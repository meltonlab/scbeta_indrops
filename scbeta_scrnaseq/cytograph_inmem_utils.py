import logging
import numpy as np
import scipy as sp
from typing import *
import loompy
import cytograph as cg
import pandas as pd

import igraph as ig
from sklearn import preprocessing

class CountsNormalizer:
    """
    Normalize and optionally standardize a dataset, dealing properly with edge cases such as division by zero.
    """
    def __init__(self,
            mu: np.ndarray = None,
            sd: np.ndarray = None,
            totals = None,
            mean_center: bool = True,
            normalize_variance: bool = True,
            level: int = 10000) -> None:

        self.sd = sd  # type: np.ndarray
        self.mu = mu  # type: np.ndarray
        self.totals = totals  # type: np.ndarray
        self.level = level
        self.mean_center = mean_center
        self.normalize_variance = normalize_variance

    def norm_factor_from_totals(self, totals):
        norm_factor = (self.level/totals)
        norm_factor[~np.isfinite(norm_factor)] = 0
        norm_factor = sp.sparse.diags(norm_factor).tocsr()
        return norm_factor

    def fit(self, vals: sp.sparse.csr_matrix, cells: np.ndarray = None) -> None:
        totals = vals.sum(0).A.reshape(-1)
        norm_factor = self.norm_factor_from_totals(totals)

        vals = vals.dot(norm_factor).T
        if cells is not None:
            vals = vals[cells, :]

        #Standard scaler is particularly fast and works with sparse arrays (it computes the mean even if it is not used in rescaling).
        scaler = preprocessing.StandardScaler(with_mean=False)
        scaler.fit(vals)

        self.mu = scaler.mean_
        self.sd = np.sqrt(scaler.var_)
        self.totals = totals

    def transform(self, vals: sp.sparse.csr_matrix, cells: np.ndarray = None, genes: np.ndarray = None) -> np.ndarray:
        """
        Normalize a matrix using the previously calculated aggregate statistics

        Args:
            vals (sp.sparse.csr_matrix):     Matrix of shape (n_genes, n_cells)
            cells (ndarray):    Optional indices of the cells that are represented in vals
            cells (ndarray):    Optional indices of the genes that are represented in vals

        Returns:
            vals_adjusted (ndarray):    The normalized values
        """

        if genes is not None:
            vals = vals[genes, :]
        else:
            genes = np.arange(vals.shape[0])

        norm_factor = norm_factor = self.norm_factor_from_totals(self.totals)
        if cells is not None:
            vals = vals[:, cells]
            norm_factor = norm_factor[cells, :][:, cells]

        vals = vals.dot(norm_factor).A

        if self.mean_center:
            vals = vals - self.mu[genes][:, None]
        if self.normalize_variance:
            vals = div0(vals.T, self.sd[genes]).T
            
        return vals
    
def div0(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """ ignore / 0, div0( [-1, 0, 1], 0 ) -> [0, 0, 0] """
    with np.errstate(divide='ignore', invalid='ignore'):
        c = np.true_divide(a, b)
        c[~np.isfinite(c)] = 0  # -inf inf NaN
    return c


from sklearn.svm import SVR

class FeatureSelection:
    def __init__(self) -> None:
        self.genes = None  # type: np.ndarray
        self.mu = None  # type: np.ndarray
        self.sd = None  # type: np.ndarray
        self.totals = None  # type: np.ndarray

    def fit(self, ds: loompy.LoomConnection, vals: sp.sparse.csr_matrix,
            cells: np.ndarray = None,
            mu: np.ndarray = None,
            sd: np.ndarray = None) -> np.ndarray:
        """
        Fits a noise model (CV vs mean)

        Args:
            ds (LoomConnection):    Dataset
            n_genes (int):  number of genes to include
            cells (ndarray): cells to include when computing mean and CV (or None)
            mu, std:        Precomputed mean and standard deviations (optional)

        Returns:
            ndarray of selected genes (list of ints)
        """
        if cells is not None:
            vals = vals[:, cells]
        self.fit_cells = cells

        if mu is None or sd is None:
            scaler = preprocessing.StandardScaler(with_mean=False)
            scaler.fit(vals.T)
            mu = scaler.mean_
            sd = np.sqrt(scaler.var_)

        # if "_Valid" in ds.ra:
        #     valid = ds.ra._Valid == 1
        # else:
        #     valid = np.ones(ds.shape[0], dtype='bool')
        # valid = ((valid) & (mu > self.min_expr)).astype(int)
        ok = (mu > 0) & (sd > 0)

        self.mu = mu
        self.sd = sd
        self.ok = ok
        # self.valid = valid

        cv = sd[ok] / mu[ok]
        log2_m = np.log2(mu[ok])
        log2_cv = np.log2(cv)

        svr_gamma = 1000. / len(mu[ok])
        clf = SVR(gamma=svr_gamma)
        clf.fit(log2_m[:, np.newaxis], log2_cv)
        fitted_fun = clf.predict
        # Score is the relative position with respect of the fitted curve
        fitted_val = fitted_fun(log2_m[:, np.newaxis])
        score = log2_cv - fitted_val
        # score = score * valid[ok]

        self.fitted_val = fitted_val
        self.score = score
        self.log2_cv = log2_cv
        self.log2_m = log2_m

    def select_genes(self, n_genes=1000, min_expr=0.001, valid_genes=None):
        _valid = (self.mu > min_expr)
        if valid_genes is not None:
            _valid = _valid & valid_genes
        _valid_score = self.score * _valid[self.ok].astype(float)

        picked_genes = np.where(self.ok)[0][np.argsort(_valid_score)][-n_genes: ]
        return picked_genes


from sklearn.decomposition import PCA

class PCAProjection:
    """
    Project a dataset into a reduced feature space using PCA. The projection can be fit
    to one dataset then used to project another. To work properly, both datasets must be normalized in the same 
    way prior to projection.
    """
    def __init__(self, genes: np.ndarray, max_n_components: int = 50) -> None:
        """
        Args:
            genes:              The genes to use for the projection
            max_n_components:   The maximum number of projected components
            nng                 Non-neuronal genes, to be zeroed in neurons (where TaxonomyRank1 == "Neurons")
        """
        self.genes = genes
        self.n_components = max_n_components

        self.cells = None  # type: np.ndarray
        self.pca = None  # type: IncrementalPCA
        self.sigs = None  # type: np.ndarray
        # self.scan_batch_size = scan_batch_size

    def fit(self, vals: sp.sparse.csr_matrix, normalizer: cg.Normalizer, cells: np.ndarray = None) -> None:
        n_cells = vals.shape[1] if cells is None else cells.shape[0]
        n_genes = self.genes.shape[0]

        self.pca = PCA(n_components=self.n_components)
        norm_vals = normalizer.transform(vals, genes=self.genes, cells=cells)
        self.pca.fit(norm_vals.T)

    def transform(self, vals: sp.sparse.csr_matrix, normalizer: cg.Normalizer, cells: np.ndarray = None) -> np.ndarray:
        
        n_cells = vals.shape[1] if cells is None else cells.shape[0]
        n_genes = self.genes.shape[0]

        norm_vals = normalizer.transform(vals, genes=self.genes, cells=cells)
        transformed = self.pca.transform(norm_vals.T)
        return transformed

    def fit_transform(self, vals: sp.sparse.csr_matrix, normalizer: cg.Normalizer, cells: np.ndarray = None) -> np.ndarray:
        self.fit(vals, normalizer, cells)
        return self.transform(vals, normalizer, cells)


import logging
import numpy as np
from typing import *
from sklearn.svm import SVR
from statsmodels.sandbox.stats.multicomp import multipletests
import loompy
import cytograph as cg


class MarkerSelection:
    def __init__(self, n_markers: int = 10) -> None:
        self.alpha = 0.1

    def fit(self, vals: sp.sparse.csr_matrix, labels: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Finds n_markers genes per cluster using enrichment score

        Args:
            vals:    Dataset, as a sparse array loaded in memory

        Returns:
            ndarray of selected genes (list of ints)
            ndarray of enrichment scores
            ndarray of FDR-corrected P values (i.e. q values)
        """
        # Get the observed enrichment statistics

        enrichment = self._fit(vals, labels)
        
        rand_labels = np.random.permutation(labels)
        null_enrichment = self._fit(vals, rand_labels)
        qvals = np.zeros_like(enrichment)
    
        for ix in range(enrichment.shape[1]):
            null_values = null_enrichment[:, ix]
            null_values.sort()
            values = enrichment[:, ix]
            pvals = 1 - np.searchsorted(null_values, values) / values.shape[0]
            (_, q, _, _) = multipletests(pvals, self.alpha, method="fdr_bh")
            qvals[:, ix] = q

        self.enrichment = enrichment
        self.qvals = qvals

    def _fit(self, vals: sp.sparse.csr_matrix, labels: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """
        Finds n_markers genes per cluster using enrichment score

        Args:
            vals:    Dataset, as a sparse array loaded in memory

        Returns:
            ndarray of selected genes (list of ints)
            ndarray of enrichment scores
        """
        n_labels = max(labels) + 1
        n_cells = vals.shape[1]
        labels_mat = sp.sparse.coo_matrix((np.ones_like(labels), (np.arange(n_cells), labels)))

        # Number of cells per cluster
        sizes = np.bincount(labels, minlength=n_labels)

        means = (vals.dot(labels_mat).A)/sizes
        nnz = (vals > 0).dot(labels_mat).A

        means_overall  = vals.mean(1).A.ravel()
        nnz_overall = (vals > 0).sum(1).A.ravel()

        f_nnz = nnz / sizes
        f_nnz_overall = nnz_overall / n_cells

        means_other = ((means_overall * n_cells)[None].T - (means * sizes)) / (n_cells - sizes)
        f_nnz_other = ((f_nnz_overall * n_cells)[None].T - (f_nnz * sizes)) / (n_cells - sizes)

        enrichment = (f_nnz + 0.1) / (f_nnz_other + 0.1) * (means + 0.01) / (means_other + 0.01)
        return enrichment

    def select_markers(self, n_markers=50, qval_threshold=0.001, valid_genes = None):
        is_marker = (pd.DataFrame(self.enrichment).rank(ascending=False).values <= n_markers) & (self.qvals < qval_threshold)
        if valid_genes is not None:
            is_marker = is_marker & (np.tile(valid_genes, (self.enrichment.shape[1], 1)).T)
        return is_marker



def clean_labels(labels, min_cluster_size=10, zero_is_null=None):
    bigs = np.where(np.bincount(labels) >= min_cluster_size)[0]
    mapping = {k: v for v, k in enumerate(bigs)}
    labels = np.array([mapping[x]+1 if x in bigs else 0 for x in labels])
    return labels


import numpy_groupies as npg
from sklearn import preprocessing
class CellLabels():
    def __init__(self, labels, null_label=0):
        
        self.original = labels
        self.null_label = null_label
        self.is_labelled = np.where(self.original != null_label)[0]
        self.is_unlabelled = np.where(self.original == null_label)[0]
        
        self.le = preprocessing.LabelEncoder()
        self.encoded = self.le.fit_transform(self.original[self.is_labelled])


    def pseudobulk_counts(self, vals, level=10000, gene_names=None):

        label_marker_counts = npg.aggregate(self.encoded,
                    vals[:, :].A,
                    func='sum', axis=1)
        label_total_counts = npg.aggregate(self.encoded,
                     vals.sum(0).A.ravel(),
                     func='sum')
        label_norm_counts = ((label_marker_counts/label_total_counts) * level).T

        original_grp_names = self.le.inverse_transform(np.arange(label_norm_counts.shape[0]))

        label_norm_counts = pd.DataFrame(label_norm_counts,
            index = original_grp_names, 
            columns = gene_names,
            )


        return label_norm_counts



from sklearn.neighbors import NearestNeighbors        
def get_mknn(X, n_neighbors=100, correct_disconnected_nodes=False, seed=0):
    np.random.seed(seed)
    nn = NearestNeighbors(n_neighbors=n_neighbors, algorithm="ball_tree", n_jobs=4)
    nn.fit(X)
    knn = nn.kneighbors_graph(mode='connectivity')
    knn = knn.tocoo()
    mknn = knn.minimum(knn.transpose())

    if correct_disconnected_nodes:
        is_isolated_in_mknn = (mknn.sum(1).A.ravel() == 0).astype(float)
        nearest_single_neighbor = nn.kneighbors_graph(mode='connectivity', n_neighbors=1)
        # Keep only edges from those isolated notes
        isolation_correction = sp.sparse.diags(is_isolated_in_mknn).dot(nearest_single_neighbor)
        mknn = mknn + isolation_correction

    return mknn.tocoo()

try:
    import igraph as ig
    import louvain
except:
    pass
def get_louvain(mknn, min_cluster_size=10, resolution_parameter=1.0, seed=0):
    g = ig.Graph(n=mknn.shape[0], edges=list(zip(mknn.row, mknn.col)), directed=False)

    # Louvain clustering over the mKNN graph
    louvain.set_rng_seed(seed)
    part = louvain.find_partition(g,
                louvain.RBConfigurationVertexPartition,
                resolution_parameter=resolution_parameter)
    
    return CellLabels(clean_labels(part.membership, min_cluster_size=min_cluster_size)) 

try:
    import leidenalg
except:
    pass
def get_leiden(mknn, min_cluster_size=10, resolution_parameter=1.0, seed=0, n_iterations=5):

    g = ig.Graph(n=mknn.shape[0], edges=list(zip(mknn.row, mknn.col)), directed=False)

    part = leidenalg.find_partition(g, leidenalg.RBConfigurationVertexPartition,
                seed = seed, n_iterations = n_iterations,
                resolution_parameter=resolution_parameter,
                )

    return CellLabels(clean_labels(part.membership, min_cluster_size=min_cluster_size))   


def update_cluster_based_filter(initial_filter, cluster_labels,
                                force_false=[], force_true = []):
    new_filter = initial_filter.copy().astype(int)
    for lb in force_false:
        new_filter[cluster_labels==lb] = 0
    for lb in force_true:
        new_filter[cluster_labels==lb] = 1
    return new_filter

def merge_split_dataset_filter(ds_dict, merged_ds, filter_attr):

    passing_cell_ids = []
    for di, _ds in ds_dict.items():
        passing_cell_ids += list(_ds.ca.CellID[np.where(_ds.ca[filter_attr] > 0)[0]])

    merged_filter = np.isin(merged_ds.ca.CellID, passing_cell_ids)
    merged_ds.ca[filter_attr] = merged_filter
