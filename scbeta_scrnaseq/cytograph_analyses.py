import logging
from typing import *
import numpy as np
import scipy as sp
import pandas as pd

import loompy
import cytograph as cg
import igraph as ig

import louvain

from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import NearestNeighbors

from scbeta_scrnaseq import utils
import scbeta_scrnaseq.cytograph_inmem_utils as cgm


def highvar_pca(ds, vals=None, 
            n_highvar_pcs=50,
            n_highvar_genes=2000,
            train_cells=None,
            namespace='', seed=0):
    # nn_k=200
    # n_highvar_pcs=30
    # n_marker_pcs=30
    # seed=0
    # n_highvar_genes = 2500
    # n_markers = 25
    # classifier_prob_threshold=0.66

    np.random.seed(seed)

    # Load counts to memory
    if vals is None:
        vals = sp.sparse.csr_matrix(ds.layers[""][:, :])

    # Flag genes detected in fewer than 10 cells as invalid.
    # nnz = (vals > 0).sum(1).A.ravel()
    # ds.ra._Valid = nnz > 10

    # Normalize each cells' counts to 10000.
    # Then, mean-center and variance-stabilize each gene.
    normalizer = cgm.CountsNormalizer()
    normalizer.fit(vals, cells=train_cells)

    # Identify high-variance genes. 
    features = cgm.FeatureSelection()
    features.fit(ds, vals, cells=train_cells)
    genes = features.select_genes(n_genes=n_highvar_genes, valid_genes=ds.ra._Valid.astype(bool))
    n_cells = len(normalizer.totals)

    # Fit and transform all cells using PCA
    high_var_pca = cgm.PCAProjection(genes=genes, max_n_components=n_highvar_pcs)
    high_var_pca.fit(vals, normalizer=normalizer, cells=train_cells)
    high_var_pca_transformed = high_var_pca.transform(vals, normalizer=normalizer)

    ds.ca[f"{namespace}NormalizerTotals"] = normalizer.totals
    ds.ra[f"{namespace}NormalizerMu"] = normalizer.mu
    ds.ra[f"{namespace}NormalizerSd"] = normalizer.sd

    is_high_var_gene = np.zeros(ds.ra._Valid.shape)
    is_high_var_gene[genes] = 1
    ds.ra[f"{namespace}HighVarGene"] = is_high_var_gene

    ds.ca[f"{namespace}HighVarPCA"] = high_var_pca_transformed

def compute_tsne(ds, tsne_input, perplexity=100, namespace="", seed=0):
    np.random.seed(seed)
    ds.ca[r"{namespace}TSNE"] = cg.TSNE(perplexity=perplexity).layout(tsne_input)


def mknn_leiden_labels(input_space, valid_cells=None,
        nn_k=250, clustering_resolution=3.0, min_cluster_size=10,
        seed=0,
        ):

    mknn = cgm.get_mknn(input_space, n_neighbors=nn_k, correct_disconnected_nodes=False, seed=seed)
    leiden_labels = cgm.get_leiden(mknn, resolution_parameter=clustering_resolution,
                                      min_cluster_size=min_cluster_size,
                                      seed=seed)

    if valid_cells is not None:
        _labels = np.zeros_like(valid_cells.astype(int))
        _labels[valid_cells.astype(bool)] = leiden_labels.original
        leiden_labels = cgm.CellLabels(_labels.astype(int))

    return leiden_labels

def recluster_label(input_space, initial_labels, labels_to_refine, 
        nn_k=200, clustering_resolution=0.5, min_cluster_size = 50, seed=0):

    refining_filter = np.isin(initial_labels, labels_to_refine).astype(bool)

    refined_labels = mknn_leiden_labels(input_space[refining_filter],
                                           refining_filter, min_cluster_size=min_cluster_size,
                                           nn_k=nn_k, clustering_resolution=clustering_resolution, seed=seed)

    def merge_initial_refined_labels(initial, refined, labels_to_refine):
        merged = initial.copy()
        _refined = refined + initial.max() + 1
        for lb in labels_to_refine:
            merged[initial==lb] = _refined[initial==lb]

        label_size = np.bincount(merged)

        mapping = {v:k for k,v in enumerate([0,] + list(np.argsort(-label_size[1:])+1))}
        merged = np.array([mapping[i] for i in merged])

        return merged

    return cgm.CellLabels(merge_initial_refined_labels(initial_labels, refined_labels.original, labels_to_refine))




def low_qual_cluster_detector(ds, vals=None, input_space_attr='RawHighVarPCA',
            n_pcs=50, nn_k=200, 
            clustering_resolution = 3.0,
            n_markers = 100, seed=0,
            ):
    np.random.seed(seed)
    if vals is None:
        vals = sp.sparse.csr_matrix(ds.layers[""][:, :])

    input_space = ds.ca[input_space_attr]
    input_space = input_space[:, :min(n_pcs, input_space.shape[1])]

    cluster_labels = mknn_leiden_labels(input_space, nn_k=nn_k, clustering_resolution=clustering_resolution, seed=seed)

    label_has_markers = np.zeros(len(cluster_labels.le.classes_))
    marker_selection_rounds = 10
    marker_selection_results = []
    for __ in range(marker_selection_rounds):
        marker_selection = cgm.MarkerSelection()
        marker_selection.fit(vals[:, cluster_labels.is_labelled], cluster_labels.encoded)

        is_label_marker = marker_selection.select_markers(n_markers = n_markers,
                                                          valid_genes = ds.ra._Valid,
                                                          qval_threshold=0.001)
        label_has_markers += (is_label_marker.sum(0) > n_markers/2)
        marker_selection_results.append(is_label_marker.sum(0))

    ds.attrs['LabelsMarkerCounts'] = label_has_markers
    
    marker_selection_results = np.array(marker_selection_results)
    ds.attrs['MarkerSelectionResults'] = marker_selection_results

    cells_passing_qual_filter = np.ones_like(ds.ca._Valid)
    cells_passing_qual_filter[cluster_labels.original == 0] = 0

    log_totals = np.log10(vals.sum(0).A.reshape(-1))
    for encoded_label in np.where(label_has_markers == 0)[0]:
        orig_label = cluster_labels.le.classes_[encoded_label]
        in_cluster = cluster_labels.original == orig_label
        cells_passing_qual_filter[in_cluster] = 0

    ds.ca["QualFilterLabels"] = cluster_labels.original
    ds.ca["__Filter1__Automated_quality"] = cells_passing_qual_filter


import numpy_groupies as npg

def pseudobulk_from_label(ds, agg_labels, norm_total = 10000):
    label_marker_counts = npg.aggregate(agg_labels.encoded,
                 ds.vals[:, agg_labels.is_labelled].A,
                 func='sum', axis=1)
    label_total_counts = npg.aggregate(agg_labels.encoded,
                 ds.vals[:, agg_labels.is_labelled].sum(0).A.ravel(),
                 func='sum')

    label_norm_counts = ((label_marker_counts/label_total_counts) * norm_total).T
    label_norm_counts = pd.DataFrame(label_norm_counts, columns=ds.ra.Gene, index=agg_labels.le.classes_)
    return label_norm_counts


def expressed_fraction_from_label(ds, agg_labels, frac_of_max=0.01):

    expr_threshold = ds.vals.max(1).A * frac_of_max
    expr_threshold[expr_threshold < 1] = 0


    detected_frac = npg.aggregate(agg_labels.encoded,
                 ds.vals[:, agg_labels.is_labelled].A > expr_threshold,
                 func='mean', axis=1)

    detected_frac = pd.DataFrame(detected_frac.T, columns=ds.ra.Gene, index=agg_labels.le.classes_)
    return detected_frac


def update_labels_with_classifier(ds, labels, input_space_attr="FinalHighVarPCA",
                                  do_not_classify_filter=None, 
                                  recovery_classification_threshold=0.66, seed=0):
    classifier = RandomForestClassifier(max_depth=20, n_estimators=100, oob_score=True, random_state=seed)
    classifier.fit(ds.ca[input_space_attr][labels.is_labelled], labels.encoded)
    pred = classifier.predict_proba(ds.ca[input_space_attr])
    
    classifier_labels = np.zeros_like(labels.original)

    # Use out-of-bootstrap prediction for cells in training set. 
    classifier_labels[labels.is_labelled] = labels.le.inverse_transform(classifier.oob_decision_function_.argmax(1))

    recovered = (pred.max(1)[labels.is_unlabelled] > recovery_classification_threshold)
    if do_not_classify_filter is not None:
        recovered = recovered & (do_not_classify_filter[labels.is_unlabelled]==0)
    recovered = labels.is_unlabelled[recovered]
    
    classifier_labels[recovered] = labels.le.inverse_transform(pred.argmax(1)[recovered])
    return cgm.CellLabels(classifier_labels, null_label="")
