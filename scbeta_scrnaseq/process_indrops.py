import subprocess
import os, sys

import numpy as np
import scipy as sp
import pandas as pd

from collections import OrderedDict, defaultdict

from scbeta_scrnaseq import utils

def fast_counts_parse(filename, cache_suffix='.parsed.npz', verbose=False, save_cached_file=True):
    preparsed_filename = filename + cache_suffix
    if os.path.isfile(preparsed_filename):
        # Read data from cache
        with np.load(preparsed_filename) as npz:
            return {
                'sparse_counts': sp.sparse.coo_matrix((npz['c_'], (npz['i_'], npz['j_'])),
                    shape=(len(npz['barcodes']), len(npz['genes']))).tocsr(),
                'genes' : npz['genes'].copy(),
                'barcodes' : npz['barcodes'].copy()
            }

    else:
        # Extract from the counts file
        genes = None
        barcodes = []
        counts = []
        i_ = []
        j_ = []
        c_ = []
        for line_i,line in enumerate(subprocess.Popen(["gzip", "--stdout", "-d", filename], stdout=subprocess.PIPE).stdout):
            line = line.decode()
            line = line[:-1].split('\t')
            if not genes:
                assert line[0] == 'barcode'
                genes = line[1:]
                i=0
                continue
            i = line_i-1
            barcodes.append(line[0])
            for j, c in enumerate(line[1:]):
                if c != '0':
                    i_.append(i)
                    j_.append(j)
                    c_.append(float(c))
            if i%100==0 and verbose:
                print(i, 'rows parsed.')

        i_ = np.array(i_)
        j_ = np.array(j_)
        c_ = np.array(c_)
        genes = np.array(genes)
        barcodes = np.array(barcodes)
        if save_cached_file:
            np.savez(filename+'.parsed.npz',
                i_=i_, j_=j_, c_=c_,
                genes=genes, barcodes=barcodes)
        return {
                'sparse_counts':  sp.sparse.coo_matrix((c_, (i_, j_)),
                    shape=(len(barcodes), len(genes))).tocsr(),
                'genes' : genes,
                'barcodes' : barcodes,
            }

def fast_counts_parse_to_df(filename, verbose=False):
    data = fast_counts_parse(filename, verbose=verbose)
    return pd.DataFrame(data['sparse_counts'].toarray(), index=data['barcodes'], columns=data['genes'])

class SparsifiedIndropsAnalysis(utils.CachedAnalysis):
    cached_attributes = OrderedDict([
        ('counts', 'csr.df.npz'),
        ('raw_counts_summary',  'df.npz'),
        ])

    other_samples_for_gene_exclusion = ()
    maximal_mito_reads_fraction = 0.33
    log10_local_lib_size_threshold = 3.40
    expected_fano_threshold = False

    cell_quality_thr = 1000
    gene_quality_thr = 10

    counts_input_files = []

    def process(self):
        np.random.seed(0)
        (counts, gene_names, cell_names), mito_counts, raw_counts_summary = self.get_final_counts()
        counts = utils.SparseCSRDataFrame(counts, columns=gene_names, index=cell_names)
        
        return locals()

    def get_raw_counts(self):
        genes = None
        cells = []
        sparse_counts = []

        for input_prefix, cell_prefix in self.counts_input_files:

            c = fast_counts_parse(os.path.join(self.root_path, '..', 'indrops_raw', '%s.counts.tsv.gz' % input_prefix))
            
            cells += ['%s.%s'%(cell_prefix, bc) for bc in c['barcodes']]
            
            if genes is not None:
                assert (c['genes'] == genes).all()
            else:
                genes = c['genes']
            
            sparse_counts.append(c['sparse_counts'])

        counts = sp.sparse.vstack(sparse_counts)
        return counts, np.array(cells), genes
    
    def get_final_counts(self):
        counts, cell_names, gene_names = self.get_raw_counts()

        not_junk_genes = [(not gene.startswith('MT-')) &
                 ('MTRNR' not in gene) &
                  ('.' not in gene) & 
                  (not gene.startswith('LINC')) &
                   ('-AS' not in gene) &
                    (not (gene.startswith('RP') &
                     ('-' in gene) )) for gene in gene_names]
        not_junk_genes = np.array(not_junk_genes)

        final_counts = counts[:, not_junk_genes]
        
        library_sizes = final_counts.sum(axis=1)

        mito_genes = [(gene.startswith('MT-') or ('MTRNR' in gene)) for gene in gene_names]
        mito_genes = np.array(mito_genes)

        mito_counts = counts[:, mito_genes].sum(axis=1)
        mito_ratio = mito_counts/library_sizes

        gene_detection_count = (final_counts>0).sum(0)
        cell_gene_count = (final_counts>0).sum(1)

        cell_quality_thr = self.cell_quality_thr
        gene_quality_thr = self.gene_quality_thr

        cell_filter = (library_sizes > cell_quality_thr) & (mito_ratio < self.maximal_mito_reads_fraction)
        gene_filter = gene_detection_count > gene_quality_thr

        cell_filter = np.array(cell_filter).ravel()
        gene_filter = np.array(gene_filter).ravel()

        # Compute genes to exclude
        max_frac = (sp.sparse.diags(np.array(1./library_sizes[cell_filter]).ravel()) * final_counts[cell_filter, :]).max(0)
        exclude_if_above_fr = 0.02

        genes_to_exclude_from_lib_size = pd.Series(max_frac.toarray().ravel() > exclude_if_above_fr, index=gene_names[not_junk_genes])
        genes_to_exclude_from_lib_size.to_csv(os.path.join(self.root_path, '%s.exclude_from_lib_size.csv'%self.__class__.__name__))


        _counts_tuple = (final_counts[cell_filter, :][:, gene_filter],
                            gene_names[not_junk_genes][gene_filter],
                            cell_names[cell_filter])
        _mito_counts = pd.Series(mito_counts[cell_filter].A.ravel(), cell_names[cell_filter])

        raw_counts_summary = pd.DataFrame([library_sizes.A.ravel(), mito_counts.A.ravel(), cell_gene_count.A.ravel()],
            index=['all_counts', 'mito_counts', 'n_genes_detected'],
            columns = cell_names).T

        return _counts_tuple, _mito_counts, raw_counts_summary



