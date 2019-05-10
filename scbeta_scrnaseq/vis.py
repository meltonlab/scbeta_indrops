
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

import logging
# set DEBUG for everything
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger('matplotlib')
# set WARNING for Matplotlib
logger.setLevel(logging.WARNING)

import palettable
import scbeta_scrnaseq.cytograph_inmem_utils as cgm

def get_fig_and_gs(n_plots, n_cols, square_size=4):
    n_rows = int(np.ceil(n_plots/n_cols))

    fig = plt.figure(figsize=(square_size*n_cols, square_size*n_rows))
    gs = gridspec.GridSpec(n_rows, n_cols)
    return fig, gs


def scatter_cont(pos, cont, title='', ax=None, **kwargs):
    if ax is None:
        fig = plt.figure(figsize=(4,4))
        ax = fig.add_subplot(111)
    if title:
        ax.set_title(title)
    
    x = pos[:, 0]
    y = pos[:, 1]
    
    default_args = dict(edgecolor='none', s=10, alpha=0.8, cmap='RdBu_r')
    default_args.update(kwargs)
        
    ax.scatter(x, y, c=cont, 
                 **default_args)
        
    return ax

import matplotlib.patheffects as PathEffects
def label_cat_on_ax(ax, pos, cat, min_elements_for_label=0):
    for _c in sorted(set(cat))[::-1]:
        _f = (cat == _c)
        _x = pos[:, 0][_f]
        _y = pos[:, 1][_f]

        if len(_x) > min_elements_for_label:
            txt = ax.text(np.median(_x), np.median(_y), str(_c), fontweight='bold', ha='center', va='center')
            txt.set_path_effects([PathEffects.withStroke(linewidth=3, foreground='w')])

def scatter_cat(pos, cat, title='', ax=None,
                min_elements_for_label=0,
                int_categories = False, fix_max_category=0,
                color_sequence = None,
                **kwargs):

    if color_sequence is None:
        color_sequence = palettable.colorbrewer.qualitative.Set1_9.colors + \
            palettable.colorbrewer.qualitative.Pastel1_9.colors + \
            palettable.colorbrewer.qualitative.Set3_12.colors + \
            palettable.colorbrewer.qualitative.Set2_8.colors

    if ax is None:
        fig = plt.figure(figsize=(4,4))
        ax = fig.add_subplot(111)
    if title:
        ax.set_title(title)
    
    if int_categories:
        _cat_elements = range(max(fix_max_category, max(cat)))
        cat = cat.astype(int)
    else:
        _cat_elements = sorted(set(cat))
    
    default_args = dict(edgecolor='none', s=5, alpha=1.0)
    default_args.update(kwargs)
    
    for _c, color in zip(_cat_elements, color_sequence):
        _f = (cat == _c)
        _x = pos[:, 0][_f]
        _y = pos[:, 1][_f]
        
        ax.scatter(_x, _y, c=(np.array(color)/255.).reshape(1,-1), 
                  **default_args)
        
    label_cat_on_ax(ax, pos, cat, min_elements_for_label=min_elements_for_label)
    return ax




import matplotlib.gridspec as gridspec
def plot_genes(display_genes, tsne_proj, indrops_obj):        
    x = tsne_proj[0]
    y = tsne_proj[1]

    cmap = palettable.cartocolors.sequential.Sunset_7.mpl_colormap

    n_plots = len(display_genes)
    n_cols = 5
    n_rows = int(np.ceil(n_plots/n_cols))

    square_size = 2.5
    fig = plt.figure(figsize=(square_size*n_cols, square_size*n_rows))
    gs = gridspec.GridSpec(n_rows, n_cols)

    for ci, g in enumerate(display_genes):
        ax = fig.add_subplot(gs[ci], xlabel='', ylabel='',
                             xticklabels=[], yticklabels=[], title=g)
        ax.set_title(g, fontsize=15)
        if g not in indrops_obj.tpm.columns:
            continue

        c = indrops_obj.log10_tpm_g(g)[x.index]
        ax.scatter(x, y, c=c, edgecolor='none', s=5, alpha=0.8, cmap=cmap)
    return fig


def plot_auto_qual_filter_by_time_point(ds_dict, square_size = 5, show_genes = [], 
        tsne_attr = 'RawHighVarTSNE', ):

    n_cols = len(ds_dict)

    n_meta_plots = 4
    n_rows = len(show_genes) + n_meta_plots


    fig = plt.figure(figsize=(square_size*n_cols, square_size*n_rows))
    gs = gridspec.GridSpec(n_rows, n_cols)

    for ci,k in enumerate(ds_dict.keys()):
        _vds = ds_dict[k]

        tsne_pos = _vds.ca[tsne_attr]

        # Plot the subclusters
        ax = fig.add_subplot(gs[0, ci], xticks=[], yticks=[], title=f'timepoint {k}', ylabel='Quality Filter Cluster')
        scatter_cat(tsne_pos,  _vds.ca.QualFilterLabels, ax=ax, s=5)
        
        # Plot the automated quality filter results
        ax = fig.add_subplot(gs[1, ci], xticks=[], yticks=[], title=f'timepoint {k}', ylabel='Autofilter results')
        scatter_cont(tsne_pos, _vds.ca.__Filter1__Automated_quality, ax=ax, s=5)
        
        # Plot the depth
        ax = fig.add_subplot(gs[2, ci], xticks=[], yticks=[], title=f'timepoint {k}', ylabel='log10(counts) per cell')
        scatter_cont(_vds.ca.RawHighVarTSNE, np.log10(_vds.ca.RawNormalizerTotals),
                        ax=ax, cmap='RdBu', s=5)
        
        # Plot the batch
        ax = fig.add_subplot(gs[3, ci], xticks=[], yticks=[], title=f'timepoint {k}', ylabel='batch')
        scatter_cat(_vds.ca.RawHighVarTSNE, _vds.ca.CellBatch,
                        ax=ax, color_sequence=palettable.colorbrewer.qualitative.Paired_10.colors, s=5)
        

        # Plot the gene expression
        _v_norm = cgm.CountsNormalizer(totals=_vds.ca.RawNormalizerTotals, mean_center=False, normalize_variance=False)
        for gi, g in enumerate(show_genes):
            _gi = np.where(_vds.ra.Gene==g)[0]
            if len(_gi)==0:
                continue
            g_expr = np.log10(_v_norm.transform(_vds.vals, genes=_gi) + 0.1).ravel()
            vmax = max(max(g_expr), 2)
            _co = np.argsort(g_expr)
            
            ax = fig.add_subplot(gs[gi+n_meta_plots, ci], xticks=[], yticks=[], title=f'timepoint {k}', ylabel=f'{g} expr')
                
            scatter_cont(_vds.ca.RawHighVarTSNE[_co], g_expr[_co], s=5,
                     cmap=palettable.cartocolors.sequential.Sunset_7.mpl_colormap, ax=ax,
                             vmin = -1, vmax=vmax)
            label_cat_on_ax(ax, tsne_pos, _vds.ca.QualFilterLabels)


def plot_gene_by_time_point(ds_dict,  show_genes = [], square_size = 5,
        tsne_attr = 'RawHighVarTSNE', ):

    n_cols = len(ds_dict)

    n_rows = len(show_genes)


    fig = plt.figure(figsize=(square_size*n_cols, square_size*n_rows))
    gs = gridspec.GridSpec(n_rows, n_cols)

    for ci,k in enumerate(ds_dict.keys()):
        _vds = ds_dict[k]

        tsne_pos = _vds.ca[tsne_attr]

        # Plot the gene expression
        _v_norm = cgm.CountsNormalizer(totals=_vds.ca.RawNormalizerTotals, mean_center=False, normalize_variance=False)
        for gi, g in enumerate(show_genes):
            _gi = np.where(_vds.ra.Gene==g)[0]
            if len(_gi)==0:
                continue
            g_expr = np.log10(_v_norm.transform(_vds.vals, genes=_gi) + 0.1).ravel()
            vmax = max(max(g_expr), 2)
            _co = np.argsort(g_expr)
            
            ax = fig.add_subplot(gs[gi, ci], xticks=[], yticks=[], title=f'timepoint {k}', ylabel=f'{g} expr')
                
            scatter_cont(tsne_pos, g_expr, s=5,
                     cmap=palettable.cartocolors.sequential.Sunset_7.mpl_colormap, ax=ax,
                             vmin = -1, vmax=vmax)

def plot_cat_by_time_point(ds_dict, cat, square_size = 5,
        tsne_attr = 'RawHighVarTSNE', min_elements_for_label=50, 
        int_categories = True):

    n_cols = len(ds_dict)
    n_rows = 1

    fig = plt.figure(figsize=(square_size*n_cols, square_size*n_rows))
    gs = gridspec.GridSpec(n_rows, n_cols)

    for ci,k in enumerate(ds_dict.keys()):
        _vds = ds_dict[k]

        tsne_pos = _vds.ca[tsne_attr]

        v_f = np.isin(_vds.ca.CellID, list(cat.index))
        v_cells = _vds.ca.CellID[v_f]
        v_pos = _vds.ca.RawHighVarTSNE[v_f]

        fix_max_category = max(cat)+1 if int_categories else None

        ax = fig.add_subplot(gs[ci], xticks=[], yticks=[], title=f'timepoint {k}')
        scatter_cat(v_pos, cat[v_cells], min_elements_for_label=min_elements_for_label, 
                int_categories = int_categories, fix_max_category=fix_max_category, ax=ax)


def plot_cont_by_time_point(ds_dict, cont, square_size = 5,
        tsne_attr = 'RawHighVarTSNE', **attrs):

    n_cols = len(ds_dict)
    n_rows = 1

    fig = plt.figure(figsize=(square_size*n_cols, square_size*n_rows))
    gs = gridspec.GridSpec(n_rows, n_cols)

    for ci,k in enumerate(ds_dict.keys()):
        _vds = ds_dict[k]

        tsne_pos = _vds.ca[tsne_attr]

        v_f = np.isin(_vds.ca.CellID, list(cont.index))
        v_cells = _vds.ca.CellID[v_f]
        v_pos = _vds.ca.RawHighVarTSNE[v_f]

        ax = fig.add_subplot(gs[ci], xticks=[], yticks=[], title=f'timepoint {k}')
        scatter_cont(v_pos, cont[v_cells], ax=ax, **attrs)


def plot_attr_by_time_point(ds_dict, attr, attr_type='cat', square_size = 5,
        tsne_attr = 'RawHighVarTSNE', min_elements_for_label=50, 
        int_categories = True, dpi=100, **attrs):

    n_cols = len(ds_dict)
    n_rows = 1

    fig = plt.figure(figsize=(square_size*n_cols, square_size*n_rows), dpi=dpi)
    gs = gridspec.GridSpec(n_rows, n_cols)

    for ci,k in enumerate(ds_dict.keys()):
        _vds = ds_dict[k]

        tsne_pos = _vds.ca[tsne_attr]

        plot_value = _vds.ca[attr]

        ax = fig.add_subplot(gs[ci], xticks=[], yticks=[])
        if attr_type == 'cat':
            fix_max_category = max(plot_value)+1 if int_categories else None
            scatter_cat(tsne_pos, plot_value, min_elements_for_label=min_elements_for_label, 
                    int_categories = int_categories, fix_max_category=fix_max_category, ax=ax)
        elif attr_type == 'cont':
            scatter_cont(tsne_pos, plot_value, ax=ax, **attrs)

