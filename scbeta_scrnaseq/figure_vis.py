
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
from collections import defaultdict

def to_array_col(color):
    return np.array(color)/255.


tsne_legend_ms = 15
heatmap_legend_ms = 15


class CoreColors():

    def __init__(self):

        cnames = ['red', 'blue', 'green', 'purple', 'orange', 'yellow', 'brown', 'pink', 'grey']


        for cname,c in zip(cnames, palettable.colorbrewer.qualitative.Set1_9.colors):
            setattr(self, cname, np.array(c)/255.)
            
        for cname, c in zip(['blue', 'green', 'red', 'orange', 'purple'],
                            palettable.colorbrewer.qualitative.Paired_10.colors[::2]):
            setattr(self, 'pale_'+cname, np.array(c)/255.)
            
        self.pale_brown = to_array_col([193, 113, 68])
        self.pale_brown = to_array_col([193, 128, 108])
        self.teal = to_array_col(palettable.colorbrewer.qualitative.Set2_3.colors[0])
        self.brown_red = to_array_col(palettable.colorbrewer.qualitative.Dark2_3.colors[1])

        self.light_grey = to_array_col(palettable.colorbrewer.qualitative.Set2_8.colors[-1])

        self.dark_green = to_array_col([52, 117, 50])
        self.dark_yellow = to_array_col(palettable.colorbrewer.qualitative.Dark2_6.colors[-1])
        self.halfdark_green = (self.green + self.dark_green)/2

        self.dark_grey = np.array([0.5, 0.5, 0.5])
core_colors = CoreColors()

thin_legend_params = {'axes.linewidth': 0.35,
                     'xtick.major.width': 0.35,
                     'ytick.major.width': 0.35}


def setup_matplotlib_params():
    from matplotlib import rcParams
    rcParams['font.family'] = 'sans-serif'
    # rcParams['font.sans-serif'] = ['Myriad Pro']
    rcParams['font.sans-serif'] = ['Helvetica Neue']
    rcParams['axes.titlesize'] = 5
    rcParams['axes.labelsize'] = 5
    rcParams['xtick.labelsize'] = 5
    rcParams['ytick.labelsize'] = 5
    rcParams['mathtext.default'] = 'regular'
    rcParams['axes.linewidth'] = 0.5
    rcParams['xtick.minor.width'] = 0.5
    rcParams['xtick.major.width'] = 0.5
    rcParams['xtick.major.size'] = 2.0
    rcParams['xtick.minor.size'] = 1.0
    rcParams['xtick.major.pad'] = 1.0
    rcParams['ytick.minor.width'] = 0.5
    rcParams['ytick.major.width'] = 0.5
    rcParams['ytick.major.size'] = 2.0
    rcParams['ytick.minor.size'] = 1.0
    rcParams['ytick.major.pad'] = 1.0

    rcParams['pdf.fonttype'] = 42

def labels_to_color(labels, labels_data):
    labels_colors = np.ones((len(labels), 3)) * 0.5
    for cl in sorted(set(labels)):
        cl_filter = labels==cl
        cl_color = labels_data[cl]['color']
        labels_colors[cl_filter, :] = cl_color
    return labels_colors

def prepare_for_scatter(proj, labels, labels_data):
    
    defined_labels = list(labels_data.keys())

    cell_label_undefined = np.zeros(labels.shape).astype(bool)

    labels_colors = np.ones((len(labels), 3)) * 0.5
    for cl in sorted(set(labels)):
        cl_filter = labels==cl
        cl_color = labels_data[cl]['color']
        labels_colors[cl_filter, :] = cl_color
        cell_label_undefined[cl_filter] = (cl not in defined_labels)
    
    from sklearn.neighbors import NearestNeighbors
    nbrs = NearestNeighbors(n_neighbors=10, algorithm='ball_tree').fit(proj)
    nbrs_labels = labels[nbrs.kneighbors(proj, return_distance=False)]
    sort_val = (nbrs_labels.T == nbrs_labels[:, 0].T).mean(0)

    sort_val[cell_label_undefined] = sort_val[cell_label_undefined] - 1
    _o = np.argsort(sort_val)

    return proj[_o].copy(), labels_colors[_o].copy()

def dataset_label_params(dataset):

    common_label_data = {
        'scbeta': dict(color=core_colors.purple,
                     short_label=r'SC-beta',
                     long_label=r'SC-beta'),
        'ph': dict(color=core_colors.red,
                     short_label=r'SC-alpha',
                     long_label=r'SC-alpha'),
        'ec': dict(color=core_colors.blue,
                     short_label=r'SC-EC',
                     long_label=r'SC-EC'),
        'exo': dict(color=core_colors.green,
                     short_label=r'Non-endo.',
                     long_label=r'Non-endocrine'),
        'pdx1': dict(color=core_colors.pale_green,
                     short_label=r'PDX1⁺ prog.',
                     long_label=r'PDX1⁺ progenitor'),
        'nkx61': dict(color=core_colors.teal,
                     short_label=r'NKX6.1⁺ prog.',
                     long_label=r'NKX6.1⁺ progenitor'),
        'foxj1': dict(color=core_colors.pale_blue,
                     short_label=r'CHGA⁺/FOXJ1⁺',
                     long_label=r'CHGA⁺/FOXJ1⁺'),
        'sst_hhex': dict(color=core_colors.pale_red,
                     short_label=r'SST⁺/HHEX⁺',
                     long_label=r'SST⁺/HHEX⁺'),
        'fev_high_isl_neg': dict(color=core_colors.pale_brown,
                     short_label=r'FEV$^{high}$/ISL$^{low}$',
                     long_label=r'FEV$^{high}$/ISL$^{low}$'),
        '': dict(color=core_colors.grey),
    }

    dataset_specific_data = defaultdict(dict, {
        'x1_S3c': {
            'repl': dict(color=core_colors.dark_grey,
                     short_label=r'Repl. PDX1⁺',
                     long_label=r'Repl. PDX1⁺'),
        },
        'x1_S4c': {
            'repl': dict(color=core_colors.dark_grey,
                     short_label=r'Repl. NKX6.1⁺',
                     long_label=r'Repl. NKX6.1⁺'),
            'neurog3': dict(color=core_colors.orange,
                     short_label=r'NEUROG3⁺ prog.',
                     long_label=r'NEUROG3⁺ prog.'),
        },
        'x1_S5c': {
            'repl': dict(color=core_colors.dark_grey,
                     short_label=r'Repl. n.-e.',
                     long_label=r'Repl. n.-e.'),
        },
        'x1_S6c': {
            'repl': dict(color=core_colors.dark_grey,
                     short_label=r'Repl. n.-e.',
                     long_label=r'Repl. n.-e.'),
        },
        'x2_S3c': {
            'repl': dict(color=core_colors.dark_grey,
                     short_label=r'Repl. PDX1⁺',
                     long_label=r'Repl. PDX1⁺'),
        },
        'x2_S4c': {
            'pre_ph': dict(color=core_colors.red,
                     short_label=r'Repl. NKX6.1⁺',
                     long_label=r'Repl. NKX6.1⁺'),
            'repl': dict(color=core_colors.grey,
                     short_label=r'Repl. NKX6.1⁺',
                     long_label=r'Repl. NKX6.1⁺'),
            'neurog3': dict(color=core_colors.orange,
                     short_label=r'NEUROG3⁺ prog.',
                     long_label=r'NEUROG3⁺ prog.'),
        },
        'x2_S5c': {
            'repl': dict(color=core_colors.dark_grey,
                     short_label=r'Repl. n.-e.',
                     long_label=r'Repl. n.-e.'),
        },
        'x2_S6c': {
            'repl': dict(color=core_colors.dark_grey,
                     short_label=r'Repl. n.-e.',
                     long_label=r'Repl. n.-e.'),
        },
        'stage6': {
            'sst': common_label_data['sst_hhex'],
            'neurog3': dict(color=core_colors.brown,
                     short_label=r'NEUROG3⁺ prog.',
                     long_label=r'NEUROG3⁺ prog.'),
            'repl': dict(color=core_colors.dark_grey,
                     short_label=r'Repl.',
                     long_label=r'Repl.'),
            'other': dict(color=core_colors.light_grey,
                     short_label=r'Other',
                     long_label=r'Other'),
            'other__phox2a': dict(color=core_colors.light_grey,
                     short_label=r'PHOX2A⁺',
                     long_label=r'PHOX2A⁺'),
            'other__onecut3': dict(color=core_colors.light_grey,
                     short_label=r'ONECUT3⁺',
                     long_label=r'ONECUT3⁺'),
            'other__gap43': dict(color=core_colors.light_grey,
                     short_label=r'GAP43⁺',
                     long_label=r'GAP43⁺'),
        },
        'stage5': {
            'prog_nkx61': common_label_data['nkx61'],
            'prog_sox2': common_label_data['exo'],
            'fev_high_isl_low' : common_label_data['fev_high_isl_neg'],
            'neurog3_late': dict(color=core_colors.brown,
                     short_label=r'Late endo. ind.',
                     long_label=r'Late endo. ind.',),
            'neurog3_mid': dict(color=core_colors.orange,
                     short_label=r'Mid endo. ind.',
                     long_label=r'Mid endo. ind.',),
            'neurog3_early': dict(color=core_colors.pale_orange,
                     short_label=r'Early endo. ind.',
                     long_label=r'Early endo. ind.',),

            'phox2a': dict(color=core_colors.light_grey,
                     short_label=r'PHOX2A⁺',
                     long_label=r'PHOX2A⁺'),

            'exo__lgals3': dict(color=core_colors.pale_brown,
                     short_label=r'Non-endo CDX2⁺',
                     long_label=r'Non-endo CDX2⁺'),

            # 'other': dict(color=core_colors.pale_brown,
            #          short_label=r'Other',
            #          long_label=r'Other'),
        },
        'human_islets': {
            'beta': dict(color=core_colors.purple,
                         short_label=r'Islet beta',
                         long_label=r'Beta'),
            'alpha': dict(color=core_colors.red,
                         short_label=r'Islet alpha',
                         long_label=r'Alpha'),
            'gamma': dict(color=core_colors.blue,
                         short_label=r'Islet gamma',
                         long_label=r'Gamma'),
            'delta': dict(color=core_colors.green,
                         short_label=r'Islet delta',
                         long_label=r'Delta'),
            'acinar': dict(color=core_colors.pale_green,
                         short_label=r'Pancreatic acinar',
                         long_label=r'Pancreatic acinar'),
            'ductal': dict(color=core_colors.teal,
                         short_label=r'Pancreatic ductal',
                         long_label=r'Pancreatic ductal'),
        }
    })

    output_ = {**common_label_data, **dataset_specific_data[dataset]}
    for cl in output_.keys():
        output_[cl]['color_vec'] = output_[cl]['color'].reshape(1, -1)

    return output_

from sklearn.neighbors import RadiusNeighborsRegressor

def proj_regressor(proj, vals, **init_args):
    default_args = dict(
        radius = 1.0)
    default_args.update(init_args)
    
    neigh = RadiusNeighborsRegressor(**default_args)
    neigh.fit(proj, vals)
    return neigh

def scatter_as_path(X, buf=0.5, **patch_args):

    from shapely.geometry import Point
    import geopandas as gp

    pts = gp.GeoSeries([Point(x, y) for x, y in X])
    mp = pts.buffer(buf).unary_union
    return mp
    
    return PolygonPatch(mp, **patch_args)



