import os, sys
import pandas as pd
import numpy as np
import scipy as sp

from scbeta_scrnaseq.process_indrops import SparsifiedIndropsAnalysis

class Stage6_timecourse(SparsifiedIndropsAnalysis):

    sample = 'stage6_tc'

    maximal_mito_reads_fraction = 0.5
    log10_local_lib_size_threshold = 2.0
    cell_quality_thr = 750
    expected_fano_threshold = 1.2

    counts_input_files = [
        ('diff1_wk0_b1', 'diff1_wk0_b1'),
        ('diff1_wk1_b1', 'diff1_wk1_b1'),
        ('diff1_wk2_b1', 'diff1_wk2_b1'),
        ('diff1_wk3_b1', 'diff1_wk3_b1'),
        ('diff1_wk4_b1', 'diff1_wk4_b1'),
        ('diff1_wk5_b1', 'diff1_wk5_b1'),
        ('diff2_wk0_b1', 'diff2_wk0_b1'),
        ('diff2_wk1_b1', 'diff2_wk1_b1'),
        ('diff2_wk2_b1', 'diff2_wk2_b1'),
        ('diff2_wk3_b1', 'diff2_wk3_b1'),
        ('diff2_wk4_b1', 'diff2_wk4_b1'),
        ('diff2_wk5_b1', 'diff2_wk5_b1'),
        ('diff3_wk0_b1', 'diff3_wk0_b1'),
        ('diff3_wk1_b1', 'diff3_wk1_b1'),
        ('diff3_wk2_b1', 'diff3_wk2_b1'),
        ('diff3_wk3_b1', 'diff3_wk3_b1'),
        ('diff3_wk4_b1', 'diff3_wk4_b1'),
        ('diff3_wk5_b1', 'diff3_wk5_b1'),
        ('diff1_wk0_b2', 'diff1_wk0_b2'),
        ('diff1_wk1_b2', 'diff1_wk1_b2'),
        ('diff1_wk2_b2', 'diff1_wk2_b2'),
        ('diff1_wk3_b2', 'diff1_wk3_b2'),
        ('diff1_wk4_b2', 'diff1_wk4_b2'),
        ('diff1_wk5_b2', 'diff1_wk5_b2'),
        ('diff2_wk1_b2', 'diff2_wk1_b2'),
        ('diff2_wk2_b2', 'diff2_wk2_b2'),
        ('diff2_wk3_b2', 'diff2_wk3_b2'),
        ('diff2_wk4_b2', 'diff2_wk4_b2'),
        ('diff2_wk5_b2', 'diff2_wk5_b2'),
        ('diff3_wk0_b2', 'diff3_wk0_b2'),
        ('diff3_wk1_b2', 'diff3_wk1_b2'),
        ('diff3_wk2_b2', 'diff3_wk2_b2'),
        ('diff3_wk3_b2', 'diff3_wk3_b2'),
        ('diff3_wk4_b2', 'diff3_wk4_b2'),
        ('diff3_wk5_b2', 'diff3_wk5_b2'),
        ]
