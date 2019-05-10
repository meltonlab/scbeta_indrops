import os, sys
import pandas as pd
import numpy as np
import scipy as sp

from scbeta_scrnaseq.process_indrops import SparsifiedIndropsAnalysis

class Stage5_timecourse(SparsifiedIndropsAnalysis):

    sample = 'stage5_tc'

    maximal_mito_reads_fraction = 0.5
    log10_local_lib_size_threshold = 2.0
    cell_quality_thr = 750
    expected_fano_threshold = 1.2

    counts_input_files = [
        ('stg5diff1_S5d0_b1', 'stg5diff1_S5d0_b1'),
        ('stg5diff1_S5d1_b1', 'stg5diff1_S5d1_b1'),
        ('stg5diff1_S5d2_b1', 'stg5diff1_S5d2_b1'),
        ('stg5diff1_S5d3_b1', 'stg5diff1_S5d3_b1'),
        ('stg5diff1_S5d4_b1', 'stg5diff1_S5d4_b1'),
        ('stg5diff1_S5d5_b1', 'stg5diff1_S5d5_b1'),
        ('stg5diff1_S5d6_b1', 'stg5diff1_S5d6_b1'),
        ('stg5diff1_S5d7_b1', 'stg5diff1_S5d7_b1'),
        ('stg5diff2_S5d0_b1', 'stg5diff2_S5d0_b1'),
        ('stg5diff2_S5d1_b1', 'stg5diff2_S5d1_b1'),
        ('stg5diff2_S5d2_b1', 'stg5diff2_S5d2_b1'),
        ('stg5diff2_S5d3_b1', 'stg5diff2_S5d3_b1'),
        ('stg5diff2_S5d4_b1', 'stg5diff2_S5d4_b1'),
        ('stg5diff2_S5d5_b1', 'stg5diff2_S5d5_b1'),
        ('stg5diff2_S5d6_b1', 'stg5diff2_S5d6_b1'),
        ('stg5diff2_S5d7_b1', 'stg5diff2_S5d7_b1'),
        ('stg5diff1_S5d0_b2', 'stg5diff1_S5d0_b2'),
        ('stg5diff1_S5d1_b2', 'stg5diff1_S5d1_b2'),
        ('stg5diff1_S5d2_b2', 'stg5diff1_S5d2_b2'),
        ('stg5diff1_S5d3_b2', 'stg5diff1_S5d3_b2'),
        ('stg5diff1_S5d4_b2', 'stg5diff1_S5d4_b2'),
        ('stg5diff1_S5d5_b2', 'stg5diff1_S5d5_b2'),
        ('stg5diff1_S5d6_b2', 'stg5diff1_S5d6_b2'),
        ('stg5diff1_S5d7_b2', 'stg5diff1_S5d7_b2'),
        ('stg5diff2_S5d0_b2', 'stg5diff2_S5d0_b2'),
        ('stg5diff2_S5d1_b2', 'stg5diff2_S5d1_b2'),
        ('stg5diff2_S5d2_b2', 'stg5diff2_S5d2_b2'),
        ('stg5diff2_S5d3_b2', 'stg5diff2_S5d3_b2'),
        ('stg5diff2_S5d4_b2', 'stg5diff2_S5d4_b2'),
        ('stg5diff2_S5d5_b2', 'stg5diff2_S5d5_b2'),
        ('stg5diff2_S5d6_b2', 'stg5diff2_S5d6_b2'),
        ('stg5diff2_S5d7_b2', 'stg5diff2_S5d7_b2'),
        ]