import os, sys
import pandas as pd
import numpy as np
import scipy as sp

from scbeta_scrnaseq.process_indrops import SparsifiedIndropsAnalysis

        
class HumanIslets(SparsifiedIndropsAnalysis):

    data_dir = None
    sample = None

    log10_local_lib_size_threshold = 2.0

    counts_input_files = [
        ('human1_lib1', 'human1_lib1'),
        ('human1_lib2', 'human1_lib2'),
        ('human1_lib3', 'human1_lib3'),
        ('human2_lib1', 'human2_lib1'),
        ('human2_lib2', 'human2_lib2'),
        ('human2_lib3', 'human2_lib3'),
        ('human3_lib1', 'human3_lib1'),
        ('human3_lib2', 'human3_lib2'),
        ('human3_lib3', 'human3_lib3'),
        ('human3_lib4', 'human3_lib4'),
        ('human4_lib1', 'human4_lib1'),
        ('human4_lib2', 'human4_lib2'),
        ]
        
class PSCsComparison(SparsifiedIndropsAnalysis):

    data_dir = None
    sample = None

    log10_local_lib_size_threshold = 2.0

    counts_input_files = [
        ('hues8_v4_b1', 'hues8_v4_b1'),
        ('hues8_v4_b2', 'hues8_v4_b2'),
        ('ips101631_v4_b1', 'ips101631_v4_b1'),
        ('ips101631_v4_b2', 'ips101631_v4_b2'),
        ('hues8_x3_b2', 'hues8_x3_b2'),
        ]


class ReaggComparison(SparsifiedIndropsAnalysis):

    data_dir = None
    sample = None
    maximal_mito_reads_fraction = 1.0
    log10_local_lib_size_threshold = 2.0

    counts_input_files = [
        ('Native_b1', 'Native_b1'),
        ('Native_b2', 'Native_b2'),
        ('Reagg_b1', 'Reagg_b1'),
        ('Reagg_b2', 'Reagg_b2'),
        ]

