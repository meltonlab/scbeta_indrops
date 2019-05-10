import os, sys
import pandas as pd
import numpy as np
import scipy as sp

from scbeta_scrnaseq.process_indrops import SparsifiedIndropsAnalysis

        
class x1_x2_Base(SparsifiedIndropsAnalysis):

    data_dir = None
    sample = None

    log10_local_lib_size_threshold = 2.0

class x1_S3c(x1_x2_Base):
    data_dir=""
    sample="x1_S3c"
    flask='x1'
    stage="S3c"

    counts_input_files = [('x1_S3c_b1', 'x1_S3c_b1'), ('x1_S3c_b2', 'x1_S3c_b2')]

class x1_S4c(x1_x2_Base):
    data_dir=""
    sample="x1_S4c"
    flask='x1'
    stage="S4c"

    counts_input_files = [('x1_S4c_b1', 'x1_S4c_b1'), ('x1_S4c_b2', 'x1_S4c_b2')]

class x1_S5c(x1_x2_Base):
    data_dir=""
    sample="x1_S5c"
    flask='x1'
    stage="S5c"

    counts_input_files = [('x1_S5c_b1', 'x1_S5c_b1'), ('x1_S5c_b2', 'x1_S5c_b2')]

class x1_S6c(x1_x2_Base):
    data_dir=""
    sample="x1_S6c"
    flask='x1'
    stage="S6c"

    counts_input_files = [('x2_S3c_b1', 'x2_S3c_b1'), ('x2_S3c_b2', 'x2_S3c_b2')]

class x2_S3c(x1_x2_Base):
    data_dir=""
    sample="x2_S3c"
    flask='x2'
    stage="S3c"

    counts_input_files = [('x2_S3c_b1', 'x2_S3c_b1'), ('x2_S3c_b2', 'x2_S3c_b2')]

class x2_S4c(x1_x2_Base):
    data_dir=""
    sample="x2_S4c"
    flask='x2'
    stage="S4c"

    counts_input_files = [('x2_S4c_b1', 'x2_S4c_b1'), ('x2_S4c_b2', 'x2_S4c_b2')]

class x2_S5c(x1_x2_Base):
    data_dir=""
    sample="x2_S5c"
    flask='x2'
    stage="S5c"

    counts_input_files = [('x2_S5c_b1', 'x2_S5c_b1'), ('x2_S5c_b2', 'x2_S5c_b2')]

class x2_S6c(x1_x2_Base):
    data_dir=""
    sample="x2_S6c"
    flask='x2'
    stage="S6c"

    counts_input_files = [('x2_S6c_b1', 'x2_S6c_b1'), ('x2_S6c_b2', 'x2_S6c_b2')]
