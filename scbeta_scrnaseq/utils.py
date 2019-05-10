import abc
import os
import pandas as pd
import numpy as np
import scipy as sp
import pickle

import palettable

def preload_ds_vals(ds):
    ds.vals = sp.sparse.csr_matrix(ds.layers[""][:, :])


def save_df(obj, filename):
    np.savez_compressed(filename, data=obj.values, index=obj.index.values, columns=obj.columns.values)
    
def load_df(filename):
    with np.load(filename) as f:
        obj = pd.DataFrame(**f)
    return obj


class SparseCSRDataFrame():
    def __init__(self, values, columns=None, index=None):
        self.values = values
        self.columns = columns
        self.index = index
        self.shape = values.shape

    def save(self, filename):
        np.savez_compressed(filename, sparse_data=self.values.data, sparse_indices=self.values.indices, sparse_indptr=self.values.indptr,
            index=self.index, columns=self.columns)

def sparse_load_df(filename):
    with np.load(filename) as f:
        obj = utils.SparseCSRDataFrame(sp.sparse.csr_matrix((f['sparse_data'], f['sparse_indices'], f['sparse_indptr']), shape=(len(f['index']), len(f['columns']))),
            index=f['index'], columns=f['columns'])
    return obj

class SavedDataFrame():
    
    def __init__(self, fn):    
        self.fn = fn
        
    def __getitem__(self, prefix):
        obj = load_df(self.fn%prefix)
        if obj.shape[1] == 1:
            return obj.iloc[:, 0]
        else:
            return obj

class SavedArray():
    
    def __init__(self, fn):    
        self.fn = fn
        self.cache = {}
        
    def __getitem__(self, prefix):
        obj = np.load(self.fn%prefix)
        return obj

class CachedAnalysis(metaclass=abc.ABCMeta):
    """
    Base class for cached analysis data. 

    Needs to define a `cached_attributes` dict of (attribute name, caching filetype) pairs
    And a `process` method that generates all the cached data.

    Upon instantiation, this class will check that all cached items exist and if any are missing,
    will run the analysis (by calling `process()`) and cache new items. 

    Cached data is not loaded until the first time it is needed (to save on memory usage).

    A list of available data attributes is present in `self.attrs` and these call be accessed using `self.[attr_name]`

    """

    def __getattr__(self, name):
        if name in self.cached_attributes.keys():
            cached_name = '_%s' % name
            if cached_name not in self.__dict__.keys():
                setattr(self, cached_name, self.read_attr(name))
            return self.__dict__[cached_name]

        else:
            raise AttributeError(name)
                

    def __init__(self, output_and_cache_dir):

        self.root_path = output_and_cache_dir
        os.makedirs(self.root_path, exist_ok=True)

        self.attrs = self.cached_attributes.keys()

        # Check that all cached files exist. If any are missing, re-run the full analysis.
        cache_incomplete = False
        for attr_name in self.cached_attributes.keys():
            attr_filename = self.cached_attr_filename(attr_name)
            if not os.path.isfile(attr_filename):
                if not cache_incomplete:
                    print('Missing (at least) file %s.' % attr_filename)
                cache_incomplete = True

        if cache_incomplete:
            print('Cache invalid for %s.' % (self.__class__.__name__))
            self.update_cache()

    def update_cache(self):
        print('Recomputing cache for %s.' % (self.__class__.__name__))
        self.save_object_dict(self.process())

    @abc.abstractmethod
    def process(self):
        pass

    @abc.abstractmethod
    def cached_attributes(self):
        pass

    def cached_attr_filename(self, attr_name):
        filename = '%s.%s.%s' % (self.__class__.__name__, attr_name, self.cached_attributes[attr_name])
        return os.path.join(self.root_path, filename)

    def save_object_dict(self, attr_dict):
        for attr_name, attr_type in self.cached_attributes.items():
            print('Saving %s to cache' % attr_name)
            self.save_attr(attr_dict[attr_name], attr_name)


    def save_attr(self, obj, attr_name):
        obj_type = self.cached_attributes[attr_name]
        filename = self.cached_attr_filename(attr_name)

        if obj_type == 'df.csv.gz':
            obj.to_csv(filename, compression = 'gzip')
        elif obj_type == 'df.h5':
            obj.to_hdf(filename, key=attr_name, complevel=1, complib='zlib')
        elif obj_type == 'df.npz':
            np.savez_compressed(filename, data=obj.values, index=obj.index.values, columns=obj.columns.values)
        elif obj_type == 'csr.df.npz':
            np.savez_compressed(filename, sparse_data=obj.values.data, sparse_indices=obj.values.indices, sparse_indptr=obj.values.indptr,
                index=obj.index, columns=obj.columns)
        # elif obj_type == 'df.npy':
        #     np.save(filename, data=obj.values, index=obj.index.values, columns=obj.columns.values)
        elif obj_type == 'series.npz':
            np.savez_compressed(filename, data=obj.values, index=obj.index.values)
        elif obj_type == 'npy':
            np.save(filename, obj)
        elif obj_type == 'pickle':
            with open(filename, 'wb') as f:
                pickle.dump(obj, f)
        elif obj_type == 'pn.npz':
            np.savez_compressed(filename, data=obj.values, items=obj.items.values,
                minor_axis=obj.minor_axis.values, major_axis=obj.major_axis.values)
        else: 
            raise('Unsupported object type')


    def read_attr(self, attr_name):
        obj_type = self.cached_attributes[attr_name]
        filename = self.cached_attr_filename(attr_name)
        if obj_type == 'df.csv.gz':
            obj = pd.read_csv(filename, index_col = 0, compression = 'gzip')
        elif obj_type == 'df.h5':
            obj = pd.read_hdf(filename, key=attr_name) 
        elif obj_type == 'df.npz':
            with np.load(filename) as f:
                obj = pd.DataFrame(**f)
        elif obj_type == 'csr.df.npz':
            with np.load(filename) as f:
                obj = SparseCSRDataFrame(sp.sparse.csr_matrix((f['sparse_data'], f['sparse_indices'], f['sparse_indptr']), shape=(len(f['index']), len(f['columns']))),
                    index=f['index'], columns=f['columns'])
        # elif obj_type == 'df.npy':
        #     with np.load(filename) as f:
        #         obj = pd.DataFrame(**f)
        elif obj_type == 'series.npz':
            with np.load(filename) as f:
                obj = pd.Series(**f)
        elif obj_type == 'npy':
            obj = np.load(filename)
        elif obj_type == 'pickle':
            with open(filename, 'rb') as f:
                obj = pickle.load(f)
        elif obj_type == 'pn.npz':
            with np.load(filename) as f:
                obj = pd.Panel(**f)
        else: 
            raise('Unsupported object type')
        return obj


    def __repr__(self):
        rep = """<%s - Cached Analysis:\n""" % self.__class__.__name__
        for key in self.cached_attributes.keys():
            rep += "  - %s\n" % key
        rep += ">"
        return rep

def two_array_corr(A,B):
    # Get number of rows in either A or B
    N = B.shape[0]

    # Store columnw-wise in A and B, as they would be used at few places
    sA = A.sum(0)
    sB = B.sum(0)

    # Basically there are four parts in the formula. We would compute them one-by-one
    p1 = N*np.einsum('ij,ik->kj',A,B)
    p2 = sA*sB[:,None]
    p3 = N*((B**2).sum(0)) - (sB**2)
    p4 = N*((A**2).sum(0)) - (sA**2)

    # Finally compute Pearson Correlation Coefficient as 2D array 
    pcorr = ((p1 - p2)/np.sqrt(p4*p3[:,None]))
    return pcorr


def where_idx(ser):
    return ser[ser].index


from collections import OrderedDict
def combine_rows(dfs, rows, na=0):
    new_df = OrderedDict({})
    for inrow in rows:
        if len(inrow) == 2:
            df, row = inrow
            label = f"{df} {row}"
        elif len(inrow) == 3:
            df, row, label = inrow
            
        new_df[label] = dfs[df].loc[row]
    new_df = pd.DataFrame(new_df).fillna(na).T
    return new_df