import sys
import numpy as np
import h5py
import pandas as pd
from joblib import dump, load
from sklearn.decomposition import IncrementalPCA

def read_slice(chunks, idx_start, idx_end):
    """ Read column slice (idx_start:idx_end) in to a numpy array.
    """
    arr = list()
    for chunk in chunks:
        with h5py.File(chunk, 'r') as f:
            arr.append(f['rare'][:,idx_start:idx_end])
    arr = np.vstack(arr).T
    return(arr)

chunks = ['pca_chunk_{0:02d}.hdf5'.format(i) for i in range(60)]

ipca = IncrementalPCA(n_components=50, whiten=True)

idx = np.linspace(0, 469835, 400, dtype=np.int32)
for idx_start, idx_end in zip(idx[:-1], idx[1:]):
    print(idx_start, idx_end)
    dat = read_slice(chunks, idx_start, idx_end)
    ipca.partial_fit(dat)
    dump(ipca, "./ipca_model.joblib")

top_pcs = list()
for idx_start, idx_end in zip(idx[:-1], idx[1:]):
    print(idx_start, idx_end)
    dat = read_slice(chunks, idx_start, idx_end)
    top_pcs.append(ipca.transform(dat))

top_pcs = np.vstack(top_pcs)