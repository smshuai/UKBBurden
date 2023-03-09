import sys
import numpy as np
import h5py
import pandas as pd
from joblib import dump, load
from sklearn.decomposition import IncrementalPCA


def read_slice(chunks, key, idx_start, idx_end):
    """ Read column slice (idx_start:idx_end) in to a numpy array.
    """
    arr = list()
    for chunk in chunks:
        with h5py.File(chunk, 'r') as f:
            arr.append(f[key][:,idx_start:idx_end])
    arr = np.vstack(arr).T
    return(arr)


## arguments
key = sys.argv[1] # key for hdf5
out1 = sys.argv[2] # joblib name for IPCA
out2 = sys.argv[3] # csv name for top PCs
x = int(sys.argv[4]) # run in x batches

chunks = ['/nfs/research/birney/users/shimin/ukbb/wes/200k/pvcf/7_pca_chunk/chunk_{0:02d}.hdf5'.format(i) for i in range(20)]

ipca = IncrementalPCA(n_components=50, whiten=True)

idx = np.linspace(0, 200643, x, dtype=np.int32)
for idx_start, idx_end in zip(idx[:-1], idx[1:]):
    print(idx_start, idx_end)
    dat = read_slice(chunks, key, idx_start, idx_end)
    ipca.partial_fit(dat)

dump(ipca, out1)

top_pcs = list()
for idx_start, idx_end in zip(idx[:-1], idx[1:]):
    print(idx_start, idx_end)
    dat = read_slice(chunks, key, idx_start, idx_end)
    top_pcs.append(ipca.transform(dat))

top_pcs = np.vstack(top_pcs)

sample_names = np.loadtxt("/nfs/research/birney/users/shimin/ukbb/wes/200k/pvcf/sample_names.txt", dtype=str)
top_pcs = pd.DataFrame(top_pcs, columns=['PC' + str(i) for i in range(1, 51)], index=sample_names)
top_pcs.to_csv(out2)
