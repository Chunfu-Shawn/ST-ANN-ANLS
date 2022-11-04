import sys
import scanpy as sc
import scipy.sparse
import pandas as pd


def read_sp(input, use_raw=True):
    adata = sc.read_h5ad(input)

    # use raw count
    if use_raw:
        adata = adata.raw.to_adata()

    # remove mito genes
    if 'mt' in adata.var.columns:
        adata = adata[:, ~adata.var['mt']].copy()

    adata.X = scipy.sparse.csr_matrix(adata.X)
    return adata


def read_sc(count, meta, celltype_key='cell_type', batch_key=None, categorical_covariate_key=None,
            continuous_covariate_key=None):
    count = sc.read_csv(count)
    count.X = scipy.sparse.csc_matrix(count.X)

    meta = pd.read_csv(meta, index_col=0)

    if all(count.obs.index == meta.index):
        count.obs = pd.concat([count.obs, meta], axis=1)

    if celltype_key not in count.obs.columns:
        raise KeyError('celltype_key ' + celltype_key + ' not found in ' + ', '.join(count.obs.columns) + '!')
    if batch_key is not None and batch_key not in count.obs.columns:
        raise KeyError('batch_key ' + batch_key + ' not found in ' + ', '.join(count.obs.columns) + '!')
    if categorical_covariate_key is not None and not all(x in count.obs.columns for x in categorical_covariate_key):
        raise KeyError('categorical_covariate_key ' + categorical_covariate_key + 'is not a subset of ' + ', '.join(
            count.obs.columns) + '!')
    if continuous_covariate_key is not None and not all(x in count.obs.columns for x in continuous_covariate_key):
        raise KeyError('continuous_covariate_key ' + continuous_covariate_key + 'is not a subset of ' + ', '.join(
            count.obs.columns) + '!')

    return count


def warn(message):
    sys.stderr.write(message + "\n")


def info(message):
    sys.stdout.write(message + "\n")
