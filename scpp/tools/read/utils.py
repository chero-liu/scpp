from typing import Optional, Union
import os
import pandas as pd
from pathlib import Path
import scanpy as sc
from scipy import sparse
import sys

def getSingleDataFormat(filename):
    if os.path.isdir(filename):
        return "10x-format"
    elif filename.endswith(".h5ad"):
        return "h5ad-format"
    else:
        return "normal-format"


def read(
    filename: Union[str, Path],
    var_names: Optional[str] = "gene_symbols",
    make_unique: bool = True,
    cache: bool = False,
    delimiter: str = "\t",
    dtype: str = None,
    prefix: str = None,
):
    if dtype == "10x-format":
        adata = sc.read_10x_mtx(
            filename, var_names=var_names, cache=cache, make_unique=make_unique
        )
    elif dtype == "h5ad-format":
        adata = sc.read(filename)
        if "raw" in adata.layers:
            adata.X = adata.layers["raw"]
        else:
            sys.exit("No raw counts matrix can be found in Anndata object")
    elif dtype == "normal-format":
        try:
            adata = sc.read_csv(filename, delimiter=delimiter).T
        except ValueError:  # BD matrix
            adata = sc.read_csv(filename, delimiter=delimiter, first_column_names=True)
    else:
        raise ValueError("Not a valid matrix format, plz check.")

    gex_rows = list(map(lambda x: x.replace("_", "-"), adata.var.index))
    adata.var.index = gex_rows

    if prefix:
        gex_cols = list(map(lambda x: prefix + "_" + x, adata.obs.index))
        adata.obs.index = gex_cols

    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    return adata


def ensureSparseMatrix(adata):
    if not sparse.issparse(adata.X):
        adata.X = sparse.csr_matrix(adata.X)
    adata.layers["raw"] = adata.X
