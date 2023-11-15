import pandas as pd
from anndata import AnnData


def transform_cluster_labels(adata: AnnData, cluster_method: str):
    """
    Transform and rename cluster labels in AnnData object.

    Parameters:
    - adata: AnnData
        An AnnData object containing single-cell data.
    - cluster_method: str
        The key indicating the original cluster labels in adata.obs.

    Returns:
    - None (in-place modification of adata)
    """
    # Copy the original cluster labels to a new key
    adata.obs["raw_cluster"] = adata.obs[cluster_method]

    # Convert cluster labels to numeric and find the maximum
    max_num = max(pd.to_numeric(adata.obs["raw_cluster"]))

    # Rename categories in "raw_cluster" to consecutive integers starting from 1
    adata.rename_categories("raw_cluster", list(map(str, list(range(1, max_num + 2)))))


def umap_argue(adata) -> int:
    celln = len(adata.obs.index)
    if celln < 1000:
        min_dist = 0.5
    elif (celln >= 1000) and (celln < 100000):
        min_dist = 0.3
    elif (celln >= 100000) and (celln < 300000):
        min_dist = 0.2
    else:
        min_dist = 0.1
    return min_dist
