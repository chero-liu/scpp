from anndata import AnnData
from typing import Union, Optional, Literal
from numpy.random import RandomState
import scanpy as sc
from typing import Union, Optional, Any, Mapping, Callable
from types import MappingProxyType
import numpy as np

_Method = Literal["umap", "gauss", "rapids"]
_MetricFn = Callable[[np.ndarray, np.ndarray], float]
# from sklearn.metrics.pairwise_distances.__doc__:
_MetricSparseCapable = Literal[
    "cityblock", "cosine", "euclidean", "l1", "l2", "manhattan"
]
_MetricScipySpatial = Literal[
    "braycurtis",
    "canberra",
    "chebyshev",
    "correlation",
    "dice",
    "hamming",
    "jaccard",
    "kulsinski",
    "mahalanobis",
    "minkowski",
    "rogerstanimoto",
    "russellrao",
    "seuclidean",
    "sokalmichener",
    "sokalsneath",
    "sqeuclidean",
    "yule",
]
_Metric = Union[_MetricSparseCapable, _MetricScipySpatial]


def neighbors(
    adata: AnnData,
    n_neighbors: int = 15,
    n_pcs: Optional[int] = None,
    use_rep: Optional[str] = None,
    knn: bool = True,
    random_state: Optional[Union[int, RandomState]] = 0,
    method: Optional[_Method] = "umap",
    metric: Union[_Metric, _MetricFn] = "euclidean",
    metric_kwds: Mapping[str, Any] = MappingProxyType({}),
    copy: bool = False,
    **kwargs,
) -> Optional[AnnData]:
    """\
    Compute a neighborhood graph of observations [McInnes18]_.
    The neighbor search efficiency of this heavily relies on UMAP [McInnes18]_,
    which also provides a method for estimating connectivities of data points -
    the connectivity of the manifold (`method=='umap'`). If `method=='gauss'`,
    connectivities are computed according to [Coifman05]_, in the adaption of
    [Haghverdi16]_.
    Parameters
    ----------
    adata
        Annotated data matrix.
    n_neighbors
        The size of local neighborhood (in terms of number of neighboring data
        points) used for manifold approximation. Larger values result in more
        global views of the manifold, while smaller values result in more local
        data being preserved. In general values should be in the range 2 to 100.
        If `knn` is `True`, number of nearest neighbors to be searched. If `knn`
        is `False`, a Gaussian kernel width is set to the distance of the
        `n_neighbors` neighbor.
    n_pcs
        Use this many PCs. If n_pcs==0 use .X if use_rep is None.
    use_rep
        Use the indicated representation. 'X' or any key for .obsm is valid.
        If None, the representation is chosen automatically: For .n_vars < 50, .X is used, otherwise ‘X_pca’ is used.
        If ‘X_pca’ is not present, it’s computed with default parameters.
    knn
        If `True`, use a hard threshold to restrict the number of neighbors to
        `n_neighbors`, that is, consider a knn graph. Otherwise, use a Gaussian
        Kernel to assign low weights to neighbors more distant than the
        `n_neighbors` nearest neighbor.
    random_state
        A numpy random seed.
    method
        Use ‘umap’ [McInnes18] or ‘gauss’ (Gauss kernel following [Coifman05] with adaptive width
        [Haghverdi16]) for computing connectivities. Use ‘rapids’ for the RAPIDS implementation of UMAP
        (experimental, GPU only).
    copy
        Return a copy instead of writing to adata.
    Returns
    -------
    Depending on `copy`, updates or returns `adata` with the following:
    **connectivities** : sparse matrix (`.uns['neighbors']`, dtype `float32`)
        Weighted adjacency matrix of the neighborhood graph of data
        points. Weights should be interpreted as connectivities.
    **distances** : sparse matrix (`.uns['neighbors']`, dtype `float32`)
        Instead of decaying weights, this stores distances for each pair of
        neighbors.
    """

    sc.pp.neighbors(
        adata,
        n_neighbors=n_neighbors,
        n_pcs=n_pcs,
        use_rep=use_rep,
        knn=knn,
        method=method,
        random_state=random_state,
        metric=metric,
        metric_kwds=metric_kwds,
        copy=copy,
        **kwargs,
    )

    print("Created k-Nearest-Neighbor graph in adata.uns['neighbors'] ")


def clustering(
    adata: AnnData,
    method: Union["louvain", "leiden"] = "leiden",
    resolution: float = 1.2,
    **kwargs,
):
    """
    Perform community detection on single-cell data using Louvain or Leiden algorithm.

    Parameters:
    - adata: AnnData
        An AnnData object containing single-cell data.
    - method: str, optional (default: 'louvain')
        The community detection algorithm to use. Options: 'louvain' or 'leiden'.
    - resolution: float, optional
        Resolution parameter for the community detection algorithm. If None, the default algorithm value is used.
    - louvain_key: str, optional (default: 'louvain')
        The key to store Louvain communities in adata.obs.
    - leiden_key: str, optional (default: 'leiden')
        The key to store Leiden communities in adata.obs.
    - **kwargs:
        Additional keyword arguments to be passed to the underlying community detection function.

    Returns:
    - None (in-place modification of adata)
    """
    if method == "louvain":
        sc.tl.louvain(adata, key_added=method, resolution=resolution, **kwargs)
    elif method == "leiden":
        sc.tl.leiden(adata, key_added=method, resolution=resolution, **kwargs)
    else:
        raise ValueError(f"Unsupported {method} method. Choose 'louvain' or 'leiden'.")
