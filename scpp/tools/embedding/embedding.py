from typing import Union, Optional
from anndata import AnnData
from numpy.random.mtrand import RandomState
import scanpy as sc


def run_tsne(
    adata: AnnData,
    n_pcs: Optional[int] = None,
    use_rep: Optional[str] = None,
    perplexity: Union[float, int] = 30,
    metric: str = "euclidean",
    early_exaggeration: Union[float, int] = 12,
    learning_rate: Union[float, int] = 1000,
    random_state: Optional[Union[int, RandomState]] = 0,
    copy: bool = False,
    output: Optional[str] = None,
) -> Optional[AnnData]:
    """\
    Wrap function scanpy.pp.tsne
    t-SNE [Maaten08] [Amir13] [Pedregosa11].
    t-distributed stochastic neighborhood embedding (tSNE) [Maaten08] has been proposed for visualizating single-cell
    data by [Amir13]. Here, by default, we use the implementation of scikit-learn [Pedregosa11].
    You can achieve a huge speedup and better convergence if you install Multicore-tSNE by [Ulyanov16],
    which will be automatically detected by Scanpy.
    Parameters
    ----------
    adata
        Annotated data matrix.
    n_pcs
        Use this many PCs. If n_pcs==0 use .X if use_rep is None.
    use_rep
        Use the indicated representation. 'X' or any key for .obsm is valid.
        If None, the representation is chosen automatically: For .n_vars < 50, .X is used,
        otherwise ‘X_pca’ is used. If ‘X_pca’ is not present, it’s computed with default parameters.
    perplexity
        The perplexity is related to the number of nearest neighbors that is used in other manifold learning algorithms.
        Larger datasets usually require a larger perplexity. Consider selecting a value between 5 and 50.
        The choice is not extremely critical since t-SNE is quite insensitive to this parameter.
    metric
        Distance metric calculate neighbors on.
    early_exaggeration
        Controls how tight natural clusters in the original space are in the embedded space and
        how much space will be between them. For larger values, the space between natural clusters will be
        larger in the embedded space. Again, the choice of this parameter is not very critical.
        If the cost function increases during initial optimization,
        the early exaggeration factor or the learning rate might be too high.
    learning_rate
        Note that the R-package “Rtsne” uses a default of 200. The learning rate can be a critical parameter.
        It should be between 100 and 1000. If the cost function increases during initial optimization,
        the early exaggeration factor or the learning rate might be too high.
        If the cost function gets stuck in a bad local minimum increasing the learning rate helps sometimes.
    random_state
        Change this to use different intial states for the optimization. If None, the initial state is not reproducible.
    copy
        Return a copy instead of writing to adata.
    output
        Save tsne embedding locations.
    Returns
    -------
    Depending on `copy`, returns or updates `adata` with the following fields.
    `X_tsne` : :class:`numpy.ndarray` (`adata.obsm`)
        Independent Component Analysis representation of data.

    """
    sc.tl.tsne(
        adata,
        n_pcs=n_pcs,
        use_rep=use_rep,
        perplexity=perplexity,
        metric=metric,
        early_exaggeration=early_exaggeration,
        learning_rate=learning_rate,
        random_state=random_state,
        copy=copy,
    )

    print("tSNE is done! Generated in adata.obsm['X_tsne'] and adata.uns['tsne']")


def run_umap(
    adata: AnnData,
    min_dist: float = 0.5,
    spread: float = 1.0,
    n_components: int = 2,
    maxiter: Optional[int] = None,
    alpha: float = 1.0,
    gamma: float = 1.0,
    negative_sample_rate: int = 5,
    init_pos: Optional[str] = "spectral",
    random_state: Optional[Union[int, RandomState]] = 0,
    a: Optional[float] = None,
    b: Optional[float] = None,
    copy: bool = False,
    method: Optional["str"] = "umap",
    output: Optional[str] = None,
) -> Optional[AnnData]:
    """\
    Wrap function scanpy.pp.umap
    Embed the neighborhood graph using UMAP [McInnes18]_.
    UMAP (Uniform Manifold Approximation and Projection) is a manifold learning
    technique suitable for visualizing high-dimensional data. Besides tending to
    be faster than tSNE, it optimizes the embedding such that it best reflects
    the topology of the data, which we represent throughout Scanpy using a
    neighborhood graph. tSNE, by contrast, optimizes the distribution of
    nearest-neighbor distances in the embedding such that these best match the
    distribution of distances in the high-dimensional space.  We use the
    implementation of `umap-learn <https://github.com/lmcinnes/umap>`__
    [McInnes18]_. For a few comparisons of UMAP with tSNE, see this `preprint
    <https://doi.org/10.1101/298430>`__.
    Parameters
    ----------
    adata
        Annotated data matrix.
    min_dist
        The effective minimum distance between embedded points.
        Smaller values will result in a more clustered/clumped embedding
        where nearby points on the manifold are drawn closer together,
        while larger values will result on a more even dispersal of points.
        The value should be set relative to the spread value,
        which determines the scale at which embedded points will be spread out.
        The default of in the umap-learn package is 0.1.
    spread
        The effective scale of embedded points.
        In combination with min_dist this determines how clustered/clumped the embedded points are.
    n_components
        The number of dimensions of the embedding.
    maxiter
        The number of iterations (epochs) of the optimization. Called n_epochs in the original UMAP.
    alpha
        The initial learning rate for the embedding optimization.
    gamma
        Weighting applied to negative samples in low dimensional embedding optimization.
        Values higher than one will result in greater weight being given to negative samples.
    negative_sample_rate
        The number of negative edge/1-simplex samples to use per positive edge/1-simplex sample in
        optimizing the low dimensional embedding.
    init_pos
        How to initialize the low dimensional embedding. Called init in the original UMAP. Options are:
        Any key for adata.obsm.
        ’paga’: positions from paga().
        ’spectral’: use a spectral embedding of the graph.
        ’random’: assign initial embedding positions at random.
        A numpy array of initial embedding positions.
    random_state
        If `int`, `random_state` is the seed used by the random number generator;
        If `RandomState`, `random_state` is the random number generator;
        If `None`, the random number generator is the `RandomState` instance used
        by `np.random`.
    a
        More specific parameters controlling the embedding.
        If None these values are set automatically as determined by min_dist and spread.
    b
        More specific parameters controlling the embedding.
        If None these values are set automatically as determined by min_dist and spread.
    method
        Use the original ‘umap’ implementation, or ‘rapids’ (experimental, GPU only)
    copy
        Return a copy instead of writing to adata.
    output
        Save umap embedding locations.
    Returns
    -------
    Depending on `copy`, returns or updates `adata` with the following fields.
    `X_umap` : :class:`numpy.ndarray` (`adata.obsm`)
        Independent Component Analysis representation of data.

    """

    sc.tl.umap(
        adata,
        min_dist=min_dist,
        spread=spread,
        n_components=n_components,
        maxiter=maxiter,
        alpha=alpha,
        gamma=gamma,
        negative_sample_rate=negative_sample_rate,
        init_pos=init_pos,
        random_state=random_state,
        a=a,
        b=b,
        copy=copy,
        method=method,
    )

    print("UMAP is done! Generated in adata.obsm['X_umap'] and adata.uns['umap']")
