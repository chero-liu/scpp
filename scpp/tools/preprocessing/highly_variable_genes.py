import scanpy as sc
from anndata import AnnData
from typing import Optional, Literal


def variable_filter(
    adata: AnnData,
    layer: Optional[str] = None,
    n_top_genes: Optional[int] = None,
    min_disp=0.5,
    max_disp=float("inf"),
    min_mean=0.0125,
    max_mean=3,
    span=0.3,
    n_bins=20,
    flavor: Optional[Literal["seurat", "cell_ranger", "seurat_v3"]] = "seurat",
    subset: bool = False,
    inplace: bool = True,
    batch_key: Optional[str] = None,
    check_values=True,
):
    """
    Identify highly variable genes in single-cell data using customizable parameters.

    Parameters:
    - adata: AnnData
        An AnnData object containing single-cell data.
    - layer: str, optional (default: None)
        The layer of adata to use. If None, the X or raw data is used.
    - n_top_genes: int, optional (default: None)
        Number of highly-variable genes to keep. If None, all are kept.
    - min_disp: float, optional (default: 0.5)
        Minimum dispersion for a gene to be considered highly variable.
    - max_disp: float, optional (default: float('inf'))
        Maximum dispersion for a gene to be considered highly variable.
    - min_mean: float, optional (default: 0.0125)
        Minimum mean expression for a gene to be considered highly variable.
    - max_mean: float, optional (default: 3)
        Maximum mean expression for a gene to be considered highly variable.
    - span: float, optional (default: 0.3)
        Fraction of the data range to use for binning.
    - n_bins: int, optional (default: 20)
        Number of bins.
    - flavor: {'seurat', 'cell_ranger'}, optional (default: 'seurat')
        Choose the flavor for computing the normalization. Use 'seurat' for Seurat-like normalization,
        'cell_ranger' for Cell Ranger-like normalization.
    - subset: bool, optional (default: False)
        Restrict the analysis to a subset of genes.
    - inplace: bool, optional (default: True)
        Whether to perform the operation in-place or return a new AnnData object.
    - batch_key: str, optional (default: None)
        Key in .obs that stores the batch origin.
    - check_values: bool, optional (default: True)
        Check the value range of each gene for the conditions min_disp and max_disp.

    Returns:
    - None or AnnData (depending on inplace parameter)
    """
    sc.pp.highly_variable_genes(
        adata,
        layer=layer,
        n_top_genes=n_top_genes,
        min_disp=min_disp,
        max_disp=max_disp,
        min_mean=min_mean,
        max_mean=max_mean,
        span=span,
        n_bins=n_bins,
        flavor=flavor,
        subset=subset,
        inplace=inplace,
        batch_key=batch_key,
        check_values=check_values,
    )

    return adata
