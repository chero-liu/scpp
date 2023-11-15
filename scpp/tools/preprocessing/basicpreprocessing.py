from anndata import AnnData
from typing import Optional, List, Union
import scanpy as sc
import numpy as np
import pandas as pd


def _filter(
    adata: AnnData,
    mingenes: Optional[float] = None,
    minumis: Optional[float] = None,
    maxgenes: Optional[float] = None,
    maxumis: Optional[float] = None,
    mincells: Optional[int] = None,
    save_path: str = None,
    dpi: Optional[int] = None,
) -> AnnData:
    sc.pp.filter_cells(adata, min_genes=mingenes)
    sc.pp.filter_cells(adata, min_counts=minumis)

    maxgenes_new = (
        maxgenes
        if maxgenes > 1
        else np.percentile(adata.obs["n_genes"], 100 * maxgenes) + 1
    )
    maxumis_new = (
        maxumis
        if maxumis > 1
        else np.percentile(adata.obs["n_counts"], 100 * maxumis) + 1
    )

    sc.pp.filter_cells(adata, max_genes=maxgenes_new)
    sc.pp.filter_cells(adata, max_counts=maxumis_new)

    sc.pp.filter_genes(adata, min_cells=mincells)

    return adata


def findMtThreshold(
    adata: AnnData,
    species: str,
    mt_thresholds: List[float],
    mtfilter: str = "default",
) -> float:
    """
    Calculate mt threshold for each sample (generic version)
    """
    if species in ["human", "mouse"]:
        adata.var["mt"] = adata.var_names.str.contains("^[Mm][Tt]-")
    else:
        adata.var["mt"] = False
        print("WARNING: can not find any mt genes.")

    mt_df = sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=False
    )[0]["pct_counts_mt"]

    mt_count = []
    for ts in mt_thresholds:
        mt_count.append(len(mt_df[mt_df > ts]) / len(mt_df))

    mt_loc = mt_count.index(min(mt_count, key=lambda x: abs(x - 0.05)))
    if mt_count[mt_loc] > 0.3:
        print("WARNING: Filter too much cells according to MT threshold")

    if mtfilter != "default":
        return float(mtfilter)
    else:
        mtfilter = mt_thresholds[mt_loc]
        return float(mtfilter)


def calculate_mt_common(mtfilter: str, mt_list: List[int]) -> int:
    if mtfilter != "default":
        try:
            mt_common = int(mtfilter)
        except (ValueError, TypeError):
            print("Warning: mtfilter is not a valid integer. Using the default value.")
            mt_common = "default"  # You can set a default value or handle it as needed
    else:
        df = pd.DataFrame(mt_list)
        mt_common = df.mode().sort_values(by=0, ascending=False)[0].iloc[0]

    return mt_common


def normalize_total(
    adata: AnnData,
    target_sum: Optional[float] = None,
    exclude_highly_expressed: bool = False,
    max_fraction: float = 0.05,
    key_added: Optional[str] = None,
    layer: Optional[str] = None,
    layers: Optional[List[str]] = None,
    layer_norm: Optional[str] = None,
    inplace: bool = True,
    copy: bool = False,
) -> None:
    """
    Normalize total counts per cell.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    target_sum : float, optional (default: None)
        If None, after normalization, each observation (cell) has a total count equal to the median of the
        *counts_per_cell* before normalization.
    exclude_highly_expressed : bool, optional (default: False)
        Exclude genes that are expressed above `max_fraction` in a cell.
    max_fraction : float, optional (default: 0.05)
        If `exclude_highly_expressed` is True, exclude genes that are expressed above `max_fraction` in a cell.
    key_added : str, optional (default: None)
        Name of the field in `adata.obs` where the total counts per cell are stored.
    layer : str, optional (default: None)
        Name of the layer where the normalized data should be stored. By default, the normalized data is stored in
        `adata.X`.
    layers : List[str], optional (default: None)
        Names of the layers where the normalized data should be stored. If `layers` is not None, `layer` should be None.
    layer_norm : str, optional (default: None)
        If not None, renormalize the values of the layer such that the values sum to `layer_norm`.
    inplace : bool, optional (default: True)
        Whether to place the result back into the `AnnData` object or return it.
    copy : bool, optional (default: False)
        If True, return a copy of `adata` instead of updating it in place.

    Returns
    -------
    None
    """
    sc.pp.normalize_total(
        adata,
        target_sum=target_sum,
        exclude_highly_expressed=exclude_highly_expressed,
        max_fraction=max_fraction,
        key_added=key_added,
        layer=layer,
        layers=layers,
        layer_norm=layer_norm,
        inplace=inplace,
        copy=copy,
    )
    print("Normalization step is finished in adata.X")


def log1p(
    adata: AnnData,
    base: Optional[float] = None,
    copy: bool = False,
    chunked: Optional[bool] = None,
    chunk_size: Optional[int] = None,
    layer: Optional[str] = None,
    obsm: Optional[str] = None,
) -> None:
    """
    Logarithmize the data matrix.

    Computes X = log(X + 1), where log denotes the natural logarithm unless a different base is given.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    base : float, optional (default: None)
        Base of the logarithm. Natural logarithm is used by default.
    copy : bool, optional (default: False)
        If an AnnData is passed, determines whether a copy is returned.
    chunked : bool, optional (default: None)
        Process the data matrix in chunks, which will save memory. Applies only to AnnData.
    chunk_size : int, optional (default: None)
        n_obs of the chunks to process the data in.
    layer : str, optional (default: None)
        Entry of layers to transform.
    obsm : str, optional (default: None)
        Entry of obsm to transform.

    Returns
    -------
    None
    """
    sc.pp.log1p(
        adata,
        base=base,
        copy=copy,
        chunked=chunked,
        chunk_size=chunk_size,
        layer=layer,
        obsm=obsm,
    )
    print("Log transformation step is finished in adata.X")


def normalize(adata: AnnData, target_sum: float = 1e4) -> None:
    """
    Normalize the data matrix by total counts and logarithmize the result.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    target_sum : float, optional (default: 1e4)
        If None, after normalization, each observation (cell) has a total count equal to the median of the
        *counts_per_cell* before normalization.

    Returns
    -------
    None
    """
    normalize_total(adata, target_sum=target_sum)
    log1p(adata)


def scale(
    adata: Union[AnnData],
    zero_center: bool = True,
    max_value: Optional[float] = None,
    copy: bool = False,
) -> Optional[AnnData]:
    """\
    Wrap function of scanpy.pp.scale

    Scale data to unit variance and zero mean.
    .. note::
        Variables (genes) that do not display any variation (are constant across
        all observations) are retained and set to 0 during this operation. In
        the future, they might be set to NaNs.
    Parameters
    ----------
    data
        The (annotated) data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    zero_center
        If `False`, omit zero-centering variables, which allows to handle sparse
        input efficiently.
    max_value
        Clip (truncate) to this value after scaling. If `None`, do not clip.
    copy
        If an :class:`~anndata.AnnData` is passed,
        determines whether a copy is returned.
    Returns
    -------
    Depending on `copy` returns or updates `adata` with a scaled `adata.X`.
    """

    sc.pp.scale(adata, zero_center=zero_center, max_value=max_value, copy=copy)

    print("Scale step is finished in adata.X")
