from typing import Literal, Optional,Union,Iterable
from anndata import AnnData
import scanpy as sc

def rank_genes_groups(
    adata: AnnData,
    groupby: str = None,
    use_raw: Optional[bool] = None,
    layer: Optional[str] = None,
    groups: Union[Literal["all"], Iterable[str]] = "all",
    reference: str = "rest",
    n_genes: Optional[int] = None,
    method: Literal["logreg", "t-test", "wilcoxon", "t-test_overestim_var"] = "wilcoxon",
    corr_method: Literal["benjamini-hochberg", "bonferroni"] = "benjamini-hochberg",
    tie_correct: bool = False,
    rankby_abs: bool = False,
    pts: bool = False,
    key_added: Optional[str] = None
) -> Optional[AnnData]:
    """\
    Wrap function scanpy.tl.rank_genes_groups
    Rank genes for characterizing groups.
    Expects logarithmized data.
    """

    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        use_raw=use_raw,
        layer=layer,
        groups=groups,
        reference=reference,
        n_genes=n_genes,
        method=method,
        corr_method=corr_method,
        tie_correct=tie_correct,
        rankby_abs=rankby_abs,
        pts=pts,
        key_added=key_added,
    )