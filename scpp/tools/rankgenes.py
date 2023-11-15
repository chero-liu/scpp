from anndata import AnnData
from typing import Optional, Literal
from scpp.tools.markergenes.rank_genes_groups import rank_genes_groups
from scpp.tools.markergenes.utils import write_diff_genes
import scanpy as sc


class RankGenes:
    """
    Find diff genes
    """

    def __init__(
        self,
        adata: AnnData,
        outdir: str = None,
        prefix: str = None,
        groupby: str = None,
        use_raw: Optional[bool] = False,
        pts: bool = True,
        method: Literal[
            "logreg", "t-test", "wilcoxon", "t-test_overestim_var"
        ] = "wilcoxon",
        logfc: float = 0.25,
        minpct: float = 0.1,
        top: int = 3,
    ):
        self.adata = adata
        self.outdir = outdir
        self.prefix = prefix
        self.method = method
        self.logfc = logfc
        self.minpct = minpct
        self.groupby = groupby
        self.use_raw = use_raw
        self.pts = pts
        self.top = top

    def run_RankGeneDefault(self):
        print("Performing rank genes ...")
        rank_genes_groups(
            self.adata,
            groupby=self.groupby,
            method=self.method,
            pts=self.pts,
            use_raw=self.use_raw,
            layer="normalised",
        )

        write_diff_genes(self.adata, self.outdir, self.prefix, self.logfc, self.minpct)
        return self.adata

    def run_RankGeneSpecified():
        pass

    @property
    def topn(self):
        if "rank_genes_groups" not in self.adata.uns:
            raise ValueError(
                "'rank_genes_groups' not in adata.uns, please rank genes first"
            )

        groups = self.adata.uns["rank_genes_groups"]["names"].dtype.names
        topn = {}

        for i in groups:
            de_i = sc.get.rank_genes_groups_df(self.adata, group=i)
            topn[i] = de_i["names"][: self.top].tolist()
        return topn
