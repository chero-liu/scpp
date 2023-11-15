from scpp.tools.utils import s_common, Step
from scpp.tools.rankgenes import RankGenes
from scpp.tools.plot import Plot
import unittest
import scanpy as sc
import omicverse as ov
from scpp.__init__ import ROOT_PATH


class AutoAnno(Step):
    def __init__(self, args):
        Step.__init__(self, args)
        self.adata = sc.read(args.h5ad)
        self.foldchange = args.foldchange
        self.pvalue = args.pvalue
        self.pySCSA_celltype = args.pySCSA_celltype
        self.target = args.target
        self.tissue = args.tissue
        self.groupby = "cluster"

        self.cal_logfc = args.cal_logfc
        self.method = args.method
        self.logfc = args.logfc
        self.minpct = args.minpct
        self.dpi = args.dpi

    def anno(self):
        scsa = ov.single.pySCSA(
            adata=self.adata,
            foldchange=self.foldchange,
            pvalue=self.pvalue,
            celltype=self.pySCSA_celltype,
            target=self.target,
            tissue=self.tissue,
            model_path=f"{ROOT_PATH}/integration/ref/pySCSA_2023_v2_plus.db",
        )
        self.adata.uns["log1p"]["base"] = 10
        scsa.cell_anno(clustertype="raw_cluster", cluster="all", rank_rep=False)
        scsa.cell_auto_anno(
            self.adata, clustertype="raw_cluster", key="scsa_celltype_cellmarker"
        )
        self.adata.obs.rename(
            columns={"scsa_celltype_cellmarker": "cluster"}, inplace=True
        )
        return self.adata

    def rankgenes(self, adata):
        rank_obj = RankGenes(
            adata=adata,
            outdir=self.outdir,
            prefix=self.prefix,
            species=self.species,
            method=self.method,
            cal_logfc=self.cal_logfc,
            logfc=self.logfc,
            minpct=self.minpct,
            groupby=self.groupby,
        )
        adata, topn, top2d = rank_obj.run()

        return adata, topn, top2d

    def plot(self, adata, cluster, topn, top2d):
        plot_obj = Plot(
            adata=adata,
            cluster=cluster,
            topn=topn,
            top2d=top2d,
            outdir=self.outdir,
            prefix=self.prefix,
            species=self.species,
            dpi=self.dpi,
        )

        plot_obj.run()

    def run(self):
        print("Start AutoAnno pipeline ...")

        adata = self.anno()

        adata, topn, top2d = self.rankgenes(adata)

        self.plot(adata, cluster=self.groupby, topn=topn, top2d=top2d)


def autoanno(args):
    with AutoAnno(args) as runner:
        runner.run()


def get_opts_autoanno(parser, sub_program=True):
    parser.add_argument(
        "--mtfilter",
        type=str,
        default=20,
        help="if default, automatically choose a mt threshold for all samples. Recommend",
    )
    parser.add_argument(
        "--mingene", type=float, default=200, help="minimum gene content cutoff"
    )
    parser.add_argument(
        "--maxgene", type=float, default=0.98, help="maximum gene content cutoff"
    )
    parser.add_argument(
        "--maxumi", type=float, default=0.98, help="maximum umi content cutoff"
    )
    parser.add_argument(
        "--minumi", type=float, default=0, help="minimum umi content cutoff"
    )
    parser.add_argument(
        "--mincell",
        type=int,
        default=5,
        help="filter genes if exists in less than an exact number of cells",
    )
    parser.add_argument(
        "--max_cell",
        type=int,
        default=500000,
        help="maximum cells cutoff to perform some analysis",
    )
    parser.add_argument(
        "--flavor", type=str, default="default", help="use all genes or hvg"
    )
    parser.add_argument(
        "--n_pcs",
        type=int,
        default=10,
        help="calculate principle components automatically or use an exact pc number",
    )
    parser.add_argument(
        "--batch_method",
        type=str,
        default="No",
        help="whether remove batch, and which method used",
    )
    parser.add_argument(
        "--clust_method",
        type=str,
        default="leiden",
        help="cluster method, louvain or leiden",
    )
    parser.add_argument(
        "--min_dist", type=str, default="default", help="a umap argument"
    )
    parser.add_argument(
        "--resolution",
        type=float,
        default=1.2,
        help="a parameter value controlling the coarseness of the clustering",
    )
    parser.add_argument(
        "--method", type=str, default="wilcoxon", help="rank genes method."
    )
    parser.add_argument(
        "--cal_logfc", type=str, default="seurat", help="method to calculate logfc."
    )
    parser.add_argument(
        "--logfc", type=float, default=0.25, help="logfoldchanges cutoff"
    )
    parser.add_argument("--minpct", type=float, default=0.1, help="min pct cutoff")
    parser.add_argument(
        "--legend_fontsize", type=int, default=10, help="legend font size"
    )
    parser.add_argument("--dpi", type=int, default=300, help="image dpi")
    parser.add_argument("--hvg_num", type=int, default=2000, help="hvg genes number")

    if sub_program:
        parser = s_common(parser)
        parser.add_argument(
            "--h5ad",
            type=str,
            help="Path to the input .h5ad file containing single-cell transcriptome data",
        )
        parser.add_argument(
            "--pySCSA_celltype",
            type=str,
            default="normal",
            choices=["normal", "cancer"],
            help="Specify the cell type for annotation (options: 'normal' or 'cancer')",
        )
        parser.add_argument(
            "--foldchange",
            type=float,
            default=1.5,
            help="Fold change threshold for identifying marker genes (default: 1.5). Higher values result in fewer markers.",
        )
        parser.add_argument(
            "--pvalue",
            type=float,
            default=0.01,
            help="P-value threshold for statistical significance in differential expression analysis (default: 0.01).",
        )
        parser.add_argument(
            "--target",
            type=str,
            default="cellmarker",
            help="Specify the target database for cell annotation (options: 'cellmarker', 'cancersea', 'panglaodb')",
        )

        parser.add_argument(
            "--tissue",
            type=str,
            default="All",
            help="Specify the tissue of interest for annotation (default: 'All')",
        )
    return parser


if __name__ == "__main__":
    unittest.main()
