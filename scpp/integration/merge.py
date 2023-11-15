from scpp.tools.utils import s_common, Step
from scpp.tools.read.create import Create
from scpp.tools.cluster import Cluster
from scpp.tools.rankgenes import RankGenes
from scpp.tools.plot import Plot
import unittest
import matplotlib

matplotlib.use("Agg")


class Merge(Step):
    def __init__(self, args):
        Step.__init__(self, args)
        self.info = args.info
        # self.outdir = args.outdir
        # self.prefi = args.prefix
        self.mtfilter = args.mtfilter

        self.mingene = args.mingene
        self.minumi = args.minumi
        self.maxgene = args.maxgene
        self.maxumi = args.maxumi
        self.mincell = args.mincell

        self.max_cell = args.max_cell
        self.flavor = args.flavor
        self.n_pcs = args.n_pcs
        self.batch_method = args.batch_method

        self.clust_method = args.clust_method
        self.min_dist = args.min_dist
        self.resolution = args.resolution

        self.method = args.method
        self.cal_logfc = args.cal_logfc
        self.logfc = args.logfc
        self.minpct = args.minpct
        self.hvg_num = args.hvg_num

        self.legend_fontsize = args.legend_fontsize
        self.dpi = args.dpi
        self.groupby = "raw_cluster"

    def create(self):
        create_obj = Create(
            info=self.info,
            outdir=self.outdir,
            prefix=self.prefix,
            species=self.species,
            mtfilter=self.mtfilter,
            mingene=self.mingene,
            minumi=self.minumi,
            maxumi=self.maxumi,
            maxgene=self.maxgene,
            mincell=self.mincell,
            dpi=self.dpi,
        )
        adata = create_obj.run()

        return adata

    def cluster(self, adata):
        cluster_obj = Cluster(
            adata=adata,
            hvg_num=self.hvg_num,
            n_pcs=self.n_pcs,
            batch_method=self.batch_method,
            clust_method=self.clust_method,
            min_dist=self.min_dist,
            resolution=self.resolution,
        )
        adata = cluster_obj.run()

        return adata

    def rankgenes(self, adata):
        rank_obj = RankGenes(
            adata=adata,
            outdir=self.outdir,
            prefix=self.prefix,
            method=self.method,
            logfc=self.logfc,
            minpct=self.minpct,
            groupby=self.groupby,
        )
        adata = rank_obj.run_RankGeneDefault()

        return adata

    def plot(self, adata):
        plot_obj = Plot(
            adata,
            self.groupby,
            "gname",
            "sample",
            self.outdir,
            self.prefix,
        )
        plot_obj.run()

    def run(self):
        print("Start Merge pipeline ...")

        adata = self.create()

        adata = self.cluster(adata)

        adata = self.rankgenes(adata)

        adata.write_h5ad(f"{self.outdir}/{self.prefix}_PRO_diff.h5ad")

        self.plot(adata)


def merge(args):
    with Merge(args) as runner:
        runner.run()


def get_opts_merge(parser, sub_program=True):
    parser.add_argument(
        "--mtfilter",
        type=str,
        default="default",
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
            "--info",
            help="a project description config file, include rawname, samplename, groupname, etc",
        )
    return parser


if __name__ == "__main__":
    unittest.main()
