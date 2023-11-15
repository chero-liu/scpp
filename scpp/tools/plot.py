from anndata import AnnData
import scanpy as sc
from scpp.tools.plotting.showmarker import save_CellMarkerPlot, PLOT_TYPE
from scpp.tools.plotting.embeddings import save_EmbeddingsPlot
from scpp.tools.plotting.normalised import plot_dual_histograms
from scpp.tools.plotting.qc import qcplot
from scpp.tools.rankgenes import RankGenes
from scpp.tools.utils import check_mkdir
from matplotlib import pyplot as plt


class Plot:
    def __init__(
        self,
        adata: AnnData,
        cluster: str,
        group: str,
        sample: str,
        outdir: str,
        prefix: str,
        dpi: int = 300,
    ):
        self.adata = adata
        self.cluster = cluster
        self.outdir = outdir
        self.prefix = prefix
        self.dpi = dpi
        self.group = group
        self.sample = sample

    def run(self):
        rankgenes_obj = RankGenes(self.adata)
        topn = rankgenes_obj.topn
        for plot_type in PLOT_TYPE.__args__:
            save_CellMarkerPlot(
                self.adata,
                topn,
                self.cluster,
                self.outdir,
                self.prefix,
                plot_type,
                self.dpi,
            )

        colors = [self.cluster, self.sample, self.group]
        for color in colors:
            save_EmbeddingsPlot(
                self.adata,
                color,
                "umap",
                self.outdir,
                self.prefix,
                self.dpi,
            )
            save_EmbeddingsPlot(
                self.adata,
                color,
                "tsne",
                self.outdir,
                self.prefix,
                self.dpi,
            )

        qcplot(self.adata, self.outdir, self.prefix, self.dpi)

        plot_dual_histograms(
            self.adata.layers["raw"].sum(1),
            self.adata.layers["normalised"].sum(1),
            "nCount_RNA",
            "Shifted logarithm",
            save_path=f"{self.outdir}/2.normalised/normalised",
        )

        check_mkdir(f"{self.outdir}/3.preprocessing/")
        sc.pl.highly_variable_genes(self.adata)
        plt.savefig(
            f"{self.outdir}/3.preprocessing/hvg.pdf",
            bbox_inches="tight",
        )
        plt.savefig(
            f"{self.outdir}/3.preprocessing/hvg.png",
            bbox_inches="tight",
            dpi=self.dpi,
        )
        plt.close()

        check_mkdir(f"{self.outdir}/4.pca/")
        sc.pl.pca_variance_ratio(self.adata)
        plt.savefig(
            f"{self.outdir}/4.pca/pca.pdf",
            bbox_inches="tight",
        )
        plt.savefig(
            f"{self.outdir}/4.pca/pca.png",
            bbox_inches="tight",
            dpi=self.dpi,
        )
        plt.close()
