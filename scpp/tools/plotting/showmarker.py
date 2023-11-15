import scanpy as sc
import matplotlib.pyplot as plt
from anndata import AnnData
from typing import Dict, Any, Literal, Optional
from scpp.tools.utils import check_mkdir

PLOT_TYPE = Literal[
    "matrixplot",
    "heatmap",
    "dotplot",
    "stacked_violin",
    "rank_genes_groups_dotplot",
]


def save_CellMarkerPlot(
    adata: AnnData,
    top_genes: Dict[str, Any],
    groupby: str,
    outdir: str,
    prefix: str,
    plot_type: Optional[PLOT_TYPE],
    dpi: int = 300,
    dendrogram: bool = False,
    **kwargs,
) -> None:
    """
    Generate and save a specified type of cell marker plot.

    Parameters:
    - adata (AnnData): Annotated data matrix.
      An AnnData object containing single-cell RNA-seq data.
    - top_genes (Dict[str, Any]): Top genes for the plot.
      A dictionary containing information about the top genes, e.g., {gene_key: gene_values}.
    - groupby (str): The column in adata.obs to group by.
      The variable by which the data will be grouped.
    - outdir (str): Output directory path.
      The path to the directory where the output files will be saved.
    - prefix (str): Prefix for the output files.
      A string used as a prefix for the names of the output files.
    - plot_type (str): Type of plot to generate ('matrixplot', 'heatmap', 'dotplot', 'stacked_violin', 'rank_genes_groups_dotplot').
      The specific type of plot to generate.
    - dpi (int, optional): Dots per inch of the saved figure. Default is 300.
      The resolution of the saved figure.
    - **kwargs: Additional keyword arguments specific to the chosen plot type.

    Returns:
    None
    """
    if plot_type == "matrixplot":
        sc.pl.matrixplot(
            adata,
            top_genes,
            groupby,
            dendrogram=dendrogram,
            log=False,
            use_raw=False,
            cmap="RdBu_r",
            **kwargs,
        )
    elif plot_type == "heatmap":
        sc.pl.heatmap(
            adata,
            top_genes,
            groupby=groupby,
            dendrogram=dendrogram,
            swap_axes=True,
            show_gene_labels=True,
            use_raw=False,
            standard_scale="var",
            **kwargs,
        )
    elif plot_type == "dotplot":
        sc.pl.dotplot(
            adata,
            top_genes,
            groupby,
            dendrogram=dendrogram,
            use_raw=False,
            swap_axes=False,
            **kwargs,
        )
    elif plot_type == "stacked_violin":
        sc.pl.stacked_violin(
            adata,
            top_genes,
            groupby=groupby,
            dendrogram=dendrogram,
            use_raw=False,
            **kwargs,
        )
    elif plot_type == "rank_genes_groups_dotplot":
        sc.pl.rank_genes_groups_dotplot(
            adata,
            n_genes=3,
            dendrogram=dendrogram,
            values_to_plot="logfoldchanges",
            vmax=7,
            vmin=-7,
            cmap="RdBu_r",
            **kwargs,
        )

    check_mkdir(f"{outdir}/5.cellmarker/")
    plt.savefig(
        f"{outdir}/5.cellmarker/{prefix}_{plot_type}.pdf",
        bbox_inches="tight",
    )
    plt.savefig(
        f"{outdir}/5.cellmarker/{prefix}_{plot_type}.png",
        bbox_inches="tight",
        dpi=dpi,
    )
    plt.close()
