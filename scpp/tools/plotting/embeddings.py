import scanpy as sc
import matplotlib.pyplot as plt
from scpp.tools.utils import check_mkdir
from anndata import AnnData


def save_EmbeddingsPlot(
    adata: AnnData,
    color: str,
    plot_type: str,
    outdir: str,
    prefix: str,
    dpi: int,
):
    """
    Plot and save embeddings (UMAP or tSNE) for the specified color variable.

    Parameters:
    - adata (AnnData): Annotated Data object.
    - color (str): Variable used for coloring the plot.
    - prefix (str): Prefix for the saved plot files.
    - outdir (str): Output directory for saving the plot files.
    - plot_type (str, optional): Type of plot to generate. Valid options are 'umap' or 'tsne'.
    Defaults to 'umap'.
    - dpi (int, optional): Dots per inch for saving the plot. Defaults to 300.
    """
    # Create output directories if they don't exist
    save_path = f"{outdir}/6.{plot_type}"
    check_mkdir(save_path)

    # Define a mapping of plot types to Scanpy plotting functions
    plot_functions = {
        "umap": sc.pl.umap,
        "tsne": sc.pl.tsne,
    }

    # Check if the specified plot_type is valid
    if plot_type not in plot_functions:
        raise ValueError("Invalid plot_type. Choose either 'umap', or 'tsne'.")

    # Call the appropriate plotting function
    plot_functions[plot_type](adata, color=color)

    plt.gca().set_box_aspect(1)

    # Save plots in both PDF and PNG formats
    plt.savefig(
        f"{save_path}/{prefix}_{color}_{plot_type}.pdf",
        bbox_inches="tight",
    )
    plt.savefig(
        f"{save_path}/{prefix}_{color}_{plot_type}.png",
        bbox_inches="tight",
        dpi=dpi,
    )

    plt.close()
