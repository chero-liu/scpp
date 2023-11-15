import scanpy as sc
from matplotlib import pyplot as plt
from anndata import AnnData
from scpp.tools.utils import check_mkdir


def qcplot(
    adata: AnnData,
    outdir: str,
    prefix: str,
    dpi: int = 300,
):
    """
    Generate quality control (QC) plots for RNA-seq data.

    This function creates scatter plots and a multi-panel violin plot to visualize quality control metrics,
    including the number of RNA counts, the number of detected features, and the percentage of mitochondrial reads.

    Parameters:
    - adata (AnnData): Annotated Data object containing RNA-seq data.
    - outdir (str): Output directory for saving QC plots.
    - prefix (str): Prefix for the saved plot files.
    - dpi (int, optional): Dots per inch for saving the plot. Defaults to 300.

    Returns:
    None
    """
    # Create QC directory if it doesn't exist
    save_path = f"{outdir}/1.qc/"
    check_mkdir(save_path)

    # Define the QC plots to generate
    qc_plots = [
        ("nCount_RNA", "nFeature_RNA", f"{prefix}_nCount_RNA_nFeature_RNA"),
        ("nCount_RNA", "percent_mt", f"{prefix}_nCount_RNA_percent_mt"),
        (["nCount_RNA", "nFeature_RNA", "percent_mt"], "y", "QCviolin"),
    ]

    for x, y, filename in qc_plots:
        # Scatter plot or violin plot depending on the structure of the input
        if isinstance(x, str) and isinstance(y, str):
            # Scatter plot
            sc.pl.scatter(adata, x, y)
        elif isinstance(x, list) and isinstance(y, str):
            # Violin plot
            sc.pl.violin(adata, x, jitter=0.4, multi_panel=True)

        # Save the plot
        plt.savefig(
            f"{save_path}/{filename}.pdf",
            bbox_inches="tight",
        )
        plt.savefig(
            f"{save_path}/{filename}.png",
            bbox_inches="tight",
            dpi=dpi,
        )
        plt.close()
