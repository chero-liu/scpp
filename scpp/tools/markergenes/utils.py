import os
import scanpy as sc
from scpp.tools.utils import check_mkdir
from anndata import AnnData


def write_diff_genes(
    adata: AnnData, outdir: str, prefix: str, logfc: float = 0.25, minpct: float = 0.1
):
    """
    Write differential genes for specified groups to separate files and create a summary list.

    Parameters:
    - adata (AnnData): Annotated data matrix.
    - outdir (str): Output directory path.
    - prefix (str): Prefix for the output files.
    - groups (list): List of groups for which differential genes will be analyzed.
    - logfc (float, optional): Log-fold change threshold. Default is 0.25.
    - minpct (float, optional): Minimum percentage threshold. Default is 0.1.

    Returns:
    None
    """
    check_mkdir(f"{outdir}/diffgene/")
    diff_genes_summary = []
    groups = adata.uns["rank_genes_groups"]["names"].dtype.names
    for group in groups:
        de_i = sc.get.rank_genes_groups_df(adata, group=group)
        de_i = de_i[de_i.logfoldchanges.abs() >= logfc]
        de_i = de_i[
            (de_i["pct_nz_group"] >= minpct) | (de_i["pct_nz_reference"] >= minpct)
        ]

        output_path = f"{outdir}/diffgene/{prefix}_cluster{group}_diffgenes.xls"
        de_i.to_csv(output_path, index=False, sep="\t")

        diff_genes_summary.append((group, os.path.abspath(output_path)))

    # Write summary list
    summary_path = f"{outdir}/diffgene/{prefix}_diffgenes.list"
    with open(summary_path, "w") as output:
        for entry in diff_genes_summary:
            output.write(entry[0] + "\t" + entry[1] + "\n")
