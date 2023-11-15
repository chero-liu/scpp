import scanpy as sc
import pandas as pd
import numpy as np
from collections import Counter
from scpp.tools.read.utils import getSingleDataFormat, read, ensureSparseMatrix
from scpp.tools.utils import check_mkdir, save_parameters_to_anndata
from scpp.tools.preprocessing.basicpreprocessing import (
    _filter,
    findMtThreshold,
    calculate_mt_common,
)

MT_THRESHOLD = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50]


class Create:
    def __init__(
        self,
        info: None,
        outdir: str,
        prefix: str,
        species: str,
        mtfilter: str = "default",
        mingene: float = 200,
        minumi: float = 0,
        maxumi: float = 0.98,
        maxgene: float = 0.98,
        mincell: int = 5,
        dpi: int = 300,
    ):
        self.info = info
        self.outdir = outdir
        self.prefix = prefix
        self.species = species
        self.mtfilter = mtfilter
        self.mingene = mingene
        self.minumi = minumi
        self.maxgene = maxgene
        self.maxumi = maxumi
        self.mincell = mincell
        self.dpi = dpi

    def run(self):
        self.info = pd.read_csv(self.info, sep="\t")
        adatas, mt_list, genes_list = [], [], []
        raw_sample, raw_cells, raw_genes, raw_umi, raw_gene = [], [], [], [], []
        filter_sample, filter_cells, filter_genes, filter_umi, filter_gene = (
            [],
            [],
            [],
            [],
            [],
        )

        for index, line in self.info.iterrows():
            print(f"Start processing {line.spname} ...")
            data_type = getSingleDataFormat(line.path)
            adata = read(line.path, prefix=line.spname, dtype=data_type)
            ensureSparseMatrix(adata)

            df = sc.pp.calculate_qc_metrics(
                adata, percent_top=None, log1p=False, inplace=False
            )[0]
            raw_sample.append(line.spname)
            raw_cells.append(adata.shape[0])
            raw_genes.append(adata.shape[1])
            raw_umi.append(np.median(df["total_counts"]))
            raw_gene.append(np.median(df["n_genes_by_counts"]))

            _filter(
                adata,
                mingenes=self.mingene,
                minumis=self.minumi,
                maxgenes=self.maxgene,
                maxumis=self.maxumi,
                mincells=self.mincell,
            )

            add_names = line.index.drop(["path"])
            for tag in add_names:
                adata.obs[tag] = str(line[tag])
            mt_ = findMtThreshold(
                adata,
                mt_thresholds=MT_THRESHOLD,
                species=self.species,
                mtfilter=self.mtfilter,
            )
            mt_list.append(mt_)
            genes_list.append(adata.shape[1])
            adatas.append(adata)

        adata = sc.concat(adatas, join="outer")
        adata.var_names_make_unique()
        adata.obs_names_make_unique()
        adata.X = np.nan_to_num(adata.X)

        mt_ = findMtThreshold(
            adata,
            mt_thresholds=MT_THRESHOLD,
            species=self.species,
            mtfilter=self.mtfilter,
        )
        sc.pp.calculate_qc_metrics(
            adata,
            qc_vars=["mt"],
            percent_top=None,
            log1p=False,
            inplace=True,
        )

        mt_common = calculate_mt_common(self.mtfilter, mt_list)

        mt_dict = {}
        mt_raw = Counter(adata.obs["spname"])
        adata = adata[adata.obs.pct_counts_mt < int(mt_common), :]
        mt_filter = Counter(adata.obs["spname"])
        for key, value in mt_raw.items():
            mt_dict[key] = mt_filter[key] / mt_raw[key]
        mt_df = pd.DataFrame(
            {"SampleID": mt_dict.keys(), "mt_filtered_percent": mt_dict.values()}
        )

        # re-order obs slot
        reorder_slot = ["rawname", "spname", "gname"]
        for slot in reorder_slot:
            adata.obs[slot] = adata.obs[slot].astype("category")
            if slot == "gname":
                try:
                    order = pd.DataFrame(self.info.gname.tolist()).drop_duplicates()
                except AttributeError:
                    order = pd.DataFrame(self.info.spname.tolist()).drop_duplicates()
                order = order[0].tolist()
            else:
                order = self.info[slot].tolist()

            order = [str(x) for x in order]
            adata.obs[slot] = adata.obs[slot].cat.set_categories(order)

        # rename obs slot name
        adata.obs.rename(
            columns={
                "n_genes_by_counts": "nFeature_RNA",
                "total_counts": "nCount_RNA",
                "pct_counts_mt": "percent_mt",
            },
            inplace=True,
        )

        adata.obs.rename(columns={"spname": "sample"}, inplace=True)

        del adata.obs["n_genes"]
        del adata.obs["n_counts"]

        for sample in set(adata.obs["sample"]):
            tmp = adata[adata.obs["sample"] == sample]
            filter_sample.append(sample)
            filter_cells.append(tmp.shape[0])
            filter_genes.append(tmp.shape[1])
            filter_umi.append(np.median(tmp.obs["nCount_RNA"]))
            filter_gene.append(np.median(tmp.obs["nFeature_RNA"]))

        raw_df = pd.DataFrame(
            dict(
                SampleID=raw_sample,
                raw_cells=raw_cells,
                raw_genes=raw_genes,
                raw_median_umi=raw_umi,
                raw_median_gene=raw_gene,
            )
        )
        filter_df = pd.DataFrame(
            dict(
                SampleID=filter_sample,
                filtered_cells=filter_cells,
                filtered_genes=filter_genes,
                filtered_median_umi=filter_umi,
                filtered_median_gene=filter_gene,
            )
        )
        qc_df = pd.merge(raw_df, filter_df, on="SampleID")
        qc_df = pd.merge(qc_df, mt_df, on="SampleID")
        check_mkdir(self.outdir)
        qc_df.to_csv(
            self.outdir + "/" + "{0}_qc_config.xls".format(self.prefix),
            sep="\t",
            index=False,
        )

        argu = {
            "mingene": str(self.mingene),
            "minumi": str(self.minumi),
            "maxgene": str(self.maxgene),
            "maxumi": str(self.maxumi),
            "mincell": str(self.mincell),
            "mtfilter": str(mt_common),
        }
        adata = save_parameters_to_anndata(adata, argu)

        return adata
