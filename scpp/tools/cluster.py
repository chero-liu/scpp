from anndata import AnnData
from scpp.tools.utils import check_adata, retrieve_layers, save_parameters_to_anndata
from scpp.tools.preprocessing.highly_variable_genes import variable_filter
from scpp.tools.embedding._batch import batch_correction
from scpp.tools.preprocessing.graph import neighbors, clustering
from scpp.tools.preprocessing.utils import transform_cluster_labels, umap_argue
from scpp.tools.embedding.embedding import run_umap, run_tsne
from scpp.tools.preprocessing.basicpreprocessing import normalize


class Cluster:
    def __init__(
        self,
        adata: AnnData,
        hvg_num: int = 2000,
        n_pcs: int = 10,
        batch_method: str = "No",
        batch_name: str = "sample",
        clust_method: str = "leiden",
        min_dist: str = "default",
        resolution: float = 1.2,
        target_sum: int = 10000,
    ):
        self.adata = adata
        self.hvg_num = hvg_num
        self.n_pcs = n_pcs
        self.batch_method = batch_method
        self.batch_name = batch_name
        self.clust_method = clust_method
        self.min_dist = min_dist
        self.resolution = resolution
        self.target_sum = target_sum

    def run(self):
        print(f"Start Clustering ...")

        self.adata = check_adata(self.adata)

        self.adata.raw = self.adata
        print("Raw adata is already saved in adata.raw")

        normalize(self.adata, target_sum=self.target_sum)

        self.adata.layers["normalised"] = self.adata.X.copy()
        print("Normalization step is finished ,copy to adata.layers['normalised']")

        data = self.adata.copy()

        self.adata = variable_filter(
            self.adata,
            "normalised",
            self.hvg_num,
            flavor="seurat_v3",
        )
        variable_filter_result = self.adata.var.copy()
        self.adata = self.adata[:, self.adata.var.highly_variable]

        self.adata, use_rep = batch_correction(
            self.adata,
            batch_key=self.batch_name,
            methods=self.batch_method,
        )

        neighbors(
            self.adata,
            n_neighbors=15,
            n_pcs=self.n_pcs,
            use_rep=use_rep,
        )

        clustering(
            self.adata,
            self.clust_method,
            self.resolution,
        )

        transform_cluster_labels(
            self.adata,
            self.clust_method,
        )

        if self.min_dist == "default":
            self.min_dist = umap_argue(self.adata)
        else:
            self.min_dist = float(self.min_dist)

        run_umap(
            self.adata,
            min_dist=self.min_dist,
        )

        run_tsne(
            self.adata,
            n_pcs=self.n_pcs,
        )

        self.adata = retrieve_layers(
            self.adata,
            data,
        )
        self.adata.var = variable_filter_result
        argu = {
            "dim": str(self.n_pcs),
            "resolution": str(self.resolution),
            "batch_method": self.batch_method,
            "clust_method": self.clust_method,
            "min_dist": str(self.min_dist),
        }
        self.adata = save_parameters_to_anndata(self.adata, argu)

        return self.adata
