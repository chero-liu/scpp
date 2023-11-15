import scanpy as sc
import numpy as np
import anndata
from typing import Union, Optional, Literal
import scanpy.external as sce
from anndata import AnnData
from scpp.tools.preprocessing.basicpreprocessing import scale
from scpp.tools.embedding.pca import run_pca

_Method = Literal["harmony", "combat", "scanorama", "No"]


def batch_correction(
    adata: AnnData,
    batch_key: str,
    methods: Optional[_Method] = "harmony",
    n_pcs: int = 50,
    **kwargs,
) -> AnnData:
    """
    Batch correction for single-cell data

    Arguments:
        adata: AnnData object
        batch_key: batch key
        methods: harmony,combat,scanorama
        n_pcs: number of PCs
        kwargs: other parameters for harmony`harmonypy.run_harmony()`,combat`sc.pp.combat()`,scanorama`scanorama.integrate_scanpy()`

    Returns:
        adata: AnnData object

    """

    print(f"...Begin using {methods} to correct batch effect")

    if methods == "harmony":
        try:
            import harmonypy

            # print('mofax have been install version:',mfx.__version__)
        except ImportError:
            raise ImportError("Please install the harmonypy: `pip install harmonypy`.")

        adata3 = adata.copy()
        scale(adata3)
        run_pca(adata3, n_comps=n_pcs)
        sce.pp.harmony_integrate(
            adata3, batch_key, basis="X_pca", adjusted_basis="X_pca_harmony", *kwargs
        )

        return adata3, "X_pca_harmony"
    elif methods == "combat":
        adata2 = adata.copy()
        sc.pp.combat(adata2, key=batch_key, *kwargs)
        scale(adata2)
        run_pca(adata2, n_comps=n_pcs)

        adata2.obsm["X_combat"] = adata2.X.copy()

        return adata2, "X_combat"
    elif methods == "scanorama":
        try:
            import scanorama

            # print('mofax have been install version:',mfx.__version__)
        except ImportError:
            raise ImportError("Please install the scanorama: `pip install scanorama`.")

        batches = adata.obs[batch_key].cat.categories.tolist()
        alldata = {}
        for batch in batches:
            alldata[batch] = adata[adata.obs[batch_key] == batch,]
        alldata2 = dict()
        for ds in alldata.keys():
            print(ds)
            alldata2[ds] = alldata[ds]

        # convert to list of AnnData objects
        adatas = list(alldata2.values())

        # run scanorama.integrate
        scanorama.integrate_scanpy(adatas, dimred=n_pcs, *kwargs)
        scanorama_int = [ad.obsm["X_scanorama"] for ad in adatas]

        # make into one matrix.
        all_s = np.concatenate(scanorama_int)
        print(all_s.shape)

        # add to the AnnData object, create a new object first
        adata.obsm["X_scanorama"] = all_s
        return adata, "X_scanorama"
    elif methods == "No":
        scale(adata)
        run_pca(adata, n_comps=n_pcs)
        return adata, "X_pca"
    else:
        print("Not supported")
