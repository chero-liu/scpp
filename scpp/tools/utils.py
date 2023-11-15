import importlib
import os
import abc
import sys
import logging
import pandas as pd
import matplotlib.pyplot as plt
from functools import wraps
from anndata import AnnData
from typing import Literal, Dict, Any
import json

logging.basicConfig(level=logging.INFO)
from scpp.__init__ import ROOT_PATH


def find_assay_init(assay):
    init_module = importlib.import_module(f"scpp.{assay}.__init__")
    return init_module


def write_embedding(adata, key, embed_fn, n_comp=None, sep="\t"):
    """Export cell embeddings as a txt table"""

    if key not in adata.obsm.keys():
        raise KeyError(f"{key} is not a valid `.obsm` key")
    mat = adata.obsm[key].copy()
    if n_comp is not None and mat.shape[1] >= n_comp:
        mat = mat[:, 0:n_comp]
    df = pd.DataFrame(mat, index=adata.obs_names)
    df.columns = [
        f"{key.split('_')[1].upper()}_" + str(i + 1) for i in range(mat.shape[1])
    ]
    df.to_csv(embed_fn, sep=sep, header=True, index=True)


def check_mkdir(path):
    """
    Ensure that the specified directory exists. If it doesn't exist, create it.

    Parameters:
    - path (str): The path of the directory to be checked/created.

    Returns:
    - bool: True if the directory exists or was successfully created, False otherwise.
    """
    try:
        # Check if the directory exists
        if os.path.exists(path):
            print(f"The directory '{path}' already exists.")
            return True

        # Create the directory if it doesn't exist
        os.makedirs(path)
        print(f"Directory '{path}' created successfully.")
        return True

    except Exception as e:
        print(f"Error: {e}")
        return False


def find_step_module(assay, step):
    file_path_dict = {
        "assay": f"{ROOT_PATH}/{assay}/{step}.py",
        "tools": f"{ROOT_PATH}/tools/{step}.py",
    }

    init_module = find_assay_init(assay)

    if os.path.exists(file_path_dict["assay"]):
        step_module = importlib.import_module(f"scpp.{assay}.{step}")
    elif hasattr(init_module, "IMPORT_DICT") and step in init_module.IMPORT_DICT:
        module_path = init_module.IMPORT_DICT[step]
        step_module = importlib.import_module(f"{module_path}.{step}")
    elif os.path.exists(file_path_dict["tools"]):
        step_module = importlib.import_module(f"scpp.tools.{step}")
    else:
        raise ModuleNotFoundError(f"No module found for {assay}.{step}")

    return step_module


def save_parameters_to_anndata(adata: AnnData, parameters: Dict[str, Any]) -> AnnData:
    """
    Save parameters to the 'uns' attribute of an AnnData object.

    Parameters:
    - adata (AnnData): The AnnData object where parameters will be saved.
    - parameters (Dict[str, Any]): A dictionary containing parameters to be saved.

    Returns:
    AnnData: Updated AnnData object with saved parameters.
    """
    if "preprocess_para" in adata.uns:
        para = json.loads(adata.uns["preprocess_para"])
    else:
        para = {}
    for key, value in parameters.items():
        para[key] = value
    adata.uns["preprocess_para"] = json.dumps(para)

    return adata


def retrieve_layers(
    adata: AnnData,
    data: AnnData,
    layer_keys: Literal["raw", "normalised"] = ["raw", "normalised"],
    use_layer: str = "normalised",
) -> AnnData:
    """
    Retrieve specified layers from 'data' and update 'adata' object.

    Parameters:
    - adata (AnnData): The AnnData object to be updated with retrieved layers.
    - data (AnnData): The AnnData object containing layers to be retrieved.
    - layer_keys (Literal["raw", "normalised"]): List of layer keys to be retrieved. Default is ["raw", "normalised"].
    - use_layer (str): The layer to be used for 'adata.X'. Default is "normalised".

    Returns:
    AnnData: Updated 'adata' object.
    """
    # Create a new AnnData object for raw data
    raw_adata = adata.raw.to_adata()

    # Set raw data to itself
    raw_adata.raw = raw_adata

    # Update layers with specified data
    for layer in layer_keys:
        raw_adata.layers[layer] = data.layers[layer]

    # Update the original 'adata' object with the preprocessed data
    raw_adata.X = data.layers[use_layer]

    return raw_adata


def check_adata(adata: AnnData) -> AnnData:
    """
    Update `adata.X` based on the existence of a valid matrix named `raw` in `adata.layers`.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.

    Returns
    -------
    AnnData
    """
    if "raw" in adata.layers:
        adata.X = adata.layers["raw"].copy()
        return adata
    else:
        print("adata has no valid matrix raw")


def save_or_show(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        save_path = kwargs.pop("save_path", None)
        dpi = kwargs.pop("dpi", 300)
        func(*args, **kwargs)

        if save_path:
            directory = os.path.dirname(save_path)
            if directory and not os.path.exists(directory):
                os.makedirs(directory)

            plt.savefig(f"{save_path}.png", dpi=dpi)
            plt.savefig(f"{save_path}.pdf")
            plt.close()
        else:
            plt.show()

    return wrapper


def s_common(parser):
    """subparser common arguments"""
    parser.add_argument("--outdir", help="Output diretory.", required=True)
    parser.add_argument("--prefix", help="Prefix of all output files.", required=True)
    parser.add_argument("--species", help="Species", required=True)
    # parser.add_argument("--thread", help="", default=4)

    return parser


class Step:
    """
    Step class
    """

    def __init__(self, args):
        self.args = args
        self.outdir = args.outdir
        self.prefix = args.prefix
        self.species = args.species
        self.assay = args.subparser_assay
        # self.thread = int(args.thread)

        check_mkdir(self.outdir)

    @abc.abstractmethod
    def run(self):
        sys.exit("Please implement run() method.")

    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        print("done")
