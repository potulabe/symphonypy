# pylint: disable=C0103, C0116, C0114, C0115, W0511

import logging
from typing import List, Union
import numpy as np

from anndata import AnnData

logger = logging.getLogger("symphonypy")


# TODO:
# def _compute_confidence()


def _harmony_integrate_R(
    adata: AnnData,
    key: Union[List[str], str],
    basis: str = "X_pca",
    adjusted_basis: str = "X_pca_harmony",
    **kwargs,
) -> None:
    """
    Function description.
    """
    import shutil

    if not shutil.which("R"):
        raise Exception(
            "R installation is necessary."
        )
    try:
        from rpy2.robjects.packages import importr
    except ImportError:
        raise ImportError("\nPlease install rpy2:\n\n\tpip install rpy2")
    try:
        harmony = importr("harmony")
    except Exception as e:
        raise Exception(
            'R package "Harmony" is necessary.\n'
            'Please install it from https://github.com/immunogenomics/harmony and try again'
        )
        
    from rpy2.robjects import numpy2ri, pandas2ri, default_converter
    from rpy2.robjects.conversion import localconverter
    import rpy2.rinterface_lib.callbacks
    rpy2.rinterface_lib.callbacks.consolewrite_warnerror = lambda x: print(x, end="")

    dollar = importr("base").__dict__["$"]

    with localconverter(
        default_converter + numpy2ri.converter + pandas2ri.converter
    ) as cv:
        ho = harmony.HarmonyMatrix(
            adata.obsm[basis],
            adata.obs,
            key,
            do_pca=False,
            return_object=True,
            **kwargs,
        )
        R = dollar(ho, "R")
        Z_corr = dollar(ho, "Z_corr").T
        K = dollar(ho, "K")
        Y = dollar(ho, "Y")
        sigma = dollar(ho, "sigma")

        adata.uns["harmony"] = {
            "Nr": R.sum(axis=1),
            "C": R @ Z_corr,
            "Y": Y.T,
            "K": K[0],
            "sigma": sigma.squeeze(1),
            "ref_basis_loadings": "PCs",
            "ref_basis_adjusted": adjusted_basis,
        }
        adata.obsm[adjusted_basis] = Z_corr

def _assign_clusters(X: np.array, sigma: np.array, Y: np.array) -> np.array:
    """_summary_
    Args:
        X (np.array): _description_
        sigma (np.array): _description_
        Y (np.array): _description_
    Returns:
        np.array: _description_
    """
    # it's made so in harmonypy,
    # maybe to prevent overflow during L2 normalization?
    X_cos = X / X.max(axis=1, keepdims=True)
    # L2 normalization for cosine distance
    X_cos /= np.linalg.norm(X_cos, ord=2, axis=1, keepdims=True)

    # [K, N] = [K, d] x [Nq, d].T
    R = -2 * (1 - np.dot(Y, X_cos.T)) / sigma[..., np.newaxis]
    R = np.exp(R)
    R /= R.sum(axis=0, keepdims=True)

    return R


def _correct_query(
    X: np.array,
    phi_: np.array,
    R: np.array,
    K: int,
    Nr: np.array,
    C: np.array,
    lamb: np.array,
):
    # [d, N] = [N, d].T
    X_corr = X.copy().T

    for i in range(K):
        # [B + 1, N] = [B + 1, N] * [N]
        Phi_Rk = np.multiply(phi_, R[i, :])

        # [B + 1, B + 1] = [B + 1, N] x [N, B + 1]
        x = np.dot(Phi_Rk, phi_.T)
        # += [1]
        x[0, 0] += Nr[i]

        # [B + 1, d] = [B + 1, N] x [N, d]
        y = np.dot(Phi_Rk, X)
        y[0, :] += C[i]

        # [B + 1, d] = [B + 1, B + 1] x [B + 1, d]
        W = np.dot(np.linalg.inv(x + lamb), y)
        W[0, :] = 0  # do not remove the intercept

        # [d, N] -= [B + 1, d].T x [B + 1, N]
        X_corr -= np.dot(W.T, Phi_Rk)

    return X_corr.T


def _map_query_to_ref(
    adata_ref: AnnData,
    adata_query: AnnData,
    query_basis_ref: str = "X_pca_reference",
    ref_basis_loadings: str = "PCs",
    max_value: Union[float, None] = 10.0,
    use_genes_column: str = "highly_variable",
):
    use_genes_list = np.array(adata_ref.var_names[adata_ref.var[use_genes_column]])

    stds = np.array(adata_ref.var["std"][adata_ref.var[use_genes_column]])
    means = np.array(adata_ref.var["mean"][adata_ref.var[use_genes_column]])[stds != 0]

    use_genes_list_present = np.isin(use_genes_list, adata_query.var_names)

    # Adjusting for missing genes.
    if not all(use_genes_list_present):
        logger.warning(
            "%i out of %i "
            "genes from reference are missing in query dataset, "
            "their expressions will be set to zero ",
            (~use_genes_list_present).sum(),
            use_genes_list.shape[0],
        )
        # anyway after scaling data wouldn't be sparse efficient,
        # so will convert it to array
        t = np.zeros((adata_query.shape[0], use_genes_list.shape[0]))
        t[:, use_genes_list_present] = adata_query[
            :, use_genes_list[use_genes_list_present]
        ].X.A
        t = t[:, stds != 0]

    else:
        t = adata_query[:, use_genes_list].X[:, stds != 0].copy()

    stds = stds[stds != 0]

    t = (t - means[np.newaxis]) / stds[np.newaxis]
    t = np.clip(t, None, max_value)

    # set zero expression to missing genes after scaling
    t[:, ~use_genes_list_present] = 0

    # map query to reference's PCA coords
    # [cells, n_comps] = [cells, genes] * [genes, n_comps]
    adata_query.obsm[query_basis_ref] = np.array(
        t @ adata_ref.varm[ref_basis_loadings][adata_ref.var_names.isin(use_genes_list)]
    )
