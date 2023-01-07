# pylint: disable=C0103, C0116, C0114, C0115, W0511
from __future__ import annotations

import logging

import numpy as np

from anndata import AnnData
from harmonypy import run_harmony

logger = logging.getLogger("symphonypy")


# TODO:
# def _compute_confidence()


def _harmony_integrate_python(
    adata: AnnData,
    key: list[str] | str,
    ref_basis_source: str = "X_pca",
    ref_basis_adjusted: str = "X_pca_harmony",
    ref_basis_loadings: str = "PCs",
    verbose: bool = False,
    **harmony_kwargs,
):
    ref_ho = run_harmony(
        adata.obsm[ref_basis_source],
        meta_data=adata.obs,
        vars_use=key,
        verbose=verbose,
        **harmony_kwargs,
    )

    adata.obsm[ref_basis_adjusted] = ref_ho.Z_corr.T

    converged = ref_ho.check_convergence(1)

    # [K, d] = [K, Nref] x [d, N_ref].T
    C = ref_ho.R @ ref_ho.Z_corr.T
    Y = C / np.linalg.norm(C, ord=2, axis=1, keepdims=True)

    adata.uns["harmony"] = {
        # [K] the number of cells softly belonging to each cluster
        "Nr": ref_ho.R.sum(axis=1),
        # [K, d]
        "C": C,
        # ref cluster centroids L2 normalized
        # [K, d]
        "Y": Y,
        # number of clusters
        "K": ref_ho.K,
        # sigma [K] (cluster cross enthropy regularization coef)
        "sigma": ref_ho.sigma,
        "ref_basis_loadings": ref_basis_loadings,
        "ref_basis_adjusted": ref_basis_adjusted,
        "vars_use": key,
        "harmony_kwargs": harmony_kwargs,
        "converged": converged,
    }

    if not converged:
        logger.warning(
            "Harmony didn't converge. Consider increasing max_iter_harmony parameter value"
        )


def _harmony_integrate_R(
    adata: AnnData,
    key: list[str] | str,
    basis: str = "X_pca",
    adjusted_basis: str = "X_pca_harmony",
    ref_basis_loadings: str = "PCs",
    random_seed: int = 1,
    verbose: bool = False,
    **harmony_kwargs,
) -> None:
    """
    Function description. Exhaustively.
    """

    import shutil

    if not shutil.which("R"):
        raise Exception("R installation is necessary.")
    try:
        from rpy2.robjects.packages import importr
    except ImportError as e:
        raise ImportError("\nPlease install rpy2:\n\n\tpip install rpy2") from e
    try:
        harmony = importr("harmony")
    except Exception as e:
        raise Exception(
            'R package "Harmony" is necessary.\n'
            "Please install it from https://github.com/immunogenomics/harmony and try again"
        ) from e

    from rpy2.robjects import numpy2ri, pandas2ri, default_converter
    from rpy2.robjects.conversion import localconverter
    import rpy2.rinterface_lib.callbacks

    rpy2.rinterface_lib.callbacks.consolewrite_warnerror = lambda x: print(x, end="")

    dollar = importr("base").__dict__["$"]

    with localconverter(
        default_converter + numpy2ri.converter + pandas2ri.converter
    ) as cv:

        importr("base").set_seed(random_seed)

        ho = harmony.HarmonyMatrix(
            adata.obsm[basis],
            adata.obs,
            key,
            do_pca=False,
            return_object=True,
            verbose=verbose,
            **harmony_kwargs,
        )

        converged = dollar(ho, "check_convergence")(1)[0]

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
            "ref_basis_loadings": ref_basis_loadings,
            "ref_basis_adjusted": adjusted_basis,
            "converged": converged,
        }
        adata.obsm[adjusted_basis] = Z_corr

    if not converged:
        logger.warning(
            "Harmony didn't converge. Consider increasing max.iter.harmony parameter value"
        )


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
    R = -2 * (1 - Y @ X_cos.T) / sigma[..., np.newaxis]
    R = np.exp(R)
    R /= R.sum(axis=0, keepdims=True)

    return R


def _correct_query(
    X: np.array,
    phi_: np.array,
    R: np.array,
    Nr: np.array,
    C: np.array,
    lamb: np.array,
):
    # [d, N] = [N, d].T
    X_corr = X.copy().T
    lamb[0, 0] = 0
    K = R.shape[0]

    for i in range(K):
        # [B + 1, N] = [B + 1, N] * [N]
        Phi_Rk = np.multiply(phi_, R[i, :])

        # [B + 1, B + 1] = [B + 1, N] x [N, B + 1]
        x = Phi_Rk @ phi_.T
        # += [1]
        x[0, 0] += Nr[i]

        # [B + 1, d] = [B + 1, N] x [N, d]
        y = Phi_Rk @ X
        y[0, :] += C[i]

        # [B + 1, d] = [B + 1, B + 1] x [B + 1, d]
        W = np.linalg.inv(x + lamb) @ y
        W[0, :] = 0  # do not remove the intercept

        # [d, N] -= [B + 1, d].T x [B + 1, N]
        X_corr -= W.T @ Phi_Rk

    return X_corr.T


def _map_query_to_ref(
    adata_ref: AnnData,
    adata_query: AnnData,
    query_basis_ref: str = "X_pca_reference",
    ref_basis_loadings: str = "PCs",
    max_value: float | None = 10.0,
    use_genes_column: str = "highly_variable",
):
    assert (
        "mean" in adata_ref.var
    ), "Gene expression means are expected to be saved in adata_ref.var"
    assert (
        "std" in adata_ref.var
    ), "Gene expression stds are expected to be saved in adata_ref.var"

    use_genes_list = np.array(adata_ref.var_names[adata_ref.var[use_genes_column]])

    # [N_genes]
    stds = np.array(adata_ref.var["std"][adata_ref.var[use_genes_column]])

    use_genes_list_present = np.isin(use_genes_list, adata_query.var_names) & (
        stds != 0
    )

    # Adjusting for missing genes.
    if not all(use_genes_list_present):
        logger.warning(
            "%i out of %i "
            "genes from reference are missing in query dataset or have zero std in reference,"
            "their expressions in query will be set to zero",
            (~use_genes_list_present).sum(),
            use_genes_list.shape[0],
        )
        # anyway after scaling data wouldn't be sparse efficient,
        # so will convert it to array
        t = np.zeros((adata_query.shape[0], use_genes_list.shape[0]))
        t[:, use_genes_list_present] = adata_query[
            :, use_genes_list[use_genes_list_present]
        ].X.A
    else:
        t = adata_query[:, use_genes_list].X.A

    means = np.array(
        adata_ref.var["mean"][adata_ref.var[use_genes_column]][use_genes_list_present]
    )
    stds = stds[use_genes_list_present]

    t[:, use_genes_list_present] = (
        t[:, use_genes_list_present] - means[np.newaxis]
    ) / stds[np.newaxis]

    if max_value is not None:
        t = np.clip(t, -max_value, max_value)

    # map query to reference's PCA coords
    # [cells, n_comps] = [cells, genes] * [genes, n_comps]
    adata_query.obsm[query_basis_ref] = np.array(
        t @ adata_ref.varm[ref_basis_loadings][adata_ref.var_names.isin(use_genes_list)]
    )
