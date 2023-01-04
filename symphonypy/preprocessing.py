# pylint: disable=C0103, W0511, C0114
from __future__ import annotations

import logging

import numpy as np
from anndata import AnnData
from harmonypy import run_harmony


logger = logging.getLogger("symphonypy")


def harmony_integrate(
    adata: AnnData,
    key: list[str] | str,
    *harmony_args,
    ref_basis_source: str = "X_pca",
    ref_basis_adjusted: str = "X_pca_harmony",
    ref_basis_loadings: str = "PCs",
    verbose: bool = False,
    **harmony_kwargs,
):
    """
    Run Harmony batch correction on adata,
    save corrected output to adata.obsm,
    save all necessary to Symphony mapping
        algorithm parameters to adata.uns

    Args:
        adata (ad.AnnData): reference adata object
        basis_source (str): adata.obsm[basis_source] will be used
            as input embedding to Harmony
        basis_adjusted (str): at adata.obsm[basis_adjusted]
            corrected embedding will be saved
        key (list[str] | str): which columns from adata.obs
            to use as batch keys (`vars_use` parameter of Harmony)
    """
    ref_ho = run_harmony(
        adata.obsm[ref_basis_source],
        meta_data=adata.obs,
        vars_use=key,
        verbose=verbose,
        *harmony_args,
        **harmony_kwargs,
    )

    adata.obsm[ref_basis_adjusted] = ref_ho.Z_corr.T

    converged = ref_ho.check_convergence(1)

    adata.uns["harmony"] = {
        # [K] the number of cells softly belonging to each cluster
        "Nr": ref_ho.R.sum(axis=1),
        # [K, d] = [K, Nref] x [d, N_ref].T
        "C": ref_ho.R @ ref_ho.Z_corr.T,
        # ref cluster centroids L2 normalized
        # [K, d] = [d, K].T
        "Y": ref_ho.Y.T,
        # number of clusters
        "K": ref_ho.K,
        # sigma [K] (cluster cross enthropy regularization coef)
        "sigma": ref_ho.sigma,
        "ref_basis_loadings": ref_basis_loadings,
        "ref_basis_adjusted": ref_basis_adjusted,
        "vars_use": key,
        "harmony_args": harmony_args,
        "harmony_kwargs": harmony_kwargs,
        "converged": converged,
    }

    if not converged:
        logger.warning(
            "Harmony didn't converge. Consider increasing max_iter_harmony parameter value"
        )
