# pylint: disable=C0103, W0511, C0114
from __future__ import annotations

import logging

from anndata import AnnData
from ._utils import _harmony_integrate_python, _harmony_integrate_R


logger = logging.getLogger("symphonypy")


def harmony_integrate(
    adata: AnnData,
    key: list[str] | str,
    flavor: str = "python",
    ref_basis_source: str = "X_pca",
    ref_basis_adjusted: str = "X_pca_harmony",
    ref_basis_loadings: str = "PCs",
    verbose: bool = False,
    random_seed: int = 1,
    **harmony_kwargs,
):
    """
    Run Harmony batch correction on adata,
    save corrected output to adata.obsm,
    save all the necessary to Symphony mapping
        algorithm parameters to adata.uns

    Args:
        adata (AnnData): adata object with batch
        key (list[str] | str): which columns from adata.obs
            to use as batch keys (`vars_use` parameter of Harmony)
        flavor (str, optional): if to run harmonypy or Harmony via rpy2. Defaults to "python".
        ref_basis_source (str, optional): adata.obsm[ref_basis_source] will be used
            as input embedding to Harmony. Defaults to "X_pca".
        ref_basis_adjusted (str, optional): at adata.obsm[ref_basis_adjusted]
            corrected embedding will be saved. Defaults to "X_pca_harmony".
        ref_basis_loadings (str, optional): gene loadings of ref_basis_source. Defaults to "PCs".
        verbose (bool, optional): verbosity level of harmony. Defaults to False.
        random_seed (int, optional): random_seed for harmony. Defaults to 1.

    """

    if flavor == "python":
        if verbose:
            print("Harmony integration with harmonypy is preforming.")
        _harmony_integrate_python(
            adata=adata,
            key=key,
            ref_basis_source=ref_basis_source,
            ref_basis_adjusted=ref_basis_adjusted,
            ref_basis_loadings=ref_basis_loadings,
            verbose=verbose,
            random_state=random_seed,
            **harmony_kwargs,
        )
    elif flavor == "R":
        if verbose:
            print("Harmony integration with R Harmony is preforming.")

        _harmony_integrate_R(
            adata=adata,
            key=key,
            ref_basis_source=ref_basis_source,
            ref_basis_adjusted=ref_basis_adjusted,
            ref_basis_loadings=ref_basis_loadings,
            verbose=verbose,
            random_seed=random_seed,
            **harmony_kwargs,
        )
    else:
        raise Exception("`flavor` argument should be `python` or `R`.")
