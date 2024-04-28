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
    Run Harmony batch correction on adata, save corrected output to ``adata.obsm``,
    save all the necessary to Symphony mapping algorithm parameters to ``adata.uns``

    :param adata: adata object with batch
    :type adata: AnnData
    :param key: which columns from ``adata.obs`` to use as batch keys (``vars_use`` parameter of Harmony)
    :type key: list[str] | str
    :param flavor: if to run harmonypy or Harmony via ``rpy2``, defaults to "python"
    :type flavor: str, optional
    :param ref_basis_source: ``adata.obsm[ref_basis_source]`` will be used as input embedding to Harmony, defaults to "X_pca"
    :type ref_basis_source: str, optional
    :param ref_basis_adjusted: slot where to put corrected coordinates, defaults to "X_pca_harmony"
    :type ref_basis_adjusted: str, optional
    :param ref_basis_loadings: slot with feature loadings to the original embedding, defaults to "PCs"
    :type ref_basis_loadings: str, optional
    :param verbose: if to print logs of steps of integration, defaults to False
    :type verbose: bool, optional
    :param random_seed: random seed, defaults to 1
    :type random_seed: int, optional
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
