# pylint: disable=C0103, C0116, C0114, C0115, W0511
from __future__ import annotations

import logging
import warnings

import numpy as np
import pandas as pd

from anndata import AnnData
from harmonypy import run_harmony
from scipy.sparse import issparse
from sklearn.cluster import KMeans

from scanpy.tools._ingest import Ingest, _DimDict

logger = logging.getLogger("symphonypy")


def _harmony_integrate_python(
    adata: AnnData,
    key: list[str] | str,
    ref_basis_source: str = "X_pca",
    ref_basis_adjusted: str = "X_pca_harmony",
    ref_basis_loadings: str = "PCs",
    verbose: bool = False,
    **harmony_kwargs,
) -> None:
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

    adata.uns["harmony"] = {
        # [K] the number of cells softly belonging to each cluster
        "Nr": ref_ho.R.sum(axis=1),
        # ref cluster centroids
        # [K, d]
        "C": C,
        # number of clusters
        "K": ref_ho.K,
        # sigma [K] (cluster cross enthropy regularization coef)
        "sigma": ref_ho.sigma,
        "ref_basis_loadings": ref_basis_loadings,
        "ref_basis_adjusted": ref_basis_adjusted,
        "vars_use": key,
        "harmony_kwargs": harmony_kwargs,
        "converged": converged,
        # [K, Nref]
        "R": ref_ho.R,
    }

    if not converged:
        logger.warning(
            "Harmony didn't converge. "
            "Consider increasing max_iter_harmony parameter value"
        )


def _harmony_integrate_R(
    adata: AnnData,
    key: list[str] | str,
    ref_basis_source: str = "X_pca",
    ref_basis_adjusted: str = "X_pca_harmony",
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
            adata.obsm[ref_basis_source],
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
        sigma = dollar(ho, "sigma")

        # [K, d] = [K, Nref] x [d, N_ref].T
        C = R @ Z_corr

        adata.uns["harmony"] = {
            "Nr": R.sum(axis=1),
            "C": C,
            "K": K[0],
            "sigma": sigma.squeeze(1),
            "ref_basis_loadings": ref_basis_loadings,
            "ref_basis_adjusted": ref_basis_adjusted,
            "converged": converged,
            "R": R,
        }
        adata.obsm[ref_basis_adjusted] = Z_corr

    if not converged:
        logger.warning(
            "Harmony didn't converge. Consider increasing max.iter.harmony parameter value"
        )


def _assign_clusters(
    X: np.array, sigma: float | np.array, Y: np.array, K: int
) -> np.array:

    if isinstance(sigma, (int, float)):
        sigma = np.array([sigma], dtype=np.float32)
    else:
        sigma = np.array(sigma, dtype=np.float32)
        assert (
            len(sigma) == K
        ), "sigma paramater must be either a single float or an array of length equal to number of clusters"

    # it's made so in harmonypy,
    # maybe to prevent overflow during L2 normalization?
    X_cos = X / X.max(axis=1, keepdims=True)
    # L2 normalization for cosine distance
    X_cos /= np.linalg.norm(X_cos, ord=2, axis=1, keepdims=True)

    # [K, N] = [K, d] x [Nq, d].T
    R = -2 * (1 - Y @ X_cos.T) / sigma[..., np.newaxis]
    R -= np.max(R, axis=0)  # maybe for numerical stability
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


def _adjust_for_missing_genes(
    adata: AnnData, use_genes_list: np.array, use_genes_list_present: np.array
) -> np.matrix:
    """
    Sets zero expression to missing genes, returns non-sparse matrix.
    Usefull when you will still do non-sparse-efficient normalization.
    Args:
        adata (AnnData): adata
        use_genes_list (np.array): which genes expressions to be left
        use_genes_list_present (np.array): which genes from use_genes_list are indeed present in adata

    Returns:
        np.matrix: non-sparse array of expressions
        of all the genes from use_genes_list
        with expressions of missing genes set to zero
    """
    logger.warning(
        "%i out of %i "
        "genes from the reference are missing in the query dataset or have zero std in the reference, "
        "their expressions in the query will be set to zero",
        (~use_genes_list_present).sum(),
        use_genes_list.shape[0],
    )
    t = np.zeros((adata.shape[0], use_genes_list.shape[0]))

    X = adata[:, use_genes_list[use_genes_list_present]].X
    t[:, use_genes_list_present] = X.A if issparse(X) else X.copy()

    return t


def _map_query_to_ref(
    adata_ref: AnnData,
    adata_query: AnnData,
    transferred_primary_basis: str = "X_pca_reference",
    ref_basis_loadings: str = "PCs",
    max_value: float | None = 10.0,
    use_genes_column: str | None = "highly_variable",
):
    assert (
        "mean" in adata_ref.var
    ), "Gene expression means are expected to be saved in adata_ref.var"
    assert (
        "std" in adata_ref.var
    ), "Gene expression stds are expected to be saved in adata_ref.var"

    use_genes_list = (
        adata_ref.var_names
        if use_genes_column is None
        else adata_ref.var_names[adata_ref.var[use_genes_column]]
    )

    # [N_genes]
    stds = (
        np.array(adata_ref.var["std"])
        if use_genes_column is None
        else np.array(adata_ref.var["std"][adata_ref.var[use_genes_column]])
    )

    use_genes_list_present = use_genes_list.isin(adata_query.var_names) & (stds != 0)

    # Adjusting for missing genes.
    if not all(use_genes_list_present):
        t = _adjust_for_missing_genes(
            adata_query, use_genes_list, use_genes_list_present
        )
    else:
        X = adata_query[:, use_genes_list].X
        t = X.A if issparse(X) else X.copy()

    means = (
        np.array(adata_ref.var["mean"][use_genes_list_present])
        if use_genes_column is None
        else np.array(
            adata_ref.var["mean"][adata_ref.var[use_genes_column]][
                use_genes_list_present
            ]
        )
    )
    stds = stds[use_genes_list_present]

    t[:, use_genes_list_present] -= means[np.newaxis]
    t[:, use_genes_list_present] /= stds[np.newaxis]

    if max_value is not None:
        t = np.clip(t, -max_value, max_value)

    # map query to reference's PCA coords
    # [cells, n_comps] = [cells, genes] * [genes, n_comps]
    adata_query.obsm[transferred_primary_basis] = np.array(
        t @ adata_ref.varm[ref_basis_loadings][adata_ref.var_names.isin(use_genes_list)]
    )


def _run_soft_kmeans(
    adata_ref: AnnData,
    ref_basis: str = "X_pca",
    K: int | None = None,
    ref_basis_loadings: str = "PCs",
    sigma: float = 0.1,
):
    N = adata_ref.shape[0]

    if K is None:
        K = np.min([np.round(N / 30.0), 100]).astype(int)

    model = KMeans(n_clusters=K, init="k-means++", n_init=10, max_iter=25)
    model.fit(adata_ref.obsm[ref_basis])
    C = model.cluster_centers_

    Z_cos = adata_ref.obsm[ref_basis]
    Z_cos /= Z_cos.max(axis=1, keepdims=True)
    Z_cos /= np.linalg.norm(Z_cos, ord=2, axis=1, keepdims=True)

    # (1) Normalize
    Y = C / np.linalg.norm(C, ord=2, axis=1, keepdims=True)
    # (2) Assign cluster probabilities
    R = _assign_clusters(Z_cos, sigma, Y, K)

    adata_ref.uns["harmony"] = {
        # [K] the number of cells softly belonging to each cluster
        "Nr": R.sum(axis=1),
        # ref cluster centroids
        # [K, d]
        "C": C,
        # number of clusters
        "K": K,
        # sigma [K] (cluster cross enthropy regularization coef)
        "sigma": None,
        "ref_basis_loadings": ref_basis_loadings,
        "ref_basis_adjusted": ref_basis,
        "vars_use": None,
        "harmony_kwargs": {},
        "converged": True,
        "R": R,
    }


class Ingest_sp(Ingest):
    def __init__(self, *args, **kwargs):

        self._force_use_rep = kwargs.pop("use_representation", None)

        super().__init__(*args, **kwargs)

    def fit(self, adata_new: AnnData):
        self._obs = pd.DataFrame(index=adata_new.obs.index)
        self._obsm = _DimDict(adata_new.n_obs, axis=0)

        self._adata_new = adata_new
        self._obsm["rep"] = self._same_rep()

    def _same_rep(self):
        adata = self._adata_new

        if self._force_use_rep is not None:
            if not self._force_use_rep == self._use_rep:
                warnings.warn(
                    f"'{self._use_rep}' adata_reference's representation was used for reference, "
                    f"while '{self._force_use_rep}' was used for query. Be sure that they correspond."
                )
            return adata.obsm[self._force_use_rep]

        if self._n_pcs is not None:
            return self._pca(self._n_pcs)

        if self._use_rep == "X":
            assert all(adata.var_names == self._adata_ref.var_names), (
                "Can't use 'X' as a representation because var_names "
                "do not match between reference and query"
            )
            return adata.X

        if self._use_rep in adata.obsm.keys():
            return adata.obsm[self._use_rep]

        return adata.X

    def _pca(self, n_pcs=None):

        if self._pca_use_hvg:
            use_genes_list = self._adata_ref.var_names[
                self._adata_ref.var["highly_variable"]
            ]
        else:
            use_genes_list = self._adata_ref.var_names

        use_genes_list_present = use_genes_list.isin(self._adata_new.var_names)

        if not all(use_genes_list_present):
            X = _adjust_for_missing_genes(
                self._adata_new, use_genes_list, use_genes_list_present
            )
        else:
            X = self._adata_new[:, use_genes_list].X
            X = X.toarray() if issparse(X) else X.copy()

        if self._pca_centered:
            X -= X.mean(axis=0)

        X_pca = np.dot(X, self._pca_basis[:, :n_pcs])

        return X_pca


def _cluster_maha_dist(
    adata_query_cluster_obs: pd.DataFrame,
    adata_query: AnnData,
    transferred_primary_basis: str,
    reference_cluster_centroids: np.array,
    reference_cluster_centroids_norm: np.array,
    u: float,
    lamb: float,
):
    adata_cluster = adata_query[adata_query_cluster_obs.index, :]
    d = reference_cluster_centroids.shape[1]
    # [Nq, d]
    query_coords = adata_cluster.obsm[transferred_primary_basis]

    # Calculate query cluster centroid and covariances in PC space
    query_cluster_centroid = query_coords.mean(axis=0)
    query_cluster_centered = query_coords - query_cluster_centroid
    query_cluster_cov = (
        query_cluster_centered.T
        @ query_cluster_centered
        / (query_cluster_centered.shape[0] - 1)
    )
    query_cluster_centroid_norm = query_cluster_centroid

    # Find nearest reference cluster centroid
    ref_centroid_closest_idx = np.argmax(
        reference_cluster_centroids_norm @ query_cluster_centroid_norm
    )
    ref_centroid_closest = reference_cluster_centroids[ref_centroid_closest_idx]

    # Calculate Mahalanobis distance from query cluster to nearest reference centroid
    cluster_size = adata_query_cluster_obs.shape[0]
    if cluster_size < u * d:
        logger.info(
            "cluster contains too few cells to estimate confidence: ', query_cluster_labels_unique[c])"
        )
        dist = None
    else:
        query_cluster_cov += np.diag([lamb] * d)
        inv_cluster_cov = np.linalg.inv(query_cluster_cov)
        dif = ref_centroid_closest - query_cluster_centroid
        dist = float((dif @ inv_cluster_cov @ dif) ** 0.5)

    return dist


def _cluster_covs(X_ref, R, K):
    R_ = R[:, np.newaxis]

    # [d, N_ref]
    X_ref_T = X_ref.T
    # [K, d, N_ref] = [K, 1, Nref] X [d, N_ref]
    X_weighted = np.multiply(R_, X_ref_T)
    # [K, 1]
    v1 = R.sum(axis=1, keepdims=True)
    v2 = (R**2).sum(axis=1, keepdims=True)
    v1_2 = v1**2
    # [K, d]
    X_weighted_mean = X_weighted.sum(axis=2) / v1
    X_weighted_mean = X_weighted_mean[..., np.newaxis]
    # [1, d, N_ref]
    X_ref_T_ = X_ref_T[np.newaxis]
    # [K, d, N_ref]
    X_ref_K = np.repeat(X_ref_T_, K, axis=0)
    # [K, d, N_ref] = [K, d, N_ref] - [K, d, 1]
    X_centered = X_ref_K - X_weighted_mean
    # [K, d, N_ref] = [K, 1, Nref] X [K, d, N_ref]
    X_centered_weighted = np.multiply(R_, X_centered)

    # [K, d, d] = [K, d, N_ref] * [K, d, N_ref]
    cluster_covs = np.einsum("npq,nrq->npr", X_centered_weighted, X_centered)
    cluster_covs = cluster_covs * v1[..., np.newaxis] / (v1_2 - v2)[..., np.newaxis]
    inv_cluster_covs = np.linalg.inv(cluster_covs)

    return inv_cluster_covs, X_weighted_mean
