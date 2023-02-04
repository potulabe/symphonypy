# pylint: disable=C0103, C0116, C0114, C0115, W0511
from __future__ import annotations

import logging
import warnings

from typing import Iterable, List, Optional, Union

import numpy as np
import pandas as pd

from anndata import AnnData
from packaging import version
from sklearn.neighbors import KNeighborsClassifier

from scanpy._compat import pkg_version

from ._utils import (
    _assign_clusters,
    _correct_query,
    _map_query_to_ref,
    Ingest_sp,
    _run_soft_kmeans,
    _cluster_maha_dist,
    _cluster_covs,
)


ANNDATA_MIN_VERSION = version.parse("0.7rc1")
logger = logging.getLogger("symphonypy")


def per_cell_confidence(
    adata_query: AnnData,
    adata_ref: AnnData,
    ref_basis_adjusted: str = "X_pca_harmony",
    query_basis_adjusted: str = "X_pca_harmony",
    transferred_primary_basis: str = "X_pca_reference",
    obs: str = "symphony_per_cell_dist",
):
    """
    Calculates the weighted Mahalanobis distance for query cells to reference clusters.
    Higher distance metric indicates less confidence.
    Saves the metric to `adata_query.obs[obs]`

    Args:
        adata_query (AnnData): query adata object mapped to `adata_ref` with Symphony
        adata_ref (AnnData): reference adata object (with Harmony object in adata_ref.uns)
        ref_basis_adjusted (str, optional): adata_ref.obsm[ref_basis_adjusted] should contain resulting (harmony integrated if batch was present) reference representation. Defaults to "X_pca_harmony".
        query_basis_adjusted (str, optional): adata_query.obsm[query_basis_adjusted] should contain symphony adjusted query representation. Defaults to "X_pca_harmony".
        transferred_primary_basis (str, optional): adata_query.obsm[transferred_primary_basis] should contain pre-Symphony reference PC query representation. Defaults to "X_pca_reference".
        obs (str, optional): at adata_query.obs[obs] confidence metric will be saved. Defaults to "symphony_confidence".
    """

    assert (
        "harmony" in adata_ref.uns
    ), "Harmony object not found in adata_ref.uns['harmony']. First, run standard symphony query mapping."

    harmony = adata_ref.uns["harmony"]

    # [K, Nref]
    R = harmony["R"]
    K, _ = R.shape

    if ("inv_cluster_covs" not in harmony) or ("cluster_centers" not in harmony):
        # count reference's covariance matrices

        # [N_ref, d]
        X_ref = adata_ref.obsm[ref_basis_adjusted]

        inv_cluster_covs, X_weighted_mean = _cluster_covs(X_ref, R, K)
        # save
        harmony["inv_cluster_covs"] = inv_cluster_covs
        harmony["cluster_centers"] = X_weighted_mean

    else:
        logger.info(
            "Using precomputed reference's clusters covariance matrices from the Harmony object"
        )
        inv_cluster_covs = harmony["inv_cluster_covs"]
        X_weighted_mean = harmony["cluster_centers"]

    # count distance from each query cell to each reference's cluster centroid
    # [Nq, d]
    X_q = adata_query.obsm[transferred_primary_basis]
    # [1, Nq, d]
    X_q_ = X_q[np.newaxis]
    # [K, Nq, d]
    X_q_K = np.repeat(X_q_, K, axis=0)
    X_weighted_mean = X_weighted_mean.squeeze(2)[:, np.newaxis]
    # [K, Nq, d] = [K, Nq, d] - [K, 1, d]
    adata_query_centered = X_q_K - X_weighted_mean
    # [K, Nq] = np.sum(([K, Nq, d] * [d, d]) X [K, Nq, d], axis=2)
    maha_dists = (
        np.sum(
            np.multiply(adata_query_centered @ inv_cluster_covs, adata_query_centered),
            axis=2,
        )
        ** 0.5
    )

    # [K, Nq]
    Rq = adata_query.obsm[f"{query_basis_adjusted}_symphony_R"].T
    # average distance weighted by cluster membership
    # [Nq] = ([K, Nq] X [K, Nq]).sum()
    adata_query.obs[obs] = np.sum(np.multiply(maha_dists, Rq), axis=0)


def per_cluster_confidence(
    adata_query: AnnData,
    adata_ref: AnnData,
    cluster_key: str,
    u: float = 2,
    lamb: float = 0,
    transferred_primary_basis: str = "X_pca_reference",
    obs: str | None = "symphony_per_cluster_dist",
    uns: str | None = "symphony_per_cluster_dist",
):
    """
    Calculates the Mahalanobis distance from user-defined query clusters to their nearest
    reference centroid after initial projection into reference PCA space.
    All query cells in a cluster get the same score. Higher distance indicates less confidence.
    Due to the instability of estimating covariance with small numbers of cells, we do not assign a
    score to clusters smaller than u * d, where d is the dimensionality of the embedding and u is specified.

    Args:
        adata_query (AnnData): query adata object mapped to `adata_ref` with Symphony
        adata_ref (AnnData): reference adata object (with Harmony object in adata_ref.uns)
        cluster_key (str, [str]): which keys from adata_query.obs to use as a cluster label (if list, adata_query will be grouped by them)
        u (float, optional): at least u * d cells are to be assigned to a cluster, where d is a dimensionality of representation. Defaults to 2.
        lamb (float, optional): ridge regression like coef for covariance matrix inversion numerical stability. Defaults to 0 (no ridge).
        transferred_primary_basis (str, optional): adata_query.obsm[transferred_primary_basis] should contain pre-Symphony reference PC query representation. Defaults to "X_pca_reference".
        obs (str | None, optional): If not None, resulted dists would be written to adata_query.obs[obs] for each cell (just the same value for each cluster)
        uns (str | None, optional): If not None, resulted dists would be written to adata_query.uns[uns] for each cluster
    """
    assert (
        "harmony" in adata_ref.uns
    ), "Harmony object not found in adata_ref.uns['harmony']. First, run standard symphony query mapping."

    harmony = adata_ref.uns["harmony"]

    # [K, d]
    C_ = harmony["C"]  # not normalised
    N = harmony["Nr"]
    reference_cluster_centroids = C_ / N[..., np.newaxis]
    reference_cluster_centroids_norm = reference_cluster_centroids / np.linalg.norm(
        reference_cluster_centroids, ord=2, axis=1, keepdims=True
    )

    dists = adata_query.obs.groupby(cluster_key).apply(
        _cluster_maha_dist,
        adata_query,
        transferred_primary_basis=transferred_primary_basis,
        reference_cluster_centroids=reference_cluster_centroids,
        reference_cluster_centroids_norm=reference_cluster_centroids_norm,
        u=u,
        lamb=lamb,
    )

    if uns is not None:
        adata_query.uns[uns] = {"key": cluster_key, "dist": dists}
    if obs is not None:
        grouped = adata_query.obs.groupby(cluster_key)[cluster_key]
        adata_query.obs[obs] = grouped.transform(lambda x: dists[x[0]])


def ingest(
    adata_query: AnnData,
    adata_ref: AnnData,
    obs: Optional[Union[str, Iterable[str]]] = None,
    embedding_method: Union[str, Iterable[str]] = "umap",
    labeling_method: str = "knn",
    neighbors_key: Optional[str] = None,
    inplace: bool = True,
    use_rep: str | None = None,
    **kwargs,
):
    """
    copied from https://github.com/scverse/scanpy/blob/master/scanpy/tools/_ingest.py
    with little change that var_names equality between adata and adata_new wouldn't be check if needless,
    and additional parameter `use_rep` is added.

    Args:
        adata_query (AnnData): target adata object.
        adata_ref (AnnData): source adata object.
        obs (Optional[Union[str, Iterable[str]]], optional): which columns from adata_red.obs to transfer. Defaults to None.
        embedding_method (Union[str, Iterable[str]], optional): which adata_ref's embeddings to transfer. Defaults to ('umap', 'pca').
        labeling_method (str, optional): which method to use for labeling transferring. Defaults to 'knn'.
        neighbors_key (Optional[str], optional): which key from adata_ref.uns to use as source of neighbors info. Defaults to None.
        inplace (bool, optional): if to write directly to adata_query or return an adjusted copy of adata_query. Defaults to True.
        use_rep (str, None): which adata_query's representation to use for embedding mappings. If None, it will be decided depending on circumstances. Defaults to None.
        kwargs: will be forwarded to the `sc.tl.Ingest.neighbors` function.

    Returns:
        if inplace is False returns a copy of adata_query with additional slots,
        otherwise adds to adata_query.
    """
    anndata_version = pkg_version("anndata")
    if anndata_version < ANNDATA_MIN_VERSION:
        raise ValueError(
            f"ingest only works correctly with anndata>={ANNDATA_MIN_VERSION} "
            f"(you have {anndata_version}) as prior to {ANNDATA_MIN_VERSION}, "
            "`AnnData.concatenate` did not concatenate `.obsm`."
        )

    start = logger.info("running ingest")
    obs = [obs] if isinstance(obs, str) else obs
    embedding_method = (
        [embedding_method] if isinstance(embedding_method, str) else embedding_method
    )
    labeling_method = (
        [labeling_method] if isinstance(labeling_method, str) else labeling_method
    )

    if len(labeling_method) == 1 and len(obs or []) > 1:
        labeling_method = labeling_method * len(obs)

    # dirty hard-code
    if neighbors_key is None:
        neighbors_key = "neighbors"
    if neighbors_key in adata_ref.uns:
        if use_rep is None:
            if "use_rep" not in adata_ref.uns[neighbors_key]["params"]:
                warnings.warn(
                    "'X_pca' representation will be used for neighbors search in adata_query"
                )
                adata_ref.uns[neighbors_key]["params"]["use_rep"] = "X_pca"
            else:
                logger.info(
                    "'%s' representation will be used for neighbors search in adata_query",
                    adata_ref.uns[neighbors_key]["params"]["use_rep"],
                )

    ing = Ingest_sp(adata_ref, neighbors_key, use_representation=use_rep)
    ing.fit(adata_query)

    for method in embedding_method:
        ing.map_embedding(method)

    if obs is not None:
        ing.neighbors(**kwargs)
        for i, col in enumerate(obs):
            ing.map_labels(col, labeling_method[i])

    logger.info("    finished", time=start)

    if "rep" in ing._obsm:
        del ing._obsm["rep"]

    return ing.to_adata(inplace)


def map_embedding(
    adata_query: AnnData,
    adata_ref: AnnData,
    key: list[str] | str | None = None,
    lamb: float | np.array | None = None,
    sigma: float | np.array = 0.1,
    use_genes_column: str | None = "highly_variable",
    transferred_adjusted_basis: str = "X_pca_harmony",
    transferred_primary_basis: str = "X_pca_reference",
    ref_basis_loadings: str = "PCs",
    K: int | None = None,
    reference_primary_basis: str = "X_pca",
) -> None:
    """
    Actually runs Symphony algorithm for mapping adata_query to adata_ref.
    Will use Harmony object from adata_ref.uns if present,
    otherwise will firstly run k-means clusterization Harmony step
    without batch correction.

    Adds mapping of query cells to reference coords
    to adata_query.obsm[transferred_primary_basis] and
    symphony-corrected coords to adata_query.obsm[transferred_adjusted_basis],
    and query to the reference's clusters membership to adata_query.obsm[transferred_adjusted_basis + "_R"]

    Args:
        adata_query (AnnData): query adata object.
        adata_ref (AnnData): reference adata object
            (to account for batch effect in reference first run harmony_integrate).
        key (list[str] | str | None, optional): which of the columns
            from adata_query.obs to consider as batch keys. Defaults to None.
        lamb (float | np.array | None, optional): Entropy regularization parameter for soft k-means. Defaults to None.
        sigma (float | np.array, optional): Ridge regularization parameter for the linear model. Defaults to 0.1.
        use_genes_column (str | None, optional): adata_ref.var[use_genes_column] genes will
            be used to map query embeddings to reference. Defaults to "highly_variable".
        transferred_adjusted_basis (str, optional): in adata_query.obsm[transferred_adjusted_basis]
            symphony-adjusted coords will be saved. Defaults to "X_pca_harmony".
        transferred_primary_basis (str, optional): in adata_query.obsm[transferred_primary_basis]
            adata_query mapping to reference coords will be saved. Defaults to "X_pca_reference".
        ref_basis_loadings (str, optional): adata_ref.varm[ref_basis_loadings] will be used
            as gene loadings to map adata_query to adata_ref coords. Defaults to "PCs".
        K (int | None, optional): Number of clusters to use for k-means clustering.
            Only used if harmony integration was not performed on adata_ref. Defaults to None.
        reference_primary_basis (str): which reference embedding to use for k-means clustering.
            Only used if harmony integration was not performed on adata_ref. Defaults to "X_pca".
    """
    # Errors
    assert (
        "mean" in adata_ref.var
    ), "Gene expression means are expected to be saved in adata_ref.var"
    assert (
        "std" in adata_ref.var
    ), "Gene expression stds are expected to be saved in adata_ref.var"

    if "harmony" not in adata_ref.uns:
        warnings.warn(
            f"Not found `harmony` object in adata_ref.uns.\n"
            f"Assuming that adata_ref doesn't have any batches, "
            f"and using '{reference_primary_basis}' representation of adata_ref for clustering.\n"
            f"Otherwise, firstly run symphonypy.pp.harmony_integrate "
            f"on adata_ref to account for them."
        )

        _run_soft_kmeans(
            adata_ref,
            ref_basis=reference_primary_basis,
            ref_basis_loadings=ref_basis_loadings,
            K=K,
            sigma=sigma,
        )

    if use_genes_column is not None:
        assert (
            use_genes_column in adata_ref.var
        ), f"Column `{use_genes_column}` not found in adata_ref.var. Set `use_genes_column` parameter properly"

    # Warning
    if "log1p" not in adata_query.uns:
        warnings.warn("Gene expressions in adata_query should be log1p-transformed")

    harmony_ref = adata_ref.uns["harmony"]
    ref_basis_loadings = harmony_ref["ref_basis_loadings"]

    # 1. map query to ref initial embedding
    _map_query_to_ref(
        adata_ref=adata_ref,
        adata_query=adata_query,
        transferred_primary_basis=transferred_primary_basis,
        ref_basis_loadings=ref_basis_loadings,
        use_genes_column=use_genes_column,
    )

    # 2. assign clusters
    # [Nq, d]
    X = adata_query.obsm[transferred_primary_basis]
    C = harmony_ref["C"]
    Y = C / np.linalg.norm(C, ord=2, axis=1, keepdims=True)

    R = _assign_clusters(X, sigma, Y, harmony_ref["K"])

    # 3. correct query embeddings
    # likewise harmonypy
    if key is None:
        batch_data = pd.DataFrame({"batch": ["1"] * len(adata_query)})
    elif isinstance(key, str):
        batch_data = pd.DataFrame(adata_query.obs[key]).astype(str)
    else:
        batch_data = adata_query.obs[key].astype(str)

    # [B, N] = [N, B].T  (B -- batch num)
    phi = pd.get_dummies(batch_data).to_numpy().T
    # [B + 1, N]
    phi_ = np.concatenate([np.ones((1, phi.shape[1])), phi], axis=0)

    phi_n = batch_data.describe().loc["unique"].to_numpy().astype(int)
    # lambda (ridge regularization coef)
    if lamb is None:
        lamb = np.repeat([1] * len(phi_n), phi_n)
    elif isinstance(lamb, (float, int)):
        lamb = np.repeat([lamb] * len(phi_n), phi_n)
    elif len(lamb) == len(phi_n):
        lamb = np.repeat([lamb], phi_n)
    else:
        assert len(lamb) == np.sum(phi_n), "each batch variable must have a lambda"

    # [B + 1, B + 1]
    lamb = np.diag(np.insert(lamb, 0, 0))

    adata_query.obsm[transferred_adjusted_basis] = _correct_query(
        X, phi_, R, harmony_ref["Nr"], harmony_ref["C"], lamb
    )
    adata_query.obsm[f"{transferred_adjusted_basis}_symphony_R"] = R.T


def transfer_labels_kNN(
    adata_query: AnnData,
    adata_ref: AnnData,
    ref_labels: list[str] | str,
    *kNN_args,
    query_labels: list[str] | str | None = None,
    ref_basis: str = "X_pca_harmony",
    query_basis: str = "X_pca_harmony",
    **kNN_kwargs,
) -> None:
    """
    Run sklearn kNN classificator for label transferring.
    Args:
        adata_ref (AnnData): adata object to use for train.
        adata_query (AnnData): adata object to use for prediction.
        kNN_args: will be passed to kNN class init function.
        ref_labels (list[str]): columns from adata_ref.obs to use as labels for model training. Default: "X_pca_harmony".
        query_labels (list[str] | None): keys in adata_query.obs where to save transferred ref_labels.
            (in corresponding to ref_labels order). If not provided, same as ref_labels will be used.
        ref_basis (str): adata_ref.obsm[ref_basis] will be used as features for kNN training.
        query_basis (str): adata_query.obsm[query_basis] will be used as features for prediction.
        kNN_kwargs: will be passed to kNN class init function.
    """
    knn = KNeighborsClassifier(*kNN_args, **kNN_kwargs)

    knn.fit(adata_ref.obsm[ref_basis], adata_ref.obs[ref_labels])

    if query_labels is None:
        query_labels = ref_labels

    # TODO: predict_proba
    adata_query.obs[query_labels] = knn.predict(adata_query.obsm[query_basis])


def tsne(
    adata: AnnData,
    use_rep: str = "X_pca",
    t_sne_slot: str = "X_tsne",
    use_model: "openTSNE.TSNEEmbedding" | str | None = None,
    save_path: str | None = None,
    use_raw: bool | None = None,
    return_model: bool = False,
    **kwargs,
) -> None | "openTSNE.TSNEEmbedding":
    """
    Run openTSNE dimension reduction on adata if use_model is None,
    or ingest adata.obsm[use_rep] to existing embedding, saved in use_model.

    Args:
        adata (Anndata): adata object
        use_rep (str): adata.obsm[use_rep] will be used as features for `openTSNE` model.
        t_sne_slot (str): to adata.obsm[t_sne_slot] embedding will be saved.
        use_model (`openTSNE.TSNEEmbedding`, str, None): `openTSNE` model object or path to pickle dumped model to use. Defaults to None.
        save_path (str, None): Filepath to save pickle of the `openTSNE` model. Defaults to None.
        use_raw (bool, None): If to use adata.raw.X as features for `openTSNE`. Defaults to None.
        return_model (bool): If to return `openTSNE` model. Defaults to False.
        kwargs: will be forwarded to the `openTSNE.TSNE` init function.

    Returns:
        if return_model is True, returns `openTSNE` model
    """

    import pickle

    try:
        from openTSNE import TSNE
        from openTSNE import TSNEEmbedding
    except ImportError as exc:
        raise ImportError(
            "\nPlease install openTSNE:\n\n\tpip install openTSNE"
        ) from exc

    if not (use_model is None) and not (save_path is None):
        logger.warning("The model that will be saved is a `PartialTSNEEmbedding`")

    if use_model is None:

        tsne_obj = TSNE(**kwargs)

        if use_rep != "X":
            model = tsne_obj.fit(adata.obsm[use_rep])
        elif use_raw:
            model = tsne_obj.fit(adata.raw.X)
        else:
            model = tsne_obj.fit(adata.X)

        adata.obsm[t_sne_slot] = np.array(model)
    else:

        if isinstance(use_model, str):
            with open(use_model, "rb") as model_file:
                model = pickle.load(model_file)
        elif isinstance(use_model, TSNEEmbedding):
            model = use_model
        else:
            raise Exception(
                "`use_model` should be a path to the model or the model itself."
            )

        if use_rep != "X":
            model = model.transform(adata.obsm[use_rep])
        elif use_raw:
            model = model.transform(adata.raw.X)
        else:
            model = model.transform(adata.X)
        adata.obsm[t_sne_slot] = np.array(model)

    if save_path:
        with open(save_path, "wb") as model_file:
            pickle.dump(model, model_file, protocol=pickle.HIGHEST_PROTOCOL)
            logger.info("Model is saved in %s", save_path)

    if return_model:
        return model
