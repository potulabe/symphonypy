# pylint: disable=C0103, C0116, C0114, C0115, W0511
from __future__ import annotations

import logging
import warnings

from typing import Iterable, Optional, Union

import numpy as np
import pandas as pd

from anndata import AnnData
from packaging import version
from sklearn.neighbors import KNeighborsClassifier
from scipy.sparse import issparse

from scanpy.tools._ingest import Ingest, _DimDict
from scanpy._compat import pkg_version

from ._utils import (
    _assign_clusters,
    _correct_query,
    _map_query_to_ref,
    _adjust_for_missing_genes,
    _run_soft_kmeans,
)


ANNDATA_MIN_VERSION = version.parse("0.7rc1")
logger = logging.getLogger("symphonypy")


class Ingest_sp(Ingest):
    def fit(self, adata_new: AnnData):
        self._obs = pd.DataFrame(index=adata_new.obs.index)
        self._obsm = _DimDict(adata_new.n_obs, axis=0)

        self._adata_new = adata_new
        self._obsm["rep"] = self._same_rep()

    def _same_rep(self):
        adata = self._adata_new

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

        use_genes_list_present = np.isin(use_genes_list, self._adata_new.var_names)

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


def ingest(
    adata_query: AnnData,
    adata_ref: AnnData,
    obs: Optional[Union[str, Iterable[str]]] = None,
    embedding_method: Union[str, Iterable[str]] = "umap",
    labeling_method: str = "knn",
    neighbors_key: Optional[str] = None,
    inplace: bool = True,
    **kwargs,
):
    """
    copied from https://github.com/scverse/scanpy/blob/master/scanpy/tools/_ingest.py
    with little change that var_names equality between adata and adata_new wouldn't be check if needless

    Args:
        adata_query (AnnData): _description_
        adata_ref (AnnData): _description_
        obs (Optional[Union[str, Iterable[str]]], optional): _description_. Defaults to None.
        embedding_method (Union[str, Iterable[str]], optional): _description_. Defaults to ('umap', 'pca').
        labeling_method (str, optional): _description_. Defaults to 'knn'.
        neighbors_key (Optional[str], optional): _description_. Defaults to None.
        inplace (bool, optional): _description_. Defaults to True.

    Raises:
        ValueError: _description_

    Returns:
        _type_: _description_
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

    ing = Ingest_sp(adata_ref, neighbors_key)
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
    adjusted_basis_query: str = "X_pca_harmony",
    query_basis_ref: str = "X_pca_reference",
    ref_basis_loadings: str = "PCs",
    K: int | None = None,
) -> None:
    """
    Actually runs Symphony algorithm for mapping adata_query to adata_ref.
    Will use Harmony object from adata_ref.uns if present,
    otherwise will firstly run k-means clusterization Harmony step
    without batch correction.

    Adds mapping of query cells to reference coords
    to adata_query.obsm[query_basis_ref] and
    symphony-corrected coords to adata_query.obsm[adjusted_basis_query].

    Args:
        adata_ref (AnnData): reference adata,
            to account for batch effect in reference
            first run harmony_integrate
        adata_query (AnnData): query adata
        key (list[str] | str | None, optional): which of the columns
            from adata_query.obs to consider as batch keys.
            Defaults to None.
        lamb (float | np.array | None, optional): Entropy regularization parameter for soft k-means. Defaults to None.
        sigma (float | np.array, optional): Ridge regularization parameter for the linear model. Defaults to 0.1.
        use_genes_column (str | None, optional): adata_ref.var[use_genes_column] genes will
            be used to map query embeddings to reference. Defaults to "highly_variable".
        adjusted_basis_query (str, optional): in adata_query.obsm[adjusted_basis_query]
            symphony-adjusted coords will be saved. Defaults to "X_pca_harmony".
        query_basis_ref (str, optional): in adata_query.obsm[query_basis_ref]
            adata_query mapping to reference coords will be saved. Defaults to "X_pca_reference".
        ref_basis_loadings (str, optional): adata_ref.varm[ref_basis_loadings] will be used
            as gene loadings to map adata_query to adata_ref coords. Defaults to "PCs".
        K (int | None, optional): Number of clusters to use for k-means clustering.
            Only used if harmony integration was not performed on adata_ref. Defaults to None.
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
            "Not found `harmony` object in adata_ref.uns. "
            "Assuming that adata_ref doesn't have any batches. "
            "Otherwise, firstly run symphonypy.pp.harmony_integrate "
            "on adata_ref to account for them."
        )
        _run_soft_kmeans(
            adata_ref,
            ref_basis=adjusted_basis_query,
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
        query_basis_ref=query_basis_ref,
        ref_basis_loadings=ref_basis_loadings,
        use_genes_column=use_genes_column,
    )

    # 2. assign clusters
    # [Nq, d]
    X = adata_query.obsm[query_basis_ref]

    R = _assign_clusters(X, sigma, harmony_ref["Y"], harmony_ref["K"])

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

    adata_query.obsm[adjusted_basis_query] = _correct_query(
        X, phi_, R, harmony_ref["Nr"], harmony_ref["C"], lamb
    )


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
    Args:
        adata_ref (AnnData): _description_
        adata_query (AnnData): _description_
        ref_basis (str): _description_
        query_basis (str): _description_
        ref_labels (list[str]): keys from adata_ref.obs to transfer to adata_query. Default: "X_pca_harmony"
        query_labels (list[str] | None): keys in adata_query.obs where to save transferred ref_labels.
            (in corresponding to ref_labels order). If not provided, same as ref_labels will be used
        n_neighbors (int): kNN parameter
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
            print(f"Model is saved in {save_path}")

    if return_model:
        return model
