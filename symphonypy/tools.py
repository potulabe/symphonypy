# pylint: disable=C0103, C0116, C0114, C0115, W0511

import logging
from typing import List, Union
import numpy as np
import pandas as pd
import scanpy as sc

from sklearn.neighbors import KNeighborsClassifier
from anndata import AnnData

from ._utils import _assign_clusters, _correct_query, _map_query_to_ref


logger = logging.getLogger("symphonypy")


def map_embedding(
    adata_ref: AnnData,
    adata_query: AnnData,
    key: Union[List[str], str, None] = None,
    lamb: Union[float, np.array, None] = None,
    use_genes_column: str = "highly_variable",
    adjusted_basis_query: str = "X_pca_harmony",
    query_basis_ref: str = "X_pca_ref",
) -> None:

    assert (
        "mean" in adata_ref.var
    ), "Gene expression means are expected to be saved in adata_ref.var"
    assert (
        "std" in adata_ref.var
    ), "Gene expression stds are expected to be saved in adata_ref.var"
    assert (
        "harmony" in adata_ref.uns
    ), "Firstly run symphonypy.pp.harmony_integrate on adata_ref"
    assert (
        use_genes_column in adata_ref.var
    ), f"Column `{use_genes_column}` not found in adata_ref.var"

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
    R = _assign_clusters(X, harmony_ref["sigma"], harmony_ref["Y"])

    # 3. correct query embeddings
    # likewise harmonypy
    if key is None:
      batch_data = pd.DataFrame({"batch" : ["1"] * len(adata_query)})
    elif type(key) == str:
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
        X, phi_, R, harmony_ref["K"], harmony_ref["Nr"], harmony_ref["C"], lamb
    )


def transfer_labels_kNN(
    adata_ref: AnnData,
    adata_query: AnnData,
    ref_labels: List[str],
    query_labels: Union[List[str], None],
    *kNN_args,
    ref_basis: str = "X_pca_harmony",
    query_basis: str = "X_pca_harmony",
    **kNN_kwargs,
) -> None:
    """

    Args:
        adata_ref (ad.AnnData): _description_
        adata_query (ad.AnnData): _description_
        ref_basis (str): _description_
        query_basis (str): _description_
        ref_labels (list[str]): keys from adata_ref.obs to transfer to adata_query
        query_labels (list[str] | None): keys in adata_query.obs where to save transferred ref_labels
            (in corresponding to ref_labels order)
            If not provided, same as ref_labels will be used
        n_neighbors (int): kNN parameter

    """
    knn = KNeighborsClassifier(*kNN_args, **kNN_kwargs)

    knn.fit(adata_ref.obsm[ref_basis], adata_ref.obs[ref_labels])

    # TODO: predict_proba
    adata_query.obs[query_labels] = knn.predict(adata_query.obsm[query_basis])
