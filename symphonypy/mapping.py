# pylint: disable=C0103, C0116, C0114, C0115, W0511

from typing import List
import anndata as ad
import numpy as np
import pandas as pd

from harmonypy import run_harmony, Harmony
from sklearn.neighbors import KNeighborsClassifier


class Symphony:
    def __init__(
        self,
        adata_ref: ad.AnnData,
        batch_keys: List[str],
        harmony_kwargs: dict,
    ) -> None:

        self.adata_ref = adata_ref
        self.batch_keys = batch_keys
        self.harmony_kwargs = harmony_kwargs

        # harmony area
        self.Nr = None
        self.C = None
        self.Y = None
        self.sigma = None
        self.K = None

        # symphony area
        self.phi_ = None
        self.lamb = None
        self.R = None

    def fit(
        self, ref_basis_source="X_pca", ref_basis_adjusted="X_pca_adjusted"
    ) -> None:
        ref_ho = run_harmony(
            self.adata_ref.obsm[ref_basis_source],
            meta_data=self.adata_ref.obs,
            vars_use=self.batch_keys,
            **self.harmony_kwargs
        )  # type: Harmony

        self.adata_ref.obsm[ref_basis_adjusted] = ref_ho.Z_corr.T

        # [K] the number of cells softly belonging to each cluster
        self.Nr = ref_ho.R.sum(axis=1)

        # [K, d] = [K, Nref] x [d, N_ref].T
        self.C = np.dot(ref_ho.R, ref_ho.Z_corr.T)

        # ref cluster centroids L2 normalized
        # [K, d] = [d, K].T
        self.Y = ref_ho.Y.T

        # number of clusters
        self.K = ref_ho.K

        # sigma [K] (cluster cross enthropy regularization coef)
        self.sigma = ref_ho.sigma

    def transform(
        self,
        adata_query: ad.AnnData,
        basis: str,
        adjusted_basis: str,
        batch_keys: List[str],
        lamb: float,
    ) -> ad.AnnData:

        # [Nq, d]
        X = adata_query.obsm[basis]
        self.assign_clusters(X)

        # likewise harmonypy
        batch_data = adata_query.obs[batch_keys]
        # [B, N] = [N, B].T  (B -- batch num)
        phi = pd.get_dummies(batch_data).to_numpy().T
        # [B + 1, N]
        self.phi_ = np.concatenate([np.ones((1, phi.shape[1])), phi], axis=0)

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
        self.lamb = np.diag(np.insert(lamb, 0, 0))

        adata_query.obsm[adjusted_basis] = self.correct_query(X)

    def assign_clusters(self, X):
        # it's made so in harmonypy, maybe to prevent overflow during L2 normalization?
        X_cos = X / X.max(axis=1, keepdims=True)
        # L2 normalization for cosine distance
        X_cos = X_cos / np.linalg.norm(X_cos, ord=2, axis=1, keepdims=True)

        # [K, N] = [K, d] x [Nq, d].T
        self.R = -2 * (1 - np.dot(self.Y, X_cos.T)) / self.sigma[..., np.newaxis]
        self.R = np.exp(self.R)
        self.R /= self.R.sum(axis=0, keepdims=True)

    def correct_query(self, X):
        # [d, N] = [N, d].T
        X_corr = X.copy().T

        for i in range(self.K):
            # [B + 1, N] = [B + 1, N] * [N]
            Phi_Rk = np.multiply(self.phi_, self.R[i, :])

            # [B + 1, B + 1] = [B + 1, N] x [N, B + 1]
            x = np.dot(Phi_Rk, self.phi_.T)
            # += [1]
            x[0, 0] += self.Nr[i]

            # [B + 1, d] = [B + 1, N] x [N, d]
            y = np.dot(Phi_Rk, X)
            y[0, :] += self.C[i]

            # [B + 1, d] = [B + 1, B + 1] x [B + 1, d]
            W = np.dot(np.linalg.inv(x + self.lamb), y)
            W[0, :] = 0  # do not remove the intercept

            # [d, N] -= [B + 1, d].T x [B + 1, N]
            X_corr -= np.dot(W.T, Phi_Rk)

        return X_corr.T

    def transfer_labels(
        self,
        adata_query: ad.AnnData,
        ref_basis: str,
        query_basis: str,
        labels: List[str],
        k_neighbours: int,
    ) -> None:
        knn = KNeighborsClassifier(k_neighbours, weights="distance")
        knn.fit(self.adata_ref.obsm[ref_basis], self.adata_ref.obs[labels])
        adata_query.obs[labels] = knn.predict(adata_query.obsm[query_basis])

    # def compute_confidence()
