# pylint: disable=E1123, W0621, C0116, W0511, E1121

"""
Symphony algorithm:

1. Reference building:
    - log(CP10K + 1) library size normalization of the cells
    - subset by the top g variable genes
        by the variance stabilizing transform (VST) method (as provided in Seurat18)
    - scaling of the genes to have mean 0 and variance 1 (saving μ and σ for each gene)
    - PCA (by default, d=20) to embed the reference cells
        in a low-dimensional space, saving the gene loadings (U)
    - Harmony (by default, Symphony uses the default parameters
        for the cluster diversity enforcement (θ = 2),
        the entropy regularization hyperparameter for soft k-means (s = 0.1),
        and the number of clusters k = min(100; n/30) )
    savings:
        - gene means (μ) and standard deviations (σ) used to scale the genes
        - PCA gene loadings
        - clusters' centroids
        - Nr and C

2. Soft clustering of query
    - map query to ref (just count coords in PCAs)
    - assign clusters' memberships (soft k-means clustering with regularisation on enthropy)

3. Mixture of experts correction
    - count coefs for linear mixture of experts (matrix inverse)
    - substract linear corrections for embedding by each cluster

4. Label transfering
    Use any unsupervised classification algorithm, e.g. kNN, to transfer labels from ref to query

5*. Count mapping confidence metrics
"""

from typing import List
import anndata as ad
import scanpy as sc

from symphonypy.preprocessing import preprocess_query_PCA, preprocess_ref_PCA
from symphonypy.mapping import Symphony


def run_symphony(
    adata_ref: ad.AnnData,
    adata_query: ad.AnnData,
    batch_keys: List[str],
    n_comps: int,
    harmony_kwargs: dict,
    lamb: float,
    labels: List[str],
    k_neighbours: int,
    raw_counts: bool,
    n_top_genes: int,
) -> None:

    # TODO: adata_ref and adata_query are mutated inside these functions
    preprocess_ref_PCA(
        adata_ref,
        n_comps=n_comps,
        batch_keys=batch_keys,
        raw_counts=raw_counts,
        n_top_genes=n_top_genes,
    )
    preprocess_query_PCA(adata_ref, adata_query, raw_counts=raw_counts)

    # hard-coded for a while
    query_basis = "X_pca_ref"
    adjusted_basis = "X_pca_adjusted"

    # TODO: add an ability to save object to disk
    so = Symphony(adata_ref, batch_keys=batch_keys, harmony_kwargs=harmony_kwargs)

    so.fit()
    so.transform(
        adata_query,
        basis=query_basis,
        adjusted_basis=adjusted_basis,
        batch_keys=batch_keys,
        lamb=lamb,
    )

    so.transfer_labels(
        adata_query, adjusted_basis, adjusted_basis, labels, k_neighbours
    )


if __name__ == "__main__":

    # adata_ref = sc.read_10x_mtx(
    #     "data/filtered_gene_bc_matrices/hg19/",  # the directory with the `.mtx` file
    #     var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
    #     cache=True,
    # )

    # adata_ref = sc.read("data/pancreas.h5ad")  # somehow normalized data

    # adata_ref = sc.read("data/PBMC_Satija.h5ad")  # not normalized, with batch
    # adata_query = adata_ref.copy()

    adata = sc.read("data/exprs_norm_all.h5ad")

    adata_query = adata[adata.obs.donor == "5'"].copy()
    adata_query.obs["ref_query"] = "query"

    adata_ref = adata[~(adata.obs.donor == "5'")].copy()
    adata_ref.obs["ref_query"] = "ref"

    raw_counts = False
    n_comps = 20
    batch_keys = [
        "donor",
    ]
    harmony_kwargs = {"sigma": 0.1}
    lamb = 1
    n_top_genes = 2000
    k_neighbours = 10
    labels = ["cell_type", "cell_subtype"]

    run_symphony(
        adata_ref,
        adata_query,
        batch_keys,
        n_comps,
        harmony_kwargs,
        lamb=lamb,
        labels=labels,
        k_neighbours=k_neighbours,
        raw_counts=raw_counts,
        n_top_genes=n_top_genes,
    )
