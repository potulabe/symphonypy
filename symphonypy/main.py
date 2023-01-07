# pylint: disable=E1123, W0621, C0116, W0511, E1121

from __future__ import annotations

import scanpy as sc
from anndata import AnnData

import symphonypy as sp


def run_symphony(
    adata_ref: AnnData,
    adata_query: AnnData,
    batch_keys: list[str] | str,
    n_comps: int,
    harmony_kwargs: dict,
    lamb: float,
    labels: list[str],
    n_neighbours: int,
    raw_counts: bool,
    n_top_genes: int,
    use_genes_column: str,
) -> None:
    """
    This function is supposed to be used mostly for debugging
    1. preprocessing
        - preprocessing_ref(adata_ref)
             -> adata_ref.uns["harmony"]
        - preprocessing_query(adata_query, adata_ref)
             -> adata_query["PCA_ref"]
    2. run symphony
        - adjust PCA
        - transfer labels
        - transfer UMAP
    """
    search_highly_variable = (
        use_genes_column == "highly_variable" and "highly_variable" not in adata_ref.var
    )

    if isinstance(batch_keys, str):
        batch_keys = [
            batch_keys,
        ]

    # HVG, PCA
    adata_ref.obs["batch_symphonypy"] = (
        (adata_ref.obs[batch_keys]).astype(str).agg("_".join, axis=1)
    )

    if raw_counts:
        if search_highly_variable:
            sc.pp.highly_variable_genes(
                adata_ref,
                batch_key="batch_symphonypy",
                n_top_genes=n_top_genes,
                flavor="seurat_v3",
            )
        sc.pp.normalize_total(adata_ref, target_sum=1e5)
        sc.pp.log1p(adata_ref)

    elif search_highly_variable:
        sc.pp.highly_variable_genes(
            adata_ref, batch_key="batch_symphonypy", n_top_genes=n_top_genes
        )

    max_value = 10
    sc.pp.scale(adata_ref, zero_center=True, max_value=max_value)
    adata_ref.X[adata_ref.X < -max_value] = -max_value

    sc.tl.pca(adata_ref, n_comps=n_comps, zero_center=False)

    ref_basis_source = "X_pca"
    basis_adjusted = "X_pca_harmony"

    ref_basis_loadings = "PCs"
    query_basis_ref = "X_pca_reference"

    # preprocess query
    if raw_counts:
        sc.pp.normalize_total(adata_query, target_sum=1e5)
        sc.pp.log1p(adata_query)

    sp.pp.harmony_integrate(
        adata_ref,
        ref_basis_source=ref_basis_source,
        ref_basis_adjusted=basis_adjusted,
        ref_basis_loadings=ref_basis_loadings,
        flavor="R",
        key=batch_keys,
        **harmony_kwargs,
    )

    sp.tl.map_embedding(
        adata_ref,
        adata_query,
        key=batch_keys,
        lamb=lamb,
        use_genes_column=use_genes_column,
        adjusted_basis_query=basis_adjusted,
        query_basis_ref=query_basis_ref,
    )

    sp.tl.transfer_labels_kNN(
        adata_ref,
        adata_query,
        labels,
        # kNN args
        n_neighbours,
        ref_basis=basis_adjusted,
        query_basis=basis_adjusted,
        # kNN kwargs
        weights="distance",
    )

    sc.pp.neighbors(
        adata_ref, n_pcs=n_comps, n_neighbors=20, knn=True, use_rep=basis_adjusted
    )
    sc.tl.umap(adata_ref)

    sc.tl.ingest(adata=adata_query, adata_ref=adata_ref, embedding_method="umap")


if __name__ == "__main__":

    # adata_ref = sc.read_10x_mtx(
    #     "data/filtered_gene_bc_matrices/hg19/",  # the directory with the `.mtx` file
    #     var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
    #     cache=True,
    # )

    # adata_ref = sc.read("data/pancreas.h5ad")  # somehow normalized data

    adata_ref = sc.read("data/PBMC_Satija.h5ad")  # not normalized, with batch
    adata_query = adata_ref.copy()
    raw_counts = True
    labels = [
        "celltype.l1",
    ]

    # adata = sc.read("data/exprs_norm_all.h5ad")

    # adata_query = adata[adata.obs.donor == "5'"].copy()
    # adata_query.obs["ref_query"] = "query"
    # adata_ref = adata[~(adata.obs.donor == "5'")].copy()
    # adata_ref.obs["ref_query"] = "ref"
    # raw_counts = False
    # labels = ["cell_type", "cell_subtype"]

    n_comps = 20
    batch_keys = "donor"
    harmony_kwargs = {"sigma": 0.1}
    lamb = 1
    n_top_genes = 2000
    n_neighbours = 10
    use_genes_column = "highly_variable"

    run_symphony(
        adata_ref=adata_ref,
        adata_query=adata_query,
        batch_keys=batch_keys,
        n_comps=n_comps,
        lamb=lamb,
        labels=labels,
        n_neighbours=n_neighbours,
        raw_counts=raw_counts,
        n_top_genes=n_top_genes,
        harmony_kwargs=harmony_kwargs,
        use_genes_column=use_genes_column,
    )
