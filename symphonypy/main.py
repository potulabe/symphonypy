# pylint: disable=E1123, W0621, C0116, W0511, E1121

from typing import List
import anndata as ad
import scanpy as sc

import symphonypy as sp


def run_symphony(
    adata_ref: ad.AnnData,
    adata_query: ad.AnnData,
    batch_keys: List[str],
    n_comps: int,
    harmony_args: List,
    harmony_kwargs: dict,
    lamb: float,
    labels: List[str],
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
    # HVG, PCA
    sp.utils.preprocess_ref_PCA(
        adata_ref,
        n_comps=n_comps,
        batch_keys=batch_keys,
        raw_counts=raw_counts,
        n_top_genes=n_top_genes,
        search_highly_variable=search_highly_variable,
    )

    ref_basis_source = "X_pca"
    basis_adjusted = "X_pca_adjusted"

    ref_basis_loadings = "PCs"
    query_basis_ref = "X_pca_ref"

    # preprocess query
    if raw_counts:
        sc.pp.normalize_total(adata_query, target_sum=1e5)
        sc.pp.log1p(adata_query)

    sp.pp.harmony_integrate(
        adata_ref,
        ref_basis_source=ref_basis_source,
        ref_basis_adjusted=basis_adjusted,
        ref_basis_loadings=ref_basis_loadings,
        vars_use=batch_keys,
        *harmony_args,  # TODO: test if it works
        **harmony_kwargs
    )

    sp.tl.map_embedding(
        adata_ref,
        adata_query,
        batch_keys=batch_keys,
        lamb=lamb,
        use_genes_column=use_genes_column,
        adjusted_basis_query=basis_adjusted,
        query_basis_ref=query_basis_ref,
    )

    sp.tl.transfer_labels_kNN(
        adata_ref,
        adata_query,
        labels,
        labels,
        # kNN args
        n_neighbours,
        ref_basis=basis_adjusted,
        query_basis=basis_adjusted,
        # kNN kwargs
        weights="distance",
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
    n_neighbours = 10
    labels = ["cell_type", "cell_subtype"]
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
        harmony_args=[],
        harmony_kwargs=harmony_kwargs,
        use_genes_column=use_genes_column,
    )
