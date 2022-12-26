# pylint: disable=C0103, W0511, C0114

from typing import List
import anndata as ad
import numpy as np
import scanpy as sc


def preprocess_ref_PCA(
    adata, n_comps: int, batch_keys: List[str], raw_counts: bool, n_top_genes: int
) -> ad.AnnData:
    # TODO: align HVG choice to original Symphony's HVG methods
    # https://github.com/immunogenomics/symphony/blob/main/R/findVariableGenes.R
    # TODO: n_top_genes

    adata.obs["batch_symphonypy"] = (
        (adata.obs[batch_keys]).astype(str).agg("_".join, axis=1)
    )

    if raw_counts:
        adata.layers["counts"] = adata.X.copy()
        sc.pp.normalize_total(adata, target_sum=1e5)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(
            adata,
            batch_key="batch_symphonypy",
            n_top_genes=n_top_genes,
            flavor="seurat_v3",
            layer="counts",
        )

    else:
        sc.pp.highly_variable_genes(
            adata, batch_key="batch_symphonypy", n_top_genes=n_top_genes
        )

    sc.pp.scale(adata, zero_center=True)

    sc.tl.pca(adata, n_comps=n_comps)

    return adata


def preprocess_query_PCA(adata_ref, adata_query, raw_counts: bool) -> ad.AnnData:

    if raw_counts:
        sc.pp.normalize_total(adata_query, target_sum=1e5)
        sc.pp.log1p(adata_query)

    # TODO: adjust for var_names compatibility
    # TODO: check for non-zero std
    # TODO: HVG here intersected by batches?
    t = adata_query[:, adata_ref.var_names[adata_ref.var.highly_variable]].X.copy()
    t = (
        t - np.array(adata_ref.var["mean"][adata_ref.var.highly_variable])[np.newaxis]
    ) / np.array(adata_ref.var["std"][adata_ref.var.highly_variable])[np.newaxis]

    # map query to reference's PCA coords
    # [cells, n_comps] = [cells, genes] * [genes, n_comps]
    adata_query.obsm["X_pca_ref"] = np.array(t * adata_ref.varm["PCs"][adata_ref.var.highly_variable])
