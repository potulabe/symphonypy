<!-- omit in toc -->
# Symphonypy
[![pytest workflow](https://github.com/potulabe/symphonypy/actions/workflows/test.yaml/badge.svg)](https://github.com/potulabe/symphonypy/actions/workflows/test.yaml)

Porting of [Symphony R](https://github.com/immunogenomics/symphony) package to Python

- [Installation](#installation)
- [Examples](#examples)
- [Instructions](#instructions)
  - [Step 1: Reference building](#step-1-reference-building)
  - [Step 2: Query preprocessing and Symphony](#step-2-query-preprocessing-and-symphony)
  - [Step 3: Label transfer](#step-3-label-transfer)
  - [Step 4 (optional): Dimensionality reduction](#step-4-optional-dimensionality-reduction)
- [Benchmarking](#benchmarking)
- [Short API description](#short-api-description)


## Installation
Symphonypy package might be installed via pip:
```
pip install symphonypy
```

## Examples
Here are Jupyter-notebooks with simple examples of how to use symphonypy
- with [pre-built references](notebooks/Symphonypy_precomputed.ipynb) from original Symphony ([nbviewer](https://nbviewer.org/github/potulabe/symphonypy/blob/main/notebooks/Symphonypy_precomputed.ipynb)) — *references contain only data that is necessary for label transfer and don't include gene expressions themselves*,
- with [reference building](notebooks/Symphonypy_simple_tutorial.ipynb) from scratch ([nbviewer](https://nbviewer.org/github/potulabe/symphonypy/blob/main/notebooks/Symphonypy_simple_tutorial.ipynb)),
- for [mapping to reference without harmony step](notebooks/Symphonypy_without_harmony_tutorial.ipynb) ([nbviewer](https://nbviewer.org/github/potulabe/symphonypy/blob/main/notebooks/Symphonypy_without_harmony_tutorial.ipynb)).



## Instructions
### Step 1: Reference building
```python
import scanpy as sc
import symphonypy as sp

sc.pp.normalize_total(adata_ref, target_sum=1e5)
sc.pp.log1p(adata_ref)
sc.pp.highly_variable_genes(
    adata_ref,
    batch_key=batch_key_ref,
    n_top_genes=n_top_genes,
)
adata_ref.raw = adata_ref
adata_ref = adata_ref[:, adata_ref.var.highly_variable]
sc.pp.scale(adata_ref, max_value=10)
sc.pp.pca(adata_ref, n_comps=30, zero_center=False)

# You can skip Harmony if you have only one batch in reference
sp.pp.harmony_integrate(adata_ref, key=batch_key_ref)  
# -> adata_ref.obsm["X_pca_harmony"] <- Harmony adjusted "X_pca"
# -> adata_ref.uns["harmony"] <- Harmony object for Symphony
```

### Step 2: Query preprocessing and Symphony
```python
# target_sum for normalize_total() should be the same as in reference
sc.pp.normalize_total(adata_query, target_sum=1e5)
sc.pp.log1p(adata_query)

# Symphony
sp.tl.map_embedding(adata_query, adata_ref, key=batch_key_query)
# -> adata_query.obsm["X_pca_harmony"] <- Symphony adjusted query's PCA
sp.tl.per_cell_confidence(adata_query, adata_ref)
# -> adata_query.obs["symphony_per_cell_dist"] <- Symphony mapping score per cell
sp.tl.per_cluster_confidence(adata_query, adata_ref, query_clusters)
# -> adata_query.uns["symphony_per_cluster_dist"] <- Symphony mapping score per cluster
```

### Step 3: Label transfer
```python
sp.tl.transfer_labels_kNN(adata_query, adata_ref, labels)
# -> adata_query.obs[labels] <- transferred labels (via sklearn kNN)
```

### Step 4 (optional): Dimensionality reduction
With UMAP:
```python
sc.pp.neighbors(adata_ref, use_rep="X_pca_harmony")
sc.tl.umap(adata_ref)
sp.tl.ingest(adata_query, adata_ref)
# -> adata_query.obsm["X_umap"] <- mapped to the reference's UMAP coords
```

With t-SNE (`openTSNE` should be installed, `pip install openTSNE`):
```python
tSNE_model = sp.tl.tsne(adata_ref, use_rep="X_pca_harmony", return_model=True)
sp.tl.tsne(adata_query, use_rep="X_pca_harmony", use_model=tSNE_model)
# -> adata_query.obsm["X_tsne"] <- mapped to the reference's tSNE coords
```

## Benchmarking
- Harmony (R) vs harmonypy benchmarking: [benchmarking/Benchmarking_harmony_PBMC_Satija.ipynb](benchmarking/Benchmarking_harmony_PBMC_Satija_CITEseq.ipynb)
- Symphony (R) vs symphonypy benchmarking: [benchmarking/Benchmarking_symphony_PBMC.ipynb](benchmarking/Benchmarking_symphony_PBMC.ipynb)

## Short API description
Arguments for each function are listed in a source code.

| Function | Description |
|-|-|
|[`datasets`](symphonypy/datasets.py)|Preprocessed datasets for label transfer|
|[`pp.harmony_integrate`](preprocessing.py#L13)|Run Harmony batch correction on adata, save corrected output to `adata.obsm`, save all the necessary to Symphony mapping algorithm parameters to `adata.uns`|
|[`tl.map_embedding`](symphonypy/tools.py#L257)|Actually runs Symphony algorithm for mapping `adata_query` to `adata_ref`.|
|[`tl.transfer_labels_kNN`](symphonypy/tools.py#L390)|Run sklearn kNN classificator for label transferring.|
|[`tl.per_cell_confidence`](symphonypy/tools.py#L33)|Calculates the weighted Mahalanobis distance for query cells to reference clusters.|
|[`tl.per_cluster_confidence`](symphonypy/tools.py#L109)|Calculates the Mahalanobis distance from user-defined query clusters to their nearest reference centroid after initial projection into reference PCA space.|
|[`tl.ingest`](symphonypy/tools.py#L171)|Copied from https://github.com/scverse/scanpy/blob/master/scanpy/tools/_ingest.py with little change that var_names equality between adata and adata_new wouldn't be check if needless, and additional parameter `use_rep` is added.|
|[`tl.tsne`](symphonypy/tools.py#L427)|Run openTSNE dimension reduction on adata if `use_model` is None, or ingest `adata.obsm[use_rep]` to existing embedding, saved in `use_model`.|
