<!-- omit in toc -->
# Symphonypy
Porting of [Symphony R](https://github.com/immunogenomics/symphony) package to Python

- [Installation](#installation)
- [Examples](#examples)
- [Instructions](#instructions)
  - [Step 1: Reference building](#step-1-reference-building)
  - [Step 2: Query preprocessing and Symphony](#step-2-query-preprocessing-and-symphony)
  - [Step 3: Label transfer](#step-3-label-transfer)
  - [Step 4 (optional): Dimensionality reduction](#step-4-optional-dimensionality-reduction)
- [Benchmarking](#benchmarking)


## Installation
Symphonypy package might be installed via pip:
```
pip install symphonypy
```

## Examples
Here are Jupyter-notebooks with simple examples of how to use symphonypy
- with [pre-built references](notebooks/Symphonypy_precomputed.ipynb) from original Symphony,
- with [reference building](notebooks/Symphonypy_simple_tutorial.ipynb) from scratch,
- for [mapping to reference without harmony step](notebooks/Symphonypy_without_harmony_tutorial.ipynb).

## Instructions
### Step 1: Reference building
```python
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
