<!-- omit in toc -->
# Symphonypy
Porting of [Symphony R](https://github.com/immunogenomics/symphony) package to Python

- [Installation](#installation)
- [Usage examples](#examples)
- [Instructions](#instructions)
- [Benchmarking](#benchmarking)


## Installation
Symphonypy package might be installed via pip:
```
pip install symphonypy
```

## Examples
Here are Jupyter-notebooks with simple examples of how to use symphonypy
- with [pre-built references](notebooks/Symphonypy_precomputed.ipynb) from original Symphony,
- with [reference building](notebooks/Symphonypy_simple_tutorial.ipynb) from scratch.

## Instructions
### Step 1: reference building
```
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
sp.pp.harmony_integrate(adata_ref, key=batch_key_ref)
```

### Step 2: Query preprocessing and Symphony
```
# target_sum for normalize_total() should be the same as in reference
sc.pp.normalize_total(adata_query, target_sum=1e5)
sc.pp.log1p(adata_query)
sp.tl.map_embedding(adata_query, adata_ref, key=batch_key_query)
```

### Step 3: Label transfer
```
sp.tl.transfer_labels_kNN(adata_query, adata_ref, labels)
```

### Step 4 (optional): Dimensionality reduction
With UMAP:
```
sc.pp.neighbors(adata_ref, use_rep="X_pca_harmony")
sc.tl.umap(adata_ref)
sp.tl.ingest(adata_query, adata_ref)
```
With t-SNE (`openTSNE` should be installed, `pip install openTSNE`):
```
tSNE_model = sp.tl.tsne(adata_ref, use_rep="X_pca_harmony", return_model=True)
sp.tl.tsne(adata_query, use_rep="X_pca_harmony", use_model=tSNE_model)
```

## Benchmarking
- Harmony (R) vs harmonypy benchmarking:
  - [benchmarking/Benchmarking_harmony_PBMC_Satija.ipynb](benchmarking/Benchmarking_harmony_PBMC_Satija_CITEseq.ipynb)
- Symphony (R) vs symphonypy benchmarking: [benchmarking/Benchmarking_symphony_PBMC.ipynb](benchmarking/Benchmarking_symphony_PBMC.ipynb)
- PBMC example from the Symphony repo: [benchmarking/validation_PBMC_example.ipynb](benchmarking/validation_PBMC_example.ipynb)

Download data used in benchmarking: [benchmarking/data_download.ipynb](benchmarking/data_download.ipynb)
