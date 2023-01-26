<!-- omit in toc -->
# Symphonypy
Porting of [Symphony R](https://github.com/immunogenomics/symphony) package to Python

- [Installation](#installation)
- [Usage](#usage)
  - [Preprocessing](#preprocessing)
  - [Run symphony](#run-symphony)
  - [Transfer labels](#transfer-labels)
  - [Map UMAP](#map-umap)
  - [Map tSNE with `openTSNE`](#map-tsne-with-opentsne)
- [Benchmarking](#benchmarking)

> Currently under development:
> - evaluation of confidence metrics
> - precomputed symphony reference datasets
> - Python package


## Installation
Symphonypy package might be installed via pip:
```
pip install symphonypy
```

  
## Usage
### Preprocessing
```python
import scanpy as sc
import symphonypy as sp


n_comps = 20
batch_key = "donor"
n_top_genes = 2000


# preprocess reference, e.g. HVG, normalize, log1p:
sc.pp.highly_variable_genes(
    adata_ref,
    batch_key=batch_key,
    n_top_genes=n_top_genes,
    flavor="seurat_v3",
)
sc.pp.normalize_total(adata_ref, target_sum=1e5)
sc.pp.log1p(adata_ref)
adata.raw = adata


# preprocess query in the same way as reference:
sc.pp.normalize_total(adata_query, target_sum=1e5)
sc.pp.log1p(adata_query)


# create reference embedding, e.g. PCA:
sc.pp.scale(adata_ref, zero_center=True, max_value=10)
adata_ref.X[adata_ref.X < -10] = -10 # for R Symphony-like preprocessing
sc.tl.pca(adata_ref, n_comps=n_comps, zero_center=False)  # do not center again
```
### Run symphony
```python
harmony_kwargs = {"sigma": 0.1}
lamb = 1


# run Harmonypy or Harmony R,
# if you want to integrate the reference dataset
# (skip this, if there is no necessity 
# of batch correction in the reference):
sp.pp.harmony_integrate(
    adata_ref,
    ref_basis_source="X_pca",
    ref_basis_adjusted="X_pca_harmony",
    ref_basis_loadings="PCs",
    flavor="python",  # might be "R" (will run Harmony via rpy2)
    key=batch_keys,
    **harmony_kwargs,
)
# -> adata_ref.obsm["X_pca_harmony"]
# -> adata_ref.uns["harmony"] (needed for the next step)


# run symphonypy to map query to the reference's embedding:
sp.tl.map_embedding(
    adata_query,
    adata_ref,
    key=batch_key,  # could be list of batch keys
    lamb=lamb,
    use_genes_column="highly_variable",
    adjusted_basis_query="X_pca_harmony",
    query_basis_ref="X_pca_reference",
)
# -> adata_query.obsm["X_pca_harmony"]
```
### Transfer labels
```python
n_neighbours = 10
labels = ["cell_type", "cell_subtype"]  # any columns from adata_ref.obs
transferred_labels = ["predicted_cell_type", "predicted_cell_subtype"]

# transfer labels via scipy kNN
sp.tl.transfer_labels_kNN(
    adata_query,
    adata_ref,
    labels,
    # kNN args
    n_neighbours,
    query_labels=transferred_labels,
    ref_basis="X_pca_harmony",
    query_basis="X_pca_harmony",
    # kNN kwargs
    weights="distance",
)
# -> adata_query.obs[:, transferred_labels]
```
### Map UMAP
```python
# map query to the reference's UMAP
# build reference UMAP
sc.pp.neighbors(
    adata_ref,
    n_pcs=n_comps,
    n_neighbors=n_neighbours,
    knn=True,
    use_rep="X_pca_harmony"
)
sc.tl.umap(adata_ref)
# run ingest (same as sc.tl.ingest, but with setting to zero expressions of var_names missed in query)
sp.tl.ingest(adata_query=adata_query, adata_ref=adata_ref, embedding_method="umap")
# -> adata_query.obsm["X_umap"]
```
### Map tSNE with `openTSNE`
```python
# map query to the reference's tSNE
tSNE = sp.tl.tsne(adata_ref, use_rep="X_pca_harmony", return_model=True)
sp.tl.tsne(adata_query, use_model=tSNE, use_rep="X_pca_harmony")
# -> adata_query.obsm["X_tsne"]
```

## Benchmarking
- Harmony (R) vs harmonypy benchmarking:
  - [benchmarking/Benchmarking_harmony_PBMC_Satija.ipynb](benchmarking/Benchmarking_harmony_PBMC_Satija.ipynb)
  - [benchmarking/Benchmarking_harmony_PBMC_Satija_CITEseq.ipynb](benchmarking/Benchmarking_harmony_PBMC_Satija_CITEseq.ipynb)
- Symphony (R) vs symphonypy benchmarking: [benchmarking/Benchmarking_symphony_PBMC.ipynb](benchmarking/Benchmarking_symphony_PBMC.ipynb)
- PBMC example from the Symphony repo: [benchmarking/validation_PBMC_example.ipynb](benchmarking/validation_PBMC_example.ipynb)

Download data used in benchmarking: [benchmarking/data_download.ipynb](benchmarking/data_download.ipynb)
