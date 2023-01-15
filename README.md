<!-- omit in toc -->
# Symphonypy
Porting of [Symphony R](https://github.com/immunogenomics/symphony) package to Python

- [Usage](#usage)
  - [Preprocessing](#preprocessing)
  - [Run symphony](#run-symphony)
  - [Transfer labels](#transfer-labels)
  - [Map UMAP](#map-umap)
  - [Map Open tSNE](#map-open-tsne)
- [Benchmarking](#benchmarking)

> Currently under development:
> - evaluation of confidence metrics
> - building reference without running Harmony
> - Symphony R via rpy2 envelop
> - precomputed symphony reference datasets

  
## Usage
### Preprocessing
```python
import scanpy as sc
import symphonypy as sp


n_comps = 20
batch_key = "donor"
n_top_genes = 2000
n_neighbours = 10


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


# create reference embedding, e.g. PCA:
sc.pp.scale(adata_ref, zero_center=True, max_value=10)
adata_ref.X[adata_ref.X < -10] = -10 # for R Symphony-like processing
sc.tl.pca(adata_ref, n_comps=n_comps, zero_center=False)


# preprocess query in the same way as reference:
sc.pp.normalize_total(adata_query, target_sum=1e5)
sc.pp.log1p(adata_query)
```
### Run symphony
```python
harmony_kwargs = {"sigma": 0.1}
lamb = 1


# run Harmonypy or Harmony R on the reference:
sp.pp.harmony_integrate(
    adata_ref,
    ref_basis_source="X_pca",
    ref_basis_adjusted="X_pca_harmony",
    ref_basis_loadings="PCs",
    flavor="python",  # might be "R" (will run Harmony via rpy2)
    key=batch_keys,
    **harmony_kwargs,
)

# run symphonypy to map query to the reference's embedding:
sp.tl.map_embedding(
    adata_ref,
    adata_query,
    key=batch_key,  # could be list of batch keys
    lamb=lamb,
    use_genes_column="highly_variable",
    adjusted_basis_query="X_pca_harmony",
    query_basis_ref="X_pca_reference",
)
```
### Transfer labels
```python
labels = ["cell_type", "cell_subtype"]  # any columns from adata_ref.obs


# transfer labels via scipy kNN
sp.tl.transfer_labels_kNN(
    adata_ref,
    adata_query,
    labels,
    # kNN args
    n_neighbours,
    ref_basis="X_pca_harmony",
    query_basis="X_pca_harmony",
    # kNN kwargs
    weights="distance",
)
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
sp.tl.ingest(adata=adata_query, adata_ref=adata_ref, embedding_method="umap")
```
### Map Open tSNE
```python
# map query to the reference's tSNE
tSNE = sc.tl.tsne(adata_ref, use_rep="X_pca_harmony", return_model=True)
sc.tl.tsne(adata_query, use_model=tSNE)
```

## Benchmarking
- Harmony (R) vs harmonypy benchmarking: [benchmarking/Benchmarking_harmony.ipynb](benchmarking/Benchmarking_harmony.ipynb)
- Symphony (R) vs symphonypy benchmarking: [benchmarking/Benchmarking_symphony.ipynb](benchmarking/Benchmarking_symphony.ipynb)
- PBMC example from the Symphony repo: [benchmarking/validation_PBMC_example.ipynb](benchmarking/validation_PBMC_example.ipynb)

Download data used in benchmarking: [benchmarking/data_download.ipynb](benchmarking/data_download.ipynb)
