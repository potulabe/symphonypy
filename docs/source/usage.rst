Usage
=====

.. _installation:

Installation
------------

Symphonypy package might be installed via pip:

.. code-block:: console

   pip install symphonypy

Instructions
------------

Step 1: Reference building
**************************

.. code-block:: python

   import scanpy as sc
   import symphonypy as sp

   # adata_ref is a filtered dataset that contains counts
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
   # zero_center=False for sc.pp.pca() is recommended, but even
   # without it the results will be pretty same
   
   # You can skip Harmony if you have only one batch in reference
   sp.pp.harmony_integrate(adata_ref, key=batch_key_ref)  
   # -> adata_ref.obsm["X_pca_harmony"] <- Harmony adjusted "X_pca"
   # -> adata_ref.uns["harmony"] <- Harmony object for Symphony

Step 2: Query preprocessing and Symphony
****************************************

.. code-block:: python

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

Step 3: Label transfer
**********************
.. code-block:: python

   sp.tl.transfer_labels_kNN(adata_query, adata_ref, labels)
   # -> adata_query.obs[labels] <- transferred labels (via sklearn kNN)

Step 4 (optional): Dimensionality reduction
*******************************************

With UMAP:

.. code-block:: python

   sc.pp.neighbors(adata_ref, use_rep="X_pca_harmony")
   sc.tl.umap(adata_ref)
   sp.tl.ingest(adata_query, adata_ref)
   # -> adata_query.obsm["X_umap"] <- mapped to the reference's UMAP coords

With t-SNE (``openTSNE`` should be installed, ``pip install openTSNE``):

.. code-block:: python

   tSNE_model = sp.tl.tsne(adata_ref, use_rep="X_pca_harmony", return_model=True)
   sp.tl.tsne(adata_query, use_rep="X_pca_harmony", use_model=tSNE_model)
   # -> adata_query.obsm["X_tsne"] <- mapped to the reference's tSNE coords
