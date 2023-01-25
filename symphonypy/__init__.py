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

from . import preprocessing as pp
from . import tools as tl
from . import datasets
