import os
import sys

import scanpy as sc

SYMPHONY_DIR = os.path.dirname(os.path.abspath("__file__"))
SYMPHONYPY_DIR = os.path.join(SYMPHONY_DIR, "symphonypy")
print(SYMPHONY_DIR)
sys.path.append(SYMPHONY_DIR)
import symphonypy as sp


def subsample_adata(adata, group: str, n_sample: int, random_seed: int):
    subsample = adata.obs.groupby(group, group_keys=False).apply(
        lambda x: x.sample(n_sample, random_state=random_seed)
    )
    return adata[adata.obs.index.isin(subsample.index)]


def prepare_test_sample(adata_fp):
    adata = sc.read_h5ad(adata_fp)

    adata_ref = adata[adata.obs["orig.ident"] == "P1_0"]
    adata_ref = subsample_adata(adata_ref, "orig.ident", n_sample_ref, random_state)

    adata_query = adata[adata.obs["orig.ident"].isin(["P2_3", "P3_3"])]
    adata_query = subsample_adata(
        adata_query, "orig.ident", n_sample_query, random_state
    )

    sc.pp.normalize_total(adata_ref, target_sum=1e5)
    sc.pp.log1p(adata_ref)
    sc.pp.highly_variable_genes(
        adata_ref,
        batch_key="orig.ident",
        n_top_genes=3000,
    )
    adata_ref.raw = adata_ref
    adata_ref = adata_ref[:, adata_ref.var.highly_variable]
    sc.pp.scale(adata_ref, max_value=10)
    sc.pp.pca(adata_ref, n_comps=30)
    sp.pp.harmony_integrate(
        adata_ref,
        key="orig.ident",
        verbose=True,
        max_iter_harmony=20,
        random_seed=random_state,
    )
    adata_ref.write("tests/data/PBMC_Satija.harmony.h5ad")

    sc.pp.normalize_total(adata_query, target_sum=1e5)
    sc.pp.log1p(adata_query)
    adata_query.raw = adata_query

    sc.pp.scale(adata_query, max_value=10)
    sc.tl.pca(adata_query, n_comps=30)
    sc.pp.neighbors(adata_query)
    sc.tl.leiden(adata_query)

    adata_query.X = adata_query.raw.X.copy()
    adata_query.write("tests/data/PBMC_Satija.query_subsample.h5ad")

    sp.tl.map_embedding(adata_query=adata_query, adata_ref=adata_ref, key="orig.ident")
    sp.tl.per_cell_confidence(adata_query, adata_ref)
    sp.tl.per_cluster_confidence(adata_query, adata_ref, "leiden")
    adata_query.write("tests/data/PBMC_Satija.symphony.h5ad")


if __name__ == "__main__":
    n_sample_ref = 4096
    n_sample_query = 512
    random_state = 42
    prepare_test_sample("data/PBMC_Satija.h5ad")
