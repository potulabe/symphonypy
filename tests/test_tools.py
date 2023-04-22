import scanpy as sc
import symphonypy as sp


class TestTools:
    adata_query_fp = "tests/data/PBMC_Satija.query_subsample.h5ad"
    adata_query_res_fp = "tests/data/PBMC_Satija.symphony.h5ad"
    adata_ref_fp = "tests/data/PBMC_Satija.harmony.h5ad"
    query_batch_key = "orig.ident"
    query_clusters = "leiden"

    @staticmethod
    def assert_equals(f, s, threshold=1e-7):
        assert (abs(f - s) < threshold).all()

    def test_map_embedding(self):
        adata_query = sc.read_h5ad(self.adata_query_fp)
        adata_ref = sc.read_h5ad(self.adata_ref_fp)
        sp.tl.map_embedding(
            adata_query=adata_query, adata_ref=adata_ref, key=self.query_batch_key
        )

        adata_query_res = sc.read_h5ad(self.adata_query_res_fp)

        self.assert_equals(
            adata_query.obsm["X_pca_reference"], adata_query_res.obsm["X_pca_reference"]
        )
        self.assert_equals(
            adata_query.obsm["X_pca_harmony"], adata_query_res.obsm["X_pca_harmony"]
        )

        sp.tl.per_cell_confidence(adata_query, adata_ref)

        self.assert_equals(
            adata_query.obs["symphony_per_cell_dist"],
            adata_query_res.obs["symphony_per_cell_dist"],
        )

        sp.tl.per_cluster_confidence(adata_query, adata_ref, self.query_clusters)

        self.assert_equals(
            adata_query.obs["symphony_per_cluster_dist"],
            adata_query_res.obs["symphony_per_cluster_dist"],
        )

        self.assert_equals(
            adata_query.uns["symphony_per_cluster_dist"]["dist"],
            adata_query_res.uns["symphony_per_cluster_dist"]["dist"],
        )

        assert (
            adata_query.uns["symphony_per_cluster_dist"]["cluster_labels"]
            == adata_query_res.uns["symphony_per_cluster_dist"]["cluster_labels"]
        ).all()
