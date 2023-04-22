import scanpy as sc
import symphonypy as sp


class TestPreprocessing:
    adata_fp = "tests/data/PBMC_Satija.query_subsample.h5ad"
    batch_key = "orig.ident"

    def assert_harmony_object(self, adata):
        assert "X_pca_harmony" in adata.obsm
        assert "harmony" in adata.uns
        assert "Nr" in adata.uns["harmony"]
        assert "C" in adata.uns["harmony"]
        assert "K" in adata.uns["harmony"]
        assert "sigma" in adata.uns["harmony"]
        assert "ref_basis_loadings" in adata.uns["harmony"]
        assert "ref_basis_adjusted" in adata.uns["harmony"]
        assert "vars_use" in adata.uns["harmony"]
        assert "harmony_kwargs" in adata.uns["harmony"]
        assert "converged" in adata.uns["harmony"]
        assert "R" in adata.uns["harmony"]

    def test_harmony_integrate_python(self):
        adata = sc.read_h5ad(self.adata_fp)
        sp.pp.harmony_integrate(adata, key=self.batch_key, flavor="python")

        self.assert_harmony_object(adata)
        

    # def test_harmony_integrate_R(self):
    #     adata = sc.read_h5ad(self.adata_fp)
    #     sp.pp.harmony_integrate(adata, key=self.batch_key, flavor="R")

    #     self.assert_harmony_object(adata)
