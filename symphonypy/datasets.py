from __future__ import annotations

from scanpy import read, AnnData
from pathlib import Path


def pbmcs_10x_reference(
    file_path: str | Path = "data/symphony_ref/pbmcs_10x_reference.h5ad",
) -> AnnData:
    url = "https://zenodo.org/record/7607565/files/pbmcs_10x_reference.h5ad"
    adata = read(file_path, backup_url=url, sparse=True, cache=True)
    return adata

  
def pancreatic_atlas(
    file_path: str | Path = "data/symphony_ref/pancreas_plate-based_reference.h5ad",
) -> AnnData:
    url = "https://zenodo.org/record/7607565/files/pancreas_plate-based_reference.h5ad"
    adata = read(file_path, backup_url=url, sparse=True, cache=True)
    return adata
  
  
def fetal_liver(
    file_path: str | Path = "data/symphony_ref/fetal_liver_reference_3p.h5ad",
) -> AnnData:
    url = "https://zenodo.org/record/7607565/files/fetal_liver_reference_3p.h5ad"
    adata = read(file_path, backup_url=url, sparse=True, cache=True)
    return adata
  
  
def kidney(
    file_path: str | Path = "data/symphony_ref/kidney_healthy_fetal_reference.h5ad",
) -> AnnData:
    url = "https://zenodo.org/record/7607565/files/kidney_healthy_fetal_reference.h5ad"
    adata = read(file_path, backup_url=url, sparse=True, cache=True)
    return adata
  
  
def t_cells(
    file_path: str | Path = "data/symphony_ref/tbru_ref.h5ad",
) -> AnnData:
    url = "https://zenodo.org/record/7607565/files/tbru_ref.h5ad"
    adata = read(file_path, backup_url=url, sparse=True, cache=True)
    return adata
  

def inflammatory_atlas(
    file_path: str | Path = "data/symphony_ref/inflammatory_reference.h5ad",
) -> AnnData:
    url = "https://zenodo.org/record/7607565/files/zhang_reference.h5ad"
    adata = read(file_path, backup_url=url, sparse=True, cache=True)
    return adata
  
  
def TMS(
    file_path: str | Path = "data/symphony_ref/TMS_facs_reference.h5ad",
) -> AnnData:
    url = "https://zenodo.org/record/7607565/files/TMS_facs_reference.h5ad"
    adata = read(file_path, backup_url=url, sparse=True, cache=True)
    return adata
