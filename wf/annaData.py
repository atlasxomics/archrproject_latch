import logging
import pandas as pd
import scanpy as sc

logging.basicConfig(
    format="%(levelname)s - %(asctime)s - %(message)s", level=logging.INFO
)

logging.info("Reading combined matrix...")
counts = pd.read_csv("combined_mat.csv", index_col=0).T

logging.info("Creating AnnData object...")
adata = sc.AnnData(counts)

logging.info("Normalizing counts...")
sc.pp.normalize_total(adata, target_sum=10000)
sc.pp.log1p(adata)

logging.info("Identifying variable genes...")
sc.pp.highly_variable_genes(
    adata, flavor="seurat_v3", n_top_genes=2000, batch_key=None, inplace=True
)

logging.info("Writing h5ad to disk...")
adata.write("/root/combined_hvg.h5ad")
