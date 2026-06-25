import logging
import gc
import glob
import math

from pathlib import Path
from typing import Dict, List, Optional, Tuple

import anndata
import numpy as np
import scanpy as sc
import scipy.sparse as sp

logging.basicConfig(
    format="%(levelname)s - %(asctime)s - %(message)s", level=logging.INFO
)


def add_spatial_offset(
    adata: anndata.AnnData,
    sample_keys: Tuple[str, ...] = ("sample", "Sample", "sample_name", "SampleName"),
    spatial_key: str = "spatial",
    new_obsm_key: str = "spatial_offset",
    tile_spacing: float = 300.0,
) -> None:
    """Tile per-sample spatial coordinates into one plotting coordinate system."""
    if spatial_key not in adata.obsm or new_obsm_key in adata.obsm:
        return

    sample_key = next((key for key in sample_keys if key in adata.obs), None)
    if sample_key is None:
        return

    sample_values = adata.obs[sample_key].astype(str)
    samples = sorted(sample_values.unique().tolist())
    n_samples = len(samples)
    if n_samples == 0:
        return

    spatial = np.asarray(adata.obsm[spatial_key])
    if spatial.ndim != 2 or spatial.shape[1] < 2:
        return

    n_cols = min(2, max(1, n_samples))
    n_rows = math.ceil(n_samples / n_cols)
    offset = np.empty_like(spatial, dtype=float)
    grid_bounds = {}
    sample_positions = {}

    for idx, sample_name in enumerate(samples):
        row = idx // n_cols
        col = idx % n_cols
        sample_positions[sample_name] = (row, col)

        mask = (sample_values == sample_name).to_numpy()
        sample_spatial = spatial[mask]
        if sample_spatial.size == 0:
            continue

        max_coords = sample_spatial.max(axis=0)
        min_coords = sample_spatial.min(axis=0)
        grid_bounds[(row, col)] = {
            "width": float(max_coords[0] - min_coords[0]),
            "height": float(max_coords[1] - min_coords[1]),
            "min_x": float(min_coords[0]),
            "max_y": float(max_coords[1]),
        }

    row_heights = [
        max(
            (
                grid_bounds[(row, col)]["height"]
                for col in range(n_cols)
                if (row, col) in grid_bounds
            ),
            default=0.0,
        )
        for row in range(n_rows)
    ]
    col_widths = [
        max(
            (
                grid_bounds[(row, col)]["width"]
                for row in range(n_rows)
                if (row, col) in grid_bounds
            ),
            default=0.0,
        )
        for col in range(n_cols)
    ]

    row_y_offsets = [0.0]
    for idx in range(n_rows - 1):
        row_y_offsets.append(row_y_offsets[-1] - row_heights[idx] - tile_spacing)

    col_x_offsets = [0.0]
    for idx in range(n_cols - 1):
        col_x_offsets.append(col_x_offsets[-1] + col_widths[idx] + tile_spacing)

    for sample_name in samples:
        row, col = sample_positions[sample_name]
        bounds = grid_bounds.get((row, col))
        if bounds is None:
            continue

        mask = (sample_values == sample_name).to_numpy()
        sample_spatial = spatial[mask].copy().astype(float)
        sample_spatial[:, 0] += col_x_offsets[col] - bounds["min_x"]
        sample_spatial[:, 1] += row_y_offsets[row] - bounds["max_y"]
        offset[mask] = sample_spatial

    adata.obsm[new_obsm_key] = offset


def clean_adata(adata: anndata.AnnData) -> anndata.AnnData:
    """Prepare the AnnData object for plotting, particularly with Latch Plots.
    We remove unnesscary obs and reduce the size of the data stored in X.
    """
    obs = [
        "sample", "condition", "sample_name", "cluster",  # atx_snap groups
        "Sample", "SampleName", "Condition",  # ArchRProject groups
        "Clusters", "condition_1",
        "TSSEnrichment", "n_fragment", "log10_nFrags",  # ArchRProject
        "frac_dup", "frac_mito", "on_off", "row",  # atx_snap
        "col", "xcor", "ycor" "tsse"
    ]

    rm_obs = [o for o in adata.obs.keys() if o not in obs]

    adata.obs.drop(rm_obs, axis=1, inplace=True)

    adata.varm.clear()
    adata.layers.clear()

    if adata.raw:
        adata.raw = None

    add_spatial_offset(adata)

    rm_obsm = [obsm for obsm in adata.obsm.keys()
               if obsm not in ["spatial_offset", "X_umap"]]

    for obsm in rm_obsm:
        del adata.obsm[obsm]

    # Sparse to dense conversion
    if sp.issparse(adata.X):
        try:
            adata.X = adata.X.toarray()
        except Exception as e:
            logging.warning(f"Could not convert sparse .X to dense: {e}")

    try:
        adata.X = adata.X.astype(np.float16)
    except Exception as e:
        logging.warning(f"Cannot convert .X to float16: {e}")

    return adata


def clean_index_columns(*adatas: anndata.AnnData) -> None:
    """Remove _index columns from AnnData objects if they exist."""
    for adata in adatas:
        if hasattr(adata, 'raw') and adata.raw is not None:
            if "_index" in adata.raw.var.columns:
                adata.raw.var.drop(columns=['_index'], inplace=True)


def load_and_combine_data() -> Tuple[anndata.AnnData, anndata.AnnData]:
    """Load and combine gene and motif AnnData objects."""
    logging.info("Reading and combining gene AnnData...")
    gene_files = glob.glob("*g_converted.h5ad")
    gene_adatas = [anndata.read_h5ad(file) for file in gene_files]
    adata_gene = sc.concat(gene_adatas)

    logging.info("Reading and combining motif AnnData...")
    motif_files = glob.glob("*m_converted.h5ad")
    motif_adatas = [anndata.read_h5ad(file) for file in motif_files]
    adata_motif = sc.concat(motif_adatas)

    # Clean up memory
    del gene_adatas, motif_adatas
    gc.collect()

    # Clean up index columns if they exist
    clean_index_columns(adata_gene, adata_motif)

    return adata_gene, adata_motif


def save_anndata_objects(
    adata_gene: anndata.AnnData,
    adata_motif: anndata.AnnData,
    base_dir: Path
) -> None:
    """Save full and reduced AnnData objects."""
    logging.info("Saving full adata...")

    add_spatial_offset(adata_gene)
    add_spatial_offset(adata_motif)

    # Save full objects
    adata_gene.X = adata_gene.X.astype(np.float32)
    adata_gene.write(f"{base_dir}/combined_ge.h5ad")
    adata_motif.write(f"{base_dir}/combined_motifs.h5ad")

    # Create and save reduced objects
    logging.info("Making reduced gene adata...")
    sm_adata_gene = clean_adata(adata_gene)
    sm_adata_motif = clean_adata(adata_motif)

    logging.info("Saving gene adata...")
    sm_adata_gene.write(f"{base_dir}/combined_sm_ge.h5ad")

    logging.info("Saving motif adata...")
    sm_adata_motif.write(f"{base_dir}/combined_sm_motifs.h5ad")


def load_analysis_results(
    adata_gene: anndata.AnnData,
    adata_motif: anndata.AnnData,
    groups: List[str]
) -> None:
    """Load various analysis results into AnnData objects."""
    # Load differential analysis results
    load_marker_genes(adata_gene)
    load_enriched_motifs(adata_motif)

    # Load heatmap data
    load_heatmaps(adata_gene, adata_motif)

    # Load volcano plots if condition analysis was performed
    if "condition" in groups:
        load_volcano_plots(adata_gene, adata_motif)


def load_csv_files_to_uns(
    pattern: str,
    target_uns: Dict,
    dtype_spec: Optional[Dict] = None,
    index_col: Optional[int] = None,
    name_transform: Optional[callable] = None
) -> None:
    import pandas as pd
    """Generic function to load CSV files matching a pattern into uns
    dictionary."""
    try:
        files = glob.glob(pattern)
        for file in files:
            name = Path(file).stem
            if name_transform:
                name = name_transform(name, file)
            df = pd.read_csv(file, dtype=dtype_spec, index_col=index_col)
            target_uns[name] = df
    except Exception as e:
        logging.warning(f"Error loading files matching {pattern}: {e}")


def load_enriched_motifs(adata_motif: anndata.AnnData) -> None:
    """Load enriched motifs data."""
    logging.info("Adding enriched motifs...")
    load_csv_files_to_uns(
        "enrichedMotifs_*.csv", adata_motif.uns, dtype_spec={"group_name": str}
    )


def load_heatmaps(
    adata_gene: anndata.AnnData, adata_motif: anndata.AnnData
) -> None:
    """Load heatmap data for both genes and motifs."""
    logging.info("Adding heatmaps...")

    # Gene heatmaps
    load_csv_files_to_uns(
        "genes_per_*_hm.csv",
        adata_gene.uns,
        dtype_spec={"cluster": str},
        index_col=0
    )

    # Motif heatmaps
    load_csv_files_to_uns(
        "motif_per_*_hm.csv",
        adata_motif.uns,
        index_col=0
    )


def load_marker_genes(adata_gene: anndata.AnnData) -> None:
    """Load ranked genes data."""
    logging.info("Adding ranked genes...")
    load_csv_files_to_uns(
        "ranked_genes_*.csv", adata_gene.uns, dtype_spec={"group_name": str}
    )


def load_volcano_plots(
    adata_gene: anndata.AnnData, adata_motif: anndata.AnnData
) -> None:
    """Load volcano plot data for genes and motifs."""
    logging.info("Adding gene volcanos...")
    load_csv_files_to_uns(
        "volcanoMarkers_genes_*.csv",
        adata_gene.uns,
        dtype_spec={"cluster": str},
        name_transform=lambda name, file: volcano_name_transform(name, file, "genes")
    )

    logging.info("Adding motif volcanos...")
    load_csv_files_to_uns(
        "volcanoMarkers_motifs_*.csv",
        adata_motif.uns,
        dtype_spec={"cluster": str},
        name_transform=lambda name,
        file: volcano_name_transform(name, file, "motifs")
    )


def transfer_embedding_data(
    adata_gene: anndata.AnnData,
    adata_motif: anndata.AnnData,
    data_path: str,
    obsm_key: str
) -> None:
    """Transfer embedding data (UMAP/spatial) to AnnData objects."""
    import pandas as pd

    logging.info(f"Transferring {obsm_key}...")

    try:
        df = pd.read_csv(data_path, index_col=0)
        aligned_data_g = df.loc[adata_gene.obs_names].values
        aligned_data_m = df.loc[adata_motif.obs_names].values

        adata_gene.obsm[obsm_key] = aligned_data_g
        adata_motif.obsm[obsm_key] = aligned_data_m

    except (FileNotFoundError, KeyError) as e:
        logging.warning(f"Error loading {obsm_key} data: {e}")


def volcano_name_transform(name: str, file: str, data_type: str) -> str:
    """Transform volcano plot file names."""
    file_name = Path(file).name
    treatment = file_name.replace(f"volcanoMarkers_{data_type}_", "").replace(".csv", "")
    return f"volcano_{treatment}"
