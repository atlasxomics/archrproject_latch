import anndata
import logging
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from matplotlib.backends.backend_pdf import PdfPages
from pathlib import Path
from typing import Iterable, List, Optional

from atx_common import filter_anndata


GROUPED_NHOOD_ENRICHMENT_KEY = "cluster_nhood_enrichment_by_group"
GROUPED_NHOOD_SCHEMA_VERSION = 1


def _ensure_obs_categorical(
    adata: anndata.AnnData, key: Optional[str]
) -> None:
    if key is None or key not in adata.obs.columns:
        return

    if adata.obs[key].dtype.name == "category":
        adata.obs[key] = adata.obs[key].cat.remove_unused_categories()
    else:
        adata.obs[key] = adata.obs[key].astype("category")


def _nhood_result_for_subgroup(
    adata: anndata.AnnData,
    mask: np.ndarray,
    cluster_key: str,
    sample_key: Optional[str],
    spatial_key: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute enrichment without copying the feature matrix from ``adata``."""
    obs_keys = [cluster_key]
    if sample_key is not None and sample_key != cluster_key:
        obs_keys.append(sample_key)

    subset = anndata.AnnData(obs=adata.obs.loc[mask, obs_keys].copy())
    subset.obsm[spatial_key] = np.asarray(adata.obsm[spatial_key])[mask].copy()
    _ensure_obs_categorical(subset, cluster_key)
    _ensure_obs_categorical(subset, sample_key)

    categories = subset.obs[cluster_key].cat.categories.astype(str).to_numpy()
    if subset.n_obs < 2:
        shape = (len(categories), len(categories))
        return categories, np.zeros(shape, dtype=float), np.zeros(shape, dtype=float)

    squidpy_analysis(
        subset,
        cluster_key=cluster_key,
        sample_key=sample_key,
        spatial_key=spatial_key,
    )
    result = subset.uns[f"{cluster_key}_nhood_enrichment"]
    return categories, np.asarray(result["zscore"]), np.asarray(result["count"])


def precompute_grouped_nhood_enrichment(
    adata: anndata.AnnData,
    group_keys: Iterable[str],
    cluster_key: str = "cluster",
    sample_key: Optional[str] = None,
    spatial_key: Optional[str] = None,
) -> None:
    """Store subgroup neighborhood results in a backed-AnnData-friendly schema.

    The result contains only small enrichment matrices and their labels.  It is
    deliberately indexed with numeric HDF5-safe keys because user-provided group
    values can contain slashes and other characters that are unsafe as HDF5 paths.
    """
    if spatial_key is None:
        spatial_key = next(
            (key for key in ("spatial_offset", "spatial") if key in adata.obsm),
            None,
        )
    if spatial_key is None:
        raise KeyError(
            "Spatial coordinates were not found in `adata.obsm`; expected "
            "`spatial_offset` or `spatial`."
        )

    stored_groups = {}
    seen = set()
    for group_key in group_keys:
        if group_key in seen or group_key == cluster_key:
            continue
        seen.add(group_key)
        if group_key not in adata.obs:
            logging.warning(
                "Skipping neighborhood precomputation for missing obs key '%s'.",
                group_key,
            )
            continue

        values = pd.unique(adata.obs[group_key].dropna())
        stored_subgroups = {}
        for subgroup_index, group_value in enumerate(values):
            mask = (adata.obs[group_key] == group_value).to_numpy()
            logging.info(
                "Precomputing neighborhood enrichment for %s=%s (%d cells).",
                group_key,
                group_value,
                int(mask.sum()),
            )
            categories, zscore, count = _nhood_result_for_subgroup(
                adata,
                mask,
                cluster_key,
                sample_key,
                spatial_key,
            )
            stored_subgroups[str(subgroup_index)] = {
                "group_value": str(group_value),
                "cluster_categories": categories,
                "zscore": zscore,
                "count": count,
            }

        stored_groups[str(len(stored_groups))] = {
            "group_key": str(group_key),
            "subgroups": stored_subgroups,
        }

    adata.uns[GROUPED_NHOOD_ENRICHMENT_KEY] = {
        "schema_version": GROUPED_NHOOD_SCHEMA_VERSION,
        "cluster_key": cluster_key,
        "groups": stored_groups,
    }


def plot_neighborhoods(
    adata: anndata.AnnData,
    group: str,
    subgroups: Optional[List[str]],
    outdir: str = "figures"
):
    from squidpy.pl import nhood_enrichment

    if group != "all" and subgroups:
        filtered_adatas = {}
        for sg in subgroups:
            filtered_adata = filter_anndata(adata, group, sg)
            squidpy_analysis(filtered_adata)
            filtered_adatas[sg] = filtered_adata

    plt.rcParams.update({'figure.autolayout': True})
    with PdfPages(f"{outdir}/{group}_neighborhoods.pdf") as pdf:

        if subgroups:
            for sg in subgroups:
                fig = nhood_enrichment(
                    filtered_adatas[sg],
                    cluster_key="cluster",
                    method="single",
                    title=f"{group} {sg}: Neighborhood enrichment",
                    cmap="bwr",
                    vmin=-50,
                    vmax=50,
                )
                pdf.savefig(fig, bbox_inches="tight")
                plt.close(fig)

        elif group == "all":
            fig = nhood_enrichment(
                adata,
                cluster_key="cluster",
                method="single",
                title="All cells: Neighborhood enrichment",
                cmap="bwr",
                vmin=-50,
                vmax=50,
            )
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)


def run_squidpy_analysis(
    adata_gene: anndata.AnnData,
    figures_dir: Path,
    sample_key: Optional[str] = None,
    group_keys: Optional[Iterable[str]] = None,
) -> anndata.AnnData:
    """Run Squidpy analysis, including results consumed by backed Plots."""
    from squidpy.gr import ripley as compute_ripley
    from squidpy.pl import ripley as plot_ripley
    import scanpy as sc

    logging.info("Running squidpy...")
    adata_gene = squidpy_analysis(adata_gene, sample_key=sample_key)
    precompute_grouped_nhood_enrichment(
        adata_gene,
        group_keys=group_keys or (),
        sample_key=sample_key,
    )

    # Generate neighborhood plots
    logging.info("Making neighborhood plots...")
    group_dict = {"all": None}

    for group_name, group_value in group_dict.items():
        plot_neighborhoods(
            adata_gene, group_name, group_value, outdir=str(figures_dir)
        )

    # Generate Ripley's plot
    logging.info("Running ripley...")
    compute_ripley(
        adata_gene,
        cluster_key="cluster",
        mode="L",
        max_dist=500,
    )
    old_figdir = sc.settings.figdir
    try:
        sc.settings.figdir = str(figures_dir)
        plot_ripley(
            adata_gene,
            cluster_key="cluster",
            mode="L",
            save="ripleys_L.pdf",
        )
    finally:
        sc.settings.figdir = old_figdir

    return adata_gene


def run_spatial_autocorr(
    adata: anndata.AnnData,
    n_jobs: int = 4,
) -> pd.DataFrame:
    """Compute Moran's I per feature using squidpy. Returns DataFrame sorted by I descending."""
    import squidpy as sq

    if "spatial_connectivities" not in adata.obsp:
        sq.gr.spatial_neighbors(adata, coord_type="grid", n_neighs=4, n_rings=1)

    layer = None
    for candidate in ["log1p", "lognorm", "normalized"]:
        if candidate in adata.layers:
            layer = candidate
            break

    features_upper = pd.Index(adata.var_names.astype(str)).str.upper()
    keep = ~(
        features_upper.str.startswith("MT-")
        | features_upper.str.startswith("RPS")
        | features_upper.str.startswith("RPL")
        | features_upper.str.startswith("MTRNR")
    )
    if "highly_variable" in adata.var.columns:
        keep = keep & adata.var["highly_variable"].to_numpy()

    test_features = adata.var_names[keep].tolist()
    if not test_features:
        raise ValueError("No features remain after filtering for spatial autocorrelation.")

    logging.info(
        "Running spatial autocorrelation (Moran's I) on %d features (layer=%s).",
        len(test_features),
        layer or "X",
    )
    sq.gr.spatial_autocorr(
        adata, mode="moran", genes=test_features, layer=layer, n_perms=None, n_jobs=n_jobs
    )
    return adata.uns["moranI"].sort_values("I", ascending=False)


def _get_page_saver(output_path: str):
    import os
    ext = os.path.splitext(output_path)[1].lower()
    page_idx = [1]
    pdf = PdfPages(output_path) if ext == ".pdf" else None
    output_stem = os.path.splitext(output_path)[0] if ext else output_path
    output_ext = ext if ext else ".png"
    saved_paths: List[str] = []

    def save_page(fig):
        if pdf is not None:
            pdf.savefig(fig)
            return
        image_path = f"{output_stem}_{page_idx[0]:03d}{output_ext}"
        fig.savefig(image_path, dpi=200, bbox_inches="tight")
        saved_paths.append(image_path)
        page_idx[0] += 1

    def close():
        if pdf is not None:
            pdf.close()

    return save_page, close, pdf is not None, saved_paths


def _write_html_gallery(output_path: str, title: str, image_paths: List[str], captions: Optional[List[str]] = None, html_output_path: Optional[str] = None) -> None:
    import os
    if not image_paths:
        return
    html_path = html_output_path if html_output_path is not None else f"{os.path.splitext(output_path)[0]}.html"
    html_dir = os.path.dirname(html_path) or "."
    relative_paths = [os.path.relpath(p, start=html_dir) for p in image_paths]
    blocks = []
    for idx, rel_path in enumerate(relative_paths, start=1):
        caption = captions[idx - 1] if captions and idx - 1 < len(captions) else f"Page {idx}"
        blocks.append(f"<section class=\"page\"><h2>{caption}</h2><img loading=\"lazy\" src=\"{rel_path}\" alt=\"{caption}\" /></section>")
    body = "\n".join(blocks)
    html = f"""<!doctype html><html lang="en"><head><meta charset="utf-8" /><title>{title}</title>
<style>body{{margin:24px auto;max-width:1600px;background:#f4f4f4;font-family:Arial,sans-serif}}.page{{margin:0 0 24px;background:#fff;border:1px solid #ddd;border-radius:8px;padding:12px}}img{{width:100%;height:auto;display:block}}</style>
</head><body><h1>{title}</h1>{body}</body></html>"""
    with open(html_path, "w", encoding="utf-8") as fh:
        fh.write(html)


def plot_svg_spatial(
    adata: anndata.AnnData,
    svg_df: pd.DataFrame,
    figures_dir: Path,
    filename: str,
    modality: str,
    top_n: int = 10,
    html_output_path: Optional[str] = None,
) -> None:
    """Spatial scatter plots for the top N spatially variable features."""
    import scipy.sparse as sparse_mod
    from squidpy.pl import spatial_scatter

    if svg_df.empty or "spatial" not in adata.obsm:
        return

    top_features = [f for f in svg_df.index[:top_n] if f in adata.var_names]
    if not top_features:
        return

    layer = None
    for candidate in ["log1p", "lognorm", "normalized"]:
        if candidate in adata.layers:
            layer = candidate
            break

    feat_indices = adata.var_names.get_indexer(top_features)
    X_src = adata.layers[layer] if layer is not None else adata.X
    if sparse_mod.issparse(X_src):
        X_sub = X_src[:, feat_indices].toarray().astype(np.float32)
    else:
        X_sub = np.asarray(X_src[:, feat_indices], dtype=np.float32)

    plot_adata = anndata.AnnData(
        X=X_sub,
        obs=adata.obs.copy(),
        var=adata.var.iloc[feat_indices].copy(),
        obsm={"spatial": adata.obsm["spatial"].copy()},
    )

    # Resolve sample column (ArchR uses "Sample", atx_snap uses "sample")
    sample_key = next(
        (k for k in ["sample", "Sample", "sample_name", "SampleName"] if k in adata.obs.columns),
        None,
    )
    samples = (
        sorted(adata.obs[sample_key].astype(str).unique())
        if sample_key is not None
        else ["all"]
    )

    n_cols = 5
    n_rows = max(1, -(-len(top_features) // n_cols))

    output_path = str(figures_dir / filename)
    page_captions: List[str] = []
    save_page, close, is_pdf, image_paths = _get_page_saver(output_path)
    try:
        for sample in samples:
            if sample_key is not None:
                sample_adata = plot_adata[plot_adata.obs[sample_key].astype(str) == sample].copy()
            else:
                sample_adata = plot_adata.copy()

            fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 4, n_rows * 4), squeeze=False)
            axes_flat = axes.ravel()

            for i, feat in enumerate(top_features):
                ax = axes_flat[i]
                moran_i = float(svg_df.at[feat, "I"]) if feat in svg_df.index else float("nan")
                try:
                    spatial_scatter(
                        sample_adata,
                        color=feat,
                        size=75,
                        shape=None,
                        ax=ax,
                        title=f"{feat}  I={moran_i:.3f}",
                        colorbar=False,
                    )
                    if ax.collections:
                        fig.colorbar(ax.collections[0], ax=ax, shrink=0.6)
                    ax.axis("off")
                except Exception as feat_err:
                    logging.warning("Could not plot %s: %s", feat, feat_err)
                    ax.set_title(feat)
                    ax.axis("off")

            for ax in axes_flat[len(top_features):]:
                ax.axis("off")

            plt.suptitle(
                f"{sample} {modality} — top {len(top_features)} spatially variable features",
                fontsize=13,
                fontweight="bold",
                y=1.02,
            )
            plt.tight_layout()
            page_captions.append(f"{sample}: top {len(top_features)} SVFs")
            save_page(fig)
            plt.close(fig)
    finally:
        close()

    if not is_pdf:
        _write_html_gallery(
            output_path,
            title=f"Top Spatially Variable {modality}",
            image_paths=image_paths,
            captions=page_captions,
            html_output_path=html_output_path,
        )


def squidpy_analysis(
    adata: anndata.AnnData,
    cluster_key: str = "cluster",
    sample_key: Optional[str] = None,
    spatial_key: Optional[str] = None,
) -> anndata.AnnData:
    """Perform squidpy Neighbors enrichment analysis.
    """
    from squidpy.gr import nhood_enrichment, spatial_neighbors

    _ensure_obs_categorical(adata, cluster_key)
    _ensure_obs_categorical(adata, sample_key)

    n_clusters = len(adata.obs[cluster_key].cat.categories)
    if spatial_key is None:
        spatial_key = next(
            (key for key in ("spatial_offset", "spatial") if key in adata.obsm),
            None,
        )
    if spatial_key is None:
        raise KeyError(
            "Spatial coordinates were not found in `adata.obsm`; expected "
            "`spatial_offset` or `spatial`."
        )

    spatial_neighbors(
        adata,
        spatial_key=spatial_key,
        coord_type="grid",
        n_neighs=4,
        n_rings=1,
        library_key=sample_key,
    )
    if n_clusters < 2:
        logging.warning(
            "Skipping Squidpy neighborhood enrichment because only "
            f"{n_clusters} cluster is present."
        )
        adata.uns[f"{cluster_key}_nhood_enrichment"] = {
            "zscore": np.zeros((n_clusters, n_clusters), dtype=float),
            "count": np.zeros((n_clusters, n_clusters), dtype=float),
            "skipped": "fewer_than_two_clusters",
        }
        return adata

    nhood_enrichment(
        adata,
        cluster_key=cluster_key,
        library_key=sample_key,
        seed=42,
    )

    return adata
