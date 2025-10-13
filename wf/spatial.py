import anndata
import logging
import matplotlib.pyplot as plt

from matplotlib.backends.backend_pdf import PdfPages
from pathlib import Path
from typing import List, Optional

from wf.utils import filter_anndata


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
    adata_gene: anndata.AnnData, figures_dir: Path
) -> anndata.AnnData:
    """Run Squidpy analysis and generate plots."""
    from squidpy.pl import ripley

    logging.info("Running squidpy...")
    adata_gene = squidpy_analysis(adata_gene)

    # Generate neighborhood plots
    logging.info("Making neighborhood plots...")
    group_dict = {"all": None}

    for group_name, group_value in group_dict.items():
        plot_neighborhoods(
            adata_gene, group_name, group_value, outdir=str(figures_dir)
        )

    # Generate Ripley's plot
    logging.info("Running ripley...")
    ripley(adata_gene, cluster_key="cluster", mode="L", save="ripleys_L.pdf")

    return adata_gene


def squidpy_analysis(
    adata: anndata.AnnData, cluster_key: str = "cluster"
) -> anndata.AnnData:
    """Perform squidpy Neighbors enrichment analysis.
    """
    from squidpy.gr import nhood_enrichment, ripley, spatial_neighbors

    if not adata.obs["cluster"].dtype.name == "category":
        adata.obs["cluster"] = adata.obs["cluster"].astype("category")

    spatial_neighbors(adata, coord_type="grid", n_neighs=4, n_rings=1)
    nhood_enrichment(adata, cluster_key=cluster_key)
    ripley(adata, cluster_key="cluster", mode="L", max_dist=500)

    return adata

