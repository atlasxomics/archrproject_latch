import anndata

from typing import List

from wf.upload_to_registry import Run


def filter_anndata(
    adata: anndata.AnnData, group: str, subgroup: List[str]
) -> anndata.AnnData:
    return adata[adata.obs[group] == subgroup]


def get_groups(runs: List[Run]):
    """Set 'groups' list for differential analysis"""

    samples = [run.run_id for run in runs]
    conditions = list({run.condition for run in runs})

    groups = ["cluster"]
    if len(samples) > 1:
        groups.append("sample")
    if len(conditions) > 1:
        groups.append("condition")

    return groups


def rename_obs_columns(adata: anndata):
    rename_map = {
        "Clusters": "cluster",
        "Condition": "condition",
        "Sample": "sample",
        "nFrags": "n_fragment"
    }

    for old_name, new_name in rename_map.items():
        if old_name in adata.obs.columns:
            adata.obs.rename(columns={old_name: new_name}, inplace=True)

    return adata
