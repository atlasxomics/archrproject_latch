import logging

from wf.upload_to_registry import Run

# Shared AtlasXomics utilities — re-exported for `import wf.utils as utils`
from atx_common import (  # noqa: F401
    Genome,
    copy_peak_files,
    filter_anndata,
    get_groups,
    rename_obs_columns,
    sanitize_condition,
)


logging.basicConfig(
    format="%(levelname)s - %(asctime)s - %(message)s", level=logging.INFO
)


# ---------------------------------------------------------------------------
# Workflow-specific reference data
# ---------------------------------------------------------------------------
coverage_dict = {
    "cluster": "Clusters",
    "sample": "Sample",
    "condition": "condition_1",
}
