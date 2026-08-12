import logging

from latch import small_task
from latch.types import LatchDir, LatchFile
from latch.registry.table import Table
from dataclasses import dataclass
from typing import List, Optional
from latch.functions.messages import message

from atx_common import get_LatchFile  # noqa: F401 — re-exported

logging.basicConfig(format="%(levelname)s - %(asctime)s - %(message)s")


@dataclass
class Run:
    run_id: str
    sample_name: str
    fragments_file: LatchFile
    spatial_dir: LatchDir
    condition: str = "None"


@small_task(retries=0)
def upload_to_registry(
    runs: List[Run],
    archr_project: LatchDir,
    run_table_id: str = "761",
    project_table_id: str = "917",
) -> LatchDir:
    run_table = Table(run_table_id)
    project_table = Table(project_table_id)
    try:
        with run_table.update() as updater:
            for run in runs:
                message(
                    "info",
                    {
                        "title": f"Updating run {run.run_id} in registry table ID {run_table_id}",
                        "body": f"Run {run.run_id}, condition {run.condition}",
                    },
                )

                updater.upsert_record(
                    run.run_id,
                    condition=run.condition,
                    spatial_directory=run.spatial_dir,
                    archrproject_outs=archr_project
                )

        for run in runs:  # loop through projects with linked uns
            for page in project_table.list_records():
                for project_id, record in page.items():
                    project = record.get_values()
                    project_name = record.get_name()
                    print(project)
                    try:
                        if len(project['Runs']) > 0:
                            for project_run in project['Runs']:
                                record_name = project_run.get_name()
                                if record_name == run:
                                    print(record_name)
                                    with project_table.update() as updater:
                                        message(
                                            "info",
                                            {
                                                "title": f"Updating project {project_name} in registry table ID 917"
                                            },
                                        )
                                        updater.upsert_record(
                                            project_name,
                                            test_string="DONE"
                                        )

                                        break
                    except:
                        break

    except Exception as err:
        print(f"Unexpected {err=}, {type(err)=}")
    finally:
        return archr_project

