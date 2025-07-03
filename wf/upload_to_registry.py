import logging

from pathlib import Path

from latch import small_task
from latch.types import LatchDir, LatchFile
from latch.registry.table import Table
from dataclasses import dataclass
from typing import List
from latch.functions.messages import message

logging.basicConfig(format="%(levelname)s - %(asctime)s - %(message)s")


@dataclass
class Run:
    run_id: str
    sample_name: str
    fragments_file: LatchFile
    spatial_dir: LatchDir
    condition: str = "None"


def get_LatchFile(directory: LatchDir, file_name: str) -> LatchFile:
    try:
        files = [file for file in directory.iterdir()
                 if isinstance(file, LatchFile) and
                 Path(file.path).name == file_name]
        if len(files) == 1:
            return files[0]
        elif len(files) == 0:
            raise FileNotFoundError(
                f"No file {file_name} found in {directory.remote_path}"
            )
        elif len(files) > 1:
            raise FileNotFoundError(
                f"Multiple files {file_name} found in {directory.remote_path}"
            )
    except Exception as e:
        logging.error(f"Failed to find file '{file_name}'; error {e}")
        return None


@small_task(retries=0)
def upload_to_registry(
    runs: List[Run],
    archr_project: LatchDir,
    run_table_id: str = "761",
    project_table_id: str = "917",
):
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
        return
    finally:
        return


if __name__ == "__main__":
    upload_to_registry(
        runs=[
            Run(
                'D1291',
                LatchFile('latch:///atac_outs/6bp_D01291_NG02620/outs/6bp_D01291_NG02620_fragments.tsv.gz'),
                'old',
                LatchDir('latch:///atx-illumina-1682977469.0200825/Images_spatial/D1291/spatial'),
            ),
            Run(
                'D1292',
                LatchFile('latch:///cleaned/Babayev_cleaned/cleaned_D01292_NG02621_fragments.tsv.gz'),
                'old',
                LatchDir('latch:///atx-illumina-1682977469.0200825/Images_spatial/D1292/spatial'),
            ),
            Run(
                'D1293',
                LatchFile('latch:///atac_outs/6bp_D01293_NG02622/outs/6bp_D01293_NG02622_fragments.tsv.gz'),
                'young',
                LatchDir('latch:///atx-illumina-1682977469.0200825/Images_spatial/D1293/spatial'),
            ),
            Run(
                'D1294',
                LatchFile('latch:///atac_outs/D01294_NG02624/outs/D01294_NG02624_fragments.tsv.gz'),
                'young',
                LatchDir('latch:///atx-illumina-1682977469.0200825/Images_spatial/D1294/spatial'),
            )
        ],
        archr_project=LatchDir("latch://13502.account/ArchRProjects/Babeyev"),
        run_table_id="761",
        project_table_id="917",
    )
