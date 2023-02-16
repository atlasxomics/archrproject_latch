''' Short workflow for converting CellRanger output (fragments.tss.gz) into ArchR
objects for downstream analysis.
'''
import subprocess

from dataclasses import dataclass
from dataclasses_json import dataclass_json
from enum import Enum
from typing import List

from latch import medium_task, workflow
from latch.resources.launch_plan import LaunchPlan
from latch.types import (
    LatchAuthor,
    LatchDir,
    LatchFile,
    LatchMetadata,
    LatchParameter
)

@dataclass_json
@dataclass
class Run:
    run_id: str
    condition: str
    fragments_file: LatchFile

class Genome(Enum):
    mm10 = 'mm10'
    hg38 = 'hg38'

@medium_task
def archr_task(
    fragments_files: List[Run],
    project_name: str,
    genome: Genome,
    threads: int,
    tile_size: int,
    min_TSS: float,
    min_frags: int
) -> LatchDir:
    
    _archr_cmd = [
        'Rscript',
        '/root/wf/archr_objs.R',
        project_name,
        genome.value,
        f'{threads}',
        f'{tile_size}',
        f'{min_TSS}',
        f'{min_frags}',
    ]

    runs = [
    f'{run.run_id},{run.fragments_file.local_path},{run.condition}'
    for run in fragments_files
    ]
    
    _archr_cmd.extend(runs)
    subprocess.run(_archr_cmd)

    out_dir = f'{project_name}_ArchRProject'
    subprocess.run(
            [
        'mv',
        f'{out_dir}/Save-ArchR-Project.rds',
        f'{out_dir}/{project_name}.rds'
        ]
    )

    return LatchDir(
        f'/root/{out_dir}',
        f'latch:///archr_outs/{out_dir}'
    )


metadata = LatchMetadata(
    display_name='archr',
    author=LatchAuthor(
        name='James McGann',
        email='jpaulmcgann@gmail.com',
        github='github.com/jpmcga',
    ),
    repository='https://github.com/jpmcga/archr_latch/',
    license='MIT',
    parameters={
        'fragments_files': LatchParameter(
            display_name='fragment files',
            description='List mapping Run ID and condition to file path to fragments.tsv file.',
            batch_table_column=True, 
        ),
        'project_name' : LatchParameter(
            display_name='project name',
            description='Name prefix of output ArchRProject folder.',
            batch_table_column=True,
        ),
        'genome': LatchParameter(
            display_name='genome',
            description='Reference genome to be used for geneAnnotation and \
                        genomeAnnotation',
            batch_table_column=True,
        ),
        'threads': LatchParameter( # Might want to set a rule here
            display_name='threads',
            description='The number of threads to be used for parallel \
                        computing; max 24',
            batch_table_column=True,
            hidden=True
        ),
        'tile_size': LatchParameter(
            display_name='tile size', 
            description='The size of the tiles used for binning counts in the \
                        'TileMatrix'.',
            batch_table_column=True,
            hidden=True
        ),
        'min_TSS': LatchParameter(
            display_name='minimum TSS',
            description='The minimum numeric transcription start site (TSS) \
                    enrichment score required for a cell to pass filtering.',
            batch_table_column=True,
            hidden=True
        ),
        'min_frags': LatchParameter(
            display_name='minumum fragments',
            description='The minimum number of mapped ATAC-seq fragments \
                        required per cell to pass filtering.',
            batch_table_column=True,
            hidden=True
        ),        
    },
    tags=[],
)


@workflow(metadata)
def archr_workflow(
    fragments_files: List[Run],
    genome: Genome,
    project_name: str,
    threads: int=24,
    tile_size: int=5000,
    min_TSS: float=2.0,
    min_frags: int=0,
) -> LatchDir:
    '''Pipeline for converting fragment.tsv.gz files from 10x cellranger to \
    ArchR .arrow files and ArchRProject folders.

    For data from DBiT-seq for spatially-resolved epigenomics.

    - See Deng, Y. et al 2022.
    '''

    return archr_task(
        fragments_files=fragments_files,
        project_name=project_name,
        genome=genome,
        threads=threads,
        tile_size=tile_size,
        min_TSS=min_TSS,
        min_frags=min_frags,
    )

LaunchPlan(
    archr_workflow,
    'Test Data',
    {
        'fragments_files' : [
                Run(
                    'dev',
                    'control',
                    LatchFile('/cr_outs/ds_D01033_NG01681/outs/ds_D01033_NG01681_fragments.tsv.gz')
                )
            ],
        'project_name' : 'dev',
        'genome' : Genome.hg38
    },
)

if __name__ == '__main__':
    archr_workflow(
        fragments_files=[
            Run(
                'dev',
                'control',
                LatchFile('/cr_outs/ds_D01033_NG01681/outs/ds_D01033_NG01681_fragments.tsv.gz')
            )
        ],
        project_name='dev',
        genome=Genome.hg38
    )

