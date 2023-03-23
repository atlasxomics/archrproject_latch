''' Workflow for converting CellRanger output (fragments.tss.gz) into ArchR
objects for downstream analysis; additionally, generates UMAP and
SpatialDimPlots for a list of lsi_varfeatures.
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
    fragments_file: LatchFile
    condition: str
    positions_file: LatchFile = LatchFile(
        'latch:///position_files/all_tissue_positions_list.csv'
    )
    spatial_dir: LatchDir

class Genome(Enum):
    mm10 = 'mm10'
    hg38 = 'hg38'

@medium_task
def archr_task(
    runs: List[Run],
    project_name: str,
    genome: Genome,
    threads: int,
    tile_size: int,
    min_TSS: float,
    min_frags: int,
    lsi_iterations: int,
    lsi_resolution: float,
    lsi_varfeatures: List[int],
    clustering_resolution: int,
    umap_mindist: float
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
        f'{lsi_iterations}',
        f'{lsi_resolution}',
        f'{",".join(str(i) for i in lsi_varfeatures)}',
        f'{clustering_resolution}',
        f'{umap_mindist}',
    ]

    runs = [
        (
        f'{run.run_id},'
        f'{run.fragments_file.local_path},'
        f'{run.condition},'
        f'{run.positions_file.local_path},'
        f'{run.spatial_dir.local_path},'
    )
    for run in runs
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
        'runs': LatchParameter(
            display_name='runs',
            description='List of runs to be analyzed; each run must contain a \
                         run_id and fragments.tsv file; optional: condition, \
                         tissue position file for filtering on/off tissue, \
                         spatial folder for SpatialDimPlot.',
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
                        TileMatrix.',
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
        'lsi_iterations': LatchParameter(
            display_name='LSI iterations',
            description='iterations parameter from addIterativeLSI function.',
            batch_table_column=True,
            hidden=True
        ),                
        'lsi_resolution': LatchParameter(
            display_name='LSI resolution',
            description='resolution parameter from \
                        addIterativeLSI/clusterParams function.',
            batch_table_column=True
        ),                
        'lsi_varfeatures': LatchParameter(
            display_name='LSI varFeatures',
            description='varFeatures parameter from addIterativeLSI function; \
                        each will correspond to a umap.pdf, the last in the \
                        will be used to make the RDS object.',
            batch_table_column=True
        ),              
        'clustering_resolution': LatchParameter(
            display_name='clustering resolution',
            description='resolution parameter from addClusters function.',
            batch_table_column=True,
            hidden=True
        ),              
        'umap_mindist': LatchParameter(
            display_name='UMAP minimum distance',
            description='minDist parameter from addUMAP function.',
            batch_table_column=True,
            hidden=True
        ),                 
    },
    tags=[],
)


@workflow(metadata)
def archr_workflow(
    runs: List[Run],
    genome: Genome,
    project_name: str,
    threads: int=24,
    tile_size: int=5000,
    min_TSS: float=2.0,
    min_frags: int=0,
    lsi_iterations: int=2,
    lsi_resolution: float=0.5,
    lsi_varfeatures: List[int]=[25000],
    clustering_resolution: int=1,
    umap_mindist: float=0.0
) -> LatchDir:
    '''Pipeline for converting fragment.tsv.gz files from 10x cellranger to \
    ArchR .arrow files and ArchRProject folders.

    For data from DBiT-seq for spatially-resolved epigenomics.

    - See Deng, Y. et al 2022.
    '''

    return archr_task(
        runs=runs,
        project_name=project_name,
        genome=genome,
        threads=threads,
        tile_size=tile_size,
        min_TSS=min_TSS,
        min_frags=min_frags,
        lsi_iterations=lsi_iterations,
        lsi_resolution=lsi_resolution,
        lsi_varfeatures=lsi_varfeatures,
        clustering_resolution=clustering_resolution,
        umap_mindist=umap_mindist
    )

LaunchPlan(
    archr_workflow,
    'Test Data',
    {
        'runs' : [
            Run(
                'dev',
                LatchFile(
                    'latch:///cr_outs/ds_D01033_NG01681/outs/ds_D01033_NG01681_fragments.tsv.gz'
                ),
                'control',
                LatchFile(
                    'latch:///position_files/all_tissue_positions_list.csv'
                ),
                LatchDir(
                    'latch:////spatials/D00866/spatial'
                )
            )
        ],
        'project_name' : 'dev',
        'genome' : Genome.hg38
    },
)

if __name__ == '__main__':
    archr_workflow(
        runs=[
            Run(
                'dev',
                LatchFile(
                    'latch:///cr_outs/ds_,v_NG01681/outs/ds_D01033_NG01681_fragments.tsv.gz'
                ),
                'control',
                LatchFile(
                    'latch:///position_files/D01033/test.csv'
                ),
                LatchDir(
                    'latch:////spatials/D00866/spatial'
                )
            )
        ],
        project_name='dev',
        genome=Genome.hg38,
        lsi_varfeatures=[25000, 10000]
    )

