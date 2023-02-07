''' Short workflow for converting CellRanger output (fragments.tss.gz) into ArchR
objects for downstream analysis.
'''
import subprocess
from enum import Enum

from latch import medium_task, workflow
from latch.types import (
    LatchAuthor,
    LatchDir,
    LatchFile,
    LatchMetadata,
    LatchParameter
)


class Genome(Enum):
    mm10 = 'mm10'
    hg38 = 'hg38'

@medium_task
def archr_task(
    fragments_file: LatchFile,
    run_id: str,
    genome: Genome,
    threads: int,
    tile_size: int,
    min_TSS: float,
    min_frags: int
) -> LatchDir:
    
    _archr_cmd = [
        'Rscript',
        '/root/wf/archr_objs.R',
        fragments_file.local_path,
        run_id,
        genome.value,
        f'{threads}',
        f'{tile_size}',
        f'{min_TSS}',
        f'{min_frags}',
    ]
    
    subprocess.run(_archr_cmd)

    subprocess.run(['mkdir', f'{run_id}_archr_outs'])
    _mv_cmd = [
        'mv',
        f'/root/{run_id}.arrow',
        f'/root/{run_id}_ArchRProject',
        f'{run_id}_archr_outs'
    ]

    subprocess.run(_mv_cmd)

    return LatchDir(
        f'/root/{run_id}_archr_outs',
        f'latch:///archr_outs/{run_id}_archr_outs'
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
        'fragments_file': LatchParameter(
            display_name='fragments file',
            description='fragments.tsv.gz from CellRanger outs folder.',
            batch_table_column=True, 
        ),
        'run_id': LatchParameter(
            display_name='run id',
            description='Name of run with optional prefix, default to \
                        Dxxxxx_NGxxxxx.',
            batch_table_column=True,
            placeholder='Dxxxxx_NGxxxxx'
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
                        "TileMatrix".',
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
    fragments_file: LatchFile,
    run_id: str,
    genome: Genome,
    threads: int=24,
    tile_size: int=5000,
    min_TSS: float=2.0,
    min_frags: int=0,
) -> LatchDir:
    '''
    Pipeline for converting fragment.tsv.gz files from 10x cellranger to 
    ArchR .arrow files and ArchRProject folders.

    For data from DBiT-seq for spatially-resolved epigenomics.

    - See Deng, Y. et al 2022.
    '''

    return archr_task(
        fragments_file=fragments_file,
        run_id=run_id,
        genome=genome,
        threads=threads,
        tile_size=tile_size,
        min_TSS=min_TSS,
        min_frags=min_frags
    )

if __name__ == '__main__':
    archr_workflow(
        fragments_file=LatchFile('latch:///runs/ds_D01033_NG01681/outs/fragments.tsv.gz'),
        run_id='dev',
        genome=Genome.hg38
        )

