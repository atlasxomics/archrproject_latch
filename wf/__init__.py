''' Workflow for converting ATAC fragments (fragments.tsv.gz) into ArchR
objects for downstream analysis; additionally, generates UMAP and
SpatialDimPlots for a list of lsi_varfeatures.
'''
import glob
import subprocess

from enum import Enum
from pathlib import Path
from typing import List

from latch import custom_task, workflow
from latch.resources.launch_plan import LaunchPlan
from latch.types import (
    LatchAuthor,
    LatchDir,
    LatchFile,
    LatchMetadata,
    LatchParameter,
    LatchRule
)

from wf.upload_to_registry import upload_to_registry, Run


class Genome(Enum):
    mm10 = 'mm10'
    hg38 = 'hg38'
    rnor6 = 'rnor6'


@custom_task(cpu=62, memory=384, storage_gib=4949)
def archr_task(
    runs: List[Run],
    project_name: str,
    genome: Genome,
    tile_size: int,
    min_TSS: float,
    min_frags: int,
    lsi_iterations: int,
    lsi_resolution: float,
    lsi_varfeatures: int,
    clustering_resolution: float,
    umap_mindist: float,
    num_threads: int,
    min_cells_cluster: int,
    max_clusters: int
) -> LatchDir:

    out_dir = project_name
    subprocess.run(['mkdir', f'{out_dir}'])

    _archr_cmd = [
        'Rscript',
        '/root/wf/archr_objs.R',
        project_name,
        genome.value,
        f'{tile_size}',
        f'{min_TSS}',
        f'{min_frags}',
        f'{lsi_iterations}',
        f'{lsi_resolution}',
        f'{lsi_varfeatures}',
        f'{clustering_resolution}',
        f'{umap_mindist}',
        f'{num_threads}',
        f'{min_cells_cluster}',
        f'{max_clusters}'
    ]

    runs = [
        (
            f'{run.run_id},'
            f'{run.sample_name},'
            f'{run.fragments_file.local_path},'
            f'{run.condition},'
            f'{run.positions_file.local_path},'
            f'{run.spatial_dir.local_path},'
        )
        for run in runs
    ]

    _archr_cmd.extend(runs)
    subprocess.run(_archr_cmd)

    project_dirs = glob.glob(f'{project_name}_*')
    www = glob.glob('www')
    seurat_objs = glob.glob('*.rds')
    h5_files = glob.glob('*.h5')
    R_files = glob.glob('*.R')
    image = glob.glob('.RData')

    _mv_cmd = (
        ['mv'] +
        project_dirs +
        www +
        seurat_objs +
        h5_files +
        R_files +
        image +
        [out_dir]
    )

    subprocess.run(_mv_cmd)

    csv_tables = glob.glob('*.csv')
    volcanos = glob.glob('*.txt')

    tables_dir = Path(f'/root/{out_dir}/tables')
    tables_dir.mkdir(parents=True, exist_ok=True)

    _mv_tables_cmd = ['mv'] + csv_tables + volcanos + [str(tables_dir)]

    subprocess.run(_mv_tables_cmd)

    # Move figures into subfolder
    figures = [fig for fig in glob.glob('*.pdf') if fig != 'Rplots.pdf']
    figures_dir = Path(f'/root/{out_dir}/figures')
    figures_dir.mkdir(parents=True, exist_ok=True)

    _mv_figures_cmd = ['mv'] + figures + [str(figures_dir)]

    subprocess.run(_mv_figures_cmd)

    return LatchDir(
        f'/root/{out_dir}',
        f'latch:///ArchRProjects/{out_dir}'
    )


metadata = LatchMetadata(
    display_name='create ArchRProject',
    author=LatchAuthor(
        name='AtlasXomics, Inc.',
        email='jamesm@atlasxomics.com',
        github='github.com/atlasxomics',
    ),
    repository='https://github.com/atlasxomics/archrproject_latch',
    license='MIT',
    parameters={
        'runs': LatchParameter(
            display_name='runs',
            description='List of runs to be analyzed; each run must contain a \
                        run_id and fragments.tsv file; optional: condition, \
                        tissue position file for filtering on/off tissue, \
                        spatial folder for SpatialDimPlot. Note that multiple \
                        Coditions must be separted by -, for example: \
                        Non_responder-post_treatment-One_month.',
            batch_table_column=True,
            samplesheet=True
        ),
        'project_name': LatchParameter(
            display_name='project name',
            description='Name of output directory in archr_outs/',
            batch_table_column=True,
            rules=[
                LatchRule(
                    regex="^[^/].*",
                    message="project name cannot start with a '/'"
                )
            ]
        ),
        'genome': LatchParameter(
            display_name='genome',
            description='Reference genome to be used for geneAnnotation and \
                        genomeAnnotation',
            batch_table_column=True,
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
                        enrichment score required for a cell to pass \
                        filtering.',
            batch_table_column=True,
            hidden=True
        ),
        'min_frags': LatchParameter(
            display_name='minimum fragments',
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
            description='varFeatures parameter from addIterativeLSI function.',
            batch_table_column=True
        ),
        'clustering_resolution': LatchParameter(
            display_name='clustering resolution',
            description='resolution parameter from addClusters function.',
            batch_table_column=True
        ),
        'umap_mindist': LatchParameter(
            display_name='UMAP minimum distance',
            description='minDist parameter from addUMAP function.',
            batch_table_column=True,
            hidden=True
        ),
        'num_threads': LatchParameter(
            display_name='number of threads',
            description='This is for peak calling step with MACS2. If you get \
                        OOMKilled Error set it to 1!',
            batch_table_column=True,
            hidden=True
        ),
        'min_cells_cluster': LatchParameter(
            display_name='minimum cells per cluster',
            description='Minimum number of cells in a cluster; passed to the \
                        nOutlier parameter of ArchR::addClusters. If a \
                        cluster falls below the minimum, it is merged into a \
                        neighboring cluster',
            batch_table_column=True,
            hidden=True
        ),
        'max_clusters': LatchParameter(
            display_name='maximum clusters',
            description='Maximum number of clusters allow for the project; \
                        passed to the maxClusters parameter of \
                        ArchR::addClusters. If above maximum clusters, \
                        clusters are merged.',
            batch_table_column=True,
            hidden=True
        ),
        'run_table_id': LatchParameter(
            display_name='Registry Table ID',
            description='The runs will be updated in Registry with its \
                        corresponding condition, spatial directory, \
                        condition, and location of the output archR project.'
        ),
        'project_table_id': LatchParameter(
            display_name='The ID of the SOWs Registry table',
            description='The ArchR project will be inserted into the SOW \
                        table for the corresponding runs.'
        )
    },
    tags=[],
)


@workflow(metadata)
def archrproject_workflow(
    runs: List[Run],
    genome: Genome,
    project_name: str,
    tile_size: int = 5000,
    min_TSS: float = 2.0,
    min_frags: int = 0,
    lsi_iterations: int = 2,
    lsi_resolution: float = 0.5,
    lsi_varfeatures: int = 25000,
    clustering_resolution: float = 1.0,
    umap_mindist: float = 0.0,
    num_threads: int = 50,
    min_cells_cluster: int = 20,
    max_clusters: int = 25,
    run_table_id: str = "761",
    project_table_id: str = "779"
) -> LatchDir:
    '''Workflow for converting fragment.tsv.gz files to ArchRProjects.

    # create ArchRProject

    **create ArchRProject** is a [latch.bio](https://latch.bio/) workflow for
    generating R objects and data for downstream analysis of epigenomic
    [DBiT-seq](https://www.nature.com/articles/s41586-022-05094-1) experiments.
    Provided fragments from a single-cell ATAC-seq preprocessing and alignment
    workflow and spatial information, **create ArchRProjectr** performs the
    heavy computational steps of ArchR and Seurat and returns files that can
    be easily input into custom scripts for more neuanced analysis without the
    need to perform intensive computation.  The workflow utilizes
    [ArchR](https://www.archrproject.com/articles/Articles/tutorial.html) to
    perform epigenomic single-cell analysis and
    [Seurat](https://satijalab.org/seurat/) to spatially align the data.
    The workflow can take data from either a single tissue-sample analyzed via
    DBiT-seq or multiple tissue-samples; in ATX parlance, tissue-samples
    analyzed via DBIT-seq are termed 'Runs'.  All Runs input to **create
    ArchRProject** are merged into a single ArchRProject for analysis.

    ## Inputs
    All input files for **create ArchRProject** must be on the latch.bio
    [file system](https://wiki.latch.bio/wiki/data/overview).  Each run in the
    workflow takes the following parameters,
    * [fragments.tsv.gz file](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/fragments): A BED-like, tab-delimited file in which each row contains an ATAC-seq fragment
    * [tissue_positions_list.csv](https://docs.atlasxomics.com/projects/AtlasXbrowser/en/latest/SpatialFolder.html): A comma-separated file in which each row contains a unique barcode, an indicator for whether the tixel is 'on-tissue' (1, 0), and a row/column index
    * [Spatial folder](https://docs.atlasxomics.com/projects/AtlasXbrowser/en/latest/SpatialFolder.html): A directory containing tissue images and experiment metadata
    * Run ID: An identifier for the run
    * Condition (_optional_):  An experimental Condition descriptor
    (ie. 'control', 'diseased')
    Individual runs are batched in a Project with the following global
    parameters,
    * Project Name: A name for the output folder
    * Genome: A reference genome to be used for alignment
    * Upload to SLIMS _(ATX-internal only)_: A T/F toggle for whether to push
    QC results to LIMS
    * LSI resolution: A decimal value used as input to the `clusterParams`
    parameter of the `addIterativeLSI` function in
    [ArchR](https://www.archrproject.com/reference/addIterativeLSI.html);
    * LSI varFeatures: An integer used as input to the `varFeatures` parameter
    of the `addIterativeLSI` function in
    [ArchR](https://www.archrproject.com/reference/addIterativeLSI.html);
    * clustering resolution: A decimal value used as input to the `resolution`
    parameter of the `addClusters` function in
    [ArchR](https://www.archrproject.com/reference/addClusters.html).
    > The Project also takes a series of single-value parameters that can be
    found under the 'Hidden Parameters' dropdown; these parameters are less
    commonly varied inputs to ArchR functions.

    ## Running the workflow (_in progress_)
    The **create ArchRProject** workflow can be found in the
    [Workflows](https://wiki.latch.bio/workflows/overview) module in your
    latch.bio workspace. For access to an ATX-collaborator workspace, please
    contact your AtlasXomics Support Scientist or email
    support@atlasxomics.com.  See
    [here](https://wiki.latch.bio/workflows/overview) for general instructions
    for running workflows in latch.bio.
    1. Navigate to the **optimize archr** workflow in the Workflows module in
    your latch.bio workspace.  Ensure you are on the 'Parameters' tab of the
    workflow.
    2. To add Runs to the Project, select the '+ runs' icon.  Add values for
    the Run parameters described above; repeat for each Run in the Project.
    3. Scroll to the bottom of the page and input values for global project
    parameters.
    4. Click the 'Hidden Parameters' button and change the global parameters
    as needed.
    5. Click the 'Launch Workflow' button on the bottom-right of the parameters
      page.  This will automatically navigate you to the Executions tab of the
      workflow.
    6. From the Executions tab, you can view the status of the launched
    workflow. Once the workflow has completed running, the status will change
    to 'Succeeded'; if the workflow has the status 'Failed', please contact an
    AtlasXomics Support Scientist.  You can click on the workflow execution to
    view a more granular workflow status and see output logs.
    7. Workflow outputs are loaded into the latch.bio
    [Data module](https://wiki.latch.bio/wiki/data/overview) in the
    `ArchRProjects` directory.

    ## Outputs (_in progress_)
    Outputs from **create ArchRProject** are loaded into latch.bio
    [Data module](https://wiki.latch.bio/wiki/data/overview) in the
    `ArchRProjects` directory.
    * ArchRProject/
    * SeuratObj.rds
    * SeuratObjMotif.rds
    * UMAPHarmony.csv
    * enrichMotifs_clusters.rds
    * enrichMotifs_sample.rds
    * enrichMotifs_condition.rds
    * genes_per_cluster_hm.csv
    * genes_per_sample_hm.csv
    * genes_per_condition_hm.csv
    * inpMarkers.txt
    * inpMarkers_motif.txt
    * markersGS_clusters.rds
    * markersGS_sample.rds
    * markersGS_condition.rds
    * motif_per_cluster_hm.csv
    * motif_per_sample_hm.csv
    * motif_per_condition_hm.csv
    * req_genes1.csv
    * req_genes2.csv
    * req_genes3.csv
    * req_motifs1.csv
    * req_motifs2.csv
    * req_motifs3.csv
    * seqlogo.rds
    ## Next Steps
    Analysis can be performed locally or in a latch.bio
    [Pod](https://wiki.latch.bio/wiki/pods/overview).  For access to
    ATX-specific Pods, please contact your AtlasXomics Support Scientist.
    Output from **create ArchRProject** can be provided to the **atlasShiny**
    workflow to create input for a DBiT-seq-specific R Shiny App.  For access
    to this workflow and app, please contact your AtlasXomics Support
    Scientist.

    ## Support
    Questions? Comments?  Contact support@atlasxomics.com or post in
    AtlasXomics
    [Discord](https://discord.com/channels/1004748539827597413/1005222888384770108).

    '''
    archr_project = archr_task(
        runs=runs,
        project_name=project_name,
        genome=genome,
        tile_size=tile_size,
        min_TSS=min_TSS,
        min_frags=min_frags,
        lsi_iterations=lsi_iterations,
        lsi_resolution=lsi_resolution,
        lsi_varfeatures=lsi_varfeatures,
        clustering_resolution=clustering_resolution,
        umap_mindist=umap_mindist,
        num_threads=num_threads,
        min_cells_cluster=min_cells_cluster,
        max_clusters=max_clusters
    )

    upload_to_registry(
        runs=runs,
        archr_project=archr_project,
        run_table_id=run_table_id,
        project_table_id=project_table_id
    )

    return archr_project


LaunchPlan(
    archrproject_workflow,
    'defaults',
    {
        'runs': [
            Run(
                'default',
                'default',
                LatchFile(
                    'latch:///chromap_outputs/demo/chromap_output/fragments.tsv.gz'
                ),
                'demo',
                LatchDir('latch:///spatials/demo/spatial'),
                LatchFile(
                    'latch:///spatials/demo/spatial/tissue_positions_list.csv'
                ),
                )
        ],
        'project_name': 'demo',
        'genome': Genome.hg38,
        'run_table_id': '761',
        'project_table_id': '779'
    },
)
