library(ArchR)
library(harmony)
library(Seurat)

args = commandArgs(trailingOnly=TRUE)

project_name <- args[1]
genome <- args[2]
threads <- as.integer(args[3])
tile_size <- as.integer(args[4])
min_TSS <- as.numeric(args[5])
min_frags <- as.integer(args[6])
lsi_iterations <- as.integer(args[7])
lsi_resolution <- as.numeric(args[8])
for (i in strsplit(args[9], ',')) {lsi_varfeatures <- as.integer(i)}
clustering_resolution <- as.integer(args[10])
umap_mindist <- as.numeric(args[11])


runs = strsplit(args[12:length(args)], ',')
inputs = c()
for (run in runs) {inputs[run[1]] = run[2]}

out_dir <- paste0(project_name, '_ArchRProject')
addArchRGenome(genome)
addArchRThreads(threads = threads)

ArrowFiles <- createArrowFiles(
   inputFiles = inputs,
   sampleNames = names(inputs),
   minTSS = min_TSS,
   minFrags = min_frags,
   maxFrags = 1e+07,
   addTileMat = TRUE,
   addGeneScoreMat = TRUE,
   offsetPlus = 0,
   offsetMinus = 0,
   TileMatParams = list(tileSize=tile_size)
)

proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = out_dir
)

# Add an additional Conditions column
for (run in runs) {proj$Condition[proj$Sample==run[1]] <- run[3]}

# Filter on-tissue
all_ontissue = c()
for (run in runs) {
  positions <- read.csv(run[4], header = FALSE)
  positions$V1 <- paste(run[1], '#', positions$V1,'-1', sep = "")
  on_tissue <- positions$V1 [which(positions$V2 == 1)]
  all_ontissue <- c(all_ontissue, on_tissue)
}
proj <- proj[proj$cellNames %in% all_ontissue]

for (i in 1:length(lsi_varfeatures)) {

  varfeatures = lsi_varfeatures[i]
  proj <- addIterativeLSI(
    ArchRProj = proj,
    useMatrix = 'TileMatrix',
    name = 'IterativeLSI',
    iterations = lsi_iterations, 
    clusterParams = list(
      resolution = c(lsi_resolution), 
      sampleCells = 10000, 
      n.start = 10
    ), 
    varFeatures = varfeatures, 
    dimsToUse = 1:30,
    force = TRUE
  )

  if (length(runs) > 1) {
    proj <- addHarmony(
      ArchRProj = proj,
      reducedDims = 'IterativeLSI',
      name = 'Harmony',
      groupBy = 'Sample',
      force = TRUE
    )
    name = 'Harmony'
  } else {
    name = 'IterativeLSI'
  }

  proj <- addClusters(
    input = proj,
    F = name,
    method = 'Seurat',
    name = 'Clusters',
    resolution = c(clustering_resolution), 
    force = TRUE
  )

  proj <- addUMAP(
    ArchRProj = proj, 
    reducedDims = name, 
    name = 'UMAP', 
    nNeighbors = 30, 
    minDist = umap_mindist,  
    metric = 'cosine',
    force = TRUE
  )

  p1 <- plotEmbedding(
    ArchRProj = proj,
    colorBy = "cellColData",
    name = "Sample",
    embedding = "UMAP"
  )

  p2 <- plotEmbedding(
    ArchRProj = proj,
    colorBy = "cellColData",
    name = "Clusters",
    embedding = "UMAP"
  )

  ggsave(
    paste0(out_dir, '/umap_', varfeatures, '.pdf'),
    p1 + p2,
    width = 10,
    height = 10
  )
}

saveArchRProject(ArchRProj = proj)
