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

runs = strsplit(args[7:length(args)], ',')
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

proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = 'TileMatrix',
  name = 'IterativeLSI',
  iterations = 2,
  clusterParams = list(
    resolution = c(0.5), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 50000, 
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
}else{
  name = 'IterativeLSI'
}

proj <- addClusters(
  input = proj,
  F = name,
  method = 'Seurat',
  name = 'Clusters',
  resolution = c(1),
  force = TRUE
)

proj <- addUMAP(
  ArchRProj = proj, 
  reducedDims = name, 
  name = 'UMAP', 
  nNeighbors = 30, 
  minDist = 0.0, 
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

ggsave(paste0(out_dir, '/umap.pdf'), p1+p2, width = 10, height = 10)

saveArchRProject(ArchRProj = proj)
