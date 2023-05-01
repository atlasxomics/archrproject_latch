library(ArchR)
library(ggplot2)
library(harmony)
library(Seurat)

# globals ---------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

project_name <- args[1]
genome <- args[2]
threads <- as.integer(args[3])
tile_size <- as.integer(args[4])
min_tss <- as.numeric(args[5])
min_frags <- as.integer(args[6])
lsi_iterations <- as.integer(args[7])
lsi_resolution <- as.numeric(args[8])
for (i in strsplit(args[9], ",")) {
  lsi_varfeatures <- as.integer(i)
  }
clustering_resolution <- as.integer(args[10])
umap_mindist <- as.numeric(args[11])

runs <- strsplit(args[12:length(args)], ",")
inputs <- c()
for (run in runs) {
  inputs[run[1]] <- run[2]
  }

out_dir <- paste0(project_name, "_ArchRProject")

build_atlas_seurat_object <- function(
  run_id,
  matrix,
  metadata,
  spatial_path) {
  # Prepare and combine gene matrix, metadata, and image for seurat object
  # for runs within a project.

  image <- Read10X_Image(
    image.dir = spatial_path,
    filter.matrix = TRUE
  )
  metadata <- metadata[metadata$Sample == run_id, ]

  matrix <- matrix[, c(grep(pattern = run_id, colnames(matrix)))]
  matrix@Dimnames[[2]] <- metadata@rownames

  object <- CreateSeuratObject(
    counts = matrix,
    assay  = "Spatial",
    meta.data = as.data.frame(metadata)
  )
  image <- image[Cells(x = object)]
  DefaultAssay(object = image) <- "Spatial"
  object[["slice1"]] <- image
  return(object)
}

spatial_plot <- function(seurat_object, name) {
  clusters <- sort(unique(seurat_object$Clusters))
  colors <- ArchRPalettes$stallion2[seq_len(length(clusters))]
  names(colors) <- clusters
  SpatialDimPlot(
    seurat_object,
    group.by = "Clusters",
    pt.size.factor = 1,
    cols = colors,
    stroke = 0) +
  ggtitle(name) +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 21))
}

feature_plot <- function(seurat_obj, feature, name) {
  Seurat::Spatialfeature_plot(
    object = seurat_obj,
    features = feature,
    alpha = c(0.2, 1),
    pt.size.factor = 1
    ) +
  ggtitle(paste0(feature, " : ", name)) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 15)
  )
}

# create archr project --------------------------------------------------------

addArchRGenome(genome)
addArchRThreads(threads = threads)

arrow_files <- createArrowFiles(
   inputFiles = inputs,
   sampleNames = names(inputs),
   minTSS = min_tss,
   minFrags = min_frags,
   maxFrags = 1e+07,
   addTileMat = TRUE,
   addGeneScoreMat = TRUE,
   offsetPlus = 0,
   offsetMinus = 0,
   TileMatParams = list(tileSize = tile_size)
)

proj <- ArchRProject(
  ArrowFiles = arrow_files,
  outputDirectory = out_dir
)

# Add an additional Conditions column
for (run in runs) {
  proj$Condition[proj$Sample == run[1]] <- run[3]
}

# Filter on-tissue
all_ontissue <- c()
for (run in runs) {
  positions <- read.csv(run[4], header = FALSE)
  positions$V1 <- paste(run[1], "#", positions$V1, "-1", sep = "")
  on_tissue <- positions$V1 [which(positions$V2 == 1)]
  all_ontissue <- c(all_ontissue, on_tissue)
}
proj <- proj[proj$cellNames %in% all_ontissue]

# iterate plotting ------------------------------------------------------------

for (i in seq_along((lsi_varfeatures))) {

  varfeatures <- lsi_varfeatures[i]
  proj <- addIterativeLSI(
    ArchRProj = proj,
    useMatrix = "TileMatrix",
    name = "IterativeLSI",
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
      reducedDims = "IterativeLSI",
      name = "Harmony",
      groupBy = "Sample",
      force = TRUE
    )
    name <- "Harmony"
  } else {
    name <- "IterativeLSI"
  }
  proj <- addClusters(
    input = proj,
    F = name,
    method = "Seurat",
    name = "Clusters",
    resolution = c(clustering_resolution),
    force = TRUE
  )
  proj <- addUMAP(
    ArchRProj = proj,
    reducedDims = name,
    name = "UMAP",
    nNeighbors = 30,
    minDist = umap_mindist,
    metric = "cosine",
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
    paste0(out_dir, "/umap_", varfeatures, ".pdf"),
    p1 + p2,
    width = 10,
    height = 10
  )

  proj <- addImputeWeights(proj)

  # Create metadata object for Seurat object
  metadata <- getCellColData(ArchRProj = proj)
  rownames(metadata) <- str_split_fixed(
  str_split_fixed(
    row.names(metadata),
    "#",
    2)[, 2],
  "-",
  2)[, 1]
  metadata["log10_nFrags"] <- log(metadata$nFrags)

  # Create gene matrix for Seurat object.
  gene_matrix <- getMatrixFromProject(
    ArchRProj = proj,
    useMatrix = "GeneScoreMatrix"
  )
  matrix <- imputeMatrix(
    mat = assay(gene_matrix),
    imputeWeights = getImputeWeights(proj)
  )
  gene_row_names <- gene_matrix@elementMetadata$name
  rownames(matrix) <- gene_row_names

  seurat_objs <- c()
  for (run in runs) {

    obj <- build_atlas_seurat_object(
      run_id = run[1],
      matrix = matrix,
      metadata = metadata,
      spatial_path = run[5]
    )

    saveRDS(
      obj,
      file = paste0(out_dir, "/", run[1], "_SeuratObj_", varfeatures, ".rds")
      )

    p1 <- spatial_plot(
      obj,
      name = paste(run[1], varfeatures)
      )

    ggsave(
      paste0(out_dir, "/", run[1], "_spatialdim_", varfeatures, ".pdf"),
      p1,
      width = 10,
      height = 10
    )
    seurat_objs <- c(seurat_objs, obj)
    }
}

for (obj in seurat_objs) {
  name <- unique(obj@meta.data[["Sample"]])

  nfrags_plot <- feature_plot(obj, "log10_nFrags", name)
  tss_plot <- feature_plot(obj, "TSSEnrichment", name)

  pdf(paste0(out_dir, "/", name, "_qc_plots.pdf"))
    print(nfrags_plot)
    par(newpage = TRUE)
    print(tss_plot)
  dev.off()

}

saveArchRProject(ArchRProj = proj)
