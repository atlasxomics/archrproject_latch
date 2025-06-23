library("ArchR")
library("BSgenome")
library("BSgenome.Hsapiens.UCSC.hg38")
library("BSgenome.Mmusculus.UCSC.mm10")
library("BSgenome.Rnorvegicus.UCSC.rn6")
library("ComplexHeatmap")
library("chromVARmotifs")
library("circlize")
library("data.table")
library("dplyr")
library("GenomicRanges")
library("gridExtra")
library("harmony")
library("plyr")
library("qdap")
library("readr")
library("Seurat")
library("seqLogo")
library("ShinyCell")
library("tidyverse")

source("/root/getDeviation_ArchR.R")
source("/root/makeShinyFiles.R")
source("/root/wf/convert.R")
source("/root/wf/utils.R")


# globals ---------------------------------------------------------------------
set.seed(42)

args <- commandArgs(trailingOnly = TRUE)
print(args)

project_name <- args[1]
genome <- args[2]
tile_size <- as.integer(args[3])
min_tss <- as.numeric(args[4])
min_frags <- as.integer(args[5])
lsi_iterations <- as.integer(args[6])
lsi_resolution <- as.numeric(args[7])
lsi_varfeatures <- as.integer(args[8])
clustering_resolution <- as.numeric(args[9])
umap_mindist <- as.numeric(args[10])
num_threads <- as.integer(args[11])
min_cells_cluster <- as.integer(args[12])
max_clusters <- as.integer(args[13])
print(paste("Number of threads:", num_threads))

runs <- strsplit(args[14:length(args)], ",")
runs

inputs <- c()
for (run in runs) {
  inputs[run[1]] <- run[3]
}
inputs

out_dir <- paste0(project_name, "_ArchRProject")

# save input metrics in csv
metrics <- as.list(args[1:13])
names(metrics) <- c(
  "project_name",
  "genome",
  "tile_size",
  "minimum_tss",
  "minimum_fragments",
  "lsi_iterations",
  "lsi_resolution",
  "lsi_varFeatures",
  "clustering_resolution",
  "umap_minimum_distance",
  "number_threads",
  "min_cells_cluster",
  "max_clusters"
)
write.csv(metrics, file = "metadata.csv", row.names = FALSE)


# create archr project --------------------------------------------------------

addArchRThreads(threads = 50)

if (genome != "rnor6") {

  addArchRGenome(genome)
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
} else if (genome == "rnor6") {

  load(
    "/root/custom_ArchR_genomes_and_annotations/rn6/rn6_liftoff_mm10NcbiRefSeq_ArchR_annotations.rda"
  )
  arrow_files <- createArrowFiles(
    inputFiles = inputs,
    geneAnnotation = geneAnnotation,
    genomeAnnotation = genomeAnnotation,
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
    outputDirectory = out_dir,
    geneAnnotation = geneAnnotation,
    genomeAnnotation = genomeAnnotation
  )
}

# Add an additional Conditions column
for (run in runs) {
  proj$Condition[proj$Sample == run[1]] <- run[4]
  proj$SampleName[proj$Sample == run[1]] <- run[2]
}

# Filter on-tissue
all_ontissue <- c()
for (run in runs) {
  positions <- read.csv(run[5], header = FALSE)
  positions$V1 <- paste0(run[1], "#", positions$V1, "-1")
  on_tissue <- positions$V1 [which(positions$V2 == 1)]
  all_ontissue <- c(all_ontissue, on_tissue)
}

proj <- proj[proj$cellNames %in% all_ontissue]

saveArchRProject(
  ArchRProj = proj, outputDirectory = paste0(project_name, "_ArchRProject")
)

# dimension reduction and clustering ------------------------------------------

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
  varFeatures = lsi_varfeatures,
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

proj <- ArchR::addClusters(
  input = proj,
  reducedDims = name,
  method = "Seurat",
  name = "Clusters",
  resolution = c(clustering_resolution),
  nOutlier = min_cells_cluster,
  maxClusters = max_clusters,
  force = TRUE
)
print(table(proj$Clusters))

proj <- ArchR::addUMAP(
  ArchRProj = proj,
  reducedDims = name,
  name = "UMAP",
  nNeighbors = 30,
  minDist = umap_mindist,
  metric = "cosine",
  force = TRUE
)

proj <- addImputeWeights(proj)

saveArchRProject(
  ArchRProj = proj, outputDirectory = paste0(project_name, "_ArchRProject")
)

conds <- strsplit(proj$Condition, split = "\\s|-")

fun1 <- function(lst, n) {
  sapply(lst, `[`, n)
}

for (i in seq_along(conds[[1]])) {
  proj <- addCellColData(
    proj,
    data = fun1(conds, i),
    name = paste0("condition_", i),
    cells = proj$cellNames,
    force = TRUE
  )
}
treatment <- names(getCellColData(proj))[
  grep("condition_", names(getCellColData(proj)))
]

print("++++ What treatments are in this project? ++++")
treatment

umap_plots <- c()
for (name in c("Clusters", "Sample")) {
  umap <- plot_umap(proj, name) # from utils
  umap_plots[[name]] <- umap
}

for (i in seq_along(treatment)) {
  umap <- plot_umap(proj, treatment[i]) # from utils
  umap_plots[[treatment[i]]] <- umap
}

pdf("umap_plots.pdf")
grid.arrange(
  grobs = umap_plots,
  ncol = 2
)
dev.off()

# create seurat objects -------------------------------------------------------

# create metadata object for Seurat object
metadata <- getCellColData(ArchRProj = proj)
rownames(metadata) <- str_split_fixed(
  str_split_fixed(
    row.names(metadata),
    "#",
    2)[, 2],
  "-",
  2
)[, 1]
metadata["log10_nFrags"] <- log(metadata$nFrags)

# create gene matrix for Seurat object
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

print("++++ creating seurat objs ++++")

seurat_objs <- c()
for (run in runs) {

  obj <- build_atlas_seurat_object(
    run_id = run[1],
    matrix = matrix,
    metadata = metadata,
    spatial_path = run[6]
  )

  saveRDS(obj, file = paste0(run[1], "_SeuratObj.rds"))
  seurat_objs <- c(seurat_objs, obj)
}

print("++++ creating spatial plots ++++")

spatial_cluster_plots <- list()
for (i in seq_along(seurat_objs)) {
  plot <- plot_spatial(seurat_objs[[i]], unique(seurat_objs[[i]]$Sample))
  spatial_cluster_plots[[i]] <- plot
}

spatial_lists <- split(
  spatial_cluster_plots, ceiling(seq_along(spatial_cluster_plots) / 4)
)

pdf("spatial_plots.pdf")
for (i in seq_along(spatial_lists)) {
  print(CombinePlots(spatial_lists[[i]], legend = "bottom"))
}
dev.off()

print("++++ creating qc plots ++++")

# feature plot from utils
metrics <- c("TSSEnrichment", "nFrags", "log10_nFrags")
all_qc_plots <- list()
for (i in seq_along(metrics)) {
  spatial_qc_plots <- list()
  for (j in seq_along(seurat_objs)) {
    plot <- plot_feature(
      seurat_objs[[j]], metrics[i], unique(seurat_objs[[j]]$Sample)
    )
    spatial_qc_plots[[j]] <- plot
  }
  all_qc_plots[[i]] <- spatial_qc_plots
}

print("These are the available SeuratObjects:")
seurat_objs

all <-  list()
for (i in seq_along(seurat_objs)) {
  all[[i]] <- seurat_objs[[i]]
  all[[i]] <- RenameCells(
    all[[i]],
    new.names = paste0(
      unique(all[[i]]@meta.data$Sample),
      "#",
      colnames(all[[i]]),
      "-1"
    )
  )
}

pdf("qc_plots.pdf")
for (i in seq_along(metrics)) {
  lists <- split(all_qc_plots[[i]], ceiling(seq_along(all_qc_plots[[i]]) / 6))
  for (list in lists) {
    grid.arrange(
      grobs = list,
      ncol = 2
    )
  }
}
dev.off()

# save ArchR object
saveArchRProject(
  ArchRProj = proj, outputDirectory = paste0(project_name, "_ArchRProject")
)

############-------------Identifying Marker Genes------------##################
# per cluster

markersGS <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "ttest"
)

saveRDS(markersGS, "markersGS_clusters.rds")

marker_list <- getMarkers(markersGS, cutOff = "FDR <= 0.02")
write.csv(
  marker_list,
  file = "marker_genes_per_cluster.csv",
  row.names = FALSE
)

for (rd in names(proj@reducedDims)){
  if (rd == "Harmony") {
    proj <- addImputeWeights(proj, reducedDims = "Harmony")
  } else {
    proj <- addImputeWeights(proj)
  }
}
# save for shiny

heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS,
  cutOff = "Pval <= 0.05 & Log2FC >= 0.10",
  plotLog2FC = TRUE,
  transpose = FALSE,
  returnMatrix = TRUE
)

write.csv(heatmapGS, "genes_per_cluster_hm.csv")

######################### Plot marker genes ###################################

heatmaps <- list()

# recompute for heatmap plot
gene_cutoff <- "Pval <= 0.05 & Log2FC >= 0.1"
heatmap_gs_plotting <- plotMarkerHeatmap(
  seMarker = markersGS,
  cutOff = gene_cutoff,
  transpose = TRUE
)
# save for plotting with peaks and motifs
gene_hm <- ComplexHeatmap::draw(
  heatmap_gs_plotting,
  heatmap_legend_side = "bot",
  annotation_legend_side = "bot",
  column_title = paste0("Marker genes (", gene_cutoff, ")"),
  column_title_gp = gpar(fontsize = 12)
)
heatmaps[[1]] <- gene_hm

##############################################################################

# per sample

######-----------Identifying Marker Genes-------------#######

if (length(unique(proj$Sample)) > 1) {

  markersGS <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "GeneScoreMatrix",
    groupBy = "Sample",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "ttest",
  )

  # save for shiny app
  saveRDS(markersGS, "markersGS_sample.rds")

  marker_list <- getMarkers(markersGS, cutOff = "FDR <= 0.02")
  write.csv(
    marker_list,
    file = "marker_genes_per_sample.csv",
    row.names = FALSE
  )
  for (rd in names(proj@reducedDims)) {
    if (rd == "Harmony") {
      proj <- addImputeWeights(proj, reducedDims = "Harmony")
    } else {
      proj <- addImputeWeights(proj)
    }
  }
  # save for shiny
  heatmapGS <- plotMarkerHeatmap(
    seMarker = markersGS,
    cutOff = "Pval <= 0.05 & Log2FC >= 0.10",
    plotLog2FC = TRUE,
    transpose = FALSE,
    returnMatrix = TRUE
  )
  write.csv(heatmapGS, "genes_per_sample_hm.csv")

} else {
  heatmapGS <- "There are not enough samples to be compared with!"
}

# per treatment
######---------------------Identifying Marker Genes----------------------######
if (length(unique(proj$Condition)) > 1) {

  for (i in seq_along(treatment)) {

    markersGS <- getMarkerFeatures(
      ArchRProj = proj,
      useMatrix = "GeneScoreMatrix",
      groupBy = treatment[i],
      bias = c("TSSEnrichment", "log10(nFrags)"),
      testMethod = "ttest",
    )
    # save for shiny app
    saveRDS(markersGS, paste0("markersGS_condition_", i, ".rds"))

    marker_list <- getMarkers(markersGS, cutOff = "FDR <= 0.02")
    write.csv(
      marker_list,
      file = paste0("marker_genes_per_condition_", i, ".csv"),
      row.names = FALSE
    )

    for (rd in names(proj@reducedDims)) {
      if (rd == "Harmony") {
        proj <- addImputeWeights(proj, reducedDims = "Harmony")
      } else {
        proj <- addImputeWeights(proj)
      }
    }

    # save for shiny
    heatmapGS <- plotMarkerHeatmap(
      seMarker = markersGS,
      cutOff = "Pval <= 0.05 & Log2FC >= 0.10",
      plotLog2FC = TRUE,
      transpose = FALSE,
      returnMatrix = TRUE
    )
    write.csv(heatmapGS, paste0("genes_per_condition_", i, "_hm.csv"))
  }
} else {
  heatmapGS <- "There are not enough Conditions to be compared with!"
}

# Volcano plots for genes
if (length(unique(proj$Condition)) > 1) {
  for (j in seq_along(treatment)) {

    ncells <- length(proj$cellNames)
    markerList <- getMarkerFeatures(
      ArchRProj = proj,
      useMatrix = "GeneScoreMatrix",
      groupBy = treatment[j],
      bias = c("TSSEnrichment", "log10(nFrags)"),
      maxCells = ncells,
      normBy = "none",
      testMethod = "ttest"
    )

    req_DF <- as.data.frame(getCellColData(proj))
    df1 <- table(req_DF$Clusters, req_DF[, treatment[j]])
    distr <- as.data.frame.matrix(round(prop.table(as.matrix(df1), 1), 2))

    lst <- list()
    for (i in 1:nrow(distr)) {
      row <- distr[i, ]
      if (
        sum(unname(unlist(row)) >= 0.90) == 1
      ) {
        lst[[i]] <- rownames(row)
      }
    }
    not_req_list <- unlist(lst)

    req_clusters <- unique(proj$Clusters)
    req_clusters <- req_clusters[
      order(as.numeric(gsub("C", "", req_clusters)))
    ]
    req_clusters <- req_clusters[which(!req_clusters %in% not_req_list)]
    markerList_C <- list()
    proj_C <- list()

    for (i in seq_along(req_clusters)) {

      idxSample <- BiocGenerics::which(proj$Clusters == req_clusters[i])
      cellsSample <- proj$cellNames[idxSample]
      proj_C[i] <- proj[cellsSample, ]
      ncells[i] <- length(proj_C[[i]]$cellNames)

      # per each cluster separately
      markerList_C[[i]] <- getMarkerFeatures(
        ArchRProj = proj_C[[i]],
        useMatrix = "GeneScoreMatrix",
        groupBy = treatment[j],
        bias = c("TSSEnrichment", "log10(nFrags)"),
        maxCells = ncells[[i]],
        normBy = "none",
        testMethod = "ttest"
      )
    }
    names(markerList_C) <- req_clusters

    gsm <- getMatrixFromProject(proj)
    gsm_mat <- assay(getMatrixFromProject(proj), "GeneScoreMatrix")

    which(rowSums(is.na(gsm_mat)) > 0)
    any(rowSums((gsm_mat)) == 0)

    empty_gene_idx <- which(rowSums((gsm_mat)) == 0)
    empty_gene <- rowData(gsm)$name[empty_gene_idx]

    markerList_df1 <- assay(markerList, "Log2FC")
    markerList_df2 <- assay(markerList, "Pval")
    markerList_df3 <- assay(markerList, "FDR")

    conditions <- sort(unique(proj@cellColData[treatment[j]][, 1]))
    markerList_df <- list()

    for (conds in conditions){
      markerList_df[[conds]] <- DataFrame(
        markerList_df1[, conds],
        markerList_df2[, conds],
        markerList_df3[, conds]
      )
      markerList_df[[conds]] <- as.data.frame(markerList_df[[conds]])
      markerList_df[[conds]]$genes <- rowData(markerList)$name
      markerList_df[[conds]]$cluster <- rep(
        "All", length(rownames(markerList_df[[conds]]))
      )
      colnames(markerList_df[[conds]]) <- c(
        "avg_log2FC", "p_val", "p_val_adj", "gene", "cluster"
      )
    }

    # per cluster
    markerList_df1_C <- list()
    markerList_df2_C <- list()
    markerList_df3_C <- list()
    markerList_df_C <- list()

    for (i in seq_along(req_clusters)) {

      cluster <- req_clusters[i]
      markerList_df1_C[[i]] <- assay(markerList_C[[i]], "Log2FC")
      markerList_df2_C[[i]] <- assay(markerList_C[[i]], "Pval")
      markerList_df3_C[[i]] <- assay(markerList_C[[i]], "FDR")

      conditions <- sort(unique(proj@cellColData[treatment[j]][, 1]))
      markerList_df_C[[i]] <- list()
      for (conds in conditions) {

        markerList_df_C[[i]][[conds]] <- DataFrame(
          markerList_df1_C[[i]][, conds],
          markerList_df2_C[[i]][, conds],
          markerList_df3_C[[i]][, conds]
        )

        markerList_df_C[[i]][[conds]] <- as.data.frame(
          markerList_df_C[[i]][[conds]]
        )
        markerList_df_C[[i]][[conds]]$genes<- rowData(markerList_C[[i]])$name
        markerList_df_C[[i]][[conds]]$cluster <- rep(
          cluster, dim(markerList_df_C[[i]][[conds]])[1]
        )
        colnames(markerList_df_C[[i]][[conds]]) <- c(
          "avg_log2FC", "p_val", "p_val_adj", "gene", "cluster"
        )
      }
    }

    names(markerList_df_C) <- req_clusters
    markersGS_merged_df <- do.call(Map, c(f = rbind, markerList_df_C))

    # also data frame for all clusters together needs to be added
    for (conds in conditions) {
      others <- paste(
        colnames(markerList)[conds != (colnames(markerList))], collapse = "|"
      )
      markersGS_merged_df[[conds]] <- rbind(
        markerList_df[[conds]], markersGS_merged_df[[conds]]
      )

      # remove empty genes
      markersGS_merged_df[[conds]] <- markersGS_merged_df[[conds]][which(
        !markersGS_merged_df[[conds]]$gene %in% empty_gene
      ), ]

      # remove na values
      markersGS_merged_df[[conds]] <- na.omit(markersGS_merged_df[[conds]])

      # remove FDR equal to 0
      markersGS_merged_df[[conds]] <- markersGS_merged_df[[conds]][
        which(!markersGS_merged_df[[conds]]$p_val_adj == 0),
      ]

      # make logfc limiation between 1 and -1
      markersGS_merged_df[[conds]] <- markersGS_merged_df[[conds]][
        which(abs(markersGS_merged_df[[conds]]$avg_log2FC) < 1.2),
      ]

      markersGS_merged_df[[conds]]$Significance <- ifelse(
        markersGS_merged_df[[conds]]$p_val < 10^-2,
        ifelse(
          markersGS_merged_df[[conds]]$avg_log2FC > 0.0,
          conds,
          others
        ),
        "Not siginficant"
      )

      de <- list()
      de[[conds]] <- markersGS_merged_df[[conds]]

      write.table(
        de[[conds]],
        paste0(
          "volcanoMarkers_genes_",
          j,
          "_",
          conds,
          ".csv"
        ),
        sep = ",",
        quote = FALSE,
        row.names = FALSE
      )

      print(
        paste0(
          "writing volcanoMarkers_genes_", j, "_", conds, ".csv is done!"
        )
      )

      features <- unique(de[[conds]]$cluster)
      volcano_plots <- list()
      for (i in seq_along(features)) {
        volcano_plots[[i]] <- scvolcano(
          de[[conds]], conds, others, features[[i]]
        )
      }

      pdf(paste0("volcano_plots_", conds, ".pdf"))
      for (plot in volcano_plots) {
        print(plot)
      }
      dev.off()
    }
  }

} else {
  de <- "There are not enough conditions to be compared with!"
}

tempdir <- "/root"
genes_per_cluster_hm <- find_func(tempdir, "genes_per_cluster_hm.csv")
hm_per_clust <- read.csv(genes_per_cluster_hm)

nClust <- length(unique(proj$Clusters))

df <- list()

for (i in seq_along(1: nClust)) {
  df[[i]] <- hm_per_clust[, c(1, i + 1)]

  #select top 5 values by group
  df[[i]] <- df[[i]][order(df[[i]][, 2], decreasing = TRUE), ][1:5, 1]
}

final <- do.call(rbind, df)
req_genes1 <- unlist(df)
req_genes1 <- req_genes1[!duplicated(req_genes1)]
req_genes1 <- na.omit(req_genes1)

# save the genes for default values in shiny app
write.csv(req_genes1, "req_genes1.csv")

# per sample
if (length(unique(proj$Sample)) > 1) {

  genes_per_sample_hm <- find_func(tempdir, "genes_per_sample_hm.csv")
  hm_per_sample <- read.csv(genes_per_sample_hm)

  nSamples <- length(unique(proj$Sample))

  df <- list()

  for (i in seq_along(1:nSamples)) {
    df[[i]] <- hm_per_sample[, c(1, i + 1)]

    #select top 10 values by group
    df[[i]] <- df[[i]][order(df[[i]][, 2], decreasing = TRUE), ][1:10, 1]
  }
  final <- do.call(rbind, df)

  req_genes3 <- unlist(df)
  req_genes3 <- req_genes3[!duplicated(req_genes3)]
  req_genes3 <- na.omit(req_genes3)
  write.csv(req_genes3, "req_genes3.csv")

} else {
  req_genes3 <- "There are not enough samples to be compared with!"
}

# per treatment
if (length(unique(proj$Condition)) > 1) {
  genes_per_cond_hm <- find_func(tempdir, "genes_per_condition_*")

  for (j in seq_along(genes_per_cond_hm)) {

    hm_per_cond <- read.csv(genes_per_cond_hm[j])
    nConds <- 2
    df <- list()
    for (i in seq_along(1:nConds)) {
      df[[i]] <- hm_per_cond[, c(1, i + 1)]

      #select top 20 values by group
      df[[i]] <- df[[i]][order(df[[i]][, 2], decreasing = TRUE), ][1:20, 1]
    }
    final <- do.call(rbind, df)

    req_genes2 <- unlist(df)
    req_genes2 <- req_genes2[!duplicated(req_genes2)]
    req_genes2 <- na.omit(req_genes2)
    write.csv(req_genes2, paste0("req_genes2_", j, ".csv"))
  }

} else {
  req_genes2 <- "There are not enough condition to be compared with!"
}

UMAPHarmony <- getEmbedding(
  ArchRProj = proj, embedding = "UMAP", returnDF = TRUE
)
write.csv(UMAPHarmony, "UMAPHarmony.csv")

saveArchRProject(
  ArchRProj = proj, outputDirectory = paste0(project_name, "_ArchRProject")
)

# sometimes we get outofMemory error; number of threads here:
addArchRThreads(threads = num_threads)

# peak calling with MACS2 for clusters-----------------------------------------

proj <- addGroupCoverages(
  ArchRProj = proj,
  groupBy = "Clusters",
  maxCells = 1500,
  force = TRUE
)

# get genome size
if (genome == "hg38") {
  genome_size <- 3.3e+09
} else if (genome == "mm10") {
  genome_size <- 3.0e+09
} else if (genome == "rnor6") {
  genome_size <- 2.9e+09
}

pathToMacs2 <- findMacs2()
proj <- addReproduciblePeakSet(
  ArchRProj = proj,
  groupBy = "Clusters",
  pathToMacs2 = pathToMacs2,
  genomeSize = genome_size,
  maxPeaks = 300000,
  force = TRUE
)
proj <- addPeakMatrix(proj, force = TRUE)

# save run data in medians.csv
metadata <- getCellColData(ArchRProj = proj)
tss <- aggregate(
  metadata@listData$TSSEnrichment,
  by = list(metadata@listData$Sample),
  FUN = median
)
nfrags <- aggregate(
  metadata@listData$nFrags,
  by = list(metadata@listData$Sample),
  FUN = median
)
conditions <- aggregate(
  metadata@listData$Condition,
  by = list(metadata@listData$Sample),
  FUN = max
)
frip <- aggregate(
  metadata@listData$FRIP,
  by = list(metadata@listData$Sample),
  FUN = median
)
frip$x <- round(frip$x, 4)
list_dfs <- list(tss, nfrags, frip, conditions)
merged_df <- Reduce(
  function(x, y) merge(x, y, by = "Group.1", all = TRUE), list_dfs
)
names(merged_df) <- c(
  "run_id", "median_TSS", "median_fragments", "median_FRIP", "condition"
)
write.csv(merged_df, file = "medians.csv", row.names = FALSE)

# Motif enrichment (Deviation)

proj <- add_motif_annotations(proj, genome) # from utils

saveArchRProject(
  ArchRProj = proj, outputDirectory = paste0(project_name, "_ArchRProject")
)

######################### get marker peaks, save ##############################

peak_data <- data.frame(proj@peakSet@ranges, proj@peakSet@elementMetadata)

# clusters
markers_peaks_c <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  k = 100,
  testMethod = "wilcoxon"
)

peak_marker_list_c <- getMarkers(
  markers_peaks_c, cutOff = "Pval <= 0.05 & Log2FC >= 0.1"
)
write.csv(
  peak_marker_list_c,
  file = "marker_peaks_per_cluster.csv",
  row.names = FALSE
)

total_peaks_c <- merge(peak_data, peak_marker_list_c, by = c("start", "end"))
write.csv(
  total_peaks_c, file = "complete_peak_list_cluster.csv", row.names = FALSE
)

# samples
markers_peaks_s <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  groupBy = "Sample",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  k = 100,
  testMethod = "wilcoxon"
)

peak_marker_list_s <- getMarkers(
  markers_peaks_s, cutOff = "Pval <= 0.05 & Log2FC >= 0.1"
)

write.csv(
  peak_marker_list_s,
  file = "marker_peaks_per_sample.csv",
  row.names = FALSE
)

total_peaks_s <- merge(peak_data, peak_marker_list_s, by = c("start", "end"))
write.csv(
  total_peaks_s, file = "complete_peak_list_sample.csv", row.names = FALSE
)

# treatment
if (length(unique(proj$Condition)) > 1) {
  for (i in seq_along(treatment)) {

    marker_peaks_t <- getMarkerFeatures(
      ArchRProj = proj,
      useMatrix = "PeakMatrix",
      groupBy = treatment[i],
      bias = c("TSSEnrichment", "log10(nFrags)"),
      k = 100,
      testMethod = "wilcoxon"
    )

    peak_marker_list_t <- getMarkers(
      marker_peaks_t, cutOff = "Pval <= 0.05 & Log2FC >= 0.1"
    )
    write.csv(
      peak_marker_list_t,
      file = paste0("marker_peaks_per_condition-", i, ".csv"),
      row.names = FALSE
    )

    total_peaks_t <- merge(
      peak_data, peak_marker_list_t, by = c("start", "end")
    )
    write.csv(
      total_peaks_t,
      file = paste0("complete_peak_list_condition-", i, ".csv"),
      row.names = FALSE
    )
  }
}

peak_cutoff <- "Pval <= 0.05 & Log2FC >= 0.1"
heatmap_peaks <- plotMarkerHeatmap(
  seMarker = markers_peaks_c,
  cutOff = peak_cutoff,
  transpose = TRUE
)

peak_hm <- ComplexHeatmap::draw(
  heatmap_peaks,
  heatmap_legend_side = "bot",
  annotation_legend_side = "bot",
  column_title = paste0("Marker peaks (", peak_cutoff, ")"),
  column_title_gp = gpar(fontsize = 12)
)

heatmaps[[2]] <- peak_hm

motifs_cutoff <- "Pval <= 0.05 & Log2FC >= 0.1"
enrichMotifs <- peakAnnoEnrichment(
  seMarker = markers_peaks_c,
  ArchRProj = proj,
  peakAnnotation = "Motif",
  cutOff = motifs_cutoff
)

motifs_df_c <- data.frame(enrichMotifs@assays@data)
write.csv(motifs_df_c, file = "enrichedMotifs_cluster.csv")

heatmapEM <- plotEnrichHeatmap(
  enrichMotifs,
  transpose = TRUE,
  n = 50,
  cutOff = 2
)

heatmap_motifs <- ComplexHeatmap::draw(
  heatmapEM,
  heatmap_legend_side = "bottom",
  column_title = paste0("Marker motifs (", motifs_cutoff),
  column_title_gp = gpar(fontsize = 12)
)

heatmaps[[3]] <- heatmap_motifs

print("+++++++++++creating heatmap plots++++++++++++++")

pdf("heatmaps_all.pdf")
for (i in seq_along(heatmaps)) {
  print(heatmaps[[i]])
}
dev.off()

###############################################################################

motif_lst <- unique(rownames(enrichMotifs))
split_string <- strsplit(motif_lst, split = "\\(")
fun1 <- function(list, nth) {
  sapply(list, `[`, 1)
}
req_motifs1 <- gsub("_", "-", fun1(split_string))
req_motifs1 <- gsub(" ", "", req_motifs1)

rownames(enrichMotifs) <- req_motifs1

saveRDS(enrichMotifs, "enrichMotifs_clusters.rds")

# cutOff A numeric cutOff that indicates the minimum P-adj enrichment to be
# included in the heatmap. default is 20 but we decrease that!

heatmapEM <- plotEnrichHeatmap(
  enrichMotifs,
  n = 50,
  transpose = FALSE,
  returnMatrix = TRUE,
  cutOff = 2
)

motif_lst <- unique(rownames(heatmapEM))
split_string <- strsplit(motif_lst, split = "\\(")
fun1 <- function(list, nth) {
  sapply(list, `[`, 1)
}
req_motifs1 <- gsub("_", "-", fun1(split_string))
req_motifs1 <- gsub(" ", "", req_motifs1)

rownames(heatmapEM) <- req_motifs1

write.csv(heatmapEM, "motif_per_cluster_hm.csv")

tempdir <- "/root"
motifs_per_cluster_hm <- find_func(tempdir, "motif_per_cluster_hm.csv")
hm_per_clust <- read.csv(motifs_per_cluster_hm)

nClust <- length(unique(proj$Clusters))

df <- list()

for (i in seq_along(1:nClust)) {
  df[[i]] <- hm_per_clust[, c(1, i + 1)]

  #select top 5 values by group
  df[[i]] <- df[[i]][order(df[[i]][, 2], decreasing = TRUE), ][1:5, 1]
}
final <- do.call(rbind, df)
req_motifs1 <- unlist(df)
req_motifs1 <- req_motifs1[!duplicated(req_motifs1)]
req_motifs1 <- na.omit(req_motifs1)

# save the motifs for default values in shiny app
write.csv(req_motifs1, "req_motifs1.csv")

######################## Add Motifs Matrix and Projections ####################

proj <- addBgdPeaks(proj, force = TRUE)
proj <- addDeviationsMatrix(
  ArchRProj = proj,
  peakAnnotation = "Motif",
  force = TRUE
)
markersMotifs <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "MotifMatrix",
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  useSeqnames = "z"
)
# create deviation score

markerMotifsList <- getMarkers(
  markersMotifs, cutOff = "FDR < 0.9 & MeanDiff >= 0"
)

motifs <- list()
for (i in seq_len(length(markerMotifsList))) {
  if (length(markerMotifsList[[i]]$name) > 1) {
    motifs <- c(motifs, markerMotifsList[[i]]$name)
  }
}

if (length(motifs) > 1) {
  motifs <- unlist(motifs)
  motifs <- paste0("z:", motifs)
  motifs <- unique(motifs)

  proj <- addImputeWeights(proj)

  dev_score <- getDeviation_ArchR(
    ArchRProj = proj,
    name = motifs,
    imputeWeights = getImputeWeights(proj)
  )

  dev_score[is.na(dev_score)] <- 0
}

dev_score2 <- dev_score[!is.infinite(rowSums(dev_score)), ]
colnames(dev_score2) <- rownames(getCellColData(proj))

# remove 0 deviations per All samples
all_zero <- names(which(rowSums(dev_score2) == 0))
dev_score2 <- dev_score2[which(!rownames(dev_score2) %in% c(all_zero)), ]

# convert to dgCmatrix
dev_score3 <- Matrix(as.matrix(dev_score2), sparse = TRUE)

# create metadata object for Seurat object
metadata <- getCellColData(ArchRProj = proj)
rownames(metadata) <- str_split_fixed(
  str_split_fixed(row.names(metadata), "#", 2)[, 2],
  "-",
  2
)[, 1]
metadata["log10_nFrags"] <- log(metadata$nFrags)

# create motif seurat objects -------------

seurat_objs <- c()

for (run in runs) {

  obj <- build_atlas_seurat_object(
    run_id = run[1],
    matrix = dev_score3,
    metadata = metadata,
    spatial_path = run[6]
  )

  saveRDS(obj, file = paste0(run[1], "_SeuratObjMotif.rds"))
  seurat_objs <- c(seurat_objs, obj)
}

print("These are the available seurat_objMotifs.")
seurat_objs

all_m <-  list()
for (i in seq_along(seurat_objs)) {
  all_m[[i]] <- seurat_objs[[i]]

  all_m[[i]] <- Seurat::RenameCells(
    all_m[[i]],
    new.names = paste0(
      unique(all_m[[i]]@meta.data$Sample),
      "#",
      colnames(all_m[[i]]),
      "-1"
    )
  )
}

# peak calling with MACS2 for Sample -------------------------------------------

if (length(unique(proj$Sample)) > 1) {

  proj <- addGroupCoverages(
    ArchRProj = proj,
    groupBy = "Sample",
    maxCells = 1500,
    force = TRUE
  )

  proj <- addReproduciblePeakSet(
    ArchRProj = proj,
    groupBy = "Sample",
    pathToMacs2 = pathToMacs2,
    genomeSize = genome_size,
    maxPeaks = 300000,
    force = TRUE
  )
  proj <- addPeakMatrix(proj, force = TRUE)

  proj <- add_motif_annotations(proj, genome) # from utils

  saveArchRProject(
    ArchRProj = proj, outputDirectory = paste0(project_name, "_ArchRProject")
  )

  markersPeaks <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix",
    groupBy = "Sample",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    k = 100,
    testMethod = "wilcoxon"
  )
  enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = proj,
    peakAnnotation = "Motif",
    cutOff = "Pval <= 0.05 & Log2FC >= 0.1"
  )

  motifs_df_s <- data.frame(enrichMotifs@assays@data)
  write.csv(motifs_df_s, file = "enrichedMotifs_sample.csv")

  motif_lst <- unique(rownames(enrichMotifs))
  split_string <- strsplit(motif_lst, split = "\\(")
  fun1 <- function(list, nth) {
    sapply(list, `[`, 1)
  }
  req_motifs3 <- gsub("_", "-", fun1(split_string))
  req_motifs3 <- gsub(" ", "", req_motifs3)

  rownames(enrichMotifs) <- req_motifs3

  saveRDS(enrichMotifs, "enrichMotifs_sample.rds")

  # cutOff: A numeric cutOff that indicates the minimum P-adj enrichment to be
  # included in the heatmap. default is 20 but we decrease that!
  heatmapEM <- plotEnrichHeatmap(
    enrichMotifs, n = 50, transpose = FALSE, returnMatrix = TRUE, cutOff = 2
  )

  motif_lst <- unique(rownames(heatmapEM))
  split_string <- strsplit(motif_lst, split = "\\(")
  fun1 <- function(list, nth) {
    sapply(list, `[`, 1)
  }
  req_motifs3 <- gsub("_", "-", fun1(split_string))
  req_motifs3 <- gsub(" ", "", req_motifs3)

  rownames(heatmapEM) <- req_motifs3
  write.csv(heatmapEM, "motif_per_sample_hm.csv")

  nSamples <- length(unique(proj$Sample))

  df <- list()
  tempdir <- "/root"
  motifs_per_sample_hm <- find_func(tempdir, "motif_per_sample_hm.csv")
  hm_per_sample <- read.csv(motifs_per_sample_hm)

  for (i in seq_along(1:nSamples)){
    df[[i]] <- hm_per_sample[, c(1, i + 1)]

    #select top 5 values by group
    df[[i]] <- df[[i]][order(df[[i]][, 2], decreasing = TRUE), ][1:10, 1]
  }
  final <- do.call(rbind, df)
  req_motifs3 <- unlist(df)

  req_motifs3 <- req_motifs3[!duplicated(req_motifs3)]

  req_motifs3 <- na.omit(req_motifs3)

  write.csv(req_motifs3, "req_motifs3.csv")

} else {
  enrichMotifs <- "There are not enough samples to be compared with!"
  heatmapEM <- "There are not enough samples to be compared with!"
  req_motifs3 <- "There are not enough samples to be compared with!"
}

# peak calling with MACS2 for treatment ---------------------------------------
if (length(unique(proj$Condition)) > 1) {
  for (i in seq_along(treatment)) {

    proj <- addGroupCoverages(
      ArchRProj = proj,
      groupBy = treatment[i],
      maxCells = 1500,
      force = TRUE
    )

    proj <- addReproduciblePeakSet(
      ArchRProj = proj,
      groupBy = treatment[i],
      pathToMacs2 = pathToMacs2,
      genomeSize = genome_size,
      maxPeaks = 300000,
      force = TRUE
    )
    proj <- addPeakMatrix(proj, force = TRUE)

    proj <- add_motif_annotations(proj, genome) # from utils

    # save ArchR object
    saveArchRProject(
      ArchRProj = proj, outputDirectory = paste0(project_name, "_ArchRProject")
    )
  }

  for (i in seq_along(treatment)) {

    markersPeaks <- getMarkerFeatures(
      ArchRProj = proj,
      useMatrix = "PeakMatrix",
      groupBy = treatment[i],
      bias = c("TSSEnrichment", "log10(nFrags)"),
      k = 100,
      testMethod = "wilcoxon"
    )

    enrichMotifs <- peakAnnoEnrichment(
      seMarker = markersPeaks,
      ArchRProj = proj,
      peakAnnotation = "Motif",
      cutOff = "Pval <= 0.05 & Log2FC >= 0.1"
    )

    motifs_df_t <- data.frame(enrichMotifs@assays@data)
    write.csv(
      motifs_df_t, file = paste0("enrichedMotifs_condition_", i, ".csv")
    )

    motif_lst <- unique(rownames(enrichMotifs))
    split_string <- strsplit(motif_lst, split = "\\(")
    fun1 <- function(list, nth) {
      sapply(list, `[`, 1)
    }
    req_motifs2 <- gsub("_", "-", fun1(split_string))
    req_motifs2 <- gsub(" ", "", req_motifs2)

    rownames(enrichMotifs) <- req_motifs2
    saveRDS(enrichMotifs, paste0("enrichMotifs_condition_", i, ".rds"))

    # cutOff A numeric cutOff that indicates the minimum P-adj enrichment to
    # be included in the heatmap. Default is 20 but we decrease that!

    heatmapEM <- plotEnrichHeatmap(
      enrichMotifs, n = 50, transpose = FALSE, returnMatrix = TRUE, cutOff = 2
    )

    motif_lst <- unique(rownames(heatmapEM))
    split_string <- strsplit(motif_lst, split = "\\(")
    fun1 <- function(list, nth) {
      sapply(list, `[`, 1)
    }
    req_motifs2 <- gsub("_", "-", fun1(split_string))
    req_motifs2 <- gsub(" ", "", req_motifs2)

    rownames(heatmapEM) <- req_motifs2
    write.csv(heatmapEM, paste0("motif_per_condition_", i, "_hm.csv"))
  }

  nConds <- 2
  df <- list()
  tempdir <- "/root"
  motifs_per_cond_hm <- find_func(tempdir, "motif_per_condition_*")

  for (j in seq_along(motifs_per_cond_hm)) {
    hm_per_cond <- read.csv(motifs_per_cond_hm[j])

    for (i in seq_along(1:nConds)) {
      df[[i]] <- hm_per_cond[, c(1, i + 1)]

      #select top 5 values by group
      df[[i]] <- df[[i]][order(df[[i]][, 2], decreasing = TRUE), ][1:5, 1]
    }

    final <- do.call(rbind, df)
    req_motifs2 <- unlist(df)

    req_motifs2 <- req_motifs2[!duplicated(req_motifs2)]
    req_motifs2 <- na.omit(req_motifs2)
    write.csv(req_motifs2, paste0("req_motifs2_", j, ".csv"))
  }
} else {
  enrichMotifs <- "There are not enough conditions to be compared with!"
  heatmapEM <- "There are not enough conditions to be compared with!"
  req_motifs2 <- "There are not enough conditions to be compared with!"
}

# Volcano plots for motifs -----------------------------------------------------
if (length(unique(proj$Condition)) > 1) {
  for (j in seq_along(treatment)) {

    ncells <- length(proj$cellNames)

    # all clusters together
    markersMotifs <- getMarkerFeatures(
      ArchRProj = proj,
      useMatrix = "MotifMatrix",
      groupBy = treatment[j],
      bias = c("TSSEnrichment", "log10(nFrags)"),
      testMethod = "wilcoxon",
      useSeqnames = "z",
      maxCells = 5000,
      normBy = "none"
    )

    req_DF <- as.data.frame(getCellColData(proj))
    df1 <- table(req_DF$Clusters, req_DF[, treatment[j]])
    distr <- as.data.frame.matrix(round(prop.table(as.matrix(df1), 1), 2))
    lst <- list()

    for (i in 1:nrow(distr)) {
      row <- distr[i, ]
      if (
        sum(unname(unlist(row)) >= 0.90) == 1) {
        rownames(row) -> lst[[i]]
      }
    }
    not_req_list <- unlist(lst)

    req_clusters <- unique(proj$Clusters)
    req_clusters <- req_clusters[order(as.numeric(gsub("C", "", req_clusters)))]
    req_clusters <- req_clusters[which(!req_clusters %in% not_req_list)]

    markersMotifs_C <- list()
    proj_C <- list()

    for (i in seq_along(req_clusters)) {

      idxSample <- BiocGenerics::which(proj$Clusters == req_clusters[i])

      cellsSample <- proj$cellNames[idxSample]
      proj_C[i] <- proj[cellsSample, ]

      ncells[i] <- length(proj_C[[i]]$cellNames)

      # per each cluster separately
      markersMotifs_C[[i]] <- getMarkerFeatures(
        ArchRProj = proj_C[[i]],
        useMatrix = "MotifMatrix",
        groupBy = treatment[j],
        bias = c("TSSEnrichment", "log10(nFrags)"),
        maxCells = ncells[[i]],
        normBy = "none",
        testMethod = "ttest"
      )
    }
    names(markersMotifs_C) <- req_clusters
    dev_score <- getDeviation_ArchR(
      ArchRProj = proj,
      name = motifs,
      imputeWeights = getImputeWeights(proj)
    )
    empty_motif_idx <- which(rowSums((dev_score)) == 0)
    empty_motif <- rownames(dev_score)[empty_motif_idx] 

    markersMotifs_df1 <- assay(markersMotifs, "MeanDiff")
    markersMotifs_df2 <- assay(markersMotifs, "Pval")
    markersMotifs_df3 <- assay(markersMotifs, "FDR")

    conditions <- sort(unique(proj@cellColData[treatment[j]][, 1]))
    markersMotifs_df <- list()

    for (conds in conditions){

      markersMotifs_df[[conds]] <- DataFrame(
        markersMotifs_df1[[conds]],
        markersMotifs_df2[[conds]],
        markersMotifs_df3[[conds]]
      )
      markersMotifs_df[[conds]] <- as.data.frame(markersMotifs_df[[conds]])
      markersMotifs_df[[conds]]$genes <- rowData(markersMotifs)$name
      markersMotifs_df[[conds]]$cluster <- rep(
        "All",
        length(rownames(markersMotifs_df[[conds]]))
      )
      colnames(markersMotifs_df[[conds]]) <- c(
        "avg_log2FC", "p_val", "p_val_adj", "gene", "cluster"
      )
    }

    # percluster
    markersMotifs_df1_C <- list()
    markersMotifs_df2_C <- list()
    markersMotifs_df3_C <- list()
    markersMotifs_df_C <- list()

    for (i in seq_along(req_clusters)) {

      cluster <- req_clusters[i]

      markersMotifs_df1_C[[i]] <- assay(markersMotifs_C[[i]], "MeanDiff")
      markersMotifs_df2_C[[i]] <- assay(markersMotifs_C[[i]], "Pval")
      markersMotifs_df3_C[[i]] <- assay(markersMotifs_C[[i]], "FDR")

      conditions <- sort(unique(proj@cellColData[treatment[j]][, 1]))
      markersMotifs_df_C[[i]] <- list()
      for (conds in conditions){

        markersMotifs_df_C[[i]][[conds]] <- DataFrame(
          markersMotifs_df1_C[[i]][, conds],
          markersMotifs_df2_C[[i]][, conds],
          markersMotifs_df3_C[[i]][, conds]
        )
        markersMotifs_df_C[[i]][[conds]] <- as.data.frame(
          markersMotifs_df_C[[i]][[conds]]
        )
        markersMotifs_df_C[[i]][[conds]]$genes <- rowData(
          markersMotifs_C[[i]]
        )$name
        markersMotifs_df_C[[i]][[conds]]$cluster <- rep(
          cluster, dim(markersMotifs_df_C[[i]][[conds]])[1]
        )
        colnames(markersMotifs_df_C[[i]][[conds]]) <- c(
          "avg_log2FC", "p_val", "p_val_adj", "gene", "cluster"
        )
      }
    }
    names(markersMotifs_df_C) <- req_clusters
    markersMotifs_merged_df <- do.call(Map, c(f = rbind, markersMotifs_df_C))

    # also data frame for all clusters together needs to be added
    for (conds in conditions) {
      others <- paste(
        colnames(markersMotifs)[conds != (colnames(markersMotifs))],
        collapse = "|"
      )

      markersMotifs_merged_df[[conds]] <- rbind(
        markersMotifs_df[[conds]], markersMotifs_merged_df[[conds]]
      )

      # remove empty genes
      markersMotifs_merged_df[[conds]] <- markersMotifs_merged_df[[conds]][
        which(!markersMotifs_merged_df[[conds]]$gene %in% empty_motif),
      ]

      # remove na values
      markersMotifs_merged_df[[conds]] <- na.omit(
        markersMotifs_merged_df[[conds]]
      )

      # remove FDR equal to 0
      markersMotifs_merged_df[[conds]] <- markersMotifs_merged_df[[conds]][
        which(!markersMotifs_merged_df[[conds]]$p_val_ad == 0),
      ]

      # make logfc limiation between 1 and -1
      markersMotifs_merged_df[[conds]] <- markersMotifs_merged_df[[conds]][
        which(abs(markersMotifs_merged_df[[conds]]$avg_log2FC) < 1.2),
      ]

      markersMotifs_merged_df[[conds]]$Significance <- ifelse(
        markersMotifs_merged_df[[conds]]$p_val < 10^-2,
        ifelse(
          markersMotifs_merged_df[[conds]]$avg_log2FC > 0.0,
          conds,
          others
        ),
        "Not siginficant"
      )
      de <- list()
      de[[conds]] <- markersMotifs_merged_df[[conds]]

      write.table(
        de[[conds]],
        paste0("volcanoMarkers_motifs_", j, "_", conds, ".csv"),
        sep = ",",
        quote = FALSE,
        row.names = FALSE
      )
      print(
        paste0(
          "writing volcanoMarkers_motifs_", j, "_", conds, ".csv is done!"
        )
      )

      features_m <- unique(de[[conds]]$cluster)
      volcano_plots_m <- list()
      for (i in seq_along(features_m)) {
        volcano_plots_m[[i]] <- scvolcano(
          de[[conds]],  conds, others, features_m[[i]]
        )
      }

      pdf(paste0("volcano_plots_motifs_", j, "_", conds, ".pdf"))
      for (plot in volcano_plots_m) {
        print(plot)
      }
      dev.off()
    }
  }
} else {
  de <- "There are not enough conditions to be compared with!"
}

################-------------- save bigwig files -------- ######################

req_conditions <- c("Clusters", treatment)

for (i in req_conditions) {
  bws <- getGroupBW(
    ArchRProj = proj,
    groupBy = i,
    normMethod = "ReadsInTSS",
    tileSize = 100,
    maxCells = 1000,
    ceiling = 4,
    verbose = TRUE,
    threads = getArchRThreads(),
    logFile = createLogFile("getGroupBW")
  )
}

saveArchRProject(
  ArchRProj = proj, outputDirectory = paste0(project_name, "_ArchRProject")
)

#################------------- Motif Logo ---------------#######################

PWMs <- getPeakAnnotation(proj, "Motif")$motifs

PWMatrixToProbMatrix <- function(x) {
  if (class(x) != "PWMatrix") stop("x must be a TFBSTools::PWMatrix object")
  (exp(as(x, "matrix"))) * TFBSTools::bg(x) / sum(TFBSTools::bg(x))
}

ProbMatrices <- lapply(PWMs, PWMatrixToProbMatrix)
lapply(ProbMatrices, colSums) %>% range

PWMatrixToProbMatrix <- function(x) {
  if (class(x) != "PWMatrix") stop("x must be a TFBSTools::PWMatrix object")
  m <- (exp(as(x, "matrix"))) * TFBSTools::bg(x) / sum(TFBSTools::bg(x))
  m <- t(t(m) / colSums(m))
  m
}

ProbMatrices <- lapply(PWMs, PWMatrixToProbMatrix)
lapply(ProbMatrices, colSums) %>% range

saveRDS(ProbMatrices, "seqlogo.rds")

##############------------- Combined SeuratObjs ---------------################

samples <- find_samples_name(all)

# extract image coordinates as -(imagecols) | imagerow
spatial <- lapply(all, function(x) {
  df <- as.data.frame(x@images[[1]]@coordinates[, c(5, 4)])
  colnames(df) <- paste0("Spatial_", 1:2)
  df$Spatial_2 <- -df$Spatial_2
  df
})

combined <- combine_objs(all, UMAPHarmony, samples, spatial, project_name)
combined_m <- combine_objs(all_m, UMAPHarmony, samples, spatial, project_name)

saveRDS(combined, "combined.rds", compress = FALSE)
saveRDS(combined_m, "combined_m.rds", compress = FALSE)

##############------------- Shiny app files ---------------################
print("Shiny App starting...")

scConf1 <- createConfig(combined)
makeShinyFiles(
  combined,
  scConf1,
  gex.assay = "scATAC",
  gex.slot = "data",
  gene.mapping = TRUE,
  shiny.prefix = "sc1",
  default.gene1 = "Tiam1",
  default.gene2 = "Ccbe1",
  default.dimred = c("UMAP_1", "UMAP_2")
)

scConf2 <- createConfig(combined_m)
makeShinyFiles(
  combined_m,
  scConf2,
  gex.assay = "scATAC",
  gex.slot = "counts",
  gene.mapping = TRUE,
  shiny.prefix = "sc2",
  default.gene1 = "RFX3-1018",
  default.gene2 = "NEUROG2-1580",
  default.dimred = c("UMAP_1", "UMAP_2")
)
citation <- list(
  title <- paste0(project_name, " Data Analysis")
)
makeShinyCodesMulti(
  shiny.title = paste0(project_name, "_Lab Data Analysis"),
  shiny.footnotes = citation,
  shiny.prefix = c("sc1", "sc2"),
  shiny.headers = c("Gene Accessibility", "Peak/Motifs"),
  shiny.dir = "./shinyApp"
)

# edit some of the prepared data
sc1def <- readRDS("/root/shinyApp/sc1def.rds")
sc2def <- readRDS("/root/shinyApp/sc2def.rds")

xlim <- lapply(spatial, function(x) {
  xlim <- c(min(x[, 1]), max(x[, 1]))
  xlim
})
ylim <- lapply(spatial, function(y) {
  ylim <- c(min(y[, 2]), max(y[, 2]))
  ylim
})

sc1def$limits <- list()
for (i in seq_along(samples)) {
  sc1def[["limits"]][[samples[i]]] <- c(
    min(xlim[[i]]), max(xlim[[i]]), min(ylim[[i]]), max(ylim[[i]])
  )
}

sc2def$limits <- list()
for (i in seq_along(samples)) {
  sc2def[["limits"]][[samples[i]]] <- c(
    min(xlim[[i]]), max(xlim[[i]]), min(ylim[[i]]), max(ylim[[i]])
  )
}

sc1def$meta1 <- "Clusters"
sc1def$meta2 <- "Sample"
sc1def$meta3 <- "SampleName"

sc2def$meta1 <- "Clusters"
sc2def$meta2 <- "Sample"
sc2def$meta3 <- "SampleName"

sc1def$Clusters <- req_genes1
sc1def$Sample <- req_genes3

sc2def$Clusters <- req_motifs1
sc2def$Sample <- req_motifs3

sc1def$dimred[3] <- paste0(names(combined@reductions)[1], "_1")
sc1def$dimred[4] <- paste0(names(combined@reductions)[1], "_2")
sc1def$dimred[5] <- paste0(names(combined@reductions)[2], "_1")
sc1def$dimred[6] <- paste0(names(combined@reductions)[2], "_2")

sc2def$dimred[3] <- paste0(names(combined_m@reductions)[1], "_1")
sc2def$dimred[4] <- paste0(names(combined_m@reductions)[1], "_2")
sc2def$dimred[5] <- paste0(names(combined_m@reductions)[2], "_1")
sc2def$dimred[6] <- paste0(names(combined_m@reductions)[2], "_2")

if (length(unique(proj$Condition)) > 1) {
  for (i in seq_along(treatment)) {

    sc1def[[paste0("meta", (2 + i))]] <- treatment[i]
    sc1def[[treatment[i]]] <- read.csv(
      find_func(tempdir, paste0("req_genes2_", i, ".csv"))
    )$x

    sc1def[[paste0(treatment[i], "_1")]] <- sort(
      unique(combined@meta.data[[treatment[i]]])
    )[1]

    sc1def[[paste0(treatment[i], "_2")]] <- sort(
      unique(combined@meta.data[[treatment[i]]])
    )[2]

    sc2def[[paste0("meta", (2 + i))]] <- treatment[i]

    sc2def[[treatment[i]]] <- read.csv(
      find_func(tempdir, paste0("req_motifs2_", i, ".csv"))
    )$x

    sc2def[[paste0(treatment[i], "_1")]] <- sort(
      unique(combined_m@meta.data[[treatment[i]]])
    )[1]

    sc2def[[paste0(treatment[i], "_2")]] <- sort(
      unique(combined_m@meta.data[[treatment[i]]])
    )[2]
  }
} else {
  print("There are not enough conditions to add to sc1def.rds or sc2def.rds")
}

saveRDS(sc1def, "/root/shinyApp/sc1def.rds")
saveRDS(sc2def, "/root/shinyApp/sc2def.rds")

sc1conf <- readRDS("/root/shinyApp/sc1conf.rds")
sc2conf <- readRDS("/root/shinyApp/sc2conf.rds")

fav <- which(sc1conf$ID %in% c("Clusters", "Sample", "SampleName", treatment))
rest <- which(!sc1conf$ID %in% c("Clusters", "Sample", "SampleName", treatment))
sc1conf <- sc1conf[c(fav, rest), ]

fav <- which(sc2conf$ID %in% c("Clusters", "Sample", "SampleName", treatment))
rest <- which(!sc2conf$ID %in% c("Clusters", "Sample", "SampleName", treatment))
sc2conf <- sc2conf[c(fav, rest), ]

saveRDS(sc1conf, "/root/shinyApp/sc1conf.rds")
saveRDS(sc2conf, "/root/shinyApp/sc2conf.rds")

number_of_pixle <- c()
for (run in runs) {
  positions <- read.csv(run[5], header = FALSE)
  number_of_pixle <- c(number_of_pixle, length(positions$V1))
}
max_number_of_pixle <- max(number_of_pixle)
print(paste0("Max number of pixels: ", max_number_of_pixle))

print("Make sure there is no ui.R or server.R file in /root")
file.remove(
  list.files(
    path = "/root",
    pattern = "^ui.*\\.R$",
    all.files = FALSE,
    full.names = FALSE,
    recursive = FALSE,
    include.dirs = FALSE
  )
)
file.remove(
  list.files(
    path = "/root",
    pattern = "^server.*\\.R$",
    all.files = FALSE,
    full.names = FALSE,
    recursive = FALSE,
    include.dirs = FALSE
  )
)

if (max_number_of_pixle <= 2500) {
  print("use ui/server from 50by50 folder")
  unlink("/root/uiserver96by96", recursive = TRUE) # will delete directory

  if (length(unique(proj$Sample)) == 1) {
    file.copy("/root/uiserver50by50/ui_3.R", "/root/ui.R", overwrite = TRUE)
    file.copy(
      "/root/uiserver50by50/server_3.R", "/root/server.R", overwrite = TRUE
    )
    unlink("/root/uiserver50by50", recursive = TRUE) # will delete directory
  } else if (length(unique(proj$Condition)) <= 1) {
    file.copy("/root/uiserver50by50/ui_2.R", "/root/ui.R", overwrite = TRUE)
    file.copy(
      "/root/uiserver50by50/server_2.R", "/root/server.R", overwrite = TRUE
    )
    unlink("/root/uiserver50by50", recursive = TRUE) # will delete directory
  } else if (
    length(unique(proj$Condition)) > 1 & length(unique(proj$condition_1)) > 2 |
      length(unique(proj$Condition)) > 1 & length(unique(proj$condition_2)) > 2
  ) {
    file.copy("/root/uiserver50by50/ui_4.R", "/root/ui.R", overwrite = TRUE)
    file.copy(
      "/root/uiserver50by50/server_4.R", "/root/server.R", overwrite = TRUE
    )
    unlink("/root/uiserver50by50", recursive = TRUE) # will delete directory
  } else {
    file.copy("/root/uiserver50by50/ui.R", "/root/ui.R", overwrite = TRUE)
    file.copy(
      "/root/uiserver50by50/server.R", "/root/server.R", overwrite = TRUE
    )
    unlink("/root/uiserver50by50", recursive = TRUE) # will delete directory
  }

} else { #if it is more than 50 by 50 use ui/server from flowgel folder
  print("use ui/server from 96by96 folder")
  unlink("/root/uiserver50by50", recursive = TRUE) # will delete directory

  if (length(unique(proj$Sample)) == 1) {
    file.copy("/root/uiserver96by96/ui_3.R", "/root/ui.R", overwrite = TRUE)
    file.copy(
      "/root/uiserver96by96/server_3.R", "/root/server.R", overwrite = TRUE
    )
    unlink("/root/uiserver96by96", recursive = TRUE) # will delete directory
  } else if (length(unique(proj$Condition)) <= 1) {
    file.copy("/root/uiserver96by96/ui_2.R", "/root/ui.R", overwrite = TRUE)
    file.copy(
      "/root/uiserver96by96/server_2.R", "/root/server.R", overwrite = TRUE
    )
    unlink("/root/uiserver96by96", recursive = TRUE) # will delete directory
  } else if (
    length(unique(proj$Condition)) > 1 & length(unique(proj$condition_1)) > 2 |
    length(unique(proj$Condition)) > 1 & length(unique(proj$condition_2)) > 2
  ) {
    file.copy("/root/uiserver96by96/ui_4.R", "/root/ui.R", overwrite = TRUE)
    file.copy(
      "/root/uiserver96by96/server_4.R", "/root/server.R", overwrite = TRUE
    )
    unlink("/root/uiserver96by96", recursive = TRUE) # will delete directory
  } else {
    file.copy("/root/uiserver96by96/ui.R", "/root/ui.R", overwrite = TRUE)
    file.copy(
      "/root/uiserver96by96/server.R", "/root/server.R", overwrite = TRUE
    )
    unlink("/root/uiserver96by96", recursive = TRUE) # will delete directory
  }
}

# copy everything in shiny app to the root

rawPath <- "/root/shinyApp/"
dataPath <- "/root/"

dataFiles <- dir(rawPath, "*.rds$", ignore.case = TRUE, all.files = TRUE)
file.copy(file.path(rawPath, dataFiles), dataPath, overwrite = TRUE)

dataFiles <- dir(rawPath, "*.h5$", ignore.case = TRUE, all.files = TRUE)
file.copy(file.path(rawPath, dataFiles), dataPath, overwrite = TRUE)

# Convert Seurat to h5ad and save ----
for (obj in all) {
  seurat_to_h5ad(obj, FALSE, paste0(unique(obj$Sample), "_g"))  # from seurat.R
}

# Convert Seurat to h5ad and save ----
for (obj in all_m) {
  seurat_to_h5ad(obj, FALSE, paste0(unique(obj$Sample), "_m"))  # from seurat.R
}
