library("ArchR")
library("BSgenome")
library("BSgenome.Hsapiens.UCSC.hg38")
library("BSgenome.Mmusculus.UCSC.mm10")
library("BSgenome.Rnorvegicus.UCSC.rn6")
library("chromVARmotifs")
library("circlize")
library("ComplexHeatmap")
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

source("/root/wf/archr.R")
source("/root/wf/seurat.R")
source("/root/wf/utils.R")

source("/root/getDeviation_ArchR.R")
source("/root/makeShinyFiles.R")


# Globals ---------------------------------------------------------------------
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

print(paste("Number of threads for peak calling:", num_threads))

runs <- strsplit(args[14:length(args)], ",")
runs

inputs <- c()
for (run in runs) {
  inputs[run[1]] <- run[3]
}
inputs

out_dir <- paste0(project_name, "_ArchRProject")

# Save input metrics in csv -----
metrics_to_csv(as.list(args[1:13])) # from utils.R

# Set genome size for peak calling -----
genome_sizes <- list("hg38" = 3.3e+09, "mm10" = 3.0e+09, "rnor6" = 2.9e+09)
genome_size <- genome_sizes[[genome]]

tempdir <- "/root"

# Create ArchRProject ---------------------------------------------------------
addArchRThreads(threads = 50)

proj <- create_archrproject( # from archr.R
  inputs, genome, min_tss, min_frags, tile_size, out_dir
)

# Add an additional Conditions column -----
for (run in runs) {
  proj$Condition[proj$Sample == run[1]] <- run[4]
  proj$SampleName[proj$Sample == run[1]] <- run[2]
}

# Filter on-tissue -----
all_ontissue <- c()
for (run in runs) {
  positions <- read.csv(run[5], header = FALSE)
  positions$V1 <- paste0(run[1], "#", positions$V1, "-1")
  on_tissue <- positions$V1 [which(positions$V2 == 1)]
  all_ontissue <- c(all_ontissue, on_tissue)
}

proj <- proj[proj$cellNames %in% all_ontissue]

# Parse conditions into 'treatments', add as columns to CellColData -----
conds <- strsplit(proj$Condition, split = "\\s|-")

for (i in seq_along(conds[[1]])) {
  proj <- ArchR::addCellColData(
    proj,
    data = extract_nth_ele(conds, i),
    name = paste0("condition_", i),
    cells = proj$cellNames,
    force = TRUE
  )
}
treatment <- names(ArchR::getCellColData(proj))[
  grep("condition_", names(ArchR::getCellColData(proj)))
]

print("++++ What treatments are in this project? ++++")
treatment

# Set global values for project data -----
n_samples <- length(unique(proj$Sample))
n_cond <- length(unique(proj$Condition))
n_cells <- length(proj$cellNames)

saveArchRProject(
  ArchRProj = proj, outputDirectory = paste0(project_name, "_ArchRProject")
)

# Dimensionality reduction and clustering -------------------------------------
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

# Batch correction with Harmony if more than one run -----
if (length(runs) > 1) {
  proj <- addHarmony(
    ArchRProj = proj,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample",
    force = TRUE
  )
  rd_name <- "Harmony"
} else {
  rd_name <- "IterativeLSI"
}

proj <- ArchR::addClusters(
  input = proj,
  reducedDims = rd_name,
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
  reducedDims = rd_name,
  name = "UMAP",
  nNeighbors = 30,
  minDist = umap_mindist,
  metric = "cosine",
  force = TRUE
)

proj <- ArchR::addImputeWeights(proj)

saveArchRProject(
  ArchRProj = proj, outputDirectory = paste0(project_name, "_ArchRProject")
)

# Plot UMAP embedding colored by Cluster, Sample, Treatment ----- *plot*
save_umap(proj, c("Clusters", "Sample", treatment))

# Create SeuratObjects for gene matrix ----------------------------------------

# create metadata object for Seurat object -----
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

# Extract gene matrix for Seurat object -----
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

print("Creating SeuratObjects...")

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

print("Available SeuratObjects...")
seurat_objs

# Plot clusters ontop of spatial coordinates -----
print("Creating spatial cluster plots...")
save_spatial_cluster_plots(seurat_objs)

# Create QC plots for TSSEnrichment, nFrags, log10_nFrags -----
print("Creating QC plots...")
save_qc_plots(seurat_objs, c("TSSEnrichment", "nFrags", "log10_nFrags"))

# Save ArchR object -----
saveArchRProject(
  ArchRProj = proj, outputDirectory = paste0(project_name, "_ArchRProject")
)

# Identify marker genes -------------------------------------------------------

# Marker genes per cluster, save marker gene rds, csv, heatmap.csv-----
cluster_marker_genes <- get_marker_genes( # from archr.R
  proj,
  group_by = "Clusters",
  markers_cutoff = "FDR <= 0.02",
  heatmap_cutoff = "Pval <= 0.05 & Log2FC >= 0.10",
  rd_name = rd_name
)

saveRDS(cluster_marker_genes$markers_gs, "markersGS_clusters.rds")
write.csv(
  cluster_marker_genes$marker_list,
  "marker_genes_per_cluster.csv",
  row.names = FALSE
)
write.csv(cluster_marker_genes$heatmap_gs, "genes_per_cluster_hm.csv")

# Marker genes per sample, save marker gene rds, csv, heatmap.csv-----
if (n_samples > 1) {

  sample_marker_genes <- get_marker_genes(
    proj,
    group_by = "Sample",
    markers_cutoff = "FDR <= 0.02",
    heatmap_cutoff = "Pval <= 0.05 & Log2FC >= 0.10",
    rd_name = rd_name
  )

  saveRDS(sample_marker_genes$markers_gs, "markersGS_sample.rds")
  write.csv(
    sample_marker_genes$marker_list,
    "marker_genes_per_sample.csv",
    row.names = FALSE
  )
  write.csv(sample_marker_genes$heatmap_gs, "genes_per_sample_hm.csv")

} else {
  print("There are not enough samples to be compared with!")
}

# Marker genes per treatment, save marker gene rds, csv, heatmap.csv-----
if (n_cond > 1) {

  for (i in seq_along(treatment)) {

    treatment_marker_genes <- get_marker_genes(
      proj,
      group_by = treatment[i],
      markers_cutoff = "FDR <= 0.02",
      heatmap_cutoff = "Pval <= 0.05 & Log2FC >= 0.10",
      rd_name = rd_name
    )

    saveRDS(
      treatment_marker_genes$markers_gs,
      paste0("markersGS_condition_", i, ".rds")
    )
    write.csv(
      treatment_marker_genes$marker_list,
      paste0("marker_genes_per_condition", i, ".csv"),
      row.names = FALSE
    )
    write.csv(
      treatment_marker_genes$heatmap_gs,
      paste0("genes_per_conditions", i, "_hm.csv")
    )
  }
} else {
  print("There are not enough Conditions to be compared with!")
}

# Initiate heatmaps list and plot marker gene per cluster heatmap -----
heatmaps <- list()

# Recompute cluster gs heatmap for plotting (transpose, plotLog2FC different)
gene_cutoff <- "Pval <= 0.05 & Log2FC >= 0.1"
heatmap_gs_plotting <- plotMarkerHeatmap(
  seMarker = cluster_marker_genes$markers_gs,
  cutOff = gene_cutoff,
  transpose = TRUE
)

# Save for plotting with peaks and motifs -----
gene_hm <- ComplexHeatmap::draw(
  heatmap_gs_plotting,
  heatmap_legend_side = "bot",
  annotation_legend_side = "bot",
  column_title = paste0("Marker genes (", gene_cutoff, ")"),
  column_title_gp = gpar(fontsize = 12)
)
heatmaps[[1]] <- gene_hm

# Volcano plots for genes -----------------------------------------------------
if (n_cond > 1) {
  for (j in seq_along(treatment)) {

    # Get gene markers for all clusters together -----
    marker_genes_df <- get_marker_df(
      proj = proj,
      group_by = treatment[j],
      matrix = "GeneScoreMatrix",
      seq_names = NULL,
      max_cells = n_cells,
      test_method = "ttest"
    )

    # Create a DF from marker genes for clusters for which no condition is
    # >90% of all cells -----
    req_clusters <- get_required_clusters(proj, treatment[j])
    marker_genes_by_cluster_df <- get_marker_df_clusters(
      proj, req_clusters, treatment[j]
    )

    # Merge and cleanup data per condition -----
    conditions <- sort(unique(proj@cellColData[treatment[j]][, 1]))
    for (cond in conditions) {

      volcano_table <- get_volcano_table( # from archr.R
        marker_genes_df, marker_genes_by_cluster_df, cond, "gene", gene_matrix
      )

      write.table(
        volcano_table,
        paste0(
          "volcanoMarkers_genes_", j, "_", cond, ".txt"
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
      )

      print(
        paste0(
          "writing volcanoMarkers_genes_", j, "_", cond, ".txt is done!"
        )
      )

      features <- unique(volcano_table$cluster)
      others <- paste(conditions[conditions != cond], collapse = "|")
      volcano_plots <- list()
      for (i in seq_along(features)) {
        volcano_plots[[i]] <- scvolcano(
          volcano_table, cond, others, features[[i]]
        )
      }

      pdf(paste0("volcano_plots_", cond, ".pdf"))
      for (plot in volcano_plots) {
        print(plot)
      }
      dev.off()
    }
  }
} else {
  print("There are not enough conditions to be compared with!")
}

# Save top marker genes for cluster, sample, condition for ShinyApp  ----------

# Top 5 marker genes by cluster -----
n_clust <- length(unique(proj$Clusters))

# Read hm back into memory as df with feats as column
hm_per_cluster <- read.csv("genes_per_cluster_hm.csv")

req_genes1 <- get_topn_hm_feats(hm_per_cluster, n_clust, 5)
write.csv(req_genes1, "req_genes1.csv")

# Top 10 marker genes per sample -----
if (n_samples > 1) {

  # Read hm back into memory as df with feats as column
  hm_per_sample <- read.csv("genes_per_sample_hm.csv")

  req_genes3 <- get_topn_hm_feats(hm_per_sample, n_samples, 10)
  write.csv(req_genes3, "req_genes3.csv")

} else {
  req_genes3 <- "There are not enough samples to be compared with!"
}

# Top 20 marker genes per treatment -----
if (n_cond > 1) {

  # Get vector of all gene condition heatmaps
  hms_per_cond <- find_func(tempdir, "genes_per_condition_*")

  for (j in seq_along(hms_per_cond)) {

    # Read hm back into memory as df with feats as column
    hm_per_cond <- read.csv(hms_per_cond[j])

    req_genes2 <- get_topn_hm_feats(hm_per_cond, n_cond, 20)
    write.csv(req_genes2, paste0("req_genes2_", j, ".csv"))
  }

} else {
  req_genes2 <- "There are not enough condition to be compared with!"
}

# Save UMAP embedding to csv, save ArchRProject -------------------------------
UMAPHarmony <- getEmbedding(
  ArchRProj = proj, embedding = "UMAP", returnDF = TRUE
)
write.csv(UMAPHarmony, "UMAPHarmony.csv")

saveArchRProject(
  ArchRProj = proj, outputDirectory = paste0(project_name, "_ArchRProject")
)

# Peak calling and motif enrichment for clusters ------------------------------

# Due to mcapply bug, too many threads in peak calling can lead to OOM error;
# decrease theads here, per specification from input parameters
addArchRThreads(threads = num_threads)

proj <- addGroupCoverages(
  ArchRProj = proj,
  groupBy = "Clusters",
  maxCells = 1500,
  force = TRUE
)

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

# Save run metrics in medians.csv ---------------------------------------------
medians <- getCellColData(ArchRProj = proj)
tss <- aggregate(
  medians@listData$TSSEnrichment,
  by = list(medians@listData$Sample),
  FUN = median
)
nfrags <- aggregate(
  medians@listData$nFrags,
  by = list(medians@listData$Sample),
  FUN = median
)
conditions <- aggregate(
  medians@listData$Condition,
  by = list(medians@listData$Sample),
  FUN = max
)
frip <- aggregate(
  medians@listData$FRIP,
  by = list(medians@listData$Sample),
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

# Add motif annotations -------------------------------------------------------
proj <- add_motif_annotations(proj, genome) # from utils

saveArchRProject(
  ArchRProj = proj, outputDirectory = paste0(project_name, "_ArchRProject")
)

# Get marker peaks for clusters, genes, treatments; save as csv ---------------

# Initialize  base data frame -----
peak_data <- data.frame(proj@peakSet@ranges, proj@peakSet@elementMetadata)

# Marker peaks per clusters -----
markers_peaks_c <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  k = 100,
  testMethod = "wilcoxon"
)

# Write significant marker peaks to csv -----
peak_marker_list_c <- getMarkers(
  markers_peaks_c, cutOff = "Pval <= 0.05 & Log2FC >= 0.1"
)
write.csv(
  peak_marker_list_c,
  file = "marker_peaks_per_cluster.csv",
  row.names = FALSE
)

# Merge all peaks with significant marker peaks, write to csv -----
total_peaks_c <- merge(peak_data, peak_marker_list_c, by = c("start", "end"))
write.csv(
  total_peaks_c, file = "complete_peak_list_cluster.csv", row.names = FALSE
)

# Marker peaks per sample -----
markers_peaks_s <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  groupBy = "Sample",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  k = 100,
  testMethod = "wilcoxon"
)

# Write significant marker peaks to csv -----
peak_marker_list_s <- getMarkers(
  markers_peaks_s, cutOff = "Pval <= 0.05 & Log2FC >= 0.1"
)
write.csv(
  peak_marker_list_s,
  file = "marker_peaks_per_sample.csv",
  row.names = FALSE
)

# Merge all peaks with significant marker peaks, write to csv -----
total_peaks_s <- merge(peak_data, peak_marker_list_s, by = c("start", "end"))
write.csv(
  total_peaks_s, file = "complete_peak_list_sample.csv", row.names = FALSE
)

# Marker peaks per treatment -----
if (n_cond > 1) {
  for (i in seq_along(treatment)) {

    marker_peaks_t <- getMarkerFeatures(
      ArchRProj = proj,
      useMatrix = "PeakMatrix",
      groupBy = treatment[i],
      bias = c("TSSEnrichment", "log10(nFrags)"),
      k = 100,
      testMethod = "wilcoxon"
    )

    # Write significant marker peaks to csv -----
    peak_marker_list_t <- getMarkers(
      marker_peaks_t, cutOff = "Pval <= 0.05 & Log2FC >= 0.1"
    )
    write.csv(
      peak_marker_list_t,
      file = paste0("marker_peaks_per_condition-", i, ".csv"),
      row.names = FALSE
    )

    # Merge all peaks with significant marker peaks, write to csv -----
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

# Plot peak, motif heatmaps by cluster, generate heatmaps .pdf ----------------

# Create peak heatmap by cluster -----
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

# Find enriched motifs per cluster, plot as heatmap -----
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

# Save heatmaps to disk as pdf (per cluster: genes, peaks, motifs) -----
print("+++++++++++creating heatmap plots++++++++++++++")

pdf("heatmaps_all.pdf")
for (i in seq_along(heatmaps)) {
  print(heatmaps[[i]])
}
dev.off()

# Per cluster: save enriched motifs, motif heatmap as .rds, .csv's -----
motif_lst <- unique(rownames(enrichMotifs))
split_string <- strsplit(motif_lst, split = "\\(")

req_motifs1 <- gsub("_", "-", extract_nth_ele(split_string))
req_motifs1 <- gsub(" ", "", req_motifs1)

rownames(enrichMotifs) <- req_motifs1

saveRDS(enrichMotifs, "enrichMotifs_clusters.rds")

heatmapEM <- plotEnrichHeatmap(
  enrichMotifs,
  n = 50,
  transpose = FALSE,
  returnMatrix = TRUE,
  cutOff = 2
)

motif_lst <- unique(rownames(heatmapEM))
split_string <- strsplit(motif_lst, split = "\\(")

req_motifs1 <- gsub("_", "-", extract_nth_ele(split_string))
req_motifs1 <- gsub(" ", "", req_motifs1)

rownames(heatmapEM) <- req_motifs1

write.csv(heatmapEM, "motif_per_cluster_hm.csv")

# Save top 5 marker motifs for cluster for ShinyApp  -----

# Read hm back into memory as df with feats as column
hm_motif_per_clust <- read.csv("motif_per_cluster_hm.csv")

req_motifs1 <- get_topn_hm_feats(hm_motif_per_clust, n_clust, 5)
write.csv(req_motifs1, "req_motifs1.csv")

# Create motif SeuratObjects --------------------------------------------------

# Create motif count matrix -----
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

# Get deviation matrix for significant motifs -----
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

# remove 0 deviations per All samples -----
all_zero <- names(which(rowSums(dev_score2) == 0))
dev_score2 <- dev_score2[which(!rownames(dev_score2) %in% c(all_zero)), ]

# convert to dgCmatrix -----
dev_score3 <- Matrix(as.matrix(dev_score2), sparse = TRUE)

# create motif seurat objects -----
seurat_objs_m <- c()
for (run in runs) {

  obj <- build_atlas_seurat_object(
    run_id = run[1],
    matrix = dev_score3,
    metadata = metadata,
    spatial_path = run[6]
  )

  saveRDS(obj, file = paste0(run[1], "_SeuratObjMotif.rds"))
  seurat_objs_m <- c(seurat_objs_m, obj)
}

print("These are the available seurat_objMotifs.")
seurat_objs_m

# Peak calling and motifs for Sample ------------------------------------------
if (n_samples > 1) {

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

  # Add motif annotation, save proj -----
  proj <- add_motif_annotations(proj, genome) # from utils

  saveArchRProject(
    ArchRProj = proj, outputDirectory = paste0(project_name, "_ArchRProject")
  )

  # Get marker peaks and Enriched motifs per sample -----
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

  req_motifs3 <- gsub("_", "-", extract_nth_ele(split_string))
  req_motifs3 <- gsub(" ", "", req_motifs3)

  rownames(enrichMotifs) <- req_motifs3

  saveRDS(enrichMotifs, "enrichMotifs_sample.rds")

  heatmapEM <- plotEnrichHeatmap(
    enrichMotifs, n = 50, transpose = FALSE, returnMatrix = TRUE, cutOff = 2
  )

  motif_lst <- unique(rownames(heatmapEM))
  split_string <- strsplit(motif_lst, split = "\\(")

  req_motifs3 <- gsub("_", "-", extract_nth_ele(split_string))
  req_motifs3 <- gsub(" ", "", req_motifs3)

  rownames(heatmapEM) <- req_motifs3
  write.csv(heatmapEM, "motif_per_sample_hm.csv")

  # Save top 10 marker motifs for sample for ShinyApp  -----

  # Read hm back into memory as df with feats as column
  hm_per_sample <- read.csv("motif_per_sample_hm.csv")

  req_motifs3 <- get_topn_hm_feats(hm_motif_per_clust, n_samples, 10)
  write.csv(req_motifs3, "req_motifs3.csv")

} else {
  req_motifs3 <- "There are not enough samples to be compared with!"
}

# Peak calling and motif enrichment per treatment ----------------------
if (n_cond > 1) {
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

    # Add motif annotation, save proj -----
    proj <- add_motif_annotations(proj, genome) # from utils

    # save ArchR object
    saveArchRProject(
      ArchRProj = proj, outputDirectory = paste0(project_name, "_ArchRProject")
    )
  }

  # Get marker peaks and Enriched motifs per treatment -----
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

    req_motifs2 <- gsub("_", "-", extract_nth_ele(split_string))
    req_motifs2 <- gsub(" ", "", req_motifs2)

    rownames(enrichMotifs) <- req_motifs2
    saveRDS(enrichMotifs, paste0("enrichMotifs_condition_", i, ".rds"))

    heatmapEM <- plotEnrichHeatmap(
      enrichMotifs, n = 50, transpose = FALSE, returnMatrix = TRUE, cutOff = 2
    )

    motif_lst <- unique(rownames(heatmapEM))
    split_string <- strsplit(motif_lst, split = "\\(")

    req_motifs2 <- gsub("_", "-", extract_nth_ele(split_string))
    req_motifs2 <- gsub(" ", "", req_motifs2)

    rownames(heatmapEM) <- req_motifs2
    write.csv(heatmapEM, paste0("motif_per_condition_", i, "_hm.csv"))
  }

  # Save top 5 marker motifs for treatment for ShinyApp  -----
  motifs_per_cond_hm <- find_func(tempdir, "motif_per_condition_*")

  for (j in seq_along(motifs_per_cond_hm)) {

    # Read hm back into memory as df with feats as column
    hm_per_cond <- read.csv(motifs_per_cond_hm[j])

    req_motifs2 <- get_topn_hm_feats(hm_per_cond, n_cond, 5)
    write.csv(req_motifs2, paste0("req_motifs2_", j, ".csv"))
  }
} else {
  req_motifs2 <- "There are not enough conditions to be compared with!"
}

# Volcano plots for motifs -----------------------------------------------------
if (n_cond > 1) {
  for (j in seq_along(treatment)) {

    # Get motif markers for all clusters together -----
    marker_motifs_df <- get_marker_df( # from archr.R
      proj = proj,
      group_by = treatment[j],
      matrix = "MotifMatrix",
      seq_names = "z",
      max_cells = 5000,
      test_method = "wilcoxon"
    )

    # Create a data from of marker genes for clusters for which no condition is
    # >90% of all cells
    req_clusters <- get_required_clusters(proj, treatment[j]) # from archr.R
    marker_motifs_by_cluster_df <- get_marker_df_clusters(
      proj, req_clusters, treatment[j]
    )

    # Merge and cleanup data -----
    conditions <- sort(unique(proj@cellColData[treatment[j]][, 1]))
    for (cond in conditions) {

      volcano_table <- get_volcano_table( # from archr.R
        marker_genes_df, marker_genes_by_cluster_df, cond, "gene", gene_matrix
      )

      write.table(
        volcano_table,
        paste0("volcanoMarkers_motifs_", j, "_", cond, ".txt"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
      )
      print(
        paste0("writing volcanoMarkers_motifs_", j, "_", cond, ".txt is done!")
      )

      features_m <- unique(volcano_table$cluster)
      others <- paste(conditions[conditions != cond], collapse = "|")
      volcano_plots_m <- list()
      for (i in seq_along(features_m)) {
        volcano_plots_m[[i]] <- scvolcano(
          volcano_table,  cond, others, features_m[[i]]
        )
      }

      pdf(paste0("volcano_plots_motifs_", j, "_", cond, ".pdf"))
      for (plot in volcano_plots_m) {
        print(plot)
      }
      dev.off()
    }
  }
} else {
  print("There are not enough conditions to be compared with!")
}

# Save bigwig files for clustersXtreatment -------------------------------------
req_conditions <- c("Clusters", treatment)

for (i in req_conditions) {
  bws <- getGroupBW(ArchRProj = proj, groupBy = i)
}

saveArchRProject(
  ArchRProj = proj, outputDirectory = paste0(project_name, "_ArchRProject")
)

# Motif Logo ------------------------------------------------------------------
pwms <- getPeakAnnotation(proj, "Motif")$motifs

pw_matrix_to_prob_matrix <- function(x) {
  if (class(x) != "PWMatrix") stop("x must be a TFBSTools::PWMatrix object")
  m <- (exp(as(x, "matrix"))) * TFBSTools::bg(x) / sum(TFBSTools::bg(x))
  m <- t(t(m) / colSums(m))
  m
}

prob_matrices <- lapply(pwms, pw_matrix_to_prob_matrix)
lapply(prob_matrices, colSums) %>% range

saveRDS(prob_matrices, "seqlogo.rds")

# combined SeuratObjs ---------------------------------------------------------

# Rename SeuratObj cells with run_id to ensure unique when combining -----
all <- rename_cells(seurat_objs)
all_m <- rename_cells(seurat_objs_m)

samples <- find_sample_names(all)

# extract image coordinates as -(imagecols) | imagerow -----
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

# create Shiny app files ------------------------------------------------------
print("Shiny App starting...")

# Create ShinyApp/ files for genes (sc1) and motifs (sc2) -----
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

# Edit def files -------------------------------------------------------
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

if (n_cond > 1) {
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

# Edit conf files -------------------------------------------------------------
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

# Select correct ui/server.R files ---------------------------------------------
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

  if (n_samples == 1) {
    file.copy("/root/uiserver50by50/ui_3.R", "/root/ui.R", overwrite = TRUE)
    file.copy(
      "/root/uiserver50by50/server_3.R", "/root/server.R", overwrite = TRUE
    )
    unlink("/root/uiserver50by50", recursive = TRUE) # will delete directory
  } else if (n_cond <= 1) {
    file.copy("/root/uiserver50by50/ui_2.R", "/root/ui.R", overwrite = TRUE)
    file.copy(
      "/root/uiserver50by50/server_2.R", "/root/server.R", overwrite = TRUE
    )
    unlink("/root/uiserver50by50", recursive = TRUE) # will delete directory
  } else if (n_cond > 1 && length(unique(proj$condition_1)) > 2 ||
      n_cond > 1 && length(unique(proj$condition_2)) > 2
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

  if (n_samples == 1) {
    file.copy("/root/uiserver96by96/ui_3.R", "/root/ui.R", overwrite = TRUE)
    file.copy(
      "/root/uiserver96by96/server_3.R", "/root/server.R", overwrite = TRUE
    )
    unlink("/root/uiserver96by96", recursive = TRUE) # will delete directory
  } else if (n_cond <= 1) {
    file.copy("/root/uiserver96by96/ui_2.R", "/root/ui.R", overwrite = TRUE)
    file.copy(
      "/root/uiserver96by96/server_2.R", "/root/server.R", overwrite = TRUE
    )
    unlink("/root/uiserver96by96", recursive = TRUE) # will delete directory
  } else if (
    n_cond > 1 && length(unique(proj$condition_1)) > 2 ||
      n_cond > 1 && length(unique(proj$condition_2)) > 2
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

# Copy everything in shinyApp/ subfolder to the /root/ -----
rawPath <- "/root/shinyApp/"
dataPath <- "/root/"

dataFiles <- dir(rawPath, "*.rds$", ignore.case = TRUE, all.files = TRUE)
file.copy(file.path(rawPath, dataFiles), dataPath, overwrite = TRUE)

dataFiles <- dir(rawPath, "*.h5$", ignore.case = TRUE, all.files = TRUE)
file.copy(file.path(rawPath, dataFiles), dataPath, overwrite = TRUE)
