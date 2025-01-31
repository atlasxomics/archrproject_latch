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

inputs <- c() # Inputs for ArrowFiles (run_id : fragment_file path)
for (run in runs) {
  inputs[run[1]] <- run[3]
}
inputs

out_dir <- paste0(project_name, "_ArchRProject")

# Save input metrics to csv ----
metrics_to_csv(as.list(args[1:13])) # from utils.R

# Set genome size for peak calling ----
genome_sizes <- list("hg38" = 3.3e+09, "mm10" = 3.0e+09, "rnor6" = 2.9e+09)
genome_size <- genome_sizes[[genome]]

tempdir <- "/root"
archrproj_dir <- paste0(project_name, "_ArchRProject")

# Create ArchRProject ---------------------------------------------------------
addArchRThreads(threads = 50)

proj <- create_archrproject( # from archr.R
  inputs, genome, min_tss, min_frags, tile_size, out_dir
)

# Add Conditions, SampleName to CellColData ----
for (run in runs) {
  proj$Condition[proj$Sample == run[1]] <- run[4]
  proj$SampleName[proj$Sample == run[1]] <- run[2]
}

# Filter to keep on-tissue only ----
all_ontissue <- c()
for (run in runs) {
  positions <- read.csv(run[5], header = FALSE)
  positions$V1 <- paste0(run[1], "#", positions$V1, "-1")
  on_tissue <- positions$V1 [which(positions$V2 == 1)]
  all_ontissue <- c(all_ontissue, on_tissue)
}

proj <- proj[proj$cellNames %in% all_ontissue]

# Parse conditions into 'treatments', add as columns to CellColData ----
conds <- strsplit(proj$Condition, split = "\\s|-")

for (i in seq_along(conds[[1]])) {
  proj <- ArchR::addCellColData(
    proj,
    data = extract_nth_ele(conds, i), # from utils.R
    name = paste0("condition_", i),
    cells = proj$cellNames,
    force = TRUE
  )
}
treatment <- names(ArchR::getCellColData(proj))[
  grep("condition_", names(ArchR::getCellColData(proj)))
]

print(paste("Treatments:", treatment))

# Set global values for project data ----
n_samples <- length(unique(proj$Sample))
n_cond <- length(unique(proj$Condition))
n_cells <- length(proj$cellNames)

saveArchRProject(ArchRProj = proj, outputDirectory = archrproj_dir)

# Dimensionality reduction and clustering ----
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

# Batch correction with Harmony if more than one run --
if (n_samples > 1) {
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

# Plot UMAP embedding colored by Cluster, Sample, Treatment --
save_umap(proj, c("Clusters", "Sample", treatment))

# Save UMAP embedding to csv --
umap_harmony <- getEmbedding(proj)
write.csv(umap_harmony, "UMAPHarmony.csv")

saveArchRProject(ArchRProj = proj, outputDirectory = archrproj_dir)

# Gene Expression Analysis ----------------------------------------------------

# Create SeuratObjects for gene matrix ----
# Extract metadata for Seurat object --
metadata <- getCellColData(ArchRProj = proj)

# Set metadata rownames to barcodes
rownames(metadata) <- str_split_fixed(
  str_split_fixed(
    row.names(metadata),
    "#",
    2
  )[, 2],
  "-",
  2
)[, 1]

# Create col for log10 of fragment counts
metadata["log10_nFrags"] <- log(metadata$nFrags)

# Extract gene matrix for SeuratObject --
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

# Create and save SeuratObjects --
print("Creating SeuratObjects...")

seurat_objs <- c()
for (run in runs) {

  obj <- build_atlas_seurat_object( # from seurat.R
    run_id = run[1],
    matrix = matrix,
    metadata = metadata,
    spatial_path = run[6]
  )

  saveRDS(obj, file = paste0(run[1], "_SeuratObj.rds"))
  seurat_objs <- c(seurat_objs, obj)
}

print("Available SeuratObjects:")
seurat_objs

# SeuratObj plots ----
# Plot clusters on top of spatial coordinates --
print("Creating spatial cluster plots...")
save_spatial_cluster_plots(seurat_objs)

# Create QC plots for TSSEnrichment, nFrags, log10_nFrags --
print("Creating QC plots...")
save_qc_plots(seurat_objs, c("TSSEnrichment", "nFrags", "log10_nFrags"))

saveArchRProject(ArchRProj = proj, outputDirectory = archrproj_dir)

# Identify marker genes ----
# Marker genes per cluster, save marker gene rds, csv, heatmap.csv --
n_clust <- length(unique(proj$Clusters))

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

# Save top 5 marker genes by cluster for ShinyApp (req_genes1.csv)
# Read hm back into memory as df with feats as column
hm_per_cluster <- read.csv("genes_per_cluster_hm.csv")
req_genes1 <- get_topn_hm_feats(hm_per_cluster, n_clust, 5) # archr.R
write.csv(req_genes1, "req_genes1.csv")

# Marker genes per sample, save marker gene rds, csv, heatmap.csv --
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

  # Save top 10 marker genes per sample for ShinyApp (req_genes3.csv)
  # Read hm back into memory as df with feats as column
  hm_per_sample <- read.csv("genes_per_sample_hm.csv")
  req_genes3 <- get_topn_hm_feats(hm_per_sample, n_samples, 10)
  write.csv(req_genes3, "req_genes3.csv")

} else {
  req_genes3 <- "There are not enough samples to be compared with!"
}

# Marker genes per treatment, save marker gene rds, csv, heatmap.csv --
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

    # Save top 20 marker genes per treatment for ShinyApp (req_genes3.csv)
    # Get vector of all gene condition heatmap csv's
    hms_per_cond <- find_func(tempdir, "genes_per_condition_*")

    for (j in seq_along(hms_per_cond)) {

      # Read hm back into memory as df with feats as column
      hm_per_cond <- read.csv(hms_per_cond[j])

      req_genes2 <- get_topn_hm_feats(hm_per_cond, n_cond, 20)
      write.csv(req_genes2, paste0("req_genes2_", j, ".csv"))
    }
  }
} else {
  req_genes2 <- "There are not enough condition to be compared with!"
}

# Volcano plots for genes ----
if (n_cond > 1) {
  for (j in seq_along(treatment)) {

    # Get gene markers df for all clusters together --
    marker_genes_df <- get_marker_df(
      proj = proj,
      group_by = treatment[j],
      matrix = "GeneScoreMatrix",
      seq_names = NULL,
      max_cells = n_cells,  # Equals total cells in project
      test_method = "ttest"
    )

    # Create a merged marker genes df for clusters for which no condition is
    # >90% of all cells --
    req_clusters <- get_required_clusters(proj, treatment[j])
    marker_genes_by_cluster_df <- get_marker_df_clusters(
      proj, req_clusters, treatment[j]
    )

    # Per condition, merge dfs and cleanup data --
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
      print(paste0("volcanoMarkers_genes_", j, "_", cond, ".txt is done!"))

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

saveArchRProject(ArchRProj = proj, outputDirectory = archrproj_dir)

# Peak and motif analysis -----------------------------------------------------

# Due to mcapply bug, too many threads in peak calling can lead to OOM error;
# decrease theads here, per specification from input parameters ----
addArchRThreads(threads = num_threads)

# Peak calling and motif enrichment for clusters ----
proj <- get_annotated_peaks(proj, "Clusters", genome_size, genome)

saveArchRProject(ArchRProj = proj, archrproj_dir)

# Save run metrics in medians.csv ----
medians <- get_proj_medians(proj)
write.csv(medians, file = "medians.csv", row.names = FALSE)

# Get marker peaks for clusters, samples, treatments; save as csv ----

# Initialize base data frame, set significance cutoff --
peak_data <- data.frame(proj@peakSet@ranges, proj@peakSet@elementMetadata)
cut_off <- "Pval <= 0.05 & Log2FC >= 0.1"

# Marker peaks per clusters --
marker_peaks_c <- get_marker_peaks(proj, "Clusters", peak_data, cut_off)

write.csv(
  marker_peaks_c$marker_peak_list,
  "marker_peaks_per_cluster.csv",
  row.names = FALSE
)

write.csv(
  marker_peaks_c$total_peaks,
  "complete_peak_list_cluster.csv",
  row.names = FALSE
)

# Marker peaks per sample --
marker_peaks_s <- get_marker_peaks(proj, "Sample", peak_data, cut_off)

write.csv(
  marker_peaks_s$marker_peak_list,
  "marker_peaks_per_sample.csv",
  row.names = FALSE
)

write.csv(
  marker_peaks_s$total_peaks,
  "complete_peak_list_sample.csv",
  row.names = FALSE
)

# Marker peaks per treatment --
if (n_cond > 1) {

  for (i in seq_along(treatment)) {

    marker_peaks_t <- get_marker_peaks(proj, treatment[i], peak_data, cut_off)

    write.csv(
      marker_peaks_t$marker_peak_list,
      file = paste0("marker_peaks_per_condition-", i, ".csv"),
      row.names = FALSE
    )

    write.csv(
      marker_peaks_t$total_peaks,
      file = paste0("complete_peak_list_condition-", i, ".csv"),
      row.names = FALSE

    )
  }
}

# Get enriched motif data for clusters, write to disk ----
enriched_motifs_c <- get_enriched_motifs(
  proj, marker_peaks_c$marker_peaks, cut_off
)

write.csv(enriched_motifs_c$enrich_df, "enrichedMotifs_cluster.csv")
saveRDS(enriched_motifs_c$enrich_motifs, "enrichMotifs_clusters.rds")
write.csv(enriched_motifs_c$heatmap_em, "motif_per_cluster_hm.csv")

# Save top 5 marker motifs for cluster for ShinyApp  --
# Read hm back into memory as df with feats as column
hm_motif_per_clust <- read.csv("motif_per_cluster_hm.csv")

req_motifs1 <- get_topn_hm_feats(hm_motif_per_clust, n_clust, 5)
write.csv(req_motifs1, "req_motifs1.csv")

# Create motif SeuratObjects ----
# Create motif count matrix --
proj <- addBgdPeaks(proj, force = TRUE)

proj <- addDeviationsMatrix(
  ArchRProj = proj,
  peakAnnotation = "Motif",
  force = TRUE
)

markers_motifs <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "MotifMatrix",
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  useSeqnames = "z"
)

# Get deviation matrix for significant motifs --
marker_motifs_list <- getMarkers(
  markers_motifs, cutOff = "FDR < 0.9 & MeanDiff >= 0"
)

motifs <- list()
for (i in seq_len(length(marker_motifs_list))) {
  if (length(marker_motifs_list[[i]]$name) > 1) {
    motifs <- c(motifs, marker_motifs_list[[i]]$name)
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

# Remove 0 deviations per All samples --
all_zero <- names(which(rowSums(dev_score2) == 0))
dev_score2 <- dev_score2[which(!rownames(dev_score2) %in% c(all_zero)), ]

# Convert to dgCmatrix --
dev_score3 <- Matrix(as.matrix(dev_score2), sparse = TRUE)

# Create motif seurat objects --
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

print("Available seurat_objMotifs:")
seurat_objs_m

# Peak calling and motifs for Sample ----
if (n_samples > 1) {

  proj <- get_annotated_peaks(proj, "Sample", genome_size, genome)

  saveArchRProject(ArchRProj = proj, outputDirectory = archrproj_dir)

  # Repeat getMarkerPeaks for new peak-set, enriched motifs per sample --
  sample_marker_peaks <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix",
    groupBy = "Sample",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    k = 100,
    testMethod = "wilcoxon"
  )

  enriched_motifs_s <- get_enriched_motifs(
    proj, sample_marker_peaks, cut_off
  )

  write.csv(enriched_motifs_s$enrich_df, "enrichedMotifs_sample.csv")
  saveRDS(enriched_motifs_s$enrich_motifs, "enrichMotifs_sample.rds")
  write.csv(enriched_motifs_s$heatmap_em, "motif_per_sample_hm.csv")

  # Save top 10 marker motifs for sample for ShinyApp  --
  # Read hm back into memory as df with feats as column
  hm_per_sample <- read.csv("motif_per_sample_hm.csv")
  req_motifs3 <- get_topn_hm_feats(hm_motif_per_clust, n_samples, 10)
  write.csv(req_motifs3, "req_motifs3.csv")

} else {
  req_motifs3 <- "There are not enough samples to be compared with!"
}

# Peak calling and motif enrichment per treatment ----
if (n_cond > 1) {

  for (i in seq_along(treatment)) {

    proj <- get_annotated_peaks(proj, treatment[i], genome_size, genome)

    saveArchRProject(ArchRProj = proj, outputDirectory = archrproj_dir)

    # Get marker peaks and Enriched motifs per treatment --
    treatment_marker_peaks <- getMarkerFeatures(
      ArchRProj = proj,
      useMatrix = "PeakMatrix",
      groupBy = treatment[i],
      bias = c("TSSEnrichment", "log10(nFrags)"),
      k = 100,
      testMethod = "wilcoxon"
    )
    enriched_motifs_t <- get_enriched_motifs(
      proj, treatment_marker_peaks, cut_off
    )

    write.csv(
      enriched_motifs_t$enrich_df,
      paste0("enrichedMotifs_condition_", i, ".csv")
    )
    saveRDS(
      enriched_motifs_t$enrich_motifs,
      paste0("enrichMotifs_condition_", i, ".rds")
    )
    write.csv(
      enriched_motifs_t$heatmap_em,
      paste0("motif_per_condition_", i, "_hm.csv")
    )
  }

  # Save top 5 marker motifs for treatment for ShinyApp  --
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

# Plot gene, peak, motif heatmaps by cluster, generate heatmaps .pdf ----
# Initiate heatmaps list and plot marker gene per cluster heatmap --
heatmaps <- list()

# Recompute gene heatmap for plotting (transpose, plotLog2FC different) --
heatmap_gs_plotting <- plotMarkerHeatmap(
  seMarker = cluster_marker_genes$markers_gs,
  cutOff = cut_off,
  transpose = TRUE
)

# Save for plotting with peaks and motifs --
gene_hm <- ComplexHeatmap::draw(
  heatmap_gs_plotting,
  heatmap_legend_side = "bot",
  annotation_legend_side = "bot",
  column_title = paste0("Marker genes (", cut_off, ")"),
  column_title_gp = gpar(fontsize = 12)
)
heatmaps[[1]] <- gene_hm

# Create peak heatmap by cluster --
heatmap_peaks <- plotMarkerHeatmap(
  seMarker = marker_peaks_c$marker_peaks,
  cutOff = cut_off,
  transpose = TRUE
)

peak_hm <- ComplexHeatmap::draw(
  heatmap_peaks,
  heatmap_legend_side = "bot",
  annotation_legend_side = "bot",
  column_title = paste0("Marker peaks (", cut_off, ")"),
  column_title_gp = gpar(fontsize = 12)
)

heatmaps[[2]] <- peak_hm

# Create motif heatmap by cluster --
heatmap_plot <- plotEnrichHeatmap(
  enriched_motifs_c$enrich_motifs,
  transpose = TRUE,
  n = 50,
  cutOff = 2
)

heatmap_motifs <- ComplexHeatmap::draw(
  heatmap_plot,
  heatmap_legend_side = "bot",
  column_title = paste0("Marker motifs (", cut_off, ")"),
  column_title_gp = gpar(fontsize = 12)
)

heatmaps[[3]] <- heatmap_motifs

# Save heatmaps to disk as pdf (per cluster: genes, peaks, motifs) --
print("Saving heatmap plots...")

pdf("heatmaps_all.pdf")
for (i in seq_along(heatmaps)) {
  print(heatmaps[[i]])
}
dev.off()

# Volcano plots for motifs ----
if (n_cond > 1) {
  for (j in seq_along(treatment)) {

    # Get motif markers for all clusters together --
    marker_motifs_df <- get_marker_df( # from archr.R
      proj = proj,
      group_by = treatment[j],
      matrix = "MotifMatrix",
      seq_names = "z",
      max_cells = 5000,
      test_method = "wilcoxon"
    )

    # Create a data from of marker genes for clusters for which no condition is
    # >90% of all cells --
    req_clusters <- get_required_clusters(proj, treatment[j]) # from archr.R
    marker_motifs_by_cluster_df <- get_marker_df_clusters(
      proj, req_clusters, treatment[j]
    )

    # Merge and cleanup data --
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

# Save bigwig files for clustersXtreatment ----
req_conditions <- c("Clusters", treatment)

for (i in req_conditions) {
  bws <- getGroupBW(ArchRProj = proj, groupBy = i)
}

saveArchRProject(ArchRProj = proj, outputDirectory = archrproj_dir)

# Motif Logo ----
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

# Create Shiny app files ------------------------------------------------------
print("Shiny App starting...")

# Create combined SeuratObjs ----
# Rename SeuratObj cells with run_id to ensure unique when combining --
all <- rename_cells(seurat_objs)
all_m <- rename_cells(seurat_objs_m)

samples <- find_sample_names(all)

# extract image coordinates as -(imagecols) | imagerow --
spatial <- lapply(all, function(x) {
  df <- as.data.frame(x@images[[1]]@coordinates[, c(5, 4)])
  colnames(df) <- paste0("Spatial_", 1:2)
  df$Spatial_2 <- -df$Spatial_2
  df
})

combined <- combine_objs(all, umap_harmony, samples, spatial, project_name)
combined_m <- combine_objs(all_m, umap_harmony, samples, spatial, project_name)

saveRDS(combined, "combined.rds", compress = FALSE)
saveRDS(combined_m, "combined_m.rds", compress = FALSE)

# Create ShinyApp/ files for genes (sc1) and motifs (sc2) ----
sc_conf1 <- createConfig(combined)
makeShinyFiles(
  combined,
  sc_conf1,
  gex.assay = "scATAC",
  gex.slot = "data",
  gene.mapping = TRUE,
  shiny.prefix = "sc1",
  default.gene1 = "Tiam1",
  default.gene2 = "Ccbe1",
  default.dimred = c("UMAP_1", "UMAP_2")
)

sc_conf2 <- createConfig(combined_m)
makeShinyFiles(
  combined_m,
  sc_conf2,
  gex.assay = "scATAC",
  gex.slot = "counts",
  gene.mapping = TRUE,
  shiny.prefix = "sc2",
  default.gene1 = "RFX3-1018",
  default.gene2 = "NEUROG2-1580",
  default.dimred = c("UMAP_1", "UMAP_2")
)

# Set Shiny App defaults ----
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

xlim <- lapply(spatial, function(x) {
  xlim <- c(min(x[, 1]), max(x[, 1]))
  xlim
})
ylim <- lapply(spatial, function(y) {
  ylim <- c(min(y[, 2]), max(y[, 2]))
  ylim
})

# Edit ShinyApp files ----
# scXdef.rds --
sc1def <- edit_shiny_def(
  "sc1", samples, combined, req_genes1, req_genes3, xlim, ylim
)
sc2def <- edit_shiny_def(
  "sc2", samples, combined_m, req_motifs1, req_motifs3, xlim, ylim
)

if (n_cond > 1) {
  sc1def <- add_def_treatments(sc1def, treatment, "genes", combined)
  sc2def <- add_def_treatments(sc2def, treatment, "motifs", combined_m)
} else {
  print("There are not enough conditions to add to sc1def.rds or sc2def.rds")
}

saveRDS(sc1def, "/root/shinyApp/sc1def.rds")
saveRDS(sc2def, "/root/shinyApp/sc2def.rds")

# scXconf.rds --
sc1conf <- edit_shiny_conf("sc1", treatment)
sc2conf <- edit_shiny_conf("sc2", treatment)

saveRDS(sc1conf, "/root/shinyApp/sc1conf.rds")
saveRDS(sc2conf, "/root/shinyApp/sc2conf.rds")

# Select correct ui/server.R files ----
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
  print("Using ui/server from 50by50 folder...")
  folder <- "/root/uiserver50by50"
} else {
  print("Using ui/server from 96by96 folder...")
  folder <- "/root/uiserver96by96"
}

files <- determine_files_to_copy(n_samples, n_cond, proj)
copy_ui_server_files(folder, files$ui_file, files$server_file)

# Copy everything in shinyApp/ subfolder to the /root/ -----
raw_path <- "/root/shinyApp/"
data_path <- "/root/"

for (ext in c("*.rds$", "*.h5$")) {
  data_files <- dir(raw_path, ext, ignore.case = TRUE, all.files = TRUE)
  file.copy(file.path(raw_path, data_files), data_path)
}
