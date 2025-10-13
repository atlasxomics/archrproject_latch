#' Generic rountines

library("dplyr")
library("Seurat")


build_atlas_seurat_object <- function(
  run_id,
  matrix,
  metadata,
  spatial_path
) {
  # Prepare and combine gene matrix, metadata, and image for seurat object
  # for runs within a project.

  for (i in seq_along(treatment)) {
    def[[paste0("meta", (2 + i))]] <- treatment[i]

    def[[treatment[i]]] <- read.csv(
      paste0("req_", feature, "2_", i, ".csv")
    )$x

    def[[paste0(treatment[i], "_1")]] <- sort(
      unique(combined@meta.data[[treatment[i]]])
    )[1]

    def[[paste0(treatment[i], "_2")]] <- sort(
      unique(combined@meta.data[[treatment[i]]])
    )[2]
  }

  return(def)
}

copy_ui_server_files <- function(folder, ui_file, server_file) {

  file.copy(file.path(folder, ui_file), "/root/ui.R")
  file.copy(file.path(folder, server_file), "/root/server.R")
}

determine_files_to_copy <- function(n_samples, n_cond, proj) {

  if (n_samples == 1) {
    return(list(ui_file = "ui_3.R", server_file = "server_3.R"))

  } else if (n_cond <= 1) {
    return(list(ui_file = "ui_2.R", server_file = "server_2.R"))

  } else if (n_cond > 1 && (
    length(unique(proj$condition_1)) > 2 || length(unique(proj$condition_2)) > 2
  )) {
    return(list(ui_file = "ui_4.R", server_file = "server_4.R"))

  } else {
    return(list(ui_file = "ui.R", server_file = "server.R"))
  }
}

scvolcano <- function(
  inpMarkers, condition1, condition2, feature = "All", fc_col = "Log2FC"
) {

  # Subset by cluster
  ggData <- inpMarkers[inpMarkers$cluster == feature, ]
  minfdr <- 0.09
  minfdr1 <- 10^-(1 / 6 * (-log10(min(ggData$p_val_adj))))

  # Add Significance column
  ggData$Significance <- ifelse(
    ggData$p_val_adj < minfdr,
    ifelse(
      ggData[[fc_col]] > 0.0,
      condition1,
      condition2
    ),
    "Not siginficant"
  )

  ggData$Significance <- factor(
    ggData$Significance,
    levels = c(condition1, condition2, "Not siginficant")
  )

  # Avoid log10(0)
  ggData[ggData$p_val_adj < 1e-300, "p_val_adj"] <- 1e-300
  ggData$log10fdr <- -log10(ggData$p_val_adj)

  # Make volcano plot
  ggOut <- ggplot(ggData, aes(x = .data[[fc_col]], y = log10fdr)) +
    geom_point() +
    sctheme() +
    ylab("-log10(FDR)") +
    geom_point(aes(color = Significance)) +
    scale_color_manual(values = c("#F8766D", "#619CFF", "gray")) +
    geom_text_repel(
      data = subset(ggData, p_val_adj < minfdr1), aes(label = gene)
    ) +
    ggtitle(paste("Markers:", feature)) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20),
      legend.text = element_text(size = 15),
      legend.title = element_text(size = 18)
    )

  return(ggOut)
}

  # Add coordinate limits --
  def$limits <- list()
  for (i in seq_along(samples)) {
    def[["limits"]][[samples[i]]] <- c(
      min(xlim[[i]]), max(xlim[[i]]), min(ylim[[i]]), max(ylim[[i]])
    )
  }

  # Set defaults for meta --
  def$meta1 <- "Clusters"
  def$meta2 <- "Sample"
  def$meta3 <- "SampleName"

  # Add default feats for Clusters, Sample --
  def$Clusters <- req_feats1
  def$Sample <- req_feats3

  # Set reductions values --
  def$dimred[3] <- paste0(names(combined[3]@reductions)[1], "_1")
  def$dimred[4] <- paste0(names(combined[3]@reductions)[1], "_2")
  def$dimred[5] <- paste0(names(combined[3]@reductions)[2], "_1")
  def$dimred[6] <- paste0(names(combined[3]@reductions)[2], "_2")

  return(def)
}

find_func <- function(tempdir, pattern) {
  list.files(
    path = tempdir,     #' Replace tempdir with the directory you want, followed
    pattern = pattern,  #' by 0 or more characters, then ".csv", and then
    full.names = TRUE,  #' nothing else ($) include the directory in the result.
    recursive = TRUE
  )
}

find_sample_names <- function(seurat_lst) {
  #' Extract list of sample names from list of SeuratObjs.
  sapply(
    seq_along(seurat_lst), function(i) {
      unique(seurat_lst[[i]]@meta.data$Sample)
    }
  )
}

extract_nth_ele <- function(lst, n = 1) {
  #' Extract nth element for each item in an array; return as a vector.
  sapply(lst, `[`, n)
}

metrics_to_csv <- function(metric_list) {
  #' Write inpute metrics to csv
  names(metric_list) <- c(
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
  write.csv(metric_list, file = "exe_metadata.csv", row.names = FALSE)
}