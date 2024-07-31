#' Generic rountines

library("dplyr")
library("Seurat")

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
  write.csv(metric_list, file = "metadata.csv", row.names = FALSE)
}