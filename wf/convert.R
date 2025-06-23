#' Rountines for manipulation of and plotting with Seurat

library("BPCells")
library("ggplot2")
library("ggrepel")
library("gridExtra")
library("harmony")
library("Matrix")
library("patchwork")
library("purrr")
library("Seurat")
library("SeuratDisk")

source("/root/wf/utils.R")


convert_spatial <- function(seurat_obj) {

  # Try to extract spatial coordinates from the 'images' slot
  if (length(seurat_obj@images) > 0) {
    message("Attempting spatial coords from image slot...")

    # Check if the image is Visium format
    if (inherits(seurat_obj@images[[1]], "VisiumV1")) {
      # Extract coordinates
      coords <- data.frame(seurat_obj@images[[1]]@coordinates[, c(5, 4)])
      coords[[2]] <- -coords[[2]]
      colnames(coords) <- c("spatial_1", "spatial_2")

      # Create and add the spatial reduction
      spatial_embedding <- Seurat::CreateDimReducObject(
        embeddings = as.matrix(coords),
        key = "spatial",
        assay = Seurat::DefaultAssay(seurat_obj)
      )
      seurat_obj@reductions[["spatial"]] <- spatial_embedding

      message("Successfully created spatial coordinates from image slot.")
      return(seurat_obj)
    } else {
      message("Image slots must be class VisiumV1; trying alternative method.")
      # Continue to next method (fall through)
    }
  }

  # Try to extract spatial coordinates from reductions if Sample metadata exists
  if ("Sample" %in% colnames(seurat_obj@meta.data)) {

    message("Attempting spatial coords from reductions...")

    samples <- unique(seurat_obj@meta.data[["Sample"]])
    reduction_names <- names(seurat_obj@reductions)

    # Check if all samples have corresponding reductions
    if (all(samples %in% reduction_names)) {
      # Extract and combine spatial coordinates for each sample
      reductions <- list()
      for (sample in samples) {
        mat <- seurat_obj@reductions[[sample]]@cell.embeddings
        sample_cells <- grepl(paste0("^", sample), rownames(mat))
        if (sum(sample_cells) > 0) {
          sample_mat <- mat[sample_cells, ]
          reductions[[sample]] <- sample_mat
        } else {
          message(paste0("No matching cells found for sample: ", sample))
        }
      }

      # Only proceed if we have data to work with
      if (length(reductions) > 0) {
        spatial <- do.call(rbind, reductions)
        spatial_embedding <- Seurat::CreateDimReducObject(
          embeddings = spatial,
          key = "spatial",
          assay = Seurat::DefaultAssay(seurat_obj)
        )
        seurat_obj@reductions[["spatial"]] <- spatial_embedding
        message("Successfully reformatted spatial SeuratObject coordinates.")
        return(seurat_obj)
      }
    } else {
      message("Mismatch between sample names and
              embeddings; skipping conversion.")
      # Continue to final return
    }
  }

  # If we reach here, we couldn't find spatial coordinates
  message("Unable to find spatial coordinates in SeuratObject;
          aborting spatial conversion")
  return(seurat_obj)
}

seurat_to_h5ad <- function(seuratobj, spatial, prefix) {

  if (spatial) {
    # Reformat spatial info for squidpy compatibility
    seuratobj <- convert_spatial(seuratobj)  # From utils.R
  }

  SaveH5Seurat(seuratobj, filename = "temp.h5Seurat", overwrite = TRUE)
  Convert("temp.h5Seurat", dest = "h5ad", overwrite = TRUE)
  file.rename("temp.h5ad", paste0(prefix, "_converted.h5ad"))

  message("Seurat object converted to h5ad file.")
}