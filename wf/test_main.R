library(ArchR)
library(dplyr)
library(purrr)
library(Seurat)
library(stringr)

source("utils.R")
source("../getDeviation_ArchR.R")

setwd("/Users/jamesm/latch/archrproject/wf")

proj <- loadArchRProject("demo_ArchRProject/")
seurat_objs <- c(
  readRDS("default_SeuratObj.rds"),
  readRDS("test_D01116_NG01973_SeuratObjMotif.rds"),
  readRDS("test_D01117_NG01974_SeuratObjMotif.rds")
)

all <-  list()
for (i in seq_along(seurat_objs)) {
  all[[i]] <- seurat_objs[[i]]
  all[[i]] <- Seurat::RenameCells(
    all[[i]],
    new.names = paste0(
      unique(all[[i]]@meta.data$Sample),
      "#",
      colnames(all[[i]]),
      "-1"
    )
  )
}

UMAPHarmony <- getEmbedding(
  ArchRProj = proj, embedding = "UMAP", returnDF = TRUE
)

find_samples_name <- function(seurat_lst) {
  # Extract list of sample names from list of SeuratObjs.
  sapply(
    seq_along(seurat_lst), function(i) {
      unique(seurat_lst[[i]]@meta.data$Sample)
    }
  )
}

main_func <- function(seurat_lst, umap_embedding) {
  
  # ----------- Prepare data -----------
  
  samples <- find_samples_name(seurat_lst)
  
  # Ensure that cell names are unique, if not append run_id
  seurat_lst <- lapply(seq_along(seurat_lst), function(i){
    if (!str_ends(Cells(seurat_lst[[i]])[[1]], "-1")) {
      new_cellnames <- paste0(samples[i], "#", Cells(seurat_lst[[i]]), "-1")
      seurat_lst[[i]] <- RenameCells(seurat_lst[[i]], new.names = new_cellnames)
    }
  })
  
  # ----------- Prepare feature matrices -----------
  # Remove features with 0 counts from objs
  filtered <- lapply(seurat_lst, function(x) {
    to_keep <- rowSums(x@assays[[1]]@counts) > 0
    x <- subset(x, features = rownames(x)[to_keep])
  })
  
  # Convert list of SeuratObjs to list of Assay counts as dataframes.
  filtered_dfs <- lapply(filtered, function(x) {
    df <- as.data.frame(x@assays[[1]]@counts)
    colnames(df) <- Seurat::Cells(x)
    df$region <- rownames(df)
    df
  })

  # ----------- Prepare coordinates -----------
  # extract image coordinates as -(imagecols) | imagerow
  spatial <- lapply(seurat_lst, function(x) {
    df <- as.data.frame(x@images[[1]]@coordinates[, c(5, 4)])
    colnames(df) <- paste0("Spatial_", 1:2)
    df$Spatial_2 <- -df$Spatial_2
    df
  })
  
  # Make combinations of spatial coordinates, return as matrix; if n=1, simply
  # return df as matix.
  spatial_all <- lapply(seq_along(spatial), function(i) {
    tmp <- bind_rows(spatial[-i])
    tmp$Spatial_1 <- 0
    tmp$Spatial_2 <- 0
    as.matrix(rbind(spatial[[i]], tmp))
  })
  
  # ----------- Combine feature matrices -----------
  # Merge counts dfs, drop region; this could be where things break?
  combined_mat <- purrr::reduce(filtered_dfs, full_join, by = "region")
  rownames(combined_mat) <- combined_mat$region
  combined_mat$region <- NULL

  # remove cells not in intersection of combined metadata and spatial cords
  extra <- setdiff(colnames(combined_mat), rownames(spatial_all[[1]]))
  combined_mat <- combined_mat[, which(!colnames(combined_mat) %in% extra)]

  # ----------- Prepare Metadata -----------
  # remove assay-specific suffixes in obj meta.data (nCount_Spatial -> nCount)
  filtered <- lapply(filtered, function(x) {
    colnames(x@meta.data) <- gsub(
      paste0("_", Seurat::Assays(x)),
      "",
      colnames(x@meta.data)
    )
    x
  })

  # get the list of metadata from seurat objects
  list_of_metadata <- lapply(filtered, function(x) {
    x@meta.data
  })

  # combine meta data
  meta.data <- do.call("rbind", list_of_metadata)
  write.csv(meta.data, "req_meta_data.csv", row.names = TRUE)

  # ----------- Create Combined SeuratObj -----------
  combined <- Seurat::CreateSeuratObject(
    counts = combined_mat,
    assay = "scATAC",
    meta.data = meta.data
  )

  # make clusters factors
  combined@meta.data$Clusters <- factor(
    combined@meta.data$Clusters,
    levels = c(paste0("C", seq_along(unique(combined@meta.data$Clusters))))
  )

  # Add spatial coordinates as reductions
  for (i in seq_along(samples)) {
    spatial_all[[i]] <- spatial_all[[i]][colnames(combined), ]
    combined[[samples[i]]] <- Seurat::CreateDimReducObject(
      embeddings = spatial_all[[i]],
      key = samples[i],
      assay = Seurat::DefaultAssay(combined)
    )
  }

  # Get Variable Features
  combined <- Seurat::NormalizeData(
    combined, normalization.method = "LogNormalize", scale.factor = 10000
  )

  combined <- Seurat::FindVariableFeatures(
    combined, selection.method = "vst", nfeatures = 2000
  )

  # Add UMAP calculated previously to as a reduction
  combined[["UMAP"]] <- Seurat::CreateDimReducObject(
    embeddings = as.matrix(umap_embedding),
    key = "UMAP",
    assay = Seurat::DefaultAssay(combined)
  )

  # remove nFeature and nCounts
  combined@meta.data$nCount_scATAC <- NULL
  combined@meta.data$nCount <- NULL
  combined@meta.data$nFeature_scATAC <- NULL
  combined@meta.data$nFeature <- NULL

  return(combined)
}

main_func(seurat_objs, UMAPHarmony)

h5ad_file <- "output_data.h5ad"

# Convert the .h5ad file to a .h5seurat file
Convert(h5ad_file, dest = "h5seurat", assay = "scRNA", overwrite = TRUE)

# Load the .h5seurat file as a Seurat object
combined_h5 <- LoadH5Seurat("output_data.h5seurat")

# Print the Seurat object to verify
print(seurat_object)
