library("ArchR")
library("BPCells")
library("dplyr")
library("ggplot2")
library("ggrepel")
library("harmony")
library("patchwork")
library("purrr")
library("Seurat")
library("SeuratDisk")

find_func <- function(tempdir, pattern) {
  list.files(
    path = tempdir,     # replace tempdir with the directory you want
    pattern = pattern,  # has "test", followed by 0 or more characters,
    full.names = TRUE,  # then ".csv", and then nothing else ($)
    recursive = TRUE    # include the directory in the result
  )
}

build_atlas_seurat_object <- function(
  run_id,
  matrix,
  metadata,
  spatial_path
) {
  #' Prepare and combine gene matrix, metadata, and image for SeuratObject
  #' for runs within a project.

  # First garbage collect to free up memory before starting
  gc(verbose = FALSE, full = TRUE)

  # Read the image first, as it's likely smaller than the matrix manipulation
  image <- Seurat::Read10X_Image(
    image.dir = spatial_path,
    filter.matrix = TRUE
  )

  # Filter metadata
  metadata <- metadata[metadata$Sample == run_id, ]

  # Find indices for subsetting instead of creating a temporary vector with grep
  col_indices <- grep(pattern = run_id, colnames(matrix))

  # Create sparse matrix directly from the subset to minimize memory usage
  message("Subsetting matrix and converting to sparse format...")

  # Option 1: If matrix is already in memory and very large
  matrix_sparse <- as(matrix[, col_indices], "dgCMatrix")

  # Immediately remove references to free memory
  rm(col_indices)
  # Don't remove matrix here as it's an input argument

  # Force garbage collection
  gc(verbose = FALSE, full = TRUE)

  # Set column names
  colnames(matrix_sparse) <- rownames(metadata)

  # Create assay object and immediately remove the source matrix
  message("Creating Seurat assay object...")
  matrix_assay <- Seurat::CreateAssayObject(counts = matrix_sparse)
  rm(matrix_sparse)
  gc(verbose = FALSE, full = TRUE)

  # Create Seurat object
  message("Creating Seurat object...")
  object <- Seurat::CreateSeuratObject(
    counts = matrix_assay,
    assay = "Spatial",
    meta.data = as.data.frame(metadata)
  )

  # Clean up
  rm(matrix_assay)
  gc(verbose = FALSE, full = TRUE)

  # Add image
  message("Adding spatial image...")
  image <- image[Seurat::Cells(x = object)]
  Seurat::DefaultAssay(object = image) <- "Spatial"
  object[["slice1"]] <- image

  # Final cleanup
  rm(image)
  gc(verbose = FALSE, full = TRUE)

  return(object)
}

plot_feature <- function(seurat_obj, feature, name) {
  # Wrapper of Seurat's SpatialFeaturePlot with specific aesthetics

  Seurat::SpatialFeaturePlot(
    object = seurat_obj,
    features = feature,
    alpha = c(0.2, 1),
    pt.size.factor = 1,
    crop = FALSE
  ) +
    ggtitle(name) +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 15, hjust = 0.5),
      text = element_text(size = 10)
    )
}

plot_spatial <- function(seurat_object, name) {
  # Wrapper of Seurat's SpatialDimPlot with specific aesthetics

  clusters <- sort(unique(seurat_object$Clusters))
  colors <- ArchR::ArchRPalettes$stallion2[seq_len(length(clusters))]
  names(colors) <- clusters
  Seurat::SpatialDimPlot(
    seurat_object,
    group.by = "Clusters",
    pt.size.factor = 1,
    image.alpha = 0,
    crop = FALSE,
    cols = colors,
    stroke = 0
    ) +
      ggtitle(name) +
      theme(
        plot.title = element_text(size = 15),
        text = element_text(size = 10),
        legend.position = "bottom"
      )
}

plot_geneset <- function(seurat_obj, marker_genes, name, title) {
  #' Return a Seurat SpatialFeaturePlot of the average expression score for a
  #' set of genes.

  seurat_obj <- Seurat::AddModuleScore(
    object = seurat_obj,
    features = marker_genes,
    name = name
  )
  plot <- Seurat::SpatialFeaturePlot(
    object = seurat_obj,
    pt = 1,
    features = paste0(name, 1),
    crop = FALSE
  ) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5))
}

plot_umap <- function(archrproj, name) {
  p <- ArchR::plotEmbedding(
    ArchRProj = archrproj,
    colorBy = "cellColData",
    name = name,
    embedding = "UMAP"
  )
  return(p)
}

sctheme <- function(base_size = 24, XYval = TRUE, Xang = 0, XjusH = 0.5) {
  oupTheme <- theme(
    text = element_text(size = base_size, family = "Helvetica"),
    panel.background = element_rect(fill = "white", colour = NA),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black", size = base_size / 20),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = base_size),
    axis.text.x = element_text(angle = Xang, hjust = XjusH),
    legend.position = "bottom",
    legend.key = element_rect(colour = NA, fill = NA),
  )
  if (!XYval) {
    oupTheme <- oupTheme + theme(
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }
  return(oupTheme)
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

add_motif_annotations <- function(proj, genome) {
  if (genome == "hg38" || genome == "mm10") {
    proj <- ArchR::addMotifAnnotations(
      ArchRProj = proj,
      motifSet = "cisbp",
      name = "Motif",
      force = TRUE
    )
  } else {
    proj <- ArchR::addMotifAnnotations(
      ArchRProj = proj,
      motifSet = "encode",
      name = "Motif",
      force = TRUE,
      species = ArchR::getGenome(ArchRProj = proj)
    )
  }
  return(proj)
}

addclust <- function(archrproj, resolution, reduceddims_name, min_cells) {
  #' Ensure all clusters >= min_cells; if less than min_cells, decrease
  #' clustering resolution by 0.1 and repeat clustering until resoultion = 0.
  #'  Return an ArhcRProject.

  cluster_df <- as.data.frame(table(archrproj$Clusters))

  while (min(cluster_df$Freq) <= min_cells && resolution > 0) {

    # Update the value in each step
    resolution <- resolution - 0.1
    print(paste("Changing the resolution to", resolution))

    archrproj <- ArchR::addClusters(
      input = archrproj,
      reducedDims = reduceddims_name,
      method = "Seurat",
      name = "Clusters",
      resolution = resolution,
      force = TRUE
    )

    cluster_df <- as.data.frame(table(archrproj$Clusters))

    print(paste("With resolution equal to", resolution))
    print(cluster_df)
  }
  return(archrproj)
}

find_samples_name <- function(seurat_lst) {
  # Extract list of sample names from list of SeuratObjs.
  sapply(
    seq_along(seurat_lst), function(i) {
      unique(seurat_lst[[i]]@meta.data$Sample)
    }
  )
}

combine_objs <- function(
  seurat_lst, umap_embedding, samples, spatial, project_name
) {

  # ----------- Filter SeuratObjs  -----------
  # Remove cell with NA counts from each SeuratObj
  to_remove <- lapply(seurat_lst, function(x) {
    names(which(colSums(is.na(x@assays[[1]]@counts)) > 0))
  })
  filtered <- mapply(
    function(x, y) x[, !colnames(x) %in% y], seurat_lst, to_remove
  )

  # ----------- Prepare coordinates -----------
  # Make combinations of spatial coordinates, return as matrix; if n=1, simply
  # return df as matix.
  spatial_all <- lapply(seq_along(spatial), function(i) {
    tmp <- dplyr::bind_rows(spatial[-i])
    tmp$Spatial_1 <- 0
    tmp$Spatial_2 <- 0
    as.matrix(rbind(spatial[[i]], tmp))
  })

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

  ncells  <- as.double(sum(sapply(filtered, ncol)))
  nfeats <- as.double(nrow(seurat_lst[[1]]))
  matrix_size <- ncells * nfeats

  if (matrix_size < 2^31 - 1) {

    print("Feature matrix less than 2^31 -1...")

    # Convert list of SeuratObjs to list of feat counts as data frames.
    filtered_dfs <- lapply(filtered, function(x) {
      df <- as.data.frame(x@assays[[1]]@counts)
      colnames(df) <- Seurat::Cells(x)
      df$region <- rownames(df)
      df
    })

    # Merge counts dfs, drop region
    combined_mat <- purrr::reduce(filtered_dfs, full_join, by = "region")
    rownames(combined_mat) <- combined_mat$region
    combined_mat$region <- NULL

    combined <- Seurat::CreateSeuratObject(
      counts = combined_mat,
      assay = "scATAC",
      meta.data = meta.data
    )

  } else {

    print("Feature matrix greater than 2^31 -1; using BPCells...")

    # Use BPCells to handle large matrices
    # Convert SeuratObjs to BPCells files
    counts_mat <- c()
    for (i in seq_along(filtered)) {

      path <- paste0(project_name, "/", samples[[i]], "_BP")

      BPCells::write_matrix_dir(
        mat = filtered[[i]]@assays[["Spatial"]]@counts,
        dir = path,
        overwrite = TRUE
      )
      counts_mat[[i]] <- BPCells::open_matrix_dir(dir = path)
    }

    # Create combiend SeuratObject
    combined <- CreateSeuratObject(
      counts = counts_mat,
      assay = "scATAC",
      meta.data = meta.data
    )

    combined <- SeuratObject::JoinLayers(combined)
  }

  # Normalize and calculate variable features
  combined <- Seurat::NormalizeData(
    combined, normalization.method = "LogNormalize", scale.factor = 10000
  )

  combined <- Seurat::FindVariableFeatures(
    combined, selection.method = "vst", nfeatures = 2000, slot = "counts"
  )

  # Make clusters factors
  combined@meta.data$Clusters <- factor(
    combined@meta.data$Clusters,
    levels = c(paste0("C", seq_along(unique(combined@meta.data$Clusters))))
  )

  # Add spatial coordinates as embeddings
  # Remove cells that remain in coordinates but not in combined (had NA counts)
  include <- intersect(Seurat::Cells(combined), rownames(spatial_all[[1]]))
  combined <- subset(combined, cells = include)

  # Add embeddings to combined
  for (i in seq_along(samples)) {
    spatial_all[[i]] <- spatial_all[[i]][colnames(combined), ]
    combined[[samples[i]]] <- Seurat::CreateDimReducObject(
      embeddings = spatial_all[[i]],
      key = samples[i],
      assay = Seurat::DefaultAssay(combined)
    )
  }

  # Add UMAP calculated previously to as a reduction
  combined[["UMAP"]] <- Seurat::CreateDimReducObject(
    embeddings = as.matrix(umap_embedding),
    key = "UMAP",
    assay = Seurat::DefaultAssay(combined)
  )

  # Remove unused metadata columns
  combined@meta.data$nCount_scATAC <- NULL
  combined@meta.data$nCount <- NULL
  combined@meta.data$nFeature_scATAC <- NULL
  combined@meta.data$nFeature <- NULL

  return(combined)
}
