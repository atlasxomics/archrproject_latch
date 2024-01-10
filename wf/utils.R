library(ArchR)
library(ggplot2)
library(ggrepel)
library(harmony)
library(patchwork)
library(Seurat)

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
  matrix <- CreateAssayObject(matrix)

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

plot_feature <- function(seurat_obj, feature, name) {
  # Wrapper of Seurat's SpatialFeaturePlot with specific aesthetics

  SpatialFeaturePlot(
    object = seurat_obj,
    features = feature,
    alpha = c(0.2, 1),
    pt.size.factor = 1) +
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
  colors <- ArchRPalettes$stallion2[seq_len(length(clusters))]
  names(colors) <- clusters
  SpatialDimPlot(
    seurat_object,
    group.by = "Clusters",
    pt.size.factor = 1,
    image.alpha = 0,
    crop = FALSE,
    cols = colors,
    stroke = 0) +
    ggtitle(name) +
    theme(
      plot.title = element_text(size= 15),
      text = element_text(size = 10),
      legend.position = "bottom"
    )
}

plot_geneset <- function(seurat_obj, marker_genes, name, title) {
  # Return a Seurat SpatialFeaturePlot of the average expression score for a
  # set of genes.

  seurat_obj <- AddModuleScore(
    object = seurat_obj,
    features = marker_genes,
    name = name
  )
  plot <- SpatialFeaturePlot(
    object = seurat_obj,
    pt = 1,
    features = paste0(name, 1)) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5)
    )
}

plot_umap <- function(archrproj, name) {
  p <- plotEmbedding(
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
  if(!XYval) {
    oupTheme <- oupTheme + theme(
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }
  return(oupTheme)
}

scvolcano <- function(inpMarkers, feature = "All") {

  # Prepare ggData
  ggData <- inpMarkers[which(inpMarkers$cluster == feature), ]
  minfdr <- 0.09
  minfdr1 <- 10^-(1 / 6 * (-log10(min(ggData$p_val_adj))))

  minfdr2 <- 10^-(2 / 3 * (-log10(min(ggData$p_val_adj))))

  ggData$Significance <- ifelse(
    ggData$p_val_adj < minfdr,
    ifelse(
      ggData$avg_log2FC > 0.0,
      markerList@colData@rownames[[1]],
      markerList@colData@rownames[[2]]),
    "Not siginficant"
  )

  ggData$Significance <- factor(
    ggData$Significance,
    levels = c(
      markerList@colData@rownames[[1]],
      markerList@colData@rownames[[2]],
      "Not siginficant")
  )

  ggData[ggData$p_val_adj < 1e-300, "p_val_adj"] <- 1e-300
  ggData$log10fdr <- -log10(ggData$p_val_adj)

  # Actual ggplot
  ggOut <-
    ggplot(ggData, aes(avg_log2FC, log10fdr)) +
    geom_point() +
    sctheme() +
    ylab("-log10(FDR)") +
    geom_point(aes(color = Significance)) +
    scale_color_manual(values = c("#F8766D", "#619CFF", "gray")) +
    geom_text_repel(
      data = subset(ggData, p_val_adj < minfdr1),
      aes(label = gene)) +
    ggtitle(paste("Marker genes:", feature)) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20),
      legend.text = element_text(size = 15),
      legend.title = element_text(size = 18)
    )
  return(ggOut)
}