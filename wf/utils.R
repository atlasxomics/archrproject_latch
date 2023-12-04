library(ArchR)
library(ggplot2)
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

scvolcano <- function(inpMarkers){
  
  # Prepare ggData
  ggData = inpMarkers[cluster == input$sc1de1inp]
  minfdr = 0.09
  minfdr1 = 10^-(1/6 *(-log10(min(inpMarkers[cluster == input$sc1de1inp]$p_val_adj))))
  
  minfdr2 = 10^-(2/3 *(-log10(min(inpMarkers[cluster == input$sc1de1inp]$p_val_adj))))
  
  ggData$Significance = ifelse(ggData$p_val_adj < minfdr,
                               # & abs(ggData$avg_log2FC) >= 0.58,
                               ifelse(ggData$avg_log2FC > 0.0
                                      ,sc1def$Condition1,sc1def$Condition2)
                               ,'Not siginficant'
  )
  ggData$Significance <- factor(ggData$Significance
                                , levels = c(sc1def$Condition1,sc1def$Condition2,'Not siginficant'))
  
  ggData[p_val_adj < 1e-300]$p_val_adj = 1e-300
  ggData$log10fdr = -log10(ggData$p_val_adj)
  
  # Actual ggplot
  ggOut =
    ggplot(ggData, aes(avg_log2FC, log10fdr)) +
    geom_point() + sctheme() + ylab("-log10(FDR)")+  geom_point(aes(color = Significance)) +
    scale_color_manual(values = c("#F8766D","#619CFF","gray")) +
    geom_text_repel(data = subset(ggData, p_val_adj < minfdr1 ),aes(label = gene))
  return(ggOut)
}


