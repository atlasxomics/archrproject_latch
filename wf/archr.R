#' Rountines for manipulation of and plotting with ArchR

library("ArchR")
library("ggplot2")
library("ggplot2")
library("gridExtra")
library("harmony")
library("patchwork")
library("SummarizedExperiment")
library("S4Vectors")

source("/root/getDeviation_ArchR.R")
source("/root/wf/utils.R")


add_motif_annotations <- function(proj, genome) {
  #' Wrapper for ArchR::addMotifAnnotations()

  motif_set <- list("mm10" = "cisbp", "hg38" = "cisbp", "rnor6" = "encode")
  species <- list(
    "mm10" = NULL, "hg38" = NULL, "rnor6" = ArchR::getGenome(ArchRProj = proj)
  )

  proj <- ArchR::addMotifAnnotations(
    ArchRProj = proj,
    motifSet = motif_set[[genome]],
    name = "Motif",
    force = TRUE,
    species = species[[genome]]
  )
  return(proj)
}

create_archrproject <- function(
  inputs, genome, min_tss, min_frags, tile_size, out_dir
) {
  #' Create ArrowFiles and ArchRProject; handles mm10, hg38, rnor6.  Inputs are
  #' a named vector mapping run_id to local fragment files path.

  if (genome %in% c("mm10", "hg38")) {

    ArchR::addArchRGenome(genome)
    geneAnnotation <- ArchR::getGeneAnnotation()
    genomeAnnotation <- ArchR::getGenomeAnnotation()

  } else if (genome == "rnor6") {

    # load in geneAnnotation and genomeAnotation from custom_ArchR repo
    load(
    "/root/custom_ArchR_genomes_and_annotations/rn6/rn6_liftoff_mm10NcbiRefSeq_ArchR_annotations.rda"
    )

  } else {

    stop(
      "Genome not one for 'mm10', 'hg38', 'rnor6'; please supply correct
      genome."
    )
  }

  arrow_files <- ArchR::createArrowFiles(
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

  proj <- ArchR::ArchRProject(
    ArrowFiles = arrow_files,
    outputDirectory = out_dir,
    geneAnnotation = geneAnnotation,
    genomeAnnotation = genomeAnnotation
  )

  return(proj)
}

get_annotated_peaks <- function(proj, group_by, genome_size, genome) {

  proj <- ArchR::addGroupCoverages(
    ArchRProj = proj,
    groupBy = group_by,
    maxCells = 1500,
    force = TRUE
  )

  proj <- ArchR::addReproduciblePeakSet(
    ArchRProj = proj,
    groupBy = group_by,
    pathToMacs2 = ArchR::findMacs2(),
    genomeSize = genome_size,
    maxPeaks = 300000,
    force = TRUE
  )
  proj <- ArchR::addPeakMatrix(proj, force = TRUE)

  # Add motif annotations -----
  proj <- add_motif_annotations(proj, genome)

  return(proj)
}

get_enriched_motifs <- function(proj, marker_peaks, cutoff) {

  enrich_motifs <- ArchR::peakAnnoEnrichment(
    seMarker = marker_peaks,
    ArchRProj = proj,
    peakAnnotation = "Motif",
    cutOff = cutoff
  )
  enrich_df <- data.frame(enrich_motifs@assays@data)

  motif_lst <- unique(rownames(enrich_motifs))
  split_string <- strsplit(motif_lst, split = "\\(")

  req_motifs <- gsub("_", "-", extract_nth_ele(split_string)) # from utils.R
  req_motifs <- gsub(" ", "", req_motifs)

  rownames(enrich_motifs) <- req_motifs

  heatmap_em <- ArchR::plotEnrichHeatmap(
    enrich_motifs, n = 50, transpose = FALSE, returnMatrix = TRUE, cutOff = 2
  )

  motif_lst <- unique(rownames(heatmap_em))
  split_string <- strsplit(motif_lst, split = "\\(")

  req_motifs <- gsub("_", "-", extract_nth_ele(split_string))  # from utils.R
  req_motifs <- gsub(" ", "", req_motifs)

  rownames(heatmap_em) <- req_motifs

  return(list(
    enrich_df = enrich_df,
    enrich_motifs = enrich_motifs,
    heatmap_em = heatmap_em
  ))
}

get_marker_df <- function(
  proj, group_by, matrix, seq_names, max_cells, test_method
) {
  #' Return data frame of ArhcR marker features with:
  #' "avg_log2FC", "p_val", "p_val_adj", "gene", "cluster"

  # Get markers -----
  markers <- ArchR::getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = matrix,
    groupBy = group_by,
    bias = c("TSSEnrichment", "log10(nFrags)"),
    useSeqnames = seq_names,
    maxCells = max_cells,
    normBy = "none",
    testMethod = test_method
  )

  # Create data frame with for MeanDiff, Pval, FDR -----
  markers_df1 <- SummarizedExperiment::assay(markers, "MeanDiff")
  markers_df2 <- SummarizedExperiment::assay(markers, "Pval")
  markers_df3 <- SummarizedExperiment::assay(markers, "FDR")

  conditions <- sort(unique(proj@cellColData[group_by][, 1]))

  markers_df <- list()
  for (conds in conditions) {

    markers_df[[conds]] <- S4Vectors::DataFrame(
      markers_df1[[conds]],
      markers_df2[[conds]],
      markers_df3[[conds]]
    )
    markers_df[[conds]] <- as.data.frame(markers_df[[conds]])
    markers_df[[conds]]$genes <- SummarizedExperiment::rowData(markers)$name
    markers_df[[conds]]$cluster <- rep(
      "All", length(rownames(markers_df[[conds]]))
    )
    colnames(markers_df[[conds]]) <- c(
      "avg_log2FC", "p_val", "p_val_adj", "gene", "cluster"
    )
  }
  return(markers_df)
}

get_marker_df_clusters <- function(proj, clusters, group_by) {

  markers_by_cluster <- list()
  for (i in seq_along(clusters)) {

    proj_filter <- BiocGenerics::which(proj$Clusters == clusters[i])
    cells_subset <- proj$cellNames[proj_filter]
    proj_subset <- proj[cells_subset, ]
    n_cells <- length(proj_subset$cellNames)

    # Get gene score matrix each cluster separately -----
    markers_by_cluster[[i]] <- ArchR::getMarkerFeatures(
      ArchRProj = proj_subset,
      useMatrix = "GeneScoreMatrix",
      groupBy = group_by,
      bias = c("TSSEnrichment", "log10(nFrags)"),
      maxCells = n_cells,
      normBy = "none",
      testMethod = "ttest"
    )
  }
  names(markers_by_cluster) <- clusters

  # create markList_df for each cluster -----
  markerlist_df1 <- list()
  markerlist_df2 <- list()
  markerlist_df3 <- list()
  markerlist_df <- list()

  for (i in seq_along(clusters)) {

    cluster <- clusters[i]
    markerlist_df1[[i]] <- SummarizedExperiment::assay(
      markers_by_cluster[[i]], "Log2FC"
    )
    markerlist_df2[[i]] <- SummarizedExperiment::assay(
      markers_by_cluster[[i]], "Pval"
    )
    markerlist_df3[[i]] <- SummarizedExperiment::assay(
      markers_by_cluster[[i]], "FDR"
    )

    conditions <- sort(unique(proj@cellColData[group_by][, 1]))

    markerlist_df[[i]] <- list()
    for (conds in conditions) {

      markerlist_df[[i]][[conds]] <- S4Vectors::DataFrame(
        markerlist_df1[[i]][, conds],
        markerlist_df2[[i]][, conds],
        markerlist_df3[[i]][, conds]
      )

      markerlist_df[[i]][[conds]] <- as.data.frame(
        markerlist_df[[i]][[conds]]
      )
      markerlist_df[[i]][[conds]]$genes <- SummarizedExperiment::rowData(
        markers_by_cluster[[i]]
      )$name
      markerlist_df[[i]][[conds]]$cluster <- rep(
        cluster, dim(markerlist_df[[i]][[conds]])[1]
      )
      colnames(markerlist_df[[i]][[conds]]) <- c(
        "avg_log2FC", "p_val", "p_val_adj", "gene", "cluster"
      )
    }
  }

  names(markerlist_df) <- clusters
  markers_df_by_cluster <- do.call(Map, c(f = rbind, markerlist_df))

  return(markers_df_by_cluster)
}

get_marker_genes <- function(
  proj, group_by, markers_cutoff, heatmap_cutoff, rd_name
) {

  markers_gs <- ArchR::getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "GeneScoreMatrix",
    groupBy = group_by,
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "ttest"
  )
  marker_list <- ArchR::getMarkers(markers_gs, cutOff = markers_cutoff)

  proj <- ArchR::addImputeWeights(proj, reducedDims = rd_name)

  heatmap_gs <- ArchR::plotMarkerHeatmap(
    seMarker = markers_gs,
    cutOff = heatmap_cutoff,
    plotLog2FC = TRUE,
    returnMatrix = TRUE
  )
  return(
    list(
      markers_gs = markers_gs,
      marker_list = marker_list,
      heatmap_gs = heatmap_gs
    )
  )
}

get_marker_peaks <- function(proj, group_by, peak_data, cut_off) {

  marker_peaks <- ArchR::getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix",
    groupBy = group_by,
    bias = c("TSSEnrichment", "log10(nFrags)"),
    k = 100,
    testMethod = "wilcoxon"
  )

  marker_peak_list <- ArchR::getMarkers(marker_peaks, cutOff = cut_off)

  # Merge all peaks with significant marker peaks, write to csv -----
  total_peaks <- merge(peak_data, marker_peak_list, by = c("start", "end"))

  return(list(
    marker_peaks = marker_peaks,
    marker_peak_list = marker_peak_list,
    total_peaks = total_peaks
  ))
}

get_proj_medians <- function(proj) {
  #' Get dataframe with median TSS, nFrags, FRIP value per run.

  medians <- ArchR::getCellColData(ArchRProj = proj)
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

  return(merged_df)
}

get_required_clusters <- function(proj, group_by) {

  req_df <- as.data.frame(ArchR::getCellColData(proj))
  df1 <- table(req_df$Clusters, req_df[, group_by])
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
  req_clusters <- req_clusters[order(as.numeric(gsub("C", "", req_clusters)))]
  req_clusters <- req_clusters[which(!req_clusters %in% not_req_list)]

  return(req_clusters)
}

get_topn_hm_feats <- function(heatmap, n_groups, n_feats) {
  #' Return a vector with top n features from heatmap for each grouping
  #' (cluster, samples, conditions); vector will have length n_groups * nfeats.

  lst <- list()
  for (i in seq_along(1: n_groups)) {
    lst[[i]] <- heatmap[, c(1, i + 1)]
    lst[[i]] <- lst[[i]][
      order(lst[[i]][, 2], decreasing = TRUE),
    ][1:n_feats, 1]
  }

  # Flatten to vector, remove duplicated features, NAs -----
  req_feats <- unlist(lst)
  req_feats <- req_feats[!duplicated(req_feats)]
  req_feats <- na.omit(req_feats)

  return(req_feats)
}

get_volcano_table <- function(
  markers_df, markers_by_cluster_df, condition, feature, matrix
) {

  # Merge df with all clusters with df for each cluster -----
  merged_df <- rbind(
    markers_df[[condition]], markers_by_cluster_df[[condition]]
  )

  # Remove empty features -----
  if (feature == "gene") {
    gsm_mat <- SummarizedExperiment::assay(matrix, "GeneScoreMatrix")
    empty_feat_idx <- which(rowSums((gsm_mat)) == 0)
    empty_feat <- SummarizedExperiment::rowData(matrix)$name[empty_feat_idx]

  } else if (feature == "motif") {
    empty_feat_idx <- which(rowSums(matrix) == 0)
    empty_feat <- rownames(matrix)[empty_feat_idx]
  }
  merged_df <- merged_df[which(!merged_df$gene %in% empty_feat), ]

  # Remove na values -----
  merged_df <- na.omit(merged_df)

  # Remove FDR equal to 0 -----
  merged_df <- merged_df[which(!merged_df$p_val_adj == 0), ]

  # Make logfc limiation between 1 and -1 -----
  merged_df <- merged_df[which(abs(merged_df$avg_log2FC) < 1.2), ]

  # Get string of conditions != cond
  others <- paste(
    names(markers_df)[condition != names(markers_df)], collapse = "|"
  )

  # Significant if p-val < 0.01; if avg log(fold change) positive then
  # condition, else other codition
  merged_df$Significance <- ifelse(
    merged_df$p_val < 10^-2,
    ifelse(
      merged_df$avg_log2FC > 0.0,
      condition,
      others
    ),
    "Not siginficant"
  )

  return(merged_df)
}

plot_umap <- function(proj, name) {
  plot <- ArchR::plotEmbedding(
    ArchRProj = proj,
    colorBy = "cellColData",
    name = name,
    embedding = "UMAP"
  )
  return(plot)
}

save_umap <- function(proj, color_by) {
  #' Save a plots of UMAP embeddings as a single .pdf, color by each feature in
  #' character vector 'color_by'.

  umap_plots <- c()
  for (feature in color_by) {
    umap <- plot_umap(proj, feature)
    umap_plots[[feature]] <- umap
  }

  pdf("umap_plots.pdf")
  gridExtra::grid.arrange(grobs = umap_plots, ncol = 2)
  dev.off()
}