library("ArchR")
library("BSgenome")
library("BSgenome.Hsapiens.UCSC.hg38")
library("BSgenome.Mmusculus.UCSC.mm10")
library("gridExtra")
library("harmony")
library("purrr")
library("Seurat")
library("GenomicRanges")
library("dplyr")                                    # Load dplyr package
library("plyr")                                     # Load plyr package
library("readr")
library("qdap")
library("ShinyCell")
library("seqLogo")
require("ggseqlogo")
library("chromVARmotifs")
suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(library("ComplexHeatmap"))
suppressPackageStartupMessages(library("circlize"))
suppressPackageStartupMessages(library(data.table))

source("/root/makeShinyFiles.R")
source("/root/getDeviation_ArchR.R")
source("/root/wf/utils.R")

find_func <- function(tempdir,pattern){
  
  list.files(
    path = tempdir, # replace with the directory you want
    pattern = pattern, # has "test", followed by 0 or more characters,
    # then ".csv", and then nothing else ($)
    full.names = TRUE # include the directory in the result
    , recursive = TRUE
  )
}

# globals ---------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

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
print(num_threads)

runs <- strsplit(args[12:length(args)], ",")
runs

inputs <- c()
for (run in runs) {
  inputs[run[1]] <- run[2]
}
inputs

out_dir <- paste0(project_name, "_ArchRProject")

# create archr project --------------------------------------------------------

addArchRGenome(genome)
addArchRThreads(threads = 50)

arrow_files <- createArrowFiles(
  inputFiles = inputs,
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

proj <- ArchRProject(
  ArrowFiles = arrow_files,
  outputDirectory = out_dir
)

# Add an additional Conditions column
for (run in runs) {
  proj$Condition[proj$Sample == run[1]] <- run[3]
}

# Filter on-tissue
all_ontissue <- c()
for (run in runs) {
  positions <- read.csv(run[4], header = FALSE)
  positions$V1 <- paste(run[1], "#", positions$V1, "-1", sep = "")
  on_tissue <- positions$V1 [which(positions$V2 == 1)]
  all_ontissue <- c(all_ontissue, on_tissue)
}

proj <- proj[proj$cellNames %in% all_ontissue]

# dimension reduction and clustering ------------------------------------------

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

if (length(runs) > 1) {
  proj <- addHarmony(
    ArchRProj = proj,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample",
    force = TRUE
  )
  name <- "Harmony"
} else {
  name <- "IterativeLSI"
}

proj <- addClusters(
  input = proj,
  reducedDims = name,
  method = "Seurat",
  name = "Clusters",
  resolution = c(clustering_resolution),
  force = TRUE
)

######check how many cells per clusters are
cluster_df <- as.data.frame(table(proj$Clusters))
resolution = c(clustering_resolution)

addclust <- function(x) {
  while (min(cluster_df$Freq) <= 20) {
    print(paste0('with resolution equal to ',resolution ))
    cluster_df <- as.data.frame(table(x$Clusters))
    print(cluster_df)
    # Update the value in each step
    
    resolution <- resolution - 0.1
    print(paste0('changing the resolution to ',resolution )) 
    
    x <- addClusters(
      input = x,
      reducedDims = name,
      method = "Seurat",
      name = "Clusters",
      resolution = resolution ,
      force = TRUE
    )   
    cluster_df <- as.data.frame(table(x$Clusters))
    print(cluster_df)
  }
  return(x)
}

proj_2 <- addclust(proj)
table(proj_2$Clusters)

proj <- proj_2 
##################

proj <- addUMAP(
  ArchRProj = proj,
  reducedDims = name,
  name = "UMAP",
  nNeighbors = 30,
  minDist = umap_mindist,
  metric = "cosine",
  force = TRUE
)

umap_plots <- c()
for (name in c("Clusters", "Condition", "Sample")) {
  umap <- plot_umap(proj, name) # from utils
  umap_plots[[name]] <- umap
}

pdf("umap_plots.pdf")
grid.arrange(
  umap_plots[["Clusters"]],
  umap_plots[["Condition"]],
  umap_plots[["Sample"]],
  ncol = 2
)
dev.off()

proj <- addImputeWeights(proj)

conds <- strsplit(proj$Condition,split = "\\s|-")

fun1 <- function(lst,n){
  sapply(lst, `[`, n)}

for (i in seq_along(conds[[1]])){
  proj <- addCellColData(proj
                         ,data = fun1(conds,i)
                         ,name = paste0('condition_',i)
                         , cells = proj$cellNames
                         ,force = TRUE)
}
treatment <- names(getCellColData(proj))[grep('condition_',names(getCellColData(proj)))]


print('++++ what is the treatments in this project +++++++')
treatment

# create seurat objects -------------------------------------------------------

# create metadata object for Seurat object
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

# create gene matrix for Seurat object
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

print("+++++++++++creating seurat objs++++++++++++++")

seurat_objs <- c()
for (run in runs) {

  obj <- build_atlas_seurat_object(
    run_id = run[1],
    matrix = matrix,
    metadata = metadata,
    spatial_path = run[5]
  )

  saveRDS(obj, file = paste0(run[1], "_SeuratObj.rds"))
  seurat_objs <- c(seurat_objs, obj)
}

print("+++++++++++creating spatial plots++++++++++++++")

run_ids <- unique(proj$Sample)
spatial_cluster_plots <- list()
for (i in seq_along(run_ids)){
  plot <- plot_spatial(seurat_objs[[i]], run_ids[i])
  spatial_cluster_plots[[i]] <- plot
}

pdf("spatial_plots.pdf")
CombinePlots(spatial_cluster_plots, legend = "bottom")
dev.off()

print("+++++++++++creating qc plots++++++++++++++")

# feature plot from utils
metrics <- c("TSSEnrichment", "nFrags", "log10_nFrags")
all_qc_plots <- list()
for (i in seq_along(metrics)) {
  spatial_qc_plots <- list()
  for (j in seq_along(run_ids)) {
    plot <- plot_feature(seurat_objs[[j]], metrics[i], run_ids[j])
    spatial_qc_plots[[j]] <- plot
  }
  all_qc_plots[[i]] <- spatial_qc_plots
}

print('this is available seurat_objs')
seurat_objs

all <-  list()
for (i in seq_along(seurat_objs)){
  all[[i]] <- seurat_objs[[i]]
  
  all[[i]] <- RenameCells(all[[i]], new.names = paste0(unique(all[[i]]@meta.data$Sample),"#",colnames(all[[i]]),"-1"))
}

pdf("qc_plots.pdf")
for (i in seq_along(metrics)) {
  grid.arrange(
    grobs = all_qc_plots[[i]],
    ncol = 2
  )
}
dev.off()


# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# save ArchR object
saveArchRProject(
  ArchRProj = proj,
  outputDirectory = paste0(project_name, "_ArchRProject")
)
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@



############-------------Identifying Marker Genes------------###################
# per cluster

markersGS <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "GeneScoreMatrix",
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "ttest",#"wilcoxon"
  
)

saveRDS(markersGS,"markersGS_clusters.rds")

marker_list <- getMarkers(markersGS, cutOff = "FDR <= 0.02")
write.csv(
  marker_list,
  file = "marker_genes_per_cluster.csv",
  row.names = FALSE
)

for (rd in names(proj@reducedDims)){
  if (rd == "Harmony") {
    proj <- addImputeWeights(proj, reducedDims = "Harmony")
  } else {
    proj <- addImputeWeights(proj)
  }
}
# save for shiny

heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS,
  cutOff = "FDR <= 0.05 & Log2FC >= 0.20",
  plotLog2FC = TRUE,
  transpose = FALSE,
  returnMatrix = TRUE
)

write.csv(heatmapGS, "genes_per_cluster_hm.csv")

######################### Plot marker genes ###################################

heatmaps <- list()

# recompute for heatmap plot
gene_cutoff <- "FDR <= 0.05 & Log2FC >= 0.2"
heatmap_gs_plotting <- plotMarkerHeatmap(
  seMarker = markersGS,
  cutOff = gene_cutoff,
  transpose = TRUE
)
# save for plotting with peaks and motifs
gene_hm <- ComplexHeatmap::draw(
  heatmap_gs_plotting,
  heatmap_legend_side = "bot",
  annotation_legend_side = "bot",
  column_title = paste0("Marker genes (", gene_cutoff, ")"),
  column_title_gp = gpar(fontsize = 12)
)
heatmaps[[1]] <- gene_hm

###############################################################################

# per sample

######-----------Identifying Marker Genes-------------#######

if (length(unique(proj$Sample))>1){
  
markersGS <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Sample",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "ttest",#"wilcoxon"
)

# save for shiny app
saveRDS(markersGS,"markersGS_sample.rds")

marker_list <- getMarkers(markersGS, cutOff = "FDR <= 0.02")
write.csv(
  marker_list,
  file = "marker_genes_per_sample.csv",
  row.names = FALSE
)


for (rd in names(proj@reducedDims)){
  if (rd=="Harmony"){
    proj <- addImputeWeights(proj,reducedDims = "Harmony")
    
  } else {
    proj <- addImputeWeights(proj)
    
  }
}


  
  # save for shiny
  
  heatmapGS <- plotMarkerHeatmap(
    seMarker = markersGS, 
    cutOff = "FDR <= 0.05 & Log2FC >= 0.20", plotLog2FC = TRUE,
    #   labelMarkers = seleceted_markers,
    transpose = F,  returnMatrix = TRUE
  )
  write.csv(heatmapGS,"genes_per_sample_hm.csv")
  
} else {
  
  heatmapGS <- "there is not enough samples to be compared with!"  
}   



# per treatment

######---------------------Identifying Marker Genes----------------------#######
if (length(unique(proj$Condition))>1){
  
for (i in seq_along(treatment)){
  
  markersGS <- getMarkerFeatures(
    ArchRProj = proj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = treatment[i],
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "ttest",#"wilcoxon"
  )
  # save for shiny app
  saveRDS(markersGS,paste0("markersGS_treatment_",i,".rds"))

  marker_list <- getMarkers(markersGS, cutOff = "FDR <= 0.02")
  write.csv(
    marker_list,
    file = paste0("marker_genes_per_treatment_", i, ".csv"),
    row.names = FALSE
  )
  
  
  for (rd in names(proj@reducedDims)){
    if (rd=="Harmony"){
      proj <- addImputeWeights(proj,reducedDims = "Harmony")
      
    } else {
      proj <- addImputeWeights(proj)
      
    }
  }
  
  
  # save for shiny
  
 
    
    heatmapGS <- plotMarkerHeatmap(
      seMarker = markersGS, 
      cutOff = "FDR <= 0.05 & Log2FC >= 0.20", 
      plotLog2FC = TRUE,
      transpose = F,  
      returnMatrix = TRUE
    )
    write.csv(heatmapGS,paste0("genes_per_treatment_",i,"_hm.csv")  )
}
    
  } else {
    
    heatmapGS <- "there is not enough Conditions to be compared with!"  
  }   
  
  
                    

# Volcano plots for genes
# Volcano plots for genes
if (length(unique(proj$Condition))>1){
  for (j in seq_along(treatment)){
    
    ncells <- length(proj$cellNames)
    
    markerList <- getMarkerFeatures(
      ArchRProj = proj,
      useMatrix = "GeneScoreMatrix", 
      groupBy = treatment[j],
      bias = c("TSSEnrichment", "log10(nFrags)"),maxCells = ncells ,normBy = "none",
      testMethod = "ttest")
    
    req_DF <- as.data.frame(getCellColData(proj))
    df1 <- table(req_DF$Clusters,req_DF[,treatment[j]])
    distr <- as.data.frame.matrix(round(prop.table(as.matrix(df1),1),2))
    lst <- list()
    
    for(i in 1:nrow(distr)) {
      row <- distr[i,]
      if (
        sum(unname(unlist(row))>= 0.90) == 1) {
        rownames(row) -> lst[[i]]
      }
    }
    not_req_list <- unlist(lst)
    
    req_clusters <- unique(proj$Clusters)
    req_clusters <- req_clusters[order(as.numeric(gsub("C","",req_clusters)))]
    req_clusters <- req_clusters[which(!req_clusters%in%not_req_list)]
    
    markerList_C <- list()
    proj_C <- list()
    
    for (i in seq_along(req_clusters)) {
      
      idxSample <- BiocGenerics::which(proj$Clusters == req_clusters[i])
      
      cellsSample <- proj$cellNames[idxSample]
      proj_C[i] <- proj[cellsSample,]
      
      ncells[i] <- length(proj_C[[i]]$cellNames)
      
      # per each cluster separately
      markerList_C[[i]] <- getMarkerFeatures(
        ArchRProj = proj_C[[i]],
        useMatrix = "GeneScoreMatrix", 
        groupBy = treatment[j],
        bias = c("TSSEnrichment", "log10(nFrags)")
        ,maxCells = ncells[[i]] ,
        normBy = "none",
        testMethod = "ttest")
    }
    names(markerList_C) <- req_clusters
    
    
    gsm <- getMatrixFromProject(proj)
    gsm_mat <- assay(getMatrixFromProject(proj),"GeneScoreMatrix")
    
    which(rowSums(is.na(gsm_mat))>0)
    any(rowSums((gsm_mat))==0)
    
    empty_gene_idx <- which(rowSums((gsm_mat))==0)
    empty_gene <- rowData(gsm)$name[empty_gene_idx]
    
    
    
    markerList_df1 <- assay(markerList, "Log2FC")
    markerList_df2 <- assay(markerList, "Pval")
    markerList_df3 <- assay(markerList, "FDR")
    markerList_df <- cbind(markerList_df1,markerList_df2,markerList_df3)
    markerList_df$genes<- rowData(markerList)$name
    markerList_df$cluster <- rep("All",length(rownames(markerList_df)))
    
    # # we only want to see results of one set , say just sham
    markerList_df <- markerList_df[,c(1,3,5,7,8)]
    colnames(markerList_df) <- c("avg_log2FC","p_val","p_val_adj","gene","cluster")
    
    # percluster
    
    
    markerList_df1_C <- list()
    markerList_df2_C <- list()
    markerList_df3_C <- list()
    markerList_df_C <- list()
    
    # for (i in (1:nClust)){
    for (i in seq_along(req_clusters)){
      
      cluster <- req_clusters[i]
      
      markerList_df1_C[[i]] <- assay(markerList_C[[i]], "Log2FC")
      markerList_df2_C[[i]] <- assay(markerList_C[[i]], "Pval")
      markerList_df3_C[[i]] <- assay(markerList_C[[i]], "FDR")
      
      markerList_df_C[[i]] <- cbind(markerList_df1_C[[i]]
                                    ,markerList_df2_C[[i]]
                                    ,markerList_df3_C[[i]])
      
      markerList_df_C[[i]]$genes<- rowData(markerList_C[[i]])$name
      markerList_df_C[[i]]$cluster <- rep(cluster,length(rownames(markerList_df_C[[i]])))
      
      # # we only want to see results of one set , say just sham
      markerList_df_C[[i]] <- markerList_df_C[[i]][,c(1,3,5,7,8)]
      colnames(markerList_df_C[[i]]) <- c("avg_log2FC","p_val","p_val_adj","gene","cluster")
    }
    
    
    names(markerList_df_C) <- req_clusters
    
    markersGS_merged_df <- do.call("rbind", markerList_df_C)
    
    # also data frame for all clusters together needs to be added
    
    markersGS_merged_df <- rbind(markerList_df,markersGS_merged_df)
    
    # remove empty genes
    markersGS_merged_df <- markersGS_merged_df[which(!markersGS_merged_df$gene%in%empty_gene),]
    
    # remove na values
    markersGS_merged_df <- na.omit(markersGS_merged_df)
    
    # remove FDR equal to 0
    markersGS_merged_df <- markersGS_merged_df[which(!markersGS_merged_df$p_val_adj== 0),]
    
    
    
    # make logfc limiation between 1 and -1
    
    markersGS_merged_df <- markersGS_merged_df[which(abs(markersGS_merged_df$avg_log2FC)< 1.2),]
    
    markersGS_merged_df$Significance = ifelse(markersGS_merged_df$p_val < 10^-2 , 
                                              #                                           abs(markersGS_merged_df$avg_log2FC) >= 0.58, 
                                              ifelse(markersGS_merged_df$avg_log2FC> 0.0 
                                                     ,colnames(markerList)[1],colnames(markerList)[2]),
                                              'Not siginficant')
    
    de <- markersGS_merged_df

    print(paste0("volcanoMarkers_genes_", j, ".txt is writing!"))
    write.table(
      de,
      paste0("volcanoMarkers_genes_", j, ".txt"),
      sep = '\t',
      quote = FALSE,
      row.names = FALSE
    )
    print(paste0("writing volcanoMarkers_genes_", j, ".txt is done!"))

    features <- unique(de$cluster)
    volcano_plots <- list()
    for (i in seq_along(features)) {
      volcano_plots[[i]] <- scvolcano(de, markerList, features[[i]])
    }

    pdf(paste0("volcano_plots_", j, ".pdf"))
    for (plot in volcano_plots) {
      print(plot)
    }
    dev.off()
  }

} else {
  
  de <- "there is not enough conditions to be compared with!"  
}   


tempdir <- "/root"
genes_per_cluster_hm <- find_func(tempdir,"genes_per_cluster_hm.csv")
hm_per_clust <- read.csv(genes_per_cluster_hm)




nClust = length(unique(proj$Clusters))


df = list()


for (i in seq_along(1:nClust)){
  df[[i]] <- hm_per_clust[,c(1,i+1)]
  
  #select top 5 values by group
  
  df[[i]] <- df[[i]][order(df[[i]][,2], decreasing = T),][1:5,1]
  
}
final <- do.call(rbind, df)

req_genes1 <- unlist(df)

req_genes1<- req_genes1[!duplicated(req_genes1)]

req_genes1<- na.omit(req_genes1)


# save the genes for default values in shiny app

write.csv(req_genes1,"req_genes1.csv")




# per sample

if (length(unique(proj$Sample))>1){
  
  genes_per_sample_hm <- find_func(tempdir,"genes_per_sample_hm.csv")
  hm_per_sample <- read.csv(genes_per_sample_hm)
  
  
  nSamples = length(unique(proj$Sample))
  
  
  df = list()
  
  for (i in seq_along(1:nSamples)){
    df[[i]] <- hm_per_sample[,c(1,i+1)]
    
    #select top 10 values by group
    
    df[[i]] <- df[[i]][order(df[[i]][,2], decreasing = T),][1:10,1]
    
  }
  final <- do.call(rbind, df)
  # final
  req_genes3 <- unlist(df)
  
  req_genes3<- req_genes3[!duplicated(req_genes3)]
  
  req_genes3<- na.omit(req_genes3)
  
  write.csv(req_genes3,"req_genes3.csv")
} else {
  
  req_genes3 <- "there is not enough samples to be compared with!"  
}





# per treatment

if (length(unique(proj$Condition))>1){
  genes_per_cond_hm <- find_func(tempdir,"genes_per_treatment_*")
  
  for (j in seq_along(genes_per_cond_hm)){
    hm_per_cond <- read.csv(genes_per_cond_hm[j])
    
    nConds = 2
    
    
    df = list()
    
    for (i in seq_along(1:nConds)){
      df[[i]] <- hm_per_cond[,c(1,i+1)]
      
      #select top 20 values by group
      
      df[[i]] <- df[[i]][order(df[[i]][,2], decreasing = T),][1:20,1]
      
    }
    final <- do.call(rbind, df)
    # final
    req_genes2 <- unlist(df)
    
    req_genes2<- req_genes2[!duplicated(req_genes2)]
    
    req_genes2<- na.omit(req_genes2)
    write.csv(req_genes2,paste0("req_genes2_",j,".csv"))
  }
  
} else {
  
  req_genes2 <- "there is not enough condition to be compared with!"  
}   



UMAPHarmony <- getEmbedding(ArchRProj = proj, embedding = "UMAP", returnDF = TRUE)
write.csv(UMAPHarmony,"UMAPHarmony.csv")


# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# save ArchR object
saveArchRProject(
  ArchRProj = proj,
  outputDirectory = paste0(project_name, "_ArchRProject")
)
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# since sometimes we get outofMemory error we decrease number of threads here:
addArchRThreads(threads = num_threads)

# peak calling with MACS2 for clusters------------------------------------------


proj <- addGroupCoverages(
  ArchRProj = proj,
  groupBy = "Clusters",
  maxCells = 1500,
  force = TRUE
)
# get genome size
species <- getGenome(ArchRProj = proj)

if (species == "BSgenome.Hsapiens.UCSC.hg38"){
  genome_size <- 3.3e+09
} else if (species == "BSgenome.Mmusculus.UCSC.mm10") {
  genome_size = 3.0e+09
}

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

# Motif enrichment (Deviation)


if("Motif" %ni% names(proj@peakAnnotation)){
  if (
    species == "BSgenome.Hsapiens.UCSC.hg38" || species == "BSgenome.Mmusculus.UCSC.mm10") {
    proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif", force = TRUE)
  } else {
    proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "encode", name = "Motif", force = TRUE
                                 , species = getGenome(ArchRProj = proj))
  }
}


# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# save ArchR object
saveArchRProject(
  ArchRProj = proj,
  outputDirectory = paste0(project_name, "_ArchRProject")
)
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

######################### get marker peaks, save ##############################

markersPeaks <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  k = 100,
  testMethod = "wilcoxon"
)

peak_marker_list <- getMarkers(markersPeaks, cutOff = "FDR <= 0.02")
write.csv(
  peak_marker_list,
  file = "marker_peaks_per_cluster.csv",
  row.names = FALSE
)

markersPeaks_sample <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  groupBy = "Sample",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  k = 100,
  testMethod = "wilcoxon"
)

peak_marker_list_s <- getMarkers(markersPeaks_sample, cutOff = "FDR <= 0.02")
write.csv(
  peak_marker_list_s,
  file = "marker_peaks_per_sample.csv",
  row.names = FALSE
)

######################### Plot marker peaks, motifs ###########################

peak_cutoff <- "Pval <= 0.05 & Log2FC >= 0.1"
heatmap_peaks <- plotMarkerHeatmap(
  seMarker = markersPeaks,
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

motifs_cutoff <- "Pval <= 0.05 & Log2FC >= 0.1"
enrichMotifs <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = proj,
  peakAnnotation = "Motif",
  cutOff = motifs_cutoff
)

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

print("+++++++++++creating heatmap plots++++++++++++++")

pdf("heatmaps_all.pdf")
for (i in seq_along(heatmaps)) {
  print(heatmaps[[i]])
}
dev.off()

###############################################################################

motif_lst <- unique(rownames(enrichMotifs))
split_string <- strsplit(motif_lst, split = "\\(")
fun1 <- function(list, nth){
  sapply(list, `[` , 1)
}
req_motifs1 <- gsub("_","-",fun1(split_string))
req_motifs1 <- gsub(" ","",req_motifs1)

rownames(enrichMotifs) <- req_motifs1

saveRDS(enrichMotifs,"enrichMotifs_clusters.rds")

# cutOff A numeric cutOff that indicates the minimum P-adj enrichment to be included in the heatmap. default is 20 but we decrease that!

heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 50, transpose = F,returnMatrix = TRUE, cutOff = 2)

motif_lst <- unique(rownames(heatmapEM))
split_string <- strsplit(motif_lst, split = "\\(")
fun1 <- function(list, nth){
  sapply(list, `[` , 1)
}
req_motifs1 <- gsub("_","-",fun1(split_string))
req_motifs1 <- gsub(" ","",req_motifs1)

rownames(heatmapEM) <- req_motifs1

write.csv(heatmapEM,"motif_per_cluster_hm.csv")

tempdir <- "/root"
motifs_per_cluster_hm <- find_func(tempdir,"motif_per_cluster_hm.csv")
hm_per_clust <- read.csv(motifs_per_cluster_hm)


nClust = length(unique(proj$Clusters))

df = list()

for (i in seq_along(1:nClust)){
  df[[i]] <- hm_per_clust[,c(1,i+1)]
  
  #select top 5 values by group
  
  df[[i]] <- df[[i]][order(df[[i]][,2], decreasing = T),][1:5,1]
  
}
final <- do.call(rbind, df)
req_motifs1 <- unlist(df)

req_motifs1<- req_motifs1[!duplicated(req_motifs1)]

req_motifs1<- na.omit(req_motifs1)


# save the motifs for default values in shiny app

write.csv(req_motifs1,"req_motifs1.csv")

######################### Add Motifs Matrix and Projections ####################

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
# create deviation score 

markerMotifsList <- getMarkers(markersMotifs, cutOff = "FDR < 0.9 & MeanDiff >= 0")

motifs <- list()
for (i in seq_len(length(markerMotifsList))) {
  if (length(markerMotifsList[[i]]$name)>1) {
    motifs <- c(motifs, markerMotifsList[[i]]$name)
    #     motifs <- c(motifs, c(cluster1112_mm$name,cluster10_mm$name))
  }
}

if (length(motifs)>1) {
  motifs <- unlist(motifs)
  motifs <- paste0('z:', motifs)
  motifs <- unique(motifs)
  
  
  proj <- addImputeWeights(proj)
  
  dev_score <- getDeviation_ArchR(ArchRProj = proj
                                  , name = motifs
                                  , imputeWeights = getImputeWeights(proj))
  

  dev_score[is.na(dev_score)] <- 0 #min(dev_score, na.rm = TRUE)
  
}

dev_score2 <- dev_score[!is.infinite(rowSums(dev_score)),]
colnames(dev_score2) <- rownames(getCellColData(proj))
# remove 0 deviations per All samples 
all_zero <- names(which(rowSums(dev_score2)==0))
dev_score2 <- dev_score2[which(!rownames(dev_score2)%in%c(all_zero)),]
# convert to dgCmatrix
dev_score3 <- Matrix(as.matrix(dev_score2), sparse = TRUE) 

# create metadata object for Seurat object
metadata <- getCellColData(ArchRProj = proj)
rownames(metadata) <- str_split_fixed(
  str_split_fixed(
    row.names(metadata),
    "#",
    2)[, 2],
  "-",
  2)[, 1]
metadata["log10_nFrags"] <- log(metadata$nFrags)

# create seurat objects -------------------------------------------------------

seurat_objs <- c()

for (run in runs) {
  
  obj <- build_atlas_seurat_object(
    run_id = run[1],
    matrix = dev_score3,
    metadata = metadata,
    spatial_path = run[5]
  )
  
  saveRDS(obj, file = paste0(run[1], "_SeuratObjMotif.rds"))
  seurat_objs <- c(seurat_objs, obj)
}

print('this is available seurat_objMotifs')
seurat_objs

all_m <-  list()
for (i in seq_along(seurat_objs)){
  all_m[[i]] <- seurat_objs[[i]]
  
  all_m[[i]] <- RenameCells(all_m[[i]], new.names = paste0(unique(all_m[[i]]@meta.data$Sample),"#",colnames(all_m[[i]]),"-1"))
}




# peak calling with MACS2 for Sample -------------------------------------------
if (length(unique(proj$Sample)) > 1) {

  proj <- addGroupCoverages(
    ArchRProj = proj,
    groupBy = "Sample",
    maxCells = 1500,
    force = TRUE
  )
  # get genome size
  species <- getGenome(ArchRProj = proj)
  if (species == "BSgenome.Hsapiens.UCSC.hg38"){
    genome_size <- 3.3e+09
  } else if (species == "BSgenome.Mmusculus.UCSC.mm10") {
    genome_size = 3.0e+09
  }
  pathToMacs2 <- findMacs2()
  proj <- addReproduciblePeakSet(
    ArchRProj = proj,
    groupBy = "Sample",
    pathToMacs2 = pathToMacs2,
    genomeSize = genome_size,
    maxPeaks = 300000,
    force = TRUE 
  )
  proj <- addPeakMatrix(proj, force = TRUE)

      proj <- addMotifAnnotations(ArchRProj = proj
                                  , motifSet = "cisbp"
                                  , name = "Motif"
                                  , force = TRUE
                                  )

      # save ArchR object
      saveArchRProject(
        ArchRProj = proj,
        outputDirectory = paste0(project_name, "_ArchRProject")
      )
  # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      
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

  motif_lst <- unique(rownames(enrichMotifs))
  split_string <- strsplit(motif_lst, split = "\\(")
  fun1 <- function(list, nth){
    sapply(list, `[` , 1)
  }
  req_motifs3 <- gsub("_","-",fun1(split_string))
  req_motifs3 <- gsub(" ","",req_motifs3)

  rownames(enrichMotifs) <- req_motifs3
  saveRDS(enrichMotifs,"enrichMotifs_sample.rds")

  # cutOff A numeric cutOff that indicates the minimum P-adj enrichment to be included in the heatmap. default is 20 but we decrease that!

  heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 50, transpose = F,returnMatrix = TRUE, cutOff= 2)

  motif_lst <- unique(rownames(heatmapEM))
  split_string <- strsplit(motif_lst, split = "\\(")
  fun1 <- function(list, nth){
    sapply(list, `[` , 1)
  }
  req_motifs3 <- gsub("_","-",fun1(split_string))
  req_motifs3 <- gsub(" ","",req_motifs3)

  rownames(heatmapEM) <- req_motifs3
  write.csv(heatmapEM,"motif_per_sample_hm.csv")


  nSamples = length(unique(proj$Sample))

  df = list()
  tempdir <- "/root"
  motifs_per_sample_hm <- find_func(tempdir,"motif_per_sample_hm.csv")
  hm_per_sample <- read.csv(motifs_per_sample_hm)


  for (i in seq_along(1:nSamples)){
    df[[i]] <- hm_per_sample[,c(1,i+1)]
    
    #select top 5 values by group
    
    df[[i]] <- df[[i]][order(df[[i]][,2], decreasing = T),][1:10,1]
    
  }
  final <- do.call(rbind, df)
  req_motifs3 <- unlist(df)

  req_motifs3<- req_motifs3[!duplicated(req_motifs3)]

  req_motifs3 <- na.omit(req_motifs3)

  write.csv(req_motifs3,"req_motifs3.csv")

  } else {
    enrichMotifs <- "there is not enough samples to be compared with!"  
    heatmapEM <- "there is not enough samples to be compared with!"  
    req_motifs3 <- "there is not enough samples to be compared with!"  
    
  }

  # peak calling with MACS2 for treatment ---------------------------------------

  if (length(unique(proj$Condition))>1){
    for (i in seq_along(treatment)){
    
    proj <- addGroupCoverages(
      ArchRProj = proj,
      groupBy = treatment[i],
      maxCells = 1500,
      force = TRUE
    )
    # get genome size
    species <- getGenome(ArchRProj = proj)
    if (species == "BSgenome.Hsapiens.UCSC.hg38"){
      genome_size <- 3.3e+09
    } else if (species == "BSgenome.Mmusculus.UCSC.mm10") {
      genome_size = 3.0e+09
    }
    pathToMacs2 <- findMacs2()
    proj <- addReproduciblePeakSet(
      ArchRProj = proj,
      groupBy = treatment[i],
      pathToMacs2 = pathToMacs2,
      genomeSize = genome_size,
      maxPeaks = 300000,
      force = TRUE 
    )
    proj <- addPeakMatrix(proj, force = TRUE)
    
    proj <- addMotifAnnotations(ArchRProj = proj
                                , motifSet = "cisbp"
                                , name = "Motif"
                                , force = TRUE
    )
    

    # save ArchR object
    saveArchRProject(
      ArchRProj = proj,
      outputDirectory = paste0(project_name, "_ArchRProject")
    )
    
    }

    # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    for (i in seq_along(treatment)){
      
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
    
    motif_lst <- unique(rownames(enrichMotifs))
    split_string <- strsplit(motif_lst, split = "\\(")
    fun1 <- function(list, nth){
      sapply(list, `[` , 1)
    }
    req_motifs2 <- gsub("_","-",fun1(split_string))
    req_motifs2 <- gsub(" ","",req_motifs2)
    
    rownames(enrichMotifs) <- req_motifs2
    saveRDS(enrichMotifs,paste0("enrichMotifs_treatment_",i,".rds"))
    
    # cutOff A numeric cutOff that indicates the minimum P-adj enrichment to be included in the heatmap. default is 20 but we decrease that!
    
    heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 50, transpose = F,returnMatrix = TRUE, cutOff= 2)
    
    motif_lst <- unique(rownames(heatmapEM))
    split_string <- strsplit(motif_lst, split = "\\(")
    fun1 <- function(list, nth){
      sapply(list, `[` , 1)
    }
    req_motifs2 <- gsub("_","-",fun1(split_string))
    req_motifs2 <- gsub(" ","",req_motifs2)
    
    rownames(heatmapEM) <- req_motifs2
    write.csv(heatmapEM,paste0("motif_per_treatment_",i,"_hm.csv"))
    
    }
    
    
    
    nConds = 2
    
    df = list()
    
    tempdir <- "/root"
    
    motifs_per_cond_hm <- find_func(tempdir,"motif_per_treatment_*")
    
    for (j in seq_along(motifs_per_cond_hm)){
    hm_per_cond <- read.csv(motifs_per_cond_hm[j])
    
    
    for (i in seq_along(1:nConds)){
      df[[i]] <- hm_per_cond[,c(1,i+1)]
      
      #select top 5 values by group
      
      df[[i]] <- df[[i]][order(df[[i]][,2], decreasing = T),][1:5,1]
      
    }
    final <- do.call(rbind, df)
    req_motifs2 <- unlist(df)
    
    req_motifs2 <- req_motifs2[!duplicated(req_motifs2)]
    
    req_motifs2 <- na.omit(req_motifs2)
    
    write.csv(req_motifs2,paste0("req_motifs2_",j,".csv"))
  }
    
  
} else {
  enrichMotifs <- "there is not enough conditions to be compared with!"  
  heatmapEM <- "there is not enough conditions to be compared with!"  
  req_motifs2 <- "there is not enough conditions to be compared with!"  
  
  
}

# Volcano plots for motifs -----------------------------------------------------
if (length(unique(proj$Condition))>1){
  for (j in seq_along(treatment)){
  
  ncells <- length(proj$cellNames)
  
  # all clusters together
  
  markersMotifs <- getMarkerFeatures(
    ArchRProj = proj, 
    useMatrix = "MotifMatrix", 
    groupBy = treatment[j],
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon",
    useSeqnames = "z",maxCells = 5000,
    normBy = "none"
  )
  
  req_DF <- as.data.frame(getCellColData(proj))
  df1 <- table(req_DF$Clusters,req_DF[,treatment[j]])
  distr <- as.data.frame.matrix(round(prop.table(as.matrix(df1),1),2))
  lst <- list()
  
  for(i in 1:nrow(distr)) {
    row <- distr[i,]
    if (
      sum(unname(unlist(row))>= 0.90) == 1) {
      rownames(row) -> lst[[i]]
    }
  }
  not_req_list <- unlist(lst)
  
  req_clusters <- unique(proj$Clusters)
  req_clusters <- req_clusters[order(as.numeric(gsub("C","",req_clusters)))]
  req_clusters <- req_clusters[which(!req_clusters%in%not_req_list)]
  
  markersMotifs_C <- list()
  proj_C <- list()
  
  for (i in seq_along(req_clusters)) {
    
    idxSample <- BiocGenerics::which(proj$Clusters == req_clusters[i])
    
    cellsSample <- proj$cellNames[idxSample]
    proj_C[i] <- proj[cellsSample,]
    
    ncells[i] <- length(proj_C[[i]]$cellNames)
    
    # per each cluster separately
    markersMotifs_C[[i]] <- getMarkerFeatures(
      ArchRProj = proj_C[[i]],
      useMatrix = "MotifMatrix", 
      groupBy = treatment[j],
      bias = c("TSSEnrichment", "log10(nFrags)"),maxCells = ncells[[i]] ,normBy = "none",
      testMethod = "ttest")
  }
  names(markersMotifs_C) <- req_clusters
  
  
  dev_score <- getDeviation_ArchR(ArchRProj = proj, name = motifs
                                  , imputeWeights = getImputeWeights(proj))
  
  # which(rowSums(is.na(dev_score))>0)
  # any(rowSums((dev_score))==0)
  
  empty_motif_idx <- which(rowSums((dev_score))==0)
  empty_motif <- rownames(dev_score)[empty_motif_idx] 
  
  
  markersMotifs_df1 <- assay(markersMotifs, "MeanDiff")
  markersMotifs_df2 <- assay(markersMotifs, "Pval")
  markersMotifs_df3 <- assay(markersMotifs, "FDR")
  markersMotifs_df <- cbind(markersMotifs_df1,markersMotifs_df2,markersMotifs_df3)
  markersMotifs_df$genes<- rowData(markersMotifs)$name
  markersMotifs_df$cluster <- rep("All",length(rownames(markersMotifs_df)))
  
  # # we only want to see results of one set , say just sham
  markersMotifs_df <- markersMotifs_df[,c(1,3,5,7,8)]
  colnames(markersMotifs_df) <- c("avg_log2FC","p_val","p_val_adj","gene","cluster")
  
  
  # percluster
  
  markersMotifs_df1_C <- list()
  markersMotifs_df2_C <- list()
  markersMotifs_df3_C <- list()
  markersMotifs_df_C <- list()
  
  # for (i in (1:nClust)){
  for (i in seq_along(req_clusters)){
    
    cluster <- req_clusters[i]
    
    markersMotifs_df1_C[[i]] <- assay(markersMotifs_C[[i]], "MeanDiff")
    markersMotifs_df2_C[[i]] <- assay(markersMotifs_C[[i]], "Pval")
    markersMotifs_df3_C[[i]] <- assay(markersMotifs_C[[i]], "FDR")
    markersMotifs_df_C[[i]] <- cbind(markersMotifs_df1_C[[i]],markersMotifs_df2_C[[i]],markersMotifs_df3_C[[i]])
    markersMotifs_df_C[[i]]$genes<- rowData(markersMotifs_C[[i]])$name
    markersMotifs_df_C[[i]]$cluster <- rep(cluster,length(rownames(markersMotifs_df_C[[i]])))
    
    # # we only want to see results of one set , say just sham
    markersMotifs_df_C[[i]] <- markersMotifs_df_C[[i]][,c(1,3,5,7,8)]
    colnames(markersMotifs_df_C[[i]]) <- c("avg_log2FC","p_val","p_val_adj","gene","cluster")
  }
  
  
  names(markersMotifs_df_C) <- req_clusters
  
  #### merge all data frames
  
  markersMotifs_merged_df <- do.call("rbind", markersMotifs_df_C)
  
  # also data frame for all clusters together needs to be added
  
  markersMotifs_merged_df <- rbind(markersMotifs_df,markersMotifs_merged_df)  
  # remove empty genes
  markersMotifs_merged_df <- markersMotifs_merged_df[which(!markersMotifs_merged_df$gene%in%empty_motif),]
  
  # remove na values
  markersMotifs_merged_df <- na.omit(markersMotifs_merged_df)
  
  # remove FDR equal to 0
  markersMotifs_merged_df <- markersMotifs_merged_df[which(!markersMotifs_merged_df$p_val_adj== 0),]
  
  
  # make logfc limiation between 1 and -1
  
  markersMotifs_merged_df <- markersMotifs_merged_df[which(abs(markersMotifs_merged_df$avg_log2FC)< 1.2),]
  
  markersMotifs_merged_df$Significance = ifelse(markersMotifs_merged_df$p_val_adj < 10^-1 , 
                                                ifelse(markersMotifs_merged_df$avg_log2FC> 0.0 
                                                       ,colnames(markersMotifs)[1],colnames(markersMotifs)[2]),
                                                'Not siginficant')
  
  de <- markersMotifs_merged_df
  print(paste0("volcanoMarkers_motifs_", j, ".txt is writing!"))
  write.table(
    de,
    paste0("volcanoMarkers_motifs_",j, ".txt"),
    sep = '\t',
    quote = FALSE,
    row.names = FALSE
  )
  print(paste0("writing volcanoMarkers_motifs_", j, ".txt is done!"))
  }

  features_m <- unique(de$cluster)
  volcano_plots_m <- list()
  for (i in seq_along(features_m)) {
    volcano_plots_m[[i]] <- scvolcano(de, markersMotifs, features_m[[i]])
  }

  pdf("volcano_plots_motifs.pdf")
  for (plot in volcano_plots_m) {
    print(plot)
  }
  dev.off()
  
} else {
  
  de <- "there is not enough conditions to be compared with!"  
}   




# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# save ArchR object
saveArchRProject(
  ArchRProj = proj,
  outputDirectory = paste0(project_name, "_ArchRProject")
)

#################------------- Motif Logo ---------------#######################


#data("human_pwms_v1")

PWMs <- getPeakAnnotation(proj, "Motif")$motifs

PWMatrixToProbMatrix <- function(x){
  if (class(x) != "PWMatrix") stop("x must be a TFBSTools::PWMatrix object")
  (exp(as(x, "matrix"))) * TFBSTools::bg(x)/sum(TFBSTools::bg(x))
}

ProbMatrices <- lapply(PWMs, PWMatrixToProbMatrix)
lapply(ProbMatrices, colSums) %>% range


PWMatrixToProbMatrix <- function(x){
  if (class(x) != "PWMatrix") stop("x must be a TFBSTools::PWMatrix object")
  m <- (exp(as(x, "matrix"))) * TFBSTools::bg(x)/sum(TFBSTools::bg(x))
  m <- t(t(m)/colSums(m))
  m
}

ProbMatrices <- lapply(PWMs, PWMatrixToProbMatrix)
lapply(ProbMatrices, colSums) %>% range


saveRDS(ProbMatrices,"seqlogo.rds")

######creating combimed.rds files ########

main_func <- function(seurat_lst, umap_embedding) {
  
  find_samples_name <- function(seurat_lst) {
    # Extract list of sample names from list of SeuratObjs.
    sapply(seq_along(seurat_lst), function(i) {
      unique(seurat_lst[[i]]@meta.data$Sample)
    })
  }
  
  samples <- find_samples_name(seurat_lst)
  
  D00_fun <- function(seurat_lst) {
    # Remove samples without "counts" from list of SeuratObjs
    toRemove <- lapply(seurat_lst, function(x) {
      names(which(colSums(is.na(x@assays[[1]]@counts)) > 0))
      }) 
    mapply(function(x, y) x[, !colnames(x) %in% y], seurat_lst, toRemove)
  }
  
  D00 <- D00_fun(seurat_lst)
  
  Spatial_D00_fun <- function(D00){
    
    Spatial_D00 <- lapply(D00, function(x) {
      as.data.frame(x@images[[1]]@coordinates[,c(5,4)])
    })
    Spatial_D00 <- lapply(Spatial_D00, function(x) {
      colnames(x) <- paste0("Spatial_", 1:2)
      x
    })  
    lapply(Spatial_D00, function(x) {
      x$Spatial_2 <- -(x$Spatial_2) 
      x
    })
  }
  
  Spatial_D00 <- Spatial_D00_fun(D00)
  
  
  Spatial_D00_all_fun <- function(Spatial_D00) {
    
    tmp <- lapply(seq_along(Spatial_D00), function(i) {
      bind_rows(Spatial_D00[-i])
    })
    tmp <- lapply(tmp, function(x) {
      x$Spatial_1<- 0
      x
    })
    tmp <- lapply(tmp, function(x) {
      x$Spatial_2<- 0
      x
    })
    tmp <- lapply(seq_along(Spatial_D00), function(i) {
      as.matrix(rbind(Spatial_D00[[i]], tmp[[i]]))
    })
  }
  
  Spatial_D00_all <- Spatial_D00_all_fun(Spatial_D00)
  
  temp_fun <- function(D00){
    # Convert list of SeuratObjs to list of Assay counts as dataframes.
    temp <- lapply(D00, function(x) {
      df <- as.data.frame(x@assays[[1]]@counts)
      colnames(df) <- Cells(x)
      return(df)a
    })
    temp <- lapply(temp, function(x) {
      x$region <- rownames(x)
      x
    })  
  }
  
  temp <- temp_fun(D00)
  
  # merge seurat objects
  combined_mat <- purrr::reduce(temp, full_join, by = "region")
  
  rownames(combined_mat) <- combined_mat$region
  combined_mat$region<- NULL
  
  # remove extra cells
  extra_cells <- setdiff(colnames(combined_mat), rownames(Spatial_D00_all[[1]]))
  combined_mat <- combined_mat[, which(!colnames(combined_mat) %in% extra_cells)]
  combined_mat <- as.matrix(combined_mat)
  
  # clean columns of metadata per sample that attached sample's name before rbind
  l <- D00
  l <- lapply(l, function(x) {
    colnames(x@meta.data) <- gsub(paste0("_", Assays(x)), "", colnames(x@meta.data))
    x
  })
  D00 <- l
  
  # first get the list of meta data
  list_of_metadata <- lapply(D00, function(x) {
    x@meta.data
  })
  
  # # rbind meta data per samples
  meta.data <- do.call("rbind", list_of_metadata)
  write.csv(meta.data, "req_meta_data.csv", row.names = TRUE)
  
  combined <- CreateSeuratObject(
    counts = as.data.frame(combined_mat),
    assay = "scATAC",
    meta.data = meta.data
  )
  
  combined@meta.data$Clusters <- factor(
    combined@meta.data$Clusters,
    levels = c(paste0("C", seq_along(unique(combined@meta.data$Clusters))))
  )
  
  Spatial_D00 <- list()
  for (i in seq_along(samples)) {
    Spatial_D00[[i]] <- Spatial_D00_all[[i]][colnames(combined), ]
    combined[[samples[i]]] <- CreateDimReducObject(
                                embeddings = Spatial_D00[[i]],
                                key = samples[i],
                                assay = DefaultAssay(combined)
                              )
    }
  
  
  # we need to run Variable Features
  combined <- NormalizeData(
    combined, normalization.method = "LogNormalize", scale.factor = 10000
  )
  combined <- FindVariableFeatures(
    combined, selection.method = "vst", nfeatures = 2000
  )
  
  combined[["UMAP"]] <- CreateDimReducObject(
                          embeddings =  as.matrix(umap_embedding),
                          key = "UMAP",
                          assay = DefaultAssay(combined)
                        )
  
  return(combined)
}

combined <- main_func(all, UMAPHarmony)
combined_m <- main_func(all_m, UMAPHarmony)

saveRDS(combined, "combined.rds", compress = FALSE)
saveRDS(combined_m, "combined_m.rds", compress = FALSE)

# =================================
print('Shiny App starting ...')

scConf1 = createConfig(combined)
makeShinyFiles(combined, scConf1, gex.assay = "scATAC", gex.slot = "data",
               gene.mapping = TRUE, shiny.prefix = "sc1",
               default.gene1 = "Tiam1", default.gene2 = "Ccbe1",
               default.dimred = c("UMAP_1", "UMAP_2"))




scConf2 = createConfig(combined_m)
makeShinyFiles(combined_m, scConf2, gex.assay = "scATAC", gex.slot = "counts",
               gene.mapping = TRUE, shiny.prefix = "sc2",
               default.gene1 = "RFX3-1018", default.gene2 = "NEUROG2-1580",
               default.dimred = c("UMAP_1", "UMAP_2"))


citation = list(
  title   = paste0(project_name," Data Analysis")
)
makeShinyCodesMulti(
  shiny.title = paste0(project_name,"_Lab Data Analysis"), 
  shiny.footnotes = citation,
  shiny.prefix = c("sc1", "sc2"),
  
  shiny.headers = c("Gene Accessibility", "Peak/Motifs"), 
  shiny.dir = "./shinyApp"
) 

# edit some of the prepared data


sc1def <- readRDS("/root/shinyApp/sc1def.rds")
sc2def <- readRDS("/root/shinyApp/sc2def.rds")

find_samples_name <- function(lst){
  
  sapply(seq_along(lst), function(i) unique(lst[[i]]@meta.data$Sample))
}   
samples <- find_samples_name(all)


D00 <- list()
for (i in seq_along(samples)) {
  D00[[i]] <- all[[i]]
  nal_cols <- which(colSums(is.na(D00[[i]]@assays[[1]]@counts))>0)
  toRemove <- names(nal_cols)
  D00[[i]] <- D00[[i]][,!colnames(D00[[i]]) %in% toRemove]
}

Spatial_D00 <- list()
for (i in seq_along(samples)) {
  Spatial_D00[[i]] <- as.data.frame(D00[[i]]@images[[1]]@coordinates[,c(5,4)])
  colnames(Spatial_D00[[i]]) <- paste0("Spatial_", 1:2)
  Spatial_D00[[i]]$Spatial_2 <- -(Spatial_D00[[i]]$Spatial_2)
}

l <- Spatial_D00
xlim <- lapply(l, function(x) { xlim <- c(min(x[,1]),max(x[,1]));xlim})
ylim <- lapply(l, function(y) { ylim <- c(min(y[,2]),max(y[,2]));ylim})


sc1def$limits <- list()
for (i in seq_along(samples)) {
  sc1def[['limits']][[samples[i]]]<- c(
    min(xlim[[i]]),max(xlim[[i]])
    ,min(ylim[[i]]),max(ylim[[i]])
  )}

sc2def$limits <- list()
for (i in seq_along(samples)) {
  sc2def[['limits']][[samples[i]]]<- c(
    min(xlim[[i]]),max(xlim[[i]])
    ,min(ylim[[i]]),max(ylim[[i]])
  )}

sc1def$meta1<- "Clusters"

sc1def$meta2<- "Sample"

sc2def$meta1<- "Clusters"

sc2def$meta2<- "Sample"


sc1def$Clusters <- req_genes1

sc1def$Sample <- req_genes3

sc2def$Clusters <- req_motifs1

sc2def$Sample <- req_motifs3


sc1def$dimred[3] <- paste0(names(combined@reductions)[1],"_1")
sc1def$dimred[4] <- paste0(names(combined@reductions)[1],"_2")
sc1def$dimred[5] <- paste0(names(combined@reductions)[2],"_1")
sc1def$dimred[6] <- paste0(names(combined@reductions)[2],"_2")


sc2def$dimred[3] <- paste0(names(combined_m@reductions)[1],"_1")
sc2def$dimred[4] <- paste0(names(combined_m@reductions)[1],"_2")
sc2def$dimred[5] <- paste0(names(combined_m@reductions)[2],"_1")
sc2def$dimred[6] <- paste0(names(combined_m@reductions)[2],"_2")


treatment <- names(getCellColData(proj))[grep('condition_',names(getCellColData(proj)))]

if (length(unique(proj$Condition))>1) {
for(i in seq_along(treatment)){
  
  sc1def[[paste0('meta',(2+i))]] <- treatment[i]
  sc1def[[treatment[i]]] <- read.csv(find_func(tempdir, paste0('req_genes2_',i,'.csv')))$x
  sc1def[[paste0(treatment[i],"_1")]] <- sort(unique(combined@meta.data[[treatment[i]]]))[1]
  sc1def[[paste0(treatment[i],"_2")]] <- sort(unique(combined@meta.data[[treatment[i]]]))[2]
  
  
  sc2def[[paste0('meta',(2+i))]] <- treatment[i]
  sc2def[[treatment[i]]] <- read.csv(find_func(tempdir, paste0('req_motifs2_',i,'.csv')))$x
  sc2def[[paste0(treatment[i],"_1")]] <- sort(unique(combined_m@meta.data[[treatment[i]]]))[1]
  sc2def[[paste0(treatment[i],"_2")]] <- sort(unique(combined_m@meta.data[[treatment[i]]]))[2]
  
}
}else{
  
  print('there is not enough conditions to add to sc1def.rds or sc2def.rds')
}

saveRDS(sc1def,"/root/shinyApp/sc1def.rds")
saveRDS(sc2def,"/root/shinyApp/sc2def.rds")

sc1conf <- readRDS("/root/shinyApp/sc1conf.rds")
sc2conf <- readRDS("/root/shinyApp/sc2conf.rds")

fav <- which(sc1conf$ID %in% c("Clusters","Sample", treatment))
rest <- which(!sc1conf$ID %in% c("Clusters","Sample", treatment))
sc1conf <- sc1conf[c(fav,rest),]


fav <- which(sc2conf$ID %in% c("Clusters","Sample", treatment))
rest <- which(!sc2conf$ID %in% c("Clusters","Sample", treatment))
sc2conf <- sc2conf[c(fav,rest),]

saveRDS(sc1conf,"/root/shinyApp/sc1conf.rds")
saveRDS(sc2conf,"/root/shinyApp/sc2conf.rds")


  
  if (length(unique(proj$Sample)) <=1){
    file.remove("/root/ui.R")
    file.remove("/root/server.R")
    file.remove("/root/ui_2.R")
    file.remove("/root/server_2.R")
    file.rename("/root/ui_3.R","/root/ui.R")
    file.rename("/root/server_3.R","/root/server.R")
    
  } else if (length(unique(proj$Condition))<=1) {
    
    file.remove("/root/ui.R")
    file.remove("/root/server.R")
    file.remove("/root/ui_3.R")
    file.remove("/root/server_3.R")
    file.rename("/root/ui_2.R","/root/ui.R")
    file.rename("/root/server_2.R","/root/server.R")
    
    
  }else{
    file.remove("/root/ui_2.R")
    file.remove("/root/server_2.R")
    file.remove("/root/ui_3.R")
    file.remove("/root/server_3.R")
    
    
    
  }

# copy everything in shiny app to the root

rawPath <- "/root/shinyApp/"
dataPath <- "/root/"

dataFiles <- dir(rawPath, "*.rds$", ignore.case = TRUE, all.files = TRUE)
file.copy(file.path(rawPath, dataFiles), dataPath, overwrite = TRUE)

dataFiles <- dir(rawPath, "*.h5$", ignore.case = TRUE, all.files = TRUE)
file.copy(file.path(rawPath, dataFiles), dataPath, overwrite = TRUE)

  
  