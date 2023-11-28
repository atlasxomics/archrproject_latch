library("ArchR")
library("BSgenome")
library("BSgenome.Hsapiens.UCSC.hg38")
library("BSgenome.Mmusculus.UCSC.mm10")
library("gridExtra")
library("harmony")
library("Seurat")
library("GenomicRanges")
library("dplyr")                                    # Load dplyr package
library("plyr")                                     # Load plyr package
library("readr")
library("qdap")
library("ShinyCell")
require("tidyverse")
library("seqLogo")
require("ggseqlogo")
library("chromVARmotifs")
library("ComplexHeatmap")
library("circlize")
library(data.table)

source("/root/getDeviation_ArchR.R")
source("/root/wf/utils.R")



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

runs <- strsplit(args[11:length(args)], ",")
inputs <- c()
for (run in runs) {
  inputs[run[1]] <- run[2]
}

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

spatial_cluster_plots <- list()
for (i in seq_along(run_ids)){
  plot <- plot_spatial(seurat_objs[[i]], run_ids[i])
  spatial_cluster_plots[[i]] <- plot
}

pdf("spatial_plots.pdf")
wrap_plots(spatial_cluster_plots, guides = "collect") &
  theme(
    legend.position = "bottom",
    text = element_text(size = 12),
    plot.title = element_text(size = 10)
  )
dev.off()

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

pdf("qc_plots.pdf")
for (i in seq_along(metrics)) {
  grid.arrange(
    grobs = all_qc_plots[[i]],
    ncol = 2
  )
}
dev.off()


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
  heatmapGS,
  heatmap_legend_side = "bot",
  annotation_legend_side = "bot",
  column_title = paste0("Marker genes ()", gene_cutoff, ")"),
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
    
    print(paste0("inpMarkers_motif_",j,".txt is writing!"))
    
    write.table(de,paste0("inpMarkers_",j,".txt"), sep = '\t', quote = F, row.names = F)
    
    print(paste0("writing inpMarkers_",j,".txt is done!"))
    
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



UMAPHarmony <-getEmbedding(ArchRProj = proj, embedding = "UMAP", returnDF = TRUE)
write.csv(UMAPHarmony,"UMAPHarmony.csv")


# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# save ArchR object
saveArchRProject(
  ArchRProj = proj,
  outputDirectory = paste0(project_name, "_ArchRProject")
)
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


# peak calling with MACS2 for clusters-----------------------------------------------------

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

markersPeaks <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  k = 100,
  testMethod = "wilcoxon"
)

######################### Plot marker peaks, motifs ###########################

hpeak_cutoff <- "FDR <= 0.05 & Log2FC >= 1"
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

heatmapEM <- plotEnrichHeatmap(enrichMotifs, transpose = TRUE)

heatmap_motifs <- ComplexHeatmap::draw(
  heatmapEM,
  heatmap_legend_side = "bottom",
  column_title = paste0("Marker motifs (", motifs_cutoff),
  column_title_gp = gpar(fontsize = 12)
)

heatmaps[[3]] <- heatmap_motifs

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





# peak calling with MACS2 for Sample -------------------------------------------
if (length(unique(proj$Sample))>1){
  
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
  names(markerList_C) <- req_clusters
  
  
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
  print(paste0("inpMarkers_motif_",j,".txt is writing!"))
  write.table(de,paste0("inpMarkers_motif_",j,".txt"), sep = '\t', quote = F, row.names = F)
  print(paste0("writing inpMarkers_motif_",j,".txt is done!"))
  }
  
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


data("human_pwms_v1")

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

