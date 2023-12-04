setwd("~/latch/archrproject/wf/")

library(ArchR)
library(ggplot2)
library(ggrepel)
library(gridExtra)

proj <- loadArchRProject("jm_test_ArchRProject/")
run_ids <- unique(proj$Sample)

#######################some treatment stuff from Noori##########################

conds <- strsplit(proj$Condition, split = "\\s|-")
fun1 <- function(lst, n) {
  sapply(lst, "[", n)
}

for (i in seq_along(conds[[1]])) {
  proj <- addCellColData(
    proj,
    data = fun1(conds, i),
    name = paste0('condition_',i),
    cells = proj$cellNames,
    force = TRUE
  )
}

treatment <- names(getCellColData(proj))[
  grep('condition_', names(getCellColData(proj)))
]


#######################make volcanos############################################

if (length(unique(proj$Condition)) > 1) {
  for (j in seq_along(treatment)) {
    
    ncells <- length(proj$cellNames)
    
    markerList <- getMarkerFeatures(
      ArchRProj = proj,
      useMatrix = "GeneScoreMatrix",
      groupBy = treatment[j],
      bias = c("TSSEnrichment", "log10(nFrags)"),
      maxCells = ncells,
      normBy = "none",
      testMethod = "ttest"
    )
    
    req_DF <- as.data.frame(getCellColData(proj))
    df1 <- table(req_DF$Clusters, req_DF[, treatment[j]])
    distr <- as.data.frame.matrix(round(prop.table(as.matrix(df1), 1), 2))
    lst <- list()
    
    for(i in 1:nrow(distr)) {
      row <- distr[i, ]
      if (
        sum(unname(unlist(row)) >= 0.90) == 1) {
        rownames(row) -> lst[[i]]
      }
    }
    not_req_list <- unlist(lst)
    
    req_clusters <- unique(proj$Clusters)
    req_clusters <- req_clusters[order(as.numeric(gsub("C", "", req_clusters)))]
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
        bias = c("TSSEnrichment", "log10(nFrags)"),
        maxCells = ncells[[i]],
        normBy = "none",
        testMethod = "ttest"
      )
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
    
    markerList_df1_C <- list()
    markerList_df2_C <- list()
    markerList_df3_C <- list()
    markerList_df_C <- list()
    
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
    
    # gene activity score by cluster
    markersGS_merged_df <- do.call("rbind", markerList_df_C)
    
    # also data frame for all clusters together needs to be added
    markersGS_merged_df <- rbind(markerList_df, markersGS_merged_df) # this overwrites the clusters
    
    # remove empty genes
    markersGS_merged_df <- markersGS_merged_df[which(!markersGS_merged_df$gene%in%empty_gene),]
    
    # remove na values
    markersGS_merged_df <- na.omit(markersGS_merged_df)
    
    # remove FDR equal to 0
    markersGS_merged_df <- markersGS_merged_df[which(!markersGS_merged_df$p_val_adj== 0),]
    
    # make logfc limiation between 1 and -1
    
    markersGS_merged_df <- markersGS_merged_df[which(abs(markersGS_merged_df$avg_log2FC)< 1.2),]
    
    markersGS_merged_df$Significance = ifelse(
      markersGS_merged_df$p_val < 10^-2,
      ifelse(markersGS_merged_df$avg_log2FC > 0.0,
             colnames(markerList)[1],
             colnames(markerList)[2]),
      'Not siginficant')
    
    de <- markersGS_merged_df
    
    print(paste0("inpMarkers_", j, ".txt is writing!"))
    
    write.table(de,paste0("inpMarkers_",j,".txt"), sep = '\t', quote = F, row.names = F)
    
    print(paste0("writing inpMarkers_",j,".txt is done!"))
    
  }
  
} else {
  
  de <- "there is not enough conditions to be compared with!"  
}   

################################## plotting ####################################

sctheme <- function(base_size = 15, XYval = TRUE, Xang = 0, XjusH = 0.5) { 
  oupTheme = theme( 
    text = element_text(size = base_size, family = "Helvetica"), 
    panel.background = element_rect(fill = "white", colour = NA), 
    axis.line = element_line(colour = "black"), 
    axis.ticks = element_line(colour = "black", size = base_size / 20), 
    axis.title = element_text(face = "bold"), 
    axis.text = element_text(size = base_size), 
    axis.text.x = element_text(angle = Xang, hjust = XjusH), 
    legend.position = "bottom", 
    legend.key = element_rect(colour = NA, fill = NA) 
  ) 
  if(!XYval){ 
    oupTheme = oupTheme + theme( 
      axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
      axis.text.y = element_blank(), axis.ticks.y = element_blank()) 
  } 
  return(oupTheme) 
} 

scvolcano <- function(inpMarkers, feature = "All"){
  
  # Prepare ggData
  ggData <- inpMarkers[which(inpMarkers$cluster==feature),]
  minfdr = 0.09
  minfdr1 = 10^-(1/6 *(-log10(min(ggData$p_val_adj))))
  
  minfdr2 = 10^-(2/3 *(-log10(min(ggData$p_val_adj))))
  
  
  ggData$Significance = ifelse(
    ggData$p_val_adj < minfdr,
    ifelse(ggData$avg_log2FC > 0.0,
           markerList@colData@rownames[[1]],
           markerList@colData@rownames[[2]]
           ),
    'Not siginficant'
  )
  
  ggData$Significance <- factor(
    ggData$Significance,
    levels = c(
      markerList@colData@rownames[[1]],
      markerList@colData@rownames[[2]],
      'Not siginficant')
  )
  
  ggData[ggData$p_val_adj < 1e-300, "p_val_adj"] <- 1e-300
  ggData$log10fdr = -log10(ggData$p_val_adj)
  
  # Actual ggplot
  ggOut <-
    ggplot(ggData, aes(avg_log2FC, log10fdr)) +
    geom_point() +
    sctheme() +
    ylab("-log10(FDR)") +
    geom_point(aes(color = Significance)) +
    scale_color_manual(values = c("#F8766D","#619CFF","gray")) +
    geom_text_repel(
      data = subset(ggData, p_val_adj < minfdr1 ),
      aes(label = gene)) +
    ggtitle(paste("Marker genes:", feature)) +
    theme(plot.title = element_text(hjust = 0.5, size = 20))
  return(ggOut)
}

features <- unique(de$cluster)
volcano_plots <- list()
for (i in seq_along(features)) {
  volcano_plots[[i]] <- scvolcano(de, features[[i]])
}

pdf("volcano_plots.pdf")
for (plot in volcano_plots) {
  print(plot)
}
dev.off()

