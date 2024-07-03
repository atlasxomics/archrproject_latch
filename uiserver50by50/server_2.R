# this version produced to use when the project has only one or None condition!

library(shiny)
library(shinyhelper)
library(data.table)
library(Matrix)
library(DT)
library(magrittr)
library(ggplot2)
library(ggrepel)
library(hdf5r)
library(ggdendro)
library(gridExtra)
require(ggseqlogo)# Motif Logo
library("ComplexHeatmap") #heatmap
library("circlize") #heatmap
library("S4Vectors")
library("SummarizedExperiment")
library("ggplotify")
library("ArchR")
library("Seurat")
library("stringi")
library(RColorBrewer)

tempdir <- getwd()
ArchRobj <- system(paste0("find ", tempdir, " -name '*_ArchRProject' -type d"), intern = TRUE)


proj <- loadArchRProject(path = ArchRobj, force = FALSE, showLogo = TRUE)


sc1conf = readRDS("sc1conf.rds")
sc1def  = readRDS("sc1def.rds")
sc1gene = readRDS("sc1gene.rds")
sc1meta = readRDS("sc1meta.rds")

sc2conf = readRDS("sc2conf.rds")
sc2def  = readRDS("sc2def.rds")
sc2gene = readRDS("sc2gene.rds")
sc2meta = readRDS("sc2meta.rds")

# for genes heatmap
seMarker_cluster <-  readRDS("markersGS_clusters.rds")
seMarker_sample <-  readRDS("markersGS_sample.rds")

# for motif heatmap
seEnrich_cluster <-  readRDS("enrichMotifs_clusters.rds")
seEnrich_sample <-  readRDS("enrichMotifs_sample.rds")

# for conditions
treatment <- names(getCellColData(proj))[grep('condition_',names(getCellColData(proj)))]

combined <- readRDS('combined.rds')

### Useful stuff
# Colour palette 
cList = list(c("grey85","#FFF7EC","#FEE8C8","#FDD49E","#FDBB84", 
               "#FC8D59","#EF6548","#D7301F","#B30000","#7F0000"), 
             c("#4575B4","#74ADD1","#ABD9E9","#E0F3F8","#FFFFBF", 
               "#FEE090","#FDAE61","#F46D43","#D73027")[c(1,1:9,9)], 
             c("#FDE725","#AADC32","#5DC863","#27AD81","#21908C", 
               "#2C728E","#3B528B","#472D7B","#440154"),
             c("#E6E7E8","#3A97FF","#8816A7","black"),
             c('#352A86','#343DAE','#0262E0','#1389D2','#2DB7A3','#A5BE6A','#F8BA43','#F6DA23','#F8FA0D')
             
             ) 
names(cList) = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple","comet","blueYellow") 
 
# Panel sizes 
pList = c("400px", "600px", "800px") 
names(pList) = c("Small", "Medium", "Large") 
pList2 = c("500px", "700px", "900px") 
names(pList2) = c("Small", "Medium", "Large") 
pList3 = c("600px", "800px", "1000px") 
names(pList3) = c("Small", "Medium", "Large") 
sList = c(18,24,30) 
names(sList) = c("Small", "Medium", "Large") 
lList = c(5,6,7) 
names(lList) = c("Small", "Medium", "Large") 
 
# Function to extract legend 
g_legend <- function(a.gplot){  
  tmp <- ggplot_gtable(ggplot_build(a.gplot))  
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")  
  legend <- tmp$grobs[[leg]]  
  legend 
}  
 
# Plot theme 
sctheme <- function(base_size = 24, XYval = TRUE, Xang = 0, XjusH = 0.5){ 
  oupTheme = theme( 
    text =             element_text(size = base_size, family = "Helvetica"), 
    panel.background = element_rect(fill = "white", colour = NA), 
    axis.line =   element_line(colour = "black"), 
    axis.ticks =  element_line(colour = "black", size = base_size / 20), 
    axis.title =  element_text(face = "bold"), 
    axis.text =   element_text(size = base_size), 
    axis.text.x = element_text(angle = Xang, hjust = XjusH), 
    legend.position = "bottom", 
    legend.key =      element_rect(colour = NA, fill = NA) 
  ) 
  if(!XYval){ 
    oupTheme = oupTheme + theme( 
      axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
      axis.text.y = element_blank(), axis.ticks.y = element_blank()) 
  } 
  return(oupTheme) 
} 
 
### Common plotting functions 
# Plot cell information on dimred 
scDRcell <- function(inpConf, inpMeta, inpdrX, inp1, inpsub1, inpsub2, inpsub3, 
                     inpsiz, inpcol, inpord, inpfsz, inpasp, inptxt, inplab){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Prepare ggData 
  inpdrY = paste0(substring(inpdrX, 1, nchar(inpdrX)-1),2)
  
  ggData = inpMeta[, c(inpConf[UI == inpdrX]$ID, inpConf[UI == inpdrY]$ID, 
                       inpConf[UI == inp1]$ID, inpConf[UI == inpsub1]$ID),  
                   with = FALSE] 
  colnames(ggData) = c("X", "Y", "val", "sub") 
  rat = (max(ggData$X) - min(ggData$X)) / (max(ggData$Y) - min(ggData$Y)) 
  bgCells = FALSE 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    bgCells = TRUE 
    ggData2 = ggData[!sub %in% inpsub2] 
    ggData = ggData[sub %in% inpsub2] 
  } 
  if(inpord == "Max-1st"){ 
    ggData = ggData[order(val)] 
  } else if(inpord == "Min-1st"){ 
    ggData = ggData[order(-val)] 
  } else if(inpord == "Random"){ 
    ggData = ggData[sample(nrow(ggData))] 
  } 
  
  # Do factoring if required 
  if(!is.na(inpConf[UI == inp1]$fCL)){ 
    # ggCol = strsplit(inpConf[UI == inp1]$fCL, "\\|")[[1]] 
    ggCol <- inpsub3
    names(ggCol) = levels(ggData$val) 
    ggLvl = levels(ggData$val)[levels(ggData$val) %in% unique(ggData$val)] 
    ggData$val = factor(ggData$val, levels = ggLvl) 
    ggCol = ggCol[ggLvl] 
  } 
 
  # Actual ggplot 
  ggOut = ggplot(ggData, aes(X, Y, color = val)) 
  if(bgCells){ 
    ggOut = ggOut + 
      geom_point(data = ggData2, color = "snow2", size = inpsiz, shape = 16) 
  } 
  ggOut = ggOut + 
    geom_point(size = inpsiz, shape = 16) + xlab(inpdrX) + ylab(inpdrY) + 
    sctheme(base_size = sList[inpfsz], XYval = inptxt) 
  if(is.na(inpConf[UI == inp1]$fCL)){ 
    ggOut = ggOut + scale_color_gradientn("", colours = cList[[inpcol]]) + 
      guides(color = guide_colorbar(barwidth = 15)) 
  } else { 
    sListX = min(nchar(paste0(levels(ggData$val), collapse = "")), 200) 
    sListX = 0.75 * (sList - (1.5 * floor(sListX/50))) 
    ggOut = ggOut + scale_color_manual("", values = ggCol) + 
      guides(color = guide_legend(override.aes = list(size = 5),  
                                  nrow = inpConf[UI == inp1]$fRow)) + 
      theme(legend.text = element_text(size = sListX[inpfsz])) 
    if(inplab){ 
      ggData3 = ggData[, .(X = mean(X), Y = mean(Y)), by = "val"] 
      lListX = min(nchar(paste0(ggData3$val, collapse = "")), 200) 
      lListX = lList - (0.25 * floor(lListX/50)) 
      # ggOut = ggOut + 
      #   geom_text_repel(data = ggData3, aes(X, Y, label = val), 
      #                   color = "grey10", bg.color = "grey95", bg.r = 0.15, 
      #                   size = lListX[inpfsz], seed = 42) 
    } 
  } 
  if(inpasp == "Square") { 
    ggOut = ggOut + coord_fixed(ratio = rat)
  } else if(inpasp == "Fixed") { 
    ggOut = ggOut + coord_fixed()
  }
  
    
    p <- ggOut

      for(i in names(sc1def$limits)){
        
        reduction <- substring(inpdrX, 1, nchar(inpdrX)-2)
        if(reduction == i) {
          

          p <- p+ xlim(sc1def$limits[[i]][1],sc1def$limits[[i]][2])+ ylim(sc1def$limits[[i]][3],sc1def$limits[[i]][4])
          
        }  else  {
          
          p <- p 
        }
      
      }
    
  

  return(p)
}
 
scDRnum <- function(inpConf, inpMeta, inp1, inp2, inpsub1, inpsub2, 
                    inpH5, inpGene, inpsplt){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Prepare ggData 
  ggData = inpMeta[, c(inpConf[UI == inp1]$ID, inpConf[UI == inpsub1]$ID), 
                   with = FALSE] 
  colnames(ggData) = c("group", "sub") 
  h5file <- H5File$new(inpH5, mode = "r") 
  h5data <- h5file[["grp"]][["data"]] 
  ggData$val2 = h5data$read(args = list(inpGene[inp2], quote(expr=))) 
  ggData[val2 < 0]$val2 = 0 
  h5file$close_all() 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    ggData = ggData[sub %in% inpsub2] 
  } 
  
  # Split inp1 if necessary 
  if(is.na(inpConf[UI == inp1]$fCL)){ 
    if(inpsplt == "Quartile"){nBk = 4} 
    if(inpsplt == "Decile"){nBk = 10} 
    ggData$group = cut(ggData$group, breaks = nBk) 
  } 
  
  # Actual data.table 
  ggData$express = FALSE 
  ggData[val2 > 0]$express = TRUE 
  ggData1 = ggData[express == TRUE, .(nExpress = .N), by = "group"] 
  ggData = ggData[, .(nCells = .N), by = "group"] 
  ggData = ggData1[ggData, on = "group"] 
  ggData = ggData[, c("group", "nCells", "nExpress"), with = FALSE] 
  ggData[is.na(nExpress)]$nExpress = 0 
  ggData$pctExpress = 100 * ggData$nExpress / ggData$nCells 
  ggData = ggData[order(group)] 
  colnames(ggData)[3] = paste0(colnames(ggData)[3], "_", inp2) 
  return(ggData) 
} 
# Plot gene expression on dimred 
scDRgene <- function(inpConf, inpMeta, inpdrX, inp1, inpsub1, inpsub2, 
                     inpH5, inpGene, 
                     inpsiz, inpcol, inpord, inpfsz, inpasp, inptxt){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Prepare ggData 
  inpdrY = paste0(substring(inpdrX, 1, nchar(inpdrX)-1),2)
  ggData = inpMeta[, c(inpConf[UI == inpdrX]$ID, inpConf[UI == inpdrY]$ID, 
                       inpConf[UI == inpsub1]$ID),  
                   with = FALSE] 
  colnames(ggData) = c("X", "Y", "sub") 
  rat = (max(ggData$X) - min(ggData$X)) / (max(ggData$Y) - min(ggData$Y)) 
  
  h5file <- H5File$new(inpH5, mode = "r") 
  h5data <- h5file[["grp"]][["data"]] 
  ggData$val = h5data$read(args = list(inpGene[inp1], quote(expr=))) 
  ggData[val < 0]$val = 0 
  h5file$close_all() 
  bgCells = FALSE 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    bgCells = TRUE 
    ggData2 = ggData[!sub %in% inpsub2] 
    ggData = ggData[sub %in% inpsub2] 
  } 
  if(inpord == "Max-1st"){ 
    ggData = ggData[order(val)] 
  } else if(inpord == "Min-1st"){ 
    ggData = ggData[order(-val)] 
  } else if(inpord == "Random"){ 
    ggData = ggData[sample(nrow(ggData))] 
  } 
   
  # Actual ggplot 
  ggOut = ggplot(ggData, aes(X, Y, color = val)) 
  if(bgCells){ 
    ggOut = ggOut + 
      geom_point(data = ggData2, color = "snow2", size = inpsiz, shape = 16) 
  } 
  ggOut = ggOut + 
    geom_point(size = inpsiz, shape = 16) + xlab(inpdrX) + ylab(inpdrY) + 
    sctheme(base_size = sList[inpfsz], XYval = inptxt) +  
    scale_color_gradientn(inp1, colours = cList[[inpcol]]) + 
      guides(color = guide_colorbar(barwidth = 15)) 
  if(inpasp == "Square") { 
    ggOut = ggOut + coord_fixed(ratio = rat) 
  } else if(inpasp == "Fixed") { 
    ggOut = ggOut + coord_fixed() 
  } 
  
    
    p <- ggOut
    
    for(i in names(sc1def$limits)){
      reduction <- substring(inpdrX, 1, nchar(inpdrX)-2)
      if(reduction == i) {
        
        p <- p+ xlim(sc1def$limits[[i]][1],sc1def$limits[[i]][2])+ ylim(sc1def$limits[[i]][3],sc1def$limits[[i]][4])
        
      }  else  {
        
        p <- p 
      }

    
  }
  
  return(p)
  
} 
 
# Plot gene coexpression on dimred 
bilinear <- function(x,y,xy,Q11,Q21,Q12,Q22){ 
  oup = (xy-x)*(xy-y)*Q11 + x*(xy-y)*Q21 + (xy-x)*y*Q12 + x*y*Q22 
  oup = oup / (xy*xy) 
  return(oup) 
} 
scDRcoex <- function(inpConf, inpMeta, inpdrX, inp, 
                     inpsub1, inpsub2, inpH5, inpGene, 
                     inpsiz, inpcol, inpord, inpfsz, inpasp, inptxt){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Prepare ggData 
  inpdrY = paste0(substring(inpdrX, 1, nchar(inpdrX)-1),2)
  ggData = inpMeta[, c(inpConf[UI == inpdrX]$ID, inpConf[UI == inpdrY]$ID, 
                       inpConf[UI == inpsub1]$ID),  
                   with = FALSE] 
  colnames(ggData) = c("X", "Y", "sub") 
  rat = (max(ggData$X) - min(ggData$X)) / (max(ggData$Y) - min(ggData$Y)) 
  
  #####NewCodes##############
  geneList = scGeneList(inp, inpGene) 
  geneList = geneList[present == TRUE] 
  # shiny::validate(need(nrow(geneList) <= 50, "More than 50 genes to plot! Please reduce the gene list!"))
  shiny::validate(need(nrow(geneList) > 1, "Please input at least 2 genes to plot!")) 
  pbmc <- AddModuleScore(combined,
                         features = list(signature = geneList$gene),
                         name="signature")
  
  reduction <- substring(inpdrX, 1, nchar(inpdrX)-2)
  
  ggOut <- FeaturePlot(pbmc, reduction = reduction,
                       features = c("signature1"), label = FALSE, repel = TRUE, pt.size = inpsiz) +
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))+
    ggtitle("")

  ggOut = ggOut +
    #   geom_point(size = inpsiz, shape = 16, color = ggData$cMix) + 
    #   xlab(inpdrX) + ylab(inpdrY) + 
    sctheme(base_size = sList[inpfsz], XYval = inptxt) +
    # scale_color_gradientn(inp1, colours = cList[[1]]) +
    guides(color = guide_colorbar(barwidth = 15))
  if(inpasp == "Square") {
    ggOut = ggOut + coord_fixed(ratio = rat)
  } else if(inpasp == "Fixed") {
    ggOut = ggOut + coord_fixed()
  }
  
  plotAllLayers <- function(ggOut){
    
    p <- ggOut
    
    for(i in names(sc1def$limits)){
      
      reduction <- substring(inpdrX, 1, nchar(inpdrX)-2)
      if(reduction == i) {
        
        
        p <- p+ xlim(sc1def$limits[[i]][1],sc1def$limits[[i]][2])+ ylim(sc1def$limits[[i]][3],sc1def$limits[[i]][4])
        
      }  else  {
        
        p <- p 
      }
      print(p)
    }
    
  }
  
  return(plotAllLayers(ggOut))
  
  
  
  
  # return(ggOut) 
} 
 
scDRcoexLeg <- function(inp1, inp2, inpcol, inpfsz){ 
  # Generate coex color palette 
  cInp = strsplit(inpcol, "; ")[[1]] 
  if(cInp[1] == "Red (Gene1)"){ 
    c10 = c(255,0,0) 
  } else if(cInp[1] == "Orange (Gene1)"){ 
    c10 = c(255,140,0) 
  } else { 
    c10 = c(0,255,0) 
  } 
  if(cInp[2] == "Green (Gene2)"){ 
    c01 = c(0,255,0) 
  } else { 
    c01 = c(0,0,255) 
  } 
  c00 = c(217,217,217) ; c11 = c10 + c01 
  nGrid = 16; nPad = 2; nTot = nGrid + nPad * 2 
  gg = data.table(v1 = rep(0:nTot,nTot+1), v2 = sort(rep(0:nTot,nTot+1))) 
  gg$vv1 = gg$v1 - nPad ; gg[vv1 < 0]$vv1 = 0; gg[vv1 > nGrid]$vv1 = nGrid 
  gg$vv2 = gg$v2 - nPad ; gg[vv2 < 0]$vv2 = 0; gg[vv2 > nGrid]$vv2 = nGrid 
  gg$cR = bilinear(gg$vv1, gg$vv2, nGrid, c00[1], c10[1], c01[1], c11[1]) 
  gg$cG = bilinear(gg$vv1, gg$vv2, nGrid, c00[2], c10[2], c01[2], c11[2]) 
  gg$cB = bilinear(gg$vv1, gg$vv2, nGrid, c00[3], c10[3], c01[3], c11[3]) 
  gg$cMix = rgb(gg$cR, gg$cG, gg$cB, maxColorValue = 255) 
  gg = gg[, c("v1", "v2", "cMix")] 
  
  # Actual ggplot 
  ggOut = ggplot(gg, aes(v1, v2)) + 
    geom_tile(fill = gg$cMix) + 
    xlab(inp1) + ylab(inp2) + coord_fixed(ratio = 1) + 
    scale_x_continuous(breaks = c(0, nTot), label = c("low", "high")) + 
    scale_y_continuous(breaks = c(0, nTot), label = c("low", "high")) + 
    sctheme(base_size = sList[inpfsz], XYval = TRUE) 
  return(ggOut) 
} 
 
scDRcoexNum <- function(inpConf, inpMeta, inp1, inp2, 
                        inpsub1, inpsub2, inpH5, inpGene){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Prepare ggData 
  ggData = inpMeta[, c(inpConf[UI == inpsub1]$ID), with = FALSE] 
  colnames(ggData) = c("sub") 
  h5file <- H5File$new(inpH5, mode = "r") 
  h5data <- h5file[["grp"]][["data"]] 
  ggData$val1 = h5data$read(args = list(inpGene[inp1], quote(expr=))) 
  ggData[val1 < 0]$val1 = 0 
  ggData$val2 = h5data$read(args = list(inpGene[inp2], quote(expr=))) 
  ggData[val2 < 0]$val2 = 0 
  h5file$close_all() 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    ggData = ggData[sub %in% inpsub2] 
  } 
  
  # Actual data.table 
  ggData$express = "none" 
  ggData[val1 > 0]$express = inp1 
  ggData[val2 > 0]$express = inp2 
  ggData[val1 > 0 & val2 > 0]$express = "both" 
  ggData$express = factor(ggData$express, levels = unique(c("both", inp1, inp2, "none"))) 
  ggData = ggData[, .(nCells = .N), by = "express"] 
  ggData$percent = 100 * ggData$nCells / sum(ggData$nCells) 
  ggData = ggData[order(express)] 
  colnames(ggData)[1] = "expression > 0" 
  return(ggData) 
} 
 
# Plot violin / boxplot 
scVioBox <- function(inpConf, inpMeta, inp1, inp2, 
                     inpsub1, inpsub2, inpsub3, inpH5, inpGene, 
                     inptyp, inppts, inpsiz, inpfsz){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Prepare ggData 
  ggData = inpMeta[, c(inpConf[UI == inp1]$ID, inpConf[UI == inpsub1]$ID), 
                   with = FALSE] 
  colnames(ggData) = c("X", "sub") 
  
  # Load in either cell meta or gene expr
  if(inp2 %in% inpConf$UI){ 
    ggData$val = inpMeta[[inpConf[UI == inp2]$ID]] 
  } else { 
    h5file <- H5File$new(inpH5, mode = "r") 
    h5data <- h5file[["grp"]][["data"]] 
    ggData$val = h5data$read(args = list(inpGene[inp2], quote(expr=))) 
    ggData[val < 0]$val = 0 
    set.seed(42) 
    tmpNoise = rnorm(length(ggData$val)) * diff(range(ggData$val)) / 1000 
    ggData$val = ggData$val + tmpNoise 
    h5file$close_all() 
  } 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    ggData = ggData[sub %in% inpsub2] 
  } 
  
  # Do factoring 
  # ggCol = strsplit(inpConf[UI == inp1]$fCL, "\\|")[[1]] 
  ggCol <- inpsub3
  names(ggCol) = levels(ggData$X) 
  ggLvl = levels(ggData$X)[levels(ggData$X) %in% unique(ggData$X)] 
  ggData$X = factor(ggData$X, levels = ggLvl) 
  ggCol = ggCol[ggLvl] 
  
  # Actual ggplot 
  if(inptyp == "violin"){ 
    ggOut = ggplot(ggData, aes(X, val, fill = X)) + geom_violin(scale = "width") 
  } else { 
    ggOut = ggplot(ggData, aes(X, val, fill = X)) + geom_boxplot() 
  } 
  if(inppts){ 
    ggOut = ggOut + geom_jitter(size = inpsiz, shape = 16) 
  } 
  ggOut = ggOut + xlab(inp1) + ylab(inp2) + 
    sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +  
    scale_fill_manual("", values = ggCol) +
    theme(legend.position = "none")
  return(ggOut) 
} 
 
# Plot proportion plot 
scProp <- function(inpConf, inpMeta, inp1, inp2, inpsub1, inpsub2,inpsub3, 
                   inptyp, inpflp, inpfsz){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Prepare ggData 
  ggData = inpMeta[, c(inpConf[UI == inp1]$ID, inpConf[UI == inp2]$ID, 
                       inpConf[UI == inpsub1]$ID),  
                   with = FALSE] 
  colnames(ggData) = c("X", "grp", "sub") 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    ggData = ggData[sub %in% inpsub2] 
  } 
  ggData = ggData[, .(nCells = .N), by = c("X", "grp")] 
  ggData = ggData[, {tot = sum(nCells) 
                      .SD[,.(pctCells = 100 * sum(nCells) / tot, 
                             nCells = nCells), by = "grp"]}, by = "X"] 
  
  # Do factoring 
  # ggCol = strsplit(inpConf[UI == inp2]$fCL, "\\|")[[1]] 
  ggCol <- inpsub3
  names(ggCol) = levels(ggData$grp) 
  ggLvl = levels(ggData$grp)[levels(ggData$grp) %in% unique(ggData$grp)] 
  ggData$grp = factor(ggData$grp, levels = ggLvl) 
  ggCol = ggCol[ggLvl] 
  
  # Actual ggplot 
  if(inptyp == "Proportion"){ 
    ggOut = ggplot(ggData, aes(X, pctCells, fill = grp)) + 
      geom_col() + ylab("Cell Proportion (%)") 
  } else { 
    ggOut = ggplot(ggData, aes(X, nCells, fill = grp)) + 
      geom_col() + ylab("Number of Cells") 
  } 
  if(inpflp){ 
    ggOut = ggOut + coord_flip() 
  } 
  ggOut = ggOut + xlab(inp1) + 
    sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +  
    scale_fill_manual("", values = ggCol) + 
    theme(legend.position = "right") 
  return(ggOut) 
} 
 
# Get gene list 
scGeneList <- function(inp, inpGene){ 
  geneList = data.table(gene = unique(trimws(strsplit(inp, ",|;|
")[[1]])), 
                        present = TRUE) 
  geneList[!gene %in% names(inpGene)]$present = FALSE 
  return(geneList) 
} 


###################### Plot gene expression bubbleplot / heatmap for genes #############################################################
# create heatmap matrix

Creat_matrix <-  function(seMarker){
  log2Norm = TRUE
  plotLog2FC = TRUE
  scaleTo = 10^4
  scaleRows = TRUE
  limits = c(-2, 2)
  
  #Evaluate AssayNames
  assayNames <- names(SummarizedExperiment::assays(seMarker))
  for(an in assayNames){
    eval(parse(text=paste0(an, " <- ", "SummarizedExperiment::assays(seMarker)[['", an, "']]")))
  }
  # cutOff = "FDR <= 0.01 & Log2FC >= 0.5"
  cutOff = "FDR <= 1"
  passMat <- eval(parse(text=cutOff))
  for(an in assayNames){
    eval(parse(text=paste0("rm(",an,")")))
  }
  #Now Get Values
  if(ncol(seMarker) <= 2){
  if(!plotLog2FC){
    stop("Must use plotLog2FC = TRUE when ncol(seMarker) <= 2!")
  }
  }
  #Get Matrix
  if(plotLog2FC){
    mat <- as.matrix(SummarizedExperiment::assays(seMarker)[["Log2FC"]])
  }else{
    mat <- as.matrix(SummarizedExperiment::assays(seMarker)[["Mean"]])
    if(log2Norm){
      mat <- log2(t(t(mat)/colSums(mat)) * scaleTo + 1)
    }
    if(scaleRows){
      mat <- sweep(mat - rowMeans(mat), 1, matrixStats::rowSds(mat), `/`)
    }
  }
  mat[mat > max(limits)] <- max(limits)
  mat[mat < min(limits)] <- min(limits)    
  
  if(ncol(mat) == 1){
    idx <- which(rowSums(passMat, na.rm = TRUE) > 0)
  }else{
    idx <- which(rowSums(passMat, na.rm = TRUE) > 0 & matrixStats::rowVars(mat) != 0 & !is.na(matrixStats::rowVars(mat)))
  }
  subsetMarkers = NULL
  if(!is.null(subsetMarkers)) {
    if(length(which(subsetMarkers %ni% 1:nrow(mat))) == 0){
      idx <- subsetMarkers
    } else {
      stop("Rownames / indices provided to the subsetMarker parameter are outside of the boundaries of seMarker.")
    }
    
  }    
  mat <- mat[idx,,drop=FALSE]
  passMat <- passMat[idx,,drop=FALSE]
  
  if(nrow(mat) == 0){
    stop("No Makers Found!")
  }   
  #add rownames
  rd <- SummarizedExperiment::rowData(seMarker)[idx,]
  if(is.null(rd$name)){
    rn <- paste0(rd$seqnames,":",rd$start,"-",rd$end)
  }else{
    if(sum(duplicated(rd$name)) > 0){
      rn <- paste0(rd$seqnames,":",rd$name)
    }else{
      rn <- rd$name
    }
  }
  rownames(mat) <- rn
  rownames(passMat) <- rn
  return(mat)
  
}


scBubbHeat <- function(inpConf, inpMeta, inp, inpGrp, inpPlt, 
                        inpsub1, inpsub2, inpH5, inpGene, inpScl, inpRow, inpCol, 
                        inpcols, inpfsz, save = FALSE){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Identify genes that are in our dataset 
  geneList = scGeneList(inp, inpGene) 
  geneList = geneList[present == TRUE] 
  # shiny::validate(need(nrow(geneList) <= 50, "More than 50 genes to plot! Please reduce the gene list!"))
  shiny::validate(need(nrow(geneList) > 1, "Please input at least 2 genes to plot!")) 
  
  
  # Prepare ggData 
  # ======================================
  treatment <- names(getCellColData(proj))[grep('condition_',names(getCellColData(proj)))]
  
  if(inpGrp=="Clusters"){
    seMarker <- seMarker_cluster


    for(iGene in geneList$gene){
    seMarker <- seMarker[which(rowData(seMarker)$name%in%geneList$gene),]
    }
    mat = Creat_matrix(seMarker)
    ggData1 = as.data.frame(mat)
    # add gene names into the first col
    ggData1 = data.frame(X = rownames(ggData1), ggData1)

    nClust <- ncol(mat)
    # 
    req_meta_data <- read.csv("./tables/req_meta_data.csv")[,c(
      'X'
      , 'Clusters'
      ,treatment
      ,"Sample")]

    d <-  list()
    out <-  list()
    n1 <-  list()
    n2 = nrow(ggData1)
    cluster <-  list()

    for (i in seq_along(1:nClust)){

      d[[i]] <- ggData1[,c(1,i+1)]
      cluster[[i]] <- paste0("C",i)
      n1[[i]] = nrow(req_meta_data[which(req_meta_data$Clusters==cluster[[i]]),])

    }
    # 

    out <- mapply(function(x,y) DataFrame(geneName=rep(x[,1],y),val=rep(x[,2],y)),d,n1)


    out <- mapply(function(x,y) DataFrame(x,sampleID=rep(req_meta_data[which(req_meta_data$Clusters==y),]$X,1,each=n2)), out,cluster)


    out <- lapply(out, function(x) DataFrame(x,grpBy=rep("Clusters",nrow(x))))


    out <- mapply(function(x,y) DataFrame(x,sub=rep(y,nrow(x))),out,cluster)
    # 
    # 
    h5data <- as.data.frame(do.call("rbind", out))
    #  
    
    
  
  } else if (inpGrp == "Sample" || inpGrp == "SampleName") {
    
      seMarker <- seMarker_sample 
      
      for(iGene in geneList$gene){
        seMarker <- seMarker[which(rowData(seMarker)$name%in%geneList$gene),]
      }
      mat = Creat_matrix(seMarker)
      
      ggData1 = as.data.frame(mat)
      
      # add gene names into the first col
      ggData1 = data.frame(X = rownames(ggData1), ggData1)
      nSample <- ncol(mat) 
      
      req_meta_data <- read.csv("./tables/req_meta_data.csv")[,c(
        'X'
        , 'Clusters'
        ,treatment
        ,"Sample")]
      
      d <-  list()
      out <-  list()
      n1 <-  list()
      n2 = nrow(ggData1)
      Allsamples <- colnames(mat)

      samples <- list()
      
      for (i in seq_along(1:nSample)){

        d[[i]] <- ggData1[,c(1,i+1)]
        samples[[i]] <- Allsamples[i]
        n1[[i]] = nrow(req_meta_data[which(req_meta_data$Sample==samples[[i]]),])

      }

      out <- mapply(function(x,y) DataFrame(geneName=rep(x[,1],y),val=rep(x[,2],y)),d,n1)


      out <- mapply(function(x,y) DataFrame(x,sampleID=rep(req_meta_data[which(req_meta_data$Sample==y),]$X,1,each=n2)), out,samples)


      out <- lapply(out, function(x) DataFrame(x,grpBy=rep("Sample",nrow(x))))


      out <- mapply(function(x,y) DataFrame(x,sub=rep(y,nrow(x))),out,samples)

      
      h5data <- as.data.frame(do.call("rbind", out))
    
    
    
  } else {
  
    
    idx = as.integer(unlist(strsplit(inpGrp,"_"))[2])
    
    
    seMarker <- seMarker_treatment[[idx]]

    for(iGene in geneList$gene){
      seMarker <- seMarker[which(rowData(seMarker)$name%in%geneList$gene),]
    }
    mat = Creat_matrix(seMarker)
    ggData1 = as.data.frame(mat)
    # add gene names into the first col
    ggData1 = data.frame(X = rownames(ggData1), ggData1)

    nTreatment <- ncol(mat)

    req_meta_data <- read.csv("./tables/req_meta_data.csv")[,c(
      'X'
      , 'Clusters'
      ,treatment
      ,"Sample")]

    d <-  list()
    out <-  list()
    n1 <-  list()
    n2 = nrow(ggData1)
    AllTreatment <-  colnames(mat)
    
    Treatment <- list()
    for (i in seq_along(1:nTreatment)){

      d[[i]] <- ggData1[,c(1,i+1)]
      Treatment[[i]] <- AllTreatment[i]
      n1[[i]] = nrow(req_meta_data[which(req_meta_data[[inpGrp]]==Treatment[[i]]),])

    }


    out <- mapply(function(x,y) DataFrame(geneName=rep(x[,1],y),val=rep(x[,2],y)),d,n1)


    out <- mapply(function(x,y) DataFrame(x,sampleID=rep(req_meta_data[which(req_meta_data[[inpGrp]]==y),]$X,1,each=n2)), out,Treatment)


    out <- lapply(out, function(x) DataFrame(x,grpBy=rep(inpGrp,nrow(x))))


    out <- mapply(function(x,y) DataFrame(x,sub=rep(y,nrow(x))),out,Treatment)



    h5data <- as.data.frame(do.call("rbind", out))

     
  }

  ggData = data.table()
  for(iGene in geneList$gene){
    tmp = inpMeta[, c("sampleID", inpConf[UI == inpsub1]$ID), with = FALSE]
    colnames(tmp) = c("sampleID", "sub")
    tmp$grpBy = inpMeta[[inpConf[UI == inpGrp]$ID]]
    tmp$geneName = iGene
    h5data_tmp =  h5data[which(h5data$geneName==iGene),] 
    tmp$val = h5data_tmp[match(tmp$sampleID,h5data_tmp$sampleID),]$val
    ggData = rbindlist(list(ggData, tmp))
  }
   
 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){
    ggData = ggData[sub %in% inpsub2]
     
  }
  shiny::validate(need(uniqueN(ggData$grpBy) > 1, "Only 1 group present, unable to plot!"))
  colRange = c(-max(abs(range(ggData$val))), max(abs(range(ggData$val))))

  if(inpRow){
    clusterRows = TRUE
  } else {
    clusterRows = FALSE
  }
  if(inpCol){
    clusterCols = TRUE
  } else {
    clusterCols = FALSE
  }
  
  cmplxData = data.table()
  for(iGene in geneList$gene){
    for(iclst in as.character(unique(ggData$grpBy))){
      temp <-  ggData[which(ggData$geneName==iGene & ggData$grpBy==iclst),][1,]
      cmplxData = rbindlist(list(cmplxData, temp))
    }
  }
  cmplxData2 = data.table()
  for(iGene in geneList$gene){
    
    temp = as.data.frame(t(setDT(cmplxData)[order(-(geneName==iGene)), .SD[1L] ,.(grpBy)][,c("geneName","val")]))[-c(1),]
    colnames(temp)= as.character(cmplxData[which(cmplxData$geneName==iGene),]$grpBy)
    cmplxData2 = rbindlist(list(cmplxData2, temp))
  }
  cmplxData2 = cmplxData2[ , lapply(.SD, as.numeric)]
  
  cmplxData3 = as.data.frame(cmplxData2)
  cmplxData3 = cmplxData3[,order(as.numeric(gsub("C","",colnames(cmplxData3))))]
  rownames(cmplxData3) = geneList$gene
  
  
  
  # Actual plot according to plottype 
  if(inpPlt == "Bubbleplot"){
    
 #we don't use bubble plot
    
  } else {
    # Heatmap 
    customRowLabel <- rownames(cmplxData3)
    customRowLabelIDs <- as.numeric(rownames(as.data.frame(rownames(cmplxData3))))
    fontSizeLabels <- 7
    ggOut = draw(Heatmap(

      #Main Stuff
      matrix = as.matrix(cmplxData3),
      #     name = name,
      col = cList[[inpcols]],
      show_column_names = T,
      cluster_columns = clusterCols,
      show_column_dend = T,
      show_row_names = T,
      cluster_rows = clusterRows)+
    rowAnnotation(foo = anno_mark(at= customRowLabelIDs, labels = customRowLabel
                                     ,labels_gp = gpar(fontsize = fontSizeLabels)
    )))
    
  }
  

  return(ggOut)
}
 
# Plot gene expression bubbleplot / heatmap for motifs ########################################################################################

# new functions required here!
.rowScale <- function(mat = NULL, min = NULL, max = NULL){
  if(!is.null(min)){
    rMin <- min
  }else{
    rMin <- matrixStats::rowMins(mat)
  }
  if(!is.null(max)){
    rMax <- max
  }else{
    rMax <- matrixStats::rowMaxs(mat)
  }
  rScale <- rMax - rMin
  matDiff <- mat - rMin
  matScale <- matDiff/rScale
  out <- list(mat=matScale, min=rMin, max=rMax)
  return(out)
}

.binarySort <- function(m = NULL, scale = FALSE, cutOff = 1, lmat = NULL, clusterCols = TRUE){
  
  if(is.null(lmat)){
    #Compute Row-Zscores
    if(scale){
      lmat <- sweep(m - rowMeans(m), 1, matrixStats::rowSds(m), `/`)
    }else{
      lmat <- m
    }
    lmat <- lmat >= cutOff
  }
  
  #Transpose
  m <- t(m)
  lmat <- t(lmat)
  
  #Identify Column Ordering
  if(clusterCols){
    hc <- hclust(dist(m))
    colIdx <- hc$order
    m <- t(m[colIdx,])
    lmat <- t(lmat[colIdx,])
  }else{
    m <- t(m)
    lmat <- t(lmat)
    hc <- NULL
  }
  
  #Identify Row Ordering
  rowIdx <- do.call("order", c(as.data.frame(lmat)[seq_len(ncol(lmat))], list(decreasing = TRUE)))
  m <- t(m[rowIdx,])
  lmat <- t(lmat[rowIdx,])
  
  #Transpose
  m <- t(m)
  lmat <- t(lmat)
  
  return(list(mat = m, hclust = hc))
  
}
Creat_matrix_motif <- function(
    seEnrich = NULL,
    pal = paletteContinuous(set = "comet", n = 100),
    n = 100,
    cutOff = 0,
    pMax = Inf,
    clusterCols = TRUE,
    binaryClusterRows = TRUE,
    labelRows = TRUE,
    rastr = TRUE,
    transpose = FALSE,
    returnMatrix = TRUE
){
  
  
  mat <- assays(seEnrich)[["mlog10Padj"]]
  
  keep <- lapply(seq_len(ncol(mat)), function(x){
    idx <- head(order(mat[, x], decreasing = TRUE), n)
    rownames(mat)[idx[which(mat[idx,x] >= cutOff)]]
  }) %>% unlist %>% unique
  mat <- mat[keep, ,drop = FALSE]
  
  if(nrow(mat)==0){
    stop("No enrichments found for your cutoff!")
  }
  
  passMat <- lapply(seq_len(nrow(mat)), function(x){
    (mat[x, ] >= 0.9*max(mat[x, ])) * 1
  }) %>% Reduce("rbind", .) %>% data.frame
  colnames(passMat) <- colnames(mat)
  
  mat[mat > pMax] <- pMax
  
  if(nrow(mat)==0){
    stop("No enrichments found for your cutofff!!!")
  }
  
  mat <- .rowScale(as.matrix(mat), min = 0)
  if(pMax != 100){
    rownames(mat[[1]]) <- paste0(rownames(mat[[1]]), " (",round(mat$max),")")
    rownames(passMat) <- rownames(mat[[1]])
  }
  
  mat2 <- mat[[1]] * 100
  
  
  motif_lst <- unique(rownames(mat2))
  split_string <- strsplit(motif_lst, split = "\\(")
  fun1 <- function(list, nth){
    sapply(list, `[` , 1)
  }
  req_motifs1 <- gsub("_","-",fun1(split_string))
  req_motifs1 <- gsub(" ","",req_motifs1)
  
  rownames(mat2) <- req_motifs1
  
  
  if(nrow(mat2) > 1 & ncol(mat2) > 1){
    if(binaryClusterRows){
      #       cn <- order(colMeans(mat2), decreasing=TRUE)
      bS <- .binarySort(mat2, lmat = passMat[rownames(mat2), colnames(mat2)], clusterCols = clusterCols)
      mat2 <- bS[[1]][,colnames(mat2)]
      clusterRows <- FALSE
      clusterCols <- bS[[2]]
    }else{
      clusterRows <- TRUE
    }
  }else{
    clusterCols <- NULL
    clusterRows <- FALSE
  }
  
  if(nrow(mat2) > 100){
    borderColor <- FALSE
  }else{
    borderColor <- TRUE
  }
  
  #   if(transpose){
  
  #     if(!is.null(clusterCols)){
  #       mat2 <- t(mat2[,clusterCols$order,drop=FALSE])
  #     }else{
  #       mat2 <- t(mat2)
  #     }
  
  if(returnMatrix){
    return(mat2)
    #     }
  }
}
# =========================



scBubbHeat2 <- function(inpConf, inpMeta, inp, inpGrp, inpPlt, 
                       inpsub1, inpsub2, inpH5, inpGene, inpScl, inpRow, inpCol, 
                       inpcols, inpfsz, save = FALSE){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Identify motifs that are in our dataset 
  geneList = scGeneList(inp, inpGene) 
  geneList = geneList[present == TRUE] 
  #shiny::validate(need(nrow(geneList) <= 50, "More than 50 genes to plot! Please reduce the gene list!"))
  shiny::validate(need(nrow(geneList) > 1, "Please input at least 2 genes to plot!"))
  
  
  # Prepare ggData 
  # ======================================
  treatment <- names(getCellColData(proj))[grep('condition_',names(getCellColData(proj)))]
  
  if(inpGrp=="Clusters"){
    seEnrich <- seEnrich_cluster
    for(iGene in geneList$gene){
      seEnrich <-seEnrich[which(rownames(seEnrich)%in%geneList$gene),]
    }
    mat = Creat_matrix_motif(seEnrich)
    ggData1 = as.data.frame(mat)
    # add gene names into the first col
    ggData1 = data.frame(X = rownames(ggData1), ggData1)
    
    nClust <- ncol(mat)
    # 
    req_meta_data <- read.csv("./tables/req_meta_data.csv")[,c(
      'X'
      , 'Clusters'
      ,treatment
      ,"Sample")]
    
    d <-  list()
    out <-  list()
    n1 <-  list()
    n2 = nrow(ggData1)
    cluster <-  list()
    
    for (i in seq_along(1:nClust)){
      
      d[[i]] <- ggData1[,c(1,i+1)]
      cluster[[i]] <- paste0("C",i)
      n1[[i]] = nrow(req_meta_data[which(req_meta_data$Clusters==cluster[[i]]),])
    }
    out <- mapply(function(x,y) DataFrame(geneName=rep(x[,1],y),val=rep(x[,2],y)),d,n1)

    out <- mapply(function(x,y) DataFrame(x,sampleID=rep(req_meta_data[which(req_meta_data$Clusters==y),]$X,1,each=n2)), out,cluster)
    
    out <- lapply(out, function(x) DataFrame(x,grpBy=rep("Clusters",nrow(x))))
    
    out <- mapply(function(x,y) DataFrame(x,sub=rep(y,nrow(x))),out,cluster)

    h5data <- as.data.frame(do.call("rbind", out))

  } else if (inpGrp == "Sample" || inpGrp == "SampleName") {

    seEnrich <- seEnrich_sample
    
    for(iGene in geneList$gene){
      seEnrich <- seEnrich[which(rownames(seEnrich)%in%geneList$gene), ]
    }
    mat = Creat_matrix_motif(seEnrich)
    
    ggData1 = as.data.frame(mat)
    
    # add gene names into the first col
    ggData1 = data.frame(X = rownames(ggData1), ggData1)
    nSample <- ncol(mat) 
    
    req_meta_data <- read.csv("./tables/req_meta_data.csv")[, c(
      'X'
      , 'Clusters'
      ,treatment
      ,"Sample")]
    
    d <-  list()
    out <-  list()
    n1 <-  list()
    n2 = nrow(ggData1)
    Allsamples <- colnames(mat)
    
    samples <- list()
    
    for (i in seq_along(1:nSample)){
      
      d[[i]] <- ggData1[,c(1,i+1)]
      samples[[i]] <- Allsamples[i]
      n1[[i]] = nrow(req_meta_data[which(req_meta_data$Sample==samples[[i]]),])
      
    }
    
    out <- mapply(function(x,y) DataFrame(geneName=rep(x[,1],y),val=rep(x[,2],y)),d,n1)
    
    out <- mapply(function(x,y) DataFrame(x,sampleID=rep(req_meta_data[which(req_meta_data$Sample==y),]$X,1,each=n2)), out,samples)
    
    out <- lapply(out, function(x) DataFrame(x,grpBy=rep("Sample",nrow(x))))
    
    out <- mapply(function(x,y) DataFrame(x,sub=rep(y,nrow(x))),out,samples)
    
    h5data <- as.data.frame(do.call("rbind", out))
    
  } else {
    idx = as.integer(unlist(strsplit(inpGrp,"_"))[2])
    
    seEnrich <- seEnrich_treatment[[idx]]
    
    for(iGene in geneList$gene){
      seEnrich <- seEnrich[which(rownames(seEnrich)%in%geneList$gene),]
    }
    mat = Creat_matrix_motif(seEnrich)
    ggData1 = as.data.frame(mat)
    # add gene names into the first col
    ggData1 = data.frame(X = rownames(ggData1), ggData1)
    
    nTreatment <- ncol(mat)
    
    req_meta_data <- read.csv("./tables/req_meta_data.csv")[,c(
      'X'
      , 'Clusters'
      ,treatment
      ,"Sample")]
    
    d <-  list()
    out <-  list()
    n1 <-  list()
    n2 = nrow(ggData1)
    AllTreatment <-  colnames(mat)
    
    Treatment <- list()
    for (i in seq_along(1:nTreatment)){
      
      d[[i]] <- ggData1[,c(1,i+1)]
      Treatment[[i]] <- AllTreatment[i]
      n1[[i]] = nrow(req_meta_data[which(req_meta_data[[inpGrp]]==Treatment[[i]]),])
      
    }
    #
    #
    #
    out <- mapply(function(x,y) DataFrame(geneName=rep(x[,1],y),val=rep(x[,2],y)),d,n1)

    out <- mapply(function(x,y) DataFrame(x,sampleID=rep(req_meta_data[which(req_meta_data[[inpGrp]]==y),]$X,1,each=n2)), out,Treatment)

    out <- lapply(out, function(x) DataFrame(x,grpBy=rep(inpGrp,nrow(x))))
    
    out <- mapply(function(x,y) DataFrame(x,sub=rep(y,nrow(x))),out,Treatment)
    #
    #
    #
    h5data <- as.data.frame(do.call("rbind", out))
    #
    
  }
  
  ggData = data.table()
  for(iGene in geneList$gene){
    tmp = inpMeta[, c("sampleID", inpConf[UI == inpsub1]$ID), with = FALSE]
    colnames(tmp) = c("sampleID", "sub")
    tmp$grpBy = inpMeta[[inpConf[UI == inpGrp]$ID]]
    tmp$geneName = iGene
    h5data_tmp =  h5data[which(h5data$geneName==iGene),] 
    tmp$val = h5data_tmp[match(tmp$sampleID,h5data_tmp$sampleID),]$val
    ggData = rbindlist(list(ggData, tmp))
  }
  # 
  
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){
    ggData = ggData[sub %in% inpsub2]
    #   
  }
  shiny::validate(need(uniqueN(ggData$grpBy) > 1, "Only 1 group present, unable to plot!"))
  colRange = c(-max(abs(range(ggData$val))), max(abs(range(ggData$val))))
  
  if(inpRow){
    clusterRows = TRUE
  } else {
    clusterRows = FALSE
  }
  if(inpCol){
    clusterCols = TRUE
  } else {
    clusterCols = FALSE
  }
  
  cmplxData = data.table()
  for(iGene in geneList$gene){
    for(iclst in as.character(unique(ggData$grpBy))){
      temp <-  ggData[which(ggData$geneName==iGene & ggData$grpBy==iclst),][1,]
      cmplxData = rbindlist(list(cmplxData, temp))
    }
  }
  cmplxData2 = data.table()
  for(iGene in geneList$gene){
    
    temp = as.data.frame(t(setDT(cmplxData)[order(-(geneName==iGene)), .SD[1L] ,.(grpBy)][,c("geneName","val")]))[-c(1),]
    colnames(temp)= as.character(cmplxData[which(cmplxData$geneName==iGene),]$grpBy)
    cmplxData2 = rbindlist(list(cmplxData2, temp))
  }
  cmplxData2 = cmplxData2[ , lapply(.SD, as.numeric)]
  
  cmplxData3 = as.data.frame(cmplxData2)
  cmplxData3 = cmplxData3[,order(as.numeric(gsub("C","",colnames(cmplxData3))))]
  rownames(cmplxData3) = geneList$gene
  
  
  
  # Actual plot according to plottype 
  if(inpPlt == "Bubbleplot"){
    
    #we don't use bubble plot
    
  } else {
    # Heatmap 
    customRowLabel <- rownames(cmplxData3)
    customRowLabelIDs <- as.numeric(rownames(as.data.frame(rownames(cmplxData3))))
    fontSizeLabels <- 7
    ggOut = draw(Heatmap(
      
      #Main Stuff
      matrix = as.matrix(cmplxData3),
      #     name = name,
      col = cList[[inpcols]],
      show_column_names = T,
      cluster_columns = clusterCols,
      show_column_dend = T,
      show_row_names = T,
      cluster_rows = clusterRows)+
        rowAnnotation(foo = anno_mark(at= customRowLabelIDs, labels = customRowLabel
                                      ,labels_gp = gpar(fontsize = fontSizeLabels)
        )))
    
  } 
    
  
  
  

  return(ggOut) ########################################################################################
}
 
 
 
 
 
### Start server code 
shinyServer(function(input, output, session) { 
  ### For all tags and Server-side selectize 
  observe_helpers() 
 optCrt="{ option_create: function(data,escape) {return('<div class=\"create\"><strong>' + '</strong></div>');} }" 
  updateSelectizeInput(session, "sc1a1inp2", choices = names(sc1gene), server = TRUE, 
                       selected = sc1def$gene1, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "sc1a3inp1", choices = names(sc1gene), server = TRUE, 
                       selected = sc1def$gene2, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "sc1a3inp2", choices = names(sc1gene), server = TRUE, 
                       selected = sc1def$gene2, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "sc1b2inp1", choices = names(sc1gene), server = TRUE, 
                       selected = sc1def$gene1, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "sc1b2inp2", choices = names(sc1gene), server = TRUE, 
                       selected = sc1def$gene2, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "sc1c1inp2", server = TRUE, 
                       choices = c(sc1conf[is.na(fID)]$UI,names(sc1gene)), 
                       selected = sc1conf[is.na(fID)]$UI[1], options = list( 
                         maxOptions = length(sc1conf[is.na(fID)]$UI) + 3, 
                         create = TRUE, persist = TRUE, render = I(optCrt))) 
  
  ##### START OF NEW CODE TO ADD ####################
  ### Tab1.x1: New tab for volcano

   
  ##### END OF NEW CODE TO ADD #########################  
 
  ### Plots for tab a1 
  output$sc1a1sub1.ui <- renderUI({ 
    sub = strsplit(sc1conf[UI == input$sc1a1sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("sc1a1sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$sc1a1sub1non, { 
    sub = strsplit(sc1conf[UI == input$sc1a1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1a1sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$sc1a1sub1all, { 
    sub = strsplit(sc1conf[UI == input$sc1a1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1a1sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$sc1a1oup1 <- renderPlot({ 
    scDRcell(sc1conf, sc1meta, input$sc1a1drX, input$sc1a1inp1,  
             input$sc1a1sub1, input$sc1a1sub2, 
             input$sc1a1siz, input$sc1a1col1, input$sc1a1ord1, 
             input$sc1a1fsz, input$sc1a1asp, input$sc1a1txt, input$sc1a1lab1) 
  }) 
  output$sc1a1oup1.ui <- renderUI({ 
    plotOutput("sc1a1oup1", height = pList[input$sc1a1psz]) 
  }) 
  output$sc1a1oup1.pdf <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1a1drX,"_",input$sc1a1drY,"_",  
                                   input$sc1a1inp1,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$sc1a1oup1.h, width = input$sc1a1oup1.w, useDingbats = FALSE, 
      plot = scDRcell(sc1conf, sc1meta, input$sc1a1drX, input$sc1a1inp1,   
                      input$sc1a1sub1, input$sc1a1sub2, 
                      input$sc1a1siz, input$sc1a1col1, input$sc1a1ord1,  
                      input$sc1a1fsz, input$sc1a1asp, input$sc1a1txt, input$sc1a1lab1) ) 
  }) 
  output$sc1a1oup1.png <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1a1drX,"_",input$sc1a1drY,"_",  
                                   input$sc1a1inp1,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$sc1a1oup1.h, width = input$sc1a1oup1.w, 
      plot = scDRcell(sc1conf, sc1meta, input$sc1a1drX, input$sc1a1inp1,   
                      input$sc1a1sub1, input$sc1a1sub2, 
                      input$sc1a1siz, input$sc1a1col1, input$sc1a1ord1,  
                      input$sc1a1fsz, input$sc1a1asp, input$sc1a1txt, input$sc1a1lab1) ) 
  }) 
  output$sc1a1.dt <- renderDataTable({ 
    ggData = scDRnum(sc1conf, sc1meta, input$sc1a1inp1, input$sc1a1inp2, 
                     input$sc1a1sub1, input$sc1a1sub2, 
                     "sc1gexpr.h5", sc1gene, input$sc1a1splt) 
    datatable(ggData, rownames = FALSE, extensions = "Buttons", 
              options = list(pageLength = -1, dom = "tB", buttons = c("copy", "csv", "excel"))) %>% 
      formatRound(columns = c("pctExpress"), digits = 2) 
  }) 
   
  output$sc1a1oup2 <- renderPlot({ 
    scDRgene(sc1conf, sc1meta, input$sc1a1drX2, input$sc1a1inp2,  
             input$sc1a1sub1, input$sc1a1sub2, 
             "sc1gexpr.h5", sc1gene, 
             input$sc1a1siz, input$sc1a1col2, input$sc1a1ord2, 
             input$sc1a1fsz, input$sc1a1asp, input$sc1a1txt) 
  }) 
  output$sc1a1oup2.ui <- renderUI({ 
    plotOutput("sc1a1oup2", height = pList[input$sc1a1psz]) 
  }) 
  output$sc1a1oup2.pdf <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1a1drX2,"_",input$sc1a1drY2,"_",  
                                   input$sc1a1inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$sc1a1oup2.h, width = input$sc1a1oup2.w, useDingbats = FALSE, 
      plot = scDRgene(sc1conf, sc1meta, input$sc1a1drX2, input$sc1a1inp2,  
                      input$sc1a1sub1, input$sc1a1sub2, 
                      "sc1gexpr.h5", sc1gene, 
                      input$sc1a1siz, input$sc1a1col2, input$sc1a1ord2, 
                      input$sc1a1fsz, input$sc1a1asp, input$sc1a1txt) ) 
  }) 
  output$sc1a1oup2.png <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1a1drX2,"_",input$sc1a1drY2,"_",  
                                   input$sc1a1inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$sc1a1oup2.h, width = input$sc1a1oup2.w, 
      plot = scDRgene(sc1conf, sc1meta, input$sc1a1drX2, input$sc1a1inp2,  
                      input$sc1a1sub1, input$sc1a1sub2, 
                      "sc1gexpr.h5", sc1gene, 
                      input$sc1a1siz, input$sc1a1col2, input$sc1a1ord2, 
                      input$sc1a1fsz, input$sc1a1asp, input$sc1a1txt) ) 
  }) 
   
   
  ### Plots for tab a2 
  output$sc1a2sub1.ui <- renderUI({ 
    sub = strsplit(sc1conf[UI == input$sc1a2sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("sc1a2sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$sc1a2sub1non, { 
    sub = strsplit(sc1conf[UI == input$sc1a2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1a2sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$sc1a2sub1all, { 
    sub = strsplit(sc1conf[UI == input$sc1a2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1a2sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$sc1a2oup1 <- renderPlot({ 
    scDRcell(sc1conf, sc1meta, input$sc1a2drX, input$sc1a2inp1,  
             input$sc1a2sub1, input$sc1a2sub2, 
             input$sc1a2siz, input$sc1a2col1, input$sc1a2ord1, 
             input$sc1a2fsz, input$sc1a2asp, input$sc1a2txt, input$sc1a2lab1) 
  }) 
  output$sc1a2oup1.ui <- renderUI({ 
    plotOutput("sc1a2oup1", height = pList[input$sc1a2psz]) 
  }) 
  output$sc1a2oup1.pdf <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1a2drX,"_",input$sc1a2drY,"_",  
                                   input$sc1a2inp1,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$sc1a2oup1.h, width = input$sc1a2oup1.w, useDingbats = FALSE, 
      plot = scDRcell(sc1conf, sc1meta, input$sc1a2drX, input$sc1a2inp1,   
                      input$sc1a2sub1, input$sc1a2sub2, 
                      input$sc1a2siz, input$sc1a2col1, input$sc1a2ord1,  
                      input$sc1a2fsz, input$sc1a2asp, input$sc1a2txt, input$sc1a2lab1) ) 
  }) 
  output$sc1a2oup1.png <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1a2drX,"_",input$sc1a2drY,"_",  
                                   input$sc1a2inp1,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$sc1a2oup1.h, width = input$sc1a2oup1.w, 
      plot = scDRcell(sc1conf, sc1meta, input$sc1a2drX, input$sc1a2inp1,   
                      input$sc1a2sub1, input$sc1a2sub2, 
                      input$sc1a2siz, input$sc1a2col1, input$sc1a2ord1,  
                      input$sc1a2fsz, input$sc1a2asp, input$sc1a2txt, input$sc1a2lab1) ) 
  }) 
   
  output$sc1a2oup2 <- renderPlot({ 
    scDRcell(sc1conf, sc1meta, input$sc1a2drX2, input$sc1a2inp2,  
             input$sc1a2sub1, input$sc1a2sub2, 
             input$sc1a2siz, input$sc1a2col2, input$sc1a2ord2, 
             input$sc1a2fsz, input$sc1a2asp, input$sc1a2txt, input$sc1a2lab2) 
  }) 
  output$sc1a2oup2.ui <- renderUI({ 
    plotOutput("sc1a2oup2", height = pList[input$sc1a2psz]) 
  }) 
  output$sc1a2oup2.pdf <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1a2drX2,"_",input$sc1a2drY2,"_",  
                                   input$sc1a2inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$sc1a2oup2.h, width = input$sc1a2oup2.w, useDingbats = FALSE, 
      plot = scDRcell(sc1conf, sc1meta, input$sc1a2drX2, input$sc1a2inp2,   
                      input$sc1a2sub1, input$sc1a2sub2, 
                      input$sc1a2siz, input$sc1a2col2, input$sc1a2ord2,  
                      input$sc1a2fsz, input$sc1a2asp, input$sc1a2txt, input$sc1a2lab2) ) 
  }) 
  output$sc1a2oup2.png <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1a2drX2,"_",input$sc1a2drY2,"_",  
                                   input$sc1a2inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$sc1a2oup2.h, width = input$sc1a2oup2.w, 
      plot = scDRcell(sc1conf, sc1meta, input$sc1a2drX2, input$sc1a2inp2,   
                      input$sc1a2sub1, input$sc1a2sub2, 
                      input$sc1a2siz, input$sc1a2col2, input$sc1a2ord2,  
                      input$sc1a2fsz, input$sc1a2asp, input$sc1a2txt, input$sc1a2lab2) ) 
  }) 
   
   
  ### Plots for tab a3 
  output$sc1a3sub1.ui <- renderUI({ 
    sub = strsplit(sc1conf[UI == input$sc1a3sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("sc1a3sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$sc1a3sub1non, { 
    sub = strsplit(sc1conf[UI == input$sc1a3sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1a3sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$sc1a3sub1all, { 
    sub = strsplit(sc1conf[UI == input$sc1a3sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1a3sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$sc1a3oup1 <- renderPlot({ 
    scDRgene(sc1conf, sc1meta, input$sc1a3drX, input$sc1a3inp1,  
             input$sc1a3sub1, input$sc1a3sub2, 
             "sc1gexpr.h5", sc1gene, 
             input$sc1a3siz, input$sc1a3col1, input$sc1a3ord1, 
             input$sc1a3fsz, input$sc1a3asp, input$sc1a3txt) 
  }) 
  output$sc1a3oup1.ui <- renderUI({ 
    plotOutput("sc1a3oup1", height = pList[input$sc1a3psz]) 
  }) 
  output$sc1a3oup1.pdf <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1a3drX,"_",input$sc1a3drY,"_",  
                                   input$sc1a3inp1,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$sc1a3oup1.h, width = input$sc1a3oup1.w, useDingbats = FALSE, 
      plot = scDRgene(sc1conf, sc1meta, input$sc1a3drX, input$sc1a3inp1,  
                      input$sc1a3sub1, input$sc1a3sub2, 
                      "sc1gexpr.h5", sc1gene, 
                      input$sc1a3siz, input$sc1a3col1, input$sc1a3ord1, 
                      input$sc1a3fsz, input$sc1a3asp, input$sc1a3txt) ) 
  }) 
  output$sc1a3oup1.png <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1a3drX,"_",input$sc1a3drY,"_",  
                                   input$sc1a3inp1,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$sc1a3oup1.h, width = input$sc1a3oup1.w, 
      plot = scDRgene(sc1conf, sc1meta, input$sc1a3drX, input$sc1a3inp1,  
                      input$sc1a3sub1, input$sc1a3sub2, 
                      "sc1gexpr.h5", sc1gene, 
                      input$sc1a3siz, input$sc1a3col1, input$sc1a3ord1, 
                      input$sc1a3fsz, input$sc1a3asp, input$sc1a3txt) ) 
  }) 
   
  output$sc1a3oup2 <- renderPlot({ 
    scDRgene(sc1conf, sc1meta, input$sc1a3drX2, input$sc1a3inp2,  
             input$sc1a3sub1, input$sc1a3sub2, 
             "sc1gexpr.h5", sc1gene, 
             input$sc1a3siz, input$sc1a3col2, input$sc1a3ord2, 
             input$sc1a3fsz, input$sc1a3asp, input$sc1a3txt) 
  }) 
  output$sc1a3oup2.ui <- renderUI({ 
    plotOutput("sc1a3oup2", height = pList[input$sc1a3psz]) 
  }) 
  output$sc1a3oup2.pdf <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1a3drX2,"_",input$sc1a3drY2,"_",  
                                   input$sc1a3inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$sc1a3oup2.h, width = input$sc1a3oup2.w, useDingbats = FALSE, 
      plot = scDRgene(sc1conf, sc1meta, input$sc1a3drX2, input$sc1a3inp2,  
                      input$sc1a3sub1, input$sc1a3sub2, 
                      "sc1gexpr.h5", sc1gene, 
                      input$sc1a3siz, input$sc1a3col2, input$sc1a3ord2, 
                      input$sc1a3fsz, input$sc1a3asp, input$sc1a3txt) ) 
  }) 
  output$sc1a3oup2.png <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1a3drX2,"_",input$sc1a3drY2,"_",  
                                   input$sc1a3inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$sc1a3oup2.h, width = input$sc1a3oup2.w, 
      plot = scDRgene(sc1conf, sc1meta, input$sc1a3drX2, input$sc1a3inp2,  
                      input$sc1a3sub1, input$sc1a3sub2, 
                      "sc1gexpr.h5", sc1gene, 
                      input$sc1a3siz, input$sc1a3col2, input$sc1a3ord2, 
                      input$sc1a3fsz, input$sc1a3asp, input$sc1a3txt) ) 
  }) 
     
   
  ### Plots for tab b2 
  output$sc1b2sub1.ui <- renderUI({ 
    sub = strsplit(sc1conf[UI == input$sc1b2sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("sc1b2sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$sc1b2sub1non, { 
    sub = strsplit(sc1conf[UI == input$sc1b2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1b2sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$sc1b2sub1all, { 
    sub = strsplit(sc1conf[UI == input$sc1b2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1b2sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  
  
  output$sc1b2oupTxt <- renderUI({ 
    textAreaInput("sc1b2inp", HTML("List of gene names (genes list, separated by , or ; or newline):"),
                  height = "110px",width = "600px",
                  value = paste0(sc1def$genes, collapse = ", ")
    ) %>%
      helper(type = "inline", size = "m", fade = TRUE,
             title = "List of genes to plot on selected spatial/UMAP plot",
             content = c("Input genes to plot",
                         "",
                         "- Genes should be separated by comma, semicolon or newline"))
    
  })
  
  output$sc1b2oup1 <- renderPlot({ 
    scDRcoex(sc1conf, sc1meta, input$sc1b2drX,    
             input$sc1b2inp, input$sc1b2sub1, input$sc1b2sub2, 
             "sc1gexpr.h5", sc1gene, 
             input$sc1b2siz, input$sc1b2col1, input$sc1b2ord1, 
             input$sc1b2fsz, input$sc1b2asp, input$sc1b2txt) 
  }) 
  output$sc1b2oup1.ui <- renderUI({ 
    plotOutput("sc1b2oup1", height = pList2[input$sc1b2psz]) 
  }) 
  output$sc1b2oup1.pdf <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1b2drX,"_",input$sc1b2drY,"_",  
                                    input$sc1b2inp,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$sc1b2oup1.h, width = input$sc1b2oup1.w, useDingbats = FALSE, 
      plot = scDRcoex(sc1conf, sc1meta, input$sc1b2drX,  
                      input$sc1b2inp, input$sc1b2sub1, input$sc1b2sub2, 
                      "sc1gexpr.h5", sc1gene, 
                      input$sc1b2siz, input$sc1b2col1, input$sc1b2ord1, 
                      input$sc1b2fsz, input$sc1b2asp, input$sc1b2txt) ) 
  }) 
  output$sc1b2oup1.png <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1b2drX,"_",input$sc1b2drY,"_",  
                                    input$sc1b2inp,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$sc1b2oup1.h, width = input$sc1b2oup1.w, 
      plot = scDRcoex(sc1conf, sc1meta, input$sc1b2drX, 
                      input$sc1b2inp, input$sc1b2sub1, input$sc1b2sub2, 
                      "sc1gexpr.h5", sc1gene, 
                      input$sc1b2siz, input$sc1b2col1, input$sc1b2ord1, 
                      input$sc1b2fsz, input$sc1b2asp, input$sc1b2txt) ) 
  }) 
  
     
   
  ### Plots for tab c1 
  output$sc1c1sub1.ui <- renderUI({ 
    sub = strsplit(sc1conf[UI == input$sc1c1sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("sc1c1sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$sc1c1sub1non, { 
    sub = strsplit(sc1conf[UI == input$sc1c1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1c1sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$sc1c1sub1all, { 
    sub = strsplit(sc1conf[UI == input$sc1c1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1c1sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$sc1c1oup <- renderPlot({ 
    scVioBox(sc1conf, sc1meta, input$sc1c1inp1, input$sc1c1inp2, 
             input$sc1c1sub1, input$sc1c1sub2, 
             "sc1gexpr.h5", sc1gene, input$sc1c1typ, input$sc1c1pts, 
             input$sc1c1siz, input$sc1c1fsz) 
  }) 
  output$sc1c1oup.ui <- renderUI({ 
    plotOutput("sc1c1oup", height = pList2[input$sc1c1psz]) 
  }) 
  output$sc1c1oup.pdf <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1c1typ,"_",input$sc1c1inp1,"_",  
                                   input$sc1c1inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$sc1c1oup.h, width = input$sc1c1oup.w, useDingbats = FALSE, 
      plot = scVioBox(sc1conf, sc1meta, input$sc1c1inp1, input$sc1c1inp2, 
                      input$sc1c1sub1, input$sc1c1sub2, 
                      "sc1gexpr.h5", sc1gene, input$sc1c1typ, input$sc1c1pts, 
                      input$sc1c1siz, input$sc1c1fsz) ) 
  }) 
  output$sc1c1oup.png <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1c1typ,"_",input$sc1c1inp1,"_",  
                                   input$sc1c1inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$sc1c1oup.h, width = input$sc1c1oup.w, 
      plot = scVioBox(sc1conf, sc1meta, input$sc1c1inp1, input$sc1c1inp2, 
                      input$sc1c1sub1, input$sc1c1sub2, 
                      "sc1gexpr.h5", sc1gene, input$sc1c1typ, input$sc1c1pts, 
                      input$sc1c1siz, input$sc1c1fsz) ) 
  }) 
     
   
### Plots for tab c2 
  output$sc1c2sub1.ui <- renderUI({
    sub = strsplit(sc1conf[UI == input$sc1c2sub1]$fID, "\\|")[[1]]
    checkboxGroupInput("sc1c2sub2", "Select which cells to show", inline = TRUE,
                       choices = sub, selected = sub)
  })
  observeEvent(input$sc1c2sub1non, { 
    sub = strsplit(sc1conf[UI == input$sc1c2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1c2sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$sc1c2sub1all, { 
    sub = strsplit(sc1conf[UI == input$sc1c2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1c2sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
output$sc1c2oup <- renderPlot({
  scProp(sc1conf, sc1meta, input$sc1c2inp1, input$sc1c2inp2,
         input$sc1c2sub1, input$sc1c2sub2,
         input$sc1c2typ, input$sc1c2flp, input$sc1c2fsz)
})
output$sc1c2oup.ui <- renderUI({
  plotOutput("sc1c2oup", height = pList2[input$sc1c2psz])
})
output$sc1c2oup.pdf <- downloadHandler( 
  filename = function() { paste0("sc1",input$sc1c2typ,"_",input$sc1c2inp1,"_",  
                                 input$sc1c2inp2,".pdf") }, 
  content = function(file) { ggsave( 
    file, device = "pdf", height = input$sc1c2oup.h, width = input$sc1c2oup.w, useDingbats = FALSE, 
    plot = scProp(sc1conf, sc1meta, input$sc1c2inp1, input$sc1c2inp2,  
                  input$sc1c2sub1, input$sc1c2sub2, 
                  input$sc1c2typ, input$sc1c2flp, input$sc1c2fsz) ) 
  }) 
output$sc1c2oup.png <- downloadHandler( 
  filename = function() { paste0("sc1",input$sc1c2typ,"_",input$sc1c2inp1,"_",  
                                 input$sc1c2inp2,".png") }, 
  content = function(file) { ggsave( 
    file, device = "png", height = input$sc1c2oup.h, width = input$sc1c2oup.w, 
    plot = scProp(sc1conf, sc1meta, input$sc1c2inp1, input$sc1c2inp2,  
                  input$sc1c2sub1, input$sc1c2sub2, 
                  input$sc1c2typ, input$sc1c2flp, input$sc1c2fsz) ) 
  }) 
     
   
  ### Plots for tab d1 
  output$sc1d1sub1.ui <- renderUI({ 
    sub = strsplit(sc1conf[UI == input$sc1d1sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("sc1d1sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$sc1d1sub1non, { 
    sub = strsplit(sc1conf[UI == input$sc1d1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1d1sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$sc1d1sub1all, { 
    sub = strsplit(sc1conf[UI == input$sc1d1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1d1sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  

  output$sc1d1oupTxt <- renderUI({ 
    x <- input$sc1d1grp
    # # updateTextAreaInput(session, "sc1d1inp", value = paste0(sc1def[[x]], collapse = ", "))
    textAreaInput("sc1d1inp", HTML("List of gene names (Max 50 genes, separated by , or ; or newline):"),
                  height = "110px",width = "600px",
                  value = paste0(sc1def[[x]], collapse = ", ")
                  # value = paste0(sc1def$genes, collapse = ", ")
    ) %>%
          helper(type = "inline", size = "m", fade = TRUE,
                 title = "List of genes to plot on bubbleplot / heatmap",
                 content = c("Input genes to plot",
                             "- Maximum 50 genes (due to ploting space limitations)",
                             "- Genes should be separated by comma, semicolon or newline"))
    # geneList = scGeneList(input$sc1d1inp, sc1gene)
    # if(nrow(geneList) > 50){
    #   HTML("More than 50 input genes! Please reduce the gene list!")
    # } else {
    #   oup = paste0(nrow(geneList[present == TRUE]), "genes OK and will be plotted")
    #   if(nrow(geneList[present == FALSE]) > 0){
    #     oup = paste0(oup, "<br/>",
    #                  nrow(geneList[present == FALSE]), "genes not found (",
    #                  paste0(geneList[present == FALSE]$gene, collapse = ", "), ")")
    #   }
    # HTML(oup)
    # }
  })
  
  output$sc1d1oup <- renderPlot({ 
    scBubbHeat(sc1conf, sc1meta, input$sc1d1inp, input$sc1d1grp, input$sc1d1plt, 
               input$sc1d1sub1, input$sc1d1sub2, "sc1gexpr.h5", sc1gene, 
               input$sc1d1scl, input$sc1d1row, input$sc1d1col, 
               input$sc1d1cols, input$sc1d1fsz) 
  }) 
  output$sc1d1oup.ui <- renderUI({
    plotOutput("sc1d1oup", height = pList3[input$sc1d1psz])
  })
  output$sc1d1oup.pdf <- downloadHandler(
    filename = function() { paste0("sc1",input$sc1d1plt,"_",input$sc1d1grp,".pdf") }, 
    content = function(file) { 
       
      pdf(file, height = input$sc1d1oup.h, width = input$sc1d1oup.w) 
      print(scBubbHeat(sc1conf, sc1meta, input$sc1d1inp, input$sc1d1grp, input$sc1d1plt, 
                        input$sc1d1sub1, input$sc1d1sub2, "sc1gexpr.h5", sc1gene, 
                        input$sc1d1scl, input$sc1d1row, input$sc1d1col, 
                        input$sc1d1cols, input$sc1d1fsz, save = TRUE))
      dev.off()
  }) 
  output$sc1d1oup.png <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1d1plt,"_",input$sc1d1grp,".png") }, 
    content = function(file) { 
      
      png(file,  width = 465, height = 225, units='mm', res = 300) 
      print(scBubbHeat(sc1conf, sc1meta, input$sc1d1inp, input$sc1d1grp, input$sc1d1plt, 
                        input$sc1d1sub1, input$sc1d1sub2, "sc1gexpr.h5", sc1gene, 
                        input$sc1d1scl, input$sc1d1row, input$sc1d1col, 
                        input$sc1d1cols, input$sc1d1fsz, save = TRUE) )
      dev.off()
  }) 
   
  ################################################################################  
  
  
  # ### Tab1.x2: New tab for genome tracks
 

  choices <- unname(getGenes(ArchRProj = proj)$symbol)


  updateSelectizeInput(session, "sc1trackinp", choices = choices, server = TRUE,
                       selected = choices[1], options = list(
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))


  grps <- c("Clusters","Sample",treatment)
  
  updateSelectizeInput(session, "sc1trackgrp", choices =grps , server = TRUE,
                       selected = grps[1], options = list(
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))
  
  

  observeEvent(input$range_1 , {

    updateNumericInput(session, "range_min_1", value = min(input$range_1))
    updateNumericInput(session, "range_max_1", value = max(input$range_1))

                               }, priority = 200)
  
  
  
  
  
  sctrack <- function(x,y,z,g){
    # Prepare tracks
    upstream <- -min(z)*1000
    downstream <- max(g)*1000
                           
    p <- plotBrowserTrack(
      ArchRProj = proj,
      groupBy = y,
      geneSymbol = x,
      upstream = upstream,
      downstream = downstream,
      loops = getPeak2GeneLinks(proj)
    )
    # Observe the inputs for ATAC-Seq  Co-accessibility
    grid::grid.newpage()

    grid::grid.draw(p[[x]])

  }
  
  
  output$sc1trackoup <- renderPlot({

  sctrack(input$sc1trackinp
          ,input$sc1trackgrp
          ,input$range_min_1,input$range_max_1
            )

  }, height = 550, width = 850)
 
  
  
  output$sc1trackoup.ui <- renderUI({
    plotOutput("sc1trackoup", height = pList[input$sc1a1psz])
  })
  
  

  output$sc1trackoup.pdf <- downloadHandler(
    filename = function() { paste0("sc1","_","tracks",".pdf") },
    content = function(file) {

      pdf(file, height = input$sc1trackoup.h, width = input$sc1trackoup.w)
      print(
        
        
        sctrack(input$sc1trackinp
                ,input$sc1trackgrp
                ,input$range_min_1,input$range_max_1
        )
        
    
        

      )
      dev.off()
    })
  output$sc1trackoup.png <- downloadHandler(
    filename = function() { paste0("sc1","_","tracks",".png") },
    content = function(file) {


      png(file,  width = 465, height = 225, units='mm', res = 300)
      print(
        sctrack(input$sc1trackinp
                ,input$sc1trackgrp
                ,input$range_min_1,input$range_max_1
        )

      )
      dev.off()

    })



  
  

  
  
  
   
   optCrt="{ option_create: function(data,escape) {return('<div class=\"create\"><strong>' + '</strong></div>');} }" 
  updateSelectizeInput(session, "sc2a1inp2", choices = names(sc2gene), server = TRUE,
                       selected = sc2def$gene1, options = list(
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))
  updateSelectizeInput(session, "sc2a3inp1", choices = names(sc2gene), server = TRUE,
                       selected = sc2def$gene1, options = list(
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))
  updateSelectizeInput(session, "sc2a3inp2", choices = names(sc2gene), server = TRUE,
                       selected = sc2def$gene2, options = list(
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))
  updateSelectizeInput(session, "sc2b2inp1", choices = names(sc2gene), server = TRUE,
                       selected = sc2def$gene1, options = list(
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))
  updateSelectizeInput(session, "sc2b2inp2", choices = names(sc2gene), server = TRUE,
                       selected = sc2def$gene2, options = list(
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))
  updateSelectizeInput(session, "sc2c1inp2", server = TRUE,
                       choices = c(sc2conf[is.na(fID)]$UI,names(sc2gene)),
                       selected = sc2conf[is.na(fID)]$UI[1], options = list(
                         maxOptions = length(sc2conf[is.na(fID)]$UI) + 3,
                         create = TRUE, persist = TRUE, render = I(optCrt)))

  
  
  ##### START OF NEW CODE TO ADD ####################
  ### Tab1.x1: New tab for volcano
  
 
  ##############################  
  
  
  ### Tab1.x2: New tab for Motif logo
  inplogos = readRDS("seqlogo.rds")

  updateSelectizeInput(session, "sc2logoinp", choices = names(inplogos), server = TRUE,
                       selected = names(inplogos)[1], options = list(
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt)))


  sclogo <- function(inplogos){
    # Prepare logo
    seqlogo = ggseqlogo(inplogos[c(grep(paste0("^",input$sc2logoinp,"$"), names(inplogos)))])
    return(seqlogo)
  }
  output$sc2logooup <- renderPlot({
    sclogo(inplogos)

  }, height = 450, width = 750)

  output$sc2logooup.ui <- renderUI({
    plotOutput("sc2logooup", height = pList[input$sc1a1psz])
  })
  output$sc2logooup.pdf <- downloadHandler(
    filename = function() { paste0("sc1","_","logo",".pdf") },
    content = function(file) { ggsave(
      file, device = "pdf", height = input$sc2logooup.h, width = input$sc2logooup.w, useDingbats = FALSE,
      plot = sclogo(inplogos) )
    })
  output$sc2logooup.png <- downloadHandler(
    filename = function() { paste0("sc1","_","logo",".png") },
    content = function(file) { ggsave(
      file, device = "png", height = input$sc2logooup.h, width = input$sc2logooup.w,
      plot = sclogo(inplogos) )
    })
  
  ##### END OF NEW CODE TO ADD #########################
  
  
  ### Plots for tab a1 
  output$sc1a1sub1.ui <- renderUI({ 
    sub = strsplit(sc1conf[UI == input$sc1a1sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("sc1a1sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$sc1a1sub1non, { 
    sub = strsplit(sc1conf[UI == input$sc1a1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1a1sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$sc1a1sub1all, { 
    sub = strsplit(sc1conf[UI == input$sc1a1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1a1sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  
  # add color change for clusters
  # the group levels
  dcn <- reactive({strsplit(sc1conf[UI == input$sc1a1inp1]$fID, "\\|")[[1]] })
  # default colors
  defCol <- reactive({
    if(file.exists(paste0('colorset_',input$sc1a1inp1,'.csv'))==TRUE){
      read.csv(paste0('colorset_',input$sc1a1inp1,'.csv'))$x
    }else{
      write.csv(colorRampPalette(brewer.pal(12, "Paired"))(length(dcn())),paste0('colorset_',input$sc1a1inp1,'.csv'))
      
      toMatch <- c("Clusters", "Sample", "condition_")
      choices <- unique (grep(paste(toMatch,collapse="|"), colnames(getCellColData(ArchRProj = proj)), value=TRUE))
      
      for(choice in choices){
        write.csv(colorRampPalette(brewer.pal(12, "Paired"))(length(unique(getCellColData(ArchRProj = proj)[[choice]]))),paste0('colorset_',choice,'.csv'))
      }
      colorRampPalette(brewer.pal(12, "Paired"))(length(dcn()))
      
    }
  })
  
  
  
  
  output$sc1a1sub3.ui <- renderUI({
    
    colInput <- function(vecFeats) { # vecFeats = vector of feature names for colors
      pickers <- (lapply(1:length(vecFeats), function(k) {
        colourInput(
          inputId = paste0("col", k)
          ,label = dcn()[k]      # color sel label for user
          ,if(file.exists(paste0('colorset_',input$sc1a1inp1,'.csv'))==TRUE){
            value = read.csv(paste0('colorset_',input$sc1a1inp1,'.csv'))$x[k]
          }else{
            value = defCol()[k]}#  this should match initial plot
          ,showColour = "both"  # show hex and color itself
          # ,width = 2
        )})) }
    
    
    #   create color selectors for plot
    
    dashboardSidebar(  # create mask for user interaction
      collapsed = TRUE
      , title = "Choose colors for the plot."  # sidebar title
      # ,.list = colInput(dcn())   # if you want to have color pallet in one column
      
      # if you want to have color pallet in 4 columns, the below code is temporarily
      
      ,tags$style(HTML(".main-sidebar { font-size: 0px; }"))
      
      ,splitLayout( tryCatch(colInput(dcn())[[1]],error = function(e){''})
                    ,tryCatch(colInput(dcn())[[2]],error = function(e){''})
                    ,tryCatch(colInput(dcn())[[3]],error = function(e){''})
                    ,tryCatch(colInput(dcn())[[4]],error = function(e){''})
                    ,tryCatch(colInput(dcn())[[5]],error = function(e){''})
                    ,tryCatch(colInput(dcn())[[6]],error = function(e){''})
                    ,tryCatch(colInput(dcn())[[7]],error = function(e){''})
                    ,tryCatch(colInput(dcn())[[8]],error = function(e){''})
                    ,tryCatch(colInput(dcn())[[9]],error = function(e){''})
                    ,tryCatch(colInput(dcn())[[10]],error = function(e){''})
                    ,tryCatch(colInput(dcn())[[11]],error = function(e){''})
                    ,tryCatch(colInput(dcn())[[12]],error = function(e){''})
                    ,cellArgs = list (style = "overflow:visible; width: 88px") #; padding: 15px
                    , align = "left")
      
      ,splitLayout( tryCatch(colInput(dcn())[[13]],error = function(e){''})
                    ,tryCatch(colInput(dcn())[[14]],error = function(e){''})
                    ,tryCatch(colInput(dcn())[[15]],error = function(e){''})
                    ,tryCatch(colInput(dcn())[[16]],error = function(e){''})
                    ,tryCatch(colInput(dcn())[[17]],error = function(e){''})
                    ,tryCatch(colInput(dcn())[[18]],error = function(e){''})
                    ,tryCatch(colInput(dcn())[[19]],error = function(e){''})
                    ,tryCatch(colInput(dcn())[[20]],error = function(e){''})
                    ,tryCatch(colInput(dcn())[[21]],error = function(e){''})
                    ,tryCatch(colInput(dcn())[[22]],error = function(e){''})
                    ,tryCatch(colInput(dcn())[[23]],error = function(e){''})
                    ,tryCatch(colInput(dcn())[[24]],error = function(e){''})
                    ,cellArgs = list (style = "overflow:visible; width: 88px")
                    , align = "left")
      
      ,splitLayout( tryCatch(colInput(dcn())[[25]],error = function(e){''})
                    ,tryCatch(colInput(dcn())[[26]],error = function(e){''})
                    ,tryCatch(colInput(dcn())[[27]],error = function(e){''})
                    ,tryCatch(colInput(dcn())[[28]],error = function(e){''})
                    ,tryCatch(colInput(dcn())[[29]],error = function(e){''})
                    ,tryCatch(colInput(dcn())[[30]],error = function(e){''})
                    ,tryCatch(colInput(dcn())[[31]],error = function(e){''})
                    ,tryCatch(colInput(dcn())[[32]],error = function(e){''})
                    ,tryCatch(colInput(dcn())[[33]],error = function(e){''})
                    ,tryCatch(colInput(dcn())[[34]],error = function(e){''})
                    ,tryCatch(colInput(dcn())[[35]],error = function(e){''})
                    ,tryCatch(colInput(dcn())[[36]],error = function(e){''})
                    ,cellArgs = list (style = "overflow:visible; width: 88px")
                    , align = "left")
      
      ,splitLayout( tryCatch(colInput(dcn())[[37]],error = function(e){''})
                    ,tryCatch(colInput(dcn())[[38]],error = function(e){''})
                    ,tryCatch(colInput(dcn())[[39]],error = function(e){''})
                    ,tryCatch(colInput(dcn())[[40]],error = function(e){''})
                    ,cellArgs = list (style = "overflow:visible; width: 88px")
                    , align = "left")
      ,actionButton("KeepMyColor", "KeepMyColor")
      ,actionButton("DefaultColor", "DefaultColor")
      
      , write.csv(
        colorRampPalette(brewer.pal(12, "Paired"))(length(dcn()))
        ,paste0('colorset_',input$sc1a1inp1,'.csv'))
      
    ) #dashboardSidebar
    
    
  }) #renderUI
  
  
  observeEvent(input$KeepMyColor, {
    write.csv(
      paste0("input$col", 1:length(dcn())) %>%   # use data color names
        map(., function(i) {eval(parse(text = i))}) %>% # convert strings to obj
        unlist() %>% setNames(dcn())
      ,paste0('colorset_',input$sc1a1inp1,'.csv'))
    
    lapply(1:length(dcn()), function(k) {
      updateColourInput(session,
                        inputId = paste0("col", k)
                        ,label = dcn()[k]      # color sel label for user
                        ,value = read.csv(paste0('colorset_',input$sc1a1inp1,'.csv'))$x[k] #  this should match initial plot
                        ,showColour = "both"  # show hex and color itself
                        # ,width = 2
      )})
    
    
  }, ignoreInit=TRUE
  )
  
  observeEvent(input$DefaultColor, {
    write.csv(
      colorRampPalette(brewer.pal(12, "Paired"))(length(dcn()))
      ,paste0('colorset_',input$sc1a1inp1,'.csv'))
    
    # update colors to default colors
    
    lapply(1:length(dcn()), function(k) {
      updateColourInput(session,
                        inputId = paste0("col", k)
                        ,label = dcn()[k]      # color sel label for user
                        ,value = colorRampPalette(brewer.pal(12, "Paired"))(length(dcn()))[k] #  this should match initial plot
                        ,showColour = "both"  # show hex and color itself
                        # ,width = 2
      )})
    
    
  }, ignoreInit=TRUE
  )
  
  output$sc1a1oup1 <- renderPlot({ 
    scDRcell(sc1conf, sc1meta, input$sc1a1drX, input$sc1a1inp1,  
             input$sc1a1sub1, input$sc1a1sub2,
             paste0("input$col", 1:length(dcn())) %>%   # use data color names
               map(., function(i) {eval(parse(text = i))}) %>% # convert strings to obj
               unlist() %>%  setNames(dcn()),
             input$sc1a1siz, input$sc1a1col1, input$sc1a1ord1, 
             input$sc1a1fsz, input$sc1a1asp, input$sc1a1txt, input$sc1a1lab1) 
  }) 
  output$sc1a1oup1.ui <- renderUI({ 
    plotOutput("sc1a1oup1", height = pList[input$sc1a1psz]) 
  }) 
  output$sc1a1oup1.pdf <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1a1drX,"_",input$sc1a1drY,"_",  
                                   input$sc1a1inp1,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$sc1a1oup1.h, width = input$sc1a1oup1.w, useDingbats = FALSE, 
      plot = scDRcell(sc1conf, sc1meta, input$sc1a1drX, input$sc1a1inp1,   
                      input$sc1a1sub1, input$sc1a1sub2, 
                      paste0("input$col", 1:length(dcn())) %>%   # use data color names
                        map(., function(i) {eval(parse(text = i))}) %>% # convert strings to obj
                        unlist() %>%  setNames(dcn()),
                      input$sc1a1siz, input$sc1a1col1, input$sc1a1ord1,  
                      input$sc1a1fsz, input$sc1a1asp, input$sc1a1txt, input$sc1a1lab1) ) 
    }) 
  output$sc1a1oup1.png <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1a1drX,"_",input$sc1a1drY,"_",  
                                   input$sc1a1inp1,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$sc1a1oup1.h, width = input$sc1a1oup1.w, 
      plot = scDRcell(sc1conf, sc1meta, input$sc1a1drX, input$sc1a1inp1,   
                      input$sc1a1sub1, input$sc1a1sub2,
                      paste0("input$col", 1:length(dcn())) %>%   # use data color names
                        map(., function(i) {eval(parse(text = i))}) %>% # convert strings to obj
                        unlist() %>%  setNames(dcn()),
                      input$sc1a1siz, input$sc1a1col1, input$sc1a1ord1,  
                      input$sc1a1fsz, input$sc1a1asp, input$sc1a1txt, input$sc1a1lab1) ) 
    }) 
  output$sc1a1.dt <- renderDataTable({ 
    ggData = scDRnum(sc1conf, sc1meta, input$sc1a1inp1, input$sc1a1inp2, 
                     input$sc1a1sub1, input$sc1a1sub2, 
                     "sc1gexpr.h5", sc1gene, input$sc1a1splt) 
    datatable(ggData, rownames = FALSE, extensions = "Buttons", 
              options = list(pageLength = -1, dom = "tB", buttons = c("copy", "csv", "excel"))) %>% 
      formatRound(columns = c("pctExpress"), digits = 2) 
  }) 
  
  output$sc1a1oup2 <- renderPlot({ 
    scDRgene(sc1conf, sc1meta, input$sc1a1drX2, input$sc1a1inp2,  
             input$sc1a1sub1, input$sc1a1sub2, 
             "sc1gexpr.h5", sc1gene, 
             input$sc1a1siz, input$sc1a1col2, input$sc1a1ord2, 
             input$sc1a1fsz, input$sc1a1asp, input$sc1a1txt) 
  }) 
  output$sc1a1oup2.ui <- renderUI({ 
    plotOutput("sc1a1oup2", height = pList[input$sc1a1psz]) 
  }) 
  output$sc1a1oup2.pdf <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1a1drX2,"_",input$sc1a1drY2,"_",  
                                   input$sc1a1inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$sc1a1oup2.h, width = input$sc1a1oup2.w, useDingbats = FALSE, 
      plot = scDRgene(sc1conf, sc1meta, input$sc1a1drX2, input$sc1a1inp2,  
                      input$sc1a1sub1, input$sc1a1sub2, 
                      "sc1gexpr.h5", sc1gene, 
                      input$sc1a1siz, input$sc1a1col2, input$sc1a1ord2, 
                      input$sc1a1fsz, input$sc1a1asp, input$sc1a1txt) ) 
    }) 
  output$sc1a1oup2.png <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1a1drX2,"_",input$sc1a1drY2,"_",  
                                   input$sc1a1inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$sc1a1oup2.h, width = input$sc1a1oup2.w, 
      plot = scDRgene(sc1conf, sc1meta, input$sc1a1drX2, input$sc1a1inp2,  
                      input$sc1a1sub1, input$sc1a1sub2, 
                      "sc1gexpr.h5", sc1gene, 
                      input$sc1a1siz, input$sc1a1col2, input$sc1a1ord2, 
                      input$sc1a1fsz, input$sc1a1asp, input$sc1a1txt) ) 
    }) 
  
  
  ### Plots for tab a2 
  output$sc1a2sub1.ui <- renderUI({ 
    sub = strsplit(sc1conf[UI == input$sc1a2sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("sc1a2sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$sc1a2sub1non, { 
    sub = strsplit(sc1conf[UI == input$sc1a2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1a2sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$sc1a2sub1all, { 
    sub = strsplit(sc1conf[UI == input$sc1a2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1a2sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  
  
  observeEvent(input$sc1a2togL, {
    
    lapply(1:length(dcn()), function(k) {
      updateColourInput(session,
                        inputId = paste0("col", k)
                        ,label = dcn()[k]      # color sel label for user
                        ,value = read.csv(paste0('colorset_',input$sc1a2inp1,'.csv'))$x[k] #  this should match initial plot
                        ,showColour = "both"  # show hex and color itself
                        # ,width = 2
      )})
    # print(read.csv(paste0('colorset_',input$sc1a2inp1,'.csv'))$x)
    
    
  }, ignoreInit=TRUE
  )
  
  dcn2 <- reactive({strsplit(sc1conf[UI == input$sc1a2inp1]$fID, "\\|")[[1]] })
  output$sc1a2oup1 <- renderPlot({ 
    
    
    scDRcell(sc1conf, sc1meta, input$sc1a2drX, input$sc1a2inp1,  
             input$sc1a2sub1, input$sc1a2sub2, 
             # read.csv(paste0('colorset_',input$sc1a2inp1,'.csv'))$x,
             paste0("input$col", 1:length(dcn2())) %>%   # use data color names
               map(., function(i) {eval(parse(text = i))}) %>% # convert strings to obj
               unlist() %>%  setNames(dcn2()),
             input$sc1a2siz, input$sc1a2col1, input$sc1a2ord1, 
             input$sc1a2fsz, input$sc1a2asp, input$sc1a2txt, input$sc1a2lab1) 
  }) 
  output$sc1a2oup1.ui <- renderUI({ 
    plotOutput("sc1a2oup1", height = pList[input$sc1a2psz]) 
  }) 
  output$sc1a2oup1.pdf <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1a2drX,"_",input$sc1a2drY,"_",  
                                   input$sc1a2inp1,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$sc1a2oup1.h, width = input$sc1a2oup1.w, useDingbats = FALSE, 
      plot = scDRcell(sc1conf, sc1meta, input$sc1a2drX, input$sc1a2inp1,   
                      input$sc1a2sub1, input$sc1a2sub2, 
                      paste0("input$col", 1:length(dcn2())) %>%   # use data color names
                        map(., function(i) {eval(parse(text = i))}) %>% # convert strings to obj
                        unlist() %>%  setNames(dcn2()),
                      input$sc1a2siz, input$sc1a2col1, input$sc1a2ord1,  
                      input$sc1a2fsz, input$sc1a2asp, input$sc1a2txt, input$sc1a2lab1) ) 
    }) 
  output$sc1a2oup1.png <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1a2drX,"_",input$sc1a2drY,"_",  
                                   input$sc1a2inp1,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$sc1a2oup1.h, width = input$sc1a2oup1.w, 
      plot = scDRcell(sc1conf, sc1meta, input$sc1a2drX, input$sc1a2inp1,   
                      input$sc1a2sub1, input$sc1a2sub2,
                      paste0("input$col", 1:length(dcn2())) %>%   # use data color names
                        map(., function(i) {eval(parse(text = i))}) %>% # convert strings to obj
                        unlist() %>%  setNames(dcn2()),
                      input$sc1a2siz, input$sc1a2col1, input$sc1a2ord1,  
                      input$sc1a2fsz, input$sc1a2asp, input$sc1a2txt, input$sc1a2lab1) ) 
    }) 
  
  dcn3 <- reactive({strsplit(sc1conf[UI == input$sc1a2inp2]$fID, "\\|")[[1]] })
  output$sc1a2oup2 <- renderPlot({ 
    scDRcell(sc1conf, sc1meta, input$sc1a2drX2, input$sc1a2inp2,  
             input$sc1a2sub1, input$sc1a2sub2, 
             paste0("input$col", 1:length(dcn3())) %>%   # use data color names
               map(., function(i) {eval(parse(text = i))}) %>% # convert strings to obj
               unlist() %>%  setNames(dcn3()),
             input$sc1a2siz, input$sc1a2col2, input$sc1a2ord2, 
             input$sc1a2fsz, input$sc1a2asp, input$sc1a2txt, input$sc1a2lab2) 
  }) 
  output$sc1a2oup2.ui <- renderUI({ 
    plotOutput("sc1a2oup2", height = pList[input$sc1a2psz]) 
  }) 
  output$sc1a2oup2.pdf <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1a2drX2,"_",input$sc1a2drY2,"_",  
                                   input$sc1a2inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$sc1a2oup2.h, width = input$sc1a2oup2.w, useDingbats = FALSE, 
      plot = scDRcell(sc1conf, sc1meta, input$sc1a2drX2, input$sc1a2inp2,   
                      input$sc1a2sub1, input$sc1a2sub2,
                      paste0("input$col", 1:length(dcn3())) %>%   # use data color names
                        map(., function(i) {eval(parse(text = i))}) %>% # convert strings to obj
                        unlist() %>%  setNames(dcn3()),
                      input$sc1a2siz, input$sc1a2col2, input$sc1a2ord2,  
                      input$sc1a2fsz, input$sc1a2asp, input$sc1a2txt, input$sc1a2lab2) ) 
    }) 
  output$sc1a2oup2.png <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1a2drX2,"_",input$sc1a2drY2,"_",  
                                   input$sc1a2inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$sc1a2oup2.h, width = input$sc1a2oup2.w, 
      plot = scDRcell(sc1conf, sc1meta, input$sc1a2drX2, input$sc1a2inp2,   
                      input$sc1a2sub1, input$sc1a2sub2, 
                      paste0("input$col", 1:length(dcn3())) %>%   # use data color names
                        map(., function(i) {eval(parse(text = i))}) %>% # convert strings to obj
                        unlist() %>%  setNames(dcn3()),
                      input$sc1a2siz, input$sc1a2col2, input$sc1a2ord2,  
                      input$sc1a2fsz, input$sc1a2asp, input$sc1a2txt, input$sc1a2lab2) ) 
    }) 
  
   
  ### Plots for tab a3 
  output$sc2a3sub1.ui <- renderUI({
    sub = strsplit(sc2conf[UI == input$sc2a3sub1]$fID, "\\|")[[1]]
    checkboxGroupInput("sc2a3sub2", "Select which cells to show", inline = TRUE,
                       choices = sub, selected = sub)
  })
  observeEvent(input$sc2a3sub1non, {
    sub = strsplit(sc2conf[UI == input$sc2a3sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "sc2a3sub2", label = "Select which cells to show",
                             choices = sub, selected = NULL, inline = TRUE)
  })
  observeEvent(input$sc2a3sub1all, {
    sub = strsplit(sc2conf[UI == input$sc2a3sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "sc2a3sub2", label = "Select which cells to show",
                             choices = sub, selected = sub, inline = TRUE)
  })
  output$sc2a3oup1 <- renderPlot({
    scDRgene(sc2conf, sc2meta, input$sc2a3drX, input$sc2a3inp1,
             input$sc2a3sub1, input$sc2a3sub2,
             "sc2gexpr.h5", sc2gene,
             input$sc2a3siz, input$sc2a3col1, input$sc2a3ord1,
             input$sc2a3fsz, input$sc2a3asp, input$sc2a3txt)
  })
  output$sc2a3oup1.ui <- renderUI({
    plotOutput("sc2a3oup1", height = pList[input$sc2a3psz])
  })
  output$sc2a3oup1.pdf <- downloadHandler(
    filename = function() { paste0("sc2",input$sc2a3drX,"_",input$sc2a3drY,"_",
                                   input$sc2a3inp1,".pdf") },
    content = function(file) { ggsave(
      file, device = "pdf", height = input$sc2a3oup1.h, width = input$sc2a3oup1.w, useDingbats = FALSE,
      plot = scDRgene(sc2conf, sc2meta, input$sc2a3drX, input$sc2a3inp1,
                      input$sc2a3sub1, input$sc2a3sub2,
                      "sc2gexpr.h5", sc2gene,
                      input$sc2a3siz, input$sc2a3col1, input$sc2a3ord1,
                      input$sc2a3fsz, input$sc2a3asp, input$sc2a3txt) )
  })
  output$sc2a3oup1.png <- downloadHandler(
    filename = function() { paste0("sc2",input$sc2a3drX,"_",input$sc2a3drY,"_",
                                   input$sc2a3inp1,".png") },
    content = function(file) { ggsave(
      file, device = "png", height = input$sc2a3oup1.h, width = input$sc2a3oup1.w,
      plot = scDRgene(sc2conf, sc2meta, input$sc2a3drX, input$sc2a3inp1,
                      input$sc2a3sub1, input$sc2a3sub2,
                      "sc2gexpr.h5", sc2gene,
                      input$sc2a3siz, input$sc2a3col1, input$sc2a3ord1,
                      input$sc2a3fsz, input$sc2a3asp, input$sc2a3txt) )
  })

  output$sc2a3oup2 <- renderPlot({
    scDRgene(sc2conf, sc2meta, input$sc2a3drX2, input$sc2a3inp2,
             input$sc2a3sub1, input$sc2a3sub2,
             "sc2gexpr.h5", sc2gene,
             input$sc2a3siz, input$sc2a3col2, input$sc2a3ord2,
             input$sc2a3fsz, input$sc2a3asp, input$sc2a3txt)
  })
  output$sc2a3oup2.ui <- renderUI({
    plotOutput("sc2a3oup2", height = pList[input$sc2a3psz])
  })
  output$sc2a3oup2.pdf <- downloadHandler(
    filename = function() { paste0("sc2",input$sc2a3drX2,"_",input$sc2a3drY2,"_",
                                   input$sc2a3inp2,".pdf") },
    content = function(file) { ggsave(
      file, device = "pdf", height = input$sc2a3oup2.h, width = input$sc2a3oup2.w, useDingbats = FALSE,
      plot = scDRgene(sc2conf, sc2meta, input$sc2a3drX2, input$sc2a3inp2,
                      input$sc2a3sub1, input$sc2a3sub2,
                      "sc2gexpr.h5", sc2gene,
                      input$sc2a3siz, input$sc2a3col2, input$sc2a3ord2,
                      input$sc2a3fsz, input$sc2a3asp, input$sc2a3txt) )
  })
  output$sc2a3oup2.png <- downloadHandler(
    filename = function() { paste0("sc2",input$sc2a3drX2,"_",input$sc2a3drY2,"_",
                                   input$sc2a3inp2,".png") },
    content = function(file) { ggsave(
      file, device = "png", height = input$sc2a3oup2.h, width = input$sc2a3oup2.w,
      plot = scDRgene(sc2conf, sc2meta, input$sc2a3drX2, input$sc2a3inp2,
                      input$sc2a3sub1, input$sc2a3sub2,
                      "sc2gexpr.h5", sc2gene,
                      input$sc2a3siz, input$sc2a3col2, input$sc2a3ord2,
                      input$sc2a3fsz, input$sc2a3asp, input$sc2a3txt) )
  })
  #    
   
  ### Plots for tab b2 
  output$sc2b2sub1.ui <- renderUI({
    sub = strsplit(sc2conf[UI == input$sc2b2sub1]$fID, "\\|")[[1]]
    checkboxGroupInput("sc2b2sub2", "Select which cells to show", inline = TRUE,
                       choices = sub, selected = sub)
  })
  observeEvent(input$sc2b2sub1non, {
    sub = strsplit(sc2conf[UI == input$sc2b2sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "sc2b2sub2", label = "Select which cells to show",
                             choices = sub, selected = NULL, inline = TRUE)
  })
  observeEvent(input$sc2b2sub1all, {
    sub = strsplit(sc2conf[UI == input$sc2b2sub1]$fID, "\\|")[[1]]
    updateCheckboxGroupInput(session, inputId = "sc2b2sub2", label = "Select which cells to show",
                             choices = sub, selected = sub, inline = TRUE)
  })
  output$sc2b2oup1 <- renderPlot({
    scDRcoex(sc2conf, sc2meta, input$sc2b2drX, 
             input$sc2b2inp1, input$sc2b2inp2, input$sc2b2sub1, input$sc2b2sub2,
             "sc2gexpr.h5", sc2gene,
             input$sc2b2siz, input$sc2b2col1, input$sc2b2ord1,
             input$sc2b2fsz, input$sc2b2asp, input$sc2b2txt)
  })
  output$sc2b2oup1.ui <- renderUI({
    plotOutput("sc2b2oup1", height = pList2[input$sc2b2psz])
  })
  output$sc2b2oup1.pdf <- downloadHandler(
    filename = function() { paste0("sc2",input$sc2b2drX,"_",input$sc2b2drY,"_",
                                    input$sc2b2inp1,"_",input$sc2b2inp2,".pdf") },
    content = function(file) { ggsave(
      file, device = "pdf", height = input$sc2b2oup1.h, width = input$sc2b2oup1.w, useDingbats = FALSE,
      plot = scDRcoex(sc2conf, sc2meta, input$sc2b2drX, 
                      input$sc2b2inp1, input$sc2b2inp2, input$sc2b2sub1, input$sc2b2sub2,
                      "sc2gexpr.h5", sc2gene,
                      input$sc2b2siz, input$sc2b2col1, input$sc2b2ord1,
                      input$sc2b2fsz, input$sc2b2asp, input$sc2b2txt) )
  })
  output$sc2b2oup1.png <- downloadHandler(
    filename = function() { paste0("sc2",input$sc2b2drX,"_",input$sc2b2drY,"_",
                                    input$sc2b2inp1,"_",input$sc2b2inp2,".png") },
    content = function(file) { ggsave(
      file, device = "png", height = input$sc2b2oup1.h, width = input$sc2b2oup1.w,
      plot = scDRcoex(sc2conf, sc2meta, input$sc2b2drX, 
                      input$sc2b2inp1, input$sc2b2inp2, input$sc2b2sub1, input$sc2b2sub2,
                      "sc2gexpr.h5", sc2gene,
                      input$sc2b2siz, input$sc2b2col1, input$sc2b2ord1,
                      input$sc2b2fsz, input$sc2b2asp, input$sc2b2txt) )
  })
  output$sc2b2oup2 <- renderPlot({
    scDRcoexLeg(input$sc2b2inp1, input$sc2b2inp2, input$sc2b2col1, input$sc2b2fsz)
  })
  output$sc2b2oup2.ui <- renderUI({
    plotOutput("sc2b2oup2", height = "300px")
  })
  output$sc2b2oup2.pdf <- downloadHandler(
    filename = function() { paste0("sc2",input$sc2b2drX,"_",input$sc2b2drY,"_",
                                    input$sc2b2inp1,"_",input$sc2b2inp2,"_leg.pdf") },
    content = function(file) { ggsave(
      file, device = "pdf", height = 3, width = 4, useDingbats = FALSE,
      plot = scDRcoexLeg(input$sc2b2inp1, input$sc2b2inp2, input$sc2b2col1, input$sc2b2fsz) )
  })
  output$sc2b2oup2.png <- downloadHandler(
    filename = function() { paste0("sc2",input$sc2b2drX,"_",input$sc2b2drY,"_",
                                    input$sc2b2inp1,"_",input$sc2b2inp2,"_leg.png") },
    content = function(file) { ggsave(
      file, device = "png", height = 3, width = 4,
      plot = scDRcoexLeg(input$sc2b2inp1, input$sc2b2inp2, input$sc2b2col1, input$sc2b2fsz) )
  })
  output$sc2b2.dt <- renderDataTable({
    ggData = scDRcoexNum(sc2conf, sc2meta, input$sc2b2inp1, input$sc2b2inp2,
                         input$sc2b2sub1, input$sc2b2sub2, "sc2gexpr.h5", sc2gene)
    datatable(ggData, rownames = FALSE, extensions = "Buttons",
              options = list(pageLength = -1, dom = "tB", buttons = c("copy", "csv", "excel"))) %>%
      formatRound(columns = c("percent"), digits = 2)
  })

   
  ### Plots for tab c1 
  output$sc1c1sub1.ui <- renderUI({ 
    sub = strsplit(sc1conf[UI == input$sc1c1sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("sc1c1sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$sc1c1sub1non, { 
    sub = strsplit(sc1conf[UI == input$sc1c1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1c1sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$sc1c1sub1all, { 
    sub = strsplit(sc1conf[UI == input$sc1c1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1c1sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  
  observeEvent(input$sc1c1togL, {
    output$sc1c1oup <- renderPlot({ 
      scVioBox(sc1conf, sc1meta, input$sc1c1inp1, input$sc1c1inp2, 
               input$sc1c1sub1, input$sc1c1sub2,
               read.csv(paste0('colorset_',input$sc1c1inp1,'.csv'))$x,
               "sc1gexpr.h5", sc1gene, input$sc1c1typ, input$sc1c1pts, 
               input$sc1c1siz, input$sc1c1fsz) 
    })  
    
    
  }, ignoreInit=TRUE
  )
  
  output$sc1c1oup <- renderPlot({ 
    scVioBox(sc1conf, sc1meta, input$sc1c1inp1, input$sc1c1inp2, 
             input$sc1c1sub1, input$sc1c1sub2, 
             read.csv(paste0('colorset_',input$sc1c1inp1,'.csv'))$x,
             "sc1gexpr.h5", sc1gene, input$sc1c1typ, input$sc1c1pts, 
             input$sc1c1siz, input$sc1c1fsz) 
  }) 
  output$sc1c1oup.ui <- renderUI({ 
    plotOutput("sc1c1oup", height = pList2[input$sc1c1psz]) 
  }) 
  output$sc1c1oup.pdf <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1c1typ,"_",input$sc1c1inp1,"_",  
                                   input$sc1c1inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$sc1c1oup.h, width = input$sc1c1oup.w, useDingbats = FALSE, 
      plot = scVioBox(sc1conf, sc1meta, input$sc1c1inp1, input$sc1c1inp2, 
                      input$sc1c1sub1, input$sc1c1sub2,
                      read.csv(paste0('colorset_',input$sc1c1inp1,'.csv'))$x,
                      "sc1gexpr.h5", sc1gene, input$sc1c1typ, input$sc1c1pts, 
                      input$sc1c1siz, input$sc1c1fsz) ) 
    }) 
  output$sc1c1oup.png <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1c1typ,"_",input$sc1c1inp1,"_",  
                                   input$sc1c1inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$sc1c1oup.h, width = input$sc1c1oup.w, 
      plot = scVioBox(sc1conf, sc1meta, input$sc1c1inp1, input$sc1c1inp2, 
                      input$sc1c1sub1, input$sc1c1sub2, 
                      read.csv(paste0('colorset_',input$sc1c1inp1,'.csv'))$x,
                      "sc1gexpr.h5", sc1gene, input$sc1c1typ, input$sc1c1pts, 
                      input$sc1c1siz, input$sc1c1fsz) ) 
    }) 
  
  
  ### Plots for tab c2 
  
  output$sc1c2sub1.ui <- renderUI({
    sub = strsplit(sc1conf[UI == input$sc1c2sub1]$fID, "\\|")[[1]]
    checkboxGroupInput("sc1c2sub2", "Select which cells to show", inline = TRUE,
                       choices = sub, selected = sub)
  })
  observeEvent(input$sc1c2sub1non, { 
    sub = strsplit(sc1conf[UI == input$sc1c2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1c2sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$sc1c2sub1all, { 
    sub = strsplit(sc1conf[UI == input$sc1c2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "sc1c2sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  
  observeEvent(input$sc1c2togL, {
    
    output$sc1c2oup <- renderPlot({
      scProp(sc1conf, sc1meta, input$sc1c2inp1, input$sc1c2inp2,
             input$sc1c2sub1, input$sc1c2sub2,
             read.csv(paste0('colorset_',input$sc1c2inp2,'.csv'))$x,
             input$sc1c2typ, input$sc1c2flp, input$sc1c2fsz)
    })
    
    
    
  }, ignoreInit=TRUE
  )
  output$sc1c2oup <- renderPlot({
    scProp(sc1conf, sc1meta, input$sc1c2inp1, input$sc1c2inp2,
           input$sc1c2sub1, input$sc1c2sub2,
           read.csv(paste0('colorset_',input$sc1c2inp2,'.csv'))$x,
           input$sc1c2typ, input$sc1c2flp, input$sc1c2fsz)
  })
  output$sc1c2oup.ui <- renderUI({
    plotOutput("sc1c2oup", height = pList2[input$sc1c2psz])
  })
  output$sc1c2oup.pdf <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1c2typ,"_",input$sc1c2inp1,"_",  
                                   input$sc1c2inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$sc1c2oup.h, width = input$sc1c2oup.w, useDingbats = FALSE, 
      plot = scProp(sc1conf, sc1meta, input$sc1c2inp1, input$sc1c2inp2,
                    input$sc1c2sub1, input$sc1c2sub2,
                    read.csv(paste0('colorset_',input$sc1c2inp2,'.csv'))$x,
                    input$sc1c2typ, input$sc1c2flp, input$sc1c2fsz) ) 
    }) 
  output$sc1c2oup.png <- downloadHandler( 
    filename = function() { paste0("sc1",input$sc1c2typ,"_",input$sc1c2inp1,"_",  
                                   input$sc1c2inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$sc1c2oup.h, width = input$sc1c2oup.w, 
      plot = scProp(sc1conf, sc1meta, input$sc1c2inp1, input$sc1c2inp2,  
                    input$sc1c2sub1, input$sc1c2sub2,
                    read.csv(paste0('colorset_',input$sc1c2inp2,'.csv'))$x,
                    input$sc1c2typ, input$sc1c2flp, input$sc1c2fsz) ) 
    }) 
  

#   ### Plots for tab d1 
output$sc2d1sub1.ui <- renderUI({
  sub = strsplit(sc2conf[UI == input$sc2d1sub1]$fID, "\\|")[[1]]
  checkboxGroupInput("sc2d1sub2", "Select which cells to show", inline = TRUE,
                     choices = sub, selected = sub)
})
observeEvent(input$sc2d1sub1non, {
  sub = strsplit(sc2conf[UI == input$sc2d1sub1]$fID, "\\|")[[1]]
  updateCheckboxGroupInput(session, inputId = "sc2d1sub2", label = "Select which cells to show",
                           choices = sub, selected = NULL, inline = TRUE)
})
observeEvent(input$sc2d1sub1all, {
  sub = strsplit(sc2conf[UI == input$sc2d1sub1]$fID, "\\|")[[1]]
  updateCheckboxGroupInput(session, inputId = "sc2d1sub2", label = "Select which cells to show",
                           choices = sub, selected = sub, inline = TRUE)
})


output$sc2d1oupTxt <- renderUI({
  x <- input$sc2d1grp
  # # updateTextAreaInput(session, "sc1d1inp", value = paste0(sc1def[[x]], collapse = ", "))
  textAreaInput("sc2d1inp", HTML("List of motif names (Max 50 motifs, separated by , or ; or newline):"),
                height = "110px",width = "600px",
                value = paste0(sc2def[[x]], collapse = ", ")
                # value = paste0(sc1def$genes, collapse = ", ")
  ) %>%
    helper(type = "inline", size = "m", fade = TRUE,
           title = "List of motifs to plot on heatmap",
           content = c("Input motifs to plot",
                       "- Maximum 50 motifs (due to ploting space limitations)",
                       "- Motifs should be separated by comma, semicolon or newline"))
#   # geneList = scGeneList(input$sc1d1inp, sc1gene)
#   # if(nrow(geneList) > 50){
#   #   HTML("More than 50 input genes! Please reduce the gene list!")
#   # } else {
#   #   oup = paste0(nrow(geneList[present == TRUE]), "genes OK and will be plotted")
#   #   if(nrow(geneList[present == FALSE]) > 0){
#   #     oup = paste0(oup, "<br/>",
#   #                  nrow(geneList[present == FALSE]), "genes not found (",
#   #                  paste0(geneList[present == FALSE]$gene, collapse = ", "), ")")
#   #   }
#   # HTML(oup)
#   # }
})

output$sc2d1oup <- renderPlot({
  scBubbHeat2(sc2conf, sc2meta, input$sc2d1inp, input$sc2d1grp, input$sc2d1plt,
             input$sc2d1sub1, input$sc2d1sub2, "sc2gexpr.h5", sc2gene,
             input$sc2d1scl, input$sc2d1row, input$sc2d1col,
             input$sc2d1cols, input$sc2d1fsz)
})
output$sc2d1oup.ui <- renderUI({
  plotOutput("sc2d1oup", height = pList3[input$sc2d1psz])
})
output$sc2d1oup.pdf <- downloadHandler(
  filename = function() { paste0("sc2",input$sc2d1plt,"_",input$sc2d1grp,".pdf") },
  content = function(file) {

   
    pdf(file, height = input$sc2d1oup.h, width = input$sc2d1oup.w)
    print(scBubbHeat2(sc2conf, sc2meta, input$sc2d1inp, input$sc2d1grp, input$sc2d1plt,
                      input$sc2d1sub1, input$sc2d1sub2, "sc2gexpr.h5", sc2gene,
                      input$sc2d1scl, input$sc2d1row, input$sc2d1col,
                      input$sc2d1cols, input$sc2d1fsz, save = TRUE))
    dev.off()
    #  
  })
output$sc2d1oup.png <- downloadHandler(
  filename = function() { paste0("sc2",input$sc2d1plt,"_",input$sc2d1grp,".png") },
  content = function(file) {
   
    png(file,  width = 465, height = 225, units='mm', res = 300) 
    print(scBubbHeat2(sc2conf, sc2meta, input$sc2d1inp, input$sc2d1grp, input$sc2d1plt,
                      input$sc2d1sub1, input$sc2d1sub2, "sc2gexpr.h5", sc2gene,
                      input$sc2d1scl, input$sc2d1row, input$sc2d1col,
                      input$sc2d1cols, input$sc2d1fsz, save = TRUE) )
dev.off()
   
  })
   
###########################################################  


      
}) 
 
 
 
 
