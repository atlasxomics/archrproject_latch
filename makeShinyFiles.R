#' Generate data files required for shiny app
#'
#' Generate data files required for shiny app. Five files will be generated,
#' namely (i) the shinycell config \code{prefix_conf.rds}, (ii) the gene
#' mapping object config \code{prefix_gene.rds}, (iii) the single-cell gene
#' expression \code{prefix_gexpr.h5}, (iv) the single-cell metadata
#' \code{prefix_meta.rds} and (v) the defaults for the Shiny app
#' \code{prefix_def.rds}. A prefix is specified in each file to allow for the
#' use of multiple single-cell datasets in one Shiny app. Note that both
#' \code{makeShinyFiles} and \code{makeShinyCodes} functions are ran when
#' running the wrapper function \code{makeShinyApp}.
#'
#' @param obj input single-cell object for Seurat (v3+) / SingleCellExperiment
#'   data or input file path for h5ad / loom files
#' @param scConf shinycell config data.table
#' @param gex.assay assay in single-cell data object to use for plotting
#'   gene expression, which must match one of the following:
#'   \itemize{
#'     \item{Seurat objects}: "RNA" or "integrated" assay,
#'       default is "RNA"
#'     \item{SCE objects}: "logcounts" or "normcounts" or "counts",
#'       default is "logcounts"
#'     \item{h5ad files}: "X" or any assay in "layers",
#'       default is "X"
#'     \item{loom files}: "matrix" or any assay in "layers",
#'       default is "matrix"
#'   }
#' @param gex.slot slot in single-cell assay to plot. This is only used
#'   for Seurat objects (v3+). Default is to use the "data" slot
#' @param gene.mapping specifies whether to convert human / mouse Ensembl gene
#'   IDs (e.g. ENSG000xxx / ENSMUSG000xxx) into "user-friendly" gene symbols.
#'   Set this to \code{TRUE} if you are using Ensembl gene IDs. Default is
#'   \code{FALSE} which is not to perform any conversion. Alternatively, users
#'   can supply a named vector where \code{names(gene.mapping)} correspond
#'   to the actual gene identifiers in the gene expression matrix and
#'   \code{gene.mapping} correspond to new identifiers to map to
#' @param shiny.prefix specify file prefix
#' @param shiny.dir specify directory to create the shiny app in
#' @param default.gene1 specify primary default gene to show
#' @param default.gene2 specify secondary default gene to show
#' @param default.multigene character vector specifying default genes to
#'   show in bubbleplot / heatmap
#' @param default.dimred character vector specifying the two default dimension
#'   reductions. Default is to use UMAP if not TSNE embeddings
#' @param chunkSize number of genes written to h5file at any one time. Lower
#'   this number to reduce memory consumption. Should not be less than 10
#'
#' @return data files required for shiny app
#'
#' @author John F. Ouyang
#'
#' @import data.table hdf5r reticulate
#'
#' @examples
#' makeShinyFiles(
#'   seu,
#'   scConf,
#'   gex.assay = "RNA",
#'   gex.slot = "data",
#'   shiny.prefix = "sc1",
#'   shiny.dir = "shinyApp/",
#'   default.gene1 = "GATA3",
#'   default.gene2 = "DNMT3L",
#'   default.multigene = c(
#'     "ANPEP","NANOG","ZIC2","NLGN4X","DNMT3L",
#'     "DPPA5","SLC7A2","GATA3","KRT19"
#'    ),
#'   default.dimred = c("UMAP_1", "UMAP_2")
#' )
#'
#' @export

makeShinyFiles <- function(
  obj,
  scConf,
  gex.assay = NA,
  gex.slot = c("data", "scale.data", "counts"),
  gene.mapping = FALSE,
  shiny.prefix = "sc1",
  shiny.dir = "shinyApp/",
  default.gene1 = NA,
  default.gene2 = NA,
  default.multigene = NA,
  default.dimred = NA,
  chunkSize = 500
) {

  ### Preprocessing and checks
  # Generate defaults for gex.assay / gex.slot

  gex.matdim <- dim(obj)
  gex.rownm <- rownames(obj)
  gex.colnm <- colnames(obj)
  defGenes <- Seurat::VariableFeatures(obj)[1:10]
  if (is.na(defGenes[1])) {
    warning(paste0(
      "Variable genes for seurat object not found! Have you ",
      "ran `FindVariableFeatures` or `SCTransform`?"
    ))
    defGenes <- gex.rownm[1:10]
  }
  sc1meta <- data.table(sampleID = rownames(obj@meta.data), obj@meta.data)

  gene.mapping <- defGenes

  # Check default.gene1 / default.gene2 / default.multigene
  default.gene1 <- default.gene1[1]
  default.gene2 <- default.gene2[1]

  if (default.gene1 %in% gene.mapping) {
    default.gene1 <- default.gene1
  } else {
    warning("default.gene1 doesn't exist in gene expression, using defaults...")
    default.gene1 <- defGenes[1]
  }

  if (default.gene2 %in% gene.mapping) {
    default.gene2 <- default.gene2
  } else {
    warning("default.gene2 doesn't exist in gene expression, using defaults...")
    default.gene2 <- defGenes[2]
  }

  if (all(default.multigene %in% gene.mapping)) {
    default.multigene <- default.multigene
  } else {
    warning(paste0(
      "default.multigene doesn't exist in gene expression using defaults..."
    ))
    default.multigene <- defGenes[1:9]
  }

  ### Actual object generation
  # Make XXXmeta.rds and XXXconf.rds (updated with dimred info)
  sc1conf <- scConf
  sc1conf$dimred <- FALSE
  sc1meta <- sc1meta[, c("sampleID", as.character(sc1conf$ID)), with = FALSE]

  # Factor metadata again
  for (i in as.character(sc1conf[!is.na(fID)]$ID)) {
    sc1meta[[i]] <- factor(
      sc1meta[[i]], levels = strsplit(sc1conf[ID == i]$fID, "\\|")[[1]]
    )
    levels(sc1meta[[i]]) <- strsplit(sc1conf[ID == i]$fUI, "\\|")[[1]]
    sc1conf[ID == i]$fID <- sc1conf[ID == i]$fUI
  }

  # Extract dimred and append to both XXXmeta.rds and XXXconf.rds...

  for (iDR in names(obj@reductions)) {
    drMat <- obj@reductions[[iDR]]@cell.embeddings
    colnames(drMat) <- c(paste0(iDR, "_1"), paste0(iDR, "_2"))
    if (ncol(drMat) > 5) {
      drMat <- drMat[, 1:5]  # Take first 5 components only
    }
    drMat <- drMat[sc1meta$sampleID, ]  # Ensure ordering
    drMat <- as.data.table(drMat)
    sc1meta <- cbind(sc1meta, drMat)

    # Update sc1conf accordingly
    tmp <- data.table(
      ID = colnames(drMat),
      UI = colnames(drMat),
      fID = NA,
      fUI = NA,
      fCL = NA,
      fRow = NA,
      default = 0,
      grp = FALSE,
      dimred = TRUE
    )
    sc1conf <- rbindlist(list(sc1conf, tmp))
  }
  sc1conf$ID <- as.character(sc1conf$ID)  # Remove levels

  # Make XXXgexpr.h5
  if (!dir.exists(shiny.dir)) {
    dir.create(shiny.dir)
  }
  filename <- paste0(shiny.dir, "/", shiny.prefix, "gexpr.h5")
  sc1gexpr <- H5File$new(filename, mode = "w")
  sc1gexpr.grp <- sc1gexpr$create_group("grp")
  sc1gexpr.grp.data <- sc1gexpr.grp$create_dataset(
    "data",
    dtype = h5types$H5T_NATIVE_FLOAT,
    space = H5S$new("simple", dims = gex.matdim, maxdims = gex.matdim),
    chunk_dims = c(1, gex.matdim[2])
  )
  chk <- chunkSize
  while (chk > (gex.matdim[1] - 8)) {
    chk <- floor(chk / 2)  # Account for cases where nGene < chunkSize
  }

  for (i in 1:floor((gex.matdim[1] - 8) / chk)) {
    sc1gexpr.grp.data[((i - 1) * chk + 1) : (i * chk), ] <- as.matrix(
      obj[[gex.assay[1]]]$data[((i - 1) * chk + 1) : (i * chk), ]
    )
  }
  sc1gexpr.grp.data[(i * chk + 1) : gex.matdim[1], ] <- as.matrix(
    obj[[gex.assay[1]]]$data[(i * chk + 1) : gex.matdim[1], ]
  )

  sc1gexpr$close_all()


  # Make XXXgenes.rds
  sc1gene <- seq(gex.matdim[1])
  names(gene.mapping) <- NULL
  names(sc1gene) <- gene.mapping
  sc1gene <- sc1gene[order(names(sc1gene))]
  sc1gene <- sc1gene[order(nchar(names(sc1gene)))]

  # Make XXXdef.rds (list of defaults)
  if (all(default.dimred %in% sc1conf[dimred == TRUE]$ID)) {
    default.dimred[1] <- sc1conf[ID == default.dimred[1]]$UI
    default.dimred[2] <- sc1conf[ID == default.dimred[2]]$UI
  } else if (all(default.dimred %in% sc1conf[dimred == TRUE]$UI)) {
    default.dimred <- default.dimred   # Nothing happens
  } else {
    warn <- TRUE
    if (is.na(default.dimred[1])) {
      default.dimred <- "umap"
      warn <- FALSE
    }

    # Try to guess... and give a warning
    guess <- gsub("[0-9]", "", default.dimred[1])
    if (
      length(grep(guess, sc1conf[dimred == TRUE]$UI, ignore.case = TRUE)) >= 2
    ) {
      default.dimred <- sc1conf[dimred == TRUE]$UI[
        grep(guess, sc1conf[dimred == TRUE]$UI, ignore.case = TRUE)[1:2]
      ]
    } else {
      nDR <- length(sc1conf[dimred == TRUE]$UI)
      default.dimred <- sc1conf[dimred == TRUE]$UI[(nDR - 1):nDR]
    }
    if (warn) {
      warning(paste0(
        "default.dimred not found, switching to ",
        default.dimred[1],
        " and ",
        default.dimred[1])
      )
    }  # Warn if user-supplied default.dimred is not found
  }

  # Note that we stored the display name here
  sc1def <- list()
  sc1def$meta1 <- sc1conf[default == 1]$UI   # Use display name
  sc1def$meta2 <- sc1conf[default == 2]$UI   # Use display name
  sc1def$gene1 <- default.gene1              # Actual == Display name
  sc1def$gene2 <- default.gene2              # Actual == Display name
  sc1def$genes <- default.multigene          # Actual == Display name
  sc1def$dimred <- default.dimred            # Use display name
  tmp <- nrow(sc1conf[default != 0 & grp == TRUE])
  if (tmp == 2) {
    sc1def$grp1 <- sc1def$meta1
    sc1def$grp2 <- sc1def$meta2
  } else if (tmp == 1) {
    sc1def$grp1 <- sc1conf[default != 0 & grp == TRUE]$UI
    if (nrow(sc1conf[default == 0 & grp == TRUE]) == 0) {
      sc1def$grp2 <- sc1def$grp1
    } else {
      sc1def$grp2 <- sc1conf[default == 0 & grp == TRUE]$UI[1]
    }
  } else {
    sc1def$grp1 <- sc1conf[default == 0 & grp == TRUE]$UI[1]
    if (nrow(sc1conf[default == 0 & grp == TRUE]) < 2) {
      sc1def$grp2 <- sc1def$grp1
    } else {
      sc1def$grp2 <- sc1conf[default == 0 & grp == TRUE]$UI[2]
    }
  }
  sc1conf <- sc1conf[, -c("fUI", "default"), with = FALSE]

  ### Saving objects
  saveRDS(sc1conf, file = paste0(shiny.dir, "/", shiny.prefix, "conf.rds"))
  saveRDS(sc1meta, file = paste0(shiny.dir, "/", shiny.prefix, "meta.rds"))
  saveRDS(sc1gene, file = paste0(shiny.dir, "/", shiny.prefix, "gene.rds"))
  saveRDS(sc1def,  file = paste0(shiny.dir, "/", shiny.prefix, "def.rds"))
  return(sc1conf)
}
