library("data.table")

scGeneList <- function(inp, inpGene) {
  geneList <- data.table::data.table(
    gene = unique(trimws(strsplit(inp, ",|;| ")[[1]])),
    present = TRUE
  )
  geneList[!gene %in% names(inpGene)]$present <- FALSE

  print(geneList)

  return(geneList)
}