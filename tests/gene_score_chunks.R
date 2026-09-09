# Dependency-free regression test of orchestration; ArchR calls are mocked.
# Run from the repository root: Rscript --vanilla tests/gene_score_chunks.R
sandbox <- new.env(parent = baseenv())
strip_namespace <- function(x) {
  if (is.call(x)) {
    if (identical(x[[1]], as.name("::"))) return(x[[3]])
    return(as.call(lapply(as.list(x), strip_namespace)))
  }
  x
}
for (expr in parse("wf/gene_score_chunks.R")) {
  eval(strip_namespace(expr), sandbox)
}
cells <- c("sample_extra#c-1", "sample#a-1", "sample#b-1")
raw <- matrix(seq_len(503 * 3) / 100, nrow = 503,
              dimnames = list(paste0("g", seq_len(503)), cells))
weights <- matrix(c(.5, .3, .2, .2, .6, .2, .1, .2, .7), 3)
proj <- list(cellNames = cells, Sample = c("sample_extra", "sample", "sample"))
sandbox$getArchRThreads <- function() 1
sandbox$getImputeWeights <- function(proj) weights
sandbox$getSeqnames <- function(...) c("chr1", "chr2")
reads <- character()
sandbox$getMatrixFromProject <- function(ArchRProj, useMatrix, useSeqnames,
                                        threads, asMatrix) {
  stopifnot(isTRUE(asMatrix), length(useSeqnames) == 1)
  reads <<- c(reads, useSeqnames)
  idx <- if (useSeqnames == "chr1") 1:501 else 502:503
  list(mat = raw[idx, , drop = FALSE])
}
sandbox$assay <- function(x, ...) x$mat
sandbox$rowData <- function(x) list(name = rownames(x$mat))
blocks <- integer()
sandbox$imputeMatrix <- function(mat, imputeWeights, threads) {
  stopifnot(identical(imputeWeights, weights), ncol(mat) == 3)
  blocks <<- c(blocks, nrow(mat))
  result <- mat %*% imputeWeights
  if (nrow(mat) == 1) as.vector(result) else result
}
root <- tempfile()
tryCatch({
  dirs <- sandbox$write_gene_score_chunks(
    proj, list(c("sample"), c("sample_extra")), root
  )
  expected <- raw %*% weights
  dimnames(expected) <- dimnames(raw)
  for (run in names(dirs)) {
    parts <- lapply(sort(list.files(dirs[[run]], full.names = TRUE)), readRDS)
    actual <- do.call(rbind, parts)
    stopifnot(is.matrix(actual))
    stopifnot(isTRUE(all.equal(actual,
      expected[, proj$Sample == run, drop = FALSE])))
  }
  stopifnot(identical(reads, c("chr1", "chr2")))
  stopifnot(identical(blocks, c(500L, 1L, 2L)))
  cat("PASS: chunked scores equal full multiplication; cell/gene ordering,\n",
      "sample prefixes, global weights, and single-gene chunks preserved.\n")
}, finally = unlink(root, recursive = TRUE))
