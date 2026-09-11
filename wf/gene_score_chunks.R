# Export gene scores without assembling a project-wide sparse matrix.
write_gene_score_chunks <- function(
  proj,
  runs,
  chunk_root,
  feature_chunk_size = 500,
  threads = ArchR::getArchRThreads()
) {
  if (!"asMatrix" %in% names(formals(ArchR::getMatrixFromProject))) {
    stop("Chunked gene export requires the pinned ArchR fork with asMatrix support.")
  }
  dir.create(chunk_root, recursive = TRUE, showWarnings = FALSE)

  run_ids <- vapply(runs, function(run) run[[1]], character(1))
  run_chunk_dirs <- file.path(
    chunk_root,
    sprintf("run_%03d", seq_along(run_ids))
  )
  names(run_chunk_dirs) <- run_ids
  for (run_chunk_dir in run_chunk_dirs) {
    dir.create(run_chunk_dir, recursive = TRUE, showWarnings = FALSE)
  }

  run_cells <- lapply(run_ids, function(run_id) {
    proj$cellNames[as.character(proj$Sample) == run_id]
  })
  names(run_cells) <- run_ids

  empty_runs <- names(run_cells)[lengths(run_cells) == 0]
  if (length(empty_runs) > 0) {
    stop(
      "No ArchR cells found for run(s): ",
      paste(empty_runs, collapse = ", ")
    )
  }

  impute_weights <- ArchR::getImputeWeights(proj)
  if (is.null(impute_weights) || length(impute_weights) == 0) {
    stop("No imputation weights found in the checkpointed ArchR project.")
  }

  matrix_seqnames <- ArchR::getSeqnames(
    ArchRProj = proj,
    useMatrix = "GeneScoreMatrix"
  )
  if (length(matrix_seqnames) == 0) {
    stop("No GeneScoreMatrix seqnames found in the ArchR project.")
  }

  chunk_index <- 0
  for (seq_index in seq_along(matrix_seqnames)) {
    seqname <- matrix_seqnames[[seq_index]]
    message(
      "Reading GeneScoreMatrix seqname ",
      seq_index,
      " of ",
      length(matrix_seqnames),
      ": ",
      seqname
    )

    # The pinned ArchR fork densifies each Arrow matrix when asMatrix is TRUE.
    # Restricting each call to one seqname bounds that allocation instead of
    # materializing every gene across all cells at once.
    gene_matrix <- ArchR::getMatrixFromProject(
      ArchRProj = proj,
      useMatrix = "GeneScoreMatrix",
      useSeqnames = seqname,
      threads = threads,
      asMatrix = TRUE
    )
    gene_assay <- SummarizedExperiment::assay(
      gene_matrix,
      "GeneScoreMatrix"
    )
    feature_names <- as.character(SummarizedExperiment::rowData(gene_matrix)$name)

    if (length(feature_names) != nrow(gene_assay)) {
      stop("Gene feature names do not match matrix rows for seqname ", seqname)
    }

    run_col_indices <- lapply(run_cells, function(cells) {
      match(cells, colnames(gene_assay))
    })
    missing_runs <- names(run_col_indices)[vapply(
      run_col_indices,
      anyNA,
      logical(1)
    )]
    if (length(missing_runs) > 0) {
      stop(
        "GeneScoreMatrix is missing cells for run(s): ",
        paste(missing_runs, collapse = ", ")
      )
    }

    n_feature_chunks <- ceiling(nrow(gene_assay) / feature_chunk_size)
    for (feature_index in seq_len(n_feature_chunks)) {
      chunk_index <- chunk_index + 1
      start_idx <- (feature_index - 1) * feature_chunk_size + 1
      end_idx <- min(feature_index * feature_chunk_size, nrow(gene_assay))
      feature_idx <- start_idx:end_idx

      message(
        "Imputing feature chunk ",
        feature_index,
        " of ",
        n_feature_chunks,
        " for ",
        seqname,
        " (global chunk ",
        chunk_index,
        ")"
      )

      mat_chunk <- gene_assay[feature_idx, , drop = FALSE]
      imputed_chunk <- ArchR::imputeMatrix(
        mat = mat_chunk,
        imputeWeights = impute_weights,
        threads = threads
      )

      # The pinned ArchR imputeMatrix implementation subsets its final result
      # without drop = FALSE. A one-feature chunk is therefore returned as a
      # dimensionless vector. Restore the feature-by-cell shape explicitly and
      # validate all other chunks before persisting them.
      expected_dim <- dim(mat_chunk)
      if (is.null(dim(imputed_chunk))) {
        if (length(imputed_chunk) != prod(expected_dim)) {
          stop(
            "Dimensionless imputed chunk has ",
            length(imputed_chunk),
            " values; expected ",
            prod(expected_dim)
          )
        }
        imputed_chunk <- matrix(
          imputed_chunk,
          nrow = expected_dim[[1]],
          ncol = expected_dim[[2]]
        )
      }
      if (!identical(as.integer(dim(imputed_chunk)), as.integer(expected_dim))) {
        stop(
          "Imputed chunk dimensions ",
          paste(dim(imputed_chunk), collapse = " x "),
          " do not match expected dimensions ",
          paste(expected_dim, collapse = " x ")
        )
      }
      rownames(imputed_chunk) <- feature_names[feature_idx]
      colnames(imputed_chunk) <- colnames(mat_chunk)

      for (run_id in run_ids) {
        # Imputed gene scores are effectively dense. Keep each bounded,
        # per-sample chunk dense instead of paying the larger dgCMatrix
        # overhead or introducing its nonzero-element limit.
        run_chunk <- imputed_chunk[
          , run_col_indices[[run_id]], drop = FALSE
        ]
        if (!is.matrix(run_chunk) || inherits(run_chunk, "sparseMatrix")) {
          run_chunk <- as.matrix(run_chunk)
        }
        chunk_path <- file.path(
          run_chunk_dirs[[run_id]],
          sprintf("chunk_%05d.rds", chunk_index)
        )
        saveRDS(run_chunk, file = chunk_path, compress = FALSE)
        rm(run_chunk)
      }

      rm(mat_chunk, imputed_chunk)
      gc(verbose = FALSE, full = TRUE)
    }

    rm(gene_matrix, gene_assay, feature_names, run_col_indices)
    gc(verbose = FALSE, full = TRUE)
  }

  rm(impute_weights)
  gc(verbose = FALSE, full = TRUE)
  run_chunk_dirs
}


# Only retain names of empty genes, never the full project-wide matrix.
find_empty_gene_features <- function(proj) {
  seqnames <- ArchR::getSeqnames(proj, useMatrix = "GeneScoreMatrix")
  if (length(seqnames) == 0) stop("No GeneScoreMatrix seqnames found.")
  empty_genes <- character(0)
  for (seqname in seqnames) {
    message("Scanning empty gene features on ", seqname)
    gene_matrix <- ArchR::getMatrixFromProject(
      ArchRProj = proj, useMatrix = "GeneScoreMatrix",
      useSeqnames = seqname, threads = 1, asMatrix = TRUE
    )
    gene_assay <- SummarizedExperiment::assay(gene_matrix, "GeneScoreMatrix")
    gene_names <- as.character(SummarizedExperiment::rowData(gene_matrix)$name)
    if (length(gene_names) != nrow(gene_assay)) {
      stop("Gene names do not match matrix rows on ", seqname)
    }
    empty_genes <- c(empty_genes, gene_names[which(Matrix::rowSums(gene_assay) == 0)])
    rm(gene_matrix, gene_assay, gene_names)
    gc(verbose = FALSE)
  }
  empty_genes
}
