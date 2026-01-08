# R/helpers/probe_id_tools.R
# ------------------------------------------------------------
# 🧬 Probe ID sanity tools (EPIC / EPICv2 friendly) :))
#
# Why this exists:
# - EPICv2 sometimes appends suffixes to CpG IDs
#     cg00000029_TC21  →  cg00000029
# - Older references / CNV tools often expect the trimmed version
#
# If you don’t standardize probe IDs:
# - intersections become empty
# - CNV reference breaks silently
# - you lose hours debugging ._.
#
# Design principle:
# Make probe alignment boring, explicit, and deterministic.
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(minfi)
})

# -----------------------------
# Small IO helper
# -----------------------------
ensure_dir <- function(path) {
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  invisible(path)
}

# ============================================================
# 1A) Trim CpG probe IDs (EPICv2 → cg########)
# ============================================================
# ✅ Trims only IDs that START with "cg" digits
# ✅ Leaves non-cg probes untouched (rs / ch / etc.) to avoid breaking them
#
# Examples:
#   "cg00000029_TC21" -> "cg00000029"
#   "rs12345"         -> "rs12345" (unchanged)

trim_cg_ids <- function(ids, verbose = TRUE) {
  ids <- as.character(ids)
  
  out <- sub("^(cg\\d+).*", "\\1", ids)
  
  # If an id doesn't start with cg####, keep it unchanged (safe default)
  keep_same <- !grepl("^cg\\d+", ids)
  out[keep_same] <- ids[keep_same]
  
  if (verbose) {
    n_changed <- sum(out != ids, na.rm = TRUE)
    message("✂️ trim_cg_ids(): changed ", n_changed, " / ", length(ids), " IDs")
  }
  
  out
}

# ============================================================
# 1B) Standardize cg IDs (alias) + generic duplicate dropper
# ============================================================

# Alias
standardize_cg_ids <- function(ids, verbose = TRUE) {
  trim_cg_ids(ids, verbose = verbose)
}

# Generic helper: drop duplicated probes AFTER any standardization step
# - For vectors: uses the vector itself as the ID list
# - For matrices/data.frames: choose rownames or colnames
drop_duplicate_probes_keep_first <- function(x, use = c("rownames", "colnames"), verbose = TRUE) {
  use <- match.arg(use)
  
  if (is.atomic(x) && is.null(dim(x))) {
    ids <- as.character(x)
    keep <- !duplicated(ids)
    if (verbose && any(!keep)) message("⚠️ Dropping ", sum(!keep), " duplicated IDs (keep-first).")
    return(ids[keep])
  }
  
  ids <- if (use == "rownames") rownames(x) else colnames(x)
  if (is.null(ids)) stop("drop_duplicate_probes_keep_first(): ", use, " are NULL.")
  
  keep <- !duplicated(ids)
  if (verbose && any(!keep)) message("⚠️ Dropping ", sum(!keep), " duplicated ", use, " (keep-first).")
  
  if (use == "rownames") x[keep, , drop = FALSE] else x[, keep, drop = FALSE]
}

# ============================================================
# 2) Trim rownames of a matrix/data.frame (probes in rows)
# ============================================================
# Useful for CNV-style objects where probes are in rows.
# Optionally drops duplicates created by trimming (keep-first policy).
#
# Requirements:
# - rownames(mat) must exist

trim_matrix_rownames_cg <- function(mat, verbose = TRUE, drop_duplicates = TRUE) {
  if (is.null(rownames(mat))) stop("trim_matrix_rownames_cg(): mat has NULL rownames.")
  
  old <- rownames(mat)
  new <- trim_cg_ids(old, verbose = FALSE)
  
  if (drop_duplicates) {
    dup <- duplicated(new)
    if (any(dup)) {
      if (verbose) message("⚠️ Dropping ", sum(dup), " duplicated probes after trim.")
      mat <- mat[!dup, , drop = FALSE]
      new <- new[!dup]
    }
  }
  
  rownames(mat) <- new
  if (verbose) message("✅ trim_matrix_rownames_cg(): done (rows=", nrow(mat), ")")
  
  mat
}

# ============================================================
# 3) Trim probes inside minfi methylation objects
# ============================================================
# Works for typical minfi objects that store probes in rownames:
# - MethylSet
# - GenomicMethylSet
# - RatioSet
# - GenomicRatioSet (often shows up)
#
# Behavior:
# - trims cg IDs
# - can drop duplicates created by trimming

trim_methylation_object_probes <- function(mset, verbose = TRUE, drop_duplicates = TRUE) {
  ok <- is(mset, "MethylSet") ||
    is(mset, "GenomicMethylSet") ||
    is(mset, "RatioSet") ||
    is(mset, "GenomicRatioSet")
  
  if (!ok) {
    stop("trim_methylation_object_probes(): unsupported object class: ", class(mset)[1])
  }
  
  old <- rownames(mset)
  new <- trim_cg_ids(old, verbose = FALSE)
  
  if (drop_duplicates) {
    dup <- duplicated(new)
    if (any(dup)) {
      if (verbose) message("⚠️ Dropping ", sum(dup), " duplicated probes after trim.")
      mset <- mset[!dup, ]
      new <- new[!dup]
    }
  }
  
  rownames(mset) <- new
  if (verbose) message("✅ trim_methylation_object_probes(): done (probes=", nrow(mset), ")")
  
  mset
}

# ============================================================
# 4) Assert shared probe universe across inputs (early scream)
# ============================================================
# Purpose:
# - Prevent silent mismatches (e.g., EPICv2 not trimmed)
# - Return common probe IDs (invisible) for immediate use
#
# Usage patterns:
#   common <- assert_same_probe_universe(betas1, betas2)  # expects colnames overlap
#   common <- assert_same_probe_universe(matA, matB, use = "rownames") # expects rownames overlap
#
# NOTE:
# - By default uses colnames (common in betas matrices: samples x probes)

assert_same_probe_universe <- function(..., use = c("colnames", "rownames")) {
  use <- match.arg(use)
  mats <- list(...)
  
  if (length(mats) < 2) stop("assert_same_probe_universe(): need at least 2 inputs.")
  
  get_ids <- switch(
    use,
    colnames = function(x) colnames(x),
    rownames = function(x) rownames(x)
  )
  
  ids_list <- lapply(mats, get_ids)
  
  if (any(vapply(ids_list, is.null, logical(1)))) {
    stop("assert_same_probe_universe(): some inputs have NULL ", use, ".")
  }
  if (any(vapply(ids_list, length, integer(1)) == 0)) {
    stop("assert_same_probe_universe(): some inputs have 0-length ", use, ".")
  }
  
  common <- Reduce(intersect, ids_list)
  
  if (length(common) == 0) {
    stop(
      "🚨 Probe overlap is EMPTY using ", use, ".\n",
      "Most common cause: ID mismatch (did you trim EPICv2 cg IDs?)."
    )
  }
  
  message("✅ Common probes (", use, "): ", length(common))
  invisible(common)
}
