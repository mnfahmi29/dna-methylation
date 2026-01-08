# R/helpers/batch_tools.R
# ------------------------------------------------------------
# Batch correction helpers :))
#
# Philosophy:
# - Best batch correction is done on INTENSITIES (meth/unmeth)
# - Beta-level correction exists, but must be used consciously
#
# This file provides:
# 1) intensity-level correction (recommended)
# 2) beta-level correction (fallback, clearly warned)
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(limma)
})

# -----------------------------
# 1) Intensity-level batch correction (preferred)
# -----------------------------
batch_correct_meth_unmeth <- function(
    methy,
    unmethy,
    material,
    ffpe_label   = "FFPE",
    frozen_label = "Frozen"
) {
  material <- trimws(as.character(material))
  
  has_ffpe   <- any(material == ffpe_label,   na.rm = TRUE)
  has_frozen <- any(material == frozen_label, na.rm = TRUE)
  
  if (!(has_ffpe && has_frozen)) {
    message("Batch correction skipped (need BOTH FFPE and Frozen samples).")
    return(list(
      methy.ba   = methy,
      unmethy.ba = unmethy,
      batch      = NULL
    ))
  }
  
  # simple 2-level batch coding
  batch <- ifelse(material == ffpe_label, 2, 1)
  
  message("Running intensity-level batch correction (meth/unmeth) ...")
  
  methy.ba   <- 2 ^ removeBatchEffect(log2(methy + 1), batch = batch)
  unmethy.ba <- 2 ^ removeBatchEffect(log2(unmethy + 1), batch = batch)
  
  list(
    methy.ba   = methy.ba,
    unmethy.ba = unmethy.ba,
    batch      = batch
  )
}

# -----------------------------
# 2) Beta-level batch correction (use with caution)
# -----------------------------
batch_correct_betas <- function(
    betas,
    batch,
    ref_level = NULL
) {
  stopifnot(is.matrix(betas) || is.data.frame(betas))
  betas <- as.matrix(betas)
  
  batch <- as.factor(trimws(as.character(batch)))
  
  if (nlevels(batch) <= 1) {
    message("Beta batch correction skipped (only one batch present).")
    return(betas)
  }
  
  if (!is.null(ref_level) && ref_level %in% levels(batch)) {
    batch <- relevel(batch, ref = ref_level)
  }
  
  message("⚠️  Running BETA-level batch correction (interpret carefully) ._.")
  
  # Work on M-values (standard approach)
  mvals <- log2(betas / (1 - betas))
  mvals_cor <- removeBatchEffect(mvals, batch = batch)
  
  betas_cor <- 2^mvals_cor / (1 + 2^mvals_cor)
  
  # clamp numerical noise
  betas_cor[betas_cor < 0] <- 0
  betas_cor[betas_cor > 1] <- 1
  
  betas_cor
}
