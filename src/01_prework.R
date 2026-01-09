# src/01_prework.R
# ============================================================
# 🧬 DNA Methylation Pipeline — From Raw IDATs to Clean Betas
# ============================================================
#
# What does this script do? 🤔
# --------------------------------
# It builds a single, harmonized beta matrix (samples × probes)
# by combining:
#
#   1️⃣ Public reference data:
#      - GSE90496 raw IDAT files (Illumina 450k / EPICv1)
#      - Read in chunks so your RAM doesn’t cry 🥲
#
#   2️⃣ Your own cases:
#      - EPICv2 IDATs (here: GBM + ICM, just as examples)
#      - EPICv2 probe IDs are trimmed to stay compatible
#      - You could also using another typer below EPICv2, 
#        you just need to ignore the probe trimming section
#
# Final outputs 🎁
# ----------------
# ✔ results/betas_all.RData
#     - betas_all     : numeric matrix (samples × probes)
#     - anno_combined : sample annotations aligned to betas_all
#
# ✔ (optional) results/pca_rspectra.rds
#     - fast PCA sanity check (RSpectra-based)
#
# ------------------------------------------------------------
# 🚦 Reproducibility rules (VERY IMPORTANT)
# ------------------------------------------------------------
# ❌ NO install.packages() here
# ❌ NO BiocManager::install() here
#
# ✅ The R environment is managed by `renv`
#    → Run `renv::restore()` ONCE before running this script
#
# ✅ This script assumes:
#    - packages already exist
#    - helpers are versioned inside this repo
#
# Think of it like a recipe 🍳:
#   Ingredients (packages) are prepared beforehand.
#   This script only cooks.
# ------------------------------------------------------------

rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(minfi)     # DNA methylation backbone 🧬
  library(GEOquery)  # Download GEO metadata 📦
  library(limma)     # Batch correction & linear algebra 🧮
})

# ------------------------------------------------------------
# 0) Load helpers (our “toolbox”) 🧰
# ------------------------------------------------------------
# These helpers keep the main script readable and boring
# (boring = reproducible = good science)._.

# Probe handling, EPICv2 trimming, beta extraction
source(file.path("helpers", "probe_id_tools.R"))

# Batch correction helpers (explicit, opt-in only)
source(file.path("helpers", "batch_tools.R"))

# ------------------------------------------------------------
# Vendor helpers (external, but pinned in this repo) 📦
# ------------------------------------------------------------
# MNPprocessIDAT_functions.R
#   - patched preprocessing for EPICv2 compatibility
#   - overrides preprocessIllumina() when available
#
# RSpectra_pca.R
#   - fast PCA for sanity checking large matrices
if (file.exists(file.path("helpers", "MNPTraining", "MNPprocessIDAT_functions.R"))) {
  source(file.path("helpers", "MNPTraining", "MNPprocessIDAT_functions.R"))
}
if (file.exists(file.path("helpers", "MNPTraining", "RSpectra_pca.R"))) {
  source(file.path("helpers", "MNPTraining", "RSpectra_pca.R"))
}

ensure_dir("results")  # create output folder if needed

# ------------------------------------------------------------
# 1) Paths — tell the script where the data lives 🧭
# ------------------------------------------------------------
GSE_RAW_DIR  <- "data/GSE90496_RAW"   
# you need first download that 22GB dataset then unpacked it, 
# the GSE is set to GSE90496 due to this example is for analyze Brain Tumor Idat Files

CASE_GBM_DIR <- "data/Case GBM_MES"
CASE_ICM_DIR <- "data/Intracranial Mesenchymal FET_CREB Fusion"
# GBM and ICM are example for this code, 
# you could use as many idat files as you like

# Check if your file really in the path
stopifnot(
  dir.exists(GSE_RAW_DIR),
  dir.exists(CASE_GBM_DIR),
  dir.exists(CASE_ICM_DIR)
)

# ------------------------------------------------------------
# 2) Preprocessing wrapper 🧠
# ------------------------------------------------------------
# If the MNP helper exists → use it
# Otherwise → fall back to vanilla minfi
preprocess_any <- function(rgset) {
  if (exists("MNPpreprocessIllumina", mode = "function")) {
    MNPpreprocessIllumina(rgset)
  } else {
    preprocessIllumina(rgset)
  }
}

# ------------------------------------------------------------
# 3) Download GEO metadata (GSE90496) 📥
# ------------------------------------------------------------
message("📥 Fetching GEO metadata: GSE90496 ...")

gse  <- getGEO("GSE90496", GSEMatrix = TRUE, getGPL = FALSE)
anno <- pData(gse$GSE90496_series_matrix.txt.gz)

# Keep all reference classes for now
anno_ref <- anno

# Convert GEO supplementary file URLs → IDAT basenames
fname <- gsub("_Grn.*", "", gsub(".*suppl/", "", anno_ref$supplementary_file))
filepath_gse_all <- file.path(GSE_RAW_DIR, fname)

# ------------------------------------------------------------
# 4) Read EPICv2 cases 🧪
# ------------------------------------------------------------
message("🧫 Reading EPICv2 cases ...")

RG_gbm <- read.metharray.exp(CASE_GBM_DIR, verbose = TRUE)
RG_icm <- read.metharray.exp(CASE_ICM_DIR, verbose = TRUE)

colnames(RG_gbm) <- "Case_GBM"
colnames(RG_icm) <- "Case_ICM"

anno_cases <- data.frame(
  row.names = c("Case_GBM", "Case_ICM"),
  `methylation class:ch1` = c("Case_GBM", "Case_ICM"),
  `material:ch1` = c("Frozen", "Frozen"),  # change if FFPE
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 5) Chunked read + preprocess of GSE data 🧱
# ------------------------------------------------------------
# Why chunking?
# Because loading hundreds of IDATs at once would lead to computer's 'poof'💀
# Well, I've tried it, if you feel curious just skip this chunking section
# Anywy, if you have a computer spec that compatible, you could also skip this section

chunk_size <- 25
chunks <- split(filepath_gse_all,
                ceiling(seq_along(filepath_gse_all) / chunk_size))

betas_chunks <- vector("list", length(chunks))
sample_names_chunks <- vector("list", length(chunks))

message("🧩 Reading + preprocessing GSE in chunks ...")
for (i in seq_along(chunks)) {
  message(sprintf("  - Chunk %d / %d", i, length(chunks)))
  
  RG_chunk <- read.metharray(chunks[[i]], verbose = FALSE)
  M_chunk  <- preprocess_any(RG_chunk)
  
  # Always standardize orientation: samples × probes
  beta_chunk <- betas_from_mset(M_chunk)
  
  betas_chunks[[i]] <- beta_chunk
  sample_names_chunks[[i]] <- rownames(beta_chunk)
  
  rm(RG_chunk, M_chunk, beta_chunk)
  gc()
}

betas_gse <- do.call(rbind, betas_chunks)
rownames(betas_gse) <- unlist(sample_names_chunks)

# ------------------------------------------------------------
# 6) Preprocess cases + EPICv2 probe trimming ✂️ (for EPICv2)
# ------------------------------------------------------------
# EPICv2 adds suffixes to CpG IDs.
# We trim them back to cg######## so references still match.

M_gbm <- preprocess_any(RG_gbm)
M_icm <- preprocess_any(RG_icm)

M_gbm <- trim_methylation_object_probes(M_gbm, verbose = TRUE)
M_icm <- trim_methylation_object_probes(M_icm, verbose = TRUE)

M_gbm <- filter_probes_mset(M_gbm)
M_icm <- filter_probes_mset(M_icm)

beta_gbm <- betas_from_mset(M_gbm)
beta_icm <- betas_from_mset(M_icm)

# ------------------------------------------------------------
# 7) Harmonize probe universe 🤝
# ------------------------------------------------------------
# Only probes shared by ALL datasets survive.
# This avoids silent NA explosions later.

common_probes <- intersect_probe_universe(list(
  gse = betas_gse,
  gbm = beta_gbm,
  icm = beta_icm
), verbose = TRUE)

betas_all <- rbind(
  betas_gse[, common_probes, drop = FALSE],
  beta_gbm[,  common_probes, drop = FALSE],
  beta_icm[,  common_probes, drop = FALSE]
)

# ------------------------------------------------------------
# 8) Build aligned annotation 🏷️
# ------------------------------------------------------------
anno_ref2 <- data.frame(
  row.names = rownames(betas_gse),
  `methylation class:ch1` = anno_ref$`methylation class:ch1`,
  `material:ch1`         = anno_ref$`material:ch1`,
  stringsAsFactors = FALSE
)

anno_combined <- rbind(anno_ref2, anno_cases)
anno_combined <- anno_combined[rownames(betas_all), , drop = FALSE]
stopifnot(nrow(anno_combined) == nrow(betas_all))

# ------------------------------------------------------------
# 9) Batch correction policy 🧼
# ------------------------------------------------------------
message("🧼 Batch correction is OFF by default.")
message("    Transparency > magic ✨")

# Uncomment ONLY if you know why you need it:
# betas_all <- batch_correct_betas(
#   betas_all,
#   batch     = anno_combined$`material:ch1`,
#   ref_level = "Frozen"
# )

# ------------------------------------------------------------
# 10) Optional PCA sanity check 🧠
# ------------------------------------------------------------
if (exists("run_pca_rspectra", mode = "function")) {
  message("🧠 Running optional RSpectra PCA sanity check ...")
  
  pca_out <- run_pca_rspectra(
    betas_all,
    n_components = 20,
    center = TRUE,
    scale = FALSE
  )
  
  saveRDS(
    list(pca = pca_out, anno = anno_combined),
    file = file.path("results", "pca_rspectra.rds")
  )
  
  message("Saved: results/pca_rspectra.rds ✅")
} else {
  message("🧠 PCA skipped (RSpectra helper not found).")
}

# ------------------------------------------------------------
# 11) Save final artifact 💾
# ------------------------------------------------------------
save(betas_all, anno_combined,
     file = file.path("results", "betas_all.RData"))

message("Saved: results/betas_all.RData ✅")
message("🧃 Hydration reminder: go drink some water.")
