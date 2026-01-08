rm(list=ls())
options(stringsAsFactors = FALSE)

if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")

# ---- packages ----
pkgs_bioc <- c(
  "minfi","GEOquery","limma",
  "IlluminaHumanMethylation450kmanifest",
  "IlluminaHumanMethylationEPICmanifest"
)
pkgs_cran <- c("Rtsne","RSpectra","ggplot2","data.table")
install.packages("renv")

BiocManager::install(setdiff(pkgs_bioc, rownames(installed.packages())), ask=FALSE)
install.packages(setdiff(pkgs_cran, rownames(installed.packages())), dependencies=TRUE)

library(minfi)
library(GEOquery)
library(limma)
library(Rtsne)
library(RSpectra)
library(ggplot2)
library(data.table)
library(renv)

# Your helper (optional)
if (file.exists(file.path("R","MNPprocessIDAT_functions.R"))) {
  source(file.path("R","MNPprocessIDAT_functions.R"))
}

dir.create("results", showWarnings=FALSE)

# ------------------ PATHS (EDIT) ------------------
GSE_RAW_DIR  <- "data/GSE90496_RAW"  # extracted *.idat pairs
CASE_GBM_DIR <- "data/Case GBM_MES"
CASE_ICM_DIR <- "data/Intracranial Mesenchymal FET_CREB Fusion"

# ------------------ Helpers ------------------
trim_cg_ids <- function(Mset) {
  old_ids <- rownames(Mset)
  new_ids <- sub("(^cg[0-9]+).*", "\\1", old_ids)
  dup <- duplicated(new_ids)
  if (any(dup)) {
    warning(sprintf("Dropping %d duplicated probes after EPICv2 trim.", sum(dup)))
    Mset <- Mset[!dup, ]
    new_ids <- new_ids[!dup]
  }
  rownames(Mset) <- new_ids
  Mset
}

preprocess_any <- function(RGset) {
  if (exists("MNPpreprocessIllumina")) {
    MNPpreprocessIllumina(RGset)
  } else {
    preprocessIllumina(RGset)
  }
}

batch_correct_meth_unmeth <- function(methy, unmethy, material) {
  material <- trimws(as.character(material))
  has_ffpe   <- any(material == "FFPE",   na.rm=TRUE)
  has_frozen <- any(material == "Frozen", na.rm=TRUE)
  
  if (!(has_ffpe && has_frozen)) {
    message("Batch correction skipped (need BOTH FFPE and Frozen).")
    return(list(methy.ba=methy, unmethy.ba=unmethy, batch=NULL))
  }
  
  batch <- ifelse(material == "FFPE", 2, 1)
  methy.ba   <- 2^removeBatchEffect(log2(methy + 1), batch)
  unmethy.ba <- 2^removeBatchEffect(log2(unmethy + 1), batch)
  
  list(methy.ba=methy.ba, unmethy.ba=unmethy.ba, batch=batch)
}

# ------------------ 2) GEO annotation ------------------
gse  <- getGEO("GSE90496", GSEMatrix=TRUE, getGPL=FALSE)
anno <- pData(gse$GSE90496_series_matrix.txt.gz)

# Choose reference subset (RECOMMENDED to reduce RAM):
# Option A: GBM only (as you were doing)
# anno_ref <- anno[grep("^GBM", anno$`methylation class:ch1`), , drop=FALSE]
# USE ALL CLASSES
anno_ref <- anno


# (If you want all classes later, set: anno_ref <- anno)

fname <- gsub("_Grn.*", "", gsub(".*suppl/", "", anno_ref$supplementary_file))
filepath_gse_all <- file.path(GSE_RAW_DIR, fname)

# ------------------ 3) Read your EPICv2 cases ------------------
RG_gbm <- read.metharray.exp(CASE_GBM_DIR, verbose=TRUE)
RG_icm <- read.metharray.exp(CASE_ICM_DIR, verbose=TRUE)
colnames(RG_gbm) <- "Case_GBM"
colnames(RG_icm) <- "Case_ICM"

# Minimal annotation (ONLY what you need)
anno_gse_small <- data.frame(
  row.names = fname,  # IMPORTANT: will align after we read GSE in chunks
  `methylation class:ch1` = anno_ref$`methylation class:ch1`,
  `material:ch1`          = anno_ref$`material:ch1`,
  stringsAsFactors = FALSE
)
anno_cases <- data.frame(
  row.names = c("Case_GBM","Case_ICM"),
  `methylation class:ch1` = c("Case_GBM","Case_ICM"),
  `material:ch1`          = c("Frozen","Frozen"),  # change if FFPE
  stringsAsFactors = FALSE
)

# ------------------ 4) CHUNKING: read + preprocess GSE safely ------------------
# Why chunk? Reading 300+ IDATs at once can blow RAM and crash RStudio.
chunk_size <- 25
chunks <- split(filepath_gse_all, ceiling(seq_along(filepath_gse_all) / chunk_size))

betas_chunks <- vector("list", length(chunks))
sample_names_chunks <- vector("list", length(chunks))

message("Reading + preprocessing GSE in chunks...")
for (i in seq_along(chunks)) {
  message(sprintf("Chunk %d / %d", i, length(chunks)))
  
  RG_chunk <- read.metharray(chunks[[i]], verbose=FALSE)
  M_chunk  <- preprocess_any(RG_chunk)
  
  # Reference is 450k/EPICv1 so no suffix trim needed here
  
  meth <- getMeth(M_chunk)
  unm  <- getUnmeth(M_chunk)
  beta <- meth / (meth + unm + 100)
  
  # store as samples x probes
  betas_chunks[[i]] <- t(beta)
  sample_names_chunks[[i]] <- colnames(meth)
  
  rm(RG_chunk, M_chunk, meth, unm, beta)
  gc()
}

betas_gse <- do.call(rbind, betas_chunks)
rownames(betas_gse) <- unlist(sample_names_chunks)

# ------------------ 5) Preprocess cases and compute betas ------------------
M_gbm <- trim_cg_ids(preprocess_any(RG_gbm))
M_icm <- trim_cg_ids(preprocess_any(RG_icm))

beta_gbm <- t(getMeth(M_gbm) / (getMeth(M_gbm) + getUnmeth(M_gbm) + 100))
beta_icm <- t(getMeth(M_icm) / (getMeth(M_icm) + getUnmeth(M_icm) + 100))

# ------------------ 6) Harmonize probes + bind samples ------------------
common_probes <- Reduce(intersect, list(
  colnames(betas_gse),
  colnames(beta_gbm),
  colnames(beta_icm)
))

betas_all <- rbind(
  betas_gse[, common_probes, drop=FALSE],
  beta_gbm[,  common_probes, drop=FALSE],
  beta_icm[,  common_probes, drop=FALSE]
)

# Build combined annotation aligned to betas_all rows
anno_ref2 <- data.frame(
  row.names = rownames(betas_gse),
  `methylation class:ch1` = anno_ref$`methylation class:ch1`,
  `material:ch1`          = anno_ref$`material:ch1`,
  stringsAsFactors = FALSE
)
anno_combined <- rbind(anno_ref2, anno_cases)
anno_combined <- anno_combined[rownames(betas_all), , drop=FALSE]

# ------------------ 7) Batch correction on meth/unmeth? ------------------
# If you only have betas, do NOT try to do the original meth/unmeth batch correction.
# (That correction is defined on meth/unmeth intensities.)
# For now: skip if you are running betas-only reference.
message("NOTE: using betas-only reference; skipping meth/unmeth batch correction.")

save(betas_all, anno_combined, file=file.path("results","betas_all.RData"))