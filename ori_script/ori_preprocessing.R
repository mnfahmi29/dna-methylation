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

if (file.exists(file.path("R","MNPprocessIDAT_functions.R"))) {
  source(file.path("R","MNPprocessIDAT_functions.R"))
}

dir.create("results", showWarnings=FALSE)

# ------------------ PATHS  ------------------
GSE_RAW_DIR  <- "data/GSE90496_RAW" 
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

# ------------------ Probe filtering ------------------
filter_probes <- function(Mset, filter_dir="filter") {
  fl_amb  <- file.path(filter_dir, "amb_3965probes.vh20151030.txt")
  fl_epic <- file.path(filter_dir, "epicV1B2_32260probes.vh20160325.txt")
  fl_snp  <- file.path(filter_dir, "snp_7998probes.vh20151030.txt")
  fl_xy   <- file.path(filter_dir, "xy_11551probes.vh20151030.txt")
  
  if (!all(file.exists(c(fl_amb, fl_epic, fl_snp, fl_xy)))) {
    warning("Probe filtering skipped: missing one or more filter files under '", filter_dir, "'.")
    return(Mset)
  }
  
  message("probe filtering ...", Sys.time())
  amb.filter  <- read.table(fl_amb,  header=FALSE)
  epic.filter <- read.table(fl_epic, header=FALSE)
  snp.filter  <- read.table(fl_snp,  header=FALSE)
  xy.filter   <- read.table(fl_xy,   header=FALSE)
  
  rs.filter <- grep("^rs", rownames(Mset))
  ch.filter <- grep("^ch", rownames(Mset))
  
  remove <- unique(c(
    match(amb.filter[,1],  rownames(Mset)),
    match(epic.filter[,1], rownames(Mset)),
    match(snp.filter[,1],  rownames(Mset)),
    match(xy.filter[,1],   rownames(Mset)),
    rs.filter,
    ch.filter
  ))
  remove <- remove[!is.na(remove)]
  
  message(sprintf("Removing %d / %d probes (%.2f%%).",
                  length(remove), nrow(Mset), 100*length(remove)/max(1,nrow(Mset))))
  Mset[-remove, ]
}

# ------------------ 2) GEO annotation ------------------
gse  <- getGEO("GSE90496", GSEMatrix=TRUE, getGPL=FALSE)
anno_ref <- pData(gse$GSE90496_series_matrix.txt.gz)

fname <- gsub("_Grn.*", "", gsub(".*suppl/", "", anno_ref$supplementary_file))
filepath_gse_all <- file.path(GSE_RAW_DIR, fname)

# ------------------ 3) Read EPICv2 cases ------------------
RG_gbm <- read.metharray.exp(CASE_GBM_DIR, verbose=TRUE)
RG_icm <- read.metharray.exp(CASE_ICM_DIR, verbose=TRUE)
colnames(RG_gbm) <- "Case_GBM"
colnames(RG_icm) <- "Case_ICM"

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
chunk_size <- 25
chunks <- split(filepath_gse_all, ceiling(seq_along(filepath_gse_all) / chunk_size))

betas_chunks <- vector("list", length(chunks))
sample_names_chunks <- vector("list", length(chunks))

message("Reading + preprocessing GSE in chunks...")
for (i in seq_along(chunks)) {
  message(sprintf("Chunk %d / %d", i, length(chunks)))
  
  RG_chunk <- read.metharray(chunks[[i]], verbose=FALSE)
  M_chunk  <- preprocess_any(RG_chunk)
  M_chunk  <- filter_probes(M_chunk)
  
  # Reference is 450k/EPICv1 so no suffix trimming
  
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
M_gbm <- filter_probes(trim_cg_ids(preprocess_any(RG_gbm)))
M_icm <- filter_probes(trim_cg_ids(preprocess_any(RG_icm)))

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

# ------------------ 7) Save to Directory ------------------
message("NOTE: using betas-only reference; skipping meth/unmeth batch correction.")

save(betas_all, anno_combined, file=file.path("results","betas_all.RData"))