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

Sys.which("make")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
options(repos = BiocManager::repositories())

BiocManager::install(version = "3.20", ask = FALSE)
options(timeout = 600)

BiocManager::install(c(
  "DNAcopy",
  "InteractionSet",
  "nullranges",
  "org.Hs.eg.db",
  "FDb.InfiniumMethylation.hg19",
  "RnBeads"
), ask = FALSE, force = TRUE)

if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("hovestadtlab/conumee2", subdir = "conumee2")

library(conumee2)
library("minfi")
library(GEOquery)
library(tidyverse)
library(dplyr)
library(stringr)
library(IlluminaHumanMethylationEPICv2manifest)
library(magrittr)

GBM <- "data/Case GBM_MES"
ICM <- "data/Intracranial Mesenchymal FET_CREB Fusion"
GSE  <- "data/GSE90496_RAW"

# ---- reference classes from GSE90496 ----
ctrl_classes <- c("REACT","PONS","PINEAL","INFLAM","HYPTHAL","HEMI","CEBM","WM","ADENOPIT")
pattern <- paste0("\\b(", paste(ctrl_classes, collapse="|"), ")\\b")

# ---- 1) Load GSE matrix ----
gse1 <- getGEO("GSE90496", GSEMatrix = TRUE, getGPL = FALSE)
eset <- gse1[[1]]
anno_gse <- pData(eset)

# ---- 2) Build a small metadata ----
mclass_col <- grep("^methylation class", colnames(anno_gse), value = TRUE)
if (length(mclass_col) == 0) stop("No 'methylation class' column found in pData(). Check colnames(anno_gse).")

anno_gse_small <- data.frame(
  gsm = rownames(anno_gse),  # GSM IDs
  methylation_class = anno_gse[[mclass_col[1]]],
  stringsAsFactors = FALSE
)

# ---- 3) GREP FILTER FIRST ----
keep_gsm <- anno_gse_small$gsm[grep(pattern, anno_gse_small$methylation_class, ignore.case = TRUE)]
keep_gsm <- unique(keep_gsm)

cat("GSM kept by class filter:", length(keep_gsm), "\n")

idat_files <- list.files(GSE, pattern="\\.(idat|idat\\.gz)$", full.names=TRUE)

grn_files <- idat_files[grepl("_Grn\\.idat(\\.gz)?$", idat_files)]
red_files <- idat_files[grepl("_Red\\.idat(\\.gz)?$", idat_files)]

grn_base <- sub("_Grn\\.idat(\\.gz)?$", "", grn_files)
red_base <- sub("_Red\\.idat(\\.gz)?$", "", red_files)

get_gsm <- function(x) regmatches(x, regexpr("GSM[0-9]+", x))

grn_df <- data.frame(gsm=get_gsm(grn_base), base=grn_base)
red_df <- data.frame(gsm=get_gsm(red_base), base=red_base)

both <- merge(grn_df, red_df, by=c("gsm","base"))
both_keep <- both[both$gsm %in% keep_gsm, ]

basenames_keep <- both_keep$base
basenames_keep <- basenames_keep[!duplicated(basenames_keep)]

# ---- 6) Read only the controls ----
RG_gse_ctrl <- read.metharray(basenames_keep, verbose = TRUE)
cat("Loaded control samples:", ncol(RG_gse_ctrl), "\n")

# Read Metharray Case
RG_gbm <- read.metharray.exp(GBM, verbose = TRUE)
RG_icm <- read.metharray.exp(ICM,  verbose = TRUE)
colnames(RG_gbm) <- "Case_GBM"
colnames(RG_icm) <- "Case_ICM"

# For CNV, stay on raw intensities
MSet_gbm <- preprocessRaw(RG_gbm)
MSet_icm <- preprocessRaw(RG_icm)
Mset_ctrl <- preprocessRaw(RG_gse_ctrl)

#=========================================

# -------------------------------
# Create conumee2 annotation
# -------------------------------
data(exclude_regions)
data(detail_regions)  # hg19 example detail regions (works for 450k/EPIC, generally OK for overlap)

anno <- CNV.create_anno(
  array_type      = c("450k", "EPICv2"),
  exclude_regions = exclude_regions,
  detail_regions  = detail_regions
)

# -------------------------------
# EPICv2 probe naming
# -------------------------------
# If you see cg00000029_TC21 style rownames, conumee2 overlap may still work,
# but trimming can prevent accidental probe drop. We'll do a safe trim only if needed.

trim_cg_ids_if_needed <- function(Mset) {
  ids <- rownames(Mset)
   #Detect suffix patterns like cg########_...
  if (any(grepl("^cg[0-9]+_", ids))) {
    new_ids <- sub("(^cg[0-9]+).*", "\\1", ids)
    dup <- duplicated(new_ids)
    if (any(dup)) {
      warning(sprintf("Dropping %d duplicated probes after EPICv2 trim", sum(dup)))
      Mset <- Mset[!dup, , drop=FALSE]
      new_ids <- new_ids[!dup]
    }
    rownames(Mset) <- new_ids
  }
  Mset
}

MSet_gbm  <- trim_cg_ids_if_needed(MSet_gbm)
MSet_icm  <- trim_cg_ids_if_needed(MSet_icm)
Mset_ctrl <- trim_cg_ids_if_needed(Mset_ctrl)

# -------------------------------
# Build CNV data objects
# -------------------------------

get_total_intensity <- function(MSet, sample_name = NULL) {
  M <- getMeth(MSet)
  U <- getUnmeth(MSet)
  tot <- M + U
  
  # safety net if only one sample
  tot <- as.matrix(tot)
  if (!is.null(sample_name)) {
    colnames(tot) <- sample_name
  } else if (is.null(colnames(tot))) {
    colnames(tot) <- paste0("S", seq_len(ncol(tot)))
  }
  tot
}

I_gbm  <- get_total_intensity(MSet_gbm,  "Case_GBM")
I_icm  <- get_total_intensity(MSet_icm,  "Case_ICM")
I_ctrl <- get_total_intensity(Mset_ctrl)

# -------------------------------
# Align probe
# -------------------------------
common_ids <- Reduce(intersect, list(rownames(I_gbm), rownames(I_icm), rownames(I_ctrl)))
cat("Shared probes (GBM/ICM/CTRL):", length(common_ids), "\n")

I_gbm  <- I_gbm [common_ids, , drop = FALSE]
I_icm  <- I_icm [common_ids, , drop = FALSE]
I_ctrl <- I_ctrl[common_ids, , drop = FALSE]

# -------------------------------
# NOW create CNV.data objects
# -------------------------------
data_gbm  <- CNV.load(as.data.frame(I_gbm))
data_icm  <- CNV.load(as.data.frame(I_icm))
data_ctrl <- CNV.load(as.data.frame(I_ctrl))

# -------------------------------
# Run CNV pipeline (GBM)
# -------------------------------
x_gbm <- CNV.fit(data_gbm, data_ctrl, anno)
x_gbm <- CNV.bin(x_gbm)
x_gbm <- CNV.detail(x_gbm)
x_gbm <- CNV.segment(x_gbm)
# focal calling:
x_gbm_foc <- CNV.focal(x_gbm)

# -------------------------------
# Run CNV pipeline (ICM)
# -------------------------------
x_icm <- CNV.fit(data_icm, data_ctrl, anno)
x_icm <- CNV.bin(x_icm)
x_icm <- CNV.detail(x_icm)
x_icm <- CNV.segment(x_icm)
# focal calling:
x_icm_foc <- CNV.focal(x_icm)

save(x_gbm, x_gbm_foc, x_icm, x_icm_foc, file=file.path("results","cnv.RData"))

# -------------------------------
# Output: plots + segments
# -------------------------------
dir.create("results", showWarnings = FALSE)
load(file.path("results","cnv.RData"))

# Segments
pdf("results/CNV_GBM_seg.pdf", width = 12, height = 4)
CNV.genomeplot(x_gbm[1])
dev.off()

pdf("results/CNV_ICM_seg.pdf", width = 12, height = 4)
CNV.genomeplot(x_icm[1])
dev.off()

# List all available focal genes
unique(CNV.write(x_gbm_foc, what = "focal")$name)
unique(CNV.write(x_icm_foc, what = "focal")$name)

# Focal
pdf("results/CNV_GBM_foc.pdf", width = 12, height = 4)
CNV.genomeplot(x_gbm_foc[1])
dev.off()

pdf("results/CNV_ICM_foc.pdf", width = 12, height = 4)
CNV.genomeplot(x_icm_foc[1])
dev.off()

# Focal Save
safe_filename <- function(x) {
  x %>%
    gsub("/", "_", .) %>%
    gsub("\\s+", "_", .) %>%
    gsub("[^A-Za-z0-9_\\-]", "", .)
}


genes_to_check <- c(
  "EGFR",
  "CDKN2A/B",
  "CDK4",
  "MDM2",
  "MET",
  "CDK6"
)

dir.create("results/focal", showWarnings = FALSE, recursive = TRUE)

for (g in genes_to_check) {
  g_file <- safe_filename(g)
  
  pdf(file.path("results/focal", paste0("GBM_", g_file, ".pdf")),
      width = 8, height = 5)
  
  try(
    CNV.detailplot(x_gbm_foc, name = g),
    silent = TRUE
  )
  
  dev.off()
}



# Focal EGFR
CNV.detailplot(x_gbm_foc, name = "EGFR")
CNV.write(x_gbm_foc, what = "focal")
pdf("results/GBM_EGFR_focal.pdf", width = 8, height = 5)
CNV.detailplot(i_gbm_foc, name = "EGFR")
dev.off()

CNV.detailplot(x_icm_foc, name = "EGFR")
CNV.write(x_icm_foc, what = "focal")
pdf("results/ICM_EGFR_focal.pdf", width = 8, height = 5)
CNV.detailplot(i_icm_foc, name = "EGFR")
dev.off()

# Focal CDKN 2A/2B
CNV.detailplot(x_gbm_foc, name = "CDKN 2")
CNV.write(x_gbm_foc, what = "focal")
pdf("results/GBM_EGFR_focal.pdf", width = 8, height = 5)
CNV.detailplot(i_gbm_foc, name = "EGFR")
dev.off()


CNV.detailplot(x_icm_foc, name = "EGFR")
CNV.write(x_icm_foc, what = "focal")
pdf("results/ICM_EGFR_focal.pdf", width = 8, height = 5)
CNV.detailplot(i_icm_foc, name = "EGFR")
dev.off()



# Save segments tables
seg_gbm <- CNV.write(x_gbm, what = "segments")
seg_icm <- CNV.write(x_icm, what = "segments")
write.csv(seg_gbm, "results/segments_GBM.csv", row.names = FALSE)
write.csv(seg_icm, "results/segments_ICM.csv", row.names = FALSE)

focal_gbm <- CNV.write(x_gbm, what = "focal")
focal_icm <- CNV.write(x_icm, what = "focal")
write.csv(focal_gbm, "results/focal_GBM.csv", row.names = FALSE)
write.csv(focal_icm, "results/focal_ICM.csv", row.names = FALSE)

