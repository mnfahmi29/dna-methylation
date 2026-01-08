# R/helpers/cnv_tools.R
# ------------------------------------------------------------
# 🧬 CNV helper toolbox (GEO controls + IDAT matching + conumee2)
# ------------------------------------------------------------
# Philosophy:
# - Make the annoying parts deterministic and testable
# - Keep src/03_cnv.R short and readable
#
# You’ll typically use these helpers in this order:
#   1) pick_controls_from_geo()
#   2) match_gse_controls_to_idats()
#   3) run_conumee2_case()
#   4) export_cnv_outputs()
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(GEOquery)
  library(minfi)
  library(dplyr)
})

# -----------------------------
# 0) Small utilities
# -----------------------------
ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ============================================================
# 1) Control selection from GEO
# ============================================================
# pick_controls_from_geo(GSE_ID, ctrl_classes)
#
# Returns:
# - anno: pData from GEO series matrix
# - keep_gsm: GSM ids that match your control class labels
#
# Notes:
# - Different GSE have different column names.
# - We allow you to specify which column holds the "class".
#
pick_controls_from_geo <- function(GSE_ID,
                                   ctrl_classes,
                                   class_col = c("methylation class:ch1", "methylation.class.ch1"),
                                   getGPL = FALSE) {
  stopifnot(is.character(GSE_ID), length(GSE_ID) == 1)
  stopifnot(is.character(ctrl_classes), length(ctrl_classes) >= 1)
  
  message("📥 GEO: downloading metadata for ", GSE_ID, " ...")
  gse <- getGEO(GSE_ID, GSEMatrix = TRUE, getGPL = getGPL)
  
  # Usually the first element is the series matrix object
  gse_obj <- gse[[1]]
  anno <- pData(gse_obj)
  
  # find an existing class column
  class_col <- class_col[class_col %in% colnames(anno)]
  if (length(class_col) == 0) {
    stop("Cannot find any class_col in GEO pData(). Available cols include: ",
         paste(colnames(anno)[1:min(30, ncol(anno))], collapse = ", "),
         if (ncol(anno) > 30) " ...")
  }
  class_col <- class_col[1]
  
  cls <- as.character(anno[[class_col]])
  keep <- which(cls %in% ctrl_classes)
  
  keep_gsm <- rownames(anno)[keep]
  message("✅ Controls selected: ", length(keep_gsm), " GSMs (class_col = '", class_col, "')")
  
  list(
    anno = anno,
    class_col = class_col,
    keep_gsm = keep_gsm
  )
}

# ============================================================
# 2) IDAT pairing + GSM matching
# ============================================================
# match_gse_controls_to_idats(GSE_RAW_DIR, keep_gsm)
#
# Goal:
# - Given a folder with unpacked IDATs, find basenames that correspond to GSMs.
#
# Returns:
# - idat_basenames: character vector of basenames (without _Grn.idat/_Red.idat)
# - idat_files: data.frame with basename + paths + status
#
# Assumptions (common GEO convention):
# - IDAT basenames contain GSM somewhere (e.g., GSM12345_....)
#
match_gse_controls_to_idats <- function(
    GSE_RAW_DIR,
    keep_gsm,
    pattern_grn = "_Grn\\.idat(\\.gz)?$",
    pattern_red = "_Red\\.idat(\\.gz)?$"
) {
  stopifnot(dir.exists(GSE_RAW_DIR))
  stopifnot(is.character(keep_gsm), length(keep_gsm) >= 1)
  
  grn <- list.files(GSE_RAW_DIR, pattern = pattern_grn, full.names = TRUE, ignore.case = TRUE)
  red <- list.files(GSE_RAW_DIR, pattern = pattern_red, full.names = TRUE, ignore.case = TRUE)
  
  if (length(grn) == 0 || length(red) == 0) {
    stop("No IDAT pairs found in ", GSE_RAW_DIR,
         ". Expect _Grn.idat(.gz) and _Red.idat(.gz)")
  }
  
  base_grn <- sub(pattern_grn, "", basename(grn), ignore.case = TRUE)
  base_red <- sub(pattern_red, "", basename(red), ignore.case = TRUE)
  base_both <- intersect(base_grn, base_red)
  
  # GSM match (vector-safe): any GSM substring present
  gsm_pat <- paste(unique(keep_gsm), collapse = "|")
  kept_basenames <- base_both[grepl(gsm_pat, base_both)]
  
  message("✅ IDAT basenames with both channels: ", length(base_both))
  message("✅ Matched to keep_gsm: ", length(kept_basenames))
  
  # Detect whether .gz exists for each basename by checking the real files
  # (prefer .gz if present)
  grn_path <- function(b) {
    p1 <- file.path(GSE_RAW_DIR, paste0(b, "_Grn.idat.gz"))
    p2 <- file.path(GSE_RAW_DIR, paste0(b, "_Grn.idat"))
    if (file.exists(p1)) p1 else p2
  }
  red_path <- function(b) {
    p1 <- file.path(GSE_RAW_DIR, paste0(b, "_Red.idat.gz"))
    p2 <- file.path(GSE_RAW_DIR, paste0(b, "_Red.idat"))
    if (file.exists(p1)) p1 else p2
  }
  
  idat_files <- data.frame(
    basename = kept_basenames,
    grn = vapply(kept_basenames, grn_path, character(1)),
    red = vapply(kept_basenames, red_path, character(1)),
    stringsAsFactors = FALSE
  )
  
  ok <- file.exists(idat_files$grn) & file.exists(idat_files$red)
  if (!all(ok)) {
    warning("⚠️ Some matched basenames are missing a channel. Dropping them.")
    idat_files <- idat_files[ok, , drop = FALSE]
  }
  
  list(idat_basenames = idat_files$basename, idat_files = idat_files)
}

# ============================================================
# 3) Run conumee2 per case
# ============================================================
# run_conumee2_case(data_cases, data_ctrl, anno, case_name)
#
# Inputs:
# - data_cases: RGset (or path) for a single case, already read
# - data_ctrl : RGset for controls (reference)
# - anno      : annotation object required by conumee2 (e.g., conumee2::EPIC)
# - case_name : label
#
# Returns:
# - cnv object (whatever conumee2 returns)
#
run_conumee2_case <- function(data_case_rgset,
                              data_ctrl_rgset,
                              conumee_anno,
                              case_name = "Case_1",
                              do_fit = TRUE,
                              do_bin = TRUE) {
  if (missing(data_case_rgset) || missing(data_ctrl_rgset)) {
    stop("Provide RGsets for case + controls.")
  }
  
  if (!requireNamespace("conumee2", quietly = TRUE)) {
    stop("conumee2 not installed in this renv environment.")
  }
  
  message("🧬 conumee2: running case = ", case_name)
  
  # Build CNV object
  cnv <- conumee2::CNV.load(
    case = data_case_rgset,
    control = data_ctrl_rgset,
    annotation = conumee_anno
  )
  
  if (do_fit) cnv <- conumee2::CNV.fit(cnv)
  if (do_bin) cnv <- conumee2::CNV.bin(cnv)
  
  cnv
}

# ============================================================
# 4) Export bundle (segments + focal exploration hooks)
# ============================================================
# export_cnv_outputs(cnv_obj, out_dir, case_name, genes_to_check)
#
# Saves:
# - cnv_obj RDS
# - segments table (if available)
# - simple focal query table (genes_to_check)
#
export_cnv_outputs <- function(cnv_obj,
                               out_dir,
                               case_name = "Case_1",
                               genes_to_check = character(0)) {
  ensure_dir(out_dir)
  
  # 1) Save full object
  saveRDS(cnv_obj, file.path(out_dir, paste0(case_name, "_cnv_obj.rds")))
  
  # 2) Try to pull segments / bins if conumee2 exposes them
  seg <- NULL
  bin <- NULL
  
  # Many conumee-style objects store results in slots or list entries.
  # We keep this defensive and non-breaking.
  if (!is.null(cnv_obj$segments)) seg <- cnv_obj$segments
  if (!is.null(cnv_obj$bin))      bin <- cnv_obj$bin
  
  if (!is.null(seg)) {
    utils::write.csv(seg, file.path(out_dir, paste0(case_name, "_segments.csv")), row.names = FALSE)
  }
  if (!is.null(bin)) {
    utils::write.csv(bin, file.path(out_dir, paste0(case_name, "_bins.csv")), row.names = FALSE)
  }
  
  # 3) Focal exploration scaffold (gene-based lookup)
  # We only create a "template" table here, because gene mapping depends on your annotation choice.
  if (length(genes_to_check) > 0) {
    focal_tbl <- data.frame(
      case = case_name,
      gene = genes_to_check,
      note = "Fill this using your chosen gene->locus mapping (hg19/hg38) + overlap with segments/bins",
      stringsAsFactors = FALSE
    )
    utils::write.csv(focal_tbl, file.path(out_dir, paste0(case_name, "_focal_genes_template.csv")),
                     row.names = FALSE)
  }
  
  message("✅ CNV bundle exported to: ", out_dir)
  invisible(TRUE)
}

# ============================================================
# ADD-ON BLOCK (for src/03_cnv_conumee2.R compatibility) 🧩
# ============================================================
# These wrappers make src/03 clean:
# - read_controls_rgset()
# - read_cases_rgset()
# - trim_methylset_probes_if_needed()
# - run_conumee2_case_pipeline()
# - export_cnv_bundle()
#
# They reuse your existing helpers:
# - match_gse_controls_to_idats()
# - run_conumee2_case()
# - export_cnv_outputs()
# ============================================================

# 2.5) Small helper: attach raw dir to idat basenames (so we can read later)
attach_gse_raw_dir <- function(idat_basenames, GSE_RAW_DIR) {
  attr(idat_basenames, "GSE_RAW_DIR") <- GSE_RAW_DIR
  idat_basenames
}

# ------------------------------------------------------------
# A) Read controls RGset from basenames
# ------------------------------------------------------------
# Accepts:
# - character vector of basenames (preferred), with attr("GSE_RAW_DIR")
# - OR list output from match_gse_controls_to_idats() (idat_files present)
read_controls_rgset <- function(idat_basenames_or_list, verbose = TRUE) {
  if (is.list(idat_basenames_or_list) && !is.null(idat_basenames_or_list$idat_files)) {
    files <- sub("(_Grn\\.idat(\\.gz)?)$", "", idat_basenames_or_list$idat_files$grn, ignore.case = TRUE)
  } else {
    basenames <- idat_basenames_or_list$idat_files$basename
    GSE_RAW_DIR <- dirname(idat_basenames_or_list$idat_files$grn[1])
    files <- file.path(GSE_RAW_DIR, basenames)
  }
  
  if (verbose) message("🧊 Reading controls: ", length(files), " paired IDAT basenames")
  minfi::read.metharray(files, verbose = verbose)
}


# ------------------------------------------------------------
# B) Read cases RGset from CASE_DIRS + CASE_NAMES
# ------------------------------------------------------------
# CASE_DIRS: vector of directories (each can contain 1+ samples)
# CASE_NAMES: labels; if a folder has >1 sample we auto-prefix
read_cases_rgset <- function(CASE_DIRS, CASE_NAMES = NULL, verbose = TRUE) {
  stopifnot(length(CASE_DIRS) >= 1)
  if (is.null(CASE_NAMES)) CASE_NAMES <- basename(CASE_DIRS)
  stopifnot(length(CASE_DIRS) == length(CASE_NAMES))
  
  rg_list <- vector("list", length(CASE_DIRS))
  
  for (i in seq_along(CASE_DIRS)) {
    d <- CASE_DIRS[[i]]
    nm <- CASE_NAMES[[i]]
    
    if (!dir.exists(d)) stop("🚨 CASE_DIR does not exist: ", d)
    
    if (verbose) message("🧫 Reading case dir: ", nm, " -> ", d)
    rg <- minfi::read.metharray.exp(d, verbose = verbose)
    
    # If folder contains multiple samples, prefix with case name
    if (ncol(rg) >= 1) {
      colnames(rg) <- if (ncol(rg) == 1) {
        nm
      } else {
        paste0(nm, "__", colnames(rg))
      }
    }
    
    rg_list[[i]] <- rg
  }
  
  # combine into one RGChannelSet
  RG_cases <- do.call(cbind, rg_list)
  RG_cases
}

# ------------------------------------------------------------
# C) EPICv2 “probe suffix” safety for MethylSet
# ------------------------------------------------------------
# IMPORTANT:
# - This does NOT change intensities
# - It only standardizes probe rownames (cg suffix trimming)
trim_methylset_probes_if_needed <- function(Mset, verbose = TRUE) {
  # only run if we detect suffix-like patterns
  ids <- rownames(Mset)
  has_suffix <- any(grepl("^cg\\d+_.+", ids))
  if (!has_suffix) {
    if (verbose) message("✅ Probe IDs look clean (no EPICv2 cg suffix detected).")
    return(Mset)
  }
  
  if (!file.exists(file.path("R", "helpers", "probe_id_tools.R"))) {
    stop("🚨 probe_id_tools.R not found. Needed for trimming.")
  }
  source(file.path("R", "helpers", "probe_id_tools.R"))
  
  if (!exists("trim_methylation_object_probes", mode = "function")) {
    stop("🚨 trim_methylation_object_probes() not found after sourcing probe_id_tools.R")
  }
  
  if (verbose) message("✂️ EPICv2 suffix detected → trimming cg IDs for compatibility.")
  trim_methylation_object_probes(Mset, verbose = verbose)
}

# ------------------------------------------------------------
# D) conumee2 pipeline per case name
# ------------------------------------------------------------
# Inputs:
# - Mset_cases: MethylSet with ALL cases (columns are samples)
# - Mset_ctrl : MethylSet with controls (columns are controls)
# - conumee_anno: CNV annotation object
# Returns:
# - list(cnv = <object>, segments = <df or NULL>, focal = <df or NULL>)
run_conumee2_case_pipeline <- function(Mset_cases,
                                       Mset_ctrl,
                                       conumee_anno,
                                       case_name,
                                       do_fit = TRUE,
                                       do_bin = TRUE,
                                       do_detail = TRUE,
                                       do_segment = TRUE,
                                       do_focal = TRUE) {
  if (!requireNamespace("conumee2", quietly = TRUE)) {
    stop("conumee2 not installed in this renv environment.")
  }
  
  if (!case_name %in% colnames(Mset_cases)) {
    stop("🚨 case_name not found in Mset_cases colnames(): ", case_name)
  }
  
  # subset a single case as a one-column MethylSet
  case_one <- Mset_cases[, case_name, drop = FALSE]
  
  cnv <- conumee2::CNV.load(case = case_one, control = Mset_ctrl, annotation = conumee_anno)
  if (do_fit) cnv <- conumee2::CNV.fit(cnv)
  if (do_bin) cnv <- conumee2::CNV.bin(cnv)
  
  # optional expansions (version-safe with tryCatch)
  detail <- NULL
  seg <- NULL
  focal <- NULL
  
  if (do_detail) {
    detail <- tryCatch(conumee2::CNV.detail(cnv), error = function(e) NULL)
  }
  if (do_segment) {
    seg <- tryCatch(conumee2::CNV.segment(cnv), error = function(e) NULL)
  }
  if (do_focal) {
    focal <- tryCatch(conumee2::CNV.focal(cnv), error = function(e) NULL)
  }
  
  list(cnv = cnv, detail = detail, segments = seg, focal = focal)
}

# ------------------------------------------------------------
# E) Export a per-case CNV bundle (plots + tables)
# ------------------------------------------------------------
export_cnv_bundle <- function(cnv_res,
                              out_dir,
                              case_name,
                              genes_to_check = character(0)) {
  ensure_dir(out_dir)
  
  # Save the whole result object
  saveRDS(cnv_res, file.path(out_dir, paste0(case_name, "_cnv_res.rds")))
  
  # Segments table (if present)
  if (!is.null(cnv_res$segments)) {
    utils::write.csv(cnv_res$segments,
                     file.path(out_dir, paste0(case_name, "_segments.csv")),
                     row.names = FALSE)
  }
  
  # Focal table (if present)
  if (!is.null(cnv_res$focal)) {
    utils::write.csv(cnv_res$focal,
                     file.path(out_dir, paste0(case_name, "_focal.csv")),
                     row.names = FALSE)
  }
  
  # Genome plot (if available)
  if (requireNamespace("conumee2", quietly = TRUE)) {
    pdf(file.path(out_dir, paste0(case_name, "_genomeplot.pdf")), width = 12, height = 4)
    try(conumee2::CNV.genomeplot(cnv_res$cnv), silent = TRUE)
    dev.off()
    
    pdf(file.path(out_dir, paste0(case_name, "_detailplot.pdf")), width = 12, height = 6)
    try(conumee2::CNV.detailplot(cnv_res$cnv), silent = TRUE)
    dev.off()
  }
  
  # Optional: gene list “note” template (dataset-specific mapping)
  if (length(genes_to_check) > 0) {
    focal_tpl <- data.frame(
      case = case_name,
      gene = genes_to_check,
      note = "Map genes → genomic loci (hg19/hg38), then overlap with segments/focal output.",
      stringsAsFactors = FALSE
    )
    utils::write.csv(focal_tpl,
                     file.path(out_dir, paste0(case_name, "_focal_genes_template.csv")),
                     row.names = FALSE)
  }
  
  message("✅ Exported CNV bundle for ", case_name, " → ", out_dir)
  invisible(TRUE)
}

