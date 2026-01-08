# R/helpers/config_tools.R
# ============================================================
# 🧭 Config Tools — make dataset-aware choices explicit
# ============================================================
# Why this helper exists 🤔
# ------------------------------------------------------------
# CNV pipelines *cannot* be fully “one-click universal” because:
#   ✅ CNV needs a biological reference (controls)
#   ✅ control definitions depend on *the dataset* (GSE labels vary)
#   ✅ your cases vary (1 case, 2 cases, 10 cases… all valid)
#
# So we do the sane thing:
#   1) Keep workflow code stable (src/03)
#   2) Put dataset-specific knobs in ONE place (R/config_cnv.R)
#   3) Validate config early (fail fast, fail loud) 🚨
# ============================================================

# ------------------------------------------------------------
# 0) Tiny utilities 🧰
# ------------------------------------------------------------

`%||%` <- function(x, y) if (!is.null(x)) x else y

is_named_character <- function(x) is.character(x) && !is.null(names(x)) && all(nzchar(names(x)))

pretty_vec <- function(x, max_n = 12) {
  x <- as.character(x)
  if (length(x) <= max_n) return(paste0(x, collapse = ", "))
  paste0(paste0(x[1:max_n], collapse = ", "), ", ... (+", length(x) - max_n, ")")
}

# ------------------------------------------------------------
# 1) Load config safely 📦
# ------------------------------------------------------------
# Usage:
#   cfg <- load_cnv_config()
#   cfg <- load_cnv_config("R/config_cnv.R")  # explicit path
#
# Notes:
# - We “source into a sandbox env” so it doesn’t pollute global env.
# - Returns a list (cfg) with only expected config fields.

load_cnv_config <- function(path = file.path("R", "config_cnv.R"), verbose = TRUE) {
  if (!file.exists(path)) {
    stop("🚨 CNV config not found: ", path, "\n",
         "Create it first at R/config_cnv.R")
  }
  
  env <- new.env(parent = emptyenv())
  sys.source(path, envir = env)
  
  # Pull only known fields (others ignored on purpose)
  cfg <- list(
    # required
    GSE_ID       = env$GSE_ID       %||% NULL,
    GSE_RAW_DIR  = env$GSE_RAW_DIR  %||% NULL,
    CTRL_CLASSES = env$CTRL_CLASSES %||% character(0),
    CASE_DIRS    = env$CASE_DIRS    %||% NULL,
    
    # optional goodies
    ARRAY_TYPE     = env$ARRAY_TYPE     %||% NULL,
    ZOOM_FAMILY    = env$ZOOM_FAMILY    %||% NULL,
    ZOOM_WINDOW    = env$ZOOM_WINDOW    %||% NULL,
    MIN_CONTROLS   = env$MIN_CONTROLS   %||% 10,
    GEO_CLASS_COL  = env$GEO_CLASS_COL  %||% NULL,
    GEO_MATERIAL_COL = env$GEO_MATERIAL_COL %||% NULL,
    
    OUT_DIR       = env$OUT_DIR       %||% "results/cnv",
    GENES_TO_CHECK = env$GENES_TO_CHECK %||% character(0),
    GENOME_BUILD  = env$GENOME_BUILD  %||% "hg19",
    SAVE_INTERMEDIATES = env$SAVE_INTERMEDIATES %||% TRUE
  )
  
  if (verbose) {
    message("📦 Loaded CNV config: ", path)
  }
  
  cfg
}

# ------------------------------------------------------------
# 2) Validate config
# ------------------------------------------------------------
# This is your “contract check” for dataset-aware choices.
#
# It checks:
# - required fields exist
# - directories exist (or warns)
# - cases are named and folders exist
# - zoom settings are well-formed
# - control classes provided (CNV should not run without them)

validate_cnv_config <- function(cfg, strict = FALSE) {
  # Required fields
  if (!is.character(cfg$GSE_ID) || !nzchar(cfg$GSE_ID)) {
    stop("🚨 GSE_ID must be a non-empty string. Example: GSE_ID <- \"GSE90496\"")
  }
  if (!is.character(cfg$GSE_RAW_DIR) || !nzchar(cfg$GSE_RAW_DIR)) {
    stop("🚨 GSE_RAW_DIR must be a non-empty path string.")
  }
  if (!is_named_character(cfg$CASE_DIRS) || length(cfg$CASE_DIRS) < 1) {
    stop("🚨 CASE_DIRS must be a named character vector.\n",
         "Example:\n  CASE_DIRS <- c(\"CaseA\"=\"data/caseA\", \"CaseB\"=\"data/caseB\")")
  }
  
  # Controls (CNV should scream if missing)
  if (!is.character(cfg$CTRL_CLASSES)) {
    stop("🚨 CTRL_CLASSES must be a character vector (can be empty while drafting).")
  }
  if (length(cfg$CTRL_CLASSES) == 0) {
    msg <- paste0(
      "🚨 CTRL_CLASSES is empty.\n",
      "CNV needs biological controls. Define them in R/config_cnv.R.\n",
      "Tip: start with control-like labels (often beginning with 'CONTR')."
    )
    if (strict) stop(msg) else warning(msg)
  }
  
  # Dir checks
  if (!dir.exists(cfg$GSE_RAW_DIR)) {
    msg <- paste0("⚠️ GSE_RAW_DIR not found (yet): ", cfg$GSE_RAW_DIR)
    if (strict) stop(msg) else warning(msg)
  }
  
  missing_cases <- cfg$CASE_DIRS[!dir.exists(cfg$CASE_DIRS)]
  if (length(missing_cases) > 0) {
    msg <- paste0(
      "⚠️ Some CASE_DIRS folders do not exist:\n",
      paste(names(missing_cases), "->", missing_cases, collapse = "\n")
    )
    if (strict) stop(msg) else warning(msg)
  }
  
  # Optional: zoom window format
  if (!is.null(cfg$ZOOM_WINDOW)) {
    ok <- is.list(cfg$ZOOM_WINDOW) &&
      all(c("xlim", "ylim") %in% names(cfg$ZOOM_WINDOW)) &&
      length(cfg$ZOOM_WINDOW$xlim) == 2 &&
      length(cfg$ZOOM_WINDOW$ylim) == 2
    
    if (!ok) {
      stop("🚨 ZOOM_WINDOW must be like:\n",
           "  ZOOM_WINDOW <- list(xlim = c(a,b), ylim = c(c,d))")
    }
  }
  
  invisible(TRUE)
}

# ------------------------------------------------------------
# 3) Print a friendly summary
# ------------------------------------------------------------
# This is optional but great for logs:
# “What dataset am I running? How many cases? Which controls?”

summarize_cnv_config <- function(cfg) {
  message("============================================================")
  message("🧬 CNV CONFIG SUMMARY (dataset-aware knobs) ")
  message("============================================================")
  message("📌 GSE_ID      : ", cfg$GSE_ID)
  message("📁 GSE_RAW_DIR : ", cfg$GSE_RAW_DIR)
  
  message("🧪 Cases (", length(cfg$CASE_DIRS), "):")
  for (nm in names(cfg$CASE_DIRS)) {
    message("   - ", nm, " -> ", cfg$CASE_DIRS[[nm]])
  }
  
  message("🛡️ Controls:")
  if (length(cfg$CTRL_CLASSES) == 0) {
    message("   - CTRL_CLASSES: (empty) ⚠️  (define before CNV)")
  } else {
    message("   - CTRL_CLASSES (", length(cfg$CTRL_CLASSES), "): ", pretty_vec(cfg$CTRL_CLASSES))
  }
  message("   - MIN_CONTROLS: ", cfg$MIN_CONTROLS)
  
  if (!is.null(cfg$ARRAY_TYPE))   message("🧫 ARRAY_TYPE  : ", cfg$ARRAY_TYPE)
  if (!is.null(cfg$ZOOM_FAMILY))  message("🔍 ZOOM_FAMILY : ", cfg$ZOOM_FAMILY)
  if (!is.null(cfg$ZOOM_WINDOW))  message("🪟 ZOOM_WINDOW : xlim=", paste(cfg$ZOOM_WINDOW$xlim, collapse = ", "),
                                          " | ylim=", paste(cfg$ZOOM_WINDOW$ylim, collapse = ", "))
  
  message("📦 OUT_DIR     : ", cfg$OUT_DIR)
  message("🧬 GENOME      : ", cfg$GENOME_BUILD)
  message("============================================================")
  invisible(cfg)
}

# ------------------------------------------------------------
# 4) CNV Check 🚨
# ------------------------------------------------------------
# Purpose:
# - Catch the classic CNV footguns BEFORE conumee2 runs
# - Prevent “silent wrong CNV” (the worst kind 😭)
#
# What it checks:
# ✅ config basic validity (calls validate_cnv_config)
# ✅ IDAT directories exist
# ✅ GEO metadata has the control label column
# ✅ CTRL_CLASSES exist in GEO metadata
# ✅ enough controls exist (>= MIN_CONTROLS)
# ✅ GSM ↔ paired IDAT matching is non-empty
# ✅ case folders contain paired IDATs (light check)
#
# Requires:
# - `anno` = pData(getGEO(...)) from GEOquery (GSEMatrix)
#
check_cnv_contract <- function(
    cfg,
    anno,
    strict = TRUE,
    verbose = TRUE
) {
  # 0) Basic config validation (your existing checker)
  validate_cnv_config(cfg, strict = strict)
  
  if (!is.data.frame(anno) || nrow(anno) == 0) {
    stop("🚨 `anno` must be a non-empty data.frame (from GEO pData).")
  }
  
  # 1) Pick GEO class column
  class_col <- cfg$GEO_CLASS_COL %||% "methylation class:ch1"
  if (!class_col %in% colnames(anno)) {
    stop(
      "🚨 GEO_CLASS_COL not found in anno: ", class_col, "\n",
      "Tip: set GEO_CLASS_COL in R/config_cnv.R\n",
      "Available columns (head): ", paste(head(colnames(anno), 30), collapse = ", ")
    )
  }
  
  class_vec <- as.character(anno[[class_col]])
  
  # 2) Controls exist in GEO metadata?
  if (length(cfg$CTRL_CLASSES) == 0) {
    stop("🚨 CTRL_CLASSES is empty. CNV cannot run without controls.")
  }
  
  ctrl_hits <- class_vec %in% cfg$CTRL_CLASSES
  n_ctrl <- sum(ctrl_hits, na.rm = TRUE)
  
  if (n_ctrl == 0) {
    stop(
      "🚨 No controls found in GEO metadata using CTRL_CLASSES.\n",
      "Fix: update CTRL_CLASSES (and/or GEO_CLASS_COL) in R/config_cnv.R.\n",
      "Debug tip: inspect available labels via:\n",
      "  sort(unique(anno[['", class_col, "']]))"
    )
  }
  
  min_controls <- cfg$MIN_CONTROLS %||% 10
  if (n_ctrl < min_controls) {
    msg <- paste0(
      "⚠️ Controls found = ", n_ctrl,
      " (< MIN_CONTROLS=", min_controls, ").\n",
      "CNV may be noisy/unstable. Consider broadening CTRL_CLASSES."
    )
    if (strict) warning(msg) else message(msg)
  }
  
  # 3) GSM ↔ IDAT matching (paired _Grn/_Red basenames)
  # We use `supplementary_file` to recover GSM basenames (common GEO pattern).
  if (!"supplementary_file" %in% colnames(anno)) {
    warning("⚠️ anno has no `supplementary_file` column → skipping GSM↔IDAT matching check.")
  } else {
    # Extract basenames like "GSMxxxxxxx" from supplementary_file entries
    gsm_all <- gsub("_Grn.*", "", gsub(".*suppl/", "", anno$supplementary_file))
    gsm_ctrl <- gsm_all[ctrl_hits]
    
    # Find paired IDAT basenames in GSE_RAW_DIR (must have BOTH Grn & Red)
    grn <- list.files(cfg$GSE_RAW_DIR, pattern = "_Grn\\.idat(\\.gz)?$", ignore.case = TRUE)
    red <- list.files(cfg$GSE_RAW_DIR, pattern = "_Red\\.idat(\\.gz)?$", ignore.case = TRUE)
    
    base_grn <- sub("_Grn\\.idat(\\.gz)?$", "", grn, ignore.case = TRUE)
    base_red <- sub("_Red\\.idat(\\.gz)?$", "", red, ignore.case = TRUE)
    paired_bases <- intersect(base_grn, base_red)
    
    overlap <- intersect(gsm_ctrl, paired_bases)
    
    if (length(paired_bases) == 0) {
      stop(
        "🚨 No paired IDAT basenames detected in GSE_RAW_DIR.\n",
        "Expected to find matched *_Grn.idat(.gz) AND *_Red.idat(.gz) files.\n",
        "Check: did you extract the raw IDAT archive correctly?"
      )
    }
    
    if (length(overlap) == 0) {
      stop(
        "🚨 GSM↔IDAT matching is EMPTY.\n",
        "Meaning: controls exist in GEO metadata, but their paired IDAT files are not found.\n",
        "Fix checklist:\n",
        " - Confirm you downloaded RAW IDATs for ", cfg$GSE_ID, "\n",
        " - Confirm GSE_RAW_DIR points to the extracted IDAT folder\n",
        " - Confirm files are paired (_Grn + _Red)\n",
        " - Confirm basenames match GSM IDs"
      )
    }
    
    if (verbose) {
      message("🧾 Contract: controls in GEO = ", n_ctrl)
      message("🧬 Contract: paired IDAT basenames in dir = ", length(paired_bases))
      message("🤝 Contract: control GSMs with paired IDATs = ", length(overlap), " ✅")
    }
  }
  
  # 4) Light check: cases contain IDATs (paired)
  for (nm in names(cfg$CASE_DIRS)) {
    d <- cfg$CASE_DIRS[[nm]]
    grn <- list.files(d, pattern = "_Grn\\.idat(\\.gz)?$", ignore.case = TRUE)
    red <- list.files(d, pattern = "_Red\\.idat(\\.gz)?$", ignore.case = TRUE)
    
    base_grn <- sub("_Grn\\.idat(\\.gz)?$", "", grn, ignore.case = TRUE)
    base_red <- sub("_Red\\.idat(\\.gz)?$", "", red, ignore.case = TRUE)
    paired <- intersect(base_grn, base_red)
    
    if (length(paired) == 0) {
      stop(
        "🚨 Case folder has no paired IDATs: ", nm, " -> ", d, "\n",
        "Expected *_Grn.idat(.gz) + *_Red.idat(.gz) pairs."
      )
    }
    
    if (verbose) message("🧪 Case ", nm, ": paired IDAT basenames = ", length(paired), " ✅")
  }
  
  if (verbose) message("✅ CNV CONTRACT PASSED. You may summon conumee2 now 😈➡️😇")
  invisible(TRUE)
}
