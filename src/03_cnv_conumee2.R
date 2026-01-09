# src/03_cnv_conumee2.R
# ============================================================
# 🧬🔍 CNV calling with conumee2 — Segments + Focals
# ============================================================
#
# What this script does 🧠
# ------------------------
# 1) Pick controls from GEO metadata 🧊
# 2) Match those GSM controls to IDAT pairs in your GSE_RAW_DIR 🧩
# 3) Read controls + your case IDAT(s) 🧫
# 4) preprocessRaw() because CNV wants intensities ⚡
# 5) Run conumee2 per case:
#       fit → bin → detail → segment → focal
# 6) Export per-case bundle:
#       segments CSV + focal CSV + genome plots + focal gene PDFs 🎁
#
# Repro rules 🚦
# --------------
# ❌ NO install.packages() here
# ✅ renv already prepared the kitchen 🍳
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(minfi)
  library(conumee2)
})

# ------------------------------------------------------------
# 0) Load helpers 🧰
# ------------------------------------------------------------
source(file.path("helpers", "cnv_tools.R"))
source(file.path("helpers", "probe_id_tools.R"))  # for trimming cg suffixes safely

# ------------------------------------------------------------
# 1) Dataset knobs 🧭🧠
# ------------------------------------------------------------
GSE_ID      <- "GSE90496"
GSE_RAW_DIR <- "data/GSE90496_RAW"

CASE_DIRS <- c(
  "data/Case GBM_MES",
  "data/Intracranial Mesenchymal FET_CREB Fusion"
)

CASE_NAMES <- c("Case_GBM", "Case_ICM")

# CNV output base
OUT_DIR <- file.path("results", "cnv")

# ✅ user-owned: what counts as “control” in THIS dataset?
CTRL_CLASSES <- c("REACT","PONS","PINEAL","INFLAM","HYPTHAL","HEMI","CEBM","WM","ADENOPIT")

# Optional: focal genes you want quick plots for
GENES_TO_CHECK <- c("EGFR", "CDKN2A/B", "CDK4", "MDM2", "MET", "CDK6")

# Optional: how strict you want the pipeline to be
MIN_CONTROLS <- 10

# ------------------------------------------------------------
# 2) Checkpoint screening (fail fast, fail loud) 🚨
# ------------------------------------------------------------
ensure_dir(OUT_DIR)

stopifnot(dir.exists(GSE_RAW_DIR))
stopifnot(length(CASE_DIRS) == length(CASE_NAMES))
stopifnot(all(dir.exists(CASE_DIRS)))

# “Scream early” checks that prevent silent wrong CNV
check_cnv_contract(
  GSE_ID = GSE_ID,
  GSE_RAW_DIR = GSE_RAW_DIR,
  CTRL_CLASSES = CTRL_CLASSES,
  CASE_DIRS = CASE_DIRS,
  CASE_NAMES = CASE_NAMES,
  min_controls = MIN_CONTROLS
)

message("📁 GSE_RAW_DIR : ", GSE_RAW_DIR)
message("🧪 Cases      : ", paste(CASE_NAMES, collapse = ", "))
message("🛡️ Controls   : ", length(CTRL_CLASSES), " labels (user-owned)")

# ------------------------------------------------------------
# 3) Pick controls from GEO metadata 🧊
# ------------------------------------------------------------
geo <- pick_controls_from_geo(
  GSE_ID = GSE_ID,
  ctrl_classes = CTRL_CLASSES
)

# ------------------------------------------------------------
# 4) Match controls → paired IDAT basenames 🧩
# ------------------------------------------------------------
idats <- match_gse_controls_to_idats(
  GSE_RAW_DIR = GSE_RAW_DIR,
  keep_gsm = geo$keep_gsm
)

if (length(idats$idat_basenames) < MIN_CONTROLS) {
  warning("⚠️ Controls with paired IDATs: ", length(idats$idat_basenames),
          " (< ", MIN_CONTROLS, "). CNV reference may be noisy.")
}

# ------------------------------------------------------------
# 5) Read IDATs (controls + cases) 🧫
# ------------------------------------------------------------
RG_ctrl  <- read_controls_rgset(idats$idat_basenames, verbose = TRUE)
RG_cases <- read_cases_rgset(CASE_DIRS, CASE_NAMES, verbose = TRUE)

# ------------------------------------------------------------
# 6) preprocessRaw() because CNV wants intensities ⚡
# ------------------------------------------------------------
message("⚡ preprocessRaw() for controls + cases ...")
Mset_ctrl  <- preprocessRaw(RG_ctrl)
Mset_cases <- preprocessRaw(RG_cases)

# EPICv2 safety: trim cg########_suffix → cg######## (ID-only, no intensity changes)
Mset_ctrl  <- trim_methylation_object_probes(Mset_ctrl,  verbose = TRUE, drop_duplicates = TRUE)
Mset_cases <- trim_methylation_object_probes(Mset_cases, verbose = TRUE, drop_duplicates = TRUE)

# ------------------------------------------------------------
# 7) conumee2 annotation 🗺️
# ------------------------------------------------------------
data(exclude_regions)
data(detail_regions)

anno <- CNV.create_anno(
  array_type      = c("450k", "EPIC", "EPICv2"),
  exclude_regions = exclude_regions,
  detail_regions  = detail_regions
)

# ------------------------------------------------------------
# 8) Run CNV per case + export bundle 🎁
# ------------------------------------------------------------
cnv_all <- list()

for (nm in CASE_NAMES) {
  message("\n============================================================")
  message("🚀 Running CNV for: ", nm)
  message("============================================================")
  
  res <- run_conumee2_case_pipeline(
    Mset_cases   = Mset_cases,
    Mset_ctrl    = Mset_ctrl,
    conumee_anno = anno,
    case_name    = nm,
    do_detail    = TRUE,
    do_segment   = TRUE,
    do_focal     = TRUE
  )
  
  cnv_all[[nm]] <- res
  
  export_cnv_bundle(
    cnv_res = res,
    out_dir = OUT_DIR,
    case_name = nm,
    genes_to_check = GENES_TO_CHECK
  )
}

saveRDS(cnv_all, file.path(OUT_DIR, "cnv_all_cases.rds"))
message("\n🎉 CNV done! Bundles exported to: ", OUT_DIR)
message("🧃 Go hydrate. Your CNVs will still be there when you return.")