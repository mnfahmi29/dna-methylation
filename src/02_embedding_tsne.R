# src/02_embedding_tsne.R
# ============================================================
# 🧠✨ Embed methylation betas into 2D (t-SNE) for visualization
# ============================================================
#
# TL;DR 🫶
#   We take betas (samples × probes) → compress signal (PCA) → paint it in 2D (t-SNE).
#
# Input 🍱
#   - results/betas_all.RData  (from src/01_prework.R)
#       betas_all     : matrix (samples × probes)
#       anno_combined : data.frame aligned to betas_all rows
#
# Output 🎁
#   - results/tsne/tsne_results.rds   (embedding + params)
#   - results/tsne/tsne_plot.pdf     (publication-ish)
#   - results/tsne/tsne_plot.png     (GitHub-friendly)
#
# Important brain note 🧠⚠️
#   - t-SNE is NOT a classifier.
#   - t-SNE is basically “pretty gossip” for high-dimensional data 😭
#   - Use it to *inspect* structure, not to *prove* truth.
#
# Repro rules 🚦
#   - NO install.packages() here. (don’t summon dependency demons 👹)
#   - renv manages packages → run renv::restore() once (outside scripts!)
# ------------------------------------------------------------

rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Rtsne)    # the 2D magic portal 🌀
  library(ggplot2)  # pretty points = happy brain 🥹
})

# Helpers
source(file.path("R", "helpers", "probe_id_tools.R"))
source(file.path("R", "helpers", "plot_tsne_tools.R"))

# Optional: fast PCA helper (turbo mode 🚀)
if (file.exists(file.path("R", "MNPTraining", "RSpectra_pca.R"))) {
  source(file.path("R", "MNPTraining", "RSpectra_pca.R"))
}

ensure_dir("results")
ensure_dir(file.path("results", "tsne"))

# ------------------------------------------------------------
# 1) Load betas + annotation 🧬
# ------------------------------------------------------------
infile <- file.path("results", "betas_all.RData")
stopifnot(file.exists(infile))

load(infile)  # loads betas_all, anno_combined
stopifnot(exists("betas_all"), exists("anno_combined"))
stopifnot(nrow(betas_all) == nrow(anno_combined))

message("🍱 Loaded betas_all: ", nrow(betas_all), " samples × ", ncol(betas_all), " probes")
message("🏷️ Annotation rows: ", nrow(anno_combined), " (should match samples) ✅")

# ------------------------------------------------------------
# 2) Filter probes ✂️ (remove boring probes)
# ------------------------------------------------------------
# Why filter?
# - Probes with near-zero variance contribute almost nothing.
# - Too many probes = slow + noisy = laptop screaming 😵‍💫
#
# Strategy:
# - Keep top variable probes (common and effective).
n_top_probes <- min(20000, ncol(betas_all))  # tweak if needed 🛠️

message("✂️ Choosing top variable probes... (top ", n_top_probes, " by variance)")

probe_var <- apply(betas_all, 2, var, na.rm = TRUE)
keep_idx  <- order(probe_var, decreasing = TRUE)[seq_len(n_top_probes)]
betas_f   <- betas_all[, keep_idx, drop = FALSE]

message("✅ Kept probes: ", ncol(betas_f), " (signal-rich snacks 🍿)")

# t-SNE hates NA like vampires hate sunlight 🧛☀️
na_cols <- colSums(is.na(betas_f)) > 0
if (any(na_cols)) {
  message("⚠️ NA detected: dropping ", sum(na_cols), " probe columns (anti-drama)")
  betas_f <- betas_f[, !na_cols, drop = FALSE]
}

# ------------------------------------------------------------
# 3) PCA first 🧠➡️📦 (compress before teleportation)
# ------------------------------------------------------------
# PCA is like:
# “Let’s pack 20,000 probes into 50 suitcases” 🧳🧳🧳
#
# Then t-SNE runs on suitcases, not on the whole house.
n_pcs <- 50
n_pcs <- min(n_pcs, ncol(betas_f) - 1)

if (n_pcs < 2) stop("Not enough probes for PCA. Something is off 😭")

if (exists("run_pca_rspectra", mode = "function")) {
  message("⚡ PCA via RSpectra helper (fast & furious 🚗💨)")
  pca_out <- run_pca_rspectra(
    betas_f,
    n_components = n_pcs,
    center = TRUE,
    scale = FALSE
  )
  # Assumption: helper returns samples × PCs in $x (common pattern)
  X <- pca_out$x
} else {
  message("🐢 PCA via prcomp (still okay, just slower)")
  pca_out <- prcomp(betas_f, center = TRUE, scale. = FALSE)
  X <- pca_out$x[, seq_len(n_pcs), drop = FALSE]
}

stopifnot(nrow(X) == nrow(betas_f))
message("📦 PCA matrix ready: ", nrow(X), " samples × ", ncol(X), " PCs ✅")

# ------------------------------------------------------------
# 4) t-SNE 🌀 (the “2D vibe check”)
# ------------------------------------------------------------
# Parameter vibes 🎛️:
# - perplexity: neighborhood size (must be < (n-1)/3)
# - theta: speed/accuracy tradeoff (0=exact, 0.5=typical)
# - max_iter: stability (more = better, but slower)
set.seed(42)  # reproducible chaos 🎲✅

n <- nrow(X)
perplexity <- min(30, floor((n - 1) / 3) - 1)
perplexity <- max(perplexity, 5)  # don’t go too tiny unless dataset tiny

message("🌀 Launching t-SNE portal...")
message("🎛️ Params: perplexity=", perplexity, " | dims=2 | theta=0.5 | max_iter=2000")

tsne_out <- Rtsne(
  X,
  dims = 2,
  perplexity = perplexity,
  theta = 0.5,
  max_iter = 2000,
  verbose = TRUE,
  pca = FALSE,              # we already did PCA
  check_duplicates = FALSE
)

emb <- as.data.frame(tsne_out$Y)
colnames(emb) <- c("TSNE1", "TSNE2")
emb$sample <- rownames(betas_all)

# Attach annotation (safe join by sample name)
emb <- cbind(emb, anno_combined[emb$sample, , drop = FALSE])

message("✅ Embedding done! Your samples now live in 2D 🥹✨")

# ------------------------------------------------------------
# 5) Save embedding artifact 💾
# ------------------------------------------------------------
out_rds <- file.path("results", "tsne", "tsne_results.rds")
saveRDS(
  list(
    embedding = emb,
    params = list(
      n_top_probes = ncol(betas_f),
      n_pcs = ncol(X),
      perplexity = perplexity,
      seed = 42
    )
  ),
  out_rds
)
message("💾 Saved embedding: ", out_rds, " ✅")

# ------------------------------------------------------------
# 6) Plot 🎨 (GLOBAL + ZOOM) using helper tools 🧰✨
# ------------------------------------------------------------
# This section uses:
#   - check_tsne_contract()
#   - build_tsne_df()
#   - plot_tsne_global()
#   - plot_tsne_zoom()
#   - save_plot_png_pdf()
#
# All are defined in R/helpers/plot_tsne_tools.R :contentReference[oaicite:3]{index=3} :contentReference[oaicite:4]{index=4}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

# -----------------------------
# 6A) Decide which annotation column is the “class label” 🏷️
# -----------------------------
# Supports common GEO-ish naming patterns:
# - "methylation class:ch1" (your prework default)
# - "methylation.class.ch1" (dot-version)
class_col <- if ("methylation class:ch1" %in% colnames(anno_combined)) {
  "methylation class:ch1"
} else if ("methylation.class.ch1" %in% colnames(anno_combined)) {
  "methylation.class.ch1"
} else {
  stop("❌ Cannot find methylation class column in anno_combined. ",
       "Expected something like 'methylation class:ch1'. ",
       "Check colnames(anno_combined).")
}

# -----------------------------
# 6B) Define which samples are “cases” 🧪
# -----------------------------
# Works with:
# - 0 cases (pure reference)
# - 1 case
# - 2 cases
# - many cases
#
# ✅ IMPORTANT: cases are detected by *sample names* (rownames of betas_all).
# So if you have 5 cases, just list all 5 names here.
case_names <- c("Case_GBM", "Case_ICM")  # <-- EDIT ME (can be character(0))

# -----------------------------
# 6C) Optional: dataset-agnostic grouping (“macro labels”) 🗺️
# -----------------------------
# If your dataset is NOT GSE90496 (or the class system differs),
# you SHOULD build this yourself (or skip it entirely).
#
# ✅ group_map format:
#   - named character vector
#   - names = exact class labels in your annotation
#   - values = group label you want on the map
#
# Example (toy):
# group_map <- c(
#   "GBM, MES" = "Diffuse glioma",
#   "GBM, RTK I" = "Diffuse glioma",
#   "EPN, PF A" = "Ependymal"
# )
# If group_map = NULL → helper will just use class (or leave group NA for cases).
group_map <- NULL  # <-- USER DECIDES (recommended for non-GSE90496)
# If you use GSE90496, you could use this, this is my grouping:
group_map <- c(
  # Embryonal
  "ETMR" = "Embryonal",
  "MB, WNT" = "Embryonal",
  "MB, G3" = "Embryonal",
  "MB, G4" = "Embryonal",
  "MB, SHH CHL AD" = "Embryonal",
  "MB, SHH INF" = "Embryonal",
  "ATRT, MYC" = "Embryonal",
  "ATRT, SHH" = "Embryonal",
  "ATRT, TYR" = "Embryonal",
  "CNS NB, FOXR2" = "Embryonal",
  "HGNET, BCOR" = "Embryonal",
  "HGNET, MN1" = "Embryonal",
  
  # Glioblastoma / Diffuse glioma
  "DMG, K27" = "Glioblastoma / Diffuse glioma",
  "GBM, G34" = "Glioblastoma / Diffuse glioma",
  "GBM, MES" = "Glioblastoma / Diffuse glioma",
  "GBM, RTK I" = "Glioblastoma / Diffuse glioma",
  "GBM, RTK II" = "Glioblastoma / Diffuse glioma",
  "GBM, RTK III" = "Glioblastoma / Diffuse glioma",
  "GBM, MID" = "Glioblastoma / Diffuse glioma",
  "GBM, MYCN" = "Glioblastoma / Diffuse glioma",
  "IHG" = "Glioblastoma / Diffuse glioma",
  "PXA" = "Glioblastoma / Diffuse glioma",
  
  # Glioneuronal
  "CN" = "Glioneuronal",
  "DLGNT" = "Glioneuronal",
  "LIPN" = "Glioneuronal",
  "LGG, DIG/DIA" = "Glioneuronal",
  "LGG, DNT" = "Glioneuronal",
  "LGG, RGNT" = "Glioneuronal",
  "LGG, GG" = "Glioneuronal",
  "RETB" = "Glioneuronal",
  "LGG, PA PF" = "Glioneuronal",
  "LGG, PA MID" = "Glioneuronal",
  "LGG, PA/GG ST" = "Glioneuronal",
  "LGG, SEGA" = "Glioneuronal",
  "LGG, MYB" = "Glioneuronal",
  
  # Ependymal
  "EPN, RELA" = "Ependymal",
  "EPN, PF A" = "Ependymal",
  "EPN, PF B" = "Ependymal",
  "EPN, SPINE" = "Ependymal",
  "EPN, YAP" = "Ependymal",
  "EPN, MPE" = "Ependymal",
  "SUBEPN, PF" = "Ependymal",
  "SUBEPN, SPINE" = "Ependymal",
  "SUBEPN, ST" = "Ependymal",
  
  # Mesenchymal
  "CHORDM" = "Mesenchymal",
  "EWS" = "Mesenchymal",
  "HMB" = "Mesenchymal",
  "MNG" = "Mesenchymal",
  "SFT HMPC" = "Mesenchymal",
  "EFT, CIC" = "Mesenchymal",
  "CPH, ADM" = "Mesenchymal",
  "CPH, PAP" = "Mesenchymal",
  
  # Plexus
  "PLEX, AD" = "Plexus",
  "PLEX, PED A" = "Plexus",
  "PLEX, PED B" = "Plexus",
  
  # Glioma IDH
  "A IDH" = "Glioma IDH",
  "A IDH, HG" = "Glioma IDH",
  "O IDH" = "Glioma IDH",
  
  # Melanocytic
  "MELAN" = "Melanocytic",
  "MELCYT" = "Melanocytic",
  
  # Hematopoietic
  "LYMPHO" = "Hematopoietic",
  "PLASMA" = "Hematopoietic",
  
  # Pituitary / pineal-ish
  "PITUI" = "Other",
  "PGG, nC" = "Other",
  "CHGL" = "Other",
  "ENB, A" = "Other",
  "ENB, B" = "Other",
  "ANA PA" = "Other",
  "PTPR, A" = "Other",
  "PTPR, B" = "Other",
  "PIN T,  PB A" = "Other",
  "PIN T,  PB B" = "Other",
  "PIN T, PPT" = "Other",
  "PITAD, ACTH" = "Other",
  "PITAD, PRL" = "Other",
  "PITAD, FSH LH" = "Other",
  "PITAD, TSH" = "Other",
  "PITAD, STH DNS A" = "Other",
  "PITAD, STH DNS B" = "Other",
  "PITAD, STH SPA" = "Other",
  "SCHW" = "Other",
  "SCHW, MEL" = "Other"
)


# -----------------------------
# 6D) Contract check (early screaming 😱✅)
# -----------------------------
# Ensures:
# - emb has sample column
# - anno has class_col
# - all emb samples exist in anno rownames
check_tsne_contract(emb, anno_combined, class_col)  # :contentReference[oaicite:5]{index=5}

# -----------------------------
# 6E) Build plotting df (coords + class + case flags + optional group) 🧩
# -----------------------------
df_plot <- build_tsne_df(
  emb        = emb,
  anno       = anno_combined,
  class_col  = class_col,
  case_names = case_names,
  group_map  = group_map,
  group_col_name = "group"
)  # :contentReference[oaicite:6]{index=6}

message("🧪 Case count check:"); print(table(df_plot$is_case))

# -----------------------------
# 6F) GLOBAL plot 🌍
# -----------------------------
# label_mode options in helper:
# - "none"
# - "group_onmap"      (nice when group_map exists)
# - "class_centroids"  (labels classes; can be clutter)
label_mode <- if (!is.null(group_map)) "group_onmap" else "none"

p_global <- plot_tsne_global(
  df         = df_plot,
  label_mode = label_mode,
  label_col  = "group",     # used only when label_mode="group_onmap"
  label_min_n = 25,
  title      = "t-SNE (reference + your cases) 🌍🧬"
)

save_plot_png_pdf(
  p_global,
  out_prefix = file.path("results", "tsne", "tsne_global"),
  width = 20, height = 8, dpi = 600
)  # :contentReference[oaicite:7]{index=7}

message("🖼️ Saved GLOBAL plot: results/tsne/tsne_global.(png|pdf) ✅")

# -----------------------------
# 6G) ZOOM plot 🔍
# -----------------------------
# Zoom strategies (pick ONE):
# A) around cases (default, if cases exist)
# B) by specific classes (set zoom_classes)
#
# If you have no cases, you should set zoom_classes manually.
zoom_classes <- NULL  # e.g. c("GBM, MES", "GBM, RTK I")

do_zoom_by_cases <- any(df_plot$is_case)
do_zoom_by_class <- !is.null(zoom_classes)

if (do_zoom_by_cases || do_zoom_by_class) {
  p_zoom <- plot_tsne_zoom(
    df                 = df_plot,
    zoom_classes        = zoom_classes,   # used if not NULL
    zoom_cases          = do_zoom_by_cases,
    pad_frac            = 0.03,
    label_classes_in_zoom = TRUE,
    class_label_min_n   = 10,
    title               = "t-SNE zoom 🔍"
  )
  
  save_plot_png_pdf(
    p_zoom,
    out_prefix = file.path("results", "tsne", "tsne_zoom"),
    width = 20, height = 8, dpi = 600
  )  # :contentReference[oaicite:8]{index=8}
  
  message("🖼️ Saved ZOOM plot: results/tsne/tsne_zoom.(png|pdf) ✅")
} else {
  message("🧠 Zoom skipped: no cases found and zoom_classes is NULL. ",
          "Set zoom_classes if you want a reference-only zoom.")
}

message("🎉 Plotting complete! Your repo now has GLOBAL + ZOOM maps 🧬✨")