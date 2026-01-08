# R/helpers/plot_tsne_tools.R
# ============================================================
# 🎨 t-SNE plotting tools (dataset-agnostic)
# ============================================================
# You provide:
#  - emb: data.frame with TSNE1, TSNE2, sample
#  - anno: data.frame with rownames = sample IDs
#
# This helper builds:
#  - a plotting df (coords + class + case flags + optional group)
#  - global plot
#  - zoom plot (by bbox, by selected classes, or around cases)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ------------------------------------------------------------
# 0) Contract check: embedding ↔ annotation sanity 🛑
# ------------------------------------------------------------
check_tsne_contract <- function(emb, anno, class_col) {
  
  stopifnot(is.data.frame(emb), is.data.frame(anno))
  
  if (!"sample" %in% colnames(emb))
    stop("❌ emb must contain a 'sample' column.")
  
  if (!class_col %in% colnames(anno))
    stop("❌ class_col not found in anno_combined: ", class_col)
  
  if (anyDuplicated(emb$sample))
    stop("❌ Duplicate sample IDs in embedding.")
  
  missing <- setdiff(emb$sample, rownames(anno))
  if (length(missing) > 0)
    stop("❌ Samples in embedding but not in annotation: ",
         paste(missing, collapse = ", "))
  
  message("✅ Contract check passed: embedding ↔ annotation are aligned.")
}

# -----------------------------
# 1) Build plotting dataframe 🧩
# -----------------------------
build_tsne_df <- function(
    emb,
    anno,
    class_col,
    case_names = character(),
    group_map = NULL,
    group_col_name = "group"
) {
  stopifnot(all(c("TSNE1", "TSNE2", "sample") %in% colnames(emb)))
  stopifnot(class_col %in% colnames(anno))
  
  # join by sample (robust)
  stopifnot(all(emb$sample %in% rownames(anno)))
  class_vec <- as.character(anno[emb$sample, class_col])
  
  df <- data.frame(
    x = emb$TSNE1,
    y = emb$TSNE2,
    sample = as.character(emb$sample),
    class = class_vec,
    stringsAsFactors = FALSE
  )
  
  # mark cases
  case_names <- as.character(case_names)
  df$case_id <- ifelse(df$sample %in% case_names, df$sample, NA_character_)
  df$is_case <- !is.na(df$case_id)
  
  # do NOT let cases pollute reference labels
  df$class[df$is_case] <- NA
  
  # optional grouping (dataset-agnostic)
  # group_map: named character vector, names=class, values=group
  if (!is.null(group_map)) {
    stopifnot(is.character(group_map), !is.null(names(group_map)))
    df[[group_col_name]] <- unname(group_map[df$class])
    df[[group_col_name]][is.na(df[[group_col_name]]) & !is.na(df$class)] <- "Other"
  } else {
    df[[group_col_name]] <- ifelse(is.na(df$class), NA_character_, "Ungrouped")
  }
  
  df
}

# -----------------------------
# 2) Auto grouping helper 🧠
# -----------------------------
# Makes a simple group_map by taking top-N most frequent classes as themselves,
# everything else -> "Other". (No biology assumptions.)
make_topn_group_map <- function(class_vec, top_n = 8, other_label = "Other") {
  class_vec <- as.character(class_vec)
  class_vec <- class_vec[!is.na(class_vec)]
  tab <- sort(table(class_vec), decreasing = TRUE)
  top_classes <- names(tab)[seq_len(min(top_n, length(tab)))]
  gm <- setNames(top_classes, top_classes)  # group = class
  list(group_map = gm, top_classes = top_classes, counts = tab)
}

# -----------------------------
# 3) Label positions (median centroids) 📍
# -----------------------------
compute_label_centroids <- function(df, label_col, min_n = 25) {
  stopifnot(label_col %in% colnames(df))
  df %>%
    filter(!is_case, !is.na(.data[[label_col]])) %>%
    group_by(.data[[label_col]]) %>%
    summarise(
      x = median(x, na.rm = TRUE),
      y = median(y, na.rm = TRUE),
      n = dplyr::n(),
      .groups = "drop"
    ) %>%
    filter(n >= min_n) %>%
    rename(label = .data[[label_col]])
}

# -----------------------------
# 4) Global plot 🌍
# -----------------------------
plot_tsne_global <- function(
    df,
    label_mode = c("none", "group_onmap", "class_centroids"),
    label_col = "group",
    label_min_n = 25,
    case_colors = NULL,
    title = "t-SNE (reference + cases)"
) {
  label_mode <- match.arg(label_mode)
  
  # cases: assign colors automatically if not supplied
  if (is.null(case_colors)) {
    case_ids <- sort(unique(df$case_id[!is.na(df$case_id)]))
    if (length(case_ids) > 0) {
      # default ggplot palette, but stable order
      case_colors <- setNames(scales::hue_pal()(length(case_ids)), case_ids)
    } else {
      case_colors <- character()
    }
  }
  
  p <- ggplot(df, aes(x, y)) +
    geom_point(
      data = subset(df, !is_case),
      aes(color = class),
      size = 1.1, alpha = 0.85
    ) +
    guides(color = "none") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    ggtitle(title)
  
  if (any(df$is_case)) {
    p <- p +
      geom_point(
        data = subset(df, is_case),
        aes(fill = case_id),
        shape = 21, color = "black", stroke = 1.2, size = 5.5
      ) +
      scale_fill_manual(values = case_colors, guide = "none") +
      geom_text_repel(
        data = subset(df, is_case),
        aes(label = case_id),
        size = 4, max.overlaps = Inf
      )
  }
  
  if (label_mode == "group_onmap") {
    lab <- compute_label_centroids(df, label_col = label_col, min_n = label_min_n)
    p <- p +
      geom_label_repel(
        data = lab,
        aes(x = x, y = y, label = label),
        size = 4,
        alpha = 0.9,
        label.size = 0.25,
        box.padding = 0.6,
        point.padding = 0.2,
        min.segment.length = 0,
        max.overlaps = Inf
      )
  }
  
  if (label_mode == "class_centroids") {
    lab <- compute_label_centroids(df, label_col = "class", min_n = label_min_n)
    p <- p +
      geom_text_repel(
        data = lab,
        aes(x = x, y = y, label = label),
        size = 3,
        color = "black",
        box.padding = 0.25,
        point.padding = 0.15,
        min.segment.length = 0,
        max.overlaps = Inf
      )
  }
  
  p
}

# -----------------------------
# 5) Zoom plot 🔍
# -----------------------------
# zoom can be defined by:
#   A) bbox: xlim/ylim
#   B) classes: zoom to classes set (uses their bounding box)
#   C) cases: zoom around cases with padding
plot_tsne_zoom <- function(
    df,
    zoom_classes = NULL,
    zoom_cases = TRUE,
    xlim = NULL,
    ylim = NULL,
    pad_frac = 0.03,
    label_classes_in_zoom = TRUE,
    class_label_min_n = 10,
    case_colors = NULL,
    title = "t-SNE zoom"
) {
  # figure out bbox
  if (!is.null(xlim) && !is.null(ylim)) {
    xlim_use <- xlim
    ylim_use <- ylim
  } else if (!is.null(zoom_classes)) {
    region <- df %>% filter(!is_case, class %in% zoom_classes)
    if (nrow(region) < 5) warning("Zoom classes have too few points.")
    xlim_use <- range(region$x, na.rm = TRUE)
    ylim_use <- range(region$y, na.rm = TRUE)
  } else if (zoom_cases && any(df$is_case)) {
    region <- df %>% filter(is_case)
    xlim_use <- range(region$x, na.rm = TRUE)
    ylim_use <- range(region$y, na.rm = TRUE)
  } else {
    stop("Provide xlim/ylim OR zoom_classes OR set zoom_cases=TRUE with cases present.")
  }
  
  # padding
  dx <- diff(xlim_use); dy <- diff(ylim_use)
  if (dx == 0) dx <- 1
  if (dy == 0) dy <- 1
  xlim_use <- xlim_use + c(-1, 1) * dx * pad_frac
  ylim_use <- ylim_use + c(-1, 1) * dy * pad_frac
  
  # case colors auto
  if (is.null(case_colors)) {
    case_ids <- sort(unique(df$case_id[!is.na(df$case_id)]))
    if (length(case_ids) > 0) case_colors <- setNames(scales::hue_pal()(length(case_ids)), case_ids)
    else case_colors <- character()
  }
  
  # build plot
  p <- ggplot(df, aes(x, y)) +
    geom_point(
      data = subset(df, !is_case),
      color = "grey85", size = 1.1, alpha = 0.6
    ) +
    geom_point(
      data = subset(df, !is_case,
                    x >= xlim_use[1], x <= xlim_use[2],
                    y >= ylim_use[1], y <= ylim_use[2]),
      aes(color = class),
      size = 1.3, alpha = 0.9
    ) +
    guides(color = "none") +
    coord_cartesian(xlim = xlim_use, ylim = ylim_use) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    ggtitle(title)
  
  if (any(df$is_case)) {
    p <- p +
      geom_point(
        data = subset(df, is_case),
        aes(fill = case_id),
        shape = 21, color = "black", stroke = 1.2, size = 5.5
      ) +
      scale_fill_manual(values = case_colors, guide = "none") +
      geom_text_repel(
        data = subset(df, is_case),
        aes(label = case_id),
        size = 4, max.overlaps = Inf
      )
  }
  
  if (label_classes_in_zoom) {
    class_lab <- df %>%
      filter(!is_case, !is.na(class)) %>%
      group_by(class) %>%
      summarise(
        x = median(x, na.rm = TRUE),
        y = median(y, na.rm = TRUE),
        n = n(),
        .groups = "drop"
      ) %>%
      filter(
        n >= class_label_min_n,
        x >= xlim_use[1], x <= xlim_use[2],
        y >= ylim_use[1], y <= ylim_use[2]
      )
    
    p <- p +
      geom_text_repel(
        data = class_lab,
        aes(x = x, y = y, label = class),
        size = 3,
        box.padding = 0.25,
        point.padding = 0.15,
        min.segment.length = 0,
        max.overlaps = Inf
      )
  }
  
  p
}

# -----------------------------
# 6) Save helpers 💾
# -----------------------------
save_plot_png_pdf <- function(p, out_prefix, width = 20, height = 8, dpi = 600) {
  ggsave(paste0(out_prefix, ".png"), p, width = width, height = height, dpi = dpi)
  ggsave(paste0(out_prefix, ".pdf"), p, width = width, height = height)
  invisible(TRUE)
}
