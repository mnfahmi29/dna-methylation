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

library(RSpectra)
source(file.path("R","RSpectra_pca.R"))

y <- as.factor(anno_combined$`methylation class:ch1`)

k_probes <- min(32000, ncol(betas_all))
ord <- order(-apply(betas_all, 2, sd))
betas_f <- betas_all[, ord[1:k_probes], drop=FALSE]

k_pc <- min(94, nrow(betas_f)-1, ncol(betas_f)-1)
pca <- prcomp_svds(betas_f, k=k_pc)

# t-SNE on PCs
set.seed(42)
res <- Rtsne(pca$x, pca=FALSE, max_iter=2500, theta=0, verbose=TRUE)

save(pca, res, file=file.path("results","pca_tsne.RData"))


# ------------------ 9) Plot (colored classes + styled cases) ------------------
library(ggplot2)

class_vec <- as.character(anno_combined$methylation.class.ch1)
stopifnot(length(class_vec) == nrow(res$Y))

df <- data.frame(
  x      = res$Y[,1],
  y      = res$Y[,2],
  class  = class_vec,
  sample = rownames(betas_all),
  stringsAsFactors = FALSE
)
df$is_case <- df$class %in% c("Case_GBM", "Case_ICM")
table(df$is_case)


my_colors <- c(
  "Case_GBM" = "#D62728",
  "Case_ICM" = "#1F77B4"
)

p <- ggplot(df, aes(x = x, y = y)) +
  
  # Reference samples
  geom_point(
    data = subset(df, !is_case),
    aes(color = class),
    size = 1.4,
    alpha = 0.7
  ) +
  
  # Your cases (big, black outline)
  geom_point(
    data = subset(df, is_case),
    aes(fill = class),
    shape  = 21,
    color  = "black",
    stroke = 1.2,
    size   = 6
  ) +
  
  # Labels for your cases
  geom_text(
    data = subset(df, is_case),
    aes(label = sample),
    vjust = -1.1,
    size  = 4
  ) +
  
  scale_fill_manual(values = my_colors, drop = FALSE) +
  theme_bw() +
  ggtitle("t-SNE (GSE90496 reference + EPICv2 cases)") +
  theme(
    legend.position = "right",
    panel.grid = element_blank()
  )

print(p)

ggsave(
  file.path("results", "tsne_cases_styled.png"),
  p, width = 12, height = 8, dpi = 200
)


# ------------------ 9) Zoom In ------------------

library(ggplot2)
library(cowplot)

# df already has: x, y, class, sample, is_case
# my_colors already defined for Case_GBM / Case_ICM

# ---- main plot (keep legend) ----
p_main <- ggplot(df, aes(x, y)) +
  geom_point(data=subset(df, !is_case),
             aes(color=class), size=1.4, alpha=0.7) +
  geom_point(data=subset(df, is_case),
             aes(fill=class), shape=21, color="black",
             stroke=1.2, size=6) +
  geom_text(data=subset(df, is_case),
            aes(label=sample), vjust=-1.1, size=4) +
  scale_fill_manual(values=my_colors, drop=FALSE) +
  theme_bw() +
  ggtitle("t-SNE (GSE90496 reference + EPICv2 cases)") +
  theme(panel.grid = element_blank())

# ---- compute zoom window around your 2 cases ----
cases_xy <- subset(df, is_case)
pad_x <- diff(range(cases_xy$x)) * 4; if (pad_x == 0) pad_x <- 10
pad_y <- diff(range(cases_xy$y)) * 4; if (pad_y == 0) pad_y <- 10

xlim_zoom <- range(cases_xy$x) + c(-pad_x, pad_x)
ylim_zoom <- range(cases_xy$y) + c(-pad_y, pad_y)

# ---- zoom plot (NO legends!) ----
p_zoom <- ggplot(df, aes(x, y)) +
  geom_point(data=subset(df, !is_case),
             color="grey80", size=1.2, alpha=0.7) +
  geom_point(data=subset(df, is_case),
             aes(fill=class), shape=21, color="black",
             stroke=1.2, size=6) +
  geom_text(data=subset(df, is_case),
            aes(label=sample), vjust=-1.1, size=4) +
  scale_fill_manual(values=my_colors, drop=FALSE) +
  coord_cartesian(xlim=xlim_zoom, ylim=ylim_zoom) +
  theme_bw() +
  theme(
    legend.position = "none",
    plot.title = element_text(size=10),
    plot.margin = margin(2,2,2,2)
  ) +
  ggtitle("Zoom: your cases")

# ---- inset it (draw once) ----
p_inset <- ggdraw(p_main) +
  draw_plot(p_zoom, x = 0.05, y = 0.62, width = 0.38, height = 0.33)

ggsave("results/tsne_with_zoom_inset.png", p_inset, width=12, height=8, dpi=200)
print(p_inset)

#=============================================================


class_vec <- as.character(anno_combined$methylation.class.ch1)
stopifnot(length(class_vec) == nrow(res$Y))

df <- data.frame(
  x      = res$Y[,1],
  y      = res$Y[,2],
  class  = class_vec,
  sample = rownames(betas_all),
  stringsAsFactors = FALSE
)
df$is_case <- df$class %in% c("Case_GBM", "Case_ICM")
table(df$is_case)

# df$class = biological methylation class for reference samples

df$case_id <- NA_character_
df$case_id[df$sample == "Case_GBM" | df$sample == "Case_GBM"] <- "Case_GBM"
df$case_id[df$sample == "Case_ICM" | df$sample == "Case_ICM"] <- "Case_ICM"
df$is_case <- !is.na(df$case_id)

# IMPORTANT: keep class clean (do not let cases pollute class labels)
df$class[df$is_case] <- NA

library(dplyr)
library(ggplot2)
library(ggrepel)

family_map <- function(cls) {
  cls <- as.character(cls)
  
  case_when(
    # Embryonal
    cls %in% c("ETMR","MB, WNT","MB, G3","MB, G4","MB, SHH CHL AD","MB, SHH INF",
               "ATRT, MYC","ATRT, SHH","ATRT, TYR","CNS NB, FOXR2","HGNET, BCOR",
               "HGNET, MN1") ~ "Embryonal",
    
    # Glioblastoma / Diffuse glioma
    cls %in% c("DMG, K27","GBM, G34","GBM, MES","GBM, RTK I","GBM, RTK II","GBM, RTK III",
               "GBM, MID","GBM, MYCN","IHG","PXA") ~ "Glioblastoma / Diffuse glioma",
    
    # Glioneuronal + low-grade glial-ish from your missing list
    cls %in% c("CN","DLGNT","LIPN","LGG, DIG/DIA","LGG, DNT","LGG, RGNT","LGG, GG","RETB",
               "LGG, PA PF","LGG, PA MID","LGG, PA/GG ST","LGG, SEGA","LGG, MYB") ~ "Glioneuronal",
    
    # Ependymal
    cls %in% c("EPN, RELA","EPN, PF A","EPN, PF B","EPN, SPINE","EPN, YAP","EPN, MPE",
               "SUBEPN, PF","SUBEPN, SPINE","SUBEPN, ST") ~ "Ependymal",
    
    # Mesenchymal
    cls %in% c("CHORDM","EWS","HMB","MNG","SFT HMPC","EFT, CIC",
               "CPH, ADM","CPH, PAP") ~ "Mesenchymal",
    
    # Plexus
    cls %in% c("PLEX, AD","PLEX, PED A","PLEX, PED B") ~ "Plexus",
    
    # Glioma IDH
    cls %in% c("A IDH","A IDH, HG","O IDH") ~ "Glioma IDH",
    
    # Melanocytic
    cls %in% c("MELAN","MELCYT") ~ "Melanocytic",
    
    # Hematopoietic
    cls %in% c("LYMPHO","PLASMA") ~ "Hematopoietic",
    
    # Controls
    grepl("^CONTR", cls) ~ "Controls",
    
    # Pituitary + pineal / pituicyte-type-ish from your missing list
    cls %in% c("PITUI","PGG, nC","CHGL","ENB, A","ENB, B","ANA PA",
               "PTPR, A","PTPR, B",
               "PIN T,  PB A","PIN T,  PB B","PIN T, PPT",
               "PITAD, ACTH","PITAD, PRL","PITAD, FSH LH","PITAD, TSH",
               "PITAD, STH DNS A","PITAD, STH DNS B","PITAD, STH SPA",
               "SCHW","SCHW, MEL") ~ "Other (Pituitary/Pineal/Peripheral)",
    
    TRUE ~ "Other"
  )
}

df$family <- family_map(df$class)

label_df <- df %>%
  filter(!is_case, !is.na(class)) %>%
  group_by(family) %>%
  summarise(
    x = median(x, na.rm = TRUE),
    y = median(y, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  filter(n >= 25)   # avoid tiny families cluttering

case_colors <- c("Case_GBM"="#D62728", "Case"="#1F77B4")

family_box_colors <- c(
  "Embryonal" = "#D7ECFF",
  "Glioblastoma / Diffuse glioma" = "#D6F5B6",
  "Glioneuronal" = "#F3D6C6",
  "Ependymal" = "#F6C6D7",
  "Mesenchymal" = "#E9D5FF",
  "Plexus" = "#E6E6E6",
  "Glioma IDH" = "#FFF2B3",
  "Melanocytic" = "#DDDDDD",
  "Hematopoietic" = "#EAD1F2",
  "Controls" = "#F0F0F0",
  "Other (Pituitary/Pineal/Peripheral)" = "#FFFFFF",
  "Other" = "#FFFFFF"
)

p_global_onmap <- ggplot(df, aes(x, y)) +
  # reference
  geom_point(
    data = subset(df, !is_case),
    aes(color = class),
    size = 1.2, alpha = 0.85
  ) +
  scale_color_manual(values = class_colors, na.value = "grey80", guide = "none") +
  
  # cases
  geom_point(
    data = subset(df, is_case),
    aes(fill = case_id),
    shape = 21, color = "black", stroke = 1.3, size = 6
  ) +
  scale_fill_manual(values = case_colors, guide = "none") +
  
  geom_text_repel(
    data = subset(df, is_case),
    aes(label = case_id),
    size = 4, max.overlaps = Inf
  ) +
  
  # family callouts (THIS is the “legend on the map”)
  geom_label_repel(
    data = label_df,
    aes(x = x, y = y, label = family, fill = family),
    label.size = 0.25,
    size = 4,
    alpha = 0.92,
    box.padding = 0.6,
    point.padding = 0.2,
    min.segment.length = 0,
    max.overlaps = Inf
  ) +
  scale_fill_manual(values = family_box_colors, guide = "none") +
  
  theme_bw() +
  theme(panel.grid = element_blank()) +
  ggtitle("t-SNE (GSE90496 reference + EPICv2 cases)")

p_global_onmap
ggsave("results/tsne_global_onmap_labels.png", p_global_onmap, width=20, height=8, dpi=600)

label_df <- tibble::tribble(
  ~family, ~x, ~y,
  "Medulloblastoma", -90, -60,
  "Ependymal",       -10,  40,
  "Glioblastoma / Diffuse glioma",  40, -20,
  "Mesenchymal",     -60,  10,
  "Controls",         20,  30
)

# ========================================================

library(dplyr)
library(ggrepel)

class_label_df <- df %>%
  filter(!is_case, !is.na(class)) %>%
  group_by(class) %>%
  summarise(
    x = median(x, na.rm = TRUE),
    y = median(y, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )


p_global_outmap <- ggplot(df, aes(x, y)) +
  # reference
  geom_point(
    data = subset(df, !is_case),
    aes(color = class),
    size = 1.2, alpha = 0.85
  ) +
  scale_color_manual(values = class_colors, na.value = "grey80", guide = "none") +
  
  # cases
  geom_point(
    data = subset(df, is_case),
    aes(fill = case_id),
    shape = 21, color = "black", stroke = 1.3, size = 6
  ) +
  scale_fill_manual(values = case_colors, guide = "none") +
  
  geom_text_repel(
    data = subset(df, is_case),
    aes(label = case_id),
    size = 4, max.overlaps = Inf
  ) +
  
  theme_bw() +
  theme(panel.grid = element_blank()) +
  ggtitle("t-SNE (GSE90496 reference + EPICv2 cases)")



p2 <- p_global_outmap +
  geom_text_repel(
    data = class_label_df,
    aes(x = x, y = y, label = class),
    size = 3,
    color = "black",
    box.padding = 0.25,
    point.padding = 0.15,
    min.segment.length = 0,
    max.overlaps = Inf
  )
p2
ggsave("results/tsne_global_outmap_labels.png", p2, width=20, height=8, dpi=600)

# ===========================================================================

gbm_region <- df %>%
  dplyr::filter(!is.na(family), family == "Glioblastoma / Diffuse glioma")

pad_x <- diff(range(gbm_region$x)) * 0.01; if (pad_x == 0) pad_x <- 1
pad_y <- diff(range(gbm_region$y)) * 0.01; if (pad_y == 0) pad_y <- 1

xlim_gbm <- range(gbm_region$x) + c(-pad_x, pad_x)
ylim_gbm <- range(gbm_region$y) + c(-pad_y, pad_y)

# class centroids (you already built class_label_df; reusing is fine)
class_label_df <- df %>%
  dplyr::filter(!is_case, !is.na(class)) %>%
  dplyr::group_by(class) %>%
  dplyr::summarise(
    x = median(x, na.rm=TRUE),
    y = median(y, na.rm=TRUE),
    n = dplyr::n(),
    .groups="drop"
  )

# keep only labels that fall inside zoom window
class_label_zoom <- class_label_df %>%
  dplyr::filter(
    x >= xlim_gbm[1], x <= xlim_gbm[2],
    y >= ylim_gbm[1], y <= ylim_gbm[2]
  ) %>%
  dplyr::filter(n >= 10)   # tune: 5/10/20 depending clutter

p_zoom_gbm <- ggplot(df, aes(x, y)) +
  # background points (grey, faint)
  geom_point(
    data = dplyr::filter(df, !is_case),
    color = "grey85", size = 1.1, alpha = 0.6
  ) +
  # highlight the points in the zoom window with true colors
  geom_point(
    data = dplyr::filter(df, !is_case,
                         x >= xlim_gbm[1], x <= xlim_gbm[2],
                         y >= ylim_gbm[1], y <= ylim_gbm[2]),
    aes(color = class),
    size = 1.4, alpha = 0.9
  ) +
  scale_color_manual(values = class_colors, na.value = "grey80", guide = "none") +
  
  # cases
  geom_point(
    data = dplyr::filter(df, is_case),
    aes(fill = case_id),
    shape = 21, color = "black", stroke = 1.3, size = 6
  ) +
  scale_fill_manual(values = case_colors, guide = "none") +
  
  ggrepel::geom_text_repel(
    data = dplyr::filter(df, is_case),
    aes(label = case_id),
    size = 4, max.overlaps = Inf
  ) +
  
  # methylation class labels (ONLY inside zoom)
  ggrepel::geom_text_repel(
    data = class_label_zoom,
    aes(label = class),
    size = 3,
    box.padding = 0.25,
    point.padding = 0.15,
    min.segment.length = 0,
    max.overlaps = Inf
  ) +
  
  coord_cartesian(xlim = xlim_gbm, ylim = ylim_gbm) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  ggtitle("Zoom: GBM / diffuse glioma region")

p_zoom_gbm
ggsave("results/tsne_zoom_labels.png", p_zoom_gbm, width=20, height=8, dpi=600)

# ============================================================================



