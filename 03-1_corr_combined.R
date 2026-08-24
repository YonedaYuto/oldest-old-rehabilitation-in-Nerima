library(data.table)
library(here)
library(ggplot2)

need_pkg <- function(p) if (!requireNamespace(p, quietly = TRUE))
  stop(sprintf("パッケージ '%s' が必要です。install.packages('%s')", p, p))
need_pkg("ggcorrplot")
need_pkg("patchwork")

LABELS_PATH <- here::here("00_labels.R")
if (!exists("relabel_vars")) {
  if (!file.exists(LABELS_PATH)) stop("00_labels.R が見つかりません: ", LABELS_PATH)
  source(LABELS_PATH)
}

OUT_DIR <- here::here("data")
FIG_DIR <- here::here("figures")
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

if (!exists("centered_pipelines"))
  centered_pipelines <- readRDS(file.path(OUT_DIR, "bnb_pipelines_centered.rds"))
if (!exists("selected_meta"))
  selected_meta <- readRDS(file.path(OUT_DIR, "bnb_selected_meta.rds"))

SYS_LABELS <- c(A = "System A (worst-value)", B = "System B (missing-category)")
sys_of  <- function(nm) sub(".*_", "", nm)
subset_of <- function(nm) sub("_.*", "", nm)

build_numeric_mat <- function(dt, meta, use_centered = TRUE, include_outcome = TRUE) {
  cols <- list()
  for (v in intersect(meta$continuous, names(dt))) {
    nm <- if (use_centered && paste0(v, "_c") %in% names(dt)) paste0(v, "_c") else v
    cols[[v]] <- as.numeric(dt[[nm]])
  }
  for (v in intersect(meta$categorical, names(dt))) {
    x <- dt[[v]]
    if (is.numeric(x) || is.integer(x)) {
      cols[[v]] <- as.numeric(x)
    } else {
      f <- as.factor(x)
      if (nlevels(f) <= 2L) {
        cols[[v]] <- as.numeric(f) - 1
      } else {
        for (lv in levels(f)[-1]) cols[[paste0(v, ":", lv)]] <- as.numeric(f == lv)
      }
    }
  }
  if (include_outcome && "good" %in% names(dt)) cols[["good"]] <- as.numeric(dt[["good"]])
  as.matrix(as.data.frame(cols, check.names = FALSE))
}

make_panels <- function(dt, meta, nm) {
  sys_lab <- SYS_LABELS[[sys_of(nm)]]

  contv <- intersect(meta$continuous, names(dt))
  matc  <- sapply(contv, function(v) {
    cn <- if (paste0(v, "_c") %in% names(dt)) paste0(v, "_c") else v
    as.numeric(dt[[cn]])
  })
  colnames(matc) <- contv
  cor_p <- cor(matc, use = "pairwise.complete.obs", method = "pearson")
  dimnames(cor_p) <- list(relabel_vars(rownames(cor_p)), relabel_vars(colnames(cor_p)))

  matm  <- build_numeric_mat(dt, meta, use_centered = TRUE, include_outcome = TRUE)
  cor_s <- cor(matm, use = "pairwise.complete.obs", method = "spearman")
  dimnames(cor_s) <- list(relabel_vars(rownames(cor_s)), relabel_vars(colnames(cor_s)))

  pear <- suppressWarnings(ggcorrplot::ggcorrplot(
    cor_p, type = "lower", lab = TRUE, lab_size = 2.2,
    colors = c("#2166AC", "white", "#B2182B"), tl.cex = 7.5
  )) + ggtitle(sprintf("%s\nPearson (continuous)", sys_lab)) +
    theme(plot.title = element_text(size = 10, face = "bold"))

  spear <- suppressWarnings(ggcorrplot::ggcorrplot(
    cor_s, type = "lower", lab = FALSE,
    colors = c("#2166AC", "white", "#B2182B"), tl.cex = 6.5
  )) + ggtitle(sprintf("%s\nSpearman (all variables)", sys_lab)) +
    theme(plot.title = element_text(size = 10, face = "bold"))

  list(pearson = pear, spearman = spear)
}

build_group_canvas <- function(members, subset_title, outfile) {
  pa <- make_panels(centered_pipelines[[members[1]]], selected_meta[[members[1]]], members[1])
  pb <- make_panels(centered_pipelines[[members[2]]], selected_meta[[members[2]]], members[2])

  canvas <- patchwork::wrap_plots(
    pa$pearson, pa$spearman,
    pb$pearson, pb$spearman,
    ncol = 2, nrow = 2
  ) + patchwork::plot_annotation(
        title = subset_title,
        theme = ggplot2::theme(plot.title = element_text(size = 14, face = "bold"))
      )

  ggsave(file.path(FIG_DIR, outfile), canvas, width = 13, height = 11, dpi = 150)
  message("保存: ", file.path(FIG_DIR, outfile))
  invisible(canvas)
}

build_group_canvas(
  members      = c("severe_A", "severe_B"),
  subset_title = "Severe subset (mFIM at Adm <= 26) — Correlation matrices",
  outfile      = "fig1_corr_severe_combined.png"
)

build_group_canvas(
  members      = c("elderly_A", "elderly_B"),
  subset_title = "Elderly subset (age >= 90) — Correlation matrices",
  outfile      = "fig1_corr_elderly_combined.png"
)

message("§3-1 完了。重症・高齢それぞれ A/B をまとめた相関図を ",
        normalizePath(FIG_DIR), " に保存しました。")
