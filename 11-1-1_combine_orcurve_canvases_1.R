library(data.table)
library(here)
library(ggplot2)

need_pkg <- function(p) if (!requireNamespace(p, quietly = TRUE))
  stop(sprintf("パッケージ '%s' が必要です。install.packages('%s')", p, p))
need_pkg("patchwork")

LABELS_PATH <- here::here("00_labels.R")
if (!exists("relabel_pipeline")) {
  if (!file.exists(LABELS_PATH)) stop("00_labels.R が見つかりません: ", LABELS_PATH)
  source(LABELS_PATH)
}

OUT_DIR <- here::here("data")
FIG_DIR <- here::here("figures")
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

BASE_SIZE   <- 15
STRIP_SIZE  <- 18
AXIS_TITLE  <- 15
AXIS_TEXT   <- 13
PANEL_TITLE <- 18
GROUP_TITLE <- 22

FACET_NCOL  <- 2
CURVE_COL   <- "#2166AC"
RIBBON_COL  <- "#2166AC"

GROUPS <- list(
  list(title = paste0("Adjusted odds-ratio curves for continuous admission prognostic factors ",
                      "in the oldest-old subset"),
       subtitle = paste0("Age 90 or more. Main-effect model; odds ratios are relative to the median ",
                         "of each factor (reference OR = 1) and bands are 95% confidence bands.\n",
                         "Left, primary analysis (worst-value system); right, sensitivity analysis ",
                         "(missing-category system)."),
       pipes = c("elderly_A", "elderly_B"),
       out   = "fig4_orcurve_provisional_elderly_combined.png"),
  list(title = paste0("Adjusted odds-ratio curves for continuous admission prognostic factors ",
                      "in the severe subset"),
       subtitle = paste0("Admission motor FIM 26 or less. Main-effect model; odds ratios are relative ",
                         "to the median of each factor (reference OR = 1) and bands are 95% confidence ",
                         "bands.\nLeft, primary analysis (worst-value system); right, sensitivity ",
                         "analysis (missing-category system)."),
       pipes = c("severe_A", "severe_B"),
       out   = "fig4_orcurve_provisional_severe_combined.png")
)
read_safe <- function(f) if (file.exists(f)) readRDS(f) else NULL

plain_num <- function(x) format(x, scientific = FALSE, trim = TRUE,
                                drop0trailing = TRUE, big.mark = "")

curve_path <- file.path(OUT_DIR, "bnb_orcurve_provisional.rds")
if (!file.exists(curve_path))
  stop("曲線データが見つかりません: ", curve_path, "（先に §11-1 を実行してください）")
ALL <- as.data.table(readRDS(curve_path))

build_plot <- function(nm) {
  d <- ALL[pipeline == nm & is.finite(OR) & is.finite(lcl) & is.finite(ucl)]
  if (nrow(d) == 0) {
    message(sprintf("  [%s] OR カーブのデータが無い（連続変数なし or §11-1 未実行）。", nm))
    return(NULL)
  }
  nvar  <- length(unique(d$variable_lab))
  ncolf <- min(FACET_NCOL, nvar)

  prov <- read_safe(file.path(OUT_DIR, sprintf("bnb_provisional_%s.rds", nm)))
  med_df <- NULL
  if (!is.null(prov) && !is.null(prov$centers)) {
    med_df <- unique(d[, .(variable, variable_lab)])
    med_df[, med := vapply(variable, function(v) {
      z <- prov$centers[[v]]; if (is.null(z)) NA_real_ else as.numeric(z) }, numeric(1))]
    med_df <- med_df[is.finite(med)]
  }

  ggplot(d, aes(x = x_raw)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    { if (!is.null(med_df) && nrow(med_df) > 0)
        geom_vline(data = med_df, aes(xintercept = med),
                   linetype = "longdash", colour = "grey55", linewidth = 0.4)
      else NULL } +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), fill = RIBBON_COL, alpha = 0.16) +
    geom_line(aes(y = OR), colour = CURVE_COL, linewidth = 1.0) +
    facet_wrap(~ variable_lab, scales = "free", ncol = ncolf) +
    scale_y_log10(labels = plain_num) +
    labs(title = relabel_pipeline(nm),
         x = "Predictor value (original scale)",
         y = "Adjusted odds ratio (log scale)") +
    theme_bw(base_size = BASE_SIZE) +
    theme(plot.title   = element_text(size = PANEL_TITLE, face = "bold"),
          strip.text   = element_text(size = STRIP_SIZE, face = "bold"),
          axis.title   = element_text(size = AXIS_TITLE),
          axis.text    = element_text(size = AXIS_TEXT),
          panel.grid.minor = element_blank(),
          panel.spacing = grid::unit(1.0, "lines"))
}

combine_group <- function(grp) {
  message(sprintf("== %s ==", grp$title))
  built <- lapply(grp$pipes, build_plot)
  names(built) <- grp$pipes
  built <- Filter(Negate(is.null), built)
  if (length(built) == 0) { message("  → 描画対象が無いためスキップ。"); return(invisible(NULL)) }

  combo <- patchwork::wrap_plots(built, ncol = length(built)) +
    patchwork::plot_annotation(
      title    = grp$title,
      subtitle = grp$subtitle,
      theme = theme(plot.title    = element_text(size = GROUP_TITLE, face = "bold"),
                    plot.subtitle = element_text(size = GROUP_TITLE - 3, colour = "grey30")))

  nvar_each  <- vapply(grp$pipes, function(nm) length(unique(ALL[pipeline == nm]$variable_lab)), integer(1))
  ncolf_each <- pmin(FACET_NCOL, pmax(nvar_each, 1L))
  nrow_each  <- ceiling(nvar_each / ncolf_each)
  nrow_each[!is.finite(nrow_each) | nrow_each < 1] <- 1L
  facet_cols_total <- sum(ncolf_each[nvar_each > 0])
  per_panel_w <- 4.6; per_panel_h <- 3.9
  width  <- max(facet_cols_total, 1L) * per_panel_w + 0.8
  height <- max(nrow_each) * per_panel_h + 1.4

  outfile <- file.path(FIG_DIR, grp$out)
  ggsave(outfile, combo, width = width, height = height, dpi = 150, limitsize = FALSE)
  message(sprintf("  → 保存: %s（%.1f x %.1f in, パネル %d 列）",
                  outfile, width, height, length(built)))
  invisible(outfile)
}

invisible(lapply(GROUPS, combine_group))

message("\n§11-1-1 完了。重症(A|B)・高齢(A|B)の OR カーブを各 1 枚に結合し、",
        normalizePath(FIG_DIR), " に保存しました。",
        " 変数名（ファセット見出し）はフォントを拡大して読みやすくしています。")
