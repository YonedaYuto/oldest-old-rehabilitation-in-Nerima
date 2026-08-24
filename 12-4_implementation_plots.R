library(data.table)
library(here)
library(ggplot2)

LABELS_PATH <- here::here("00_labels.R")
if (!exists("relabel_vars")) {
  if (!file.exists(LABELS_PATH)) stop("00_labels.R が見つかりません: ", LABELS_PATH)
  source(LABELS_PATH)
}

OUT_DIR <- here::here("data")
FIG_DIR <- here::here("figures")
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

PIPELINES <- names(PIPE_LABELS)

IMPL_MODEL <- c(severe_A = "Final", severe_B = "Final",
                elderly_A = "Provisional", elderly_B = "Provisional")

BASE_SIZE   <- 18
STRIP_SIZE  <- 17
AXIS_T_SIZE <- 18
AXIS_X_SIZE <- 14
ANNOT_SIZE  <- 7
CAL_BINS <- 10L; N_GRID <- 100L; CI_LEVEL <- 0.95

CAL_COL <- "#2166AC"; ROC_COL <- "#B2182B"; REF_COL <- "grey45"; RIBBON_A <- 0.15

read_safe <- function(f) if (file.exists(f)) readRDS(f) else NULL

calib_curve_data <- function(y, p, bins = CAL_BINS, n_grid = N_GRID) {
  d <- data.table(y = y, p = pmin(pmax(p, 0), 1))
  lo <- tryCatch(stats::loess(y ~ p, data = d, span = 0.9, degree = 1, family = "gaussian",
                              control = stats::loess.control(surface = "direct")),
                 error = function(e) NULL)
  grid <- seq(min(d$p), max(d$p), length.out = n_grid)
  curve <- if (!is.null(lo)) {
    pr <- stats::predict(lo, newdata = data.frame(p = grid), se = TRUE)
    z  <- stats::qnorm(1 - (1 - CI_LEVEL) / 2)
    data.table(p = grid, obs = pmin(pmax(pr$fit, 0), 1),
               lcl = pmin(pmax(pr$fit - z * pr$se.fit, 0), 1),
               ucl = pmin(pmax(pr$fit + z * pr$se.fit, 0), 1))
  } else data.table()
  br <- unique(stats::quantile(d$p, probs = seq(0, 1, length.out = bins + 1), na.rm = TRUE))
  binpts <- data.table()
  if (length(br) >= 3) {
    d[, grp := cut(p, breaks = br, include.lowest = TRUE)]
    binpts <- d[, {
      n <- .N; k <- sum(y); ph <- mean(p); oh <- k / n
      zc <- stats::qnorm(1 - (1 - CI_LEVEL) / 2); den <- 1 + zc^2 / n
      ctr <- (oh + zc^2 / (2 * n)) / den
      hw  <- zc * sqrt(oh * (1 - oh) / n + zc^2 / (4 * n^2)) / den
      .(mean_pred = ph, obs_prop = oh, lcl = max(ctr - hw, 0), ucl = min(ctr + hw, 1), n = n)
    }, by = grp]
  }
  list(curve = curve, bins = binpts)
}
roc_curve_data <- function(y, p) {
  thr <- sort(unique(c(-Inf, p, Inf)), decreasing = TRUE)
  np <- sum(y == 1); nn <- sum(y == 0)
  tpr <- vapply(thr, function(t) sum(p >  t & y == 1) / np, numeric(1))
  fpr <- vapply(thr, function(t) sum(p >  t & y == 0) / nn, numeric(1))
  data.table(fpr = fpr, tpr = tpr)[order(fpr, tpr)]
}

path_tab <- read_safe(file.path(OUT_DIR, "bnb_path_summary.rds"))
if (!is.null(path_tab)) { path_tab <- as.data.table(path_tab)
  path_tab[, pipeline := as.character(pipeline)]; path_tab[, model := as.character(model)] }

get_cal_roc <- function(nm, impl) {
  po <- read_safe(file.path(OUT_DIR, sprintf("bnb_performance_%s.rds", nm)))
  if (is.null(po)) { message(sprintf("[%s] §12 出力なし。スキップ。", nm)); return(NULL) }

  auc <- NA_real_
  if (!is.null(po$performance)) {
    pf <- as.data.table(po$performance); pf[, model := as.character(model)]
    a <- pf[model == impl, AUC]; if (length(a)) auc <- a[1]
  }
  est <- NA_character_
  if (!is.null(path_tab)) { pt <- path_tab[pipeline == nm & model == impl]; if (nrow(pt)) est <- pt$estimator[1] }

  curve <- bins <- data.table(); roc <- data.table()
  cal <- if (!is.null(po$calibration)) as.data.table(po$calibration) else data.table()
  if (nrow(cal)) {
    cal[, model := as.character(model)]
    if (!"kind" %in% names(cal)) cal[, kind := NA_character_]
    curve <- cal[model == impl & kind == "curve",
                 intersect(c("p", "obs", "lcl", "ucl"), names(cal)), with = FALSE]
    bins  <- cal[model == impl & is.na(kind),
                 intersect(c("mean_pred", "obs_prop", "lcl", "ucl", "n"), names(cal)), with = FALSE]
  }
  rc <- if (!is.null(po$roc)) as.data.table(po$roc) else data.table()
  if (nrow(rc)) { rc[, model := as.character(model)]; roc <- rc[model == impl, .(fpr, tpr)] }

  if ((nrow(curve) == 0 && nrow(bins) == 0) || nrow(roc) == 0) {
    pp <- if (!is.null(po$pooled_pred)) as.data.table(po$pooled_pred) else data.table()
    if (nrow(pp)) {
      pp[, model := as.character(model)]; sub <- pp[model == impl]
      if (nrow(sub)) {
        y <- as.integer(sub$y); p <- sub$p
        if (nrow(curve) == 0 && nrow(bins) == 0) {
          cc <- calib_curve_data(y, p); curve <- cc$curve; bins <- cc$bins
        }
        if (nrow(roc) == 0) roc <- roc_curve_data(y, p)
        if (!is.finite(auc) && length(unique(y)) == 2) {
          pos <- p[y == 1]; neg <- p[y == 0]; rr <- rank(c(pos, neg))
          auc <- (sum(rr[seq_along(pos)]) - length(pos) * (length(pos) + 1) / 2) /
                 (length(pos) * length(neg))
        }
      }
    }
  }
  list(curve = curve, bins = bins, roc = roc, auc = auc, model = impl, estimator = est)
}

curve_all <- list(); bins_all <- list(); roc_all <- list(); ann_all <- list()
for (nm in PIPELINES) {
  impl <- IMPL_MODEL[[nm]]
  if (is.null(impl) || is.na(impl)) { message(sprintf("[%s] 実装モデル未定義。スキップ。", nm)); next }
  cr <- get_cal_roc(nm, impl); if (is.null(cr)) next
  message(sprintf("[%s] 実装モデル=%s, AUC=%.3f%s", nm, impl, cr$auc,
                  if (!is.na(cr$estimator)) sprintf(" (%s)", cr$estimator) else ""))
  if (nrow(cr$curve)) curve_all[[nm]] <- copy(cr$curve)[, pipeline := nm]
  if (nrow(cr$bins))  bins_all[[nm]]  <- copy(cr$bins)[, pipeline := nm]
  if (nrow(cr$roc))   roc_all[[nm]]   <- copy(cr$roc)[, pipeline := nm]
  ann_all[[nm]] <- data.table(pipeline = nm, auc = cr$auc, model = impl, estimator = cr$estimator)
}
if (length(ann_all) == 0)
  stop("較正/ROC を組み立てられませんでした。先に §12 を実行してください。")

plabs <- vapply(PIPELINES, relabel_pipeline, character(1))
mk_fac <- function(dt) { dt[, pipeline := factor(pipeline, levels = PIPELINES, labels = plabs)]; dt }
curve_dt <- if (length(curve_all)) mk_fac(rbindlist(curve_all, fill = TRUE)) else data.table()
bins_dt  <- if (length(bins_all))  mk_fac(rbindlist(bins_all,  fill = TRUE)) else data.table()
roc_dt   <- if (length(roc_all))   mk_fac(rbindlist(roc_all,   fill = TRUE)) else data.table()
ann_dt   <- mk_fac(rbindlist(ann_all, fill = TRUE))
ann_dt[, label := ifelse(is.finite(auc), sprintf("AUC %.3f", auc), "AUC —")]

big_theme <- theme_bw(base_size = BASE_SIZE) +
  theme(legend.position = "none",
        strip.text = element_text(size = STRIP_SIZE, face = "bold"),
        strip.background = element_rect(fill = "grey92", colour = NA),
        axis.title = element_text(size = AXIS_T_SIZE),
        axis.text  = element_text(size = AXIS_X_SIZE),
        plot.title = element_text(size = BASE_SIZE + 2, face = "bold"),
        panel.grid.minor = element_blank())

gcal <- ggplot() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = REF_COL, linewidth = 0.7) +
  { if (nrow(curve_dt)) geom_ribbon(data = curve_dt, aes(x = p, ymin = lcl, ymax = ucl),
                                    fill = CAL_COL, alpha = RIBBON_A) } +
  { if (nrow(curve_dt)) geom_line(data = curve_dt, aes(x = p, y = obs),
                                  colour = CAL_COL, linewidth = 1.2) } +
  { if (nrow(bins_dt))  geom_point(data = bins_dt, aes(x = mean_pred, y = obs_prop),
                                   colour = CAL_COL, size = 2.8) } +
  { if (nrow(bins_dt))  geom_errorbar(data = bins_dt, aes(x = mean_pred, ymin = lcl, ymax = ucl),
                                      colour = CAL_COL, width = 0.02, linewidth = 0.6) } +
  facet_wrap(~ pipeline, ncol = 2) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(x = "Predicted probability", y = "Observed proportion") +
  big_theme
ggsave(file.path(FIG_DIR, "fig12-4_calibration_4panel.png"),
       gcal, width = 10, height = 10.4, dpi = 200)
message("\n較正プロット（4パネル）を保存: ",
        file.path(FIG_DIR, "fig12-4_calibration_4panel.png"))

groc <- ggplot() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = REF_COL, linewidth = 0.7) +
  { if (nrow(roc_dt)) geom_line(data = roc_dt, aes(x = fpr, y = tpr),
                                colour = ROC_COL, linewidth = 1.3) } +
  geom_text(data = ann_dt, aes(x = 0.62, y = 0.08, label = label),
            size = ANNOT_SIZE, hjust = 0, colour = "grey15") +
  facet_wrap(~ pipeline, ncol = 2) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(x = "1 - Specificity", y = "Sensitivity") +
  big_theme
ggsave(file.path(FIG_DIR, "fig12-4_roc_4panel.png"),
       groc, width = 10, height = 10.4, dpi = 200)
message("ROC（4パネル）を保存: ", file.path(FIG_DIR, "fig12-4_roc_4panel.png"))

message("\n§12-4 完了。")
message("  実装モデル: ", paste(sprintf("%s=%s", names(IMPL_MODEL), IMPL_MODEL), collapse = ", "))
message("  ・較正 4 パネル : figures/fig12-4_calibration_4panel.png")
message("  ・ROC  4 パネル : figures/fig12-4_roc_4panel.png")
message("注記: 較正/ROC は §12 が保存した MI プール予測確率ベースの素材を流用（§12 図5/図6 と同一）。")
message("      フォントは BASE_SIZE 等で調整可能。凡例なし・各パネル見出し = パイプライン名。")
