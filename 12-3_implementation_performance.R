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

DIGITS_AUC <- 3L
DIGITS_CAL <- 2L
DIGITS_GEN <- 3L

read_safe <- function(f) if (file.exists(f)) readRDS(f) else NULL
pipe_disp <- function(x) vapply(as.character(x), relabel_pipeline, character(1))

subset_disp <- function(nm) {
  if (grepl("^severe", nm)) "Severe (mFIM_in <= 26)"
  else if (grepl("^elderly", nm)) "Elderly (age >= 90)"
  else nm
}
system_disp <- function(nm) {
  s <- sub(".*_", "", nm)
  if (s == "A") "A (worst-value)" else if (s == "B") "B (missing-category)" else s
}

fmt_ci <- function(est, lcl, ucl, d) {
  ifelse(is.finite(est),
         ifelse(is.finite(lcl) & is.finite(ucl),
                sprintf(paste0("%.", d, "f [%.", d, "f, %.", d, "f]"), est, lcl, ucl),
                sprintf(paste0("%.", d, "f"), est)),
         "—")
}
fmt_num <- function(x, d) ifelse(is.finite(x), sprintf(paste0("%.", d, "f"), x), "—")

perf <- read_safe(file.path(OUT_DIR, "bnb_performance_summary.rds"))
if (is.null(perf)) {
  message("bnb_performance_summary.rds が無いため per-pipeline から再構成します。")
  rows <- lapply(PIPELINES, function(nm) {
    p <- read_safe(file.path(OUT_DIR, sprintf("bnb_performance_%s.rds", nm)))
    if (is.null(p)) NULL else p$performance
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0)
    stop("§12 の性能出力（bnb_performance_*.rds）が見つかりません。先に §12 を実行してください。")
  perf <- rbindlist(rows, fill = TRUE)
}
perf <- as.data.table(perf)
perf[, pipeline := as.character(pipeline)]
perf[, model    := as.character(model)]

gvif_tab <- read_safe(file.path(OUT_DIR, "bnb_gvif_table.rds"))
if (!is.null(gvif_tab)) { gvif_tab <- as.data.table(gvif_tab)
  gvif_tab[, pipeline := as.character(pipeline)]; gvif_tab[, model := as.character(model)] }
path_tab <- read_safe(file.path(OUT_DIR, "bnb_path_summary.rds"))
if (!is.null(path_tab)) { path_tab <- as.data.table(path_tab)
  path_tab[, pipeline := as.character(pipeline)]; path_tab[, model := as.character(model)] }
stab_tab <- read_safe(file.path(OUT_DIR, "bnb_var_stability.rds"))
if (!is.null(stab_tab)) { stab_tab <- as.data.table(stab_tab)
  stab_tab[, pipeline := as.character(pipeline)]; stab_tab[, model := as.character(model)] }

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a[1])) b else a
get_df_events <- function(nm, impl) {
  f <- if (impl == "Final") sprintf("bnb_final_%s.rds", nm) else sprintf("bnb_provisional_%s.rds", nm)
  o <- read_safe(file.path(OUT_DIR, f))
  if (is.null(o) || is.null(o$selection)) return(list(df = NA_real_, events = NA_real_))
  list(df = o$selection$final_df %||% NA_real_, events = o$selection$n_events %||% NA_real_)
}

rows <- list()
for (nm in PIPELINES) {
  impl <- IMPL_MODEL[[nm]]
  if (is.null(impl) || is.na(impl)) { message(sprintf("[%s] 実装モデル未定義。スキップ。", nm)); next }
  pr <- perf[pipeline == nm & model == impl]
  if (!nrow(pr)) { message(sprintf("[%s] §12 に %s モデルの性能が無い。スキップ。", nm, impl)); next }
  pr <- pr[1]

  est <- NA_character_; cim <- NA_character_; pth <- pr$path
  if (!is.null(path_tab)) {
    pt <- path_tab[pipeline == nm & model == impl]
    if (nrow(pt)) { pth <- pt$path[1]; est <- pt$estimator[1]; cim <- pt$ci_method[1] }
  }
  if (is.na(est)) est <- if (identical(pth, "firth")) "logistf (Firth)" else "glm"
  if (is.na(cim)) cim <- if (identical(pth, "firth")) "CLIP profile (MI-pooled)" else "Wald/Rubin"

  gvadj <- pr$max_GVIF_adj %||% NA_real_; gvar <- NA_character_
  if (!is.null(gvif_tab)) {
    gt <- gvif_tab[pipeline == nm & model == impl & is.finite(GVIF_adj)]
    if (nrow(gt)) { i <- which.max(gt$GVIF_adj); gvadj <- gt$GVIF_adj[i]; gvar <- gt$variable[i] }
  }

  de  <- get_df_events(nm, impl)
  ev  <- if (is.finite(pr$events)) pr$events else de$events
  epv <- if (is.finite(ev) && is.finite(de$df) && de$df > 0) ev / de$df else NA_real_

  nfrag <- NA_integer_
  if (!is.null(stab_tab)) {
    st <- stab_tab[pipeline == nm & model == impl]
    if (nrow(st)) nfrag <- sum(st$verdict == "fragile", na.rm = TRUE)
  }

  rows[[nm]] <- data.table(
    pipeline = nm, pipeline_disp = relabel_pipeline(nm),
    subset = subset_disp(nm), system = system_disp(nm),
    impl_model = impl, path = pth, estimator = est, ci_method = cim,
    n = pr$n, events = ev, df = de$df, EPV = epv,
    AUC = pr$AUC, AUC_lcl = pr$AUC_lcl, AUC_ucl = pr$AUC_ucl,
    AUC_corrected = pr$AUC_corrected,
    cal_slope = pr$cal_slope, cal_slope_lcl = pr$cal_slope_lcl, cal_slope_ucl = pr$cal_slope_ucl,
    cal_slope_corrected = pr$cal_slope_corrected, cal_intercept = pr$cal_intercept,
    Brier = pr$Brier, Nagelkerke_R2 = pr$Nagelkerke_R2,
    HL_stat = pr$HL_stat, HL_df = pr$HL_df, HL_p = pr$HL_p,
    GVIF_adj_max = gvadj, GVIF_adj_var = gvar, n_fragile = nfrag)
}
if (length(rows) == 0) stop("実装モデルの性能を組み立てられませんでした。§12 を確認してください。")

impl_dt <- rbindlist(rows, fill = TRUE)
impl_dt[, pipeline := factor(pipeline, levels = PIPELINES)]
setorder(impl_dt, pipeline)

saveRDS(impl_dt, file.path(OUT_DIR, "bnb_implementation_performance.rds"))
fwrite(impl_dt, file.path(OUT_DIR, "bnb_implementation_performance.csv"))

fmt_dt <- impl_dt[, .(
  Pipeline = pipeline_disp,
  Subset = subset, System = system,
  `AUC (95% CI)` = fmt_ci(AUC, AUC_lcl, AUC_ucl, DIGITS_AUC),
  `Optimism-corrected AUC` = fmt_num(AUC_corrected, DIGITS_AUC),
  `Calibration slope (95% CI)` = fmt_ci(cal_slope, cal_slope_lcl, cal_slope_ucl, DIGITS_CAL),
  `Corrected slope` = fmt_num(cal_slope_corrected, DIGITS_CAL),
  `Calibration intercept` = fmt_num(cal_intercept, DIGITS_CAL),
  `Hosmer-Lemeshow` = ifelse(is.finite(HL_p),
                             sprintf("chi2=%.1f, df=%s, p=%.3f", HL_stat,
                                     ifelse(is.finite(HL_df), as.character(as.integer(HL_df)), "—"), HL_p),
                             "—"),
  `Max GVIF^(1/2Df)` = ifelse(is.finite(GVIF_adj_max),
                              ifelse(!is.na(GVIF_adj_var),
                                     sprintf("%.2f (%s)", GVIF_adj_max, GVIF_adj_var),
                                     sprintf("%.2f", GVIF_adj_max)),
                              "—")
)]
fwrite(fmt_dt, file.path(OUT_DIR, "bnb_implementation_performance_formatted.csv"))

message("\n==== §12-3 実装モデル 性能サマリ ====")
message("（実装モデル: 重症=最終 Final / 高齢=暫定 Provisional）\n")
print(fmt_dt, row.names = FALSE)

if (requireNamespace("gridExtra", quietly = TRUE) &&
    requireNamespace("grid", quietly = TRUE)) {
  metric_cols <- setdiff(names(fmt_dt), c("Pipeline", "Subset", "System"))
  tb <- data.table(Metric = metric_cols)
  for (i in seq_len(nrow(fmt_dt))) {
    colname <- as.character(fmt_dt$Pipeline[i])
    tb[[colname]] <- unlist(fmt_dt[i, ..metric_cols], use.names = FALSE)
  }
  th <- gridExtra::ttheme_minimal(
    core = list(fg_params = list(hjust = 0, x = 0.02, fontsize = 8),
                bg_params = list(fill = c("grey97", "white"))),
    colhead = list(fg_params = list(fontsize = 8.5, fontface = "bold"),
                   bg_params = list(fill = "grey90")),
    rowhead = list(fg_params = list(fontsize = 8)))
  g <- gridExtra::tableGrob(tb, rows = NULL, theme = th)
  ncol_tab <- ncol(tb); nrow_tab <- nrow(tb)
  ggsave(file.path(FIG_DIR, "fig_implementation_performance.png"),
         g, width = 2.4 + 2.0 * (ncol_tab - 1), height = 0.34 * nrow_tab + 0.8,
         dpi = 200, limitsize = FALSE)
  message("\n報告用の表画像を保存: ", file.path(FIG_DIR, "fig_implementation_performance.png"))
} else {
  message("\n（gridExtra/grid が無いため表画像はスキップ。CSV をご利用ください。")
  message("  install.packages('gridExtra') で表画像も生成されます。）")
}

message("\n§12-3 完了。")
message("  ・数値テーブル : data/bnb_implementation_performance.rds, .csv")
message("  ・整形テーブル : data/bnb_implementation_performance_formatted.csv")
message("  ・表画像       : figures/fig_implementation_performance.png（gridExtra があれば）")
message("注記: 本表は §12（性能）/§12-1（GVIF）/§12-2（経路・安定性）の保存値を流用して組み立てたもので、")
message("      再計算は行っていない（値を変えるには該当スクリプトを再実行）。AUC/較正勾配の補正は Harrell 流")
message("      optimism 補正（§12）、Hosmer–Lemeshow は分位群 D2 プール（calibration plot を主・H&L を補助）、")
message("      GVIF^(1/(2Df)) は自由度調整済み GVIF（§12-1, Fox–Monette）。EPV は events / final_df（§8/§10）。")
