library(data.table)
library(here)

need_pkg <- function(p) if (!requireNamespace(p, quietly = TRUE))
  stop(sprintf("パッケージ '%s' が必要です。install.packages('%s')", p, p))

LABELS_PATH <- here::here("00_labels.R")
if (!exists("relabel_vars")) {
  if (!file.exists(LABELS_PATH)) stop("00_labels.R が見つかりません: ", LABELS_PATH)
  source(LABELS_PATH)
}

OUT_DIR   <- here::here("data")
PIPELINES <- names(PIPE_LABELS)
CI_LEVEL  <- 0.95
DIGITS_B  <- 3L
DIGITS_OR <- 2L

MANUSCRIPT_TERMS <- list(
  severe_A  = NULL,
  severe_B  = NULL,
  elderly_A = NULL,
  elderly_B = NULL
)

`%||%` <- function(a, b) if (is.null(a)) b else a

fmt_ci <- function(est, lo, hi, digits = DIGITS_OR) {
  sprintf("%.*f (%.*f to %.*f)", digits, est, digits, lo, digits, hi)
}

fmt_p <- function(p) {
  ifelse(is.na(p), "—",
         ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

pick_col <- function(dt, candidates) {
  hit <- candidates[candidates %in% names(dt)]
  if (!length(hit)) return(NULL)
  hit[1]
}

extract_coef_table <- function(fin, pipeline) {
  cand <- c("coef_table", "pooled", "pooled_table", "estimates", "coefs", "summary")
  tab <- NULL
  for (nm in cand) if (!is.null(fin[[nm]])) { tab <- fin[[nm]]; break }
  if (is.null(tab))
    stop(sprintf("[%s] 係数表が見つかりません。bnb_final_%s.rds の構造を確認してください: %s",
                 pipeline, pipeline, paste(names(fin), collapse = ", ")))
  tab <- as.data.table(tab)

  c_term <- pick_col(tab, c("term", "variable", "coef", "name", "rowname"))
  c_est  <- pick_col(tab, c("estimate", "beta", "coef_est", "logOR", "est"))
  c_se   <- pick_col(tab, c("std.error", "se", "SE", "std_error"))
  c_lo   <- pick_col(tab, c("conf.low", "lo", "lower", "ci_lo", "CI_low"))
  c_hi   <- pick_col(tab, c("conf.high", "hi", "upper", "ci_hi", "CI_high"))
  c_p    <- pick_col(tab, c("p.value", "p", "pval", "p_value"))
  if (is.null(c_term) || is.null(c_est))
    stop(sprintf("[%s] term / estimate 列を特定できません: %s",
                 pipeline, paste(names(tab), collapse = ", ")))

  out <- data.table(
    pipeline = pipeline,
    term     = as.character(tab[[c_term]]),
    beta     = as.numeric(tab[[c_est]]),
    se       = if (is.null(c_se)) NA_real_ else as.numeric(tab[[c_se]]),
    lo_beta  = if (is.null(c_lo)) NA_real_ else as.numeric(tab[[c_lo]]),
    hi_beta  = if (is.null(c_hi)) NA_real_ else as.numeric(tab[[c_hi]]),
    p        = if (is.null(c_p))  NA_real_ else as.numeric(tab[[c_p]])
  )

  if (all(is.finite(out$lo_beta)) && all(out$lo_beta > 0) &&
      mean(abs(log(out$lo_beta) - out$beta), na.rm = TRUE) <
      mean(abs(out$lo_beta      - out$beta), na.rm = TRUE)) {
    out[, `:=`(lo_beta = log(lo_beta), hi_beta = log(hi_beta))]
  }
  z <- stats::qnorm(1 - (1 - CI_LEVEL) / 2)
  out[is.na(lo_beta) & is.finite(se), `:=`(lo_beta = beta - z * se,
                                           hi_beta = beta + z * se)]
  out[, `:=`(OR = exp(beta), OR_lo = exp(lo_beta), OR_hi = exp(hi_beta))]
  out[]
}

order_terms <- function(dt) {
  is_int   <- grepl("^\\(?Intercept\\)?$", dt$term, ignore.case = TRUE)
  is_inter <- grepl(":", dt$term, fixed = TRUE)
  dt[order(!is_int, is_inter, dt$term)]
}

cat("\n================ Table 2 : 係数パネル ================\n")

coef_all <- rbindlist(lapply(PIPELINES, function(nm) {
  f <- file.path(OUT_DIR, sprintf("bnb_final_%s.rds", nm))
  if (!file.exists(f)) stop("見つかりません: ", f, "（§10-1 を先に実行してください）")
  fin <- readRDS(f)
  tb  <- extract_coef_table(fin, nm)
  tb  <- order_terms(tb)

  mt <- MANUSCRIPT_TERMS[[nm]]
  if (!is.null(mt)) {
    fitted_terms <- setdiff(tb$term, "(Intercept)")
    only_fit <- setdiff(fitted_terms, mt)
    only_ms  <- setdiff(mt, fitted_terms)
    if (length(only_fit) || length(only_ms))
      warning(sprintf("[%s] 原稿と保持項が食い違います。モデルのみ: {%s} / 原稿のみ: {%s}",
                      nm, paste(only_fit, collapse = ", "), paste(only_ms, collapse = ", ")))
  }
  tb
}), use.names = TRUE)

coef_all[, label := relabel_vars(term)]
coef_all[grepl("^\\(?Intercept\\)?$", term, ignore.case = TRUE),
         label := "Intercept (log odds at the reference profile)"]

table2_coef <- coef_all[, .(
  Pipeline      = PIPE_LABELS[pipeline],
  Term          = label,
  `Coefficient (log odds)` = sprintf("%.*f", DIGITS_B, beta),
  `SE`          = ifelse(is.finite(se), sprintf("%.*f", DIGITS_B, se), "—"),
  `Odds ratio (95% CI)` = ifelse(grepl("^Intercept", label), "—",
                                 fmt_ci(OR, OR_lo, OR_hi)),
  `p`           = fmt_p(p)
)]
print(table2_coef, nrows = 200)

fwrite(table2_coef, file.path(OUT_DIR, "table2_model_specification.csv"))

cat("\n================ Table 2 : 仕様パネル ================\n")

centers <- if (file.exists(file.path(OUT_DIR, "bnb_centers.rds")))
  readRDS(file.path(OUT_DIR, "bnb_centers.rds")) else NULL
knots   <- if (file.exists(file.path(OUT_DIR, "bnb_rcs_knots.rds")))
  readRDS(file.path(OUT_DIR, "bnb_rcs_knots.rds")) else NULL
rcsmeta <- if (file.exists(file.path(OUT_DIR, "bnb_rcs_meta.rds")))
  readRDS(file.path(OUT_DIR, "bnb_rcs_meta.rds")) else NULL
findex  <- if (file.exists(file.path(OUT_DIR, "bnb_final_index.rds")))
  as.data.table(readRDS(file.path(OUT_DIR, "bnb_final_index.rds"))) else NULL

const_rows <- rbindlist(lapply(PIPELINES, function(nm) {
  rows <- list()

  cn <- centers[[nm]]
  if (!is.null(cn)) {
    rows[[length(rows) + 1L]] <- data.table(
      pipeline = nm, item = "Median-centring constant (original scale)",
      variable = relabel_vars(names(cn)),
      value    = sprintf("%.4g", as.numeric(cn))
    )
  }

  kn <- knots[[nm]]
  rcs_cont <- rcsmeta[[nm]]$rcs_cont %||% names(kn)
  if (!is.null(kn) && length(rcs_cont)) {
    for (v in intersect(names(kn), rcs_cont)) {
      k_c <- as.numeric(kn[[v]])
      if (!length(k_c) || all(is.na(k_c))) next
      k_o <- if (!is.null(cn) && v %in% names(cn)) k_c + as.numeric(cn[[v]]) else rep(NA_real_, length(k_c))
      rows[[length(rows) + 1L]] <- data.table(
        pipeline = nm, item = "RCS knots, centred scale",
        variable = relabel_vars(v), value = paste(sprintf("%.4g", k_c), collapse = ", "))
      rows[[length(rows) + 1L]] <- data.table(
        pipeline = nm, item = "RCS knots, original scale",
        variable = relabel_vars(v), value = paste(sprintf("%.4g", k_o), collapse = ", "))
    }
    rows[[length(rows) + 1L]] <- data.table(
      pipeline = nm, item = "RCS specification", variable = "—",
      value = sprintf("%d knots at quantiles %s of the observed distribution; basis = Hmisc::rcspline.eval(x, knots, inclx = FALSE)",
                      rcsmeta[[nm]]$n_knots %||% NA_integer_,
                      paste(rcsmeta[[nm]]$knot_probs %||% NA, collapse = "/")))
  }

  rows[[length(rows) + 1L]] <- data.table(
    pipeline = nm, item = "Reference categories", variable = "—",
    value = paste("Disease class = cerebrovascular disease;",
                  "pre-admission long-term care = independent;",
                  "JCS = alert (0-3); sex = male;",
                  "history of cancer = no; history of CNS disease = no"))

  if (!is.null(findex) && "pipeline" %in% names(findex)) {
    fr <- findex[pipeline == nm]
    if (nrow(fr)) {
      est <- fr[[pick_col(fr, c("path", "estimator", "method"))]][1] %||% NA
      mm  <- fr[[pick_col(fr, c("m", "n_imp", "M"))]][1] %||% NA
      nn  <- fr[[pick_col(fr, c("n", "n_used", "N"))]][1] %||% NA
      rows[[length(rows) + 1L]] <- data.table(
        pipeline = nm, item = "Estimator / imputations / n analysed", variable = "—",
        value = sprintf("%s / m = %s / n = %s", est, mm, nn))
    }
  }

  rbindlist(rows, use.names = TRUE, fill = TRUE)
}), use.names = TRUE, fill = TRUE)

const_rows[, Pipeline := PIPE_LABELS[pipeline]]
table2_const <- const_rows[, .(Pipeline, Item = item, Variable = variable, Value = value)]
print(table2_const, nrows = 300)

fwrite(table2_const, file.path(OUT_DIR, "table2_model_constants.csv"))

if (requireNamespace("officer", quietly = TRUE) &&
    requireNamespace("flextable", quietly = TRUE)) {
  ft1 <- flextable::autofit(flextable::flextable(as.data.frame(table2_coef)))
  ft2 <- flextable::autofit(flextable::flextable(as.data.frame(table2_const)))
  doc <- officer::read_docx()
  doc <- officer::body_add_par(doc, "Table 2. Complete specification of the four final models", style = "heading 1")
  doc <- officer::body_add_par(doc, "Panel A. Pooled coefficients of every retained term")
  doc <- flextable::body_add_flextable(doc, ft1)
  doc <- officer::body_add_par(doc, "")
  doc <- officer::body_add_par(doc, "Panel B. Constants needed to reconstruct the linear predictor")
  doc <- flextable::body_add_flextable(doc, ft2)
  print(doc, target = file.path(OUT_DIR, "table2_model_specification.docx"))
  cat("\n貼り付け用 Word 表を書き出しました: data/table2_model_specification.docx\n")
} else {
  cat("\n[注] officer / flextable が無いため .docx は出力していません。CSV を使ってください。\n")
}

cat("\n完了: data/table2_model_specification.csv, data/table2_model_constants.csv\n")
cat("本文 Table 2 の ‡ を付した空欄に、上記 2 ファイルの値をそのまま転記してください。\n")
