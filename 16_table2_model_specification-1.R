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

CI_SOURCE <- "pooled"
if (!CI_SOURCE %in% c("pooled", "wald")) stop("CI_SOURCE は 'pooled' か 'wald'。")

WRITE_CHANGE_LOG <- TRUE

MAIN_PIPELINES <- grep("_A$", PIPELINES, value = TRUE)
SENS_PIPELINES <- grep("_B$", PIPELINES, value = TRUE)

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

ci_kind_of <- function(fin, ci_source) {
  path <- fin$path %||% NA_character_
  clip <- isTRUE(fin$method_detail$clip_ok)
  if (identical(ci_source, "wald"))
    return(list(short = "Wald (normal)",
                long  = "Wald 型（正規近似）。beta ± 1.96 × SE を再計算したもの"))
  if (identical(path, "firth")) {
    if (clip)
      return(list(short = "CLIP profile",
                  long  = paste0("Firth 罰則付き推定に対する combination-of-likelihood-profile ",
                                 "区間（logistf::CLIP.confint）")))
    return(list(short = "Wald (normal, Firth SE)",
                long  = paste0("Firth 罰則付き推定の Rubin プール SE による正規近似区間。",
                               "CLIP.confint が失敗したため profile 区間は得られていない")))
  }
  list(short = "Rubin t",
       long  = paste0("Rubin の規則でプールし、Barnard–Rubin 自由度の t 分布を用いた区間",
                      "（mice::pool）"))
}

p_kind_of <- function(fin) {
  path <- fin$path %||% NA_character_
  clip <- isTRUE(fin$method_detail$clip_ok)
  if (identical(path, "firth")) {
    if (clip) return(list(short = "CLIP profile",
                          long = "CLIP.confint によるプロファイル尤度 p 値"))
    return(list(short = "none", long = "p 値は得られていない（CLIP.confint が失敗）"))
  }
  list(short = "Rubin t (Wald)",
       long = "Rubin プール後の Wald 型 t 検定（mice::pool）")
}

CAND_TERM   <- c("term", "variable", "coef", "name", "rowname")
CAND_EST    <- c("estimate", "beta", "coef_est", "logOR", "est")
CAND_SE     <- c("std.error", "se", "SE", "std_error")
CAND_P      <- c("p.value", "p", "pval", "p_value")
CAND_LO_LOG <- c("conf.low",  "lo", "lower", "ci_lo", "CI_low")
CAND_HI_LOG <- c("conf.high", "hi", "upper", "ci_hi", "CI_high")
CAND_LO_OR  <- c("lcl", "OR_lcl", "or_lcl", "ci_lo_or")
CAND_HI_OR  <- c("ucl", "OR_ucl", "or_ucl", "ci_hi_or")

extract_coef_table <- function(fin, pipeline, ci_source = CI_SOURCE) {
  cand <- c("coef_table", "pooled", "pooled_table", "estimates", "coefs", "summary")
  tab <- NULL
  for (nm in cand) if (!is.null(fin[[nm]])) { tab <- fin[[nm]]; break }
  if (is.null(tab))
    stop(sprintf("[%s] 係数表が見つかりません。bnb_final_%s.rds の構造を確認してください: %s",
                 pipeline, pipeline, paste(names(fin), collapse = ", ")))
  tab <- as.data.table(tab)

  c_term <- pick_col(tab, CAND_TERM)
  c_est  <- pick_col(tab, CAND_EST)
  c_se   <- pick_col(tab, CAND_SE)
  c_p    <- pick_col(tab, CAND_P)
  if (is.null(c_term) || is.null(c_est))
    stop(sprintf("[%s] term / estimate 列を特定できません: %s",
                 pipeline, paste(names(tab), collapse = ", ")))

  out <- data.table(
    pipeline = pipeline,
    term     = as.character(tab[[c_term]]),
    beta     = as.numeric(tab[[c_est]]),
    se       = if (is.null(c_se)) NA_real_ else as.numeric(tab[[c_se]]),
    lo_beta  = NA_real_,
    hi_beta  = NA_real_,
    p        = if (is.null(c_p))  NA_real_ else as.numeric(tab[[c_p]])
  )

  c_lo_log <- pick_col(tab, CAND_LO_LOG); c_hi_log <- pick_col(tab, CAND_HI_LOG)
  c_lo_or  <- pick_col(tab, CAND_LO_OR);  c_hi_or  <- pick_col(tab, CAND_HI_OR)
  ci_col_source <- NA_character_
  heuristic_fired <- FALSE

  if (identical(ci_source, "wald")) {
    ci_col_source <- "recomputed"
  } else if (!is.null(c_lo_log) && !is.null(c_hi_log)) {
    out[, `:=`(lo_beta = as.numeric(tab[[c_lo_log]]),
               hi_beta = as.numeric(tab[[c_hi_log]]))]
    ci_col_source <- sprintf("stored, log-odds scale (%s / %s)", c_lo_log, c_hi_log)
    if (all(is.finite(out$lo_beta)) && all(out$lo_beta > 0) &&
        mean(abs(log(out$lo_beta) - out$beta), na.rm = TRUE) <
        mean(abs(out$lo_beta      - out$beta), na.rm = TRUE)) {
      out[, `:=`(lo_beta = log(lo_beta), hi_beta = log(hi_beta))]
      heuristic_fired <- TRUE
      ci_col_source <- paste0(ci_col_source, " [heuristically re-scaled from OR]")
      warning(sprintf(paste0("[%s] '%s' は log-odds 尺度の名前だが値は OR 尺度と判断して ",
                             "log を取った。列名と保存尺度の対応を確認すること。"),
                      pipeline, c_lo_log))
    }
  } else if (!is.null(c_lo_or) && !is.null(c_hi_or)) {
    lo_or <- as.numeric(tab[[c_lo_or]]); hi_or <- as.numeric(tab[[c_hi_or]])
    if (any(lo_or <= 0, na.rm = TRUE) || any(hi_or <= 0, na.rm = TRUE))
      stop(sprintf("[%s] '%s' / '%s' は OR 尺度のはずだが非正の値を含む。",
                   pipeline, c_lo_or, c_hi_or))
    out[, `:=`(lo_beta = log(lo_or), hi_beta = log(hi_or))]
    ci_col_source <- sprintf("stored, odds-ratio scale (%s / %s), log taken",
                             c_lo_or, c_hi_or)
  } else {
    ci_col_source <- "recomputed"
  }

  z <- stats::qnorm(1 - (1 - CI_LEVEL) / 2)
  n_wald <- sum(is.na(out$lo_beta) & is.finite(out$se))
  out[is.na(lo_beta) & is.finite(se), `:=`(lo_beta = beta - z * se,
                                           hi_beta = beta + z * se)]
  if (identical(ci_col_source, "recomputed"))
    ci_col_source <- sprintf("recomputed as beta +/- %.4f * SE (normal Wald)", z)
  else if (n_wald > 0)
    ci_col_source <- sprintf("%s; %d term(s) filled by normal Wald", ci_col_source, n_wald)

  out[, `:=`(OR = exp(beta), OR_lo = exp(lo_beta), OR_hi = exp(hi_beta))]

  data.table::setattr(out, "ci_col_source",   ci_col_source)
  data.table::setattr(out, "n_wald_filled",   n_wald)
  data.table::setattr(out, "heuristic_fired", heuristic_fired)
  out[]
}

PROV_ATTRS <- c("ci_col_source", "n_wald_filled", "heuristic_fired")
order_terms <- function(dt) {
  is_int   <- grepl("^\\(?Intercept\\)?$", dt$term, ignore.case = TRUE)
  is_inter <- grepl(":", dt$term, fixed = TRUE)
  out <- dt[order(!is_int, is_inter, dt$term)]
  for (a in PROV_ATTRS) data.table::setattr(out, a, attr(dt, a, exact = TRUE))
  out
}

cat("\n================ Table 2 : 係数パネル ================\n")
cat(sprintf("CI_SOURCE = '%s'\n", CI_SOURCE))

fins <- lapply(PIPELINES, function(nm) {
  f <- file.path(OUT_DIR, sprintf("bnb_final_%s.rds", nm))
  if (!file.exists(f)) stop("見つかりません: ", f, "（§10-1 を先に実行してください）")
  readRDS(f)
})
names(fins) <- PIPELINES

prov_rows <- list()

coef_all <- rbindlist(lapply(PIPELINES, function(nm) {
  fin <- fins[[nm]]
  tb  <- extract_coef_table(fin, nm)
  tb  <- order_terms(tb)

  ck <- ci_kind_of(fin, CI_SOURCE)
  pk <- p_kind_of(fin)
  prov_rows[[nm]] <<- data.table(
    pipeline          = nm,
    label             = unname(PIPE_LABELS[[nm]]),
    table             = if (nm %in% MAIN_PIPELINES) "Table 2" else "Supplementary Table S8",
    estimator         = fin$path %||% NA_character_,
    clip_ok           = isTRUE(fin$method_detail$clip_ok),
    n_terms           = nrow(tb),
    ci_source_setting = CI_SOURCE,
    ci_columns_used   = attr(tb, "ci_col_source"),
    ci_kind_short     = ck$short,
    ci_kind_long      = ck$long,
    p_kind_short      = pk$short,
    p_kind_long       = pk$long,
    n_wald_filled     = attr(tb, "n_wald_filled"),
    heuristic_fired   = attr(tb, "heuristic_fired"))

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

make_panel <- function(dt) dt[, .(
  Pipeline      = PIPE_LABELS[pipeline],
  Term          = label,
  `Coefficient (log odds)` = sprintf("%.*f", DIGITS_B, beta),
  `SE`          = ifelse(is.finite(se), sprintf("%.*f", DIGITS_B, se), "—"),
  `Odds ratio (95% CI)` = ifelse(grepl("^Intercept", label), "—",
                                 fmt_ci(OR, OR_lo, OR_hi)),
  `p`           = fmt_p(p)
)]

table2_coef <- make_panel(coef_all)
print(table2_coef, nrows = 200)
fwrite(table2_coef, file.path(OUT_DIR, "table2_model_specification.csv"))

tab_main <- make_panel(coef_all[pipeline %in% MAIN_PIPELINES])
tab_sens <- make_panel(coef_all[pipeline %in% SENS_PIPELINES])
fwrite(tab_main, file.path(OUT_DIR, "table2_main_analysis.csv"))
fwrite(tab_sens, file.path(OUT_DIR, "table_s8_sensitivity.csv"))
cat(sprintf("\n本文 Table 2 = %s / Supplementary Table S8 = %s\n",
            paste(MAIN_PIPELINES, collapse = " + "),
            paste(SENS_PIPELINES, collapse = " + ")))

prov <- rbindlist(prov_rows, use.names = TRUE, fill = TRUE)
fwrite(prov, file.path(OUT_DIR, "table2_ci_provenance.csv"))

cat("\n================ 区間と p 値の出所 ================\n")
print(prov[, .(pipeline, table, estimator, clip_ok, ci_kind_short, p_kind_short,
               n_wald_filled)])
if (any(prov$heuristic_fired))
  cat("\n[警告] 尺度判定のヒューリスティックが発動したパイプラインがあります。",
      "列名と保存尺度の対応を確認してください。\n")

if (WRITE_CHANGE_LOG) {
  other <- if (identical(CI_SOURCE, "pooled")) "wald" else "pooled"
  alt <- rbindlist(lapply(PIPELINES, function(nm)
    order_terms(extract_coef_table(fins[[nm]], nm, ci_source = other))),
    use.names = TRUE)

  cmp <- merge(coef_all[, .(pipeline, term, OR, OR_lo, OR_hi)],
               alt[, .(pipeline, term, OR_lo_alt = OR_lo, OR_hi_alt = OR_hi)],
               by = c("pipeline", "term"), all.x = TRUE, sort = FALSE)
  cmp[, `:=`(
    label      = relabel_vars(term),
    used       = fmt_ci(OR, OR_lo,     OR_hi),
    alternative= fmt_ci(OR, OR_lo_alt, OR_hi_alt),
    rel_change_lcl = 100 * (OR_lo / OR_lo_alt - 1),
    rel_change_ucl = 100 * (OR_hi / OR_hi_alt - 1))]
  cmp[, moves_at_printed_precision := used != alternative]
  cmp[, `:=`(ci_source_used = CI_SOURCE, ci_source_alternative = other)]
  setcolorder(cmp, c("pipeline", "term", "label", "ci_source_used",
                     "ci_source_alternative", "used", "alternative",
                     "moves_at_printed_precision"))
  fwrite(cmp, file.path(OUT_DIR, "table2_ci_change_log.csv"))

  nmv <- sum(cmp$moves_at_printed_precision, na.rm = TRUE)
  cat(sprintf("\n================ '%s' と '%s' の差分 ================\n", CI_SOURCE, other))
  cat(sprintf("小数第%d位で表示が変わる項: %d / %d\n", DIGITS_OR, nmv, nrow(cmp)))
  if (nmv > 0) {
    print(cmp[moves_at_printed_precision == TRUE,
              .(pipeline, term = label, used, alternative,
                lcl_pct = round(rel_change_lcl, 3),
                ucl_pct = round(rel_change_ucl, 3))])
    cat("\n上の行が、原稿本文・Table 2・S8・S3 で書き換えが必要になるセルです。\n")
    cat("Results と抄録が引用している区間も同じ値なので、あわせて点検してください。\n")
  } else {
    cat("表示上の変化はありません。丸めの境界にある項がなかったということです。\n")
  }
}

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

  pr <- prov[pipeline == nm]
  if (nrow(pr)) {
    rows[[length(rows) + 1L]] <- data.table(
      pipeline = nm, item = "Confidence interval", variable = "—", value = pr$ci_kind_long[1])
    rows[[length(rows) + 1L]] <- data.table(
      pipeline = nm, item = "p value", variable = "—", value = pr$p_kind_long[1])
  }

  rbindlist(rows, use.names = TRUE, fill = TRUE)
}), use.names = TRUE, fill = TRUE)

const_rows[, Pipeline := PIPE_LABELS[pipeline]]
table2_const <- const_rows[, .(Pipeline, Item = item, Variable = variable, Value = value)]
print(table2_const, nrows = 300)

fwrite(table2_const, file.path(OUT_DIR, "table2_model_constants.csv"))

footnote_for <- function(pipes, table_name) {
  pr <- prov[pipeline %in% pipes]
  if (!nrow(pr)) return(NULL)
  kinds <- unique(pr$ci_kind_long)
  est   <- unique(pr$estimator)
  s <- sprintf("%s. ", table_name)
  s <- paste0(s, "すべての連続項は中央値でセンタリングしてあるため、各係数は",
              "センタリング後の変数1単位あたりの対数オッズの変化である。")
  if (length(est) == 1L && identical(est, "glm")) {
    s <- paste0(s, "本表のモデルはいずれも通常のロジスティック回帰（最尤法）で推定した。")
  } else {
    s <- paste0(s, "推定器はパイプラインにより異なる：",
                paste(sprintf("%s は %s", pr$label, pr$estimator), collapse = "、"), "。")
  }
  if (length(kinds) == 1L) {
    s <- paste0(s, "信頼区間は", kinds, "である。")
  } else {
    s <- paste0(s, "信頼区間の種類はモデルにより異なる：",
                paste(sprintf("%s は%s", pr$label, pr$ci_kind_long), collapse = "、"), "。")
  }
  pk <- unique(pr$p_kind_long)
  s <- paste0(s, if (length(pk) == 1L) paste0("p 値は", pk, "である。")
              else paste0("p 値：", paste(sprintf("%s は%s", pr$label, pr$p_kind_long),
                                          collapse = "、"), "。"))
  if (any(pr$p_kind_short != "none"))
    s <- paste0(s, "本表の p 値は各係数に対する検定であり、",
                "交互作用項の採用を判定した §9 のプール尤度比検定の p 値とは異なる。")
  s <- paste0(s, "本表は Supplementary Material のスクリプト16が当てはめ済みモデル",
              "オブジェクトから直接書き出したものである。当該スクリプトは Figure 4 および ",
              "Supplementary Table S3 が読むのと同一の区間列（lcl / ucl）を用いる。")
  s
}

fn_main <- footnote_for(MAIN_PIPELINES, "Table 2")
fn_sens <- footnote_for(SENS_PIPELINES, "Supplementary Table S8")
fn_txt  <- paste(c(fn_main, "", fn_sens), collapse = "\n")
fn_con <- file(file.path(OUT_DIR, "table2_footnotes.txt"), open = "wt", encoding = "UTF-8")
writeLines(enc2utf8(fn_txt), fn_con, useBytes = TRUE)
close(fn_con)

cat("\n================ 生成した脚注文 ================\n")
cat(fn_txt, "\n")

if (requireNamespace("officer", quietly = TRUE) &&
    requireNamespace("flextable", quietly = TRUE)) {
  ftt <- function(x) flextable::autofit(flextable::flextable(as.data.frame(x)))
  doc <- officer::read_docx()
  doc <- officer::body_add_par(doc, "Table 2. Complete specification of the final models",
                               style = "heading 1")
  if (nrow(tab_main)) {
    doc <- officer::body_add_par(doc, "Main analysis (worst-value system)")
    doc <- flextable::body_add_flextable(doc, ftt(tab_main))
    doc <- officer::body_add_par(doc, "")
    if (!is.null(fn_main)) doc <- officer::body_add_par(doc, fn_main)
    doc <- officer::body_add_par(doc, "")
  }
  if (nrow(tab_sens)) {
    doc <- officer::body_add_par(doc, "Supplementary Table S8. Sensitivity analysis (missing-category system)",
                                 style = "heading 1")
    doc <- flextable::body_add_flextable(doc, ftt(tab_sens))
    doc <- officer::body_add_par(doc, "")
    if (!is.null(fn_sens)) doc <- officer::body_add_par(doc, fn_sens)
    doc <- officer::body_add_par(doc, "")
  }
  doc <- officer::body_add_par(doc, "Constants needed to reconstruct the linear predictor (Supplementary Table S9)",
                               style = "heading 1")
  doc <- flextable::body_add_flextable(doc, ftt(table2_const))
  print(doc, target = file.path(OUT_DIR, "table2_model_specification.docx"))
  cat("\n貼り付け用 Word 表を書き出しました: data/table2_model_specification.docx\n")
} else {
  cat("\n[注] officer / flextable が無いため .docx は出力していません。CSV を使ってください。\n")
}

cat("\n完了（改訂版 2026-08-16）:\n")
cat("  data/table2_model_specification.csv   … 4本まとめた係数パネル（従来互換）\n")
cat("  data/table2_main_analysis.csv         … 本文 Table 2\n")
cat("  data/table_s8_sensitivity.csv         … Supplementary Table S8\n")
cat("  data/table2_model_constants.csv       … 仕様パネル（S9）\n")
cat("  data/table2_ci_provenance.csv         … 区間と p 値の出所\n")
if (WRITE_CHANGE_LOG)
  cat("  data/table2_ci_change_log.csv         … 書き換えが必要なセルの一覧\n")
cat("  data/table2_footnotes.txt             … 生成した脚注文\n")
cat("\n次の作業:\n")
cat("  1. table2_ci_change_log.csv の moves_at_printed_precision = TRUE の行を、\n")
cat("     Table 2・S8・S3・Results 本文・抄録で書き換える。\n")
cat("  2. S8 の脚注を table2_footnotes.txt の文に差し替える。\n")
cat("     CLIP が失敗していた場合、従来の「combination-of-likelihood-profile 区間」は誤りである。\n")
cat("  3. 15_reporting_supplements.2026.08.16.R の ADJ_CI_SOURCE を 'stored' に変える。\n")
cat("     これで S3・Table 2・S8・Figure 4 が同一の区間を指すようになる。\n")
