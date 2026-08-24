library(data.table)
library(here)
library(ggplot2)

need_pkg <- function(p) if (!requireNamespace(p, quietly = TRUE))
  stop(sprintf("パッケージ '%s' が必要です。install.packages('%s')", p, p))
need_pkg("logistf")
HAS_PROC <- requireNamespace("pROC", quietly = TRUE)
have_tic <- requireNamespace("tictoc", quietly = TRUE)

LABELS_PATH <- here::here("00_labels.R")
if (!exists("relabel_vars")) {
  if (!file.exists(LABELS_PATH)) stop("00_labels.R が見つかりません: ", LABELS_PATH)
  source(LABELS_PATH)
}

OUT_DIR <- here::here("data")
FIG_DIR <- here::here("figures")
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

PIPELINES <- names(PIPE_LABELS)
CI_LEVEL  <- 0.95
PROB_EPS  <- 1e-8
DIGITS    <- 3L

DO_LOVO  <- TRUE
N_BOOT   <- 200L
M_BOOT   <- 5L
BOOT_IMP_IDS <- NULL

BETA_MAX     <- 15
LOGOR_CIW_MAX<- 8
FMI_HIGH     <- 0.5

DAUC_ROBUST  <- 0.005
DAUC_NULL    <- 0.000
OPT_MARGIN   <- 0.010

set.seed(20240601)

read_safe <- function(f) if (file.exists(f)) readRDS(f) else NULL
pipe_disp <- function(x) vapply(as.character(x), relabel_pipeline, character(1))

reduce_imps <- function(obj, n0, outlier_ids) {
  imps_full <- obj$smcfcs$impDatasets
  bad <- sort(unique(outlier_ids)); bad <- bad[bad %in% seq_len(n0)]
  if (length(bad) > 0) lapply(imps_full, function(d) d[-bad, , drop = FALSE]) else imps_full
}
pick_imps <- function(m, k, ids = NULL) {
  if (!is.null(ids)) return(sort(unique(ids[ids >= 1 & ids <= m])))
  if (is.null(k) || k >= m) return(seq_len(m))
  unique(round(seq(1, m, length.out = k)))
}

col_to_varkey <- function(s) {
  s <- as.character(s)
  if (s == "(Intercept)") return("(Intercept)")
  if (grepl(":", s, fixed = TRUE)) {
    parts <- strsplit(s, ":", fixed = TRUE)[[1]]
    ks <- vapply(parts, col_to_varkey, character(1))
    return(paste(sort(ks), collapse = ":"))
  }
  if (grepl("_measurable", s))  return(s)
  if (grepl("_c[0-9]+$", s))    return(sub("_c[0-9]+$", "", s))
  if (grepl("_c$", s))          return(sub("_c$", "", s))
  hit <- names(VAR_LABELS)[vapply(names(VAR_LABELS),
                                  function(v) startsWith(s, v), logical(1))]
  if (length(hit) > 0) return(hit[which.max(nchar(hit))])
  s
}
label_var_keys <- function(label) {
  parts <- strsplit(as.character(label), ":", fixed = TRUE)[[1]]
  unique(vapply(parts, function(p) {
    if (grepl("_measurable", p)) return(p)
    if (grepl("_c[0-9]+$", p))   return(sub("_c[0-9]+$", "", p))
    if (grepl("_c$", p))         return(sub("_c$", "", p))
    hit <- names(VAR_LABELS)[vapply(names(VAR_LABELS),
                                    function(v) startsWith(p, v), logical(1))]
    if (length(hit) > 0) return(hit[which.max(nchar(hit))])
    p
  }, character(1)))
}
varkey_disp <- function(k) vapply(k, function(s) {
  if (grepl(":", s, fixed = TRUE)) {
    parts <- strsplit(s, ":", fixed = TRUE)[[1]]
    return(paste(vapply(parts, relabel_vars, character(1)), collapse = " x "))
  }
  relabel_vars(s)
}, character(1), USE.NAMES = FALSE)

fit_firth_coef <- function(smf, d) {
  fit <- tryCatch(logistf::logistf(smf, data = d, pl = FALSE), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  b <- coef(fit); if (is.null(names(b))) names(b) <- rownames(vcov(fit)); b
}
predict_lp <- function(coef_named, smf, newdata) {
  mm <- tryCatch(stats::model.matrix(stats::delete.response(stats::terms(smf)), data = newdata),
                 error = function(e) NULL)
  if (is.null(mm)) return(rep(NA_real_, nrow(newdata)))
  common <- intersect(colnames(mm), names(coef_named))
  if (length(common) == 0) return(rep(NA_real_, nrow(newdata)))
  as.numeric(mm[, common, drop = FALSE] %*% coef_named[common])
}
prob_from_lp <- function(lp) { p <- stats::plogis(lp); pmin(pmax(p, PROB_EPS), 1 - PROB_EPS) }

auc_one <- function(y, p) {
  pos <- p[y == 1]; neg <- p[y == 0]
  if (length(pos) == 0 || length(neg) == 0) return(NA_real_)
  if (HAS_PROC) {
    r <- tryCatch(pROC::roc(y, p, quiet = TRUE, levels = c(0, 1), direction = "<"),
                  error = function(e) NULL)
    if (!is.null(r)) return(as.numeric(pROC::auc(r)))
  }
  rr <- rank(c(pos, neg))
  (sum(rr[seq_along(pos)]) - length(pos) * (length(pos) + 1) / 2) / (length(pos) * length(neg))
}
cal_slope <- function(y, lp) {
  fs <- tryCatch(suppressWarnings(glm(y ~ lp, family = binomial())), error = function(e) NULL)
  if (is.null(fs) || !("lp" %in% names(coef(fs)))) return(NA_real_)
  unname(coef(fs)["lp"])
}

bootstrap_optimism_firth <- function(imps, smf, n_boot = N_BOOT, m_use = M_BOOT, ids = BOOT_IMP_IDS) {
  m <- length(imps); use_idx <- pick_imps(m, m_use, ids)
  tl <- attr(stats::terms(smf), "term.labels")
  app_a <- app_s <- opt_a <- opt_s <- numeric(0)
  for (i in use_idx) {
    d <- imps[[i]]; y <- as.integer(d$good); n <- nrow(d)
    b0 <- fit_firth_coef(smf, d); if (is.null(b0)) next
    p0 <- prob_from_lp(predict_lp(b0, smf, d))
    a_app <- auc_one(y, p0); s_app <- cal_slope(y, stats::qlogis(p0))
    if (!is.finite(a_app)) next
    oa <- os <- numeric(0)
    for (bnum in seq_len(n_boot)) {
      idx <- sample.int(n, n, replace = TRUE); db <- d[idx, , drop = FALSE]
      bb <- fit_firth_coef(smf, db); if (is.null(bb)) next
      pbb <- prob_from_lp(predict_lp(bb, smf, db))
      pbo <- prob_from_lp(predict_lp(bb, smf, d))
      a_b <- auc_one(as.integer(db$good), pbb); a_t <- auc_one(y, pbo)
      s_b <- cal_slope(as.integer(db$good), stats::qlogis(pbb)); s_t <- cal_slope(y, stats::qlogis(pbo))
      if (is.finite(a_b) && is.finite(a_t)) oa <- c(oa, a_b - a_t)
      if (is.finite(s_b) && is.finite(s_t)) os <- c(os, s_b - s_t)
    }
    app_a <- c(app_a, a_app); app_s <- c(app_s, s_app)
    opt_a <- c(opt_a, mean(oa, na.rm = TRUE)); opt_s <- c(opt_s, mean(os, na.rm = TRUE))
  }
  if (length(app_a) == 0)
    return(list(AUC_apparent = NA_real_, AUC_optimism = NA_real_, AUC_corrected = NA_real_,
                slope_apparent = NA_real_, slope_optimism = NA_real_, slope_corrected = NA_real_))
  list(AUC_apparent = mean(app_a), AUC_optimism = mean(opt_a, na.rm = TRUE),
       AUC_corrected = mean(app_a) - mean(opt_a, na.rm = TRUE),
       slope_apparent = mean(app_s, na.rm = TRUE), slope_optimism = mean(opt_s, na.rm = TRUE),
       slope_corrected = mean(app_s, na.rm = TRUE) - mean(opt_s, na.rm = TRUE))
}

estimate_stability <- function(pooled, unit_key, unit_type) {
  if (is.null(pooled) || !nrow(pooled)) return(NULL)
  pl <- as.data.table(copy(pooled))
  pl <- pl[term != "(Intercept)"]
  pl[, vk := vapply(term, col_to_varkey, character(1))]
  rows <- if (unit_type == "interaction") pl[vk == unit_key] else pl[vk == unit_key]
  if (!nrow(rows)) return(NULL)
  has_fmi <- "fmi" %in% names(rows)
  logOR <- rows$estimate
  lcl <- rows$lcl; ucl <- rows$ucl
  ciw  <- log(pmax(ucl, .Machine$double.eps)) - log(pmax(lcl, .Machine$double.eps))
  rep_i <- which.max(abs(logOR))
  list(
    OR = rows$OR[rep_i], lcl = lcl[rep_i], ucl = ucl[rep_i],
    logOR = logOR[rep_i], ci_logwidth = max(ciw, na.rm = TRUE),
    ci_finite = all(is.finite(lcl) & is.finite(ucl) & lcl > 0),
    excludes_1 = any((lcl > 1 & ucl > 1) | (lcl < 1 & ucl < 1)),
    fmi = if (has_fmi) max(rows$fmi, na.rm = TRUE) else NA_real_,
    beta_extreme = any(abs(logOR) > BETA_MAX),
    n_terms = nrow(rows))
}

model_units <- function(md, obj) {
  if (md == "Provisional") {
    keys <- obj$prov_vars; groups <- obj$groups
    keep <- if (!is.null(obj$keep_vars)) obj$keep_vars else character(0)
    lapply(keys, function(k) list(type = "main", key = k,
      drop_labels = groups[[k]], forced = k %in% keep))
  } else {
    groups <- obj$groups; main_keys <- obj$main_vars
    int_groups <- obj$selected_int_groups; int_keys <- obj$selected_int_keys
    if (is.null(int_keys)) int_keys <- character(0)
    int_comp <- lapply(int_keys, function(ik)
      unique(unlist(lapply(int_groups[[ik]], label_var_keys))))
    names(int_comp) <- int_keys
    units <- lapply(main_keys, function(k) {
      involving <- int_keys[vapply(int_keys, function(ik) k %in% int_comp[[ik]], logical(1))]
      drop <- c(groups[[k]], unlist(int_groups[involving], use.names = FALSE))
      list(type = "main", key = k, drop_labels = unique(drop), forced = TRUE)
    })
    units2 <- lapply(int_keys, function(ik)
      list(type = "interaction", key = paste(sort(int_comp[[ik]]), collapse = ":"),
           int_key = ik, drop_labels = int_groups[[ik]], forced = FALSE))
    c(units, units2)
  }
}

analyze_firth_model <- function(nm, md, obj, imps) {
  smf_full <- stats::as.formula(if (md == "Provisional") obj$smformula_prov else obj$smformula_final)
  all_labels <- attr(stats::terms(smf_full), "term.labels")
  if (length(all_labels) == 0) { message(sprintf("  [%s] 切片のみ → スキップ。", md)); return(NULL) }

  if (have_tic) tictoc::tic(sprintf("  [%s] full-model bootstrap (B=%d)", md, N_BOOT))
  full <- bootstrap_optimism_firth(imps, smf_full)
  if (have_tic) tictoc::toc()

  units <- model_units(md, obj)
  rows <- list()
  for (u in units) {
    red_labels <- setdiff(all_labels, u$drop_labels)
    smf_red <- if (length(red_labels) == 0) stats::as.formula("good ~ 1")
               else stats::reformulate(red_labels, response = "good")
    red <- bootstrap_optimism_firth(imps, smf_red)

    es <- estimate_stability(obj$pooled, u$key, u$type)
    dauc <- full$AUC_corrected - red$AUC_corrected
    dopt <- full$AUC_optimism  - red$AUC_optimism

    rows[[length(rows) + 1L]] <- data.table(
      pipeline = nm, model = md, type = u$type,
      var = u$key, variable = varkey_disp(u$key), forced_keep = isTRUE(u$forced),
      AUC_corr_full = full$AUC_corrected, AUC_corr_reduced = red$AUC_corrected,
      dAUC_corrected = dauc,
      optimism_full = full$AUC_optimism, optimism_reduced = red$AUC_optimism, d_optimism = dopt,
      slope_corr_full = full$slope_corrected, slope_corr_reduced = red$slope_corrected,
      OR = if (!is.null(es)) es$OR else NA_real_,
      lcl = if (!is.null(es)) es$lcl else NA_real_,
      ucl = if (!is.null(es)) es$ucl else NA_real_,
      ci_logwidth = if (!is.null(es)) es$ci_logwidth else NA_real_,
      ci_finite = if (!is.null(es)) es$ci_finite else NA,
      excludes_1 = if (!is.null(es)) es$excludes_1 else NA,
      fmi = if (!is.null(es)) es$fmi else NA_real_,
      beta_extreme = if (!is.null(es)) es$beta_extreme else NA)
  }
  dt <- rbindlist(rows, fill = TRUE)

  dt[, identifiable := (is.na(ci_finite) | ci_finite) &
                       (is.na(beta_extreme) | !beta_extreme) &
                       (is.na(ci_logwidth) | ci_logwidth <= LOGOR_CIW_MAX) &
                       (is.na(fmi) | fmi <= FMI_HIGH)]
  dt[, adds_optimism := is.finite(d_optimism) & d_optimism > OPT_MARGIN]
  dt[, contributes  := is.finite(dAUC_corrected) & dAUC_corrected >= DAUC_ROBUST]
  dt[, no_contrib   := is.finite(dAUC_corrected) & dAUC_corrected <= DAUC_NULL]

  dt[, verdict := fifelse(
        contributes & identifiable & !adds_optimism, "robust",
      fifelse(
        (no_contrib) | adds_optimism | (identifiable == FALSE), "fragile",
        "uncertain"))]
  setorder(dt, -dAUC_corrected)
  dt
}

path_rows <- list(); stab_all <- list()

for (nm in PIPELINES) {
  prov <- read_safe(file.path(OUT_DIR, sprintf("bnb_provisional_%s.rds", nm)))
  fin  <- read_safe(file.path(OUT_DIR, sprintf("bnb_final_%s.rds",       nm)))
  obj  <- read_safe(file.path(OUT_DIR, sprintf("bnb_imp_%s.rds",         nm)))

  models <- list()
  if (!is.null(prov)) models$Provisional <- prov
  if (!is.null(fin))  models$Final       <- fin
  if (length(models) == 0) { message(sprintf("[%s] 暫定/最終モデルなし。スキップ。", nm)); next }

  message(sprintf("\n================ [%s] §12-2 経路判定・変数安定性 ================", nm))

  for (md in names(models)) {
    mo <- models[[md]]
    path <- mo$path
    path_rows[[paste0(nm, "_", md)]] <- data.table(
      pipeline = nm, model = md, path = path,
      n0 = mo$n0, n = mo$n,
      estimator = if (path == "firth") "logistf (Firth)" else "glm",
      ci_method = if (path == "firth") "CLIP profile (MI-pooled)" else "Wald/Rubin")
    message(sprintf("  [%s] 採用経路 = %s", md, toupper(path)))

    if (path != "firth") {
      message("    → glm 経路。LOVO(Firth) はスキップ（glm では §12 の選択頻度・optimism が有効）。")
      next
    }
    if (!DO_LOVO) { message("    DO_LOVO=FALSE のため LOVO をスキップ。"); next }
    if (is.null(obj)) { message("    代入データなし → LOVO 不可。"); next }

    imps <- reduce_imps(obj, mo$n0, mo$outlier_ids)
    res  <- analyze_firth_model(nm, md, mo, imps)
    if (is.null(res)) next
    stab_all[[paste0(nm, "_", md)]] <- res
    saveRDS(res, file.path(OUT_DIR, sprintf("bnb_var_stability_%s_%s.rds", nm, tolower(md))))

    message(sprintf("    変数安定性（ΔAUC_corrected 降順, B=%d, 代入 %s 本）:",
                    N_BOOT, ifelse(is.null(M_BOOT), "all", as.character(M_BOOT))))
    print(res[, .(variable, type,
                  dAUC_corr = round(dAUC_corrected, DIGITS),
                  d_optim = round(d_optimism, DIGITS),
                  OR = round(OR, 2),
                  CI = ifelse(is.finite(lcl) & is.finite(ucl),
                              sprintf("[%.2f, %.2f]", lcl, ucl), "—"),
                  fmi = round(fmi, 2),
                  ident = identifiable, verdict)], row.names = FALSE)

    d <- copy(res)[is.finite(dAUC_corrected)]
    if (nrow(d)) {
      d[, ylab := factor(variable, levels = d[order(dAUC_corrected), variable])]
      g <- ggplot(d, aes(dAUC_corrected, ylab, colour = verdict)) +
        geom_vline(xintercept = c(DAUC_NULL, DAUC_ROBUST),
                   linetype = c("solid", "dashed"), colour = c("grey60", "#2166AC")) +
        geom_segment(aes(x = 0, xend = dAUC_corrected, y = ylab, yend = ylab),
                     colour = "grey75", linewidth = 0.4) +
        geom_point(size = 2.6) +
        scale_colour_manual(values = c(robust = "#2166AC", uncertain = "#F4A582",
                                       fragile = "#B2182B"), name = NULL, drop = FALSE) +
        labs(title = sprintf("%s - %s model: variable stability (Firth)",
                             relabel_pipeline(nm), md),
             subtitle = "LOVO optimism-corrected AUC drop. Right = removing the variable hurts generalization (robust contributor).",
             x = expression(Delta * "AUC"[corrected] * "  (full - without variable)"), y = NULL) +
        theme_bw(base_size = 9) +
        theme(legend.position = "top",
              plot.title = element_text(size = 11, face = "bold"),
              plot.subtitle = element_text(size = 7.5, colour = "grey30"),
              panel.grid.minor = element_blank())
      ggsave(file.path(FIG_DIR, sprintf("fig_var_stability_%s_%s.png", nm, tolower(md))),
             g, width = 7.2, height = max(2.6, 0.4 * nrow(d) + 1.4), dpi = 150)
    }
  }
}

if (length(path_rows)) {
  path_dt <- rbindlist(path_rows, fill = TRUE)
  path_dt[, pipeline_disp := pipe_disp(pipeline)]
  setcolorder(path_dt, c("pipeline", "pipeline_disp", "model", "path", "estimator", "ci_method"))
  saveRDS(path_dt, file.path(OUT_DIR, "bnb_path_summary.rds"))
  fwrite(path_dt, file.path(OUT_DIR, "bnb_path_summary.csv"))
} else path_dt <- data.table()

if (length(stab_all)) {
  stab_dt <- rbindlist(stab_all, fill = TRUE)
  stab_dt[, pipeline_disp := pipe_disp(pipeline)]
  saveRDS(stab_dt, file.path(OUT_DIR, "bnb_var_stability.rds"))
  fwrite(stab_dt, file.path(OUT_DIR, "bnb_var_stability.csv"))
} else stab_dt <- data.table()

message("\n==== §12-2 採用経路の一覧 ====")
if (nrow(path_dt))
  print(path_dt[, .(Pipeline = pipeline_disp, Model = model, Path = toupper(path),
                    Estimator = estimator)], row.names = FALSE)

message("\n==== §12-2 Firth モデルの変数安定性 verdict 集計 ====")
if (nrow(stab_dt)) {
  print(stab_dt[, .N, by = .(Pipeline = pipe_disp(pipeline), Model = model, verdict)][
        order(Pipeline, Model, verdict)], row.names = FALSE)
  message("\n  ■ fragile（不安定／過適合誘発）と判定された変数:")
  fr <- stab_dt[verdict == "fragile"]
  if (nrow(fr))
    print(fr[, .(Pipeline = pipe_disp(pipeline), Model = model, Variable = variable,
                 dAUC_corr = round(dAUC_corrected, DIGITS),
                 d_optim = round(d_optimism, DIGITS),
                 fmi = round(fmi, 2), identifiable)], row.names = FALSE)
  else message("    （なし）")
} else message("  Firth 経路のモデルがない、または LOVO 未実行のため安定性結果はありません。")

message("\n§12-2 完了。")
message("  ・採用経路一覧               : data/bnb_path_summary.csv")
message("  ・変数安定性（Firth, 統合）   : data/bnb_var_stability.csv（モデル別 .rds も保存）")
message("  ・図                         : figures/fig_var_stability_<pipeline>_<model>.png")
message("判定の考え方:")
message("  予測側(汎化): LOVO で当該変数(＋それを含む交互作用=周辺性)を除き、Harrell 流 optimism 補正の")
message("              ΔAUC_corrected と Δoptimism を算出。除くと補正AUCが下がる(寄与)・入れると optimism")
message("              が増える(過適合誘発)で評価。selection 頻度の飽和に依らない予測ベースの安定性。")
message("  推定側(同定): §8/§10 保存の pooled から CLIP CI(MIプール・プロファイル)・FMI を読み、有限で")
message("              常識的幅か / OR CI が 1 を跨がないか / FMI が高くないかで同定性を評価。")
message(sprintf("  総合 verdict: robust=寄与あり&同定良&optimism増やさない / fragile=寄与なし or 過適合誘発 or 同定不良 / それ以外 uncertain。"))
message("  ※閾値(DAUC_ROBUST/OPT_MARGIN/FMI_HIGH/LOGOR_CIW_MAX/BETA_MAX)はヒューリスティック。EPV 不足下では")
message("    『方向は信頼できるが大きさは不確実』に留まる変数が出るのが正常で、verdict は補助的な目安。")
