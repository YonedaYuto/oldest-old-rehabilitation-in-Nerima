library(data.table)
library(here)

need_pkg <- function(p) if (!requireNamespace(p, quietly = TRUE))
  stop(sprintf("パッケージ '%s' が必要です。install.packages('%s')", p, p))
need_pkg("mice")
need_pkg("logistf")

HAS_DETECTSEP <- requireNamespace("detectseparation", quietly = TRUE)
HAS_BRGLM2    <- requireNamespace("brglm2",          quietly = TRUE)
if (!HAS_DETECTSEP && !HAS_BRGLM2)
  stop("分離検出に detectseparation か brglm2 のいずれかが必要です。",
       "install.packages('detectseparation')")

LABELS_PATH <- here::here("00_labels.R")
if (!exists("relabel_vars")) {
  if (!file.exists(LABELS_PATH)) stop("00_labels.R が見つかりません: ", LABELS_PATH)
  source(LABELS_PATH)
}

OUT_DIR <- here::here("data")
have_tic <- requireNamespace("tictoc", quietly = TRUE)

CI_LEVEL <- 0.95

GLM_LRT_METHOD <- "D3"

ALPHA <- 0.05
BONFERRONI_M <- "candidates"

BETA_MAX <- 15
SE_MAX   <- 30
PROB_EPS <- 1e-6
FIRTH_LOCK_AFTER_SWITCH <- FALSE
AUTO_FALLBACK_FIRTH <- TRUE

EPV_PER_VAR <- 10
EPV_EVENTS  <- "min"
ENFORCE_EPV <- TRUE

FORCE_PATH <- list(
  severe_A  = NULL, severe_B  = NULL,
  elderly_A = NULL, elderly_B = NULL
)

KEEP_VARS <- list(
  severe_A  = character(0), severe_B  = character(0),
  elderly_A = character(0), elderly_B = character(0)
)
FORCE_KEEP_SUBSET_VAR <- FALSE

M_SELECT_MAX <- NULL

CLIP_VARS <- list(
  severe_A  = NULL, severe_B  = NULL,
  elderly_A = NULL, elderly_B = NULL
)
LOGISTF_PL <- FALSE

if (!exists("est_index"))
  est_index <- readRDS(file.path(OUT_DIR, "bnb_estimation_index.rds"))

SUBSET_VAR <- list(severe_A = "mFIM_in", severe_B = "mFIM_in",
                   elderly_A = "age",    elderly_B = "age")

detect_one <- function(smf, dat) {
  detector <- if (HAS_DETECTSEP) detectseparation::detect_separation
              else                brglm2::detect_separation
  ds <- tryCatch(
    stats::glm(smf, data = dat, family = binomial("logit"), method = detector),
    error = function(e) e)
  if (inherits(ds, "error"))
    return(list(separated = NA, inf_terms = character(0), ok = FALSE,
                msg = conditionMessage(ds)))
  betas <- if (!is.null(ds$betas)) ds$betas
           else if (!is.null(ds$coefficients)) ds$coefficients else numeric(0)
  inf_terms <- names(betas)[is.infinite(betas)]
  sep_flag  <- if (!is.null(ds$separation)) isTRUE(ds$separation)
               else if (!is.null(ds$outcome)) isTRUE(ds$outcome)
               else length(inf_terms) > 0
  list(separated = sep_flag, inf_terms = inf_terms, ok = TRUE, msg = "")
}

glm_unstable_one <- function(f) {
  cf <- stats::coef(f); se <- suppressWarnings(sqrt(diag(stats::vcov(f)))); pr <- stats::fitted(f)
  cf_f <- cf[is.finite(cf)]; se_f <- se[is.finite(se)]
  any(!is.finite(cf)) || any(!is.finite(se)) ||
    (length(cf_f) > 0 && any(abs(cf_f) > BETA_MAX)) ||
    (length(se_f) > 0 && any(se_f > SE_MAX)) ||
    any(pr < PROB_EPS | pr > 1 - PROB_EPS) || !isTRUE(f$converged)
}

sep_or_unstable <- function(smf, imps) {
  exact_any <- FALSE; det_err <- FALSE; inf_terms <- character(0)
  for (d in imps) {
    z <- detect_one(smf, d)
    if (!isTRUE(z$ok)) { det_err <- TRUE; next }
    if (isTRUE(z$separated)) { exact_any <- TRUE; inf_terms <- union(inf_terms, z$inf_terms) }
  }
  if (exact_any)
    return(list(sep = TRUE, kind = "exact", n_unstable = NA_integer_, m = length(imps),
                inf_terms = inf_terms, det_error = det_err, max_abeta = NA_real_, max_se = NA_real_))
  n_unstable <- 0L; mab <- 0; mse <- 0
  for (d in imps) {
    f <- suppressWarnings(stats::glm(smf, binomial(), d))
    if (isTRUE(glm_unstable_one(f))) n_unstable <- n_unstable + 1L
    cf <- stats::coef(f); se <- suppressWarnings(sqrt(diag(stats::vcov(f))))
    cf <- cf[is.finite(cf)]; se <- se[is.finite(se)]
    if (length(cf)) mab <- max(mab, max(abs(cf))); if (length(se)) mse <- max(mse, max(se))
  }
  list(sep = n_unstable > 0, kind = if (n_unstable > 0) "near/unstable" else "none",
       n_unstable = n_unstable, m = length(imps), inf_terms = inf_terms,
       det_error = det_err, max_abeta = mab, max_se = mse)
}

pred_df <- function(smf, imps) {
  nc <- as.integer(round(stats::median(vapply(imps, function(d)
    length(stats::coef(suppressWarnings(glm(smf, binomial(), d)))), numeric(1)))))
  max(nc - 1L, 0L)
}

group_key <- function(lab) {
  if (grepl("_measurable$", lab)) return(lab)
  if (grepl("_c[0-9]+$", lab))    return(sub("_c[0-9]+$", "", lab))
  if (grepl("_c$", lab))          return(sub("_c$", "", lab))
  lab
}

make_formula <- function(term_labels) {
  if (length(term_labels) == 0) return(stats::as.formula("good ~ 1"))
  stats::reformulate(termlabels = term_labels, response = "good")
}

firth_pll <- function(fit) {
  ll <- fit$loglik
  if (!is.null(names(ll)) && "full" %in% names(ll)) return(unname(ll[["full"]]))
  ll[length(ll)]
}

pool_lrt_D2 <- function(stats, k) {
  stats <- pmax(as.numeric(stats), 0); stats <- stats[is.finite(stats)]
  m <- length(stats)
  if (k <= 0 || m == 0) return(list(p = NA_real_, D = NA_real_, df1 = k, df2 = NA_real_, r = NA_real_))
  mbar <- mean(stats)
  if (m == 1) return(list(p = stats::pchisq(mbar, df = k, lower.tail = FALSE),
                          D = mbar / k, df1 = k, df2 = Inf, r = 0))
  r  <- (1 + 1/m) * stats::var(sqrt(stats))
  if (!is.finite(r) || r < 0) r <- 0
  D  <- (mbar/k - (m + 1)/(m - 1) * r) / (1 + r)
  if (!is.finite(D) || D < 0) D <- 0
  df2 <- if (r <= 0) 1e6 else k^(-3/m) * (m - 1) * (1 + 1/r)^2
  if (!is.finite(df2) || df2 <= 0) df2 <- 1e6
  list(p = stats::pf(D, df1 = k, df2 = df2, lower.tail = FALSE),
       D = D, df1 = k, df2 = df2, r = r)
}

pool_lrt_D3_glm <- function(fits_full, fits_red) {
  out <- tryCatch({
    res <- mice::D3(mice::as.mira(fits_full), mice::as.mira(fits_red))
    tab <- as.data.frame(res$result)
    pcol <- grep("^P", names(tab))
    p <- if (length(pcol)) as.numeric(tab[1, pcol[1]]) else NA_real_
    list(p = p, D = as.numeric(tab[1, 1]), df1 = NA_real_, df2 = NA_real_, r = NA_real_)
  }, error = function(e) NULL, warning = function(w) NULL)
  out
}

pooled_lrt <- function(full_labels, red_labels, imps, path, method = GLM_LRT_METHOD) {
  smf_full <- make_formula(full_labels)
  smf_red  <- make_formula(red_labels)

  if (path == "glm") {
    fits_full <- lapply(imps, function(d) suppressWarnings(glm(smf_full, binomial(), d)))
    fits_red  <- lapply(imps, function(d) suppressWarnings(glm(smf_red,  binomial(), d)))
    lr <- mapply(function(ff, fr) 2 * (as.numeric(logLik(ff)) - as.numeric(logLik(fr))),
                 fits_full, fits_red)
    kvec <- mapply(function(ff, fr) length(coef(ff)) - length(coef(fr)), fits_full, fits_red)
    k <- as.integer(round(stats::median(kvec)))
    if (identical(method, "D3")) {
      d3 <- pool_lrt_D3_glm(fits_full, fits_red)
      if (!is.null(d3) && is.finite(d3$p)) return(c(d3, list(k = k, engine = "D3")))
    }
    res <- pool_lrt_D2(lr, k); return(c(res, list(k = k, engine = "D2")))

  } else {
    pll_full <- vapply(imps, function(d) {
      fit <- tryCatch(logistf::logistf(smf_full, data = d, pl = LOGISTF_PL),
                      error = function(e) NULL)
      if (is.null(fit)) NA_real_ else firth_pll(fit)
    }, numeric(1))
    pll_red <- vapply(imps, function(d) {
      fit <- tryCatch(logistf::logistf(smf_red, data = d, pl = LOGISTF_PL),
                      error = function(e) NULL)
      if (is.null(fit)) NA_real_ else firth_pll(fit)
    }, numeric(1))
    lr <- 2 * (pll_full - pll_red)
    fit1f <- tryCatch(logistf::logistf(smf_full, data = imps[[1]], pl = FALSE), error = function(e) NULL)
    fit1r <- tryCatch(logistf::logistf(smf_red,  data = imps[[1]], pl = FALSE), error = function(e) NULL)
    k <- if (!is.null(fit1f) && !is.null(fit1r))
           length(coef(fit1f)) - length(coef(fit1r)) else NA_integer_
    res <- pool_lrt_D2(lr, k); return(c(res, list(k = k, engine = "D2(Firth-PLR)")))
  }
}

pool_rubin <- function(fitlist, n_complete, ci_level = CI_LEVEL) {
  m <- length(fitlist); betas <- lapply(fitlist, coef); terms <- names(betas[[1]])
  Vdg <- lapply(fitlist, function(f) diag(as.matrix(vcov(f))))
  getv <- function(lst, t) vapply(lst, function(z) unname(z[t]), numeric(1))
  Qbar <- vapply(terms, function(t) mean(getv(betas, t)), numeric(1))
  B    <- vapply(terms, function(t) stats::var(getv(betas, t)), numeric(1))
  Ubar <- vapply(terms, function(t) mean(getv(Vdg, t)), numeric(1))
  Tvar <- Ubar + (1 + 1/m) * B
  riv  <- (1 + 1/m) * B / Ubar
  lambda <- (1 + 1/m) * B / Tvar
  k <- length(terms); dfc <- max(n_complete - k, 1)
  df_old <- ifelse(lambda > 0, (m - 1) / lambda^2, Inf)
  df_obs <- ((dfc + 1) / (dfc + 3)) * dfc * (1 - lambda)
  df  <- df_old * df_obs / (df_old + df_obs)
  fmi <- (riv + 2 / (df + 3)) / (riv + 1)
  se <- sqrt(Tvar); z <- stats::qnorm(1 - (1 - ci_level) / 2)
  data.table(term = terms, estimate = Qbar, se = se, OR = exp(Qbar),
             lcl = exp(Qbar - z * se), ucl = exp(Qbar + z * se), riv = riv, fmi = fmi)
}

estimate_glm <- function(smf, imps, n, ci_level = CI_LEVEL) {
  sep_warn <- FALSE
  fits <- lapply(imps, function(d) withCallingHandlers(
    glm(smf, family = binomial(), data = d),
    warning = function(w) {
      if (grepl("fitted probabilities", conditionMessage(w))) sep_warn <<- TRUE
      invokeRestart("muffleWarning") }))
  pooled <- mice::pool(fits)
  ps <- as.data.table(summary(pooled, conf.int = TRUE, conf.level = ci_level))
  ci_cols <- grep("%", names(ps), value = TRUE)
  if (length(ci_cols) >= 2) { lo <- ps[[ci_cols[1]]]; hi <- ps[[ci_cols[2]]] }
  else { zz <- stats::qnorm(1 - (1 - ci_level)/2); lo <- ps$estimate - zz*ps$std.error; hi <- ps$estimate + zz*ps$std.error }
  list(pooled = data.table(term = ps$term, estimate = ps$estimate, se = ps$std.error,
                           OR = exp(ps$estimate), lcl = exp(lo), ucl = exp(hi), p = ps$p.value),
       glm_sep_warn = sep_warn)
}

estimate_firth <- function(smf, imps, n, clip_vars = NULL, ci_level = CI_LEVEL) {
  fitlist <- lapply(imps, function(d) logistf::logistf(smf, data = d, pl = LOGISTF_PL, dataout = TRUE))
  rb <- pool_rubin(fitlist, n_complete = n, ci_level = ci_level)
  if (have_tic) tictoc::tic("  CLIP.confint(provisional)")
  cc <- tryCatch(logistf::CLIP.confint(obj = fitlist, variable = clip_vars,
                                       ci.level = ci_level, pvalue = TRUE),
                 error = function(e) e)
  if (have_tic) tictoc::toc()
  pooled <- copy(rb)[, .(term, estimate, se, OR, lcl, ucl, fmi)]
  clip_ok <- !inherits(cc, "error")
  if (clip_ok) {
    ci_mat <- if (!is.null(cc$ci)) as.matrix(cc$ci) else if (!is.null(cc$confint)) as.matrix(cc$confint) else NULL
    vnm <- if (!is.null(cc$variable)) cc$variable else if (!is.null(rownames(ci_mat))) rownames(ci_mat) else NULL
    if (!is.null(ci_mat) && !is.null(vnm)) {
      clip_dt <- data.table(term = vnm, lcl_clip = ci_mat[, 1], ucl_clip = ci_mat[, 2],
                            p = if (!is.null(cc$pvalue)) cc$pvalue else NA_real_)
      pooled <- merge(pooled, clip_dt, by = "term", all.x = TRUE, sort = FALSE)
      hv <- is.finite(pooled$lcl_clip) & is.finite(pooled$ucl_clip)
      pooled[hv, `:=`(lcl = exp(lcl_clip), ucl = exp(ucl_clip))]
      pooled[, c("lcl_clip", "ucl_clip") := NULL]
    } else clip_ok <- FALSE
  }
  if (!clip_ok) { warning("CLIP.confint 失敗/未対応 → CI は罰則付き Rubin 正規近似。"); pooled[, p := NA_real_] }
  list(pooled = pooled, clip_ok = clip_ok)
}

process_pipeline <- function(nm) {
  est <- readRDS(file.path(OUT_DIR, sprintf("bnb_estimation_%s.rds", nm)))
  obj <- readRDS(file.path(OUT_DIR, sprintf("bnb_imp_%s.rds",        nm)))
  odf <- readRDS(file.path(OUT_DIR, sprintf("bnb_outlier_def_%s.rds", nm)))

  imps_full   <- obj$smcfcs$impDatasets
  smf_full    <- stats::as.formula(est$smformula)
  outlier_ids <- sort(unique(odf$outlier_ids))
  n0 <- est$n; m <- length(imps_full)

  message(sprintf("\n================ [%s] §8 変数選択 (n0=%d, m=%d) ================",
                  nm, n0, m))
  message("  フルモデル smformula: ", est$smformula)

  if (length(outlier_ids) > 0) {
    bad <- outlier_ids[outlier_ids %in% seq_len(n0)]
    imps <- lapply(imps_full, function(d) d[-bad, , drop = FALSE])
    n <- n0 - length(bad)
    message(sprintf("  外れ値除外: %d 件除去（.id: %s） → n=%d",
                    length(bad), paste(bad, collapse = ", "), n))
  } else {
    imps <- imps_full; n <- n0
    message("  外れ値なし（§7-1 で 0 件）→ 全 n を使用。")
  }

  su_after <- sep_or_unstable(smf_full, imps)
  any_sep_after <- isTRUE(su_after$sep); inf_all <- su_after$inf_terms

  sep_before <- isTRUE(est$separation$any)
  change <- if (sep_before && !any_sep_after)      "解消（除外前=分離あり → 除外後=なし）"
            else if (!sep_before && any_sep_after)  "★誘発（除外前=なし → 除外後=分離相当）"
            else if (sep_before && any_sep_after)   "不変（前後とも分離相当）"
            else                                    "不変（前後とも分離なし）"
  message(sprintf("  分離 再判定（除外後フルモデル）: %s%s",
                  su_after$kind,
                  if (su_after$kind == "near/unstable")
                    sprintf("（不安定 %d/%d 代入, |β|max=%.1f, SEmax=%.1f）",
                            su_after$n_unstable, su_after$m, su_after$max_abeta, su_after$max_se)
                  else ""))
  message(sprintf("    §6 除外前=%s, §8 除外後=%s  → %s",
                  if (sep_before) "分離あり" else "分離なし",
                  if (any_sep_after) "分離相当" else "分離なし", change))
  if (any_sep_after && length(inf_all))
    message("    除外後に無限化した係数（いずれかの代入で）: ", paste(inf_all, collapse = ", "))
  if (isTRUE(su_after$det_error))
    message("    ※detect_separation が一部代入でエラー/NA → glm 直接診断で判定（握りつぶさず）。")

  forced <- FORCE_PATH[[nm]]
  base_path <- if (!is.null(forced)) {
    if (!forced %in% c("glm", "firth")) stop("FORCE_PATH は 'glm'/'firth' のみ。")
    message("  推定経路を手動指定: ", forced); forced
  } else if (any_sep_after) "firth" else "glm"
  message(sprintf("  → §8 選択の基準経路: %s（除外後の分離=%s に基づく。各ステップで逐次再判定）",
                  toupper(base_path), if (any_sep_after) "あり" else "なし"))

  all_labels <- attr(stats::terms(smf_full), "term.labels")
  keys  <- vapply(all_labels, group_key, character(1))
  groups <- split(all_labels, keys)

  keep <- KEEP_VARS[[nm]]; if (is.null(keep)) keep <- character(0)
  if (FORCE_KEEP_SUBSET_VAR && !is.null(SUBSET_VAR[[nm]]))
    keep <- union(keep, SUBSET_VAR[[nm]])
  keep <- intersect(keep, names(groups))
  if (length(keep) > 0)
    message("  選択中も保持する変数（keep）: ", paste(keep, collapse = ", "))

  sel_imps <- if (!is.null(M_SELECT_MAX) && M_SELECT_MAX < m) imps[seq_len(M_SELECT_MAX)] else imps
  if (length(sel_imps) < m)
    message(sprintf("  選択探索は m=%d 代入に間引き（最終推定は全 m=%d）。", length(sel_imps), m))

  removable0 <- setdiff(names(groups), keep)
  m_test <- if (is.numeric(BONFERRONI_M)) as.numeric(BONFERRONI_M) else length(removable0)
  m_test <- max(as.integer(round(m_test)), 1L)
  p_thresh <- ALPHA / m_test
  y_out  <- imps[[1]][["good"]]; ytab <- table(factor(y_out))
  n_events <- if (identical(EPV_EVENTS, "events")) sum(y_out == 1, na.rm = TRUE) else as.integer(min(ytab))
  epv_cap <- n_events / EPV_PER_VAR
  message(sprintf("  保持基準: Bonferroni 閾値 = ALPHA/m_test = %.3g/%d = %.3g（p>=閾値で除去）",
                  ALPHA, m_test, p_thresh))
  message(sprintf("  EPV 上限: events=%d（%s）, 上限 df=events/%g=%.1f（df 超過の間は最も非有意な変数を除去）",
                  n_events, if (identical(EPV_EVENTS, "events")) "good=1" else "min(2群)", EPV_PER_VAR, epv_cap))

  current <- names(groups)
  trace <- data.table()
  step <- 0L
  firth_locked <- (base_path == "firth")
  repeat {
    cand <- setdiff(current, keep)
    current_labels <- unlist(groups[current], use.names = FALSE)
    smf_cur <- make_formula(current_labels)
    cur_df  <- pred_df(smf_cur, sel_imps)

    if (length(cand) == 0) {
      if (ENFORCE_EPV && cur_df > epv_cap)
        message(sprintf("  EPV 未達（df=%d>上限%.1f）だが除去可能な変数なし（残りは全て keep）→ 満たせないまま停止。",
                        cur_df, epv_cap))
      else message("  これ以上落とせる変数なし（全て keep）→ 停止。")
      break
    }

    if (firth_locked) {
      step_path <- "firth"
    } else {
      cur_su <- sep_or_unstable(smf_cur, sel_imps)
      step_path <- if (isTRUE(cur_su$sep)) "firth" else "glm"
      if (isTRUE(cur_su$sep))
        message(sprintf("  [step %d] 現モデルが分離相当(%s) → 本ステップの検定は Firth%s。",
                        step + 1L, cur_su$kind, if (FIRTH_LOCK_AFTER_SWITCH) "（以降固定）" else ""))
      if (FIRTH_LOCK_AFTER_SWITCH && step_path == "firth") firth_locked <- TRUE
    }

    res_dt <- rbindlist(lapply(cand, function(g) {
      red_labels <- unlist(groups[setdiff(current, g)], use.names = FALSE)
      lr <- pooled_lrt(current_labels, red_labels, sel_imps, path = step_path)
      data.table(var = g, p = lr$p, k = lr$k, engine = lr$engine)
    }))
    res_dt <- res_dt[order(-p)]
    worst <- res_dt[1]

    epv_violated <- ENFORCE_EPV && (cur_df > epv_cap)
    message(sprintf("  [step %d] 候補LRT（p降順, %s | 現df=%d, 上限df=%.1f, 保持閾値=%.3g）:",
                    step + 1L, worst$engine, cur_df, epv_cap, p_thresh))
    print(res_dt[, .(variable = relabel_vars(var), df = k, p = signif(p, 4),
                     bonf_p = signif(pmin(p * m_test, 1), 4))])

    do_remove <- FALSE; reason <- ""
    if (epv_violated) {
      do_remove <- TRUE; reason <- sprintf("EPV(df=%d>上限%.1f)", cur_df, epv_cap)
    } else if (is.finite(worst$p) && worst$p >= p_thresh) {
      do_remove <- TRUE; reason <- sprintf("Bonferroni(p=%.4g>=%.3g)", worst$p, p_thresh)
    }

    if (do_remove) {
      step <- step + 1L
      trace <- rbind(trace, data.table(step = step, removed = worst$var,
                                       p = worst$p, bonf_p = min(worst$p * m_test, 1),
                                       df = worst$k, cur_df = cur_df, epv_cap = epv_cap,
                                       engine = worst$engine, step_path = step_path,
                                       reason = reason, n_remaining = length(current) - 1L))
      message(sprintf("    → 除去: %s [%s]", relabel_vars(worst$var), reason))
      current <- setdiff(current, worst$var)
      if (length(current) == 0) { message("  全変数除去 → 切片のみ。"); break }
    } else {
      message(sprintf("    → 停止（EPV 充足 df=%d<=%.1f かつ 最大 生p=%.4g < 保持閾値 %.3g）。暫定モデル確定。",
                      cur_df, epv_cap, worst$p, p_thresh))
      break
    }
  }

  prov_labels <- unlist(groups[current], use.names = FALSE)
  smf_prov <- make_formula(prov_labels)
  final_df <- pred_df(smf_prov, sel_imps)
  message("  暫定モデル: ", paste(deparse(smf_prov), collapse = ""))
  message(sprintf("  暫定モデルの変数: %s | 予測子df=%d, EPV上限df=%.1f（events=%d）",
                  if (length(current)) paste(relabel_vars(current), collapse = ", ") else "(切片のみ)",
                  final_df, epv_cap, n_events))
  if (ENFORCE_EPV && final_df > epv_cap)
    message("  ※EPV 上限を満たせませんでした（keep 変数の df が上限超過の可能性）。")

  su_final <- sep_or_unstable(smf_prov, imps)
  any_sep_final <- isTRUE(su_final$sep)
  path <- if (!is.null(forced)) forced else if (any_sep_final) "firth" else "glm"
  message(sprintf("  暫定モデル 分離再判定: %s → 最終推定の経路: %s（選択基準経路=%s）",
                  su_final$kind, toupper(path), toupper(base_path)))

  sel_desc <- sprintf("backward selection; retain rule: Bonferroni alpha/m_test=%.3g/%d=%.3g; EPV cap df<=events/%g=%.1f; per-step adaptive Firth on (quasi-)separation",
                      ALPHA, m_test, p_thresh, EPV_PER_VAR, epv_cap)
  if (have_tic) tictoc::tic(paste0(nm, " final estimate"))
  if (path == "glm") {
    fin <- estimate_glm(smf_prov, imps, n, ci_level = CI_LEVEL)
    pooled <- fin$pooled
    pooled_unstable <- any(!is.finite(pooled$estimate)) || any(abs(pooled$estimate) > BETA_MAX) ||
                       any(!is.finite(pooled$se)) || any(pooled$se > SE_MAX) || isTRUE(fin$glm_sep_warn)
    if (pooled_unstable && AUTO_FALLBACK_FIRTH && is.null(forced)) {
      message(sprintf("  ※glm 最終推定が不安定（|β|max=%.1f, SEmax=%.1f%s）→ Firth に自動フォールバック。",
                      max(abs(pooled$estimate[is.finite(pooled$estimate)]), 0),
                      max(pooled$se[is.finite(pooled$se)], 0),
                      if (isTRUE(fin$glm_sep_warn)) ", fitted 0/1 警告" else ""))
      path <- "firth"; any_sep_final <- TRUE
      fin <- estimate_firth(smf_prov, imps, n, clip_vars = CLIP_VARS[[nm]], ci_level = CI_LEVEL)
      pooled <- fin$pooled
      method_detail <- list(estimator = "logistf (Firth)", pool = "Rubin (point) + CLIP.confint (CI/p)",
                            selection = sel_desc, clip_ok = fin$clip_ok, auto_fallback = TRUE)
    } else {
      method_detail <- list(estimator = "glm", pool = "mice::pool (Rubin)", selection = sel_desc)
    }
  } else {
    fin <- estimate_firth(smf_prov, imps, n, clip_vars = CLIP_VARS[[nm]], ci_level = CI_LEVEL)
    pooled <- fin$pooled
    method_detail <- list(estimator = "logistf (Firth)", pool = "Rubin (point) + CLIP.confint (CI/p)",
                          selection = sel_desc, clip_ok = fin$clip_ok)
  }
  if (have_tic) tictoc::toc()

  show <- copy(pooled)[, `:=`(OR = round(OR, 3), lcl = round(lcl, 3), ucl = round(ucl, 3))]
  if ("p" %in% names(show)) show[, p := signif(p, 4)]
  show[, label := relabel_vars(sub("(運動器|廃用|有|無|女|男|あり|なし|impaired|Yes|No|Female|Male|Needed|Independent)$", "", term))]
  message("  暫定モデル プール推定:")
  print(show[, intersect(c("term", "label", "OR", "lcl", "ucl", "p", "fmi"), names(show)), with = FALSE])

  saveRDS(list(
    pipeline       = nm,
    path           = path,
    smformula_full = est$smformula,
    smformula_prov = paste(deparse(smf_prov), collapse = ""),
    prov_vars      = current,
    prov_labels    = prov_labels,
    groups         = groups,
    keep_vars      = keep,
    selection      = list(method = "backward", direction = "decrease",
                          retain_rule = "bonferroni", alpha = ALPHA,
                          m_test = m_test, m_test_basis = if (is.numeric(BONFERRONI_M)) "numeric" else BONFERRONI_M,
                          p_threshold = p_thresh,
                          enforce_epv = ENFORCE_EPV, epv_per_var = EPV_PER_VAR,
                          n_events = n_events, epv_cap = epv_cap, final_df = final_df,
                          sequential_separation = TRUE,
                          firth_lock_after_switch = FIRTH_LOCK_AFTER_SWITCH,
                          lrt = if (path == "glm") GLM_LRT_METHOD else "D2",
                          trace = trace,
                          note = "Overfitting control: Bonferroni retain threshold + EPV cap (df<=events/10) + per-step adaptive Firth on (quasi-)separation."),
    pooled         = pooled,
    method_detail  = method_detail,
    outlier_ids    = outlier_ids,
    n0 = n0, n = n, m = m,
    base_path = base_path,
    separation_after = list(any = any_sep_after, kind = su_after$kind,
                            n_unstable = su_after$n_unstable, m = m, inf_terms = inf_all,
                            max_abeta = su_after$max_abeta, max_se = su_after$max_se,
                            det_error = su_after$det_error,
                            before = sep_before, change = change),
    separation_final = list(any = any_sep_final, kind = su_final$kind,
                            n_unstable = su_final$n_unstable, inf_terms = su_final$inf_terms,
                            auto_fallback = isTRUE(method_detail$auto_fallback)),
    mnar_system = obj$mnar_system,
    selected_cont = obj$selected_cont, selected_cat = obj$selected_cat,
    rcs_cont = obj$rcs_cont, rcs_knots = obj$rcs_knots, n_knots = obj$n_knots,
    centers = obj$centers, mnar_indicator = obj$mnar_indicator,
    basis_cols = obj$basis_cols
  ), file.path(OUT_DIR, sprintf("bnb_provisional_%s.rds", nm)))

  message(sprintf("[%s] §8 完了。暫定モデル(%d変数)を bnb_provisional_%s.rds に保存。",
                  nm, length(current), nm))

  data.table(pipeline = nm, mnar_system = obj$mnar_system, n0 = n0, n = n, m = m,
             n_excluded = length(outlier_ids),
             sep_before = sep_before, sep_after = any_sep_after, sep_change = change,
             path = path, n_prov_vars = length(current),
             m_test = m_test, p_threshold = signif(p_thresh, 3),
             n_events = n_events, epv_cap = round(epv_cap, 1), final_df = final_df,
             prov_vars = paste(current, collapse = " + "))
}

prov_index <- rbindlist(lapply(est_index$pipeline, process_pipeline))
saveRDS(prov_index, file.path(OUT_DIR, "bnb_provisional_index.rds"))

message("\n==== §8 サマリ（暫定モデル・分離再判定）====")
print(prov_index[, .(pipeline, n0, n, n_excluded, sep_before, sep_after, path, n_prov_vars,
                     m_test, p_threshold, n_events, epv_cap, final_df)])

induced <- prov_index[sep_change %like% "誘発"]
if (nrow(induced) > 0) {
  message("\n★外れ値除外により分離が誘発され、§8 で Firth 経路に切り替えたパイプライン:")
  print(induced[, .(pipeline, sep_before, sep_after, path)])
} else {
  message("\n外れ値除外による分離の新規誘発はありませんでした（前後で経路の切替なし or 既存の分離のみ）。")
}

message("\n§8 完了。各パイプラインの暫定モデルを bnb_provisional_<pipeline>.rds に保存しました。",
        " §9 はこの prov_vars/groups から任意2変数の交互作用を作り、プールLRT(p<.05) で選抜します。",
        " 経路フラグ(path)は『外れ値除外後』の分離再判定に基づく点に注意（§6 とは異なりうる）。")
