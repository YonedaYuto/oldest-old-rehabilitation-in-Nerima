library(data.table)
library(here)

need_pkg <- function(p) if (!requireNamespace(p, quietly = TRUE))
  stop(sprintf("パッケージ '%s' が必要です。install.packages('%s')", p, p))
need_pkg("mice")
need_pkg("logistf")

HAS_DETECTSEP <- requireNamespace("detectseparation", quietly = TRUE)
HAS_BRGLM2    <- requireNamespace("brglm2",          quietly = TRUE)
if (!HAS_DETECTSEP && !HAS_BRGLM2)
  stop("分離検出に detectseparation か brglm2 のいずれかが必要です。install.packages('detectseparation')")

LABELS_PATH <- here::here("00_labels.R")
if (!exists("relabel_vars")) {
  if (!file.exists(LABELS_PATH)) stop("00_labels.R が見つかりません: ", LABELS_PATH)
  source(LABELS_PATH)
}

OUT_DIR <- here::here("data")
have_tic <- requireNamespace("tictoc", quietly = TRUE)

CI_LEVEL   <- 0.95
GLM_LRT_METHOD <- "D3"
RECHECK_SEP_PER_INTERACTION <- TRUE

FIRTH_LOCK_AFTER_SWITCH <- TRUE

ALPHA <- 0.05
BONFERRONI_M <- "candidates"

EPV_PER_VAR <- 10
EPV_EVENTS  <- "min"
ENFORCE_EPV <- TRUE

FORCE_PATH <- list(
  severe_A  = NULL, severe_B  = NULL,
  elderly_A = NULL, elderly_B = NULL
)

M_SELECT_MAX <- NULL

CLIP_VARS <- list(
  severe_A  = NULL, severe_B  = NULL,
  elderly_A = NULL, elderly_B = NULL
)
LOGISTF_PL <- FALSE

BETA_MAX <- 15
SE_MAX   <- 30
PROB_EPS <- 1e-6
AUTO_FALLBACK_FIRTH <- TRUE

if (!exists("int_index"))
  int_index <- readRDS(file.path(OUT_DIR, "bnb_interactions_index.rds"))

detect_one <- function(smf, dat) {
  detector <- if (HAS_DETECTSEP) detectseparation::detect_separation
              else                brglm2::detect_separation
  ds <- tryCatch(stats::glm(smf, data = dat, family = binomial("logit"), method = detector),
                 error = function(e) e)
  if (inherits(ds, "error")) return(list(separated = NA, inf_terms = character(0), ok = FALSE))
  betas <- if (!is.null(ds$betas)) ds$betas
           else if (!is.null(ds$coefficients)) ds$coefficients else numeric(0)
  inf_terms <- names(betas)[is.infinite(betas)]
  sep_flag  <- if (!is.null(ds$separation)) isTRUE(ds$separation)
               else if (!is.null(ds$outcome)) isTRUE(ds$outcome)
               else length(inf_terms) > 0
  list(separated = sep_flag, inf_terms = inf_terms, ok = TRUE)
}

glm_unstable_one <- function(f) {
  cf <- stats::coef(f); se <- suppressWarnings(sqrt(diag(stats::vcov(f)))); pr <- stats::fitted(f)
  cf_f <- cf[is.finite(cf)]; se_f <- se[is.finite(se)]
  bad <- any(!is.finite(cf)) || any(!is.finite(se)) ||
         (length(cf_f) > 0 && any(abs(cf_f) > BETA_MAX)) ||
         (length(se_f) > 0 && any(se_f > SE_MAX)) ||
         any(pr < PROB_EPS | pr > 1 - PROB_EPS) || !isTRUE(f$converged)
  list(bad = bad,
       max_abeta = if (length(cf_f)) max(abs(cf_f)) else NA_real_,
       max_se    = if (length(se_f)) max(se_f) else NA_real_)
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
    u <- glm_unstable_one(f)
    if (isTRUE(u$bad)) n_unstable <- n_unstable + 1L
    if (is.finite(u$max_abeta)) mab <- max(mab, u$max_abeta)
    if (is.finite(u$max_se))    mse <- max(mse, u$max_se)
  }
  list(sep = n_unstable > 0, kind = if (n_unstable > 0) "near/unstable" else "none",
       n_unstable = n_unstable, m = length(imps), inf_terms = inf_terms,
       det_error = det_err, max_abeta = mab, max_se = mse)
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
  if (is.na(k) || k <= 0 || m == 0) return(list(p = NA_real_, D = NA_real_, df1 = k, df2 = NA_real_, r = NA_real_))
  mbar <- mean(stats)
  if (m == 1) return(list(p = stats::pchisq(mbar, df = k, lower.tail = FALSE), D = mbar/k, df1 = k, df2 = Inf, r = 0))
  r <- (1 + 1/m) * stats::var(sqrt(stats)); if (!is.finite(r) || r < 0) r <- 0
  D <- (mbar/k - (m + 1)/(m - 1) * r) / (1 + r); if (!is.finite(D) || D < 0) D <- 0
  df2 <- if (r <= 0) 1e6 else k^(-3/m) * (m - 1) * (1 + 1/r)^2
  if (!is.finite(df2) || df2 <= 0) df2 <- 1e6
  list(p = stats::pf(D, df1 = k, df2 = df2, lower.tail = FALSE), D = D, df1 = k, df2 = df2, r = r)
}

pool_lrt_D3_glm <- function(fits_full, fits_red) {
  tryCatch({
    res <- mice::D3(mice::as.mira(fits_full), mice::as.mira(fits_red))
    tab <- as.data.frame(res$result); pcol <- grep("^P", names(tab))
    list(p = if (length(pcol)) as.numeric(tab[1, pcol[1]]) else NA_real_,
         D = as.numeric(tab[1, 1]), df1 = NA_real_, df2 = NA_real_, r = NA_real_)
  }, error = function(e) NULL, warning = function(w) NULL)
}

pooled_lrt <- function(full_labels, red_labels, imps, path, method = GLM_LRT_METHOD) {
  smf_full <- make_formula(full_labels); smf_red <- make_formula(red_labels)
  if (path == "glm") {
    fits_full <- lapply(imps, function(d) suppressWarnings(glm(smf_full, binomial(), d)))
    fits_red  <- lapply(imps, function(d) suppressWarnings(glm(smf_red,  binomial(), d)))
    lr <- mapply(function(ff, fr) 2 * (as.numeric(logLik(ff)) - as.numeric(logLik(fr))), fits_full, fits_red)
    kvec <- mapply(function(ff, fr) length(coef(ff)) - length(coef(fr)), fits_full, fits_red)
    k <- as.integer(round(stats::median(kvec)))
    if (identical(method, "D3")) {
      d3 <- pool_lrt_D3_glm(fits_full, fits_red)
      if (!is.null(d3) && is.finite(d3$p)) return(c(d3, list(k = k, engine = "D3")))
    }
    return(c(pool_lrt_D2(lr, k), list(k = k, engine = "D2")))
  } else {
    pll_full <- vapply(imps, function(d) { fit <- tryCatch(logistf::logistf(smf_full, data = d, pl = LOGISTF_PL), error = function(e) NULL)
      if (is.null(fit)) NA_real_ else firth_pll(fit) }, numeric(1))
    pll_red  <- vapply(imps, function(d) { fit <- tryCatch(logistf::logistf(smf_red, data = d, pl = LOGISTF_PL), error = function(e) NULL)
      if (is.null(fit)) NA_real_ else firth_pll(fit) }, numeric(1))
    lr <- 2 * (pll_full - pll_red)
    f1f <- tryCatch(logistf::logistf(smf_full, data = imps[[1]], pl = FALSE), error = function(e) NULL)
    f1r <- tryCatch(logistf::logistf(smf_red,  data = imps[[1]], pl = FALSE), error = function(e) NULL)
    k <- if (!is.null(f1f) && !is.null(f1r)) length(coef(f1f)) - length(coef(f1r)) else NA_integer_
    return(c(pool_lrt_D2(lr, k), list(k = k, engine = "D2(Firth-PLR)")))
  }
}

pool_rubin <- function(fitlist, n_complete, ci_level = CI_LEVEL) {
  m <- length(fitlist); betas <- lapply(fitlist, coef); terms <- names(betas[[1]])
  Vdg <- lapply(fitlist, function(f) diag(as.matrix(vcov(f))))
  getv <- function(lst, t) vapply(lst, function(z) unname(z[t]), numeric(1))
  Qbar <- vapply(terms, function(t) mean(getv(betas, t)), numeric(1))
  B    <- vapply(terms, function(t) stats::var(getv(betas, t)), numeric(1))
  Ubar <- vapply(terms, function(t) mean(getv(Vdg, t)), numeric(1))
  Tvar <- Ubar + (1 + 1/m) * B; riv <- (1 + 1/m) * B / Ubar; lambda <- (1 + 1/m) * B / Tvar
  k <- length(terms); dfc <- max(n_complete - k, 1)
  df_old <- ifelse(lambda > 0, (m - 1) / lambda^2, Inf)
  df_obs <- ((dfc + 1)/(dfc + 3)) * dfc * (1 - lambda)
  df <- df_old * df_obs / (df_old + df_obs); fmi <- (riv + 2/(df + 3)) / (riv + 1)
  se <- sqrt(Tvar); z <- stats::qnorm(1 - (1 - ci_level)/2)
  data.table(term = terms, estimate = Qbar, se = se, OR = exp(Qbar),
             lcl = exp(Qbar - z*se), ucl = exp(Qbar + z*se), riv = riv, fmi = fmi)
}

estimate_glm <- function(smf, imps, n, ci_level = CI_LEVEL) {
  sep_warn <- FALSE
  fits <- lapply(imps, function(d) withCallingHandlers(glm(smf, binomial(), d),
    warning = function(w) { if (grepl("fitted probabilities", conditionMessage(w))) sep_warn <<- TRUE
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
  if (have_tic) tictoc::tic("  CLIP.confint(final)")
  cc <- tryCatch(logistf::CLIP.confint(obj = fitlist, variable = clip_vars, ci.level = ci_level, pvalue = TRUE),
                 error = function(e) e)
  if (have_tic) tictoc::toc()
  pooled <- copy(rb)[, .(term, estimate, se, OR, lcl, ucl, fmi)]; clip_ok <- !inherits(cc, "error")
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

inter_disp <- function(key) {
  pr <- strsplit(key, ":", fixed = TRUE)[[1]]
  paste(relabel_vars(pr[1]), "x", relabel_vars(pr[2]))
}

process_pipeline <- function(nm) {
  ix  <- readRDS(file.path(OUT_DIR, sprintf("bnb_interactions_%s.rds", nm)))
  obj <- readRDS(file.path(OUT_DIR, sprintf("bnb_imp_%s.rds",          nm)))

  prov_labels     <- ix$prov_labels
  groups          <- ix$groups
  kept_int_groups <- ix$kept_int_groups
  outlier_ids     <- ix$outlier_ids
  n0 <- ix$n0; m <- length(obj$smcfcs$impDatasets)

  base_path <- ix$path

  message(sprintf("\n================ [%s] §10 最終モデル（前進選択）================", nm))
  message("  主効果（固定）: ",
          if (length(ix$prov_vars)) paste(relabel_vars(ix$prov_vars), collapse = ", ") else "(切片のみ)")
  message(sprintf("  §9 残存交互作用（候補）: %d 件%s",
                  length(kept_int_groups),
                  if (length(kept_int_groups))
                    paste0("（", paste(vapply(names(kept_int_groups), inter_disp, character(1)), collapse = ", "), "）")
                  else ""))

  bad <- sort(unique(outlier_ids)); bad <- bad[bad %in% seq_len(n0)]
  imps <- if (length(bad) > 0) lapply(obj$smcfcs$impDatasets, function(d) d[-bad, , drop = FALSE])
          else obj$smcfcs$impDatasets
  n <- n0 - length(bad)
  message(sprintf("  縮小データ: n=%d（除外 %d 件）, m=%d, 選択基準経路=%s",
                  n, length(bad), m, toupper(base_path)))

  sel_imps <- if (!is.null(M_SELECT_MAX) && M_SELECT_MAX < m) imps[seq_len(M_SELECT_MAX)] else imps

  n_cand <- length(kept_int_groups)
  m_test <- if (is.numeric(BONFERRONI_M)) as.numeric(BONFERRONI_M)
            else if (identical(BONFERRONI_M, "pairs")) {
              np <- tryCatch(nrow(ix$candidates), error = function(e) NA_integer_)
              if (is.null(np) || is.na(np) || np < 1) n_cand else np
            } else n_cand
  m_test <- max(as.integer(round(m_test)), 1L)
  p_thresh <- ALPHA / m_test
  message(sprintf("  投入基準: Bonferroni 調整 閾値 = ALPHA/m_test = %.3g/%d = %.3g（m_test=%s）",
                  ALPHA, m_test, p_thresh,
                  if (is.numeric(BONFERRONI_M)) "指定値" else BONFERRONI_M))

  y_out  <- imps[[1]][["good"]]
  ytab   <- table(factor(y_out))
  n_events <- if (identical(EPV_EVENTS, "events")) sum(y_out == 1, na.rm = TRUE)
              else as.integer(min(ytab))
  epv_cap <- n_events / EPV_PER_VAR
  main_df <- as.integer(round(stats::median(vapply(sel_imps, function(d)
               length(stats::coef(suppressWarnings(glm(make_formula(prov_labels), binomial(), d)))),
               numeric(1))))) - 1L
  main_df <- max(main_df, 0L)
  current_df <- main_df
  message(sprintf("  EPV 上限: events=%d（%s）, 上限 df=events/%g=%.1f / 主効果 df=%d → 交互作用に使える df=%.1f",
                  n_events, if (identical(EPV_EVENTS, "events")) "good=1" else "min(2群)",
                  EPV_PER_VAR, epv_cap, main_df, max(epv_cap - main_df, 0)))
  if (ENFORCE_EPV && main_df > epv_cap)
    message("  ※主効果のみで既に EPV 上限を超過 → 交互作用は追加不可（0 件）。")

  current <- prov_labels
  remaining <- names(kept_int_groups)
  selected_keys <- character(0)
  trace <- data.table()
  step <- 0L

  firth_locked <- (base_path == "firth")
  if (firth_locked)
    message("  §9 経路が Firth → 最初から Firth 固定（分離判定をスキップ、計算量削減）。")

  if (length(remaining) == 0) {
    message("  §9 残存交互作用が 0 → 最終モデル = 暫定モデル（主効果のみ）。")
  } else repeat {
    if (length(remaining) == 0) { message("  全候補を投入 → 停止。"); break }

    if (!firth_locked) {
      cur_su <- sep_or_unstable(make_formula(current), sel_imps)
      if (isTRUE(cur_su$sep)) {
        if (FIRTH_LOCK_AFTER_SWITCH) firth_locked <- TRUE
        message(sprintf("  [step %d] 現モデルが分離相当(%s) → 本ステップ以降は Firth%s。",
                        step + 1L, cur_su$kind,
                        if (FIRTH_LOCK_AFTER_SWITCH) " 固定（判定スキップ）" else " で評価"))
      }
    }

    res_dt <- rbindlist(lapply(remaining, function(g) {
      int_labels  <- kept_int_groups[[g]]
      full_labels <- c(current, int_labels)
      if (firth_locked) {
        cand_path <- "firth"; sep_induced <- TRUE
      } else if (RECHECK_SEP_PER_INTERACTION) {
        aug_su <- sep_or_unstable(make_formula(full_labels), sel_imps)
        cand_path <- if (isTRUE(aug_su$sep)) "firth" else "glm"
        sep_induced <- cand_path == "firth"
      } else {
        cand_path <- "glm"; sep_induced <- FALSE
      }
      lr <- pooled_lrt(full_labels, current, sel_imps, path = cand_path)
      p_adj <- if (is.finite(lr$p)) min(lr$p * m_test, 1) else NA_real_
      data.table(ikey = g, df = lr$k, engine = lr$engine, sep_induced = sep_induced,
                 p = lr$p, p_adj = p_adj)
    }))
    res_dt <- res_dt[order(p)]
    res_dt[, eligible := if (ENFORCE_EPV) (current_df + df <= epv_cap) else TRUE]

    message(sprintf("  [step %d] 未投入候補LRT（生p昇順, 閾値=%.3g, Bonferroni m_test=%d, %s | 現df=%d, 上限df=%.1f）:",
                    step + 1L, p_thresh, m_test,
                    if (firth_locked) "Firth 固定" else "glm/Firth 逐次判定", current_df, epv_cap))
    print(res_dt[, .(interaction = vapply(ikey, inter_disp, character(1)), df, engine, sep_induced,
                     p = signif(p, 4), p_bonf = signif(p_adj, 4), eligible)])

    elig <- res_dt[eligible == TRUE]
    add_now <- nrow(elig) > 0 && is.finite(elig[1]$p) && elig[1]$p < p_thresh

    if (add_now) {
      best <- elig[1]
      step <- step + 1L
      current <- c(current, kept_int_groups[[best$ikey]])
      current_df <- current_df + best$df
      selected_keys <- c(selected_keys, best$ikey)
      remaining <- setdiff(remaining, best$ikey)
      newly_locked <- FALSE
      if (!firth_locked && FIRTH_LOCK_AFTER_SWITCH && isTRUE(best$sep_induced)) {
        firth_locked <- TRUE; newly_locked <- TRUE
      }
      trace <- rbind(trace, data.table(step = step, added = best$ikey, interaction = inter_disp(best$ikey),
                                       p = best$p, p_adj = best$p_adj, p_thresh = p_thresh,
                                       df = best$df, total_df = current_df, epv_cap = epv_cap,
                                       engine = best$engine,
                                       sep_induced = best$sep_induced, firth_locked = firth_locked,
                                       n_int = length(selected_keys)))
      message(sprintf("    → 追加: %s (生p=%.4g < %.3g; Bonferroni p=%.4g < %.2f; df+%d→総df=%d/上限%.1f)%s%s",
                      inter_disp(best$ikey), best$p, p_thresh, best$p_adj, ALPHA,
                      best$df, current_df, epv_cap,
                      if (best$sep_induced) " ［分離相当→Firth評価］" else "",
                      if (newly_locked) " → 以降 Firth 固定（計算量削減）" else ""))
    } else {
      sig_inelig <- res_dt[eligible == FALSE & is.finite(p) & p < p_thresh]
      if (nrow(elig) == 0) {
        message(sprintf("    → 停止（EPV 上限により追加可能な候補なし; 残り df 予算=%.1f）。",
                        max(epv_cap - current_df, 0)))
      } else {
        message(sprintf("    → 停止（適格候補の最小 生p=%.4g >= 閾値 %.3g; Bonferroni p=%.4g >= %.2f）。",
                        elig[1]$p, p_thresh, elig[1]$p_adj, ALPHA))
      }
      if (nrow(sig_inelig) > 0)
        message(sprintf("      （有意だが EPV 超過で不採用: %s）",
                        paste(vapply(sig_inelig$ikey, inter_disp, character(1)), collapse = ", ")))
      break
    }
  }

  final_labels <- current
  final_df <- current_df
  selected_int_labels <- if (length(selected_keys)) unlist(kept_int_groups[selected_keys], use.names = FALSE) else character(0)
  smf_final <- make_formula(final_labels)
  message("  最終モデル: ", paste(deparse(smf_final), collapse = ""))
  message(sprintf("  追加された交互作用: %d / %d 件 | 最終 予測子df=%d（主効果%d＋交互作用%d）, EPV上限df=%.1f（events=%d）",
                  length(selected_keys), length(kept_int_groups),
                  final_df, main_df, final_df - main_df, epv_cap, n_events))

  su <- sep_or_unstable(smf_final, imps)
  any_sep_final <- isTRUE(su$sep); inf_all <- su$inf_terms
  if (su$kind == "exact") {
    message(sprintf("  最終モデル 判定: 厳密(準)完全分離あり | 無限化: %s",
                    paste(inf_all, collapse = ", ")))
  } else if (su$kind == "near/unstable") {
    message(sprintf("  最終モデル 判定: near-separation/不安定 %d/%d 代入（|β|max=%.1f, SEmax=%.1f, |β|>%g or SE>%g or fitted∉(%g,1-%g) or 非収束）",
                    su$n_unstable, su$m, su$max_abeta, su$max_se, BETA_MAX, SE_MAX, PROB_EPS, PROB_EPS))
  } else {
    message(sprintf("  最終モデル 判定: 分離なし・glm 安定（|β|max=%.1f, SEmax=%.1f）",
                    su$max_abeta, su$max_se))
  }
  if (isTRUE(su$det_error))
    message("  ※detect_separation が一部代入でエラー/NA → glm 直接診断で判定（握りつぶさず）。")

  forced <- FORCE_PATH[[nm]]
  path <- if (!is.null(forced)) { if (!forced %in% c("glm","firth")) stop("FORCE_PATH は 'glm'/'firth'。"); forced }
          else if (any_sep_final) "firth" else "glm"
  message(sprintf("  → §10 最終推定の経路: %s（§8/§9 base=%s, §9 path_hint=%s）",
                  toupper(path), toupper(base_path), toupper(ix$path_hint)))

  sel_desc <- sprintf("forward interaction selection (main fixed); entry: Bonferroni alpha/m_test=%.3g/%d=%.3g; per-step adaptive Firth on (quasi-)separation",
                      ALPHA, m_test, p_thresh)
  if (have_tic) tictoc::tic(paste0(nm, " final estimate"))
  if (path == "glm") {
    fin <- estimate_glm(smf_final, imps, n, ci_level = CI_LEVEL); pooled <- fin$pooled
    pooled_unstable <- any(!is.finite(pooled$estimate)) || any(abs(pooled$estimate) > BETA_MAX) ||
                       any(!is.finite(pooled$se)) || any(pooled$se > SE_MAX) || isTRUE(fin$glm_sep_warn)
    if (pooled_unstable && AUTO_FALLBACK_FIRTH && is.null(FORCE_PATH[[nm]])) {
      message(sprintf("  ※glm 最終推定が不安定（|β|max=%.1f, SEmax=%.1f%s）→ Firth に自動フォールバック。",
                      max(abs(pooled$estimate[is.finite(pooled$estimate)]), 0),
                      max(pooled$se[is.finite(pooled$se)], 0),
                      if (isTRUE(fin$glm_sep_warn)) ", fitted 0/1 警告" else ""))
      path <- "firth"; any_sep_final <- TRUE
      if (identical(su$kind, "none")) su$kind <- "near/unstable(pooled)"
      su$sep <- TRUE
      fin <- estimate_firth(smf_final, imps, n, clip_vars = CLIP_VARS[[nm]], ci_level = CI_LEVEL)
      pooled <- fin$pooled
      method_detail <- list(estimator = "logistf (Firth)", pool = "Rubin (point) + CLIP.confint (CI/p)",
                            selection = sel_desc,
                            clip_ok = fin$clip_ok, auto_fallback = TRUE)
    } else {
      if (isTRUE(fin$glm_sep_warn))
        message("  ※glm で fitted prob 0/1 警告（準完全分離の兆候）。FORCE_PATH='firth' で感度確認を推奨。")
      method_detail <- list(estimator = "glm", pool = "mice::pool (Rubin)",
                            selection = sel_desc)
    }
  } else {
    fin <- estimate_firth(smf_final, imps, n, clip_vars = CLIP_VARS[[nm]], ci_level = CI_LEVEL); pooled <- fin$pooled
    method_detail <- list(estimator = "logistf (Firth)", pool = "Rubin (point) + CLIP.confint (CI/p)",
                          selection = sel_desc,
                          clip_ok = fin$clip_ok)
  }
  if (have_tic) tictoc::toc()

  show <- copy(pooled)[, `:=`(OR = round(OR, 3), lcl = round(lcl, 3), ucl = round(ucl, 3))]
  if ("p" %in% names(show)) show[, p := signif(p, 4)]
  message("  最終モデル プール推定:")
  print(show[, intersect(c("term", "OR", "lcl", "ucl", "p", "fmi"), names(show)), with = FALSE])

  saveRDS(list(
    pipeline        = nm,
    path            = path,
    smformula_final = paste(deparse(smf_final), collapse = ""),
    smformula_prov  = ix$smformula_prov,
    main_labels     = prov_labels,
    main_vars       = ix$prov_vars,
    selected_int_keys   = selected_keys,
    selected_int_labels = selected_int_labels,
    selected_int_groups = if (length(selected_keys)) setNames(kept_int_groups[selected_keys], selected_keys) else list(),
    candidate_int_keys  = names(kept_int_groups),
    final_labels    = final_labels,
    groups          = groups,
    selection       = list(method = "forward", direction = "increase",
                           main_effects_fixed = TRUE,
                           entry_rule = "bonferroni", alpha = ALPHA,
                           m_test = m_test, m_test_basis = if (is.numeric(BONFERRONI_M)) "numeric" else BONFERRONI_M,
                           p_threshold = p_thresh,
                           sequential_separation = TRUE,
                           firth_lock_after_switch = FIRTH_LOCK_AFTER_SWITCH,
                           enforce_epv = ENFORCE_EPV, epv_per_var = EPV_PER_VAR,
                           n_events = n_events, epv_cap = epv_cap,
                           main_df = main_df, final_df = final_df, int_df = final_df - main_df,
                           lrt = if (path == "glm") GLM_LRT_METHOD else "D2",
                           trace = trace,
                           note = "Plan §10 backward → forward (user); overfitting control: Bonferroni entry + per-step adaptive Firth on (quasi-)separation + EPV cap (main+interactions df <= events/10)."),
    pooled          = pooled,
    method_detail   = method_detail,
    separation_final = list(any = any_sep_final, kind = su$kind, n_unstable = su$n_unstable,
                            max_abeta = su$max_abeta, max_se = su$max_se,
                            m = m, inf_terms = inf_all, det_error = su$det_error,
                            beta_max = BETA_MAX, se_max = SE_MAX, prob_eps = PROB_EPS,
                            auto_fallback = isTRUE(method_detail$auto_fallback)),
    outlier_ids = outlier_ids, n0 = n0, n = n, m = m,
    mnar_system = ix$mnar_system,
    selected_cont = ix$selected_cont, selected_cat = ix$selected_cat,
    rcs_cont = ix$rcs_cont, rcs_knots = ix$rcs_knots, n_knots = ix$n_knots,
    centers = ix$centers, mnar_indicator = ix$mnar_indicator, basis_cols = ix$basis_cols
  ), file.path(OUT_DIR, sprintf("bnb_final_%s.rds", nm)))

  message(sprintf("[%s] §10 完了。最終モデル（主効果%d＋交互作用%d）を bnb_final_%s.rds に保存。",
                  nm, length(ix$prov_vars), length(selected_keys), nm))

  data.table(pipeline = nm, mnar_system = ix$mnar_system, n = n, m = m,
             path = path, sep_final = any_sep_final,
             n_main = length(ix$prov_vars),
             n_int_candidate = length(kept_int_groups), n_int_selected = length(selected_keys),
             m_test = m_test, p_threshold = signif(p_thresh, 3),
             n_events = n_events, epv_cap = round(epv_cap, 1),
             main_df = main_df, final_df = final_df,
             selected_int = if (length(selected_keys)) paste(selected_keys, collapse = ", ") else "")
}

final_index <- rbindlist(lapply(int_index$pipeline, process_pipeline), fill = TRUE)
saveRDS(final_index, file.path(OUT_DIR, "bnb_final_index.rds"))

message("\n==== §10 サマリ（最終モデル・前進選択）====")
print(final_index[, .(pipeline, n, path, n_main, n_int_candidate, n_int_selected,
                      m_test, p_threshold, n_events, epv_cap, main_df, final_df, selected_int)])

message("\n§10 完了。各パイプラインの最終モデルを bnb_final_<pipeline>.rds に保存しました。",
        " 方式: 主効果固定＋交互作用の前進(変数増加)選択（プラン §10 の後退から変更）。",
        " 過適合対策として、(1) 各ステップで(準)分離・near-separation を逐次判定し分離相当の候補は",
        " Firth(罰則付きLRT)で評価、(2) 投入基準を Bonferroni 調整閾値 ALPHA/m_test に変更、",
        " (3) EPV 制約（主効果＋追加交互作用の予測子df ≤ events/", EPV_PER_VAR, "）を適用。",
        " §11 は暫定/最終モデルの変数-OR図、§12 は性能評価、§13 は感度統合に進みます。")
