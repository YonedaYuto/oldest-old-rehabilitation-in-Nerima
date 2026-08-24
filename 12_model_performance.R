library(data.table)
library(here)
library(ggplot2)

need_pkg <- function(p) if (!requireNamespace(p, quietly = TRUE))
  stop(sprintf("パッケージ '%s' が必要です。install.packages('%s')", p, p))
need_pkg("patchwork")
HAS_LOGISTF <- requireNamespace("logistf",         quietly = TRUE)
HAS_PROC    <- requireNamespace("pROC",            quietly = TRUE)
HAS_HL      <- requireNamespace("ResourceSelection", quietly = TRUE)
HAS_CAR     <- requireNamespace("car",             quietly = TRUE)
have_tic    <- requireNamespace("tictoc",          quietly = TRUE)

LABELS_PATH <- here::here("00_labels.R")
if (!exists("relabel_vars")) {
  if (!file.exists(LABELS_PATH)) stop("00_labels.R が見つかりません: ", LABELS_PATH)
  source(LABELS_PATH)
}

OUT_DIR <- here::here("data")
FIG_DIR <- here::here("figures")
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

CI_LEVEL <- 0.95
HL_GROUPS <- 10L
CAL_BINS  <- 10L
N_GRID    <- 100L
PROB_EPS  <- 1e-8
PIPELINES <- names(PIPE_LABELS)

DO_BOOTSTRAP <- TRUE
N_BOOT       <- 200L
M_BOOT       <- 10L

DO_STABILITY <- TRUE
B_STAB       <- 200L
P_CRIT_STAB  <- 0.05

STAB_N_IMP  <- 5L
STAB_IMP_IDS <- NULL
VIF_N_IMP   <- 5L
VIF_IMP_IDS <- NULL

BOOT_SEP_STRATEGY <- "adaptive"

BETA_MAX   <- 15
SE_MAX     <- 30
FITTED_EPS <- 1e-6

set.seed(20240601)

CAL_COL  <- "#2166AC"; CAL_RIBBON <- "#2166AC"; REF_COL <- "grey40"
ROC_COL  <- "#B2182B"; BIN_COL    <- "#444444"

read_safe <- function(f) if (file.exists(f)) readRDS(f) else NULL

reduce_imps <- function(obj, n0, outlier_ids) {
  imps_full <- obj$smcfcs$impDatasets
  bad <- sort(unique(outlier_ids)); bad <- bad[bad %in% seq_len(n0)]
  if (length(bad) > 0) lapply(imps_full, function(d) d[-bad, , drop = FALSE]) else imps_full
}

fit_get_coef <- function(smf, d, path) {
  if (path == "firth") {
    if (!HAS_LOGISTF) stop("Firth 経路だが logistf が無い。install.packages('logistf')")
    fit <- tryCatch(logistf::logistf(smf, data = d, pl = FALSE), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    b <- coef(fit); if (is.null(names(b))) names(b) <- rownames(vcov(fit))
    b
  } else {
    fit <- tryCatch(suppressWarnings(glm(smf, family = binomial(), data = d)),
                    error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    coef(fit)
  }
}

glm_unstable <- function(fit) {
  if (is.null(fit)) return(TRUE)
  if (!isTRUE(fit$converged)) return(TRUE)
  b <- coef(fit); if (any(!is.finite(b)) || any(abs(b) > BETA_MAX)) return(TRUE)
  se <- tryCatch(sqrt(diag(vcov(fit))), error = function(e) NA_real_)
  if (any(!is.finite(se)) || any(se > SE_MAX)) return(TRUE)
  fv <- fit$fitted.values
  if (any(fv < FITTED_EPS | fv > 1 - FITTED_EPS)) return(TRUE)
  FALSE
}

fit_robust <- function(smf, d, base_path, strategy = BOOT_SEP_STRATEGY) {
  use_firth_always <- identical(base_path, "firth") || identical(strategy, "firth")
  if (use_firth_always) {
    b <- fit_get_coef(smf, d, "firth")
    if (is.null(b)) return(NULL)
    return(list(coef = b, estimator = "firth"))
  }
  fit <- tryCatch(suppressWarnings(glm(smf, family = binomial(), data = d)),
                  error = function(e) NULL)
  if (!is.null(fit) && !glm_unstable(fit))
    return(list(coef = coef(fit), estimator = "glm"))
  b <- fit_get_coef(smf, d, "firth")
  if (is.null(b)) {
    if (!is.null(fit)) return(list(coef = coef(fit), estimator = "glm(unstable)"))
    return(NULL)
  }
  list(coef = b, estimator = "firth")
}

firth_pll <- function(fit) {
  ll <- fit$loglik
  if (!is.null(names(ll)) && "full" %in% names(ll)) return(unname(ll[["full"]]))
  ll[length(ll)]
}

pick_imps <- function(m, k, ids = NULL) {
  if (!is.null(ids)) return(sort(unique(ids[ids >= 1 & ids <= m])))
  if (is.null(k) || k >= m) return(seq_len(m))
  unique(round(seq(1, m, length.out = k)))
}

predict_lp <- function(coef_named, smf, newdata) {
  mm <- tryCatch(stats::model.matrix(stats::delete.response(stats::terms(smf)),
                                      data = newdata),
                 error = function(e) NULL)
  if (is.null(mm)) return(rep(NA_real_, nrow(newdata)))
  common <- intersect(colnames(mm), names(coef_named))
  if (length(common) == 0) return(rep(NA_real_, nrow(newdata)))
  as.numeric(mm[, common, drop = FALSE] %*% coef_named[common])
}

prob_from_lp <- function(lp) {
  p <- stats::plogis(lp)
  pmin(pmax(p, PROB_EPS), 1 - PROB_EPS)
}

auc_one <- function(y, p) {
  pos <- p[y == 1]; neg <- p[y == 0]
  if (length(pos) == 0 || length(neg) == 0)
    return(list(auc = NA_real_, se = NA_real_))
  if (HAS_PROC) {
    r <- tryCatch(pROC::roc(response = y, predictor = p, quiet = TRUE,
                            levels = c(0, 1), direction = "<"),
                  error = function(e) NULL)
    if (!is.null(r)) {
      a  <- as.numeric(pROC::auc(r))
      ci <- tryCatch(as.numeric(pROC::ci.auc(r, method = "delong",
                                             conf.level = CI_LEVEL)),
                     error = function(e) NULL)
      se <- if (!is.null(ci) && length(ci) == 3)
              (ci[3] - ci[1]) / (2 * stats::qnorm(1 - (1 - CI_LEVEL) / 2)) else NA_real_
      return(list(auc = a, se = se))
    }
  }
  rr <- rank(c(pos, neg))
  a  <- (sum(rr[seq_along(pos)]) - length(pos) * (length(pos) + 1) / 2) /
        (length(pos) * length(neg))
  list(auc = a, se = NA_real_)
}

brier_one <- function(y, p) mean((p - y)^2)

nagelkerke_one <- function(y, p) {
  n <- length(y); pbar <- mean(y)
  if (pbar <= 0 || pbar >= 1) return(NA_real_)
  ll_m <- sum(y * log(p) + (1 - y) * log(1 - p))
  ll_0 <- sum(y * log(pbar) + (1 - y) * log(1 - pbar))
  cox  <- 1 - exp((2 / n) * (ll_0 - ll_m))
  cox / (1 - exp((2 / n) * ll_0))
}

cal_one <- function(y, lp) {
  fs <- tryCatch(suppressWarnings(glm(y ~ lp, family = binomial())),
                 error = function(e) NULL)
  fi <- tryCatch(suppressWarnings(glm(y ~ 1, family = binomial(), offset = lp)),
                 error = function(e) NULL)
  list(
    slope    = if (!is.null(fs)) unname(coef(fs)["lp"]) else NA_real_,
    se_slope = if (!is.null(fs)) sqrt(vcov(fs)["lp", "lp"]) else NA_real_,
    intercept= if (!is.null(fi)) unname(coef(fi)[1]) else NA_real_,
    se_int   = if (!is.null(fi)) sqrt(vcov(fi)[1, 1]) else NA_real_
  )
}

hl_one <- function(y, p, g = HL_GROUPS) {
  if (HAS_HL) {
    r <- tryCatch(ResourceSelection::hoslem.test(y, p, g = g), error = function(e) NULL)
    if (!is.null(r)) return(list(stat = unname(r$statistic), df = unname(r$parameter),
                                 p = unname(r$p.value)))
  }
  br <- unique(stats::quantile(p, probs = seq(0, 1, length.out = g + 1), na.rm = TRUE))
  if (length(br) < 3) return(list(stat = NA_real_, df = NA_real_, p = NA_real_))
  grp <- cut(p, breaks = br, include.lowest = TRUE)
  obs <- tapply(y, grp, sum); exp <- tapply(p, grp, sum)
  n_g <- tapply(y, grp, length)
  ok  <- is.finite(exp) & exp > 0 & exp < n_g
  chi <- sum(((obs - exp)^2 / (exp * (1 - exp / n_g)))[ok], na.rm = TRUE)
  df  <- sum(ok) - 2
  list(stat = chi, df = df, p = if (df > 0) stats::pchisq(chi, df, lower.tail = FALSE) else NA_real_)
}

rubin_scalar <- function(est, se, ci = CI_LEVEL) {
  keep <- is.finite(est); est <- est[keep]; se <- se[keep]
  m <- length(est)
  if (m == 0) return(list(est = NA_real_, se = NA_real_, lcl = NA_real_, ucl = NA_real_))
  Qbar <- mean(est)
  Ubar <- if (all(is.finite(se))) mean(se^2) else NA_real_
  B    <- if (m >= 2) stats::var(est) else 0
  Tvar <- if (is.finite(Ubar)) Ubar + (1 + 1 / m) * B else (1 + 1 / m) * B
  SE   <- sqrt(max(Tvar, 0))
  z    <- stats::qnorm(1 - (1 - ci) / 2)
  list(est = Qbar, se = SE, lcl = Qbar - z * SE, ucl = Qbar + z * SE)
}

pool_auc <- function(aucs, ses, ci = CI_LEVEL) {
  a <- aucs; ok <- is.finite(a) & a > 0 & a < 1
  a <- a[ok]; s <- ses[ok]
  if (length(a) == 0) return(list(auc = NA_real_, lcl = NA_real_, ucl = NA_real_))
  la <- stats::qlogis(a)
  ls <- ifelse(is.finite(s), s / (a * (1 - a)), NA_real_)
  pr <- rubin_scalar(la, ls, ci)
  list(auc = stats::plogis(pr$est),
       lcl = stats::plogis(pr$lcl), ucl = stats::plogis(pr$ucl))
}

pool_hl_D2 <- function(stats_vec, df_vec) {
  s <- stats_vec[is.finite(stats_vec) & stats_vec >= 0]
  k <- as.integer(round(stats::median(df_vec[is.finite(df_vec)], na.rm = TRUE)))
  m <- length(s)
  if (m == 0 || is.na(k) || k <= 0) return(list(p = NA_real_, D = NA_real_, df1 = k))
  mbar <- mean(s)
  if (m == 1) return(list(p = stats::pchisq(mbar, df = k, lower.tail = FALSE),
                          D = mbar / k, df1 = k))
  r  <- (1 + 1 / m) * stats::var(sqrt(s)); if (!is.finite(r) || r < 0) r <- 0
  D  <- (mbar / k - (m + 1) / (m - 1) * r) / (1 + r); if (!is.finite(D) || D < 0) D <- 0
  df2 <- if (r <= 0) 1e6 else k^(-3 / m) * (m - 1) * (1 + 1 / r)^2
  list(p = stats::pf(D, df1 = k, df2 = df2, lower.tail = FALSE), D = D, df1 = k)
}

eval_model <- function(imps, smf, path) {
  m <- length(imps); y <- as.integer(imps[[1]]$good)
  n <- length(y)
  P  <- matrix(NA_real_, nrow = n, ncol = m)
  coefs <- vector("list", m)
  auc_v <- se_v <- bri_v <- nag_v <- numeric(m)
  cs_v  <- cs_se <- ci_v <- ci_se <- numeric(m)
  hl_s  <- hl_df <- hl_p <- numeric(m)

  for (i in seq_len(m)) {
    b <- fit_get_coef(smf, imps[[i]], path)
    coefs[[i]] <- b
    if (is.null(b)) { auc_v[i] <- se_v[i] <- bri_v[i] <- nag_v[i] <- NA;
                      cs_v[i] <- cs_se[i] <- ci_v[i] <- ci_se[i] <- NA;
                      hl_s[i] <- hl_df[i] <- hl_p[i] <- NA; next }
    lp <- predict_lp(b, smf, imps[[i]]); p <- prob_from_lp(lp); P[, i] <- p
    a  <- auc_one(y, p);  auc_v[i] <- a$auc; se_v[i] <- a$se
    bri_v[i] <- brier_one(y, p); nag_v[i] <- nagelkerke_one(y, p)
    cc <- cal_one(y, lp); cs_v[i] <- cc$slope; cs_se[i] <- cc$se_slope
    ci_v[i] <- cc$intercept; ci_se[i] <- cc$se_int
    hh <- hl_one(y, p); hl_s[i] <- hh$stat; hl_df[i] <- hh$df; hl_p[i] <- hh$p
  }

  auc_pool   <- pool_auc(auc_v, se_v)
  slope_pool <- rubin_scalar(cs_v, cs_se)
  int_pool   <- rubin_scalar(ci_v, ci_se)
  hl_pool    <- pool_hl_D2(hl_s, hl_df)

  perf <- data.table(
    n = n, events = sum(y == 1), nonevents = sum(y == 0), m = m,
    AUC = auc_pool$auc, AUC_lcl = auc_pool$lcl, AUC_ucl = auc_pool$ucl,
    Brier = mean(bri_v, na.rm = TRUE),
    Nagelkerke_R2 = mean(nag_v, na.rm = TRUE),
    cal_slope = slope_pool$est, cal_slope_lcl = slope_pool$lcl, cal_slope_ucl = slope_pool$ucl,
    cal_intercept = int_pool$est, cal_intercept_lcl = int_pool$lcl, cal_intercept_ucl = int_pool$ucl,
    HL_stat = hl_pool$D * hl_pool$df1, HL_df = hl_pool$df1, HL_p = hl_pool$p,
    HL_p_median = stats::median(hl_p, na.rm = TRUE),
    HL_reject_frac = mean(hl_p < 0.05, na.rm = TRUE)
  )
  list(perf = perf, p_pool = rowMeans(P, na.rm = TRUE), y = y, coefs = coefs)
}

calib_curve_data <- function(y, p, bins = CAL_BINS, n_grid = N_GRID) {
  d <- data.table(y = y, p = pmin(pmax(p, 0), 1))
  lo <- tryCatch(stats::loess(y ~ p, data = d, span = 0.9, degree = 1,
                              family = "gaussian", control = stats::loess.control(surface = "direct")),
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
      zc <- stats::qnorm(1 - (1 - CI_LEVEL) / 2)
      den <- 1 + zc^2 / n; ctr <- (oh + zc^2 / (2 * n)) / den
      hw  <- zc * sqrt(oh * (1 - oh) / n + zc^2 / (4 * n^2)) / den
      .(mean_pred = ph, obs_prop = oh, lcl = max(ctr - hw, 0), ucl = min(ctr + hw, 1), n = n)
    }, by = grp]
  }
  list(curve = curve, bins = binpts)
}

roc_curve_data <- function(y, p) {
  if (HAS_PROC) {
    r <- tryCatch(pROC::roc(response = y, predictor = p, quiet = TRUE,
                            levels = c(0, 1), direction = "<"),
                  error = function(e) NULL)
    if (!is.null(r))
      return(data.table(fpr = 1 - r$specificities, tpr = r$sensitivities)[order(fpr, tpr)])
  }
  thr <- sort(unique(c(-Inf, p, Inf)), decreasing = TRUE)
  np <- sum(y == 1); nn <- sum(y == 0)
  tpr <- vapply(thr, function(t) sum(p >  t & y == 1) / np, numeric(1))
  fpr <- vapply(thr, function(t) sum(p >  t & y == 0) / nn, numeric(1))
  data.table(fpr = fpr, tpr = tpr)[order(fpr, tpr)]
}

vif_one <- function(smf, d) {
  if (!HAS_CAR) return(NULL)
  tl <- attr(stats::terms(smf), "term.labels")
  if (length(tl) < 2) return(NULL)
  fit <- tryCatch(suppressWarnings(glm(smf, family = binomial(), data = d)),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  v <- tryCatch(car::vif(fit), error = function(e) NULL)
  if (is.null(v)) return(NULL)
  if (is.matrix(v)) {
    dt <- as.data.table(v, keep.rownames = "term")
    setnames(dt, old = names(dt)[-1],
             new = sub("GVIF\\^\\(1/\\(2\\*Df\\)\\)", "GVIF_adj", names(dt)[-1]))
    if ("GVIF_adj" %in% names(dt)) dt[, GVIF2 := GVIF_adj^2]
    dt[]
  } else {
    data.table(term = names(v), GVIF = as.numeric(v), Df = 1, GVIF_adj = sqrt(as.numeric(v)),
               GVIF2 = as.numeric(v))
  }
}

bootstrap_optimism <- function(imps, smf, path, n_boot = N_BOOT, m_use = M_BOOT,
                               strategy = BOOT_SEP_STRATEGY) {
  m <- length(imps)
  use_idx <- pick_imps(m, m_use)
  opt_auc <- opt_slope <- app_auc <- app_slope <- numeric(0)
  n_firth_fits <- 0L; n_total_fits <- 0L

  fit_rb <- function(d) {
    r <- fit_robust(smf, d, base_path = path, strategy = strategy)
    if (!is.null(r)) {
      n_total_fits <<- n_total_fits + 1L
      if (grepl("firth", r$estimator)) n_firth_fits <<- n_firth_fits + 1L
    }
    r
  }

  for (i in use_idx) {
    d <- imps[[i]]; y <- as.integer(d$good); n <- nrow(d)
    r0 <- fit_rb(d); if (is.null(r0)) next
    p0  <- prob_from_lp(predict_lp(r0$coef, smf, d))
    a_app <- auc_one(y, p0)$auc; s_app <- cal_one(y, qlogis(p0))$slope
    if (!is.finite(a_app)) next
    o_a <- o_s <- numeric(0)
    for (bnum in seq_len(n_boot)) {
      idx <- sample.int(n, n, replace = TRUE)
      db  <- d[idx, , drop = FALSE]
      rb  <- fit_rb(db); if (is.null(rb)) next
      bb  <- rb$coef
      pbb <- prob_from_lp(predict_lp(bb, smf, db))
      a_b <- auc_one(as.integer(db$good), pbb)$auc
      s_b <- cal_one(as.integer(db$good), qlogis(pbb))$slope
      pbo <- prob_from_lp(predict_lp(bb, smf, d))
      a_t <- auc_one(y, pbo)$auc
      s_t <- cal_one(y, qlogis(pbo))$slope
      if (is.finite(a_b) && is.finite(a_t)) o_a <- c(o_a, a_b - a_t)
      if (is.finite(s_b) && is.finite(s_t)) o_s <- c(o_s, s_b - s_t)
    }
    app_auc   <- c(app_auc, a_app);  app_slope <- c(app_slope, s_app)
    opt_auc   <- c(opt_auc, mean(o_a, na.rm = TRUE))
    opt_slope <- c(opt_slope, mean(o_s, na.rm = TRUE))
  }
  firth_frac <- if (n_total_fits > 0) n_firth_fits / n_total_fits else NA_real_
  if (length(app_auc) == 0)
    return(data.table(AUC_apparent = NA_real_, AUC_optimism = NA_real_, AUC_corrected = NA_real_,
                      slope_apparent = NA_real_, slope_optimism = NA_real_, slope_corrected = NA_real_,
                      n_boot = n_boot, m_used = length(use_idx),
                      sep_strategy = strategy, firth_fit_frac = firth_frac))
  data.table(
    AUC_apparent  = mean(app_auc),  AUC_optimism  = mean(opt_auc, na.rm = TRUE),
    AUC_corrected = mean(app_auc) - mean(opt_auc, na.rm = TRUE),
    slope_apparent  = mean(app_slope), slope_optimism = mean(opt_slope, na.rm = TRUE),
    slope_corrected = mean(app_slope) - mean(opt_slope, na.rm = TRUE),
    n_boot = n_boot, m_used = length(use_idx),
    sep_strategy = strategy, firth_fit_frac = firth_frac)
}

group_key <- function(lab) {
  if (grepl("_measurable$", lab)) return(lab)
  if (grepl("_c[0-9]+$", lab))    return(sub("_c[0-9]+$", "", lab))
  if (grepl("_c$", lab))          return(sub("_c$", "", lab))
  lab
}
backward_lrt <- function(full_labels, d, p_crit = P_CRIT_STAB, keep = character(0),
                         base_path = "glm", strategy = BOOT_SEP_STRATEGY) {
  keys   <- vapply(full_labels, group_key, character(1))
  groups <- split(full_labels, keys)
  current <- names(groups)
  mk <- function(labs) if (length(labs) == 0) stats::as.formula("good ~ 1")
                       else stats::reformulate(labs, response = "good")
  fit_glm   <- function(labs) tryCatch(suppressWarnings(glm(mk(labs), binomial(), d)),
                                       error = function(e) NULL)
  fit_firth <- function(labs) tryCatch(logistf::logistf(mk(labs), data = d, pl = FALSE),
                                       error = function(e) NULL)
  always_firth <- identical(base_path, "firth") || identical(strategy, "firth")

  repeat {
    cand <- setdiff(current, keep)
    if (length(cand) == 0) break
    cur_labs <- unlist(groups[current], use.names = FALSE)

    mode_firth <- always_firth; gf <- NULL
    if (!always_firth) {
      gf <- fit_glm(cur_labs)
      if (is.null(gf) || glm_unstable(gf)) mode_firth <- TRUE
    }

    if (!mode_firth) {
      ll_full <- as.numeric(stats::logLik(gf)); k_full <- length(coef(gf))
      pv <- vapply(cand, function(g) {
        red <- unlist(groups[setdiff(current, g)], use.names = FALSE)
        gr <- fit_glm(red); if (is.null(gr)) return(NA_real_)
        lr <- 2 * (ll_full - as.numeric(stats::logLik(gr))); k <- k_full - length(coef(gr))
        if (k <= 0) return(NA_real_); stats::pchisq(lr, df = k, lower.tail = FALSE)
      }, numeric(1))
    } else {
      ff <- fit_firth(cur_labs); if (is.null(ff)) break
      ll_full <- firth_pll(ff); k_full <- length(coef(ff))
      pv <- vapply(cand, function(g) {
        red <- unlist(groups[setdiff(current, g)], use.names = FALSE)
        fr <- fit_firth(red); if (is.null(fr)) return(NA_real_)
        lr <- 2 * (ll_full - firth_pll(fr)); k <- k_full - length(coef(fr))
        if (k <= 0) return(NA_real_); stats::pchisq(lr, df = k, lower.tail = FALSE)
      }, numeric(1))
    }
    worst <- which.max(pv)
    if (length(worst) == 0 || !is.finite(pv[worst]) || pv[worst] < p_crit) break
    current <- setdiff(current, cand[worst])
  }
  current
}
model_stability <- function(full_labels, imps, imp_ids, b_stab = B_STAB,
                            p_crit = P_CRIT_STAB, keep = character(0),
                            base_path = "glm", strategy = BOOT_SEP_STRATEGY) {
  keys     <- vapply(full_labels, group_key, character(1))
  all_keys <- unique(keys)
  sel_count <- setNames(numeric(length(all_keys)), all_keys)
  total_runs <- 0L
  for (ii in imp_ids) {
    d <- imps[[ii]]; n <- nrow(d)
    for (b in seq_len(b_stab)) {
      idx <- sample.int(n, n, replace = TRUE)
      sel <- tryCatch(backward_lrt(full_labels, d[idx, , drop = FALSE], p_crit, keep,
                                   base_path, strategy),
                      error = function(e) character(0))
      sel_count[intersect(sel, all_keys)] <- sel_count[intersect(sel, all_keys)] + 1
      total_runs <- total_runs + 1L
    }
  }
  data.table(var = all_keys, variable = relabel_vars(all_keys),
             selected_freq = as.integer(round(sel_count[all_keys])),
             selected_pct = 100 * sel_count[all_keys] / max(total_runs, 1L))
}

process_pipeline <- function(nm) {
  prov <- read_safe(file.path(OUT_DIR, sprintf("bnb_provisional_%s.rds", nm)))
  fin  <- read_safe(file.path(OUT_DIR, sprintf("bnb_final_%s.rds",       nm)))
  obj  <- read_safe(file.path(OUT_DIR, sprintf("bnb_imp_%s.rds",         nm)))
  if (is.null(obj)) { message(sprintf("[%s] 代入データが無い。スキップ。", nm)); return(NULL) }
  if (is.null(prov) && is.null(fin)) {
    message(sprintf("[%s] 暫定/最終モデルが無い。スキップ。", nm)); return(NULL)
  }
  message(sprintf("\n================ [%s] §12 性能評価 ================", nm))

  models <- list()
  if (!is.null(prov)) models$Provisional <-
    list(smf = stats::as.formula(prov$smformula_prov),  path = prov$path,
         n0 = prov$n0, outlier_ids = prov$outlier_ids,
         full_labels = attr(stats::terms(stats::as.formula(prov$smformula_full)), "term.labels"),
         keep = prov$keep_vars)
  if (!is.null(fin)) models$Final <-
    list(smf = stats::as.formula(fin$smformula_final), path = fin$path,
         n0 = fin$n0, outlier_ids = fin$outlier_ids,
         full_labels = NULL, keep = NULL)

  perf_rows <- list(); calib_list <- list(); roc_list <- list()
  boot_list <- list(); vif_list <- list(); stab_obj <- NULL
  pp_store  <- list()

  for (md in names(models)) {
    mo  <- models[[md]]
    imps <- reduce_imps(obj, mo$n0, mo$outlier_ids)
    m <- length(imps)
    tl <- attr(stats::terms(mo$smf), "term.labels")
    message(sprintf("  [%s] path=%s, 項数=%d, m=%d", md, mo$path, length(tl), m))

    if (length(tl) == 0) {
      message("    切片のみモデル → 性能評価をスキップ。"); next
    }

    if (have_tic) tictoc::tic(sprintf("  [%s] eval_model", md))
    ev <- eval_model(imps, mo$smf, mo$path)
    if (have_tic) tictoc::toc()
    perf <- copy(ev$perf)[, `:=`(pipeline = nm, pipeline_disp = relabel_pipeline(nm),
                                 model = md, path = mo$path)]
    perf_rows[[md]] <- perf
    pp_store[[md]]  <- data.table(pipeline = nm, model = md, y = ev$y, p = ev$p_pool)

    cc <- calib_curve_data(ev$y, ev$p_pool)
    if (nrow(cc$curve) > 0) calib_list[[paste0(md, "_curve")]] <-
      copy(cc$curve)[, `:=`(pipeline = nm, model = md, kind = "curve")]
    if (nrow(cc$bins)  > 0) calib_list[[paste0(md, "_bins")]]  <-
      copy(cc$bins)[, `:=`(pipeline = nm, model = md)]
    rc <- roc_curve_data(ev$y, ev$p_pool)
    roc_list[[md]] <- copy(rc)[, `:=`(pipeline = nm, model = md,
                                      auc = perf$AUC)]

    vif_ids <- pick_imps(m, VIF_N_IMP, VIF_IMP_IDS)
    vlist <- lapply(vif_ids, function(ii) {
      v <- vif_one(mo$smf, imps[[ii]]); if (!is.null(v)) v[, .imp := ii]; v
    })
    vlist <- vlist[!vapply(vlist, is.null, logical(1))]
    if (length(vlist) > 0) {
      vraw <- rbindlist(vlist, fill = TRUE)
      vf <- vraw[, .(GVIF = mean(GVIF, na.rm = TRUE),
                     Df = stats::median(Df, na.rm = TRUE),
                     GVIF_adj = mean(GVIF_adj, na.rm = TRUE),
                     n_imp = .N), by = term]
      vf[, `:=`(pipeline = nm, model = md)]; vif_list[[md]] <- vf
      perf_rows[[md]][, `:=`(max_GVIF = max(vf$GVIF, na.rm = TRUE),
                             max_GVIF_adj = max(vf$GVIF_adj, na.rm = TRUE),
                             vif_n_imp = length(vif_ids))]
    } else perf_rows[[md]][, `:=`(max_GVIF = NA_real_, max_GVIF_adj = NA_real_,
                                  vif_n_imp = 0L)]

    if (DO_BOOTSTRAP) {
      if (have_tic) tictoc::tic(sprintf("  [%s] bootstrap optimism (B=%d)", md, N_BOOT))
      bo <- bootstrap_optimism(imps, mo$smf, mo$path, N_BOOT, M_BOOT)
      if (have_tic) tictoc::toc()
      bo[, `:=`(pipeline = nm, model = md, path = mo$path)]; boot_list[[md]] <- bo
      perf_rows[[md]][, `:=`(AUC_corrected = bo$AUC_corrected,
                             cal_slope_corrected = bo$slope_corrected)]
    } else perf_rows[[md]][, `:=`(AUC_corrected = NA_real_, cal_slope_corrected = NA_real_)]
  }

  if (DO_STABILITY && !is.null(prov) && length(models$Provisional$full_labels) > 0) {
    imps     <- reduce_imps(obj, prov$n0, prov$outlier_ids)
    m        <- length(imps)
    stab_ids <- pick_imps(m, STAB_N_IMP, STAB_IMP_IDS)
    keep     <- if (!is.null(prov$keep_vars)) prov$keep_vars else character(0)
    base_p   <- prov$path
    strat    <- if (identical(base_p, "firth")) "firth" else BOOT_SEP_STRATEGY
    if (have_tic) tictoc::tic(sprintf("  model stability (B=%d x %d imp, %s)",
                                      B_STAB, length(stab_ids), strat))
    stab <- model_stability(models$Provisional$full_labels, imps, stab_ids,
                            B_STAB, P_CRIT_STAB, keep, base_p, strat)
    if (have_tic) tictoc::toc()
    stab[, in_provisional := var %in% prov$prov_vars]
    setorder(stab, -selected_pct)
    stab_obj <- list(pipeline = nm, b_stab = B_STAB, p_crit = P_CRIT_STAB,
                     imp_ids = stab_ids, n_imp = length(stab_ids),
                     base_path = base_p, sep_strategy = strat,
                     total_runs = B_STAB * length(stab_ids), freq = stab,
                     note = paste0("後退選択を ", length(stab_ids), " 本の代入 × B=", B_STAB,
                                   " のブートストラップで繰り返し、選択頻度を平均。分離対策=", strat,
                                   "（adaptive: 分離相当のみ Firth / firth: 全ステップ Firth）。",
                                   " MIプール(D3/D2)選択そのものの完全再現ではなく代入ごとの単一データ選択の近似。"))
    saveRDS(stab_obj, file.path(OUT_DIR, sprintf("bnb_model_stability_%s.rds", nm)))
    message(sprintf("  モデル安定性（変数選択頻度, 上位 / 代入%d本×B=%d, %s）:",
                    length(stab_ids), B_STAB, strat))
    print(stab[, .(variable, selected_pct = round(selected_pct, 1), in_provisional)][1:min(.N, 8)],
          row.names = FALSE)
  }

  perf_dt <- rbindlist(perf_rows, fill = TRUE)
  saveRDS(list(
    pipeline = nm,
    performance = perf_dt,
    pooled_pred = rbindlist(pp_store, fill = TRUE),
    calibration = if (length(calib_list)) rbindlist(calib_list, fill = TRUE) else data.table(),
    roc         = if (length(roc_list))   rbindlist(roc_list,   fill = TRUE) else data.table(),
    bootstrap   = if (length(boot_list))  rbindlist(boot_list,  fill = TRUE) else data.table(),
    vif         = if (length(vif_list))   rbindlist(vif_list,   fill = TRUE) else data.table(),
    stability   = stab_obj,
    settings    = list(ci_level = CI_LEVEL, hl_groups = HL_GROUPS, cal_bins = CAL_BINS,
                       n_boot = N_BOOT, m_boot = M_BOOT, b_stab = B_STAB)
  ), file.path(OUT_DIR, sprintf("bnb_performance_%s.rds", nm)))

  show <- copy(perf_dt)
  message(sprintf("  [%s] プール性能サマリ:", nm))
  print(show[, .(model, path,
                 AUC = sprintf("%.3f [%.3f, %.3f]", AUC, AUC_lcl, AUC_ucl),
                 AUC_corr = round(AUC_corrected, 3),
                 Brier = round(Brier, 3), R2 = round(Nagelkerke_R2, 3),
                 cal_slope = sprintf("%.2f [%.2f, %.2f]", cal_slope, cal_slope_lcl, cal_slope_ucl),
                 slope_corr = round(cal_slope_corrected, 2),
                 HL_p = signif(HL_p, 3), maxGVIF = round(max_GVIF, 1))],
        row.names = FALSE)

  list(perf = perf_dt,
       calib = if (length(calib_list)) rbindlist(calib_list, fill = TRUE) else data.table(),
       roc   = if (length(roc_list))   rbindlist(roc_list,   fill = TRUE) else data.table())
}

results <- lapply(PIPELINES, process_pipeline)
names(results) <- PIPELINES
results <- results[!vapply(results, is.null, logical(1))]
if (length(results) == 0)
  stop("性能評価できるパイプラインがありません。§8/§10 を先に実行してください。")

perf_all <- rbindlist(lapply(results, `[[`, "perf"), fill = TRUE)
setcolorder(perf_all, c("pipeline", "pipeline_disp", "model", "path"))
perf_all[, pipeline := factor(pipeline, levels = PIPELINES)]
perf_all[, model    := factor(model, levels = c("Provisional", "Final"))]
setorder(perf_all, pipeline, model)
saveRDS(perf_all, file.path(OUT_DIR, "bnb_performance_summary.rds"))
fwrite(perf_all,  file.path(OUT_DIR, "bnb_performance_summary.csv"))

message("\n==== §12 性能サマリ（全パイプライン×モデル）====")
print(perf_all[, .(Pipeline = pipeline_disp, Model = as.character(model), Path = path,
                   AUC = sprintf("%.3f", AUC), AUC_corr = sprintf("%.3f", AUC_corrected),
                   Brier = sprintf("%.3f", Brier), R2 = sprintf("%.3f", Nagelkerke_R2),
                   cal_slope = sprintf("%.2f", cal_slope), slope_corr = sprintf("%.2f", cal_slope_corrected),
                   cal_int = sprintf("%.2f", cal_intercept), HL_p = signif(HL_p, 3),
                   maxGVIF = round(max_GVIF, 1))], row.names = FALSE)

calib_all <- rbindlist(lapply(results, `[[`, "calib"), fill = TRUE)
roc_all   <- rbindlist(lapply(results, `[[`, "roc"),   fill = TRUE)
if (nrow(calib_all) && !"kind" %in% names(calib_all)) calib_all[, kind := NA_character_]
if (nrow(calib_all)) calib_all[, model := factor(model, levels = c("Provisional", "Final"))]
if (nrow(roc_all))   roc_all[,   model := factor(model, levels = c("Provisional", "Final"))]

make_calib <- function(curve_dt, bin_dt, title, subtitle) {
  p <- ggplot() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = REF_COL) +
    { if (nrow(curve_dt))
        geom_ribbon(data = curve_dt, aes(x = p, ymin = lcl, ymax = ucl, fill = model),
                    alpha = 0.15, show.legend = FALSE) else NULL } +
    { if (nrow(curve_dt))
        geom_line(data = curve_dt, aes(x = p, y = obs, colour = model), linewidth = 0.8) else NULL } +
    { if (nrow(bin_dt))
        geom_point(data = bin_dt, aes(x = mean_pred, y = obs_prop, colour = model), size = 1.8) else NULL } +
    { if (nrow(bin_dt))
        geom_errorbar(data = bin_dt, aes(x = mean_pred, ymin = lcl, ymax = ucl, colour = model),
                      width = 0.012, linewidth = 0.4) else NULL } +
    scale_colour_manual(values = c(Provisional = CAL_COL, Final = ROC_COL), name = NULL) +
    scale_fill_manual(values  = c(Provisional = CAL_COL, Final = ROC_COL), name = NULL) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(title = title, subtitle = subtitle,
         x = "Predicted probability", y = "Observed proportion") +
    theme_bw(base_size = 10) +
    theme(legend.position = "top",
          plot.title = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 8.5, colour = "grey30"),
          panel.grid.minor = element_blank())
  p
}

if (nrow(calib_all)) {
  cal_plots <- list()
  for (pp in PIPELINES) {
    cu <- calib_all[pipeline == pp & kind == "curve"]
    bn <- calib_all[pipeline == pp & is.na(kind)]
    if (nrow(cu) == 0 && nrow(bn) == 0) next
    g <- make_calib(cu, bn, relabel_pipeline(pp),
                    "Pooled calibration (MI-averaged predicted prob; loess + decile bins, 95% CI)")
    cal_plots[[pp]] <- g
    ggsave(file.path(FIG_DIR, sprintf("fig5_calibration_%s.png", pp)),
           g, width = 6.0, height = 5.6, dpi = 150)
  }
  if (length(cal_plots)) {
    panels <- lapply(names(cal_plots), function(pp)
      cal_plots[[pp]] + labs(subtitle = NULL))
    combo <- patchwork::wrap_plots(panels, ncol = 2, guides = "collect") +
      patchwork::plot_annotation(
        title = "Figure 5. Calibration plots (pooled across imputations)",
        subtitle = "MI-averaged predicted probability vs observed proportion. Dashed = ideal; points = decile bins (95% Wilson CI).",
        theme = theme(plot.title = element_text(size = 12, face = "bold"),
                      plot.subtitle = element_text(size = 9, colour = "grey30"))) &
      theme(legend.position = "top")
    nr <- ceiling(length(panels) / 2)
    ggsave(file.path(FIG_DIR, "fig5_calibration_all.png"),
           combo, width = 11, height = 5.0 * nr + 0.6, dpi = 150)
  }
  message("\nCalibration plot を保存: ", normalizePath(FIG_DIR),
          "/fig5_calibration_<pipeline>.png（および ..._all.png）")
}

make_roc <- function(d, title) {
  lab <- d[, .(auc = auc[1]), by = model]
  lab[, lab := sprintf("%s (AUC %.3f)", model, auc)]
  labmap <- stats::setNames(lab$lab, as.character(lab$model))
  ggplot(d, aes(x = fpr, y = tpr, colour = model)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = REF_COL) +
    geom_line(linewidth = 0.8) +
    scale_colour_manual(values = c(Provisional = CAL_COL, Final = ROC_COL),
                        labels = labmap, name = NULL) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(title = title,
         subtitle = "Pooled ROC (MI-averaged predicted probability); AUC = Rubin-pooled across imputations",
         x = "1 - Specificity (FPR)", y = "Sensitivity (TPR)") +
    theme_bw(base_size = 10) +
    theme(legend.position = "top",
          plot.title = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 8, colour = "grey30"),
          panel.grid.minor = element_blank())
}

if (nrow(roc_all)) {
  roc_all <- merge(roc_all,
                   perf_all[, .(pipeline = as.character(pipeline), model, auc_pool = AUC)],
                   by = c("pipeline", "model"), all.x = TRUE)
  roc_all[is.finite(auc_pool), auc := auc_pool]
  roc_plots <- list()
  for (pp in PIPELINES) {
    d <- roc_all[pipeline == pp]
    if (nrow(d) == 0) next
    g <- make_roc(copy(d), relabel_pipeline(pp))
    roc_plots[[pp]] <- g
    ggsave(file.path(FIG_DIR, sprintf("fig6_roc_%s.png", pp)),
           g, width = 5.8, height = 5.6, dpi = 150)
  }
  if (length(roc_plots)) {
    panels <- lapply(names(roc_plots), function(pp)
      roc_plots[[pp]] + labs(subtitle = NULL))
    combo <- patchwork::wrap_plots(panels, ncol = 2) +
      patchwork::plot_annotation(
        title = "Figure 6. ROC curves with pooled AUC",
        subtitle = "Pooled ROC from MI-averaged predicted probabilities; AUC = Rubin-pooled across imputations (logit scale).",
        theme = theme(plot.title = element_text(size = 12, face = "bold"),
                      plot.subtitle = element_text(size = 9, colour = "grey30")))
    nr <- ceiling(length(panels) / 2)
    ggsave(file.path(FIG_DIR, "fig6_roc_all.png"),
           combo, width = 11, height = 5.0 * nr + 0.6, dpi = 150)
  }
  message("ROC plot を保存: ", normalizePath(FIG_DIR),
          "/fig6_roc_<pipeline>.png（および ..._all.png）")
}

message("\n§12 完了。暫定・最終モデルのプール性能を評価しました:")
message("  ・図5 calibration plot（fig5_calibration_*.png）、図6 ROC+AUC（fig6_roc_*.png）")
message("  ・プール AUC（logit-Rubin, 95%CI）/ Brier / Nagelkerke R² / 較正切片・勾配 / Hosmer–Lemeshow(D2プール)")
message(sprintf("  ・VIF（GVIF）: 代入 %s 本を平均",
                ifelse(is.null(VIF_N_IMP), "all", as.character(VIF_N_IMP))))
message(sprintf("  ・ブートストラップ内部妥当性: optimism 補正 AUC・較正勾配（B=%d, 代入 m_used≤%s, 分離対策=%s）",
                N_BOOT, ifelse(is.null(M_BOOT), "all", as.character(M_BOOT)), BOOT_SEP_STRATEGY))
message(sprintf("  ・モデル安定性: フルモデルからの後退選択の変数選択頻度（B=%d × 代入 %s 本, 分離対策=%s）",
                B_STAB, ifelse(is.null(STAB_N_IMP), "all", as.character(STAB_N_IMP)), BOOT_SEP_STRATEGY))
message("  保存: data/bnb_performance_<pipeline>.rds, bnb_performance_summary.rds/.csv, ",
        "bnb_model_stability_<pipeline>.rds")
message("【変更1】安定性評価・VIF とも代入 1 本固定をやめ、STAB_N_IMP / VIF_N_IMP（または *_IMP_IDS）で本数を調整可能。")
message("【変更2】外れ値除外後ブートストラップの分離リスクに対し BOOT_SEP_STRATEGY=", BOOT_SEP_STRATEGY,
        " を適用（adaptive: 各当てはめで(準)分離を評価し分離相当のみ Firth / firth: 全当てはめ・全LRTを Firth）。",
        " 元モデルが Firth 経路の場合は本設定に依らず常に Firth。")
message("注記（プラン付録C）: H&L は分位群依存のため calibration plot を主・H&L を補助とする。",
        " 交互作用/スプライン項を含むモデルの GVIF は本質的でない高値を取りうる（中央値centeringで緩和）。",
        " 後退選択後の性能は楽観バイアスを持つため optimism 補正値を併記。",
        " 安定性評価は代入ごとの単一データ後退選択の近似で、§8 の MI プール選択(D3/D2)そのものの完全再現ではない。")
