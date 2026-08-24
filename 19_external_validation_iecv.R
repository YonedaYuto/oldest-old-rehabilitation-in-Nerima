library(data.table)
library(here)

OUT_DIR <- here::here("data")
FIG_DIR <- here::here("figures")
for (d in c(OUT_DIR, FIG_DIR)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

LABELS_PATH <- here::here("00_labels.R")
if (!exists("relabel_vars")) {
  if (!file.exists(LABELS_PATH)) stop("00_labels.R not found: ", LABELS_PATH)
  source(LABELS_PATH)
}

HAS_LOGISTF <- requireNamespace("logistf", quietly = TRUE)
HAS_PROC    <- requireNamespace("pROC",    quietly = TRUE)
HAS_GG      <- requireNamespace("ggplot2", quietly = TRUE)

PIPELINES <- names(PIPE_LABELS)
CI_LEVEL  <- 0.95
read_safe <- function(f) if (file.exists(f)) readRDS(f) else NULL
`%||%`    <- function(a, b) if (is.null(a)) b else a
set.seed(20240601)

CLUSTER_VAR <- "class"

CLUSTER_DERIVE <- NULL

CLUSTER_ORDER <- NULL

CLUSTER_IN_MODEL <- "drop"

MIN_CLUSTER_N        <- 30L
MIN_CLUSTER_EVENTS   <- 5L
MIN_TRAIN_EVENTS     <- 10L

IECV_N_IMP <- 25L

IECV_SEP_STRATEGY <- "adaptive"
BETA_MAX   <- 15; SE_MAX <- 30; FITTED_EPS <- 1e-6
PROB_EPS   <- 1e-8

META_METHOD <- "DL"

COL_POINT <- "#2166AC"; COL_SUMMARY <- "#B2182B"; COL_REF <- "grey40"

reduce_imps <- function(obj, n0, outlier_ids) {
  imps_full <- obj$smcfcs$impDatasets
  bad <- sort(unique(outlier_ids)); bad <- bad[bad %in% seq_len(n0)]
  if (length(bad) > 0) lapply(imps_full, function(d) d[-bad, , drop = FALSE]) else imps_full
}

pick_imps <- function(m, k) {
  if (is.null(k) || k >= m) return(seq_len(m))
  unique(round(seq(1, m, length.out = k)))
}

prob_from_lp <- function(lp) pmin(pmax(stats::plogis(lp), PROB_EPS), 1 - PROB_EPS)

glm_unstable_fit <- function(b, se, fitted) {
  if (any(!is.finite(b)) || any(abs(b) > BETA_MAX)) return(TRUE)
  if (any(!is.finite(se)) || any(se > SE_MAX))      return(TRUE)
  any(fitted < FITTED_EPS | fitted > 1 - FITTED_EPS)
}

fit_train <- function(x, y, base_path, strategy = IECV_SEP_STRATEGY) {
  cn <- colnames(x)
  firth_fit <- function() {
    if (!HAS_LOGISTF) return(NULL)
    keep <- setdiff(cn, "(Intercept)")
    df <- as.data.frame(x[, keep, drop = FALSE])
    safe <- make.names(keep, unique = TRUE)
    names(df) <- safe
    df$..y.. <- y
    ff <- tryCatch(logistf::logistf(stats::as.formula("..y.. ~ ."), data = df, pl = FALSE),
                   error = function(e) NULL)
    if (is.null(ff)) return(NULL)
    b <- coef(ff)
    if (is.null(names(b))) names(b) <- c("(Intercept)", safe)
    out <- stats::setNames(rep(NA_real_, length(cn)), cn)
    out["(Intercept)"] <- unname(b[1])
    idx <- match(safe, names(b))
    out[keep] <- unname(b[idx])
    if (any(!is.finite(out))) return(NULL)
    list(coef = out, estimator = "firth")
  }

  if (identical(base_path, "firth") || identical(strategy, "firth")) return(firth_fit())

  fit <- tryCatch(suppressWarnings(
    stats::glm.fit(x = x, y = y, family = stats::binomial())), error = function(e) NULL)
  if (!is.null(fit) && isTRUE(fit$converged)) {
    b <- fit$coefficients
    names(b) <- cn
    se <- tryCatch({
      xtwx <- crossprod(x * sqrt(fit$weights))
      sqrt(diag(solve(xtwx)))
    }, error = function(e) rep(NA_real_, length(b)))
    if (!any(is.na(b)) && !glm_unstable_fit(b, se, fit$fitted.values))
      return(list(coef = b, estimator = "glm"))
  }
  ff <- firth_fit()
  if (!is.null(ff)) ff$estimator <- "firth (fallback)"
  ff
}

auc_one <- function(y, p) {
  pos <- p[y == 1]; neg <- p[y == 0]
  if (!length(pos) || !length(neg)) return(list(auc = NA_real_, se = NA_real_))
  if (HAS_PROC) {
    r <- tryCatch(pROC::roc(response = y, predictor = p, quiet = TRUE,
                            levels = c(0, 1), direction = "<"), error = function(e) NULL)
    if (!is.null(r)) {
      a  <- as.numeric(pROC::auc(r))
      ci <- tryCatch(as.numeric(pROC::ci.auc(r, method = "delong", conf.level = CI_LEVEL)),
                     error = function(e) NULL)
      se <- if (!is.null(ci) && length(ci) == 3)
        (ci[3] - ci[1]) / (2 * stats::qnorm(1 - (1 - CI_LEVEL) / 2)) else NA_real_
      return(list(auc = a, se = se))
    }
  }
  rr <- rank(c(pos, neg))
  n1 <- length(pos); n0 <- length(neg)
  a  <- (sum(rr[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  q1 <- a / (2 - a); q2 <- 2 * a^2 / (1 + a)
  se <- sqrt((a * (1 - a) + (n1 - 1) * (q1 - a^2) + (n0 - 1) * (q2 - a^2)) / (n1 * n0))
  list(auc = a, se = se)
}

perf_one <- function(y, lp) {
  p <- prob_from_lp(lp)
  a <- auc_one(y, p)
  fs <- tryCatch(suppressWarnings(glm(y ~ lp, family = stats::binomial())),
                 error = function(e) NULL)
  fi <- tryCatch(suppressWarnings(glm(y ~ 1, family = stats::binomial(), offset = lp)),
                 error = function(e) NULL)
  o <- sum(y); e <- sum(p); n <- length(y)
  list(
    c_stat        = a$auc,        se_c_stat        = a$se,
    cal_slope     = if (!is.null(fs)) unname(coef(fs)["lp"]) else NA_real_,
    se_cal_slope  = if (!is.null(fs)) sqrt(vcov(fs)["lp", "lp"]) else NA_real_,
    cal_intercept = if (!is.null(fi)) unname(coef(fi)[1]) else NA_real_,
    se_cal_int    = if (!is.null(fi)) sqrt(vcov(fi)[1, 1]) else NA_real_,
    log_oe        = if (o > 0 && e > 0) log(o / e) else NA_real_,
    se_log_oe     = if (o > 0) sqrt((1 - o / n) / o) else NA_real_,
    brier         = mean((p - y)^2),
    o_events      = o, e_events = e, n = n)
}

rubin1 <- function(est, se) {
  ok <- is.finite(est) & is.finite(se); est <- est[ok]; se <- se[ok]
  m <- length(est)
  if (m == 0) return(c(NA_real_, NA_real_, 0))
  Q <- mean(est); U <- mean(se^2); B <- if (m >= 2) stats::var(est) else 0
  c(Q, sqrt(U + (1 + 1 / m) * B), m)
}
logit  <- function(x) log(x / (1 - x))
expit  <- function(x) 1 / (1 + exp(-x))

meta_dl <- function(est, se, method = META_METHOD, level = CI_LEVEL) {
  ok <- is.finite(est) & is.finite(se) & se > 0
  est <- est[ok]; se <- se[ok]; k <- length(est)
  if (k == 0) return(list(k = 0L, mu = NA_real_, se = NA_real_, lcl = NA_real_,
                          ucl = NA_real_, tau2 = NA_real_, i2 = NA_real_,
                          pi_lo = NA_real_, pi_hi = NA_real_))
  w  <- 1 / se^2
  mu_fe <- sum(w * est) / sum(w)
  Q  <- sum(w * (est - mu_fe)^2)
  df <- k - 1
  C  <- sum(w) - sum(w^2) / sum(w)
  tau2 <- if (identical(method, "FE") || df <= 0 || C <= 0) 0 else max(0, (Q - df) / C)
  wr <- 1 / (se^2 + tau2)
  mu <- sum(wr * est) / sum(wr)
  se_mu <- sqrt(1 / sum(wr))
  z <- stats::qnorm(1 - (1 - level) / 2)
  i2 <- if (df > 0 && Q > 0) max(0, (Q - df) / Q) else 0
  if (k >= 3) {
    tq <- stats::qt(1 - (1 - level) / 2, df = k - 2)
    pi_lo <- mu - tq * sqrt(tau2 + se_mu^2)
    pi_hi <- mu + tq * sqrt(tau2 + se_mu^2)
  } else { pi_lo <- NA_real_; pi_hi <- NA_real_ }
  list(k = k, mu = mu, se = se_mu, lcl = mu - z * se_mu, ucl = mu + z * se_mu,
       tau2 = tau2, i2 = i2, pi_lo = pi_lo, pi_hi = pi_hi, Q = Q, df = df)
}

run_pipeline <- function(nm) {
  fin <- read_safe(file.path(OUT_DIR, sprintf("bnb_final_%s.rds", nm)))
  imp <- read_safe(file.path(OUT_DIR, sprintf("bnb_imp_%s.rds",   nm)))
  if (is.null(fin) || is.null(imp)) {
    message(sprintf("[%s] skipped: bnb_final_ or bnb_imp_ is missing.", nm)); return(NULL)
  }

  imps_all <- reduce_imps(imp, fin$n0, fin$outlier_ids)
  ids      <- pick_imps(length(imps_all), IECV_N_IMP)
  imps     <- imps_all[ids]
  d1       <- as.data.frame(imps[[1]])

  if (!is.null(CLUSTER_DERIVE)) {
    cl <- factor(CLUSTER_DERIVE(d1))
  } else {
    if (!CLUSTER_VAR %in% names(d1)) {
      message(sprintf("[%s] skipped: cluster variable '%s' is not in the imputed data.",
                      nm, CLUSTER_VAR)); return(NULL)
    }
    cl <- factor(d1[[CLUSTER_VAR]])
  }
  cl_disp <- relabel_levels(CLUSTER_VAR, levels(cl))

  smf_txt <- fin$smformula_final
  smf     <- stats::as.formula(smf_txt)
  in_model <- CLUSTER_VAR %in% all.vars(smf)
  dropped  <- FALSE
  if (in_model) {
    if (identical(CLUSTER_IN_MODEL, "skip")) {
      message(sprintf(paste0("[%s] skipped: '%s' is in the final model, so its coefficient ",
                             "is not estimable when that cluster is held out ",
                             "(CLUSTER_IN_MODEL = 'skip')."), nm, CLUSTER_VAR))
      return(NULL)
    }
    smf <- stats::update(smf, stats::as.formula(sprintf(". ~ . - %s", CLUSTER_VAR)))
    dropped <- TRUE
    message(sprintf(paste0("[%s] '%s' removed from the model for this validation; ",
                           "its coefficient cannot be estimated when the cluster is ",
                           "held out. The validated model therefore has %d fewer terms ",
                           "than Table 2."), nm, CLUSTER_VAR, nlevels(cl) - 1L))
  }

  mats <- lapply(imps, function(d) {
    d <- as.data.frame(d)
    stats::model.matrix(stats::delete.response(stats::terms(smf, data = d)), data = d)
  })
  ys <- lapply(imps, function(d) as.integer(as.data.frame(d)$good))
  if (!all(vapply(mats, function(m) identical(colnames(m), colnames(mats[[1]])), logical(1))))
    stop(sprintf("[%s] the model matrix differs between imputations; cannot proceed.", nm))
  ncol_mm <- ncol(mats[[1]])

  pred_store <- list()
  rows <- rbindlist(lapply(levels(cl), function(lv) {
    idx  <- which(cl == lv)
    n_k  <- length(idx)
    y1   <- ys[[1]]
    ev_k <- sum(y1[idx]);  nev_k <- n_k - ev_k
    ev_t <- sum(y1[-idx]); nev_t <- length(y1) - n_k - ev_t

    base <- data.table(
      pipeline = nm, label = unname(PIPE_LABELS[[nm]]),
      cluster_variable = CLUSTER_VAR,
      cluster = lv, cluster_label = cl_disp[match(lv, levels(cl))],
      n_holdout = n_k, events_holdout = ev_k, nonevents_holdout = nev_k,
      n_train = length(y1) - n_k, events_train = ev_t, nonevents_train = nev_t,
      n_parameters = ncol_mm - 1L,
      cluster_var_dropped = dropped, m_used = length(imps))

    if (n_k < MIN_CLUSTER_N ||
        min(ev_k, nev_k) < MIN_CLUSTER_EVENTS ||
        min(ev_t, nev_t) < MIN_TRAIN_EVENTS) {
      message(sprintf("  [%s / %s] not validated: n = %d, events = %d, non-events = %d.",
                      nm, lv, n_k, ev_k, nev_k))
      return(cbind(base, data.table(validated = FALSE, estimators = NA_character_)))
    }

    per <- lapply(seq_along(imps), function(i) {
      x <- mats[[i]]; y <- ys[[i]]
      fit <- fit_train(x[-idx, , drop = FALSE], y[-idx], fin$path)
      if (is.null(fit)) return(NULL)
      lp <- as.numeric(x[idx, , drop = FALSE] %*% fit$coef[colnames(x)])
      c(perf_one(y[idx], lp), list(estimator = fit$estimator, lp = lp))
    })
    per <- Filter(Negate(is.null), per)
    if (!length(per))
      return(cbind(base, data.table(validated = FALSE, estimators = "no fit converged")))

    grab <- function(f) vapply(per, function(z) z[[f]], numeric(1))
    lp_bar <- Reduce(`+`, lapply(per, `[[`, "lp")) / length(per)
    pred_store[[lv]] <<- data.table(pipeline = nm, cluster = lv,
                                    cluster_label = cl_disp[match(lv, levels(cl))],
                                    row_in_subset = idx,
                                    y = ys[[1]][idx], lp = lp_bar,
                                    p = prob_from_lp(lp_bar))

    cst <- grab("c_stat"); cse <- grab("se_c_stat")
    lg_est <- logit(pmin(pmax(cst, 1e-6), 1 - 1e-6))
    lg_se  <- cse / (cst * (1 - cst))
    rc <- rubin1(lg_est, lg_se)
    rs <- rubin1(grab("cal_slope"),     grab("se_cal_slope"))
    ri <- rubin1(grab("cal_intercept"), grab("se_cal_int"))
    ro <- rubin1(grab("log_oe"),        grab("se_log_oe"))
    z  <- stats::qnorm(1 - (1 - CI_LEVEL) / 2)

    cbind(base, data.table(
      validated = TRUE,
      estimators = paste(sort(unique(vapply(per, `[[`, character(1), "estimator"))),
                         collapse = " / "),
      c_stat = expit(rc[1]), c_lcl = expit(rc[1] - z * rc[2]), c_ucl = expit(rc[1] + z * rc[2]),
      c_logit = rc[1], c_logit_se = rc[2],
      cal_slope = rs[1], cal_slope_lcl = rs[1] - z * rs[2], cal_slope_ucl = rs[1] + z * rs[2],
      cal_slope_se = rs[2],
      cal_intercept = ri[1], cal_int_lcl = ri[1] - z * ri[2], cal_int_ucl = ri[1] + z * ri[2],
      cal_int_se = ri[2],
      oe_ratio = exp(ro[1]), oe_lcl = exp(ro[1] - z * ro[2]), oe_ucl = exp(ro[1] + z * ro[2]),
      log_oe = ro[1], log_oe_se = ro[2],
      brier = mean(grab("brier")),
      observed_risk = mean(ys[[1]][idx]),
      expected_risk = mean(vapply(per, function(z) mean(prob_from_lp(z$lp)), numeric(1)))))
  }), fill = TRUE)

  if (length(pred_store))
    saveRDS(rbindlist(pred_store, fill = TRUE),
            file.path(OUT_DIR, sprintf("bnb_iecv_predictions_%s.rds", nm)))
  rows
}

message("==== Section 19: internal-external cross-validation ====")
message(sprintf("Clusters defined by '%s'; %s; estimator strategy '%s'.",
                CLUSTER_VAR,
                if (identical(CLUSTER_IN_MODEL, "drop"))
                  "the cluster variable is dropped from the model when it appears in it"
                else "pipelines containing the cluster variable are skipped",
                IECV_SEP_STRATEGY))

by_cluster <- rbindlist(Filter(Negate(is.null), lapply(PIPELINES, run_pipeline)), fill = TRUE)
if (!nrow(by_cluster))
  stop("nothing was validated; check that Sections 05 and 10-1 have been run and that ",
       "CLUSTER_VAR exists in the imputed data.")

saveRDS(by_cluster, file.path(OUT_DIR, "table_s10_iecv_by_cluster.rds"))
fwrite(by_cluster,  file.path(OUT_DIR, "table_s10_iecv_by_cluster.csv"))

message("\n---- 19.1 Performance in each held-out cluster ----")
print(by_cluster[validated == TRUE, .(
  pipeline, cluster = cluster_label, n = n_holdout, events = events_holdout,
  C = sprintf("%.3f (%.3f to %.3f)", c_stat, c_lcl, c_ucl),
  slope = sprintf("%.2f (%.2f to %.2f)", cal_slope, cal_slope_lcl, cal_slope_ucl),
  `O:E` = sprintf("%.2f (%.2f to %.2f)", oe_ratio, oe_lcl, oe_ucl),
  Brier = round(brier, 3), estimators)])

summarise_metric <- function(dt, est_col, se_col, back, metric_name) {
  m <- meta_dl(dt[[est_col]], dt[[se_col]])
  data.table(metric = metric_name, k_clusters = m$k,
             summary = back(m$mu), summary_lcl = back(m$lcl), summary_ucl = back(m$ucl),
             pi_lcl = back(m$pi_lo), pi_ucl = back(m$pi_hi),
             tau2 = m$tau2, I2 = m$i2)
}

summ <- rbindlist(lapply(split(by_cluster[validated == TRUE], by = "pipeline", drop = TRUE),
                         function(dt) {
  if (!nrow(dt)) return(NULL)
  cbind(data.table(pipeline = dt$pipeline[1], label = dt$label[1],
                   cluster_variable = CLUSTER_VAR,
                   cluster_var_dropped = dt$cluster_var_dropped[1]),
        rbindlist(list(
          summarise_metric(dt, "c_logit",      "c_logit_se",   expit,    "C statistic"),
          summarise_metric(dt, "cal_slope",    "cal_slope_se", identity, "Calibration slope"),
          summarise_metric(dt, "cal_intercept","cal_int_se",   identity, "Calibration-in-the-large"),
          summarise_metric(dt, "log_oe",       "log_oe_se",    exp,      "Observed : expected"))))
}), fill = TRUE)

saveRDS(summ, file.path(OUT_DIR, "table_s10_iecv_summary.rds"))
fwrite(summ,  file.path(OUT_DIR, "table_s10_iecv_summary.csv"))

message("\n---- 19.2 Random-effects summary and 95% prediction interval ----")
print(summ[, .(pipeline, metric, k_clusters,
               summary = sprintf("%.3f (%.3f to %.3f)", summary, summary_lcl, summary_ucl),
               prediction_interval = ifelse(is.finite(pi_lcl),
                                            sprintf("%.3f to %.3f", pi_lcl, pi_ucl), "—"),
               tau2 = signif(tau2, 3), I2 = round(I2, 3))])
message("The prediction interval, not the confidence interval, is the quantity to quote: ",
        "it is where performance in a new cluster of this kind is expected to fall. ",
        "It is undefined for fewer than three validated clusters.")

if (HAS_GG && nrow(by_cluster[validated == TRUE])) {
  library(ggplot2)

  mk <- function(metric_name, col_est, col_lo, col_hi, xintercept, xlab) {
    d <- by_cluster[validated == TRUE]
    dd <- data.table(pipeline = d$label,
                     row = sprintf("%s  (n = %d)", d$cluster_label, d$n_holdout),
                     est = d[[col_est]], lo = d[[col_lo]], hi = d[[col_hi]],
                     kind = "cluster")
    s <- summ[metric == metric_name]
    if (nrow(s)) {
      ds <- data.table(pipeline = s$label, row = "Random-effects summary",
                       est = s$summary, lo = s$summary_lcl, hi = s$summary_ucl,
                       kind = "summary")
      dp <- data.table(pipeline = s$label, row = "95% prediction interval",
                       est = s$summary, lo = s$pi_lcl, hi = s$pi_ucl,
                       kind = "prediction")
      dd <- rbindlist(list(dd, ds, dp), fill = TRUE)
    }
    dd[, row := factor(row, levels = rev(unique(row)))]
    ggplot(dd[is.finite(est)], aes(x = est, y = row, colour = kind, shape = kind)) +
      geom_vline(xintercept = xintercept, linetype = "dashed", colour = COL_REF) +
      geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.18, na.rm = TRUE) +
      geom_point(size = 2.1, na.rm = TRUE) +
      facet_wrap(~ pipeline, ncol = 2, scales = "free_y") +
      scale_colour_manual(values = c(cluster = COL_POINT, summary = COL_SUMMARY,
                                     prediction = COL_SUMMARY), guide = "none") +
      scale_shape_manual(values = c(cluster = 16, summary = 18, prediction = 124),
                         guide = "none") +
      labs(x = xlab, y = NULL,
           title = sprintf("Internal-external cross-validation: %s", tolower(metric_name)),
           subtitle = sprintf("Clusters defined by %s; each estimate comes from patients not used to fit the model",
                              relabel_vars(CLUSTER_VAR))) +
      theme_bw(base_size = 11) +
      theme(panel.grid.minor = element_blank(),
            plot.title = element_text(face = "bold", size = 12),
            strip.background = element_rect(fill = "grey92", colour = NA))
  }

  p_c <- mk("C statistic", "c_stat", "c_lcl", "c_ucl", 0.5, "C statistic (95% CI)")
  p_s <- mk("Calibration slope", "cal_slope", "cal_slope_lcl", "cal_slope_ucl", 1,
            "Calibration slope (95% CI)")

  n_rows <- length(unique(by_cluster[validated == TRUE]$cluster)) + 2
  h <- max(4.2, 1.2 + 0.42 * n_rows * 2)
  ggsave(file.path(FIG_DIR, "fig_iecv_forest.png"),
         p_c, width = 9.5, height = h, dpi = 300)
  ggsave(file.path(FIG_DIR, "fig_iecv_forest_slope.png"),
         p_s, width = 9.5, height = h, dpi = 300)
  message("\nSaved figures/fig_iecv_forest.png and figures/fig_iecv_forest_slope.png")

  preds <- rbindlist(lapply(PIPELINES, function(nm)
    read_safe(file.path(OUT_DIR, sprintf("bnb_iecv_predictions_%s.rds", nm)))), fill = TRUE)
  if (!is.null(preds) && nrow(preds)) {
    preds[, label := factor(unname(PIPE_LABELS[pipeline]),
                            levels = unname(PIPE_LABELS[PIPELINES]))]
    if (!"cluster_label" %in% names(preds)) {
      preds[, cluster_label := cluster]
      message("bnb_iecv_predictions_*.rds has no cluster_label column; the facet strips will ",
              "show the raw levels. Re-run section 19.1 to regenerate them.")
    }
    lab_lv <- if (!is.null(CLUSTER_ORDER))
      intersect(CLUSTER_ORDER, unique(as.character(preds$cluster_label)))
    else unique(as.character(preds$cluster_label))
    lab_lv <- c(lab_lv, setdiff(unique(as.character(preds$cluster_label)), lab_lv))
    preds[, cluster_label := factor(as.character(cluster_label), levels = lab_lv)]

    preds[, bin := cut(p, breaks = unique(stats::quantile(p, seq(0, 1, length.out = 6),
                                                          na.rm = TRUE)),
                       include.lowest = TRUE), by = .(pipeline, cluster_label)]
    bins <- preds[, .(x = mean(p), y = mean(y), n = .N), by = .(label, cluster_label, bin)]
    p_cal <- ggplot(preds, aes(x = p, y = y)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = COL_REF) +
      geom_smooth(method = "loess", se = TRUE, span = 1, colour = COL_POINT,
                  fill = COL_POINT, alpha = 0.15, linewidth = 0.7, na.rm = TRUE) +
      geom_point(data = bins, aes(x = x, y = y, size = n), inherit.aes = FALSE,
                 shape = 21, fill = "white", colour = "grey20", na.rm = TRUE) +
      scale_size_continuous(range = c(1.2, 3.4), guide = "none") +
      coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
      facet_grid(label ~ cluster_label) +
      labs(x = "Predicted probability of a good functional outcome",
           y = "Observed proportion",
           title = "Calibration in held-out clusters (internal-external cross-validation)",
           subtitle = sprintf(paste0("Clusters defined by %s. Predictions come from models ",
                                     "fitted without the patients shown in each panel."),
                              relabel_vars(CLUSTER_VAR))) +
      theme_bw(base_size = 10) +
      theme(panel.grid.minor = element_blank(),
            plot.title = element_text(face = "bold", size = 12),
            strip.background = element_rect(fill = "grey92", colour = NA))
    ggsave(file.path(FIG_DIR, "fig_iecv_calibration.png"), p_cal,
           width = 10, height = 8.5, dpi = 300)
    message("Saved figures/fig_iecv_calibration.png")
  }
} else if (!HAS_GG) {
  message("\nggplot2 is not installed; tables were written but no figure was drawn.")
}

message("\n==== Section 19 complete ====")
message("  data/table_s10_iecv_by_cluster.csv, .rds")
message("  data/table_s10_iecv_summary.csv, .rds")
message("  data/bnb_iecv_predictions_<pipeline>.rds")
message("  figures/fig_iecv_forest.png, fig_iecv_forest_slope.png, fig_iecv_calibration.png")
message("")
message("For the manuscript, the three C statistics must be labelled distinctly:")
message("  apparent            - fitted and evaluated on the same patients (Table S6)")
message("  optimism-corrected  - bootstrap, same patients, same case mix (Table S6)")
message("  internal-external   - evaluated on patients whose cluster was not used (this table)")
message("")
message("Limits to state in the legend, in this order:")
message("  1. This is not external validation. Every patient comes from one ward in one")
message("     hospital, so between-cluster variation understates between-hospital variation.")
message("  2. Variable selection, spline knots and centring constants were fixed on the")
message("     whole subset before this validation, so selection-induced optimism is not")
message("     removed; only the coefficients are re-estimated in each fold.")
if (any(by_cluster$cluster_var_dropped, na.rm = TRUE))
  message(sprintf(paste0("  3. '%s' was removed from the model wherever it was retained, because a",
                         " coefficient for the held-out cluster is not estimable. The validated",
                         " model is therefore not identical to the model in Table 2."),
                  CLUSTER_VAR))
message("  4. Clusters are unequal in size; the prediction interval reflects that and is",
        " the interval to quote.")
