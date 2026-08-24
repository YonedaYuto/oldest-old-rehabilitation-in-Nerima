library(data.table)
library(here)
library(ggplot2)

need_pkg <- function(p) if (!requireNamespace(p, quietly = TRUE))
  stop(sprintf("package '%s' is required. install.packages('%s')", p, p))
HAS_LOGISTF <- requireNamespace("logistf", quietly = TRUE)
HAS_PROC    <- requireNamespace("pROC",    quietly = TRUE)

LABELS_PATH <- here::here("00_labels.R")
if (!exists("relabel_vars")) {
  if (!file.exists(LABELS_PATH)) stop("00_labels.R not found: ", LABELS_PATH)
  source(LABELS_PATH)
}

OUT_DIR <- here::here("data")
FIG_DIR <- here::here("figures")
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

CI_LEVEL   <- 0.95
PROB_EPS   <- 1e-8
N_BOOT     <- 200L
M_BOOT     <- 10L
P_CRIT_MAIN <- 0.05
P_CRIT_INT  <- 0.05
BOOT_SEP_STRATEGY <- "adaptive"
BETA_MAX   <- 15; SE_MAX <- 30; FITTED_EPS <- 1e-6

PIPELINES  <- names(PIPE_LABELS)
set.seed(20240601)

read_safe <- function(f) if (file.exists(f)) readRDS(f) else NULL

read_fixed_correction <- function() {
  perf <- read_safe(file.path(OUT_DIR, "bnb_performance_summary.rds"))
  if (is.null(perf)) {
    message("Section 12 summary (data/bnb_performance_summary.rds) not found; ",
            "the formula-fixed correction cannot be shown. Run Section 12 first.")
    return(NULL)
  }
  pf <- as.data.table(perf)
  if ("model" %in% names(pf)) {
    sel <- tolower(as.character(pf$model)) == "final"
    if (!any(sel)) {
      message("No rows with model == 'Final' in bnb_performance_summary.rds; ",
              "models present: ", paste(unique(as.character(pf$model)), collapse = ", "))
      return(NULL)
    }
    pf <- pf[sel]
  }
  if ("pipeline" %in% names(pf)) pf[, pipeline := as.character(pipeline)]
  slope_col <- intersect(c("cal_slope_corrected", "slope_corrected"), names(pf))
  auc_col   <- intersect(c("AUC_corrected"), names(pf))
  if (!length(auc_col)) {
    message("AUC_corrected is absent from bnb_performance_summary.rds ",
            "(Section 12 was probably run with DO_BOOTSTRAP = FALSE).")
    return(NULL)
  }
  out <- pf[, .(pipeline,
                AUC_corrected_fixed   = as.numeric(get(auc_col[1])),
                slope_corrected_fixed = if (length(slope_col))
                  as.numeric(get(slope_col[1])) else NA_real_)]
  if (all(!is.finite(out$AUC_corrected_fixed)))
    message("AUC_corrected in Section 12 is all NA (DO_BOOTSTRAP = FALSE?); ",
            "the formula-fixed correction will be blank.")
  unique(out, by = "pipeline")
}

reduce_imps <- function(obj, n0, outlier_ids) {
  imps_full <- obj$smcfcs$impDatasets
  bad <- sort(unique(outlier_ids)); bad <- bad[bad %in% seq_len(n0)]
  if (length(bad) > 0) lapply(imps_full, function(d) d[-bad, , drop = FALSE]) else imps_full
}

pick_imps <- function(m, k) {
  if (is.null(k) || k >= m) return(seq_len(m))
  unique(round(seq(1, m, length.out = k)))
}

glm_unstable <- function(fit) {
  if (is.null(fit)) return(TRUE)
  if (!isTRUE(fit$converged)) return(TRUE)
  b <- coef(fit); if (any(!is.finite(b)) || any(abs(b) > BETA_MAX)) return(TRUE)
  se <- tryCatch(sqrt(diag(vcov(fit))), error = function(e) NA_real_)
  if (any(!is.finite(se)) || any(se > SE_MAX)) return(TRUE)
  fv <- fit$fitted.values
  any(fv < FITTED_EPS | fv > 1 - FITTED_EPS)
}

firth_pll <- function(fit) {
  ll <- fit$loglik
  if (!is.null(names(ll)) && "full" %in% names(ll)) return(unname(ll[["full"]]))
  ll[length(ll)]
}

mk_formula <- function(labs) {
  if (length(labs) == 0) return(stats::as.formula("good ~ 1"))
  stats::reformulate(labs, response = "good")
}

fit_robust <- function(smf, d, base_path, strategy = BOOT_SEP_STRATEGY) {
  always_firth <- identical(base_path, "firth") || identical(strategy, "firth")
  if (always_firth) {
    if (!HAS_LOGISTF) stop("Firth path requested but logistf is not installed.")
    fit <- tryCatch(logistf::logistf(smf, data = d, pl = FALSE), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    b <- coef(fit); if (is.null(names(b))) names(b) <- rownames(vcov(fit))
    return(list(coef = b, estimator = "firth", ll = firth_pll(fit), k = length(b)))
  }
  fit <- tryCatch(suppressWarnings(glm(smf, family = binomial(), data = d)),
                  error = function(e) NULL)
  if (!is.null(fit) && !glm_unstable(fit))
    return(list(coef = coef(fit), estimator = "glm",
                ll = as.numeric(stats::logLik(fit)), k = length(coef(fit))))
  if (!HAS_LOGISTF) return(if (is.null(fit)) NULL else
    list(coef = coef(fit), estimator = "glm(unstable)",
         ll = as.numeric(stats::logLik(fit)), k = length(coef(fit))))
  ff <- tryCatch(logistf::logistf(smf, data = d, pl = FALSE), error = function(e) NULL)
  if (is.null(ff)) return(NULL)
  b <- coef(ff); if (is.null(names(b))) names(b) <- rownames(vcov(ff))
  list(coef = b, estimator = "firth", ll = firth_pll(ff), k = length(b))
}

predict_lp <- function(coef_named, smf, newdata) {
  mm <- tryCatch(stats::model.matrix(stats::delete.response(stats::terms(smf)),
                                     data = newdata), error = function(e) NULL)
  if (is.null(mm)) return(rep(NA_real_, nrow(newdata)))
  common <- intersect(colnames(mm), names(coef_named))
  if (length(common) == 0) return(rep(NA_real_, nrow(newdata)))
  as.numeric(mm[, common, drop = FALSE] %*% coef_named[common])
}

prob_from_lp <- function(lp) pmin(pmax(stats::plogis(lp), PROB_EPS), 1 - PROB_EPS)

auc_one <- function(y, p) {
  pos <- p[y == 1]; neg <- p[y == 0]
  if (length(pos) == 0 || length(neg) == 0) return(NA_real_)
  if (HAS_PROC) {
    r <- tryCatch(pROC::roc(response = y, predictor = p, quiet = TRUE,
                            levels = c(0, 1), direction = "<"), error = function(e) NULL)
    if (!is.null(r)) return(as.numeric(pROC::auc(r)))
  }
  rr <- rank(c(pos, neg))
  (sum(rr[seq_along(pos)]) - length(pos) * (length(pos) + 1) / 2) /
    (length(pos) * length(neg))
}

cal_slope <- function(y, lp) {
  fs <- tryCatch(suppressWarnings(glm(y ~ lp, family = binomial())), error = function(e) NULL)
  if (is.null(fs)) NA_real_ else unname(coef(fs)["lp"])
}

group_key <- function(lab) {
  if (grepl("_measurable$", lab)) return(lab)
  if (grepl("_c[0-9]+$", lab))    return(sub("_c[0-9]+$", "", lab))
  if (grepl("_c$", lab))          return(sub("_c$", "", lab))
  lab
}

backward_lrt <- function(groups, current, d, p_crit = P_CRIT_MAIN,
                         keep = character(0), base_path = "glm") {
  repeat {
    cand <- setdiff(current, keep)
    if (length(cand) == 0) break
    full_labs <- unlist(groups[current], use.names = FALSE)
    rf <- fit_robust(mk_formula(full_labs), d, base_path); if (is.null(rf)) break
    pv <- vapply(cand, function(g) {
      red <- unlist(groups[setdiff(current, g)], use.names = FALSE)
      rr  <- fit_robust(mk_formula(red), d, base_path); if (is.null(rr)) return(NA_real_)
      lr  <- 2 * (rf$ll - rr$ll); k <- rf$k - rr$k
      if (k <= 0 || !is.finite(lr)) return(NA_real_)
      stats::pchisq(max(lr, 0), df = k, lower.tail = FALSE)
    }, numeric(1))
    worst <- which.max(pv)
    if (length(worst) == 0 || !is.finite(pv[worst]) || pv[worst] < p_crit) break
    current <- setdiff(current, cand[worst])
  }
  current
}

int_labels <- function(groups, a, b)
  as.vector(outer(groups[[a]], groups[[b]], paste, sep = ":"))

forward_int <- function(groups, main_keys, pairs, d, p_crit = P_CRIT_INT,
                        base_path = "glm") {
  selected <- list()
  base_labs <- unlist(groups[main_keys], use.names = FALSE)
  repeat {
    cand <- pairs[vapply(pairs, function(pr) all(pr %in% main_keys) &&
                           !paste(pr, collapse = ":") %in% names(selected), logical(1))]
    if (length(cand) == 0) break
    cur_labs <- c(base_labs, unlist(selected, use.names = FALSE))
    r0 <- fit_robust(mk_formula(cur_labs), d, base_path); if (is.null(r0)) break
    pv <- vapply(cand, function(pr) {
      add <- int_labels(groups, pr[1], pr[2])
      r1  <- fit_robust(mk_formula(c(cur_labs, add)), d, base_path)
      if (is.null(r1)) return(NA_real_)
      lr <- 2 * (r1$ll - r0$ll); k <- r1$k - r0$k
      if (k <= 0 || !is.finite(lr)) return(NA_real_)
      stats::pchisq(max(lr, 0), df = k, lower.tail = FALSE)
    }, numeric(1))
    best <- which.min(pv)
    if (length(best) == 0 || !is.finite(pv[best]) || pv[best] >= p_crit) break
    pr <- cand[[best]]
    selected[[paste(pr, collapse = ":")]] <- int_labels(groups, pr[1], pr[2])
  }
  c(base_labs, unlist(selected, use.names = FALSE))
}

optimism_with_selection <- function(imps, groups, full_keys, keep_keys, pairs,
                                    smf_final, path, n_boot = N_BOOT, m_use = M_BOOT) {
  m <- length(imps); use_idx <- pick_imps(m, m_use)
  per_imp <- vector("list", length(use_idx))

  for (j in seq_along(use_idx)) {
    i <- use_idx[j]
    d <- imps[[i]]; y <- as.integer(d$good); n <- nrow(d)

    r_org <- fit_robust(smf_final, d, path); if (is.null(r_org)) next
    p_org <- prob_from_lp(predict_lp(r_org$coef, smf_final, d))
    a_app <- auc_one(y, p_org); s_app <- cal_slope(y, stats::qlogis(p_org))

    o_a <- o_s <- numeric(0); n_terms <- integer(0)
    for (b in seq_len(n_boot)) {
      idx <- sample.int(n, n, replace = TRUE)
      db  <- d[idx, , drop = FALSE]

      keys_b <- tryCatch(backward_lrt(groups, full_keys, db, P_CRIT_MAIN,
                                      keep_keys, path), error = function(e) NULL)
      if (is.null(keys_b) || length(keys_b) == 0) next
      labs_b <- tryCatch(forward_int(groups, keys_b, pairs, db, P_CRIT_INT, path),
                         error = function(e) NULL)
      if (is.null(labs_b) || length(labs_b) == 0) next
      smf_b <- mk_formula(labs_b)

      rb <- fit_robust(smf_b, db, path); if (is.null(rb)) next
      pbb <- prob_from_lp(predict_lp(rb$coef, smf_b, db))
      a_b <- auc_one(as.integer(db$good), pbb)
      s_b <- cal_slope(as.integer(db$good), stats::qlogis(pbb))
      pbo <- prob_from_lp(predict_lp(rb$coef, smf_b, d))
      a_t <- auc_one(y, pbo)
      s_t <- cal_slope(y, stats::qlogis(pbo))
      if (is.finite(a_b) && is.finite(a_t)) o_a <- c(o_a, a_b - a_t)
      if (is.finite(s_b) && is.finite(s_t)) o_s <- c(o_s, s_b - s_t)
      n_terms <- c(n_terms, length(labs_b))
    }
    per_imp[[j]] <- data.table(
      .imp = i, AUC_apparent = a_app, slope_apparent = s_app,
      AUC_optimism = mean(o_a, na.rm = TRUE), slope_optimism = mean(o_s, na.rm = TRUE),
      n_terms_mean = mean(n_terms), n_ok = length(o_a))
  }

  pi <- rbindlist(Filter(Negate(is.null), per_imp), fill = TRUE)
  if (nrow(pi) == 0) return(NULL)

  mc <- function(x) if (length(x) > 1) stats::sd(x) / sqrt(length(x)) else NA_real_
  auc_corr_i   <- pi$AUC_apparent   - pi$AUC_optimism
  slope_corr_i <- pi$slope_apparent - pi$slope_optimism

  list(per_imp = pi, summary = data.table(
    AUC_apparent      = mean(pi$AUC_apparent),
    AUC_optimism      = mean(pi$AUC_optimism,   na.rm = TRUE),
    AUC_corrected     = mean(auc_corr_i,        na.rm = TRUE),
    AUC_corrected_mcse   = mc(auc_corr_i[is.finite(auc_corr_i)]),
    slope_apparent    = mean(pi$slope_apparent),
    slope_optimism    = mean(pi$slope_optimism, na.rm = TRUE),
    slope_corrected   = mean(slope_corr_i,      na.rm = TRUE),
    slope_corrected_mcse = mc(slope_corr_i[is.finite(slope_corr_i)]),
    n_boot = n_boot, m_used = nrow(pi),
    mean_terms_selected = mean(pi$n_terms_mean, na.rm = TRUE)))
}

process_pipeline <- function(nm) {
  prov <- read_safe(file.path(OUT_DIR, sprintf("bnb_provisional_%s.rds", nm)))
  ix   <- read_safe(file.path(OUT_DIR, sprintf("bnb_interactions_%s.rds", nm)))
  fin  <- read_safe(file.path(OUT_DIR, sprintf("bnb_final_%s.rds", nm)))
  imp  <- read_safe(file.path(OUT_DIR, sprintf("bnb_imp_%s.rds", nm)))
  if (is.null(prov) || is.null(fin) || is.null(imp)) {
    message(sprintf("[%s] inputs missing; skipped.", nm)); return(NULL)
  }
  imps <- reduce_imps(imp, fin$n0, fin$outlier_ids)

  groups    <- fin$groups
  full_labs <- attr(stats::terms(stats::as.formula(prov$smformula_full)), "term.labels")
  full_keys <- unique(vapply(full_labs, group_key, character(1)))
  full_keys <- intersect(full_keys, names(groups))
  keep_keys <- intersect(prov$keep_vars, names(groups))

  pairs <- list()
  if (!is.null(ix) && !is.null(ix$candidates)) {
    cd <- as.data.table(ix$candidates)
    pc <- intersect(c("var1", "var2"), names(cd))
    if (length(pc) == 2) {
      pairs <- lapply(seq_len(nrow(cd)),
                      function(i) as.character(unlist(cd[i, ..pc], use.names = FALSE)))
    } else if ("key" %in% names(cd)) {
      pairs <- strsplit(as.character(cd$key), ":", fixed = TRUE)
    }
  }
  if (length(pairs) == 0 && length(fin$candidate_int_keys))
    pairs <- lapply(strsplit(fin$candidate_int_keys, ":", fixed = TRUE), identity)
  pairs <- Filter(function(pr) length(pr) == 2 && all(pr %in% names(groups)), pairs)

  message(sprintf("[%s] candidates: %d main-effect groups, %d interaction pairs; path = %s",
                  nm, length(full_keys), length(pairs), fin$path))

  res <- optimism_with_selection(imps, groups, full_keys, keep_keys, pairs,
                                 stats::as.formula(fin$smformula_final), fin$path)
  if (is.null(res)) { message(sprintf("[%s] optimism could not be estimated.", nm)); return(NULL) }

  saveRDS(c(res, list(pipeline = nm)),
          file.path(OUT_DIR, sprintf("bnb_optimism_selection_%s.rds", nm)))
  cbind(data.table(pipeline = nm, label = PIPE_LABELS[[nm]], path = fin$path), res$summary)
}

summ <- rbindlist(Filter(Negate(is.null), lapply(PIPELINES, process_pipeline)), fill = TRUE)
if (!nrow(summ))
  stop("No pipeline produced an optimism estimate. Run Sections 8, 9, 10 and 5 first.")

fixed <- read_fixed_correction()
if (!is.null(fixed) && nrow(summ)) {
  summ <- merge(summ, fixed, by = "pipeline", all.x = TRUE, sort = FALSE)
} else {
  summ[, `:=`(AUC_corrected_fixed = NA_real_, slope_corrected_fixed = NA_real_)]
}
summ[, `:=`(
  AUC_optimism_from_selection   = AUC_corrected_fixed   - AUC_corrected,
  slope_optimism_from_selection = slope_corrected_fixed - slope_corrected)]
setcolorder(summ, c("pipeline", "label", "path",
                    "AUC_apparent", "AUC_corrected_fixed", "AUC_corrected",
                    "AUC_corrected_mcse", "AUC_optimism_from_selection"))

saveRDS(summ, file.path(OUT_DIR, "bnb_optimism_selection_summary.rds"))
data.table::fwrite(summ, file.path(OUT_DIR, "bnb_optimism_selection_summary.csv"))

message("---- Optimism correction: formula fixed (Section 12) vs selection included (12-2-1) ----")
print(summ[, .(pipeline,
               AUC_apparent    = round(AUC_apparent, 3),
               AUC_fixed       = round(AUC_corrected_fixed, 3),
               AUC_selection   = round(AUC_corrected, 3),
               AUC_mcse        = signif(AUC_corrected_mcse, 2),
               dAUC_selection  = round(AUC_optimism_from_selection, 3),
               slope_fixed     = round(slope_corrected_fixed, 3),
               slope_selection = round(slope_corrected, 3),
               slope_mcse      = signif(slope_corrected_mcse, 2))])
message("AUC_fixed      = Section 12, model formula held fixed (read from ",
        "data/bnb_performance_summary.rds; not recomputed here).")
message("AUC_selection  = this script, backward elimination and interaction ",
        "screening repeated inside every bootstrap replicate.")
message("dAUC_selection = AUC_fixed - AUC_selection = the part of the optimism ",
        "attributable to variable and interaction selection.")

if (nrow(summ)) {
  pl <- melt(summ[, .(pipeline, label, AUC_apparent,
                      AUC_corrected_fixed, AUC_corrected_selection = AUC_corrected)],
             id.vars = c("pipeline", "label"),
             measure.vars = c("AUC_apparent", "AUC_corrected_fixed", "AUC_corrected_selection"),
             variable.name = "kind", value.name = "AUC")
  pl[, kind := factor(kind, levels = c("AUC_apparent", "AUC_corrected_fixed",
                                       "AUC_corrected_selection"),
                      labels = c("Apparent", "Corrected (formula fixed)",
                                 "Corrected (selection included)"))]
  if (all(is.na(pl[kind == "Corrected (formula fixed)", AUC])))
    warning("The formula-fixed correction is unavailable; the figure will show ",
            "two bars per pipeline. Run Section 12 with DO_BOOTSTRAP = TRUE.")
  pl[, lab_txt := fifelse(is.finite(AUC), sprintf("%.3f", AUC), "n/a")]
  g <- ggplot(pl, aes(x = label, y = AUC, fill = kind)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7, na.rm = TRUE) +
    geom_text(aes(label = lab_txt, y = pmax(AUC, 0.5, na.rm = TRUE)),
              position = position_dodge(width = 0.8), vjust = -0.4, size = 3, na.rm = TRUE) +
    scale_fill_manual(values = c("#999999", "#4393C3", "#B2182B"), name = NULL,
                      drop = FALSE) +
    coord_cartesian(ylim = c(0.5, 1.02)) +
    labs(x = NULL, y = "AUC",
         title = "Optimism correction with and without the selection process") +
    theme_bw(base_size = 11) + theme(legend.position = "top")
  ggsave(file.path(FIG_DIR, "fig_optimism_selection.png"), g,
         width = 9, height = 5, dpi = 300)
}

message("Section 12-2-1 complete.")
message("  data/bnb_optimism_selection_<pipeline>.rds")
message("  data/bnb_optimism_selection_summary.rds, .csv")
message("  figures/fig_optimism_selection.png")
