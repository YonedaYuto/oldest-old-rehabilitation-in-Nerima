library(data.table)
library(here)

OUT_DIR <- here::here("data")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

LABELS_PATH <- here::here("00_labels.R")
if (!exists("relabel_vars")) {
  if (!file.exists(LABELS_PATH)) stop("00_labels.R not found: ", LABELS_PATH)
  source(LABELS_PATH)
}

HAS_LOGISTF <- requireNamespace("logistf", quietly = TRUE)
PIPELINES   <- names(PIPE_LABELS)

MC_REL_SE_TOL <- 0.05
M_RULE        <- function(fmi) ceiling(100 * fmi)

LOGISTF_PL <- FALSE

BETA_MAX <- 20; SE_MAX <- 10; FITTED_EPS <- 1e-8

read_safe <- function(f) if (file.exists(f)) readRDS(f) else NULL

reduce_imps <- function(obj, n0, outlier_ids) {
  imps_full <- obj$smcfcs$impDatasets
  bad <- sort(unique(outlier_ids)); bad <- bad[bad %in% seq_len(n0)]
  if (length(bad) > 0) lapply(imps_full, function(d) d[-bad, , drop = FALSE]) else imps_full
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

fit_with_se <- function(smf, d, base_path) {
  firth_fit <- function() {
    if (!HAS_LOGISTF) stop("logistf is required for the Firth path but is not installed.")
    ff <- tryCatch(logistf::logistf(smf, data = d, pl = LOGISTF_PL), error = function(e) NULL)
    if (is.null(ff)) return(NULL)
    b <- coef(ff); V <- as.matrix(vcov(ff))
    if (is.null(names(b))) names(b) <- rownames(V)
    list(coef = b, se = sqrt(diag(V)), estimator = "firth")
  }
  if (identical(base_path, "firth")) return(firth_fit())

  fit <- tryCatch(suppressWarnings(glm(smf, family = binomial(), data = d)),
                  error = function(e) NULL)
  if (!is.null(fit) && !glm_unstable(fit))
    return(list(coef = coef(fit), se = sqrt(diag(as.matrix(vcov(fit)))), estimator = "glm"))
  ff <- firth_fit()
  if (!is.null(ff)) ff$estimator <- "firth (fallback)"
  ff
}

pool_full <- function(coefs, ses, n_complete) {
  m     <- length(coefs)
  terms <- Reduce(intersect, lapply(coefs, names))
  if (length(terms) == 0) stop("no coefficient is present in every imputation")
  getb <- function(t) vapply(coefs, function(z) unname(z[[t]]), numeric(1))
  getu <- function(t) vapply(ses,   function(z) unname(z[[t]])^2, numeric(1))

  Qbar <- vapply(terms, function(t) mean(getb(t)), numeric(1))
  B    <- vapply(terms, function(t) stats::var(getb(t)), numeric(1))
  Ubar <- vapply(terms, function(t) mean(getu(t)), numeric(1))

  Tvar   <- Ubar + (1 + 1 / m) * B
  riv    <- (1 + 1 / m) * B / Ubar
  lambda <- (1 + 1 / m) * B / Tvar

  k      <- length(terms)
  dfc    <- max(n_complete - k, 1)
  df_old <- ifelse(lambda > 0, (m - 1) / lambda^2, Inf)
  df_obs <- ((dfc + 1) / (dfc + 3)) * dfc * (1 - lambda)
  df     <- df_old * df_obs / (df_old + df_obs)
  fmi    <- (riv + 2 / (df + 3)) / (riv + 1)

  se       <- sqrt(Tvar)
  mcse     <- sqrt(B / m)
  mcse_rel <- mcse / se

  data.table(term = terms, estimate = Qbar, se = se, B = B, Ubar = Ubar,
             riv = riv, lambda = lambda, df = df, fmi = fmi,
             mcse = mcse, mcse_rel = mcse_rel, m_needed = M_RULE(fmi))
}

terms_all <- list()
summ_all  <- list()

for (nm in PIPELINES) {
  fin <- read_safe(file.path(OUT_DIR, sprintf("bnb_final_%s.rds", nm)))
  imp <- read_safe(file.path(OUT_DIR, sprintf("bnb_imp_%s.rds",   nm)))
  if (is.null(fin) || is.null(imp)) {
    message(sprintf("[%s] skipped: bnb_final_ or bnb_imp_ is missing.", nm))
    next
  }

  imps <- reduce_imps(imp, fin$n0, fin$outlier_ids)
  smf  <- stats::as.formula(fin$smformula_final)

  coefs <- list(); ses <- list(); estimators <- character(0)
  for (d in imps) {
    r <- fit_with_se(smf, as.data.frame(d), fin$path)
    if (is.null(r)) next
    coefs[[length(coefs) + 1L]] <- r$coef
    ses[[length(ses)   + 1L]]   <- r$se
    estimators <- c(estimators, r$estimator)
  }
  if (length(coefs) < 2L) {
    message(sprintf("[%s] skipped: fewer than two imputations could be fitted.", nm))
    next
  }

  n_complete <- nrow(as.data.frame(imps[[1]]))
  tt <- pool_full(coefs, ses, n_complete)
  tt[, `:=`(pipeline = nm, label = unname(PIPE_LABELS[[nm]]), path = fin$path)]
  setcolorder(tt, c("pipeline", "label", "path", "term"))
  terms_all[[nm]] <- tt

  m_used   <- length(coefs)
  m_needed <- max(tt$m_needed, na.rm = TRUE)
  summ_all[[nm]] <- data.table(
    pipeline     = nm,
    label        = unname(PIPE_LABELS[[nm]]),
    path         = fin$path,
    estimators   = paste(sort(unique(estimators)), collapse = " / "),
    n            = n_complete,
    n0           = fin$n0,
    n_excluded   = fin$n0 - n_complete,
    m            = m_used,
    m_available  = length(imps),
    n_terms      = nrow(tt),
    fmi_min      = min(tt$fmi,      na.rm = TRUE),
    fmi_max      = max(tt$fmi,      na.rm = TRUE),
    mcse_max     = max(tt$mcse,     na.rm = TRUE),
    mcse_rel_max = max(tt$mcse_rel, na.rm = TRUE),
    m_needed     = m_needed,
    m_sufficient = (m_used >= m_needed) &&
                   (max(tt$mcse_rel, na.rm = TRUE) <= MC_REL_SE_TOL))

  message(sprintf("[%s] %s | n = %d (%d excluded), m = %d | FMI %.3f-%.3f | MCSE max %.4f (rel %.4f)",
                  nm, fin$path, n_complete, fin$n0 - n_complete, m_used,
                  min(tt$fmi, na.rm = TRUE), max(tt$fmi, na.rm = TRUE),
                  max(tt$mcse, na.rm = TRUE), max(tt$mcse_rel, na.rm = TRUE)))
}

if (length(summ_all) == 0) stop("nothing was computed; check that Sections 05 and 10 have been run.")

terms_dt <- rbindlist(terms_all, fill = TRUE)
summ_dt  <- rbindlist(summ_all,  fill = TRUE)

saveRDS(terms_dt, file.path(OUT_DIR, "bnb_fmi_mcse_terms.rds"))
saveRDS(summ_dt,  file.path(OUT_DIR, "bnb_fmi_mcse_summary.rds"))
data.table::fwrite(terms_dt, file.path(OUT_DIR, "bnb_fmi_mcse_terms.csv"))
data.table::fwrite(summ_dt,  file.path(OUT_DIR, "bnb_fmi_mcse_summary.csv"))

message("\n==== Section 05-1-1 summary ====")
print(summ_dt[, .(pipeline, path, n, m, n_terms,
                  fmi_min = round(fmi_min, 3), fmi_max = round(fmi_max, 3),
                  mcse_max = signif(mcse_max, 3), mcse_rel_max = round(mcse_rel_max, 4),
                  m_needed, m_sufficient)])

fmi_lo <- min(summ_dt$fmi_min, na.rm = TRUE)
fmi_hi <- max(summ_dt$fmi_max, na.rm = TRUE)
mcse_hi <- max(summ_dt$mcse_max, na.rm = TRUE)

message("\nFor the manuscript (Methods, Statistical analysis; legend of Supplementary Figure 3):")
message(sprintf("  fraction of missing information ranged from %.2f to %.2f", fmi_lo, fmi_hi))
message(sprintf("  Monte Carlo standard error of the pooled log odds ratios was at most %.3f", mcse_hi))
message(sprintf("  (largest relative Monte Carlo error %.3f, against the 0.05 convention)",
                max(summ_dt$mcse_rel_max, na.rm = TRUE)))

message("\nSection 05-1-1 complete. Written to:")
message("  data/bnb_fmi_mcse_terms.csv, .rds")
message("  data/bnb_fmi_mcse_summary.csv, .rds")
message("Section 05-1 remains the diagnostic for imputation convergence and for the",
        " observed-versus-imputed distributions; this script supplies the FMI and the",
        " Monte Carlo standard error of the FINAL models, which Section 05-1 does not.")
