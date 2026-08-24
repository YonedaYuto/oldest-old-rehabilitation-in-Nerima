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
CI_LEVEL   <- 0.95
THR_GOOD   <- 65; THR_SEVERE <- 26; THR_ELDER <- 90
PIPELINES  <- names(PIPE_LABELS)
read_safe  <- function(f) if (file.exists(f)) readRDS(f) else NULL
`%||%`     <- function(a, b) if (is.null(a)) b else a
set.seed(20240601)

S5B_SAMPLE <- "analysed"

ADJ_CI_SOURCE <- "stored"

DIGITS_OR <- 2L

if (!S5B_SAMPLE %in% c("analysed", "all"))
  stop("S5B_SAMPLE must be 'analysed' or 'all'.")
if (!ADJ_CI_SOURCE %in% c("table2_compatible", "stored"))
  stop("ADJ_CI_SOURCE must be 'table2_compatible' or 'stored'.")

if (!exists("BNB_small")) stop("BNB_small not found. Load it first.")
dat <- as.data.table(BNB_small)
if (!"good" %in% names(dat)) dat[, good := as.integer(as.numeric(mFIM_out) >= THR_GOOD)]
if (!"severe"  %in% names(dat)) dat[, severe  := as.numeric(mFIM_in) <= THR_SEVERE]
if (!"elderly" %in% names(dat)) dat[, elderly := as.numeric(age)     >= THR_ELDER]

fmt_or_ci <- function(est, lo, hi, digits = DIGITS_OR)
  sprintf("%.*f (%.*f to %.*f)", digits, est, digits, lo, digits, hi)

reduce_imps <- function(obj, n0, outlier_ids) {
  imps_full <- obj$smcfcs$impDatasets
  bad <- sort(unique(outlier_ids)); bad <- bad[bad %in% seq_len(n0)]
  if (length(bad) > 0) lapply(imps_full, function(d) d[-bad, , drop = FALSE]) else imps_full
}

C_ANTICIPATED <- 0.75

r2cs_from_c <- function(cstat, phi, n = 2e5) {
  sigma <- sqrt(2) * stats::qnorm(cstat)
  a <- stats::uniroot(function(a0) mean(stats::plogis(a0 + stats::rnorm(2e4, 0, sigma))) - phi,
                      interval = c(-15, 15))$root
  lp <- a + stats::rnorm(n, 0, sigma)
  p  <- stats::plogis(lp)
  y  <- stats::rbinom(n, 1, p)
  ll_m <- sum(y * log(p) + (1 - y) * log(1 - p))
  pbar <- mean(y)
  ll_0 <- sum(y * log(pbar) + (1 - y) * log(1 - pbar))
  1 - exp((2 / n) * (ll_0 - ll_m))
}
riley_min_n <- function(p_params, phi, r2cs) {
  max_r2cs <- 1 - (phi^phi * (1 - phi)^(1 - phi))^2
  n1 <- p_params / ((0.9 - 1) * log(1 - r2cs / 0.9))
  s2 <- r2cs / (r2cs + 0.05 * max_r2cs)
  n2 <- p_params / ((s2 - 1) * log(1 - r2cs / s2))
  n3 <- (1.96 / 0.05)^2 * phi * (1 - phi)
  n  <- ceiling(max(n1, n2, n3))
  data.table(n_criterion_i = ceiling(n1), n_criterion_ii = ceiling(n2),
             n_criterion_iii = ceiling(n3), n_required = n,
             events_required = ceiling(n * phi),
             epp_required = round(n * phi / p_params, 1),
             R2_cs_anticipated = round(r2cs, 4), max_R2_cs = round(max_r2cs, 4))
}

ss_rows <- lapply(PIPELINES, function(nm) {
  fin <- read_safe(file.path(OUT_DIR, sprintf("bnb_final_%s.rds", nm)))
  if (is.null(fin)) return(NULL)
  sub  <- if (grepl("severe", nm)) dat[severe == TRUE] else dat[elderly == TRUE]
  phi  <- mean(sub$good, na.rm = TRUE)
  p_par <- if (!is.null(fin$pooled) && "term" %in% names(as.data.table(fin$pooled)))
    sum(as.data.table(fin$pooled)$term != "(Intercept)")
  else length(fin$final_labels)
  if (p_par != length(fin$final_labels))
    message(sprintf("  [%s] %d parameters from %d formula terms (multi-level factor present).",
                    nm, p_par, length(fin$final_labels)))
  r2cs <- r2cs_from_c(C_ANTICIPATED, phi)

  n_av    <- nrow(sub)
  ev_av   <- sum(sub$good, na.rm = TRUE)
  nev_av  <- n_av - ev_av
  min_av  <- min(ev_av, nev_av)

  cbind(data.table(pipeline = nm, label = PIPE_LABELS[[nm]],
                   n_available = n_av,
                   events_available = ev_av,
                   nonevents_available = nev_av,
                   minority_class_available = min_av,
                   minority_class_is = if (ev_av <= nev_av) "events (good outcome)"
                                       else "non-events (poor outcome)",
                   outcome_proportion = round(phi, 3),
                   p_parameters = p_par,
                   epp_events_available   = round(ev_av  / p_par, 1),
                   epp_minority_available = round(min_av / p_par, 1),
                   epp_denominator_note = paste(
                     "epp_events_available uses all outcome events;",
                     "epp_minority_available uses the smaller outcome class,",
                     "which is the convention quoted in the manuscript and in",
                     "Supplementary Table S2.")),
        riley_min_n(p_par, phi, r2cs))
})
tab_ss <- rbindlist(Filter(Negate(is.null), ss_rows), fill = TRUE)
fwrite(tab_ss, file.path(OUT_DIR, "table_s2_sample_size.csv"))
message("---- 15.1 Minimum sample size (Riley et al.) ----")
print(tab_ss[, .(pipeline, n_available, events_available, nonevents_available,
                 minority_class_available, p_parameters,
                 epp_events_available, epp_minority_available,
                 n_required, epp_required)])
message("EPP: the manuscript and Supplementary Table S2 quote epp_minority_available.")

rubin <- function(est, se, ci = CI_LEVEL) {
  ok <- is.finite(est) & is.finite(se); est <- est[ok]; se <- se[ok]
  m <- length(est); if (m == 0) return(c(NA_real_, NA_real_, NA_real_, NA_real_))
  Q <- mean(est); U <- mean(se^2); B <- if (m >= 2) stats::var(est) else 0
  Tv <- U + (1 + 1 / m) * B; S <- sqrt(max(Tv, 0))
  z <- stats::qnorm(1 - (1 - ci) / 2)
  fmi <- if (Tv > 0) ((1 + 1 / m) * B) / Tv else 0
  c(Q, Q - z * S, Q + z * S, fmi)
}

group_key <- function(lab) {
  lab <- sub("^(.*?):.*$", "\\1", lab)
  if (grepl("_measurable$", lab)) return(lab)
  if (grepl("_c[0-9]+$", lab))    return(sub("_c[0-9]+$", "", lab))
  if (grepl("_c$", lab))          return(sub("_c$", "", lab))
  lab
}

part_key <- function(part, keys) {
  part <- as.character(part)
  keys <- as.character(keys)
  if (!length(keys)) return(group_key(part))
  hit <- keys[startsWith(part, keys)]
  if (length(hit)) hit[which.max(nchar(hit))] else group_key(part)
}
term_keys <- function(term, keys)
  vapply(strsplit(as.character(term), ":", fixed = TRUE)[[1]],
         part_key, character(1), keys = keys)
term_key <- function(term, keys) paste(term_keys(term, keys), collapse = ":")

crude_one <- function(var, imps) {
  est <- se <- numeric(0)
  for (d in imps) {
    if (!var %in% names(d)) next
    f <- stats::as.formula(sprintf("good ~ %s", var))
    fit <- tryCatch(suppressWarnings(glm(f, binomial(), d)), error = function(e) NULL)
    unstable <- is.null(fit) || !isTRUE(fit$converged) ||
      any(!is.finite(coef(fit))) || any(abs(coef(fit)) > 15)
    if (unstable && HAS_LOGISTF) {
      ff <- tryCatch(logistf::logistf(f, data = d, pl = FALSE), error = function(e) NULL)
      if (is.null(ff)) next
      b <- coef(ff)[-1]; s <- sqrt(diag(vcov(ff)))[-1]
    } else {
      if (is.null(fit)) next
      b <- coef(fit)[-1]; s <- sqrt(diag(vcov(fit)))[-1]
    }
    if (length(b) == 0) next
    est <- c(est, b[1]); se <- c(se, s[1])
  }
  r <- rubin(est, se)
  data.table(term = var, variable = relabel_vars(var),
             crude_OR = exp(r[1]), crude_lcl = exp(r[2]), crude_ucl = exp(r[3]),
             fmi = r[4], m_used = length(est))
}

T2_CI_CANDIDATES_LO <- c("conf.low", "lo", "lower", "ci_lo", "CI_low")
T2_CI_CANDIDATES_HI <- c("conf.high", "hi", "upper", "ci_hi", "CI_high")

adj_or_ci <- function(pooled_dt) {
  p <- copy(as.data.table(pooled_dt))
  if (!"estimate" %in% names(p) && "OR" %in% names(p)) p[, estimate := log(OR)]
  if (!"se" %in% names(p)) p[, se := NA_real_]

  p[, `:=`(OR_point   = if ("OR" %in% names(p)) OR else exp(estimate),
           lcl_stored = if ("lcl" %in% names(p)) lcl else NA_real_,
           ucl_stored = if ("ucl" %in% names(p)) ucl else NA_real_)]

  has_t2_ci <- any(c(T2_CI_CANDIDATES_LO) %in% names(p))
  z <- stats::qnorm(1 - (1 - CI_LEVEL) / 2)
  if (has_t2_ci) {
    lo_col <- intersect(T2_CI_CANDIDATES_LO, names(p))[1]
    hi_col <- intersect(T2_CI_CANDIDATES_HI, names(p))[1]
    p[, `:=`(lcl_t2 = exp(as.numeric(get(lo_col))),
             ucl_t2 = exp(as.numeric(get(hi_col))))]
  } else {
    p[, `:=`(lcl_t2 = exp(estimate - z * se),
             ucl_t2 = exp(estimate + z * se))]
  }
  p[]
}

crude_rows <- lapply(PIPELINES, function(nm) {
  fin <- read_safe(file.path(OUT_DIR, sprintf("bnb_final_%s.rds", nm)))
  imp <- read_safe(file.path(OUT_DIR, sprintf("bnb_imp_%s.rds", nm)))
  if (is.null(fin) || is.null(imp)) return(NULL)
  imps <- imp$smcfcs$impDatasets
  gkeys <- names(fin$groups)
  keys <- unique(vapply(fin$final_labels, function(l) term_keys(l, gkeys)[1], character(1)))
  vars <- vapply(keys, function(k)
    if (k %in% names(imps[[1]])) k
    else if (paste0(k, "_c") %in% names(imps[[1]])) paste0(k, "_c")
    else k, character(1))
  ct <- rbindlist(lapply(seq_along(vars), function(i) {
    r <- crude_one(vars[i], imps)
    if (is.null(r) || !nrow(r)) return(NULL)
    r[, key := keys[i]][]
  }), fill = TRUE)
  if (!nrow(ct)) return(NULL)

  adj <- adj_or_ci(fin$pooled)
  if ("term" %in% names(adj)) {
    adj[, term := as.character(term)]
    adj <- adj[term != "(Intercept)"]
    adj[, key := vapply(term, function(t) term_keys(t, gkeys)[1], character(1))]

    lo_use <- if (ADJ_CI_SOURCE == "stored") adj$lcl_stored else adj$lcl_t2
    hi_use <- if (ADJ_CI_SOURCE == "stored") adj$ucl_stored else adj$ucl_t2
    lo_alt <- if (ADJ_CI_SOURCE == "stored") adj$lcl_t2 else adj$lcl_stored
    hi_alt <- if (ADJ_CI_SOURCE == "stored") adj$ucl_t2 else adj$ucl_stored
    lo_use[!is.finite(lo_use)] <- lo_alt[!is.finite(lo_use)]
    hi_use[!is.finite(hi_use)] <- hi_alt[!is.finite(hi_use)]
    adj[, `:=`(lcl_used = lo_use, ucl_used = hi_use)]

    adj_c <- adj[, .(adjusted = paste(sprintf("%s: %s", term,
                                              fmt_or_ci(OR_point, lcl_used, ucl_used)),
                                      collapse = "; ")), by = key]
    ct <- merge(ct, adj_c, by = "key", all.x = TRUE)
    unmatched <- setdiff(adj$key, ct$key)
    if (length(unmatched))
      message(sprintf("  [%s] adjusted terms with no crude counterpart: %s",
                      nm, paste(unmatched, collapse = ", ")))
  }
  attr_qc <- adj[, .(pipeline = nm, label = PIPE_LABELS[[nm]], path = fin$path %||% NA,
                     term, OR = OR_point,
                     lcl_stored, ucl_stored, lcl_t2, ucl_t2)]
  list(tab = cbind(data.table(pipeline = nm, label = PIPE_LABELS[[nm]]), ct), qc = attr_qc)
})

crude_rows <- Filter(Negate(is.null), crude_rows)
tab_crude <- rbindlist(lapply(crude_rows, `[[`, "tab"), fill = TRUE)
saveRDS(tab_crude, file.path(OUT_DIR, "table_s3_crude_vs_adjusted.rds"))
fwrite(tab_crude, file.path(OUT_DIR, "table_s3_crude_vs_adjusted.csv"))
message("---- 15.2 Crude vs adjusted odds ratios ----")
print(tab_crude[, .(pipeline, variable,
                    crude = fmt_or_ci(crude_OR, crude_lcl, crude_ucl),
                    adjusted)])

qc <- rbindlist(lapply(crude_rows, `[[`, "qc"), fill = TRUE)
if (nrow(qc)) {
  qc[, `:=`(
    fmt_stored = fmt_or_ci(OR, lcl_stored, ucl_stored),
    fmt_table2 = fmt_or_ci(OR, lcl_t2,     ucl_t2))]
  qc[, comparable := is.finite(lcl_stored) & is.finite(ucl_stored) &
                     is.finite(lcl_t2)     & is.finite(ucl_t2)]
  qc[, differs_after_rounding := comparable & (fmt_stored != fmt_table2)]
  fwrite(qc, file.path(OUT_DIR, "qc_s3_vs_table2_ci.csv"))
  nd <- sum(qc$differs_after_rounding, na.rm = TRUE)
  message(sprintf("---- A-7 (a) interval audit: %d of %d terms differ at %d decimal places ----",
                  nd, nrow(qc), DIGITS_OR))
  if (nd > 0) {
    print(qc[differs_after_rounding == TRUE,
             .(pipeline, term, stored = fmt_stored, table2 = fmt_table2)])
    message("`stored`  = the interval Section 10 saved (Rubin t for glm; for Firth, a CLIP")
    message("            profile interval only when method_detail$clip_ok is TRUE).")
    message("`table2`  = the normal Wald interval the ORIGINAL script 16 rebuilt, because its")
    message("            pick_col() candidates omitted 'lcl' / 'ucl'.")
    message(sprintf("S3 was written using ADJ_CI_SOURCE = '%s'.", ADJ_CI_SOURCE))
    message("These are the cells that move between the two conventions; ",
            "16_table2_model_specification-1.R writes the same list to ",
            "data/table2_ci_change_log.csv.")
  }
}

prov_path <- file.path(OUT_DIR, "table2_ci_provenance.csv")
if (file.exists(prov_path)) {
  prov16 <- data.table::fread(prov_path)
  setting16 <- unique(as.character(prov16$ci_source_setting))
  expected  <- if (ADJ_CI_SOURCE == "stored") "pooled" else "wald"
  message("---- Table 2 / S8 on disk was produced with CI_SOURCE = '",
          paste(setting16, collapse = "/"), "' ----")
  if (!all(setting16 == expected))
    warning(sprintf(paste0("S3 is being written with ADJ_CI_SOURCE = '%s', which pairs with ",
                           "CI_SOURCE = '%s', but data/table2_ci_provenance.csv says script 16 ",
                           "was last run with '%s'. S3 and Table 2 will disagree. Re-run one of ",
                           "the two."),
                   ADJ_CI_SOURCE, expected, paste(setting16, collapse = "/")))
  fir <- prov16[estimator == "firth"]
  if (nrow(fir))
    for (i in seq_len(nrow(fir)))
      message(sprintf("  [%s] Firth, clip_ok = %s -> %s",
                      fir$pipeline[i], fir$clip_ok[i],
                      if (isTRUE(as.logical(fir$clip_ok[i])))
                        "combination-of-likelihood-profile interval"
                      else paste0("Rubin-pooled NORMAL approximation, NOT a profile interval; ",
                                  "correct S8's footnote")))
} else {
  message("---- data/table2_ci_provenance.csv not found ----")
  message("Run 16_table2_model_specification-1.R first, so that S3 and Table 2 are known to ",
          "use the same interval convention.")
}

contrast <- function(d, flag, lab, stratum = "All") {
  a <- d[get(flag) == TRUE]; b <- d[get(flag) == FALSE]
  n1 <- nrow(a); n0 <- nrow(b)
  if (n1 == 0 || n0 == 0) return(NULL)
  e1 <- sum(a$good, na.rm = TRUE); e0 <- sum(b$good, na.rm = TRUE)
  p1 <- e1 / n1; p0 <- e0 / n0
  or  <- (p1 / (1 - p1)) / (p0 / (1 - p0))
  rr  <- p1 / p0
  rd  <- p1 - p0
  se_lrr  <- sqrt(1 / max(e1, .5) - 1 / n1 + 1 / max(e0, .5) - 1 / n0)
  se_rd   <- sqrt(p1 * (1 - p1) / n1 + p0 * (1 - p0) / n0)
  se_lor  <- sqrt(1 / max(e1, .5) + 1 / max(n1 - e1, .5) +
                  1 / max(e0, .5) + 1 / max(n0 - e0, .5))
  z <- stats::qnorm(1 - (1 - CI_LEVEL) / 2)
  data.table(subset = lab, stratum = stratum,
             n_subset = n1, events_subset = e1, risk_subset = round(p1, 3),
             n_reference = n0, events_reference = e0, risk_reference = round(p0, 3),
             OR = round(or, 3),
             OR_lcl = round(exp(log(or) - z * se_lor), 3),
             OR_ucl = round(exp(log(or) + z * se_lor), 3),
             RR = round(rr, 3), RR_lcl = round(exp(log(rr) - z * se_lrr), 3),
             RR_ucl = round(exp(log(rr) + z * se_lrr), 3),
             RD = round(rd, 3), RD_lcl = round(rd - z * se_rd, 3),
             RD_ucl = round(rd + z * se_rd, 3),
             OR_over_RR = round(or / rr, 2))
}

dat[, sev_and_eld := severe & elderly]
abs_rows <- list(contrast(dat, "severe", "Severe (admission motor FIM <=26)"),
                 contrast(dat, "elderly", "Elderly (age >=90)"),
                 contrast(dat, "sev_and_eld", "Severe and elderly"))
if ("class" %in% names(dat)) {
  for (cl in levels(factor(dat$class))) {
    in_class <- as.character(dat$class) == cl
    if (!sum(in_class, na.rm = TRUE)) next

    dat[, "__f_sev" := dat$severe  & in_class]
    dat[, "__f_eld" := dat$elderly & in_class]
    abs_rows <- c(abs_rows,
                  list(contrast(dat, "__f_sev", "Severe (admission motor FIM <=26)", cl),
                       contrast(dat, "__f_eld", "Elderly (age >=90)", cl)))
    dat[, c("__f_sev", "__f_eld") := NULL]

    sub <- dat[class == cl]
    abs_rows <- c(abs_rows,
                  list(contrast(sub, "severe",  "Severe (admission motor FIM <=26)",
                                paste0(cl, " (within-disease reference)")),
                       contrast(sub, "elderly", "Elderly (age >=90)",
                                paste0(cl, " (within-disease reference)"))))
  }
}
tab_abs <- rbindlist(Filter(Negate(is.null), abs_rows), fill = TRUE)
fwrite(tab_abs, file.path(OUT_DIR, "table_s4_absolute_risk.csv"))
message("---- 15.3 Absolute risks, risk ratios and risk differences ----")
print(tab_abs[, .(subset, stratum, risk_subset, risk_reference, OR, RR, RD, OR_over_RR)])

int_rows <- lapply(PIPELINES, function(nm) {
  fin <- read_safe(file.path(OUT_DIR, sprintf("bnb_final_%s.rds", nm)))
  ix  <- read_safe(file.path(OUT_DIR, sprintf("bnb_interactions_%s.rds", nm)))
  if (is.null(fin) || !length(fin$selected_int_keys)) return(NULL)
  pooled <- as.data.table(fin$pooled)
  if ("term" %in% names(pooled)) pooled[, term := as.character(term)]
  gkeys  <- names(fin$groups)
  cd <- if (!is.null(ix) && !is.null(ix$candidates)) as.data.table(ix$candidates) else NULL
  if (!is.null(cd) && !"key" %in% names(cd) && all(c("v1", "v2") %in% names(cd)))
    cd[, key := paste0(v1, ":", v2)]
  cd_get <- function(k, col) {
    if (is.null(cd) || !"key" %in% names(cd) || !col %in% names(cd)) return(NA)
    v <- cd[key == k][[col]]
    if (length(v)) v[1] else NA
  }
  rbindlist(lapply(fin$selected_int_keys, function(k) {
    est <- pooled[vapply(term, function(t) identical(term_key(t, gkeys), k), logical(1))]
    if (!nrow(est)) {
      est <- pooled[term %in% fin$selected_int_groups[[k]]]
      if (!nrow(est)) {
        message(sprintf(paste0("  [%s] interaction '%s' was selected in Section 10 ",
                               "but no matching coefficient is in the pooled table."),
                        nm, k))
        return(NULL)
      }
    }
    data.table(pipeline = nm, label = PIPE_LABELS[[nm]], interaction = k,
               interaction_disp = paste(relabel_vars(strsplit(k, ":", fixed = TRUE)[[1]]),
                                        collapse = " x "),
               term = est$term, OR = round(est$OR, 3),
               lcl = round(est$lcl, 3), ucl = round(est$ucl, 3),
               LRT_p = signif(as.numeric(cd_get(k, "p")), 3),
               LRT_df = as.integer(cd_get(k, "df")))
  }), fill = TRUE)
})
tab_int <- rbindlist(Filter(Negate(is.null), int_rows), fill = TRUE)
fwrite(tab_int, file.path(OUT_DIR, "table_s5_interaction_estimates.csv"))
message("---- 15.4 Interaction estimates ----"); print(tab_int)

FOCUS_INT <- list(c("age", "support_in"))

fit_pooled_one <- function(smf, imps, want, path = "glm") {
  est <- se <- numeric(0)
  for (d in imps) {
    b <- s <- NULL
    use_firth <- identical(path, "firth")
    fit <- NULL
    if (!use_firth) {
      fit <- tryCatch(suppressWarnings(glm(smf, binomial(), d)), error = function(e) NULL)
      use_firth <- is.null(fit) || !isTRUE(fit$converged) ||
        any(!is.finite(coef(fit))) || any(abs(coef(fit)) > 15)
    }
    if (use_firth && HAS_LOGISTF) {
      ff <- tryCatch(logistf::logistf(smf, data = d, pl = FALSE), error = function(e) NULL)
      if (is.null(ff)) next
      b <- coef(ff); s <- suppressWarnings(sqrt(diag(vcov(ff))))
      if (is.null(names(b))) names(b) <- rownames(vcov(ff))
    } else if (!is.null(fit)) {
      b <- coef(fit); s <- suppressWarnings(sqrt(diag(vcov(fit))))
    }
    if (is.null(b) || !want %in% names(b)) next
    est <- c(est, unname(b[want])); se <- c(se, unname(s[which(names(b) == want)]))
  }
  r <- rubin(est, se)
  data.table(term = want, OR = exp(r[1]), lcl = exp(r[2]), ucl = exp(r[3]),
             fmi = r[4], m_used = length(est))
}

focus_one_sample <- function(nm, fin, imps, sample_id, sample_note) {
  keys_here <- vapply(FOCUS_INT, paste, character(1), collapse = ":")
  keys_here <- intersect(keys_here, fin$selected_int_keys)
  if (!length(keys_here)) return(NULL)
  smf <- stats::as.formula(fin$smformula_final)
  rbindlist(lapply(keys_here, function(k) {
    pr <- strsplit(k, ":", fixed = TRUE)[[1]]
    cont_key <- pr[1]; fac_key <- pr[2]
    cont_labs <- fin$groups[[cont_key]]
    if (length(cont_labs) != 1) {
      message(sprintf(paste0("  [%s] %s enters the model as %d spline basis terms, ",
                             "so a single simple slope is not defined; see the ",
                             "odds-ratio curve of Section 11 instead."),
                      nm, cont_key, length(cont_labs)))
      return(NULL)
    }
    if (!fac_key %in% names(imps[[1]]) ||
        !is.factor(as.data.frame(imps[[1]])[[fac_key]])) return(NULL)
    lvs <- levels(as.data.frame(imps[[1]])[[fac_key]])
    rbindlist(lapply(lvs, function(lv) {
      di <- lapply(imps, function(d) {
        d <- as.data.frame(d)
        d[[fac_key]] <- stats::relevel(factor(d[[fac_key]]), ref = lv)
        d
      })
      r <- fit_pooled_one(smf, di, cont_labs, fin$path)
      data.table(pipeline = nm, label = PIPE_LABELS[[nm]],
                 sample = sample_id, sample_note = sample_note,
                 interaction = k,
                 stratum_variable = unname(relabel_vars(fac_key)), stratum = lv,
                 n_stratum = sum(as.character(as.data.frame(imps[[1]])[[fac_key]]) == lv),
                 n_total = nrow(as.data.frame(imps[[1]])),
                 slope_of = unname(relabel_vars(cont_key)), term = r$term,
                 OR = round(r$OR, 3), lcl = round(r$lcl, 3), ucl = round(r$ucl, 3),
                 fmi = round(r$fmi, 3), m_used = r$m_used)
    }), fill = TRUE)
  }), fill = TRUE)
}

focus_rows <- lapply(PIPELINES, function(nm) {
  fin <- read_safe(file.path(OUT_DIR, sprintf("bnb_final_%s.rds", nm)))
  imp <- read_safe(file.path(OUT_DIR, sprintf("bnb_imp_%s.rds", nm)))
  if (is.null(fin) || is.null(imp)) return(NULL)

  imps_all <- imp$smcfcs$impDatasets
  imps_red <- reduce_imps(imp, fin$n0, fin$outlier_ids)
  n_drop   <- fin$n0 - nrow(as.data.frame(imps_red[[1]]))

  rbindlist(list(
    focus_one_sample(nm, fin, imps_red, "analysed",
                     sprintf(paste0("Sections 10-12 analysis sample: %d of %d cases, ",
                                    "%d influential observations excluded. Same sample ",
                                    "as Table 2."),
                             fin$n0 - n_drop, fin$n0, n_drop)),
    focus_one_sample(nm, fin, imps_all, "all_imputed",
                     sprintf(paste0("All imputed cases: %d. Influential observations ",
                                    "NOT excluded; differs from Table 2."), fin$n0))
  ), fill = TRUE)
})
tab_focus_both <- rbindlist(Filter(Negate(is.null), focus_rows), fill = TRUE)

if (nrow(tab_focus_both)) {
  fwrite(tab_focus_both, file.path(OUT_DIR, "table_s5b_focus_interaction_bysample.csv"))

  primary_id <- if (S5B_SAMPLE == "analysed") "analysed" else "all_imputed"
  tab_focus  <- tab_focus_both[sample == primary_id]
  fwrite(tab_focus, file.path(OUT_DIR, "table_s5b_focus_interaction.csv"))

  message("---- 15.4b Age slope within each care group (severe subset) ----")
  message(sprintf("Primary sample written to table_s5b_focus_interaction.csv: '%s'.", primary_id))
  print(tab_focus[, .(pipeline, stratum_variable, stratum, n_stratum, n_total, slope_of,
                      OR = sprintf("%.3f (%.3f to %.3f)", OR, lcl, ucl), m_used)])

  message("\n---- 15.4b both samples, for the S5b footnote ----")
  print(tab_focus_both[, .(pipeline, sample, stratum, n_stratum, n_total,
                           OR = sprintf("%.3f (%.3f to %.3f)", OR, lcl, ucl))])
  chk <- tab_focus_both[, .(strata_sum = sum(n_stratum), n_total = n_total[1]),
                        by = .(pipeline, sample)]
  message("Reconciliation of the stratum counts with the analysed n:")
  print(chk)
  message("The strata of the primary sample must sum to the n quoted in the Methods ",
          "and in Table 2; if they do not, S5b and Table 2 are fitted on different data.")
  message("OR is per one-unit increase in the continuous variable, inside the ",
          "stratum named; the two rows come from one model refitted with the ",
          "reference level switched, so they share the interaction term.")
} else {
  tab_focus <- data.table()
  message("---- 15.4b no focus interaction present in any final model ----")
  message("FOCUS_INT is set to: ",
          paste(vapply(FOCUS_INT, paste, character(1), collapse = " x "), collapse = "; "),
          ". Selected interactions are: ",
          paste(unique(tab_int$interaction), collapse = "; "))
}

message("\nSection 15 complete (revision 2026-08-16).")
message("  data/table_s2_sample_size.csv                  (S2; both EPP denominators)")
message("  data/table_s3_crude_vs_adjusted.csv, .rds      (S3)")
message("  data/table_s4_absolute_risk.csv                (S4; also feeds Figure 3)")
message("  data/table_s5_interaction_estimates.csv        (S5)")
if (nrow(tab_focus)) {
  message("  data/table_s5b_focus_interaction.csv           (S5b, primary sample)")
  message("  data/table_s5b_focus_interaction_bysample.csv  (S5b, audit trail)")
}
message("  data/qc_s3_vs_table2_ci.csv                    (A-7a interval audit)")
message("\nSettings used: S5B_SAMPLE = '", S5B_SAMPLE,
        "', ADJ_CI_SOURCE = '", ADJ_CI_SOURCE, "'.")
