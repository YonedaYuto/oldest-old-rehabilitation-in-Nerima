library(data.table)
library(here)

OUT_DIR <- here::here("data")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

if (!exists("relabel_vars")) {
  LABELS_CANDIDATES <- c(
    here::here("00_labels.R"),
    here::here("scripts", "00_labels.R"),
    here::here("00_labels_1.R"),
    here::here("scripts", "00_labels_1.R")
  )
  LABELS_PATH <- LABELS_CANDIDATES[file.exists(LABELS_CANDIDATES)][1]
  if (is.na(LABELS_PATH))
    stop("00_labels.R not found. Looked in:\n  ",
         paste(LABELS_CANDIDATES, collapse = "\n  "))
  source(LABELS_PATH)
}

HAS_LOGISTF <- requireNamespace("logistf", quietly = TRUE)
CI_LEVEL    <- 0.95
PIPELINES   <- names(PIPE_LABELS)
read_safe   <- function(f) if (file.exists(f)) readRDS(f) else NULL
`%||%`      <- function(a, b) if (is.null(a)) b else a
set.seed(20240601)

ADJ_CI_SOURCE <- "stored"
DIGITS_OR     <- 2L

OUT_SUFFIX <- "_multilevel"

if (!ADJ_CI_SOURCE %in% c("table2_compatible", "stored"))
  stop("ADJ_CI_SOURCE must be 'table2_compatible' or 'stored'.")

out_path <- function(stem, ext) file.path(OUT_DIR, paste0(stem, OUT_SUFFIX, ".", ext))

fmt_or_ci <- function(est, lo, hi, digits = DIGITS_OR)
  sprintf("%.*f (%.*f to %.*f)", digits, est, digits, lo, digits, hi)

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

crude_one_multi <- function(var, imps) {
  est_by <- list(); se_by <- list()
  order_seen <- character(0)
  n_firth <- 0L; n_fit <- 0L
  ref_level <- NA_character_; var_levels <- character(0)

  for (d in imps) {
    if (!var %in% names(d)) next
    if (is.na(ref_level) && (is.factor(d[[var]]) || is.character(d[[var]]))) {
      lv <- if (is.factor(d[[var]])) levels(d[[var]]) else sort(unique(as.character(d[[var]])))
      var_levels <- lv
      ref_level  <- lv[1]
    }
    f   <- stats::as.formula(sprintf("good ~ %s", var))
    fit <- tryCatch(suppressWarnings(glm(f, binomial(), d)), error = function(e) NULL)
    unstable <- is.null(fit) || !isTRUE(fit$converged) ||
      any(!is.finite(coef(fit))) || any(abs(coef(fit)) > 15)
    if (unstable && HAS_LOGISTF) {
      ff <- tryCatch(logistf::logistf(f, data = d, pl = FALSE), error = function(e) NULL)
      if (is.null(ff)) next
      b <- coef(ff)[-1]; s <- sqrt(diag(vcov(ff)))[-1]
      n_firth <- n_firth + 1L
    } else {
      if (is.null(fit)) next
      b <- coef(fit)[-1]; s <- sqrt(diag(vcov(fit)))[-1]
    }
    if (length(b) == 0) next
    nmb <- names(b)
    if (is.null(nmb) || any(!nzchar(nmb)))
      nmb <- if (length(b) == 1L) var else paste0(var, seq_along(b))
    for (j in seq_along(b)) {
      k <- nmb[j]
      if (!k %in% order_seen) order_seen <- c(order_seen, k)
      est_by[[k]] <- c(est_by[[k]], unname(b[j]))
      se_by[[k]]  <- c(se_by[[k]],  unname(s[j]))
    }
    n_fit <- n_fit + 1L
  }
  if (!length(order_seen)) return(NULL)

  rbindlist(lapply(order_seen, function(k) {
    r <- rubin(est_by[[k]], se_by[[k]])
    lev <- if (length(var_levels) && startsWith(k, var)) {
      cand <- substring(k, nchar(var) + 1L)
      if (nzchar(cand) && cand %in% var_levels)
        cand else NA_character_
    } else NA_character_
    term_label <- if (!is.na(lev)) paste0(var, ":", lev) else k
    contrast <- if (!is.na(lev) && !is.na(ref_level))
      sprintf("%s vs %s", relabel_levels(var, lev), relabel_levels(var, ref_level))
    else NA_character_
    data.table(
      term       = k,
      term_label = term_label,
      level      = lev,
      ref_level  = ref_level,
      contrast   = contrast,
      variable   = relabel_vars(var),
      term_variable = relabel_vars(term_label),
      crude_OR = exp(r[1]), crude_lcl = exp(r[2]), crude_ucl = exp(r[3]),
      fmi = r[4], m_used = length(est_by[[k]]),
      n_imp_fitted = n_fit, n_imp_firth = n_firth)
  }), fill = TRUE)
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
  imp <- read_safe(file.path(OUT_DIR, sprintf("bnb_imp_%s.rds",   nm)))
  if (is.null(fin) || is.null(imp)) return(NULL)
  imps  <- imp$smcfcs$impDatasets
  gkeys <- names(fin$groups)
  keys  <- unique(vapply(fin$final_labels, function(l) term_keys(l, gkeys)[1], character(1)))
  vars  <- vapply(keys, function(k)
    if (k %in% names(imps[[1]])) k
    else if (paste0(k, "_c") %in% names(imps[[1]])) paste0(k, "_c")
    else k, character(1))
  ct <- rbindlist(lapply(seq_along(vars), function(i) {
    r <- crude_one_multi(unname(vars[i]), imps)
    if (is.null(r) || !nrow(r)) return(NULL)
    r[, `:=`(key = unname(keys[i]), fitted_on = unname(vars[i]))][]
  }), fill = TRUE)
  if (!nrow(ct)) return(NULL)

  adj_c <- NULL
  adj   <- adj_or_ci(fin$pooled)
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
    adj[, adj_fmt := fmt_or_ci(OR_point, lcl_used, ucl_used)]

    adj_c <- adj[, .(adjusted = paste(sprintf("%s: %s", term, adj_fmt), collapse = "; "),
                     n_adj_terms = .N), by = key]
    ct <- merge(ct, adj_c, by = "key", all.x = TRUE)

    ct <- merge(ct, adj[, .(key, term, adjusted_term = adj_fmt)],
                by = c("key", "term"), all.x = TRUE, sort = FALSE)

    unmatched <- setdiff(adj$key, ct$key)
    if (length(unmatched))
      message(sprintf("  [%s] adjusted terms with no crude counterpart: %s",
                      nm, paste(unmatched, collapse = ", ")))
  } else {
    ct[, `:=`(adjusted = NA_character_, n_adj_terms = NA_integer_,
              adjusted_term = NA_character_)]
  }

  ct[, crude := fmt_or_ci(crude_OR, crude_lcl, crude_ucl)]
  front <- intersect(c("key", "variable", "term", "term_label", "level",
                       "ref_level", "contrast", "term_variable",
                       "fitted_on", "crude_OR", "crude_lcl", "crude_ucl",
                       "crude", "fmi", "m_used", "n_imp_fitted",
                       "n_imp_firth", "adjusted_term", "n_adj_terms",
                       "adjusted"), names(ct))
  setcolorder(ct, c(front, setdiff(names(ct), front)))

  attr_qc <- if ("term" %in% names(adj)) {
    adj[, .(pipeline = nm, label = PIPE_LABELS[[nm]], path = fin$path %||% NA,
            term, OR = OR_point, lcl_stored, ucl_stored, lcl_t2, ucl_t2)]
  } else data.table()

  n_adj_by_key <- if (!is.null(adj_c)) {
    adj_c[, .(key, n_adj_terms)]
  } else {
    as.data.table(data.frame(key = character(0), n_adj_terms = integer(0),
                             stringsAsFactors = FALSE))
  }
  cov_qc <- merge(ct[, .(n_crude_terms = .N), by = key], n_adj_by_key,
                  by = "key", all.x = TRUE)
  cov_qc[, `:=`(pipeline = nm, label = PIPE_LABELS[[nm]],
                variable = relabel_vars(key))]
  cov_qc[, crude_terms_lt_adjusted :=
           is.finite(n_adj_terms) & n_crude_terms < n_adj_terms]

  list(tab = cbind(data.table(pipeline = nm, label = PIPE_LABELS[[nm]]), ct),
       qc  = attr_qc, cov = cov_qc)
})

crude_rows <- Filter(Negate(is.null), crude_rows)
tab_crude  <- rbindlist(lapply(crude_rows, `[[`, "tab"), fill = TRUE)
saveRDS(tab_crude, out_path("table_s3_crude_vs_adjusted", "rds"))
fwrite(  tab_crude, out_path("table_s3_crude_vs_adjusted", "csv"))

message("---- 15.2 Crude vs adjusted odds ratios (all factor contrasts) ----")
print(tab_crude[, .(pipeline,
                    variable = ifelse(is.na(contrast), variable,
                                      paste0(variable, " [", contrast, "]")),
                    crude, adjusted_term)])

covq <- rbindlist(lapply(crude_rows, `[[`, "cov"), fill = TRUE)
if (nrow(covq)) {
  setcolorder(covq, c("pipeline", "label", "key", "variable",
                      "n_crude_terms", "n_adj_terms", "crude_terms_lt_adjusted"))
  fwrite(covq, out_path("qc_s3_crude_terms", "csv"))
  short <- covq[crude_terms_lt_adjusted == TRUE]
  message(sprintf("---- coverage audit: %d of %d variable keys have fewer crude than adjusted terms ----",
                  nrow(short), nrow(covq)))
  if (nrow(short)) {
    print(short[, .(pipeline, variable, n_crude_terms, n_adj_terms)])
    message("A key can legitimately show fewer crude than adjusted terms when the ",
            "adjusted side expands into an RCS basis or an interaction; a k-level ",
            "FACTOR showing 1 crude term against k-1 adjusted terms is the defect ",
            "this script repairs and should now be absent.")
  } else {
    message("Every variable key now contributes at least as many crude terms as ",
            "adjusted terms of the same key.")
  }
  multi <- covq[n_crude_terms > 1]
  if (nrow(multi)) {
    message("---- multi-contrast predictors now reported in full ----")
    print(multi[, .(pipeline, variable, n_crude_terms)])
  }
}

qc <- rbindlist(lapply(crude_rows, `[[`, "qc"), fill = TRUE)
if (nrow(qc)) {
  qc[, `:=`(
    fmt_stored = fmt_or_ci(OR, lcl_stored, ucl_stored),
    fmt_table2 = fmt_or_ci(OR, lcl_t2,     ucl_t2))]
  qc[, comparable := is.finite(lcl_stored) & is.finite(ucl_stored) &
                     is.finite(lcl_t2)     & is.finite(ucl_t2)]
  qc[, differs_after_rounding := comparable & (fmt_stored != fmt_table2)]
  fwrite(qc, out_path("qc_s3_vs_table2_ci", "csv"))
  nd <- sum(qc$differs_after_rounding, na.rm = TRUE)
  message(sprintf("---- A-7 (a) interval audit: %d of %d terms differ at %d decimal places ----",
                  nd, nrow(qc), DIGITS_OR))
  if (nd > 0)
    print(qc[differs_after_rounding == TRUE,
             .(pipeline, term, stored = fmt_stored, table2 = fmt_table2)])
}

prov_path <- file.path(OUT_DIR, "table2_ci_provenance.csv")
if (file.exists(prov_path)) {
  prov16    <- data.table::fread(prov_path)
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
}

message("\nSection 15-1 complete (corrected Supplementary Table S3 only).")
message("  ", out_path("table_s3_crude_vs_adjusted", "csv"), "   (S3, corrected)")
message("  ", out_path("table_s3_crude_vs_adjusted", "rds"))
if (nrow(covq))
  message("  ", out_path("qc_s3_crude_terms",  "csv"), "   (crude vs adjusted term counts)")
if (nrow(qc))
  message("  ", out_path("qc_s3_vs_table2_ci", "csv"), "   (A-7a interval audit)")
message("Settings used: ADJ_CI_SOURCE = '", ADJ_CI_SOURCE,
        "', DIGITS_OR = ", DIGITS_OR, ", OUT_SUFFIX = '", OUT_SUFFIX, "'.")
message("15.1 (S2), 15.3 (S4) and 15.4 (S5 / S5b) were NOT re-run; ",
        "their files from 15_reporting_supplements-1.R are unchanged.")
