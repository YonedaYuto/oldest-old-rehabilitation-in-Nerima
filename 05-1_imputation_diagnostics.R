library(data.table)
library(here)
library(ggplot2)

need_pkg <- function(p) if (!requireNamespace(p, quietly = TRUE))
  stop(sprintf("パッケージ '%s' が必要です。install.packages('%s')", p, p))

LABELS_PATH <- here::here("00_labels.R")
if (!exists("relabel_vars")) {
  if (!file.exists(LABELS_PATH)) stop("00_labels.R が見つかりません: ", LABELS_PATH)
  source(LABELS_PATH)
}

OUT_DIR <- here::here("data")
FIG_DIR <- here::here("figures")
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

RERUN_TRACE_CHAINS <- 0L
RERUN_RJLIMIT      <- 5000L
MC_REL_SE_TOL      <- 0.05
VALID_RANGE <- list(MMSE_in = c(0, 30), BBS_in = c(0, 56))

imp_index <- readRDS(file.path(OUT_DIR, "bnb_imp_index.rds"))

pool_rubin <- function(fitlist, n_complete) {
  m      <- length(fitlist)
  betas  <- lapply(fitlist, coef)
  terms  <- names(betas[[1]])
  Vdiag  <- lapply(fitlist, function(f) diag(as.matrix(vcov(f))))

  getv <- function(lst, t) vapply(lst, function(z) unname(z[t]), numeric(1))
  Qbar <- vapply(terms, function(t) mean(getv(betas, t)), numeric(1))
  B    <- vapply(terms, function(t) stats::var(getv(betas, t)), numeric(1))
  Ubar <- vapply(terms, function(t) mean(getv(Vdiag, t)), numeric(1))

  Tvar   <- Ubar + (1 + 1/m) * B
  riv    <- (1 + 1/m) * B / Ubar
  lambda <- (1 + 1/m) * B / Tvar

  k          <- length(terms)
  dfc        <- max(n_complete - k, 1)
  df_old     <- ifelse(lambda > 0, (m - 1) / lambda^2, Inf)
  df_obs     <- ((dfc + 1) / (dfc + 3)) * dfc * (1 - lambda)
  df         <- df_old * df_obs / (df_old + df_obs)
  fmi        <- (riv + 2 / (df + 3)) / (riv + 1)

  se        <- sqrt(Tvar)
  mc_rel_se <- sqrt(B / m) / se
  m_needed  <- ceiling(100 * fmi)

  data.table(term = terms, estimate = Qbar, se = se,
             riv = riv, lambda = lambda, df = df,
             fmi = fmi, mc_rel_se = mc_rel_se, m_needed = m_needed)
}

trace_dt <- function(sci, coef_names) {
  d <- dim(sci)
  if (is.null(d) || length(d) != 3) return(NULL)
  nch <- d[1]; nco <- d[2]; nit <- d[3]
  if (length(coef_names) != nco) coef_names <- paste0("b", seq_len(nco))
  out <- vector("list", nch * nco)
  z <- 1L
  for (c0 in seq_len(nch)) for (j in seq_len(nco)) {
    out[[z]] <- data.table(chain = c0, coef = coef_names[j],
                           iter = seq_len(nit), est = as.numeric(sci[c0, j, ]))
    z <- z + 1L
  }
  rbindlist(out)
}

for (nm in imp_index) {
  obj  <- readRDS(file.path(OUT_DIR, sprintf("bnb_imp_%s.rds", nm)))
  imps <- obj$smcfcs$impDatasets
  m    <- obj$m
  smf  <- stats::as.formula(obj$smformula)
  long <- as.data.table(obj$long)
  orig <- long[.imp == 0]
  n    <- obj$n
  imputed_vars <- obj$imputed

  message(sprintf("================ [%s] 代入診断 (n=%d, m=%d) ================", nm, n, m))

  fit1 <- glm(smf, family = binomial(), data = imps[[1]])
  cf_names <- names(coef(fit1))

  sci <- obj$smcfcs$smCoefIter
  tdt <- trace_dt(sci, cf_names)

  if (RERUN_TRACE_CHAINS > 0L) {
    need_pkg("smcfcs"); need_pkg("Hmisc")
    od <- as.data.frame(orig[, setdiff(names(orig), c(".imp", ".id")), with = FALSE])
    od <- od[, names(obj$meth), drop = FALSE]
    for (cl in names(od)) if (is.factor(imps[[1]][[cl]]))
      od[[cl]] <- factor(od[[cl]], levels = levels(imps[[1]][[cl]]))
    nit <- if (!is.null(dim(sci))) dim(sci)[3] else 100L
    set.seed(20240L)
    re <- lapply(seq_len(RERUN_TRACE_CHAINS), function(i)
      smcfcs::smcfcs(od, smtype = "logistic", smformula = obj$smformula,
                     method = obj$meth, predictorMatrix = obj$predmat,
                     m = 1, numit = nit, rjlimit = RERUN_RJLIMIT))
    tdt <- rbindlist(lapply(seq_along(re), function(i) {
      t0 <- trace_dt(re[[i]]$smCoefIter, cf_names); t0[, chain := i]; t0
    }))
  }

  if (!is.null(tdt)) {
    p_conv <- ggplot(tdt, aes(iter, est, group = factor(chain), colour = factor(chain))) +
      geom_line(linewidth = 0.4, alpha = 0.8) +
      facet_wrap(~ coef, scales = "free_y") +
      labs(title = sprintf("%s  -  smcfcs convergence trace", relabel_pipeline(nm)),
           subtitle = "Estimates should stabilise (no trend) by the final iteration",
           x = "Iteration", y = "Coefficient estimate", colour = "Chain") +
      theme_bw(base_size = 9) +
      theme(legend.position = if (max(tdt$chain) > 1) "top" else "none")
    nco <- length(unique(tdt$coef)); ncol_f <- min(4L, nco)
    ggsave(file.path(FIG_DIR, sprintf("fig_imp_convergence_%s.png", nm)),
           p_conv, width = 3.2 * ncol_f, height = 2.6 * ceiling(nco / ncol_f) + 0.6, dpi = 150)
  } else {
    message("  収束トレース: smCoefIter が取得できませんでした。")
  }

  sep_flag <- FALSE
  fits <- lapply(imps, function(d) {
    withCallingHandlers(
      glm(smf, family = binomial(), data = d),
      warning = function(w) {
        if (grepl("fitted probabilities", conditionMessage(w))) sep_flag <<- TRUE
        invokeRestart("muffleWarning")
      })
  })
  pool <- pool_rubin(fits, n_complete = n)
  pool_print <- copy(pool)[, `:=`(estimate = round(estimate, 3), se = round(se, 3),
                                  riv = round(riv, 3), lambda = round(lambda, 3),
                                  df = round(df, 1), fmi = round(fmi, 3),
                                  mc_rel_se = round(mc_rel_se, 3))]
  print(pool_print)

  max_fmi   <- max(pool$fmi, na.rm = TRUE)
  m_needed  <- max(pool$m_needed, na.rm = TRUE)
  max_mcse  <- max(pool$mc_rel_se, na.rm = TRUE)
  m_ok      <- (m >= m_needed) && (max_mcse <= MC_REL_SE_TOL)

  message(sprintf("  最大FMI=%.3f → 目安の必要m≈%d（現在m=%d）", max_fmi, m_needed, m))
  message(sprintf("  最大の相対MC誤差=%.3f（許容%.2f）→ m は %s",
                  max_mcse, MC_REL_SE_TOL, ifelse(m_ok, "十分", "増やすことを推奨")))
  if (sep_flag)
    message("  ※完全分離の兆候（fitted prob 0/1）を検出。FMI/SE は不安定な可能性。§6 で Firth 等に分岐。")

  range_tab <- data.table()
  dens_dt   <- data.table()
  cont_imp  <- intersect(imputed_vars, names(VALID_RANGE))

  for (v in cont_imp) {
    miss_id  <- orig[is.na(get(v)), .id]
    obs_val  <- orig[!is.na(get(v)), as.numeric(get(v))]
    imp_val  <- long[.imp >= 1 & .id %in% miss_id, as.numeric(get(v))]
    if (length(imp_val) == 0) next
    dens_dt <- rbind(dens_dt,
      data.table(variable = relabel_vars(v), source = "Observed", value = obs_val),
      data.table(variable = relabel_vars(v), source = "Imputed",  value = imp_val))
    lo <- VALID_RANGE[[v]][1]; hi <- VALID_RANGE[[v]][2]
    range_tab <- rbind(range_tab, data.table(
      variable = v, n_missing = length(miss_id), n_imp_values = length(imp_val),
      below = sum(imp_val < lo), above = sum(imp_val > hi),
      pct_out = round(100 * mean(imp_val < lo | imp_val > hi), 2)))
  }

  if (nrow(dens_dt) > 0) {
    p_dist <- ggplot(dens_dt, aes(value, colour = source, fill = source)) +
      geom_density(alpha = 0.25, linewidth = 0.5) +
      facet_wrap(~ variable, scales = "free") +
      scale_colour_manual(values = c(Observed = "#2166AC", Imputed = "#B2182B")) +
      scale_fill_manual(values   = c(Observed = "#2166AC", Imputed = "#B2182B")) +
      labs(title = sprintf("%s  -  Observed vs imputed (continuous)", relabel_pipeline(nm)),
           x = NULL, y = "Density", colour = NULL, fill = NULL) +
      theme_bw(base_size = 10) + theme(legend.position = "top")
    nfac <- length(unique(dens_dt$variable))
    ggsave(file.path(FIG_DIR, sprintf("fig_imp_distributions_%s.png", nm)),
           p_dist, width = 4.2 * min(2L, nfac), height = 3.4, dpi = 150)
  }
  if (nrow(range_tab) > 0) {
    message("  連続代入値の範囲外チェック:"); print(range_tab)
    if (any(range_tab$pct_out > 0))
      message("  ※範囲外の代入値あり（norm 法の特性）。境界切り詰めや別法の検討余地。")
  }

  if ("JCS_bin" %in% imputed_vars) {
    lv <- levels(imps[[1]]$JCS_bin); target <- lv[length(lv)]
    miss_id <- orig[is.na(JCS_bin), .id]
    p_obs <- mean(orig[!is.na(JCS_bin), JCS_bin] == target)
    p_imp <- mean(long[.imp >= 1 & .id %in% miss_id, JCS_bin] == target)
    message(sprintf("  JCS_bin '%s' 割合: 観測=%.3f / 代入=%.3f（欠損%d例）",
                    target, p_obs, p_imp, length(miss_id)))
  }

  saveRDS(list(pool = pool, max_fmi = max_fmi, m_needed = m_needed,
               max_mc_rel_se = max_mcse, m_sufficient = m_ok,
               separation = sep_flag, range_check = range_tab),
          file.path(OUT_DIR, sprintf("bnb_imp_diag_%s.rds", nm)))
}

message("§5-1 完了。収束・分布の図を ", normalizePath(FIG_DIR),
        " に、判定結果を data/bnb_imp_diag_<pipeline>.rds に保存しました。",
        " FMI 由来の必要 m を満たさない／相対MC誤差が大きい場合は 05 の N_IMPS を増やして再実行してください。")
