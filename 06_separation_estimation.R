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
FIG_DIR <- here::here("figures")
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

CI_LEVEL <- 0.95

FORCE_PATH <- list(
  severe_A  = NULL, severe_B  = NULL,
  elderly_A = NULL, elderly_B = NULL
)

CLIP_VARS <- list(
  severe_A  = NULL, severe_B  = NULL,
  elderly_A = NULL, elderly_B = NULL
)

GEN_CLIP_PROFILE <- TRUE
LOGISTF_PL       <- FALSE

have_tic <- requireNamespace("tictoc", quietly = TRUE)

if (!exists("imp_index"))
  imp_index <- readRDS(file.path(OUT_DIR, "bnb_imp_index.rds"))

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

pool_rubin <- function(fitlist, n_complete, ci_level = CI_LEVEL) {
  m     <- length(fitlist)
  betas <- lapply(fitlist, coef)
  terms <- names(betas[[1]])
  Vdg   <- lapply(fitlist, function(f) diag(as.matrix(vcov(f))))

  getv <- function(lst, t) vapply(lst, function(z) unname(z[t]), numeric(1))
  Qbar <- vapply(terms, function(t) mean(getv(betas, t)), numeric(1))
  B    <- vapply(terms, function(t) stats::var(getv(betas, t)), numeric(1))
  Ubar <- vapply(terms, function(t) mean(getv(Vdg, t)), numeric(1))

  Tvar   <- Ubar + (1 + 1/m) * B
  riv    <- (1 + 1/m) * B / Ubar
  lambda <- (1 + 1/m) * B / Tvar
  k      <- length(terms); dfc <- max(n_complete - k, 1)
  df_old <- ifelse(lambda > 0, (m - 1) / lambda^2, Inf)
  df_obs <- ((dfc + 1) / (dfc + 3)) * dfc * (1 - lambda)
  df     <- df_old * df_obs / (df_old + df_obs)
  fmi    <- (riv + 2 / (df + 3)) / (riv + 1)

  se <- sqrt(Tvar)
  z  <- stats::qnorm(1 - (1 - ci_level) / 2)
  data.table(term = terms, estimate = Qbar, se = se,
             OR = exp(Qbar), lcl = exp(Qbar - z * se), ucl = exp(Qbar + z * se),
             riv = riv, fmi = fmi)
}

estimate_glm <- function(smf, imps, n, ci_level = CI_LEVEL) {
  sep_warn <- FALSE
  fits <- lapply(imps, function(d) withCallingHandlers(
    glm(smf, family = binomial(), data = d),
    warning = function(w) {
      if (grepl("fitted probabilities", conditionMessage(w))) sep_warn <<- TRUE
      invokeRestart("muffleWarning")
    }))

  pooled <- mice::pool(fits)
  ps <- summary(pooled, conf.int = TRUE, conf.level = ci_level)
  ps <- as.data.table(ps)

  ci_cols <- grep("%", names(ps), value = TRUE)
  if (length(ci_cols) >= 2) {
    lo <- ps[[ci_cols[1]]]; hi <- ps[[ci_cols[2]]]
  } else {
    zz <- stats::qnorm(1 - (1 - ci_level) / 2)
    lo <- ps$estimate - zz * ps$std.error; hi <- ps$estimate + zz * ps$std.error
  }
  out <- data.table(
    term = ps$term, estimate = ps$estimate, se = ps$std.error,
    OR = exp(ps$estimate), lcl = exp(lo), ucl = exp(hi),
    p = ps$p.value)
  list(pooled = out, fits = fits, glm_sep_warn = sep_warn)
}

estimate_firth <- function(smf, imps, n, clip_vars = NULL, ci_level = CI_LEVEL) {
  fitlist <- lapply(imps, function(d)
    logistf::logistf(smf, data = d, pl = LOGISTF_PL, dataout = TRUE))

  rb <- pool_rubin(fitlist, n_complete = n, ci_level = ci_level)

  if (have_tic) tictoc::tic("  CLIP.confint")
  cc <- tryCatch(
    logistf::CLIP.confint(obj = fitlist, variable = clip_vars,
                          ci.level = ci_level, pvalue = TRUE),
    error = function(e) e)
  if (have_tic) tictoc::toc()

  pooled <- copy(rb)[, .(term, estimate, se, OR, lcl, ucl, fmi)]
  clip_ok <- !inherits(cc, "error")

  if (clip_ok) {
    ci_mat <- if (!is.null(cc$ci)) as.matrix(cc$ci)
              else if (!is.null(cc$confint)) as.matrix(cc$confint) else NULL
    vnm <- if (!is.null(cc$variable)) cc$variable
           else if (!is.null(rownames(ci_mat))) rownames(ci_mat) else NULL
    if (!is.null(ci_mat) && !is.null(vnm)) {
      clip_dt <- data.table(term = vnm,
                            lcl_clip = ci_mat[, 1], ucl_clip = ci_mat[, 2],
                            p = if (!is.null(cc$pvalue)) cc$pvalue else NA_real_)
      pooled <- merge(pooled, clip_dt, by = "term", all.x = TRUE, sort = FALSE)
      have_clip <- is.finite(pooled$lcl_clip) & is.finite(pooled$ucl_clip)
      pooled[have_clip, `:=`(lcl = exp(lcl_clip), ucl = exp(ucl_clip))]
      pooled[, c("lcl_clip", "ucl_clip") := NULL]
    } else {
      clip_ok <- FALSE
    }
  }
  if (!clip_ok) {
    warning("CLIP.confint が失敗/未対応のため、CIは罰則付きRubinの正規近似を使用します。")
    pooled[, p := NA_real_]
  }

  list(pooled = pooled, fits = fitlist, clip = cc, clip_ok = clip_ok)
}

save_clip_profile <- function(fitlist, variable, nm) {
  pr <- tryCatch(
    logistf::CLIP.profile(obj = fitlist, variable = variable),
    error = function(e) e)
  if (inherits(pr, "error")) {
    message("  CLIP.profile 図はスキップ（", conditionMessage(pr), "）")
    return(invisible(NULL))
  }
  fp <- file.path(FIG_DIR, sprintf("fig_clip_profile_%s.png", nm))
  grDevices::png(fp, width = 1200, height = 900, res = 150)
  on.exit(grDevices::dev.off(), add = TRUE)
  tryCatch(
    plot(pr, main = sprintf("%s  -  CLIP profile: %s",
                            relabel_pipeline(nm), relabel_vars(variable))),
    error = function(e) message("  CLIP.profile の plot に失敗: ", conditionMessage(e)))
  message("  CLIP プロファイル図を保存: ", fp)
}

est_index <- data.table()

for (nm in imp_index) {
  obj  <- readRDS(file.path(OUT_DIR, sprintf("bnb_imp_%s.rds", nm)))
  imps <- obj$smcfcs$impDatasets
  smf  <- stats::as.formula(obj$smformula)
  n    <- obj$n
  m    <- length(imps)

  message(sprintf("================ [%s] 分離判定→推定分岐 (n=%d, m=%d) ================",
                  nm, n, m))
  message("  smformula: ", obj$smformula)

  det <- lapply(imps, function(d) detect_one(smf, d))
  sep_vec   <- vapply(det, function(z) isTRUE(z$separated), logical(1))
  ok_vec    <- vapply(det, function(z) isTRUE(z$ok),        logical(1))
  inf_all   <- sort(unique(unlist(lapply(det, function(z) z$inf_terms))))
  n_sep     <- sum(sep_vec, na.rm = TRUE)
  any_sep   <- n_sep > 0

  if (any(!ok_vec))
    message(sprintf("  ※分離検出に失敗した代入が %d 件（モデル特異等の可能性）。",
                    sum(!ok_vec)))
  message(sprintf("  分離あり代入: %d / %d  (%.0f%%)", n_sep, m, 100 * n_sep / m))
  if (any_sep)
    message("  無限化した係数（いずれかの代入で）: ",
            paste(inf_all, collapse = ", "))

  sep_table <- data.table(
    imp        = seq_len(m),
    separated  = sep_vec,
    detect_ok  = ok_vec,
    inf_terms  = vapply(det, function(z) paste(z$inf_terms, collapse = "; "), character(1))
  )

  forced <- FORCE_PATH[[nm]]
  path <- if (!is.null(forced)) {
    if (!forced %in% c("glm", "firth")) stop("FORCE_PATH は 'glm'/'firth' のみ。")
    message("  推定経路を手動指定: ", forced,
            if (forced == "glm" && any_sep) "（分離ありだが glm を強制）" else "")
    forced
  } else if (any_sep) "firth" else "glm"

  message(sprintf("  → 採用する推定経路: %s", toupper(path)))

  if (have_tic) tictoc::tic(nm)
  if (path == "glm") {
    est <- estimate_glm(smf, imps, n, ci_level = CI_LEVEL)
    pooled <- est$pooled
    if (isTRUE(est$glm_sep_warn))
      message("  ※glm で fitted prob 0/1 の警告あり（準完全分離の兆候）。",
              "FORCE_PATH=\"firth\" での感度確認を推奨。")
    method_detail <- list(estimator = "glm", pool = "mice::pool (Rubin)")
  } else {
    est <- estimate_firth(smf, imps, n, clip_vars = CLIP_VARS[[nm]], ci_level = CI_LEVEL)
    pooled <- est$pooled
    method_detail <- list(estimator = "logistf (Firth)",
                          pool = "Rubin (point) + CLIP.confint (CI/p)",
                          clip_ok = est$clip_ok)
    if (GEN_CLIP_PROFILE) {
      cand <- setdiff(inf_all, "(Intercept)")
      if (length(cand) == 0)
        cand <- setdiff(pooled$term, "(Intercept)")
      if (length(cand) > 0) save_clip_profile(est$fits, cand[1], nm)
    }
  }
  if (have_tic) tictoc::toc()

  show <- copy(pooled)[, `:=`(
    OR = round(OR, 3), lcl = round(lcl, 3), ucl = round(ucl, 3),
    estimate = round(estimate, 3), se = round(se, 3))]
  if ("p" %in% names(show)) show[, p := signif(p, 4)]
  show[, label := relabel_vars(sub("(運動器|廃用|有|無|女|男|あり|なし|impaired|Yes|No|Female|Male|Needed|Independent)$", "", term))]
  message("  フルモデル（選択前）プール推定:")
  print(show[, intersect(c("term", "OR", "lcl", "ucl", "p", "fmi"), names(show)), with = FALSE])

  saveRDS(list(
    pipeline      = nm,
    path          = path,
    pooled        = pooled,
    smformula     = obj$smformula,
    separation    = list(any = any_sep, n_sep = n_sep, m = m,
                         inf_terms = inf_all, table = sep_table),
    method_detail = method_detail,
    n = n, m = m, mnar_system = obj$mnar_system,
    selected_cont = obj$selected_cont, selected_cat = obj$selected_cat,
    rcs_cont = obj$rcs_cont, rcs_knots = obj$rcs_knots, n_knots = obj$n_knots,
    centers = obj$centers, mnar_indicator = obj$mnar_indicator
  ), file.path(OUT_DIR, sprintf("bnb_estimation_%s.rds", nm)))

  est_index <- rbind(est_index, data.table(
    pipeline = nm, path = path, mnar_system = obj$mnar_system,
    n = n, m = m, n_sep = n_sep, any_separation = any_sep,
    inf_terms = paste(inf_all, collapse = "; ")))

  message(sprintf("[%s] 完了。経路=%s を bnb_estimation_%s.rds に保存。", nm, path, nm))
}

saveRDS(est_index, file.path(OUT_DIR, "bnb_estimation_index.rds"))

message("\n==== §6 サマリ（推定経路）====")
print(est_index[, .(pipeline, mnar_system, n, m, n_sep, path)])
message("§6 完了。各パイプラインの推定経路を確定しました。",
        " GLM 経路は §8 で psfmi(D3) により選択、",
        " FIRTH 経路は §8 で logistf の罰則付き尤度比を手動プール（D2/D3・Median-P）します。",
        " 分離時の CLIP プロファイル図は ", normalizePath(FIG_DIR), " に保存されます。")
