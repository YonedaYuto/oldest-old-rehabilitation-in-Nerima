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

INT_P_CRIT <- 0.05
INTERACTION_RCS_MODE <- "linear"
SKIP_MEASURABLE_SELF_PAIR <- TRUE
RECHECK_SEP_PER_INTERACTION <- TRUE
GLM_LRT_METHOD <- "D3"

BETA_MAX <- 15
SE_MAX   <- 30
PROB_EPS <- 1e-6

FORCE_PATH <- list(
  severe_A  = NULL, severe_B  = NULL,
  elderly_A = NULL, elderly_B = NULL
)

M_SCREEN_MAX <- NULL

LOGISTF_PL <- FALSE

if (!exists("prov_index"))
  prov_index <- readRDS(file.path(OUT_DIR, "bnb_provisional_index.rds"))

detect_one <- function(smf, dat) {
  detector <- if (HAS_DETECTSEP) detectseparation::detect_separation
              else                brglm2::detect_separation
  ds <- tryCatch(stats::glm(smf, data = dat, family = binomial("logit"), method = detector),
                 error = function(e) e)
  if (inherits(ds, "error"))
    return(list(separated = NA, inf_terms = character(0), ok = FALSE))
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
                inf_terms = inf_terms, det_error = det_err))
  n_unstable <- 0L
  for (d in imps) {
    f <- suppressWarnings(stats::glm(smf, binomial(), d))
    if (isTRUE(glm_unstable_one(f))) n_unstable <- n_unstable + 1L
  }
  list(sep = n_unstable > 0, kind = if (n_unstable > 0) "near/unstable" else "none",
       n_unstable = n_unstable, m = length(imps), inf_terms = inf_terms, det_error = det_err)
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
  if (m == 1) return(list(p = stats::pchisq(mbar, df = k, lower.tail = FALSE),
                          D = mbar / k, df1 = k, df2 = Inf, r = 0))
  r <- (1 + 1/m) * stats::var(sqrt(stats)); if (!is.finite(r) || r < 0) r <- 0
  D <- (mbar/k - (m + 1)/(m - 1) * r) / (1 + r); if (!is.finite(D) || D < 0) D <- 0
  df2 <- if (r <= 0) 1e6 else k^(-3/m) * (m - 1) * (1 + 1/r)^2
  if (!is.finite(df2) || df2 <= 0) df2 <- 1e6
  list(p = stats::pf(D, df1 = k, df2 = df2, lower.tail = FALSE), D = D, df1 = k, df2 = df2, r = r)
}

pool_lrt_D3_glm <- function(fits_full, fits_red) {
  tryCatch({
    res <- mice::D3(mice::as.mira(fits_full), mice::as.mira(fits_red))
    tab <- as.data.frame(res$result)
    pcol <- grep("^P", names(tab))
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
    pll_full <- vapply(imps, function(d) {
      fit <- tryCatch(logistf::logistf(smf_full, data = d, pl = LOGISTF_PL), error = function(e) NULL)
      if (is.null(fit)) NA_real_ else firth_pll(fit) }, numeric(1))
    pll_red <- vapply(imps, function(d) {
      fit <- tryCatch(logistf::logistf(smf_red, data = d, pl = LOGISTF_PL), error = function(e) NULL)
      if (is.null(fit)) NA_real_ else firth_pll(fit) }, numeric(1))
    lr <- 2 * (pll_full - pll_red)
    f1f <- tryCatch(logistf::logistf(smf_full, data = imps[[1]], pl = FALSE), error = function(e) NULL)
    f1r <- tryCatch(logistf::logistf(smf_red,  data = imps[[1]], pl = FALSE), error = function(e) NULL)
    k <- if (!is.null(f1f) && !is.null(f1r)) length(coef(f1f)) - length(coef(f1r)) else NA_integer_
    return(c(pool_lrt_D2(lr, k), list(k = k, engine = "D2(Firth-PLR)")))
  }
}

inter_side_labels <- function(g, groups, rcs_mode) {
  labs <- groups[[g]]
  if (identical(rcs_mode, "full")) return(labs)
  lin <- labs[grepl("_c$", labs)]
  if (length(lin) >= 1) lin else labs
}

make_inter_labels <- function(g1, g2, groups, rcs_mode) {
  s1 <- inter_side_labels(g1, groups, rcs_mode)
  s2 <- inter_side_labels(g2, groups, rcs_mode)
  as.vector(outer(s1, s2, function(a, b) paste0(a, ":", b)))
}

base_underlying <- function(g) sub("_measurable$", "", g)

process_pipeline <- function(nm) {
  prov <- readRDS(file.path(OUT_DIR, sprintf("bnb_provisional_%s.rds", nm)))
  obj  <- readRDS(file.path(OUT_DIR, sprintf("bnb_imp_%s.rds",        nm)))

  groups      <- prov$groups
  prov_vars   <- prov$prov_vars
  prov_labels <- prov$prov_labels
  outlier_ids <- prov$outlier_ids
  n0 <- prov$n0; m <- length(obj$smcfcs$impDatasets)

  forced <- FORCE_PATH[[nm]]
  path <- if (!is.null(forced)) { if (!forced %in% c("glm","firth")) stop("FORCE_PATH は 'glm'/'firth'。"); forced }
          else prov$path

  message(sprintf("\n================ [%s] §9 交互作用スクリーニング ================", nm))
  message("  暫定モデル: ", prov$smformula_prov)
  message("  暫定モデルの変数: ",
          if (length(prov_vars)) paste(relabel_vars(prov_vars), collapse = ", ") else "(切片のみ)")
  message(sprintf("  §8 経路=%s → §9 採用経路=%s | 交互作用形式=%s",
                  toupper(prov$path), toupper(path), INTERACTION_RCS_MODE))

  bad <- sort(unique(outlier_ids)); bad <- bad[bad %in% seq_len(n0)]
  imps <- if (length(bad) > 0) lapply(obj$smcfcs$impDatasets, function(d) d[-bad, , drop = FALSE])
          else obj$smcfcs$impDatasets
  n <- n0 - length(bad)
  message(sprintf("  縮小データ: n=%d（除外 %d 件）, m=%d", n, length(bad), m))

  screen_imps <- if (!is.null(M_SCREEN_MAX) && M_SCREEN_MAX < m) imps[seq_len(M_SCREEN_MAX)] else imps
  if (length(screen_imps) < m)
    message(sprintf("  探索は m=%d 代入に間引き。", length(screen_imps)))

  base_firth <- (path == "firth")
  if (!base_firth) {
    base_su <- sep_or_unstable(make_formula(prov_labels), screen_imps)
    if (isTRUE(base_su$sep)) {
      base_firth <- TRUE
      message(sprintf("  暫定(ベース)モデルが分離相当(%s) → 全候補を Firth で評価（判定スキップ）。",
                      base_su$kind))
    }
  }
  path_eff <- if (base_firth) "firth" else "glm"

  if (length(prov_vars) < 2L) {
    message("  暫定モデルの変数が 2 未満 → 交互作用候補なし。")
    cand_dt <- data.table(); kept <- list(); path9 <- path_eff
  } else {
    pairs <- utils::combn(prov_vars, 2L, simplify = FALSE)
    if (SKIP_MEASURABLE_SELF_PAIR)
      pairs <- Filter(function(pr) base_underlying(pr[1]) != base_underlying(pr[2]), pairs)
    message(sprintf("  交互作用候補: %d ペア", length(pairs)))

    if (have_tic) tictoc::tic(paste0(nm, " interaction screen"))
    rows <- lapply(pairs, function(pr) {
      g1 <- pr[1]; g2 <- pr[2]
      int_labels <- make_inter_labels(g1, g2, groups, INTERACTION_RCS_MODE)
      aug_labels <- c(prov_labels, int_labels)

      if (base_firth) {
        cand_path <- "firth"; sep_induced <- TRUE
      } else if (RECHECK_SEP_PER_INTERACTION) {
        su <- sep_or_unstable(make_formula(aug_labels), screen_imps)
        sep_induced <- isTRUE(su$sep)
        cand_path   <- if (sep_induced) "firth" else "glm"
      } else {
        cand_path <- "glm"; sep_induced <- FALSE
      }
      lr <- pooled_lrt(aug_labels, prov_labels, screen_imps, path = cand_path)
      data.table(v1 = g1, v2 = g2, n_terms = length(int_labels), df = lr$k,
                 engine = lr$engine, sep_induced = sep_induced, p = lr$p,
                 int_labels = paste(int_labels, collapse = " + "))
    })
    cand_dt <- rbindlist(rows)[order(p)]
    cand_dt[, kept := is.finite(p) & p < INT_P_CRIT]
    if (have_tic) tictoc::toc()

    message("  候補LRT（p昇順）:")
    print(cand_dt[, .(interaction = paste(relabel_vars(v1), "x", relabel_vars(v2)),
                      df, engine, sep_induced, p = signif(p, 4), kept)])

    keep_rows <- cand_dt[kept == TRUE]
    kept <- lapply(seq_len(nrow(keep_rows)), function(i) {
      r <- keep_rows[i]
      g1 <- r$v1; g2 <- r$v2
      list(key = paste0(g1, ":", g2), pair = c(g1, g2),
           labels = make_inter_labels(g1, g2, groups, INTERACTION_RCS_MODE),
           p = r$p, df = r$df, engine = r$engine, sep_induced = r$sep_induced)
    })
    names(kept) <- vapply(kept, function(z) z$key, character(1))

    path9 <- if (base_firth || any(keep_rows$sep_induced)) "firth" else "glm"

    message(sprintf("  → 残す交互作用: %d 件%s",
                    length(kept),
                    if (length(kept)) paste0("（", paste(vapply(kept, function(z)
                      paste(relabel_vars(z$pair[1]), "x", relabel_vars(z$pair[2])), character(1)),
                      collapse = ", "), "）") else ""))
    if (any(keep_rows$sep_induced))
      message("  ※残った交互作用に、拡大モデルで分離を誘発したものあり → §10 は Firth 経路を考慮。")
  }

  kept_labels <- if (length(kept)) unlist(lapply(kept, function(z) z$labels), use.names = FALSE) else character(0)
  kept_int_groups <- if (length(kept)) setNames(lapply(kept, function(z) z$labels), names(kept)) else list()

  saveRDS(list(
    pipeline   = nm,
    path       = path_eff,
    path_nominal = path,
    path_hint  = path9,
    prov_vars  = prov_vars, prov_labels = prov_labels, groups = groups,
    smformula_prov = prov$smformula_prov,
    candidates = cand_dt,
    kept       = kept,
    kept_keys  = names(kept),
    kept_labels = kept_labels,
    kept_int_groups = kept_int_groups,
    settings   = list(int_p_crit = INT_P_CRIT, rcs_mode = INTERACTION_RCS_MODE,
                      skip_measurable_self = SKIP_MEASURABLE_SELF_PAIR,
                      recheck_sep = RECHECK_SEP_PER_INTERACTION, glm_lrt = GLM_LRT_METHOD,
                      beta_max = BETA_MAX, se_max = SE_MAX, prob_eps = PROB_EPS),
    outlier_ids = outlier_ids, n0 = n0, n = n, m = m,
    mnar_system = prov$mnar_system,
    selected_cont = prov$selected_cont, selected_cat = prov$selected_cat,
    rcs_cont = prov$rcs_cont, rcs_knots = prov$rcs_knots, n_knots = prov$n_knots,
    centers = prov$centers, mnar_indicator = prov$mnar_indicator, basis_cols = prov$basis_cols
  ), file.path(OUT_DIR, sprintf("bnb_interactions_%s.rds", nm)))

  n_cand <- if (is.data.table(cand_dt)) nrow(cand_dt) else 0L
  message(sprintf("[%s] §9 完了。候補 %d / 残存 %d を bnb_interactions_%s.rds に保存。",
                  nm, n_cand, length(kept), nm))

  data.table(pipeline = nm, mnar_system = prov$mnar_system, n = n, m = m,
             path = path_eff, path_nominal = path, base_upgraded = (path_eff != path),
             n_prov_vars = length(prov_vars),
             n_candidates = n_cand, n_kept = length(kept),
             any_sep_induced = if (n_cand > 0) any(cand_dt$kept & cand_dt$sep_induced) else FALSE,
             kept = if (length(kept))
               paste(vapply(kept, function(z) paste0(z$pair[1], "*", z$pair[2]), character(1)), collapse = ", ")
               else "")
}

int_index <- rbindlist(lapply(prov_index$pipeline, process_pipeline), fill = TRUE)
saveRDS(int_index, file.path(OUT_DIR, "bnb_interactions_index.rds"))

message("\n==== §9 サマリ（交互作用スクリーニング）====")
print(int_index[, .(pipeline, n, path, base_upgraded, n_prov_vars, n_candidates, n_kept, kept)])

message("\n§9 完了。各パイプラインの残存交互作用を bnb_interactions_<pipeline>.rds に保存しました。",
        " 候補探索では §10 と同じ統合判定 sep_or_unstable（厳密な完全/準完全分離 ＋ glm 直接診断による",
        " near-separation）で(準)分離を評価し、分離相当の候補は Firth で尤度比検定を行いました。",
        " §10 は kept_int_groups を候補として、主効果固定のまま前進選択で交互作用を選びます。",
        " path（=path_eff）は near-separation 評価でベースが firth に昇格した結果を反映します。")
