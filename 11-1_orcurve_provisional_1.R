library(data.table)
library(here)
library(ggplot2)

need_pkg <- function(p) if (!requireNamespace(p, quietly = TRUE))
  stop(sprintf("パッケージ '%s' が必要です。install.packages('%s')", p, p))
need_pkg("Hmisc")
HAS_LOGISTF <- requireNamespace("logistf", quietly = TRUE)

LABELS_PATH <- here::here("00_labels.R")
if (!exists("relabel_vars")) {
  if (!file.exists(LABELS_PATH)) stop("00_labels.R が見つかりません: ", LABELS_PATH)
  source(LABELS_PATH)
}

OUT_DIR <- here::here("data")
FIG_DIR <- here::here("figures")
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

CI_LEVEL  <- 0.95
N_GRID    <- 200L
ORCURVE_CAP <- 1e3
PIPELINES <- names(PIPE_LABELS)

CURVE_COL  <- "#2166AC"
RIBBON_COL <- "#2166AC"
REF_COL    <- "grey40"
KNOT_COL   <- "grey70"

read_safe <- function(f) if (file.exists(f)) readRDS(f) else NULL

plain_num <- function(x) format(x, scientific = FALSE, trim = TRUE,
                                drop0trailing = TRUE, big.mark = "")

is_continuous_key <- function(g, groups) {
  paste0(g, "_c") %in% groups[[g]]
}

var_basis <- function(v, xc_grid, knots, n_nl, model_terms) {
  lin <- paste0(v, "_c")
  nlt <- if (n_nl >= 1L) paste0(v, "_c", seq_len(n_nl)) else character(0)
  want <- intersect(c(lin, nlt), model_terms)
  if (length(want) == 0) return(NULL)

  xc_all <- c(xc_grid, 0)
  cols <- list()
  if (lin %in% want) cols[[lin]] <- xc_all
  if (length(nlt) > 0 && any(nlt %in% want)) {
    if (is.null(knots) || length(knots) < 3L)
      stop(sprintf("変数 %s は RCS 項を持つがノットが不正です。", v))
    b <- as.matrix(Hmisc::rcspline.eval(xc_all, knots = knots, inclx = FALSE))
    for (j in seq_len(ncol(b))) {
      nm <- paste0(v, "_c", j)
      if (nm %in% want) cols[[nm]] <- b[, j]
    }
  }
  M <- do.call(cbind, cols[want])
  colnames(M) <- want
  list(B  = M[seq_along(xc_grid), , drop = FALSE],
       b0 = as.numeric(M[length(xc_all), ]),
       terms = want)
}

fit_coef_vcov <- function(smf, d, path) {
  if (path == "firth") {
    if (!HAS_LOGISTF) stop("Firth 経路だが logistf が無い。install.packages('logistf')")
    fit <- logistf::logistf(smf, data = d, pl = FALSE)
    b <- coef(fit); V <- vcov(fit)
  } else {
    fit <- suppressWarnings(glm(smf, family = binomial(), data = d))
    b <- coef(fit); V <- vcov(fit)
  }
  if (is.null(names(b)) && !is.null(rownames(V))) names(b) <- rownames(V)
  list(b = b, V = as.matrix(V))
}

or_curve_one <- function(v, imps, smf, path, center, knots, n_nl, ci_level) {
  vc <- paste0(v, "_c")
  rng <- range(unlist(lapply(imps, function(d) range(d[[vc]], na.rm = TRUE))), finite = TRUE)
  if (!all(is.finite(rng)) || rng[1] == rng[2]) return(NULL)
  xc_grid <- seq(rng[1], rng[2], length.out = N_GRID)

  m <- length(imps); G <- length(xc_grid)
  est_im <- matrix(NA_real_, nrow = m, ncol = G)
  var_im <- matrix(NA_real_, nrow = m, ncol = G)

  for (i in seq_len(m)) {
    cf <- tryCatch(fit_coef_vcov(smf, imps[[i]], path), error = function(e) NULL)
    if (is.null(cf)) next
    bb <- var_basis(v, xc_grid, knots, n_nl, names(cf$b))
    if (is.null(bb)) return(NULL)
    terms <- bb$terms
    b_sub <- cf$b[terms]
    V_sub <- cf$V[terms, terms, drop = FALSE]
    C <- sweep(bb$B, 2, bb$b0, FUN = "-")
    est_im[i, ] <- as.numeric(C %*% b_sub)
    CV <- C %*% V_sub
    var_im[i, ] <- pmax(rowSums(CV * C), 0)
  }

  ok <- colSums(is.finite(est_im)) > 0
  if (!any(ok)) return(NULL)

  Qbar <- colMeans(est_im, na.rm = TRUE)
  Ubar <- colMeans(var_im, na.rm = TRUE)
  B    <- apply(est_im, 2, function(z) { z <- z[is.finite(z)]
                if (length(z) < 2) 0 else stats::var(z) })
  Tvar <- Ubar + (1 + 1 / m) * B
  SE   <- sqrt(pmax(Tvar, 0))
  z    <- stats::qnorm(1 - (1 - ci_level) / 2)

  dt <- data.table(
    x_raw = xc_grid + center,
    logOR = Qbar, se = SE,
    OR  = exp(Qbar),
    lcl = exp(Qbar - z * SE),
    ucl = exp(Qbar + z * SE)
  )
  if (!is.null(ORCURVE_CAP) && is.finite(ORCURVE_CAP)) {
    dt[, `:=`(lcl = pmax(lcl, OR / ORCURVE_CAP),
              ucl = pmin(ucl, OR * ORCURVE_CAP))]
  }
  dt[is.finite(OR) & is.finite(lcl) & is.finite(ucl)]
}

process_pipeline <- function(nm) {
  prov <- read_safe(file.path(OUT_DIR, sprintf("bnb_provisional_%s.rds", nm)))
  obj  <- read_safe(file.path(OUT_DIR, sprintf("bnb_imp_%s.rds",        nm)))
  if (is.null(prov)) { message(sprintf("[%s] 暫定モデルが無い。スキップ。", nm)); return(NULL) }
  if (is.null(obj))  { message(sprintf("[%s] 代入データが無い。スキップ。", nm)); return(NULL) }

  smf  <- stats::as.formula(prov$smformula_prov)
  path <- prov$path
  groups <- prov$groups
  prov_vars <- prov$prov_vars

  cont_keys <- prov_vars[vapply(prov_vars, is_continuous_key, logical(1), groups = groups)]
  if (length(cont_keys) == 0) {
    message(sprintf("[%s] 暫定モデルに連続変数なし → OR カーブはスキップ。", nm))
    return(NULL)
  }
  message(sprintf("[%s] 連続変数 %d 個の OR カーブを作成: %s（経路=%s）",
                  nm, length(cont_keys),
                  paste(relabel_vars(cont_keys), collapse = ", "), path))

  imps_full   <- obj$smcfcs$impDatasets
  n0          <- prov$n0
  outlier_ids <- sort(unique(prov$outlier_ids))
  if (length(outlier_ids) > 0) {
    bad  <- outlier_ids[outlier_ids %in% seq_len(n0)]
    imps <- lapply(imps_full, function(d) d[-bad, , drop = FALSE])
  } else imps <- imps_full

  n_nl   <- max(0L, prov$n_knots - 2L)
  curves <- list()
  knot_dt <- data.table()

  for (v in cont_keys) {
    center  <- as.numeric(prov$centers[[v]])
    is_rcs  <- !is.null(prov$rcs_cont) && v %in% prov$rcs_cont
    knots   <- if (is_rcs) prov$rcs_knots[[v]] else NULL
    nnl_v   <- if (is_rcs) n_nl else 0L

    cv <- or_curve_one(v, imps, smf, path, center, knots, nnl_v, CI_LEVEL)
    if (is.null(cv) || nrow(cv) == 0) {
      message(sprintf("  [%s] %s: 曲線を計算できずスキップ。", nm, relabel_vars(v)))
      next
    }
    cv[, `:=`(pipeline = nm, variable = v,
              variable_lab = relabel_vars(v),
              type = if (is_rcs) "RCS" else "Linear")]
    curves[[v]] <- cv
    if (is_rcs && !is.null(knots))
      knot_dt <- rbind(knot_dt, data.table(variable_lab = relabel_vars(v),
                                            knot_raw = as.numeric(knots) + center))
  }
  if (length(curves) == 0) return(NULL)
  cdt <- rbindlist(curves)

  med_dt <- unique(cdt[, .(variable_lab,
                           med_raw = NA_real_)])
  for (v in cont_keys) {
    if (!v %in% names(curves)) next
    med_dt[variable_lab == relabel_vars(v), med_raw := as.numeric(prov$centers[[v]])]
  }

  p <- ggplot(cdt, aes(x = x_raw)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = REF_COL, linewidth = 0.4) +
    geom_vline(data = med_dt, aes(xintercept = med_raw),
               linetype = "longdash", colour = REF_COL, linewidth = 0.4) +
    { if (nrow(knot_dt) > 0)
        geom_vline(data = knot_dt, aes(xintercept = knot_raw),
                   linetype = "dotted", colour = KNOT_COL, linewidth = 0.3)
      else NULL } +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), fill = RIBBON_COL, alpha = 0.18) +
    geom_line(aes(y = OR), colour = CURVE_COL, linewidth = 0.8) +
    facet_wrap(~ variable_lab, scales = "free") +
    scale_y_log10(labels = plain_num) +
    labs(title = sprintf("%s  -  Provisional model: adjusted OR curves", relabel_pipeline(nm)),
         subtitle = paste0("OR vs. continuous predictor (reference = median, OR=1). ",
                           "Band = ", round(100 * CI_LEVEL), "% CI (Rubin-pooled, delta method). ",
                           "Long-dash = median; dotted = RCS knots."),
         x = "Predictor value (original scale)",
         y = "Adjusted odds ratio (log scale)") +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 8, colour = "grey30"),
          panel.grid.minor = element_blank())

  n_facet <- length(curves)
  ncol_f  <- min(3L, n_facet)
  nrow_f  <- ceiling(n_facet / ncol_f)
  ggsave(file.path(FIG_DIR, sprintf("fig4_orcurve_provisional_%s.png", nm)),
         p, width = 4.4 * ncol_f, height = 3.5 * nrow_f + 0.8, dpi = 150)
  message(sprintf("  [%s] OR カーブを保存: fig4_orcurve_provisional_%s.png（%d 変数）",
                  nm, nm, n_facet))

  cdt[]
}

all_curves <- rbindlist(lapply(PIPELINES, process_pipeline), fill = TRUE)

if (nrow(all_curves) == 0) {
  message("\n描画できる連続変数が暫定モデルにありませんでした。")
} else {
  saveRDS(all_curves, file.path(OUT_DIR, "bnb_orcurve_provisional.rds"))
  message("\n曲線データを保存: ", file.path(OUT_DIR, "bnb_orcurve_provisional.rds"))
}

message("\n§11-1 完了。暫定モデルの連続変数について OR カーブ＋95%CI を ",
        normalizePath(FIG_DIR), "/fig4_orcurve_provisional_<pipeline>.png に出力。",
        " CI は §6/§8 と同じ Rubin プール（デルタ法）。RCS 変数は非線形カーブ、",
        "線形変数は対数 OR 軸上の直線として描画。")
