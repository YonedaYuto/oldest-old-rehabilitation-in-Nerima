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

CI_LEVEL    <- 0.95
N_GRID      <- 200L
MOD_PROBS   <- c(0.10, 0.50, 0.90)
ORCURVE_CAP <- 1e3
PIPELINES   <- names(PIPE_LABELS)

STRATA_PAL <- c("#1B7837", "#762A83", "#E08214", "#2166AC", "#B2182B", "#666666")

read_safe <- function(f) if (file.exists(f)) readRDS(f) else NULL

classify_term <- function(term, var_names) {
  if (term == "(Intercept)")
    return(list(var = NA_character_, type = "intercept"))
  if (grepl("_measurable$", term))
    return(list(var = term, type = "indicator"))
  if (grepl("_c[0-9]+$", term))
    return(list(var = sub("_c[0-9]+$", "", term), type = "cont_nl",
                j = as.integer(sub(".*_c", "", term))))
  if (grepl("_c$", term))
    return(list(var = sub("_c$", "", term), type = "cont_lin"))
  hit <- var_names[vapply(var_names, function(v) startsWith(term, v), logical(1))]
  if (length(hit) == 0) return(list(var = term, type = "factor", level = ""))
  v <- hit[which.max(nchar(hit))]
  list(var = v, type = "factor", level = sub(paste0("^", v), "", term))
}

component_value <- function(cls, assigned, centers, knots) {
  switch(cls$type,
    cont_lin  = as.numeric(assigned) - centers[[cls$var]],
    cont_nl   = {
      vc <- as.numeric(assigned) - centers[[cls$var]]
      kn <- knots[[cls$var]]
      if (is.null(kn)) stop(sprintf("変数 %s の RCS ノットがありません。", cls$var))
      as.numeric(Hmisc::rcspline.eval(vc, knots = kn, inclx = FALSE)[, cls$j])
    },
    indicator = as.numeric(assigned),
    factor    = as.numeric(as.character(assigned) == cls$level),
    intercept = rep(0, length(assigned)),
    stop("未知の項タイプ: ", cls$type)
  )
}

build_A <- function(coef_names, focal, focal_grid, focal_ref,
                    mod, mod_assign, ref_assign, var_names, centers, knots) {
  G <- length(focal_grid)
  A <- matrix(0, nrow = G, ncol = length(coef_names),
              dimnames = list(NULL, coef_names))
  for (cn in coef_names) {
    comps <- strsplit(cn, ":", fixed = TRUE)[[1]]
    cls   <- lapply(comps, classify_term, var_names = var_names)
    fpos  <- which(vapply(cls, function(z) identical(z$var, focal), logical(1)))
    if (length(fpos) == 0) next
    fc <- cls[[fpos[1]]]
    fval_grid <- component_value(fc, focal_grid, centers, knots)
    fval_ref  <- component_value(fc, focal_ref,  centers, knots)
    oval <- 1
    others <- setdiff(seq_along(comps), fpos[1])
    for (oi in others) {
      oc <- cls[[oi]]
      assigned <- if (identical(oc$var, mod)) mod_assign else ref_assign(oc$var)
      oval <- oval * component_value(oc, assigned, centers, knots)
    }
    A[, cn] <- (fval_grid - fval_ref) * oval
  }
  A
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

pool_curve <- function(A_list, imps, smf, path, ci_level) {
  m <- length(imps)
  coef_names <- colnames(A_list[[1]])
  fits <- vector("list", m)
  for (i in seq_len(m)) {
    cf <- tryCatch(fit_coef_vcov(smf, imps[[i]], path), error = function(e) NULL)
    if (is.null(cf)) next
    b <- cf$b[coef_names]; V <- cf$V[coef_names, coef_names, drop = FALSE]
    fits[[i]] <- list(b = b, V = V)
  }
  fits <- Filter(Negate(is.null), fits)
  m_ok <- length(fits)
  if (m_ok == 0) return(NULL)
  z <- stats::qnorm(1 - (1 - ci_level) / 2)

  out <- vector("list", length(A_list))
  for (s in seq_along(A_list)) {
    A <- A_list[[s]]; G <- nrow(A)
    est <- matrix(NA_real_, m_ok, G); var <- matrix(NA_real_, m_ok, G)
    for (i in seq_len(m_ok)) {
      b <- fits[[i]]$b; V <- fits[[i]]$V
      est[i, ] <- as.numeric(A %*% b)
      AV <- A %*% V
      var[i, ] <- pmax(rowSums(AV * A), 0)
    }
    Qbar <- colMeans(est, na.rm = TRUE); Ubar <- colMeans(var, na.rm = TRUE)
    B <- apply(est, 2, function(zz) { zz <- zz[is.finite(zz)]
               if (length(zz) < 2) 0 else stats::var(zz) })
    SE <- sqrt(pmax(Ubar + (1 + 1 / m_ok) * B, 0))
    out[[s]] <- data.table(OR = exp(Qbar),
                           lcl = exp(Qbar - z * SE),
                           ucl = exp(Qbar + z * SE))
  }
  out
}

build_interaction <- function(key, fin, imps, smf, path, coef_names,
                               var_names, centers, knots) {
  parts <- strsplit(key, ":", fixed = TRUE)[[1]]
  if (length(parts) != 2) return(NULL)
  g1 <- parts[1]; g2 <- parts[2]; groups <- fin$groups

  kind <- function(g) {
    if (paste0(g, "_c") %in% groups[[g]]) "cont"
    else if (grepl("_measurable$", g))   "indicator"
    else                                  "factor"
  }
  k1 <- kind(g1); k2 <- kind(g2)

  ref_assign <- function(v) {
    if (v %in% names(centers))            return(centers[[v]])
    if (grepl("_measurable$", v))         return(0)
    lv <- levels(imps[[1]][[v]]); if (is.null(lv)) lv <- sort(unique(as.character(imps[[1]][[v]])))
    lv[1]
  }
  fac_levels <- function(v) { lv <- levels(imps[[1]][[v]])
    if (is.null(lv)) sort(unique(as.character(imps[[1]][[v]]))) else lv }
  cont_range <- function(v) { vc <- paste0(v, "_c")
    range(unlist(lapply(imps, function(d) range(d[[vc]], na.rm = TRUE))), finite = TRUE) }
  cont_strata <- function(v) { vc <- paste0(v, "_c")
    x <- unlist(lapply(imps, function(d) d[[vc]]))
    qs <- as.numeric(stats::quantile(x, MOD_PROBS, na.rm = TRUE))
    data.table(raw = qs + centers[[v]], lab = sprintf("%s = %.4g", relabel_vars(v), qs + centers[[v]]))
  }

  pick <- function() {
    if (k1 == "cont" && k2 != "cont") return(list(f = g1, m = g2, fk = "cont"))
    if (k2 == "cont" && k1 != "cont") return(list(f = g2, m = g1, fk = "cont"))
    if (k1 == "cont" && k2 == "cont") {
      f <- if (!is.null(fin$rcs_cont) && g1 %in% fin$rcs_cont) g1
           else if (!is.null(fin$rcs_cont) && g2 %in% fin$rcs_cont) g2 else g1
      list(f = f, m = setdiff(c(g1, g2), f), fk = "cont")
    } else {
      n1 <- length(fac_levels(g1)); n2 <- length(fac_levels(g2))
      f <- if (n1 >= n2) g1 else g2
      list(f = f, m = setdiff(c(g1, g2), f), fk = "cat")
    }
  }
  pk <- pick(); focal <- pk$f; mod <- pk$m

  mk <- kind(mod)
  strata <- if (mk == "cont") {
    cs <- cont_strata(mod); lapply(seq_len(nrow(cs)), function(i) list(assign = cs$raw[i], lab = cs$lab[i]))
  } else if (mk == "indicator") {
    list(list(assign = 0, lab = paste0(relabel_vars(mod), ": No")),
         list(assign = 1, lab = paste0(relabel_vars(mod), ": Yes")))
  } else {
    lv <- fac_levels(mod); lapply(lv, function(l) list(assign = l, lab = relabel_levels(mod, l)))
  }

  disp <- sprintf("%s x %s", relabel_vars(g1), relabel_vars(g2))

  if (pk$fk == "cont") {
    rng <- cont_range(focal)
    if (!all(is.finite(rng)) || rng[1] == rng[2]) return(NULL)
    xc <- seq(rng[1], rng[2], length.out = N_GRID)
    focal_grid <- xc + centers[[focal]]
    focal_ref  <- centers[[focal]]
    A_list <- lapply(strata, function(st)
      build_A(coef_names, focal, focal_grid, focal_ref, mod, st$assign,
              ref_assign, var_names, centers, knots))
    res <- pool_curve(A_list, imps, smf, path, CI_LEVEL)
    if (is.null(res)) return(NULL)
    dt <- rbindlist(lapply(seq_along(strata), function(s)
      cbind(res[[s]], data.table(x_raw = focal_grid, stratum = strata[[s]]$lab))))
    list(type = "curve", focal = focal, mod = mod, disp = disp,
         med_raw = centers[[focal]],
         x_lab = relabel_vars(focal), data = dt)

  } else {
    levs <- fac_levels(focal); ref_lev <- levs[1]; show_levs <- setdiff(levs, ref_lev)
    if (length(show_levs) == 0) return(NULL)
    focal_grid <- show_levs
    focal_ref  <- ref_lev
    A_list <- lapply(strata, function(st)
      build_A(coef_names, focal, focal_grid, focal_ref, mod, st$assign,
              ref_assign, var_names, centers, knots))
    res <- pool_curve(A_list, imps, smf, path, CI_LEVEL)
    if (is.null(res)) return(NULL)
    dt <- rbindlist(lapply(seq_along(strata), function(s)
      cbind(res[[s]], data.table(focal_level = vapply(show_levs, function(l) relabel_levels(focal, l), character(1)),
                                 stratum = strata[[s]]$lab))))
    list(type = "points", focal = focal, mod = mod, disp = disp,
         ref_lab = relabel_levels(focal, ref_lev),
         x_lab = relabel_vars(focal), data = dt)
  }
}

clamp_ci <- function(dt) {
  if (!is.null(ORCURVE_CAP) && is.finite(ORCURVE_CAP))
    dt[, `:=`(lcl = pmax(lcl, OR / ORCURVE_CAP), ucl = pmin(ucl, OR * ORCURVE_CAP))]
  dt[is.finite(OR) & is.finite(lcl) & is.finite(ucl)]
}

plot_curve <- function(info, pipe_disp) {
  dt <- clamp_ci(copy(info$data))
  ggplot(dt, aes(x = x_raw, colour = stratum, fill = stratum)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    geom_vline(xintercept = info$med_raw, linetype = "longdash", colour = "grey60", linewidth = 0.35) +
    geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.13, colour = NA) +
    geom_line(aes(y = OR), linewidth = 0.8) +
    scale_colour_manual(values = STRATA_PAL, name = relabel_vars(info$mod)) +
    scale_fill_manual(values = STRATA_PAL, guide = "none") +
    scale_y_log10() +
    labs(title = sprintf("%s  -  Interaction: %s", pipe_disp, info$disp),
         subtitle = sprintf("Adjusted OR of %s by strata of %s (ref = %s median; OR=1). Band = %d%% CI (Rubin, delta).",
                            info$x_lab, relabel_vars(info$mod), info$x_lab, round(100 * CI_LEVEL)),
         x = sprintf("%s (original scale)", info$x_lab),
         y = "Adjusted odds ratio (log scale)") +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 8, colour = "grey30"),
          legend.position = "top", panel.grid.minor = element_blank())
}

plot_points <- function(info, pipe_disp) {
  dt <- clamp_ci(copy(info$data))
  pd <- position_dodge(width = 0.5)
  ggplot(dt, aes(x = focal_level, y = OR, colour = stratum)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0.18, position = pd, linewidth = 0.5) +
    geom_point(size = 2.2, position = pd) +
    scale_colour_manual(values = STRATA_PAL, name = relabel_vars(info$mod)) +
    scale_y_log10() +
    labs(title = sprintf("%s  -  Interaction: %s", pipe_disp, info$disp),
         subtitle = sprintf("Adjusted OR of %s (vs. %s) by strata of %s; OR=1. Bars = %d%% CI (Rubin, delta).",
                            info$x_lab, info$ref_lab, relabel_vars(info$mod), round(100 * CI_LEVEL)),
         x = info$x_lab, y = "Adjusted odds ratio (log scale)") +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 8, colour = "grey30"),
          legend.position = "top", panel.grid.minor = element_blank())
}

sanitize <- function(s) gsub("[^A-Za-z0-9]+", "_", s)

process_pipeline <- function(nm) {
  fin <- read_safe(file.path(OUT_DIR, sprintf("bnb_final_%s.rds", nm)))
  obj <- read_safe(file.path(OUT_DIR, sprintf("bnb_imp_%s.rds",   nm)))
  if (is.null(fin)) { message(sprintf("[%s] 最終モデルが無い。スキップ。", nm)); return(NULL) }
  if (is.null(obj)) { message(sprintf("[%s] 代入データが無い。スキップ。", nm)); return(NULL) }

  int_keys <- fin$selected_int_keys
  if (length(int_keys) == 0) {
    message(sprintf("[%s] 最終モデルに交互作用なし → スキップ。", nm)); return(NULL)
  }
  message(sprintf("[%s] 交互作用 %d 件を図示: %s（経路=%s）", nm, length(int_keys),
                  paste(int_keys, collapse = ", "), fin$path))

  smf  <- stats::as.formula(fin$smformula_final)
  path <- fin$path
  centers <- fin$centers; knots <- fin$rcs_knots

  imps_full   <- obj$smcfcs$impDatasets
  n0          <- fin$n0
  outlier_ids <- sort(unique(fin$outlier_ids))
  if (length(outlier_ids) > 0) {
    bad  <- outlier_ids[outlier_ids %in% seq_len(n0)]
    imps <- lapply(imps_full, function(d) d[-bad, , drop = FALSE])
  } else imps <- imps_full

  ref_fit <- tryCatch(fit_coef_vcov(smf, imps[[1]], path), error = function(e) NULL)
  if (is.null(ref_fit)) { message(sprintf("[%s] 参照フィット失敗。スキップ。", nm)); return(NULL) }
  coef_names <- names(ref_fit$b)
  var_names  <- unique(c(names(centers), fin$selected_cat, names(fin$groups)))
  var_names  <- var_names[order(-nchar(var_names))]

  collected <- list()
  for (key in int_keys) {
    info <- tryCatch(build_interaction(key, fin, imps, smf, path, coef_names,
                                       var_names, centers, knots),
                     error = function(e) { message(sprintf("  [%s] %s: %s", nm, key, conditionMessage(e))); NULL })
    if (is.null(info)) { message(sprintf("  [%s] %s: 図を作れずスキップ。", nm, key)); next }
    p <- if (info$type == "curve") plot_curve(info, relabel_pipeline(nm))
         else                       plot_points(info, relabel_pipeline(nm))
    fn <- file.path(FIG_DIR, sprintf("fig4_interaction_%s_%s.png", nm, sanitize(key)))
    h  <- if (info$type == "curve") 4.2 else max(3.6, 0.0 + 4.0)
    ggsave(fn, p, width = 7.6, height = h, dpi = 150)
    message(sprintf("  [%s] %s を保存: %s", nm, info$disp, basename(fn)))
    d <- copy(info$data); d[, `:=`(pipeline = nm, interaction = info$disp, type = info$type)]
    collected[[key]] <- d
  }
  if (length(collected) == 0) return(NULL)
  rbindlist(collected, fill = TRUE)
}

all_int <- rbindlist(lapply(PIPELINES, process_pipeline), fill = TRUE)

if (nrow(all_int) == 0) {
  message("\nどのパイプラインの最終モデルにも図示すべき交互作用がありませんでした。")
} else {
  saveRDS(all_int, file.path(OUT_DIR, "bnb_interaction_plots.rds"))
  message("\n交互作用の図データを保存: ", file.path(OUT_DIR, "bnb_interaction_plots.rds"))
}

message("\n§11-2 完了。最終モデルに残った交互作用のみを ",
        normalizePath(FIG_DIR), "/fig4_interaction_<pipeline>_<key>.png に出力。",
        " 焦点が連続→層別 OR カーブ、カテゴリ→層別 OR 点。各層内の焦点参照点を OR=1 に揃え、",
        " 層間の開きが交互作用（効果修飾）を表す。CI は §11-1 と同じ Rubin プール（デルタ法）。")
