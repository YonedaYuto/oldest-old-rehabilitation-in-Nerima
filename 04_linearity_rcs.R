library(data.table)
library(here)
library(ggplot2)

need_pkg <- function(p) if (!requireNamespace(p, quietly = TRUE))
  stop(sprintf("パッケージ '%s' が必要です。install.packages('%s')", p, p))
need_pkg("Hmisc")

LABELS_PATH <- here::here("00_labels.R")
if (!exists("relabel_vars")) {
  if (!file.exists(LABELS_PATH)) stop("00_labels.R が見つかりません: ", LABELS_PATH)
  source(LABELS_PATH)
}

OUT_DIR <- here::here("data")
FIG_DIR <- here::here("figures")
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

N_KNOTS <- 3L

KNOT_PROBS <- list(
  "3" = c(0.10, 0.50, 0.90),
  "4" = c(0.05, 0.35, 0.65, 0.95),
  "5" = c(0.05, 0.275, 0.50, 0.725, 0.95)
)
if (!as.character(N_KNOTS) %in% names(KNOT_PROBS))
  stop("N_KNOTS は 3, 4, 5 のいずれかにしてください。")
PROBS <- KNOT_PROBS[[as.character(N_KNOTS)]]

N_EMP_BIN <- 10L

RCS_VARS <- list(
  severe_A  = character(0),
  severe_B  = character(0),
  elderly_A = character(0),
  elderly_B = character(0)
)

RCS_VARS <- list(
  severe_A  = c("STEF_better_in", "BBS_in"),
  severe_B  = c("STEF_better_in", "BBS_in"),
  elderly_A = c("STEF_better_in"),
  elderly_B = c("STEF_better_in")
)

if (!exists("centered_pipelines"))
  centered_pipelines <- readRDS(file.path(OUT_DIR, "bnb_pipelines_centered.rds"))
if (!exists("selected_meta"))
  selected_meta      <- readRDS(file.path(OUT_DIR, "bnb_selected_meta.rds"))

assess_rows <- function(dt, v, meta) {
  cvar <- paste0(v, "_c")
  use  <- is.finite(dt[[cvar]])
  ind  <- paste0(v, "_measurable")
  if (identical(meta$mnar_system, "B") &&
      v %in% meta$mnar_vars && ind %in% names(dt)) {
    use <- use & (dt[[ind]] == 1L)
  }
  use
}

rcs_nl_basis <- function(x, knots) {
  if (length(unique(knots)) < length(knots) || length(knots) < 3L) return(NULL)
  b <- tryCatch(Hmisc::rcspline.eval(x, knots = knots, inclx = FALSE),
                error = function(e) NULL)
  if (is.null(b) || NCOL(b) < 1L) return(NULL)
  b <- as.matrix(b)
  colnames(b) <- paste0("nl", seq_len(ncol(b)))
  b
}

empirical_logit <- function(x, y, nbin = N_EMP_BIN) {
  br <- unique(stats::quantile(x, probs = seq(0, 1, length.out = nbin + 1), na.rm = TRUE))
  if (length(br) < 3L) return(NULL)
  g  <- cut(x, breaks = br, include.lowest = TRUE)
  dt <- data.table(x = x, y = y, g = g)[!is.na(g)]
  agg <- dt[, .(xm = median(x), ev = sum(y), n = .N), by = g]
  agg[, logit := log((ev + 0.5) / (n - ev + 0.5))]
  agg[, .(x = xm, logit, n)]
}

assess_one <- function(dt, v, meta) {
  cvar <- paste0(v, "_c")
  use  <- assess_rows(dt, v, meta)
  x <- as.numeric(dt[[cvar]][use])
  y <- as.integer(dt[["good"]][use])
  ok <- is.finite(x) & !is.na(y)
  x <- x[ok]; y <- y[ok]
  n <- length(x)

  out <- list(var = v, n = n, n_events = sum(y),
              knots = NA_real_, feasible = FALSE,
              lr_chisq = NA_real_, df = NA_integer_, p = NA_real_,
              pred = NULL, emp = NULL, knot_x = NULL)

  if (n < 10L || length(unique(x)) < 4L) {
    warning(sprintf("[%s] 変数 %s: 有効例/ユニーク値が少なく評価不能（n=%d）。",
                    meta$mnar_system, v, n))
    return(out)
  }

  knots <- as.numeric(stats::quantile(x, probs = PROBS, na.rm = TRUE))
  out$knots  <- knots
  out$knot_x <- knots

  xg <- seq(min(x), max(x), length.out = 200L)

  fit_lin   <- glm(y ~ x, family = binomial())
  lin_logit <- as.numeric(predict(fit_lin,
                                  newdata = data.frame(x = xg), type = "link"))
  pred <- data.table(x = xg, logit = lin_logit, method = "Linear")

  basis <- rcs_nl_basis(x, knots)
  if (!is.null(basis)) {
    dfit    <- data.frame(y = y, x = x, basis)
    nlnames <- colnames(basis)
    form    <- stats::as.formula(paste("y ~ x +", paste(nlnames, collapse = " + ")))
    fit_rcs <- glm(form, data = dfit, family = binomial())

    an <- anova(fit_lin, fit_rcs, test = "LRT")
    out$feasible <- TRUE
    out$df       <- as.integer(an$Df[2])
    out$lr_chisq <- as.numeric(an$Deviance[2])
    out$p        <- as.numeric(an$`Pr(>Chi)`[2])

    basis_g  <- Hmisc::rcspline.eval(xg, knots = knots, inclx = FALSE)
    basis_g  <- as.matrix(basis_g); colnames(basis_g) <- nlnames
    rcs_logit <- as.numeric(predict(fit_rcs,
                  newdata = data.frame(x = xg, basis_g), type = "link"))
    pred <- rbind(pred,
                  data.table(x = xg, logit = rcs_logit, method = "RCS"))
  } else {
    warning(sprintf("[%s] 変数 %s: ノット分位点が縮退し RCS 不可。線形のみ表示。",
                    meta$mnar_system, v))
  }

  out$pred <- pred
  out$emp  <- empirical_logit(x, y)
  out
}

make_linearity_figure <- function(results, meta, label) {
  disp <- relabel_pipeline(label)

  pred_all <- data.table(); emp_all <- data.table(); knot_all <- data.table()
  for (r in results) {
    vlab <- relabel_vars(r$var)
    if (!is.null(r$pred) && nrow(r$pred) > 0)
      pred_all <- rbind(pred_all, cbind(r$pred, variable = vlab))
    if (!is.null(r$emp) && nrow(r$emp) > 0)
      emp_all  <- rbind(emp_all,  cbind(r$emp,  variable = vlab))
    if (!is.null(r$knot_x))
      knot_all <- rbind(knot_all,
                        data.table(variable = vlab, knot_x = r$knot_x))
  }
  if (nrow(pred_all) == 0) return(NULL)

  p <- ggplot() +
    geom_vline(data = knot_all, aes(xintercept = knot_x),
               linetype = "dotted", colour = "grey70", linewidth = 0.3) +
    geom_point(data = emp_all, aes(x = x, y = logit, size = n),
               colour = "grey35", alpha = 0.7) +
    geom_line(data = pred_all,
              aes(x = x, y = logit, colour = method, linetype = method),
              linewidth = 0.7) +
    facet_wrap(~ variable, scales = "free") +
    scale_colour_manual(values = c(Linear = "#2166AC", RCS = "#B2182B")) +
    scale_linetype_manual(values = c(Linear = "dashed", RCS = "solid")) +
    scale_size_continuous(range = c(0.8, 3), guide = "none") +
    labs(title = sprintf("%s  -  Linearity check (logit vs centered predictor)", disp),
         subtitle = sprintf("Knots fixed at quantiles %s (N_KNOTS=%d); points = continuity-corrected empirical logit",
                            paste(PROBS, collapse = ", "), N_KNOTS),
         x = "Centered predictor (0 = median)", y = "logit (log-odds of Good Outcome)",
         colour = "Fit", linetype = "Fit") +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(size = 11),
          plot.subtitle = element_text(size = 8, colour = "grey30"),
          legend.position = "top")
  p
}

linearity_results <- list()
rcs_knots_all      <- list()
rcs_meta           <- list()

for (nm in names(centered_pipelines)) {
  dt   <- centered_pipelines[[nm]]
  meta <- selected_meta[[nm]]

  centered_present <- sub("_c$", "", grep("_c$", names(dt), value = TRUE))
  contv <- intersect(meta$selected_cont, centered_present)

  res_list <- lapply(contv, function(v) assess_one(dt, v, meta))
  names(res_list) <- contv

  tab <- rbindlist(lapply(res_list, function(r) data.table(
    variable = r$var, label = relabel_vars(r$var),
    n = r$n, events = r$n_events,
    df = r$df, LR_chisq = round(r$lr_chisq, 3), p = signif(r$p, 4),
    nonlinear = ifelse(is.na(r$p), NA, r$p < 0.05),
    feasible = r$feasible
  )))
  linearity_results[[nm]] <- tab

  rcs_knots_all[[nm]] <- lapply(res_list, function(r) r$knots)

  fig <- make_linearity_figure(res_list, meta, nm)
  if (!is.null(fig)) {
    n_facet <- length(contv)
    ncol_f  <- min(3L, n_facet)
    nrow_f  <- ceiling(n_facet / ncol_f)
    ggsave(file.path(FIG_DIR, sprintf("fig2_linearity_%s.png", nm)),
           fig, width = 4.2 * ncol_f, height = 3.4 * nrow_f + 0.6, dpi = 150)
  }

  message(sprintf("==== [%s] 非線形性 LRT（complete-case 単変量 screening） ====", nm))
  print(tab)
  flagged <- tab[nonlinear == TRUE, variable]
  if (length(flagged) > 0)
    message(sprintf("[%s] p<0.05 で非線形の示唆: %s  → 図2を見て RCS_VARS への記入を検討。",
                    nm, paste(flagged, collapse = ", ")))

  rcs_cont <- intersect(RCS_VARS[[nm]], contv)
  feasible_vars <- tab[feasible == TRUE, variable]
  bad <- setdiff(rcs_cont, feasible_vars)
  if (length(bad) > 0) {
    warning(sprintf("[%s] RCS 指定だが RCS 不可のため除外: %s",
                    nm, paste(bad, collapse = ", ")))
    rcs_cont <- intersect(rcs_cont, feasible_vars)
  }

  m <- meta
  m$n_knots   <- N_KNOTS
  m$knot_probs <- PROBS
  m$rcs_knots <- rcs_knots_all[[nm]]
  m$rcs_cont  <- rcs_cont
  m$linearity_lrt <- tab
  rcs_meta[[nm]] <- m

  if (length(rcs_cont) > 0)
    message(sprintf("[%s] RCS 化する変数: %s（%d ノット）",
                    nm, paste(rcs_cont, collapse = ", "), N_KNOTS))
  else
    message(sprintf("[%s] RCS 化する変数なし（全連続変数を線形で投入）", nm))
}

saveRDS(rcs_knots_all,     file.path(OUT_DIR, "bnb_rcs_knots.rds"))
saveRDS(linearity_results, file.path(OUT_DIR, "bnb_linearity_results.rds"))
saveRDS(rcs_meta,          file.path(OUT_DIR, "bnb_rcs_meta.rds"))

message("§4 完了。図2を ", normalizePath(FIG_DIR),
        " に保存。図2を見て RCS_VARS（変数ごと）を記入し再実行すると、",
        "RCS 化の判定が bnb_rcs_meta.rds に反映され §5 の代入式に渡されます。")
