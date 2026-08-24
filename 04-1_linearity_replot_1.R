library(data.table)
library(here)
library(ggplot2)

need_pkg <- function(p) if (!requireNamespace(p, quietly = TRUE))
  stop(sprintf("パッケージ '%s' が必要です。install.packages('%s')", p, p))
need_pkg("Hmisc")
need_pkg("patchwork")

LABELS_PATH <- here::here("00_labels.R")
if (!exists("relabel_vars")) {
  if (!file.exists(LABELS_PATH)) stop("00_labels.R が見つかりません: ", LABELS_PATH)
  source(LABELS_PATH)
}

OUT_DIR <- here::here("data")
FIG_DIR <- here::here("figures")
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

BASE_SIZE  <- 15
ANNO_SIZE  <- 4.3
LINE_W     <- 0.9
N_EMP_BIN  <- 10L

BL_VARS    <- c("age")

SUBSETS <- list(
  severe  = c("severe_A",  "severe_B"),
  elderly = c("elderly_A", "elderly_B")
)
SUBSET_TITLE <- c(severe = "Severe subset", elderly = "Elderly subset")

if (!exists("centered_pipelines"))
  centered_pipelines <- readRDS(file.path(OUT_DIR, "bnb_pipelines_centered.rds"))
if (!exists("rcs_meta"))
  rcs_meta <- readRDS(file.path(OUT_DIR, "bnb_rcs_meta.rds"))

assess_rows <- function(dt, v, meta) {
  cvar <- paste0(v, "_c")
  use  <- is.finite(dt[[cvar]])
  ind  <- paste0(v, "_measurable")
  if (identical(meta$mnar_system, "B") &&
      v %in% meta$mnar_vars && ind %in% names(dt)) use <- use & (dt[[ind]] == 1L)
  use
}

rcs_nl_basis <- function(x, knots) {
  if (is.null(knots) || length(unique(knots)) < length(knots) || length(knots) < 3L)
    return(NULL)
  b <- tryCatch(Hmisc::rcspline.eval(x, knots = knots, inclx = FALSE),
                error = function(e) NULL)
  if (is.null(b) || NCOL(b) < 1L) return(NULL)
  b <- as.matrix(b); colnames(b) <- paste0("nl", seq_len(ncol(b))); b
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

recompute_one <- function(dt, v, meta) {
  cvar <- paste0(v, "_c")
  use  <- assess_rows(dt, v, meta)
  x <- as.numeric(dt[[cvar]][use]); y <- as.integer(dt[["good"]][use])
  ok <- is.finite(x) & !is.na(y); x <- x[ok]; y <- y[ok]
  out <- list(pred = NULL, emp = NULL, knot_x = NULL)
  if (length(x) < 10L || length(unique(x)) < 4L) return(out)

  knots  <- meta$rcs_knots[[v]]
  out$knot_x <- if (!is.null(knots)) knots[is.finite(knots)] else NULL
  xg <- seq(min(x), max(x), length.out = 200L)

  fit_lin <- glm(y ~ x, family = binomial())
  pred <- data.table(x = xg,
                     logit = as.numeric(predict(fit_lin, data.frame(x = xg), type = "link")),
                     method = "Linear")

  basis <- rcs_nl_basis(x, knots)
  if (!is.null(basis)) {
    dfit <- data.frame(y = y, x = x, basis); nl <- colnames(basis)
    fit_rcs <- glm(stats::as.formula(paste("y ~ x +", paste(nl, collapse = " + "))),
                   data = dfit, family = binomial())
    bg <- as.matrix(Hmisc::rcspline.eval(xg, knots = knots, inclx = FALSE)); colnames(bg) <- nl
    pred <- rbind(pred, data.table(x = xg,
                  logit = as.numeric(predict(fit_rcs, data.frame(x = xg, bg), type = "link")),
                  method = "RCS"))
  }
  out$pred <- pred
  out$emp  <- empirical_logit(x, y)
  out
}

fmt_p <- function(p) ifelse(is.na(p), "n/a",
                     ifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = 3)))

build_pipeline_plot <- function(nm, dt, meta) {
  centered_present <- sub("_c$", "", grep("_c$", names(dt), value = TRUE))
  contv <- intersect(meta$selected_cont, centered_present)

  res <- lapply(contv, function(v) recompute_one(dt, v, meta)); names(res) <- contv

  pred_all <- data.table(); emp_all <- data.table(); knot_all <- data.table()
  for (v in contv) {
    r <- res[[v]]; vlab <- relabel_vars(v)
    if (!is.null(r$pred)) pred_all <- rbind(pred_all, cbind(r$pred, variable = vlab))
    if (!is.null(r$emp))  emp_all  <- rbind(emp_all,  cbind(r$emp,  variable = vlab))
    if (!is.null(r$knot_x))
      knot_all <- rbind(knot_all, data.table(variable = vlab, knot_x = r$knot_x))
  }

  lrt <- as.data.table(meta$linearity_lrt)
  anno <- data.table()
  for (v in contv) {
    row <- lrt[variable == v]
    tag <- if (v %in% meta$rcs_cont) "  [RCS applied]" else ""
    txt <- if (nrow(row) == 0 || !isTRUE(row$feasible) || is.na(row$p))
      paste0("LRT: not estimable", tag)
    else
      sprintf("LRT: chi-sq=%.1f (df=%d)\np=%s%s",
              row$LR_chisq, row$df, fmt_p(row$p), tag)
    bl <- v %in% BL_VARS
    anno <- rbind(anno, data.table(
      variable = relabel_vars(v), label = txt,
      ypos = if (bl) -Inf else Inf,
      vj   = if (bl) -0.25 else 1.15))
  }

  ggplot() +
    geom_vline(data = knot_all, aes(xintercept = knot_x),
               linetype = "dotted", colour = "grey70", linewidth = 0.35) +
    geom_point(data = emp_all, aes(x = x, y = logit, size = n),
               colour = "grey35", alpha = 0.7) +
    geom_line(data = pred_all,
              aes(x = x, y = logit, colour = method, linetype = method),
              linewidth = LINE_W) +
    geom_text(data = anno, aes(x = -Inf, y = ypos, label = label, vjust = vj),
              hjust = -0.04, size = ANNO_SIZE, lineheight = 0.95,
              colour = "grey15") +
    facet_wrap(~ variable, scales = "free") +
    scale_colour_manual(values = c(Linear = "#2166AC", RCS = "#B2182B")) +
    scale_linetype_manual(values = c(Linear = "dashed", RCS = "solid")) +
    scale_size_continuous(range = c(1.2, 4), guide = "none") +
    labs(title = relabel_pipeline(nm), x = "Centered predictor (0 = median)",
         y = "logit (log-odds of Good Outcome)", colour = "Fit", linetype = "Fit") +
    theme_bw(base_size = BASE_SIZE) +
    theme(plot.title = element_text(size = BASE_SIZE + 1, face = "bold"),
          strip.text = element_text(size = BASE_SIZE - 1),
          legend.text = element_text(size = BASE_SIZE - 1),
          legend.title = element_text(size = BASE_SIZE - 1))
}

for (sub in names(SUBSETS)) {
  pipes <- SUBSETS[[sub]]
  plots <- list()
  nfac  <- 0L
  for (nm in pipes) {
    if (is.null(centered_pipelines[[nm]]) || is.null(rcs_meta[[nm]])) {
      warning("入力が見つかりません: ", nm, " をスキップ"); next
    }
    dt <- centered_pipelines[[nm]]; meta <- rcs_meta[[nm]]
    plots[[nm]] <- build_pipeline_plot(nm, dt, meta)
    centered_present <- sub("_c$", "", grep("_c$", names(dt), value = TRUE))
    nfac <- max(nfac, length(intersect(meta$selected_cont, centered_present)))
  }
  if (length(plots) == 0) next

  ncol_f <- min(4L, max(1L, nfac))
  nrow_f <- ceiling(nfac / ncol_f)

  probs_txt <- paste(rcs_meta[[pipes[1]]]$knot_probs, collapse = ", ")
  combo <- patchwork::wrap_plots(plots, ncol = 1) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      title = sprintf("%s  -  linearity check (logit vs centered predictor)",
                      SUBSET_TITLE[[sub]]),
      subtitle = sprintf(
        "Dashed = linear, solid = RCS; grey points = continuity-corrected empirical logit (size = bin n); dotted = knots at quantiles %s",
        probs_txt),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(size = BASE_SIZE + 4, face = "bold"),
        plot.subtitle = ggplot2::element_text(size = BASE_SIZE - 2, colour = "grey30"))) &
    ggplot2::theme(legend.position = "top")

  per_h <- nrow_f * 2.9 + 0.7
  ggsave(file.path(FIG_DIR, sprintf("fig2_linearity_%s.png", sub)),
         combo, width = 3.6 * ncol_f, height = per_h * length(plots) + 1.0,
         dpi = 150, limitsize = FALSE)
  message(sprintf("[%s] %s を保存（%d パイプライン × %d 変数）",
                  sub, sprintf("fig2_linearity_%s.png", sub), length(plots), nfac))
}

message("§4-1 完了。サブセット統合版の線形性図を ", normalizePath(FIG_DIR), " に保存しました。")
