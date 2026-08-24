library(data.table)
library(here)
library(ggplot2)

need_pkg <- function(p) if (!requireNamespace(p, quietly = TRUE))
  stop(sprintf("パッケージ '%s' が必要です。install.packages('%s')", p, p))
need_pkg("patchwork")

LABELS_PATH <- here::here("00_labels.R")
if (!exists("relabel_vars")) {
  if (!file.exists(LABELS_PATH)) stop("00_labels.R が見つかりません: ", LABELS_PATH)
  source(LABELS_PATH)
}

OUT_DIR <- here::here("data")
FIG_DIR <- here::here("figures")
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

BASE_SIZE <- 15
LINE_W    <- 0.5

SUBSETS <- list(severe  = c("severe_A",  "severe_B"),
                elderly = c("elderly_A", "elderly_B"))
SUBSET_TITLE <- c(severe = "Severe subset", elderly = "Elderly subset")

imp_index <- readRDS(file.path(OUT_DIR, "bnb_imp_index.rds"))

pretty_term <- function(cn, factor_vars) {
  fv <- factor_vars[order(nchar(factor_vars), decreasing = TRUE)]
  lab <- function(b) if (b %in% names(VAR_LABELS)) unname(VAR_LABELS[b]) else b
  vapply(cn, function(s) {
    if (s == "(Intercept)") return("Intercept")
    m <- regmatches(s, regexec("^(.*)_c([0-9]+)$", s))[[1]]
    if (length(m) == 3) return(paste0(lab(m[2]), " (spline ", m[3], ")"))
    if (grepl("_c$", s)) return(lab(sub("_c$", "", s)))
    mm <- regmatches(s, regexec("^(.*)_measurable[0-9]*$", s))[[1]]
    if (length(mm) == 2) return(paste0(lab(mm[2]), " (measurable)"))
    for (b in fv) if (startsWith(s, b)) {
      lev <- substring(s, nchar(b) + 1)
      return(paste0(lab(b), ": ", relabel_levels(b, lev)))
    }
    if (s %in% names(VAR_LABELS)) return(lab(s))
    s
  }, character(1), USE.NAMES = FALSE)
}

trace_dt <- function(sci, term_labels) {
  d <- dim(sci); if (is.null(d) || length(d) != 3) return(NULL)
  nch <- d[1]; nco <- d[2]; nit <- d[3]
  if (length(term_labels) != nco) term_labels <- paste0("b", seq_len(nco))
  out <- rbindlist(lapply(seq_len(nch), function(c0)
    rbindlist(lapply(seq_len(nco), function(j)
      data.table(chain = c0, term = term_labels[j],
                 iter = seq_len(nit), est = as.numeric(sci[c0, j, ]))))))
  out[, term := factor(term, levels = unique(term_labels))]
  out
}

build_conv_plot <- function(nm) {
  obj  <- readRDS(file.path(OUT_DIR, sprintf("bnb_imp_%s.rds", nm)))
  sci  <- obj$smcfcs$smCoefIter
  if (is.null(dim(sci))) return(NULL)
  imps <- obj$smcfcs$impDatasets
  fit1 <- suppressWarnings(glm(stats::as.formula(obj$smformula),
                               family = binomial(), data = imps[[1]]))
  cn   <- names(coef(fit1))
  fvars <- names(Filter(is.factor, imps[[1]]))
  td   <- trace_dt(sci, pretty_term(cn, fvars))
  if (is.null(td)) return(NULL)
  multi <- max(td$chain) > 1

  ggplot(td, aes(iter, est, group = factor(chain), colour = factor(chain))) +
    geom_line(linewidth = LINE_W, alpha = 0.85) +
    facet_wrap(~ term, scales = "free_y") +
    labs(title = relabel_pipeline(nm), x = "Iteration",
         y = "Coefficient estimate", colour = "Chain") +
    theme_bw(base_size = BASE_SIZE) +
    theme(plot.title  = element_text(size = BASE_SIZE + 1, face = "bold"),
          strip.text  = element_text(size = BASE_SIZE - 2),
          axis.text   = element_text(size = BASE_SIZE - 3),
          legend.position = if (multi) "top" else "none")
}

build_dist_plot <- function(nm) {
  obj  <- readRDS(file.path(OUT_DIR, sprintf("bnb_imp_%s.rds", nm)))
  long <- as.data.table(obj$long); orig <- long[.imp == 0]
  imps <- obj$smcfcs$impDatasets
  cont <- Filter(function(v) v %in% names(orig) && !is.factor(imps[[1]][[v]]),
                 obj$imputed)

  dd <- data.table()
  for (v in cont) {
    miss <- orig[is.na(get(v)), .id]
    obs  <- orig[!is.na(get(v)), as.numeric(get(v))]
    imp  <- long[.imp >= 1 & .id %in% miss, as.numeric(get(v))]
    if (length(imp) == 0) next
    dd <- rbind(dd,
      data.table(variable = relabel_vars(v), source = "Observed", value = obs),
      data.table(variable = relabel_vars(v), source = "Imputed",  value = imp))
  }
  if ("JCS_bin" %in% obj$imputed && is.factor(imps[[1]]$JCS_bin)) {
    lv <- levels(imps[[1]]$JCS_bin); tg <- lv[length(lv)]
    miss <- orig[is.na(JCS_bin), .id]
    message(sprintf("  [%s] JCS '%s' 割合 観測=%.3f / 代入=%.3f",
                    nm, relabel_levels("JCS_bin", tg),
                    mean(orig[!is.na(JCS_bin), JCS_bin] == tg),
                    mean(long[.imp >= 1 & .id %in% miss, JCS_bin] == tg)))
  }
  if (nrow(dd) == 0) return(NULL)

  ggplot(dd, aes(value, colour = source, fill = source)) +
    geom_density(alpha = 0.25, linewidth = 0.7) +
    facet_wrap(~ variable, scales = "free") +
    scale_colour_manual(values = c(Observed = "#2166AC", Imputed = "#B2182B")) +
    scale_fill_manual(values   = c(Observed = "#2166AC", Imputed = "#B2182B")) +
    labs(title = relabel_pipeline(nm), x = NULL, y = "Density",
         colour = NULL, fill = NULL) +
    theme_bw(base_size = BASE_SIZE) +
    theme(plot.title = element_text(size = BASE_SIZE + 1, face = "bold"),
          strip.text = element_text(size = BASE_SIZE - 1),
          legend.position = "top")
}

n_panels <- function(p) {
  b <- ggplot_build(p)
  max(as.integer(b$data[[length(b$data)]]$PANEL))
}

combine_and_save <- function(sub, pipes, builder, fname,
                             subtitle, ncol_facet, per_facet, panel_unit_h) {
  plots <- list(); max_fac <- 1L
  for (nm in pipes) {
    if (!file.exists(file.path(OUT_DIR, sprintf("bnb_imp_%s.rds", nm)))) {
      warning("入力なし: ", nm, " をスキップ"); next
    }
    p <- builder(nm)
    if (is.null(p)) { warning(nm, " は描画不能"); next }
    plots[[nm]] <- p
    max_fac <- max(max_fac, n_panels(p))
  }
  if (length(plots) == 0) return(invisible(NULL))

  nrow_f <- ceiling(max_fac / ncol_facet)
  combo <- patchwork::wrap_plots(plots, ncol = 1) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      title = sprintf("%s  -  %s", SUBSET_TITLE[[sub]], subtitle),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(size = BASE_SIZE + 4, face = "bold"))) &
    ggplot2::theme(legend.position = "top")

  per_h <- nrow_f * panel_unit_h + 0.7
  ggsave(file.path(FIG_DIR, fname), combo,
         width = per_facet * ncol_facet, height = per_h * length(plots) + 1.0,
         dpi = 150, limitsize = FALSE)
  message(sprintf("[%s] %s を保存（%d パイプライン）", sub, fname, length(plots)))
}

for (sub in names(SUBSETS)) {
  pipes <- intersect(SUBSETS[[sub]], imp_index)
  message(sprintf("================ [%s] 診断図 ================", sub))

  combine_and_save(sub, pipes, build_conv_plot,
                   sprintf("fig_imp_convergence_%s.png", sub),
                   subtitle = "smcfcs convergence trace (estimate vs iteration)",
                   ncol_facet = 5L, per_facet = 3.0, panel_unit_h = 2.5)

  combine_and_save(sub, pipes, build_dist_plot,
                   sprintf("fig_imp_distributions_%s.png", sub),
                   subtitle = "Observed vs imputed (post-truncation), continuous MAR variables",
                   ncol_facet = 2L, per_facet = 4.4, panel_unit_h = 3.4)
}

message("§5-1 完了。サブセット統合の convergence / distribution 図を ",
        normalizePath(FIG_DIR), " に保存しました。")
