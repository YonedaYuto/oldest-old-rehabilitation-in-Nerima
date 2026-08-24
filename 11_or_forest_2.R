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

CI_LEVEL    <- 0.95
DIGITS_OR   <- 3
DIGITS_BETA <- 3

FOREST_OR_LIM <- c(0.05, 20)

PIPELINES <- names(PIPE_LABELS)

relabel_coef <- function(cn) {
  cn <- as.character(cn)
  one <- function(s) {
    s <- as.character(s)
    if (is.na(s) || !nzchar(s)) return(s)
    if (s == "(Intercept)") return("(Intercept)")
    if (grepl(":", s, fixed = TRUE)) {
      parts <- strsplit(s, ":", fixed = TRUE)[[1]]
      return(paste(vapply(parts, one, character(1)), collapse = " x "))
    }
    if (grepl("_c[0-9]+$", s)) {
      base <- sub("_c[0-9]+$", "", s); j <- sub(".*_c", "", s)
      return(paste0(relabel_vars(base), " (nl", j, ")"))
    }
    if (grepl("_c$", s))         return(relabel_vars(s))
    if (grepl("_measurable", s)) return(relabel_vars(s))
    hit <- names(VAR_LABELS)[vapply(names(VAR_LABELS),
                                    function(v) startsWith(s, v), logical(1))]
    if (length(hit) > 0) {
      v   <- hit[which.max(nchar(hit))]
      lev <- sub(paste0("^", v), "", s)
      if (nzchar(lev)) return(paste0(relabel_vars(v), ": ", relabel_levels(v, lev)))
      return(relabel_vars(v))
    }
    s
  }
  vapply(cn, one, character(1), USE.NAMES = FALSE)
}

tidy_pooled <- function(pooled, pipeline, model, path) {
  if (is.null(pooled) || nrow(pooled) == 0) return(data.table())
  dt <- as.data.table(pooled)
  dt[, term := as.character(term)]
  getcol <- function(nm) if (nm %in% names(dt)) dt[[nm]] else NA_real_
  out <- data.table(
    pipeline      = pipeline,
    pipeline_disp = relabel_pipeline(pipeline),
    model         = model,
    path          = path,
    term          = dt$term,
    variable      = relabel_coef(dt$term),
    beta          = as.numeric(getcol("estimate")),
    se            = as.numeric(getcol("se")),
    OR            = as.numeric(getcol("OR")),
    lcl           = as.numeric(getcol("lcl")),
    ucl           = as.numeric(getcol("ucl")),
    p             = as.numeric(getcol("p")),
    fmi           = as.numeric(getcol("fmi"))
  )
  out[, sig := is.finite(lcl) & is.finite(ucl) & (lcl > 1 | ucl < 1) &
        term != "(Intercept)"]
  out[]
}

read_safe <- function(f) if (file.exists(f)) readRDS(f) else NULL

coef_list <- list()

for (nm in PIPELINES) {
  prov <- read_safe(file.path(OUT_DIR, sprintf("bnb_provisional_%s.rds", nm)))
  fin  <- read_safe(file.path(OUT_DIR, sprintf("bnb_final_%s.rds",       nm)))

  if (is.null(prov))
    message(sprintf("[%s] 暫定モデル(bnb_provisional_%s.rds)が見つかりません。スキップ。", nm, nm))
  else
    coef_list[[paste0(nm, "_prov")]] <-
      tidy_pooled(prov$pooled, nm, "Provisional", prov$path)

  if (is.null(fin))
    message(sprintf("[%s] 最終モデル(bnb_final_%s.rds)が見つかりません。スキップ。", nm, nm))
  else
    coef_list[[paste0(nm, "_final")]] <-
      tidy_pooled(fin$pooled, nm, "Final", fin$path)
}

coef_dt <- rbindlist(coef_list, fill = TRUE)
if (nrow(coef_dt) == 0)
  stop("暫定/最終モデルの保存結果が読めませんでした。§8・§10 を先に実行してください。")

coef_dt[, pipeline := factor(pipeline, levels = PIPELINES)]
coef_dt[, model    := factor(model,    levels = c("Provisional", "Final"))]
coef_dt[, term_order := seq_len(.N), by = .(pipeline, model)]
setorder(coef_dt, pipeline, model, term_order)

fmt_or  <- function(x) ifelse(is.finite(x), formatC(x, format = "f", digits = DIGITS_OR), "NA")
fmt_b   <- function(x) ifelse(is.finite(x), formatC(x, format = "f", digits = DIGITS_BETA), "NA")
fmt_p   <- function(x) ifelse(is.finite(x), ifelse(x < 1e-4, "<0.0001", signif(x, 3)), "NA")

show_dt <- coef_dt[, .(
  Pipeline = pipeline_disp,
  Model    = as.character(model),
  Path     = path,
  Variable = variable,
  beta     = fmt_b(beta),
  SE       = fmt_b(se),
  OR       = fmt_or(OR),
  `95%CI`  = sprintf("[%s, %s]", fmt_or(lcl), fmt_or(ucl)),
  p        = fmt_p(p),
  FMI      = ifelse(is.finite(fmi), formatC(fmi, format = "f", digits = 2), "")
)]

message("\n==== §11.1 暫定・最終モデル 係数・CI 統合表 ====")
message("（係数 beta = log-OR、OR は exp(beta)、95%CI は OR スケール。",
        "Path: glm=Rubin / firth=Firth+CLIP）\n")
for (pp in PIPELINES) {
  for (md in c("Provisional", "Final")) {
    sub <- show_dt[Pipeline == relabel_pipeline(pp) & Model == md]
    if (nrow(sub) == 0) next
    message(sprintf("---- %s | %s model (path=%s) ----",
                    relabel_pipeline(pp), md, sub$Path[1]))
    print(sub[, .(Variable, beta, SE, OR, `95%CI`, p, FMI)], row.names = FALSE)
    cat("\n")
  }
}

saveRDS(coef_dt, file.path(OUT_DIR, "bnb_coef_table.rds"))
fwrite(coef_dt[, .(pipeline, pipeline_disp, model, path, term, variable,
                   beta, se, OR, lcl, ucl, p, fmi, sig)],
       file.path(OUT_DIR, "bnb_coef_table.csv"))
message("係数・CI 統合表を保存: ",
        file.path(OUT_DIR, "bnb_coef_table.rds"), " / .csv")

HILITE <- "#B2182B"

forest_dt <- coef_dt[model == "Provisional" & term != "(Intercept)" &
                       is.finite(OR) & is.finite(lcl) & is.finite(ucl)]

if (nrow(forest_dt) == 0) {
  message("暫定モデルに描画できる項がありません（切片のみ等）。forest plot はスキップ。")
} else {
  lo <- FOREST_OR_LIM[1]; hi <- FOREST_OR_LIM[2]
  forest_dt[, `:=`(
    lcl_disp   = pmax(lcl, lo),
    ucl_disp   = pmin(ucl, hi),
    OR_disp    = pmin(pmax(OR, lo), hi),
    trunc_low  = lcl < lo,
    trunc_high = ucl > hi
  )]

  FS_BASE   <- 14
  FS_TITLE  <- 17
  FS_SUB    <- 11
  FS_AXIS   <- 14
  FS_AXIST  <- 14
  FS_CITEXT <- 4.8
  PT_SIZE   <- 2.8
  EB_LW     <- 0.7

  fmt1 <- function(x) ifelse(!is.finite(x), "Inf",
            ifelse(abs(x) >= 1000, formatC(x, format = "g", digits = 3),
                   formatC(x, format = "f", digits = 2)))
  forest_dt[, or_text := sprintf("%s [%s, %s]", fmt1(OR), fmt1(lcl), fmt1(ucl))]

  make_forest <- function(d, title, show_legend = TRUE, show_sub = TRUE) {
    setorder(d, OR_disp)
    d[, ylab := factor(variable, levels = variable)]

    forest <- ggplot(d, aes(x = OR_disp, y = ylab, colour = sig)) +
      geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
      geom_errorbarh(aes(xmin = lcl_disp, xmax = ucl_disp),
                     height = 0.24, linewidth = EB_LW) +
      geom_point(size = PT_SIZE) +
      { if (any(d$trunc_low))
          geom_segment(data = d[trunc_low == TRUE],
                       aes(x = lcl_disp * 1.02, xend = lo, y = ylab, yend = ylab),
                       arrow = grid::arrow(length = grid::unit(0.045, "npc")),
                       linewidth = EB_LW, inherit.aes = FALSE, colour = "grey30")
        else NULL } +
      { if (any(d$trunc_high))
          geom_segment(data = d[trunc_high == TRUE],
                       aes(x = ucl_disp * 0.98, xend = hi, y = ylab, yend = ylab),
                       arrow = grid::arrow(length = grid::unit(0.045, "npc")),
                       linewidth = EB_LW, inherit.aes = FALSE, colour = "grey30")
        else NULL } +
      scale_colour_manual(values = c(`FALSE` = "grey45", `TRUE` = HILITE),
                          labels = c(`FALSE` = "CI includes 1",
                                     `TRUE`  = "CI excludes 1"),
                          name = NULL, drop = FALSE) +
      scale_x_log10(limits = FOREST_OR_LIM) +
      labs(title = title,
           subtitle = if (show_sub) "Provisional model: adjusted OR (95% CI), reference = median" else NULL,
           x = "Odds ratio (log scale)", y = NULL) +
      theme_bw(base_size = FS_BASE) +
      theme(legend.position = if (show_legend) "top" else "none",
            plot.title    = element_text(size = FS_TITLE, face = "bold"),
            plot.subtitle = element_text(size = FS_SUB, colour = "grey30"),
            axis.text.y   = element_text(size = FS_AXIS),
            axis.text.x   = element_text(size = FS_AXIS - 1),
            axis.title.x  = element_text(size = FS_AXIST),
            legend.text   = element_text(size = FS_SUB),
            panel.grid.minor = element_blank())

    tab <- ggplot(d, aes(x = 0.5, y = ylab)) +
      geom_text(aes(label = or_text), hjust = 0.5, size = FS_CITEXT) +
      scale_x_continuous(limits = c(0, 1)) +
      coord_cartesian(clip = "off") +
      labs(title = "OR [95% CI]") +
      theme_void(base_size = FS_BASE) +
      theme(plot.title  = element_text(size = FS_AXIST, face = "bold", hjust = 0.5),
            plot.margin = margin(5.5, 5.5, 5.5, 0))

    forest + tab + patchwork::plot_layout(widths = c(3, 1.35))
  }

  done_pipes <- character(0)
  for (pp in PIPELINES) {
    d <- forest_dt[pipeline == pp]
    if (nrow(d) == 0) next
    g <- make_forest(copy(d), relabel_pipeline(pp))
    done_pipes <- c(done_pipes, pp)
    h <- max(2.8, 0.55 * nrow(d) + 1.8)
    ggsave(file.path(FIG_DIR, sprintf("fig4_forest_provisional_%s.png", pp)),
           g, width = 9.6, height = h, dpi = 150, limitsize = FALSE)
  }

  if (length(done_pipes) > 0) {
    panels <- lapply(done_pipes, function(pp)
      make_forest(copy(forest_dt[pipeline == pp]), relabel_pipeline(pp),
                  show_legend = FALSE, show_sub = FALSE))
    combo <- patchwork::wrap_plots(panels, ncol = 2) +
      patchwork::plot_annotation(
        title = "Main effect models - adjusted odds ratios (forest plots)",
        subtitle = paste0("OR (95% CI) annotated; reference = median; OR=1 dashed. ",
                          "Red = CI excludes 1. Arrows = CI beyond display limits ",
                          sprintf("[%g, %g].", FOREST_OR_LIM[1], FOREST_OR_LIM[2])),
        theme = theme(plot.title = element_text(size = 18, face = "bold"),
                      plot.subtitle = element_text(size = 12, colour = "grey30")))
    n_rows   <- ceiling(length(done_pipes) / 2)
    max_nrow <- max(vapply(done_pipes, function(pp) nrow(forest_dt[pipeline == pp]), integer(1)))
    ggsave(file.path(FIG_DIR, "fig4_forest_provisional_all.png"),
           combo, width = 19, height = n_rows * (0.5 * max_nrow + 2.4) + 0.8,
           dpi = 150, limitsize = FALSE)
  }

  message("\n暫定モデル forest plot を保存: ",
          normalizePath(FIG_DIR), "/fig4_forest_provisional_<pipeline>.png ",
          "（および ..._all.png）")
  trunc_n <- forest_dt[trunc_low | trunc_high, .N]
  if (trunc_n > 0)
    message(sprintf("  ※ CI が表示域 [%g, %g] を超え矢印で打ち切った項: %d（広い CI は Firth/分離由来の可能性）。",
                    FOREST_OR_LIM[1], FOREST_OR_LIM[2], trunc_n))
}

message("\n§11 完了。",
        "(1) 暫定+最終モデルの係数・CI 統合表（bnb_coef_table.csv/.rds）、",
        "(2) 暫定モデルの forest plot（fig4_forest_provisional_*.png）を出力。",
        " 最終モデルの OR カーブ/層別/±MDC 図は本スクリプトの対象外（ユーザー指示）。")
