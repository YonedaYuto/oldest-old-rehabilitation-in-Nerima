library(data.table)
library(here)
library(ggplot2)

need_pkg <- function(p) if (!requireNamespace(p, quietly = TRUE))
  stop(sprintf("パッケージ '%s' が必要です。install.packages('%s')", p, p))
need_pkg("ggcorrplot")
need_pkg("patchwork")

LABELS_PATH <- here::here("00_labels.R")
if (!exists("relabel_vars")) {
  if (!file.exists(LABELS_PATH)) stop("00_labels.R が見つかりません: ", LABELS_PATH)
  source(LABELS_PATH)
}

OUT_DIR <- here::here("data")
FIG_DIR <- here::here("figures")
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

COR_THRESH <- 0.80
COR_NOTE   <- 0.70

DROP_VARS <- list(
  severe_A  = character(0),
  severe_B  = character(0),
  elderly_A = character(0),
  elderly_B = character(0)
)

if (!exists("pipelines"))     pipelines     <- readRDS(file.path(OUT_DIR, "bnb_pipelines_mnar.rds"))
if (!exists("pipeline_meta")) pipeline_meta <- readRDS(file.path(OUT_DIR, "bnb_pipeline_meta.rds"))

center_continuous <- function(dt, meta) {
  out <- copy(dt)
  contv <- intersect(meta$continuous, names(out))
  centers <- numeric(0)
  for (v in contv) {
    med <- median(out[[v]], na.rm = TRUE)
    centers[v] <- med
    out[, paste0(v, "_c") := get(v) - med]
  }
  list(dt = out, centers = centers)
}

build_numeric_mat <- function(dt, meta, use_centered = TRUE, include_outcome = TRUE) {
  cols <- list()

  for (v in intersect(meta$continuous, names(dt))) {
    nm <- if (use_centered && paste0(v, "_c") %in% names(dt)) paste0(v, "_c") else v
    cols[[v]] <- as.numeric(dt[[nm]])
  }
  for (v in intersect(meta$categorical, names(dt))) {
    x <- dt[[v]]
    if (is.numeric(x) || is.integer(x)) {
      cols[[v]] <- as.numeric(x)
    } else {
      f <- as.factor(x)
      if (nlevels(f) <= 2L) {
        cols[[v]] <- as.numeric(f) - 1
      } else {
        for (lv in levels(f)[-1]) {
          cols[[paste0(v, ":", lv)]] <- as.numeric(f == lv)
        }
      }
    }
  }
  if (include_outcome && "good" %in% names(dt)) cols[["good"]] <- as.numeric(dt[["good"]])

  as.matrix(as.data.frame(cols, check.names = FALSE))
}

high_cor_pairs <- function(cmat, thresh) {
  cmat[lower.tri(cmat, diag = TRUE)] <- NA
  idx <- which(abs(cmat) >= thresh, arr.ind = TRUE)
  if (nrow(idx) == 0) return(data.table(var1 = character(0), var2 = character(0), r = numeric(0)))
  dt <- data.table(
    var1 = rownames(cmat)[idx[, 1]],
    var2 = colnames(cmat)[idx[, 2]],
    r    = round(cmat[idx], 3)
  )
  dt[order(-abs(r))]
}

make_corr_figure <- function(dt, meta, label) {
  disp <- relabel_pipeline(label)

  contv <- intersect(meta$continuous, names(dt))
  matc  <- sapply(contv, function(v) {
    nm <- if (paste0(v, "_c") %in% names(dt)) paste0(v, "_c") else v
    as.numeric(dt[[nm]])
  })
  colnames(matc) <- contv
  cor_p <- cor(matc, use = "pairwise.complete.obs", method = "pearson")
  dimnames(cor_p) <- list(relabel_vars(rownames(cor_p)), relabel_vars(colnames(cor_p)))

  matm  <- build_numeric_mat(dt, meta, use_centered = TRUE, include_outcome = TRUE)
  cor_s <- cor(matm, use = "pairwise.complete.obs", method = "spearman")
  dimnames(cor_s) <- list(relabel_vars(rownames(cor_s)), relabel_vars(colnames(cor_s)))

  pA <- suppressWarnings(ggcorrplot::ggcorrplot(
    cor_p, type = "lower", lab = TRUE, lab_size = 2.4,
    colors = c("#2166AC", "white", "#B2182B"), tl.cex = 8
  )) + ggtitle(sprintf("%s\nA: Continuous variables (Pearson)", disp)) +
    theme(plot.title = element_text(size = 10))

  pB <- suppressWarnings(ggcorrplot::ggcorrplot(
    cor_s, type = "lower", lab = FALSE,
    colors = c("#2166AC", "white", "#B2182B"), tl.cex = 7
  )) + ggtitle(sprintf("%s\nB: All variables incl. categorical (Spearman)", disp)) +
    theme(plot.title = element_text(size = 10))

  fig <- patchwork::wrap_plots(pA, pB, ncol = 2)
  list(fig = fig, cor_pearson = cor_p, cor_spearman = cor_s)
}

centered_pipelines <- list()
centers_list       <- list()
selected_meta      <- list()

for (nm in names(pipelines)) {
  dt0  <- pipelines[[nm]]
  meta <- pipeline_meta[[nm]]

  cc <- center_continuous(dt0, meta)
  centered_pipelines[[nm]] <- cc$dt
  centers_list[[nm]]       <- cc$centers

  cf <- make_corr_figure(cc$dt, meta, nm)
  ggsave(file.path(FIG_DIR, sprintf("fig1_corr_%s.png", nm)),
         cf$fig, width = 12, height = 6, dpi = 150)

  hp_p <- high_cor_pairs(cf$cor_pearson,  COR_NOTE)
  hp_s <- high_cor_pairs(cf$cor_spearman, COR_NOTE)
  message(sprintf("==== [%s] 高相関ペア (|r| >= %.2f) ====", nm, COR_NOTE))
  message("-- Pearson(連続) --");  print(hp_p)
  message("-- Spearman(全体) --"); print(hp_s)

  drop <- DROP_VARS[[nm]]
  cont_sel <- setdiff(meta$continuous,  drop)
  cat_sel  <- setdiff(meta$categorical, drop)
  if (length(drop) > 0)
    message(sprintf("[%s] 事前選択で除外: %s", nm, paste(drop, collapse = ", ")))

  m <- meta
  m$centers          <- cc$centers
  m$selected_cont    <- cont_sel
  m$selected_cat     <- cat_sel
  m$dropped_preselect <- drop
  selected_meta[[nm]] <- m
}

saveRDS(centered_pipelines, file.path(OUT_DIR, "bnb_pipelines_centered.rds"))
saveRDS(centers_list,       file.path(OUT_DIR, "bnb_centers.rds"))
saveRDS(selected_meta,      file.path(OUT_DIR, "bnb_selected_meta.rds"))

message("§3 完了。図1を ", normalizePath(FIG_DIR),
        " に保存。DROP_VARS を図1に基づき記入して再実行すると事前選択が反映されます。")
