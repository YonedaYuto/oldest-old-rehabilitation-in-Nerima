library(data.table)
library(here)
library(ggplot2)

LABELS_PATH <- here::here("00_labels.R")
if (!exists("relabel_vars")) {
  if (!file.exists(LABELS_PATH)) stop("00_labels.R が見つかりません: ", LABELS_PATH)
  source(LABELS_PATH)
}

OUT_DIR <- here::here("data")
FIG_DIR <- here::here("figures")
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

PIPELINES <- names(PIPE_LABELS)
MODELS    <- c("Provisional", "Final")

GVIF_ADJ_HIGH <- sqrt(10)
GVIF_ADJ_MOD  <- sqrt(5)

COR_N_IMP   <- 5L
COR_IMP_IDS <- NULL

TOP_PAIRS   <- 20L
TOP_DRIVERS <- 5L
DIGITS      <- 3L

read_safe <- function(f) if (file.exists(f)) readRDS(f) else NULL

pipe_disp <- function(x) vapply(as.character(x), relabel_pipeline, character(1))

reduce_imps <- function(obj, n0, outlier_ids) {
  imps_full <- obj$smcfcs$impDatasets
  bad <- sort(unique(outlier_ids)); bad <- bad[bad %in% seq_len(n0)]
  if (length(bad) > 0) lapply(imps_full, function(d) d[-bad, , drop = FALSE]) else imps_full
}
pick_imps <- function(m, k, ids = NULL) {
  if (!is.null(ids)) return(sort(unique(ids[ids >= 1 & ids <= m])))
  if (is.null(k) || k >= m) return(seq_len(m))
  unique(round(seq(1, m, length.out = k)))
}

col_to_varkey <- function(s) {
  s <- as.character(s)
  if (s == "(Intercept)") return("(Intercept)")
  if (grepl(":", s, fixed = TRUE)) {
    parts <- strsplit(s, ":", fixed = TRUE)[[1]]
    return(paste(vapply(parts, col_to_varkey, character(1)), collapse = ":"))
  }
  if (grepl("_measurable", s))  return(s)
  if (grepl("_c[0-9]+$", s))    return(sub("_c[0-9]+$", "", s))
  if (grepl("_c$", s))          return(sub("_c$", "", s))
  hit <- names(VAR_LABELS)[vapply(names(VAR_LABELS),
                                  function(v) startsWith(s, v), logical(1))]
  if (length(hit) > 0) return(hit[which.max(nchar(hit))])
  s
}

varkey_disp <- function(k) {
  vapply(k, function(s) {
    if (grepl(":", s, fixed = TRUE)) {
      parts <- strsplit(s, ":", fixed = TRUE)[[1]]
      return(paste(vapply(parts, function(p) relabel_vars(p), character(1)), collapse = " x "))
    }
    relabel_vars(s)
  }, character(1), USE.NAMES = FALSE)
}

aux_r2 <- function(R, j, S) {
  S <- setdiff(S, j); if (length(S) == 0) return(0)
  RSS <- R[S, S, drop = FALSE]
  inv <- tryCatch(solve(RSS), error = function(e) MASS::ginv(RSS))
  rjS <- R[j, S, drop = TRUE]
  val <- as.numeric(crossprod(rjS, inv %*% rjS))
  min(max(val, 0), 1 - 1e-12)
}

gvif_det <- function(R, idx) {
  others <- setdiff(seq_len(nrow(R)), idx)
  dv <- det(R[idx, idx, drop = FALSE])
  do <- if (length(others) > 0) det(R[others, others, drop = FALSE]) else 1
  da <- det(R)
  if (!is.finite(da) || abs(da) < 1e-12) return(NA_real_)
  (dv * do) / da
}

design_cor <- function(smf, d) {
  mm <- tryCatch(stats::model.matrix(stats::delete.response(stats::terms(smf)), data = d),
                 error = function(e) NULL)
  if (is.null(mm)) return(NULL)
  keep <- setdiff(colnames(mm), "(Intercept)")
  mm <- mm[, keep, drop = FALSE]
  sds <- apply(mm, 2, stats::sd)
  ok  <- is.finite(sds) & sds > 0
  mm  <- mm[, ok, drop = FALSE]
  if (ncol(mm) < 2) return(NULL)
  list(R = stats::cor(mm), cols = colnames(mm))
}

analyze_model <- function(nm, md, smf, imps, ids) {
  cors <- lapply(ids, function(i) design_cor(smf, imps[[i]]))
  cors <- cors[!vapply(cors, is.null, logical(1))]
  if (length(cors) == 0) return(NULL)
  common <- Reduce(intersect, lapply(cors, `[[`, "cols"))
  if (length(common) < 2) return(NULL)
  Rsum <- matrix(0, length(common), length(common), dimnames = list(common, common))
  for (cc in cors) Rsum <- Rsum + cc$R[common, common]
  R <- Rsum / length(cors)
  diag(R) <- 1

  cols   <- common
  varkey <- vapply(cols, col_to_varkey, character(1))
  vars   <- unique(varkey)
  col_idx_by_var <- split(seq_along(cols), varkey)[vars]

  Rinv <- tryCatch(solve(R), error = function(e) MASS::ginv(R))
  col_vif <- if (is.matrix(Rinv) && all(dim(Rinv) == dim(R))) diag(Rinv) else rep(NA_real_, length(cols))
  names(col_vif) <- cols

  gvif_dt <- rbindlist(lapply(vars, function(v) {
    idx <- col_idx_by_var[[v]]
    g <- gvif_det(R, idx); dfv <- length(idx)
    data.table(pipeline = nm, model = md, var = v, variable = varkey_disp(v),
               Df = dfv, GVIF = g,
               GVIF_adj = if (is.finite(g) && g > 0) g^(1 / (2 * dfv)) else NA_real_,
               max_col_VIF = max(col_vif[idx], na.rm = TRUE))
  }))
  gvif_dt[, level := fifelse(is.finite(GVIF_adj) & GVIF_adj > GVIF_ADJ_HIGH, "high",
                             fifelse(is.finite(GVIF_adj) & GVIF_adj > GVIF_ADJ_MOD, "moderate", "ok"))]
  setorder(gvif_dt, -GVIF_adj)

  pair_list <- list()
  for (a in seq_along(vars)) for (b in seq_len(a - 1)) {
    ia <- col_idx_by_var[[vars[a]]]; ib <- col_idx_by_var[[vars[b]]]
    sub <- abs(R[ia, ib, drop = FALSE])
    wm  <- which(sub == max(sub), arr.ind = TRUE)[1, ]
    pair_list[[length(pair_list) + 1L]] <- data.table(
      pipeline = nm, model = md,
      var_a = vars[a], var_b = vars[b],
      variable_a = varkey_disp(vars[a]), variable_b = varkey_disp(vars[b]),
      max_abs_cor = max(sub),
      col_a = cols[ia][wm[1]], col_b = cols[ib][wm[2]],
      signed_cor = R[ia[wm[1]], ib[wm[2]]])
  }
  pairs_dt <- if (length(pair_list)) rbindlist(pair_list) else data.table()
  if (nrow(pairs_dt)) setorder(pairs_dt, -max_abs_cor)

  high_vars <- gvif_dt[level == "high", var]
  attr_list <- list()
  for (v in high_vars) {
    idx <- col_idx_by_var[[v]]
    jstar <- idx[which.max(col_vif[idx])]
    O <- setdiff(seq_along(cols), idx)
    r2_full <- aux_r2(R, jstar, O)
    other_vars <- setdiff(vars, v)
    contrib <- vapply(other_vars, function(k) {
      idk <- col_idx_by_var[[k]]
      r2_full - aux_r2(R, jstar, setdiff(O, idk))
    }, numeric(1))
    maxcor_k <- vapply(other_vars, function(k) {
      idk <- col_idx_by_var[[k]]; max(abs(R[jstar, idk]))
    }, numeric(1))
    ad <- data.table(pipeline = nm, model = md,
                     high_var = v, high_variable = varkey_disp(v),
                     GVIF_adj = gvif_dt[var == v, GVIF_adj][1],
                     rep_col = cols[jstar], aux_R2 = r2_full,
                     driver_var = other_vars, driver_variable = varkey_disp(other_vars),
                     contrib_R2 = contrib, max_abs_cor = maxcor_k)
    setorder(ad, -contrib_R2)
    attr_list[[v]] <- ad
  }
  attr_dt <- if (length(attr_list)) rbindlist(attr_list) else data.table()

  list(R = R, cols = cols, varkey = varkey, vars = vars,
       gvif = gvif_dt, pairs = pairs_dt, attribution = attr_dt)
}

plot_corr_heatmap <- function(res, title) {
  vars <- res$vars
  A <- matrix(NA_real_, length(vars), length(vars), dimnames = list(vars, vars))
  ci <- split(seq_along(res$cols), res$varkey)[vars]
  for (a in seq_along(vars)) for (b in seq_along(vars)) {
    ia <- ci[[vars[a]]]; ib <- ci[[vars[b]]]
    A[a, b] <- if (a == b) 1 else max(abs(res$R[ia, ib, drop = FALSE]))
  }
  vlab <- varkey_disp(vars)
  hv <- res$gvif[level == "high", var]
  mod <- res$gvif[level == "moderate", var]
  mark <- ifelse(vars %in% hv, "**", ifelse(vars %in% mod, "*", ""))
  vlab2 <- paste0(vlab, ifelse(mark == "", "", paste0(" ", mark)))

  dt <- as.data.table(as.table(A)); setnames(dt, c("a", "b", "absr"))
  dt[, a := factor(a, levels = vars, labels = vlab2)]
  dt[, b := factor(b, levels = vars, labels = vlab2)]
  ggplot(dt, aes(a, b, fill = absr)) +
    geom_tile(colour = "grey90") +
    geom_text(aes(label = ifelse(is.finite(absr), formatC(absr, format = "f", digits = 2), "")),
              size = 2.6, colour = "grey15") +
    scale_fill_gradient2(low = "white", mid = "#FDDBC7", high = "#B2182B",
                         midpoint = 0.5, limits = c(0, 1), name = "max |r|") +
    labs(title = title,
         subtitle = "Variable-block max |correlation| of design columns. ** = high GVIF, * = moderate.",
         x = NULL, y = NULL) +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 8, colour = "grey30"),
          panel.grid = element_blank())
}

plot_gvif_bar <- function(gvif_dt, title) {
  d <- copy(gvif_dt)[is.finite(GVIF_adj)]
  if (nrow(d) == 0) return(NULL)
  setorder(d, GVIF_adj)
  d[, ylab := factor(variable, levels = variable)]
  ggplot(d, aes(GVIF_adj, ylab, fill = level)) +
    geom_col(width = 0.7) +
    geom_vline(xintercept = GVIF_ADJ_MOD,  linetype = "dotted", colour = "grey40") +
    geom_vline(xintercept = GVIF_ADJ_HIGH, linetype = "dashed", colour = "#B2182B") +
    scale_fill_manual(values = c(ok = "grey60", moderate = "#F4A582", high = "#B2182B"),
                      name = NULL, drop = FALSE) +
    labs(title = title,
         subtitle = sprintf("GVIF^(1/(2*Df)); dotted = moderate (VIF~%.0f), dashed = high (VIF~%.0f).",
                            GVIF_ADJ_MOD^2, GVIF_ADJ_HIGH^2),
         x = "GVIF^(1/(2*Df))", y = NULL) +
    theme_bw(base_size = 9) +
    theme(legend.position = "top",
          plot.title = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 8, colour = "grey30"),
          panel.grid.minor = element_blank())
}

gvif_all <- list(); pairs_all <- list(); attr_all <- list()

for (nm in PIPELINES) {
  prov <- read_safe(file.path(OUT_DIR, sprintf("bnb_provisional_%s.rds", nm)))
  fin  <- read_safe(file.path(OUT_DIR, sprintf("bnb_final_%s.rds",       nm)))
  obj  <- read_safe(file.path(OUT_DIR, sprintf("bnb_imp_%s.rds",         nm)))
  if (is.null(obj)) { message(sprintf("[%s] 代入データなし。スキップ。", nm)); next }

  defs <- list()
  if (!is.null(prov)) defs$Provisional <-
    list(smf = stats::as.formula(prov$smformula_prov), n0 = prov$n0, outlier_ids = prov$outlier_ids)
  if (!is.null(fin))  defs$Final <-
    list(smf = stats::as.formula(fin$smformula_final), n0 = fin$n0, outlier_ids = fin$outlier_ids)
  if (length(defs) == 0) { message(sprintf("[%s] 暫定/最終モデルなし。スキップ。", nm)); next }

  message(sprintf("\n================ [%s] §12-1 GVIF 共線性分解 ================", nm))

  for (md in names(defs)) {
    de   <- defs[[md]]
    imps <- reduce_imps(obj, de$n0, de$outlier_ids)
    m    <- length(imps)
    ids  <- pick_imps(m, COR_N_IMP, COR_IMP_IDS)
    tl   <- attr(stats::terms(de$smf), "term.labels")
    if (length(tl) < 2) { message(sprintf("  [%s] 予測子 <2 → 共線性なし。スキップ。", md)); next }

    res <- analyze_model(nm, md, de$smf, imps, ids)
    if (is.null(res)) { message(sprintf("  [%s] 設計行列の相関を作れずスキップ。", md)); next }

    gvif_all[[paste0(nm, "_", md)]]  <- res$gvif
    if (nrow(res$pairs))       pairs_all[[paste0(nm, "_", md)]] <- res$pairs
    if (nrow(res$attribution)) attr_all[[paste0(nm, "_", md)]]  <- res$attribution

    n_high <- res$gvif[level == "high", .N]
    n_mod  <- res$gvif[level == "moderate", .N]
    message(sprintf("  [%s] 変数 %d / 設計列 %d（代入 %d 本平均）。高GVIF=%d, 中=%d",
                    md, length(res$vars), length(res$cols), length(ids), n_high, n_mod))

    g_show <- copy(res$gvif)[, .(Variable = variable, Df,
                                 GVIF = round(GVIF, DIGITS),
                                 `GVIF^(1/2Df)` = round(GVIF_adj, DIGITS),
                                 level)]
    print(head(g_show, 10), row.names = FALSE)

    if (n_high > 0) {
      message(sprintf("  → 高GVIF変数の共線性分解（代表列の補助R²を他変数別の低下量で分解, 上位%d）:", TOP_DRIVERS))
      for (v in res$gvif[level == "high", var]) {
        sub <- res$attribution[high_var == v]
        setorder(sub, -contrib_R2)
        head_sub <- head(sub, TOP_DRIVERS)
        message(sprintf("    ■ %s  (GVIF^(1/2Df)=%.2f, 代表列の補助R²=%.3f)",
                        varkey_disp(v), res$gvif[var == v, GVIF_adj][1],
                        sub$aux_R2[1]))
        print(head_sub[, .(driver = driver_variable,
                           contrib_R2 = round(contrib_R2, DIGITS),
                           max_abs_cor = round(max_abs_cor, DIGITS))], row.names = FALSE)
      }
    }

    if (n_high > 0 || n_mod > 0) {
      ttl <- sprintf("%s - %s model", relabel_pipeline(nm), md)
      gh <- plot_corr_heatmap(res, paste0(ttl, ": collinearity"))
      sz <- max(4.5, 0.42 * length(res$vars) + 2)
      ggsave(file.path(FIG_DIR, sprintf("fig_gvif_corr_%s_%s.png", nm, tolower(md))),
             gh, width = sz, height = sz - 0.4, dpi = 150)
      gb <- plot_gvif_bar(res$gvif, paste0(ttl, ": GVIF"))
      if (!is.null(gb))
        ggsave(file.path(FIG_DIR, sprintf("fig_gvif_bar_%s_%s.png", nm, tolower(md))),
               gb, width = 6.2, height = max(2.6, 0.34 * nrow(res$gvif) + 1.2), dpi = 150)
    }
  }
}

if (length(gvif_all)) {
  gvif_dt <- rbindlist(gvif_all, fill = TRUE)
  gvif_dt[, pipeline_disp := vapply(as.character(pipeline), relabel_pipeline, character(1))]
  saveRDS(gvif_dt, file.path(OUT_DIR, "bnb_gvif_table.rds"))
  fwrite(gvif_dt, file.path(OUT_DIR, "bnb_gvif_table.csv"))
} else gvif_dt <- data.table()

if (length(pairs_all)) {
  pairs_dt <- rbindlist(pairs_all, fill = TRUE)
  pairs_dt[, pipeline_disp := vapply(as.character(pipeline), relabel_pipeline, character(1))]
  saveRDS(pairs_dt, file.path(OUT_DIR, "bnb_collinearity_pairs.rds"))
  fwrite(pairs_dt, file.path(OUT_DIR, "bnb_collinearity_pairs.csv"))
} else pairs_dt <- data.table()

if (length(attr_all)) {
  attr_dt <- rbindlist(attr_all, fill = TRUE)
  attr_dt[, pipeline_disp := vapply(as.character(pipeline), relabel_pipeline, character(1))]
  saveRDS(attr_dt, file.path(OUT_DIR, "bnb_gvif_attribution.rds"))
  fwrite(attr_dt, file.path(OUT_DIR, "bnb_gvif_attribution.csv"))
} else attr_dt <- data.table()

message("\n==== §12-1 高 GVIF 変数の一覧（全モデル）====")
if (nrow(gvif_dt)) {
  hi <- gvif_dt[level %in% c("high", "moderate")]
  if (nrow(hi)) {
    setorder(hi, pipeline, model, -GVIF_adj)
    print(hi[, .(Pipeline = pipeline_disp, Model = model, Variable = variable, Df,
                 `GVIF^(1/2Df)` = round(GVIF_adj, DIGITS), level)], row.names = FALSE)
  } else message("  高/中程度の GVIF を持つ変数はありませんでした。")
}

message("\n==== §12-1 高 GVIF を押し上げている変数対（寄与上位）====")
if (nrow(attr_dt)) {
  topdrv <- attr_dt[, .SD[order(-contrib_R2)][1:min(.N, TOP_DRIVERS)],
                    by = .(pipeline, model, high_var)]
  print(topdrv[, .(Pipeline = vapply(as.character(pipeline), relabel_pipeline, character(1)),
                   Model = model,
                   `High-GVIF var` = high_variable, `Driven by` = driver_variable,
                   `contrib R2` = round(contrib_R2, DIGITS),
                   `max|r|` = round(max_abs_cor, DIGITS))], row.names = FALSE)
} else message("  高 GVIF の変数がないため、寄与分解はありません。")

message("\n§12-1 完了。")
message("  ・GVIF 表（高/中フラグ）       : data/bnb_gvif_table.csv")
message("  ・共線性ペア（最大|相関|, 降順）: data/bnb_collinearity_pairs.csv")
message("  ・高GVIF変数 -> 寄与する他変数  : data/bnb_gvif_attribution.csv")
message("  ・図                           : figures/fig_gvif_corr_<pipeline>_<model>.png, fig_gvif_bar_*.png")
message("注記: GVIF は『変数 vs 残り全変数』の指標でペアワイズではない。本スクリプトは設計行列の相関 R から")
message("      GVIF を det 公式で再計算（§12/car と整合）し、(a) 変数ブロック間の最大|相関|、")
message("      (b) 補助R²(=1-1/VIF)を他変数別の低下量に分解した『一意寄与 contrib_R2』で、どの変数が")
message("      高GVIFを押し上げているかを特定する。寄与は最後に抜く方式のため順序非依存の一意成分。")
message("      因子・RCS は変数ブロック（複数設計列）として集約。交互作用項は 1 変数として扱う。")
message("      相関は外れ値除外後の縮小代入データ ", ifelse(is.null(COR_N_IMP), "全", as.character(COR_N_IMP)),
        " 本で平均（§12 の VIF と同様の代入横断平均）。")
