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

DIAG_MODE <- "average"
REP_IMP   <- 1L

PEARSON_CUT <- 3
LEV_MULT    <- 4
COOK_NUM    <- 6
DFB_MULT    <- 5

COOK_CONTOUR_EXTRA <- numeric(0)

OUTLIER_IDS <- list(
  severe_A  = integer(0),
  severe_B  = integer(0),
  elderly_A = integer(0),
  elderly_B = integer(0)
)

if (!exists("est_index"))
  est_index <- readRDS(file.path(OUT_DIR, "bnb_estimation_index.rds"))

relabel_coef <- function(cn) {
  vapply(cn, function(s) {
    if (s == "(Intercept)") return("(Intercept)")
    if (grepl("_c[0-9]+$", s)) {
      base <- sub("_c[0-9]+$", "", s); j <- sub(".*_c", "", s)
      return(paste0(relabel_vars(base), " (nl", j, ")"))
    }
    if (grepl("_c$", s))          return(relabel_vars(s))
    if (grepl("_measurable", s))  return(relabel_vars(s))
    hit <- names(VAR_LABELS)[vapply(names(VAR_LABELS),
                                    function(v) startsWith(s, v), logical(1))]
    if (length(hit) > 0) {
      v   <- hit[which.max(nchar(hit))]
      lev <- sub(paste0("^", v), "", s)
      return(paste0(relabel_vars(v), ": ", relabel_levels(v, lev)))
    }
    s
  }, character(1), USE.NAMES = FALSE)
}

outliers_index <- data.table()

for (nm in est_index$pipeline) {
  est  <- readRDS(file.path(OUT_DIR, sprintf("bnb_estimation_%s.rds", nm)))
  obj  <- readRDS(file.path(OUT_DIR, sprintf("bnb_imp_%s.rds", nm)))
  imps <- obj$smcfcs$impDatasets
  smf  <- stats::as.formula(est$smformula)
  n    <- est$n
  m    <- length(imps)
  ids  <- seq_len(n)

  if (identical(est$path, "firth"))
    message(sprintf("[%s] §6 は Firth 経路。診断は glm による近似である点に留意。", nm))

  if (DIAG_MODE == "single") {
    fit <- glm(smf, family = binomial(), data = imps[[REP_IMP]])
    d   <- diag_one(fit)
    p   <- length(coef(fit))
    pearson_v <- d$pearson; leverage_v <- d$leverage; cook_v <- d$cook
    dfb_m <- d$dfb
    overlay_long <- NULL
    mode_note <- sprintf("Single imputation (#%d)", REP_IMP)

  } else {
    fits <- lapply(imps, function(dd) glm(smf, family = binomial(), data = dd))
    p    <- length(coef(fits[[1]]))
    PRS  <- vapply(fits, function(f) residuals(f, type = "pearson"), numeric(n))
    LEV  <- vapply(fits, hatvalues,      numeric(n))
    CKS  <- vapply(fits, cooks.distance, numeric(n))
    dfbl <- lapply(fits, function(f) as.matrix(dfbetas(f)))
    cn   <- colnames(dfbl[[1]])
    DFB  <- array(unlist(dfbl), dim = c(n, length(cn), m))

    pearson_v  <- rowMeans(PRS); leverage_v <- rowMeans(LEV); cook_v <- rowMeans(CKS)
    dfb_m      <- apply(DFB, c(1, 2), mean); colnames(dfb_m) <- cn

    mode_note <- sprintf("Averaged over m=%d imputations (per .id)", m)
    overlay_long <- if (DIAG_MODE == "overlay") {
      rbindlist(lapply(seq_len(m), function(k) data.table(
        .id = ids, .imp = k,
        Pearson = PRS[, k], Leverage = LEV[, k], Cook = CKS[, k])))
    } else NULL
  }

  lev_cut  <- LEV_MULT * p / n
  cook_cut <- COOK_NUM / n
  dfb_cut  <- DFB_MULT / sqrt(n)

  dfb_max      <- apply(abs(dfb_m), 1, max)
  dfb_exceed   <- abs(dfb_m) > dfb_cut
  n_dfb_exceed <- rowSums(dfb_exceed)
  flag_pearson <- abs(pearson_v) > PEARSON_CUT
  flag_lev     <- leverage_v     > lev_cut
  flag_cook    <- cook_v         > cook_cut
  flag_dfb     <- dfb_max        > dfb_cut
  flag_any     <- flag_pearson | flag_lev | flag_cook | flag_dfb

  diag_dt <- data.table(
    .id = ids,
    pearson = pearson_v, leverage = leverage_v, cook = cook_v, dfb_max = dfb_max,
    n_dfb_exceed = n_dfb_exceed,
    flag_pearson, flag_lev, flag_cook, flag_dfb, flag_any)

  flag_A <- flag_pearson | flag_lev | flag_cook
  panelA <- data.table(.id = ids, pearson = pearson_v, leverage = leverage_v,
                       flag = flag_A)

  x_pad <- 0.5
  xlim  <- range(c(panelA$pearson, -PEARSON_CUT - x_pad, PEARSON_CUT + x_pad))
  ylim  <- c(0, max(panelA$leverage, lev_cut) * 1.10)

  cook_levels <- c(cook_cut, COOK_CONTOUR_EXTRA)
  cook_dt <- make_cook_contour(cook_levels, p = p, hmax = ylim[2])

  dfb_dt <- as.data.table(dfb_m); dfb_dt[, .id := ids]
  panelB <- melt(dfb_dt, id.vars = ".id", variable.name = "coef", value.name = "value")
  panelB[, coef := factor(relabel_coef(as.character(coef)))]
  panelB[, flag := abs(value) > dfb_cut]

  sub <- sprintf(paste0("%s | n=%d, p=%d | A: dotted purple = Cook contour(s) {%s}; ",
                        "dashed = |Pearson|=%g (vert) & hat=3p/n=%.3f (horiz). ",
                        "B: |dfbetas|=2/sqrt(n)=%.3f"),
                 mode_note, n, p,
                 paste(sprintf("%.4g", cook_levels), collapse = ", "),
                 PEARSON_CUT, lev_cut, dfb_cut)

  fig <- make_influence_figure(panelA, cook_dt, panelB,
                               pearson_cut = PEARSON_CUT, lev_cut = lev_cut,
                               dfb_cut = dfb_cut, xlim = xlim, ylim = ylim,
                               label = nm, sub = sub)
  n_coef <- length(unique(panelB$coef)); ncolB <- min(4L, n_coef)
  ggsave(file.path(FIG_DIR, sprintf("fig3_influence_%s.png", nm)),
         fig, width = 6 + 1.8 * ncolB,
         height = max(6, 1.6 * ceiling(n_coef / ncolB) + 1.8), dpi = 150)

  message(sprintf("==== [%s] 外れ値診断 (%s) | n=%d, p=%d ====", nm, mode_note, n, p))
  message(sprintf("  閾値  |Pearson|>%g, hat>%.4f(=3p/n), Cook>%.5f(=4/n), |dfbetas|>%.4f(=2/√n)",
                  PEARSON_CUT, lev_cut, cook_cut, dfb_cut))
  message(sprintf("  超過数 Pearson=%d, Leverage=%d, Cook=%d, dfbetas=%d, いずれか=%d",
                  sum(flag_pearson), sum(flag_lev), sum(flag_cook), sum(flag_dfb), sum(flag_any)))
  if (sum(flag_any) > 0) {
    flagged <- diag_dt[flag_any == TRUE][order(-pmax(abs(pearson), cook, dfb_max))]
    message("  閾値を超えた観測（上位20、影響度順）:")
    print(head(flagged[, .(.id, pearson = round(pearson, 2), leverage = round(leverage, 3),
                           cook = round(cook, 4), dfb_max = round(dfb_max, 3),
                           flag_pearson, flag_lev, flag_cook, flag_dfb)], 20))
  } else {
    message("  → 既定閾値を超える観測なし。")
  }

  excl <- OUTLIER_IDS[[nm]]
  excl <- excl[excl %in% ids]
  if (length(excl) > 0)
    message(sprintf("  除外指定 .id（パイプライン共通で全代入から除去）: %s",
                    paste(sort(excl), collapse = ", ")))
  else
    message("  除外指定なし（図3を見て OUTLIER_IDS に記入し再実行すると §8 へ反映）。")

  saveRDS(list(
    pipeline = nm, n = n, p = p, m = m, diag_mode = DIAG_MODE,
    thresholds = list(pearson = PEARSON_CUT, leverage = lev_cut,
                      cook = cook_cut, dfbetas = dfb_cut),
    cook_contour_levels = cook_levels,
    diagnostics = diag_dt, overlay = overlay_long,
    dfbetas = cbind(data.table(.id = ids), as.data.table(dfb_m)),
    exclude_ids = sort(unique(excl)),
    smformula = est$smformula, path = est$path
  ), file.path(OUT_DIR, sprintf("bnb_influence_%s.rds", nm)))

  outliers_index <- rbind(outliers_index, data.table(
    pipeline = nm, n = n, p = p,
    cut_lev = round(lev_cut, 4), cut_cook = round(cook_cut, 5), cut_dfb = round(dfb_cut, 4),
    n_flag_pearson = sum(flag_pearson), n_flag_lev = sum(flag_lev),
    n_flag_cook = sum(flag_cook), n_flag_dfb = sum(flag_dfb),
    n_flag_any = sum(flag_any), n_excluded = length(excl)))
}

saveRDS(outliers_index, file.path(OUT_DIR, "bnb_outliers_index.rds"))

message("\n==== §7 サマリ（閾値超過・除外）====")
print(outliers_index)
message("§7 完了。図3を ", normalizePath(FIG_DIR),
        " に保存。パネルAは横軸=ピアソン残差・縦軸=てこ比、Cook距離は等高線で重畳。",
        " 図3を見て OUTLIER_IDS（.id）を記入し再実行すると、除外 .id が ",
        "bnb_influence_<pipeline>.rds に保存され、§8 が全代入から共通除去します。")
