library(data.table)
library(here)

OUT_DIR <- here::here("data")

LABELS_PATH <- here::here("00_labels.R")
if (!exists("relabel_vars")) {
  if (!file.exists(LABELS_PATH)) stop("00_labels.R が見つかりません: ", LABELS_PATH)
  source(LABELS_PATH)
}

DFB_MIN_VARS      <- 2L
COUNT_LEVEL       <- "coef"
INCLUDE_INTERCEPT <- FALSE

if (!exists("est_index"))
  est_index <- readRDS(file.path(OUT_DIR, "bnb_estimation_index.rds"))

coef_to_var <- function(cn) {
  vapply(cn, function(s) {
    if (s == "(Intercept)")        return("(Intercept)")
    if (grepl("_measurable$", s))  return(s)
    if (grepl("_c[0-9]+$", s))     return(sub("_c[0-9]+$", "", s))
    if (grepl("_c$", s))           return(sub("_c$", "", s))
    hit <- names(VAR_LABELS)[vapply(names(VAR_LABELS),
                                    function(v) startsWith(s, v), logical(1))]
    if (length(hit) > 0)           return(hit[which.max(nchar(hit))])
    s
  }, character(1), USE.NAMES = FALSE)
}

recompute_dfbetas <- function(nm, diag_mode) {
  est  <- readRDS(file.path(OUT_DIR, sprintf("bnb_estimation_%s.rds", nm)))
  obj  <- readRDS(file.path(OUT_DIR, sprintf("bnb_imp_%s.rds", nm)))
  imps <- obj$smcfcs$impDatasets
  smf  <- stats::as.formula(est$smformula)
  n    <- est$n; m <- length(imps)
  if (identical(diag_mode, "single")) {
    fit <- glm(smf, family = binomial(), data = imps[[1]])
    dfb_m <- as.matrix(dfbetas(fit))
  } else {
    dfbl <- lapply(imps, function(dd) as.matrix(dfbetas(glm(smf, binomial(), data = dd))))
    cn   <- colnames(dfbl[[1]])
    DFB  <- array(unlist(dfbl), dim = c(n, length(cn), m))
    dfb_m <- apply(DFB, c(1, 2), mean); colnames(dfb_m) <- cn
  }
  cbind(data.table(.id = seq_len(n)), as.data.table(dfb_m))
}

def_index <- data.table()

for (nm in est_index$pipeline) {
  inf  <- readRDS(file.path(OUT_DIR, sprintf("bnb_influence_%s.rds", nm)))
  thr  <- inf$thresholds
  diag <- as.data.table(inf$diagnostics)
  n <- inf$n; p <- inf$p
  ids <- diag$.id

  crit1 <- diag$flag_pearson | diag$flag_lev | diag$flag_cook

  if (COUNT_LEVEL == "coef" && "n_dfb_exceed" %in% names(diag) && !INCLUDE_INTERCEPT) {
    exceed_count <- diag$n_dfb_exceed
    count_unit <- "coefficients"
  } else {
    dfb_dt <- recompute_dfbetas(nm, inf$diag_mode)
    setkey(dfb_dt, .id); setkey(diag, .id)
    stopifnot(identical(diag$.id, dfb_dt$.id))

    coefcols <- setdiff(names(dfb_dt), ".id")
    if (!INCLUDE_INTERCEPT) coefcols <- setdiff(coefcols, "(Intercept)")
    exceed <- abs(as.matrix(dfb_dt[, ..coefcols])) > thr$dfbetas

    if (COUNT_LEVEL == "variable") {
      var_of <- coef_to_var(coefcols)
      uv <- unique(var_of)
      exceed_var <- vapply(uv, function(v) {
        cols <- which(var_of == v)
        rowSums(exceed[, cols, drop = FALSE]) > 0
      }, logical(nrow(exceed)))
      exceed_count <- rowSums(exceed_var)
      count_unit <- "variables"
    } else {
      exceed_count <- rowSums(exceed)
      count_unit <- "coefficients"
    }
  }
  crit2 <- exceed_count >= DFB_MIN_VARS

  outlier <- crit1 & crit2

  def_dt <- data.table(
    .id = ids,
    pearson = round(diag$pearson, 3), leverage = round(diag$leverage, 4),
    cook = round(diag$cook, 5), dfb_exceed = exceed_count,
    crit1_resid_lev_cook = crit1, crit2_dfbetas_2plus = crit2, outlier = outlier)

  outlier_ids <- sort(ids[outlier])
  crit1_ids   <- sort(ids[crit1])
  crit2_ids   <- sort(ids[crit2])

  message(sprintf("==== [%s] 外れ値定義 (n=%d, p=%d) | dfbetas は%s単位で数える ====",
                  nm, n, p, count_unit))
  message(sprintf("  閾値（より大きい）: |Pearson|>%g, hat>%.4f, Cook>%.5f, |dfbetas|>%.4f (>=%d %s)",
                  thr$pearson, thr$leverage, thr$cook, thr$dfbetas, DFB_MIN_VARS, count_unit))
  message(sprintf("  基準1(残差/てこ比/Cook いずれか): %d 件", length(crit1_ids)))
  message(sprintf("  基準2(dfbetas %d%s以上)         : %d 件", DFB_MIN_VARS, count_unit, length(crit2_ids)))
  message(sprintf("  両基準に該当                    : %d 件", sum(crit1 & crit2)))
  message(sprintf("  外れ値(基準1 ∩ 基準2)           : %d 件 (%.1f%%)",
                  length(outlier_ids), 100 * length(outlier_ids) / n))
  if (length(outlier_ids) > 0) {
    message("  外れ値 .id: ", paste(outlier_ids, collapse = ", "))
    print(def_dt[outlier == TRUE][order(-pmax(abs(pearson), cook), -dfb_exceed)])
  } else {
    message("  → 外れ値なし。")
  }

  saveRDS(list(
    pipeline = nm, n = n, p = p,
    thresholds = thr,
    rule = list(crit1 = "any of |Pearson|/hat/Cook > threshold",
                crit2 = sprintf("dfbetas > threshold in >= %d %s", DFB_MIN_VARS, count_unit),
                combine = "AND  (outlier = crit1 & crit2)",
                count_level = COUNT_LEVEL, dfb_min_vars = DFB_MIN_VARS,
                include_intercept = INCLUDE_INTERCEPT, comparison = ">"),
    table = def_dt,
    outlier_ids = outlier_ids, crit1_ids = crit1_ids, crit2_ids = crit2_ids
  ), file.path(OUT_DIR, sprintf("bnb_outlier_def_%s.rds", nm)))

  def_index <- rbind(def_index, data.table(
    pipeline = nm, n = n, p = p, count_unit = count_unit,
    n_crit1 = length(crit1_ids), n_crit2 = length(crit2_ids),
    n_both = sum(crit1 & crit2), n_outlier = length(outlier_ids),
    pct_outlier = round(100 * length(outlier_ids) / n, 1)))
}

saveRDS(def_index, file.path(OUT_DIR, "bnb_outlier_def_index.rds"))

message("\n==== §7-1 サマリ（外れ値定義）====")
print(def_index)
message("§7-1 完了。外れ値 .id を bnb_outlier_def_<pipeline>.rds に保存しました。",
        " §8 はこの outlier_ids を全代入データから共通除去して変数選択に進みます。")
