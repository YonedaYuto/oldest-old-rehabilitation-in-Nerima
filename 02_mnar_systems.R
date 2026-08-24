library(data.table)
library(here)

OUT_DIR <- here::here("data")

if (!exists("df_severe"))  df_severe  <- readRDS(file.path(OUT_DIR, "bnb_preprocessed_severe.rds"))
if (!exists("df_elderly")) df_elderly <- readRDS(file.path(OUT_DIR, "bnb_preprocessed_elderly.rds"))
if (!exists("var_meta"))   var_meta   <- readRDS(file.path(OUT_DIR, "bnb_var_meta.rds"))

MNAR_VARS <- var_meta$mnar_vars

apply_mnar_system <- function(dt, system = c("A", "B"),
                              worst = c("zero", "obs_min"),
                              fillB = c("median", "zero"),
                              mnar_vars = MNAR_VARS) {
  system <- match.arg(system)
  worst  <- match.arg(worst)
  fillB  <- match.arg(fillB)
  out <- copy(dt)

  ind_cols <- character(0)

  for (v in mnar_vars) {
    if (!v %in% names(out)) stop(sprintf("MNAR変数 %s が見つかりません。", v))
    n_miss <- out[is.na(get(v)), .N]

    if (system == "A") {
      floor_val <- if (worst == "zero") 0 else min(out[[v]], na.rm = TRUE)
      out[is.na(get(v)), (v) := floor_val]

    } else {
      ind <- paste0(v, "_measurable")
      out[, (ind) := as.integer(!is.na(get(v)))]
      ind_cols <- c(ind_cols, ind)

      const_val <- if (fillB == "median") median(out[[v]], na.rm = TRUE) else 0
      out[is.na(get(v)), (v) := const_val]

      if (out[, uniqueN(get(ind))] < 2L) {
        warning(sprintf("指標 %s が無分散です（欠損 %d 件）。モデルから除外を検討してください。",
                        ind, n_miss))
      }
    }
  }

  remaining <- out[, vapply(.SD, function(x) sum(is.na(x)), integer(1)),
                    .SDcols = mnar_vars]
  if (any(remaining > 0)) {
    stop("系統 ", system, " 変換後も MNAR 変数に欠損が残っています: ",
         paste(names(remaining)[remaining > 0], collapse = ", "))
  }

  setattr(out, "mnar_system",    system)
  setattr(out, "mnar_value",     mnar_vars)
  setattr(out, "mnar_indicator", ind_cols)
  out[]
}

df_severe_A  <- apply_mnar_system(df_severe,  system = "A")
df_severe_B  <- apply_mnar_system(df_severe,  system = "B")
df_elderly_A <- apply_mnar_system(df_elderly, system = "A")
df_elderly_B <- apply_mnar_system(df_elderly, system = "B")

pipelines <- list(
  severe_A  = df_severe_A,
  severe_B  = df_severe_B,
  elderly_A = df_elderly_A,
  elderly_B = df_elderly_B
)

make_meta <- function(base_meta, dt) {
  m <- base_meta
  sys <- attr(dt, "mnar_system")
  m$mnar_system    <- sys
  m$mnar_indicator <- attr(dt, "mnar_indicator")
  if (sys == "B") {
    m$categorical <- c(m$categorical, attr(dt, "mnar_indicator"))
  }
  m$impute_targets <- m$mar_vars
  m
}

meta_severe_A  <- make_meta(var_meta, df_severe_A)
meta_severe_B  <- make_meta(var_meta, df_severe_B)
meta_elderly_A <- make_meta(var_meta, df_elderly_A)
meta_elderly_B <- make_meta(var_meta, df_elderly_B)

pipeline_meta <- list(
  severe_A  = meta_severe_A,
  severe_B  = meta_severe_B,
  elderly_A = meta_elderly_A,
  elderly_B = meta_elderly_B
)

report_missing <- function(dt, label, meta) {
  cols <- c(meta$mnar_vars, meta$mnar_indicator, meta$mar_vars)
  na   <- dt[, vapply(.SD, function(x) sum(is.na(x)), integer(1)), .SDcols = cols]
  message(sprintf("---- [%s] (系統%s) n=%d, good=%d ----",
                  label, meta$mnar_system, nrow(dt), dt[, sum(good)]))
  print(na)
}
for (nm in names(pipelines)) {
  report_missing(pipelines[[nm]], nm, pipeline_meta[[nm]])
}

saveRDS(pipelines,     file.path(OUT_DIR, "bnb_pipelines_mnar.rds"))
saveRDS(pipeline_meta, file.path(OUT_DIR, "bnb_pipeline_meta.rds"))

message("§2 完了。4本のパイプライン(severe_A/severe_B/elderly_A/elderly_B)を生成・保存しました。")
