library(data.table)
library(forcats)
library(here)

OUT_DIR <- here::here("data")

THR_GOOD   <- 65
THR_SEVERE <- 26
THR_ELDER  <- 90

if (!exists("BNB_small")) {
  stop("BNB_small が見つかりません。先に読み込んでください。")
}
dfo <- as.data.table(BNB_small)

expected_cols <- c(
  "age", "BBS_in", "cFIM_in", "class", "GP_better_in",
  "his_car", "his_CNS", "JCS_in", "mFIM_in", "MMSE_in",
  "sex", "STEF_better_in", "support_in", "walk_speed_in", "X6MD_in",
  "mFIM_out"
)
missing_cols <- setdiff(expected_cols, colnames(dfo))
if (length(missing_cols) > 0) {
  stop("入力データに想定列がありません: ", paste(missing_cols, collapse = ", "))
}

df <- copy(dfo)

df[, JCS_num := suppressWarnings(as.numeric(as.character(JCS_in)))]
df[JCS_num == 300, JCS_num := 30]
unexpected_jcs <- df[!is.na(JCS_in) & is.na(JCS_num), .N]
if (unexpected_jcs > 0) {
  warning(sprintf("JCS_in を数値化できなかった非欠損セルが %d 件あります。", unexpected_jcs))
}

df[, JCS_bin := fcase(
  is.na(JCS_num),                 NA_character_,
  JCS_num %in% c(0, 1, 2, 3),     "clear",
  default =                       "impaired"
)]
df[, JCS_bin := factor(JCS_bin, levels = c("clear", "impaired"))]

df[, class      := fct_relevel(factor(class),      "脳血管")]
df[, his_car    := fct_relevel(factor(his_car),    "無")]
df[, his_CNS    := fct_relevel(factor(his_CNS),    "無")]
df[, sex        := fct_relevel(factor(sex),        "男")]
df[, support_in := fct_relevel(factor(support_in), "なし")]

check_ref <- function(x, ref, name) {
  if (!ref %in% levels(x)) {
    stop(sprintf("変数 %s に基準カテゴリ '%s' が存在しません。実際の水準: %s",
                 name, ref, paste(levels(x), collapse = ", ")))
  }
}
check_ref(df$class,      "脳血管", "class")
check_ref(df$his_car,    "無",     "his_car")
check_ref(df$his_CNS,    "無",     "his_CNS")
check_ref(df$sex,        "男",     "sex")
check_ref(df$support_in, "なし",   "support_in")

if (df[is.na(mFIM_out), .N] > 0) {
  warning("mFIM_out に欠損があります（readme では欠損なしの想定）。")
}
df[, good := as.integer(mFIM_out >= THR_GOOD)]

df_severe  <- df[mFIM_in <= THR_SEVERE]
df_elderly <- df[age      >= THR_ELDER]

var_meta <- list(
  outcome = "good",
  raw_outcome = "mFIM_out",

  continuous = c("age", "BBS_in", "cFIM_in", "GP_better_in", "mFIM_in",
                 "MMSE_in", "STEF_better_in", "walk_speed_in", "X6MD_in"),
  categorical = c("class", "his_car", "his_CNS", "JCS_bin", "sex", "support_in"),

  mar_vars  = c("JCS_in", "MMSE_in", "BBS_in"),
  mnar_vars = c("GP_better_in", "STEF_better_in", "walk_speed_in", "X6MD_in"),

  mar_aux = list(JCS_in = "cFIM_in", MMSE_in = "cFIM_in", BBS_in = "mFIM_in"),

  subset_def = list(severe = "mFIM_in", elderly = "age")
)

report_size <- function(d, label, exp_n, exp_good) {
  n    <- nrow(d)
  ngd  <- d[, sum(good)]
  flag <- if (n == exp_n && ngd == exp_good) "OK" else "*** 不一致 ***"
  message(sprintf("[%s] n=%d (期待%d), good=%d (期待%d)  %s",
                  label, n, exp_n, ngd, exp_good, flag))
}
message("---- サイズ検証（readme との突き合わせ）----")
report_size(df,        "全体", 2400, 1771)
report_size(df_severe, "重症",  576,  144)
report_size(df_elderly,"高齢",  337,  209)

message("---- 欠損数（MAR / MNAR 変数）----")
na_counts <- df[, lapply(.SD, function(x) sum(is.na(x))),
                .SDcols = c(var_meta$mar_vars, var_meta$mnar_vars)]
print(na_counts)

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
saveRDS(df,         file.path(OUT_DIR, "bnb_preprocessed_full.rds"))
saveRDS(df_severe,  file.path(OUT_DIR, "bnb_preprocessed_severe.rds"))
saveRDS(df_elderly, file.path(OUT_DIR, "bnb_preprocessed_elderly.rds"))
saveRDS(var_meta,   file.path(OUT_DIR, "bnb_var_meta.rds"))

message("前処理完了。保存先: ", normalizePath(OUT_DIR))
