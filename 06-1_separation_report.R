library(data.table)
library(here)

OUT_DIR <- here::here("data")

if (!exists("est_index"))
  est_index <- readRDS(file.path(OUT_DIR, "bnb_estimation_index.rds"))

summary_tab <- est_index[, .(
  pipeline, mnar_system, n, m, n_sep,
  pct_sep   = round(100 * n_sep / m, 1),
  any_separation,
  inf_terms = fifelse(inf_terms == "", "(none)", inf_terms),
  path
)]

message("==== §6 分離判定サマリ（全パイプライン）====")
print(summary_tab)

for (nm in est_index$pipeline) {
  obj <- readRDS(file.path(OUT_DIR, sprintf("bnb_estimation_%s.rds", nm)))
  sp  <- obj$separation

  message(sprintf("\n---- [%s] (系統%s)  n=%d, m=%d ----",
                  nm, obj$mnar_system, obj$n, sp$m))
  message(sprintf("  分離あり代入: %d / %d  (%.1f%%)",
                  sp$n_sep, sp$m, 100 * sp$n_sep / sp$m))
  message(sprintf("  無限化した係数: %s",
                  if (length(sp$inf_terms) == 0) "(none)" else paste(sp$inf_terms, collapse = ", ")))
  message(sprintf("  判定: %s  → 採用経路: %s",
                  if (isTRUE(sp$any)) "分離あり" else "分離なし", toupper(obj$path)))

  st <- as.data.table(sp$table)
  if (any(st$separated, na.rm = TRUE)) {
    message("  分離が出た代入（imp 番号と無限化項）:")
    print(st[separated == TRUE, .(imp, inf_terms)])
  } else {
    message("  → 全代入データで分離なし。")
  }
  if (any(!st$detect_ok)) {
    message(sprintf("  ※分離検出に失敗した代入が %d 件あります（要確認）:", sum(!st$detect_ok)))
    print(st[detect_ok == FALSE, .(imp)])
  }
}

all_glm <- all(est_index$path == "glm")
message(sprintf("\n結論: %s",
  if (all_glm)
    "全サブセット・全系統で分離なし → 4本とも GLM + Rubin's rule 経路を採用。"
  else
    paste0("分離ありのパイプラインが存在 → Firth 経路: ",
           paste(est_index[path == "firth", pipeline], collapse = ", "))))
