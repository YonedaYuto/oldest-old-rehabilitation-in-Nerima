library(data.table)
library(here)

OUT_DIR   <- here::here("data")
PIPELINES <- c("severe_A", "severe_B", "elderly_A", "elderly_B")
LABELS    <- c(severe_A  = "Severe / worst-value",
               severe_B  = "Severe / missing-category",
               elderly_A = "Elderly / worst-value",
               elderly_B = "Elderly / missing-category")

rd <- function(f) if (file.exists(f)) readRDS(f) else NULL
p_ <- function(...) file.path(OUT_DIR, sprintf(...))
`%||%` <- function(a, b) if (is.null(a)) b else a

describe <- function(obj, kind) {
  if (is.null(obj)) return(NULL)
  pooled <- as.data.table(obj$pooled)
  terms  <- setdiff(pooled$term, "(Intercept)")
  ints   <- grep(":", terms, value = TRUE)
  data.table(
    kind        = kind,
    curated     = isTRUE(obj$curated),
    path        = obj$path,
    n0          = obj$n0,
    n           = obj$n,
    n_excluded  = length(obj$outlier_ids %||% integer(0)),
    m           = obj$m,
    n_param     = length(terms),
    n_int       = length(ints),
    interactions= if (length(ints)) paste(ints, collapse = ", ") else "-",
    formula     = obj$smformula_final %||% obj$smformula_prov)
}

summ <- rbindlist(lapply(PIPELINES, function(nm) {
  rows <- rbindlist(list(
    describe(rd(p_("bnb_provisional_%s.rds",      nm)), "main-effect (current)"),
    describe(rd(p_("bnb_provisional_auto_%s.rds", nm)), "main-effect (auto §8)"),
    describe(rd(p_("bnb_final_%s.rds",            nm)), "final (current)"),
    describe(rd(p_("bnb_final_auto_%s.rds",       nm)), "final (auto §10)")), fill = TRUE)
  if (!nrow(rows)) return(NULL)
  cbind(pipeline = nm, label = unname(LABELS[nm]), rows)
}), fill = TRUE)

message("==== 保持されているモデルの一覧 ====")
print(summ[, .(pipeline, kind, curated, path, n0, n, n_excluded, m, n_param, n_int)])

message("\n==== モデル式 ====")
for (i in seq_len(nrow(summ)))
  cat(sprintf("[%s] %-22s %s\n", summ$pipeline[i], summ$kind[i], summ$formula[i]))

message("\n==== 保持されている交互作用 ====")
print(summ[n_int > 0, .(pipeline, kind, interactions)])

message("\n==== 最終モデルの係数（現在保持されているもの）====")
for (nm in PIPELINES) {
  fin <- rd(p_("bnb_final_%s.rds", nm)); if (is.null(fin)) next
  pl  <- as.data.table(fin$pooled)
  cols <- intersect(c("term", "estimate", "se", "OR", "lcl", "ucl", "p"), names(pl))
  cat(sprintf("\n--- %s (%s) | n = %d, m = %d, %s ---\n",
              nm, LABELS[nm], fin$n, fin$m, fin$path))
  print(pl[, ..cols])
}

message("\n==== チェック ====")
for (nm in PIPELINES) {
  fin <- rd(p_("bnb_final_%s.rds", nm)); if (is.null(fin)) next
  terms <- setdiff(as.data.table(fin$pooled)$term, "(Intercept)")
  cat(sprintf("[%s] 推定パラメータ数(切片を除く) = %d / 交互作用 = %d %s\n",
              nm, length(terms), sum(grepl(":", terms)),
              if (any(grepl("BBS_in.*:.*MMSE_in|MMSE_in.*:.*BBS_in", terms)))
                "  <-- BBS×MMSE が保持されている" else ""))
}
idx_f <- rd(file.path(OUT_DIR, "bnb_final_index.rds"))
idx_p <- rd(file.path(OUT_DIR, "bnb_provisional_index.rds"))
if (!is.null(idx_p)) { message("\n-- bnb_provisional_index.rds --"); print(idx_p) }
if (!is.null(idx_f)) { message("\n-- bnb_final_index.rds --");       print(idx_f) }

message("\n§18 完了（読み取りのみ。いかなるファイルも上書きしていない）。")
