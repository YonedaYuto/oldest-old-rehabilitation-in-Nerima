library(data.table)
library(here)

need_pkg <- function(p) if (!requireNamespace(p, quietly = TRUE))
  stop(sprintf("パッケージ '%s' が必要です。install.packages('%s')", p, p))
need_pkg("mice")
need_pkg("logistf")

HAS_DETECTSEP <- requireNamespace("detectseparation", quietly = TRUE)
HAS_BRGLM2    <- requireNamespace("brglm2",          quietly = TRUE)
if (!HAS_DETECTSEP && !HAS_BRGLM2)
  stop("分離検出に detectseparation か brglm2 のいずれかが必要です。install.packages('detectseparation')")

LABELS_PATH <- here::here("00_labels.R")
if (!exists("relabel_vars")) {
  if (!file.exists(LABELS_PATH)) stop("00_labels.R が見つかりません: ", LABELS_PATH)
  source(LABELS_PATH)
}

OUT_DIR  <- here::here("data")
have_tic <- requireNamespace("tictoc", quietly = TRUE)
PIPELINES <- names(PIPE_LABELS)

USER_MODEL <- list(
  severe_A  = NULL,
  severe_B  = NULL,
  elderly_A = "BBS_in + cFIM_in + STEF_better_in + mFIM_in + support_in + MMSE_in",
  elderly_B = "BBS_in + cFIM_in + class + mFIM_in + MMSE_in + STEF_better_in + support_in"
)

CI_LEVEL <- 0.95
LOGISTF_PL <- FALSE
FORCE_PATH <- list(severe_A = NULL, severe_B = NULL, elderly_A = NULL, elderly_B = NULL)
CLIP_VARS  <- list(severe_A = NULL, severe_B = NULL, elderly_A = NULL, elderly_B = NULL)

BETA_MAX <- 15
SE_MAX   <- 30
PROB_EPS <- 1e-6
AUTO_FALLBACK_FIRTH <- TRUE

EPV_PER_VAR <- 10

auto_path  <- function(nm) file.path(OUT_DIR, sprintf("bnb_final_auto_%s.rds", nm))
final_path <- function(nm) file.path(OUT_DIR, sprintf("bnb_final_%s.rds",      nm))
AUTO_INDEX  <- file.path(OUT_DIR, "bnb_final_index_auto.rds")
FINAL_INDEX <- file.path(OUT_DIR, "bnb_final_index.rds")

detect_one <- function(smf, dat) {
  detector <- if (HAS_DETECTSEP) detectseparation::detect_separation
              else                brglm2::detect_separation
  ds <- tryCatch(stats::glm(smf, data = dat, family = binomial("logit"), method = detector),
                 error = function(e) e)
  if (inherits(ds, "error")) return(list(separated = NA, inf_terms = character(0), ok = FALSE))
  betas <- if (!is.null(ds$betas)) ds$betas
           else if (!is.null(ds$coefficients)) ds$coefficients else numeric(0)
  inf_terms <- names(betas)[is.infinite(betas)]
  sep_flag  <- if (!is.null(ds$separation)) isTRUE(ds$separation)
               else if (!is.null(ds$outcome)) isTRUE(ds$outcome)
               else length(inf_terms) > 0
  list(separated = sep_flag, inf_terms = inf_terms, ok = TRUE)
}

glm_unstable_one <- function(f) {
  cf <- stats::coef(f); se <- suppressWarnings(sqrt(diag(stats::vcov(f)))); pr <- stats::fitted(f)
  cf_f <- cf[is.finite(cf)]; se_f <- se[is.finite(se)]
  bad <- any(!is.finite(cf)) || any(!is.finite(se)) ||
         (length(cf_f) > 0 && any(abs(cf_f) > BETA_MAX)) ||
         (length(se_f) > 0 && any(se_f > SE_MAX)) ||
         any(pr < PROB_EPS | pr > 1 - PROB_EPS) || !isTRUE(f$converged)
  list(bad = bad,
       max_abeta = if (length(cf_f)) max(abs(cf_f)) else NA_real_,
       max_se    = if (length(se_f)) max(se_f) else NA_real_)
}

sep_or_unstable <- function(smf, imps) {
  exact_any <- FALSE; det_err <- FALSE; inf_terms <- character(0)
  for (d in imps) {
    z <- detect_one(smf, d)
    if (!isTRUE(z$ok)) { det_err <- TRUE; next }
    if (isTRUE(z$separated)) { exact_any <- TRUE; inf_terms <- union(inf_terms, z$inf_terms) }
  }
  if (exact_any)
    return(list(sep = TRUE, kind = "exact", n_unstable = NA_integer_, m = length(imps),
                inf_terms = inf_terms, det_error = det_err, max_abeta = NA_real_, max_se = NA_real_))
  n_unstable <- 0L; mab <- 0; mse <- 0
  for (d in imps) {
    f <- suppressWarnings(stats::glm(smf, binomial(), d))
    u <- glm_unstable_one(f)
    if (isTRUE(u$bad)) n_unstable <- n_unstable + 1L
    if (is.finite(u$max_abeta)) mab <- max(mab, u$max_abeta)
    if (is.finite(u$max_se))    mse <- max(mse, u$max_se)
  }
  list(sep = n_unstable > 0, kind = if (n_unstable > 0) "near/unstable" else "none",
       n_unstable = n_unstable, m = length(imps), inf_terms = inf_terms,
       det_error = det_err, max_abeta = mab, max_se = mse)
}

make_formula <- function(term_labels) {
  if (length(term_labels) == 0) return(stats::as.formula("good ~ 1"))
  stats::reformulate(termlabels = term_labels, response = "good")
}

pool_rubin <- function(fitlist, n_complete, ci_level = CI_LEVEL) {
  m <- length(fitlist); betas <- lapply(fitlist, coef); terms <- names(betas[[1]])
  Vdg <- lapply(fitlist, function(f) diag(as.matrix(vcov(f))))
  getv <- function(lst, t) vapply(lst, function(z) unname(z[t]), numeric(1))
  Qbar <- vapply(terms, function(t) mean(getv(betas, t)), numeric(1))
  B    <- vapply(terms, function(t) stats::var(getv(betas, t)), numeric(1))
  Ubar <- vapply(terms, function(t) mean(getv(Vdg, t)), numeric(1))
  Tvar <- Ubar + (1 + 1/m) * B; riv <- (1 + 1/m) * B / Ubar; lambda <- (1 + 1/m) * B / Tvar
  k <- length(terms); dfc <- max(n_complete - k, 1)
  df_old <- ifelse(lambda > 0, (m - 1) / lambda^2, Inf)
  df_obs <- ((dfc + 1)/(dfc + 3)) * dfc * (1 - lambda)
  df <- df_old * df_obs / (df_old + df_obs); fmi <- (riv + 2/(df + 3)) / (riv + 1)
  se <- sqrt(Tvar); z <- stats::qnorm(1 - (1 - ci_level)/2)
  data.table(term = terms, estimate = Qbar, se = se, OR = exp(Qbar),
             lcl = exp(Qbar - z*se), ucl = exp(Qbar + z*se), riv = riv, fmi = fmi)
}

estimate_glm <- function(smf, imps, n, ci_level = CI_LEVEL) {
  sep_warn <- FALSE
  fits <- lapply(imps, function(d) withCallingHandlers(glm(smf, binomial(), d),
    warning = function(w) { if (grepl("fitted probabilities", conditionMessage(w))) sep_warn <<- TRUE
      invokeRestart("muffleWarning") }))
  pooled <- mice::pool(fits)
  ps <- as.data.table(summary(pooled, conf.int = TRUE, conf.level = ci_level))
  ci_cols <- grep("%", names(ps), value = TRUE)
  if (length(ci_cols) >= 2) { lo <- ps[[ci_cols[1]]]; hi <- ps[[ci_cols[2]]] }
  else { zz <- stats::qnorm(1 - (1 - ci_level)/2); lo <- ps$estimate - zz*ps$std.error; hi <- ps$estimate + zz*ps$std.error }
  list(pooled = data.table(term = ps$term, estimate = ps$estimate, se = ps$std.error,
                           OR = exp(ps$estimate), lcl = exp(lo), ucl = exp(hi), p = ps$p.value),
       glm_sep_warn = sep_warn)
}

estimate_firth <- function(smf, imps, n, clip_vars = NULL, ci_level = CI_LEVEL) {
  fitlist <- lapply(imps, function(d) logistf::logistf(smf, data = d, pl = LOGISTF_PL, dataout = TRUE))
  rb <- pool_rubin(fitlist, n_complete = n, ci_level = ci_level)
  if (have_tic) tictoc::tic("  CLIP.confint(final)")
  cc <- tryCatch(logistf::CLIP.confint(obj = fitlist, variable = clip_vars, ci.level = ci_level, pvalue = TRUE),
                 error = function(e) e)
  if (have_tic) tictoc::toc()
  pooled <- copy(rb)[, .(term, estimate, se, OR, lcl, ucl, fmi)]; clip_ok <- !inherits(cc, "error")
  if (clip_ok) {
    ci_mat <- if (!is.null(cc$ci)) as.matrix(cc$ci) else if (!is.null(cc$confint)) as.matrix(cc$confint) else NULL
    vnm <- if (!is.null(cc$variable)) cc$variable else if (!is.null(rownames(ci_mat))) rownames(ci_mat) else NULL
    if (!is.null(ci_mat) && !is.null(vnm)) {
      clip_dt <- data.table(term = vnm, lcl_clip = ci_mat[, 1], ucl_clip = ci_mat[, 2],
                            p = if (!is.null(cc$pvalue)) cc$pvalue else NA_real_)
      pooled <- merge(pooled, clip_dt, by = "term", all.x = TRUE, sort = FALSE)
      hv <- is.finite(pooled$lcl_clip) & is.finite(pooled$ucl_clip)
      pooled[hv, `:=`(lcl = exp(lcl_clip), ucl = exp(ucl_clip))]
      pooled[, c("lcl_clip", "ucl_clip") := NULL]
    } else clip_ok <- FALSE
  }
  if (!clip_ok) { warning("CLIP.confint 失敗/未対応 → CI は罰則付き Rubin 正規近似。"); pooled[, p := NA_real_] }
  list(pooled = pooled, clip_ok = clip_ok)
}

inter_side_labels <- function(g, groups, rcs_mode) {
  labs <- groups[[g]]
  if (identical(rcs_mode, "full")) return(labs)
  lin <- labs[grepl("_c$", labs)]
  if (length(lin) >= 1) lin else labs
}
make_inter_labels <- function(g1, g2, groups, rcs_mode) {
  s1 <- inter_side_labels(g1, groups, rcs_mode)
  s2 <- inter_side_labels(g2, groups, rcs_mode)
  as.vector(outer(s1, s2, function(a, b) paste0(a, ":", b)))
}

inter_disp <- function(key) {
  pr <- strsplit(key, ":", fixed = TRUE)[[1]]
  paste(relabel_vars(pr[1]), "x", relabel_vars(pr[2]))
}

read_safe <- function(f) if (file.exists(f)) readRDS(f) else NULL

parse_user_spec <- function(spec, groups, nm) {
  spec <- trimws(spec)
  if (!nzchar(spec)) stop(sprintf("[%s] USER_MODEL が空文字です。", nm))
  tt <- tryCatch(stats::terms(stats::as.formula(paste("~", spec)), keep.order = FALSE),
                 error = function(e)
                   stop(sprintf("[%s] USER_MODEL の式を解釈できません: '%s'（%s）",
                                nm, spec, conditionMessage(e))))
  labs <- attr(tt, "term.labels")
  if (length(labs) == 0)
    stop(sprintf("[%s] USER_MODEL に項がありません（切片のみは不可）: '%s'", nm, spec))
  is_int    <- grepl(":", labs, fixed = TRUE)
  main_keys <- labs[!is_int]
  int_pairs <- lapply(labs[is_int], function(s) strsplit(s, ":", fixed = TRUE)[[1]])
  bad_order <- vapply(int_pairs, function(p) length(p) != 2L, logical(1))
  if (any(bad_order))
    stop(sprintf("[%s] 3 変数以上の交互作用は非対応です: %s",
                 nm, paste(labs[is_int][bad_order], collapse = ", ")))
  used    <- unique(c(main_keys, unlist(int_pairs)))
  unknown <- setdiff(used, names(groups))
  if (length(unknown) > 0)
    stop(sprintf("[%s] 未知の変数キー: %s\n  使用可能なキー: %s",
                 nm, paste(unknown, collapse = ", "), paste(names(groups), collapse = ", ")))
  ord <- names(groups)
  norm_pair <- function(p) { p2 <- p[order(match(p, ord))]; list(pair = p2, key = paste0(p2[1], ":", p2[2])) }
  np      <- lapply(int_pairs, norm_pair)
  int_pairs_ord <- lapply(np, `[[`, "pair")
  int_keys      <- vapply(np, `[[`, character(1), "key")
  for (p in int_pairs_ord) {
    miss <- setdiff(p, main_keys)
    if (length(miss) > 0)
      message(sprintf("  ※[%s] 交互作用 %s の主効果 %s が主効果に含まれていません（先行研究上の意図でなければ追加を検討）。",
                      nm, paste(relabel_vars(p), collapse = " x "),
                      paste(relabel_vars(miss), collapse = ", ")))
  }
  list(main_keys = main_keys[order(match(main_keys, ord))],
       int_pairs = int_pairs_ord, int_keys = int_keys)
}

print_key_legend <- function(nm, auto) {
  groups <- auto$groups
  keys   <- names(groups)
  message(sprintf("  使用可能な変数キー（USER_MODEL に書くトークン）: %d 個", length(keys)))
  leg <- data.table(token = keys, label = relabel_vars(keys),
                    in_provisional = keys %in% auto$main_vars)
  print(leg, row.names = FALSE)
  message("  §8 暫定モデルの主効果: ",
          if (length(auto$main_vars)) paste(relabel_vars(auto$main_vars), collapse = ", ") else "(なし)")
  message("  §10 最終で選ばれた交互作用: ",
          if (length(auto$selected_int_keys))
            paste(vapply(auto$selected_int_keys, inter_disp, character(1)), collapse = ", ")
          else "(なし)")
}

process_pipeline <- function(nm) {
  if (!file.exists(auto_path(nm))) {
    cur <- read_safe(final_path(nm))
    if (is.null(cur)) { message(sprintf("[%s] bnb_final_%s.rds が無い。§10 を先に実行。スキップ。", nm, nm)); return(NULL) }
    saveRDS(cur, auto_path(nm))
    message(sprintf("[%s] §10 自動最終モデルを bnb_final_auto_%s.rds に退避。", nm, nm))
  }
  auto <- readRDS(auto_path(nm))

  message(sprintf("\n================ [%s] §10-1 手動キュレーション ================", nm))
  print_key_legend(nm, auto)

  spec <- USER_MODEL[[nm]]

  if (is.null(spec)) {
    saveRDS(auto, final_path(nm))
    message(sprintf("  USER_MODEL[['%s']]=NULL → §10 の最終モデルをそのまま使用（bnb_final_%s.rds を auto から復元）。", nm, nm))
    return(data.table(pipeline = nm, mnar_system = auto$mnar_system, n = auto$n, m = auto$m,
                      curated = FALSE, path = auto$path,
                      n_main = length(auto$main_vars),
                      n_int = length(auto$selected_int_keys),
                      formula = auto$smformula_final))
  }

  groups <- auto$groups
  ix_int  <- read_safe(file.path(OUT_DIR, sprintf("bnb_interactions_%s.rds", nm)))
  rcs_mode <- if (!is.null(ix_int) && !is.null(ix_int$settings$rcs_mode)) ix_int$settings$rcs_mode else "linear"

  parsed <- parse_user_spec(spec, groups, nm)

  main_labels <- unlist(groups[parsed$main_keys], use.names = FALSE)
  selected_int_groups <- list(); selected_int_labels <- character(0)
  for (i in seq_along(parsed$int_pairs)) {
    p   <- parsed$int_pairs[[i]]; key <- parsed$int_keys[i]
    lab <- make_inter_labels(p[1], p[2], groups, rcs_mode)
    selected_int_groups[[key]] <- lab
    selected_int_labels <- c(selected_int_labels, lab)
  }
  final_labels <- c(main_labels, selected_int_labels)
  smf_final    <- make_formula(final_labels)
  message("  確定モデル（ユーザー指定）: ", paste(deparse(smf_final), collapse = ""))
  message(sprintf("  主効果 %d 群 / 交互作用 %d 件%s",
                  length(parsed$main_keys), length(parsed$int_keys),
                  if (length(parsed$int_keys))
                    paste0("（", paste(vapply(parsed$int_keys, inter_disp, character(1)), collapse = ", "), "）")
                  else ""))

  obj <- read_safe(file.path(OUT_DIR, sprintf("bnb_imp_%s.rds", nm)))
  if (is.null(obj)) stop(sprintf("[%s] bnb_imp_%s.rds が無い（§5 出力）。", nm, nm))
  n0 <- auto$n0; m <- length(obj$smcfcs$impDatasets)
  bad <- sort(unique(auto$outlier_ids)); bad <- bad[bad %in% seq_len(n0)]
  imps <- if (length(bad) > 0) lapply(obj$smcfcs$impDatasets, function(d) d[-bad, , drop = FALSE])
          else obj$smcfcs$impDatasets
  n <- n0 - length(bad)
  message(sprintf("  縮小データ: n=%d（除外 %d 件）, m=%d", n, length(bad), m))

  y_out  <- imps[[1]][["good"]]
  n_events <- as.integer(min(table(factor(y_out))))
  pred_df  <- length(stats::coef(suppressWarnings(glm(smf_final, binomial(), imps[[1]])))) - 1L
  message(sprintf("  参考 EPV: 予測子df=%d, events=%d（少数側）, 目安上限 events/%g=%.1f%s",
                  pred_df, n_events, EPV_PER_VAR, n_events / EPV_PER_VAR,
                  if (pred_df > n_events / EPV_PER_VAR) "  ※目安超過（過適合・分離に注意）" else ""))

  su <- sep_or_unstable(smf_final, imps)
  any_sep <- isTRUE(su$sep)
  if (su$kind == "exact")
    message(sprintf("  分離 判定: 厳密(準)完全分離あり | 無限化: %s", paste(su$inf_terms, collapse = ", ")))
  else if (su$kind == "near/unstable")
    message(sprintf("  分離 判定: near-separation/不安定 %d/%d 代入（|β|max=%.1f, SEmax=%.1f）",
                    su$n_unstable, su$m, su$max_abeta, su$max_se))
  else
    message(sprintf("  分離 判定: 分離なし・glm 安定（|β|max=%.1f, SEmax=%.1f）", su$max_abeta, su$max_se))

  forced <- FORCE_PATH[[nm]]
  path <- if (!is.null(forced)) { if (!forced %in% c("glm","firth")) stop("FORCE_PATH は 'glm'/'firth'。"); forced }
          else if (any_sep) "firth" else "glm"
  message(sprintf("  → 推定経路: %s", toupper(path)))

  sel_desc <- "user-curated model (variables fixed by prior knowledge; §10-1)"
  if (have_tic) tictoc::tic(paste0(nm, " §10-1 estimate"))
  if (path == "glm") {
    fin <- estimate_glm(smf_final, imps, n, ci_level = CI_LEVEL); pooled <- fin$pooled
    pooled_unstable <- any(!is.finite(pooled$estimate)) || any(abs(pooled$estimate) > BETA_MAX) ||
                       any(!is.finite(pooled$se)) || any(pooled$se > SE_MAX) || isTRUE(fin$glm_sep_warn)
    if (pooled_unstable && AUTO_FALLBACK_FIRTH && is.null(FORCE_PATH[[nm]])) {
      message(sprintf("  ※glm 推定が不安定（|β|max=%.1f, SEmax=%.1f%s）→ Firth に自動フォールバック。",
                      max(abs(pooled$estimate[is.finite(pooled$estimate)]), 0),
                      max(pooled$se[is.finite(pooled$se)], 0),
                      if (isTRUE(fin$glm_sep_warn)) ", fitted 0/1 警告" else ""))
      path <- "firth"; any_sep <- TRUE
      if (identical(su$kind, "none")) su$kind <- "near/unstable(pooled)"
      su$sep <- TRUE
      fin <- estimate_firth(smf_final, imps, n, clip_vars = CLIP_VARS[[nm]], ci_level = CI_LEVEL)
      pooled <- fin$pooled
      method_detail <- list(estimator = "logistf (Firth)", pool = "Rubin (point) + CLIP.confint (CI/p)",
                            selection = sel_desc, clip_ok = fin$clip_ok, auto_fallback = TRUE)
    } else {
      if (isTRUE(fin$glm_sep_warn))
        message("  ※glm で fitted prob 0/1 警告（準完全分離の兆候）。FORCE_PATH='firth' で感度確認を推奨。")
      method_detail <- list(estimator = "glm", pool = "mice::pool (Rubin)", selection = sel_desc)
    }
  } else {
    fin <- estimate_firth(smf_final, imps, n, clip_vars = CLIP_VARS[[nm]], ci_level = CI_LEVEL); pooled <- fin$pooled
    method_detail <- list(estimator = "logistf (Firth)", pool = "Rubin (point) + CLIP.confint (CI/p)",
                          selection = sel_desc, clip_ok = fin$clip_ok)
  }
  if (have_tic) tictoc::toc()

  show <- copy(pooled)[, `:=`(OR = round(OR, 3), lcl = round(lcl, 3), ucl = round(ucl, 3))]
  if ("p" %in% names(show)) show[, p := signif(p, 4)]
  message("  確定モデル プール推定:")
  print(show[, intersect(c("term", "OR", "lcl", "ucl", "p", "fmi"), names(show)), with = FALSE])

  saveRDS(list(
    pipeline        = nm,
    path            = path,
    smformula_final = paste(deparse(smf_final), collapse = ""),
    smformula_prov  = auto$smformula_prov,
    main_labels     = main_labels,
    main_vars       = parsed$main_keys,
    selected_int_keys   = parsed$int_keys,
    selected_int_labels = selected_int_labels,
    selected_int_groups = selected_int_groups,
    candidate_int_keys  = auto$candidate_int_keys,
    final_labels    = final_labels,
    groups          = groups,
    selection       = list(method = "user_curated", main_effects_fixed = TRUE,
                           note = "Variables hand-selected from prior research (§10-1); estimation/separation handling identical to §10."),
    pooled          = pooled,
    method_detail   = method_detail,
    separation_final = list(any = any_sep, kind = su$kind, n_unstable = su$n_unstable,
                            max_abeta = su$max_abeta, max_se = su$max_se,
                            m = m, inf_terms = su$inf_terms, det_error = su$det_error,
                            beta_max = BETA_MAX, se_max = SE_MAX, prob_eps = PROB_EPS,
                            auto_fallback = isTRUE(method_detail$auto_fallback)),
    outlier_ids = auto$outlier_ids, n0 = n0, n = n, m = m,
    mnar_system = auto$mnar_system,
    selected_cont = auto$selected_cont, selected_cat = auto$selected_cat,
    rcs_cont = auto$rcs_cont, rcs_knots = auto$rcs_knots, n_knots = auto$n_knots,
    centers = auto$centers, mnar_indicator = auto$mnar_indicator, basis_cols = auto$basis_cols,
    curated = TRUE, curated_spec = spec, rcs_mode = rcs_mode
  ), final_path(nm))

  message(sprintf("[%s] §10-1 完了。ユーザー確定モデル（主効果%d＋交互作用%d）を bnb_final_%s.rds に上書き保存。",
                  nm, length(parsed$main_keys), length(parsed$int_keys), nm))

  data.table(pipeline = nm, mnar_system = auto$mnar_system, n = n, m = m,
             curated = TRUE, path = path,
             n_main = length(parsed$main_keys), n_int = length(parsed$int_keys),
             formula = paste(deparse(smf_final), collapse = ""))
}

if (file.exists(FINAL_INDEX) && !file.exists(AUTO_INDEX)) {
  saveRDS(readRDS(FINAL_INDEX), AUTO_INDEX)
  message("§10 自動最終インデックスを bnb_final_index_auto.rds に退避。")
}

curated_index <- rbindlist(lapply(PIPELINES, process_pipeline), fill = TRUE)
saveRDS(curated_index, FINAL_INDEX)

message("\n==== §10-1 サマリ（手動キュレーション後の最終モデル）====")
if (nrow(curated_index))
  print(curated_index[, .(pipeline, mnar_system, n, m, curated, path, n_main, n_int, formula)],
        row.names = FALSE)

message("\n§10-1 完了。")
message("  ・USER_MODEL に書いたパイプラインは、その変数で最終モデルを当てはめ直し bnb_final_<pipeline>.rds を上書き。")
message("  ・USER_MODEL=NULL のパイプラインは §10 の最終モデルをそのまま使用（auto バックアップから復元）。")
message("  ・§10 自動版は bnb_final_auto_<pipeline>.rds / bnb_final_index_auto.rds に保全（再実行は常にこの auto を素材に再構築）。")
message("  ・出力スキーマは §10 と同一なので、§11 以降（11 / 11-1 / 11-2 / 12 / 12-1 …）は無改変で実行できます。")
