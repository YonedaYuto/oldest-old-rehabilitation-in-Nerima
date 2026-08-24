library(data.table)
library(here)

need_pkg <- function(p) if (!requireNamespace(p, quietly = TRUE))
  stop(sprintf("パッケージ '%s' が必要です。install.packages('%s')", p, p))
need_pkg("smcfcs")
need_pkg("Hmisc")

OUT_DIR <- here::here("data")

N_IMPS_BY <- list(
  severe_A  = 100L, severe_B  = 100L,
  elderly_A = 50L,  elderly_B = 50L
)
N_IMPS_DEFAULT <- 50L
NUMIT    <- 100L
RJLIMIT  <- 5000L
IMP_SEED <- 2024L
USE_PARALLEL <- TRUE

IMP_BOUNDS   <- list(MMSE_in = c(0, 30), BBS_in = c(0, 56))
ROUND_TO_INT <- FALSE

EXTRA_SM_TERMS <- list(
  severe_A = character(0), severe_B = character(0),
  elderly_A = character(0), elderly_B = character(0)
)

if (!exists("centered_pipelines"))
  centered_pipelines <- readRDS(file.path(OUT_DIR, "bnb_pipelines_centered.rds"))
if (!exists("rcs_meta"))
  rcs_meta <- readRDS(file.path(OUT_DIR, "bnb_rcs_meta.rds"))

fmt <- function(x) formatC(x, format = "g", digits = 15)
knots_to_c <- function(k) paste0("c(", paste(fmt(k), collapse = ", "), ")")

build_smcfcs_inputs <- function(dt, meta, extra_terms = character(0)) {
  sel_cont <- meta$selected_cont
  sel_cat  <- meta$selected_cat
  rcs_cont <- intersect(meta$rcs_cont, sel_cont)
  centers  <- meta$centers
  n_knots  <- meta$n_knots
  n_nl     <- max(0L, n_knots - 2L)

  cols_raw <- intersect(sel_cont, names(dt))
  cols_c   <- paste0(cols_raw, "_c")
  miss_c   <- setdiff(cols_c, names(dt))
  if (length(miss_c) > 0)
    stop("centered 列が見つかりません（§3 を再実行してください）: ",
         paste(miss_c, collapse = ", "))

  dfi <- copy(dt)[, c("good", cols_raw, cols_c, sel_cat), with = FALSE]
  dfi[, rowid := .I]
  setcolorder(dfi, c("rowid", "good", cols_raw, cols_c, sel_cat))

  basis_cols <- character(0)
  for (v in rcs_cont) {
    vc <- paste0(v, "_c"); kn <- meta$rcs_knots[[v]]
    if (is.null(kn) || n_nl < 1L) next
    b <- Hmisc::rcspline.eval(dfi[[vc]], knots = kn, inclx = FALSE)
    b <- as.matrix(b)
    for (j in seq_len(ncol(b))) {
      nm <- paste0(v, "_c", j)
      dfi[, (nm) := b[, j]]
      basis_cols <- c(basis_cols, nm)
    }
  }

  vars <- names(dfi)

  meth <- setNames(rep("", length(vars)), vars)

  na_any <- vapply(dfi, function(x) anyNA(x), logical(1))
  imp_raw <- intersect(cols_raw, names(na_any)[na_any])
  imp_cat <- intersect(sel_cat,  names(na_any)[na_any])

  for (v in imp_raw) meth[v] <- "norm"
  for (v in imp_cat) {
    f <- dfi[[v]]
    nlev <- if (is.factor(f)) nlevels(f) else length(unique(stats::na.omit(f)))
    meth[v] <- if (nlev <= 2L) "logreg" else if (is.ordered(f)) "podds" else "mlogit"
  }
  unexpected <- setdiff(c(imp_raw, imp_cat),
                        c("MMSE_in", "BBS_in", "JCS_bin"))
  if (length(unexpected) > 0)
    warning("[", meta$mnar_system, "] 想定外の欠損変数を代入対象に検出: ",
            paste(unexpected, collapse = ", "))

  for (v in imp_raw) {
    vc <- paste0(v, "_c")
    if (vc %in% vars) meth[vc] <- sprintf("%s - (%s)", v, fmt(centers[[v]]))
  }
  for (v in intersect(rcs_cont, imp_raw)) {
    vc <- paste0(v, "_c"); kn <- meta$rcs_knots[[v]]
    for (j in seq_len(n_nl)) {
      nm <- paste0(v, "_c", j)
      if (nm %in% vars)
        meth[nm] <- sprintf("Hmisc::rcspline.eval(%s, knots = %s)[, %d]",
                            vc, knots_to_c(kn), j)
    }
  }

  predmat <- matrix(1L, length(vars), length(vars), dimnames = list(vars, vars))
  diag(predmat) <- 0L
  derived <- c("good", "rowid", cols_c, basis_cols)
  predmat[, intersect(derived, vars)] <- 0L

  cont_terms <- unlist(lapply(sel_cont, function(v) {
    vc <- paste0(v, "_c")
    if (v %in% rcs_cont && n_nl >= 1L) c(vc, paste0(v, "_c", seq_len(n_nl))) else vc
  }), use.names = FALSE)
  rhs <- paste(c(cont_terms, sel_cat, extra_terms), collapse = " + ")
  smformula <- paste("good ~", rhs)

  dfi[, good := as.integer(good)]
  list(dfi = as.data.frame(dfi), meth = meth, predmat = predmat,
       smformula = smformula, imputed = c(imp_raw, imp_cat),
       basis_cols = basis_cols)
}

run_one_smcfcs <- function(dat, smf, meth, pm, numit, rjlimit) {
  smcfcs::smcfcs(originaldata = dat, smtype = "logistic", smformula = smf,
                 method = meth, predictorMatrix = pm,
                 m = 1, numit = numit, rjlimit = rjlimit)
}

stack_long <- function(orig_df, imp_obj) {
  m <- length(imp_obj$impDatasets)
  base0 <- data.table(.imp = 0L, .id = seq_len(nrow(orig_df)), as.data.table(orig_df))
  parts <- lapply(seq_len(m), function(k) {
    d <- imp_obj$impDatasets[[k]]
    data.table(.imp = k, .id = seq_len(nrow(d)), as.data.table(d))
  })
  rbindlist(c(list(base0), parts), use.names = TRUE, fill = TRUE)
}

clamp_one <- function(d, bounds, centers, rcs_cont, rcs_knots, n_nl, round_int) {
  for (v in names(bounds)) {
    if (!v %in% names(d)) next
    lo <- bounds[[v]][1]; hi <- bounds[[v]][2]
    x  <- pmin(pmax(as.numeric(d[[v]]), lo), hi)
    if (round_int) x <- round(x)
    d[[v]] <- x
    vc <- paste0(v, "_c")
    if (vc %in% names(d)) d[[vc]] <- x - centers[[v]]
    if (v %in% rcs_cont && n_nl >= 1L && !is.null(rcs_knots[[v]])) {
      b <- as.matrix(Hmisc::rcspline.eval(d[[vc]], knots = rcs_knots[[v]], inclx = FALSE))
      for (j in seq_len(min(ncol(b), n_nl))) d[[paste0(v, "_c", j)]] <- b[, j]
    }
  }
  d
}

have_furrr <- USE_PARALLEL &&
  requireNamespace("furrr", quietly = TRUE) &&
  requireNamespace("future", quietly = TRUE)
have_tic   <- requireNamespace("tictoc", quietly = TRUE)

if (have_furrr) {
  max_m <- max(unlist(N_IMPS_BY), N_IMPS_DEFAULT)
  N_WORKERS <- min(max_m, max(1L, parallel::detectCores() - 1L))
  future::plan(future::multisession, workers = N_WORKERS)
  message(sprintf("並列実行: %d ワーカ", N_WORKERS))
} else {
  message("逐次実行（furrr/future 未使用）")
}

imp_index <- character(0)

for (nm in names(centered_pipelines)) {
  dt   <- centered_pipelines[[nm]]
  meta <- rcs_meta[[nm]]
  extra <- if (!is.null(EXTRA_SM_TERMS[[nm]])) EXTRA_SM_TERMS[[nm]] else character(0)

  inp <- build_smcfcs_inputs(dt, meta, extra_terms = extra)

  n_imps <- if (!is.null(N_IMPS_BY[[nm]])) N_IMPS_BY[[nm]] else N_IMPS_DEFAULT

  message(sprintf("==== [%s] 代入開始  n=%d, m=%d, numit=%d ====",
                  nm, nrow(inp$dfi), n_imps, NUMIT))
  message("  smformula: ", inp$smformula)
  message("  代入対象 : ", paste(inp$imputed, collapse = ", "))

  if (have_tic) tictoc::tic(nm)

  if (have_furrr) {
    imps <- furrr::future_map(
      seq_len(n_imps),
      function(i) run_one_smcfcs(inp$dfi, inp$smformula, inp$meth, inp$predmat,
                                 NUMIT, RJLIMIT),
      .options = furrr::furrr_options(seed = IMP_SEED, packages = c("smcfcs", "Hmisc"))
    )
  } else {
    set.seed(IMP_SEED)
    imps <- lapply(seq_len(n_imps),
                   function(i) run_one_smcfcs(inp$dfi, inp$smformula, inp$meth,
                                              inp$predmat, NUMIT, RJLIMIT))
  }

  combined <- imps[[1]]
  combined$impDatasets <- do.call(c, lapply(imps, function(z) z$impDatasets))

  if (length(IMP_BOUNDS) > 0) {
    n_nl <- max(0L, meta$n_knots - 2L)
    combined$impDatasets <- lapply(
      combined$impDatasets, clamp_one,
      bounds = IMP_BOUNDS, centers = meta$centers,
      rcs_cont = intersect(meta$rcs_cont, meta$selected_cont),
      rcs_knots = meta$rcs_knots, n_nl = n_nl, round_int = ROUND_TO_INT)
  }

  if (have_tic) tictoc::toc()
  wn <- warnings(); if (length(wn) > 0) print(wn)

  long_dt <- stack_long(inp$dfi, combined)

  saveRDS(list(
    smcfcs    = combined,
    long      = long_dt,
    smformula = inp$smformula,
    meth      = inp$meth,
    predmat   = inp$predmat,
    imputed   = inp$imputed,
    basis_cols = inp$basis_cols,
    rcs_cont  = meta$rcs_cont,
    n_knots   = meta$n_knots,
    rcs_knots = meta$rcs_knots,
    centers   = meta$centers,
    selected_cont = meta$selected_cont,
    selected_cat  = meta$selected_cat,
    mnar_system   = meta$mnar_system,
    mnar_indicator = meta$mnar_indicator,
    n = nrow(inp$dfi), m = n_imps
  ), file.path(OUT_DIR, sprintf("bnb_imp_%s.rds", nm)))

  imp_index <- c(imp_index, nm)
  message(sprintf("[%s] 代入完了。m=%d データセットを保存。", nm, n_imps))
}

if (have_furrr) future::plan(future::sequential)

saveRDS(imp_index, file.path(OUT_DIR, "bnb_imp_index.rds"))
message("§5 完了。代入結果を ", normalizePath(OUT_DIR),
        " に bnb_imp_<pipeline>.rds として保存しました（§6 で分離判定→推定分岐に進む）。")
