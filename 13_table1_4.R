library(data.table)
library(here)

need_pkg <- function(p) if (!requireNamespace(p, quietly = TRUE))
  stop(sprintf("Package '%s' is required. install.packages('%s')", p, p))

LABELS_PATH <- here::here("00_labels.R")
if (!exists("relabel_vars")) {
  if (!file.exists(LABELS_PATH)) stop("00_labels.R not found: ", LABELS_PATH)
  source(LABELS_PATH)
}

OUT_DIR <- here::here("data")
FIG_DIR <- here::here("figures")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

THR_GOOD   <- 65
THR_SEVERE <- 26
THR_ELDER  <- 90

TERM_COL   <- "term"
PCT_DIGITS <- 1
SHOW_MEAN  <- FALSE
SHOW_MINMAX<- FALSE

DERIVED_VARS <- list(
  walk_speed_in = function(b) 10 / b[["X1m_time_in"]]
)
DERIVED_SRC <- list(
  walk_speed_in = "X1m_time_in"
)

TERM_LABEL <- "Length of stay"

if (!exists("BNB_small")) {
  stop("BNB_small not found. Load it first (used as the variable-set reference).")
}
if (!exists("BNB")) {
  stop("BNB not found (the source data containing term). Load it first.")
}

bnb        <- as.data.table(BNB)
small_cols <- names(BNB_small)

direct_cols  <- intersect(small_cols, names(bnb))
derived_cols <- setdiff(small_cols, names(bnb))

no_rule <- setdiff(derived_cols, names(DERIVED_VARS))
if (length(no_rule) > 0) {
  stop("Variables absent from BNB and without a DERIVED_VARS rule: ",
       paste(no_rule, collapse = ", "),
       "\n(Add the transformation used when BNB_small was created to DERIVED_VARS, ",
       "e.g. walk_speed_in = function(b) 10 / b[['X1m_time']])")
}
for (dv in intersect(derived_cols, names(DERIVED_SRC))) {
  miss_src <- setdiff(DERIVED_SRC[[dv]], names(bnb))
  if (length(miss_src) > 0) {
    stop(sprintf("BNB lacks source column(s) needed to reconstruct '%s': %s",
                 dv, paste(miss_src, collapse = ", ")))
  }
}
if (!TERM_COL %in% names(bnb)) {
  stop(sprintf("BNB has no length-of-stay column '%s'. Set TERM_COL. Columns: %s",
               TERM_COL, paste(names(bnb), collapse = ", ")))
}
if (length(derived_cols) > 0) {
  message("Reconstructed from BNB: ",
          paste(sprintf("%s (<- %s)", derived_cols,
                        ifelse(derived_cols %in% names(DERIVED_SRC),
                               vapply(derived_cols, function(x)
                                 paste(DERIVED_SRC[[x]], collapse = ","), character(1)),
                               "rule")), collapse = ", "))
}

BNB_smallplus <- bnb[, direct_cols, with = FALSE]
for (dv in derived_cols) {
  vals <- DERIVED_VARS[[dv]](bnb)
  if (length(vals) != nrow(bnb))
    stop(sprintf("Derived '%s' length (%d) does not match BNB rows (%d).",
                 dv, length(vals), nrow(bnb)))
  nbad <- sum(!is.finite(suppressWarnings(as.numeric(vals))) & !is.na(vals))
  if (nbad > 0)
    message(sprintf("Note: derived '%s' has %d non-finite value(s) (e.g. from 0 in source).",
                    dv, nbad))
  BNB_smallplus[, (dv) := vals]
}
BNB_smallplus[, (TERM_COL) := bnb[[TERM_COL]]]
setcolorder(BNB_smallplus, unique(c(small_cols, TERM_COL)))

n_total   <- nrow(BNB_smallplus)
n_term_na <- sum(is.na(BNB_smallplus[[TERM_COL]]))
message(sprintf("BNB_smallplus built: n=%d, vars=%d (direct %d + derived %d + term), term missing=%d (%.1f%%)",
                n_total, ncol(BNB_smallplus), length(direct_cols),
                length(derived_cols), n_term_na, 100 * n_term_na / max(n_total, 1)))
if (nrow(BNB_smallplus) != nrow(BNB_small)) {
  message(sprintf(paste0("Note: BNB_smallplus rows (%d) differ from BNB_small rows (%d). ",
                         "If BNB_small also filtered rows (not only columns), check intent."),
                  nrow(BNB_smallplus), nrow(BNB_small)))
}

saveRDS(BNB_smallplus, file.path(OUT_DIR, "BNB_smallplus.rds"))
data.table::fwrite(BNB_smallplus, file.path(OUT_DIR, "BNB_smallplus.csv"))

derive_table_vars <- function(d) {
  d <- copy(as.data.table(d))
  if (!"JCS_bin" %in% names(d) && "JCS_in" %in% names(d)) {
    jn <- suppressWarnings(as.numeric(as.character(d$JCS_in)))
    jn[jn == 300] <- 30
    d[, JCS_bin := fcase(
      is.na(jn),               NA_character_,
      jn %in% c(0, 1, 2, 3),   "clear",
      default =                "impaired")]
    d[, JCS_bin := factor(JCS_bin, levels = c("clear", "impaired"))]
  }
  if (!"good" %in% names(d) && "mFIM_out" %in% names(d)) {
    d[, good := as.integer(mFIM_out >= THR_GOOD)]
  }
  d
}

work <- derive_table_vars(BNB_smallplus)

work[, severe  := mFIM_in <= THR_SEVERE]
work[, elderly := age      >= THR_ELDER]

message(sprintf("Subsets: severe=%d, elderly=%d, severe&elderly=%d",
                sum(work$severe, na.rm = TRUE),
                sum(work$elderly, na.rm = TRUE),
                sum(work$severe & work$elderly, na.rm = TRUE)))

saveRDS(work, file.path(OUT_DIR, "bnb_smallplus_table_input.rds"))

as_flag <- function(x) !is.na(x) & x
strata <- list(
  Overall          = rep(TRUE, nrow(work)),
  Severe           = as_flag(work$severe),
  Elderly          = as_flag(work$elderly),
  `Severe&Elderly` = as_flag(work$severe & work$elderly)
)
col_keys  <- names(strata)
col_disp  <- c(Overall = "Overall", Severe = "Severe", Elderly = "Elderly",
               `Severe&Elderly` = "Severe & Elderly")

var_spec <- list(
  list(var = "age",            type = "cont", group = "Demographics"),
  list(var = "sex",            type = "bin",  group = "Demographics",        level = "\u5973"),
  list(var = "class",          type = "cat",  group = "Diagnosis"),
  list(var = "his_car",        type = "bin",  group = "Comorbidity",         level = "\u6709"),
  list(var = "his_CNS",        type = "bin",  group = "Comorbidity",         level = "\u6709"),
  list(var = "support_in",     type = "bin",  group = "Status",              level = "\u3042\u308a"),
  list(var = "mFIM_in",        type = "cont", group = "Function (admission)"),
  list(var = "cFIM_in",        type = "cont", group = "Function (admission)"),
  list(var = "MMSE_in",        type = "cont", group = "Function (admission)"),
  list(var = "JCS_bin",        type = "bin",  group = "Function (admission)", level = "impaired"),
  list(var = "BBS_in",         type = "cont", group = "Function (admission)"),
  list(var = "GP_better_in",   type = "cont", group = "Function (admission)"),
  list(var = "STEF_better_in", type = "cont", group = "Function (admission)"),
  list(var = "walk_speed_in",  type = "cont", group = "Function (admission)"),
  list(var = "X6MD_in",        type = "cont", group = "Function (admission)"),
  list(var = TERM_COL,         type = "cont", group = "Course"),
  list(var = "mFIM_out",       type = "cont", group = "Outcome"),
  list(var = "good",           type = "bin",  group = "Outcome",             level = 1L)
)

fmt_pct <- function(x) formatC(x, format = "f", digits = PCT_DIGITS)

num_fmt <- function(x) {
  if (all(is.na(x))) return("NA")
  if (isTRUE(all.equal(x, round(x))) && all(abs(x) < 1e6))
    formatC(round(x), format = "d", big.mark = "")
  else
    formatC(x, format = "f", digits = 1)
}

label_for <- function(v) {
  if (identical(v, TERM_COL)) return(TERM_LABEL)
  relabel_vars(v)
}

miss_tag <- function(k, n) {
  if (n == 0) return("\u2014")
  if (k == 0) return("0")
  sprintf("%d (%s%%)", k, fmt_pct(100 * k / n))
}

cell_cont <- function(x) {
  v  <- x[!is.na(x)]
  nN <- length(x)
  if (nN == 0) return(list(stat = "\u2014", miss = "\u2014"))
  if (length(v) == 0) return(list(stat = "NA", miss = miss_tag(nN, nN)))
  q <- stats::quantile(v, c(.25, .5, .75), names = FALSE, type = 7)
  stat <- sprintf("%s [%s, %s]", num_fmt(q[2]), num_fmt(q[1]), num_fmt(q[3]))
  if (SHOW_MEAN)
    stat <- sprintf("%s; %s (%s)", stat, num_fmt(mean(v)), num_fmt(stats::sd(v)))
  if (SHOW_MINMAX)
    stat <- sprintf("%s; [%s, %s]", stat, num_fmt(min(v)), num_fmt(max(v)))
  list(stat = stat, miss = miss_tag(sum(is.na(x)), nN))
}

cell_cat_level <- function(x, lvl) {
  nN    <- length(x)
  valid <- sum(!is.na(x))
  if (nN == 0) return("\u2014")
  k <- sum(!is.na(x) & as.character(x) == as.character(lvl))
  if (valid == 0) return("NA")
  sprintf("%d (%s%%)", k, fmt_pct(100 * k / valid))
}

rows <- list()
push <- function(group, label, desc, indent, cells, kind) {
  rows[[length(rows) + 1]] <<- data.table(
    group = group, label = label, desc = desc, indent = indent, kind = kind,
    Overall          = cells[["Overall"]],
    Severe           = cells[["Severe"]],
    Elderly          = cells[["Elderly"]],
    `Severe&Elderly` = cells[["Severe&Elderly"]]
  )
}

cont_desc <- function() {
  d <- "median [Q1, Q3]"
  if (SHOW_MEAN)   d <- paste0(d, "; mean (SD)")
  if (SHOW_MINMAX) d <- paste0(d, "; [min, max]")
  d
}
empty_cells <- function() setNames(as.list(rep("", length(col_keys))), col_keys)

for (sp in var_spec) {
  v <- sp$var
  if (!v %in% names(work)) {
    warning(sprintf("Variable '%s' not in `work`; excluded from Table 1.", v))
    next
  }
  lab <- label_for(v)

  if (sp$type == "cont") {
    stat_cells <- setNames(vector("list", length(col_keys)), col_keys)
    miss_cells <- setNames(vector("list", length(col_keys)), col_keys)
    for (ck in col_keys) {
      cc <- cell_cont(work[[v]][strata[[ck]]])
      stat_cells[[ck]] <- cc$stat
      miss_cells[[ck]] <- cc$miss
    }
    push(sp$group, lab, cont_desc(), 0, stat_cells, "stat")
    push(sp$group, "Missing", "n (%)", 1, miss_cells, "miss")

  } else if (sp$type == "bin") {
    lv     <- sp$level
    lv_lab <- relabel_levels(v, as.character(lv))
    stat_cells <- setNames(vector("list", length(col_keys)), col_keys)
    miss_cells <- setNames(vector("list", length(col_keys)), col_keys)
    for (ck in col_keys) {
      xx <- work[[v]][strata[[ck]]]
      stat_cells[[ck]] <- cell_cat_level(xx, lv)
      miss_cells[[ck]] <- miss_tag(sum(is.na(xx)), length(xx))
    }
    push(sp$group, sprintf("%s = %s", lab, lv_lab), "n (%)", 0, stat_cells, "stat")
    push(sp$group, "Missing", "n (%)", 1, miss_cells, "miss")

  } else if (sp$type == "cat") {
    miss_cells <- setNames(vector("list", length(col_keys)), col_keys)
    for (ck in col_keys) {
      xx <- work[[v]][strata[[ck]]]
      miss_cells[[ck]] <- miss_tag(sum(is.na(xx)), length(xx))
    }
    push(sp$group, lab, "n (%)", 0, empty_cells(), "head")
    push(sp$group, "Missing", "n (%)", 1, miss_cells, "miss")
    lvls <- if (is.factor(work[[v]])) levels(work[[v]]) else
              sort(unique(as.character(work[[v]][!is.na(work[[v]])])))
    for (lv in lvls) {
      lv_lab <- relabel_levels(v, lv)
      cells  <- setNames(vector("list", length(col_keys)), col_keys)
      for (ck in col_keys) cells[[ck]] <- cell_cat_level(work[[v]][strata[[ck]]], lv)
      push(sp$group, lv_lab, "", 1, cells, "level")
    }
  }
}

tab <- rbindlist(rows, use.names = TRUE)

header_label <- function(ck) sprintf("%s (N=%d)", col_disp[[ck]], sum(strata[[ck]]))
col_header <- vapply(col_keys, header_label, character(1))

disp <- copy(tab)
pad_csv <- ifelse(disp$indent == 1, "    ", "")
disp[, Variable := fifelse(
        desc == "", paste0(pad_csv, label),
        fifelse(kind == "miss", paste0(pad_csv, label, ", ", desc),
                paste0(pad_csv, label, "\n", desc)))]
disp_out <- disp[, .(Group = group, Variable,
                     Overall, Severe, Elderly, `Severe&Elderly`)]
setnames(disp_out, c("Overall", "Severe", "Elderly", "Severe&Elderly"),
         unname(col_header))
data.table::fwrite(disp_out, file.path(OUT_DIR, "table1.csv"))

long <- melt(tab, id.vars = c("group", "label", "desc", "indent", "kind"),
             measure.vars = col_keys, variable.name = "stratum",
             value.name = "value")
data.table::fwrite(long, file.path(OUT_DIR, "table1_long.csv"))

esc <- function(s) {
  s <- gsub("&", "&amp;", s, fixed = TRUE)
  s <- gsub("<", "&lt;",  s, fixed = TRUE)
  s <- gsub(">", "&gt;",  s, fixed = TRUE)
  s
}
lab_html <- function(label, desc, indent, kind) {
  pad <- if (indent == 1) "&nbsp;&nbsp;&nbsp;&nbsp;" else ""
  if (identical(kind, "miss"))
    return(paste0("<span class='ml'>", pad, esc(label), ", ", esc(desc), "</span>"))
  out <- paste0("<span class='vname'>", pad, esc(label), "</span>")
  if (nzchar(desc))
    out <- paste0(out, "<br><span class='vdesc'>", esc(desc), "</span>")
  out
}
build_html <- function(tab, col_header, col_keys) {
  th <- paste0("<th>", esc(col_header), "</th>", collapse = "")
  body <- ""
  prev_group <- ""
  for (i in seq_len(nrow(tab))) {
    r <- tab[i]
    if (r$group != prev_group) {
      body <- paste0(body, sprintf("<tr class='grp'><td colspan='%d'>%s</td></tr>\n",
                                   length(col_keys) + 1, esc(r$group)))
      prev_group <- r$group
    }
    cls <- switch(r$kind, miss = "miss", head = "head", level = "level", "stat")
    labcell <- lab_html(r$label, r$desc, r$indent, r$kind)
    tds <- paste0(vapply(col_keys, function(ck)
      sprintf("<td>%s</td>", esc(as.character(r[[ck]]))), character(1)), collapse = "")
    body <- paste0(body, sprintf("<tr class='%s'><td class='lab'>%s</td>%s</tr>\n",
                                 cls, labcell, tds))
  }
  css <- paste0(
    "body{font-family:-apple-system,Segoe UI,Helvetica,Arial,sans-serif;",
    "color:#1a1a1a;margin:28px;}",
    "h2{font-size:22px;margin:0 0 4px;} .sub{color:#555;font-size:14px;margin:0 0 16px;}",
    "table{border-collapse:collapse;font-size:16px;min-width:860px;}",
    "th,td{padding:7px 16px;text-align:left;vertical-align:top;}",
    "thead th{border-bottom:2px solid #333;white-space:nowrap;font-size:16px;}",
    "th:not(:first-child),td:not(:first-child){text-align:right;}",
    "tr.grp td{font-weight:700;background:#f0f1f3;border-top:1px solid #ddd;",
    "padding-top:9px;padding-bottom:5px;font-size:16px;}",
    "td.lab{text-align:left;}",
    ".vname{font-weight:500;} .vdesc{color:#666;font-size:13px;}",
    "tr.miss td{color:#777;font-size:13px;}",
    "tr.head .vname{font-weight:700;}",
    "tfoot td{font-size:12.5px;color:#555;border-top:1px solid #ddd;",
    "padding-top:10px;line-height:1.5;}")
  foot <- paste0(
    "Continuous variables: median [Q1, Q3]. Categorical variables: n (%); ",
    "percentages are computed among non-missing observations. ",
    "&ldquo;Missing&rdquo; rows give the count (and % of the column N) with a missing value; ",
    "a dash (&mdash;) marks a stratum with no eligible observations. ",
    "Severe: admission mFIM &le; 26. Elderly: age &ge; 90. Good outcome: discharge mFIM &ge; 65. ",
    "missingness reporting (3.1).")
  paste0("<!doctype html><html><head><meta charset='utf-8'><style>", css,
         "</style></head><body>",
         "<h2>Table 1. Baseline characteristics of the study sample</h2>",
         "<p class='sub'>Overall and by severity / age strata</p>",
         "<table><thead><tr><th>Characteristic</th>", th, "</tr></thead><tbody>",
         body, "</tbody><tfoot><tr><td colspan='", length(col_keys) + 1, "'>",
         foot, "</td></tr></tfoot></table></body></html>")
}
writeLines(build_html(tab, col_header, col_keys), file.path(FIG_DIR, "table1.html"))

if (requireNamespace("gridExtra", quietly = TRUE) &&
    requireNamespace("grid", quietly = TRUE)) {
  png_df <- as.data.frame(disp_out)
  png_df[[2]] <- gsub("\n", ", ", png_df[[2]], fixed = TRUE)
  th <- gridExtra::ttheme_default(
    core    = list(fg_params = list(hjust = 0, x = 0.02, fontsize = 11)),
    colhead = list(fg_params = list(fontsize = 12, fontface = "bold")))
  g <- gridExtra::tableGrob(png_df, rows = NULL, theme = th)
  grDevices::png(file.path(FIG_DIR, "table1.png"),
                 width = 13, height = 0.34 * nrow(png_df) + 1.4,
                 units = "in", res = 200)
  grid::grid.newpage(); grid::grid.draw(g)
  grDevices::dev.off()
  message("Wrote PNG: figures/table1.png")
} else {
  message("gridExtra/grid not available; PNG skipped (HTML written to figures/table1.html).")
}

message("---- Table 1 column headers ----")
for (i in seq_along(col_keys)) message("  ", col_header[i])
message("---- Table 1 (first rows) ----")
print(disp_out[1:min(.N, 16)])

message("Table 1 complete.")
message("  data/BNB_smallplus.rds, data/BNB_smallplus.csv")
message("  data/table1.csv, data/table1_long.csv")
message("  figures/table1.html", if (file.exists(file.path(FIG_DIR, "table1.png")))
        ", figures/table1.png" else "")
