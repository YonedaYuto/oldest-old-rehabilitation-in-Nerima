library(data.table)
library(here)
library(ggplot2)

OUT_DIR <- here::here("data")
FIG_DIR <- here::here("figures")
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

LABELS_PATH <- here::here("00_labels.R")
if (!exists("relabel_vars")) {
  if (!file.exists(LABELS_PATH)) stop("00_labels.R not found: ", LABELS_PATH)
  source(LABELS_PATH)
}

SOURCE_ALL    <- "Dataset"
STATUS_COL    <- "disposition"
STATUS_EXCL   <- c("医療機関", "病院・診療所へ転院", "終了（死亡等）")
OUTCOME_COL   <- "mFIM_out"

THR_GOOD <- 65; THR_SEVERE <- 26; THR_ELDER <- 90

if (!exists(SOURCE_ALL))
  stop(sprintf(paste0("'%s' not found. Run load(\"…/BNB_all.rda\") first, then ",
                      "check the restored name with ls() and set SOURCE_ALL to it."),
               SOURCE_ALL))
if (!exists("BNB_small")) stop("BNB_small (the analysed cohort) not found. Load it first.")

src <- get(SOURCE_ALL)

allx <- copy(as.data.table(src))
ana  <- copy(as.data.table(BNB_small))

only_all   <- setdiff(names(allx), names(ana))
only_small <- setdiff(names(ana), names(allx))
if (length(only_all) || length(only_small))
  message(sprintf("Column differences - only in %s: [%s]; only in BNB_small: [%s]",
                  SOURCE_ALL, paste(only_all, collapse = ", "),
                  paste(only_small, collapse = ", ")))
message(sprintf("%s: %d rows x %d columns; BNB_small: %d rows x %d columns.",
                SOURCE_ALL, nrow(allx), ncol(allx), nrow(ana), ncol(ana)))

harmonise <- function(d) {
  if ("JCS_in" %in% names(d) && !"JCS_bin" %in% names(d)) {
    d[, .jcs_num := suppressWarnings(as.numeric(as.character(JCS_in)))]
    d[.jcs_num == 300, .jcs_num := 30]
    d[, JCS_bin := factor(fcase(is.na(.jcs_num),                  NA_character_,
                                .jcs_num %in% c(0, 1, 2, 3),      "clear",
                                default =                         "impaired"),
                          levels = c("clear", "impaired"))]
    d[, .jcs_num := NULL]
  }
  d[]
}
allx <- harmonise(allx); ana <- harmonise(ana)

if (!STATUS_COL %in% names(allx))
  stop(sprintf(paste0("Column '%s' is not in %s; it carries the primary exclusion ",
                      "criterion, so Section 14 cannot proceed without it.\n",
                      "  Columns present: %s"),
               STATUS_COL, SOURCE_ALL, paste(names(allx), collapse = ", ")))
if (!OUTCOME_COL %in% names(allx))
  stop(sprintf("Column '%s' is not in %s; it carries the second exclusion criterion.",
               OUTCOME_COL, SOURCE_ALL))

allx[, .disp := trimws(as.character(get(STATUS_COL)))]

disp_tab <- allx[, .N, by = .disp][order(-N)]
disp_tab[, excludes := .disp %in% STATUS_EXCL]
message("---- disposition values in ", SOURCE_ALL, " ----")
print(disp_tab)

hit <- intersect(STATUS_EXCL, disp_tab$.disp)
if (!length(hit))
  stop(sprintf(paste0("None of the excluding dispositions [%s] occurs in '%s'. ",
                      "The labels differ from the data - check for full-width vs ",
                      "half-width characters against the values printed above."),
               paste(STATUS_EXCL, collapse = " / "), STATUS_COL))
if (length(hit) < length(STATUS_EXCL))
  warning(sprintf("Excluding dispositions not observed in the data: %s",
                  paste(setdiff(STATUS_EXCL, hit), collapse = " / ")))

allx[, .excl_disp := .disp %in% STATUS_EXCL]
allx[, .excl_out  := !.excl_disp &
       is.na(suppressWarnings(as.numeric(as.character(get(OUTCOME_COL)))))]
allx[, .excluded  := .excl_disp | .excl_out]
allx[, .group := fifelse(.excluded, "Excluded", "Included")]

n_inc_flag <- allx[.group == "Included", .N]; n_exc_flag <- allx[.group == "Excluded", .N]
message(sprintf(paste0("Excluded %d = %d by disposition + %d with %s missing; ",
                       "Included %d (of %d admissions)."),
                n_exc_flag, allx[, sum(.excl_disp)], allx[, sum(.excl_out)],
                OUTCOME_COL, n_inc_flag, nrow(allx)))
if (n_inc_flag != nrow(ana))
  warning(sprintf(paste0("Rows flagged as included (%d) differ from BNB_small (%d). ",
                         "Check STATUS_EXCL against the disposition values printed above."),
                  n_inc_flag, nrow(ana)))

reason_tab <- rbind(
  allx[.excl_disp == TRUE, .(n = .N), by = .(reason = .disp)],
  allx[.excl_out  == TRUE, .(reason = sprintf("%s missing (outcome unevaluable)",
                                              OUTCOME_COL), n = .N)][n > 0])[order(-n)]
if (nrow(reason_tab)) {
  reason_tab[, pct_of_all_admissions := round(100 * n / nrow(allx), 1)]
  fwrite(reason_tab, file.path(OUT_DIR, "table_s1_exclusion_reasons.csv"))
  message("---- reason for exclusion ----"); print(reason_tab)
}

allx[, (STATUS_COL) := NULL]

CONT_VARS <- intersect(c("age", "mFIM_in", "cFIM_in", "MMSE_in", "BBS_in",
                         "GP_better_in", "STEF_better_in", "walk_speed_in",
                         "X6MD_in", "term"), names(allx))
CAT_VARS  <- intersect(c("sex", "class", "his_car", "his_CNS", "support_in",
                         "JCS_bin"), names(allx))
missed <- setdiff(names(allx), c(CONT_VARS, CAT_VARS, OUTCOME_COL, "JCS_in",
                                 ".disp", ".excl_disp", ".excl_out",
                                 ".excluded", ".group"))
if (length(missed))
  message("Columns present but not compared: ", paste(missed, collapse = ", "))

fmt_med <- function(x) {
  q <- stats::quantile(x, c(0.5, 0.25, 0.75), na.rm = TRUE)
  sprintf("%s [%s, %s]", signif(q[1], 3), signif(q[2], 3), signif(q[3], 3))
}
smd_cont <- function(x, g) {
  a <- x[g == "Included"]; b <- x[g == "Excluded"]
  s <- sqrt((stats::var(a, na.rm = TRUE) + stats::var(b, na.rm = TRUE)) / 2)
  if (!is.finite(s) || s == 0) return(NA_real_)
  (mean(a, na.rm = TRUE) - mean(b, na.rm = TRUE)) / s
}
smd_prop <- function(p1, p2) {
  s <- sqrt((p1 * (1 - p1) + p2 * (1 - p2)) / 2)
  if (!is.finite(s) || s == 0) return(NA_real_)
  (p1 - p2) / s
}

rows <- list()
rows[[length(rows) + 1]] <- data.table(
  Characteristic = "n", Included = as.character(allx[.group == "Included", .N]),
  Excluded = as.character(allx[.group == "Excluded", .N]), SMD = NA_character_)

for (v in CONT_VARS) {
  x <- suppressWarnings(as.numeric(allx[[v]])); g <- allx$.group
  rows[[length(rows) + 1]] <- data.table(
    Characteristic = sprintf("%s, median [Q1, Q3]", relabel_vars(v)),
    Included = fmt_med(x[g == "Included"]),
    Excluded = fmt_med(x[g == "Excluded"]),
    SMD = sprintf("%.2f", smd_cont(x, g)))
  rows[[length(rows) + 1]] <- data.table(
    Characteristic = "    Missing, n (%)",
    Included = sprintf("%d (%.1f%%)", sum(is.na(x[g == "Included"])),
                       100 * mean(is.na(x[g == "Included"]))),
    Excluded = sprintf("%d (%.1f%%)", sum(is.na(x[g == "Excluded"])),
                       100 * mean(is.na(x[g == "Excluded"]))),
    SMD = NA_character_)
}

for (v in CAT_VARS) {
  f <- factor(allx[[v]]); g <- allx$.group
  for (lv in levels(f)) {
    p1 <- mean(f[g == "Included"] == lv, na.rm = TRUE)
    p2 <- mean(f[g == "Excluded"] == lv, na.rm = TRUE)
    rows[[length(rows) + 1]] <- data.table(
      Characteristic = sprintf("%s = %s, n (%%)", relabel_vars(v), lv),
      Included = sprintf("%d (%.1f%%)", sum(f[g == "Included"] == lv, na.rm = TRUE), 100 * p1),
      Excluded = sprintf("%d (%.1f%%)", sum(f[g == "Excluded"] == lv, na.rm = TRUE), 100 * p2),
      SMD = sprintf("%.2f", smd_prop(p1, p2)))
  }
}

tab_s1 <- rbindlist(rows, fill = TRUE)
saveRDS(tab_s1, file.path(OUT_DIR, "table_s1_excluded_vs_included.rds"))
fwrite(tab_s1, file.path(OUT_DIR, "table_s1_excluded_vs_included.csv"))
message("---- Supplementary Table S1 (first rows) ----"); print(head(tab_s1, 15))

if (all(c("mFIM_in", "age") %in% names(allx))) {
  allx[, severe  := as.numeric(mFIM_in) <= THR_SEVERE]
  allx[, elderly := as.numeric(age)     >= THR_ELDER]
  allx[, sev_and_eld := severe & elderly]
  good_rate <- function(d) {
    if (!OUTCOME_COL %in% names(d)) return(NA_real_)
    y <- as.numeric(d[[OUTCOME_COL]])
    if (all(is.na(y))) return(NA_real_)
    100 * mean(y >= THR_GOOD, na.rm = TRUE)
  }
  by_str <- rbindlist(lapply(c("severe", "elderly", "sev_and_eld"), function(s) {
    sub <- allx[get(s) == TRUE]
    inc <- sub[.group == "Included"]; exc <- sub[.group == "Excluded"]
    data.table(stratum = s,
               n_total = nrow(sub),
               n_included = nrow(inc), n_excluded = nrow(exc),
               pct_excluded = round(100 * nrow(exc) / max(nrow(sub), 1), 1),
               age_med_included = stats::median(as.numeric(inc$age), na.rm = TRUE),
               age_med_excluded = stats::median(as.numeric(exc$age), na.rm = TRUE),
               mFIMin_med_included = stats::median(as.numeric(inc$mFIM_in), na.rm = TRUE),
               mFIMin_med_excluded = stats::median(as.numeric(exc$mFIM_in), na.rm = TRUE),
               good_pct_included = round(good_rate(inc), 1),
               good_pct_worst_case = round(good_rate(inc) * nrow(inc) /
                                             max(nrow(sub), 1), 1))
  }))
  by_str <- rbind(by_str, data.table(
    stratum = "whole cohort", n_total = nrow(allx),
    n_included = allx[.group == "Included", .N],
    n_excluded = allx[.group == "Excluded", .N],
    pct_excluded = round(100 * allx[.group == "Excluded", .N] / nrow(allx), 1),
    age_med_included = stats::median(as.numeric(allx[.group == "Included"]$age), na.rm = TRUE),
    age_med_excluded = stats::median(as.numeric(allx[.group == "Excluded"]$age), na.rm = TRUE),
    mFIMin_med_included = stats::median(as.numeric(allx[.group == "Included"]$mFIM_in), na.rm = TRUE),
    mFIMin_med_excluded = stats::median(as.numeric(allx[.group == "Excluded"]$mFIM_in), na.rm = TRUE),
    good_pct_included = round(good_rate(allx[.group == "Included"]), 1),
    good_pct_worst_case = round(good_rate(allx[.group == "Included"]) *
                                  allx[.group == "Included", .N] / nrow(allx), 1)),
    fill = TRUE)
  fwrite(by_str, file.path(OUT_DIR, "table_s1_by_stratum.csv"))
  message("---- Exclusion by stratum ----"); print(by_str)
}

n_inc  <- allx[.group == "Included", .N]
n_exc  <- allx[.group == "Excluded", .N]
good_n <- if ("mFIM_out" %in% names(ana)) {
  sum(as.numeric(ana$mFIM_out) >= THR_GOOD, na.rm = TRUE)
} else if ("good" %in% names(ana)) {
  sum(as.integer(ana$good), na.rm = TRUE)
} else NA_integer_
if (is.na(good_n)) stop("Neither mFIM_out nor good is present in BNB_small.")

tip <- data.table(assumed_good_rate_in_excluded = seq(0, 1, by = 0.05))
tip[, overall_good_pct := 100 * (good_n + assumed_good_rate_in_excluded * n_exc) / (n_inc + n_exc)]
tip[, note := fifelse(assumed_good_rate_in_excluded == 0,
                      "worst case (all excluded assumed not to reach the outcome)",
                      fifelse(abs(assumed_good_rate_in_excluded - good_n / n_inc) < 0.026,
                              "missing-at-random-like (same rate as the analysed cohort)", ""))]
fwrite(tip, file.path(OUT_DIR, "table_s1_tipping_point.csv"))
message(sprintf(paste0("Good-outcome proportion: %.1f%% as analysed (%d/%d); ",
                       "%.1f%% if none of the %d excluded patients had reached it."),
                100 * good_n / n_inc, good_n, n_inc, tip[1]$overall_good_pct, n_exc))

g <- ggplot(tip, aes(x = 100 * assumed_good_rate_in_excluded, y = overall_good_pct)) +
  geom_line(colour = "#2166AC", linewidth = 0.9) +
  geom_hline(yintercept = 100 * good_n / n_inc, linetype = 2, colour = "grey40") +
  annotate("text", x = 5, y = 100 * good_n / n_inc + 1.5, hjust = 0, size = 3,
           label = sprintf("as analysed: %.1f%%", 100 * good_n / n_inc)) +
  labs(x = "Assumed good-outcome rate among the excluded patients (%)",
       y = "Good-outcome proportion over all admissions (%)",
       title = "Sensitivity of the good-outcome proportion to the excluded patients") +
  theme_bw(base_size = 11)
ggsave(file.path(FIG_DIR, "fig_s_tipping_point.png"), g, width = 7, height = 4.5, dpi = 300)

message("Section 14 complete.")
message("  data/table_s1_excluded_vs_included.csv, table_s1_exclusion_reasons.csv,")
message("  table_s1_by_stratum.csv, table_s1_tipping_point.csv")
message("  figures/fig_s_tipping_point.png")
