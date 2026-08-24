library(data.table)
library(here)
library(ggplot2)

OUT_DIR <- here::here("data")
FIG_DIR <- here::here("figures")
for (d in c(OUT_DIR, FIG_DIR)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

LABELS_PATH <- here::here("00_labels.R")
if (!exists("relabel_vars")) {
  if (!file.exists(LABELS_PATH)) stop("00_labels.R not found: ", LABELS_PATH)
  source(LABELS_PATH)
}

CI_LEVEL   <- 0.95
THR_GOOD   <- 65
THR_SEVERE <- 26
THR_ELDER  <- 90
set.seed(20240601)

DISEASE_LONG <- c(MSD = "Musculoskeletal disease",
                  CVD = "Cerebrovascular disease",
                  DS  = "Disuse syndrome")

DISEASE_ORDER <- c("MSD", "CVD", "DS")

COL_POINT <- "#2166AC"; COL_REF <- "grey40"; COL_TEXT <- "grey15"

theme_paper <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linewidth = 0.25, colour = "grey90"),
          strip.background = element_rect(fill = "grey92", colour = NA),
          strip.text       = element_text(face = "bold", size = base_size - 0.5),
          plot.title       = element_text(face = "bold", size = base_size + 1),
          plot.subtitle    = element_text(size = base_size - 1, colour = COL_TEXT),
          plot.title.position = "plot",
          legend.position  = "none")
}

if (!exists("BNB_small")) stop("BNB_small not found. Load it first.")
dat <- as.data.table(BNB_small)

dat[, age      := as.numeric(age)]
dat[, mFIM_in  := as.numeric(mFIM_in)]
dat[, mFIM_out := as.numeric(mFIM_out)]
if (!"good"    %in% names(dat)) dat[, good    := as.integer(mFIM_out >= THR_GOOD)]
if (!"severe"  %in% names(dat)) dat[, severe  := mFIM_in <= THR_SEVERE]
if (!"elderly" %in% names(dat)) dat[, elderly := age     >= THR_ELDER]
dat[, sev_and_eld := severe & elderly]

if (!"class" %in% names(dat)) stop("BNB_small has no 'class' column.")
dat[, class := factor(class)]

lv_raw <- levels(dat$class)
lv_eng <- relabel_levels("class", lv_raw)
if (anyNA(lv_eng) || any(lv_eng == lv_raw))
  warning("Some disease levels were not translated by relabel_levels(): ",
          paste(lv_raw[lv_eng == lv_raw], collapse = ", "))
dat[, class_eng := factor(lv_eng[match(as.character(class), lv_raw)],
                          levels = intersect(DISEASE_ORDER, lv_eng))]

qc_lab <- dat[, .(n = .N), by = .(raw = as.character(class))]
qc_lab[, english := lv_eng[match(raw, lv_raw)]]
qc_lab[, percent_of_cohort := round(100 * n / nrow(dat), 1)]
qc_lab[, level_position := match(raw, lv_raw)]
qc_lab[, table1_expected := c("脳血管" = 45.5, "運動器" = 40.3, "廃用" = 14.1)[raw]]
setcolorder(qc_lab, c("level_position", "raw", "english", "n",
                      "percent_of_cohort", "table1_expected"))
fwrite(qc_lab, file.path(OUT_DIR, "qc_disease_class_labelling.csv"))

message("---- 20.0 Disease-class labelling, reconciled against Table 1 ----")
print(qc_lab)
bad <- qc_lab[is.finite(table1_expected) &
              abs(percent_of_cohort - table1_expected) > 0.5]
if (nrow(bad)) {
  warning("The class percentages do not match Table 1. Do not use these figures ",
          "until this is resolved: ", paste(bad$raw, collapse = ", "))
} else {
  message("Each Japanese level maps to the English label whose Table 1 percentage ",
          "it reproduces. Note that 脳血管 is level 1 because it is the reference ",
          "category, which is what made a positional mapping produce CVD/MSD the ",
          "wrong way round in the published Figure 2.")
}

SEV_LAB <- c("TRUE"  = "Admission motor FIM <= 26",
             "FALSE" = "Admission motor FIM >= 27")

fig2 <- rbindlist(list(
  copy(dat)[, disease_facet := "All"],
  copy(dat)[, disease_facet := as.character(class_eng)]
), fill = TRUE)
fig2 <- fig2[!is.na(disease_facet) & is.finite(age) & is.finite(mFIM_out)]
fig2[, disease_facet := factor(disease_facet, levels = c("All", DISEASE_ORDER))]
fig2[, sev_facet := factor(SEV_LAB[as.character(severe)], levels = unname(SEV_LAB))]

quad <- fig2[, {
  n <- .N
  tl <- sum(age <  THR_ELDER & mFIM_out >= THR_GOOD)
  tr <- sum(age >= THR_ELDER & mFIM_out >= THR_GOOD)
  bl <- sum(age <  THR_ELDER & mFIM_out <  THR_GOOD)
  br <- sum(age >= THR_ELDER & mFIM_out <  THR_GOOD)
  .(n = n,
    pct_young_good = 100 * tl / n, pct_old_good = 100 * tr / n,
    pct_young_poor = 100 * bl / n, pct_old_poor = 100 * br / n,
    pct_good = 100 * (tl + tr) / n)
}, by = .(disease_facet, sev_facet)]
fwrite(quad, file.path(OUT_DIR, "fig2_quadrant_percentages.csv"))
message("\n---- 20.1 Figure 2 quadrant percentages ----")
print(quad[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 1) else x)])

x_lo <- floor(min(fig2$age, na.rm = TRUE) / 10) * 10
x_hi <- ceiling(max(fig2$age, na.rm = TRUE) / 10) * 10
y_hi <- max(fig2$mFIM_out, na.rm = TRUE)

quad_lab <- rbindlist(list(
  quad[, .(disease_facet, sev_facet, x = x_lo + 0.30 * (x_hi - x_lo), y = y_hi * 1.06,
           lab = sprintf("%.1f%%", pct_young_good), vj = 1)],
  quad[, .(disease_facet, sev_facet, x = x_lo + 0.93 * (x_hi - x_lo), y = y_hi * 1.06,
           lab = sprintf("%.1f%%", pct_old_good), vj = 1)],
  quad[, .(disease_facet, sev_facet, x = x_lo + 0.30 * (x_hi - x_lo), y = -1,
           lab = sprintf("%.1f%%", pct_young_poor), vj = 0)],
  quad[, .(disease_facet, sev_facet, x = x_lo + 0.93 * (x_hi - x_lo), y = -1,
           lab = sprintf("%.1f%%", pct_old_poor), vj = 0)]
), fill = TRUE)
n_lab <- quad[, .(disease_facet, sev_facet, x = x_lo + 0.30 * (x_hi - x_lo),
                  y = -1, lab = sprintf("n = %s", format(n, big.mark = ",")))]

p2 <- ggplot(fig2, aes(x = age, y = mFIM_out)) +
  geom_hline(yintercept = THR_GOOD,  linetype = "dashed", linewidth = 0.4, colour = COL_REF) +
  geom_vline(xintercept = THR_ELDER, linetype = "dashed", linewidth = 0.4, colour = COL_REF) +
  geom_point(position = position_jitter(width = 0.45, height = 0.9, seed = 20240601),
             shape = 16, size = 0.55, alpha = 0.5, colour = COL_POINT) +
  geom_text(data = quad_lab, aes(x = x, y = y, label = lab, vjust = vj),
            inherit.aes = FALSE, size = 2.7, fontface = "bold", colour = COL_TEXT) +
  geom_text(data = n_lab, aes(x = x, y = y, label = lab),
            inherit.aes = FALSE, size = 2.7, vjust = 1.7, colour = COL_TEXT) +
  facet_grid(sev_facet ~ disease_facet, switch = "y") +
  scale_x_continuous(breaks = seq(x_lo, x_hi, by = 10), limits = c(x_lo, x_hi)) +
  scale_y_continuous(breaks = c(0, 20, 40, THR_GOOD, 80, 100),
                     limits = c(-9, y_hi * 1.12), expand = c(0, 0)) +
  labs(x = "Age at admission (years)",
       y = "Motor FIM at discharge",
       title = "Motor FIM at discharge versus age, by disease class and severity stratum",
       subtitle = paste0("Dashed lines: the good-outcome threshold (motor FIM ", THR_GOOD,
                         ") and age ", THR_ELDER,
                         ". Percentages are the share of each panel falling in each quadrant.")) +
  theme_paper() +
  theme(strip.placement = "outside")

ggsave(file.path(FIG_DIR, "fig2_outcome_by_age_disease.png"), p2,
       width = 11.5, height = 6.4, dpi = 400)
ggsave(file.path(FIG_DIR, "fig2_outcome_by_age_disease.pdf"), p2,
       width = 11.5, height = 6.4)
message("Saved figures/fig2_outcome_by_age_disease.png and .pdf")

contrast <- function(d, flag, lab, stratum = "All") {
  a <- d[get(flag) == TRUE]; b <- d[get(flag) == FALSE]
  n1 <- nrow(a); n0 <- nrow(b)
  if (n1 == 0 || n0 == 0) return(NULL)
  e1 <- sum(a$good, na.rm = TRUE); e0 <- sum(b$good, na.rm = TRUE)
  p1 <- e1 / n1; p0 <- e0 / n0
  or  <- (p1 / (1 - p1)) / (p0 / (1 - p0))
  rr  <- p1 / p0
  rd  <- p1 - p0
  se_lor <- sqrt(1 / max(e1, .5) + 1 / max(n1 - e1, .5) +
                 1 / max(e0, .5) + 1 / max(n0 - e0, .5))
  se_lrr <- sqrt(1 / max(e1, .5) - 1 / n1 + 1 / max(e0, .5) - 1 / n0)
  se_rd  <- sqrt(p1 * (1 - p1) / n1 + p0 * (1 - p0) / n0)
  z <- stats::qnorm(1 - (1 - CI_LEVEL) / 2)
  data.table(subset = lab, stratum = stratum,
             n_subset = n1, events_subset = e1, risk_subset = p1,
             n_reference = n0, events_reference = e0, risk_reference = p0,
             OR = or, OR_lcl = exp(log(or) - z * se_lor), OR_ucl = exp(log(or) + z * se_lor),
             RR = rr, RR_lcl = exp(log(rr) - z * se_lrr), RR_ucl = exp(log(rr) + z * se_lrr),
             RD = rd, RD_lcl = rd - z * se_rd, RD_ucl = rd + z * se_rd)
}

SUBSETS <- list(
  list(flag = "elderly",     lab = sprintf("Elderly (age >= %d)", THR_ELDER)),
  list(flag = "severe",      lab = sprintf("Severe (admission motor FIM <= %d)", THR_SEVERE))
)

fig3 <- rbindlist(lapply(SUBSETS, function(s) {
  overall <- contrast(dat, s$flag, s$lab, "All")
  by_cls  <- rbindlist(lapply(DISEASE_ORDER, function(cl) {
    if (!sum(dat$class_eng == cl, na.rm = TRUE)) return(NULL)
    flag_vec <- dat[[s$flag]] & (as.character(dat$class_eng) == cl)
    dat[, "__f_cl" := flag_vec]
    out <- contrast(dat, "__f_cl", s$lab, cl)
    dat[, "__f_cl" := NULL]
    out
  }), fill = TRUE)
  rbindlist(list(overall, by_cls), fill = TRUE)
}), fill = TRUE)

fig3[, subset := factor(subset, levels = vapply(SUBSETS, `[[`, character(1), "lab"))]
fig3[, is_overall := stratum == "All"]
fig3[, row_label := ifelse(is_overall, as.character(subset),
                           paste0("    ", DISEASE_LONG[stratum]))]
fig3[, stratum_ord := ifelse(is_overall, 0L, match(stratum, DISEASE_ORDER))]
setorder(fig3, subset, stratum_ord)
fig3[, row_id := .I]
fig3[, ypos := nrow(fig3) - row_id + 1L]
fig3[, count_label := sprintf("%d / %d", events_subset, n_subset)]
fig3[, or_label := sprintf("%.2f (%.2f to %.2f)", OR, OR_lcl, OR_ucl)]

fwrite(fig3[, .(subset, stratum, n_subset, events_subset,
                risk_subset = round(risk_subset, 3),
                n_reference, events_reference,
                risk_reference = round(risk_reference, 3),
                OR = round(OR, 3), OR_lcl = round(OR_lcl, 3), OR_ucl = round(OR_ucl, 3),
                RR = round(RR, 3), RD = round(RD, 3))],
       file.path(OUT_DIR, "fig3_subset_odds.csv"))
message("\n---- 20.2 Figure 3 odds ratios (these must equal Supplementary Table S4) ----")
print(fig3[, .(subset, stratum, count_label, or_label)])

x_min <- min(fig3$OR_lcl, na.rm = TRUE)
x_max <- max(fig3$OR_ucl, na.rm = TRUE)
brks  <- c(0.01, 0.03125, 0.0625, 0.125, 0.25, 0.5, 1, 2, 4)
brks  <- brks[brks >= x_min / 1.5 & brks <= x_max * 1.5]

Y_LIM  <- c(0.4, nrow(fig3) + 1.6)
HDR_Y  <- nrow(fig3) + 1.0

FIG3_TITLE <- "Odds of a good functional outcome in each subset relative to the rest of the cohort"
FIG3_SUB   <- paste0("Good outcome = motor FIM ", THR_GOOD,
                     " or more at discharge. ")

X_TITLE <- "Odds ratio for a good functional outcome (log scale)"

theme_col <- function() {
  theme_paper(11) +
    theme(panel.border     = element_blank(),
          panel.background = element_blank(),
          panel.grid       = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.y      = element_blank(),
          axis.ticks.y     = element_blank(),
          axis.text.x      = element_text(colour = NA),
          axis.title.x     = element_text(colour = NA),
          axis.ticks.x     = element_line(colour = NA),
          plot.margin      = margin(0, 2, 0, 2))
}

forest_panel <- function(show_y_labels, y_labels = fig3$row_label) {
  p <- ggplot(fig3, aes(x = OR, y = ypos)) +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.4, colour = COL_REF) +
    geom_errorbarh(aes(xmin = OR_lcl, xmax = OR_ucl), height = 0.18,
                   linewidth = 0.45, colour = COL_POINT) +
    geom_point(aes(size = n_subset), shape = 15, colour = COL_POINT) +
    scale_size_area(max_size = 3.4, guide = "none") +
    scale_x_log10(breaks = brks,
                  labels = function(v) format(v, drop0trailing = TRUE, scientific = FALSE),
                  limits = c(x_min / 1.6, x_max * 1.6), expand = c(0, 0)) +
    scale_y_continuous(limits = Y_LIM, expand = c(0, 0),
                       breaks = fig3$ypos, labels = y_labels) +
    labs(x = X_TITLE, y = NULL) +
    theme_paper() +
    theme(panel.grid.major.y = element_blank(),
          axis.ticks.y = element_blank())
  if (show_y_labels) p + theme(axis.text.y = element_text(hjust = 0, colour = COL_TEXT))
  else               p + theme(axis.text.y = element_blank())
}

text_panel <- function(value_col, header, hjust = 0.5, x_at = 0, x_lim = c(-1, 1)) {
  bold  <- fig3[is_overall  == TRUE]
  plain <- fig3[is_overall  == FALSE]
  ggplot() +
    geom_text(data = bold,  aes(x = x_at, y = ypos, label = .data[[value_col]]),
              hjust = hjust, size = 2.95, colour = COL_TEXT, fontface = "bold") +
    geom_text(data = plain, aes(x = x_at, y = ypos, label = .data[[value_col]]),
              hjust = hjust, size = 2.95, colour = COL_TEXT, fontface = "plain") +
    annotate("text", x = x_at, y = HDR_Y, label = header, hjust = hjust,
             size = 2.95, fontface = "bold", colour = COL_TEXT) +
    scale_x_continuous(limits = x_lim, expand = c(0, 0), breaks = x_at, labels = " ") +
    scale_y_continuous(limits = Y_LIM, expand = c(0, 0)) +
    labs(x = X_TITLE, y = NULL) +
    theme_col()
}

if (requireNamespace("patchwork", quietly = TRUE)) {
  lab_panel <- text_panel("row_label", "Subgroup", hjust = 0, x_at = 0, x_lim = c(0, 1))
  cnt_panel <- text_panel("count_label", "Good outcome\nn / N")
  est_panel <- text_panel("or_label",    "Odds ratio\n(95% CI)")

  p3 <- patchwork::wrap_plots(lab_panel, cnt_panel, forest_panel(FALSE), est_panel,
                              nrow = 1, widths = c(0.30, 0.12, 0.36, 0.22)) +
    patchwork::plot_annotation(
      title = FIG3_TITLE, subtitle = FIG3_SUB,
      theme = theme(plot.title = element_text(face = "bold", size = 12),
                    plot.subtitle = element_text(size = 10, colour = COL_TEXT),
                    plot.title.position = "plot"))
  fig3_w <- 11.0
} else {
  message("patchwork is not installed; drawing Figure 3 as a single panel. ",
          "install.packages('patchwork') gives the aligned-column layout.")
  fig3[, row_label_n := sprintf("%s   (%s)", row_label, count_label)]
  p3 <- forest_panel(TRUE, y_labels = fig3$row_label_n) +
    geom_text(aes(x = Inf, label = or_label), hjust = -0.12,
              size = 2.95, colour = COL_TEXT) +
    annotate("text", x = Inf, y = HDR_Y, label = "Odds ratio (95% CI)",
             hjust = -0.05, size = 2.95, fontface = "bold", colour = COL_TEXT) +
    coord_cartesian(clip = "off") +
    labs(title = FIG3_TITLE, subtitle = FIG3_SUB) +
    theme(plot.margin = margin(6, 96, 6, 6))
  fig3_w <- 11.5
}

ggsave(file.path(FIG_DIR, "fig3_subset_forest.png"), p3,
       width = fig3_w, height = 5.8, dpi = 400)
ggsave(file.path(FIG_DIR, "fig3_subset_forest.pdf"), p3,
       width = fig3_w, height = 5.8)
message("Saved figures/fig3_subset_forest.png and .pdf")

message("\n==== Section 20 complete ====")
message("  figures/fig2_outcome_by_age_disease.png, .pdf")
message("  figures/fig3_subset_forest.png, .pdf")
message("  data/fig2_quadrant_percentages.csv")
message("  data/fig3_subset_odds.csv")
message("  data/qc_disease_class_labelling.csv")
message("")
message("Before replacing the images in the manuscript, check two things:")
message("  1. qc_disease_class_labelling.csv reproduces Table 1 (CVD 45.5, MSD 40.3, DS 14.1).")
message("     The published Figure 2 does not: its MSD column holds the CVD patients.")
message("  2. fig3_subset_odds.csv equals Supplementary Table S4 row for row.")
message("The Results text needs no change on account of Figure 2: the only disease-specific")
message("sentence about that figure concerns disuse syndrome, which was labelled correctly.")
