library(data.table)
library(here)

OUT_DIR <- here::here("data")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

THR_ELDER <- 90  # age defining the oldest-old group, as everywhere else

# ---- data -------------------------------------------------------------------
if (!exists("BNB_small")) stop("BNB_small not found. Load it first.")
dat <- as.data.table(BNB_small)

need <- c("age", "mFIM_out", "cFIM_out")
if (length(setdiff(need, names(dat))) > 0) {
  stop("BNB_small is missing: ", paste(setdiff(need, names(dat)), collapse = ", "))
}

dat[, `:=`(age      = as.numeric(age),
           mFIM_out = as.numeric(mFIM_out),
           cFIM_out = as.numeric(cFIM_out))]

dat[, FIM_total_out := mFIM_out + cFIM_out]
dat[, age_group := factor(fifelse(age >= THR_ELDER, "elder", "young"),
                          levels = c("elder", "young"))]

n_drop <- dat[is.na(FIM_total_out) | is.na(age), .N]
if (n_drop > 0) message(sprintf("Dropped %d record(s) with a missing age or FIM total.", n_drop))
d <- dat[!is.na(FIM_total_out) & !is.na(age)]

# ---- 21.1  Welch's t test ---------------------------------------------------
tt <- t.test(FIM_total_out ~ age_group, data = d, var.equal = FALSE)

res <- data.table(
  group      = c("elder (age >= 90)", "young (age < 90)"),
  n          = d[, .N, by = age_group][order(age_group), N],
  mean       = round(as.numeric(tt$estimate), 1),
  sd         = round(d[, sd(FIM_total_out), by = age_group][order(age_group), V1], 1)
)
res[, `:=`(mean_difference = round(diff(rev(as.numeric(tt$estimate))), 2),
           ci_low  = round(tt$conf.int[1], 2),
           ci_high = round(tt$conf.int[2], 2),
           t       = round(as.numeric(tt$statistic), 3),
           df      = round(as.numeric(tt$parameter), 1),
           p_value = signif(tt$p.value, 3),
           test    = "Welch two-sample t test (two-sided)")]

print(res)
fwrite(res, file.path(OUT_DIR, "table_welch_discharge_fim_total.csv"))