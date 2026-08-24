VAR_LABELS <- c(
  age            = "Age",
  BBS_in         = "BBS",
  cFIM_in        = "cFIM",
  class          = "Disease",
  GP_better_in   = "Grip Strength",
  his_car        = "History of Cancer",
  his_CNS        = "History of diseases of CNS",
  JCS_in         = "JCS",
  JCS_bin        = "JCS",
  mFIM_in        = "mFIM",
  MMSE_in        = "MMSE",
  sex            = "Sex",
  STEF_better_in = "STEF",
  support_in     = "Care",
  walk_speed_in  = "Walk Speed",
  X6MD_in        = "6MWD",
  mFIM_out       = "mFIM at Discharge",
  good           = "Good Outcome"
)

LEVEL_LABELS <- c(
  "脳血管"   = "CVD",
  "運動器"   = "MSD",
  "廃用"     = "DS",
  "有"       = "Yes",
  "無"       = "No",
  "女"       = "Female",
  "男"       = "Male",
  "あり"     = "Needed",
  "なし"     = "Independent",
  "clear"    = "Clear (0-3)",
  "impaired" = "Impaired (>=10)"
)

PIPE_LABELS <- c(
  severe_A  = "Severe / System A (worst-value)",
  severe_B  = "Severe / System B (missing-category)",
  elderly_A = "Elderly / System A (worst-value)",
  elderly_B = "Elderly / System B (missing-category)"
)

relabel_levels <- function(var, lv) {
  out <- ifelse(lv %in% names(LEVEL_LABELS), LEVEL_LABELS[lv], lv)
  unname(out)
}

relabel_vars <- function(x) {
  vapply(x, function(s) {
    if (grepl(":", s, fixed = TRUE)) {
      parts <- strsplit(s, ":", fixed = TRUE)[[1]]
      base  <- parts[1]; lev <- paste(parts[-1], collapse = ":")
      blab  <- if (base %in% names(VAR_LABELS)) unname(VAR_LABELS[base]) else base
      llab  <- if (lev %in% names(LEVEL_LABELS)) unname(LEVEL_LABELS[lev]) else lev
      return(paste0(blab, ": ", llab))
    }
    if (grepl("_measurable$", s)) {
      base <- sub("_measurable$", "", s)
      blab <- if (base %in% names(VAR_LABELS)) unname(VAR_LABELS[base]) else base
      return(paste0(blab, " (measurable)"))
    }
    if (grepl("_c$", s)) {
      base <- sub("_c$", "", s)
      if (base %in% names(VAR_LABELS)) return(unname(VAR_LABELS[base]))
    }
    if (s %in% names(VAR_LABELS)) return(unname(VAR_LABELS[s]))
    s
  }, character(1), USE.NAMES = FALSE)
}

relabel_pipeline <- function(label) {
  if (label %in% names(PIPE_LABELS)) unname(PIPE_LABELS[label]) else label
}
