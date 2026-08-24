library(here)

`%||%` <- function(a, b) if (is.null(a) || !length(a) || all(is.na(a))) b else a

FORCE_LOAD_ANALYSIS_PKGS <- TRUE

SESSION_DIR <- here::here("session")
dir.create(SESSION_DIR, showWarnings = FALSE, recursive = TRUE)

KEY_PKGS <- c("smcfcs", "logistf", "rms", "mice", "Hmisc", "data.table",
              "ggplot2", "detectseparation", "brglm2", "here",
              "officer", "flextable", "pROC", "boot")

COPYRIGHT_HOLDER <- "Yuto Yoneda"
COPYRIGHT_YEAR   <- format(Sys.Date(), "%Y")

if (FORCE_LOAD_ANALYSIS_PKGS) {
  for (p in KEY_PKGS)
    if (requireNamespace(p, quietly = TRUE)) suppressMessages(library(p, character.only = TRUE))
}

si <- utils::sessionInfo()

si_path <- file.path(SESSION_DIR, "sessionInfo.txt")
con <- file(si_path, open = "wt", encoding = "UTF-8")
writeLines(c(
  "# ---------------------------------------------------------------------------",
  "# sessionInfo() for the analysis reported in this manuscript",
  sprintf("# recorded on %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "# ---------------------------------------------------------------------------",
  ""), con)
writeLines(utils::capture.output(print(si)), con)
close(con)
cat("書き出しました: session/sessionInfo.txt\n")

collect <- function(lst, origin) {
  if (!length(lst)) return(NULL)
  data.frame(
    package = vapply(lst, function(x) x$Package, character(1)),
    version = vapply(lst, function(x) x$Version, character(1)),
    origin  = origin,
    stringsAsFactors = FALSE
  )
}

pkg_df <- do.call(rbind, list(
  data.frame(package = "R", version = paste(si$R.version$major, si$R.version$minor, sep = "."),
             origin = "base R", stringsAsFactors = FALSE),
  collect(si$otherPkgs,  "attached"),
  collect(si$loadedOnly, "loaded via a namespace")
))
pkg_df <- pkg_df[!duplicated(pkg_df$package), ]
pkg_df$key_analysis_package <- pkg_df$package %in% KEY_PKGS
pkg_df <- pkg_df[order(!pkg_df$key_analysis_package, pkg_df$package), ]

utils::write.csv(pkg_df, file.path(SESSION_DIR, "package_versions.csv"),
                 row.names = FALSE, fileEncoding = "UTF-8")
cat("書き出しました: session/package_versions.csv (", nrow(pkg_df), " rows)\n", sep = "")

missing_key <- setdiff(intersect(KEY_PKGS, rownames(utils::installed.packages())), pkg_df$package)
if (length(missing_key))
  message("[注] 未読み込みのため記録されなかった主要パッケージ: ",
          paste(missing_key, collapse = ", "),
          "  → 解析と同じセッションで実行するか FORCE_LOAD_ANALYSIS_PKGS <- TRUE にしてください。")

if (requireNamespace("renv", quietly = TRUE)) {
  ok <- tryCatch({
    renv::snapshot(project = here::here(), prompt = FALSE, force = TRUE)
    TRUE
  }, error = function(e) { message("renv::snapshot に失敗: ", conditionMessage(e)); FALSE })
  if (ok) {
    lock <- here::here("renv.lock")
    if (file.exists(lock)) {
      file.copy(lock, file.path(SESSION_DIR, "renv.lock"), overwrite = TRUE)
      cat("書き出しました: session/renv.lock\n")
    }
  }
} else {
  message("[注] renv が未導入のため renv.lock は作成していません。",
          " install.packages('renv') 後に再実行すると、",
          "第三者がパッケージ群を厳密に復元できるようになります。")
}

mit <- c(
  sprintf("MIT License"),
  "",
  sprintf("Copyright (c) %s %s", COPYRIGHT_YEAR, COPYRIGHT_HOLDER),
  "",
  "Permission is hereby granted, free of charge, to any person obtaining a copy",
  "of this software and associated documentation files (the \"Software\"), to deal",
  "in the Software without restriction, including without limitation the rights",
  "to use, copy, modify, merge, publish, distribute, sublicense, and/or sell",
  "copies of the Software, and to permit persons to whom the Software is",
  "furnished to do so, subject to the following conditions:",
  "",
  "The above copyright notice and this permission notice shall be included in all",
  "copies or substantial portions of the Software.",
  "",
  "THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR",
  "IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,",
  "FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE",
  "AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER",
  "LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,",
  "OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE",
  "SOFTWARE."
)
writeLines(mit, file.path(here::here(), "LICENSE"), useBytes = TRUE)

docs <- c(
  "Licence for the documentation, figures, tables and reporting checklists",
  "in this repository (i.e. everything except the source code)",
  "",
  sprintf("Copyright (c) %s %s", COPYRIGHT_YEAR, COPYRIGHT_HOLDER),
  "",
  "These materials are licensed under the Creative Commons Attribution 4.0",
  "International Licence (CC BY 4.0). You are free to share and adapt them for",
  "any purpose, including commercially, provided that appropriate credit is",
  "given, a link to the licence is provided, and any changes are indicated.",
  "",
  "Full licence text: https://creativecommons.org/licenses/by/4.0/legalcode",
  "",
  "NOTE ON THE INDIVIDUAL PATIENT DATA",
  "The individual patient records underlying this analysis are NOT covered by",
  "either licence and are not redistributed. They are held under institutional",
  "and privacy restrictions and may be available from the corresponding author",
  "on reasonable request and with the necessary approvals."
)
writeLines(docs, file.path(here::here(), "LICENSE-docs.txt"), useBytes = TRUE)
cat("書き出しました: LICENSE (MIT), LICENSE-docs.txt (CC BY 4.0)\n")

key_line <- paste(
  vapply(intersect(KEY_PKGS, pkg_df$package), function(p)
    sprintf("%s %s", p, pkg_df$version[match(p, pkg_df$package)]), character(1)),
  collapse = ", ")

stmt <- c(
  sprintf("The analysis was carried out in R %s on %s.",
          paste(si$R.version$major, si$R.version$minor, sep = "."),
          si$running %||% "the platform recorded in sessionInfo.txt"),
  sprintf("The packages used, with the versions actually loaded, were: %s.", key_line),
  "The complete output of sessionInfo(), the full package version table and, where",
  "available, the renv lockfile that pins every dependency are provided as",
  "Supplementary Material (session/sessionInfo.txt, session/package_versions.csv,",
  "session/renv.lock; script 17).",
  "The analysis code is released under the MIT licence and the supplementary",
  "documents, figures, tables and reporting checklists under the Creative Commons",
  "Attribution 4.0 International licence (CC BY 4.0), so that both may be reused",
  "and adapted with attribution. The individual patient data are not redistributed."
)
writeLines(stmt, file.path(SESSION_DIR, "reproducibility_statement.txt"), useBytes = TRUE)

cat("\n---------------- 原稿貼り付け用 ----------------\n")
cat(paste(stmt, collapse = "\n"), "\n")
cat("------------------------------------------------\n")
cat("\n完了: session/ 以下と LICENSE / LICENSE-docs.txt を確認してください。\n")
