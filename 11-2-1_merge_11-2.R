library(here)

need_pkg <- function(p) if (!requireNamespace(p, quietly = TRUE))
  stop(sprintf("パッケージ '%s' が必要です。install.packages('%s')", p, p))
need_pkg("magick")
library(magick)

LABELS_PATH <- here::here("00_labels.R")
if (!exists("relabel_pipeline")) {
  if (!file.exists(LABELS_PATH)) stop("00_labels.R が見つかりません: ", LABELS_PATH)
  source(LABELS_PATH)
}

FIG_DIR <- here::here("figures")

BG        <- "white"
GAP_X     <- 24L
HEADER_H  <- 64L
TITLE_H   <- 84L
HEADER_SZ <- 34
TITLE_SZ  <- 52

GROUPS <- list(
  list(title = "Severe subset (mFIM_in <= 26)",
       pipes = c("severe_A",  "severe_B"),
       out   = "fig4_interaction_severe_combined.png"),
  list(title = "Elderly subset (age >= 90)",
       pipes = c("elderly_A", "elderly_B"),
       out   = "fig4_interaction_elderly_combined.png")
)

list_pngs <- function(nm) {
  pat <- file.path(FIG_DIR, sprintf("fig4_interaction_%s_*.png", nm))
  sort(Sys.glob(pat))
}

pad_width <- function(img, w) {
  h <- magick::image_info(img)$height
  magick::image_extent(img, geometry = sprintf("%dx%d", w, h),
                       gravity = "west", color = BG)
}
pad_height <- function(img, h) {
  w <- magick::image_info(img)$width
  magick::image_extent(img, geometry = sprintf("%dx%d", w, h),
                       gravity = "north", color = BG)
}
label_strip <- function(text, w, ht, size) {
  magick::image_annotate(
    magick::image_blank(width = w, height = ht, color = BG),
    text, gravity = "center", size = size, color = "black", weight = 700)
}

build_column <- function(nm) {
  files <- list_pngs(nm)
  if (length(files) == 0) {
    message(sprintf("  [%s] 交互作用 PNG が見つからない（交互作用なし or §11-2 未実行）。", nm))
    return(NULL)
  }
  imgs <- lapply(files, magick::image_read)
  w    <- max(vapply(imgs, function(im) magick::image_info(im)$width, integer(1)))
  imgs <- lapply(imgs, pad_width, w = w)
  body <- magick::image_append(do.call(c, imgs), stack = TRUE)
  hdr  <- label_strip(relabel_pipeline(nm), w, HEADER_H, HEADER_SZ)
  message(sprintf("  [%s] %d 図を結合（列幅=%dpx）。", nm, length(files), w))
  magick::image_append(c(hdr, body), stack = TRUE)
}

combine_group <- function(grp) {
  message(sprintf("== %s ==", grp$title))
  cols <- Filter(Negate(is.null), lapply(grp$pipes, build_column))
  if (length(cols) == 0) {
    message("  → 図示すべき交互作用が無いためスキップ。"); return(invisible(NULL))
  }
  Hmax <- max(vapply(cols, function(im) magick::image_info(im)$height, integer(1)))
  cols <- lapply(cols, pad_height, h = Hmax)
  if (length(cols) > 1L) {
    gap  <- magick::image_blank(width = GAP_X, height = Hmax, color = BG)
    seq2 <- vector("list", 2L * length(cols) - 1L)
    seq2[seq(1, length(seq2), by = 2L)] <- cols
    if (length(cols) > 1L) seq2[seq(2, length(seq2), by = 2L)] <- rep(list(gap), length(cols) - 1L)
    row <- magick::image_append(do.call(c, seq2), stack = FALSE)
  } else row <- cols[[1]]

  W     <- magick::image_info(row)$width
  title <- label_strip(grp$title, W, TITLE_H, TITLE_SZ)
  canvas <- magick::image_append(c(title, row), stack = TRUE)

  outfile <- file.path(FIG_DIR, grp$out)
  magick::image_write(canvas, outfile, format = "png")
  message(sprintf("  → 保存: %s（%dx%d px）",
                  outfile, magick::image_info(canvas)$width, magick::image_info(canvas)$height))
  invisible(outfile)
}

invisible(lapply(GROUPS, combine_group))

message("\n§11-2-1 完了。重症(A|B)・高齢(A|B)の交互作用プロットを各 1 枚に結合し、",
        normalizePath(FIG_DIR), " に保存しました（fig4_interaction_severe_combined.png /",
        " fig4_interaction_elderly_combined.png）。各パイプラインを 1 列に縦積みし、列ヘッダーと",
        " グループ見出しを付与。鍵の数がパイプライン間で異なっても自動対応します。")
