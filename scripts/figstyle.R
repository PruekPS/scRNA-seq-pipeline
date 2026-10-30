# Shared Nature-style figure settings for the R steps (07, 11); the Python twin is
# figstyle.py. Arial or Helvetica at 5-7 pt, 0.5 pt lines, 89/120/183 mm widths, no
# titles inside figures, editable text in the PDF and SVG. Colours are colour-blind safe.
#
# Not a pipeline step: nothing to run. Other R scripts source() this file.

NATURE_SINGLE_MM <- 89
NATURE_DOUBLE_MM <- 183
NATURE_INK <- "#1a1a1a"
NATURE_FONT <- "sans"  # R maps this to Helvetica in PDF and to the system sans font in PNG/TIFF

# ggplot line widths are in mm: 0.5 pt = 0.176 mm
theme_nature <- function(base_size = 6) {
  ggplot2::theme_classic(base_size = base_size, base_family = NATURE_FONT) +
    ggplot2::theme(
      text = ggplot2::element_text(colour = NATURE_INK, size = base_size),
      axis.text = ggplot2::element_text(colour = NATURE_INK, size = 5.5),
      axis.title = ggplot2::element_text(colour = NATURE_INK, size = 6),
      axis.line = ggplot2::element_line(colour = NATURE_INK, linewidth = 0.176),
      axis.ticks = ggplot2::element_line(colour = NATURE_INK, linewidth = 0.176),
      axis.ticks.length = grid::unit(2, "pt"),
      legend.text = ggplot2::element_text(size = 5.5),
      legend.title = ggplot2::element_text(size = 5.5),
      legend.key.size = grid::unit(6, "pt"),
      legend.background = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(2, 2, 2, 2, "pt")
    )
}

# Write a ggplot as PDF + SVG + 600 dpi PNG + TIFF, honouring FIG_FORMATS (see figstyle.py).
# `path` has no extension. SVG is skipped with a note if svglite is not installed.
save_nature <- function(plot, path, width_mm = NATURE_SINGLE_MM, height_mm = 70) {
  formats <- trimws(strsplit(Sys.getenv("FIG_FORMATS", "pdf,svg,png,tiff"), ",", fixed = TRUE)[[1]])
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  for (ext in formats) {
    if (ext == "svg" && !requireNamespace("svglite", quietly = TRUE)) {
      cat("  note: svglite is not installed, skipping", paste0(path, ".svg"), "\n")
      next
    }
    args <- list(filename = paste0(path, ".", ext), plot = plot, width = width_mm, height = height_mm,
                 units = "mm", dpi = 600)
    if (ext == "tiff") args$compression <- "lzw"
    if (ext == "pdf") args$useDingbats <- FALSE
    do.call(ggplot2::ggsave, args)
  }
}

# Open a base-graphics device for one figure in several formats. Call `draw_fn()` for each.
save_base_nature <- function(draw_fn, path, width_mm = NATURE_SINGLE_MM, height_mm = 70) {
  formats <- trimws(strsplit(Sys.getenv("FIG_FORMATS", "pdf,svg,png,tiff"), ",", fixed = TRUE)[[1]])
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  for (ext in formats) {
    file <- paste0(path, ".", ext)
    if (ext == "pdf") {
      grDevices::cairo_pdf(file, width = w, height = h, pointsize = 6)
    } else if (ext == "svg") {
      grDevices::svg(file, width = w, height = h, pointsize = 6)
    } else if (ext == "png") {
      grDevices::png(file, width = w, height = h, units = "in", res = 600, pointsize = 6, type = "cairo")
    } else if (ext == "tiff") {
      grDevices::tiff(file, width = w, height = h, units = "in", res = 600, pointsize = 6, compression = "lzw",
                      type = "cairo")
    } else {
      next
    }
    tryCatch(draw_fn(), finally = grDevices::dev.off())
  }
}
