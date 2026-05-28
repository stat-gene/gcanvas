# Miami plot: a pair of mirrored Manhattan panels for comparing two
# association results.

#' Miami plot (two stacked, mirrored Manhattan panels)
#'
#' Combines two `manhattan()` objects (or two summary-statistic tables) into
#' a back-to-back layout with shared x-axis and configurable y-axis ratios.
#'
#' @param p1,p2 The two manhattan plots (or summary-stat tables).
#' @param gap Absolute gap (in y units) between the panels.
#' @param gap.ratio Gap expressed as a fraction of the combined panel height.
#' @param y.ratio Length-2 numeric controlling the relative heights of the
#'   top and bottom panels.
#' @param panel.space,panel.space.top,panel.space.bottom,panel.space.left,panel.box
#'   Panel-spacing and panel-border controls.
#' @param x.text.size,y.text.size,tick.size Axis text / tick sizing.
#' @param direction.reverse Logical. Mirror the direction overlay between panels.
#' @param show.legend Logical. Draw the legend.
#' @param drop.internal Logical. Drop internal book-keeping columns from the
#'   returned plot data.
#' @param silent Logical. Suppress progress notes.
#'
#' @return A `ggplot` (or patchwork) object.
#' @export
miami <- function(p1,
                   p2,
                   gap = NULL,
                   gap.ratio = 0.1,
                   y.ratio = c(1, 1),
                   panel.space = NULL,
                   panel.space.top = NULL,
                   panel.space.bottom = NULL,
                   panel.space.left = NULL,
                   panel.box = TRUE,
                   x.text.size = NULL,
                   y.text.size = NULL,
                   tick.size = NULL,
                   direction.reverse = FALSE,
                   show.legend = FALSE,
                   drop.internal = TRUE,
                   silent = FALSE) {
  require_pkg(c("ggplot2", "data.table"))
  silent <- isTRUE(silent)
  panel.box <- isTRUE(panel.box)
  direction.reverse <- isTRUE(direction.reverse)
  y.ratio <- suppressWarnings(as.numeric(y.ratio))
  if (!length(y.ratio)) y.ratio <- c(1, 1)
  if (length(y.ratio) == 1L) y.ratio <- rep(y.ratio, 2L)
  if (length(y.ratio) >= 2L) y.ratio <- y.ratio[1:2]
  bad_ratio <- !is.finite(y.ratio) | is.na(y.ratio) | (y.ratio <= 0)
  if (any(bad_ratio)) y.ratio[bad_ratio] <- 1
  panel.space <- suppressWarnings(as.numeric(panel.space))
  if (!length(panel.space)) panel.space <- c(NA_real_, NA_real_)
  if (length(panel.space) == 1L) panel.space <- rep(panel.space, 2L)
  if (length(panel.space) >= 2L) panel.space <- panel.space[1:2]
  bad_space <- !is.finite(panel.space) | is.na(panel.space)
  if (any(bad_space)) panel.space[bad_space] <- NA_real_
  panel.space[!is.na(panel.space) & panel.space < 0] <- 0
  panel.space.top <- suppressWarnings(as.numeric(panel.space.top))[1]
  if (!is.finite(panel.space.top) || is.na(panel.space.top)) panel.space.top <- NA_real_
  panel.space.bottom <- suppressWarnings(as.numeric(panel.space.bottom))[1]
  if (!is.finite(panel.space.bottom) || is.na(panel.space.bottom)) panel.space.bottom <- NA_real_
  panel.space.left <- suppressWarnings(as.numeric(panel.space.left))[1]
  if (!is.finite(panel.space.left) || is.na(panel.space.left)) panel.space.left <- NA_real_
  x.text.size <- suppressWarnings(as.numeric(x.text.size))[1]
  y.text.size <- suppressWarnings(as.numeric(y.text.size))[1]
  tick.size <- suppressWarnings(as.numeric(tick.size))[1]
  show.legend <- isTRUE(show.legend)
  drop.internal <- isTRUE(drop.internal)
  panel_space_unit_pt <- 10

  if (!inherits(p1, "ggplot") || !inherits(p2, "ggplot")) {
    stop("p1 and p2 must be ggplot objects returned by manhattan().", call. = FALSE)
  }
  meta1 <- attr(p1, "gcanvas_meta")
  meta2 <- attr(p2, "gcanvas_meta")
  if (!is.list(meta1) || !is.list(meta2) ||
      !identical(as.character(meta1$type %||% ""), "manhattan") ||
      !identical(as.character(meta2$type %||% ""), "manhattan")) {
    stop("p1 and p2 must be ggplot objects returned by manhattan().", call. = FALSE)
  }

  .gcanvas_note("gcanvas::miami", .miami_summary_line("p1", p1, meta1), silent = silent)
  .gcanvas_note("gcanvas::miami", .miami_summary_line("p2", p2, meta2), silent = silent)
  .gcanvas_note("gcanvas::miami", "Extracting and merging built Manhattan layers", silent = silent)

  ctx <- .miami_prepare_context(
    p1 = p1, p2 = p2, meta1 = meta1, meta2 = meta2,
    gap = gap, gap.ratio = gap.ratio, y.ratio = y.ratio,
    drop.internal = drop.internal, direction.reverse = direction.reverse
  )
  base <- .miami_compose_plot(
    p1 = p1, p2 = p2, meta1 = meta1, meta2 = meta2, ctx = ctx,
    show.legend = show.legend, panel.box = panel.box,
    x.text.size = x.text.size, y.text.size = y.text.size, tick.size = tick.size
  )
  .miami_finalize_plot(
    base = base, ctx = ctx, meta1 = meta1, meta2 = meta2,
    y.ratio = y.ratio, direction.reverse = direction.reverse,
    panel.space = panel.space, panel.space.top = panel.space.top,
    panel.space.bottom = panel.space.bottom, panel.space.left = panel.space.left,
    panel_space_unit_pt = panel_space_unit_pt, silent = silent
  )
}

