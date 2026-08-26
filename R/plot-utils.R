# Layout helpers shared across plotting functions: extracting legend grobs
# and shrinking long category labels.

#' Extract the legend grob from a ggplot
#'
#' Pulls the legend (or all legends) out of a `ggplot` object as a `gtable`,
#' optionally wrapping it back into a ggplot for placement in a layout.
#'
#' @param x A `ggplot` object.
#' @param which `"auto"`, `"first"`, or `"all"`: which legend(s) to return.
#' @param drop.empty Logical. Drop legends with no entries.
#' @param as.ggplot Logical. If `TRUE`, return as `ggplot`; otherwise the
#'   underlying `gtable`.
#'
#' @return A legend grob, ggplot wrapper, or list of those.
#' @export
get.legend <- function(x, which = c("auto", "first", "all"), drop.empty = TRUE, as.ggplot = TRUE) {
  require_pkg(c("ggplot2", "gtable", "grid"))
  which <- match.arg(which)

  if (!inherits(x, "ggplot")) {
    stop("x must be a ggplot object.", call. = FALSE)
  }

  gt <- ggplot2::ggplotGrob(x)
  idx <- grepl("^guide-box", gt$layout$name)
  if (!any(idx)) return(NULL)

  legend_names <- gt$layout$name[idx]
  legends <- gt$grobs[idx]

  is_empty_legend <- function(g) {
    if (inherits(g, "zeroGrob")) return(TRUE)
    if (inherits(g, "gtable")) {
      if (!length(g$grobs)) return(TRUE)
      return(all(vapply(g$grobs, inherits, logical(1), what = "zeroGrob")))
    }
    FALSE
  }

  if (isTRUE(drop.empty)) {
    keep <- !vapply(legends, is_empty_legend, logical(1))
    legends <- legends[keep]
    legend_names <- legend_names[keep]
  }

  if (!length(legends)) return(NULL)

  .get_legend_to_plot <- function(lg) {
    p_leg <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::annotation_custom(lg, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
      ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE, clip = "off") +
      ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0))
    attr(p_leg, "legend_grob") <- lg
    p_leg
  }

  if (which == "all") {
    out <- stats::setNames(legends, legend_names)
    if (!isTRUE(as.ggplot)) return(out)
    return(stats::setNames(lapply(out, .get_legend_to_plot), names(out)))
  }
  if (which == "first") {
    if (!isTRUE(as.ggplot)) return(legends[[1]])
    return(.get_legend_to_plot(legends[[1]]))
  }

  pref <- c("guide-box-right", "guide-box-left", "guide-box-bottom", "guide-box-top", "guide-box-inside", "guide-box")
  hit <- match(pref, legend_names, nomatch = 0L)
  hit <- hit[hit > 0L]
  lg <- if (length(hit)) legends[[hit[1]]] else legends[[1]]
  if (!isTRUE(as.ggplot)) return(lg)
  .get_legend_to_plot(lg)
}

# Split a ggplot into point layers vs non-point layers.
# The two returned plots keep the same scales/coords/theme so saved images
# with the same device size can be perfectly overlaid.
#' Split a plot into point-only and label-only layers
#'
#' Splits a `ggplot` into two plots that share an identical layout: one holding
#' only the point layers (the heavy scatter data, e.g. a full genome-wide cloud)
#' and one holding the labels, axes, reference lines, and legend as lightweight
#' vector graphics. Render the point plot as a raster and overlay the label plot
#' as vector to keep text crisp while keeping the file small. Non-point
#' decorations on the point plot are made transparent (not removed) so the two
#' plots align exactly when overlaid at the same device size.
#'
#' `gcanvas` performs no rasterisation itself; this function only separates the
#' layers so you can rasterise the point plot with the tool of your choice (for
#' example `ggrastr::rasterise()`, or by saving it to PNG) and then combine it
#' with the vector label plot.
#'
#' This is a plain function (the dot in its name does not denote an S3 method
#' for [base::split()]).
#'
#' @details
#' For memory efficiency the returned plots hold independent *data* but delegate
#' each layer's *structure* (geom/stat/mapping) to `x` by reference (via
#' [ggplot2::ggproto()]), avoiding a full serialize of the environment-heavy
#' layer objects. They render identically to a full deep copy (verified
#' pixel-for-pixel). The only consequence is that saving a returned plot with
#' [saveRDS()] pulls in `x`'s layer structure and is therefore large; save the
#' rendered image (PNG/PDF) instead, which is unaffected.
#'
#' @usage split.plot_label(x, label.point = c("legend.only", "thin", "keep"), keep.n = 96L)
#'
#' @param x A `ggplot` object, for example the return value of [manhattan()].
#' @param label.point How much of the point-layer data to keep in the label
#'   plot: `"legend.only"` (default) retains only a minimal representative subset
#'   needed to reproduce the legend and hides the points themselves; `"thin"`
#'   retains up to `keep.n` representative points; `"keep"` retains all point
#'   data.
#' @param keep.n Maximum number of points to retain when `label.point = "thin"`.
#'
#' @return A list of two `ggplot` objects: `point` (point layers only, with
#'   non-point decorations drawn transparent) and `label` (all layers, with the
#'   point-layer data reduced according to `label.point`, and panel/plot/legend
#'   backgrounds made transparent so it can be overlaid on the rasterised
#'   `point` plot). Save `label` with a transparent device background
#'   (e.g. `ggsave(..., bg = "transparent")`) when compositing.
#' @rawNamespace export("split.plot_label")
split.plot_label <- function(x, label.point = c("legend.only", "thin", "keep"), keep.n = 96L) {
  require_pkg("ggplot2")
  require_pkg("data.table")
  require_pkg("rlang")
  label.point <- match.arg(label.point)
  keep.n <- as_int(keep.n)
  if (!is.finite(keep.n) || is.na(keep.n) || keep.n < 8L) keep.n <- 8L
  if (!inherits(x, "ggplot")) {
    stop("x must be a ggplot object.", call. = FALSE)
  }

  # Build lightweight skeletons (shared theme/scales/coord/labels) WITHOUT the
  # heavy per-layer data: replacing $layers on a shallow copy does not mutate x,
  # so serialising the skeleton copies almost nothing. This replaces the previous
  # six serialize round-trips (three full-plot + three layer copies) of the
  # ~1e6-row point cloud, which was the time/memory blow-up.
  skel <- x
  skel$layers <- list()
  p_point <- unserialize(serialize(skel, NULL))
  p_label <- unserialize(serialize(skel, NULL))

  if (!length(x$layers)) {
    return(list(point = p_point, label = p_label))
  }

  is_point_layer <- vapply(
    x$layers,
    function(lyr) inherits(lyr$geom, "GeomPoint"),
    logical(1)
  )

  # Copy a layer without serialising its (environment-heavy) structure: give it an
  # independent data copy and delegate geom/stat/mapping/params to the source layer
  # via ggproto(). Serialising a ggplot2 layer expands to ~5x its data because the
  # mapping and geom/stat chains capture data-sized environments; delegation copies
  # only the data. Renders identically (verified pixel-for-pixel); outputs are
  # data-independent but share x's layer *structure* by reference (see @details).
  copy_layer <- function(L) {
    d <- L$data
    if (data.table::is.data.table(d)) d <- data.table::copy(d)
    ggplot2::ggproto(NULL, L, data = d)
  }

  # p_label keeps every layer; its point-layer data is thinned below.
  p_label$layers <- lapply(x$layers, copy_layer)

  # $label is meant to be overlaid on the (rasterised) $point plot, so make its
  # panel/plot/legend backgrounds transparent; an opaque panel fill would hide the
  # points underneath when the two are composited.
  p_label <- p_label + ggplot2::theme(
    panel.background      = ggplot2::element_rect(fill = NA, colour = NA),
    plot.background       = ggplot2::element_rect(fill = NA, colour = NA),
    legend.background     = ggplot2::element_rect(fill = NA, colour = NA),
    legend.box.background = ggplot2::element_rect(fill = NA, colour = NA)
  )

  # Keep plot geometry identical while hiding non-point decorations.
  # Use transparent text/lines (not dropping elements) to preserve layout for overlay.
  p_point <- p_point + ggplot2::theme(
    axis.title.x = ggplot2::element_text(color = "transparent"),
    axis.title.y = ggplot2::element_text(color = "transparent"),
    axis.text.x = ggplot2::element_text(color = "transparent"),
    axis.text.y = ggplot2::element_text(color = "transparent"),
    axis.ticks.x = ggplot2::element_line(color = "transparent"),
    axis.ticks.y = ggplot2::element_line(color = "transparent"),
    axis.line.x = ggplot2::element_line(color = "transparent"),
    axis.line.y = ggplot2::element_line(color = "transparent"),
    panel.border = ggplot2::element_rect(color = "transparent", fill = NA),
    panel.grid.major.x = ggplot2::element_line(color = "transparent"),
    panel.grid.major.y = ggplot2::element_line(color = "transparent"),
    panel.grid.minor.x = ggplot2::element_line(color = "transparent"),
    panel.grid.minor.y = ggplot2::element_line(color = "transparent"),
    plot.title = ggplot2::element_text(color = "transparent"),
    plot.subtitle = ggplot2::element_text(color = "transparent"),
    plot.caption = ggplot2::element_text(color = "transparent"),
    plot.tag = ggplot2::element_text(color = "transparent"),
    strip.text.x = ggplot2::element_text(color = "transparent"),
    strip.text.y = ggplot2::element_text(color = "transparent"),
    strip.background = ggplot2::element_rect(color = "transparent", fill = "transparent")
  )

  # Point-only view: remove legend entirely (outside + inside).
  p_point <- p_point + ggplot2::theme(legend.position = "none")

  # Label view: keep guide legend while shrinking heavy point data
  # to a minimal representative subset.
  .split_label_mapped_vars_for_legend <- function(layer_map, plot_map, data_names) {
    aes_names <- c("colour", "color", "fill", "shape", "linetype", "size", "alpha")
    .vars_from_expr <- function(ex) {
      if (is.null(ex)) return(character())
      e0 <- tryCatch(rlang::get_expr(ex), error = function(e) ex)
      vv <- tryCatch(all.vars(e0), error = function(e) character())
      vv <- as.character(vv)
      vv <- vv[!is.na(vv) & nzchar(vv)]
      vv[vv %in% data_names]
    }
    .extract <- function(mp) {
      if (is.null(mp) || !length(mp)) return(character())
      out <- character()
      for (an in aes_names) {
        ex <- mp[[an]]
        if (is.null(ex)) next
        out <- c(out, .vars_from_expr(ex))
        nm <- tryCatch(rlang::as_name(ex), error = function(e) NA_character_)
        if (is.na(nm) || !nzchar(nm) || !(nm %in% data_names)) {
          lb <- tryCatch(rlang::as_label(ex), error = function(e) NA_character_)
          lb <- trimws(sub("^~", "", as.character(lb)))
          nm <- if (!is.na(lb) && nzchar(lb) && (lb %in% data_names)) lb else NA_character_
        }
        if (!is.na(nm) && nzchar(nm)) out <- c(out, nm)
      }
      unique(out)
    }
    unique(c(.extract(layer_map), .extract(plot_map)))
  }
  .split_label_mapped_vars_all <- function(layer_map, plot_map, data_names) {
    .vars_from_map <- function(mp) {
      if (is.null(mp) || !length(mp)) return(character())
      out <- character()
      for (nm in names(mp)) {
        ex <- mp[[nm]]
        if (is.null(ex)) next
        e0 <- tryCatch(rlang::get_expr(ex), error = function(e) ex)
        vv <- tryCatch(all.vars(e0), error = function(e) character())
        vv <- as.character(vv)
        vv <- vv[!is.na(vv) & nzchar(vv) & vv %in% data_names]
        if (length(vv)) out <- c(out, vv)
      }
      unique(out)
    }
    unique(c(.vars_from_map(layer_map), .vars_from_map(plot_map)))
  }
  .split_label_mapped_var1 <- function(layer_map, plot_map, aes_name, data_names) {
    .get1 <- function(mp) {
      if (is.null(mp) || !length(mp) || is.null(mp[[aes_name]])) return(NA_character_)
      ex <- mp[[aes_name]]
      nm <- tryCatch(rlang::as_name(ex), error = function(e) NA_character_)
      if (!is.na(nm) && nzchar(nm) && (nm %in% data_names)) return(nm)
      lb <- tryCatch(rlang::as_label(ex), error = function(e) NA_character_)
      lb <- trimws(sub("^~", "", as.character(lb)))
      if (!is.na(lb) && nzchar(lb) && (lb %in% data_names)) return(lb)
      NA_character_
    }
    v <- .get1(layer_map)
    if (!is.na(v)) return(v)
    .get1(plot_map)
  }
  .split_label_legend_only_point_data <- function(src, layer_map, plot_map, max_levels = 160L) {
    if (!is.data.frame(src) || !nrow(src)) return(src)
    dt <- if (data.table::is.data.table(src)) data.table::copy(src) else data.table::as.data.table(src)
    map_vars <- .split_label_mapped_vars_for_legend(layer_map, plot_map, names(dt))
    map_vars_all <- .split_label_mapped_vars_all(layer_map, plot_map, names(dt))
    if (!length(map_vars)) {
      keep_cols <- unique(c(map_vars_all, names(dt)))
      out <- dt[1, ..keep_cols]
    } else {
      keep_cols <- unique(c(map_vars, map_vars_all))
      out <- unique(dt[, ..keep_cols])
      if (nrow(out) > as_int(max_levels)) {
        idx <- unique(as_int(round(seq(1, nrow(out), length.out = as_int(max_levels)))))
        out <- out[idx]
      }
    }
    x_var <- .split_label_mapped_var1(layer_map, plot_map, "x", names(dt))
    y_var <- .split_label_mapped_var1(layer_map, plot_map, "y", names(dt))
    if (!is.na(x_var) && nzchar(x_var) && !(x_var %in% names(out))) out[, (x_var) := NA_real_]
    if (!is.na(y_var) && nzchar(y_var) && !(y_var %in% names(out))) out[, (y_var) := NA_real_]
    if (!is.na(x_var) && nzchar(x_var) && (x_var %in% names(out))) out[, (x_var) := NA_real_]
    if (!is.na(y_var) && nzchar(y_var) && (y_var %in% names(out))) out[, (y_var) := NA_real_]
    out[]
  }
  .split_label_thin_point_data <- function(df, map_vars, keep_n = 96L) {
    if (!is.data.frame(df) || !nrow(df)) return(df)
    n <- as_int(nrow(df))
    keep_n <- as_int(keep_n)
    if (!is.finite(keep_n) || is.na(keep_n) || keep_n < 8L) keep_n <- 8L
    if (n <= keep_n) return(if (data.table::is.data.table(df)) data.table::copy(df) else df)

    idx <- c(1L, n)
    for (v in unique(map_vars)) {
      if (!(v %in% names(df))) next
      xv <- df[[v]]
      if (is.numeric(xv)) {
        ok <- which(is.finite(xv) & !is.na(xv))
        if (!length(ok)) next
        xok <- xv[ok]
        qs <- unique(stats::quantile(xok, probs = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1), na.rm = TRUE, names = FALSE))
        for (qv in qs) {
          j <- ok[which.min(abs(xok - qv))]
          if (is.finite(j) && !is.na(j)) idx <- c(idx, as_int(j))
        }
      } else {
        ch <- as.character(xv)
        lev <- unique(ch[!is.na(ch) & nzchar(ch)])
        if (!length(lev)) next
        if (length(lev) > 48L) lev <- lev[unique(as_int(round(seq(1, length(lev), length.out = 48L))))]
        for (lv in lev) {
          j <- which(ch == lv)[1]
          if (is.finite(j) && !is.na(j)) idx <- c(idx, as_int(j))
        }
      }
    }
    idx <- unique(as_int(idx))
    idx <- idx[is.finite(idx) & !is.na(idx) & idx >= 1L & idx <= n]
    target <- min(keep_n, n)
    if (length(idx) < target) {
      rem <- setdiff(seq_len(n), idx)
      need <- target - length(idx)
      if (need > 0L && length(rem)) {
        add_idx <- rem[unique(as_int(round(seq(1, length(rem), length.out = need))))]
        idx <- c(idx, add_idx)
      }
    }
    idx <- sort(unique(as_int(idx)))
    if (data.table::is.data.table(df)) return(data.table::copy(df[idx]))
    df[idx, , drop = FALSE]
  }

  if (length(p_label$layers)) {
    for (i in seq_along(p_label$layers)) {
      lyr <- p_label$layers[[i]]
      if (!inherits(lyr$geom, "GeomPoint")) next
      src <- NULL
      # Read-only source; the thinning helpers subset/copy internally and the
      # "keep" branch copies explicitly, so no defensive copy is needed here.
      if (is.data.frame(lyr$data) && nrow(lyr$data)) {
        src <- lyr$data
      } else if (is.data.frame(p_label$data) && nrow(p_label$data)) {
        src <- p_label$data
      }
      if (!is.null(src) && nrow(src)) {
        map_vars <- .split_label_mapped_vars_for_legend(lyr$mapping, p_label$mapping, names(src))
        if (identical(label.point, "keep")) {
          lyr$data <- if (data.table::is.data.table(src)) data.table::copy(src) else src
        } else if (identical(label.point, "thin")) {
          lyr$data <- .split_label_thin_point_data(src, map_vars = map_vars, keep_n = keep.n)
        } else {
          lyr$data <- .split_label_legend_only_point_data(src, lyr$mapping, p_label$mapping, max_levels = max(32L, keep.n))
          lyr$geom_params$na.rm <- TRUE
          lyr$stat_params$na.rm <- TRUE
        }
      } else {
        lyr$data <- data.frame(.gcanvas_dummy = 1)
      }
      p_label$layers[[i]] <- lyr
    }
  }

  # Point-only plot: point layers with independent data, structure delegated.
  p_point$layers <- lapply(x$layers[is_point_layer], copy_layer)

  attr(p_point, "gcanvas_meta") <- attr(x, "gcanvas_meta")
  attr(p_label, "gcanvas_meta") <- attr(x, "gcanvas_meta")
  attr(p_point, "split.plot_label") <- "point"
  attr(p_label, "split.plot_label") <- "label"

  list(point = p_point, label = p_label)
}

