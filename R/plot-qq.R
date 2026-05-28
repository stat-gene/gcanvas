# Q-Q plot of association p-values, with optional category stratification.

#' Q-Q plot of GWAS p-values
#'
#' Plots observed `-log10(p)` against the uniform expectation, optionally
#' stratified by a category column and with optional `lambda_GC` annotation.
#'
#' @param data A `data.frame`/`data.table` of summary statistics.
#' @param p.col Column name for the p-value.
#' @param y.col Optional alternative y-column override.
#' @param category Optional column for category-stratified Q-Q curves.
#' @param legend.position Where the legend is placed.
#' @param title Optional plot title.
#' @param x.title,y.title Axis titles (defaults use expected/observed -log10P).
#' @param ci Logical. Draw the uniform-expectation confidence band.
#' @param lambda.gc Logical. Annotate the genomic inflation factor.
#' @param y.breaks Optional y-axis break vector.
#' @param point.size,point.color,category.color Point styling.
#' @param line.color,line.linewidth Diagonal reference-line styling.
#' @param ci.color,ci.alpha Confidence-band styling.
#' @param panel.box Logical. Draw the panel border.
#' @param grid,grid.major,grid.minor,grid.major.x,grid.major.y,grid.minor.x,grid.minor.y
#'   Grid-line toggles.
#' @param silent Logical. Suppress progress notes.
#'
#' @return A `ggplot` object.
#' @export
qq <- function(data,
               p.col = "P", y.col = NULL,
               category = NULL,
               legend.position = c("bottom", "top", "left", "right", "left.outside", "right.outside"),
               title = NULL,
               x.title = bquote(paste("Expected ", -log[10], " (", italic(P), ")")),
               y.title = bquote(paste("Observed ", -log[10], " (", italic(P), ")")),
               ci = TRUE,
               lambda.gc = FALSE,
               y.breaks = NULL,
               point.size = 1.2,
               point.color = NULL,
               category.color = NULL,
               line.color = "#F94144",
               ci.color = "grey80",
               ci.alpha = 0.3,
               silent = FALSE,
               line.linewidth = 0.5,
               panel.box = TRUE,
               grid = FALSE,
               grid.major = FALSE,
               grid.minor = FALSE,
               grid.major.x = FALSE,
               grid.major.y = TRUE,
               grid.minor.x = FALSE,
               grid.minor.y = FALSE) {
  require_pkg(c("ggplot2", "data.table", "scales"))
  silent <- isTRUE(silent)
  panel.box <- isTRUE(panel.box)
  grid <- isTRUE(grid)
  grid.major <- isTRUE(grid.major) || grid
  grid.minor <- isTRUE(grid.minor) || grid
  grid.major.x <- isTRUE(grid.major.x) || grid.major
  grid.major.y <- isTRUE(grid.major.y) || grid.major
  grid.minor.x <- isTRUE(grid.minor.x) || grid.minor
  grid.minor.y <- isTRUE(grid.minor.y) || grid.minor
  line.linewidth <- suppressWarnings(as.numeric(line.linewidth))[1]
  if (!is.finite(line.linewidth) || is.na(line.linewidth) || line.linewidth <= 0) line.linewidth <- 0.5
  legend.position <- match.arg(tolower(as.character(legend.position)[1]), choices = c("bottom", "top", "left", "right", "left.outside", "right.outside"))
  .qq_numvec <- function(x) suppressWarnings(as.numeric(x))
  .qq_as_char1 <- function(x) {
    if (is.null(x) || !length(x)) return(NULL)
    v <- as.character(x)[1]
    if (is.na(v) || !nzchar(v)) NULL else v
  }
  .qq_as_ref_flag <- function(x) {
    if (is.null(x)) return(logical())
    if (is.logical(x)) return(as.logical(x))
    xc <- tolower(trimws(as.character(x)))
    out <- rep(NA, length(xc))
    out[xc %in% c("true", "t", "1", "yes", "y", "ref", "reference")] <- TRUE
    out[xc %in% c("false", "f", "0", "no", "n", "target")] <- FALSE
    as.logical(out)
  }
  .qq_clean_values <- function(v, use_p = TRUE) {
    if (use_p) {
      v <- .gcanvas_p_filter(v)
      v <- v[!is.na(v) & nzchar(v)]
      return(v)
    }
    v <- .qq_numvec(v)
    v[is.finite(v) & !is.na(v)]
  }
  .qq_normalize_point_colors <- function(point.color, cat_names) {
    cat_names <- as.character(cat_names)
    n_cat <- length(cat_names)
    if (n_cat <= 0L) return(character())
    .qq_is_valid_color <- function(z) {
      if (is.null(z) || length(z) == 0L) return(FALSE)
      zz <- as.character(z)[1]
      if (is.na(zz) || !nzchar(zz)) return(FALSE)
      !is.null(tryCatch(grDevices::col2rgb(zz), error = function(e) NULL))
    }
    .palette_to_colors <- function(mode_name) {
      out <- tryCatch(
        get.colors(mode_name, n_cat, discrete = FALSE, far = TRUE, random = FALSE, plot = FALSE, silent = TRUE),
        error = function(e) NULL
      )
      if (is.null(out) || !length(out)) return(NULL)
      as.character(out)
    }
    if (is.null(point.color) || length(point.color) == 0L) {
      if (n_cat == 1L) return(stats::setNames(rep_len("grey40", n_cat), cat_names))
      return(stats::setNames(
        get.colors("darkrainbow", n_cat, far = FALSE, plot = FALSE, silent = TRUE),
        cat_names
      ))
    }
    pc <- point.color
    if (is.list(pc)) pc <- unlist(pc, use.names = TRUE)
    pc <- as.character(pc)
    pc <- pc[!is.na(pc) & nzchar(pc)]
    if (!length(pc)) {
      return(stats::setNames(
        get.colors("darkrainbow", n_cat, far = FALSE, plot = FALSE, silent = TRUE),
        cat_names
      ))
    }
    if (length(pc) == 1L && !.qq_is_valid_color(pc[1])) {
      pal_cols <- .palette_to_colors(pc[1])
      if (!is.null(pal_cols) && length(pal_cols)) {
        out <- rep_len(pal_cols, n_cat)
        names(out) <- cat_names
        return(out)
      }
    }
    nms <- names(pc)
    if (!is.null(nms) && any(nzchar(nms))) {
      out <- rep(NA_character_, n_cat)
      names(out) <- cat_names
      hit <- match(cat_names, nms)
      ok <- !is.na(hit)
      out[ok] <- pc[hit[ok]]
      miss <- which(!ok)
      if (length(miss)) {
        fill <- pc[is.na(nms) | !nzchar(nms)]
        if (!length(fill)) fill <- pc
        out[miss] <- rep_len(fill, length(miss))
      }
      out[is.na(out) | !nzchar(out)] <- rep_len(
        get.colors("darkrainbow", n_cat, far = FALSE, plot = FALSE, silent = TRUE),
        n_cat
      )[is.na(out) | !nzchar(out)]
      return(out)
    }
    out <- rep_len(pc, n_cat)
    names(out) <- cat_names
    out
  }

  use_p <- is.null(y.col)
  category_col <- .qq_as_char1(category)
  if (!is.null(category.color) && length(category.color) > 0L) {
    point.color <- category.color
  }
  category_levels <- NULL
  q_raw <- data.table::data.table(category = character(), val = numeric())
  p_vals_all <- numeric()

  if (is.list(data) && !is.data.frame(data) && !data.table::is.data.table(data)) {
    if (length(data) == 1L) {
      data <- data[[1]]
    } else {
      nm <- names(data)
      if (is.null(nm)) nm <- character(length(data))
      nm[is.na(nm) | !nzchar(nm)] <- paste0("group", seq_len(length(data)))[is.na(nm) | !nzchar(nm)]
      category_levels <- nm
      for (i in seq_along(data)) {
        vv <- .qq_clean_values(data[[i]], use_p = use_p)
        if (!length(vv)) next
        q_raw <- data.table::rbindlist(list(q_raw, data.table::data.table(category = nm[i], val = vv)), use.names = TRUE)
      }
    }
  }

  if (!nrow(q_raw)) {
    if (is.data.frame(data) || data.table::is.data.table(data)) {
      dt <- if (data.table::is.data.table(data)) data else data.table::as.data.table(data)
      if (use_p) {
        if (!(p.col %in% names(dt))) stop("p.col not found in data.", call. = FALSE)
        vv <- .qq_clean_values(dt[[p.col]], use_p = TRUE)
      } else {
        if (!(y.col %in% names(dt))) stop("y.col not found in data.", call. = FALSE)
        vv <- .qq_clean_values(dt[[y.col]], use_p = FALSE)
      }
      if (!length(vv)) stop("No valid values.", call. = FALSE)
      if (!is.null(category_col) && (category_col %in% names(dt))) {
        if (is.factor(dt[[category_col]])) {
          category_levels <- levels(dt[[category_col]])
        }
        catv <- as.character(dt[[category_col]])
        keep <- if (use_p) {
          x <- .gcanvas_p_filter(dt[[p.col]])
          !is.na(x) & nzchar(x)
        } else {
          x <- .qq_numvec(dt[[y.col]])
          is.finite(x) & !is.na(x)
        }
        catv <- catv[keep]
        valv <- vv
        catv[is.na(catv) | !nzchar(catv)] <- "NA"
        q_raw <- data.table::data.table(category = catv, val = valv)
      } else {
        q_raw <- data.table::data.table(category = "all", val = vv)
      }
    } else {
      vv <- .qq_clean_values(data, use_p = use_p)
      if (!length(vv)) stop("No valid values.", call. = FALSE)
      q_raw <- data.table::data.table(category = "all", val = vv)
    }
  }

  if (!nrow(q_raw)) stop("No valid values.", call. = FALSE)
  q_raw[, category := as.character(category)]
  q_raw[is.na(category) | !nzchar(category), category := "all"]
  # Category order policy: sort(unique(category)) first.
  cat_names <- sort(unique(as.character(q_raw$category)))
  cat_names <- cat_names[!is.na(cat_names) & nzchar(cat_names)]
  if (!length(cat_names)) cat_names <- "all"
  if (use_p) {
    min_idx <- suppressWarnings(which.max(-log10c(q_raw$val)))[1]
    min_p <- if (is.finite(min_idx) && !is.na(min_idx) && min_idx >= 1L) .gcanvas_format_minp(q_raw$val[min_idx], digits = 3L, cutoff = 1e-3) else "NA"
    .gcanvas_note("gcanvas::qq",
                  sprintf("n_variants=%d | n_categories=%d | minP=%s", as_int(nrow(q_raw)), as_int(length(cat_names)), min_p),
                  silent = silent)
  } else {
    .gcanvas_note("gcanvas::qq",
                  sprintf("n_variants=%d | n_categories=%d | maxY=%.4g", as_int(nrow(q_raw)), as_int(length(cat_names)), max(as.numeric(q_raw$val), na.rm = TRUE)),
                  silent = silent)
  }
  if (use_p) p_vals_all <- q_raw$val

  qdt <- q_raw[, {
    o <- sort(if (use_p) -log10c(val) else as.numeric(val), decreasing = TRUE)
    n0 <- length(o)
    e <- -log10(seq_len(n0) / (n0 + 1))
    .(e = e, o = o, n_cat = n0)
  }, by = category]

  qdt[, category := factor(category, levels = cat_names)]
  multi_cat <- length(cat_names) > 1L
  point_cols <- .qq_normalize_point_colors(point.color, cat_names)
  point_cols <- point_cols[cat_names]
  if (isTRUE(multi_cat)) {
    qdt[, cat_draw_order := as_int(category)]
    data.table::setorder(qdt, -cat_draw_order)
    qdt[, cat_draw_order := NULL]
  }

  ci_dt <- NULL
  if (isTRUE(ci)) {
    n_ci <- max(qdt$n_cat, na.rm = TRUE)
    idx <- seq_len(n_ci)
    ci_dt <- data.table::data.table(
      e = -log10(idx / (n_ci + 1)),
      c025 = -log10(stats::qbeta(0.025, idx, n_ci - idx + 1)),
      c975 = -log10(stats::qbeta(0.975, idx, n_ci - idx + 1))
    )
  }

  if (is.null(y.breaks) || (is.character(y.breaks) && length(y.breaks) == 1L &&
                            tolower(trimws(as.character(y.breaks)[1])) == "auto")) {
    y_breaks_draw <- scales::pretty_breaks()(range(c(0, qdt$o), na.rm = TRUE))
  } else {
    y_breaks_draw <- suppressWarnings(as.numeric(y.breaks))
  }
  y_breaks_draw <- unique(y_breaks_draw[is.finite(y_breaks_draw) & y_breaks_draw >= 0])
  if (!length(y_breaks_draw)) y_breaks_draw <- 0

  y_max <- max(qdt$o, if (!is.null(ci_dt)) ci_dt$c975 else qdt$o, na.rm = TRUE) * 1.05
  x_max <- max(qdt$e, na.rm = TRUE) * 1.05
  if (!is.finite(y_max) || y_max <= 0) y_max <- 1
  if (!is.finite(x_max) || x_max <= 0) x_max <- 1
  x_axis_title <- if (is.character(x.title)) as.character(x.title)[1] else x.title
  y_axis_title <- if (is.character(y.title)) as.character(y.title)[1] else y.title
  axis_line_col <- "black"
  axis_line_lwd <- 0.35
  panel_box_lwd <- 0.35
  panel_border_elem <- if (isTRUE(panel.box)) ggplot2::element_rect(fill = NA, colour = "grey20", linewidth = panel_box_lwd) else ggplot2::element_blank()
  grid_major_x_elem <- if (isTRUE(grid.major.x)) ggplot2::element_line(color = "grey90", linewidth = 0.35) else ggplot2::element_blank()
  grid_major_y_elem <- if (isTRUE(grid.major.y)) ggplot2::element_line(color = "grey90", linewidth = 0.35) else ggplot2::element_blank()
  grid_minor_x_elem <- if (isTRUE(grid.minor.x)) ggplot2::element_line(color = "grey90", linewidth = 0.25) else ggplot2::element_blank()
  grid_minor_y_elem <- if (isTRUE(grid.minor.y)) ggplot2::element_line(color = "grey90", linewidth = 0.25) else ggplot2::element_blank()

  p <- ggplot2::ggplot(data = qdt) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.line.x = ggplot2::element_line(colour = axis_line_col, linewidth = axis_line_lwd),
      axis.line.y = ggplot2::element_line(colour = axis_line_col, linewidth = axis_line_lwd),
      panel.border = panel_border_elem,
      panel.grid.major.x = grid_major_x_elem,
      panel.grid.major.y = grid_major_y_elem,
      panel.grid.minor.x = grid_minor_x_elem,
      panel.grid.minor.y = grid_minor_y_elem,
      axis.title = ggplot2::element_text(size = 16),
      axis.text = ggplot2::element_text(size = 14),
      legend.title = ggplot2::element_text(size = 14),
      legend.text = ggplot2::element_text(size = 12)
    ) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linewidth = line.linewidth, color = line.color) +
    ggplot2::coord_cartesian(xlim = c(0, x_max), ylim = c(0, y_max)) +
    ggplot2::scale_x_continuous(name = x_axis_title, expand = c(0, 0)) +
    ggplot2::scale_y_continuous(name = y_axis_title, breaks = y_breaks_draw, expand = c(0, 0)) +
    ggplot2::ggtitle(title)

  if (!is.null(ci_dt)) {
    p <- p + ggplot2::geom_ribbon(
      data = ci_dt,
      ggplot2::aes(x = e, ymin = c025, ymax = c975),
      inherit.aes = FALSE,
      alpha = ci.alpha, fill = ci.color
    )
  }

  if (isTRUE(multi_cat)) {
    p <- p +
      ggplot2::geom_point(ggplot2::aes(x = e, y = o, color = category), size = point.size) +
      ggplot2::scale_color_manual(
        values = point_cols, breaks = cat_names, limits = cat_names, drop = FALSE, name = "Category",
        guide = ggplot2::guide_legend(override.aes = list(size = 3, alpha = 1))
      )
  } else {
    p <- p + ggplot2::geom_point(
      ggplot2::aes(x = e, y = o),
      color = unname(point_cols[1]),
      size = point.size
    )
  }

  lambda_gc_val <- NA_real_
  if (isTRUE(lambda.gc)) {
    if (use_p && length(p_vals_all)) lambda_gc_val <- .gcanvas_lambda_gc(p_vals_all)
  } else if (is.numeric(lambda.gc) && length(lambda.gc) > 0L) {
    lambda_gc_val <- suppressWarnings(as.numeric(lambda.gc))[1]
  }
  lambda_x <- if (identical(legend.position, "left")) x_max * 0.02 else x_max * 0.05
  if (is.finite(lambda_gc_val)) {
    p <- p + ggplot2::annotate(
      "text",
      x = lambda_x, y = y_max * 0.95,
      label = as.expression(bquote(lambda[GC] == .(formatC(lambda_gc_val, format = "f", digits = 3)))),
      size = 6, hjust = 0
    )
  }
  legend_left_y <- if (is.finite(lambda_gc_val)) 0.88 else 0.98

  if (legend.position == "bottom") {
    p <- p + ggplot2::theme(
      legend.position = "bottom",
      legend.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0),
      legend.box.margin = ggplot2::margin(t = -4, r = 0, b = 0, l = 0),
      legend.spacing.y = grid::unit(0.05, "cm")
    )
  } else if (legend.position == "top") {
    p <- p + ggplot2::theme(
      legend.position = "top",
      legend.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0),
      legend.box.margin = ggplot2::margin(t = 0, r = 0, b = -4, l = 0),
      legend.spacing.y = grid::unit(0.05, "cm")
    )
  } else if (legend.position == "left.outside") {
    p <- p + ggplot2::theme(legend.position = "left")
  } else if (legend.position == "right.outside") {
    p <- p + ggplot2::theme(legend.position = "right")
  } else if (legend.position == "left") {
    p <- p + ggplot2::theme(
      legend.position = c(0.02, legend_left_y),
      legend.justification = c(0, 1),
      legend.background = ggplot2::element_rect(fill = scales::alpha("white", 0.4), color = "grey20")
    )
  } else if (legend.position == "right") {
    p <- p + ggplot2::theme(
      legend.position = c(0.98, 0.98),
      legend.justification = c(1, 1),
      legend.background = ggplot2::element_rect(fill = scales::alpha("white", 0.4), color = "grey20")
    )
  }

  attr(p, "gcanvas_meta") <- list(
    type = "qq",
    n_variants = as_int(nrow(q_raw)),
    n_categories = as_int(length(cat_names)),
    categories = cat_names,
    columns = list(p.col = p.col, y.col = y.col, category = category_col),
    ci = isTRUE(ci),
    lambda_gc = if (is.finite(lambda_gc_val)) as.numeric(lambda_gc_val) else NA_real_,
    point.color = point_cols
  )
  p
}

