# Genome-wide Manhattan plot of association p-values.

#' Manhattan plot of GWAS p-values
#'
#' Draws a genome-wide Manhattan plot with chromosome banding, configurable
#' significance thresholds, lead-variant annotation, and optional split-track
#' overlays.
#'
#' @param data A `data.frame`/`data.table` of summary statistics.
#' @param snp.col,chrom.col,pos.col,p.col Column names in `data`.
#' @param y.col Optional column to use for the y axis (default `-log10(p)`).
#' @param lead.label.col Optional column whose values label lead variants.
#' @param build Genome build (`37` or `38`).
#' @param chroms,chroms.drop Which chromosomes to keep / drop.
#' @param alpha,point.size,chrom.color,dark Base point styling and theme mode.
#' @param lead,lead.color,lead.size,lead.stroke Lead-variant highlighting.
#' @param lead.flank,lead.flank.snp,lead.flank.color,lead.flank.alpha
#'   Flanking-region highlighting around each lead variant.
#' @param lead.label,lead.label.color,lead.label.size,lead.label.angle
#'   Lead-variant text labels.
#' @param drag.label,drag.label.arrow,drag.label.line.color,drag.label.linewidth
#'   Repulsion / leader-line behavior for drag-style labels.
#' @param seed Integer seed for the label layout RNG.
#' @param y.max,y.breaks Y-axis controls.
#' @param y.rescale.at,y.rescale.ratio,y.rescale.breaks,y.rescale.line,y.rescale.line.color,y.rescale.line.type,y.rescale.line.linewidth
#'   Y-axis break-and-rescale controls for compressing the top of the y-axis.
#' @param y.ceiling,y.ceiling.line,y.ceiling.line.color,y.ceiling.line.type,y.ceiling.line.linewidth
#'   Y-axis ceiling clipping.
#' @param panel.space,panel.space.top,panel.space.bottom,panel.space.left,panel.box
#'   Panel-spacing and panel-border controls.
#' @param line.x Logical. Draw the horizontal x-axis line.
#' @param drag.chrom,drag.chrom.tick Compact / stagger small chromosome labels.
#' @param threshold,threshold.color,threshold.type,threshold.linewidth,line.alpha
#'   Significance line controls.
#' @param direction,beta.col,a1.col,a2.col,direction.a1,direction.color,direction.color.lead,direction.size
#'   Direction-of-effect glyph controls.
#' @param x.text.size,y.text.size,tick.size Axis text / tick sizing.
#' @param grid,grid.major,grid.minor,grid.major.x,grid.major.y,grid.minor.x,grid.minor.y
#'   Grid-line toggles.
#' @param title,x.title,y.title,x.title.nudge Plot/axis titles.
#' @param silent Logical. Suppress progress notes.
#'
#' @return A `ggplot` object.
#' @export
manhattan <- function(data,
                      snp.col = "SNP", chrom.col = "CHR", pos.col = "POS", p.col = "P",
                      y.col = NULL, lead.label.col = NULL,
                      build = 38L,
                      chroms = "auto",
                      chroms.drop = FALSE,
                      alpha = 1,
                      point.size = 1.4,
                      chrom.color = c("grey60", "grey80"),
                      dark = FALSE,
                      lead = NULL,
                      lead.color = "#E63946",
                      lead.size = NULL,
                      lead.stroke = 1,
                      lead.flank = NULL,
                      lead.flank.snp = NULL,
                      lead.flank.color = NULL,
                      lead.flank.alpha = 1,
                      lead.label = FALSE,
                      lead.label.color = "grey20",
                      drag.label = FALSE,
                      drag.label.arrow = TRUE,
                      drag.label.line.color = "grey70",
                      drag.label.linewidth = 0.35,
                      lead.label.size = 4,
                      lead.label.angle = 0,
                      seed = 23L,
                      y.max = NULL,
                      y.breaks = NULL,
                      y.rescale.at = NULL,
                      y.rescale.ratio = 0.25,
                      y.rescale.breaks = NULL,
                      y.rescale.line = TRUE,
                      y.rescale.line.color = "grey80",
                      y.rescale.line.type = "solid",
                      y.rescale.line.linewidth = 0.5,
                      y.ceiling = NULL,
                      y.ceiling.line = TRUE,
                      y.ceiling.line.color = "grey80",
                      y.ceiling.line.type = "solid",
                      y.ceiling.line.linewidth = 0.5,
                      panel.space = NULL,
                      panel.space.top = NULL,
                      panel.space.bottom = NULL,
                      panel.space.left = NULL,
                      panel.box = TRUE,
                      line.x = TRUE,
                      drag.chrom = c(19, 21),
                      drag.chrom.tick = NULL,
                      threshold = 5e-8,
                      threshold.color = "grey20",
                      threshold.type = NULL,
                      threshold.linewidth = 0.7,
                      line.alpha = 0.8,
                      direction = FALSE,
                      beta.col = "BETA",
                      a1.col = NULL,
                      a2.col = NULL,
                      direction.a1 = NULL,
                      direction.color = c("#EF476F", "#258AB2"),
                      direction.color.lead = NULL,
                      direction.size = NULL,
                      x.text.size = 14,
                      y.text.size = 14,
                      tick.size = NULL,
                      grid = FALSE,
                      grid.major = FALSE,
                      grid.minor = FALSE,
                      grid.major.x = FALSE,
                      grid.major.y = TRUE,
                      grid.minor.x = TRUE,
                      grid.minor.y = FALSE,
                      title = NULL,
                      x.title = "Chromosome",
                      x.title.nudge = 0,
                      y.title = NULL,
                      silent = FALSE) {
  require_pkg(c("ggplot2", "data.table", "scales"))
  if (!is.data.frame(data) && !data.table::is.data.table(data)) stop("data must be a data.frame/data.table.", call. = FALSE)
  old_dt_verbose <- getOption("datatable.verbose")
  options(datatable.verbose = FALSE)
  on.exit(options(datatable.verbose = old_dt_verbose), add = TRUE)
  list2env(.manhattan_resolve_options(
    silent = silent,
    panel.box = panel.box,
    x.text.size = x.text.size,
    y.text.size = y.text.size,
    x.title.nudge = x.title.nudge,
    tick.size = tick.size,
    panel.space = panel.space,
    panel.space.top = panel.space.top,
    panel.space.bottom = panel.space.bottom,
    panel.space.left = panel.space.left,
    grid = grid,
    grid.major = grid.major,
    grid.minor = grid.minor,
    grid.major.x = grid.major.x,
    grid.major.y = grid.major.y,
    grid.minor.x = grid.minor.x,
    grid.minor.y = grid.minor.y,
    direction = direction,
    direction.color = direction.color,
    direction.color.lead = direction.color.lead,
    line.alpha = line.alpha,
    alpha = alpha,
    lead.flank.alpha = lead.flank.alpha,
    point.size = point.size,
    direction.size = direction.size,
    drag.chrom = drag.chrom,
    drag.chrom.tick = drag.chrom.tick,
    drag.label = drag.label,
    drag.label.arrow = drag.label.arrow,
    chroms = chroms,
    chroms.drop = chroms.drop,
    drag.label.linewidth = drag.label.linewidth,
    lead.stroke = lead.stroke,
    lead.label.size = lead.label.size,
    lead.label.angle = lead.label.angle,
    seed = seed,
    line.x = line.x,
    y.ceiling.line = y.ceiling.line,
    y.rescale.line = y.rescale.line,
    y.ceiling.line.color = y.ceiling.line.color,
    y.ceiling.line.type = y.ceiling.line.type,
    y.ceiling.line.linewidth = y.ceiling.line.linewidth,
    y.rescale.line.color = y.rescale.line.color,
    y.rescale.line.type = y.rescale.line.type,
    y.rescale.line.linewidth = y.rescale.line.linewidth,
    lead.color = lead.color
  ), environment())

  list2env(.manhattan_prepare_input(
    data = data,
    snp.col = snp.col,
    chrom.col = chrom.col,
    pos.col = pos.col,
    p.col = p.col,
    y.col = y.col,
    lead.label.col = lead.label.col,
    direction = direction,
    beta.col = beta.col,
    a1.col = a1.col,
    a2.col = a2.col,
    direction.a1 = direction.a1
  ), environment())
  .gcanvas_note("gcanvas::manhattan", summary_msg, silent = silent)

  list2env(.manhattan_prepare_layout(
    dt = dt,
    build = build,
    chroms_auto = chroms_auto,
    chroms.drop = chroms.drop,
    chroms_req = chroms_req,
    threshold = threshold,
    y.col = y.col,
    threshold.color = threshold.color,
    threshold.type = threshold.type,
    threshold.linewidth = threshold.linewidth,
    y.rescale.at = y.rescale.at,
    y.rescale.ratio = y.rescale.ratio,
    y.rescale.breaks = y.rescale.breaks,
    y.ceiling = y.ceiling,
    dark = dark,
    chrom.color = chrom.color,
    alpha = alpha,
    silent = silent
  ), environment())

  lead_for_tbl <- lead
  pass_mode <- FALSE
  if (is.atomic(lead) || (is.list(lead) && !is.data.frame(lead) && !data.table::is.data.table(lead))) {
    lead_raw <- .gcanvas_as_snp_vector(lead)
    lead_tok <- tolower(lead_raw)
    pass_mode <- any(lead_tok %in% c("pass", "threshold"))
    if (isTRUE(pass_mode)) {
      lead_keep <- lead_raw[!(lead_tok %in% c("pass", "threshold"))]
      lead_for_tbl <- if (length(lead_keep)) lead_keep else FALSE
    }
  }
  if (isTRUE(direction) && isTRUE(pass_mode) && !("beta" %in% names(dt))) {
    stop("direction=TRUE with lead='pass'/'threshold' requires beta.col in data.", call. = FALSE)
  }

  lead_tbl <- .manhattan_normalize_lead_table(
    lead = lead_for_tbl,
    lead.color = lead.color,
    lead.size = lead.size,
    lead.stroke = lead.stroke,
    lead.flank.color = lead.flank.color,
    lead.label = lead.label,
    lead.label.color = lead.label.color
  )
  flank_lead_input <- if (is.null(lead.flank.snp) || length(lead.flank.snp) == 0L ||
                          (is.logical(lead.flank.snp) && length(lead.flank.snp) == 1L && !isTRUE(lead.flank.snp))) {
    lead_for_tbl
  } else {
    lead.flank.snp
  }
  flank_tbl <- .manhattan_normalize_lead_table(
    lead = flank_lead_input,
    lead.color = lead.color,
    lead.size = lead.size,
    lead.stroke = lead.stroke,
    lead.flank.color = lead.flank.color,
    lead.label = FALSE,
    lead.label.color = lead.label.color
  )
  flank_tbl <- .manhattan_resolve_lead_positions(flank_tbl, dt)

  lead_derived <- .manhattan_derive_leads(
    lead_tbl = lead_tbl,
    dt = dt,
    data = data,
    snp_col_use = snp_col_use,
    lead_label_col_std = lead_label_col_std,
    lead.label = lead.label,
    pass_mode = pass_mode,
    threshold_zero_mode = threshold_zero_mode,
    threshold_lines = threshold_lines,
    lead_color_enabled = lead_color_enabled,
    lead_col0 = lead_col0
  )
  lead_tbl <- lead_derived$lead_tbl
  dt <- lead_derived$dt
  lead_draw_map <- lead_derived$lead_draw_map
  .gcanvas_note("gcanvas::manhattan", "Lead SNP and flank annotation", silent = silent)

  flank_bp <- .gcanvas_parse_bp_span(lead.flank, arg_name = "lead.flank")
  lead_flank_color_supplied <- {
    x <- lead.flank.color
    if (is.list(x) && !is.data.frame(x) && !data.table::is.data.table(x)) x <- unlist(x, use.names = TRUE)
    x <- as.character(x)
    any(!is.na(x) & nzchar(x))
  }
  suppress_explicit_lead <- (!isTRUE(lead_color_enabled)) &&
    is.finite(flank_bp) && !is.na(flank_bp) && flank_bp > 0 &&
    isTRUE(lead_flank_color_supplied)
  lead_hit_rows <- if (isTRUE(suppress_explicit_lead)) .manhattan_lead_row_idx(dt, lead_tbl) else integer()

  if (nrow(flank_tbl) && is.finite(flank_bp) && !is.na(flank_bp) && flank_bp > 0) {
    win <- flank_tbl[is.finite(POS) & !is.na(CHR) & nzchar(CHR),
                    .(lead_id, CHR, start = pmax(1, POS - flank_bp), end = POS + flank_bp, lead_pos = POS, flank_color)]
    if (nrow(win)) {
      q_dt <- dt[, .(row_id = .I, CHR, POS)]
      hits <- win[q_dt, on = .(CHR, start <= POS, end >= POS), nomatch = 0L, allow.cartesian = TRUE]
      if (nrow(hits)) {
        hits[, dist := abs(start - lead_pos)]
        # When a point falls in several lead flanks, the LAST lead in the table
        # (largest lead_id) wins, matching the lead-point z-order where later leads
        # are drawn on top; distance only breaks ties within the same lead. (Was
        # nearest-lead-wins, which conflicted with the point ordering.)
        data.table::setorder(hits, row_id, -lead_id, dist)
        hits <- hits[!duplicated(row_id)]
        if (length(lead_hit_rows)) hits <- hits[!(row_id %in% lead_hit_rows)]
        dt[hits$row_id, `:=`(plot_color = hits$flank_color, plot_alpha = lead.flank.alpha, flank_highlight = TRUE)]
      }
    }
  }

  dt <- .manhattan_mark_ceiling_candidates(dt = dt, lead_tbl = lead_tbl, use_ceiling = use_ceiling)

  .gcanvas_note("gcanvas::manhattan", "Building plot layers", silent = silent)

  y_top <- suppressWarnings(as.numeric(y.max))[1]
  y_max_auto <- is.null(y.max) || length(y.max) == 0L ||
    is.na(as.character(y.max)[1]) || !nzchar(as.character(y.max)[1]) ||
    identical(tolower(as.character(y.max)[1]), "auto")
  if (isTRUE(use_ceiling) && is.finite(y_ceiling0)) {
    y_top <- .manhattan_y_map(y_ceiling0)
  } else {
    if (!is.finite(y_top) || isTRUE(y_max_auto)) {
      y_top <- max(c(dt$ydraw, line_dt$y), na.rm = TRUE)
      if (!is.finite(y_top) || y_top <= 0) y_top <- 1
      y_top <- y_top * 1.05
    }
    has_lead_label <- isTRUE(nrow(lead_tbl) > 0L) && any(!is.na(lead_tbl$lead_label) & nzchar(lead_tbl$lead_label))
    if (isTRUE(has_lead_label)) {
      y_top <- y_top * if (isTRUE(is_angle_zero)) 1.10 else 1.24
    } else if (!isTRUE(is_angle_zero) && isTRUE(nrow(lead_tbl) > 0L)) {
      y_top <- y_top * 1.18
    }
    if (!is.finite(y_top) || y_top <= 0) y_top <- 1
  }
  y_span <- max(1, y_top)
  tick_len_auto <- max(0.06, y_span * 0.06)
  .manhattan_ceiling_nudge <- function(size_val) {
    sz <- suppressWarnings(as.numeric(size_val))
    base_sz <- suppressWarnings(as.numeric(point.size))[1]
    if (!is.finite(base_sz) || is.na(base_sz) || base_sz <= 0) base_sz <- 1.4
    if (!length(sz)) return(numeric())
    sz[!is.finite(sz) | is.na(sz) | sz <= 0] <- base_sz
    base_nudge <- max(0.02, y_span * 0.006)
    as.numeric(base_nudge * (sz / base_sz))
  }
  label_dist_mult <- 2
  label_nudge_base <- max(0.04, min(0.28, y_span * 0.01)) * label_dist_mult
  label_nudge_step <- max(0.01, min(0.10, y_span * 0.0035)) * label_dist_mult
  dt[, y_plot := ydraw]
  lead_size0 <- suppressWarnings(as.numeric(lead.size))
  if (!length(lead_size0) || all(!is.finite(lead_size0))) {
    lead_size0 <- c(point.size * 1.6, point.size * 2.2)
  } else if (length(lead_size0) == 1L) {
    lead_size0 <- c(lead_size0[1], lead_size0[1] * 1.35)
  } else {
    lead_size0 <- lead_size0[1:2]
  }
  lead_size0[!is.finite(lead_size0) | lead_size0 <= 0] <- point.size * c(1.6, 2.2)

  dir_derived <- .manhattan_derive_direction_leads(
    dt = dt,
    lead_draw_map = lead_draw_map,
    direction = direction,
    y.col = y.col,
    threshold = threshold,
    threshold_zero_mode = threshold_zero_mode,
    threshold_lines = threshold_lines,
    direction_color_inherit = direction_color_inherit,
    dir_col_pos = dir_col_pos,
    dir_col_neg = dir_col_neg,
    direction_size0 = direction_size0,
    lead_size0 = lead_size0,
    lead_color_enabled = lead_color_enabled,
    lead_col0 = lead_col0,
    direction_lead_fill_override = direction_lead_fill_override,
    dir_col_lead_pos = dir_col_lead_pos,
    dir_col_lead_neg = dir_col_lead_neg,
    use_ceiling = use_ceiling,
    ceiling_nudge_fun = .manhattan_ceiling_nudge
  )
  dt <- dir_derived$dt
  dir_dt <- dir_derived$dir_dt
  dir_row_idx <- dir_derived$dir_row_idx

  # Prevent double-drawing of lead points: draw them only in explicit lead layer.
  lead_draw_idx <- integer()
  if (nrow(lead_tbl)) {
    lead_draw_snp <- unique(lead_tbl$snp[!is.na(lead_tbl$snp) & nzchar(lead_tbl$snp)])
    if (length(lead_draw_snp)) {
      idx_snp <- dt[snp %in% lead_draw_snp, which = TRUE]
      lead_draw_idx <- c(lead_draw_idx, idx_snp[!is.na(idx_snp)])
    }
    lead_draw_pos <- unique(lead_tbl[is.finite(POS_data) & !is.na(CHR_data) & nzchar(CHR_data), .(CHR = CHR_data, POS = POS_data)])
    if (nrow(lead_draw_pos)) {
      idx_pos <- dt[lead_draw_pos, on = .(CHR, POS), which = TRUE]
      lead_draw_idx <- c(lead_draw_idx, idx_pos[!is.na(idx_pos)])
    }
  }
  lead_draw_idx <- unique(as_int(lead_draw_idx))
  dt[, lead_draw := FALSE]
  if (length(lead_draw_idx)) dt[lead_draw_idx, lead_draw := TRUE]

  if (is.null(y.breaks) || (is.character(y.breaks) && length(y.breaks) == 1L &&
                            tolower(trimws(as.character(y.breaks)[1])) == "auto")) {
    if (isTRUE(use_y_rescale)) {
      raw_max <- max(c(dt$ydraw_raw, line_dt$y_raw), na.rm = TRUE)
      if (!is.finite(raw_max) || raw_max <= 0) raw_max <- max(dt$ydraw_raw, na.rm = TRUE)
      if (!is.finite(raw_max) || raw_max <= 0) raw_max <- y_rescale_at
      br_low <- suppressWarnings(as.numeric(base::pretty(c(0, y_rescale_at), n = 3)))
      br_low <- br_low[is.finite(br_low) & br_low <= y_rescale_at]
      if (length(y_rescale_breaks)) {
        br_high <- sort(unique(y_rescale_breaks[y_rescale_breaks >= y_rescale_at]))
      } else {
        br_high <- suppressWarnings(as.numeric(base::pretty(c(y_rescale_at, raw_max), n = 2)))
        br_high <- br_high[is.finite(br_high) & br_high >= y_rescale_at]
      }
      y_breaks_raw <- sort(unique(c(0, br_low, y_rescale_at, br_high)))
      # Keep rescaled axis readable: aggressively thin to a small set of representative ticks.
      if (length(y_breaks_raw) > 6L) {
        lo <- y_breaks_raw[y_breaks_raw < y_rescale_at]
        hi <- y_breaks_raw[y_breaks_raw > y_rescale_at]
        lo_keep <- if (length(lo) > 2L) lo[unique(as_int(round(seq(1, length(lo), length.out = 2))))] else lo
        hi_keep <- if (length(hi) > 2L) hi[unique(as_int(round(seq(1, length(hi), length.out = 2))))] else hi
        y_breaks_raw <- sort(unique(c(0, lo_keep, y_rescale_at, hi_keep)))
      }
      if (!length(y_breaks_raw)) y_breaks_raw <- c(0, y_rescale_at)
      y_breaks_draw <- .manhattan_y_map(y_breaks_raw)
      y_breaks_labels <- as.character(signif(y_breaks_raw, 4))
    } else {
      y_breaks_draw <- scales::pretty_breaks()(c(0, y_top))
      y_breaks_labels <- ggplot2::waiver()
    }
  } else {
    y_breaks_raw <- suppressWarnings(as.numeric(y.breaks))
    if (isTRUE(use_y_rescale)) {
      y_breaks_draw <- .manhattan_y_map(y_breaks_raw)
      y_breaks_labels <- as.character(signif(y_breaks_raw, 4))
    } else {
      y_breaks_draw <- y_breaks_raw
      y_breaks_labels <- ggplot2::waiver()
    }
  }
  if (!inherits(y_breaks_labels, "waiver")) {
    br_dt <- data.table::data.table(draw = y_breaks_draw, lab = y_breaks_labels)
    br_dt <- br_dt[is.finite(draw) & draw >= 0]
    if (nrow(br_dt)) {
      data.table::setorder(br_dt, draw)
      br_dt[, draw_key := round(draw, 6)]
      br_dt <- br_dt[!duplicated(draw_key)]
      br_dt[, draw_key := NULL]
      br_dt <- br_dt[!duplicated(draw)]
      if (isTRUE(use_y_rescale) && is.finite(y_rescale_at) && !is.na(y_rescale_at) && nrow(br_dt) >= 2L) {
        y_at_draw <- .manhattan_y_map(y_rescale_at)
        min_gap_at <- max(1.0, y_span * 0.04)
        idx_at <- which.min(abs(br_dt$draw - y_at_draw))
        if (length(idx_at) == 1L && is.finite(br_dt$draw[idx_at])) {
          if (idx_at > 1L) {
            gap_below <- y_at_draw - br_dt$draw[idx_at - 1L]
            if (is.finite(gap_below) && gap_below < min_gap_at) {
              br_dt <- br_dt[-(idx_at - 1L)]
            }
          }
        }
      }
      if (!any(abs(br_dt$draw) < .Machine$double.eps)) {
        br_dt <- data.table::rbindlist(list(data.table::data.table(draw = 0, lab = "0"), br_dt), use.names = TRUE)
        br_dt <- br_dt[!duplicated(draw)]
      }
      y_breaks_draw <- br_dt$draw
      y_breaks_labels <- br_dt$lab
    } else {
      y_breaks_draw <- 0
      y_breaks_labels <- "0"
    }
  } else {
    y_breaks_draw <- unique(y_breaks_draw[is.finite(y_breaks_draw) & y_breaks_draw >= 0])
    if (!length(y_breaks_draw)) y_breaks_draw <- 0
    if (!isTRUE(any(abs(y_breaks_draw) < .Machine$double.eps))) y_breaks_draw <- c(0, y_breaks_draw)
    y_breaks_draw <- sort(unique(y_breaks_draw))
  }
  if (!exists("y_breaks_labels", inherits = FALSE)) y_breaks_labels <- ggplot2::waiver()

  x_axis_title <- if (is.null(x.title)) {
    "Chromosome"
  } else if (is.character(x.title)) {
    as.character(x.title)[1]
  } else {
    x.title
  }
  y_axis_title <- if (!is.null(y.title)) {
    if (is.character(y.title)) as.character(y.title)[1] else y.title
  } else if (!is.null(y.col)) {
    y.col
  } else {
    bquote(paste(-log[10], " (", italic(P), ")"))
  }
  shape_map <- if (use_ceiling) c("0" = 16, "1" = 17) else c("0" = 16)
  size_map <- if (use_ceiling) c("0" = point.size, "1" = point.size * 2.2) else c("0" = point.size)
  base_shape <- suppressWarnings(as.numeric(shape_map["0"]))[1]
  if (!is.finite(base_shape) || is.na(base_shape)) base_shape <- 16
  ceil_shape <- suppressWarnings(as.numeric(shape_map["1"]))[1]
  if (!is.finite(ceil_shape) || is.na(ceil_shape)) ceil_shape <- base_shape
  base_size <- suppressWarnings(as.numeric(size_map["0"]))[1]
  if (!is.finite(base_size) || is.na(base_size) || base_size <= 0) base_size <- point.size
  ceil_size <- suppressWarnings(as.numeric(size_map["1"]))[1]
  if (!is.finite(ceil_size) || is.na(ceil_size) || ceil_size <= 0) ceil_size <- base_size
  special_chr_ticks <- intersect(drag_chrom_vec, chr_map$CHR)
  x_breaks <- chr_map$tick
  x_labels <- chr_map$CHR
  if (length(special_chr_ticks)) {
    x_labels[chr_map$CHR %in% special_chr_ticks] <- paste0("\n", x_labels[chr_map$CHR %in% special_chr_ticks])
  }

  x_max <- max(chr_map$chr_offset + chr_map$chr_width, na.rm = TRUE)
  if (!is.finite(x_max) || x_max <= 0) x_max <- max(dt$x, na.rm = TRUE)
  if (!is.finite(x_max) || x_max <= 0) x_max <- 1
  x_pad <- x_max * 0.01
  panel_space_unit_pt <- 10
  .manhattan_panel_space_to_pt <- function(x) {
    x <- suppressWarnings(as.numeric(x))[1]
    if (!is.finite(x) || is.na(x) || x < 0) return(NA_real_)
    as.numeric(x * panel_space_unit_pt)
  }
  bottom_margin_pt <- 8
  top_margin_pt <- if (isTRUE(drag.label)) 88 else 12
  if (!is.na(panel.space[1])) top_margin_pt <- .manhattan_panel_space_to_pt(panel.space[1])
  if (!is.na(panel.space[2])) bottom_margin_pt <- .manhattan_panel_space_to_pt(panel.space[2])
  if (!is.na(panel.space.top)) top_margin_pt <- .manhattan_panel_space_to_pt(panel.space.top)
  if (!is.na(panel.space.bottom)) bottom_margin_pt <- .manhattan_panel_space_to_pt(panel.space.bottom)
  top_margin_user_set <- (!is.na(panel.space[1]) || !is.na(panel.space.top))
  label_ang0 <- suppressWarnings(as.numeric(lead.label.angle))[1]
  if (!is.finite(label_ang0) || is.na(label_ang0)) label_ang0 <- 0
  ang_fac0 <- abs(sin(label_ang0 * pi / 180))
  ceiling_label_n0 <- if (isTRUE(use_ceiling) && isTRUE(nrow(lead_tbl) > 0L)) {
    sum(as_int(lead_tbl$ceiling_data) == 1L & !is.na(lead_tbl$lead_label) & nzchar(lead_tbl$lead_label))
  } else 0L
  has_ceiling_label0 <- ceiling_label_n0 > 0L
  if (isTRUE(has_ceiling_label0) && !isTRUE(drag.label) && !isTRUE(top_margin_user_set)) {
    extra_n <- min(6, as_int(ceiling_label_n0))
    top_margin_pt <- max(top_margin_pt, 34 + 22 * ang_fac0 + 4 * extra_n)
  }
  left_margin_pt <- if (!is.na(panel.space.left)) .manhattan_panel_space_to_pt(panel.space.left) else 8
  axis_line_col <- "black"
  axis_line_lwd <- 0.35
  panel_box_lwd <- 0.35
  axis_line_x_elem <- if (isTRUE(line.x)) ggplot2::element_line(colour = axis_line_col, linewidth = axis_line_lwd) else ggplot2::element_blank()
  axis_line_y_elem <- ggplot2::element_line(colour = axis_line_col, linewidth = axis_line_lwd)
  panel_border_elem <- if (isTRUE(panel.box)) ggplot2::element_rect(fill = NA, colour = "grey20", linewidth = panel_box_lwd) else ggplot2::element_blank()
  grid_major_x_elem <- if (isTRUE(grid.major.x)) ggplot2::element_line(color = "grey90", linewidth = 0.35) else ggplot2::element_blank()
  grid_major_y_elem <- if (isTRUE(grid.major.y)) ggplot2::element_line(color = "grey90", linewidth = 0.35) else ggplot2::element_blank()
  grid_minor_x_elem <- if (isTRUE(grid.minor.x)) ggplot2::element_line(color = "grey90", linewidth = 0.25) else ggplot2::element_blank()
  grid_minor_y_elem <- if (isTRUE(grid.minor.y)) ggplot2::element_line(color = "grey90", linewidth = 0.25) else ggplot2::element_blank()

  draw_parts <- .manhattan_partition_draw_data(
    dt = dt,
    use_ceiling = use_ceiling,
    direction = direction,
    dir_row_idx = dir_row_idx,
    base_shape = base_shape,
    ceil_shape = ceil_shape,
    base_size = base_size,
    ceil_size = ceil_size,
    ceiling_nudge_fun = .manhattan_ceiling_nudge
  )
  dt_base <- draw_parts$dt_base
  dt_highlight <- draw_parts$dt_highlight
  dt_ceiling_plain <- draw_parts$dt_ceiling_plain
  .manhattan_is_valid_color <- function(v) {
    vv <- as.character(v)
    ok <- !is.na(vv) & nzchar(vv)
    if (any(ok)) {
      ok_idx <- which(ok)
      ok_val <- vapply(vv[ok_idx], function(cc) {
        !inherits(try(grDevices::col2rgb(cc), silent = TRUE), "try-error")
      }, logical(1))
      ok[ok_idx] <- ok_val
    }
    ok
  }

  p <- .manhattan_init_plot(
    dt_base = dt_base,
    x_axis_title = x_axis_title,
    x_breaks = x_breaks,
    x_labels = x_labels,
    y_axis_title = y_axis_title,
    y_breaks_draw = y_breaks_draw,
    y_breaks_labels = y_breaks_labels,
    x_pad = x_pad,
    x_max = x_max,
    y_top = y_top,
    title = title,
    axis_line_x_elem = axis_line_x_elem,
    axis_line_y_elem = axis_line_y_elem,
    tick.size = tick.size,
    panel_border_elem = panel_border_elem,
    grid_major_x_elem = grid_major_x_elem,
    grid_major_y_elem = grid_major_y_elem,
    grid_minor_x_elem = grid_minor_x_elem,
    grid_minor_y_elem = grid_minor_y_elem,
    x.text.size = x.text.size,
    y.text.size = y.text.size,
    x.title.nudge = x.title.nudge,
    top_margin_pt = top_margin_pt,
    bottom_margin_pt = bottom_margin_pt,
    left_margin_pt = left_margin_pt
  )
  p <- .manhattan_draw_chr_ticks(
    p = p,
    chr_map = chr_map,
    special_chr_ticks = special_chr_ticks,
    drag_chrom_tick = drag_chrom_tick,
    tick_len_auto = tick_len_auto,
    linewidth = 0.5,
    color = "grey20"
  )
  p <- .manhattan_draw_highlight_points(p, dt_highlight)
  p <- .manhattan_draw_ceiling_points(
    p = p,
    dt_ceiling_plain = dt_ceiling_plain,
    direction = direction,
    direction_color_inherit = direction_color_inherit,
    dir_col_pos = dir_col_pos,
    dir_col_neg = dir_col_neg,
    ceil_shape = ceil_shape,
    ceil_size = ceil_size
  )
  p <- .manhattan_draw_reference_lines(
    p = p,
    line_dt = line_dt,
    line.alpha = line.alpha,
    y.rescale.line = y.rescale.line,
    use_y_rescale = use_y_rescale,
    y_rescale_at = y_rescale_at,
    y_map_fun = .manhattan_y_map,
    y.rescale.line.linewidth = y.rescale.line.linewidth,
    y.rescale.line.color = y.rescale.line.color,
    y.rescale.line.type = y.rescale.line.type,
    y.ceiling.line = y.ceiling.line,
    use_ceiling = use_ceiling,
    y_ceiling0 = y_ceiling0,
    panel.box = panel.box,
    y.ceiling.line.linewidth = y.ceiling.line.linewidth,
    y.ceiling.line.color = y.ceiling.line.color,
    y.ceiling.line.type = y.ceiling.line.type
  )

  p <- .manhattan_draw_leads(
    p = p,
    lead_tbl = lead_tbl,
    dir_dt = dir_dt,
    dir_row_idx = dir_row_idx,
    suppress_explicit_lead = suppress_explicit_lead,
    lead_color_enabled = lead_color_enabled,
    lead_col0 = lead_col0,
    lead.stroke = lead.stroke,
    ceiling_nudge_fun = .manhattan_ceiling_nudge,
    validate_color_fn = .manhattan_is_valid_color,
    x_max = x_max,
    x_pad = x_pad,
    y_top = y_top,
    y_span = y_span,
    drag.label = drag.label,
    drag.label.line.color = drag.label.line.color,
    drag.label.arrow = drag.label.arrow,
    drag.label.linewidth = drag.label.linewidth,
    lead.label.angle = lead.label.angle,
    lead.label.size = lead.label.size,
    label_hjust = label_hjust,
    label_vjust = label_vjust,
    is_angle_90 = is_angle_90,
    is_angle_270 = is_angle_270,
    is_angle_repel_safe = is_angle_repel_safe,
    is_angle_zero = is_angle_zero,
    label_nudge_base = label_nudge_base,
    label_nudge_step = label_nudge_step,
    seed = seed
  )

  p <- .manhattan_draw_direction_points(
    p = p,
    dir_dt = dir_dt,
    validate_color_fn = .manhattan_is_valid_color,
    direction_size0 = direction_size0,
    lead.stroke = lead.stroke
  )
  p <- .manhattan_draw_overlays(
    p = p,
    y_span = y_span,
    point.size = point.size,
    chr_map = chr_map,
    special_chr_ticks = special_chr_ticks,
    drag_chrom_tick = drag_chrom_tick,
    tick_len_auto = tick_len_auto
  )

  lead_missing <- 0L
  if (nrow(lead_tbl)) {
    lead_missing <- as_int(sum((is.na(lead_tbl$x_data) | !is.finite(lead_tbl$x_data)) & !is.na(lead_tbl$snp)))
    if (lead_missing > 0) .gcanvas_warn_msg(sprintf("Lead SNPs not found in data: n=%d", lead_missing))
  }

  attr(p, "gcanvas_meta") <- list(
    type = "manhattan",
    build = as_int(build),
    columns = list(snp.col = snp_col_use, chrom.col = chrom_col_use, pos.col = pos_col_use, p.col = p.col, y.col = y.col),
    n_variants = as_int(nrow(dt)),
    n_chromosomes = as_int(nrow(chr_map)),
    chromosome_bounds = chr_map[, .(chr = CHR, start = chr_start, end = chr_end)],
    highlights = list(
      n_lead = as_int(nrow(lead_tbl)),
      n_lead_missing = as_int(lead_missing),
      lead_flank_bp = if (is.finite(flank_bp)) as.numeric(flank_bp) else NA_real_,
      n_flank_highlighted = as_int(sum(dt$plot_color != dt$chr_plot_color, na.rm = TRUE))
    )
  )
  .gcanvas_note("gcanvas::manhattan", "Done", silent = silent)
  p
}

.circos_has_color <- function(x) {
  x <- as.character(x)
  vapply(x, function(z) {
    !is.na(z) && nzchar(z) && !is.null(tryCatch(grDevices::col2rgb(z), error = function(e) NULL))
  }, logical(1))
}

.circos_add_alpha <- function(col, alpha = 1) {
  col <- as.character(col)
  alpha <- suppressWarnings(as.numeric(alpha))
  if (!length(alpha)) alpha <- 1
  alpha <- rep_len(alpha, length(col))
  alpha[!is.finite(alpha) | is.na(alpha)] <- 1
  alpha <- pmax(0, pmin(1, alpha))
  out <- col
  ok <- .circos_has_color(col)
  if (any(ok)) {
    out[ok] <- vapply(which(ok), function(i) {
      grDevices::adjustcolor(col[i], alpha.f = alpha[i])
    }, character(1))
  }
  out
}

.circos_map_values_by_name <- function(keys, values, default = NA) {
  keys <- as.character(keys)
  out <- rep_len(default, length(keys))
  if (is.null(values) || length(values) == 0L) return(out)
  vals <- as.character(values)
  nms <- names(values)
  if (!is.null(nms) && any(!is.na(nms) & nzchar(nms))) {
    hit <- match(keys, nms)
    ok <- is.finite(hit) & !is.na(hit) & hit > 0L
    out[ok] <- vals[hit[ok]]
    return(out)
  }
  if (length(vals) == 1L) return(rep_len(vals, length(keys)))
  if (length(vals) == length(keys)) return(vals)
  rep_len(vals, length(keys))
}

.circos_track_match_index <- function(key, tracks) {
  tracks <- as.character(tracks)
  key <- as.character(key)[1]
  if (!length(tracks) || is.na(key) || !nzchar(key)) return(NA_integer_)
  hit <- match(key, tracks)
  if (!is.na(hit) && hit > 0L) return(as_int(hit))
  idx <- suppressWarnings(as.integer(key))[1]
  if (is.finite(idx) && !is.na(idx) && idx >= 1L && idx <= length(tracks)) return(as_int(idx))
  NA_integer_
}

.circos_track_scalar_map <- function(value, tracks, default = NULL, numeric = FALSE) {
  tracks <- as.character(tracks)
  n <- length(tracks)
  if (!n) return(setNames(vector(if (isTRUE(numeric)) "numeric" else "character", 0L), character()))
  base_default <- if (length(default)) rep_len(default, n) else rep_len(NA, n)
  if (isTRUE(numeric)) {
    base_default <- suppressWarnings(as.numeric(base_default))
    out <- base_default
  } else {
    out <- as.character(base_default)
  }
  names(out) <- tracks
  if (is.null(value) || length(value) == 0L) return(out)

  if (is.list(value) && !is.data.frame(value) && !data.table::is.data.table(value)) {
    nms <- names(value)
    has_names <- !is.null(nms) && any(!is.na(nms) & nzchar(nms))
    if (has_names) {
      for (i in seq_along(value)) {
        idx <- .circos_track_match_index(nms[i], tracks)
        if (is.na(idx) || idx <= 0L) next
        vi <- value[[i]]
        if (is.null(vi) || !length(vi)) next
        out[idx] <- if (isTRUE(numeric)) suppressWarnings(as.numeric(vi[1])) else as.character(vi[1])
      }
      return(out)
    }
    value <- unlist(value, use.names = FALSE)
  }

  vals <- if (isTRUE(numeric)) suppressWarnings(as.numeric(value)) else as.character(value)
  nms <- names(value)
  if (!is.null(nms) && any(!is.na(nms) & nzchar(nms))) {
    for (i in seq_along(vals)) {
      idx <- .circos_track_match_index(nms[i], tracks)
      if (is.na(idx) || idx <= 0L) next
      out[idx] <- vals[i]
    }
    return(out)
  }
  if (length(vals) == 1L) {
    out[] <- vals[1]
  } else if (length(vals) == n) {
    out[] <- vals
  } else {
    out[] <- rep_len(vals, n)
  }
  names(out) <- tracks
  out
}

.circos_track_list_map <- function(value, tracks, default = NULL) {
  tracks <- as.character(tracks)
  out <- stats::setNames(vector("list", length(tracks)), tracks)
  if (!length(tracks)) return(out)
  if (is.null(value) || length(value) == 0L) {
    for (i in seq_along(tracks)) out[[i]] <- default
    return(out)
  }
  if (is.list(value) && !is.data.frame(value) && !data.table::is.data.table(value)) {
    nms <- names(value)
    if (!is.null(nms) && any(!is.na(nms) & nzchar(nms))) {
      for (i in seq_along(tracks)) out[[i]] <- default
      for (i in seq_along(value)) {
        idx <- .circos_track_match_index(nms[i], tracks)
        if (is.na(idx) || idx <= 0L) next
        out[[idx]] <- value[[i]]
      }
      return(out)
    }
    if (length(value) == length(tracks)) {
      for (i in seq_along(tracks)) out[[i]] <- value[[i]]
      names(out) <- tracks
      return(out)
    }
    if (length(value) == 1L) {
      for (i in seq_along(tracks)) out[[i]] <- value[[1]]
      return(out)
    }
  }
  nms <- names(value)
  if (!is.null(nms) && any(!is.na(nms) & nzchar(nms))) {
    for (i in seq_along(tracks)) out[[i]] <- default
    for (i in seq_along(value)) {
      idx <- .circos_track_match_index(nms[i], tracks)
      if (is.na(idx) || idx <= 0L) next
      out[[idx]] <- value[[i]]
    }
    return(out)
  }
  for (i in seq_along(tracks)) out[[i]] <- value
  out
}

.circos_track_bp_map <- function(value, tracks, default = NA_real_, arg_name = "value") {
  if (is.null(value) || length(value) == 0L || (is.logical(value) && length(value) == 1L)) {
    out0 <- rep_len(default, length(tracks))
    names(out0) <- as.character(tracks)
    return(out0)
  }
  raw_map <- .circos_track_scalar_map(value, tracks, default = default, numeric = FALSE)
  out <- suppressWarnings(as.numeric(raw_map))
  bad <- !is.finite(out) | is.na(out)
  if (any(bad)) {
    out[bad] <- vapply(raw_map[bad], function(x) .gcanvas_parse_bp_span(x, arg_name = arg_name), numeric(1))
  }
  out[!is.finite(out) | is.na(out) | out <= 0] <- NA_real_
  names(out) <- names(raw_map)
  out
}

.circos_apply_y_rescale <- function(v, at = NA_real_, ratio = 0.25) {
  vv <- suppressWarnings(as.numeric(v))
  at <- suppressWarnings(as.numeric(at))[1]
  ratio <- suppressWarnings(as.numeric(ratio))[1]
  if (!is.finite(at) || is.na(at) || at <= 0) return(vv)
  if (!is.finite(ratio) || is.na(ratio) || ratio <= 0 || ratio > 1) ratio <- 0.25
  out <- vv
  idx <- is.finite(vv) & !is.na(vv) & vv > at
  out[idx] <- at + (vv[idx] - at) * ratio
  out
}

.circos_nice_y_max <- function(y_max_raw, user_specified = FALSE) {
  x <- suppressWarnings(as.numeric(y_max_raw))[1]
  if (!is.finite(x) || is.na(x) || x <= 0) return(1)
  if (isTRUE(user_specified)) {
    return(max(1, ceiling(x)))
  }
  max(1, ceiling(x * 1.05))
}

.circos_compute_y_breaks <- function(y_raw,
                                     y_max_raw,
                                     y_breaks = NULL,
                                     y_rescale.at = NULL,
                                     y_rescale.ratio = 0.25,
                                     y_rescale.breaks = NULL) {
  y_raw <- suppressWarnings(as.numeric(y_raw))
  y_raw <- y_raw[is.finite(y_raw) & !is.na(y_raw)]
  y_max_raw <- suppressWarnings(as.numeric(y_max_raw))[1]
  if (!is.finite(y_max_raw) || is.na(y_max_raw) || y_max_raw <= 0) {
    y_max_raw <- if (length(y_raw)) max(y_raw, na.rm = TRUE) else 1
  }
  y_rescale.at <- suppressWarnings(as.numeric(y_rescale.at))[1]
  use_rescale <- is.finite(y_rescale.at) && !is.na(y_rescale.at) && y_rescale.at > 0
  if (is.null(y_breaks) || (is.character(y_breaks) && length(y_breaks) == 1L &&
                            tolower(trimws(as.character(y_breaks)[1])) == "auto")) {
    if (isTRUE(use_rescale)) {
      br_low <- suppressWarnings(as.numeric(base::pretty(c(0, y_rescale.at), n = 3)))
      br_low <- br_low[is.finite(br_low) & br_low <= y_rescale.at]
      br_hi_src <- suppressWarnings(as.numeric(y_rescale.breaks))
      br_hi_src <- br_hi_src[is.finite(br_hi_src) & br_hi_src >= y_rescale.at]
      if (!length(br_hi_src)) {
        br_hi_src <- suppressWarnings(as.numeric(base::pretty(c(y_rescale.at, y_max_raw), n = 2)))
        br_hi_src <- br_hi_src[is.finite(br_hi_src) & br_hi_src >= y_rescale.at]
      }
      br_raw <- sort(unique(c(0, br_low, y_rescale.at, br_hi_src)))
      if (length(br_raw) > 6L) {
        lo <- br_raw[br_raw < y_rescale.at]
        hi <- br_raw[br_raw > y_rescale.at]
        lo_keep <- if (length(lo) > 2L) lo[unique(as_int(round(seq(1, length(lo), length.out = 2L))))] else lo
        hi_keep <- if (length(hi) > 2L) hi[unique(as_int(round(seq(1, length(hi), length.out = 2L))))] else hi
        br_raw <- sort(unique(c(0, lo_keep, y_rescale.at, hi_keep)))
      }
    } else {
      br_raw <- suppressWarnings(as.numeric(base::pretty(c(0, y_max_raw), n = 4)))
      br_raw <- br_raw[is.finite(br_raw) & br_raw >= 0]
      if (!length(br_raw)) br_raw <- c(0, y_max_raw)
    }
  } else {
    br_raw <- suppressWarnings(as.numeric(y_breaks))
    br_raw <- br_raw[is.finite(br_raw) & br_raw >= 0]
    if (!length(br_raw)) br_raw <- 0
  }
  br_raw <- sort(unique(c(0, br_raw)))
  br_raw <- br_raw[br_raw <= y_max_raw]
  if (!length(br_raw)) br_raw <- c(0, y_max_raw)
  if (!isTRUE(use_rescale) && length(br_raw) > 5L) {
    keep <- unique(as_int(round(seq(1, length(br_raw), length.out = 5L))))
    br_raw <- br_raw[keep]
  }
  br_draw <- .circos_apply_y_rescale(br_raw, at = y_rescale.at, ratio = y_rescale.ratio)
  data.table::data.table(
    y_raw = br_raw,
    y_draw = br_draw,
    label = as.character(signif(br_raw, 4))
  )
}

.circos_arc_xy <- function(start.degree, end.degree, radius, n = 120L) {
  ang <- seq(as.numeric(start.degree), as.numeric(end.degree), length.out = max(2L, as_int(n)))
  rad <- ang * pi / 180
  data.frame(
    x = radius * cos(rad),
    y = radius * sin(rad)
  )
}

.circos_ring_polygon <- function(start.degree, end.degree, r.inner, r.outer, n = 120L) {
  a1 <- seq(as.numeric(start.degree), as.numeric(end.degree), length.out = max(2L, as_int(n)))
  a2 <- rev(a1)
  r1 <- a1 * pi / 180
  r2 <- a2 * pi / 180
  data.frame(
    x = c(r.outer * cos(r1), r.inner * cos(r2)),
    y = c(r.outer * sin(r1), r.inner * sin(r2))
  )
}

.circos_build_shadow_regions <- function(shadow_tbl,
                                         chr_map,
                                         track_meta,
                                         tracks,
                                         shadow.alpha.map) {
  require_pkg("data.table")
  empty <- data.table::data.table(
    track = character(),
    CHR = character(),
    start.degree = numeric(),
    end.degree = numeric(),
    r.inner = numeric(),
    r.outer = numeric(),
    fill = character(),
    alpha = numeric()
  )
  if (!nrow(shadow_tbl)) return(empty)
  shadow_use <- data.table::copy(shadow_tbl)
  shadow_use <- shadow_use[is.finite(flank_bp) & flank_bp > 0 & !is.na(CHR) & nzchar(CHR) & is.finite(POS)]
  if (!nrow(shadow_use)) return(empty)
  shadow_use <- merge(
    shadow_use,
    chr_map[, .(CHR, chr_end, start.degree.chr = start.degree, width.degree)],
    by = "CHR",
    all.x = TRUE,
    sort = FALSE
  )
  shadow_use <- shadow_use[is.finite(chr_end) & is.finite(start.degree.chr) & is.finite(width.degree)]
  if (!nrow(shadow_use)) return(empty)

  parts <- vector("list", nrow(shadow_use))
  for (i in seq_len(nrow(shadow_use))) {
    li <- shadow_use[i]
    tr_use <- if (!is.na(li$track[[1]]) && nzchar(li$track[[1]])) as.character(li$track[[1]]) else as.character(tracks)
    tr_use <- tr_use[tr_use %in% track_meta$track]
    if (!length(tr_use)) next
    start.bp <- max(1, li$POS[[1]] - li$flank_bp[[1]])
    end.bp <- min(li$chr_end[[1]], li$POS[[1]] + li$flank_bp[[1]])
    if (!is.finite(start.bp) || !is.finite(end.bp) || end.bp < start.bp) next
    frac.start <- if (li$chr_end[[1]] <= 1) 0 else (start.bp - 1) / li$chr_end[[1]]
    frac.end <- if (li$chr_end[[1]] <= 1) 0 else (end.bp - 1) / li$chr_end[[1]]
    frac.start <- pmax(0, pmin(1, frac.start))
    frac.end <- pmax(0, pmin(1, frac.end))
    start.deg <- li$start.degree.chr[[1]] - frac.start * li$width.degree[[1]]
    end.deg <- li$start.degree.chr[[1]] - frac.end * li$width.degree[[1]]
    tm <- track_meta[match(tr_use, track)]
    alpha.now <- as.numeric(shadow.alpha.map[match(tr_use, names(shadow.alpha.map))])
    alpha.now[!is.finite(alpha.now) | is.na(alpha.now)] <- 1
    alpha.now <- pmax(0, pmin(1, alpha.now))
    parts[[i]] <- data.table::data.table(
      track = tr_use,
      CHR = as.character(li$CHR[[1]]),
      start.degree = as.numeric(start.deg),
      end.degree = as.numeric(end.deg),
      r.inner = tm$r.inner,
      r.outer = tm$r.outer,
      fill = as.character(li$flank_color[[1]]),
      alpha = alpha.now
    )
  }
  parts <- Filter(Negate(is.null), parts)
  if (!length(parts)) return(empty)
  out <- data.table::rbindlist(parts, use.names = TRUE, fill = TRUE)
  out[is.na(fill) | !nzchar(fill), fill := "#E63946"]
  out
}

.circos_prepare_shadow_table <- function(lead.shadow,
                                         lead_tbl,
                                         tracks,
                                         lead.color = "#E63946",
                                         lead.flank = NULL,
                                         lead.flank.color = NULL,
                                         lead.shadow.flank = NULL,
                                         lead.shadow.color = NULL) {
  require_pkg("data.table")
  empty <- data.table::data.table(
    track = character(),
    CHR = character(),
    POS = numeric(),
    flank_bp = numeric(),
    flank_color = character()
  )
  if (is.null(lead.shadow) || length(lead.shadow) == 0L || (is.logical(lead.shadow) && length(lead.shadow) == 1L && !isTRUE(lead.shadow))) {
    return(empty)
  }

  shadow_tbl <- NULL
  if (is.logical(lead.shadow) && length(lead.shadow) == 1L && isTRUE(lead.shadow)) {
    if (!nrow(lead_tbl)) return(empty)
    shadow_tbl <- data.table::copy(lead_tbl[, .(track, CHR, POS, flank_bp, flank_color, lead_color)])
  } else {
    shadow_tbl <- .circos_normalize_lead_table(
      lead = lead.shadow,
      tracks = tracks,
      lead.color = lead.shadow.color %||% lead.color,
      lead.size = NULL,
      lead.stroke = NULL,
      lead.flank = lead.shadow.flank %||% lead.flank,
      lead.flank.color = lead.shadow.color %||% lead.flank.color,
      lead.label = FALSE,
      lead.label.color = "grey20"
    )[, .(track, CHR, POS, flank_bp, flank_color, lead_color)]
  }
  if (!nrow(shadow_tbl)) return(empty)

  shadow_flank_map <- .circos_track_bp_map(lead.shadow.flank, tracks, default = NA_real_, arg_name = "lead.shadow.flank")
  lead_flank_map <- .circos_track_bp_map(lead.flank, tracks, default = NA_real_, arg_name = "lead.flank")
  shadow_color_map <- .circos_track_scalar_map(lead.shadow.color, tracks, default = NA_character_)
  lead_flank_color_map <- .circos_track_scalar_map(lead.flank.color, tracks, default = NA_character_)
  lead_color_map <- .circos_track_scalar_map(lead.color, tracks, default = "#E63946")

  map_track_default <- function(track_vec, mp, fallback = NA) {
    out0 <- rep_len(fallback, length(track_vec))
    hit <- match(as.character(track_vec), names(mp))
    ok <- is.finite(hit) & !is.na(hit) & hit > 0L
    out0[ok] <- mp[hit[ok]]
    if (length(mp)) out0[!ok] <- rep_len(mp[1], sum(!ok))
    out0
  }

  shadow_tbl[is.na(flank_bp), flank_bp := as.numeric(map_track_default(track, shadow_flank_map, NA_real_))]
  shadow_tbl[is.na(flank_bp), flank_bp := as.numeric(map_track_default(track, lead_flank_map, NA_real_))]
  shadow_tbl[is.na(flank_color) | !nzchar(flank_color), flank_color := as.character(map_track_default(track, shadow_color_map, NA_character_))]
  shadow_tbl[is.na(flank_color) | !nzchar(flank_color), flank_color := as.character(map_track_default(track, lead_flank_color_map, NA_character_))]
  shadow_tbl[is.na(flank_color) | !nzchar(flank_color), flank_color := as.character(map_track_default(track, lead_color_map, "#E63946"))]
  shadow_tbl[!is.finite(flank_bp) | flank_bp <= 0, flank_bp := NA_real_]
  shadow_tbl[!is.na(CHR) & nzchar(CHR) & is.finite(POS) & is.finite(flank_bp) & flank_bp > 0]
}

.circos_prepare_outer_band <- function(dt,
                                       lead_tbl,
                                       chr_map,
                                       outer.band = FALSE,
                                       outer.band.bin.size = "1Mb",
                                       outer.band.metric = "auto",
                                       outer.band.color = "auto",
                                       outer.band.zero.color = "grey90",
                                       band.name = "outer.band") {
  require_pkg("data.table")
  bin_arg_name <- paste0(band.name, ".bin.size")
  metric_arg_name <- paste0(band.name, ".metric")
  resolve_group_value_map <- function(value, groups, default = NULL, bp = FALSE) {
    out <- rep_len(default, length(groups))
    names(out) <- groups
    if (is.null(value) || length(value) == 0L) {
      if (isTRUE(bp)) {
        out_num <- suppressWarnings(as.numeric(out))
        names(out_num) <- groups
        return(out_num)
      }
      return(out)
    }
    assign_one <- function(key, val) {
      idx <- suppressWarnings(as.integer(key))
      if (is.finite(idx) && !is.na(idx) && idx >= 1L && idx <= length(groups)) {
        out[idx] <<- val
      } else if (!is.na(key) && nzchar(key)) {
        hit <- match(as.character(key), groups)
        if (is.finite(hit) && !is.na(hit) && hit > 0L) out[hit] <<- val
      }
    }
    if (is.list(value) && !is.data.frame(value) && !data.table::is.data.table(value)) {
      vn <- names(value)
      if (is.null(vn) || !any(!is.na(vn) & nzchar(vn))) {
        for (i in seq_len(min(length(value), length(groups)))) out[i] <- unlist(value[i], recursive = TRUE, use.names = FALSE)[1]
      } else {
        for (i in seq_along(value)) assign_one(vn[i], unlist(value[i], recursive = TRUE, use.names = FALSE)[1])
      }
    } else {
      vv <- value
      vn <- names(vv)
      if (!is.null(vn) && any(!is.na(vn) & nzchar(vn))) {
        for (i in seq_along(vv)) assign_one(vn[i], vv[i])
      } else if (length(vv) == 1L) {
        out[] <- vv[1]
      } else {
        out[seq_len(min(length(vv), length(groups)))] <- vv[seq_len(min(length(vv), length(groups)))]
      }
    }
    if (isTRUE(bp)) {
      out_num <- vapply(as.list(out), function(z) .gcanvas_parse_bp_span(z, arg_name = bin_arg_name), numeric(1))
      out_num[!is.finite(out_num) | is.na(out_num) | out_num <= 0] <- .gcanvas_parse_bp_span(default, arg_name = bin_arg_name)
      names(out_num) <- groups
      return(out_num)
    }
    out
  }
  parse_outer_band_groups <- function(spec, default.group.colors, default.bin.size) {
    if (!(is.list(spec) && !is.data.frame(spec) && !data.table::is.data.table(spec))) return(NULL)
    parts <- spec
    if (!length(parts)) return(NULL)
    gnames <- names(parts)
    if (is.null(gnames) || !length(gnames)) gnames <- rep("", length(parts))
    blank <- which(is.na(gnames) | !nzchar(gnames))
    if (length(blank)) gnames[blank] <- as.character(blank)
    gnames <- make.unique(gnames, sep = ".")
    snp_map <- unique(dt[, .(snp, CHR, POS)])
    group_color_map <- resolve_group_value_map(default.group.colors, gnames, default = NA_character_, bp = FALSE)
    group_bin_map <- resolve_group_value_map(default.bin.size, gnames, default = default.bin.size, bp = TRUE)
    default_pal <- tryCatch(get.colors(mode = "lagoon", n = length(gnames), plot = FALSE, silent = TRUE), error = function(e) NULL)
    if (!length(default_pal)) default_pal <- grDevices::rainbow(length(gnames))

    parse_group_df <- function(x) {
      x <- data.table::as.data.table(x)
      nms <- names(x)
      low <- tolower(nms)
      pick <- function(cands) {
        hit <- match(cands, low, nomatch = 0L)
        hit <- hit[hit > 0L]
        if (!length(hit)) return(NA_character_)
        nms[hit[1]]
      }
      snp_col <- pick(c("snp", "snpid", "rsid", "variant", "variant_id"))
      chr_col <- pick(c("chr", "chrom", "chromosome"))
      pos_col <- pick(c("pos", "bp", "position"))
      data.table::data.table(
        snp = if (!is.na(snp_col)) as.character(x[[snp_col]]) else NA_character_,
        CHR = if (!is.na(chr_col)) normalize.chrom(x[[chr_col]]) else NA_character_,
        POS = if (!is.na(pos_col)) suppressWarnings(as.numeric(x[[pos_col]])) else NA_real_
      )
    }

    out_parts <- vector("list", length(parts))
    for (i in seq_along(parts)) {
      gi <- gnames[i]
      xi <- parts[[i]]
      if (is.null(xi) || length(xi) == 0L) next
      if (is.data.frame(xi) || data.table::is.data.table(xi)) {
        dt0 <- parse_group_df(xi)
      } else {
        snps <- .gcanvas_as_snp_vector(xi)
        if (!length(snps)) next
        dt0 <- data.table::data.table(snp = snps, CHR = NA_character_, POS = NA_real_)
      }
      if (!nrow(dt0)) next
      dt0 <- merge(dt0, snp_map, by = "snp", all.x = TRUE, sort = FALSE)
      if ("CHR.x" %in% names(dt0)) {
        dt0[is.na(CHR.x) | !nzchar(CHR.x), CHR.x := CHR.y]
        dt0[, CHR := normalize.chrom(CHR.x)]
        dt0[, c("CHR.x", "CHR.y") := NULL]
      } else if ("CHR" %in% names(dt0)) {
        dt0[, CHR := normalize.chrom(CHR)]
      }
      if ("POS.x" %in% names(dt0)) {
        dt0[!is.finite(POS.x) | is.na(POS.x), POS.x := POS.y]
        dt0[, POS := suppressWarnings(as.numeric(POS.x))]
        dt0[, c("POS.x", "POS.y") := NULL]
      } else if ("POS" %in% names(dt0)) {
        dt0[, POS := suppressWarnings(as.numeric(POS))]
      }
      dt0 <- dt0[!is.na(CHR) & nzchar(CHR) & is.finite(POS) & !is.na(POS) & POS > 0]
      if (!nrow(dt0)) next
      dt0[, `:=`(
        group = gi,
        group_color = as.character(group_color_map[gi]),
        flank_bp = as.numeric(group_bin_map[gi])
      )]
      if (is.na(dt0$group_color[1]) || !nzchar(dt0$group_color[1])) dt0[, group_color := default_pal[i]]
      out_parts[[i]] <- unique(dt0[, .(group, snp, CHR, POS, group_color, flank_bp)])
    }
    out_parts <- Filter(Negate(is.null), out_parts)
    if (!length(out_parts)) return(NULL)
    data.table::rbindlist(out_parts, use.names = TRUE, fill = TRUE)
  }
  resolve_outer_palette <- function(x, n = 100L) {
    x0 <- as.character(x)
    x0 <- x0[!is.na(x0) & nzchar(x0)]
    if (!length(x0)) return(character())
    if (length(x0) == 1L) {
      key <- tolower(trimws(x0[1]))
      chrom_pal <- .gcanvas_chrom_palette_list()
      if (key %in% names(chrom_pal)) {
        return(grDevices::colorRampPalette(as.character(chrom_pal[[key]]), space = "Lab")(n))
      }
      pal_try <- tryCatch(get.colors(mode = x0[1], n = n, plot = FALSE, silent = TRUE), error = function(e) NULL)
      if (length(pal_try)) return(as.character(pal_try))
    }
    x0
  }
  map_counts_to_palette <- function(n_vec, pal, zero_col = "grey85") {
    out <- rep_len(zero_col, length(n_vec))
    pos <- which(is.finite(n_vec) & !is.na(n_vec) & n_vec > 0)
    if (!length(pos) || !length(pal)) return(out)
    n_pos <- n_vec[pos]
    uq <- sort(unique(n_pos))
    if (length(uq) == 1L) {
      idx0 <- max(1L, min(length(pal), as.integer(round(length(pal) * 0.6))))
      out[pos] <- pal[idx0]
      return(out)
    }
    idx_map <- round(seq(1, length(pal), length.out = length(uq)))
    idx_map <- pmax(1L, pmin(length(pal), as.integer(idx_map)))
    out[pos] <- pal[idx_map[match(n_pos, uq)]]
    out
  }
  empty <- data.table::data.table(
    CHR = character(),
    bin.id = integer(),
    start = numeric(),
    end = numeric(),
    n = integer(),
    start.degree = numeric(),
    end.degree = numeric(),
    fill = character()
  )
  outer.band.enabled <- !(is.null(outer.band) || length(outer.band) == 0L || (is.logical(outer.band) && length(outer.band) == 1L && !isTRUE(outer.band)))
  if (!isTRUE(outer.band.enabled)) return(empty)

  bin_bp <- .gcanvas_parse_bp_span(outer.band.bin.size, arg_name = bin_arg_name)
  if (!is.finite(bin_bp) || is.na(bin_bp) || bin_bp <= 0) {
    stop(sprintf("%s must be numeric (bp) or string like '250kb'/'1mb'.", bin_arg_name), call. = FALSE)
  }

  metric0 <- tolower(trimws(as.character(outer.band.metric)[1]))
  if (is.na(metric0) || !nzchar(metric0)) metric0 <- "auto"
  if (!metric0 %in% c("auto", "variant", "lead")) {
    stop(sprintf("%s must be one of 'auto', 'variant', or 'lead'.", metric_arg_name), call. = FALSE)
  }
  if (identical(metric0, "auto")) metric0 <- if (nrow(lead_tbl)) "lead" else "variant"
  color_mode <- tolower(trimws(as.character(outer.band.color)[1]))
  if (is.na(color_mode) || !nzchar(color_mode)) color_mode <- "auto"
  use_lead_color <- identical(color_mode, "auto") && identical(metric0, "lead")
  pal <- outer.band.color
  if (identical(color_mode, "auto")) {
    pal <- .gcanvas_ld_palette(100L)
  } else if (identical(color_mode, "ld")) {
    pal <- .gcanvas_ld_palette(100L)
  } else {
    pal <- resolve_outer_palette(pal, n = 100L)
  }
  pal <- as.character(pal)
  pal <- pal[!is.na(pal) & nzchar(pal)]
  if (!length(pal)) pal <- .gcanvas_ld_palette(100L)
  zero_col <- as.character(outer.band.zero.color)[1]
  if (is.na(zero_col) || !nzchar(zero_col)) zero_col <- "grey85"

  explicit_group_dt <- parse_outer_band_groups(
    spec = outer.band,
    default.group.colors = outer.band.color,
    default.bin.size = outer.band.bin.size
  )
  if (data.table::is.data.table(explicit_group_dt) && nrow(explicit_group_dt)) {
    chr_len <- data.table::copy(chr_map[, .(CHR, chr_end, chr_start.degree = start.degree, chr_end.degree = end.degree, width.degree)])
    explicit_group_dt <- merge(explicit_group_dt, chr_len, by = "CHR", all.x = TRUE, sort = FALSE)
    explicit_group_dt <- explicit_group_dt[is.finite(chr_end) & !is.na(chr_end) & chr_end > 0]
    if (!nrow(explicit_group_dt)) return(empty)
    explicit_group_dt[POS > chr_end, POS := chr_end]
    explicit_group_dt[!is.finite(flank_bp) | is.na(flank_bp) | flank_bp <= 0, flank_bp := bin_bp]
    explicit_group_dt[, `:=`(
      start = pmax(1, POS - floor(flank_bp / 2)),
      end = pmin(chr_end, POS + ceiling(flank_bp / 2))
    )]
    explicit_group_dt[, bin.id := .I]
    explicit_group_dt[, `:=`(
      start.degree = chr_start.degree - ((start - 1) / chr_end) * width.degree,
      end.degree = chr_start.degree - ((end / chr_end) * width.degree),
      n = 1L,
      fill = as.character(group_color)
    )]
    explicit_group_dt[is.na(fill) | !nzchar(fill), fill := zero_col]
    bg_dt <- chr_len[, .(
      CHR = CHR,
      bin.id = 0L,
      start = 1,
      end = chr_end,
      n = 0L,
      start.degree = chr_start.degree,
      end.degree = chr_end.degree,
      fill = zero_col
    )]
    out <- data.table::rbindlist(
      list(bg_dt, explicit_group_dt[, .(CHR, bin.id, start, end, n, start.degree, end.degree, fill)]),
      use.names = TRUE,
      fill = TRUE
    )
    return(out)
  }

  lead_color_dt <- NULL
  if (identical(metric0, "lead")) {
    if (!nrow(lead_tbl)) return(empty)
    src <- data.table::copy(lead_tbl)
    if ("snp" %in% names(src)) {
      src <- unique(src[, .(snp, CHR, POS, lead_color)], by = c("snp", "CHR", "POS", "lead_color"))
    } else {
      src <- unique(src[, .(CHR, POS, lead_color)])
      src[, snp := NA_character_]
    }
  } else {
    src <- data.table::copy(dt)
    if ("snp" %in% names(src)) {
      src <- unique(src[, .(snp, CHR, POS)], by = c("snp", "CHR", "POS"))
    } else {
      src <- unique(src[, .(CHR, POS)])
      src[, snp := NA_character_]
    }
  }
  src <- src[!is.na(CHR) & nzchar(CHR) & is.finite(POS) & !is.na(POS) & POS > 0]
  if (!nrow(src)) return(empty)

  chr_len <- data.table::copy(chr_map[, .(CHR, chr_end, chr_start.degree = start.degree, chr_end.degree = end.degree, width.degree)])
  src <- merge(src, chr_len, by = "CHR", all.x = TRUE, sort = FALSE)
  src <- src[is.finite(chr_end) & !is.na(chr_end) & chr_end > 0]
  if (!nrow(src)) return(empty)
  src[POS > chr_end, POS := chr_end]

  if (identical(metric0, "lead")) {
    half_left <- floor(bin_bp / 2)
    half_right <- ceiling(bin_bp / 2)
    src[, `:=`(
      start = pmax(1, POS - half_left),
      end = pmin(chr_end, POS + half_right)
    )]
    win_dt <- src[, .(
      n = .N,
      lead_color = {
        v <- lead_color[!is.na(lead_color) & nzchar(lead_color)]
        if (!length(v)) NA_character_ else names(sort(table(v), decreasing = TRUE))[1]
      },
      chr_end = chr_end[1],
      chr_start.degree = chr_start.degree[1],
      chr_end.degree = chr_end.degree[1],
      width.degree = width.degree[1]
    ), by = .(CHR, start, end)]
    win_dt[, bin.id := .I]
    win_dt[, `:=`(
      start.degree = chr_start.degree - ((start - 1) / chr_end) * width.degree,
      end.degree = chr_start.degree - ((end / chr_end) * width.degree)
    )]
    bg_dt <- chr_len[, .(
      CHR = CHR,
      bin.id = 0L,
      start = 1,
      end = chr_end,
      n = 0L,
      lead_color = NA_character_,
      start.degree = chr_start.degree,
      end.degree = chr_end.degree
    )]
    if (isTRUE(use_lead_color)) {
      win_dt[, fill := ifelse(!is.na(lead_color) & nzchar(lead_color), lead_color, zero_col)]
    } else {
      win_dt[, fill := map_counts_to_palette(n, pal = pal, zero_col = zero_col)]
    }
    bg_dt[, fill := zero_col]
    out <- data.table::rbindlist(list(bg_dt, win_dt[, .(CHR, bin.id, start, end, n, lead_color, start.degree, end.degree, fill)]), use.names = TRUE, fill = TRUE)
    return(out[, .(CHR, bin.id, start, end, n, start.degree, end.degree, fill)])
  }

  src[, bin.id := as.integer(floor((POS - 1) / bin_bp) + 1L)]

  count_dt <- src[, .(n = .N), by = .(CHR, bin.id)]
  if (identical(metric0, "lead") && ("lead_color" %in% names(src))) {
    lead_color_dt <- src[
      !is.na(lead_color) & nzchar(lead_color),
      .N,
      by = .(CHR, bin.id, lead_color)
    ]
    if (nrow(lead_color_dt)) {
      data.table::setorderv(lead_color_dt, c("CHR", "bin.id", "N", "lead_color"), c(1L, 1L, -1L, 1L), na.last = TRUE)
      lead_color_dt <- lead_color_dt[, .SD[1], by = .(CHR, bin.id)]
      lead_color_dt[, N := NULL]
    }
  }

  bin_parts <- vector("list", nrow(chr_len))
  for (i in seq_len(nrow(chr_len))) {
    chr_i <- chr_len$CHR[i]
    chr_end_i <- chr_len$chr_end[i]
    n_bin_i <- max(1L, as.integer(ceiling(chr_end_i / bin_bp)))
    start_i <- ((seq_len(n_bin_i) - 1L) * bin_bp) + 1
    end_i <- pmin(seq_len(n_bin_i) * bin_bp, chr_end_i)
    bin_parts[[i]] <- data.table::data.table(
      CHR = chr_i,
      bin.id = seq_len(n_bin_i),
      start = start_i,
      end = end_i
    )
  }
  bins <- data.table::rbindlist(bin_parts, use.names = TRUE, fill = TRUE)
  bins <- merge(bins, count_dt, by = c("CHR", "bin.id"), all.x = TRUE, sort = FALSE)
  bins[is.na(n), n := 0L]
  if (data.table::is.data.table(lead_color_dt) && nrow(lead_color_dt)) {
    bins <- merge(bins, lead_color_dt, by = c("CHR", "bin.id"), all.x = TRUE, sort = FALSE)
  } else {
    bins[, lead_color := NA_character_]
  }
  bins <- merge(bins, chr_len, by = "CHR", all.x = TRUE, sort = FALSE)
  bins[, `:=`(
    start.degree = chr_start.degree - ((start - 1) / chr_end) * width.degree,
    end.degree = chr_start.degree - ((end / chr_end) * width.degree
    )
  )]
  bins[, fill := zero_col]
  bins[, fill := map_counts_to_palette(n, pal = pal, zero_col = zero_col)]
  bins[]
}

.circos_lty <- function(x) {
  x0 <- tolower(trimws(as.character(x)[1]))
  if (is.na(x0) || !nzchar(x0)) return("solid")
  if (x0 %in% c("dottd", "dot", "dott")) return("dotted")
  if (x0 %in% c("dash", "dashed")) return("dashed")
  if (x0 %in% c("solid", "dotted", "dotdash", "longdash", "twodash")) return(x0)
  x0
}

.circos_open_device <- function(file, width = 20, height = 20, units = "cm", dpi = 300) {
  if (is.null(file) || length(file) == 0L) return(NULL)
  f0 <- as.character(file)[1]
  if (is.na(f0) || !nzchar(f0)) return(NULL)
  ext <- tolower(tools::file_ext(f0))
  if (!nzchar(ext) && tolower(f0) %in% c("png", "pdf", "jpg", "jpeg", "tiff", "bmp")) {
    ext <- tolower(f0)
    f0 <- paste0("circos.", ext)
  }
  if (!nzchar(ext)) stop("file must be NULL, a file path, or one of png/pdf/jpg/jpeg/tiff/bmp.", call. = FALSE)
  width <- suppressWarnings(as.numeric(width))[1]
  height <- suppressWarnings(as.numeric(height))[1]
  units <- tolower(trimws(as.character(units)[1]))
  dpi <- suppressWarnings(as.numeric(dpi))[1]
  if (!is.finite(width) || is.na(width) || width <= 0) width <- 20
  if (!is.finite(height) || is.na(height) || height <= 0) height <- 20
  if (is.na(units) || !nzchar(units)) units <- "cm"
  if (!(units %in% c("cm", "mm", "in", "px"))) units <- "cm"
  if (!is.finite(dpi) || is.na(dpi) || dpi <= 0) dpi <- 300
  to_inches <- function(x, unit, dpi) {
    if (identical(unit, "in")) return(x)
    if (identical(unit, "cm")) return(x / 2.54)
    if (identical(unit, "mm")) return(x / 25.4)
    x / dpi
  }
  to_pixels <- function(x, unit, dpi) {
    if (identical(unit, "px")) return(as_int(round(x)))
    as_int(round(to_inches(x, unit, dpi) * dpi))
  }
  ext <- if (ext == "jpg") "jpeg" else ext
  width_in <- to_inches(width, units, dpi)
  height_in <- to_inches(height, units, dpi)
  width_px <- to_pixels(width, units, dpi)
  height_px <- to_pixels(height, units, dpi)
  if (ext == "pdf") {
    grDevices::pdf(f0, width = width_in, height = height_in, onefile = TRUE)
  } else if (ext == "png") {
    grDevices::png(f0, width = width_px, height = height_px, res = dpi)
  } else if (ext == "jpeg") {
    grDevices::jpeg(f0, width = width_px, height = height_px, res = dpi, quality = 100)
  } else if (ext == "tiff") {
    grDevices::tiff(f0, width = width_px, height = height_px, res = dpi, compression = "lzw")
  } else if (ext == "bmp") {
    grDevices::bmp(f0, width = width_px, height = height_px, res = dpi)
  } else {
    stop("Unsupported file extension: ", ext, call. = FALSE)
  }
  abs_path(f0)
}

.circos_bind_named_list <- function(data, track.col = NULL, track.order = NULL) {
  require_pkg("data.table")
  if (!(is.list(data) && !is.data.frame(data) && !data.table::is.data.table(data))) {
    return(list(data = data, track.col = track.col, from.named.list = FALSE))
  }
  ok_df <- vapply(data, function(x) is.data.frame(x) || data.table::is.data.table(x), logical(1))
  if (!all(ok_df)) {
    stop("If data is a list for circos(), every element must be a data.frame/data.table.", call. = FALSE)
  }
  nms <- names(data)
  if (is.null(nms) || !length(nms)) nms <- rep("", length(data))
  fill_idx <- which(is.na(nms) | !nzchar(nms))
  if (length(fill_idx)) {
    tr_order <- as.character(track.order)
    tr_order <- tr_order[!is.na(tr_order) & nzchar(tr_order)]
    if (length(tr_order) == length(data)) {
      nms[fill_idx] <- tr_order[fill_idx]
    } else {
      nms[fill_idx] <- as.character(fill_idx)
    }
  }
  nms <- make.unique(nms, sep = ".")
  track.col.use <- as.character(track.col)[1]
  if (is.null(track.col) || !length(track.col.use) || is.na(track.col.use) || !nzchar(track.col.use)) {
    track.col.use <- "..gcanvas.track.."
  }
  while (any(vapply(data, function(x) track.col.use %in% names(x), logical(1)))) {
    track.col.use <- paste0(track.col.use, ".")
  }
  parts <- vector("list", length(data))
  for (i in seq_along(data)) {
    di <- if (data.table::is.data.table(data[[i]])) data.table::copy(data[[i]]) else data.table::as.data.table(data[[i]])
    di[, (track.col.use) := as.character(nms[i])]
    parts[[i]] <- di
  }
  list(
    data = data.table::rbindlist(parts, use.names = TRUE, fill = TRUE),
    track.col = track.col.use,
    from.named.list = TRUE
  )
}

.circos_prepare_input <- function(data,
                                  snp.col,
                                  chrom.col,
                                  pos.col,
                                  p.col,
                                  y.col,
                                  track.col,
                                  track.order = NULL,
                                  lead.label.col = NULL,
                                  silent = FALSE) {
  require_pkg("data.table")
  list_input <- .circos_bind_named_list(data = data, track.col = track.col, track.order = track.order)
  data <- list_input$data
  if (!is.data.frame(data) && !data.table::is.data.table(data)) {
    stop("data must be a data.frame/data.table, or a list of data.frame/data.table objects.", call. = FALSE)
  }
  if (isTRUE(list_input$from.named.list)) {
    track.col <- list_input$track.col
  }
  if (missing(track.col) || is.null(track.col) || !length(track.col) ||
      is.na(as.character(track.col)[1]) || !nzchar(as.character(track.col)[1])) {
    stop("track.col is required for circos() when data is a data.frame/data.table.", call. = FALSE)
  }
  if (!is.null(y.col) && length(y.col) > 1L) stop("y.col must be a single column name.", call. = FALSE)
  if (is.null(y.col) && length(p.col) > 1L) stop("p.col must be a single column name.", call. = FALSE)

  snp_col_use <- .gcanvas_resolve_colname(names(data), snp.col, required = TRUE, arg_label = "snp.col")
  chrom_col_use <- .gcanvas_resolve_colname(names(data), chrom.col, required = TRUE, arg_label = "chrom.col")
  pos_col_use <- .gcanvas_resolve_colname(names(data), pos.col, required = TRUE, arg_label = "pos.col")
  track_col_use <- .gcanvas_resolve_colname(names(data), track.col, required = TRUE, arg_label = "track.col")
  y_col_use <- if (is.null(y.col)) NULL else .gcanvas_resolve_colname(names(data), y.col, required = TRUE, arg_label = "y.col")
  p_col_use <- if (is.null(y.col)) .gcanvas_resolve_colname(names(data), p.col, required = TRUE, arg_label = "p.col") else NULL
  lead_label_col_use <- NULL
  if (!is.null(lead.label.col) && length(lead.label.col) && !is.na(as.character(lead.label.col)[1]) &&
      nzchar(as.character(lead.label.col)[1])) {
    lead_label_col_use <- .gcanvas_resolve_colname(names(data), lead.label.col, required = TRUE, arg_label = "lead.label.col")
  }

  cols <- unique(c(snp_col_use, chrom_col_use, pos_col_use, track_col_use, y_col_use, p_col_use, lead_label_col_use))
  dt <- if (data.table::is.data.table(data)) data.table::copy(data[, ..cols]) else data.table::as.data.table(data[, cols, drop = FALSE])
  data.table::setnames(dt, snp_col_use, "snp")
  data.table::setnames(dt, chrom_col_use, "CHR")
  data.table::setnames(dt, pos_col_use, "POS")
  data.table::setnames(dt, track_col_use, "track")
  if (!is.null(y_col_use)) data.table::setnames(dt, y_col_use, "Y")
  if (!is.null(p_col_use)) data.table::setnames(dt, p_col_use, "P")
  if (!is.null(lead_label_col_use) && lead_label_col_use %in% names(dt)) data.table::setnames(dt, lead_label_col_use, "lead_label_src")

  dt[, `:=`(
    snp = as.character(snp),
    CHR = normalize.chrom(CHR),
    POS = suppressWarnings(as.numeric(POS)),
    track = as.character(track)
  )]
  dt[is.na(track) | !nzchar(track), track := NA_character_]
  if (is.null(y.col)) {
    dt[, P := .gcanvas_p_filter(P)]
    dt[, `:=`(value.raw = as.character(P), y.raw = -log10c(P))]
  } else {
    dt[, `:=`(Y = suppressWarnings(as.numeric(Y)), value.raw = as.character(Y), y.raw = as.numeric(Y))]
  }
  dt <- dt[!is.na(snp) & nzchar(snp) & !is.na(CHR) & nzchar(CHR) & is.finite(POS) & !is.na(track) & nzchar(track)]
  dt <- dt[is.finite(y.raw) & !is.na(y.raw)]
  if (!nrow(dt)) stop("No valid rows after parsing.", call. = FALSE)
  dt[, .row_id := .I]

  tracks_seen <- unique(dt$track)
  if (is.null(track.order) || !length(track.order)) {
    track.order.use <- tracks_seen
  } else {
    track.order.use <- as.character(track.order)
    track.order.use <- track.order.use[!is.na(track.order.use) & nzchar(track.order.use)]
    miss <- setdiff(track.order.use, tracks_seen)
    if (length(miss)) .gcanvas_warn_msg(sprintf("track.order values not found in data and ignored: %s", paste(miss, collapse = ", ")))
    track.order.use <- track.order.use[track.order.use %in% tracks_seen]
    track.order.use <- c(track.order.use, setdiff(tracks_seen, track.order.use))
  }

  .gcanvas_note(
    "gcanvas::circos",
    sprintf(
      "n_rows=%d | n_tracks=%d | n_chromosomes=%d%s",
      as_int(nrow(dt)),
      as_int(length(track.order.use)),
      as_int(data.table::uniqueN(dt$CHR)),
      if (isTRUE(list_input$from.named.list)) " | input=named.list" else ""
    ),
    silent = silent
  )

  list(
    dt = dt,
    track.order = track.order.use,
    input.from.named.list = isTRUE(list_input$from.named.list),
    columns = list(
      snp.col = snp_col_use,
      chrom.col = chrom_col_use,
      pos.col = pos_col_use,
      p.col = p_col_use,
      y.col = y_col_use,
      track.col = track_col_use,
      lead.label.col = lead_label_col_use
    )
  )
}

.circos_prepare_layout <- function(dt,
                                   build = 38L,
                                   chroms = "auto",
                                   chroms.drop = FALSE,
                                   chrom.gap = NULL,
                                   gap.degree = NULL,
                                   last.chrom.gap = 9.5,
                                   start.degree = 90,
                                   track.order,
                                   track.height = 1,
                                   track.gap = 0.15,
                                   inner.radius = 0.6,
                                   outward = TRUE) {
  require_pkg("data.table")
  chr_obs <- .gcanvas_sort_chr_unique(dt$CHR)
  chroms_auto <- is.null(chroms) || length(chroms) == 0L ||
    (length(chroms) == 1L && is.character(chroms) && tolower(trimws(as.character(chroms)[1])) == "auto")
  if (isTRUE(chroms_auto)) {
    chr_keep <- chr_obs
  } else {
    chr_req <- normalize.chrom(chroms)
    chr_req <- chr_req[!is.na(chr_req) & nzchar(chr_req)]
    chr_req <- unique(chr_req)
    if (isTRUE(chroms.drop)) {
      chr_keep <- intersect(chr_req, chr_obs)
    } else {
      chr_keep <- chr_req
    }
  }
  if (!length(chr_keep)) stop("No chromosomes remain after applying chroms/chroms.drop.", call. = FALSE)
  dt <- dt[CHR %in% chr_keep]
  if (!nrow(dt)) stop("No rows remain after chromosome filtering.", call. = FALSE)

  obs_bounds <- dt[, .(obs_end = max(POS, na.rm = TRUE)), by = CHR]
  ref_bounds <- .gcanvas_hg_chr_bounds(build = build)
  chr_map <- merge(data.table::data.table(CHR = chr_keep), ref_bounds, by = "CHR", all.x = TRUE, sort = FALSE)
  chr_map <- merge(chr_map, obs_bounds, by = "CHR", all.x = TRUE, sort = FALSE)
  chr_map[, chr_end := pmax(as.numeric(end), as.numeric(obs_end), 1, na.rm = TRUE)]
  chr_map[!is.finite(chr_end) | is.na(chr_end) | chr_end <= 0, chr_end := 1]

  gap_in <- if (!is.null(gap.degree) && length(gap.degree)) gap.degree else chrom.gap
  if (is.null(gap_in) || length(gap_in) == 0L) {
    gap_vec <- rep(1.5, nrow(chr_map))
    if (length(gap_vec)) {
      last.chrom.gap0 <- suppressWarnings(as.numeric(unlist(last.chrom.gap, use.names = FALSE))[1])
      if (is.finite(last.chrom.gap0) && !is.na(last.chrom.gap0) && last.chrom.gap0 >= 0) {
        gap_vec[length(gap_vec)] <- last.chrom.gap0
      } else {
        gap_vec[length(gap_vec)] <- 9.5
      }
    }
  } else {
    gap_vec <- suppressWarnings(as.numeric(gap_in))
    gap_vec <- gap_vec[is.finite(gap_vec) & !is.na(gap_vec) & gap_vec >= 0]
    if (!length(gap_vec)) gap_vec <- rep(1.5, nrow(chr_map))
    if (length(gap_vec) == 1L) gap_vec <- rep(gap_vec, nrow(chr_map))
    if (length(gap_vec) != nrow(chr_map)) stop("gap.degree/chrom.gap must be length 1 or equal to the number of chromosomes.", call. = FALSE)
    last.chrom.gap0 <- suppressWarnings(as.numeric(unlist(last.chrom.gap, use.names = FALSE))[1])
    if (length(gap_vec) && is.finite(last.chrom.gap0) && !is.na(last.chrom.gap0) && last.chrom.gap0 >= 0) {
      gap_vec[length(gap_vec)] <- last.chrom.gap0
    }
  }
  if (sum(gap_vec) >= 360) stop("Sum of chromosome gaps must be < 360 degrees.", call. = FALSE)

  start.degree <- suppressWarnings(as.numeric(start.degree))[1]
  if (!is.finite(start.degree) || is.na(start.degree)) start.degree <- 90
  span.total <- 360 - sum(gap_vec)
  chr_map[, width.degree := span.total * (chr_end / sum(chr_end))]
  chr_map[, gap.after := gap_vec]
  chr_map[, start.degree := NA_real_]
  chr_map[, end.degree := NA_real_]
  cur <- start.degree
  for (i in seq_len(nrow(chr_map))) {
    chr_map$start.degree[i] <- cur
    chr_map$end.degree[i] <- cur - chr_map$width.degree[i]
    cur <- chr_map$end.degree[i] - chr_map$gap.after[i]
  }
  chr_map[, center.degree := (start.degree + end.degree) / 2]

  dt <- merge(dt, chr_map[, .(CHR, chr_end, start.degree, end.degree, width.degree)], by = "CHR", all.x = TRUE, sort = FALSE)
  frac <- ifelse(dt$chr_end <= 1, 0, (dt$POS - 1) / dt$chr_end)
  frac[!is.finite(frac) | is.na(frac)] <- 0
  frac <- pmax(0, pmin(1, frac))
  dt[, angle.degree := start.degree - frac * width.degree]

  track.height.default <- suppressWarnings(as.numeric(unlist(track.height, use.names = FALSE))[1])
  if (!is.finite(track.height.default) || is.na(track.height.default) || track.height.default <= 0) track.height.default <- 1
  track.gap.default <- suppressWarnings(as.numeric(unlist(track.gap, use.names = FALSE))[1])
  if (!is.finite(track.gap.default) || is.na(track.gap.default) || track.gap.default < 0) track.gap.default <- 0.15
  inner.radius <- suppressWarnings(as.numeric(inner.radius))[1]
  if (!is.finite(inner.radius) || is.na(inner.radius) || inner.radius < 0.05) inner.radius <- 0.6

  track.plot.order <- if (isTRUE(outward)) as.character(track.order) else rev(as.character(track.order))
  track_height_map <- .circos_track_scalar_map(track.height, track.order, default = track.height.default, numeric = TRUE)
  track_gap_map <- .circos_track_scalar_map(track.gap, track.order, default = track.gap.default, numeric = TRUE)
  track_height_map[!is.finite(track_height_map) | is.na(track_height_map) | track_height_map <= 0] <- track.height.default
  track_gap_map[!is.finite(track_gap_map) | is.na(track_gap_map) | track_gap_map < 0] <- track.gap.default
  track_meta <- data.table::data.table(
    track = track.plot.order,
    plot_index = seq_along(track.plot.order),
    track_height = as.numeric(track_height_map[match(track.plot.order, names(track_height_map))]),
    track_gap = as.numeric(track_gap_map[match(track.plot.order, names(track_gap_map))])
  )
  cur_r <- inner.radius
  track_meta[, `:=`(r.inner = NA_real_, r.outer = NA_real_)]
  for (i in seq_len(nrow(track_meta))) {
    h0 <- track_meta$track_height[i]
    g0 <- track_meta$track_gap[i]
    track_meta$r.inner[i] <- cur_r
    track_meta$r.outer[i] <- cur_r + h0
    cur_r <- track_meta$r.outer[i] + g0
  }
  track_meta <- track_meta[match(track.order, track)]

  outer.radius <- max(track_meta$r.outer, na.rm = TRUE)
  chr.band.outer <- max(0.14, inner.radius - max(0.05, track.gap.default * 0.55))
  chr.band.inner <- max(0.04, chr.band.outer - min(0.14, max(0.08, track.height.default * 0.28)))
  if (chr.band.outer <= chr.band.inner) {
    chr.band.outer <- max(0.14, inner.radius * 0.8)
    chr.band.inner <- max(0.04, chr.band.outer - 0.1)
  }

  list(
    dt = dt,
    chr.map = chr_map,
    track.meta = track_meta,
    outer.radius = outer.radius,
    chr.band.inner = chr.band.inner,
    chr.band.outer = chr.band.outer
  )
}

.circos_normalize_lead_table <- function(lead,
                                         tracks,
                                         lead.color = "#E63946",
                                         lead.size = NULL,
                                         lead.stroke = 1,
                                         lead.flank = NULL,
                                         lead.flank.color = NULL,
                                         lead.label = FALSE,
                                         lead.label.color = "grey20") {
  require_pkg("data.table")
  empty <- data.table::data.table(
    lead_id = integer(),
    track = character(),
    snp = character(),
    CHR = character(),
    POS = numeric(),
    lead_color = character(),
    lead_size = numeric(),
    lead_stroke = numeric(),
    flank_bp = numeric(),
    flank_color = character(),
    lead_label = character(),
    lead_label_color = character()
  )
  if (is.null(lead) || length(lead) == 0L || (is.logical(lead) && length(lead) == 1L && !isTRUE(lead))) {
    return(empty)
  }
  if (is.character(lead) && length(lead) == 1L && !is.na(lead[1]) && nzchar(lead[1]) && file.exists(lead[1])) {
    lead <- data.table::fread(lead[1], data.table = TRUE, showProgress = FALSE)
  }

  parse_lead_dt <- function(x) {
    x <- data.table::as.data.table(x)
    nms <- names(x)
    low <- tolower(nms)
    pick <- function(cands) {
      hit <- match(cands, low, nomatch = 0L)
      hit <- hit[hit > 0L]
      if (!length(hit)) return(NA_character_)
      nms[hit[1]]
    }
    track_col <- pick(c("track", "trait", "phenotype"))
    snp_col <- pick(c("snp", "snpid", "rsid", "variant", "variant_id"))
    chr_col <- pick(c("chr", "chrom", "chromosome"))
    pos_col <- pick(c("pos", "bp", "position"))
    color_col <- pick(c("lead_color", "color"))
    size_col <- pick(c("lead_size", "size"))
    stroke_col <- pick(c("lead_stroke", "stroke"))
    flank_col <- pick(c("flank", "lead_flank"))
    flank_color_col <- pick(c("flank_color", "lead_flank_color", "region_color"))
    label_col <- pick(c("lead_label", "label", "text"))
    label_color_col <- pick(c("lead_label_color", "label_color", "text_color"))
    out <- data.table::data.table(
      track = if (!is.na(track_col)) as.character(x[[track_col]]) else NA_character_,
      snp = if (!is.na(snp_col)) as.character(x[[snp_col]]) else NA_character_,
      CHR = if (!is.na(chr_col)) normalize.chrom(x[[chr_col]]) else NA_character_,
      POS = if (!is.na(pos_col)) suppressWarnings(as.numeric(x[[pos_col]])) else NA_real_,
      lead_color = if (!is.na(color_col)) as.character(x[[color_col]]) else NA_character_,
      lead_size = if (!is.na(size_col)) suppressWarnings(as.numeric(x[[size_col]])) else NA_real_,
      lead_stroke = if (!is.na(stroke_col)) suppressWarnings(as.numeric(x[[stroke_col]])) else NA_real_,
      flank_bp = if (!is.na(flank_col)) vapply(x[[flank_col]], .gcanvas_parse_bp_span, numeric(1), arg_name = "lead.flank") else NA_real_,
      flank_color = if (!is.na(flank_color_col)) as.character(x[[flank_color_col]]) else NA_character_,
      lead_label = if (!is.na(label_col)) as.character(x[[label_col]]) else NA_character_,
      lead_label_color = if (!is.na(label_color_col)) as.character(x[[label_color_col]]) else NA_character_
    )
    out
  }

  if (is.data.frame(lead) || data.table::is.data.table(lead)) {
    out <- parse_lead_dt(lead)
  } else if (is.list(lead) && !is.data.frame(lead) && !data.table::is.data.table(lead) &&
             !is.null(names(lead)) && any(!is.na(names(lead)) & nzchar(names(lead)))) {
    parts <- vector("list", length(lead))
    for (i in seq_along(lead)) {
      gi <- as.character(names(lead)[i])
      li <- lead[[i]]
      if (is.null(li) || length(li) == 0L) next
      if (is.data.frame(li) || data.table::is.data.table(li)) {
        dt0 <- parse_lead_dt(li)
        dt0[is.na(track) | !nzchar(track), track := gi]
        parts[[i]] <- dt0
      } else {
        snps <- .gcanvas_as_snp_vector(li)
        if (!length(snps)) next
        parts[[i]] <- data.table::data.table(track = gi, snp = snps)
      }
    }
    parts <- Filter(Negate(is.null), parts)
    if (!length(parts)) return(empty)
    out <- data.table::rbindlist(parts, use.names = TRUE, fill = TRUE)
  } else {
    snps <- .gcanvas_as_snp_vector(lead)
    out <- data.table::data.table(track = NA_character_, snp = snps)
  }
  if (!nrow(out)) return(empty)
  need_cols <- c("track", "snp", "CHR", "POS", "lead_color", "lead_size", "lead_stroke", "flank_bp", "flank_color", "lead_label", "lead_label_color")
  col_defaults <- list(
    track = NA_character_,
    snp = NA_character_,
    CHR = NA_character_,
    POS = NA_real_,
    lead_color = NA_character_,
    lead_size = NA_real_,
    lead_stroke = NA_real_,
    flank_bp = NA_real_,
    flank_color = NA_character_,
    lead_label = NA_character_,
    lead_label_color = NA_character_
  )
  miss_cols <- setdiff(need_cols, names(out))
  if (length(miss_cols)) {
    for (nm in miss_cols) out[, (nm) := col_defaults[[nm]]]
  }

  out[, `:=`(
    track = as.character(track),
    snp = as.character(snp),
    CHR = normalize.chrom(CHR),
    POS = suppressWarnings(as.numeric(POS)),
    lead_size = suppressWarnings(as.numeric(lead_size)),
    lead_stroke = suppressWarnings(as.numeric(lead_stroke)),
    flank_bp = suppressWarnings(as.numeric(flank_bp))
  )]
  out[!nzchar(track), track := NA_character_]
  out[!nzchar(snp), snp := NA_character_]
  out[!is.finite(POS), POS := NA_real_]
  out[!is.finite(lead_size) | lead_size <= 0, lead_size := NA_real_]
  out[!is.finite(lead_stroke) | lead_stroke < 0, lead_stroke := NA_real_]
  out[!is.finite(flank_bp) | flank_bp <= 0, flank_bp := NA_real_]
  out <- out[(!is.na(snp) & nzchar(snp)) | (!is.na(CHR) & nzchar(CHR) & is.finite(POS))]
  if (!nrow(out)) return(empty)

  track_color_map <- .circos_track_scalar_map(lead.color, tracks, default = "#E63946")
  track_size_map <- .circos_track_scalar_map(lead.size, tracks, default = 1.05, numeric = TRUE)
  track_stroke_map <- .circos_track_scalar_map(lead.stroke, tracks, default = 1.4, numeric = TRUE)
  track_flank_map <- .circos_track_bp_map(lead.flank, tracks, default = NA_real_, arg_name = "lead.flank")
  track_flank_color_map <- .circos_track_scalar_map(lead.flank.color, tracks, default = NA_character_)
  track_label_color_map <- .circos_track_scalar_map(lead.label.color, tracks, default = "grey20")

  map_track_default <- function(track_vec, mp, fallback = NA) {
    out0 <- rep_len(fallback, length(track_vec))
    hit <- match(as.character(track_vec), names(mp))
    ok <- is.finite(hit) & !is.na(hit) & hit > 0L
    out0[ok] <- mp[hit[ok]]
    if (length(mp)) out0[!ok] <- rep_len(mp[1], sum(!ok))
    out0
  }

  out[is.na(lead_color) | !nzchar(lead_color), lead_color := as.character(map_track_default(track, track_color_map, "#E63946"))]
  out[is.na(lead_color) | !nzchar(lead_color), lead_color := "#E63946"]
  out[is.na(lead_size), lead_size := as.numeric(map_track_default(track, track_size_map, 1.05))]
  out[!is.finite(lead_size) | lead_size <= 0, lead_size := 1.05]
  out[is.na(lead_stroke), lead_stroke := as.numeric(map_track_default(track, track_stroke_map, 1.4))]
  out[!is.finite(lead_stroke) | lead_stroke < 0, lead_stroke := 1.4]
  out[is.na(flank_bp), flank_bp := as.numeric(map_track_default(track, track_flank_map, NA_real_))]
  out[is.na(flank_color) | !nzchar(flank_color), flank_color := as.character(map_track_default(track, track_flank_color_map, NA_character_))]
  out[is.na(flank_color) | !nzchar(flank_color), flank_color := lead_color]
  out[is.na(lead_label_color) | !nzchar(lead_label_color), lead_label_color := as.character(map_track_default(track, track_label_color_map, "grey20"))]
  out[is.na(lead_label_color) | !nzchar(lead_label_color), lead_label_color := "grey20"]
  if (isTRUE(lead.label)) {
    out[is.na(lead_label) | !nzchar(lead_label), lead_label := snp]
  }
  out[, lead_id := .I]
  data.table::setcolorder(out, c(
    "lead_id", "track", "snp", "CHR", "POS", "lead_color", "lead_size", "lead_stroke",
    "flank_bp", "flank_color", "lead_label", "lead_label_color"
  ))
  out[]
}

.circos_apply_leads <- function(dt,
                                lead_tbl,
                                base.point.color = NULL,
                                alpha = 1,
                                lead.flank.alpha = 1,
                                lead.label.col = NULL,
                                lead.label = FALSE) {
  require_pkg("data.table")
  if (!is.null(base.point.color) && length(base.point.color) == 1L &&
      is.character(base.point.color) && base.point.color %in% names(dt)) {
    dt[, plot.color := as.character(get(base.point.color))]
  } else if (!is.null(base.point.color) && length(base.point.color) == nrow(dt)) {
    dt[, plot.color := as.character(base.point.color)]
  } else {
    dt[, plot.color := "grey60"]
  }
  alpha_map <- .circos_track_scalar_map(alpha, unique(dt$track), default = 1, numeric = TRUE)
  alpha_map[!is.finite(alpha_map) | is.na(alpha_map)] <- 1
  alpha_map <- pmax(0, pmin(1, alpha_map))
  flank_alpha_map <- .circos_track_scalar_map(lead.flank.alpha, unique(dt$track), default = alpha_map, numeric = TRUE)
  flank_alpha_map[!is.finite(flank_alpha_map) | is.na(flank_alpha_map)] <- alpha_map[match(names(flank_alpha_map), names(alpha_map))]
  flank_alpha_map <- pmax(0, pmin(1, flank_alpha_map))
  dt[, `:=`(
    plot.alpha = as.numeric(alpha_map[match(track, names(alpha_map))]),
    lead.hit = FALSE,
    lead.flank.hit = FALSE
  )]
  dt[!is.finite(plot.alpha) | is.na(plot.alpha), plot.alpha := 1]
  dt[is.na(plot.color) | !nzchar(plot.color), plot.color := "grey60"]

  if (!nrow(lead_tbl)) {
    return(list(dt = dt, lead.tbl = lead_tbl, lead.draw = dt[0], lead.missing = 0L))
  }

  snp_map <- dt[, .(CHR_data = CHR[1], POS_data = POS[1], track_first = track[1]), by = snp]
  lead_tbl <- merge(lead_tbl, snp_map, by = "snp", all.x = TRUE, sort = FALSE)
  lead_tbl[is.na(CHR) & !is.na(CHR_data), CHR := CHR_data]
  lead_tbl[is.na(POS) & is.finite(POS_data), POS := POS_data]
  if (!("track_first" %in% names(lead_tbl))) lead_tbl[, track_first := NA_character_]
  if (is.null(lead.label.col) || !(lead.label.col %in% names(dt))) {
    lead_lab_map <- NULL
  } else {
    lead_lab_map <- dt[!is.na(snp) & nzchar(snp), .SD[1], by = .(snp, track), .SDcols = lead.label.col]
    data.table::setnames(lead_lab_map, lead.label.col, "lead_label_data")
  }
  if (!is.null(lead_lab_map)) {
    lead_tbl <- merge(lead_tbl, lead_lab_map, by = c("snp", "track"), all.x = TRUE, sort = FALSE)
    lead_tbl[(is.na(lead_label) | !nzchar(lead_label)) & !is.na(lead_label_data), lead_label := as.character(lead_label_data)]
    lead_tbl[, lead_label_data := NULL]
  }
  if (isTRUE(lead.label)) lead_tbl[(is.na(lead_label) | !nzchar(lead_label)) & !is.na(snp), lead_label := snp]

  lead_draw_parts <- list()
  flank_hits <- list()
  lead_missing <- 0L
  for (i in seq_len(nrow(lead_tbl))) {
    li <- lead_tbl[i]
    track_mask <- if (!is.na(li$track[[1]]) && nzchar(li$track[[1]])) dt$track == li$track[[1]] else rep(TRUE, nrow(dt))
    exact_idx <- integer()
    if (!is.na(li$snp[[1]]) && nzchar(li$snp[[1]])) {
      exact_idx <- which(track_mask & dt$snp == li$snp[[1]])
    }
    if (!length(exact_idx) && !is.na(li$CHR[[1]]) && nzchar(li$CHR[[1]]) && is.finite(li$POS[[1]])) {
      exact_idx <- which(track_mask & dt$CHR == li$CHR[[1]] & dt$POS == li$POS[[1]])
    }
    if (!length(exact_idx)) {
      lead_missing <- lead_missing + 1L
    } else {
      dt[exact_idx, `:=`(lead.hit = TRUE, plot.color = as.character(li$lead_color[[1]]))]
      lead_draw_parts[[length(lead_draw_parts) + 1L]] <- data.table::copy(dt[exact_idx])[, `:=`(
        lead_id = as_int(li$lead_id[[1]]),
        lead_color = as.character(li$lead_color[[1]]),
        lead_size = as.numeric(li$lead_size[[1]]),
        lead_stroke = as.numeric(li$lead_stroke[[1]]),
        lead_label = as.character(li$lead_label[[1]]),
        lead_label_color = as.character(li$lead_label_color[[1]])
      )]
    }
    if (is.finite(li$flank_bp[[1]]) && !is.na(li$CHR[[1]]) && nzchar(li$CHR[[1]]) && is.finite(li$POS[[1]])) {
      idx <- which(track_mask & dt$CHR == li$CHR[[1]] & dt$POS >= (li$POS[[1]] - li$flank_bp[[1]]) & dt$POS <= (li$POS[[1]] + li$flank_bp[[1]]))
      if (length(idx)) {
        alpha_now <- flank_alpha_map[match(as.character(li$track[[1]]), names(flank_alpha_map))]
        if (!length(alpha_now) || !is.finite(alpha_now) || is.na(alpha_now)) alpha_now <- 1
        flank_hits[[length(flank_hits) + 1L]] <- data.table::data.table(
          .row_id = dt$.row_id[idx],
          lead_id = as_int(li$lead_id[[1]]),
          track = as.character(dt$track[idx]),
          dist = abs(dt$POS[idx] - li$POS[[1]]),
          flank_color = as.character(li$flank_color[[1]]),
          flank_alpha = as.numeric(alpha_now)
        )
      }
    }
  }
  if (length(flank_hits)) {
    flank_dt <- data.table::rbindlist(flank_hits, use.names = TRUE, fill = TRUE)
    flank_dt[, flank_alpha := as.numeric(flank_alpha_map[match(track, names(flank_alpha_map))])]
    flank_dt[!is.finite(flank_alpha) | is.na(flank_alpha), flank_alpha := 1]
    data.table::setorderv(flank_dt, c(".row_id", "dist", "lead_id"), c(1L, 1L, 1L), na.last = TRUE)
    flank_dt <- flank_dt[!duplicated(.row_id)]
    dt[flank_dt$.row_id, `:=`(
      lead.flank.hit = TRUE,
      plot.color = flank_dt$flank_color,
      plot.alpha = flank_dt$flank_alpha
    )]
  }
  lead_draw <- if (length(lead_draw_parts)) {
    x <- data.table::rbindlist(lead_draw_parts, use.names = TRUE, fill = TRUE)
    data.table::setorderv(x, c(".row_id", "lead_id"), c(1L, 1L), na.last = TRUE)
    x[!duplicated(.row_id)]
  } else {
    dt[0]
  }
  list(dt = dt, lead.tbl = lead_tbl, lead.draw = lead_draw, lead.missing = lead_missing)
}

