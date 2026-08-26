# Circos-style genome-wide plot.
#
# NOTE: this function is experimental and is being actively developed.
# Interfaces and defaults may change in a future release.

#' Circos plot (experimental)
#'
#' @description
#' **Experimental / work in progress.** Interfaces and defaults may change.
#'
#' Draws a circos-style overview of genome-wide association data: one ring
#' per chromosome, with points placed by position and colored / sized by the
#' chosen statistic.
#'
#' @param data A `data.frame`/`data.table` of summary statistics.
#' @param snp.col,chrom.col,pos.col,p.col Column names in `data`.
#' @param y.col Optional alternative y-column override.
#' @param track.col,track.order Optional column / order for stacking multiple
#'   data tracks in the circle.
#' @param build Genome build (`37` or `38`).
#' @param chroms,chroms.drop Which chromosomes to include / drop.
#' @param chrom.gap,start.degree,gap.degree,last.chrom.gap Per-chromosome
#'   layout angles.
#' @param alpha,point.size,chrom.color,point.shape Base point styling.
#' @param track.height,track.gap,inner.radius Track geometry.
#' @param inner.band,inner.band.bin.size,inner.band.metric,inner.band.color,inner.band.zero.color
#'   Inner heat-band controls.
#' @param outward Logical. Direction the data fills (outward vs inward).
#' @param threshold,threshold.color,threshold.type,threshold.linewidth
#'   Significance-line controls.
#' @param lead,lead.color,lead.size,lead.stroke Lead-variant highlighting.
#' @param lead.flank,lead.flank.color,lead.flank.alpha,lead.line Flank shading
#'   and connectors around lead variants.
#' @param lead.shadow,lead.shadow.flank,lead.shadow.color,lead.shadow.alpha
#'   Drop-shadow effect around the lead variant.
#' @param lead.label,lead.label.col,lead.label.color,lead.label.size,lead.label.angle
#'   Lead-variant text labels.
#' @param seed Integer seed for the label layout RNG.
#' @param y.max,y.breaks,y.rescale.at,y.rescale.ratio,y.rescale.breaks,y.rescale.line,y.rescale.line.color,y.rescale.line.type,y.rescale.line.linewidth
#'   Y-axis range and break-and-rescale controls.
#' @param chrom.label,chrom.label.col,chrom.label.size Chromosome-label styling.
#' @param outer.band,outer.band.bin.size,outer.band.metric,outer.band.color,outer.band.zero.color,outer.band.height,outer.band.gap
#'   Outer heat-band controls.
#' @param axis,axis.col,axis.linewidth,axis.text.size,axis.text.offset Axis styling.
#' @param title Optional plot title.
#' @param file,width,height,units,dpi Output file controls (`NULL` to display only).
#' @param return.data Logical. Return the prepared layout data alongside the plot.
#' @param silent Logical. Suppress progress notes.
#'
#' @return A plot object (currently a base-graphics render; subject to change).
#' @export
circos <- function(data,
                   snp.col = "SNP",
                   chrom.col = "CHR",
                   pos.col = "POS",
                   p.col = "P",
                   y.col = NULL,
                   track.col = NULL,
                   track.order = NULL,
                   build = 38L,
                   chroms = "auto",
                   chroms.drop = FALSE,
                   chrom.gap = NULL,
                   start.degree = 90,
                   gap.degree = NULL,
                   last.chrom.gap = 10,
                   alpha = 1,
                   point.size = 0.7,
                   chrom.color = c("grey60", "grey80"),
                   point.shape = 16,
                   track.height = 1,
                   track.gap = 0.15,
                   inner.radius = 0.6,
                   inner.band = FALSE,
                   inner.band.bin.size = "1Mb",
                   inner.band.metric = "chromosome",
                   inner.band.color = NULL,
                   inner.band.zero.color = "grey90",
                   outward = TRUE,
                   threshold = 5e-8,
                   threshold.color = "#E63946",
                   threshold.type = "auto",
                   threshold.linewidth = 0.7,
                   lead = NULL,
                   lead.color = "#E63946",
                   lead.size = 1.05,
                   lead.stroke = 1.4,
                   lead.flank = NULL,
                   lead.flank.color = NULL,
                   lead.flank.alpha = 1,
                   lead.line = FALSE,
                   lead.shadow = FALSE,
                   lead.shadow.flank = NULL,
                   lead.shadow.color = NULL,
                   lead.shadow.alpha = 0.28,
                   lead.label = FALSE,
                   lead.label.col = NULL,
                   lead.label.color = "grey20",
                   lead.label.size = 3.5,
                   lead.label.angle = 0,
                   seed = 23L,
                   y.max = NULL,
                   y.breaks = NULL,
                   y.rescale.at = NULL,
                   y.rescale.ratio = 0.05,
                   y.rescale.breaks = NULL,
                   y.rescale.line = TRUE,
                   y.rescale.line.color = "grey50",
                   y.rescale.line.type = "solid",
                   y.rescale.line.linewidth = 0.5,
                   chrom.label = TRUE,
                   chrom.label.col = "grey20",
                   chrom.label.size = 3,
                   outer.band = FALSE,
                   outer.band.bin.size = "1Mb",
                   outer.band.metric = "auto",
                   outer.band.color = "auto",
                   outer.band.zero.color = "grey90",
                   outer.band.height = 0.12,
                   outer.band.gap = 0.04,
                   axis = TRUE,
                   axis.col = "grey20",
                   axis.linewidth = 0.6,
                   axis.text.size = 0.88,
                   axis.text.offset = 0.05,
                   title = NULL,
                   file = NULL,
                   width = 25,
                   height = 25,
                   units = "cm",
                   dpi = 300,
                   return.data = FALSE,
                   silent = FALSE) {
  require_pkg("data.table")
  point.size.in <- point.size
  point.shape.in <- point.shape
  alpha.in <- alpha
  lead.flank.alpha.in <- lead.flank.alpha
  lead.shadow.alpha.in <- lead.shadow.alpha
  lead.label.size.in <- lead.label.size
  lead.label.angle.in <- lead.label.angle
  axis.col.in <- axis.col
  axis.linewidth.in <- axis.linewidth
  axis.text.size.in <- axis.text.size
  axis.text.offset.in <- axis.text.offset
  point.size.default <- suppressWarnings(as.numeric(unlist(point.size.in, use.names = FALSE))[1])
  if (!is.finite(point.size.default) || is.na(point.size.default) || point.size.default <= 0) point.size.default <- 0.9
  point.shape.default <- suppressWarnings(as.integer(unlist(point.shape.in, use.names = FALSE))[1])
  if (!is.finite(point.shape.default) || is.na(point.shape.default)) point.shape.default <- 16L
  alpha.default <- suppressWarnings(as.numeric(unlist(alpha.in, use.names = FALSE))[1])
  if (!is.finite(alpha.default) || is.na(alpha.default)) alpha.default <- 1
  alpha.default <- max(0, min(1, alpha.default))
  lead.flank.alpha.default <- suppressWarnings(as.numeric(unlist(lead.flank.alpha.in, use.names = FALSE))[1])
  if (!is.finite(lead.flank.alpha.default) || is.na(lead.flank.alpha.default)) lead.flank.alpha.default <- alpha.default
  lead.flank.alpha.default <- max(0, min(1, lead.flank.alpha.default))
  lead.shadow.alpha.default <- suppressWarnings(as.numeric(unlist(lead.shadow.alpha.in, use.names = FALSE))[1])
  if (!is.finite(lead.shadow.alpha.default) || is.na(lead.shadow.alpha.default)) lead.shadow.alpha.default <- 0.28
  lead.shadow.alpha.default <- max(0, min(1, lead.shadow.alpha.default))
  lead.label.size.default <- suppressWarnings(as.numeric(unlist(lead.label.size.in, use.names = FALSE))[1])
  if (!is.finite(lead.label.size.default) || is.na(lead.label.size.default) || lead.label.size.default <= 0) lead.label.size.default <- 3.5
  lead.label.angle.default <- suppressWarnings(as.numeric(unlist(lead.label.angle.in, use.names = FALSE))[1])
  if (!is.finite(lead.label.angle.default) || is.na(lead.label.angle.default)) lead.label.angle.default <- 0
  chrom.label.size <- suppressWarnings(as.numeric(chrom.label.size))[1]
  if (!is.finite(chrom.label.size) || is.na(chrom.label.size) || chrom.label.size <= 0) chrom.label.size <- 3
  axis.linewidth.default <- suppressWarnings(as.numeric(unlist(axis.linewidth.in, use.names = FALSE))[1])
  if (!is.finite(axis.linewidth.default) || is.na(axis.linewidth.default) || axis.linewidth.default <= 0) axis.linewidth.default <- 0.5
  axis.col.default <- as.character(unlist(axis.col.in, use.names = FALSE))[1]
  if (is.na(axis.col.default) || !nzchar(axis.col.default)) axis.col.default <- "grey20"
  axis.text.size.default <- suppressWarnings(as.numeric(unlist(axis.text.size.in, use.names = FALSE))[1])
  if (!is.finite(axis.text.size.default) || is.na(axis.text.size.default) || axis.text.size.default <= 0) axis.text.size.default <- 0.88
  axis.text.offset.default <- suppressWarnings(as.numeric(unlist(axis.text.offset.in, use.names = FALSE))[1])
  if (!is.finite(axis.text.offset.default) || is.na(axis.text.offset.default) || axis.text.offset.default < 0) axis.text.offset.default <- 0.05
  outer.band.height <- suppressWarnings(as.numeric(outer.band.height))[1]
  if (!is.finite(outer.band.height) || is.na(outer.band.height) || outer.band.height < 0) outer.band.height <- 0.12
  outer.band.gap <- suppressWarnings(as.numeric(outer.band.gap))[1]
  if (!is.finite(outer.band.gap) || is.na(outer.band.gap) || outer.band.gap < 0) outer.band.gap <- 0.04
  inner.band.enabled <- !(is.null(inner.band) || length(inner.band) == 0L || (is.logical(inner.band) && length(inner.band) == 1L && !isTRUE(inner.band)))
  outer.band.enabled <- !(is.null(outer.band) || length(outer.band) == 0L || (is.logical(outer.band) && length(outer.band) == 1L && !isTRUE(outer.band)))
  seed <- as_int(seed)
  if (!is.finite(seed) || is.na(seed)) seed <- 23L
  set.seed(seed)

  prep <- .circos_prepare_input(
    data = data,
    snp.col = snp.col,
    chrom.col = chrom.col,
    pos.col = pos.col,
    p.col = p.col,
    y.col = y.col,
    track.col = track.col,
    track.order = track.order,
    lead.label.col = lead.label.col,
    silent = silent
  )
  dt <- prep$dt
  track.order <- prep$track.order
  col_meta <- prep$columns

  layout <- .circos_prepare_layout(
    dt = dt,
    build = build,
    chroms = chroms,
    chroms.drop = chroms.drop,
    chrom.gap = chrom.gap,
    gap.degree = gap.degree,
    last.chrom.gap = last.chrom.gap,
    start.degree = start.degree,
    track.order = track.order,
    track.height = track.height,
    track.gap = track.gap,
    inner.radius = inner.radius,
    outward = outward
  )
  dt <- layout$dt
  chr_map <- layout$chr.map
  track_meta <- layout$track.meta
  outer.radius <- layout$outer.radius
  chr.band.inner <- layout$chr.band.inner
  chr.band.outer <- layout$chr.band.outer
  axis.angle <- suppressWarnings(as.numeric(chr_map$start.degree[1]))
  if (!is.finite(axis.angle) || is.na(axis.angle)) axis.angle <- 90
  axis.rad <- axis.angle * pi / 180
  tick.length.major <- 0.06
  tick.length.minor <- 0.035
  zero.label.nudge <- 0.028

  chrom_palette_key <- NA_character_
  chrom.cols <- as.character(chrom.color)
  if (length(chrom.cols) == 1L && !is.na(chrom.cols) && nzchar(chrom.cols)) {
    key <- tolower(trimws(chrom.cols[1]))
    pal_list <- .gcanvas_chrom_palette_list()
    if (key %in% names(pal_list)) {
      chrom.cols <- as.character(pal_list[[key]])
      chrom_palette_key <- key
    }
  }
  chrom.cols <- chrom.cols[!is.na(chrom.cols) & nzchar(chrom.cols)]
  if (!length(chrom.cols)) chrom.cols <- c("grey60", "grey80")
  if (identical(chrom_palette_key, "random")) {
    .gcanvas_note("gcanvas::circos", sprintf("chrom.color=random -> %s", paste(chrom.cols, collapse = ", ")), silent = silent)
  }
  chr_map[, chr_plot_color := rep(chrom.cols, length.out = .N)]
  dt[, chr_plot_color := chr_map$chr_plot_color[match(CHR, chr_map$CHR)]]
  dt[is.na(chr_plot_color) | !nzchar(chr_plot_color), chr_plot_color := rep(chrom.cols, length.out = .N)]

  point.size.map <- .circos_track_scalar_map(point.size.in, track.order, default = point.size.default, numeric = TRUE)
  point.size.map[!is.finite(point.size.map) | is.na(point.size.map) | point.size.map <= 0] <- point.size.default
  point.shape.map <- .circos_track_scalar_map(point.shape.in, track.order, default = point.shape.default, numeric = TRUE)
  point.shape.map[!is.finite(point.shape.map) | is.na(point.shape.map)] <- point.shape.default
  point.shape.map <- as.integer(round(point.shape.map))
  alpha.map <- .circos_track_scalar_map(alpha.in, track.order, default = alpha.default, numeric = TRUE)
  alpha.map[!is.finite(alpha.map) | is.na(alpha.map)] <- alpha.default
  alpha.map <- pmax(0, pmin(1, alpha.map))
  lead.flank.alpha.map <- .circos_track_scalar_map(lead.flank.alpha.in, track.order, default = lead.flank.alpha.default, numeric = TRUE)
  lead.flank.alpha.map[!is.finite(lead.flank.alpha.map) | is.na(lead.flank.alpha.map)] <- lead.flank.alpha.default
  lead.flank.alpha.map <- pmax(0, pmin(1, lead.flank.alpha.map))
  lead.shadow.alpha.map <- .circos_track_scalar_map(lead.shadow.alpha.in, track.order, default = lead.shadow.alpha.default, numeric = TRUE)
  lead.shadow.alpha.map[!is.finite(lead.shadow.alpha.map) | is.na(lead.shadow.alpha.map)] <- lead.shadow.alpha.default
  lead.shadow.alpha.map <- pmax(0, pmin(1, lead.shadow.alpha.map))
  lead.label.size.map <- .circos_track_scalar_map(lead.label.size.in, track.order, default = lead.label.size.default, numeric = TRUE)
  lead.label.size.map[!is.finite(lead.label.size.map) | is.na(lead.label.size.map) | lead.label.size.map <= 0] <- lead.label.size.default
  lead.label.angle.map <- .circos_track_scalar_map(lead.label.angle.in, track.order, default = lead.label.angle.default, numeric = TRUE)
  lead.label.angle.map[!is.finite(lead.label.angle.map) | is.na(lead.label.angle.map)] <- lead.label.angle.default
  axis.col.map <- .circos_track_scalar_map(axis.col.in, track.order, default = axis.col.default)
  axis.col.map[is.na(axis.col.map) | !nzchar(axis.col.map)] <- axis.col.default
  axis.linewidth.map <- .circos_track_scalar_map(axis.linewidth.in, track.order, default = axis.linewidth.default, numeric = TRUE)
  axis.linewidth.map[!is.finite(axis.linewidth.map) | is.na(axis.linewidth.map) | axis.linewidth.map <= 0] <- axis.linewidth.default
  axis.text.size.map <- .circos_track_scalar_map(axis.text.size.in, track.order, default = axis.text.size.default, numeric = TRUE)
  axis.text.size.map[!is.finite(axis.text.size.map) | is.na(axis.text.size.map) | axis.text.size.map <= 0] <- axis.text.size.default
  axis.text.offset.map <- .circos_track_scalar_map(axis.text.offset.in, track.order, default = axis.text.offset.default, numeric = TRUE)
  axis.text.offset.map[!is.finite(axis.text.offset.map) | is.na(axis.text.offset.map) | axis.text.offset.map < 0] <- axis.text.offset.default

  lead_tbl <- .circos_normalize_lead_table(
    lead = lead,
    tracks = track.order,
    lead.color = lead.color,
    lead.size = lead.size,
    lead.stroke = lead.stroke,
    lead.flank = lead.flank,
    lead.flank.color = lead.flank.color,
    lead.label = lead.label,
    lead.label.color = lead.label.color
  )
  lead_applied <- .circos_apply_leads(
    dt = dt,
    lead_tbl = lead_tbl,
    base.point.color = "chr_plot_color",
    alpha = alpha.map,
    lead.flank.alpha = lead.flank.alpha.map,
    lead.label.col = "lead_label_src",
    lead.label = lead.label
  )
  dt <- lead_applied$dt
  lead_tbl <- lead_applied$lead.tbl
  lead_draw <- lead_applied$lead.draw
  lead_missing <- lead_applied$lead.missing
  shadow_tbl <- .circos_prepare_shadow_table(
    lead.shadow = lead.shadow,
    lead_tbl = lead_tbl,
    tracks = track.order,
    lead.color = lead.color,
    lead.flank = lead.flank,
    lead.flank.color = lead.flank.color,
    lead.shadow.flank = lead.shadow.flank,
    lead.shadow.color = lead.shadow.color
  )
  shadow_regions <- .circos_build_shadow_regions(
    shadow_tbl = shadow_tbl,
    chr_map = chr_map,
    track_meta = track_meta,
    tracks = track.order,
    shadow.alpha.map = lead.shadow.alpha.map
  )
  inner.band.metric0 <- tolower(trimws(as.character(inner.band.metric)[1]))
  if (is.na(inner.band.metric0) || !nzchar(inner.band.metric0)) inner.band.metric0 <- "chromosome"
  inner.band.use.default <- isTRUE(inner.band.enabled) &&
    isTRUE(is.logical(inner.band) && length(inner.band) == 1L && isTRUE(inner.band)) &&
    identical(inner.band.metric0, "chromosome") &&
    is.null(inner.band.color)
  inner_band_tbl <- if (isTRUE(inner.band.enabled) && !isTRUE(inner.band.use.default)) {
    .circos_prepare_outer_band(
      dt = dt,
      lead_tbl = lead_tbl,
      chr_map = chr_map,
      outer.band = inner.band,
      outer.band.bin.size = inner.band.bin.size,
      outer.band.metric = if (identical(inner.band.metric0, "chromosome")) "variant" else inner.band.metric,
      outer.band.color = if (is.null(inner.band.color)) "auto" else inner.band.color,
      outer.band.zero.color = inner.band.zero.color,
      band.name = "inner.band"
    )
  } else {
    data.table::data.table()
  }
  outer_band_tbl <- .circos_prepare_outer_band(
    dt = dt,
    lead_tbl = lead_tbl,
    chr_map = chr_map,
    outer.band = outer.band,
    outer.band.bin.size = outer.band.bin.size,
    outer.band.metric = outer.band.metric,
    outer.band.color = outer.band.color,
    outer.band.zero.color = outer.band.zero.color,
    band.name = "outer.band"
  )
  outer.band.inner <- outer.radius + outer.band.gap
  outer.band.outer <- outer.band.inner + outer.band.height

  threshold_map <- .circos_track_list_map(threshold, track.order, default = NULL)
  threshold_color_map <- .circos_track_list_map(threshold.color, track.order, default = "grey20")
  threshold_type_map <- .circos_track_list_map(threshold.type, track.order, default = NULL)
  threshold_lwd_map <- .circos_track_list_map(threshold.linewidth, track.order, default = 0.7)
  y_max_map <- .circos_track_scalar_map(y.max, track.order, default = NA_real_, numeric = TRUE)
  y_breaks_map <- .circos_track_list_map(y.breaks, track.order, default = NULL)
  y_rescale_at_map <- .circos_track_scalar_map(y.rescale.at, track.order, default = NA_real_, numeric = TRUE)
  y_rescale_ratio_map <- .circos_track_scalar_map(y.rescale.ratio, track.order, default = 0.25, numeric = TRUE)
  y_rescale_breaks_map <- .circos_track_list_map(y.rescale.breaks, track.order, default = NULL)
  y_rescale_line_color_map <- .circos_track_scalar_map(y.rescale.line.color, track.order, default = "grey20")
  y_rescale_line_type_map <- .circos_track_scalar_map(y.rescale.line.type, track.order, default = "solid")
  y_rescale_line_lwd_map <- .circos_track_scalar_map(y.rescale.line.linewidth, track.order, default = 0.5, numeric = TRUE)
  y_rescale_line_lwd_map[!is.finite(y_rescale_line_lwd_map) | is.na(y_rescale_line_lwd_map) | y_rescale_line_lwd_map <= 0] <- 0.5

  track_state <- vector("list", length(track.order))
  names(track_state) <- track.order
  dt[, y.draw := NA_real_]
  dt[, y.clip := NA_real_]
  dt[, radius := NA_real_]
  if (nrow(lead_draw)) {
    lead_draw[, `:=`(y.draw = NA_real_, y.clip = NA_real_, radius = NA_real_)]
  }
  axis_text_draw <- vector("list", 0L)
  major.grid.col <- "grey90"
  ring.boundary.col <- "grey20"
  .circos_radius_from_y <- function(y_draw, tm, y_scale_top) {
    yy <- suppressWarnings(as.numeric(y_draw))
    top <- suppressWarnings(as.numeric(y_scale_top))[1]
    if (!is.finite(top) || is.na(top) || top <= 0) top <- 1
    yy[!is.finite(yy) | is.na(yy)] <- 0
    yy <- pmax(0, pmin(yy, top))
    tm$r.inner[[1]] + (yy / top) * (tm$r.outer[[1]] - tm$r.inner[[1]])
  }

  for (tr in track.order) {
    idx <- which(dt$track == tr)
    if (!length(idx)) next
    y_at <- y_rescale_at_map[tr]
    y_ratio <- y_rescale_ratio_map[tr]
    y_map_fun <- function(v) .circos_apply_y_rescale(v, at = y_at, ratio = y_ratio)
    thr_lines <- .gcanvas_threshold_lines(
      threshold = threshold_map[[tr]],
      y.col = y.col,
      threshold.color = threshold_color_map[[tr]],
      threshold.type = threshold_type_map[[tr]],
      threshold.linewidth = threshold_lwd_map[[tr]],
      map_fun = y_map_fun,
      default_color = "grey20"
    )
    y_raw <- dt$y.raw[idx]
    y_draw <- y_map_fun(y_raw)
    y_max_raw_in <- y_max_map[tr]
    y_max_user <- is.finite(y_max_raw_in) && !is.na(y_max_raw_in) && y_max_raw_in > 0
    if (!isTRUE(y_max_user)) {
      y_max_raw_in <- max(c(y_raw, thr_lines$y_raw), na.rm = TRUE)
      if (!is.finite(y_max_raw_in) || is.na(y_max_raw_in) || y_max_raw_in <= 0) y_max_raw_in <- max(y_raw, na.rm = TRUE)
      if (!is.finite(y_max_raw_in) || is.na(y_max_raw_in) || y_max_raw_in <= 0) y_max_raw_in <- 1
    }
    y_max_raw_i <- .circos_nice_y_max(y_max_raw_in, user_specified = y_max_user)
    y_max_draw_i <- y_map_fun(y_max_raw_i)
    if (!is.finite(y_max_draw_i) || is.na(y_max_draw_i) || y_max_draw_i <= 0) y_max_draw_i <- max(y_draw, na.rm = TRUE)
    if (!is.finite(y_max_draw_i) || is.na(y_max_draw_i) || y_max_draw_i <= 0) y_max_draw_i <- 1
    y_scale_top_i <- y_max_draw_i
    if (!is.finite(y_scale_top_i) || is.na(y_scale_top_i) || y_scale_top_i <= 0) y_scale_top_i <- y_max_raw_i
    br <- .circos_compute_y_breaks(
      y_raw = y_raw,
      y_max_raw = y_max_raw_i,
      y_breaks = y_breaks_map[[tr]],
      y_rescale.at = y_at,
      y_rescale.ratio = y_ratio,
      y_rescale.breaks = y_rescale_breaks_map[[tr]]
    )
    tm <- track_meta[track == tr]
    y_clip <- pmin(y_draw, y_scale_top_i)
    radius_now <- .circos_radius_from_y(y_clip, tm, y_scale_top_i)
    dt[idx, `:=`(y.draw = y_draw, y.clip = y_clip, radius = radius_now)]
    if (nrow(lead_draw)) {
      ld_idx <- which(lead_draw$track == tr)
      if (length(ld_idx)) {
        y_raw_ld <- lead_draw$y.raw[ld_idx]
        y_draw_ld <- y_map_fun(y_raw_ld)
        y_clip_ld <- pmin(y_draw_ld, y_scale_top_i)
        radius_ld <- .circos_radius_from_y(y_clip_ld, tm, y_scale_top_i)
        lead_draw[ld_idx, `:=`(y.draw = y_draw_ld, y.clip = y_clip_ld, radius = radius_ld)]
      }
    }
    if (nrow(thr_lines)) {
      thr_lines[, `:=`(
        y_max_raw = y_max_raw_i,
        y_max_draw = y_max_draw_i,
        y_scale_top = y_scale_top_i,
        radius = .circos_radius_from_y(pmin(y, y_scale_top_i), tm, y_scale_top_i),
        track = tr
      )]
    }
    track_state[[tr]] <- list(
      track = tr,
      y_max_raw = y_max_raw_i,
      y_max_draw = y_max_draw_i,
      y_scale_top = y_scale_top_i,
      y_rescale.at = y_at,
      y_rescale.ratio = y_ratio,
      y_breaks = br,
      threshold_lines = thr_lines
    )
  }

  out_path <- .circos_open_device(file = file, width = width, height = height, units = units, dpi = dpi)
  if (!is.null(out_path)) on.exit(grDevices::dev.off(), add = TRUE)
  op <- graphics::par(no.readonly = TRUE)
  on.exit(try(graphics::par(op), silent = TRUE), add = TRUE)
  gap_pad <- if ("track_gap" %in% names(track_meta) && any(is.finite(track_meta$track_gap), na.rm = TRUE)) {
    max(track_meta$track_gap[is.finite(track_meta$track_gap)], na.rm = TRUE)
  } else {
    0.15
  }
  outer.extra <- if (isTRUE(outer.band.enabled) && nrow(outer_band_tbl)) outer.band.gap + outer.band.height else 0
  lim <- outer.radius + outer.extra + max(0.35, gap_pad * 2.2, 0.22 + 0.04 * length(track.order))
  graphics::par(mar = c(1, 1, if (is.null(title)) 1 else 2.5, 1), xpd = NA)
  graphics::plot.new()
  graphics::plot.window(xlim = c(-lim, lim), ylim = c(-lim, lim), asp = 1)

  if (isTRUE(inner.band.use.default)) {
    for (i in seq_len(nrow(chr_map))) {
      fill_col <- chr_map$chr_plot_color[i]
      poly <- .circos_ring_polygon(chr_map$start.degree[i], chr_map$end.degree[i], chr.band.inner, chr.band.outer, n = 90L)
      graphics::polygon(poly$x, poly$y, col = fill_col, border = "white", lwd = 0.6)
    }
  } else if (isTRUE(inner.band.enabled) && nrow(inner_band_tbl)) {
    for (i in seq_len(nrow(inner_band_tbl))) {
      poly <- .circos_ring_polygon(
        inner_band_tbl$start.degree[i],
        inner_band_tbl$end.degree[i],
        chr.band.inner,
        chr.band.outer,
        n = 20L
      )
      graphics::polygon(poly$x, poly$y, col = inner_band_tbl$fill[i], border = NA)
    }
    for (i in seq_len(nrow(chr_map))) {
      xy0 <- .circos_arc_xy(chr_map$start.degree[i], chr_map$end.degree[i], chr.band.inner, n = 90L)
      xy1 <- .circos_arc_xy(chr_map$start.degree[i], chr_map$end.degree[i], chr.band.outer, n = 90L)
      graphics::lines(xy0$x, xy0$y, col = "grey70", lwd = 0.4)
      graphics::lines(xy1$x, xy1$y, col = "grey70", lwd = 0.4)
    }
  }

  if (nrow(shadow_regions)) {
    for (i in seq_len(nrow(shadow_regions))) {
      poly <- .circos_ring_polygon(
        shadow_regions$start.degree[i],
        shadow_regions$end.degree[i],
        shadow_regions$r.inner[i],
        shadow_regions$r.outer[i],
        n = 60L
      )
      fill_now <- .circos_add_alpha(shadow_regions$fill[i], shadow_regions$alpha[i])
      graphics::polygon(poly$x, poly$y, col = fill_now, border = NA)
    }
  }

  if (isTRUE(outer.band.enabled) && nrow(outer_band_tbl)) {
    for (i in seq_len(nrow(outer_band_tbl))) {
      poly <- .circos_ring_polygon(
        outer_band_tbl$start.degree[i],
        outer_band_tbl$end.degree[i],
        outer.band.inner,
        outer.band.outer,
        n = 20L
      )
      graphics::polygon(poly$x, poly$y, col = outer_band_tbl$fill[i], border = NA)
    }
    for (i in seq_len(nrow(chr_map))) {
      xy0 <- .circos_arc_xy(chr_map$start.degree[i], chr_map$end.degree[i], outer.band.inner, n = 90L)
      xy1 <- .circos_arc_xy(chr_map$start.degree[i], chr_map$end.degree[i], outer.band.outer, n = 90L)
      graphics::lines(xy0$x, xy0$y, col = "grey70", lwd = 0.4)
      graphics::lines(xy1$x, xy1$y, col = "grey70", lwd = 0.4)
    }
  }

  if (nrow(lead_draw)) {
    lead_draw <- lead_draw[is.finite(radius) & is.finite(angle.degree)]
    if (isTRUE(lead.line) && nrow(lead_draw)) {
      lead_line_draw <- unique(lead_draw[, .(CHR, POS, angle.degree)])
      if (nrow(lead_line_draw)) {
        line.r.inner <- min(track_meta$r.inner, na.rm = TRUE)
        line.r.outer <- max(track_meta$r.outer, na.rm = TRUE)
        if (is.finite(line.r.inner) && is.finite(line.r.outer) && line.r.outer > line.r.inner) {
          ang_line <- lead_line_draw$angle.degree * pi / 180
          x0 <- line.r.inner * cos(ang_line)
          y0 <- line.r.inner * sin(ang_line)
          x1 <- line.r.outer * cos(ang_line)
          y1 <- line.r.outer * sin(ang_line)
          graphics::segments(x0, y0, x1, y1, col = "grey50", lty = 2, lwd = 0.8)
        }
      }
    }
  }

  if (isTRUE(axis)) {
    for (tr in track.order) {
      tm <- track_meta[track == tr]
      if (!nrow(tm)) next
      br <- track_state[[tr]]$y_breaks
      for (i in seq_len(nrow(chr_map))) {
        xy0 <- .circos_arc_xy(chr_map$start.degree[i], chr_map$end.degree[i], tm$r.inner[[1]], n = 90L)
        xy1 <- .circos_arc_xy(chr_map$start.degree[i], chr_map$end.degree[i], tm$r.outer[[1]], n = 90L)
        graphics::lines(xy0$x, xy0$y, col = ring.boundary.col, lwd = 0.6)
        graphics::lines(xy1$x, xy1$y, col = ring.boundary.col, lwd = 0.6)
      }
      if (is.data.frame(br) && nrow(br)) {
        for (j in seq_len(nrow(br))) {
          rj <- .circos_radius_from_y(br$y_draw[j], tm, track_state[[tr]]$y_scale_top)
          if (j > 1L && j < nrow(br)) {
            for (i in seq_len(nrow(chr_map))) {
              arc <- .circos_arc_xy(chr_map$start.degree[i], chr_map$end.degree[i], rj, n = 90L)
              graphics::lines(arc$x, arc$y, col = major.grid.col, lwd = 0.45)
            }
          }
        }
      }
    }
  }

  for (tr in track.order) {
    dsub <- dt[track == tr]
    if (!nrow(dsub)) next
    point.size.now <- point.size.map[tr]
    if (!length(point.size.now) || !is.finite(point.size.now) || is.na(point.size.now) || point.size.now <= 0) point.size.now <- point.size.default
    point.shape.now <- point.shape.map[tr]
    if (!length(point.shape.now) || !is.finite(point.shape.now) || is.na(point.shape.now)) point.shape.now <- point.shape.default
    draw_point_layer <- function(xdt) {
      if (!nrow(xdt)) return(invisible(NULL))
      ang <- xdt$angle.degree * pi / 180
      x <- xdt$radius * cos(ang)
      y <- xdt$radius * sin(ang)
      pcols <- .circos_add_alpha(xdt$plot.color, xdt$plot.alpha)
      if (point.shape.now %in% 21:25) {
        graphics::points(x, y, pch = point.shape.now, cex = point.size.now, bg = pcols, col = pcols)
      } else {
        graphics::points(x, y, pch = point.shape.now, cex = point.size.now, col = pcols)
      }
      invisible(NULL)
    }
    draw_point_layer(dsub[!(lead.flank.hit %in% TRUE) & !(lead.hit %in% TRUE)])
    draw_point_layer(dsub[(lead.flank.hit %in% TRUE) & !(lead.hit %in% TRUE)])
    thr <- track_state[[tr]]$threshold_lines
    if (nrow(thr)) {
      for (k in seq_len(nrow(thr))) {
        for (i in seq_len(nrow(chr_map))) {
          arc <- .circos_arc_xy(chr_map$start.degree[i], chr_map$end.degree[i], thr$radius[k], n = 90L)
          graphics::lines(arc$x, arc$y, col = thr$color[k], lwd = thr$linewidth[k], lty = .circos_lty(thr$type[k]))
        }
      }
    }
      if (isTRUE(y.rescale.line)) {
        y_at <- track_state[[tr]]$y_rescale.at
        if (is.finite(y_at) && !is.na(y_at) && y_at > 0 && y_at < track_state[[tr]]$y_max_raw) {
        r_rs <- .circos_radius_from_y(
          .circos_apply_y_rescale(y_at, at = y_at, ratio = track_state[[tr]]$y_rescale.ratio),
          track_meta[track == tr],
          track_state[[tr]]$y_scale_top
        )
          for (i in seq_len(nrow(chr_map))) {
            arc <- .circos_arc_xy(chr_map$start.degree[i], chr_map$end.degree[i], r_rs, n = 90L)
            graphics::lines(
              arc$x, arc$y,
            col = y_rescale_line_color_map[tr],
            lwd = y_rescale_line_lwd_map[tr],
            lty = .circos_lty(y_rescale_line_type_map[tr])
            )
          }
      }
    }
  }

  if (isTRUE(axis)) {
    for (tr in track.order) {
      tm <- track_meta[track == tr]
      if (!nrow(tm)) next
      axis.col.now <- axis.col.map[tr]
      if (!length(axis.col.now) || is.na(axis.col.now) || !nzchar(axis.col.now)) axis.col.now <- axis.col.default
      axis.linewidth.now <- axis.linewidth.map[tr]
      if (!length(axis.linewidth.now) || !is.finite(axis.linewidth.now) || is.na(axis.linewidth.now) || axis.linewidth.now <= 0) axis.linewidth.now <- axis.linewidth.default
      axis.text.size.now <- axis.text.size.map[tr]
      if (!length(axis.text.size.now) || !is.finite(axis.text.size.now) || is.na(axis.text.size.now) || axis.text.size.now <= 0) axis.text.size.now <- axis.text.size.default
      axis.text.offset.now <- axis.text.offset.map[tr]
      if (!length(axis.text.offset.now) || !is.finite(axis.text.offset.now) || is.na(axis.text.offset.now) || axis.text.offset.now < 0) axis.text.offset.now <- axis.text.offset.default
      graphics::segments(
        tm$r.inner[[1]] * cos(axis.rad),
        tm$r.inner[[1]] * sin(axis.rad),
        tm$r.outer[[1]] * cos(axis.rad),
        tm$r.outer[[1]] * sin(axis.rad),
        col = axis.col.now,
        lwd = axis.linewidth.now
      )
      br <- track_state[[tr]]$y_breaks
      if (is.data.frame(br) && nrow(br)) {
        minor_draw <- numeric()
        if (nrow(br) >= 2L) {
          minor_draw <- br$y_draw[-nrow(br)] + diff(br$y_draw) / 2
          minor_draw <- minor_draw[is.finite(minor_draw) & !is.na(minor_draw)]
        }
        if (length(minor_draw)) {
          for (rm in minor_draw) {
            rj <- .circos_radius_from_y(rm, tm, track_state[[tr]]$y_scale_top)
            x0 <- rj * cos(axis.rad)
            y0 <- rj * sin(axis.rad)
            x1 <- x0 - tick.length.minor
            y1 <- y0
            graphics::segments(x0, y0, x1, y1, col = axis.col.now, lwd = max(0.4, axis.linewidth.now * 0.9))
          }
        }
        for (j in seq_len(nrow(br))) {
          rj <- .circos_radius_from_y(br$y_draw[j], tm, track_state[[tr]]$y_scale_top)
          x0 <- rj * cos(axis.rad)
          y0 <- rj * sin(axis.rad)
          x1 <- x0 - tick.length.major
          y1 <- y0
          graphics::segments(x0, y0, x1, y1, col = axis.col.now, lwd = axis.linewidth.now)
          y_lab <- if (j == 1L) y1 + zero.label.nudge else y1
          axis_text_draw[[length(axis_text_draw) + 1L]] <- list(
            x = x1 - axis.text.offset.now,
            y = y_lab,
            label = br$label[j],
            col = axis.col.now,
            cex = axis.text.size.now
          )
        }
      }
    }
  }

  if (nrow(lead_draw)) {
    if (nrow(lead_draw)) {
      ang <- lead_draw$angle.degree * pi / 180
      x <- lead_draw$radius * cos(ang)
      y <- lead_draw$radius * sin(ang)
      border_col <- ifelse(is.finite(lead_draw$lead_stroke) & lead_draw$lead_stroke > 0, "grey20", lead_draw$lead_color)
      graphics::points(
        x, y,
        pch = 21,
        cex = lead_draw$lead_size,
        bg = lead_draw$lead_color,
        col = border_col,
        lwd = pmax(0.5, lead_draw$lead_stroke)
      )
      lab_dt <- lead_draw[!is.na(lead_label) & nzchar(lead_label)]
      if (nrow(lab_dt)) {
        ang_lab <- lab_dt$angle.degree * pi / 180
        r_lab <- lab_dt$radius + 0.08
        x_lab <- r_lab * cos(ang_lab)
        y_lab <- r_lab * sin(ang_lab)
        for (i in seq_len(nrow(lab_dt))) {
          cex_now <- lead.label.size.map[as.character(lab_dt$track[i])]
          if (!length(cex_now) || !is.finite(cex_now) || is.na(cex_now) || cex_now <= 0) cex_now <- lead.label.size.default
          angle_now <- lead.label.angle.map[as.character(lab_dt$track[i])]
          if (!length(angle_now) || !is.finite(angle_now) || is.na(angle_now)) angle_now <- lead.label.angle.default
          srt_now <- lab_dt$angle.degree[i] + angle_now - 90
          srt_now <- ifelse(srt_now < -90, srt_now + 180, ifelse(srt_now > 90, srt_now - 180, srt_now))
          graphics::text(
            x_lab[i], y_lab[i],
            labels = lab_dt$lead_label[i],
            cex = cex_now / 3.5,
            col = lab_dt$lead_label_color[i],
            srt = srt_now
          )
        }
      }
    }
  }

  if (length(axis_text_draw)) {
    for (lab in axis_text_draw) {
      graphics::text(
        lab$x,
        lab$y,
        labels = lab$label,
        cex = lab$cex,
        adj = c(1, 0.5),
        col = lab$col
      )
    }
  }

  if (isTRUE(chrom.label)) {
    lab_r <- outer.radius + outer.extra + 0.22
    chrom.label.col <- as.character(chrom.label.col)[1]
    if (is.na(chrom.label.col) || !nzchar(chrom.label.col)) chrom.label.col <- "grey20"
    for (i in seq_len(nrow(chr_map))) {
      ang0 <- chr_map$center.degree[i]
      rad0 <- ang0 * pi / 180
      x <- lab_r * cos(rad0)
      y <- lab_r * sin(rad0)
      srt <- ang0 - 90
      srt <- ((srt + 180) %% 360) - 180
      graphics::text(
        x, y,
        labels = chr_map$CHR[i],
        cex = chrom.label.size / 3,
        col = chrom.label.col,
        srt = srt,
        adj = c(0.5, 0.5),
        xpd = TRUE
      )
    }
  }

  if (!is.null(title)) graphics::title(main = title)

  if (!is.null(out_path) && nzchar(as.character(out_path)[1])) {
    .gcanvas_note("gcanvas::circos", paste0("saved: ", abs_path(out_path)), silent = silent)
  }

  meta <- list(
    type = "circos",
    engine = "base",
    build = as_int(build),
    columns = col_meta,
    n_rows = as_int(nrow(dt)),
    n_tracks = as_int(length(track.order)),
    n_chromosomes = as_int(nrow(chr_map)),
    track.order = track.order,
    chrom.colors = chrom.cols,
    inner.band = list(
      enabled = isTRUE(inner.band.enabled),
      metric = if (!isTRUE(inner.band.enabled)) {
        "none"
      } else if (isTRUE(inner.band.use.default)) {
        "chromosome"
      } else if (is.list(inner.band) && !is.data.frame(inner.band) && !data.table::is.data.table(inner.band) &&
                 !(is.logical(inner.band) && length(inner.band) == 1L)) {
        "custom"
      } else {
        inner.band.metric0
      },
      bin.size = inner.band.bin.size,
      n_bins = as_int(nrow(inner_band_tbl))
    ),
    chromosome.bounds = chr_map[, .(chr = CHR, end = chr_end, start.degree, end.degree, width.degree, chr_plot_color)],
    outer.band = list(
      enabled = isTRUE(outer.band.enabled),
      metric = if (isTRUE(outer.band.enabled)) {
        if (is.list(outer.band) && !is.data.frame(outer.band) && !data.table::is.data.table(outer.band) &&
            !(is.logical(outer.band) && length(outer.band) == 1L)) {
          "custom"
        } else {
          x <- tolower(trimws(as.character(outer.band.metric)[1]))
          if (is.na(x) || !nzchar(x) || identical(x, "auto")) if (nrow(lead_tbl)) "lead" else "variant" else x
        }
      } else {
        "none"
      },
      bin.size = outer.band.bin.size,
      n_bins = as_int(nrow(outer_band_tbl))
    ),
    highlights = list(
      n_lead = as_int(nrow(lead_tbl)),
      n_lead_missing = as_int(lead_missing),
      n_exact = as_int(sum(dt$lead.hit %in% TRUE, na.rm = TRUE)),
      n_flank = as_int(sum(dt$lead.flank.hit %in% TRUE, na.rm = TRUE))
    ),
    output = list(
      file = out_path %||% NULL,
      width = suppressWarnings(as.numeric(width))[1],
      height = suppressWarnings(as.numeric(height))[1],
      units = as.character(units)[1],
      dpi = suppressWarnings(as.numeric(dpi))[1]
    )
  )
  if (isTRUE(return.data)) {
    meta$data <- dt
    meta$lead <- lead_tbl
    meta$track.meta <- track_meta
    meta$track.state <- track_state
  }
  .gcanvas_note("gcanvas::circos", "Done", silent = silent)
  invisible(meta)
}

.miami_layer_col <- function(p, col) {
  out <- vector("list", length(p$layers))
  k <- 0L
  for (i in seq_along(p$layers)) {
    d <- p$layers[[i]]$data
    if (!is.data.frame(d) || !nrow(d) || !(col %in% names(d))) next
    k <- k + 1L
    out[[k]] <- d[[col]]
  }
  if (k == 0L) return(NULL)
  unlist(out[seq_len(k)], use.names = FALSE)
}

.miami_chr_summary <- function(meta, p) {
  chr_now <- character()
  cb <- meta$chromosome_bounds
  if (is.data.frame(cb) && nrow(cb) && ("chr" %in% names(cb))) chr_now <- cb$chr
  if (!length(chr_now)) chr_now <- .miami_layer_col(p, "CHR")
  chr_now <- .gcanvas_sort_chr_unique(chr_now)
  if (!length(chr_now)) return("n=0")
  auto_chr <- as.character(1:22)
  chr_extra <- c("X", "Y", "MT")
  if (all(auto_chr %in% chr_now)) {
    extra <- chr_extra[chr_extra %in% chr_now]
    if (length(extra)) return(paste0("1-22,", paste(extra, collapse = ",")))
    return("1-22")
  }
  if (length(chr_now) <= 12L) return(paste(chr_now, collapse = ","))
  sprintf("n=%d", as_int(length(chr_now)))
}

.miami_n_variants <- function(meta, p) {
  n0 <- as_int(meta$n_variants)[1]
  if (is.finite(n0) && !is.na(n0) && n0 > 0L) return(as_int(n0))
  snp <- as.character(.miami_layer_col(p, "snp"))
  snp <- snp[!is.na(snp) & nzchar(snp)]
  if (length(snp)) return(as_int(length(unique(snp))))
  n1 <- 0L
  for (i in seq_along(p$layers)) {
    d <- p$layers[[i]]$data
    if (is.data.frame(d) && nrow(d)) n1 <- max(n1, nrow(d))
  }
  as_int(n1)
}

.miami_minp <- function(p) {
  pv <- as.character(.miami_layer_col(p, "P"))
  if (!length(pv)) return("NA")
  # Find the minimum p-value in log10 space (underflow-safe via log10c) and
  # format the original string, matching manhattan(); converting to numeric
  # first would underflow extreme-small p-values to 0/NA and cap min-P near
  # the denormal floor (~1e-324).
  lp <- suppressWarnings(log10c(pv))
  ok <- is.finite(lp) & !is.na(lp) & lp <= 0
  if (!any(ok)) return("NA")
  i <- which(ok)[which.min(lp[ok])]
  .gcanvas_format_minp(pv[i], digits = 3L, cutoff = 1e-3)
}

.miami_maxy <- function(p) {
  yv <- suppressWarnings(as.numeric(.miami_layer_col(p, "yval")))
  yv <- yv[is.finite(yv) & !is.na(yv)]
  if (!length(yv)) {
    yv <- suppressWarnings(as.numeric(.miami_layer_col(p, "Y")))
    yv <- yv[is.finite(yv) & !is.na(yv)]
  }
  if (!length(yv)) return(NA_real_)
  max(yv, na.rm = TRUE)
}

.miami_summary_line <- function(tag, p, meta) {
  nvar <- .miami_n_variants(meta, p)
  chr_sum <- .miami_chr_summary(meta, p)
  ycol <- if (is.list(meta$columns)) meta$columns$y.col else NULL
  use_p <- is.null(ycol) || length(ycol) == 0L || is.na(ycol)[1]
  if (isTRUE(use_p)) {
    return(sprintf("%s: n_variants=%d | chromosomes=%s | minP=%s", tag, as_int(nvar), chr_sum, .miami_minp(p)))
  }
  maxy <- .miami_maxy(p)
  maxy_msg <- if (is.finite(maxy) && !is.na(maxy)) sprintf("%.4g", maxy) else "NA"
  sprintf("%s: n_variants=%d | chromosomes=%s | maxY=%s", tag, as_int(nvar), chr_sum, maxy_msg)
}

.miami_remove_scales <- function(p, aesthetics) {
  if (!inherits(p, "ggplot")) return(p)
  p <- unserialize(serialize(p, NULL))
  if (!length(p$scales$scales)) return(p)
  keep <- vapply(
    p$scales$scales,
    function(s) {
      aes <- s$aesthetics %||% character()
      !any(aes %in% aesthetics)
    },
    logical(1)
  )
  p$scales$scales <- p$scales$scales[keep]
  p
}

.miami_get_axis_info <- function(p, axis = c("x", "y"), built = NULL) {
  axis <- match.arg(axis)
  b <- built %||% suppressWarnings(ggplot2::ggplot_build(p))
  pp <- b$layout$panel_params[[1]]
  sc <- pp[[axis]]
  if (is.null(sc)) return(list(breaks = NULL, labels = NULL, minor_breaks = NULL, range = NULL))
  get_safe <- function(fn, default = NULL) tryCatch(fn(), error = function(e) default)
  breaks <- get_safe(function() sc$get_breaks(), NULL)
  if (is.null(breaks)) breaks <- sc$breaks %||% NULL
  minor_breaks <- get_safe(function() sc$get_minor_breaks(), NULL)
  if (is.null(minor_breaks)) minor_breaks <- sc$minor_breaks %||% NULL
  labels <- NULL
  if (!is.null(breaks) && is.function(sc$get_labels)) labels <- get_safe(function() sc$get_labels(breaks), NULL)
  if (is.null(labels)) labels <- sc$labels %||% NULL
  rng <- pp[[paste0(axis, ".range")]] %||% (sc$range$range %||% NULL)
  rng <- if (!is.null(rng)) suppressWarnings(as.numeric(rng)) else NULL
  list(breaks = breaks, labels = labels, minor_breaks = minor_breaks, range = rng)
}

.miami_point_ymax <- function(p, built = NULL) {
  b <- built %||% suppressWarnings(ggplot2::ggplot_build(p))
  vals <- numeric()
  for (i in seq_along(b$data)) {
    if (i > length(b$plot$layers)) next
    geom_name <- class(b$plot$layers[[i]]$geom)[1]
    if (!grepl("GeomPoint", geom_name, fixed = TRUE)) next
    d <- b$data[[i]]
    if ("y" %in% names(d)) {
      v <- suppressWarnings(as.numeric(d$y))
      v <- v[is.finite(v) & !is.na(v)]
      if (length(v)) vals <- c(vals, v)
    }
    if (all(c("ymin", "ymax") %in% names(d))) {
      v <- c(suppressWarnings(as.numeric(d$ymin)), suppressWarnings(as.numeric(d$ymax)))
      v <- v[is.finite(v) & !is.na(v)]
      if (length(v)) vals <- c(vals, v)
    }
  }
  if (!length(vals)) return(NA_real_)
  max(vals, na.rm = TRUE)
}

.miami_axis_from_meta <- function(meta) {
  cb <- meta$chromosome_bounds
  if (is.null(cb) || !is.data.frame(cb) || !all(c("chr", "start", "end") %in% names(cb))) return(NULL)
  dt <- data.table::as.data.table(cb)
  dt[, chr := normalize.chrom(chr)]
  dt[, start := suppressWarnings(as.numeric(start))]
  dt[, end := suppressWarnings(as.numeric(end))]
  dt <- dt[!is.na(chr) & nzchar(chr) & is.finite(start) & is.finite(end)]
  if (!nrow(dt)) return(NULL)
  sw <- which(dt$end < dt$start)
  if (length(sw)) {
    tmp <- dt$start[sw]
    dt$start[sw] <- dt$end[sw]
    dt$end[sw] <- tmp
  }
  dt[, width := pmax(1, end - start + 1)]
  dt[, chr_order := rank.chrom(chr)]
  data.table::setorderv(dt, c("chr_order", "chr"), c(1L, 1L), na.last = TRUE)
  dt[, offset := data.table::shift(cumsum(width), fill = 0)]
  dt[, tick := offset + width / 2]
  dt[]
}

.miami_axis_union <- function(a, b) {
  if (is.null(a) || !nrow(a)) return(NULL)
  if (is.null(b) || !nrow(b)) return(NULL)
  ab <- data.table::rbindlist(list(
    a[, .(chr = as.character(chr), start = as.numeric(start), end = as.numeric(end))],
    b[, .(chr = as.character(chr), start = as.numeric(start), end = as.numeric(end))]
  ), use.names = TRUE, fill = TRUE)
  ab[, chr := normalize.chrom(chr)]
  ab <- ab[!is.na(chr) & nzchar(chr) & is.finite(start) & is.finite(end)]
  if (!nrow(ab)) return(NULL)
  out <- ab[, .(start = min(start, na.rm = TRUE), end = max(end, na.rm = TRUE)), by = chr]
  sw <- which(out$end < out$start)
  if (length(sw)) {
    tmp <- out$start[sw]
    out$start[sw] <- out$end[sw]
    out$end[sw] <- tmp
  }
  out[, width := pmax(1, end - start + 1)]
  out[, chr_order := rank.chrom(chr)]
  data.table::setorderv(out, c("chr_order", "chr"), c(1L, 1L), na.last = TRUE)
  out[, offset := data.table::shift(cumsum(width), fill = 0)]
  out[, tick := offset + width / 2]
  out[]
}

.miami_is_drag_tick_layer <- function(layer) {
  inherits(layer$geom, "GeomSegment") &&
    is.data.frame(layer$data) &&
    setequal(names(layer$data), c("CHR", "x_tick"))
}

.miami_is_mask_rect_layer <- function(layer) {
  if (!inherits(layer$geom, "GeomRect")) return(FALSE)
  if (!is.data.frame(layer$data) || !all(c("xmin", "xmax", "ymin", "ymax") %in% names(layer$data))) return(FALSE)
  d <- layer$data
  if (!(is.infinite(d$xmin[1]) && is.infinite(d$xmax[1]))) return(FALSE)
  fill_val <- d$fill[1] %||% layer$params$fill %||% NA_character_
  fill_val <- tolower(as.character(fill_val)[1])
  fill_val %in% c("white", "#ffffff")
}

.miami_observed_chr <- function(p) {
  vals <- vector("list", length(p$layers))
  k <- 0L
  for (i in seq_along(p$layers)) {
    lyr <- p$layers[[i]]
    if (.miami_is_drag_tick_layer(lyr) || .miami_is_mask_rect_layer(lyr)) next
    d <- lyr$data
    if (!is.data.frame(d) || !nrow(d) || !("CHR" %in% names(d))) next
    k <- k + 1L
    vals[[k]] <- d$CHR
  }
  if (k == 0L) return(character())
  chr <- normalize.chrom(unlist(vals[seq_len(k)], use.names = FALSE))
  .gcanvas_sort_chr_unique(chr)
}

.miami_map_x <- function(chr, pos, axis_map) {
  chr0 <- normalize.chrom(chr)
  pos0 <- suppressWarnings(as.numeric(pos))
  key <- axis_map[, .(chr, start, offset)]
  ridx <- match(chr0, key$chr)
  out <- rep(NA_real_, length(pos0))
  ok <- !is.na(ridx) & is.finite(pos0)
  if (any(ok)) out[ok] <- pos0[ok] - key$start[ridx[ok]] + 1 + key$offset[ridx[ok]]
  out
}

.miami_theme_text_pt <- function(p, key = c("axis.text.x", "axis.title.x"), default = 14) {
  key <- match.arg(key)
  th <- p$theme %||% ggplot2::theme_get()
  elt <- switch(
    key,
    "axis.text.x" = (th$axis.text.x %||% th$axis.text),
    "axis.title.x" = (th$axis.title.x %||% th$axis.title)
  )
  sz <- suppressWarnings(as.numeric(elt$size))[1]
  if (!is.finite(sz) || is.na(sz) || sz <= 0) sz <- default
  sz
}

.miami_theme_color <- function(p, key = c("axis.text.x", "axis.line.x"), default = "grey20") {
  key <- match.arg(key)
  th <- p$theme %||% ggplot2::theme_get()
  elt <- switch(
    key,
    "axis.text.x" = (th$axis.text.x %||% th$axis.text),
    "axis.line.x" = (th$axis.line.x %||% th$axis.line)
  )
  col <- elt$colour %||% elt$color %||% default
  col <- as.character(col)[1]
  if (is.na(col) || !nzchar(col)) col <- default
  col
}

.miami_theme_linewidth <- function(p, key = c("axis.line.x", "panel.border"), default = 0.35) {
  key <- match.arg(key)
  th <- p$theme %||% ggplot2::theme_get()
  elt <- switch(
    key,
    "axis.line.x" = (th$axis.line.x %||% th$axis.line),
    "panel.border" = th$panel.border
  )
  lw <- suppressWarnings(as.numeric(elt$linewidth %||% elt$size))[1]
  if (!is.finite(lw) || is.na(lw) || lw <= 0) lw <- default
  lw
}

.miami_append_extra_cols <- function(built_df, layer, plot_data) {
  src <- NULL
  if (is.data.frame(layer$data)) src <- layer$data else if (is.data.frame(plot_data)) src <- plot_data
  if (is.null(src) || nrow(src) != nrow(built_df)) return(built_df)
  extra <- setdiff(names(src), names(built_df))
  if (!length(extra)) return(built_df)
  src_df <- as.data.frame(src)
  built_df[extra] <- src_df[, extra, drop = FALSE]
  built_df
}

.miami_transform_y <- function(df, mode = c("top", "bottom"), gap0, y_max0, panel_height0) {
  mode <- match.arg(mode)
  half <- gap0 / 2
  y_max0 <- suppressWarnings(as.numeric(y_max0))[1]
  if (!is.finite(y_max0) || is.na(y_max0) || y_max0 <= 0) y_max0 <- 1
  panel_height0 <- suppressWarnings(as.numeric(panel_height0))[1]
  if (!is.finite(panel_height0) || is.na(panel_height0) || panel_height0 <= 0) panel_height0 <- y_max0
  flip_val <- function(v) {
    if (!is.numeric(v)) return(v)
    if (identical(mode, "top")) return((v / y_max0) * panel_height0 + half)
    -((v / y_max0) * panel_height0) - half
  }
  if ("y" %in% names(df)) df$y <- flip_val(df$y)
  if ("yend" %in% names(df)) df$yend <- flip_val(df$yend)
  if ("yintercept" %in% names(df)) df$yintercept <- flip_val(df$yintercept)
  if (all(c("ymin", "ymax") %in% names(df)) && is.numeric(df$ymin) && is.numeric(df$ymax)) {
    new_min <- flip_val(df$ymin)
    new_max <- flip_val(df$ymax)
    df$ymin <- pmin(new_min, new_max)
    df$ymax <- pmax(new_min, new_max)
  } else {
    if ("ymin" %in% names(df)) df$ymin <- flip_val(df$ymin)
    if ("ymax" %in% names(df)) df$ymax <- flip_val(df$ymax)
  }
  if (identical(mode, "bottom") && "angle" %in% names(df)) {
    ang <- suppressWarnings(as.numeric(as.character(df$angle)))
    ok <- is.finite(ang) & !is.na(ang)
    if (any(ok)) {
      ang[ok] <- -1 * ang[ok]
      df$angle <- ang
    }
  }
  df
}

.miami_transform_x_using_meta <- function(df, axis_map) {
  if (is.null(axis_map) || !nrow(axis_map)) return(df)
  if (all(c("CHR", "POS", "x") %in% names(df))) {
    x_new <- .miami_map_x(df$CHR, df$POS, axis_map)
    ok <- is.finite(x_new) & !is.na(x_new)
    if (any(ok)) df$x[ok] <- x_new[ok]
  }
  if (all(c("CHR_data", "POS_data") %in% names(df))) {
    x_data_old <- if ("x_data" %in% names(df)) suppressWarnings(as.numeric(df$x_data)) else rep(NA_real_, nrow(df))
    x_data_new <- .miami_map_x(df$CHR_data, df$POS_data, axis_map)
    if ("x_data" %in% names(df)) {
      ok <- is.finite(x_data_new) & !is.na(x_data_new)
      if (any(ok)) df$x_data[ok] <- x_data_new[ok]
    }
    if ("x_label" %in% names(df) && "x_data" %in% names(df)) {
      dx <- suppressWarnings(as.numeric(df$x_label) - x_data_old)
      ok <- is.finite(x_data_new) & !is.na(x_data_new) & is.finite(dx)
      if (any(ok)) df$x_label[ok] <- x_data_new[ok] + dx[ok]
    }
  }
  df
}

.miami_make_identity_mapping <- function(df) {
  candidates <- c(
    "x", "y", "xend", "yend",
    "xmin", "xmax", "ymin", "ymax",
    "xintercept", "yintercept",
    "label", "group",
    "colour", "fill", "alpha",
    "size", "linewidth",
    "linetype", "shape", "stroke",
    "angle", "hjust", "vjust",
    "fontface", "family"
  )
  present <- intersect(candidates, names(df))
  if (!length(present)) return(ggplot2::aes())
  aes_args <- setNames(vector("list", length(present)), present)
  for (nm in present) aes_args[[nm]] <- as.name(nm)
  do.call(ggplot2::aes, aes_args)
}

.miami_extract_layers_built <- function(p, mode = c("top", "bottom"), axis_map, gap0,
                                                y_max0, panel_height0, direction.reverse = FALSE,
                                                drop.internal0 = TRUE, built = NULL) {
  mode <- match.arg(mode)
  b <- built %||% suppressWarnings(ggplot2::ggplot_build(p))
  layers_out <- list()
  dfs_out <- list()
  for (i in seq_along(b$data)) {
    orig_layer <- b$plot$layers[[i]]
    if (isTRUE(drop.internal0)) {
      if (.miami_is_drag_tick_layer(orig_layer)) next
      if (.miami_is_mask_rect_layer(orig_layer)) next
    }
    df <- b$data[[i]]
    df <- .miami_append_extra_cols(df, orig_layer, b$plot$data)
    df <- .miami_transform_x_using_meta(df, axis_map)
    df <- .miami_transform_y(df, mode = mode, gap0 = gap0, y_max0 = y_max0, panel_height0 = panel_height0)
    if (identical(mode, "bottom") && isTRUE(direction.reverse) && ("shape" %in% names(df))) {
      shp <- suppressWarnings(as.numeric(as.character(df$shape)))
      if ("direction_sign" %in% names(df)) {
        is_dir <- suppressWarnings(as.numeric(df$direction_sign))
        is_dir <- is.finite(is_dir) & !is.na(is_dir) & is_dir != 0
      } else {
        is_dir <- rep(TRUE, length(shp))
      }
      sw24 <- is_dir & is.finite(shp) & (shp == 24)
      sw25 <- is_dir & is.finite(shp) & (shp == 25)
      if (any(sw24 | sw25)) {
        shp[sw24] <- 25
        shp[sw25] <- 24
        df$shape <- shp
      }
    }
    if ("size" %in% names(df) && !is.numeric(df$size)) df$size <- suppressWarnings(as.numeric(as.character(df$size)))
    if ("linewidth" %in% names(df) && !is.numeric(df$linewidth)) df$linewidth <- suppressWarnings(as.numeric(as.character(df$linewidth)))
    if ("alpha" %in% names(df) && !is.numeric(df$alpha)) df$alpha <- suppressWarnings(as.numeric(as.character(df$alpha)))
    geom_obj <- orig_layer$geom
    geom_name <- class(geom_obj)[1]
    keep_geom_params <- TRUE
    if (grepl("Repel", geom_name, fixed = TRUE)) {
      geom_obj <- if (grepl("Label", geom_name, fixed = TRUE)) ggplot2::GeomLabel else ggplot2::GeomText
      keep_geom_params <- FALSE
    }
    lyr <- ggplot2::layer(
      data = df,
      mapping = .miami_make_identity_mapping(df),
      geom = geom_obj,
      stat = "identity",
      position = "identity",
      inherit.aes = FALSE,
      params = orig_layer$params
    )
    if (!is.null(orig_layer$aes_params)) lyr$aes_params <- orig_layer$aes_params
    if (!is.null(orig_layer$stat_params)) lyr$stat_params <- orig_layer$stat_params
    if (isTRUE(keep_geom_params) && !is.null(orig_layer$geom_params)) lyr$geom_params <- orig_layer$geom_params
    if (identical(mode, "bottom") && isTRUE(direction.reverse)) {
      .swap_shape_24_25 <- function(v) {
        s <- suppressWarnings(as.numeric(v))
        if (!length(s)) return(v)
        i24 <- is.finite(s) & !is.na(s) & s == 24
        i25 <- is.finite(s) & !is.na(s) & s == 25
        if (any(i24 | i25)) {
          s[i24] <- 25
          s[i25] <- 24
        }
        s
      }
      if (!is.null(lyr$aes_params$shape)) lyr$aes_params$shape <- .swap_shape_24_25(lyr$aes_params$shape)
      if (!is.null(lyr$geom_params$shape)) lyr$geom_params$shape <- .swap_shape_24_25(lyr$geom_params$shape)
    }
    if (identical(mode, "bottom")) {
      .flip_just01 <- function(v) {
        x <- suppressWarnings(as.numeric(v))
        ok <- is.finite(x) & !is.na(x) & x >= 0 & x <= 1
        if (any(ok)) x[ok] <- 1 - x[ok]
        x
      }
      if (!is.null(lyr$aes_params$angle)) {
        ang0 <- suppressWarnings(as.numeric(lyr$aes_params$angle))
        if (length(ang0)) lyr$aes_params$angle <- -1 * ang0
      }
      if (!is.null(lyr$geom_params$angle)) {
        ang1 <- suppressWarnings(as.numeric(lyr$geom_params$angle))
        if (length(ang1)) lyr$geom_params$angle <- -1 * ang1
      }
      if (!is.null(lyr$aes_params$vjust)) {
        vj0 <- .flip_just01(lyr$aes_params$vjust)
        if (length(vj0)) lyr$aes_params$vjust <- vj0
      }
      if (!is.null(lyr$geom_params$vjust)) {
        vj1 <- .flip_just01(lyr$geom_params$vjust)
        if (length(vj1)) lyr$geom_params$vjust <- vj1
      }
    }
    layers_out[[length(layers_out) + 1L]] <- lyr
    dfs_out[[length(dfs_out) + 1L]] <- df
  }
  list(layers = layers_out, dfs = dfs_out)
}

.miami_add_identity_scales <- function(p, dfs, show.legend = FALSE) {
  cols <- unique(unlist(lapply(dfs, names), use.names = FALSE))
  guide_opt <- if (isTRUE(show.legend)) "legend" else "none"
  if ("colour" %in% cols) p <- p + ggplot2::scale_colour_identity(guide = guide_opt)
  if ("fill" %in% cols) p <- p + ggplot2::scale_fill_identity(guide = guide_opt)
  if ("alpha" %in% cols) p <- p + ggplot2::scale_alpha_identity(guide = guide_opt)
  if ("linetype" %in% cols) p <- p + ggplot2::scale_linetype_identity(guide = guide_opt)
  if ("shape" %in% cols) p <- p + ggplot2::scale_shape_identity(guide = guide_opt)
  if ("size" %in% cols) p <- p + ggplot2::scale_size_identity(guide = guide_opt)
  if ("linewidth" %in% cols) p <- p + ggplot2::scale_linewidth_identity(guide = guide_opt)
  p
}

.miami_axis_max_from_range <- function(rng) {
  rng <- suppressWarnings(as.numeric(rng))
  if (length(rng) != 2L || !all(is.finite(rng))) return(NA_real_)
  max(rng)
}

.miami_extract_breaks_labels <- function(sc, ymax, source_breaks, source_labels) {
  br_src <- source_breaks
  if (!is.null(sc) && is.function(sc$get_breaks)) {
    br_sc <- tryCatch(sc$get_breaks(), error = function(e) NULL)
    if (!is.null(br_sc)) br_src <- br_sc
  }
  br <- suppressWarnings(as.numeric(br_src))
  br <- br[is.finite(br) & !is.na(br)]
  if (is.finite(ymax)) br <- br[br >= 0 & br <= (ymax + 1e-10)]
  br <- sort(unique(br))
  if (!length(br)) br <- c(0, ymax)
  br <- suppressWarnings(as.numeric(br))
  br <- br[is.finite(br) & !is.na(br)]
  br <- sort(unique(c(0, br)))
  if (!length(br)) br <- 0
  lb <- NULL
  if (!is.null(sc) && is.function(sc$get_labels)) {
    lb <- tryCatch(sc$get_labels(br), error = function(e) NULL)
  }
  if ((is.null(lb) || length(lb) != length(br)) && !is.null(source_labels)) {
    src_br <- suppressWarnings(as.numeric(source_breaks))
    if (length(src_br) && length(source_labels) == length(src_br)) {
      map <- stats::setNames(as.character(source_labels), as.character(src_br))
      lb <- unname(map[as.character(br)])
    }
  }
  if (is.null(lb) || length(lb) != length(br)) lb <- as.character(signif(br, 5))
  lb <- as.character(lb)
  lb[is.na(lb) | !nzchar(lb)] <- as.character(signif(br[is.na(lb) | !nzchar(lb)], 5))
  list(breaks = br, labels = lb)
}

.miami_min_x_from_layers <- function(dfs) {
  vals <- numeric()
  for (df in dfs) {
    if (!is.data.frame(df) || !"x" %in% names(df)) next
    xv <- suppressWarnings(as.numeric(df$x))
    xv <- xv[is.finite(xv) & !is.na(xv)]
    if (length(xv)) vals <- c(vals, xv)
  }
  if (!length(vals)) return(NA_real_)
  min(vals, na.rm = TRUE)
}

.miami_max_x_from_layers <- function(dfs) {
  vals <- numeric()
  for (df in dfs) {
    if (!is.data.frame(df) || !"x" %in% names(df)) next
    xv <- suppressWarnings(as.numeric(df$x))
    xv <- xv[is.finite(xv) & !is.na(xv)]
    if (length(xv)) vals <- c(vals, xv)
  }
  if (!length(vals)) return(NA_real_)
  max(vals, na.rm = TRUE)
}

.miami_digit_count <- function(lbl) {
  if (is.null(lbl)) return(0L)
  x <- as.character(lbl)
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) return(0L)
  x <- gsub("[^0-9]", "", x)
  x <- x[nzchar(x)]
  if (!length(x)) return(0L)
  max(nchar(x), na.rm = TRUE)
}

.miami_y_extent_from_dfs <- function(dfs) {
  vals <- numeric()
  ycols <- c("y", "yend", "ymin", "ymax", "yintercept")
  for (df in dfs) {
    if (!is.data.frame(df) || !nrow(df)) next
    for (nm in ycols) {
      if (!nm %in% names(df)) next
      v <- suppressWarnings(as.numeric(df[[nm]]))
      v <- v[is.finite(v) & !is.na(v)]
      if (length(v)) vals <- c(vals, v)
    }
  }
  if (!length(vals)) return(c(NA_real_, NA_real_))
  c(min(vals, na.rm = TRUE), max(vals, na.rm = TRUE))
}

.miami_ytitle_text <- function(lbl, meta = NULL) {
  if (is.language(lbl) || inherits(lbl, "expression")) {
    return(list(label = paste(deparse(lbl), collapse = ""), parse = TRUE))
  }
  if (is.character(lbl) && length(lbl) > 0L && !is.na(lbl[1]) && nzchar(lbl[1])) {
    return(list(label = as.character(lbl)[1], parse = FALSE))
  }
  if (is.list(meta) && is.list(meta$columns) && !is.null(meta$columns$y.col)) {
    ycol <- as.character(meta$columns$y.col)[1]
    if (!is.na(ycol) && nzchar(ycol)) return(list(label = ycol, parse = FALSE))
  }
  list(label = "paste(-log[10], ' (', italic(P), ')')", parse = TRUE)
}

.miami_space_to_pt <- function(x, unit_pt) {
  x <- suppressWarnings(as.numeric(x))[1]
  if (!is.finite(x) || is.na(x) || x < 0) return(NA_real_)
  as.numeric(x * unit_pt)
}

.miami_prepare_context <- function(p1, p2, meta1, meta2, gap, gap.ratio, y.ratio,
                                           drop.internal, direction.reverse) {
  axis_map <- .miami_axis_union(.miami_axis_from_meta(meta1), .miami_axis_from_meta(meta2))
  obs_chr_union <- unique(c(.miami_observed_chr(p1), .miami_observed_chr(p2)))
  obs_chr_union <- obs_chr_union[!is.na(obs_chr_union) & nzchar(obs_chr_union)]
  if (!is.null(axis_map) && nrow(axis_map) && length(obs_chr_union)) {
    axis_map <- axis_map[chr %in% obs_chr_union]
    if (nrow(axis_map)) {
      axis_map[, chr_order := rank.chrom(chr)]
      data.table::setorderv(axis_map, c("chr_order", "chr"), c(1L, 1L), na.last = TRUE)
      axis_map[, width := pmax(1, end - start + 1)]
      axis_map[, offset := data.table::shift(cumsum(width), fill = 0)]
      axis_map[, tick := offset + width / 2]
      axis_map[, chr_order := NULL]
    }
  }

  # Build each 1e6-point plot once and reuse; the miami helpers below otherwise
  # call ggplot_build() again (3x per plot), which dominates miami() runtime.
  b1 <- suppressWarnings(ggplot2::ggplot_build(p1))
  b2 <- suppressWarnings(ggplot2::ggplot_build(p2))
  y1 <- .miami_get_axis_info(p1, "y", built = b1)
  y2 <- .miami_get_axis_info(p2, "y", built = b2)
  sc1 <- b1$layout$panel_params[[1]]$y
  sc2 <- b2$layout$panel_params[[1]]$y

  y1_max <- .miami_axis_max_from_range(b1$layout$panel_params[[1]]$y.range %||% y1$range)
  y2_max <- .miami_axis_max_from_range(b2$layout$panel_params[[1]]$y.range %||% y2$range)
  if (!is.finite(y1_max) || is.na(y1_max)) y1_max <- .miami_point_ymax(p1, built = b1)
  if (!is.finite(y2_max) || is.na(y2_max)) y2_max <- .miami_point_ymax(p2, built = b2)
  if (!is.finite(y1_max) || is.na(y1_max) || y1_max <= 0) y1_max <- 1
  if (!is.finite(y2_max) || is.na(y2_max) || y2_max <= 0) y2_max <- 1

  if (is.null(gap)) {
    gap_ratio0 <- suppressWarnings(as.numeric(gap.ratio))[1]
    if (!is.finite(gap_ratio0) || is.na(gap_ratio0) || gap_ratio0 <= 0 || gap_ratio0 > 0.5) gap_ratio0 <- 0.08
    gap0 <- max(y1_max, y2_max) * gap_ratio0
  } else {
    gap0 <- suppressWarnings(as.numeric(gap))[1]
  }
  if (!is.finite(gap0) || is.na(gap0) || gap0 <= 0) gap0 <- max(y1_max, y2_max) * 0.08
  half0 <- gap0 / 2

  panel_height_top <- y1_max * y.ratio[1]
  panel_height_bottom <- y2_max * y.ratio[2]
  if (!is.finite(panel_height_top) || is.na(panel_height_top) || panel_height_top <= 0) panel_height_top <- y1_max
  if (!is.finite(panel_height_bottom) || is.na(panel_height_bottom) || panel_height_bottom <= 0) panel_height_bottom <- y2_max

  top <- .miami_extract_layers_built(
    p1, mode = "top", axis_map = axis_map, gap0 = gap0,
    y_max0 = y1_max, panel_height0 = panel_height_top,
    drop.internal0 = drop.internal, direction.reverse = direction.reverse, built = b1
  )
  bot <- .miami_extract_layers_built(
    p2, mode = "bottom", axis_map = axis_map, gap0 = gap0,
    y_max0 = y2_max, panel_height0 = panel_height_bottom,
    drop.internal0 = drop.internal, direction.reverse = direction.reverse, built = b2
  )

  bl1 <- .miami_extract_breaks_labels(sc1, y1_max, y1$breaks, y1$labels)
  bl2 <- .miami_extract_breaks_labels(sc2, y2_max, y2$breaks, y2$labels)
  br1 <- bl1$breaks
  br2 <- bl2$breaks
  top_breaks <- (br1 / y1_max) * panel_height_top + half0
  bot_breaks <- -((br2 / y2_max) * panel_height_bottom) - half0
  y_breaks <- c(bot_breaks, top_breaks)
  y_labels <- c(bl2$labels, bl1$labels)
  ord <- order(y_breaks)
  y_breaks <- y_breaks[ord]
  y_labels <- y_labels[ord]
  y_minor <- numeric()

  x_label <- p1$labels$x %||% "Chromosome"
  x_lim <- NULL
  axis_dt <- NULL
  if (!is.null(axis_map) && nrow(axis_map)) {
    axis_dt <- data.table::data.table(CHR = as.character(axis_map$chr), x_tick = as.numeric(axis_map$tick))
    x_max <- max(axis_map$offset + axis_map$width, na.rm = TRUE)
    x_pad <- max(1, x_max * 0.01)
    x_lim <- c(0, x_max + x_pad)
  } else {
    x_info1 <- .miami_get_axis_info(p1, "x", built = b1)
    x_lim <- x_info1$range
    if (is.null(x_lim) || length(x_lim) != 2 || !all(is.finite(x_lim))) x_lim <- c(0, 1)
    if (!is.finite(x_lim[2]) || is.na(x_lim[2]) || x_lim[2] <= 0) x_lim[2] <- 1
    x_lim[1] <- 0
  }
  x_span0 <- diff(x_lim)
  if (!is.finite(x_span0) || is.na(x_span0) || x_span0 <= 0) x_span0 <- 1

  x_data_min <- .miami_min_x_from_layers(c(top$dfs, bot$dfs))
  x_data_max <- .miami_max_x_from_layers(c(top$dfs, bot$dfs))
  if (!is.finite(x_data_min) || is.na(x_data_min)) {
    if (!is.null(axis_map) && nrow(axis_map)) {
      x_data_min <- min(axis_map$offset + 1, na.rm = TRUE)
    } else {
      x_data_min <- x_lim[1] + max(0.5, x_span0 * 0.005)
    }
  }
  if (!is.finite(x_data_min) || is.na(x_data_min)) x_data_min <- 1
  if (!is.finite(x_data_max) || is.na(x_data_max) || x_data_max <= x_data_min) x_data_max <- x_data_min + x_span0
  x_data_span <- x_data_max - x_data_min
  if (!is.finite(x_data_span) || is.na(x_data_span) || x_data_span <= 0) x_data_span <- x_span0
  axis_left_gap <- max(0.5, x_data_span * 0.01)
  y_axis_x <- x_data_min - axis_left_gap

  y_digit_max <- max(
    .miami_digit_count(bl1$labels),
    .miami_digit_count(bl2$labels),
    nchar(as.character(max(0, suppressWarnings(as.integer(floor(y1_max))), suppressWarnings(as.integer(floor(y2_max))), na.rm = TRUE))),
    na.rm = TRUE
  )
  if (!is.finite(y_digit_max) || is.na(y_digit_max) || y_digit_max < 1) y_digit_max <- 1
  y_title_x <- y_axis_x
  y_title_hjust <- 0.5
  y_title_vjust <- -(0.9 + 0.45 * y_digit_max)
  x_lim_draw <- c(min(x_lim[1], y_axis_x), x_lim[2])
  if (!all(is.finite(x_lim_draw)) || diff(x_lim_draw) <= 0) x_lim_draw <- x_lim

  y_lim_base <- c(-half0 - panel_height_bottom, half0 + panel_height_top)
  y_ext <- .miami_y_extent_from_dfs(c(top$dfs, bot$dfs))

  list(
    axis_map = axis_map,
    axis_dt = axis_dt,
    y1 = y1,
    y2 = y2,
    y1_max = y1_max,
    y2_max = y2_max,
    gap = gap0,
    half = half0,
    panel_height_top = panel_height_top,
    panel_height_bottom = panel_height_bottom,
    top = top,
    bot = bot,
    y_breaks = y_breaks,
    y_labels = y_labels,
    y_minor = y_minor,
    x_label = x_label,
    x_lim = x_lim,
    x_lim_draw = x_lim_draw,
    x_span = x_span0,
    x_data_min = x_data_min,
    x_data_max = x_data_max,
    y_axis_x = y_axis_x,
    y_digit_max = y_digit_max,
    y_title_x = y_title_x,
    y_title_hjust = y_title_hjust,
    y_title_vjust = y_title_vjust,
    y_lim = y_lim_base,
    y_ext = y_ext
  )
}

.miami_compose_plot <- function(p1, p2, meta1, meta2, ctx, show.legend,
                                        panel.box, x.text.size, y.text.size, tick.size) {
  # Strip the heavy (~1e6-row) layers/data BEFORE deep-copying. Replacing these
  # slots does not mutate p1, and .miami_remove_scales() then serialises only a
  # lightweight skeleton instead of the full point cloud (the layers are discarded
  # here anyway). Avoids ~15s of serialize/unserialize on large plots.
  base <- p1
  base$layers <- list()
  base$data <- data.frame()
  base$mapping <- ggplot2::aes()
  base <- .miami_remove_scales(base, c("x", "y", "shape", "size", "colour", "fill", "alpha", "linetype", "linewidth"))
  for (ly in c(ctx$top$layers, ctx$bot$layers)) base <- base + ly
  base <- .miami_add_identity_scales(base, dfs = c(ctx$top$dfs, ctx$bot$dfs), show.legend = show.legend)

  y_title_top <- p1$labels$y %||% NULL
  y_title_bottom <- p2$labels$y %||% NULL
  base <- base + ggplot2::scale_y_continuous(
    name = NULL,
    breaks = ctx$y_breaks,
    labels = ctx$y_labels,
    minor_breaks = ctx$y_minor,
    expand = ggplot2::expansion(mult = c(0, 0), add = c(0, 0))
  )
  if (!is.null(ctx$axis_dt) && nrow(ctx$axis_dt)) {
    base <- base + ggplot2::scale_x_continuous(
      name = NULL,
      breaks = ctx$axis_dt$x_tick,
      labels = rep("", nrow(ctx$axis_dt)),
      limits = ctx$x_lim_draw,
      expand = c(0, 0)
    )
  }
  base$coordinates <- ggplot2::coord_cartesian(xlim = ctx$x_lim_draw, ylim = ctx$y_lim, clip = "off")

  line_col <- .miami_theme_color(p1, "axis.line.x", default = "grey20")
  line_lw <- .miami_theme_linewidth(p1, "axis.line.x", default = 0.35)
  y_text_pt <- if (is.finite(y.text.size) && !is.na(y.text.size) && y.text.size > 0) y.text.size else .miami_theme_text_pt(p1, "axis.text.x", default = 14)
  midline_dt <- data.table::data.table(
    x = c(ctx$y_axis_x, ctx$y_axis_x),
    xend = c(ctx$x_lim[2], ctx$x_lim[2]),
    y = c(+ctx$half, -ctx$half),
    yend = c(+ctx$half, -ctx$half)
  )
  base <- base +
    ggplot2::geom_segment(data = midline_dt, ggplot2::aes(x = x, y = y, xend = xend, yend = yend), inherit.aes = FALSE, colour = line_col, linewidth = line_lw) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = y_text_pt),
      axis.ticks.x = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank()
    )

  if (isTRUE(panel.box) && length(ctx$x_lim) == 2L && all(is.finite(ctx$x_lim)) && length(ctx$y_lim) == 2L && all(is.finite(ctx$y_lim))) {
    box_top <- data.table::data.table(
      x = c(ctx$y_axis_x, ctx$y_axis_x, ctx$y_axis_x, ctx$x_lim[2]),
      y = c(ctx$half, ctx$half + ctx$panel_height_top, ctx$half, ctx$half),
      xend = c(ctx$x_lim[2], ctx$x_lim[2], ctx$y_axis_x, ctx$x_lim[2]),
      yend = c(ctx$half, ctx$half + ctx$panel_height_top, ctx$half + ctx$panel_height_top, ctx$half + ctx$panel_height_top)
    )
    box_bottom <- data.table::data.table(
      x = c(ctx$y_axis_x, ctx$y_axis_x, ctx$y_axis_x, ctx$x_lim[2]),
      y = c(-ctx$half - ctx$panel_height_bottom, -ctx$half, -ctx$half - ctx$panel_height_bottom, -ctx$half - ctx$panel_height_bottom),
      xend = c(ctx$x_lim[2], ctx$x_lim[2], ctx$y_axis_x, ctx$x_lim[2]),
      yend = c(-ctx$half - ctx$panel_height_bottom, -ctx$half, -ctx$half, -ctx$half)
    )
    base <- base + ggplot2::geom_segment(
      data = data.table::rbindlist(list(box_top, box_bottom), use.names = TRUE),
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      inherit.aes = FALSE,
      colour = line_col,
      linewidth = line_lw
    )
  }

  if (!is.null(ctx$axis_dt) && nrow(ctx$axis_dt)) {
    y_axis_dt <- data.table::data.table(
      x = c(ctx$y_axis_x, ctx$y_axis_x),
      xend = c(ctx$y_axis_x, ctx$y_axis_x),
      y = c(ctx$half, -ctx$half - ctx$panel_height_bottom),
      yend = c(ctx$half + ctx$panel_height_top, -ctx$half)
    )
    base <- base + ggplot2::geom_segment(
      data = y_axis_dt,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      inherit.aes = FALSE,
      colour = line_col,
      linewidth = line_lw
    )
  }

  if (!is.null(ctx$axis_dt) && nrow(ctx$axis_dt)) {
    txt_pt <- if (is.finite(x.text.size) && !is.na(x.text.size) && x.text.size > 0) x.text.size else .miami_theme_text_pt(p1, "axis.text.x", default = 14)
    ttl_pt <- .miami_theme_text_pt(p1, "axis.title.x", default = 16)
    txt_mm <- txt_pt / 2.845
    ttl_mm <- (ttl_pt * 0.7) / 2.845
    tick_len <- if (is.finite(tick.size) && !is.na(tick.size) && tick.size >= 0) tick.size else max(ctx$gap * 0.16, 0.16)
    y_chr <- 0
    base <- base +
      ggplot2::geom_text(data = ctx$axis_dt, ggplot2::aes(x = x_tick, y = y_chr, label = CHR), inherit.aes = FALSE, colour = line_col, size = txt_mm, vjust = 0.5) +
      ggplot2::annotate("text", x = ctx$y_axis_x, y = y_chr, label = as.character(ctx$x_label), hjust = 0.75, vjust = 0.5, colour = line_col, size = ttl_mm)
    if (is.finite(tick_len) && !is.na(tick_len) && tick_len > 0) {
      base <- base +
        ggplot2::geom_segment(data = ctx$axis_dt, ggplot2::aes(x = x_tick, xend = x_tick, y = ctx$half, yend = ctx$half - tick_len), inherit.aes = FALSE, colour = line_col, linewidth = line_lw) +
        ggplot2::geom_segment(data = ctx$axis_dt, ggplot2::aes(x = x_tick, xend = x_tick, y = -ctx$half, yend = -ctx$half + tick_len), inherit.aes = FALSE, colour = line_col, linewidth = line_lw)
    }
  }

  top_ttl <- .miami_ytitle_text(y_title_top, meta = meta1)
  bot_ttl <- .miami_ytitle_text(y_title_bottom, meta = meta2)
  if (!is.null(top_ttl$label) && nzchar(top_ttl$label)) {
    base <- base + ggplot2::annotate("text", x = ctx$y_title_x, y = ctx$half + (ctx$panel_height_top / 2), label = top_ttl$label, parse = isTRUE(top_ttl$parse), angle = 90, hjust = ctx$y_title_hjust, vjust = ctx$y_title_vjust, colour = line_col, size = .miami_theme_text_pt(p1, "axis.title.x", default = 16) / 2.845)
  }
  if (!is.null(bot_ttl$label) && nzchar(bot_ttl$label)) {
    base <- base + ggplot2::annotate("text", x = ctx$y_title_x, y = -ctx$half - (ctx$panel_height_bottom / 2), label = bot_ttl$label, parse = isTRUE(bot_ttl$parse), angle = 90, hjust = ctx$y_title_hjust, vjust = ctx$y_title_vjust, colour = line_col, size = .miami_theme_text_pt(p2, "axis.title.x", default = 16) / 2.845)
  }
  base
}

.miami_finalize_plot <- function(base, ctx, meta1, meta2, y.ratio, direction.reverse,
                                         panel.space, panel.space.top, panel.space.bottom,
                                         panel.space.left, panel_space_unit_pt, silent) {
  top_over <- 0
  bottom_over <- 0
  if (all(is.finite(ctx$y_ext))) {
    top_over <- max(0, ctx$y_ext[2] - (ctx$half + ctx$panel_height_top))
    bottom_over <- max(0, (-ctx$half - ctx$panel_height_bottom) - ctx$y_ext[1])
  }
  top_ratio <- top_over / max(1, ctx$panel_height_top)
  bottom_ratio <- bottom_over / max(1, ctx$panel_height_bottom)
  top_extra_pt <- min(320, 6 + top_ratio * 180)
  bottom_extra_pt <- min(320, 6 + bottom_ratio * 180)
  if (!is.na(panel.space[1])) top_extra_pt <- .miami_space_to_pt(panel.space[1], panel_space_unit_pt)
  if (!is.na(panel.space[2])) bottom_extra_pt <- .miami_space_to_pt(panel.space[2], panel_space_unit_pt)
  if (!is.na(panel.space.top)) top_extra_pt <- .miami_space_to_pt(panel.space.top, panel_space_unit_pt)
  if (!is.na(panel.space.bottom)) bottom_extra_pt <- .miami_space_to_pt(panel.space.bottom, panel_space_unit_pt)
  left_margin_pt <- max(44, 24 + 7 * ctx$y_digit_max)
  if (!is.na(panel.space.left)) left_margin_pt <- .miami_space_to_pt(panel.space.left, panel_space_unit_pt)
  base <- base + ggplot2::theme(plot.margin = ggplot2::margin(10 + top_extra_pt, 10, 10 + bottom_extra_pt, left_margin_pt))
  attr(base, "gcanvas_meta") <- list(
    type = "miami",
    source = list(p1 = meta1, p2 = meta2),
    gap = as.numeric(ctx$gap),
    y.ratio = as.numeric(y.ratio),
    direction.reverse = as.logical(direction.reverse),
    panel_height = c(top = as.numeric(ctx$panel_height_top), bottom = as.numeric(ctx$panel_height_bottom)),
    x_label = as.character(ctx$x_label),
    n_chromosomes = if (!is.null(ctx$axis_dt)) as_int(nrow(ctx$axis_dt)) else NA_integer_
  )
  .gcanvas_note("gcanvas::miami", "Done", silent = silent)
  base
}

