# Regional association plot: the lead-variant panel + LD coloring + gene
# track layout (via `regional.track()`) + optional credible-set markers.

.regional_draw_ld_points <- function(p,
                                             df0,
                                             use_snp2,
                                             legend.position,
                                             lead_snp,
                                             lead_snp2,
                                             lead1_color,
                                             lead2_color,
                                             lead_shape,
                                             lead_shape2,
                                             ldcolorset,
                                             alpha,
                                             shared_mode,
                                             shared.color,
                                             shared.cutoff,
                                             shared.alpha) {
  ld_is_horiz <- legend.position %in% c("bottom", "top")
  ld_barwidth <- if (ld_is_horiz) grid::unit(5.5, "cm") else grid::unit(0.45, "cm")
  ld_barheight <- if (ld_is_horiz) grid::unit(0.45, "cm") else grid::unit(3, "cm")
  if (isTRUE(use_snp2)) {
    require_pkg(c("ggnewscale", "colorspace"))
    .mix_hex <- function(col_a, col_b, w = 0.5) {
      w <- pmin(pmax(suppressWarnings(as.numeric(w)), 0), 1)
      a <- grDevices::col2rgb(col_a)[, 1]
      b <- grDevices::col2rgb(col_b)[, 1]
      r <- round(a[1] * (1 - w) + b[1] * w)
      g <- round(a[2] * (1 - w) + b[2] * w)
      bb <- round(a[3] * (1 - w) + b[3] * w)
      grDevices::rgb(r, g, bb, maxColorValue = 255)
    }
    .ld_pal <- function(lead_col) {
      base_pal <- grDevices::colorRampPalette(c(.mix_hex(lead_col, "#F2F2F2", w = 0.97), lead_col))(5)
      c(.mix_hex(lead_col, "#F2F2F2", w = 0.985), .mix_hex(lead_col, "#F2F2F2", w = 0.935), base_pal[3], base_pal[4], base_pal[5])
    }
    df0[, `:=`(
      ld1_plot = pmin(pmax(suppressWarnings(as.numeric(ld1)), 0), 1),
      ld2_plot = pmin(pmax(suppressWarnings(as.numeric(ld2)), 0), 1)
    )]
    df0[, ld_group := 1L]
    df0[is.finite(ld2_plot) & !is.finite(ld1_plot), ld_group := 2L]
    df0[is.finite(ld1_plot) & is.finite(ld2_plot) & ld2_plot > ld1_plot, ld_group := 2L]
    d1 <- df0[ld_group == 1L]
    d2 <- df0[ld_group == 2L]
    d_shared <- df0[0]
    pal1 <- .ld_pal(lead1_color)
    pal2 <- .ld_pal(lead2_color)
    pal_shared <- .ld_pal(shared.color)
    ld_breaks <- seq(0, 1, by = 0.2)
    if (!identical(shared_mode, "off")) {
      d_shared <- df0[is.finite(ld1_plot) & is.finite(ld2_plot)]
      if (nrow(d_shared)) {
        d_shared[, shared_r2 := pmin(ld1_plot, ld2_plot)]
        d_shared <- d_shared[shared_r2 >= shared.cutoff]
        if (identical(shared_mode, "split")) {
          shared_snp <- unique(d_shared$snp)
          d1 <- d1[!(snp %in% shared_snp)]
          d2 <- d2[!(snp %in% shared_snp)]
        }
      }
    }
    data.table::setorderv(d1, c("yval", "xvalue", "snp"), c(1L, 1L, 1L), na.last = TRUE)
    data.table::setorderv(d2, c("yval", "xvalue", "snp"), c(1L, 1L, 1L), na.last = TRUE)
    if (nrow(d_shared)) data.table::setorderv(d_shared, c("yval", "xvalue", "snp"), c(1L, 1L, 1L), na.last = TRUE)
    legend_dummy <- data.table::data.table(x = min(df0$xvalue, na.rm = TRUE), y = min(df0$yval, na.rm = TRUE), ld_legend = c(0, 1))
    p <- p +
      ggplot2::geom_point(data = d1, ggplot2::aes(x = xvalue, y = yval, fill = ld1_plot, shape = point_shape), color = "grey20", size = 3, alpha = alpha, show.legend = FALSE) +
      ggplot2::geom_point(data = legend_dummy, ggplot2::aes(x = x, y = y, fill = ld_legend), inherit.aes = FALSE, shape = 21, size = 0, alpha = 0, stroke = 0, show.legend = TRUE) +
      ggplot2::scale_fill_stepsn(colors = pal1, limits = c(0, 1), breaks = ld_breaks, show.limits = TRUE, na.value = "grey35",
        guide = ggplot2::guide_coloursteps(title = if (is.null(lead_snp) || is.na(lead_snp)) bquote(r[1]^2) else bquote(r[1]^2 ~ "(" * .(lead_snp) * ")"),
          title.theme = ggplot2::element_text(face = "bold", size = 14),
          label.theme = ggplot2::element_text(size = 11),
          title.position = "top",
          title.hjust = 0.5,
          direction = if (ld_is_horiz) "horizontal" else "vertical",
          barwidth = ld_barwidth,
          barheight = ld_barheight,
          frame.colour = "black",
          frame.linewidth = 0.4,
          ticks.colour = "black",
          ticks.linewidth = 0.4,
          order = 1)) +
      ggnewscale::new_scale_fill() +
      ggplot2::geom_point(data = d2, ggplot2::aes(x = xvalue, y = yval, fill = ld2_plot, shape = point_shape), color = "grey20", size = 3, alpha = alpha, show.legend = FALSE) +
      ggplot2::geom_point(data = legend_dummy, ggplot2::aes(x = x, y = y, fill = ld_legend), inherit.aes = FALSE, shape = 21, size = 0, alpha = 0, stroke = 0, show.legend = TRUE) +
      ggplot2::scale_fill_stepsn(colors = pal2, limits = c(0, 1), breaks = ld_breaks, show.limits = TRUE, na.value = "grey35",
        guide = ggplot2::guide_coloursteps(title = if (is.null(lead_snp2) || is.na(lead_snp2)) bquote(r[2]^2) else bquote(r[2]^2 ~ "(" * .(lead_snp2) * ")"),
          title.theme = ggplot2::element_text(face = "bold", size = 14),
          label.theme = ggplot2::element_text(size = 11),
          title.position = "top",
          title.hjust = 0.5,
          direction = if (ld_is_horiz) "horizontal" else "vertical",
          barwidth = ld_barwidth,
          barheight = ld_barheight,
          frame.colour = "black",
          frame.linewidth = 0.4,
          ticks.colour = "black",
          ticks.linewidth = 0.4,
          order = 2))
    if (!identical(shared_mode, "off")) {
      if (identical(shared_mode, "split")) {
        split_size <- 4.2
        split_dx <- (max(df0$xvalue, na.rm = TRUE) - min(df0$xvalue, na.rm = TRUE)) * 0.0002
        .ld_bin5 <- function(x) {
          out <- as.integer(cut(x, breaks = ld_breaks, include.lowest = TRUE, right = TRUE, labels = FALSE))
          out[is.na(out)] <- 1L
          pmax(1L, pmin(5L, out))
        }
        d_shared[, `:=`(ld1_bin = .ld_bin5(ld1_plot), ld2_bin = .ld_bin5(ld2_plot), draw_id = seq_len(.N))]
        d_shared_split <- data.table::rbindlist(list(
          d_shared[, .(xvalue = xvalue + split_dx, yval = yval, draw_id = draw_id, part = 1L, glyph = "\u25D0", col_hex = pal1[ld1_bin])],
          d_shared[, .(xvalue = xvalue - split_dx, yval = yval, draw_id = draw_id, part = 2L, glyph = "\u25D1", col_hex = pal2[ld2_bin])],
          d_shared[, .(xvalue = xvalue, yval = yval, draw_id = draw_id, part = 3L, glyph = "\u25CB", col_hex = "grey20")]
        ), use.names = TRUE)
        data.table::setorderv(d_shared_split, c("draw_id", "part"), c(1L, 1L), na.last = TRUE)
        p <- p +
          ggplot2::geom_text(data = d_shared_split, ggplot2::aes(x = xvalue, y = yval, label = glyph, color = col_hex), family = "Apple Symbols", size = split_size, alpha = shared.alpha, show.legend = FALSE) +
          ggplot2::scale_colour_identity(guide = "none")
      } else {
        p <- p +
          ggnewscale::new_scale_fill() +
          ggplot2::geom_point(data = d_shared, ggplot2::aes(x = xvalue, y = yval, fill = shared_r2, shape = point_shape), color = "grey20", size = 3, alpha = shared.alpha, show.legend = FALSE) +
          ggplot2::geom_point(data = legend_dummy, ggplot2::aes(x = x, y = y, fill = ld_legend), inherit.aes = FALSE, shape = 21, size = 0, alpha = 0, stroke = 0, show.legend = TRUE) +
          ggplot2::scale_fill_stepsn(colors = pal_shared, limits = c(0, 1), breaks = ld_breaks, show.limits = TRUE, na.value = "grey35",
            guide = ggplot2::guide_coloursteps(title = bquote(paste("Shared LD: min(", r[1]^2, ", ", r[2]^2, ") > ", .(shared.cutoff))),
              title.theme = ggplot2::element_text(face = "bold", size = 14),
              label.theme = ggplot2::element_text(size = 11),
              title.position = "top",
              title.hjust = 0.5,
              direction = if (ld_is_horiz) "horizontal" else "vertical",
              barwidth = ld_barwidth,
              barheight = ld_barheight,
              frame.colour = "black",
              frame.linewidth = 0.4,
              ticks.colour = "black",
              ticks.linewidth = 0.4,
              order = 3))
      }
    }
  } else {
    p <- p +
      ggplot2::geom_point(data = df0, ggplot2::aes(x = xvalue, y = yval, fill = ld, shape = point_shape), color = "grey20", size = 3, alpha = alpha) +
      ggplot2::scale_fill_gradientn(colors = ldcolorset, limits = c(0, 1), breaks = seq(0, 1, by = 0.2), na.value = "grey80",
        guide = ggplot2::guide_colourbar(title = expression(r^2),
          title.theme = ggplot2::element_text(face = "bold", size = 14),
          label.theme = ggplot2::element_text(size = 11),
          title.position = "top",
          title.hjust = 0.5,
          direction = if (ld_is_horiz) "horizontal" else "vertical",
          barwidth = ld_barwidth,
          barheight = ld_barheight,
          frame.colour = "black",
          frame.linewidth = 0.4,
          ticks.colour = "black",
          ticks.linewidth = 0.4))
  }
  p <- p + ggplot2::geom_point(data = df0[snp == lead_snp], ggplot2::aes(x = xvalue, y = yval), fill = "white", color = lead1_color, size = 4, shape = lead_shape, stroke = 2)
  if (isTRUE(use_snp2)) {
    p <- p + ggplot2::geom_point(data = df0[snp == lead_snp2], ggplot2::aes(x = xvalue, y = yval), fill = "white", color = lead2_color, size = 4, shape = lead_shape2, stroke = 2)
  }
  p
}

.regional_apply_legend_theme <- function(p, threshold_lines, legend.position, use_snp2, title) {
  if (nrow(threshold_lines)) {
    data.table::setorderv(threshold_lines, c("idx"), c(1L), na.last = TRUE)
    for (i in seq_len(nrow(threshold_lines))) {
      p <- p + ggplot2::geom_hline(yintercept = threshold_lines$y[i], linetype = threshold_lines$type[i], color = threshold_lines$color[i], linewidth = threshold_lines$linewidth[i], alpha = 0.8)
    }
  }
  if (legend.position == "bottom") {
    p <- p + ggplot2::theme(legend.position = "bottom", legend.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0), legend.box.margin = ggplot2::margin(t = -4, r = 0, b = 0, l = 0), legend.spacing.y = grid::unit(0.05, "cm"))
  } else if (legend.position == "top") {
    p <- p + ggplot2::theme(legend.position = "top", legend.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0), legend.box.margin = ggplot2::margin(t = 0, r = 0, b = -4, l = 0), legend.spacing.y = grid::unit(0.05, "cm"))
  } else if (legend.position == "left.outside") {
    p <- p + ggplot2::theme(legend.position = "left")
  } else if (legend.position == "right.outside") {
    p <- p + ggplot2::theme(legend.position = "right")
  } else if (legend.position == "left") {
    p <- p + ggplot2::theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1), legend.background = ggplot2::element_rect(fill = scales::alpha("white", 0.4), color = "grey20"))
  } else if (legend.position == "right") {
    p <- p + ggplot2::theme(legend.position = c(0.98, 0.98), legend.justification = c(1, 1), legend.background = ggplot2::element_rect(fill = scales::alpha("white", 0.4), color = "grey20"))
  }
  if (isTRUE(use_snp2)) p <- p + ggplot2::theme(legend.box = "vertical")
  if (!is.null(title)) p <- p + ggplot2::ggtitle(title)
  p
}

.regional_draw_tracks <- function(p, show.gene, tracks, hl_map, show.exon, show.gene.label, show.omit_label, gene_meta, reduc, xbound, y.min, y.max0, label_max_y, hl.gene.size, .compute_y_breaks_reg, exon.size) {
  ylims <- NULL
  y.breaks_draw <- NULL
  y.breaks_labels <- NULL
  if (isTRUE(show.gene) && !is.null(tracks) && nrow(tracks$gene)) {
    gene_draw <- tracks$gene
    exon_draw <- tracks$exon
    y_anchor <- if (nrow(gene_draw)) gene_draw$gene_y else tracks$gene$gene_y
    y.min2 <- min(y.min, min(y_anchor, na.rm = TRUE))
    gap <- diff(c(y.min, y.max0)) * 0.05
    y_top <- y.max0 * 1.1
    if (is.finite(label_max_y)) y_top <- max(y_top, label_max_y * 1.05)
    ylims <- c(y.min2 - gap, y_top)
    bumper <- diff(range(xbound * reduc)) * 1e-3 * exon.size
    if (nrow(gene_draw)) p <- p + ggplot2::geom_segment(data = gene_draw, ggplot2::aes(x = start / reduc, xend = end / reduc, y = gene_y, yend = gene_y), linewidth = 0.7, color = "grey20")
    if (isTRUE(show.exon) && !is.null(exon_draw) && nrow(exon_draw)) {
      p <- p + ggplot2::geom_segment(data = exon_draw, ggplot2::aes(x = (start - bumper) / reduc, xend = (end + bumper) / reduc, y = gene_y, yend = gene_y), linewidth = 2.1, color = "grey20")
    }
    if (!is.null(hl_map) && nrow(gene_draw) && any(gene_draw$hl, na.rm = TRUE)) {
      p <- p + ggplot2::geom_segment(data = gene_draw[hl %in% TRUE], ggplot2::aes(x = start / reduc, xend = end / reduc, y = gene_y, yend = gene_y, color = I(hl_col)), linewidth = 0.9)
      if (isTRUE(show.exon) && !is.null(exon_draw) && nrow(exon_draw) && any(exon_draw$hl, na.rm = TRUE)) {
        p <- p + ggplot2::geom_segment(data = exon_draw[hl %in% TRUE], ggplot2::aes(x = (start - bumper) / reduc, xend = (end + bumper) / reduc, y = gene_y, yend = gene_y, color = I(hl_col)), linewidth = 2.2)
      }
    }
    if (isTRUE(show.gene.label)) {
      tracks$gene$label <- .gcanvas_as_char_no_null(tracks$gene[["label"]], empty = "")
      tracks$gene$gene_name <- .gcanvas_as_char_no_null(tracks$gene[["gene_name"]], empty = "")
      if ("gene_id" %in% names(tracks$gene)) {
        tracks$gene$gene_id <- .gcanvas_as_char_no_null(tracks$gene[["gene_id"]], empty = "")
      } else {
        tracks$gene$gene_id <- ""
      }
      tracks$gene <- tracks$gene[nzchar(gene_name) | nzchar(gene_id)]
      tr_text <- tracks$gene[nzchar(label)]
      if (nrow(tr_text)) {
        if (!is.null(hl_map) && nrow(tracks$gene) && any(tracks$gene$hl, na.rm = TRUE)) {
          p <- p + ggplot2::geom_text(data = tr_text[!(hl %in% TRUE)], ggplot2::aes(x = text_x / reduc, y = text_y, label = label, hjust = hjust), vjust = "bottom", size = 3, fontface = "bold.italic", color = "grey20", na.rm = TRUE)
          p <- p + ggplot2::geom_text(data = tr_text[hl %in% TRUE], ggplot2::aes(x = text_x / reduc, y = text_y, label = label, hjust = hjust, color = I(hl_col)), vjust = "bottom", size = hl.gene.size, fontface = "bold.italic", na.rm = TRUE)
        } else {
          p <- p + ggplot2::geom_text(data = tr_text, ggplot2::aes(x = text_x / reduc, y = text_y, label = label, hjust = hjust), vjust = "bottom", size = 3, fontface = "bold.italic", color = "grey20", na.rm = TRUE)
        }
      }
    }
    y.breaks_obj <- .compute_y_breaks_reg(ylims)
    y.breaks_draw <- y.breaks_obj$draw
    y.breaks_labels <- y.breaks_obj$labels
    if (isTRUE(show.omit_label) && !is.null(gene_meta) && isTRUE(gene_meta$n_omitted > 0)) {
      omit_txt <- sprintf("+%d genes omitted", gene_meta$n_omitted)
      y_omit <- min(tracks$gene$gene_y, na.rm = TRUE) - (0.15 * y.max0) * 0.8
      p <- p + ggplot2::annotate("text", x = xbound[1], y = y_omit, label = omit_txt, hjust = 0, vjust = 1, size = 3, color = "grey30")
    }
  }
  list(p = p, ylims = ylims, y.breaks_draw = y.breaks_draw, y.breaks_labels = y.breaks_labels, tracks = tracks)
}

.regional_draw_labels <- function(p, label.mode, label_df, lead_snps, lead_snp, lead_snp2, use_snp2, label.seed, label.size, lead1_color, lead2_color) {
  if (label.mode != "off" && nrow(label_df)) {
    lab_other <- label_df[!(snp %in% lead_snps)]
    lab_lead <- label_df[snp == lead_snp]
    lab_lead2 <- if (isTRUE(use_snp2)) label_df[snp == lead_snp2] else label_df[0]
    curved <- if (label.mode == "smart") 0.15 else 0
    repel_dir <- if (label.mode == "top") "y" else "both"
    if (nrow(lab_other)) {
      p <- p + ggrepel::geom_text_repel(data = lab_other, ggplot2::aes(x = xvalue, y = yval, label = label_text), nudge_x = lab_other$nudge_x, nudge_y = lab_other$nudge_y, direction = repel_dir, seed = label.seed, size = label.size, color = "grey15", box.padding = 0.35, point.padding = 0.25, min.segment.length = 0, segment.color = "grey45", segment.size = 0.35, segment.curvature = curved, segment.ncp = if (label.mode == "smart") 3 else 1, max.overlaps = Inf)
    }
    if (nrow(lab_lead)) {
      p <- p + ggrepel::geom_text_repel(data = lab_lead, ggplot2::aes(x = xvalue, y = yval, label = label_text), nudge_x = lab_lead$nudge_x, nudge_y = lab_lead$nudge_y, direction = repel_dir, seed = label.seed, size = label.size * 1.1, fontface = "bold", color = lead1_color, box.padding = 0.4, point.padding = 0.25, min.segment.length = 0, segment.color = lead1_color, segment.size = 0.45, segment.curvature = curved, segment.ncp = if (label.mode == "smart") 3 else 1, max.overlaps = Inf)
    }
    if (nrow(lab_lead2)) {
      p <- p + ggrepel::geom_text_repel(data = lab_lead2, ggplot2::aes(x = xvalue, y = yval, label = label_text), nudge_x = lab_lead2$nudge_x, nudge_y = lab_lead2$nudge_y, direction = repel_dir, seed = label.seed, size = label.size * 1.1, fontface = "bold", color = lead2_color, box.padding = 0.4, point.padding = 0.25, min.segment.length = 0, segment.color = lead2_color, segment.size = 0.45, segment.curvature = curved, segment.ncp = if (label.mode == "smart") 3 else 1, max.overlaps = Inf)
    }
  }
  p
}

.regional_attach_meta <- function(p, df0, use_snp2, lead_snps, beta.col, direction.a1, label.mode, label.col, label.size, label.seed, label.top.n, label_df, ld_meta, ld_rds0, ld_bfile0, build, build.gwas, chrom, pos.range, pos.unit, gene_meta) {
  if (isTRUE(use_snp2)) {
    ld_r2 <- df0[, .(SNP = snp, CHR = CHR, POS = POS, r2 = suppressWarnings(as.numeric(ld)), r2_1 = suppressWarnings(as.numeric(ld1)), r2_2 = suppressWarnings(as.numeric(ld2)))]
  } else {
    ld_r2 <- df0[, .(SNP = snp, CHR = CHR, POS = POS, r2 = suppressWarnings(as.numeric(ld)))]
  }
  ld_r2[, is_lead := (SNP %in% lead_snps)]
  direction_meta <- list(
    enabled = isTRUE(any(df0$point_shape != 21)),
    beta_col = beta.col,
    direction_a1_input = if (is.null(direction.a1)) NULL else class(direction.a1)[1],
    n_shape_up = as_int(sum(df0$point_shape == 24, na.rm = TRUE)),
    n_shape_down = as_int(sum(df0$point_shape == 25, na.rm = TRUE)),
    n_shape_circle = as_int(sum(df0$point_shape == 21, na.rm = TRUE))
  )
  label_meta <- list(
    mode = label.mode,
    label_col = label.col,
    label_size = as.numeric(label.size),
    label_seed = as_int(label.seed),
    label_top_n = as_int(label.top.n),
    n_requested = as_int(length(unique(as.character(label_df$snp)))),
    n_drawn = as_int(nrow(label_df))
  )
  ld_meta$ld_rds <- if (!is.null(ld_rds0) && nzchar(as.character(ld_rds0)[1])) abs_path(ld_rds0) else NA_character_
  ld_meta$ld_bfile <- if (!is.null(ld_bfile0) && nzchar(as.character(ld_bfile0)[1])) abs_path(ld_bfile0) else NA_character_
  ld_meta$ld_path <- if (identical(ld_meta$mode, "rds")) ld_meta$ld_rds else if (identical(ld_meta$mode, "bigsnpr")) ld_meta$ld_bfile else if (identical(ld_meta$mode, "matrix")) "ld.matrix" else NA_character_
  attr(p, "gcanvas_meta") <- list(
    build = list(plot = build, gwas = build.gwas),
    region = list(
      chromosome = chrom,
      position = c(start = as.numeric(pos.range[1]), end = as.numeric(pos.range[2])),
      pos_unit = pos.unit,
      n_variants = as_int(nrow(df0))
    ),
    ld = ld_meta,
    ld_r2 = ld_r2,
    direction = direction_meta,
    labels = label_meta,
    genes = gene_meta
  )
  p
}

#' Regional association plot
#'
#' Draws a regional Manhattan-style plot for one locus: variants colored by
#' LD r-squared to the lead variant, optional second-track overlay (e.g. eQTL),
#' direction-of-effect glyphs, credible-set highlighting, and gene/exon tracks
#' supplied by [regional.track()].
#'
#' The lead variant can be specified via `snp`, by `chrom` + `pos`, or
#' omitted entirely. When neither a SNP nor coordinates pin a lead variant
#' (or when the requested `snp` is not in `data`), the function falls back
#' to a no-lead mode: it draws the region defined by `pos.range` without a
#' highlighted lead point, lead-driven LD computation
#' (`ld.matrix` / `ld.rds` / `ld.bfile`) is skipped with a warning, and
#' point coloring uses precomputed `ld.col` / `ld2.col` columns if present.
#'
#' @param data A `data.frame`/`data.table` of summary statistics.
#' @param snp,snp2 Lead variant id(s) — `snp2` enables the second-track overlay.
#' @param chrom,pos Optional coordinates if `snp` is not provided.
#' @param snp.col,chrom.col,pos.col,p.col Column names in `data`.
#' @param flank Window size (bp) when only a lead variant is given.
#' @param pos.range Explicit `(start, end)` window override.
#' @param pos.unit X-axis units: `"mb"`, `"kb"`, or `"base"`.
#' @param alpha Point alpha.
#' @param y.col,x.title,y.title,title Optional axis / title overrides.
#' @param ld.col,ld2.col LD columns to look up.
#' @param beta.col,a1.col,a2.col,direction,direction.a1 Direction-of-effect
#'   options.
#' @param label.snp,label.top.n,label.col,label.size,label.seed,label.mode
#'   Variant-labeling controls.
#' @param y.max,y.breaks,y.rescale.at,y.rescale.ratio,y.rescale.breaks Y-axis
#'   range and break-and-rescale controls.
#' @param threshold,threshold.color,threshold.type,threshold.linewidth
#'   Significance-line controls.
#' @param plot.width,units Plot-size hints for label placement.
#' @param build,build.gwas Genome build of the GTF (and of `data` if different).
#' @param liftover.dir,liftover.chain Liftover controls when builds differ.
#' @param ld.matrix,ld.rds,ld.bfile,ld.cache.dir,threads LD reference inputs.
#' @param ld.color,snp.color,snp2.color Color scales for LD and points.
#' @param shared.mode,shared.color,shared.cutoff,shared.alpha Shared-signal
#'   overlay controls (when both `snp` and `snp2` are given).
#' @param legend.position Where the legend is placed.
#' @param gtf Optional path to a GTF / `gtf2rds` cache.
#' @param hl.gene,hl.gene.color,hl.gene.size Gene-highlight controls.
#' @param show.gene,show.gene.label,show.exon,show.omit_label Visibility toggles
#'   for the gene track.
#' @param gene.max_row,gene.max_n,gene.add,gene.force,gene.force.label,gene.force.ignore_max_row
#'   Gene-track layout controls (passed to [regional.track()]).
#' @param biotype.keep,biotype.priority,biotype.keep.all Biotype filters.
#' @param exon.size Exon-thickness in the gene track.
#' @param snp.color.missing,snp2.color.missing Logical. Treat missing LD as
#'   "color missing" instead of dropping the point.
#' @param biotype.keep.default Logical. Use the default biotype keep set.
#'
#' @return A `ggplot` object combining the association panel and gene tracks.
#' @export
regional <- function(data,
                     snp = NULL, snp2 = NULL, chrom = NULL, pos = NULL,
                     snp.col = "SNP", chrom.col = "CHR", pos.col = "POS", p.col = "P",
                     flank = 5e5, pos.range = NULL, pos.unit = c("mb", "kb", "base"),
                     alpha = 1,
                     y.col = NULL, x.title = NULL, y.title = NULL,
                     title = NULL,
                     ld.col = "ld", ld2.col = "ld2",
                     beta.col = NULL, a1.col = NULL, a2.col = NULL,
                     direction = FALSE, direction.a1 = NULL,
                     label.snp = NULL, label.top.n = 0L, label.col = "SNP",
                     label.size = 3, label.seed = 23L,
                     label.mode = c("smart", "top", "off"),
                     y.max = NULL, y.breaks = NULL,
                     y.rescale.at = NULL, y.rescale.ratio = 0.25, y.rescale.breaks = NULL,
                     threshold = 5e-08,
                     threshold.color = "grey20",
                     threshold.type = NULL,
                     threshold.linewidth = 0.7,
                     plot.width = NULL, units = "cm",
                     build = 38L,
                     build.gwas = NULL,
                     liftover.dir = NULL,
                     liftover.chain = NULL,
                     ld.matrix = NULL,
                     ld.rds = NULL,
                     ld.bfile = NULL,
                     ld.cache.dir = NULL,
                     threads = 4L,
                     ld.color = NULL,
                     snp.color = "#CC2936",
                     snp2.color = "#08415C",
                     shared.mode = c("off", "split", "single"),
                     shared.color = "#702B9E",
                     shared.cutoff = 0.6,
                     shared.alpha = 1,
                     legend.position = c("bottom", "top", "left", "right", "left.outside", "right.outside"),
                     gtf = NULL,
                     hl.gene = NULL,
                     hl.gene.color = NULL,
                     hl.gene.size = 3,
                     show.gene = TRUE,
                     show.gene.label = TRUE,
                     show.exon = TRUE,
                     show.omit_label = TRUE,
                     gene.max_row = 8L,
                     gene.max_n = 50L,
                     gene.add = NULL,
                     gene.force = character(),
                     gene.force.label = NULL,
                     gene.force.ignore_max_row = TRUE,
                     biotype.keep = c("protein_coding", "lncRNA"),
                     biotype.priority = c("protein_coding", "lncRNA"),
                     biotype.keep.all = FALSE,
                     exon.size = 1,
                     snp.color.missing = FALSE,
                     snp2.color.missing = FALSE,
                     biotype.keep.default = TRUE) {
  df <- data
  require_pkg(c("ggplot2", "RColorBrewer", "scales", "data.table"))

  resolve_options <- function() {
    pos.unit <- match.arg(tolower(as.character(pos.unit)[1]), choices = c("mb", "kb", "base"))
    legend.position <- match.arg(tolower(as.character(legend.position)[1]), choices = c("bottom", "top", "left", "right", "left.outside", "right.outside"))
    label.mode <- match.arg(tolower(as.character(label.mode)[1]), choices = c("smart", "top", "off"))
    label.top.n <- as_int(label.top.n)
    if (is.na(label.top.n) || label.top.n < 0L) label.top.n <- 0L
    label.size <- suppressWarnings(as.numeric(label.size))[1]
    if (!is.finite(label.size) || label.size <= 0) label.size <- 3
    label.seed <- as_int(label.seed)
    if (is.na(label.seed)) label.seed <- 42L
    alpha <- suppressWarnings(as.numeric(alpha))[1]
    if (!is.finite(alpha)) alpha <- 1
    alpha <- max(0, min(1, alpha))
    default_snp_color <- as.character(formals(regional)$snp.color)[1]
    default_snp2_color <- as.character(formals(regional)$snp2.color)[1]
    default_shared_color <- as.character(formals(regional)$shared.color)[1]
    shared_mode_in <- shared.mode[1]
    if (is.logical(shared_mode_in)) {
      shared_mode <- if (isTRUE(shared_mode_in)) "split" else "off"
    } else {
      sm <- tolower(as.character(shared_mode_in))
      if (sm %in% c("true", "t", "1", "yes", "y")) {
        shared_mode <- "split"
      } else if (sm %in% c("false", "f", "0", "no", "n")) {
        shared_mode <- "off"
      } else if (sm %in% c("off", "none", "single", "split")) {
        shared_mode <- if (sm == "none") "off" else sm
      } else {
        shared_mode <- "off"
      }
    }
    shared.cutoff <- suppressWarnings(as.numeric(shared.cutoff))[1]
    if (!is.finite(shared.cutoff)) shared.cutoff <- 0.2
    shared.cutoff <- max(0, min(1, shared.cutoff))
    shared.alpha <- suppressWarnings(as.numeric(shared.alpha))[1]
    if (!is.finite(shared.alpha)) shared.alpha <- 1
    shared.alpha <- max(0, min(1, shared.alpha))
    shared.color <- as.character(shared.color)[1]
    if (is.na(shared.color) || !nzchar(shared.color)) shared.color <- default_shared_color
    if (is.null(build.gwas) || length(build.gwas) == 0L || is.na(build.gwas)) build.gwas <- build
    ld_color0 <- as.character(ld.color)[1]
    if (is.null(ld.color) || length(ld.color) == 0L || is.na(ld_color0) || !nzchar(ld_color0) || identical(tolower(ld_color0), "auto")) {
      ldcolorset <- .gcanvas_ld_palette(100)
    } else {
      ldcolorset <- ld.color
    }
    list(
      pos.unit = pos.unit,
      legend.position = legend.position,
      label.mode = label.mode,
      label.top.n = label.top.n,
      label.size = label.size,
      label.seed = label.seed,
      alpha = alpha,
      default_snp_color = default_snp_color,
      default_snp2_color = default_snp2_color,
      default_shared_color = default_shared_color,
      shared_mode = shared_mode,
      shared.cutoff = shared.cutoff,
      shared.alpha = shared.alpha,
      shared.color = shared.color,
      build.gwas = build.gwas,
      ldcolorset = ldcolorset
    )
  }

  prepare_input <- function() {
    snp_col_use <- .gcanvas_resolve_colname(names(df), snp.col, aliases = character(), required = TRUE, arg_label = "snp.col")
    chrom_col_use <- .gcanvas_resolve_colname(names(df), chrom.col, aliases = character(), required = TRUE, arg_label = "chrom.col")
    pos_col_use <- .gcanvas_resolve_colname(names(df), pos.col, aliases = character(), required = TRUE, arg_label = "pos.col")
    cols <- c(snp_col_use, chrom_col_use, pos_col_use)
    if (is.null(y.col)) cols <- c(cols, p.col)
    if (!is.null(y.col)) cols <- c(cols, y.col)
    if (label.mode != "off" && !is.null(label.col) && nzchar(as.character(label.col)[1]) && as.character(label.col)[1] != "SNP") {
      cols <- c(cols, as.character(label.col)[1])
    }
    cols <- unique(cols)
    miss <- setdiff(cols, names(df))
    if (length(miss)) stop("Missing columns in data: ", paste(miss, collapse = ", "), call. = FALSE)
    if (data.table::is.data.table(df)) {
      df0 <- data.table::copy(df[, ..cols])
    } else {
      df0 <- data.table::as.data.table(df[, cols, drop = FALSE])
    }
    if (snp_col_use %in% names(df0)) data.table::setnames(df0, snp_col_use, "snp")
    if (chrom_col_use %in% names(df0)) data.table::setnames(df0, chrom_col_use, "CHR")
    if (pos_col_use %in% names(df0)) data.table::setnames(df0, pos_col_use, "POS")
    if (is.null(y.col)) {
      if (p.col %in% names(df0)) data.table::setnames(df0, p.col, "P")
    } else {
      if (y.col %in% names(df0)) data.table::setnames(df0, y.col, "Y")
    }
    allele_a1_col <- a1.col
    allele_a2_col <- a2.col
    if (!is.null(allele_a1_col) && !is.null(allele_a2_col)) {
      if (!(allele_a1_col %in% names(df)) || !(allele_a2_col %in% names(df))) {
        stop("a1.col/a2.col not found in df.", call. = FALSE)
      }
      df0$A1 <- toupper(as.character(df[[allele_a1_col]]))
      df0$A2 <- toupper(as.character(df[[allele_a2_col]]))
    }
    if (ld.col %in% names(df)) df0$ld_precomputed <- suppressWarnings(as.numeric(df[[ld.col]]))
    if (ld2.col %in% names(df)) df0$ld2_precomputed <- suppressWarnings(as.numeric(df[[ld2.col]]))
    if (!is.null(beta.col)) {
      if (!(beta.col %in% names(df))) stop("beta.col not found in df.", call. = FALSE)
      df0$beta <- suppressWarnings(as.numeric(df[[beta.col]]))
    }
    if (is.character(direction.a1) && length(direction.a1) == 1L && !is.na(direction.a1) && (direction.a1 %in% names(df))) {
      df0$direction_a1_input <- toupper(as.character(df[[direction.a1]]))
    }
    label_col_src <- as.character(label.col)[1]
    if (is.na(label_col_src) || !nzchar(label_col_src)) label_col_src <- "SNP"
    if (identical(label_col_src, "SNP") || identical(label_col_src, snp_col_use)) {
      label_col_src <- "snp"
    } else if (identical(label_col_src, chrom_col_use)) {
      label_col_src <- "CHR"
    } else if (identical(label_col_src, pos_col_use)) {
      label_col_src <- "POS"
    } else if (identical(label_col_src, p.col)) {
      label_col_src <- "P"
    } else if (!is.null(y.col) && identical(label_col_src, y.col)) {
      label_col_src <- "Y"
    }
    df0$snp <- as.character(df0$snp)
    df0$CHR <- normalize.chrom(df0$CHR)
    df0$POS <- suppressWarnings(as.numeric(df0$POS))
    if (is.null(y.col)) {
      df0$P <- .gcanvas_p_filter(df0$P)
    } else {
      df0$Y <- suppressWarnings(as.numeric(df0$Y))
    }
    if (is.null(y.col)) {
      df0 <- df0[!is.na(P) & nzchar(P) & is.finite(POS) & !is.na(CHR) & nzchar(CHR)]
    } else {
      df0 <- df0[is.finite(Y) & is.finite(POS) & !is.na(CHR) & nzchar(CHR)]
    }
    if (!nrow(df0)) stop("No valid rows after parsing.", call. = FALSE)
    flank_use <- .gcanvas_parse_bp_span(flank, arg_name = "flank")
    if (!is.finite(flank_use) || is.na(flank_use) || flank_use <= 0) flank_use <- 5e5
    df0$yval <- if (is.null(y.col)) -log10c(df0$P) else as.numeric(df0$Y)
    df0 <- df0[is.finite(yval) & !is.na(yval)]
    if (!nrow(df0)) stop("No valid rows after parsing.", call. = FALSE)
    y_rescale_at0 <- suppressWarnings(as.numeric(y.rescale.at))[1]
    if (!is.finite(y_rescale_at0) || is.na(y_rescale_at0) || y_rescale_at0 <= 0) y_rescale_at0 <- NA_real_
    y_rescale_ratio0 <- suppressWarnings(as.numeric(y.rescale.ratio))[1]
    if (!is.finite(y_rescale_ratio0) || is.na(y_rescale_ratio0) || y_rescale_ratio0 <= 0 || y_rescale_ratio0 > 1) y_rescale_ratio0 <- 0.25
    y_rescale_breaks0 <- suppressWarnings(as.numeric(y.rescale.breaks))
    y_rescale_breaks0 <- y_rescale_breaks0[is.finite(y_rescale_breaks0)]
    use_y_rescale <- is.finite(y_rescale_at0) && !is.na(y_rescale_at0)
    .manhattan_y_map_reg <- function(v) {
      vv <- suppressWarnings(as.numeric(v))
      if (!isTRUE(use_y_rescale)) return(vv)
      out <- vv
      idx <- is.finite(vv) & vv > y_rescale_at0
      out[idx] <- y_rescale_at0 + (vv[idx] - y_rescale_at0) * y_rescale_ratio0
      out
    }
    .compute_y_breaks_reg <- function(ylim_vec = NULL) {
      if (is.null(ylim_vec) || length(ylim_vec) < 2L || any(!is.finite(ylim_vec))) {
        ylim_vec <- range(df0$yval, na.rm = TRUE)
      }
      if (is.null(y.breaks) || (is.character(y.breaks) && length(y.breaks) == 1L && tolower(trimws(as.character(y.breaks)[1])) == "auto")) {
        if (isTRUE(use_y_rescale)) {
          raw_max <- max(df0$yval_raw, na.rm = TRUE)
          if (!is.finite(raw_max) || raw_max <= 0) raw_max <- y_rescale_at0
          br_low <- suppressWarnings(as.numeric(base::pretty(c(0, y_rescale_at0), n = 3)))
          br_low <- br_low[is.finite(br_low) & br_low <= y_rescale_at0]
          if (length(y_rescale_breaks0)) {
            br_high <- sort(unique(y_rescale_breaks0[y_rescale_breaks0 >= y_rescale_at0]))
          } else {
            br_high <- suppressWarnings(as.numeric(base::pretty(c(y_rescale_at0, raw_max), n = 2)))
            br_high <- br_high[is.finite(br_high) & br_high >= y_rescale_at0]
          }
          y_raw <- sort(unique(c(0, br_low, y_rescale_at0, br_high)))
          if (length(y_raw) > 6L) {
            lo <- y_raw[y_raw < y_rescale_at0]
            hi <- y_raw[y_raw > y_rescale_at0]
            lo_keep <- if (length(lo) > 2L) lo[unique(as_int(round(seq(1, length(lo), length.out = 2))))] else lo
            hi_keep <- if (length(hi) > 2L) hi[unique(as_int(round(seq(1, length(hi), length.out = 2))))] else hi
            y_raw <- sort(unique(c(0, lo_keep, y_rescale_at0, hi_keep)))
          }
          if (!length(y_raw)) y_raw <- c(0, y_rescale_at0)
          y_draw <- .manhattan_y_map_reg(y_raw)
          br <- data.table::data.table(draw = y_draw, lab = as.character(signif(y_raw, 4)))
          br <- br[is.finite(draw) & draw >= 0]
          if (nrow(br)) {
            data.table::setorder(br, draw)
            br[, draw_key := round(draw, 6)]
            br <- br[!duplicated(draw_key)]
            br[, draw_key := NULL]
            if (!any(abs(br$draw) < .Machine$double.eps)) {
              br <- data.table::rbindlist(list(data.table::data.table(draw = 0, lab = "0"), br), use.names = TRUE)
              br <- br[!duplicated(draw)]
            }
            y_at_draw <- .manhattan_y_map_reg(y_rescale_at0)
            min_gap_at <- max(1.0, diff(range(ylim_vec, na.rm = TRUE)) * 0.04)
            idx_at <- which.min(abs(br$draw - y_at_draw))
            if (length(idx_at) == 1L && idx_at > 1L) {
              gap_below <- y_at_draw - br$draw[idx_at - 1L]
              if (is.finite(gap_below) && gap_below < min_gap_at) br <- br[-(idx_at - 1L)]
            }
            return(list(draw = br$draw, labels = br$lab))
          }
          return(list(draw = 0, labels = "0"))
        }
        yb <- suppressWarnings(as.numeric(scales::pretty_breaks()(ylim_vec)))
        yb <- yb[is.finite(yb) & yb >= 0]
        if (!length(yb)) yb <- 0
        if (!any(abs(yb) < .Machine$double.eps)) yb <- c(0, yb)
        yb <- sort(unique(yb))
        return(list(draw = yb, labels = ggplot2::waiver()))
      }
      yb_raw <- suppressWarnings(as.numeric(y.breaks))
      yb_raw <- yb_raw[is.finite(yb_raw)]
      if (!length(yb_raw)) yb_raw <- 0
      if (isTRUE(use_y_rescale)) {
        yb_draw <- .manhattan_y_map_reg(yb_raw)
        br <- data.table::data.table(draw = yb_draw, lab = as.character(signif(yb_raw, 4)))
        br <- br[is.finite(draw) & draw >= 0]
        data.table::setorder(br, draw)
        br[, draw_key := round(draw, 6)]
        br <- br[!duplicated(draw_key)]
        br[, draw_key := NULL]
        if (!nrow(br)) return(list(draw = 0, labels = "0"))
        if (!any(abs(br$draw) < .Machine$double.eps)) {
          br <- data.table::rbindlist(list(data.table::data.table(draw = 0, lab = "0"), br), use.names = TRUE)
          br <- br[!duplicated(draw)]
        }
        return(list(draw = br$draw, labels = br$lab))
      }
      yb <- yb_raw[yb_raw >= 0]
      if (!length(yb)) yb <- 0
      if (!any(abs(yb) < .Machine$double.eps)) yb <- c(0, yb)
      yb <- sort(unique(yb))
      list(draw = yb, labels = ggplot2::waiver())
    }
    df0$yval_raw <- df0$yval
    df0$yval <- .manhattan_y_map_reg(df0$yval_raw)
    if (label_col_src %in% names(df0)) {
      df0$label_text <- as.character(df0[[label_col_src]])
    } else {
      df0$label_text <- as.character(df0$snp)
    }
    df0$label_text[is.na(df0$label_text)] <- ""
    df0$point_shape <- 21
    snp_is_auto <- is.null(snp) || (is.character(snp) && length(snp) == 1L && tolower(trimws(as.character(snp)[1])) == "auto")
    chrom_is_auto <- is.null(chrom) || (is.character(chrom) && length(chrom) == 1L && tolower(trimws(as.character(chrom)[1])) == "auto")
    pos_is_auto <- is.null(pos) || (is.character(pos) && length(pos) == 1L && tolower(trimws(as.character(pos)[1])) == "auto")
    posrange_given <- !is.null(pos.range) &&
      !(length(pos.range) == 1L && is.character(pos.range) && tolower(trimws(as.character(pos.range)[1])) == "auto")

    has_lead <- FALSE
    .resolve_chrom_from_data <- function() {
      uc <- unique(normalize.chrom(df0$CHR))
      uc <- uc[!is.na(uc) & nzchar(uc)]
      if (length(uc) != 1L) {
        stop("No lead SNP / chrom given. Either supply `chrom`, or filter `data` to a single chromosome.", call. = FALSE)
      }
      uc[1]
    }

    if (!snp_is_auto) {
      snp_query <- as.character(snp)[1]
      ii <- which(df0$snp == snp_query)
      if (length(ii)) {
        if (!chrom_is_auto || !pos_is_auto) {
          warning("snp matched in data; chrom/pos are ignored and lead coordinates are taken from the matched SNP row.", call. = FALSE)
        }
        snp <- snp_query
        chrom <- normalize.chrom(df0$CHR[ii[1]])[1]
        pos <- as.numeric(df0$POS[ii[1]])
        has_lead <- TRUE
      } else {
        # SNP not found in data: try chrom+pos fallback, else fall back to no-lead pos.range mode.
        .gcanvas_warn_msg(sprintf(
          "Lead SNP '%s' not found in data; proceeding without highlighting a lead variant.",
          snp_query
        ))
        if (!chrom_is_auto && !pos_is_auto) {
          chrom <- normalize.chrom(chrom)[1]
          pos <- suppressWarnings(as.numeric(pos))
          if (!is.finite(pos)) stop("Invalid pos.", call. = FALSE)
          snp <- NA_character_
        } else if (posrange_given) {
          chrom <- if (chrom_is_auto) .resolve_chrom_from_data() else normalize.chrom(chrom)[1]
          snp <- NA_character_
          pos <- NA_real_
        } else {
          stop("Lead SNP not found in data and no fallback (chrom+pos or pos.range) supplied.", call. = FALSE)
        }
      }
    } else if (!chrom_is_auto && !pos_is_auto) {
      chrom <- normalize.chrom(chrom)[1]
      pos <- suppressWarnings(as.numeric(pos))
      if (!is.finite(pos)) stop("Invalid pos.", call. = FALSE)
      ii <- which(df0$CHR == chrom & df0$POS == pos)
      if (length(ii)) {
        snp <- df0$snp[ii[1]]
        has_lead <- TRUE
      } else {
        .gcanvas_warn_msg(sprintf(
          "No SNP at chr%s:%s found in data; proceeding without highlighting a lead variant.",
          as.character(chrom), format(pos, scientific = FALSE)
        ))
        snp <- NA_character_
      }
    } else {
      # Neither snp nor (chrom + pos) given -- require pos.range to know the window.
      if (!posrange_given) {
        stop(
          "regional() needs at least one of: (1) snp; (2) chrom + pos; or (3) chrom + pos.range. None were provided.",
          call. = FALSE
        )
      }
      chrom <- if (chrom_is_auto) .resolve_chrom_from_data() else normalize.chrom(chrom)[1]
      snp <- NA_character_
      pos <- NA_real_
    }

    use_snp2 <- FALSE
    pos2 <- NA_real_
    if (has_lead && !is.null(snp2) && length(snp2) > 0L) {
      snp2_in <- as.character(snp2)[1]
      if (!is.na(snp2_in) && nzchar(snp2_in) && !(tolower(snp2_in) %in% c("false", "auto"))) {
        snp2 <- snp2_in
        if (!identical(snp2, snp)) {
          ii2 <- which(df0$snp == snp2)
          if (!length(ii2)) stop("Lead SNP2 not found in data.", call. = FALSE)
          chrom2 <- normalize.chrom(df0$CHR[ii2[1]])[1]
          pos2 <- as.numeric(df0$POS[ii2[1]])
          if (!identical(chrom2, chrom)) stop("snp2 must be on the same chromosome as snp.", call. = FALSE)
          use_snp2 <- TRUE
        }
      }
    } else if (!has_lead && !is.null(snp2) && length(snp2) > 0L) {
      snp2_in <- as.character(snp2)[1]
      if (!is.na(snp2_in) && nzchar(snp2_in) && !(tolower(snp2_in) %in% c("false", "auto"))) {
        .gcanvas_warn_msg("snp2 is ignored when no primary lead SNP is resolved.")
      }
    }

    if (is.null(pos.range) || (length(pos.range) == 1 && tolower(as.character(pos.range)) == "auto")) {
      if (!has_lead) {
        stop("pos.range is required when no lead SNP is resolved.", call. = FALSE)
      }
      if (isTRUE(use_snp2) && is.finite(pos2)) {
        pos.range <- c(min(pos, pos2) - flank_use, max(pos, pos2) + flank_use)
      } else {
        pos.range <- c(pos - flank_use, pos + flank_use)
      }
    }
    pos.range <- .gcanvas_as_num2(pos.range)
    if (length(pos.range) != 2 || any(!is.finite(pos.range))) stop("pos.range must be numeric length-2.", call. = FALSE)
    pos.range <- c(min(pos.range), max(pos.range))
    build <- as_int(build)
    build.gwas <- as_int(build.gwas)
    if (!is.na(build) && !is.na(build.gwas) && build != build.gwas) {
      has_liftover_dir <- !is.null(liftover.dir) && length(liftover.dir) > 0L && !is.na(liftover.dir[1]) && nzchar(as.character(liftover.dir[1]))
      if (!has_liftover_dir) {
        stop("build and build.gwas differ, but liftover.dir is NULL/empty. Provide liftover.dir or set build.gwas = build.", call. = FALSE)
      }
      liftover.dir <- as.character(liftover.dir)[1]
      .gcanvas_note("gcanvas::regioanl", sprintf("Liftover GWAS: b%s -> b%s", build.gwas, build))
      tmp <- data.frame(SNP = df0$snp, CHR = df0$CHR, POS = df0$POS, stringsAsFactors = FALSE)
      tmp2 <- liftover(tmp, from = build.gwas, to = build, liftover.dir = liftover.dir, liftover.chain = liftover.chain, SNP = "SNP", CHR = "CHR", POS = "POS", silent = TRUE)
      pos_new <- tmp2[[paste0("POS_b", build)]]
      drop_n <- sum(is.na(pos_new), na.rm = TRUE)
      if (drop_n) .gcanvas_warn_msg(sprintf("Dropped %d variants (liftover NA).", drop_n))
      df0$POS <- pos_new
      df0 <- df0[!is.na(POS)]
      if (!nrow(df0)) stop("No variants left after liftover.", call. = FALSE)
      if (has_lead) {
        pos <- liftover_positions(chrom, pos, from = build.gwas, to = build, liftover.dir = liftover.dir, liftover.chain = liftover.chain)
      }
      if (isTRUE(use_snp2)) pos2 <- liftover_positions(chrom, pos2, from = build.gwas, to = build, liftover.dir = liftover.dir, liftover.chain = liftover.chain)
      pos.range <- liftover_positions(chrom, pos.range, from = build.gwas, to = build, liftover.dir = liftover.dir, liftover.chain = liftover.chain)
      if (any(is.na(pos.range)) || (has_lead && is.na(pos)) || (isTRUE(use_snp2) && is.na(pos2))) stop("Lead/pos.range liftover failed.", call. = FALSE)
      pos.range <- c(min(pos.range), max(pos.range))
    }
    df0 <- df0[CHR == chrom & POS >= pos.range[1] & POS <= pos.range[2]]
    if (!nrow(df0)) stop("No variants in region after filtering.", call. = FALSE)
    min_p_msg <- "NA"
    if (is.null(y.col) && ("P" %in% names(df0))) {
      min_idx <- suppressWarnings(which.max(df0$yval))[1]
      if (is.finite(min_idx) && !is.na(min_idx) && min_idx >= 1L) {
        min_p_msg <- .gcanvas_format_minp(df0$P[min_idx], digits = 3L, cutoff = 1e-3)
      }
    }
    .gcanvas_note("gcanvas::regioanl", sprintf("n_variants=%d | chromosome=%s | minP=%s", as_int(nrow(df0)), as.character(chrom)[1], min_p_msg))
    direction <- isTRUE(direction)
    if (direction) {
      if (is.null(beta.col) || !("beta" %in% names(df0))) stop("direction=TRUE requires beta.col.", call. = FALSE)
      beta_aligned <- as.numeric(df0$beta)
      beta_aligned[!is.finite(beta_aligned)] <- NA_real_
      if (!is.null(direction.a1)) {
        dir_obj <- .gcanvas_resolve_direction_a1(direction.a1, df = df, snp.col = snp_col_use, a1.col = allele_a1_col %||% "A1")
        dir_ref <- rep(NA_character_, nrow(df0))
        if (identical(dir_obj$type, "column")) {
          dir_ref <- if ("direction_a1_input" %in% names(df0)) df0$direction_a1_input else rep(NA_character_, nrow(df0))
        } else if (identical(dir_obj$type, "map")) {
          dir_ref <- unname(dir_obj$map[df0$snp])
        }
        have_a1a2 <- ("A1" %in% names(df0)) && ("A2" %in% names(df0))
        if (!have_a1a2 && !identical(dir_obj$type, "column")) {
          stop("direction.a1 requires a1.col and a2.col unless direction.a1 is a df column name.", call. = FALSE)
        }
        if (have_a1a2) {
          m_a1 <- !is.na(dir_ref) & nzchar(dir_ref) & !is.na(df0$A1) & (dir_ref == df0$A1)
          m_a2 <- !is.na(dir_ref) & nzchar(dir_ref) & !is.na(df0$A2) & (dir_ref == df0$A2)
          beta_aligned[m_a2] <- -1 * beta_aligned[m_a2]
          aligned_ok <- m_a1 | m_a2
          beta_aligned[!aligned_ok] <- 0
        } else {
          has_ref <- !is.na(dir_ref) & nzchar(dir_ref)
          beta_aligned[!has_ref] <- 0
        }
      }
      sgn <- sign(beta_aligned)
      sgn[!is.finite(sgn)] <- 0
      df0$point_shape <- ifelse(sgn > 0, 24, ifelse(sgn < 0, 25, 21))
    }
    y.min <- min(0, df0$yval, na.rm = TRUE)
    y.max0 <- if (is.null(y.max) || (is.character(y.max) && length(y.max) == 1L && tolower(trimws(as.character(y.max)[1])) == "auto")) {
      ceiling(max(df0$yval, na.rm = TRUE) * 1.05)
    } else {
      .manhattan_y_map_reg(as.numeric(y.max))
    }
    y.breaks_obj <- .compute_y_breaks_reg(c(y.min, y.max0))
    y.axis.title <- if (!is.null(y.title)) {
      if (is.character(y.title)) as.character(y.title)[1] else y.title
    } else if (!is.null(y.col)) {
      y.col
    } else {
      bquote(paste(-log[10], " (", italic(P), ")"))
    }
    threshold_lines <- .gcanvas_threshold_lines(
      threshold = threshold,
      y.col = y.col,
      threshold.color = threshold.color,
      threshold.type = threshold.type,
      threshold.linewidth = threshold.linewidth,
      map_fun = .manhattan_y_map_reg,
      default_color = "grey20"
    )
    list(
      df0 = df0,
      snp = snp,
      snp2 = snp2,
      chrom = chrom,
      pos = pos,
      pos2 = pos2,
      use_snp2 = use_snp2,
      has_lead = has_lead,
      pos.range = pos.range,
      build = build,
      build.gwas = build.gwas,
      y.min = y.min,
      y.max0 = y.max0,
      y.breaks_draw = y.breaks_obj$draw,
      y.breaks_labels = y.breaks_obj$labels,
      y.axis.title = y.axis.title,
      threshold_lines = threshold_lines,
      .manhattan_y_map_reg = .manhattan_y_map_reg,
      .compute_y_breaks_reg = .compute_y_breaks_reg
    )
  }

  compute_ld <- function(df0, chrom, pos, pos2, pos.range, snp, snp2, use_snp2, has_lead, ldcolorset) {
    ld_mode <- "none"
    ld_matrix0 <- NULL
    ld_rds0 <- NULL
    ld_bfile0 <- NULL
    has_ld1_pre <- "ld_precomputed" %in% names(df0)
    has_ld2_pre <- "ld2_precomputed" %in% names(df0)
    has_valid_ld_matrix <- FALSE
    if (!is.null(ld.matrix)) {
      ld_matrix0 <- .gcanvas_ld_extract_named_matrix(ld.matrix, max_depth = 8L)
      has_valid_ld_matrix <- !is.null(ld_matrix0)
      if (!has_valid_ld_matrix) {
        .gcanvas_warn_msg("ld.matrix provided but no named square numeric LD matrix found; ignoring ld.matrix.")
        ld_matrix0 <- NULL
      }
    }
    if (!is.null(ld.rds) && length(ld.rds) > 0L && !is.na(ld.rds[1]) && nzchar(as.character(ld.rds[1]))) ld_rds0 <- abs_path(as.character(ld.rds)[1])
    if (!is.null(ld.bfile) && length(ld.bfile) > 0L && !is.na(ld.bfile[1]) && nzchar(as.character(ld.bfile[1]))) ld_bfile0 <- .gcanvas_normalize_bfile_prefix(ld.bfile)

    # Lead-required LD backends (matrix / rds / bfile) need a lead SNP to
    # compute r2 against. If there's no lead, fall through to precomputed
    # columns (ld.col / ld2.col) or "none".
    if (!isTRUE(has_lead)) {
      if (has_valid_ld_matrix || !is.null(ld_rds0) || !is.null(ld_bfile0)) {
        .gcanvas_warn_msg("LD reference (ld.matrix/ld.rds/ld.bfile) provided but no lead SNP resolved; falling back to ld.col / ld2.col columns if present.")
      }
      has_valid_ld_matrix <- FALSE
      ld_matrix0 <- NULL
      ld_rds0 <- NULL
      ld_bfile0 <- NULL
      if (has_ld1_pre || (isTRUE(use_snp2) && has_ld2_pre)) {
        ld_mode <- "precomputed"
      } else {
        ld_mode <- "none"
      }
    } else if (has_valid_ld_matrix) {
      if (!is.null(ld_rds0) || !is.null(ld_bfile0) || has_ld1_pre || (isTRUE(use_snp2) && has_ld2_pre)) {
        .gcanvas_warn_msg("ld.matrix provided; using ld.matrix and ignoring other LD sources.")
      }
      ld_mode <- "matrix"
    } else if (!is.null(ld_rds0) && file.exists(ld_rds0)) {
      if (!is.null(ld_bfile0)) .gcanvas_warn_msg("Both ld.rds and ld.bfile provided; using ld.rds.")
      ld_mode <- "rds"
    } else if (!is.null(ld_rds0) && !file.exists(ld_rds0)) {
      .gcanvas_warn_msg(paste0("ld.rds not found: ", ld_rds0, " -> trying remaining LD sources"))
      if (!is.null(ld_bfile0)) {
        if (.gcanvas_ld_bfile_ready(ld_bfile0)) ld_mode <- "bigsnpr" else {
          .gcanvas_warn_msg(paste0("ld.bfile found but incomplete PLINK set (.bed/.bim/.fam): ", ld_bfile0, " -> using LD none"))
          ld_mode <- "none"
        }
      } else if (has_ld1_pre || (isTRUE(use_snp2) && has_ld2_pre)) {
        ld_mode <- "precomputed"
      }
    } else if (!is.null(ld_bfile0)) {
      if (.gcanvas_ld_bfile_ready(ld_bfile0)) ld_mode <- "bigsnpr" else {
        .gcanvas_warn_msg(paste0("ld.bfile found but incomplete PLINK set (.bed/.bim/.fam): ", ld_bfile0, " -> using LD none"))
        ld_mode <- "none"
      }
    } else if (has_ld1_pre || (isTRUE(use_snp2) && has_ld2_pre)) {
      ld_mode <- "precomputed"
    }
    ld_meta <- list(mode = ld_mode)
    if (ld_mode == "none") {
      df0$ld1 <- 1
      df0$ld2 <- if (isTRUE(use_snp2)) 1 else NA_real_
      df0$ld <- 1
      ldcolorset <- c("grey40", "grey40")
      if (!is.null(ld_bfile0) && !.gcanvas_ld_bfile_ready(ld_bfile0)) {
        ld_meta$reason <- "invalid_ld_bfile"
      } else if (!is.null(ld.matrix) && is.null(ld_matrix0)) {
        ld_meta$reason <- "invalid_ld_matrix"
      } else if (!is.null(ld_rds0) && !file.exists(ld_rds0)) {
        ld_meta$reason <- "invalid_ld_rds"
      } else {
        ld_meta$reason <- "no_ld_source"
      }
    } else if (ld_mode == "matrix") {
      r2_res <- ld_r2_from_matrix(ld_matrix = ld_matrix0, snp_ids = df0$snp, lead_snp = snp, source_label = "ld.matrix")
      df0$ld1 <- r2_res$r2
      if (isTRUE(use_snp2)) {
        r2_res2 <- ld_r2_from_matrix(ld_matrix = ld_matrix0, snp_ids = df0$snp, lead_snp = snp2, source_label = "ld.matrix")
        df0$ld2 <- r2_res2$r2
        ldmax <- pmax(df0$ld1, df0$ld2, na.rm = TRUE)
        ldmax[is.na(df0$ld1) & is.na(df0$ld2)] <- NA_real_
        df0$ld <- ldmax
        ld_meta$secondary <- c(list(snp2 = snp2), r2_res2$meta)
      } else {
        df0$ld2 <- NA_real_
        df0$ld <- df0$ld1
      }
      df0$ld[df0$snp == snp] <- 1
      if (isTRUE(use_snp2)) df0$ld[df0$snp == snp2] <- 1
      ld_meta <- c(ld_meta, r2_res$meta, list(ld_matrix = "in_memory"))
    } else if (ld_mode == "rds") {
      r2_res <- ld_r2_from_rds(ld_rds = ld_rds0, snp_ids = df0$snp, lead_snp = snp)
      df0$ld1 <- r2_res$r2
      if (isTRUE(use_snp2)) {
        r2_res2 <- ld_r2_from_rds(ld_rds = ld_rds0, snp_ids = df0$snp, lead_snp = snp2)
        df0$ld2 <- r2_res2$r2
        ldmax <- pmax(df0$ld1, df0$ld2, na.rm = TRUE)
        ldmax[is.na(df0$ld1) & is.na(df0$ld2)] <- NA_real_
        df0$ld <- ldmax
        ld_meta$secondary <- c(list(snp2 = snp2), r2_res2$meta)
      } else {
        df0$ld2 <- NA_real_
        df0$ld <- df0$ld1
      }
      df0$ld[df0$snp == snp] <- 1
      if (isTRUE(use_snp2)) df0$ld[df0$snp == snp2] <- 1
      ld_meta <- c(ld_meta, r2_res$meta, list(ld_rds = ld_rds0))
    } else if (ld_mode == "precomputed") {
      ld1 <- if (has_ld1_pre) suppressWarnings(as.numeric(df0$ld_precomputed)) else rep(NA_real_, nrow(df0))
      ld2 <- if (isTRUE(use_snp2) && has_ld2_pre) suppressWarnings(as.numeric(df0$ld2_precomputed)) else rep(NA_real_, nrow(df0))
      df0$ld1 <- ld1
      df0$ld2 <- if (isTRUE(use_snp2)) ld2 else NA_real_
      if (isTRUE(use_snp2)) {
        if (!has_ld1_pre) .gcanvas_warn_msg("ld.col not found; using ld2.col only for LD coloring.")
        if (!has_ld2_pre) .gcanvas_warn_msg("ld2.col not found; using ld.col only for LD coloring.")
        ldmax <- pmax(ld1, ld2, na.rm = TRUE)
        ldmax[is.na(ld1) & is.na(ld2)] <- NA_real_
        df0$ld <- ldmax
      } else {
        df0$ld <- ld1
      }
      ld_meta$ld_col <- ld.col
      ld_meta$ld2_col <- if (isTRUE(use_snp2)) ld2.col else NULL
      ld_meta$n_missing <- as_int(sum(!is.finite(df0$ld) | is.na(df0$ld)))
    } else {
      cache_dir <- as.character(ld.cache.dir)[1]
      if (is.null(ld.cache.dir) || length(ld.cache.dir) == 0L || is.na(cache_dir) || !nzchar(cache_dir) || tolower(cache_dir) == "auto") {
        cache_dir <- .gcanvas_default_ld_cache_dir(ld_bfile0)
      }
      cache_dir <- abs_path(cache_dir)
      dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
      .gcanvas_note("gcanvas::regioanl", paste0("LD mode=plink (bigsnpr), cache: ", cache_dir))
      ld_obj <- attach_ld_ref(bfile = ld_bfile0, ld_rds = NULL, cache_dir = cache_dir)
      ref_dt <- ref_bim_table(ld_obj)
      r2_res <- ld_r2_bigsnpr(ld_obj = ld_obj, ref_dt = ref_dt, chr = chrom, start = pos.range[1], end = pos.range[2], df_pos = df0$POS, lead_pos = pos, df_ea = NULL, df_nea = NULL, threads = as_int(threads))
      df0$ld1 <- r2_res$r2
      if (isTRUE(use_snp2)) {
        r2_res2 <- ld_r2_bigsnpr(ld_obj = ld_obj, ref_dt = ref_dt, chr = chrom, start = pos.range[1], end = pos.range[2], df_pos = df0$POS, lead_pos = pos2, df_ea = NULL, df_nea = NULL, threads = as_int(threads))
        df0$ld2 <- r2_res2$r2
        ldmax <- pmax(df0$ld1, df0$ld2, na.rm = TRUE)
        ldmax[is.na(df0$ld1) & is.na(df0$ld2)] <- NA_real_
        df0$ld <- ldmax
        ld_meta$secondary <- c(list(snp2 = snp2), r2_res2$meta)
      } else {
        df0$ld2 <- NA_real_
        df0$ld <- df0$ld1
      }
      df0$ld[df0$snp == snp] <- 1
      if (isTRUE(use_snp2)) df0$ld[df0$snp == snp2] <- 1
      ld_meta <- c(ld_meta, r2_res$meta, list(cache_dir = cache_dir, ld_bfile = ld_bfile0))
    }
    ld_meta$lead_snp <- snp
    ld_meta$lead_snp2 <- if (isTRUE(use_snp2)) snp2 else NULL
    list(df0 = df0, ld_meta = ld_meta, ldcolorset = ldcolorset, ld_rds0 = ld_rds0, ld_bfile0 = ld_bfile0)
  }

  prepare_tracks <- function() {
    tracks <- NULL
    gene_meta <- NULL
    use_gtf <- !is.null(gtf) && length(gtf) && !is.na(gtf[1]) && nzchar(as.character(gtf[1]))
    if (isTRUE(show.gene)) {
      biotype.keep.use <- biotype.keep
      include_other_biotypes <- isTRUE(biotype.keep.all)
      if (!is.null(biotype.keep.use)) {
        bk <- unique(tolower(as.character(biotype.keep.use)))
        bk <- bk[!is.na(bk) & nzchar(bk)]
        if ("all" %in% bk) {
          biotype.keep.use <- NULL
          include_other_biotypes <- TRUE
        }
      }
      lncrna_symbol_only_use <- isTRUE(biotype.keep.default) &&
        !isTRUE(include_other_biotypes) &&
        !is.null(biotype.keep.use) &&
        any(tolower(as.character(biotype.keep.use)) %in% "lncrna", na.rm = TRUE)
      track_src_bgz <- NULL
      track_src_obj <- NULL
      if (isTRUE(use_gtf)) {
        track_src_bgz <- if (grepl("\\.bgz$", gtf, ignore.case = TRUE)) gtf else gtf_prepare_tabix(gtf, sort = "auto", chr_order = "natural")
      } else {
        # No GTF supplied: fall back to the package-bundled annotation tracks
        # (tracks.b37 / tracks.b38), selected by `build`.
        track_src_obj <- tryCatch(
          .gcanvas_bundled_tracks(build = build),
          error = function(e) {
            .gcanvas_note("gcanvas::regional", sprintf("gene track skipped: %s", conditionMessage(e)))
            NULL
          }
        )
      }
      lead_pos_for_tracks <- if (isTRUE(has_lead) && is.finite(pos)) pos else mean(pos.range)
      if (!is.null(track_src_bgz) || !is.null(track_src_obj)) {
      tracks <- regional.track(
        gtf_bgz = track_src_bgz,
        tracks = track_src_obj,
        chrom = chrom,
        pos.range = pos.range,
        y.max = y.max0,
        lead_pos = lead_pos_for_tracks,
        keep_biotype = biotype.keep.use,
        gene_max_levels = gene.max_row,
        gene_max_n = gene.max_n,
        gene_add = gene.add,
        gene_force = gene.force,
        gene_force_label = gene.force.label,
        gene_force_ignore_max_levels = gene.force.ignore_max_row,
        biotype_priority = biotype.priority,
        lncrna_symbol_only = lncrna_symbol_only_use,
        include_other_biotypes = include_other_biotypes,
        plot.width = plot.width,
        units = units
      )
      gene_meta <- tracks$meta
      }
    }
    hl_map <- .gcanvas_normalize_highlight_map(hl.gene, hl.gene.color)
    if (!is.null(hl_map) && isTRUE(show.gene) && !is.null(tracks) && nrow(tracks$gene)) {
      tracks$gene$hl_col <- mapply(.gcanvas_pick_hl_color, MoreArgs = list(map = hl_map), gene_name = tracks$gene$gene_name, gene_id = tracks$gene$gene_id, USE.NAMES = FALSE)
      tracks$gene$hl <- !is.na(tracks$gene$hl_col)
      if (!is.null(tracks$exon) && nrow(tracks$exon)) {
        tracks$exon$hl_col <- mapply(.gcanvas_pick_hl_color, MoreArgs = list(map = hl_map), gene_name = tracks$exon$gene_name, gene_id = tracks$exon$gene_id, USE.NAMES = FALSE)
        tracks$exon$hl <- !is.na(tracks$exon$hl_col)
      }
    }
    if (!is.null(gene_meta) && !is.null(gene_meta$n_dropped_by_biotype_keep) && isTRUE(gene_meta$n_dropped_by_biotype_keep > 0)) {
      cnt <- gene_meta$biotype_filtered_counts
      if (!is.null(cnt) && length(cnt)) {
        kv <- paste0(names(cnt), "=", as.integer(cnt))
        .gcanvas_note("gcanvas::regioanl", sprintf("Biotype filtered genes: total=%d (%s)", as.integer(gene_meta$n_dropped_by_biotype_keep), paste(kv, collapse = ", ")))
      } else {
        .gcanvas_note("gcanvas::regioanl", sprintf("Biotype filtered genes: total=%d", as.integer(gene_meta$n_dropped_by_biotype_keep)))
      }
    }
    list(tracks = tracks, gene_meta = gene_meta, hl_map = hl_map)
  }

  opt <- resolve_options()
  list2env(opt, environment())
  prep <- prepare_input()
  list2env(prep, environment())
  ld_res <- compute_ld(df0 = df0, chrom = chrom, pos = pos, pos2 = pos2, pos.range = pos.range, snp = snp, snp2 = snp2, use_snp2 = use_snp2, has_lead = has_lead, ldcolorset = ldcolorset)
  df0 <- ld_res$df0
  ld_meta <- ld_res$ld_meta
  ldcolorset <- ld_res$ldcolorset
  ld_rds0 <- ld_res$ld_rds0
  ld_bfile0 <- ld_res$ld_bfile0
  trk <- prepare_tracks()
  tracks <- trk$tracks
  gene_meta <- trk$gene_meta
  hl_map <- trk$hl_map

  if (pos.unit == "mb") {
    reduc <- 1e06
    xlab_suffix <- " (Mb)"
  } else if (pos.unit == "kb") {
    reduc <- 1e03
    xlab_suffix <- " (kb)"
  } else {
    reduc <- 1
    xlab_suffix <- ""
  }
  df0$xvalue <- df0$POS / reduc
  xbound <- pos.range / reduc
  x.axis.title <- {
    x_auto <- is.null(x.title) ||
      (is.character(x.title) &&
        (length(x.title) == 0L || is.na(as.character(x.title)[1]) || !nzchar(as.character(x.title)[1]) ||
          tolower(trimws(as.character(x.title)[1])) == "auto"))
    if (isTRUE(x_auto)) paste0("Position on chr", chrom, xlab_suffix) else if (is.character(x.title)) as.character(x.title)[1] else x.title
  }
  lead_snp <- snp
  lead_snp2 <- if (isTRUE(use_snp2)) snp2 else NULL
  lead_snps <- if (isTRUE(use_snp2)) c(lead_snp, lead_snp2) else lead_snp
  lead1_color <- as.character(snp.color)[1]
  if (!isTRUE(use_snp2) && snp.color.missing) {
    lead1_color <- tail(ldcolorset, 1)
  } else if (is.na(lead1_color) || !nzchar(lead1_color) || tolower(lead1_color) == "auto") {
    lead1_color <- if (isTRUE(use_snp2)) default_snp_color else tail(ldcolorset, 1)
  }
  lead2_color <- as.character(snp2.color)[1]
  if (snp2.color.missing) {
    lead2_color <- default_snp2_color
  } else if (is.na(lead2_color) || !nzchar(lead2_color) || tolower(lead2_color) == "auto") {
    lead2_color <- default_snp2_color
  }
  lead_shape <- if (direction) {
    shp <- df0$point_shape[df0$snp == snp]
    if (length(shp) && is.finite(shp[1])) shp[1] else 21
  } else {
    23
  }
  lead_shape2 <- if (isTRUE(use_snp2) && direction) {
    shp2 <- df0$point_shape[df0$snp == snp2]
    if (length(shp2) && is.finite(shp2[1])) shp2[1] else 21
  } else {
    23
  }
  label_df <- df0[0L, .(snp, xvalue, yval, label_text)]
  label_max_y <- NA_real_
  if (label.mode != "off") {
    label_idx <- integer(0)
    if (label.top.n > 0L) {
      ord <- order(df0$yval, decreasing = TRUE, na.last = NA)
      label_idx <- unique(c(label_idx, ord[seq_len(min(label.top.n, length(ord)))]))
    }
    if (!is.null(label.snp) && length(label.snp)) {
      want <- unique(as.character(label.snp))
      label_idx <- unique(c(label_idx, which(df0$snp %in% want)))
    }
    if (length(label_idx)) {
      label_df <- df0[label_idx, .(snp, xvalue, yval, label_text)]
      label_df <- label_df[!is.na(label_text) & nzchar(label_text)]
    }
    if (nrow(label_df)) {
      require_pkg("ggrepel")
      data.table::setorder(label_df, -yval)
      y_span <- diff(range(df0$yval, na.rm = TRUE))
      if (!is.finite(y_span) || y_span <= 0) y_span <- 1
      x_span <- diff(xbound)
      if (!is.finite(x_span) || x_span <= 0) x_span <- 1
      step_y <- y_span * 0.06
      if (!is.finite(step_y) || step_y <= 0) step_y <- max(0.5, abs(max(df0$yval, na.rm = TRUE)) * 0.03)
      n_lab <- nrow(label_df)
      rank_desc <- seq_len(n_lab)
      label_df$target_y <- max(df0$yval, na.rm = TRUE) + (n_lab - rank_desc + 1) * step_y * 0.25
      label_df$nudge_x <- 0
      label_df$nudge_y <- pmax(step_y * 0.4, label_df$target_y - label_df$yval)
      if (label.mode == "smart") {
        dx <- x_span * 0.05
        dy <- step_y * 2
        side_step <- x_span * 0.06
        if (!is.finite(side_step) || side_step <= 0) side_step <- 0.2
        nx <- numeric(n_lab)
        ny <- numeric(n_lab)
        for (i in seq_len(n_lab)) {
          xi <- label_df$xvalue[i]
          yi <- label_df$yval[i]
          left_space <- xi - xbound[1]
          right_space <- xbound[2] - xi
          left_ct <- sum(df0$xvalue >= (xi - dx) & df0$xvalue < xi & abs(df0$yval - yi) <= dy, na.rm = TRUE)
          right_ct <- sum(df0$xvalue <= (xi + dx) & df0$xvalue > xi & abs(df0$yval - yi) <= dy, na.rm = TRUE)
          side <- if (right_space >= left_space) 1 else -1
          pref_ct <- if (side > 0) right_ct else left_ct
          alt_ct <- if (side > 0) left_ct else right_ct
          if (pref_ct > 8 && alt_ct + 2 < pref_ct) side <- -side
          both_crowded <- (left_ct > 10 && right_ct > 10)
          nx[i] <- if (both_crowded) 0 else side * side_step
          rank_boost <- (n_lab - i + 1) * step_y * 0.08
          ny[i] <- if (both_crowded) {
            max(step_y * 0.9, label_df$target_y[i] - yi) + rank_boost
          } else {
            max(step_y * 0.45, label_df$target_y[i] - yi) + rank_boost
          }
        }
        label_df$nudge_x <- nx
        label_df$nudge_y <- ny
      }
      label_max_y <- max(label_df$yval + label_df$nudge_y, na.rm = TRUE)
    }
  }
  ld_is_horiz <- legend.position %in% c("bottom", "top")
  ld_barwidth <- if (ld_is_horiz) grid::unit(5.5, "cm") else grid::unit(0.45, "cm")
  ld_barheight <- if (ld_is_horiz) grid::unit(0.45, "cm") else grid::unit(3, "cm")
  p <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, color = "grey20", linewidth = 0.5) +
    ggplot2::scale_x_continuous(name = x.axis.title, breaks = scales::pretty_breaks(), expand = c(0, 0)) +
    ggplot2::scale_shape_identity() +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(color = "grey20"),
      axis.text = ggplot2::element_text(size = 14),
      axis.title = ggplot2::element_text(size = 16),
      plot.title = ggplot2::element_text(size = 18),
      legend.title = ggplot2::element_text(face = "bold", size = 12),
      legend.text = ggplot2::element_text(size = 10)
    )
  p <- .regional_draw_ld_points(p = p, df0 = df0, use_snp2 = use_snp2, legend.position = legend.position, lead_snp = lead_snp, lead_snp2 = lead_snp2, lead1_color = lead1_color, lead2_color = lead2_color, lead_shape = lead_shape, lead_shape2 = lead_shape2, ldcolorset = ldcolorset, alpha = alpha, shared_mode = shared_mode, shared.color = shared.color, shared.cutoff = shared.cutoff, shared.alpha = shared.alpha)
  p <- .regional_apply_legend_theme(p = p, threshold_lines = threshold_lines, legend.position = legend.position, use_snp2 = use_snp2, title = title)
  track_draw <- .regional_draw_tracks(
    p = p, show.gene = show.gene, tracks = tracks, hl_map = hl_map, show.exon = show.exon, show.gene.label = show.gene.label,
    show.omit_label = show.omit_label, gene_meta = gene_meta, reduc = reduc, xbound = xbound, y.min = y.min, y.max0 = y.max0,
    label_max_y = label_max_y, hl.gene.size = hl.gene.size, .compute_y_breaks_reg = .compute_y_breaks_reg, exon.size = exon.size
  )
  p <- track_draw$p
  tracks <- track_draw$tracks
  ylims <- track_draw$ylims
  if (!is.null(track_draw$y.breaks_draw)) {
    y.breaks_draw <- track_draw$y.breaks_draw
    y.breaks_labels <- track_draw$y.breaks_labels
  }
  p <- p + ggplot2::scale_y_continuous(name = y.axis.title, breaks = y.breaks_draw, labels = y.breaks_labels)
  p <- .regional_draw_labels(p = p, label.mode = label.mode, label_df = label_df, lead_snps = lead_snps, lead_snp = lead_snp, lead_snp2 = lead_snp2, use_snp2 = use_snp2, label.seed = label.seed, label.size = label.size, lead1_color = lead1_color, lead2_color = lead2_color)
  if (is.null(ylims)) {
    p <- p + ggplot2::coord_cartesian(xlim = xbound)
  } else {
    p <- p + ggplot2::coord_cartesian(xlim = xbound, ylim = ylims)
  }
  p <- .regional_attach_meta(
    p = p, df0 = df0, use_snp2 = use_snp2, lead_snps = lead_snps, beta.col = beta.col, direction.a1 = direction.a1,
    label.mode = label.mode, label.col = label.col, label.size = label.size, label.seed = label.seed, label.top.n = label.top.n,
    label_df = label_df, ld_meta = ld_meta, ld_rds0 = ld_rds0, ld_bfile0 = ld_bfile0, build = build, build.gwas = build.gwas,
    chrom = chrom, pos.range = pos.range, pos.unit = pos.unit, gene_meta = gene_meta
  )
  attr(p, "gcanvas_tracks") <- tracks
  p
}

.gcanvas_map_values_by_name <- function(keys, values, default = NA) {
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

.manhattan_map_lead_values <- function(snps, groups, values, default = NA) {
  snps <- as.character(snps)
  groups <- as.character(groups)
  out <- rep_len(default, length(snps))
  if (is.null(values) || length(values) == 0L) return(out)
  if (is.list(values) && !is.data.frame(values) && !data.table::is.data.table(values)) {
    nms <- names(values)
    has_names <- !is.null(nms) && any(!is.na(nms) & nzchar(nms))
    if (has_names) {
      for (i in seq_along(values)) {
        gi <- as.character(nms[i])
        idx <- which(!is.na(groups) & groups == gi)
        if (!length(idx)) next
        vi <- values[[i]]
        if (is.null(vi) || length(vi) == 0L) next
        out[idx] <- .gcanvas_map_values_by_name(snps[idx], vi, default = out[idx])
      }
      return(out)
    }
  }
  .gcanvas_map_values_by_name(snps, values, default = out)
}

.manhattan_map_lead_numeric <- function(snps, groups, values, default = NA_real_) {
  out <- suppressWarnings(as.numeric(.manhattan_map_lead_values(snps, groups, values, default = default)))
  out[!is.finite(out)] <- suppressWarnings(as.numeric(rep_len(default, length(out)))[!is.finite(out)])
  out
}

.manhattan_map_lead_sizes <- function(snps,
                                      groups,
                                      values,
                                      default_main = 1.6,
                                      default_ceiling = 2.2) {
  lead_size_main <- .manhattan_map_lead_numeric(snps, groups, values, default = default_main)
  lead_size_main[!is.finite(lead_size_main) | is.na(lead_size_main) | lead_size_main <= 0] <- default_main
  lead_size_ceiling <- lead_size_main * 1.35
  bad_ceiling <- !is.finite(lead_size_ceiling) | is.na(lead_size_ceiling) | lead_size_ceiling <= 0
  lead_size_ceiling[bad_ceiling] <- default_ceiling
  list(
    lead_size_main = lead_size_main,
    lead_size_ceiling = lead_size_ceiling
  )
}

.manhattan_normalize_lead_table <- function(lead,
                                          lead.color = "#E63946",
                                          lead.size = NULL,
                                          lead.stroke = 1,
                                          lead.flank.color = NULL,
                                          lead.label = FALSE,
                                          lead.label.color = "grey20") {
  require_pkg("data.table")
  empty <- data.table::data.table(
    lead_id = integer(),
    lead_group = character(),
    snp = character(),
    CHR = character(),
    POS = numeric(),
    lead_color = character(),
    lead_size_main = numeric(),
    lead_size_ceiling = numeric(),
    lead_stroke = numeric(),
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

  .manhattan_parse_lead_dt <- function(dt, default_group = NA_character_) {
    dt <- data.table::as.data.table(dt)
    nms <- names(dt)
    nms_low <- tolower(nms)
    pick_col <- function(cands) {
      hit <- match(cands, nms_low, nomatch = 0L)
      hit <- hit[hit > 0L]
      if (!length(hit)) return(NA_character_)
      nms[hit[1]]
    }
    snp_col <- pick_col(c("snp", "snpid", "rsid", "variant", "variant_id"))
    chr_col <- pick_col(c("chr", "chrom", "chromosome"))
    pos_col <- pick_col(c("pos", "bp", "position"))
    color_col <- pick_col(c("lead_color", "color"))
    size_col <- pick_col(c("lead_size", "size"))
    stroke_col <- pick_col(c("lead_stroke", "stroke"))
    flank_color_col <- pick_col(c("lead_flank_color", "flank_color", "region_color"))
    label_col <- pick_col(c("lead_label", "label", "text"))
    label_color_col <- pick_col(c("lead_label_color", "label_color", "text_color"))
    group_col <- pick_col(c("lead_group", "group", "set", "cluster"))

    out0 <- data.table::data.table(
      lead_group = if (!is.na(group_col)) as.character(dt[[group_col]]) else as.character(default_group),
      snp = if (!is.na(snp_col)) as.character(dt[[snp_col]]) else NA_character_,
      CHR = if (!is.na(chr_col)) normalize.chrom(dt[[chr_col]]) else NA_character_,
      POS = if (!is.na(pos_col)) suppressWarnings(as.numeric(dt[[pos_col]])) else NA_real_,
      lead_color = if (!is.na(color_col)) as.character(dt[[color_col]]) else NA_character_,
      lead_size = if (!is.na(size_col)) suppressWarnings(as.numeric(dt[[size_col]])) else NA_real_,
      lead_stroke = if (!is.na(stroke_col)) suppressWarnings(as.numeric(dt[[stroke_col]])) else NA_real_,
      flank_color = if (!is.na(flank_color_col)) as.character(dt[[flank_color_col]]) else NA_character_,
      lead_label = if (!is.na(label_col)) as.character(dt[[label_col]]) else NA_character_,
      lead_label_color = if (!is.na(label_color_col)) as.character(dt[[label_color_col]]) else NA_character_
    )
    out0
  }

  if (is.data.frame(lead) || data.table::is.data.table(lead)) {
    dt <- if (data.table::is.data.table(lead)) data.table::copy(lead) else data.table::as.data.table(lead)
    out <- .manhattan_parse_lead_dt(dt)
  } else if (is.list(lead) && !is.data.frame(lead) && !data.table::is.data.table(lead) &&
             !is.null(names(lead)) && any(!is.na(names(lead)) & nzchar(names(lead)))) {
    lead_nm <- names(lead)
    parts <- lapply(seq_along(lead), function(i) {
      li <- lead[[i]]
      gi <- as.character(lead_nm[i])
      if (is.na(gi) || !nzchar(gi)) gi <- sprintf("group%d", i)
      if (is.null(li) || length(li) == 0L) return(NULL)
      if (is.character(li) && length(li) == 1L && !is.na(li[1]) && nzchar(li[1]) && file.exists(li[1])) {
        li <- data.table::fread(li[1], data.table = TRUE, showProgress = FALSE)
      }
      if (is.data.frame(li) || data.table::is.data.table(li)) {
        dt0 <- if (data.table::is.data.table(li)) data.table::copy(li) else data.table::as.data.table(li)
        return(.manhattan_parse_lead_dt(dt0, default_group = gi))
      }
      snps <- .gcanvas_as_snp_vector(li)
      if (!length(snps)) return(NULL)
      data.table::data.table(
        lead_group = gi,
        snp = snps,
        CHR = NA_character_,
        POS = NA_real_,
        lead_color = NA_character_,
        lead_size = NA_real_,
        lead_stroke = NA_real_,
        flank_color = NA_character_,
        lead_label = NA_character_,
        lead_label_color = NA_character_
      )
    })
    parts <- Filter(Negate(is.null), parts)
    if (!length(parts)) return(empty)
    out <- data.table::rbindlist(parts, use.names = TRUE, fill = TRUE)
  } else {
    snps <- .gcanvas_as_snp_vector(lead)
    out <- data.table::data.table(
      lead_group = NA_character_,
      snp = snps,
      CHR = NA_character_,
      POS = NA_real_,
      lead_color = NA_character_,
      lead_size = NA_real_,
      lead_stroke = NA_real_,
      flank_color = NA_character_,
      lead_label = NA_character_,
      lead_label_color = NA_character_
    )
  }

  if (!nrow(out)) return(empty)
  out[, snp := as.character(snp)]
  out[, lead_group := as.character(lead_group)]
  out[, CHR := normalize.chrom(CHR)]
  out[, POS := suppressWarnings(as.numeric(POS))]
  out[!nzchar(snp), snp := NA_character_]
  out[!is.finite(POS), POS := NA_real_]
  if ("lead_size" %in% names(out)) out[!is.finite(lead_size), lead_size := NA_real_]
  if ("lead_stroke" %in% names(out)) out[!is.finite(lead_stroke) | lead_stroke < 0, lead_stroke := NA_real_]

  out <- out[(!is.na(snp) & nzchar(snp)) | (!is.na(CHR) & nzchar(CHR) & is.finite(POS))]
  if (!nrow(out)) return(empty)

  out[, lead_id := .I]
  lead_col_in <- lead.color
  lead_col_disabled <- is.null(lead_col_in) ||
    (is.logical(lead_col_in) && length(lead_col_in) == 1L && !isTRUE(lead_col_in)) ||
    (is.character(lead_col_in) && length(lead_col_in) == 1L &&
       tolower(trimws(as.character(lead_col_in)[1])) %in% c("false", "f", "none", "null", "na"))
  if (isTRUE(lead_col_disabled)) {
    out[, lead_color_arg := NA_character_]
  } else {
    out[, lead_color_arg := .manhattan_map_lead_values(snp, lead_group, lead_col_in, default = "#E63946")]
  }
  out[is.na(lead_color) | !nzchar(lead_color), lead_color := lead_color_arg]
  if (!isTRUE(lead_col_disabled)) {
    out[is.na(lead_color) | !nzchar(lead_color), lead_color := "#E63946"]
  }
  out[, lead_color_arg := NULL]

  size_map <- .manhattan_map_lead_sizes(
    snps = out$snp,
    groups = out$lead_group,
    values = lead.size,
    default_main = 1.6,
    default_ceiling = 2.2
  )
  out[, lead_size_main := size_map$lead_size_main]
  out[, lead_size_ceiling := size_map$lead_size_ceiling]
  if ("lead_size" %in% names(out)) {
    out[is.finite(lead_size) & !is.na(lead_size) & lead_size > 0, `:=`(
      lead_size_main = lead_size,
      lead_size_ceiling = lead_size * 1.35
    )]
    out[, lead_size := NULL]
  }

  out[, lead_stroke_arg := .manhattan_map_lead_numeric(snp, lead_group, lead.stroke, default = 1)]
  out[is.na(lead_stroke), lead_stroke := lead_stroke_arg]
  out[!is.finite(lead_stroke) | lead_stroke < 0, lead_stroke := 1]
  out[, lead_stroke_arg := NULL]

  flank_default <- if (is.null(lead.flank.color) || !length(lead.flank.color)) {
    out$lead_color
  } else {
    .manhattan_map_lead_values(out$snp, out$lead_group, lead.flank.color, default = out$lead_color)
  }
  out[is.na(flank_color) | !nzchar(flank_color), flank_color := flank_default]
  out[is.na(flank_color) | !nzchar(flank_color), flank_color := out$lead_color]

  out[, label_color_arg := .manhattan_map_lead_values(snp, lead_group, lead.label.color, default = "grey20")]
  out[is.na(lead_label_color) | !nzchar(lead_label_color), lead_label_color := label_color_arg]
  out[is.na(lead_label_color) | !nzchar(lead_label_color), lead_label_color := "grey20"]
  out[, label_color_arg := NULL]

  if (isTRUE(lead.label)) {
    out[!is.na(snp), lead_label := snp]
  } else if (is.list(lead.label) && !is.data.frame(lead.label) && !data.table::is.data.table(lead.label)) {
    out[, lead_label := .manhattan_map_lead_values(snp, lead_group, lead.label, default = lead_label)]
  }
  out
}

.manhattan_derive_leads <- function(lead_tbl,
                                            dt,
                                            data,
                                            snp_col_use,
                                            lead_label_col_std,
                                            lead.label,
                                            pass_mode,
                                            threshold_zero_mode,
                                            threshold_lines,
                                            lead_color_enabled,
                                            lead_col0) {
  require_pkg("data.table")
  empty_map <- data.table::data.table(
    .row_id = integer(),
    lead_size_main = numeric(),
    lead_size_ceiling = numeric(),
    lead_color = character(),
    lead_stroke = numeric()
  )

  pass_cut <- if (isTRUE(threshold_zero_mode)) {
    -Inf
  } else if (nrow(threshold_lines)) {
    max(threshold_lines$y_raw, na.rm = TRUE)
  } else {
    NA_real_
  }
  if (isTRUE(pass_mode) && (isTRUE(threshold_zero_mode) || is.finite(pass_cut))) {
    pass_idx <- if (isTRUE(threshold_zero_mode)) {
      rep(TRUE, nrow(dt))
    } else {
      dt$yval >= pass_cut
    }
    if (isTRUE(lead_color_enabled) && length(lead_col0)) {
      dt[pass_idx, plot_color := as.character(lead_col0[1])]
    }
    dt[pass_idx, `:=`(plot_alpha = plot_alpha, lead_highlight = TRUE, pass_highlight = TRUE)]
  }

  if (!nrow(lead_tbl)) {
    return(list(lead_tbl = lead_tbl, dt = dt, lead_draw_map = empty_map, pass_cut = pass_cut))
  }

  lead_map <- dt[!is.na(snp) & nzchar(snp), .SD[which.max(yval)], by = snp, .SDcols = c(".row_id", "CHR", "POS", "x", "ydraw", "ceiling")]
  data.table::setnames(
    lead_map,
    c(".row_id", "CHR", "POS", "x", "ydraw", "ceiling"),
    c(".row_id_data", "CHR_data", "POS_data", "x_data", "ydraw_data", "ceiling_data")
  )
  lead_tbl <- merge(lead_tbl, lead_map, by = "snp", all.x = TRUE, sort = FALSE)
  lead_tbl[is.na(CHR) & !is.na(CHR_data), CHR := CHR_data]
  lead_tbl[is.na(POS) & is.finite(POS_data), POS := POS_data]
  if (!(".row_id_data" %in% names(lead_tbl))) lead_tbl[, .row_id_data := NA_integer_]
  if (!("x_data" %in% names(lead_tbl))) lead_tbl[, x_data := NA_real_]
  if (!("ydraw_data" %in% names(lead_tbl))) lead_tbl[, ydraw_data := NA_real_]
  if (!("ceiling_data" %in% names(lead_tbl))) lead_tbl[, ceiling_data := 0L]

  lead_col_map <- dt[!is.na(snp) & nzchar(snp), .SD[which.max(yval)], by = snp, .SDcols = "plot_color"]
  data.table::setnames(lead_col_map, "plot_color", "plot_color_data")
  lead_tbl <- merge(lead_tbl, lead_col_map, by = "snp", all.x = TRUE, sort = FALSE)
  if (!("plot_color_data" %in% names(lead_tbl))) lead_tbl[, plot_color_data := NA_character_]

  if (!is.null(lead_label_col_std) && lead_label_col_std %in% names(dt)) {
    if (identical(lead_label_col_std, "snp")) {
      lead_tbl[(is.na(lead_label) | !nzchar(lead_label)) & !is.na(snp), lead_label := snp]
    } else {
      map_lab <- dt[!is.na(snp) & nzchar(snp), .SD[1], by = snp, .SDcols = lead_label_col_std]
      data.table::setnames(map_lab, lead_label_col_std, "lead_label_from_col")
      lead_tbl <- merge(lead_tbl, map_lab, by = "snp", all.x = TRUE, sort = FALSE)
      lead_tbl[(is.na(lead_label) | !nzchar(lead_label)) & !is.na(lead_label_from_col), lead_label := as.character(lead_label_from_col)]
      lead_tbl[, lead_label_from_col := NULL]
    }
  }

  if (is.character(lead.label) && length(lead.label)) {
    if (length(lead.label) == 1L && !is.na(lead.label[1]) && nzchar(lead.label[1])) {
      lb <- as.character(lead.label[1])
      if (lb %in% names(data)) {
        dt_raw <- if (data.table::is.data.table(data)) data else data.table::as.data.table(data)
        if (all(c(snp_col_use, lb) %in% names(dt_raw))) {
          map_lb <- dt_raw[, .SD[1], by = snp_col_use, .SDcols = lb]
          data.table::setnames(map_lb, c(snp_col_use, lb), c("snp", "lead_label_arg"))
          lead_tbl <- merge(lead_tbl, map_lb, by = "snp", all.x = TRUE, sort = FALSE)
          lead_tbl[!is.na(lead_label_arg), lead_label := as.character(lead_label_arg)]
          lead_tbl[, lead_label_arg := NULL]
        }
      } else {
        vals <- .gcanvas_as_snp_vector(lead.label)
        if (length(vals) && all(vals %in% lead_tbl$snp)) {
          lead_tbl[snp %in% vals, lead_label := snp]
        }
      }
    } else {
      lbv <- as.character(lead.label)
      lbn <- names(lbv)
      if (!is.null(lbn) && any(nzchar(lbn))) {
        lead_tbl[, lead_label := .gcanvas_map_values_by_name(snp, lbv, default = lead_label)]
      } else if (length(lbv) == nrow(lead_tbl)) {
        lead_tbl[, lead_label := lbv]
      } else {
        lead_tbl[, lead_label := rep(lbv, length.out = .N)]
      }
    }
  }

  lead_tbl[, lead_label := as.character(lead_label)]
  lead_tbl[is.na(lead_label) | !nzchar(lead_label), lead_label := NA_character_]

  lead_match_snp <- unique(lead_tbl$snp[!is.na(lead_tbl$snp) & nzchar(lead_tbl$snp)])
  if (length(lead_match_snp)) dt[snp %in% lead_match_snp, lead_highlight := TRUE]

  lead_pos_dt <- unique(lead_tbl[is.finite(POS) & !is.na(CHR) & nzchar(CHR), .(CHR, POS)])
  if (nrow(lead_pos_dt)) {
    idx_lead <- dt[lead_pos_dt, on = .(CHR, POS), which = TRUE]
    idx_lead <- idx_lead[!is.na(idx_lead)]
    if (length(idx_lead)) dt[idx_lead, lead_highlight := TRUE]
  }

  lead_draw_map <- lead_tbl[
    is.finite(.row_id_data) & !is.na(.row_id_data),
    .(
      .row_id = as_int(.row_id_data),
      lead_size_main,
      lead_size_ceiling,
      lead_color,
      lead_stroke
    )
  ]
  data.table::setorderv(lead_draw_map, c(".row_id"), c(1L), na.last = TRUE)
  lead_draw_map <- lead_draw_map[!duplicated(.row_id)]

  list(lead_tbl = lead_tbl, dt = dt, lead_draw_map = lead_draw_map, pass_cut = pass_cut)
}

.manhattan_resolve_lead_positions <- function(lead_tbl, dt) {
  require_pkg("data.table")
  if (!nrow(lead_tbl) || !nrow(dt)) return(lead_tbl)
  lead_map <- dt[!is.na(snp) & nzchar(snp), .SD[which.max(yval)], by = snp, .SDcols = c("CHR", "POS")]
  data.table::setnames(lead_map, c("CHR", "POS"), c("CHR_data", "POS_data"))
  out <- merge(lead_tbl, lead_map, by = "snp", all.x = TRUE, sort = FALSE)
  out[is.na(CHR) & !is.na(CHR_data), CHR := CHR_data]
  out[is.na(POS) & is.finite(POS_data), POS := POS_data]
  out[, c("CHR_data", "POS_data") := NULL]
  out
}

.manhattan_lead_row_idx <- function(dt, lead_tbl) {
  require_pkg("data.table")
  if (!nrow(lead_tbl) || !nrow(dt)) return(integer())
  idx <- integer()
  snp_hit <- unique(as.character(lead_tbl$snp))
  snp_hit <- snp_hit[!is.na(snp_hit) & nzchar(snp_hit)]
  if (length(snp_hit)) idx <- c(idx, dt[snp %in% snp_hit, which = TRUE])
  pos_hit <- unique(lead_tbl[is.finite(POS) & !is.na(CHR) & nzchar(CHR), .(CHR, POS)])
  if (nrow(pos_hit)) idx <- c(idx, dt[pos_hit, on = .(CHR, POS), which = TRUE])
  unique(as.integer(idx[!is.na(idx)]))
}

.manhattan_derive_direction_leads <- function(dt,
                                                      lead_draw_map,
                                                      direction,
                                                      y.col,
                                                      threshold,
                                                      threshold_zero_mode,
                                                      threshold_lines,
                                                      direction_color_inherit,
                                                      dir_col_pos,
                                                      dir_col_neg,
                                                      direction_size0,
                                                      lead_size0,
                                                      lead_color_enabled,
                                                      lead_col0,
                                                      direction_lead_fill_override,
                                                      dir_col_lead_pos,
                                                      dir_col_lead_neg,
                                                      use_ceiling,
                                                      ceiling_nudge_fun) {
  require_pkg("data.table")
  dt[, direction_hit := FALSE]
  dir_dt <- dt[0]
  if (!(isTRUE(direction) && ("direction_sign" %in% names(dt)))) {
    dt[, direction_draw := FALSE]
    return(list(dt = dt, dir_dt = dir_dt, dir_row_idx = integer()))
  }

  direction_p_cut <- NA_real_
  if (is.null(y.col)) {
    th_num <- suppressWarnings(as.numeric(threshold))
    th_num <- th_num[is.finite(th_num) & !is.na(th_num)]
    th_p <- th_num[th_num > 0 & th_num <= 1]
    if (length(th_p)) direction_p_cut <- min(th_p, na.rm = TRUE)
  }
  direction_cut <- if (isTRUE(threshold_zero_mode)) {
    -Inf
  } else if (nrow(threshold_lines)) {
    max(threshold_lines$y_raw, na.rm = TRUE)
  } else {
    -Inf
  }
  if (isTRUE(threshold_zero_mode)) {
    dt[, direction_hit := (direction_sign != 0)]
  } else if (is.null(y.col) && is.finite(direction_p_cut)) {
    dt[, direction_hit := (!is.na(P) & nzchar(P) & .gcanvas_p_to_num(P) <= direction_p_cut & direction_sign != 0)]
  } else if (is.finite(direction_cut)) {
    dt[, direction_hit := (yval >= direction_cut & direction_sign != 0)]
  } else {
    dt[, direction_hit := FALSE]
  }
  dt[, direction_hit := (direction_hit %in% TRUE & lead_highlight %in% TRUE)]
  dir_dt <- dt[direction_hit %in% TRUE, .(.row_id, snp, CHR, POS, x, y_plot, ceiling, direction_sign, lead_highlight, pass_highlight, lead_border_color, plot_color)]
  if (nrow(dir_dt)) {
    if (nrow(lead_draw_map)) {
      dir_dt <- merge(dir_dt, lead_draw_map, by = ".row_id", all.x = TRUE, sort = FALSE)
    } else {
      dir_dt[, `:=`(
        lead_size_main = NA_real_,
        lead_size_ceiling = NA_real_,
        lead_color = NA_character_,
        lead_stroke = NA_real_
      )]
    }
    dir_dt[, `:=`(
      dir_shape = ifelse(direction_sign > 0, 24, 25),
      dir_fill = if (isTRUE(direction_color_inherit)) as.character(plot_color) else ifelse(direction_sign > 0, dir_col_pos, dir_col_neg),
      dir_color = if (isTRUE(direction_color_inherit)) NA_character_ else ifelse(direction_sign > 0, dir_col_pos, dir_col_neg),
      dir_stroke = 0.45,
      dir_size = direction_size0
    )]
    dir_dt[is.na(dir_fill) | !nzchar(dir_fill), dir_fill := as.character(plot_color)]
    dir_dt[is.na(dir_fill) | !nzchar(dir_fill), dir_fill := "grey70"]
    dir_dt[
      lead_highlight %in% TRUE,
      `:=`(
        dir_stroke = data.table::fifelse(is.finite(lead_stroke) & !is.na(lead_stroke) & lead_stroke >= 0, lead_stroke, 1),
        dir_color = "grey20"
      )
    ]
    dir_dt[
      lead_highlight %in% TRUE & pass_highlight %in% FALSE &
        is.finite(lead_size_main) & !is.na(lead_size_main) & lead_size_main > 0,
      dir_size := lead_size_main
    ]
    dir_dt[
      lead_highlight %in% TRUE & pass_highlight %in% FALSE &
        (!is.finite(dir_size) | is.na(dir_size) | dir_size <= 0),
      dir_size := lead_size0[1]
    ]
    if (isTRUE(direction_color_inherit) && isTRUE(lead_color_enabled) && length(lead_col0)) {
      dir_dt[
        lead_highlight %in% TRUE & !is.na(lead_color) & nzchar(lead_color),
        dir_fill := as.character(lead_color)
      ]
      dir_dt[
        lead_highlight %in% TRUE & (is.na(dir_fill) | !nzchar(dir_fill)),
        dir_fill := as.character(lead_col0[1])
      ]
    }
    if (isTRUE(direction_lead_fill_override)) {
      dir_dt[lead_highlight %in% TRUE & direction_sign > 0, dir_fill := dir_col_lead_pos]
      dir_dt[lead_highlight %in% TRUE & direction_sign < 0, dir_fill := dir_col_lead_neg]
    }
    if (!isTRUE(lead_color_enabled)) {
      dir_dt[lead_highlight %in% TRUE, dir_fill := as.character(plot_color)]
    }
    dir_dt <- dir_dt[is.finite(x) & is.finite(y_plot)]
    if (isTRUE(use_ceiling)) {
      dir_dt[, y_dir := {
        n0 <- ceiling_nudge_fun(dir_size)
        n1 <- n0 * 0.28
        y_plot + data.table::fifelse(
          as_int(ceiling) == 1L,
          n0 + data.table::fifelse(direction_sign < 0, n1, 0),
          0
        )
      }]
    } else {
      dir_dt[, y_dir := y_plot]
    }
    data.table::setorderv(dir_dt, c("y_plot", "x"), c(1L, 1L), na.last = TRUE)
  }
  dt[, direction_draw := (direction_hit %in% TRUE)]
  dir_row_idx <- unique(as_int(dir_dt$.row_id[is.finite(dir_dt$.row_id)]))
  list(dt = dt, dir_dt = dir_dt, dir_row_idx = dir_row_idx)
}

.manhattan_draw_leads <- function(p,
                                          lead_tbl,
                                          dir_dt,
                                          dir_row_idx,
                                          suppress_explicit_lead,
                                          lead_color_enabled,
                                          lead_col0,
                                          lead.stroke,
                                          ceiling_nudge_fun,
                                          validate_color_fn,
                                          x_max,
                                          x_pad,
                                          y_top,
                                          y_span,
                                          drag.label,
                                          drag.label.line.color,
                                          drag.label.arrow,
                                          drag.label.linewidth,
                                          lead.label.angle,
                                          lead.label.size,
                                          label_hjust,
                                          label_vjust,
                                          is_angle_90,
                                          is_angle_270,
                                          is_angle_repel_safe,
                                          is_angle_zero,
                                          label_nudge_base,
                                          label_nudge_step,
                                          seed) {
  require_pkg("data.table")
  .manhattan_add_split_lead_points <- function(p, dat, shape_val, size_col) {
    if (!nrow(dat)) return(p)
    stroke_vals <- sort(unique(dat$lead_stroke[is.finite(dat$lead_stroke) & !is.na(dat$lead_stroke) & dat$lead_stroke >= 0]))
    if (!length(stroke_vals)) stroke_vals <- suppressWarnings(as.numeric(lead.stroke))[1]
    if (!length(stroke_vals) || !is.finite(stroke_vals) || is.na(stroke_vals) || stroke_vals < 0) stroke_vals <- 1
    for (st in stroke_vals) {
      dsub <- dat[is.finite(lead_stroke) & !is.na(lead_stroke) & abs(lead_stroke - st) < 1e-8]
      if (!nrow(dsub) && length(stroke_vals) == 1L) {
        dsub <- dat
      }
      if (!nrow(dsub)) next
      p <- p + ggplot2::geom_point(
        data = dsub,
        ggplot2::aes(x = x_data, y = ydraw_plot, fill = I(lead_fill), size = .data[[size_col]]),
        shape = shape_val, color = "grey20", stroke = as.numeric(st), show.legend = FALSE
      )
    }
    p
  }
  lead_pts <- data.table::copy(lead_tbl)
  lead_label_pts <- data.table::copy(lead_tbl)
  if (!("x_data" %in% names(lead_pts))) lead_pts[, x_data := NA_real_]
  if (!("ydraw_data" %in% names(lead_pts))) lead_pts[, ydraw_data := NA_real_]
  if (!("ceiling_data" %in% names(lead_pts))) lead_pts[, ceiling_data := 0L]
  if (!("plot_color_data" %in% names(lead_pts))) lead_pts[, plot_color_data := NA_character_]
  if (!("x_data" %in% names(lead_label_pts))) lead_label_pts[, x_data := NA_real_]
  if (!("ydraw_data" %in% names(lead_label_pts))) lead_label_pts[, ydraw_data := NA_real_]
  if (!("ceiling_data" %in% names(lead_label_pts))) lead_label_pts[, ceiling_data := 0L]

  if (nrow(lead_pts)) {
    dir_lead_snp <- unique(dir_dt[!is.na(snp) & nzchar(snp), snp])
    if (length(dir_lead_snp)) lead_pts <- lead_pts[!(snp %in% dir_lead_snp)]
    dir_lead_pos <- unique(dir_dt[!is.na(CHR) & nzchar(CHR) & is.finite(POS), .(CHR, POS)])
    if (nrow(dir_lead_pos)) {
      lead_pts <- lead_pts[!dir_lead_pos, on = .(CHR_data = CHR, POS_data = POS)]
    }
    if (length(dir_row_idx)) {
      lead_pts <- lead_pts[!(as_int(.row_id_data) %in% dir_row_idx)]
    }
    lead_pts <- lead_pts[is.finite(x_data) & is.finite(ydraw_data)]
    if (nrow(lead_pts) && !isTRUE(suppress_explicit_lead)) {
      lead_pts[, ydraw_plot := ydraw_data + data.table::fifelse(
        as_int(ceiling_data) == 1L,
        ceiling_nudge_fun(lead_size_ceiling),
        0
      )]
      if (isTRUE(lead_color_enabled)) {
        lead_pts[is.na(lead_color) | !nzchar(lead_color), lead_color := if (length(lead_col0)) as.character(lead_col0[1]) else "#E63946"]
        lead_pts[, lead_fill := lead_color]
      } else {
        lead_pts[, lead_fill := as.character(plot_color_data)]
      }
      lead_pts[is.na(lead_fill) | !nzchar(lead_fill), lead_fill := as.character(plot_color_data)]
      lead_pts[is.na(lead_fill) | !nzchar(lead_fill), lead_fill := "grey70"]
      lead_pts <- lead_pts[validate_color_fn(lead_fill)]
      lp0 <- lead_pts[as_int(ceiling_data) == 0L]
      lp1 <- lead_pts[as_int(ceiling_data) == 1L]

      if (nrow(lp0)) p <- .manhattan_add_split_lead_points(p, lp0, 21, "lead_size_main")
      if (nrow(lp1)) p <- .manhattan_add_split_lead_points(p, lp1, 24, "lead_size_ceiling")
    }
  }

  lead_label_pts <- lead_label_pts[is.finite(x_data) & is.finite(ydraw_data)]
  if (nrow(lead_label_pts)) {
    lead_label_pts[, ydraw_plot := ydraw_data + data.table::fifelse(
      as_int(ceiling_data) == 1L,
      ceiling_nudge_fun(lead_size_ceiling),
      0
    )]
  }
  label_dt <- lead_label_pts[!is.na(lead_label) & nzchar(lead_label)]
  if (!nrow(label_dt)) return(p)

  label_ang0 <- suppressWarnings(as.numeric(lead.label.angle))[1]
  if (!is.finite(label_ang0) || is.na(label_ang0)) label_ang0 <- 0
  if (isTRUE(drag.label)) {
    label_drag <- data.table::copy(label_dt)
    data.table::setorder(label_drag, x_data, -ydraw_plot)
    min_dx <- x_max * 0.03
    if (!is.finite(min_dx) || min_dx <= 0) min_dx <- 1
    label_drag[, x_label := as.numeric(x_data)]
    if (nrow(label_drag) >= 2L) {
      for (i in 2:nrow(label_drag)) {
        label_drag$x_label[i] <- max(label_drag$x_label[i], label_drag$x_label[i - 1L] + min_dx)
      }
    }
    x_hi <- x_max + x_pad * 0.9
    if (max(label_drag$x_label, na.rm = TRUE) > x_hi) {
      shift_all <- max(label_drag$x_label, na.rm = TRUE) - x_hi
      label_drag[, x_label := x_label - shift_all]
    }
    label_drag[, x_label := pmax(x_pad * 0.5, pmin(x_label, x_hi))]
    bend_y <- y_top + y_span * 0.02
    label_y <- y_top + y_span * 0.06
    text_gap <- max(0.02, y_span * 0.004)
    tip_gap <- max(0.06, y_span * 0.02)
    tri_extra <- max(0.03, y_span * 0.01)
    label_drag[, y_tip := ydraw_plot + tip_gap + ifelse(as_int(ceiling_data) == 1L, tri_extra, 0)]
    label_drag[, `:=`(y_bend = bend_y, y_label = label_y, y_text = label_y + text_gap)]
    label_drag[, hjust_drag := ifelse(x_label >= x_data, 0, 1)]
    label_drag[, drag_seg_col := .gcanvas_map_values_by_name(snp, drag.label.line.color, default = "grey80")]
    label_drag[is.na(drag_seg_col) | !nzchar(drag_seg_col), drag_seg_col := "grey80"]
    label_drag[, drag_arrow_col := drag_seg_col]
    arrow_spec <- if (isTRUE(drag.label.arrow)) grid::arrow(length = grid::unit(0.12, "cm"), type = "closed") else NULL

    p <- p +
      ggplot2::geom_segment(
        data = label_drag,
        ggplot2::aes(x = x_data, y = y_bend, xend = x_data, yend = y_tip, color = I(drag_arrow_col)),
        linewidth = drag.label.linewidth,
        alpha = 0.8,
        arrow = arrow_spec,
        show.legend = FALSE
      ) +
      ggplot2::geom_segment(
        data = label_drag,
        ggplot2::aes(x = x_data, y = y_bend, xend = x_label, yend = y_label, color = I(drag_seg_col)),
        linewidth = 0.35,
        alpha = 0.8,
        show.legend = FALSE
      ) +
      ggplot2::geom_text(
        data = label_drag,
        ggplot2::aes(x = x_label, y = y_text, label = lead_label, color = I(lead_label_color), hjust = hjust_drag),
        angle = lead.label.angle,
        vjust = 0.5,
        size = lead.label.size,
        show.legend = FALSE
      )
    return(p)
  }

  label_dt <- data.table::copy(label_dt)
  label_dt_ceiling <- label_dt[as_int(ceiling_data) == 1L]
  label_dt_main <- label_dt[as_int(ceiling_data) != 1L | is.na(ceiling_data)]

  if (nrow(label_dt_main)) {
    seg_col <- NA_character_
    seg_min_len <- Inf
    if ((isTRUE(is_angle_90) || isTRUE(is_angle_270))) {
      data.table::setorder(label_dt_main, -ydraw_plot, x_data)
      base_y_shift <- label_nudge_base * 1.2
      step_y_shift <- label_nudge_step
      label_dt_main[, y_label := ydraw_plot + base_y_shift + (.N - seq_len(.N)) * step_y_shift]
      x_shift <- x_max * 0.0015
      label_dt_main[, x_label := x_data + if (isTRUE(is_angle_90)) x_shift else -x_shift]
      p <- p + ggplot2::geom_text(
        data = label_dt_main,
        ggplot2::aes(x = x_label, y = y_label, label = lead_label, color = I(lead_label_color)),
        size = lead.label.size, angle = lead.label.angle, hjust = label_hjust, vjust = label_vjust, show.legend = FALSE
      )
    } else if (requireNamespace("ggrepel", quietly = TRUE) && isTRUE(is_angle_repel_safe)) {
      nudge_y0 <- if (isTRUE(is_angle_zero)) label_nudge_base else label_nudge_base * 0.8
      p <- p + ggrepel::geom_text_repel(
        data = label_dt_main,
        ggplot2::aes(x = x_data, y = ydraw_plot, label = lead_label, color = I(lead_label_color)),
        seed = seed,
        size = lead.label.size,
        angle = lead.label.angle,
        hjust = label_hjust,
        vjust = label_vjust,
        box.padding = 0.35,
        point.padding = 0.2,
        nudge_y = nudge_y0,
        direction = "y",
        min.segment.length = seg_min_len,
        segment.color = seg_col,
        show.legend = FALSE
      )
    } else {
      if (!requireNamespace("ggrepel", quietly = TRUE)) {
        .gcanvas_warn_msg("Package 'ggrepel' is not installed. Falling back to geom_text for lead labels.")
      }
      data.table::setorder(label_dt_main, -ydraw_plot, x_data)
      base_y_shift <- if (isTRUE(is_angle_zero)) label_nudge_base else label_nudge_base * 1.1
      step_y_shift <- label_nudge_step
      label_dt_main[, y_label := ydraw_plot + base_y_shift + (.N - seq_len(.N)) * step_y_shift]
      x_shift <- x_max * 0.002
      label_dt_main[, x_label := x_data + if (label_hjust <= 0.05) x_shift else if (label_hjust >= 0.95) -x_shift else 0]
      p <- p + ggplot2::geom_text(
        data = label_dt_main,
        ggplot2::aes(x = x_label, y = y_label, label = lead_label, color = I(lead_label_color)),
        size = lead.label.size, angle = lead.label.angle, hjust = label_hjust, vjust = label_vjust, show.legend = FALSE
      )
    }
  }

  if (nrow(label_dt_ceiling)) {
    data.table::setorder(label_dt_ceiling, x_data, -ydraw_plot)
    ang_fac <- abs(sin(label_ang0 * pi / 180))
    base_out <- max(0.05, y_span * (0.03 + 0.03 * ang_fac))
    step_out <- max(0.010, y_span * (0.007 + 0.010 * ang_fac))
    out_start <- y_top + max(0.03, y_span * (0.015 + 0.02 * ang_fac))
    label_dt_ceiling[, y_anchor := ydraw_plot + base_out]
    label_dt_ceiling[, y_outside := out_start + ((seq_len(.N) - 1L) %% 3L) * step_out]
    label_dt_ceiling[, y_label := pmax(y_anchor, y_outside)]
    x_shift_ce <- x_max * 0.002
    label_dt_ceiling[, x_label := x_data + if (isTRUE(is_angle_90)) x_shift_ce * 2 else if (isTRUE(is_angle_270)) -x_shift_ce * 2 else 0]
    p <- p + ggplot2::geom_text(
      data = label_dt_ceiling,
      ggplot2::aes(x = x_label, y = y_label, label = lead_label, color = I(lead_label_color)),
      size = lead.label.size, angle = lead.label.angle, hjust = label_hjust, vjust = label_vjust, show.legend = FALSE
    )
  }

  p
}

.manhattan_draw_chr_ticks <- function(p,
                                              chr_map,
                                              special_chr_ticks,
                                              drag_chrom_tick,
                                              tick_len_auto,
                                              linewidth = 0.5,
                                              color = "grey20") {
  if (!length(special_chr_ticks)) return(p)
  tick_dt <- chr_map[CHR %in% special_chr_ticks, .(CHR, x_tick = tick)]
  tick_dt[, chr_num := as_int(CHR)]
  data.table::setorder(tick_dt, chr_num)
  tick_dt[, chr_num := NULL]
  tick_len <- if (is.finite(drag_chrom_tick) && !is.na(drag_chrom_tick) && drag_chrom_tick > 0) {
    as.numeric(drag_chrom_tick)
  } else {
    tick_len_auto
  }
  p + ggplot2::geom_segment(
    data = tick_dt,
    ggplot2::aes(x = x_tick, xend = x_tick, y = 0, yend = -tick_len),
    inherit.aes = FALSE,
    linewidth = linewidth,
    color = color
  )
}

.manhattan_partition_draw_data <- function(dt,
                                                   use_ceiling,
                                                   direction,
                                                   dir_row_idx,
                                                   base_shape,
                                                   ceil_shape,
                                                   base_size,
                                                   ceil_size,
                                                   ceiling_nudge_fun) {
  dt_main <- dt[(ceiling != 1L) & !(lead_draw %in% TRUE)]
  dt_ceiling <- dt[(ceiling == 1L) & !(lead_draw %in% TRUE)]
  dt_highlight <- dt_main[(flank_highlight | lead_highlight) & !(direction_draw %in% TRUE)]
  dt_base <- dt_main[!(flank_highlight | lead_highlight) & !(direction_draw %in% TRUE)]
  dt_ceiling_plain <- dt_ceiling[(ceiling_marker %in% TRUE) & !(direction_draw %in% TRUE)]
  if (nrow(dt_ceiling_plain)) {
    n0 <- ceiling_nudge_fun(ceil_size)
    n1 <- n0 * 0.28
    if (isTRUE(direction) && ("direction_sign" %in% names(dt_ceiling_plain))) {
      dt_ceiling_plain[, y_ceiling_plot := y_plot + n0 + data.table::fifelse(direction_sign < 0, n1, 0)]
    } else {
      dt_ceiling_plain[, y_ceiling_plot := y_plot + n0]
    }
  }
  if (length(dir_row_idx)) {
    dt_ceiling_plain <- dt_ceiling_plain[!(.row_id %in% dir_row_idx)]
  }
  dt_base[, shape_plot := ifelse(ceiling == 1, ceil_shape, base_shape)]
  dt_base[, size_plot := ifelse(ceiling == 1, ceil_size, base_size)]
  dt_highlight[, shape_plot := ifelse(ceiling == 1, ceil_shape, base_shape)]
  dt_highlight[, size_plot := ifelse(ceiling == 1, ceil_size, base_size)]
  list(
    dt_base = dt_base,
    dt_highlight = dt_highlight,
    dt_ceiling_plain = dt_ceiling_plain
  )
}

.manhattan_init_plot <- function(dt_base,
                                         x_axis_title,
                                         x_breaks,
                                         x_labels,
                                         y_axis_title,
                                         y_breaks_draw,
                                         y_breaks_labels,
                                         x_pad,
                                         x_max,
                                         y_top,
                                         title,
                                         axis_line_x_elem,
                                         axis_line_y_elem,
                                         tick.size,
                                         panel_border_elem,
                                         grid_major_x_elem,
                                         grid_major_y_elem,
                                         grid_minor_x_elem,
                                         grid_minor_y_elem,
                                         x.text.size,
                                         y.text.size,
                                         x.title.nudge,
                                         top_margin_pt,
                                         bottom_margin_pt,
                                         left_margin_pt) {
  ggplot2::ggplot() +
    ggplot2::geom_point(
      data = dt_base,
      ggplot2::aes(
        x = x, y = y_plot,
        shape = shape_plot,
        size = size_plot,
        color = I(plot_color),
        alpha = I(plot_alpha)
      ),
      show.legend = FALSE
    ) +
    ggplot2::scale_shape_identity() +
    ggplot2::scale_size_identity() +
    ggplot2::scale_x_continuous(name = x_axis_title, breaks = x_breaks, labels = x_labels, expand = c(0, 0)) +
    ggplot2::scale_y_continuous(
      name = y_axis_title,
      breaks = y_breaks_draw,
      labels = y_breaks_labels,
      expand = ggplot2::expansion(mult = c(0, 0), add = c(0, 0))
    ) +
    ggplot2::coord_cartesian(xlim = c(0 - x_pad, x_max + x_pad), ylim = c(0, y_top), clip = "off") +
    ggplot2::ggtitle(title) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.line.x = axis_line_x_elem,
      axis.line.y = axis_line_y_elem,
      axis.ticks = if (is.finite(tick.size) && !is.na(tick.size) && tick.size <= 0) ggplot2::element_blank() else ggplot2::element_line(),
      axis.ticks.length = if (is.finite(tick.size) && !is.na(tick.size) && tick.size > 0) grid::unit(tick.size, "pt") else grid::unit(2.75, "pt"),
      panel.border = panel_border_elem,
      panel.grid.major.x = grid_major_x_elem,
      panel.grid.major.y = grid_major_y_elem,
      panel.grid.minor.x = grid_minor_x_elem,
      panel.grid.minor.y = grid_minor_y_elem,
      axis.title = ggplot2::element_text(size = 16),
      axis.title.x = ggplot2::element_text(size = 16, margin = ggplot2::margin(t = 6 - x.title.nudge)),
      axis.text.x = ggplot2::element_text(size = x.text.size),
      axis.text.y = ggplot2::element_text(size = y.text.size),
      plot.margin = ggplot2::margin(t = top_margin_pt, r = 8, b = bottom_margin_pt, l = left_margin_pt)
    )
}

.manhattan_draw_highlight_points <- function(p, dt_highlight) {
  if (!nrow(dt_highlight)) return(p)
  p + ggplot2::geom_point(
    data = dt_highlight,
    ggplot2::aes(
      x = x, y = y_plot,
      shape = shape_plot,
      size = size_plot,
      color = I(plot_color),
      alpha = I(plot_alpha)
    ),
    show.legend = FALSE
  )
}

.manhattan_draw_ceiling_points <- function(p,
                                                   dt_ceiling_plain,
                                                   direction,
                                                   direction_color_inherit,
                                                   dir_col_pos,
                                                   dir_col_neg,
                                                   ceil_shape,
                                                   ceil_size) {
  if (!nrow(dt_ceiling_plain)) return(p)
  if (isTRUE(direction) && ("direction_sign" %in% names(dt_ceiling_plain))) {
    ce_dir <- dt_ceiling_plain[is.finite(direction_sign) & direction_sign != 0]
    ce_plain <- dt_ceiling_plain[!(is.finite(direction_sign) & direction_sign != 0)]
    if (nrow(ce_plain)) {
      p <- p + ggplot2::geom_point(
        data = ce_plain,
        ggplot2::aes(
          x = x, y = y_ceiling_plot,
          color = I(plot_color),
          alpha = I(plot_alpha)
        ),
        shape = ceil_shape,
        size = ceil_size,
        show.legend = FALSE
      )
    }
    if (nrow(ce_dir)) {
      ce_dir[, dir_fill_ce := if (isTRUE(direction_color_inherit)) {
        as.character(plot_color)
      } else ifelse(direction_sign > 0, dir_col_pos, dir_col_neg)]
      ce_dir[is.na(dir_fill_ce) | !nzchar(dir_fill_ce), dir_fill_ce := as.character(plot_color)]
      ce_dir[is.na(dir_fill_ce) | !nzchar(dir_fill_ce), dir_fill_ce := "grey70"]
      ce_pos <- ce_dir[direction_sign > 0]
      ce_neg <- ce_dir[direction_sign < 0]
      if (nrow(ce_pos)) {
        p <- p + ggplot2::geom_point(
          data = ce_pos,
          ggplot2::aes(
            x = x, y = y_ceiling_plot,
            fill = I(dir_fill_ce),
            alpha = I(plot_alpha)
          ),
          shape = 24,
          size = ceil_size,
          color = NA,
          stroke = 0,
          show.legend = FALSE
        )
      }
      if (nrow(ce_neg)) {
        p <- p + ggplot2::geom_point(
          data = ce_neg,
          ggplot2::aes(
            x = x, y = y_ceiling_plot,
            fill = I(dir_fill_ce),
            alpha = I(plot_alpha)
          ),
          shape = 25,
          size = ceil_size,
          color = NA,
          stroke = 0,
          show.legend = FALSE
        )
      }
    }
    return(p)
  }
  p + ggplot2::geom_point(
    data = dt_ceiling_plain,
    ggplot2::aes(
      x = x, y = y_ceiling_plot,
      color = I(plot_color),
      alpha = I(plot_alpha)
    ),
    shape = ceil_shape,
    size = ceil_size,
    show.legend = FALSE
  )
}

.manhattan_draw_reference_lines <- function(p,
                                                    line_dt,
                                                    line.alpha,
                                                    y.rescale.line,
                                                    use_y_rescale,
                                                    y_rescale_at,
                                                    y_map_fun,
                                                    y.rescale.line.linewidth,
                                                    y.rescale.line.color,
                                                    y.rescale.line.type,
                                                    y.ceiling.line,
                                                    use_ceiling,
                                                    y_ceiling0,
                                                    panel.box,
                                                    y.ceiling.line.linewidth,
                                                    y.ceiling.line.color,
                                                    y.ceiling.line.type) {
  if (nrow(line_dt)) {
    for (i in seq_len(nrow(line_dt))) {
      p <- p + ggplot2::geom_hline(
        yintercept = line_dt$y[i],
        linewidth = line_dt$linewidth[i],
        colour = line_dt$color[i],
        linetype = line_dt$type[i],
        alpha = line.alpha
      )
    }
  }
  y_ref_line_alpha <- 0.9
  if (isTRUE(y.rescale.line) && isTRUE(use_y_rescale) && is.finite(y_rescale_at) && !is.na(y_rescale_at)) {
    y_rescale_plot <- y_map_fun(y_rescale_at)
    p <- p + ggplot2::geom_hline(
      yintercept = y_rescale_plot,
      linewidth = y.rescale.line.linewidth,
      colour = y.rescale.line.color,
      linetype = y.rescale.line.type,
      alpha = y_ref_line_alpha
    )
  }
  if (isTRUE(y.ceiling.line) && isTRUE(use_ceiling) && is.finite(y_ceiling0) && !isTRUE(panel.box)) {
    y_ceiling_plot <- y_map_fun(y_ceiling0)
    p <- p + ggplot2::geom_hline(
      yintercept = y_ceiling_plot,
      linewidth = y.ceiling.line.linewidth,
      colour = y.ceiling.line.color,
      linetype = y.ceiling.line.type,
      alpha = y_ref_line_alpha
    )
  }
  p
}

.manhattan_draw_direction_points <- function(p,
                                                     dir_dt,
                                                     validate_color_fn,
                                                     direction_size0,
                                                     lead.stroke) {
  if (!nrow(dir_dt)) return(p)
  dir_dt <- dir_dt[validate_color_fn(dir_fill)]
  dir_dt[is.na(dir_color) | !nzchar(dir_color), dir_color := "grey20"]
  dir_dt <- dir_dt[validate_color_fn(dir_color)]
  dir_nonlead <- dir_dt[lead_highlight %in% FALSE]
  dir_lead <- dir_dt[lead_highlight %in% TRUE]
  if (nrow(dir_nonlead)) {
    p <- p + ggplot2::geom_point(
      data = dir_nonlead[direction_sign > 0],
      ggplot2::aes(x = x, y = y_dir, fill = I(dir_fill)),
      shape = 24, size = direction_size0,
      color = NA, stroke = 0, show.legend = FALSE
    ) +
      ggplot2::geom_point(
        data = dir_nonlead[direction_sign < 0],
        ggplot2::aes(x = x, y = y_dir, fill = I(dir_fill)),
        shape = 25, size = direction_size0,
        color = NA, stroke = 0, show.legend = FALSE
      )
  }
  if (nrow(dir_lead)) {
    for (sgn in c(1, -1)) {
      dsub <- dir_lead[direction_sign == sgn]
      if (!nrow(dsub)) next
      shp <- if (sgn > 0) 24 else 25
      for (sz in sort(unique(dsub$dir_size))) {
        if (!is.finite(sz) || is.na(sz) || sz <= 0) next
        dsz <- dsub[is.finite(dir_size) & abs(dir_size - sz) < 1e-8]
        if (!nrow(dsz)) next
        stroke_vals <- sort(unique(dsz$dir_stroke[is.finite(dsz$dir_stroke) & !is.na(dsz$dir_stroke) & dsz$dir_stroke >= 0]))
        if (!length(stroke_vals)) stroke_vals <- suppressWarnings(as.numeric(lead.stroke))[1]
        if (!length(stroke_vals) || !is.finite(stroke_vals) || is.na(stroke_vals) || stroke_vals < 0) stroke_vals <- 1
        for (st in stroke_vals) {
          dst <- dsz[is.finite(dir_stroke) & !is.na(dir_stroke) & abs(dir_stroke - st) < 1e-8]
          if (!nrow(dst) && length(stroke_vals) == 1L) dst <- dsz
          if (!nrow(dst)) next
          p <- p + ggplot2::geom_point(
            data = dst,
            ggplot2::aes(x = x, y = y_dir, fill = I(dir_fill), color = I(dir_color)),
            shape = shp, size = as.numeric(sz),
            stroke = as.numeric(st), show.legend = FALSE
          )
        }
      }
    }
  }
  p
}

.manhattan_draw_overlays <- function(p,
                                             y_span,
                                             point.size,
                                             chr_map,
                                             special_chr_ticks,
                                             drag_chrom_tick,
                                             tick_len_auto) {
  mask_depth <- max(0.05, y_span * 0.012 * (point.size / 1.4))
  mask_top <- -max(1e-4, y_span * 2e-4)
  p <- p + ggplot2::annotate(
    "rect",
    xmin = -Inf, xmax = Inf,
    ymin = -mask_depth, ymax = mask_top,
    fill = "white", color = NA
  )
  .manhattan_draw_chr_ticks(
    p = p,
    chr_map = chr_map,
    special_chr_ticks = special_chr_ticks,
    drag_chrom_tick = drag_chrom_tick,
    tick_len_auto = tick_len_auto,
    linewidth = 0.5,
    color = "grey20"
  )
}

.manhattan_resolve_options <- function(silent,
                                               panel.box,
                                               x.text.size,
                                               y.text.size,
                                               x.title.nudge,
                                               tick.size,
                                               panel.space,
                                               panel.space.top,
                                               panel.space.bottom,
                                               panel.space.left,
                                               grid,
                                               grid.major,
                                               grid.minor,
                                               grid.major.x,
                                               grid.major.y,
                                               grid.minor.x,
                                               grid.minor.y,
                                               direction,
                                               direction.color,
                                               direction.color.lead,
                                               line.alpha,
                                               alpha,
                                               lead.flank.alpha,
                                               point.size,
                                               direction.size,
                                               drag.chrom,
                                               drag.chrom.tick,
                                               drag.label,
                                               drag.label.arrow,
                                               chroms,
                                               chroms.drop,
                                               drag.label.linewidth,
                                               lead.stroke,
                                               lead.label.size,
                                               lead.label.angle,
                                               seed,
                                               line.x,
                                               y.ceiling.line,
                                               y.rescale.line,
                                               y.ceiling.line.color,
                                               y.ceiling.line.type,
                                               y.ceiling.line.linewidth,
                                               y.rescale.line.color,
                                               y.rescale.line.type,
                                               y.rescale.line.linewidth,
                                               lead.color) {
  silent <- isTRUE(silent)
  panel.box <- isTRUE(panel.box)
  x.text.size <- suppressWarnings(as.numeric(x.text.size))[1]
  if (!is.finite(x.text.size) || is.na(x.text.size) || x.text.size <= 0) x.text.size <- 14
  y.text.size <- suppressWarnings(as.numeric(y.text.size))[1]
  if (!is.finite(y.text.size) || is.na(y.text.size) || y.text.size <= 0) y.text.size <- 14
  x.title.nudge <- suppressWarnings(as.numeric(x.title.nudge))[1]
  if (!is.finite(x.title.nudge) || is.na(x.title.nudge)) x.title.nudge <- 0
  tick.size <- suppressWarnings(as.numeric(tick.size))[1]
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
  grid <- isTRUE(grid)
  grid.major <- isTRUE(grid.major) || grid
  grid.minor <- isTRUE(grid.minor) || grid
  grid.major.x <- isTRUE(grid.major.x) || grid.major
  grid.major.y <- isTRUE(grid.major.y) || grid.major
  grid.minor.x <- isTRUE(grid.minor.x) || grid.minor
  grid.minor.y <- isTRUE(grid.minor.y) || grid.minor
  direction <- isTRUE(direction)
  line.alpha <- suppressWarnings(as.numeric(line.alpha))[1]
  if (!is.finite(line.alpha) || is.na(line.alpha)) line.alpha <- 0.8
  line.alpha <- max(0, min(1, line.alpha))
  direction_color_inherit <- is.null(direction.color) ||
    (is.logical(direction.color) && length(direction.color) == 1L && !isTRUE(direction.color)) ||
    (is.character(direction.color) && length(direction.color) == 1L &&
       tolower(trimws(as.character(direction.color)[1])) %in% c("false", "f", "none", "null", "na"))
  dc <- direction.color
  if (is.list(dc) && !is.data.frame(dc) && !data.table::is.data.table(dc)) dc <- unlist(dc, use.names = TRUE)
  dc <- as.character(dc)
  dc <- dc[!is.na(dc) & nzchar(dc)]
  if (!length(dc)) {
    direction_color_inherit <- TRUE
    dir_col_pos <- NA_character_
    dir_col_neg <- NA_character_
  } else {
    dcn <- names(dc)
    if (!is.null(dcn) && any(nzchar(dcn))) {
      pos_hit <- match("+", dcn)
      neg_hit <- match("-", dcn)
      if (is.na(pos_hit)) pos_hit <- match("pos", tolower(dcn))
      if (is.na(neg_hit)) neg_hit <- match("neg", tolower(dcn))
      dir_col_pos <- if (!is.na(pos_hit)) dc[pos_hit] else dc[1]
      dir_col_neg <- if (!is.na(neg_hit)) dc[neg_hit] else if (length(dc) >= 2L) dc[2] else dc[1]
    } else {
      dir_col_pos <- dc[1]
      dir_col_neg <- if (length(dc) >= 2L) dc[2] else dc[1]
    }
  }
  dcl <- direction.color.lead
  if (is.list(dcl) && !is.data.frame(dcl) && !data.table::is.data.table(dcl)) dcl <- unlist(dcl, use.names = TRUE)
  dcl <- as.character(dcl)
  dcl <- dcl[!is.na(dcl) & nzchar(dcl)]
  direction_lead_fill_override <- (!isTRUE(direction_color_inherit) || length(dcl) > 0L)
  if (length(dcl)) {
    dcln <- names(dcl)
    if (!is.null(dcln) && any(nzchar(dcln))) {
      pos_hit <- match("+", dcln)
      neg_hit <- match("-", dcln)
      if (is.na(pos_hit)) pos_hit <- match("pos", tolower(dcln))
      if (is.na(neg_hit)) neg_hit <- match("neg", tolower(dcln))
      dir_col_lead_pos <- if (!is.na(pos_hit)) dcl[pos_hit] else dcl[1]
      dir_col_lead_neg <- if (!is.na(neg_hit)) dcl[neg_hit] else if (length(dcl) >= 2L) dcl[2] else dcl[1]
    } else {
      dir_col_lead_pos <- dcl[1]
      dir_col_lead_neg <- if (length(dcl) >= 2L) dcl[2] else dcl[1]
    }
  } else {
    dir_col_lead_pos <- dir_col_pos
    dir_col_lead_neg <- dir_col_neg
  }
  alpha <- suppressWarnings(as.numeric(alpha))[1]
  if (!is.finite(alpha)) alpha <- 1
  alpha <- max(0, min(1, alpha))
  lead.flank.alpha <- suppressWarnings(as.numeric(lead.flank.alpha))[1]
  if (!is.finite(lead.flank.alpha)) lead.flank.alpha <- alpha
  lead.flank.alpha <- max(0, min(1, lead.flank.alpha))
  point.size <- suppressWarnings(as.numeric(point.size))[1]
  if (!is.finite(point.size) || point.size <= 0) point.size <- 1.4
  direction_size0 <- suppressWarnings(as.numeric(direction.size))[1]
  if (!is.finite(direction_size0) || is.na(direction_size0) || direction_size0 <= 0) direction_size0 <- point.size * 1.45
  drag_chrom_vec <- character()
  if (is.null(drag.chrom) || length(drag.chrom) == 0L) {
    drag_chrom_vec <- character()
  } else if (is.logical(drag.chrom) && length(drag.chrom) == 1L) {
    if (isTRUE(drag.chrom)) {
      drag_chrom_vec <- c("19", "21")
    } else {
      drag_chrom_vec <- character()
    }
  } else {
    drag_chrom_vec <- normalize.chrom(as.character(drag.chrom))
    drag_chrom_vec <- unique(drag_chrom_vec[!is.na(drag_chrom_vec) & nzchar(drag_chrom_vec)])
  }
  drag_chrom_tick <- {
    dct0 <- drag.chrom.tick
    if (is.null(dct0) || length(dct0) == 0L) {
      NA_real_
    } else if (is.character(dct0) && tolower(trimws(as.character(dct0)[1])) %in% c("auto", "default")) {
      NA_real_
    } else if (is.logical(dct0)) {
      if (isTRUE(dct0[1])) NA_real_ else NA_real_
    } else {
      dct <- suppressWarnings(as.numeric(dct0))[1]
      if (is.finite(dct) && !is.na(dct) && dct > 0) dct else NA_real_
    }
  }
  drag.label <- isTRUE(drag.label)
  drag.label.arrow <- isTRUE(drag.label.arrow)
  chroms.drop <- isTRUE(chroms.drop)
  chroms_auto <- is.null(chroms) || length(chroms) == 0L ||
    (length(chroms) == 1L && is.character(chroms) && tolower(trimws(as.character(chroms)[1])) == "auto")
  chroms_req <- character()
  if (!isTRUE(chroms_auto)) {
    chroms_req <- normalize.chrom(chroms)
    chroms_req <- chroms_req[!is.na(chroms_req) & nzchar(chroms_req)]
    if (length(chroms_req)) {
      chroms_req <- unique(chroms_req)
      chroms_req <- chroms_req[order(rank.chrom(chroms_req), chroms_req)]
    }
  }
  drag.label.linewidth <- suppressWarnings(as.numeric(drag.label.linewidth))[1]
  if (!is.finite(drag.label.linewidth) || is.na(drag.label.linewidth) || drag.label.linewidth <= 0) {
    drag.label.linewidth <- 0.35
  }
  lead.stroke <- suppressWarnings(as.numeric(lead.stroke))[1]
  if (!is.finite(lead.stroke) || lead.stroke < 0) lead.stroke <- 1
  lead.label.size <- suppressWarnings(as.numeric(lead.label.size))[1]
  if (!is.finite(lead.label.size) || lead.label.size <= 0) lead.label.size <- 3.5
  lead.label.angle <- suppressWarnings(as.numeric(lead.label.angle))[1]
  if (!is.finite(lead.label.angle) || is.na(lead.label.angle)) lead.label.angle <- 0
  label_angle_norm <- lead.label.angle %% 360
  if (!is.finite(label_angle_norm) || is.na(label_angle_norm)) label_angle_norm <- 0
  if (label_angle_norm < 0) label_angle_norm <- label_angle_norm + 360
  is_angle_90 <- abs(label_angle_norm - 90) < 1e-8
  is_angle_270 <- abs(label_angle_norm - 270) < 1e-8
  is_angle_zero <- abs(label_angle_norm) < 1e-8 || abs(label_angle_norm - 360) < 1e-8
  is_angle_repel_safe <- abs((label_angle_norm %% 90)) < 1e-8 || abs((label_angle_norm %% 90) - 90) < 1e-8
  label_hjust <- if (is_angle_zero) {
    0.5
  } else if (is_angle_90 || is_angle_270) {
    if (is_angle_90) 0 else 1
  } else if (label_angle_norm <= 90 || label_angle_norm >= 270) {
    0
  } else {
    1
  }
  label_vjust <- if (is_angle_zero) {
    -0.6
  } else if (is_angle_90 || is_angle_270) {
    0.5
  } else {
    0
  }
  seed <- as_int(seed)
  if (is.na(seed)) seed <- 23L
  line.x <- isTRUE(line.x)
  y.ceiling.line <- isTRUE(y.ceiling.line)
  y.rescale.line <- isTRUE(y.rescale.line)
  y.ceiling.line.color <- as.character(y.ceiling.line.color)[1]
  if (is.na(y.ceiling.line.color) || !nzchar(y.ceiling.line.color)) y.ceiling.line.color <- "grey80"
  y.ceiling.line.type <- tolower(trimws(as.character(y.ceiling.line.type)[1]))
  if (is.na(y.ceiling.line.type) || !nzchar(y.ceiling.line.type)) y.ceiling.line.type <- "solid"
  y.ceiling.line.type <- ifelse(y.ceiling.line.type %in% c("dottd", "dot", "dott"), "dotted", y.ceiling.line.type)
  y.ceiling.line.linewidth <- suppressWarnings(as.numeric(y.ceiling.line.linewidth))[1]
  if (!is.finite(y.ceiling.line.linewidth) || is.na(y.ceiling.line.linewidth) || y.ceiling.line.linewidth <= 0) y.ceiling.line.linewidth <- 0.5
  y.rescale.line.color <- as.character(y.rescale.line.color)[1]
  if (is.na(y.rescale.line.color) || !nzchar(y.rescale.line.color)) y.rescale.line.color <- "grey80"
  y.rescale.line.type <- tolower(trimws(as.character(y.rescale.line.type)[1]))
  if (is.na(y.rescale.line.type) || !nzchar(y.rescale.line.type)) y.rescale.line.type <- "solid"
  y.rescale.line.type <- ifelse(y.rescale.line.type %in% c("dottd", "dot", "dott"), "dotted", y.rescale.line.type)
  y.rescale.line.linewidth <- suppressWarnings(as.numeric(y.rescale.line.linewidth))[1]
  if (!is.finite(y.rescale.line.linewidth) || is.na(y.rescale.line.linewidth) || y.rescale.line.linewidth <= 0) y.rescale.line.linewidth <- 0.5
  lead_color_enabled <- !(is.null(lead.color) ||
    (is.logical(lead.color) && length(lead.color) == 1L && !isTRUE(lead.color)) ||
    (is.character(lead.color) && length(lead.color) == 1L &&
       tolower(trimws(as.character(lead.color)[1])) %in% c("false", "f", "none", "null", "na")))
  lead_col0 <- if (isTRUE(lead_color_enabled)) {
    if (is.list(lead.color) && !is.data.frame(lead.color) && !data.table::is.data.table(lead.color)) {
      unlist(lead.color, use.names = FALSE)
    } else {
      lead.color
    }
  } else character(0)
  lead_col0 <- as.character(lead_col0)
  lead_col0 <- lead_col0[!is.na(lead_col0) & nzchar(lead_col0)]
  if (isTRUE(lead_color_enabled) && !length(lead_col0)) lead_col0 <- "#E63946"
  list(
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
    line.alpha = line.alpha,
    direction_color_inherit = direction_color_inherit,
    dir_col_pos = dir_col_pos,
    dir_col_neg = dir_col_neg,
    direction_lead_fill_override = direction_lead_fill_override,
    dir_col_lead_pos = dir_col_lead_pos,
    dir_col_lead_neg = dir_col_lead_neg,
    alpha = alpha,
    lead.flank.alpha = lead.flank.alpha,
    point.size = point.size,
    direction_size0 = direction_size0,
    drag_chrom_vec = drag_chrom_vec,
    drag_chrom_tick = drag_chrom_tick,
    drag.label = drag.label,
    drag.label.arrow = drag.label.arrow,
    chroms.drop = chroms.drop,
    chroms_auto = chroms_auto,
    chroms_req = chroms_req,
    drag.label.linewidth = drag.label.linewidth,
    lead.stroke = lead.stroke,
    lead.label.size = lead.label.size,
    lead.label.angle = lead.label.angle,
    label_angle_norm = label_angle_norm,
    is_angle_90 = is_angle_90,
    is_angle_270 = is_angle_270,
    is_angle_zero = is_angle_zero,
    is_angle_repel_safe = is_angle_repel_safe,
    label_hjust = label_hjust,
    label_vjust = label_vjust,
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
    lead_color_enabled = lead_color_enabled,
    lead_col0 = lead_col0
  )
}

.manhattan_prepare_input <- function(data,
                                             snp.col,
                                             chrom.col,
                                             pos.col,
                                             p.col,
                                             y.col,
                                             lead.label.col,
                                             direction,
                                             beta.col,
                                             a1.col,
                                             a2.col,
                                             direction.a1) {
  snp_col_use <- .gcanvas_resolve_colname(names(data), snp.col, aliases = character(), required = TRUE, arg_label = "snp.col")
  chrom_col_use <- .gcanvas_resolve_colname(names(data), chrom.col, aliases = character(), required = TRUE, arg_label = "chrom.col")
  pos_col_use <- .gcanvas_resolve_colname(names(data), pos.col, aliases = character(), required = TRUE, arg_label = "pos.col")
  cols <- c(snp_col_use, chrom_col_use, pos_col_use)
  if (is.null(y.col)) cols <- c(cols, p.col) else cols <- c(cols, y.col)
  if (!is.null(lead.label.col) && nzchar(as.character(lead.label.col)[1])) cols <- c(cols, as.character(lead.label.col)[1])
  if (isTRUE(direction) && !is.null(beta.col) && nzchar(as.character(beta.col)[1])) cols <- c(cols, as.character(beta.col)[1])
  if (!is.null(a1.col) && nzchar(as.character(a1.col)[1])) cols <- c(cols, as.character(a1.col)[1])
  if (!is.null(a2.col) && nzchar(as.character(a2.col)[1])) cols <- c(cols, as.character(a2.col)[1])
  if (is.character(direction.a1) && length(direction.a1) == 1L && !is.na(direction.a1) && (direction.a1 %in% names(data))) {
    cols <- c(cols, as.character(direction.a1)[1])
  }
  cols <- unique(cols)
  miss <- setdiff(cols, names(data))
  if (length(miss)) stop("Missing columns in data: ", paste(miss, collapse = ", "), call. = FALSE)

  dt <- if (data.table::is.data.table(data)) {
    data.table::copy(data[, ..cols])
  } else {
    data.table::as.data.table(data[, cols, drop = FALSE])
  }
  data.table::setnames(dt, snp_col_use, "snp")
  data.table::setnames(dt, chrom_col_use, "CHR")
  data.table::setnames(dt, pos_col_use, "POS")
  if (is.null(y.col)) {
    data.table::setnames(dt, p.col, "P")
  } else {
    data.table::setnames(dt, y.col, "Y")
  }
  if (isTRUE(direction) && !is.null(beta.col) && nzchar(as.character(beta.col)[1]) && (beta.col %in% names(dt))) {
    data.table::setnames(dt, beta.col, "beta")
  }

  allele_a1_col <- a1.col
  allele_a2_col <- a2.col
  if (!is.null(allele_a1_col) && !is.null(allele_a2_col)) {
    if (!(allele_a1_col %in% names(data)) || !(allele_a2_col %in% names(data))) {
      stop("a1.col/a2.col not found in data.", call. = FALSE)
    }
    if (allele_a1_col %in% names(dt)) data.table::setnames(dt, allele_a1_col, "A1")
    if (allele_a2_col %in% names(dt)) data.table::setnames(dt, allele_a2_col, "A2")
  }
  if (is.character(direction.a1) && length(direction.a1) == 1L && !is.na(direction.a1) && (direction.a1 %in% names(dt))) {
    data.table::setnames(dt, direction.a1, "direction_a1_input")
  }

  lead_label_col_std <- NULL
  if (!is.null(lead.label.col) && nzchar(as.character(lead.label.col)[1])) {
    c0 <- as.character(lead.label.col)[1]
    if (identical(c0, snp_col_use)) {
      lead_label_col_std <- "snp"
    } else if (identical(c0, chrom_col_use)) {
      lead_label_col_std <- "CHR"
    } else if (identical(c0, pos_col_use)) {
      lead_label_col_std <- "POS"
    } else if (is.null(y.col) && identical(c0, p.col)) {
      lead_label_col_std <- "P"
    } else if (!is.null(y.col) && identical(c0, y.col)) {
      lead_label_col_std <- "Y"
    } else if (c0 %in% names(dt)) {
      data.table::setnames(dt, c0, "lead_label_src")
      lead_label_col_std <- "lead_label_src"
    }
  }

  dt[, snp := as.character(snp)]
  dt[, CHR := normalize.chrom(CHR)]
  dt[, POS := suppressWarnings(as.numeric(POS))]
  if (is.null(y.col)) {
    dt[, P := .gcanvas_p_filter(P)]
    dt <- dt[!is.na(P) & nzchar(P) & is.finite(POS) & !is.na(CHR) & nzchar(CHR)]
    dt[, yval := -log10c(P)]
  } else {
    dt[, Y := suppressWarnings(as.numeric(Y))]
    dt <- dt[is.finite(Y) & is.finite(POS) & !is.na(CHR) & nzchar(CHR)]
    dt[, yval := as.numeric(Y)]
  }
  if (!nrow(dt)) stop("No valid rows after parsing.", call. = FALSE)

  if (isTRUE(direction) && ("beta" %in% names(dt))) {
    dt[, beta := suppressWarnings(as.numeric(beta))]
    beta_aligned <- as.numeric(dt$beta)
    beta_aligned[!is.finite(beta_aligned)] <- NA_real_
    if (!is.null(direction.a1)) {
      dir_obj <- .gcanvas_resolve_direction_a1(direction.a1, df = data, snp.col = snp_col_use, a1.col = allele_a1_col %||% "A1")
      dir_ref <- rep(NA_character_, nrow(dt))
      if (identical(dir_obj$type, "column")) {
        dir_ref <- if ("direction_a1_input" %in% names(dt)) toupper(as.character(dt$direction_a1_input)) else rep(NA_character_, nrow(dt))
      } else if (identical(dir_obj$type, "map")) {
        dir_ref <- unname(dir_obj$map[dt$snp])
      }
      have_a1a2 <- ("A1" %in% names(dt)) && ("A2" %in% names(dt))
      if (!have_a1a2 && !identical(dir_obj$type, "column")) {
        stop("direction.a1 requires a1.col and a2.col unless direction.a1 is a data column name.", call. = FALSE)
      }
      if (have_a1a2) {
        dt[, A1 := toupper(as.character(A1))]
        dt[, A2 := toupper(as.character(A2))]
        m_a1 <- !is.na(dir_ref) & nzchar(dir_ref) & !is.na(dt$A1) & (dir_ref == dt$A1)
        m_a2 <- !is.na(dir_ref) & nzchar(dir_ref) & !is.na(dt$A2) & (dir_ref == dt$A2)
        beta_aligned[m_a2] <- -1 * beta_aligned[m_a2]
        aligned_ok <- m_a1 | m_a2
        beta_aligned[!aligned_ok] <- 0
      } else {
        has_ref <- !is.na(dir_ref) & nzchar(dir_ref)
        beta_aligned[!has_ref] <- 0
      }
    }
    dt[, beta_aligned := beta_aligned]
    dt[, direction_sign := sign(beta_aligned)]
    dt[!is.finite(direction_sign) | is.na(direction_sign), direction_sign := 0]
  }

  chr_now <- .gcanvas_sort_chr_unique(dt$CHR)
  auto_chr <- as.character(1:22)
  chr_extra <- c("X", "Y", "MT")
  chr_summary <- if (all(auto_chr %in% chr_now)) {
    extra <- chr_extra[chr_extra %in% chr_now]
    if (length(extra)) paste0("1-22,", paste(extra, collapse = ",")) else "1-22"
  } else if (length(chr_now) <= 12L) {
    paste(chr_now, collapse = ",")
  } else {
    sprintf("n=%d", as_int(length(chr_now)))
  }
  if (is.null(y.col)) {
    min_idx <- suppressWarnings(which.max(dt$yval))[1]
    summary_msg <- sprintf(
      "n_variants=%d | chromosomes=%s | minP=%s",
      as_int(nrow(dt)),
      chr_summary,
      if (is.finite(min_idx) && !is.na(min_idx) && min_idx >= 1L) .gcanvas_format_minp(dt$P[min_idx], digits = 3L, cutoff = 1e-3) else "NA"
    )
  } else {
    summary_msg <- sprintf("n_variants=%d | chromosomes=%s | maxY=%.4g", as_int(nrow(dt)), chr_summary, max(dt$yval, na.rm = TRUE))
  }

  list(
    dt = dt,
    snp_col_use = snp_col_use,
    chrom_col_use = chrom_col_use,
    pos_col_use = pos_col_use,
    allele_a1_col = allele_a1_col,
    allele_a2_col = allele_a2_col,
    lead_label_col_std = lead_label_col_std,
    chr_summary = chr_summary,
    summary_msg = summary_msg
  )
}

.manhattan_prepare_layout <- function(dt,
                                              build,
                                              chroms_auto,
                                              chroms.drop,
                                              chroms_req,
                                              threshold,
                                              y.col,
                                              threshold.color,
                                              threshold.type,
                                              threshold.linewidth,
                                              y.rescale.at,
                                              y.rescale.ratio,
                                              y.rescale.breaks,
                                              y.ceiling,
                                              dark,
                                              chrom.color,
                                              alpha,
                                              silent) {
  y_rescale_at <- suppressWarnings(as.numeric(y.rescale.at))[1]
  if (!is.finite(y_rescale_at) || is.na(y_rescale_at) || y_rescale_at <= 0) y_rescale_at <- NA_real_
  y_rescale_ratio <- suppressWarnings(as.numeric(y.rescale.ratio))[1]
  if (!is.finite(y_rescale_ratio) || is.na(y_rescale_ratio) || y_rescale_ratio <= 0 || y_rescale_ratio > 1) y_rescale_ratio <- 0.25
  y_rescale_breaks <- suppressWarnings(as.numeric(y.rescale.breaks))
  y_rescale_breaks <- y_rescale_breaks[is.finite(y_rescale_breaks)]
  threshold_num_in <- suppressWarnings(as.numeric(threshold))
  threshold_num_in <- threshold_num_in[is.finite(threshold_num_in) & !is.na(threshold_num_in)]
  threshold_zero_mode <- (length(threshold_num_in) == 1L && abs(threshold_num_in[1]) < .Machine$double.eps)
  use_y_rescale <- is.finite(y_rescale_at) && !is.na(y_rescale_at)
  .manhattan_y_map <- function(v) {
    vv <- suppressWarnings(as.numeric(v))
    if (!isTRUE(use_y_rescale)) return(vv)
    out <- vv
    idx <- is.finite(vv) & vv > y_rescale_at
    out[idx] <- y_rescale_at + (vv[idx] - y_rescale_at) * y_rescale_ratio
    out
  }
  y_ceiling0 <- suppressWarnings(as.numeric(y.ceiling))[1]
  use_ceiling <- is.finite(y_ceiling0)
  if (use_ceiling) {
    dt[, ceiling := as_int(yval > y_ceiling0)]
    dt[, ydraw := pmin(yval, y_ceiling0)]
  } else {
    dt[, ceiling := 0L]
    dt[, ydraw := yval]
  }
  dt[, ydraw_raw := ydraw]
  dt[, ydraw := .manhattan_y_map(ydraw_raw)]
  threshold_lines <- .gcanvas_threshold_lines(
    threshold = threshold,
    y.col = y.col,
    threshold.color = threshold.color,
    threshold.type = threshold.type,
    threshold.linewidth = threshold.linewidth,
    map_fun = .manhattan_y_map,
    default_color = "grey20"
  )
  if (isTRUE(threshold_zero_mode)) threshold_lines <- threshold_lines[0]
  chr_obs <- dt[, .(obs_min = min(POS, na.rm = TRUE), obs_max = max(POS, na.rm = TRUE)), by = CHR]
  chr_obs[, chr_order := rank.chrom(CHR)]
  data.table::setorder(chr_obs, chr_order, CHR)
  chr_ref <- .gcanvas_hg_chr_bounds(build = build)
  chr_obs_vals <- as.character(chr_obs$CHR)
  chr_obs_vals <- chr_obs_vals[!is.na(chr_obs_vals) & nzchar(chr_obs_vals)]
  auto_chr_keep <- as.character(1:22)
  chr_keep <- if (isTRUE(chroms_auto)) {
    if (isTRUE(chroms.drop)) unique(chr_obs_vals) else unique(c(auto_chr_keep, chr_obs_vals))
  } else {
    unique(c(chr_obs_vals, chroms_req))
  }
  chr_keep <- chr_keep[!is.na(chr_keep) & nzchar(chr_keep)]
  chr_keep <- unique(chr_keep)
  chr_keep <- chr_keep[order(rank.chrom(chr_keep), chr_keep)]
  if (!length(chr_keep)) stop("No chromosomes available after applying chroms/chroms.drop.", call. = FALSE)
  chr_keep_dt <- data.table::data.table(CHR = chr_keep)
  chr_map <- merge(chr_keep_dt, chr_ref, by = "CHR", all.x = TRUE, sort = FALSE)
  chr_map <- merge(chr_map, chr_obs[, .(CHR, obs_min, obs_max)], by = "CHR", all.x = TRUE, sort = FALSE)
  chr_map[, chr_order := rank.chrom(CHR)]
  data.table::setorderv(chr_map, c("chr_order", "CHR"), c(1L, 1L), na.last = TRUE)
  drop_unresolved <- chr_map[
    ((!is.finite(start) | is.na(start)) & (!is.finite(obs_min) | is.na(obs_min))) |
      ((!is.finite(end) | is.na(end)) & (!is.finite(obs_max) | is.na(obs_max))),
    as.character(CHR)
  ]
  drop_unresolved <- unique(drop_unresolved[!is.na(drop_unresolved) & nzchar(drop_unresolved)])
  if (length(drop_unresolved)) {
    .gcanvas_warn_msg(sprintf(
      "Skipping unresolved chromosomes (no reference bounds and no observed POS): %s",
      paste(drop_unresolved, collapse = ",")
    ))
    chr_map <- chr_map[!(CHR %in% drop_unresolved)]
  }
  if (!nrow(chr_map)) stop("No chromosomes available after resolving chromosome bounds.", call. = FALSE)
  chr_map[, chr_start := data.table::fifelse(
    is.finite(obs_min),
    data.table::fifelse(is.finite(start), pmin(start, obs_min), obs_min),
    start
  )]
  chr_map[!is.finite(chr_start) | is.na(chr_start), chr_start := 1]
  chr_map[, chr_start := pmax(1, chr_start)]
  chr_map[, chr_end := data.table::fifelse(
    is.finite(obs_max),
    data.table::fifelse(is.finite(end), pmax(end, obs_max), obs_max),
    end
  )]
  chr_map[!is.finite(chr_end) | is.na(chr_end), chr_end := chr_start]
  chr_map[chr_end < chr_start, chr_end := chr_start]
  chr_map[, chr_width := pmax(1, chr_end - chr_start + 1)]
  chr_map[, chr_offset := data.table::shift(cumsum(chr_width), fill = 0)]
  chr_map[, tick := chr_offset + chr_width / 2]
  chrom_palette_key <- NA_character_
  if (isTRUE(dark)) {
    chrom_cols <- c("grey25", "grey50")
  } else {
    chrom_cols <- as.character(chrom.color)
    if (length(chrom_cols) == 1L && !is.na(chrom_cols) && nzchar(chrom_cols)) {
      key <- tolower(trimws(chrom_cols[1]))
      pal_list <- .gcanvas_chrom_palette_list()
      if (key %in% names(pal_list)) {
        chrom_cols <- as.character(pal_list[[key]])
        chrom_palette_key <- key
      }
    }
  }
  chrom_cols <- chrom_cols[!is.na(chrom_cols) & nzchar(chrom_cols)]
  if (!length(chrom_cols)) chrom_cols <- c("grey60", "grey80")
  if (identical(chrom_palette_key, "random")) {
    .gcanvas_note("gcanvas::manhattan", sprintf("chrom.color=random -> %s", paste(chrom_cols, collapse = ", ")), silent = silent)
  }
  chr_map[, chr_plot_color := rep(chrom_cols, length.out = .N)]
  dt <- chr_map[dt, on = .(CHR)]
  dt[, x := chr_offset + (POS - chr_start + 1)]
  dt[, .row_id := .I]
  dt[, `:=`(
    plot_color = chr_plot_color,
    plot_alpha = alpha,
    flank_highlight = FALSE,
    lead_highlight = FALSE,
    pass_highlight = FALSE,
    lead_border_color = "grey20"
  )]
  line_dt <- data.table::copy(threshold_lines)
  if (nrow(line_dt)) {
    data.table::setorderv(line_dt, c("idx"), c(1L), na.last = TRUE)
    line_dt <- line_dt[, .(idx, y_raw, y, color, type, linewidth)]
  } else {
    line_dt <- data.table::data.table(
      idx = integer(),
      y_raw = numeric(),
      y = numeric(),
      color = character(),
      type = character(),
      linewidth = numeric()
    )
  }
  list(
    dt = dt,
    y_rescale_at = y_rescale_at,
    y_rescale_ratio = y_rescale_ratio,
    y_rescale_breaks = y_rescale_breaks,
    threshold_zero_mode = threshold_zero_mode,
    use_y_rescale = use_y_rescale,
    .manhattan_y_map = .manhattan_y_map,
    y_ceiling0 = y_ceiling0,
    use_ceiling = use_ceiling,
    threshold_lines = threshold_lines,
    line_dt = line_dt,
    chr_map = chr_map
  )
}

.manhattan_mark_ceiling_candidates <- function(dt, lead_tbl, use_ceiling) {
  dt[, ceiling_marker := FALSE]
  if (!isTRUE(use_ceiling)) return(dt)
  .manhattan_pick_ceiling_chr <- function(pos, score, flank_bp) {
    n <- length(pos)
    if (!n) return(integer())
    ord <- order(-score, pos, na.last = TRUE)
    picked <- integer()
    sel_pos <- numeric()
    for (ii in ord) {
      x <- pos[ii]
      if (!length(sel_pos)) {
        picked <- c(picked, ii)
        sel_pos <- x
        next
      }
      sp <- sort(sel_pos)
      k <- findInterval(x, sp)
      d <- Inf
      if (k > 0L) d <- min(d, abs(x - sp[k]))
      if (k < length(sp)) d <- min(d, abs(sp[k + 1L] - x))
      if (d > flank_bp) {
        picked <- c(picked, ii)
        sel_pos <- c(sel_pos, x)
      }
    }
    picked
  }

  ceiling_pool <- dt[ceiling == 1L & is.finite(POS) & !is.na(CHR) & nzchar(CHR)]
  keep_rows <- integer()
  if (nrow(ceiling_pool)) {
    pick_rows <- ceiling_pool[, {
      kk <- .manhattan_pick_ceiling_chr(pos = POS, score = yval, flank_bp = 1e6)
      .(.row_id = .row_id[kk])
    }, by = .(CHR)]
    keep_rows <- unique(as_int(pick_rows$.row_id))
  }
  if (nrow(lead_tbl)) {
    force_rows <- integer()
    rid <- as_int(lead_tbl$.row_id_data)
    rid <- rid[is.finite(rid) & !is.na(rid)]
    if (length(rid)) force_rows <- c(force_rows, rid)
    lead_snp_keep <- unique(lead_tbl$snp[!is.na(lead_tbl$snp) & nzchar(lead_tbl$snp)])
    if (length(lead_snp_keep)) {
      idx_snp <- dt[snp %in% lead_snp_keep, which = TRUE]
      force_rows <- c(force_rows, idx_snp[!is.na(idx_snp)])
    }
    lead_pos_keep <- unique(lead_tbl[is.finite(POS) & !is.na(CHR) & nzchar(CHR), .(CHR, POS)])
    if (nrow(lead_pos_keep)) {
      idx_pos <- dt[lead_pos_keep, on = .(CHR, POS), which = TRUE]
      force_rows <- c(force_rows, idx_pos[!is.na(idx_pos)])
    }
    if (length(force_rows)) keep_rows <- c(keep_rows, as_int(force_rows))
  }
  keep_rows <- unique(keep_rows[is.finite(keep_rows) & !is.na(keep_rows)])
  if (length(keep_rows)) {
    ceiling_rows <- dt[ceiling == 1L, as_int(.row_id)]
    keep_rows <- intersect(keep_rows, ceiling_rows)
    if (length(keep_rows)) dt[.row_id %in% keep_rows, ceiling_marker := TRUE]
  }
  dt
}

.gcanvas_lambda_gc <- function(p) {
  p <- .gcanvas_p_to_num(p)
  p <- p[is.finite(p) & !is.na(p) & p > 0]
  if (!length(p)) return(NA_real_)
  chisq <- suppressWarnings(stats::qchisq(1 - p, df = 1))
  chisq <- chisq[is.finite(chisq) & !is.na(chisq)]
  if (!length(chisq)) return(NA_real_)
  stats::median(chisq) / stats::qchisq(0.5, df = 1)
}

.gcanvas_p_to_num <- function(p) {
  p_chr <- .gcanvas_p_filter(p)
  p_num <- suppressWarnings(as.numeric(p_chr))
  p_num[!is.finite(p_num) | is.na(p_num) | p_num <= 0 | p_num > 1] <- NA_real_
  p_num
}

.gcanvas_format_minp <- function(p, digits = 3L, cutoff = 1e-3) {
  d <- as_int(digits)[1]
  if (!is.finite(d) || is.na(d) || d < 1L) d <- 3L
  c0 <- as_num(cutoff)[1]
  if (!is.finite(c0) || is.na(c0) || c0 <= 0 || c0 >= 1) c0 <- 1e-3
  p_num <- .gcanvas_p_to_num(p)[1]
  if (!is.finite(p_num) || is.na(p_num)) return("NA")
  p_sig <- suppressWarnings(signif(p_num, d))
  if (!is.finite(p_sig) || is.na(p_sig) || p_sig <= 0) return("NA")
  if (p_num < c0) {
    return(formatC(p_sig, format = "e", digits = max(0L, d - 1L)))
  }
  exp10 <- floor(log10(abs(p_sig)))
  dec <- max(0L, as_int(d - exp10 - 1L))
  formatC(p_sig, format = "f", digits = dec)
}

.gcanvas_p_filter <- function(p) {
  p_chr <- trimws(as.character(p))
  p_chr[is.na(p_chr) | !nzchar(p_chr)] <- NA_character_
  if (!length(p_chr)) return(p_chr)
  lp <- suppressWarnings(log10c(p_chr))
  ok <- is.finite(lp) & !is.na(lp) & lp <= 0
  out <- rep(NA_character_, length(p_chr))
  out[ok] <- p_chr[ok]
  out
}

.gcanvas_parse_bp_span <- function(x, arg_name = "value") {
  if (is.null(x) || length(x) == 0L) return(NA_real_)
  if (is.logical(x) && length(x) == 1L && !isTRUE(x)) return(NA_real_)

  x0 <- x[1]
  if (is.numeric(x0)) {
    v <- suppressWarnings(as.numeric(x0))
    if (!is.finite(v) || is.na(v) || v <= 0) return(NA_real_)
    return(v)
  }

  s <- tolower(trimws(as.character(x0)))
  if (!nzchar(s) || s %in% c("na", "nan", "null", "false", "f", "no", "none")) return(NA_real_)
  s <- gsub(",", "", s, fixed = TRUE)
  s <- gsub("\\s+", "", s, perl = TRUE)

  rx <- regexec("^([0-9]*\\.?[0-9]+)(bp|kb|mb)?$", s, perl = TRUE)
  mt <- regmatches(s, rx)[[1]]
  if (!length(mt)) {
    stop(sprintf("%s must be numeric (bp) or string like '250kb'/'1mb'.", arg_name), call. = FALSE)
  }
  val <- suppressWarnings(as.numeric(mt[2]))
  if (!is.finite(val) || is.na(val) || val <= 0) return(NA_real_)
  unit <- mt[3]
  mult <- switch(unit, bp = 1, kb = 1e3, mb = 1e6, 1)
  as.numeric(val * mult)
}

.gcanvas_threshold_lines <- function(threshold,
                                     y.col = NULL,
                                     threshold.color = NULL,
                                     threshold.type = NULL,
                                     threshold.linewidth = 0.7,
                                     map_fun = NULL,
                                     default_color = "grey20") {
  th_in <- suppressWarnings(as.numeric(threshold))
  if (!length(th_in)) {
    return(data.table::data.table(
      idx = integer(), threshold_input = numeric(),
      y_raw = numeric(), y = numeric(),
      color = character(), type = character(), linewidth = numeric()
    ))
  }
  idx_all <- seq_along(th_in)
  ok <- is.finite(th_in) & !is.na(th_in)
  if (!any(ok)) {
    return(data.table::data.table(
      idx = integer(), threshold_input = numeric(),
      y_raw = numeric(), y = numeric(),
      color = character(), type = character(), linewidth = numeric()
    ))
  }
  y_raw_all <- th_in
  if (is.null(y.col)) {
    i_p <- is.finite(y_raw_all) & !is.na(y_raw_all) & y_raw_all > 0 & y_raw_all <= 1
    y_raw_all[i_p] <- -log10(y_raw_all[i_p])
  }
  dt <- data.table::data.table(
    idx = as_int(idx_all[ok]),
    threshold_input = as.numeric(th_in[ok]),
    y_raw = as.numeric(y_raw_all[ok])
  )
  dt <- dt[is.finite(y_raw) & !is.na(y_raw)]
  if (!nrow(dt)) {
    return(data.table::data.table(
      idx = integer(), threshold_input = numeric(),
      y_raw = numeric(), y = numeric(),
      color = character(), type = character(), linewidth = numeric()
    ))
  }

  colv <- as.character(threshold.color)
  colv <- colv[!is.na(colv) & nzchar(colv)]
  if (!length(colv)) colv <- as.character(default_color)[1]
  col_all <- rep_len(colv, length(th_in))
  dt[, color := as.character(col_all[idx])]
  dt[is.na(color) | !nzchar(color), color := as.character(default_color)[1]]

  y_max <- max(dt$y_raw, na.rm = TRUE)
  type_default_all <- rep("dotted", length(th_in))
  top_idx <- dt$idx[is.finite(dt$y_raw) & dt$y_raw == y_max]
  if (length(top_idx)) type_default_all[top_idx] <- "dashed"

  typv <- as.character(threshold.type)
  typv <- typv[!is.na(typv) & nzchar(typv)]
  typv <- typv[!tolower(trimws(typv)) %in% "auto"]
  if (length(typv)) {
    type_all <- rep_len(typv, length(th_in))
  } else {
    type_all <- type_default_all
  }
  dt[, type := as.character(type_all[idx])]
  dt[, type := tolower(trimws(type))]
  dt[type %in% c("dottd", "dot", "dott"), type := "dotted"]
  dt[is.na(type) | !nzchar(type), type := "dotted"]

  lwv <- suppressWarnings(as.numeric(threshold.linewidth))
  lwv <- lwv[is.finite(lwv) & !is.na(lwv) & lwv > 0]
  if (!length(lwv)) lwv <- 0.7
  lw_all <- rep_len(lwv, length(th_in))
  dt[, linewidth := as.numeric(lw_all[idx])]
  dt[!is.finite(linewidth) | is.na(linewidth) | linewidth <= 0, linewidth := as.numeric(lwv[1])]

  if (is.null(map_fun)) {
    dt[, y := as.numeric(y_raw)]
  } else {
    dt[, y := suppressWarnings(as.numeric(map_fun(y_raw)))]
  }
  dt <- dt[is.finite(y) & !is.na(y)]
  dt[]
}

