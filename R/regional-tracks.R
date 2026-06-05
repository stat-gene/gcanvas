# Gene/exon track layout and the user-facing `regional.track()` builder used
# by `regional()` and callable standalone for custom regional figures.

# Load a package-bundled annotation track object (`tracks.b37` / `tracks.b38`)
# by genome build. Works both when installed and under devtools::load_all().
.gcanvas_bundled_tracks <- function(build = 38L) {
  build <- as_int(build)[1]
  if (is.na(build)) build <- 38L
  if (build == 19L) build <- 37L
  if (!(build %in% c(37L, 38L))) {
    stop("Bundled tracks are only available for build 37 or 38.", call. = FALSE)
  }
  nm <- sprintf("tracks.b%s", build)
  env <- new.env(parent = emptyenv())
  ok <- tryCatch({
    utils::data(list = nm, package = "gcanvas", envir = env)
    TRUE
  }, error = function(e) FALSE, warning = function(w) FALSE)
  if (!isTRUE(ok) || !exists(nm, envir = env, inherits = FALSE)) {
    stop(sprintf("Bundled annotation tracks '%s' not found in package 'gcanvas'.", nm), call. = FALSE)
  }
  get(nm, envir = env, inherits = FALSE)
}

# Query an in-memory tracks object (same schema as gtf2rds()/tracks.b3x) for a
# region, returning the same shape as gtf_query_gene_exon().
.gcanvas_tracks_query_gene_exon <- function(tracks, chrom, start, end) {
  require_pkg("data.table")
  empty_gene <- data.table::data.table(
    gene_name = character(), gene_id = character(), biotype = character(), biotype_raw = character(),
    strand = character(), start = numeric(), end = numeric()
  )
  empty_exon <- data.table::data.table(
    gene_name = character(), gene_id = character(), transcript_id = character(),
    exon_number = character(), strand = character(), start = numeric(), end = numeric()
  )
  chrom <- normalize.chrom(chrom)[1]
  start <- as_int(start); end <- as_int(end)
  if (is.na(start) || is.na(end) || end < start) stop("Invalid start/end.", call. = FALSE)

  if (!is.list(tracks) || is.null(tracks$gene)) {
    stop("In-memory tracks object must be a list with a $gene table.", call. = FALSE)
  }
  gene <- data.table::as.data.table(tracks$gene)
  exon <- if (!is.null(tracks$exon)) data.table::as.data.table(tracks$exon) else data.table::copy(empty_exon)

  gene_chr <- normalize.chrom(gene[["CHR"]])
  gsel <- !is.na(gene_chr) & gene_chr == chrom &
    gene[["end"]] >= start & gene[["start"]] <= end
  gene <- gene[gsel]
  if (nrow(gene)) {
    biotype_raw <- if ("biotype_raw" %in% names(gene)) gene[["biotype_raw"]] else gene[["biotype"]]
    gene <- data.table::data.table(
      gene_name = as.character(gene[["gene_name"]]),
      gene_id = as.character(gene[["gene_id"]]),
      biotype = as.character(gene[["biotype"]]),
      biotype_raw = as.character(biotype_raw),
      strand = as.character(gene[["strand"]]),
      start = as.numeric(gene[["start"]]),
      end = as.numeric(gene[["end"]])
    )
  } else {
    gene <- data.table::copy(empty_gene)
  }

  if (nrow(exon) && "CHR" %in% names(exon)) {
    exon_chr <- normalize.chrom(exon[["CHR"]])
    esel <- !is.na(exon_chr) & exon_chr == chrom &
      exon[["end"]] >= start & exon[["start"]] <= end
    exon <- exon[esel]
    if (nrow(exon)) {
      exon <- data.table::data.table(
        gene_name = as.character(exon[["gene_name"]]),
        gene_id = as.character(exon[["gene_id"]]),
        transcript_id = as.character(exon[["transcript_id"]]),
        exon_number = as.character(exon[["exon_number"]]),
        strand = as.character(exon[["strand"]]),
        start = as.numeric(exon[["start"]]),
        end = as.numeric(exon[["end"]])
      )
    } else {
      exon <- data.table::copy(empty_exon)
    }
  } else {
    exon <- data.table::copy(empty_exon)
  }

  list(gene = gene, exon = exon, chrom = chrom)
}

.gcanvas_calc_label_space <- function(label, x_center, pos.range,
                              unit.text = 0.025,
                              margin_frac = 0.01) {
  pos.range <- as.numeric(pos.range)
  span <- diff(pos.range)
  if (!is.finite(span) || span <= 0) span <- 1

  margin <- span * margin_frac

  # Conservative label width estimate:
  # previous formula tended to overestimate text span and blocked upper rows too aggressively.
  label0 <- as.character(label)
  if (is.na(label0) || !nzchar(label0)) label0 <- ""
  label_core <- gsub("[\u2190\u2192\u2026]", "", label0, perl = TRUE)
  n_all <- nchar(label0, type = "width")
  n_core <- nchar(label_core, type = "width")
  arrow_bonus <- max(0, n_all - n_core) * 0.35
  eff_chars <- max(1, n_core + arrow_bonus)
  htext <- (eff_chars + 1) * unit.text * span * 0.65
  if (!is.finite(htext) || htext < 0) htext <- 0

  pick_space <- function(x, hj) {
    if (hj <= 0) return(c(x, x + htext))
    if (hj >= 1) return(c(x - htext, x))
    c(x - htext/2, x + htext/2)
  }

  x0 <- as.numeric(x_center)[1]
  if (!is.finite(x0)) x0 <- mean(pos.range)

  hjust <- 0.5
  if (x0 - htext/2 < pos.range[1] + margin) hjust <- 0
  if (x0 + htext/2 > pos.range[2] - margin) hjust <- 1

  if (hjust <= 0) {
    lo <- pos.range[1] + margin
    hi <- pos.range[2] - margin - htext
  } else if (hjust >= 1) {
    lo <- pos.range[1] + margin + htext
    hi <- pos.range[2] - margin
  } else {
    lo <- pos.range[1] + margin + htext/2
    hi <- pos.range[2] - margin - htext/2
  }

  if (is.finite(lo) && is.finite(hi) && hi >= lo) {
    x0 <- min(max(x0, lo), hi)
  } else {
    x0 <- min(max(x0, pos.range[1] + margin), pos.range[2] - margin)
  }

  sp <- pick_space(x0, hjust)
  list(text_x = x0, tstart = sp[1], tend = sp[2], hjust = hjust, htext = htext)
}


layout_gene_tracks <- function(gene_df, exon_df, pos.range, lead_pos,
                               ystep,
                               unit.text = 0.025,
                               max_levels = 6L,
                               max_genes = 50L,
                               force_genes = character(),
                               force_biotype_override = character(),
                               force_label_map = NULL,
                               force_ignore_max_levels = TRUE,
                               biotype_priority = c("protein_coding","lncRNA","miRNA"),
                               biotype_keep = c("protein_coding","lncRNA"),
                               lncrna_symbol_only = FALSE) {
  require_pkg("data.table")
  # --- gene data.table ---
  gene_dt <- data.table::as.data.table(gene_df)
  # sanitize gene identifiers early (prevents NULL/list poisoning)
  gene_dt[, gene_name := .gcanvas_as_char_no_null(gene_name, empty = "")]
  if ("gene_id" %in% names(gene_dt)) {
    gene_dt[, gene_id := .gcanvas_as_char_no_null(gene_id, empty = "")]
  } else {
    gene_dt[, gene_id := ""]
  }
  gene_dt[, gene_name := .gcanvas_normalize_ensembl_gene_id(gene_name)]
  gene_dt[, gene_id := .gcanvas_normalize_ensembl_gene_id(gene_id)]
  # if gene_name missing but gene_id exists, use gene_id as name
  gene_dt[!nzchar(gene_name) & nzchar(gene_id), gene_name := gene_id]
  # drop genes that have no usable identifier -> bar/exon/label should disappear
  gene_dt <- gene_dt[nzchar(gene_name) | nzchar(gene_id)]

  # --- exon data.table ---
  exon_dt <- data.table::as.data.table(exon_df)
  exon_dt[, gene_name := .gcanvas_as_char_no_null(gene_name, empty = "")]
  if ("gene_id" %in% names(exon_dt)) {
    exon_dt[, gene_id := .gcanvas_as_char_no_null(gene_id, empty = "")]
  } else {
    exon_dt[, gene_id := ""]
  }
  exon_dt[, gene_name := .gcanvas_normalize_ensembl_gene_id(gene_name)]
  exon_dt[, gene_id := .gcanvas_normalize_ensembl_gene_id(gene_id)]
  exon_dt[!nzchar(gene_name) & nzchar(gene_id), gene_name := gene_id]
  exon_dt <- exon_dt[nzchar(gene_name) | nzchar(gene_id)]

  n_total_all0 <- nrow(gene_dt)
  if (!n_total_all0) {
    meta <- list(
      n_total = 0L, n_kept = 0L, n_omitted = 0L,
      n_omitted_by_max_genes = 0L, n_omitted_by_max_levels = 0L,
      n_dropped_by_biotype_keep = 0L,
      biotype_filtered_counts = integer(0),
      n_total_all = 0L
    )
    return(list(gene = gene_dt[0], exon = exon_dt[0], meta = meta))
  }

  pos.range <- .gcanvas_as_num2(pos.range)
  pos.range <- pos.range[is.finite(pos.range)]
  if (length(pos.range) != 2) stop("layout_gene_tracks: pos.range must be numeric length-2.", call. = FALSE)
  pos.range <- c(min(pos.range), max(pos.range))

  lead_pos <- suppressWarnings(as.numeric(lead_pos))
  if (!is.finite(lead_pos)) lead_pos <- mean(pos.range)

  force_genes <- unique(.gcanvas_normalize_ensembl_gene_id(as.character(force_genes)))
  force_biotype_override <- unique(.gcanvas_normalize_ensembl_gene_id(as.character(force_biotype_override)))
  gene_dt[, forced := gene_name %in% force_genes]

  if (!("biotype" %in% names(gene_dt))) gene_dt[, biotype := "other"]
  gene_dt[, biotype := as.character(biotype)]
  gene_dt[is.na(biotype) | !nzchar(biotype), biotype := "other"]

  biotype_keep <- unique(as.character(biotype_keep))
  if (length(biotype_keep) == 0L || all(!nzchar(biotype_keep))) biotype_keep <- NULL

  n_drop_keep <- 0L
  dropped_counts <- integer(0)

  if (!is.null(biotype_keep)) {
    # User-forced genes can bypass biotype filter.
    keep_mask <- (gene_dt$biotype %in% biotype_keep) | (gene_dt$gene_name %in% force_biotype_override)
    if (isTRUE(lncrna_symbol_only)) {
      is_lnc <- (!is.na(gene_dt$biotype) & gene_dt$biotype == "lncRNA")
      is_id_like <- .gcanvas_is_ensembl_gene_id(gene_dt$gene_name)
      keep_mask <- keep_mask & (!(is_lnc & is_id_like) | (gene_dt$gene_name %in% force_biotype_override))
    }
    n_drop_keep <- sum(!keep_mask, na.rm = TRUE)

    if (n_drop_keep > 0L) {
      dropped_dt <- gene_dt[!keep_mask]
      dropped_counts <- as.integer(table(dropped_dt$biotype))
      names(dropped_counts) <- names(table(dropped_dt$biotype))
    }

    gene_dt <- gene_dt[keep_mask]
  }

  n_total_sel0 <- nrow(gene_dt)

  if (!n_total_sel0) {
    meta <- list(
      n_total = 0L, n_kept = 0L, n_omitted = 0L,
      n_omitted_by_max_genes = 0L, n_omitted_by_max_levels = 0L,
      n_dropped_by_biotype_keep = as_int(n_drop_keep),
      biotype_filtered_counts = dropped_counts,
      max_levels = as_int(max_levels),
      max_genes = as_int(max_genes),
      n_forced = as_int(sum(gene_df$gene_name %in% force_genes, na.rm = TRUE)),
      biotype_priority = unique(as.character(biotype_priority)),
      biotype_keep = biotype_keep,
      n_total_all = as_int(n_total_all0)
    )
    return(list(
      gene = gene_dt[0],
      exon = exon_dt[0],
      meta = meta
    ))
  }

  biotype_priority <- unique(as.character(biotype_priority))
  if (length(biotype_priority) == 0L || all(!nzchar(biotype_priority))) biotype_priority <- character()

  if (length(biotype_priority)) {
    gene_dt[, biotype_rank := match(biotype, biotype_priority)]
    gene_dt[is.na(biotype_rank), biotype_rank := length(biotype_priority) + 1L]
  } else {
    gene_dt[, biotype_rank := 1L]
  }

  gene_dt[, mid := (start + end) / 2]
  gene_dt[, dist := pmin(abs(start - lead_pos), abs(end - lead_pos), abs(mid - lead_pos), na.rm = TRUE)]
  gene_dt[lead_pos >= start & lead_pos <= end, dist := 0]

  forced_dt <- gene_dt[forced == TRUE]
  nonforced_dt <- gene_dt[forced == FALSE][order(biotype_rank, dist, abs(end - start), na.last = TRUE)]

  max_genes <- as_int(max_genes)
  if (is.na(max_genes) || max_genes < 1L) max_genes <- nrow(nonforced_dt)
  kept_nonforced <- nonforced_dt[seq_len(min(max_genes, nrow(nonforced_dt)))]

  cand <- data.table::rbindlist(list(forced_dt, kept_nonforced), fill = TRUE)
  cand <- unique(cand, by = c("gene_name","gene_id","start","end","strand"))

  n_omitted_by_max_genes <- nrow(gene_dt) - nrow(cand)

  data.table::setorder(cand, start, end)

  n_label_omitted <- 0L
  levels_end <- numeric(0)
  cand[, `:=`(level = NA_integer_, label = NA_character_, text_x = NA_real_, hstart = NA_real_, hend = NA_real_, hjust = NA_real_)]

  span0 <- diff(pos.range)
  if (!is.finite(span0) || span0 <= 0) span0 <- 1
  gene_margin <- 0.01 * span0

  for (i in seq_len(nrow(cand))) {
    gstart <- cand$start[i] - gene_margin
    gend   <- cand$end[i] + gene_margin

    if (!is.finite(gstart)) gstart <- cand$start[i]
    if (!is.finite(gend))   gend   <- cand$end[i]

    placed <- FALSE
    if (length(levels_end)) {
      for (k in seq_along(levels_end)) {
        if (gstart >= levels_end[k]) {
          cand$level[i] <- k
          levels_end[k] <- gend
          placed <- TRUE
          break
        }
      }
    }

    if (!placed) {
      can_new <- length(levels_end) < as_int(max_levels) ||
        (isTRUE(force_ignore_max_levels) && isTRUE(cand$forced[i]))
      if (can_new) {
        levels_end <- c(levels_end, gend)
        cand$level[i] <- length(levels_end)
      }
    }
  }

  kept <- cand[!is.na(level)]

  if (nrow(kept)) {
    label_name <- .gcanvas_as_char_no_null(kept$gene_name, empty = "")
    if (!is.null(force_label_map) && length(force_label_map)) {
      map_keys <- names(force_label_map)
      idx_force <- which((kept$forced %in% TRUE) & (kept$gene_name %in% map_keys))
      if (length(idx_force)) {
        label_name[idx_force] <- unname(force_label_map[kept$gene_name[idx_force]])
      }
    }

    kept[, label := ifelse(
      strand == "+", paste0(label_name, "\u2192"),
      ifelse(strand == "-", paste0("\u2190", label_name), label_name)
    )]

    max_chars <- max(5L, as.integer(floor((1 - 2 * 0.01) / unit.text) - 2L))

    kept[, label := vapply(label, function(s) {
      if (is.na(s) || !nzchar(s)) return("")
      pre <- ""; suf <- ""
      if (startsWith(s, "\u2190")) { pre <- "\u2190"; s <- substring(s, 2) }
      if (endsWith(s, "\u2192")) { suf <- "\u2192"; s <- substring(s, 1, nchar(s)-1) }
      keep <- max_chars - nchar(pre) - nchar(suf)
      keep <- max(3L, keep)
      if (nchar(s) > keep) s <- paste0(substr(s, 1, keep-1), "\u2026")
      paste0(pre, s, suf)
    }, character(1))]

    kept[, mid := (start + end) / 2]
    tmp <- lapply(seq_len(nrow(kept)), function(i) {
      .gcanvas_calc_label_space(kept$label[i], kept$mid[i], pos.range, unit.text = unit.text)
    })

    kept[, `:=`(
      text_x = vapply(tmp, `[[`, numeric(1), "text_x"),
      tstart = vapply(tmp, `[[`, numeric(1), "tstart"),
      tend   = vapply(tmp, `[[`, numeric(1), "tend"),
      hjust  = vapply(tmp, `[[`, numeric(1), "hjust"),
      htext  = vapply(tmp, `[[`, numeric(1), "htext")
    )]

    # Left-to-right interval packing by genomic coordinates:
    # place each gene on the highest possible row (smallest level index),
    # and reuse freed rows when intervals no longer overlap.
    # Keep start/end first so row order follows actual gene position, not label arrow direction.
    data.table::setorder(kept, start, end, tstart, tend)

    kept[, label_ok := TRUE]
    kept[is.na(label) | !nzchar(label), label_ok := FALSE]
    kept[!is.finite(tstart) | !is.finite(tend), label_ok := FALSE]

    # Use gene + label span for collision packing so labels do not overlap.
    kept[, span_start := ifelse(label_ok, pmin(start, tstart), start)]
    kept[, span_end := ifelse(label_ok, pmax(end, tend), end)]

    row_end <- numeric(0)
    new_level <- integer(nrow(kept))
    for (j in seq_len(nrow(kept))) {
      s <- kept$span_start[j]
      e <- kept$span_end[j]
      if (!is.finite(s)) s <- kept$start[j]
      if (!is.finite(e)) e <- kept$end[j]
      if (!is.finite(s)) s <- -Inf
      if (!is.finite(e)) e <- s

      if (length(row_end)) {
        free <- which(s >= row_end)
      } else {
        free <- integer(0)
      }

      if (length(free)) {
        lv <- free[1]
        new_level[j] <- lv
        row_end[lv] <- e
      } else {
        row_end <- c(row_end, e)
        new_level[j] <- length(row_end)
      }
    }

    kept[, level := new_level]
    n_label_omitted <- sum(!kept$label_ok, na.rm = TRUE)
    kept[!(label_ok), label := NA_character_]

    kept[, c("label_ok","tstart","tend","htext","mid","span_start","span_end") := NULL]

  }


  n_omitted_by_max_levels <- nrow(cand) - nrow(kept)
  n_kept <- nrow(kept)

  n_omitted <- n_total_sel0 - n_kept

  if (n_kept) {
    kept[, gene_y := -1 * (level + 0.5) * ystep]
    kept[, text_y := gene_y + ystep * 0.2]
    kept[, gene_name := .gcanvas_as_char_no_null(kept$gene_name, empty = "")]
    kept[, label := .gcanvas_as_char_no_null(kept$label, empty = "")]
  }

  if (nrow(exon_dt) && n_kept) {
    exon_dt <- exon_dt[gene_name %in% kept$gene_name]
    exon_dt <- merge(exon_dt, kept[, .(gene_name, gene_y, text_x, text_y)], by = "gene_name", all.x = TRUE)
  } else {
    exon_dt <- exon_dt[0]
  }

  meta <- list(
    n_total = as_int(n_total_sel0),
    n_kept = as_int(n_kept),
    n_omitted = as_int(n_omitted),
    n_omitted_by_max_genes = as_int(max(0L, n_omitted_by_max_genes)),
    n_omitted_by_max_levels = as_int(max(0L, n_omitted_by_max_levels)),
    n_label_omitted = as_int(n_label_omitted),
    n_dropped_by_biotype_keep = as_int(n_drop_keep),
    biotype_filtered_counts = dropped_counts,
    max_levels = as_int(max_levels),
    max_genes = as_int(max_genes),
    n_forced = as_int(sum(gene_dt$gene_name %in% force_genes, na.rm = TRUE)),
    biotype_priority = biotype_priority,
    biotype_keep = biotype_keep,
    n_total_all = as_int(n_total_all0)
  )

  list(gene = kept, exon = exon_dt, meta = meta)
}


#' Build gene/exon tracks for a regional plot
#'
#' Queries a `gtf2rds()` annotation cache for the requested region, applies
#' biotype filters and gene-force / gene-add rules, lays out genes across
#' available levels, and returns the layout tables that [regional()] (or a
#' bespoke regional figure) can render below the association track.
#'
#' @param gtf_bgz Path to the bgzipped GTF (the `.rds` cache is located
#'   alongside it). Optional when `tracks` is supplied.
#' @param tracks Optional in-memory annotation object (the same list shape as
#'   produced by [gtf2rds()], e.g. the bundled [tracks.b37] / [tracks.b38]).
#'   When supplied, gene/exon records are taken from it directly and `gtf_bgz`
#'   is ignored (no tabix query is performed).
#' @param chrom Chromosome of the region.
#' @param pos.range Numeric length-2 vector `(start, end)` in bp.
#' @param y.max Y-axis maximum from the association panel (for proportioning).
#' @param lead_pos Position of the lead variant (used for label spacing).
#' @param keep_biotype Character vector of biotypes to keep (default
#'   `c("protein_coding", "lncRNA")` unless `include_other_biotypes` is set).
#' @param gene_max_levels Maximum number of vertical levels for gene layout.
#' @param gene_max_n Maximum number of genes drawn.
#' @param gene_force Symbols always drawn, even past `gene_max_n`/levels.
#' @param gene_force_label Optional label override for `gene_force`.
#' @param gene_add Additional gene metadata to overlay.
#' @param gene_force_ignore_max_levels Logical. Allow forced genes past the
#'   level cap.
#' @param biotype_priority Sort order applied when picking which genes win.
#' @param lncrna_symbol_only Logical. Restrict lncRNAs to symbol-named entries.
#' @param include_other_biotypes Logical. Allow non-default biotypes.
#' @param plot.width,units Plot-size hints used for label placement.
#' @param force_edge_genes Logical. Keep genes whose edges touch the region.
#'
#' @return A list of layout tables (gene, exon, label) consumed by the
#'   regional plotter.
#' @export
regional.track <- function(gtf_bgz = NULL, chrom, pos.range, y.max, lead_pos,
                                   tracks = NULL,
                                   keep_biotype = NULL,
                                   gene_max_levels = 6L,
                                   gene_max_n = 50L,
                                   gene_force = character(),
                                   gene_force_label = NULL,
                                   gene_add = NULL,
                                   gene_force_ignore_max_levels = TRUE,
                                   biotype_priority = c("protein_coding","lncRNA","miRNA"),
                                   lncrna_symbol_only = FALSE,
                                   include_other_biotypes = FALSE,
                                   plot.width = NULL,
                                   units = "cm",
                                   force_edge_genes = TRUE) {
  require_pkg("data.table")
  pos.range <- .gcanvas_as_num2(pos.range)
  if (length(pos.range) != 2 || any(!is.finite(pos.range))) stop("regional.track: invalid pos.range.", call. = FALSE)
  pos.range <- c(min(pos.range), max(pos.range))

  if (isFALSE(include_other_biotypes) && is.null(keep_biotype)) {
    keep_biotype <- c("protein_coding", "lncRNA")
  }
  biotype_keep_for_layout <- if (isTRUE(include_other_biotypes) && is.null(keep_biotype)) NULL else keep_biotype
  if (is.null(biotype_keep_for_layout)) {
    lncrna_symbol_only <- FALSE
  } else {
    keep_lower <- tolower(as.character(biotype_keep_for_layout))
    keep_lower <- keep_lower[!is.na(keep_lower) & nzchar(keep_lower)]
    if (!("lncrna" %in% keep_lower)) lncrna_symbol_only <- FALSE
  }

  if (!is.null(tracks)) {
    anno <- .gcanvas_tracks_query_gene_exon(
      tracks, chrom = chrom,
      start = min(pos.range), end = max(pos.range)
    )
  } else {
    if (is.null(gtf_bgz) || !length(gtf_bgz) || is.na(gtf_bgz[1]) || !nzchar(as.character(gtf_bgz[1]))) {
      stop("regional.track: either `tracks` or `gtf_bgz` must be supplied.", call. = FALSE)
    }
    anno <- gtf_query_gene_exon(
      gtf_bgz, chrom = chrom,
      start = min(pos.range), end = max(pos.range),
      features = c("gene","exon"),
      keep_biotype = NULL
    )
  }

  gene <- data.table::as.data.table(anno$gene)
  exon <- data.table::as.data.table(anno$exon)

  add_obj <- .gcanvas_extract_gene_add(gene_add)
  if (nrow(add_obj$gene)) {
    gene <- data.table::rbindlist(list(data.table::as.data.table(gene), data.table::as.data.table(add_obj$gene)), fill = TRUE)
    gene <- unique(gene, by = c("gene_name", "gene_id", "strand", "start", "end"))
  }
  if (nrow(add_obj$exon)) {
    exon <- data.table::rbindlist(list(data.table::as.data.table(exon), data.table::as.data.table(add_obj$exon)), fill = TRUE)
    exon <- unique(exon, by = c("gene_name", "gene_id", "strand", "start", "end", "transcript_id", "exon_number"))
  }

  gene <- gene[end >= pos.range[1] & start <= pos.range[2]]
  if (nrow(gene)) {
    gene[, `:=`(start = pmax(start, pos.range[1]), end = pmin(end, pos.range[2]))]
  }

  exon <- exon[end >= pos.range[1] & start <= pos.range[2]]
  if (nrow(exon)) {
    exon[, `:=`(start = pmax(start, pos.range[1]), end = pmin(end, pos.range[2]))]
  }

  gene_force_user <- unique(.gcanvas_normalize_ensembl_gene_id(as.character(gene_force)))
  gene_force <- gene_force_user
  force_label_map <- NULL
  if (!is.null(gene_force_label)) {
    gf_label <- gene_force_label
    if (is.list(gf_label)) gf_label <- unlist(gf_label, use.names = TRUE)
    gf_label <- as.character(gf_label)
    if (length(gf_label)) {
      nm <- names(gf_label)
      if (!is.null(nm) && any(!is.na(nm) & nzchar(nm))) {
        nm <- .gcanvas_normalize_ensembl_gene_id(nm)
        ok <- !is.na(nm) & nzchar(nm) & !is.na(gf_label) & nzchar(gf_label)
        if (any(ok)) {
          force_label_map <- stats::setNames(gf_label[ok], nm[ok])
          force_label_map <- force_label_map[!duplicated(names(force_label_map))]
        }
      } else if (length(gf_label) == length(gene_force)) {
        ok <- !is.na(gene_force) & nzchar(gene_force) & !is.na(gf_label) & nzchar(gf_label)
        if (any(ok)) {
          force_label_map <- stats::setNames(gf_label[ok], gene_force[ok])
          force_label_map <- force_label_map[!duplicated(names(force_label_map))]
        }
      }
    }
  }
  if (isTRUE(force_edge_genes) && nrow(gene)) {
    li <- which.min(gene$start)
    ri <- which.max(gene$end)
    edge_force <- unique(gene$gene_name[c(li, ri)])
    edge_force <- edge_force[!is.na(edge_force) & nzchar(edge_force)]
    gene_force <- unique(c(gene_force, edge_force))
  }

  ystep <- 0.15 * y.max
  unit.text <- 0.025

  plot.width <- suppressWarnings(as.numeric(plot.width))[1]
  units <- tolower(as.character(units)[1])

  if (is.finite(plot.width) && plot.width > 0) {
    if (units %in% c("in", "inch", "inches")) plot.width <- plot.width * 2.54
    baseline_cm <- 16
    unit.text <- unit.text * (baseline_cm / plot.width)
  }

  laid <- layout_gene_tracks(
    gene_df = gene, exon_df = exon,
    pos.range = pos.range, lead_pos = lead_pos,
    ystep = ystep,
    unit.text = unit.text,
    max_levels = gene_max_levels,
    max_genes = gene_max_n,
    force_genes = gene_force,
    force_biotype_override = gene_force_user,
    force_label_map = force_label_map,
    force_ignore_max_levels = gene_force_ignore_max_levels,
    biotype_priority = biotype_priority,
    biotype_keep = biotype_keep_for_layout,
    lncrna_symbol_only = lncrna_symbol_only
  )

  list(gene = laid$gene, exon = laid$exon, chrom = normalize.chrom(chrom)[1], meta = laid$meta)
}

