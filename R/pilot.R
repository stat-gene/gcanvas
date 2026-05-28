# Pilot subsampling: pull a compact, lead-preserving subset of variants
# suitable for quick exploratory plotting before running the full pipeline.

#' Downsample summary statistics for a pilot plot
#'
#' Returns a stratified subset of `data` that preserves all variants near
#' lead signals and significance peaks while thinning the rest, suitable for
#' fast iteration on plot styling before running the full dataset.
#'
#' @param data A `data.frame`/`data.table` of summary statistics.
#' @param snp.col,chrom.col,pos.col,p.col Column names in `data`.
#' @param y.col Optional alternative y-column.
#' @param lead Optional table of pre-identified lead variants.
#' @param lead.flank Window kept around each lead variant (bp).
#' @param threshold Significance threshold for automatic peak retention.
#' @param threshold.flank Window around each above-threshold variant.
#' @param n.variant Approximate total number of variants in the output.
#' @param per.chrom.min Minimum variants kept per chromosome.
#' @param method One of `"stratified"` (default) or `"uniform"`.
#' @param seed Integer seed for reproducibility.
#' @param keep.boundary Logical. Always include the chromosome-boundary rows.
#' @param keep.chr.minp Logical. Always include the per-chromosome minimum-p row.
#' @param return.type Output container: `"auto"`, `"data.table"`, or `"data.frame"`.
#' @param silent Logical. Suppress progress notes.
#'
#' @return A `data.table` of the sampled rows.
#' @export
get.pilot <- function(data,
                  snp.col = "SNP", chrom.col = "CHR", pos.col = "POS", p.col = "P",
                  y.col = NULL,
                  lead = NULL,
                  lead.flank = 5e5,
                  threshold = 5e-8,
                  threshold.flank = 5e5,
                  n.variant = 30000L,
                  per.chrom.min = 50L,
                  method = c("stratified", "uniform"),
                  seed = 23L,
                  keep.boundary = TRUE,
                  keep.chr.minp = TRUE,
                  return.type = c("auto", "data.table", "data.frame"),
                  silent = FALSE) {
  require_pkg(c("data.table"))
  if (!is.data.frame(data) && !data.table::is.data.table(data)) stop("data must be a data.frame/data.table.", call. = FALSE)
  method <- match.arg(method)
  return.type <- match.arg(return.type)
  silent <- isTRUE(silent)

  snp_col_use <- .gcanvas_resolve_colname(names(data), snp.col, aliases = character(), required = TRUE, arg_label = "snp.col")
  chrom_col_use <- .gcanvas_resolve_colname(names(data), chrom.col, aliases = character(), required = TRUE, arg_label = "chrom.col")
  pos_col_use <- .gcanvas_resolve_colname(names(data), pos.col, aliases = character(), required = TRUE, arg_label = "pos.col")
  req <- c(snp_col_use, chrom_col_use, pos_col_use, if (is.null(y.col)) p.col else y.col)
  miss <- setdiff(req, names(data))
  if (length(miss)) stop("Missing columns in data: ", paste(miss, collapse = ", "), call. = FALSE)
  .gcanvas_note("gcanvas::get.pilot", "Parsing input and scoring variants", silent = silent)

  dt <- if (data.table::is.data.table(data)) data.table::copy(data) else data.table::as.data.table(data)
  dt[, row_id := seq_len(.N)]
  dt[, snp_tmp := as.character(get(snp_col_use))]
  dt[, CHR_tmp := normalize.chrom(get(chrom_col_use))]
  dt[, POS_tmp := suppressWarnings(as.numeric(get(pos_col_use)))]
  if (is.null(y.col)) {
    dt[, P_tmp := .gcanvas_p_filter(get(p.col))]
    dt <- dt[!is.na(P_tmp) & nzchar(P_tmp) & is.finite(POS_tmp) & !is.na(CHR_tmp) & nzchar(CHR_tmp)]
    dt[, y_tmp := -log10c(P_tmp)]
    thr_y <- if (is.finite(suppressWarnings(as.numeric(threshold))[1]) && suppressWarnings(as.numeric(threshold))[1] > 0 && suppressWarnings(as.numeric(threshold))[1] <= 1) {
      -log10(suppressWarnings(as.numeric(threshold))[1])
    } else {
      suppressWarnings(as.numeric(threshold))[1]
    }
  } else {
    dt[, y_tmp := suppressWarnings(as.numeric(get(y.col)))]
    dt <- dt[is.finite(y_tmp) & is.finite(POS_tmp) & !is.na(CHR_tmp) & nzchar(CHR_tmp)]
    thr0 <- suppressWarnings(as.numeric(threshold))[1]
    thr_y <- if (is.finite(thr0) && thr0 > 0 && thr0 <= 1) -log10(thr0) else thr0
  }
  if (!nrow(dt)) stop("No valid rows after parsing.", call. = FALSE)
  dt[, chr_order := rank.chrom(CHR_tmp)]
  data.table::setorderv(dt, c("chr_order", "POS_tmp", "row_id"), c(1L, 1L, 1L), na.last = TRUE)

  lead.flank <- .gcanvas_parse_bp_span(lead.flank, arg_name = "lead.flank")
  threshold.flank <- .gcanvas_parse_bp_span(threshold.flank, arg_name = "threshold.flank")
  if (!is.finite(lead.flank) || lead.flank <= 0) lead.flank <- 5e5
  if (!is.finite(threshold.flank) || threshold.flank <= 0) threshold.flank <- 2e5

  keep <- rep(FALSE, nrow(dt))
  anchor_idx <- integer()
  .get_pilot_keep_i <- function(i) {
    i <- as_int(i)
    i <- i[is.finite(i) & !is.na(i) & i >= 1L & i <= nrow(dt)]
    if (length(i)) keep[i] <<- TRUE
    invisible(NULL)
  }
  .get_pilot_add_anchor <- function(i) {
    i <- as_int(i)
    i <- i[is.finite(i) & !is.na(i) & i >= 1L & i <= nrow(dt)]
    if (length(i)) {
      anchor_idx <<- unique(c(anchor_idx, i))
      keep[i] <<- TRUE
    }
    invisible(NULL)
  }
  .get_pilot_sample_spread <- function(idx, n_take, min_per_chr = 0L, prefer_chr = NULL) {
    idx <- unique(as_int(idx))
    idx <- idx[is.finite(idx) & !is.na(idx) & idx >= 1L & idx <= nrow(dt)]
    n_take <- as_int(n_take)
    min_per_chr <- as_int(min_per_chr)
    if (!length(idx) || !is.finite(n_take) || n_take <= 0L) return(integer())
    if (n_take >= length(idx)) return(idx)
    if (!is.finite(min_per_chr) || min_per_chr < 0L) min_per_chr <- 0L

    cand <- dt[idx, .(idx = idx, CHR_tmp, POS_tmp, chr_order)]
    data.table::setorderv(cand, c("chr_order", "POS_tmp", "idx"), c(1L, 1L, 1L), na.last = TRUE)
    cnt <- cand[, .N, by = .(CHR_tmp, chr_order)]
    data.table::setorderv(cnt, c("chr_order"), c(1L), na.last = TRUE)
    if (!is.null(prefer_chr) && length(prefer_chr)) {
      pref <- as.character(prefer_chr)
      cnt[, pref := as_int(CHR_tmp %in% pref)]
      data.table::setorderv(cnt, c("pref", "chr_order"), c(-1L, 1L), na.last = TRUE)
      cnt[, pref := NULL]
    }

    n_chr <- nrow(cnt)
    cnt[, alloc := 0L]
    cnt[, cap := N]
    if (n_chr > 0L && min_per_chr > 0L) {
      k <- min(min_per_chr, as_int(floor(n_take / n_chr)))
      if (k > 0L) cnt[, alloc := pmin(cap, k)]
    }
    left <- n_take - sum(cnt$alloc)
    if (left > 0L) {
      cap1 <- pmax(0L, cnt$cap - cnt$alloc)
      if (sum(cap1) > 0L) {
        add <- pmin(cap1, as_int(floor(left * (cap1 / sum(cap1)))))
        cnt[, alloc := alloc + add]
        left <- n_take - sum(cnt$alloc)
      }
    }
    if (left > 0L) {
      cap2 <- pmax(0L, cnt$cap - cnt$alloc)
      ord <- order(cap2, cnt$N, decreasing = TRUE)
      for (k in ord) {
        if (left <= 0L) break
        if (cap2[k] <= 0L) next
        cnt$alloc[k] <- cnt$alloc[k] + 1L
        cap2[k] <- cap2[k] - 1L
        left <- left - 1L
      }
    }
    if (sum(cnt$alloc) > n_take) {
      over <- sum(cnt$alloc) - n_take
      ord <- order(cnt$alloc, decreasing = TRUE)
      for (k in ord) {
        if (over <= 0L) break
        if (cnt$alloc[k] <= 0L) next
        dec <- min(over, cnt$alloc[k])
        cnt$alloc[k] <- cnt$alloc[k] - dec
        over <- over - dec
      }
    }

    picked <- integer()
    for (i in seq_len(nrow(cnt))) {
      n_i <- as_int(cnt$alloc[i])
      if (!is.finite(n_i) || n_i <= 0L) next
      sub <- cand[CHR_tmp == cnt$CHR_tmp[i]]
      n_i <- min(n_i, nrow(sub))
      take <- unique(as_int(round(seq(1, nrow(sub), length.out = n_i))))
      picked <- c(picked, sub$idx[take])
    }
    unique(as_int(picked))
  }

  if (isTRUE(keep.boundary)) {
    .get_pilot_add_anchor(dt[, .I[which.min(POS_tmp)], by = CHR_tmp]$V1)
    .get_pilot_add_anchor(dt[, .I[which.max(POS_tmp)], by = CHR_tmp]$V1)
  }
  if (isTRUE(keep.chr.minp)) {
    .get_pilot_add_anchor(dt[, .I[which.max(y_tmp)], by = CHR_tmp]$V1)
  }

  lead_ids <- .gcanvas_as_snp_vector(lead)
  if (!length(lead_ids)) {
    idx0 <- suppressWarnings(which.max(dt$y_tmp))[1]
    if (is.finite(idx0) && !is.na(idx0) && idx0 >= 1L) lead_ids <- dt$snp_tmp[idx0]
  }
  lead_ids <- unique(as.character(lead_ids))
  lead_pos_dt <- dt[snp_tmp %in% lead_ids, .(CHR_tmp, POS_tmp, row_idx = .I)]
  if (nrow(lead_pos_dt)) {
    .get_pilot_add_anchor(lead_pos_dt$row_idx)
    win <- lead_pos_dt[, .(CHR_tmp, start = pmax(1, POS_tmp - lead.flank), end = POS_tmp + lead.flank)]
    q <- dt[, .(i = .I, CHR_tmp, POS_tmp)]
    hit <- win[q, on = .(CHR_tmp, start <= POS_tmp, end >= POS_tmp), nomatch = 0L, allow.cartesian = TRUE]
    if (nrow(hit)) .get_pilot_keep_i(unique(hit$i))
  }

  .gcanvas_note("gcanvas::get.pilot", "Including significant regions and anchor variants", silent = silent)
  sig_idx <- which(dt$y_tmp >= thr_y)
  if (length(sig_idx)) {
    .get_pilot_keep_i(sig_idx)
    sig_win <- dt[sig_idx, .(CHR_tmp, start = pmax(1, POS_tmp - threshold.flank), end = POS_tmp + threshold.flank)]
    q <- dt[, .(i = .I, CHR_tmp, POS_tmp)]
    hit <- sig_win[q, on = .(CHR_tmp, start <= POS_tmp, end >= POS_tmp), nomatch = 0L, allow.cartesian = TRUE]
    if (nrow(hit)) .get_pilot_keep_i(unique(hit$i))
  }

  .gcanvas_note("gcanvas::get.pilot", "Sampling sparse background points", silent = silent)
  n.variant <- as_int(n.variant); if (!is.finite(n.variant) || n.variant < 1L) n.variant <- 30000L
  per.chrom.min <- as_int(per.chrom.min); if (!is.finite(per.chrom.min) || per.chrom.min < 0L) per.chrom.min <- 50L
  selected_idx0 <- which(keep)
  if (length(selected_idx0) > n.variant) {
    .gcanvas_note("gcanvas::get.pilot", "Downsampling prioritized set with chromosome spread", silent = silent)
    all_chr <- unique(dt$CHR_tmp)
    sel_chr <- unique(dt$CHR_tmp[selected_idx0])
    miss_chr <- setdiff(all_chr, sel_chr)
    bg_pool <- which(!keep)

    reserve_bg <- 0L
    if (length(bg_pool)) {
      reserve_bg <- max(as_int(floor(n.variant * 0.2)), length(miss_chr))
      reserve_bg <- min(reserve_bg, length(bg_pool))
      reserve_bg <- max(0L, reserve_bg)
    }
    main_budget <- max(0L, n.variant - reserve_bg)

    forced <- unique(as_int(anchor_idx[anchor_idx %in% selected_idx0]))
    if (length(forced) > main_budget) {
      forced <- .get_pilot_sample_spread(forced, main_budget, min_per_chr = 1L)
    }
    main_left <- max(0L, main_budget - length(forced))
    main_pool <- setdiff(selected_idx0, forced)
    picked_main <- .get_pilot_sample_spread(main_pool, main_left, min_per_chr = if (main_left >= length(all_chr)) 1L else 0L)
    final_idx <- unique(as_int(c(forced, picked_main)))

    if (reserve_bg > 0L) {
      picked_bg <- .get_pilot_sample_spread(bg_pool, reserve_bg, min_per_chr = if (reserve_bg >= length(all_chr)) 1L else 0L, prefer_chr = miss_chr)
      final_idx <- unique(as_int(c(final_idx, picked_bg)))
    }
    if (length(final_idx) < n.variant) {
      extra_pool <- setdiff(seq_len(nrow(dt)), final_idx)
      extra_take <- .get_pilot_sample_spread(extra_pool, n.variant - length(final_idx), min_per_chr = 0L)
      final_idx <- unique(as_int(c(final_idx, extra_take)))
    }
    keep <- rep(FALSE, nrow(dt))
    keep[final_idx] <- TRUE
  }

  need <- max(0L, n.variant - sum(keep, na.rm = TRUE))
  rem_idx <- which(!keep)

  if (need > 0L && length(rem_idx)) {
    set.seed(as_int(seed %||% 1L))
    rem <- dt[rem_idx, .(idx = rem_idx, CHR_tmp, POS_tmp, chr_order)]
    if (identical(method, "uniform")) {
      data.table::setorderv(rem, c("chr_order", "POS_tmp", "idx"), c(1L, 1L, 1L), na.last = TRUE)
      pick <- unique(as_int(round(seq(1, nrow(rem), length.out = min(need, nrow(rem))))))
      .get_pilot_keep_i(rem$idx[pick])
    } else {
      .get_pilot_keep_i(.get_pilot_sample_spread(rem$idx, need, min_per_chr = if (need >= length(unique(rem$CHR_tmp))) 1L else 0L))
    }
  }

  out_idx <- dt$row_id[keep]
  out_idx <- out_idx[is.finite(out_idx) & !is.na(out_idx)]
  out_idx <- unique(as_int(out_idx))
  out <- if (data.table::is.data.table(data)) data[out_idx] else data[out_idx, , drop = FALSE]

  attr(out, "gcanvas_meta") <- list(
    type = "get.pilot",
    n_input = as_int(nrow(data)),
    n_output = as_int(nrow(out)),
    sampling = list(
      n_variant = as_int(n.variant),
      threshold = as.numeric(threshold),
      lead_flank = as.numeric(lead.flank),
      threshold_flank = as.numeric(threshold.flank),
      method = method
    ),
    retained = list(
      n_lead = as_int(length(lead_ids)),
      n_significant = as_int(length(sig_idx)),
      n_priority = as_int(sum(keep, na.rm = TRUE))
    )
  )

  if (identical(return.type, "data.frame") && data.table::is.data.table(out)) return(as.data.frame(out))
  if (identical(return.type, "data.table") && !data.table::is.data.table(out)) return(data.table::as.data.table(out))
  out
}

