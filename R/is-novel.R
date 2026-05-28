# Novelty annotation: which hits in `x` are >= `distance` away from every
# entry in `reported` (and, optionally, also distant from a `novel` list).

#' Mark genomic hits as novel relative to a reference set
#'
#' For each entry of `x`, returns whether it sits at least `distance` bp from
#' every entry of `reported` on the same chromosome. Optionally compares to a
#' second `novel` set as well.
#'
#' @param x A `data.frame`/`data.table` of hits with chromosome and position columns.
#' @param reported A reference set of previously reported variants.
#' @param novel Optional second reference set.
#' @param distance Numeric (bp) or string (e.g. `"1Mb"`) defining the
#'   proximity window that counts as "not novel".
#' @param ... Further options (column overrides, etc.).
#'
#' @return A logical or annotated vector / table parallel to `x`.
#' @export
is.novel <- function(x, reported = NULL, novel = NULL, distance = "1Mb", ...) {
  require_pkg("data.table")

  dots <- list(...)
  dot_nms <- names(dots)
  if (length(dots) && (is.null(dot_nms) || any(is.na(dot_nms) | !nzchar(dot_nms)))) {
    stop("All ... arguments must be named.", call. = FALSE)
  }
  allowed <- c(
    "chrom", "pos", "start", "end",
    "chrom.col", "pos.col", "start.col", "end.col",
    "reported.chrom", "reported.pos", "reported.start", "reported.end",
    "reported.chrom.col", "reported.pos.col", "reported.start.col", "reported.end.col",
    "novel.chrom", "novel.pos", "novel.start", "novel.end",
    "novel.chrom.col", "novel.pos.col", "novel.start.col", "novel.end.col"
  )
  bad <- setdiff(dot_nms %||% character(), allowed)
  if (length(bad)) stop("Unknown ... argument(s): ", paste(bad, collapse = ", "), call. = FALSE)

  dist_bp <- .gcanvas_parse_bp_span(distance, arg_name = "distance")
  if (!is.finite(dist_bp) || is.na(dist_bp) || dist_bp < 0) {
    stop("distance must be non-negative bp span (numeric or string like '1Mb').", call. = FALSE)
  }

  chrom_col <- as.character(dots[["chrom.col"]] %||% "CHR")[1]
  pos_col <- as.character(dots[["pos.col"]] %||% "POS")[1]
  start_col <- as.character(dots[["start.col"]] %||% "START")[1]
  end_col <- as.character(dots[["end.col"]] %||% "END")[1]

  rep_chrom_col <- as.character(dots[["reported.chrom.col"]] %||% chrom_col)[1]
  rep_pos_col <- as.character(dots[["reported.pos.col"]] %||% pos_col)[1]
  rep_start_col <- as.character(dots[["reported.start.col"]] %||% start_col)[1]
  rep_end_col <- as.character(dots[["reported.end.col"]] %||% end_col)[1]
  novel_chrom_col <- as.character(dots[["novel.chrom.col"]] %||% chrom_col)[1]
  novel_pos_col <- as.character(dots[["novel.pos.col"]] %||% pos_col)[1]
  novel_start_col <- as.character(dots[["novel.start.col"]] %||% start_col)[1]
  novel_end_col <- as.character(dots[["novel.end.col"]] %||% end_col)[1]

  .is_novel_pick_ci <- function(nms, cands) {
    nms <- as.character(nms %||% character())
    cands <- as.character(cands %||% character())
    cands <- cands[!is.na(cands) & nzchar(cands)]
    if (!length(nms) || !length(cands)) return(NA_character_)
    nms_low <- tolower(nms)
    for (cc in cands) {
      i <- match(tolower(cc), nms_low, nomatch = 0L)
      if (i > 0L) return(nms[i])
    }
    NA_character_
  }

  .is_novel_to_tabular_list <- function(obj) {
    if (!is.list(obj) || is.data.frame(obj) || data.table::is.data.table(obj)) return(NULL)
    nm <- names(obj)
    if (is.null(nm) || !length(nm)) return(NULL)
    if (any(is.na(nm) | !nzchar(nm))) return(NULL)
    ln <- vapply(obj, length, integer(1))
    if (!length(ln) || any(ln == 0L) || length(unique(ln)) != 1L) return(NULL)
    data.table::as.data.table(obj)
  }

  .is_novel_parse_chr_pos_vectors <- function(chrom, pos, who = "x") {
    if (is.null(chrom) || is.null(pos)) return(NULL)
    chr <- as.character(chrom)
    p <- suppressWarnings(as.numeric(pos))
    n <- length(chr)
    if (length(p) != n) {
      stop(sprintf("%s: chrom and pos must have the same length.", who), call. = FALSE)
    }
    data.table::data.table(
      idx = seq_len(n),
      CHR = normalize.chrom(chr),
      START = p,
      END = p
    )
  }

  .is_novel_parse_chr_range_vectors <- function(chrom, start, end, who = "x") {
    if (is.null(chrom) || is.null(start) || is.null(end)) return(NULL)
    chr <- as.character(chrom)
    s <- suppressWarnings(as.numeric(start))
    e <- suppressWarnings(as.numeric(end))
    n <- length(chr)
    if (length(s) != n || length(e) != n) {
      stop(sprintf("%s: chrom/start/end must have the same length.", who), call. = FALSE)
    }
    data.table::data.table(
      idx = seq_len(n),
      CHR = normalize.chrom(chr),
      START = s,
      END = e
    )
  }

  .is_novel_parse_char_tokens <- function(v) {
    vv <- as.character(v)
    n <- length(vv)
    out <- data.table::data.table(
      idx = seq_len(n),
      CHR = rep(NA_character_, n),
      START = rep(NA_real_, n),
      END = rep(NA_real_, n)
    )
    if (!n) return(out)
    s <- trimws(vv)
    s[is.na(s) | !nzchar(s)] <- NA_character_
    s2 <- gsub(",", "", s, fixed = TRUE)

    is_locus <- !is.na(s2) & grepl("^[^:]+:[0-9]+-[0-9]+$", s2)
    if (any(is_locus)) {
      m <- regexec("^([^:]+):([0-9]+)-([0-9]+)$", s2[is_locus], perl = TRUE)
      mm <- regmatches(s2[is_locus], m)
      out$CHR[is_locus] <- normalize.chrom(vapply(mm, `[`, character(1), 2))
      out$START[is_locus] <- suppressWarnings(as.numeric(vapply(mm, `[`, character(1), 3)))
      out$END[is_locus] <- suppressWarnings(as.numeric(vapply(mm, `[`, character(1), 4)))
    }

    left <- !is_locus & !is.na(s2)
    if (any(left)) {
      sp <- strsplit(s2[left], ":", fixed = TRUE)
      chr <- vapply(sp, function(z) if (length(z) >= 1L) z[1] else NA_character_, character(1))
      pos <- suppressWarnings(as.numeric(vapply(sp, function(z) if (length(z) >= 2L) z[2] else NA_character_, character(1))))
      out$CHR[left] <- normalize.chrom(chr)
      out$START[left] <- pos
      out$END[left] <- pos
    }
    out
  }

  .is_novel_parse_table <- function(dt, chrom_col, pos_col, start_col, end_col) {
    nms <- names(dt)
    chr_use <- .is_novel_pick_ci(nms, c(chrom_col, "CHR", "chrom", "chr"))
    pos_use <- .is_novel_pick_ci(nms, c(pos_col, "POS", "pos", "BP", "bp"))
    start_use <- .is_novel_pick_ci(nms, c(start_col, "START", "start"))
    end_use <- .is_novel_pick_ci(nms, c(end_col, "END", "end"))

    if (!is.na(chr_use) && !is.na(start_use) && !is.na(end_use)) {
      return(data.table::data.table(
        idx = seq_len(nrow(dt)),
        CHR = normalize.chrom(dt[[chr_use]]),
        START = suppressWarnings(as.numeric(dt[[start_use]])),
        END = suppressWarnings(as.numeric(dt[[end_use]]))
      ))
    }
    if (!is.na(chr_use) && !is.na(pos_use)) {
      pp <- suppressWarnings(as.numeric(dt[[pos_use]]))
      return(data.table::data.table(
        idx = seq_len(nrow(dt)),
        CHR = normalize.chrom(dt[[chr_use]]),
        START = pp,
        END = pp
      ))
    }

    if (ncol(dt) >= 3L) {
      return(data.table::data.table(
        idx = seq_len(nrow(dt)),
        CHR = normalize.chrom(dt[[1]]),
        START = suppressWarnings(as.numeric(dt[[2]])),
        END = suppressWarnings(as.numeric(dt[[3]]))
      ))
    }
    if (ncol(dt) == 2L) {
      pp <- suppressWarnings(as.numeric(dt[[2]]))
      return(data.table::data.table(
        idx = seq_len(nrow(dt)),
        CHR = normalize.chrom(dt[[1]]),
        START = pp,
        END = pp
      ))
    }
    stop("Could not parse tabular input: need CHR+POS or CHR+START+END (or 2/3 columns).", call. = FALSE)
  }

  .is_novel_parse_any <- function(obj,
                         chrom_vec = NULL, pos_vec = NULL, start_vec = NULL, end_vec = NULL,
                         chrom_col = "CHR", pos_col = "POS", start_col = "START", end_col = "END",
                         who = "x") {
    if (!is.null(pos_vec) && (!is.null(start_vec) || !is.null(end_vec))) {
      stop(sprintf("%s: provide either chrom+pos or chrom+start+end, not both.", who), call. = FALSE)
    }
    dt_vec <- .is_novel_parse_chr_pos_vectors(chrom_vec, pos_vec, who = who)
    if (is.null(dt_vec)) dt_vec <- .is_novel_parse_chr_range_vectors(chrom_vec, start_vec, end_vec, who = who)
    if (!is.null(dt_vec)) {
      out <- dt_vec
      nm0 <- names(chrom_vec) %||% names(pos_vec) %||% names(start_vec)
      if (is.null(nm0) || !length(nm0) || !any(!is.na(nm0) & nzchar(nm0))) nm0 <- NULL
      attr(out, "orig_names") <- nm0
      return(out)
    }

    tab_obj <- .is_novel_to_tabular_list(obj)
    if (!is.null(tab_obj)) obj <- tab_obj

    if (is.data.frame(obj) || data.table::is.data.table(obj)) {
      dt0 <- if (data.table::is.data.table(obj)) data.table::copy(obj) else data.table::as.data.table(obj)
      out <- .is_novel_parse_table(dt0, chrom_col = chrom_col, pos_col = pos_col, start_col = start_col, end_col = end_col)
      attr(out, "orig_names") <- NULL
      return(out)
    }

    vals <- .gcanvas_as_snp_vector(obj)
    out <- .is_novel_parse_char_tokens(vals)
    nm0 <- names(vals)
    if (is.null(nm0) || !length(nm0) || !any(!is.na(nm0) & nzchar(nm0))) nm0 <- NULL
    attr(out, "orig_names") <- nm0
    out
  }

  x_dt <- .is_novel_parse_any(
    obj = x,
    chrom_vec = dots[["chrom"]] %||% NULL,
    pos_vec = dots[["pos"]] %||% NULL,
    start_vec = dots[["start"]] %||% NULL,
    end_vec = dots[["end"]] %||% NULL,
    chrom_col = chrom_col, pos_col = pos_col, start_col = start_col, end_col = end_col,
    who = "x"
  )
  novel_mode <- !is.null(novel) || !is.null(dots[["novel.chrom"]]) || !is.null(dots[["novel.pos"]]) || !is.null(dots[["novel.start"]]) || !is.null(dots[["novel.end"]])
  if (!novel_mode && is.null(reported) && is.null(dots[["reported.chrom"]]) && is.null(dots[["reported.pos"]]) &&
      is.null(dots[["reported.start"]]) && is.null(dots[["reported.end"]])) {
    stop("Provide either reported or novel.", call. = FALSE)
  }

  rep_dt <- NULL
  novel_dt <- NULL
  if (isTRUE(novel_mode)) {
    novel_dt <- .is_novel_parse_any(
      obj = novel,
      chrom_vec = dots[["novel.chrom"]] %||% NULL,
      pos_vec = dots[["novel.pos"]] %||% NULL,
      start_vec = dots[["novel.start"]] %||% NULL,
      end_vec = dots[["novel.end"]] %||% NULL,
      chrom_col = novel_chrom_col, pos_col = novel_pos_col, start_col = novel_start_col, end_col = novel_end_col,
      who = "novel"
    )
  } else {
    rep_dt <- .is_novel_parse_any(
      obj = reported,
      chrom_vec = dots[["reported.chrom"]] %||% NULL,
      pos_vec = dots[["reported.pos"]] %||% NULL,
      start_vec = dots[["reported.start"]] %||% NULL,
      end_vec = dots[["reported.end"]] %||% NULL,
      chrom_col = rep_chrom_col, pos_col = rep_pos_col, start_col = rep_start_col, end_col = rep_end_col,
      who = "reported"
    )
  }

  if (!nrow(x_dt)) return(logical())
  if (isTRUE(novel_mode) && !nrow(novel_dt)) {
    out0 <- rep(FALSE, nrow(x_dt))
    nm <- attr(x_dt, "orig_names")
    if (!is.null(nm) && length(nm) == length(out0)) names(out0) <- nm
    return(out0)
  }
  if (!isTRUE(novel_mode) && !nrow(rep_dt)) {
    out0 <- rep(TRUE, nrow(x_dt))
    nm <- attr(x_dt, "orig_names")
    if (!is.null(nm) && length(nm) == length(out0)) names(out0) <- nm
    return(out0)
  }

  x_dt[, `:=`(
    CHR = normalize.chrom(CHR),
    START = suppressWarnings(as.numeric(START)),
    END = suppressWarnings(as.numeric(END))
  )]
  if (isTRUE(novel_mode)) {
    novel_dt[, `:=`(
      CHR = normalize.chrom(CHR),
      START = suppressWarnings(as.numeric(START)),
      END = suppressWarnings(as.numeric(END))
    )]
  } else {
    rep_dt[, `:=`(
      CHR = normalize.chrom(CHR),
      START = suppressWarnings(as.numeric(START)),
      END = suppressWarnings(as.numeric(END))
    )]
  }

  swap_x <- is.finite(x_dt$START) & is.finite(x_dt$END) & x_dt$START > x_dt$END
  if (any(swap_x)) {
    tmp <- x_dt$START[swap_x]
    x_dt$START[swap_x] <- x_dt$END[swap_x]
    x_dt$END[swap_x] <- tmp
  }
  if (isTRUE(novel_mode)) {
    swap_n <- is.finite(novel_dt$START) & is.finite(novel_dt$END) & novel_dt$START > novel_dt$END
    if (any(swap_n)) {
      tmp <- novel_dt$START[swap_n]
      novel_dt$START[swap_n] <- novel_dt$END[swap_n]
      novel_dt$END[swap_n] <- tmp
    }
  } else {
    swap_r <- is.finite(rep_dt$START) & is.finite(rep_dt$END) & rep_dt$START > rep_dt$END
    if (any(swap_r)) {
      tmp <- rep_dt$START[swap_r]
      rep_dt$START[swap_r] <- rep_dt$END[swap_r]
      rep_dt$END[swap_r] <- tmp
    }
  }

  x_valid <- !is.na(x_dt$CHR) & nzchar(x_dt$CHR) &
    is.finite(x_dt$START) & is.finite(x_dt$END) &
    !is.na(x_dt$START) & !is.na(x_dt$END)
  if (isTRUE(novel_mode)) {
    novel_dt <- novel_dt[!is.na(CHR) & nzchar(CHR) & is.finite(START) & is.finite(END) & !is.na(START) & !is.na(END)]
    novel_dt <- unique(novel_dt, by = c("CHR", "START", "END"))
  } else {
    rep_dt <- rep_dt[!is.na(CHR) & nzchar(CHR) & is.finite(START) & is.finite(END) & !is.na(START) & !is.na(END)]
    rep_dt <- unique(rep_dt, by = c("CHR", "START", "END"))
  }

  out <- rep(NA, nrow(x_dt))
  if (!any(x_valid)) {
    nm <- attr(x_dt, "orig_names")
    if (!is.null(nm) && length(nm) == length(out)) names(out) <- nm
    return(out)
  }
  if (isTRUE(novel_mode) && !nrow(novel_dt)) {
    out[x_valid] <- FALSE
    nm <- attr(x_dt, "orig_names")
    if (!is.null(nm) && length(nm) == length(out)) names(out) <- nm
    return(out)
  }
  if (!isTRUE(novel_mode) && !nrow(rep_dt)) {
    out[x_valid] <- TRUE
    nm <- attr(x_dt, "orig_names")
    if (!is.null(nm) && length(nm) == length(out)) names(out) <- nm
    return(out)
  }

  x_work <- x_dt[x_valid, .(idx = idx, CHR, START, END)]
  x_work[, chr_order := rank.chrom(CHR)]
  data.table::setorderv(x_work, c("chr_order", "CHR", "START", "END", "idx"), c(1L, 1L, 1L, 1L, 1L), na.last = TRUE)
  x_work[, chr_order := NULL]
  chr_levels <- unique(as.character(x_work$CHR))

  for (chr_i in chr_levels) {
    x_chr <- x_work[CHR == chr_i]
    r_chr <- if (isTRUE(novel_mode)) novel_dt[CHR == chr_i, .(START, END)] else rep_dt[CHR == chr_i, .(START, END)]
    if (!nrow(x_chr)) next
    if (!nrow(r_chr)) {
      out[x_chr$idx] <- if (isTRUE(novel_mode)) FALSE else TRUE
      next
    }

    if (isTRUE(novel_mode)) {
      x_ov <- data.table::copy(x_chr[, .(qidx = seq_len(.N), START, END)])
      r_ov <- data.table::copy(r_chr[, .(START, END)])
      data.table::setkey(x_ov, START, END)
      data.table::setkey(r_ov, START, END)
      ov <- data.table::foverlaps(x_ov, r_ov, by.x = c("START", "END"), by.y = c("START", "END"), nomatch = 0L, type = "any")
      hit <- rep(FALSE, nrow(x_chr))
      if (nrow(ov)) hit[unique(ov$qidx)] <- TRUE
      out[x_chr$idx] <- hit
      next
    }

    d <- rep(Inf, nrow(x_chr))

    x_ov <- data.table::copy(x_chr[, .(qidx = seq_len(.N), START, END)])
    r_ov <- data.table::copy(r_chr[, .(START, END)])
    data.table::setkey(x_ov, START, END)
    data.table::setkey(r_ov, START, END)
    ov <- data.table::foverlaps(x_ov, r_ov, by.x = c("START", "END"), by.y = c("START", "END"), nomatch = 0L, type = "any")
    if (nrow(ov)) d[unique(ov$qidx)] <- 0

    r_end <- sort(as.numeric(r_chr$END))
    r_start <- sort(as.numeric(r_chr$START))
    xs <- as.numeric(x_chr$START)
    xe <- as.numeric(x_chr$END)

    i_left <- findInterval(xs, r_end)
    d_left <- rep(Inf, length(xs))
    ok_left <- i_left > 0L
    d_left[ok_left] <- xs[ok_left] - r_end[i_left[ok_left]]

    i_right <- findInterval(xe, r_start)
    d_right <- rep(Inf, length(xe))
    ok_right <- i_right < length(r_start)
    d_right[ok_right] <- r_start[i_right[ok_right] + 1L] - xe[ok_right]

    d <- pmin(d, d_left, d_right)
    d[d < 0] <- 0
    out[x_chr$idx] <- (d >= dist_bp)
  }

  nm <- attr(x_dt, "orig_names")
  if (!is.null(nm) && length(nm) == length(out)) names(out) <- nm
  out
}

