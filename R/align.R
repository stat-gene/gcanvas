# Allele-aware alignment helpers for effect estimates, frequencies, and LD r.
# `ld.align()` lives further down (defined in this file alongside its
# matching beta/maf counterparts so all the alignment surface lives together).

#' Align effect sizes (betas) to a reference allele
#'
#' Flips the sign of `x` where the dataset's A2 matches the reference A1, sets
#' entries with no matching allele to `NA`, and returns the aligned vector.
#'
#' @param x Numeric vector of effect sizes.
#' @param a1.x Effect allele for `x` (the allele the beta is reported against).
#' @param a2.x Non-effect (other) allele for `x`.
#' @param a1.ref Reference effect allele to align against.
#'
#' @return Numeric vector of `x` aligned to `a1.ref`. Unresolvable entries
#'   become `NA_real_`.
#' @export
beta.align <- function(x, a1.x, a2.x, a1.ref) {
  n <- length(x)
  if (length(a1.x) != n || length(a2.x) != n || length(a1.ref) != n) {
    stop("x, a1.x, a2.x, and a1.ref must have the same length.", call. = FALSE)
  }
  out <- as_num(x)
  ax1 <- toupper(as.character(a1.x))
  ax2 <- toupper(as.character(a2.x))
  ar <- toupper(as.character(a1.ref))

  no_ref <- is.na(ar) | !nzchar(ar)
  hit_a1 <- !no_ref & !is.na(ax1) & nzchar(ax1) & (ar == ax1)
  hit_a2 <- !no_ref & !is.na(ax2) & nzchar(ax2) & (ar == ax2)
  unknown <- !no_ref & !(hit_a1 | hit_a2)

  out[hit_a2] <- -1 * out[hit_a2]
  out[no_ref | unknown] <- NA_real_
  out
}

#' Align minor allele frequencies to a reference allele
#'
#' Returns `1 - x` where the dataset's A2 matches the reference A1, sets
#' entries with no matching allele to `NA`, and preserves the rest.
#'
#' @inheritParams beta.align
#' @return Numeric vector of frequencies aligned to `a1.ref`. Unresolvable
#'   entries become `NA_real_`.
#' @export
maf.align <- function(x, a1.x, a2.x, a1.ref) {
  n <- length(x)
  if (length(a1.x) != n || length(a2.x) != n || length(a1.ref) != n) {
    stop("x, a1.x, a2.x, and a1.ref must have the same length.", call. = FALSE)
  }
  out <- as_num(x)
  ax1 <- toupper(as.character(a1.x))
  ax2 <- toupper(as.character(a2.x))
  ar <- toupper(as.character(a1.ref))

  no_ref <- is.na(ar) | !nzchar(ar)
  hit_a1 <- !no_ref & !is.na(ax1) & nzchar(ax1) & (ar == ax1)
  hit_a2 <- !no_ref & !is.na(ax2) & nzchar(ax2) & (ar == ax2)
  unknown <- !no_ref & !(hit_a1 | hit_a2)

  out[hit_a2] <- 1 - out[hit_a2]
  out[no_ref | unknown] <- NA_real_
  out
}

# Backward-compatible alias (kept internal; not exported).
freq.align <- maf.align

#' Align LD r values to a reference allele
#'
#' Sign-flips entries of a long-form LD r pair table or a square r matrix so
#' that values are consistent with a reference effect allele. Expects signed
#' r (not r-squared).
#'
#' @param x A pair `data.frame`/`data.table` with `SNP`, `SNP2`, `r` columns,
#'   a list with a `data` or `R` slot from [calcld()], or a numeric matrix
#'   with `dimnames`.
#' @param reference A `data.frame`/`data.table` with `SNP` and `A1` columns
#'   describing the desired reference allele per variant.
#' @param snp.col Name of the SNP id column in `reference` / `x`.
#' @param a1.col Name of the A1 column in `reference` / `x`.
#' @param a1.x Optional A1 mapping for `x` (named character vector, list, or
#'   table). Required when `x` is a matrix or a bare pair table.
#'
#' @return An object matching the input shape with r values aligned, plus
#'   bookkeeping about which variants were retained.
#' @export
ld.align <- function(x, reference, snp.col = "SNP", a1.col = "A1", a1.x = NULL) {
  require_pkg("data.table")
  .gcanvas_note("gcanvas::ld.align", "IMPORTANT: ld.align expects signed LD r values (not r-squared, r2).", silent = FALSE)

  if (!is.data.frame(reference) && !data.table::is.data.table(reference)) {
    stop("reference must be a data.frame/data.table.", call. = FALSE)
  }
  ref_dt <- if (data.table::is.data.table(reference)) data.table::copy(reference) else data.table::as.data.table(reference)
  ref_snp_col <- .gcanvas_resolve_colname(names(ref_dt), snp.col, aliases = character(), required = TRUE, arg_label = "reference snp.col")
  ref_a1_col <- .gcanvas_resolve_colname(names(ref_dt), a1.col, aliases = c("a1", "A1"), required = TRUE, arg_label = "reference a1.col")
  ref_map_dt <- ref_dt[, .(
    snp = as.character(get(ref_snp_col)),
    a1_ref = toupper(as.character(get(ref_a1_col)))
  )]
  ref_map_dt <- ref_map_dt[!is.na(snp) & nzchar(snp) & !is.na(a1_ref) & nzchar(a1_ref)]
  ref_map_dt <- ref_map_dt[!duplicated(snp)]
  ref_map <- stats::setNames(ref_map_dt$a1_ref, ref_map_dt$snp)

  .ld_align_as_a1_map <- function(obj) {
    if (is.null(obj) || length(obj) == 0L) return(NULL)
    if (is.data.frame(obj) || data.table::is.data.table(obj)) {
      dt0 <- if (data.table::is.data.table(obj)) data.table::copy(obj) else data.table::as.data.table(obj)
      snp_use <- .gcanvas_resolve_colname(names(dt0), snp.col, aliases = character(), required = TRUE, arg_label = "a1.x snp.col")
      a1_use <- .gcanvas_resolve_colname(names(dt0), a1.col, aliases = character(), required = TRUE, arg_label = "a1.x a1.col")
      map_dt <- dt0[, .(snp = as.character(get(snp_use)), a1 = toupper(as.character(get(a1_use))))]
      map_dt <- map_dt[!is.na(snp) & nzchar(snp)]
      map_dt <- map_dt[!duplicated(snp)]
      return(stats::setNames(map_dt$a1, map_dt$snp))
    }
    if (is.list(obj) && !is.data.frame(obj) && !data.table::is.data.table(obj)) {
      obj <- unlist(obj, use.names = TRUE)
    }
    nms <- names(obj)
    vv <- as.character(obj)
    if (is.null(nms) || !length(nms)) {
      stop("a1.x must be a named vector/list or a data.frame/data.table with snp/a1 columns.", call. = FALSE)
    }
    vv <- toupper(vv)
    vv[is.na(vv) | !nzchar(vv)] <- NA_character_
    nms <- as.character(nms)
    keep <- !is.na(nms) & nzchar(nms)
    vv <- vv[keep]; nms <- nms[keep]
    keep2 <- !duplicated(nms)
    stats::setNames(vv[keep2], nms[keep2])
  }

  .ld_align_pair_dt <- function(dt_in, a1_map_in) {
    dt <- if (data.table::is.data.table(dt_in)) data.table::copy(dt_in) else data.table::as.data.table(dt_in)
    s1_col <- .gcanvas_resolve_colname(names(dt), "SNP", aliases = c("markerID", "snp1", "id1"), required = TRUE, arg_label = "pair SNP")
    s2_col <- .gcanvas_resolve_colname(names(dt), "SNP2", aliases = c("markerID2", "snp2", "id2"), required = TRUE, arg_label = "pair SNP2")
    r_col <- .gcanvas_resolve_colname(names(dt), "r", aliases = c("R", "ld_r", "cor", "corr"), required = FALSE)
    if (is.null(r_col)) {
      if (!is.null(.gcanvas_resolve_colname(names(dt), "r2", aliases = c("R2", "ld_r2"), required = FALSE))) {
        stop("ld.align expects LD r (signed), but detected only r2 column.", call. = FALSE)
      }
      stop("Could not find LD r column in x. Expected one of: r, R, ld_r, cor, corr.", call. = FALSE)
    }
    if (tolower(r_col) == "r2") stop("ld.align expects LD r (signed), not r2.", call. = FALSE)

    dt[, `:=`(
      .snp1 = as.character(get(s1_col)),
      .snp2 = as.character(get(s2_col)),
      .r = suppressWarnings(as.numeric(get(r_col)))
    )]
    vars_in <- unique(c(dt$.snp1, dt$.snp2))
    vars_in <- vars_in[!is.na(vars_in) & nzchar(vars_in)]
    keep_vars <- intersect(vars_in, names(ref_map))
    dt <- dt[.snp1 %in% keep_vars & .snp2 %in% keep_vars]
    if (!nrow(dt)) {
      dt[, c(".snp1", ".snp2", ".r") := NULL]
      return(list(data = dt[0], a1_ref = stats::setNames(character(), character())))
    }

    a1_map <- .ld_align_as_a1_map(a1_map_in)
    if (is.null(a1_map)) {
      stop("a1.x is required unless x contains named a1 information (e.g., calcld list$a1).", call. = FALSE)
    }
    common <- intersect(keep_vars, names(a1_map))
    dt <- dt[.snp1 %in% common & .snp2 %in% common]
    if (!nrow(dt)) {
      dt[, c(".snp1", ".snp2", ".r") := NULL]
      return(list(data = dt[0], a1_ref = stats::setNames(character(), character())))
    }
    a1a <- toupper(as.character(a1_map[dt$.snp1]))
    a1b <- toupper(as.character(a1_map[dt$.snp2]))
    refa <- toupper(as.character(ref_map[dt$.snp1]))
    refb <- toupper(as.character(ref_map[dt$.snp2]))

    d1 <- a1a != refa
    d2 <- a1b != refb
    invalid <- is.na(d1) | is.na(d2)
    flip <- !invalid & xor(d1, d2)
    r_new <- dt$.r
    r_new[flip] <- -1 * r_new[flip]
    r_new[invalid] <- NA_real_
    dt[, (r_col) := r_new]
    dt[, c(".snp1", ".snp2", ".r") := NULL]
    used <- unique(c(dt[[s1_col]], dt[[s2_col]]))
    used <- as.character(used[!is.na(used) & nzchar(as.character(used))])
    list(data = dt[], a1_ref = ref_map[used])
  }

  .ld_align_matrix <- function(M, a1_map_in) {
    M <- suppressWarnings(as.matrix(M))
    if (!is.matrix(M) || !is.numeric(M)) stop("LD matrix must be a numeric matrix.", call. = FALSE)
    rn <- rownames(M); cn <- colnames(M)
    if (is.null(rn) || is.null(cn)) stop("LD matrix must have rownames and colnames.", call. = FALSE)
    vars <- intersect(rn, cn)
    vars <- vars[!is.na(vars) & nzchar(vars)]
    if (!length(vars)) return(list(R = M[0, 0, drop = FALSE], a1_ref = stats::setNames(character(), character())))

    a1_map <- .ld_align_as_a1_map(a1_map_in)
    if (is.null(a1_map)) {
      stop("For matrix input, a1.x must be provided (or list$a1 must exist for list input).", call. = FALSE)
    }
    vars <- intersect(vars, names(ref_map))
    vars <- intersect(vars, names(a1_map))
    if (!length(vars)) return(list(R = M[0, 0, drop = FALSE], a1_ref = stats::setNames(character(), character())))
    M2 <- M[vars, vars, drop = FALSE]
    d <- toupper(as.character(a1_map[vars])) != toupper(as.character(ref_map[vars]))
    ok <- !is.na(d)
    vars <- vars[ok]
    d <- d[ok]
    if (!length(vars)) return(list(R = M2[0, 0, drop = FALSE], a1_ref = stats::setNames(character(), character())))
    M2 <- M2[vars, vars, drop = FALSE]
    s <- ifelse(d, -1, 1)
    M2 <- outer(s, s, `*`) * M2
    list(R = M2, a1_ref = ref_map[vars])
  }

  # Input type 1: calcld-style list
  if (is.list(x) && !is.data.frame(x) && !data.table::is.data.table(x)) {
    if ("data" %in% names(x) && (is.data.frame(x$data) || data.table::is.data.table(x$data))) {
      a1_src <- if (!is.null(x$a1) && length(x$a1)) x$a1 else a1.x
      res <- .ld_align_pair_dt(x$data, a1_src)
      out <- x
      out$data <- res$data
      out$a1 <- res$a1_ref
      if (is.list(out$meta)) {
        out$meta$aligned <- TRUE
        out$meta$n_variant <- as_int(length(res$a1_ref))
      }
      return(out)
    }
    if ("R" %in% names(x) && is.matrix(x$R)) {
      a1_src <- if (!is.null(x$a1) && length(x$a1)) x$a1 else a1.x
      res <- .ld_align_matrix(x$R, a1_src)
      out <- x
      out$R <- res$R
      out$a1 <- res$a1_ref
      if (is.list(out$meta)) {
        out$meta$aligned <- TRUE
        out$meta$n_variant <- as_int(length(res$a1_ref))
      }
      return(out)
    }
  }

  # Input type 3: pair data.table/data.frame
  if (is.data.frame(x) || data.table::is.data.table(x)) {
    return(.ld_align_pair_dt(x, a1.x)$data)
  }

  # Input type 2: single named matrix
  if (is.matrix(x)) {
    res <- .ld_align_matrix(x, a1.x)
    return(res$R)
  }

  stop("Unsupported x type. Use calcld result/list, pair table, or named LD matrix.", call. = FALSE)
}

