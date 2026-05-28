# Round-trip conversion between long-form pair tables and symmetric matrices.

#' Pair table to symmetric matrix
#'
#' Pivots a 3-column `(id1, id2, value)` table to a symmetric matrix. Columns
#' are auto-detected (case-insensitive aliases such as `SNP`/`SNP2`/`r2`/`r`)
#' or, failing that, the first two columns are treated as ids and the third
#' as the value.
#'
#' @param x A `data.frame`/`data.table` (or a list with a `data` slot).
#' @param diag Value placed on the matrix diagonal (default `1`).
#' @param na Value used for missing pairs (default `NA`).
#'
#' @return A symmetric numeric matrix with row/column names from the ids.
#' @export
dt2mat <- function(x, diag = 1, na = NA) {
  require_pkg("data.table")
  if (is.list(x) && !is.data.frame(x) && !data.table::is.data.table(x) && ("data" %in% names(x))) {
    x <- x$data
  }
  if (!is.data.frame(x) && !data.table::is.data.table(x)) {
    stop("x must be a data.frame/data.table.", call. = FALSE)
  }
  dt <- if (data.table::is.data.table(x)) data.table::copy(x) else data.table::as.data.table(x)
  if (!nrow(dt)) return(matrix(na, nrow = 0L, ncol = 0L))

  nms <- names(dt)
  nms_low <- tolower(nms)

  pick_col <- function(cands) {
    hit <- which(nms_low %in% cands)
    if (!length(hit)) return(NA_character_)
    nms[hit[1]]
  }

  id1_col <- pick_col(c("snp", "markerid", "snp1", "id1", "variant1", "var1", "rsid1"))
  id2_col <- pick_col(c("snp2", "markerid2", "id2", "variant2", "var2", "rsid2"))

  if (is.na(id1_col) || is.na(id2_col)) {
    if (ncol(dt) == 3L) {
      id1_col <- nms[1]
      id2_col <- nms[2]
    } else if (!is.na(id1_col) && is.na(id2_col) && ncol(dt) >= 2L) {
      id2_col <- setdiff(nms, id1_col)[1]
    } else if (is.na(id1_col) && !is.na(id2_col) && ncol(dt) >= 2L) {
      id1_col <- setdiff(nms, id2_col)[1]
    } else {
      stop("Could not detect pair columns. Provide markerID/markerID2 or 3-column pair table.", call. = FALSE)
    }
  }

  val_col <- pick_col(c("r2", "r", "value", "val", "ld"))
  if (is.na(val_col)) {
    if (ncol(dt) == 3L) {
      val_col <- setdiff(nms, c(id1_col, id2_col))[1]
    } else {
      cand <- setdiff(nms, c(id1_col, id2_col))
      if (!length(cand)) stop("Could not detect value column.", call. = FALSE)
      val_col <- cand[1]
    }
  }

  dt2 <- dt[, .(
    snp1 = as.character(get(id1_col)),
    snp2 = as.character(get(id2_col)),
    val = get(val_col)
  )]
  dt2 <- dt2[!is.na(snp1) & nzchar(snp1) & !is.na(snp2) & nzchar(snp2)]
  if (!nrow(dt2)) return(matrix(na, nrow = 0L, ncol = 0L))

  ids <- unique(c(dt2$snp1, dt2$snp2))
  ids <- ids[!is.na(ids) & nzchar(ids)]
  ids <- as.character(ids)
  if (!length(ids)) return(matrix(na, nrow = 0L, ncol = 0L))

  mat <- matrix(na, nrow = length(ids), ncol = length(ids), dimnames = list(ids, ids))
  ii <- match(dt2$snp1, ids)
  jj <- match(dt2$snp2, ids)
  ok <- !is.na(ii) & !is.na(jj)
  if (any(ok)) {
    mat[cbind(ii[ok], jj[ok])] <- dt2$val[ok]
    mat[cbind(jj[ok], ii[ok])] <- dt2$val[ok]
  }
  diag(mat) <- diag
  mat
}

#' Matrix to pair table
#'
#' Melts a (typically symmetric) matrix into a long-form `data.table` with one
#' row per `(id1, id2)` pair.
#'
#' @param x A matrix, or a list with an `R` slot containing a matrix.
#' @param columns Character vector of length 3 giving the output column names
#'   for `(id1, id2, value)`.
#'
#' @return A `data.table` with the named columns and one row per matrix cell.
#' @export
mat2dt <- function(x, columns = c("snp1", "snp2", "value")) {
  require_pkg("data.table")

  if (is.list(x) && !is.data.frame(x) && !data.table::is.data.table(x) && ("R" %in% names(x))) {
    x <- x$R
  }
  if (!is.matrix(x)) {
    stop("x must be a matrix (or a list containing matrix in x$R).", call. = FALSE)
  }

  col.names <- as.character(columns)
  if (length(col.names) != 3L || any(is.na(col.names) | !nzchar(col.names))) {
    stop("columns must be a character vector of length 3.", call. = FALSE)
  }
  if (length(unique(col.names)) != 3L) {
    stop("columns must contain 3 unique names.", call. = FALSE)
  }

  M <- as.matrix(x)
  nr <- nrow(M)
  nc <- ncol(M)
  if (nr == 0L || nc == 0L) {
    out0 <- data.table::data.table(snp1 = character(), snp2 = character(), value = numeric())
    data.table::setnames(out0, c("snp1", "snp2", "value"), col.names)
    return(out0)
  }

  rn <- rownames(M)
  cn <- colnames(M)
  if (is.null(rn)) rn <- as.character(seq_len(nr))
  if (is.null(cn)) cn <- as.character(seq_len(nc))

  ii <- rep.int(seq_len(nr), times = nc)
  jj <- rep(seq_len(nc), each = nr)
  out <- data.table::data.table(
    snp1 = as.character(rn[ii]),
    snp2 = as.character(cn[jj]),
    value = as.vector(M)
  )
  data.table::setnames(out, c("snp1", "snp2", "value"), col.names)
  out[]
}

