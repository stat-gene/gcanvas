# Chromosome name normalization and ordering, plus a handful of small
# internal helpers (Ensembl id, column-resolution, tabix shell-out) that were
# grouped here in the original script for convenience.

#' Normalize chromosome names to canonical form
#'
#' Strips `chr`/`Chr` prefixes, uppercases the rest, and remaps numeric
#' aliases (`23 -> X`, `24 -> Y`, `25 -> XY`, `26/M/MT -> MT`).
#'
#' @param x Character or numeric vector of chromosome identifiers.
#' @return Character vector of canonical chromosome names.
#' @export
normalize.chrom <- function(x) {
  x <- as.character(x)
  x <- sub("^chr", "", x, ignore.case = TRUE)
  x <- toupper(x)
  x[x == "23"] <- "X"
  x[x == "24"] <- "Y"
  x[x == "25"] <- "XY"
  x[x %in% c("26", "M", "MT")] <- "MT"
  x
}

#' Integer rank for chromosome ordering
#'
#' Maps autosomes to their integer (1-22), then sex chromosomes (X=23, Y=24,
#' XY=25, MT=26), and anything unrecognized to 1000 so it sorts last.
#'
#' @param chr Character or numeric chromosome identifiers.
#' @return Integer vector suitable for `order()` / sorting.
#' @export
rank.chrom <- function(chr) {
  x <- normalize.chrom(chr)
  out <- suppressWarnings(as.integer(x))
  out[is.na(out) & x == "X"] <- 23L
  out[is.na(out) & x == "Y"] <- 24L
  out[is.na(out) & x == "XY"] <- 25L
  out[is.na(out) & x == "MT"] <- 26L
  out[is.na(out)] <- 1000L
  out
}

.gcanvas_sort_chr_unique <- function(x) {
  x <- normalize.chrom(x)
  x <- x[!is.na(x) & nzchar(x)]
  x <- unique(x)
  if (!length(x)) return(character())
  sort.chrom(x)
}

#' Sort by chromosome order
#'
#' Reorders a vector, factor, or data.frame/data.table by natural chromosome
#' order (1, 2, ..., 22, X, Y, XY, MT, then anything else).
#'
#' This is a plain function (the dot in its name does not denote an S3
#' method for [base::sort()]).
#'
#' @param x Vector, factor, or table to sort.
#' @return Object of the same type, reordered.
#' @rawNamespace export("sort.chrom")
sort.chrom <- function(x) {
  x0 <- x
  x_chr <- normalize.chrom(x)
  o <- order(rank.chrom(x_chr), seq_along(x_chr), na.last = TRUE)

  if (is.factor(x0)) {
    return(factor(x0[o], levels = unique(x0[o])))
  }
  if (is.data.frame(x0) || data.table::is.data.table(x0)) {
    return(x0[o, , drop = FALSE])
  }
  x0[o]
}

.gcanvas_hg_chr_bounds <- function(build = 38L) {
  build <- as_int(build)[1]
  if (is.na(build)) build <- 38L
  if (!(build %in% c(37L, 38L))) stop("build must be 37 or 38.", call. = FALSE)
  if (identical(build, 37L)) {
    chr <- c(as.character(1:22), "X", "Y", "MT")
    end <- c(
      249250621, 243199373, 198022430, 191154276, 180915260, 171115067,
      159138663, 146364022, 141213431, 135534747, 135006516, 133851895,
      115169878, 107349540, 102531392, 90354753, 81195210, 78077248,
      59128983, 63025520, 48129895, 51304566, 155270560, 59373566, 16569
    )
  } else {
    chr <- c(as.character(1:22), "X", "Y", "MT")
    end <- c(
      248956422, 242193529, 198295559, 190214555, 181538259, 170805979,
      159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
      114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
      58617616, 64444167, 46709983, 50818468, 156040895, 57227415, 16569
    )
  }
  data.table::data.table(
    CHR = chr,
    start = rep(1, length(chr)),
    end = as.numeric(end)
  )
}

chr_to_ucsc1 <- function(chrom) {
  x <- normalize.chrom(chrom)[1]
  if (is.na(x) || !nzchar(x)) return(NA_character_)
  if (x == "MT") return("chrM")
  paste0("chr", x)
}

.gcanvas_normalize_ensembl_gene_id <- function(x) {
  x <- as.character(x)
  ok <- !is.na(x) & grepl("^ENSG", x, ignore.case = TRUE)
  x[ok] <- sub("\\..*$", "", x[ok])
  x[ok] <- toupper(x[ok])
  x
}

.gcanvas_is_ensembl_gene_id <- function(x) {
  x <- as.character(x)
  !is.na(x) & grepl("^ENSG", x, ignore.case = TRUE)
}

.gcanvas_select_cols_df <- function(df, cols) {
  require_pkg("data.table")
  cols <- as.character(cols)
  miss <- setdiff(cols, names(df))
  if (length(miss)) stop("Missing columns in df: ", paste(miss, collapse = ", "), call. = FALSE)

  if (data.table::is.data.table(df)) {
    out <- df[, ..cols]
    return(as.data.frame(out, stringsAsFactors = FALSE))
  }
  as.data.frame(df[, cols, drop = FALSE], stringsAsFactors = FALSE)
}

.gcanvas_resolve_colname <- function(nms, primary, aliases = character(), required = TRUE, arg_label = "column") {
  nms <- as.character(nms %||% character())
  cand <- unique(c(as.character(primary)[1], as.character(aliases)))
  cand <- cand[!is.na(cand) & nzchar(cand)]
  if (!length(cand)) {
    if (isTRUE(required)) stop(sprintf("Could not resolve %s: no candidate names.", arg_label), call. = FALSE)
    return(NULL)
  }
  hit_exact <- cand[cand %in% nms]
  if (length(hit_exact)) return(hit_exact[1])
  if (isTRUE(required)) {
    stop(sprintf("Could not resolve %s. tried: %s", arg_label, paste(cand, collapse = ", ")), call. = FALSE)
  }
  NULL
}

.gcanvas_tabix_bin <- function() {
  x <- Sys.which("tabix")
  if (!nzchar(x)) stop("tabix binary not found in PATH.", call. = FALSE)
  x
}

.gcanvas_tabix_query_lines <- function(bgzip_path, seqname, start, end) {
  start <- as_int(start); end <- as_int(end)
  if (is.na(start) || is.na(end) || end < start) return(character())

  tabix <- .gcanvas_tabix_bin()
  region <- sprintf("%s:%d-%d", seqname, start, end)

  errf <- tempfile(fileext = ".tabix.err")
  on.exit(unlink(errf), add = TRUE)

  out <- system2(tabix, args = c(bgzip_path, region), stdout = TRUE, stderr = errf)

  out <- out[!grepl("^#", out)]
  out
}

