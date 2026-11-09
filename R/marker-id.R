# Marker IDs (chrom:pos:allele:allele) that do not depend on allele order,
# allele case, position type or the session locale.

#' Marker ID with alphabetically sorted alleles
#'
#' Builds `chrom:pos:a:b` IDs in which `a` and `b` are the two alleles in upper
#' case and in byte (C-locale) order, so a variant gets the same ID whichever
#' allele is listed first and whatever the session locale.
#'
#' Positions may be integer, numeric or character and give the same ID: whole
#' numbers are written in fixed notation (`100000`, never `1e+05`), and a
#' character position is read as a number (`"1e+05"`, `" 100000"` and
#' `"100000"` all give `100000`; `""`, `"."` and `"NA"` are missing). A
#' non-whole or negative position is an error, as are logical alleles
#' (`read.table()` reads a column of `T` alleles as `TRUE`). A missing
#' chromosome, position or allele gives `NA`.
#'
#' The ID is not strand-free: `A/G` and its complement `C/T` differ, and a
#' palindromic SNP (`A/T`, `C/G`) keeps its ID on either strand although its
#' alleles flip. After [liftover()], pass its `STRAND_b<to>` column as `strand`
#' to reverse-complement the alleles of rows on the `-` strand. Only `A/C/G/T/N`
#' alleles are complemented, and indels flipped that way are not
#' left-normalized.
#'
#' @param chrom Chromosome.
#' @param pos Base position: integer, numeric or character.
#' @param a1,a2 The two alleles, in any order.
#' @param strand Optional `"+"`/`"-"` per row; the alleles of `"-"` rows are
#'   reverse-complemented first. `NA` counts as `"+"`.
#' @param normalize Logical. Normalize chromosome names with
#'   [normalize.chrom()] (`chr1` -> `1`, `23` -> `X`) so IDs from sources with
#'   different conventions match.
#' @param sep Separator between the four parts.
#'
#' @return Character vector of IDs, one per row. Arguments of length 1 are
#'   recycled.
#' @examples
#' marker.id(1, c(100000, 100000L), c("G", "a"), c("A", "g"))
#' marker.id(c("chr2", "23"), c(5, 6), c("T", "C"), c("C", "A"), normalize = TRUE)
#' # alleles on the "-" strand of the new build, after liftover()
#' marker.id(1, 100, "A", "G", strand = "-")
#' @export
marker.id <- function(chrom, pos, a1, a2, strand = NULL, normalize = FALSE, sep = ":") {
  args <- list(chrom = chrom, pos = pos, a1 = a1, a2 = a2)
  if (!is.null(strand)) args$strand <- strand
  len <- lengths(args)
  n <- max(len)
  if (any(!(len %in% c(1L, n)))) {
    stop(sprintf("%s must have the same length (or length 1).", paste(names(args), collapse = ", ")), call. = FALSE)
  }
  x1 <- .gcanvas_marker_allele(a1, "a1")
  x2 <- .gcanvas_marker_allele(a2, "a2")
  p <- .gcanvas_marker_pos(pos)
  chr <- .gcanvas_marker_chrom(chrom)
  if (isTRUE(normalize)) chr <- .gcanvas_on_unique(chr, normalize.chrom)
  if (!n) return(character())
  chr <- rep_len(chr, n)
  p <- rep_len(p, n)
  x1 <- rep_len(x1, n)
  x2 <- rep_len(x2, n)

  if (!is.null(strand)) {
    s <- rep_len(as.character(strand), n)
    if (!all(s %in% c("+", "-", NA))) stop("strand must be \"+\", \"-\" or NA.", call. = FALSE)
    flip <- which(s == "-")
    x1[flip] <- .gcanvas_revcomp(x1[flip])
    x2[flip] <- .gcanvas_revcomp(x2[flip])
  }

  # byte order of the distinct alleles (order(method = "radix") is C-locale)
  u <- unique(c(x1, x2))
  r <- integer(length(u))
  r[order(u, method = "radix")] <- seq_along(u)
  swap <- which(r[data.table::chmatch(x2, u)] < r[data.table::chmatch(x1, u)])
  lo <- x1
  lo[swap] <- x2[swap]
  hi <- x2
  hi[swap] <- x1[swap]

  # one sprintf() call: several times faster than paste() on millions of rows
  sp <- gsub("%", "%%", as.character(sep)[1], fixed = TRUE)
  fmt <- paste0("%s", sp, if (is.integer(p)) "%d" else "%.0f", sp, "%s", sp, "%s")
  out <- sprintf(fmt, chr, p, lo, hi)
  out[is.na(chr) | is.na(p) | is.na(lo) | is.na(hi)] <- NA_character_
  out
}

# f(x) evaluated once per distinct value of a character vector.
.gcanvas_on_unique <- function(x, f) {
  u <- unique(x)
  f(u)[data.table::chmatch(x, u)]
}

.gcanvas_marker_allele <- function(x, what) {
  if (is.logical(x) && !all(is.na(x))) {
    stop(sprintf("%s is logical (an allele T read as TRUE?); read alleles as character.", what), call. = FALSE)
  }
  .gcanvas_on_unique(as.character(x), toupper)
}

# Chromosome names; whole numbers in fixed notation (`as.character(1)` is
# "1e+00" under a negative `scipen`).
.gcanvas_marker_chrom <- function(chrom) {
  if (!is.numeric(chrom)) return(as.character(chrom))
  u <- unique(chrom)
  v <- as.character(u)
  whole <- !is.na(u) & u == floor(u) & abs(u) <= .Machine$integer.max
  v[whole] <- as.character(as.integer(u[whole]))
  v[match(chrom, u)]
}

# Positions as whole numbers >= 0 (integer, or double above the integer
# range), NA where missing.
.gcanvas_marker_pos <- function(pos) {
  if (is.factor(pos) || inherits(pos, "integer64")) pos <- as.character(pos)
  if (is.character(pos)) {
    s <- trimws(pos)
    s[s %in% c("", ".", "NA")] <- NA_character_
    pos <- suppressWarnings(as.numeric(s))
    bad <- !is.na(s) & is.na(pos)
    if (any(bad)) {
      stop(sprintf("pos is not a number in %d row(s), e.g. \"%s\".", sum(bad), s[which(bad)[1]]), call. = FALSE)
    }
  }
  if (is.logical(pos) && all(is.na(pos))) pos <- as.integer(pos)
  if (!is.numeric(pos)) stop("pos must be integer, numeric or character.", call. = FALSE)
  if (is.integer(pos)) {
    if (any(pos < 0L, na.rm = TRUE)) stop("pos must be whole numbers >= 0.", call. = FALSE)
    return(pos)
  }
  pos[is.nan(pos)] <- NA_real_
  bad <- !is.na(pos) & (!is.finite(pos) | pos < 0 | pos != floor(pos))
  if (any(bad)) {
    stop(sprintf("pos must be whole numbers >= 0 (%d row(s), e.g. %s).", sum(bad), format(pos[which(bad)[1]], digits = 15)),
         call. = FALSE)
  }
  if (all(is.na(pos) | pos <= .Machine$integer.max)) as.integer(pos) else pos
}

# Reverse complement of A/C/G/T/N alleles; other alleles (`-`, `*`, `<DEL>`,
# `I`, `D`) are kept.
.gcanvas_revcomp <- function(x) {
  .gcanvas_on_unique(x, function(u) {
    dna <- !is.na(u) & grepl("^[ACGTN]+$", u)
    v <- chartr("ACGTN", "TGCAN", u[dna])
    long <- nchar(v) > 1L
    v[long] <- vapply(strsplit(v[long], "", fixed = TRUE), function(z) paste(rev(z), collapse = ""), character(1))
    u[dna] <- v
    u
  })
}
