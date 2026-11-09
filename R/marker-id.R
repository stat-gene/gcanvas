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
#' `"100000"` all give `100000`). A non-whole or negative position is an error,
#' as are logical alleles (`read.table()` reads a column of `T` alleles as
#' `TRUE`).
#'
#' A row gets `NA` when its chromosome, position or either allele is missing:
#' `NA`, or a string that is `""`, `"."` or `"NA"` once leading and trailing
#' spaces are removed (spaces are removed from every chromosome and allele), or
#' an allele `"0"` (PLINK's missing allele).
#'
#' Unusual alleles are reported with one warning per kind, giving the number of
#' rows and up to three of them as examples: the missing codes above (`0`, `.`,
#' `""`, `NA`; the row gets `NA`), `-` (the empty allele of an indel, as in
#' ANNOVAR), `*` (VCF spanning deletion), numeric codes (PLINK `1/2` or `1234`
#' coding), `I`/`D`, symbolic alleles such as `<DEL>`, `N`, and anything else
#' that is not `A/C/G/T`. Only the missing codes change the result; the other
#' alleles are kept but do not match the same variant written in `A/C/G/T` or as
#' a VCF-style indel (e.g. `AT/A`). Use [suppressWarnings()] to silence them.
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
  al <- .gcanvas_allele_check(x1, x2, function(i) sprintf(if (is.integer(p)) "%s:%d" else "%s:%.0f", chr[i], p[i]))
  x1 <- al[[1]]
  x2 <- al[[2]]

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

#' Sort the two alleles of each row
#'
#' Returns the two alleles of each row in upper case and in byte (C-locale)
#' order, the same order [marker.id()] writes them in. With these two columns
#' and the chromosome and position, two tables can be joined on the variant
#' whichever allele each lists first, without building ID strings: on ten
#' million rows such a join is several times faster and uses a fraction of the
#' memory of a join on [marker.id()] strings.
#'
#' A row with a missing allele (`NA`, or `""`, `"."`, `"NA"` or `"0"` as in
#' [marker.id()]) gets `NA` in both columns. data.table joins match `NA` keys
#' with each other, so drop those rows before joining. Unusual alleles are
#' reported as in [marker.id()], and logical alleles are an error
#' (`read.table()` reads a column of `T` alleles as `TRUE`).
#'
#' @inheritParams marker.id
#' @return A list of two character vectors, `lo` and `hi`, one element per row.
#'   Arguments of length 1 are recycled.
#' @examples
#' allele.sort(c("G", "a", "T"), c("A", "G", "C"))
#' # join two tables on position and alleles, whichever allele is listed first
#' x <- data.table::data.table(CHR = 1, POS = c(100L, 200L), A1 = c("A", "C"), A2 = c("G", "T"))
#' y <- data.table::data.table(CHR = 1, POS = c(100L, 200L), A1 = c("G", "C"), A2 = c("A", "G"))
#' x[, c("A_lo", "A_hi") := allele.sort(A1, A2)]
#' y[, c("A_lo", "A_hi") := allele.sort(A1, A2)]
#' y[x, on = c("CHR", "POS", "A_lo", "A_hi"), nomatch = NULL]
#' @export
allele.sort <- function(a1, a2, strand = NULL) {
  args <- list(a1 = a1, a2 = a2)
  if (!is.null(strand)) args$strand <- strand
  len <- lengths(args)
  n <- max(len)
  if (any(!(len %in% c(1L, n)))) {
    stop(sprintf("%s must have the same length (or length 1).", paste(names(args), collapse = ", ")), call. = FALSE)
  }
  al <- .gcanvas_allele_check(rep_len(.gcanvas_marker_allele(a1, "a1"), n), rep_len(.gcanvas_marker_allele(a2, "a2"), n))
  x1 <- al[[1]]
  x2 <- al[[2]]
  if (!is.null(strand)) {
    s <- rep_len(as.character(strand), n)
    if (!all(s %in% c("+", "-", NA))) stop("strand must be \"+\", \"-\" or NA.", call. = FALSE)
    flip <- which(s == "-")
    x1[flip] <- .gcanvas_revcomp(x1[flip])
    x2[flip] <- .gcanvas_revcomp(x2[flip])
  }
  miss <- is.na(x1) | is.na(x2)
  x1[miss] <- NA_character_
  x2[miss] <- NA_character_
  u <- unique(c(x1, x2))
  r <- integer(length(u))
  r[order(u, method = "radix")] <- seq_along(u)
  swap <- which(r[data.table::chmatch(x2, u)] < r[data.table::chmatch(x1, u)])
  lo <- x1
  lo[swap] <- x2[swap]
  hi <- x2
  hi[swap] <- x1[swap]
  list(lo = lo, hi = hi)
}

# Strings that count as missing in chromosome, position and allele columns;
# alleles also have PLINK's "0".
.gcanvas_marker_na <- c("", ".", "NA")
.gcanvas_allele_na <- c(.gcanvas_marker_na, "0")

# Kinds of unusual allele reported by .gcanvas_allele_check(), in that order.
.gcanvas_allele_kinds <- c(
  missing = "missing allele code (0, ., \"\", NA), set to NA",
  dash = "allele \"-\" (empty allele of an indel), kept; VCF-style indels such as AT/A do not match it",
  star = "allele \"*\" (VCF spanning deletion), kept",
  number = "numeric allele code (PLINK 1/2 or 1234 coding?), kept; it does not match A/C/G/T",
  indel = "allele I/D, kept; VCF-style indels do not match it",
  symbolic = "symbolic allele (e.g. <DEL>), kept",
  n = "allele with N (unknown base), kept",
  other = "allele not made of A/C/G/T, kept"
)

# Kind of each (upper-cased, trimmed) allele; NA for A/C/G/T and missing values.
.gcanvas_allele_kind <- function(u) {
  k <- rep(NA_character_, length(u))
  set <- function(hit, kind) k[is.na(k) & !is.na(u) & hit] <<- kind
  set(u %in% .gcanvas_allele_na, "missing")
  set(u == "-", "dash")
  set(u == "*", "star")
  set(grepl("^[0-9]+$", u), "number")
  set(u %in% c("I", "D"), "indel")
  set(grepl("^<.*>$", u), "symbolic")
  set(grepl("^[ACGTN]+$", u), "n")
  set(!grepl("^[ACGT]+$", u), "other")
  k[k %in% "n" & !grepl("N", u, fixed = TRUE)] <- NA_character_
  k
}

# Warns once per kind of unusual allele with up to three rows as examples
# (`label(i)` names rows `i`, e.g. "1:100"), and returns the two alleles with
# the missing codes set to NA.
.gcanvas_allele_check <- function(x1, x2, label = function(i) paste("row", i)) {
  u <- unique(c(x1, x2))
  kind <- .gcanvas_allele_kind(u)
  if (!all(is.na(kind))) {
    k1 <- kind[data.table::chmatch(x1, u)]
    k2 <- kind[data.table::chmatch(x2, u)]
    show <- function(a) ifelse(is.na(a), "NA", ifelse(nzchar(a), a, "\"\""))
    for (k in intersect(names(.gcanvas_allele_kinds), kind)) {
      rows <- which(k1 %in% k | k2 %in% k)
      ex <- rows[seq_len(min(3L, length(rows)))]
      warning(sprintf("%s: %s row(s), e.g. %s", .gcanvas_allele_kinds[[k]], format(length(rows), big.mark = ","),
                      paste(sprintf("%s %s/%s", label(ex), show(x1[ex]), show(x2[ex])), collapse = ", ")), call. = FALSE)
    }
    x1[k1 %in% "missing"] <- NA_character_
    x2[k2 %in% "missing"] <- NA_character_
  }
  list(x1, x2)
}

# Spaces trimmed, missing-value strings set to NA.
.gcanvas_marker_clean <- function(x) {
  x <- trimws(x)
  x[x %in% .gcanvas_marker_na] <- NA_character_
  x
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
  .gcanvas_on_unique(as.character(x), function(u) trimws(toupper(u)))
}

# Chromosome names; whole numbers in fixed notation (`as.character(1)` is
# "1e+00" under a negative `scipen`).
.gcanvas_marker_chrom <- function(chrom) {
  if (!is.numeric(chrom)) return(.gcanvas_on_unique(as.character(chrom), .gcanvas_marker_clean))
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
    s <- .gcanvas_marker_clean(pos)
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
