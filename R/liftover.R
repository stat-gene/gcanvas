# Coordinate liftover between genome builds. The default path reads the chain
# file in R (no external binary); `method = "binary"` runs UCSC `liftOver
# -multiple` instead. Both collect every hit of a position and pass them through
# `.gcanvas_lift_resolve()`, so a new position is attached only when it lies on
# the input's own chromosome at a single place, and every row gets a status.

# Canonical assembly key for a build argument (accepts 38/hg38/GRCh38,
# 37/19/hg19/GRCh37, 36/18/hg18/NCBI36). Returns NA for anything unknown.
.gcanvas_build_key <- function(b) {
  s <- toupper(trimws(as.character(b)[1]))
  if (s %in% c("38", "HG38", "GRCH38")) return("GRCh38")
  if (s %in% c("37", "19", "HG19", "GRCH37")) return("GRCh37")
  if (s %in% c("36", "18", "HG18", "NCBI36")) return("hg18")
  NA_character_
}

# Numeric build label used for output column suffixes (POS_b38 etc.).
.gcanvas_build_num <- function(b) {
  k <- .gcanvas_build_key(b)
  if (is.na(k)) return(NA_character_)
  switch(k, GRCh38 = "38", GRCh37 = "37", hg18 = "18", NA_character_)
}

# Candidate chain basenames for a (from -> to) build pair, in preference order.
.gcanvas_chain_candidates <- function(from_key, to_key) {
  pair <- paste(from_key, to_key, sep = "=>")
  switch(pair,
    "GRCh37=>GRCh38" = c("GRCh37_to_GRCh38.chain.gz", "hg19ToHg38.over.chain.gz"),
    "GRCh38=>GRCh37" = c("GRCh38_to_GRCh37.chain.gz", "hg38ToHg19.over.chain.gz"),
    "hg18=>GRCh37"   = c("hg18ToHg19.over.chain.gz", "NCBI36_to_GRCh37.chain.gz"),
    character()
  )
}

# UCSC chain files that chain.download() can fetch (not bundled: UCSC allows
# its chain files for non-commercial use only).
.gcanvas_ucsc_chains <- list(
  "hg18=>GRCh37" = list(
    file = "hg18ToHg19.over.chain.gz",
    url = "https://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz",
    readme = "https://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/README.txt",
    md5 = "34e6ad1eaeb83af301617f2460641933"
  )
)

# Directories searched by `liftover.chain = "auto"`: the chains bundled in
# `inst/extdata`, then the user data directory that chain.download() fills.
.gcanvas_chain_dirs <- function() {
  d <- c(system.file("extdata", package = "gcanvas"), tools::R_user_dir("gcanvas", which = "data"))
  d[nzchar(d)]
}

.gcanvas_chain_is_auto <- function(liftover.chain) {
  chain0 <- as.character(liftover.chain)[1]
  is.null(liftover.chain) || length(liftover.chain) == 0L ||
    is.na(chain0) || !nzchar(chain0) || toupper(chain0) == "AUTO"
}

# Resolve a chain file path: an explicit `liftover.chain` must exist; "auto"
# picks the first candidate for the build pair found in the chain dirs. Errors
# (with what to do next) when nothing is found.
.gcanvas_resolve_chain <- function(from, to, liftover.chain) {
  if (!.gcanvas_chain_is_auto(liftover.chain)) {
    chain0 <- as.character(liftover.chain)[1]
    if (!file.exists(chain0)) stop(sprintf("chain file not found: %s", chain0), call. = FALSE)
    return(abs_path(chain0))
  }
  from_key <- .gcanvas_build_key(from)
  to_key <- .gcanvas_build_key(to)
  if (is.na(from_key) || is.na(to_key)) {
    stop(sprintf("unknown build %s -> %s (use 37/38/18, hg19/hg38/hg18 or GRCh37/GRCh38), or pass `liftover.chain`.",
                 as.character(from)[1], as.character(to)[1]), call. = FALSE)
  }
  names <- .gcanvas_chain_candidates(from_key, to_key)
  if (!length(names)) {
    stop(sprintf("no chain file is known for %s -> %s; pass `liftover.chain`.", from_key, to_key), call. = FALSE)
  }
  for (d in .gcanvas_chain_dirs()) {
    for (nm in names) {
      p <- file.path(d, nm)
      if (file.exists(p)) return(abs_path(p))
    }
  }
  hint <- if (!is.null(.gcanvas_ucsc_chains[[paste(from_key, to_key, sep = "=>")]])) {
    sprintf(" Download it once with chain.download(from = %s, to = %s) (UCSC chain: non-commercial use only),",
            as.character(from)[1], as.character(to)[1])
  } else {
    ""
  }
  stop(sprintf("no chain file found for %s -> %s (looked for %s).%s or pass `liftover.chain`.",
               from_key, to_key, paste(names, collapse = ", "), hint), call. = FALSE)
}

# Locate the UCSC liftOver binary for `method = "binary"`. A configured value
# wins -- the `liftover` argument, then `options(gcanvas.liftover)`, then the
# `GCANVAS_LIFTOVER` environment variable -- and must resolve (a file path or a
# command on PATH). Otherwise `liftOver` is looked up on PATH, then in the bin
# directory of the active conda env and of the conda env R itself runs from.
.gcanvas_resolve_liftover_bin <- function(liftover = NULL) {
  conf <- c(
    as.character(liftover %||% "")[1],
    as.character(getOption("gcanvas.liftover", ""))[1],
    Sys.getenv("GCANVAS_LIFTOVER", unset = "")
  )
  conf <- conf[!is.na(conf) & nzchar(conf)]
  if (length(conf)) {
    cmd <- conf[1]
    if (file.exists(cmd) && !dir.exists(cmd)) return(abs_path(cmd))
    w <- Sys.which(cmd)
    return(if (nzchar(w)) unname(w) else structure(NA_character_, configured = cmd))
  }
  w <- Sys.which("liftOver")
  if (nzchar(w)) return(unname(w))
  envs <- unique(c(Sys.getenv("CONDA_PREFIX", unset = ""), dirname(dirname(R.home()))))
  envs <- envs[nzchar(envs) & dir.exists(file.path(envs, "conda-meta"))]
  for (e in envs) {
    p <- file.path(e, "bin", "liftOver")
    if (file.exists(p)) return(abs_path(p))
  }
  NA_character_
}

.gcanvas_lift_std <- c(as.character(1:22), "X", "Y", "MT")

# TRUE when a chain spells the standard chromosomes UCSC-style (chr1, chrM):
# none of its target names is Ensembl-style (1, MT) and one is UCSC-style.
# (Ensembl chains also carry a "chrM" alias next to "MT".)
.gcanvas_chain_ucsc <- function(tName) {
  !any(tName %in% .gcanvas_lift_std) && any(tName %in% .gcanvas_lift_seqname(.gcanvas_lift_std, TRUE))
}

# Chromosome naming convention of a chain file from its header lines:
# "ucsc" (chr-prefixed) or "plain" (Ensembl-style, no prefix).
.gcanvas_chain_seqstyle <- function(chain_path) {
  con <- gzfile(chain_path, open = "rt")
  x <- tryCatch(readLines(con, warn = FALSE), finally = close(con))
  hdr <- x[startsWith(x, "chain")]
  tName <- vapply(strsplit(trimws(hdr), "[ \t]+"), function(f) f[3], character(1))
  if (.gcanvas_chain_ucsc(tName)) "ucsc" else "plain"
}

# Standard chromosome (1-22, X, Y, MT) of a name, or NA. Strips a chr prefix
# and maps 23/24/26/M; 25/XY (the X pseudo-autosomal region, in X
# coordinates) is looked up on X. Patch, alt and unplaced contigs (HG*_PATCH,
# CHR_*, *_alt, GL*, KI*, chrUn_*) are NA.
.gcanvas_lift_chr <- function(x) {
  x <- sub("^CHR", "", toupper(trimws(as.character(x))))
  x[x == "23"] <- "X"
  x[x == "24"] <- "Y"
  x[x %in% c("25", "XY")] <- "X"
  x[x %in% c("26", "M")] <- "MT"
  x[!(x %in% .gcanvas_lift_std)] <- NA_character_
  x
}

# Standard chromosome names as a chain file spells them.
.gcanvas_lift_seqname <- function(chr, ucsc) {
  if (!isTRUE(ucsc)) return(chr)
  out <- paste0("chr", chr)
  out[chr %in% "MT"] <- "chrM"
  out
}

.gcanvas_chain_cache <- new.env(parent = emptyenv())

# Gap-free aligned blocks of a chain file whose target is a standard
# chromosome, cached per file for the session. Columns: tName (standard name),
# tStart/tEnd (0-based half-open on the target), qName (as in the chain),
# qStd (its standard name, NA for contigs), qSize, qStrand, qStart (0-based on
# the query strand) and chain (header order; Ensembl chain ids repeat).
.gcanvas_chain_blocks <- function(path) {
  info <- file.info(path)
  key <- paste(normalizePath(path, winslash = "/", mustWork = TRUE), info$size, as.numeric(info$mtime), sep = "|")
  hit <- .gcanvas_chain_cache[[key]]
  if (!is.null(hit)) return(hit)
  blocks <- .gcanvas_chain_read(path)
  assign(key, blocks, envir = .gcanvas_chain_cache)
  blocks
}

.gcanvas_chain_read <- function(path) {
  bad <- function(what) stop(sprintf("malformed chain file (%s): %s", what, path), call. = FALSE)
  con <- gzfile(path, open = "rt")
  x <- tryCatch(readLines(con, warn = FALSE), finally = close(con))
  x <- trimws(x)
  x <- x[nzchar(x) & !startsWith(x, "#")]
  hdr <- startsWith(x, "chain")
  if (!length(x) || !hdr[1]) bad("no leading 'chain' header")

  # chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
  h <- data.table::tstrsplit(x[hdr], "[ \t]+")
  if (length(h) < 12L) bad("short header line")
  tName <- h[[3]]
  qName <- h[[8]]
  qStrand <- h[[10]]
  num <- function(v) suppressWarnings(as.numeric(v))
  tStart0 <- num(h[[6]]); tEnd0 <- num(h[[7]])
  qSize <- num(h[[9]]); qStart0 <- num(h[[11]]); qEnd0 <- num(h[[12]])
  if (anyNA(c(tStart0, tEnd0, qSize, qStart0, qEnd0)) || !all(h[[5]] == "+") || !all(qStrand %in% c("+", "-"))) {
    bad("header fields")
  }

  # "size dt dq" lines, then a last "size" line per chain
  chain <- cumsum(hdr)[!hdr]
  if (!length(chain) || any(tabulate(chain, nbins = length(tName)) == 0L)) bad("chain without blocks")
  b <- data.table::tstrsplit(x[!hdr], "[ \t]+", fill = NA_character_)
  size <- num(b[[1]])
  dt <- if (length(b) >= 3L) num(b[[2]]) else rep(NA_real_, length(size))
  dq <- if (length(b) >= 3L) num(b[[3]]) else rep(NA_real_, length(size))
  last <- c(chain[-1L] != chain[-length(chain)], TRUE)
  if (anyNA(size) || anyNA(dt[!last]) || anyNA(dq[!last]) || !all(is.na(dt[last]))) bad("block lines")
  dt[last] <- 0
  dq[last] <- 0

  # block starts: chain start + sizes and gaps of the earlier blocks of the chain
  step_t <- size + dt
  step_q <- size + dq
  off_t <- cumsum(step_t) - step_t
  off_q <- cumsum(step_q) - step_q
  first <- match(chain, chain)
  tS <- tStart0[chain] + off_t - off_t[first]
  qS <- qStart0[chain] + off_q - off_q[first]
  if (any(tS[last] + size[last] != tEnd0[chain[last]]) || any(qS[last] + size[last] != qEnd0[chain[last]])) {
    bad("blocks do not add up to the header range")
  }

  t_std <- .gcanvas_lift_std[match(tName, .gcanvas_lift_seqname(.gcanvas_lift_std, .gcanvas_chain_ucsc(tName)))]
  keep <- !is.na(t_std[chain])
  i <- chain[keep]
  blocks <- data.table::data.table(
    tName = t_std[i],
    tStart = as.integer(tS[keep]),
    tEnd = as.integer(tS[keep] + size[keep]),
    qName = qName[i],
    qStd = .gcanvas_lift_chr(qName)[i],
    qSize = as.integer(qSize[i]),
    qStrand = qStrand[i],
    qStart = as.integer(qS[keep]),
    chain = i
  )
  data.table::setkeyv(blocks, c("tName", "tStart"))
  blocks
}

.gcanvas_lift_hits_empty <- function() {
  data.table::data.table(row = integer(), qName = character(), qStd = character(), pos = integer(), strand = character())
}

# Every hit of the 0-based positions `p0` on standard chromosomes `chr` (NA =
# not looked up): the chain blocks covering the base. Returns input row, new
# chromosome (as named in the chain), its standard name (NA for contigs),
# 1-based new position and strand. A base on a "-" query block at q0 (0-based
# on the reverse strand) is qSize - q0 on the + strand (1-based).
.gcanvas_lift_hits_chain <- function(chr, p0, blocks) {
  q <- data.table::data.table(row = seq_along(chr), tName = chr, p0 = p0)
  q <- q[!is.na(tName) & !is.na(p0)]
  if (!nrow(q) || !nrow(blocks)) return(.gcanvas_lift_hits_empty())
  h <- blocks[q, on = c("tName", "tStart<=p0", "tEnd>p0"), nomatch = NULL, allow.cartesian = TRUE,
              list(row = i.row, qName = x.qName, qStd = x.qStd, strand = x.qStrand,
                   q0 = x.qStart + (i.p0 - x.tStart), qSize = x.qSize)]
  h[, pos := data.table::fifelse(strand == "-", qSize - q0, q0 + 1L)]
  h[, list(row, qName, qStd, pos, strand)]
}

# The same hits from UCSC `liftOver -multiple` (BED named by row number).
.gcanvas_lift_hits_binary <- function(chr, p0, chain, bin, rm.tmp = TRUE) {
  rows <- which(!is.na(chr) & !is.na(p0))
  if (!length(rows)) return(.gcanvas_lift_hits_empty())
  ucsc <- identical(.gcanvas_chain_seqstyle(chain), "ucsc")
  f <- paste0(tempfile("gcanvas-liftover-"), c(".bed", ".lifted.bed", ".unlifted.bed"))
  if (isTRUE(rm.tmp)) on.exit(unlink(f), add = TRUE)
  bed <- data.table::data.table(.gcanvas_lift_seqname(chr[rows], ucsc), p0[rows], p0[rows] + 1L, rows, 0L, "+")
  data.table::fwrite(bed, f[1], sep = "\t", col.names = FALSE, quote = FALSE)
  log <- suppressWarnings(system2(bin, c("-multiple", shQuote(f[1]), shQuote(chain), shQuote(f[2]), shQuote(f[3])),
                                  stdout = TRUE, stderr = TRUE))
  rc <- attr(log, "status")
  if ((!is.null(rc) && rc != 0L) || !file.exists(f[2])) {
    stop(sprintf("liftOver failed (exit status %s): %s", rc %||% "?",
                 paste(utils::tail(log, 5L), collapse = " | ")), call. = FALSE)
  }
  if (!file.size(f[2])) return(.gcanvas_lift_hits_empty())
  # chrom start end name serial strand
  out <- data.table::fread(f[2], header = FALSE, sep = "\t", colClasses = "character")
  data.table::data.table(
    row = as.integer(out[[4]]),
    qName = out[[1]],
    qStd = .gcanvas_lift_chr(out[[1]]),
    pos = as.integer(out[[3]]),
    strand = out[[ncol(out)]]
  )
}

# Rows from hits: identical hits are merged; exactly one hit on the row's own
# chromosome gives "ok" with its position and strand. Otherwise no hit is
# "unmapped", several distinct same-chromosome hits "multiple", hits only on
# other standard chromosomes "other_chrom", and hits only on patch / alt /
# unplaced contigs "contig".
.gcanvas_lift_resolve <- function(hits, chr, bad) {
  n <- length(chr)
  hits <- unique(hits, by = c("row", "qName", "pos", "strand"))
  same <- !is.na(hits$qStd) & hits$qStd == chr[hits$row]
  other <- !is.na(hits$qStd) & !same
  n_same <- tabulate(hits$row[same], nbins = n)
  status <- rep("unmapped", n)
  status[tabulate(hits$row, nbins = n) > 0L] <- "contig"
  status[tabulate(hits$row[other], nbins = n) > 0L] <- "other_chrom"
  status[n_same > 1L] <- "multiple"
  status[n_same == 1L] <- "ok"
  status[bad] <- "bad_input"
  pos <- rep(NA_integer_, n)
  strand <- rep(NA_character_, n)
  one <- same & n_same[hits$row] == 1L
  pos[hits$row[one]] <- hits$pos[one]
  strand[hits$row[one]] <- hits$strand[one]
  data.table::data.table(pos = pos, strand = strand, status = status)
}

# Lift 1-based positions `pos` on chromosomes `chrom` from `from` to `to`.
# Returns list(result = one row per input with pos/strand/status, chain,
# method). Errors when the liftover cannot run at all (no chain, no binary).
.gcanvas_lift_core <- function(chrom, pos, from, to, method = "chain", liftover = NULL,
                               liftover.chain = "auto", rm.tmp = TRUE) {
  chr <- .gcanvas_lift_chr(chrom)
  p <- suppressWarnings(as.numeric(pos))
  bad <- is.na(chr) | !is.finite(p) | p < 1 | p > .Machine$integer.max | p != floor(p)
  chr[bad] <- NA_character_
  p0 <- rep(NA_integer_, length(p))
  p0[!bad] <- as.integer(p[!bad]) - 1L

  from_key <- .gcanvas_build_key(from)
  if (!is.na(from_key) && identical(from_key, .gcanvas_build_key(to)) && .gcanvas_chain_is_auto(liftover.chain)) {
    res <- data.table::data.table(
      pos = data.table::fifelse(bad, NA_integer_, p0 + 1L),
      strand = data.table::fifelse(bad, NA_character_, "+"),
      status = data.table::fifelse(bad, "bad_input", "ok")
    )
    return(list(result = res, chain = NA_character_, method = "identity"))
  }

  chain <- .gcanvas_resolve_chain(from, to, liftover.chain)
  if (identical(method, "binary")) {
    bin <- .gcanvas_resolve_liftover_bin(liftover)
    if (is.na(bin) && !is.null(attr(bin, "configured"))) {
      stop(sprintf("liftOver binary not found: %s", attr(bin, "configured")), call. = FALSE)
    }
    if (is.na(bin)) {
      stop("method = \"binary\" needs the UCSC liftOver binary: pass `liftover`, set options(gcanvas.liftover) or ",
           "GCANVAS_LIFTOVER, or put liftOver on PATH.", call. = FALSE)
    }
    hits <- .gcanvas_lift_hits_binary(chr, p0, chain, bin, rm.tmp = rm.tmp)
  } else {
    hits <- .gcanvas_lift_hits_chain(chr, p0, .gcanvas_chain_blocks(chain))
  }
  list(result = .gcanvas_lift_resolve(hits, chr, bad), chain = chain, method = method)
}

.gcanvas_lift_status_levels <- c("ok", "unmapped", "other_chrom", "contig", "multiple", "bad_input", "failed")

# Build label for messages: the canonical key when known, else as given.
.gcanvas_build_label <- function(b) {
  k <- .gcanvas_build_key(b)
  if (is.na(k)) as.character(b)[1] else k
}

# "ok=9,950 | unmapped=40 | contig=10" (non-zero statuses, fixed order).
.gcanvas_lift_summary <- function(status) {
  tab <- table(factor(status, levels = .gcanvas_lift_status_levels))
  tab <- tab[tab > 0L]
  if (!length(tab)) return("no rows")
  paste(sprintf("%s=%s", names(tab), format(as.integer(tab), big.mark = ",", trim = TRUE)), collapse = " | ")
}

#' Lift variant coordinates between genome builds
#'
#' Translates each row's chromosome and position from build `from` to build
#' `to`, and attaches the new position only when it can be trusted. By default
#' the chain file is read directly in R, so no external program is needed and
#' it works on every OS; `method = "binary"` runs the UCSC `liftOver`
#' executable with `-multiple` instead and applies the same rules to its hits.
#' `from`/`to` accept `38`/`hg38`/`GRCh38`, `37`/`19`/`hg19`/`GRCh37`, and
#' `36`/`18`/`hg18`/`NCBI36`.
#'
#' @section Rules:
#' 1. **One output row per input row.** Row count and order are unchanged and
#'    results are attached by row, never by SNP id, so duplicated or missing
#'    ids cannot mix positions up.
#' 2. **Same chromosome only.** A new position is kept only on the input's own
#'    chromosome. Names are compared after normalization (`chr` prefix;
#'    `23`/`24`/`26`/`M` for X/Y/MT; `25`/`XY`, the pseudo-autosomal region,
#'    is looked up on X). A position that lands on another chromosome, or only
#'    on a patch, alt or unplaced contig (`HG*_PATCH`, `CHR_*`, `*_alt`, `GL*`,
#'    `KI*`, `chrUn_*`), is `NA`; a contig coordinate is never attached as a
#'    position on the original chromosome.
#' 3. **Several hits** (what `liftOver` reports as "Duplicated in new"): every
#'    hit is examined. Hits at the same chromosome and position are merged;
#'    exactly one hit on the input's chromosome is used, none or two or more
#'    distinct ones give `NA`.
#' 4. **Every row gets a status** (`lift_status`), and the counts per status
#'    are summarized at the end (unless `silent`).
#'
#' `lift_status` is one of
#' \describe{
#'   \item{`ok`}{lifted to a single position on the same chromosome.}
#'   \item{`unmapped`}{no aligned block covers the base (deleted in the new
#'     build, or the chromosome is not in the chain).}
#'   \item{`other_chrom`}{hits only on other standard chromosomes.}
#'   \item{`contig`}{hits only on patch, alt or unplaced contigs.}
#'   \item{`multiple`}{two or more distinct positions on the same chromosome.}
#'   \item{`bad_input`}{chromosome not 1-22/X/Y/MT, or position not a positive
#'     whole number.}
#'   \item{`failed`}{the liftover could not run at all (no chain file,
#'     unreadable chain, missing or failing binary); a warning is raised.}
#' }
#'
#' @section Chain files:
#' With `liftover.chain = "auto"` the chain is chosen by build pair. GRCh37
#' <-> GRCh38 use the Ensembl assembly-mapping chains bundled with gcanvas
#' (Ensembl places no restrictions on the use of its data). hg18 -> hg19 uses
#' the UCSC chain, which is not bundled because UCSC allows its liftOver chain
#' files for non-commercial use only; fetch it once with [chain.download()].
#' Any other chain file can be given as a path. A chain is parsed once per
#' session and cached. When `from` and `to` are the same build, positions are
#' returned unchanged (status `ok`).
#'
#' @section liftOver binary:
#' Only `method = "binary"` uses the UCSC `liftOver` executable; it gives the
#' same result as the default and is mainly a cross-check. It is not bundled:
#' UCSC licenses it free for non-commercial use only (commercial use needs a
#' UCSC license), and excludes it from the tools that may be redistributed via
#' package managers. If you are eligible, `conda install bioconda::ucsc-liftover`
#' or a download from <https://hgdownload.soe.ucsc.edu/admin/exe/> works. It
#' is found through the `liftover` argument, `options(gcanvas.liftover)`, the
#' `GCANVAS_LIFTOVER` environment variable, `PATH`, and the `bin` directory
#' of the active conda environment, in that order.
#'
#' @param df A `data.frame`/`data.table` of variants.
#' @param from,to Source and destination genome builds (see Description).
#' @param liftover Path to the `liftOver` executable (or a command name on
#'   `PATH`) for `method = "binary"`; ignored otherwise.
#' @param liftover.chain Explicit chain file path, or `"auto"` (the default;
#'   `NULL` is treated the same) to choose one by build pair.
#' @param snp.col Ignored; kept for backward compatibility. Results are
#'   attached by row, so no SNP id column is needed.
#' @param chrom.col,pos.col Column names in `df` for the chromosome and the
#'   1-based position.
#' @param method `"auto"` (the default) and `"chain"` read the chain file in R;
#'   `"binary"` runs UCSC `liftOver -multiple`.
#' @param strict Logical. When the liftover cannot run at all, raise an error
#'   instead of a warning. Per-row failures never raise; they are reported in
#'   `lift_status`.
#' @param rm.tmp Logical. With `method = "binary"`, remove the temporary BED
#'   files (under [tempdir()]) when done.
#' @param silent Logical. Suppress the summary note. Notes go to stderr via
#'   [message()], so [suppressMessages()] silences them as well.
#'
#' @return `df` (a `data.table` if `df` was one, otherwise a `data.frame` with
#'   its row names) with these columns added or replaced, `<to>` being the
#'   numeric build label (`38`, `37`):
#'   \describe{
#'     \item{`<pos.col>_b<to>`}{integer new position; `NA` unless
#'       `lift_status` is `"ok"`.}
#'     \item{`<chrom.col>_b<to>`}{the input chromosome, in its own notation,
#'       when `"ok"`; `NA` otherwise.}
#'     \item{`STRAND_b<to>`}{`"+"` or `"-"` when `"ok"`. `"-"` means the base
#'       is reverse-complemented in the new build, so alleles must be
#'       complemented.}
#'     \item{`lift_status`}{see the list above.}
#'   }
#'
#' @seealso [chain.download()]
#' @examples
#' df <- data.frame(
#'   SNP = c("rs1", "rs2", "rs3", "rs3"),
#'   CHR = c(1, 1, 1, 1),
#'   POS = c(1005000, 143185000, 50005000, 60005000)
#' )
#' # rs2 lands only on a GRCh37 patch contig, so it gets NA and "contig";
#' # the two rs3 rows keep their own positions.
#' liftover(df, from = 38, to = 37)
#' @export
liftover <- function(df, from = 37, to = 38,
                     liftover = NULL,
                     liftover.chain = "auto",
                     snp.col = "SNP", chrom.col = "CHR", pos.col = "POS",
                     method = c("auto", "chain", "binary"),
                     strict = FALSE,
                     rm.tmp = TRUE, silent = FALSE) {
  require_pkg("data.table")
  method <- match.arg(method)
  if (identical(method, "auto")) method <- "chain"

  to_label <- .gcanvas_build_num(to)
  if (is.na(to_label)) to_label <- as.character(to)[1]

  is_dt <- data.table::is.data.table(df)
  rn0 <- if (!is_dt) rownames(df) else NULL
  dt <- if (is_dt) data.table::copy(df) else data.table::as.data.table(df)
  chr_use <- .gcanvas_resolve_colname(names(dt), chrom.col, aliases = character(), required = TRUE, arg_label = "chrom.col")
  pos_use <- .gcanvas_resolve_colname(names(dt), pos.col, aliases = character(), required = TRUE, arg_label = "pos.col")
  finalize <- function(xdt) {
    if (is_dt) return(xdt)
    out <- as.data.frame(xdt, stringsAsFactors = FALSE)
    if (!is.null(rn0) && length(rn0) == nrow(out)) rownames(out) <- rn0
    out
  }

  run <- tryCatch(
    .gcanvas_lift_core(dt[[chr_use]], dt[[pos_use]], from = from, to = to, method = method,
                       liftover = liftover, liftover.chain = liftover.chain, rm.tmp = rm.tmp),
    error = function(e) e
  )
  if (inherits(run, "error")) {
    msg <- sprintf("liftover %s -> %s failed: %s", .gcanvas_build_label(from), .gcanvas_build_label(to), conditionMessage(run))
    if (isTRUE(strict)) stop(msg, call. = FALSE)
    warning(msg, call. = FALSE)
    n <- nrow(dt)
    run <- list(
      result = data.table::data.table(pos = rep(NA_integer_, n), strand = rep(NA_character_, n), status = rep("failed", n)),
      chain = NA_character_, method = method
    )
  }
  res <- run$result

  chr_new <- dt[[chr_use]]
  if (is.factor(chr_new)) chr_new <- as.character(chr_new)
  chr_new[res$status != "ok"] <- NA
  data.table::set(dt, j = paste0(pos_use, "_b", to_label), value = res$pos)
  data.table::set(dt, j = paste0(chr_use, "_b", to_label), value = chr_new)
  data.table::set(dt, j = paste0("STRAND_b", to_label), value = res$strand)
  data.table::set(dt, j = "lift_status", value = res$status)

  if (!isTRUE(silent)) {
    how <- if (is.na(run$chain)) run$method else sprintf("%s, %s", run$method, basename(run$chain))
    .gcanvas_note("gcanvas::liftover", sprintf("%s -> %s (%s): n=%s | %s",
                                               .gcanvas_build_label(from), .gcanvas_build_label(to), how,
                                               format(nrow(dt), big.mark = ","), .gcanvas_lift_summary(res$status)))
  }
  finalize(dt)
}

# Lifted positions and status (one row per element of `pos`) on chromosome
# `chrom`; regional() lifts its variants, lead(s) and window ends in one call.
# Errors when the liftover cannot run.
liftover_positions <- function(chrom, pos, from, to, liftover = NULL, liftover.chain = "auto", method = "auto") {
  pos <- .gcanvas_as_num2(pos)
  method <- match.arg(method, c("auto", "chain", "binary"))
  if (identical(method, "auto")) method <- "chain"
  .gcanvas_lift_core(rep_len(as.character(chrom), length(pos)), pos, from = from, to = to, method = method,
                     liftover = liftover, liftover.chain = liftover.chain)$result
}

#' Download a UCSC liftOver chain file
#'
#' Fetches a UCSC chain file that gcanvas does not bundle into the directory
#' where [liftover()] looks with `liftover.chain = "auto"`. Currently hg18 ->
#' hg19 (`hg18ToHg19.over.chain.gz`); GRCh37 <-> GRCh38 chains ship with the
#' package.
#'
#' UCSC makes its liftOver chain files the exception to its otherwise free
#' data: they are free only for non-commercial use (independent researchers
#' and nonprofit organizations), any other use needs a commercial license, and
#' downloading or using them means accepting the UCSC EULA. See
#' <https://genome.ucsc.edu/license/> and the README next to the chain file.
#' Nothing is downloaded until you accept: pass `accept.license = TRUE`, or
#' answer the prompt in an interactive session.
#'
#' @param from,to Builds of the chain (aliases as in [liftover()]).
#' @param dir Destination directory. Default: `tools::R_user_dir("gcanvas",
#'   "data")`, which `liftover.chain = "auto"` searches.
#' @param accept.license Logical. `TRUE` confirms that your use is
#'   non-commercial and that you accept the UCSC terms.
#' @param overwrite Logical. Download again even if the file exists.
#'
#' @return Invisibly, the path of the chain file.
#' @examples
#' \dontrun{
#' chain.download(18, 19, accept.license = TRUE)
#' liftover(df, from = 18, to = 19)
#' }
#' @export
chain.download <- function(from = 18, to = 19, dir = NULL, accept.license = FALSE, overwrite = FALSE) {
  from_key <- .gcanvas_build_key(from)
  to_key <- .gcanvas_build_key(to)
  ent <- .gcanvas_ucsc_chains[[paste(from_key, to_key, sep = "=>")]]
  if (is.null(ent)) {
    stop(sprintf("no downloadable chain for %s -> %s (available: hg18 -> hg19). GRCh37 <-> GRCh38 chains are bundled.",
                 as.character(from)[1], as.character(to)[1]), call. = FALSE)
  }
  dir <- dir %||% tools::R_user_dir("gcanvas", which = "data")
  dest <- file.path(dir, ent$file)
  if (file.exists(dest) && !isTRUE(overwrite)) return(invisible(abs_path(dest)))

  notice <- paste0(
    "UCSC liftOver chain files are free for non-commercial use only (independent researchers and nonprofit ",
    "organizations); any other use needs a UCSC commercial license. Downloading or using them means accepting ",
    "the UCSC EULA (https://genome.ucsc.edu/license/EULA.pdf). Terms: ", ent$readme
  )
  if (!isTRUE(accept.license)) {
    ok <- interactive() && isTRUE(utils::askYesNo(paste0(notice, "\nDownload ", ent$file, " under these terms?"), default = FALSE))
    if (!ok) stop(notice, "\nNot downloaded; pass accept.license = TRUE to accept the terms.", call. = FALSE)
  } else {
    message(notice)
  }

  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(fileext = ".chain.gz")
  on.exit(unlink(tmp), add = TRUE)
  utils::download.file(ent$url, tmp, mode = "wb", quiet = TRUE)
  md5 <- unname(tools::md5sum(tmp))
  if (!identical(md5, ent$md5)) stop(sprintf("checksum mismatch for %s (md5 %s, expected %s)", ent$url, md5, ent$md5), call. = FALSE)
  if (!file.copy(tmp, dest, overwrite = TRUE)) stop("could not write ", dest, call. = FALSE)
  .gcanvas_note("gcanvas::chain.download", sprintf("saved %s", dest))
  invisible(abs_path(dest))
}
