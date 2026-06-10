# Coordinate liftover via the UCSC `liftOver` binary and chain files.

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

# Resolve a chain file path. Honours an explicit `liftover.chain`; otherwise
# searches (in order) `liftover.dir/chain`, `liftover.dir`, and the package's
# bundled `inst/extdata` chains for a name matching the requested build pair.
.gcanvas_resolve_chain <- function(from, to, liftover.dir, liftover.chain) {
  chain0 <- as.character(liftover.chain)[1]
  auto <- is.null(liftover.chain) || length(liftover.chain) == 0L ||
    is.na(chain0) || !nzchar(chain0) || toupper(chain0) == "AUTO"
  if (!auto) return(chain0)

  from_key <- .gcanvas_build_key(from)
  to_key <- .gcanvas_build_key(to)
  if (is.na(from_key) || is.na(to_key)) return(NA_character_)
  names <- .gcanvas_chain_candidates(from_key, to_key)
  # Legacy generic naming honoured in user-supplied dirs.
  legacy <- sprintf("%s_to_%s.chain.gz", from_key, to_key)
  names <- unique(c(names, legacy))
  if (!length(names)) return(NA_character_)

  dirs <- character()
  if (!is.null(liftover.dir) && length(liftover.dir) && !is.na(liftover.dir[1]) && nzchar(liftover.dir[1])) {
    dirs <- c(file.path(liftover.dir, "chain"), liftover.dir)
  }
  bundle_dir <- system.file("extdata", package = "gcanvas")
  if (nzchar(bundle_dir)) dirs <- c(dirs, bundle_dir)

  for (d in dirs) {
    for (nm in names) {
      p <- file.path(d, nm)
      if (file.exists(p)) return(p)
    }
  }
  NA_character_
}

# Locate the UCSC liftOver binary: prefer `liftover.dir/liftOver`, then PATH.
.gcanvas_resolve_liftover_bin <- function(liftover.dir) {
  if (!is.null(liftover.dir) && length(liftover.dir) && !is.na(liftover.dir[1]) && nzchar(liftover.dir[1])) {
    p <- file.path(liftover.dir, "liftOver")
    if (file.exists(p)) return(p)
  }
  w <- Sys.which("liftOver")
  if (nzchar(w)) return(unname(w))
  NA_character_
}

# Detect the chromosome naming convention of a chain file from its first
# `chain` header line (target/tName field). Returns "ucsc" (chr-prefixed) or
# "plain" (Ensembl-style, no prefix).
.gcanvas_chain_seqstyle <- function(chain_path) {
  con <- tryCatch(
    if (grepl("\\.gz$", chain_path, ignore.case = TRUE)) gzfile(chain_path, "rt") else file(chain_path, "rt"),
    error = function(e) NULL
  )
  if (is.null(con)) return("plain")
  on.exit(close(con), add = TRUE)
  style <- "plain"
  repeat {
    ln <- tryCatch(readLines(con, n = 1L, warn = FALSE), error = function(e) character())
    if (!length(ln)) break
    if (startsWith(ln, "chain")) {
      f <- strsplit(trimws(ln), "\\s+")[[1]]
      if (length(f) >= 3L) {
        style <- if (grepl("^chr", f[3], ignore.case = TRUE)) "ucsc" else "plain"
      }
      break
    }
  }
  style
}

#' Lift coordinates between genome builds
#'
#' Wraps the UCSC `liftOver` binary to translate variant coordinates from one
#' build to another. `from`/`to` accept `38`/`hg38`/`GRCh38`, `37`/`19`/`hg19`/
#' `GRCh37`, and `36`/`18`/`hg18`/`NCBI36`.
#'
#' When `liftover.chain = "auto"` (the default), the chain file is resolved
#' automatically: first from `liftover.dir` (a `chain/` subdir or the dir
#' itself), then from the chain files bundled with the package
#' (`GRCh37<->GRCh38` from Ensembl and `hg18->hg19` from UCSC). The chromosome
#' naming convention of the chosen chain (UCSC `chr1` vs Ensembl `1`) is
#' detected from its header and the input coordinates are reformatted to match,
#' so mixed-convention inputs lift correctly. The `liftOver` binary itself is
#' not bundled; it is located via `liftover.dir` or the system `PATH`.
#'
#' @param df A `data.frame`/`data.table` of variants.
#' @param from,to Source and destination genome builds. Accepts numeric or
#'   string aliases (see Description). Supported pairs: 37<->38 and 18->19.
#' @param liftover.dir Optional directory containing the `liftOver` binary and
#'   chain files. Chains are looked up in `liftover.dir/chain` and `liftover.dir`.
#' @param liftover.chain Explicit chain file path, or `"auto"` (the default;
#'   `NULL` is treated the same) to resolve from `liftover.dir` / the bundled
#'   chains.
#' @param SNP,CHR,POS Column names in `df` for the SNP id, chromosome, and position.
#' @param rm.tmp Logical. Remove temporary BED files when done.
#' @param silent Logical. Suppress progress notes.
#'
#' @return A `data.table` with updated `CHR`/`POS` and (where applicable)
#'   columns indicating which variants failed to lift.
#' @export
liftover <- function(df, from = 37, to = 38,
                     liftover.dir = NULL,
                     liftover.chain = "auto",
                     SNP = "SNP", CHR = "CHR", POS = "POS",
                     rm.tmp = TRUE, silent = FALSE) {
  require_pkg("data.table")

  liftover.chain <- .gcanvas_resolve_chain(from, to, liftover.dir, liftover.chain)
  if (is.na(liftover.chain) || !nzchar(liftover.chain) || !file.exists(liftover.chain)) {
    if (!silent) {
      message("No chain file found for build ", from, " -> ", to,
              " (looked in liftover.dir and the bundled chains).")
    }
    return(df)
  }
  chain_style <- .gcanvas_chain_seqstyle(liftover.chain)

  lift_bin <- .gcanvas_resolve_liftover_bin(liftover.dir)
  if (is.na(lift_bin)) {
    if (!silent) message("No liftOver binary found (set liftover.dir or put `liftOver` on PATH).")
    return(df)
  }

  to_label <- .gcanvas_build_num(to)
  if (is.na(to_label)) to_label <- as.character(to)[1]

  is_dt <- data.table::is.data.table(df)
  rn0 <- if (!is_dt) rownames(df) else NULL
  dt <- if (is_dt) data.table::copy(df) else data.table::as.data.table(df)
  snp_use <- .gcanvas_resolve_colname(names(dt), SNP, aliases = character(), required = TRUE, arg_label = "SNP")
  chr_use <- .gcanvas_resolve_colname(names(dt), CHR, aliases = character(), required = TRUE, arg_label = "CHR")
  pos_use <- .gcanvas_resolve_colname(names(dt), POS, aliases = character(), required = TRUE, arg_label = "POS")
  finalize <- function(xdt) {
    if (is_dt) return(xdt)
    out <- as.data.frame(xdt, stringsAsFactors = FALSE)
    if (!is.null(rn0) && length(rn0) == nrow(out)) rownames(out) <- rn0
    out
  }

  fmt_chr <- if (identical(chain_style, "ucsc")) {
    function(z) chr_to_ucsc1(z)
  } else {
    function(z) {
      v <- normalize.chrom(z)[1]
      if (is.na(v) || !nzchar(v)) NA_character_ else v
    }
  }

  bed <- dt[, c(pos_use, snp_use), with = FALSE]
  bed[, CHR := vapply(dt[[chr_use]], fmt_chr, character(1))]
  bed[, POSm1 := suppressWarnings(as.numeric(get(pos_use))) - 1]
  data.table::setcolorder(bed, c("CHR", "POSm1", pos_use, snp_use))
  bed <- bed[is.finite(POSm1) & !is.na(CHR)]
  if (!nrow(bed)) return(finalize(dt))

  bed[, c("POSm1", pos_use) := lapply(.SD, function(x) format(x, scientific = FALSE)), .SDcols = c("POSm1", pos_use)]
  bed <- unique(bed, by = c("CHR", "POSm1", pos_use, snp_use))

  randn <- round(stats::runif(1) * 1e9)
  bedf <- sprintf("tmp.%s.liftover.bed", randn)
  liftedf <- sprintf("tmp.%s.lifted", randn)
  unliftedf <- sprintf("tmp.%s.unlifted", randn)
  write.table(as.data.frame(bed), file = bedf, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

  cmd <- sprintf('bash -c "%s %s %s %s %s"',
                 shQuote(lift_bin), shQuote(bedf), shQuote(liftover.chain), shQuote(liftedf), shQuote(unliftedf))
  try(system(cmd, intern = TRUE), silent = TRUE)

  lifted <- tryCatch(
    data.table::fread(liftedf, header = FALSE, colClasses = "character", data.table = TRUE, verbose = FALSE),
    error = function(e) data.table::data.table()
  )

  if (rm.tmp) {
    suppressWarnings(try(unlink(c(bedf, liftedf, unliftedf)), silent = TRUE))
  }

  if (!nrow(lifted)) return(finalize(dt))

  lifted <- lifted[, .(tmp_snp = V4, tmp_pos = suppressWarnings(as.numeric(V3)))]
  data.table::setnames(lifted, c("tmp_snp", "tmp_pos"), c(snp_use, paste0(pos_use, "_b", to_label)))

  dt[, ..ord := .I]
  dt[lifted, on = setNames(snp_use, snp_use), (paste0(pos_use, "_b", to_label)) := get(paste0("i.", paste0(pos_use, "_b", to_label)))]
  data.table::setorder(dt, ..ord)
  dt[, ..ord := NULL]

  finalize(dt)
}

liftover_positions <- function(chrom, pos, from, to, liftover.dir, liftover.chain = "auto") {
  pos <- .gcanvas_as_num2(pos)
  if (!length(pos)) return(pos)
  key <- paste0("K", seq_along(pos))
  tmp <- data.frame(SNP = key, CHR = rep(normalize.chrom(chrom)[1], length(pos)), POS = pos, stringsAsFactors = FALSE)
  out <- liftover(tmp, from = from, to = to, liftover.dir = liftover.dir,
                  liftover.chain = liftover.chain, SNP = "SNP", CHR = "CHR", POS = "POS",
                  silent = TRUE)
  to_label <- .gcanvas_build_num(to)
  if (is.na(to_label)) to_label <- as.character(to)[1]
  pos2 <- out[[paste0("POS_b", to_label)]]
  names(pos2) <- out$SNP
  unname(pos2[key])
}

