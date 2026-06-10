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
# resolves from the package's bundled `inst/extdata` chains by build pair.
.gcanvas_resolve_chain <- function(from, to, liftover.chain) {
  chain0 <- as.character(liftover.chain)[1]
  auto <- is.null(liftover.chain) || length(liftover.chain) == 0L ||
    is.na(chain0) || !nzchar(chain0) || toupper(chain0) == "AUTO"
  if (!auto) return(chain0)

  from_key <- .gcanvas_build_key(from)
  to_key <- .gcanvas_build_key(to)
  if (is.na(from_key) || is.na(to_key)) return(NA_character_)
  names <- .gcanvas_chain_candidates(from_key, to_key)
  if (!length(names)) return(NA_character_)

  bundle_dir <- system.file("extdata", package = "gcanvas")
  if (!nzchar(bundle_dir)) return(NA_character_)
  for (nm in names) {
    p <- file.path(bundle_dir, nm)
    if (file.exists(p)) return(p)
  }
  NA_character_
}

# Locate the UCSC liftOver binary. Uses the user-supplied `liftover` value when
# given (a file path or a command name), otherwise defaults to looking up
# `liftOver` on the environment `PATH`. The search mirrors the PLINK resolver:
# direct file -> `Sys.which` -> common bin dirs -> login-shell `command -v`
# (the last step finds conda-env installs that R's own PATH may not expose).
.gcanvas_resolve_liftover_bin <- function(liftover = NULL) {
  cmd <- if (!is.null(liftover) && length(liftover) && !is.na(liftover[1]) && nzchar(liftover[1])) {
    as.character(liftover[1])
  } else {
    "liftOver"
  }
  if (file.exists(cmd)) return(abs_path(cmd))
  w <- Sys.which(cmd)
  if (nzchar(w)) return(unname(w))
  cand_dir <- c(path.expand("~/bin"), "/opt/homebrew/bin", "/usr/local/bin", "/usr/bin", "/bin")
  cand <- file.path(cand_dir, basename(cmd))
  hit <- cand[file.exists(cand)]
  if (length(hit)) return(abs_path(hit[1]))
  shs <- unique(c(Sys.getenv("SHELL", unset = "/bin/sh"), "/bin/zsh", "/bin/bash", "/bin/sh"))
  shs <- shs[file.exists(shs)]
  for (sh in shs) {
    out <- tryCatch(
      suppressWarnings(system2(sh, c("-lc", sprintf("command -v %s 2>/dev/null", shQuote(cmd))), stdout = TRUE, stderr = TRUE)),
      error = function(e) character(0)
    )
    out <- trimws(as.character(out))
    out <- out[nzchar(out)]
    if (!length(out)) next
    last <- out[length(out)]
    if (file.exists(last)) return(abs_path(last))
    p2 <- Sys.which(last)
    if (nzchar(p2)) return(unname(p2))
  }
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
#' automatically from the chain files bundled with the package
#' (`GRCh37<->GRCh38` from Ensembl and `hg18->hg19` from UCSC). The chromosome
#' naming convention of the chosen chain (UCSC `chr1` vs Ensembl `1`) is
#' detected from its header and the input coordinates are reformatted to match,
#' so mixed-convention inputs lift correctly. For build pairs that are not
#' bundled, pass an explicit chain file via `liftover.chain`.
#'
#' The `liftOver` binary itself is not bundled. When the `liftover` argument is
#' `NULL`, the executable is resolved automatically by searching the system
#' `PATH` (and common install locations, including conda environments), so a
#' `liftOver` reachable from your shell is used without any extra configuration.
#' The recommended way to install it is via conda:
#' \preformatted{
#' conda install bioconda::ucsc-liftover
#' }
#' Alternatively, download the UCSC binary directly (Linux x86_64 example):
#' \preformatted{
#' wget https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
#' chmod +x ./liftOver
#' }
#' Binaries for other platforms are under
#' <https://hgdownload.soe.ucsc.edu/admin/exe/>.
#'
#' @param df A `data.frame`/`data.table` of variants.
#' @param from,to Source and destination genome builds. Accepts numeric or
#'   string aliases (see Description). Supported pairs: 37<->38 and 18->19.
#' @param liftover Optional path to the `liftOver` executable (a file path or a
#'   command name on `PATH`). When `NULL` (the default), `liftOver` is resolved
#'   automatically from the system `PATH` and common install locations
#'   (including conda environments).
#' @param liftover.chain Explicit chain file path, or `"auto"` (the default;
#'   `NULL` is treated the same) to resolve from the bundled chains.
#' @param snp.col,chrom.col,pos.col Column names in `df` for the SNP id,
#'   chromosome, and position.
#' @param rm.tmp Logical. Remove temporary BED files when done.
#' @param silent Logical. Suppress progress notes.
#'
#' @return A `data.table` with updated `CHR`/`POS` and (where applicable)
#'   columns indicating which variants failed to lift.
#' @export
liftover <- function(df, from = 37, to = 38,
                     liftover = NULL,
                     liftover.chain = "auto",
                     snp.col = "SNP", chrom.col = "CHR", pos.col = "POS",
                     rm.tmp = TRUE, silent = FALSE) {
  require_pkg("data.table")

  liftover.chain <- .gcanvas_resolve_chain(from, to, liftover.chain)
  if (is.na(liftover.chain) || !nzchar(liftover.chain) || !file.exists(liftover.chain)) {
    if (!silent) {
      message("No chain file found for build ", from, " -> ", to,
              " (no bundled chain; pass `liftover.chain` for this build pair).")
    }
    return(df)
  }
  chain_style <- .gcanvas_chain_seqstyle(liftover.chain)

  lift_bin <- .gcanvas_resolve_liftover_bin(liftover)
  if (is.na(lift_bin)) {
    if (!silent) message("No liftOver binary found. Set `liftover` to the liftOver executable path or put `liftOver` on PATH (download: https://hgdownload.soe.ucsc.edu/admin/exe/).")
    return(df)
  }

  to_label <- .gcanvas_build_num(to)
  if (is.na(to_label)) to_label <- as.character(to)[1]

  is_dt <- data.table::is.data.table(df)
  rn0 <- if (!is_dt) rownames(df) else NULL
  dt <- if (is_dt) data.table::copy(df) else data.table::as.data.table(df)
  snp_use <- .gcanvas_resolve_colname(names(dt), snp.col, aliases = character(), required = TRUE, arg_label = "snp.col")
  chr_use <- .gcanvas_resolve_colname(names(dt), chrom.col, aliases = character(), required = TRUE, arg_label = "chrom.col")
  pos_use <- .gcanvas_resolve_colname(names(dt), pos.col, aliases = character(), required = TRUE, arg_label = "pos.col")
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

  lifted <- lifted[, .(tmp_snp = V4, tmp_pos = suppressWarnings(as.integer(as.numeric(V3))))]
  data.table::setnames(lifted, c("tmp_snp", "tmp_pos"), c(snp_use, paste0(pos_use, "_b", to_label)))

  dt[, ..ord := .I]
  dt[lifted, on = setNames(snp_use, snp_use), (paste0(pos_use, "_b", to_label)) := get(paste0("i.", paste0(pos_use, "_b", to_label)))]
  data.table::setorder(dt, ..ord)
  dt[, ..ord := NULL]

  finalize(dt)
}

liftover_positions <- function(chrom, pos, from, to, liftover = NULL, liftover.chain = "auto") {
  pos <- .gcanvas_as_num2(pos)
  if (!length(pos)) return(pos)
  key <- paste0("K", seq_along(pos))
  tmp <- data.frame(SNP = key, CHR = rep(normalize.chrom(chrom)[1], length(pos)), POS = pos, stringsAsFactors = FALSE)
  out <- liftover(tmp, from = from, to = to, liftover = liftover,
                  liftover.chain = liftover.chain, snp.col = "SNP", chrom.col = "CHR", pos.col = "POS",
                  silent = TRUE)
  to_label <- .gcanvas_build_num(to)
  if (is.na(to_label)) to_label <- as.character(to)[1]
  pos2 <- out[[paste0("POS_b", to_label)]]
  if (is.null(pos2)) return(rep(NA_real_, length(pos)))
  names(pos2) <- out$SNP
  unname(pos2[key])
}

