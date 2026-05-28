# Coordinate liftover via the UCSC `liftOver` binary and chain files.

#' Lift coordinates between genome builds
#'
#' Wraps the UCSC `liftOver` binary to translate variant coordinates from one
#' build (e.g. `37`) to another (e.g. `38`). The chain file is auto-resolved
#' from `liftover.dir` when not provided.
#'
#' @param df A `data.frame`/`data.table` of variants.
#' @param from,to Source and destination genome builds (`37` or `38`).
#' @param liftover.dir Optional directory containing `liftOver` and chain files.
#' @param liftover.chain Optional explicit chain file path.
#' @param SNP,CHR,POS Column names in `df` for the SNP id, chromosome, and position.
#' @param rm.tmp Logical. Remove temporary BED files when done.
#' @param silent Logical. Suppress progress notes.
#'
#' @return A `data.table` with updated `CHR`/`POS` and (where applicable)
#'   columns indicating which variants failed to lift.
#' @export
liftover <- function(df, from = 37, to = 38,
                     liftover.dir = NULL,
                     liftover.chain = NULL,
                     SNP = "SNP", CHR = "CHR", POS = "POS",
                     rm.tmp = TRUE, silent = FALSE) {
  require_pkg("data.table")

  chain0 <- as.character(liftover.chain)[1]
  chain_is_auto <- is.null(liftover.chain) || length(liftover.chain) == 0L ||
    is.na(chain0) || !nzchar(chain0) || toupper(chain0) == "AUTO"
  if (chain_is_auto) {
    if (from == 19) from <- 37
    if (to == 19) to <- 37
    liftover.chain <- file.path(liftover.dir, "chain", sprintf("GRCh%s_to_GRCh%s.chain.gz", from, to))
  }
  if (!file.exists(liftover.chain)) {
    if (!silent) message("No chain file: ", liftover.chain)
    return(df)
  }

  lift_bin <- file.path(liftover.dir, "liftOver")
  if (!file.exists(lift_bin)) {
    if (!silent) message("No liftOver binary: ", lift_bin)
    return(df)
  }

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

  bed <- dt[, c(pos_use, snp_use), with = FALSE]
  bed[, CHR := vapply(dt[[chr_use]], chr_to_ucsc1, character(1))]
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
  data.table::setnames(lifted, c("tmp_snp", "tmp_pos"), c(snp_use, paste0(pos_use, "_b", to)))

  dt[, ..ord := .I]
  dt[lifted, on = setNames(snp_use, snp_use), (paste0(pos_use, "_b", to)) := get(paste0("i.", paste0(pos_use, "_b", to)))]
  data.table::setorder(dt, ..ord)
  dt[, ..ord := NULL]

  finalize(dt)
}

liftover_positions <- function(chrom, pos, from, to, liftover.dir, liftover.chain = NULL) {
  pos <- .gcanvas_as_num2(pos)
  if (!length(pos)) return(pos)
  key <- paste0("K", seq_along(pos))
  tmp <- data.frame(SNP = key, CHR = rep(normalize.chrom(chrom)[1], length(pos)), POS = pos, stringsAsFactors = FALSE)
  out <- liftover(tmp, from = from, to = to, liftover.dir = liftover.dir,
                  liftover.chain = liftover.chain, SNP = "SNP", CHR = "CHR", POS = "POS",
                  silent = TRUE)
  pos2 <- out[[paste0("POS_b", to)]]
  names(pos2) <- out$SNP
  unname(pos2[key])
}

