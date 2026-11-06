# Cross-check of liftover() against UCSC `liftOver -multiple`, row by row, in
# both directions (GRCh37 <-> GRCh38), for `auto` (NCBI chain, then Ensembl
# chain) and for the Ensembl chain alone. The expected values come from the raw
# `-multiple` output with the rules written out again here, independently of
# the package code:
#   hits at the same (chromosome, position, strand) merge; exactly one hit on
#   the input's own chromosome -> that position; none -> NA (other_chrom when a
#   hit is on another standard chromosome, else contig, or unmapped without
#   hits); two or more distinct ones -> NA (multiple).
#   `auto`: a row the NCBI chain gave no position on its chromosome takes the
#   Ensembl chain's position; an unmapped row takes the Ensembl status.
# Position sets: 100k random positions, every block boundary (+-1), a random
# base in every reverse-strand block, the middle of every target gap, a 20 kb
# grid, and (GRCh38 -> GRCh37) the regression cases and windows of
# tests/testthat/test-liftover.R. Also times both methods and writes
# tests/testthat/fixtures/liftover_b38_b37_windows.tsv.
#
#   conda create -p <dir> -c conda-forge -c bioconda ucsc-liftover
#   Rscript data-raw/liftover_validate.R <dir>/bin/liftOver

suppressPackageStartupMessages({
  library(data.table)
  pkgload::load_all(".", quiet = TRUE)
})
args <- commandArgs(trailingOnly = TRUE)
bin <- if (length(args)) args[1] else unname(Sys.which("liftOver"))
stopifnot(nzchar(bin), file.exists(bin))
set.seed(20261030)

std <- c(as.character(1:22), "X", "Y", "MT")
len <- list(
  "37" = c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022,
           141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753,
           81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560, 59373566, 16569),
  "38" = c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636,
           138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345,
           83257441, 80373285, 58617616, 64444167, 46709983, 50818468, 156040895, 57227415, 16569)
)
chains_of <- list(
  "37" = c(ncbi = system.file("extdata", "GRCh37_to_GRCh38.ncbi.chain.gz", package = "gcanvas"),
           ensembl = system.file("extdata", "GRCh37_to_GRCh38.chain.gz", package = "gcanvas")),
  "38" = c(ncbi = system.file("extdata", "GRCh38_to_GRCh37.ncbi.chain.gz", package = "gcanvas"),
           ensembl = system.file("extdata", "GRCh38_to_GRCh37.chain.gz", package = "gcanvas"))
)

# Target blocks of the chains on standard chromosomes, parsed chain by chain.
read_blocks <- function(path) {
  x <- readLines(gzfile(path))
  x <- x[nzchar(trimws(x))]
  h <- which(startsWith(x, "chain"))
  last <- c(h[-1] - 1L, length(x))
  res <- vector("list", length(h))
  for (k in seq_along(h)) {
    f <- strsplit(x[h[k]], "[ \t]+")[[1]]
    if (!(f[3] %in% std)) next
    m <- do.call(rbind, lapply(strsplit(x[(h[k] + 1L):last[k]], "[ \t]+"),
                               function(v) as.numeric(c(v, 0, 0)[1:3])))
    t0 <- as.numeric(f[6]) + c(0, cumsum(m[, 1] + m[, 2]))[seq_len(nrow(m))]
    res[[k]] <- data.table(chr = f[3], tStart = t0, tEnd = t0 + m[, 1], dt = m[, 2], qStrand = f[10])
  }
  rbindlist(res)
}

grid <- rbindlist(lapply(1:24, function(i) data.table(chr = std[i], pos = seq(5000, len[["37"]][i], by = 20000))))
cases <- data.table(chr = "1", pos = c(1005000, 143185000, 3925000, 50005000, 60005000))
# 500 positions per window at the middle of each 1 kb bin (142,900,500, 142,901,500, ...)
windows <- data.table(
  window = rep(c("chr1:142.9-143.4Mb", "chr10:46.9-47.4Mb", "chr17:21.0-21.5Mb"), each = 500L),
  chr = rep(c("1", "10", "17"), each = 500L),
  pos = c(142.9e6, 46.9e6, 21.0e6)[rep(1:3, each = 500L)] + 500 + 1000 * rep(0:499, 3L)
)

position_sets <- function(from) {
  b <- rbindlist(lapply(chains_of[[from]], read_blocks))
  L <- len[[from]]
  i <- sample.int(25L, 1e5, replace = TRUE, prob = L)
  out <- list(
    random = data.table(chr = std[i], pos = ceiling(stats::runif(1e5) * L[i])),
    block_edges = unique(rbind(b[, .(chr, pos = tStart)], b[, .(chr, pos = tStart + 1)],
                               b[, .(chr, pos = tEnd)], b[, .(chr, pos = tEnd + 1)])[pos >= 1]),
    reverse_blocks = b[qStrand == "-", .(chr, pos = tStart + 1 + floor(stats::runif(.N) * (tEnd - tStart)))],
    gap_middles = b[dt > 0, .(chr, pos = tEnd + 1 + floor(dt / 2))],
    grid_20kb = grid
  )
  if (from == "38") {
    out$cases <- cases
    out$windows <- windows[, .(chr, pos)]
  }
  out
}

# Raw `liftOver -multiple` hits for rows (chr, pos).
binary_hits <- function(chr, pos, chain) {
  f <- tempfile(c("in", "out", "un"), fileext = ".bed")
  on.exit(unlink(f))
  pos <- as.integer(pos)
  fwrite(data.table(chr, pos - 1L, pos, seq_along(pos), 0L, "+"), f[1], sep = "\t", col.names = FALSE)
  rc <- system2(bin, c("-multiple", f[1], chain, f[2], f[3]), stdout = FALSE, stderr = FALSE)
  stopifnot(rc == 0L)
  if (!file.size(f[2])) return(data.table(row = integer(), chr = character(), pos = integer(), strand = character()))
  out <- fread(f[2], header = FALSE, sep = "\t", colClasses = "character")
  data.table(row = as.integer(out$V4), chr = out$V1, pos = as.integer(out$V3), strand = out$V6)
}

expected <- function(chr, pos, chain) {
  n <- length(pos)
  h <- unique(binary_hits(chr, pos, chain))
  conflicts <- h[, uniqueN(strand), by = .(row, chr, pos)][V1 > 1L, .N]
  h[, same := chr == ..chr[row]]
  per <- h[, .(n_same = sum(same), n_other = sum(!same & chr %in% std),
               pos = if (sum(same) == 1L) pos[same] else NA_integer_,
               strand = if (sum(same) == 1L) strand[same] else NA_character_), by = row]
  e <- data.table(pos = rep(NA_integer_, n), strand = NA_character_, status = "unmapped")
  e[per$row, `:=`(pos = per$pos, strand = per$strand,
                  status = with(per, fifelse(n_same == 1L, "ok", fifelse(n_same > 1L, "multiple",
                                                                        fifelse(n_other > 0L, "other_chrom", "contig")))))]
  list(e = e, conflicts = conflicts, n_multi_hit = h[, .N, by = row][N > 1L, .N])
}

# `auto` from the expected rows of the NCBI (a) and Ensembl (b) chains.
combine <- function(a, b) {
  open <- a$status %in% c("unmapped", "other_chrom", "contig")
  use <- (open & b$status == "ok") | (a$status == "unmapped" & b$status != "unmapped")
  out <- copy(a)
  out[use, `:=`(pos = b$pos[use], strand = b$strand[use], status = b$status[use])]
  out[]
}

neq <- function(a, b) !((is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b))

compare <- function(got, e, to) {
  p <- got[[paste0("POS_b", to)]]
  cc <- got[[paste0("CHR_b", to)]]
  s <- got[[paste0("STRAND_b", to)]]
  bad <- neq(p, e$pos) | neq(s, e$strand) | neq(got$lift_status, e$status) |
    neq(cc, ifelse(e$status == "ok", got$CHR, NA))
  c(rows = nrow(got), mismatch = sum(bad))
}

report <- list()
for (from in c("37", "38")) {
  to <- if (from == "37") "38" else "37"
  ens <- chains_of[[from]][["ensembl"]]
  sets <- position_sets(from)
  for (nm in names(sets)) {
    d <- sets[[nm]]
    d[, pos := as.integer(pos)]
    ex_n <- expected(d$chr, d$pos, chains_of[[from]][["ncbi"]])
    ex_e <- expected(d$chr, d$pos, ens)
    e <- combine(ex_n$e, ex_e$e)
    in_df <- data.frame(CHR = d$chr, POS = d$pos)
    g_chain <- liftover(in_df, from = from, to = to, method = "chain", silent = TRUE)
    g_bin <- liftover(in_df, from = from, to = to, method = "binary", liftover = bin, silent = TRUE)
    g_ens <- liftover(in_df, from = from, to = to, method = "chain", liftover.chain = ens, silent = TRUE)
    st <- table(factor(e$status, levels = c("ok", "unmapped", "other_chrom", "contig", "multiple")))
    report[[length(report) + 1L]] <- data.table(
      direction = sprintf("%s->%s", from, to), set = nm, n = nrow(d),
      chain_vs_expected = compare(g_chain, e, to)[["mismatch"]],
      binary_vs_expected = compare(g_bin, e, to)[["mismatch"]],
      ensembl_vs_expected = compare(g_ens, ex_e$e, to)[["mismatch"]],
      rows_out_equal_in = nrow(g_chain) == nrow(d) && nrow(g_bin) == nrow(d),
      ok = st[["ok"]], unmapped = st[["unmapped"]], other_chrom = st[["other_chrom"]],
      contig = st[["contig"]], multiple = st[["multiple"]],
      ncbi_multi_hit_rows = ex_n$n_multi_hit, strand_conflicts = ex_n$conflicts + ex_e$conflicts,
      minus_strand_ok = sum(g_chain[[paste0("STRAND_b", to)]] == "-", na.rm = TRUE)
    )
    if (from == "38" && nm == "windows") {
      fx <- cbind(windows[, .(window, chr, pos)], pos_b37 = g_ens$POS_b37, strand_b37 = g_ens$STRAND_b37,
                  status = g_ens$lift_status)
      dir.create("tests/testthat/fixtures", recursive = TRUE, showWarnings = FALSE)
      fwrite(fx, "tests/testthat/fixtures/liftover_b38_b37_windows.tsv", sep = "\t", na = "NA", quote = FALSE)
    }
  }
}
report <- rbindlist(report)
print(report, nrows = 100)

# The 20 kb grid in liftOver's default (single) mode on the Ensembl chain, as in
# the original measurement.
single_mode <- function(from) {
  f <- tempfile(c("in", "out", "un"), fileext = ".bed")
  on.exit(unlink(f))
  fwrite(grid[, .(chr, as.integer(pos) - 1L, as.integer(pos), .I)], f[1], sep = "\t", col.names = FALSE)
  system2(bin, c(f[1], chains_of[[from]][["ensembl"]], f[2], f[3]), stdout = FALSE, stderr = FALSE)
  out <- fread(f[2], header = FALSE, sep = "\t", colClasses = "character")
  un <- readLines(f[3])
  data.table(direction = from, lifted = nrow(out), lifted_to_contig = sum(!(out$V1 %in% std)),
             lifted_to_other_chrom = sum(out$V1 %in% std & out$V1 != grid$chr[as.integer(out$V4)]),
             deleted = sum(un == "#Deleted in new"), duplicated = sum(un == "#Duplicated in new"))
}
cat("\n20 kb grid, liftOver single mode:\n")
print(rbind(single_mode("37"), single_mode("38")))

w <- fread("tests/testthat/fixtures/liftover_b38_b37_windows.tsv")
cat("\nGRCh38 -> GRCh37 windows on the Ensembl chain, rows lifted to the same chromosome:\n")
print(w[, .(n = .N, ok = sum(status == "ok")), by = window])

# Timing: chain (cold = parse + lookup, warm = cached chain) vs binary.
cat("\nTiming (seconds):\n")
for (n in c(1e5, 1e6)) {
  for (from in c("37", "38")) {
    to <- if (from == "37") "38" else "37"
    L <- len[[from]]
    i <- sample.int(25L, n, replace = TRUE, prob = L)
    d <- data.frame(CHR = std[i], POS = ceiling(stats::runif(n) * L[i]))
    rm(list = ls(envir = .gcanvas_chain_cache), envir = .gcanvas_chain_cache)
    t_cold <- system.time(liftover(d, from = from, to = to, silent = TRUE))[["elapsed"]]
    t_warm <- system.time(liftover(d, from = from, to = to, silent = TRUE))[["elapsed"]]
    t_bin <- system.time(liftover(d, from = from, to = to, method = "binary", liftover = bin, silent = TRUE))[["elapsed"]]
    cat(sprintf("  n=%-9s %s->%s  chain cold %.2f  chain warm %.2f  binary %.2f\n",
                format(n, big.mark = ",", scientific = FALSE), from, to, t_cold, t_warm, t_bin))
  }
}
