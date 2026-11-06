# Builds inst/extdata/GRCh37_to_GRCh38.ncbi.chain.gz and
# GRCh38_to_GRCh37.ncbi.chain.gz from the NCBI Remap alignments between the
# GRCh37 and GRCh38 primary assemblies. Only first-pass alignments are kept
# (reciprocity=3: best for both assemblies), so a base has one placement at
# most. Only alignments from a standard chromosome are written; the files have
# no MT, which liftover() takes from the Ensembl chain.
#
# In these GFFs column 1 is the second assembly of the file name and `Target`
# the first, so `Target` is the chain target. Gap ops: M both, I target only,
# D column 1 only.
#
#   Rscript data-raw/ncbi_remap_chains.R

suppressPackageStartupMessages(library(data.table))

remap <- "https://ftp.ncbi.nlm.nih.gov/pub/remap/Homo_sapiens/3.0"
genomes <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405"
reports <- c(
  "37" = file.path(genomes, "GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_assembly_report.txt"),
  "38" = file.path(genomes, "GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_report.txt")
)
gffs <- c(
  "37" = file.path(remap, "GCF_000001305.13_GRCh37_Primary_Assembly/GCF_000001305.16_GRCh38.p14_Primary_Assembly",
                   "GCF_000001305.13-GCF_000001305.16.gff"),
  "38" = file.path(remap, "GCF_000001305.16_GRCh38.p14_Primary_Assembly/GCF_000001305.13_GRCh37_Primary_Assembly",
                   "GCF_000001305.16-GCF_000001305.13.gff")
)

# Name and length per RefSeq accession: chromosomes get their plain name
# (1-22, X, Y), every other sequence keeps its accession.
read_report <- function(url) {
  x <- readLines(url, warn = FALSE)
  r <- fread(text = x[!startsWith(x, "#") | startsWith(x, "# Sequence-Name")], sep = "\t", header = TRUE)
  setnames(r, 1L, "Sequence-Name")
  chrom <- r$`Sequence-Role` == "assembled-molecule"
  data.table(acc = r$`RefSeq-Accn`, name = ifelse(chrom, r$`Sequence-Name`, r$`RefSeq-Accn`),
             size = as.numeric(r$`Sequence-Length`), chrom = chrom)
}

read_gff <- function(url) {
  x <- readLines(url, warn = FALSE)
  x <- x[!startsWith(x, "#")]
  f <- tstrsplit(x, "\t")
  tag <- function(k) {
    out <- rep(NA_character_, length(x))
    i <- grepl(paste0("(^|;)", k, "="), f[[9]])
    out[i] <- sub(paste0("^(.*;)?", k, "=([^;]*).*$"), "\\2", f[[9]][i])
    out
  }
  g <- data.table(q = f[[1]], q0 = as.numeric(f[[4]]) - 1, q1 = as.numeric(f[[5]]), q_strand = f[[7]],
                  target = tag("Target"), gap = tag("Gap"), reciprocity = tag("reciprocity"))
  g[, c("t", "t0", "t1", "t_strand") := tstrsplit(target, " ")]
  g[, `:=`(t0 = as.numeric(t0) - 1, t1 = as.numeric(t1), target = NULL)]
  g[]
}

# Gapless blocks of one alignment, ascending on the target. Positions are
# 0-based starts on the + strand of each sequence.
alignment_blocks <- function(a) {
  w <- if (is.na(a$gap)) paste0("M", a$q1 - a$q0) else strsplit(a$gap, " ")[[1]]
  op <- substr(w, 1, 1)
  n <- as.numeric(substring(w, 2))
  stopifnot(op %in% c("M", "I", "D"))
  q_step <- ifelse(op %in% c("M", "D"), n, 0)
  t_step <- ifelse(op %in% c("M", "I"), n, 0)
  stopifnot(sum(q_step) == a$q1 - a$q0, sum(t_step) == a$t1 - a$t0)
  m <- op == "M"
  q_off <- (cumsum(q_step) - q_step)[m]
  t_off <- (cumsum(t_step) - t_step)[m]
  size <- n[m]
  b <- data.table(
    size = size,
    t = if (a$t_strand == "-") a$t1 - t_off - size else a$t0 + t_off,
    q = if (a$q_strand == "-") a$q1 - q_off - size else a$q0 + q_off
  )
  b[order(t)]
}

chain_lines <- function(a, id, t_info, q_info) {
  b <- alignment_blocks(a)
  flip <- a$t_strand != a$q_strand
  q_size <- q_info$size
  q <- if (flip) q_size - (b$q + b$size) else b$q
  k <- nrow(b)
  dt <- b$t[-1] - (b$t[-k] + b$size[-k])
  dq <- q[-1] - (q[-k] + b$size[-k])
  stopifnot(dt >= 0, dq >= 0)
  c(
    sprintf("chain %.0f %s %.0f + %.0f %.0f %s %.0f %s %.0f %.0f %d", sum(b$size), t_info$name, t_info$size,
            b$t[1], b$t[k] + b$size[k], q_info$name, q_size, if (flip) "-" else "+", q[1], q[k] + b$size[k], id),
    if (k > 1L) sprintf("%.0f %.0f %.0f", b$size[-k], dt, dq),
    sprintf("%.0f", b$size[k]),
    ""
  )
}

write_chain <- function(from, to) {
  t_seq <- read_report(reports[[from]])
  q_seq <- read_report(reports[[to]])
  g <- read_gff(gffs[[from]])
  g <- g[reciprocity == "3" & t %in% t_seq[chrom == TRUE, acc]]
  stopifnot(g$q %in% q_seq$acc)
  o <- match(g$t, t_seq$acc)
  g <- g[order(o, t0)]
  lines <- unlist(lapply(seq_len(nrow(g)), function(i) {
    chain_lines(g[i], i, t_seq[acc == g$t[i]], q_seq[acc == g$q[i]])
  }))
  out <- sprintf("inst/extdata/GRCh%s_to_GRCh%s.ncbi.chain.gz", from, to)
  con <- gzfile(out, "wb", compression = 9L)
  on.exit(close(con))
  writeLines(lines, con)
  cat(sprintf("%s: %d alignments, %d lines\n", out, nrow(g), length(lines)))
}

write_chain("37", "38")
write_chain("38", "37")
