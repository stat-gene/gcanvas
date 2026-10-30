# Small chain file whose blocks exercise each rule of liftover(). Target and
# query names are Ensembl-style; every block is 100 bp unless noted.
#   1: t[100,200) -> 1:q[200,300), gap (t 50, q 50), t[250,300) -> 1:q[350,400)
#   2: t[0,100) -> 2 on the "-" strand (qSize 1000): pos p -> 1001 - p
#   3: t[0,100) -> 4:q[500,600)                       (other chromosome)
#   5: t[0,100) -> 5:q[0,100) and 5:q[500,600)        (two positions)
#   6: t[0,100) -> 6:q[0,100) twice                   (identical hits)
#   7: t[0,100) -> 7:q[10,110) and HG99_PATCH:q[0,100)
#   8: t[0,100) -> GL000192.1:q[0,100)                (contig only)
#   X: t[0,100) -> X:q[1000,1100);  MT: t[0,100) -> MT:q[0,100)
#   HG1_PATCH: t[0,100) -> 9:q[0,100)                 (non-standard target)
toy_chain <- function() {
  path <- tempfile(fileext = ".chain")
  hdr <- function(tn, ts, te, qn, qs, qstr, q0, q1, id) {
    sprintf("chain 1 %s 1000 + %d %d %s %d %s %d %d %d", tn, ts, te, qn, qs, qstr, q0, q1, id)
  }
  writeLines(c(
    hdr("1", 100, 300, "1", 1000, "+", 200, 400, 1), "100 50 50", "50", "",
    hdr("2", 0, 100, "2", 1000, "-", 0, 100, 2), "100", "",
    hdr("3", 0, 100, "4", 1000, "+", 500, 600, 3), "100", "",
    hdr("5", 0, 100, "5", 1000, "+", 0, 100, 4), "100", "",
    hdr("5", 0, 100, "5", 1000, "+", 500, 600, 5), "100", "",
    hdr("6", 0, 100, "6", 1000, "+", 0, 100, 6), "100", "",
    hdr("6", 0, 100, "6", 1000, "+", 0, 100, 6), "100", "",
    hdr("7", 0, 100, "7", 1000, "+", 10, 110, 7), "100", "",
    hdr("7", 0, 100, "HG99_PATCH", 1000, "+", 0, 100, 8), "100", "",
    hdr("8", 0, 100, "GL000192.1", 1000, "+", 0, 100, 9), "100", "",
    hdr("X", 0, 100, "X", 2000, "+", 1000, 1100, 10), "100", "",
    hdr("MT", 0, 100, "MT", 1000, "+", 0, 100, 11), "100", "",
    hdr("HG1_PATCH", 0, 100, "9", 1000, "+", 0, 100, 12), "100", ""
  ), path)
  path
}

# UCSC-style names: chr1 -> chr1 (+10), chrM -> chrM on the "-" strand, and a
# non-standard target.
toy_chain_ucsc <- function() {
  path <- tempfile(fileext = ".chain")
  writeLines(c(
    "chain 1 chr1 1000 + 0 100 chr1 1000 + 10 110 1", "100", "",
    "chain 1 chrM 1000 + 0 100 chrM 1000 - 0 100 2", "100", "",
    "chain 1 chr1_gl000191_random 1000 + 0 100 chr2 1000 + 0 100 3", "100", ""
  ), path)
  path
}

lift_toy <- function(chr, pos, ...) {
  liftover(data.frame(CHR = chr, POS = pos, stringsAsFactors = FALSE), from = 38, to = 37,
           liftover.chain = toy_chain(), silent = TRUE, ...)
}

lift_bundled <- function(chr, pos, from = 38, to = 37, ...) {
  liftover(data.frame(CHR = chr, POS = pos, stringsAsFactors = FALSE), from = from, to = to, silent = TRUE, ...)
}

liftover_bin <- function() {
  bin <- .gcanvas_resolve_liftover_bin(NULL)
  if (is.na(bin)) "" else bin
}
