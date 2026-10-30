# Expected values were measured with UCSC `liftOver -multiple` on the bundled
# Ensembl chains (data-raw/liftover_validate.R, which also writes the window
# fixture).

test_that("GRCh38 -> GRCh37 regression cases", {
  out <- lift_bundled("1", c(1005000, 143185000, 3925000, 50005000, 60005000))
  # 143,185,000 hits only HG1287_PATCH:143,871,415; 3,925,000 also hits HG1472_PATCH
  expect_identical(out$POS_b37, c(940380L, NA, 3841564L, 50470672L, 60470672L))
  expect_identical(out$lift_status, c("ok", "contig", "ok", "ok", "ok"))
  expect_identical(out$CHR_b37, c("1", NA, "1", "1", "1"))
  expect_identical(out$STRAND_b37, c("+", NA, "+", "+", "+"))
})

test_that("GRCh38 -> GRCh37 windows: same-chromosome counts and positions", {
  fx <- utils::read.delim(test_path("fixtures", "liftover_b38_b37_windows.tsv"),
                          colClasses = c(window = "character", chr = "character", pos = "integer",
                                         pos_b37 = "integer", strand_b37 = "character", status = "character"))
  out <- lift_bundled(fx$chr, fx$pos)
  expect_identical(out$POS_b37, fx$pos_b37)
  expect_identical(out$STRAND_b37, fx$strand_b37)
  expect_identical(out$lift_status, fx$status)
  n_ok <- tapply(out$lift_status == "ok", fx$window, sum)
  expect_identical(as.integer(n_ok[c("chr1:142.9-143.4Mb", "chr10:46.9-47.4Mb", "chr17:21.0-21.5Mb")]),
                   c(1L, 175L, 500L))
})

test_that("GRCh37 -> GRCh38 basics", {
  out <- lift_bundled(c("1", "1", "chrX", "23"), c(1005000, 3925000, 2500000, 2500000), from = 37, to = 38)
  expect_identical(out$POS_b38, c(1069620L, NA, 2581959L, 2581959L))
  expect_identical(out$lift_status, c("ok", "unmapped", "ok", "ok"))
  expect_identical(out$CHR_b38, c("1", NA, "chrX", "23"))
})

test_that("rule 1: one row per input row, attached by row and never by SNP id", {
  df <- data.frame(SNP = c("rs1", "rs1", NA, "rs1"), CHR = 1, POS = c(50005000, 60005000, 1005000, 50005000))
  out <- liftover(df, from = 38, to = 37, silent = TRUE)
  expect_identical(nrow(out), 4L)
  expect_identical(out$SNP, df$SNP)
  expect_identical(out$POS, df$POS)
  expect_identical(out$POS_b37, c(50470672L, 60470672L, 940380L, 50470672L))
  expect_identical(liftover(df[c("CHR", "POS")], from = 38, to = 37, silent = TRUE)$POS_b37, out$POS_b37)
  o <- c(4L, 2L, 3L, 1L)
  expect_identical(liftover(df[o, ], from = 38, to = 37, silent = TRUE)$POS_b37, out$POS_b37[o])
})

test_that("input class, row names and the input itself are preserved", {
  df <- data.frame(CHR = c(1, 1), POS = c(1005000, 143185000), row.names = c("a", "b"))
  out <- liftover(df, from = 38, to = 37, silent = TRUE)
  expect_s3_class(out, "data.frame")
  expect_false(data.table::is.data.table(out))
  expect_identical(rownames(out), c("a", "b"))
  expect_identical(class(out$CHR_b37), "numeric")

  dt <- data.table::data.table(SNP = "x", CHR = "1", POS = 1005000)
  out_dt <- liftover(dt, from = 38, to = 37, silent = TRUE)
  expect_true(data.table::is.data.table(out_dt))
  expect_false("lift_status" %in% names(dt))
  # lifting again replaces the columns instead of adding new ones
  again <- liftover(out_dt, from = 38, to = 37, silent = TRUE)
  expect_identical(names(again), names(out_dt))
})

test_that("custom column names give matching output names", {
  df <- data.frame(chrom = "1", bp = 1005000)
  out <- liftover(df, from = 38, to = 37, chrom.col = "chrom", pos.col = "bp", silent = TRUE)
  expect_identical(names(out), c("chrom", "bp", "bp_b37", "chrom_b37", "STRAND_b37", "lift_status"))
  expect_identical(out$bp_b37, 940380L)
})

test_that("block boundaries, gaps and the reverse strand (toy chain)", {
  out <- lift_toy("1", c(100, 101, 200, 201, 250, 251, 300, 301))
  expect_identical(out$POS_b37, c(NA, 201L, 300L, NA, NA, 351L, 400L, NA))
  expect_identical(out$lift_status, c("unmapped", "ok", "ok", "unmapped", "unmapped", "ok", "ok", "unmapped"))

  rev <- lift_toy("2", c(1, 50, 100))
  expect_identical(rev$POS_b37, c(1000L, 951L, 901L))
  expect_identical(rev$STRAND_b37, c("-", "-", "-"))
})

test_that("rule 2: other chromosomes and contigs are NA (toy chain)", {
  out <- lift_toy(c("3", "8", "7"), c(10, 10, 10))
  expect_identical(out$lift_status, c("other_chrom", "contig", "ok"))
  # a same-chromosome hit wins over a patch hit
  expect_identical(out$POS_b37, c(NA, NA, 20L))
  # non-standard input chromosomes are not looked up, even when the chain has them
  bad <- lift_toy(c("HG1_PATCH", "chrUn_gl000220", NA, ""), 10)
  expect_identical(bad$lift_status, rep("bad_input", 4))
})

test_that("rule 3: identical hits merge, distinct same-chromosome hits are NA (toy chain)", {
  out <- lift_toy(c("6", "5"), c(10, 10))
  expect_identical(out$POS_b37, c(10L, NA))
  expect_identical(out$lift_status, c("ok", "multiple"))
})

test_that("chromosome names are normalised, output keeps the input notation (toy chain)", {
  x <- lift_toy(c("X", "chrX", "23", "25", "XY", "MT", "M", "chrM", "26"), 10)
  expect_identical(x$POS_b37, c(rep(1010L, 5), rep(10L, 4)))
  expect_identical(x$CHR_b37, c("X", "chrX", "23", "25", "XY", "MT", "M", "chrM", "26"))
})

test_that("UCSC-style chains (chr1, chrM) are read too", {
  out <- liftover(data.frame(CHR = c("1", "chr1", "MT", "26", "chr1_gl000191_random"), POS = 5), from = 18, to = 19,
                  liftover.chain = toy_chain_ucsc(), silent = TRUE)
  expect_identical(out$POS_b37, c(15L, 15L, 996L, 996L, NA))
  expect_identical(out$STRAND_b37, c("+", "+", "-", "-", NA))
  expect_identical(out$lift_status, c(rep("ok", 4), "bad_input"))
})

test_that("zero rows give the output columns", {
  out <- lift_bundled(character(), numeric())
  expect_identical(nrow(out), 0L)
  expect_identical(names(out), c("CHR", "POS", "POS_b37", "CHR_b37", "STRAND_b37", "lift_status"))
})

test_that("bad positions are reported, not lifted", {
  out <- lift_toy("1", c(0, -5, 150.5, NA, 1e12, 150))
  expect_identical(out$lift_status, c(rep("bad_input", 5), "ok"))
  expect_identical(out$POS_b37, c(rep(NA_integer_, 5), 250L))
  chr_pos <- lift_toy("1", c("150", "abc"))
  expect_identical(chr_pos$lift_status, c("ok", "bad_input"))
})

test_that("rule 4: status column and summary note", {
  df <- data.frame(CHR = 1, POS = c(1005000, 143185000, 1e12))
  expect_message(liftover(df, from = 38, to = 37), "ok=1 \\| contig=1 \\| bad_input=1")
  expect_silent(liftover(df, from = 38, to = 37, silent = TRUE))
  expect_silent(suppressMessages(liftover(df, from = 38, to = 37)))
  expect_output(suppressMessages(liftover(df, from = 38, to = 37)), NA)
})

test_that("the same build is the identity", {
  out <- lift_bundled(c("1", "1"), c(123, -1), from = 38, to = "hg38")
  expect_identical(out$POS_b38, c(123L, NA))
  expect_identical(out$lift_status, c("ok", "bad_input"))
})

test_that("a liftover that cannot run warns (or errors with strict) instead of passing df through", {
  df <- data.frame(SNP = "a", CHR = 1, POS = 1005000)
  missing_chain <- file.path(tempdir(), "no-such.chain.gz")
  expect_warning(out <- liftover(df, from = 38, to = 37, liftover.chain = missing_chain, silent = TRUE),
                 "chain file not found")
  expect_identical(out$lift_status, "failed")
  expect_true(is.na(out$POS_b37))
  expect_error(liftover(df, from = 38, to = 37, liftover.chain = missing_chain, strict = TRUE, silent = TRUE),
               "chain file not found")
  expect_warning(liftover(df, from = 38, to = 99, silent = TRUE), "unknown build")

  broken <- tempfile(fileext = ".chain")
  writeLines(c("chain 1 1 1000 + 0 100 1 1000 + 0 200 1", "100"), broken)
  expect_error(liftover(df, from = 38, to = 37, liftover.chain = broken, strict = TRUE, silent = TRUE),
               "malformed chain")

  expect_warning(liftover(df, from = 38, to = 37, method = "binary", liftover = "/no/such/liftOver", silent = TRUE),
                 "liftOver binary not found: /no/such/liftOver")
})

test_that("hg18 -> hg19 is not bundled and points to chain.download()", {
  skip_if(file.exists(file.path(tools::R_user_dir("gcanvas", "data"), "hg18ToHg19.over.chain.gz")))
  expect_warning(liftover(data.frame(CHR = 1, POS = 1e6), from = 18, to = 19, silent = TRUE), "chain.download")
})

test_that("chain.download() needs the UCSC terms accepted and knows its pairs", {
  d <- tempfile("chains-")
  expect_error(chain.download(18, 19, dir = d), "non-commercial")
  expect_false(dir.exists(d))
  expect_error(chain.download(37, 38, dir = d, accept.license = TRUE), "no downloadable chain")
})

test_that("rcoloc call: binary path given, chain used", {
  df <- data.frame(SNP = c("a", ".start", ".end"), CHR = "1", POS = c(1005000, 905000, 1105000))
  out <- liftover(df, from = 38, to = 37, liftover = "/no/such/liftOver", liftover.chain = "auto", silent = TRUE)
  expect_identical(out$POS_b37, c(940380L, 840380L, 1040380L))
  expect_true(all(out$lift_status == "ok"))
})

test_that("the binary path agrees with the chain path", {
  bin <- liftover_bin()
  skip_if(!nzchar(bin), "UCSC liftOver not found")
  toy <- toy_chain()
  chr <- rep(c("1", "2", "3", "5", "6", "7", "8", "X", "25", "MT", "HG1_PATCH"), each = 6)
  pos <- rep(c(1, 50, 100, 101, 251, 300), 11)
  a <- liftover(data.frame(CHR = chr, POS = pos), 38, 37, liftover.chain = toy, silent = TRUE)
  b <- liftover(data.frame(CHR = chr, POS = pos), 38, 37, liftover.chain = toy, method = "binary",
                liftover = bin, silent = TRUE)
  expect_identical(a, b)
  u <- data.frame(CHR = c("1", "MT", "2"), POS = c(5, 5, 5))
  expect_identical(liftover(u, 18, 19, liftover.chain = toy_chain_ucsc(), silent = TRUE),
                   liftover(u, 18, 19, liftover.chain = toy_chain_ucsc(), method = "binary", liftover = bin, silent = TRUE))

  set.seed(1)
  n <- 5000
  d <- data.frame(CHR = sample(c(1:22, "X", "Y"), n, replace = TRUE), POS = sample.int(1.5e8, n))
  for (dir in list(c(38, 37), c(37, 38))) {
    a <- liftover(d, dir[1], dir[2], silent = TRUE)
    b <- liftover(d, dir[1], dir[2], method = "binary", liftover = bin, silent = TRUE)
    expect_identical(a, b)
  }
})
