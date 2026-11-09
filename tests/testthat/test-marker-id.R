test_that("integer, numeric and character positions give the same ID", {
  pos <- c(123456L, 100000L, 15000000L, 248956422L)
  want <- c("1:123456:A:G", "1:100000:A:G", "1:15000000:A:G", "1:248956422:A:G")
  expect_identical(marker.id(1, pos, "G", "A"), want)
  expect_identical(marker.id(1, as.numeric(pos), "G", "A"), want)
  expect_identical(marker.id("1", as.character(pos), "G", "A"), want)
  expect_identical(marker.id(1, factor(as.character(pos)), "G", "A"), want)
  expect_identical(marker.id(1, c("1e+05", " 100000", "0100000"), "A", "G"), rep("1:100000:A:G", 3))
  op <- options(scipen = -9)
  on.exit(options(op))
  expect_identical(marker.id(1, as.numeric(pos), "G", "A"), want)
})

test_that("allele order and case do not matter, and the order is byte-wise", {
  expect_identical(marker.id(1, 5, c("G", "a", "A", "g"), c("A", "G", "g", "A")), rep("1:5:A:G", 4))
  expect_identical(marker.id(1, 5, c("<DEL>", "-", "AT", "I"), c("*", "A", "A", "D")),
                   c("1:5:*:<DEL>", "1:5:-:A", "1:5:A:AT", "1:5:D:I"))
  expect_identical(marker.id(1, 5, factor("T"), factor("C")), "1:5:C:T")
})

test_that("missing parts give NA, bad input is an error", {
  expect_identical(marker.id(c(1, NA, 1, 1, 1), c(5, 5, NA, 5, 5), c("A", "A", "A", NA, "A"), c("G", "G", "G", "G", NA)),
                   c("1:5:A:G", NA, NA, NA, NA))
  expect_identical(marker.id(1, c("", ".", "NA"), "A", "G"), rep(NA_character_, 3))
  expect_identical(marker.id(character(), integer(), character(), character()), character())
  expect_error(marker.id(1, 1234.5, "A", "G"), "whole numbers")
  expect_error(marker.id(1, -1L, "A", "G"), "whole numbers")
  expect_error(marker.id(1, "12a", "A", "G"), "not a number")
  t_alleles <- utils::read.table(text = "A1 A2\nT C\nT G", header = TRUE)
  expect_error(marker.id(1, 5:6, t_alleles$A1, t_alleles$A2), "logical")
  expect_error(marker.id(1:2, 1:3, "A", "G"), "same length")
})

test_that("strand, chromosome normalization and separator", {
  expect_identical(marker.id(1, 100, c("A", "AT", "-", "A"), c("G", "A", "C", "T"), strand = "-"),
                   c("1:100:C:T", "1:100:AT:T", "1:100:-:G", "1:100:A:T"))
  expect_identical(marker.id(1, 100, "A", "G", strand = c(NA)), "1:100:A:G")
  expect_error(marker.id(1, 100, "A", "G", strand = "x"), "strand")
  expect_identical(marker.id(c("chr1", "23", "chrX", "M"), 5, "A", "G", normalize = TRUE),
                   c("1:5:A:G", "X:5:A:G", "X:5:A:G", "MT:5:A:G"))
  expect_identical(marker.id(c("chr1", "23"), 5, "A", "G"), c("chr1:5:A:G", "23:5:A:G"))
  expect_identical(marker.id(1, 5, "A", "G", sep = "_"), "1_5_A_G")
  expect_identical(marker.id(1, 5, "A", "G", sep = "%"), "1%5%A%G")
  expect_identical(marker.id(1, 3e9, "A", "G"), "1:3000000000:A:G")
})

test_that("IDs after liftover() use the lifted position and strand", {
  d <- data.frame(CHR = c(17, 15), POS = c(36098040, 22755185), A1 = c("A", "A"), A2 = c("G", "C"))
  out <- liftover(d, from = 37, to = 38, silent = TRUE)
  expect_identical(out$STRAND_b38, c("+", "-"))
  expect_identical(marker.id(out$CHR, out$POS_b38, out$A1, out$A2, strand = out$STRAND_b38),
                   c("17:37738049:A:G", "15:23117883:G:T"))
})

test_that("allele.sort() gives the alleles in marker.id() order", {
  out <- allele.sort(c("G", "a", "T", "<DEL>", NA, "A"), c("A", "G", "C", "*", "C", NA))
  expect_identical(out, list(lo = c("A", "A", "C", "*", NA, NA), hi = c("G", "G", "T", "<DEL>", NA, NA)))
  expect_identical(allele.sort("A", c("G", "C"), strand = c("-", "+")), list(lo = c("C", "A"), hi = c("T", "C")))
  expect_identical(allele.sort(factor("T"), factor("C")), list(lo = "C", hi = "T"))
  expect_identical(allele.sort(character(), character()), list(lo = character(), hi = character()))
  expect_error(allele.sort(TRUE, "C"), "logical")
  expect_error(allele.sort(c("A", "C"), c("G", "T", "A")), "same length")

  set.seed(1)
  al <- c("A", "C", "G", "T", "a", "AT", "-")
  d <- data.frame(CHR = sample(1:3, 2000, TRUE), POS = sample.int(1e6, 2000), A1 = sample(al, 2000, TRUE),
                  A2 = sample(al, 2000, TRUE), S = sample(c("+", "-", NA), 2000, TRUE))
  k <- allele.sort(d$A1, d$A2, strand = d$S)
  expect_identical(paste(d$CHR, d$POS, k$lo, k$hi, sep = ":"), marker.id(d$CHR, d$POS, d$A1, d$A2, strand = d$S))
})

test_that("allele.sort() columns join two tables like marker.id() strings do", {
  set.seed(2)
  al <- c("A", "C", "G", "T")
  x <- data.table::data.table(CHR = sample(1:2, 500, TRUE), POS = sample.int(1e4, 500), A1 = sample(al, 500, TRUE))
  x[, A2 := al[(match(A1, al) + sample(0:2, .N, TRUE)) %% 4 + 1]]
  y <- data.table::copy(x)[sample(.N, 300)][, c("A1", "A2") := .(A2, A1)]
  x[, c("A_lo", "A_hi") := allele.sort(A1, A2)]
  y[, c("A_lo", "A_hi") := allele.sort(A1, A2)]
  by_cols <- y[x, on = c("CHR", "POS", "A_lo", "A_hi"), nomatch = NULL]
  by_id <- intersect(marker.id(x$CHR, x$POS, x$A1, x$A2), marker.id(y$CHR, y$POS, y$A1, y$A2))
  expect_identical(nrow(by_cols), 300L)
  expect_setequal(marker.id(by_cols$CHR, by_cols$POS, by_cols$A_lo, by_cols$A_hi), by_id)
})
