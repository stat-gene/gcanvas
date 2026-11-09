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
