# regional(build.gwas != build) lifts the variants, lead(s) and window ends in
# one liftover() call and draws only rows with lift_status "ok".

regional_quiet <- function(...) suppressMessages(regional(..., show.gene = FALSE))

toy_region <- function(chr, from, to, n, seed = 1) {
  set.seed(seed)
  data.frame(SNP = sprintf("v%d", seq_len(n)), CHR = chr, POS = sort(sample(seq(from, to), n)), P = stats::runif(n))
}

test_that("variants, lead and window are lifted without a liftOver binary", {
  d <- toy_region(1, 0.9e6, 1.1e6, 200)
  d$P[100] <- 1e-10
  g <- regional_quiet(d, snp = "v100", flank = 5e4, build = 38, build.gwas = 37, liftover = "/no/such/liftOver")
  m <- attr(g, "gcanvas_meta")
  ends <- lift_bundled("1", d$POS[100] + c(-5e4, 5e4), from = 37, to = 38)$POS_b38
  expect_equal(unname(m$region$position), as.numeric(ends))
  lifted <- liftover(d, from = 37, to = 38, silent = TRUE)
  keep <- lifted$lift_status == "ok" & lifted$POS_b38 >= ends[1] & lifted$POS_b38 <= ends[2]
  expect_equal(sort(m$ld_r2$POS), sort(as.numeric(lifted$POS_b38[keep])))
  expect_identical(m$ld_r2$SNP[m$ld_r2$is_lead], "v100")
  expect_identical(m$build$liftover$status[["ok"]], 200L)
})

test_that("a region missing from the Ensembl chain keeps its variants", {
  # GRCh37 chr17:35.9-36.3 Mb around rs4430796 (HNF1B)
  d <- toy_region(17, 35.9e6, 36.3e6, 200)
  d$POS[100] <- 36098040
  d$SNP[100] <- "rs4430796"
  d$P[100] <- 1e-12
  g <- regional_quiet(d, snp = "rs4430796", pos.range = c(35.9e6, 36.3e6), build = 38, build.gwas = 37)
  m <- attr(g, "gcanvas_meta")
  expect_identical(m$build$liftover$status[["ok"]], 200L)
  expect_identical(nrow(m$ld_r2), 200L)
  expect_identical(m$ld_r2$POS[m$ld_r2$is_lead], 37738049)
  # back from GRCh38 too
  d38 <- toy_region(17, 37.6e6, 37.9e6, 100)
  d38$POS[50] <- 37738049
  d38$SNP[50] <- "rs4430796"
  d38$P[50] <- 1e-12
  m <- attr(regional_quiet(d38, snp = "rs4430796", build = 37, build.gwas = 38), "gcanvas_meta")
  expect_identical(m$ld_r2$POS[m$ld_r2$is_lead], 36098040)
})

test_that("duplicated SNP ids keep their own lifted positions", {
  d <- toy_region(1, 0.95e6, 1.05e6, 50)
  d$SNP[c(10, 40)] <- "dup"
  d$P[25] <- 1e-9
  g <- regional_quiet(d, snp = d$SNP[25], flank = 1e5, build = 38, build.gwas = 37)
  m <- attr(g, "gcanvas_meta")
  want <- lift_bundled("1", d$POS[c(10, 40)], from = 37, to = 38)$POS_b38
  expect_equal(sort(m$ld_r2$POS[m$ld_r2$SNP == "dup"]), sort(as.numeric(want)))
})

test_that("untrusted variants are not drawn and an unliftable window end is replaced", {
  # GRCh38 chr1:143.30-143.40 Mb maps only to GRCh37 patch contigs
  d <- toy_region("chr1", 143.30e6, 143.80e6, 400)
  lead <- d$SNP[which.min(abs(d$POS - 143.7e6))]
  d$P[d$SNP == lead] <- 1e-15
  suppressMessages(expect_message(
    g <- regional(d, snp = lead, pos.range = c(143.30e6, 143.80e6), build = 37, build.gwas = 38, show.gene = FALSE),
    "pos.range end\\(s\\) could not be lifted"
  ))
  m <- attr(g, "gcanvas_meta")
  lifted <- liftover(d, from = 38, to = 37, silent = TRUE)
  ok <- lifted$lift_status == "ok"
  expect_gt(sum(!ok), 0L)
  expect_equal(sort(m$ld_r2$POS), sort(as.numeric(lifted$POS_b37[ok])))
  ends <- lift_bundled("chr1", c(143.30e6, 143.80e6), from = 38, to = 37)
  expect_identical(ends$lift_status, c("contig", "ok"))
  expect_equal(unname(m$region$position), range(as.numeric(c(ends$POS_b37[2], lifted$POS_b37[ok]))))
  expect_identical(m$build$liftover$status[["contig"]], sum(lifted$lift_status == "contig"))
})

test_that("a lead that cannot be lifted falls back to no-lead mode", {
  d <- toy_region("1", 143.30e6, 143.80e6, 400)
  lead <- d$SNP[which.min(abs(d$POS - 143.35e6))]
  suppressMessages(expect_message(
    g <- regional(d, snp = lead, pos.range = c(143.30e6, 143.80e6), build = 37, build.gwas = 38, show.gene = FALSE),
    "could not be lifted to b37 \\(contig\\); proceeding without highlighting"
  ))
  expect_false(any(attr(g, "gcanvas_meta")$ld_r2$is_lead))
})

test_that("regional() notes go to stderr only", {
  d <- toy_region(1, 0.9e6, 1.1e6, 50)
  suppressMessages(expect_message(regional(d, snp = "v25", build = 38, build.gwas = 37, show.gene = FALSE),
                                  "\\[gcanvas::regional\\]"))
  expect_output(regional_quiet(d, snp = "v25", build = 38, build.gwas = 37), NA)
})
