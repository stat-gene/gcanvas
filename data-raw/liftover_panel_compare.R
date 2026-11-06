# Compares the chains on a variant panel: `auto` (NCBI, then Ensembl) against
# the Ensembl and UCSC chains alone and UCSC filled by Ensembl. The UCSC chains
# are not bundled (non-commercial use only), so their directory is an argument.
#
#   Rscript data-raw/liftover_panel_compare.R <panel> <ucsc chain dir>
#
# <panel> needs columns CHR, SNP, BP on GRCh37.

suppressPackageStartupMessages({
  library(data.table)
  pkgload::load_all(".", quiet = TRUE)
})
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L, file.exists(args[1]), dir.exists(args[2]))
panel <- fread(args[1], select = c("CHR", "SNP", "BP"))
ucsc <- file.path(args[2], "hg19ToHg38.over.chain.gz")
ucsc_back <- file.path(args[2], "hg38ToHg19.over.chain.gz")
ens <- system.file("extdata", "GRCh37_to_GRCh38.chain.gz", package = "gcanvas")

up <- function(chain = "auto") liftover(panel, 37, 38, pos.col = "BP", liftover.chain = chain, silent = TRUE)
a <- up()
e <- up(ens)
u <- up(ucsc)
d <- data.table(panel, pa = a$BP_b38, pe = e$BP_b38, pu = u$BP_b38, sa = a$lift_status, se = e$lift_status, su = u$lift_status)

cat("Ensembl vs UCSC:\n")
d[, cat := fcase(se == "ok" & su == "ok" & pe == pu, "same position",
                 se == "ok" & su == "ok", "different position",
                 se == "ok", "Ensembl only", su == "ok", "UCSC only", default = "neither")]
print(d[, .N, by = cat][order(-N)])

cat("\nRows not lifted to GRCh38:\n")
print(data.table(
  chains = c("Ensembl", "UCSC", "UCSC, then Ensembl", "auto (NCBI, then Ensembl)"),
  n = c(d[se != "ok", .N], d[su != "ok", .N], d[su != "ok" & se != "ok", .N], d[sa != "ok", .N])
))
cat("\nauto status:\n")
print(table(d$sa))
cat("\nauto vs UCSC, both lifted: same", d[sa == "ok" & su == "ok" & pa == pu, .N],
    "| different", d[sa == "ok" & su == "ok" & pa != pu, .N], "\n")

cat("\nBack to GRCh37 from the auto position:\n")
ok <- d[sa == "ok"]
back <- function(chain = "auto") {
  liftover(ok[, .(CHR, POS = pa)], 38, 37, liftover.chain = chain, silent = TRUE)
}
b <- back()
bu <- back(ucsc_back)
print(data.table(
  chains = c("auto", "UCSC"), n = nrow(ok),
  original = c(sum(b$POS_b37 == ok$BP, na.rm = TRUE), sum(bu$POS_b37 == ok$BP, na.rm = TRUE)),
  other = c(sum(b$POS_b37 != ok$BP, na.rm = TRUE), sum(bu$POS_b37 != ok$BP, na.rm = TRUE)),
  not_lifted = c(sum(is.na(b$POS_b37)), sum(is.na(bu$POS_b37)))
))
