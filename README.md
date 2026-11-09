# gcanvas

<!-- badges: start -->
[![R-CMD-check](https://github.com/stat-gene/gcanvas/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/stat-gene/gcanvas/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

A visualization and analysis toolkit for GWAS workflows:

Manhattan / Q-Q
/ regional / Miami / circos / heterozygosity-missingness plots, PCA with
reference projection and ancestry estimation, LD computation / clumping /
proxy lookup, GTF-driven gene tracks, liftover, and a set of helpers for
working with summary statistics.

> **Status:** v0.1.2 — early release. The API is stable enough to use day to
> day; `circos()` is experimental and being actively iterated. Detailed
> documentation and vignettes will follow.

## Install

![Version](https://img.shields.io/badge/Version-0.1.2-blue) ![Last Updated](https://img.shields.io/badge/Last%20Updated-Nov%2006%202026-red)

```r
# Latest from main
remotes::install_github("stat-gene/gcanvas")

# Pin a tagged release
remotes::install_github("stat-gene/gcanvas@v0.1.2")
```

Heavy R dependencies (`bigsnpr`, `Rsamtools`, `rcartocolor`, ...) are listed
as `Suggests`; install them when you need the functions that touch them.

### External binaries

Some functions shell out to standard genomics tools. Install whichever you
need for your workflow:

| Tool | Used by | macOS | Linux | Windows |
|---|---|---|---|---|
| `plink2` | `pca`, `pca.projection`, `plink.extract`, `hetmiss` | `brew install plink2` | apt / conda / native binary | Native binary on `PATH` |
| `plink` (1.9) | `plink.extract` (multi-input merge), `hetmiss` (`plink.version="plink"`) | `brew install plink` | apt / conda / native binary | Native binary on `PATH` |
| `tabix` (htslib) | `gtf2rds`, `regional` (GTF mode), `geneinfo` (GTF mode) | `brew install htslib` | `apt install tabix` / conda | **WSL required** |
| UCSC `liftOver` (optional) | `liftover(method = "binary")` only | conda / UCSC binary | conda / UCSC binary | not available |

The package itself loads on every OS (CI tests macOS / Linux / Windows). On
Windows native, the LD / PCA / Manhattan / Q-Q / regional (LD-only) / circos
and liftover paths work; GTF paths need WSL or a Linux/macOS shell.
`liftover()` reads chain files in R and needs no external program (see
[Liftover](#liftover)).

#### Downloading the binaries

**PLINK.** `gcanvas` uses **PLINK 2** for most genotype operations and
additionally needs **PLINK 1.9** for the multi-input merge in `plink.extract()`.
Download the build for your platform and put it on your `PATH` (or pass the
path via the function's `plink` argument):

- PLINK 2: <https://www.cog-genomics.org/plink/2.0/>
- PLINK 1.9: <https://www.cog-genomics.org/plink/1.9/>

### Liftover

`liftover()` reads chain files directly in R: no `liftOver` binary, no WSL.
The GRCh37 ⟷ GRCh38 chains ship with the package and are picked by build pair:

```r
lifted <- liftover(sumstats, from = 37, to = 38)   # adds POS_b38, CHR_b38, STRAND_b38, lift_status
lifted[lifted$lift_status == "ok", ]
```

It attaches only positions it can trust:

1. **One output row per input row**, in the same order, attached by row (never
   by SNP id), so duplicated or missing ids cannot mix positions up.
2. **Same chromosome only.** A position that moves to another chromosome, or
   only to a patch / alt / unplaced contig (`HG*_PATCH`, `*_alt`, `GL*`, `KI*`,
   ...), gets `NA`; a contig coordinate is never attached to the original
   chromosome.
3. **Several hits** (liftOver's "Duplicated in new"): identical hits merge;
   exactly one hit on the same chromosome is used, none or several distinct
   ones give `NA`.
4. **Every row gets a `lift_status`** (`ok`, `unmapped`, `other_chrom`,
   `contig`, `multiple`, `not_reversible`, `bad_input`), and the counts are
   summarized.

`STRAND_b38 == "-"` means the base is reverse-complemented in the new build,
so alleles must be complemented. `regional(..., build.gwas = 37, build = 38)`
uses the same rules and draws only `ok` variants.

GRCh37 ⟷ GRCh38 use two bundled chains, both free of use restrictions, tried in
order: one built from the NCBI Remap first-pass alignments, then the Ensembl
assembly-mapping chain for rows the first left without a position on their
chromosome. The Ensembl chain alone has no blocks for several regions (e.g.
GRCh37 chr17:34-37 Mb, which holds *HNF1B*, and chr10:47.5-48.7 Mb): on the
1.12M-SNP 1000 Genomes HapMap3 panel it leaves 1,431 SNPs unlifted, against 20
for the two chains. A row filled by the Ensembl chain must also lift back to
where it came from, otherwise it gets `NA` (`not_reversible`: two places of the
old build land on one base), so every `ok` row returns to its input position.
The summary note shows the rows each chain lifted and points out rows still
`unmapped`. To use one chain only, or a chain of your own, pass `liftover.chain`
(paths are tried in order, without the lift-back check).

hg18 → hg19 uses a UCSC chain. UCSC chain files are free for non-commercial use
only, so it is not bundled; download it once (you are asked to accept the UCSC
terms):

```r
chain.download(18, 19)
liftover(df, from = 18, to = 19)
```

The UCSC `liftOver` binary is optional. `liftover(method = "binary")` runs
`liftOver -multiple` and applies the same rules, which is useful as a
cross-check (the results are identical). The binary is free for non-commercial
use only (commercial use needs a UCSC license) and is not bundled; if you are
eligible, install it yourself (`conda install bioconda::ucsc-liftover`, or
<https://hgdownload.soe.ucsc.edu/admin/exe/>) and pass `liftover = "<path>"`,
set `options(gcanvas.liftover = "<path>")` or the `GCANVAS_LIFTOVER`
environment variable, or put it on `PATH`. See
`inst/extdata/CHAIN_FILES_LICENSE.md` for the chain-file terms.

## What's inside

### Plots
| Function | Purpose |
|---|---|
| `manhattan()` | Genome-wide Manhattan plot of association p-values |
| `qq()` | Q-Q plot of p-values with optional category stratification |
| `regional()` | Locus zoom plot (variants + LD coloring + gene track) |
| `regional.track()` | Build the gene/exon layout used inside `regional()` (also standalone) |
| `miami()` | Two stacked, mirrored Manhattan panels for comparing two results |
| `circos()` | Circos-style genome-wide plot (**experimental**) |
| `hetmiss()` | Per-sample heterozygosity vs. missingness QC scatter |

### Plot utilities
| Function | Purpose |
|---|---|
| `get.legend()` | Extract the legend grob from a `ggplot` |
| `split.plot_label()` | Split a plot into point-only + label-only layers (rasterise the points, keep labels vector) |

### PCA
| Function | Purpose |
|---|---|
| `pca()` | PCA via PLINK2 on a bfile/pfile |
| `pca.projection()` | Project a target dataset onto a reference PCA |
| `pca.plot()` | 2D / 3D / facet-grid PC scatter, with categories and projection overlay |
| `pca.screeplot()` | Eigenvalue / proportion-variance scree plot |

### PCA ancestry estimation
| Function | Purpose |
|---|---|
| `pca.refqc()` | Reference-PCA QC: drop within-population outliers (adaptive `k`) |
| `pca.ancestry()` | Assign / score target samples against QC'd reference labels |

### LD / cache
| Function | Purpose |
|---|---|
| `calcld()` | LD (signed r and r-squared) for single / pairwise / set / lead-vs-region modes |
| `ldproxy()` | Proxy SNPs in LD with a query at a given r-squared threshold |
| `ldclump()` | PLINK-style clumping of summary statistics |
| `cache.list()` | Inspect the in-memory / on-disk LD reference cache |
| `cache.reset()` | Clear the LD reference cache |

### Genomic analysis helpers
| Function | Purpose |
|---|---|
| `liftover()` | Coordinate liftover between builds from chain files (same-chromosome hits only, per-row `lift_status`) |
| `marker.id()` | Marker ID `chrom:pos:a:b` with upper-cased, sorted alleles; the same ID for integer, numeric or character positions; warns about missing or unusual allele codes |
| `allele.sort()` | The two alleles of each row in `marker.id()` order, for joining tables on `CHR`, `POS` and alleles without ID strings |
| `chain.download()` | Fetch the UCSC hg18 → hg19 chain (non-commercial terms) for `liftover()` |
| `is.novel()` | Mark hits as novel vs. a reported set (with distance threshold) |
| `get.lead()` | Greedy lead-variant detection over summary stats |
| `beta.align()` | Align effect sizes to a reference allele |
| `maf.align()` | Align minor allele frequencies to a reference allele |
| `ld.align()` | Sign-flip LD r values (pair tables or matrices) to a reference allele |
| `plink.extract()` | PLINK2 subset-and-filter wrapper |
| `get.pilot()` | Lead-preserving downsample of summary stats for fast pilot plots |
| `geneinfo()` | Gene metadata lookup by symbol / Ensembl id |

### Utilities
| Function | Purpose |
|---|---|
| `get.colors()` | Flexible palette generator (random / Brewer / Carto / curated bases) |
| `show.pal()` | Display a palette by name or vector |
| `sort.chrom()` | Sort vectors / tables in natural chromosome order |
| `rank.chrom()` | Integer rank for chromosome ordering |
| `normalize.chrom()` | Canonicalize chromosome names (strip `chr`, remap `23 → X` ...) |
| `pvalue()` | Two-/one-tailed p-values from z (or beta/se); switches to scientific-notation character for values below `tiny.threshold` so tiny p never underflows |
| `format_pvalue()` | Format p-values (numeric / char / `log10(p)`) as `geom_richtext`/Markdown HTML (default) or plain `"3.20 x 10^-8"` (`html = FALSE`) |
| `log10c()` | `log10()` that accepts character p-values in scientific notation (handles values smaller than `5e-324`) |
| `pow10c()` | Inverse of `log10c()`: turns `log10(p)` into character scientific notation, never underflows (`pow10c(-12012) -> "1e-12012"`) |
| `zabs()` | `\|beta / se\|` (or `\|beta / sqrt(varbeta)\|`) |
| `quantvec()` | Quantile bucket labels (factor of `"1st", "2nd", ...`) |
| `credibleset()` | Credible-set labels from posterior inclusion probabilities |
| `invnorm()` | Rank-based inverse normal transform |
| `logger()` | Lightweight console + file logger |
| `dt2mat()` | Pair table → symmetric matrix |
| `mat2dt()` | Matrix → long-form pair table |
| `gtf2rds()` | Convert bgzipped + tabix-indexed GTF into a `gcanvas` annotation cache |

### Toy data generators
| Function | Purpose |
|---|---|
| `toy.gwas()` | Synthetic GWAS summary stats + matching PLINK1 bfile (single- or multi-ancestry); used in every example below |
| `toy.eqtl()` | Synthetic eQTL summary stats; can piggyback on a `toy.gwas()` result to share variants and bfile |

Built-in datasets: `tracks.b37`, `tracks.b38` (lazy-loaded annotation tracks).

## Minimal example

All examples below are self-contained — `toy.gwas()` writes a small PLINK1
bfile on disk so the LD / PCA / het-miss steps work without external data.

```r
library(gcanvas)

# ---------------------------------------------------------------------------
# 1. Generate a toy GWAS (summary stats + matching PLINK1 bfile "toy")
# ---------------------------------------------------------------------------
sumstats <- toy.gwas(
  n.sample = 5000, n.snp = 10000, n.causal = 3,
  seed = 23L, bfile = "toy"
)
leads <- get.lead(sumstats)$leadSNP

# ---------------------------------------------------------------------------
# 2. Genome-wide diagnostic plots
# ---------------------------------------------------------------------------
manhattan(sumstats, lead = leads, lead.label.col = "SNP")
qq(sumstats, lambda.gc = TRUE)

# ---------------------------------------------------------------------------
# 3. Regional locus zoom (LD coloring from the toy bfile)
# ---------------------------------------------------------------------------
regional(sumstats, snp = leads[1], ld.bfile = "toy", flank = 2.5e5)

# ---------------------------------------------------------------------------
# 4. LD utilities
# ---------------------------------------------------------------------------
ld_pairs <- calcld(snp = leads[1], flank = "10kb", ld.bfile = "toy")
clumps   <- ldclump(sumstats, ld.bfile = "toy", rsq = 0.2)
proxies  <- ldproxy(leads[1], ld.bfile = "toy", rsq = 0.8)
is.novel(leads, reported = leads[1], distance = "1Mb")

# ---------------------------------------------------------------------------
# 5. Two-study comparison: Miami plot
# ---------------------------------------------------------------------------
sumstats2 <- toy.gwas(
  n.sample = 5000, n.snp = 10000, n.causal = 2,
  seed = 42L, bfile = "toy2"
)
m1 <- manhattan(sumstats,  lead = leads,        lead.label.col = "SNP")
m2 <- manhattan(sumstats2, lead.label.col = "SNP")
miami(m1, m2)

# ---------------------------------------------------------------------------
# 6. eQTL view of the same locus (toy.eqtl reuses the GWAS scaffold)
# ---------------------------------------------------------------------------
eqtl <- toy.eqtl(base = sumstats, n.sample = 500, seed = 23L)
regional(
  data = eqtl[eqtl$gene_is_focal],
  snp  = leads[1], ld.bfile = "toy", flank = 2.5e5
)

# ---------------------------------------------------------------------------
# 7. Sample-level QC and PCA on the toy bfile
# ---------------------------------------------------------------------------
hetmiss(bfile = "toy")$plot
pc <- pca(bfile = "toy")
pca.plot(pc, pc.use = 1:4)$plot

# ---------------------------------------------------------------------------
# 8. Pretty p-value labels for ggtext / ggplot annotations
# ---------------------------------------------------------------------------
format_pvalue(c(0.5, 5e-8, "1.2e-400"))  # HTML for geom_richtext()
format_pvalue(5e-8, html = FALSE)        # "5 x 10^-8" plain text
```

Per-function help: `?manhattan`, `?regional`, `?pca`, `?toy.gwas`, ...

For the multi-ancestry / PCA-projection / ancestry-estimation workflow
(`pca.refqc()` + `pca.ancestry()`), see `?pca.refqc` — that pipeline assumes
a merged reference/target PLINK fileset and is briefly sketched in the
package vignettes (forthcoming).

## License

MIT © 2026 Beomsu Kim
