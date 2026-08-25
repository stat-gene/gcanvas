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

> **Status:** v0.1.0 — early release. The API is stable enough to use day to
> day; `circos()` is experimental and being actively iterated. Detailed
> documentation and vignettes will follow.

## Install

```r
# Latest from main
remotes::install_github("stat-gene/gcanvas")

# Pin a tagged release
remotes::install_github("stat-gene/gcanvas@v0.1.0")
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
| UCSC `liftOver` | `liftover` | UCSC binary | UCSC binary | **WSL required** |

The package itself loads on every OS (CI tests macOS / Linux / Windows). On
Windows native, the LD / PCA / Manhattan / Q-Q / regional (LD-only) / circos
paths work; GTF and liftover paths need WSL or a Linux/macOS shell.

#### Downloading the binaries

**PLINK.** `gcanvas` uses **PLINK 2** for most genotype operations and
additionally needs **PLINK 1.9** for the multi-input merge in `plink.extract()`.
Download the build for your platform and put it on your `PATH` (or pass the
path via the function's `plink` argument):

- PLINK 2: <https://www.cog-genomics.org/plink/2.0/>
- PLINK 1.9: <https://www.cog-genomics.org/plink/1.9/>

**UCSC `liftOver`.** The chain files for the common build pairs
(GRCh37 ⟷ GRCh38, hg18 → hg19) ship with the package, so `liftover()` finds them
automatically. You only need the `liftOver` executable, which is **not**
bundled. The recommended way to install it is via conda:

```sh
conda install bioconda::ucsc-liftover
```

Alternatively, download the UCSC binary directly (Linux x86_64 example):

```sh
wget https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
chmod +x ./liftOver
```

Binaries for other platforms are under
<https://hgdownload.soe.ucsc.edu/admin/exe/>.

When you don't pass the `liftover` argument, the executable is resolved
automatically from your `PATH` and common install locations (including conda
environments), so a conda-installed `liftOver` is picked up with no extra setup:

```r
# liftOver auto-resolved from PATH / conda env:
liftover(df, from = 37, to = 38)
# or point at the executable explicitly:
liftover(df, from = 37, to = 38, liftover = "./liftOver")
```

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
| `liftover()` | Coordinate liftover via UCSC `liftOver` |
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
