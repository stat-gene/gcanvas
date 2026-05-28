# gcanvas

<!-- badges: start -->
[![R-CMD-check](https://github.com/stat-gene/gcanvas/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/stat-gene/gcanvas/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

Genomic visualization and analysis toolkit for GWAS workflows: Manhattan / Q-Q
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

Heavy dependencies (`bigsnpr`, `Rsamtools`, `rcartocolor`, ...) are listed as
`Suggests`; install them when you need the functions that touch them.

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
| `split.label()` | Wrap / thin long category labels so high-cardinality legends fit |

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
| `pvalue()` | Two-/one-tailed p-values from z (or beta/se), with tiny-p string fallback |
| `zabs()` | `|beta / se|` (or `|beta / sqrt(varbeta)|`) |
| `quantvec()` | Quantile bucket labels (factor of `"1st", "2nd", ...`) |
| `credibleset()` | Credible-set labels from posterior inclusion probabilities |
| `invnorm()` | Rank-based inverse normal transform |
| `logger()` | Lightweight console + file logger |
| `dt2mat()` | Pair table → symmetric matrix |
| `mat2dt()` | Matrix → long-form pair table |
| `gtf2rds()` | Convert bgzipped + tabix-indexed GTF into a `gcanvas` annotation cache |

Built-in datasets: `tracks.b37`, `tracks.b38` (lazy-loaded annotation tracks).

## Minimal example

```r
library(gcanvas)

# Manhattan + Q-Q on a summary-stat data.table with SNP/CHR/POS/P columns
manhattan(sumstats)
qq(sumstats)

# Regional locus zoom around a lead SNP
regional(sumstats, snp = "rs12345", flank = 5e5)

# PCA + ancestry estimation workflow
ref  <- pca(bfile = "ref")
proj <- pca.projection(ref, target = "target")
qcd  <- pca.refqc(ref, label.data = ref_labels, label.col = "ancestry")
anc  <- pca.ancestry(proj, label.data = qcd$label.data, label.col = "ancestry")
```

Per-function help is available via `?manhattan`, `?regional`, etc.

## License

MIT © 2026 Beomsu Kim
