# Bundled liftOver chain files — sources and licenses

The chain files in this directory are redistributed third-party data. They are
**not** covered by the gcanvas package's MIT license; each retains the license
of its original source, listed below. They are used by `gcanvas::liftover()`
(and `liftover.chain = "auto"`) to translate coordinates between genome builds.

| File | Build mapping | Chromosome style | Source | License |
|---|---|---|---|---|
| `GRCh37_to_GRCh38.chain.gz` | GRCh37 → GRCh38 (hg19 → hg38) | Ensembl (`1`, `X`, `MT`) | Ensembl assembly mapping | Apache License 2.0 |
| `GRCh38_to_GRCh37.chain.gz` | GRCh38 → GRCh37 (hg38 → hg19) | Ensembl (`1`, `X`, `MT`) | Ensembl assembly mapping | Apache License 2.0 |
| `hg18ToHg19.over.chain.gz`  | hg18 → hg19 (NCBI36 → GRCh37) | UCSC (`chr1`, `chrX`, `chrM`) | UCSC Genome Browser | UCSC data terms (see below) |

## Ensembl (GRCh37 <-> GRCh38)

Source: Ensembl assembly-mapping chain files
(`ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/`).

Ensembl data and code are released under the Apache License 2.0, with no
additional restrictions on use. See <https://www.ensembl.org/info/about/legal/>.
A copy of the Apache 2.0 license is available at
<https://www.apache.org/licenses/LICENSE-2.0>.

## UCSC (hg18 -> hg19)

Source: UCSC Genome Browser liftOver chain files
(`hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/`).

The UCSC genome annotation and assembly data (including liftOver chain files)
are freely available for any use, including commercial use, for the human
assemblies. The commercial-license requirement applies to the UCSC Genome
Browser / Blat **software** (e.g. the `liftOver` binary), which is *not*
bundled here. When redistributing UCSC data, please credit the UCSC Genome
Browser. See <https://genome.ucsc.edu/license/> and
<https://genome.ucsc.edu/conditions.html>.

> Kent WJ, et al. The Human Genome Browser at UCSC. Genome Res. 2002.

## Note on the liftOver binary

The `liftOver` executable is **not** distributed with this package. Install it
yourself (e.g. from UCSC) and either place it in `liftover.dir` or make it
available on your `PATH`. Commercial users should review the UCSC software
license for the binary.
