# Bundled liftOver chain files — sources and licenses

The chain files in this directory are third-party data redistributed with
gcanvas. They are **not** covered by the gcanvas package's MIT license; each
keeps the terms of its source, listed below. `gcanvas::liftover()` reads them
directly (`liftover.chain = "auto"`) to translate coordinates between genome
builds; no external program is needed.

| File | Build mapping | Chromosome style | Source | Terms |
|---|---|---|---|---|
| `GRCh37_to_GRCh38.chain.gz` | GRCh37 → GRCh38 (hg19 → hg38) | Ensembl (`1`, `X`, `MT`; patches `HG*_PATCH`, `CHR_*`) | Ensembl assembly mapping | No restrictions on use |
| `GRCh38_to_GRCh37.chain.gz` | GRCh38 → GRCh37 (hg38 → hg19) | Ensembl (`1`, `X`, `MT`; patches `HG*_PATCH`, `CHR_*`) | Ensembl assembly mapping | No restrictions on use |

## Ensembl (GRCh37 <-> GRCh38)

Source: Ensembl assembly-mapping chain files
(`https://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/`).

Ensembl's disclaimer states that Ensembl places no restrictions on access to,
or use of, the data it provides
(<https://www.ensembl.org/info/about/legal/disclaimer.html>). The Apache
License 2.0 mentioned on the same page covers Ensembl's **code**, not these
data files.

## Not bundled: UCSC chain files (e.g. hg18 -> hg19)

UCSC's license page (<https://genome.ucsc.edu/license/>) makes the liftOver
chain files the one exception to its otherwise free data: they may be linked,
downloaded, used and redistributed for non-commercial use only. The README in
each UCSC liftOver directory limits free use to independent researchers and
nonprofit organizations, requires a commercial license for any other use, and
treats downloading or using the files as accepting the UCSC EULA
(<https://genome.ucsc.edu/license/EULA.pdf>).

gcanvas therefore does not ship UCSC chains. `gcanvas::chain.download(18, 19)`
fetches `hg18ToHg19.over.chain.gz` from
`https://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/` into
`tools::R_user_dir("gcanvas", "data")` (where `liftover.chain = "auto"`
looks) after you accept those terms (`accept.license = TRUE`, or the prompt in
an interactive session).

> Kent WJ, et al. The Human Genome Browser at UCSC. Genome Res. 2002.

## Not bundled: the liftOver binary

The UCSC `liftOver` executable is not needed by default and is **not**
distributed with gcanvas. Only `liftover(method = "binary")` uses it, as a
cross-check of the built-in chain reader. UCSC licenses it free for
non-commercial use; commercial use requires a UCSC license, and liftOver (with
BLAT) is excluded from the UCSC tools that may be redistributed through package
managers. If you are eligible, install it yourself (e.g.
`conda install bioconda::ucsc-liftover`, or
<https://hgdownload.soe.ucsc.edu/admin/exe/>) and point gcanvas at it with the
`liftover` argument, `options(gcanvas.liftover = ...)` or the
`GCANVAS_LIFTOVER` environment variable.
