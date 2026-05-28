# GTF preparation (bgzip + tabix), conversion to RDS, gene metadata lookup,
# and region queries. Exports `gtf2rds()` and `geneinfo()`.

.gcanvas_decompress_gz_to_file <- function(gz_path, out_path, chunk_size = 1024L * 1024L * 8L) {
  in_con <- gzfile(gz_path, open = "rb")
  on.exit(close(in_con), add = TRUE)
  out_con <- file(out_path, open = "wb")
  on.exit(close(out_con), add = TRUE)
  repeat {
    buf <- readBin(in_con, what = "raw", n = chunk_size)
    if (!length(buf)) break
    writeBin(buf, out_con)
  }
  invisible(out_path)
}

.gtf.sort_plain_system <- function(in_plain, out_plain, chr_order = c("natural", "lex")) {
  chr_order <- match.arg(chr_order)
  bash <- Sys.which("bash")
  if (!nzchar(bash)) stop("bash not found; cannot system-sort GTF.", call. = FALSE)

  in_plain <- normalizePath(in_plain, winslash = "/", mustWork = TRUE)
  out_plain <- normalizePath(out_plain, winslash = "/", mustWork = FALSE)

  if (chr_order == "lex") {
    cmd <- sprintf(
      'set -euo pipefail; (grep "^#" %s || true) > %s; (grep -v "^#" %s || true) | LC_ALL=C sort -k1,1 -k4,4n >> %s',
      shQuote(in_plain), shQuote(out_plain),
      shQuote(in_plain), shQuote(out_plain)
    )
  } else {
    cmd <- sprintf(
      'set -euo pipefail; (grep "^#" %s || true) > %s; ' %+%
        '(grep -v "^#" %s || true) ' %+%
        '| awk \'BEGIN{OFS="\\t"} {c=$1; sub(/^chr/i,"",c); u=toupper(c); ' %+%
        'r=1000; if(u=="X") r=23; else if(u=="Y") r=24; else if(u=="MT"||u=="M") r=25; ' %+%
        'else if(c ~ /^[0-9]+$/) r=c+0; print r,$0}\' ' %+%
        '| LC_ALL=C sort -k1,1n -k2,2 -k5,5n ' %+%
        '| cut -f2- >> %s',
      shQuote(in_plain), shQuote(out_plain),
      shQuote(in_plain), shQuote(out_plain)
    )
  }

  system2(bash, c("-lc", cmd), stdout = TRUE, stderr = TRUE)
  invisible(out_plain)
}

gtf_prepare_tabix <- function(gtf_gz,
                              out_bgz = NULL,
                              overwrite = FALSE,
                              keep_tmp = FALSE,
                              sort = c("auto", "never", "always"),
                              chr_order = c("natural", "lex")) {
  require_pkg("Rsamtools")
  sort <- match.arg(sort)
  chr_order <- match.arg(chr_order)

  if (!file.exists(gtf_gz)) stop("GTF not found: ", gtf_gz, call. = FALSE)

  if (is.null(out_bgz)) {
    base <- sub("\\.gz$", "", gtf_gz, ignore.case = TRUE)
    out_bgz <- paste0(base, ".bgz")
  }
  tbi <- paste0(out_bgz, ".tbi")

  if (!overwrite && file.exists(out_bgz) && file.exists(tbi)) return(abs_path(out_bgz))

  tmp_plain <- tempfile(fileext = ".gtf")
  if (grepl("\\.gz$", gtf_gz, ignore.case = TRUE)) {
    .gcanvas_decompress_gz_to_file(gtf_gz, tmp_plain)
  } else {
    file.copy(gtf_gz, tmp_plain, overwrite = TRUE)
  }
  tmp_sorted <- tempfile(fileext = ".sorted.gtf")

  attempt_index <- function(plain_path) {
    Rsamtools::bgzip(plain_path, dest = out_bgz, overwrite = TRUE)
    Rsamtools::indexTabix(out_bgz, format = "gff")
    TRUE
  }

  ok <- FALSE
  err_msg <- NULL

  if (sort != "always") {
    ok <- tryCatch(attempt_index(tmp_plain), error = function(e) { err_msg <<- conditionMessage(e); FALSE })
  }

  if (!ok) {
    if (sort == "never") stop(err_msg %||% "indexTabix failed.", call. = FALSE)
    message("Sorting GTF for tabix indexing...")
    .gtf.sort_plain_system(tmp_plain, tmp_sorted, chr_order = chr_order)
    ok2 <- tryCatch(attempt_index(tmp_sorted), error = function(e) { err_msg <<- conditionMessage(e); FALSE })
    if (!ok2) stop(err_msg %||% "indexTabix failed after sorting.", call. = FALSE)
  }

  if (!keep_tmp) unlink(c(tmp_plain, tmp_sorted))
  message("Prepared .bgz and .tbi.")
  abs_path(out_bgz)
}

#' Convert a bgzipped + tabix-indexed GTF into an `.rds` annotation cache
#'
#' Reads gene and exon records from `gtf.bgz` for each contig listed in the
#' tabix index, normalizes chromosome names, and writes a single `.rds` file
#' containing a list with `genes` and `exons` data.tables. The resulting cache
#' is what [regional.track()] and [geneinfo()] consume.
#'
#' @param gtf.bgz Path to a bgzipped GTF with a `.tbi` alongside.
#' @param build Genome build (`37` or `38`).
#' @param out.rds Output `.rds` path. If `NULL`, derived from `gtf.bgz`.
#' @param overwrite Logical. Overwrite an existing output file.
#' @param compress Compression used by [saveRDS()].
#'
#' @return Invisibly returns the path to the written `.rds`.
#' @export
gtf2rds <- function(gtf.bgz,
                    build = 38L,
                    out.rds = NULL,
                    overwrite = FALSE,
                    compress = "xz") {
  require_pkg("data.table")

  if (!file.exists(gtf.bgz)) stop("gtf.bgz not found: ", gtf.bgz, call. = FALSE)
  gtf.bgz <- abs_path(gtf.bgz)

  tbi1 <- paste0(gtf.bgz, ".tbi")
  tbi2 <- sub("\\.bgz$", ".tbi", gtf.bgz, ignore.case = TRUE)
  if (!file.exists(tbi1) && !file.exists(tbi2)) {
    stop("Tabix index (.tbi) not found for: ", gtf.bgz, call. = FALSE)
  }

  build <- as_int(build)
  if (is.na(build)) build <- 38L
  if (build == 19L) build <- 37L

  if (is.null(out.rds) || !length(out.rds) || is.na(out.rds[1]) || !nzchar(as.character(out.rds[1]))) {
    out.rds <- file.path(dirname(gtf.bgz), sprintf("gcanvas.tracks.b%s.rds", build))
  }
  out.rds <- .gcanvas_pretty_path(abs_path(out.rds))
  if (file.exists(out.rds) && !isTRUE(overwrite)) {
    stop("Output already exists: ", out.rds, ". Set overwrite=TRUE to replace.", call. = FALSE)
  }
  dir.create(dirname(out.rds), recursive = TRUE, showWarnings = FALSE)

  tabix <- .gcanvas_tabix_bin()
  seqnames <- system2(tabix, args = c("-l", gtf.bgz), stdout = TRUE, stderr = TRUE)
  seqnames <- unique(seqnames[!is.na(seqnames) & nzchar(seqnames)])
  if (!length(seqnames)) stop("No sequence names found via tabix -l: ", gtf.bgz, call. = FALSE)

  .gcanvas_note("gcanvas::gtf2rds", sprintf("reading %d sequence(s) from %s", length(seqnames), gtf.bgz))

  parse_one_seq <- function(seqname) {
    lines <- .gcanvas_tabix_query_lines(gtf.bgz, seqname, 1L, 2147483647L)
    if (!length(lines)) {
      return(list(
        gene = data.table::data.table(
          CHR = character(), gene_name = character(), gene_id = character(),
          biotype = character(), biotype_raw = character(), strand = character(),
          start = numeric(), end = numeric()
        ),
        exon = data.table::data.table(
          CHR = character(), gene_name = character(), gene_id = character(),
          transcript_id = character(), exon_number = character(),
          strand = character(), start = numeric(), end = numeric()
        )
      ))
    }

    dt <- data.table::fread(text = lines, sep = "\t", header = FALSE, data.table = TRUE, quote = "")
    if (ncol(dt) < 9) {
      return(list(
        gene = data.table::data.table(
          CHR = character(), gene_name = character(), gene_id = character(),
          biotype = character(), biotype_raw = character(), strand = character(),
          start = numeric(), end = numeric()
        ),
        exon = data.table::data.table(
          CHR = character(), gene_name = character(), gene_id = character(),
          transcript_id = character(), exon_number = character(),
          strand = character(), start = numeric(), end = numeric()
        )
      ))
    }

    data.table::setnames(dt, 1:9, c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"))
    dt <- dt[feature %in% c("gene", "exon")]
    if (!nrow(dt)) {
      return(list(
        gene = data.table::data.table(
          CHR = character(), gene_name = character(), gene_id = character(),
          biotype = character(), biotype_raw = character(), strand = character(),
          start = numeric(), end = numeric()
        ),
        exon = data.table::data.table(
          CHR = character(), gene_name = character(), gene_id = character(),
          transcript_id = character(), exon_number = character(),
          strand = character(), start = numeric(), end = numeric()
        )
      ))
    }

    gene_name <- .gcanvas_parse_attr_one(dt$attribute, "gene_name")
    gene_id <- .gcanvas_parse_attr_one(dt$attribute, "gene_id")
    gene_name <- .gcanvas_normalize_ensembl_gene_id(gene_name)
    gene_id <- .gcanvas_normalize_ensembl_gene_id(gene_id)
    gene_name[is.na(gene_name) | !nzchar(gene_name)] <- gene_id[is.na(gene_name) | !nzchar(gene_name)]

    biotype_raw <- .gcanvas_parse_attr_one(dt$attribute, "gene_biotype")
    biotype2 <- .gcanvas_parse_attr_one(dt$attribute, "gene_type")
    biotype_raw[is.na(biotype_raw)] <- biotype2[is.na(biotype_raw)]

    tx_id <- .gcanvas_parse_attr_one(dt$attribute, "transcript_id")
    ex_no <- .gcanvas_parse_attr_one(dt$attribute, "exon_number")
    chr_std <- normalize.chrom(dt$seqname)

    dt[, `:=`(
      CHR = chr_std,
      start = suppressWarnings(as.numeric(start)),
      end = suppressWarnings(as.numeric(end)),
      gene_name = gene_name,
      gene_id = gene_id,
      biotype_raw = biotype_raw,
      biotype = .gcanvas_biotype_group_vec(biotype_raw),
      transcript_id = tx_id,
      exon_number = ex_no
    )]
    dt <- dt[is.finite(start) & is.finite(end) & end >= start]
    dt <- dt[!is.na(CHR) & nzchar(CHR)]

    gene_dt <- unique(
      dt[feature == "gene" & (!is.na(gene_name) | !is.na(gene_id)),
         .(CHR, gene_name, gene_id, biotype, biotype_raw, strand, start, end)],
      by = c("CHR", "gene_name", "gene_id", "strand", "start", "end")
    )
    exon_dt <- unique(
      dt[feature == "exon" & (!is.na(gene_name) | !is.na(gene_id)),
         .(CHR, gene_name, gene_id, transcript_id, exon_number, strand, start, end)],
      by = c("CHR", "gene_name", "gene_id", "transcript_id", "exon_number", "strand", "start", "end")
    )
    list(gene = gene_dt, exon = exon_dt)
  }

  parsed <- lapply(seqnames, parse_one_seq)
  gene_all <- data.table::rbindlist(lapply(parsed, `[[`, "gene"), use.names = TRUE, fill = TRUE)
  exon_all <- data.table::rbindlist(lapply(parsed, `[[`, "exon"), use.names = TRUE, fill = TRUE)

  if (nrow(gene_all)) {
    gene_all[, gene_name := .gcanvas_as_char_no_null(gene_name, empty = "")]
    gene_all[, gene_id := .gcanvas_as_char_no_null(gene_id, empty = "")]
    gene_all <- gene_all[nzchar(gene_name) | nzchar(gene_id)]
    gene_all[!nzchar(gene_name) & nzchar(gene_id), gene_name := gene_id]
    gene_all <- unique(gene_all, by = c("CHR", "gene_name", "gene_id", "strand", "start", "end"))
    gene_all[, chr_order := rank.chrom(CHR)]
    data.table::setorderv(gene_all, c("chr_order", "CHR", "start", "end"), c(1L, 1L, 1L, 1L), na.last = TRUE)
    gene_all[, chr_order := NULL]
    data.table::setkey(gene_all, CHR, start, end)
  }
  if (nrow(exon_all)) {
    exon_all[, gene_name := .gcanvas_as_char_no_null(gene_name, empty = "")]
    exon_all[, gene_id := .gcanvas_as_char_no_null(gene_id, empty = "")]
    exon_all <- exon_all[nzchar(gene_name) | nzchar(gene_id)]
    exon_all[!nzchar(gene_name) & nzchar(gene_id), gene_name := gene_id]
    exon_all <- unique(exon_all, by = c("CHR", "gene_name", "gene_id", "transcript_id", "exon_number", "strand", "start", "end"))
    exon_all[, chr_order := rank.chrom(CHR)]
    data.table::setorderv(exon_all, c("chr_order", "CHR", "start", "end"), c(1L, 1L, 1L, 1L), na.last = TRUE)
    exon_all[, chr_order := NULL]
    data.table::setkey(exon_all, CHR, start, end)
  }

  out <- list(
    gene = gene_all,
    exon = exon_all,
    meta = list(
      source = gtf.bgz,
      build = build,
      created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      n_seq = as_int(length(seqnames)),
      n_gene = as_int(nrow(gene_all)),
      n_exon = as_int(nrow(exon_all))
    )
  )
  class(out) <- c("gcanvas_tracks_rds", class(out))

  saveRDS(out, file = out.rds, compress = compress)
  .gcanvas_note("gcanvas::gtf2rds", sprintf("saved %s (gene=%d, exon=%d)", out.rds, nrow(gene_all), nrow(exon_all)))

  invisible(list(
    file = out.rds,
    build = build,
    n_seq = as_int(length(seqnames)),
    n_gene = as_int(nrow(gene_all)),
    n_exon = as_int(nrow(exon_all))
  ))
}

# ----- GTF: biotype grouping -----
.gcanvas_biotype_group1 <- function(x) {
  x0 <- as.character(x)
  x <- tolower(x0)

  if (!is.na(x) && x == "protein_coding") return("protein_coding")

  if (!is.na(x) && (x == "mirna" || x == "microrna")) return("miRNA")

  lnc_set <- c(
    "lncrna", "lincrna",
    "antisense", "sense_intronic", "sense_overlapping", "processed_transcript",
    "3prime_overlapping_ncrna", "bidirectional_promoter_lncrna", "macro_lncrna",
    "non_coding", "noncoding"
  )
  if (!is.na(x) && x %in% lnc_set) return("lncRNA")

  "other"
}

.gcanvas_biotype_group_vec <- function(x) vapply(x, .gcanvas_biotype_group1, character(1))

.gcanvas_pick_gene_representative <- function(dt, key_col) {
  require_pkg("data.table")
  x <- data.table::as.data.table(data.table::copy(dt))
  key_col <- as.character(key_col)[1]
  if (is.na(key_col) || !nzchar(key_col) || !(key_col %in% names(x))) {
    stop("key_col must be a column in dt.", call. = FALSE)
  }
  if (!nrow(x)) return(x)

  if (!("biotype" %in% names(x))) x[, biotype := ""]
  if (!("start" %in% names(x))) x[, start := NA_real_]
  if (!("end" %in% names(x))) x[, end := NA_real_]

  x[, (key_col) := .gcanvas_as_char_no_null(get(key_col), empty = "")]
  x <- x[nzchar(get(key_col))]
  if (!nrow(x)) return(x)

  x[, start := suppressWarnings(as.numeric(start))]
  x[, end := suppressWarnings(as.numeric(end))]
  x[, width := data.table::fifelse(is.finite(start) & is.finite(end), end - start, -Inf)]
  x[, biotype_rank := data.table::fifelse(.gcanvas_as_char_no_null(biotype, empty = "") == "protein_coding", 0L, 1L)]
  data.table::setorderv(x, c(key_col, "biotype_rank", "width", "start", "end"), c(1L, 1L, -1L, 1L, 1L), na.last = TRUE)
  x <- unique(x, by = key_col)
  x[, c("width", "biotype_rank") := NULL]
  x
}

# ----- GTF: query gene/exon -----
.gcanvas_parse_attr_one <- function(attr, key) {
  x <- as.character(attr)
  out <- rep(NA_character_, length(x))

  ok <- !is.na(x) & nzchar(x)
  if (!any(ok)) return(out)

  pat <- sprintf('%s "([^"]+)"', key)
  m <- regexec(pat, x[ok], perl = TRUE)
  mm <- regmatches(x[ok], m)

  out_ok <- vapply(mm, function(z) {
    if (length(z) >= 2L) z[2] else NA_character_
  }, character(1))

  out[ok] <- out_ok
  out
}

.gcanvas_find_upward_file <- function(filename, start = getwd(), max_depth = 5L) {
  filename <- as.character(filename)[1]
  if (is.na(filename) || !nzchar(filename)) return(NULL)
  dir0 <- abs_path(start)
  max_depth <- as_int(max_depth)
  if (is.na(max_depth) || max_depth < 0L) max_depth <- 0L
  for (i in 0:max_depth) {
    cand <- file.path(dir0, filename)
    if (file.exists(cand)) return(abs_path(cand))
    par0 <- dirname(dir0)
    if (identical(par0, dir0)) break
    dir0 <- par0
  }
  NULL
}

.gcanvas_default_tracks_rds <- function(build = 38L) {
  build <- as_int(build)[1]
  if (is.na(build)) build <- 38L
  if (build == 19L) build <- 37L
  if (!(build %in% c(37L, 38L))) stop("build must be 37 or 38.", call. = FALSE)
  fname <- sprintf("gcanvas.tracks.b%s.rds", build)
  path0 <- .gcanvas_find_upward_file(fname, start = getwd(), max_depth = 6L)
  if (!is.null(path0) && nzchar(path0)) return(path0)
  stop(sprintf("Default track RDS not found for build %s: %s", build, fname), call. = FALSE)
}

.geneinfo_load_reference <- function(build = 38L, gtf = NULL) {
  require_pkg("data.table")
  rds_path <- NULL
  if (!is.null(gtf) && length(gtf) && !is.na(gtf[1]) && nzchar(as.character(gtf[1]))) {
    gtf0 <- abs_path(as.character(gtf[1]))
    if (!file.exists(gtf0)) stop("gtf not found: ", gtf0, call. = FALSE)
    if (grepl("\\.rds$", gtf0, ignore.case = TRUE)) {
      rds_path <- gtf0
    } else {
      require_pkg("digest")
      fi <- file.info(gtf0)
      key <- paste(
        gtf0,
        sprintf("%.0f", as.numeric(fi$size)[1]),
        sprintf("%.0f", as.numeric(fi$mtime)[1]),
        sprintf("%.0f", as.numeric(fi$ctime)[1]),
        sep = "|"
      )
      cache_dir <- .gcanvas_default_cache_dir(scope = "geneinfo", anchor = gtf0)
      dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
      stamp <- digest::digest(key, algo = "xxhash64", serialize = FALSE)
      bgz_path <- file.path(cache_dir, paste0(stamp, ".gtf.bgz"))
      rds_path <- file.path(cache_dir, paste0(stamp, ".rds"))
      if (!file.exists(rds_path)) {
        gtf_bgz <- if (grepl("\\.bgz$", gtf0, ignore.case = TRUE)) {
          gtf0
        } else {
          gtf_prepare_tabix(
            gtf0,
            out_bgz = bgz_path,
            overwrite = FALSE,
            sort = "auto",
            chr_order = "natural"
          )
        }
        gtf2rds(gtf.bgz = gtf_bgz, build = build, out.rds = rds_path, overwrite = TRUE)
      }
    }
  } else {
    rds_path <- .gcanvas_default_tracks_rds(build = build)
  }
  obj <- readRDS(rds_path)
  if (!is.list(obj) || is.null(obj$gene)) {
    stop("Track RDS must contain $gene table.", call. = FALSE)
  }
  gene_ref <- data.table::as.data.table(obj$gene)
  req <- c("CHR", "gene_name", "gene_id", "start", "end", "strand")
  miss <- setdiff(req, names(gene_ref))
  if (length(miss)) stop("Track RDS gene table missing columns: ", paste(miss, collapse = ", "), call. = FALSE)
  biotype_col <- if ("biotype_raw" %in% names(gene_ref)) "biotype_raw" else if ("biotype" %in% names(gene_ref)) "biotype" else NULL
  if (is.null(biotype_col)) stop("Track RDS gene table missing biotype/biotype_raw column.", call. = FALSE)

  gene_ref <- gene_ref[, .(
    gene_name = .gcanvas_normalize_ensembl_gene_id(.gcanvas_as_char_no_null(gene_name, empty = "")),
    gene_id = .gcanvas_normalize_ensembl_gene_id(.gcanvas_as_char_no_null(gene_id, empty = "")),
    chr = normalize.chrom(CHR),
    start = suppressWarnings(as.numeric(start)),
    end = suppressWarnings(as.numeric(end)),
    strand = .gcanvas_as_char_no_null(strand, empty = ""),
    biotype = .gcanvas_as_char_no_null(get(biotype_col), empty = "")
  )]
  gene_ref <- gene_ref[is.finite(start) & is.finite(end) & end >= start]
  gene_ref <- gene_ref[nzchar(gene_name) | nzchar(gene_id)]
  gene_ref[!nzchar(gene_name) & nzchar(gene_id), gene_name := gene_id]
  gene_ref[!nzchar(gene_id) & nzchar(gene_name) & .gcanvas_is_ensembl_gene_id(gene_name), gene_id := gene_name]
  sym_ref <- gene_ref[nzchar(gene_name)]
  if (nrow(sym_ref)) {
    sym_ref <- .gcanvas_pick_gene_representative(sym_ref, "gene_name")
    data.table::setkey(sym_ref, gene_name)
  }

  id_ref <- gene_ref[nzchar(gene_id)]
  if (nrow(id_ref)) {
    id_ref <- .gcanvas_pick_gene_representative(id_ref, "gene_id")
    data.table::setkey(id_ref, gene_id)
  }

  list(
    path = rds_path,
    gene = gene_ref,
    symbol = sym_ref,
    ensg = id_ref
  )
}

.geneinfo_as_type <- function(type) {
  if (is.list(type)) return(type)
  type0 <- tolower(as.character(type))
  type0[!nzchar(type0) | is.na(type0)] <- "auto"
  type0[type0 %in% c("ensgid", "ensembl", "ensembl_gene_id")] <- "ensg"
  if (!all(type0 %in% c("auto", "symbol", "ensg"))) {
    stop("type must be one of auto, symbol, ensg, ensgid, a character vector of those values, or list(symbol=..., ensg=...).", call. = FALSE)
  }
  if (length(type0) == 1L) return(type0)
  type0
}

.geneinfo_guess_symbol_col <- function(dt) {
  nms <- names(dt)
  cand <- c("gene", "Gene", "gene_name", "GENE", "symbol", "SYMBOL", "hgnc", "HGNC")
  hit <- cand[cand %in% nms]
  if (length(hit)) return(hit[1])
  NULL
}

.geneinfo_guess_ensg_col <- function(dt) {
  nms <- names(dt)
  cand <- c("ensg", "ENSG", "ensgid", "ENSGID", "gene_id", "GENE_ID", "ensembl_gene_id", "ENSEMBL_GENE_ID")
  hit <- cand[cand %in% nms]
  if (length(hit)) return(hit[1])
  NULL
}

.geneinfo_detect_two_cols <- function(dt) {
  nms <- names(dt)
  if (length(nms) != 2L) return(list(symbol = NULL, ensg = NULL))
  c1 <- nms[1]
  c2 <- nms[2]
  v1 <- .gcanvas_normalize_ensembl_gene_id(.gcanvas_as_char_no_null(dt[[c1]], empty = ""))
  v2 <- .gcanvas_normalize_ensembl_gene_id(.gcanvas_as_char_no_null(dt[[c2]], empty = ""))
  p1 <- mean(.gcanvas_is_ensembl_gene_id(v1) & nzchar(v1), na.rm = TRUE)
  p2 <- mean(.gcanvas_is_ensembl_gene_id(v2) & nzchar(v2), na.rm = TRUE)
  if (is.na(p1)) p1 <- 0
  if (is.na(p2)) p2 <- 0
  if (p1 > p2) return(list(symbol = c2, ensg = c1))
  if (p2 > p1) return(list(symbol = c1, ensg = c2))
  list(symbol = c1, ensg = c2)
}

.geneinfo_prepare_input <- function(gene, type = c("auto", "symbol", "ensg")) {
  require_pkg("data.table")
  type0 <- .geneinfo_as_type(type)
  .infer_symbol_from_vec <- function(v) {
    v0 <- .gcanvas_as_char_no_null(v, empty = "")
    out <- v0
    out[.gcanvas_is_ensembl_gene_id(out)] <- ""
    out
  }
  .infer_ensg_from_vec <- function(v) {
    v0 <- .gcanvas_as_char_no_null(v, empty = "")
    out <- .gcanvas_normalize_ensembl_gene_id(v0)
    out[!.gcanvas_is_ensembl_gene_id(v0)] <- NA_character_
    out
  }

  if (is.atomic(gene) && !is.data.frame(gene)) {
    out <- data.table::data.table(Gene = as.character(gene))
    if (identical(type0, "symbol")) {
      out[, `.geneinfo_symbol` := Gene]
      out[, `.geneinfo_ensg` := NA_character_]
    } else if (identical(type0, "ensg")) {
      out[, `.geneinfo_symbol` := NA_character_]
      out[, `.geneinfo_ensg` := .gcanvas_normalize_ensembl_gene_id(Gene)]
    } else {
      out[, `.geneinfo_ensg` := data.table::fifelse(.gcanvas_is_ensembl_gene_id(Gene), .gcanvas_normalize_ensembl_gene_id(Gene), NA_character_)]
      out[, `.geneinfo_symbol` := data.table::fifelse(is.na(.geneinfo_ensg), Gene, NA_character_)]
    }
    return(list(data = out, gene_col = "Gene"))
  }

  if (!(is.data.frame(gene) || data.table::is.data.table(gene))) {
    stop("gene must be a vector, data.frame, or data.table.", call. = FALSE)
  }

  dt <- data.table::as.data.table(data.table::copy(gene))
  if (!nrow(dt)) {
    if (!length(names(dt))) data.table::setnames(dt, "V1", "Gene")
    dt[, `.geneinfo_symbol` := character()]
    dt[, `.geneinfo_ensg` := character()]
    return(list(data = dt, gene_col = names(dt)[1]))
  }

  gene_col <- names(dt)[1]
  symbol_col <- NULL
  ensg_col <- NULL

  if (is.list(type0)) {
    symbol_col <- as.character(type0$symbol)[1] %||% NULL
    ensg_col <- as.character(type0$ensg)[1] %||% NULL
    if (!is.null(symbol_col) && !(symbol_col %in% names(dt))) stop("symbol column not found: ", symbol_col, call. = FALSE)
    if (!is.null(ensg_col) && !(ensg_col %in% names(dt))) stop("ensg column not found: ", ensg_col, call. = FALSE)
  } else if (length(type0) > 1L) {
    cols_use <- names(dt)[seq_len(min(length(type0), ncol(dt)))]
    for (i in seq_along(cols_use)) {
      ti <- type0[i]
      if (identical(ti, "symbol") && is.null(symbol_col)) symbol_col <- cols_use[i]
      if (identical(ti, "ensg") && is.null(ensg_col)) ensg_col <- cols_use[i]
    }
  } else if (length(names(dt)) == 2L) {
    symbol_col <- .geneinfo_guess_symbol_col(dt)
    ensg_col <- .geneinfo_guess_ensg_col(dt)
    if (is.null(symbol_col) || is.null(ensg_col)) {
      det <- .geneinfo_detect_two_cols(dt)
      if (is.null(symbol_col)) symbol_col <- det$symbol
      if (is.null(ensg_col)) ensg_col <- det$ensg
    }
  } else {
    symbol_col <- .geneinfo_guess_symbol_col(dt)
    ensg_col <- .geneinfo_guess_ensg_col(dt)
  }

  if (identical(type0, "symbol")) {
    if (is.null(symbol_col)) symbol_col <- gene_col
    dt[, `.geneinfo_symbol` := .gcanvas_as_char_no_null(get(symbol_col), empty = "")]
    dt[, `.geneinfo_ensg` := NA_character_]
    gene_col <- symbol_col
  } else if (identical(type0, "ensg")) {
    if (is.null(ensg_col)) ensg_col <- gene_col
    dt[, `.geneinfo_symbol` := NA_character_]
    dt[, `.geneinfo_ensg` := .gcanvas_normalize_ensembl_gene_id(get(ensg_col))]
    gene_col <- ensg_col
  } else {
    if (!is.null(symbol_col) && !is.null(ensg_col) && !identical(symbol_col, ensg_col)) {
      s0 <- .gcanvas_as_char_no_null(dt[[symbol_col]], empty = "")
      e0 <- .gcanvas_as_char_no_null(dt[[ensg_col]], empty = "")
      s1 <- .infer_symbol_from_vec(s0)
      e1 <- .infer_ensg_from_vec(e0)
      s2 <- .infer_symbol_from_vec(e0)
      e2 <- .infer_ensg_from_vec(s0)
      s1[!nzchar(s1) & nzchar(s2)] <- s2[!nzchar(s1) & nzchar(s2)]
      e1[(is.na(e1) | !nzchar(e1)) & !is.na(e2) & nzchar(e2)] <- e2[(is.na(e1) | !nzchar(e1)) & !is.na(e2) & nzchar(e2)]
      dt[, `.geneinfo_symbol` := s1]
      dt[, `.geneinfo_ensg` := e1]
      gene_col <- symbol_col
    } else {
      if (is.null(gene_col) || !nzchar(gene_col)) gene_col <- names(dt)[1]
      vals0 <- .gcanvas_as_char_no_null(get(gene_col), empty = "")
      dt[, `.geneinfo_ensg` := .infer_ensg_from_vec(vals0)]
      dt[, `.geneinfo_symbol` := .infer_symbol_from_vec(vals0)]
    }
  }

  list(data = dt, gene_col = gene_col, symbol_col = symbol_col, ensg_col = ensg_col)
}

.geneinfo_attach_map <- function(dt, ref) {
  require_pkg("data.table")
  out <- data.table::copy(dt)
  out[, `:=`(
    gene_name = NA_character_,
    ensg_id = NA_character_,
    chr = NA_character_,
    start = NA_real_,
    end = NA_real_,
    strand = NA_character_,
    biotype = NA_character_
  )]

  sym_query <- .gcanvas_as_char_no_null(out$.geneinfo_symbol, empty = "")
  ensg_query <- .gcanvas_normalize_ensembl_gene_id(out$.geneinfo_ensg)

  sym_map <- NULL
  ensg_map <- NULL

  if (nrow(ref$symbol) && any(nzchar(sym_query))) {
    sym_map <- ref$symbol[.(sym_query)]
  }
  if (nrow(ref$ensg) && any(!is.na(ensg_query) & nzchar(ensg_query))) {
    ensg_map <- ref$ensg[.(ensg_query)]
  }

  .coalesce_chr <- function(a, b) {
    a0 <- .gcanvas_as_char_no_null(a, empty = "")
    b0 <- .gcanvas_as_char_no_null(b, empty = "")
    data.table::fifelse(nzchar(a0), a0, b0)
  }
  .coalesce_num <- function(a, b) {
    a0 <- suppressWarnings(as.numeric(a))
    b0 <- suppressWarnings(as.numeric(b))
    out0 <- a0
    miss <- !is.finite(out0) | is.na(out0)
    out0[miss] <- b0[miss]
    out0
  }

  if (!is.null(sym_map)) {
    out[, `:=`(
      gene_name = sym_map$gene_name,
      ensg_id = sym_map$gene_id,
      chr = sym_map$chr,
      start = sym_map$start,
      end = sym_map$end,
      strand = sym_map$strand,
      biotype = sym_map$biotype
    )]
  }
  if (!is.null(ensg_map)) {
    out[, `:=`(
      gene_name = .coalesce_chr(gene_name, ensg_map$gene_name),
      ensg_id = .coalesce_chr(ensg_id, ensg_map$gene_id),
      chr = .coalesce_chr(chr, ensg_map$chr),
      start = .coalesce_num(start, ensg_map$start),
      end = .coalesce_num(end, ensg_map$end),
      strand = .coalesce_chr(strand, ensg_map$strand),
      biotype = .coalesce_chr(biotype, ensg_map$biotype)
    )]
    out[!nzchar(gene_name), gene_name := NA_character_]
    out[!nzchar(ensg_id), ensg_id := NA_character_]
    out[!nzchar(chr), chr := NA_character_]
    out[!nzchar(strand), strand := NA_character_]
    out[!nzchar(biotype), biotype := NA_character_]
  }

  out[]
}

.geneinfo_fill_possible <- function(dt) {
  require_pkg("data.table")
  out <- data.table::copy(dt)
  anno_chr <- c("gene_name", "ensg_id", "chr", "strand", "biotype")
  anno_num <- c("start", "end")

  .is_source <- function(d) {
    (!is.na(d$ensg_id) & nzchar(d$ensg_id)) |
      (!is.na(d$chr) & nzchar(d$chr)) |
      is.finite(d$start) |
      is.finite(d$end)
  }
  .needs_fill <- function(d) {
    is.na(d$gene_name) | !nzchar(d$gene_name) |
      is.na(d$ensg_id) | !nzchar(d$ensg_id) |
      is.na(d$chr) | !nzchar(d$chr) |
      !is.finite(d$start) |
      !is.finite(d$end) |
      is.na(d$strand) | !nzchar(d$strand) |
      is.na(d$biotype) | !nzchar(d$biotype)
  }
  .fill_from_key <- function(d, key_col) {
    key <- .gcanvas_as_char_no_null(d[[key_col]], empty = "")
    ok_key <- nzchar(key)
    source_idx <- .is_source(d) & ok_key
    if (!any(source_idx)) return(d)

    ref_dt <- data.table::copy(d[source_idx, c(key_col, anno_chr, anno_num), with = FALSE])
    ref_dt[, .completeness := rowSums(cbind(
      !is.na(gene_name) & nzchar(gene_name),
      !is.na(ensg_id) & nzchar(ensg_id),
      !is.na(chr) & nzchar(chr),
      is.finite(start),
      is.finite(end),
      !is.na(strand) & nzchar(strand),
      !is.na(biotype) & nzchar(biotype)
    ))]
    data.table::setorderv(ref_dt, c(key_col, ".completeness"), c(1L, -1L), na.last = TRUE)
    ref_dt <- unique(ref_dt, by = key_col)
    data.table::setkeyv(ref_dt, key_col)

    idx <- which(.needs_fill(d) & ok_key)
    if (!length(idx)) return(d)
    m <- ref_dt[.(key[idx])]
    for (nm in anno_chr) {
      cur <- .gcanvas_as_char_no_null(d[[nm]], empty = "")
      add <- .gcanvas_as_char_no_null(m[[nm]], empty = "")
      fill <- idx[!nzchar(cur[idx]) & nzchar(add)]
      if (length(fill)) data.table::set(d, i = fill, j = nm, value = add[!nzchar(cur[idx]) & nzchar(add)])
    }
    for (nm in anno_num) {
      cur <- suppressWarnings(as.numeric(d[[nm]]))
      add <- suppressWarnings(as.numeric(m[[nm]]))
      fill_mask <- (!is.finite(cur[idx]) | is.na(cur[idx])) & is.finite(add)
      fill <- idx[fill_mask]
      if (length(fill)) data.table::set(d, i = fill, j = nm, value = add[fill_mask])
    }
    d
  }

  out <- .fill_from_key(out, ".geneinfo_symbol")
  out <- .fill_from_key(out, ".geneinfo_ensg")
  out[]
}

.geneinfo_prepare_update <- function(info.update) {
  require_pkg("data.table")
  if (is.null(info.update)) return(NULL)
  if (!(is.data.frame(info.update) || data.table::is.data.table(info.update))) {
    stop("info.update must be a data.frame/data.table when provided.", call. = FALSE)
  }
  dt <- data.table::as.data.table(data.table::copy(info.update))
  if (!nrow(dt) || !ncol(dt)) return(NULL)

  canon <- function(x) tolower(gsub("[^a-z0-9]+", "", as.character(x)))
  pick_col <- function(targets) {
    nms <- names(dt)
    cn <- canon(nms)
    tg <- canon(targets)
    hit <- match(tg, cn, nomatch = 0L)
    hit <- hit[hit > 0L]
    if (length(hit)) return(nms[hit[1]])
    hit <- vapply(cn, function(z) any(vapply(tg, function(t) nzchar(t) && (grepl(t, z, fixed = TRUE) || grepl(z, t, fixed = TRUE)), logical(1))), logical(1))
    if (any(hit)) return(nms[which(hit)[1]])
    NULL
  }

  symbol_col <- .geneinfo_guess_symbol_col(dt)
  if (is.null(symbol_col)) symbol_col <- pick_col(c("Gene", "gene", "gene_name", "symbol"))
  ensg_col <- .geneinfo_guess_ensg_col(dt)
  if (is.null(ensg_col)) ensg_col <- pick_col(c("ensg_id", "gene_id", "ENSGID", "ensg"))

  anno_map <- list(
    gene_name = pick_col(c("gene_name", "gene")),
    ensg_id = pick_col(c("ensg_id", "gene_id", "ENSGID", "ensg")),
    chr = pick_col(c("chr", "chrom", "chromosome")),
    start = pick_col(c("start", "pos_start", "tx_start")),
    end = pick_col(c("end", "pos_end", "tx_end")),
    strand = pick_col(c("strand")),
    biotype = pick_col(c("biotype", "gene_biotype", "gene_type"))
  )
  anno_map <- anno_map[!vapply(anno_map, is.null, logical(1))]

  if (is.null(symbol_col) && is.null(ensg_col) && !length(anno_map)) return(NULL)

  out <- data.table::data.table(
    .info_symbol = if (!is.null(symbol_col)) .gcanvas_as_char_no_null(dt[[symbol_col]], empty = "") else rep("", nrow(dt)),
    .info_ensg = if (!is.null(ensg_col)) .gcanvas_normalize_ensembl_gene_id(dt[[ensg_col]]) else rep(NA_character_, nrow(dt))
  )
  for (nm in names(anno_map)) {
    src <- anno_map[[nm]]
    if (nm %in% c("start", "end")) {
      out[[nm]] <- suppressWarnings(as.numeric(dt[[src]]))
    } else if (nm == "chr") {
      out[[nm]] <- normalize.chrom(dt[[src]])
    } else if (nm == "ensg_id") {
      out[[nm]] <- .gcanvas_normalize_ensembl_gene_id(dt[[src]])
    } else {
      out[[nm]] <- .gcanvas_as_char_no_null(dt[[src]], empty = "")
    }
  }
  out[]
}

.geneinfo_apply_update <- function(dt, info.update, replace = FALSE) {
  require_pkg("data.table")
  info0 <- .geneinfo_prepare_update(info.update)
  if (is.null(info0) || !nrow(info0)) return(dt)

  out <- data.table::copy(dt)
  anno_chr <- intersect(c("gene_name", "ensg_id", "chr", "strand", "biotype"), names(info0))
  anno_num <- intersect(c("start", "end"), names(info0))
  if (!length(anno_chr) && !length(anno_num)) return(out)

  rank_info <- function(x, key_col) {
    ref_dt <- data.table::copy(x)
    ref_dt <- ref_dt[
      (!is.na(get(key_col)) & nzchar(get(key_col))) |
        ("start" %in% names(ref_dt) & is.finite(start)) |
        ("end" %in% names(ref_dt) & is.finite(end))
    ]
    if (!nrow(ref_dt)) return(ref_dt[0])
    chr_ok <- if ("chr" %in% names(ref_dt)) !is.na(ref_dt$chr) & nzchar(ref_dt$chr) else rep(FALSE, nrow(ref_dt))
    start_ok <- if ("start" %in% names(ref_dt)) is.finite(ref_dt$start) else rep(FALSE, nrow(ref_dt))
    end_ok <- if ("end" %in% names(ref_dt)) is.finite(ref_dt$end) else rep(FALSE, nrow(ref_dt))
    ref_dt[, .completeness := rowSums(cbind(
      if ("gene_name" %in% names(ref_dt)) !is.na(gene_name) & nzchar(gene_name) else rep(FALSE, .N),
      if ("ensg_id" %in% names(ref_dt)) !is.na(ensg_id) & nzchar(ensg_id) else rep(FALSE, .N),
      chr_ok, start_ok, end_ok,
      if ("strand" %in% names(ref_dt)) !is.na(strand) & nzchar(strand) else rep(FALSE, .N),
      if ("biotype" %in% names(ref_dt)) !is.na(biotype) & nzchar(biotype) else rep(FALSE, .N)
    ))]
    data.table::setorderv(ref_dt, c(key_col, ".completeness"), c(1L, -1L), na.last = TRUE)
    ref_dt <- unique(ref_dt, by = key_col)
    data.table::setkeyv(ref_dt, key_col)
    ref_dt
  }

  sym_ref <- rank_info(info0[nzchar(.info_symbol)], ".info_symbol")
  ensg_ref <- rank_info(info0[!is.na(.info_ensg) & nzchar(.info_ensg)], ".info_ensg")

  coalesce_chr <- function(cur, add, do_replace = FALSE) {
    cur0 <- .gcanvas_as_char_no_null(cur, empty = "")
    add0 <- .gcanvas_as_char_no_null(add, empty = "")
    if (isTRUE(do_replace)) {
      out0 <- cur0
      idx <- nzchar(add0)
      out0[idx] <- add0[idx]
      return(out0)
    }
    data.table::fifelse(nzchar(cur0), cur0, add0)
  }
  coalesce_num <- function(cur, add, do_replace = FALSE) {
    cur0 <- suppressWarnings(as.numeric(cur))
    add0 <- suppressWarnings(as.numeric(add))
    if (isTRUE(do_replace)) {
      idx <- is.finite(add0) & !is.na(add0)
      cur0[idx] <- add0[idx]
      return(cur0)
    }
    miss <- !is.finite(cur0) | is.na(cur0)
    cur0[miss] <- add0[miss]
    cur0
  }

  if (nrow(sym_ref) && any(nzchar(.gcanvas_as_char_no_null(out$.geneinfo_symbol, empty = "")))) {
    sym_map <- sym_ref[.(.gcanvas_as_char_no_null(out$.geneinfo_symbol, empty = ""))]
    for (nm in anno_chr) out[[nm]] <- coalesce_chr(out[[nm]], sym_map[[nm]], do_replace = replace)
    for (nm in anno_num) out[[nm]] <- coalesce_num(out[[nm]], sym_map[[nm]], do_replace = replace)
  }
  if (nrow(ensg_ref) && any(!is.na(out$.geneinfo_ensg) & nzchar(out$.geneinfo_ensg))) {
    ensg_map <- ensg_ref[.(.gcanvas_normalize_ensembl_gene_id(out$.geneinfo_ensg))]
    for (nm in anno_chr) out[[nm]] <- coalesce_chr(out[[nm]], ensg_map[[nm]], do_replace = replace)
    for (nm in anno_num) out[[nm]] <- coalesce_num(out[[nm]], ensg_map[[nm]], do_replace = replace)
  }
  out[]
}

#' Look up gene metadata by symbol or Ensembl id
#'
#' Returns gene-level annotation (chromosome, start/end, strand, biotype)
#' for the requested gene(s) by symbol or Ensembl gene id, using the bundled
#' or user-supplied GTF cache.
#'
#' @param gene Character vector of gene symbols or Ensembl gene ids, or a
#'   `data.frame`/`data.table` containing one such column.
#' @param build Genome build (`37` or `38`).
#' @param gtf Optional path to a `gtf2rds()`-produced cache to use instead of
#'   the package-bundled tracks.
#' @param type Input type detector: `"auto"`, `"symbol"`, or `"ensg"`.
#' @param all.possible Logical. Return every matching annotation row instead
#'   of one row per query.
#' @param info.update Optional named list / table of overrides to merge in.
#' @param replace Logical. If `TRUE`, overrides replace existing annotation
#'   values; otherwise they only fill missing ones.
#'
#' @return A `data.table` with one row per resolved annotation.
#' @export
geneinfo <- function(gene,
                     build = 38L,
                     gtf = NULL,
                     type = c("auto", "symbol", "ensg"),
                     all.possible = TRUE,
                     info.update = NULL,
                     replace = FALSE) {
  require_pkg("data.table")
  build <- as_int(build)[1]
  if (is.na(build)) build <- 38L
  if (build == 19L) build <- 37L
  if (!(build %in% c(37L, 38L))) stop("build must be 37 or 38.", call. = FALSE)

  prep <- .geneinfo_prepare_input(gene, type = type)
  anno_cols <- c("gene_name", "ensg_id", "chr", "start", "end", "strand", "biotype")
  input_cols <- NULL
  rename_map <- character()
  if (is.data.frame(gene) || data.table::is.data.table(gene)) {
    input_cols <- setdiff(names(prep$data), c(".geneinfo_symbol", ".geneinfo_ensg"))
    conflicts <- intersect(input_cols, anno_cols)
    if (length(conflicts)) {
      tmp_names <- paste0(".geneinfo_input_", conflicts)
      data.table::setnames(prep$data, conflicts, tmp_names)
      rename_map <- stats::setNames(conflicts, tmp_names)
    }
  }
  ref <- .geneinfo_load_reference(build = build, gtf = gtf)
  out <- .geneinfo_attach_map(prep$data, ref = ref)
  if (isTRUE(all.possible)) {
    out <- .geneinfo_fill_possible(out)
  }
  if (!is.null(info.update)) {
    out <- .geneinfo_apply_update(out, info.update = info.update, replace = replace)
  }

  if (length(rename_map)) {
    rename_anno <- stats::setNames(paste0(anno_cols, "_gtf"), anno_cols)
    rename_anno <- rename_anno[names(rename_anno) %in% names(out) & names(rename_anno) %in% unname(rename_map)]
    if (length(rename_anno)) {
      data.table::setnames(out, names(rename_anno), unname(rename_anno))
    }
    data.table::setnames(out, names(rename_map), unname(rename_map))
    anno_cols_final <- data.table::fifelse(
      anno_cols %in% unname(rename_map),
      paste0(anno_cols, "_gtf"),
      anno_cols
    )
  } else {
    anno_cols_final <- anno_cols
  }

  if (!is.data.frame(gene) && !data.table::is.data.table(gene)) {
    data.table::setcolorder(out, c("Gene", "gene_name", "ensg_id", "chr", "start", "end", "strand", "biotype"))
  } else {
    data.table::setcolorder(out, unique(c(input_cols, anno_cols_final)))
  }
  int_cols <- intersect(c("start", "end", "start_gtf", "end_gtf"), names(out))
  for (nm in int_cols) {
    out[[nm]] <- suppressWarnings(as.integer(out[[nm]]))
  }
  out[, c(".geneinfo_symbol", ".geneinfo_ensg") := NULL]
  out[]
}

gtf_query_gene_exon <- function(gtf_bgz, chrom, start, end,
                                features = c("gene", "exon"),
                                keep_biotype = NULL) {
  require_pkg(c("Rsamtools", "GenomicRanges", "IRanges", "data.table"))
  empty_gene <- data.table::data.table(
    gene_name = character(), gene_id = character(), biotype = character(), biotype_raw = character(),
    strand = character(), start = numeric(), end = numeric()
  )
  empty_exon <- data.table::data.table(
    gene_name = character(), gene_id = character(), transcript_id = character(),
    exon_number = character(), strand = character(), start = numeric(), end = numeric()
  )

  chrom <- normalize.chrom(chrom)[1]
  start <- as_int(start); end <- as_int(end)
  if (is.na(start) || is.na(end) || end < start) stop("Invalid start/end.", call. = FALSE)
  if (start < 1L) start <- 1L

  gtf_path <- abs_path(gtf_bgz)

  seq_candidates <- unique(c(chrom, chr_to_ucsc1(chrom)))
  seq_candidates <- seq_candidates[!is.na(seq_candidates) & nzchar(seq_candidates)]

  lines <- character()
  for (seqname in seq_candidates) {
    got <- .gcanvas_tabix_query_lines(gtf_path, seqname, start, end)
    if (length(got)) {
      lines <- got
      break
    }
  }


  if (!length(lines)) {
    return(list(
      gene = data.table::copy(empty_gene),
      exon = data.table::copy(empty_exon),
      chrom = chrom
    ))
  }

  dt <- data.table::fread(text = lines, sep = "\t", header = FALSE, data.table = TRUE, quote = "")
  if (ncol(dt) < 9) {
    return(list(
      gene = data.table::copy(empty_gene),
      exon = data.table::copy(empty_exon),
      chrom = chrom
    ))
  }

  data.table::setnames(dt, 1:9, c("seqname","source","feature","start","end","score","strand","frame","attribute"))
  dt <- dt[feature %in% features]
  if (!nrow(dt)) return(list(gene = data.table::copy(empty_gene), exon = data.table::copy(empty_exon), chrom = chrom))

  gene_name <- .gcanvas_parse_attr_one(dt$attribute, "gene_name")
  gene_id <- .gcanvas_parse_attr_one(dt$attribute, "gene_id")
  gene_name <- .gcanvas_normalize_ensembl_gene_id(gene_name)
  gene_id <- .gcanvas_normalize_ensembl_gene_id(gene_id)

  biotype_raw <- .gcanvas_parse_attr_one(dt$attribute, "gene_biotype")
  biotype2 <- .gcanvas_parse_attr_one(dt$attribute, "gene_type")
  biotype_raw[is.na(biotype_raw)] <- biotype2[is.na(biotype_raw)]

  tx_id <- .gcanvas_parse_attr_one(dt$attribute, "transcript_id")
  ex_no <- .gcanvas_parse_attr_one(dt$attribute, "exon_number")

  gene_name[is.na(gene_name) | !nzchar(gene_name)] <- gene_id[is.na(gene_name) | !nzchar(gene_name)]

  dt[, `:=`(
    gene_name = gene_name,
    gene_id = gene_id,
    biotype_raw = biotype_raw,
    biotype = .gcanvas_biotype_group_vec(biotype_raw),
    transcript_id = tx_id,
    exon_number = ex_no,
    start = suppressWarnings(as.numeric(start)),
    end = suppressWarnings(as.numeric(end))
  )]

  if (!is.null(keep_biotype)) {
    keep_biotype <- unique(as.character(keep_biotype))
    dt <- dt[biotype %in% keep_biotype | is.na(biotype)]
  }

  gene_dt <- dt[feature == "gene"]
  if (nrow(gene_dt)) {
    gene_df <- gene_dt[, .(gene_name, gene_id, biotype, biotype_raw, strand, start, end)]
    gene_df[, .gene_key := data.table::fifelse(!is.na(gene_id) & nzchar(gene_id), gene_id, gene_name)]
    gene_df <- .gcanvas_pick_gene_representative(gene_df, ".gene_key")
    gene_df[, .gene_key := NULL]
  } else {
    exon_dt0 <- dt[feature == "exon" & !is.na(gene_name)]
    if (!nrow(exon_dt0)) {
      gene_df <- data.table::copy(empty_gene)
    } else {
      mn <- exon_dt0[, .(start = min(start, na.rm = TRUE)), by = .(gene_name, gene_id, strand)]
      mx <- exon_dt0[, .(end = max(end, na.rm = TRUE)), by = .(gene_name, gene_id, strand)]
      gene_df <- merge(mn, mx, by = c("gene_name","gene_id","strand"), all = TRUE)
      gene_df[, `:=`(biotype = NA_character_, biotype_raw = NA_character_)]
      data.table::setcolorder(gene_df, c("gene_name","gene_id","biotype","biotype_raw","strand","start","end"))
    }
  }

  exon_dt <- dt[feature == "exon" & !is.na(gene_name)]
  exon_df <- unique(exon_dt[, .(gene_name, gene_id, transcript_id, exon_number, strand, start, end)])
  list(gene = gene_df, exon = exon_df, chrom = chrom)
}

.gcanvas_normalize_added_gene_table <- function(x, default_biotype = "protein_coding") {
  require_pkg("data.table")
  empty_gene <- data.table::data.table(
    gene_name = character(), gene_id = character(),
    biotype = character(), biotype_raw = character(),
    strand = character(), start = numeric(), end = numeric()
  )
  if (is.null(x)) {
    return(data.table::copy(empty_gene))
  }
  dt <- data.table::as.data.table(x)
  if (!nrow(dt)) {
    return(data.table::copy(empty_gene))
  }

  nm <- names(dt)
  name_col <- nm[nm %in% c("gene_name", "label", "gene", "name")][1]
  if (is.na(name_col) || !nzchar(name_col)) {
    if ("gene_id" %in% nm) {
      dt[, gene_name := .gcanvas_as_char_no_null(gene_id, empty = "")]
    } else {
      stop("gene.add gene table must include one of: gene_name, label, gene, name.", call. = FALSE)
    }
  } else {
    dt[, gene_name := .gcanvas_as_char_no_null(get(name_col), empty = "")]
  }

  if ("gene_id" %in% nm) {
    dt[, gene_id := .gcanvas_as_char_no_null(gene_id, empty = "")]
  } else {
    dt[, gene_id := gene_name]
  }
  if ("strand" %in% nm) {
    dt[, strand := .gcanvas_as_char_no_null(strand, empty = "*")]
  } else {
    dt[, strand := "*"]
  }
  if ("biotype" %in% nm) {
    dt[, biotype := .gcanvas_as_char_no_null(biotype, empty = default_biotype)]
  } else {
    dt[, biotype := default_biotype]
  }
  if ("biotype_raw" %in% nm) {
    dt[, biotype_raw := .gcanvas_as_char_no_null(biotype_raw, empty = "")]
  } else {
    dt[, biotype_raw := biotype]
  }

  if (!("start" %in% nm) || !("end" %in% nm)) {
    stop("gene.add gene table must include numeric start and end columns.", call. = FALSE)
  }
  dt[, start := suppressWarnings(as.numeric(start))]
  dt[, end := suppressWarnings(as.numeric(end))]

  dt[, gene_name := .gcanvas_normalize_ensembl_gene_id(gene_name)]
  dt[, gene_id := .gcanvas_normalize_ensembl_gene_id(gene_id)]
  dt <- dt[!is.na(start) & !is.na(end) & is.finite(start) & is.finite(end) & end >= start]
  dt <- dt[nzchar(gene_name) | nzchar(gene_id)]
  if (!nrow(dt)) {
    return(data.table::copy(empty_gene))
  }
  dt[, .(gene_name, gene_id, biotype, biotype_raw, strand, start, end)]
}

.gcanvas_normalize_added_exon_table <- function(x, gene_table = NULL) {
  require_pkg("data.table")
  empty_exon <- data.table::data.table(
    gene_name = character(), gene_id = character(),
    transcript_id = character(), exon_number = character(),
    strand = character(), start = numeric(), end = numeric()
  )
  if (is.null(x)) {
    return(data.table::copy(empty_exon))
  }
  dt <- data.table::as.data.table(x)
  if (!nrow(dt)) {
    return(data.table::copy(empty_exon))
  }

  nm <- names(dt)
  name_col <- nm[nm %in% c("gene_name", "label", "gene", "name")][1]
  if (is.na(name_col) || !nzchar(name_col)) {
    if ("gene_id" %in% nm) {
      dt[, gene_name := .gcanvas_as_char_no_null(gene_id, empty = "")]
    } else {
      stop("gene.add exon table must include one of: gene_name, label, gene, name.", call. = FALSE)
    }
  } else {
    dt[, gene_name := .gcanvas_as_char_no_null(get(name_col), empty = "")]
  }
  if ("gene_id" %in% nm) {
    dt[, gene_id := .gcanvas_as_char_no_null(gene_id, empty = "")]
  } else {
    dt[, gene_id := gene_name]
  }
  if ("strand" %in% nm) {
    dt[, strand := .gcanvas_as_char_no_null(strand, empty = "*")]
  } else {
    dt[, strand := "*"]
  }
  if ("transcript_id" %in% nm) {
    dt[, transcript_id := .gcanvas_as_char_no_null(transcript_id, empty = "")]
  } else {
    dt[, transcript_id := ""]
  }
  if ("exon_number" %in% nm) {
    dt[, exon_number := .gcanvas_as_char_no_null(exon_number, empty = "")]
  } else {
    dt[, exon_number := ""]
  }
  if (!("start" %in% nm) || !("end" %in% nm)) {
    stop("gene.add exon table must include numeric start and end columns.", call. = FALSE)
  }
  dt[, start := suppressWarnings(as.numeric(start))]
  dt[, end := suppressWarnings(as.numeric(end))]
  dt[, gene_name := .gcanvas_normalize_ensembl_gene_id(gene_name)]
  dt[, gene_id := .gcanvas_normalize_ensembl_gene_id(gene_id)]
  dt <- dt[!is.na(start) & !is.na(end) & is.finite(start) & is.finite(end) & end >= start]
  dt <- dt[nzchar(gene_name) | nzchar(gene_id)]

  if (!is.null(gene_table) && nrow(gene_table) && nrow(dt)) {
    gdt <- data.table::as.data.table(gene_table)[, .(gene_name, gene_id)]
    gdt[, gene_name := .gcanvas_normalize_ensembl_gene_id(gene_name)]
    gdt[, gene_id := .gcanvas_normalize_ensembl_gene_id(gene_id)]
    dt <- merge(dt, unique(gdt), by = c("gene_name", "gene_id"), all.x = FALSE, all.y = FALSE)
  }

  dt[, .(gene_name, gene_id, transcript_id, exon_number, strand, start, end)]
}

.gcanvas_extract_gene_add <- function(gene_add) {
  empty_gene <- data.table::data.table(
    gene_name = character(), gene_id = character(),
    biotype = character(), biotype_raw = character(),
    strand = character(), start = numeric(), end = numeric()
  )
  empty_exon <- data.table::data.table(
    gene_name = character(), gene_id = character(),
    transcript_id = character(), exon_number = character(),
    strand = character(), start = numeric(), end = numeric()
  )
  if (is.null(gene_add)) return(list(gene = empty_gene, exon = empty_exon))

  if (is.data.frame(gene_add) || data.table::is.data.table(gene_add)) {
    g <- .gcanvas_normalize_added_gene_table(gene_add)
    ex <- empty_exon
    if (all(c("exon_start", "exon_end") %in% names(gene_add))) {
      ex0 <- data.table::as.data.table(gene_add)
      nm <- names(ex0)
      name_col <- nm[nm %in% c("gene_name", "label", "gene", "name")][1]
      if (is.na(name_col) || !nzchar(name_col)) name_col <- "gene_id"
      ex <- .gcanvas_normalize_added_exon_table(data.table::data.table(
        gene_name = ex0[[name_col]],
        gene_id = if ("gene_id" %in% nm) ex0$gene_id else ex0[[name_col]],
        strand = if ("strand" %in% nm) ex0$strand else "*",
        start = ex0$exon_start,
        end = ex0$exon_end
      ), gene_table = g)
    }
    return(list(gene = g, exon = ex))
  }

  if (is.list(gene_add)) {
    g <- .gcanvas_normalize_added_gene_table(gene_add$gene %||% gene_add$genes)
    ex <- .gcanvas_normalize_added_exon_table(gene_add$exon %||% gene_add$exons, gene_table = g)
    return(list(gene = g, exon = ex))
  }

  stop("gene.add must be data.frame/data.table or list(gene=..., exon=...).", call. = FALSE)
}

