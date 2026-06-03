# Toy/synthetic data generators used for examples and quick demos.
# Kept unexported in v0.1.0 (will graduate to exported once tests / vignettes
# are in place).

.gcanvas_toy_sample_alleles <- function(n) {
  n <- as_int(n)
  if (!is.finite(n) || n <= 0L) {
    return(data.table::data.table(A1 = character(), A2 = character()))
  }
  base <- c("A", "C", "G", "T")
  a1 <- character(n)
  a2 <- character(n)
  for (i in seq_len(n)) {
    aa <- sort(sample(base, size = 2L, replace = FALSE))
    a1[i] <- aa[1]
    a2[i] <- aa[2]
  }
  data.table::data.table(A1 = a1, A2 = a2)
}

.gcanvas_toy_sim_haplotype <- function(n, maf, copy_prob) {
  n <- as_int(n)
  maf <- suppressWarnings(as.numeric(maf))
  m <- length(maf)
  if (!is.finite(n) || n <= 0L || m <= 0L) {
    return(matrix(0L, nrow = max(1L, n), ncol = max(1L, m)))
  }
  h <- matrix(0L, nrow = n, ncol = m)
  h[, 1] <- stats::rbinom(n, size = 1L, prob = maf[1])
  if (m >= 2L) {
    for (j in 2:m) {
      cp <- max(0.02, min(0.995, as.numeric(copy_prob[j])))
      copied <- stats::runif(n) < cp
      draw <- stats::rbinom(n, size = 1L, prob = maf[j])
      hj <- ifelse(copied, h[, j - 1L], draw)
      flip <- stats::runif(n) < 0.003
      if (any(flip)) hj[flip] <- 1L - hj[flip]
      h[, j] <- as_int(hj)
    }
  }
  h
}

.gcanvas_toy_parse_bfile_option <- function(bfile, prefix_pattern = "toy_gwas_") {
  bfile_prefix <- NULL
  b0 <- if (is.null(bfile) || length(bfile) == 0L) TRUE else bfile[1]
  write_bfile <- TRUE
  if (is.logical(b0)) {
    write_bfile <- isTRUE(b0)
    if (isTRUE(write_bfile)) bfile_prefix <- tempfile(pattern = prefix_pattern)
  } else {
    btxt <- trimws(as.character(b0))
    blow <- tolower(btxt)
    if (is.na(btxt) || !nzchar(btxt) || blow %in% c("true", "t", "yes", "y", "1", "auto")) {
      write_bfile <- TRUE
      bfile_prefix <- tempfile(pattern = prefix_pattern)
    } else if (blow %in% c("false", "f", "no", "n", "0", "none", "null", "na")) {
      write_bfile <- FALSE
    } else {
      write_bfile <- TRUE
      bfile_prefix <- .gcanvas_normalize_bfile_prefix(btxt)
      if (is.null(bfile_prefix)) bfile_prefix <- abs_path(btxt)
    }
  }
  list(write_bfile = isTRUE(write_bfile), bfile_prefix = bfile_prefix)
}

.gcanvas_toy_split_samples <- function(n.sample, n.ancestry) {
  n.ancestry <- as_int(n.ancestry)[1]
  if (!is.finite(n.ancestry) || is.na(n.ancestry) || n.ancestry < 1L) n.ancestry <- 1L
  ns <- as_int(n.sample)
  ns <- ns[is.finite(ns) & !is.na(ns)]
  if (!length(ns)) ns <- as_int(10000L)

  if (length(ns) == 1L) {
    total <- as_int(ns[1])
    if (!is.finite(total) || is.na(total) || total < n.ancestry) total <- n.ancestry * 100L
    min_each <- if (total >= as_int(n.ancestry * 10L)) 10L else 1L
    w <- rep(1, n.ancestry) + stats::runif(n.ancestry, min = -0.15, max = 0.15)
    w <- pmax(0.2, w)
    w <- w / sum(w)
    alloc <- as_int(floor(total * w))
    alloc <- pmax(min_each, alloc)

    diff <- as_int(total - sum(alloc))
    if (diff > 0L) {
      add_idx <- sample.int(n.ancestry, size = diff, replace = TRUE)
      tab <- table(add_idx)
      alloc[as_int(names(tab))] <- alloc[as_int(names(tab))] + as_int(tab)
    } else if (diff < 0L) {
      need <- as_int(abs(diff))
      for (k in seq_len(need)) {
        can <- which(alloc > min_each)
        if (!length(can)) can <- seq_len(n.ancestry)
        j <- sample(can, size = 1L)
        alloc[j] <- alloc[j] - 1L
      }
    }
    if (sum(alloc) != total) {
      tail_diff <- as_int(total - sum(alloc))
      alloc[1] <- alloc[1] + tail_diff
    }
    alloc[alloc < 1L] <- 1L
    return(list(per = as_int(alloc), total = as_int(sum(alloc)), by_vector = FALSE))
  }

  if (length(ns) != n.ancestry) {
    stop("For n.ancestry > 1, n.sample must be scalar(total sample size) or length == n.ancestry.", call. = FALSE)
  }
  ns <- as_int(ns)
  ns[!is.finite(ns) | is.na(ns) | ns < 1L] <- 1L
  list(per = as_int(ns), total = as_int(sum(ns)), by_vector = TRUE)
}

.gcanvas_toy_parse_var_list <- function(var.list, chr_allow = NULL) {
  dt <- data.table::data.table(CHR = character(), POS = numeric(), A1 = character(), A2 = character(), SNP = character())
  if (is.null(var.list) || length(var.list) == 0L) return(dt)

  out <- NULL
  if (data.table::is.data.table(var.list) || is.data.frame(var.list)) {
    x <- data.table::as.data.table(var.list)
    if ("SNP" %in% names(x)) {
      vx <- as.character(x$SNP)
      vx <- vx[!is.na(vx) & nzchar(vx)]
      if (length(vx)) out <- vx
    } else if (all(c("CHR", "POS", "A1", "A2") %in% names(x))) {
      y <- x[, .(
        CHR = normalize.chrom(CHR),
        POS = suppressWarnings(as.numeric(POS)),
        A1 = toupper(as.character(A1)),
        A2 = toupper(as.character(A2))
      )]
      y <- y[!is.na(CHR) & nzchar(CHR) & is.finite(POS) & POS >= 1 & !is.na(A1) & !is.na(A2) & nzchar(A1) & nzchar(A2) & A1 != A2]
      if (!nrow(y)) return(dt)
      y[, POS := as_int(round(POS))]
      y[, SNP := paste0(CHR, ":", POS, ":", A1, ":", A2)]
      if (!is.null(chr_allow)) {
        keep <- normalize.chrom(chr_allow)
        y <- y[CHR %in% keep]
      }
      y[, chr_order := rank.chrom(CHR)]
      data.table::setorder(y, chr_order, POS, SNP)
      y[, chr_order := NULL]
      return(unique(y[, .(CHR, POS, A1, A2, SNP)], by = "SNP"))
    }
  } else {
    out <- as.character(var.list)
  }

  out <- as.character(out %||% character())
  out <- trimws(out)
  out <- out[!is.na(out) & nzchar(out)]
  if (!length(out)) return(dt)

  sp <- strsplit(out, ":", fixed = TRUE)
  keep <- lengths(sp) >= 4L
  if (!any(keep)) return(dt)
  sp <- sp[keep]
  y <- data.table::rbindlist(lapply(sp, function(v) {
    data.table::data.table(
      CHR = normalize.chrom(v[1]),
      POS = suppressWarnings(as.numeric(v[2])),
      A1 = toupper(as.character(v[3])),
      A2 = toupper(as.character(v[4]))
    )
  }), use.names = TRUE, fill = TRUE)
  y <- y[!is.na(CHR) & nzchar(CHR) & is.finite(POS) & POS >= 1 & !is.na(A1) & !is.na(A2) & nzchar(A1) & nzchar(A2) & A1 != A2]
  if (!nrow(y)) return(dt)
  y[, POS := as_int(round(POS))]
  if (!is.null(chr_allow)) {
    keep_chr <- normalize.chrom(chr_allow)
    y <- y[CHR %in% keep_chr]
  }
  if (!nrow(y)) return(dt)
  y[, SNP := paste0(CHR, ":", POS, ":", A1, ":", A2)]
  y[, chr_order := rank.chrom(CHR)]
  data.table::setorder(y, chr_order, POS, SNP)
  y[, chr_order := NULL]
  unique(y[, .(CHR, POS, A1, A2, SNP)], by = "SNP")
}

.gcanvas_toy_normalize_iid_prefix <- function(iid.prefix) {
  if (is.null(iid.prefix) || length(iid.prefix) == 0L) return(NULL)
  p <- trimws(as.character(iid.prefix)[1])
  if (is.na(p) || !nzchar(p)) return(NULL)
  if (tolower(p) %in% c("null", "na", "none", "false")) return(NULL)
  p
}

.gcanvas_toy_apply_iid_prefix <- function(ids, iid.prefix = NULL) {
  out <- as.character(ids)
  pref <- .gcanvas_toy_normalize_iid_prefix(iid.prefix)
  if (is.null(pref)) return(out)
  paste0(pref, "_", out)
}

.toy_gwas_pick_centers <- function(n.causal, chr_bounds, causal.window) {
  n.causal <- as_int(n.causal)
  if (!is.finite(n.causal) || n.causal < 1L) {
    return(data.table::data.table(locus_id = integer(), CHR = character(), center_pos = numeric()))
  }
  centers <- data.table::data.table(
    locus_id = seq_len(n.causal),
    CHR = character(n.causal),
    center_pos = numeric(n.causal)
  )
  chr_prob <- as.numeric(chr_bounds$end)
  chr_prob <- chr_prob / sum(chr_prob)
  for (i in seq_len(n.causal)) {
    picked <- FALSE
    for (k in seq_len(800L)) {
      chr_i <- sample(chr_bounds$CHR, size = 1L, prob = chr_prob)
      chr_len <- as_int(chr_bounds[CHR == chr_i, end][1])
      margin <- as_int(min(causal.window, floor(chr_len / 4)))
      lo <- max(2L, margin + 1L)
      hi <- max(lo, chr_len - margin - 1L)
      pos_i <- if (hi > lo) sample.int(hi - lo + 1L, 1L) + lo - 1L else sample.int(chr_len, 1L)
      if (i == 1L) {
        centers[i, `:=`(CHR = chr_i, center_pos = as.numeric(pos_i))]
        picked <- TRUE
        break
      }
      prev <- centers[seq_len(i - 1L)]
      ov <- prev[CHR == chr_i & abs(center_pos - pos_i) < (2 * causal.window)]
      if (!nrow(ov)) {
        centers[i, `:=`(CHR = chr_i, center_pos = as.numeric(pos_i))]
        picked <- TRUE
        break
      }
    }
    if (!picked) {
      chr_i <- sample(chr_bounds$CHR, size = 1L, prob = chr_prob)
      chr_len <- as_int(chr_bounds[CHR == chr_i, end][1])
      centers[i, `:=`(CHR = chr_i, center_pos = as.numeric(sample.int(chr_len, 1L)))]
      .gcanvas_warn_msg(sprintf("Locus %d overlap guard failed; placed anyway.", i))
    }
  }
  centers
}

.toy_gwas_sample_bg_positions <- function(n, chr_bounds, block_windows = NULL) {
  require_pkg("data.table")
  n <- as_int(n)
  if (!is.finite(n) || n <= 0L) {
    return(data.table::data.table(CHR = character(), POS = numeric()))
  }
  chr_prob <- as.numeric(chr_bounds$end)
  chr_prob <- chr_prob / sum(chr_prob)

  blk <- data.table::as.data.table(block_windows %||% data.table::data.table())
  if (nrow(blk)) {
    blk <- blk[, .(CHR = as.character(CHR), start = suppressWarnings(as.numeric(start)), end = suppressWarnings(as.numeric(end)))]
    blk <- blk[!is.na(CHR) & nzchar(CHR) & is.finite(start) & is.finite(end) & end >= start]
    data.table::setkey(blk, CHR, start, end)
  } else {
    blk <- data.table::data.table(CHR = character(), start = numeric(), end = numeric())
  }

  out <- data.table::data.table(CHR = character(), POS = numeric())
  out_key <- character()
  iter <- 0L
  while (nrow(out) < n && iter < 100L) {
    iter <- iter + 1L
    need <- as_int(n - nrow(out))
    k <- max(2000L, need * 4L)
    chr_s <- sample(chr_bounds$CHR, size = k, replace = TRUE, prob = chr_prob)
    len_s <- as.numeric(chr_bounds$end[match(chr_s, chr_bounds$CHR)])
    pos_s <- as.numeric(floor(stats::runif(k, min = 1, max = pmax(2, len_s + 1))))
    cand <- data.table::data.table(CHR = chr_s, POS = as.numeric(round(pos_s)))
    cand <- cand[!is.na(CHR) & nzchar(CHR) & is.finite(POS) & POS >= 1]
    cand <- unique(cand, by = c("CHR", "POS"))
    if (!nrow(cand)) next

    if (nrow(blk)) {
      cand_iv <- cand[, .(CHR, start = POS, end = POS, POS)]
      data.table::setkey(cand_iv, CHR, start, end)
      hit <- data.table::foverlaps(cand_iv, blk, nomatch = 0L)
      if (nrow(hit)) {
        drop_key <- unique(paste0(hit$CHR, ":", as_int(hit$POS)))
        cand <- cand[!(paste0(CHR, ":", as_int(POS)) %in% drop_key)]
      }
    }
    if (!nrow(cand)) next
    if (length(out_key)) {
      cand <- cand[!(paste0(CHR, ":", as_int(POS)) %in% out_key)]
    }
    if (!nrow(cand)) next

    add <- cand[seq_len(min(need, nrow(cand)))]
    out <- data.table::rbindlist(list(out, add), use.names = TRUE)
    out_key <- c(out_key, paste0(add$CHR, ":", as_int(add$POS)))
  }

  if (nrow(out) < n) {
    stop("Failed to sample enough background positions. Reduce n.snp or causal.window.", call. = FALSE)
  }
  out[seq_len(n), .(CHR, POS)]
}

#' Simulate a toy GWAS summary-statistics table
#'
#' Generates a synthetic GWAS for quick examples and tests. Produces a
#' `data.table` of summary stats (`SNP`, `CHR`, `POS`, `A1`, `A2`, `EAF`,
#' `BETA`, `SE`, `Z`, `P`, `N`, `is_lead`, ...) and -- when `bfile` is set --
#' a matching PLINK1 bfile suitable for downstream LD / PCA workflows.
#'
#' Supports single-ancestry (default, returns a `data.table`) and multi-
#' ancestry mode (`n.ancestry > 1`, returns a named list of `data.table`s
#' plus a merged PLINK fileset).
#'
#' @param n.sample Sample size used to scale standard errors.
#' @param n.ancestry Number of ancestries to simulate. `>1` returns a list.
#' @param n.snp Approximate number of variants.
#' @param n.causal Number of causal loci.
#' @param iid.prefix Prefix for synthetic sample IDs.
#' @param build Genome build (`37` or `38`).
#' @param chr Chromosomes to populate.
#' @param var.list Optional vector of SNP IDs to use as the variant scaffold.
#' @param causal.window,causal.nsnp,causal.min.sig Causal-locus geometry.
#' @param p.threshold,p.floor P-value threshold for "significant" tagging and
#'   the floor used when truncating extreme p-values.
#' @param eaf.range Effect-allele frequency draw range.
#' @param h2 Per-locus heritability.
#' @param causal.mix Mixing weights for effect-direction patterns at causal loci.
#' @param ld.decay.bp Approximate LD decay distance in bp.
#' @param bfile Logical or path prefix. `TRUE` (default) writes a default
#'   PLINK1 bfile next to the working directory; a character path writes to
#'   that prefix; `FALSE` skips bfile generation.
#' @param overwrite Logical. Overwrite existing bfile / output.
#' @param show.progress Logical. Show a progress bar.
#' @param seed Integer seed for reproducibility.
#'
#' @return A `data.table` of summary stats (single ancestry) or a named
#'   list of such tables (multi-ancestry). Carries `attr(., "gcanvas_meta")`
#'   with run parameters.
#'
#' @seealso [toy.eqtl()] for matching eQTL-style data,
#'   [get.lead()] to extract lead variants from the result.
#' @export
toy.gwas <- function(n.sample = 10000L,
                     n.ancestry = 1L,
                     n.snp = 10000L,
                     n.causal = 3L,
                     iid.prefix = NULL,
                     build = 38L,
                     chr = 1:22,
                     var.list = NULL,
                     causal.window = 2e5L,
                     causal.nsnp = 120L,
                     causal.min.sig = 8L,
                     p.threshold = 5e-8,
                     p.floor = 1e-50,
                     eaf.range = c(0.01, 0.5),
                     h2 = 0.2,
                     causal.mix = c(0.55, 0.30, 0.15),
                     ld.decay.bp = 3e4,
                     bfile = TRUE,
                     overwrite = TRUE,
                     show.progress = TRUE,
                     seed = NULL) {
  require_pkg(c("data.table", "bigsnpr", "bigstatsr"))
  show_progress <- isTRUE(show.progress)
  .toy_gwas_progress_new <- function(total, label = NULL) {
    if (!show_progress) return(NULL)
    .gcanvas_note("gcanvas::toy.gwas", label, silent = !show_progress)
    utils::txtProgressBar(min = 0, max = max(1, as.numeric(total)), initial = 0, style = 3)
  }
  .toy_gwas_progress_set <- function(pb, value) {
    if (is.null(pb)) return(invisible(NULL))
    utils::setTxtProgressBar(pb, as.numeric(value))
    invisible(NULL)
  }
  .toy_gwas_progress_close <- function(pb) {
    if (is.null(pb)) return(invisible(NULL))
    close(pb)
    invisible(NULL)
  }

  seed_use <- .gcanvas_seed_resolve(seed)
  if (!is.null(seed_use)) set.seed(seed_use)
  iid_prefix <- .gcanvas_toy_normalize_iid_prefix(iid.prefix)
  n.ancestry <- as_int(n.ancestry)[1]
  if (!is.finite(n.ancestry) || is.na(n.ancestry) || n.ancestry < 1L) n.ancestry <- 1L
  seed_label <- .gcanvas_seed_label(seed = seed, seed_use = seed_use)
  .gcanvas_note(
    "gcanvas::toy.gwas",
    sprintf(
      "Start: n.sample=%s | n.ancestry=%d | n.snp=%d | n.causal=%d | build=%s | seed=%s",
      paste(as_int(n.sample), collapse = ","), as_int(n.ancestry), as_int(n.snp), as_int(n.causal), as.character(build)[1], seed_label
    ),
    silent = FALSE
  )

  n.causal <- as_int(n.causal)
  n.snp <- as_int(n.snp)
  causal.window <- as_int(causal.window)
  causal.nsnp <- as_int(causal.nsnp)
  causal.min.sig <- as_int(causal.min.sig)
  n.sample <- as_int(n.sample)
  ld.decay.bp <- suppressWarnings(as.numeric(ld.decay.bp))[1]
  h2 <- suppressWarnings(as.numeric(h2))[1]
  overwrite <- isTRUE(overwrite)

  if (!is.finite(n.causal) || n.causal < 1L) stop("n.causal must be >= 1.", call. = FALSE)
  if (!is.finite(n.snp) || n.snp < 200L) n.snp <- 200L
  if (!is.finite(causal.window) || causal.window < 1000L) causal.window <- 1000L
  if (!is.finite(causal.nsnp) || causal.nsnp < 10L) causal.nsnp <- 10L
  if (!is.finite(causal.min.sig) || causal.min.sig < 1L) causal.min.sig <- 1L
  if (!length(n.sample) || all(!is.finite(n.sample) | is.na(n.sample))) n.sample <- as_int(10000L)
  if (isTRUE(as_int(n.ancestry) > 1L)) {
    n.sample[!is.finite(n.sample) | is.na(n.sample) | n.sample < 1L] <- 1L
  } else {
    n.sample <- as_int(n.sample[1])
    if (!is.finite(n.sample) || is.na(n.sample) || n.sample < 50L) n.sample <- 50L
  }
  if (!is.finite(ld.decay.bp) || ld.decay.bp <= 100) ld.decay.bp <- 3e4
  if (!is.finite(h2) || h2 <= 0 || h2 >= 1) h2 <- 0.2
  h2_target <- as.numeric(h2)

  p.threshold <- suppressWarnings(as.numeric(p.threshold))[1]
  if (!is.finite(p.threshold) || p.threshold <= 0 || p.threshold >= 1) p.threshold <- 5e-8
  p.floor <- suppressWarnings(as.numeric(p.floor))[1]
  if (!is.finite(p.floor) || p.floor <= 0 || p.floor >= p.threshold) p.floor <- 1e-50

  eaf.range <- suppressWarnings(as.numeric(eaf.range))
  if (length(eaf.range) < 2L || any(!is.finite(eaf.range[1:2]))) eaf.range <- c(0.01, 0.5)
  eaf.range <- sort(eaf.range[1:2])
  eaf.range[1] <- max(1e-4, eaf.range[1])
  eaf.range[2] <- min(0.999, eaf.range[2])
  if (eaf.range[2] <= eaf.range[1]) eaf.range <- c(0.01, 0.5)

  mix <- suppressWarnings(as.numeric(causal.mix))
  mix <- mix[is.finite(mix) & mix >= 0]
  if (length(mix) < 3L) mix <- c(0.55, 0.30, 0.15)
  mix <- mix[1:3]
  if (sum(mix) <= 0) mix <- c(0.55, 0.30, 0.15)
  mix <- mix / sum(mix)

  if (isTRUE(as_int(n.ancestry) > 1L)) {
    ns_alloc <- .gcanvas_toy_split_samples(n.sample = n.sample, n.ancestry = as_int(n.ancestry))
    n.sample.each <- as_int(ns_alloc$per)
    n.sample.total <- as_int(ns_alloc$total)

    seed_anchor <- seed_use
    if (is.null(seed_anchor) || !is.finite(seed_anchor) || is.na(seed_anchor)) {
      seed_anchor <- as_int(sample.int(.Machine$integer.max - 1L, size = 1L))
    }
    seed_anchor <- as_int(abs(seed_anchor))
    if (!is.finite(seed_anchor) || is.na(seed_anchor) || seed_anchor == 0L) seed_anchor <- 1L

    bopt_multi <- .gcanvas_toy_parse_bfile_option(bfile = bfile, prefix_pattern = "toy_gwas_multi_")
    bfile_base <- NULL
    if (isTRUE(bopt_multi$write_bfile)) {
      bfile_base <- .gcanvas_pretty_path(abs_path(bopt_multi$bfile_prefix))
      dir.create(dirname(bfile_base), recursive = TRUE, showWarnings = FALSE)
    }

    .toy_gwas_merge_bfiles <- function(prefixes, out_prefix, overwrite = TRUE) {
      prefixes <- as.character(prefixes)
      prefixes <- prefixes[!is.na(prefixes) & nzchar(prefixes)]
      prefixes <- unique(prefixes)
      if (!length(prefixes)) return(NULL)
      out_prefix <- .gcanvas_pretty_path(abs_path(out_prefix))
      targets <- paste0(out_prefix, c(".bed", ".bim", ".fam"))
      if (any(file.exists(targets))) {
        if (isTRUE(overwrite)) unlink(targets, force = TRUE) else {
          stop("Merged bfile already exists. Set overwrite=TRUE or use another prefix.", call. = FALSE)
        }
      }
      if (length(prefixes) == 1L) {
        src <- paste0(prefixes[1], c(".bed", ".bim", ".fam"))
        if (!all(file.exists(src))) stop("Cannot copy merged bfile: missing source files.", call. = FALSE)
        ok <- file.copy(src, targets, overwrite = isTRUE(overwrite))
        if (!all(ok)) stop("Failed to write merged bfile.", call. = FALSE)
        return(out_prefix)
      }

      tmp_backing <- character()
      fbm_base <- tempfile(pattern = "toy_merge_fbm_")
      on.exit({
        rmv <- c(
          tmp_backing,
          paste0(tmp_backing, ".bk"),
          paste0(tmp_backing, ".rds"),
          fbm_base,
          paste0(fbm_base, ".bk"),
          paste0(fbm_base, ".rds")
        )
        rmv <- rmv[file.exists(rmv)]
        if (length(rmv)) unlink(rmv, force = TRUE)
      }, add = TRUE)

      obj_list <- vector("list", length(prefixes))
      fam_list <- vector("list", length(prefixes))
      map0 <- NULL
      n_total <- 0L
      for (i in seq_along(prefixes)) {
        pref <- prefixes[i]
        bed <- paste0(pref, ".bed")
        bim <- paste0(pref, ".bim")
        fam <- paste0(pref, ".fam")
        if (!all(file.exists(c(bed, bim, fam)))) {
          stop("Missing ancestry bfile component(s): ", pref, call. = FALSE)
        }
        bb <- tempfile(pattern = sprintf("toy_merge_%d_", as_int(i)))
        tmp_backing <- c(tmp_backing, bb)
        bigsnpr::snp_readBed(bed, backingfile = bb)
        obj <- bigsnpr::snp_attach(paste0(bb, ".rds"))
        obj_list[[i]] <- obj
        fam_list[[i]] <- data.table::as.data.table(obj$fam)
        n_total <- n_total + as_int(nrow(obj$fam))

        mapi <- data.table::as.data.table(obj$map)
        if (is.null(map0)) {
          map0 <- mapi
        } else {
          if (!identical(as.character(map0$marker.ID), as.character(mapi$marker.ID))) {
            stop("Ancestry bfiles do not share identical marker.ID order; cannot merge.", call. = FALSE)
          }
          if (!identical(as.character(map0$allele1), as.character(mapi$allele1)) ||
              !identical(as.character(map0$allele2), as.character(mapi$allele2))) {
            stop("Ancestry bfiles have allele mismatch; cannot merge.", call. = FALSE)
          }
        }
      }
      if (is.null(map0) || !nrow(map0) || n_total < 1L) return(NULL)

      m <- as_int(nrow(map0))
      G_fbm <- bigstatsr::FBM.code256(
        nrow = as_int(n_total),
        ncol = as_int(m),
        code = c(0, 1, 2, 3, rep(NA_real_, 252)),
        init = 0,
        backingfile = fbm_base
      )

      row_cursor <- 1L
      chunk <- 600L
      for (i in seq_along(obj_list)) {
        obj <- obj_list[[i]]
        ni <- as_int(nrow(obj$fam))
        if (ni < 1L) next
        rr <- row_cursor:(row_cursor + ni - 1L)
        for (from in seq.int(1L, m, by = chunk)) {
          to <- min(m, from + chunk - 1L)
          G_fbm[rr, from:to] <- obj$genotypes[, from:to, drop = FALSE]
        }
        row_cursor <- row_cursor + ni
      }

      fam_all <- data.table::rbindlist(fam_list, use.names = TRUE, fill = TRUE)
      req_fam <- c("family.ID", "sample.ID", "paternal.ID", "maternal.ID", "sex", "affection")
      for (nm in req_fam) if (!(nm %in% names(fam_all))) fam_all[, (nm) := 0]
      fam_all <- as.data.frame(fam_all[, ..req_fam], stringsAsFactors = FALSE)
      map_out <- as.data.frame(map0, stringsAsFactors = FALSE)
      obj_out <- structure(list(genotypes = G_fbm, fam = fam_all, map = map_out), class = "bigSNP")
      bigsnpr::snp_writeBed(obj_out, paste0(out_prefix, ".bed"))
      out_prefix
    }

    .toy_gwas_write_ancestry_bfile <- function(scaffold, present, eaf_vec, n_i, local_G, out_prefix, ancestry_idx, iid.prefix = NULL, apply_fid_prefix = FALSE, overwrite = TRUE) {
      out_prefix <- .gcanvas_pretty_path(abs_path(out_prefix))
      dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)
      targets <- paste0(out_prefix, c(".bed", ".bim", ".fam"))
      if (any(file.exists(targets))) {
        if (isTRUE(overwrite)) {
          unlink(targets, force = TRUE)
        } else {
          stop("bfile already exists. Set overwrite=TRUE or provide another bfile prefix.", call. = FALSE)
        }
      }

      n_total <- as_int(nrow(scaffold))
      if ((as.numeric(n_total) * as.numeric(n_i)) > 3e8) {
        stop("Requested bfile is too large. Reduce n.snp/n.sample or set bfile=FALSE.", call. = FALSE)
      }

      fbm_base <- tempfile(pattern = "toy_gwas_multi_fbm_")
      on.exit({
        tmp_files <- c(fbm_base, paste0(fbm_base, ".bk"), paste0(fbm_base, ".rds"))
        unlink(tmp_files[file.exists(tmp_files)], force = TRUE)
      }, add = TRUE)

      G_fbm <- bigstatsr::FBM.code256(
        nrow = as_int(n_i), ncol = as_int(n_total),
        code = c(0, 1, 2, 3, rep(NA_real_, 252)),
        init = 0, backingfile = fbm_base
      )

      chunk <- 600L
      for (from in seq.int(1L, n_total, by = chunk)) {
        to <- min(n_total, from + chunk - 1L)
        idx <- from:to
        block <- matrix(0L, nrow = as_int(n_i), ncol = length(idx))
        present_idx <- idx[present[idx]]
        if (length(present_idx)) {
          loc_mask <- !is.na(scaffold$is_local[present_idx]) &
            as.logical(scaffold$is_local[present_idx]) &
            is.finite(scaffold$local_col_idx[present_idx])
          loc_idx <- present_idx[loc_mask]
          if (length(loc_idx) && !is.null(local_G) && ncol(local_G) > 0L) {
            lcol <- as_int(scaffold$local_col_idx[loc_idx])
            block[, match(loc_idx, idx)] <- local_G[, lcol, drop = FALSE]
          }
          bg_idx <- setdiff(present_idx, loc_idx)
          if (length(bg_idx)) {
            p_bg <- pmin(pmax(as.numeric(eaf_vec[bg_idx]), 1e-4), 0.999)
            block[, match(bg_idx, idx)] <- matrix(
              stats::rbinom(as_int(n_i * length(bg_idx)), size = 2L, prob = rep(p_bg, each = as_int(n_i))),
              nrow = as_int(n_i), ncol = length(bg_idx)
            )
          }
        }
        G_fbm[, idx] <- block
      }

      fam <- data.frame(
        family.ID = rep(sprintf("anc%d", as_int(ancestry_idx)), as_int(n_i)),
        sample.ID = sprintf("anc%d_%06d", as_int(ancestry_idx), seq_len(as_int(n_i))),
        paternal.ID = 0,
        maternal.ID = 0,
        sex = 0,
        affection = -9,
        stringsAsFactors = FALSE
      )
      if (isTRUE(apply_fid_prefix)) {
        fam$family.ID <- .gcanvas_toy_apply_iid_prefix(fam$family.ID, iid.prefix = iid.prefix)
      }
      fam$sample.ID <- .gcanvas_toy_apply_iid_prefix(fam$sample.ID, iid.prefix = iid.prefix)
      map <- data.frame(
        chromosome = as_int(rank.chrom(scaffold$CHR)),
        marker.ID = as.character(scaffold$SNP),
        genetic.dist = 0,
        physical.pos = as_int(scaffold$POS),
        allele1 = as.character(scaffold$A1),
        allele2 = as.character(scaffold$A2),
        stringsAsFactors = FALSE
      )
      obj <- structure(list(genotypes = G_fbm, fam = fam, map = map), class = "bigSNP")
      bigsnpr::snp_writeBed(obj, paste0(out_prefix, ".bed"))
      out_prefix
    }

    .gcanvas_note(
      "gcanvas::toy.gwas",
      sprintf(
        "Multi-ancestry mode: n.ancestry=%d | n.sample.total=%d | n.sample.each=[%s] | seed.anchor=%d",
        as_int(n.ancestry),
        as_int(n.sample.total),
        paste(as_int(n.sample.each), collapse = ","),
        as_int(seed_anchor)
      ),
      silent = !show_progress
    )

    .gcanvas_note("gcanvas::toy.gwas", "Generating shared SNP scaffold", silent = !show_progress)
    forced_vars <- .gcanvas_toy_parse_var_list(var.list = var.list, chr_allow = chr)
    forced_snp <- as.character(forced_vars$SNP)
    scaffold_dt <- toy.gwas(
      n.sample = as_int(max(200L, max(n.sample.each))),
      n.ancestry = 1L,
      n.snp = as_int(n.snp),
      n.causal = as_int(n.causal),
      build = build,
      chr = chr,
      var.list = forced_snp,
      causal.window = as_int(causal.window),
      causal.nsnp = as_int(causal.nsnp),
      causal.min.sig = as_int(causal.min.sig),
      p.threshold = as.numeric(p.threshold),
      p.floor = as.numeric(p.floor),
      eaf.range = as.numeric(eaf.range),
      h2 = as.numeric(h2),
      causal.mix = as.numeric(mix),
      ld.decay.bp = as.numeric(ld.decay.bp),
      bfile = FALSE,
      overwrite = TRUE,
      show.progress = FALSE,
      seed = as_int(seed_anchor)
	    )
	    scaffold <- data.table::copy(scaffold_dt)
	    scaffold[, chr_order := rank.chrom(CHR)]
	    data.table::setorderv(scaffold, c("chr_order", "CHR", "POS", "SNP"), c(1L, 1L, 1L, 1L), na.last = TRUE)
	    scaffold[, chr_order := NULL]
	    scaffold <- scaffold[, .(
	      SNP, CHR, POS, A1, A2, EAF,
	      locus_id, is_local, is_causal, center_pos, true_beta, local_col_idx
	    )]
    scaffold[, row_id := .I]
    n_total_var <- as_int(nrow(scaffold))
    if (n_total_var < 1L) stop("Failed to generate shared scaffold.", call. = FALSE)

    local_rows <- which(!is.na(scaffold$is_local) & as.logical(scaffold$is_local) & is.finite(scaffold$local_col_idx))
    m_local <- if (length(local_rows)) as_int(max(scaffold$local_col_idx[local_rows])) else 0L
    locus_ids <- sort(unique(as_int(scaffold$locus_id[is.finite(scaffold$locus_id) & scaffold$locus_id > 0L])))
    locus_rows <- setNames(vector("list", length(locus_ids)), as.character(locus_ids))
    for (k in seq_along(locus_ids)) {
      locus_rows[[k]] <- which(as_int(scaffold$locus_id) == as_int(locus_ids[k]))
    }
    locus_active <- matrix(TRUE, nrow = length(locus_ids), ncol = as_int(n.ancestry))
    colnames(locus_active) <- paste0("ancestry", seq_len(as_int(n.ancestry)))
    rownames(locus_active) <- as.character(locus_ids)
    if (length(locus_ids) && as_int(n.ancestry) > 1L) {
      k_opts <- seq_len(as_int(n.ancestry))
      lambda <- max(1.2, min(2.2, as.numeric(n.ancestry) / 1.6))
      k_prob <- stats::dpois(k_opts, lambda = lambda) + 0.02
      k_prob[k_opts == as_int(n.ancestry)] <- k_prob[k_opts == as_int(n.ancestry)] * 0.65
      k_prob <- k_prob / sum(k_prob)
      for (k in seq_along(locus_ids)) {
        k_use <- sample(k_opts, size = 1L, prob = k_prob)
        pick <- sample.int(as_int(n.ancestry), size = as_int(k_use), replace = FALSE)
        locus_active[k, ] <- FALSE
        locus_active[k, pick] <- TRUE
      }
      for (i in seq_len(as_int(n.ancestry))) {
        if (!any(locus_active[, i])) {
          locus_active[sample(seq_along(locus_ids), size = 1L), i] <- TRUE
        }
      }
    }

    ld_mult <- if (as_int(n.ancestry) > 1L) seq(0.80, 1.20, length.out = as_int(n.ancestry)) else 1
    h2_mult <- if (as_int(n.ancestry) > 1L) seq(0.90, 1.10, length.out = as_int(n.ancestry)) else 1
    maf_shift <- if (as_int(n.ancestry) > 1L) seq(-0.05, 0.05, length.out = as_int(n.ancestry)) else 0

    out <- vector("list", as_int(n.ancestry))
    ancestry_bfiles <- rep(NA_character_, as_int(n.ancestry))
    for (i in seq_len(as_int(n.ancestry))) {
      anc_tag <- paste0("ancestry", as_int(i))
      ld_i <- max(500, as.numeric(ld.decay.bp) * as.numeric(ld_mult[i]))
      h2_i <- max(0.02, min(0.95, as.numeric(h2) * as.numeric(h2_mult[i])))
      lo_i <- max(1e-4, min(0.49, as.numeric(eaf.range[1]) + as.numeric(maf_shift[i])))
      hi_i <- min(0.999, max(lo_i + 1e-3, as.numeric(eaf.range[2]) + as.numeric(maf_shift[i])))
      eaf_i <- c(lo_i, hi_i)
      n_i <- as_int(n.sample.each[i])
      active_loci_i <- if (length(locus_ids)) as_int(locus_ids[as.logical(locus_active[, i])]) else integer(0)
      active_rows <- as_int(scaffold$locus_id) %in% as_int(active_loci_i)
      bfile_i <- if (isTRUE(bopt_multi$write_bfile)) sprintf("%s.%s", bfile_base, anc_tag) else FALSE

      .gcanvas_note(
        "gcanvas::toy.gwas",
        sprintf(
          "Generating %s: n.sample=%d | ld.decay.bp=%.0f | h2=%.3f | eaf.range=[%.3f, %.3f]",
          anc_tag, as_int(n_i), as.numeric(ld_i), as.numeric(h2_i), as.numeric(eaf_i[1]), as.numeric(eaf_i[2])
        ),
        silent = !show_progress
      )

      base_eaf <- pmin(pmax(as.numeric(scaffold$EAF), 1e-4), 0.999)
      eaf_vec <- base_eaf + as.numeric(maf_shift[i]) + stats::rnorm(n_total_var, mean = 0, sd = 0.02)
      eaf_vec <- pmin(pmax(eaf_vec, as.numeric(eaf_i[1])), as.numeric(eaf_i[2]))

      base_drop <- min(0.25, max(0.01, 0.02 + abs(as.numeric(maf_shift[i])) * 0.60 + stats::runif(1, min = 0, max = 0.02)))
      low_maf_penalty <- pmax(0, 0.03 - eaf_vec) * 2.0
      drop_prob <- pmin(0.85, pmax(0, base_drop + low_maf_penalty))
      present <- stats::runif(n_total_var) > drop_prob
      if (length(forced_snp)) {
        present[as.character(scaffold$SNP) %in% forced_snp] <- TRUE
      }

      for (k in seq_along(locus_rows)) {
        idx <- locus_rows[[k]]
        if (!length(idx)) next
        if (!any(present[idx])) present[sample(idx, size = 1L)] <- TRUE
        if (length(active_loci_i) && as_int(locus_ids[k]) %in% as_int(active_loci_i)) {
          cidx <- idx[
            !is.na(scaffold$is_causal[idx]) &
              as.logical(scaffold$is_causal[idx]) &
              active_rows[idx]
          ]
          if (length(cidx) && !any(present[cidx])) present[sample(cidx, size = 1L)] <- TRUE
        }
      }

      local_G <- if (m_local > 0L) matrix(0L, nrow = as_int(n_i), ncol = as_int(m_local)) else NULL
      if (m_local > 0L) {
        for (k in seq_along(locus_rows)) {
          idx <- locus_rows[[k]]
          if (!length(idx)) next
          idx <- idx[
            !is.na(scaffold$is_local[idx]) &
              as.logical(scaffold$is_local[idx]) &
              is.finite(scaffold$local_col_idx[idx]) &
              present[idx]
          ]
          if (!length(idx)) next
          ord <- order(as_int(scaffold$local_col_idx[idx]))
          idx <- idx[ord]
          lcol <- as_int(scaffold$local_col_idx[idx])
          pos <- as.numeric(scaffold$POS[idx])
          maf <- pmin(pmax(as.numeric(eaf_vec[idx]), 1e-4), 0.999)
          cp <- c(0, exp(-diff(pos) / max(1000, as.numeric(ld_i))))
          hap1 <- .gcanvas_toy_sim_haplotype(n = as_int(n_i), maf = maf, copy_prob = cp)
          hap2 <- .gcanvas_toy_sim_haplotype(n = as_int(n_i), maf = maf, copy_prob = cp)
          local_G[, lcol] <- hap1 + hap2
        }
      }

      if (m_local > 0L && as_int(n_i) > 1L) {
        local_present_rows <- which(
          present &
            !is.na(scaffold$is_local) &
            as.logical(scaffold$is_local) &
            is.finite(scaffold$local_col_idx)
        )
        if (length(local_present_rows)) {
          local_present_cols <- as_int(scaffold$local_col_idx[local_present_rows])
          ord <- order(local_present_cols)
          local_present_rows <- local_present_rows[ord]
          local_present_cols <- local_present_cols[ord]
          vv <- apply(local_G[, local_present_cols, drop = FALSE], 2, stats::var)
          bad <- local_present_rows[!is.finite(vv) | vv <= 0]
          if (length(bad)) present[bad] <- FALSE
        }
      }

      var_g <- pmax(2 * as.numeric(eaf_vec) * (1 - as.numeric(eaf_vec)), 1e-8)
      se <- rep(NA_real_, n_total_var)
      se[present] <- 1 / sqrt(as.numeric(n_i) * var_g[present])
      beta_eff <- as.numeric(scaffold$true_beta) * as.numeric(active_rows) * sqrt(max(1e-6, as.numeric(h2_i) / 0.2))
      z <- rep(NA_real_, n_total_var)

      if (m_local > 0L && as_int(n_i) > 1L) {
        local_rows_present <- which(
          present &
            !is.na(scaffold$is_local) &
            as.logical(scaffold$is_local) &
            is.finite(scaffold$local_col_idx)
        )
        if (length(local_rows_present)) {
          cols <- as_int(scaffold$local_col_idx[local_rows_present])
          ord <- order(cols)
          local_rows_present <- local_rows_present[ord]
          cols <- cols[ord]
          G_use <- local_G[, cols, drop = FALSE]
          vv <- apply(G_use, 2, stats::var)
          keep <- is.finite(vv) & vv > 0
          if (!all(keep)) {
            drop_rows <- local_rows_present[!keep]
            if (length(drop_rows)) {
              present[drop_rows] <- FALSE
              se[drop_rows] <- NA_real_
            }
            local_rows_present <- local_rows_present[keep]
            cols <- cols[keep]
            G_use <- G_use[, keep, drop = FALSE]
          }
          if (length(local_rows_present)) {
            R_use <- suppressWarnings(stats::cor(G_use))
            R_use[!is.finite(R_use)] <- 0
            R_use <- (R_use + t(R_use)) / 2
            diag(R_use) <- 1
            eig <- eigen(R_use, symmetric = TRUE)
            eig$values[eig$values < 1e-8] <- 1e-8
            R_use <- eig$vectors %*% (eig$values * t(eig$vectors))
            R_use <- (R_use + t(R_use)) / 2
            diag(R_use) <- 1
            mu <- sqrt(as.numeric(n_i)) * as.numeric(R_use %*% beta_eff[local_rows_present])
            zn <- as.numeric(t(chol(R_use)) %*% stats::rnorm(length(local_rows_present)))
            z[local_rows_present] <- mu + zn
          }
        }
      }

      rem_idx <- which(present & is.na(z) & is.finite(se) & !is.na(se))
      if (length(rem_idx)) {
        z[rem_idx] <- stats::rnorm(length(rem_idx), mean = beta_eff[rem_idx] / se[rem_idx], sd = 1)
      }

      p <- rep(NA_real_, n_total_var)
      ok_z <- which(present & is.finite(z))
      if (length(ok_z)) {
        p[ok_z] <- pmin(pmax(2 * stats::pnorm(-abs(z[ok_z])), as.numeric(p.floor)), 1)
      }
      beta_hat <- rep(NA_real_, n_total_var)
      ok_b <- which(present & is.finite(z) & is.finite(se))
      if (length(ok_b)) beta_hat[ok_b] <- z[ok_b] * se[ok_b]

      for (k in seq_along(locus_rows)) {
        if (!(length(active_loci_i) && as_int(locus_ids[k]) %in% as_int(active_loci_i))) next
        idx_all <- locus_rows[[k]]
        idx <- idx_all[present[idx_all] & is.finite(p[idx_all]) & !is.na(p[idx_all])]
        if (!length(idx)) next
        sig_target_i <- min(as_int(causal.min.sig), as_int(length(idx)))
        if (sig_target_i < 1L) next
        sig_now <- as_int(sum(p[idx] < p.threshold, na.rm = TRUE))
        need <- as_int(sig_target_i - sig_now)
        if (need <= 0L) next
        ord <- idx[order(p[idx])]
        cand <- ord[p[ord] >= p.threshold]
        pick <- cand[seq_len(min(need, length(cand)))]
        if (!length(pick)) next
        p_forced <- p.threshold * stats::runif(length(pick), min = 1e-3, max = 0.8)
        p_forced <- pmax(p_forced, as.numeric(p.floor))
        z_abs <- stats::qnorm(p_forced / 2, lower.tail = FALSE)
        sgn <- sign(z[pick])
        sgn[!is.finite(sgn) | sgn == 0] <- sample(c(-1, 1), size = sum(!is.finite(sgn) | sgn == 0), replace = TRUE)
        p[pick] <- p_forced
        z[pick] <- sgn * z_abs
        beta_hat[pick] <- z[pick] * se[pick]
      }

      neglog <- rep(NA_real_, n_total_var)
      ok_p <- which(is.finite(p) & !is.na(p))
      if (length(ok_p)) neglog[ok_p] <- -log10(pmax(p[ok_p], as.numeric(p.floor)))

      is_lead <- rep(FALSE, n_total_var)
      for (k in seq_along(locus_rows)) {
        idx_all <- locus_rows[[k]]
        idx <- idx_all[present[idx_all] & is.finite(p[idx_all]) & !is.na(p[idx_all])]
        if (!length(idx)) next
        is_lead[idx[which.min(p[idx])]] <- TRUE
      }

      dt_i <- data.table::copy(scaffold)
      dt_i[, `:=`(
        is_causal = as.logical(is_causal & active_rows),
        true_beta = as.numeric(true_beta) * as.numeric(active_rows),
        EAF = ifelse(present, as.numeric(eaf_vec), NA_real_),
        P = as.numeric(p),
        negLog10P = as.numeric(neglog),
        BETA = as.numeric(beta_hat),
        SE = as.numeric(se),
        Z = as.numeric(z),
        N = as_int(ifelse(present, as_int(n_i), 0L)),
        is_lead = as.logical(is_lead)
      )]
      data.table::setcolorder(dt_i, c(
        "SNP", "CHR", "POS", "A1", "A2", "EAF",
        "P", "negLog10P", "BETA", "SE", "Z", "N",
        "locus_id", "is_local", "is_causal", "is_lead", "center_pos", "true_beta", "local_col_idx", "row_id"
      ))
      dt_i[, row_id := NULL]

      if (isTRUE(bopt_multi$write_bfile) && !isFALSE(bfile_i)) {
        ancestry_bfiles[i] <- .toy_gwas_write_ancestry_bfile(
          scaffold = scaffold,
          present = as.logical(present),
          eaf_vec = as.numeric(eaf_vec),
          n_i = as_int(n_i),
          local_G = local_G,
          out_prefix = bfile_i,
          ancestry_idx = as_int(i),
          iid.prefix = iid_prefix,
          apply_fid_prefix = FALSE,
          overwrite = isTRUE(overwrite)
        )
        ancestry_bfiles[i] <- .gcanvas_pretty_path(abs_path(ancestry_bfiles[i]))
      }

      attr(dt_i, "gcanvas_meta") <- list(
        type = "toy.gwas",
        ancestry = anc_tag,
        ancestry_index = as_int(i),
        n_ancestry = as_int(n.ancestry),
        n_samples = as_int(n_i),
        n_active_loci = as_int(length(active_loci_i)),
        active_loci = as_int(active_loci_i),
        n_variants = as_int(nrow(dt_i)),
        n_present = as_int(sum(dt_i$N > 0, na.rm = TRUE)),
        bfile = if (isTRUE(bopt_multi$write_bfile)) ancestry_bfiles[i] else NULL
      )
      out[[i]] <- dt_i
    }
    names(out) <- paste0("ancestry", seq_len(as_int(n.ancestry)))

    merged_bfile <- NULL
    if (isTRUE(bopt_multi$write_bfile)) {
      merged_bfile <- .toy_gwas_merge_bfiles(
        prefixes = ancestry_bfiles,
        out_prefix = sprintf("%s.merged", bfile_base),
        overwrite = isTRUE(overwrite)
      )
      if (!is.null(merged_bfile) && nzchar(as.character(merged_bfile)[1])) {
        merged_bfile <- .gcanvas_pretty_path(abs_path(merged_bfile))
        .gcanvas_note("gcanvas::toy.gwas", paste0("merged bfile written: ", merged_bfile), silent = !show_progress)
      }
    }

    attr(out, "gcanvas_meta") <- list(
      type = "toy.gwas.multi",
      n_ancestry = as_int(n.ancestry),
      n_sample_each = as_int(n.sample.each),
      n_samples_total = as_int(n.sample.total),
      forced_variants = as.character(forced_snp),
      seed_anchor = as_int(seed_anchor),
      ancestry_bfile = ancestry_bfiles,
      merged_bfile = merged_bfile,
      ancestry_parameters = data.table::data.table(
        ancestry = names(out),
        n_samples = as_int(n.sample.each),
        n_active_loci = as_int(colSums(locus_active)),
        ld_decay_bp = as.numeric(max(500, as.numeric(ld.decay.bp) * as.numeric(ld_mult))),
        h2 = as.numeric(pmax(0.02, pmin(0.95, as.numeric(h2) * as.numeric(h2_mult)))),
        eaf_lo = as.numeric(pmax(1e-4, pmin(0.49, as.numeric(eaf.range[1]) + as.numeric(maf_shift)))),
        eaf_hi = as.numeric(pmin(0.999, pmax(
          pmax(1e-4, pmin(0.49, as.numeric(eaf.range[1]) + as.numeric(maf_shift))) + 1e-3,
          as.numeric(eaf.range[2]) + as.numeric(maf_shift)
        )))
      )
    )
    return(out)
  }

  chr_bounds <- .gcanvas_hg_chr_bounds(build = build)
  chr_in <- normalize.chrom(chr)
  chr_in <- unique(chr_in[!is.na(chr_in) & nzchar(chr_in)])
  if (!length(chr_in)) chr_in <- as.character(1:22)
  forced_vars <- .gcanvas_toy_parse_var_list(var.list = var.list, chr_allow = chr_in)
  chr_bounds <- chr_bounds[CHR %in% chr_in]
  if (!nrow(chr_bounds)) stop("No valid chromosomes after filtering.", call. = FALSE)
  chr_bounds[, chr_order := rank.chrom(CHR)]
  data.table::setorder(chr_bounds, chr_order, CHR)

  min_local_per_locus <- 20L
  if (n.snp < as_int(n.causal * min_local_per_locus)) {
    stop(sprintf("n.snp=%d is too small for n.causal=%d. Increase n.snp or reduce n.causal.", n.snp, n.causal), call. = FALSE)
  }
  causal.nsnp_eff <- min(as_int(causal.nsnp), as_int(floor(n.snp / n.causal)))
  causal.nsnp_eff <- max(min_local_per_locus, as_int(causal.nsnp_eff))
  if (causal.nsnp_eff < causal.nsnp) {
    .gcanvas_warn_msg(sprintf("causal.nsnp reduced from %d to %d to satisfy n.snp.", causal.nsnp, causal.nsnp_eff))
  }
  local_total <- as_int(n.causal * causal.nsnp_eff)
  n_bg <- max(0L, as_int(n.snp - local_total))
  sig_target <- as_int(min(causal.min.sig, causal.nsnp_eff))
  if (!is.finite(sig_target) || sig_target < 1L) sig_target <- 1L
  bopt <- .gcanvas_toy_parse_bfile_option(bfile = bfile, prefix_pattern = "toy_gwas_")
  has_ld_matrix_expected <- isTRUE(local_total > 1L)
  total_steps_expected <- 2L + as_int(has_ld_matrix_expected) + as_int(isTRUE(bopt$write_bfile))

  centers <- .toy_gwas_pick_centers(n.causal = n.causal, chr_bounds = chr_bounds, causal.window = causal.window)
  causal_nvar <- sample(c(1L, 2L, 3L), size = n.causal, replace = TRUE, prob = mix)
  centers[, n_causal_var := causal_nvar]

  local_G_list <- vector("list", n.causal)
  local_dt_list <- vector("list", n.causal)
  causal_cols <- integer(0)
  causal_locus_of_causal <- integer(0)
  beta_base <- numeric(0)
  col_cursor <- 0L

  pb_locus <- .toy_gwas_progress_new(
    n.causal,
    sprintf("1/%d generating local loci (n=%d)", as_int(total_steps_expected), as_int(n.causal))
  )
  for (i in seq_len(n.causal)) {
    chr_i <- centers$CHR[i]
    center_i <- as_int(centers$center_pos[i])
    chr_len <- as_int(chr_bounds[CHR == chr_i, end][1])
    left <- max(1L, center_i - causal.window)
    right <- min(chr_len, center_i + causal.window)
    width <- as_int(right - left + 1L)
    if (width < causal.nsnp_eff) {
      stop(sprintf("Locus %d has window width smaller than causal.nsnp.", i), call. = FALSE)
    }

    pos_i <- sort(sample.int(width, size = causal.nsnp_eff, replace = FALSE) + left - 1L)
    maf_i <- pmin(pmax(stats::rbeta(causal.nsnp_eff, shape1 = 0.9, shape2 = 5), eaf.range[1]), eaf.range[2])
    ld_decay_i <- max(1000, ld.decay.bp * stats::runif(1, min = 0.7, max = 1.3))
    cp <- c(0, exp(-diff(pos_i) / ld_decay_i))

    hap1 <- .gcanvas_toy_sim_haplotype(n = n.sample, maf = maf_i, copy_prob = cp)
    hap2 <- .gcanvas_toy_sim_haplotype(n = n.sample, maf = maf_i, copy_prob = cp)
    geno_i <- hap1 + hap2

    eaf_emp <- colMeans(geno_i) / 2
    mono <- which(!is.finite(eaf_emp) | eaf_emp <= 0 | eaf_emp >= 1)
    if (length(mono)) {
      for (j in mono) {
        ix <- sample.int(n.sample, size = 2L, replace = FALSE)
        geno_i[ix, j] <- c(0L, 2L)
      }
      eaf_emp <- colMeans(geno_i) / 2
    }

    alleles <- .gcanvas_toy_sample_alleles(causal.nsnp_eff)
    local_col_idx <- seq.int(col_cursor + 1L, col_cursor + causal.nsnp_eff)
    col_cursor <- col_cursor + causal.nsnp_eff
    local_G_list[[i]] <- geno_i

    w <- exp(-abs(pos_i - center_i) / max(1, as.numeric(causal.window) * 0.2))
    k_i <- min(causal_nvar[i], causal.nsnp_eff)
    c_local <- sort(sample(seq_len(causal.nsnp_eff), size = k_i, replace = FALSE, prob = w))
    c_global <- local_col_idx[c_local]
    causal_cols <- c(causal_cols, c_global)
    causal_locus_of_causal <- c(causal_locus_of_causal, rep(i, length(c_global)))

    eff_sd <- if (k_i == 1L) {
      stats::runif(1, min = 0.030, max = 0.055)
    } else if (k_i == 2L) {
      stats::runif(1, min = 0.020, max = 0.040)
    } else {
      stats::runif(1, min = 0.015, max = 0.030)
    }
    eff <- stats::rnorm(length(c_global), mean = 0, sd = eff_sd)
    if (!length(eff) || !any(abs(eff) > 1e-8)) eff <- rep(eff_sd, length(c_global))
    beta_base <- c(beta_base, eff)

    local_dt_list[[i]] <- data.table::data.table(
      locus_id = as_int(i),
      CHR = as.character(chr_i),
      POS = as.numeric(pos_i),
      A1 = alleles$A1,
      A2 = alleles$A2,
      EAF = as.numeric(eaf_emp),
      center_pos = as.numeric(center_i),
      local_col_idx = as_int(local_col_idx),
      is_local = TRUE,
      is_causal = local_col_idx %in% c_global,
      true_beta = 0
    )
    .toy_gwas_progress_set(pb_locus, i)
  }
  .toy_gwas_progress_close(pb_locus)

  var_local <- data.table::rbindlist(local_dt_list, use.names = TRUE, fill = TRUE)
  if (!nrow(var_local)) stop("No local variants generated.", call. = FALSE)
  if (!length(local_G_list)) stop("No local genotype matrix generated.", call. = FALSE)
  G_local <- do.call(cbind, local_G_list)
  mode(G_local) <- "numeric"

  var_local[match(causal_cols, local_col_idx), true_beta := beta_base]
  locus_col_list <- split(var_local$local_col_idx, var_local$locus_id)

  m_local <- as_int(ncol(G_local))
  has_ld_matrix_step <- isTRUE(m_local > 1L)
  total_steps <- 2L + as_int(has_ld_matrix_step) + as_int(isTRUE(bopt$write_bfile))
  local_ord <- order(var_local$local_col_idx)
  if (!all(var_local$local_col_idx[local_ord] == seq_len(m_local))) {
    stop("Internal error: local_col_idx is not aligned to local genotype matrix columns.", call. = FALSE)
  }
  local_info <- var_local[local_ord]

  if (has_ld_matrix_step) {
    .gcanvas_note("gcanvas::toy.gwas", sprintf("2/%d building LD matrix (local m=%d)", as_int(total_steps), as_int(m_local)), silent = !show_progress)
    R_local <- suppressWarnings(stats::cor(G_local))
    R_local[!is.finite(R_local)] <- 0
    R_local <- (R_local + t(R_local)) / 2
    diag(R_local) <- 1
    eig <- eigen(R_local, symmetric = TRUE)
    eig$values[eig$values < 1e-8] <- 1e-8
    R_local <- eig$vectors %*% (eig$values * t(eig$vectors))
    R_local <- (R_local + t(R_local)) / 2
    diag(R_local) <- 1
    chol_R <- suppressWarnings(chol(R_local))
    z_noise <- as.numeric(t(chol_R) %*% stats::rnorm(m_local))
  } else {
    R_local <- matrix(1, nrow = 1L, ncol = 1L)
    z_noise <- stats::rnorm(1L)
  }

  var_g_local <- pmax(2 * as.numeric(local_info$EAF) * (1 - as.numeric(local_info$EAF)), 1e-8)
  se_local <- 1 / sqrt(n.sample * var_g_local)
  effect_scale_h2 <- sqrt(max(1e-6, h2_target / 0.2))

  assoc_once <- function(locus_scale) {
    b <- numeric(m_local)
    eff <- beta_base * effect_scale_h2 * locus_scale[causal_locus_of_causal]
    b[causal_cols] <- eff
    mu_z <- sqrt(n.sample) * as.numeric(R_local %*% b)
    z <- mu_z + z_noise
    p <- pmin(pmax(2 * stats::pnorm(-abs(z)), p.floor), 1)
    beta <- z * se_local
    list(p = p, beta = beta, se = se_local, z = z)
  }

  assoc_step <- if (has_ld_matrix_step) 3L else 2L
  .gcanvas_note("gcanvas::toy.gwas", sprintf("%d/%d computing association (local m=%d)", as_int(assoc_step), as_int(total_steps), as_int(m_local)), silent = !show_progress)
  scales <- rep(1, n.causal)
  assoc <- assoc_once(scales)
  forced_loci <- integer(0)
  for (iter in seq_len(15L)) {
    locus_sig <- vapply(seq_len(n.causal), function(i) {
      idx <- locus_col_list[[as.character(i)]]
      pi <- assoc$p[idx]
      pi <- pi[is.finite(pi) & !is.na(pi)]
      if (!length(pi)) return(0L)
      as_int(sum(pi < p.threshold, na.rm = TRUE))
    }, integer(1))
    fail <- which(locus_sig < sig_target)
    if (!length(fail)) break
    scales[fail] <- pmin(scales[fail] * 1.15, 2.5)
    assoc <- assoc_once(scales)
  }

  locus_sig <- vapply(seq_len(n.causal), function(i) {
    idx <- locus_col_list[[as.character(i)]]
    pi <- assoc$p[idx]
    pi <- pi[is.finite(pi) & !is.na(pi)]
    if (!length(pi)) return(0L)
    as_int(sum(pi < p.threshold, na.rm = TRUE))
  }, integer(1))
  fail <- which(locus_sig < sig_target)
  if (length(fail)) {
    for (i in fail) {
      idx <- locus_col_list[[as.character(i)]]
      idx <- idx[is.finite(assoc$p[idx]) & !is.na(assoc$p[idx])]
      if (!length(idx)) next
      sig_now <- as_int(sum(assoc$p[idx] < p.threshold, na.rm = TRUE))
      need <- sig_target - sig_now
      if (need <= 0L) next
      ord <- idx[order(assoc$p[idx])]
      cand <- ord[assoc$p[ord] >= p.threshold]
      pick <- cand[seq_len(min(need, length(cand)))]
      if (!length(pick)) next
      p_forced <- p.threshold * stats::runif(length(pick), min = 1e-3, max = 0.8)
      p_forced <- pmax(p_forced, p.floor)
      z_abs <- stats::qnorm(p_forced / 2, lower.tail = FALSE)
      sgn <- sign(assoc$beta[pick])
      sgn[!is.finite(sgn) | sgn == 0] <- sample(c(-1, 1), size = sum(!is.finite(sgn) | sgn == 0), replace = TRUE)
      se0 <- assoc$se[pick]
      se0[!is.finite(se0) | se0 <= 0] <- 1 / sqrt(max(10, n.sample))
      assoc$p[pick] <- p_forced
      assoc$z[pick] <- sgn * z_abs
      assoc$beta[pick] <- sgn * z_abs * se0
      assoc$se[pick] <- se0
      forced_loci <- c(forced_loci, i)
    }
  }

  map_idx <- match(var_local$local_col_idx, seq_len(m_local))
  var_local[, `:=`(
    P = assoc$p[map_idx],
    BETA = assoc$beta[map_idx],
    SE = assoc$se[map_idx],
    Z = assoc$z[map_idx]
  )]

  windows <- centers[, .(
    CHR = as.character(CHR),
    start = pmax(1, center_pos - as.numeric(causal.window)),
    end = center_pos + as.numeric(causal.window)
  )]
  bg_pos <- .toy_gwas_sample_bg_positions(n = n_bg, chr_bounds = chr_bounds, block_windows = windows)

  if (nrow(bg_pos)) {
    bg_alleles <- .gcanvas_toy_sample_alleles(nrow(bg_pos))
    eaf_bg <- pmin(pmax(stats::rbeta(nrow(bg_pos), shape1 = 0.9, shape2 = 5), eaf.range[1]), eaf.range[2])
    var_g_bg <- pmax(2 * eaf_bg * (1 - eaf_bg), 1e-8)
    se_bg <- 1 / sqrt(n.sample * var_g_bg)
    z_bg <- stats::rnorm(nrow(bg_pos))
    beta_bg <- z_bg * se_bg
    p_bg <- pmin(pmax(2 * stats::pnorm(-abs(z_bg)), p.floor), 1)
    var_bg <- data.table::data.table(
      locus_id = 0L,
      CHR = as.character(bg_pos$CHR),
      POS = as.numeric(bg_pos$POS),
      A1 = bg_alleles$A1,
      A2 = bg_alleles$A2,
      EAF = as.numeric(eaf_bg),
      center_pos = NA_real_,
      local_col_idx = NA_integer_,
      is_local = FALSE,
      is_causal = FALSE,
      true_beta = 0,
      P = p_bg,
      BETA = beta_bg,
      SE = se_bg,
      Z = z_bg
    )
  } else {
    var_bg <- data.table::data.table(
      locus_id = integer(), CHR = character(), POS = numeric(),
      A1 = character(), A2 = character(), EAF = numeric(), center_pos = numeric(),
      local_col_idx = integer(), is_local = logical(), is_causal = logical(),
      true_beta = numeric(), P = numeric(), BETA = numeric(), SE = numeric(), Z = numeric()
    )
  }

  var_dt <- data.table::rbindlist(list(var_local, var_bg), use.names = TRUE, fill = TRUE)
  var_dt[, `:=`(
    P = pmin(pmax(suppressWarnings(as.numeric(P)), p.floor), 1),
    BETA = suppressWarnings(as.numeric(BETA)),
    SE = suppressWarnings(as.numeric(SE)),
    Z = suppressWarnings(as.numeric(Z)),
    N = as_int(n.sample)
  )]
  var_dt[, negLog10P := -log10(pmax(P, p.floor))]
  var_dt[!is.finite(negLog10P) | is.na(negLog10P), negLog10P := -log10(p.floor)]
  var_dt[, SNP := paste0(CHR, ":", as_int(POS), ":", A1, ":", A2)]

  if (nrow(forced_vars)) {
    forced_snp <- as.character(forced_vars$SNP)
    present_force <- unique(as.character(var_dt$SNP))
    miss_force <- forced_vars[!(SNP %in% present_force)]
    if (nrow(miss_force)) {
      n_add <- as_int(nrow(miss_force))
      drop_idx <- var_dt[locus_id == 0L & !(SNP %in% forced_snp), .I]
      if (length(drop_idx)) {
        nd <- min(as_int(length(drop_idx)), n_add)
        pick_drop <- sample(drop_idx, size = nd)
        var_dt <- var_dt[-pick_drop]
      }
      eaf_f <- pmin(pmax(stats::rbeta(n_add, shape1 = 0.9, shape2 = 5), eaf.range[1]), eaf.range[2])
      vg_f <- pmax(2 * eaf_f * (1 - eaf_f), 1e-8)
      se_f <- 1 / sqrt(as.numeric(n.sample) * vg_f)
      z_f <- stats::rnorm(n_add)
      p_f <- pmin(pmax(2 * stats::pnorm(-abs(z_f)), p.floor), 1)
      beta_f <- z_f * se_f
      add_dt <- data.table::data.table(
        SNP = as.character(miss_force$SNP),
        CHR = as.character(miss_force$CHR),
        POS = as.numeric(miss_force$POS),
        A1 = as.character(miss_force$A1),
        A2 = as.character(miss_force$A2),
        EAF = as.numeric(eaf_f),
        P = as.numeric(p_f),
        negLog10P = as.numeric(-log10(pmax(p_f, p.floor))),
        BETA = as.numeric(beta_f),
        SE = as.numeric(se_f),
        Z = as.numeric(z_f),
        N = as_int(n.sample),
        locus_id = 0L,
        is_local = FALSE,
        is_causal = FALSE,
        is_lead = FALSE,
        center_pos = NA_real_,
        true_beta = 0,
        local_col_idx = NA_integer_
      )
      var_dt <- data.table::rbindlist(list(var_dt, add_dt), use.names = TRUE, fill = TRUE)
    }
    if (nrow(var_dt) > as_int(n.snp)) {
      need_drop <- as_int(nrow(var_dt) - as_int(n.snp))
      drop_bg <- var_dt[!(SNP %in% forced_snp) & locus_id == 0L, .I]
      drop_pool <- if (length(drop_bg) >= need_drop) {
        drop_bg
      } else {
        c(drop_bg, var_dt[!(SNP %in% forced_snp) & !(.I %in% drop_bg), .I])
      }
      if (length(drop_pool)) {
        ord <- drop_pool[order(var_dt$P[drop_pool], decreasing = TRUE, na.last = TRUE)]
        var_dt <- var_dt[-ord[seq_len(min(need_drop, length(ord)))]]
      }
    }
  }

  var_dt[, is_lead := FALSE]
  lead_row <- var_dt[locus_id > 0, .I[which.min(P)], by = locus_id]$V1
  lead_row <- lead_row[is.finite(lead_row)]
  if (length(lead_row)) var_dt[lead_row, is_lead := TRUE]

  var_dt[, chr_order := rank.chrom(CHR)]
  data.table::setorder(var_dt, chr_order, POS, SNP)
  var_dt[, chr_order := NULL]
  data.table::setcolorder(var_dt, c(
    "SNP", "CHR", "POS", "A1", "A2", "EAF",
    "P", "negLog10P", "BETA", "SE", "Z", "N",
    "locus_id", "is_local", "is_causal", "is_lead", "center_pos", "true_beta", "local_col_idx"
  ))

  bfile_prefix <- NULL
  if (isTRUE(bopt$write_bfile)) {
    bfile_step <- if (has_ld_matrix_step) 4L else 3L
    .gcanvas_note("gcanvas::toy.gwas", sprintf("%d/%d writing bfile (n_sample=%d, n_snp=%d)", as_int(bfile_step), as_int(total_steps), as_int(n.sample), as_int(nrow(var_dt))), silent = !show_progress)
    bfile_prefix <- .gcanvas_pretty_path(abs_path(bopt$bfile_prefix))
    dir.create(dirname(bfile_prefix), recursive = TRUE, showWarnings = FALSE)
    targets <- paste0(bfile_prefix, c(".bed", ".bim", ".fam"))
    if (any(file.exists(targets))) {
      if (isTRUE(overwrite)) {
        unlink(targets, force = TRUE)
      } else {
        stop("bfile already exists. Set overwrite=TRUE or provide another bfile prefix.", call. = FALSE)
      }
    }

    n_total <- as_int(nrow(var_dt))
    if ((as.numeric(n_total) * as.numeric(n.sample)) > 3e8) {
      stop("Requested bfile is too large. Reduce n.snp/n.sample or set bfile=FALSE.", call. = FALSE)
    }

    fbm_base <- tempfile(pattern = "toy_gwas_fbm_")
    G_fbm <- bigstatsr::FBM.code256(
      nrow = n.sample, ncol = n_total,
      code = c(0, 1, 2, 3, rep(NA_real_, 252)),
      init = 0, backingfile = fbm_base
    )

    chunk <- 600L
    n_chunk <- as_int(ceiling(as.numeric(n_total) / as.numeric(chunk)))
    pb_bed <- .toy_gwas_progress_new(n_chunk, "writing PLINK chunks")
    chunk_i <- 0L
    for (from in seq.int(1L, n_total, by = chunk)) {
      to <- min(n_total, from + chunk - 1L)
      idx <- from:to
      dtc <- var_dt[idx]
      block <- matrix(0L, nrow = n.sample, ncol = length(idx))
      loc <- which(!is.na(dtc$is_local) & dtc$is_local & is.finite(dtc$local_col_idx))
      if (length(loc)) {
        lc <- as_int(dtc$local_col_idx[loc])
        block[, loc] <- G_local[, lc, drop = FALSE]
      }
      bg <- setdiff(seq_along(idx), loc)
      if (length(bg)) {
        p_bg <- pmin(pmax(as.numeric(dtc$EAF[bg]), 1e-4), 0.999)
        block[, bg] <- matrix(
          stats::rbinom(n.sample * length(bg), size = 2L, prob = rep(p_bg, each = n.sample)),
          nrow = n.sample, ncol = length(bg)
        )
      }
      G_fbm[, idx] <- block
      chunk_i <- chunk_i + 1L
      .toy_gwas_progress_set(pb_bed, chunk_i)
    }
    .toy_gwas_progress_close(pb_bed)

    fam <- data.frame(
      family.ID = rep(1, n.sample),
      sample.ID = sprintf("toy_%06d", seq_len(n.sample)),
      paternal.ID = 0,
      maternal.ID = 0,
      sex = 0,
      affection = -9,
      stringsAsFactors = FALSE
    )
    fam$family.ID <- .gcanvas_toy_apply_iid_prefix(fam$family.ID, iid.prefix = iid_prefix)
    fam$sample.ID <- .gcanvas_toy_apply_iid_prefix(fam$sample.ID, iid.prefix = iid_prefix)
    map <- data.frame(
      chromosome = as_int(rank.chrom(var_dt$CHR)),
      marker.ID = var_dt$SNP,
      genetic.dist = 0,
      physical.pos = as_int(var_dt$POS),
      allele1 = var_dt$A1,
      allele2 = var_dt$A2,
      stringsAsFactors = FALSE
    )
    obj <- structure(list(genotypes = G_fbm, fam = fam, map = map), class = "bigSNP")
    bigsnpr::snp_writeBed(obj, paste0(bfile_prefix, ".bed"))
    .gcanvas_note("gcanvas::toy.gwas", paste0("bfile written: ", bfile_prefix))

    tmp_files <- c(fbm_base, paste0(fbm_base, ".bk"), paste0(fbm_base, ".rds"))
    unlink(tmp_files[file.exists(tmp_files)], force = TRUE)
  }

  locus_meta <- var_dt[locus_id > 0, .(
    CHR = CHR[1],
    center_pos = center_pos[1],
    n_variants = .N,
    n_causal = sum(is_causal, na.rm = TRUE),
    n_sig = sum(P < p.threshold, na.rm = TRUE),
    lead_snp = SNP[which.min(P)],
    lead_pos = POS[which.min(P)],
    lead_p = min(P, na.rm = TRUE)
  ), by = locus_id]
  locus_meta <- merge(locus_meta, centers[, .(locus_id, n_causal_var)], by = "locus_id", all.x = TRUE, sort = TRUE)
  locus_meta[, forced_threshold := locus_id %in% unique(forced_loci)]

  attr(var_dt, "gcanvas_meta") <- list(
    type = "toy.gwas",
    build = as_int(build),
    p_threshold = as.numeric(p.threshold),
    n_variants = as_int(nrow(var_dt)),
    n_samples = as_int(n.sample),
    n_causal_loci = as_int(n.causal),
    bfile = bfile_prefix,
      simulation = list(
        n_snp = as_int(n.snp),
        n_local = as_int(local_total),
        n_background = as_int(n_bg),
        forced_variants = as.character(forced_vars$SNP),
        causal_window = as_int(causal.window),
        causal_nsnp = as_int(causal.nsnp_eff),
      causal_min_sig = as_int(causal.min.sig),
      causal_mix = mix,
      h2 = as.numeric(h2_target),
      ld_decay_bp = as.numeric(ld.decay.bp),
      p_floor = as.numeric(p.floor),
      eaf_range = eaf.range
    ),
    causal_loci = locus_meta
  )
  var_dt
}

#' Simulate a toy eQTL summary-statistics table
#'
#' Generates a synthetic eQTL dataset compatible with the [regional()] /
#' [manhattan()] / [calcld()] surface. Adds gene-level columns
#' (`gene_id`, `gene_name`, `gene_tss`, `gene_locus_id`, `gene_is_focal`,
#' `is_lead_gene`, ...) on top of the standard GWAS layout.
#'
#' Can either build a new variant scaffold from scratch or piggyback on an
#' existing [toy.gwas()] table via `base = <toy.gwas result>` so the eQTL
#' and GWAS share variants and bfile.
#'
#' @param base Optional output of [toy.gwas()] to reuse the variant scaffold.
#' @param n.sample Sample size used to scale standard errors.
#' @param n.snp,n.causal,causal.window,causal.nsnp,causal.min.sig Causal-
#'   locus geometry.
#' @param iid.prefix,build,chr,p.threshold,p.floor,eaf.range,h2,causal.mix,ld.decay.bp
#'   Same meaning as in [toy.gwas()].
#' @param gene.per.locus Number of genes drawn per causal locus.
#' @param gene.prefix Prefix used to synthesize Ensembl-style gene ids.
#' @param gtf Optional path to a bgzipped + tabix-indexed GTF -- if supplied,
#'   gene names / TSS are pulled from real annotation instead of being made up.
#' @param gtf.rds Optional `gtf2rds()` cache (alternative to `gtf`).
#' @param cis.window,cis.decay.bp,focal.cis.scale,nonfocal.cis.scale,focal.causal.boost,nonfocal.causal.boost,trans.noise.ratio
#'   Cis-window geometry and effect-size mixing for focal vs. non-focal genes.
#' @param bfile Logical / path. See [toy.gwas()] for behavior.
#' @param overwrite,show.progress,seed Same meaning as in [toy.gwas()].
#' @param silent Logical. Suppress progress notes.
#'
#' @return A `data.table` of eQTL summary statistics with one row per
#'   (gene, variant) pair. Carries `attr(., "gcanvas_meta")` with run params.
#'
#' @seealso [toy.gwas()] for the matching GWAS-style generator.
#' @export
toy.eqtl <- function(base = NULL,
                     n.sample = 500L,
                     n.snp = 10000L,
                     n.causal = 3L,
                     iid.prefix = NULL,
                     build = 38L,
                     chr = 1:22,
                     causal.window = 2e5L,
                     causal.nsnp = 120L,
                     causal.min.sig = 12L,
                     p.threshold = 5e-8,
                     p.floor = 1e-50,
                     eaf.range = c(0.01, 0.5),
                     h2 = 0.35,
                     causal.mix = c(0.55, 0.30, 0.15),
                     ld.decay.bp = 3e4,
                     gene.per.locus = 3L,
                     gene.prefix = "ENSGTOY",
                     gtf = NULL,
                     gtf.rds = NULL,
                     cis.window = 1e6L,
                     cis.decay.bp = 5e4L,
                     focal.cis.scale = c(0.35, 0.65),
                     nonfocal.cis.scale = c(0.08, 0.28),
                     focal.causal.boost = c(2.0, 4.0),
                     nonfocal.causal.boost = c(0.3, 1.1),
                     trans.noise.ratio = 0.15,
                     bfile = FALSE,
                     overwrite = TRUE,
                     show.progress = TRUE,
                     seed = NULL,
                     silent = FALSE) {
  require_pkg("data.table")
  silent <- isTRUE(silent)

  seed_use <- .gcanvas_seed_resolve(seed)
  if (!is.null(seed_use)) set.seed(seed_use)
  seed_label <- .gcanvas_seed_label(seed = seed, seed_use = seed_use)
  .gcanvas_note(
    "gcanvas::toy.eqtl",
    sprintf(
      "Start: base=%s | n.sample=%d | n.snp=%d | n.causal=%d | build=%s | seed=%s",
      ifelse(is.null(base), "auto", "provided"),
      as_int(n.sample), as_int(n.snp), as_int(n.causal), as.character(build)[1], seed_label
    ),
    silent = silent
  )

  p.threshold <- suppressWarnings(as.numeric(p.threshold))[1]
  if (!is.finite(p.threshold) || p.threshold <= 0 || p.threshold >= 1) p.threshold <- 5e-8
  p.floor <- suppressWarnings(as.numeric(p.floor))[1]
  if (!is.finite(p.floor) || p.floor <= 0 || p.floor >= p.threshold) p.floor <- 1e-50
  gene.per.locus <- as_int(gene.per.locus)
  if (!is.finite(gene.per.locus) || gene.per.locus < 1L) gene.per.locus <- 1L
  cis.window <- as_int(cis.window)
  if (!is.finite(cis.window) || cis.window < 1L) cis.window <- 1e6L
  cis.decay.bp <- as_int(cis.decay.bp)
  if (!is.finite(cis.decay.bp) || cis.decay.bp < 500L) cis.decay.bp <- 5e4L
  trans.noise.ratio <- suppressWarnings(as.numeric(trans.noise.ratio))[1]
  if (!is.finite(trans.noise.ratio)) trans.noise.ratio <- 0.15
  trans.noise.ratio <- max(0, min(0.95, trans.noise.ratio))

  normalize_pair <- function(x, default = c(0.1, 0.3)) {
    v <- suppressWarnings(as.numeric(x))
    if (length(v) < 2L || any(!is.finite(v[1:2]))) v <- default
    v <- sort(v[1:2])
    c(v[1], max(v[1], v[2]))
  }
  focal.cis.scale <- normalize_pair(focal.cis.scale, c(0.35, 0.65))
  nonfocal.cis.scale <- normalize_pair(nonfocal.cis.scale, c(0.08, 0.28))
  focal.causal.boost <- normalize_pair(focal.causal.boost, c(2.0, 4.0))
  nonfocal.causal.boost <- normalize_pair(nonfocal.causal.boost, c(0.3, 1.1))

  generated_base <- FALSE
  if (is.null(base)) {
    generated_base <- TRUE
    .gcanvas_note(
      "gcanvas::toy.eqtl",
      sprintf("Generating base toy.gwas (show.progress=%s)", ifelse(isTRUE(show.progress), "TRUE", "FALSE")),
      silent = silent
    )
    base <- toy.gwas(
      n.sample = n.sample,
      n.snp = n.snp,
      n.causal = n.causal,
      build = build,
      chr = chr,
      causal.window = causal.window,
      causal.nsnp = causal.nsnp,
      causal.min.sig = causal.min.sig,
      p.threshold = p.threshold,
      p.floor = p.floor,
      eaf.range = eaf.range,
      h2 = h2,
      causal.mix = causal.mix,
      ld.decay.bp = ld.decay.bp,
      bfile = bfile,
      overwrite = overwrite,
      show.progress = show.progress,
      iid.prefix = iid.prefix,
      seed = seed_use %||% seed
    )
  } else {
    if (!(is.data.frame(base) || data.table::is.data.table(base))) {
      stop("base must be a toy.gwas-like data.frame/data.table.", call. = FALSE)
    }
  }

  meta0 <- attr(base, "gcanvas_meta")
  if (!isTRUE(generated_base)) {
    bmeta <- as_int(meta0$build %||% NA_integer_)
    if (is.finite(bmeta)) build <- bmeta
  }
  dt0 <- if (data.table::is.data.table(base)) data.table::copy(base) else data.table::as.data.table(base)
  req <- c("SNP", "CHR", "POS", "A1", "A2", "EAF")
  miss <- setdiff(req, names(dt0))
  if (length(miss)) stop("base is missing columns: ", paste(miss, collapse = ", "), call. = FALSE)
  dt0[, `:=`(
    CHR = normalize.chrom(CHR),
    POS = suppressWarnings(as.numeric(POS)),
    EAF = pmin(pmax(suppressWarnings(as.numeric(EAF)), 1e-4), 0.999)
  )]
  dt0 <- dt0[!is.na(CHR) & nzchar(CHR) & is.finite(POS)]
  if (!nrow(dt0)) stop("base has no valid variants.", call. = FALSE)

  if (!("locus_id" %in% names(dt0))) dt0[, locus_id := 0L]
  dt0[, locus_id := as_int(locus_id)]
  dt0[is.na(locus_id), locus_id := 0L]
  if (!("is_causal" %in% names(dt0))) dt0[, is_causal := FALSE]
  dt0[, is_causal := as.logical(is_causal)]
  dt0[is.na(is_causal), is_causal := FALSE]
  if (!("center_pos" %in% names(dt0))) dt0[, center_pos := NA_real_]
  dt0[, center_pos := suppressWarnings(as.numeric(center_pos))]

  if (!("SE" %in% names(dt0))) {
    se0 <- 1 / sqrt(max(1, as.numeric(n.sample)) * pmax(2 * dt0$EAF * (1 - dt0$EAF), 1e-8))
    dt0[, SE := as.numeric(se0)]
  } else {
    dt0[, SE := suppressWarnings(as.numeric(SE))]
    bad_se <- !is.finite(dt0$SE) | dt0$SE <= 0
    if (any(bad_se)) {
      dt0[bad_se, SE := 1 / sqrt(max(1, as.numeric(n.sample)) * pmax(2 * EAF * (1 - EAF), 1e-8))]
    }
  }
  if (!("Z" %in% names(dt0))) {
    if ("BETA" %in% names(dt0)) {
      dt0[, Z := suppressWarnings(as.numeric(BETA)) / pmax(SE, 1e-12)]
    } else if ("P" %in% names(dt0)) {
      p0 <- pmin(pmax(suppressWarnings(as.numeric(dt0$P)), p.floor), 1)
      z0 <- stats::qnorm(p0 / 2, lower.tail = FALSE)
      if ("BETA" %in% names(dt0)) {
        sgn <- sign(suppressWarnings(as.numeric(dt0$BETA)))
        sgn[!is.finite(sgn) | sgn == 0] <- sample(c(-1, 1), size = sum(!is.finite(sgn) | sgn == 0), replace = TRUE)
      } else {
        sgn <- sample(c(-1, 1), size = nrow(dt0), replace = TRUE)
      }
      dt0[, Z := as.numeric(sgn * z0)]
    } else {
      stop("base must include Z, or BETA+SE, or P (+/- BETA sign).", call. = FALSE)
    }
  } else {
    dt0[, Z := suppressWarnings(as.numeric(Z))]
  }
  dt0[!is.finite(Z), Z := 0]

  if (!("P" %in% names(dt0))) {
    dt0[, P := pmin(pmax(2 * stats::pnorm(-abs(Z)), p.floor), 1)]
  } else {
    dt0[, P := pmin(pmax(suppressWarnings(as.numeric(P)), p.floor), 1)]
  }
  if (!("BETA" %in% names(dt0))) dt0[, BETA := as.numeric(Z * SE)]
  if (!("N" %in% names(dt0))) dt0[, N := as_int(n.sample)]

  loci_meta <- data.table::as.data.table(meta0$causal_loci %||% data.frame())
  if (!nrow(loci_meta)) {
    loci_meta <- dt0[locus_id > 0, .(
      CHR = CHR[1],
      center_pos = if (any(is.finite(center_pos))) center_pos[which.max(is.finite(center_pos))] else stats::median(POS),
      lead_pos = POS[which.min(P)],
      lead_snp = SNP[which.min(P)]
    ), by = .(locus_id)]
  } else {
    loci_meta <- loci_meta[, .(
      locus_id = as_int(locus_id),
      CHR = normalize.chrom(CHR),
      center_pos = suppressWarnings(as.numeric(center_pos)),
      lead_pos = suppressWarnings(as.numeric(lead_pos)),
      lead_snp = as.character(lead_snp)
    )]
    fill <- dt0[locus_id > 0, .(
      CHR0 = CHR[1],
      center0 = if (any(is.finite(center_pos))) center_pos[which.max(is.finite(center_pos))] else stats::median(POS),
      lead0 = POS[which.min(P)],
      snp0 = SNP[which.min(P)]
    ), by = .(locus_id)]
    loci_meta <- merge(loci_meta, fill, by = "locus_id", all = TRUE)
    loci_meta[is.na(CHR) | !nzchar(CHR), CHR := CHR0]
    loci_meta[!is.finite(center_pos), center_pos := center0]
    loci_meta[!is.finite(lead_pos), lead_pos := lead0]
    loci_meta[is.na(lead_snp) | !nzchar(lead_snp), lead_snp := snp0]
    loci_meta[, c("CHR0", "center0", "lead0", "snp0") := NULL]
  }
  loci_meta <- loci_meta[is.finite(locus_id) & !is.na(CHR) & nzchar(CHR)]
  loci_meta <- unique(loci_meta, by = "locus_id")
  if (!nrow(loci_meta)) stop("toy.eqtl failed: no causal loci resolved from base.", call. = FALSE)

  if (!any(dt0$locus_id > 0)) {
    cwin <- as.numeric(meta0$simulation$causal_window %||% causal.window)
    if (!is.finite(cwin) || cwin < 1000) cwin <- as.numeric(causal.window)
    dt0[, locus_id := 0L]
    for (i in seq_len(nrow(loci_meta))) {
      row <- loci_meta[i]
      hit <- dt0$CHR == row$CHR[[1]] &
        dt0$POS >= (row$center_pos[[1]] - cwin) &
        dt0$POS <= (row$center_pos[[1]] + cwin)
      dt0[hit, locus_id := as_int(row$locus_id[[1]])]
    }
  }

  dt_local <- dt0[locus_id > 0]
  if (!nrow(dt_local)) stop("toy.eqtl failed: no locus_id > 0 variants in base.", call. = FALSE)

  gene_source <- "synthetic"
  genes_ref <- NULL
  gtf_bgz <- NULL
  if (!is.null(gtf.rds) && length(gtf.rds) && !is.na(gtf.rds[1]) && nzchar(as.character(gtf.rds[1]))) {
    gtf_rds_path <- abs_path(gtf.rds[1])
    if (!file.exists(gtf_rds_path)) stop("gtf.rds not found: ", gtf_rds_path, call. = FALSE)
    gobj <- readRDS(gtf_rds_path)
    if (is.list(gobj) && !is.null(gobj$gene)) {
      genes_ref <- data.table::as.data.table(gobj$gene)
    } else {
      genes_ref <- data.table::as.data.table(gobj)
    }
    gene_source <- "gtf.rds"
  } else if (!is.null(gtf) && length(gtf) && !is.na(gtf[1]) && nzchar(as.character(gtf[1]))) {
    gtf0 <- as.character(gtf[1])
    if (grepl("\\.rds$", gtf0, ignore.case = TRUE)) {
      gtf_rds_path <- abs_path(gtf0)
      if (!file.exists(gtf_rds_path)) stop("gtf.rds not found: ", gtf_rds_path, call. = FALSE)
      gobj <- readRDS(gtf_rds_path)
      if (is.list(gobj) && !is.null(gobj$gene)) {
        genes_ref <- data.table::as.data.table(gobj$gene)
      } else {
        genes_ref <- data.table::as.data.table(gobj)
      }
      gene_source <- "gtf.rds"
    } else {
      gtf_bgz <- if (grepl("\\.bgz$", gtf0, ignore.case = TRUE)) abs_path(gtf0) else gtf_prepare_tabix(gtf0, sort = "auto", chr_order = "natural")
      gene_source <- "gtf"
    }
  }
  if (!is.null(genes_ref)) {
    need_gene_cols <- c("CHR", "start", "end")
    miss_gene_cols <- setdiff(need_gene_cols, names(genes_ref))
    if (length(miss_gene_cols)) stop("gtf.rds gene table missing columns: ", paste(miss_gene_cols, collapse = ", "), call. = FALSE)
    genes_ref[, `:=`(
      CHR = normalize.chrom(CHR),
      start = suppressWarnings(as.numeric(start)),
      end = suppressWarnings(as.numeric(end)),
      strand = if ("strand" %in% names(genes_ref)) as.character(strand) else "*",
      gene_id = if ("gene_id" %in% names(genes_ref)) .gcanvas_as_char_no_null(gene_id, empty = "") else "",
      gene_name = if ("gene_name" %in% names(genes_ref)) .gcanvas_as_char_no_null(gene_name, empty = "") else "",
      biotype = if ("biotype" %in% names(genes_ref)) as.character(biotype) else NA_character_,
      biotype_raw = if ("biotype_raw" %in% names(genes_ref)) as.character(biotype_raw) else NA_character_
    )]
    genes_ref[, gene_name := .gcanvas_normalize_ensembl_gene_id(gene_name)]
    genes_ref[, gene_id := .gcanvas_normalize_ensembl_gene_id(gene_id)]
    genes_ref <- genes_ref[is.finite(start) & is.finite(end) & end >= start & !is.na(CHR) & nzchar(CHR)]
    genes_ref[, biotype := ifelse(!is.na(biotype) & nzchar(biotype), as.character(biotype), .gcanvas_biotype_group_vec(biotype_raw))]
    genes_ref[!nzchar(gene_name) & nzchar(gene_id), gene_name := gene_id]
    genes_ref[!nzchar(gene_id) & nzchar(gene_name), gene_id := gene_name]
    if (nrow(genes_ref)) {
      genes_ref[, .gene_key := data.table::fifelse(nzchar(gene_id), gene_id, gene_name)]
      genes_ref <- .gcanvas_pick_gene_representative(genes_ref, ".gene_key")
      genes_ref[, .gene_key := NULL]
    }
  }

  .toy_eqtl_gene_priority <- function(biotype, gene_id, gene_name) {
    bt <- tolower(as.character(biotype))
    gid <- as.character(gene_id)
    gnm <- as.character(gene_name)
    is_pc <- !is.na(bt) & bt == "protein_coding"
    is_lnc <- !is.na(bt) & bt == "lncrna"
    has_symbol <- !is.na(gnm) & nzchar(gnm) & !.gcanvas_is_ensembl_gene_id(gnm)
    has_ensgid <- !is.na(gid) & nzchar(gid) & .gcanvas_is_ensembl_gene_id(gid)
    out <- rep(4L, length(bt))
    out[is_lnc] <- 3L
    out[is_lnc & has_symbol] <- 2L
    out[is_lnc & !has_symbol & !has_ensgid] <- 2L
    out[is_pc] <- 1L
    out
  }

  fetch_locus_genes <- function(chr_i, start_i, end_i) {
    if (!is.null(genes_ref) && nrow(genes_ref)) {
      return(data.table::copy(genes_ref[CHR == chr_i & end >= start_i & start <= end_i, .(CHR, gene_id, gene_name, strand, start, end, biotype)]))
    }
    if (!is.null(gtf_bgz) && nzchar(gtf_bgz)) {
      anno <- gtf_query_gene_exon(
        gtf_bgz, chrom = chr_i, start = as_int(start_i), end = as_int(end_i),
        features = c("gene"), keep_biotype = NULL
      )
      gd <- data.table::as.data.table(anno$gene)
      if (!nrow(gd)) return(gd)
      gd[, `:=`(
        CHR = normalize.chrom(chr_i),
        gene_name = .gcanvas_normalize_ensembl_gene_id(.gcanvas_as_char_no_null(gene_name, empty = "")),
        gene_id = .gcanvas_normalize_ensembl_gene_id(.gcanvas_as_char_no_null(gene_id, empty = "")),
        biotype = as.character(biotype),
        strand = as.character(strand),
        start = suppressWarnings(as.numeric(start)),
        end = suppressWarnings(as.numeric(end))
      )]
      gd[!nzchar(gene_name) & nzchar(gene_id), gene_name := gene_id]
      gd[!nzchar(gene_id) & nzchar(gene_name), gene_id := gene_name]
      return(gd[, .(CHR, gene_id, gene_name, strand, start, end, biotype)])
    }
    data.table::data.table(CHR = character(), gene_id = character(), gene_name = character(), strand = character(), start = numeric(), end = numeric(), biotype = character())
  }

  chr_bounds <- .gcanvas_hg_chr_bounds(build = build)
  chr_len_map <- stats::setNames(as.numeric(chr_bounds$end), chr_bounds$CHR)
  gene_rows <- list()
  gidx <- 0L
  for (i in seq_len(nrow(loci_meta))) {
    l <- loci_meta[i]
    chr_i <- as.character(l$CHR[[1]])
    center_i <- as.numeric(l$center_pos[[1]])
    lead_i <- as.numeric(l$lead_pos[[1]])
    if (!is.finite(lead_i)) lead_i <- center_i
    if (!is.finite(center_i)) center_i <- lead_i
    rstart <- max(1, as_int(round(center_i - cis.window)))
    rend <- as_int(round(center_i + cis.window))
    cand <- fetch_locus_genes(chr_i, rstart, rend)

    if (nrow(cand)) {
      cand <- unique(cand, by = c("gene_id", "gene_name", "start", "end", "strand"))
      cand <- cand[!is.na(CHR) & nzchar(CHR)]
      cand[!nzchar(gene_name) & nzchar(gene_id), gene_name := gene_id]
      cand[!nzchar(gene_id) & nzchar(gene_name), gene_id := gene_name]
      cand <- cand[nzchar(gene_id) | nzchar(gene_name)]
      cand[, biotype := as.character(biotype)]
      cand[, tss := ifelse(strand == "-", end, start)]
      cand[!is.finite(tss), tss := (start + end) / 2]
      cand[, dist := abs(tss - lead_i)]
      cand[, gene_prio := .toy_eqtl_gene_priority(biotype, gene_id, gene_name)]
      data.table::setorderv(cand, c("gene_prio", "dist", "start", "end"), c(1L, 1L, 1L, 1L), na.last = TRUE)
      keep_n <- min(gene.per.locus, nrow(cand))
      cand <- cand[seq_len(keep_n)]
      cand[, `:=`(
        locus_id = as_int(l$locus_id[[1]]),
        is_focal = FALSE
      )]
      cand[1, is_focal := TRUE]
      cand <- cand[, .(gene_id, gene_name, locus_id, CHR, tss, is_focal, biotype)]
    } else {
      chr_len <- as.numeric(chr_len_map[[chr_i]])
      if (!is.finite(chr_len)) chr_len <- max(dt_local[CHR == chr_i, POS], na.rm = TRUE)
      if (!is.finite(chr_len) || chr_len < 1) chr_len <- center_i + cis.window
      kmax <- max(1L, gene.per.locus)
      cand <- data.table::data.table(
        gene_id = sprintf("%s%09d", toupper(gene.prefix), seq_len(kmax) + gidx),
        gene_name = sprintf("TOYGENE%03d", seq_len(kmax) + gidx),
        locus_id = as_int(l$locus_id[[1]]),
        CHR = chr_i,
        tss = pmin(pmax(as.numeric(round(stats::rnorm(kmax, mean = center_i, sd = max(2000, cis.window / 3)))), 1), chr_len),
        is_focal = c(TRUE, rep(FALSE, kmax - 1L)),
        biotype = "synthetic"
      )
      gidx <- gidx + kmax
    }
    gene_rows[[length(gene_rows) + 1L]] <- cand
  }
  if (!length(gene_rows)) stop("toy.eqtl failed: no genes constructed.", call. = FALSE)
  genes <- data.table::rbindlist(gene_rows, use.names = TRUE, fill = TRUE)
  genes <- genes[!is.na(CHR) & nzchar(CHR) & is.finite(tss)]
  genes <- genes[!is.na(gene_id) & nzchar(gene_id)]
  if (!nrow(genes)) stop("toy.eqtl failed: no valid genes after filtering.", call. = FALSE)

  eqtl_list <- vector("list", nrow(genes))
  n_out <- 0L
  for (i in seq_len(nrow(genes))) {
    g <- genes[i]
    sub <- dt_local[locus_id == g$locus_id[[1]]]
    if (!nrow(sub)) next
    in_cis <- abs(sub$POS - g$tss[[1]]) <= cis.window
    use <- in_cis | sub$is_causal
    sub <- sub[use]
    if (!nrow(sub)) next

    d <- abs(sub$POS - g$tss[[1]])
    cis_w <- exp(-d / max(500, cis.decay.bp))
    cis_w <- pmax(0, pmin(1, cis_w))
    base_z <- suppressWarnings(as.numeric(sub$Z))
    base_z[!is.finite(base_z)] <- 0

    if (isTRUE(g$is_focal[[1]])) {
      sc <- stats::runif(1, min = focal.cis.scale[1], max = focal.cis.scale[2])
      floor_mult <- 0.18
      noise_total <- stats::runif(1, min = 0.55, max = 0.90)
      boost_rng <- focal.causal.boost
    } else {
      sc <- stats::runif(1, min = nonfocal.cis.scale[1], max = nonfocal.cis.scale[2])
      floor_mult <- 0.03
      noise_total <- stats::runif(1, min = 0.85, max = 1.35)
      boost_rng <- nonfocal.causal.boost
    }

    cis_sd <- noise_total * (1 - trans.noise.ratio)
    trans_sd <- noise_total * trans.noise.ratio
    gene_trans <- stats::rnorm(1, mean = 0, sd = trans_sd)

    z <- base_z * (floor_mult + sc * cis_w) +
      stats::rnorm(nrow(sub), mean = 0, sd = cis_sd) +
      gene_trans +
      stats::rnorm(nrow(sub), mean = 0, sd = trans_sd * 0.35)

    if (any(sub$is_causal, na.rm = TRUE)) {
      ii <- which(sub$is_causal)
      sgn <- sign(base_z[ii])
      bad <- !is.finite(sgn) | sgn == 0
      if (any(bad)) sgn[bad] <- sample(c(-1, 1), size = sum(bad), replace = TRUE)
      z[ii] <- z[ii] + sgn * stats::runif(length(ii), min = boost_rng[1], max = boost_rng[2])
    }
    z <- pmax(pmin(z, 14), -14)

    se <- 1 / sqrt(as.numeric(n.sample) * pmax(2 * sub$EAF * (1 - sub$EAF), 1e-8))
    beta <- z * se
    p <- pmin(pmax(2 * stats::pnorm(-abs(z)), p.floor), 1)

    sub[, `:=`(
      gene_id = as.character(g$gene_id[[1]]),
      gene_name = as.character(g$gene_name[[1]]),
      gene_tss = as.numeric(g$tss[[1]]),
      gene_locus_id = as_int(g$locus_id[[1]]),
      gene_is_focal = isTRUE(g$is_focal[[1]]),
      BETA = as.numeric(beta),
      SE = as.numeric(se),
      Z = as.numeric(z),
      P = as.numeric(p),
      negLog10P = -log10(as.numeric(p)),
      N = as_int(n.sample)
    )]
    n_out <- n_out + 1L
    eqtl_list[[n_out]] <- sub
  }

  if (n_out == 0L) stop("toy.eqtl produced no cis variants around causal loci.", call. = FALSE)
  eqtl <- data.table::rbindlist(eqtl_list[seq_len(n_out)], use.names = TRUE, fill = TRUE)
  if ("is_lead" %in% names(eqtl)) eqtl[, is_lead := NULL]
  eqtl[, is_lead_gene := FALSE]
  lead_gene <- eqtl[, .I[which.min(P)], by = gene_id]$V1
  lead_gene <- lead_gene[is.finite(lead_gene)]
  if (length(lead_gene)) eqtl[lead_gene, is_lead_gene := TRUE]

  eqtl[, chr_order := rank.chrom(CHR)]
  data.table::setorder(eqtl, gene_id, chr_order, POS, SNP)
  eqtl[, chr_order := NULL]
  data.table::setcolorder(eqtl, c(
    "gene_id", "gene_name", "gene_tss", "gene_locus_id", "gene_is_focal",
    "SNP", "CHR", "POS", "A1", "A2", "EAF",
    "P", "negLog10P", "BETA", "SE", "Z", "N",
    "locus_id", "is_local", "is_causal", "is_lead_gene", "center_pos", "true_beta", "local_col_idx"
  ))

  attr(eqtl, "gcanvas_meta") <- list(
    type = "toy.eqtl",
    build = as_int(build),
    n_samples = as_int(n.sample),
    n_variants = as_int(nrow(eqtl)),
    n_genes = as_int(data.table::uniqueN(eqtl$gene_id)),
    cis_window = as_int(cis.window),
    cis_decay_bp = as_int(cis.decay.bp),
    p_floor = as.numeric(p.floor),
    trans_noise_ratio = as.numeric(trans.noise.ratio),
    gene_source = gene_source,
    base_generated = isTRUE(generated_base),
    genes = genes,
    bfile = meta0$bfile %||% NULL,
    source = meta0
  )
  .gcanvas_note(
    "gcanvas::toy.eqtl",
    sprintf(
      "Done: n_genes=%d | n_variants=%d",
      as_int(data.table::uniqueN(eqtl$gene_id)),
      as_int(nrow(eqtl))
    ),
    silent = silent
  )
  eqtl
}
