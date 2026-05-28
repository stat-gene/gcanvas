# Lead-variant detection over summary stats and the `log10c` numeric helper
# used by some plotting code.

#' Identify lead variants in summary statistics
#'
#' Greedy lead-variant selection: scans rows by increasing p-value (or
#' decreasing |Z|), accepts the most significant variant as a lead, then
#' suppresses everything within `window` bp on the same chromosome, and repeats.
#'
#' @param data A `data.frame`/`data.table` of summary statistics.
#' @param snp.col,chrom.col,pos.col Column names for SNP id, chromosome,
#'   and position.
#' @param p.col,z.col,beta.col,se.col Optional statistic columns used by the
#'   ranking method.
#' @param method One of `"p"`, `"z"`, `"beta_se"` — what to rank by.
#' @param flank Suppression window around each lead (bp or e.g. `"500kb"`).
#' @param flank.locus Optional locus-level suppression window (defaults to `flank`).
#' @param p.threshold Maximum p-value for a row to be eligible as a lead.
#' @param silent Logical. Suppress progress notes.
#'
#' @return A `data.table` of lead variants and the rows they suppress.
#' @export
get.lead <- function(data,
                     snp.col = "SNP",
                     chrom.col = "CHR",
                     pos.col = "POS",
                     p.col = "P",
                     z.col = "Z",
                     beta.col = "BETA",
                     se.col = "SE",
                     method = "p",
                     flank = 5e5,
                     flank.locus = flank,
                     p.threshold = 5e-8,
                     silent = FALSE) {
  require_pkg(c("data.table"))
  if (!is.data.frame(data) && !data.table::is.data.table(data)) {
    stop("data must be a data.frame/data.table.", call. = FALSE)
  }
  silent <- isTRUE(silent)
  method <- tolower(as.character(method)[1])
  if (!(method %in% c("p", "z"))) stop("method must be one of: 'p', 'z'.", call. = FALSE)
  snp_col_use <- .gcanvas_resolve_colname(names(data), snp.col, aliases = character(), required = TRUE, arg_label = "snp.col")
  chrom_col_use <- .gcanvas_resolve_colname(names(data), chrom.col, aliases = character(), required = TRUE, arg_label = "chrom.col")
  pos_col_use <- .gcanvas_resolve_colname(names(data), pos.col, aliases = character(), required = TRUE, arg_label = "pos.col")
  req_base <- c(snp_col_use, chrom_col_use, pos_col_use)
  miss_base <- setdiff(req_base, names(data))
  if (length(miss_base)) stop("Missing columns in data: ", paste(miss_base, collapse = ", "), call. = FALSE)
  if (identical(method, "p") && !(p.col %in% names(data))) {
    method <- "z"
    .gcanvas_note("gcanvas::get.lead", sprintf("There is no %s in columns of data. Running with method='z'.", p.col), silent = silent)
  }
  if (identical(method, "p") && !(p.col %in% names(data))) {
    stop("p.col not found for method='p'.", call. = FALSE)
  }
  if (identical(method, "z")) {
    has_z_col <- z.col %in% names(data)
    has_beta_se <- all(c(beta.col, se.col) %in% names(data))
    if (!has_z_col && !has_beta_se) {
      stop("Missing columns in data for method='z': provide z.col or both beta.col and se.col.", call. = FALSE)
    }
  }

  .get_lead_pick_index_chr <- function(pos, p, flank_bp) {
    n <- length(pos)
    if (!n) return(integer())
    ord <- order(p, pos, na.last = TRUE)
    picked <- integer()
    sel_pos <- numeric()
    for (ii in ord) {
      x <- pos[ii]
      if (!length(sel_pos)) {
        picked <- c(picked, ii)
        sel_pos <- x
        next
      }
      sp <- sort(sel_pos)
      k <- findInterval(x, sp)
      d <- Inf
      if (k > 0L) d <- min(d, abs(x - sp[k]))
      if (k < length(sp)) d <- min(d, abs(sp[k + 1L] - x))
      if (d > flank_bp) {
        picked <- c(picked, ii)
        sel_pos <- c(sel_pos, x)
      }
    }
    picked
  }

  flank <- .gcanvas_parse_bp_span(flank, arg_name = "flank")
  flank.locus <- .gcanvas_parse_bp_span(flank.locus, arg_name = "flank.locus")
  if (!is.finite(flank) || flank <= 0) stop("flank must be a positive bp span.", call. = FALSE)
  if (!is.finite(flank.locus) || flank.locus <= 0) stop("flank.locus must be a positive bp span.", call. = FALSE)

  .gcanvas_note("gcanvas::get.lead", "Preparing variants", silent = silent)
  dt <- if (data.table::is.data.table(data)) data.table::copy(data) else data.table::as.data.table(data)
  dt[, snp0 := as.character(get(snp_col_use))]
  dt[, chr0 := normalize.chrom(get(chrom_col_use))]
  dt[, pos0 := suppressWarnings(as.numeric(get(pos_col_use)))]
  dt <- dt[!is.na(snp0) & nzchar(snp0) & !is.na(chr0) & nzchar(chr0) & is.finite(pos0)]
  if (identical(method, "p")) {
    dt[, p0 := suppressWarnings(as.numeric(get(p.col)))]
    dt <- dt[is.finite(p0) & p0 > 0 & p0 <= 1]
  } else {
    if (z.col %in% names(dt)) {
      dt[, z0 := suppressWarnings(as.numeric(get(z.col)))]
      dt <- dt[is.finite(z0)]
    } else {
      dt[, beta0 := suppressWarnings(as.numeric(get(beta.col)))]
      dt[, se0 := suppressWarnings(as.numeric(get(se.col)))]
      dt <- dt[is.finite(beta0) & is.finite(se0) & se0 != 0]
      dt[, z0 := beta0 / se0]
    }
  }
  if (!nrow(dt)) {
    out0 <- data.table::data.table(
      indexSNP = character(), leadSNP = character(), locus = character(),
      CHR = character(), START = numeric(), END = numeric(), n_indexSNP = integer()
    )
    return(out0)
  }
  if (identical(method, "p")) dt[p0 < 1e-317, p0 := 1e-317]
  dt[, chr_order := rank.chrom(chr0)]

  thr_raw <- suppressWarnings(as.numeric(p.threshold))[1]
  use_thr <- is.finite(thr_raw) && !is.na(thr_raw) && thr_raw > 0 && thr_raw < 1
  if (identical(method, "p")) {
    if (use_thr) dt <- dt[p0 < thr_raw]
    dt[, score0 := -log10c(as.character(p0))]
  } else {
    if (!use_thr) stop("p.threshold should be a single numeric value in (0, 1) for method='z'.", call. = FALSE)
    z_cut <- stats::qnorm(log(thr_raw) - log(2), lower.tail = FALSE, log.p = TRUE)
    if (!is.finite(z_cut)) stop("Failed to calculate Z cutoff from p.threshold.", call. = FALSE)
    dt <- dt[abs(z0) > z_cut]
    dt[, score0 := abs(z0)]
  }
  if (!nrow(dt)) {
    out0 <- data.table::data.table(
      indexSNP = character(), leadSNP = character(), locus = character(),
      CHR = character(), START = numeric(), END = numeric(), n_indexSNP = integer()
    )
    return(out0)
  }

  .gcanvas_note("gcanvas::get.lead", "Selecting positional index SNPs", silent = silent)
  idx_rows <- dt[, {
    pick <- .get_lead_pick_index_chr(pos = pos0, p = -score0, flank_bp = flank)
    .(row_idx = .I[pick])
  }, by = .(chr0)]
  idx <- dt[idx_rows$row_idx]
  data.table::setorderv(idx, c("chr_order", "pos0", "score0", "snp0"), c(1L, 1L, -1L, 1L), na.last = TRUE)

  .gcanvas_note("gcanvas::get.lead", "Merging index SNPs into loci", silent = silent)
  idx[, win_start := pmax(1, pos0 - flank.locus)]
  idx[, win_end := pos0 + flank.locus]
  idx[, locus_id := {
    n <- .N
    out <- integer(n)
    if (!n) out else {
      cur <- 1L
      cur_end <- win_end[1L]
      out[1L] <- cur
      if (n >= 2L) {
        for (i in 2:n) {
          if (win_start[i] <= cur_end) {
            out[i] <- cur
            if (win_end[i] > cur_end) cur_end <- win_end[i]
          } else {
            cur <- cur + 1L
            out[i] <- cur
            cur_end <- win_end[i]
          }
        }
      }
      out
    }
  }, by = .(chr0)]

  locus_dt <- idx[, {
    n_idx <- .N
    lead_i <- which.max(score0)[1]
    min_p <- if ("p0" %in% names(.SD)) {
      suppressWarnings(min(as.numeric(p0), na.rm = TRUE))
    } else {
      zabs <- abs(as.numeric(score0))
      zabs <- zabs[is.finite(zabs)]
      if (length(zabs)) suppressWarnings(min(2 * stats::pnorm(zabs, lower.tail = FALSE), na.rm = TRUE)) else NA_real_
    }
    if (!is.finite(min_p)) min_p <- NA_real_
    .(
      leadSNP = snp0[lead_i],
      n_indexSNP = as_int(n_idx),
      minP = as.numeric(min_p),
      START = pmax(1, min(pos0) - flank),
      END = max(pos0) + flank
    )
  }, by = .(chr0, locus_id)]
  locus_dt[, CHR := chr0]
  locus_dt[, `:=`(START = as.numeric(START), END = as.numeric(END))]
  locus_dt[, chrom_locus := ifelse(grepl("^chr", tolower(CHR)), CHR, paste0("chr", CHR))]
  locus_dt[, locus := paste0(chrom_locus, ":", as_int(round(START)), "-", as_int(round(END)))]
  locus_dt[, chrom_locus := NULL]
  out <- merge(
    idx[, .(indexSNP = snp0, chr0, locus_id)],
    locus_dt[, .(chr0, locus_id, leadSNP, locus, CHR, START, END, n_indexSNP, minP)],
    by = c("chr0", "locus_id"),
    all.x = TRUE,
    sort = FALSE
  )
  out[, c("chr0", "locus_id") := NULL]
  out[, locus_size := as.numeric(END - START)]
  data.table::setcolorder(out, c("indexSNP", "leadSNP", "locus", "CHR", "START", "END", "n_indexSNP", "minP", "locus_size"))
  out[, chr_order := rank.chrom(CHR)]
  data.table::setorderv(out, c("chr_order", "CHR", "START", "END", "indexSNP"), c(1L, 1L, 1L, 1L, 1L), na.last = TRUE)
  out[, chr_order := NULL]
  unique(out)
}

log10c <- function(x) {
  name.ori <- names(x)
  x <- toupper(as.character(x)); names(x) <- seq_along(x)
  idx.exp <- grep("E-", x)
  if (!length(idx.exp)) {
    res <- log10(as.numeric(x)); names(res) <- name.ori; return(res)
  }
  x.exp <- x[idx.exp]; x.non <- x[-idx.exp]
  name.exp <- names(x.exp); name.non <- names(x.non)
  parts <- strsplit(x.exp, "E-")
  mat <- matrix(as.numeric(unlist(parts)), ncol = 2, byrow = TRUE)
  v <- log10(mat[,1]) - mat[,2]; names(v) <- name.exp
  w <- log10(as.numeric(x.non)); names(w) <- name.non
  out <- c(v, w)
  out <- out[order(as.numeric(names(out)))]
  names(out) <- name.ori
  out
}

.gcanvas_normalize_highlight_map <- function(genes, colors = NULL, default_color = "#9307E8") {
  genes <- unique(.gcanvas_normalize_ensembl_gene_id(as.character(genes)))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  if (!length(genes)) return(NULL)

  if (is.null(colors) || !length(colors)) {
    out <- rep(default_color, length(genes))
    names(out) <- genes
    return(out)
  }

  if (is.list(colors)) colors <- unlist(colors, use.names = TRUE)
  cols <- as.character(colors)
  nm <- names(cols)
  if (!is.null(nm) && any(nzchar(nm))) {
    nm <- .gcanvas_normalize_ensembl_gene_id(nm)
    names(cols) <- nm
    out <- rep(default_color, length(genes))
    names(out) <- genes
    hit <- intersect(genes, nm[nzchar(nm)])
    out[hit] <- cols[hit]
    return(out)
  }

  out <- rep(cols, length.out = length(genes))
  names(out) <- genes
  out
}

.gcanvas_pick_hl_color <- function(map, gene_name, gene_id) {
  if (is.null(map) || !length(map)) return(NA_character_)
  gene_name <- .gcanvas_normalize_ensembl_gene_id(gene_name)
  gene_id <- .gcanvas_normalize_ensembl_gene_id(gene_id)
  if (!is.na(gene_name) && nzchar(gene_name) && gene_name %in% names(map)) return(unname(map[gene_name]))
  if (!is.na(gene_id) && nzchar(gene_id) && gene_id %in% names(map)) return(unname(map[gene_id]))
  NA_character_
}

.gcanvas_resolve_direction_a1 <- function(direction.a1, df, snp.col = "SNP", a1.col = "A1") {
  if (is.null(direction.a1)) return(list(type = "none"))

  if (is.character(direction.a1) && length(direction.a1) == 1L && !is.na(direction.a1) && (direction.a1 %in% names(df))) {
    return(list(
      type = "column",
      col = direction.a1,
      vec = toupper(as.character(df[[direction.a1]]))
    ))
  }

  if (is.data.frame(direction.a1) || data.table::is.data.table(direction.a1)) {
    x <- as.data.frame(direction.a1, stringsAsFactors = FALSE)
    if (!nrow(x)) return(list(type = "map", map = setNames(character(), character())))

    if (ncol(x) == 2L) {
      s <- as.character(x[[1]])
      a <- toupper(as.character(x[[2]]))
    } else {
      nms <- names(x)
      snp_hit <- which(nms %in% unique(c(snp.col)))[1]
      a1_hit <- which(nms %in% unique(c(a1.col)))[1]
      if (is.na(snp_hit) || is.na(a1_hit)) {
        stop("direction.a1 data.frame/data.table must have identifiable SNP and A1 columns.", call. = FALSE)
      }
      s <- as.character(x[[snp_hit]])
      a <- toupper(as.character(x[[a1_hit]]))
    }

    ok <- !is.na(s) & nzchar(s) & !is.na(a) & nzchar(a)
    s <- s[ok]; a <- a[ok]
    if (!length(s)) return(list(type = "map", map = setNames(character(), character())))
    mp <- stats::setNames(a, s)
    mp <- mp[!duplicated(names(mp))]
    return(list(type = "map", map = mp))
  }

  if (is.list(direction.a1)) {
    nm <- names(direction.a1)
    if (is.null(nm) || !length(nm)) stop("direction.a1 list must be named with SNP ids.", call. = FALSE)
    vals <- vapply(direction.a1, function(z) {
      if (is.null(z) || !length(z)) return(NA_character_)
      toupper(as.character(z)[1])
    }, character(1))
    ok <- !is.na(nm) & nzchar(nm) & !is.na(vals) & nzchar(vals)
    mp <- stats::setNames(vals[ok], nm[ok])
    mp <- mp[!duplicated(names(mp))]
    return(list(type = "map", map = mp))
  }

  if (is.atomic(direction.a1)) {
    nm <- names(direction.a1)
    if (!is.null(nm) && length(nm)) {
      vals <- toupper(as.character(direction.a1))
      ok <- !is.na(nm) & nzchar(nm) & !is.na(vals) & nzchar(vals)
      mp <- stats::setNames(vals[ok], nm[ok])
      mp <- mp[!duplicated(names(mp))]
      return(list(type = "map", map = mp))
    }
  }

  stop("Unsupported direction.a1 format.", call. = FALSE)
}

