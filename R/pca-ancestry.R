# Reference-PCA-based ancestry estimation: outlier QC on the reference
# (`pca.refqc()`) and ancestry-fraction estimation for target samples
# (`pca.ancestry()`).

.pca_ancestry_as_char1 <- function(x) {
  if (is.null(x) || !length(x)) return(NULL)
  v <- as.character(x)[1]
  if (is.na(v) || !nzchar(v)) NULL else v
}

.pca_ancestry_canon <- function(x) {
  toupper(gsub("[^A-Za-z0-9#]+", "", as.character(x)))
}

.pca_ancestry_pick_col <- function(nm, keys, default = NA_character_) {
  cn <- .pca_ancestry_canon(nm)
  for (k in keys) {
    hit <- which(cn == .pca_ancestry_canon(k))
    if (length(hit)) return(nm[hit[1]])
  }
  default
}

.pca_ancestry_pick_id_col <- function(nm, preferred = NULL) {
  keys <- unique(c(preferred, "IID", "#IID", "ID", "sample.ID", "FID"))
  out <- .pca_ancestry_pick_col(nm, keys, default = NA_character_)
  if (is.na(out) || !nzchar(out)) return(NULL)
  out
}

.pca_ancestry_pc_cols <- function(nm) {
  nm <- as.character(nm)
  hit <- grep("^PC[0-9]+$", nm, ignore.case = FALSE, value = TRUE)
  hit[order(as_int(sub("^PC", "", hit)))]
}

.pca_ancestry_as_pc_dt <- function(x) {
  if (is.null(x)) return(NULL)
  if (data.table::is.data.table(x)) return(data.table::copy(x))
  if (is.data.frame(x)) return(data.table::as.data.table(x))
  if (is.matrix(x)) return(data.table::as.data.table(as.data.frame(x, stringsAsFactors = FALSE)))
  NULL
}

.pca_ancestry_parse_pc_use <- function(x, pc_cols) {
  idx_all <- as_int(sub("^PC", "", pc_cols))
  if (is.null(x) || length(x) == 0L) return(idx_all)
  if (length(x) == 1L) {
    x1 <- tolower(trimws(as.character(x)[1]))
    if (x1 %in% c("all", "auto")) return(idx_all)
  }
  if (is.numeric(x)) {
    y <- as_int(x)
    y <- y[is.finite(y) & !is.na(y) & y > 0]
    return(y)
  }
  xc <- as.character(x)
  xc <- xc[!is.na(xc) & nzchar(xc)]
  if (!length(xc)) return(idx_all)
  y <- suppressWarnings(as.integer(xc))
  need_parse <- is.na(y)
  if (any(need_parse)) {
    z <- toupper(trimws(xc[need_parse]))
    z <- sub("^PC", "", z)
    y2 <- suppressWarnings(as.integer(z))
    y[need_parse] <- y2
  }
  y <- y[is.finite(y) & !is.na(y) & y > 0]
  y
}

.pca_ancestry_read_eval <- function(eval_in) {
  if (is.null(eval_in) || length(eval_in) == 0L) return(NULL)
  if (data.table::is.data.table(eval_in) || is.data.frame(eval_in)) {
    dt <- if (data.table::is.data.table(eval_in)) data.table::copy(eval_in) else data.table::as.data.table(eval_in)
    if (!nrow(dt) && !ncol(dt)) stop("Empty eval table.", call. = FALSE)
    cn <- names(dt)
    ev_col <- .pca_ancestry_pick_col(cn, c("EIGENVAL", "EVAL", "VALUE"))
    if (is.na(ev_col)) {
      num_cols <- cn[vapply(dt, function(v) {
        x <- suppressWarnings(as.numeric(v))
        any(is.finite(x) & !is.na(x))
      }, logical(1))]
      ev_col <- if (length(num_cols)) num_cols[1] else cn[1]
    }
    ev <- suppressWarnings(as.numeric(dt[[ev_col]]))
    ev <- ev[is.finite(ev) & !is.na(ev)]
    if (!length(ev)) stop("No numeric eigenvalue found in eval table.", call. = FALSE)
    return(ev)
  }
  if (is.character(eval_in) && length(eval_in) == 1L) {
    pth <- abs_path(eval_in)[1]
    if (!is.na(pth) && nzchar(pth) && file.exists(pth)) {
      dt <- data.table::fread(pth, data.table = TRUE, check.names = FALSE, showProgress = FALSE)
      if (!nrow(dt) && !ncol(dt)) stop("Empty eval file: ", pth, call. = FALSE)
      cn <- names(dt)
      ev_col <- .pca_ancestry_pick_col(cn, c("EIGENVAL", "EVAL", "VALUE"))
      if (is.na(ev_col)) {
        num_cols <- cn[vapply(dt, function(v) {
          x <- suppressWarnings(as.numeric(v))
          any(is.finite(x) & !is.na(x))
        }, logical(1))]
        ev_col <- if (length(num_cols)) num_cols[1] else cn[1]
      }
      ev <- suppressWarnings(as.numeric(dt[[ev_col]]))
      ev <- ev[is.finite(ev) & !is.na(ev)]
      if (!length(ev)) stop("No numeric eigenvalue found in eval file: ", pth, call. = FALSE)
      return(ev)
    }
  }
  ev <- suppressWarnings(as.numeric(eval_in))
  if (!length(ev)) return(NULL)
  names(ev) <- names(eval_in)
  ev
}

.pca_ancestry_make_eval_named <- function(eval_in, pc_cols) {
  ev <- .pca_ancestry_read_eval(eval_in)
  if (is.null(ev) || !length(ev)) return(NULL)
  nm <- names(ev)
  if (is.null(nm) || !length(nm) || all(is.na(nm) | !nzchar(nm))) {
    nm <- pc_cols[seq_len(min(length(pc_cols), length(ev)))]
    if (length(nm) < length(ev)) {
      nm <- c(nm, paste0("PC", seq_len(length(ev) - length(nm)) + length(nm)))
    }
  }
  nm <- as.character(nm)
  bad_nm <- is.na(nm) | !nzchar(nm)
  if (any(bad_nm)) nm[bad_nm] <- paste0("PC", seq_len(sum(bad_nm)))
  out <- as.numeric(ev)
  names(out) <- nm
  out
}

.pca_ancestry_is_missingish <- function(x) {
  if (is.character(x)) {
    return(is.na(x) | !nzchar(trimws(x)))
  }
  is.na(x)
}

.pca_ancestry_as_ref_flag <- function(x) {
  if (is.null(x)) return(logical())
  if (is.logical(x)) return(as.logical(x))
  xc <- tolower(trimws(as.character(x)))
  out <- rep(NA, length(xc))
  out[xc %in% c("true", "t", "1", "yes", "y", "ref", "reference")] <- TRUE
  out[xc %in% c("false", "f", "0", "no", "n", "target")] <- FALSE
  raw_lgl <- suppressWarnings(as.logical(x))
  fill_idx <- is.na(out) & !is.na(raw_lgl)
  if (any(fill_idx)) out[fill_idx] <- raw_lgl[fill_idx]
  as.logical(out)
}

.pca_ancestry_resolve_input <- function(pc, eval = NULL, pc.use = 1:10, id.col = "IID") {
  require_pkg("data.table")
  eval_use <- eval
  dt <- NULL
  pc_in <- pc
  if (is.list(pc_in) && !data.table::is.data.table(pc_in) && !is.data.frame(pc_in)) {
    cand_keys <- c("pc", "PC", "data", "scores", "eigenvec")
    for (kk in cand_keys) {
      if (is.null(pc_in[[kk]])) next
      dt0 <- .pca_ancestry_as_pc_dt(pc_in[[kk]])
      if (is.null(dt0)) next
      dt <- dt0
      break
    }
    if (is.null(dt) && length(pc_in)) {
      for (ii in seq_along(pc_in)) {
        dt0 <- .pca_ancestry_as_pc_dt(pc_in[[ii]])
        if (is.null(dt0)) next
        pcc0 <- .pca_ancestry_pc_cols(names(dt0))
        if (length(pcc0) >= 1L) {
          dt <- dt0
          break
        }
      }
    }
    if (is.null(dt)) stop("pc must be pca() return object or data.frame/data.table with PC columns.", call. = FALSE)
    if (is.null(eval_use) || length(eval_use) == 0L) {
      if (!is.null(pc_in$eigenvalue) && length(pc_in$eigenvalue)) eval_use <- pc_in$eigenvalue
      else if (!is.null(pc_in$prop.var) && length(pc_in$prop.var)) eval_use <- pc_in$prop.var
      else if (!is.null(pc_in$eval) && length(pc_in$eval)) eval_use <- pc_in$eval
      else if (!is.null(attr(pc_in, "eigenvalue")) && length(attr(pc_in, "eigenvalue"))) eval_use <- attr(pc_in, "eigenvalue")
      else if (!is.null(attr(pc_in, "prop.var")) && length(attr(pc_in, "prop.var"))) eval_use <- attr(pc_in, "prop.var")
    }
  } else {
    dt <- .pca_ancestry_as_pc_dt(pc_in)
  }
  if (is.null(dt) || !nrow(dt)) stop("No rows in pc.", call. = FALSE)

  id_col_use <- .pca_ancestry_as_char1(id.col)
  if (is.null(id_col_use) || !(id_col_use %in% names(dt))) {
    id_candidates <- intersect(c("IID", "#IID", "ID", "sample.ID", "FID"), names(dt))
    if (length(id_candidates)) id_col_use <- id_candidates[1]
  }
  if (is.null(id_col_use) || !(id_col_use %in% names(dt))) stop("id.col not found in pc.", call. = FALSE)

  dt[, .ROW_ID := .I]
  dt[, IID_USE := as.character(get(id_col_use))]
  dt[is.na(IID_USE) | !nzchar(IID_USE), IID_USE := paste0("row", .I)]

  pcc <- .pca_ancestry_pc_cols(names(dt))
  if (!length(pcc)) stop("At least one PC column (PC1, PC2, ...) is required.", call. = FALSE)
  for (nm in pcc) dt[, (nm) := suppressWarnings(as.numeric(get(nm)))]

  pc_idx <- .pca_ancestry_parse_pc_use(pc.use, pcc)
  pc_idx <- unique(pc_idx)
  if (!length(pc_idx)) stop("pc.use must contain at least one PC.", call. = FALSE)
  pc_names <- paste0("PC", pc_idx)
  miss_pc <- pc_names[!(pc_names %in% pcc)]
  if (length(miss_pc)) stop("Missing PC columns in data: ", paste(miss_pc, collapse = ", "), call. = FALSE)

  eval_named <- .pca_ancestry_make_eval_named(eval_use, pcc)
  if (is.null(eval_named) || !length(eval_named)) {
    stop("eval is required. Provide eval or use a pca()/pca.projection() return object with eigenvalue.", call. = FALSE)
  }
  if (any(!(pc_names %in% names(eval_named)))) {
    stop("Selected PCs are not available in eval: ", paste(pc_names[!(pc_names %in% names(eval_named))], collapse = ", "), call. = FALSE)
  }

  x <- as.matrix(dt[, ..pc_names])
  storage.mode(x) <- "double"
  if (any(!is.finite(x) | is.na(x))) stop("Selected PC columns contain NA/non-finite values.", call. = FALSE)

  list(
    dt = dt,
    x = x,
    pc_names = pc_names,
    pc_names_all = pcc,
    eval_named = eval_named,
    id_col_use = id_col_use
  )
}

.pca_ancestry_attach_label_data <- function(dt, id_col_use, label.data = NULL, label.col = NULL) {
  dt <- data.table::copy(dt)
  label_col_use <- .pca_ancestry_as_char1(label.col)
  ref_pc_col <- .pca_ancestry_pick_col(names(dt), c("reference"))
  if (is.na(ref_pc_col)) {
    dt[, .REFERENCE_PC := NA]
  } else {
    dt[, .REFERENCE_PC := get(ref_pc_col)]
  }

  if (!is.null(label.data)) {
    if (is.atomic(label.data) && !is.list(label.data) && !is.data.frame(label.data) && !data.table::is.data.table(label.data)) {
      if (length(label.data) != nrow(dt)) stop("label.data vector length must match nrow(pc).", call. = FALSE)
      dt[, .LABEL_INPUT := as.character(label.data)]
    } else {
      gd <- if (data.table::is.data.table(label.data)) data.table::copy(label.data) else if (is.data.frame(label.data)) data.table::as.data.table(label.data) else NULL
      if (is.null(gd)) stop("label.data must be vector/data.frame/data.table.", call. = FALSE)
      ref_gd_col <- .pca_ancestry_pick_col(names(gd), c("reference"))
      id_col_gd <- NULL
      label_from_gd <- NULL

      if (ncol(gd) == 2L && is.null(label_col_use)) {
        id_col_gd <- names(gd)[1]
        label_from_gd <- names(gd)[2]
      } else {
        id_col_gd <- .pca_ancestry_pick_id_col(names(gd), preferred = id_col_use)
        if (is.null(id_col_gd) || !(id_col_gd %in% names(gd))) stop("id.col not found in label.data.", call. = FALSE)
        if (!is.null(label_col_use)) {
          if (label_col_use %in% names(gd)) {
            label_from_gd <- label_col_use
          } else if (!(label_col_use %in% names(dt))) {
            stop("label.col not found in label.data or pc: ", label_col_use, call. = FALSE)
          }
        } else {
          cands <- setdiff(names(gd), c(id_col_gd, ref_gd_col))
          if (length(cands)) label_from_gd <- cands[1]
        }
      }

      keep_cols <- unique(c(id_col_gd, label_from_gd, ref_gd_col))
      keep_cols <- keep_cols[!is.na(keep_cols) & nzchar(keep_cols)]
      gd <- gd[, ..keep_cols]
      data.table::setnames(gd, id_col_gd, "IDM")
      if (!is.null(label_from_gd) && (label_from_gd %in% names(gd))) data.table::setnames(gd, label_from_gd, ".LABEL_INPUT")
      if (!is.na(ref_gd_col) && (ref_gd_col %in% names(gd))) data.table::setnames(gd, ref_gd_col, ".REFERENCE_INPUT")
      gd[, IDM := as.character(IDM)]
      gd <- gd[!is.na(IDM) & nzchar(IDM)]
      data.table::setkey(gd, IDM)
      dt <- merge(dt, gd, by.x = "IID_USE", by.y = "IDM", all.x = TRUE, sort = FALSE)

      if (!(".LABEL_INPUT" %in% names(dt)) && !is.null(label_col_use) && (label_col_use %in% names(dt))) {
        dt[, .LABEL_INPUT := as.character(get(label_col_use))]
      }
    }
  } else if (!is.null(label_col_use)) {
    if (!(label_col_use %in% names(dt))) stop("label.col not found in pc: ", label_col_use, call. = FALSE)
    dt[, .LABEL_INPUT := as.character(get(label_col_use))]
  }

  if (!(".LABEL_INPUT" %in% names(dt))) dt[, .LABEL_INPUT := NA_character_]
  if (!(".REFERENCE_INPUT" %in% names(dt))) dt[, .REFERENCE_INPUT := .REFERENCE_PC]

  ref_missing <- .pca_ancestry_is_missingish(dt$.REFERENCE_INPUT)
  if (any(ref_missing)) dt[ref_missing, .REFERENCE_INPUT := .REFERENCE_PC[ref_missing]]

  dt[, .LABEL_INPUT := as.character(.LABEL_INPUT)]
  dt[is.na(.LABEL_INPUT) | !nzchar(trimws(.LABEL_INPUT)), .LABEL_INPUT := NA_character_]
  dt[, .REFERENCE_FLAG := .pca_ancestry_as_ref_flag(.REFERENCE_INPUT)]
  dt[, .REFERENCE_PC := NULL]
  dt
}

.pca_ancestry_attach_subgroup_data <- function(dt, id_col_use, subgroup.data = NULL, subgroup.col = NULL) {
  dt <- data.table::copy(dt)
  subgroup_col_use <- .pca_ancestry_as_char1(subgroup.col)

  if (!is.null(subgroup.data)) {
    if (is.atomic(subgroup.data) && !is.list(subgroup.data) && !is.data.frame(subgroup.data) && !data.table::is.data.table(subgroup.data)) {
      if (length(subgroup.data) != nrow(dt)) stop("subgroup.data vector length must match nrow(pc).", call. = FALSE)
      dt[, .SUBGROUP_INPUT := as.character(subgroup.data)]
    } else {
      gd <- if (data.table::is.data.table(subgroup.data)) data.table::copy(subgroup.data) else if (is.data.frame(subgroup.data)) data.table::as.data.table(subgroup.data) else NULL
      if (is.null(gd)) stop("subgroup.data must be vector/data.frame/data.table.", call. = FALSE)
      id_col_gd <- NULL
      subgroup_from_gd <- NULL

      if (ncol(gd) == 2L && is.null(subgroup_col_use)) {
        id_col_gd <- names(gd)[1]
        subgroup_from_gd <- names(gd)[2]
      } else {
        id_col_gd <- .pca_ancestry_pick_id_col(names(gd), preferred = id_col_use)
        if (is.null(id_col_gd) || !(id_col_gd %in% names(gd))) stop("id.col not found in subgroup.data.", call. = FALSE)
        if (!is.null(subgroup_col_use)) {
          if (subgroup_col_use %in% names(gd)) {
            subgroup_from_gd <- subgroup_col_use
          } else if (!(subgroup_col_use %in% names(dt))) {
            stop("subgroup.col not found in subgroup.data or pc: ", subgroup_col_use, call. = FALSE)
          }
        } else {
          cands <- setdiff(names(gd), id_col_gd)
          if (length(cands)) subgroup_from_gd <- cands[1]
        }
      }

      keep_cols <- unique(c(id_col_gd, subgroup_from_gd))
      keep_cols <- keep_cols[!is.na(keep_cols) & nzchar(keep_cols)]
      gd <- gd[, ..keep_cols]
      data.table::setnames(gd, id_col_gd, "IDM")
      if (!is.null(subgroup_from_gd) && (subgroup_from_gd %in% names(gd))) data.table::setnames(gd, subgroup_from_gd, ".SUBGROUP_INPUT")
      gd[, IDM := as.character(IDM)]
      gd <- gd[!is.na(IDM) & nzchar(IDM)]
      data.table::setkey(gd, IDM)
      dt <- merge(dt, gd, by.x = "IID_USE", by.y = "IDM", all.x = TRUE, sort = FALSE)

      if (!(".SUBGROUP_INPUT" %in% names(dt)) && !is.null(subgroup_col_use) && (subgroup_col_use %in% names(dt))) {
        dt[, .SUBGROUP_INPUT := as.character(get(subgroup_col_use))]
      }
    }
  } else if (!is.null(subgroup_col_use)) {
    if (!(subgroup_col_use %in% names(dt))) stop("subgroup.col not found in pc: ", subgroup_col_use, call. = FALSE)
    dt[, .SUBGROUP_INPUT := as.character(get(subgroup_col_use))]
  }

  if (!(".SUBGROUP_INPUT" %in% names(dt))) dt[, .SUBGROUP_INPUT := NA_character_]
  dt[, .SUBGROUP_INPUT := as.character(.SUBGROUP_INPUT)]
  dt[is.na(.SUBGROUP_INPUT) | !nzchar(trimws(.SUBGROUP_INPUT)), .SUBGROUP_INPUT := NA_character_]
  dt
}

.pca_ancestry_weight_vector <- function(eval_named, pc_names, power = 1) {
  power <- as_num(power)[1]
  if (!is.finite(power) || is.na(power) || power < 0) stop("power must be a non-negative finite value.", call. = FALSE)
  w <- as.numeric(eval_named[pc_names])
  if (any(!is.finite(w) | is.na(w))) stop("Selected PC eigenvalues contain NA/non-finite values.", call. = FALSE)
  w <- w ^ power
  sw <- sum(w, na.rm = TRUE)
  if (!is.finite(sw) || is.na(sw) || sw <= 0) stop("Selected PC eigenvalues must sum to a positive value.", call. = FALSE)
  w / sw
}

.pca_ancestry_outlier_pc_spec <- function(pc_names, eval_named, min_pc_reduce = 5L, prop_target = 0.9, weight.power = 1, pc_names_all = NULL) {
  pc_names <- as.character(pc_names)
  if (is.null(pc_names_all) || !length(pc_names_all)) pc_names_all <- pc_names
  pc_names_all <- as.character(pc_names_all)
  pc_names_all <- pc_names_all[pc_names_all %in% names(eval_named)]
  eval_all <- as_num(eval_named[pc_names_all])
  if (any(!is.finite(eval_all) | is.na(eval_all))) stop("Outlier PC eigenvalues contain NA/non-finite values.", call. = FALSE)
  prop_all <- eval_all / sum(eval_all)
  keep_n <- length(pc_names_all)
  min_pc_reduce <- as_int(min_pc_reduce)[1]
  prop_target <- as_num(prop_target)[1]
  if (is.finite(min_pc_reduce) && !is.na(min_pc_reduce) && keep_n >= min_pc_reduce &&
      is.finite(prop_target) && !is.na(prop_target) && prop_target > 0 && prop_target < 1) {
    cum_prop <- cumsum(prop_all)
    hit <- which(cum_prop >= prop_target)
    if (length(hit)) keep_n <- max(1L, hit[1])
  }
  pc_keep_global <- pc_names_all[seq_len(keep_n)]
  pc_keep <- pc_names[pc_names %in% pc_keep_global]
  if (!length(pc_keep)) pc_keep <- pc_names[1]
  list(
    pc_names = pc_keep,
    pc_idx = match(pc_keep, pc_names),
    weight = .pca_ancestry_weight_vector(eval_named = eval_named, pc_names = pc_keep, power = weight.power),
    cum_prop = sum(prop_all[seq_len(keep_n)]),
    label = if (length(pc_keep) == 1L) pc_keep[1] else sprintf("%s-%s", pc_keep[1], pc_keep[length(pc_keep)])
  )
}

.pca_ancestry_posterior_pc_spec <- function(pc_names, eval_named, min_pc_reduce = 5L, prop_target = 0.95, weight.power = 1, pc_names_all = NULL) {
  .pca_ancestry_outlier_pc_spec(
    pc_names = pc_names,
    eval_named = eval_named,
    min_pc_reduce = min_pc_reduce,
    prop_target = prop_target,
    weight.power = weight.power,
    pc_names_all = pc_names_all
  )
}

.pca_ancestry_fit_scale <- function(x_ref) {
  center <- colMeans(x_ref)
  scale <- apply(x_ref, 2L, stats::sd)
  scale[!is.finite(scale) | is.na(scale) | scale <= 0] <- 1
  list(center = center, scale = scale)
}

.pca_ancestry_transform <- function(x, center, scale, weight) {
  z <- sweep(x, 2L, center, "-", check.margin = FALSE)
  z <- sweep(z, 2L, scale, "/", check.margin = FALSE)
  sweep(z, 2L, sqrt(weight), "*", check.margin = FALSE)
}

.pca_ancestry_topk_self <- function(xw_ref, kmax) {
  n <- nrow(xw_ref)
  if (!is.finite(kmax) || is.na(kmax) || kmax < 1L) stop("kmax must be >= 1.", call. = FALSE)
  kmax <- min(as_int(kmax), max(1L, n - 1L))
  rs <- rowSums(xw_ref^2)
  dmat <- outer(rs, rs, "+") - 2 * tcrossprod(xw_ref)
  dmat[!is.finite(dmat)] <- Inf
  dmat[dmat < 0] <- 0
  dmat[is.finite(dmat)] <- sqrt(dmat[is.finite(dmat)])
  diag(dmat) <- Inf
  idx <- matrix(NA_integer_, n, kmax)
  dist <- matrix(NA_real_, n, kmax)
  for (i in seq_len(n)) {
    ord <- order(dmat[i, ], decreasing = FALSE, na.last = NA)
    take <- ord[seq_len(min(kmax, length(ord)))]
    idx[i, seq_along(take)] <- take
    dist[i, seq_along(take)] <- dmat[i, take]
  }
  list(idx = idx, dist = dist)
}

.pca_ancestry_topk_query <- function(xw_query, xw_ref, kmax, chunk_size = 1024L) {
  nq <- nrow(xw_query)
  nr <- nrow(xw_ref)
  if (!is.finite(kmax) || is.na(kmax) || kmax < 1L) stop("kmax must be >= 1.", call. = FALSE)
  kmax <- min(as_int(kmax), nr)
  if (nq == 0L) {
    return(list(
      idx = matrix(NA_integer_, 0L, kmax),
      dist = matrix(NA_real_, 0L, kmax)
    ))
  }
  ref_norm <- rowSums(xw_ref^2)
  qry_norm <- rowSums(xw_query^2)
  idx <- matrix(NA_integer_, nq, kmax)
  dist <- matrix(NA_real_, nq, kmax)
  chunk_size <- as_int(chunk_size)[1]
  if (!is.finite(chunk_size) || is.na(chunk_size) || chunk_size < 1L) chunk_size <- 1024L
  starts <- seq.int(1L, nq, by = chunk_size)
  for (st in starts) {
    en <- min(nq, st + chunk_size - 1L)
    dmat <- outer(qry_norm[st:en], ref_norm, "+") - 2 * tcrossprod(xw_query[st:en, , drop = FALSE], xw_ref)
    dmat[!is.finite(dmat)] <- Inf
    dmat[dmat < 0] <- 0
    dmat[is.finite(dmat)] <- sqrt(dmat[is.finite(dmat)])
    for (ii in seq_len(nrow(dmat))) {
      ord <- order(dmat[ii, ], decreasing = FALSE, na.last = NA)
      take <- ord[seq_len(min(kmax, length(ord)))]
      idx[st + ii - 1L, seq_along(take)] <- take
      dist[st + ii - 1L, seq_along(take)] <- dmat[ii, take]
    }
  }
  list(idx = idx, dist = dist)
}

.pca_ancestry_group_topk_self <- function(xw_ref, labels_ref, groups, kmax) {
  labels_ref <- as.character(labels_ref)
  groups <- as.character(groups)
  n <- nrow(xw_ref)
  kmax <- as_int(kmax)[1]
  if (!is.finite(kmax) || is.na(kmax) || kmax < 1L) stop("kmax must be >= 1.", call. = FALSE)
  rs <- rowSums(xw_ref^2)
  dmat <- outer(rs, rs, "+") - 2 * tcrossprod(xw_ref)
  dmat[!is.finite(dmat)] <- Inf
  dmat[dmat < 0] <- 0
  diag(dmat) <- Inf

  out <- vector("list", length(groups))
  names(out) <- groups
  for (g in groups) {
    idx_g <- which(labels_ref == g)
    kg <- min(kmax, length(idx_g))
    dist_g <- matrix(NA_real_, nrow = n, ncol = max(1L, kg))
    if (!length(idx_g)) {
      out[[g]] <- dist_g[, 0L, drop = FALSE]
      next
    }
    for (i in seq_len(n)) {
      d0 <- dmat[i, idx_g]
      if (labels_ref[i] == g) {
        hit <- which(idx_g == i)
        if (length(hit)) d0[hit[1]] <- Inf
      }
      ord <- order(d0, decreasing = FALSE, na.last = NA)
      ord <- ord[is.finite(d0[ord])]
      if (!length(ord)) next
      take_n <- min(kg, length(ord))
      dist_g[i, seq_len(take_n)] <- d0[ord[seq_len(take_n)]]
    }
    out[[g]] <- dist_g
  }
  out
}

.pca_ancestry_group_topk_query <- function(xw_query, xw_ref, labels_ref, groups, kmax, chunk_size = 1024L) {
  labels_ref <- as.character(labels_ref)
  groups <- as.character(groups)
  out <- vector("list", length(groups))
  names(out) <- groups
  for (g in groups) {
    idx_g <- which(labels_ref == g)
    if (!length(idx_g)) {
      out[[g]] <- matrix(NA_real_, nrow(xw_query), 0L)
      next
    }
    kg <- min(as_int(kmax)[1], length(idx_g))
    topk <- .pca_ancestry_topk_query(
      xw_query = xw_query,
      xw_ref = xw_ref[idx_g, , drop = FALSE],
      kmax = kg,
      chunk_size = chunk_size
    )
    out[[g]] <- topk$dist
  }
  out
}

.pca_ancestry_complement_topk_self <- function(xw_ref, labels_ref, groups, kmax) {
  labels_ref <- as.character(labels_ref)
  groups <- as.character(groups)
  n <- nrow(xw_ref)
  kmax <- as_int(kmax)[1]
  if (!is.finite(kmax) || is.na(kmax) || kmax < 1L) stop("kmax must be >= 1.", call. = FALSE)
  rs <- rowSums(xw_ref^2)
  dmat <- outer(rs, rs, "+") - 2 * tcrossprod(xw_ref)
  dmat[!is.finite(dmat)] <- Inf
  dmat[dmat < 0] <- 0
  dmat[is.finite(dmat)] <- sqrt(dmat[is.finite(dmat)])
  diag(dmat) <- Inf

  out <- vector("list", length(groups))
  names(out) <- groups
  for (g in groups) {
    idx_g <- which(labels_ref != g)
    kg <- min(kmax, length(idx_g))
    dist_g <- matrix(NA_real_, nrow = n, ncol = max(1L, kg))
    if (!length(idx_g)) {
      out[[g]] <- dist_g[, 0L, drop = FALSE]
      next
    }
    for (i in seq_len(n)) {
      d0 <- dmat[i, idx_g]
      ord <- order(d0, decreasing = FALSE, na.last = NA)
      ord <- ord[is.finite(d0[ord])]
      if (!length(ord)) next
      take_n <- min(kg, length(ord))
      dist_g[i, seq_len(take_n)] <- d0[ord[seq_len(take_n)]]
    }
    out[[g]] <- dist_g
  }
  out
}

.pca_ancestry_complement_topk_query <- function(xw_query, xw_ref, labels_ref, groups, kmax, chunk_size = 1024L) {
  labels_ref <- as.character(labels_ref)
  groups <- as.character(groups)
  out <- vector("list", length(groups))
  names(out) <- groups
  for (g in groups) {
    idx_g <- which(labels_ref != g)
    if (!length(idx_g)) {
      out[[g]] <- matrix(NA_real_, nrow(xw_query), 0L)
      next
    }
    kg <- min(as_int(kmax)[1], length(idx_g))
    topk <- .pca_ancestry_topk_query(
      xw_query = xw_query,
      xw_ref = xw_ref[idx_g, , drop = FALSE],
      kmax = kg,
      chunk_size = chunk_size
    )
    out[[g]] <- topk$dist
  }
  out
}

.pca_ancestry_prior_vec <- function(groups, labels_ref, prior = "empirical", prior.user = NULL) {
  groups <- as.character(groups)
  labels_ref <- as.character(labels_ref)
  prior_mode <- .pca_ancestry_as_char1(prior)
  if (is.null(prior_mode)) prior_mode <- "empirical"
  prior_mode <- tolower(trimws(prior_mode))

  if (prior_mode == "empirical") {
    tab <- table(factor(labels_ref, levels = groups))
    out <- as_num(tab) / sum(tab)
    names(out) <- groups
    return(out)
  }

  if (prior_mode == "uniform") {
    out <- rep(1 / length(groups), length(groups))
    names(out) <- groups
    return(out)
  }

  if (prior_mode == "user") {
    if (is.null(prior.user) || !length(prior.user)) stop("prior.user must be provided when prior='user'.", call. = FALSE)
    pu <- as_num(prior.user)
    if (all(is.na(names(prior.user)) | !nzchar(names(prior.user)))) {
      if (length(pu) != length(groups)) stop("prior.user must have length equal to number of groups or be a named vector.", call. = FALSE)
      names(pu) <- groups
    }
    pu <- pu[groups]
    if (any(!is.finite(pu) | is.na(pu) | pu < 0)) stop("prior.user must contain non-negative finite values for all groups.", call. = FALSE)
    su <- sum(pu)
    if (!is.finite(su) || is.na(su) || su <= 0) stop("prior.user must sum to a positive value.", call. = FALSE)
    out <- pu / su
    names(out) <- groups
    return(out)
  }

  stop("prior must be one of: empirical, uniform, user", call. = FALSE)
}

.pca_ancestry_kth_distance <- function(dist_mat, k) {
  if (is.null(dist_mat) || !length(dist_mat)) return(rep(NA_real_, 0L))
  k <- as_int(k)[1]
  if (!is.finite(k) || is.na(k) || k < 1L) stop("k must be >= 1.", call. = FALSE)
  n <- nrow(dist_mat)
  out <- rep(NA_real_, n)
  if (!n) return(out)
  for (i in seq_len(n)) {
    d0 <- as_num(dist_mat[i, , drop = TRUE])
    d0 <- d0[is.finite(d0) & !is.na(d0)]
    if (!length(d0)) next
    take_n <- min(k, length(d0))
    out[i] <- d0[take_n]
  }
  out
}

.pca_ancestry_indep_posterior <- function(group_dist_pos,
                                          group_dist_neg,
                                          groups,
                                          labels_ref,
                                          k,
                                          dim_p,
                                          prior = "empirical",
                                          prior.user = NULL) {
  groups <- as.character(groups)
  labels_ref <- as.character(labels_ref)
  k <- as_int(k)[1]
  dim_p <- as_int(dim_p)[1]
  if (!is.finite(k) || is.na(k) || k < 1L) stop("k must be >= 1.", call. = FALSE)
  if (!is.finite(dim_p) || is.na(dim_p) || dim_p < 1L) stop("dim_p must be >= 1.", call. = FALSE)
  n <- if (length(group_dist_pos)) nrow(group_dist_pos[[1]]) else 0L
  g <- length(groups)
  pp_mat <- matrix(NA_real_, n, g)
  colnames(pp_mat) <- groups
  label_max <- rep(NA_character_, n)
  prior_vec <- .pca_ancestry_prior_vec(groups = groups, labels_ref = labels_ref, prior = prior, prior.user = prior.user)
  n_total <- length(labels_ref)
  tiny <- .Machine$double.eps

  d_pos <- lapply(group_dist_pos, .pca_ancestry_kth_distance, k = k)
  d_neg <- lapply(group_dist_neg, .pca_ancestry_kth_distance, k = k)

  for (gg in seq_along(groups)) {
    gname <- groups[gg]
    n_pos <- sum(labels_ref == gname, na.rm = TRUE)
    n_neg <- n_total - n_pos
    k_pos <- min(k, n_pos)
    k_neg <- min(k, n_neg)
    pi_g <- prior_vec[gname]
    pp_g <- rep(NA_real_, n)

    if (n_pos >= 1L && n_neg >= 1L && is.finite(pi_g) && !is.na(pi_g) && pi_g > 0 && pi_g < 1) {
      r_pos <- pmax(as_num(d_pos[[gg]]), tiny)
      r_neg <- pmax(as_num(d_neg[[gg]]), tiny)
      dens_pos <- k_pos / (n_pos * (r_pos ^ dim_p))
      dens_neg <- k_neg / (n_neg * (r_neg ^ dim_p))
      num <- pi_g * dens_pos
      den <- num + (1 - pi_g) * dens_neg
      ok <- is.finite(num) & !is.na(num) & is.finite(den) & !is.na(den) & den > 0
      pp_g[ok] <- num[ok] / den[ok]
    } else if (n_pos >= 1L && n_neg < 1L) {
      pp_g[] <- 1
    }

    pp_mat[, gg] <- pp_g
  }

  for (i in seq_len(n)) {
    x <- as_num(pp_mat[i, ])
    ok <- is.finite(x) & !is.na(x)
    if (!any(ok)) next
    im <- which.max(replace(x, !ok, -Inf))
    label_max[i] <- groups[im[1]]
  }

  list(
    score = pp_mat,
    label.max = label_max
  )
}

.pca_ancestry_stats_from_summary <- function(sum_vec, sumsq_vec, n, var_default = NULL) {
  sum_vec <- as_num(sum_vec)
  sumsq_vec <- as_num(sumsq_vec)
  n <- as_int(n)[1]
  p <- length(sum_vec)
  if (!p) return(list(center = numeric(), var = numeric()))

  center <- rep(NA_real_, p)
  var <- rep(NA_real_, p)
  if (is.finite(n) && !is.na(n) && n >= 1L) {
    center <- sum_vec / n
  }
  if (is.finite(n) && !is.na(n) && n >= 2L) {
    var <- (sumsq_vec - (sum_vec^2) / n) / (n - 1)
  }
  ok <- is.finite(var) & !is.na(var) & (var > 0)
  if (!is.null(var_default)) {
    fb <- as_num(var_default)
    if (length(fb) == 1L && p > 1L) fb <- rep(fb, p)
    fb_ok <- is.finite(fb) & !is.na(fb) & (fb > 0)
    fill <- !ok & fb_ok
    if (any(fill)) var[fill] <- fb[fill]
  }
  var[!is.finite(var) | is.na(var) | var <= 0] <- 1
  center[!is.finite(center) | is.na(center)] <- 0
  list(center = center, var = var)
}

.pca_ancestry_diag_distance_matrix <- function(x, center, var, weight = NULL) {
  if (!nrow(x)) return(numeric())
  d <- sweep(x, 2L, center, "-", check.margin = FALSE)
  d <- d^2
  if (!is.null(var)) {
    v <- as_num(var)
    if (length(v) == 1L && ncol(d) > 1L) v <- rep(v, ncol(d))
    if (length(v) != ncol(d)) stop("var length must match number of PCs.", call. = FALSE)
    v[!is.finite(v) | is.na(v) | v <= 0] <- 1
    d <- sweep(d, 2L, v, "/", check.margin = FALSE)
  }
  if (!is.null(weight)) {
    w <- as_num(weight)
    if (length(w) == 1L && ncol(d) > 1L) w <- rep(w, ncol(d))
    if (length(w) != ncol(d)) stop("weight length must match number of PCs.", call. = FALSE)
    d <- sweep(d, 2L, w, "*", check.margin = FALSE)
  }
  out <- rowSums(d)
  out[!is.finite(out) | is.na(out)] <- NA_real_
  ok <- is.finite(out) & !is.na(out) & out >= 0
  out[ok] <- sqrt(out[ok])
  out
}

.pca_ancestry_sample_int <- function(n, size, seed = NULL) {
  n <- as_int(n)[1]
  size <- as_int(size)[1]
  if (!is.finite(n) || is.na(n) || n < 1L || !is.finite(size) || is.na(size) || size < 1L) return(integer())
  size <- min(size, n)
  old_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (isTRUE(old_exists)) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (isTRUE(old_exists)) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  if (!is.null(seed) && length(seed)) set.seed(as_int(seed)[1])
  sample.int(n, size = size, replace = FALSE)
}

.pca_ancestry_lof_unit_model <- function(x_unit, weight, k, tau = 3) {
  n <- nrow(x_unit)
  p <- ncol(x_unit)
  if (n < 2L || p < 1L) return(NULL)
  tau <- as_num(tau)[1]
  if (!is.finite(tau) || is.na(tau) || tau <= 0) tau <- 3
  weight <- as_num(weight)
  if (length(weight) == 1L && p > 1L) weight <- rep(weight, p)
  if (length(weight) != p) stop("outlier weight length must match number of PCs.", call. = FALSE)
  weight[!is.finite(weight) | is.na(weight) | weight < 0] <- 0
  xw_ref <- sweep(x_unit, 2L, sqrt(weight), "*", check.margin = FALSE)
  k_use <- min(as_int(k)[1], n - 1L)
  if (!is.finite(k_use) || is.na(k_use) || k_use < 1L) return(NULL)
  topk <- .pca_ancestry_topk_self(xw_ref, kmax = k_use)
  kdist_ref <- topk$dist[, k_use]
  reach <- topk$dist
  for (j in seq_len(k_use)) {
    idx_j <- topk$idx[, j]
    reach[, j] <- pmax(reach[, j], kdist_ref[idx_j])
  }
  lrd_ref <- 1 / rowMeans(reach)
  lrd_ref[!is.finite(lrd_ref) | is.na(lrd_ref) | lrd_ref <= 0] <- NA_real_
  lof_ref <- rep(NA_real_, n)
  pdist_ref <- tau * sqrt(rowMeans(topk$dist^2))
  pdist_ref[!is.finite(pdist_ref) | is.na(pdist_ref) | pdist_ref <= 0] <- NA_real_
  for (i in seq_len(n)) {
    idx_i <- as_int(topk$idx[i, , drop = TRUE])
    ok <- !is.na(idx_i)
    if (!any(ok) || !is.finite(lrd_ref[i]) || is.na(lrd_ref[i]) || lrd_ref[i] <= 0) next
    idx_i <- idx_i[ok]
    lrd_nb <- lrd_ref[idx_i]
    ok_nb <- is.finite(lrd_nb) & !is.na(lrd_nb) & lrd_nb > 0
    if (!any(ok_nb)) next
    lof_ref[i] <- mean(lrd_nb[ok_nb] / lrd_ref[i])
  }
  plof_ref <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    idx_i <- as_int(topk$idx[i, , drop = TRUE])
    ok <- !is.na(idx_i)
    if (!any(ok) || !is.finite(pdist_ref[i]) || is.na(pdist_ref[i]) || pdist_ref[i] <= 0) next
    idx_i <- idx_i[ok]
    pdist_nb <- pdist_ref[idx_i]
    ok_nb <- is.finite(pdist_nb) & !is.na(pdist_nb) & pdist_nb > 0
    if (!any(ok_nb)) next
    plof_ref[i] <- pdist_ref[i] / mean(pdist_nb[ok_nb]) - 1
  }
  plof_finite <- plof_ref[is.finite(plof_ref) & !is.na(plof_ref)]
  nplof <- if (length(plof_finite)) tau * sqrt(mean(plof_finite^2)) else NA_real_
  if (!is.finite(nplof) || is.na(nplof) || nplof <= 0) nplof <- NA_real_
  loop_ref <- rep(NA_real_, n)
  ok_loop <- is.finite(plof_ref) & !is.na(plof_ref) & is.finite(nplof) & !is.na(nplof) & nplof > 0
  if (any(ok_loop)) {
    loop_ref[ok_loop] <- pmax(0, 2 * stats::pnorm(plof_ref[ok_loop] / nplof) - 1)
  }
  list(
    n = n,
    k = k_use,
    tau = tau,
    center = colMeans(x_unit),
    weight = weight,
    x_ref = x_unit,
    xw_ref = xw_ref,
    kdist_ref = kdist_ref,
    lrd_ref = lrd_ref,
    lof_ref = lof_ref,
    pdist_ref = pdist_ref,
    plof_ref = plof_ref,
    nplof = nplof,
    loop_ref = loop_ref
  )
}

.pca_ancestry_lof_score_query <- function(x_query, unit_model) {
  if (is.null(unit_model) || !nrow(x_query)) {
    return(list(
      lof = rep(NA_real_, nrow(x_query)),
      loop = rep(NA_real_, nrow(x_query))
    ))
  }
  k_use <- as_int(unit_model$k)[1]
  if (!is.finite(k_use) || is.na(k_use) || k_use < 1L) {
    return(list(
      lof = rep(NA_real_, nrow(x_query)),
      loop = rep(NA_real_, nrow(x_query))
    ))
  }
  topk <- .pca_ancestry_topk_query(
    xw_query = sweep(x_query, 2L, sqrt(unit_model$weight), "*", check.margin = FALSE),
    xw_ref = unit_model$xw_ref,
    kmax = min(k_use, nrow(unit_model$xw_ref))
  )
  out_lof <- rep(NA_real_, nrow(x_query))
  out_loop <- rep(NA_real_, nrow(x_query))
  for (i in seq_len(nrow(x_query))) {
    idx <- as_int(topk$idx[i, , drop = TRUE])
    dist <- as_num(topk$dist[i, , drop = TRUE])
    ok <- !is.na(idx) & is.finite(dist) & !is.na(dist)
    if (!any(ok)) next
    idx <- idx[ok]
    dist <- dist[ok]
    lrd_nb <- unit_model$lrd_ref[idx]
    kdist_nb <- unit_model$kdist_ref[idx]
    ok2 <- is.finite(lrd_nb) & !is.na(lrd_nb) & lrd_nb > 0 & is.finite(kdist_nb) & !is.na(kdist_nb)
    if (!any(ok2)) next
    lrd_nb <- lrd_nb[ok2]
    kdist_nb <- kdist_nb[ok2]
    dist <- dist[ok2]
    idx <- idx[ok2]
    reach <- pmax(dist, kdist_nb)
    if (!all(is.finite(reach) & !is.na(reach) & reach > 0)) next
    lrd_q <- 1 / mean(reach)
    if (!is.finite(lrd_q) || is.na(lrd_q) || lrd_q <= 0) next
    out_lof[i] <- mean(lrd_nb / lrd_q)

    pdist_q <- unit_model$tau * sqrt(mean(dist^2))
    pdist_nb <- unit_model$pdist_ref[idx]
    ok3 <- is.finite(pdist_nb) & !is.na(pdist_nb) & pdist_nb > 0
    if (!is.finite(pdist_q) || is.na(pdist_q) || pdist_q <= 0 || !any(ok3)) next
    plof_q <- pdist_q / mean(pdist_nb[ok3]) - 1
    nplof <- as_num(unit_model$nplof)[1]
    if (!is.finite(plof_q) || is.na(plof_q) || !is.finite(nplof) || is.na(nplof) || nplof <= 0) next
    out_loop[i] <- pmax(0, 2 * stats::pnorm(plof_q / nplof) - 1)
  }
  list(lof = out_lof, loop = out_loop)
}

.pca_ancestry_outlier_group_model <- function(x_group,
                                              subgroup_group = NULL,
                                              min_subgroup_n = 3L,
                                              weight = NULL,
                                              k = 20L) {
  n <- nrow(x_group)
  p <- ncol(x_group)
  if (n < 2L || p < 1L) {
    return(list(
      mode = "none",
      units = list(),
      weight = as_num(weight),
      k = as_int(k)[1]
    ))
  }

  subgroup_clean <- rep(NA_character_, n)
  if (!is.null(subgroup_group)) {
    subgroup_clean <- as.character(subgroup_group)
    subgroup_clean[is.na(subgroup_clean) | !nzchar(trimws(subgroup_clean))] <- NA_character_
  }

  subgroup_mode <- FALSE
  unit_id <- rep(".GROUP", n)
  if (any(!is.na(subgroup_clean))) {
    tab <- table(subgroup_clean, useNA = "no")
    if (length(tab) && all(!is.na(subgroup_clean)) && all(tab >= min_subgroup_n)) {
      subgroup_mode <- TRUE
      unit_id <- subgroup_clean
    }
  }

  units_idx <- split(seq_len(n), unit_id)
  units <- vector("list", length(units_idx))
  names(units) <- names(units_idx)
  for (nm in names(units_idx)) {
    idx <- as_int(units_idx[[nm]])
    x_unit <- x_group[idx, , drop = FALSE]
    unit_model <- .pca_ancestry_lof_unit_model(x_unit = x_unit, weight = weight, k = k)
    if (is.null(unit_model)) next
    units[[nm]] <- c(list(idx = idx), unit_model)
  }
  units <- Filter(Negate(is.null), units)

  list(
    mode = if (isTRUE(subgroup_mode)) "subgroup" else "group",
    units = units,
    weight = as_num(weight),
    k = as_int(k)[1]
  )
}

.pca_ancestry_build_outlier_model <- function(x_ref,
                                              label_ref,
                                              subgroup_ref = NULL,
                                              min_subgroup_n = 3L,
                                              weight = NULL,
                                              k = 20L) {
  label_ref <- as.character(label_ref)
  groups <- sort(unique(label_ref[!is.na(label_ref) & nzchar(label_ref)]))
  subgroup_ref <- if (is.null(subgroup_ref)) rep(NA_character_, length(label_ref)) else as.character(subgroup_ref)
  subgroup_ref[is.na(subgroup_ref) | !nzchar(trimws(subgroup_ref))] <- NA_character_

  models <- vector("list", length(groups))
  names(models) <- groups
  for (g in groups) {
    idx <- which(label_ref == g)
    models[[g]] <- .pca_ancestry_outlier_group_model(
      x_group = x_ref[idx, , drop = FALSE],
      subgroup_group = subgroup_ref[idx],
      min_subgroup_n = min_subgroup_n,
      weight = weight,
      k = k
    )
  }

  list(
    groups = groups,
    models = models
  )
}

.pca_ancestry_outlier_score <- function(x_query, label_assign, model) {
  label_assign <- as.character(label_assign)
  n <- nrow(x_query)
  out_lof <- rep(NA_real_, n)
  out_loop <- rep(NA_real_, n)
  grp_query <- unique(label_assign[!is.na(label_assign) & nzchar(label_assign)])
  for (g in grp_query) {
    mod_g <- model$models[[g]]
    if (is.null(mod_g) || !length(mod_g$units)) next
    idx <- which(label_assign == g)
    if (!length(idx)) next
    d_mat <- matrix(NA_real_, nrow = length(idx), ncol = length(mod_g$units))
    colnames(d_mat) <- names(mod_g$units)
    unit_no <- 0L
    for (nm in names(mod_g$units)) {
      unit_no <- unit_no + 1L
      unit <- mod_g$units[[nm]]
      d_mat[, unit_no] <- .pca_ancestry_diag_distance_matrix(
        x_query[idx, , drop = FALSE],
        unit$center,
        var = NULL,
        weight = mod_g$weight
      )
    }
    best_no <- max.col(-d_mat, ties.method = "first")
    res_g <- lapply(seq_len(length(idx)), function(ii) {
      unit_nm <- names(mod_g$units)[best_no[ii]]
      .pca_ancestry_lof_score_query(
        x_query = x_query[idx[ii], , drop = FALSE],
        unit_model = mod_g$units[[unit_nm]]
      )
    })
    out_lof[idx] <- vapply(res_g, function(z) as_num(z$lof)[1], numeric(1))
    out_loop[idx] <- vapply(res_g, function(z) as_num(z$loop)[1], numeric(1))
  }
  list(lof = out_lof, loop = out_loop)
}

.pca_ancestry_outlier_ref_max <- function(model) {
  vals <- unlist(
    lapply(model$models, function(m) {
      if (is.null(m) || !length(m$units)) return(numeric())
      unlist(lapply(m$units, function(u) as_num(u$lof_ref)), use.names = FALSE)
    }),
    use.names = FALSE
  )
  vals <- as_num(vals)
  vals <- vals[is.finite(vals) & !is.na(vals)]
  if (!length(vals)) return(NA_real_)
  max(vals)
}

.pca_ancestry_outlier_ref_score <- function(label_ref, model) {
  label_ref <- as.character(label_ref)
  out_lof <- rep(NA_real_, length(label_ref))
  out_loop <- rep(NA_real_, length(label_ref))
  grp_ref <- unique(label_ref[!is.na(label_ref) & nzchar(label_ref)])
  for (g in grp_ref) {
    idx_g <- which(label_ref == g)
    mod_g <- model$models[[g]]
    if (!length(idx_g) || is.null(mod_g) || !length(mod_g$units)) next
    score_lof_g <- rep(NA_real_, length(idx_g))
    score_loop_g <- rep(NA_real_, length(idx_g))
    for (unit in mod_g$units) {
      idx_u <- as_int(unit$idx)
      lof_u <- as_num(unit$lof_ref)
      loop_u <- as_num(unit$loop_ref)
      ok <- idx_u >= 1L & idx_u <= length(score_lof_g)
      if (!any(ok)) next
      score_lof_g[idx_u[ok]] <- lof_u[ok]
      score_loop_g[idx_u[ok]] <- loop_u[ok]
    }
    out_lof[idx_g] <- score_lof_g
    out_loop[idx_g] <- score_loop_g
  }
  list(lof = out_lof, loop = out_loop)
}

.pca_ancestry_binom_ci <- function(m, n, conf = 0.95) {
  m <- as_int(m)[1]
  n <- as_int(n)[1]
  conf <- as_num(conf)[1]
  if (!is.finite(n) || is.na(n) || n <= 0L) return(c(NA_real_, NA_real_))
  alpha <- 1 - conf
  lower <- if (m <= 0L) 0 else stats::qbeta(alpha / 2, m, n - m + 1)
  upper <- if (m >= n) 1 else stats::qbeta(1 - alpha / 2, m + 1, n - m)
  c(lower, upper)
}

.pca_ancestry_format_rate <- function(err, ci) {
  if (!is.finite(err) || is.na(err) || length(ci) != 2L || any(!is.finite(ci) | is.na(ci))) {
    return("NA [NA, NA]")
  }
  sprintf("%.3f [%.3f, %.3f]", err, ci[1], ci[2])
}

.pca_ancestry_safe_pp_cols <- function(groups) {
  groups <- as.character(groups)
  safe <- gsub("[^A-Za-z0-9._-]+", ".", groups)
  safe[is.na(safe) | !nzchar(safe)] <- "group"
  safe <- make.unique(safe, sep = ".")
  out <- paste0("PP.", safe)
  names(out) <- as.character(groups)
  out
}

.pca_ancestry_build_table <- function(row_id,
                                      ids,
                                      id_col_use,
                                      label_input,
                                      label_max,
                                      estimated,
                                      score_mat,
                                      groups,
                                      keep = NULL,
                                      reference = NULL,
                                      score_outlier = NULL) {
  dt <- data.table::data.table(.ROW_ID = as_int(row_id))
  dt[, (id_col_use) := as.character(ids)]
  dt[, estimated := as.character(estimated)]
  if (!is.null(keep)) dt[, keep := as.logical(keep)]
  if (!is.null(reference)) dt[, reference := as.logical(reference)]
  dt[, label.input := as.character(label_input)]
  dt[, label.max := as.character(label_max)]
  pp_max <- rep(NA_real_, nrow(dt))
  if (!is.null(score_mat) && length(score_mat)) {
    pp_max <- apply(score_mat, 1L, function(x) {
      x <- as_num(x)
      x <- x[is.finite(x) & !is.na(x)]
      if (!length(x)) return(NA_real_)
      max(x)
    })
  }
  dt[, PP.max := as_num(pp_max)]
  if (!is.null(score_outlier)) dt[, score.outlier := as_num(score_outlier)]
  rs_cols <- .pca_ancestry_safe_pp_cols(groups)
  for (g in groups) {
    dt[, (rs_cols[g]) := as_num(score_mat[, g])]
  }
  dt
}

# ----- PCA reference QC -----
# Reference ancestry QC on projected/reference PCs.
# - `pc` can be a `pca()` / `pca.projection()` return object or a data.frame/data.table
#   containing ID and PC columns.
# - `label.data` / `label.col` follow the same merge policy as `pca.plot(group.data, group.col)`.
# - Returns:
#   * `posterior`: per-sample posterior table with `keep`
#   * `id.drop`: IDs failing QC
#   * `param.best`: selected `k`
# Example:
#   ref.qc <- pca.refqc(
#     pc = pca.ref,
#     eval = pca.ref$eigenval,
#     pc.use = 1:10,
#     id.col = "IID",
#     label.data = sample.meta,
#     label.col = "ancestry"
#   )
#' Reference-PCA QC for ancestry estimation
#'
#' Removes within-population outliers from a reference PCA result so that
#' downstream [pca.ancestry()] estimates are not biased by mislabeled or
#' admixed samples in the reference. Selects the QC depth `k` adaptively.
#'
#' @param pc A reference PCA result (output of [pca()]).
#' @param eval Optional eigenvalue vector (defaults to `pc$eigenval`).
#' @param pc.use Integer vector of PCs to use (default `1:10`).
#' @param id.col Column name carrying the sample id.
#' @param label.data Optional `data.frame` of labels (one row per sample).
#' @param label.col Column name carrying the ancestry/population label.
#' @param k.grid Integer vector of `k` values to scan adaptively.
#' @param pve.cutoff,pp.cutoff Posterior thresholds for assignment / inclusion.
#' @param prior Prior mode (`"empirical"` or `"uniform"`).
#' @param prior.user Optional named prior overrides.
#' @param loop Logical. Iterate QC until convergence.
#' @param outlier.k,outlier.score Outlier-detection neighborhood size and score.
#'
#' @return A list with the QC-passed PC table, the dropped IDs, and the
#'   selected `param.best` value.
#' @export
pca.refqc <- function(pc,
                      eval = NULL,
                      pc.use = 1:10,
                      id.col = "IID",
                      label.data = NULL,
                      label.col = NULL,
                      k.grid = seq(11, 61, by = 2),
                      pve.cutoff = 0.95,
                      pp.cutoff = 0.5,
                      prior = "empirical",
                      prior.user = NULL,
                      loop = TRUE,
                      outlier.k = 11L,
                      outlier.score = 0.95) {
  require_pkg("data.table")
  conf <- 0.95
  pve.cutoff <- as_num(pve.cutoff)[1]
  pp.cutoff <- as_num(pp.cutoff)[1]
  loop <- isTRUE(loop)
  outlier.k <- as_int(outlier.k)[1]
  outlier.score_in <- outlier.score
  outlier.score <- if (is.null(outlier.score_in) || (is.character(outlier.score_in) && length(outlier.score_in) == 1L && tolower(trimws(outlier.score_in[1])) == "auto")) NA_real_ else as_num(outlier.score_in)[1]
  if (is.na(outlier.score)) outlier.score <- 0.95
  if (!is.finite(pve.cutoff) || is.na(pve.cutoff) || pve.cutoff <= 0 || pve.cutoff > 1) stop("pve.cutoff must be in (0, 1].", call. = FALSE)
  if (!is.finite(pp.cutoff) || is.na(pp.cutoff) || pp.cutoff <= 0 || pp.cutoff > 1) stop("pp.cutoff must be in (0, 1].", call. = FALSE)
  if (!is.finite(outlier.k) || is.na(outlier.k) || outlier.k < 1L) stop("outlier.k must be >= 1.", call. = FALSE)

  inp <- .pca_ancestry_resolve_input(pc = pc, eval = eval, pc.use = pc.use, id.col = id.col)
  dt <- .pca_ancestry_attach_label_data(inp$dt, id_col_use = inp$id_col_use, label.data = label.data, label.col = label.col)
  ref_idx_label <- which(!is.na(dt$.LABEL_INPUT) & nzchar(dt$.LABEL_INPUT))
  use_label_gate <- length(ref_idx_label) >= 2L
  if (use_label_gate) {
    ref_idx <- ref_idx_label
    labels_ref <- as.character(dt$.LABEL_INPUT[ref_idx])
  } else {
    ref_idx <- seq_len(nrow(dt))
    if (length(ref_idx) < 2L) stop("pca.refqc requires at least two samples.", call. = FALSE)
    labels_ref <- rep("Reference", length(ref_idx))
  }

  groups <- sort(unique(labels_ref))
  post_pc <- .pca_ancestry_posterior_pc_spec(inp$pc_names, inp$eval_named, min_pc_reduce = 5L, prop_target = pve.cutoff, pc_names_all = inp$pc_names_all)
  x_ref <- inp$x[ref_idx, post_pc$pc_idx, drop = FALSE]
  outlier_pc <- .pca_ancestry_outlier_pc_spec(inp$pc_names, inp$eval_named, min_pc_reduce = 5L, prop_target = pve.cutoff, pc_names_all = inp$pc_names_all)
  x_ref_outlier <- inp$x[ref_idx, outlier_pc$pc_idx, drop = FALSE]
  scale_fit <- .pca_ancestry_fit_scale(x_ref)
  weight <- post_pc$weight
  xw_ref <- .pca_ancestry_transform(x_ref, center = scale_fit$center, scale = scale_fit$scale, weight = weight)

  k_grid <- as_int(k.grid)
  k_grid <- sort(unique(k_grid[is.finite(k_grid) & !is.na(k_grid) & k_grid >= 1L]))
  if (!length(k_grid)) stop("k.grid must contain at least one positive integer.", call. = FALSE)
  k_grid <- unique(pmin(k_grid, length(ref_idx) - 1L))
  k_grid <- k_grid[k_grid >= 1L]
  if (!length(k_grid)) stop("k.grid is invalid after clipping to n_reference - 1.", call. = FALSE)

  knn_ref_group_pos <- .pca_ancestry_group_topk_self(xw_ref, labels_ref = labels_ref, groups = groups, kmax = max(k_grid))
  knn_ref_group_neg <- .pca_ancestry_complement_topk_self(xw_ref, labels_ref = labels_ref, groups = groups, kmax = max(k_grid))

  .gcanvas_note(
    "gcanvas::pca.refqc",
    sprintf("Start: n=%d | p=%d | groups=%d", as_int(length(ref_idx)), as_int(length(post_pc$pc_names)), as_int(length(groups)))
  )
  .gcanvas_note(
    "gcanvas::pca.refqc",
    sprintf("Posterior PCs=%s | prop.var=%.3f", post_pc$label, as_num(post_pc$cum_prop)[1])
  )
  .gcanvas_note(
    "gcanvas::pca.refqc",
    sprintf("Grid search: k.grid.n=%d | prior=%s", as_int(length(k_grid)), .pca_ancestry_as_char1(prior))
  )
  outlier_model <- NULL
  score_outlier_ref <- list(
    lof = rep(NA_real_, length(ref_idx)),
    loop = rep(NA_real_, length(ref_idx))
  )
  outlier_score_apply <- outlier.score
  if (isTRUE(loop)) {
    outlier_model <- .pca_ancestry_build_outlier_model(
      x_ref = x_ref_outlier,
      label_ref = labels_ref,
      subgroup_ref = NULL,
      weight = outlier_pc$weight,
      k = outlier.k
    )
    score_outlier_ref <- .pca_ancestry_outlier_ref_score(
      label_ref = labels_ref,
      model = outlier_model
    )
    .gcanvas_note(
      "gcanvas::pca.refqc",
      sprintf(
        "Density gate: PCs=%s | prop.var=%.3f | outlier.k=%d | outlier.score=%s",
        outlier_pc$label,
        as_num(outlier_pc$cum_prop)[1],
        as_int(outlier.k),
        if (is.na(outlier_score_apply)) "NA" else format(signif(outlier_score_apply, 6), scientific = FALSE, trim = TRUE)
      )
    )
  }

  if (!use_label_gate || length(groups) <= 1L) {
    group1 <- as.character(groups[1])
    if (!nzchar(group1)) group1 <- "Reference"
    pp_mat <- matrix(1, nrow = length(ref_idx), ncol = 1L)
    colnames(pp_mat) <- group1
    hard_drop <- rep(FALSE, length(ref_idx))
    drop_reason <- rep(NA_character_, length(ref_idx))
    if (isTRUE(loop)) {
      fail_loop <- !is.na(outlier_score_apply) & !is.na(score_outlier_ref$loop) & (score_outlier_ref$loop > outlier_score_apply)
      hard_drop <- fail_loop
      drop_reason[fail_loop] <- "LoOP_density_outlier"
    }
    ref_estimated <- rep(group1, length(ref_idx))
    ref_estimated[hard_drop] <- "Unknown"
    ref_dt <- .pca_ancestry_build_table(
      row_id = dt$.ROW_ID[ref_idx],
      ids = dt$IID_USE[ref_idx],
      id_col_use = inp$id_col_use,
      label_input = labels_ref,
      label_max = rep(group1, length(ref_idx)),
      estimated = ref_estimated,
      score_mat = pp_mat,
      groups = group1,
      keep = !hard_drop,
      score_outlier = score_outlier_ref$loop
    )
    data.table::setorderv(ref_dt, ".ROW_ID")
    ref_dt[, PP.input := rep(1, .N)]
    ref_dt[, PP.diff := rep(0, .N)]
    ref_dt[, ambiguous := FALSE]
    ref_dt[, REASON := drop_reason]
    data.table::setcolorder(
      ref_dt,
      c(
        inp$id_col_use,
        "estimated",
        "keep",
        "ambiguous",
        "REASON",
        "label.input",
        "label.max",
        "PP.input",
        "PP.max",
        "PP.diff",
        setdiff(names(ref_dt), c(".ROW_ID", inp$id_col_use, "estimated", "keep", "ambiguous", "REASON", "label.input", "label.max", "PP.input", "PP.max", "PP.diff"))
      )
    )
    id_drop <- as.character(ref_dt[get("keep") %in% FALSE, get(inp$id_col_use)])
    id_ambiguous <- character(0)
    drop_info <- ref_dt[get("keep") %in% FALSE, .(IID = as.character(get(inp$id_col_use)), REASON = as.character(REASON))]
    ref_dt[, .ROW_ID := NULL]
    .gcanvas_note(
      "gcanvas::pca.refqc",
      sprintf(
        "Shortcut: label gate skipped (%s); assigned PP.%s=1 for all samples%s",
        if (!use_label_gate) "no/insufficient group labels" else "single group",
        group1,
        if (isTRUE(loop)) ", then applied LoOP gate" else ""
      )
    )
    .gcanvas_note(
      "gcanvas::pca.refqc",
      sprintf(
        "Done: drop=%d | ambiguous=0 | error(before)=NA | error(after)=NA | k=%d",
        as_int(length(id_drop)),
        as_int(min(k_grid))
      )
    )
    dt_plot <- data.table::data.table(
      IID_USE = dt$IID_USE[ref_idx],
      .REFQC_GROUP = labels_ref
    )
    dt_plot[, (inp$id_col_use) := dt[[inp$id_col_use]][ref_idx]]
    for (nm in inp$pc_names_all) {
      dt_plot[, (nm) := dt[[nm]][ref_idx]]
    }
    data.table::setcolorder(dt_plot, c(inp$id_col_use, "IID_USE", ".REFQC_GROUP", inp$pc_names_all))
    out <- list(
      posterior = ref_dt,
      id.drop = id_drop,
      drop.info = drop_info,
      id.ambiguous = id_ambiguous,
      pc = dt_plot,
      eval = inp$eval_named,
      group.title = .pca_ancestry_as_char1(label.col) %||% "Group",
      pc.use = as_int(post_pc$pc_idx),
      param.best = list(
        k = as_int(min(k_grid)),
        pp.cutoff = as_num(pp.cutoff),
        prior = .pca_ancestry_as_char1(prior),
        loop = isTRUE(loop),
        outlier.k = as_int(outlier.k),
        outlier.score = as_num(outlier_score_apply)
      )
    )
    class(out) <- c("gcanvas_pca_refqc", "list")
    return(out)
  }

  best <- NULL

  for (k0 in k_grid) {
    loo <- .pca_ancestry_indep_posterior(
      group_dist_pos = knn_ref_group_pos,
      group_dist_neg = knn_ref_group_neg,
      groups = groups,
      labels_ref = labels_ref,
      k = k0,
      dim_p = ncol(x_ref),
      prior = prior,
      prior.user = prior.user
    )
    pp_max <- apply(loo$score, 1L, function(x) {
      x <- as_num(x)
      x <- x[is.finite(x) & !is.na(x)]
      if (!length(x)) return(NA_real_)
      max(x)
    })
    pp_own <- rep(NA_real_, length(labels_ref))
    for (gg in groups) {
      hit <- which(labels_ref == gg)
      if (!length(hit)) next
      pp_own[hit] <- as_num(loo$score[hit, gg])
    }
    pp_diff <- pp_max - pp_own
    est <- loo$label.max
    est[is.na(est) | !nzchar(est) | is.na(pp_max) | (pp_max < pp.cutoff)] <- "Unknown"
    gate1_drop <- is.na(est) |
      is.na(pp_own) |
      (pp_own < pp.cutoff) |
      (
        (est != labels_ref) &
          (is.na(pp_diff) | (pp_diff > 0.1))
      )
    gate2_drop <- rep(FALSE, length(labels_ref))
    fail_loop <- rep(FALSE, length(labels_ref))
    if (isTRUE(loop)) {
      fail_loop <- !is.na(outlier_score_apply) & !is.na(score_outlier_ref$loop) & (score_outlier_ref$loop > outlier_score_apply)
      gate2_drop <- fail_loop & !gate1_drop
    }
    hard_drop <- gate1_drop | gate2_drop
    drop_reason <- rep(NA_character_, length(labels_ref))
    drop_reason[gate1_drop] <- "label_consistency"
    drop_reason[gate2_drop & fail_loop] <- "LoOP_density_outlier"
    soft_flag <- !hard_drop & (est != labels_ref)
    miss <- hard_drop
    err <- mean(miss, na.rm = TRUE)
    ci <- .pca_ancestry_binom_ci(sum(miss, na.rm = TRUE), length(labels_ref), conf = conf)
    cand <- list(
      k = k0,
      loo = loo,
      err = err,
      ci = ci,
      keep = !miss,
      ambiguous = soft_flag,
      est = est,
      pp_max = pp_max,
      pp_own = pp_own,
      pp_diff = pp_diff,
      drop_reason = drop_reason,
      score_outlier_loop = score_outlier_ref$loop
    )
    if (is.null(best) ||
        cand$err < best$err - 1e-12 ||
        (abs(cand$err - best$err) <= 1e-12 && cand$k < best$k)) {
      best <- cand
    }
  }

  if (is.null(best)) stop("Failed to select best k.", call. = FALSE)

  ref_estimated <- best$est
  ref_estimated[is.na(ref_estimated) | !nzchar(ref_estimated)] <- "Unknown"
  ref_dt <- .pca_ancestry_build_table(
    row_id = dt$.ROW_ID[ref_idx],
    ids = dt$IID_USE[ref_idx],
    id_col_use = inp$id_col_use,
    label_input = labels_ref,
    label_max = best$loo$label.max,
      estimated = ref_estimated,
      score_mat = best$loo$score,
      groups = groups,
      keep = best$keep,
      score_outlier = best$score_outlier_loop
    )
  data.table::setorderv(ref_dt, ".ROW_ID")
  ref_dt[, PP.input := as_num(best$pp_own)]
  ref_dt[, PP.diff := as_num(best$pp_diff)]
  ref_dt[, ambiguous := as.logical(best$ambiguous)]
  ref_dt[, REASON := as.character(best$drop_reason)]
  data.table::setcolorder(
    ref_dt,
    c(
      inp$id_col_use,
      "estimated",
      "keep",
      "ambiguous",
      "REASON",
      "label.input",
      "label.max",
      "PP.input",
      "PP.max",
      "PP.diff",
      setdiff(names(ref_dt), c(".ROW_ID", inp$id_col_use, "estimated", "keep", "ambiguous", "REASON", "label.input", "label.max", "PP.input", "PP.max", "PP.diff"))
    )
  )
  id_drop <- as.character(ref_dt[get("keep") %in% FALSE, get(inp$id_col_use)])
  id_ambiguous <- as.character(ref_dt[get("ambiguous") %in% TRUE & get("keep") %in% TRUE, get(inp$id_col_use)])
  drop_info <- ref_dt[get("keep") %in% FALSE, .(IID = as.character(get(inp$id_col_use)), REASON = as.character(REASON))]
  ref_dt[, .ROW_ID := NULL]
  n_ambiguous <- sum(ref_dt$ambiguous %in% TRUE, na.rm = TRUE)

  after_keep <- sum(best$keep, na.rm = TRUE)
  after_ci <- .pca_ancestry_binom_ci(0, after_keep, conf = conf)
  before_txt <- .pca_ancestry_format_rate(best$err, best$ci)
  after_txt <- .pca_ancestry_format_rate(0, after_ci)
  msg_tail <- sprintf(
    "Done: drop=%d | ambiguous=%d | error(before)=%s | error(after)=%s | k=%d",
    as_int(length(id_drop)),
    as_int(n_ambiguous),
    before_txt,
    after_txt,
    as_int(best$k)
  )
  .gcanvas_note("gcanvas::pca.refqc", msg_tail)

  dt_plot <- data.table::data.table(
    IID_USE = dt$IID_USE[ref_idx],
    .REFQC_GROUP = labels_ref
  )
  dt_plot[, (inp$id_col_use) := dt[[inp$id_col_use]][ref_idx]]
  for (nm in inp$pc_names_all) {
    dt_plot[, (nm) := dt[[nm]][ref_idx]]
  }
  data.table::setcolorder(dt_plot, c(inp$id_col_use, "IID_USE", ".REFQC_GROUP", inp$pc_names_all))
  out <- list(
    posterior = ref_dt,
    id.drop = id_drop,
    drop.info = drop_info,
    id.ambiguous = id_ambiguous,
    pc = dt_plot,
    eval = inp$eval_named,
    group.title = .pca_ancestry_as_char1(label.col) %||% "Group",
    pc.use = as_int(post_pc$pc_idx),
      param.best = list(
        k = as_int(best$k),
        pp.cutoff = as_num(pp.cutoff),
        prior = .pca_ancestry_as_char1(prior),
        loop = isTRUE(loop),
        outlier.k = as_int(outlier.k),
        outlier.score = as_num(outlier_score_apply)
      )
  )
  class(out) <- c("gcanvas_pca_refqc", "list")
  out
}


# ----- PCA ancestry estimation -----
# Ancestry posterior estimation for reference + target samples in the same projected PC space.
# - Reference/target are inferred in this order:
#   `target.id` > `label.ref` > `reference` column in `pc`/`label.data` > missing label rows as target.
# - Posterior uses group-wise local weighted-kNN on the selected PCs.
# - Outlier gate uses assigned-group LoOP on reduced outlier PCs.
#   If subgroup information is available, the nearest subgroup is selected first and LoOP is
#   computed within that subgroup only.
#   It reports `score.outlier = LoOP`.
# - Returns:
#   * `posterior`: posterior table for reference/target rows
#   * `outlier.score`: applied LoOP cutoff
# Example:
#   anc <- pca.ancestry(
#     pc = pca.proj,
#     eval = pca.ref$eigenval,
#     pc.use = 1:10,
#     id.col = "IID",
#     label.data = sample.meta,
#     label.col = "ancestry",
#     subgroup.data = sample.meta,
#     subgroup.col = "subgroup",
#     k = ref.qc$param.best$k,
#     outlier.k = 21,
#     outlier.score = NULL
#   )
# Example with explicit target IDs:
#   anc.target <- pca.ancestry(
#     pc = pca.proj,
#     eval = pca.ref$eigenval,
#     label.data = sample.meta,
#     label.col = "ancestry",
#     target.id = c("S1", "S2", "S3"),
#     k = ref.qc$param.best$k,
#     target.only = TRUE
#   )
#' Estimate ancestry membership from a reference-aligned PCA
#'
#' Assigns each target sample to (or scores its membership across) the
#' reference ancestry labels by comparing its position in PC space to the
#' QC'd reference populations.
#'
#' @param pc A PCA result containing both reference and target samples
#'   (typically the output of [pca.projection()]).
#' @param eval Optional eigenvalue vector (defaults to `pc$eigenval`).
#' @param pc.use Integer vector of PCs to use (default `1:10`).
#' @param id.col Column name carrying the sample id.
#' @param label.data A `data.frame` of `(id, ancestry)` for reference samples.
#' @param label.col Column name carrying the ancestry label in `label.data`.
#' @param subgroup.data,subgroup.col Optional finer-grained subgroup labels.
#' @param label.ref Optional fallback label for unmatched reference rows.
#' @param target.id Optional explicit target-sample ids (otherwise inferred
#'   from missing labels).
#' @param target.only Logical. Return only target-sample rows in the result.
#' @param k Neighborhood size used for posterior estimation.
#' @param pve.cutoff,pp.cutoff Posterior thresholds for assignment / inclusion.
#' @param prior Prior mode (`"empirical"` or `"uniform"`).
#' @param prior.user Optional named prior overrides.
#' @param outlier.k,outlier.score Outlier-detection neighborhood size and score.
#'
#' @return A list with per-sample ancestry assignments / probabilities and
#'   any QC diagnostics.
#' @export
pca.ancestry <- function(pc,
                         eval = NULL,
                         pc.use = 1:10,
                         id.col = "IID",
                         label.data = NULL,
                         label.col = NULL,
                         subgroup.data = NULL,
                         subgroup.col = NULL,
                         label.ref = NULL,
                         target.id = NULL,
                         target.only = FALSE,
                         k = 11L,
                         pve.cutoff = 0.95,
                         pp.cutoff = 0.5,
                         prior = "empirical",
                         prior.user = NULL,
                         outlier.k = 11L,
                         outlier.score = 0.95) {
  require_pkg("data.table")
  target.only <- isTRUE(target.only)
  pve.cutoff <- as_num(pve.cutoff)[1]
  pp.cutoff <- as_num(pp.cutoff)[1]
  conf <- 0.95
  outlier.k <- as_int(outlier.k)[1]
  outlier.score_in <- outlier.score
  outlier.score <- if (is.null(outlier.score_in) || (is.character(outlier.score_in) && length(outlier.score_in) == 1L && tolower(trimws(outlier.score_in[1])) == "auto")) NA_real_ else as_num(outlier.score_in)[1]
  if (is.na(outlier.score)) outlier.score <- 0.95
  k_in <- as_int(k)[1]
  if (!is.finite(pve.cutoff) || is.na(pve.cutoff) || pve.cutoff <= 0 || pve.cutoff > 1) stop("pve.cutoff must be in (0, 1].", call. = FALSE)
  if (!is.finite(pp.cutoff) || is.na(pp.cutoff) || pp.cutoff <= 0 || pp.cutoff > 1) stop("pp.cutoff must be in (0, 1].", call. = FALSE)
  if (!is.finite(outlier.k) || is.na(outlier.k) || outlier.k < 1L) stop("outlier.k must be >= 1.", call. = FALSE)
  if (!is.na(outlier.score) && (!is.finite(outlier.score) || outlier.score < 0)) stop("outlier.score must be >= 0.", call. = FALSE)
  if (!is.finite(k_in) || is.na(k_in) || k_in < 1L) stop("k must be >= 1.", call. = FALSE)

  inp <- .pca_ancestry_resolve_input(pc = pc, eval = eval, pc.use = pc.use, id.col = id.col)
  dt <- .pca_ancestry_attach_label_data(inp$dt, id_col_use = inp$id_col_use, label.data = label.data, label.col = label.col)
  dt <- .pca_ancestry_attach_subgroup_data(dt, id_col_use = inp$id_col_use, subgroup.data = subgroup.data, subgroup.col = subgroup.col)

  label_input <- as.character(dt$.LABEL_INPUT)
  label_input[is.na(label_input) | !nzchar(label_input)] <- NA_character_
  subgroup_input <- as.character(dt$.SUBGROUP_INPUT)
  subgroup_input[is.na(subgroup_input) | !nzchar(subgroup_input)] <- NA_character_
  ids <- as.character(dt$IID_USE)
  ref_flag <- as.logical(dt$.REFERENCE_FLAG)
  label_ref_use <- as.character(label.ref)
  label_ref_use <- label_ref_use[!is.na(label_ref_use) & nzchar(label_ref_use)]
  target_id_use <- as.character(target.id)
  target_id_use <- target_id_use[!is.na(target_id_use) & nzchar(target_id_use)]

  target_idx <- integer()
  ref_idx <- integer()
  excluded_idx <- integer()

  if (length(target_id_use)) {
    target_idx <- which(ids %in% target_id_use)
    if (!length(target_idx)) stop("No target.id matched id.col.", call. = FALSE)
    ref_mask <- !(seq_len(nrow(dt)) %in% target_idx) & !is.na(label_input)
    if (length(label_ref_use)) ref_mask <- ref_mask & (label_input %in% label_ref_use)
    ref_idx <- which(ref_mask)
    excluded_idx <- setdiff(seq_len(nrow(dt)), c(ref_idx, target_idx))
  } else if (length(label_ref_use)) {
    ref_idx <- which(!is.na(label_input) & (label_input %in% label_ref_use))
    target_idx <- which(is.na(label_input) | !(label_input %in% label_ref_use))
  } else if (any(!is.na(ref_flag))) {
    ref_idx <- which(ref_flag %in% TRUE & !is.na(label_input))
    target_idx <- which(!(ref_flag %in% TRUE & !is.na(label_input)))
  } else {
    ref_idx <- which(!is.na(label_input))
    target_idx <- which(is.na(label_input))
  }

  ref_idx <- sort(unique(ref_idx))
  target_idx <- sort(unique(target_idx))
  excluded_idx <- sort(unique(setdiff(excluded_idx, c(ref_idx, target_idx))))

  if (length(ref_idx) < 2L) stop("pca.ancestry requires at least two labeled reference samples.", call. = FALSE)

  if (length(excluded_idx)) {
    .gcanvas_note(
      "gcanvas::pca.ancestry",
      sprintf("Excluded rows outside reference/target rules: n=%d", as_int(length(excluded_idx)))
    )
  }

  groups <- sort(unique(label_input[ref_idx]))
  post_pc <- .pca_ancestry_posterior_pc_spec(inp$pc_names, inp$eval_named, min_pc_reduce = 5L, prop_target = pve.cutoff, pc_names_all = inp$pc_names_all)
  x_ref <- inp$x[ref_idx, post_pc$pc_idx, drop = FALSE]
  outlier_pc <- .pca_ancestry_outlier_pc_spec(inp$pc_names, inp$eval_named, min_pc_reduce = 5L, prop_target = pve.cutoff, pc_names_all = inp$pc_names_all)
  x_ref_outlier <- inp$x[ref_idx, outlier_pc$pc_idx, drop = FALSE]
  outlier_model <- .pca_ancestry_build_outlier_model(
    x_ref = x_ref_outlier,
    label_ref = label_input[ref_idx],
    subgroup_ref = subgroup_input[ref_idx],
    weight = outlier_pc$weight,
    k = outlier.k
  )
  score_outlier_ref <- .pca_ancestry_outlier_ref_score(
    label_ref = label_input[ref_idx],
    model = outlier_model
  )
  outlier_score_apply <- outlier.score
  scale_fit <- .pca_ancestry_fit_scale(x_ref)
  weight <- post_pc$weight
  xw_ref <- .pca_ancestry_transform(x_ref, center = scale_fit$center, scale = scale_fit$scale, weight = weight)

  k_ref <- min(k_in, length(ref_idx) - 1L)
  if (k_ref < 1L) stop("Reference size is too small for LOO posterior.", call. = FALSE)
  if (k_ref != k_in) {
    .gcanvas_note("gcanvas::pca.ancestry", sprintf("k clipped for reference LOO: %d -> %d", as_int(k_in), as_int(k_ref)))
  }
  k_target <- min(k_in, length(ref_idx))
  if (k_target < 1L) stop("Reference size is too small for target posterior.", call. = FALSE)

  .gcanvas_note(
    "gcanvas::pca.ancestry",
    sprintf("Start: n=%d | reference=%d | target=%d | outlier.k=%d | outlier.score=%s | prior=%s", as_int(nrow(dt)), as_int(length(ref_idx)), as_int(length(target_idx)), as_int(outlier.k), if (is.na(outlier_score_apply)) "NA" else format(signif(outlier_score_apply, 6), scientific = FALSE, trim = TRUE), .pca_ancestry_as_char1(prior))
  )
  .gcanvas_note(
    "gcanvas::pca.ancestry",
    sprintf(
      "Posterior PCs=%s | prop.var=%.3f | k=%d",
      post_pc$label,
      as_num(post_pc$cum_prop)[1],
      as_int(k_target)
    )
  )
  .gcanvas_note(
    "gcanvas::pca.ancestry",
    sprintf(
      "Outlier score: PCs=%s | prop.var=%.3f | k=%d | score=LoOP",
      outlier_pc$label,
      as_num(outlier_pc$cum_prop)[1],
      as_int(outlier.k)
    )
  )

  knn_ref_pos <- .pca_ancestry_group_topk_self(xw_ref, labels_ref = label_input[ref_idx], groups = groups, kmax = k_ref)
  knn_ref_neg <- .pca_ancestry_complement_topk_self(xw_ref, labels_ref = label_input[ref_idx], groups = groups, kmax = k_ref)
  loo_ref <- .pca_ancestry_indep_posterior(
    group_dist_pos = knn_ref_pos,
    group_dist_neg = knn_ref_neg,
    groups = groups,
    labels_ref = label_input[ref_idx],
    k = k_ref,
    dim_p = ncol(x_ref),
    prior = prior,
    prior.user = prior.user
  )
  ref_pp_max <- apply(loo_ref$score, 1L, function(x) {
    x <- as_num(x)
    x <- x[is.finite(x) & !is.na(x)]
    if (!length(x)) return(NA_real_)
    max(x)
  })
  ref_est <- loo_ref$label.max
  ref_est[is.na(ref_est) | !nzchar(ref_est) | is.na(ref_pp_max) | (ref_pp_max < pp.cutoff)] <- "Unknown"

  dt_out_list <- list()
  if (!target.only) {
    dt_ref <- .pca_ancestry_build_table(
      row_id = dt$.ROW_ID[ref_idx],
      ids = ids[ref_idx],
      id_col_use = inp$id_col_use,
      label_input = label_input[ref_idx],
      label_max = loo_ref$label.max,
      estimated = ref_est,
      score_mat = loo_ref$score,
      groups = groups,
      reference = TRUE,
      score_outlier = score_outlier_ref$loop
    )
    dt_out_list[[length(dt_out_list) + 1L]] <- dt_ref
  }

  unknown_n <- 0L
  unknown_assign_n <- 0L
  unknown_outlier_n <- 0L
  if (length(target_idx)) {
    x_target <- inp$x[target_idx, post_pc$pc_idx, drop = FALSE]
    x_target_outlier <- inp$x[target_idx, outlier_pc$pc_idx, drop = FALSE]
    xw_target <- .pca_ancestry_transform(x_target, center = scale_fit$center, scale = scale_fit$scale, weight = weight)
    knn_target_pos <- .pca_ancestry_group_topk_query(xw_target, xw_ref, labels_ref = label_input[ref_idx], groups = groups, kmax = k_target)
    knn_target_neg <- .pca_ancestry_complement_topk_query(xw_target, xw_ref, labels_ref = label_input[ref_idx], groups = groups, kmax = k_target)
    post_target <- .pca_ancestry_indep_posterior(
      group_dist_pos = knn_target_pos,
      group_dist_neg = knn_target_neg,
      groups = groups,
      labels_ref = label_input[ref_idx],
      k = k_target,
      dim_p = ncol(x_ref),
      prior = prior,
      prior.user = prior.user
    )
    pp_max_target <- apply(post_target$score, 1L, function(x) {
      x <- as_num(x)
      x <- x[is.finite(x) & !is.na(x)]
      if (!length(x)) return(NA_real_)
      max(x)
    })
    score_outlier_target <- .pca_ancestry_outlier_score(
      x_query = x_target_outlier,
      label_assign = post_target$label.max,
      model = outlier_model
    )
    assign_fail <- is.na(post_target$label.max) | !nzchar(post_target$label.max) | is.na(pp_max_target) | (pp_max_target < pp.cutoff)
    outlier_fail <- !assign_fail & !is.na(outlier_score_apply) & !is.na(score_outlier_target$loop) & (score_outlier_target$loop > outlier_score_apply)
    target_est <- post_target$label.max
    target_est[is.na(target_est) | !nzchar(target_est)] <- "Unknown"
    target_est[assign_fail | outlier_fail] <- "Unknown"
    unknown_n <- sum(target_est == "Unknown", na.rm = TRUE)
    unknown_assign_n <- sum(assign_fail, na.rm = TRUE)
    unknown_outlier_n <- sum(outlier_fail, na.rm = TRUE)
    dt_target <- .pca_ancestry_build_table(
      row_id = dt$.ROW_ID[target_idx],
      ids = ids[target_idx],
      id_col_use = inp$id_col_use,
      label_input = label_input[target_idx],
      label_max = post_target$label.max,
      estimated = target_est,
      score_mat = post_target$score,
      groups = groups,
      reference = FALSE,
      score_outlier = score_outlier_target$loop
    )
    dt_out_list[[length(dt_out_list) + 1L]] <- dt_target
  }

  if (!length(dt_out_list)) {
    posterior_dt <- data.table::data.table()
  } else {
    posterior_dt <- data.table::rbindlist(dt_out_list, use.names = TRUE, fill = TRUE)
    data.table::setorderv(posterior_dt, ".ROW_ID")
    posterior_dt[, .ROW_ID := NULL]
  }

  ref_miss <- is.na(ref_est) | (ref_est != label_input[ref_idx])
  ref_n <- length(ref_idx)
  ref_err <- if (ref_n > 0L) sum(ref_miss, na.rm = TRUE) / ref_n else NA_real_
  ref_ci <- .pca_ancestry_binom_ci(
    sum(ref_miss, na.rm = TRUE),
    ref_n,
    conf = conf
  )
  msg_done <- sprintf(
    "Done: reference error=%s | unknown=%d/%d (PP.max<%s=%d, outlier=%d) | k=%d | outlier.k=%d | outlier.score=%s",
    .pca_ancestry_format_rate(ref_err, ref_ci),
    as_int(unknown_n),
    as_int(length(target_idx)),
    format(signif(pp.cutoff, 6), scientific = FALSE, trim = TRUE),
    as_int(unknown_assign_n),
    as_int(unknown_outlier_n),
    as_int(k_ref),
    as_int(outlier.k),
    if (is.na(outlier_score_apply)) "NA" else format(signif(outlier_score_apply, 6), scientific = FALSE, trim = TRUE)
  )
  .gcanvas_note("gcanvas::pca.ancestry", msg_done)

  list(
    posterior = posterior_dt,
    pc.use = as_int(post_pc$pc_idx),
    outlier.score = as_num(outlier_score_apply)
  )
}

