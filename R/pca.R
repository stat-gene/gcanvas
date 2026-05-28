# PCA toolkit: compute principal components from genotype data (`pca()`),
# project new samples onto a reference PCA (`pca.projection()`), and visualize
# the result (`pca.plot()`, `pca.screeplot()`).

.gcanvas_pca_pc_cols <- function(nm) {
  nm <- as.character(nm)
  hit <- grep("^PC[0-9]+$", nm, ignore.case = FALSE, value = TRUE)
  hit[order(as_int(sub("^PC", "", hit)))]
}

.gcanvas_pca_parse_pc_use <- function(x, pc_cols) {
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

.gcanvas_pca_as_pc_dt <- function(x) {
  if (is.null(x)) return(NULL)
  if (data.table::is.data.table(x)) return(data.table::copy(x))
  if (is.data.frame(x)) return(data.table::as.data.table(x))
  if (is.matrix(x)) return(data.table::as.data.table(as.data.frame(x, stringsAsFactors = FALSE)))
  NULL
}

.gcanvas_pca_extract_input <- function(x) {
  dt <- NULL
  eval_use <- NULL
  if (is.list(x) && !data.table::is.data.table(x) && !is.data.frame(x)) {
    cand_keys <- c("pc", "PC", "data", "scores", "eigenvec")
    for (kk in cand_keys) {
      if (!is.null(x[[kk]])) {
        dt0 <- .gcanvas_pca_as_pc_dt(x[[kk]])
        if (!is.null(dt0)) {
          dt <- dt0
          break
        }
      }
    }
    if (is.null(dt) && length(x)) {
      for (ii in seq_along(x)) {
        dt0 <- .gcanvas_pca_as_pc_dt(x[[ii]])
        if (is.null(dt0)) next
        pcc0 <- .gcanvas_pca_pc_cols(names(dt0))
        if (length(pcc0) >= 1L) {
          dt <- dt0
          break
        }
      }
    }
    if (!is.null(x$eigenvalue) && length(x$eigenvalue)) eval_use <- x$eigenvalue
    else if (!is.null(x$prop.var) && length(x$prop.var)) eval_use <- x$prop.var
    else if (!is.null(x$eval) && length(x$eval)) eval_use <- x$eval
    else if (!is.null(attr(x, "eigenvalue")) && length(attr(x, "eigenvalue"))) eval_use <- attr(x, "eigenvalue")
    else if (!is.null(attr(x, "prop.var")) && length(attr(x, "prop.var"))) eval_use <- attr(x, "prop.var")
  } else {
    dt <- .gcanvas_pca_as_pc_dt(x)
    if (!is.null(attr(x, "eigenvalue")) && length(attr(x, "eigenvalue"))) eval_use <- attr(x, "eigenvalue")
    else if (!is.null(attr(x, "prop.var")) && length(attr(x, "prop.var"))) eval_use <- attr(x, "prop.var")
  }
  list(pc = dt, eval = eval_use)
}

.gcanvas_pca_normalize_eval <- function(eval_in, pc_cols = NULL) {
  if (is.null(eval_in) || length(eval_in) == 0L) return(NULL)
  ev <- suppressWarnings(as.numeric(eval_in))
  if (!length(ev)) return(NULL)
  nm <- names(eval_in)
  if (is.null(nm) || !length(nm) || all(is.na(nm) | !nzchar(nm))) {
    if (!is.null(pc_cols) && length(pc_cols)) {
      nm <- pc_cols[seq_len(min(length(pc_cols), length(ev)))]
      if (length(nm) < length(ev)) nm <- c(nm, paste0("PC", seq_len(length(ev) - length(nm)) + length(nm)))
    } else {
      nm <- paste0("PC", seq_along(ev))
    }
  }
  nm <- as.character(nm)
  nm[is.na(nm) | !nzchar(nm)] <- paste0("PC", seq_len(sum(is.na(nm) | !nzchar(nm))))
  names(ev) <- nm
  ev
}

.gcanvas_pca_propvar_from_eval <- function(eval_named) {
  ev <- .gcanvas_pca_normalize_eval(eval_named)
  if (is.null(ev) || !length(ev)) return(NULL)
  denom <- sum(as.numeric(ev), na.rm = TRUE)
  if (!is.finite(denom) || is.na(denom) || denom <= 0) return(ev * NA_real_)
  out <- as.numeric(ev) / denom
  names(out) <- names(ev)
  out
}

.gcanvas_pca_build_screeplot <- function(prop.var) {
  pv <- as_num(prop.var)
  nm <- names(prop.var)
  if (is.null(nm) || !length(nm)) nm <- paste0("PC", seq_along(pv))
  if (length(nm) != length(pv)) nm <- paste0("PC", seq_along(pv))
  dt <- data.table::data.table(
    PC = seq_along(pv),
    PC_label = factor(nm, levels = nm),
    prop.var = pv
  )
  dt[, cum.prop := cumsum(prop.var)]
  dt[!is.finite(cum.prop) | is.na(cum.prop), cum.prop := NA_real_]
  dt[, cum.prop := pmin(pmax(cum.prop, 0), 1)]
  y_left_min <- suppressWarnings(min(dt$prop.var, na.rm = TRUE))
  y_left_max <- suppressWarnings(max(dt$prop.var, na.rm = TRUE))
  if (!is.finite(y_left_min) || is.na(y_left_min)) y_left_min <- 0
  if (!is.finite(y_left_max) || is.na(y_left_max) || y_left_max <= 0) y_left_max <- 1
  y_span <- y_left_max - y_left_min
  if (!is.finite(y_span) || is.na(y_span) || y_span <= 0) y_span <- 1
  dt[, cum.y := y_left_min + cum.prop * y_span]

  ggplot2::ggplot(dt, ggplot2::aes(x = PC_label, y = prop.var, group = 1)) +
    ggplot2::geom_line(ggplot2::aes(y = cum.y, group = 1), color = "#0A9396", linewidth = 0.45) +
    ggplot2::geom_point(ggplot2::aes(y = cum.y), color = "#0A9396", size = 1.4) +
    ggplot2::geom_line(color = "grey20", linewidth = 0.5) +
    ggplot2::geom_point(color = "grey20", size = 1.5) +
    ggplot2::xlab("Principal component") +
    ggplot2::ylab("Proportion of variance") +
    ggplot2::scale_y_continuous(
      sec.axis = ggplot2::sec_axis(
        ~ (. - y_left_min) / y_span,
        name = "Cumulative proportion of variance",
        breaks = seq(0, 1, by = 0.2)
      )
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = 14),
      axis.text = ggplot2::element_text(size = 12),
      axis.title.y.right = ggplot2::element_text(color = "#0A9396", size = 14),
      axis.text.y.right = ggplot2::element_text(color = "#0A9396", size = 12),
      axis.ticks.y.right = ggplot2::element_line(color = "#0A9396", linewidth = 0.3),
      panel.border = ggplot2::element_rect(fill = NA, color = "grey20", linewidth = 0.3),
      axis.line = ggplot2::element_line(color = "grey20", linewidth = 0.3),
      axis.ticks = ggplot2::element_line(color = "grey20", linewidth = 0.3),
      panel.grid.major = ggplot2::element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_line(color = "grey92", linewidth = 0.2)
    )
}

#' Scree plot of PCA eigenvalues or PC variance
#'
#' Plots eigenvalues (or proportion of variance) against PC index.
#'
#' @param pc A PCA result list (output of [pca()]).
#' @param eval Optional eigenvalue vector if `pc` is not provided.
#' @param mode Either `"eval"` (eigenvalues) or `"pc"` (proportion variance).
#' @param pc.use Subset of PCs to plot: `"all"` or an integer vector.
#'
#' @return A `ggplot` object.
#' @export
pca.screeplot <- function(pc = NULL, eval = NULL, mode = "eval", pc.use = "all") {
  require_pkg(c("data.table", "ggplot2"))
  mode0 <- tolower(trimws(as.character(mode)[1]))
  if (!nzchar(mode0) || is.na(mode0)) mode0 <- "eval"
  if (!(mode0 %in% c("eval", "pc"))) stop("mode must be one of 'eval' or 'pc'.", call. = FALSE)

  inp <- .gcanvas_pca_extract_input(pc)
  dt <- inp$pc
  eval_use <- eval %||% inp$eval
  if (is.null(eval_use) &&
      !is.null(pc) &&
      !is.list(pc) &&
      !is.data.frame(pc) &&
      !data.table::is.data.table(pc) &&
      !is.matrix(pc)) {
    eval_use <- pc
  }
  pc_cols <- if (is.null(dt)) character() else .gcanvas_pca_pc_cols(names(dt))

  if (identical(mode0, "eval")) {
    eval_named <- .gcanvas_pca_normalize_eval(eval_use, pc_cols = pc_cols)
    if (is.null(eval_named) || !length(eval_named)) {
      stop("eval is required for mode='eval'. Provide eval or use a pca()/pca.projection() object with eigenvalue/prop.var.", call. = FALSE)
    }
    pc_idx <- .gcanvas_pca_parse_pc_use(pc.use, names(eval_named))
    pc_names <- paste0("PC", unique(pc_idx))
    pc_names <- pc_names[pc_names %in% names(eval_named)]
    if (!length(pc_names)) stop("No requested PCs found in eval.", call. = FALSE)
    prop_all <- .gcanvas_pca_propvar_from_eval(eval_named)
    total <- sum(prop_all, na.rm = TRUE)
    if (!is.finite(total) || is.na(total) || total <= 0) stop("eval must sum to a positive value.", call. = FALSE)
    prop.var <- prop_all[pc_names]
    names(prop.var) <- pc_names
    return(.gcanvas_pca_build_screeplot(prop.var))
  }

  if (is.null(dt) || !nrow(dt)) stop("pc data is required for mode='pc'.", call. = FALSE)
  if (!length(pc_cols)) stop("No PC columns found in pc.", call. = FALSE)
  for (nm in pc_cols) dt[, (nm) := suppressWarnings(as.numeric(get(nm)))]
  pc_idx <- .gcanvas_pca_parse_pc_use(pc.use, pc_cols)
  pc_names <- paste0("PC", unique(pc_idx))
  pc_names <- pc_names[pc_names %in% pc_cols]
  if (!length(pc_names)) stop("No requested PCs found in pc.", call. = FALSE)
  var_all <- vapply(pc_cols, function(nm) stats::var(dt[[nm]], na.rm = TRUE), numeric(1))
  names(var_all) <- pc_cols
  var_all[!is.finite(var_all) | is.na(var_all)] <- 0
  total <- sum(var_all, na.rm = TRUE)
  if (!is.finite(total) || is.na(total) || total <= 0) stop("PC variances must sum to a positive value.", call. = FALSE)
  prop.var <- var_all[pc_names] / total
  names(prop.var) <- pc_names
  .gcanvas_pca_build_screeplot(prop.var)
}

#' Principal component analysis of PLINK genotype data
#'
#' Runs PCA via PLINK2 (`--pca`) on a bfile or pfile, returning a list with
#' the PC scores, eigenvalues, and the resolved PLINK invocation. Provides
#' threading, output-directory, and QC-pass-through controls.
#'
#' @param bfile,pfile PLINK1 bfile prefix or PLINK2 pfile prefix.
#' @param plink Path to the `plink2` binary (default `"plink2"`).
#' @param plink.out Logical or path. Persist the PLINK log/outputs.
#' @param show.plinklog Logical. Stream PLINK's log to the R console.
#' @param threads Integer thread count passed to PLINK.
#' @param chrom Chromosomes to include (default autosomes `1:22`).
#' @param keep Optional sample-id file for `--keep`.
#' @param maf,geno Variant QC thresholds.
#' @param indep.window,indep.step,indep.r2 LD-pruning controls
#'   (`--indep-pairwise`).
#' @param pca.mode One of `"approx"` (default) or `"exact"`.
#' @param pca.count Number of PCs to extract.
#' @param plink.version PLINK version selector (`"auto"`, `"plink"`, `"plink2"`).
#' @param silent Logical. Suppress progress notes.
#'
#' @return A list with `pc` (sample scores), `eval` (eigenvalues), `meta`
#'   (run metadata), and other auxiliary slots.
#' @export
pca <- function(bfile = NULL,
                pfile = NULL,
                plink = "plink2",
                plink.out = FALSE,
                show.plinklog = TRUE,
                threads = 4L,
                chrom = 1:22,
                keep = NULL,
                maf = 0.05,
                geno = 0.02,
                indep.window = 50L,
                indep.step = 5L,
                indep.r2 = 0.1,
                pca.mode = "approx",
                pca.count = 20L,
                plink.version = "auto",
                silent = FALSE) {
  require_pkg("data.table")
  silent <- isTRUE(silent)
  show.plinklog <- isTRUE(show.plinklog)
  if (isTRUE(silent)) show.plinklog <- FALSE
  thr <- .gcanvas_resolve_threads(threads)
  pv <- tolower(trimws(as.character(plink.version)[1]))
  if (!nzchar(pv) || is.na(pv)) pv <- "auto"

  .pca_strip_quotes <- function(s) {
    s <- trimws(as.character(s))
    if (!nzchar(s)) return(s)
    gsub("^[\"']|[\"']$", "", s)
  }
  .pca_alias_target_from_rc <- function(cmd, rc_files = c("~/.zshrc", "~/.zprofile", "~/.bashrc", "~/.bash_profile")) {
    cmd <- as.character(cmd)[1]
    if (is.na(cmd) || !nzchar(cmd)) return(NA_character_)
    pat <- sprintf("^\\s*alias\\s+%s\\s*=\\s*(.+)\\s*$", gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", cmd))
    for (rf in rc_files) {
      f <- path.expand(rf)
      if (!file.exists(f)) next
      ln <- tryCatch(readLines(f, warn = FALSE), error = function(e) character(0))
      if (!length(ln)) next
      hit <- grep(pat, ln, value = TRUE, perl = TRUE)
      if (!length(hit)) next
      rhs <- sub(pat, "\\1", hit[length(hit)], perl = TRUE)
      rhs <- sub("\\s+#.*$", "", rhs)
      rhs <- .pca_strip_quotes(rhs)
      if (!nzchar(rhs)) next
      tok <- strsplit(rhs, "\\s+")[[1]]
      if (!length(tok)) next
      bin <- .pca_strip_quotes(tok[1])
      if (!nzchar(bin)) next
      w <- Sys.which(bin)
      if (nzchar(w)) return(w)
      if (file.exists(path.expand(bin))) return(abs_path(path.expand(bin)))
    }
    NA_character_
  }
  .pca_resolve_exec <- function(cmd) {
    cmd <- as.character(cmd)[1]
    if (is.na(cmd) || !nzchar(cmd)) return(NA_character_)
    if (file.exists(cmd)) return(abs_path(cmd))
    w <- Sys.which(cmd)
    if (nzchar(w)) return(w)
    cand_dir <- c(path.expand("~/bin"), "/opt/homebrew/bin", "/usr/local/bin", "/usr/bin", "/bin")
    cand <- file.path(cand_dir, cmd)
    hit <- cand[file.exists(cand)]
    if (length(hit)) return(abs_path(hit[1]))
    a <- .pca_alias_target_from_rc(cmd)
    if (!is.na(a) && nzchar(a)) return(a)
    shs <- unique(c(Sys.getenv("SHELL", unset = "/bin/sh"), "/bin/zsh", "/bin/bash", "/bin/sh"))
    shs <- shs[file.exists(shs)]
    for (sh in shs) {
      out <- tryCatch(
        suppressWarnings(system2(sh, c("-lc", sprintf("command -v %s 2>/dev/null", shQuote(cmd))), stdout = TRUE, stderr = TRUE)),
        error = function(e) character(0)
      )
      out <- trimws(as.character(out))
      out <- out[nzchar(out)]
      if (!length(out)) next
      last <- out[length(out)]
      if (file.exists(last)) return(abs_path(last))
      p3 <- Sys.which(last)
      if (nzchar(p3)) return(p3)
    }
    NA_character_
  }
  .pca_run_plink <- function(exec, args) {
    cmd_line <- paste(c(shQuote(exec), vapply(args, shQuote, character(1))), collapse = " ")
    args_run <- gsub("\\$", "\\\\\\$", as.character(args))
    if (isTRUE(show.plinklog)) {
      st <- tryCatch(
        suppressWarnings(system2(command = exec, args = args_run, stdout = "", stderr = "")),
        error = function(e) {
          stop(sprintf("Failed to run PLINK command: %s\nCommand: %s", conditionMessage(e), cmd_line), call. = FALSE)
        }
      )
      status <- as_int(st)
      if (!identical(as_int(status), 0L)) {
        stop("PLINK command failed.\nCommand: ", cmd_line, call. = FALSE)
      }
    } else {
      rc <- tryCatch(
        suppressWarnings(system2(command = exec, args = args_run, stdout = TRUE, stderr = TRUE)),
        error = function(e) {
          stop(sprintf("Failed to run PLINK command: %s\nCommand: %s", conditionMessage(e), cmd_line), call. = FALSE)
        }
      )
      status <- attr(rc, "status") %||% 0L
      if (!identical(as_int(status), 0L)) {
        msg <- if (length(rc)) paste(rc, collapse = "\n") else "PLINK command failed."
        stop("PLINK command failed.\n", msg, "\nCommand: ", cmd_line, call. = FALSE)
      }
    }
    invisible(TRUE)
  }

  if (!is.null(bfile) && !is.null(pfile)) stop("Use either bfile or pfile, not both.", call. = FALSE)
  if (is.null(bfile) && is.null(pfile)) stop("Provide bfile or pfile.", call. = FALSE)

  cache_dir0 <- .gcanvas_default_cache_dir(scope = "pca", anchor = if (!is.null(bfile)) bfile else pfile)
  dir.create(cache_dir0, recursive = TRUE, showWarnings = FALSE)
  .gcanvas_register_cache_dir(cache_dir0)
  base_pref <- file.path(cache_dir0, sprintf("pca_%s_%s", format(Sys.time(), "%Y%m%d_%H%M%S"), as.integer(stats::runif(1, 1, 1e9))))
  keep_plink <- FALSE
  if (!is.null(plink.out) && length(plink.out) > 0L) {
    if (is.logical(plink.out) && length(plink.out) == 1L) {
      keep_plink <- isTRUE(plink.out)
    } else {
      out_chr <- as.character(plink.out)[1]
      if (!is.na(out_chr) && nzchar(out_chr) && !(tolower(trimws(out_chr)) %in% c("false", "null", "na"))) {
        out_chr <- abs_path(path.expand(out_chr))
        if (dir.exists(out_chr) || grepl("[/\\\\]$", as.character(plink.out)[1])) {
          dir.create(out_chr, recursive = TRUE, showWarnings = FALSE)
          .gcanvas_register_cache_dir(out_chr)
          base_pref <- file.path(out_chr, sprintf("pca_%s_%s", format(Sys.time(), "%Y%m%d_%H%M%S"), as.integer(stats::runif(1, 1, 1e9))))
        } else {
          dir.create(dirname(out_chr), recursive = TRUE, showWarnings = FALSE)
          .gcanvas_register_cache_dir(dirname(out_chr))
          base_pref <- out_chr
        }
        keep_plink <- TRUE
      }
    }
  }
  prune_pref <- paste0(base_pref, ".prune")
  pca_pref <- paste0(base_pref, ".pca")
  keep_file <- NA_character_

  .pca_pick_col <- function(nm, keys) {
    cn <- toupper(gsub("[^A-Z0-9#]+", "", as.character(nm)))
    for (k in keys) {
      kk <- toupper(gsub("[^A-Z0-9#]+", "", as.character(k)))
      hit <- which(cn == kk)
      if (length(hit)) return(nm[hit[1]])
    }
    NA_character_
  }
  .pca_prepare_keep <- function(x) {
    if (is.null(x) || length(x) == 0L) return(NULL)
    if (data.table::is.data.table(x) || is.data.frame(x)) {
      dt <- if (data.table::is.data.table(x)) data.table::copy(x) else data.table::as.data.table(x)
      if (!nrow(dt)) return(NULL)
      c_fid <- .pca_pick_col(names(dt), c("FID", "#FID", "family.ID", "familyID", "FAMILY_ID"))
      c_iid <- .pca_pick_col(names(dt), c("IID", "#IID", "ID", "sample.ID", "sampleID", "INDV", "INDIVIDUAL"))
      if (is.na(c_iid)) {
        if (ncol(dt) >= 2L) {
          c_fid <- names(dt)[1]
          c_iid <- names(dt)[2]
        } else {
          c_iid <- names(dt)[1]
        }
      }
      if (is.na(c_fid)) {
        out <- dt[, .(IID = as.character(get(c_iid)))]
        out <- out[!is.na(IID) & nzchar(IID)]
        out <- unique(out, by = "IID")
        return(out)
      }
      out <- dt[, .(FID = as.character(get(c_fid)), IID = as.character(get(c_iid)))]
      out <- out[!is.na(IID) & nzchar(IID)]
      out <- out[is.na(FID) | !nzchar(FID), FID := IID]
      out <- unique(out, by = c("FID", "IID"))
      out
    } else {
      iid <- as.character(x)
      iid <- iid[!is.na(iid) & nzchar(iid)]
      if (!length(iid)) return(NULL)
      iid <- unique(iid)
      data.table::data.table(IID = iid)
    }
  }
  .pca_load_sample_map <- function() {
    if (!is.null(bfile)) {
      fam <- paste0(as.character(bfile)[1], ".fam")
      if (!file.exists(fam)) return(NULL)
      dtf <- data.table::fread(fam, data.table = TRUE, header = FALSE, showProgress = FALSE)
      if (ncol(dtf) < 2L) return(NULL)
      out <- data.table::data.table(FID = as.character(dtf[[1]]), IID = as.character(dtf[[2]]))
      out <- out[!is.na(IID) & nzchar(IID)]
      return(unique(out, by = c("FID", "IID")))
    }
    if (!is.null(pfile)) {
      psam <- paste0(as.character(pfile)[1], ".psam")
      if (!file.exists(psam)) return(NULL)
      dts <- data.table::fread(psam, data.table = TRUE, check.names = FALSE, showProgress = FALSE)
      if (!nrow(dts)) return(NULL)
      c_iid <- .pca_pick_col(names(dts), c("#IID", "IID", "ID", "sample.ID", "sampleID"))
      c_fid <- .pca_pick_col(names(dts), c("#FID", "FID", "family.ID", "familyID"))
      if (is.na(c_iid)) return(NULL)
      out <- if (is.na(c_fid)) {
        data.table::data.table(IID = as.character(dts[[c_iid]]))
      } else {
        data.table::data.table(FID = as.character(dts[[c_fid]]), IID = as.character(dts[[c_iid]]))
      }
      out <- out[!is.na(IID) & nzchar(IID)]
      if (!("FID" %in% names(out))) out[, FID := IID]
      out[is.na(FID) | !nzchar(FID), FID := IID]
      return(unique(out[, .(FID, IID)], by = c("FID", "IID")))
    }
    NULL
  }

  args_base <- character()
  if (!is.null(bfile)) args_base <- c(args_base, "--bfile", as.character(bfile)[1])
  if (!is.null(pfile)) {
    pfx <- as.character(pfile)[1]
    use_vzs <- file.exists(paste0(pfx, ".pvar.zst")) && !file.exists(paste0(pfx, ".pvar"))
    args_base <- c(args_base, "--pfile", pfx, if (isTRUE(use_vzs)) "vzs")
  }
  keep_dt <- .pca_prepare_keep(keep)
  keep_n <- 0L
  if (!is.null(keep_dt) && nrow(keep_dt)) {
    sample_map <- .pca_load_sample_map()
    if (!is.null(sample_map) && nrow(sample_map)) {
      if (!("FID" %in% names(keep_dt))) {
        keep_dt <- merge(
          keep_dt[, .(IID = as.character(IID))],
          sample_map[, .(FID, IID)],
          by = "IID",
          all.x = FALSE,
          all.y = FALSE,
          sort = FALSE
        )
      } else {
        keep_dt <- keep_dt[, .(FID = as.character(FID), IID = as.character(IID))]
        keep_dt <- keep_dt[!is.na(IID) & nzchar(IID)]
        keep_dt[is.na(FID) | !nzchar(FID), FID := IID]
        data.table::setkey(sample_map, FID, IID)
        data.table::setkey(keep_dt, FID, IID)
        keep_dt <- sample_map[keep_dt, nomatch = 0L]
      }
      keep_dt <- unique(keep_dt, by = c("FID", "IID"))
      if (!nrow(keep_dt)) {
        stop("keep has no matching samples in input data (FID/IID after matching).", call. = FALSE)
      }
    } else {
      if (!("FID" %in% names(keep_dt))) keep_dt[, FID := IID]
      keep_dt <- unique(keep_dt[, .(FID, IID)], by = c("FID", "IID"))
    }
    keep_file <- paste0(base_pref, ".keep.txt")
    data.table::fwrite(keep_dt[, .(FID, IID)], keep_file, sep = "\t", col.names = FALSE, quote = FALSE)
    keep_n <- as_int(nrow(keep_dt))
    args_base <- c(args_base, "--keep", keep_file)
  }
  .pca_normalize_chrom <- function(x) {
    if (is.null(x) || length(x) == 0L) return(character())
    xc <- as.character(x)
    xc <- xc[!is.na(xc) & nzchar(xc)]
    if (!length(xc)) return(character())
    out <- character()
    for (v in xc) {
      vv <- trimws(v)
      if (!nzchar(vv)) next
      # Keep range/string tokens as-is (e.g., "1-22", "X,Y")
      if (grepl("[-,]", vv)) {
        parts <- unlist(strsplit(vv, ",", fixed = TRUE))
        parts <- trimws(parts)
        parts <- parts[nzchar(parts)]
        out <- c(out, parts)
      } else {
        out <- c(out, vv)
      }
    }
    out <- toupper(gsub("^CHR", "", out, ignore.case = TRUE))
    out[out == "M"] <- "MT"
    out
  }
  chr_vec <- .pca_normalize_chrom(chrom)
  if (!length(chr_vec)) chr_vec <- as.character(1:22)
  chr0 <- paste(chr_vec, collapse = ",")
  maf0 <- as_num(maf)[1]; if (!is.finite(maf0) || is.na(maf0) || maf0 < 0 || maf0 > 1) maf0 <- 0.05
  geno0 <- as_num(geno)[1]; if (!is.finite(geno0) || is.na(geno0) || geno0 < 0 || geno0 > 1) geno0 <- 0.02
  indep_window0 <- as_int(indep.window)[1]; if (!is.finite(indep_window0) || is.na(indep_window0) || indep_window0 < 1L) indep_window0 <- 50L
  indep_step0 <- as_int(indep.step)[1]; if (!is.finite(indep_step0) || is.na(indep_step0) || indep_step0 < 1L) indep_step0 <- 5L
  indep_r20 <- as_num(indep.r2)[1]; if (!is.finite(indep_r20) || is.na(indep_r20) || indep_r20 <= 0 || indep_r20 >= 1) indep_r20 <- 0.1
  mode0 <- tolower(trimws(as.character(pca.mode)[1])); if (is.na(mode0) || !nzchar(mode0)) mode0 <- "approx"

  plink_cmd <- as.character(plink)[1]
  if (is.na(plink_cmd) || !nzchar(plink_cmd)) stop("PLINK command is empty.", call. = FALSE)
  plink_exec <- .pca_resolve_exec(plink_cmd)
  if (is.na(plink_exec) || !nzchar(plink_exec)) {
    stop(sprintf("PLINK executable not found for '%s'. Set `plink` to an executable path.", plink_cmd), call. = FALSE)
  }

  .gcanvas_note(
    "gcanvas::pca",
    sprintf("Start: mode=%s | plink=%s | plink.version=%s | threads=%d | keep.n=%d",
            if (!is.null(bfile)) "bfile" else "pfile",
            plink_cmd, pv, as_int(thr), as_int(keep_n)),
    silent = silent
  )
  if (isTRUE(keep_plink)) {
    .gcanvas_note("gcanvas::pca", sprintf("PLINK outputs kept: %s", base_pref), silent = silent)
  }

  args_prune <- c(
    args_base,
    "--chr", chr0,
    "--maf", format(maf0, scientific = FALSE),
    "--geno", format(geno0, scientific = FALSE),
    "--indep-pairwise", as.character(indep_window0), as.character(indep_step0), format(indep_r20, scientific = FALSE),
    "--threads", as.character(as_int(thr)),
    "--out", prune_pref
  )
  .pca_run_plink(plink_exec, args_prune)

  args_pca <- c(
    args_base,
    "--extract", paste0(prune_pref, ".prune.in"),
    "--freq",
    "--pca", mode0, "allele-wts",
    "--threads", as.character(as_int(thr)),
    "--out", pca_pref
  )
  if (!is.null(pca.count) && length(pca.count) > 0L && !is.na(pca.count[1])) {
    k <- as_int(pca.count)[1]
    if (is.finite(k) && !is.na(k) && k > 0L) {
      args_pca <- c(
        args_base,
        "--extract", paste0(prune_pref, ".prune.in"),
        "--freq",
        "--pca", mode0, as.character(k), "allele-wts",
        "--threads", as.character(as_int(thr)),
        "--out", pca_pref
      )
    }
  }
  .pca_run_plink(plink_exec, args_pca)

  pc_file <- paste0(pca_pref, ".eigenvec")
  ev_file <- paste0(pca_pref, ".eigenval")
  wt_file <- paste0(pca_pref, ".eigenvec.allele")
  if (!file.exists(pc_file) || !file.exists(ev_file) || !file.exists(wt_file)) {
    stop("PLINK PCA finished but expected .eigenvec/.eigenval/.eigenvec.allele not found.", call. = FALSE)
  }

  .pca_pick_aux_col <- function(nm, keys) {
    cn <- toupper(gsub("[^A-Z0-9#]+", "", as.character(nm)))
    for (k in keys) {
      kk <- toupper(gsub("[^A-Z0-9#]+", "", as.character(k)))
      hit <- which(cn == kk)
      if (length(hit)) return(nm[hit[1]])
    }
    NA_character_
  }
  .pca_add_a1_freq <- function(dt) {
    if (is.null(dt) || !nrow(dt)) return(dt)
    if ("A1.freq" %in% names(dt)) {
      dt[, A1.freq := suppressWarnings(as.numeric(dt[["A1.freq"]]))]
      return(dt)
    }
    c_a1 <- .pca_pick_aux_col(names(dt), c("A1", "ALLELE1", "ALT", "EFFECTALLELE"))
    c_ref <- .pca_pick_aux_col(names(dt), c("REF"))
    c_alt <- .pca_pick_aux_col(names(dt), c("ALT"))
    c_af <- .pca_pick_aux_col(names(dt), c("ALT_FREQS", "ALTFREQS", "ALTFREQ", "AF"))
    if (is.na(c_a1) || is.na(c_ref) || is.na(c_alt) || is.na(c_af)) {
      if (!("A1.freq" %in% names(dt))) dt[, A1.freq := NA_real_]
      return(dt)
    }
    a1v <- toupper(as.character(dt[[c_a1]]))
    refv <- toupper(as.character(dt[[c_ref]]))
    altv <- toupper(as.character(dt[[c_alt]]))
    afv <- suppressWarnings(as.numeric(dt[[c_af]]))
    a1f <- ifelse(!is.na(a1v) & !is.na(altv) & a1v == altv, afv,
                  ifelse(!is.na(a1v) & !is.na(refv) & a1v == refv, 1 - afv, NA_real_))
    dt[, A1.freq := a1f]
    dt
  }
  .pca_load_variant_pos <- function() {
    if (!is.null(bfile)) {
      bim <- paste0(as.character(bfile)[1], ".bim")
      if (!file.exists(bim)) return(NULL)
      dtb <- data.table::fread(bim, data.table = TRUE, header = FALSE, showProgress = FALSE)
      if (ncol(dtb) < 4L) return(NULL)
      out <- data.table::data.table(
        CHR = as.character(dtb[[1]]),
        ID = as.character(dtb[[2]]),
        POS = suppressWarnings(as.numeric(dtb[[4]]))
      )
      return(out)
    }
    if (!is.null(pfile)) {
      pfx <- as.character(pfile)[1]
      pv <- if (file.exists(paste0(pfx, ".pvar"))) paste0(pfx, ".pvar") else if (file.exists(paste0(pfx, ".pvar.zst"))) paste0(pfx, ".pvar.zst") else NA_character_
      if (is.na(pv) || !file.exists(pv)) return(NULL)
      dtp <- tryCatch(suppressWarnings(data.table::fread(pv, data.table = TRUE, showProgress = FALSE, quote = "")), error = function(e) NULL)
      if (is.null(dtp) || !nrow(dtp)) return(NULL)
      cn <- names(dtp)
      c_chr <- .pca_pick_aux_col(cn, c("#CHROM", "CHROM", "CHR"))
      c_pos <- .pca_pick_aux_col(cn, c("POS", "BP"))
      c_id <- .pca_pick_aux_col(cn, c("ID", "SNP"))
      if (is.na(c_chr) || is.na(c_pos) || is.na(c_id)) return(NULL)
      out <- data.table::data.table(
        CHR = as.character(dtp[[c_chr]]),
        ID = as.character(dtp[[c_id]]),
        POS = suppressWarnings(as.numeric(dtp[[c_pos]]))
      )
      return(out)
    }
    NULL
  }
  .pca_fread_quiet <- function(...) {
    withCallingHandlers(
      data.table::fread(...),
      warning = function(w) {
        msg <- conditionMessage(w)
        if (grepl("Previous fread\\(\\) session was not cleaned up properly", msg, fixed = FALSE)) {
          invokeRestart("muffleWarning")
        }
      }
    )
  }

  pc <- .pca_fread_quiet(pc_file, data.table = TRUE, check.names = FALSE, showProgress = FALSE)
  if ("#FID" %in% names(pc)) data.table::setnames(pc, "#FID", "FID")
  if ("#IID" %in% names(pc)) data.table::setnames(pc, "#IID", "IID")
  if ("FID" %in% names(pc)) {
    fidv <- as.character(pc$FID)
    no_fid <- all(is.na(fidv) | !nzchar(fidv) | fidv %in% c("0", "NA"))
    if (isTRUE(no_fid)) pc[, FID := NULL]
  }

  eigenvalue_dt <- .pca_fread_quiet(ev_file, data.table = TRUE, header = FALSE, showProgress = FALSE)
  if (ncol(eigenvalue_dt) == 1L) data.table::setnames(eigenvalue_dt, "V1", "EIGENVAL")
  if (!("EIGENVAL" %in% names(eigenvalue_dt))) {
    data.table::setnames(eigenvalue_dt, names(eigenvalue_dt)[1], "EIGENVAL")
  }
  eigenvalue_dt[, EIGENVAL := suppressWarnings(as.numeric(EIGENVAL))]
  eig_vec <- as.numeric(eigenvalue_dt$EIGENVAL)
  names(eig_vec) <- paste0("PC", seq_along(eig_vec))
  prop.var <- .gcanvas_pca_propvar_from_eval(eig_vec)
  prop.var <- as.numeric(prop.var)
  names(prop.var) <- paste0("PC", seq_along(prop.var))
  afreq <- NULL
  afreq_file <- paste0(pca_pref, ".afreq")
  if (file.exists(afreq_file)) {
    afreq <- .pca_fread_quiet(afreq_file, data.table = TRUE, check.names = FALSE, showProgress = FALSE)
    if ("#CHROM" %in% names(afreq)) data.table::setnames(afreq, "#CHROM", "CHR")
    c_id_af <- .pca_pick_aux_col(names(afreq), c("ID", "SNP", "VAR"))
    if (!is.na(c_id_af) && c_id_af != "ID") data.table::setnames(afreq, c_id_af, "ID")
    if ("ID" %in% names(afreq)) {
      afreq <- afreq[!is.na(ID) & nzchar(ID)]
      if (!("POS" %in% names(afreq))) {
        pos_ref0 <- .pca_load_variant_pos()
        if (!is.null(pos_ref0) && nrow(pos_ref0)) {
          pos_ref0 <- unique(pos_ref0, by = "ID")
          afreq <- merge(afreq, pos_ref0[, .(ID, POS_ref = POS)], by = "ID", all.x = TRUE, sort = FALSE)
          if (!("POS" %in% names(afreq))) afreq[, POS := POS_ref]
          afreq[, POS_ref := NULL]
        }
      }
      if ("POS" %in% names(afreq)) afreq[, POS := suppressWarnings(as.numeric(POS))]
      afreq <- unique(afreq, by = "ID")
    } else {
      afreq <- NULL
    }
  }

  weight <- .pca_fread_quiet(wt_file, data.table = TRUE, check.names = FALSE, showProgress = FALSE)
  if ("#CHROM" %in% names(weight)) data.table::setnames(weight, "#CHROM", "CHR")
  c_id_w <- .pca_pick_aux_col(names(weight), c("ID", "SNP"))
  if (is.na(c_id_w)) c_id_w <- names(weight)[1]
  if (c_id_w != "ID") data.table::setnames(weight, c_id_w, "ID")
  pos_ref <- .pca_load_variant_pos()
  if (!is.null(pos_ref) && nrow(pos_ref)) {
    pos_ref <- unique(pos_ref, by = "ID")
    pos_ref <- pos_ref[, .(ID, CHR_ref = CHR, POS_ref = POS)]
    weight <- merge(weight, pos_ref, by = "ID", all.x = TRUE, sort = FALSE)
  } else {
    if (!("CHR_ref" %in% names(weight))) weight[, CHR_ref := NA_character_]
    if (!("POS_ref" %in% names(weight))) weight[, POS_ref := NA_real_]
  }
  if ("CHR.x" %in% names(weight) || "CHR.y" %in% names(weight)) {
    chrx <- if ("CHR.x" %in% names(weight)) as.character(weight[["CHR.x"]]) else rep(NA_character_, nrow(weight))
    chry <- if ("CHR.y" %in% names(weight)) as.character(weight[["CHR.y"]]) else rep(NA_character_, nrow(weight))
    weight[, CHR := data.table::fifelse(!is.na(chrx) & nzchar(chrx), chrx, chry)]
    rm_chrxy <- intersect(c("CHR.x", "CHR.y"), names(weight))
    if (length(rm_chrxy)) weight[, (rm_chrxy) := NULL]
  }
  if (!("CHR" %in% names(weight))) weight[, CHR := NA_character_]
  if ("CHR_ref" %in% names(weight)) {
    wchr <- as.character(weight$CHR)
    cref <- as.character(weight$CHR_ref)
    fill_chr <- (is.na(wchr) | !nzchar(wchr)) & !is.na(cref) & nzchar(cref)
    if (any(fill_chr)) weight[fill_chr, CHR := CHR_ref]
    weight[, CHR_ref := NULL]
  }
  if (!("POS" %in% names(weight))) weight[, POS := NA_real_]
  if ("POS_ref" %in% names(weight)) {
    p0 <- suppressWarnings(as.numeric(weight$POS))
    pr <- suppressWarnings(as.numeric(weight$POS_ref))
    fill_pos <- !is.finite(p0) & is.finite(pr)
    if (any(fill_pos)) weight[fill_pos, POS := POS_ref]
    weight[, POS_ref := NULL]
  }
  weight[, POS := suppressWarnings(as.numeric(POS))]
  if (!is.null(afreq) && nrow(afreq)) {
    keep_af_cols <- intersect(c("ID", "ALT_FREQS", "OBS_CT"), names(afreq))
    if (length(keep_af_cols) >= 2L) {
      weight <- merge(weight, afreq[, ..keep_af_cols], by = "ID", all.x = TRUE, sort = FALSE)
    }
  }
  weight <- .pca_add_a1_freq(weight)
  if ("ALT_FREQS" %in% names(weight)) weight[, ALT_FREQS := NULL]
  lead_cols <- c("CHR", "POS", "ID", "REF", "ALT", "PROVISIONAL_REF?", "A1", "A1.freq", "OBS_CT")
  rest_cols <- setdiff(names(weight), lead_cols)
  data.table::setcolorder(weight, c(lead_cols, rest_cols))
  attr(weight, "gcanvas_pca_eigenvalue") <- eig_vec
  attr(weight, "gcanvas_pca_prop.var") <- prop.var
  attr(weight, "gcanvas_pca_afreq") <- afreq

  scree <- .gcanvas_pca_build_screeplot(prop.var)

  if (!isTRUE(keep_plink)) {
    rm_files <- c(
      paste0(prune_pref, ".prune.in"),
      paste0(prune_pref, ".prune.out"),
      paste0(prune_pref, ".log"),
      paste0(prune_pref, ".nosex"),
      keep_file,
      pc_file, ev_file, wt_file,
      afreq_file,
      paste0(pca_pref, ".acf"),
      paste0(pca_pref, ".log"),
      paste0(pca_pref, ".nosex")
    )
    rm_files <- rm_files[file.exists(rm_files)]
    if (length(rm_files)) unlink(rm_files, force = TRUE)
  }
  if (isTRUE(show.plinklog) && !isTRUE(keep_plink)) {
    .gcanvas_note(
      "gcanvas::pca",
      "PLINK intermediate files were removed (plink.out=FALSE). Set plink.out=TRUE or a path to keep them.",
      silent = silent
    )
  }

  .gcanvas_note(
    "gcanvas::pca",
    sprintf("Done: n_samples=%d | n_pcs=%d | n_weights=%d",
            as_int(nrow(pc)),
            as_int(max(0L, ncol(pc) - as_int(sum(toupper(names(pc)) %in% c("FID", "IID", "#FID", "#IID", "ID"))))),
            as_int(nrow(weight))),
    silent = silent
  )

  out <- list(
    pc = pc,
    eigenvalue = eig_vec,
    weight = weight,
    prop.var = prop.var,
    screeplot = scree
  )
  class(out) <- c("gcanvas_pca", "list")
  out
}

#' Project new samples onto a reference PCA
#'
#' Uses the variant loadings from a reference [pca()] run to compute PC
#' scores for samples in a separate genotype file (target), enabling apples-
#' to-apples comparison to the reference population's PC space.
#'
#' @param ref A reference PCA result (output of [pca()]).
#' @param target A PLINK bfile/pfile prefix or precomputed target genotype matrix.
#' @param plink Path to the `plink2` binary.
#' @param plink.format `"pfile"` or `"bfile"`.
#' @param plink.out Logical or path. Persist PLINK outputs.
#' @param id.aligned Logical. Assume target sample ids already align with `ref`.
#' @param force Logical. Force re-projection even when cached intermediates exist.
#' @param show.plinklog Logical. Stream PLINK's log to the R console.
#' @param threads Integer thread count.
#' @param chrom Chromosomes to include (default autosomes).
#' @param maf,geno Variant QC thresholds.
#' @param indep.window,indep.step,indep.r2 LD-pruning controls.
#' @param pca.mode One of `"approx"` (default) or `"exact"`.
#' @param pca.count Number of PCs to project.
#' @param plink.version PLINK version selector.
#' @param silent Logical. Suppress progress notes.
#'
#' @return A list mirroring the [pca()] output with `pc` containing the
#'   projected scores for `target`.
#' @export
pca.projection <- function(ref,
                           target,
                           plink = "plink2",
                           plink.format = "pfile",
                           plink.out = FALSE,
                           id.aligned = FALSE,
                           force = FALSE,
                           show.plinklog = TRUE,
                           threads = 4L,
                           chrom = 1:22,
                           maf = 0.05,
                           geno = 0.02,
                           indep.window = 50L,
                           indep.step = 5L,
                           indep.r2 = 0.1,
                           pca.mode = "approx",
                           pca.count = 20L,
                           plink.version = "auto",
                           silent = FALSE) {
  require_pkg(c("data.table", "ggplot2"))
  silent <- isTRUE(silent)
  show.plinklog <- isTRUE(show.plinklog)
  if (isTRUE(silent)) show.plinklog <- FALSE
  id.aligned <- isTRUE(id.aligned)
  force <- isTRUE(force)
  thr <- .gcanvas_resolve_threads(threads)

  .pca_projection_canon <- function(nm) toupper(gsub("[^A-Z0-9#]+", "", as.character(nm)))
  .pca_projection_pick_col <- function(nm, keys, default = NA_character_) {
    cn <- .pca_projection_canon(nm)
    for (k in keys) {
      hit <- which(cn == .pca_projection_canon(k))
      if (length(hit)) return(nm[hit[1]])
    }
    default
  }
  .pca_projection_strip_quotes <- function(s) {
    s <- trimws(as.character(s))
    if (!nzchar(s)) return(s)
    gsub("^[\"']|[\"']$", "", s)
  }
  .pca_projection_alias_target_from_rc <- function(cmd, rc_files = c("~/.zshrc", "~/.zprofile", "~/.bashrc", "~/.bash_profile")) {
    cmd <- as.character(cmd)[1]
    if (is.na(cmd) || !nzchar(cmd)) return(NA_character_)
    pat <- sprintf("^\\s*alias\\s+%s\\s*=\\s*(.+)\\s*$", gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", cmd))
    for (rf in rc_files) {
      f <- path.expand(rf)
      if (!file.exists(f)) next
      ln <- tryCatch(readLines(f, warn = FALSE), error = function(e) character(0))
      if (!length(ln)) next
      hit <- grep(pat, ln, value = TRUE, perl = TRUE)
      if (!length(hit)) next
      rhs <- sub(pat, "\\1", hit[length(hit)], perl = TRUE)
      rhs <- sub("\\s+#.*$", "", rhs)
      rhs <- .pca_projection_strip_quotes(rhs)
      if (!nzchar(rhs)) next
      tok <- strsplit(rhs, "\\s+")[[1]]
      if (!length(tok)) next
      bin <- .pca_projection_strip_quotes(tok[1])
      if (!nzchar(bin)) next
      w <- Sys.which(bin)
      if (nzchar(w)) return(w)
      if (file.exists(path.expand(bin))) return(abs_path(path.expand(bin)))
    }
    NA_character_
  }
  .pca_projection_resolve_exec <- function(cmd) {
    cmd <- as.character(cmd)[1]
    if (is.na(cmd) || !nzchar(cmd)) return(NA_character_)
    if (file.exists(cmd)) return(abs_path(cmd))
    w <- Sys.which(cmd)
    if (nzchar(w)) return(w)
    cand_dir <- c(path.expand("~/bin"), "/opt/homebrew/bin", "/usr/local/bin", "/usr/bin", "/bin")
    cand <- file.path(cand_dir, cmd)
    hit <- cand[file.exists(cand)]
    if (length(hit)) return(abs_path(hit[1]))
    a <- .pca_projection_alias_target_from_rc(cmd)
    if (!is.na(a) && nzchar(a)) return(a)
    shs <- unique(c(Sys.getenv("SHELL", unset = "/bin/sh"), "/bin/zsh", "/bin/bash", "/bin/sh"))
    shs <- shs[file.exists(shs)]
    for (sh in shs) {
      out <- tryCatch(
        suppressWarnings(system2(sh, c("-lc", sprintf("command -v %s 2>/dev/null", shQuote(cmd))), stdout = TRUE, stderr = TRUE)),
        error = function(e) character(0)
      )
      out <- trimws(as.character(out))
      out <- out[nzchar(out)]
      if (!length(out)) next
      last <- out[length(out)]
      if (file.exists(last)) return(abs_path(last))
      p2 <- Sys.which(last)
      if (nzchar(p2)) return(p2)
    }
    NA_character_
  }
  .pca_projection_run_plink <- function(exec, args) {
    cmd_line <- paste(c(shQuote(exec), vapply(args, shQuote, character(1))), collapse = " ")
    args_run <- gsub("\\$", "\\\\\\$", as.character(args))
    if (isTRUE(show.plinklog)) {
      st <- tryCatch(
        suppressWarnings(system2(command = exec, args = args_run, stdout = "", stderr = "")),
        error = function(e) {
          stop(sprintf("Failed to run PLINK command: %s\nCommand: %s", conditionMessage(e), cmd_line), call. = FALSE)
        }
      )
      status <- as_int(st)
      if (!identical(as_int(status), 0L)) {
        stop("PLINK command failed.\nCommand: ", cmd_line, call. = FALSE)
      }
    } else {
      rc <- tryCatch(
        suppressWarnings(system2(command = exec, args = args_run, stdout = TRUE, stderr = TRUE)),
        error = function(e) {
          stop(sprintf("Failed to run PLINK command: %s\nCommand: %s", conditionMessage(e), cmd_line), call. = FALSE)
        }
      )
      status <- attr(rc, "status") %||% 0L
      if (!identical(as_int(status), 0L)) {
        msg <- if (length(rc)) paste(rc, collapse = "\n") else "PLINK command failed."
        stop("PLINK command failed.\n", msg, "\nCommand: ", cmd_line, call. = FALSE)
      }
    }
    invisible(TRUE)
  }
  .pca_projection_norm_prefix <- function(x) {
    s <- as.character(x)[1]
    if (is.na(s) || !nzchar(s)) return(NA_character_)
    s <- trimws(s)
    s <- sub("\\.(bed|bim|fam|pgen|psam|pvar|pvar\\.zst)$", "", s, ignore.case = TRUE)
    abs_path(s)
  }
  .pca_projection_dataset_exists <- function(prefix) {
    is_b <- all(file.exists(paste0(prefix, c(".bed", ".bim", ".fam"))))
    is_p <- file.exists(paste0(prefix, ".pgen")) &&
      file.exists(paste0(prefix, ".psam")) &&
      (file.exists(paste0(prefix, ".pvar")) || file.exists(paste0(prefix, ".pvar.zst")))
    list(bfile = is_b, pfile = is_p)
  }
  .pca_projection_detect_dataset <- function(x, tag, prefer = "pfile") {
    pfx <- .pca_projection_norm_prefix(x)
    if (is.na(pfx) || !nzchar(pfx)) stop(sprintf("%s prefix is empty.", tag), call. = FALSE)
    ex <- .pca_projection_dataset_exists(pfx)
    if (!isTRUE(ex$bfile) && !isTRUE(ex$pfile)) {
      stop(sprintf("%s not found as bfile/pfile prefix: %s", tag, pfx), call. = FALSE)
    }
    fmt <- NA_character_
    if (isTRUE(ex$bfile) && isTRUE(ex$pfile)) {
      if (prefer %in% c("bfile", "pfile")) fmt <- prefer else fmt <- "pfile"
    } else if (isTRUE(ex$bfile)) {
      fmt <- "bfile"
    } else {
      fmt <- "pfile"
    }
    use_vzs <- isTRUE(fmt == "pfile") && file.exists(paste0(pfx, ".pvar.zst")) && !file.exists(paste0(pfx, ".pvar"))
    list(prefix = pfx, format = fmt, use_vzs = use_vzs, tag = tag)
  }
  .pca_projection_dataset_args <- function(ds) {
    if (identical(ds$format, "bfile")) {
      return(c("--bfile", ds$prefix))
    }
    c("--pfile", ds$prefix, if (isTRUE(ds$use_vzs)) "vzs")
  }
  .pca_projection_normalize_chrom <- function(x) {
    if (is.null(x) || length(x) == 0L) return(character())
    xc <- as.character(x)
    xc <- xc[!is.na(xc) & nzchar(xc)]
    if (!length(xc)) return(character())
    out <- character()
    for (v in xc) {
      vv <- trimws(v)
      if (!nzchar(vv)) next
      if (grepl("[-,]", vv)) {
        parts <- unlist(strsplit(vv, ",", fixed = TRUE))
        parts <- trimws(parts)
        parts <- parts[nzchar(parts)]
        out <- c(out, parts)
      } else {
        out <- c(out, vv)
      }
    }
    out <- toupper(gsub("^CHR", "", out, ignore.case = TRUE))
    out[out == "M"] <- "MT"
    out
  }
  .pca_projection_read_variant_ids <- function(ds) {
    if (identical(ds$format, "bfile")) {
      bim <- paste0(ds$prefix, ".bim")
      if (!file.exists(bim)) stop("Missing .bim: ", bim, call. = FALSE)
      dt <- data.table::fread(bim, data.table = TRUE, header = FALSE, showProgress = FALSE)
      if (ncol(dt) < 2L) stop("Invalid .bim (requires >=2 columns): ", bim, call. = FALSE)
      ids <- as.character(dt[[2]])
    } else {
      pvar <- if (file.exists(paste0(ds$prefix, ".pvar"))) paste0(ds$prefix, ".pvar") else paste0(ds$prefix, ".pvar.zst")
      if (!file.exists(pvar)) stop("Missing .pvar/.pvar.zst for prefix: ", ds$prefix, call. = FALSE)
      dt <- suppressWarnings(data.table::fread(pvar, data.table = TRUE, check.names = FALSE, showProgress = FALSE, quote = ""))
      c_id <- .pca_projection_pick_col(names(dt), c("ID", "SNP", "VAR"))
      if (is.na(c_id)) stop("Could not find ID column in pvar: ", pvar, call. = FALSE)
      ids <- as.character(dt[[c_id]])
    }
    ids <- ids[!is.na(ids) & nzchar(ids)]
    unique(ids)
  }
  .pca_projection_compress_colnums <- function(idx) {
    idx <- sort(unique(as_int(idx)))
    idx <- idx[is.finite(idx) & !is.na(idx) & idx > 0L]
    if (!length(idx)) return(NA_character_)
    if (length(idx) == 1L) return(as.character(idx))
    d <- c(Inf, diff(idx))
    grp <- cumsum(d != 1L)
    sp <- split(idx, grp)
    parts <- vapply(sp, function(v) {
      if (length(v) == 1L) as.character(v[1]) else sprintf("%d-%d", v[1], v[length(v)])
    }, character(1))
    paste(parts, collapse = ",")
  }
  .pca_projection_parse_score_spec <- function(weight_file, pca_count = NULL) {
    hdr <- data.table::fread(weight_file, nrows = 0L, data.table = TRUE, check.names = FALSE, showProgress = FALSE)
    nm <- names(hdr)
    if (!length(nm)) stop("Failed to read header from weight file: ", weight_file, call. = FALSE)
    c_id <- .pca_projection_pick_col(nm, c("ID", "SNP", "VAR", "VARIANT"))
    c_a1 <- .pca_projection_pick_col(nm, c("A1", "ALLELE1", "ALT", "EFFECTALLELE"))
    if (is.na(c_id) || is.na(c_a1)) {
      stop("Weight file must contain ID and A1-like columns for --score.", call. = FALSE)
    }
    pc_nm <- nm[grepl("^PC[0-9]+$", toupper(nm))]
    if (!length(pc_nm)) pc_nm <- nm[grepl("^PC[0-9]+", toupper(nm))]
    if (!length(pc_nm)) stop("No PC weight columns found in weight file.", call. = FALSE)
    ord <- order(as_int(gsub("[^0-9]", "", toupper(pc_nm))), na.last = TRUE)
    pc_nm <- pc_nm[ord]
    if (!is.null(pca_count) && length(pca_count) > 0L && !is.na(pca_count[1])) {
      k <- as_int(pca_count)[1]
      if (is.finite(k) && !is.na(k) && k > 0L) {
        pc_nm <- pc_nm[seq_len(min(length(pc_nm), as_int(k)))]
      }
    }
    idx <- match(pc_nm, nm)
    spec <- .pca_projection_compress_colnums(idx)
    list(
      id_col = as_int(match(c_id, nm)),
      a1_col = as_int(match(c_a1, nm)),
      pc_names = as.character(pc_nm),
      pc_indices = as_int(idx),
      score_spec = as.character(spec)
    )
  }
  .pca_projection_read_sscore <- function(prefix, pc_names, eig_vals, is_ref) {
    f <- paste0(prefix, ".sscore")
    if (!file.exists(f)) stop("Missing projected score file: ", f, call. = FALSE)
    dt <- data.table::fread(f, data.table = TRUE, check.names = FALSE, showProgress = FALSE)
    if ("#FID" %in% names(dt)) data.table::setnames(dt, "#FID", "FID")
    if ("#IID" %in% names(dt)) data.table::setnames(dt, "#IID", "IID")
    c_iid <- .pca_projection_pick_col(names(dt), c("IID", "ID", "#IID"))
    if (is.na(c_iid)) stop("Could not find IID column in .sscore: ", f, call. = FALSE)
    if (c_iid != "IID") data.table::setnames(dt, c_iid, "IID")
    c_fid <- .pca_projection_pick_col(names(dt), c("FID", "#FID"))
    out <- data.table::data.table(
      IID = as.character(dt$IID),
      reference = as.logical(is_ref)
    )
    if (!is.na(c_fid)) out[, FID := as.character(dt[[c_fid]])]

    nm <- names(dt)
    pc_avg <- paste0(pc_names, "_AVG")
    if (all(pc_avg %in% nm)) {
      src <- pc_avg
    } else if (all(pc_names %in% nm)) {
      src <- pc_names
    } else {
      cand <- nm[grepl("^SCORE[0-9]+_AVG$", toupper(nm))]
      if (!length(cand)) cand <- nm[grepl("_AVG$", toupper(nm))]
      if (length(cand) < length(pc_names)) {
        stop("Could not map projected PC columns from .sscore: ", f, call. = FALSE)
      }
      ord <- order(as_int(gsub("[^0-9]", "", toupper(cand))), na.last = TRUE)
      src <- cand[ord][seq_len(length(pc_names))]
    }

    eig_use <- as_num(eig_vals)
    eig_use <- eig_use[seq_len(min(length(eig_use), length(pc_names)))]
    if (length(eig_use) < length(pc_names)) {
      eig_use <- c(eig_use, rep(NA_real_, length(pc_names) - length(eig_use)))
    }
    # PLINK2 PCA projection via --score variance-standardize returns the
    # correct projected coordinates up to a global sign convention. For the
    # current allele-wts output, matching back to ref .eigenvec requires an
    # additional global sign flip; the sqrt(eigenvalue)/2 factor keeps the
    # projected scale aligned with the original reference PCs.
    div <- sqrt(eig_use) / 2
    div[!is.finite(div) | is.na(div) | div == 0] <- NA_real_
    for (i in seq_along(pc_names)) {
      v <- suppressWarnings(as.numeric(dt[[src[i]]]))
      out[[pc_names[i]]] <- -1 * v / div[i]
    }
    out <- out[!is.na(IID) & nzchar(IID)]
    if ("FID" %in% names(out)) {
      fidv <- as.character(out$FID)
      no_fid <- all(is.na(fidv) | !nzchar(fidv) | fidv %in% c("0", "NA"))
      if (isTRUE(no_fid)) out[, FID := NULL]
    }
    if ("FID" %in% names(out)) data.table::setcolorder(out, c("FID", "IID", "reference", pc_names))
    else data.table::setcolorder(out, c("IID", "reference", pc_names))
    out
  }
  .pca_projection_load_variant_pos <- function(ds) {
    if (identical(ds$format, "bfile")) {
      bim <- paste0(ds$prefix, ".bim")
      if (!file.exists(bim)) return(NULL)
      dtb <- data.table::fread(bim, data.table = TRUE, header = FALSE, showProgress = FALSE)
      if (ncol(dtb) < 4L) return(NULL)
      return(data.table::data.table(
        CHR = as.character(dtb[[1]]),
        ID = as.character(dtb[[2]]),
        POS = suppressWarnings(as.numeric(dtb[[4]]))
      ))
    }
    pvar <- if (file.exists(paste0(ds$prefix, ".pvar"))) paste0(ds$prefix, ".pvar") else if (file.exists(paste0(ds$prefix, ".pvar.zst"))) paste0(ds$prefix, ".pvar.zst") else NA_character_
    if (is.na(pvar) || !file.exists(pvar)) return(NULL)
    dtp <- tryCatch(suppressWarnings(data.table::fread(pvar, data.table = TRUE, showProgress = FALSE, quote = "")), error = function(e) NULL)
    if (is.null(dtp) || !nrow(dtp)) return(NULL)
    cn <- names(dtp)
    c_chr <- .pca_projection_pick_col(cn, c("#CHROM", "CHROM", "CHR"))
    c_pos <- .pca_projection_pick_col(cn, c("POS", "BP"))
    c_id <- .pca_projection_pick_col(cn, c("ID", "SNP"))
    if (is.na(c_chr) || is.na(c_pos) || is.na(c_id)) return(NULL)
    data.table::data.table(
      CHR = as.character(dtp[[c_chr]]),
      ID = as.character(dtp[[c_id]]),
      POS = suppressWarnings(as.numeric(dtp[[c_pos]]))
    )
  }
  .pca_projection_add_a1_freq <- function(dt) {
    if (is.null(dt) || !nrow(dt)) return(dt)
    if ("A1.freq" %in% names(dt)) {
      dt[, A1.freq := suppressWarnings(as.numeric(dt[["A1.freq"]]))]
      return(dt)
    }
    c_a1 <- .pca_projection_pick_col(names(dt), c("A1", "ALLELE1", "ALT", "EFFECTALLELE"))
    c_ref <- .pca_projection_pick_col(names(dt), c("REF"))
    c_alt <- .pca_projection_pick_col(names(dt), c("ALT"))
    c_af <- .pca_projection_pick_col(names(dt), c("ALT_FREQS", "ALTFREQS", "ALTFREQ", "AF"))
    if (is.na(c_a1) || is.na(c_ref) || is.na(c_alt) || is.na(c_af)) {
      dt[, A1.freq := NA_real_]
      return(dt)
    }
    a1v <- toupper(as.character(dt[[c_a1]]))
    refv <- toupper(as.character(dt[[c_ref]]))
    altv <- toupper(as.character(dt[[c_alt]]))
    afv <- suppressWarnings(as.numeric(dt[[c_af]]))
    a1f <- ifelse(!is.na(a1v) & !is.na(altv) & a1v == altv, afv,
                  ifelse(!is.na(a1v) & !is.na(refv) & a1v == refv, 1 - afv, NA_real_))
    dt[, A1.freq := a1f]
    dt
  }
  .pca_projection_merge_pos_ref <- function(dt, pos_ref, add_chr = FALSE) {
    if (is.null(dt) || !nrow(dt) || is.null(pos_ref) || !nrow(pos_ref)) return(dt)
    if (!all(c("ID", "POS") %in% names(pos_ref))) return(dt)
    pos_use <- if (isTRUE(add_chr) && ("CHR" %in% names(pos_ref))) {
      unique(pos_ref[, .(ID, CHR_ref = CHR, POS_ref = POS)], by = "ID")
    } else {
      unique(pos_ref[, .(ID, POS_ref = POS)], by = "ID")
    }
    merge(dt, pos_use, by = "ID", all.x = TRUE, sort = FALSE)
  }
  .pca_projection_normalize_weight <- function(x) {
    dt <- if (data.table::is.data.table(x)) data.table::copy(x) else data.table::as.data.table(x)
    if (!nrow(dt)) stop("Precomputed weight table is empty.", call. = FALSE)
    if ("#CHROM" %in% names(dt)) data.table::setnames(dt, "#CHROM", "CHR")
    c_id <- .pca_projection_pick_col(names(dt), c("ID", "SNP", "VAR", "VARIANT"))
    if (is.na(c_id)) stop("Precomputed weight must contain an ID column.", call. = FALSE)
    if (c_id != "ID") data.table::setnames(dt, c_id, "ID")
    c_a1 <- .pca_projection_pick_col(names(dt), c("A1", "ALLELE1", "ALT", "EFFECTALLELE"))
    if (is.na(c_a1)) stop("Precomputed weight must contain an A1-like column.", call. = FALSE)
    dt <- .pca_projection_add_a1_freq(dt)
    pc_nm <- names(dt)[grepl("^PC[0-9]+$", toupper(names(dt)))]
    if (!length(pc_nm)) pc_nm <- names(dt)[grepl("^PC[0-9]+", toupper(names(dt)))]
    if (!length(pc_nm)) stop("Precomputed weight must contain PC weight columns.", call. = FALSE)
    lead_cols <- intersect(c("CHR", "POS", "ID", "REF", "ALT", "A1", "A1.freq", "OBS_CT", "ALT_FREQS"), names(dt))
    rest_cols <- setdiff(names(dt), lead_cols)
    data.table::setcolorder(dt, c(lead_cols, rest_cols))
    dt
  }
  .pca_projection_normalize_afreq <- function(x) {
    if (is.null(x)) return(NULL)
    dt <- if (data.table::is.data.table(x)) data.table::copy(x) else data.table::as.data.table(x)
    if (!nrow(dt)) return(NULL)
    if ("#CHROM" %in% names(dt)) data.table::setnames(dt, "#CHROM", "CHR")
    c_chr <- .pca_projection_pick_col(names(dt), c("CHR", "#CHROM", "CHROM"))
    c_pos <- .pca_projection_pick_col(names(dt), c("POS", "BP"))
    c_id <- .pca_projection_pick_col(names(dt), c("ID", "SNP", "VAR", "VARIANT"))
    c_ref <- .pca_projection_pick_col(names(dt), c("REF"))
    c_alt <- .pca_projection_pick_col(names(dt), c("ALT"))
    c_af <- .pca_projection_pick_col(names(dt), c("ALT_FREQS", "ALTFREQS", "ALTFREQ", "AF"))
    c_a1 <- .pca_projection_pick_col(names(dt), c("A1", "ALLELE1", "EFFECTALLELE"))
    c_a1f <- .pca_projection_pick_col(names(dt), c("A1.FREQ", "A1_FREQ", "A1FREQ"))
    c_obs <- .pca_projection_pick_col(names(dt), c("OBS_CT", "OBSCT", "N", "NOBS"))
    if (is.na(c_af) && !is.na(c_id) && !is.na(c_ref) && !is.na(c_alt) && !is.na(c_a1) && !is.na(c_a1f)) {
      a1v <- toupper(as.character(dt[[c_a1]]))
      refv <- toupper(as.character(dt[[c_ref]]))
      altv <- toupper(as.character(dt[[c_alt]]))
      a1fv <- suppressWarnings(as.numeric(dt[[c_a1f]]))
      altf <- ifelse(!is.na(a1v) & !is.na(altv) & a1v == altv, a1fv,
                     ifelse(!is.na(a1v) & !is.na(refv) & a1v == refv, 1 - a1fv, NA_real_))
      dt[, ALT_FREQS := altf]
      c_af <- "ALT_FREQS"
    }
    if (is.na(c_id) || is.na(c_ref) || is.na(c_alt) || is.na(c_af)) return(NULL)
    out <- data.table::data.table(
      ID = as.character(dt[[c_id]]),
      REF = as.character(dt[[c_ref]]),
      ALT = as.character(dt[[c_alt]]),
      ALT_FREQS = suppressWarnings(as.numeric(dt[[c_af]]))
    )
    if (!is.na(c_chr)) out[, CHR := as.character(dt[[c_chr]])]
    if (!is.na(c_pos)) out[, POS := suppressWarnings(as.numeric(dt[[c_pos]]))]
    out[, OBS_CT := if (!is.na(c_obs)) suppressWarnings(as.numeric(dt[[c_obs]])) else NA_real_]
    out <- out[!is.na(ID) & nzchar(ID)]
    out <- unique(out, by = "ID")
    lead_cols <- intersect(c("CHR", "POS", "ID", "REF", "ALT", "ALT_FREQS", "OBS_CT"), names(out))
    rest_cols <- setdiff(names(out), lead_cols)
    data.table::setcolorder(out, c(lead_cols, rest_cols))
    out
  }
  .pca_projection_afreq_from_weight <- function(weight) {
    .pca_projection_normalize_afreq(weight)
  }
  .pca_projection_normalize_ref_pc <- function(x, pc_names) {
    if (!(data.table::is.data.table(x) || is.data.frame(x))) return(NULL)
    dt <- if (data.table::is.data.table(x)) data.table::copy(x) else data.table::as.data.table(x)
    if (!nrow(dt)) return(NULL)
    if ("#FID" %in% names(dt)) data.table::setnames(dt, "#FID", "FID")
    if ("#IID" %in% names(dt)) data.table::setnames(dt, "#IID", "IID")
    c_iid <- .pca_projection_pick_col(names(dt), c("IID", "ID", "#IID"))
    if (is.na(c_iid)) return(NULL)
    if (c_iid != "IID") data.table::setnames(dt, c_iid, "IID")
    c_fid <- .pca_projection_pick_col(names(dt), c("FID", "#FID"))
    if ("reference" %in% names(dt)) {
      ref_flag <- as.logical(dt$reference)
      if (any(ref_flag %in% TRUE, na.rm = TRUE)) {
        dt <- dt[ref_flag %in% TRUE]
      }
    }
    if (!nrow(dt)) return(NULL)
    missing_pc <- setdiff(pc_names, names(dt))
    if (length(missing_pc)) {
      stop(
        "Precomputed ref PC table is missing projected PCs: ",
        paste(missing_pc, collapse = ", "),
        call. = FALSE
      )
    }
    out <- data.table::data.table(IID = as.character(dt$IID), reference = TRUE)
    if (!is.na(c_fid)) out[, FID := as.character(dt[[c_fid]])]
    for (pcn in pc_names) out[[pcn]] <- suppressWarnings(as.numeric(dt[[pcn]]))
    out <- out[!is.na(IID) & nzchar(IID)]
    if ("FID" %in% names(out)) {
      fidv <- as.character(out$FID)
      no_fid <- all(is.na(fidv) | !nzchar(fidv) | fidv %in% c("0", "NA"))
      if (isTRUE(no_fid)) out[, FID := NULL]
    }
    if ("FID" %in% names(out)) data.table::setcolorder(out, c("FID", "IID", "reference", pc_names))
    else data.table::setcolorder(out, c("IID", "reference", pc_names))
    out
  }
  .pca_projection_extract_ref_input <- function(x, prefer = "pfile") {
    if (data.table::is.data.table(x) || is.data.frame(x)) {
      weight_dt <- .pca_projection_normalize_weight(x)
      eig <- attr(x, "gcanvas_pca_eigenvalue", exact = TRUE)
      prop <- attr(x, "gcanvas_pca_prop.var", exact = TRUE)
      afreq <- attr(x, "gcanvas_pca_afreq", exact = TRUE)
      if (is.null(afreq)) afreq <- .pca_projection_afreq_from_weight(weight_dt)
      return(list(
        type = "precomputed",
        weight = weight_dt,
        eigenvalue = eig,
        prop.var = prop,
        afreq = afreq,
        pc = NULL
      ))
    }
    if (is.list(x) && !is.null(x$weight) && (data.table::is.data.table(x$weight) || is.data.frame(x$weight))) {
      weight_dt <- .pca_projection_normalize_weight(x$weight)
      eig <- x$eigenvalue %||% attr(x$weight, "gcanvas_pca_eigenvalue", exact = TRUE)
      prop <- x$prop.var %||% attr(x$weight, "gcanvas_pca_prop.var", exact = TRUE)
      afreq <- x$afreq %||% attr(x$weight, "gcanvas_pca_afreq", exact = TRUE)
      if (is.null(afreq)) afreq <- .pca_projection_afreq_from_weight(weight_dt)
      return(list(
        type = "precomputed",
        weight = weight_dt,
        eigenvalue = eig,
        prop.var = prop,
        afreq = afreq,
        pc = x$pc %||% NULL
      ))
    }
    list(type = "dataset", dataset = .pca_projection_detect_dataset(x, tag = "ref", prefer = prefer))
  }
  .pca_projection_build_scree <- function(prop.var) {
    .gcanvas_pca_build_screeplot(prop.var)
  }

  pf <- tolower(trimws(as.character(plink.format)[1]))
  if (!nzchar(pf) || is.na(pf)) pf <- "pfile"
  if (!(pf %in% c("pfile", "bfile", "auto"))) {
    stop("plink.format must be one of 'pfile', 'bfile', or 'auto'.", call. = FALSE)
  }

  pv <- tolower(trimws(as.character(plink.version)[1]))
  if (!nzchar(pv) || is.na(pv)) pv <- "auto"
  if (pv %in% c("plink1", "1", "1.9")) {
    stop("pca.projection requires PLINK2.", call. = FALSE)
  }

  ref_input <- .pca_projection_extract_ref_input(ref, prefer = pf)
  target_ds <- .pca_projection_detect_dataset(target, tag = "target", prefer = pf)
  ref_ds <- if (identical(ref_input$type, "dataset")) ref_input$dataset else NULL

  chr_vec <- .pca_projection_normalize_chrom(chrom)
  if (!length(chr_vec)) chr_vec <- as.character(1:22)
  chr0 <- paste(chr_vec, collapse = ",")
  maf0 <- as_num(maf)[1]; if (!is.finite(maf0) || is.na(maf0) || maf0 < 0 || maf0 > 1) maf0 <- 0.05
  geno0 <- as_num(geno)[1]; if (!is.finite(geno0) || is.na(geno0) || geno0 < 0 || geno0 > 1) geno0 <- 0.02
  indep_window0 <- as_int(indep.window)[1]; if (!is.finite(indep_window0) || is.na(indep_window0) || indep_window0 < 1L) indep_window0 <- 50L
  indep_step0 <- as_int(indep.step)[1]; if (!is.finite(indep_step0) || is.na(indep_step0) || indep_step0 < 1L) indep_step0 <- 5L
  indep_r20 <- as_num(indep.r2)[1]; if (!is.finite(indep_r20) || is.na(indep_r20) || indep_r20 <= 0 || indep_r20 >= 1) indep_r20 <- 0.1
  mode0 <- tolower(trimws(as.character(pca.mode)[1])); if (is.na(mode0) || !nzchar(mode0)) mode0 <- "approx"

  plink_cmd <- as.character(plink)[1]
  if (is.na(plink_cmd) || !nzchar(plink_cmd)) stop("PLINK command is empty.", call. = FALSE)
  plink_exec <- .pca_projection_resolve_exec(plink_cmd)
  if (is.na(plink_exec) || !nzchar(plink_exec)) {
    stop(sprintf("PLINK executable not found for '%s'. Set `plink` to an executable path.", plink_cmd), call. = FALSE)
  }

  cache_dir0 <- .gcanvas_default_cache_dir(
    scope = "pca_projection",
    anchor = if (!is.null(ref_ds)) ref_ds$prefix else target_ds$prefix
  )
  dir.create(cache_dir0, recursive = TRUE, showWarnings = FALSE)
  .gcanvas_register_cache_dir(cache_dir0)
  base_pref <- file.path(cache_dir0, sprintf("pca_projection_%s_%s", format(Sys.time(), "%Y%m%d_%H%M%S"), as.integer(stats::runif(1, 1, 1e9))))
  keep_plink <- FALSE
  if (!is.null(plink.out) && length(plink.out) > 0L) {
    if (is.logical(plink.out) && length(plink.out) == 1L) {
      keep_plink <- isTRUE(plink.out)
    } else {
      out_chr <- as.character(plink.out)[1]
      if (!is.na(out_chr) && nzchar(out_chr) && !(tolower(trimws(out_chr)) %in% c("false", "null", "na"))) {
        out_chr <- abs_path(path.expand(out_chr))
        if (dir.exists(out_chr) || grepl("[/\\\\]$", as.character(plink.out)[1])) {
          dir.create(out_chr, recursive = TRUE, showWarnings = FALSE)
          .gcanvas_register_cache_dir(out_chr)
          base_pref <- file.path(out_chr, sprintf("pca_projection_%s_%s", format(Sys.time(), "%Y%m%d_%H%M%S"), as.integer(stats::runif(1, 1, 1e9))))
        } else {
          dir.create(dirname(out_chr), recursive = TRUE, showWarnings = FALSE)
          .gcanvas_register_cache_dir(dirname(out_chr))
          base_pref <- out_chr
        }
        keep_plink <- TRUE
      }
    }
  }

  .gcanvas_note(
    "gcanvas::pca.projection",
    sprintf("Start: ref=%s | target=%s | plink=%s | plink.version=%s | threads=%d",
            if (!is.null(ref_ds)) ref_ds$format else "precomputed", target_ds$format, plink_cmd, pv, as_int(thr)),
    silent = silent
  )
  if (isTRUE(keep_plink)) {
    .gcanvas_note("gcanvas::pca.projection", sprintf("PLINK outputs kept: %s", base_pref), silent = silent)
  }

  if (!identical(ref_input$type, "dataset")) {
    target_eff <- target_ds
    target_align_pref <- paste0(base_pref, ".target.align")
    if (!id.aligned) {
      .gcanvas_note("gcanvas::pca.projection", "Aligning target variant IDs with --set-all-var-ids '@:#:$1:$2'", silent = silent)
      args_target_align <- c(
        .pca_projection_dataset_args(target_ds),
        "--set-all-var-ids", "@:#:$1:$2",
        "--new-id-max-allele-len", "1000", "missing",
        "--make-pgen",
        "--threads", as.character(as_int(thr)),
        "--out", target_align_pref
      )
      .pca_projection_run_plink(plink_exec, args_target_align)
      target_eff <- list(prefix = target_align_pref, format = "pfile", use_vzs = FALSE, tag = "target")
    }

    weight_full <- data.table::copy(ref_input$weight)
    weight_ids_all <- unique(as.character(weight_full$ID))
    weight_ids_all <- weight_ids_all[!is.na(weight_ids_all) & nzchar(weight_ids_all)]
    if (!length(weight_ids_all)) stop("Precomputed weight has no valid variant IDs.", call. = FALSE)

    eig_all <- as_num(ref_input$eigenvalue)
    eig_all <- eig_all[is.finite(eig_all) & !is.na(eig_all)]
    if (!length(eig_all)) {
      stop(
        "Precomputed ref must provide eigenvalue. Use a pca()/pca.projection() object, or pass weight with gcanvas_pca_eigenvalue attribute.",
        call. = FALSE
      )
    }
    if (is.null(names(eig_all)) || !length(names(eig_all)) || any(is.na(names(eig_all)) | !nzchar(names(eig_all)))) {
      names(eig_all) <- paste0("PC", seq_along(eig_all))
    }

    prop_all <- ref_input$prop.var
    if (!is.null(prop_all)) {
      prop_all <- as_num(prop_all)
      if (length(prop_all) >= 1L) {
        if (is.null(names(prop_all)) || any(is.na(names(prop_all)) | !nzchar(names(prop_all)))) {
          names(prop_all) <- names(eig_all)[seq_len(min(length(prop_all), length(eig_all)))]
        }
      }
    }

    afreq_full <- .pca_projection_normalize_afreq(ref_input$afreq)
    if (is.null(afreq_full)) afreq_full <- .pca_projection_afreq_from_weight(weight_full)
    if (is.null(afreq_full) || !nrow(afreq_full)) {
      stop(
        "Precomputed ref must provide reference allele frequencies. Rerun pca()/pca.projection() with the current gcanvas version, or supply raw reference data.",
        call. = FALSE
      )
    }

    target_ids <- .pca_projection_read_variant_ids(target_eff)
    missing_ids <- setdiff(weight_ids_all, target_ids)
    if (length(missing_ids) && !force) {
      stop(
        sprintf(
          "Target is missing %d/%d weight variants. Set force=TRUE to drop missing weight variants and continue.",
          as_int(length(missing_ids)),
          as_int(length(weight_ids_all))
        ),
        call. = FALSE
      )
    }
    if (length(missing_ids) && force) {
      warning(
        sprintf(
          "Target is missing %d/%d weight variants; dropping them because force=TRUE.",
          as_int(length(missing_ids)),
          as_int(length(weight_ids_all))
        ),
        call. = FALSE
      )
    }
    score_ids <- setdiff(weight_ids_all, missing_ids)
    if (!length(score_ids)) stop("No weight variants remain after matching to target.", call. = FALSE)

    weight_use <- weight_full[ID %chin% score_ids]
    afreq_use <- afreq_full[ID %chin% score_ids]
    afreq_ids <- unique(as.character(afreq_use$ID))
    miss_af <- setdiff(score_ids, afreq_ids)
    if (length(miss_af)) {
      stop(
        sprintf(
          "Reference allele frequency table is missing %d/%d projected variants.",
          as_int(length(miss_af)),
          as_int(length(score_ids))
        ),
        call. = FALSE
      )
    }

    weight_file <- paste0(base_pref, ".ref.weight.tsv")
    afreq_file <- paste0(base_pref, ".ref.afreq.tsv")
    extract_file <- paste0(base_pref, ".weight.extract.txt")
    c_a1_w <- .pca_projection_pick_col(names(weight_use), c("A1", "ALLELE1", "ALT", "EFFECTALLELE"))
    if (is.na(c_a1_w)) stop("Precomputed weight is missing an A1-like column for scoring.", call. = FALSE)
    pc_write <- names(weight_use)[grepl("^PC[0-9]+$", toupper(names(weight_use)))]
    if (!length(pc_write)) pc_write <- names(weight_use)[grepl("^PC[0-9]+", toupper(names(weight_use)))]
    if (!length(pc_write)) stop("Precomputed weight has no PC columns for scoring.", call. = FALSE)
    pc_write <- pc_write[order(as_int(gsub("[^0-9]", "", toupper(pc_write))), na.last = TRUE)]
    weight_write <- data.table::data.table(
      ID = as.character(weight_use$ID),
      A1 = as.character(weight_use[[c_a1_w]])
    )
    for (pcn in pc_write) weight_write[[pcn]] <- suppressWarnings(as.numeric(weight_use[[pcn]]))
    data.table::fwrite(weight_write, weight_file, sep = "\t")
    afreq_write <- data.table::copy(afreq_use)
    if (!("CHR" %in% names(afreq_write))) afreq_write[, CHR := NA_character_]
    if (!("OBS_CT" %in% names(afreq_write))) afreq_write[, OBS_CT := NA_real_]
    keep_afreq_cols <- intersect(c("CHR", "ID", "REF", "ALT", "ALT_FREQS", "OBS_CT"), names(afreq_write))
    afreq_write <- afreq_write[, ..keep_afreq_cols]
    data.table::setcolorder(afreq_write, c(intersect("CHR", names(afreq_write)), setdiff(names(afreq_write), "CHR")))
    if ("CHR" %in% names(afreq_write)) data.table::setnames(afreq_write, "CHR", "#CHROM")
    data.table::fwrite(afreq_write, afreq_file, sep = "\t")
    data.table::fwrite(data.table::data.table(SNP = score_ids), extract_file, sep = "\t", col.names = FALSE, quote = FALSE)

    score_spec <- .pca_projection_parse_score_spec(weight_file, pca_count = pca.count)
    if (is.na(score_spec$score_spec) || !nzchar(score_spec$score_spec)) {
      stop("Failed to build --score-col-nums from precomputed weight.", call. = FALSE)
    }
    pc_names <- as.character(score_spec$pc_names)
    if (!all(pc_names %in% names(eig_all))) {
      if (length(eig_all) < length(pc_names)) {
        stop("Precomputed eigenvalue does not cover all PC columns in weight.", call. = FALSE)
      }
      eig <- eig_all[seq_len(length(pc_names))]
      names(eig) <- pc_names
    } else {
      eig <- eig_all[pc_names]
    }
    if (!is.null(prop_all) && length(prop_all) && all(pc_names %in% names(prop_all))) {
      prop.var <- .gcanvas_pca_propvar_from_eval(prop_all)[pc_names]
    } else {
      prop.var <- .gcanvas_pca_propvar_from_eval(eig_all)[pc_names]
      names(prop.var) <- pc_names
    }

    target_proj_pref <- paste0(base_pref, ".target.proj")
    args_score_common <- c(
      "--extract", extract_file,
      "--read-freq", afreq_file,
      "--score", weight_file, as.character(score_spec$id_col), as.character(score_spec$a1_col),
      "header-read", "no-mean-imputation", "variance-standardize", "list-variants",
      "--score-col-nums", score_spec$score_spec,
      "--threads", as.character(as_int(thr))
    )
    .pca_projection_run_plink(plink_exec, c(.pca_projection_dataset_args(target_eff), args_score_common, "--out", target_proj_pref))

    pc_target <- .pca_projection_read_sscore(target_proj_pref, pc_names = pc_names, eig_vals = eig, is_ref = FALSE)
    pc_ref <- .pca_projection_normalize_ref_pc(ref_input$pc, pc_names = pc_names)
    pc <- if (is.null(pc_ref) || !nrow(pc_ref)) {
      pc_target
    } else {
      data.table::rbindlist(list(pc_ref, pc_target), use.names = TRUE, fill = TRUE)
    }
    if ("FID" %in% names(pc)) data.table::setcolorder(pc, c("FID", "IID", "reference", pc_names))
    else data.table::setcolorder(pc, c("IID", "reference", pc_names))

    weight_use <- .pca_projection_add_a1_freq(weight_use)
    if ("ALT_FREQS" %in% names(weight_use)) weight_use[, ALT_FREQS := NULL]
    lead_cols <- intersect(c("CHR", "POS", "ID", "REF", "ALT", "PROVISIONAL_REF?", "A1", "A1.freq", "OBS_CT"), names(weight_use))
    rest_cols <- setdiff(names(weight_use), lead_cols)
    data.table::setcolorder(weight_use, c(lead_cols, rest_cols))
    attr(weight_use, "gcanvas_pca_eigenvalue") <- eig
    attr(weight_use, "gcanvas_pca_prop.var") <- prop.var
    attr(weight_use, "gcanvas_pca_afreq") <- afreq_use

    scree <- .pca_projection_build_scree(prop.var)

    all_temp <- unique(c(
      weight_file,
      afreq_file,
      extract_file,
      Sys.glob(paste0(target_proj_pref, ".*")),
      Sys.glob(paste0(target_align_pref, ".*"))
    ))
    all_temp <- all_temp[file.exists(all_temp)]
    keep_files <- character()
    if (isTRUE(keep_plink)) {
      keep_files <- unique(c(
        paste0(target_proj_pref, ".sscore"),
        weight_file,
        afreq_file
      ))
      keep_files <- keep_files[file.exists(keep_files)]
      if (length(keep_files)) .gcanvas_register_cache_files(keep_files)
    }
    rm_files <- setdiff(all_temp, keep_files)
    if (length(rm_files)) unlink(rm_files, force = TRUE)

    if (isTRUE(show.plinklog) && !isTRUE(keep_plink)) {
      .gcanvas_note(
        "gcanvas::pca.projection",
        "PLINK intermediate files were removed (plink.out=FALSE). Set plink.out=TRUE or a path to keep projection outputs.",
        silent = silent
      )
    }

    .gcanvas_note(
      "gcanvas::pca.projection",
      sprintf(
        "Done: n_ref=%d | n_target=%d | n_variants=%d | n_pcs=%d",
        as_int(sum(as.logical(pc$reference) %in% TRUE, na.rm = TRUE)),
        as_int(sum(as.logical(pc$reference) %in% FALSE, na.rm = TRUE)),
        as_int(length(score_ids)),
        as_int(length(pc_names))
      ),
      silent = silent
    )

    out <- list(
      pc = pc,
      eigenvalue = eig,
      weight = weight_use,
      prop.var = prop.var,
      screeplot = scree
    )
    class(out) <- c("gcanvas_pca_projection", "gcanvas_pca", "list")
    return(out)
  }

  ref_eff <- ref_ds
  target_eff <- target_ds
  ref_align_pref <- paste0(base_pref, ".ref.align")
  target_align_pref <- paste0(base_pref, ".target.align")
  if (!id.aligned) {
    .gcanvas_note("gcanvas::pca.projection", "Aligning variant IDs with --set-all-var-ids '@:#:$1:$2'", silent = silent)
    args_ref_align <- c(
      .pca_projection_dataset_args(ref_ds),
      "--set-all-var-ids", "@:#:$1:$2",
      "--new-id-max-allele-len", "1000", "missing",
      "--make-pgen",
      "--threads", as.character(as_int(thr)),
      "--out", ref_align_pref
    )
    args_target_align <- c(
      .pca_projection_dataset_args(target_ds),
      "--set-all-var-ids", "@:#:$1:$2",
      "--new-id-max-allele-len", "1000", "missing",
      "--make-pgen",
      "--threads", as.character(as_int(thr)),
      "--out", target_align_pref
    )
    .pca_projection_run_plink(plink_exec, args_ref_align)
    .pca_projection_run_plink(plink_exec, args_target_align)
    ref_eff <- list(prefix = ref_align_pref, format = "pfile", use_vzs = FALSE, tag = "ref")
    target_eff <- list(prefix = target_align_pref, format = "pfile", use_vzs = FALSE, tag = "target")
  }

  ref_ids <- .pca_projection_read_variant_ids(ref_eff)
  target_ids <- .pca_projection_read_variant_ids(target_eff)
  common_ids <- intersect(ref_ids, target_ids)
  if (!length(common_ids)) stop("No intersected variants between ref and target.", call. = FALSE)
  snplist_file <- paste0(base_pref, ".intersect.snplist")
  data.table::fwrite(data.table::data.table(SNP = common_ids), snplist_file, sep = "\t", col.names = FALSE, quote = FALSE)
  .gcanvas_note("gcanvas::pca.projection", sprintf("Intersected variants: n=%d", as_int(length(common_ids))), silent = silent)

  prune_pref <- paste0(base_pref, ".prune")
  args_prune <- c(
    .pca_projection_dataset_args(ref_eff),
    "--extract", snplist_file,
    "--chr", chr0,
    "--maf", format(maf0, scientific = FALSE),
    "--geno", format(geno0, scientific = FALSE),
    "--indep-pairwise", as.character(indep_window0), as.character(indep_step0), format(indep_r20, scientific = FALSE),
    "--threads", as.character(as_int(thr)),
    "--out", prune_pref
  )
  .pca_projection_run_plink(plink_exec, args_prune)
  prune_in <- paste0(prune_pref, ".prune.in")
  if (!file.exists(prune_in)) stop("Missing prune.in file after indep-pairwise.", call. = FALSE)

  pca_pref <- paste0(base_pref, ".ref.pca")
  args_pca <- c(
    .pca_projection_dataset_args(ref_eff),
    "--extract", prune_in,
    "--freq",
    "--pca", mode0, "allele-wts",
    "--threads", as.character(as_int(thr)),
    "--out", pca_pref
  )
  if (!is.null(pca.count) && length(pca.count) > 0L && !is.na(pca.count[1])) {
    k <- as_int(pca.count)[1]
    if (is.finite(k) && !is.na(k) && k > 0L) {
      args_pca <- c(
        .pca_projection_dataset_args(ref_eff),
        "--extract", prune_in,
        "--freq",
        "--pca", mode0, as.character(k), "allele-wts",
        "--threads", as.character(as_int(thr)),
        "--out", pca_pref
      )
    }
  }
  .pca_projection_run_plink(plink_exec, args_pca)

  ev_file <- paste0(pca_pref, ".eigenval")
  wt_file <- if (file.exists(paste0(pca_pref, ".eigenvec.allele"))) {
    paste0(pca_pref, ".eigenvec.allele")
  } else if (file.exists(paste0(pca_pref, ".eigenvec.var"))) {
    paste0(pca_pref, ".eigenvec.var")
  } else {
    NA_character_
  }
  afreq_file <- paste0(pca_pref, ".afreq")
  if (!file.exists(ev_file) || is.na(wt_file) || !file.exists(wt_file) || !file.exists(afreq_file)) {
    stop("Reference PCA step finished but required files (.eigenval/.eigenvec.allele(or .var)/.afreq) are missing.", call. = FALSE)
  }

  score_spec <- .pca_projection_parse_score_spec(wt_file, pca_count = pca.count)
  if (is.na(score_spec$score_spec) || !nzchar(score_spec$score_spec)) {
    stop("Failed to build --score-col-nums from weight file.", call. = FALSE)
  }

  ref_proj_pref <- paste0(base_pref, ".ref.proj")
  target_proj_pref <- paste0(base_pref, ".target.proj")
  args_score_common <- c(
    "--extract", prune_in,
    "--read-freq", afreq_file,
    "--score", wt_file, as.character(score_spec$id_col), as.character(score_spec$a1_col),
    "header-read", "no-mean-imputation", "variance-standardize", "list-variants",
    "--score-col-nums", score_spec$score_spec,
    "--threads", as.character(as_int(thr))
  )
  .pca_projection_run_plink(plink_exec, c(.pca_projection_dataset_args(ref_eff), args_score_common, "--out", ref_proj_pref))
  .pca_projection_run_plink(plink_exec, c(.pca_projection_dataset_args(target_eff), args_score_common, "--out", target_proj_pref))

  eig_dt <- data.table::fread(ev_file, data.table = TRUE, header = FALSE, showProgress = FALSE)
  if (!ncol(eig_dt)) stop("Empty eigenvalue file: ", ev_file, call. = FALSE)
  eig <- as_num(eig_dt[[1]])
  eig <- eig[is.finite(eig) & !is.na(eig)]
  if (!length(eig)) stop("No finite eigenvalues found.", call. = FALSE)
  pc_names <- as.character(score_spec$pc_names)
  k_use <- min(length(pc_names), length(eig))
  pc_names <- pc_names[seq_len(k_use)]
  eig <- eig[seq_len(k_use)]
  names(eig) <- pc_names
  prop.var <- .gcanvas_pca_propvar_from_eval(eig)
  names(prop.var) <- pc_names

  pc_ref <- .pca_projection_read_sscore(ref_proj_pref, pc_names = pc_names, eig_vals = eig, is_ref = TRUE)
  pc_target <- .pca_projection_read_sscore(target_proj_pref, pc_names = pc_names, eig_vals = eig, is_ref = FALSE)
  pc <- data.table::rbindlist(list(pc_ref, pc_target), use.names = TRUE, fill = TRUE)
  if ("FID" %in% names(pc)) data.table::setcolorder(pc, c("FID", "IID", "reference", pc_names))
  else data.table::setcolorder(pc, c("IID", "reference", pc_names))

  weight <- data.table::fread(wt_file, data.table = TRUE, check.names = FALSE, showProgress = FALSE)
  if ("#CHROM" %in% names(weight)) data.table::setnames(weight, "#CHROM", "CHR")
  c_id_w <- .pca_projection_pick_col(names(weight), c("ID", "SNP", "VAR"))
  if (is.na(c_id_w)) c_id_w <- names(weight)[1]
  if (c_id_w != "ID") data.table::setnames(weight, c_id_w, "ID")
  afreq <- data.table::fread(afreq_file, data.table = TRUE, check.names = FALSE, showProgress = FALSE)
  if ("#CHROM" %in% names(afreq)) data.table::setnames(afreq, "#CHROM", "CHR")
  c_id_af <- .pca_projection_pick_col(names(afreq), c("ID", "SNP", "VAR"))
  if (!is.na(c_id_af) && c_id_af != "ID") data.table::setnames(afreq, c_id_af, "ID")
  if ("ID" %in% names(afreq)) afreq <- afreq[!is.na(ID) & nzchar(ID)]
  pos_ref <- .pca_projection_load_variant_pos(ref_eff)
  if (!is.null(pos_ref) && nrow(pos_ref)) {
    weight <- .pca_projection_merge_pos_ref(weight, pos_ref, add_chr = TRUE)
    if (!is.null(afreq) && nrow(afreq)) {
      afreq <- .pca_projection_merge_pos_ref(afreq, pos_ref, add_chr = FALSE)
      if (!("POS" %in% names(afreq))) afreq[, POS := POS_ref]
      if ("POS_ref" %in% names(afreq)) afreq[, POS_ref := NULL]
    }
  } else {
    if (!("CHR_ref" %in% names(weight))) weight[, CHR_ref := NA_character_]
    if (!("POS_ref" %in% names(weight))) weight[, POS_ref := NA_real_]
  }
  if ("CHR.x" %in% names(weight) || "CHR.y" %in% names(weight)) {
    chrx <- if ("CHR.x" %in% names(weight)) as.character(weight[["CHR.x"]]) else rep(NA_character_, nrow(weight))
    chry <- if ("CHR.y" %in% names(weight)) as.character(weight[["CHR.y"]]) else rep(NA_character_, nrow(weight))
    weight[, CHR := data.table::fifelse(!is.na(chrx) & nzchar(chrx), chrx, chry)]
    rm_chrxy <- intersect(c("CHR.x", "CHR.y"), names(weight))
    if (length(rm_chrxy)) weight[, (rm_chrxy) := NULL]
  }
  if (!("CHR" %in% names(weight))) weight[, CHR := NA_character_]
  if ("CHR_ref" %in% names(weight)) {
    wchr <- as.character(weight$CHR)
    cref <- as.character(weight$CHR_ref)
    fill_chr <- (is.na(wchr) | !nzchar(wchr)) & !is.na(cref) & nzchar(cref)
    if (any(fill_chr)) weight[fill_chr, CHR := CHR_ref]
    weight[, CHR_ref := NULL]
  }
  if (!("POS" %in% names(weight))) weight[, POS := NA_real_]
  if ("POS_ref" %in% names(weight)) {
    p0 <- suppressWarnings(as.numeric(weight$POS))
    pr <- suppressWarnings(as.numeric(weight$POS_ref))
    fill_pos <- !is.finite(p0) & is.finite(pr)
    if (any(fill_pos)) weight[fill_pos, POS := POS_ref]
    weight[, POS_ref := NULL]
  }
  weight[, POS := suppressWarnings(as.numeric(POS))]
  if (!is.null(afreq) && nrow(afreq)) {
    if ("POS" %in% names(afreq)) afreq[, POS := suppressWarnings(as.numeric(POS))]
    afreq <- unique(afreq, by = "ID")
    keep_af_cols <- intersect(c("ID", "ALT_FREQS", "OBS_CT"), names(afreq))
    if (length(keep_af_cols) >= 2L) {
      weight <- merge(weight, afreq[, ..keep_af_cols], by = "ID", all.x = TRUE, sort = FALSE)
    }
  }
  weight <- .pca_projection_add_a1_freq(weight)
  if ("ALT_FREQS" %in% names(weight)) weight[, ALT_FREQS := NULL]
  lead_cols <- c("CHR", "POS", "ID", "REF", "ALT", "PROVISIONAL_REF?", "A1", "A1.freq", "OBS_CT")
  rest_cols <- setdiff(names(weight), lead_cols)
  data.table::setcolorder(weight, c(lead_cols, rest_cols))
  attr(weight, "gcanvas_pca_eigenvalue") <- eig
  attr(weight, "gcanvas_pca_prop.var") <- prop.var
  attr(weight, "gcanvas_pca_afreq") <- afreq

  scree <- .pca_projection_build_scree(prop.var)

  all_temp <- unique(c(
    snplist_file,
    Sys.glob(paste0(prune_pref, ".*")),
    Sys.glob(paste0(pca_pref, ".*")),
    Sys.glob(paste0(ref_proj_pref, ".*")),
    Sys.glob(paste0(target_proj_pref, ".*")),
    Sys.glob(paste0(ref_align_pref, ".*")),
    Sys.glob(paste0(target_align_pref, ".*"))
  ))
  all_temp <- all_temp[file.exists(all_temp)]
  keep_files <- character()
  if (isTRUE(keep_plink)) {
    keep_files <- unique(c(
      paste0(ref_proj_pref, ".sscore"),
      paste0(target_proj_pref, ".sscore"),
      ev_file,
      wt_file
    ))
    keep_files <- keep_files[file.exists(keep_files)]
    if (length(keep_files)) .gcanvas_register_cache_files(keep_files)
  }
  rm_files <- setdiff(all_temp, keep_files)
  if (length(rm_files)) unlink(rm_files, force = TRUE)

  if (isTRUE(show.plinklog) && !isTRUE(keep_plink)) {
    .gcanvas_note(
      "gcanvas::pca.projection",
      "PLINK intermediate files were removed (plink.out=FALSE). Set plink.out=TRUE or a path to keep projection outputs.",
      silent = silent
    )
  }

  .gcanvas_note(
    "gcanvas::pca.projection",
    sprintf("Done: n_ref=%d | n_target=%d | n_variants=%d | n_pcs=%d",
            as_int(sum(as.logical(pc$reference) %in% TRUE, na.rm = TRUE)),
            as_int(sum(as.logical(pc$reference) %in% FALSE, na.rm = TRUE)),
            as_int(length(common_ids)),
            as_int(length(pc_names))),
    silent = silent
  )

  out <- list(
    pc = pc,
    eigenvalue = eig,
    weight = weight,
    prop.var = prop.var,
    screeplot = scree
  )
  class(out) <- c("gcanvas_pca_projection", "gcanvas_pca", "list")
  out
}

#' Plot PCA results
#'
#' Renders 2D / 3D / facet-grid PC scatter plots with optional category
#' coloring, density overlays, projection comparison, and reference shading.
#'
#' @param pc A PCA result list (output of [pca()] or [pca.projection()]).
#' @param eval Optional eigenvalue vector for variance-explained axis labels.
#' @param pc.use Which PCs to plot (`"all"` or integer vector).
#' @param id.col Sample-id column name in `pc`.
#' @param group.data,group.col,group.title,group.color,group.shape,group.size,group.stroke
#'   Category / grouping controls.
#' @param point.color,point.size,point.stroke,point.shape,alpha Base point styling.
#' @param outlier,ambiguous,fade Optional outlier / ambiguous-sample annotations.
#' @param line.sd,line.sd.color,line.sd.type,line.sd.linewidth Standard-deviation
#'   ellipse / line controls.
#' @param legend.title,legend.position Legend controls.
#' @param equal.lim One of `"panel"`, `"global"`, `"none"` — axis-limit policy.
#' @param projection.helper,projection.group Helpers for plotting projected
#'   target samples on top of a reference PCA.
#'
#' @return A `ggplot` object (or list of them for multi-panel modes).
#' @export
pca.plot <- function(pc,
                     eval = NULL,
                     pc.use = "all",
                     id.col = "IID",
                     group.data = NULL,
                     group.col = NULL,
                     group.title = NULL,
                     group.color = NULL,
                     group.shape = NULL,
                     group.size = NULL,
                     group.stroke = NULL,
                     point.color = "grey20",
                     point.size = 1.5,
                     point.stroke = 1,
                     point.shape = 16,
                     alpha = 0.65,
                     outlier = NULL,
                     ambiguous = NULL,
                     fade = NULL,
                     line.sd = FALSE,
                     line.sd.color = NULL,
                     line.sd.type = "dashed",
                     line.sd.linewidth = 0.3,
                     legend.title = NULL,
                     legend.position = c("bottom", "top", "left", "right"),
                     equal.lim = c("panel", "global", "none"),
                     projection.helper = FALSE,
                     projection.group = NULL) {
  require_pkg(c("data.table", "ggplot2", "scales"))
  legend.position <- match.arg(legend.position)
  equal.lim <- match.arg(equal.lim)
  projection.helper <- isTRUE(projection.helper)
  refqc.helper <- FALSE
  refqc_drop_ids <- character()
  refqc_ambiguous_ids <- character()
  refqc_group_title <- NULL
  refqc_drop_label <- "Dropped"
  point_color_missing <- missing(point.color)
  user_group_style_input <- (!is.null(group.color) && length(group.color) > 0L) ||
    (!is.null(group.shape) && length(group.shape) > 0L) ||
    (!is.null(group.size) && length(group.size) > 0L) ||
    (!is.null(group.stroke) && length(group.stroke) > 0L)

  .pca_plot_numvec <- function(x) suppressWarnings(as.numeric(x))
  .pca_plot_as_char1 <- function(x) {
    if (is.null(x) || !length(x)) return(NULL)
    v <- as.character(x)[1]
    if (is.na(v) || !nzchar(v)) NULL else v
  }
  .pca_plot_as_flag1 <- function(x, default = NA) {
    if (is.null(x) || !length(x)) return(default)
    if (is.logical(x)) {
      y <- isTRUE(x[1])
      if (is.na(x[1])) return(default)
      return(y)
    }
    xx <- tolower(trimws(as.character(x)[1]))
    if (!nzchar(xx) || is.na(xx)) return(default)
    if (xx %in% c("t", "true", "1", "yes", "y")) return(TRUE)
    if (xx %in% c("f", "false", "0", "no", "n")) return(FALSE)
    default
  }
  .pca_plot_legend_title <- function(group_title, grp_col_final, group_col_input = NULL) {
    .normalize_group_title <- function(x) {
      x1 <- .pca_plot_as_char1(x)
      if (is.null(x1)) return(NULL)
      x_key <- toupper(trimws(x1))
      if (x_key %in% c(".GROUP_TMP", ".PROJ_GROUP")) return("Group")
      x1
    }
    gt <- .pca_plot_as_char1(group_title)
    if (!is.null(gt)) return(gt)
    gci <- .normalize_group_title(group_col_input)
    if (!is.null(gci)) return(gci)
    gc <- .normalize_group_title(grp_col_final)
    if (is.null(gc)) return("Group")
    gc
  }
  .pca_plot_projection_labels <- function(x) {
    if (is.null(x) || length(x) == 0L) return(c("Reference", "Target"))
    y <- as.character(x)
    y <- y[!is.na(y) & nzchar(trimws(y))]
    if (!length(y)) return(c("Reference", "Target"))
    if (length(y) == 1L) return(c(y[1], "Target"))
    c(y[1], y[2])
  }
  .pca_plot_as_ref_flag <- function(x) {
    if (is.null(x)) return(logical())
    if (is.logical(x)) return(as.logical(x))
    xc <- tolower(trimws(as.character(x)))
    out <- rep(NA, length(xc))
    out[xc %in% c("true", "t", "1", "yes", "y", "ref", "reference")] <- TRUE
    out[xc %in% c("false", "f", "0", "no", "n", "target")] <- FALSE
    raw_lgl <- suppressWarnings(as.logical(x))
    out[is.na(out) & !is.na(raw_lgl)] <- raw_lgl[is.na(out) & !is.na(raw_lgl)]
    as.logical(out)
  }
  .pca_plot_is_valid_color <- function(z) {
    if (is.null(z) || length(z) == 0L) return(FALSE)
    zz <- as.character(z)[1]
    if (is.na(zz) || !nzchar(zz)) return(FALSE)
    !is.null(tryCatch(grDevices::col2rgb(zz), error = function(e) NULL))
  }
  .pca_plot_lighten_hex <- function(cols, factor = 0.55) {
    nms <- names(cols)
    cols <- as.character(cols)
    out <- cols
    ok <- !is.na(cols) & nzchar(cols)
    if (!any(ok)) {
      names(out) <- nms
      return(out)
    }
    rgb <- tryCatch(grDevices::col2rgb(cols[ok]), error = function(e) NULL)
    if (is.null(rgb)) {
      names(out) <- nms
      return(out)
    }
    rgb_new <- round(rgb + (255 - rgb) * factor)
    rgb_new[rgb_new < 0] <- 0
    rgb_new[rgb_new > 255] <- 255
    out[ok] <- grDevices::rgb(rgb_new[1, ], rgb_new[2, ], rgb_new[3, ], maxColorValue = 255)
    names(out) <- nms
    out
  }
  .pca_plot_lookup_color <- function(col_map, key, default = "grey40") {
    if (is.null(col_map) || !length(col_map) || is.null(key) || !length(key)) return(default)
    key1 <- trimws(as.character(key)[1])
    nms <- names(col_map)
    if (is.null(nms) || !length(nms)) return(default)
    hit <- match(key1, nms)
    if (is.na(hit)) {
      hit <- match(tolower(key1), tolower(trimws(nms)))
    }
    if (is.na(hit)) return(default)
    val <- as.character(col_map[hit])[1]
    if (is.na(val) || !nzchar(val)) default else val
  }
  proj_labels <- .pca_plot_projection_labels(projection.group)
  proj_ref_label <- proj_labels[1]
  proj_tar_label <- proj_labels[2]
  .pca_plot_normalize_group_colors <- function(group_color, levs) {
    levs <- as.character(levs)
    n_cat <- length(levs)
    if (n_cat <= 0L) return(character())
    unknown_idx <- which(tolower(trimws(levs)) == "unknown")
    non_unknown_idx <- setdiff(seq_len(n_cat), unknown_idx)
    n_non_unknown <- length(non_unknown_idx)
    .apply_unknown_default <- function(cols, force = FALSE) {
      if (!length(unknown_idx)) return(cols)
      if (isTRUE(force)) cols[unknown_idx] <- "#4D4D4D"
      cols
    }
    .palette_to_colors <- function(mode_name, n_use = n_non_unknown) {
      if (!is.finite(n_use) || is.na(n_use) || n_use <= 0L) return(character())
      out <- tryCatch(
        get.colors(mode_name, as_int(n_use), discrete = FALSE, far = TRUE, random = FALSE, plot = FALSE, silent = TRUE),
        error = function(e) NULL
      )
      if (is.null(out) || !length(out)) return(NULL)
      as.character(out)
    }
    .fill_non_unknown <- function(out, vals) {
      if (!n_non_unknown) return(out)
      out[non_unknown_idx] <- rep_len(as.character(vals), n_non_unknown)
      out
    }
    if (is.null(group_color) || length(group_color) == 0L) {
      if (n_non_unknown <= 1L) {
        out <- stats::setNames(rep(NA_character_, n_cat), levs)
        out <- .fill_non_unknown(out, get.colors("darkrainbow", max(1L, n_non_unknown), far = TRUE, random = FALSE, plot = FALSE, silent = TRUE))
        return(.apply_unknown_default(out, force = TRUE))
      }
      out <- stats::setNames(rep(NA_character_, n_cat), levs)
      out <- .fill_non_unknown(out, get.colors("darkrainbow", n_non_unknown, far = TRUE, random = FALSE, plot = FALSE, silent = TRUE))
      return(.apply_unknown_default(out, force = TRUE))
    }
    cc <- group_color
    if (is.list(cc) && !is.data.frame(cc) && !data.table::is.data.table(cc)) cc <- unlist(cc, use.names = TRUE)
    cc <- as.character(cc)
    cc <- cc[!is.na(cc) & nzchar(cc)]
    if (!length(cc)) {
      if (n_non_unknown <= 1L) {
        out <- stats::setNames(rep(NA_character_, n_cat), levs)
        out <- .fill_non_unknown(out, get.colors("darkrainbow", max(1L, n_non_unknown), far = TRUE, random = FALSE, plot = FALSE, silent = TRUE))
        return(.apply_unknown_default(out, force = TRUE))
      }
      out <- stats::setNames(rep(NA_character_, n_cat), levs)
      out <- .fill_non_unknown(out, get.colors("darkrainbow", n_non_unknown, far = TRUE, random = FALSE, plot = FALSE, silent = TRUE))
      return(.apply_unknown_default(out, force = TRUE))
    }
    if (length(cc) == 1L && !.pca_plot_is_valid_color(cc[1])) {
      pal_cols <- .palette_to_colors(cc[1], n_use = n_non_unknown)
      if (!is.null(pal_cols) && length(pal_cols)) {
        out <- stats::setNames(rep(NA_character_, n_cat), levs)
        out <- .fill_non_unknown(out, pal_cols)
        return(.apply_unknown_default(out, force = TRUE))
      }
    }
    nms <- names(cc)
    if (!is.null(nms) && any(nzchar(nms))) {
      out <- rep(NA_character_, n_cat); names(out) <- levs
      hit <- match(levs, nms)
      ok <- !is.na(hit)
      out[ok] <- cc[hit[ok]]
      miss <- which(!ok & !(seq_len(n_cat) %in% unknown_idx))
      if (length(miss)) {
        fill <- cc[is.na(nms) | !nzchar(nms)]
        if (!length(fill)) fill <- cc
        out[miss] <- rep_len(fill, length(miss))
      }
      need <- which((is.na(out) | !nzchar(out)) & !(seq_len(n_cat) %in% unknown_idx))
      if (length(need)) {
        out[need] <- rep_len(get.colors("darkrainbow", length(need), far = TRUE, random = FALSE, plot = FALSE, silent = TRUE), length(need))
      }
      unknown_named <- any(tolower(trimws(nms)) == "unknown", na.rm = TRUE)
      if (!isTRUE(unknown_named)) out <- .apply_unknown_default(out, force = TRUE)
      return(out)
    }
    out <- stats::setNames(rep(NA_character_, n_cat), levs)
    out <- .fill_non_unknown(out, cc)
    .apply_unknown_default(out, force = TRUE)
  }
  .pca_plot_normalize_group_shapes <- function(shape_spec, levs, default_shape = 16) {
    levs <- as.character(levs)
    n_cat <- length(levs)
    if (n_cat <= 0L) return(numeric())
    default_shape_num <- suppressWarnings(as.numeric(unlist(default_shape, use.names = FALSE))[1])
    if (!is.finite(default_shape_num) || is.na(default_shape_num)) default_shape_num <- 16
    if (is.null(shape_spec) || length(shape_spec) == 0L) {
      out <- rep(default_shape_num, n_cat)
      names(out) <- levs
      return(out)
    }
    sh <- shape_spec
    if (is.list(sh) && !is.data.frame(sh) && !data.table::is.data.table(sh)) sh <- unlist(sh, use.names = TRUE)
    sh_num <- suppressWarnings(as.numeric(sh))
    nms <- names(sh)
    if (!is.null(nms) && any(nzchar(nms))) {
      out <- rep(NA_real_, n_cat); names(out) <- levs
      hit <- match(levs, nms)
      ok <- !is.na(hit)
      out[ok] <- sh_num[hit[ok]]
      miss <- which(!ok | !is.finite(out))
      if (length(miss)) {
        fill <- sh_num[is.na(nms) | !nzchar(nms)]
        if (!length(fill)) fill <- sh_num
        fill <- fill[is.finite(fill)]
        if (!length(fill)) fill <- default_shape_num
        out[miss] <- rep_len(fill, length(miss))
      }
      out[!is.finite(out) | is.na(out)] <- default_shape_num
      return(out)
    }
    sh_num <- sh_num[is.finite(sh_num)]
    if (!length(sh_num)) sh_num <- default_shape_num
    out <- rep_len(sh_num, n_cat)
    names(out) <- levs
    out
  }
  .pca_plot_normalize_group_sizes <- function(size_spec, levs, default_size = 1.2) {
    levs <- as.character(levs)
    n_cat <- length(levs)
    if (n_cat <= 0L) return(numeric())
    if (is.null(size_spec) || length(size_spec) == 0L) {
      out <- rep(as.numeric(default_size), n_cat)
      names(out) <- levs
      return(out)
    }
    sz <- size_spec
    if (is.list(sz) && !is.data.frame(sz) && !data.table::is.data.table(sz)) sz <- unlist(sz, use.names = TRUE)
    sz_num <- suppressWarnings(as.numeric(sz))
    nms <- names(sz)
    if (!is.null(nms) && any(nzchar(nms))) {
      out <- rep(NA_real_, n_cat); names(out) <- levs
      hit <- match(levs, nms)
      ok <- !is.na(hit)
      out[ok] <- sz_num[hit[ok]]
      miss <- which(!ok | !is.finite(out) | out <= 0)
      if (length(miss)) {
        fill <- sz_num[is.na(nms) | !nzchar(nms)]
        if (!length(fill)) fill <- sz_num
        fill <- fill[is.finite(fill) & fill > 0]
        if (!length(fill)) fill <- as.numeric(default_size)
        out[miss] <- rep_len(fill, length(miss))
      }
      out[!is.finite(out) | is.na(out) | out <= 0] <- as.numeric(default_size)
      return(out)
    }
    sz_num <- sz_num[is.finite(sz_num) & sz_num > 0]
    if (!length(sz_num)) sz_num <- as.numeric(default_size)
    out <- rep_len(sz_num, n_cat)
    names(out) <- levs
    out
  }
  .pca_plot_normalize_group_strokes <- function(stroke_spec, levs, default_stroke = 1) {
    levs <- as.character(levs)
    n_cat <- length(levs)
    if (n_cat <= 0L) return(numeric())
    if (is.null(stroke_spec) || length(stroke_spec) == 0L) {
      out <- rep(as.numeric(default_stroke), n_cat)
      names(out) <- levs
      return(out)
    }
    st <- stroke_spec
    if (is.list(st) && !is.data.frame(st) && !data.table::is.data.table(st)) st <- unlist(st, use.names = TRUE)
    st_num <- suppressWarnings(as.numeric(st))
    nms <- names(st)
    if (!is.null(nms) && any(nzchar(nms))) {
      out <- rep(NA_real_, n_cat); names(out) <- levs
      hit <- match(levs, nms)
      ok <- !is.na(hit)
      out[ok] <- st_num[hit[ok]]
      miss <- which(!ok | !is.finite(out) | out < 0)
      if (length(miss)) {
        fill <- st_num[is.na(nms) | !nzchar(nms)]
        if (!length(fill)) fill <- st_num
        fill <- fill[is.finite(fill) & fill >= 0]
        if (!length(fill)) fill <- as.numeric(default_stroke)
        out[miss] <- rep_len(fill, length(miss))
      }
      out[!is.finite(out) | is.na(out) | out < 0] <- as.numeric(default_stroke)
      return(out)
    }
    st_num <- st_num[is.finite(st_num) & st_num >= 0]
    if (!length(st_num)) st_num <- as.numeric(default_stroke)
    out <- rep_len(st_num, n_cat)
    names(out) <- levs
    out
  }
  .pca_plot_parse_line_sd <- function(x) {
    if (is.null(x) || length(x) == 0L) return(numeric())
    if (is.logical(x) && length(x) == 1L) {
      if (isTRUE(x)) return(c(5, 3))
      return(numeric())
    }
    if (is.character(x)) {
      xs <- tolower(trimws(as.character(x)))
      xs <- xs[!is.na(xs) & nzchar(xs) & xs != "false"]
      if (!length(xs)) return(numeric())
      xv <- suppressWarnings(as.numeric(gsub("sd$", "", xs)))
      xv <- xv[is.finite(xv) & !is.na(xv) & xv > 0]
      return(unique(xv))
    }
    xv <- suppressWarnings(as.numeric(x))
    xv <- xv[is.finite(xv) & !is.na(xv) & xv > 0]
    unique(xv)
  }
  .pca_plot_pc_cols <- function(nm) {
    nm <- as.character(nm)
    hit <- grep("^PC[0-9]+$", nm, ignore.case = FALSE, value = TRUE)
    hit[order(as_int(sub("^PC", "", hit)))]
  }
  .pca_plot_make_eval_named <- function(eval_in, pc_cols) {
    if (is.null(eval_in) || length(eval_in) == 0L) return(NULL)
    ev <- NULL
    nm <- NULL
    .normalize_pc_names <- function(x) {
      if (is.null(x) || !length(x)) return(x)
      x <- as.character(x)
      x <- trimws(x)
      x[x %in% c("", "NA")] <- NA_character_
      is_num <- !is.na(x) & grepl("^[0-9]+$", x)
      x[is_num] <- paste0("PC", x[is_num])
      has_pc <- !is.na(x) & grepl("^PC[0-9]+$", toupper(x))
      x[has_pc] <- paste0("PC", sub("^PC", "", toupper(x[has_pc])))
      x
    }
    if (data.table::is.data.table(eval_in) || is.data.frame(eval_in)) {
      dt_ev <- if (data.table::is.data.table(eval_in)) data.table::copy(eval_in) else data.table::as.data.table(eval_in)
      cn_ev <- names(dt_ev)
      if (ncol(dt_ev) == 1L) {
        ev_col <- cn_ev[1]
      } else {
        key_hits <- match(c("EIGENVAL", "EVAL", "VALUE"), toupper(cn_ev))
        key_hits <- key_hits[!is.na(key_hits)]
        if (length(key_hits)) {
          ev_col <- cn_ev[key_hits[1]]
        } else {
          num_cols <- cn_ev[vapply(dt_ev, function(v) {
            x <- suppressWarnings(as.numeric(v))
            any(is.finite(x) & !is.na(x))
          }, logical(1))]
          ev_col <- if (length(num_cols)) num_cols[1] else cn_ev[1]
        }
      }
      ev <- suppressWarnings(as.numeric(dt_ev[[ev_col]]))
      keep <- is.finite(ev) & !is.na(ev)
      ev <- ev[keep]
      if (!length(ev)) return(NULL)
      pc_col <- cn_ev[toupper(cn_ev) %in% c("PC", "COMPONENT", "INDEX")]
      if (length(pc_col)) nm <- .normalize_pc_names(dt_ev[[pc_col[1]]])[keep]
    } else if (is.character(eval_in) && length(eval_in) == 1L) {
      pth <- abs_path(eval_in)[1]
      if (!is.na(pth) && nzchar(pth) && file.exists(pth)) {
        dt_ev <- data.table::fread(pth, data.table = TRUE, check.names = FALSE, showProgress = FALSE)
        cn_ev <- names(dt_ev)
        if (ncol(dt_ev) == 1L) {
          ev_col <- cn_ev[1]
        } else {
          key_hits <- match(c("EIGENVAL", "EVAL", "VALUE"), toupper(cn_ev))
          key_hits <- key_hits[!is.na(key_hits)]
          if (length(key_hits)) {
            ev_col <- cn_ev[key_hits[1]]
          } else {
            num_cols <- cn_ev[vapply(dt_ev, function(v) {
              x <- suppressWarnings(as.numeric(v))
              any(is.finite(x) & !is.na(x))
            }, logical(1))]
            ev_col <- if (length(num_cols)) num_cols[1] else cn_ev[1]
          }
        }
        ev <- suppressWarnings(as.numeric(dt_ev[[ev_col]]))
        keep <- is.finite(ev) & !is.na(ev)
        ev <- ev[keep]
        if (!length(ev)) return(NULL)
        pc_col <- cn_ev[toupper(cn_ev) %in% c("PC", "COMPONENT", "INDEX")]
        if (length(pc_col)) nm <- .normalize_pc_names(dt_ev[[pc_col[1]]])[keep]
      } else {
        ev <- suppressWarnings(as.numeric(eval_in))
      }
    } else {
      ev <- suppressWarnings(as.numeric(eval_in))
    }
    if (is.null(ev) || !length(ev)) return(NULL)
    if (is.null(nm)) nm <- .normalize_pc_names(names(eval_in))
    if (is.null(nm) || !length(nm) || all(is.na(nm) | !nzchar(nm))) {
      nm <- pc_cols[seq_len(min(length(pc_cols), length(ev)))]
      if (length(nm) < length(ev)) nm <- c(nm, paste0("PC", seq_len(length(ev) - length(nm)) + length(nm)))
    }
    nm <- as.character(nm)
    nm[is.na(nm) | !nzchar(nm)] <- paste0("PC", seq_len(sum(is.na(nm) | !nzchar(nm))))
    out <- ev
    names(out) <- nm
    out
  }
  .pca_plot_axis_label <- function(k, eval_named) {
    nm <- paste0("PC", as_int(k))
    if (is.null(eval_named) || !(nm %in% names(eval_named))) return(nm)
    ev <- suppressWarnings(as.numeric(eval_named[nm]))
    s <- sum(as.numeric(eval_named), na.rm = TRUE)
    if (!is.finite(ev) || is.na(ev) || !is.finite(s) || is.na(s) || s <= 0) return(nm)
    sprintf("%s (%.2f%%)", nm, 100 * ev / s)
  }
  .pca_plot_as_pc_dt <- function(x) {
    if (is.null(x)) return(NULL)
    if (data.table::is.data.table(x)) return(data.table::copy(x))
    if (is.data.frame(x)) return(data.table::as.data.table(x))
    if (is.matrix(x)) {
      xx <- as.data.frame(x, stringsAsFactors = FALSE)
      return(data.table::as.data.table(xx))
    }
    NULL
  }
  .pca_plot_extract_eval <- function(x) {
    if (is.null(x)) return(NULL)
    if (data.table::is.data.table(x) || is.data.frame(x) || is.matrix(x)) return(NULL)
    if (is.list(x)) {
      if (!is.null(x$eigenvalue) && length(x$eigenvalue)) return(x$eigenvalue)
      if (!is.null(x$prop.var) && length(x$prop.var)) return(x$prop.var)
      if (!is.null(x$eval) && length(x$eval)) return(x$eval)
    }
    if (!is.null(attr(x, "eigenvalue")) && length(attr(x, "eigenvalue"))) return(attr(x, "eigenvalue"))
    if (!is.null(attr(x, "prop.var")) && length(attr(x, "prop.var"))) return(attr(x, "prop.var"))
    NULL
  }
  .pca_plot_groupdata_from_list <- function(x) {
    if (is.null(x) || !is.list(x) || is.data.frame(x) || data.table::is.data.table(x)) return(NULL)
    if (!length(x)) return(data.table::data.table(IID = character(), .GROUP_LIST = character()))
    nms <- names(x)
    if (is.null(nms)) nms <- rep("", length(x))
    parts <- vector("list", length(x))
    lev <- character()
    for (i in seq_along(x)) {
      gi <- as.character(nms[i])
      if (is.na(gi) || !nzchar(gi)) gi <- as.character(i)
      xi <- x[[i]]
      if (is.null(xi) || length(xi) == 0L) next
      iid <- as.character(xi)
      iid <- iid[!is.na(iid) & nzchar(iid)]
      if (!length(iid)) next
      lev <- c(lev, gi)
      parts[[i]] <- data.table::data.table(
        IID = unique(iid),
        .GROUP_LIST = gi
      )
    }
    parts <- Filter(Negate(is.null), parts)
    out <- if (length(parts)) data.table::rbindlist(parts, use.names = TRUE, fill = TRUE) else data.table::data.table(IID = character(), .GROUP_LIST = character())
    out <- unique(out, by = c("IID", ".GROUP_LIST"))
    lev <- lev[!is.na(lev) & nzchar(lev)]
    attr(out, "group_levels") <- sort(unique(as.character(lev)))
    out
  }
  .pca_plot_extract_outlier <- function(x, id_col_pref = "IID") {
    out <- list(
      mode = "none",
      drop = character(),
      ambiguous = character(),
      group.data = NULL,
      group.col = NULL,
      group.title = NULL
    )
    if (is.null(x) || length(x) == 0L) return(out)
    if (inherits(x, "gcanvas_pca_refqc")) {
      out$mode <- "refqc"
      out$drop <- as.character(x$id.drop)
      out$ambiguous <- as.character(x$id.ambiguous)
      out$group.title <- .pca_plot_as_char1(x$group.title)
      if (!is.null(x$posterior) && (data.table::is.data.table(x$posterior) || is.data.frame(x$posterior))) {
        x <- x$posterior
      } else {
        return(out)
      }
    }
    if (data.table::is.data.table(x) || is.data.frame(x)) {
      dt0 <- if (data.table::is.data.table(x)) data.table::copy(x) else data.table::as.data.table(x)
      id0 <- .pca_ancestry_pick_id_col(names(dt0), preferred = id_col_pref)
      if (is.null(id0) || !(id0 %in% names(dt0))) return(out)
      dt0[, .ID := as.character(get(id0))]
      dt0 <- dt0[!is.na(.ID) & nzchar(.ID)]
      if (!nrow(dt0)) return(out)
      has_keep <- "keep" %in% names(dt0)
      has_amb <- "ambiguous" %in% names(dt0)
      has_reason <- "REASON" %in% names(dt0)
      has_label_input <- "label.input" %in% names(dt0)
      has_label_max <- "label.max" %in% names(dt0)
      if (has_keep || has_amb || has_reason) {
        out$mode <- "refqc"
        if (has_keep) out$drop <- as.character(dt0[get("keep") %in% FALSE, .ID])
        if (has_amb) out$ambiguous <- as.character(dt0[get("ambiguous") %in% TRUE & (!has_keep | get("keep") %in% TRUE), .ID])
        grp_col0 <- NULL
        if (has_label_input) grp_col0 <- "label.input" else if (has_label_max) grp_col0 <- "label.max"
        if (!is.null(grp_col0)) {
          gd0 <- dt0[, .(IID = .ID, .OUTLIER_GROUP = as.character(get(grp_col0)))]
          gd0 <- gd0[!is.na(.OUTLIER_GROUP) & nzchar(trimws(.OUTLIER_GROUP))]
          gd0 <- unique(gd0, by = c("IID", ".OUTLIER_GROUP"))
          if (nrow(gd0)) {
            out$group.data <- gd0
            out$group.col <- ".OUTLIER_GROUP"
          }
        }
        return(out)
      }
      out$mode <- "simple"
      out$drop <- unique(as.character(dt0$.ID))
      return(out)
    }
    if (is.atomic(x) && !is.list(x)) {
      ids <- as.character(x)
      ids <- unique(ids[!is.na(ids) & nzchar(ids)])
      if (!length(ids)) return(out)
      out$mode <- "simple"
      out$drop <- ids
      return(out)
    }
    if (is.list(x) && !is.data.frame(x) && !data.table::is.data.table(x)) {
      ids <- unlist(x, recursive = TRUE, use.names = FALSE)
      ids <- unique(as.character(ids))
      ids <- ids[!is.na(ids) & nzchar(ids)]
      if (!length(ids)) return(out)
      out$mode <- "simple"
      out$drop <- ids
      return(out)
    }
    out
  }
  .pca_plot_extract_ambiguous <- function(x, id_col_pref = "IID") {
    ids <- character()
    if (is.null(x) || length(x) == 0L) return(ids)
    if (inherits(x, "gcanvas_pca_refqc")) {
      ids <- as.character(x$id.ambiguous)
      ids <- ids[!is.na(ids) & nzchar(ids)]
      if (length(ids)) return(unique(ids))
      if (!is.null(x$posterior) && (data.table::is.data.table(x$posterior) || is.data.frame(x$posterior))) {
        x <- x$posterior
      } else {
        return(ids)
      }
    }
    if (data.table::is.data.table(x) || is.data.frame(x)) {
      dt0 <- if (data.table::is.data.table(x)) data.table::copy(x) else data.table::as.data.table(x)
      id0 <- .pca_ancestry_pick_id_col(names(dt0), preferred = id_col_pref)
      if (is.null(id0) || !(id0 %in% names(dt0))) return(ids)
      dt0[, .ID := as.character(get(id0))]
      dt0 <- dt0[!is.na(.ID) & nzchar(.ID)]
      if (!nrow(dt0)) return(ids)
      if ("ambiguous" %in% names(dt0)) {
        ids <- as.character(dt0[get("ambiguous") %in% TRUE, .ID])
        ids <- ids[!is.na(ids) & nzchar(ids)]
        return(unique(ids))
      }
      return(unique(as.character(dt0$.ID)))
    }
    if (is.atomic(x) && !is.list(x)) {
      ids <- unique(as.character(x))
      ids <- ids[!is.na(ids) & nzchar(ids)]
      return(ids)
    }
    if (is.list(x) && !is.data.frame(x) && !data.table::is.data.table(x)) {
      ids <- unlist(x, recursive = TRUE, use.names = FALSE)
      ids <- unique(as.character(ids))
      ids <- ids[!is.na(ids) & nzchar(ids)]
      return(ids)
    }
    ids
  }

  dt <- NULL
  group_levels <- NULL
  point.size <- suppressWarnings(as.numeric(point.size))[1]
  if (!is.finite(point.size) || is.na(point.size) || point.size <= 0) point.size <- 1.2
  point.stroke <- suppressWarnings(as.numeric(point.stroke))[1]
  if (!is.finite(point.stroke) || is.na(point.stroke) || point.stroke < 0) point.stroke <- 1
  alpha <- suppressWarnings(as.numeric(alpha))[1]
  if (!is.finite(alpha) || is.na(alpha)) alpha <- 0.65
  alpha <- max(0, min(1, alpha))
  eval_use <- eval
  pc_in <- pc
  if (inherits(pc_in, "gcanvas_pca_refqc")) {
    refqc.helper <- TRUE
    eval_refqc <- .pca_plot_extract_eval(pc_in)
    refqc_drop_ids <- as.character(pc_in$id.drop)
    refqc_drop_ids <- refqc_drop_ids[!is.na(refqc_drop_ids) & nzchar(refqc_drop_ids)]
    refqc_ambiguous_ids <- as.character(pc_in$id.ambiguous)
    refqc_ambiguous_ids <- refqc_ambiguous_ids[!is.na(refqc_ambiguous_ids) & nzchar(refqc_ambiguous_ids)]
    refqc_group_title <- .pca_plot_as_char1(pc_in$group.title)
    refqc_out_info <- NULL
    if (!is.null(pc_in$posterior) && (data.table::is.data.table(pc_in$posterior) || is.data.frame(pc_in$posterior))) {
      refqc_out_info <- .pca_plot_extract_outlier(pc_in$posterior, id_col_pref = id.col)
    }
    if (!is.null(pc_in$pc)) pc_in <- pc_in$pc
    if ((is.null(eval_use) || length(eval_use) == 0L) && !is.null(eval_refqc) && length(eval_refqc)) eval_use <- eval_refqc
    if (is.null(group.data) && (is.null(group.col) || !nzchar(.pca_plot_as_char1(group.col)))) {
      if (!is.null(refqc_out_info$group.data)) {
        group.data <- refqc_out_info$group.data
        group.col <- refqc_out_info$group.col %||% ".OUTLIER_GROUP"
      } else {
        group.col <- ".REFQC_GROUP"
      }
    }
    if (is.null(group.title) || !nzchar(.pca_plot_as_char1(group.title))) {
      group.title <- refqc_group_title
    }
  }
  if (!isTRUE(refqc.helper) && !is.null(outlier)) {
    out_info <- .pca_plot_extract_outlier(outlier, id_col_pref = id.col)
    if (identical(out_info$mode, "refqc")) {
      refqc.helper <- TRUE
      refqc_drop_ids <- unique(as.character(out_info$drop))
      refqc_drop_ids <- refqc_drop_ids[!is.na(refqc_drop_ids) & nzchar(refqc_drop_ids)]
      refqc_ambiguous_ids <- unique(as.character(out_info$ambiguous))
      refqc_ambiguous_ids <- refqc_ambiguous_ids[!is.na(refqc_ambiguous_ids) & nzchar(refqc_ambiguous_ids)]
      if ((is.null(group.data) || length(group.data) == 0L) &&
          (is.null(group.col) || !nzchar(.pca_plot_as_char1(group.col))) &&
          !is.null(out_info$group.data)) {
        group.data <- out_info$group.data
        group.col <- out_info$group.col %||% ".OUTLIER_GROUP"
      }
      if (is.null(group.title) || !nzchar(.pca_plot_as_char1(group.title))) {
        if (!is.null(out_info$group.title) && nzchar(.pca_plot_as_char1(out_info$group.title))) {
          group.title <- out_info$group.title
        }
      }
    } else if (identical(out_info$mode, "simple")) {
      refqc.helper <- TRUE
      refqc_drop_ids <- unique(as.character(out_info$drop))
      refqc_drop_ids <- refqc_drop_ids[!is.na(refqc_drop_ids) & nzchar(refqc_drop_ids)]
      refqc_ambiguous_ids <- character()
      refqc_drop_label <- "Outlier"
    }
  }
  if (!is.null(ambiguous)) {
    amb_ids <- .pca_plot_extract_ambiguous(ambiguous, id_col_pref = id.col)
    if (length(amb_ids)) {
      refqc.helper <- TRUE
      refqc_ambiguous_ids <- unique(c(refqc_ambiguous_ids, amb_ids))
      refqc_ambiguous_ids <- refqc_ambiguous_ids[!is.na(refqc_ambiguous_ids) & nzchar(refqc_ambiguous_ids)]
    }
  }
  fade_on <- .pca_plot_as_flag1(fade, default = NA)
  if (is.na(fade_on)) {
    fade_on <- isTRUE(refqc.helper) || length(refqc_drop_ids) > 0L || length(refqc_ambiguous_ids) > 0L
  }
  fade_on <- isTRUE(fade_on)
  if (is.list(pc_in) && !data.table::is.data.table(pc_in) && !is.data.frame(pc_in)) {
    # Accept PCA return objects and generic list containers.
    cand_keys <- c("pc", "PC", "data", "scores", "eigenvec")
    for (kk in cand_keys) {
      if (!is.null(pc_in[[kk]])) {
        dt0 <- .pca_plot_as_pc_dt(pc_in[[kk]])
        if (!is.null(dt0)) {
          dt <- dt0
          break
        }
      }
    }
    if (is.null(dt) && length(pc_in)) {
      for (ii in seq_along(pc_in)) {
        dt0 <- .pca_plot_as_pc_dt(pc_in[[ii]])
        if (is.null(dt0)) next
        pcc0 <- grep("^PC[0-9]+$", names(dt0), ignore.case = FALSE, value = TRUE)
        if (length(pcc0) >= 2L) {
          dt <- dt0
          break
        }
      }
    }
    if (is.null(dt)) stop("pc must be pca() return object or data.frame/data.table with PC columns.", call. = FALSE)
    if (is.null(eval_use) || length(eval_use) == 0L) {
      eval_use <- .pca_plot_extract_eval(pc_in)
    }
  } else if (is.data.frame(pc_in) || data.table::is.data.table(pc_in)) {
    dt <- if (data.table::is.data.table(pc_in)) data.table::copy(pc_in) else data.table::as.data.table(pc_in)
  } else if (is.matrix(pc_in)) {
    dt <- data.table::as.data.table(as.data.frame(pc_in, stringsAsFactors = FALSE))
  } else {
    stop("pc must be pca() return object or data.frame/data.table with PC columns.", call. = FALSE)
  }
  if (is.null(dt) || !nrow(dt)) stop("No rows in pc.", call. = FALSE)

  id_col_use <- .pca_plot_as_char1(id.col)
  if (is.null(id_col_use) || !(id_col_use %in% names(dt))) {
    id_candidates <- intersect(c("IID", "#IID", "ID", "sample.ID", "FID"), names(dt))
    if (length(id_candidates)) id_col_use <- id_candidates[1]
  }
  if (is.null(id_col_use) || !(id_col_use %in% names(dt))) stop("id.col not found in pc.", call. = FALSE)
  dt[, IID_USE := as.character(get(id_col_use))]
  dt[is.na(IID_USE) | !nzchar(IID_USE), IID_USE := paste0("row", .I)]

  pcc <- .pca_plot_pc_cols(names(dt))
  if (length(pcc) < 2L) stop("At least two PC columns (PC1, PC2, ...) are required.", call. = FALSE)
  for (nm in pcc) dt[, (nm) := suppressWarnings(as.numeric(get(nm)))]

  .pca_plot_parse_pc_use <- function(x, pc_cols) {
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
  pc_idx <- .pca_plot_parse_pc_use(pc.use, pcc)
  pc_idx <- unique(pc_idx)
  if (length(pc_idx) < 2L) stop("pc.use must contain at least two PCs.", call. = FALSE)
  pair_list <- lapply(seq_len(length(pc_idx) - 1L), function(i) c(pc_idx[i], pc_idx[i + 1L]))
  miss_pc <- unique(unlist(lapply(pair_list, function(v) paste0("PC", v))))
  miss_pc <- miss_pc[!(miss_pc %in% pcc)]
  if (length(miss_pc)) stop("Missing PC columns in data: ", paste(miss_pc, collapse = ", "), call. = FALSE)

  eval_named <- .pca_plot_make_eval_named(eval_use, pcc)

  # Projection helper: auto-map reference/target without extra group settings.
  if (isTRUE(projection.helper) &&
      is.null(group.data) &&
      (is.null(group.col) || !nzchar(as.character(group.col)[1])) &&
      ("reference" %in% names(dt))) {
    ref_flag <- .pca_plot_as_ref_flag(dt$reference)
    dt[, .PROJ_GROUP := data.table::fifelse(ref_flag %in% TRUE, "Reference", "Target")]
    group.col <- ".PROJ_GROUP"
    if (is.null(group.color) || length(group.color) == 0L) {
      group.color <- c(Reference = "#4D4D4D", Target = "#D7263D")
    }
    if (is.null(group.shape) || length(group.shape) == 0L) {
      group.shape <- c(Reference = 16, Target = 16)
    }
  }

  # Group handling
  grp_col_final <- NULL
  if (!is.null(group.data)) {
    if (is.list(group.data) && !is.data.frame(group.data) && !data.table::is.data.table(group.data)) {
      gd <- .pca_plot_groupdata_from_list(group.data)
      if (is.null(gd)) stop("group.data must be vector/data.frame/data.table/list.", call. = FALSE)
      group_levels <- attr(gd, "group_levels", exact = TRUE)
      data.table::setnames(gd, c("IID", ".GROUP_LIST"), c("IDM", "GRP"))
      data.table::setkey(gd, IDM)
      dt <- merge(dt, gd, by.x = "IID_USE", by.y = "IDM", all.x = TRUE, sort = FALSE)
      data.table::setnames(dt, "GRP", ".GROUP_TMP")
      grp_col_final <- ".GROUP_TMP"
      if (is.null(group.col) || !nzchar(as.character(group.col)[1])) group.col <- grp_col_final
    } else if (is.atomic(group.data) && !is.list(group.data) && !is.data.frame(group.data) && !data.table::is.data.table(group.data)) {
      if (length(group.data) == nrow(dt)) {
        if (is.factor(group.data)) {
          group_levels <- levels(group.data)
        }
        dt[, .GROUP_TMP := as.character(group.data)]
        grp_col_final <- ".GROUP_TMP"
        if (is.null(group.col) || !nzchar(as.character(group.col)[1])) group.col <- grp_col_final
      } else {
        iid_sel <- as.character(group.data)
        iid_sel <- iid_sel[!is.na(iid_sel) & nzchar(iid_sel)]
        if (!length(iid_sel)) stop("group.data IID vector has no valid IDs.", call. = FALSE)
        iid_sel <- unique(iid_sel)
        dt[, .GROUP_TMP := data.table::fifelse(IID_USE %in% iid_sel, "Target", "Others")]
        group_levels <- c("Target", "Others")
        grp_col_final <- ".GROUP_TMP"
        if (is.null(group.col) || !nzchar(as.character(group.col)[1])) group.col <- grp_col_final
        if (is.null(group.color) || length(group.color) == 0L) {
          group.color <- c(Target = "#D7263D", Others = "#4D4D4D")
        }
      }
    } else {
      gd <- if (data.table::is.data.table(group.data)) data.table::copy(group.data) else if (is.data.frame(group.data)) data.table::as.data.table(group.data) else NULL
      if (is.null(gd)) stop("group.data must be vector/data.frame/data.table.", call. = FALSE)
      gcol <- .pca_plot_as_char1(group.col)
      if (is.null(gcol) && ncol(gd) == 2L) {
        id_col_gd <- names(gd)[1]
        grp_col_gd <- names(gd)[2]
      } else {
        id_col_gd <- .pca_ancestry_pick_id_col(names(gd), preferred = id_col_use)
        if (is.null(id_col_gd) || !(id_col_gd %in% names(gd))) stop("id.col not found in group.data.", call. = FALSE)
        if (is.null(gcol)) {
          cands <- setdiff(names(gd), id_col_gd)
          if (!length(cands)) stop("group.col not found in group.data.", call. = FALSE)
          grp_col_gd <- cands[1]
        } else if (!(gcol %in% names(gd))) {
          stop("group.col not found in group.data: ", gcol, call. = FALSE)
        } else {
          grp_col_gd <- gcol
        }
      }
      if (is.factor(gd[[grp_col_gd]])) {
        group_levels <- levels(gd[[grp_col_gd]])
      }
      ref_col_gd <- NULL
      ref_hit <- names(gd)[tolower(names(gd)) %in% c("reference", "ref", "is_reference", "is.reference", "isref")]
      if (length(ref_hit)) ref_col_gd <- ref_hit[1]
      if (!is.null(ref_col_gd)) {
        ref_vec <- .pca_plot_as_ref_flag(gd[[ref_col_gd]])
        gd <- gd[, .(
          IDM = as.character(get(id_col_gd)),
          GRP = as.character(get(grp_col_gd)),
          REF = ref_vec
        )]
      } else {
        gd <- gd[, .(IDM = as.character(get(id_col_gd)), GRP = as.character(get(grp_col_gd)))]
      }
      gd <- gd[!is.na(IDM) & nzchar(IDM)]
      data.table::setkey(gd, IDM)
      dt <- merge(dt, gd, by.x = "IID_USE", by.y = "IDM", all.x = TRUE, sort = FALSE)
      data.table::setnames(dt, "GRP", ".GROUP_TMP")
      if ("REF" %in% names(dt)) {
        if (!("reference" %in% names(dt))) {
          data.table::setnames(dt, "REF", "reference")
        } else {
          dt[is.na(reference), reference := REF]
          dt[, REF := NULL]
        }
      }
      grp_col_final <- ".GROUP_TMP"
    }
  } else {
    gcol <- .pca_plot_as_char1(group.col)
    if (!is.null(gcol)) {
      if (!(gcol %in% names(dt))) stop("group.col not found in data: ", gcol, call. = FALSE)
      grp_col_final <- gcol
    }
  }
  has_group <- !is.null(grp_col_final) && (grp_col_final %in% names(dt))
  if (has_group) {
    dt[, Group := as.character(get(grp_col_final))]
    dt[is.na(Group) | !nzchar(Group), Group := "Unknown"]
    lev_raw <- unique(as.character(dt$Group))
    lev_raw <- lev_raw[!is.na(lev_raw) & nzchar(lev_raw)]
    if (isTRUE(projection.helper) && identical(grp_col_final, ".PROJ_GROUP")) {
      lev <- c("Reference", "Target")
      lev <- lev[lev %in% lev_raw]
      lev <- c(lev, setdiff(sort(lev_raw), lev))
    } else if (!is.null(group_levels) && length(group_levels)) {
      lev0 <- as.character(group_levels)
      lev0 <- lev0[!is.na(lev0) & nzchar(lev0)]
      lev <- c(lev0[lev0 %in% lev_raw], setdiff(sort(lev_raw), lev0))
    } else {
      lev <- sort(lev_raw)
    }
    if (!length(lev)) lev <- "Unknown"
    dt[, Group := factor(Group, levels = lev)]
  }

  # Outliers across all requested adjacent PC pairs
  out3 <- character()
  out5 <- character()
  for (pr in pair_list) {
    xnm <- paste0("PC", pr[1]); ynm <- paste0("PC", pr[2])
    x <- suppressWarnings(as.numeric(dt[[xnm]]))
    y <- suppressWarnings(as.numeric(dt[[ynm]]))
    ok <- is.finite(x) & !is.na(x) & is.finite(y) & !is.na(y)
    if (!any(ok)) next
    mx <- mean(x[ok], na.rm = TRUE); sx <- stats::sd(x[ok], na.rm = TRUE)
    my <- mean(y[ok], na.rm = TRUE); sy <- stats::sd(y[ok], na.rm = TRUE)
    if (!is.finite(sx) || is.na(sx)) sx <- 0
    if (!is.finite(sy) || is.na(sy)) sy <- 0
    o3 <- ok & ((x < (mx - 3 * sx)) | (x > (mx + 3 * sx)) | (y < (my - 3 * sy)) | (y > (my + 3 * sy)))
    o5 <- ok & ((x < (mx - 5 * sx)) | (x > (mx + 5 * sx)) | (y < (my - 5 * sy)) | (y > (my + 5 * sy)))
    out3 <- unique(c(out3, as.character(dt$IID_USE[o3])))
    out5 <- unique(c(out5, as.character(dt$IID_USE[o5])))
  }

  # line.sd options
  line_sd <- .pca_plot_parse_line_sd(line.sd)
  if (length(line_sd)) {
    line_sd <- sort(unique(line_sd), decreasing = TRUE)
  }
  n_line <- length(line_sd)
  line_col <- as.character(line.sd.color)
  line_col <- line_col[!is.na(line_col) & nzchar(line_col)]
  if (!length(line_col)) {
    if (n_line <= 1L) {
      line_col <- "#EF476F"
    } else {
      # Higher SD threshold -> red, lower threshold -> blue.
      line_col <- grDevices::colorRampPalette(c("#EF476F", "#258AB2"))(n_line)
    }
  } else {
    line_col <- rep(line_col, length.out = max(1L, n_line))
  }
  line_typ <- as.character(line.sd.type)
  line_typ <- line_typ[!is.na(line_typ) & nzchar(line_typ)]
  if (!length(line_typ)) line_typ <- "dashed"
  line_typ <- rep(line_typ, length.out = max(1L, n_line))
  line_lwd <- .pca_plot_numvec(line.sd.linewidth)
  line_lwd <- line_lwd[is.finite(line_lwd) & !is.na(line_lwd) & line_lwd > 0]
  if (!length(line_lwd)) line_lwd <- 0.2
  line_lwd <- rep(line_lwd, length.out = max(1L, n_line))

  plot_list <- list()
  if (has_group) {
    group_legend_title <- .pca_plot_legend_title(
      if (!is.null(legend.title)) legend.title else group.title,
      grp_col_final,
      group_col_input = group.col
    )
    if (isTRUE(refqc.helper) && !is.null(refqc_group_title) && !nzchar(.pca_plot_as_char1(group.title))) {
      group_legend_title <- refqc_group_title
    }
    lev <- levels(dt$Group)
    col_spec <- if (!is.null(group.color) && length(group.color)) {
      group.color
    } else if (!isTRUE(point_color_missing)) {
      point.color
    } else {
      NULL
    }
    grp_cols <- .pca_plot_normalize_group_colors(col_spec, lev)
    grp_cols_overlay <- grp_cols
    grp_cols_ambiguous <- grp_cols
    if (isTRUE(fade_on)) {
      grp_cols <- .pca_plot_lighten_hex(grp_cols)
      grp_cols_ambiguous <- .pca_plot_lighten_hex(grp_cols_overlay, factor = 0.35)
    }
    shp_spec <- if (!is.null(group.shape) && length(group.shape)) group.shape else point.shape
    grp_shapes <- .pca_plot_normalize_group_shapes(shp_spec, lev, default_shape = point.shape)
    siz_spec <- if (!is.null(group.size) && length(group.size)) group.size else point.size
    grp_sizes <- .pca_plot_normalize_group_sizes(siz_spec, lev, default_size = point.size)
    str_spec <- if (!is.null(group.stroke) && length(group.stroke)) group.stroke else point.stroke
    grp_strokes <- .pca_plot_normalize_group_strokes(str_spec, lev, default_stroke = point.stroke)
    has_proj_reference <- isTRUE(projection.helper) && ("reference" %in% names(dt))
    is_proj_helper_group <- has_proj_reference && identical(grp_col_final, ".PROJ_GROUP")
    use_proj_style_override <- isTRUE(has_proj_reference)
    if (isTRUE(use_proj_style_override) && isTRUE(is_proj_helper_group) && ("Target" %in% names(grp_shapes))) {
      grp_shapes["Target"] <- 15
    }
    grp_alpha <- stats::setNames(rep(alpha, length(lev)), lev)
    if (isTRUE(use_proj_style_override) && isTRUE(is_proj_helper_group) && ("Reference" %in% names(grp_alpha))) {
      grp_alpha["Reference"] <- min(1, alpha + 0.2)
    }
  } else {
    grp_cols <- character()
    grp_shapes <- numeric()
    grp_sizes <- numeric()
    grp_strokes <- numeric()
    grp_alpha <- numeric()
    use_proj_style_override <- FALSE
  }
  .pca_plot_max_abs <- function(x, y) {
    v <- suppressWarnings(as.numeric(c(x, y)))
    v <- v[is.finite(v) & !is.na(v)]
    if (!length(v)) return(NA_real_)
    max(abs(v), na.rm = TRUE)
  }
  .pca_plot_limit_with_sd <- function(x, y, sds = numeric()) {
    xv <- suppressWarnings(as.numeric(x))
    yv <- suppressWarnings(as.numeric(y))
    ok <- is.finite(xv) & !is.na(xv) & is.finite(yv) & !is.na(yv)
    if (!any(ok)) return(NA_real_)
    mx <- mean(xv[ok], na.rm = TRUE)
    my <- mean(yv[ok], na.rm = TRUE)
    sx <- stats::sd(xv[ok], na.rm = TRUE)
    sy <- stats::sd(yv[ok], na.rm = TRUE)
    if (!is.finite(sx) || is.na(sx)) sx <- 0
    if (!is.finite(sy) || is.na(sy)) sy <- 0
    base_lim <- .pca_plot_max_abs(xv[ok], yv[ok])
    if (!length(sds)) return(base_lim)
    sds <- suppressWarnings(as.numeric(sds))
    sds <- sds[is.finite(sds) & !is.na(sds) & sds > 0]
    if (!length(sds)) return(base_lim)
    lim_x_sd <- max(abs(c(mx - sds * sx, mx + sds * sx)), na.rm = TRUE)
    lim_y_sd <- max(abs(c(my - sds * sy, my + sds * sy)), na.rm = TRUE)
    max(base_lim, lim_x_sd, lim_y_sd, na.rm = TRUE)
  }
  global_lim <- NA_real_
  if (identical(equal.lim, "global")) {
    lim_pool <- numeric()
    for (pr in pair_list) {
      xnm <- paste0("PC", as_int(pr[1]))
      ynm <- paste0("PC", as_int(pr[2]))
      if (!(xnm %in% names(dt)) || !(ynm %in% names(dt))) next
      l0 <- .pca_plot_limit_with_sd(dt[[xnm]], dt[[ynm]], sds = line_sd)
      if (is.finite(l0) && !is.na(l0) && l0 > 0) lim_pool <- c(lim_pool, l0)
    }
    lim_pool <- lim_pool[is.finite(lim_pool) & !is.na(lim_pool) & lim_pool > 0]
    if (length(lim_pool)) global_lim <- max(lim_pool, na.rm = TRUE)
  }
  for (pr in pair_list) {
    xk <- as_int(pr[1]); yk <- as_int(pr[2])
    xnm <- paste0("PC", xk); ynm <- paste0("PC", yk)
    sub <- dt[is.finite(get(xnm)) & !is.na(get(xnm)) & is.finite(get(ynm)) & !is.na(get(ynm))]
    if (!nrow(sub)) next
    if (has_group && ("Group" %in% names(sub))) {
      sub[, .draw_unknown := tolower(trimws(as.character(Group))) == "unknown"]
      sub[, .draw_target := trimws(as.character(Group)) == "Target"]
      data.table::setorderv(sub, c(".draw_unknown", ".draw_target"), c(-1L, 1L), na.last = TRUE)
      sub[, c(".draw_unknown", ".draw_target") := NULL]
    }

    mx <- mean(sub[[xnm]], na.rm = TRUE); sx <- stats::sd(sub[[xnm]], na.rm = TRUE)
    my <- mean(sub[[ynm]], na.rm = TRUE); sy <- stats::sd(sub[[ynm]], na.rm = TRUE)
    if (!is.finite(sx) || is.na(sx)) sx <- 0
    if (!is.finite(sy) || is.na(sy)) sy <- 0

    xlab <- .pca_plot_axis_label(xk, eval_named)
    ylab <- .pca_plot_axis_label(yk, eval_named)
    p <- ggplot2::ggplot(sub, ggplot2::aes(x = .data[[xnm]], y = .data[[ynm]]))
    if (has_group) {
      special_ids <- character()
      if (isTRUE(refqc.helper)) {
        special_ids <- unique(c(refqc_ambiguous_ids, refqc_drop_ids))
      }
      sub_base <- if (length(special_ids)) sub[!(IID_USE %in% special_ids)] else sub
      if (isTRUE(refqc.helper)) {
        base_shape <- suppressWarnings(as.numeric(point.shape))[1]
        if (!is.finite(base_shape) || is.na(base_shape)) base_shape <- 16
        ds <- if (length(special_ids)) sub[IID_USE %in% special_ids] else sub[0]
        da <- if (nrow(ds)) ds[IID_USE %in% refqc_ambiguous_ids] else ds
        dd <- if (nrow(ds)) ds[IID_USE %in% refqc_drop_ids] else ds
        base_size <- point.size
        mark_size <- max(3, point.size * 1.45)
        lg_pt_size <- max(3, as.numeric(point.size) * 2)
        mark_stroke <- max(1.4, point.stroke)
        amb_stroke <- max(0.7, point.stroke * 0.75)
        grp_cols_base <- grp_cols
        grp_cols_mark <- grp_cols_overlay
        for (g in lev) {
          g_chr <- as.character(g)
          if (is.na(g_chr) || !nzchar(g_chr)) next
          col_base <- unname(grp_cols_base[g_chr])
          if (!is.character(col_base) || !length(col_base) || is.na(col_base) || !nzchar(col_base)) col_base <- "grey70"
          d0 <- sub_base[as.character(Group) == g_chr]
          if (nrow(d0)) {
            p <- p + ggplot2::geom_point(
              data = d0,
              mapping = ggplot2::aes(x = .data[[xnm]], y = .data[[ynm]]),
              inherit.aes = FALSE,
              colour = col_base,
              shape = base_shape,
              size = base_size,
              stroke = point.stroke,
              alpha = alpha,
              show.legend = FALSE
            )
          }
        }

        for (g in lev) {
          g_chr <- as.character(g)
          if (is.na(g_chr) || !nzchar(g_chr)) next
          col_amb <- unname(grp_cols_mark[g_chr])
          if (!is.character(col_amb) || !length(col_amb) || is.na(col_amb) || !nzchar(col_amb)) col_amb <- "grey20"
          d1 <- da[as.character(Group) == g_chr]
          if (nrow(d1)) {
            p <- p + ggplot2::geom_point(
              data = d1,
              mapping = ggplot2::aes(x = .data[[xnm]], y = .data[[ynm]]),
              inherit.aes = FALSE,
              fill = col_amb,
              colour = "grey40",
              shape = 24,
              size = mark_size,
              stroke = amb_stroke,
              alpha = 1,
              show.legend = FALSE
            )
          }
        }

        for (g in lev) {
          g_chr <- as.character(g)
          if (is.na(g_chr) || !nzchar(g_chr)) next
          col_drop <- unname(grp_cols_mark[g_chr])
          if (!is.character(col_drop) || !length(col_drop) || is.na(col_drop) || !nzchar(col_drop)) col_drop <- "grey20"
          d2 <- dd[as.character(Group) == g_chr]
          if (nrow(d2)) {
            p <- p + ggplot2::geom_point(
              data = d2,
              mapping = ggplot2::aes(x = .data[[xnm]], y = .data[[ynm]]),
              inherit.aes = FALSE,
              colour = col_drop,
              shape = 4,
              size = mark_size,
              stroke = mark_stroke,
              alpha = 1,
              show.legend = FALSE
            )
          }
        }

        d_group_legend <- data.table::data.table(
          .x = rep(sub[[xnm]][1], length(lev)),
          .y = rep(sub[[ynm]][1], length(lev)),
          .legend_group = factor(lev, levels = lev)
        )
        p <- p +
          ggplot2::geom_point(
            data = d_group_legend,
            mapping = ggplot2::aes(x = .x, y = .y, colour = .legend_group),
            inherit.aes = FALSE,
            shape = 16,
            size = lg_pt_size,
            alpha = 0,
            show.legend = TRUE
          ) +
          ggplot2::scale_color_manual(
            values = stats::setNames(as.character(grp_cols_mark[lev]), lev),
            breaks = lev,
            drop = FALSE,
            name = group_legend_title,
            guide = ggplot2::guide_legend(
              order = 1,
              byrow = TRUE,
              override.aes = list(
                shape = 16,
                alpha = 1,
                size = lg_pt_size,
                colour = as.character(grp_cols_mark[lev])
              )
            )
          )

        if (nrow(da) || nrow(dd)) {
          status_levels <- c("Ambiguous", refqc_drop_label)
          status_keep <- c(nrow(da) > 0L, nrow(dd) > 0L)
          d_status_legend <- data.table::data.table(
            .x = rep(sub[[xnm]][1], length(status_levels)),
            .y = rep(sub[[ynm]][1], length(status_levels)),
            .qc_status = factor(status_levels, levels = status_levels)
          )
          d_status_legend <- d_status_legend[status_keep]
          status_breaks <- as.character(d_status_legend$.qc_status)
          shape_vals <- ifelse(status_breaks == "Ambiguous", 24, 4)
          fill_vals <- ifelse(status_breaks == "Ambiguous", "white", NA_character_)
          stroke_vals <- ifelse(status_breaks == "Ambiguous", amb_stroke, mark_stroke)
          shape_map <- c(Ambiguous = 24)
          shape_map[refqc_drop_label] <- 4
          p <- p +
            ggplot2::geom_point(
              data = d_status_legend,
              mapping = ggplot2::aes(x = .x, y = .y, shape = .qc_status),
              inherit.aes = FALSE,
              colour = "grey20",
              fill = "grey20",
              size = mark_size,
              stroke = mark_stroke,
              alpha = 0,
              na.rm = TRUE,
              show.legend = TRUE
            ) +
            ggplot2::scale_shape_manual(
              values = shape_map,
              breaks = status_breaks,
              name = NULL,
              drop = FALSE,
              guide = ggplot2::guide_legend(
                order = 2,
                override.aes = list(
                  shape = shape_vals,
                  colour = rep("grey20", length(status_breaks)),
                  fill = fill_vals,
                  alpha = rep(1, length(status_breaks)),
                  size = rep(lg_pt_size, length(status_breaks)),
                  stroke = stroke_vals
                )
              )
            )
        }
      } else if (isTRUE(use_proj_style_override)) {
        sub[, .proj_ref := .pca_plot_as_ref_flag(reference) %in% TRUE]
        sub[, .proj_class := data.table::fifelse(.proj_ref %in% TRUE, proj_ref_label, proj_tar_label)]
        sub[, .proj_class := factor(.proj_class, levels = c(proj_ref_label, proj_tar_label))]
        sub[, .proj_fill := NA_character_]
        sub[.proj_ref %in% FALSE, .proj_fill := "white"]
        sub[, .proj_stroke := unname(grp_strokes[as.character(Group)])]
        sub[, .proj_col_ref := unname(.pca_plot_lighten_hex(grp_cols_overlay, factor = 0.6)[as.character(Group)])]
        sub[is.na(.proj_col_ref) | !nzchar(.proj_col_ref), .proj_col_ref := "grey70"]
        d_ref <- sub[.proj_ref %in% TRUE]
        d_tar <- sub[.proj_ref %in% FALSE]
        d_grp_legend <- data.table::data.table(
          .x = rep(sub[[xnm]][1], length(lev)),
          .y = rep(sub[[ynm]][1], length(lev)),
          Group = factor(lev, levels = lev)
        )
        d_proj_legend <- data.table::data.table(
          .x = rep(sub[[xnm]][1], 2L),
          .y = rep(sub[[ynm]][1], 2L),
          .proj_class = factor(c(proj_ref_label, proj_tar_label), levels = c(proj_ref_label, proj_tar_label))
        )
        if (nrow(d_ref)) {
          p <- p + ggplot2::geom_point(
            data = d_ref,
            mapping = ggplot2::aes(shape = .proj_class, size = Group, stroke = .proj_stroke),
            inherit.aes = TRUE,
            colour = d_ref$.proj_col_ref,
            fill = NA,
            alpha = alpha,
            show.legend = FALSE
          )
        }
        if (nrow(d_tar)) {
          p <- p + ggplot2::geom_point(
            data = d_tar,
            mapping = ggplot2::aes(color = Group, shape = .proj_class, size = Group, fill = .proj_fill, stroke = .proj_stroke),
            inherit.aes = TRUE,
            alpha = alpha,
            show.legend = FALSE
          )
        }
        p <- p +
          ggplot2::geom_point(
            data = d_grp_legend,
            mapping = ggplot2::aes(x = .x, y = .y, color = Group, size = Group),
            inherit.aes = FALSE,
            alpha = 0,
            shape = 16,
            show.legend = TRUE
          ) +
          ggplot2::geom_point(
            data = d_proj_legend,
            mapping = ggplot2::aes(x = .x, y = .y, shape = .proj_class),
            inherit.aes = FALSE,
            colour = "grey20",
            fill = "white",
            alpha = 0,
            size = point.size,
            stroke = point.stroke,
            show.legend = TRUE
          ) +
          ggplot2::scale_shape_manual(values = stats::setNames(c(16, 22), c(proj_ref_label, proj_tar_label)), name = NULL, drop = FALSE) +
          ggplot2::scale_fill_identity(guide = "none") +
          ggplot2::scale_continuous_identity(aesthetics = "stroke", guide = "none")
      } else {
        p <- p + ggplot2::geom_point(ggplot2::aes(color = Group, shape = Group, size = Group, stroke = Group), alpha = alpha)
      }
      if (!isTRUE(refqc.helper)) {
        p <- p + ggplot2::scale_color_manual(values = grp_cols, drop = FALSE, name = group_legend_title)
      }
      if (isTRUE(use_proj_style_override)) {
        p <- p + ggplot2::scale_size_manual(values = grp_sizes, drop = FALSE, name = group_legend_title)
      } else if (!isTRUE(refqc.helper)) {
        p <- p +
          ggplot2::scale_size_manual(values = grp_sizes, drop = FALSE, name = group_legend_title) +
          ggplot2::scale_discrete_manual(aesthetics = "stroke", values = grp_strokes, drop = FALSE, name = group_legend_title)
      }
      if (!isTRUE(use_proj_style_override) && !isTRUE(refqc.helper)) {
        p <- p + ggplot2::scale_shape_manual(values = grp_shapes, drop = FALSE, name = group_legend_title)
      }
    } else {
      # No group: allow direct per-point or scalar control by point.color/point.shape/group.size.
      n_sub <- nrow(sub)
      col_vec <- point.color
      if (is.list(col_vec) && !is.data.frame(col_vec) && !data.table::is.data.table(col_vec)) col_vec <- unlist(col_vec, use.names = TRUE)
      if (is.null(col_vec) || length(col_vec) == 0L) {
        sub[, .pt_col := "grey20"]
      } else {
        cv <- as.character(col_vec)
        cv <- cv[!is.na(cv) & nzchar(cv)]
        if (!length(cv)) {
          sub[, .pt_col := "grey20"]
        } else if (length(cv) == 1L && !.pca_plot_is_valid_color(cv[1])) {
          pal <- tryCatch(get.colors(cv[1], n_sub, discrete = FALSE, far = TRUE, random = FALSE, plot = FALSE, silent = TRUE), error = function(e) NULL)
          if (is.null(pal) || !length(pal)) pal <- "grey20"
          sub[, .pt_col := rep_len(as.character(pal), n_sub)]
        } else {
          sub[, .pt_col := rep_len(cv, n_sub)]
        }
      }
      sub[, .pt_col_orig := .pt_col]
      if (isTRUE(fade_on)) {
        sub[, .pt_col := .pca_plot_lighten_hex(.pt_col_orig)]
      }
      shp_in <- point.shape
      if (is.list(shp_in) && !is.data.frame(shp_in) && !data.table::is.data.table(shp_in)) {
        shp_in <- unlist(shp_in, use.names = TRUE)
      }
      shp_vec <- suppressWarnings(as.numeric(shp_in))
      shp_vec <- shp_vec[is.finite(shp_vec) & !is.na(shp_vec)]
      if (!length(shp_vec)) shp_vec <- 16
      sub[, .pt_shape := rep_len(shp_vec, n_sub)]
      sz_in <- if (!is.null(group.size) && length(group.size)) group.size else point.size
      if (is.list(sz_in) && !is.data.frame(sz_in) && !data.table::is.data.table(sz_in)) {
        sz_in <- unlist(sz_in, use.names = TRUE)
      }
      sz_base <- suppressWarnings(as.numeric(sz_in))
      sz_base <- sz_base[is.finite(sz_base) & !is.na(sz_base) & sz_base > 0]
      if (!length(sz_base)) sz_base <- point.size
      sub[, .pt_size := rep_len(sz_base, n_sub)]
      if (isTRUE(refqc.helper)) {
        special_ids <- unique(c(refqc_ambiguous_ids, refqc_drop_ids))
        sub_base <- if (length(special_ids)) sub[!(IID_USE %in% special_ids)] else sub
        da <- if (length(refqc_ambiguous_ids)) sub[IID_USE %in% refqc_ambiguous_ids] else sub[0]
        dd <- if (length(refqc_drop_ids)) sub[IID_USE %in% refqc_drop_ids] else sub[0]
        mark_size <- max(3, point.size * 1.45)
        mark_stroke <- max(1.4, point.stroke)
        amb_stroke <- max(0.7, point.stroke * 0.75)
        p <- p +
          ggplot2::geom_point(
            data = sub_base,
            mapping = ggplot2::aes(x = .data[[xnm]], y = .data[[ynm]], color = .pt_col, shape = .pt_shape, size = .pt_size),
            inherit.aes = FALSE,
            alpha = alpha,
            stroke = point.stroke,
            show.legend = FALSE
          )
        if (nrow(da)) {
          p <- p + ggplot2::geom_point(
            data = da,
            mapping = ggplot2::aes(x = .data[[xnm]], y = .data[[ynm]]),
            inherit.aes = FALSE,
            fill = da$.pt_col_orig,
            colour = "grey40",
            shape = 24,
            size = mark_size,
            stroke = amb_stroke,
            alpha = 1,
            show.legend = FALSE
          )
        }
        if (nrow(dd)) {
          p <- p + ggplot2::geom_point(
            data = dd,
            mapping = ggplot2::aes(x = .data[[xnm]], y = .data[[ynm]]),
            inherit.aes = FALSE,
            colour = dd$.pt_col_orig,
            shape = 4,
            size = mark_size,
            stroke = mark_stroke,
            alpha = 1,
            show.legend = FALSE
          )
        }
        p <- p +
          ggplot2::scale_color_identity() +
          ggplot2::scale_shape_identity() +
          ggplot2::scale_size_identity()
      } else {
        p <- p +
          ggplot2::geom_point(ggplot2::aes(color = .pt_col, shape = .pt_shape, size = .pt_size), alpha = alpha, stroke = point.stroke) +
          ggplot2::scale_color_identity() +
          ggplot2::scale_shape_identity() +
          ggplot2::scale_size_identity()
      }
    }
    if (n_line > 0L) {
      for (i in seq_len(n_line)) {
        s <- line_sd[i]
        p <- p +
          ggplot2::geom_vline(xintercept = c(mx - s * sx, mx + s * sx), color = line_col[i], linetype = line_typ[i], linewidth = line_lwd[i]) +
          ggplot2::geom_hline(yintercept = c(my - s * sy, my + s * sy), color = line_col[i], linetype = line_typ[i], linewidth = line_lwd[i])
      }
    }
    p <- p +
      ggplot2::xlab(xlab) +
      ggplot2::ylab(ylab)
    lim_use <- NA_real_
    if (identical(equal.lim, "global")) {
      lim_use <- global_lim
    } else if (identical(equal.lim, "panel")) {
      lim_use <- .pca_plot_limit_with_sd(sub[[xnm]], sub[[ynm]], sds = line_sd)
    }
    if (is.finite(lim_use) && !is.na(lim_use) && lim_use > 0) {
      lim_use <- lim_use * 1.02
      p <- p + ggplot2::coord_fixed(ratio = 1, xlim = c(-lim_use, lim_use), ylim = c(-lim_use, lim_use))
    } else {
      p <- p + ggplot2::coord_fixed(ratio = 1)
    }
    p <- p +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.title = ggplot2::element_text(size = 14),
        axis.text = ggplot2::element_text(size = 12),
        legend.title = ggplot2::element_text(size = 12),
        legend.text = ggplot2::element_text(size = 12),
        panel.border = ggplot2::element_rect(fill = NA, color = "grey20", linewidth = 0.35),
        axis.line = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_line(color = "grey20", linewidth = 0.35)
      )
    if (has_group) {
      lg_pt_size <- max(3, as.numeric(point.size) * 2)
      if (isTRUE(refqc.helper)) {
        p <- p +
          ggplot2::theme(
            legend.position = legend.position,
            legend.background = ggplot2::element_rect(fill = scales::alpha("white", 0.4), color = "grey20"),
            legend.box.background = ggplot2::element_blank(),
            legend.key = ggplot2::element_rect(fill = "white", color = NA)
          )
      } else if (isTRUE(use_proj_style_override)) {
        leg_cols <- as.character(unname(grp_cols[lev]))
        leg_cols[is.na(leg_cols) | !nzchar(leg_cols)] <- "grey20"
        p <- p +
          ggplot2::guides(
            color = ggplot2::guide_legend(
              byrow = TRUE,
              override.aes = list(
                shape = rep(16, length(lev)),
                size = rep(lg_pt_size, length(lev)),
                alpha = rep(1, length(lev)),
                colour = leg_cols
              )
            ),
            shape = ggplot2::guide_legend(
              title = NULL,
              override.aes = list(color = "grey20", fill = c(NA, "white"), size = lg_pt_size, alpha = 1, stroke = point.stroke)
            ),
            size = "none",
            stroke = "none"
          ) +
          ggplot2::theme(
            legend.position = legend.position,
            legend.background = ggplot2::element_rect(fill = scales::alpha("white", 0.4), color = "grey20"),
            legend.box.background = ggplot2::element_blank(),
            legend.key = ggplot2::element_rect(fill = "white", color = NA)
          )
      } else {
        p <- p +
          ggplot2::guides(
            color = ggplot2::guide_legend(
              byrow = TRUE,
              override.aes = list(
                shape = unname(grp_shapes),
                size = rep(lg_pt_size, length(grp_shapes)),
                stroke = unname(grp_strokes),
                alpha = 1
              )
            ),
            shape = "none",
            size = "none",
            stroke = "none"
          ) +
          ggplot2::theme(
            legend.position = legend.position,
            legend.background = ggplot2::element_rect(fill = scales::alpha("white", 0.4), color = "grey20"),
            legend.box.background = ggplot2::element_blank(),
            legend.key = ggplot2::element_rect(fill = "white", color = NA)
          )
      }
      p <- p +
        ggplot2::theme(
          legend.title = ggplot2::element_text(size = 12),
          legend.text = ggplot2::element_text(size = 12)
        )
    } else {
      p <- p + ggplot2::theme(legend.position = "none")
    }
    plot_list[[paste0("PC", xk, "_PC", yk)]] <- p
  }
  if (!length(plot_list)) stop("No valid PCA plot panels generated.", call. = FALSE)

  list(
    plot = plot_list,
    outlier = list(
      `3sd` = as.character(unique(out3)),
      `5sd` = as.character(unique(out5))
    )
  )
}

