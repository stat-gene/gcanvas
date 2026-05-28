# Sample-level heterozygosity vs. missingness QC scatter plot.

#' Heterozygosity vs. missingness plot
#'
#' Computes (or accepts precomputed) per-sample heterozygosity rates and
#' missingness rates from PLINK genotype files and plots them with QC
#' thresholds drawn on top.
#'
#' @param het Either a path to a PLINK heterozygosity file or a data.frame.
#' @param miss Either a path to a PLINK missingness file or a data.frame.
#' @param bfile Optional PLINK1 bfile prefix (computes `het`/`miss` if both
#'   `het` and `miss` are `NULL`).
#' @param pfile Optional PLINK2 pfile prefix (alternative to `bfile`).
#' @param plink PLINK2 binary name / path.
#' @param plink.version One of `"auto"`, `"plink2"`, `"plink"`.
#' @param plink.out Logical or path. Persist PLINK outputs.
#' @param show.plinklog Logical. Stream PLINK's log to the R console.
#' @param line.sd,line.sd.color,line.sd.type,line.sd.linewidth Threshold
#'   line style controls.
#' @param group Optional grouping column / vector for coloring.
#' @param group.color Optional color override per group.
#' @param show.legend Logical. Draw the legend.
#' @param legend.position Where the legend is placed.
#' @param silent Logical. Suppress progress notes.
#'
#' @return A `ggplot` object.
#' @export
hetmiss <- function(het = NULL,
                    miss = NULL,
                    bfile = NULL,
                    pfile = NULL,
                    plink = "plink2",
                    plink.version = "auto",
                    plink.out = FALSE,
                    show.plinklog = TRUE,
                    line.sd = c(5, 3),
                    line.sd.color = c("#EF476F", "#258AB2"),
                    line.sd.type = "dashed",
                    line.sd.linewidth = 0.2,
                    group = NULL,
                    group.color = NULL,
                    show.legend = TRUE,
                    legend.position = c("bottom", "top", "left", "right", "left.outside", "right.outside"),
                    silent = FALSE) {
  require_pkg(c("data.table", "ggplot2"))
  silent <- isTRUE(silent)
  show.plinklog <- isTRUE(show.plinklog)
  if (isTRUE(silent)) show.plinklog <- FALSE
  show.legend <- isTRUE(show.legend)
  legend.position <- match.arg(tolower(as.character(legend.position)[1]), choices = c("bottom", "top", "left", "right", "left.outside", "right.outside"))

  .hetmiss_canon <- function(nm) toupper(gsub("[^A-Z0-9#]+", "", as.character(nm)))
  .hetmiss_pick_col <- function(nm, keys, default = NA_character_) {
    cn <- .hetmiss_canon(nm)
    for (k in keys) {
      hit <- which(cn == toupper(k))
      if (length(hit)) return(nm[hit[1]])
    }
    default
  }
  .hetmiss_read_any_dt <- function(x) {
    if (is.null(x)) return(NULL)
    if (data.table::is.data.table(x)) return(data.table::copy(x))
    if (is.data.frame(x)) return(data.table::as.data.table(x))
    if (is.character(x) && length(x) == 1L && file.exists(x)) {
      return(data.table::fread(x, data.table = TRUE, check.names = FALSE, showProgress = FALSE))
    }
    stop("Input must be data.table/data.frame or existing file path.", call. = FALSE)
  }
  .hetmiss_read_group_dt <- function(g) {
    if (is.null(g)) return(NULL)
    if (data.table::is.data.table(g)) dtg <- data.table::copy(g)
    else if (is.data.frame(g)) dtg <- data.table::as.data.table(g)
    else stop("group must be data.table/data.frame.", call. = FALSE)
    if (!nrow(dtg)) return(NULL)
    nms <- names(dtg)
    i_iid <- .hetmiss_pick_col(nms, c("IID", "#IID", "ID"))
    if (is.na(i_iid)) {
      if (ncol(dtg) >= 2L) {
        i_iid <- nms[1]
      } else {
        stop("group must include IID column (or at least 2 columns: IID + group).", call. = FALSE)
      }
    }
    i_grp <- setdiff(nms, c(i_iid, .hetmiss_pick_col(nms, c("FID", "#FID"))))
    i_grp <- if (length(i_grp)) i_grp[1] else NA_character_
    if (is.na(i_grp)) stop("group must include a group column.", call. = FALSE)
    out <- data.table::data.table(
      IID = as.character(dtg[[i_iid]]),
      Group = as.character(dtg[[i_grp]])
    )
    out <- out[!is.na(IID) & nzchar(IID)]
    out[]
  }
  .hetmiss_normalize_version <- function(v) {
    if (is.null(v) || length(v) == 0L) return("auto")
    vv <- tolower(trimws(as.character(v)[1]))
    if (vv %in% c("auto", "1", "plink1", "1.9", "plink2", "2", "2.0")) {
      if (vv %in% c("1", "plink1", "1.9")) return("plink1")
      if (vv %in% c("2", "plink2", "2.0")) return("plink2")
      return("auto")
    }
    "auto"
  }
  .hetmiss_detect_plink_version <- function(dt_het, dt_miss) {
    nmh <- names(dt_het %||% data.table::data.table())
    nmm <- names(dt_miss %||% data.table::data.table())
    cnh <- .hetmiss_canon(nmh)
    cnm <- .hetmiss_canon(nmm)
    if (any(cnh %in% c("NNM")) || any(cnm %in% c("NMISS", "NGENO"))) return("plink1")
    if (any(cnh %in% c("OBSCT")) || any(cnm %in% c("MISSINGCT", "OBSCT"))) return("plink2")
    "auto"
  }
  .hetmiss_normalize_het <- function(dt) {
    if (is.null(dt) || !nrow(dt)) stop("het is empty.", call. = FALSE)
    nm <- names(dt)
    i_fid <- .hetmiss_pick_col(nm, c("FID", "#FID"))
    i_iid <- .hetmiss_pick_col(nm, c("IID", "#IID", "ID"))
    i_ohom <- .hetmiss_pick_col(nm, c("OHOM"))
    i_nnm <- .hetmiss_pick_col(nm, c("NNM", "OBSCT"))
    if (is.na(i_iid) || is.na(i_ohom) || is.na(i_nnm)) {
      stop("het columns are not recognized. Need IID/#IID(ID), O(HOM), and N(NM)/OBS_CT.", call. = FALSE)
    }
    out <- data.table::data.table(
      FID = if (!is.na(i_fid)) as.character(dt[[i_fid]]) else NA_character_,
      IID = as.character(dt[[i_iid]]),
      O_HOM = suppressWarnings(as.numeric(dt[[i_ohom]])),
      OBS_CT = suppressWarnings(as.numeric(dt[[i_nnm]]))
    )
    out <- out[!is.na(IID) & nzchar(IID)]
    out[, het_rate := (OBS_CT - O_HOM) / OBS_CT]
    out[!is.finite(het_rate) | is.na(het_rate), het_rate := NA_real_]
    out[]
  }
  .hetmiss_normalize_miss <- function(dt) {
    if (is.null(dt) || !nrow(dt)) stop("miss is empty.", call. = FALSE)
    nm <- names(dt)
    i_fid <- .hetmiss_pick_col(nm, c("FID", "#FID"))
    i_iid <- .hetmiss_pick_col(nm, c("IID", "#IID", "ID"))
    i_fmiss <- .hetmiss_pick_col(nm, c("FMISS"))
    i_nmiss <- .hetmiss_pick_col(nm, c("NMISS", "MISSINGCT"))
    i_ngeno <- .hetmiss_pick_col(nm, c("NGENO", "OBSCT"))
    if (is.na(i_iid) || is.na(i_fmiss)) {
      stop("miss columns are not recognized. Need IID/#IID(ID) and F_MISS.", call. = FALSE)
    }
    out <- data.table::data.table(
      FID = if (!is.na(i_fid)) as.character(dt[[i_fid]]) else NA_character_,
      IID = as.character(dt[[i_iid]]),
      F_MISS = suppressWarnings(as.numeric(dt[[i_fmiss]])),
      N_MISS = if (!is.na(i_nmiss)) suppressWarnings(as.numeric(dt[[i_nmiss]])) else NA_real_,
      N_GENO = if (!is.na(i_ngeno)) suppressWarnings(as.numeric(dt[[i_ngeno]])) else NA_real_
    )
    out <- out[!is.na(IID) & nzchar(IID)]
    out[!is.finite(F_MISS) | is.na(F_MISS), F_MISS := NA_real_]
    out[]
  }
  .hetmiss_run_plink <- function(bfile0, pfile0, plink0, ver0, plink_out0 = FALSE) {
    .hetmiss_strip_quotes <- function(s) {
      s <- trimws(as.character(s))
      if (!nzchar(s)) return(s)
      s <- gsub("^[\"']|[\"']$", "", s)
      s
    }
    .hetmiss_alias_target_from_rc <- function(cmd, rc_files = c("~/.zshrc", "~/.zprofile", "~/.bashrc", "~/.bash_profile")) {
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
        rhs <- .hetmiss_strip_quotes(rhs)
        if (!nzchar(rhs)) next
        tok <- strsplit(rhs, "\\s+")[[1]]
        if (!length(tok)) next
        bin <- .hetmiss_strip_quotes(tok[1])
        if (!nzchar(bin)) next
        w <- Sys.which(bin)
        if (nzchar(w)) return(w)
        if (file.exists(path.expand(bin))) return(abs_path(path.expand(bin)))
      }
      NA_character_
    }
    .hetmiss_resolve_exec <- function(cmd) {
      cmd <- as.character(cmd)[1]
      if (is.na(cmd) || !nzchar(cmd)) return(NA_character_)
      if (file.exists(cmd)) return(abs_path(cmd))

      w <- Sys.which(cmd)
      if (nzchar(w)) return(w)

      # Fallback search for common binary locations when R PATH is minimal.
      cand_dir <- c(
        path.expand("~/bin"),
        "/opt/homebrew/bin",
        "/usr/local/bin",
        "/usr/bin",
        "/bin"
      )
      cand <- file.path(cand_dir, cmd)
      hit <- cand[file.exists(cand)]
      if (length(hit)) return(abs_path(hit[1]))

      # Try simple alias definitions in shell rc files.
      a <- .hetmiss_alias_target_from_rc(cmd)
      if (!is.na(a) && nzchar(a)) return(a)

      # Final check: ask login shell where the command is.
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
        # command -v may print alias text; extract command token if possible.
        if (grepl("aliased to", last, fixed = TRUE)) {
          rhs <- sub(".*aliased to\\s*", "", last)
          rhs <- gsub("^['`\"]|['`\"]$", "", rhs)
          tok <- strsplit(rhs, "\\s+")[[1]]
          cand2 <- tok[1]
          if (!is.na(cand2) && nzchar(cand2)) {
            p2 <- Sys.which(cand2)
            if (nzchar(p2)) return(p2)
            if (file.exists(path.expand(cand2))) return(abs_path(path.expand(cand2)))
          }
        } else {
          if (file.exists(last)) return(abs_path(last))
          p3 <- Sys.which(last)
          if (nzchar(p3)) return(p3)
        }
      }
      NA_character_
    }

    if (!is.null(bfile0) && !is.null(pfile0)) {
      stop("Use either bfile or pfile, not both.", call. = FALSE)
    }
    if (is.null(bfile0) && is.null(pfile0)) {
      stop("Provide het/miss or bfile/pfile.", call. = FALSE)
    }
    keep_plink <- FALSE
    cache_dir0 <- .gcanvas_default_cache_dir(scope = "hetmiss", anchor = if (!is.null(bfile0)) bfile0 else pfile0)
    dir.create(cache_dir0, recursive = TRUE, showWarnings = FALSE)
    .gcanvas_register_cache_dir(cache_dir0)
    out_pref <- file.path(cache_dir0, sprintf("hetmiss_%s_%s", format(Sys.time(), "%Y%m%d_%H%M%S"), as.integer(stats::runif(1, 1, 1e9))))
    if (!is.null(plink_out0) && length(plink_out0) > 0L) {
      if (is.logical(plink_out0) && length(plink_out0) == 1L) {
        if (isTRUE(plink_out0)) {
          cache_dir <- .gcanvas_default_cache_dir(scope = "hetmiss", anchor = if (!is.null(bfile0)) bfile0 else pfile0)
          dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
          .gcanvas_register_cache_dir(cache_dir)
          out_pref <- file.path(cache_dir, sprintf("hetmiss_%s_%s", format(Sys.time(), "%Y%m%d_%H%M%S"), as.integer(stats::runif(1, 1, 1e9))))
          keep_plink <- TRUE
        }
      } else {
        out_chr <- as.character(plink_out0)[1]
        if (!is.na(out_chr) && nzchar(out_chr) && !(tolower(trimws(out_chr)) %in% c("false", "null", "na"))) {
          out_chr <- abs_path(path.expand(out_chr))
          if (dir.exists(out_chr) || grepl("[/\\\\]$", as.character(plink_out0)[1])) {
            dir.create(out_chr, recursive = TRUE, showWarnings = FALSE)
            .gcanvas_register_cache_dir(out_chr)
            out_pref <- file.path(out_chr, sprintf("hetmiss_%s", format(Sys.time(), "%Y%m%d_%H%M%S")))
          } else {
            dir.create(dirname(out_chr), recursive = TRUE, showWarnings = FALSE)
            .gcanvas_register_cache_dir(dirname(out_chr))
            out_pref <- out_chr
          }
          keep_plink <- TRUE
        }
      }
    }
    args <- character()
    if (!is.null(bfile0)) args <- c(args, "--bfile", as.character(bfile0)[1])
    if (!is.null(pfile0)) {
      pfx <- as.character(pfile0)[1]
      use_vzs <- file.exists(paste0(pfx, ".pvar.zst")) && !file.exists(paste0(pfx, ".pvar"))
      args <- c(args, "--pfile", pfx, if (isTRUE(use_vzs)) "vzs")
    }
    args <- c(args, "--het", "--missing", "--allow-no-sex", "--out", out_pref)
    plink_cmd <- as.character(plink0)[1]
    if (is.na(plink_cmd) || !nzchar(plink_cmd)) {
      stop("PLINK command is empty. Set `plink` to a valid executable path or name.", call. = FALSE)
    }
    plink_exec <- .hetmiss_resolve_exec(plink_cmd)
    if (is.na(plink_exec) || !nzchar(plink_exec)) {
      stop(
        sprintf(
          "PLINK executable not found for '%s'. Alias-only setup is not available to R here; set `plink` to an executable path.",
          plink_cmd
        ),
        call. = FALSE
      )
    }
    plink_run <- plink_exec
    .gcanvas_note(
      "gcanvas::hetmiss",
      sprintf("Start: mode=%s | plink=%s | plink.version=%s",
              if (!is.null(bfile0)) "bfile" else "pfile",
              as.character(plink0)[1],
              as.character(ver0)[1]),
      silent = silent
    )
    if (isTRUE(keep_plink)) {
      .gcanvas_note("gcanvas::hetmiss", sprintf("PLINK outputs kept: %s", out_pref), silent = silent)
    }
    .exec_plink <- function(args_one) {
      cmd_line <- paste(c(shQuote(plink_run), vapply(args_one, shQuote, character(1))), collapse = " ")
      args_run <- gsub("\\$", "\\\\\\$", as.character(args_one))
      if (isTRUE(show.plinklog)) {
        st <- tryCatch(
          suppressWarnings(system2(command = plink_run, args = args_run, stdout = "", stderr = "")),
          error = function(e) {
            stop(
              sprintf("Failed to run PLINK command: %s\nCommand: %s", conditionMessage(e), cmd_line),
              call. = FALSE
            )
          }
        )
        status <- as_int(st)
        if (!identical(as_int(status), 0L)) {
          stop("PLINK command failed.\nCommand: ", cmd_line, call. = FALSE)
        }
      } else {
        rc <- tryCatch(
          suppressWarnings(system2(command = plink_run, args = args_run, stdout = TRUE, stderr = TRUE)),
          error = function(e) {
            stop(
              sprintf("Failed to run PLINK command: %s\nCommand: %s", conditionMessage(e), cmd_line),
              call. = FALSE
            )
          }
        )
        status <- attr(rc, "status") %||% 0L
        if (!identical(as_int(status), 0L)) {
          msg <- if (length(rc)) paste(rc, collapse = "\n") else "PLINK command failed."
          stop(
            "PLINK command failed.\n", msg,
            "\nCommand: ", cmd_line,
            call. = FALSE
          )
        }
      }
      invisible(TRUE)
    }
    .exec_plink(args)
    het_file <- paste0(out_pref, ".het")
    miss_file <- if (file.exists(paste0(out_pref, ".imiss"))) paste0(out_pref, ".imiss") else paste0(out_pref, ".smiss")
    if (!file.exists(het_file) || !file.exists(miss_file)) {
      stop("PLINK finished but expected .het/.imiss(.smiss) not found.", call. = FALSE)
    }
    list(het_file = het_file, miss_file = miss_file, out_pref = out_pref, keep_plink = keep_plink)
  }

  pv_in <- .hetmiss_normalize_version(plink.version)
  dt_het <- .hetmiss_read_any_dt(het)
  dt_miss <- .hetmiss_read_any_dt(miss)

  if (is.null(dt_het) || is.null(dt_miss)) {
    rr <- .hetmiss_run_plink(bfile0 = bfile, pfile0 = pfile, plink0 = plink, ver0 = pv_in, plink_out0 = plink.out)
    dt_het <- .hetmiss_read_any_dt(rr$het_file)
    dt_miss <- .hetmiss_read_any_dt(rr$miss_file)
    if (!isTRUE(rr$keep_plink)) {
      rm_files <- c(
        paste0(rr$out_pref, ".het"),
        paste0(rr$out_pref, ".imiss"),
        paste0(rr$out_pref, ".smiss"),
        paste0(rr$out_pref, ".log"),
        paste0(rr$out_pref, ".nosex")
      )
      rm_files <- rm_files[file.exists(rm_files)]
      if (length(rm_files)) unlink(rm_files, force = TRUE)
    }
  }
  if ((is.null(het) || is.null(miss)) && isTRUE(show.plinklog) &&
      exists("rr", inherits = FALSE) && !isTRUE(rr$keep_plink)) {
    .gcanvas_note(
      "gcanvas::hetmiss",
      "PLINK intermediate files were removed (plink.out=FALSE). Set plink.out=TRUE or a path to keep them.",
      silent = silent
    )
  }

  pv_auto <- .hetmiss_detect_plink_version(dt_het, dt_miss)
  pv <- if (identical(pv_in, "auto")) pv_auto else pv_in
  if (identical(pv, "auto")) pv <- "plink2"

  het0 <- .hetmiss_normalize_het(dt_het)
  miss0 <- .hetmiss_normalize_miss(dt_miss)
  if ("FID" %in% names(het0) && !all(is.na(het0$FID)) && ("FID" %in% names(miss0)) && !all(is.na(miss0$FID))) {
    dt <- merge(miss0, het0, by = c("FID", "IID"), all = FALSE, sort = FALSE)
  } else {
    dt <- merge(miss0, het0, by = "IID", all = FALSE, sort = FALSE)
    if (!("FID" %in% names(dt))) dt[, FID := NA_character_]
  }
  dt <- dt[is.finite(F_MISS) & is.finite(het_rate)]
  if (!nrow(dt)) stop("No overlapping valid samples after merging het/miss.", call. = FALSE)

  allmiss <- dt[is.finite(N_MISS) & is.finite(N_GENO) & N_GENO > 0 & N_MISS >= N_GENO]
  if (nrow(allmiss)) {
    dt <- dt[!(IID %in% allmiss$IID)]
  }
  if (!nrow(dt)) stop("All samples removed after full-missing filter.", call. = FALSE)

  dt_group <- .hetmiss_read_group_dt(group)
  if (!is.null(dt_group) && nrow(dt_group)) {
    dt <- merge(dt, dt_group, by = "IID", all.x = TRUE, sort = FALSE)
  }
  has_group <- "Group" %in% names(dt)
  if (has_group) {
    dt[, Group := as.character(Group)]
    dt[is.na(Group) | !nzchar(Group), Group := "Unknown"]
    grp_levels <- unique(dt$Group)
    grp_levels <- grp_levels[order(grp_levels)]
    dt[, Group := factor(Group, levels = grp_levels)]
  }

  m_miss <- mean(dt$F_MISS, na.rm = TRUE)
  sd_miss <- stats::sd(dt$F_MISS, na.rm = TRUE)
  m_het <- mean(dt$het_rate, na.rm = TRUE)
  sd_het <- stats::sd(dt$het_rate, na.rm = TRUE)
  if (!is.finite(sd_miss) || is.na(sd_miss)) sd_miss <- 0
  if (!is.finite(sd_het) || is.na(sd_het)) sd_het <- 0
  miss_all_zero <- all(is.finite(dt$F_MISS) & !is.na(dt$F_MISS) & dt$F_MISS == 0)

  miss_3 <- dt[F_MISS > (m_miss + 3 * sd_miss), unique(IID)]
  miss_5 <- dt[F_MISS > (m_miss + 5 * sd_miss), unique(IID)]
  het_3 <- dt[het_rate < (m_het - 3 * sd_het) | het_rate > (m_het + 3 * sd_het), unique(IID)]
  het_5 <- dt[het_rate < (m_het - 5 * sd_het) | het_rate > (m_het + 5 * sd_het), unique(IID)]

  line_sd <- as_num(line.sd)
  line_sd <- line_sd[is.finite(line_sd) & !is.na(line_sd) & line_sd > 0]
  if (!length(line_sd)) line_sd <- c(5, 3)
  line_sd <- unique(line_sd)
  n_line <- length(line_sd)
  line_col <- as.character(line.sd.color)
  line_col <- line_col[!is.na(line_col) & nzchar(line_col)]
  if (!length(line_col)) line_col <- c("#EF476F", "#258AB2")
  line_col <- rep(line_col, length.out = n_line)
  line_typ <- as.character(line.sd.type)
  line_typ <- line_typ[!is.na(line_typ) & nzchar(line_typ)]
  if (!length(line_typ)) line_typ <- "dashed"
  line_typ <- rep(line_typ, length.out = n_line)
  line_lwd <- as_num(line.sd.linewidth)
  line_lwd <- line_lwd[is.finite(line_lwd) & !is.na(line_lwd) & line_lwd > 0]
  if (!length(line_lwd)) line_lwd <- 0.2
  line_lwd <- rep(line_lwd, length.out = n_line)

  p <- ggplot2::ggplot(dt, ggplot2::aes(x = het_rate, y = F_MISS))
  if (has_group) {
    p <- p + ggplot2::geom_point(ggplot2::aes(color = Group), size = 0.8, alpha = 0.6)
    if (!is.null(group.color) && length(group.color)) {
      cols <- group.color
      if (is.list(cols) && !is.data.frame(cols) && !data.table::is.data.table(cols)) cols <- unlist(cols, use.names = TRUE)
      cols <- as.character(cols)
      cols <- cols[!is.na(cols) & nzchar(cols)]
      if (length(cols)) {
        lev <- levels(dt$Group)
        if (!is.null(names(cols)) && any(nzchar(names(cols)))) {
          map <- cols
        } else {
          map <- stats::setNames(rep(cols, length.out = length(lev)), lev)
        }
        p <- p + ggplot2::scale_color_manual(values = map, drop = FALSE)
      }
    }
  } else {
    p <- p + ggplot2::geom_point(size = 0.8, color = "grey20", alpha = 0.6)
  }
  p <- p +
    ggplot2::xlab("Heterozygosity") +
    ggplot2::ylab("Missingness") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = 14),
      axis.text = ggplot2::element_text(size = 12),
      legend.title = ggplot2::element_text(size = 12),
      legend.text = ggplot2::element_text(size = 12)
    )
  if (isTRUE(miss_all_zero)) {
    p <- p + ggplot2::scale_y_continuous(limits = c(0, 0.05))
  } else {
    p <- p + ggplot2::scale_y_continuous(limits = c(0, NA_real_))
  }
  for (i in seq_len(n_line)) {
    s <- line_sd[i]
    p <- p +
      ggplot2::geom_hline(
        yintercept = m_miss + s * sd_miss,
        color = line_col[i],
        linetype = line_typ[i],
        linewidth = line_lwd[i]
      ) +
      ggplot2::geom_vline(
        xintercept = c(m_het - s * sd_het, m_het + s * sd_het),
        color = line_col[i],
        linetype = line_typ[i],
        linewidth = line_lwd[i]
      )
  }

  if (!show.legend || !has_group) {
    p <- p + ggplot2::theme(legend.position = "none")
  } else if (legend.position == "bottom") {
    p <- p + ggplot2::theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.box.just = "center"
    )
  } else if (legend.position == "top") {
    p <- p + ggplot2::theme(
      legend.position = "top",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.box.just = "center"
    )
  } else if (legend.position == "left.outside") {
    p <- p + ggplot2::theme(legend.position = "left")
  } else if (legend.position == "right.outside") {
    p <- p + ggplot2::theme(legend.position = "right")
  } else if (legend.position == "left") {
    p <- p + ggplot2::theme(
      legend.position = c(0.02, 0.98),
      legend.justification = c(0, 1),
      legend.background = ggplot2::element_rect(fill = scales::alpha("white", 0.85), color = "grey20"),
      legend.box.background = ggplot2::element_rect(fill = scales::alpha("white", 0.85), color = "grey20")
    )
  } else if (legend.position == "right") {
    p <- p + ggplot2::theme(
      legend.position = c(0.98, 0.98),
      legend.justification = c(1, 1),
      legend.background = ggplot2::element_rect(fill = scales::alpha("white", 0.85), color = "grey20"),
      legend.box.background = ggplot2::element_rect(fill = scales::alpha("white", 0.85), color = "grey20")
    )
  }

  attr(p, "gcanvas_meta") <- list(
    type = "hetmiss",
    plink_version = pv,
    n_samples = as_int(nrow(dt)),
    n_all_missing_removed = as_int(nrow(allmiss)),
    thresholds = list(
      miss_mean = m_miss, miss_sd = sd_miss,
      het_mean = m_het, het_sd = sd_het
    )
  )
  .gcanvas_note(
    "gcanvas::hetmiss",
    sprintf("Done: n_samples=%d | plink.version=%s | miss.3sd=%d | het.3sd=%d | miss.5sd=%d | het.5sd=%d",
            as_int(nrow(dt)), pv, as_int(length(miss_3)), as_int(length(het_3)), as_int(length(miss_5)), as_int(length(het_5))),
    silent = silent
  )

  list(
    plot = p,
    data = data.table::copy(dt),
    outlier = list(
      het.3sd = as.character(het_3),
      het.5sd = as.character(het_5),
      miss.3sd = as.character(miss_3),
      miss.5sd = as.character(miss_5)
    )
  )
}
