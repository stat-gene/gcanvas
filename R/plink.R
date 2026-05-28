# PLINK subset extraction: variants / samples / filters into a new bfile/pfile
# with sensible defaults.

#' Extract a PLINK sub-dataset
#'
#' Convenience wrapper around `plink2` for the common subset-and-filter
#' workflow: keep/remove samples, extract/exclude variants, apply MAF / geno /
#' HWE thresholds, and write a new bfile/pfile.
#'
#' @param bfile,pfile Input PLINK1 bfile prefix or PLINK2 pfile prefix.
#' @param keep,remove Sample-id files (one per line).
#' @param extract,exclude Variant-id files (one per line).
#' @param geno,maf,hwe Numeric per-variant QC thresholds.
#' @param mind Per-sample missingness threshold.
#' @param plink Path to the PLINK binary (`plink2` by default).
#' @param out Output prefix.
#' @param out.format Output file format (`"pgen"`, `"bed"`, ...).
#' @param vst Logical. Use PLINK's variant-set table operations where applicable.
#' @param silent Logical. Suppress progress notes.
#'
#' @return Invisibly, a list with the input/output paths and the resolved
#'   PLINK command line.
#' @export
plink.extract <- function(bfile = NULL,
                          pfile = NULL,
                          keep = NULL,
                          remove = NULL,
                          extract = NULL,
                          exclude = NULL,
                          geno = NULL,
                          maf = NULL,
                          hwe = NULL,
                          mind = NULL,
                          plink = "plink2",
                          out = "gcanvas.extracted",
                          out.format = "pgen",
                          vst = TRUE,
                          silent = FALSE) {
  require_pkg("data.table")
  silent <- isTRUE(silent)
  vst <- isTRUE(vst)

  .plink_extract_strip_quotes <- function(s) {
    s <- trimws(as.character(s))
    if (!nzchar(s)) return(s)
    gsub("^[\"']|[\"']$", "", s)
  }
  .plink_extract_alias_target_from_rc <- function(cmd, rc_files = c("~/.zshrc", "~/.zprofile", "~/.bashrc", "~/.bash_profile")) {
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
      rhs <- .plink_extract_strip_quotes(rhs)
      if (!nzchar(rhs)) next
      tok <- strsplit(rhs, "\\s+")[[1]]
      if (!length(tok)) next
      bin <- .plink_extract_strip_quotes(tok[1])
      if (!nzchar(bin)) next
      w <- Sys.which(bin)
      if (nzchar(w)) return(w)
      if (file.exists(path.expand(bin))) return(abs_path(path.expand(bin)))
    }
    NA_character_
  }
  .plink_extract_resolve_exec <- function(cmd) {
    cmd <- as.character(cmd)[1]
    if (is.na(cmd) || !nzchar(cmd)) return(NA_character_)
    if (file.exists(cmd)) return(abs_path(cmd))
    w <- Sys.which(cmd)
    if (nzchar(w)) return(w)
    cand_dir <- c(path.expand("~/bin"), "/opt/homebrew/bin", "/usr/local/bin", "/usr/bin", "/bin")
    cand <- file.path(cand_dir, cmd)
    hit <- cand[file.exists(cand)]
    if (length(hit)) return(abs_path(hit[1]))
    a <- .plink_extract_alias_target_from_rc(cmd)
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
  .plink_extract_resolve_plink1_exec <- function(plink2_exec) {
    cand <- unique(c(
      file.path(dirname(plink2_exec), "plink"),
      "/Users/bsu/software/plink/plink",
      Sys.which("plink")
    ))
    cand <- cand[!is.na(cand) & nzchar(cand)]
    cand <- cand[file.exists(cand) | nzchar(Sys.which(cand))]
    for (cc in cand) {
      exe <- if (file.exists(cc)) abs_path(cc) else Sys.which(cc)
      if (!nzchar(exe)) next
      vv <- tryCatch(suppressWarnings(system2(exe, "--version", stdout = TRUE, stderr = TRUE)), error = function(e) character(0))
      txt <- tolower(paste(vv, collapse = " "))
      if (grepl("plink v1\\.9", txt, fixed = FALSE)) return(exe)
    }
    NA_character_
  }
  .plink_extract_run_plink <- function(exec, args) {
    cmd_line <- paste(c(shQuote(exec), vapply(args, shQuote, character(1))), collapse = " ")
    args_run <- gsub("\\$", "\\\\\\$", as.character(args))
    rc <- tryCatch(
      suppressWarnings(system2(command = exec, args = args_run, stdout = TRUE, stderr = TRUE)),
      error = function(e) stop(sprintf("Failed to run PLINK command: %s\nCommand: %s", conditionMessage(e), cmd_line), call. = FALSE)
    )
    status <- attr(rc, "status") %||% 0L
    if (!isTRUE(silent) && length(rc)) cat(paste0(rc, collapse = "\n"), "\n")
    if (!identical(as_int(status), 0L)) {
      msg <- if (length(rc)) paste(rc, collapse = "\n") else "PLINK command failed."
      stop("PLINK command failed.\n", msg, "\nCommand: ", cmd_line, call. = FALSE)
    }
    invisible(TRUE)
  }
  .plink_extract_norm_prefix <- function(x) {
    y <- as.character(x)
    y <- y[!is.na(y) & nzchar(y)]
    if (!length(y)) return(character())
    y <- trimws(y)
    y <- sub("\\.(bed|bim|fam|pgen|psam|pvar|pvar\\.zst)$", "", y, ignore.case = TRUE)
    unique(abs_path(path.expand(y)))
  }
  .plink_extract_dataset_exists <- function(prefix) {
    is_b <- all(file.exists(paste0(prefix, c(".bed", ".bim", ".fam"))))
    is_p <- file.exists(paste0(prefix, ".pgen")) &&
      file.exists(paste0(prefix, ".psam")) &&
      (file.exists(paste0(prefix, ".pvar")) || file.exists(paste0(prefix, ".pvar.zst")))
    list(bfile = is_b, pfile = is_p)
  }
  .plink_extract_dataset_args <- function(ds) {
    if (identical(ds$format, "bfile")) return(c("--bfile", ds$prefix))
    c("--pfile", ds$prefix, if (isTRUE(ds$use_vzs)) "vzs")
  }
  .plink_extract_pick_col <- function(nm, keys) {
    cn <- toupper(gsub("[^A-Z0-9#]+", "", as.character(nm)))
    for (k in keys) {
      kk <- toupper(gsub("[^A-Z0-9#]+", "", as.character(k)))
      hit <- which(cn == kk)
      if (length(hit)) return(nm[hit[1]])
    }
    NA_character_
  }
  .plink_extract_prepare_sample_dt <- function(x) {
    if (is.null(x) || length(x) == 0L) return(NULL)
    if (data.table::is.data.table(x) || is.data.frame(x)) {
      dt <- if (data.table::is.data.table(x)) data.table::copy(x) else data.table::as.data.table(x)
      if (!nrow(dt)) return(NULL)
      c_fid <- .plink_extract_pick_col(names(dt), c("FID", "#FID", "family.ID", "familyID", "FAMILY_ID"))
      c_iid <- .plink_extract_pick_col(names(dt), c("IID", "#IID", "ID", "sample.ID", "sampleID", "INDV", "INDIVIDUAL"))
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
        out[, FID := IID]
        out <- out[, .(FID, IID)]
      } else {
        out <- dt[, .(FID = as.character(get(c_fid)), IID = as.character(get(c_iid)))]
        out <- out[!is.na(IID) & nzchar(IID)]
        out[is.na(FID) | !nzchar(FID), FID := IID]
        out <- unique(out, by = c("FID", "IID"))
      }
      out
    } else {
      iid <- as.character(x)
      iid <- iid[!is.na(iid) & nzchar(iid)]
      if (!length(iid)) return(NULL)
      iid <- unique(iid)
      data.table::data.table(IID = iid)
    }
  }
  .plink_extract_load_psam_map <- function(prefix) {
    psam <- paste0(prefix, ".psam")
    if (!file.exists(psam)) return(NULL)
    dts <- data.table::fread(psam, data.table = TRUE, check.names = FALSE, showProgress = FALSE)
    if (!nrow(dts)) return(NULL)
    c_iid <- .plink_extract_pick_col(names(dts), c("#IID", "IID", "ID", "sample.ID", "sampleID"))
    c_fid <- .plink_extract_pick_col(names(dts), c("#FID", "FID", "family.ID", "familyID"))
    if (is.na(c_iid)) return(NULL)
    out <- if (is.na(c_fid)) {
      data.table::data.table(IID = as.character(dts[[c_iid]]))
    } else {
      data.table::data.table(FID = as.character(dts[[c_fid]]), IID = as.character(dts[[c_iid]]))
    }
    out <- out[!is.na(IID) & nzchar(IID)]
    if (!("FID" %in% names(out))) out[, FID := IID]
    out[is.na(FID) | !nzchar(FID), FID := IID]
    unique(out[, .(FID, IID)], by = c("FID", "IID"))
  }
  .plink_extract_align_sample_filter <- function(x_dt, sample_map, label = "keep", strict = TRUE) {
    if (is.null(x_dt) || !nrow(x_dt)) return(NULL)
    if (is.null(sample_map) || !nrow(sample_map)) {
      if (!("FID" %in% names(x_dt))) x_dt[, FID := IID]
      x_dt <- x_dt[!is.na(IID) & nzchar(IID)]
      x_dt[is.na(FID) | !nzchar(FID), FID := IID]
      return(unique(x_dt[, .(FID, IID)], by = c("FID", "IID")))
    }
    if (!("FID" %in% names(x_dt))) {
      y <- merge(
        x_dt[, .(IID = as.character(IID))],
        sample_map[, .(FID, IID)],
        by = "IID",
        all.x = FALSE,
        all.y = FALSE,
        sort = FALSE
      )
      y <- unique(y[, .(FID, IID)], by = c("FID", "IID"))
    } else {
      y <- x_dt[, .(FID = as.character(FID), IID = as.character(IID))]
      y <- y[!is.na(IID) & nzchar(IID)]
      y[is.na(FID) | !nzchar(FID), FID := IID]
      data.table::setkey(sample_map, FID, IID)
      data.table::setkey(y, FID, IID)
      y <- sample_map[y, nomatch = 0L]
      y <- unique(y[, .(FID, IID)], by = c("FID", "IID"))
    }
    if (!nrow(y) && isTRUE(strict)) {
      stop(sprintf("%s has no matching samples in merged input (FID/IID after matching).", label), call. = FALSE)
    }
    y
  }
  .plink_extract_prepare_variant_dt <- function(x) {
    if (is.null(x) || length(x) == 0L) return(NULL)
    if (data.table::is.data.table(x) || is.data.frame(x)) {
      dt <- if (data.table::is.data.table(x)) data.table::copy(x) else data.table::as.data.table(x)
      if (!nrow(dt)) return(NULL)
      c_id <- .plink_extract_pick_col(names(dt), c("SNP", "ID", "marker.ID", "VARIANT", "RSID"))
      if (is.na(c_id)) c_id <- names(dt)[1]
      out <- dt[, .(ID = as.character(get(c_id)))]
      out <- out[!is.na(ID) & nzchar(ID)]
      unique(out, by = "ID")
    } else {
      id <- as.character(x)
      id <- id[!is.na(id) & nzchar(id)]
      if (!length(id)) return(NULL)
      data.table::data.table(ID = unique(id))
    }
  }
  .plink_extract_load_sample_map <- function(ds) {
    if (identical(ds$format, "bfile")) {
      fam <- paste0(ds$prefix, ".fam")
      if (!file.exists(fam)) return(NULL)
      dtf <- data.table::fread(fam, data.table = TRUE, header = FALSE, showProgress = FALSE)
      if (ncol(dtf) < 2L) return(NULL)
      out <- data.table::data.table(FID = as.character(dtf[[1]]), IID = as.character(dtf[[2]]))
      out <- out[!is.na(IID) & nzchar(IID)]
      return(unique(out, by = c("FID", "IID")))
    }
    .plink_extract_load_psam_map(ds$prefix)
  }
  .plink_extract_load_variant_ids <- function(ds) {
    if (identical(ds$format, "bfile")) {
      bim <- paste0(ds$prefix, ".bim")
      if (!file.exists(bim)) return(character())
      dt <- data.table::fread(bim, data.table = TRUE, header = FALSE, showProgress = FALSE)
      if (ncol(dt) < 2L) return(character())
      ids <- as.character(dt[[2]])
    } else {
      pvar <- if (file.exists(paste0(ds$prefix, ".pvar"))) paste0(ds$prefix, ".pvar") else paste0(ds$prefix, ".pvar.zst")
      if (!file.exists(pvar)) return(character())
      dt <- suppressWarnings(data.table::fread(pvar, data.table = TRUE, check.names = FALSE, showProgress = FALSE, quote = ""))
      c_id <- .plink_extract_pick_col(names(dt), c("ID", "SNP", "VAR", "marker.ID", "RSID"))
      if (is.na(c_id)) return(character())
      ids <- as.character(dt[[c_id]])
    }
    ids <- ids[!is.na(ids) & nzchar(ids)]
    unique(ids)
  }

  fmt_input <- as.character(out.format)[1]
  fmt <- tolower(trimws(fmt_input))
  if (is.na(fmt) || !nzchar(fmt)) fmt <- "pgen"
  out_mode <- if (fmt %in% c("p", "pgen", "pfile")) {
    "pgen"
  } else if (fmt %in% c("b", "bed", "bfile")) {
    "bed"
  } else {
    stop("out.format must be one of: p/pgen/pfile or b/bed/bfile.", call. = FALSE)
  }
  geno0 <- as_num(geno)[1]
  use_geno <- is.finite(geno0) && !is.na(geno0) && geno0 >= 0 && geno0 <= 1
  maf0 <- as_num(maf)[1]
  use_maf <- is.finite(maf0) && !is.na(maf0) && maf0 > 0 && maf0 < 1
  hwe_in <- as.character(hwe)[1]
  use_hwe <- !is.na(hwe_in) && nzchar(trimws(hwe_in)) && !(tolower(trimws(hwe_in)) %in% c("false", "null", "na"))
  if (use_hwe) {
    hwe_num <- suppressWarnings(as.numeric(hwe_in))
    if (is.finite(hwe_num) && !is.na(hwe_num)) {
      if (hwe_num <= 0 || hwe_num >= 1) use_hwe <- FALSE else hwe_in <- format(hwe_num, scientific = TRUE)
    }
  }
  mind0 <- as_num(mind)[1]
  use_mind <- is.finite(mind0) && !is.na(mind0) && mind0 >= 0 && mind0 <= 1

  bvec <- .plink_extract_norm_prefix(bfile)
  pvec <- .plink_extract_norm_prefix(pfile)
  if (!length(bvec) && !length(pvec)) stop("Provide at least one bfile or pfile prefix.", call. = FALSE)

  inputs <- list()
  for (x in bvec) inputs[[length(inputs) + 1L]] <- list(prefix = x, format = "bfile", use_vzs = FALSE)
  for (x in pvec) {
    use_vzs <- file.exists(paste0(x, ".pvar.zst")) && !file.exists(paste0(x, ".pvar"))
    inputs[[length(inputs) + 1L]] <- list(prefix = x, format = "pfile", use_vzs = use_vzs)
  }
  if (!length(inputs)) stop("No valid inputs detected.", call. = FALSE)

  for (i in seq_along(inputs)) {
    ex <- .plink_extract_dataset_exists(inputs[[i]]$prefix)
    if (identical(inputs[[i]]$format, "bfile") && !isTRUE(ex$bfile)) {
      stop("bfile not found: ", inputs[[i]]$prefix, call. = FALSE)
    }
    if (identical(inputs[[i]]$format, "pfile") && !isTRUE(ex$pfile)) {
      stop("pfile not found: ", inputs[[i]]$prefix, call. = FALSE)
    }
  }

  plink_cmd <- as.character(plink)[1]
  if (is.na(plink_cmd) || !nzchar(plink_cmd)) stop("PLINK command is empty.", call. = FALSE)
  plink_exec <- .plink_extract_resolve_exec(plink_cmd)
  if (is.na(plink_exec) || !nzchar(plink_exec)) {
    stop(sprintf("PLINK executable not found for '%s'. Set `plink` to an executable path.", plink_cmd), call. = FALSE)
  }

  out_prefix <- abs_path(path.expand(as.character(out)[1]))
  out_prefix <- sub("\\.(bed|bim|fam|pgen|psam|pvar|pvar\\.zst)$", "", out_prefix, ignore.case = TRUE)
  dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)

  cache_dir <- .gcanvas_default_cache_dir(scope = "plink.extract", anchor = out_prefix)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  .gcanvas_register_cache_dir(cache_dir)
  work_pref <- file.path(cache_dir, sprintf("plink_extract_%s_%d", format(Sys.time(), "%Y%m%d_%H%M%S"), as_int(stats::runif(1, 1, 1e9))))
  dir.create(work_pref, recursive = TRUE, showWarnings = FALSE)
  .gcanvas_register_cache_dir(work_pref)

  .gcanvas_note(
    "gcanvas::plink.extract",
    sprintf("Start: n_input=%d | out=%s | format=%s | vst=%s | geno=%s | maf=%s | hwe=%s | mind=%s",
            as_int(length(inputs)), out_prefix, out_mode, ifelse(isTRUE(vst), "TRUE", "FALSE"),
            if (use_geno) format(geno0, scientific = FALSE) else "NULL",
            if (use_maf) format(maf0, scientific = FALSE) else "NULL",
            if (use_hwe) hwe_in else "NULL",
            if (use_mind) format(mind0, scientific = FALSE) else "NULL"),
    silent = silent
  )

  keep_dt0 <- .plink_extract_prepare_sample_dt(keep)
  remove_dt0 <- .plink_extract_prepare_sample_dt(remove)
  extract_dt0 <- .plink_extract_prepare_variant_dt(extract)
  exclude_dt0 <- .plink_extract_prepare_variant_dt(exclude)

  # Pass 1: sample/variant availability scan per input.
  scan_list <- vector("list", length(inputs))
  n_skipped <- 0L
  for (i in seq_along(inputs)) {
    dsi <- inputs[[i]]
    smap_i <- .plink_extract_load_sample_map(dsi)
    keep_i <- .plink_extract_align_sample_filter(keep_dt0, smap_i, label = sprintf("keep (input %d)", i), strict = FALSE)
    remove_i <- .plink_extract_align_sample_filter(remove_dt0, smap_i, label = sprintf("remove (input %d)", i), strict = FALSE)
    if (!is.null(keep_i) && nrow(keep_i) == 0L) {
      n_skipped <- n_skipped + 1L
      next
    }
    if (!is.null(remove_i) && nrow(remove_i) && !is.null(smap_i) && nrow(smap_i)) {
      remain_iid <- setdiff(unique(smap_i$IID), unique(remove_i$IID))
      if (length(remain_iid) == 0L) {
        n_skipped <- n_skipped + 1L
        next
      }
    }
    var_ids <- .plink_extract_load_variant_ids(dsi)
    if (!length(var_ids)) {
      n_skipped <- n_skipped + 1L
      next
    }
    keep_var_i <- var_ids
    if (!is.null(extract_dt0) && nrow(extract_dt0)) keep_var_i <- intersect(keep_var_i, as.character(extract_dt0$ID))
    if (!is.null(exclude_dt0) && nrow(exclude_dt0)) keep_var_i <- setdiff(keep_var_i, as.character(exclude_dt0$ID))
    if (!length(keep_var_i)) {
      n_skipped <- n_skipped + 1L
      next
    }
    scan_list[[i]] <- list(ds = dsi, keep = keep_i, remove = remove_i, keep_var = unique(keep_var_i))
  }
  scan_list <- Filter(Negate(is.null), scan_list)
  if (!length(scan_list)) {
    stop("No dataset remained after keep/remove/extract/exclude pre-filtering.", call. = FALSE)
  }
  # Force concatenating merge in plink2 by extracting common variant set across active inputs.
  shared_var <- Reduce(intersect, lapply(scan_list, `[[`, "keep_var"))
  shared_var <- unique(as.character(shared_var))
  shared_var <- shared_var[!is.na(shared_var) & nzchar(shared_var)]
  if (!length(shared_var)) {
    stop("No common variants remained across filtered inputs; cannot merge.", call. = FALSE)
  }
  shared_extract_file <- file.path(work_pref, "extract_shared.txt")
  data.table::fwrite(data.table::data.table(ID = shared_var), shared_extract_file, sep = "\t", col.names = FALSE, quote = FALSE)

  # Pass 2: apply filters and write per-input bed with common variant scaffold.
  filtered_beds <- character()
  for (i in seq_along(scan_list)) {
    si <- scan_list[[i]]
    keep_file_i <- NA_character_
    remove_file_i <- NA_character_
    if (!is.null(si$keep) && nrow(si$keep)) {
      keep_file_i <- file.path(work_pref, sprintf("keep_%03d.txt", as_int(i)))
      data.table::fwrite(si$keep[, .(FID, IID)], keep_file_i, sep = "\t", col.names = FALSE, quote = FALSE)
    }
    if (!is.null(si$remove) && nrow(si$remove)) {
      remove_file_i <- file.path(work_pref, sprintf("remove_%03d.txt", as_int(i)))
      data.table::fwrite(si$remove[, .(FID, IID)], remove_file_i, sep = "\t", col.names = FALSE, quote = FALSE)
    }
    bed_i <- file.path(work_pref, sprintf("flt_%03d", as_int(i)))
    args_i <- c(.plink_extract_dataset_args(si$ds))
    if (!is.na(keep_file_i) && nzchar(keep_file_i) && file.exists(keep_file_i)) args_i <- c(args_i, "--keep", keep_file_i)
    if (!is.na(remove_file_i) && nzchar(remove_file_i) && file.exists(remove_file_i)) args_i <- c(args_i, "--remove", remove_file_i)
    args_i <- c(args_i, "--extract", shared_extract_file, "--make-bed", "--out", bed_i)
    .plink_extract_run_plink(plink_exec, args_i)
    if (all(file.exists(paste0(bed_i, c(".bed", ".bim", ".fam"))))) filtered_beds <- c(filtered_beds, bed_i) else n_skipped <- n_skipped + 1L
  }
  filtered_beds <- unique(filtered_beds)
  if (!length(filtered_beds)) stop("No filtered bed files produced.", call. = FALSE)

  merged_bed <- filtered_beds[1]
  if (length(filtered_beds) > 1L) {
    plink1_exec <- .plink_extract_resolve_plink1_exec(plink_exec)
    if (is.na(plink1_exec) || !nzchar(plink1_exec)) {
      stop("Multiple-input merge requires PLINK 1.9 (plink). Could not find a PLINK 1.9 executable.", call. = FALSE)
    }
    merge_list <- file.path(work_pref, "merge_list.txt")
    data.table::fwrite(data.table::data.table(prefix = filtered_beds), merge_list, sep = "\t", col.names = FALSE, quote = FALSE)
    merged_bed <- file.path(work_pref, "merged")
    args_merge <- c(
      "--merge-list", merge_list,
      "--make-bed",
      "--out", merged_bed
    )
    .plink_extract_run_plink(plink1_exec, args_merge)
  }

  if (identical(out_mode, "pgen")) {
    args_final <- c("--bfile", merged_bed)
    if (use_geno) args_final <- c(args_final, "--geno", format(geno0, scientific = FALSE))
    if (use_maf) args_final <- c(args_final, "--maf", format(maf0, scientific = FALSE))
    if (use_hwe) args_final <- c(args_final, "--hwe", hwe_in)
    if (use_mind) args_final <- c(args_final, "--mind", format(mind0, scientific = FALSE))
    args_final <- c(args_final, "--make-pgen", if (isTRUE(vst)) "vzs", "--out", out_prefix)
  } else {
    args_final <- c("--bfile", merged_bed)
    if (use_geno) args_final <- c(args_final, "--geno", format(geno0, scientific = FALSE))
    if (use_maf) args_final <- c(args_final, "--maf", format(maf0, scientific = FALSE))
    if (use_hwe) args_final <- c(args_final, "--hwe", hwe_in)
    if (use_mind) args_final <- c(args_final, "--mind", format(mind0, scientific = FALSE))
    args_final <- c(args_final, "--make-bed", "--out", out_prefix)
  }
  .plink_extract_run_plink(plink_exec, args_final)

  out_files <- if (identical(out_mode, "pgen")) {
    c(paste0(out_prefix, ".pgen"), paste0(out_prefix, ".psam"),
      if (isTRUE(vst)) paste0(out_prefix, ".pvar.zst") else paste0(out_prefix, ".pvar"))
  } else {
    c(paste0(out_prefix, ".bed"), paste0(out_prefix, ".bim"), paste0(out_prefix, ".fam"))
  }

  .gcanvas_note(
    "gcanvas::plink.extract",
    sprintf("Done: out=%s | format=%s | n_input=%d | used=%d | skipped=%d",
            out_prefix, out_mode, as_int(length(inputs)), as_int(length(filtered_beds)), as_int(n_skipped)),
    silent = silent
  )

  list(
    out = out_prefix,
    format = out_mode,
    file = out_files[file.exists(out_files)]
  )
}

