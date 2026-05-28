# LD reference attachment (bigsnpr backend), r/r-squared computation, proxy
# lookup, and clumping. Exports `ldproxy()` and `calcld()`. (`ldclump()` lives
# in this file too, appended at the end.)

canon_a1 <- function(a1, a2) {
  a1 <- toupper(as.character(a1)); a2 <- toupper(as.character(a2))
  ifelse(a1 <= a2, a1, a2)
}
canon_a2 <- function(a1, a2) {
  a1 <- toupper(as.character(a1)); a2 <- toupper(as.character(a2))
  ifelse(a1 <= a2, a2, a1)
}

.gcanvas_ld_with_file_lock <- function(lock_path, expr, timeout_sec = 300, poll_sec = 0.2) {
  lock_path <- as.character(lock_path)[1]
  t0 <- Sys.time()
  repeat {
    ok <- FALSE
    try({ ok <- file.create(lock_path) }, silent = TRUE)
    if (isTRUE(ok)) break
    if (as.numeric(difftime(Sys.time(), t0, units = "secs")) > timeout_sec) {
      stop("Timeout acquiring lock: ", lock_path, call. = FALSE)
    }
    Sys.sleep(poll_sec)
  }
  on.exit(unlink(lock_path), add = TRUE)
  force(expr)
}

.gcanvas_ld_thread_lock <- function(n = 1L) {
  n <- as_int(n); if (is.na(n) || n < 1L) n <- 1L
  st <- list(
    env = Sys.getenv(c(
      "OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS",
      "VECLIB_MAXIMUM_THREADS","BLIS_NUM_THREADS","NUMEXPR_NUM_THREADS"
    ), unset = NA_character_),
    blas = NA_integer_,
    omp  = NA_integer_,
    have_rhpc = requireNamespace("RhpcBLASctl", quietly = TRUE)
  )
  if (isTRUE(st$have_rhpc)) {
    st$blas <- suppressWarnings(tryCatch(RhpcBLASctl::blas_get_num_procs(), error = function(e) NA_integer_))
    st$omp  <- suppressWarnings(tryCatch(RhpcBLASctl::omp_get_num_procs(),  error = function(e) NA_integer_))
    suppressWarnings(try(RhpcBLASctl::blas_set_num_threads(n), silent = TRUE))
    suppressWarnings(try(RhpcBLASctl::omp_set_num_threads(n),  silent = TRUE))
  }
  Sys.setenv(
    OMP_NUM_THREADS = as.character(n),
    OPENBLAS_NUM_THREADS = as.character(n),
    MKL_NUM_THREADS = as.character(n),
    VECLIB_MAXIMUM_THREADS = as.character(n),
    BLIS_NUM_THREADS = as.character(n),
    NUMEXPR_NUM_THREADS = as.character(n)
  )
  st
}

.gcanvas_max_threads <- function() {
  n <- suppressWarnings(tryCatch(parallel::detectCores(logical = TRUE), error = function(e) NA_integer_))
  if (!is.finite(n) || is.na(n) || n < 1L) n <- 1L
  as_int(n)
}

.gcanvas_resolve_threads <- function(threads = 1L) {
  req <- as_int(threads)
  n <- req
  if (!is.finite(n) || is.na(n) || n < 1L) n <- 1L
  nmax <- .gcanvas_max_threads()
  if (!is.finite(nmax) || is.na(nmax) || nmax < 1L) nmax <- 1L
  if (n > nmax) n <- nmax
  req_show <- if (is.finite(req) && !is.na(req)) as_int(req) else 1L
  if (is.finite(req_show) && !is.na(req_show) && req_show != as_int(n)) {
    .gcanvas_note("gcanvas", sprintf("threads requested=%d, using=%d", as_int(req_show), as_int(n)))
  }
  as_int(n)
}

.gcanvas_ld_thread_restore <- function(st) {
  if (!is.list(st)) return(invisible(FALSE))
  if (isTRUE(st$have_rhpc) && requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    if (is.finite(st$blas)) suppressWarnings(try(RhpcBLASctl::blas_set_num_threads(as_int(st$blas)), silent = TRUE))
    if (is.finite(st$omp))  suppressWarnings(try(RhpcBLASctl::omp_set_num_threads(as_int(st$omp)),  silent = TRUE))
  }
  env0 <- st$env
  if (is.null(env0)) return(invisible(TRUE))
  for (nm in names(env0)) {
    val <- env0[[nm]]
    if (is.na(val)) Sys.unsetenv(nm) else do.call(Sys.setenv, setNames(list(val), nm))
  }
  invisible(TRUE)
}

ldref_signature <- function(bfile, bed_subset_bytes = 1e6) {
  require_pkg("digest")
  bfile <- as.character(bfile)[1]
  bed <- paste0(bfile, ".bed")
  bim <- paste0(bfile, ".bim")
  fam <- paste0(bfile, ".fam")
  if (!file.exists(bed)) stop("Missing .bed: ", bed, call. = FALSE)
  if (!file.exists(bim)) stop("Missing .bim: ", bim, call. = FALSE)
  if (!file.exists(fam)) stop("Missing .fam: ", fam, call. = FALSE)

  h_bim <- digest::digest(file = bim, algo = "xxhash64")
  h_fam <- digest::digest(file = fam, algo = "xxhash64")
  bed_size <- as.numeric(file.info(bed)$size)
  if (!is.finite(bed_size)) bed_size <- NA_real_

  h_bed_subs <- NA_character_
  subset_n <- as_int(bed_subset_bytes)
  if (is.finite(subset_n) && subset_n > 0L) {
    con <- file(bed, "rb")
    on.exit(close(con), add = TRUE)
    head_raw <- readBin(con, what = "raw", n = subset_n)
    tail_raw <- raw(0)
    if (is.finite(bed_size) && bed_size > 0) {
      tail_n <- min(subset_n, as_int(bed_size))
      seek(con, where = max(0, bed_size - tail_n), origin = "start")
      tail_raw <- readBin(con, what = "raw", n = tail_n)
    }
    h_bed_subs <- digest::digest(c(head_raw, tail_raw), algo = "xxhash64", serialize = FALSE)
  }

  key_str <- paste0(
    "bim=", h_bim,
    "|fam=", h_fam,
    "|bed_size=", sprintf("%.0f", bed_size),
    "|bed_subset=", h_bed_subs
  )
  digest::digest(key_str, algo = "xxhash64", serialize = FALSE)
}

attach_ld_ref <- function(bfile = NULL, ld_rds = NULL, cache_dir) {
  require_pkg(c("bigsnpr", "digest"))
  cache_dir <- as.character(cache_dir)[1]
  cache_dir_was_new <- !dir.exists(cache_dir)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  if (isTRUE(cache_dir_was_new) && dir.exists(cache_dir)) .gcanvas_register_cache_dir(cache_dir)

  .gcanvas_ld_rds_stamp <- function(path) {
    p <- abs_path(path)
    fi <- file.info(p)
    size <- as.numeric(fi$size)[1]
    mt <- as.numeric(fi$mtime)[1]
    ct <- as.numeric(fi$ctime)[1]
    hs <- NA_character_
    if (file.exists(p) && is.finite(size) && !is.na(size) && size > 0) {
      con <- file(p, "rb")
      on.exit(close(con), add = TRUE)
      n <- min(262144L, as_int(size))
      head_raw <- readBin(con, what = "raw", n = n)
      mid_raw <- raw(0)
      tail_raw <- raw(0)
      mid_start <- max(0, as_int(floor(size / 2) - n / 2))
      seek(con, where = mid_start, origin = "start")
      mid_raw <- readBin(con, what = "raw", n = n)
      tail_start <- max(0, as_int(size - n))
      seek(con, where = tail_start, origin = "start")
      tail_raw <- readBin(con, what = "raw", n = n)
      hs <- digest::digest(c(head_raw, mid_raw, tail_raw), algo = "xxhash64", serialize = FALSE)
    }
    paste0("size=", sprintf("%.0f", size), "|mtime=", sprintf("%.0f", mt), "|ctime=", sprintf("%.0f", ct), "|sample3=", hs)
  }

  if (!is.null(ld_rds) && nzchar(ld_rds) && file.exists(ld_rds)) {
    .gcanvas_note("gcanvas", paste0("LD ref: snp_attach(", abs_path(ld_rds), ")"))
    obj <- bigsnpr::snp_attach(ld_rds)
    attr(obj, "gcanvas_ldref_id") <- paste0("rds:", abs_path(ld_rds), "|", .gcanvas_ld_rds_stamp(ld_rds))
    return(obj)
  }

  if (is.null(bfile) || !nzchar(bfile)) stop("Provide bfile or ld_rds.", call. = FALSE)
  bed <- paste0(bfile, ".bed")
  if (!file.exists(bed)) stop("Missing bed: ", bed, call. = FALSE)

  sig <- ldref_signature(bfile, bed_subset_bytes = 1e6)
  backing_base <- file.path(cache_dir, paste0(sig, ".ldref"))
  rds_path <- paste0(backing_base, ".rds")
  bk_path <- paste0(backing_base, ".bk")
  lk <- file.path(cache_dir, paste0(sig, ".ldref.lock"))

  .gcanvas_ld_with_file_lock(lk, {
    if (file.exists(rds_path) && file.exists(bk_path)) {
      .gcanvas_note("gcanvas", paste0("LD ref cached: ", abs_path(rds_path)))
      obj <- bigsnpr::snp_attach(rds_path)
      attr(obj, "gcanvas_ldref_id") <- paste0("sig:", sig)
      return(obj)
    }
    .gcanvas_note("gcanvas", paste0("Converting PLINK ---> bigsnpr backing: ", abs_path(backing_base)))
    bigsnpr::snp_readBed(bed, backingfile = backing_base)
    if (!file.exists(rds_path) || !file.exists(bk_path)) stop("snp_readBed failed to create .rds/.bk", call. = FALSE)
    .gcanvas_register_cache_files(c(rds_path, bk_path))
    obj <- bigsnpr::snp_attach(rds_path)
    attr(obj, "gcanvas_ldref_id") <- paste0("sig:", sig)
    obj
  })
}

.gcanvas_ld_bfile_ready <- function(bfile) {
  x <- bfile[1]
  if (is.null(x) || length(x) == 0L || is.na(x)) return(FALSE)
  bfile <- as.character(x)
  if (!nzchar(bfile)) return(FALSE)
  need <- paste0(bfile, c(".bed", ".bim", ".fam"))
  all(file.exists(need))
}

.gcanvas_normalize_bfile_prefix <- function(x) {
  if (is.null(x) || length(x) == 0L || is.na(x[1])) return(NULL)
  s <- as.character(x[1])
  if (!nzchar(s)) return(NULL)
  s <- sub("\\.(bed|bim|fam)$", "", s, ignore.case = TRUE)
  abs_path(s)
}

.gcanvas_default_cache_dir <- function(scope = "default", anchor = NULL) {
  sc <- as.character(scope)[1]
  if (is.na(sc) || !nzchar(sc)) sc <- "default"
  sc <- gsub("[^A-Za-z0-9._-]+", "_", sc)
  if (!nzchar(sc)) sc <- "default"

  x <- if (is.null(anchor)) NULL else anchor[1]
  if (!is.null(x) && !is.na(x) && nzchar(as.character(x))) {
    ap <- abs_path(as.character(x))
    if (dir.exists(ap)) {
      return(file.path(ap, ".gcanvas_cache", sc))
    }
    return(file.path(dirname(ap), ".gcanvas_cache", sc))
  }
  if ("R_user_dir" %in% getNamespaceExports("tools")) {
    return(file.path(tools::R_user_dir("gcanvas", "cache"), sc))
  }
  file.path(tempdir(), "gcanvas", sc)
}

.gcanvas_default_ld_cache_dir <- function(ld.bfile = NULL) {
  x <- if (is.null(ld.bfile)) NULL else ld.bfile[1]
  .gcanvas_default_cache_dir(scope = "ldref", anchor = x)
}

.gcanvas_session_cache_registry <- local({
  e <- new.env(parent = emptyenv())
  e$created_files <- character(0)
  e$created_dirs <- character(0)
  e
})

.gcanvas_register_cache_files <- function(paths) {
  if (is.null(paths) || !length(paths)) return(invisible(FALSE))
  p <- abs_path(paths)
  p <- unique(as.character(p))
  p <- p[!is.na(p) & nzchar(p)]
  if (!length(p)) return(invisible(FALSE))
  e <- .gcanvas_session_cache_registry
  e$created_files <- unique(c(e$created_files, p))
  e$created_dirs <- unique(c(e$created_dirs, dirname(p)))
  invisible(TRUE)
}

.gcanvas_register_cache_dir <- function(path) {
  if (is.null(path) || !length(path)) return(invisible(FALSE))
  p <- abs_path(path)[1]
  if (is.na(p) || !nzchar(p)) return(invisible(FALSE))
  e <- .gcanvas_session_cache_registry
  e$created_dirs <- unique(c(e$created_dirs, p))
  invisible(TRUE)
}

.gcanvas_ld_mem_cache <- local({
  e <- new.env(parent = emptyenv())
  e$entries <- new.env(parent = emptyenv())
  e$order <- character(0)
  e$total_bytes <- 0
  e
})

.gcanvas_ld_mem_cache_max_entries <- function() {
  n <- suppressWarnings(as.integer(getOption("gcanvas.ld.mem.cache.max_entries", 32L)))
  if (!is.finite(n) || is.na(n) || n < 0L) n <- 32L
  n
}

.gcanvas_ld_mem_cache_max_bytes <- function() {
  mb <- suppressWarnings(as.numeric(getOption("gcanvas.ld.mem.cache.max_mb", 128)))
  if (!is.finite(mb) || is.na(mb) || mb < 0) mb <- 128
  as.numeric(mb) * 1024^2
}

.gcanvas_ld_mem_cache_enabled <- function() {
  .gcanvas_ld_mem_cache_max_entries() > 0L && .gcanvas_ld_mem_cache_max_bytes() > 0
}

.gcanvas_ld_mem_cache_key <- function(ld_obj, chr, lead_pos, ridx_sorted, pos_sorted) {
  ld_id <- attr(ld_obj, "gcanvas_ldref_id", exact = TRUE)
  if (is.null(ld_id) || !nzchar(as.character(ld_id)[1])) {
    map <- ld_obj$map
    ld_id <- paste0("map:", nrow(map), ":", sum(.gcanvas_as_num2(map$physical.pos), na.rm = TRUE))
  }
  n <- length(ridx_sorted)
  ridx_sig <- if (!n) "0" else paste0(
    n, ":",
    ridx_sorted[1], ":",
    ridx_sorted[n], ":",
    sum(ridx_sorted), ":",
    sum((seq_len(n) %% 97L) * ridx_sorted)
  )
  pos_sig <- if (!length(pos_sorted)) "0" else paste0(
    length(pos_sorted), ":",
    as_int(min(pos_sorted, na.rm = TRUE)), ":",
    as_int(max(pos_sorted, na.rm = TRUE))
  )
  paste(as.character(ld_id)[1], chr, as_int(lead_pos), ridx_sig, pos_sig, sep = "|")
}

.gcanvas_ld_mem_cache_get <- function(key) {
  if (!.gcanvas_ld_mem_cache_enabled()) return(NULL)
  key <- as.character(key)[1]
  e <- .gcanvas_ld_mem_cache
  if (!exists(key, envir = e$entries, inherits = FALSE)) return(NULL)
  x <- get(key, envir = e$entries, inherits = FALSE)
  e$order <- c(e$order[e$order != key], key)
  x$value
}

.gcanvas_ld_mem_cache_prune <- function() {
  e <- .gcanvas_ld_mem_cache
  max_entries <- .gcanvas_ld_mem_cache_max_entries()
  max_bytes <- .gcanvas_ld_mem_cache_max_bytes()
  while (length(e$order) > max_entries || e$total_bytes > max_bytes) {
    if (!length(e$order)) break
    victim <- e$order[1]
    e$order <- e$order[-1]
    if (exists(victim, envir = e$entries, inherits = FALSE)) {
      x <- get(victim, envir = e$entries, inherits = FALSE)
      e$total_bytes <- e$total_bytes - as.numeric(x$bytes)
      rm(list = victim, envir = e$entries)
    }
  }
  if (!is.finite(e$total_bytes) || e$total_bytes < 0) e$total_bytes <- 0
  invisible(TRUE)
}

.gcanvas_ld_mem_cache_set <- function(key, value) {
  if (!.gcanvas_ld_mem_cache_enabled()) return(invisible(FALSE))
  key <- as.character(key)[1]
  e <- .gcanvas_ld_mem_cache
  bytes <- as.numeric(utils::object.size(value))
  if (!is.finite(bytes) || is.na(bytes) || bytes <= 0) return(invisible(FALSE))
  if (bytes > .gcanvas_ld_mem_cache_max_bytes()) return(invisible(FALSE))

  if (exists(key, envir = e$entries, inherits = FALSE)) {
    old <- get(key, envir = e$entries, inherits = FALSE)
    e$total_bytes <- e$total_bytes - as.numeric(old$bytes)
    rm(list = key, envir = e$entries)
  }

  assign(key, list(value = value, bytes = bytes), envir = e$entries)
  e$total_bytes <- e$total_bytes + bytes
  e$order <- c(e$order[e$order != key], key)
  .gcanvas_ld_mem_cache_prune()
  invisible(TRUE)
}

.gcanvas_cache_counts <- function() {
  mem_n <- 0L
  e <- tryCatch(.gcanvas_ld_mem_cache, error = function(err) NULL)
  if (is.environment(e) && is.environment(e$entries)) {
    mem_n <- as_int(length(ls(envir = e$entries, all.names = TRUE)))
  }

  reg <- .gcanvas_session_cache_registry
  created_files <- unique(as.character(reg$created_files))
  created_files <- created_files[!is.na(created_files) & nzchar(created_files)]
  files_exist <- created_files[file.exists(created_files)]

  created_dirs <- unique(as.character(reg$created_dirs))
  created_dirs <- created_dirs[!is.na(created_dirs) & nzchar(created_dirs)]
  if (length(created_files)) {
    created_dirs <- unique(c(created_dirs, dirname(created_files)))
  }
  created_dirs <- created_dirs[dir.exists(created_dirs)]

  lock_files <- character(0)
  if (length(created_dirs)) {
    for (d in created_dirs) {
      lf <- list.files(d, pattern = "\\.ldref\\.lock$", full.names = TRUE)
      if (length(lf)) lock_files <- c(lock_files, lf)
    }
  }
  lock_files <- unique(lock_files[file.exists(lock_files)])

  list(
    memory = list(ld_entries_removed = as_int(mem_n)),
    disk = list(
      files_removed = as_int(length(files_exist)),
      lock_files_removed = as_int(length(lock_files)),
      dirs_removed = as_int(length(created_dirs))
    )
  )
}

ref_bim_table <- function(ld_obj) {
  require_pkg("data.table")
  map <- ld_obj$map
  chr <- normalize.chrom(map$chromosome)
  pos <- as_int(map$physical.pos)
  a1_bim <- toupper(as.character(map$allele1))
  a2_bim <- toupper(as.character(map$allele2))
  ra1 <- canon_a1(a1_bim, a2_bim)
  ra2 <- canon_a2(a1_bim, a2_bim)

  dt <- data.table::data.table(chr = chr, pos = pos, ra1 = ra1, ra2 = ra2, a1_bim = a1_bim, a2_bim = a2_bim)
  dt[, ridx := .I]
  data.table::setkey(dt, chr, pos, ra1, ra2)
  dt <- dt[!duplicated(dt, by = c("chr","pos","ra1","ra2"))]
  dt
}

ld_r2_bigsnpr <- function(ld_obj, ref_dt, chr, start, end,
                          df_pos, lead_pos,
                          df_ea = NULL, df_nea = NULL,
                          threads = 1L) {
  require_pkg(c("bigsnpr", "data.table"))
  threads <- .gcanvas_resolve_threads(threads)

  chr_val <- normalize.chrom(chr)[1]
  start_val <- max(1L, as_int(start))
  end_val <- as_int(end)
  if (is.na(end_val) || end_val < start_val) stop("Invalid start/end.", call. = FALSE)

  df_pos0 <- .gcanvas_as_num2(df_pos)
  dt <- data.table::data.table(chr = rep(chr_val, length(df_pos0)), pos = as_int(df_pos0))
  dt <- dt[is.finite(pos) & !is.na(pos)]
  if (!nrow(dt)) {
    return(list(r2 = rep(NA_real_, length(df_pos0)),
                meta = list(n_mapped = 0L, n_missing = length(df_pos0), lead_in_ref = FALSE)))
  }

  have_allele <- !is.null(df_ea) && !is.null(df_nea)
  if (have_allele) {
    ea0 <- toupper(as.character(df_ea))
    nea0 <- toupper(as.character(df_nea))
    dt[, ea := ea0[seq_len(.N)]]
    dt[, nea := nea0[seq_len(.N)]]
    dt[, ra1 := canon_a1(ea, nea)]
    dt[, ra2 := canon_a2(ea, nea)]
    data.table::setkey(dt, chr, pos, ra1, ra2)
    hit <- ref_dt[dt, nomatch = 0L]
    miss <- dt[!dt$pos %in% hit$pos]
    if (nrow(miss)) {
      hit2 <- ref_dt[miss[, .(chr, pos)], on = .(chr, pos), nomatch = 0L, mult = "first"]
      hit <- data.table::rbindlist(list(hit, hit2), fill = TRUE)
    }
  } else {
    hit <- ref_dt[dt, on = .(chr, pos), nomatch = 0L, mult = "first"]
  }

  hit <- unique(hit[, .(chr, pos, ridx)])
  hit <- hit[pos >= start_val & pos <= end_val]
  if (nrow(hit) < 3) {
    r2_out <- rep(NA_real_, length(df_pos0))
    return(list(r2 = r2_out, meta = list(n_mapped = nrow(hit), n_missing = sum(is.na(r2_out)), lead_in_ref = FALSE)))
  }

  lead_pos <- as_int(lead_pos)
  if (!(lead_pos %in% hit$pos)) {
    lead_row <- ref_dt[data.table::data.table(chr = chr_val, pos = lead_pos),
                       on = .(chr, pos), nomatch = 0L, mult = "first"]
    if (nrow(lead_row) == 1) hit <- data.table::rbindlist(list(hit, lead_row[, .(chr, pos, ridx)]), fill = TRUE)
  }
  hit <- unique(hit, by = c("chr","pos","ridx"))

  pos_all <- ld_obj$map$physical.pos
  hit[, pos_ref := pos_all[as_int(ridx)]]
  hit <- hit[is.finite(pos_ref) & !is.na(pos_ref)]
  data.table::setorder(hit, pos_ref, ridx)

  lead_i <- match(lead_pos, hit$pos)
  if (is.na(lead_i)) {
    r2_out <- rep(NA_real_, length(df_pos0))
    return(list(r2 = r2_out, meta = list(n_mapped = nrow(hit), n_missing = sum(is.na(r2_out)), lead_in_ref = FALSE)))
  }

  ridx_sorted <- as_int(hit$ridx)
  pos_sorted <- hit$pos_ref
  size_pos <- max(pos_sorted, na.rm = TRUE) - min(pos_sorted, na.rm = TRUE) + 1
  if (!is.finite(size_pos) || size_pos <= 0) size_pos <- 1e12

  lead_tag <- if (is.finite(lead_pos)) paste0(chr_val, ":", as_int(lead_pos)) else paste0(chr_val, ":NA")
  cache_key <- .gcanvas_ld_mem_cache_key(
    ld_obj = ld_obj, chr = chr_val, lead_pos = lead_pos,
    ridx_sorted = ridx_sorted, pos_sorted = pos_sorted
  )
  map_r2 <- .gcanvas_ld_mem_cache_get(cache_key)
  cache_hit <- !is.null(map_r2)

  if (cache_hit) {
    .gcanvas_note("gcanvas", sprintf("LD memory-cache hit for %s: n_variants=%d", lead_tag, length(ridx_sorted)))
  } else {
    .gcanvas_note("gcanvas", sprintf("LD compute for %s: n_variants=%d", lead_tag, length(ridx_sorted)))
    st <- .gcanvas_ld_thread_lock(1L)
    on.exit(.gcanvas_ld_thread_restore(st), add = TRUE)

    R_sp <- bigsnpr::snp_cor(ld_obj$genotypes, ind.col = ridx_sorted, infos.pos = pos_sorted,
                             size = size_pos, ncores = as_int(threads))
    r <- as.numeric(R_sp[lead_i, ])
    r2 <- r * r
    r2[!is.finite(r2)] <- NA_real_
    r2[lead_i] <- 1

    map_r2 <- data.table::data.table(pos = hit$pos, r2 = r2)
    .gcanvas_ld_mem_cache_set(cache_key, map_r2)
  }

  data.table::setkey(map_r2, pos)

  q_dt <- data.table::data.table(pos = as_int(df_pos0))
  r2_out <- map_r2[q_dt, on = "pos"]$r2

  list(r2 = r2_out, meta = list(
    n_mapped = as_int(nrow(hit)),
    n_missing = as_int(sum(is.na(r2_out))),
    lead_in_ref = TRUE,
    ld_mem_cache_hit = cache_hit
  ))
}

.gcanvas_ld_extract_named_matrix <- function(x, max_depth = 6L) {
  max_depth <- as_int(max_depth)
  if (is.na(max_depth) || max_depth < 0L) max_depth <- 0L

  .gcanvas_is_named_square_numeric <- function(m) {
    if (!is.matrix(m) || !is.numeric(m)) return(FALSE)
    if (nrow(m) != ncol(m) || nrow(m) < 1L) return(FALSE)
    rn <- rownames(m); cn <- colnames(m)
    if (is.null(rn) || is.null(cn)) return(FALSE)
    if (!length(rn) || !length(cn)) return(FALSE)
    TRUE
  }

  if (.gcanvas_is_named_square_numeric(x)) return(x)
  if (max_depth == 0L || !is.list(x)) return(NULL)

  nms <- names(x)
  if (!is.null(nms) && length(nms)) {
    pri <- c("ld", "LD", "r", "R", "corr", "cor", "matrix", "mat")
    idx_pri <- which(nms %in% pri)
    if (length(idx_pri)) {
      for (i in idx_pri) {
        out <- .gcanvas_ld_extract_named_matrix(x[[i]], max_depth = max_depth - 1L)
        if (!is.null(out)) return(out)
      }
    }
  }

  for (i in seq_along(x)) {
    out <- .gcanvas_ld_extract_named_matrix(x[[i]], max_depth = max_depth - 1L)
    if (!is.null(out)) return(out)
  }
  NULL
}

ld_r2_from_rds <- function(ld_rds, snp_ids, lead_snp) {
  ld_rds <- abs_path(as.character(ld_rds)[1])
  if (!file.exists(ld_rds)) stop("ld.rds file not found: ", ld_rds, call. = FALSE)

  obj <- readRDS(ld_rds)
  mat <- .gcanvas_ld_extract_named_matrix(obj, max_depth = 8L)
  if (is.null(mat)) {
    stop("No named square numeric LD matrix found in ld.rds.", call. = FALSE)
  }

  rn <- rownames(mat)
  cn <- colnames(mat)
  if (!identical(rn, cn)) {
    shared <- intersect(rn, cn)
    if (!length(shared)) stop("ld.rds matrix row/col names do not overlap.", call. = FALSE)
    mat <- mat[shared, shared, drop = FALSE]
    rn <- rownames(mat)
    cn <- colnames(mat)
  }
  if (!(lead_snp %in% rn)) {
    return(list(
      r2 = rep(NA_real_, length(snp_ids)),
      meta = list(n_mapped = 0L, n_missing = length(snp_ids), lead_in_ref = FALSE, source = ld_rds)
    ))
  }

  lead_tag <- as.character(lead_snp)[1]
  .gcanvas_note("gcanvas::ld_r2_from_rds", sprintf("LD compute for %s: n_variants=%d", lead_tag, nrow(mat)))

  q <- as.character(snp_ids)
  hit <- match(q, rn)
  r <- rep(NA_real_, length(q))
  ok <- !is.na(hit)
  if (any(ok)) r[ok] <- as.numeric(mat[lead_snp, rn[hit[ok]], drop = TRUE])
  r2 <- r * r
  r2[!is.finite(r2)] <- NA_real_
  r2 <- pmax(0, pmin(1, r2))
  r2[q == lead_snp] <- 1

  list(
    r2 = r2,
    meta = list(
      n_mapped = as_int(sum(ok)),
      n_missing = as_int(sum(is.na(r2))),
      lead_in_ref = TRUE,
      source = ld_rds
    )
  )
}

ld_r2_from_matrix <- function(ld_matrix, snp_ids, lead_snp, source_label = "ld.matrix") {
  mat <- .gcanvas_ld_extract_named_matrix(ld_matrix, max_depth = 8L)
  if (is.null(mat)) {
    stop("No named square numeric LD matrix found in ld.matrix.", call. = FALSE)
  }

  rn <- rownames(mat)
  cn <- colnames(mat)
  if (!identical(rn, cn)) {
    shared <- intersect(rn, cn)
    if (!length(shared)) stop("ld.matrix row/col names do not overlap.", call. = FALSE)
    mat <- mat[shared, shared, drop = FALSE]
    rn <- rownames(mat)
  }
  if (!(lead_snp %in% rn)) {
    return(list(
      r2 = rep(NA_real_, length(snp_ids)),
      meta = list(n_mapped = 0L, n_missing = length(snp_ids), lead_in_ref = FALSE, source = source_label)
    ))
  }

  lead_tag <- as.character(lead_snp)[1]
  .gcanvas_note("gcanvas::ld_r2_from_matrix", sprintf("LD compute for %s: n_variants=%d", lead_tag, nrow(mat)))

  q <- as.character(snp_ids)
  hit <- match(q, rn)
  r <- rep(NA_real_, length(q))
  ok <- !is.na(hit)
  if (any(ok)) r[ok] <- as.numeric(mat[lead_snp, rn[hit[ok]], drop = TRUE])
  r2 <- r * r
  r2[!is.finite(r2)] <- NA_real_
  r2 <- pmax(0, pmin(1, r2))
  r2[q == lead_snp] <- 1

  list(
    r2 = r2,
    meta = list(
      n_mapped = as_int(sum(ok)),
      n_missing = as_int(sum(is.na(r2))),
      lead_in_ref = TRUE,
      source = source_label
    )
  )
}

#' LD proxy variants for a query SNP
#'
#' Returns variants in linkage disequilibrium with the query SNP above the
#' requested r-squared threshold within `window` bp, using the configured
#' bigsnpr LD reference (or a precomputed RDS / matrix).
#'
#' @param x Query SNP id (character of length 1) or a table of queries.
#' @param rsq Minimum r-squared threshold (default 0.8).
#' @param window LD search window. Numeric (bp) or a string like `"500kb"`.
#' @param ld.rds,ld.bfile Reference panel pointers (RDS or PLINK bfile).
#' @param snp.col,chrom.col,pos.col Column names in the query table.
#' @param ld.cache.dir Optional directory for the LD reference cache.
#' @param threads Integer thread count.
#'
#' @return A `data.table` of proxies with id, chromosome, position, alleles,
#'   r/r-squared, and distance.
#' @export
ldproxy <- function(x,
                    rsq = 0.8,
                    window = "500kb",
                    ld.rds = NULL,
                    ld.bfile = NULL,
                    snp.col = "SNP",
                    chrom.col = "CHR",
                    pos.col = "POS",
                    ld.cache.dir = NULL,
                    threads = 4L) {
  require_pkg(c("data.table"))

  rsq_cut <- suppressWarnings(as.numeric(rsq))[1]
  if (!is.finite(rsq_cut) || is.na(rsq_cut) || rsq_cut < 0 || rsq_cut > 1) {
    stop("rsq must be a numeric value in [0, 1].", call. = FALSE)
  }
  window_bp <- .gcanvas_parse_bp_span(window, arg_name = "window")
  if (!is.finite(window_bp) || is.na(window_bp) || window_bp <= 0) {
    stop("window must be a positive bp span.", call. = FALSE)
  }
  ld_threads <- .gcanvas_resolve_threads(threads)

  empty_out <- data.table::data.table(
    SNP = character(), CHR = character(), POS = numeric(), proxyID = character(), r2 = numeric()
  )

  marker_dt <- NULL
  if (is.data.frame(x) || data.table::is.data.table(x)) {
    dt0 <- if (data.table::is.data.table(x)) data.table::copy(x) else data.table::as.data.table(x)
    snp_col_use <- .gcanvas_resolve_colname(names(dt0), snp.col, aliases = character(), required = TRUE, arg_label = "snp.col")
    chrom_col_use <- .gcanvas_resolve_colname(names(dt0), chrom.col, aliases = character(), required = TRUE, arg_label = "chrom.col")
    pos_col_use <- .gcanvas_resolve_colname(names(dt0), pos.col, aliases = character(), required = TRUE, arg_label = "pos.col")
    marker_dt <- dt0[, .(
      markerID = as.character(get(snp_col_use)),
      chrom = normalize.chrom(get(chrom_col_use)),
      pos = suppressWarnings(as.numeric(get(pos_col_use)))
    )]
  } else {
    ids <- .gcanvas_as_snp_vector(x)
    if (!length(ids)) return(empty_out)
    marker_dt <- data.table::data.table(markerID = as.character(ids))
    sp <- strsplit(marker_dt$markerID, ":", fixed = TRUE)
    marker_dt[, chrom := normalize.chrom(vapply(sp, function(v) if (length(v) >= 1L) v[1] else NA_character_, character(1)))]
    marker_dt[, pos := suppressWarnings(as.numeric(vapply(sp, function(v) if (length(v) >= 2L) v[2] else NA_character_, character(1))))]
  }

  marker_dt <- marker_dt[!is.na(markerID) & nzchar(markerID)]
  if (!nrow(marker_dt)) return(empty_out)
  marker_dt[, chrom := normalize.chrom(chrom)]
  marker_dt[, pos := suppressWarnings(as.numeric(pos))]
  marker_dt <- unique(marker_dt, by = "markerID")

  ld_rds0 <- NULL
  if (!is.null(ld.rds) && length(ld.rds) > 0L && !is.na(ld.rds[1]) && nzchar(as.character(ld.rds[1]))) {
    ld_rds0 <- abs_path(as.character(ld.rds)[1])
    if (!file.exists(ld_rds0)) {
      .gcanvas_warn_msg(paste0("ld.rds not found: ", ld_rds0, " -> trying ld.bfile"))
      ld_rds0 <- NULL
    }
  }
  ld_bfile0 <- NULL
  if (!is.null(ld.bfile) && length(ld.bfile) > 0L && !is.na(ld.bfile[1]) && nzchar(as.character(ld.bfile[1]))) {
    ld_bfile0 <- .gcanvas_normalize_bfile_prefix(ld.bfile)
    if (!is.null(ld_bfile0) && !.gcanvas_ld_bfile_ready(ld_bfile0)) {
      .gcanvas_warn_msg(paste0("ld.bfile incomplete (.bed/.bim/.fam): ", ld_bfile0))
      ld_bfile0 <- NULL
    }
  }
  if (is.null(ld_rds0) && is.null(ld_bfile0)) {
    stop("Provide valid ld.rds or ld.bfile.", call. = FALSE)
  }

  # ---- ld.rds path ----
  if (!is.null(ld_rds0)) {
    obj <- readRDS(ld_rds0)
    mat <- .gcanvas_ld_extract_named_matrix(obj, max_depth = 8L)
    if (is.null(mat)) stop("No named square numeric LD matrix found in ld.rds.", call. = FALSE)
    rn <- rownames(mat); cn <- colnames(mat)
    if (!identical(rn, cn)) {
      shared <- intersect(rn, cn)
      if (!length(shared)) stop("ld.rds matrix row/col names do not overlap.", call. = FALSE)
      mat <- mat[shared, shared, drop = FALSE]
      rn <- rownames(mat)
    }
    var_dt <- data.table::data.table(proxy = rn)
    sp <- strsplit(var_dt$proxy, ":", fixed = TRUE)
    var_dt[, chrom := normalize.chrom(vapply(sp, function(v) if (length(v) >= 1L) v[1] else NA_character_, character(1)))]
    var_dt[, pos := suppressWarnings(as.numeric(vapply(sp, function(v) if (length(v) >= 2L) v[2] else NA_character_, character(1))))]

    out_list <- vector("list", nrow(marker_dt))
    for (i in seq_len(nrow(marker_dt))) {
      mk <- marker_dt$markerID[i]
      chr_i <- marker_dt$chrom[i]
      pos_i <- marker_dt$pos[i]

      lead_id <- mk
      lead_idx <- match(lead_id, rn)
      if (is.na(lead_idx) && !is.na(chr_i) && nzchar(chr_i) && is.finite(pos_i)) {
        hit <- which(var_dt$chrom == chr_i & var_dt$pos == pos_i)
        if (length(hit)) {
          lead_id <- var_dt$proxy[hit[1]]
          lead_idx <- match(lead_id, rn)
        }
      }
      if (is.na(lead_idx)) next
      lead_chr <- var_dt$chrom[lead_idx]
      lead_pos <- var_dt$pos[lead_idx]

      r <- as.numeric(mat[lead_idx, , drop = TRUE])
      r2 <- r * r
      r2[!is.finite(r2)] <- NA_real_
      r2 <- pmax(0, pmin(1, r2))

      dt_out <- data.table::data.table(
        markerID = mk,
        chrom = var_dt$chrom,
        pos = var_dt$pos,
        proxyID = var_dt$proxy,
        r2 = r2
      )
      dt_out <- dt_out[proxyID != lead_id & is.finite(r2) & r2 > rsq_cut]
      if (!is.na(chr_i) && nzchar(chr_i) && is.finite(pos_i)) {
        dt_out <- dt_out[chrom == chr_i & is.finite(pos) & abs(pos - pos_i) <= window_bp]
      } else if (!is.na(lead_chr) && nzchar(lead_chr) && is.finite(lead_pos)) {
        dt_out <- dt_out[chrom == lead_chr & is.finite(pos) & abs(pos - lead_pos) <= window_bp]
      }
      if (nrow(dt_out)) out_list[[i]] <- dt_out
    }
    out <- data.table::rbindlist(out_list, use.names = TRUE, fill = TRUE)
    if (!nrow(out)) return(empty_out)
    out <- unique(out, by = c("markerID", "proxyID"))
    data.table::setnames(out, c("markerID", "chrom", "pos"), c("SNP", "CHR", "POS"))
    data.table::setcolorder(out, c("SNP", "CHR", "POS", "proxyID", "r2"))
    return(out[])
  }

  # ---- ld.bfile path ----
  cache_dir <- {
    x0 <- as.character(ld.cache.dir)[1]
    if (is.null(ld.cache.dir) || length(ld.cache.dir) == 0L || is.na(x0) || !nzchar(x0) || tolower(x0) == "auto") {
      .gcanvas_default_ld_cache_dir(ld_bfile0)
    } else {
      abs_path(ld.cache.dir)
    }
  }
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  ld_obj <- attach_ld_ref(bfile = ld_bfile0, ld_rds = NULL, cache_dir = cache_dir)
  map <- ld_obj$map
  var_dt <- data.table::data.table(
    ridx = seq_len(nrow(map)),
    chrom = normalize.chrom(map$chromosome),
    pos = suppressWarnings(as.numeric(map$physical.pos)),
    a1 = toupper(as.character(map$allele1)),
    a2 = toupper(as.character(map$allele2))
  )
  if ("marker.ID" %in% names(map)) {
    mid <- as.character(map$marker.ID)
    mid[is.na(mid) | !nzchar(mid)] <- NA_character_
    var_dt[, proxy := mid]
  } else {
    var_dt[, proxy := NA_character_]
  }
  var_dt[is.na(proxy) | !nzchar(proxy), proxy := paste0(chrom, ":", as_int(round(pos)), ":", a1, ":", a2)]
  var_dt <- var_dt[!is.na(chrom) & nzchar(chrom) & is.finite(pos) & !is.na(proxy) & nzchar(proxy)]
  if (!nrow(var_dt)) return(empty_out)
  data.table::setkey(var_dt, chrom, pos)

  th <- .gcanvas_ld_thread_lock(ld_threads)
  on.exit(.gcanvas_ld_thread_restore(th), add = TRUE)

  out_list <- vector("list", nrow(marker_dt))
  for (i in seq_len(nrow(marker_dt))) {
    mk <- marker_dt$markerID[i]
    chr_i <- marker_dt$chrom[i]
    pos_i <- marker_dt$pos[i]

    lead_row <- var_dt[proxy == mk][1]
    if (!nrow(lead_row) && !is.na(chr_i) && nzchar(chr_i) && is.finite(pos_i)) {
      lead_row <- var_dt[.(chr_i, pos_i)][1]
    }
    if (!nrow(lead_row)) next
    chr0 <- lead_row$chrom[1]
    pos0 <- as.numeric(lead_row$pos[1])
    ridx0 <- as_int(lead_row$ridx[1])

    sub <- var_dt[chrom == chr0 & abs(pos - pos0) <= window_bp]
    if (nrow(sub) < 2L) next
    data.table::setorder(sub, pos, ridx)
    ridx_sorted <- as_int(sub$ridx)
    lead_local <- match(ridx0, ridx_sorted)
    if (is.na(lead_local)) next
    pos_sorted <- suppressWarnings(as.numeric(ld_obj$map$physical.pos[ridx_sorted]))
    size_pos <- max(pos_sorted, na.rm = TRUE) - min(pos_sorted, na.rm = TRUE) + 1
    if (!is.finite(size_pos) || size_pos <= 0) size_pos <- window_bp * 2 + 1

    R_sp <- bigsnpr::snp_cor(
      ld_obj$genotypes,
      ind.col = ridx_sorted,
      infos.pos = pos_sorted,
      size = size_pos,
      ncores = ld_threads,
      alpha = 1,
      thr_r2 = 0
    )
    r <- as.numeric(R_sp[lead_local, ])
    r2 <- r * r
    r2[!is.finite(r2)] <- NA_real_
    r2 <- pmax(0, pmin(1, r2))

    sub_out <- data.table::data.table(
      markerID = mk,
      chrom = sub$chrom,
      pos = as.numeric(sub$pos),
      proxyID = sub$proxy,
      r2 = r2
    )
    sub_out <- sub_out[!(chrom == chr0 & pos == pos0 & proxyID == lead_row$proxy[1]) & is.finite(r2) & r2 > rsq_cut]
    if (nrow(sub_out)) out_list[[i]] <- sub_out
  }
  out <- data.table::rbindlist(out_list, use.names = TRUE, fill = TRUE)
  if (!nrow(out)) return(empty_out)
  out <- unique(out, by = c("markerID", "proxyID"))
  data.table::setnames(out, c("markerID", "chrom", "pos"), c("SNP", "CHR", "POS"))
  data.table::setcolorder(out, c("SNP", "CHR", "POS", "proxyID", "r2"))
  out[]
}

#' Compute LD (r and r-squared) between variants
#'
#' Flexible LD calculator covering single-variant lookup, pairwise mode, set
#' mode (full square matrix), and lead-vs-region mode. Resolves variants via
#' SNP id or coordinates, attaches the configured bigsnpr LD reference, and
#' optionally writes the long-form pair table to disk.
#'
#' @param snp Character vector of SNP ids (or a table).
#' @param chrom,pos Coordinate-based input alternative to `snp`.
#' @param data Optional `data.frame` providing both lookups and metadata.
#' @param a1 Optional A1 mapping for sign-aware alignment.
#' @param snp.col,chrom.col,pos.col,a1.col Column names in `data`.
#' @param flank Window size (numeric bp or e.g. `"250kb"`) for flank mode.
#' @param r2 Logical. Return r-squared instead of signed r.
#' @param matrix Logical. Return a square matrix in addition to / instead of
#'   the pair table.
#' @param pairwise One of `"auto"`, `TRUE`, `FALSE` — pairwise vs. lead-vs-set.
#' @param mode One of `"auto"`, `"flank"`, `"set"` — region resolution.
#' @param ld.bfile,ld.rds,ld.cache.dir Reference panel pointers.
#' @param threads Integer thread count.
#' @param silent Logical. Suppress progress notes.
#'
#' @return A list containing `data` (pair table) and/or `R` (matrix), the
#'   variant `a1` mapping, and a `meta` slot of run metadata.
#' @export
calcld <- function(snp = NULL,
                   chrom = NULL,
                   pos = NULL,
                   data = NULL,
                   a1 = NULL,
                   snp.col = "SNP",
                   chrom.col = "CHR",
                   pos.col = "POS",
                   a1.col = "A1",
                   flank = "250kb",
                   r2 = TRUE,
                   matrix = FALSE,
                   pairwise = c("auto", TRUE, FALSE),
                   mode = c("auto", "flank", "set"),
                   ld.bfile = NULL,
                   ld.rds = NULL,
                   ld.cache.dir = NULL,
                   threads = 4L,
                   silent = FALSE) {
  require_pkg(c("data.table"))
  use_r2 <- isTRUE(r2)
  use_matrix <- isTRUE(matrix)
  silent <- isTRUE(silent)

  flank_bp <- .gcanvas_parse_bp_span(flank, arg_name = "flank")
  if (!is.finite(flank_bp) || is.na(flank_bp) || flank_bp <= 0) {
    stop("flank must be numeric (bp) or string like '250kb'/'1mb'.", call. = FALSE)
  }

  .calcld_parse_mode <- function(x) {
    m <- tolower(trimws(as.character(x)[1]))
    if (is.na(m) || !nzchar(m)) m <- "auto"
    if (!(m %in% c("auto", "flank", "set"))) {
      stop("mode must be one of c('auto','flank','set').", call. = FALSE)
    }
    m
  }
  .calcld_parse_pairwise <- function(x, mode_use) {
    val <- NA
    if (is.logical(x) && length(x) == 1L && !is.na(x)) {
      val <- isTRUE(x)
    } else {
      x0 <- tolower(trimws(as.character(x)[1]))
      if (is.na(x0) || !nzchar(x0) || x0 == "auto") {
        val <- NA
      } else if (x0 %in% c("true", "t", "1", "yes", "y")) {
        val <- TRUE
      } else if (x0 %in% c("false", "f", "0", "no", "n")) {
        val <- FALSE
      } else {
        stop("pairwise must be one of c('auto', TRUE, FALSE).", call. = FALSE)
      }
    }
    if (is.na(val)) {
      if (identical(mode_use, "flank")) FALSE else TRUE
    } else {
      isTRUE(val)
    }
  }
  .calcld_resolve_named_or_positional <- function(v, n, keys = NULL) {
    if (is.null(v) || length(v) == 0L) return(rep(NA_character_, n))
    if (is.list(v) && !is.data.frame(v) && !data.table::is.data.table(v)) {
      v <- unlist(v, use.names = TRUE)
    }
    vv <- toupper(as.character(v))
    vv[is.na(vv) | !nzchar(vv)] <- NA_character_
    nms <- names(vv)
    if (!is.null(keys) && !is.null(nms) && any(!is.na(nms) & nzchar(nms))) {
      out <- unname(vv[as.character(keys)])
      if (length(out) == n) return(out)
    }
    if (length(vv) == 1L) return(rep(vv[1], n))
    if (length(vv) == n) return(vv)
    rep(NA_character_, n)
  }
  .calcld_empty_pair_dt <- function(value_col = "r2") {
    out <- data.table::data.table(
      markerID = character(),
      chrom = character(),
      pos = numeric(),
      markerID2 = character()
    )
    out[, (value_col) := numeric()]
    data.table::setcolorder(out, c("markerID", "chrom", "pos", "markerID2", value_col))
    out
  }
  .calcld_pair_template <- function(vdt, pairwise_use, value_col = "r2") {
    vdt <- data.table::as.data.table(vdt)
    if (!nrow(vdt) || nrow(vdt) < 2L) return(.calcld_empty_pair_dt(value_col))
    if (isTRUE(pairwise_use)) {
      cb <- utils::combn(seq_len(nrow(vdt)), 2L)
      i <- cb[1, ]; j <- cb[2, ]
    } else {
      i <- rep(1L, nrow(vdt) - 1L)
      j <- seq.int(2L, nrow(vdt))
    }
    out <- data.table::data.table(
      markerID = vdt$markerID[i],
      chrom = vdt$chrom[i],
      pos = as.numeric(vdt$pos[i]),
      markerID2 = vdt$markerID[j]
    )
    out[, (value_col) := NA_real_]
    data.table::setcolorder(out, c("markerID", "chrom", "pos", "markerID2", value_col))
    out
  }
  .calcld_a1_named <- function(vdt = NULL, out_dt = NULL) {
    if (isTRUE(use_r2)) return(NULL)
    if (is.null(vdt) || !nrow(vdt) || !("markerID" %in% names(vdt)) || !("a1_out" %in% names(vdt))) {
      return(stats::setNames(character(), character()))
    }
    vv <- data.table::as.data.table(vdt)[, .(markerID = as.character(markerID), a1_out = as.character(a1_out))]
    vv <- vv[!is.na(markerID) & nzchar(markerID)]
    vv <- vv[!duplicated(markerID)]
    if (!nrow(vv)) return(stats::setNames(character(), character()))
    out <- stats::setNames(vv$a1_out, vv$markerID)
    if (!is.null(out_dt) && nrow(out_dt)) {
      id1 <- if ("SNP" %in% names(out_dt)) "SNP" else "markerID"
      id2 <- if ("SNP2" %in% names(out_dt)) "SNP2" else "markerID2"
      ids <- unique(c(as.character(out_dt[[id1]]), as.character(out_dt[[id2]])))
      ids <- ids[!is.na(ids) & nzchar(ids)]
      out <- out[ids]
    }
    out
  }
  .calcld_meta_pack <- function(chrom = character(),
                         start = NA_real_,
                         end = NA_real_,
                         n_locus = 0L,
                         n_mapped_use = 0L) {
    list(
      chrom = chrom,
      start = start,
      end = end,
      n_locus = as_int(n_locus),
      n_input = as_int(n_input),
      n_mapped = as_int(n_mapped_use),
      mode = mode_use,
      pairwise = pairwise_use,
      source = ld_mode,
      ld_path = as.character(ld_path)[1]
    )
  }
  .calcld_ret_table <- function(dt, a1_vec = NULL, meta = NULL) {
    dt <- if (data.table::is.data.table(dt)) data.table::copy(dt) else data.table::as.data.table(dt)
    if ("a1" %in% names(dt)) dt[, a1 := NULL]
    if ("markerID" %in% names(dt)) data.table::setnames(dt, "markerID", "SNP")
    if ("markerID2" %in% names(dt)) data.table::setnames(dt, "markerID2", "SNP2")
    if ("chrom" %in% names(dt)) data.table::setnames(dt, "chrom", "CHR")
    if ("pos" %in% names(dt)) data.table::setnames(dt, "pos", "POS")
    val_col <- intersect(c("r2", "r"), names(dt))
    if ("SNP" %in% names(dt) && "CHR" %in% names(dt) && "POS" %in% names(dt) && "SNP2" %in% names(dt) && length(val_col) == 1L) {
      data.table::setcolorder(dt, c("SNP", "CHR", "POS", "SNP2", val_col))
    }
    if (is.null(a1_vec)) a1_vec <- .calcld_a1_named(NULL, NULL)
    if (is.null(meta)) meta <- .calcld_meta_pack()
    list(data = dt[], a1 = a1_vec, meta = meta)
  }

  # -------- Parse input variants --------
  seed_dt <- data.table::data.table(
    markerID = character(),
    chrom = character(),
    pos = numeric(),
    a1_input = character(),
    input_order = integer()
  )

  if (!is.null(data)) {
    dt0 <- if (data.table::is.data.table(data)) data.table::copy(data) else data.table::as.data.table(data)
    pick_col <- function(primary, secondary, required = FALSE, aliases = character()) {
      out <- .gcanvas_resolve_colname(
        names(dt0),
        primary = if (!is.null(primary) && length(primary) == 1L && is.character(primary) && !is.na(primary)) primary else NA_character_,
        aliases = c(as.character(secondary)[1], aliases),
        required = FALSE
      )
      if (isTRUE(required) && is.null(out)) {
        stop(sprintf("Required column not found in data (primary='%s', secondary='%s').",
                     as.character(primary)[1], as.character(secondary)[1]), call. = FALSE)
      }
      out
    }
    chr_use <- pick_col(chrom, chrom.col, required = TRUE, aliases = character())
    pos_use <- pick_col(pos, pos.col, required = TRUE, aliases = character())
    snp_use <- pick_col(snp, snp.col, required = FALSE, aliases = character())
    a1_use <- pick_col(a1, a1.col, required = FALSE, aliases = c("A1", "a1"))

    seed_dt <- data.table::data.table(
      markerID = if (!is.null(snp_use)) as.character(dt0[[snp_use]]) else NA_character_,
      chrom = normalize.chrom(dt0[[chr_use]]),
      pos = suppressWarnings(as.numeric(dt0[[pos_use]])),
      a1_input = if (!is.null(a1_use)) toupper(as.character(dt0[[a1_use]])) else NA_character_,
      input_order = seq_len(nrow(dt0))
    )
  } else {
    has_chr_pos <- !is.null(chrom) && length(chrom) > 0L && !is.null(pos) && length(pos) > 0L
    if (isTRUE(has_chr_pos)) {
      chr_v <- normalize.chrom(chrom)
      pos_v <- suppressWarnings(as.numeric(pos))
      n <- max(length(chr_v), length(pos_v))
      if (length(chr_v) == 1L && n > 1L) chr_v <- rep(chr_v, n)
      if (length(pos_v) == 1L && n > 1L) pos_v <- rep(pos_v, n)
      if (length(chr_v) != length(pos_v)) stop("chrom and pos lengths must match (or one must be length 1).", call. = FALSE)

      snp_v <- NULL
      if (!is.null(snp) && length(snp) > 0L) {
        snp_v <- as.character(snp)
        if (length(snp_v) == 1L && n > 1L) snp_v <- rep(snp_v, n)
        if (length(snp_v) != n) snp_v <- NULL
      }
      marker_v <- if (is.null(snp_v)) paste0(chr_v, ":", as_int(pos_v)) else snp_v
      key_v <- if (is.null(snp_v)) marker_v else snp_v
      a1_v <- .calcld_resolve_named_or_positional(a1, n = n, keys = key_v)

      seed_dt <- data.table::data.table(
        markerID = as.character(marker_v),
        chrom = chr_v,
        pos = pos_v,
        a1_input = a1_v,
        input_order = seq_len(n)
      )
    } else if (!is.null(snp) && length(snp) > 0L) {
      snp_v <- .gcanvas_as_snp_vector(snp)
      if (!length(snp_v)) stop("No valid input variants.", call. = FALSE)
      sp <- strsplit(snp_v, ":", fixed = TRUE)
      chr_v <- normalize.chrom(vapply(sp, function(z) if (length(z) >= 1L) z[1] else NA_character_, character(1)))
      pos_v <- suppressWarnings(as.numeric(vapply(sp, function(z) if (length(z) >= 2L) z[2] else NA_character_, character(1))))
      a1_v <- .calcld_resolve_named_or_positional(a1, n = length(snp_v), keys = snp_v)
      seed_dt <- data.table::data.table(
        markerID = as.character(snp_v),
        chrom = chr_v,
        pos = pos_v,
        a1_input = a1_v,
        input_order = seq_along(snp_v)
      )
    } else {
      stop("Provide input via (chrom + pos), snp, or data.", call. = FALSE)
    }
  }

  seed_dt[, markerID := as.character(markerID)]
  seed_dt[, chrom := normalize.chrom(chrom)]
  seed_dt[, pos := suppressWarnings(as.numeric(pos))]
  seed_dt[, a1_input := toupper(as.character(a1_input))]
  seed_dt[is.na(markerID) | !nzchar(markerID), markerID := ifelse(!is.na(chrom) & nzchar(chrom) & is.finite(pos), paste0(chrom, ":", as_int(pos)), NA_character_)]
  seed_dt <- seed_dt[(!is.na(markerID) & nzchar(markerID)) | (!is.na(chrom) & nzchar(chrom) & is.finite(pos))]
  if (!nrow(seed_dt)) stop("No valid input variants after parsing.", call. = FALSE)
  data.table::setorder(seed_dt, input_order)
  seed_dt[, dedup_key := ifelse(!is.na(chrom) & nzchar(chrom) & is.finite(pos), paste0(chrom, ":", as_int(pos)), markerID)]
  seed_dt <- seed_dt[!duplicated(dedup_key)]
  seed_dt[, dedup_key := NULL]
  if (!nrow(seed_dt)) stop("No valid input variants after deduplication.", call. = FALSE)
  n_input <- nrow(seed_dt)

  mode_use <- .calcld_parse_mode(mode)
  if (identical(mode_use, "auto")) {
    if (!is.null(data)) {
      mode_use <- "set"
    } else if (n_input <= 1L) {
      mode_use <- "flank"
    } else {
      mode_use <- "set"
    }
  }
  if (identical(mode_use, "set") && n_input <= 1L) {
    .gcanvas_note("gcanvas::calcld", "mode='set' requires >=2 variants; returning NULL.", silent = silent)
    return(NULL)
  }
  pairwise_use <- .calcld_parse_pairwise(pairwise, mode_use = mode_use)
  if (isTRUE(use_matrix) && !isTRUE(pairwise_use)) {
    pairwise_use <- TRUE
    .gcanvas_note("gcanvas::calcld", "matrix=TRUE -> forcing pairwise=TRUE", silent = silent)
  }

  # -------- Attach LD reference --------
  ld_rds0 <- NULL
  if (!is.null(ld.rds) && length(ld.rds) > 0L && !is.na(ld.rds[1]) && nzchar(as.character(ld.rds[1]))) {
    ld_rds0 <- abs_path(as.character(ld.rds)[1])
    if (!file.exists(ld_rds0)) {
      .gcanvas_warn_msg(paste0("ld.rds not found: ", ld_rds0, " -> trying ld.bfile"))
      ld_rds0 <- NULL
    }
  }
  ld_bfile0 <- NULL
  if (!is.null(ld.bfile) && length(ld.bfile) > 0L && !is.na(ld.bfile[1]) && nzchar(as.character(ld.bfile[1]))) {
    ld_bfile0 <- .gcanvas_normalize_bfile_prefix(ld.bfile)
    if (!is.null(ld_bfile0) && !.gcanvas_ld_bfile_ready(ld_bfile0)) {
      .gcanvas_warn_msg(paste0("ld.bfile incomplete (.bed/.bim/.fam): ", ld_bfile0))
      ld_bfile0 <- NULL
    }
  }
  if (is.null(ld_rds0) && is.null(ld_bfile0)) {
    stop("Provide valid ld.rds or ld.bfile.", call. = FALSE)
  }

  ld_mode <- if (!is.null(ld_rds0)) "rds" else "bfile"
  ld_path <- if (identical(ld_mode, "rds")) abs_path(ld_rds0) else abs_path(ld_bfile0)
  ld_threads <- .gcanvas_resolve_threads(threads)

  ref_dt <- NULL
  mat <- NULL
  ld_obj <- NULL
  if (identical(ld_mode, "rds")) {
    obj <- readRDS(ld_rds0)
    mat <- .gcanvas_ld_extract_named_matrix(obj, max_depth = 8L)
    if (is.null(mat)) stop("No named square numeric LD matrix found in ld.rds.", call. = FALSE)
    rn <- rownames(mat); cn <- colnames(mat)
    if (!identical(rn, cn)) {
      shared <- intersect(rn, cn)
      if (!length(shared)) stop("ld.rds matrix row/col names do not overlap.", call. = FALSE)
      mat <- mat[shared, shared, drop = FALSE]
      rn <- rownames(mat)
    }
    sp <- strsplit(rn, ":", fixed = TRUE)
    ref_dt <- data.table::data.table(
      ref_idx = seq_along(rn),
      ref_id = as.character(rn),
      chrom = normalize.chrom(vapply(sp, function(v) if (length(v) >= 1L) v[1] else NA_character_, character(1))),
      pos = suppressWarnings(as.numeric(vapply(sp, function(v) if (length(v) >= 2L) v[2] else NA_character_, character(1)))),
      a1_ref = toupper(vapply(sp, function(v) if (length(v) >= 3L) v[3] else NA_character_, character(1))),
      a2_ref = toupper(vapply(sp, function(v) if (length(v) >= 4L) v[4] else NA_character_, character(1)))
    )
  } else {
    lim <- suppressWarnings(as.numeric(getOption("bigstatsr.ncores.max"))[1])
    if (!is.finite(lim) || is.na(lim) || lim < 1) {
      nc <- suppressWarnings(as.integer(parallel::detectCores()))
      if (!is.finite(nc) || is.na(nc) || nc < 1L) nc <- 1L
      options(bigstatsr.ncores.max = nc)
    }
    cache_dir <- {
      x0 <- as.character(ld.cache.dir)[1]
      if (is.null(ld.cache.dir) || length(ld.cache.dir) == 0L || is.na(x0) || !nzchar(x0) || tolower(x0) == "auto") {
        .gcanvas_default_ld_cache_dir(ld_bfile0)
      } else {
        abs_path(ld.cache.dir)
      }
    }
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    ld_obj <- attach_ld_ref(bfile = ld_bfile0, ld_rds = NULL, cache_dir = cache_dir)
    map <- ld_obj$map
    mid <- if ("marker.ID" %in% names(map)) as.character(map$marker.ID) else rep(NA_character_, nrow(map))
    mid[is.na(mid) | !nzchar(mid)] <- NA_character_
    ref_dt <- data.table::data.table(
      ref_idx = seq_len(nrow(map)),
      ref_id = mid,
      chrom = normalize.chrom(map$chromosome),
      pos = suppressWarnings(as.numeric(map$physical.pos)),
      a1_ref = toupper(as.character(map$allele1)),
      a2_ref = toupper(as.character(map$allele2))
    )
    ref_dt[is.na(ref_id) | !nzchar(ref_id), ref_id := paste0(chrom, ":", as_int(round(pos)), ":", a1_ref, ":", a2_ref)]
  }
  ref_dt <- ref_dt[!is.na(ref_idx)]
  ref_dt <- ref_dt[!is.na(ref_id) & nzchar(ref_id)]
  ref_dt <- ref_dt[!is.na(chrom) & nzchar(chrom) & is.finite(pos)]
  if (!nrow(ref_dt)) stop("No valid LD reference variants.", call. = FALSE)
  ref_by_id <- unique(ref_dt, by = "ref_id")
  ref_by_pos <- unique(ref_dt, by = c("chrom", "pos"))
  data.table::setkey(ref_by_pos, chrom, pos)

  # -------- Map input variants --------
  seed_dt[, input_idx := seq_len(.N)]
  q_id <- seed_dt[!is.na(markerID) & nzchar(markerID), .(input_idx, markerID)]
  hit_id <- if (nrow(q_id)) {
    ref_by_id[q_id, on = .(ref_id = markerID), nomatch = 0L, mult = "first"][
      , .(input_idx, ref_idx, ref_id, chrom_ref = chrom, pos_ref = pos, a1_ref, a2_ref)
    ]
  } else data.table::data.table(input_idx = integer(), ref_idx = integer(), ref_id = character(), chrom_ref = character(), pos_ref = numeric(), a1_ref = character(), a2_ref = character())

  not_id <- setdiff(seed_dt$input_idx, hit_id$input_idx)
  q_pos <- seed_dt[input_idx %in% not_id & !is.na(chrom) & nzchar(chrom) & is.finite(pos), .(input_idx, chrom, pos)]
  hit_pos <- if (nrow(q_pos)) {
    ref_by_pos[q_pos, on = .(chrom, pos), nomatch = 0L, mult = "first"][
      , .(input_idx, ref_idx, ref_id, chrom_ref = chrom, pos_ref = pos, a1_ref, a2_ref)
    ]
  } else data.table::data.table(input_idx = integer(), ref_idx = integer(), ref_id = character(), chrom_ref = character(), pos_ref = numeric(), a1_ref = character(), a2_ref = character())

  hit_all <- data.table::rbindlist(list(hit_id, hit_pos), use.names = TRUE, fill = TRUE)
  hit_all <- hit_all[!duplicated(input_idx)]
  seed_map <- merge(seed_dt, hit_all, by = "input_idx", all.x = TRUE, sort = FALSE)
  seed_map[, markerID := as.character(markerID)]
  seed_map[is.na(markerID) | !nzchar(markerID), markerID := paste0(chrom, ":", as_int(pos))]
  seed_map[, a1_input := toupper(as.character(a1_input))]
  seed_map[is.na(a1_input) | !nzchar(a1_input), a1_input := NA_character_]
  mapped_seed_n <- as_int(sum(!is.na(seed_map$ref_idx)))

  if (identical(mode_use, "set")) {
    .gcanvas_note(
      "gcanvas::calcld",
      sprintf("mode=%s | pairwise=%s | n_input=%d | n_mapped=%d", mode_use, ifelse(pairwise_use, "TRUE", "FALSE"), n_input, mapped_seed_n),
      silent = silent
    )
  }

  # -------- Build locus variant set --------
  ref_keep <- ref_dt[, .(ref_idx, ref_id, chrom, pos, a1_ref, a2_ref)]
  calc_var <- NULL
  if (identical(mode_use, "set")) {
    calc_var <- seed_map[!is.na(ref_idx), .(ref_idx, markerID, a1_input, is_seed = TRUE, input_idx)]
    if (nrow(calc_var)) {
      calc_var <- ref_keep[calc_var, on = .(ref_idx)]
    } else {
      calc_var <- data.table::data.table(
        ref_idx = integer(), ref_id = character(), chrom = character(), pos = numeric(),
        a1_ref = character(), a2_ref = character(), markerID = character(), a1_input = character(),
        is_seed = logical(), input_idx = integer()
      )
    }
  } else {
    win_src <- seed_map[!is.na(chrom_ref) & nzchar(chrom_ref) & is.finite(pos_ref), .(chrom = chrom_ref, pos = pos_ref)]
    if (!nrow(win_src)) {
      win_src <- seed_map[!is.na(chrom) & nzchar(chrom) & is.finite(pos), .(chrom, pos)]
    }
    if (!nrow(win_src)) {
      stop("Could not resolve focal position for flank mode.", call. = FALSE)
    }
    win <- win_src[, .(start = pmax(1, min(pos, na.rm = TRUE) - flank_bp), end = max(pos, na.rm = TRUE) + flank_bp), by = chrom]
    .fmt_bp <- function(v) {
      vv <- suppressWarnings(as.numeric(v))
      vv[!is.finite(vv)] <- NA_real_
      out <- ifelse(is.na(vv), NA_character_, format(round(vv), scientific = FALSE, trim = TRUE, big.mark = ","))
      as.character(out)
    }
    reg_list <- lapply(seq_len(nrow(win)), function(i) {
      ww <- win[i]
      ref_keep[chrom == ww$chrom & pos >= ww$start & pos <= ww$end]
    })
    calc_var <- data.table::rbindlist(reg_list, use.names = TRUE, fill = TRUE)
    calc_var <- unique(calc_var, by = "ref_idx")
    calc_var[, `:=`(markerID = as.character(ref_id), a1_input = NA_character_, is_seed = FALSE, input_idx = NA_integer_)]
    seed_ref <- seed_map[!is.na(ref_idx)][order(input_idx)]
    seed_ref <- seed_ref[!duplicated(ref_idx)]
    if (nrow(seed_ref)) {
      calc_var[seed_ref, on = .(ref_idx), `:=`(markerID = i.markerID, a1_input = i.a1_input, is_seed = TRUE, input_idx = i.input_idx)]
    }
    n_variant_flank <- as_int(sum(!is.na(calc_var$ref_idx)))
    .gcanvas_note(
      "gcanvas::calcld",
      sprintf(
        "mode=%s | pairwise=%s | flank=%s | chrom=%s | start=%s | end=%s | n_variant=%d",
        mode_use,
        ifelse(pairwise_use, "TRUE", "FALSE"),
        .fmt_bp(flank_bp)[1],
        paste(as.character(win$chrom), collapse = ","),
        paste(.fmt_bp(win$start), collapse = ","),
        paste(.fmt_bp(win$end), collapse = ","),
        n_variant_flank
      ),
      silent = silent
    )
  }

  if (!nrow(calc_var)) {
    if (identical(mode_use, "set") && n_input == 2L && !isTRUE(use_matrix)) return(as.numeric(NA))
    if (isTRUE(use_matrix)) {
      return(list(
        R = matrix(numeric(), nrow = 0L, ncol = 0L),
        a1 = if (isTRUE(use_r2)) NULL else stats::setNames(character(), character()),
        meta = .calcld_meta_pack(chrom = character(), start = NA_real_, end = NA_real_, n_locus = 0L, n_mapped_use = 0L)
      ))
    }
    value_col0 <- if (isTRUE(use_r2)) "r2" else "r"
    return(.calcld_ret_table(
      .calcld_empty_pair_dt(value_col0),
      a1_vec = .calcld_a1_named(NULL, NULL),
      meta = .calcld_meta_pack(chrom = character(), start = NA_real_, end = NA_real_, n_locus = 0L, n_mapped_use = 0L)
    ))
  }

  data.table::setorder(calc_var, input_idx, chrom, pos, ref_idx)
  calc_var <- calc_var[!duplicated(ref_idx)]
  calc_var[, markerID := as.character(markerID)]
  calc_var[is.na(markerID) | !nzchar(markerID), markerID := as.character(ref_id)]
  calc_var[is.na(markerID) | !nzchar(markerID), markerID := paste0(chrom, ":", as_int(pos))]
  calc_var[, markerID := make.unique(markerID)]

  # Rule 8: multiple input but only <=2 mapped variants -> return NA for multi-variant output.
  if (!isTRUE(use_matrix) && identical(mode_use, "set") && n_input > 2L && mapped_seed_n <= 2L) {
    seed_out <- seed_map[order(input_idx), .(
      markerID = markerID,
      chrom = data.table::fifelse(!is.na(chrom) & nzchar(chrom), chrom, chrom_ref),
      pos = data.table::fifelse(is.finite(pos), pos, pos_ref),
      a1_out = a1_input
    )]
    seed_out <- seed_out[!is.na(markerID) & nzchar(markerID)]
    seed_out <- seed_out[!duplicated(markerID)]
    if (!nrow(seed_out)) {
      value_col0 <- if (isTRUE(use_r2)) "r2" else "r"
      return(.calcld_ret_table(
        .calcld_empty_pair_dt(value_col0),
        a1_vec = .calcld_a1_named(seed_out, NULL),
        meta = .calcld_meta_pack(chrom = character(), start = NA_real_, end = NA_real_, n_locus = 0L, n_mapped_use = mapped_seed_n)
      ))
    }
    value_col0 <- if (isTRUE(use_r2)) "r2" else "r"
    out0 <- .calcld_pair_template(seed_out, pairwise_use = pairwise_use, value_col = value_col0)
    return(.calcld_ret_table(
      out0,
      a1_vec = .calcld_a1_named(seed_out, out0),
      meta = .calcld_meta_pack(
        chrom = unique(seed_out$chrom),
        start = suppressWarnings(min(seed_out$pos, na.rm = TRUE)),
        end = suppressWarnings(max(seed_out$pos, na.rm = TRUE)),
        n_locus = nrow(seed_out),
        n_mapped_use = mapped_seed_n
      )
    ))
  }

  calc_m <- calc_var[!is.na(ref_idx)]
  if (!nrow(calc_m) || nrow(calc_m) < 2L) {
    if (identical(mode_use, "set") && n_input == 2L && !isTRUE(use_matrix)) return(as.numeric(NA))
    if (isTRUE(use_matrix)) {
      nm <- calc_var$markerID
      M <- matrix(NA_real_, nrow = length(nm), ncol = length(nm), dimnames = list(nm, nm))
      return(list(
        R = M,
        a1 = if (isTRUE(use_r2)) NULL else setNames(calc_var$a1_input, calc_var$markerID),
        meta = .calcld_meta_pack(
          chrom = unique(calc_var$chrom),
          start = suppressWarnings(min(calc_var$pos, na.rm = TRUE)),
          end = suppressWarnings(max(calc_var$pos, na.rm = TRUE)),
          n_locus = nrow(calc_var),
          n_mapped_use = nrow(calc_m)
        )
      ))
    }
    value_col0 <- if (isTRUE(use_r2)) "r2" else "r"
    v0 <- calc_var[, .(markerID, chrom, pos, a1_out = a1_input)]
    out0 <- .calcld_pair_template(v0, pairwise_use = pairwise_use, value_col = value_col0)
    return(.calcld_ret_table(
      out0,
      a1_vec = .calcld_a1_named(v0, out0),
      meta = .calcld_meta_pack(
        chrom = unique(calc_var$chrom),
        start = suppressWarnings(min(calc_var$pos, na.rm = TRUE)),
        end = suppressWarnings(max(calc_var$pos, na.rm = TRUE)),
        n_locus = nrow(calc_var),
        n_mapped_use = nrow(calc_m)
      )
    ))
  }

  # -------- Compute LD matrix --------
  if (identical(ld_mode, "rds")) {
    idx <- match(calc_m$ref_id, rownames(mat))
    keep <- !is.na(idx)
    calc_m <- calc_m[keep]
    idx <- idx[keep]
    if (length(idx) < 2L) {
      if (identical(mode_use, "set") && n_input == 2L && !isTRUE(use_matrix)) return(as.numeric(NA))
      value_col0 <- if (isTRUE(use_r2)) "r2" else "r"
      v0 <- calc_var[, .(markerID, chrom, pos, a1_out = a1_input)]
      out0 <- .calcld_pair_template(v0, pairwise_use = pairwise_use, value_col = value_col0)
      return(.calcld_ret_table(
        out0,
        a1_vec = .calcld_a1_named(v0, out0),
        meta = .calcld_meta_pack(
          chrom = unique(calc_var$chrom),
          start = suppressWarnings(min(calc_var$pos, na.rm = TRUE)),
          end = suppressWarnings(max(calc_var$pos, na.rm = TRUE)),
          n_locus = nrow(calc_var),
          n_mapped_use = nrow(calc_m)
        )
      ))
    }
    R <- as.matrix(mat[idx, idx, drop = FALSE])
  } else {
    data.table::setorder(calc_m, ref_idx)
    ridx <- as_int(calc_m$ref_idx)
    pos_sorted <- suppressWarnings(as.numeric(ld_obj$map$physical.pos[ridx]))
    size_pos <- max(pos_sorted, na.rm = TRUE) - min(pos_sorted, na.rm = TRUE) + 1
    if (!is.finite(size_pos) || size_pos <= 0) size_pos <- max(flank_bp * 2, 1)
    th <- .gcanvas_ld_thread_lock(ld_threads)
    on.exit(.gcanvas_ld_thread_restore(th), add = TRUE)
    R_sp <- bigsnpr::snp_cor(
      ld_obj$genotypes,
      ind.col = ridx,
      infos.pos = pos_sorted,
      size = size_pos,
      ncores = ld_threads,
      alpha = 1,
      thr_r2 = 0
    )
    R <- as.matrix(R_sp)
  }
  R <- suppressWarnings(as.matrix(R))
  R[!is.finite(R)] <- NA_real_
  R <- (R + t(R)) / 2
  diag(R) <- 1

  calc_m[, a1_ref := toupper(as.character(a1_ref))]
  calc_m[, a2_ref := toupper(as.character(a2_ref))]
  calc_m[, a1_use := toupper(as.character(a1_input))]
  calc_m[is.na(a1_use) | !nzchar(a1_use), a1_use := a2_ref]

  if (isTRUE(use_r2)) {
    R_use <- R * R
    R_use[!is.finite(R_use)] <- NA_real_
    R_use[R_use < 0] <- 0
    R_use[R_use > 1] <- 1
    calc_m[, a1_out := NA_character_]
  } else {
    sgn <- rep(1, nrow(calc_m))
    has_req <- !is.na(calc_m$a1_input) & nzchar(calc_m$a1_input)
    hit_a2 <- has_req & !is.na(calc_m$a2_ref) & (calc_m$a1_input == calc_m$a2_ref)
    hit_a1 <- has_req & !is.na(calc_m$a1_ref) & (calc_m$a1_input == calc_m$a1_ref)
    sgn[hit_a2] <- 1
    sgn[hit_a1] <- -1
    bad <- has_req & !(hit_a1 | hit_a2)
    sgn[bad] <- NA_real_
    S <- outer(sgn, sgn, `*`)
    R_use <- R * S
    R_use[!is.finite(R_use)] <- NA_real_
    R_use[R_use < -1] <- -1
    R_use[R_use > 1] <- 1
    calc_m[, a1_out := a1_use]
  }

  meta <- .calcld_meta_pack(
    chrom = unique(calc_m$chrom),
    start = suppressWarnings(min(calc_m$pos, na.rm = TRUE)),
    end = suppressWarnings(max(calc_m$pos, na.rm = TRUE)),
    n_locus = nrow(calc_m),
    n_mapped_use = mapped_seed_n
  )

  if (isTRUE(use_matrix)) {
    dimnames(R_use) <- list(calc_m$markerID, calc_m$markerID)
    return(list(
      R = R_use,
      a1 = if (isTRUE(use_r2)) NULL else stats::setNames(calc_m$a1_out, calc_m$markerID),
      meta = meta
    ))
  }

  # Rule 7: set mode + exactly 2 input variants -> numeric output.
  if (identical(mode_use, "set") && n_input == 2L) {
    i2 <- seed_map[order(input_idx)][1:2]
    if (any(is.na(i2$ref_idx))) return(as.numeric(NA))
    ii <- match(i2$ref_idx, calc_m$ref_idx)
    if (any(is.na(ii))) return(as.numeric(NA))
    return(as.numeric(R_use[ii[1], ii[2]]))
  }

  value_col <- if (isTRUE(use_r2)) "r2" else "r"
  if (isTRUE(pairwise_use)) {
    if (nrow(calc_m) < 2L) {
      out0 <- .calcld_empty_pair_dt(value_col = value_col)
      return(.calcld_ret_table(out0, a1_vec = .calcld_a1_named(calc_m, out0), meta = meta))
    }
    cb <- utils::combn(seq_len(nrow(calc_m)), 2L)
    i <- cb[1, ]; j <- cb[2, ]
  } else {
    if (nrow(calc_m) < 2L) {
      out0 <- .calcld_empty_pair_dt(value_col = value_col)
      return(.calcld_ret_table(out0, a1_vec = .calcld_a1_named(calc_m, out0), meta = meta))
    }
    anchor <- seed_map[order(input_idx)][!is.na(ref_idx)][1]
    anchor_idx <- if (nrow(anchor)) match(anchor$ref_idx[1], calc_m$ref_idx) else NA_integer_
    if (is.na(anchor_idx)) {
      v0 <- calc_m[, .(markerID, chrom, pos, a1_out)]
      out0 <- .calcld_pair_template(v0, pairwise_use = FALSE, value_col = value_col)
      return(.calcld_ret_table(out0, a1_vec = .calcld_a1_named(v0, out0), meta = meta))
    }
    j <- setdiff(seq_len(nrow(calc_m)), anchor_idx)
    if (!length(j)) {
      out0 <- .calcld_empty_pair_dt(value_col = value_col)
      return(.calcld_ret_table(out0, a1_vec = .calcld_a1_named(calc_m, out0), meta = meta))
    }
    i <- rep(anchor_idx, length(j))
  }

  out <- data.table::data.table(
    markerID = calc_m$markerID[i],
    chrom = calc_m$chrom[i],
    pos = as.numeric(calc_m$pos[i]),
    markerID2 = calc_m$markerID[j]
  )
  out[, (value_col) := as.numeric(R_use[cbind(i, j)])]
  data.table::setcolorder(out, c("markerID", "chrom", "pos", "markerID2", value_col))
  .calcld_ret_table(out[], a1_vec = .calcld_a1_named(calc_m[, .(markerID, a1_out)], out), meta = meta)
}

#' LD-based clumping of association summary statistics
#'
#' Performs the standard PLINK-style clumping algorithm: while there are
#' variants below `p.threshold`, pick the smallest-p one as an index variant
#' and absorb everything within `window` bp at r-squared >= `rsq` into its clump.
#'
#' @param data A `data.frame`/`data.table` of summary statistics with SNP,
#'   chromosome, position, and p-value columns.
#' @param rsq Minimum r-squared for a variant to be clumped into an index.
#' @param p.threshold Maximum p-value for an index variant.
#' @param window Search window. Numeric (bp) or a string like `"500kb"`.
#' @param ld.rds,ld.bfile Reference panel pointers (RDS or PLINK bfile).
#' @param snp.col,chrom.col,pos.col,p.col,z.col Column names in `data`.
#' @param prune.window,prune.step Window/step for variant pruning mode (bp/SNP).
#' @param mode One of `"auto"`, `"clump"`, `"prune"` — algorithm selector.
#' @param ld.cache.dir Optional directory for the LD reference cache.
#' @param threads Integer thread count.
#' @param silent Logical. Suppress progress notes.
#'
#' @return A `data.table` of clumps with the index variant and a list-column
#'   of clumped variants.
#' @export
ldclump <- function(data,
                    rsq = 0.2,
                    p.threshold = 5e-8,
                    window = "500kb",
                    ld.rds = NULL,
                    ld.bfile = NULL,
                    snp.col = "SNP",
                    chrom.col = "CHR",
                    pos.col = "POS",
                    p.col = NULL,
                    z.col = NULL,
                    prune.window = 200L,
                    prune.step = 50L,
                    mode = c("auto", "clump", "prune"),
                    ld.cache.dir = NULL,
                    threads = 4L,
                    silent = FALSE) {
  require_pkg(c("data.table"))
  if (!is.data.frame(data) && !data.table::is.data.table(data)) {
    stop("data must be a data.frame/data.table.", call. = FALSE)
  }
  mode <- match.arg(mode)
  silent <- isTRUE(silent)
  rsq_cut <- suppressWarnings(as.numeric(rsq))[1]
  if (!is.finite(rsq_cut) || is.na(rsq_cut) || rsq_cut < 0 || rsq_cut > 1) {
    stop("rsq must be a numeric value in [0, 1].", call. = FALSE)
  }
  window_bp <- .gcanvas_parse_bp_span(window, arg_name = "window")
  if (!is.finite(window_bp) || is.na(window_bp) || window_bp <= 0) {
    stop("window must be a positive bp span.", call. = FALSE)
  }
  ld_threads <- .gcanvas_resolve_threads(threads)
  prune.window <- as_int(prune.window); if (is.na(prune.window) || prune.window < 2L) prune.window <- 200L
  prune.step <- as_int(prune.step); if (is.na(prune.step) || prune.step < 1L) prune.step <- 50L
  snp_col_use <- .gcanvas_resolve_colname(names(data), snp.col, aliases = character(), required = TRUE, arg_label = "snp.col")
  chrom_col_use <- .gcanvas_resolve_colname(names(data), chrom.col, aliases = character(), required = TRUE, arg_label = "chrom.col")
  pos_col_use <- .gcanvas_resolve_colname(names(data), pos.col, aliases = character(), required = TRUE, arg_label = "pos.col")
  cols <- c(snp_col_use, chrom_col_use, pos_col_use)
  if (!is.null(p.col) && nzchar(as.character(p.col)[1])) cols <- c(cols, as.character(p.col)[1])
  if (!is.null(z.col) && nzchar(as.character(z.col)[1])) cols <- c(cols, as.character(z.col)[1])
  cols <- unique(cols)
  miss <- setdiff(cols, names(data))
  if (length(miss)) stop("Missing columns in data: ", paste(miss, collapse = ", "), call. = FALSE)

  dt_raw <- if (data.table::is.data.table(data)) data.table::copy(data) else data.table::as.data.table(data)
  dt <- dt_raw[, ..cols]
  dt[, row_id := seq_len(.N)]
  data.table::setnames(dt, snp_col_use, "markerID")
  data.table::setnames(dt, chrom_col_use, "chrom")
  data.table::setnames(dt, pos_col_use, "pos")
  if (!is.null(p.col) && (p.col %in% names(dt))) data.table::setnames(dt, p.col, "P_raw")
  if (!is.null(z.col) && (z.col %in% names(dt))) data.table::setnames(dt, z.col, "Z_raw")

  dt[, markerID := as.character(markerID)]
  dt[, chrom := normalize.chrom(chrom)]
  dt[, pos := suppressWarnings(as.numeric(pos))]
  dt <- dt[!is.na(markerID) & nzchar(markerID) & !is.na(chrom) & nzchar(chrom) & is.finite(pos)]
  .ldclump_ret <- function(kept_dt = NULL, variant_obj = NULL, mode_used = mode, prune_out = NULL) {
    if (is.null(kept_dt) || !nrow(kept_dt)) {
      out0 <- dt_raw[0]
    } else {
      out0 <- dt_raw[kept_dt$row_id]
    }
    if (is.null(variant_obj)) {
      variant_obj <- if (identical(mode_used, "prune")) character() else list()
    }
    if (identical(mode_used, "prune")) {
      if (is.null(prune_out)) prune_out <- character()
      out_list <- list(
        data = out0,
        prune.in = variant_obj,
        prune.out = prune_out,
        mode = mode_used
      )
    } else {
      out_list <- list(
        data = out0,
        variant = variant_obj,
        mode = mode_used
      )
    }
    out_list
  }
  if (!nrow(dt)) return(.ldclump_ret(NULL, NULL, mode, character()))
  dt <- dt[!duplicated(markerID)]

  have_p <- "P_raw" %in% names(dt)
  have_z <- "Z_raw" %in% names(dt)
  use_mode <- mode
  if (identical(use_mode, "auto")) {
    use_mode <- if (have_p || have_z) "clump" else "prune"
  }
  if (identical(use_mode, "clump")) {
    if (!have_p && !have_z) {
      .gcanvas_warn_msg("Neither p.col nor z.col is provided; switching to pruning mode.")
      use_mode <- "prune"
    }
  }
  .gcanvas_note(
    "gcanvas::ldclump",
    sprintf("mode=%s | rsq=%.4f | p.threshold=%s | window=%s | prune.window=%d | prune.step=%d | threads=%d",
            as.character(use_mode)[1], rsq_cut, as.character(p.threshold)[1], as.character(window)[1],
            prune.window, prune.step, ld_threads),
    silent = silent
  )
  .gcanvas_note("gcanvas::ldclump", sprintf("n_variants=%d", nrow(dt)), silent = silent)
  do_prefilter <- have_p || have_z
  if (do_prefilter) {
    if (have_p) {
      p_chr <- .gcanvas_p_filter(dt$P_raw)
      y <- -log10c(p_chr)
      thr <- suppressWarnings(as.numeric(p.threshold))[1]
      if (is.finite(thr) && !is.na(thr) && thr > 0 && thr < 1) {
        y_thr <- -log10(thr)
      } else if (is.finite(thr) && !is.na(thr) && thr > 0) {
        y_thr <- thr
      } else {
        y_thr <- -log10(5e-8)
      }
      dt[, score := y]
      dt <- dt[is.finite(score) & !is.na(score) & score >= y_thr]
      .gcanvas_note("gcanvas::ldclump", sprintf("Filtered variants by p-value threshold: n=%d", nrow(dt)), silent = silent)
    } else {
      z0 <- abs(suppressWarnings(as.numeric(dt$Z_raw)))
      thr <- suppressWarnings(as.numeric(p.threshold))[1]
      if (is.finite(thr) && !is.na(thr) && thr > 0 && thr < 1) {
        z_cut <- stats::qnorm(log(thr) - log(2), lower.tail = FALSE, log.p = TRUE)
      } else {
        z_cut <- -Inf
      }
      dt[, score := z0]
      dt <- dt[is.finite(score) & !is.na(score) & score >= z_cut]
      .gcanvas_note("gcanvas::ldclump", sprintf("Filtered variants by z-value threshold: n=%d", nrow(dt)), silent = silent)
    }
    if (!nrow(dt)) return(.ldclump_ret(NULL, NULL, use_mode, character()))
  } else {
    dt[, score := NA_real_]
    if (identical(use_mode, "prune")) {
      .gcanvas_note("gcanvas::ldclump", sprintf("Pruning mode (no p/z pre-filter): n=%d", nrow(dt)), silent = silent)
    }
  }
  dt[, chr_order := rank.chrom(chrom)]

  ld_rds0 <- NULL
  if (!is.null(ld.rds) && length(ld.rds) > 0L && !is.na(ld.rds[1]) && nzchar(as.character(ld.rds[1]))) {
    ld_rds0 <- abs_path(as.character(ld.rds)[1])
    if (!file.exists(ld_rds0)) {
      .gcanvas_warn_msg(paste0("ld.rds not found: ", ld_rds0, " -> trying ld.bfile"))
      ld_rds0 <- NULL
    }
  }
  ld_bfile0 <- NULL
  if (!is.null(ld.bfile) && length(ld.bfile) > 0L && !is.na(ld.bfile[1]) && nzchar(as.character(ld.bfile[1]))) {
    ld_bfile0 <- .gcanvas_normalize_bfile_prefix(ld.bfile)
    if (!is.null(ld_bfile0) && !.gcanvas_ld_bfile_ready(ld_bfile0)) {
      .gcanvas_warn_msg(paste0("ld.bfile incomplete (.bed/.bim/.fam): ", ld_bfile0))
      ld_bfile0 <- NULL
    }
  }
  if (is.null(ld_rds0) && is.null(ld_bfile0)) stop("Provide valid ld.rds or ld.bfile.", call. = FALSE)

  # Build reference index table.
  if (!is.null(ld_rds0)) {
    .gcanvas_note("gcanvas::ldclump", "Loading LD reference from ld.rds", silent = silent)
    obj <- readRDS(ld_rds0)
    mat <- .gcanvas_ld_extract_named_matrix(obj, max_depth = 8L)
    if (is.null(mat)) stop("No named square numeric LD matrix found in ld.rds.", call. = FALSE)
    rn <- rownames(mat); cn <- colnames(mat)
    if (!identical(rn, cn)) {
      shared <- intersect(rn, cn)
      if (!length(shared)) stop("ld.rds matrix row/col names do not overlap.", call. = FALSE)
      mat <- mat[shared, shared, drop = FALSE]
      rn <- rownames(mat)
    }
    ref_dt <- data.table::data.table(ref_id = rn)
    sp <- strsplit(ref_dt$ref_id, ":", fixed = TRUE)
    ref_dt[, chrom := normalize.chrom(vapply(sp, function(v) if (length(v) >= 1L) v[1] else NA_character_, character(1)))]
    ref_dt[, pos := suppressWarnings(as.numeric(vapply(sp, function(v) if (length(v) >= 2L) v[2] else NA_character_, character(1))))]
    ref_dt <- ref_dt[!is.na(chrom) & nzchar(chrom) & is.finite(pos)]
    ref_dt[, ref_idx := match(ref_id, rn)]
  } else {
    cache_dir <- {
      x0 <- as.character(ld.cache.dir)[1]
      if (is.null(ld.cache.dir) || length(ld.cache.dir) == 0L || is.na(x0) || !nzchar(x0) || tolower(x0) == "auto") {
        .gcanvas_default_ld_cache_dir(ld_bfile0)
      } else {
        abs_path(ld.cache.dir)
      }
    }
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    .gcanvas_note("gcanvas::ldclump", sprintf("Attaching LD reference from ld.bfile (cache: %s)", cache_dir), silent = silent)
    ld_obj <- attach_ld_ref(bfile = ld_bfile0, ld_rds = NULL, cache_dir = cache_dir)
    map <- ld_obj$map
    ref_dt <- data.table::data.table(
      ref_idx = seq_len(nrow(map)),
      chrom = normalize.chrom(map$chromosome),
      pos = suppressWarnings(as.numeric(map$physical.pos)),
      a1 = toupper(as.character(map$allele1)),
      a2 = toupper(as.character(map$allele2))
    )
    if ("marker.ID" %in% names(map)) {
      mid <- as.character(map$marker.ID)
      mid[is.na(mid) | !nzchar(mid)] <- NA_character_
      ref_dt[, ref_id := mid]
    } else {
      ref_dt[, ref_id := NA_character_]
    }
    ref_dt[is.na(ref_id) | !nzchar(ref_id), ref_id := paste0(chrom, ":", as_int(round(pos)), ":", a1, ":", a2)]
    ref_dt <- ref_dt[!is.na(chrom) & nzchar(chrom) & is.finite(pos) & !is.na(ref_id) & nzchar(ref_id)]
  }
  if (!nrow(ref_dt)) return(.ldclump_ret(NULL, NULL, use_mode, character()))
  data.table::setkey(ref_dt, chrom, pos)

  # Map input variants to LD reference.
  map_by_id <- ref_dt[dt[, .(markerID)], on = .(ref_id = markerID), nomatch = 0L][, .(markerID = ref_id, ref_idx)]
  miss_id <- setdiff(dt$markerID, map_by_id$markerID)
  map_by_pos <- data.table::data.table()
  if (length(miss_id)) {
    miss_dt <- dt[markerID %in% miss_id, .(markerID, chrom, pos)]
    map_by_pos <- ref_dt[miss_dt, on = .(chrom, pos), nomatch = 0L, mult = "first"][, .(markerID = i.markerID, ref_idx)]
  }
  map_dt <- unique(data.table::rbindlist(list(map_by_id, map_by_pos), use.names = TRUE, fill = TRUE), by = "markerID")
  dt <- merge(dt, map_dt, by = "markerID", all.x = FALSE, sort = FALSE)
  .gcanvas_note("gcanvas::ldclump", sprintf("Mapped variants to LD reference: n=%d", nrow(dt)), silent = silent)
  if (!nrow(dt)) return(.ldclump_ret(NULL, NULL, use_mode, character()))

  # Build helper for r2 query between lead and candidates.
  get_r2 <- function(lead_idx, cand_idx) {
    if (!length(cand_idx)) return(rep(NA_real_, 0L))
    if (!is.null(ld_rds0)) {
      r <- as.numeric(mat[lead_idx, cand_idx, drop = TRUE])
      r2 <- r * r
      r2[!is.finite(r2)] <- NA_real_
      pmax(0, pmin(1, r2))
    } else {
      ridx_all <- unique(as_int(c(lead_idx, cand_idx)))
      ridx_all <- ridx_all[is.finite(ridx_all) & !is.na(ridx_all)]
      if (length(ridx_all) < 2L) return(rep(NA_real_, length(cand_idx)))
      pos_all <- suppressWarnings(as.numeric(ld_obj$map$physical.pos[ridx_all]))
      size_pos <- max(pos_all, na.rm = TRUE) - min(pos_all, na.rm = TRUE) + 1
      if (!is.finite(size_pos) || size_pos <= 0) size_pos <- window_bp * 2 + 1
      data.table::setorder(data.table::data.table(ridx = ridx_all), ridx)
      ridx_sorted <- sort(ridx_all)
      lead_local <- match(lead_idx, ridx_sorted)
      if (is.na(lead_local)) return(rep(NA_real_, length(cand_idx)))
      R_sp <- bigsnpr::snp_cor(
        ld_obj$genotypes,
        ind.col = ridx_sorted,
        infos.pos = suppressWarnings(as.numeric(ld_obj$map$physical.pos[ridx_sorted])),
        size = size_pos,
        ncores = ld_threads,
        alpha = 1,
        thr_r2 = 0
      )
      col_pos <- match(cand_idx, ridx_sorted)
      r <- rep(NA_real_, length(cand_idx))
      ok <- !is.na(col_pos)
      if (any(ok)) r[ok] <- as.numeric(R_sp[lead_local, col_pos[ok]])
      r2 <- r * r
      r2[!is.finite(r2)] <- NA_real_
      pmax(0, pmin(1, r2))
    }
  }

  th <- .gcanvas_ld_thread_lock(ld_threads)
  on.exit(.gcanvas_ld_thread_restore(th), add = TRUE)

  # Fast pruning path cache: build chromosome-wise LD edges once for both ld.rds / ld.bfile.
  prune_edge_cache <- list()
  if (identical(use_mode, "prune")) {
    require_pkg("Matrix")
    dt_pr <- data.table::copy(dt)
    data.table::setorderv(dt_pr, c("chr_order", "pos", "markerID"), c(1L, 1L, 1L), na.last = TRUE)
    dt <- dt_pr

    chr_levels <- unique(dt$chrom)
    chr_levels <- chr_levels[!is.na(chr_levels) & nzchar(chr_levels)]
    chr_data <- list()
    for (chr_i in chr_levels) {
      idx_chr <- which(dt$chrom == chr_i)
      if (length(idx_chr) < 2L) next
      ridx_chr <- as_int(dt$ref_idx[idx_chr])
      pos_chr <- as.numeric(dt$pos[idx_chr])
      ok_chr <- is.finite(ridx_chr) & !is.na(ridx_chr) & is.finite(pos_chr) & !is.na(pos_chr)
      idx_chr <- idx_chr[ok_chr]
      ridx_chr <- ridx_chr[ok_chr]
      pos_chr <- pos_chr[ok_chr]
      if (length(idx_chr) < 2L) next
      ord_chr <- order(pos_chr, ridx_chr)
      idx_chr <- idx_chr[ord_chr]
      ridx_chr <- ridx_chr[ord_chr]
      pos_chr <- pos_chr[ord_chr]
      chr_data[[as.character(chr_i)]] <- list(idx = idx_chr, ridx = ridx_chr, pos = pos_chr, n = length(idx_chr))
    }
    chr_keys <- names(chr_data)
    if (length(chr_keys)) {
      workers_base <- 1L
      if (ld_threads > 1L && .Platform$OS.type != "windows" && requireNamespace("parallel", quietly = TRUE)) {
        workers_base <- as_int(ld_threads)
      }
      if (!is.finite(workers_base) || is.na(workers_base) || workers_base < 1L) workers_base <- 1L
      chunk_size <- max(128L, min(1024L, as_int(round(5e5 / max(50, prune.window)))))
      if (!is.finite(chunk_size) || is.na(chunk_size) || chunk_size < 1L) chunk_size <- 512L

      .run_chunked <- function(tasks, FUN, workers, label = "LD prune edge cache") {
        n_task <- length(tasks)
        if (!n_task) return(list())
        out <- vector("list", n_task)
        done <- 0L
        if (workers <= 1L) {
          for (i in seq_len(n_task)) {
            out[[i]] <- FUN(tasks[[i]])
            done <- i
            .gcanvas_note("gcanvas::ldclump", sprintf("%s: chunk %d/%d", label, done, n_task), silent = silent)
          }
          return(out)
        }
        b_st <- seq.int(1L, n_task, by = workers)
        for (b in b_st) {
          e <- min(n_task, b + workers - 1L)
          ids <- b:e
          res_b <- parallel::mclapply(ids, function(k) FUN(tasks[[k]]), mc.cores = length(ids))
          out[ids] <- res_b
          for (k in ids) {
            done <- done + 1L
            .gcanvas_note("gcanvas::ldclump", sprintf("%s: chunk %d/%d", label, done, n_task), silent = silent)
          }
        }
        out
      }

      if (!is.null(ld_rds0)) {
        task_list <- list()
        tid <- 0L
        for (chr_i in chr_keys) {
          n_chr <- as_int(chr_data[[chr_i]]$n)
          starts <- seq.int(1L, n_chr, by = chunk_size)
          for (s in starts) {
            tid <- tid + 1L
            task_list[[tid]] <- list(chr = chr_i, s = as_int(s), e = as_int(min(n_chr, s + chunk_size - 1L)))
          }
        }
        workers_rds <- if (workers_base > 1L) min(workers_base, as_int(length(task_list))) else 1L
        if (!is.finite(workers_rds) || is.na(workers_rds) || workers_rds < 1L) workers_rds <- 1L
        .gcanvas_note(
          "gcanvas::ldclump",
          sprintf("Building LD prune edge cache (source=ld.rds, chunks=%d, workers=%d, chunk.size=%d)", length(task_list), workers_rds, chunk_size),
          silent = silent
        )
        .task_rds <- function(task) {
          cd <- chr_data[[task$chr]]
          if (is.null(cd)) return(NULL)
          ridx_chr <- as_int(cd$ridx)
          pos_chr <- as.numeric(cd$pos)
          n_chr <- as_int(cd$n)
          s <- max(1L, as_int(task$s))
          e <- min(n_chr, as_int(task$e))
          if (!is.finite(s) || !is.finite(e) || s > e) return(NULL)

          li <- vector("list", e - s + 1L)
          lj <- vector("list", e - s + 1L)
          lr <- vector("list", e - s + 1L)
          nbuf <- 0L
          for (i_local in seq.int(s, e)) {
            if (!is.finite(ridx_chr[i_local]) || is.na(ridx_chr[i_local])) next
            j_bp <- findInterval(pos_chr[i_local] + window_bp, pos_chr)
            j_end <- min(n_chr, i_local + prune.window - 1L, j_bp)
            if (!is.finite(j_end) || is.na(j_end) || j_end <= i_local) next
            cand <- seq.int(i_local + 1L, j_end)
            r <- as.numeric(mat[ridx_chr[i_local], ridx_chr[cand], drop = TRUE])
            r2 <- r * r
            sel <- is.finite(r2) & !is.na(r2) & r2 >= rsq_cut
            if (!any(sel)) next
            nbuf <- nbuf + 1L
            li[[nbuf]] <- rep.int(i_local, sum(sel))
            lj[[nbuf]] <- cand[sel]
            lr[[nbuf]] <- r2[sel]
          }
          if (nbuf == 0L) return(NULL)
          data.table::data.table(
            chr = as.character(task$chr),
            i_local = as_int(unlist(li[seq_len(nbuf)], use.names = FALSE)),
            j_local = as_int(unlist(lj[seq_len(nbuf)], use.names = FALSE)),
            rsq = as.numeric(unlist(lr[seq_len(nbuf)], use.names = FALSE))
          )
        }
        res_chunks <- .run_chunked(task_list, .task_rds, workers = workers_rds, label = "LD prune edge cache")
        edge_dt <- data.table::rbindlist(res_chunks, use.names = TRUE, fill = TRUE)
        if (nrow(edge_dt)) {
          for (chr_i in chr_keys) {
            ed <- edge_dt[chr == chr_i]
            if (!nrow(ed)) next
            idx_chr <- as_int(chr_data[[chr_i]]$idx)
            ed[, `:=`(
              i = idx_chr[i_local],
              j = idx_chr[j_local]
            )]
            ed <- ed[is.finite(rsq) & !is.na(rsq) & rsq >= rsq_cut & (j_local - i_local) < prune.window]
            if (!nrow(ed)) next
            data.table::setorderv(ed, c("i_local", "j_local"), c(1L, 1L), na.last = TRUE)
            prune_edge_cache[[chr_i]] <- ed[, .(i, j, i_local, j_local, rsq)]
          }
        }
      } else {
        workers_bfile <- if (workers_base > 1L) min(workers_base, as_int(length(chr_keys))) else 1L
        if (!is.finite(workers_bfile) || is.na(workers_bfile) || workers_bfile < 1L) workers_bfile <- 1L
        .gcanvas_note(
          "gcanvas::ldclump",
          sprintf("Building LD prune edge cache (source=ld.bfile, chunks=%d, workers=%d)", length(chr_keys), workers_bfile),
          silent = silent
        )
        .task_bfile <- function(chr_i) {
          cd <- chr_data[[chr_i]]
          if (is.null(cd)) return(NULL)
          idx_chr <- as_int(cd$idx)
          ridx_chr <- as_int(cd$ridx)
          pos_chr <- as.numeric(cd$pos)
          if (length(idx_chr) < 2L) return(NULL)
          size_pos <- max(pos_chr, na.rm = TRUE) - min(pos_chr, na.rm = TRUE) + 1
          if (!is.finite(size_pos) || size_pos <= 0) size_pos <- window_bp * 2 + 1
          ncores_inner <- if (workers_bfile > 1L) 1L else ld_threads
          R_sp <- bigsnpr::snp_cor(
            ld_obj$genotypes,
            ind.col = ridx_chr,
            infos.pos = pos_chr,
            size = size_pos,
            ncores = as_int(ncores_inner),
            alpha = 1,
            thr_r2 = rsq_cut
          )
          sm <- Matrix::summary(R_sp)
          if (!nrow(sm)) return(NULL)
          sm <- sm[sm$i < sm$j, , drop = FALSE]
          if (!nrow(sm)) return(NULL)
          i_local <- as_int(sm$i)
          j_local <- as_int(sm$j)
          rsqv <- as.numeric(sm$x * sm$x)
          dbp <- pos_chr[j_local] - pos_chr[i_local]
          keep_e <- is.finite(rsqv) & !is.na(rsqv) & rsqv >= rsq_cut &
            (j_local - i_local) < prune.window &
            is.finite(dbp) & !is.na(dbp) & dbp <= window_bp
          if (!any(keep_e)) return(NULL)
          i_local <- i_local[keep_e]
          j_local <- j_local[keep_e]
          rsqv <- rsqv[keep_e]
          data.table::data.table(
            i = idx_chr[i_local],
            j = idx_chr[j_local],
            i_local = i_local,
            j_local = j_local,
            rsq = rsqv
          )
        }
        res_chunks <- .run_chunked(as.list(chr_keys), .task_bfile, workers = workers_bfile, label = "LD prune edge cache")
        for (ii in seq_along(chr_keys)) {
          chr_i <- chr_keys[ii]
          ed <- res_chunks[[ii]]
          if (is.null(ed) || !nrow(ed)) next
          data.table::setorderv(ed, c("i_local", "j_local"), c(1L, 1L), na.last = TRUE)
          prune_edge_cache[[chr_i]] <- ed
        }
      }
      .gcanvas_note("gcanvas::ldclump", sprintf("LD edge cache ready: %d chromosomes", length(prune_edge_cache)), silent = silent)
    }
  }

  ld_row_cache <- new.env(parent = emptyenv())
  get_r2_rows <- function(i_row, cand_rows) {
    cand_rows <- as_int(cand_rows)
    cand_rows <- cand_rows[is.finite(cand_rows) & !is.na(cand_rows) & cand_rows >= 1L & cand_rows <= nrow(dt)]
    if (!length(cand_rows)) return(rep(NA_real_, 0L))
    k <- as.character(as_int(i_row))[1]
    ent <- if (exists(k, envir = ld_row_cache, inherits = FALSE)) get(k, envir = ld_row_cache, inherits = FALSE) else NULL
    if (is.null(ent)) {
      chr_i <- dt$chrom[i_row]
      pos_i <- dt$pos[i_row]
      cand_all <- which(dt$chrom == chr_i & abs(dt$pos - pos_i) <= window_bp)
      r2_all <- get_r2(dt$ref_idx[i_row], dt$ref_idx[cand_all])
      ent <- list(rows = as_int(cand_all), r2 = as.numeric(r2_all))
      assign(k, ent, envir = ld_row_cache)
    }
    m <- match(cand_rows, ent$rows)
    r <- rep(NA_real_, length(cand_rows))
    ok <- !is.na(m)
    if (any(ok)) r[ok] <- ent$r2[m[ok]]
    r
  }

  keep <- rep(TRUE, nrow(dt))
  dep_map <- list()
  .ldclump_add_dep <- function(lead_id, dep_ids, dep_r2) {
    if (is.null(lead_id) || !length(lead_id)) return(invisible(NULL))
    lead_id <- as.character(lead_id)[1]
    dep_ids <- as.character(dep_ids)
    dep_r2 <- suppressWarnings(as.numeric(dep_r2))
    ok <- !is.na(dep_ids) & nzchar(dep_ids) & is.finite(dep_r2) & !is.na(dep_r2)
    if (!any(ok)) {
      if (is.null(dep_map[[lead_id]])) dep_map[[lead_id]] <<- setNames(numeric(), character())
      return(invisible(NULL))
    }
    dep_ids <- dep_ids[ok]
    dep_r2 <- dep_r2[ok]
    addv <- stats::setNames(dep_r2, dep_ids)
    cur <- dep_map[[lead_id]]
    if (is.null(cur)) cur <- setNames(numeric(), character())
    for (j in seq_along(addv)) {
      nm <- names(addv)[j]
      rv <- as.numeric(addv[j])
      if (nm %in% names(cur)) {
        if (is.finite(rv) && (is.na(cur[[nm]]) || rv > cur[[nm]])) cur[[nm]] <- rv
      } else {
        cur <- c(cur, stats::setNames(rv, nm))
      }
    }
    dep_map[[lead_id]] <<- cur
    invisible(NULL)
  }
  if (identical(use_mode, "clump")) {
    .gcanvas_note("gcanvas::ldclump", "Running LD clumping loop", silent = silent)
    data.table::setorderv(dt, c("chr_order", "score", "pos", "markerID"), c(1L, -1L, 1L, 1L), na.last = TRUE)
    for (i in seq_len(nrow(dt))) {
      if (!keep[i]) next
      chr_i <- dt$chrom[i]
      pos_i <- dt$pos[i]
      cand <- which(keep & dt$chrom == chr_i & abs(dt$pos - pos_i) <= window_bp)
      if (!length(cand)) next
      r2 <- get_r2_rows(i, cand)
      drop <- cand[is.finite(r2) & r2 >= rsq_cut]
      drop <- setdiff(drop, i)
      if (length(drop)) {
        keep[drop] <- FALSE
        lead_id <- dt$markerID[i]
        r2_drop <- r2[match(drop, cand)]
        .ldclump_add_dep(lead_id, dt$markerID[drop], r2_drop)
      }
      if (is.null(dep_map[[dt$markerID[i]]])) dep_map[[dt$markerID[i]]] <- setNames(numeric(), character())
    }
  } else {
    chr_all <- unique(as.character(dt$chrom))
    chr_all <- chr_all[!is.na(chr_all) & nzchar(chr_all)]
    workers_prune <- 1L
    if (ld_threads > 1L && .Platform$OS.type != "windows" && requireNamespace("parallel", quietly = TRUE)) {
      workers_prune <- min(as_int(ld_threads), as_int(length(chr_all)))
    }
    if (!is.finite(workers_prune) || is.na(workers_prune) || workers_prune < 1L) workers_prune <- 1L
    .gcanvas_note(
      "gcanvas::ldclump",
      sprintf("Running LD pruning loop | chromosomes=%d | workers=%d", as_int(length(chr_all)), as_int(workers_prune)),
      silent = silent
    )
    # PLINK-like indep-pairwise controls:
    # prune.window (variant count), prune.step (step), rsq threshold.
    data.table::setorderv(dt, c("chr_order", "pos", "markerID"), c(1L, 1L, 1L), na.last = TRUE)
    .prune_chr <- function(chr_i) {
      idx_chr <- which(dt$chrom == chr_i)
      if (!length(idx_chr)) return(list(chr = chr_i, keep_local = logical(0)))
      n_chr <- length(idx_chr)
      keep_local <- rep(TRUE, n_chr)
      if (n_chr < 2L) return(list(chr = chr_i, keep_local = keep_local))
      starts <- seq.int(1L, n_chr, by = prune.step)
      ed_chr <- prune_edge_cache[[as.character(chr_i)]]
      edge_by_i <- NULL
      if (!is.null(ed_chr) && nrow(ed_chr)) {
        if (!("i_local" %in% names(ed_chr)) || !("j_local" %in% names(ed_chr))) {
          map_local <- integer(nrow(dt))
          map_local[idx_chr] <- seq_len(n_chr)
          ed_chr[, i_local := map_local[i]]
          ed_chr[, j_local := map_local[j]]
        }
        ed_chr <- ed_chr[is.finite(rsq) & !is.na(rsq) & i_local >= 1L & j_local > i_local]
        if (nrow(ed_chr)) {
          data.table::setorderv(ed_chr, c("i_local", "j_local"), c(1L, 1L), na.last = TRUE)
          edge_by_i <- split(seq_len(nrow(ed_chr)), ed_chr$i_local)
        }
      }
      processed_until <- integer(n_chr)
      for (s in starts) {
        e <- min(n_chr, s + prune.window - 1L)
        if (s >= e) next
        for (i_local in seq.int(s, e)) {
          if (!keep_local[i_local]) next
          prev_u <- processed_until[i_local]
          if (e <= prev_u) next
          drop_local <- integer()
          if (!is.null(edge_by_i)) {
            hit_idx <- edge_by_i[[as.character(i_local)]]
            if (!is.null(hit_idx) && length(hit_idx)) {
              jloc <- as_int(ed_chr$j_local[hit_idx])
              sel_j <- jloc > max(i_local, prev_u) & jloc <= e
              if (any(sel_j)) {
                cand_loc <- jloc[sel_j]
                if (length(cand_loc)) {
                  cand_loc <- cand_loc[is.finite(cand_loc) & !is.na(cand_loc) & cand_loc >= 1L & cand_loc <= n_chr]
                  if (length(cand_loc)) {
                    cand_loc <- unique(cand_loc)
                    drop_local <- cand_loc[keep_local[cand_loc]]
                  }
                }
              }
            }
          }
          processed_until[i_local] <- max(processed_until[i_local], e)
          if (length(drop_local)) keep_local[drop_local] <- FALSE
        }
      }
      list(chr = chr_i, keep_local = keep_local)
    }

    keep[] <- FALSE
    if (workers_prune <= 1L) {
      res_pr <- lapply(chr_all, .prune_chr)
    } else {
      res_pr <- parallel::mclapply(chr_all, .prune_chr, mc.cores = workers_prune)
    }
    for (res in res_pr) {
      if (is.null(res) || is.null(res$chr) || is.null(res$keep_local)) next
      idx_chr <- which(dt$chrom == res$chr)
      if (!length(idx_chr)) next
      kk <- as.logical(res$keep_local)
      if (length(kk) != length(idx_chr)) {
        kk <- rep_len(kk, length(idx_chr))
      }
      keep[idx_chr] <- kk
    }
  }

  kept <- dt[keep]
  .gcanvas_note("gcanvas::ldclump", sprintf("Done: kept=%d, dropped=%d", nrow(kept), nrow(dt) - nrow(kept)), silent = silent)
  if (identical(use_mode, "prune")) {
    variant_keep <- unique(as.character(dt$markerID[keep]))
    variant_keep <- variant_keep[!is.na(variant_keep) & nzchar(variant_keep)]
    prune_out <- unique(as.character(dt$markerID[!keep]))
    prune_out <- prune_out[!is.na(prune_out) & nzchar(prune_out)]
    if (!nrow(kept)) return(.ldclump_ret(NULL, variant_keep, use_mode, prune_out))
    return(.ldclump_ret(kept, variant_keep, use_mode, prune_out))
  }

  if (!nrow(kept)) return(.ldclump_ret(NULL, dep_map, use_mode))
  .ldclump_ret(kept, dep_map, use_mode)
}

