# LD-reference cache inspectors. The cache itself lives inside the LD backend
# helpers; these are the user-facing entry points.

#' List currently cached LD reference entries
#'
#' Returns counts of items held in the in-memory and on-disk LD reference
#' caches used by [calcld()] / [ldproxy()] / [ldclump()].
#'
#' @return A named list of cache counts.
#' @export
cache.list <- function() {
  .gcanvas_cache_counts()
}

#' Reset the LD reference cache
#'
#' Clears both the in-memory and on-disk caches populated by the LD helpers.
#'
#' @param silent Logical. Suppress progress notes.
#' @return Invisibly, a list summarizing how many entries were removed.
#' @export
cache.reset <- function(silent = FALSE) {
  silent <- isTRUE(silent)

  # 1) Memory cache
  mem_removed <- 0L
  e <- tryCatch(.gcanvas_ld_mem_cache, error = function(err) NULL)
  if (is.environment(e) && is.environment(e$entries)) {
    keys <- ls(envir = e$entries, all.names = TRUE)
    mem_removed <- as_int(length(keys))
    if (length(keys)) rm(list = keys, envir = e$entries)
    e$order <- character(0)
    e$total_bytes <- 0
  }

  # 2) Disk cache files created in this session
  reg <- .gcanvas_session_cache_registry
  created_files <- unique(as.character(reg$created_files))
  created_files <- created_files[!is.na(created_files) & nzchar(created_files)]
  files_exist <- created_files[file.exists(created_files)]
  if (length(files_exist)) suppressWarnings(unlink(files_exist, recursive = FALSE, force = TRUE))
  files_removed <- if (length(files_exist)) as_int(sum(!file.exists(files_exist))) else 0L

  # 3) Leftover lock files in created dirs
  created_dirs <- unique(as.character(reg$created_dirs))
  created_dirs <- created_dirs[!is.na(created_dirs) & nzchar(created_dirs)]
  if (length(created_files)) {
    file_dirs <- unique(dirname(created_files))
    file_dirs <- file_dirs[!is.na(file_dirs) & nzchar(file_dirs)]
    created_dirs <- unique(c(created_dirs, file_dirs))
  }
  lock_files <- character(0)
  if (length(created_dirs)) {
    for (d in created_dirs) {
      if (!dir.exists(d)) next
      lf <- list.files(d, pattern = "\\.ldref\\.lock$", full.names = TRUE)
      if (length(lf)) lock_files <- c(lock_files, lf)
    }
  }
  lock_files <- unique(lock_files[file.exists(lock_files)])
  if (length(lock_files)) suppressWarnings(unlink(lock_files, recursive = FALSE, force = TRUE))
  lock_removed <- if (length(lock_files)) as_int(sum(!file.exists(lock_files))) else 0L

  # 4) Remove now-empty cache dirs that were created in this session
  dirs_removed <- 0L
  if (length(created_dirs)) {
    created_dirs <- created_dirs[order(nchar(created_dirs), decreasing = TRUE)]
    for (d in created_dirs) {
      if (!dir.exists(d)) next
      n_inside <- length(list.files(d, all.files = TRUE, no.. = TRUE))
      if (n_inside == 0L) {
        ok <- suppressWarnings(unlink(d, recursive = TRUE, force = TRUE))
        if (identical(ok, 0L) || !dir.exists(d)) dirs_removed <- dirs_removed + 1L
      }
    }
  }

  # 5) Keep only remaining tracked paths
  reg$created_files <- created_files[file.exists(created_files)]
  reg$created_dirs <- created_dirs[dir.exists(created_dirs)]

  removed <- list(
    memory = list(ld_entries_removed = as_int(mem_removed)),
    disk = list(
      files_removed = as_int(files_removed),
      lock_files_removed = as_int(lock_removed),
      dirs_removed = as_int(dirs_removed)
    )
  )
  current <- .gcanvas_cache_counts()

  if (!isTRUE(silent)) {
    .gcanvas_note("gcanvas::cache.reset", sprintf(
      "memory_entries_removed=%d | disk_files_removed=%d | lock_removed=%d | dirs_removed=%d",
      as_int(removed$memory$ld_entries_removed),
      as_int(removed$disk$files_removed),
      as_int(removed$disk$lock_files_removed),
      as_int(removed$disk$dirs_removed)
    ), silent = FALSE)
  }
  invisible(list(removed = removed, current = current))
}

