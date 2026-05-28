# Package-internal helpers: dependency checks, coercion, path normalization,
# input parsing, seed handling, and the shared message sink used by all
# verbose-progress callers (`.gcanvas_note`).

require_pkg <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Missing packages: ", paste(miss, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

as_int <- function(x) suppressWarnings(as.integer(x))
as_num <- function(x) suppressWarnings(as.numeric(x))
abs_path <- function(x) {
  if (is.null(x) || length(x) == 0) return(x)
  x <- as.character(x)
  x[is.na(x) | !nzchar(x)] <- NA_character_
  y <- path.expand(x)
  is_abs <- grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)", y)
  y[!is_abs & !is.na(y)] <- file.path(getwd(), y[!is_abs & !is.na(y)])
  normalizePath(y, winslash = "/", mustWork = FALSE)
}

.gcanvas_as_snp_vector <- function(x) {
  if (is.null(x) || length(x) == 0L) return(character())
  if (is.data.frame(x) || data.table::is.data.table(x)) {
    if (!nrow(x) || !ncol(x)) return(character())
    out <- as.character(x[[1L]])
    out[!is.na(out) & nzchar(out)]
  } else if (is.list(x)) {
    out <- unlist(x, recursive = TRUE, use.names = TRUE)
    out <- as.character(out)
    out[!is.na(out) & nzchar(out)]
  } else {
    if (is.character(x) && length(x) == 1L && !is.na(x[1]) && nzchar(x[1]) && file.exists(x[1])) {
      out <- tryCatch({
        dt0 <- data.table::fread(x[1], data.table = FALSE, showProgress = FALSE)
        if (!nrow(dt0) || !ncol(dt0)) character() else as.character(dt0[[1L]])
      }, error = function(e) {
        readLines(x[1], warn = FALSE)
      })
    } else {
      out <- as.character(x)
    }
    out[!is.na(out) & nzchar(out)]
  }
}

.gcanvas_seed_resolve <- function(seed = NULL) {
  if (is.null(seed) || length(seed) == 0L) return(NULL)
  s0 <- seed[1]
  if (is.character(s0)) {
    key <- tolower(trimws(as.character(s0)))
    if (!is.na(key) && key %in% c("random", "rand")) {
      return(as_int(sample.int(.Machine$integer.max - 1L, size = 1L)))
    }
  }
  s <- as_int(s0)
  if (!is.finite(s) || is.na(s)) return(NULL)
  s <- abs(s)
  if (!is.finite(s) || is.na(s) || s == 0L) s <- 1L
  as_int(s)
}

.gcanvas_seed_label <- function(seed = NULL, seed_use = NULL) {
  if (is.null(seed) || length(seed) == 0L) {
    return("NULL (session RNG state; not reset)")
  }
  s0 <- seed[1]
  if (is.character(s0)) {
    key <- tolower(trimws(as.character(s0)))
    if (!is.na(key) && key %in% c("random", "rand")) {
      if (!is.null(seed_use) && is.finite(seed_use) && !is.na(seed_use)) {
        return(sprintf("random -> %d", as_int(seed_use)))
      }
      return("random")
    }
  }
  s_num <- as_int(s0)
  if (is.finite(s_num) && !is.na(s_num)) {
    return(sprintf("%d", as_int(abs(s_num))))
  }
  "invalid -> NULL (session RNG state; not reset)"
}

.gcanvas_msg <- function(..., sep = "", appendLF = TRUE) {
  cat(..., sep = sep)
  if (isTRUE(appendLF)) cat("\n")
  invisible(NULL)
}

.gcanvas_warn_msg <- function(..., sep = "", appendLF = TRUE) {
  cat("[WARN] ", ..., sep = sep)
  if (isTRUE(appendLF)) cat("\n")
  invisible(NULL)
}

.gcanvas_note <- function(module = "gcanvas", msg = NULL, silent = FALSE) {
  if (isTRUE(silent)) return(invisible(NULL))
  msg0 <- as.character(msg)[1]
  if (is.na(msg0) || !nzchar(msg0)) return(invisible(NULL))
  module0 <- as.character(module)[1]
  if (is.na(module0) || !nzchar(module0)) module0 <- "gcanvas"
  .gcanvas_msg(sprintf("[%s] %s", module0, msg0))
  invisible(NULL)
}
