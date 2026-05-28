# Lightweight logger used by long-running gcanvas calls. Returns a list of
# closures (info / warn / error / debug / flat / close) sharing one optional
# file connection.

#' Create a logger object
#'
#' Returns a list of logging closures (`info`, `warn`, `error`, `debug`,
#' `flat`, `close`) that print to the console and optionally tee to a log
#' file. Timestamps and verbosity are configurable.
#'
#' @param when Logical, `"now"`, or a `POSIXt` object. Controls whether
#'   messages are timestamped, and pins the timestamp if a `POSIXt` is given.
#' @param log_path Optional path to a log file; parent directories are created.
#' @param append Logical. Append to an existing log file instead of truncating.
#' @param verbose Logical. If `FALSE`, suppress `DEBUG` console messages.
#'
#' @return A named list of logger closures.
#' @export
logger <- function(when = TRUE, log_path = NULL, append = FALSE, verbose = TRUE) {
  con <- NULL
  use_ts <- FALSE
  ts_fixed <- NULL

  if (inherits(when, "POSIXt")) {
    use_ts <- TRUE
    ts_fixed <- format(when, "%Y-%m-%d %H:%M:%S")
  } else if (isTRUE(when) || identical(tolower(as.character(when)[1]), "now")) {
    use_ts <- TRUE
  } else {
    use_ts <- FALSE
  }

  if (!is.null(log_path) && nzchar(log_path)) {
    dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)
    if (isTRUE(append) && file.exists(log_path)) {
      con <- file(log_path, open = "at")
      writeLines("----------------------------------------------------------------", con)
      writeLines("", con)
    } else {
      con <- file(log_path, open = "wt")
    }
  }

  logf <- function(level, msg) {
    ts <- if (use_ts) (ts_fixed %||% format(Sys.time(), "%Y-%m-%d %H:%M:%S")) else NULL
    line <- if (!is.null(ts)) sprintf("[%s] [%s] %s", ts, level, msg) else sprintf("[%s] %s", level, msg)

    if (!is.null(con)) {
      writeLines(line, con)
      flush(con)
    }
    if (isTRUE(verbose) || level != "DEBUG") message(line)
    invisible(line)
  }

  list(
    info  = function(msg)  logf("INFO", msg),
    warn  = function(msg)  logf("WARN", msg),
    error = function(msg)  logf("ERROR", msg),
    debug = function(msg)  logf("DEBUG", msg),
    flat  = function(msg)  if (!is.null(con)) { writeLines(msg, con); flush(con) } else if (isTRUE(verbose)) message(msg),
    close = function() { if (!is.null(con)) tryCatch(close(con), error = function(e) invisible(NULL)) }
  )
}
