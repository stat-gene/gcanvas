# Statistical helpers used across the package: z-scores, p-values, credible
# sets, quantile-based binning, and inverse-normal transformation.

#' Absolute z-score from beta and SE (or variance)
#'
#' Computes `|beta / se|` (or `|beta / sqrt(varbeta)|`) elementwise, returning
#' `NA_real_` where inputs are missing, non-finite, or have a non-positive
#' standard error.
#'
#' @param beta Numeric vector of effect sizes.
#' @param se Numeric vector of standard errors (mutually exclusive with `varbeta`).
#' @param varbeta Numeric vector of variances of `beta` (mutually exclusive with `se`).
#'
#' @return Numeric vector of absolute z-scores, same length as `beta`.
#' @export
zabs <- function(beta, se = NULL, varbeta = NULL) {
  beta <- as_num(beta)
  se_use <- NULL
  if (!is.null(varbeta)) {
    vb <- as_num(varbeta)
    vb[!is.na(vb) & vb <= 0] <- NA_real_
    se_use <- sqrt(vb)
  } else if (!is.null(se)) {
    se_use <- as_num(se)
  }
  if (is.null(se_use)) return(rep(NA_real_, length(beta)))
  se_use[!is.na(se_use) & se_use <= 0] <- NA_real_
  z <- suppressWarnings(beta / se_use)
  out <- abs(z)
  out[!is.finite(out)] <- NA_real_
  out
}

#' `10^x` as a scientific-notation character string
#'
#' Inverse companion to [log10c()]: given a (possibly very negative) `log10`
#' value, returns its scientific-notation character form without underflow.
#' For example, `pow10c(-12012.3)` returns `"5.012e-12013"` even though the
#' actual numeric value would underflow to `0` in IEEE-754 doubles.
#'
#' @param x Numeric vector of `log10` values.
#' @param digits Integer mantissa digits to keep (default `8`).
#'
#' @return Character vector the same length as `x`. `NA_character_` for
#'   non-finite entries. Round-trips with [log10c()]:
#'   `log10c(pow10c(x)) ≈ x` and `pow10c(log10c(s)) ≈ s` (up to mantissa rounding).
#'
#' @seealso [log10c()] for the reverse direction.
#' @export
pow10c <- function(x, digits = 8L) {
  x <- as_num(x)
  out <- rep(NA_character_, length(x))
  ok <- is.finite(x) & !is.na(x)
  if (!any(ok)) return(out)

  d <- as_int(digits)[1]
  if (!is.finite(d) || is.na(d) || d < 1L) d <- 8L

  expo <- floor(x[ok])
  mant <- 10^(x[ok] - expo)

  # Carry mantissa into next decade if rounding pushes it over 10.
  mant_rounded <- round(mant, d)
  carry <- mant_rounded >= 10
  if (any(carry)) {
    mant[carry] <- mant[carry] / 10
    expo[carry] <- expo[carry] + 1L
  }

  mtxt <- formatC(mant, format = "f", digits = d)
  mtxt <- sub("\\.?0+$", "", mtxt)
  out[ok] <- paste0(mtxt, "e", ifelse(expo >= 0, "+", ""), as_int(expo))
  out
}

#' Two-tailed (or one-tailed) p-value from z or (beta, SE)
#'
#' Computes p-values on the log scale internally. For mid-range p-values the
#' return is numeric; for entries below `tiny.threshold` — where `exp(logp)`
#' would underflow to zero — the return is automatically switched to a
#' character scientific-notation vector via [pow10c()] so precision is
#' preserved (e.g. `"1.23e-500"`).
#'
#' This inherits the precision-without-underflow philosophy of [log10c()] /
#' [pow10c()]: the function never silently rounds a tiny p to `0`.
#'
#' @param z Numeric vector of z-scores. If supplied, `beta`/`se`/`varbeta` are ignored.
#' @param beta,se,varbeta Alternative inputs passed to [zabs()] when `z` is `NULL`.
#' @param two.tailed Logical. Two-tailed p-value if `TRUE` (default).
#' @param log.p Logical. Return natural-log p-values if `TRUE` (matches
#'   `stats::pnorm(..., log.p = TRUE)` convention).
#' @param tiny.threshold Numeric `(0, 1)`. P-values below this are emitted
#'   as character strings (default `1e-317`, the IEEE-754 normal-range edge).
#' @param tiny.as.character Logical. If `FALSE`, force a numeric return even
#'   when tiny entries exist (those underflow to `0`).
#' @param digits Integer. Mantissa digits used in the character
#'   representation (default `8`).
#'
#' @return Numeric vector when no entry falls below `tiny.threshold` or
#'   `tiny.as.character = FALSE`. Otherwise a character vector with
#'   scientific notation throughout (mid-range entries are also rendered
#'   as character so the result is type-consistent).
#'
#' @seealso [pow10c()] which underlies the character formatting, and
#'   [format_pvalue()] for HTML / Markdown labels.
#' @export
pvalue <- function(z = NULL,
                   beta = NULL,
                   se = NULL,
                   varbeta = NULL,
                   two.tailed = TRUE,
                   log.p = FALSE,
                   tiny.threshold = 1e-317,
                   tiny.as.character = TRUE,
                   digits = 8L) {
  if (!is.null(z)) {
    absz0 <- abs(as_num(z))
  } else {
    absz0 <- zabs(beta = beta, se = se, varbeta = varbeta)
  }
  absz0[!is.finite(absz0)] <- NA_real_

  two.tailed <- isTRUE(two.tailed)
  logp <- if (two.tailed) {
    log(2) + stats::pnorm(-absz0, log.p = TRUE)
  } else {
    stats::pnorm(-absz0, log.p = TRUE)
  }
  logp[!is.finite(logp)] <- -Inf
  if (isTRUE(log.p)) return(logp)

  tiny.threshold <- as_num(tiny.threshold)[1]
  if (!is.finite(tiny.threshold) || is.na(tiny.threshold) || tiny.threshold <= 0 || tiny.threshold >= 1) {
    tiny.threshold <- 1e-317
  }
  p_num <- exp(logp)
  tiny <- is.finite(logp) & !is.na(logp) & (logp < log(tiny.threshold))

  if (!isTRUE(tiny.as.character) || !any(tiny, na.rm = TRUE)) {
    p_num[!is.finite(p_num)] <- NA_real_
    return(p_num)
  }

  # Convert the whole vector to character via pow10c so output is type-stable.
  out <- pow10c(logp / log(10), digits = digits)
  out
}

#' Format p-values for `ggtext::geom_richtext()` / Markdown labels
#'
#' Converts numeric, character, or `log10(p)` inputs into a typeset
#' scientific-notation string like `"3.2 &times; 10<sup>-8</sup>"` (HTML,
#' default — suitable for [ggtext::geom_richtext()] /
#' [ggtext::element_markdown()]) or `"3.2 x 10^-8"` (plain text).
#'
#' Tiny p-values that would underflow under `as.numeric()` are routed through
#' [log10c()] when the input is character, or supplied directly via
#' `log10p = TRUE`. Either way precision is preserved.
#'
#' @param x P-values (numeric or character) or, when `log10p = TRUE`,
#'   base-10-log p-values (i.e. `log10(p)`, **not** natural log).
#' @param log10p Logical. If `TRUE`, `x` is treated as `log10(p)`. Default `FALSE`.
#' @param digits Mantissa digits to keep (default `2`).
#' @param plain.cutoff Numeric `(0, 1]` or `NULL`. P-values at or above this
#'   value are rendered as plain decimals (e.g. `"0.05"`) instead of
#'   scientific notation. `NULL` (default) means "always use scientific".
#' @param html Logical. If `TRUE` (default), emit HTML-rich output using
#'   `<sup>` and the `times` glyph -- ready for
#'   [ggtext::geom_richtext()]. If `FALSE`, emit a plain ASCII form,
#'   `"<mantissa> x 10^<exp>"` (e.g. `"3.20 x 10^-8"`).
#' @param times Character. The multiplication glyph in HTML mode. Defaults
#'   to the HTML entity `"&times;"`. Ignored when `html = FALSE`.
#'
#' @return Character vector of the same length as `x`.
#'
#' @examples
#' \dontrun{
#'   library(ggplot2); library(ggtext)
#'   ggplot(df, aes(x, y, label = format_pvalue(P))) +
#'     geom_richtext()
#'   # log10(p) input (e.g. from -negLog10P in a Manhattan plot)
#'   format_pvalue(-7.3, log10p = TRUE)
#'   # plain-text fallback (e.g. for base graphics or copy-paste)
#'   format_pvalue(5e-8, html = FALSE)
#' }
#' @export
format_pvalue <- function(x,
                          log10p = FALSE,
                          digits = 2L,
                          plain.cutoff = NULL,
                          html = TRUE,
                          times = "&times;") {
  if (is.null(x) || length(x) == 0L) return(character())

  if (isTRUE(log10p)) {
    log10v_all <- as_num(x)
  } else {
    if (is.character(x)) {
      log10v_all <- suppressWarnings(log10c(x))
    } else {
      log10v_all <- suppressWarnings(log10(as_num(x)))
    }
  }

  d <- as_int(digits)[1]
  if (!is.finite(d) || is.na(d) || d < 1L) d <- 2L

  pc_use <- NULL
  if (!is.null(plain.cutoff) && length(plain.cutoff)) {
    pc_val <- suppressWarnings(as.numeric(plain.cutoff)[1])
    if (is.finite(pc_val) && !is.na(pc_val) && pc_val > 0 && pc_val <= 1) {
      pc_use <- pc_val
    }
  }

  html_use <- isTRUE(html)
  times <- as.character(times)[1]
  if (is.na(times) || !nzchar(times)) times <- "&times;"

  out <- rep(NA_character_, length(log10v_all))
  ok <- is.finite(log10v_all) & !is.na(log10v_all)
  if (!any(ok)) return(out)

  log10v <- log10v_all[ok]
  expo <- floor(log10v)
  mant <- 10^(log10v - expo)

  mant_rounded <- round(mant, d)
  carry <- mant_rounded >= 10
  if (any(carry)) {
    mant[carry] <- mant[carry] / 10
    expo[carry] <- expo[carry] + 1L
  }

  mtxt <- formatC(mant, format = "f", digits = d)
  mtxt <- sub("\\.?0+$", "", mtxt)

  res <- if (html_use) {
    sprintf("%s %s 10<sup>%d</sup>", mtxt, times, as_int(expo))
  } else {
    sprintf("%s x 10^%d", mtxt, as_int(expo))
  }

  if (!is.null(pc_use)) {
    p_val <- mant * 10^expo
    plain_idx <- p_val >= pc_use
    if (any(plain_idx)) {
      plain_txt <- formatC(p_val[plain_idx], format = "fg", digits = d + 1L, flag = "#")
      plain_txt <- sub("\\.?0+$", "", plain_txt)
      plain_txt <- sub("^[[:space:]]+", "", plain_txt)
      res[plain_idx] <- plain_txt
    }
  }

  out[ok] <- res
  out
}

#' Credible set labels from posterior inclusion probabilities
#'
#' Greedy cumulative-mass selection: variants are added in decreasing PIP
#' order until the running sum reaches each requested confidence threshold.
#' The returned vector labels each input with the smallest confidence level it
#' belongs to (0 means not in any credible set).
#'
#' @param pip Numeric vector of posterior inclusion probabilities.
#' @param confidence Numeric vector of confidence levels in `(0, 1]`.
#'
#' @return Numeric vector the same length as `pip`. The attribute `"cutoff"`
#'   carries the PIP cutoff at each level.
#' @export
credibleset <- function(pip, confidence = c(0.95, 0.99)) {
  pip_num <- as_num(pip)
  name_ori <- names(pip)
  n <- length(pip_num)
  if (!n) return(numeric())

  conf <- as_num(confidence)
  conf <- conf[is.finite(conf) & !is.na(conf) & conf > 0 & conf <= 1]
  conf <- sort(unique(conf))
  if (!length(conf)) stop("confidence must contain numeric values in (0, 1].", call. = FALSE)

  valid <- is.finite(pip_num) & !is.na(pip_num) & pip_num >= 0
  out <- rep(0, n)
  if (all(!valid)) {
    names(out) <- name_ori
    return(out)
  }

  pip_work <- pip_num[valid]
  ord <- order(pip_work, decreasing = TRUE, na.last = TRUE)
  pip_sorted <- pip_work[ord]
  total <- sum(pip_sorted)
  if (!is.finite(total) || total <= 0) {
    names(out) <- name_ori
    return(out)
  }

  csum <- cumsum(pip_sorted)
  bounds <- total * conf
  cutoffs <- vapply(bounds, function(b) {
    i <- which(csum >= b)[1]
    if (is.na(i) || !is.finite(i)) min(pip_sorted) else pip_sorted[i]
  }, numeric(1))

  for (i in seq_along(conf)) {
    j <- length(conf) - i + 1L
    out[valid & pip_num >= cutoffs[j]] <- conf[j]
  }
  names(out) <- name_ori
  names(cutoffs) <- as.character(conf)
  attr(out, "cutoff") <- cutoffs
  out
}

#' Quantile bucket labels
#'
#' Assigns each value of `x` to one of `n` equal-probability bins, returning a
#' factor with levels `"1st", "2nd", "3rd", "4th", ...`.
#'
#' @param x Numeric vector.
#' @param n Integer, number of quantile bins (default 4 = quartiles).
#'
#' @return A factor of length `length(x)`.
#' @export
quantvec <- function(x, n = 4L) {
  x_num <- as_num(x)
  n <- as_int(n)[1]
  if (!is.finite(n) || is.na(n) || n < 1L) stop("n must be a positive integer.", call. = FALSE)

  qnames <- c("1st", "2nd", "3rd", paste0(4:(4 + max(0L, n - 4L)), "th"))[seq_len(n)]
  out_chr <- rep(NA_character_, length(x_num))

  ok <- is.finite(x_num) & !is.na(x_num)
  if (!any(ok)) {
    out <- factor(out_chr, levels = qnames)
    names(out) <- names(x)
    return(out)
  }

  qv <- stats::quantile(x_num[ok], probs = seq(0, 1, length.out = n + 1L), na.rm = TRUE, names = FALSE)
  for (i in seq_len(n - 1L)) {
    idx <- ok & x_num >= qv[i] & x_num < qv[i + 1L]
    out_chr[idx] <- qnames[i]
  }
  idx_last <- ok & x_num >= qv[n]
  out_chr[idx_last] <- qnames[n]

  out <- factor(out_chr, levels = qnames)
  names(out) <- names(x)
  out
}

#' Rank-based inverse normal transform
#'
#' Maps `x` to standard-normal quantiles via `qnorm((rank - c) / (n - 2c + 1))`.
#'
#' @param x Numeric vector to transform.
#' @param c Blom-style constant (default `3/8`).
#' @param stochastic Logical. If `TRUE`, breaks ties by a deterministic random
#'   shuffle (seeded internally for reproducibility).
#'
#' @return Numeric vector of transformed values, same length as `x`.
#' @export
invnorm <- function(x, c = 3.0 / 8, stochastic = FALSE) {
  names_ori <- names(x)
  x_num <- as_num(x)
  n <- length(x_num)
  out <- rep(NA_real_, n)
  if (!n) {
    names(out) <- names_ori
    return(out)
  }

  c0 <- suppressWarnings(as.numeric(c))[1]
  if (!is.finite(c0) || is.na(c0)) c0 <- 3.0 / 8

  ok <- !is.na(x_num) & is.finite(x_num)
  if (!any(ok)) {
    names(out) <- names_ori
    return(out)
  }

  xv <- x_num[ok]
  if (isTRUE(stochastic)) {
    set.seed(3L)
    ord <- sample.int(length(xv))
    xv_shuf <- xv[ord]
    rk_shuf <- rank(xv_shuf, ties.method = "first")
    rk <- numeric(length(xv))
    rk[ord] <- rk_shuf
  } else {
    rk <- rank(xv, ties.method = "average")
  }

  denom <- length(xv) - 2 * c0 + 1
  if (!is.finite(denom) || denom <= 0) {
    names(out) <- names_ori
    return(out)
  }
  p <- (rk - c0) / denom
  p <- pmax(.Machine$double.eps, pmin(1 - .Machine$double.eps, p))
  out[ok] <- stats::qnorm(p)
  names(out) <- names_ori
  out
}
