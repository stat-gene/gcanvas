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

#' Two-tailed (or one-tailed) p-value from z or (beta, SE)
#'
#' Computes p-values on the log scale internally and switches to a string
#' representation in scientific notation for entries that would underflow to
#' zero under `exp()`.
#'
#' @param z Numeric vector of z-scores. If supplied, `beta`/`se`/`varbeta` are ignored.
#' @param beta,se,varbeta Alternative inputs passed to [zabs()] when `z` is `NULL`.
#' @param two.tailed Logical. Two-tailed p-value if `TRUE` (default).
#' @param log.p Logical. Return natural-log p-values if `TRUE`.
#' @param tiny.threshold Numeric. P-values below this are emitted as character
#'   strings (preserving precision instead of underflowing to 0).
#' @param tiny.as.character Logical. Disable the character fallback if `FALSE`.
#' @param digits Integer. Mantissa digits used in the character representation.
#'
#' @return Numeric (or character) vector of p-values, same length as the input.
#' @export
pvalue <- function(z = NULL,
                   beta = NULL,
                   se = NULL,
                   varbeta = NULL,
                   two.tailed = TRUE,
                   log.p = FALSE,
                   tiny.threshold = 1e-317,
                   tiny.as.character = TRUE,
                   digits = 6L) {
  .pvalue_logp_to_sci <- function(lp, digits = 6L) {
    lp <- as_num(lp)
    out <- rep(NA_character_, length(lp))
    ok <- is.finite(lp) & !is.na(lp)
    if (!any(ok)) return(out)
    log10v <- lp[ok] / log(10)
    expo <- floor(log10v)
    mant <- exp(lp[ok] - expo * log(10))
    shift <- mant >= 10
    if (any(shift)) {
      mant[shift] <- mant[shift] / 10
      expo[shift] <- expo[shift] + 1
    }
    d <- as_int(digits)[1]
    if (!is.finite(d) || is.na(d) || d < 1L) d <- 6L
    mtxt <- formatC(mant, format = "f", digits = d)
    mtxt <- sub("\\.?0+$", "", mtxt)
    out[ok] <- paste0(mtxt, "e", ifelse(expo >= 0, "+", ""), as_int(expo))
    out
  }

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

  out <- rep(NA_character_, length(logp))
  ok <- is.finite(logp) & !is.na(logp)
  out[ok] <- format(signif(p_num[ok], 16), scientific = TRUE, trim = TRUE)
  out[tiny] <- .pvalue_logp_to_sci(logp[tiny], digits = digits)
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
