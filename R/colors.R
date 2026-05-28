# Color palette generators and inspector for the gcanvas plotting layer.
# Exposes `get.colors()` (random + curated + palette modes) and `show.pal()`
# (visualize a palette by name or vector). All other helpers are internal.

.gcanvas_ld_palette <- function(n = 100L) {
  n <- as_int(n)[1]
  if (!is.finite(n) || is.na(n) || n < 1L) return(character())
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    stop("mode='ld' requires package 'RColorBrewer'.", call. = FALSE)
  }
  cols <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")))(n)
  cols[1] <- "#333399"
  cols
}

#' Generate a color palette
#'
#' Flexible palette generator supporting random spread, perceptual-distance
#' selection, curated palettes (Brewer, Carto, viridis, named bases), tonal
#' variation, and pinning of must-include colors.
#'
#' @param mode Palette mode (e.g. `"random"`, `"ld"`, a palette name, or one of
#'   the curated bases). See body for the full list.
#' @param n Number of colors to return. When `NULL`, a sensible default for
#'   the chosen mode is used.
#' @param discrete Logical. If `TRUE`, snap to discrete palette entries instead
#'   of interpolating.
#' @param far Logical. If `TRUE`, choose colors that maximize perceptual
#'   distance in Lab space.
#' @param random Logical. If `FALSE`, force deterministic palette ordering
#'   (`seed` ignored, no tonal shifts).
#' @param silent Logical. Suppress progress notes.
#' @param include.hue,include.color Optional hue / color anchors to pin.
#' @param tone,ash Tonal-shift / desaturation modifiers.
#' @param hue.shift,satur.mult,light.offset HSL adjustment knobs applied to
#'   the chosen palette.
#' @param alpha Output alpha (0-1).
#' @param seed Integer / `"random"` seed (only honored in random / tonal modes).
#' @param plot Logical. Also display the palette swatches.
#' @param label One of `"top"`, `"axis"`, `"none"` — label placement on the swatch.
#' @param title Optional title for the swatch.
#' @param label.size Relative label text size.
#'
#' @return Character vector of hex color codes of length `n`.
#' @export
get.colors <- function(mode = "random",
                       n = NULL,
                       discrete = FALSE,
                       far = TRUE,
                       random = TRUE,
                       silent = FALSE,
                       include.hue = NULL,
                       include.color = NULL,
                       tone = FALSE,
                       ash = FALSE,
                       hue.shift = 0,
                       satur.mult = 1,
                       light.offset = 0,
                       alpha = 1,
                       seed = NULL,
                       plot = TRUE,
                       label = c("top", "axis", "none"),
                       title = NULL,
                       label.size = 1) {
  n_missing <- is.null(n) || length(n) == 0L || !is.finite(suppressWarnings(as.numeric(n))[1])
  mode_num_in <- suppressWarnings(as.numeric(mode))[1]
  # If the first positional argument is numeric (e.g., get.colors(10)),
  # treat it as n and use random mode.
  if (length(mode) == 1L && n_missing && is.finite(mode_num_in) && !is.na(mode_num_in)) {
    n <- mode_num_in
    mode <- "random"
  }

  custom_mode_base_raw <- NULL
  mode_chr_all <- as.character(mode)
  mode_chr_all <- mode_chr_all[!is.na(mode_chr_all) & nzchar(mode_chr_all)]
  is_color_vec <- FALSE
  if (length(mode_chr_all) > 0L) {
    chk <- vapply(mode_chr_all, function(z) {
      !is.null(tryCatch(grDevices::col2rgb(z), error = function(e) NULL))
    }, logical(1))
    is_color_vec <- all(chk)
  }
  if (isTRUE(is_color_vec)) {
    custom_mode_base_raw <- mode_chr_all
    mode_raw <- "__custom__"
  } else {
    mode_raw <- as.character(mode)[1]
    mode_num <- suppressWarnings(as.numeric(mode_raw))
    if (is.finite(mode_num) && (is.null(n) || length(n) == 0L || !is.finite(suppressWarnings(as.numeric(n))[1]))) {
      n <- mode_num
      mode_raw <- "random"
    }
  }
  mode_raw <- tolower(trimws(mode_raw))
  if (mode_raw %in% c("continuous", "discrete")) {
    discrete <- identical(mode_raw, "discrete")
    mode_raw <- "random"
  }
  if (!nzchar(mode_raw) || is.na(mode_raw)) mode_raw <- "random"
  mode <- mode_raw
  discrete <- isTRUE(discrete)
  far <- isTRUE(far)
  random <- isTRUE(random)
  .get_colors_carto_palette_canonical <- function(mode_name) {
    if (!requireNamespace("rcartocolor", quietly = TRUE)) return(NULL)
    m <- tolower(as.character(mode_name)[1])
    if (is.na(m) || !nzchar(m)) return(NULL)
    nm <- as.character(rcartocolor::cartocolors$Name)
    hit <- nm[tolower(nm) == m]
    if (!length(hit)) return(NULL)
    as.character(hit[1])
  }
  .get_colors_carto_palette_meta <- function(pal_name) {
    if (is.null(pal_name) || !nzchar(as.character(pal_name)[1])) return(NULL)
    if (!requireNamespace("rcartocolor", quietly = TRUE)) return(NULL)
    meta <- rcartocolor::metacartocolors
    i <- which(as.character(meta$Name) == as.character(pal_name)[1])[1]
    if (!is.finite(i) || is.na(i)) return(NULL)
    list(
      min_n = as_int(meta$Min_n[i]),
      max_n = as_int(meta$Max_n[i]),
      type = tolower(as.character(meta$Type[i]))
    )
  }
  .get_colors_default_n_from_mode <- function(mode_name) {
    m <- tolower(as.character(mode_name)[1])
    if (identical(m, "__custom__")) return(length(custom_mode_base_raw))
    if (!nzchar(m) || is.na(m) || identical(m, "random")) return(10L)
    if (identical(m, "rainbow")) return(10L)
    if (identical(m, "darkrainbow")) return(10L)
    if (identical(m, "lightrainbow")) return(10L)
    if (identical(m, "ld")) return(100L)
    if (identical(m, "ggplot")) return(6L)
    carto_name <- .get_colors_carto_palette_canonical(m)
    if (!is.null(carto_name)) {
      meta <- .get_colors_carto_palette_meta(carto_name)
      if (!is.null(meta) && is.finite(meta$max_n) && !is.na(meta$max_n)) return(as_int(meta$max_n))
      return(10L)
    }
    if (m %in% c("viridis", "magma", "plasma", "inferno", "cividis", "rocket", "mako", "turbo")) return(256L)
    if (requireNamespace("RColorBrewer", quietly = TRUE)) {
      info <- RColorBrewer::brewer.pal.info
      rn <- rownames(info)
      idx <- match(m, tolower(rn))
      if (!is.na(idx)) return(as_int(info[idx, "maxcolors"]))
    }
    10L
  }
  if (is.null(n) || length(n) == 0L || !is.finite(suppressWarnings(as.numeric(n))[1])) {
    n <- .get_colors_default_n_from_mode(mode)
  }
  palette_base_n <- max(as_int(n), as_int(.get_colors_default_n_from_mode(mode)))

  if (is.logical(tone) && length(tone) > 0L && !isTRUE(tone)) {
    tone_level <- NA_real_
  } else if (is.logical(tone) && length(tone) > 0L && isTRUE(tone)) {
    tone_level <- 0.5
  } else if (is.numeric(tone) && length(tone) > 0L) {
    tone_level <- suppressWarnings(as.numeric(tone))[1]
    if (!is.finite(tone_level) || is.na(tone_level)) tone_level <- 0
    tone_level <- max(0, min(1, tone_level))
  } else {
    tone_chr <- tolower(as.character(tone)[1])
    if (!nzchar(tone_chr) || is.na(tone_chr)) tone_chr <- "soft"
    if (!(tone_chr %in% c("soft", "vivid"))) stop("tone must be 'soft', 'vivid', or numeric in [0, 1].", call. = FALSE)
    tone_level <- if (identical(tone_chr, "vivid")) 1 else 0
  }
  if (is.logical(ash) && length(ash) > 0L && !isTRUE(ash)) {
    ash_level <- NA_real_
  } else if (is.logical(ash) && length(ash) > 0L && isTRUE(ash)) {
    ash_level <- 0.6
  } else if (is.numeric(ash) && length(ash) > 0L) {
    ash_level <- suppressWarnings(as.numeric(ash))[1]
    if (!is.finite(ash_level) || is.na(ash_level)) {
      ash_level <- NA_real_
    } else {
      ash_level <- max(0, min(1, ash_level))
    }
  } else {
    ash_level <- NA_real_
  }
  tone_enabled <- is.finite(tone_level) && !is.na(tone_level)
  ash_enabled <- is.finite(ash_level) && !is.na(ash_level)
  transform_enabled <- isTRUE(tone_enabled || ash_enabled)

  label <- match.arg(label)
  n <- as_int(n)[1]
  if (!is.finite(n) || is.na(n) || n < 1L) return(character())
  alpha <- suppressWarnings(as.numeric(alpha))[1]
  if (!is.finite(alpha) || is.na(alpha)) alpha <- 1
  alpha <- max(0, min(1, alpha))
  hue.shift <- suppressWarnings(as.numeric(hue.shift))[1]
  if (!is.finite(hue.shift) || is.na(hue.shift)) hue.shift <- 0
  satur.mult <- suppressWarnings(as.numeric(satur.mult))[1]
  if (!is.finite(satur.mult) || is.na(satur.mult) || satur.mult < 0) satur.mult <- 1
  light.offset <- suppressWarnings(as.numeric(light.offset))[1]
  if (!is.finite(light.offset) || is.na(light.offset)) light.offset <- 0
  label.size <- suppressWarnings(as.numeric(label.size))[1]
  if (!is.finite(label.size) || is.na(label.size) || label.size <= 0) label.size <- 1

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- NULL
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)

  .get_colors_canonical_seed <- function(x) {
    m <- 2147483647
    x <- suppressWarnings(as.numeric(x))[1]
    if (!is.finite(x) || is.na(x)) return(NA_integer_)
    s <- suppressWarnings(x %% m)
    if (!is.finite(s) || is.na(s)) return(NA_integer_)
    if (s <= 0) s <- s + m
    s <- as.integer(floor(s))
    if (!is.finite(s) || is.na(s) || s <= 0L) s <- 1L
    s
  }
  seed0 <- .get_colors_canonical_seed(seed)
  if (!is.finite(seed0) || is.na(seed0)) seed0 <- .get_colors_canonical_seed((as.numeric(Sys.time()) * 1e6 + Sys.getpid()))
  seed_use <- as_int(seed0)
  deterministic_palette_mode <- (!identical(mode, "random")) && !isTRUE(transform_enabled)
  random <- if (isTRUE(deterministic_palette_mode)) FALSE else isTRUE(random)
  use_rng <- identical(mode, "random") || isTRUE(random)
  if (isTRUE(use_rng)) {
    set.seed(as_int(seed_use))
    .gcanvas_note("gcanvas::get.colors", sprintf("seed=%d", as_int(seed_use)), silent = silent)
  } else {
    if (isTRUE(deterministic_palette_mode)) {
      .gcanvas_note("gcanvas::get.colors", "deterministic palette mode (tone=FALSE, ash=FALSE; seed ignored)", silent = silent)
    } else {
      .gcanvas_note("gcanvas::get.colors", "random=FALSE (deterministic palette mode)", silent = silent)
    }
  }

  .get_colors_fix_hex <- function(x) {
    x <- as.character(x)
    x[is.na(x) | !nzchar(x)] <- "#808080"
    x
  }
  .get_colors_as_hex <- function(x) {
    x <- as.character(x)
    x <- x[!is.na(x) & nzchar(x)]
    if (!length(x)) return(character())
    out <- vapply(x, function(z) {
      rgb <- tryCatch(grDevices::col2rgb(z), error = function(e) NULL)
      if (is.null(rgb)) return(NA_character_)
      grDevices::rgb(rgb[1, 1], rgb[2, 1], rgb[3, 1], maxColorValue = 255, alpha = as_int(alpha * 255))
    }, character(1))
    unique(out[!is.na(out) & nzchar(out)])
  }
  custom_mode_base <- if (!is.null(custom_mode_base_raw)) .get_colors_as_hex(custom_mode_base_raw) else character()
  if (identical(mode, "__custom__") && length(custom_mode_base) == 0L) {
    stop("mode vector is provided but no valid colors were found.", call. = FALSE)
  }
  .get_colors_as_hue <- function(x) {
    if (is.null(x) || !length(x)) return(numeric())
    if (is.numeric(x)) return(suppressWarnings(as.numeric(x))[is.finite(suppressWarnings(as.numeric(x)))] %% 360)
    xx <- tolower(trimws(as.character(x)))
    map <- c(red = 15, orange = 35, yellow = 60, lime = 95, green = 120, teal = 165, cyan = 190, sky = 210, blue = 240, indigo = 265, purple = 285, magenta = 315, pink = 335)
    h <- suppressWarnings(as.numeric(xx)); out <- h[is.finite(h)]
    if (any(!is.finite(h))) {key <- xx[!is.finite(h)]; out <- c(out, unname(map[key[key %in% names(map)]]))}
    out[is.finite(out)] %% 360
  }
  .get_colors_text_col_for_bg <- function(hex) {
    rgb <- grDevices::col2rgb(hex) / 255
    y <- 0.2126 * rgb[1, ] + 0.7152 * rgb[2, ] + 0.0722 * rgb[3, ]
    ifelse(y < 0.55, "white", "black")
  }
  .get_colors_lerp <- function(a, b, t) a + (b - a) * t
  .get_colors_circ_dist <- function(x, y) {d <- abs(x - y); pmin(d, 360 - d)}
  .get_colors_sample_hues_far <- function(k, base_h = numeric()) {
    if (k <= 0L) return(numeric())
    base_h <- unique(suppressWarnings(as.numeric(base_h))[is.finite(suppressWarnings(as.numeric(base_h)))] %% 360)
    cand <- 0:359; chosen <- base_h; out <- numeric(k)
    for (i in seq_len(k)) {
      if (!length(chosen)) pick <- sample(cand, 1L) else {
        dmin <- vapply(cand, function(h) min(.get_colors_circ_dist(h, chosen)), numeric(1))
        pick <- cand[which.max(dmin + stats::runif(length(dmin), 0, 2.5))]
      }
      out[i] <- pick; chosen <- c(chosen, pick)
    }
    out
  }
  .get_colors_pick_far_hex <- function(pool_hex, k, randomize = TRUE) {
    pool_hex <- unique(.get_colors_fix_hex(pool_hex))
    if (k <= 0L || !length(pool_hex)) return(character())
    if (k >= length(pool_hex)) return(pool_hex[seq_len(min(k, length(pool_hex)))])
    rgb <- t(grDevices::col2rgb(pool_hex))
    lab <- grDevices::convertColor(rgb, from = "sRGB", to = "Lab", scale.in = 255)
    picked <- if (isTRUE(randomize)) sample.int(nrow(lab), 1L) else 1L
    while (length(picked) < k) {
      rem <- setdiff(seq_len(nrow(lab)), picked)
      dmin <- vapply(rem, function(i) min(sqrt(rowSums((lab[picked, , drop = FALSE] - matrix(lab[i, ], nrow = length(picked), ncol = 3, byrow = TRUE))^2))), numeric(1))
      picked <- c(picked, rem[which.max(dmin)])
    }
    pool_hex[picked]
  }
  .get_colors_pick_stride_hex <- function(pool_hex, k) {
    pool_hex <- .get_colors_fix_hex(pool_hex)
    m <- length(pool_hex)
    if (k <= 0L || m == 0L) return(character())
    if (k == 1L) return(pool_hex[1])
    idx <- unique(as_int(round(seq(1, m, length.out = k))))
    if (length(idx) < k) {
      add <- setdiff(seq_len(m), idx)
      idx <- c(idx, head(add, k - length(idx)))
    }
    pool_hex[idx[seq_len(min(k, length(idx)))]]
  }
  .get_colors_apply_hsv_tone <- function(cols, sat_rng, val_rng, hue_shift = 0, satur_mult = 1, light_offset = 0, randomize = TRUE) {
    rgb_mat <- grDevices::col2rgb(cols)
    hsv_mat <- grDevices::rgb2hsv(rgb_mat[1, ], rgb_mat[2, ], rgb_mat[3, ])
    if (isTRUE(randomize)) {
      s_mul <- stats::runif(length(cols), sat_rng[1], sat_rng[2])
      v_add <- stats::runif(length(cols), val_rng[1], val_rng[2])
    } else {
      t <- if (length(cols) <= 1L) 0.5 else seq(0, 1, length.out = length(cols))
      s_mul <- sat_rng[1] + (sat_rng[2] - sat_rng[1]) * t
      v_add <- val_rng[1] + (val_rng[2] - val_rng[1]) * t
    }
    s_new <- pmin(1, pmax(0, hsv_mat["s", ] * s_mul * satur_mult))
    v_new <- pmin(1, pmax(0, hsv_mat["v", ] + v_add + light_offset))
    h_new <- (hsv_mat["h", ] + (hue_shift / 360)) %% 1
    .get_colors_fix_hex(grDevices::hsv(h_new, s_new, v_new, alpha = alpha))
  }

  include_cols <- .get_colors_as_hex(include.color)
  include_h <- .get_colors_as_hue(include.hue)
  tv <- if (!isTRUE(tone_enabled)) NA_real_ else tone_level^1.35
  tv_rand <- if (is.na(tv)) 0.5 else tv
  l_hi <- c(.get_colors_lerp(88, 70, tv_rand), .get_colors_lerp(97, 82, tv_rand))
  l_mid <- c(.get_colors_lerp(72, 50, tv_rand), .get_colors_lerp(84, 64, tv_rand))
  l_lo <- c(.get_colors_lerp(56, 32, tv_rand), .get_colors_lerp(70, 46, tv_rand))
  c_lo <- c(.get_colors_lerp(6, 24, tv_rand), .get_colors_lerp(18, 42, tv_rand))
  c_mid <- c(.get_colors_lerp(14, 48, tv_rand), .get_colors_lerp(30, 72, tv_rand))
  c_hi <- c(.get_colors_lerp(28, 72, tv_rand), .get_colors_lerp(48, 96, tv_rand))
  disc_c_min <- .get_colors_lerp(18, 66, tv_rand); disc_c_max <- .get_colors_lerp(45, 98, tv_rand)
  disc_l_min <- .get_colors_lerp(64, 34, tv_rand); disc_l_max <- .get_colors_lerp(90, 74, tv_rand)
  sat_rng <- if (is.na(tv)) c(1, 1) else c(.get_colors_lerp(0.55, 1.10, tv), .get_colors_lerp(0.90, 1.45, tv))
  val_rng <- if (is.na(tv)) c(0, 0) else c(.get_colors_lerp(0.10, -0.08, tv), .get_colors_lerp(0.24, 0.04, tv))
  if (isTRUE(ash_enabled)) {
    sat_rng <- sat_rng * (1 - 0.65 * ash_level)
    sat_rng <- pmax(0, sat_rng)
    val_rng <- val_rng + c(-0.18 * ash_level, -0.04 * ash_level)
  }
  if (!isTRUE(transform_enabled)) {
    sat_rng <- c(1, 1)
    val_rng <- c(0, 0)
  }
  hue_shift_eff <- if (isTRUE(transform_enabled)) hue.shift + if (isTRUE(ash_enabled)) (-8 * ash_level) else 0 else 0
  satur_mult_eff <- if (isTRUE(transform_enabled)) satur.mult else 1
  light_offset_eff <- if (isTRUE(transform_enabled)) light.offset else 0
  rb_base <- c("#FF595E", "#FF924C", "#FFCA3A", "#C5CA30", "#8AC926", "#36949D", "#1982C4", "#4267AC", "#565AA0", "#6A4C93")
  drb_base <- c("#CD434D", "#D47240", "#D2A646", "#C5C145", "#8BB631", "#39A6AA", "#2F91C3", "#4B77B7", "#4B539F", "#6C49A4")
  lrb_base <- c("#FA7B7F", "#FAA670", "#F9D060", "#CBCF50", "#9BD441", "#4EA5AD", "#3394D1", "#5A7BBB", "#6E72AD", "#7E62A4")
  viridis_opts <- c(
    viridis = "D", magma = "A", inferno = "B", plasma = "C",
    cividis = "E", rocket = "F", mako = "G", turbo = "H"
  )
  .get_colors_palette_from_mode <- function(mode_name, n_out) {
    n_out <- as_int(n_out)
    if (!is.finite(n_out) || n_out < 1L) return(character())
    m <- tolower(as.character(mode_name)[1])
    if (identical(m, "__custom__")) {
      if (!length(custom_mode_base)) return(character())
      if (n_out <= length(custom_mode_base)) return(custom_mode_base[seq_len(n_out)])
      return(.get_colors_fix_hex(grDevices::colorRampPalette(custom_mode_base, space = "Lab")(n_out)))
    }
    if (identical(m, "rainbow")) {
      if (n_out <= length(rb_base)) return(rb_base[seq_len(n_out)])
      return(.get_colors_fix_hex(grDevices::colorRampPalette(rb_base, space = "Lab")(n_out)))
    }
    if (identical(m, "darkrainbow")) {
      if (n_out <= length(drb_base)) return(drb_base[seq_len(n_out)])
      return(.get_colors_fix_hex(grDevices::colorRampPalette(drb_base, space = "Lab")(n_out)))
    }
    if (identical(m, "lightrainbow")) {
      if (n_out <= length(lrb_base)) return(lrb_base[seq_len(n_out)])
      return(.get_colors_fix_hex(grDevices::colorRampPalette(lrb_base, space = "Lab")(n_out)))
    }
    if (identical(m, "ld")) {
      return(.gcanvas_ld_palette(n_out))
    }
    if (identical(m, "ggplot")) {
      # ggplot-like qualitative hue sequence:
      # hues <- seq(15, 375, length = n + 1); hcl(h = hues, l = 65, c = 100)[1:n]
      hues <- seq(15, 375, length.out = n_out + 1L)
      return(.get_colors_fix_hex(grDevices::hcl(h = hues[seq_len(n_out)], c = 100, l = 65, alpha = alpha)))
    }
    carto_name <- .get_colors_carto_palette_canonical(m)
    if (!is.null(carto_name)) {
      meta <- .get_colors_carto_palette_meta(carto_name)
      if (is.null(meta)) {
        return(.get_colors_fix_hex(rcartocolor::carto_pal(n = n_out, name = carto_name)))
      }
      min_n <- if (is.finite(meta$min_n) && !is.na(meta$min_n)) as_int(meta$min_n) else 2L
      max_n <- if (is.finite(meta$max_n) && !is.na(meta$max_n)) as_int(meta$max_n) else max(n_out, 2L)
      if (n_out <= 1L) {
        base_n <- max(min_n, min(max_n, 3L))
        base_cols <- .get_colors_fix_hex(rcartocolor::carto_pal(n = base_n, name = carto_name))
        return(base_cols[as_int(ceiling(length(base_cols) / 2))])
      }
      if (n_out >= min_n && n_out <= max_n) {
        return(.get_colors_fix_hex(rcartocolor::carto_pal(n = n_out, name = carto_name)))
      }
      base_cols <- .get_colors_fix_hex(rcartocolor::carto_pal(n = max_n, name = carto_name))
      if (n_out <= length(base_cols)) return(base_cols[seq_len(n_out)])
      return(.get_colors_fix_hex(grDevices::colorRampPalette(base_cols, space = "Lab")(n_out)))
    }
    if (m %in% names(viridis_opts)) {
      opt <- unname(viridis_opts[m])
      if (requireNamespace("viridisLite", quietly = TRUE)) {
        return(viridisLite::viridis(n_out, option = opt, alpha = alpha))
      }
      if (requireNamespace("viridis", quietly = TRUE)) {
        return(viridis::viridis(n_out, option = opt, alpha = alpha))
      }
      stop("mode='", m, "' requires package 'viridisLite' or 'viridis'.", call. = FALSE)
    }
    if (requireNamespace("RColorBrewer", quietly = TRUE)) {
      info <- RColorBrewer::brewer.pal.info
      rn <- rownames(info)
      idx <- match(m, tolower(rn))
      if (!is.na(idx)) {
        pal_name <- rn[idx]
        maxc <- as_int(info[idx, "maxcolors"])
        maxc <- max(3L, maxc)
        pal <- RColorBrewer::brewer.pal(maxc, pal_name)
        if (n_out <= length(pal)) return(pal[seq_len(n_out)])
        return(.get_colors_fix_hex(grDevices::colorRampPalette(pal, space = "Lab")(n_out)))
      }
    }
    stop("Unsupported mode: ", mode_name, ". Supported: random, rainbow, darkrainbow, lightrainbow, ld, ggplot, all RColorBrewer palettes, all rcartocolor palettes, viridis/magma/plasma/inferno/cividis/rocket/mako/turbo.", call. = FALSE)
  }
  .get_colors_mode_is_continuous_palette <- function(mode_name) {
    m <- tolower(as.character(mode_name)[1])
    if (m %in% c("rainbow", "darkrainbow", "lightrainbow", "ld", "ggplot", names(viridis_opts))) return(TRUE)
    carto_name <- .get_colors_carto_palette_canonical(m)
    if (!is.null(carto_name)) {
      meta <- .get_colors_carto_palette_meta(carto_name)
      if (!is.null(meta) && identical(meta$type, "qualitative")) return(FALSE)
      return(TRUE)
    }
    if (m %in% c("dark2")) return(FALSE)
    if (requireNamespace("RColorBrewer", quietly = TRUE)) {
      info <- RColorBrewer::brewer.pal.info
      rn <- rownames(info)
      idx <- match(m, tolower(rn))
      if (!is.na(idx)) {
        typ <- tolower(as.character(info[idx, "category"]))
        return(!identical(typ, "qual"))
      }
    }
    FALSE
  }
  .get_colors_mode_is_explicit_discrete_palette <- function(mode_name) {
    m <- tolower(as.character(mode_name)[1])
    brewer_discrete <- tolower(c("Dark2", "Paired", "Pastel1", "Pastel2", "Set1", "Set2", "Set3"))
    carto_discrete <- tolower(c("Antique", "Bold", "Pastel", "Prism", "Safe", "Vivid"))
    m %in% c(brewer_discrete, carto_discrete)
  }

  cols <- character()
  if (!discrete) {
    if (mode == "random") {
      if (n == 1L) {
        h <- if (length(include_h)) include_h[1] else stats::runif(1, 0, 360)
        cols <- .get_colors_fix_hex(grDevices::hcl(h, c = stats::runif(1, c_mid[1], c_mid[2]), l = stats::runif(1, l_mid[1], l_mid[2]), alpha = alpha))
      } else {
        h1 <- if (length(include_h)) include_h[1] else stats::runif(1, 0, 360)
        h2 <- if (length(include_h) >= 2L) include_h[2] else (h1 + stats::runif(1, 70, 170)) %% 360
        hm <- (h1 + h2) / 2
        anchors <- c(
          grDevices::hcl(h1, c = stats::runif(1, c_lo[1], c_lo[2]), l = stats::runif(1, l_hi[1], l_hi[2]), alpha = alpha),
          grDevices::hcl(hm, c = stats::runif(1, c_mid[1], c_mid[2]), l = stats::runif(1, l_mid[1], l_mid[2]), alpha = alpha),
          grDevices::hcl(h2, c = stats::runif(1, c_hi[1], c_hi[2]), l = stats::runif(1, l_lo[1], l_lo[2]), alpha = alpha)
        )
        if (length(include_cols)) anchors <- c(anchors[1], include_cols, anchors[length(anchors)])
        cols <- .get_colors_fix_hex(grDevices::colorRampPalette(unique(anchors), space = "Lab")(n))
      }
    } else {
      mode_native_n <- as_int(.get_colors_default_n_from_mode(mode))
      explicit_discrete_mode <- .get_colors_mode_is_explicit_discrete_palette(mode)
      if (isTRUE(explicit_discrete_mode) && is.finite(mode_native_n) && !is.na(mode_native_n) && mode_native_n > 0L) {
        base_native <- .get_colors_palette_from_mode(mode, n_out = mode_native_n)
        if (n <= mode_native_n) {
          cols <- base_native[seq_len(n)]
          cols <- .get_colors_apply_hsv_tone(cols, sat_rng, val_rng, hue_shift = hue_shift_eff, satur_mult = satur_mult_eff, light_offset = light_offset_eff, randomize = random)
        } else {
          pool <- .get_colors_fix_hex(grDevices::colorRampPalette(base_native, space = "Lab")(max(100L, n * 6L)))
          pool <- .get_colors_apply_hsv_tone(pool, sat_rng, val_rng, hue_shift = hue_shift_eff, satur_mult = satur_mult_eff, light_offset = light_offset_eff, randomize = random)
          cols <- .get_colors_pick_stride_hex(pool, n)
        }
      } else {
      base_n_use <- if (identical(mode, "ggplot")) n else palette_base_n
      far_use <- isTRUE(far) && !identical(mode, "ggplot")
      base <- .get_colors_palette_from_mode(mode, n_out = base_n_use)
      if (isTRUE(far_use)) {
        if (isTRUE(random)) {
          pool <- .get_colors_fix_hex(grDevices::colorRampPalette(base, space = "Lab")(max(100L, n * 6L)))
          pool <- .get_colors_apply_hsv_tone(pool, sat_rng, val_rng, hue_shift = hue_shift_eff, satur_mult = satur_mult_eff, light_offset = light_offset_eff, randomize = random)
          cols <- .get_colors_pick_stride_hex(pool, n)
        } else {
          base2 <- .get_colors_apply_hsv_tone(base, sat_rng, val_rng, hue_shift = hue_shift_eff, satur_mult = satur_mult_eff, light_offset = light_offset_eff, randomize = FALSE)
          cols <- .get_colors_pick_stride_hex(base2, n)
        }
      } else {
        if (n <= length(base)) {
          cols <- base[seq_len(n)]
        } else {
          cols <- .get_colors_fix_hex(grDevices::colorRampPalette(base, space = "Lab")(n))
        }
        cols <- .get_colors_apply_hsv_tone(cols, sat_rng, val_rng, hue_shift = hue_shift_eff, satur_mult = satur_mult_eff, light_offset = light_offset_eff, randomize = random)
      }
      }
    }
  } else {
    if (mode == "random") {
      n_inc <- min(length(include_cols), n); n_rem <- n - n_inc
      cols <- character(n); if (n_inc > 0L) cols[seq_len(n_inc)] <- include_cols[seq_len(n_inc)]
      if (n_rem > 0L) {
        hue <- if (length(include_h)) rep_len(include_h, n_rem) + stats::runif(n_rem, -18, 18) else .get_colors_sample_hues_far(n_rem) + stats::rnorm(n_rem, 0, .get_colors_lerp(10, 4, tv))
        hue <- hue %% 360
        chroma <- stats::runif(n_rem, disc_c_min, disc_c_max)
        light <- stats::runif(n_rem, disc_l_min, disc_l_max)
        amp <- .get_colors_lerp(4, 10, tv)
        chroma <- pmin(100, pmax(0, chroma + rep_len(c(amp, -0.7 * amp, 0), n_rem)))
        light <- pmin(100, pmax(0, light + rep_len(c(-0.9 * amp, 0.8 * amp, 0), n_rem)))
        cols[(n_inc + 1L):n] <- .get_colors_fix_hex(grDevices::hcl(hue, c = chroma, l = light, alpha = alpha))
      }
      cols <- .get_colors_fix_hex(cols)[seq_len(n)]
    } else {
      mode_native_n <- as_int(.get_colors_default_n_from_mode(mode))
      explicit_discrete_mode <- .get_colors_mode_is_explicit_discrete_palette(mode)
      if (isTRUE(explicit_discrete_mode) && is.finite(mode_native_n) && !is.na(mode_native_n) && mode_native_n > 0L) {
        base_native <- .get_colors_palette_from_mode(mode, n_out = mode_native_n)
        if (n <= mode_native_n) {
          cols <- base_native[seq_len(n)]
          cols <- .get_colors_apply_hsv_tone(cols, sat_rng, val_rng, hue_shift = hue_shift_eff, satur_mult = satur_mult_eff, light_offset = light_offset_eff, randomize = random)
        } else {
          pool <- .get_colors_fix_hex(grDevices::colorRampPalette(base_native, space = "Lab")(max(length(base_native), n * 6L)))
          pool <- .get_colors_apply_hsv_tone(pool, sat_rng, val_rng, hue_shift = hue_shift_eff, satur_mult = satur_mult_eff, light_offset = light_offset_eff, randomize = random)
          cols <- .get_colors_pick_stride_hex(pool, n)
        }
      } else {
        base <- .get_colors_palette_from_mode(mode, n_out = palette_base_n)
        if (isTRUE(far) || .get_colors_mode_is_continuous_palette(mode)) {
          pool <- .get_colors_fix_hex(grDevices::colorRampPalette(base, space = "Lab")(max(length(base), n * 6L)))
          pool <- .get_colors_apply_hsv_tone(pool, sat_rng, val_rng, hue_shift = hue_shift_eff, satur_mult = satur_mult_eff, light_offset = light_offset_eff, randomize = random)
          cols <- .get_colors_pick_far_hex(pool, n, randomize = random)
        } else {
          if (n <= length(base)) {
            cols <- base[seq_len(n)]
          } else {
            cols <- .get_colors_fix_hex(grDevices::colorRampPalette(base, space = "Lab")(n))
          }
          cols <- .get_colors_apply_hsv_tone(cols, sat_rng, val_rng, hue_shift = hue_shift_eff, satur_mult = satur_mult_eff, light_offset = light_offset_eff, randomize = random)
        }
      }
    }
  }
  if (length(include_cols)) {
    idx <- unique(as_int(round(seq(1, n, length.out = min(length(include_cols), n)))))
    cols[idx] <- include_cols[seq_along(idx)]
  }
  cols <- .get_colors_fix_hex(cols)
  if (isTRUE(is.finite(alpha) && alpha >= 1)) {
    cols <- ifelse(grepl("^#[0-9A-Fa-f]{8}$", cols), substr(cols, 1, 7), cols)
  }

  if (isTRUE(plot)) {
    title0 <- if (!is.null(title) && nzchar(as.character(title)[1])) {
      as.character(title)[1]
    } else {
      sprintf("get.colors (mode=%s, discrete=%s, far=%s, n=%d)", mode, ifelse(discrete, "TRUE", "FALSE"), ifelse(far, "TRUE", "FALSE"), n)
    }
    op <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(op), add = TRUE)
    if (identical(label, "axis")) {
      graphics::par(mar = c(8, 2, 4, 1))
    } else {
      graphics::par(mar = c(2, 2, 4, 1))
    }

    graphics::plot(0, 0, type = "n", xlim = c(0, n), ylim = c(0, 1),
                   xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n", main = title0)
    for (i in seq_len(n)) {
      graphics::rect(i - 1, 0, i, 1, col = cols[i], border = NA)
    }

    if (identical(label, "top")) {
      tc <- .get_colors_text_col_for_bg(cols)
      graphics::text(x = seq_len(n) - 0.5, y = 0.5, labels = cols, col = tc, cex = label.size, srt = 90)
    } else if (identical(label, "axis")) {
      graphics::text(x = seq_len(n) - 0.5, y = -0.08, labels = cols, xpd = NA, srt = 45,
                     adj = 1, cex = label.size, col = "grey20")
    }
  }

  cols
}

.gcanvas_chrom_palette_list <- function() {
  list(
    gray = c("grey60", "grey80"),
    lightgray = c("grey70", "grey90"),
    darkgray = c("grey25", "grey50"),
    grey = c("grey60", "grey80"),
    lightgrey = c("grey70", "grey90"),
    darkgrey = c("grey25", "grey50"),
    lagoon = c("#38a3a5", "#c7f9cc"),
    lightlagoon = c("#69B5B7", "#E3F5E4"),
    darklagoon = c("#268183", "#85C88C"),
    coral = c("#f38375", "#fbc3bc"),
    soda = c("#7DA3C7", "#F2F3DF"),
    lightsoda = c("#90caf9", "#F4F5ED"),
    darksoda = c("#5F80A0", "#edeec9"),
    bbongdda = c("#7DA3C7", "#F2F3DF"),
    lightbbongdda = c("#90caf9", "#F4F5ED"),
    darkbbongdda = c("#5F80A0", "#edeec9"),
    lightcoral = c("#f8ad9d", "#fae0e4"),
    darkcoral = c("#C86254", "#D59288"),
    pink = c("#ff7096", "#f7cad0"),
    lightpink = c("#f9bec7", "#fae0e4"),
    darkpink = c("#D56D93", "#D2979F"),
    lavender = c("#a76eee", "#d6befa"),
    lightlavender = c("#dfc5fe", "#d3d4ff"),
    darklavender = c("#4a4e69", "#9a8c98"),
    random = rev(get.colors(n = 2, tone = 0.1, ash = 0.7, plot = FALSE, silent = TRUE)),
    spectral = get.colors("spectral", 22, plot = FALSE, silent = TRUE),
    viridis = get.colors("viridis", 22, plot = FALSE, silent = TRUE),
    magma = get.colors("magma", 22, plot = FALSE, silent = TRUE),
    plasma = get.colors("plasma", 22, plot = FALSE, silent = TRUE),
    inferno = get.colors("inferno", 22, plot = FALSE, silent = TRUE),
    cividis = get.colors("cividis", 22, plot = FALSE, silent = TRUE),
    rocket = get.colors("rocket", 22, plot = FALSE, silent = TRUE),
    mako = get.colors("mako", 22, plot = FALSE, silent = TRUE),
    turbo = get.colors("turbo", 22, plot = FALSE, silent = TRUE)
  )
}

#' Display a color palette
#'
#' Renders a swatch of `n` colors, either from a named palette/source
#' (`type`) or a vector of hex codes passed via `type`.
#'
#' @param type Palette name, source identifier, or a character vector of hex codes.
#' @param n Number of swatches to draw (when applicable).
#' @param silent Logical. Suppress progress notes.
#'
#' @return Invisibly returns the hex codes; called for its plotting side effect.
#' @export
show.pal <- function(type = NULL, n = 10L, silent = FALSE) {
  require_pkg(c("ggplot2", "data.table"))
  n <- as_int(n)[1]
  if (!is.finite(n) || is.na(n) || n < 2L) n <- 10L

  req <- as.character(type %||% "all")
  req <- tolower(trimws(req))
  req <- req[!is.na(req) & nzchar(req)]
  if (!length(req)) req <- "all"
  req_raw <- as.character(type %||% "all")
  req_raw <- trimws(req_raw)
  req_raw <- req_raw[!is.na(req_raw) & nzchar(req_raw)]

  cat_alias <- c(chrom = "manhattan", chromosome = "manhattan", viridis = "viridis")
  req <- ifelse(req %in% names(cat_alias), unname(cat_alias[req]), req)

  .show_pal_is_color_string <- function(x) {
    if (is.null(x) || !length(x)) return(logical())
    vapply(as.character(x), function(z) {
      !is.null(tryCatch(grDevices::col2rgb(z), error = function(e) NULL))
    }, logical(1))
  }
  .show_pal_to_hex <- function(x) {
    x <- as.character(x)
    out <- vapply(x, function(z) {
      rgb <- tryCatch(grDevices::col2rgb(z), error = function(e) NULL)
      if (is.null(rgb)) return(NA_character_)
      grDevices::rgb(rgb[1, 1], rgb[2, 1], rgb[3, 1], maxColorValue = 255)
    }, character(1))
    out[!is.na(out) & nzchar(out)]
  }

  core_modes <- c("random", "rainbow", "darkrainbow", "lightrainbow", "ld", "ggplot")
  viridis_modes <- c("viridis", "magma", "plasma", "inferno", "cividis", "rocket", "mako", "turbo")
  groups_all <- c("core", "viridis", "rcolorbrewer", "rcartocolor", "manhattan")
  req_groups <- intersect(req, c(groups_all, "all"))

  brewer_names <- character()
  if (requireNamespace("RColorBrewer", quietly = TRUE)) {
    brewer_names <- tolower(rownames(RColorBrewer::brewer.pal.info))
  }
  carto_names <- character()
  if (requireNamespace("rcartocolor", quietly = TRUE)) {
    carto_names <- tolower(unique(as.character(rcartocolor::metacartocolors$Name)))
  }
  man_names <- tolower(names(.gcanvas_chrom_palette_list()))
  pref_man <- req[grepl("^(manhattan\\.|man\\.)", req)]
  pref_suffix <- sub("^(manhattan\\.|man\\.)", "", pref_man)
  pref_suffix <- pref_suffix[pref_suffix %in% man_names]
  pal_names <- unique(c(core_modes, viridis_modes, brewer_names, carto_names, man_names))
  req_pals <- setdiff(req, c(groups_all, "all"))
  req_pals <- req_pals[req_pals %in% pal_names]
  req_pals <- unique(c(req_pals, pref_suffix))
  is_custom <- length(req_raw) > 0L && all(.show_pal_is_color_string(req_raw))

  if (!is_custom && !("all" %in% req) && length(req_raw) > 0L && length(req_groups) == 0L && length(req_pals) == 0L) {
    stop("Unknown palette/category: ", paste(req_raw, collapse = ", "), call. = FALSE)
  }

  wanted_groups <- if ("all" %in% req_groups || (!length(req_groups) && !length(req_pals))) groups_all else setdiff(req_groups, "all")
  if (is_custom) wanted_groups <- character()

  mk_dt <- function(name, cols, grp) {
    cc <- as.character(cols)
    cc <- cc[!is.na(cc) & nzchar(cc)]
    if (!length(cc)) return(NULL)
    data.table::data.table(
      palette = as.character(name)[1],
      group = as.character(grp)[1],
      idx = seq_along(cc),
      color = cc
    )
  }

  dts <- list()

  if ("core" %in% wanted_groups) {
    for (m in core_modes) {
      cols <- tryCatch(get.colors(m, n = n, plot = FALSE, silent = TRUE), error = function(e) character())
      dts[[length(dts) + 1L]] <- mk_dt(m, cols, "core")
    }
  }
  if ("viridis" %in% wanted_groups) {
    for (m in viridis_modes) {
      cols <- tryCatch(get.colors(m, n = n, plot = FALSE, silent = TRUE), error = function(e) character())
      dts[[length(dts) + 1L]] <- mk_dt(m, cols, "Viridis")
    }
  }
  if ("rcolorbrewer" %in% wanted_groups && requireNamespace("RColorBrewer", quietly = TRUE)) {
    binfo <- RColorBrewer::brewer.pal.info
    br_names <- rownames(binfo)
    for (nm in br_names) {
      nn <- min(max(3L, as_int(binfo[nm, "maxcolors"])), max(3L, n))
      cols <- tryCatch(get.colors(nm, n = nn, plot = FALSE, silent = TRUE), error = function(e) character())
      dts[[length(dts) + 1L]] <- mk_dt(nm, cols, "RColorBrewer")
    }
  }
  if ("rcartocolor" %in% wanted_groups && requireNamespace("rcartocolor", quietly = TRUE)) {
    cmeta <- rcartocolor::metacartocolors
    cnames <- unique(as.character(cmeta$Name))
    for (nm in cnames) {
      i <- which(as.character(cmeta$Name) == nm)[1]
      max_n <- as_int(cmeta$Max_n[i]); if (!is.finite(max_n) || is.na(max_n) || max_n < 2L) max_n <- n
      min_n <- as_int(cmeta$Min_n[i]); if (!is.finite(min_n) || is.na(min_n) || min_n < 2L) min_n <- 2L
      nn <- min(max_n, max(min_n, n))
      cols <- tryCatch(get.colors(nm, n = nn, plot = FALSE, silent = TRUE), error = function(e) character())
      dts[[length(dts) + 1L]] <- mk_dt(nm, cols, "rcartocolor")
    }
  }
  if ("manhattan" %in% wanted_groups) {
    man <- .gcanvas_chrom_palette_list()
    dep_pal <- c("spectral", "viridis", "magma", "plasma", "inferno", "cividis", "rocket", "mako", "turbo")
    mk <- names(man)
    for (i in seq_along(man)) {
      if (mk[i] %in% dep_pal) next
      dts[[length(dts) + 1L]] <- mk_dt(paste0("manhattan.", mk[i]), man[[i]], "manhattan")
    }
  }

  if (length(req_pals) > 0L) {
    for (pm in req_pals) {
      pal_label <- pm
      if (pm %in% man_names) {
        cols <- .gcanvas_chrom_palette_list()[[pm]]
        pal_label <- paste0("manhattan.", pm)
      } else {
        cols <- tryCatch(get.colors(pm, n = n, plot = FALSE, silent = TRUE), error = function(e) character())
      }
      dts[[length(dts) + 1L]] <- mk_dt(pal_label, cols, "selected")
    }
  }
  if (is_custom) {
    custom_cols <- .show_pal_to_hex(req_raw)
    dts[[length(dts) + 1L]] <- mk_dt("custom", custom_cols, "selected")
  }

  dts <- Filter(Negate(is.null), dts)
  if (!length(dts)) stop("No palettes available for the requested type.", call. = FALSE)
  dt <- data.table::rbindlist(dts, use.names = TRUE, fill = TRUE)
  dt <- dt[!is.na(color) & nzchar(color)]
  if (!nrow(dt)) stop("No palettes available for the requested type.", call. = FALSE)
  show_hex_overlay <- isTRUE(is_custom) || (length(req_pals) == 1L && length(req_groups) == 0L)
  .show_pal_text_col_for_bg <- function(hex) {
    rgb <- grDevices::col2rgb(hex) / 255
    y <- 0.2126 * rgb[1, ] + 0.7152 * rgb[2, ] + 0.0722 * rgb[3, ]
    ifelse(y < 0.55, "white", "black")
  }
  dt[, txt_col := .show_pal_text_col_for_bg(color)]
  dt[, group := factor(group, levels = unique(group))]
  dt[, palette_lab := as.character(palette)]
  ord_key <- unique(paste(dt$group, dt$palette_lab, sep = ":::"))
  dt[, palette_key := paste(group, palette_lab, sep = ":::")]
  dt[, palette_fac := factor(palette_key, levels = rev(ord_key))]
  lab_map <- stats::setNames(dt$palette_lab, dt$palette_key)

  if (isTRUE(show_hex_overlay)) {
    p <- ggplot2::ggplot(dt, ggplot2::aes(x = idx, y = palette_fac)) +
      ggplot2::geom_tile(ggplot2::aes(fill = color), height = 0.82, width = 0.95, color = NA) +
      ggplot2::geom_text(
        ggplot2::aes(label = toupper(color), color = txt_col),
        angle = 90,
        size = 3.5,
        fontface = "bold"
      ) +
      ggplot2::scale_y_discrete(labels = lab_map, expand = c(0, 0)) +
      ggplot2::scale_fill_identity() +
      ggplot2::scale_color_identity() +
      ggplot2::labs(x = NULL, y = NULL, title = "gcanvas palettes") +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        panel.border = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_text(size = 9),
        axis.ticks.y = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(face = "bold", size = 12),
        strip.background = ggplot2::element_blank(),
        strip.text = ggplot2::element_blank()
      )
  } else {
    p <- ggplot2::ggplot(dt, ggplot2::aes(x = idx, y = palette_fac)) +
      ggplot2::geom_tile(ggplot2::aes(fill = color), height = 0.82, width = 0.95, color = NA) +
      ggplot2::facet_grid(group ~ ., scales = "free_y", space = "free_y") +
      ggplot2::scale_y_discrete(labels = lab_map) +
      ggplot2::scale_fill_identity() +
      ggplot2::scale_color_identity() +
      ggplot2::labs(
        x = "Color Index",
        y = NULL,
        title = "gcanvas palettes"
      ) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_text(size = 8),
        plot.title = ggplot2::element_text(face = "bold", size = 12),
        axis.ticks.y = ggplot2::element_blank(),
        strip.background = ggplot2::element_rect(fill = "grey95", color = "grey80"),
        strip.text.y = ggplot2::element_text(face = "bold", size = 9)
      )
  }

  if (!isTRUE(silent)) {
    print_target <- NULL
    if (is_custom) {
      print_target <- unique(dt[group == "selected" & palette == "custom"]$color)
    } else if (length(req_pals) == 1L && length(req_groups) == 0L) {
      print_target <- unique(dt[group == "selected" & palette == req_pals[1]]$color)
    }
    if (!is.null(print_target) && length(print_target) > 0L) {
      cat(paste(toupper(print_target), collapse = ", "), "\n")
    }
  }
  p
}

.gcanvas_pretty_path <- function(x) {
  if (is.null(x) || length(x) == 0) return(x)
  x <- as.character(x)
  fix1 <- function(s) {
    if (is.na(s) || !nzchar(s)) return(s)
    if (grepl("^//", s)) {
      return(paste0("//", gsub("/{2,}", "/", substr(s, 3, nchar(s)))))
    }
    gsub("(?<!:)/{2,}", "/", s, perl = TRUE)
  }
  vapply(x, fix1, character(1))
}

.gcanvas_as_num2 <- function(x) {
  if (is.list(x)) x <- unlist(x, use.names = FALSE)
  suppressWarnings(as.numeric(x))
}

.gcanvas_as_char_no_null <- function(x, empty = "") {
  if (is.null(x)) return(character())

  if (is.list(x)) {
    return(vapply(x, function(z) {
      if (is.null(z) || length(z) == 0L) return(empty)
      z0 <- as.character(z)[1]
      if (is.na(z0) || !nzchar(z0)) empty else z0
    }, character(1)))
  }

  out <- as.character(x)
  out[is.na(out) | !nzchar(out)] <- empty
  out
}

