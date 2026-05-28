# Build lazy-loaded reference data objects from legacy .rds files.
#
# Run once interactively from the package root:
#   setwd("/Users/bsu/scripts/gcanvas")
#   source("data-raw/build_tracks.R")
#
# This produces data/tracks.b37.rda and data/tracks.b38.rda which become
# available to users via `data(tracks.b37)` / `data(tracks.b38)` once
# the package is installed. The package's `LazyData: true` setting means
# the objects are only loaded into memory on first reference.

stopifnot(file.exists("gcanvas.tracks.b37.rds"))
stopifnot(file.exists("gcanvas.tracks.b38.rds"))

tracks.b37 <- readRDS("gcanvas.tracks.b37.rds")
tracks.b38 <- readRDS("gcanvas.tracks.b38.rds")

if (requireNamespace("usethis", quietly = TRUE)) {
  usethis::use_data(tracks.b37, tracks.b38, compress = "xz", overwrite = TRUE)
} else {
  dir.create("data", showWarnings = FALSE)
  save(tracks.b37, file = "data/tracks.b37.rda", compress = "xz")
  save(tracks.b38, file = "data/tracks.b38.rda", compress = "xz")
}
