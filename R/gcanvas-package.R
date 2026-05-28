#' gcanvas: Genomic Visualization and Analysis Toolkit
#'
#' A toolkit for visualizing and analyzing genomic association data.
#' Plotting functions include [manhattan()], [qq()], [regional()], [miami()],
#' [circos()] (experimental), [hetmiss()], and the helper [regional.track()].
#' PCA workflows are provided by [pca()], [pca.projection()], [pca.plot()],
#' [pca.screeplot()], [pca.refqc()], and [pca.ancestry()]. Genomic helpers
#' include [liftover()], [calcld()], [ldproxy()], [ldclump()], [is.novel()],
#' [get.lead()], [beta.align()], [maf.align()], [ld.align()], [plink.extract()],
#' [get.pilot()], and [geneinfo()]. General utilities cover color palettes
#' ([get.colors()], [show.pal()]), chromosome handling
#' ([sort.chrom()], [rank.chrom()], [normalize.chrom()]), p-value handling
#' ([pvalue()], [zabs()], [quantvec()], [credibleset()], [invnorm()]),
#' matrix conversion ([dt2mat()], [mat2dt()]), GTF processing ([gtf2rds()]),
#' and the [logger()] sink.
#'
#' Reference annotation tracks are bundled as lazy-loaded datasets
#' ([tracks.b37], [tracks.b38]).
#'
#' @keywords internal
"_PACKAGE"

#' Null-coalescing operator
#'
#' Returns `a` if it is not `NULL`, otherwise returns `b`.
#'
#' @param a Left operand.
#' @param b Right operand (fallback when `a` is `NULL`).
#' @return `a` when `!is.null(a)`, else `b`.
#' @name grapes-or-or-grapes
#' @keywords internal
`%||%` <- function(a, b) if (!is.null(a)) a else b

#' String concatenation operator
#'
#' Convenience infix wrapper around [paste0()].
#'
#' @param a,b Character vectors to concatenate.
#' @return A character vector of `paste0(a, b)`.
#' @name grapes-plus-grapes
#' @keywords internal
`%+%` <- function(a, b) paste0(a, b)
