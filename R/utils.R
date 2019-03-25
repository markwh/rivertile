# utils.R
# misc. utility functions

#' Copied from markstats package
#'
#' @param strvec Character vector
#' @param split Character to split by
#' @param piece of split element to keep (numeric)
#' @param ... Passed to \code{strsplit()}
#'
splitPiece <- function (strvec, split, piece, ...) {
  spl <- strsplit(strvec, split = split, ...)
  out <- vapply(spl, `[`, character(1), piece)
  out
}

