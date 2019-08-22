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

#' Add variables to a data frame with given names
#'
#' Helper for conflictsToCharacters function
#'
#' @param df data.frame
#' @param names names of variables to add
#' @param verbose control messaging
#' @importFrom stats setNames
addVars <- function(df, names, verbose = FALSE) {
  if (length(names) == 0)
    return(df)

  newcols <- matrix(NA, nrow = nrow(df), ncol = length(names)) %>%
    as.data.frame() %>%
    setNames(names)

  if (verbose)
    message(sprintf("Adding the following columns: %s",
                    paste(names, collapse = ", ")))

  out <- cbind(df, newcols)
  out
}

#' Used in bind_rows2
#'
#' Converts conflicted column types to character columns and rbinds the result.
#' If one data.frame is missing columns, they are created and populated by NA's.
#'
#' @param df1 data.frame (or similar) to combine
#' @param df2 data.frame (or similar) to combine
#' @param addMissing Add columns to result if they are missing in either df1 or df2?
#' @param verbose Include messages about conversions?

conflictsToCharacter <- function (df1, df2, addMissing = TRUE, verbose = FALSE)
{
  df1 <- as.data.frame(df1, stringsAsFactors = FALSE)
  df2 <- as.data.frame(df2, stringsAsFactors = FALSE)

  allcols <- union(names(df1), names(df2))
  newTo1 <- setdiff(allcols, names(df1))
  newTo2 <- setdiff(allcols, names(df2))

  if (addMissing) {
    df1 <- addVars(df1, newTo1, verbose = verbose)[allcols]
    df2 <- addVars(df2, newTo2, verbose = verbose)[allcols]
  } else if (length(newTo1) > 0 || length(newTo2) > 0)
    stop("Column names must match unless addMissing is TRUE")

  conflicts = !unlist(Map(identical, lapply(df1, class),
                          lapply(df2, class)))
  conflicts[c(newTo1, newTo2)] <- FALSE

  if (sum(conflicts) > 0 && verbose)
    message(paste("Converting to character:", names(conflicts)[conflicts]))
  df1[conflicts] = lapply(df1[conflicts], as.character)
  df2[conflicts] = lapply(df2[conflicts], as.character)
  rbind(df1, df2, stringsAsFactors = FALSE)
}

#' A more flexible (but slower) version of dplyr::bind_rows
#'
#' @param dfList a list of data.frames
#' @param addMissing Add columns if missing from a subset of dfList?
#' @param verbose Include messages about conversions?
#' @export

bind_rows2 <- function (dfList, addMissing = TRUE, verbose = FALSE) {
  redfun <- function(x, y)
    conflictsToCharacter(x, y, addMissing = addMissing, verbose = verbose)
  out <- Reduce(redfun, dfList)
  out
}

#' Check if sf package is installed.
#'
check_sf <- function() {
  requireNamespace("sf", quietly = TRUE)
}

#' Check if geosphere package is installed.
#'
check_geosphere <- function() {
  requireNamespace("sf", quietly = TRUE)
}
