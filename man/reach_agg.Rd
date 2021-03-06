% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reach_agg.R
\name{reach_agg}
\alias{reach_agg}
\title{Aggregate node-level data to reach-scale}
\usage{
reach_agg(nodedata, usenodes = "all", weight = TRUE)
}
\arguments{
\item{nodedata}{Must have additional columns from \code{add_nodelen(), add_offset()}}

\item{usenodes}{Which nodes to use. Defaults to "all"}

\item{weight}{Use weighted regression based on height uncertainty?}
}
\description{
A slightly imprecise reproduction of RiverObs reach aggregation for a handful of
 variables.
}
