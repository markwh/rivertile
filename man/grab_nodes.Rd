% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/maps.R
\name{grab_nodes}
\alias{grab_nodes}
\alias{grab_nodebox}
\alias{rt_nodebbox}
\alias{grab_bbox}
\title{Filter a pixel cloud}
\usage{
grab_nodes(pixdf, nodeids, pcvdf = pixdf)

grab_nodebox(pixdf, nodeids, pcvdf = pixdf, dilate = 0)

rt_nodebbox(pixdf, nodeids, pcvdf = pixdf, dilate_frac = 0)

grab_bbox(pixdf, bbox)
}
\arguments{
\item{pixdf}{as returned by \code{pixc_read()} or \code{pixc_join()}}

\item{nodeids}{values of node_index to keep}

\item{pcvdf}{a data.frame that includes node and range/azimuth information,
e.g. as returned by \code{pixcvec_read()}. May be omitted if this information
is already in \code{pixdf}}

\item{dilate}{Number of pixels to extend beyond the supplied nodes
(default is 0).}

\item{dilate_frac}{fraction of lat/lon range to dilate the bounding box
(applied on both sides)}

\item{bbox}{as returned by \code{rt_nodebbox()}}
}
\description{
Only include pixels assigned to a given vector of nodes.
}
\section{Functions}{
\itemize{
\item \code{grab_nodebox}: Filter to a range/azimuth bounding box.

\item \code{rt_nodebbox}: Create a lat/lon bounding box.

\item \code{grab_bbox}: Filter to a lat/lon bounding box.
}}

