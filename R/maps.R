# Map pixel and rivertile products

#' Convert linear distance to lat/lon distance
#'
#' Helper for map_node_pixc.
#'
#' @param x distance in meters
#' @param lat latitude, used for determining degree length
to_latlon <- function(x, lat) {
  lon_m <- cos(lat * pi / 180) * 110567
  lat_m <- 111000 # approximately
  out <- x * sqrt(1 / (lat_m * lon_m))
  out
}

#' Filter a pixel cloud
#'
#' Only include pixels assigned to a given vector of nodes.
#'
#' @param pixdf as returned by \code{pixc_read()} or \code{pixc_join()}
#' @param nodeids values of node_index to keep
#' @param pcvdf a data.frame that includes node and range/azimuth information,
#'  e.g. as returned by \code{pixcvec_read()}. May be omitted if this information
#'  is already in \code{pixdf}
#' @export
grab_nodes <- function(pixdf, nodeids, pcvdf = pixdf) {
  radf <- pcvdf %>%
    dplyr::filter(node_index %in% nodeids)

  out <- pixdf %>%
    dplyr::filter((range_index * 1e5 + azimuth_index) %in%
             (radf$range_index * 1e5 + radf$azimuth_index))
  out
}

#' Filter a pixel data frame to a range/azimuth bounding box.
#'
#' @describeIn grab_nodes Filter to a range/azimuth bounding box.
#' @inheritParams grab_nodes
#' @param dilate Number of pixels to extend beyond the supplied nodes
#'   (default is 0).
#' @export
grab_nodebox <- function(pixdf, nodeids, pcvdf = pixdf, dilate = 0) {
  radf <- pcvdf %>%
    dplyr::filter(node_index %in% nodeids)

  minrange <- min(radf$range_index, na.rm = TRUE) - dilate
  maxrange <- max(radf$range_index, na.rm = TRUE) + dilate
  minazimuth <- min(radf$azimuth_index, na.rm = TRUE) - dilate
  maxazimuth <- max(radf$azimuth_index, na.rm = TRUE) + dilate

  out <- pixdf %>%
    dplyr::filter(range_index >= minrange, range_index <= maxrange,
           azimuth_index >= minazimuth, azimuth_index <= maxazimuth)
  out
}

#' Create a bounding box (lat/lon) from pixc and node ids.
#'
#' @describeIn grab_nodes Create a lat/lon bounding box.
#' @inheritParams grab_nodes
#' @param dilate_frac fraction of lat/lon range to dilate the bounding box
#'   (applied on both sides)
#' @export
rt_nodebbox <- function(pixdf, nodeids, pcvdf = pixdf, dilate = 0) {
  radf <- pcvdf %>%
    dplyr::filter(node_index %in% nodeids)

  minlat <- min(radf$latitude, na.rm = TRUE)
  maxlat <- max(radf$latitude, na.rm = TRUE)
  minlon <- min(radf$longitude, na.rm = TRUE)
  maxlon <- max(radf$longitude, na.rm = TRUE)

  latrange <- maxlat - minlat
  lonrange <- maxlon - minlon

  out <- list(minlat = minlat - latrange * dilate,
              maxlat = maxlat + latrange * dilate,
              minlon = minlon - lonrange * dilate,
              maxlon = maxlon + lonrange * dilate)
  out
}

#' Filter a pixel data frame to a lat/lon bounding box.
#'
#' @describeIn grab_nodes Filter to a lat/lon bounding box.
#' @inheritParams grab_nodes
#' @param dilate_frac fraction of lat/lon range to dilate the bounding box
#'   (applied on both sides)
#' @param bbox as returned by \code{rt_nodebbox()}
#' @export
grab_bbox <- function(pixdf, bbox) {
  out <- pixdf %>%
    dplyr::filter(latitude >= bbox$minlat, longitude >= bbox$minlon,
                  latitude <= bbox$maxlat, longitude <= bbox$maxlon)
  out
}

#' Display a pixel cloud as a map.
#'
#' @param pixdf As returned by \code{pixc_read()}, possibly joined by
#'   \code{pixc_join()}
#' @param colorby variable name in \code{pixcdf} to use for coloring
#' @param geoloc Which geolocation to use: "best" (default) chooses based on
#'   which columns are available in \code{pixcdf}.
#' @param real_area If \code{TRUE}, use display actual pixel areas. Note that
#'   this will take considerably more computation time.
#' @param water_frac Scale pixel size by \code{water_frac} column?
#' @param plot if FALSE return the plot data instead of generating the plot.
#' @param maxpoints Maximum number of points to allow. Meant to encourage
#'   filtering so as not to overwhelm the renderer.
#' @param ... passed to \code{ggplot2::geom_point} or \code{ggforce::geom_circle}
#'
#' @importFrom ggplot2 ggplot geom_point scale_size_identity labs
#' @export
pixc_map <- function(pixdf,
                     colorby = "classification",
                     geoloc = c("best", "orig", "improved"),
                     real_area = FALSE,
                     water_frac = FALSE,
                     plot = TRUE,
                     maxpoints = 2500, ...) {

  geoloc <- match.arg(geoloc)

  if (nrow(pixdf) > maxpoints)
    stop(sprintf("Number of pixels (%s) is greater than maxpoints (%s).\n",
                 nrow(pixdf), maxpoints),
         "Please filter pixdf (e.g. grab_nodes()) or increase maxpoints")

  # geolocation:
  if (geoloc == "improved" ||
      (geoloc == "best" && !is.null(pixdf$latitude_vectorproc))) {
    pixdf$latitude <- pixdf$latitude_vectorproc
    pixdf$longitude <- pixdf$longitude_vectorproc
  }

  # sizing, coloring of points/circles
  pixdf$sizescale <- 1
  if (water_frac)  pixdf$sizescale <- pixdf$water_frac
  pixdf$colorvar <- pixdf[[colorby]]


  # Construct ggplot object
  mapgg <- pixdf %>%
    ggplot()

  if (real_area) {

    # convert pixel area to radius in meters, then to lat/lon
    pixradius_m <- sqrt(pixdf$pixel_area * pixdf$sizescale / pi)

    pixdf$radius_ll <- to_latlon(pixradius_m, pixdf$latitude)

    if (!plot) return(pixdf)

    mapgg <- mapgg +
      ggforce::geom_circle(aes(x0 = longitude, y0 = latitude,
                      fill = colorvar,
                      r = radius_ll),
                  data = pixdf, n = 8, linetype = 0, ...)
  } else { # use points instead of circles
    if (!plot) return(pixdf)
    pixcsize = 2.5
    mapgg <- mapgg +
      geom_point(aes(x = longitude, y = latitude,
                     color = colorvar,
                     size = pixcsize * sqrt(sizescale)),
                 shape = 20, ...) +
      scale_size_identity() +
      labs(color = colorby)
  }

  mapgg
}

#' Slant-plane map of pixel cloud
#'
#' @describeIn pixc_map Slant-plane map of pixel cloud
#' @inheritParams pixc_map
#' @importFrom ggplot2 geom_raster aes
#' @export
pixc_slantmap <- function(pixdf,
                          colorby = "classification",
                          maxpoints = 5000, ...) {
  # coloring of points/circles
  pixdf$colorvar <- pixdf[[colorby]]

  # Construct ggplot object
  mapgg <- pixdf %>%
    ggplot()

  pixcsize = 2.5
  mapgg <- mapgg +
    geom_raster(aes(x = range_index, y = azimuth_index,
                    fill = colorvar), ...) +
    labs(fill = colorby)

  mapgg
}
