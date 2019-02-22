# Plotting functions

#' Plot a histogram of errors, possibly scaled by uncertainty estimates
#'
#' Currently only implemented for node-level data
#'
#' @param valdata data.frame as returned by code{rt_valdata()}
#' @param center Subtract the mean (bias-correct) the errors?
#' @param scale Scale the errors by 1-sigma uncertainty estimates?
#' @param curve Overlay a standard normal curve? This is the default if
#'   both \code{center} and \code{scale} are \code{TRUE}.
#' @param vars Which variables to plot? Defaults to "all"
#' @param ... Passed to \code{rt_valdata()}
#' @importFrom dplyr group_by ungroup
#' @export

rt_val_hist <- function(valdata, center = FALSE, scale = FALSE,
                        curve = center && scale,
                        vars = "all",
                        ...) {
  # Currently only implemented for nodes

  if (length(vars) > 1 || vars != "all") {
    valdata <- valdata %>%
      dplyr::filter(variable %in% vars)
  }

  out <- valdata %>%
    mutate(err = pixc_err) %>%
    group_by(variable)


  if (center) {
    out <- mutate(out, err = err - mean(err, na.rm = TRUE))
  }

  if (scale) {
    out <- mutate(out, err = err / sigma_est)
  }

  out <- out %>%
    ungroup() %>%
    ggplot(aes(x = err)) +
    geom_histogram(aes(y = ..density..), bins = 15) +
    facet_wrap(~variable, scales = "free")

  if(curve) {
    out <- out +
      stat_function(fun = dnorm, args = list(mean = 0, sd = 1), color = "blue")
  }
  out
}

#' Plot predictions and uncertainty bounds across nodes
#'
#' @param valdata data.frame as returned by code{rt_valdata()}
#' @param variable Which variable to plot?
#' @param err_only Only plot the errors (subtract off truth value)?
#' @param ... Passed to \code{rt_valdata()}
#'
#' @importFrom dplyr mutate
#' @export

rt_val_nodeseries <- function(valdata, variable = "height", err_only = TRUE, ...) {

  valdata <- valdata %>%
    `[`(.$variable == variable, ) %>%
    mutate(value = pixc_val, err = pixc_err)

  if (err_only) {
    valdata$baseval <- 0
  } else {
    valdata$baseval <- valdata$value
  }
  # return(valdata)
  yvar <- ifelse(err_only, "err", "value")

  out <- valdata %>%
    ggplot(aes(x = node_id)) +
    geom_ribbon(aes(ymin = baseval - 1.96 * sigma_est,
                    ymax = baseval + 1.96 * sigma_est), fill = "pink") +
    geom_ribbon(aes(ymin = baseval - sigma_est,
                    ymax = baseval + sigma_est), fill = "#7780ff") +
    geom_point(aes_string(y = yvar)) +
    theme_bw()
  out
}

#' Returns a vector of worst-performing nodes (by error)
#'
#' @export
badnodes <- function(valdata, variable = "width", n = 4,
                     which = c("abs", "min", "max")) {
  which <- match.arg(which)
  valdata <- valdata[valdata[["variable"]] == variable, ]
  errvec <- valdata$pixc_err
  if (which == "abs") errvec <- -abs(errvec) else
    if (which == "max") errvec <- -errvec

  badords <- order(errvec)[1:n]

  out <- valdata$node_id[badords]
  out
}

#' Map a given number of nodes' pixcvec locations--gdem versus rivertile
#'
#' @param valdata data.frame as returned by \code{rt_valdata()}
#' @param nodes A vector of node_id's, defaults to worst nodes using \code{badnodes()}
#' @param pcv1 Name of first pixcvec netcdf
#' @param pcv2 Name of second pixcvec netcdf
#' @param maxpixels Maximum number of pixels to plot, as a random sample if necessary.
#' @param leaflet Use interactive leaflet plot? Otherwise uses ggplot.
#' @param ... Passed to \code{rt_valdata}
#'
#' @importFrom fs path
#' @importFrom dplyr rename
#'
#' @export
val_map_node <- function(dir, nodes = badnodes(rt_valdata(dir)),
                         pcv1 = "pcv.nc", pcv2 = "pcv_gdem.nc",
                         maxpixels = 1000, leaflet = TRUE, ...) {

  pcvdata1 <- pixcvec_read(path(dir, pcv1)) %>%
    dplyr::filter(node_index %in% nodes) %>%
    rename(lat = latitude_vectorproc, lon = longitude_vectorproc)
  pcvdata2 <- pixcvec_read(path(dir, pcv2)) %>%
    dplyr::filter(node_index %in% nodes) %>%
    rename(lat = latitude_vectorproc, lon = longitude_vectorproc)

  if (nrow(pcvdata1) > maxpixels) {
    message("Subsampling pixcvec 1. Change using `maxpixels` argument.")
    pcvdata1 <- dplyr::sample_n(pcvdata1, maxpixels)
  }
  if (nrow(pcvdata2) > maxpixels) {
    message("Subsampling pixcvec 2. Change using `maxpixels` argument.")
    pcvdata2 <- dplyr::sample_n(pcvdata2, maxpixels)
  }

  if (leaflet) {
    out <- leaflet() %>%
      addTiles() %>%
      addCircleMarkers(data = pcvdata1, lng = ~lon, lat = ~lat,
                       radius = 2, color = "blue") %>%
      addCircleMarkers(data = pcvdata2, lng = ~lon, lat = ~lat,
                       radius = 2, color = "red")
  } else {
    bbox <- ggmap::make_bbox(pcvdata1$lon, pcvdata1$lat)
    satmap <- ggmap::get_map(location = bbox, maptype = "terrain", source = "google")
    out <- ggmap(satmap) +
      geom_point(data = pcvdata1, aes(x = lon, y = lat), color = "blue", size = 2) +
      geom_point(data = pcvdata2, aes(x = lon, y = lat), color = "red", size = 2)
    out
  }

  out
}
