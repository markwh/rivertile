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

#' Validation scatterplot
#'
#' @param valdata As returned by \code{rt_valdata()}
#' @param variables Variables to plot
#' @param xvar,yvar x- and y-axis variables
#' @param ci1,ci2 Confidence-interval levels for 2 ribbon geoms.
#'  Either or can be set to FALSE to disable plotting.
#' @param plot if FALSE, just return the data used for plotting
#' @export
rt_val_scatter <- function(valdata, variables = c("height", "width", "slope"),
                           xvar = "id", yvar = c("value", "err", "relerr"),
                           ci1 = 0.6827, ci2 = 0.95,
                           plot = TRUE) {
  varnames <- c(value = "pixc_val", err = "pixc_err", relerr = "rel_err")
  yvar <- match.arg(yvar)
  yvarname <- yvar # for plot axis label
  yvar <- varnames[yvar]

  if (xvar == "id") xvar <- ifelse(is.null(valdata$node_id),
                                   "reach_id", "node_id")
  plotdata <- valdata %>%
    dplyr::filter(variable %in% variables) %>%
    dplyr::mutate(rel_err = pixc_err / sigma_est)

  # Manually add x and y axis variables based on inputs
  plotdata[["xval"]] <- plotdata[[xvar]]
  plotdata[["yval"]] <- plotdata[[yvar]]

  if (yvar == "pixc_val") {
    plotdata[["ymiddle"]] <- plotdata[["gdem_val"]]
  } else {plotdata[["ymiddle"]] <- 0}
  if (yvar == "rel_err") {
    plotdata$ysigma <- 1
  } else {plotdata$ysigma <- plotdata$sigma_est}

  # Sigma multipliers for confidence intervals
  if (ci1 > 0) {
    stopifnot(ci1 < 1)
    sigmult1 <- qnorm((1 + ci1) / 2)
    plotdata <- plotdata %>%
      mutate(ci1_lwr = ymiddle - sigmult1 * ysigma,
             ci1_upr = ymiddle + sigmult1 * ysigma)
  }
  if (ci2 > 0) {
    stopifnot(ci2 < 1)
    sigmult2 <- qnorm((1 + ci2) / 2)
    plotdata <- plotdata %>%
      mutate(ci2_lwr = ymiddle - sigmult2 * ysigma,
             ci2_upr = ymiddle + sigmult2 * ysigma)
  }

  # Return the data or build the plot
  if (!plot) return(plotdata)

  out <- ggplot(plotdata, aes(x = xval))

  # Add ribbons
  if (ci1 > 0) {
    out <- out + geom_ribbon(aes(ymin = ci2_lwr, ymax = ci2_upr),
                             fill = "pink")
  }
  if (ci2 > 0) {
    out <- out + geom_ribbon(aes(ymin = ci1_lwr, ymax = ci1_upr),
                             fill = "#7780ff")
  }

  # Add points, wrap variables, label
  out <- out + geom_point(aes(y = yval)) +
    facet_wrap(~variable, scales = "free") +
    ylab(yvarname) + xlab(xvar)

  out
}


#' Returns a vector of worst-performing nodes (by error)
#'
#' @param valdata As returned by \code{rt_valdata()}
#' @param variable Which variable's errors define "bad" nodes?
#' @param n Number of bad nodes to return
#' @param which "abs" for worst absolute errors, "min" for worst
#'  negative errors, "max" for worst positive errors
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

val_node_data <- function(dir, nodes = badnodes(rt_valdata(dir)),
                          pcv1 = "pcv.nc", pcv2 = "pcv_gdem.nc") {
  pcvdata1 <- pixcvec_read(path(dir, pcv1)) %>%
    dplyr::filter(node_index %in% nodes) %>%
    rename(lat = latitude_vectorproc, lon = longitude_vectorproc)
  pcvdata2 <- pixcvec_read(path(dir, pcv2)) %>%
    dplyr::filter(node_index %in% nodes) %>%
    rename(lat = latitude_vectorproc, lon = longitude_vectorproc)

  out <- list(pcv1 = pcvdata1, pcv2 = pcvdata2)
  out
}

#' Map a given number of nodes' pixcvec locations--gdem versus rivertile
#'
#' @param dir A directory containing rivertile output.
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

  datalist <- val_node_data(dir, nodes, pcv1, pcv2)
  pcvdata1 <- datalist[["pcv1"]]
  pcvdata2 <- datalist[["pcv2"]]

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
                       radius = 2, color = "blue",
                       popup = ~paste0("node: ", node_index)) %>%
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



# Area plot ---------------------------------------------------------------

#' Join pixcvec to pixel cloud
#'
#' @param dir Directory containing pixel_cloud netcdfs
#' @param pcvname,pixcname Names of pixcvec and pixel cloud netcdf files
#' @export
join_pixc <- function(dir, pcvname = "pcv.nc",
                      pixcname = "pixel_cloud.nc") {

  pcvdf <- pixcvec_read(path(dir, pcvname))
  pixcdf <- pixc_read(path(dir, pixcname)) %>%
    inner_join(pcvdf, by = c("azimuth_index", "range_index"))
  pixcdf
}

#' Plot the area computation for a node
#'
#' @param pixc_joined A joined pixc data frame, as returned by \code{join_pixc()}
#' @param nodes Vector giving indices of nodes to plot
#' @param node_truth Optional truth data, as returned by \code{rt_read()}
#' @param plot if FALSE, return the plot data but don't construct the ggplot
#' @export
plot_area <- function(pixc_joined, nodes, node_truth = NULL, plot = TRUE) {
  sumrydf <- pixc_joined %>%
    filter(node_index %in% nodes) %>%
    group_by(node_index) %>%
    arrange(desc(water_frac)) %>%
    mutate(cum_area = cumsum(pixel_area),
           area_lag = dplyr::lag(cum_area, default = 0),
           classification = as.factor(classification)) %>%
    ungroup()

  if (!is.null(node_truth)) {
    joindf <- node_truth %>%
      transmute(reach_index = reach_id, node_index = node_id,
                true_area = area_total)
    sumrydf <- sumrydf %>%
      left_join(joindf, by = c("node_index", "reach_index"))
  }

  if (!plot) return(sumrydf)

  out <- ggplot(sumrydf)

  out <- out +
    geom_rect(aes(xmin = area_lag, xmax = cum_area,
                  ymin = 0, ymax = water_frac, fill = classification)) +
    xlab("Cumulative Pixel Area (m^2)") + ylab("Pixel Water Fraction") +
    facet_wrap(~node_index, scales = "free_x")

  # Add gdem truth
  if (!is.null(node_truth)) {
    out <- out +
      geom_rect(aes(xmin = 0, ymin = 0, xmax = true_area, ymax = 1),
                fill = NA, color = "gray30", linetype = 2)
  }

  out
}
