# Aggregating nodes to reaches
# These operations are simplified versions of those in RiverObs, used to do on-the-fly
# recalculations of reach products using filtered, modified node-level data.


#' Weighted linear model for height as a function of distance downstream
#'
#' Used for computing reach-level heights and slopes with corresponding uncertainty.
#'
#' @param node_h Vector of node heights
#' @param node_h_u Vector of node height uncertainty (1-sigma)
#' @param node_x Vector of node downstream distance
#' @param loc_offset From reach-level rivertile, distance from observed to prior
#'  reach center
#' @param weight Use weighted least-squares? Default is TRUE.
#' @importFrom stats lm
#' @export
reach_height_lm <- function(node_h, node_h_u, node_x, loc_offset,
                            weight = TRUE) {
  weights <- 1 / (node_h_u^2)
  if (!weight) weights <- NULL

  hxdf <- data.frame(x = node_x - mean(node_x) + loc_offset, h = node_h)
  hxmod <- lm(h ~ x, hxdf, weights = weights)

  hxmod
}


#' Aggregate node-level data to reach-scale
#'
#' A slightly imprecise reproduction of RiverObs reach aggregation for a handful of
#'  variables.
#'
#' @param nodedata Must have additional columns from \code{add_nodelen()}
#' @importFrom stats median
#' @importFrom rlang .data
#' @export
reach_agg <- function(nodedata) {

  if (is.null(nodedata$cumlen) || is.null(nodedata$nodelen)) {
    stop("nodedata must have the following precomputed: nodelen, cumlen, loc_offset")
  }

  # Make linear models for height, slope
  hxmods <- split(nodedata, f = nodedata$reach_id) %>%
    purrr::map(~reach_height_lm(node_h = .$height, node_h_u = .$height2_u,
                                node_x = .$cumlen, loc_offset = .$loc_offset, weight = TRUE))
  hxcoef <- map(hxmods, ~as.data.frame(summary(.)$coefficients, row.names = FALSE)) %>%
    map(~setNames(.[, 1:2], c("est", "std"))) %>%
    map(~mutate(., param = c("intercept", "slope"))) %>%
    bind_rows(.id = "reach_id")

  #  values and uncertainties come from param estimates and standard errors
  reach_heights <- dplyr::filter(hxcoef, param == "intercept")$est
  reach_heights_u <- dplyr::filter(hxcoef, param == "intercept")$std
  reach_slopes <- dplyr::filter(hxcoef, param == "slope")$est * 1e6
  reach_slopes_u <- dplyr::filter(hxcoef, param == "slope")$std * 1e6

  nd_agg <- nodedata %>%
    group_by(.data$reach_id) %>%
    summarize(time = median(.data$time), time_tai = median(.data$time_tai),
              area = sum(.data$area_total),
              area_u = sqrt(sum(.data$area_tot_u^2)),
              width = area / sum(.data$nodelen),
              width_u = area_u / sum(.data$nodelen)) %>%
    mutate(height = reach_heights,
           height_u = reach_heights_u,
           slope = reach_slopes,
           slope_u = reach_slopes_u) %>%
    rename(area_total = area, area_tot_u = area_u)


  nd_agg
}


#' Check whether all node IDs are sequential --i.e. none are missing.
#'
#' Very simple for now, but will be updated in future as node_id format is updated.
#'
#' @param node_ids Vector of IDs that should be sequential and complete.
all_sequential <- function(node_ids) {
  ids_adj <- node_ids - min(node_ids) + 1
  isTRUE(all.equal(ids_adj, 1:length(ids_adj)))
}

#' Add loc_offset to node data
#'
#' @param nodedata,reachdata As returned by \code{rt_read()}
#' @export

add_offset <- function(nodedata, reachdata) {
  reachadjdf <- reachdata %>%
    dplyr::select(reach_id, loc_offset)
  # join to loc_offset from reachdata
  out <- nodedata %>%
    left_join(reachadjdf, by = "reach_id")

  out
}

#' Add node length and cumulative length (distance) downstream to node data
#'
#' @param nodedata As returned by \code{rt_read()}
#' @export
add_nodelen <- function(nodedata) {

  nodeids <- nodedata$node_id
  if (!all_sequential(nodeids)) stop ("Gaps exist in node data")

  out <- nodedata %>%
    arrange(node_id) %>%
    mutate(nodelen = area_total / width,
           cumlen = cumsum(nodelen))

  out
}
