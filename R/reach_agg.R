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
#' @param x_ctr center of reach, as given by meandist in \code{add_offset()}
#' @param loc_offset From reach-level rivertile, distance from observed to prior
#'  reach center
#' @param weight Use weighted least-squares? Default is TRUE.
#' @importFrom stats lm
#' @export
reach_height_lm <- function(node_h, node_h_u, node_x, x_ctr, loc_offset,
                            weight = TRUE) {
  # browser()
  weights <- 1 / (node_h_u^2)
  if (!weight) weights <- NULL

  hxdf <- data.frame(x = node_x - x_ctr + loc_offset, h = node_h)
  hxmod <- lm(h ~ x, hxdf, weights = weights)

  hxmod
}


#' Aggregate node-level data to reach-scale
#'
#' A slightly imprecise reproduction of RiverObs reach aggregation for a handful of
#'  variables.
#'
#' @param nodedata Must have additional columns from \code{add_nodelen(), add_offset()}
#' @param usenodes Which nodes to use. Defaults to "all"
#' @param weight Use weighted regression based on height uncertainty?
#' @importFrom stats median
#' @importFrom rlang .data
#' @importFrom dplyr bind_rows rename
#' @export
reach_agg <- function(nodedata, usenodes = "all", weight = TRUE) {

  if (is.null(nodedata$cumlen) || is.null(nodedata$nodelen)) {
    stop("nodedata must have the following precomputed: nodelen, cumlen, loc_offset")
  }

  # Filter to specified nodes.
  if (usenodes != "all") {
    nodedata <- dplyr::filter(nodedata, node_id %in% usenodes)
  }

  # Make sure sufficient number of nodes are present per reach
  usereaches <- reaches_tokeep(nodedata, minnodes = 3)
  warn_reaches(nodedata, usereaches)
  nodedata <- dplyr::filter(nodedata, reach_id %in% usereaches)
  if (nrow(nodedata) == 0) return(NULL)

  # Make linear models for height, slope
  hxmods <- split(nodedata, f = nodedata$reach_id) %>%
    purrr::map(~reach_height_lm(node_h = .$wse, node_h_u = .$wse_r_u,
                                node_x = .$cumlen,
                                x_ctr = .$meandist,
                                loc_offset = .$loc_offset,
                                weight = weight))
  hxcoef <- map(hxmods, ~as.data.frame(summary(.)$coefficients,
                                       row.names = FALSE)) %>%
    map(~setNames(.[, 1:2], c("est", "std"))) %>%
    map(~mutate(., param = c("intercept", "slope"))) %>%
    bind_rows(.id = "reach_id")

  #  values and uncertainties come from param estimates and standard errors
  reach_heights <- dplyr::filter(hxcoef, param == "intercept")$est
  reach_heights_u <- dplyr::filter(hxcoef, param == "intercept")$std
  reach_slopes <- -dplyr::filter(hxcoef, param == "slope")$est
  reach_slopes_u <- dplyr::filter(hxcoef, param == "slope")$std

  nd_agg <- nodedata %>%
    group_by(.data$reach_id) %>%
    summarize(time = median(.data$time), time_tai = median(.data$time_tai),
              area = sum(.data$area_total),
              area_u = sqrt(sum(.data$area_tot_u^2)),
              width = area / sum(.data$nodelen),
              width_u = area_u / sum(.data$nodelen)) %>%
    mutate(wse = reach_heights,
           wse_r_u = reach_heights_u,
           slope = reach_slopes,
           slope_r_u = reach_slopes_u) %>%
    rename(area_total = area, area_tot_u = area_u)


  nd_agg
}

reaches_tokeep <- function(nodedata, minnodes = 3) {
  if (nrow(nodedata) == 0) {
    return(numeric(0))
  }
  nperdf <- nodedata %>%
    group_by(reach_id) %>%
    summarize(n = n())
  out <- nperdf$reach_id[nperdf$n > minnodes]
  out
}

warn_reaches <- function(nodedata, usereaches) {
  badreaches <- setdiff(nodedata$reach_id, usereaches)
  if (length(badreaches))
    warning("The following reaches had insufficient data and were ",
            "removed prior to aggregating: ",
            badreaches)
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

  # check for node distance
  if (is.null(nodedata$cumlen)) stop ("nodedata must have length data (add_nodelen())")

  nnodedf <- reachdata %>%
    dplyr::select(reach_id, n_good_nod)

  meandistdf <- nodedata %>%
    group_by(reach_id) %>%
    summarize(n = n(), meandist = mean(cumlen)) %>%
    left_join(nnodedf, by = "reach_id")

  # Check number of nodes is same used in riverobs
  badnreaches <- meandistdf$n != meandistdf$n_good_nod
  if (sum(badnreaches) > 0) {
    warning("The following reaches have the wrong number of nodes:\n",
            paste(meandistdf$reach_id[badnreaches], collapse = ", "), "\n",
            paste(meandistdf$n[badnreaches], collapse = ", "),
            " Should be: \n",
            paste(meandistdf$n_good_nod[badnreaches]),
            "\nDistances may be wrong.")
  }

  # join to loc_offset from reachdata
  out <- nodedata %>%
    left_join(reachadjdf, by = "reach_id") %>%
    left_join(meandistdf[c("reach_id", "meandist")], by = "reach_id")

  out
}

#' Add node length and cumulative length (distance) downstream to node data
#'
#' @param nodedata As returned by \code{rt_read()}
#' @param force Logical, force via estimate if gaps exist?
#'
#' @importFrom dplyr arrange
#' @export
add_nodelen <- function(nodedata, force = FALSE) {

  nodeids <- nodedata$node_id
  if (!all_sequential(nodeids) && ! force) stop ("Gaps exist in node data")

  out <- nodedata %>%
    arrange(node_id) %>%
    group_by(reach_id) %>%
    mutate(nodelen = area_total / width,
           cumlen = cumsum(nodelen),
           meanlen = mean(nodelen)) %>%
    ungroup()
  if (force) { # attempt to estimate cumlen.
    allnodes <- seq(min(nodeids), max(nodeids))
    missingnodes <- sort(setdiff(allnodes, out$node_id))

    meandist <- mean(out$nodelen)

    for (i in seq_along(missingnodes)) {
      nodeidi <- missingnodes[i]
      closestnodeind <- which.min(abs(nodeidi - out$node_id)) # TODO: fix this hack
      meanleni <- out[closestnodeind, ]$meanlen

      gtidinds <- which(out$node_id > nodeidi)

      out$cumlen[gtidinds] <- out$cumlen[gtidinds] + meanleni
    }
  }
  out$meanlen <- NULL

  out
}
