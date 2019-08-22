# Functions to read SWOT river products
# Migrated from swot-error/lib/rivertile.R



#' Read prior database for node info
#'
#' @param ncfile netcdf file containing prior info
#' @param nodeids optional vector of node indices
#' @param reachids optional vector of reach indices
#' @param as_sf Convert to spatial frame?
#'
#'
#'
#' @export
priornode_read <- function(ncfile, nodeids = NULL, reachids = NULL, as_sf = TRUE) {
  nc <- nc_open(ncfile)
  on.exit(nc_close(nc))

  getvar <- function(var, ...) as.vector(ncvar_get(nc, var, ...))
  ncnodeids <- getvar("nodes/node_id")
  ncreachids <- getvar("nodes/reach_id")

  ncinds <- 1:length(ncreachids)
  if (!is.null(nodeids)) {
    ncinds <- match(nodeids, which(ncreachids %in% nodeids))
  }
  if (!is.null(reachids)) {
    ncinds <- intersect(ncinds, which(ncreachids %in% reachids))
  }
  # browser()
  outinds <- ncinds - min(ncinds) + 1
  readstart <- min(ncinds)
  readlen <- max(ncinds) - min(ncinds) + 1

  node_df <- data.frame(
    node_id = ncnodeids[outinds],
    reach_id = ncreachids[outinds],
    latitude = getvar("nodes/y", start = readstart, count = readlen)[outinds],
    longitude = getvar("nodes/x", start = readstart, count = readlen)[outinds])

  if (as_sf) {
    if (!check_sf()) stop("as_sf argument requires sf package. Please install it.")
    out <- sf::st_as_sf(node_df, coords = c("longitude", "latitude"),
                    crs = "+proj=longlat +datum=WGS84")

  } else {
    out <- node_df
  }

  out
}


#' Read prior database for node info
#'
#' Unclear when this function would be useful.
#'
#' @param reachids vector of reach indices
#' @param ncfile netcdf file containing prior info
#'
priorreach_read <- function(ncfile, reachids = NULL) {
  nc <- nc_open(ncfile)
  on.exit(nc_close(nc))

  getvar <- function(var, ...) as.vector(ncvar_get(nc, var, ...))
  ncreachids <- getvar("reaches/reach_id")

  ncinds <- 1:length(ncreachids)
  if (!is.null(reachids)) {
    ncinds <- which(ncreachids %in% reachids)
  }

  outinds <- ncinds - min(ncinds) + 1
  readstart <- min(ncinds)
  readlen <- max(ncinds) - min(ncinds) + 1

  out <- data.frame(
    reach_id = getvar("reaches/reach_id", start = readstart, count = readlen)[outinds],
    latitude = getvar("reaches/y", start = readstart, count = readlen)[outinds],
    longitude = getvar("reaches/x", start = readstart, count = readlen)[outinds])
  out
}


#' Read prior database for centerline info
#'
#' Returns a sf object with linestring geometry.
#'
#' @param ncfile netcdf file containing prior info
#' @param nodeids optional vector of node indices
#' @param reachids optional vector of reach indices
#' @param as_sf Convert to spatial frame object?
#'
#' @importFrom dplyr group_by summarize ungroup
#'
#' @export
priorcl_read <- function(ncfile, nodeids = NULL, reachids = NULL,
                         as_sf = TRUE) {
  nc <- nc_open(ncfile)
  on.exit(nc_close(nc))

  getvar <- function(var, ...) as.vector(ncvar_get(nc, var, ...))
  ncnodeids <- ncvar_get(nc, "centerlines/node_id")[, 1]
  ncreachids <- ncvar_get(nc, "centerlines/reach_id")[, 1]

  ncinds <- 1:length(ncnodeids)
  if (!is.null(nodeids)) {
    ncinds <- match(nodeids, which(ncreachids %in% nodeids))
  }
  if (!is.null(reachids)) {
    ncinds <- intersect(ncinds, which(ncreachids %in% reachids))
  }

  outinds <- ncinds - min(ncinds) + 1
  readstart <- min(ncinds)
  readlen <- max(ncinds) - min(ncinds) + 1

  cl_df <- data.frame(
    node_id = ncnodeids[outinds],
    reach_id = ncreachids[outinds],
    latitude = getvar("centerlines/y", start = readstart, count = readlen)[outinds],
    longitude = getvar("centerlines/x", start = readstart, count = readlen)[outinds])

  if (as_sf) {
    if (!check_sf()) stop("as_sf argument requires sf package. Please install it.")
    cl_points <- sf::st_as_sf(cl_df, coords = c("longitude", "latitude"),
                          crs = "+proj=longlat +datum=WGS84")
    out <- cl_points %>%
      group_by(node_id, reach_id) %>%
      summarize(geometry = sf::st_cast(geometry, to = "LINESTRING", ids = 1)) %>%
      ungroup()

  } else {
    out <- cl_df
  }

  out
}

#' Read orbit locations
#'
#' Returns a sf object with multilinestring geometry.
#'
#' @param ncfile netcdf file containing prior info
#' @param as_sf Convert to spatial frame object?
#' @param maxpoints Maximum number of points to use for line resolution.
#'
#' @importFrom dplyr mutate
#'
#' @export
orbit_read <- function(ncfile, as_sf = TRUE, maxpoints = 1000) {
  nc <- nc_open(ncfile)
  on.exit(nc_close(nc))

  getvar <- function(var, ...) as.vector(ncvar_get(nc, var, ...))
  passlat <- getvar("latitude")
  passlon <- getvar("longitude")
  npts <- length(passlat)

  keepinds <- 1:npts
  if (npts > maxpoints) {
    keepinds <- seq(1, npts, length.out = maxpoints)
  }

  out <- data.frame(latitude = passlat, longitude = passlon) %>%
    mutate(londif = c(0, diff(longitude)),
           cutoff = abs(londif) > 300,
           splitvar = cumsum(cutoff)) %>%
    `[`(keepinds, )

  if (as_sf) {
    if (!check_sf()) stop("as_sf argument requires sf package. Please install it.")
    out <- out %>%
      sf::st_as_sf(coords = c("longitude", "latitude"),
               crs = "+proj=longlat +datum=WGS84") %>%
      group_by(splitvar) %>%
      summarize(geometry = sf::st_cast(geometry, to = "LINESTRING", ids = splitvar)) %>%
      summarize(geometry = sf::st_cast(geometry, to = "MULTILINESTRING", ids = 1))
  }

  out
}

#' Get corners of tile polygon
#'
#' @param nadir1,nadir2 length-2 vector giving longitude and latitude
#' @param heading satellite heading in degrees
#' @param xtstart,xtend Distance from nadir to start/end the swath, meters
#' @param half Which half of the tile: "L" for left, "R" for right
getTileCorners <- function(nadir1, nadir2, heading, xtstart = 4000,
                           xtend = 64000, half = c("L", "R")) {
  if (!check_geosphere()) stop("getTileCorners requires geosphere package. Please install it.")
  half <- match.arg(half)

  xtdir <- heading + ifelse(half == "R", 90, -90)
  xtdir <- ifelse(abs(xtdir) > 180, xtdir - 360, xtdir)


  points1 <- geosphere::destPoint(nadir1, b = xtdir, d = c(xtstart, xtend))
  points2 <- geosphere::destPoint(nadir2, b = xtdir, d = c(xtstart, xtend))

  out <- rbind(points1, points2[2:1, ], points1[1, ])
  out
}

# corner2sf <- function(cornermat) {
#   out <- sf::st_polygon(cornermat)
#   out
# }

#' Get tiles as spatial frames POLYGON sfc
#'
#' Returns a sfc with tile polygons
#'
#' @param nadir1,nadir2 Either a length-2 vector or a 2-column matrix giving lon, lat
#' @inheritParams getTileCorners
getTilePolygons <- function(nadir1, nadir2, heading, half, xtstart = 4000,
                            xtend = 64000) {
  if (!check_sf()) stop("getTilePolygons requires sf package. Please install it.")
  splitfun <- function(x) {
    if (is.numeric(x) && is.vector(x)) {
      stopifnot(length(x) == 2)
      x <- matrix(x, ncol = 2)
    }
    out <- split(x, f = 1:nrow(x))
  }
  nadir1 <- splitfun(nadir1)
  nadir2 <- splitfun(nadir2)
  stopifnot(length(nadir1) == length(nadir2))

  inputlist <- list(nadir1 = nadir1, nadir2 = nadir2, heading = heading, half = half)
  # browser()
  cornermats <- purrr::pmap(inputlist, getTileCorners, xtstart = xtstart,
                            xtend = xtend)
  cornerpolys <- purrr::map(cornermats, ~sf::st_polygon(list(.)))
  out <- sf::st_sfc(cornerpolys, crs = "+proj=longlat +datum=WGS84")
  out
}



