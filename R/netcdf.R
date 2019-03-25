# netcdf.R
# Process netcdf outputs of RiverObs
# modified from notebook20190204.Rmd, notebook20190318.Rmd, others?


#' Read in only a subset of a netcdf variable
#'
#' @param nc netcdf with variable of interest
#' @param varid passed to \code{ncvar_get}
#' @param inds Vector or matrix giving indices for which the variable value is desired.
#'
#' @importFrom stats df dnorm pchisq qnorm setNames
#' @export
ncvar_ss <- function(nc, varid=NA, inds) {
  if (length(dim(inds)) > 2) stop("dimensionality > 2 not supported.")
  if (length(dim(inds)) == 0) {
    return(ncvar_ss_1d(nc, varid, inds))
  }
  indlist <- split(inds[, 1], f = inds[, 2])

  pb <- progress::progress_bar$new(total = length(indlist))
  pb$tick(0)

  vals <- list()
  for (i in 1:length(indlist)) {
    pb$tick()
    vals[[i]] <- ncvar_ss_1d(nc, varid = varid, indvec = indlist[[i]],
                             inds2 = as.numeric(names(indlist)[i]))
  }
  out <- unlist(vals)
  out
}

ncvar_ss_1d <- function(nc, varid=NA, indvec, inds2 = NULL) {
  stopifnot(is.numeric(indvec) && is.vector(indvec))

  minind <- min(indvec)
  maxind <- max(indvec)
  indcnt <- maxind - minind + 1
  newinds <- indvec - minind + 1

  if (!is.null(inds2)) {
    minind <- cbind(minind, inds2)
    indcnt <- cbind(indcnt, 1)
  }

  out <- ncvar_get(nc, varid = varid, start = minind, count = indcnt)

  out <- as.vector(out[newinds])
  out
}

#' Check for inefficiencies in ncdf subset read function.
#'
#' Result should be close to 1.
#'
#' @param inds Matrix giving rows and columns of indices
#' @export
check_ineff <- function(inds) {
  num1 <- nrow(inds)
  num2 <- as.data.frame(inds) %>%
    group_by(col) %>%
    summarize(nret = max(row) - min(row) + 1) %>%
    summarize(sum = sum(nret)) %>%
    `[[`("sum")
  out <- num2 / num1
  out
}

adjust_longitude <- function(df) {
  lonname_cands <- c("longitude", "longitude_vectorproc", "p_longitud")
  lonnames <- intersect(lonname_cands, names(df))

  for (lonname in lonnames) {
    lonvec <- df[[lonname]]
    df[[lonname]][lonvec > 180] <- lonvec[lonvec > 180] - 360
  }
  df
}

#' Get data from a rivertile netcdf
#'
#' @param ncfile A rivertile netcdf file
#' @param group Which group to get data from: "nodes" or "reaches
#' @param keep_na_vars Keep variables that only contain missing values?
#'
#' @importFrom ncdf4 nc_open nc_close ncvar_get
#' @importFrom purrr map map_lgl
#' @importFrom dplyr "%>%"
#' @export

rt_read <- function(ncfile, group = c("nodes", "reaches"),
                    keep_na_vars = FALSE) {
  group <- match.arg(group)

  rt_nc <- nc_open(ncfile)
  on.exit(nc_close(rt_nc))

  grepstr <- sprintf("^%s/", group)

  grpvars <- names(rt_nc$var)[grepl(grepstr, names(rt_nc$var))]
  grpnames <- splitPiece(grpvars, "/", 2, fixed = TRUE)

  outvals_list <- map(grpvars, ~as.vector(ncvar_get(rt_nc, .))) %>%
    setNames(grpnames)

  outvals_df <- as.data.frame(outvals_list)

  # Comply with -180:180 convention used by RiverObs
  outvals_df <- adjust_longitude(outvals_df)

  if (! keep_na_vars) {
    nacols <- map_lgl(outvals_list, ~sum(!is.na(.)) == 0)
    outvals_df <- outvals_df[!nacols]
  }
  outvals_df
}


#' Read a pixcvec from a netcdf file
#'
#' @param ncfile PIXCVEC netcdf file
#' @param keep_na_vars Keep variables that only contain missing values?
#' @export
pixcvec_read <- function(ncfile, keep_na_vars = FALSE) {

  pcv_nc <- nc_open(ncfile)
  on.exit(nc_close(pcv_nc))

  pcvvars <- names(pcv_nc$var)

  outvals_list <- map(pcvvars, ~as.vector(ncvar_get(pcv_nc, .))) %>%
    setNames(pcvvars)

  outvals_df <- as.data.frame(outvals_list)
  # Comply with -180:180 convention used by RiverObs
  outvals_df <- adjust_longitude(outvals_df)

  if (! keep_na_vars) {
    nacols <- map_lgl(outvals_list, ~sum(!is.na(.)) == 0)
    outvals_df <- outvals_df[!nacols]
  }
  outvals_df
}

#' Read a gdem netcdf
#'
#' @param ncfile gdem netcdf file
#'
#' @export
gdem_read <- function(ncfile) {
  pixc_nc <- nc_open(ncfile)
  on.exit(nc_close(pixc_nc))

  ltype <- ncvar_get(pixc_nc, "landtype")
  waterpix <- which(!is.na(ltype) & ltype == 1, arr.ind = TRUE)

  lats <- ncvar_ss(pixc_nc, "latitude", inds = waterpix)
  lons <- ncvar_ss(pixc_nc, "longitude", inds = waterpix)

  out <- data.frame(latitude = lats, longitude = lons)
  out <- adjust_longitude(out)
  out
}

#' Read data from a pixel_cloud netcdf file
#'
#' @param ncfile a pixel_cloud netcdf file
#' @param group Which group to get data from: "pixel_cloud", "tvp", or "noise"
#' @param latlim,lonlim Limits of latitude and longitude, as length-2 vectors
#' @param keep_na_vars Keep variables with nothing but missing values?
#'
#' @export
pixc_read <- function(ncfile, group = c("pixel_cloud", "tvp", "noise"),
                      latlim = c(-90, 90), lonlim = c(-180, 180),
                      keep_na_vars = FALSE) {
  # browser()
  group <- match.arg(group)

  pixc_nc <- nc_open(ncfile)
  on.exit(nc_close(pixc_nc))

  grepstr <- sprintf("^%s/", group)

  grpvars <- names(pixc_nc$var)[grepl(grepstr, names(pixc_nc$var))]
  grpnames <- splitPiece(grpvars, "/", 2, fixed = TRUE)

  outvals_list <- map(grpvars, ~ncvar_get(pixc_nc, .)) %>%
    setNames(grpnames)
  # split interferogram array into real and imaginary vectors
  outvals_list$interferogram_r <- outvals_list$interferogram[1, ]
  outvals_list$interferogram_i <- outvals_list$interferogram[2, ]
  outvals_list$interferogram <- NULL
  outvals_list <- purrr::map(outvals_list, ~as.vector(.))

  outvals_df <- as.data.frame(outvals_list)

  if (! keep_na_vars) {
    nacols <- map_lgl(outvals_list, ~sum(!is.na(.)) == 0)
    outvals_df <- outvals_df[!nacols]
  }

  # Comply with -180:180 convention used by RiverObs
  outvals_df <- adjust_longitude(outvals_df)

  # Filter lat/lon that are exactly 0, restrict to bounding box
  if (group == "pixel_cloud") {
    outvals_df <- dplyr::filter(outvals_df,
                    !is.na(latitude), !is.na(longitude),
                    latitude != 0, longitude != 0,
                    longitude >= lonlim[1], longitude <= lonlim[2],
                    latitude >= latlim[1], latitude <= latlim[2])
  }

  outvals_df
}


priordb_read <- function(ncfile, group = c("reaches", "nodes", "centerlines"),
                         latlim = c(-90, 90), lonlim = c(-180, 180),
                         keep_na_vars = FALSE) {

  group <- match.arg(group)

  pdb_nc <- nc_open(ncfile)
  on.exit(nc_close(pdb_nc))

  grepstr <- sprintf("^%s/", group)

  grpvars <- names(pdb_nc$var)[grepl(grepstr, names(pdb_nc$var))]
  grpnames <- splitPiece(grpvars, "/", 2, fixed = TRUE)

  outvals_list <- map(grpvars, ~as.vector(ncvar_get(pdb_nc, .))) %>%
    setNames(grpnames)

  outvals_df <- as.data.frame(outvals_list) %>%
    dplyr::mutate(x = ifelse(x > 180, x - 360, x),
                  lon = x, lat = y) %>%
    dplyr::filter(x >= lonlim[1], x <= lonlim[2],
                  y >= latlim[1], y <= latlim[2])

  if (! keep_na_vars) {
    nacols <- map_lgl(outvals_list, ~sum(!is.na(.)) == 0)
    outvals_df <- outvals_df[!nacols]
  }
  outvals_df
}

