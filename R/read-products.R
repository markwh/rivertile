# Functions for reading prior database
# Migrated from netcdf.R
# Process netcdf outputs of RiverObs
# modified from notebook20190204.Rmd, notebook20190318.Rmd, others?

#' Helper function to coerce attribute list to tabular
#'
#' Current approach is to paste everything together--in future
#' I may use list-valued columns.
#'
#' @param attlist as returnd by \code{ncdf4::ncatt_get()}
atts_df_fun <- function(attlist) {
  out <- attlist
  lens <- purrr::map_int(attlist, length)
  out[lens > 1] <- lapply(out[lens > 1], paste, collapse = ", ")
  dplyr::as_tibble(out)
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
#' @param as_sf Convert to a spatial frame?
#'
#' @importFrom ncdf4 nc_open nc_close ncvar_get ncatt_get
#' @importFrom purrr map map_lgl
#' @importFrom dplyr "%>%"
#' @export

rt_read <- function(ncfile, group = c("nodes", "reaches"),
                    keep_na_vars = FALSE, as_sf = FALSE) {
  group <- match.arg(group)

  rt_nc <- nc_open(ncfile)
  on.exit(nc_close(rt_nc))

  grepstr <- sprintf("^%s/", group)
  grpvars <- names(rt_nc$var)[grepl(grepstr, names(rt_nc$var))]
  grpnames <- splitPiece(grpvars, "/", 2, fixed = TRUE)

  if (!as_sf) {
    ctlvar <- grepl("centerline", grpvars)
    grpvars <- grpvars[!ctlvar]
    grpnames <- grpnames[!ctlvar]
  } else {
    stop("as_sf not yet implemented.")
  }

  # Smartly vectorize arrays. If already vector-like, use as.vector.
  # Otherwise make a data.frame.
  vecfun <- function(x) {
    if (length(dim(x)) < 2) return(as.vector(x))
    if (length(dim(x)) == 2) return(as.data.frame(t(x)))
    stop("Too many dimensions")
  }
  # browser()
  outvals_list <- map(grpvars, ~vecfun(ncvar_get(rt_nc, .))) %>%
    setNames(grpnames)
  ncol_list <- purrr::map_int(outvals_list, length)
  outatts_list <- map(grpvars, ~vecfun(ncatt_get(rt_nc, .))) %>%
    setNames(grpnames) %>%
    map(atts_df_fun)
  # browser()

  outvals_df <- tidyr::unnest(as.data.frame(outvals_list))
  outatts_df <- bind_rows2(outatts_list) %>%
    mutate(name = grpnames)

  # Comply with -180:180 convention used by RiverObs
  outvals_df <- adjust_longitude(outvals_df)

  if (! keep_na_vars) {
    nacols <- map_lgl(outvals_df, ~sum(!is.na(.)) == 0)
    outvals_df <- outvals_df[!nacols]
    # outatts_df <- outatts_df[!nacols, ]
  }
  outvals_df
  out <- structure(outvals_df, atts = outatts_df)
  out
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

  outatts_list <- map(pcvvars, ~ncatt_get(pcv_nc, .)) %>%
    setNames(pcvvars) %>%
    map(atts_df_fun)
  # browser()

  outvals_df <- as.data.frame(outvals_list)
  outatts_df <- bind_rows2(outatts_list) %>%
    mutate(name = pcvvars)

  # Comply with -180:180 convention used by RiverObs
  outvals_df <- adjust_longitude(outvals_df)

  if (! keep_na_vars) {
    nacols <- map_lgl(outvals_list, ~sum(!is.na(.)) == 0)
    outvals_df <- outvals_df[!nacols]
    outatts_df <- outatts_df[!nacols, ]
  }
  out <- structure(outvals_df, atts = outatts_df)
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
  group <- match.arg(group)

  pixc_nc <- nc_open(ncfile)
  on.exit(nc_close(pixc_nc))

  grepstr <- sprintf("^%s/", group)

  grpvars <- names(pixc_nc$var)[grepl(grepstr, names(pixc_nc$var))]
  grpnames <- splitPiece(grpvars, "/", 2, fixed = TRUE)

  outvals_list <- map(grpvars, ~ncvar_get(pixc_nc, .)) %>%
    setNames(grpnames)
  outatts_list <- map(grpvars, ~ncatt_get(pixc_nc, .)) %>%
    setNames(grpnames)

  # split interferogram array into real and imaginary vectors
  if (!is.null(outvals_list$interferogram)) {
    outvals_list$interferogram_r <- outvals_list$interferogram[1, ]
    outvals_list$interferogram_i <- outvals_list$interferogram[2, ]
    outvals_list$interferogram <- NULL

    # Replace interferogram parts of attribute table
    outatts_list$interferogram_r <- outatts_list$interferogram
    outatts_list$interferogram_i <- outatts_list$interferogram
    outatts_list$interferogram <- NULL
  }

  # Convert to data.frames
  outvals_list <- purrr::map(outvals_list, ~as.vector(.))
  outvals_df <- as.data.frame(outvals_list)

  outatts_df <- attslist_to_dataframe(outatts_list)

  if (! keep_na_vars) {
    nacols <- map_lgl(outvals_list, ~sum(!is.na(.)) == 0)
    outvals_df <- outvals_df[!nacols]
    if (!is.null(outatts_df)) outatts_df <- outatts_df[!nacols, ]
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

  # reclass variables
  outvals_df$classification <- as.factor(outvals_df$classification)

  out <- structure(outvals_df, atts = outatts_df)
  out
}

#' Join pixcvec to pixel cloud
#'
#' @param pixcdf as returned by \code{pixc_read()}
#' @param pcvdf as returned by \code{pixcvec_read()}
#'
#' @param type What kind of join to perform? Default is "inner"
#' @importFrom fs path
#' @importFrom dplyr inner_join
#' @export
#'
pixc_join <- function(pixcdf, pcvdf, type = c("inner", "outer")) {

  type <- match.arg(type)
  joinfun <- if (type == "inner") dplyr::inner_join else dplyr::full_join

  outdf <- pixcdf %>%
    joinfun(pcvdf, by = c("azimuth_index", "range_index")) %>%
    dplyr::mutate(pixel_id = 1:nrow(.))
  outattdf <- bind_rows2(list(attr(pcvdf, "atts"),
                              attr(pixcdf, "atts")),
                         addMissing = TRUE)

  out <- structure(outdf, atts = outattdf)
}


#' Convert pixc attributes to data frame.
#'
#' @importFrom dplyr mutate
attslist_to_dataframe <- function(attslist) {
  if (sum(purrr::map_int(attslist, length)) == 0) return(NULL)
  smplfunc <- function(x) {
    needscollapse <- lengths(x) > 1
    x[needscollapse] <- lapply(x[needscollapse], paste, collapse = ", ")
    x
  }
  veclist <- lapply(attslist, smplfunc)

  out <- bind_rows2(veclist) %>%
    mutate(name = names(attslist))
  out
}
