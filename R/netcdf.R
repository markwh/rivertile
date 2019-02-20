# netcdf.R
# Process netcdf outputs of RiverObs
# modified from notebook20190204.Rmd.


#' Get data from a rivertile netcdf
#'
#' @param ncfile A rivertile netcdf file
#' @param group Which group to get data from: "nodes" or "reaches
#' @param keep_na_vars Keep variables that only contain missing values?
#'
#' @importFrom ncdf4 nc_open nc_close ncvar_get
#' @importFrom purrr map map_lgl
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
  if (! keep_na_vars) {
    nacols <- map_lgl(outvals_list, ~sum(!is.na(.)) == 0)
    outvals_df <- outvals_df[!nacols]
  }
  outvals_df
}

pixcvec_read <- function(ncfile, keep_na_vars = FALSE) {

  pcv_nc <- nc_open(ncfile)
  on.exit(nc_close(pcv_nc))

  pcvvars <- names(pcv_nc$var)

  outvals_list <- map(pcvvars, ~as.vector(ncvar_get(pcv_nc, .))) %>%
    setNames(pcvvars)

  outvals_df <- as.data.frame(outvals_list)
  if (! keep_na_vars) {
    nacols <- map_lgl(outvals_list, ~sum(!is.na(.)) == 0)
    outvals_df <- outvals_df[!nacols]
  }
  outvals_df
}


#' Read data from a pixel_cloud netcdf file
#'
#' @param ncfile a pixel_cloud netcdf file
#' @param group Which group to get data from: "pixel_cloud", "tvp", or "noise"
#' @param latlim
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

  outvals_list <- map(grpvars, ~as.vector(ncvar_get(pixc_nc, .))) %>%
    setNames(grpnames)

  outvals_df <- as.data.frame(outvals_list)

  if (! keep_na_vars) {
    nacols <- map_lgl(outvals_list, ~sum(!is.na(.)) == 0)
    outvals_df <- outvals_df[!nacols]
  }

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

#' Get a validation dataset from a set of RiverObs runs

rt_valdata <- function(dir, group = c("nodes", "reaches"),
                       rtname = "rt.nc", gdname = "rt_gdem.nc",
                       keep_na_vars = FALSE,
                       time_round_digits = -2) {

  group <- match.arg(group)
  rtdf <- rt_read(paste0(dir, "/", rtname), group = group,
                  keep_na_vars = keep_na_vars)
  gddf <- rt_read(paste0(dir, "/", gdname), group = group,
                  keep_na_vars = keep_na_vars)

  # ID variables for joining rivertile to gdem
  idvars <- c("reach_id", "node_id", "time", "time_tai")
  idvars <- intersect(names(rtdf), idvars)

  # time variables need to be rounded.
  timevars <- intersect(c("time", "time_tai"), names(rtdf))
  rtdf[timevars] <- round(rtdf[timevars], digits = time_round_digits)
  gddf[timevars] <- round(gddf[timevars], digits = time_round_digits)

  # variables assumed constant between rivertile and gdem, can be joined
  # separately. Or, variables only meaningful for actual rivertile data
  commonvars_rch <- c(
    "p_latitud", "p_longitud", "p_n_nodes", "xtrk_dist", "partial_f",
    "n_good_nod", "obs_frac_n", "reach_q", "geoid_height", "geoid_slop",
    "solid_tide", "pole_tide", "load_tide", "dry_trop_c", "wet_trp_c", "iono_c",
    "xover_cal_c", "p_n_nodes", "p_dist_out"
  )
  commonvars_nod <- c(
    "area_of_ht", "node_dist", "xtrk_dist", "n_good_pix", "node_q",
    "solid_tide", "pole_tide", "load_tide", "dry_trop_c", "wet_trop_c",
    "iono_c", "xover_cal_c", "p_dist_out"
  )
  commonvars <- intersect(names(rtdf), c(commonvars_rch, commonvars_nod))

  # Vector of variables to compare between rivertile and gdem
  varnames <- c("height", "height2", "slope", "width", "area_detct", "area_total",
                "latitude", "longitude")
  # Corresponding uncertainty variables
  uncnames <- setNames(
    c("height_u", "height2_u", "slope_u", "width_u", "area_det_u", "area_tot_u",
      "latitude_u", "longitud_u"),
    varnames
  )

  varnames <- intersect(names(rtdf), varnames)
  uncnames <- uncnames[varnames]

  # Make gathered data.frames
  rtdf_g <- gather(rtdf[c(idvars, varnames)],
                   key = "variable", value = "pixc_val", -!!idvars)
  gddf_g <- gather(gddf[c(idvars, varnames)],
                   key = "variable", value = "gdem_val", -!!idvars)
  uncdf_g <- rtdf[c(idvars, uncnames)] %>%
    setNames(plyr::mapvalues(names(.), from = uncnames, to = varnames)) %>%
    gather(key = "variable", value = "sigma_est", -!!idvars)

  # Join together, including "common" variables
  commondf <- rtdf[c(idvars, commonvars)]
  out <- rtdf_g %>%
    left_join(gddf_g, by = c(idvars, "variable")) %>%
    mutate(pixc_err = pixc_val - gdem_val) %>%
    left_join(uncdf_g, by = c(idvars, "variable")) %>%
    left_join(commondf, by = idvars)
  out
}

