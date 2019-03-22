# Validation functions


#' Create validation data from two node-level rivertile data.frames.
#'
#' Used internally by \code{rt_valdata()}
#'
#' @param obs data.frame with observed node-level values
#' @param truth data.frame with node-level truth values
#' @param time_round_digits how many digits to round time (secondes) to
#'  consider equal between gdem pixc and "real" pixc
#'
#' @export
rt_valdata_df <- function(obs, truth, time_round_digits = -2) {
  # ID variables for joining rivertile to gdem
  idvars <- c("reach_id", "node_id", "time", "time_tai")
  idvars <- intersect(names(obs), idvars)

  # time variables need to be rounded.
  timevars <- intersect(c("time", "time_tai"), names(obs))
  obs[timevars] <- round(obs[timevars], digits = time_round_digits)
  truth[timevars] <- round(truth[timevars], digits = time_round_digits)

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
  commonvars <- intersect(names(obs), c(commonvars_rch, commonvars_nod))

  # Vector of variables to compare between rivertile and gdem
  varnames <- c("height", "height2", "slope", "width", "area_detct", "area_total",
                "latitude", "longitude")
  # Corresponding uncertainty variables
  uncnames <- setNames(
    c("height_u", "height2_u", "slope_u", "width_u", "area_det_u", "area_tot_u",
      "latitude_u", "longitud_u"),
    varnames
  )

  varnames <- intersect(names(obs), varnames)
  uncnames <- uncnames[varnames]

  # Make gathered data.frames
  obs_g <- gather(obs[c(idvars, varnames)],
                  key = "variable", value = "pixc_val", -!!idvars)
  truth_g <- gather(truth[c(idvars, varnames)],
                    key = "variable", value = "gdem_val", -!!idvars)
  uncdf_g <- obs[c(idvars, uncnames)] %>%
    setNames(plyr::mapvalues(names(.), from = uncnames, to = varnames)) %>%
    gather(key = "variable", value = "sigma_est", -!!idvars)

  # Join together, including "common" variables
  commondf <- obs[c(idvars, commonvars)]
  out <- obs_g %>%
    left_join(truth_g, by = c(idvars, "variable")) %>%
    dplyr::mutate(pixc_err = pixc_val - gdem_val) %>%
    left_join(uncdf_g, by = c(idvars, "variable")) %>%
    left_join(commondf, by = idvars)
  out
}

#' Get a validation dataset from a set of RiverObs runs
#'
#' @param dir directory containing riverobs output including gdem truth
#' @param group Which group to get data from: "nodes" or "reaches
#' @param rtname Name of rivertile file to read
#' @param gdname Name of gdem-truth rivertile file
#' @param keep_na_vars Keep variables that only contain missing values?
#' @param time_round_digits how many digits to round time (secondes) to
#'  consider equal between gdem pixc and "real" pixc
#' @importFrom dplyr left_join
#' @importFrom tidyr gather
#' @export
rt_valdata <- function(dir, group = c("nodes", "reaches"),
                       rtname = "rt.nc", gdname = "rt_gdem.nc",
                       keep_na_vars = FALSE,
                       time_round_digits = -2,
                       flag_out_nodes = TRUE) {

  group <- match.arg(group)
  rtdf <- rt_read(paste0(dir, "/", rtname), group = group,
                  keep_na_vars = keep_na_vars)
  gddf <- rt_read(paste0(dir, "/", gdname), group = group,
                  keep_na_vars = keep_na_vars)

  out <- rt_valdata_df(obs = rtdf, truth = gddf, time_round_digits = -2)

  if (group == "nodes" && flag_out_nodes) {
    rmnodes <- flag_nodes(dir)
    out <- out[!out[["node_id"]] %in% rmnodes, ]
  }

  out
}

#' Flag nodes with questionable / unstable gdem truth
#'
#' Based on differences in "dilated" and non-dilated gdem.
#'
#' @param dir Directory containing gdem rivertile netcdfs
#' @param gdem1,gdem2 gdem rivertile netcdfs to compare
#' @param thresh Threshold for relative difference in widths, defaults to 0.15.
#'
#' @export
flag_nodes <- function(dir, gdem1 = "rt_gdem.nc", gdem2 = "rt_gdem_dil2.nc",
                       thresh = 0.15) {
  valdata <- rt_valdata(dir, group = "nodes", rtname = gdem1, gdname = gdem2,
                        keep_na_vars = TRUE, flag_out_nodes = FALSE) %>%
    dplyr::select(reach_id, node_id, variable, pixc_val, gdem_val, pixc_err) %>%
    dplyr::filter(variable == "width") %>%
    dplyr::mutate(pct_diff = - pixc_err / gdem_val)

  out <- with(valdata, node_id[pct_diff > thresh])
  out
}


#' Table of coverage for various (normal) confidence intervals.
#'
#' @param valdata as returned by \code{rt_valdata()}
#' @param ci Which confidence levels to validate?
#'
#' @export
val_coverage <- function(valdata, ci = c(68, 90, 95, 99), debias = FALSE) {
  covfun <- function(x, sigma, pctl) {
    if (debias) x <- x - mean(x, na.rm = TRUE)
    pctl <- pctl / 100
    bnd <- -qnorm((1 - pctl) / 2, mean = 0, sd = 1)
    numin <- vapply(bnd, function(y) sum(abs(x / sigma) <= y), numeric(1))
    out <- numin / length(x) * 100
    out
  }

  vallist <- split(valdata, f = valdata$variable)
  out <- purrr::map(vallist, ~covfun(x = .$pixc_err, sigma = .$sigma_est, pctl = ci)) %>%
    as.data.frame() %>%
    t() %>%
    as.data.frame() %>%
    setNames(paste0("ci", ci))

  out
}

