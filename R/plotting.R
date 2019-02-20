# Plotting functions

#' Plot a histogram of errors, possibly scaled by uncertainty estimates
#'
#' Currently only implemented for node-level data
#'
#' @param dir directory containing riverobs output including gdem truth
#' @param center Subtract the mean (bias-correct) the errors?
#' @param scale Scale the errors by 1-sigma uncertainty estimates?
#' @param curve Overlay a standard normal curve? This is the default if
#'   both \code{center} and \code{scale} are \code{TRUE}.
#' @param vars Which variables to plot? Defaults to "all"
#' @param ... Passed to \code{rt_valdata()}
#' @export

rt_val_hist <- function(dir, center = FALSE, scale = TRUE,
                        curve = center && scale,
                        vars = "all",
                        ...) {
  # Currently only implemented for nodes

  valdata <- rt_valdata(dir, group = "nodes", ...)

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
#' @param dir directory containing riverobs output including gdem truth
#' @param variable Which variable to plot?
#' @param err_only Only plot the errors (subtract off truth value)?
#' @param ... Passed to \code{rt_valdata()}
#' @export

rt_val_nodeseries <- function(dir, variable = "height", err_only = TRUE, ...) {
  valdata <- rt_valdata(dir, group = "nodes", ...) %>%
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
