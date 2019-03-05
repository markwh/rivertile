# Hypothesis testing

#' Test uncertainty estimates against chi-square distribution
#'
#' @param valdata A data.frame as returned by \code{rt_valdata()}
#' @param debias Remove effect of bias? Defaults to TRUE
#' @importFrom dplyr group_by mutate summarize
#' @export
rt_hyptest <- function(valdata, debias = TRUE) {
  out <- valdata %>%
    mutate(relerr = pixc_err / sigma_est) %>%
    group_by(variable) %>%
    mutate(meanrelerr = mean(relerr, na.rm = TRUE))

  if (!debias) out[["meanrelerr"]] <- 0

  out <- summarize(out,
                   teststat = sum((relerr - meanrelerr)^2, na.rm = TRUE),
                   df = sum(!is.na(relerr)) - 1) %>%
    mutate(pval = 2 * (1 - pchisq(teststat, df = df)))
  out
}

