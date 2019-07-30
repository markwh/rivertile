context("reading riverproduct files")

test_that("pixel cloud returns data.frame with netcdf attributes", {
  # cat(getwd())
  # print(getwd())
  expect_is(pixc1 <- pixc_read("pixel_cloud.nc"),
            "data.frame")
  expect_is(pixc1.5 <- pixc_read("pixel_cloud.nc", "pixel_cloud"),
            "data.frame")
  expect_is(pixc2 <- pixc_read("pixel_cloud.nc", group = "tvp"),
            "data.frame")
  expect_is(pixc3 <- pixc_read("pixel_cloud.nc", group = "noise"),
            "data.frame")

  expect_identical(pixc1, pixc1.5)

  expect_is(atts1 <- attr(pixc1, "atts"), "data.frame")
  expect_is(atts2 <- attr(pixc2, "atts"), "data.frame")
  expect_is(atts3 <- attr(pixc3, "atts"), "data.frame")

  expect_equivalent(names(atts1), names(atts2))
  expect_equivalent(names(atts2), names(atts3))
  expect_false(nrow(atts1) == nrow(atts2))
  expect_false(nrow(atts3) == nrow(atts2))
  expect_false(nrow(atts1) == nrow(atts3))

})


test_that("pixcvec returns data.frame with netcdf attributes", {
  expect_is(pcv1 <- pixcvec_read("pcv.nc"), "data.frame")
  expect_is(atts1 <- attr(pcv1, "atts"), "data.frame")

})


test_that("rivertile returns data.frame with netcdf attributes", {
  expect_is(rt1 <- rt_read("rt.nc"), "data.frame")
  expect_is(rt1.1 <- rt_read("rt.nc", group = "nodes"), "data.frame")
  expect_is(rt1.2 <- rt_read("rt.nc", group = "nodes", keep_na_vars = TRUE),
            "data.frame")
  # expect_is(rt1.3 <- rt_read("rt.nc", group = "nodes", as_sf = TRUE),
  #           "sf")


  expect_is(rt2 <- rt_read("rt.nc", group = "reaches"), "data.frame")
  expect_is(rt2.1 <- rt_read("rt.nc", group = "reaches", keep_na_vars = TRUE),
            "data.frame")
  # expect_is(rt2.2 <- rt_read("rt.nc", group = "reaches", as_sf = TRUE),
  #           "sf")


  expect_equivalent(rt1, rt1.1)

  expect_is(atts1 <- attr(rt1, "atts"), "data.frame")
  expect_is(atts2 <- attr(rt2, "atts"), "data.frame")

  expect_lt(ncol(rt1), ncol(rt1.2))
  expect_lt(ncol(rt2), ncol(rt2.1))


  expect_equivalent(names(atts1), names(atts2))
  expect_false(nrow(atts1) == nrow(atts2))

})
