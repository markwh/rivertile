# Preparing data for tests and extdata.

# Start with full-size netcdfs.
library(ncdf4)
library(purrr)
library(subsetnc)



# pixel cloud, pixcvec
# TODO: replace with url
raw_pixcfile <- "D:/data/riverobs-output/sacruns_20190709/90/pixel_cloud.nc"
raw_pcvfile <- "D:/data/riverobs-output/sacruns_20190709/90/pcv.nc"
raw_rtfile <- "D:/data/riverobs-output/sacruns_20190709/90/rt.nc"

test_pixcfile <- "tests/testthat/pixel_cloud.nc"
test_pcvfile <- "tests/testthat/pcv.nc"
test_rtfile <- "tests/testthat/rt.nc"

ext_pixcfile <- "extdata/pixel_cloud.nc"
ext_pcvfile <- "extdata/pcv.nc"
ext_rtfile <- "extdata/rt.nc"

pixcnc <- nc_open(raw_pixcfile)
pcvnc <- nc_open(raw_pcvfile)

pcvnc_test <- nc_subset(pcvnc, node_index %in% (node_index[1] + 0:1),
                        filename = test_pcvfile)
rangeinds_test <- ncvar_get(pcvnc_test, "range_index")
aziminds_test <- ncvar_get(pcvnc_test, "azimuth_index")
pixcnc_test <- nc_subset(pixcnc,
                         `pixel_cloud/range_index` %in% rangeinds_test,
                         `pixel_cloud/azimuth_index` %in% aziminds_test,
                         `tvp/latitude` %in% `tvp/latitude`[1:100],
                         `noise/num_lines` %in% 1:100,
                         filename = test_pixcfile)
nc_close(pcvnc)
nc_close(pixcnc)
nc_close(pixcnc_test)
nc_close(pcvnc_test)

# rivertile file is small enough not to have to subset. TODO: make this smaller
rtnc <- nc_open(raw_rtfile)
rtnc_test <- nc_subset(rtnc,
                       `reaches/reach_id` %in% (`reaches/reach_id`[1] + 0:1),
                       `nodes/reach_id` %in% (`nodes/reach_id`[1] + 0:1),
                       filename = "tests/testthat/rt.nc")
nc_close(rtnc_test)
nc_close(rtnc)


# Prior reach database ----------------------------------------------------
raw_clfile <- "~/Documents/swot-error/data/priordb-update/Sac_sample_db15.nc"
test_clfile <- "tests/testthat/centerline.nc"
file.copy(raw_clfile, test_clfile, overwrite = TRUE)

# This one is messed up (bad dimension names), so I have to do more work.
clnc <- nc_open(raw_clfile)

reachdimvars <- sapply(clnc$var,
                       function(x) sum(sapply(x$dim,
                          function(y) grepl("^reaches",
                                            y$name))))
badvars <- reachdimvars > 0 &
  sapply(clnc$var, function(x) grepl("^area_fits", x$name))
newdim <- clnc$dim$`reaches/numb_reaches`
newdim$name <- "area_fits/numb_reaches"
gooddim <- clnc$dim$`area_fits/nCoeffs`
newdim$id <- length(clnc$dim) + 1L
newdim$group_id <- gooddim$group_id
newdim$group_index <- gooddim$group_index


clnc$ndims <- newdim$id
clnc$dim[[newdim$name]] <- newdim

for (i in which(badvars)) {
  baddimi <- which(sapply(clnc$var[[i]]$dim,
                          function(x) x$name) == "reaches/numb_reaches")
  if (!length(baddimi)) print("HEY!")
  clnc$var[[i]]$dim[[baddimi]] <- newdim
}

# Again for discharge_models!

badvars <- reachdimvars > 0 &
  sapply(clnc$var, function(x) grepl("^discharge_models", x$name))
newdim <- clnc$dim$`reaches/numb_reaches`
newdim$name <- "discharge_models/numb_reaches"
gooddimind <- clnc$var$`discharge_models/MetroMan_Abar`$group_index
newdim$id <- length(clnc$dim) + 1L
newdim$group_id <- NULL # gooddim$group_id
newdim$group_index <- gooddimind


clnc$ndims <- newdim$id
clnc$dim[[newdim$name]] <- newdim

for (i in which(badvars)) {
  baddimi <- which(sapply(clnc$var[[i]]$dim,
                          function(x) x$name) == "reaches/numb_reaches")
  if (!length(baddimi)) print("HEY!")
  clnc$var[[i]]$dim[[baddimi]] <- newdim
}


clnc_test <- nc_subset(clnc, `nodes/reach_id` == `nodes/reach_id`[1],
                       `reaches/reach_id` == `reaches/reach_id`[1],
                       filename = test_clfile)
nc_close(clnc_test)


