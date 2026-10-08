# ---
# Title: Script for downloading and processing NDVI data
# Author: Anne K Scharf, Max Planck Institute for Animal Behavior, DE.
# Date: July 2025
# Description: downloading rasters from NASA earth data in quarters of the globe, 
#             quarters are merged, and processed by the pixel reliability
# ---

library("appeears")
## vignette: https://cran.r-project.org/web/packages/appeears/vignettes/appeears_vignette.html
library("sf")
library("terra")


##################################################
#### downloading rasters from NASA earth data ####
##################################################

## downloading monthly 1km ndvi:
## MOD13A3 v061 - MODIS/Terra Vegetation Indices Monthly L3 Global 1 km SIN Grid

## monthly product is calculated with a weighted temporal average of the 16day product: "In generating this monthly product, the algorithm ingests all the MOD13A2 products that overlap the month and employs a weighted temporal average."
## so layer of 2025-05-01 is the result of all layers collected in May 2025.

## download of global data had to be done in chucks, ie quarter globe tiles and than stitched together (second part of this script)

genPath <- "/path_to_folder/NDVI_monthly/"
dir.create(paste0(genPath, "quarters"))
qpath <- paste0(genPath, "quarters/")
dir.create(paste0(genPath, "global_MOD13A3.061__1_km_monthly_NDVI"))
clnPth <- paste0(genPath, "global_MOD13A3.061__1_km_monthly_NDVI/")

options(keyring_backend = "file") # at least on ubuntu this is needed. You will have to provide a pw for keyring that will be asked for every session you need to access your credentials
## this only needs to be run once, to store the credentials in keyring, after that only the pw of keyring will be asked for
rs_set_key(
  user = "earthdata_username",
  password = "earthdata_pw"
)

## find out with products are available
products <- rs_products()
# Example: MOD13A3.061 is a common global monthly NDVI product

## find out which layers are availabe for a chosen product
ndvi_layers <- rs_layers("MOD13A3.061")
print(ndvi_layers)


### define the polygon of the area of interest. You can also read in a polygon. It must be a sf object (sfc needs to be converted to sf)

## the entire globe cannot be downloaded. Half a globe also does not work, but quarter globe workes nicely
## defining one polygon per quarter globe
roiNE <- st_as_sf(data.frame(
  id = "NE",
  geom = "POLYGON((0 0, 180 0, 180 90, 0 90, 0 0))"
), wkt = "geom", crs = 4326)

roiNW <- st_as_sf(data.frame(
  id = "NW",
  geom = "POLYGON((-180 0, 0 0, 0 90, -180 90, -180 0))"
), wkt = "geom", crs = 4326)

roiSE <- st_as_sf(data.frame(
  id = "SE",
  geom = "POLYGON((0 -90, 180 -90, 180 0, 0 0, 0 -90))"
), wkt = "geom", crs = 4326)

roiSW <- st_as_sf(data.frame(
  id = "SW",
  geom = "POLYGON((-180 -90, 0 -90, 0 0, -180 0, -180 -90))"
), wkt = "geom", crs = 4326)


## defining start and end years and which layers to download
yrs <- 2000:2025
startT <- paste0(yrs, "-01-01")
startT[startT == "2000-01-01"] <- "2000-02-01" ## terra modis ndvi data starts here
endT <- paste0(yrs, "-12-31")
endT[endT == "2025-12-31"] <- "2025-05-31" ## One can also use "today" as the end date, but this way it is clearer till when the data go. Today is 9th June, and the last available layer is from 1st may.
myproduct <- "MOD13A3.061"
mylayers <- c("_1_km_monthly_NDVI", "_1_km_monthly_pixel_reliability")

## creating data.frames with details of query. A folder will be created with the same name as given in the argument "task". Make sure these are unique for each query.

dfNE_list <- list()
for (i in 1:length(yrs)) {
  dfNE_list[[paste0("dfNE_", yrs[i])]] <- data.frame(
    task = paste0("NE_", yrs[i]), ## a folder will be created with this name when the data are downloaded
    subtask = "global_ndvi",
    start = startT[i],
    end = endT[i],
    product = myproduct,
    layer = mylayers
  )
}

dfNW_list <- list()
for (i in 1:length(yrs)) {
  dfNW_list[[paste0("dfNW_", yrs[i])]] <- data.frame(
    task = paste0("NW_", yrs[i]),
    subtask = "global_ndvi",
    start = startT[i],
    end = endT[i],
    product = myproduct,
    layer = mylayers
  )
}


dfSE_list <- list()
for (i in 1:length(yrs)) {
  dfSE_list[[paste0("dfSE_", yrs[i])]] <- data.frame(
    task = paste0("SE_", yrs[i]),
    subtask = "global_ndvi",
    start = startT[i],
    end = endT[i],
    product = myproduct,
    layer = mylayers
  )
}

dfSW_list <- list()
for (i in 1:length(yrs)) {
  dfSW_list[[paste0("dfSW_", yrs[i])]] <- data.frame(
    task = paste0("SW_", yrs[i]),
    subtask = "global_ndvi",
    start = startT[i],
    end = endT[i],
    product = myproduct,
    layer = mylayers
  )
}

####################################
#### Build the task with the ROI ###
####################################
seq5yrs <- seq(1, length(yrs), by = 4) ### dividing years into chuncks of 3, -- total time 48h, for 12 workers, each hast 4h... hopefully enough

# lapply(seq5yrs, function(y){ ## keeps giving error, introducing manually "y". for some reason the looping does not work as hopped
y <- 21
task_list <- lapply(y:(y + 3), function(i) {
  list(
    rs_build_task(df = dfNE_list[[i]], roi = roiNE), # , format="netcdf4") ## netcdf4 files do not contain the date in the name, but has to be looked up in a separate table (in a not straight forward way). It is a bit lighter, but the probability of messing up dates was to high for me
    rs_build_task(df = dfNW_list[[i]], roi = roiNW), # , format="netcdf4")
    rs_build_task(df = dfSE_list[[i]], roi = roiSE), # , format="netcdf4")
    rs_build_task(df = dfSW_list[[i]], roi = roiSW) # , format="netcdf4")
  )
})
flat_task_list <- unlist(task_list, recursive = FALSE)

# t0 <- Sys.time()
rs_request_batch(
  request_list = flat_task_list, # task_l
  user = "earthdata_username",
  path = qpath,
  time_out = 14400, # 4h
  workers = 12, # max 20
  verbose = TRUE,
  total_timeout = 172800 # max 48h
)
# Sys.time()-t0

# })

## ERRORS:
# Error in if (private$status != "done") { : argument is of length zero
## this error mostly (but not always) comes if one deletes a download in the appeears webpage. Than the request will be processed, but never downloaded to your pc

## check status of request: https://appeears.earthdatacloud.nasa.gov/explore
## do not delete anything while the R script is still running!


#####################################
#### merging the global quarters #####
#####################################

### NDVI -> FLT4S (datatype)
# # 0.3-1 -- vegetation
# # 0-0.3 -- bare soil
# # neg -- no data (??)

### pixel_reliability -> INT1U (datatype)
# # 0	Good data, use with confidence
# # 1	Marginal data, Useful, but look at other QA information
# # 2	Snow/Ice Target covered with snow/ice
# # 3	Cloudy data, Target not visible, covered with cloud


yrs <- 2000:2025

# yr <- 2024
lapply(yrs, function(yr) {
  print(yr)
  dir.create(paste0(genPath, yr))

  ## ndvi
  flsNs <- list.files(paste0(qpath, "NE_", yr, "/"), pattern = "_NDVI_", full.names = F)
  fldrsNs <- paste0(c("NE_", "NW_", "SE_", "SW_"), yr)
  r_pth <- paste0(qpath, fldrsNs)

  # mth <- flsNs[1]
  lapply(flsNs, function(mth) {
    pthL <- paste0(r_pth, "/", mth)
    rstL <- lapply(pthL, rast)
    mm <- mosaic(rstL[[1]], rstL[[2]], rstL[[3]], rstL[[4]],
      fun = "first", filename = paste0(genPath, yr, "/", mth), overwrite = T,
      wopt = list(datatype = "FLT4S")
    )
    # plot(mm)
  })

  ## pixel_reliability
  R_flsNs <- list.files(paste0(qpath, "NE_", yr, "/"), pattern = "_pixel_reliability", full.names = F)
  R_fldrsNs <- paste0(c("NE_", "NW_", "SE_", "SW_"), yr)
  R_r_pth <- paste0(qpath, R_fldrsNs)

  # mth <- R_flsNs[1]
  lapply(R_flsNs, function(mth) {
    R_pthL <- paste0(R_r_pth, "/", mth)
    R_rstL <- lapply(R_pthL, rast)
    R_mm <- mosaic(R_rstL[[1]], R_rstL[[2]], R_rstL[[3]], R_rstL[[4]],
      fun = "first", filename = paste0(genPath, yr, "/", mth), overwrite = T,
      wopt = list(datatype = "INT1U")
    )
    # plot(R_mm)
  })
})


#######################################
### Script for processing NDVI data ###
#######################################

# MODIS/Terra Vegetation Indices monthly Global 1km V061 - NDVI and pixel reliability
# Downloaded from earthdata.nasa.gov for the years 2000-2025 (see script "A.earthdata_NDVI_global_download_and_stitching.R")
# Download was done in quarters, stitched to one raster for globe in yearly subfolders
# results from this script are all rasters in one folder named by date
## sample file names
## "MOD13A3.061__1_km_monthly_NDVI_doy2000061000000_aid0001.tif"
## "MOD13A3.061__1_km_monthly_pixel_reliability_doy2000061000000_aid0001.tif"
## filtered raster will be saved as:
## "20000301.tif" # date: "2000-03-01"


### pixel_reliability
# # 0	Good data, use with confidence
# # 1	Marginal data, Useful, but look at other QA information
# # 2	Snow/Ice Target covered with snow/ice
# # 3	Cloudy data, Target not visible, covered with cloud


yrs <- 2014:2025

# t0 <- Sys.time()
# yr <- 2024
lapply(yrs, function(yr) {
  print(yr)

  flsNs <- list.files(paste0(genPath, yr, "/"), pattern = "_NDVI_", full.names = F)
  ydoyVctr <- gsub(".*doy([0-9]+)_.*", "\\1", flsNs) # extracting year doy date
  # ydoy <- ydoyVctr[1]
  lapply(ydoyVctr, function(ydoy) {
    ndvi <- rast(list.files(paste0(genPath, yr, "/"), pattern = paste0("NDVI_doy", ydoy), full.names = T))
    pxreliab <- rast(list.files(paste0(genPath, yr, "/"), pattern = paste0("pixel_reliability_doy", ydoy), full.names = T))

    # mask only clouds, snow/ice kept to its original values for now
    ndvi[pxreliab == 3] <- NA # cloud

    yr_doy <- sub("^([0-9]{4})(.*)", "\\1-\\2", ydoy) ## separating year and doy
    dte <- as.Date(yr_doy, format = "%Y-%j") ## converting to date
    flname <- gsub("-", "", dte) ## converting date to filename format
    writeRaster(ndvi, paste0(clnPth, flname, ".tif"), overwrite = T, datatype = "FLT4S")
  })
})

# Sys.time()-t0 ## ca 30mins per year
