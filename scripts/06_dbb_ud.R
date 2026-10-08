# ---
# Title: Calculate monthly dBBMM and UD
# Author: Anne K Scharf, MPIAB
# Date: March 2025
# Description: This script calculates the monthly dBBMMvar; monthly dBBMM; 
#             monthly UD; table with UD extracted values (year-month, locations 
#             per month, UD size Km2, UD centroid, UD weighted mean coordinates); 
#             extracts coordinates of UD with associated dBBMM value for posterior annotation.
# Note: This script requires significant computational resources and may need
#       to be run in parts due to memory constraints. Each section in this script 
#       was run in separate jobs on a HPC (high performance computer).
# ---


library("move") ## for the dbbmm & ud functions
library("move2")
library("lubridate")
library("data.table")
library("sf")

library("doParallel")
library("plyr")
registerDoParallel(4) ## chose number of cores (maximum available-1)


genPath <- "/path_to_folder/"

pathFolder <- paste0(genPath, "5.vultureIndv_mv2_1h_outlspeed_dist/")
dir.create(paste0(genPath, "6.vultureIndv_monthlydBBvar"))
pthmonthlydBBvar <- paste0(genPath, "6.vultureIndv_monthlydBBvar/")
dir.create(paste0(genPath, "7.vultureIndv_monthlydBBL"))
pthmonthlydBBL <- paste0(genPath, "7.vultureIndv_monthlydBBL/")
dir.create(paste0(genPath, "8.vultureIndv_monthlyUDL"))
pthmonthlyUDL <- paste0(genPath, "8.vultureIndv_monthlyUDL/")
dir.create(paste0(genPath, "9.vultureIndv_monthlyUDtable"))
monthlyUDtable <- paste0(genPath, "9.vultureIndv_monthlyUDtable/")
dir.create(paste0(genPath, "10.vultureIndv_monthlyDBBcoordinatesSF"))
pthmonthlyDBBcoordinatesSF <- paste0(genPath, "10.vultureIndv_monthlyDBBcoordinatesSF/")

### setting needed for the calculations below ################
rasterLayer <- 1000 ## resolution in mts of raster
locationError <- 20 ## location error
minLocationsMonth <- 360 # min locations required per month
functionmonthlyMotionVariance <- "mean" # how to sumarize the monthly variance
UDpercentage <- 0.99 ## UD used to calculate the area, centroid, coordinates extracted
##################

##############
## monthly dBBvar
##############
month_dBBvar <- function(indvFileRDS) {
  moveObj <- readRDS(paste0(pathFolder, indvFileRDS))

  if ((nrow(moveObj) / length(unique(mt_track_id(moveObj)))) < 360) {
    stop(paste("individual contains less than 360 locations"))
  }

  yearMonth <- paste0(year(mt_time(moveObj)), "-", sprintf("%02d", month(mt_time(moveObj))))
  locPerMonthDF <- data.frame(table(yearMonth))

  monthToInclude <- as.character(locPerMonthDF$yearMonth[locPerMonthDF$Freq >= minLocationsMonth])

  if (length(monthToInclude) == 0) {
    stop(paste("individual contains 0 month with min nb of locations and no migration"))
  }

  aeqd_crs <- mt_aeqd_crs(moveObj, "center", "m")
  moveObj_c <- sf::st_transform(moveObj, aeqd_crs)
  moveObj_mv1 <- to_move(moveObj_c)

  monthlyBurst <- burst(moveObj_mv1, f = as.character(yearMonth[-length(yearMonth)]))
  monthlydBBvar <- brownian.motion.variance.dyn(monthlyBurst, location.error = locationError, margin = 11, window.size = 31)

  saveRDS(monthlydBBvar, file = paste0(pthmonthlydBBvar, indvFileRDS))
}

fls1h <- list.files(pathFolder, full.names = F)
llply(fls1h, function(x) {
  try(month_dBBvar(indvFileRDS = x))
}, .parallel = T)


#################
## monthly dBBMM ##
#################
month_dBB <- function(indvFileRDS) {
  moveObj_mv2 <- readRDS(paste0(pathFolder, indvFileRDS))
  moveObj_mv1 <- to_move(moveObj_mv2)
  monthlydBBvar <- readRDS(paste0(pthmonthlydBBvar, indvFileRDS))
  monthlydBBvar$yearMonth <- paste0(year(mt_time(moveObj_mv2)), "-", sprintf("%02d", month(mt_time(moveObj_mv2))))
  locPerMonthDF <- data.frame(table(monthlydBBvar$yearMonth))
  monthlydBBvar@interest[timeLag(moveObj_mv1, "hours") > 6] <- FALSE ## excluding segments longer than 5 hours from the dbbmm
  monthToInclude <- as.character(locPerMonthDF$Var1[locPerMonthDF$Freq >= minLocationsMonth])
  monthlydBBvar@interest[!monthlydBBvar$yearMonth %in% monthToInclude] <- FALSE # calc UD only on days with enough points and not migrating
  monthlydBBvarL <- split(monthlydBBvar)
  ## removing days with all @interest==FALSE
  monthlydBBvarL_sel <- lapply(monthlydBBvarL, function(x) {
    if (any(x@interest)) {
      x
    } else {
      x <- NULL
    }
  })
  monthlydBBvarL_sel <- monthlydBBvarL_sel[which(!sapply(monthlydBBvarL_sel, is.null))] # removing "NULL" elements from list
  mbbox <- bbox(monthlydBBvar)
  extExpansionInMts <- round((max(mbbox[, 2] - mbbox[, 1])) / 2)
  monthlydBBL0 <- lapply(monthlydBBvarL_sel, function(dBBvar) {
    try({
      db_r <- raster(ext = extent(dBBvar) + c(-extExpansionInMts, extExpansionInMts, -extExpansionInMts, extExpansionInMts), resolution = rasterLayer, crs = projection(dBBvar), vals = 1) # creating an empty raster, as this is the only option not giving constantly error because of the extent being to small
      brownian.bridge.dyn(dBBvar, raster = db_r, location.error = rep(locationError, length(dBBvar)), margin = 11, window.size = 31, time.step = 45 / 15, verbose = F)
    })
  })
  is.error <- function(x) inherits(x, "try-error")
  monthlydBBL <- monthlydBBL0[!vapply(monthlydBBL0, is.error, logical(1))]

  saveRDS(monthlydBBL, file = paste0(pthmonthlydBBL, indvFileRDS))
}

flsdbbvar <- list.files(pthmonthlydBBvar, full.names = F)
llply(flsdbbvar, function(x) {
  try(month_dBB(indvFileRDS = x))
}, .parallel = T)

################
## monthly UD ##
################
month_UD <- function(indvFileRDS) {
  monthlydBBL <- readRDS(paste0(pthmonthlydBBL, indvFileRDS))

  monthlyUDL <- lapply(monthlydBBL, function(dBB) {
    getVolumeUD(dBB)
  })
  saveRDS(monthlyUDL, file = paste0(pthmonthlyUDL, indvFileRDS))
}

flsdbb <- list.files(pthmonthlydBBL, full.names = F)
llply(flsdbb, function(x) {
  try(month_UD(indvFileRDS = x))
}, .parallel = T)


######################
## monthly UD table ##
######################
month_UDtab <- function(indvFileRDS) {
  monthlyUDL <- readRDS(paste0(pthmonthlyUDL, indvFileRDS))

  moveObj_mv2 <- readRDS(paste0(pathFolder, indvFileRDS))
  moveObj <- to_move(moveObj_mv2)
  yearMonth <- paste0(year(timestamps(moveObj)), "-", sprintf("%02d", month(timestamps(moveObj))))
  locPerMonthDF <- data.frame(table(yearMonth))

  UDsizeL <- lapply(monthlyUDL, function(ud) {
    UDsel <- ud <= UDpercentage
    UDsizem2 <- cellStats(UDsel, "sum") * rasterLayer * rasterLayer
    UDsizeKm2 <- UDsizem2 / 1000000
    return(UDsizeKm2)
  })

  dateschr <- names(monthlyUDL)
  dateschr <- gsub("X", "", dateschr)
  datesDF <- data.frame(yearMonth = gsub("\\.", "-", dateschr))
  includedDatesNbLoc <- merge(datesDF, locPerMonthDF, by = "yearMonth", all.x = T)

  #### monthly UD centroid ##
  centroidsL <- lapply(monthlyUDL, function(ud) {
    UDsel <- ud <= UDpercentage
    UDsel[UDsel == 0] <- NA
    stars_ud <- st_as_stars(UDsel)
    sf_ud <- st_as_sf(stars_ud)
    sf_ud_u <- st_union(sf_ud)
    centr <- st_centroid(sf_ud_u)
    centr_ll <- st_transform(centr, projection(moveObj))
    centr_coor <- st_coordinates(centr_ll)
    return(centr_coor)
  })
  centroids <- do.call("rbind", centroidsL)

  #### monthly dBB weighted mean coordinates ##
  monthlydBBL <- readRDS(paste0(pthmonthlydBBL, indvFileRDS))
  dBBwCoordinatesL <- lapply(1:length(monthlydBBL), function(x) {
    # print(x)
    ud <- monthlyUDL[[x]]
    dbb <- monthlydBBL[[x]]
    dbb[ud > UDpercentage] <- NA
    stars_dbb <- st_as_stars(dbb)
    sf_dbb <- st_as_sf(stars_dbb)
    sf_dbb_c <- st_centroid(sf_dbb)
    dbb_spdf_ll <- st_transform(sf_dbb_c, projection(moveObj))
    wMeanLong <- weighted.mean(st_coordinates(dbb_spdf_ll)[, 1], dbb_spdf_ll$layer)
    wMeanLat <- weighted.mean(st_coordinates(dbb_spdf_ll)[, 2], dbb_spdf_ll$layer)
    wMeanLL <- data.frame(date = names(monthlydBBL[x]), wMeanLong = wMeanLong, wMeanLat = wMeanLat)
    return(wMeanLL)
  })
  dBBwCoordinates <- do.call("rbind", dBBwCoordinatesL)

  #### table of all results of UD calculations ##
  monthlyUDcalc <- data.frame(
    fileName = rep(indvFileRDS, length(includedDatesNbLoc$yearMonth)),
    yearMonth = includedDatesNbLoc$yearMonth,
    locsPerMonth = includedDatesNbLoc$Freq,
    UDsizeKm2 = unlist(UDsizeL),
    UDcentroidsLongitude = centroids[, 1],
    UDcentroidsLatitude = centroids[, 2],
    UDwMeanLongitude = dBBwCoordinates$wMeanLong,
    UDwMeanLatitude = dBBwCoordinates$wMeanLat,
    row.names = NULL
  )

  saveRDS(monthlyUDcalc, file = paste0(monthlyUDtable, indvFileRDS))
}

flsUD <- list.files(pthmonthlyUDL, full.names = F)
llply(flsUD, function(x) {
  try(month_UDtab(indvFileRDS = x))
}, .parallel = T)

###################################
## monthly dBB coordinates as sf ##
###################################
month_DBBcoords <- function(indvFileRDS) {
  monthlyUDL <- readRDS(paste0(pthmonthlyUDL, indvFileRDS))
  monthlydBBL <- readRDS(paste0(pthmonthlydBBL, indvFileRDS))
  moveObj_mv2 <- readRDS(paste0(pathFolder, indvFileRDS))

  dBBcoordinatesL <- lapply(1:length(monthlydBBL), function(x) {
    ud <- monthlyUDL[[x]]
    dbb <- monthlydBBL[[x]]
    dbb[ud > UDpercentage] <- NA
    stars_dbb <- st_as_stars(dbb)
    sf_dbb <- st_as_sf(stars_dbb)
    sf_dbb_c <- st_centroid(sf_dbb)
    dbb_sf_ll <- st_transform(sf_dbb_c, projection(moveObj_mv2))
    dbb_sf_ll$yearMonth <- gsub("\\.", "-", (as.character(gsub("X", "", names(monthlydBBL[x])))))
    dbb_sf_ll$fileName <- indvFileRDS
    return(dbb_sf_ll)
  })

  dBBcoordinates_sf <- dplyr::bind_rows(dBBcoordinatesL)
  names(dBBcoordinates_sf)[1] <- "dBBvalue"

  saveRDS(dBBcoordinates_sf, file = paste0(pthmonthlyDBBcoordinatesSF, indvFileRDS))
}

flsUD <- list.files(pthmonthlyUDL, full.names = F)
llply(flsUD, function(x) {
  try(month_DBBcoords(indvFileRDS = x))
}, .parallel = T)
