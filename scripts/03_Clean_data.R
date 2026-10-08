# ---
# Title: Basic cleaning of data
# Author: Anne K Scharf, MPIAB
# Date: February 2025
# Description: in this script empty locations and duplicated timestamps are removed
# ---


library("move2")
library("units")
library("dplyr")
library("R.utils")


genPath <- "/path_to_folder/"
pthDownld <- "1.vultureIndv_mv2/"
filePath <- paste0(genPath, pthDownld)
dir.create(paste0(genPath, "2.vultureIndv_mv2_clean_empty_duply"))
pthClean <- "2.vultureIndv_mv2_clean_empty_duply/"
savePath <- paste0(genPath, pthClean)


flsMV <- list.files(filePath, full.names = F)
done <- list.files(savePath, full.names = F)
flsMV <- flsMV[!flsMV %in% done]

###############
## if data are updated, files will be present from the last round but have been newly downloaded. In this case this piece of code checks for this.
lastDwld <- "2024-12-13"
thisDwnl <- "2025-02-21"
tb_pth <- data.frame(pth = list.files(filePath, full.names = T), filenm = list.files(filePath, full.names = F), mbid = sub("_.*", "", list.files(filePath)))
tb_pth$lastSaved <- do.call(c, (lapply(tb_pth$pth, lastModified)))
FLSs_doneMdl <- unique(tb_pth$filenm[tb_pth$lastSaved >= as.POSIXct(thisDwnl)])

thisCalc <- "2025-02-24"
tb_pth_c <- data.frame(pth = list.files(savePath, full.names = T), filenm = list.files(savePath, full.names = F), mbid = sub("_.*", "", list.files(savePath)))
tb_pth_c$lastSaved <- do.call(c, (lapply(tb_pth_c$pth, lastModified)))
FLSs_done <- unique(tb_pth_c$filenm[tb_pth_c$lastSaved >= as.POSIXct(thisCalc)])

flsMV <- FLSs_doneMdl[!FLSs_doneMdl %in% FLSs_done]
#############


start_time <- Sys.time()
lapply(flsMV, function(indPth) {
  print(indPth)
  vultr <- readRDS(paste0(filePath, indPth))

  ## ensuring structure
  if (!mt_is_track_id_cleaved(vultr)) {
    vultr <- vultr |> dplyr::arrange(mt_track_id(vultr))
  } ## order by tracks
  if (!mt_is_time_ordered(vultr)) {
    vultr <- vultr |> dplyr::arrange(mt_track_id(vultr), mt_time(vultr))
  } # order time within tracks

  ## remove empty locs
  if (!mt_has_no_empty_points(vultr)) {
    vultr <- vultr[!sf::st_is_empty(vultr), ]
  }
  ## sometimes only lat or long are NA
  crds <- sf::st_coordinates(vultr)
  rem <- unique(c(which(is.na(crds[, 1])), which(is.na(crds[, 2]))))
  if (length(rem) > 0) {
    vultr <- vultr[-rem, ]
  }

  ## remove 0,0 coordinates
  rem0 <- which(crds[, 1] == 0 & crds[, 2] == 0)
  if (length(rem0) > 0) {
    vultr <- vultr[-rem0, ]
  }

  ## retain the duplicate entry which contains the least number of columns with NA values
  vultr <- vultr %>%
    mutate(n_na = rowSums(is.na(pick(everything())))) %>%
    arrange(n_na) %>%
    mt_filter_unique(criterion = "first") %>% # this always needs to be "first" because the duplicates get ordered according to the number of columns with NA.
    dplyr::arrange(mt_track_id()) %>%
    dplyr::arrange(mt_track_id(), mt_time())

  saveRDS(vultr, file = paste0(savePath, indPth))
  rm(vultr)
  gc()
})
end_time <- Sys.time()
end_time - start_time # 4.5h
