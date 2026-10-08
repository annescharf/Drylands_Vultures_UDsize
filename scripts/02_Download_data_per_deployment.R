# ---
# Title: Download data per deployment
# Author: Anne K Scharf, MPIAB
# Date: February 2025
# Description: this scripts downloads data from Movebank, each individual and
#             deployment are downloaded separately.
#             License terms are accepted automatically. All license terms are
#             contained in the table "license_terms_table.csv" created in this script.
#             This scripts includes solutions if data download is done several
#             times, checking which deployments have not yet been downloaded, and
#             which deployments need to be downloaded again as they contain more
#             data now (CASE_1).
#             It also contains a fix if the download got interrupted because of
#             bad internet connection, to check what has already been downloaded,
#             and what still needs to be done (CASE_2).
#
# ---


library("move2")
library("bit64")
library("units")
library("lubridate")
library("R.utils")


genPath <- "/path_to_folder/"
dir.create(paste0(genPath, "1.vultureIndv_mv2"))
pthDownld <- paste0(genPath, "1.vultureIndv_mv2/")

## specify Movebank account to use
keyring::key_list()
options("move2_movebank_key_name" = "your_movebank_account")

## studies to download
sharedStudies <- readRDS(paste0(genPath, "table_shared_studies_20250221.rds"))
publicStudies <- readRDS(paste0(genPath, "table_public_studies_20250221.rds"))
all <- rbind(sharedStudies, publicStudies) # 105

### ----####
### license terms table ###
### ----####
## save a table with license type and license terms to have them all in one place ##
license_table <- all[, c("id", "name", "taxon_ids", "principal_investigator_name", "contact_person_name", "license_type", "license_terms")]
write.csv(license_table, paste0(genPath, "license_terms_table.csv"))


### ----####
### data overview ###
### ----####
# some summary info on the available data (with public studies)
range(c(
  year(all$timestamp_first_deployed_location), # from 2006 (1994)
  year(all$timestamp_last_deployed_location)
), na.rm = T)
sum(all$number_of_individuals) # 3382 individuals
length(unique(all$id)) # 105 MB studies
length(unique(all$contact_person_name)) # 53 data owners
sort(unique(unlist(strsplit(all$taxon_ids, ",")))) # 17 species
sum(all$number_of_deployed_locations) # 691375874 locations

## check which studies are already downloaded. If last timestamp after last download date, download again, if not, jump. For 1st time downloading directly run from L77

lastDwld <- "2024-12-13" ## date of last download if dataset needs to updated (eg because new data have come in)
thisDwnl <- "2025-02-21" ## if interrupted (date of download that was interrupted, probably today)

##################
## when download is interrupted in the middle (due ie to interrupted internet connection)
tb_pth <- data.frame(pth = list.files(pthDownld, full.names = T), filenm = list.files(pthDownld, full.names = F), mbid = sub("_.*", "", list.files(pthDownld)))

tb_pth$lastSaved <- do.call(c, (lapply(tb_pth$pth, lastModified)))

IDs_doneMdl <- unique(tb_pth$mbid[tb_pth$lastSaved >= as.POSIXct(thisDwnl)])
###################


IDs_done <- unique(sub("_.*", "", list.files(pthDownld)))
compareDF <- all
compareDF$status <- NA
compareDF$status[compareDF$id %in% IDs_done] <- "done"
compareDF$status[compareDF$status == "done" & 
                   compareDF$timestamp_last_deployed_location > as.POSIXct(lastDwld, tz = "UTC")] <- "live"
table(compareDF$status)
table(is.na(compareDF$status))

Ids_toDo <- all$id ## first round
# CASE 1: and # subsequent rounds
Ids_toDo <- compareDF$id[!compareDF$status %in% c("done")]
# CASE 2: and ## interrupted in the middle
Ids_toDo <- Ids_toDo[!Ids_toDo %in% IDs_doneMdl]
# Ids_toDo_later <- as.integer64(c(560810760,232133400,66002314,81273742,149548138,1781483507)) ## these studies gave error and were downloaded separately.
# Ids_toDo <- Ids_toDo[!Ids_toDo%in%Ids_toDo_later]

start_time <- Sys.time()

### uncomment as needed
## CASE 1: subsequent download to update dataset
# lastDownload <- tb_pth$filenm[tb_pth$lastSaved >= as.POSIXct(lastDwld)] ## to download missed last time
## CASE 2: if download got interrupted in the middle
# newerToday <- tb_pth$filenm[tb_pth$lastSaved >= as.POSIXct(thisDwnl)] ## to exclude those downloaded today already

results <- lapply(Ids_toDo, function(studyId) {
  try({
    ## create table with individuals&deployments per study, to be able to download per deployment
    class(studyId) <- "integer64" ## lapply changes the class when looping though it
    print(studyId)

    ## if license terms have to be accepted, this is done here. These license terms are recorded in the license_terms_table.csv of this script
    reftb <- tryCatch(
      {
        movebank_download_deployment(study_id = studyId, omit_derived_data = F)
      },
      error = function(e) {
        movebank_download_deployment(
          study_id = studyId, omit_derived_data = F,
          "license-md5" = sub("...Alternat.*", "", sub(".*se-md5.=.", "", as.character(rlang::catch_cnd(movebank_download_study(studyId)))))
        )
      }
    )
    # vroom::problems(reftb)
    reftb <- reftb[reftb$number_of_events > units::set_units(0, "count"), ]
    indiv <- reftb$individual_local_identifier
    if (any(grepl("/", indiv) == T)) {
      indiv <- gsub("/", "-", indiv)
    } ## names can contain many symbols, "/" is the only one that surly gives errors in R
    reftb$individual_local_identifierNObsl <- indiv
    ## creating name for each downloaded file as the combination of study+indiv+deplym, this makes each file easily identifiable
    reftb$pthName <- paste0(studyId, "__", reftb$individual_local_identifierNObsl, "__", reftb$deployment_id, ".rds")
    doneIndv <- list.files(pthDownld)[grep(studyId, list.files(pthDownld))] ## indiv downloaded
    # allStInd <- reftb$pthName ## individuals in study
    # missInd <- allStInd[!allStInd%in%doneIndv] ## missing indiv
    ## identify indv that are live
    reftb$status <- NA
    reftb$status[reftb$pthName %in% doneIndv] <- "done"
    reftb$status[reftb$status == "done" & reftb$timestamp_end_individual > as.POSIXct(lastDwld, tz = "UTC")] <- "live"
    liveAndmissingIndv <- reftb$pthName[!reftb$status %in% c("done")]

    ## select accordingly
    ## CASE 1: subsequent download to update dataset
    # liveAndmissingIndv <- liveAndmissingIndv[!liveAndmissingIndv%in%lastDownload] # to download missed/not present last download
    ## CASE 2: if download got interrupted in the middle
    # liveAndmissingIndv <- liveAndmissingIndv[!liveAndmissingIndv%in%newerToday] # to exclude those downloaded today already

    print(paste0("done:", length(doneIndv), "-todo:", length(liveAndmissingIndv)))

    todoIndv <- reftb$deployment_id[reftb$pthName %in% liveAndmissingIndv] ## using deployment id to have one file per deployment and reduce mistakes downstream. Also single objects might get smaller and better to handle

    ## download each deployment separatly
    results2 <- lapply(todoIndv, function(dpy) {
      try({
        class(dpy) <- "integer64" ## lapply changes the class when looping though it
        print(paste0(studyId, "_", dpy))
        mv2 <- movebank_download_study(studyId,
          sensor_type_id = c("gps", "argos-doppler-shift", "radio-transmitter"),
          deployment_id = dpy,
          timestamp_end = as.POSIXct(Sys.time(), tz = "UTC")
        ) # to avoid locations in the future

        mt_track_id(mv2) <- "deployment_id"

        if ("individual_local_identifier" %in% names(mt_track_data(mv2))) {
          indiv <- as.character(unique(mt_track_data(mv2)$individual_local_identifier))
          if (any(grepl("/", indiv) == T)) {
            indiv <- gsub("/", "-", indiv)
          }
          print(indiv)
        } else if ("individual_local_identifier" %in% names(mv2)) {
          indiv <- as.character(unique(mv2$individual_local_identifier))
          if (any(grepl("/", indiv) == T)) {
            indiv <- gsub("/", "-", indiv)
          }
          print(indiv)
        } else {
          indiv <- mt_track_data(mv2)$individual_id
          print(indiv)
        }
        saveRDS(mv2, file = paste0(pthDownld, studyId, "__", indiv, "__", dpy, ".rds"))
      })
    })
  })
})
end_time <- Sys.time()
end_time - start_time # ~24h

is.error <- function(x) inherits(x, "try-error")
table(vapply(results, is.error, logical(1)))
names(results) <- seq_along(results)
results[vapply(results, is.error, logical(1))]
# Check studies that returned errors:
giveError <- Ids_toDo[vapply(results, is.error, logical(1))]
tberr <- all[as.numeric(all$id) %in% giveError, ]
