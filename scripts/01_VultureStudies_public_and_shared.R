# ---
# Title: Create table of studies to download from Movebank
# Author: Anne K Scharf, MPIAB
# Date: February 2025
# Description: this scripts gatherers all studies to download, those shared with
#             the specific movebank user and those that are publicly available
# ---

library("move2")
library("units")
library("dplyr")

keyring::key_list()
options("move2_movebank_key_name" = "your_movebank_account")

genPath <- "/path_to_folder/"

### ----####
## Get studies of specified user ##
### ----####
#### downloading studies to which the user "Drylands MPIAB" has been added as collaborator or manager
all <- movebank_download_study_info(study_permission = c("data_manager", "collaborator"))
saveRDS(all, paste0(genPath, "table_shared_studies_20250221.rds"))

### ----####
## Get public studies ##
### ----####
#### searching for studies with vulture name in taxon or study name, which have not been shared with Drylands user and have license terms that make them publicly available
studiesAll <- movebank_download_study_info()
studiesAll <- studiesAll[studiesAll$number_of_deployed_locations > units::set_units(0, "count"), ] # with locations
studiesAll <- studiesAll[!is.na(studiesAll$number_of_deployed_locations), ] # with locations
studiesNoTaxon <- studiesAll[is.na(studiesAll$taxon_ids), ] ## with NO taxon

#### getting vulture studies looking for genus and common name in taxon name, and in study name ####
vult <- c(
  "Gyps", "Necrosyrtes", "Trigonoceps", "Torgos", "Vultur", "Sarcogyps",
  "Neophron", "Gypaetus", "Aegypius", "Gymnogyps", "Coragyps", "Cathartes", "Sarcoramphus", "Gypohierax"
)
vultSmall <- c(
  "gyps", "necrosyrtes", "trigonoceps", "torgos", "vultur", "sarcogyps",
  "neophron", "gypaetus", "aegypius", "gymnogyps", "coragyps", "cathartes", "sarcoramphus", "gypohierax"
)
vultCom <- c("griffon", "vultures", "condor")

studiesVult <- studiesAll[grep(studiesAll$taxon_ids, pattern = paste(vult, collapse = "|")), ] # in taxon column
studiesVultNs <- studiesNoTaxon[grep(studiesNoTaxon$name, pattern = paste(vult, collapse = "|")), ] # in study name when taxon not present
studiesVultSmall <- studiesNoTaxon[grep(studiesNoTaxon$name, pattern = paste(vultSmall, collapse = "|")), ] #  in study name when taxon not present
studiesVultCom <- studiesNoTaxon[grep(studiesNoTaxon$name, pattern = paste(vultCom, collapse = "|")), ] #  in study name when taxon not present

allID <- c(studiesVult$id, studiesVultNs$id, studiesVultSmall$id, studiesVultCom$id)
table(duplicated(allID))

studiesVultures <- rbind(studiesVult, studiesVultNs, studiesVultSmall, studiesVultCom)
studiesVultures <- studiesVultures[!duplicated(studiesVultures$id), ]

#### several checks
table(studiesVultures$i_have_download_access)
table(studiesVultures$license_type)

#### table of all potential vulture studies with download access
alldwn <- studiesVultures[studiesVultures$i_have_download_access == T, ]

#### all given access to
all <- movebank_download_study_info(study_permission = c("data_manager", "collaborator"))
alldwn$id[!alldwn$id %in% all$id]

#### studies that are not shared with project
dwnNotDP <- alldwn[alldwn$id %in% alldwn$id[!alldwn$id %in% all$id], ]

tb <- dwnNotDP[, c("id", "name", "taxon_ids", "license_type", "license_terms")]

remv <- c("Poultry network Thailand 2022", "Thailand ducks 2022", "Movement syndromes across vertebrate taxa (data from Abrahms et al. 2017)", "Peru Los Amigos SigFox ") # studies that do not seem to have useful vulture data

dwnNotDP2 <- dwnNotDP[!dwnNotDP$name %in% remv, ]

#### check license and license terms
tb2 <- dwnNotDP2[, c("id", "name", "license_type", "license_terms")]
print(tb2, n = Inf)
tb2 |>
  select("id", "license_type", "license_terms") |>
  print(n = Inf)

rem1 <- c(893458555, 643943383, 16924201) # license custom with no terms (very few indiv, not worthwhile contacting)
rem2 <- c(326568799, 154820583, 186178781, 186178781, 154820583, 326568799) # release of captured individuals
rem3 <- c(217784323) # full study also available

publicStudies <- dwnNotDP2[!dwnNotDP2$id %in% c(rem1, rem2, rem3), ]
saveRDS(publicStudies, paste0(genPath, "table_public_studies_20250221.rds"))
