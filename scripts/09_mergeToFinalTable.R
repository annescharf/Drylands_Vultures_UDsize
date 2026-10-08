# ---
# Title: Creating table for models
# Author: Anne K Scharf, MPIAB
# Date: July 2025
# Description: Merges all tables, reference data tbele, tables from 
#             "9.vultureIndv_monthlyUDtable", annotates ecoregions, table with 
#             annotated NDVI and livestock values
# ---


library(sf)

genPath <- "/path_to_folder/"
annotcorrdspth <- paste0(genPath, "/11.vultureIndv_UDannotation/")
monthlyUDtablept <- paste0(genPath, "/9.vultureIndv_monthlyUDtable/")
refTb <- readRDS(paste0(genPath, "/referenceTableStudies_ALL_excludedColumn.rds"))
refTb$animal_id <- paste0(refTb$MBid, "_", refTb$individual_local_identifier)

ecoregions <- st_read(paste0(genPath, "Ecoregions2017/Ecoregions2017.shp")) # BIOME_NAME and REALM - "A biogeographic realm is the broadest biogeographic division of Earth's land surface, based on distributional patterns of terrestrial organisms. They are subdivided into bioregions, which are further subdivided into ecoregions. " (Wikipedia)
# SOURCE: Olson et al. Terrestrial Ecoregions of the World: A New Map of Life on Earth: A new global map of terrestrial ecoregions provides an innovative tool for conserving biodiversity, BioScience, Volume 51, Issue 11, November 2001, Pages 933–938, https://doi.org/10.1641/0006-3568(2001)051[0933:TEOTWA]2.0.CO;2
# MAP: https://ecoregions.appspot.com/

# fileName, species, animal_id

udtb_l <- lapply(list.files(monthlyUDtablept, full.names = T), function(x) {
  readRDS(x)
})
udtb_df <- do.call("rbind", udtb_l)
head(udtb_df)

udtb_df_sf <- st_as_sf(udtb_df[, c("fileName", "yearMonth")], coords = c("UDwMeanLongitude", "UDwMeanLatitude"), crs = "EPSG:4326")
udtb_df_sf_pj <- st_transform(udtb_df_sf, st_crs(ecoregions))

sf::sf_use_s2(FALSE)
udtb_df_sf_pj_eco <- st_join(udtb_df_sf_pj, ecoregions)

udtb_df_eco <- merge(udtb_df, udtb_df_sf_pj_eco[, c("fileName", "yearMonth", "BIOME_NAME", "REALM")])
head(udtb_df_eco)
udtb_df_eco$geometry <- NULL

annot_l <- lapply(list.files(annotcorrdspth, full.names = T), function(x) {
  readRDS(x)
})
annot_df <- do.call("rbind", annot_l)
head(annot_df)
tbudann <- merge(udtb_df_eco, annot_df, by = c("fileName", "yearMonth"))

summary(tbudann)

tb_mdls <- merge(tbudann, refTb[, c("fileName", "species", "animal_id")], by = "fileName", all.x = T, all.y = F)
head(tb_mdls)

tb_mdls$UDsizeKm2 <- NULL
tb_mdls$locsPerMonth <- NULL
tb_mdls$UDcentroidsLongitude <- NULL
tb_mdls$UDcentroidsLatitude <- NULL
head(tb_mdls)
saveRDS(tb_mdls, paste0(genPath, "/table_for_models_july2025.rds"))
