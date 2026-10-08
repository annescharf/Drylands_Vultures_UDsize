# ---
# Title: Generate the final data
# Author: Dongmin (Dennis) Kim, UMN; Harvard
# Date: May 2026
# Description: Finalize the data for the model fits.
# ---

genPath <- "/path_to_folder/"

# call the table generated from the "09_mergeToFinalTable.R"
final_df <- readRDS(paste0(genPath, "/table_for_models_july2025.rds"))

# factorize the cateogrical columns (animal_id, REALM, BIOME_NAME) 
final_df$BIOME_NAME <- as.factor(final_df$BIOME_NAME)
final_df$REALM <- as.factor(final_df$REALM)
final_df$animal_id <- as.factor(final_df$animal_id)
final_df$wmodeDrylands <- factor(final_df$wmodeDrylands, level = c("Non-dryland", "Dry subhumid", "Semiarid", "Arid", "Hyperarid"))

# remove species that has small samples (4 species removed - now 13 species)
final_df1 <- final_df %>% dplyr::filter(! species %in% c("Cathartes burrovianus", "Gyps bengalensis", "Gyps rueppellii", "Sarcoramphus papa"))

# scaling NDVI and lifestcok
final_df1$wmLifestock_scale <- scale(final_df1$wmLifestock)
final_df2 <- final_df1 %>% filter(!is.na(REALM) & REALM != "N/A") %>% droplevels()

# final data
final_df3 <- na.omit(final_df2)

# separate the 'yearMonth' column to 'year' and 'month' columns 
final_df3 <- final_df3 %>% separate(yearMonth, into = c("year", "month"), sep = "-", convert = TRUE)

# futher filterization
final_df4 <- final_df3 %>% 
  # select only valid columns for the dataframe
  dplyr::select(species, animal_id, year, month, UDsizeKm2_99, wmNDVI, wmLifestock, wmLifestock_scale, REALM, UDwMeanLongitude, UDwMeanLatitude, REALM, wmodeDrylands) %>% 
  # mutate data column for rearranging the order of year and month for further filtering
  mutate(date = as.Date(paste(year, month, "01", sep = "-"))) %>%
  arrange(animal_id, date) %>% 
  group_by(animal_id) %>% 
  mutate(dsl_UD = UDsizeKm2_99 - lag(UDsizeKm2_99),
         dsl_wmNDVI = wmNDVI - lag(wmNDVI)) %>% 
  ungroup()

# reorder the columns 
final_df5 <- final_df4 %>% dplyr::select(species, animal_id, REALM, wmodeDrylands, UDwMeanLongitude, UDwMeanLatitude, year, month, UDsizeKm2_99, dsl_UD, wmNDVI, dsl_wmNDVI, wmLifestock, wmLifestock_scale)

# load tot_number of blocks data 
tot_df <- readRDS(here::here("data/0.final_data", "Final_dryland_data_mar2026.rds"))
tot_df <- tot_df %>% dplyr::select(species, animal_id, REALM, year, month, UDsizeKm2_99, totNbLocs)

# left join
final_df6 <- final_df5 %>% left_join(., tot_df)

# create a 2-category system of Dryland (Dryland vs Non-Dryland)
final_df7 <- final_df6 %>% mutate(wmodeDrylands1 = ifelse(wmodeDrylands == "Non-dryland", "Non-dryland", "Drylands"))
final_df7$wmodeDrylands1 <-factor(final_df7$wmodeDrylands1, level = c("Non-dryland", "Drylands"))

# save the data
saveRDS(final_df7, paste0(genPath, "/Final_dryland_data_may2026.rds"))
