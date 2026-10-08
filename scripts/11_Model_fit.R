# ---
# Title: GAM fitting
# Author: Dongmin (Dennis) Kim, UMN; Harvard
# Date: April 2026
# Description: Fit the final models
# ---

# call libraries 
library(here)
library(MASS)
library(nlme)
library(mgcv)
library(dplyr)
library(tidyr)
library(ggplot2)

genPath <- "/path_to_folder/"

final_df <- readRDS(here::here(genPath, "Final_dryland_data_may2026.rds"))

# 1. model: no interaction terms
model1 <- bam(
  dsl_UD ~
    # Fixed effects:
    # Species-specific intercepts (baseline UD differences)
    # and REALM-level differences shared across species
    # and dryland-level differences shared across species
    species + REALM+ wmodeDrylands+
    
    # Environmental predictors:
    # Species-specific smooths allow each species
    # to respond differently (nonlinearly) to NDVI
    s(dsl_wmNDVI, by = species, k = 5) +
    
    # Species-specific nonlinear response to livestock density
    s(wmLifestock, by = species, k = 5) +
    
    # Sampling effort correction:
    # Nonlinear effect of number of locations per month,
    # allowed to differ by species (NOT a random effect)
    # This controls for effort-induced bias in UD estimation
    s(totNbLocs, by = species, k =5)+
    
    # Spatial structure:
    # 2D smooth over centroid longitude/latitude,
    # with separate spatial surfaces for each species
    # (accounts for broad-scale spatial autocorrelation)
    s(UDwMeanLongitude, UDwMeanLatitude, by = species, k = 100) +
    
    # Individual-level random intercept:
    # Accounts for repeated measures within animals
    # (IDs can repeat across species)
    s(animal_id, bs = "re"),
  
  data = final_df,
  family = gaussian(),
  method = "fREML",   # stable & efficient
  discrete = TRUE     # good for large datasets
)

# 2. Finalized Model format: dryland and REALM
model2 <- bam(
  dsl_UD ~
    # Fixed effects:
    # Species-specific intercepts (baseline UD differences)
    # and REALM-level differences shared across species
    # and dryland-level differences shared across species
    species + REALM+ wmodeDrylands+
    
    # Environmental predictors:
    # Species-specific smooths allow each species
    # to respond differently (nonlinearly) to NDVI
    s(dsl_wmNDVI, by = interaction(species, wmodeDrylands),  k = 5) +
    
    # Species-specific nonlinear response to livestock density
    s(wmLifestock, by = interaction(species, wmodeDrylands),  k = 5) +
    
    # Sampling effort correction:
    # Nonlinear effect of number of locations per month,
    # allowed to differ by species (NOT a random effect)
    # This controls for effort-induced bias in UD estimation
    s(totNbLocs, by = species, k =5)+
    
    # Spatial structure:
    # 2D smooth over centroid longitude/latitude,
    # with separate spatial surfaces for each species
    # (accounts for broad-scale spatial autocorrelation)
    s(UDwMeanLongitude, UDwMeanLatitude, by = species, k = 100) +
    
    # Individual-level random intercept:
    # Accounts for repeated measures within animals
    # (IDs can repeat across species)
    s(animal_id, bs = "re"),
  
  data = final_df,
  family = gaussian(),
  method = "fREML",   # stable & efficient
  discrete = TRUE     # good for large datasets
)

# 2. Finalized Model format: dryland only
model3 <- bam(
  dsl_UD ~
    # Fixed effects:
    # Species-specific intercepts (baseline UD differences)
    # and REALM-level differences shared across species
    # and dryland-level differences shared across species
    species + wmodeDrylands+
    
    # Environmental predictors:
    # Species-specific smooths allow each species
    # to respond differently (nonlinearly) to NDVI
    s(dsl_wmNDVI, by = interaction(species, wmodeDrylands),  k = 5) +
    
    # Species-specific nonlinear response to livestock density
    s(wmLifestock, by = interaction(species, wmodeDrylands),  k = 5) +
    
    # Sampling effort correction:
    # Nonlinear effect of number of locations per month,
    # allowed to differ by species (NOT a random effect)
    # This controls for effort-induced bias in UD estimation
    s(totNbLocs, by = species, k =5)+
    
    # Spatial structure:
    # 2D smooth over centroid longitude/latitude,
    # with separate spatial surfaces for each species
    # (accounts for broad-scale spatial autocorrelation)
    s(UDwMeanLongitude, UDwMeanLatitude, by = species, k = 100) +
    
    # Individual-level random intercept:
    # Accounts for repeated measures within animals
    # (IDs can repeat across species)
    s(animal_id, bs = "re"),
  
  data = final_df,
  family = gaussian(),
  method = "fREML",   # stable & efficient
  discrete = TRUE     # good for large datasets
)

# Finalized Model format: dryland only (grouping by dryland vs non-dryland)
model4 <- bam(
  dsl_UD ~
    # Fixed effects:
    # Species-specific intercepts (baseline UD differences)
    # and REALM-level differences shared across species
    # and dryland-level differences shared across species
    species + wmodeDrylands1+
    
    # Environmental predictors:
    # Species-specific smooths allow each species
    # to respond differently (nonlinearly) to NDVI
    s(dsl_wmNDVI, by = interaction(species, wmodeDrylands1),  k = 5) +
    
    # Species-specific nonlinear response to livestock density
    s(wmLifestock, by = interaction(species, wmodeDrylands1),  k = 5) +
    
    # Sampling effort correction:
    # Nonlinear effect of number of locations per month,
    # allowed to differ by species (NOT a random effect)
    # This controls for effort-induced bias in UD estimation
    s(totNbLocs, by = species, k =5)+
    
    # Spatial structure:
    # 2D smooth over centroid longitude/latitude,
    # with separate spatial surfaces for each species
    # (accounts for broad-scale spatial autocorrelation)
    s(UDwMeanLongitude, UDwMeanLatitude, by = species, k = 100) +
    
    # Individual-level random intercept:
    # Accounts for repeated measures within animals
    # (IDs can repeat across species)
    s(animal_id, bs = "re"),
  
  data = final_df1,
  family = gaussian(),
  method = "fREML",   # stable & efficient
  discrete = TRUE     # good for large datasets
)

# save the final model and modified data frame 
saveRDS(model1, file="model1_may2026.rds")
saveRDS(model2, file = "model2_may20262026.rds")
saveRDS(model3, file = "model3_may2026.rds")
saveRDS(model4, file = "model4_may2026.rds")
