
#______________________________________________
# Latest analysis, July 2026, Martina & Claude
#______________________________________________

# Two sections:
# S1 is about state: given an individual's chronic aridity, position in the 20-year trend, seasonal position, and anomaly, what's their standing UD size?
# S2 is about reactivity: does a change in NDVI provoke a change in UD, and on what timescale?

library(MASS)
library(nlme)
library(mgcv)        # bam() for large GAMMs
library(gratia)      # smooth_estimates(), derivatives() for GAM plotting
library(dplyr)       # loading it after loading gratia ensures the correct package version is loaded
library(slider)      # rolling means
library(tidyr)
library(ggplot2)
library(ggeffects)
library(emmeans)
library(patchwork)
library(trend)       # Mann-Kendall test + Sen's slope
library(lme4)        # linear mixed models (seasonal + long-term LMMs)
library(lmerTest)    # p-values via Satterthwaite df
library(MuMIn)       # r.squaredGLMM()
library(broom.mixed) # tidy model output


#_____
# Input/Output folders
genPath <- "/path_to_folder/"
outputPath <- "/path_to_folder/"
#genPath <- getwd()
#outputPath <- "/home/martina/ownCloud/Martina/ProgettiVari/Drylands/DryLands/NewAnalysis_June2026/modelResults"

final_df <- readRDS(paste0(genPath, "/data_for_models/table_for_models_july2026_w_filter_columns.rds"))

#____________
# Data checks

final_df %>%
  group_by(wmodeDrylands) %>%
  summarise(n_species = n_distinct(species),
            n_ind = n_distinct(animal_id),
            n_months = n_distinct(month),
            n_years = n_distinct(year),
            n_UDs = n()) %>%
  arrange(desc(n_ind)) %>%
  print(n = Inf)
# hyperarid is the far least represented month, we combine hyperarid and arid into one category (below)


#_______________________________________
# Grouping variables and transformations

final_df <- final_df %>% 
  tidyr::separate(yearMonth, into = c("year", "month"), sep = "-", convert = TRUE) %>%
  mutate(
    # Sampling effort: log transform then scale
    locsPerMonth_logz = as.numeric(scale(log(locsPerMonth))),
    # Center year
    year_c = year - mean(year), # year centered
    # Ensure correct order of months (ordered factor) for seasonal models
    month_f = factor(month.name[month], levels=month.name),
    # Ensure grouping variables are factors
    species = factor(species),
    animal_id = factor(animal_id),
    # Group drylands categories four classes (hyperarid and arid into one)
    wmodeDrylands_4cat = recode(wmodeDrylands, "Hyperarid" = "Arid", .default = wmodeDrylands),
    wmodeDrylands_4cat = factor(wmodeDrylands_4cat, levels = c("Non-dryland", "Dry subhumid", "Semiarid", "Arid")),
    # Group drylands categories into binary
    wmodeDrylands = factor(wmodeDrylands, levels = c("Non-dryland", "Dry subhumid", "Semiarid", "Arid", "Hyperarid")),
    wmodeDrylands_bin = ifelse(wmodeDrylands == "Non-dryland", "Non-dryland", "Drylands"),
    wmodeDrylands_bin = factor(wmodeDrylands_bin, levels=c('Non-dryland','Drylands')),
    # Create a new grouping variables species.dryland
    ssp_dry = interaction(species, wmodeDrylands_bin)
  )

#________________________
# Add a resident category column (see script 10 for details)
final_df$resident_100km <- final_df$dist_to_1st_cntr_km <= 100


#____________________________________________________
### Model long-term trend in proportion of residency
#____________________________________________________

# We model the expected proportion of months in a year that each individual remains resident, and how this proportion shifts over the 20-year period (year_c) and by dryland category. 
# A negative year_c coefficient means the probability of being resident in any given month is declining over the years.
# final_df has one row per animal_id x year x month with the residence_100km column indicating "residency" in a given month

residency_by_year <- final_df %>%
  distinct(animal_id, species, wmodeDrylands_4cat, year, month, resident_100km) %>%
  group_by(animal_id, species, wmodeDrylands_4cat, year) %>%
  summarise(
    n_resident    = sum(resident_100km, na.rm = TRUE),
    n_months      = n_distinct(month),
    prop_resident = n_resident / n_months,
    .groups = "drop"
  ) %>%
  mutate(year_c = year - mean(year))

# The response variable is, for each animal_id × year row, "out of the months we could have observed this individual that year, how many did it spend classified as resident (in that dryland/location) vs. not (migratory/away)
# The random effects (animal_id, species) account for the fact that some individuals/species are just chronically more resident than others, independent of any trend.

range(residency_by_year$n_months)


mod_residency_trend <- glmer(
  cbind(n_resident, n_months - n_resident) ~  # binary 0-1 resident or not in that month
    year_c * wmodeDrylands_4cat +
    n_months +
    (1 | animal_id) + (1 | species),
  data    = residency_by_year,
  family  = binomial,
  control = glmerControl(optimizer = "bobyqa")
)

AIC(mod_residency_trend)
summary(mod_residency_trend)
MuMIn::r.squaredGLMM(mod_residency_trend)

# Important findinf affecting the next models on long-term trends:
# residency itself is most unstable where the arid/semiarid categories are — individuals are increasingly likely to leave/become non-resident there over the same 20 years. 
# This selection-into-residency is a confound for the next models: if birds are dropping out of the resident sample disproportionately from degrading arid/semiarid sites as conditions worsen, 
# the long-term NDVI trend in those dryland types could partly reflect survivorship into residency — the resident sample in arid/semiarid categories in later years may be a shrinking, potentially unrepresentative subset of who's still willing/able to stay. 

#__________
# Plots

# predicted probability of residency by year, faceted by dryland category
pred_residency <- ggpredict(mod_residency_trend, terms = c("year_c [all]", "wmodeDrylands_4cat"))

pan_residency <- ggplot(pred_residency, aes(x = x, y = predicted, colour = group, fill = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 1.2) +
  scale_colour_manual(values = c("Non-dryland"="#008000","Dry subhumid"="#7E8C54",
                                 "Semiarid"="#795C34","Arid"="#DAA520")) +
  scale_fill_manual(values = c("Non-dryland"="#008000","Dry subhumid"="#7E8C54",
                               "Semiarid"="#795C34","Arid"="#DAA520")) +
  labs(title = "Predicted probability of monthly residency over time",
       subtitle = "By dryland category — declining residency signals increasing migratory behaviour",
       x = "Year (centred)", y = "P(resident month)", colour = NULL, fill = NULL) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom")

pan_residency
ggsave(paste0(outputPath, "/residency_trend_by_dryland.png"), pan_residency, width = 8, height = 6, dpi = 300)


#_________________________________
# Select ONLY RESIDENT individuals (displacement of UD centroid within 100 km)
#_________________________________

resident100km <- final_df[final_df$resident_100km == TRUE,]

#________________________________
# Plot residents per dryland type

world <- map_data("world")

ggplot() +
  geom_polygon(data = world, aes(long, lat, group = group),
               fill = "grey95", colour = "grey80", linewidth = 0.2) +
  geom_point(data = resident100km, aes(x = UDwMeanLongitude, y = UDwMeanLatitude, colour = wmodeDrylands),
             alpha = 0.6, size = 1.2) +
  coord_quickmap() +
  labs(x = "Longitude", y = "Latitude", colour = "Dryland type") +
  theme_bw()
world <- map_data("world")

#______________________
# Divide per hemisphere

resident100km <- resident100km %>%
  group_by(animal_id) %>%
  mutate(
    hemisphere = ifelse(mean(UDwMeanLatitude, na.rm = TRUE) >= 0, "Northern", "Southern")) %>%
  ungroup() %>%
  mutate(hemisphere = factor(hemisphere, levels = c("Northern", "Southern")))

#_____________
# Data check

# Check how many individuals and ind/species fall in each hemisphere
resident100km %>%
  group_by(species, wmodeDrylands_bin, hemisphere) %>%
  summarise(n_ind = n_distinct(animal_id),
            n_months = n_distinct(month),
            n_years = n_distinct(year)) %>%
  arrange(desc(n_ind)) %>%
  print(n = Inf)

resident100km %>%
  group_by(wmodeDrylands, hemisphere) %>%
  summarise(n_species = n_distinct(species),
            n_ind = n_distinct(animal_id),
            n_months = n_distinct(month),
            n_years = n_distinct(year),
            n_UDs = n()) %>%
  arrange(desc(n_ind)) %>%
  print(n = Inf)
# hyperarid is the far least represented month

# number of years observations per dryland
resident100km %>%
  group_by(wmodeDrylands_4cat, year) %>%
  summarise(n_species = n_distinct(species),
            n_ind = n_distinct(animal_id),
            n_monthlyUDs = n()) %>%
  arrange(year) %>%
  print(n = Inf)

# sample size per dryland/hemisphere
resident100km %>%
  group_by(wmodeDrylands_4cat, hemisphere) %>%
  summarise(n_species = n_distinct(species),
            n_ind = n_distinct(animal_id),
            n_months = n_distinct(month),
            n_years = n_distinct(year),
            n_UDs = n()) %>%
  arrange(desc(n_ind)) %>%
  print(n = Inf)

# Important check for interpretation of final results:
table(resident100km$species, resident100km$wmodeDrylands_4cat)
# any prediction for a species in a category it was never observed in (e.g. Gyps himalayensis in Semiarid) is a pure extrapolation resting entirely on the assumption that the dryland effect is the same for every species — additivity, not evidence. That's a real limitation to name explicitly (e.g., in a supplementary table listing which species×category cells have zero data, with a caveat that predictions there are model-based extrapolations, not empirically supported


#__________________________________
# Model NDVI shared expected baseline
#___________________________________

# Use a GAM to define the expected NDVI at a given location and time of year, independent of any particular behaviouralor environmental anomaly.
# Here species with few observations won't need special handling with this formulation because the seasonal expectation comes from the continuous spatial surface, so those individuals just borrow the local spatial smooth's estimate at their coordinates rather than needing enough of their own data to support an independent curve.

m_ndvi_trend <- bam(
  wmNDVI ~
    wmodeDrylands_4cat +
    s(year_c, by = wmodeDrylands_4cat, k = 10) +     # use 4-level dryland categories (arid and hyperarid grouped)
    s(year_c, species, bs = "fs", k = 6) +           # long-term trend still allowed to vary by species
    te(UDwMeanLongitude, UDwMeanLatitude, month,     # spatial-seasonal baseline, continuous
       bs = c("tp", "cc"), k = c(45, 8),             # spatial-by-seasonal tensor (check if 30 k for lat-long is enough by looking at k-check)
       d = c(2, 1)) +                                # d = c(2,1) tells mgcv to treat lon/lat as one 2D isotropic marginal and month as a separate 1D cyclic marginal, rather than a fully joint 3D penalty
    s(animal_id, bs = "re"),
  data     = resident100km,
  method   = "fREML",
  family = gaussian(),
  #family = betar(link = "logit"),
  discrete = TRUE,
  knots    = list(month = c(0.5, 12.5)))   # For a cyclic smooth (bs = "cc"), mgcv needs to know where the cycle wraps around (min(month) and max(month) by default, which for month coded 1–12 means it forces month 12 and month 1 to be the same month, with zero gap in between).
# knots month = c(0.5, 12.5) tells the cc term that the wraparound boundary sits at 0.5 and 12.5 ('December (≈12) and January (≈1) are one month apart, connected smoothly through a wrap point that sits just outside the observed range).

#AIC(m_ndvi_trend_bin) - AIC(m_ndvi_trend) # with 4 categories deltaAIC 304 lower than binary dryland definition
summary(m_ndvi_trend)
plot(m_ndvi_trend, page=1)
k.check(m_ndvi_trend)
# In kcheck a low k index and a significant pvalue (and edf/k' value ratio close to 1) indicate that the smooth is using almost every degree of freedom, which is the signature of a ceiling, not a well-supported fit. So this would mean a higher number of k is required.
# ideally we would want a k index close to 1, non significant, and a low edf/k' ratio, but probably there is too much spatial heterogeneity not captured by this data
# NA for random effect type smooths is normal

gam.check(m_ndvi_trend)

# # Run the concurvity check between species and location on a simplified model with fewer levels (remove animal_id and drylands)
# m_concurv_check <- bam(
#   wmNDVI ~
#     te(UDwMeanLongitude, UDwMeanLatitude, month_f,
#        bs = c("tp", "cc"), k = c(30, 8), d = c(2, 1)) +
#     s(species, bs = "re"),
#   data     = resident100km,
#   method   = "fREML",
#   discrete = TRUE,
#   knots    = list(month_f = c(0.5, 12.5)))
# concurvity(m_concurv_check, full = FALSE)
# # look at the pairwise concurvity between the te(lon,lat,month) term and s(species, bs="re"). 
# # the difference between observed and estimate are the important ones (para is the intercept and the others are the covariates effects)
# # te surface (row) relative to species is 0.99 in observed and 0.58 in estimate, so a large part of location variation could be reconstructed by species
# # the species row relative to the te is 0.92 in observed and 0.1 in estimate
# # so species is more redundant given location than location is redundant given species
# # plus species is also included in the year_c smooth so no need to add it as random effect too


# Results summary:
# Time, location and dryland type explain >94% variance in NDVI
# The parametric estimates show that NDVI in Arid < Semiarid < Dry Subhumid < non-dryland, all significant
# Long-term NDVI trend varies across dryland types with a significant long-term trend:
# Non-dryland show a more or less stable trend
# dry subhumid show increasing NDVI over years
# While semiarid and arid show a slight decrease 

#________________________
# Plot the smooth terms: 
# "How has NDVI changed through time in each dryland class, after accounting for seasonal, spatial, species-level and individual-level variation?"

sm <- smooth_estimates(m_ndvi_trend)
yr0 <- mean(resident100km$year)
(plotLims <- quantile(resident100km$year_c, probs = c(0.025, 0.975), na.rm = TRUE))

plotDf <- sm %>%
  filter(grepl("^s\\(year_c\\):wmodeDrylands_4cat", .smooth)) %>%
  mutate(year = .data$year_c + yr0,
         dryland = factor(sub("^s\\(year_c\\):wmodeDrylands_4cat", "", .smooth),
                          levels = c("Arid", "Semiarid", "Dry subhumid", "Non-dryland")))

ggplot(plotDf, aes(year, .estimate, colour = dryland, fill = dryland)) +
  geom_ribbon(aes(ymin = .estimate - 2 * .se,
                  ymax = .estimate + 2 * .se),
              alpha = 0.2,
              colour = NA) +
  geom_line(linewidth = 1.2) +
  scale_x_continuous(breaks = seq(min(resident100km$year), max(resident100km$year), 2)) +
  scale_color_manual(values=c("#DAA520","#795C34", "#7E8C54", "#008000"))+
  scale_fill_manual(values=c("#DAA520","#795C34", "#7E8C54", "#008000")) +
  facet_grid(~ dryland) + # dry per row
  # facet_grid(dryland ~ .) + # dry per col
  labs(x = "Year",
       y = "Partial effect of year on NDVI",
       colour = "Dryland type",
       fill = "Dryland type") +
  theme_classic()

ggsave(paste0(outputPath, "/longTerm_NDVItrends_perDryland_facet.png"), width = 20, height = 5, dpi = 300)


#______________________________________________________________
### Separate long-term, seasonal and anomaly signal from gam
#______________________________________________________________

# NDVI decomposition based on above GAM

# NDVIobserved = 
#   NDVI_baseline (where the animal lives) + 
#   NDVI_longterm (directional change over 20 years) +
#   NDVI_seasonal (predictable annual cycle) +
#   NDVI_anomaly (unusual departure) + 
#   error

# First compute the chronic individual NDVI baseline (average productivity in their area)
resident100km <- resident100km %>%
  group_by(animal_id) %>%
  mutate(ndvi_between= mean(wmNDVI, na.rm=TRUE)) %>%
  ungroup()

# Extract seasonal only and long-term only trends:
gam_terms <- predict(m_ndvi_trend, newdata = resident100km, type="terms")
colnames(gam_terms)

# LONG-TERM trend: only the predicted trend estimated by s(year_c, dryland)
# to answer "do vultures respond to gradual productivity change over 20 years?"
resident100km$ndvi_longterm <- rowSums(gam_terms[, grep("year_c.*wmodeDrylands", colnames(gam_terms))])

# Relative SEASONAL effect
# "At this location, is this month normally a high-productivity or low-productivity month? relative to the annual cycle" (it ignores the individual context)
resident100km$ndvi_seasonalRel <- gam_terms[, "te(month,UDwMeanLongitude,UDwMeanLatitude)"]

# SEASONAL deviation from the individual's chronic mean
# to answer: "do vultures respond to seasonal departures from their individual own baseline productivity?"
seasonalExp <- predict(m_ndvi_trend,  # seasonal signal including all other terms (drylands, year etc)
  newdata = resident100km %>% mutate(year_c = 0),
  exclude = c("s(animal_id)","s(year_c,species)"))
resident100km$ndvi_seasonal_dev <- seasonalExp - resident100km$ndvi_between  # NDVI_between is the long-term individual NDVI mean

# Full GAM prediction for that location, time of the year, species and dryland type
resident100km$ndvi_expected_gam <- as.numeric(predict(m_ndvi_trend, newdata = resident100km))

# ANOMALY: observed NDVI - full prediction
# how far this specific month/year/individual deviates from its own full expectation
# to answer "do vultures respond to any deviation/unusual departure from the expected pattern"
resident100km$ndvi_anomaly <- resident100km$wmNDVI - resident100km$ndvi_expected_gam


# Scale all variables:
resident100km <- resident100km %>%
  mutate(
    ndvi_between_z = as.numeric(scale(ndvi_between)),
    ndvi_longterm_z = as.numeric(scale(ndvi_longterm)),
    ndvi_seasonalRel_z = as.numeric(scale(ndvi_seasonalRel)),
    ndvi_seasonal_dev_z = as.numeric(scale(ndvi_seasonal_dev)),
    ndvi_anomaly_z = as.numeric(scale(ndvi_anomaly)))

# test collinearity among those variables
cor(
  resident100km %>% 
    select(
      ndvi_between_z,
      ndvi_longterm_z,
      #ndvi_seasonalRel_z,
      ndvi_seasonal_dev_z,
      ndvi_anomaly_z
    ),
  use="complete.obs"
)


#______________________________________________________________________
### Section 1 on state: Model UD size as a function of NDVI components
#______________________________________________________________________

# Use a smooth for seasonal pattern
ggplot(resident100km,
      aes(ndvi_seasonal_dev, log(UDsizeKm2))) +
     geom_point(alpha=0.05) +
     geom_smooth(method="loess")


mod_UDsize <- bam(
  log(UDsizeKm2) ~
    ndvi_between_z +
    ndvi_longterm_z +
    s(ndvi_seasonal_dev_z, k = 5) +
    ndvi_anomaly_z +
    wmodeDrylands_4cat +
    species +
    locsPerMonth_logz +
    s(animal_id, bs = "re"),
  data     = resident100km,
  method   = "fREML",
  family = gaussian(),
  discrete = TRUE)

summary(mod_UDsize)
plot(mod_UDsize, page = 1)

mod_UDsize_month <- bam(
  log(UDsizeKm2) ~
    ndvi_between_z +
    ndvi_longterm_z +
    s(month, by = hemisphere, bs = 'cc', k = 12) +
    s(ndvi_seasonal_dev_z, k = 8) +
    ndvi_anomaly_z +
    wmodeDrylands_4cat +
    species +
    locsPerMonth_logz +
    s(animal_id, bs = "re"),
  data     = resident100km,
  method   = "fREML",
  family   = gaussian(),
  discrete = TRUE,
  knots     = list(month = c(0.5, 12.5)))

summary(mod_UDsize_month)
k.check(mod_UDsize_month)
plot.gam(mod_UDsize_month, page = 1)
AIC(mod_UDsize_month)
AIC(mod_UDsize)

#_____________
# Model plots

# Panel A: parametric NDVI effects (between, long-term, anomaly)
coefs <- broom.mixed::tidy(mod_UDsize_month, parametric = TRUE, conf.int = TRUE) %>%
  filter(term %in% c("ndvi_between_z", "ndvi_longterm_z", "ndvi_anomaly_z")) %>%
  mutate(term = recode(term,
                       "ndvi_between_z"   = "Chronic NDVI\n(individual baseline)",
                       "ndvi_longterm_z"  = "20-yr trend position\n(year_c × dryland)",
                       "ndvi_anomaly_z"   = "Anomalous NDVI\n(full-model residual)"))

panA <- ggplot(coefs, aes(x = estimate, y = term, xmin = conf.low, xmax = conf.high)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(height = 0.2, linewidth = 0.8) +
  geom_point(size = 3.5, colour = "#2e6b8a") +
  labs(title = "NDVI component effects on log(UD size)", x = "Estimate", y = NULL) +
  theme_bw(base_size = 12)

# Panel B: smooth seasonal-deviation effect
panB <- draw(smooth_estimates(mod_UDsize_month, select = "s(ndvi_seasonal_dev_z)")) +
  labs(title = "Seasonal NDVI deviation effect", x = "Seasonal deviation (z)", y = "Partial effect") +
  theme_bw(base_size = 11)

# Panel C: hemisphere-specific month smooth
panC <- draw(smooth_estimates(mod_UDsize_month,
                              select = c("s(month):hemisphereNorthern", "s(month):hemisphereSouthern"))) +
  labs(title = "Residual seasonal (calendar) effect by hemisphere") +
  theme_bw(base_size = 11)

# Panel D: predicted UD size by dryland category
predD <- ggpredict(mod_UDsize_month, terms = "wmodeDrylands_4cat")
panD <- ggplot(predD, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high)) +
  geom_pointrange(colour = "#e67e22", linewidth = 0.9) +
  labs(title = "Predicted log(UD size) by dryland category", x = NULL, y = "Predicted log(UD size)") +
  theme_bw(base_size = 11)

(panA | panD) / (panB | panC)
ggsave(paste0(outputPath, "/UDsize_NDVIdecomposition_model.png"), width = 13, height = 9, dpi = 300)



#______________________________________________________________________
### Section 2 on reactivity: 
### Model UD change as a function of monthly change in NDVI components
#______________________________________________________________________

#________________________________
# Compute month to month response

# We calculate lag 1 change in UD and in both the seasonal (predictable change) and anomaly (unusual change) components
# For non-consecutive months these change variables are set to NA

resident100km <- resident100km %>%
  arrange(animal_id, date) %>%
  group_by(animal_id) %>%
  mutate(
    # needed to make sure changes at month differences > 1 are set to NA
    months_since_prev = as.numeric((year - lag(year)) * 12 + (month - lag(month))),
    # UD change
    log_UD_size    = log(UDsizeKm2),
    UD_logChange   = ifelse(months_since_prev == 1, log(UDsizeKm2 / lag(UDsizeKm2)), NA_real_),
    UD_absChange   = ifelse(months_since_prev == 1, UDsizeKm2 - lag(UDsizeKm2), NA_real_),
    # Predictable NDVI change. Here we use seasonalRel, which is already centered at 0 (doind this difference on seasonal_dev would give the same, as the between and all other components are constant between months and would difference themselves out)
    ndvi_seasonalRel_lag = lag(ndvi_seasonalRel),
    ndvi_seasonal_change = ifelse(months_since_prev == 1,
                                  ndvi_seasonalRel - ndvi_seasonalRel_lag, NA_real_),
    # Unpredictable NDVI change.
    ndvi_anomaly_lag    = lag(ndvi_anomaly),
    ndvi_anomaly_change = ifelse(months_since_prev == 1,
                                 ndvi_anomaly - ndvi_anomaly_lag, NA_real_)) %>%
  ungroup() %>%
  # scale both
  mutate(ndvi_seasonal_change_z = as.numeric(scale(ndvi_seasonal_change)),
         ndvi_anomaly_change_z  = as.numeric(scale(ndvi_anomaly_change)))



# Filter
resident100km_s2 <- filter(resident100km, !is.na(UD_logChange))
(spToDrop <- names(which(table(as.character(resident100km_s2$ssp_dry)) < 20)))
resident100km_s2 <- resident100km_s2 %>% dplyr::filter(! ssp_dry %in% spToDrop)

# remove species that confound the species*dryland interaction (those that only occurr in one dryland type)
confounded_spp <- resident100km_s2 %>%
  distinct(species, wmodeDrylands_bin) %>%
  count(species) %>%
  filter(n < 2) %>%
  pull(species)

resident100km_s2_clean <- resident100km_s2 %>% filter(!species %in% confounded_spp)


#_________________________________________________________________
# Test reactive response to anomalies at different temporal scales

# Beyond seasonal change, how quickly do vultures react to anomaly conditions?
# does UD change respond to this month's anomaly, a past month's anomaly, or a sustained run of anomalous conditions? 

resident100km_s2_clean <- resident100km_s2_clean %>%
  arrange(animal_id, date) %>%
  group_by(animal_id) %>%
  mutate(
    # Lagged anomaly-change (concurrent already built as ndvi_anomaly_change_z)
    ndvi_anomaly_lag2  = ifelse(months_since_prev == 1, lag(ndvi_anomaly, 2),  NA_real_),
    ndvi_anomaly_lag3  = ifelse(months_since_prev == 1, lag(ndvi_anomaly, 3),  NA_real_),
    ndvi_anomaly_lag6  = ifelse(months_since_prev == 1, lag(ndvi_anomaly, 6),  NA_real_),
    
    # Sustained/cumulative anomaly: how anomalous has it been recently, on average?
    ndvi_anomaly_roll3    = slider::slide_dbl(ndvi_anomaly, mean, .before = 2, .complete = TRUE),
    ndvi_anomaly_roll6    = slider::slide_dbl(ndvi_anomaly, mean, .before = 5, .complete = TRUE),
    # Is the current anomaly deviating from the recent run of anomalies? (acceleration of anomaly)
    ndvi_anomaly_cumdev3  = ndvi_anomaly - ndvi_anomaly_roll3,
    ndvi_anomaly_cumdev6  = ndvi_anomaly - ndvi_anomaly_roll6
  ) %>%
  ungroup() %>%
  mutate(
    ndvi_anomaly_lag2_z    = as.numeric(scale(ndvi_anomaly_lag2)),
    ndvi_anomaly_lag3_z    = as.numeric(scale(ndvi_anomaly_lag3)),
    ndvi_anomaly_lag6_z    = as.numeric(scale(ndvi_anomaly_lag6)),
    ndvi_anomaly_roll3_z   = as.numeric(scale(ndvi_anomaly_roll3)),
    ndvi_anomaly_roll6_z   = as.numeric(scale(ndvi_anomaly_roll6)),
    ndvi_anomaly_cumdev3_z = as.numeric(scale(ndvi_anomaly_cumdev3)),
    ndvi_anomaly_cumdev6_z = as.numeric(scale(ndvi_anomaly_cumdev6))
  )

# Common subset of same rows across all models
df_common_scale <- resident100km_s2_clean %>%
  filter(!is.na(ndvi_anomaly_change_z), !is.na(ndvi_anomaly_lag6_z),
         !is.na(ndvi_anomaly_roll6_z), !species %in% confounded_spp)

make_scale_model <- function(ndvi_var, df) {
  df_sub <- df[!is.na(df[[ndvi_var]]), ]
  f <- as.formula(paste0(
    "UD_logChange ~ species + wmodeDrylands_bin + ",
    "s(ndvi_seasonal_change_z, by = interaction(species, wmodeDrylands_bin), k = 5) + ",
    "s(", ndvi_var, ", by = interaction(species, wmodeDrylands_bin), k = 5) + ",
    "s(month, by = hemisphere, bs = 'cc', k = 12) + ",
    "s(locsPerMonth_logz, by = species, k = 5) + ",
    "s(UDwMeanLongitude, UDwMeanLatitude, by = species, k = 50)"
  ))
  bam(f, data = df_sub, family = gaussian(), method = "fREML",
      discrete = TRUE, nthreads = 4, knots = list(month = c(0.5, 12.5)))
}

scale_models <- list(
  "Concurrent anomaly change"    = make_scale_model("ndvi_anomaly_change_z", df_common_scale),
  "Lag-2 anomaly"                = make_scale_model("ndvi_anomaly_lag2_z",   df_common_scale),
  "Lag-3 anomaly"                = make_scale_model("ndvi_anomaly_lag3_z",   df_common_scale),
  "Lag-6 anomaly"                = make_scale_model("ndvi_anomaly_lag6_z",   df_common_scale),
  "3-month rolling anomaly"      = make_scale_model("ndvi_anomaly_roll3_z",  df_common_scale),
  "6-month rolling anomaly"      = make_scale_model("ndvi_anomaly_roll6_z",  df_common_scale),
  "3-month cumulative deviation" = make_scale_model("ndvi_anomaly_cumdev3_z",df_common_scale),
  "6-month cumulative deviation" = make_scale_model("ndvi_anomaly_cumdev6_z",df_common_scale)
)

aic_scale_comparison <- data.frame(
  model = names(scale_models),
  AIC   = sapply(scale_models, AIC)
) %>%
  mutate(delta_AIC = AIC - min(AIC)) %>%
  arrange(delta_AIC)

print(aic_scale_comparison)
write.csv(aic_scale_comparison,
          paste0(outputPath, "/ndvi_anomaly_timescale_aic_comparison.csv"),
          row.names = FALSE)


#_________________________________________
# Explore output of best performing model


# Model
mod_UDchange <- bam(
  UD_logChange ~
    species + wmodeDrylands_bin +
    s(ndvi_seasonal_change_z, by = interaction(species, wmodeDrylands_bin), k = 5) +
    s(ndvi_anomaly_cumdev6_z,  by = interaction(species, wmodeDrylands_bin), k = 5) +
    s(month, by = hemisphere, bs = "cc", k = 12) +
    s(locsPerMonth_logz, by = species, k = 5) +
    s(UDwMeanLongitude, UDwMeanLatitude, by = species, k = 50),
  knots = list(month = c(0.5, 12.5)),
  data = resident100km_s2_clean, family = gaussian(),
  method = "fREML", 
  discrete = TRUE, 
  nthreads = 4)

summary(mod_UDchange)
k.check(mod_UDchange)
plot(mod_UDchange, page = 1)

fixedTable <- summary(mod_UDchange)$p.table
fixedTable[fixedTable[,"Pr(>|t|)"] < 0.05,]
smoothTable <- summary(mod_UDchange)["s.table"]$s.table
smoothTable[smoothTable[,"p-value"] < 0.05,]
smoothTable[smoothTable[,"edf"] < 2,]

# Less change in drylands. More change in some species specifically
smooth_labels <- sapply(mod_UDchange$smooth, function(s) s$label)
seas_labels <- grep("ndvi_seasonal_change_z", smooth_labels, fixed = TRUE, value = TRUE)
anom_labels <- grep("ndvi_anomaly_cumdev6_z", smooth_labels, fixed = TRUE, value = TRUE)
seas_smooths <- smooth_estimates(mod_UDchange, select = seas_labels)
anom_smooths <- smooth_estimates(mod_UDchange, select = anom_labels)
seas_smooths$species_lab <- gsub(".*\\)|\\..*", "", seas_smooths$'interaction(species, wmodeDrylands_bin)')
anom_smooths$species_lab <- gsub(".*\\)|\\..*", "", anom_smooths$'interaction(species, wmodeDrylands_bin)')
seas_smooths$dryland <- gsub(".*\\.", "", seas_smooths$'interaction(species, wmodeDrylands_bin)')
anom_smooths$dryland <- gsub(".*\\.", "", anom_smooths$'interaction(species, wmodeDrylands_bin)')


guild_seas <- seas_smooths %>% group_by(dryland, ndvi_seasonal_change_z) %>%
  summarise(est = mean(.estimate), se = mean(.se), .groups = "drop")
guild_anom <- anom_smooths %>% group_by(dryland, ndvi_anomaly_cumdev6_z) %>%
  summarise(est = mean(.estimate), se = mean(.se), .groups = "drop")

plotLims_seas <- quantile(resident100km_s2_clean$ndvi_seasonal_change_z, c(.025,.975), na.rm=TRUE)
plotLims_anom <- quantile(resident100km_s2_clean$ndvi_anomaly_cumdev6_z,  c(.025,.975), na.rm=TRUE)

col_dry <- c("Non-dryland"="#27ae60","Drylands"="#e67e22")

#_________
# Plots

# ## Guild level plots
# pan_guild_seas <- ggplot() +
#   geom_line(data = seas_smooths, aes(ndvi_seasonal_change_z, .estimate,
#                                      group = interaction(species_lab, dryland), colour = dryland), alpha = 0.3, linewidth = 0.5) +
#   geom_ribbon(data = guild_seas, aes(ndvi_seasonal_change_z, ymin = est-1.96*se, ymax = est+1.96*se, fill = dryland), alpha = 0.25) +
#   geom_line(data = guild_seas, aes(ndvi_seasonal_change_z, est, colour = dryland), linewidth = 1.3) +
#   geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
#   geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
#   scale_colour_manual(values = col_dry) + scale_fill_manual(values = col_dry) +
#   coord_cartesian(xlim = plotLims_seas, ylim = c(-1, 1.5)) +
#   labs(title = "Guild-level: predictable seasonal NDVI change",
#        x = "Seasonal NDVI change (z)", y = "log(UD change)", colour=NULL, fill=NULL) +
#   theme_bw(base_size = 11)
# 
# pan_guild_anom <- ggplot() +
#   geom_line(data = anom_smooths, aes(ndvi_anomaly_cumdev6_z, .estimate,
#                                      group = interaction(species_lab, dryland), colour = dryland), alpha = 0.3, linewidth = 0.5) +
#   geom_ribbon(data = guild_anom, aes(ndvi_anomaly_cumdev6_z, ymin = est-1.96*se, ymax = est+1.96*se, fill = dryland), alpha = 0.25) +
#   geom_line(data = guild_anom, aes(ndvi_anomaly_cumdev6_z, est, colour = dryland), linewidth = 1.3) +
#   geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
#   geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
#   scale_colour_manual(values = col_dry) + scale_fill_manual(values = col_dry) +
#   coord_cartesian(xlim = plotLims_anom, ylim = c(-1, 1.5)) +
#   labs(title = "Guild-level: anomalous NDVI change",
#        x = "Anomalous NDVI change (z)", y = "log(UD change)", colour=NULL, fill=NULL) +
#   theme_bw(base_size = 11)
# 
# pan_guild_seas | pan_guild_anom

## Species level plots
pan_species_seas <- ggplot(seas_smooths,
                           aes(ndvi_seasonal_change_z, .estimate,
                               ymin = .estimate-1.96*.se, ymax = .estimate+1.96*.se,
                               colour = dryland, fill = dryland)) +
  geom_ribbon(alpha = 0.2, colour = NA) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  facet_wrap(~species_lab, scales = "free_y", ncol = 4) +
  scale_colour_manual(values = col_dry) + scale_fill_manual(values = col_dry) +
  coord_cartesian(xlim = plotLims_seas, ylim = c(-4, 4)) +
  labs(title = "Species-level: seasonal NDVI change response",
       x = "Seasonal NDVI change (z)", y = "log(UD change)",
       caption = "Y-axis clipped to ±4 for display", colour=NULL, fill=NULL) +
  theme_bw(base_size = 9) + theme(legend.position = "bottom", strip.text = element_text(face="italic"))

pan_species_anom <- ggplot(anom_smooths,
                           aes(ndvi_anomaly_cumdev6_z, .estimate,
                               ymin = .estimate-1.96*.se, ymax = .estimate+1.96*.se,
                               colour = dryland, fill = dryland)) +
  geom_ribbon(alpha = 0.2, colour = NA) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  facet_wrap(~species_lab, scales = "free_y", ncol = 4) +
  scale_colour_manual(values = col_dry) + scale_fill_manual(values = col_dry) +
  coord_cartesian(xlim = plotLims_anom, ylim = c(-4, 4)) +
  labs(title = "Species-level: anomalous NDVI change response",
       x = "Cumulative anomalous NDVI change over 6 months (z)", y = "log(UD change)",
       caption = "Y-axis clipped to ±4 for display", colour=NULL, fill=NULL) +
  theme_bw(base_size = 9) + theme(legend.position = "bottom", strip.text = element_text(face="italic"))

(pan_species_seas | pan_species_anom)
ggsave(paste0(outputPath, "/UDchange_monthlyModel_seasonal_anomaly_perSpecies.png"), width = 18, height = 12, dpi = 300)

#______________
# Summary table

sm <- broom.mixed::tidy(mod_UDchange) %>%
  filter(grepl("ndvi_seasonal_change_z|ndvi_anomaly_change_z", term)) %>%
  mutate(
    predictor = ifelse(grepl("seasonal", term), "Seasonal change", "Anomaly change"),
    spp_dry   = sub(".*\\)", "", term)
  ) %>%
  separate(spp_dry, into = c("species", "dryland"), sep = "\\.", extra = "merge") %>%
  mutate(sig = case_when(p.value < 0.001 ~ "***", p.value < 0.01 ~ "**",
                         p.value < 0.05 ~ "*", TRUE ~ "")) %>%
  select(species, dryland, predictor, edf, sig) %>%
  pivot_wider(names_from = c(predictor, dryland), values_from = c(edf, sig))

write.csv(sm, paste0(outputPath, "/UDchange_species_significance_table.csv"), row.names = FALSE)



# #__________________
# # Livestock context
# # For each individual, we computed the mean livestock density across all monthly UDs, log-transformed to reduce right-skew, and standardised to unit variance for use as a between-individual contextual covariate.
# # to answer "Does livestock buffer the effect of NDVI change on UD change?" "does having more livestock around dampen the response to NDVI?" as a guild-level linear test
# 
# ind_livestock <- resident100km_s2_clean %>% # livestock per individual
#   group_by(animal_id) %>%
#   summarise(livestock_mean_raw = mean(wmLifestock, na.rm = TRUE),
#             livestock_mean_log = mean(log1p(wmLifestock), na.rm = TRUE),
#             .groups = "drop") %>%
#   mutate(livestock_logmean_z = as.numeric(scale(log1p(livestock_mean_raw))),
#          livestock_meanlog_z = as.numeric(scale(livestock_mean_log)))
# 
# resident100km_s2_clean <- resident100km_s2_clean %>% 
#   left_join(ind_livestock, by = "animal_id")
# 
# 
# # Check
# table(resident100km_s2_clean$wmLifestock == 0)
# hist(resident100km_s2_clean$wmLifestock, breaks = 80)
# hist(resident100km_s2_clean$livestock_logmean_z, breaks = 30)
# hist(resident100km_s2_clean$livestock_meanlog_z, breaks = 30)
# # Should give mean ~0, sd ~1, no NaN
# 
# # three-way, dryland-conditional buffering, as in your original follow-up model
# mod_livestock_dry <- lmer(
#   UD_logChange ~
#     ndvi_seasonal_change_z * livestock_logmean_z * wmodeDrylands_bin +
#     ndvi_anomaly_change_z  * livestock_logmean_z * wmodeDrylands_bin +
#     species +
#     (1 | animal_id),
#   data    = resident100km_s2_clean,
#   REML    = TRUE,
#   control = lmerControl(optimizer = "bobyqa")
# )
# 
# summary(mod_livestock_dry)
# 
# # the random effect has 0 variance in this and simpler model, meaning animal_id can be safely dropped
# mod_livestock_dry <- lm(
#   UD_logChange ~
#     ndvi_seasonal_change_z * livestock_logmean_z * wmodeDrylands_bin +
#     ndvi_anomaly_change_z  * livestock_logmean_z * wmodeDrylands_bin +
#     species,
#   data    = resident100km_s2_clean)
# 
# summary(mod_livestock_dry)
# 
# # Both three-way terms are significant, meaning dryland status genuinely modulates the buffering effect for both NDVI components.
# # During resular seasonal greening vultures decrease their UD, and they do it more when livestock density is higher?
# # When there is an anomalous green month, vultures would be expanding their range and would do it more if livestock density is higher.
# 
# seas_slopes_live <- emtrends(mod_livestock_dry, ~ wmodeDrylands_bin | livestock_logmean_z,
#                              var = "ndvi_seasonal_change_z",
#                              at = list(livestock_logmean_z = c(-1, 0, 1)))
# anom_slopes_live <- emtrends(mod_livestock_dry, ~ wmodeDrylands_bin | livestock_logmean_z,
#                              var = "ndvi_anomaly_change_z",
#                              at = list(livestock_logmean_z = c(-1, 0, 1)))
# 
# summary(seas_slopes_live, infer = c(TRUE, TRUE))
# summary(anom_slopes_live, infer = c(TRUE, TRUE))
# 
# #_______
# # Plots
# 
# pred_seas <- ggpredict(mod_livestock_dry,
#                        terms = c("ndvi_seasonal_change_z [all]", "livestock_logmean_z [-1,0,1]", "wmodeDrylands_bin"))
# pred_anom <- ggpredict(mod_livestock_dry,
#                        terms = c("ndvi_anomaly_change_z [all]",  "livestock_logmean_z [-1,0,1]", "wmodeDrylands_bin"))
# 
# pred_seas$group <- factor(pred_seas$group, levels = c("-1","0","1"),
#                           labels = c("Low (-1 SD)", "Mean", "High (+1 SD)"))
# pred_anom$group <- factor(pred_anom$group, levels = c("-1","0","1"),
#                           labels = c("Low (-1 SD)", "Mean", "High (+1 SD)"))
# 
# col_dry <- c("Non-dryland" = "#27ae60", "Drylands" = "#e67e22")
# 
# pan_liv_seas <- ggplot(pred_seas, aes(x = x, y = predicted,
#                                       colour = facet, fill = facet,
#                                       linetype = group, group = interaction(facet, group))) +
#   geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.12, colour = NA) +
#   geom_line(linewidth = 1.1) +
#   scale_colour_manual(values = col_dry) +
#   scale_fill_manual(values = col_dry) +
#   scale_linetype_manual(values = c("dotted","solid","dashed")) +
#   labs(title = "Livestock buffering: seasonal NDVI change",
#        x = "Seasonal NDVI change (z)", y = "Predicted log(UD change)",
#        colour = "Habitat", fill = "Habitat", linetype = "Livestock density") +
#   theme_bw(base_size = 11) + theme(legend.position = "bottom")
# 
# pan_liv_anom <- ggplot(pred_anom, aes(x = x, y = predicted,
#                                       colour = facet, fill = facet,
#                                       linetype = group, group = interaction(facet, group))) +
#   geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.12, colour = NA) +
#   geom_line(linewidth = 1.1) +
#   scale_colour_manual(values = col_dry) +
#   scale_fill_manual(values = col_dry) +
#   scale_linetype_manual(values = c("dotted","solid","dashed")) +
#   labs(title = "Livestock buffering: anomalous NDVI change",
#        x = "Anomalous NDVI change (z)", y = "Predicted log(UD change)",
#        colour = "Habitat", fill = "Habitat", linetype = "Livestock density") +
#   theme_bw(base_size = 11) + theme(legend.position = "bottom")
# 
# pan_liv_seas | pan_liv_anom
# ggsave(paste0(outputPath, "/livestock_buffering_bydryland.png"), width = 13, height = 6, dpi = 300)
