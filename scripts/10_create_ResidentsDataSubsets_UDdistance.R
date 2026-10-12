# ---
# Title: Create Residents subsets
# Author: Anne Scharf, MPIAB
# Date: July 2026
# Description: Identify, for each individual, periods of residency based on the distance between the weighted centroids of all monthly UDs relative to the first UD
# ---

library(sf)
library(ggplot2)

## reading in data
df <- readRDS("./data_for_models/table_for_models_july2026.rds")
## adding date
df$date <- as.Date(paste0(df$yearMonth,"-01"),format="%Y-%m-%d")
## enrusing correct ordering
df <- df[order(df$animal_id,df$date),]
head(df)

## ---- distance to first centroid (needed for the criteria) ----------------
## transfrom into spatial object
dfst <- st_as_sf(df, coords=c("UDwMeanLongitude", "UDwMeanLatitude"))
st_crs(dfst) <- "EPSG:4326"
## calculating distance to first centroid per individual
dfst_L <- split(dfst,dfst$animal_id)
dist_to_1st_cntr <- unlist(lapply(dfst_L, function(x){
  sf::st_distance(x,x[1,])
}))
dfst$dist_to_1st_cntr_km <- dist_to_1st_cntr/1000
df$dist_to_1st_cntr_km <- dfst$dist_to_1st_cntr_km

## ---- annotate subsetting criteria ----------------------------------------
## monotonic month index: consecutive calendar months differ by exactly 1
month_index <- function(ym){
  as.integer(substr(ym, 1, 4)) * 12L + as.integer(substr(ym, 6, 7)) - 1L
}

## flag rows that belong to a run of >= 12 consecutive months for which the
## centroid stays within `thr` km of the individual's first centroid.
## A run is broken either by a missing month or by a month exceeding `thr`.
consec_within <- function(m, dist_km, thr, min_len = 12L){
  n <- length(m)
  e <- dist_km <= thr                                   # month is close enough
  ## row continues the previous run only if consecutive in time AND both close
  cont <- c(FALSE, (diff(m) == 1L) & e[-1] & e[-n])
  run_id  <- cumsum(!cont)                              # new id at each run start
  run_len <- ave(run_id, run_id, FUN = length)
  e & (run_len >= min_len)
}

## flag rows of individuals whose data within `thr` km of the first centroid
## spans at least `min_years` years. Gaps are allowed: both missing months and
## months further than `thr` km simply do not count towards the span.
span_years_within <- function(m, dist_km, thr, min_years = 5L){
  e <- dist_km <= thr                                   # month is close enough
  if (!any(e)) return(rep(FALSE, length(m)))
  me   <- m[e]
  span <- max(me) - min(me) + 1L                        # first->last close month, in months
  e & (span >= min_years * 12L)
}

## process each individual separately (rows already ordered by animal & date)
df_L <- split(df, df$animal_id)
df_L <- lapply(df_L, function(x){
  m <- month_index(x$yearMonth)          # sorted, unique months for this individual
  d <- diff(m)                           # 1 = consecutive month, >1 = gap

  ## 1. >=12 consecutive months within max 50 km of the first centroid
  x[["completeseason_50km"]]  <- consec_within(m, x$dist_to_1st_cntr_km, 50)

  ## 2. >=12 consecutive months within max 100 km of the first centroid
  x[["completeseason_100km"]] <- consec_within(m, x$dist_to_1st_cntr_km, 100)

  ## 3. >=5 years of data within max 50 km of the first centroid (gaps allowed)
  x[["resident_5yr_50km"]]    <- span_years_within(m, x$dist_to_1st_cntr_km, 50, min_years = 5L)

  ## 4. >=10 years of data within max 50 km of the first centroid (gaps allowed)
  x[["resident_10yr_50km"]]   <- span_years_within(m, x$dist_to_1st_cntr_km, 50, min_years = 10L)

  ## 5. >=5 years of data within max 50 km of the first centroid (gaps allowed)
  x[["resident_5yr_100km"]]    <- span_years_within(m, x$dist_to_1st_cntr_km, 100, min_years = 5L)
  
  ## 6. >=10 years of data within max 50 km of the first centroid (gaps allowed)
  x[["resident_10yr_100km"]]   <- span_years_within(m, x$dist_to_1st_cntr_km, 100, min_years = 10L)
  
  ## number of distinct calendar years available
  x[["tracking_years"]] <- length(unique(as.integer(substr(x$yearMonth, 1, 4))))

  ## number of gaps (a gap = one or more missing months between two present months)
  gaps <- d[d > 1L] - 1L                 # missing-month count for each gap
  x[["nb_gaps"]] <- length(gaps)

  ## longest gap duration in number of missing months (0 if no gaps)
  x[["max_gap_duration"]] <- if (length(gaps)) max(gaps) else 0L
  x
})
df <- do.call(rbind, df_L)
rownames(df) <- NULL
head(df)
## --------------------------------------------------------------------------
saveRDS(df, "./data_for_models/table_for_models_july2026_w_filter_columns.rds")


## ---- summary: completeseason_50km per species ----------------------------
## per species: number of individuals and number of individual-months that
## satisfy the completeseason_50km criterion
sub50 <- df[df$completeseason_50km == TRUE, ]
sp50   <- droplevels(factor(sub50$species))   # only species with qualifying data
summary_50km <- data.frame(
  species       = levels(sp50),
  n_individuals_50km = tapply(sub50$animal_id, sp50, function(z) length(unique(z))),
  n_months_50km      = as.integer(table(sp50)),
  row.names     = NULL
)
summary_50km <- summary_50km[order(-summary_50km$n_individuals_50km), ]
# print(summary_50km)
# write.csv(summary_50km, "./data_for_models/summary_completeseason_50km.csv", row.names = FALSE)
## --------------------------------------------------------------------------

## ---- summary: completeseason_100km per species ----------------------------
## per species: number of individuals and number of individual-months that
## satisfy the completeseason_100km criterion
sub100 <- df[df$completeseason_100km == TRUE, ]
sp100   <- droplevels(factor(sub100$species))   # only species with qualifying data
summary_100km <- data.frame(
  species       = levels(sp100),
  n_individuals_100km = tapply(sub100$animal_id, sp100, function(z) length(unique(z))),
  n_months_100km      = as.integer(table(sp100)),
  row.names     = NULL
)
summary_100km <- summary_100km[order(-summary_100km$n_individuals_100km), ]
# print(summary_100km)
# write.csv(summary_100km, "./data_for_models/summary_completeseason_100km.csv", row.names = FALSE)
## --------------------------------------------------------------------------
summary_seasonal_50_100km <- merge(summary_50km,summary_100km, by="species")
write.csv(summary_seasonal_50_100km, "./data_for_models/summary_completeseason_50_100km_july2026.csv", row.names = FALSE)

## ---- visualise data availability per criterion ---------------------------
## one data-availability plot per species (individuals on y-axis, time on x-axis),
## restricted to the rows that satisfy each criterion.
out_base <- "./exploration_plots/"
criteria <- list(
  list(col = "completeseason_50km",  folder = "max50km",
       desc = ">=12 consecutive months within max 50 km"),
  list(col = "completeseason_100km", folder = "max100km",
       desc = ">=12 consecutive months within max 100 km"),
  list(col = "resident_5yr_50km",    folder = "resident5yr_50km",
       desc = ">=5 years of data within max 50 km (gaps allowed)"),
  list(col = "resident_10yr_50km",   folder = "resident10yr_50km",
       desc = ">=10 years of data within max 50 km (gaps allowed)")
)

invisible(lapply(criteria, function(cr){
  out_dir <- file.path(out_base, cr$folder)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  sub   <- df[df[[cr$col]] == TRUE, ]
  sub_L <- split(sub, sub$species)
  sub_L <- sub_L[sapply(sub_L, nrow) > 0]   # keep only species that have data

  invisible(lapply(names(sub_L), function(sp){
    d <- sub_L[[sp]]
    ## order individuals on the y-axis by their first month of data
    ord <- names(sort(tapply(d$date, d$animal_id, min)))
    d$animal_id <- factor(d$animal_id, levels = ord)
    n_ind <- nlevels(d$animal_id)

    p <- ggplot(d, aes(x = date, y = animal_id)) +
      geom_point(shape = 15, size = 1.5, colour = "steelblue") +
      labs(title = sp, x = "Date", y = "Individual",
           subtitle = paste0(cr$desc,
                             "  |  ", n_ind, " individuals, ", nrow(d), " individual-months")) +
      theme_minimal(base_size = 9) +
      theme(plot.title = element_text(face = "italic"))

    ggsave(file.path(out_dir, paste0("data_availability_", gsub("[^A-Za-z0-9]+", "_", sp), ".png")),
           plot = p, width = 10, height = max(3, n_ind * 0.18), limitsize = FALSE)
  }))
}))
## --------------------------------------------------------------------------

## ---- variant plots: kept vs unused data ----------------------------------
## same plots as above, but for each qualifying individual ALSO show the data
## that are NOT used (do not satisfy the criterion) in a different colour.
## Written to separate "_with_unused" folders so the original plots are kept.
invisible(lapply(criteria, function(cr){
  out_dir <- file.path(out_base, paste0(cr$folder, "_with_unused"))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  qual_ids   <- unique(df$animal_id[df[[cr$col]]])   # individuals that qualify
  sub        <- df[df$animal_id %in% qual_ids, ]     # ALL of their data
  sub$used   <- factor(ifelse(sub[[cr$col]], "used", "not used"),
                       levels = c("used", "not used"))
  sub_L <- split(sub, sub$species)
  sub_L <- sub_L[sapply(sub_L, nrow) > 0]   # keep only species that have data

  invisible(lapply(names(sub_L), function(sp){
    d <- sub_L[[sp]]
    ## order individuals on the y-axis by their first month of data
    ord <- names(sort(tapply(d$date, d$animal_id, min)))
    d$animal_id <- factor(d$animal_id, levels = ord)
    n_ind  <- nlevels(d$animal_id)
    n_used <- sum(d$used == "used")

    p <- ggplot(d, aes(x = date, y = animal_id, colour = used)) +
      geom_point(shape = 15, size = 1.5) +
      scale_colour_manual(values = c("used" = "steelblue", "not used" = "grey70")) +
      labs(title = sp, x = "Date", y = "Individual", colour = NULL,
           subtitle = paste0(cr$desc, "  |  ", n_ind, " individuals, ",
                             n_used, " used / ", nrow(d), " total individual-months")) +
      theme_minimal(base_size = 9) +
      theme(plot.title = element_text(face = "italic"), legend.position = "top")

    ggsave(file.path(out_dir, paste0("data_availability_", gsub("[^A-Za-z0-9]+", "_", sp), ".png")),
           plot = p, width = 10, height = max(3, n_ind * 0.18), limitsize = FALSE)
  }))
}))
## --------------------------------------------------------------------------