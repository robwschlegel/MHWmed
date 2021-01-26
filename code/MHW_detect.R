# code/MHW_detect.R
# This script detects MHWs in the temperature data


# Setup -------------------------------------------------------------------

# Necessary packages
library(tidyverse)
library(heatwaveR)
library(tidync)
library(doParallel); registerDoParallel(cores = 7)

# The data locations
# NB: These files are not hosted on GitHub as they are too large
# Contact Robert Schlegel to receive them: robert.schlegel@imev-mer.fr
med_SST_files <- dir("data/SST", pattern = ".nc", full.names = T, recursive = T)

# The lon/lat indexes
med_lat <- tidync(med_SST_files[1]) %>% 
  activate("D1") %>% 
  hyper_tibble()
med_lon <- tidync(med_SST_files[1]) %>% 
  activate("D2") %>% 
  hyper_tibble()

# The coords with SST data
# med_sea_coords <- tidync(med_SST_files[1]) %>% 
#   hyper_filter(time = time == 31536000) %>% 
#   hyper_tibble() %>% 
#   na.omit() %>% 
#   dplyr::select(lon, lat) %>% 
#   distinct()
# save(med_sea_coords, file = "metadata/med_sea_coords.RData")
med_sea_coords <- load("metadata/med_sea_coords.RData")


# MHW pipeline ------------------------------------------------------------

# Function to load a latitude subset of a single NetCDF file
# testers...
# file_name <- med_SST_files[100]
# lat_row <- 1
load_nc_sub <- function(file_name, lat_row){
  SST_sub <- tidync(file_name) %>% 
    hyper_filter(lat = lat == med_lat$lat[lat_row]) %>% #,
                 # lon = lon == med_lon$lon[1]) %>%
    hyper_tibble() %>% 
    mutate(t = as.Date(as.POSIXct(time, origin = "1981-01-01")),
           temp = round(analysed_sst - 273.15, 2)) %>% 
    dplyr::select(lon, lat, t, temp)
}

# 378 latitude pixels, 1305 longitude
MHW_pipeline <- function(lat_row){
  
  # Begin
  lat_row_pad <- str_pad(lat_row, width = 3, pad = "0", side = "left")
  print(paste("Began run", lat_row_pad, "at", Sys.time()))
  
  # Load data
  # system.time(
  SST_prep <- map_df(.x = med_SST_files, .f = load_nc_sub, lat_row = lat_row)
  # ) # 60 seconds for 1 full lat slice
  
  # Calculate MHWs
  # system.time(
  MHW_res <- SST_prep %>%
    # filter(lat == -63.375, lon == 0.125) %>% # tester...
    group_by(lon, lat) %>%
    nest() %>% 
    mutate(clim = purrr::map(data, ts2clm, climatologyPeriod = c("1982-01-01", "2011-12-31")),
           event = purrr::map(clim, detect_event), 
           cat = purrr::map(event, category, climatology = T, season = "peak", S = F)) %>%
    select(-data, -clim)
  # ) # 55 seconds for one full lat slice
  
  # Finish
  saveRDS(MHW_res, paste0("data/MHW/MHW_calc_", lat_row_pad,".Rds"))
  rm(SST_prep, MCS_res); gc()
  # print(paste("Completed run",lat_row_pad,"at",Sys.time()))
}

# Run it
plyr::l_ply(seq_len(nrow(med_lat)), MHW_pipeline, .parallel = T) # ~8 hours on 7 cores


# Extract pixel -----------------------------------------------------------

# If there is a desire to look at a single pixel, that can be done here

# Load the monthly results to peruse pixel based results
load("data/MHW_cat_pixel_monthly.RData")

# Looking at the above dataframe we can see that the following pixel had the highest max int.
lat_round <- round(med_lat$lat, 5)
single_pixel <- map_df(.x = med_SST_files, .f = load_nc_sub, lat_row = which(lat_round == 35.72916)) %>% 
  mutate(lat = round(lat, 5),
         lon = round(lon, 7)) %>%
  filter(lon == -5.1458511)
gc()


# We then run the MHW stats for this one pixel
MHW_pixel <- single_pixel %>%
  # filter(lat == -63.375, lon == 0.125) %>% # tester...
  group_by(lon, lat) %>%
  nest() %>% 
  mutate(clim = purrr::map(data, ts2clm, climatologyPeriod = c("1982-01-01", "2011-12-31")),
         event = purrr::map(clim, detect_event), 
         cat = purrr::map(event, category, climatology = T, season = "peak")) %>%
  select(-data, -clim)

# Extract the four different data.frames
event_clim <- MHW_pixel %>% 
  dplyr::select(-cat) %>% 
  unnest(event) %>% 
  filter(row_number() %% 2 == 1) %>% 
  unnest(event) %>% 
  ungroup()
event_event <- MHW_pixel %>% 
  dplyr::select(-cat) %>% 
  unnest(event) %>% 
  filter(row_number() %% 2 == 0) %>% 
  unnest(event) %>% 
  ungroup()
cat_clim <- MHW_pixel %>% 
  dplyr::select(-event) %>% 
  unnest(cat) %>% 
  filter(row_number() %% 2 == 1) %>% 
  unnest(cat) %>% 
  ungroup()
cat_event <- MHW_pixel %>% 
  dplyr::select(-event) %>% 
  unnest(cat) %>% 
  filter(row_number() %% 2 == 0) %>% 
  unnest(cat) %>% 
  ungroup()

# Combine same shaped dataframes
all_clim <- left_join(event_clim, cat_clim, by = c("lon", "lat", "t", "event_no"))
all_event <- left_join(event_event, cat_event, 
                       by = c("lon", "lat", "event_no", "duration", 
                              "intensity_max" = "i_max", "date_peak" = "peak_date")) 

# Save extracts as .csv files for sharing across languages
write_csv(all_clim, paste0("data/extract/MHW_clim_",all_clim$lon[1],"_",all_clim$lat[1],".csv"))
write_csv(all_event, paste0("data/extract/MHW_event_",all_event$lon[1],"_",all_event$lat[1],".csv"))

# Look at the one crazy event
event_line(detect_event(ts2clm(single_pixel, climatologyPeriod = c("1982-01-01", "2011-12-31"))))
ggsave(paste0("data/extract/MHW_plot_",all_clim$lon[1],"_",all_clim$lat[1],".png"))

