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
med_SST_files <- dir("data", pattern = ".nc", full.names = T, recursive = T)

# The lon/lat indexes
med_lat <- tidync(med_SST_files[1]) %>% 
  activate("D1") %>% 
  hyper_tibble()
med_lon <- tidync(med_SST_files[1]) %>% 
  activate("D2") %>% 
  hyper_tibble()

# The 2019 data
# SST_2019 <- readMat("data/2019/L4_REP_SST_MED_2019.mat")


# MHW pipeline ------------------------------------------------------------

# Function to load a latitude subset of a single NetCDF file
# testers...
# file_name <- med_SST_files[100]
# lat_int <- 1
load_nc_sub <- function(file_name, lat_int){
  SST_sub <- tidync(file_name) %>% 
    hyper_filter(lat = lat == med_lat$lat[lat_int],
                 lon = lon == med_lon$lon[1]) %>%
    hyper_tibble() %>% 
    mutate(t = as.Date(as.POSIXct(time, origin = "1981-01-01")),
           temp = round(analysed_sst - 273.15, 2)) %>% 
    dplyr::select(lon, lat, t, temp)
}

# Convenience function for parallel processing
detect_event_cat <- function(df){
  df_event <- detect_event(ts2clm(df, climatologyPeriod = c("1982-01-01", "2011-12-31")))
  df_cat <- category(df_event, climatology = T, season = "peak")
  res <- list(event = df_event, cat = df_cat)
  return(res)
}

# 378 latitude pixels, 1305 longitude
MHW_pipeline <- function(lat_int){
  
  # Load data
  SST_prep <- plyr::ldply(med_SST_files, load_nc_sub, lat_int = lat_int, .parallel = T, .paropts = c(.inorder = FALSE))
  
  # Calculate MHWs
  # Make calculations
  MCS_res <- SST_prep %>%
    # filter(lat == -63.375, lon == 0.125) %>% # tester...
    group_by(lon, lat) %>%
    nest() %>% 
    mutate(clim = purrr::map(data, ts2clm, climatologyPeriod = c("1982-01-01", "2011-12-31")),
           event = purrr::map(clim, detect_event), 
           cat = purrr::map(event, category, climatology = T, season = "peak")) %>%
    select(-data, -clim)
}

