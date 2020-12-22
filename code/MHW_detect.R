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

# The 2019 data
# SST_2019 <- readMat("data/2019/L4_REP_SST_MED_2019.mat")


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
           cat = purrr::map(event, category, climatology = T, season = "peak")) %>%
    select(-data, -clim)
  # ) # 55 seconds for one full lat slice
  
  # Finish
  saveRDS(MHW_res, paste0("data/MHW/MHW_calc_", lat_row_pad,".Rds"))
  rm(SST_prep, MCS_res); gc()
  print(paste("Completed run",lon_row_pad,"at",Sys.time()))
}

# Run it
plyr::l_ply(seq_len(nrow(med_lat)), MHW_pipeline, .parallel = T)

