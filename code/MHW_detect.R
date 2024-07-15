# code/MHW_detect.R
# This script detects MHWs in the temperature data

# NB: This script may be run via: 
# source("code/MHW_detect.R")


# Setup -------------------------------------------------------------------

# Necessary packages
library(tidyverse)
library(heatwaveR)
library(tidync)
library(doParallel)

# Set cores
registerDoParallel(cores = detectCores()-1)

# The data locations
# NB: These files are not hosted on GitHub as they are too large
# Contact Robert Schlegel to receive them: robwschlegel@gmail.com
# Or run the code at 'code/SST_downloads.R'
med_SST_files <- dir("~/data/MED_REP", pattern = ".nc", full.names = T)

# The lon/lat indexes
med_lat <- tidync(med_SST_files[1]) |> 
  activate("D1") |> hyper_tibble() |> 
  rename(lat = latitude)
med_lon <- tidync(med_SST_files[1]) |>
  activate("D2") |> hyper_tibble() |> 
  rename(lon = longitude)

# The coords with SST data
# med_sea_coords <- tidync(med_SST_files[1]) |>
#   hyper_filter(time = time == 31536000) |>
#   hyper_tibble() |> na.omit() |>
#   rename(lon = longitude, lat = latitude) |> 
#   dplyr::select(lon, lat) |> distinct()
# save(med_sea_coords, file = "metadata/med_sea_coords.RData")
load("metadata/med_sea_coords.RData")


# MHW pipeline ------------------------------------------------------------

# Function to load a latitude subset of a single NetCDF file
# testers...
# file_name <- med_SST_files[23]; lat_row <- 1
load_nc_sub <- function(file_name, lat_row){
  SST_sub <- tidync(file_name) |>
    hyper_filter(latitude = latitude == med_lat$lat[lat_row]) |>
    hyper_tibble() |>
    mutate(t = as.Date(as.POSIXct(time, origin = "1981-01-01")),
           temp = round(analysed_sst - 273.15, 2)) |>
    dplyr::rename(lon = longitude, lat = latitude) |> 
    dplyr::select(lon, lat, t, temp)
  # rm(file_name, lat_row, SST_sub); gc()
}

# 318 latitude pixels, 1089 longitude
# lat_row <- 92
# base_years <- c(1982, 2011)
MHW_pipeline <- function(lat_row, base_years){
  
  # Begin
  lat_row_pad <- str_pad(lat_row, width = 3, pad = "0", side = "left")
  print(paste("Began run", lat_row_pad, "at", round(Sys.time())))
  
  # Create baseline period and file name
  base_line <- c(paste0(base_years[1],"-01-01"), paste0(base_years[2],"-12-31"))
  MHW_file_name <- paste0("data/MHW/MHW_calc_",lat_row_pad,"_",base_years[1],"_",base_years[2],".Rds")
  
  if(!file.exists(MHW_file_name)){
    
    # Load data
    # system.time(
    SST_prep <- map_df(.x = med_SST_files, .f = load_nc_sub, lat_row = lat_row)
    # ) # 5 seconds for 1 full lat slice
    
    # Calculate MHWs
    # system.time(
    MHW_res <- SST_prep |>
      # filter(lat == -63.375, lon == 0.125) |># tester...
      group_by(lon, lat) |>
      nest() |>
      mutate(clim = purrr::map(data, ts2clm, climatologyPeriod = base_line),
             event = purrr::map(clim, detect_event), 
             cat = purrr::map(event, category, climatology = T, season = "peak", S = F)) |>
      select(-data, -clim)
    # ) # ~2 minutes for one full lat slice
    
    # Finish
    saveRDS(MHW_res, MHW_file_name)
    rm(lat_row_pad, SST_prep, base_line, MHW_res); gc()
    # print(paste("Completed run",lat_row_pad,"at",Sys.time()))
  }
  # return()
}

# Run it for 1982-2011 baseline
# system.time(MHW_pipeline(90, base_years = c(1982, 2011))) # ~ 2-3 minutes for one, depending on latitude
plyr::l_ply(seq_len(nrow(med_lat)), MHW_pipeline, .parallel = T, base_years = c(1982, 2011)) # ~160 minutes on 15 cores
Sys.sleep(10) # Sleep for ten seconds
plyr::l_ply(seq_len(nrow(med_lat)), MHW_pipeline, .parallel = T, base_years = c(1991, 2020)) # ~160 minutes on 15 cores

