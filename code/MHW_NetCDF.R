# code/MHW_NetCDF.R
# This script contains the code used to convert MHW results to NetCDF files


# Setup -------------------------------------------------------------------

# Needed packages
library(tidyverse)
library(tidync)
library(ncdf4)

# Original SST files for reference
med_SST_files <- dir("data/SST", pattern = ".nc", full.names = T, recursive = T)

# The lon/lat indexes
med_lat <- tidync(med_SST_files[1]) %>% 
  activate("D1") %>% 
  hyper_tibble()
med_lon <- tidync(med_SST_files[1]) %>% 
  activate("D2") %>% 
  hyper_tibble()
med_coords <- expand.grid(lat = as.vector(med_lat$lat), 
                          lon = as.vector(med_lon$lon)) %>% 
  data.frame()

# File locations
MHW_files <- dir("data/MHW", full.names = T)


# Functions ---------------------------------------------------------------

# Function for creating arrays from data.frames
# df <- filter(df_step, t == 5786)
df_acast <- function(df){
  
  # Ensure correct grid size
  med_coords_sub <- med_coords %>% 
    filter(lat == df$lat[1])
  
  # Round data for massive file size reduction
  df$temp <- round(df$temp, 2)
  
  # Force grid
  res <- df %>%
    right_join(med_coords_sub, by = c("lon", "lat")) %>% 
    arrange(lon, lat)
  
  # Create array
  res_array <- base::array(res$temp, dim = c(1305,1,1))
  dimnames(res_array) <- list(lon = unique(med_coords_sub$lon),
                              lat = unique(med_coords_sub$lat),
                              t = unique(na.omit(res$t)))
  return(res_array)
}

# Wrapper function for last step before data are entered into NetCDF files
# df <- event_clim
df_proc <- function(df){
  
  # Filter NA and convert dates to integer
  df_step <- df %>% 
    mutate(temp = ifelse(is.na(temp), NA, temp),
           t = as.integer(t)) %>% 
    na.omit()
  
  # Acast
  dfa <- df_step %>%
    mutate(t2 = t) %>% 
    group_by(t2) %>%
    nest() %>%
    mutate(data2 = purrr::map(data, OISST_acast)) %>%
    select(-data)
  
  # Final form
  dfa_temp <- abind(dfa$data2, along = 3, hier.names = T)
  # dimnames(dfa_temp)
  return(dfa_temp)
}

# Load a Rds file and save it as NetCDF
Rds_to_NetCDF <- function(file_name){

  # Load R file
  MHW_res <- readRDS(file_name)
  
  # Extract the four different data.frames
  event_clim <- MHW_res %>% 
    dplyr::select(-cat) %>% 
    unnest(event) %>% 
    filter(row_number() %% 2 == 1) %>% 
    unnest(event) %>% 
    ungroup()
  event_event <- MHW_res %>% 
    dplyr::select(-cat) %>% 
    unnest(event) %>% 
    filter(row_number() %% 2 == 0) %>% 
    unnest(event) %>% 
    ungroup()
  cat_clim <- MHW_res %>% 
    dplyr::select(-event) %>% 
    unnest(cat) %>% 
    filter(row_number() %% 2 == 1) %>% 
    unnest(cat) %>% 
    ungroup()
  cat_event <- MHW_res %>% 
    dplyr::select(-event) %>% 
    unnest(cat) %>% 
    filter(row_number() %% 2 == 0) %>% 
    unnest(cat) %>% 
    ungroup()
  
  # Get file attributes
  lon <- unique(event_clim$lon)
  lat <- unique(event_clim$lat)
  time <- as.integer(unique(event_clim$t))
  tunits <- "days since 1970-01-01"
  doy <- unique(event_clim$doy)
  eventno <- unique(event_event$event_no)
  
  # Length of each attribute
  nlon <- length(lon)
  nlat <- length(lat)
  ndoy <- length(doy)
  nt <- length(t)
  nen <- length(eventno)
  
  # Convert single columns in data.frames into arrays
  
  
  
}


# Convert files -----------------------------------------------------------


