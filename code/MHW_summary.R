# code/MHW_summary.R
# This script creates annual and total summaries from the MHW Med results


# Setup -------------------------------------------------------------------

# The needed packages
library(tidyverse)
library(doParallel); registerDoParallel(cores = 7)

# The file location
res_files <- dir("data/MHW", full.names = T)

# The coords with SST
load("metadata/med_sea_coords.RData")


# Functions ---------------------------------------------------------------

# Function for loading only the daily cat results
load_cat_daily <- function(file_name){
  
  # Load data
  res_full <- readRDS(file_name)
  
  # Extract only daily cat values
  cat_clim <- res_full %>% 
    dplyr::select(-event) %>% 
    unnest(cat) %>% 
    filter(row_number() %% 2 == 1) %>% 
    unnest(cat) %>% 
    ungroup() %>% 
    mutate(category = factor(category))
  
  # Clean up and exit
  rm(res_full); gc()
  return(cat_clim)
}

# Function for finding the first date of the highest category MHW per pixel
max_event_date <- function(df){
  df %>% 
    filter(as.integer(category) == max(as.integer(category))) %>% 
    filter(t == min(t))
}

# Function to calculate the max categories per pixel
max_pixel <- function(file_name){
  
  # Load data
  MHW_cat <- load_cat_daily(file_name)
  
  # Sum of intensity per pixel
  MHW_intensity <- MHW_cat %>% 
    group_by(lon, lat) %>% 
    summarise(intensity_sum = sum(intensity), .groups = "drop")
  
  # The date of first occurrence of the highest  event
  MHW_cat_pixel <- MHW_cat %>% 
    group_by(lon, lat) %>%
    filter(as.integer(category) == max(as.integer(category), na.rm = T)) %>% 
    filter(t == min(t)) %>% 
    unique() %>%
    left_join(MHW_intensity, by = c("lon", "lat"))
  
  # Clean up and exit
}


MHW_summary_pipeline <- function(){
  
  # Load cat daily data
  MHW_cat <- plyr::ldply(res_files, load_cat_daily, .parallel = T)
  
}
