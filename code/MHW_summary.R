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
    na.omit() %>% 
    mutate(category = factor(category))
  
  # Clean up and exit
  rm(res_full); gc()
  return(cat_clim)
}

# Function for calculating the earliest date of occurrence of highest event
cat_pixel_calc <- function(file_name){
  
  # Load data
  MHW_cat <- load_cat_daily(file_name)
  
  # The sum of intensities per pixel for the year
  MHW_intensity <- MHW_cat %>% 
    group_by(lon, lat) %>% 
    summarise(intensity_sum = sum(intensity), .groups = "drop")
  
  # The earliest date of the highest category of event + the sum of intensities
  MHW_cat_pixel <- MHW_cat %>% 
    group_by(lon, lat) %>%
    filter(as.integer(category) == max(as.integer(category))) %>% 
    filter(t == min(t)) %>% 
    unique() %>%
    left_join(MHW_intensity, by = c("lon", "lat"))
  
}


# Pipeline to calculate the summary stats from each file
MHW_summary_pipeline <- function(file_name){
  

  
  # Summarise the count of how many of each category of events were experienced in each pixel
  MHW_cat_count <- MHW_cat %>% 
    group_by(lon, lat, event_no) %>% 
    summarise(max_cat = max(as.integer(category)), .groups = "drop") %>% 
    mutate(max_cat = factor(max_cat, levels = c(1:4),  labels = levels(MHW_cat$category))) %>%
    dplyr::select(-event_no) %>%
    group_by(lon, lat) %>% 
    table() %>% 
    as.data.frame() %>% 
    pivot_wider(values_from = Freq, names_from = max_cat) %>% 
    mutate(lon = as.numeric(as.character(lon)),
           lat = as.numeric(as.character(lat)))
  
  # Complete dates by categories data.frame
  full_grid <- expand_grid(t = seq(as.Date(paste0("1982-01-01")), as.Date("2019-12-31"), by = "day"), 
                           category = as.factor(levels(MHW_cat$category))) %>% 
    mutate(category = factor(category, levels = levels(MHW_cat$category)))
  
  # The count and area of the first time the largest category pixel occurs at each pixel and the cumulative values
  MHW_cat_first <- MHW_cat_pixel %>%
    # dplyr::select(lon, lat, t, category) %>% 
    group_by(t, category) %>%
    summarise(first_n = n(), .groups = "drop") %>% 
    right_join(full_grid, by = c("t", "category")) %>%
    mutate(first_n = ifelse(is.na(first_n), 0, first_n),
           first_n_prop = round(first_n/nrow(med_sea_coords), 4)) %>% 
    arrange(t, category) %>% 
    group_by(category) %>%
    mutate(first_n_cum = cumsum(first_n),
           first_n_cum_prop = round(first_n_cum/nrow(med_sea_coords), 4)) %>% 
    ungroup()
  
  # Clean up and exit
  rm(MHW_cat, MHW_intensity); gc()
  return(MHW_cat_pixel)
  
}
