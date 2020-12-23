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

# Complete annual dates by categories
full_annual_grid <- expand_grid(year = seq(1982:2019), 
                                category = as.factor(c("I Moderate", "II Strong", "III Severe", "IV Extreme")))

# Complete daily dates by categories
full_daily_grid <- expand_grid(t = seq(as.Date(paste0("1982-01-01")), as.Date("2019-12-31"), by = "day"), 
                               category = as.factor(c("I Moderate", "II Strong", "III Severe", "IV Extreme")))

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
    mutate(category = factor(category),
           year = lubridate::year(t))
  
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
    group_by(lon, lat, year) %>% 
    summarise(intensity_sum = sum(intensity), .groups = "drop")
  
  # The count of the highest category of each event in each pixel
  MHW_cat_count <- MHW_cat %>% 
    group_by(lon, lat, year, event_no) %>% 
    summarise(category = max(as.integer(category), na.rm = T), .groups = "drop") %>% 
    mutate(category = factor(category, levels = c(1:4),
                            labels = c("I Moderate", "II Strong", "III Severe", "IV Extreme"))) %>%
    dplyr::select(-event_no) %>%
    group_by(lon, lat, year, category) %>%
    summarise(count = n(), .groups = "drop") %>% 
    pivot_wider(values_from = count, names_from = category) %>% 
    replace(is.na(.), 0)
  if(!"III Severe" %in% colnames(MHW_cat_count)){
    MHW_cat_count <- MHW_cat_count %>% 
      cbind(data.frame('III Severe' = 0))
    colnames(MHW_cat_count)[colnames(MHW_cat_count) == 'III.Severe'] <- 'III Severe'
  }
  if(!"IV Extreme" %in% colnames(MHW_cat_count)){
    MHW_cat_count <- MHW_cat_count %>% 
      cbind(data.frame('IV Extreme' = 0))
    colnames(MHW_cat_count)[colnames(MHW_cat_count) == 'IV.Extreme'] <- 'IV Extreme'
  }
  
  # The earliest date of the highest category of event
  MHW_cat_pixel <- MHW_cat %>% 
    group_by(lon, lat, year) %>%
    filter(as.integer(category) == max(as.integer(category), na.rm = T)) %>% 
    filter(t == min(t)) %>% 
    unique() %>%
    data.frame() %>% 
    left_join(MHW_intensity, by = c("lon", "lat", "year")) %>% 
    left_join(MHW_cat_count, by = c("lon", "lat", "year"))
  
  # Clean up and exit
  rm(MHW_cat, MHW_intensity, MHW_cat_count); gc()
  return(MHW_cat_pixel)
}

# Function to count how many of each category of events were experienced in each pixel
cat_count_calc <- function(file_name){
  
  # Load data
  MHW_cat <- load_cat_daily(file_name)
 

  
  # Clean up and exit
  rm(MHW_cat); gc()
  return(MHW_cat_count)
}

# Pipeline to calculate the summary stats from each file
MHW_summary_pipeline <- function(file_name){

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


# Summaries ---------------------------------------------------------------

# The first highest occurrence per pixel
# system.time(
MHW_cat_pixel_all <- plyr::ldply(res_files, cat_pixel_calc, .parallel = T)
# ) # 265 seconds on 7 cores
save(MHW_cat_pixel_all, file = "data/MHW_cat_pixel_all.RData")

# The total count per pixel
# system.time(
MHW_cat_count_all <- plyr::ldply(res_files, cat_count_calc, .parallel = T)
# ) # 293 seconds on 7 cores
save(MHW_cat_count_all, file = "data/MHW_cat_count_all.RData")

# Complete dates by categories data.frame
full_grid <- expand_grid(t = seq(as.Date(paste0("1982-01-01")), as.Date("2019-12-31"), by = "day"), 
                         category = as.factor(c("I Moderate", "II Strong", "III Severe", "IV Extreme"))) #%>% 
  # mutate(category = factor(category, levels = levels(MHW_cat$category)))

# The count of the first time the largest category pixel occurs at each pixel and the cumulative values
system.time(
MHW_cat_first <- MHW_cat_pixel_all %>%
  group_by(t, category) %>%
  summarise(first_n = n(), .groups = "drop") %>% 
  right_join(full_grid, by = c("t", "category")) %>%
  mutate(first_n = ifelse(is.na(first_n), 0, first_n),
         first_n_prop = round(first_n/nrow(med_sea_coords), 4)) %>% 
  arrange(t, category) %>% 
  group_by(category) %>%
  mutate(first_n_cum = cumsum(first_n),
         first_area_cum = cumsum(first_area),
         first_n_cum_prop = round(first_n_cum/nrow(med_sea_coords), 4)) %>% 
  ungroup()
) # 1 second

# The count and area of categories of MCSs happening on a given day, and cumulatively throughout the year
# system.time(
MCS_cat_daily <- MCS_cat %>% 
  dplyr::select(lon, lat, t, category, category_correct, category_ice) %>% 
  right_join(lon_lat_OISST_area, by = c("lon", "lat")) %>% 
  pivot_longer(cols = category:category_ice, values_to = "category") %>% 
  group_by(t, name, category) %>%
  summarise(cat_n = n(),
            cat_area = sum(sq_area), .groups = "drop") %>% 
  right_join(full_grid, by = c("t", "category", "name")) %>% 
  mutate(cat_n = ifelse(is.na(cat_n), 0, cat_n),
         cat_n_prop = round(cat_n/nrow(OISST_ocean_coords), 4),
         cat_area = ifelse(is.na(cat_area), 0, cat_area),
         cat_area_prop = round(cat_area/sum(lon_lat_OISST_area$sq_area), 4)) %>% 
  arrange(t, name, category) %>% 
  group_by(name, category) %>%
  mutate(cat_n_cum = cumsum(cat_n),
         cat_area_cum = cumsum(cat_area),
         cat_n_cum_prop = round(cat_n_cum/nrow(OISST_ocean_coords), 4),
         cat_area_cum_prop = round(cat_area_cum/sum(lon_lat_OISST_area$sq_area), 4)) %>% 
  right_join(MCS_cat_first, by = c("t", "name", "category"))
# ) # 1 second
saveRDS(MCS_cat_daily, file = paste0("annual_summary_MCS/MCS_cat_daily_",chosen_year,".Rds"))

