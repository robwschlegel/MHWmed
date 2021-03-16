# MHW_other.R
# This script contains analyses performed in addition to the main project


# Setup -------------------------------------------------------------------

library(tidyverse)
library(heatwaveR)
library(tidync)
library(FNN)
library(doParallel); registerDoParallel(cores = 7)

# Load data from Israel coast
  # NB: These data are not on GitHub as they are not public
  # Contact Gil Rilov for access
coast_raw <- read_csv("~/R/forOthers/Gil/Core intertidal sites water logger data.csv", guess_max = 100000) %>% 
  filter(!is.na(Date)) %>% 
  mutate(Date = as.Date(Date, format = "%d/%m/%Y")) %>% 
  # Fix dates
  separate(Date, c("year", "month", "day")) %>% 
  mutate(year = as.numeric(year), 
         month = as.numeric(month),
         day = as.numeric(day)) %>% 
  mutate(year = case_when(year < 2000 ~ 2020,
                          TRUE ~ year)) %>% 
  unite("t", year:day, sep = "-") %>% 
  mutate(t = as.Date(t))

# Create daily means
coast_daily <- coast_raw %>% 
  select(-Time) %>% 
  group_by(t) %>% 
  summarise_if(is.numeric, mean, na.rm = T) %>% 
  pivot_longer(Achziv:Palmachim, names_to = "site", values_to = "temp") %>% 
  mutate(temp = ifelse(is.na(temp), NA, temp)) %>% 
  arrange(site, t) %>% 
  mutate(type = "in situ")

# Site coordinates
coast_coords <- data.frame(site = c("Achziv", "Shikmona", "Habonim", "Palmachim"),
                           lat = c(33.06375, 32.82467, 32.63018, 31.93024),
                           lon = c(35.10381, 34.95417, 34.9197, 34.6981))

# Complete table for forcing joins
full_annual_grid <- expand_grid(site = c("Achziv", "Shikmona", "Habonim", "Palmachim"),
                                year = seq(2009, 2020),
                                type = c("in situ", "SST"))

# The coords with SST
load("metadata/med_sea_coords.RData")
med_sea_coords$env_index <- 1:nrow(med_sea_coords)

# The data locations
# NB: These files are not hosted on GitHub as they are too large
# Contact Robert Schlegel to receive them: robert.schlegel@imev-mer.fr
med_SST_files <- dir("~/pCloudDrive/MHWmed_data/SST", pattern = ".nc", full.names = T, recursive = T)


# Find nearest SST pixels -------------------------------------------------

coast_coords_pixel <- coast_coords %>% 
  mutate(env_index = as.vector(knnx.index(as.matrix(med_sea_coords[,c("lon", "lat")]),
                                          as.matrix(.[,2:3]), k = 1))) %>% 
  left_join(med_sea_coords, by = "env_index") %>% 
  dplyr::rename(lon = lon.x, lat = lat.x, lon_pixel = lon.y, lat_pixel = lat.y) %>% 
  dplyr::select(-env_index)


# Load SST data -----------------------------------------------------------

# Function to load a latitude subset of a single NetCDF file
load_nc_pixel <- function(file_name, lon_pixel, lat_pixel){
  SST_sub <- tidync(file_name) %>% 
    hyper_filter(lat = lat == lat_pixel,
                 lon = lon == lon_pixel) %>%
    hyper_tibble() %>% 
    mutate(t = as.Date(as.POSIXct(time, origin = "1981-01-01")),
           temp = round(analysed_sst - 273.15, 2)) %>% 
    dplyr::select(lon, lat, t, temp)
}

# Convenience wrapper
load_nc_wrap <- function(site_sub){
  sub_coords <- filter(coast_coords_pixel, site == site_sub)
  res <- plyr::ldply(med_SST_files, load_nc_pixel, .parallel = T, 
                     lon_pixel = sub_coords$lon_pixel, 
                     lat_pixel = sub_coords$lat_pixel)
}

# Get the data
# pixel_dat <- plyr::ldply(coast_coords_pixel$site, load_nc_wrap, .parallel = F) %>% 
#   left_join(coast_coords_pixel, by = c("lon" = "lon_pixel", "lat" = "lat_pixel"))
# saveRDS(pixel_dat, "data/extract/levant_pixels.Rds")
# write_csv(dplyr::select(pixel_dat, site, t, temp), "data/extract/levant_pixels.csv")
pixel_dat <- readRDS("data/extract/levant_pixels.Rds")
pixel_dat <- pixel_dat %>% 
  mutate(type = "SST") %>% 
  dplyr::select(site, type, t, temp)


# Detect MHWs -------------------------------------------------------------

# All results
all_MHW <- rbind(coast_daily, pixel_dat) %>% 
  filter(t >= "2011-01-01",
         t <= "2019-12-31") %>% 
  group_by(site, type) %>%
  nest() %>% 
  mutate(clim = purrr::map(data, ts2clm, climatologyPeriod = c("2011-01-01", "2019-12-31")),
         event = purrr::map(clim, detect_event, categories = T, climatology = T, season = "peak", S = F))%>% 
  select(-data, -clim)

# Clims only
all_clim <- all_MHW %>% 
  unnest(event) %>% 
  filter(row_number() %% 2 == 1) %>% 
  unnest(event) %>% 
  ungroup()

# Events only
all_event <- all_MHW %>% 
  unnest(event) %>% 
  filter(row_number() %% 2 == 0) %>% 
  unnest(event) %>% 
  ungroup()


# Summarise results -------------------------------------------------------

summary_clim <- all_clim %>% 
  filter(event_no > 0) %>%
  mutate(year = lubridate::year(t)) %>% 
  group_by(site, type, year) %>% 
  mutate(mean_int = mean(intensity, na.rm = T),
         max_int = max(intensity, na.rm = T),
         cum_int = sum(intensity, na.rm = T),
         duration = n()) %>% 
  group_by(site, type, year, category) %>% 
  mutate(dur_cat = n()) %>% 
  ungroup() %>% 
  dplyr::select(site, type, year, category, mean_int, max_int, cum_int, duration, dur_cat) %>% 
  pivot_wider(names_from = "category", values_from = "dur_cat", values_fn = mean) %>% 
  distinct() %>% 
  dplyr::select(-`NA`) %>% 
  right_join(full_annual_grid, by = c("site", "type", "year")) %>%
  replace(is.na(.), 0) %>% 
  arrange(site, type, year)


# Compare MHW results -----------------------------------------------------

# Line and dot plot
summary_clim %>% 
  pivot_longer(mean_int:`IV Extreme`, names_to = "name", values_to = "value") %>%
  mutate(name = factor(name,
                       levels = c("mean_int", "max_int", "cum_int", "duration", 
                                  "I Moderate", "II Strong", "III Severe", "IV Extreme"))) %>% 
  ggplot(aes(x = year, y = value, colour = site)) +
  geom_point(aes(shape = type)) +
  # geom_line(alpha = 0.5) +
  geom_smooth(method = "lm", se = F, aes(linetype = type)) +
  facet_wrap(~name, scales = "free_y") +
  labs(y = NULL) +
  theme(legend.position = "bottom")

# Boxplot
summary_clim %>% 
  pivot_longer(mean_int:`IV Extreme`, names_to = "name", values_to = "value") %>%
  mutate(name = factor(name,
                       levels = c("mean_int", "max_int", "cum_int", "duration", 
                                  "I Moderate", "II Strong", "III Severe", "IV Extreme"))) %>% 
  ggplot(aes(x = site, y = value, fill = type)) +
  geom_boxplot(notch = T) +
  # geom_point(aes(shape = type)) +
  # geom_line(alpha = 0.5) +
  # geom_smooth(method = "lm", se = F, linetype = "dashed") +
  facet_wrap(~name, scales = "free_y") +
  labs(y = NULL) +
  theme(legend.position = "bottom")
