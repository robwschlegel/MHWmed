# code/MHW_summary.R
# This script creates annual and total summaries from the MHW Med results


# Setup -------------------------------------------------------------------

# The needed packages
library(tidyverse)
library(FNN)
library(geosphere)
library(sf)
library(sfheaders)
library(doParallel); registerDoParallel(cores = 7)

# Disable scientific notation
options(scipen = 999)

# The file location
res_files <- dir("data/MHW", full.names = T)

# The coords with SST
load("metadata/med_sea_coords.RData")
lat_idx <- unique(med_sea_coords$lat)
lon_idx <- unique(med_sea_coords$lon)

# MME data
mme <- read_csv("data/Collaborative_tasks_version_database_protected - MME dataset.csv", guess_max = 1000) %>% 
  dplyr::rename(lon = Longitude, lat = Latitude, year = Year) %>% 
  mutate(year = as.numeric(gsub('[.]', '', as.character(year))),
         lon = as.numeric(gsub('[,]', '.', as.character(lon))),
         lat = as.numeric(sub("(.{2})(.*)", "\\1.\\2", lat)),
         `Lower Depth` = as.numeric(gsub('[,]', '.', as.character(`Lower Depth`))),
         `Upper Depth` = as.numeric(gsub('[,]', '.', as.character(`Upper Depth`))),
         `Mortality Lower Depth` = as.numeric(gsub('[,]', '.', as.character(`Mortality Lower Depth`))),
         `Mortality Upper Depth` = as.numeric(gsub('[,]', '.', as.character(`Mortality Upper Depth`)))) %>% 
  dplyr::select(year:`Damaged qualitative`) %>%
  # filter(`Damaged qualitative` != "No", # Filter out 'No' values
         # `Upper Depth` <= 15,
         # Species != "Pinna nobilis") %>% 
  mutate(Ecoregion = case_when(Ecoregion == "Western Mediterranean" & lat >= 39 ~ "Northwestern Mediterranean",
                               Ecoregion == "Western Mediterranean" & lat < 39 ~ "Southwestern Mediterranean",
                               TRUE ~ Ecoregion))
unique(mme$Ecoregion)

# Complete annual dates by categories
full_annual_grid <- expand_grid(year = seq(1982, 2019), 
                                category = as.factor(c("I Moderate", "II Strong", "III Severe", "IV Extreme")))

# Full monthly grid
full_monthly_grid <- expand_grid(year = seq(1982, 2019), 
                                 month = lubridate::month(seq(1:12), label = T, abb = T), 
                                 category = as.factor(c("I Moderate", "II Strong", "III Severe", "IV Extreme")))

# Complete daily dates by categories
full_daily_grid <- expand_grid(t = seq(as.Date(paste0("1982-01-01")), as.Date("2019-12-31"), by = "day"), 
                               category = as.factor(c("I Moderate", "II Strong", "III Severe", "IV Extreme"))) %>% 
  mutate(year = lubridate::year(t))

# The MHW category colour palette
MHW_colours <- c(
  "I Moderate" = "#ffc866",
  "II Strong" = "#ff6900",
  "III Severe" = "#9e0000",
  "IV Extreme" = "#2d0000"
)

# The base global map
map_base <- ggplot2::fortify(maps::map(fill = TRUE, col = "grey80", plot = FALSE)) %>%
  dplyr::rename(lon = long) %>%
  mutate(group = ifelse(lon > 180, group+9999, group),
         lon = ifelse(lon > 180, lon-360, lon))

# The base Med map
med_base <- ggplot() +
  geom_polygon(data = map_base, aes(x = lon, y = lat, group = group), fill = "grey40") +
  theme_void() +
  guides(fill = guide_legend(override.aes = list(size = 10))) #+
  # theme(panel.border = element_rect(colour = "black", fill = NA),
  #       legend.position = "bottom",
  #       legend.text = element_text(size = 14),
  #       legend.title = element_text(size = 16),
  #       panel.background = element_rect(fill = "grey90"))


# Functions ---------------------------------------------------------------

# Function for finding and cleaning up points within a given region polygon
points_in_region <- function(region_in){
  region_sub <- MEOW %>% 
    filter(ECOREGION == region_in) %>% 
    dplyr::select(geometry)
  region_sub <- as.data.frame(region_sub$geometry[[1]][[1]]) %>%
    `colnames<-`(c("lon", "lat"))
    # unnest()
  coords_in <- med_sea_coords %>% 
    mutate(in_grid = sp::point.in.polygon(point.x = med_sea_coords[["lon"]], point.y = med_sea_coords[["lat"]], 
                                          pol.x = region_sub[["lon"]], pol.y = region_sub[["lat"]])) %>% 
    filter(in_grid >= 1) %>% 
    mutate(Ecoregion = region_in) %>% 
    dplyr::select(lon, lat, Ecoregion)
  return(coords_in)
}

# Function for extracting lon/lat values from sf objects
extract_coords <- function(df){
  res_sub <- as.data.frame(df$geometry[[1]][[1]]) %>% 
    `colnames<-`(c("lon", "lat"))
  lat_val <- ifelse(min(res_sub$lat) < 32, min(res_sub$lat), max(res_sub$lat))  
  lon_val <- res_sub$lon[res_sub$lat == lat_val]
  res <- data.frame(lon = lon_val, lat = lat_val)
  return(res)
}

# Function for loading only the daily cat results
load_cat_daily <- function(file_name, lon_range = NA){
  
  # Load data
  res_full <- readRDS(file_name)
  
  # Extract only daily cat values
  cat_clim <- res_full %>% 
    dplyr::select(-event) %>% 
    unnest(cat) %>% 
    filter(row_number() %% 2 == 1) %>% 
    unnest(cat) %>% 
    ungroup() %>% 
    filter(!is.na(t)) %>%
    mutate(category = factor(category),
           year = lubridate::year(t),
           month = lubridate::month(t, label = T, abb = T))
  
  # Trim lon if desired
  if(length(lon_range) > 1){
    cat_clim <- cat_clim %>% 
      filter(lon >= lon_range[1],
             lon <= lon_range[2])
  }
  
  # Clean up and exit
  rm(res_full); gc()
  return(cat_clim)
}

# Load MHW results per pixel
load_event_cat <- function(df){
  
  lat_pixel <- df$lat_sst[1] 
  lon_pixel <- df$lon_sst
  
  # Load data
  res_full <- readRDS(res_files[which(lat_idx == lat_pixel)]) %>% 
    dplyr::filter(lon %in% lon_pixel)
  gc()
  
  # Unnest event and cat results
  event_event <- res_full %>% 
    dplyr::select(-cat) %>% 
    unnest(event) %>% 
    filter(row_number() %% 2 == 0) %>% 
    unnest(event) %>% 
    ungroup()
  cat_event <- res_full %>% 
    dplyr::select(-event) %>% 
    unnest(cat) %>% 
    filter(row_number() %% 2 == 0) %>% 
    unnest(cat) %>% 
    ungroup()
  
  # Merge and exit
  all_event <- left_join(event_event, cat_event, 
                         by = c("lon", "lat", "event_no", "duration", 
                                "intensity_max" = "i_max", "date_peak" = "peak_date"))
  rm(res_full); gc()
  return(all_event)
}

# Find the nearest grid cells for each site
grid_match <- function(coords1, coords2){
  coords2$idx <- 1:nrow(coords2)
  grid_index <- data.frame(coords1,
                           idx = knnx.index(data = as.matrix(coords2[,1:2]),
                                            query = as.matrix(coords1[,1:2]), k = 1))
  grid_points <- left_join(grid_index, coords2, by = c("idx")) %>% 
    mutate(dist = round(distHaversine(cbind(lon.x, lat.x),
                                      cbind(lon.y, lat.y))/1000, 2), idx = NULL)
  return(grid_points)
}

# Function for calculating stats for each individual pixel
cat_pixel_calc <- function(file_name){
  
  # Load data
  MHW_cat <- load_cat_daily(file_name)
  
  # The sum of intensities per pixel for the year
  MHW_intensity <- MHW_cat %>% 
    group_by(lon, lat, year, month) %>% 
    summarise(duration = n(),
              cum_int = sum(intensity), .groups = "drop")
  
  # The count of the highest category of each event in each pixel
  MHW_cat_count <- MHW_cat %>% 
    group_by(lon, lat, year, month, event_no) %>% 
    summarise(category = max(as.integer(category), na.rm = T), .groups = "drop") %>% 
    mutate(category = factor(category, levels = c(1:4),
                             labels = c("I Moderate", "II Strong", "III Severe", "IV Extreme"))) %>%
    dplyr::select(-event_no) %>%
    group_by(lon, lat, year, month, category) %>%
    summarise(count = n(), .groups = "drop") %>% 
    group_by(lon, lat) %>% 
    right_join(full_monthly_grid, by = c("year", "month", "category")) %>% # This used to be a left join...
    pivot_wider(values_from = count, names_from = category) %>% 
    replace(is.na(.), 0)
  
  # The earliest date of the highest category of event
  MHW_cat_pixel <- MHW_cat %>% 
    group_by(lon, lat, year, month) %>%
    filter(as.integer(category) == max(as.integer(category), na.rm = T)) %>% 
    filter(t == min(t)) %>% 
    dplyr::rename(max_int = intensity) %>% 
    unique() %>%
    data.frame() %>% 
    left_join(MHW_intensity, by = c("lon", "lat", "year", "month")) %>% 
    left_join(MHW_cat_count, by = c("lon", "lat", "year", "month")) %>% 
    dplyr::select(lon, lat, year, month, duration, t, event_no, category, everything())
  
  # Clean up and exit
  rm(MHW_cat, MHW_intensity, MHW_cat_count); gc()
  return(MHW_cat_pixel)
}

# Function for calculating stats for each individual year
cat_pixel_annual_calc <- function(sub_months = seq(1, 12)){
  cat_pixel_annual_sum <- MHW_cat_pixel_monthly %>%
    filter(month %in% sub_months) %>% 
    dplyr::select(lon, lat, year, duration, cum_int:`IV Extreme`) %>%
    group_by(lon, lat, year) %>%
    summarise_all(sum)
  gc()
  cat_pixel_annual <- MHW_cat_pixel_monthly %>%
    filter(month %in% sub_months) %>% 
    group_by(lon, lat, year) %>%
    filter(as.integer(category) == max(as.integer(category), na.rm = T)) %>%
    filter(t == min(t)) %>%
    dplyr::select(lon:category, max_int) %>%
    unique() %>%
    left_join(cat_pixel_annual_sum, by = c("lon", "lat", "year")) %>%
    dplyr::rename(duration_max = duration.x,
                  duration_sum = duration.y)
  gc()
  return(cat_pixel_annual)
}

# Function for calculating stats for each individual day
cat_daily_calc <- function(file_name, sub_months = seq(1, 12)){
  
  # Load data
  MHW_cat <- load_cat_daily(file_name)
  
  # Calculate daily MHW occurrence across all pixels
  MHW_cat_daily <- MHW_cat %>% 
    filter(as.numeric(month) %in% sub_months) %>% 
    group_by(year, t, category) %>%
    summarise(cat_n = n(), .groups = "drop") %>% 
    right_join(full_daily_grid, by = c("year", "t", "category")) %>% 
    mutate(cat_n = ifelse(is.na(cat_n), 0, cat_n)) %>% 
    arrange(year, t, category)
  
  # Clean up and exit
  rm(MHW_cat); gc()
  return(MHW_cat_daily)
}

# Function for calculating summary stats
cat_summary_calc <- function(df_pixel, df_daily, JJASON = F){
  
  # The daily count of the first time the largest category pixel occurs over the whole Med and the cumulative values
  cat_first_annual <- df_pixel %>%
    group_by(t, year, category) %>%
    summarise(first_n = n(), .groups = "drop") %>%
    right_join(full_daily_grid, by = c("t", "year", "category")) %>%
    arrange(year, t, category) %>%
    mutate(first_n = ifelse(is.na(first_n), 0, first_n),
           first_n_prop = round(first_n/nrow(med_sea_coords), 4)) %>%
    # arrange(t, category) %>%
    group_by(year, category) %>%
    mutate(first_n_cum = cumsum(first_n),
           first_n_cum_prop = round(first_n_cum/nrow(med_sea_coords), 4)) %>%
    ungroup()
  
  # The count of categories of MHWs happening on a given day, and cumulatively throughout the year
  cat_summary_annual <- df_daily %>%
    arrange(year, t, category) %>%
    group_by(t, year, category) %>%
    summarise(cat_n = sum(cat_n), .groups = "keep") %>%
    mutate(cat_n_prop = round(cat_n/nrow(med_sea_coords), 4)) %>%
    group_by(year, category) %>%
    mutate(cat_n_cum = cumsum(cat_n),
           cat_n_cum_prop = round(cat_n_cum/nrow(med_sea_coords), 4)) %>%
    right_join(cat_first_annual, by = c("t", "year", "category"))
  
  # Filter out days for JJASON
  if(JJASON){
    cat_summary_annual <- filter(cat_summary_annual, 
                                 lubridate::month(t) >= 6,
                                 lubridate::month(t) <= 11)
  }
  
  return(cat_summary_annual)
}

# An in between function to help with RAM use
region_proc <- function(file_name, region_coords){
  
  # Load data
  cat_daily <- load_cat_daily(file_name) %>% 
    right_join(region_coords, by = c("lon", "lat")) %>% 
    filter(!is.na(t))
  
  # Number of pixels each month that experience a MHW
  cat_pixels <- cat_daily %>% 
    dplyr::select(lon, lat, year, month) %>% 
    group_by(year, month) %>% 
    distinct() %>% 
    summarise(pixels = n(), .groups = "drop")
  
  # Spatial coverage
  cat_surface <- cat_daily %>% 
    dplyr::select(lon, lat, year, month) %>% 
    distinct() %>% 
    group_by(year, month) %>% 
    summarise(surface = n()/nrow(region_coords), .groups = "drop") %>% 
    right_join(distinct(full_monthly_grid[,1:2]), by = c("year", "month")) %>% 
    replace(is.na(.), 0)
  
  # Count of category days per year, month, Ecoregion
  # Don't calculate values that account for pixels with NO MHW
  # We only want to know about temperature anomalies from MHWs
  cat_count <- cat_daily %>% 
    group_by(year, month, category) %>% 
    summarise(count = n(), .groups = "drop") %>% 
    right_join(full_monthly_grid, by = c("year", "month", "category")) %>% 
    pivot_wider(values_from = count, names_from = category) %>% 
    replace(is.na(.), 0)
  
  # Calculations for MHW metrics
  cat_calc <- cat_daily %>% 
    group_by(year, month) %>% 
    summarise(duration = n(),
              max_int = max(intensity),
              mean_int = mean(intensity),
              cum_int = sum(intensity), .groups = "drop") %>% 
    pivot_longer(duration:cum_int) %>% 
    right_join(expand.grid(year = seq(1982, 2019), 
                           month = lubridate::month(seq(1:12), label = T, abb = T),
                           name = c("duration", "max_int", "mean_int", "cum_int")), 
               by = c("year", "month", "name")) %>% 
    pivot_wider(values_from = value, names_from = name) %>% 
    arrange(year, month) %>% 
    left_join(cat_pixels, by = c("year", "month")) %>% 
    left_join(cat_count, by = c("year", "month")) %>% 
    left_join(cat_surface, by = c("year", "month")) %>% 
    replace(is.na(.), 0) %>% 
    mutate(region = region_coords$region[1]) %>% 
    dplyr::select(region, year, month, pixels, surface, max_int:cum_int, duration, everything())
  
  # Clean up and exit
  rm(cat_daily, cat_count, cat_surface); gc()
  return(cat_calc)
}

# Regional/seasonal summary calculations
region_calc <- function(region_name){
  
  print(paste0("Began run on ",region_name," at ",Sys.time()))
  
  # Find region coords
  region_coords <- filter(med_regions, region == region_name) %>% 
    distinct()
  file_sub <- data.frame(lat_index = seq_len(length(unique(med_sea_coords$lat))),
                         lat = unique(med_sea_coords$lat)) %>% 
    filter(lat %in% region_coords$lat)
  
  # load necessary files
  # registerDoParallel(cores = 7)
  cat_res <- plyr::ldply(res_files[file_sub$lat_index], region_proc,
                         .parallel = F,
                         region_coords = region_coords) %>% 
    group_by(region, year, month) %>%
    summarise(pixels = sum(pixels),
              surface = sum(surface),
              duration = sum(duration),
              max_int = mean(max_int),
              mean_int = mean(mean_int),
              cum_int = sum(cum_int),
              `I Moderate` = sum(`I Moderate`),
              `II Strong` = sum(`II Strong`),
              `III Severe` = sum(`III Severe`),
              `IV Extreme` = sum(`IV Extreme`), .groups = "drop")
  # Clean and exit
  rm(region_coords, file_sub); gc() # Free up some RAM
  return(cat_res)
}

# Event analysis of specific years
event_analysis <- function(){}

# Function for creating multi-faceted stacked barplots for ecoregion summaries
ecoregion_trend_fig <- function(region_sub, 
                                year_range = seq(1982, 2019), 
                                month_range = lubridate::month(seq(6, 11), label = T, abb = T)){
  
  # File name
  fig_name <- paste0("figures/MHW_trend_",region_sub,"_",
                     min(year_range),"-",max(year_range),"_",
                     min(month_range),"-",max(month_range),".png")
  if(region_sub == "Tunisian Plateau/Gulf of Sidra"){
    fig_name <- str_remove(fig_name, "/Gulf of Sidra")
  }
  
  # Plot
  region_plot <- MHW_cat_region %>% 
    filter(region == region_sub,
           year %in% year_range,
           month %in% month_range) %>% 
    pivot_longer(cols = surface:`IV Extreme`) %>%
    mutate(value = case_when(!name %in%  c("surface", "max_int", "mean_int") & pixels > 0 ~ value/pixels,
                             TRUE ~ value),
           name = factor(name,
                         levels = c("surface", "mean_int", "max_int", "cum_int", "duration", 
                                    "I Moderate", "II Strong", "III Severe", "IV Extreme"))) %>%
    ggplot(aes(x = year, y = value, colour = month)) +
    geom_point() +
    # geom_line(alpha = 0.5) +
    geom_smooth(method = "lm", se = F, linetype = "dashed") +
    facet_wrap(~name, scales = "free_y") +
    labs(y = NULL, title = region_sub)
  
  # Save
  ggsave(fig_name, region_plot, height = 10, width = 10)
}

# Warm season figure per ecoregion
ecoregion_summary_fig <- function(region_sub, 
                                  year_range = seq(2015, 2019), 
                                  month_range = lubridate::month(seq(6, 11), label = T, abb = T)){
  
  # File name
  fig_name <- paste0("figures/MHW_summary_",region_sub,"_",
                     min(year_range),"-",max(year_range),"_",
                     min(month_range),"-",max(month_range),".png")
  if(region_sub == "Tunisian Plateau/Gulf of Sidra"){
    fig_name <- str_remove(fig_name, "/Gulf of Sidra")
  }
  
  # Plot
  region_plot <- MHW_cat_region %>% 
    filter(region == region_sub,
           year %in% year_range,
           month %in% month_range) %>% 
    pivot_longer(cols = surface:`IV Extreme`) %>% 
    mutate(value = case_when(!name %in%  c("surface", "max_int", "mean_int") & pixels > 0 ~ value/pixels,
                             TRUE ~ value),
           name = factor(name,
                         levels = c("surface", "mean_int", "max_int", "cum_int", "duration", 
                                    "I Moderate", "II Strong", "III Severe", "IV Extreme"))) %>%
    ggplot(aes(x = year, y = value)) +
    geom_bar(aes(fill = month), 
             stat = "identity", 
             show.legend = T,
             position = "dodge",
             # position = position_stack(reverse = TRUE), 
             width = 1) +
    facet_wrap(~name, scales = "free_y") +
    labs(y = NULL, title = region_sub)
  # region_plot
  
  # Save
  ggsave(fig_name, region_plot, height = 10, width = 10)
}

# Annual summary figure
annual_summary_fig <- function(chosen_year){
  
  paste0("Started run on ",chosen_year," at ", Sys.time())
  
  # Chose the year of categories to display
  MHW_cat_filter <- filter(MHW_cat_summary_annual, year == chosen_year)
  MHW_cat_pixel_filter <- filter(MHW_cat_pixel_annual, year == chosen_year)
  gc()
  
  # Extract small data.frame for easier labeling
  MHW_cat_filter_labels <- MHW_cat_filter %>% 
    group_by(category) %>% 
    filter(t == max(t)) %>% 
    ungroup() %>% 
    mutate(label_first_n_cum = cumsum(first_n_cum_prop))
  
  # Title
  fig_title <- paste0("Mediterranean MHW categories of ",chosen_year,
                      "\nNOAA Med SST ~4km; Climatogy period: 1982-2011")
  
  ## Create figures
  # Global map of MHW occurrence
  fig_map <- ggplot(MHW_cat_pixel_filter, aes(x = lon, y = lat)) +
    # geom_tile(data = OISST_ice_coords, fill = "powderblue", colour = NA, alpha = 0.5) +
    geom_tile(aes(fill = category), colour = NA) +
    geom_polygon(data = map_base, aes(x = lon, y = lat, group = group)) +
    scale_fill_manual("Category", values = MHW_colours) +
    coord_cartesian(expand = F, 
                    xlim = c(min(med_sea_coords$lon), max(med_sea_coords$lon)),
                    ylim = c(min(med_sea_coords$lat), max(med_sea_coords$lat))) +
    theme_void() +
    guides(fill = guide_legend(override.aes = list(size = 10))) +
    theme(panel.border = element_rect(colour = "black", fill = NA),
          legend.position = "bottom",
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          panel.background = element_rect(fill = "grey90"))
  # fig_map
  
  # Stacked barplot of global daily count of MHWs by category
  fig_count <- ggplot(MHW_cat_filter, aes(x = t, y = cat_n_prop)) +
    geom_bar(aes(fill = category), stat = "identity", show.legend = F,
             position = position_stack(reverse = TRUE), width = 1) +
    scale_fill_manual("Category", values = MHW_colours) +
    scale_y_continuous(limits = c(0, 1),
                       breaks = seq(0.2, 0.8, length.out = 4),
                       labels = paste0(seq(20, 80, by = 20), "%")) +
    scale_x_date(date_breaks = "2 months", date_labels = "%Y-%m") +
    labs(y = "Daily MHWs\n(non-cumulative)", x = "Day of the year") +
    coord_cartesian(expand = F) +
    theme(panel.border = element_rect(colour = "black", fill = NA),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 13))
  # fig_count
  
  # Stacked barplot of cumulative percent of Mediterranean affected by MHWs
  fig_cum <- ggplot(MHW_cat_filter, aes(x = t, y = first_n_cum_prop)) +
    geom_bar(aes(fill = category), stat = "identity", show.legend = F,
             position = position_stack(reverse = TRUE), width = 1) +
    geom_hline(data = MHW_cat_filter_labels, show.legend = F,
               aes(yintercept = label_first_n_cum, colour = category)) +
    scale_fill_manual("Category", values = MHW_colours) +
    scale_colour_manual("Category", values = MHW_colours) +
    scale_y_continuous(limits = c(0, 1),
                       breaks = seq(0.2, 0.8, length.out = 4),
                       labels = paste0(seq(20, 80, by = 20), "%")) +
    scale_x_date(date_breaks = "2 months", date_labels = "%Y-%m") +
    labs(y = "Top MHW category per pixel\n(cumulative)", x = "Day of first occurrence") +
    coord_cartesian(expand = F) +
    theme(panel.border = element_rect(colour = "black", fill = NA),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 13))
  # fig_cum
  
  # Stacked barplot of average cumulative MHW days per pixel
  fig_prop <- ggplot(MHW_cat_filter, aes(x = t, y = cat_n_cum_prop)) +
    geom_bar(aes(fill = category), stat = "identity", show.legend = F,
             position = position_stack(reverse = TRUE), width = 1) +
    scale_fill_manual("Category", values = MHW_colours) +
    scale_y_continuous(breaks = round(seq(sum(MHW_cat_filter_labels$cat_n_cum_prop)*0.25,
                                          sum(MHW_cat_filter_labels$cat_n_cum_prop)*0.75, length.out = 3), 0)) +
    scale_x_date(date_breaks = "2 months", date_labels = "%Y-%m") +  
    labs(y = "Average MHW days per pixel\n(cumulative)", x = "Day of the year") +
    coord_cartesian(expand = F) +
    theme(panel.border = element_rect(colour = "black", fill = NA),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 13))
  # fig_prop
  
  # print("Combining figures")
  fig_ALL_sub <- ggpubr::ggarrange(fig_count, fig_cum, fig_prop, ncol = 3, align = "hv",
                                   labels = c("B)", "C)", "D)"), font.label = list(size = 16))
  fig_ALL <- ggpubr::ggarrange(fig_map, fig_ALL_sub, ncol = 1, heights = c(1, 0.6),
                               labels = c("A)"), common.legend = T, legend = "bottom",
                               font.label = list(size = 16))
  
  # Standard caption technique
  fig_ALL_cap <- grid::textGrob(fig_title, x = 0.01, just = "left", gp = grid::gpar(fontsize = 20))
  fig_ALL_cap <- ggpubr::ggarrange(fig_ALL_cap, fig_ALL, heights = c(0.07, 1), nrow = 2)
  
  # print("Saving final figure")
  ggsave(fig_ALL_cap, height = 12, width = 18, 
         filename = paste0("figures/MHW_cat_summary_",chosen_year,".png"))
  # RAM help
  rm(MHW_cat_filter, MHW_cat_pixel_filter);gc()
}

# total summary figure
total_summary_fig <- function(df){
  
  # Detect JJASON period
  if(max(lubridate::month(df$t)) == 12){
    JJASON_bit <- ""
    month_sub <- 12
    day_sub <- 31
    dd <- 1
    y2_labs <- c("5%", "10%", "15%", "20%")
  }  else {
    JJASON_bit <- " (JJASON)"
    month_sub <- 11
    day_sub <- 30
    dd <- 2
    y2_labs <- c("10%", "20%", "30%", "40%")
  }
  
  # Load OISST annual global MHW summaries
  OISST_global <- readRDS("data/OISST_cat_daily_1992-2018_total.Rds") %>% 
    group_by(t) %>% 
    mutate(cat_n_prop_stack = cumsum(cat_n_prop),
           first_n_cum_prop_stack = cumsum(first_n_cum_prop))
  
  # Total summary
  # Create mean values of daily count
  cat_daily_mean <- df %>%
    group_by(year, category) %>%
    summarise(cat_n_prop_mean = mean(cat_n_prop, na.rm = T),
              cat_n_cum_prop = max(cat_n_cum_prop, na.rm = T), .groups = "drop")
  
  # Extract only values from December 31st
  cat_daily <- df %>%
    group_by(year, category) %>%
    filter(lubridate::month(t) == month_sub, lubridate::day(t) == day_sub)
  
  # Stacked barplot of global daily count of MHWs by category
  fig_count_historic <- ggplot(cat_daily_mean, aes(x = year, y = cat_n_cum_prop)) +
    geom_bar(aes(fill = category), stat = "identity", show.legend = T,
             position = position_stack(reverse = TRUE), width = 1) +
    geom_line(data = OISST_global, aes(x = t, y = cat_n_prop_stack/dd, colour = category), 
              linetype = "dotted", show.legend = F) +
    geom_point(data = OISST_global, aes(x = t, y = cat_n_prop_stack/dd, fill = category), 
               shape = 21, show.legend = F) +
    scale_fill_manual("Category", values = MHW_colours) +
    scale_colour_manual("Category", values = MHW_colours) +
    scale_y_continuous(limits = c(0, 80),
                       breaks = seq(20, 60, length.out = 3),
                       sec.axis = sec_axis(name = "Average daily MHW coverage", trans = ~ . + 0,
                                           breaks = c(18.25, 36.5, 54.75, 73),
                                           labels = y2_labs)) +
    scale_x_continuous(breaks = seq(1984, 2019, 7)) +
    guides(pattern_colour = FALSE, colour = FALSE) +
    labs(y = "Average MHW days", x = NULL) +
    coord_cartesian(expand = F) +
    theme(panel.border = element_rect(colour = "black", fill = NA),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12))
  # fig_count_historic
  
  # Stacked barplot of cumulative percent of ocean affected by MHWs
  fig_cum_historic <- ggplot(cat_daily, aes(x = year, y = first_n_cum_prop)) +
    geom_bar(aes(fill = category), stat = "identity", show.legend = T,
             position = position_stack(reverse = TRUE), width = 1) +
    geom_line(data = OISST_global, aes(x = t, y = first_n_cum_prop_stack, colour = category), 
              linetype = "dotted", show.legend = F) +
    geom_point(data = OISST_global, aes(x = t, y = first_n_cum_prop_stack, fill = category), 
               shape = 21, show.legend = F) +
    scale_fill_manual("Category", values = MHW_colours) +
    scale_colour_manual("Category", values = MHW_colours) +
    scale_y_continuous(position = "right", 
                       limits = c(0, 1),
                       breaks = seq(0.2, 0.8, length.out = 4),
                       labels = paste0(seq(20, 80, by = 20), "%")) +
    scale_x_continuous(breaks = seq(1984, 2019, 7)) +
    labs(y = "Total MHW coverage", x = NULL) +
    coord_cartesian(expand = F) +
    theme(panel.border = element_rect(colour = "black", fill = NA),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.position = "none",
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12))
  # fig_cum_historic
  
  # Create the figure title
  min_year <- min(cat_daily_mean$year)
  max_year <- max(cat_daily_mean$year)
  fig_title <- paste0("Mediterranean MHW categories summary: ",min_year," - ", max_year, JJASON_bit,
                      "\nCMEMS Med SST ~4km; Climatogy period: 1982-2011")
  
  # Stick them together and save
  fig_ALL_historic <- ggpubr::ggarrange(fig_count_historic, fig_cum_historic,
                                        ncol = 2, align = "hv", labels = c("A)", "B)"), hjust = -0.1,
                                        font.label = list(size = 14), common.legend = T, legend = "bottom")
  fig_ALL_cap <- grid::textGrob(fig_title, x = 0.01, just = "left", gp = grid::gpar(fontsize = 18))
  fig_ALL_full <- ggpubr::ggarrange(fig_ALL_cap, fig_ALL_historic, heights = c(0.25, 1), nrow = 2)
  # ggsave(fig_ALL_full, filename = "figures/MHW_cat_historic.png", height = 4.25, width = 12)
}

# Function to create monthly maps of warm season 2015-2019
monthly_map_fig_one <- function(month_choice, year_choice, common_scales){
  
  # Filter data
  monthly_data <- MHW_cat_region %>% 
    filter(year == year_choice, month == month_choice) %>% 
    left_join(MEOW, by = c("region" = "ECOREGION"))
  
  ## Plot summaries per variable
  # Surface
  plot_surface <- med_base + 
    geom_sf(data = monthly_data, alpha = 0.9, aes(geometry = geometry, fill = surface)) +
    coord_sf(expand = F, xlim = c(-10, 45), ylim = c(25, 50)) +
    scale_fill_viridis_c("Prop. of surface", option = "E", 
                         limits = c(0, max(common_scales$surface)))
  # Duration
  plot_duration <- med_base + 
    geom_sf(data = monthly_data, alpha = 0.9, aes(geometry = geometry, fill = duration/pixels)) +
    coord_sf(expand = F, xlim = c(-10, 45), ylim = c(25, 50)) +
    scale_fill_viridis_c("Duration of MHWs (days)", option = "D", 
                         limits = c(0, max(common_scales$duration)))
  # Mean intensity
  plot_mean_int <- med_base + 
    geom_sf(data = monthly_data, alpha = 0.9, aes(geometry = geometry, fill = mean_int)) +
    coord_sf(expand = F, xlim = c(-10, 45), ylim = c(25, 50)) +
    scale_fill_viridis_c("Mean intensity (°C)", option = "A", 
                         limits = c(0, max(common_scales$mean_int)))
  # Max intensity
  plot_max_int <- med_base + 
    geom_sf(data = monthly_data, alpha = 0.9, aes(geometry = geometry, fill = max_int)) +
    coord_sf(expand = F, xlim = c(-10, 45), ylim = c(25, 50)) +
    scale_fill_viridis_c("Max intensity (°C)", option = "B", 
                         limits = c(0, max(common_scales$max_int)))
  # Cumulative intensity
  plot_cum_int <- med_base + 
    geom_sf(data = monthly_data, alpha = 0.9, aes(geometry = geometry, fill = cum_int/pixels)) +
    coord_sf(expand = F, xlim = c(-10, 45), ylim = c(25, 50)) +
    scale_fill_viridis_c("Cum. intensity (°C days)", option = "C", 
                         limits = c(0, max(common_scales$cum_int)))
  
  # Title
  fig_title <- paste0(year_choice,": ",month_choice)
  # Combine into one tall figure
  fig_ALL <- ggpubr::ggarrange(plot_surface, plot_duration, plot_mean_int, 
                               plot_max_int, plot_cum_int, ncol = 1)
  fig_ALL_cap <- grid::textGrob(fig_title, x = 0.01, just = "left", gp = grid::gpar(fontsize = 30))
  fig_ALL_cap <- ggpubr::ggarrange(fig_ALL_cap, fig_ALL, heights = c(0.05, 1), nrow = 2)
  return(fig_ALL_cap)
}

# Figure to stick one year of plots together
monthly_map_fig_full <- function(year_choice){
  
  # Create common scales
  common_scales <- MHW_cat_region %>% 
    filter(year >= 2015) %>% 
    mutate(cum_int = cum_int/pixels, 
           duration = duration/pixels) %>% 
    dplyr::select(surface, max_int, mean_int, cum_int, duration) %>% 
    na.omit() %>% 
    distinct()
  
  # Create the six months of figures
  fig_jun <- monthly_map_fig_one("Jun", year_choice, common_scales)
  fig_jul <- monthly_map_fig_one("Jul", year_choice, common_scales)
  fig_aug <- monthly_map_fig_one("Aug", year_choice, common_scales)
  fig_sep <- monthly_map_fig_one("Sep", year_choice, common_scales)
  fig_oct <- monthly_map_fig_one("Oct", year_choice, common_scales)
  fig_nov <- monthly_map_fig_one("Nov", year_choice, common_scales)
  
  # Save and exit
  fig_full <- ggpubr::ggarrange(fig_jun, fig_jul, fig_aug, fig_sep, fig_oct, fig_nov, nrow = 1)
  ggsave(paste0("figures/MHW_monthly_ecoregions_",year_choice,".png"), fig_full, height = 18, width = 42)
}

# Figure that plots the per pixel maps
monthly_map_pixel <- function(var_choice,
                              annual = T,
                              year_range = seq(2015, 2019), 
                              month_range = lubridate::month(seq(6, 11), label = F, abb = T)){
  
  # Reduce the dataframe to the desired  dimensions
  if(annual){
    MHW_cat_pixel_filter <- MHW_cat_pixel_annual_JJASON %>% 
      filter(year %in% year_range)
  } else{
    MHW_cat_pixel_filter <- MHW_cat_pixel_monthly %>% 
      filter(year %in% year_range,
             month %in% month_range)
  }
  gc()
  
  # Ecoregions for faceting
  monthly_MEOW <- MHW_cat_region %>% 
    mutate(month = as.integer(month)) %>% 
    filter(year %in% year_range, month %in% month_range) %>% 
    left_join(MEOW, by = c("region" = "ECOREGION"))
  
  # Prepare MME points
  mme_points <- mme %>% 
    filter(year %in% year_range,
           EvenStart %in% c("Summer", "Autumn"))
  
  # Complete region/year grid
  full_region_year_grid <- expand_grid(year = year_range, 
                                       Ecoregion = unique(mme$Ecoregion))
  # Prepare MME labels
  mme_labels <- mme_points %>% 
    group_by(year, Ecoregion) %>% 
    summarise(count = n(), .groups = "drop") %>% 
    right_join(full_region_year_grid, by = c("year", "Ecoregion")) %>% 
    left_join(MEOW_label_coords, by = c("Ecoregion" = "ECOREGION")) %>% 
    mutate(count = ifelse(is.na(count), 0, count)) %>% 
    arrange(year, Ecoregion)
  
  # Determine colour for MME dots
  if(var_choice != "category"){
    point_colour <- "forestgreen"
  } else{
    point_colour <- "black"
  }
  
  # Global map of MHW occurrence
  multi_panel <- med_base + 
    geom_tile(data = MHW_cat_pixel_filter, colour = NA,
              aes_string(fill = var_choice, x = "lon", y = "lat")) +
    geom_sf(data = monthly_MEOW, alpha = 1, aes(geometry = geometry), fill = NA, colour = "grey70") +
    geom_point(data = mme_points, aes(x = lon, y = lat), shape = 1, colour = point_colour) +
    geom_label(data = mme_labels, aes(x = lon, y = lat, label = count)) +
    coord_sf(expand = F, xlim = c(-10, 45), ylim = c(25, 50)) +
    theme(panel.border = element_rect(colour = "black", fill = NA),
          # legend.position = "bottom",
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 18),
          panel.background = element_rect(fill = "grey90"), 
          strip.text = element_text(size = 16))
  if(annual){
    multi_panel <- multi_panel + 
      facet_wrap(~year, ncol = 3) +
      theme(legend.direction = "horizontal", 
            legend.position = c(0.82, 0.3))
  } else{
    multi_panel <- multi_panel + 
      facet_grid(year ~ month) +
      theme(legend.position = "bottom")
  } 
  if(var_choice == "duration") multi_panel <- multi_panel + 
    scale_fill_viridis_c("Duration")
  if(var_choice == "cum_int") multi_panel <- multi_panel + 
    scale_fill_viridis_c("Cumulative\nIntensity", option = "B")
  if(var_choice == "category") multi_panel <- multi_panel + 
    scale_fill_manual("Category", values = MHW_colours) + 
    guides(fill = guide_legend(nrow = 2, byrow = TRUE))
  rm(MHW_cat_pixel_filter); gc()
  return(multi_panel)
}
#

# Ecoregions --------------------------------------------------------------

# Load MEOW
MEOW <- read_sf("metadata/MEOW/meow_ecos.shp") %>% 
  filter(PROVINCE == "Mediterranean Sea")

## Create Northwestern + Southwestern Mediterranean regions
# Extract lon/lat values
MEOW_sub <- MEOW %>% 
  filter(ECOREGION == "Western Mediterranean") %>% 
  dplyr::select(geometry)
MEOW_sub <- as.data.frame(MEOW_sub$geometry[[1]][[1]]) %>%
  `colnames<-`(c("lon", "lat"))
# Create polygons for each new region
NW_polygon <- MEOW_sub %>%
  filter(lat >= 39.1) %>% 
  sf_multipolygon()
st_crs(NW_polygon) <- 4326
SW_polygon <- MEOW_sub %>%
  filter(lat <= 39.3) %>% 
  sf_multipolygon()
st_crs(SW_polygon) <- 4326
# Create new data.frame
MEOW_new <- MEOW[c(7,7),] %>% 
  mutate(ECOREGION = c("Northwestern Mediterranean", 
                       "Southwestern Mediterranean"),
         geometry = c(NW_polygon$geometry, SW_polygon$geometry))
# Reintroduce to MEOW
MEOW <- rbind(MEOW[-7,], MEOW_new)
unique(MEOW$ECOREGION)

# Find SST pixels within Med MEOW
registerDoParallel(cores = 7)
med_regions <- plyr::ldply(unique(MEOW$ECOREGION), points_in_region, .parallel = T)

# Prepare MME labels 32
MEOW_label_coords <- MEOW %>% 
  group_by(ECOREGION) %>% 
  nest() %>% 
  mutate(coords = purrr::map(data, extract_coords)) %>% 
  dplyr::select(-data) %>% 
  unnest(coords)


# Ecoregion summaries -----------------------------------------------------

# NB: Do not run in parallel
# system.time(
# MHW_cat_region <- plyr::ldply(unique(med_regions$region), region_calc, .parallel = F)
# ) # ~150 seconds for 1 on 7 cores, ~18 minutes total
# save(MHW_cat_region, file = "data/MHW_cat_region.RData")
# readr::write_csv(MHW_cat_region, "data/MHW_cat_region.csv")
load("data/MHW_cat_region.RData")


# Ecoregion summary figures -----------------------------------------------

# plyr::l_ply(unique(med_regions$region), ecoregion_summary_fig, .parallel = T)


# Ecoregion trend figures -------------------------------------------------

# plyr::l_ply(unique(med_regions$region), ecoregion_trend_fig, .parallel = T)


# Map summary figures -----------------------------------------------------

# plyr::l_ply(2015:2019, monthly_map_fig_full, .parallel = T)


# Annual summaries --------------------------------------------------------

# The occurrences per month per pixel
# doParallel::registerDoParallel(cores = 15)
# system.time(
# MHW_cat_pixel_monthly <- plyr::ldply(res_files, cat_pixel_calc, .parallel = T)
# ) # ~15 minutes on 7 cores, 344 seconds on 15 cores
# save(MHW_cat_pixel_monthly, file = "data/MHW_cat_pixel_monthly.RData")
load("data/MHW_cat_pixel_monthly.RData") # This is very large, only load if necessary

# The occurrences per year per pixel
## NB: Requires MHW_cat_pixel_monthly
# MHW_cat_pixel_annual <- cat_pixel_annual_calc()
# save(MHW_cat_pixel_annual, file = "data/MHW_cat_pixel_annual.RData")
load("data/MHW_cat_pixel_annual.RData")

# The occurrences per year per pixel JJASON
## NB: Requires MHW_cat_pixel_monthly
# system.time(
# MHW_cat_pixel_annual_JJASON <- cat_pixel_annual_calc(sub_months = seq(6, 11))
# ) # 300 seconds
# save(MHW_cat_pixel_annual_JJASON, file = "data/MHW_cat_pixel_annual_JJASON.RData")
load("data/MHW_cat_pixel_annual_JJASON.RData")

# The occurrences per day all months
# registerDoParallel(cores = 15)
# system.time(
# MHW_cat_daily_annual <- plyr::ldply(res_files, cat_daily_calc, .parallel = T)
# ) # 227 seconds on 15 cores
# save(MHW_cat_daily_annual, file = "data/MHW_cat_daily_annual.RData")
load("data/MHW_cat_daily_annual.RData")

# The occurrences per day JJASON
# system.time(
# MHW_cat_daily_annual_JJASON <- plyr::ldply(res_files, cat_daily_calc, .parallel = T, sub_months = seq(6, 11))
# ) # 230 seconds on 15 cores
# save(MHW_cat_daily_annual_JJASON, file = "data/MHW_cat_daily_annual_JJASON.RData")
load("data/MHW_cat_daily_annual_JJASON.RData")


# Total summaries ---------------------------------------------------------

# The daily count of the first time the largest category pixel occurs over the whole Med and the cumulative values
# MHW_cat_summary_annual <- cat_summary_calc(MHW_cat_pixel_annual, MHW_cat_daily_annual)
# save(MHW_cat_summary_annual, file = "data/MHW_cat_summary_annual.RData")
load("data/MHW_cat_summary_annual.RData")

# Same as above but for JJASON
# MHW_cat_summary_annual_JJASON <- cat_summary_calc(MHW_cat_pixel_annual_JJASON, MHW_cat_daily_annual_JJASON, JJASON = T)
# save(MHW_cat_summary_annual_JJASON, file = "data/MHW_cat_summary_annual_JJASON.RData")
load("data/MHW_cat_summary_annual_JJASON.RData")


# Summary figures ---------------------------------------------------------

## NB: These require objects to be in the environment that are added by the above code
# MHW_cat_summary_annual, MHW_cat_pixel_annual

# Create annual summary figures
# NB: This is very RAM heavy
# doParallel::registerDoParallel(cores = 2)
# plyr::l_ply(1982:2019, annual_summary_fig, .parallel = T)

# Create total summary figure
total_summary <- total_summary_fig(MHW_cat_summary_annual)
ggsave("figures/MHW_cat_historic.png", total_summary, height = 4.25, width = 8)
total_summary_JJASON <- total_summary_fig(MHW_cat_summary_annual_JJASON)
ggsave("figures/MHW_cat_historic_JJASON.png", total_summary_JJASON, height = 4.25, width = 8)


# Per pixel maps with MME -------------------------------------------------

# Requires: MHW_cat_pixel_monthly.RData and MHW_cat_region.RData

# Per pixel maps
map_pixel_duration <- monthly_map_pixel("duration", annual = T)
ggsave("figures/MHW_pixel_duration.png", map_pixel_duration, height = 7, width = 20)
map_pixel_category <- monthly_map_pixel("category", annual = T)
ggsave("figures/MHW_pixel_category.png", map_pixel_category, height = 7, width = 20)
map_pixel_cum_int <- monthly_map_pixel("cum_int", annual = T)
ggsave("figures/MHW_pixel_cum_int.png", map_pixel_cum_int, height = 7, width = 20)


# Single summary map ------------------------------------------------------

# Requires: MHW_cat_pixel_monthly, MHW_cat_pixel_annual_JJASON
load("data/MHW_cat_pixel_annual_JJASON.RData")

# Prepare MME points
mme_damage <- mme %>% 
  filter(`Damaged qualitative` != "No",
         `Upper Depth` <= 15,
         Species != "Pinna nobilis",
         EvenStart %in% c("Summer", "Autumn")) %>% 
  group_by(lon, lat) %>% 
  summarise(count = n(), .groups = "drop")
mme_no <- mme %>% 
  filter(`Damaged qualitative` == "No",
         `Upper Depth` <= 15,
         Species != "Pinna nobilis",
         EvenStart %in% c("Summer", "Autumn")) %>% 
  group_by(lon, lat) %>% 
  summarise(count = n(), .groups = "drop")

# Determine historic medians per pixel
MHW_cat_pixel_JJASON_clim_median <- MHW_cat_pixel_annual_JJASON %>% 
  filter(year %in% seq(1982, 2014)) %>%
  group_by(lon, lat) %>% 
  summarise(duration_median = median(duration, na.rm = T),
            icum_median = median(cum_int, na.rm = T), .groups = "drop")

# Create study period anomalies
MHW_cat_pixel_JJASON_anom <- MHW_cat_pixel_annual_JJASON %>% 
  right_join(med_regions, by = c("lon", "lat")) %>% 
  filter(year %in% seq(2015, 2019)) %>%
  group_by(lon, lat) %>% 
  summarise(duration = mean(duration, na.rm = T),
            icum = mean(cum_int, na.rm = T), .groups = "drop") %>% 
  left_join(MHW_cat_pixel_JJASON_clim_median, by = c("lon", "lat")) %>% 
  mutate(duration_anom = duration-duration_median,
         icum_anom = icum-icum_median) %>% 
  mutate(duration_anom = ifelse(duration_anom < 0, 0, duration_anom),
         icum_anom = ifelse(icum_anom < 0, 0, icum_anom))
  
# Map of duration anom for JJASON
anom_plot_dur <- med_base + 
  geom_tile(data = MHW_cat_pixel_JJASON_anom, colour = NA,
            aes(fill = duration_anom, x = lon, y = lat)) +
  geom_contour(data = MHW_cat_pixel_JJASON_anom, colour = "black",
               aes(x = lon, y = lat, z = duration_anom), breaks = c(1)) +
  geom_sf(data = MEOW, alpha = 1, aes(geometry = geometry), fill = NA, colour = "grey70") +
  scale_fill_gradient2(low = "white", high = "forestgreen") +
  coord_sf(expand = F, xlim = c(-10, 45), ylim = c(25, 50)) +
  labs(fill = "Duration (days)",
       title = "MHW Duration",
       subtitle = "Average per year for 2015-2019 JJASON in excess of the 1982-2011 JJASON median") +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        legend.position = "bottom",
        title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        panel.background = element_rect(fill = "grey90"), 
        strip.text = element_text(size = 16))
anom_plot_dur

# Map of icum anom for JJASON
anom_plot_icum <- med_base + 
  geom_tile(data = MHW_cat_pixel_JJASON_anom, colour = NA,
            aes(fill = icum_anom, x = lon, y = lat)) +
  geom_contour(data = MHW_cat_pixel_JJASON_anom, colour = "black",
               aes(x = lon, y = lat, z = icum_anom), breaks = c(1)) +
  geom_sf(data = MEOW, alpha = 1, aes(geometry = geometry), fill = NA, colour = "grey70") +
  scale_fill_gradient2(low = "white", high = "darkorchid") +
  coord_sf(expand = F, xlim = c(-10, 45), ylim = c(25, 50)) +
  labs(fill = "Cumulative\nintensity (°C days)", 
       title = "MHW Cumulative Intensity",
       subtitle = "Average per year for 2015-2019 JJASON in excess of the 1982-2011 JJASON median") +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        legend.position = "bottom",
        title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        panel.background = element_rect(fill = "grey90"), 
        strip.text = element_text(size = 16))
anom_plot_icum

# MME damage
anom_plot_mme <- med_base + 
  geom_sf(data = MEOW, alpha = 1, aes(geometry = geometry), fill = NA, colour = "grey70") +
  geom_point(data = mme_damage, aes(x = lon, y = lat, size = count), 
             alpha = 0.7, shape = 21, fill = "yellow", colour = "red") +
  coord_sf(expand = F, xlim = c(-10, 45), ylim = c(25, 50)) +
  labs(size = "MME count (n)",
       title = "MME damage",
       subtitle = "Recorded damaging events for 2015-2019 Summer/Autumn") +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        legend.position = "bottom",
        title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        panel.background = element_rect(fill = "grey90"), 
        strip.text = element_text(size = 16))
# anom_plot_mme

# MME damage
anom_plot_mme_no <- med_base + 
  geom_sf(data = MEOW, alpha = 1, aes(geometry = geometry), fill = NA, colour = "grey70") +
  geom_point(data = mme_no, aes(x = lon, y = lat, size = count), 
             alpha = 0.7, shape = 21, fill = "yellow", colour = "red") +
  coord_sf(expand = F, xlim = c(-10, 45), ylim = c(25, 50)) +
  labs(size = "MME count (n)",
       title = "MME no damage",
       subtitle = "Recorded non-damaging events for 2015-2019 Summer/Autumn") +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        legend.position = "bottom",
        title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        panel.background = element_rect(fill = "grey90"), 
        strip.text = element_text(size = 16))
# anom_plot_mme_no

# Combine and save
anom_all <- ggpubr::ggarrange(anom_plot_dur, anom_plot_icum, anom_plot_mme, anom_plot_mme_no,
                              ncol = 2, nrow = 2, align = "hv")
ggsave("figures/MHW_pixel_median_anom.png", anom_all, height = 16, width = 22)


# MME vs MHW pixels -------------------------------------------------------

# Scale the figures so that it is the prop for the pixels in the region
# Also look at scatterplot by separating out the low, mid, high impact MME
# Also do a comparison of all depths vs top 15

# Create another scatterplot showing when no MME are recorded vs when they are
# Use size to show count, and colour to show quality of impact
# Get a curved line in there that shows how the yes and no values can be separated

# Match pixels
lon_lat_match <- grid_match(distinct(dplyr::select(mme, lon, lat)), distinct(dplyr::select(med_regions, lon, lat))) %>% 
  dplyr::rename(lon_mme = lon.x, lat_mme = lat.x, lon_sst = lon.y, lat_sst = lat.y) %>% 
  left_join(med_regions, by = c("lon_sst" = "lon", "lat_sst" = "lat"))

# Create full annual grid
lon_lat_match_full_grid <- expand_grid(year = seq(1982, 2019), lon_lat_match)

# Extract MHW results for paired pixels
doParallel::registerDoParallel(cores = 15)
site_event_cat <- plyr::ddply(lon_lat_match, c("lat_sst"), load_event_cat, .parallel = T)

# Calculate annual MHW stats per pixel
site_MHW_summary <- site_event_cat %>% 
  dplyr::rename(lon_sst = lon) %>% 
  dplyr::select(-lat) %>% 
  mutate(year = lubridate::year(date_peak),
         month = lubridate::month(date_peak)) %>%
  filter(month %in% seq(6, 11)) %>% 
  group_by(lon_sst, lat_sst, year) %>% 
  summarise(count_MHW = n(),
            duration = sum(duration, na.rm = T),
            # imean = round(mean(intensity_mean, na.rm = T), 2),
            icum = round(sum(intensity_cumulative, na.rm = T)), .groups = "drop")

# Calculate annual MME stats per site
site_MME_summary <- mme %>% 
  filter(`Damaged qualitative` != "No",
         `Upper Depth` <= 15,
         Species != "Pinna nobilis",
         EvenStart %in% c("Summer", "Autumn")) %>% 
  group_by(year, Ecoregion, Location, lon, lat) %>% 
  summarise(`Damaged percentage` = round(mean(`Damaged percentage`, na.rm = T)),
            count_MME = n(), .groups = "drop") %>% 
  dplyr::rename(lon_mme = lon, lat_mme = lat)

# Join them together
site_MME_MHW_summary <- lon_lat_match_full_grid %>% 
  left_join(site_MHW_summary, by = c("lon_sst", "lat_sst", "year")) %>% 
  replace(is.na(.), 0) %>%  # Fill in the no MHW rows with 0's
  left_join(site_MME_summary, by = c("lon_mme", "lat_mme", "Ecoregion", "year")) %>% 
  filter(year >= 2015)
write_csv(site_MME_MHW_summary, "data/site_MME_MHW_summary.csv")

# Scatterplot of MME and MHW count summaries
scatter_MME_MHW_site <- site_MME_MHW_summary %>% 
  na.omit() %>% 
  pivot_longer(count_MHW:icum) %>% 
  ggplot(aes(x = value, y = count_MME)) +
  geom_point(aes(colour = Ecoregion)) +
  geom_smooth(aes(colour = Ecoregion), method = "lm", se = F) +
  labs(y = "MME count per site/year", x = NULL,
       title = "MME count per site and year compared to MHW metrics (JJASON)") +
  facet_wrap(~name, scales = "free_x", strip.position = "bottom") +
  theme(legend.position = "bottom", 
        strip.placement = "outside", strip.background = element_blank())
scatter_MME_MHW_site

# Scatterplot of MME and MHW damage summaries
scatter_MME_MHW_site_dam <- site_MME_MHW_summary %>% 
  na.omit() %>% 
  pivot_longer(count_MHW:icum) %>% 
  ggplot(aes(x = value, y = `Damaged percentage`)) +
  geom_point(aes(colour = Ecoregion)) +
  geom_smooth(aes(colour = Ecoregion), method = "lm", se = F) +
  labs(y = "MME damage per site/year", x = NULL,
       title = "MME damage per site and year compared to MHW metrics (JJASON)") +
  facet_wrap(~name, scales = "free_x", strip.position = "bottom") +
  theme(legend.position = "bottom", 
        strip.placement = "outside", strip.background = element_blank())
scatter_MME_MHW_site_dam

# Scatterplot of MME and MHW count summaries by ecoregion
scatter_MME_MHW_ecoregion <- site_MME_MHW_summary %>% 
  na.omit() %>% 
  group_by(Ecoregion, year) %>% 
  summarise(`Damaged percentage` = mean(`Damaged percentage`),
            count_MME_sum = sum(count_MME),
            count_MME_mean = mean(count_MME),
            count_MHW_sum = sum(count_MHW),
            # count_MHW_mean = mean(count_MHW),
            duration = mean(duration),
            # imean = mean(imean),
            icum = mean(icum), .groups = "drop") %>% 
  pivot_longer(count_MHW_sum:icum) %>% 
  ggplot(aes(x = value, y = count_MME_sum)) +
  geom_point(aes(colour = Ecoregion)) +
  geom_smooth(aes(colour = Ecoregion), method = "lm", se = F) +
  labs(y = "MME count per ecoregion/year", x = NULL,
       title = "MME count per ecoregion and year compared to MHW metrics (JJASON)") +
  facet_wrap(~name, scales = "free_x", strip.position = "bottom") +
  theme(legend.position = "bottom", 
        strip.placement = "outside", strip.background = element_blank())
scatter_MME_MHW_ecoregion

# Scatterplot of MME and MHW damage summaries by ecoregion
scatter_MME_MHW_ecoregion_dam <- site_MME_MHW_summary %>% 
  na.omit() %>% 
  group_by(Ecoregion, year) %>% 
  summarise(`Damaged percentage` = mean(`Damaged percentage`),
            count_MME_sum = sum(count_MME),
            count_MME_mean = mean(count_MME),
            count_MHW_sum = sum(count_MHW),
            # count_MHW_mean = mean(count_MHW),
            duration = mean(duration),
            # imean = mean(imean),
            icum = mean(icum), .groups = "drop") %>% 
  pivot_longer(count_MHW_sum:icum) %>% 
  ggplot(aes(x = value, y = `Damaged percentage`)) +
  geom_point(aes(colour = Ecoregion)) +
  geom_smooth(aes(colour = Ecoregion), method = "lm", se = F) +
  labs(y = "MME damage (%) per ecoregion/year", x = NULL,
       title = "MME damage per ecoregion and year compared to MHW metrics (JJASON)") +
  facet_wrap(~name, scales = "free_x", strip.position = "bottom") +
  theme(legend.position = "bottom", 
        strip.placement = "outside", strip.background = element_blank())
scatter_MME_MHW_ecoregion_dam

scatter_MME_MHW <- ggpubr::ggarrange(scatter_MME_MHW_site, scatter_MME_MHW_site_dam,
                                     scatter_MME_MHW_ecoregion, scatter_MME_MHW_ecoregion_dam, 
                                     ncol = 2, nrow = 2, align = "hv", common.legend = T, legend = "bottom")
ggsave("figures/scatter_MME_MHW.png", height = 6, width = 14)


# Spatial MME vs MHW maps -------------------------------------------------

# MHW days may be the easiest to communicate
# °C is important, but may not be as accessible
# May need to have multiple maps. Faceted in one figure.

# MME filtering
# Definitely run comparisons with MHW at the smallest scale
# Pixel level at the nearest sites with MME records
# MMEs on the surface should be 15 or shallower
# The regularly monitored MME sites are going to be labeled in the main spreadsheet over the week
# For comparisons try the multiple different temporal levels: all 5 years, per year, per season per year, whole period per ecoregion
# Barplots with top pointing bar for MME and downward pointing bar for MHW

# Will need to extract a coastal stripe of pixels

# Comparisons of MME and MHW per pixel vs. whole ecoregion

# The multi-panel MHW metric maps will be better in a supplemental, so we will want to keep them in circulation
# Or maybe we do want multi-panels because MMEs change so much from year to year spatially

# A good investigation would be to see where MMEs occurred when there were no or very low MHW
# This would be done as part of the per pixel comparisons
# A second layer would be to see where MMEs do/don't occur compared against 
# the median MHW anomaly values against the climatology period

# Comparisons must be made by species
# Gorgonians are a good starting point as they are almost certainly temperature driven MME
# Certain species are better choices than others: Paramuricea clavata, Eunicella singularis, Eunicella cavolini, Corallium rubrum
# The first three can be grouped as necessary
# Gorgonians tend to start dying more often near the end of summer
# Fish are a bad choice, and Sapia

# To the maps also add a panel showing SSTmax. Could show trends, too.

# 4 maps to bring to the overall group
# Trend in SSTmax or 99th percentile JJASON shown as the total (slope * total years)
# Trend for average SST JJASON as above  
# MHW duration, median anomaly perhaps JJASON
# Cumulative intensity as for Duration JJASON

# We as the experts should chose the MHW metric to show
# Show how this metric relates to MME by species for spatial, temporal, depth ranges

# A figure somehow showing areas that were monitored but did not have mortality would be interesting
# A boxplot of some sort

# Need different figures for different genus etc.

# When did the highest MME occur in a year, and what did the MHW look like then

# Keep it simple. Figures that are a summary of all five years.
# The goal is to be able to show all of the info in one or two figures.

# Could use alpha to show count/days of MHWs over five years and colour for mean icum
# Like a density plot, sort of...
# Could use alpha with histograms to show overlay of different years for MME and MHW stats

# Time series barplots per site will be good to show as back up
# But instead of showing MHW time series, show the occurrence/days/icum as bars per season with MME

# Consider allowing for a minimum limit of MMEs in a region etc. before including it in the stats

# Show the MME rug plot bits by colour for different taxa


# Temporal summary figure -------------------------------------------------



# MHW metric time series and MME rug plot ---------------------------------



# Histograms of MME and MHW -----------------------------------------------

# Requires: site_MME_MHW_summary
site_MME_MHW_summary <- read_csv("data/site_MME_MHW_summary.csv")

# Prepare MME points
mme_points <- mme %>% 
  filter(Species != "Pinna nobilis",
         `Damaged qualitative` != "No",
         # `Upper Depth` <= 15,
         EvenStart %in% c("Summer", "Autumn"))

# Complete region/year grid
full_region_year_grid <- expand_grid(year = seq(2015, 2019), 
                                     Ecoregion = unique(mme$Ecoregion))
# Prepare MME labels
mme_labels <- mme_points %>% 
  group_by(year, Ecoregion) %>% 
  summarise(count = n(), .groups = "drop") %>% 
  right_join(full_region_year_grid, by = c("year", "Ecoregion")) %>% 
  mutate(count = ifelse(is.na(count), 0, count)) %>% 
  arrange(year, Ecoregion)

# Filter data to target time period
ecoregion_MME_MHW <- site_MME_MHW_summary %>% 
  dplyr::select(Ecoregion, year, count_MHW:count_MME, -Location) %>% 
  group_by(Ecoregion, year) %>% 
  summarise(`Damaged percentage` = mean(`Damaged percentage`, na.rm = T),
            count_MME_sum = sum(count_MME, na.rm = T),
            count_MME_mean = mean(count_MME, na.rm = T),
            count_MHW_sum = sum(count_MHW, na.rm = T),
            # count_MHW_mean = mean(count_MHW),
            duration = mean(duration, na.rm = T),
            # imean = mean(imean),
            icum = mean(icum, na.rm = T), .groups = "drop") %>% 
  left_join(mme_labels, by = c("year", "Ecoregion")) %>%
  replace(is.na(.), 0)

# Create a whole Med summary
ecoregion_MME_MHW_med <- ecoregion_MME_MHW %>% 
  group_by(year) %>% 
  summarise_all(mean, na.rm = T, .groups = "drop") %>% 
  mutate(Ecoregion = "Mediterranean",
         count_MME_sum = round(count_MME_sum))

# Combine and order factor for plotting
ecoregion_MME_MHW_all <- rbind(ecoregion_MME_MHW, ecoregion_MME_MHW_med) %>% 
  mutate(Ecoregion = factor(Ecoregion, 
                            levels = c("Mediterranean",
                                       "Alboran Sea", "Northwestern Mediterranean", 
                                       "Southwestern Mediterranean", "Adriatic Sea",
                                       "Ionian Sea", "Tunisian Plateau/Gulf of Sidra",
                                       "Aegean Sea", "Levantine Sea")),
         count = round(count))

# Barplot of durations
bar_dur <- ecoregion_MME_MHW_all %>% 
  ggplot(aes(x = year, y = duration)) +
  geom_bar(aes(fill = icum), 
           colour = "black",
           stat = "identity", 
           show.legend = T,
           position = "dodge",
           # position = position_stack(reverse = TRUE), 
           width = 1) +
  geom_label(aes(label = count)) +
  scale_fill_viridis_c("Cumulative\nIntensity (°C days)", option = "B") +
  facet_wrap(~Ecoregion) +
  scale_y_continuous(limits = c(0, 105), breaks = c(25, 50, 75, 100)) +
  scale_x_continuous(breaks = c(2015, 2017, 2019)) +
  coord_cartesian(expand = F) +
  labs(x = NULL, y = "Duration (days)", 
       title = "Total MHW days over JJASON period",
       subtitle = "Bar colour shows cumulative intensity and labels show MME count") +
  # guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5)) +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        # legend.position = c(0.83, 0.16),
        legend.position = "bottom",
        # legend.direction = "horizontal",
        legend.key.width = unit(1, "cm"),
        panel.background = element_rect(fill = "grey90"), 
        strip.text = element_text(size = 12))
bar_dur
ggsave("figures/MHW_ecoregion_summary.png", bar_dur, height = 8, width = 8)

