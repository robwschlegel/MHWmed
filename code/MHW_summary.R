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
                               category = as.factor(c("I Moderate", "II Strong", "III Severe", "IV Extreme"))) %>% 
  mutate(year = lubridate::year(t))

# The MHW category colour palette
MHW_colours <- c(
  "I Moderate" = "#ffc866",
  "II Strong" = "#ff6900",
  "III Severe" = "#9e0000",
  "IV Extreme" = "#2d0000"
)

# The base map
map_base <- ggplot2::fortify(maps::map(fill = TRUE, col = "grey80", plot = FALSE)) %>%
  dplyr::rename(lon = long) %>%
  mutate(group = ifelse(lon > 180, group+9999, group),
         lon = ifelse(lon > 180, lon-360, lon))


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

# Function for calculating stats for each individual pixel
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

# Function for calculating stats for each individual day
cat_daily_calc <- function(file_name){
  
  # Load data
  MHW_cat <- load_cat_daily(file_name)
  
  # Calculate daily MHW occurrence across all pixels
  MHW_cat_daily <- MHW_cat %>% 
    group_by(year, t, category) %>%
    summarise(cat_n = n(), .groups = "drop") %>% 
    right_join(full_daily_grid, by = c("year", "t", "category")) %>% 
    mutate(cat_n = ifelse(is.na(cat_n), 0, cat_n)) %>% 
    arrange(year, t, category)
  
  # Clean up and exit
  rm(MHW_cat); gc()
  return(MHW_cat_daily)
}

# Annual summary figure
annual_summary_fig <- function(chosen_year){
  
  # Chose the year of categories to display
  MHW_cat_filter <- filter(MHW_cat_summary_annual, year == chosen_year)
  MHW_cat_pixel_filter <- filter(MHW_cat_pixel_annual, year == chosen_year)
  
  # Extract small data.frame for easier labelling
  MHW_cat_filter_labels <- MHW_cat_filter %>% 
    group_by(category) %>% 
    filter(t == max(t)) %>% 
    ungroup() %>% 
    mutate(label_first_n_cum = cumsum(first_n_cum_prop))
  
  # Title
  fig_title <- paste0("Mediterranean MHWs: ",chosen_year)
  
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
    labs(y = "Top MHW category\n(cumulative)", x = "Day of first occurrence") +
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
    labs(y = "Average MHW days\n(cumulative)", x = "Day of the year") +
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
}

# total summary figure
total_summary_fig <- function(){
  
  # Total summary
  # Create mean values of daily count
  cat_daily_mean <- MHW_cat_summary_annual %>%
    group_by(year, category) %>%
    summarise(cat_n_prop_mean = mean(cat_n_prop, na.rm = T),
              cat_n_cum_prop = max(cat_n_cum_prop, na.rm = T), .groups = "drop")
  
  # Extract only values from December 31st
  cat_daily <- MHW_cat_first_annual %>%
    filter(lubridate::month(t) == 12, lubridate::day(t) == 31) %>%
    left_join(cat_daily_mean, by = c("year", "category"))
  
  # Stacked barplot of global daily count of MHWs by category
  fig_count_historic <- ggplot(cat_daily, aes(x = year, y = cat_n_prop_mean)) +
    geom_bar(aes(fill = category), stat = "identity", show.legend = T,
             position = position_stack(reverse = TRUE), width = 1) +
    scale_fill_manual("Category", values = MHW_colours) +
    scale_y_continuous(limits = c(0, 1),
                       breaks = seq(0.2, 0.8, length.out = 4),
                       labels = paste0(seq(20, 80, by = 20), "%")) +
    scale_x_continuous(breaks = seq(1982, 2019, 5)) +
    labs(y = "Daily MHWs\n(annual average)", x = NULL) +
    coord_cartesian(expand = F) +
    theme(panel.border = element_rect(colour = "black", fill = NA),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 16))
  # fig_count_historic
  
  # Stacked barplot of cumulative percent of ocean affected by MHWs
  fig_cum_historic <- ggplot(cat_daily, aes(x = year, y = first_n_cum_prop)) +
    geom_bar(aes(fill = category), stat = "identity", show.legend = T,
             position = position_stack(reverse = TRUE), width = 1) +
    scale_fill_manual("Category", values = MHW_colours) +
    scale_y_continuous(limits = c(0, 1),
                       breaks = seq(0.2, 0.8, length.out = 4),
                       labels = paste0(seq(20, 80, by = 20), "%")) +
    scale_x_continuous(breaks = seq(1982, 2019, 5)) +
    labs(y = "Top MHW category\n(annual total)", x = NULL) +
    coord_cartesian(expand = F) +
    theme(panel.border = element_rect(colour = "black", fill = NA),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 16))
  # fig_cum_historic
  
  # Stacked barplot of average cumulative MHW days per pixel
  fig_prop_historic <- ggplot(cat_daily, aes(x = year, y = cat_n_cum_prop)) +
    geom_bar(aes(fill = category), stat = "identity", show.legend = T,
             position = position_stack(reverse = TRUE), width = 1) +
    scale_fill_manual("Category", values = MHW_colours) +
    scale_y_continuous(limits = c(0, 80),
                       breaks = seq(20, 60, length.out = 3)) +
    scale_x_continuous(breaks = seq(1982, 2019, 5)) +
    labs(y = "Average MHW days\n(annual max)", x = NULL) +
    coord_cartesian(expand = F) +
    theme(panel.border = element_rect(colour = "black", fill = NA),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 16))
  # fig_prop_historic
  
  # Create the figure title
  fig_title <- paste0("Mediterranean MHW category summaries: 1982 - 2019")
  
  # Stick them together and save
  fig_ALL_historic <- ggpubr::ggarrange(fig_count_historic, fig_cum_historic, fig_prop_historic,
                                        ncol = 3, align = "hv", labels = c("(A)", "(B)", "(C)"), hjust = -0.1,
                                        font.label = list(size = 14), common.legend = T, legend = "bottom")
  ggsave(fig_ALL_historic, filename = "figures/MHW_cat_historic.png", height = 4.25, width = 12)
}


# Summaries ---------------------------------------------------------------

# The occurrences per pixel
# system.time(
# MHW_cat_pixel_annual <- plyr::ldply(res_files, cat_pixel_calc, .parallel = T)
# ) # 435 seconds on 7 cores
# save(MHW_cat_pixel_annual, file = "data/MHW_cat_pixel_annual.RData")
load("data/MHW_cat_pixel_annual.RData")

# The occurrences per day
# system.time(
# MHW_cat_daily_annual <- plyr::ldply(res_files, cat_daily_calc, .parallel = T)
# ) # 258 seconds on 7 cores
# save(MHW_cat_daily_annual, file = "data/MHW_cat_daily_annual.RData")
load("data/MHW_cat_daily_annual.RData")

### Summaries of total data

# The daily count of the first time the largest category pixel occurs over the whole Med and the cumulative values
MHW_cat_first_annual <- MHW_cat_pixel_annual %>%
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
MHW_cat_summary_annual <- MHW_cat_daily_annual %>% 
  arrange(year, t, category) %>% 
  group_by(t, year, category) %>% 
  summarise(cat_n = sum(cat_n), .groups = "keep") %>% 
  mutate(cat_n_prop = round(cat_n/nrow(med_sea_coords), 4)) %>% 
  group_by(year, category) %>%
  mutate(cat_n_cum = cumsum(cat_n),
         cat_n_cum_prop = round(cat_n_cum/nrow(med_sea_coords), 4)) %>% 
  right_join(MHW_cat_first_annual, by = c("t", "year", "category"))
save(MHW_cat_summary_annual, file = "data/MHW_cat_summary_annual.RData")


# Figures -----------------------------------------------------------------
# NB: These require objects to be in the environment that are added by the above code
# They are not standalone functions...

# Create annual summary figures
plyr::l_ply(1982:2019, annual_summary_fig, .parallel = T)

# Create total summary figure
total_summary_fig()

