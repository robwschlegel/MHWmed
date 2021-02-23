# talks/brief_talk.R
# The code used for talks/brief_talk.Rmd


# Setup -------------------------------------------------------------------

library(tidyverse)
library(tidync)


# NetCDF info -------------------------------------------------------------

# Latitude range
tidync("~/pCloudDrive/MHWmed_data/SST/1982_SST_MED_SST_L4_REP_OBSERVATIONS_010_021_a_1449092207440.nc") %>% activate("D1")

# Longitude range
tidync("~/pCloudDrive/MHWmed_data/SST/1982_SST_MED_SST_L4_REP_OBSERVATIONS_010_021_a_1449092207440.nc") %>% activate("D2")


# Match small spreadsheet to local MHW ------------------------------------

# Load local spreadsheet
# NB: This is not stored on GitHub as it is not public data

# Extract SST pixels based on spreadsheet lon/lat range

# Create mean SST time series and detect events

# Figure showing MHWs + MMEs


# NW Med stats ------------------------------------------------------------

# Load results
load("data/MHW_cat_region.RData")

# NW Med data
NW_Med <- filter(MHW_cat_region, region == "Northwestern Mediterranean")

# The warm periods only
NW_Med_warm <- filter(NW_Med, month %in% c("Jun", "Jul", "Aug", "Sep", "Oct", "Nov"))

# Sum of spatial cover per year
NW_med_warm_sum <- NW_Med_warm %>% 
  group_by(year) %>% 
  summarise(surface = sum(surface),
            duration = sum(duration/pixels),
            cum_int = sum(cum_int/pixels))

# Single values for plotting
stats_old <- NW_med_warm_sum %>% 
  filter(year >= 1982, year <= 2014) %>% 
  summarise(surface = mean(surface),
            duration = mean(duration, na.rm = T),
            cum_int = mean(cum_int, na.rm = T))
stats_new <- NW_med_warm_sum %>% 
  filter(year >= 2015, year <= 2019) %>% 
  summarise(surface = mean(surface),
            duration = mean(duration, na.rm = T),
            cum_int = mean(cum_int, na.rm = T))


# NW Med figures ----------------------------------------------------------

# Bar plots of surface cover
bar_surface <- ggplot(data = NW_Med_warm, aes(x = year, y = surface)) +
  # Bars
  geom_bar(aes(fill = month), 
           stat = "identity", 
           show.legend = T,
           # position = "dodge",
           position = position_stack(reverse = TRUE),
           width = 1) +
  # Old period average
  geom_segment(aes(x = 1981.5, xend = 2014.5, 
                   y = stats_old$surface, yend = stats_old$surface),
               colour = "darkviolet", size = 1.5) +
  geom_label(aes(x = 1998, y = stats_old$surface, 
                 label = round(stats_old$surface, 1)), 
             label.size = 2, colour = "darkviolet") +
  # New period average
  geom_segment(aes(x = 2014.5, xend = 2019.5, 
                   y = stats_new$surface, yend = stats_new$surface),
               colour = "hotpink", size = 1.5) +
  geom_label(aes(x = 2017, y = stats_new$surface, 
                 label = round(stats_new$surface, 1)), 
             label.size = 2, colour = "hotpink") +
  # Other
  scale_y_continuous(breaks = c(1, 2, 3, 4, 5)) +
  scale_x_continuous(breaks = c(1985, 1995, 2005, 2015)) +
  coord_cartesian(ylim = c(0, 6), expand = F) +
  labs(x = NULL, y = "Sum of region covered (proportion)", fill = "Month",
       title = "Northwestern Mediterranean: surface area of MHWs per month", 
       subtitle = "The stacked bars will add up to 6 if all six months experienced full MHW coverage")
ggsave("talks/graph/NWMed_bar_surface.png", bar_surface, width = 7)

# Bar plots for duration
bar_duration <- ggplot(data = NW_Med_warm, aes(x = year, y = duration/pixels)) +
  # Bars
  geom_bar(aes(fill = month), 
           stat = "identity", 
           show.legend = T,
           # position = "dodge",
           position = position_stack(reverse = TRUE),
           width = 1) +
  # Old period average
  geom_segment(aes(x = 1981.5, xend = 2014.5,
                   y = stats_old$duration, yend = stats_old$duration),
               colour = "darkviolet", size = 1.5) +
  geom_label(aes(x = 1998, y = stats_old$duration,
                 label = round(stats_old$duration, 1)),
             label.size = 2, colour = "darkviolet") +
  # New period average
  geom_segment(aes(x = 2014.5, xend = 2019.5,
                   y = stats_new$duration, yend = stats_new$duration),
               colour = "hotpink", size = 1.5) +
  geom_label(aes(x = 2017, y = stats_new$duration,
                 label = round(stats_new$duration, 1)),
             label.size = 2, colour = "hotpink") +
  # Other
  # scale_y_continuous(breaks = c(1, 2, 3, 4, 5)) +
  scale_x_continuous(breaks = c(1985, 1995, 2005, 2015)) +
  coord_cartesian(expand = F) +
  labs(x = NULL, y = "Sum of duration (days)", fill = "Month",
       title = "Northwestern Mediterranean: duration of MHWs per month", 
       subtitle = "The stacked bars show the sum of MHW days experienced per year")
ggsave("talks/graph/NWMed_bar_duration.png", bar_duration, width = 7)

# Bar plots for duration
bar_cum_int <- ggplot(data = NW_Med_warm, aes(x = year, y = cum_int/pixels)) +
  # Bars
  geom_bar(aes(fill = month), 
           stat = "identity", 
           show.legend = T,
           # position = "dodge",
           position = position_stack(reverse = TRUE),
           width = 1) +
  # Old period average
  geom_segment(aes(x = 1981.5, xend = 2014.5,
                   y = stats_old$cum_int, yend = stats_old$cum_int),
               colour = "darkviolet", size = 1.5) +
  geom_label(aes(x = 1998, y = stats_old$cum_int,
                 label = round(stats_old$cum_int, 1)),
             label.size = 2, colour = "darkviolet") +
  # New period average
  geom_segment(aes(x = 2014.5, xend = 2019.5,
                   y = stats_new$cum_int, yend = stats_new$cum_int),
               colour = "hotpink", size = 1.5) +
  geom_label(aes(x = 2017, y = stats_new$cum_int,
                 label = round(stats_new$cum_int, 1)),
             label.size = 2, colour = "hotpink") +
  # Other
  # scale_y_continuous(breaks = c(1, 2, 3, 4, 5)) +
  scale_x_continuous(breaks = c(1985, 1995, 2005, 2015)) +
  coord_cartesian(expand = F) +
  labs(x = NULL, y = "Sum of cum. intensity (°C days)", fill = "Month",
       title = "Northwestern Mediterranean: cum. intensity of MHWs per month", 
       subtitle = "The stacked bars show the sum of cum. intensities experienced per year")
ggsave("talks/graph/NWMed_bar_cum_int.png", bar_cum_int, width = 7)

