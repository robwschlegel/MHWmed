# talks/brief_talk.R
# The code used for talks/brief_talk.Rmd


# Setup -------------------------------------------------------------------

library(tidyverse)
library(tidync)
library(FNN)
library(geosphere)
library(heatwaveR)
library(doParallel); registerDoParallel(cores = 6)

# The data locations
# NB: These files are not hosted on GitHub as they are too large
# Contact Robert Schlegel to receive them: robert.schlegel@imev-mer.fr
med_SST_files <- dir("~/pCloudDrive/MHWmed_data/SST", pattern = ".nc", full.names = T, recursive = T)

# The lon/lat indexes
med_lat <- tidync(med_SST_files[1]) %>% 
  activate("D1") %>% 
  hyper_tibble()
med_lon <- tidync(med_SST_files[1]) %>% 
  activate("D2") %>% 
  hyper_tibble()

# Function to load a latitude subset of a single NetCDF file
# testers...
# file_name <- med_SST_files[100]
# lon_range <- range(lon_range_pixels$lon); lat_range <- range(lat_range_pixels$lat)
load_nc_range <- function(file_name, lon_range, lat_range){
  SST_sub <- tidync(file_name) %>% 
    hyper_filter(lon = between(lon, lon_range[1], lon_range[2]),
                 lat = between(lat, lat_range[1], lat_range[2])) %>% 
    # hyper_filter(lat = lat == med_lat$lat[lat_row],
                 # lon = lon == med_lon$lon[1]) %>%
    hyper_tibble() %>% 
    mutate(t = as.Date(as.POSIXct(time, origin = "1981-01-01")),
           temp = round(analysed_sst - 273.15, 2)) %>% 
    dplyr::select(lon, lat, t, temp)
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

# Set line colours
lineColCat <- c(
  "Temperature" = "black",
  "Climatology" = "grey20",
  "Threshold" = "darkgreen",
  "2x Threshold" = "darkgreen",
  "3x Threshold" = "darkgreen",
  "4x Threshold" = "darkgreen"
)

# Set category fill colours
fillColCat <- c(
  "Moderate" = "#ffc866",
  "Strong" = "#ff6900",
  "Severe" = "#9e0000",
  "Extreme" = "#2d0000"
)
fillColCatNum <- c(
  "I Moderate" = "#ffc866",
  "II Strong" = "#ff6900",
  "III Severe" = "#9e0000",
  "IV Extreme" = "#2d0000"
)


# NetCDF info -------------------------------------------------------------

# Latitude range
tidync("~/pCloudDrive/MHWmed_data/SST/1982_SST_MED_SST_L4_REP_OBSERVATIONS_010_021_a_1449092207440.nc") %>% activate("D1")

# Longitude range
tidync("~/pCloudDrive/MHWmed_data/SST/1982_SST_MED_SST_L4_REP_OBSERVATIONS_010_021_a_1449092207440.nc") %>% activate("D2")


# Load MME data -----------------------------------------------------------

# Load local spreadsheet
# NB: This is not stored on GitHub as it is not public data
MME <- read_csv("data/Cnidaria_CatalanCoast.csv") %>% 
  dplyr::select(Year:`Damaged qualitative`) %>% 
  distinct() %>% 
  dplyr::rename(lon = Longitude, lat = Latitude) %>% 
  mutate(date_start = case_when(EvenStart == "Autumn" ~ as.Date(paste0(Year,"-09-01")),
                                EvenStart == "Summer" ~ as.Date(paste0(Year,"-06-01"))),
         date_end = case_when(EvenStart == "Autumn" ~ as.Date(paste0(Year,"-11-30")),
                              EvenStart == "Summer" ~ as.Date(paste0(Year,"-08-31")))) %>% 
  mutate(Location = case_when(Location == "Falconera" ~ "Punta Falconera", 
                              Location == "caials" ~ "Caials", 
                              TRUE ~ Location))

# Mean MMEs over Catalan Coast
MME_mean <- MME %>% 
  group_by(Year, date_start, date_end) %>% 
  summarise(damage_mean = round(mean(`Damaged percentage`, na.rm = T), 0),
            count = n(), .groups = "drop")

# Lon/lat values
lon_range_MME <- range(MME$lon)
lat_range_MME <- range(MME$lat)
lon_lat_site <- select(MME, lon, lat, Location) %>% 
  distinct()

# Find nearest SST pixels to MMW data
lon_range_pixels <- med_lon %>% 
  filter(lon >= lon_range_MME[1], lon <= lon_range_MME[2])
lat_range_pixels <- med_lat %>% 
  filter(lat >= lat_range_MME[1], lat <= lat_range_MME[2])


# Load SST data -----------------------------------------------------------

# Extract SST pixels based on spreadsheet lon/lat range
system.time(
SST_prep <- plyr::ldply(med_SST_files, load_nc_range, .parallel = T, 
                        lon_range = range(lon_range_pixels$lon),
                        lat_range = range(lat_range_pixels$lat))
) # 37 seconds

# Match pixels
lon_lat_match <- grid_match(lon_lat_site, distinct(dplyr::select(SST_prep, lon, lat)))


# Detect MHWs per pixel ---------------------------------------------------

# Detect events for each pixel
MHW_res <- SST_prep %>%
  group_by(lon, lat) %>%
  nest() %>% 
  mutate(clim = purrr::map(data, ts2clm, climatologyPeriod = c("1982-01-01", "2011-12-31")),
         event = purrr::map(clim, detect_event, categories = T, climatology = T, season = "peak", S = F)) %>%
  select(-data, -clim) %>% 
  right_join(lon_lat_match, by  = c("lon" = "lon.y", "lat" = "lat.y"))

# Extract MHW clims and events separately for per pixel results
MHW_res_clim <- MHW_res %>% 
  unnest(event) %>% 
  filter(row_number() %% 2 == 1) %>% 
  unnest(event) %>% 
  ungroup() %>% 
  mutate(diff = thresh - seas,
         thresh_2x = thresh + diff,
         thresh_3x = thresh_2x + diff,
         thresh_4x = thresh_3x + diff)
MHW_res_event <- MHW_res %>% 
  unnest(event) %>% 
  filter(row_number() %% 2 == 0) %>% 
  unnest(event) %>% 
  ungroup()


# Detect mean MHWs --------------------------------------------------------

# Create mean SST time series
SST_mean <- SST_prep %>% 
  group_by(t) %>% 
  summarise(temp = mean(temp, na.rm = T))

# Detect mean events
MHW_mean <- detect_event(ts2clm(SST_mean, climatologyPeriod = c("1982-01-01", "2011-12-31")), categories = T, climatology = T)

# Extract clims and events
MHW_mean_clim <- MHW_mean$climatology %>% 
  mutate(diff = thresh - seas,
         thresh_2x = thresh + diff,
         thresh_3x = thresh_2x + diff,
         thresh_4x = thresh_3x + diff)
MHW_mean_event <- MHW_mean$event


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
            duration = sum(duration/pixels, na.rm = T),
            cum_int = sum(cum_int/pixels, na.rm = T))

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
                 label = paste0(round(stats_old$surface/6*100),"%")), 
             label.size = 2, colour = "darkviolet") +
  # New period average
  geom_segment(aes(x = 2014.5, xend = 2019.5, 
                   y = stats_new$surface, yend = stats_new$surface),
               colour = "hotpink", size = 1.5) +
  geom_label(aes(x = 2017, y = stats_new$surface, 
                 label = paste0(round(stats_new$surface/6*100),"%")), 
             label.size = 2, colour = "hotpink") +
  # Other
  scale_y_continuous(breaks = c(0, 3, 6),
                     labels = c("0%", "50%", "100%")) +
  scale_x_continuous(breaks = c(1985, 1995, 2005, 2015)) +
  coord_cartesian(ylim = c(0, 6), expand = F) +
  labs(x = NULL, y = "Surface area covered (%)", fill = "Month",
       title = "Northwestern Mediterranean: surface area covered by MHWs (JJASON)", 
       subtitle = "The stacked bars will add up to 100% if all six months experienced full MHW coverage")
ggsave("talks/graph/NWMed_bar_surface.png", bar_surface, width = 8, height = 5)

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
                 label = paste0(round(stats_old$duration, 0), " days")),
             label.size = 2, colour = "darkviolet") +
  # New period average
  geom_segment(aes(x = 2014.5, xend = 2019.5,
                   y = stats_new$duration, yend = stats_new$duration),
               colour = "hotpink", size = 1.5) +
  geom_label(aes(x = 2017, y = stats_new$duration,
                 label = paste0(round(stats_new$duration, 0), " days")),
             label.size = 2, colour = "hotpink") +
  # Other
  scale_y_continuous(limits = c(0, 183), breaks = c(30, 61, 91, 122, 152)) +
  scale_x_continuous(breaks = c(1985, 1995, 2005, 2015)) +
  coord_cartesian(expand = F) +
  labs(x = NULL, y = "Sum of MHW days", fill = "Month",
       title = "Northwestern Mediterranean: sum of MHW days (JJASON)", 
       subtitle = "If a MHW were to last from June to November it would be 183 days long")
ggsave("talks/graph/NWMed_bar_duration.png", bar_duration, width = 8, height = 5)

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
                 label = paste0(round(stats_old$cum_int, 0), "°C days")),
             label.size = 2, colour = "darkviolet") +
  # New period average
  geom_segment(aes(x = 2014.5, xend = 2019.5,
                   y = stats_new$cum_int, yend = stats_new$cum_int),
               colour = "hotpink", size = 1.5) +
  geom_label(aes(x = 2017, y = stats_new$cum_int,
                 label = paste0(round(stats_new$cum_int, 0), "°C days")),
             label.size = 2, colour = "hotpink") +
  # 100 °C days line
  geom_hline(aes(yintercept = 100), linetype = "dashed", colour = "grey20") +
  # Other
  scale_y_continuous(breaks = c(100, 200)) +
  scale_x_continuous(breaks = c(1985, 1995, 2005, 2015)) +
  coord_cartesian(expand = F) +
  labs(x = NULL, y = "Total temperature anomaly (°C days)", fill = "Month",
       title = "Northwestern Mediterranean: sum of all temperature anomalies (JJASON)", 
       subtitle = "Each bar shows the total temperature anomaly experienced for the warm period per year")
ggsave("talks/graph/NWMed_bar_cum_int.png", bar_cum_int, width = 8, height = 5)



# MME vs MHW plots --------------------------------------------------------

# PLot showing the mean MHWs and MMEs
MME_vs_MHW_mean_plot <- MHW_mean_clim %>% 
  mutate(Year = lubridate::year(t),
         month = lubridate::month(t)) %>% 
  filter(t >= "2015-01-01",
         t <= "2019-12-31",
         # month %in% 6:11
         ) %>% 
  ggplot(aes(x = t, y = temp)) +
  # MHW flames
  geom_flame(aes(y2 = thresh, fill = "Moderate")) +
  geom_flame(aes(y2 = thresh_2x, fill = "Strong")) +
  geom_flame(aes(y2 = thresh_3x, fill = "Severe")) +
  geom_flame(aes(y2 = thresh_4x, fill = "Extreme")) +
  # MHW lines
  geom_line(aes(y = seas), colour = "grey20", size = 0.7) +
  geom_line(aes(y = thresh), colour = "darkgreen", size = 0.7) +
  geom_line(aes(y = temp), colour = "black", size = 0.6) +
  # MME geoms
  geom_segment(data = MME_mean, size = 3,
               aes(x = date_start, xend = date_end,
                   y = 13, yend = 13, colour = damage_mean)) +
  geom_label(data = MME_mean, label.size = 2,
             aes(x = date_start+45, y = 13, colour = damage_mean, 
                 label = paste0("MME = ", count))) +
  # Other
  scale_fill_manual(name = NULL, values = fillColCat, guide = FALSE) +
  scale_colour_gradient(low = "blue", high = "red") +
  scale_x_date(date_labels = "%b %Y", date_breaks = "2 months") +
  scale_y_continuous(limits = c(11, 26.5)) +
  labs(y = "Temperature (°C)", x = NULL, colour = "Damage (%)",
       title = "Catalan Coast: annual time series of MHWs and MMEs",
       subtitle = "The strips at the bottom show the three month season of a recorded MME") +
  coord_cartesian(expand = F) +
  facet_wrap(~Year, scales = "free_x", ncol = 2) +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        legend.position = c(0.7, 0.1),
        legend.direction = "horizontal")
ggsave("talks/graph/MME_vs_MHW_mean_plot.png", height = 10, width = 8)

# Merge MME and MHW results
MME_vs_MHW_plot <- ggplot(data = filter(MHW_res_event,
                                        # Location == "Caials",
                                        date_start >= "2015-01-01"),
                          aes(x = date_peak, y = intensity_cumulative)) +
  geom_point(aes(fill = category), shape = 21, size = 2) +
  geom_segment(data = MME, size = 4,
               aes(x = date_start, xend = date_end,
                   y = 0, yend = 0, colour = `Damaged percentage`)) +
  geom_hline(aes(yintercept = 100), linetype = "dashed", colour = "grey20") +
  scale_fill_manual("Category", values = fillColCatNum) +
  scale_colour_gradient(low = "blue", high = "red") +
  labs(y = "Cumulative intensity (°C)", x = NULL, colour = "Damage (%)",
       title = "Co-occurrence of MME and MHWs") +
  facet_wrap(~Location) +
  theme(legend.position = c(0.65, 0.05),
        legend.direction = "horizontal")
ggsave("talks/graph/MME_vs_MHW_plot.png", height = 10, width = 12)

# Line plot showing MHWs per Location
# Note very useful...
ggplot(data = filter(MHW_res_clim, 
                     # Location == "Caials",
                     t >= "2015-06-01"),
       aes(x = t, y = temp)) +
  geom_flame(aes(y2 = thresh, fill = "Moderate")) +
  geom_flame(aes(y2 = thresh_2x, fill = "Strong")) +
  geom_flame(aes(y2 = thresh_3x, fill = "Severe")) +
  geom_flame(aes(y2 = thresh_4x, fill = "Extreme")) +
  geom_line(aes(y = thresh_2x, col = "2x Threshold"), size = 0.7, linetype = "dashed") +
  geom_line(aes(y = thresh_3x, col = "3x Threshold"), size = 0.7, linetype = "dotdash") +
  geom_line(aes(y = thresh_4x, col = "4x Threshold"), size = 0.7, linetype = "dotted") +
  geom_line(aes(y = seas, col = "Climatology"), size = 0.7) +
  geom_line(aes(y = thresh, col = "Threshold"), size = 0.7) +
  geom_line(aes(y = temp, col = "Temperature"), size = 0.6) +
  scale_colour_manual(name = NULL, values = lineColCat,
                      breaks = c("Temperature", "Climatology", "Threshold",
                                 "2x Threshold", "3x Threshold", "4x Threshold")) +
  scale_fill_manual(name = NULL, values = fillColCat, guide = FALSE) +
  scale_x_date(date_labels = "%b %Y") +
  guides(colour = guide_legend(override.aes = list(linetype = c("solid", "solid", "solid",
                                                                "dashed", "dotdash", "dotted"),
                                                   size = c(0.6, 0.7, 0.7, 0.7, 0.7, 0.7)))) +
  labs(y = "Temperature (°C)", x = NULL) +
  facet_wrap(~Location)