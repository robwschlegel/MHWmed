# code/MHW_summary.R
# This script creates annual and total summaries from the MHW Med results


# Setup -------------------------------------------------------------------

# The project-wide functions
source("code/functions.R")


# Ecoregion summaries -----------------------------------------------------

## NB: Do not run following code in parallel as the functions have internal parallelism
registerDoParallel(cores = 15)

# All pixels
# system.time(
# MHW_cat_region <- plyr::ldply(unique(med_regions$Ecoregion), region_calc, .parallel = F)
# ) # 38 seconds for 1 on 15 cores, ~6.5 minutes total
# save(MHW_cat_region, file = "data/MHW_cat_region.RData")
# readr::write_csv(MHW_cat_region, "data/MHW_cat_region.csv")
load("data/MHW_cat_region.RData")

# Coastal pixels
# system.time(
#   MHW_cat_region_coast <- plyr::ldply(unique(med_regions$Ecoregion), region_calc, pixel_sub = "coast", .parallel = F)
# ) # ~6 minutes total
# save(MHW_cat_region_coast, file = "data/MHW_cat_region_coast.RData")
# readr::write_csv(MHW_cat_region_coast, "data/MHW_cat_region_coast.csv")
load("data/MHW_cat_region_coast.RData")

# MME pixels
# NB: This filters out pixels with "No" damage and uses selected_4 rows
# system.time(
#   MHW_cat_region_pixel <- plyr::ldply(unique(med_regions$Ecoregion), region_calc, pixel_sub = "pixel", .parallel = F)
# ) # ~1 minutes total
# save(MHW_cat_region_pixel, file = "data/MHW_cat_region_pixel.RData")
# readr::write_csv(MHW_cat_region_pixel, "data/MHW_cat_region_pixel.csv")
load("data/MHW_cat_region_pixel.RData")


# Ecoregion summary figures -----------------------------------------------

# plyr::l_ply(unique(med_regions$region), ecoregion_summary_fig, .parallel = T)


# Ecoregion trend figures -------------------------------------------------

# plyr::l_ply(unique(med_regions$region), ecoregion_trend_fig, .parallel = T)


# Ecoregion pixel figures -------------------------------------------------
# These figures compare the MHW stats for different Ecoregions
# based on the different ways of filtering the pixels used

# NB: Requires - MHW_cat_region, MHW_cat_region_coast, MHW_cat_region_pixel

# Prep and combine files
# MHW_cat_region$sub <- "Ecoregion"
# MHW_cat_region_coast$sub <- "Coast"
# MHW_cat_region_pixel$sub <- "MME pixels"
# MHW_cat_region_all <- rbind(MHW_cat_region, MHW_cat_region_coast, MHW_cat_region_pixel) %>% 
#   mutate(sub = factor(sub, levels = c("Ecoregion", "Coast", "MME pixels")))

# Boxplots of statistics
# eco_box_plot <- MHW_cat_region_all %>% 
#   filter(as.numeric(month) %in% 6:11, year >= 2015) %>% 
#   mutate(region = factor(region,
#                          levels = c("Alboran Sea", "Northwestern Mediterranean", 
#                                     "Southwestern Mediterranean", "Adriatic Sea",
#                                     "Ionian Sea", "Tunisian Plateau/Gulf of Sidra",
#                                     "Aegean Sea", "Levantine Sea")),
#          dur_prop = duration/pixels) %>% 
#   select(region, sub, surface, dur_prop, mean_int) %>% 
#   pivot_longer(surface:mean_int) %>% 
#   mutate(name = case_when(name == "dur_prop" ~ "MHW days (n)",
#                           name == "mean_int" ~ "Mean intensity (°C)",
#                           name == "surface" ~ "Surface area (%)",
#                           TRUE ~ name),
#          name = factor(name, levels = c("Surface area (%)", "MHW days (n)", "Mean intensity (°C)"))) %>% 
#   ggplot(aes(x = region, y = value)) +
#   geom_boxplot(aes(fill = region, linetype = sub)) +
#   coord_flip() +
#   guides(fill = FALSE) +
#   # scale_x_reverse() +
#   facet_wrap(~name, scales = "free_x", nrow = 1, strip.position = "bottom") +
#   labs(x = NULL, y = NULL, linetype = "subset",
#        title = "Boxplots of MHW values for JJASON months from 2015-2019",
#        subtitle = "Different outlined boxplots show subsets of MHW values for whole Ecoregion (solid), coast only (dotted), or MME pixels only (dashed)") +
#   theme(legend.position = "bottom")
# ggsave("figures/MHW_eco_pixels_boxplot.png", eco_box_plot, height = 8, width = 14)


# Map summary figures -----------------------------------------------------

# plyr::l_ply(2015:2019, monthly_map_fig_full, .parallel = T)


# Annual summaries --------------------------------------------------------

# The occurrences per month per pixel
doParallel::registerDoParallel(cores = 15)
# system.time(
# MHW_cat_pixel_monthly <- plyr::ldply(res_files, cat_pixel_calc, .parallel = T)
# ) # 258 seconds on 15 cores
# save(MHW_cat_pixel_monthly, file = "data/MHW_cat_pixel_monthly.RData")
# load("data/MHW_cat_pixel_monthly.RData") # This is very large, only load if necessary

# The occurrences per year per pixel
## NB: Requires MHW_cat_pixel_monthly
# system.time(
# MHW_cat_pixel_annual <- cat_pixel_annual_calc()
# ) # 220 seconds
# save(MHW_cat_pixel_annual, file = "data/MHW_cat_pixel_annual.RData")
load("data/MHW_cat_pixel_annual.RData")

# The occurrences per year per pixel JJASON
## NB: Requires MHW_cat_pixel_monthly
# system.time(
# MHW_cat_pixel_annual_JJASON <- cat_pixel_annual_calc(sub_months = seq(6, 11))
# ) # 197 seconds
# save(MHW_cat_pixel_annual_JJASON, file = "data/MHW_cat_pixel_annual_JJASON.RData")
load("data/MHW_cat_pixel_annual_JJASON.RData")

# The occurrences per day all months
# system.time(
# MHW_cat_daily_annual <- plyr::ldply(res_files, cat_daily_calc, .parallel = T)
# ) # 170 seconds on 15 cores
# save(MHW_cat_daily_annual, file = "data/MHW_cat_daily_annual.RData")
load("data/MHW_cat_daily_annual.RData")

# The occurrences per day JJASON
# system.time(
# MHW_cat_daily_annual_JJASON <- plyr::ldply(res_files, cat_daily_calc, .parallel = T, sub_months = seq(6, 11))
# ) # 164 seconds on 15 cores
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
doParallel::registerDoParallel(cores = 15)
# plyr::l_ply(1982:2019, annual_summary_fig, .parallel = T)

# Create total summary figure
total_summary <- total_summary_fig(MHW_cat_summary_annual)
ggsave("figures/MHW_cat_historic.png", total_summary, height = 4.25, width = 8)
total_summary_JJASON <- total_summary_fig(MHW_cat_summary_annual_JJASON)
ggsave("figures/MHW_cat_historic_JJASON.png", total_summary_JJASON, height = 4.25, width = 8)

# Subset data for following two figures
MHW_cat_pixel_annual_sub <- MHW_cat_pixel_annual %>%
  filter(year %in% seq(2015, 2019))

# Med maps of duration summary values
med_map_dur <- MHW_cat_pixel_annual_sub %>%
  group_by(lon, lat) %>%
  summarise(duration_sum = sum(duration_sum, na.rm = T), .groups = "drop") %>%
  ggplot(aes(x = lon, y = lat)) +
  # geom_tile(data = OISST_ice_coords, fill = "powderblue", colour = NA, alpha = 0.5) +
  geom_tile(aes(fill = duration_sum), colour = NA) +
  geom_polygon(data = map_base, aes(x = lon, y = lat, group = group),
               fill = "grey70", colour = "black") +
  scale_fill_distiller("Days", palette = "Greens", direction = 1) +
  scale_y_continuous(breaks = NULL) +
  scale_x_continuous(breaks = NULL) +
  coord_cartesian(expand = F,
                  xlim = c(min(med_regions$lon), max(med_regions$lon)),
                  ylim = c(min(med_regions$lat), max(med_regions$lat))) +
  theme_bw() +
  labs(title = "__C)__    Sum of MHW days from 2015-2019", x = NULL, y = NULL) +
  # guides(fill = guide_legend(override.aes = list(size = 10))) +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        plot.title = ggtext::element_markdown(),
        legend.position = c(0.92, 0.80),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        panel.background = element_rect(fill = "grey90"),
        axis.text = element_blank(),
        axis.ticks = element_blank())
med_map_dur

# Med maps of duration summary values
med_map_cat <- MHW_cat_pixel_annual_sub %>%
  group_by(lon, lat) %>%
  summarise(category = max(as.numeric(category), na.rm = T), .groups = "drop") %>%
  mutate(category = factor(category, labels = c("I Moderate", "II Strong", "III Severe", "IV Extreme"))) %>%
  ggplot(aes(x = lon, y = lat)) +
  # geom_tile(data = OISST_ice_coords, fill = "powderblue", colour = NA, alpha = 0.5) +
  geom_tile(aes(fill = category), colour = NA, show.legend = F) +
  geom_polygon(data = map_base, aes(x = lon, y = lat, group = group),
               fill = "grey70", colour = "black") +
  scale_fill_manual("Category", values = MHW_colours) +
  scale_y_continuous(breaks = NULL) +
  scale_x_continuous(breaks = NULL) +
  coord_cartesian(expand = F,
                  xlim = c(min(med_regions$lon), max(med_regions$lon)),
                  ylim = c(min(med_regions$lat), max(med_regions$lat))) +
  # theme_void() +
  theme_bw() +
  labs(title = "__D)__   Highest MHW categories from 2015-2019", x = NULL, y = NULL) +
  # guides(fill = guide_legend(override.aes = list(size = 10))) +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        plot.title = ggtext::element_markdown(),
        legend.position = "bottom",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        panel.background = element_rect(fill = "grey90"),
        axis.text = element_blank(),
        axis.ticks = element_blank())
med_map_cat

# Combine with summary bar plots
med_map_combi <- ggpubr::ggarrange(med_map_dur, med_map_cat, nrow = 1)
total_summary_quad <- ggpubr::ggarrange(total_summary, med_map_combi, ncol = 1)
ggsave("figures/MHW_cat_historic_quad.png", total_summary_quad, height = 9, width = 10)


# Per pixel maps with MME -------------------------------------------------

# Requires: MHW_cat_pixel_monthly.RData and MHW_cat_region.RData

# Per pixel maps
# map_pixel_duration <- monthly_map_pixel("duration_sum", annual = T)
# ggsave("figures/MHW_pixel_duration.png", map_pixel_duration, height = 7, width = 20)
# map_pixel_category <- monthly_map_pixel("category", annual = T)
# ggsave("figures/MHW_pixel_category.png", map_pixel_category, height = 7, width = 20)
# map_pixel_cum_int <- monthly_map_pixel("cum_int", annual = T)
# ggsave("figures/MHW_pixel_cum_int.png", map_pixel_cum_int, height = 7, width = 20)


# Single summary map ------------------------------------------------------

# Requires: MHW_cat_pixel_monthly, MHW_cat_pixel_annual_JJASON
# load("data/MHW_cat_pixel_monthly.RData")
# load("data/MHW_cat_pixel_annual_JJASON.RData")

# Prepare MME points
# mme_damage <- mme_selected_4 %>% 
#   filter(`Damaged qualitative` != "No") %>% 
#   group_by(lon, lat) %>% 
#   summarise(count = n(), .groups = "drop")
# mme_no <- mme_selected_4 %>% 
#   filter(`Damaged qualitative` == "No") %>% 
#   group_by(lon, lat) %>% 
#   summarise(count = n(), .groups = "drop")

# Determine historic medians per pixel
# MHW_cat_pixel_JJASON_clim_median <- MHW_cat_pixel_annual_JJASON %>% 
#   ungroup() %>% 
#   filter(year %in% seq(1982, 2014)) %>%
#   group_by(lon, lat) %>% 
#   summarise(duration_median = median(duration_sum, na.rm = T),
#             icum_median = median(cum_int, na.rm = T), .groups = "drop")

# Create study period anomalies
# MHW_cat_pixel_JJASON_anom <- MHW_cat_pixel_annual_JJASON %>% 
#   ungroup() %>% 
#   right_join(med_regions, by = c("lon", "lat")) %>% 
#   filter(year %in% seq(2015, 2019)) %>%
#   group_by(lon, lat) %>% 
#   summarise(duration = median(duration_sum, na.rm = T),
#             icum = median(cum_int, na.rm = T), .groups = "drop") %>% 
#   left_join(MHW_cat_pixel_JJASON_clim_median, by = c("lon", "lat")) %>% 
#   mutate(duration_anom = duration-duration_median,
#          icum_anom = icum-icum_median) %>% 
#   mutate(duration_anom = ifelse(duration_anom < 0, 0, duration_anom),
#          icum_anom = ifelse(icum_anom < 0, 0, icum_anom))
  
# Map of duration anom for JJASON
# anom_plot_dur <- med_base + 
#   geom_tile(data = MHW_cat_pixel_JJASON_anom, colour = NA,
#             aes(fill = duration_anom, x = lon, y = lat)) +
#   # geom_contour(data = MHW_cat_pixel_JJASON_anom, colour = "black",
#                # aes(x = lon, y = lat, z = duration_anom), breaks = c(1)) +
#   geom_sf(data = MEOW, alpha = 1, aes(geometry = geometry), fill = NA, colour = "grey70") +
#   scale_fill_gradient2(low = "white", high = "forestgreen") +
#   coord_sf(expand = F, xlim = c(-10, 45), ylim = c(25, 50)) +
#   labs(fill = "Duration (days)",
#        title = "MHW Duration",
#        subtitle = "Median per year for 2015-2019 JJASON in excess of the 1982-2014 JJASON median") +
#   theme(panel.border = element_rect(colour = "black", fill = NA),
#         legend.position = "bottom",
#         title = element_text(size = 18),
#         legend.text = element_text(size = 16),
#         legend.title = element_text(size = 18),
#         panel.background = element_rect(fill = "grey90"), 
#         strip.text = element_text(size = 16))
# anom_plot_dur

# Map of icum anom for JJASON
# anom_plot_icum <- med_base + 
#   geom_tile(data = MHW_cat_pixel_JJASON_anom, colour = NA,
#             aes(fill = icum_anom, x = lon, y = lat)) +
#   # geom_contour(data = MHW_cat_pixel_JJASON_anom, colour = "black",
#                # aes(x = lon, y = lat, z = icum_anom), breaks = c(1)) +
#   geom_sf(data = MEOW, alpha = 1, aes(geometry = geometry), fill = NA, colour = "grey70") +
#   scale_fill_gradient2(low = "white", high = "darkorchid") +
#   coord_sf(expand = F, xlim = c(-10, 45), ylim = c(25, 50)) +
#   labs(fill = "Cumulative\nintensity (°C days)", 
#        title = "MHW Cumulative Intensity",
#        subtitle = "Median per year for 2015-2019 JJASON in excess of the 1982-2014 JJASON median") +
#   theme(panel.border = element_rect(colour = "black", fill = NA),
#         legend.position = "bottom",
#         title = element_text(size = 18),
#         legend.text = element_text(size = 16),
#         legend.title = element_text(size = 18),
#         panel.background = element_rect(fill = "grey90"), 
#         strip.text = element_text(size = 16))
# anom_plot_icum

# MME damage
# anom_plot_mme <- med_base + 
#   geom_sf(data = MEOW, alpha = 1, aes(geometry = geometry), fill = NA, colour = "grey70") +
#   geom_point(data = mme_damage, aes(x = lon, y = lat, size = count), 
#              alpha = 0.7, shape = 21, fill = "yellow", colour = "red") +
#   coord_sf(expand = F, xlim = c(-10, 45), ylim = c(25, 50)) +
#   labs(size = "MME count (n)",
#        title = "MME damage",
#        subtitle = "Recorded damaging events for 2015-2019 Summer/Autumn") +
#   theme(panel.border = element_rect(colour = "black", fill = NA),
#         legend.position = "bottom",
#         title = element_text(size = 18),
#         legend.text = element_text(size = 16),
#         legend.title = element_text(size = 18),
#         panel.background = element_rect(fill = "grey90"), 
#         strip.text = element_text(size = 16))
# anom_plot_mme

# MME damage
# anom_plot_mme_no <- med_base + 
#   geom_sf(data = MEOW, alpha = 1, aes(geometry = geometry), fill = NA, colour = "grey70") +
#   geom_point(data = mme_no, aes(x = lon, y = lat, size = count), 
#              alpha = 0.7, shape = 21, fill = "yellow", colour = "red") +
#   coord_sf(expand = F, xlim = c(-10, 45), ylim = c(25, 50)) +
#   labs(size = "MME count (n)",
#        title = "MME no damage",
#        subtitle = "Recorded non-damaging events for 2015-2019 Summer/Autumn") +
#   theme(panel.border = element_rect(colour = "black", fill = NA),
#         legend.position = "bottom",
#         title = element_text(size = 18),
#         legend.text = element_text(size = 16),
#         legend.title = element_text(size = 18),
#         panel.background = element_rect(fill = "grey90"), 
#         strip.text = element_text(size = 16))
# anom_plot_mme_no

# Combine and save
# anom_all <- ggpubr::ggarrange(anom_plot_dur, anom_plot_icum,# anom_plot_mme, anom_plot_mme_no,
#                               ncol = 2, nrow = 1, align = "hv")
# ggsave("figures/MHW_pixel_median_anom.png", anom_all, height = 8, width = 22)


# MME vs MHW pixels -------------------------------------------------------

# Match pixels
lon_lat_match <- grid_match(distinct(dplyr::select(mme_selected_4, lon, lat)), 
                            distinct(dplyr::select(med_regions, lon, lat))) %>% 
  dplyr::rename(lon_mme = lon.x, lat_mme = lat.x, lon_sst = lon.y, lat_sst = lat.y) %>% 
  left_join(med_regions, by = c("lon_sst" = "lon", "lat_sst" = "lat")) %>% 
  filter(dist < 10)

# Create full annual grid
lon_lat_match_full_grid <- expand_grid(year = seq(2015, 2019), lon_lat_match) 

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
site_MME_summary <- mme_selected_4 %>% 
  group_by(Ecoregion, year, `Monitoring series`, `Damaged qualitative`, EvenStart, Location, lon, lat) %>% 
  summarise(`Damaged percentage` = round(mean(`Damaged percentage`, na.rm = T)),
            count_MME = n(), .groups = "drop") %>% 
  dplyr::rename(lon_mme = lon, lat_mme = lat)

# Join them together
site_MME_MHW_summary <- lon_lat_match_full_grid %>% 
  left_join(site_MHW_summary, by = c("lon_sst", "lat_sst", "year")) %>% 
  replace(is.na(.), 0) %>%  # Fill in the no MHW rows with 0's
  full_join(site_MME_summary, by = c("lon_mme", "lat_mme", "Ecoregion", "year")) %>% 
  filter(year >= 2015)
write_csv(site_MME_MHW_summary, "data/site_MME_MHW_summary.csv")

# Scatterplot of MME and MHW count summaries
# scatter_MME_MHW_site <- site_MME_MHW_summary %>% 
#   na.omit() %>% 
#   pivot_longer(count_MHW:icum) %>% 
#   ggplot(aes(x = value, y = count_MME)) +
#   geom_point(aes(colour = Ecoregion)) +
#   geom_smooth(aes(colour = Ecoregion), method = "lm", se = F) +
#   labs(y = "MME count per site/year", x = NULL,
#        title = "MME count per site and year compared to MHW metrics (JJASON)") +
#   facet_wrap(~name, scales = "free_x", strip.position = "bottom") +
#   theme(legend.position = "bottom", 
#         strip.placement = "outside", strip.background = element_blank())
# scatter_MME_MHW_site

# Scatterplot of MME and MHW damage summaries
# scatter_MME_MHW_site_dam <- site_MME_MHW_summary %>% 
#   na.omit() %>% 
#   pivot_longer(count_MHW:icum) %>% 
#   ggplot(aes(x = value, y = `Damaged percentage`)) +
#   geom_point(aes(colour = Ecoregion)) +
#   geom_smooth(aes(colour = Ecoregion), method = "lm", se = F) +
#   labs(y = "MME damage per site/year", x = NULL,
#        title = "MME damage per site and year compared to MHW metrics (JJASON)") +
#   facet_wrap(~name, scales = "free_x", strip.position = "bottom") +
#   theme(legend.position = "bottom", 
#         strip.placement = "outside", strip.background = element_blank())
# scatter_MME_MHW_site_dam

# Scatterplot of MME and MHW count summaries by ecoregion
# scatter_MME_MHW_ecoregion_sum <- site_MME_MHW_summary %>% 
#   na.omit() %>% 
#   group_by(Ecoregion, year) %>% 
#   summarise(`Damaged percentage` = mean(`Damaged percentage`),
#             count_MME_sum = sum(count_MME),
#             count_MME_mean = mean(count_MME),
#             count_MHW_sum = sum(count_MHW),
#             # count_MHW_mean = mean(count_MHW),
#             duration = mean(duration),
#             # imean = mean(imean),
#             icum = mean(icum), .groups = "drop") %>% 
#   pivot_longer(count_MHW_sum:icum) %>% 
#   ggplot(aes(x = value, y = count_MME_sum)) +
#   geom_point(aes(colour = Ecoregion)) +
#   geom_smooth(aes(colour = Ecoregion), method = "lm", se = F) +
#   labs(y = "MME count per ecoregion/year (sum)", x = NULL,
#        title = "MME count (sum) per ecoregion and year compared to MHW metrics (JJASON)") +
#   facet_wrap(~name, scales = "free_x", strip.position = "bottom") +
#   theme(legend.position = "bottom", 
#         strip.placement = "outside", strip.background = element_blank())
# scatter_MME_MHW_ecoregion_sum

# Scatterplot of MME and MHW damage summaries by ecoregion
# scatter_MME_MHW_ecoregion_sum_dam <- site_MME_MHW_summary %>% 
#   na.omit() %>% 
#   group_by(Ecoregion, year) %>% 
#   summarise(`Damaged percentage` = mean(`Damaged percentage`),
#             count_MME_sum = sum(count_MME),
#             count_MME_mean = mean(count_MME),
#             count_MHW_sum = sum(count_MHW),
#             # count_MHW_mean = mean(count_MHW),
#             duration = mean(duration),
#             # imean = mean(imean),
#             icum = mean(icum), .groups = "drop") %>% 
#   pivot_longer(count_MHW_sum:icum) %>% 
#   ggplot(aes(x = value, y = `Damaged percentage`)) +
#   geom_point(aes(colour = Ecoregion)) +
#   geom_smooth(aes(colour = Ecoregion), method = "lm", se = F) +
#   labs(y = "MME damage (%) per ecoregion/year (sum)", x = NULL,
#        title = "MME damage per ecoregion and year compared to MHW metrics (JJASON)") +
#   facet_wrap(~name, scales = "free_x", strip.position = "bottom") +
#   theme(legend.position = "bottom", 
#         strip.placement = "outside", strip.background = element_blank())
# scatter_MME_MHW_ecoregion_sum_dam

# Scatterplot of MME and MHW count summaries by ecoregion
# scatter_MME_MHW_ecoregion_mean <- site_MME_MHW_summary %>% 
#   na.omit() %>% 
#   group_by(Ecoregion, year) %>% 
#   summarise(`Damaged percentage` = mean(`Damaged percentage`),
#             count_MME_sum = sum(count_MME),
#             count_MME_mean = mean(count_MME),
#             # count_MHW_sum = sum(count_MHW),
#             count_MHW_mean = mean(count_MHW),
#             duration = mean(duration),
#             # imean = mean(imean),
#             icum = mean(icum), .groups = "drop") %>% 
#   pivot_longer(count_MHW_mean:icum) %>% 
#   ggplot(aes(x = value, y = count_MME_mean)) +
#   geom_point(aes(colour = Ecoregion)) +
#   geom_smooth(aes(colour = Ecoregion), method = "lm", se = F) +
#   labs(y = "MME count per ecoregion/year (mean)", x = NULL,
#        title = "MME count per ecoregion and year compared to MHW metrics (JJASON)") +
#   facet_wrap(~name, scales = "free_x", strip.position = "bottom") +
#   theme(legend.position = "bottom", 
#         strip.placement = "outside", strip.background = element_blank())
# scatter_MME_MHW_ecoregion_mean

# Scatterplot of MME and MHW damage summaries by ecoregion
# scatter_MME_MHW_ecoregion_mean_dam <- site_MME_MHW_summary %>% 
#   na.omit() %>% 
#   group_by(Ecoregion, year) %>% 
#   summarise(`Damaged percentage` = mean(`Damaged percentage`),
#             count_MME_sum = sum(count_MME),
#             count_MME_mean = mean(count_MME),
#             # count_MHW_sum = sum(count_MHW),
#             count_MHW_mean = mean(count_MHW),
#             duration = mean(duration),
#             # imean = mean(imean),
#             icum = mean(icum), .groups = "drop") %>% 
#   pivot_longer(count_MHW_mean:icum) %>% 
#   ggplot(aes(x = value, y = `Damaged percentage`)) +
#   geom_point(aes(colour = Ecoregion)) +
#   geom_smooth(aes(colour = Ecoregion), method = "lm", se = F) +
#   labs(y = "MME damage (%) per ecoregion/year", x = NULL,
#        title = "MME damage per ecoregion and year compared to MHW metrics (JJASON)") +
#   facet_wrap(~name, scales = "free_x", strip.position = "bottom") +
#   theme(legend.position = "bottom", 
#         strip.placement = "outside", strip.background = element_blank())
# scatter_MME_MHW_ecoregion_mean_dam

# Stitch together and save
# scatter_MME_MHW <- ggpubr::ggarrange(scatter_MME_MHW_site, scatter_MME_MHW_site_dam,
#                                      scatter_MME_MHW_ecoregion_sum, scatter_MME_MHW_ecoregion_sum_dam, 
#                                      scatter_MME_MHW_ecoregion_mean, scatter_MME_MHW_ecoregion_mean_dam, 
#                                      ncol = 2, nrow = 3, align = "hv", common.legend = T, legend = "bottom")
# ggsave("figures/scatter_MME_MHW.png", scatter_MME_MHW, height = 9, width = 16)


# Yes vs No MME -----------------------------------------------------------

# TODO: Find a way to show thresholds above which a MHW metric relates to rapid increase in MME
# e.g. The threshold of 3 MHW vs 2 MHW
# Scatterplots with icum for all species and global
# Also ecoregions

# Calculate the r and p values

# Load MME MHW pairing data
site_MME_MHW_summary <- read_csv("data/site_MME_MHW_summary.csv")

# Load species grouping sheet
species_groups <- read_csv("data/MME_MHWs_relationship_species_selection.csv") %>% 
  `colnames<-`(c("species", "damage", "group", "group_single"))

# Extract only records with regular monitoring
mme_reg <- mme_selected_4 %>% 
  filter(`Monitoring series` %in% c("more.than.two.per.year", "one.per.year.monitoring")) %>% 
  left_join(site_MME_MHW_summary, by = c("lon" = "lon_mme", "lat" = "lat_mme",
                                         "year", "Ecoregion", "Location")) %>% 
  dplyr::rename(`Damaged percentage` = `Damaged percentage.x`,
                `Damaged percentage (mean)` = `Damaged percentage.y`)

# List of grouped species
spp_1 <- filter(species_groups, group == "1")
spp_2 <- filter(species_groups, group == "2")
spp_3 <- filter(species_groups, group == "3")

# Filter out regularly montiored sites by species groups
mme_reg_1 <-  filter(mme_reg, Species %in% spp_1$species)
mme_reg_2 <-  filter(mme_reg, Species %in% spp_2$species)
mme_reg_3 <-  filter(mme_reg, Species %in% spp_3$species)

# Paramuricea clavata is perhaps the best candidate for heat stress
mme_reg_Pcla <- filter(mme_reg, Species == "Paramuricea clavata")#,
         # `Damaged qualitative` != "No")

# Scatterplot of Paramuricea clavata MME vs MHW per pixel
scatter_Pcla <- species_scatter(mme_reg_Pcla, "Paramuricea clavata ")
ggsave("figures/scatter_Pcla.png", scatter_Pcla, height = 6, width = 9)

# Scatterplot for group 1 species
scatter_spp_1 <- species_scatter(mme_reg_1, "Group 1 species ")
ggsave("figures/scatter_spp_1.png", scatter_spp_1, height = 6, width = 9)

# Scatterplot for group 1+2 species
scatter_spp_1_2 <- species_scatter(rbind(mme_reg_1, mme_reg_2), "Group 1+2 species ")
ggsave("figures/scatter_spp_1_2.png", scatter_spp_1_2, height = 6, width = 9)

# Scatterplot for group 1+2 species
scatter_spp_1_2_3 <- species_scatter(rbind(mme_reg_1, mme_reg_2, mme_reg_3), "Group 1+2+3 species ")
ggsave("figures/scatter_spp_1_2_3.png", scatter_spp_1_2_3, height = 6, width = 9)

# Scatterplot for all species
scatter_spp_all <- species_scatter(mme_reg, "All species ")
ggsave("figures/scatter_spp_all.png", scatter_spp_all, height = 6, width = 9)


# Spatial MME vs MHW maps -------------------------------------------------

# MHW days may be the easiest to communicate
# °C is important, but may not be as accessible
# May need to have multiple maps. Faceted in one figure.

# MME filtering
# Definitely run comparisons with MHW at the smallest scale
# Pixel level at the nearest sites with MME records
# MMEs on the surface should be 15 or shallower
# The regularly monitored MME sites are going to be labeled in the main spreadsheet over the week
# For comparisons try the multiple different temporal levels: 
  # all 5 years, per year, per season per year, whole period per ecoregion
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
# Certain species are better choices than others: 
# Paramuricea clavata, Eunicella singularis, Eunicella cavolini, Corallium rubrum
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

# Calculate annual SST anomaly for globe
# Get the linear trend and the 2015-2019 period


# MHW metric time series and MME rug plot ---------------------------------



# Histograms of MME and MHW -----------------------------------------------

# Requires: site_MME_MHW_summary
site_MME_MHW_summary <- read_csv("data/site_MME_MHW_summary.csv")

# Complete region/year grid
full_region_grid <- expand_grid(Ecoregion = unique(mme_selected_4$Ecoregion))
full_region_year_grid <- expand_grid(year = seq(2015, 2019), 
                                     Ecoregion = unique(mme_selected_4$Ecoregion))

# Prepare MME points
mme_points <- mme_selected_4 %>% 
  filter(`Monitoring series` %in% c("more.than.two.per.year", "one.per.year.monitoring"),
         # `Upper Depth` <= 15,
         `Damaged qualitative` != "No")

# Regularly sampled sites
mme_sites_unique <- mme_points %>% 
  dplyr::select(Ecoregion, lon, lat) %>% 
  unique()

# Average MHW stats per regularly sampled site
mhw_sites_average <- site_MME_MHW_summary %>% 
  group_by(Ecoregion, year) %>% 
  right_join(mme_sites_unique, by = c("Ecoregion", "lon_mme" = "lon", "lat_mme" = "lat")) %>% 
  summarise(sites_mhw = n(), 
            duration = mean(duration, na.rm = T),
            icum = mean(icum, na.rm = T), .groups = "drop") %>% 
  right_join(full_region_year_grid, by = c("Ecoregion", "year")) %>% 
  replace(is.na(.), 0)

# Unique sites per ecoregion
mme_sites <- site_MME_MHW_summary %>% 
  filter(`Monitoring series` %in% c("more.than.two.per.year", "one.per.year.monitoring")) %>% 
  group_by(Ecoregion, year, lon_mme, lat_mme) %>%
  # dplyr::select(Ecoregion, lon, lat) %>%
  # unique(.) %>%
  summarise(sites_mme = n(), .groups = "drop") %>%
  right_join(full_region_year_grid, by = c("Ecoregion", "year")) %>%
  replace(is.na(.), 0) %>% 
  group_by(Ecoregion, year) %>%
  summarise(sites_mme = sum(sites_mme), .groups = "drop")
  
# Prepare MME labels
# mme_labels <- mme_points %>% 
#   group_by(year, Ecoregion) %>% 
#   summarise(mme_count = n(), .groups = "drop") %>% 
#   right_join(full_region_year_grid, by = c("year", "Ecoregion")) %>% 
#   mutate(mme_count = ifelse(is.na(mme_count), 0, mme_count)) %>% 
#   # left_join(mme_sites, by = "Ecoregion") %>% 
#   # mutate(mme_count_prop = round(mme_count/site_count, 1)) %>% 
#   replace(is.na(.), 0) %>% 
#   arrange(year, Ecoregion)

# Filter data to target time period
ecoregion_MME_MHW <- site_MME_MHW_summary %>% 
  filter(`Monitoring series` %in% c("more.than.two.per.year", "one.per.year.monitoring"),
         `Damaged qualitative` != "No") %>%
  # right_join(mme_sites_unique, by = c("Ecoregion", "lon_mme" = "lon", "lat_mme" = "lat")) %>% 
  # group_by(lon_mme, lat_mme) %>% 
  # ungroup() %>% 
  # dplyr::select(Ecoregion, year, count_MHW:count_MME, -Location) %>% 
  group_by(Ecoregion, year) %>% 
  summarise(`Damaged percentage` = mean(`Damaged percentage`, na.rm = T),
            # sites_mme = n(),
            count_MME_sum = sum(count_MME, na.rm = T),
            # count_MHW_sum = sum(count_MHW, na.rm = T),
            # count_MHW_mean = mean(count_MHW),
            # duration = mean(duration, na.rm = T),
            # imean = mean(imean),
            # icum = mean(icum, na.rm = T), 
            .groups = "drop") %>% 
  right_join(full_region_year_grid, by = c("Ecoregion", "year")) %>% 
  # mutate(count_MME_mean = round(count_MME_sum/sites_mme, 1)) %>% 
  left_join(mhw_sites_average, by = c("Ecoregion", "year")) %>% 
  left_join(mme_sites, by = c("Ecoregion", "year")) %>% 
  # right_join(mme_labels, by = c("year", "Ecoregion")) %>%
  replace(is.na(.), 0)

# Create a whole Med summary
ecoregion_MME_MHW_med <- ecoregion_MME_MHW %>% 
  group_by(year) %>% 
  summarise(`Damaged percentage` = mean(`Damaged percentage`, na.rm = T),
            count_MME_sum = sum(count_MME_sum, na.rm = T),
            sites_mme = sum(sites_mme, na.rm = T),
            sites_mhw = sum(sites_mhw, na.rm = T),
            duration = mean(duration, na.rm = T),
            icum = mean(icum, na.rm = T), .groups = "drop") %>% 
  mutate(Ecoregion = "Mediterranean")

# Combine and order factor for plotting
ecoregion_MME_MHW_all <- rbind(ecoregion_MME_MHW, ecoregion_MME_MHW_med) %>% 
  mutate(Ecoregion = factor(Ecoregion, 
                            levels = c("Mediterranean",
                                       "Alboran Sea", "Northwestern Mediterranean", 
                                       "Southwestern Mediterranean", "Adriatic Sea",
                                       "Ionian Sea", "Tunisian Plateau/Gulf of Sidra",
                                       "Aegean Sea", "Levantine Sea")))

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
  # geom_label(aes(label = count_MME_mean)) +
  geom_label(aes(label = paste0(count_MME_sum,"/",sites_mme))) +
  scale_fill_viridis_c("Cumulative\nIntensity (°C days)", option = "B") +
  facet_wrap(~Ecoregion) +
  scale_y_continuous(limits = c(-5, 125), breaks = c(25, 50, 75, 100)) +
  scale_x_continuous(breaks = c(2015, 2017, 2019)) +
  coord_cartesian(expand = F) +
  labs(x = NULL, y = "MHW days", 
       title = "Total MHW days and cumulative intensity (°C days) per ecoregion/year (JJASON)",
       subtitle = "Labels show count of MME with damage from regular monitoring at all depths: count/sites") +
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
# bar_dur
ggsave("figures/MHW_ecoregion_summary.png", bar_dur, height = 8, width = 8)

