# code/figures.R
# This script creates annual and total summaries from the MHW Med results


# Setup -------------------------------------------------------------------

# The project-wide functions
source("code/functions.R")


# Figure 1 ----------------------------------------------------------------

# Load med and global annual MHW summary stats
load("data/MHW_cat_summary_annual.RData")
OISST_global <- readRDS("data/OISST_cat_daily_1992-2018_total.Rds") %>% 
  group_by(t) %>% 
  mutate(cat_n_prop_stack = cumsum(cat_n_prop),
         first_n_cum_prop_stack = cumsum(first_n_cum_prop)) %>% 
  ungroup()

# Total annual Med MHW summary with global overlay
total_summary <- total_summary_fig(MHW_cat_summary_annual)
ggsave("figures/fig_1.png", total_summary, height = 4.25, width = 8)

# Prep data for stats
MHW_cat_summary <- MHW_cat_summary_annual %>% 
  group_by(t) %>% 
  mutate(cat_n_prop_stack = cumsum(cat_n_prop),
         first_n_cum_prop_stack = cumsum(first_n_cum_prop)) %>% 
  ungroup() %>% 
  filter(category == "IV Extreme",
         grepl("12-31", as.character(t)))
OISST_global_summary <- OISST_global %>% 
  filter(category == "IV Extreme")

# MHW stats
lm(cat_n_prop_stack ~ year, MHW_cat_summary)
0.006095*365
lm(cat_n_prop_stack ~ t, OISST_global_summary)


# Figure 2 ----------------------------------------------------------------

# Barplots showing duration and MME
# Use coastal pixels

# Complete region/year grid
full_region_grid <- expand_grid(Ecoregion = unique(mme_selected_4$Ecoregion))
full_region_year_grid <- expand_grid(year = seq(2015, 2019), 
                                     Ecoregion = unique(mme_selected_4$Ecoregion))
join_vals <- c("Ecoregion", "year")

# Load coastal ecoregion summaries
load("data/MHW_cat_region_coast.RData")

# MHW JJASON stats
mhw_eco_summary_JJASON <- MHW_cat_region_coast %>% 
  dplyr::rename(Ecoregion = region) %>% 
  filter(as.numeric(month) %in% 6:11) %>% 
  dplyr::select(-month) %>% 
  mutate(duration = duration/pixels,
         cum_int = cum_int/pixels) %>% 
  group_by(Ecoregion, year) %>% 
  summarise(surface = mean(surface, na.rm = T),
            duration = sum(duration, na.rm = T),
            icum = sum(cum_int, na.rm = T), .groups = "drop")

# Summary of MME per region per year
mme_eco_summary <- mme_selected_4 %>% 
  filter(`Damaged qualitative` != "No") %>% 
  group_by(Ecoregion, year) %>% 
  summarise(`Damaged percentage` = round(mean(`Damaged percentage`, na.rm = T)),
            mme_count = n(), .groups = "drop")

# Regularly sampled sites
mme_eco_sites_unique <- mme_selected_4 %>% 
  filter(`Damaged qualitative` != "No") %>% 
  dplyr::select(Ecoregion, year, lon, lat) %>% 
  unique()

# Count of sites per year
mme_eco_sites_count <- mme_eco_sites_unique %>% 
  group_by(Ecoregion, year) %>% 
  summarise(site_count = n(), .groups = "drop")

# Merge eco results
eco_MME_MHW_JJASON <- left_join(full_region_year_grid, mme_eco_summary, by = join_vals) %>% 
  left_join(mme_eco_sites_count, by = join_vals) %>% 
  left_join(mhw_eco_summary_JJASON, by = join_vals) %>% 
  replace(is.na(.), 0)

# Create a whole Med summary
eco_MME_MHW_med_JJASON <- eco_MME_MHW_JJASON %>% 
  group_by(year) %>% 
  summarise(`Damaged percentage` = mean(`Damaged percentage`, na.rm = T),
            mme_count = sum(mme_count, na.rm = T),
            site_count = sum(site_count, na.rm = T),
            surface = mean(surface, na.rm = T),
            duration = mean(duration, na.rm = T),
            icum = mean(icum, na.rm = T), .groups = "drop") %>% 
  mutate(Ecoregion = "Mediterranean")

# Combine and order factor for plotting
eco_MME_MHW_JJASON_all <- rbind(eco_MME_MHW_JJASON, eco_MME_MHW_med_JJASON) %>% 
  mutate(Ecoregion = factor(Ecoregion, 
                            levels = c("Mediterranean",
                                       "Alboran Sea", "Northwestern Mediterranean", 
                                       "Southwestern Mediterranean", "Adriatic Sea",
                                       "Ionian Sea", "Tunisian Plateau/Gulf of Sidra",
                                       "Aegean Sea", "Levantine Sea")))

# Create figures
fig_2_JJASON <- bar_dur_fig(eco_MME_MHW_JJASON_all, " (JJASON)")

# Save
ggsave("figures/fig_2.png", fig_2_JJASON, height = 8, width = 8)


# Figure 3 ----------------------------------------------------------------

# Scatterplots showing iCum and MME % damage
# Use MME pixels

# Load MME MHW pairing data
site_MME_MHW_summary <- read_csv("data/site_MME_MHW_summary.csv")

# Load species grouping sheet
species_groups <- read_csv("data/MME_MHWs_relationship_species_selection.csv") %>% 
  `colnames<-`(c("species", "damage", "group", "group_single"))

# Join MME to MHW
mme_mhw <- mme %>% 
  left_join(site_MME_MHW_summary, by = c("lon" = "lon_mme", "lat" = "lat_mme",
                                         "year", "Ecoregion", "Location", 
                                         "Monitoring series", "EvenStart", "Damaged qualitative")) %>% 
  dplyr::rename(`Damaged percentage` = `Damaged percentage.x`,
                `Damaged percentage (mean)` = `Damaged percentage.y`)

# Extract only records with regular monitoring
mme_reg <- mme_mhw %>% 
  filter(`Monitoring series` %in% c("more.than.two.per.year", "one.per.year.monitoring") | Ecoregion == "Alboran Sea") %>% 
  group_by(year, Ecoregion, lon_sst, lat_sst, Taxa, Species) %>% summarise_all(mean, na.rm = T) 
  #filter(selected_5 %in% c("2015_MHW", "2016_MHW", "2017_MHW", "2018_MHW", "2019_MHW"))

# Create data.frames based on four pre-determined filter columns
mme_Plot_1A <- filter(mme_mhw, Plot_1A %in% c("2015_MHW", "2016_MHW", "2017_MHW", "2018_MHW", "2019_MHW")) %>% 
  group_by(year, Ecoregion, lon_sst, lat_sst, Taxa, Species) %>% summarise_all(mean, na.rm = T)
mme_Plot_1B <- filter(mme_mhw, Plot_1B %in% c("2015_MHW", "2016_MHW", "2017_MHW", "2018_MHW", "2019_MHW")) %>% 
  group_by(year, Ecoregion, lon_sst, lat_sst, Taxa, Species) %>% summarise_all(mean, na.rm = T)
mme_Plot_2A <- filter(mme_mhw, Plot_2A %in% c("2015_MHW", "2016_MHW", "2017_MHW", "2018_MHW", "2019_MHW")) %>% 
  group_by(year, Ecoregion, lon_sst, lat_sst, Taxa, Species) %>% summarise_all(mean, na.rm = T)
mme_Plot_2B <- filter(mme_mhw, Plot_2B %in% c("2015_MHW", "2016_MHW", "2017_MHW", "2018_MHW", "2019_MHW")) %>% 
  group_by(year, Ecoregion, lon_sst, lat_sst, Taxa, Species) %>% summarise_all(mean, na.rm = T)
mme_Plot_3A <- filter(mme_mhw, Plot_3A %in% c("2015_MHW", "2016_MHW", "2017_MHW", "2018_MHW", "2019_MHW")) %>% 
  group_by(year, Ecoregion, lon_sst, lat_sst, Taxa, Species) %>% summarise_all(mean, na.rm = T)
mme_Plot_3B <- filter(mme_mhw, Plot_3B %in% c("2015_MHW", "2016_MHW", "2017_MHW", "2018_MHW", "2019_MHW")) %>% 
  group_by(year, Ecoregion, lon_sst, lat_sst, Taxa, Species) %>% summarise_all(mean, na.rm = T)

# List of grouped species
spp_1 <- filter(species_groups, group == "1")
spp_2 <- filter(species_groups, group == "2")
spp_3 <- filter(species_groups, group == "3")

# Filter out regularly monitored sites by species groups
mme_reg_1 <-  filter(mme_reg, Species %in% spp_1$species)
mme_reg_2 <-  filter(mme_reg, Species %in% spp_2$species)
mme_reg_3 <-  filter(mme_reg, Species %in% spp_3$species)

# Paramuricea clavata is perhaps the best candidate for heat stress
mme_reg_Pcla <- filter(mme_reg, Species == "Paramuricea clavata",
                       Ecoregion != "Ionian Sea") # Only two data points, which causes a correlation error

# Scatterplot of Paramuricea clavata MME vs MHW per pixel
scatter_Pcla <- species_scatter(mme_reg_Pcla, "Paramuricea clavata ")
ggsave("figures/scatter_Pcla.png", scatter_Pcla, height = 6, width = 8)

# Scatterplot for group 1 species
scatter_spp_1 <- species_scatter(mme_reg_1, "Group 1 species ")
ggsave("figures/scatter_spp_1.png", scatter_spp_1, height = 6, width = 8)

# Scatterplot for group 1+2 species
scatter_spp_1_2 <- species_scatter(rbind(mme_reg_1, mme_reg_2), "Group 1+2 species ")
ggsave("figures/scatter_spp_1_2.png", scatter_spp_1_2, height = 6, width = 8)

# Scatterplot for group 1+2 species
scatter_spp_1_2_3 <- species_scatter(rbind(mme_reg_1, mme_reg_2, mme_reg_3), "Group 1+2+3 species ")
ggsave("figures/scatter_spp_1_2_3.png", scatter_spp_1_2_3, height = 6, width = 8)

# Scatterplot for all species
scatter_spp_all <- species_scatter(mme_reg, "All species ")
ggsave("figures/fig_3.png", scatter_spp_all, height = 6, width = 8)

# Plots for the four pre-determined columns
scatter_Plot_1A <- species_scatter(mme_Plot_1A, "Plot 1A ")
ggsave("figures/Plot_1A.png", scatter_Plot_1A, height = 6, width = 8)
write_csv(mme_Plot_1A, "data/MME_MHW_Plot_1A.csv")
scatter_Plot_1B <- species_scatter(mme_Plot_1B, "Plot 1B ")
ggsave("figures/Plot_1B.png", scatter_Plot_1B, height = 6, width = 8)
write_csv(mme_Plot_1B, "data/MME_MHW_Plot_1B.csv")
scatter_Plot_2A <- species_scatter(mme_Plot_2A, "Plot 2A ")
ggsave("figures/Plot_2A.png", scatter_Plot_2A, height = 6, width = 8)
write_csv(mme_Plot_2A, "data/MME_MHW_Plot_2A.csv")
scatter_Plot_2B <- species_scatter(mme_Plot_2B, "Plot 2B ")
ggsave("figures/Plot_2B.png", scatter_Plot_2B, height = 6, width = 8)
write_csv(mme_Plot_2B, "data/MME_MHW_Plot_2B.csv")
scatter_Plot_3A <- species_scatter(mme_Plot_3A, "Plot 3A ")
ggsave("figures/Plot_3A.png", scatter_Plot_3A, height = 6, width = 8)
write_csv(mme_Plot_3A, "data/MME_MHW_Plot_3A.csv")
scatter_Plot_3B <- species_scatter(mme_Plot_3B, "Plot 3B ")
ggsave("figures/Plot_3B.png", scatter_Plot_3B, height = 6, width = 8)
write_csv(mme_Plot_3B, "data/MME_MHW_Plot_3B.csv")

# Correlation results
mme_reg %>% 
  group_by(Ecoregion) %>% 
  na.omit() %>% 
  summarise(r_dur = cor.test(`Damaged percentage`, duration)$estimate,
            p_dur = cor.test(`Damaged percentage`, duration)$p.value,
            r_icum = cor.test(`Damaged percentage`, icum)$estimate,
            p_icum = cor.test(`Damaged percentage`, icum)$p.value)

# Save NW and SW Med data for Quim
write_csv(mme_reg, "data/MME_MHW.csv")
mme_reg_sub <- mme_reg %>% 
  filter(Ecoregion %in% c("Northwestern Mediterranean", "Souththwestern Mediterranean"))
write_csv(mme_reg_sub, "data/MME_NW_SW.csv")


# Figure 4 ----------------------------------------------------------------

# Barplot showing relationship between MHW days and percent damage per ecoregion
## Bar colour would show one of each of the five years
##  Floating symbols would show percent mortality


# MHW stats taken from the average for all ecoregion pixels
# load("data/MHW_cat_region.RData")

# Also create a version using the co-occurring pixels
# load("data/MHW_cat_region_pixel.RData")

# Use only the species from AX column: YES
## Only take records that showed mortality
mme_selected_2_mort <- filter(mme, selected_2 == "YES") %>% #, `Damaged qualitative` != "No") %>% 
  # group_by(year, Ecoregion, lon, lat) %>% 
  # summarise(`Damaged percentage` = mean(`Damaged percentage`), .groups = "drop")
  # summarise()
  mutate(mme_damage = case_when(`Damaged qualitative` == "No" ~ 0, TRUE ~ 1)) %>% 
  # group_by(lon, lat, year) %>% 
  group_by(Ecoregion, year) %>% 
  mutate(mme_record = n(),
         mme_prop = sum(mme_damage)/mme_record) %>% 
  # dplyr::select(lon, lat, year, mme_prop, mme_record) %>% 
  dplyr::select(Ecoregion, year, mme_prop, mme_record) %>% 
  distinct()

# Match nearest pixels
mme_mhw_pixel_match <- grid_match(mme_selected_2_mort_pixel[c("lon", "lat")],
                                  MHW_pixels[c("lon", "lat")]) %>% 
  dplyr::rename(lon_mme = lon.x, lat_mme = lat.x, lon = lon.y, lat = lat.y) %>% 
  distinct() %>% 
  left_join(mme_selected_2_mort[,c("lon", "lat", "Ecoregion")], by = c("lon_mme" = "lon", "lat_mme" = "lat"))

# Get MHW file subset
file_sub <- data.frame(lat_index = seq_len(length(unique(med_sea_coords$lat))),
                       lat = unique(med_sea_coords$lat)) %>% 
  filter(lat %in% mme_mhw_pixel_match$lat)

# Calculate broad MHW stats per region/matching pixels
registerDoParallel(cores = 14)
system.time(
MHW_cat_region <- plyr::ldply(unique(med_regions$Ecoregion), region_calc, .parallel = F,
                                    mme_select = mme_selected_2_mort, pixel_sub = "full")
) # 7 minutes on 14 cores
system.time(
MHW_cat_region_pixel <- plyr::ldply(unique(med_regions$Ecoregion), region_calc, .parallel = F,
                                    mme_select = mme_selected_2_mort, pixel_sub = "pixel")
) # 68 seconds on 14 cores

# MHW stats for JJASON period
MHW_cat_region_JJASON <- MHW_cat_region %>% 
  filter(year >= 2015, month %in% c("juin", "juil", "août", "sept", "oct", "nov")) %>% 
  mutate(duration = duration/pixels) %>% 
  group_by(region, year) %>% 
  summarise(duration = sum(duration), .groups = "drop")
MHW_cat_region_pixel_JJASON <- MHW_cat_region_pixel %>% 
  filter(year >= 2015, month %in% c("juin", "juil", "août", "sept", "oct", "nov")) %>% 
  mutate(duration = duration/pixels) %>% 
  group_by(region, year) %>% 
  summarise(duration = sum(duration), .groups = "drop")

# Create region averages by all pixels
mme_region <- mme_selected_2_mort %>% 
  group_by(year, Ecoregion) %>%
  mutate(mme_prop = mme_prop*100) %>% 
  filter(Ecoregion == "Northwestern Mediterranean")
  # summarise(`Damaged percentage` = mean(`Damaged percentage`, na.rm = T), .groups = "drop")
mme_mhw_region <- left_join(mme_region, MHW_cat_region_JJASON, by = c("Ecoregion" = "region", "year"))
mme_mhw_region_pixel <- left_join(mme_region, MHW_cat_region_pixel_JJASON, by = c("Ecoregion" = "region", "year"))

# Create plot
fig_4_region_JJASON <- bar_dur_fig(mme_mhw_region, " for full ecoregion (JJASON)")
ggsave("figures/fig_4_region_JJASON.png", fig_4_region_JJASON, height = 6, width = 10)
fig_4_region_pixel_JJASON <- bar_dur_fig(mme_mhw_region_pixel, " for matching pixels (JJASON)")
ggsave("figures/fig_4_region_pixel_JJASON.png", fig_4_region_pixel_JJASON, height = 6, width = 10)


# Figure 5 ----------------------------------------------------------------

# Analysis of mortality per pixel per ecoregion
## Get MHW days and relate them to percentage of observations that have mortality
load("data/MHW_cat_pixel_annual_JJASON.RData")

## Use filter 3B for these analyses
# mme_3B <- filter(mme, Plot_3B %in% c("2015_MHW", "2016_MHW", "2017_MHW", "2018_MHW", "2019_MHW")) %>% 
mme_3B <- mme %>%
  filter(Ecoregion == "Northwestern Mediterranean") %>%
  mutate(mme_damage = case_when(`Damaged qualitative` == "No" ~ 0, TRUE ~ 1)) %>% 
  group_by(Ecoregion, `Area Monitored`, lon, lat, year) %>%
  # group_by(Ecoregion, year) %>% 
  mutate(mme_record = n(),
         mme_prop = sum(mme_damage)/mme_record) %>% 
  ungroup() %>% 
  dplyr::select(Ecoregion, `Area Monitored`, lon, lat, year, mme_prop, mme_record) %>%
  # dplyr::select(Ecoregion, year, mme_prop, mme_record) %>% 
  distinct()

# Match nearest pixels
# mme_mhw_pixel_match <- grid_match(mme_3B[c("lon", "lat")],
mme_mhw_pixel_match <- grid_match(mme[c("lon", "lat")],
                                  MHW_pixels[c("lon", "lat")]) %>% 
  dplyr::rename(lon = lon.x, lat = lat.x, lon_sst = lon.y, lat_sst = lat.y) %>% 
  distinct()
mme_mhw_3B <- mme_3B %>% 
  left_join(mme_mhw_pixel_match, by = c("lon", "lat")) %>% 
  left_join(MHW_cat_pixel_annual_JJASON, by = c("year", "lon_sst" = "lon", "lat_sst" = "lat")) %>% 
  group_by(Ecoregion, `Area Monitored`) %>% 
  summarise(duration_sum = mean(duration_sum, na.rm = T),
            mme_prop = mean(mme_prop, na.rm = T), .groups = "drop")

# Perform analysis against "severe" MME
## Or all MME above the lowest category etc.
## Create a scatterplot figure showing these results, similar to ManuFig4
## This may not work for all ecoregions
mme_mhw_3B_label_all <- mme_mhw_3B %>%
  # filter(`Damaged qualitative` != "No") %>% 
  summarise(count = n(),
            r_val = round(cor.test(mme_prop, duration_sum)$estimate, 2),
            p_val = round(cor.test(mme_prop, duration_sum)$p.value, 2),
            x_point = sum(range(duration_sum, na.rm = T))/2, .groups = "drop") %>% 
  mutate(p_val = case_when(p_val < 0.01 ~ "p < 0.01",
                           TRUE ~ paste0("p = ",p_val)))
fig_5_all <- mme_mhw_3B %>% 
  # filter(`Damaged qualitative` != "No") %>% 
  ggplot(aes(x = duration_sum, y = mme_prop)) +
  geom_smooth(method = "lm", se = F, colour = "black") +
  geom_point() +
  geom_label(data = mme_mhw_3B_label_all, alpha = 0.6, 
             aes(y = 0.6, x = x_point, label = paste0("r = ",r_val,", ",p_val, ", n = ",count))) +
  scale_y_continuous(limits = c(-0.1, 1.1), breaks = c(0.25, 0.50, 0.75, 0.1)) +
  # scale_x_continuous(limits = c(0, 125), breaks = c(30, 60, 90, 120)) +
  coord_cartesian(expand = F) +
  guides(colour = guide_legend(override.aes = list(shape = 15, size = 5))) +
  # scale_colour_brewer(palette = "Set1") +
  labs(y = "Pixel proportion moderate+", x = "MHW days") +
  theme(legend.position = "bottom")
fig_5_all
ggsave("figures/fig_5_all.png", fig_5_all, height = 6, width = 7)

# Per ecoregion
mme_mhw_3B_label_ecoregion <- mme_mhw_3B %>%
  group_by(Ecoregion) %>% 
  summarise(count = n(),
            r_val = round(cor.test(prop_mod_up, duration_sum)$estimate, 2),
            p_val = round(cor.test(prop_mod_up, duration_sum)$p.value, 2),
            x_point = sum(range(duration_sum, na.rm = T))/2, .groups = "drop") %>% 
  mutate(p_val = case_when(p_val < 0.01 ~ "p < 0.01",
                           TRUE ~ paste0("p = ",p_val)))
fig_5_ecoregion <- mme_mhw_3B %>%
  # filter(`Damaged qualitative` != "No") %>%
  ggplot(aes(x = duration_sum, y = mme_prop)) +
  geom_smooth(method = "lm", se = F, colour = "black") +
  geom_point() +
  geom_label(data = mme_mhw_3B_label_ecoregion, alpha = 0.6, 
             aes(y = 0.6, x = x_point, label = paste0("r = ",r_val,", ",p_val, ", n = ",count))) +
  scale_y_continuous(limits = c(0.4, 1.1), breaks = c(0.50, 0.75, 0.1)) +
  scale_x_continuous(limits = c(0, 125), breaks = c(30, 60, 90, 120)) +
  facet_wrap(~Ecoregion) +
  coord_cartesian(expand = F) +
  guides(colour = guide_legend(override.aes = list(shape = 15, size = 5))) +
  # scale_colour_brewer(palette = "Set1") +
  labs(y = "Pixel proportion moderate+", colour = "Taxa", x = "MHW days") +
  theme(legend.position = c(0.8, 0.2))
fig_5_ecoregion
ggsave("figures/fig_5_ecoregion.png", fig_5_ecoregion, height = 10, width = 15)


# Manuscript figure 1 -----------------------------------------------------

## A: Map of temperature difference mean 1982-1986 vs 2015-2019
# TODO: Add ecoregions to all maps
# Load data
load("data/MHW_clim_pixel_annual.RData")
# Prep data
pixel_pentad <- MHW_clim_pixel_annual %>% 
  right_join(med_regions, by = c("lon", "lat")) %>% 
  mutate(pentad = cut(year, c(1981, 1986, 1992, 1998, 2003, 2009, 2014, 2019))) %>% 
  group_by(lon, lat, pentad) %>% 
  summarise(temp = mean(temp, na.rm = T), .groups = "drop") %>% 
  pivot_wider(id_cols = c("lon", "lat"), names_from = "pentad", values_from = "temp") %>% 
  mutate(temp_diff = `(2014,2019]` - `(1981,1986]`)
# Plot data
panel_A <- ggplot(pixel_pentad, aes(x = lon, y = lat)) +
  geom_tile(aes(fill = temp_diff), colour = NA) +
  geom_polygon(data = map_base, aes(x = lon, y = lat, group = group)) +
  # scale_fill_manual("Category", values = MHW_colours) +
  scale_fill_gradient2(low = "yellow", mid = "orange", high = "red",
                       breaks = c(0.9, 1.2, 1.5), midpoint = 1.1) +#,
                       # labels = c("0", "0.5", "1.0", "1.5", "2.0")) +
                       # labels = c("0", "0.5", "1.0", "1.5", "2.0")) +
  coord_cartesian(expand = F, 
                  xlim = c(min(med_regions$lon), max(med_regions$lon)),
                  ylim = c(min(med_regions$lat), max(med_regions$lat))) +
  # theme_void() +
  labs(title = "__(a)__   Temperature difference [2015 to 2019] minus [1982 to 1986]",
       y = "Latitude (°N)", x = "Longitude (°E)", fill = "Temp. (°C)") +
  # guides(fill = guide_legend(override.aes = list(size = 10))) +
  theme(plot.title = ggtext::element_markdown(),
        panel.border = element_rect(colour = "black", fill = NA),
        # legend.position = "bottom",
        legend.position = c(0.9, 0.79),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        panel.background = element_rect(fill = "grey90"))
# panel_A

## B: Barplots of mean Med SST anom with ~5 year average segments
# Prep data
med_annual <- MHW_clim_pixel_annual %>% 
  right_join(med_regions, by = c("lon", "lat")) %>% 
  group_by(year) %>% 
  summarise(temp = mean(temp, na.rm = T), .groups = "drop") %>% 
  mutate(anom = temp - mean(temp))
med_pentad <- med_annual %>% 
  mutate(pentad = cut(year, c(1981, 1986, 1992, 1998, 2003, 2009, 2014, 2019))) %>% 
  group_by(pentad) %>% 
  summarise(anom_pentad = mean(anom, na.rm = T), .groups = "drop") %>%
  separate(pentad, into = c("start_year", "end_year"), sep = ",", remove = F) %>% 
  mutate(start_year = as.numeric(sub("[(]", "", start_year)) + 1,
         end_year = as.numeric(sub("]", "", end_year)))
# Plot data
panel_B <- ggplot(med_annual, aes(x = year, y = anom)) +
  geom_bar(aes(fill = anom), stat = "identity", width = 1, show.legend = F) +
  geom_hline(yintercept = 0) +
  geom_segment(data = med_pentad, size = 2, lineend = "round",
               aes(x = start_year, xend = end_year, 
                   y = anom_pentad, yend = anom_pentad)) +
  scale_fill_gradient2(low = "#4575b4", high = "#d73027") +
  scale_x_continuous(breaks = seq(1984, 2019, 7), expand = c(0, 0)) +
  labs(title = "__(b)__   Annual SST anomalies [1982 to 2019]",
       y = "Temperature (°C)", x = "Year") +
  theme(plot.title = ggtext::element_markdown(),
        panel.border = element_rect(colour = "black", fill = NA),
        # legend.position = "bottom",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        panel.background = element_rect(fill = "grey90"))
# panel_B

## C: Map of difference in cat 2+ days between first and last pentad
# Custom legend
MHW_colours_compare <- c(
  "Same" = "salmon",
  MHW_colours[2],
  MHW_colours[3],
  MHW_colours[4]
)

# Prep data
pixel_cat_pentad <- MHW_cat_pixel_annual %>% 
  right_join(med_regions, by = c("lon", "lat")) %>%
  group_by(lon, lat, year) %>% 
  summarise(cat2 = sum(`II Strong`, `III Severe`, `IV Extreme`),
            cat_max = max(as.numeric(category), na.rm = T), .groups = "drop") %>% 
  mutate(pentad = cut(year, c(1981, 1986, 1992, 1998, 2003, 2009, 2014, 2019))) %>% 
  group_by(lon, lat, pentad) %>% 
  summarise(cat2 = mean(cat2, na.rm = T), 
            cat_max = max(cat_max, na.rm = T), .groups = "drop") %>% 
  dplyr::select(-cat2) %>% 
  pivot_wider(id_cols = c("lon", "lat"), names_from = "pentad", values_from = "cat_max") %>%
  # mutate(`(1981,1986]` = replace_na(`(1981,1986]`, 0)) %>%
  replace(is.na(.), 0) %>% 
  mutate(cat_max_diff = case_when(`(1981,1986]` < `(2014,2019]` ~ "greater",
                                  `(1981,1986]` == `(2014,2019]` ~ "same",
                                  `(1981,1986]` > `(2014,2019]` ~ "less")) %>% 
  mutate(cat_max_diff_cat = ifelse(`(1981,1986]` <= 1 & `(1981,1986]` < `(2014,2019]`, `(2014,2019]`, NA)) %>% 
  mutate(cat_max_diff_cat = case_when(cat_max_diff_cat == 2 ~ "II Strong",
                                      cat_max_diff_cat == 3 ~ "III Severe",
                                      cat_max_diff_cat == 4 ~ "IV Extreme",
                                      TRUE ~ "Same"))

# Plot data
panel_C <- pixel_cat_pentad %>%
  na.omit() %>%
  # group_by(lon, lat) %>%
  # summarise(category = max(as.numeric(category), na.rm = T), .groups = "drop") %>%
  # mutate(cat_max_diff_cat = factor(cat_max_diff_cat, labels = c("I Moderate", "II Strong", "III Severe", "IV Extreme"))) %>%
  ggplot() +
  # geom_tile(data = OISST_ice_coords, fill = "powderblue", colour = NA, alpha = 0.5) +
  geom_tile(data = filter(pixel_cat_pentad, cat_max_diff_cat == "Same"), 
            aes(x = lon, y = lat, fill = cat_max_diff_cat), alpha = 0.3) + # This background fill colour needs tweaking
  geom_tile(data = filter(pixel_cat_pentad, cat_max_diff_cat != "Same"),
            aes(x = lon, y = lat, fill = cat_max_diff_cat), colour = NA, show.legend = T) +
  geom_polygon(data = map_base, aes(x = lon, y = lat, group = group),
               fill = "grey70", colour = "black") +
  geom_sf(data = MEOW, alpha = 1, aes(geometry = geometry), fill = NA, colour = "forestgreen") +
  # scale_fill_manual(values = c("red", "white", "blue")) +
  scale_fill_manual("Category", values = MHW_colours_compare,
                    breaks = c("Same", "II Strong", "III Severe", "IV Extreme")) +
  # scale_fill_gradient2(low = "darkorchid", mid = "white", high = "hotpink", midpoint = 0) +#,
                       # breaks = c(0.9, 1.2, 1.5), midpoint = 1.1) +
  # scale_y_continuous(breaks = NULL) +
  # scale_x_continuous(breaks = NULL) +
  coord_sf(expand = F,
           xlim = c(min(med_regions$lon), max(med_regions$lon)),
           ylim = c(min(med_regions$lat), max(med_regions$lat))) +
  # theme_void() +
  # theme_bw() +
  labs(title = "__(c)__   The highest category for [2015 to 2019] when [1982 to 1986] was 'I Moderate' or less",
       y = "Latitude (°N)", x = "Longitude (°E)", fill = "Max category") +
  # guides(fill = guide_legend(override.aes = list(size = 10))) +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        plot.title = ggtext::element_markdown(),
        legend.position = c(0.89, 0.78),
        # legend.position = "bottom",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        panel.background = element_rect(fill = "grey90"))
# panel_C

## D: Barplot of Med surface area affected by Cat 2+ MHWs 
# Shortened colour palette
MHW_colours_no_mod <- c(
  MHW_colours[2],
  MHW_colours[3],
  MHW_colours[4]
)

# Load data
load("data/MHW_cat_summary_annual.RData")
OISST_global <- readRDS("data/OISST_cat_daily_1992-2018_total.Rds") %>% 
  filter(category != "I Moderate") %>% 
  group_by(t) %>% 
  mutate(cat_n_prop_stack = cumsum(cat_n_prop),
         first_n_cum_prop_stack = cumsum(first_n_cum_prop)) %>% 
  filter(category == "IV Extreme")
# Prep data
cat_daily_mean <- MHW_cat_summary_annual %>%
  filter(category != "I Moderate") %>% 
  group_by(year, category) %>%
  summarise(first_n_cum_prop = max(first_n_cum_prop),
            cat_n_prop_mean = mean(cat_n_prop, na.rm = T),
            cat_n_cum_prop = max(cat_n_cum_prop, na.rm = T), .groups = "drop")
cat_pentad <- cat_daily_mean %>% 
  group_by(year) %>% 
  summarise(cat_n_cum_prop_sum = sum(cat_n_cum_prop, na.rm = T),
            first_n_cum_prop_sum = sum(first_n_cum_prop), .groups = "drop") %>% 
  mutate(pentad = cut(year, c(1981, 1986, 1992, 1998, 2003, 2009, 2014, 2019))) %>% 
  group_by(pentad) %>% 
  summarise(cat_n_cum_prop_pentad = mean(cat_n_cum_prop_sum, na.rm = T),
            first_n_cum_prop_pentad = mean(first_n_cum_prop_sum), .groups = "drop") %>%
  separate(pentad, into = c("start_year", "end_year"), sep = ",", remove = F) %>% 
  mutate(start_year = as.numeric(sub("[(]", "", start_year)) + 1,
         end_year = as.numeric(sub("]", "", end_year)))
# Plot data
panel_D <- ggplot(cat_daily_mean, aes(x = year, y = first_n_cum_prop)) +
  geom_bar(aes(fill = category), stat = "identity", show.legend = T,
           position = position_stack(reverse = TRUE), width = 1) +
  geom_segment(data = cat_pentad, size = 2, lineend = "round",
               aes(x = start_year, xend = end_year, 
                   y = first_n_cum_prop_pentad, yend = first_n_cum_prop_pentad)) +
  geom_point(data = OISST_global, aes(x = t, y = first_n_cum_prop_stack), 
             shape = 21, fill = "grey", show.legend = F) +
  scale_fill_manual("Category", values = MHW_colours_no_mod) +
  scale_colour_manual("Category", values = MHW_colours_no_mod) +
  # scale_y_continuous(limits = c(0, 20),
  #                    breaks = seq(5, 15, length.out = 3)) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = c(0.25, 0.5, 0.75),
                     labels = c("25", "50", "75")) +
  scale_x_continuous(breaks = seq(1984, 2019, 7)) +
  guides(pattern_colour = FALSE, colour = FALSE) +
  # labs(title = "__(d)__   Surface area affected by category 'II Strong'+ MHWs",
  labs(title = "Surface area affected by category 'II Strong'+ MHWs",
       y = "Cover (%)", x = "Year") +
  coord_cartesian(expand = F) +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        # plot.title = ggtext::element_markdown(),
        legend.position = c(0.1, 0.85),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))
# panel_D

## Combine and save
manu_fig_1 <- ggpubr::ggarrange(panel_A, panel_B, panel_C, panel_D, labels = c("(a)", "(b)", "(c)", "(d)"))
ggsave("figures/manu_fig_1.png", manu_fig_1, height = 10, width = 16)

# Summary stats for text
cat_daily_mean %>% 
  group_by(year) %>% 
  summarise(first_n_cum_prop_sum = sum(first_n_cum_prop)) %>% 
  arrange(-first_n_cum_prop_sum)

# Difference in SST from first to last pentad
med_pentad$anom_pentad[7]-med_pentad$anom_pentad[1]


# Manuscript figure 4 -----------------------------------------------------

# Scatterplot
# Keep only Med panel
# Show taxa as the colour of dots
# Use all data points. No averaging per pixel.
# Create multiple panels: One with no averaging, and others with layers of averaging

# TODO: Look into what relationship may exist between the vertical movement of the dots per year
# E.g. the range of MHW days in a year is not large, but MME damage is
# What is it that drives the dots up within a narrow x-axis range? Species? Location?

# Load MME MHW pairing data
site_MME_MHW_summary <- read_csv("data/site_MME_MHW_summary.csv")

# Load species grouping sheet
species_groups <- read_csv("data/MME_MHWs_relationship_species_selection.csv") %>% 
  `colnames<-`(c("species", "damage", "group", "group_single"))

# Join with MHW data to get pixel lon/lat
mme_mhw <- mme %>% 
  left_join(site_MME_MHW_summary, by = c("lon" = "lon_mme", "lat" = "lat_mme",
                                         "year", "Ecoregion", "Location", 
                                         "Monitoring series", "EvenStart", "Damaged qualitative")) %>% 
  dplyr::rename(`Damaged percentage` = `Damaged percentage.x`,
                `Damaged percentage (mean)` = `Damaged percentage.y`)

# Extract only records with regular monitoring
mme_reg <- mme_mhw %>% 
  filter(`Monitoring series` %in% c("more.than.two.per.year", "one.per.year.monitoring"))# | Ecoregion == "Alboran Sea") #%>% 
  # group_by(year, Ecoregion, lon_sst, lat_sst, Tax) %>% summarise_all(mean, na.rm = T) 

# The 3B filtered group
mme_3B <- mme_mhw %>% 
  filter(Plot_3B %in% c("2015_MHW", "2016_MHW", "2017_MHW", "2018_MHW", "2019_MHW"))

# Final filtering of 3B filter
mme_final <- mme_3B %>% 
  filter(Taxa != "Mollusca (Bivalvia)", Taxa != "Tracheophyta", Species != "Encrusting calcareous algae")

# Average pixel distance
mean(mme_final$dist, na.rm = T)
range(mme_final$dist, na.rm = T)

# MHW to MME occurrence stats
nrow(filter(mme_final, count_MHW == 0))/nrow(mme_final)

# Function for showing scatterplots for full med with different rounding methods
# df <- mme_final; round_type <- "species"; x_var <- "mhw_days" # Look inside the function
species_scatter_full <- function(df, round_type, x_var){
  # Round data accordingly
  if(round_type == "none"){
    df_round <- df
    plot_sub <- "Filter 3B; full records"
    write_csv(df_round, "data/MME_MHW_Plot_3B.csv")
  } else if(round_type == "species"){
    df_round <- df %>% 
      group_by(lon_sst, lat_sst, year, Taxa, Species) %>% 
      summarise_all(mean, na.rm  = T, .groups = "drop")
    plot_sub <- "Filter 3B; average per species per pixel per year"
    write_csv(df_round, "data/MME_MHW_Plot_3B_species.csv")
  } else if(round_type == "taxa"){
    df_round <- df %>% 
      group_by(lon_sst, lat_sst, year, Taxa) %>% 
      summarise_all(mean, na.rm  = T, .groups = "drop")
    plot_sub <- "Filter 3B; average per taxa per pixel per year"
    write_csv(df_round, "data/MME_MHW_Plot_3B_taxa.csv")
  } else if(round_type == "pixel"){
    df_round <- df %>% 
      group_by(lon_sst, lat_sst, year) %>% 
      summarise_all(mean, na.rm  = T, .groups = "drop") %>% 
      mutate(Taxa = "pixel")
    plot_sub <- "Filter 3B; average per pixel per year"
    write_csv(df_round, "data/MME_MHW_Plot_3B_pixel.csv")
  }
  # Set x-axis label
  if(x_var == "mhw_days"){
    x_lab <- "MHW days"
    plot_title <- "MME damage vs MHW days (JJASON)"
  } else if(x_var == "e_days"){
    x_lab <- "Days above 90th percentile threshold"
    plot_title <- "MME damage vs days above 90th perc. thresh. (JJASON)"
  } else if(x_var == "sum_anom"){
    x_lab <- "Sum of temperature anomalies (°C days)"
    plot_title <- "MME damage vs temperature anomalies. (JJASON)"
  } else if (x_var == "icum"){
    x_lab <- "MHW cumulative intensity (°C days)"
    plot_title <- "MME damage vs MHW cumulative intensity JJASON)"
  }
  # Cast long for more control
  df_long <- df_round %>%
    pivot_longer(mhw_days:icum) %>%
    filter(name == x_var)  # Pick the variable for the X-axis
  # Get label for plot
  df_label <- df_long %>%
    ungroup() %>% 
    summarise(count = n(),
              r_val = round(cor.test(`Damaged percentage`, value)$estimate, 2),
              p_val = round(cor.test(`Damaged percentage`, value)$p.value, 2),
              x_point = sum(range(value, na.rm = T))/2, .groups = "drop") %>% 
    mutate(p_val = case_when(p_val < 0.01 ~ "p < 0.01",
                             TRUE ~ paste0("p = ",p_val)))
  # The figure
  ggplot(data = df_long, aes(x = value, y = `Damaged percentage`)) +
    # geom_smooth(data = df_15, linetype = "dashed", colour = "black", method = "lm", se = F) +
    geom_smooth(colour = "black", method = "lm", se = F) +
    geom_point(aes(colour = Taxa)) +
    # geom_point() + # Simple black dots
    # geom_label(data = df_label, aes(y = 10, x = x_point, label = paste0("n = ",count)), alpha = 0.6) +
    geom_label(data = df_label, alpha = 0.6, 
               aes(y = 90, x = x_point, label = paste0("r = ",r_val,", ",p_val, ", n = ",count))) +
    guides(colour = guide_legend(override.aes = list(shape = 15, size = 5))) +
    # scale_colour_brewer(palette = "Set1") +
    labs(y = "MME damage (%)", colour = "Taxa", x = x_lab,
         # subtitle = plot_sub,
         title = plot_title) +
    # facet_wrap(~Ecoregion, scales = "free_x") +#, strip.position = "bottom") +
    scale_y_continuous(limits = c(-2, max(df_round$`Damaged percentage`, na.rm = T)+2)) +
    # scale_x_continuous(limits = c(-2, max(df$duration, na.rm = T)*1.1)) +
    theme(panel.border = element_rect(colour = "black", fill = NA),
          legend.position = "bottom", legend.box = "vertical",
          strip.placement = "outside", strip.background = element_blank())
}

## Different rounding approaches
panel_all_days <- species_scatter_full(mme_final, "none", "mhw_days")
panel_species_days <- species_scatter_full(mme_final, "species", "mhw_days")
panel_taxa_days <- species_scatter_full(mme_final, "taxa", "mhw_days")
panel_pixel_days <- species_scatter_full(mme_final, "pixel", "mhw_days")
panel_all_icum <- species_scatter_full(mme_final, "none", "icum")
panel_species_icum <- species_scatter_full(mme_final, "species", "icum")
panel_taxa_icum <- species_scatter_full(mme_final, "taxa", "icum")
panel_pixel_icum <- species_scatter_full(mme_final, "pixel", "icum")

## Combine and save
# manu_fig_4_days <- ggpubr::ggarrange(panel_all_days, panel_species_days, panel_taxa_days, panel_pixel_days, nrow = 1, align = "hv")
# manu_fig_4_icum <- ggpubr::ggarrange(panel_all_icum, panel_species_icum, panel_taxa_icum, panel_pixel_icum, nrow = 1, align = "hv")
# manu_fig_4 <- ggpubr::ggarrange(manu_fig_4_days, manu_fig_4_icum, ncol = 1)
# ggsave("figures/manu_fig_4.png", manu_fig_4, height = 10, width = 25)
manu_fig_4 <- panel_species_days
ggsave("figures/manu_fig_4.png", manu_fig_4, height = 5, width = 6)

