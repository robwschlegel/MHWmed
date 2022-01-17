# code/figures.R
# This script creates annual and total summaries from the MHW Med results


# Setup -------------------------------------------------------------------

# The project-wide functions
source("code/functions.R")
library(grid)
library(gridExtra)
library(gtable)


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
# mme_mhw_pixel_match <- grid_match(mme_selected_2_mort_pixel[c("lon", "lat")],
#                                   MHW_pixels[c("lon", "lat")]) %>% 
#   dplyr::rename(lon_mme = lon.x, lat_mme = lat.x, lon = lon.y, lat = lat.y) %>% 
#   distinct() %>% 
#   left_join(mme_selected_2_mort[,c("lon", "lat", "Ecoregion")], by = c("lon_mme" = "lon", "lat_mme" = "lat"))

# Get MHW file subset
# file_sub <- data.frame(lat_index = seq_len(length(unique(med_sea_coords$lat))),
#                        lat = unique(med_sea_coords$lat)) %>% 
#   filter(lat %in% mme_mhw_pixel_match$lat)

# Calculate broad MHW stats per region/matching pixels
# registerDoParallel(cores = 14)
# system.time(
# MHW_cat_region <- plyr::ldply(unique(med_regions$Ecoregion), region_calc, .parallel = F,
#                                     mme_select = mme_selected_2_mort, pixel_sub = "full")
# ) # 7 minutes on 14 cores
# system.time(
# MHW_cat_region_pixel <- plyr::ldply(unique(med_regions$Ecoregion), region_calc, .parallel = F,
#                                     mme_select = mme_selected_2_mort, pixel_sub = "pixel")
# ) # 68 seconds on 14 cores

# MHW stats for JJASON period
load("data/MHW_cat_region.RData")
MHW_cat_region_JJASON <- MHW_cat_region %>%
  filter(year >= 2015, month %in% c("juin", "juil", "août", "sept", "oct", "nov")) %>%
  mutate(duration = duration/pixels) %>%
  group_by(region, year) %>%
  summarise(duration = sum(duration), .groups = "drop")
# MHW_cat_region_pixel_JJASON <- MHW_cat_region_pixel %>% 
#   filter(year >= 2015, month %in% c("juin", "juil", "août", "sept", "oct", "nov")) %>% 
#   mutate(duration = duration/pixels) %>% 
#   group_by(region, year) %>% 
#   summarise(duration = sum(duration), .groups = "drop")

# Create region averages by all pixels
mme_region <- mme_selected_2_mort %>% 
  group_by(year, Ecoregion) %>%
  mutate(mme_prop = (mme_prop*100)+20) %>%
  filter(Ecoregion == "Northwestern Mediterranean")
  # summarise(`Damaged percentage` = mean(`Damaged percentage`, na.rm = T), .groups = "drop")
mme_mhw_region <- left_join(mme_region, MHW_cat_region_JJASON, by = c("Ecoregion" = "region", "year"))

# Create plot
fig_4_region_JJASON <- ggplot(mme_mhw_region, aes(x = year, y = duration)) +
  geom_bar(aes(fill = as.factor(year)), 
           colour = "black",
           stat = "identity", 
           show.legend = F,
           # position = "dodge",
           width = 0.8) +
  geom_point(aes(y = mme_prop), size = 5) +
  scale_fill_viridis_d("Year", option = "D", aesthetics = c("colour", "fill")) +
  scale_y_continuous(limits = c(0, 110), breaks = c(20, 40, 60, 80, 100), expand = c(0, 0),
                     sec.axis = sec_axis(name = "MME records\n(proportion with mortality)", 
                                         trans = ~ . + 0,
                                         breaks = c(20, 45, 70, 95),
                                         labels = c("0.00", "0.25", "0.50", "0.75"))) +
  # scale_x_continuous(breaks = c(2015, 2017, 2019)) +
  # coord_cartesian(expand = F) +
  labs(x = NULL, y = "MHW days", 
       title = "Average MHW days per year in NW Med",
       subtitle = "Points show proportion of records with mortality") +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.position = "bottom",
        legend.key.width = unit(1, "cm"),
        panel.background = element_rect(fill = "grey90"), 
        strip.text = element_text(size = 12))
fig_4_region_JJASON
ggsave("figures/fig_4_region_JJASON.png", fig_4_region_JJASON, height = 6, width = 10)


# mme_mhw_region_pixel <- left_join(mme_region, MHW_cat_region_pixel_JJASON, by = c("Ecoregion" = "region", "year"))
# fig_4_region_pixel_JJASON <- bar_dur_fig(mme_mhw_region_pixel, " for matching pixels (JJASON)")
# ggsave("figures/fig_4_region_pixel_JJASON.png", fig_4_region_pixel_JJASON, height = 6, width = 10)


# Figure 5 ----------------------------------------------------------------

# Analysis of mortality per pixel per ecoregion
## Get MHW days and relate them to percentage of observations that have mortality
load("data/MHW_cat_pixel_annual_JJASON.RData")

## Use filter 3B for these analyses
mme_3B <- filter(mme, Plot_3B %in% c("2015_MHW", "2016_MHW", "2017_MHW", "2018_MHW", "2019_MHW")) %>%
# mme_3B <- filter(mme, selected_2 == "YES") %>%
# mme_3B <- mme %>% 
  # filter(Taxa != "Tracheophyta") %>%
  filter(Ecoregion == "Northwestern Mediterranean") %>%
  mutate(mme_damage = case_when(`Damaged qualitative` %in% c("No") ~ 0, TRUE ~ 1)) %>% 
  group_by(Ecoregion, `Area Monitored`, lon, lat, year) %>%
  # group_by(Ecoregion, year) %>% 
  mutate(mme_record = n(),
         mme_prop = sum(mme_damage)/mme_record,
         damage_mean = mean(`Damaged percentage`, na.rm = T)) %>% 
  ungroup() %>% 
  dplyr::select(Ecoregion, `Area Monitored`, lon, lat, year, mme_prop, mme_record, damage_mean) %>%
  # dplyr::select(Ecoregion, year, mme_prop, mme_record) %>% 
  distinct()

# Match nearest pixels
mme_mhw_pixel_match <- grid_match(mme_3B[c("lon", "lat")],
                                  MHW_pixels[c("lon", "lat")]) %>% 
  dplyr::rename(lon = lon.x, lat = lat.x, lon_sst = lon.y, lat_sst = lat.y) %>% 
  distinct()
mme_mhw_3B <- mme_3B %>% 
  left_join(mme_mhw_pixel_match, by = c("lon", "lat")) %>% 
  left_join(MHW_cat_pixel_annual_JJASON, by = c("year", "lon_sst" = "lon", "lat_sst" = "lat")) %>% 
  mutate(duration_sum = replace_na(duration_sum, 0)) %>% 
  group_by(Ecoregion, `Area Monitored`, year) %>%
  # group_by(Ecoregion, `Area Monitored`) %>% 
  summarise(duration_sum = mean(duration_sum, na.rm = T),
            mme_prop = mean(mme_prop, na.rm = T),
            damage_mean = mean(damage_mean, na.rm = T), .groups = "drop")

# Perform analysis against "severe" MME
## Or all MME above the lowest category etc.
## Create a scatterplot figure showing these results, similar to ManuFig4
## This may not work for all ecoregions
mme_mhw_3B_label_all <- mme_mhw_3B %>%
  # group_by(year) %>%
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
  geom_point(aes(fill = as.factor(year)), shape = 21, size = 5, alpha = 0.9, show.legend = F) +
  geom_label(data = mme_mhw_3B_label_all, alpha = 0.6, 
             aes(y = 0.6, x = x_point, label = paste0("r = ",r_val,", ",p_val, ", n = ",count))) +
  scale_fill_viridis_d("Year", option = "D", aesthetics = c("colour", "fill")) +
  scale_y_continuous(limits = c(-0.1, 1.1), breaks = c(0, 0.25, 0.50, 0.75, 1), 
                     labels = c("0", "0.25", "0.5", "0.75", "1.0")) +
  # scale_x_continuous(limits = c(0, 125), breaks = c(30, 60, 90, 120)) +
  coord_cartesian(expand = F) +
  # guides(colour = guide_legend(override.aes = list(shape = 15, size = 5))) +
  # scale_colour_brewer(palette = "Set1") +
  labs(y = "MME records\n(proportion with mortality)", x = "MHW days (per year; JJASON)",
       title = "MHW days and MME records per monitoring area in the NW Med",
       subtitle = "MHW days per year during JJASON period and proportion of MME records with damage") +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.position = "bottom",
        legend.key.width = unit(1, "cm"),
        panel.background = element_rect(fill = "grey90"), 
        strip.text = element_text(size = 12))
fig_5_all
ggsave("figures/fig_5_all.png", fig_5_all, height = 6, width = 7)

# Per ecoregion
# mme_mhw_3B_label_ecoregion <- mme_mhw_3B %>%
#   group_by(Ecoregion) %>% 
#   summarise(count = n(),
#             r_val = round(cor.test(mme_prop, duration_sum)$estimate, 2),
#             p_val = round(cor.test(mme_prop, duration_sum)$p.value, 2),
#             x_point = sum(range(duration_sum, na.rm = T))/2, .groups = "drop") %>% 
#   mutate(p_val = case_when(p_val < 0.01 ~ "p < 0.01",
#                            TRUE ~ paste0("p = ",p_val)))
# fig_5_ecoregion <- mme_mhw_3B %>%
#   na.omit() %>% 
#   # filter(`Damaged qualitative` != "No") %>%
#   ggplot(aes(x = duration_sum, y = mme_prop)) +
#   geom_smooth(method = "lm", se = F, colour = "black") +
#   geom_point() +
#   geom_label(data = mme_mhw_3B_label_ecoregion, alpha = 0.6, 
#              aes(y = 0.6, x = x_point, label = paste0("r = ",r_val,", ",p_val, ", n = ",count))) +
#   scale_y_continuous(limits = c(0.4, 1.1), breaks = c(0.50, 0.75, 0.1)) +
#   # scale_x_continuous(limits = c(0, 125), breaks = c(30, 60, 90, 120)) +
#   facet_wrap(~Ecoregion) +
#   coord_cartesian(expand = F) +
#   guides(colour = guide_legend(override.aes = list(shape = 15, size = 5))) +
#   # scale_colour_brewer(palette = "Set1") +
#   labs(y = "Pixel proportion moderate+", colour = "Taxa", x = "MHW days") +
#   theme(legend.position = c(0.8, 0.2))
# fig_5_ecoregion
# ggsave("figures/fig_5_ecoregion.png", fig_5_ecoregion, height = 10, width = 15)


# Figure 6 ----------------------------------------------------------------

# The trends in SST and MHW days per ecoregion
load("data/MHW_clim_pixel_annual.RData")

# Function for SST trend calculations
# df <- filter(MHW_clim_pixel_annual, 
#              lon == MHW_clim_pixel_annual$lon[1], lat == MHW_clim_pixel_annual$lat[1])
dec_trend_calc <- function(df){
  
  # Decadal trends
  dec_trend_temp <- broom::tidy(lm(temp ~ year, df)) %>% 
    slice(2) %>% 
    mutate(dec_trend_temp = round(estimate*10, 3),
           p.value_temp = round(p.value, 4)) %>% 
    dplyr::select(dec_trend_temp, p.value_temp)
  dec_trend_dur <- broom::tidy(lm(mhw_days ~ year, df)) %>% 
    slice(2) %>% 
    mutate(dec_trend_dur = round(estimate*10, 3),
           p.value_dur = round(p.value, 4)) %>% 
    dplyr::select(dec_trend_dur, p.value_dur)
  
  # Total means
  total_temp <- df %>% 
    summarise(temp_average = round(mean(temp, na.rm = T), 2), .groups = "drop")
  total_dur <- df %>% 
    summarise(dur_average = round(mean(mhw_days, na.rm = T), 1),.groups = "drop")
  
  # Combine and exit
  res <- cbind(total_temp, dec_trend_temp, total_dur, dec_trend_dur)
  rm(df, dec_trend_temp, dec_trend_dur, total_temp, total_dur); gc()
  return(res)
}

# Convenience wrapper for summarising stats
trend_summarise <- function(df){
  df %>% 
    summarise(temp_average_mean = round(mean(temp_average, na.rm = T), 2),
              temp_average_range = max(temp_average, na.rm = T) - min(temp_average, na.rm = T),
              temp_average_sd = round(sd(temp_average, na.rm = T), 2),
              dec_trend_temp_mean = round(mean(dec_trend_temp, na.rm = T), 2),
              dec_trend_temp_range = max(dec_trend_temp, na.rm = T) - min(dec_trend_temp, na.rm = T),
              dec_trend_temp_sd = round(sd(dec_trend_temp, na.rm = T), 2),
              dur_average_mean = round(mean(dur_average, na.rm = T), 2),
              dur_average_range = max(dur_average, na.rm = T) - min(temp_average, na.rm = T),
              dur_average_sd = round(sd(dur_average, na.rm = T), 2),
              dec_trend_dur_mean = round(mean(dec_trend_dur, na.rm = T), 2),
              dec_trend_dur_range = max(dec_trend_dur, na.rm = T) - min(dec_trend_temp, na.rm = T),
              dec_trend_dur_sd = round(sd(dec_trend_dur, na.rm = T), 2), .groups = "drop") %>% 
    arrange(-dec_trend_temp_mean, -dec_trend_dur_mean) %>% 
    mutate(rank = 1:n())
}

# Prep data
MHW_clim_pixel_annual_ecoregion <- MHW_clim_pixel_annual %>% 
  left_join(med_regions, by = c("lon", "lat")) %>% 
  filter(!is.na(Ecoregion))

# Load lon files
# registerDoParallel(cores = 12)
# system.time(
#   SST_MHW_trends <- plyr::ddply(MHW_clim_pixel_annual_ecoregion, c("Ecoregion", "lon", "lat"), dec_trend_calc, .parallel = T)
# ) # xxx minutes on 12 cores
# gc()
# save(SST_MHW_trends, file = "data/SST_MHW_trends.RData")
load("data/SST_MHW_trends.RData")

## Per ecoregion stats
# Average per ecoregion
SST_MHW_trends_ecoregion <- SST_MHW_trends %>% 
  group_by(Ecoregion) %>%
  trend_summarise()

# Average per province
SST_MHW_trends_province <- SST_MHW_trends %>% 
  left_join(MEOW[,c("ECOREGION", "PROVINCE")], by = c("Ecoregion" = "ECOREGION")) %>% 
  group_by(PROVINCE) %>%
  trend_summarise()

# Create table for plotting
med_trends <- bind_rows(SST_MHW_trends_ecoregion, SST_MHW_trends_province) %>% 
  dplyr::select(PROVINCE, everything()) %>% 
  mutate(PROVINCE = "Mediterranean Sea",
         Ecoregion = ifelse(is.na(Ecoregion), "ALL", Ecoregion))

# Map of SST mean per pixel
map_SST_total <- SST_MHW_trends %>%
  na.omit() %>% 
  # Rather not remove the tails of the SST distribution
  # mutate(temp_average_05 = quantile(temp_average, probs = 0.05),
  #        temp_average_95 = quantile(temp_average, probs = 0.95),
  #        temp_average = case_when(temp_average > temp_average_95 ~ temp_average_95,
  #                                 temp_average < temp_average_05 ~ temp_average_05,
  #                                 TRUE ~ temp_average)) %>% 
  ggplot() +
  geom_tile(aes(fill = temp_average, x = lon, y = lat)) +
  geom_sf(data = MEOW, aes(colour = ECOREGION), fill = NA, show.legend = F) +
  geom_polygon(data = map_base, aes(group = group, x = lon, y = lat)) +
  scale_fill_viridis_c() +
  labs(x = NULL, y = NULL, fill = "Annual average\ntemperature (°C)") +
  coord_sf(expand = F,
           xlim = c(min(med_regions$lon), max(med_regions$lon)),
           ylim = c(min(med_regions$lat), max(med_regions$lat))) +
  theme(legend.position = "top", 
        panel.background = element_rect(colour = "black", fill = NULL))
map_SST_total

# Map of decadal trend per pixel
map_SST_trend <- SST_MHW_trends %>%
  na.omit() %>% 
  filter(p.value_temp <= 0.05) %>% 
  # Rather not remove the tails
  # mutate(dec_trend_temp_05 = quantile(dec_trend_temp, probs = 0.05),
  #        dec_trend_temp_95 = quantile(dec_trend_temp, probs = 0.95),
  #        dec_trend_temp = case_when(dec_trend_temp > dec_trend_temp_95 ~ dec_trend_temp_95,
  #                                   dec_trend_temp < dec_trend_temp_05 ~ dec_trend_temp_05,
  #                                   TRUE ~ dec_trend_temp)) %>% 
  ggplot() +
  geom_tile(aes(fill = dec_trend_temp, x = lon, y = lat)) +
  # geom_sf(data = MEOW, aes(colour = ECOREGION), fill = NA, show.legend = F) +
  geom_polygon(data = map_base, aes(group = group, x = lon, y = lat)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  labs(x = NULL, y = NULL, 
       fill = "SST trend\n(°C/dec.)", colour = "Ecoregion") +
  coord_sf(expand = F,
           xlim = c(min(med_regions$lon), max(med_regions$lon)),
           ylim = c(min(med_regions$lat), max(med_regions$lat))) +
  theme(legend.position = "top", 
        panel.background = element_rect(colour = "black", fill = NULL))
map_SST_trend

# Rank tables of ecoregions

# Table theme
t1 <- ttheme_default(core = list(bg_params = list(fill = c(rep(c("grey90", "grey95"), length.out = 8),
                                                           rep(c("grey100"), length.out = 1)))))

# Table
med_SST_table <- med_trends %>% 
  dplyr::rename(Province = PROVINCE, Rank = rank, Trend = dec_trend_temp_mean, SD = dec_trend_temp_sd) %>% 
  select(Rank, Province, Ecoregion, Trend, SD) %>% 
  tableGrob(rows = NULL, theme = t1) %>%
  gtable_add_grob(grobs = segmentsGrob(
    x0 = unit(0,"npc"),
    y0 = unit(0,"npc"),
    x1 = unit(1,"npc"),
    y1 = unit(0,"npc"),
    gp = gpar(lwd = 2.0)),
    t = 9, b = 9, l = 1, r = 5)
# plot(med_SST_table)

# Combine SST figures
fig_SST_maps <- ggpubr::ggarrange(map_SST_total, map_SST_trend, 
                                  ncol = 1, nrow = 2, align = "hv", 
                                  labels = c("A)", "B)"))
fig_SST_all <- ggpubr::ggarrange(fig_SST_maps, med_SST_table, ncol = 2, nrow = 1, 
                                 labels = c(NA, "C)"), widths = c(1, 0.7))
ggsave("figures/Med_SST.png", fig_SST_all, height = 11, width = 16.5)


# Figure 7 ---------------------------------------------------------------

# NB: This requires all the code to be run for Figure 6 section

# Map of SST mean per pixel SST 
map_MHW_dur_total <- SST_MHW_trends %>%
  na.omit() %>% 
  # Rather not remove the tails
  # mutate(dur_average_05 = quantile(dur_average, probs = 0.05),
  #        dur_average_95 = quantile(dur_average, probs = 0.95),
  #        dur_average = case_when(dur_average > dur_average_95 ~ dur_average_95,
  #                                dur_average < dur_average_05 ~ dur_average_05,
  #                                TRUE ~ dur_average)) %>% 
  ggplot() +
  geom_tile(aes(fill = dur_average, x = lon, y = lat)) +
  geom_sf(data = MEOW, aes(colour = ECOREGION), fill = NA, show.legend = F) +
  geom_polygon(data = map_base, aes(group = group, x = lon, y = lat)) +
  scale_fill_viridis_c(option = "A") +
  # scale_colour_brewer(palette = "Set1") +
  labs(x = NULL, y = NULL, fill = "MHW days\n(annual average)") +
  coord_sf(expand = F,
           xlim = c(min(med_regions$lon), max(med_regions$lon)),
           ylim = c(min(med_regions$lat), max(med_regions$lat))) +
  theme(legend.position = "top", 
        panel.background = element_rect(colour = "black", fill = NULL))
map_MHW_dur_total

# Map of decadal trend per pixel
map_MHW_dur_trend <- SST_MHW_trends %>%
  na.omit() %>% 
  filter(p.value_dur <= 0.05) %>% 
  # mutate(dec_trend_05 = quantile(dec_trend, probs = 0.05),
  #        dec_trend_95 = quantile(dec_trend, probs = 0.95),
  #        dec_trend = case_when(dec_trend > dec_trend_95 ~ dec_trend_95,
  #                              dec_trend < dec_trend_05 ~ dec_trend_05,
  #                              TRUE ~ dec_trend)) %>% 
  ggplot() +
  geom_tile(aes(fill = dec_trend_dur, x = lon, y = lat)) +
  # geom_sf(data = MEOW, aes(colour = ECOREGION), fill = NA, show.legend = F) +
  geom_polygon(data = map_base, aes(group = group, x = lon, y = lat)) +
  scale_fill_gradient2(low = "yellow4", mid = "white", high = "green4") +
  # scale_fill_viridis_c(option = "A") +
  # scale_colour_brewer(palette = "Set1") +
  labs(x = NULL, y = NULL, 
       fill = "Annual MHW\ndays (n/decade)", colour = "Ecoregion") +
  coord_sf(expand = F,
           xlim = c(min(med_regions$lon), max(med_regions$lon)),
           ylim = c(min(med_regions$lat), max(med_regions$lat))) +
  theme(legend.position = "top", 
        panel.background = element_rect(colour = "black", fill = NULL))
map_MHW_dur_trend

# Table
med_MHW_dur_table <- med_trends %>% 
  dplyr::rename(Province = PROVINCE, Rank = rank, Trend = dec_trend_dur_mean, SD = dec_trend_dur_sd) %>% 
  select(Rank, Province, Ecoregion, Trend, SD) %>% 
  tableGrob(rows = NULL, theme = t1) %>%
  gtable_add_grob(grobs = segmentsGrob(
    x0 = unit(0,"npc"),
    y0 = unit(0,"npc"),
    x1 = unit(1,"npc"),
    y1 = unit(0,"npc"),
    gp = gpar(lwd = 2.0)),
    t = 9, b = 9, l = 1, r = 5)
# plot(med_MHW_dur_table)

# Combine SST figures
fig_MHW_dur_maps <- ggpubr::ggarrange(map_MHW_dur_total, map_MHW_dur_trend, 
                                      ncol = 1, nrow = 2, align = "hv", 
                                      labels = c("A)", "B)"))
fig_MHW_dur_all <- ggpubr::ggarrange(fig_MHW_dur_maps, med_MHW_dur_table, ncol = 2, nrow = 1, 
                                     labels = c(NA, "C)"), widths = c(1, 0.7))
ggsave("figures/Med_MHW_duration.png", fig_MHW_dur_all, height = 11, width = 16.5)


# Manuscript figure 1 -----------------------------------------------------

## A: Map of temperature difference mean 1982-1986 vs 2015-2019

# Load data
load("data/MHW_clim_pixel_annual.RData")

# Custom legend
MHW_colours_compare <- c(
  "Same" = "salmon",
  MHW_colours[2],
  MHW_colours[3],
  MHW_colours[4]
)

# Manual label placements
# MEOW$ECOREGION
ecoregion_labels <- data.frame(lon = c(17.4, 33.6, 12.1, 20.8, 25.7, -2.26, 3.5, 4.1),
                               lat = c(44.7, 37.9, 31.5, 39.9, 42.1, 33.7, 44.5, 34.7),
                               Ecoregion = c("Adriatic\nSea", "Levantine\nSea", "Tunisian Plateau", "Ionian\nSea",
                                             "Aegean\nSea", "Alboran Sea", "NW Mediterranean", "SW Mediterranean"))

# Prep data
pixel_pentad <- MHW_clim_pixel_annual %>% 
  right_join(med_regions, by = c("lon", "lat")) %>% 
  mutate(pentad = cut(year, c(1981, 1986, 1992, 1998, 2003, 2009, 2014, 2019))) %>% 
  group_by(lon, lat, pentad) %>% 
  summarise(temp = mean(temp, na.rm = T), .groups = "drop") %>% 
  pivot_wider(id_cols = c("lon", "lat"), names_from = "pentad", values_from = "temp") %>% 
  mutate(temp_diff = `(2014,2019]` - `(1981,1986]`)

# Plot data
panel_A <- ggplot(pixel_pentad) +
  geom_tile(aes(fill = temp_diff, x = lon, y = lat), colour = NA) +
  geom_polygon(data = map_base, aes(x = lon, y = lat, group = group), fill = "grey60") +
  geom_sf(data = MEOW, alpha = 1, aes(geometry = geometry), fill = NA, colour = "forestgreen") +
  geom_point(data = insitu_sites, aes(x = lon, y = lat), size = 3, shape = 21, fill = "cyan") +
  geom_label(data = ecoregion_labels, aes(x = lon, y = lat, label = Ecoregion), alpha = 0.8) +
  scale_fill_gradient2(low = "yellow", mid = "orange", high = "red",
                       breaks = c(0.7, 1.2, 1.7), midpoint = 1.1) +
  coord_sf(expand = F,
           xlim = c(min(med_regions$lon), max(med_regions$lon)),
           ylim = c(min(med_regions$lat), max(med_regions$lat))) +
  scale_x_continuous(breaks = c(0, 10, 20, 30)) +
  scale_y_continuous(breaks = c(32, 36, 40, 44)) +
  labs(y = "Latitude", x = "Longitude", 
       # title = "Temperature difference [2015 to 2019] minus [1982 to 1986]",
       fill = "\U0394 Temp. (°C)") +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        panel.background = element_rect(fill = "white"),
        legend.position = c(0.91, 0.8),
        legend.background = element_rect(colour = "black"),
        legend.margin = margin(t = 5, r = 20, b = 5, l = 5),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))
panel_A

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
  labs(y = "Temperature (°C)", 
       # title = "Annual SST anomalies [1982 to 2019]",
       x = "Year") +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        panel.background = element_rect(fill = "white"),
        # legend.position = c(0.91, 0.83),
        legend.background = element_rect(colour = "black"),
        legend.margin = margin(t = 5, r = 15, b = 5, l = 5),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))
panel_B

## C: Map of difference in cat 2+ days between first and last pentad
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
  ggplot() +
  geom_tile(data = filter(pixel_cat_pentad, cat_max_diff_cat == "Same"), 
            aes(x = lon, y = lat, fill = cat_max_diff_cat), alpha = 0.3) +
  geom_tile(data = filter(pixel_cat_pentad, cat_max_diff_cat != "Same"),
            aes(x = lon, y = lat, fill = cat_max_diff_cat), colour = NA, show.legend = T) +
  geom_polygon(data = map_base, aes(x = lon, y = lat, group = group), fill = "grey60") +
  geom_sf(data = MEOW, alpha = 1, aes(geometry = geometry), fill = NA, colour = "forestgreen") +
  geom_point(data = insitu_sites, aes(x = lon, y = lat), size = 3, shape = 21, fill = "cyan") +
  geom_label(data = ecoregion_labels, aes(x = lon, y = lat, label = Ecoregion), alpha = 0.8) +
  scale_fill_manual("Category", values = MHW_colours_compare,
                    breaks = c("Same", "II Strong", "III Severe", "IV Extreme")) +
  coord_sf(expand = F,
           xlim = c(min(med_regions$lon), max(med_regions$lon)),
           ylim = c(min(med_regions$lat), max(med_regions$lat))) +
  scale_x_continuous(breaks = c(0, 10, 20, 30)) +
  scale_y_continuous(breaks = c(32, 36, 40, 44)) +
  labs(y = "Latitude", x = "Longitude", 
       # title = "The highest category for [2015 to 2019] when [1982 to 1986] was 'I Moderate' or less",
       fill = "Max category") +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        panel.background = element_rect(fill = "white"),
        legend.position = c(0.91, 0.83),
        legend.background = element_rect(colour = "black"),
        legend.margin = margin(t = 5, r = 15, b = 5, l = 5),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))
panel_C

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
  geom_bar(aes(fill = category), stat = "identity", show.legend = F,
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
  labs(y = "Cover (%)",
       # title = "Surface area affected by category 'II Strong'+ MHWs",
       x = "Year") +
  coord_cartesian(expand = F) +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        panel.background = element_rect(fill = "white"),
        legend.margin = margin(t = 5, r = 15, b = 5, l = 5),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))
panel_D

## Combine and save
manu_fig_1 <- ggpubr::ggarrange(panel_A, panel_B, panel_C, panel_D, labels = c("(a)", "(b)", "(c)", "(d)"))
ggsave("figures/manu_fig_1.png", manu_fig_1, height = 10, width = 20.3)
ggsave("figures/manu_fig_1.pdf", manu_fig_1, height = 10, width = 20.3)

# Summary stats for text
cat_daily_mean %>% 
  group_by(year) %>% 
  summarise(first_n_cum_prop_sum = sum(first_n_cum_prop)) %>% 
  arrange(-first_n_cum_prop_sum)

# Difference in SST from first to last pentad
med_pentad$anom_pentad[7]-med_pentad$anom_pentad[1]

# Proportion of surface area affected by specific categories over the study period
## Prop I
nrow(distinct(dplyr::select(filter(MHW_cat_pixel_annual, year >= 2015, year <= 2019,  `I Moderate` > 0), lon, lat)))/nrow(MHW_pixels)
nrow(distinct(dplyr::select(filter(MHW_cat_pixel_annual, year >= 2015, year <= 2019,  `II Strong` > 0), lon, lat)))/nrow(MHW_pixels)
nrow(distinct(dplyr::select(filter(MHW_cat_pixel_annual, year >= 2015, year <= 2019,  `III Severe` > 0), lon, lat)))/nrow(MHW_pixels)
nrow(distinct(dplyr::select(filter(MHW_cat_pixel_annual, year >= 2015, year <= 2019,  `IV Extreme` > 0), lon, lat)))/nrow(MHW_pixels)

# Years with Extreme events
filter(MHW_cat_pixel_annual, `IV Extreme` > 0) %>% 
  dplyr::select(year) %>% 
  distinct() %>% 
  arrange(year)


# Manuscript figure 4 -----------------------------------------------------

# Use only the species from AX column: YES
# Only take records that showed mortality
mme_selected_2_mort <- filter(mme, selected_2 == "YES") %>%
# mme_3B_mort <- filter(mme, Plot_3B %in% c("2015_MHW", "2016_MHW", "2017_MHW", "2018_MHW", "2019_MHW")) %>% # NB: This works much worse
  mutate(mme_damage = case_when(`Damaged qualitative` == "No" ~ 0, TRUE ~ 1)) %>% 
  group_by(Ecoregion, year) %>% 
  summarise(mme_record = n(),
            mme_prop = sum(mme_damage)/mme_record, .groups = "drop")

# MHW stats for JJASON period
load("data/MHW_cat_region.RData")
MHW_cat_region_JJASON <- MHW_cat_region %>%
  filter(year >= 2015, month %in% c("juin", "juil", "août", "sept", "oct", "nov")) %>%
  mutate(duration = duration/pixels) %>%
  group_by(region, year) %>%
  summarise(duration = sum(duration), .groups = "drop")

# Create region averages by all pixels
mme_region <- mme_selected_2_mort %>%
# mme_region <- mme_3B_mort %>% # This is much morse
  group_by(Ecoregion, year) %>%
  mutate(mme_prop = (mme_prop*100)+20) %>%
  filter(Ecoregion == "Northwestern Mediterranean")
manu_fig_4_a_data <- left_join(mme_region, MHW_cat_region_JJASON, by = c("Ecoregion" = "region", "year"))
write_csv(manu_fig_4_a_data, "data/manu_fig_4_a_data.csv")

# Create plot
manu_fig_4_a <- ggplot(manu_fig_4_a_data, aes(x = year, y = duration)) +
  geom_bar(aes(fill = as.factor(year)), colour = "black", stat = "identity", 
           show.legend = F, width = 0.8) +
  geom_point(aes(y = mme_prop), size = 5) +
  scale_fill_viridis_d("Year", option = "D", aesthetics = c("colour", "fill")) +
  scale_y_continuous(limits = c(0, 110), breaks = c(20, 40, 60, 80, 100), expand = c(0, 0),
                     sec.axis = sec_axis(name = "Overall proportion of records\nwith mass mortality", 
                                         trans = ~ . + 0,
                                         breaks = c(20, 45, 70, 95),
                                         labels = c("0.00", "0.25", "0.50", "0.75"))) +
  labs(x = NULL,
       title = "N.W. Mediterranean",
       # title = "Average MHW days per year in NW Med",
       # subtitle = "Points show proportion of records with mortality",
       y = "MHW days (bars)",) +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        panel.background = element_rect(fill = "white"),
        plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
        # panel.grid = element_line(colour = "black"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        axis.text = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(face = "bold", size = 14),
        axis.title = element_text(size = 14))
manu_fig_4_a

# Analysis of mortality per pixel per ecoregion
# Get MHW days and relate them to percentage of observations that have mortality
load("data/MHW_cat_pixel_annual_JJASON.RData")

## Use filter 3B in the NW Med for these analyses
mme_3B <- filter(mme, Plot_3B %in% c("2015_MHW", "2016_MHW", "2017_MHW", "2018_MHW", "2019_MHW")) %>%
  filter(Ecoregion == "Northwestern Mediterranean") %>%
  mutate(mme_damage = case_when(`Damaged qualitative` %in% c("No") ~ 0, TRUE ~ 1))

# Match nearest pixels
mme_mhw_pixel_match <- grid_match(mme_3B[c("lon", "lat")],
                                  MHW_pixels[c("lon", "lat")]) %>% 
  dplyr::rename(lon = lon.x, lat = lat.x, lon_sst = lon.y, lat_sst = lat.y) %>% 
  distinct()

# Calculate data for plotting
manu_fig_4_b_data <- mme_3B %>% 
  left_join(mme_mhw_pixel_match, by = c("lon", "lat")) %>% 
  left_join(MHW_cat_pixel_annual_JJASON, by = c("year", "lon_sst" = "lon", "lat_sst" = "lat")) %>% 
  mutate(duration_sum = replace_na(duration_sum, 0)) %>% 
  group_by(Ecoregion, `Area Monitored`, year) %>%
  summarise(NMME = sum(mme_damage),
            NEVENTS = n(),
            mme_prop = NMME/NEVENTS,
            duration_sum = mean(duration_sum, na.rm = T),
            damage_mean = mean(`Damaged percentage`, na.rm = T), .groups = "drop") %>% 
  dplyr::select(Ecoregion, `Area Monitored`, year, NMME, NEVENTS, mme_prop, damage_mean, duration_sum) %>% 
  filter(NEVENTS >= 3)
write_csv(manu_fig_4_b_data, "data/manu_fig_4_b_data.csv")

# Number of SST pixels per area monitored
mme_mhw_pixels_per_area <- mme_3B %>% 
  left_join(mme_mhw_pixel_match, by = c("lon", "lat")) %>% 
  left_join(MHW_cat_pixel_annual_JJASON, by = c("year", "lon_sst" = "lon", "lat_sst" = "lat")) %>% 
  dplyr::select(Ecoregion, `Area Monitored`, year, lon_sst, lat_sst) %>% 
  group_by(Ecoregion, `Area Monitored`, year) %>% 
  summarise(sst_pixels = n(), .groups = "drop")

# Perform analysis against all MME above the lowest category etc.
manu_fig_4_b_label <- manu_fig_4_b_data %>%
  summarise(count = n(),
            r_val = round(cor.test(mme_prop, duration_sum)$estimate, 2),
            p_val = round(cor.test(mme_prop, duration_sum)$p.value, 2),
            x_point = sum(range(duration_sum, na.rm = T))/2, .groups = "drop") %>% 
  mutate(p_val = case_when(p_val < 0.001 ~ "p < 0.001",
                           TRUE ~ paste0("p = ",p_val)))
manu_fig_4_b <- manu_fig_4_b_data %>% 
  # filter(`Damaged qualitative` != "No") %>% 
  ggplot(aes(x = duration_sum, y = mme_prop)) +
  geom_smooth(method = "lm", se = F, colour = "blue") +
  geom_point(shape = 21, size = 5, alpha = 0.9,
             # aes(fill = as.factor(year)),
             show.legend = F) +
  # Trick to get a black border around a blue text label
  geom_label(data = manu_fig_4_b_label, alpha = 0.9, colour = "black",
             aes(y = 0.2, x = 70, label = paste0("r = ",r_val,", ",p_val, ", n = ",count))) +
  # The coloured text label
  geom_label(data = manu_fig_4_b_label, alpha = 0.9, colour = "blue", label.size = 0, label.padding = unit(0, "mm"),
             aes(y = 0.2, x = 70, label = paste0("r = ",r_val,", ",p_val, ", n = ",count))) +
  scale_fill_viridis_d("Year", option = "D", aesthetics = c("colour", "fill")) +
  scale_y_continuous(limits = c(-0.1, 1.1), breaks = c(0, 0.25, 0.50, 0.75, 1), 
                     labels = c("0", "0.25", "0.5", "0.75", "1.0")) +
  coord_cartesian(expand = F) +
  labs(y = "Proportion of records with MME\nwithin each monitored area",
       title = "Within monitored areas",
       # title = "MHW days and MME records per monitoring area in the NW Med",
       # subtitle = "MHW days per year during JJASON period and proportion of MME records with damage",
       x = "MHW days (per year and monitored area)") +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        panel.background = element_rect(fill = "white"),
        plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
        # panel.grid = element_line(colour = "black"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 14))
manu_fig_4_b

manu_fig_4 <- ggpubr::ggarrange(manu_fig_4_a, manu_fig_4_b, labels = c("(a)", "(b)"), ncol = 1, align = "hv")
ggsave("figures/manu_fig_4.png", manu_fig_4, height = 10, width = 8)
ggsave("figures/manu_fig_4.pdf", manu_fig_4, height = 10, width = 8)

