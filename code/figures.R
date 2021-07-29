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

# Barplots showing duration/iCum and MME
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


# Manuscript figure 1 -----------------------------------------------------

## A: Map of temperature difference mean 1982-1986 vs 2015-2019
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
  scale_fill_gradient2(low = "blue", high = "red", 
                       breaks = c(0, 0.5, 1.0, 1.5, 2.0),
                       # labels = c("0", "0.5", "1.0", "1.5", "2.0")) +
                       labels = c("0", "0.5", "1.0", "1.5", "2.0")) +
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
        legend.position = c(0.9, 0.8),
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
  scale_fill_gradient2(low = "blue", high = "red") +
  scale_x_continuous(breaks = seq(1984, 2019, 7), expand = c(0, 0)) +
  labs(title = "__(b)__   Annual SST anomalies (1982 - 2019)",
       y = "Temperature (°C)", x = "Year") +
  theme(plot.title = ggtext::element_markdown(),
        panel.border = element_rect(colour = "black", fill = NA),
        # legend.position = "bottom",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        panel.background = element_rect(fill = "grey90"))
# panel_B

## C: Map of areas affected by most intense MHWs 2015-2019
# Prep data
MHW_cat_pixel_annual_sub <- MHW_cat_pixel_annual %>%
  filter(year %in% seq(2015, 2019))
# Plot data
panel_C <- MHW_cat_pixel_annual_sub %>%
  group_by(lon, lat) %>%
  summarise(category = max(as.numeric(category), na.rm = T), .groups = "drop") %>%
  mutate(category = factor(category, labels = c("I Moderate", "II Strong", "III Severe", "IV Extreme"))) %>%
  ggplot(aes(x = lon, y = lat)) +
  # geom_tile(data = OISST_ice_coords, fill = "powderblue", colour = NA, alpha = 0.5) +
  geom_tile(aes(fill = category), colour = NA, show.legend = F) +
  geom_polygon(data = map_base, aes(x = lon, y = lat, group = group),
               fill = "grey70", colour = "black") +
  scale_fill_manual("Category", values = MHW_colours) +
  # scale_y_continuous(breaks = NULL) +
  # scale_x_continuous(breaks = NULL) +
  coord_cartesian(expand = F,
                  xlim = c(min(med_regions$lon), max(med_regions$lon)),
                  ylim = c(min(med_regions$lat), max(med_regions$lat))) +
  # theme_void() +
  # theme_bw() +
  labs(title = "__(c)__   Highest MHW categories from 2015-2019",
       y = "Latitude (°N)", x = "Longitude (°E)") +
  # guides(fill = guide_legend(override.aes = list(size = 10))) +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        plot.title = ggtext::element_markdown(),
        legend.position = "bottom",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        panel.background = element_rect(fill = "grey90"))
# panel_C

## D: Barplot of Med surface area affected by Cat 2+ MHWs 
# Load data
load("data/MHW_cat_summary_annual.RData")
OISST_global <- readRDS("data/OISST_cat_daily_1992-2018_total.Rds") %>% 
  group_by(t) %>% 
  mutate(cat_n_prop_stack = cumsum(cat_n_prop),
         first_n_cum_prop_stack = cumsum(first_n_cum_prop)) %>% 
  filter(category == "IV Extreme")
# Prep data
cat_daily_mean <- MHW_cat_summary_annual %>%
  group_by(year, category) %>%
  summarise(cat_n_prop_mean = mean(cat_n_prop, na.rm = T),
            cat_n_cum_prop = max(cat_n_cum_prop, na.rm = T), .groups = "drop")
cat_pentad <- cat_daily_mean %>% 
  group_by(year) %>% 
  summarise(cat_n_cum_prop_sum = sum(cat_n_cum_prop, na.rm = T), .groups = "drop") %>% 
  mutate(pentad = cut(year, c(1981, 1986, 1992, 1998, 2003, 2009, 2014, 2019))) %>% 
  group_by(pentad) %>% 
  summarise(cat_n_cum_prop_pentad = mean(cat_n_cum_prop_sum, na.rm = T), .groups = "drop") %>%
  separate(pentad, into = c("start_year", "end_year"), sep = ",", remove = F) %>% 
  mutate(start_year = as.numeric(sub("[(]", "", start_year)) + 1,
         end_year = as.numeric(sub("]", "", end_year)))
# Plot data
panel_D <- ggplot(cat_daily_mean, aes(x = year, y = cat_n_cum_prop)) +
  geom_bar(aes(fill = category), stat = "identity", show.legend = T,
           position = position_stack(reverse = TRUE), width = 1) +
  geom_segment(data = cat_pentad, size = 2, lineend = "round",
               aes(x = start_year, xend = end_year, 
                   y = cat_n_cum_prop_pentad, yend = cat_n_cum_prop_pentad)) +
  geom_point(data = OISST_global, aes(x = t, y = cat_n_prop_stack), 
             shape = 21, fill = "grey", show.legend = F) +
  scale_fill_manual("Category", values = MHW_colours) +
  scale_colour_manual("Category", values = MHW_colours) +
  scale_y_continuous(limits = c(0, 100),
                     breaks = seq(20, 80, length.out = 4)) +
  scale_x_continuous(breaks = seq(1984, 2019, 7)) +
  guides(pattern_colour = FALSE, colour = FALSE) +
  labs(title = "__(d)__   Surface area of Med affected by MHWs",
       y = "Cover (%)", x = "Year") +
  coord_cartesian(expand = F) +
  theme(plot.title = ggtext::element_markdown(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.position = c(0.1, 0.8),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))
# panel_D

## Combine and save
manu_fig_1 <- ggpubr::ggarrange(panel_A, panel_B, panel_C, panel_D)
ggsave("figures/manu_fig_1.png", manu_fig_1, height = 10, width = 16)

