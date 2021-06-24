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
mme_reg <- filter(mme_mhw, selected_5 %in% c("2015_MHW", "2016_MHW", "2017_MHW", "2018_MHW", "2019_MHW")) %>% 
  filter(`Monitoring series` %in% c("more.than.two.per.year", "one.per.year.monitoring") | Ecoregion == "Alboran Sea")

# Create data.frames based on four pre-determined filter columns
mme_Plot_1A <- filter(mme_mhw, Plot_1A %in% c("2015_MHW", "2016_MHW", "2017_MHW", "2018_MHW", "2019_MHW"))
mme_Plot_1B <- filter(mme_mhw, Plot_1B %in% c("2015_MHW", "2016_MHW", "2017_MHW", "2018_MHW", "2019_MHW"))
mme_Plot_2A <- filter(mme_mhw, Plot_2A %in% c("2015_MHW", "2016_MHW", "2017_MHW", "2018_MHW", "2019_MHW"))
mme_Plot_2B <- filter(mme_mhw, Plot_2B %in% c("2015_MHW", "2016_MHW", "2017_MHW", "2018_MHW", "2019_MHW"))

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


