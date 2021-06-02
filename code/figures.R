# code/figures.R
# This script creates annual and total summaries from the MHW Med results


# Setup -------------------------------------------------------------------

# The project-wide functions
source("code/functions.R")


# Figure 1 ----------------------------------------------------------------

# Total annual Med MHW summary with global overlay
total_summary <- total_summary_fig(MHW_cat_summary_annual)
ggsave("figures/fig_1.png", total_summary, height = 4.25, width = 8)


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

# MHW annual stats
mhw_eco_summary_annual <- MHW_cat_region_coast %>% 
  dplyr::rename(Ecoregion = region) %>% 
  dplyr::select(-month) %>% 
  mutate(duration = duration/pixels,
         cum_int = cum_int/pixels) %>% 
  group_by(Ecoregion, year) %>% 
  summarise(surface = mean(surface, na.rm = T),
            duration = sum(duration, na.rm = T),
            icum = sum(cum_int, na.rm = T), .groups = "drop")

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

# Requires: site_MME_MHW_summary
# site_MME_MHW_summary <- read_csv("data/site_MME_MHW_summary.csv")

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
eco_MME_MHW_annual <- left_join(full_region_year_grid, mme_eco_summary, by = join_vals) %>% 
  left_join(mme_eco_sites_count, by = join_vals) %>% 
  left_join(mhw_eco_summary_annual, by = join_vals) %>% 
  replace(is.na(.), 0)
eco_MME_MHW_JJASON <- left_join(full_region_year_grid, mme_eco_summary, by = join_vals) %>% 
  left_join(mme_eco_sites_count, by = join_vals) %>% 
  left_join(mhw_eco_summary_JJASON, by = join_vals) %>% 
  replace(is.na(.), 0)

# Create a whole Med summary
eco_MME_MHW_med_annual <- eco_MME_MHW_annual %>% 
  group_by(year) %>% 
  summarise(`Damaged percentage` = mean(`Damaged percentage`, na.rm = T),
            mme_count = sum(mme_count, na.rm = T),
            site_count = sum(site_count, na.rm = T),
            surface = mean(surface, na.rm = T),
            duration = mean(duration, na.rm = T),
            icum = mean(icum, na.rm = T), .groups = "drop") %>% 
  mutate(Ecoregion = "Mediterranean")
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
eco_MME_MHW_annual_all <- rbind(eco_MME_MHW_annual, eco_MME_MHW_med_annual) %>% 
  mutate(Ecoregion = factor(Ecoregion, 
                            levels = c("Mediterranean",
                                       "Alboran Sea", "Northwestern Mediterranean", 
                                       "Southwestern Mediterranean", "Adriatic Sea",
                                       "Ionian Sea", "Tunisian Plateau/Gulf of Sidra",
                                       "Aegean Sea", "Levantine Sea")))
eco_MME_MHW_JJASON_all <- rbind(eco_MME_MHW_JJASON, eco_MME_MHW_med_JJASON) %>% 
  mutate(Ecoregion = factor(Ecoregion, 
                            levels = c("Mediterranean",
                                       "Alboran Sea", "Northwestern Mediterranean", 
                                       "Southwestern Mediterranean", "Adriatic Sea",
                                       "Ionian Sea", "Tunisian Plateau/Gulf of Sidra",
                                       "Aegean Sea", "Levantine Sea")))

# Barplot of durations
bar_dur_fig <- function(df, title_bit){
  res_fig <- df %>% 
    ggplot(aes(x = year, y = duration)) +
    geom_bar(aes(fill = icum), 
             colour = "black",
             stat = "identity", 
             show.legend = T,
             position = "dodge",
             # position = position_stack(reverse = TRUE), 
             width = 1) +
    # geom_label(aes(label = count_MME_mean)) +
    geom_label(aes(label = paste0(mme_count,"/",site_count))) +
    scale_fill_viridis_c("Cumulative\nIntensity (°C days)", option = "B") +
    facet_wrap(~Ecoregion) +
    # scale_y_continuous(limits = c(-5, 125), breaks = c(25, 50, 75, 100)) +
    scale_x_continuous(breaks = c(2015, 2017, 2019)) +
    coord_cartesian(expand = F) +
    labs(x = NULL, y = "MHW days", 
         title = paste0("Total MHW days and cumulative intensity (°C days) per ecoregion/year", title_bit),
         subtitle = "Labels show count of MME with damage at all depths: count/sites") +
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
  # res_fig
  return(res_fig)
}

# Create figures
fig_2_annual <- bar_dur_fig(eco_MME_MHW_annual_all, " (annual)")
fig_2_JJASON <- bar_dur_fig(eco_MME_MHW_JJASON_all, " (JJASON)")

# Save
ggsave("figures/fig_2_annual.png", fig_2_annual, height = 8, width = 8)
ggsave("figures/fig_2_JJASON.png", fig_2_JJASON, height = 8, width = 8)


# Figure 3 ----------------------------------------------------------------

# Scatterplots showing iCum and MME % damage
# Use coastal pixels


# Figure 4 ----------------------------------------------------------------


