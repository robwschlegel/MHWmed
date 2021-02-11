# code/MHW_NetCDF.R
# This script contains the code used to convert MHW results to NetCDF files


# Setup -------------------------------------------------------------------

# Needed packages
library(dplyr) # For basic data manipulation
library(tidyr) # For unnesting
library(ncdf4) # For creating NetCDF files
library(tidync) # For easily dealing with NetCDF data
library(ggplot2) # For visualising data
library(doParallel) # For parallel processing

# Original SST files for reference
med_SST_files <- dir("data/SST", pattern = ".nc", full.names = T, recursive = T)

# The lon/lat indexes
med_lat <- tidync(med_SST_files[1]) %>% 
  activate("D1") %>% 
  hyper_tibble()
med_lon <- tidync(med_SST_files[1]) %>% 
  activate("D2") %>% 
  hyper_tibble()
med_coords <- expand.grid(lat = as.vector(med_lat$lat), 
                          lon = as.vector(med_lon$lon)) %>% 
  data.frame()

# File locations
MHW_files <- dir("data/MHW", full.names = T)


# Functions ---------------------------------------------------------------

# Load an Rds MHW result file and return only event+cat data
Rds_to_event_cat <- function(file_name){
  
  # Load Rds MHW result file
  MHW_res <- readRDS(file_name)
  
  # Extract the four different data.frames
  event_event <- MHW_res %>% 
    dplyr::select(-cat) %>% 
    unnest(event) %>% 
    filter(row_number() %% 2 == 0) %>% 
    unnest(event) %>% 
    ungroup()
  cat_event <- MHW_res %>% 
    dplyr::select(-event) %>% 
    unnest(cat) %>% 
    filter(row_number() %% 2 == 0) %>% 
    unnest(cat) %>% 
    ungroup()
  rm(MHW_res); gc()
  
  # Combine
  all_event <- left_join(event_event, cat_event, 
                         by = c("lon", "lat", "event_no", "duration", 
                                "intensity_max" = "i_max", "date_peak" = "peak_date")) %>%
    dplyr::select(lon:event_no, duration:intensity_max, intensity_cumulative,
                  rate_onset, rate_decline, category) %>%
    mutate(category = factor(category, levels = c("I Moderate", "II Strong", "III Severe", "IV Extreme")))
  # filter(lon < -18) # For testing
  rm(event_event, cat_event); gc()
  return(all_event)
}

# Function for creating arrays from data.frames
df_acast <- function(df, lon_lat){
  
  # Force grid
  res <- df %>%
    right_join(lon_lat, by = c("lon", "lat")) %>%
    arrange(lon, lat)
  
  # Convert date values to integers if they are present
  if(lubridate::is.Date(res[1,4])) res[,4] <- as.integer(res[,4])
  
  # Convert factors to integers if they are present
  if(is.factor(res[1,4])) res[,4] <- as.integer(res[,4])
  
  # Create array
  res_array <- base::array(res[,4], dim = c(length(unique(lon_lat$lon)), length(unique(lon_lat$lat))))
  dimnames(res_array) <- list(lon = unique(lon_lat$lon),
                              lat = unique(lon_lat$lat))
  return(res_array)
}

# Wrapper function for last step before data are entered into NetCDF files
df_proc <- function(df, col_choice, lon_lat = NA){
  
  if(is.na(lon_lat[1,1])){
    # Determine the correct array dimensions
    lon_step <- mean(diff(sort(unique(df$lon))))
    lat_step <- mean(diff(sort(unique(df$lat))))
    lon <- seq(min(df$lon), max(df$lon), by = lon_step)
    lat <- seq(min(df$lat), max(df$lat), by = lat_step)
    
    # Create full lon/lat grid
    lon_lat <- expand.grid(lon = lon, lat = lat) %>% 
      data.frame()
  }
  
  # Acast only the desired column
  dfa <- plyr::daply(df[c("lon", "lat", "event_no", col_choice)], 
                     c("event_no"), df_acast, .parallel = T, lon_lat = lon_lat)
  return(dfa)
}

# Save event data.frame to NetCDF
# event_data <- Rds_to_event_cat(MHW_files[1])
event_to_NetCDF <- function(event_data, lon_lat = NA){
  
  # Process data for adding to NetCDF file
  # We must now run this function on each column of data we want to add to the NetCDF file
  doParallel::registerDoParallel(cores = 14)
  print(paste0("Began prepping dates at ", Sys.time()))
  prep_start <- df_proc(event_data, "date_start", lon_lat)
  prep_peak <- df_proc(event_data, "date_peak", lon_lat)
  prep_end <- df_proc(event_data, "date_end", lon_lat)
  prep_dur <- df_proc(event_data, "duration", lon_lat)
  print(paste0("Began prepping intensities at ", Sys.time()))
  prep_mean_int <- df_proc(event_data, "intensity_mean", lon_lat)
  prep_max_int <- df_proc(event_data, "intensity_max", lon_lat)
  prep_cum_int <- df_proc(event_data, "intensity_cumulative", lon_lat)
  print(paste0("Began prepping other at ", Sys.time()))
  prep_onset <- df_proc(event_data, "rate_onset", lon_lat)
  prep_decline <- df_proc(event_data, "rate_decline", lon_lat)
  prep_cat <- df_proc(event_data, "category", lon_lat)
  
  # Check array
  # dim(prep_dur)
  # dimnames(prep_cat)
  # plot(prep_cat[,,1])
  # prep_test <- data.frame(prep_mean_int[,1,1])
  
  # Get file attributes
  lon <- unique(lon_lat$lon)
  lat <- unique(lon_lat$lat)
  event_no <- unique(event_data$event_no)
  tunits <- "days since 1970-01-01"
  
  # Length of each attribute
  nlon <- length(lon)
  nlat <- length(lat)
  nen <- length(event_no)
  
  # path and file name, set dname
  ncpath <- "~/pCloudDrive/MHWmed_data/"
  ncname <- "MHW_event"
  ncfname <- paste0(ncpath, ncname, ".nc")
  
  # define dimensions
  londim <- ncdim_def("lon", "degrees_east", lon)
  latdim <- ncdim_def("lat", "degrees_north", lat)
  endim <- ncdim_def("event_no", "event_number", event_no)
  
  # define variables
  fillvalue <- 9999
  def_start <- ncvar_def(name = "date_start", units = tunits, dim = list(endim, latdim, londim),
                         missval = 0, longname = "start date of MHW", prec = "integer")
  def_peak <- ncvar_def(name = "date_peak", units = tunits, dim = list(endim, latdim, londim),
                        missval = 0, longname = "date of peak temperature anomaly during MHW", prec = "integer")
  def_end <- ncvar_def(name = "date_end", units = tunits, dim = list(endim, latdim, londim),
                        missval = 0, longname = "end date of MHW", prec = "integer")
  def_dur <- ncvar_def(name = "duration", units = "days", dim = list(endim, latdim, londim), 
                       missval = fillvalue, longname = "duration of MHW", prec = "integer")
  def_mean_int <- ncvar_def(name = "mean_int", units = "deg_C", dim = list(endim, latdim, londim), 
                            missval = fillvalue, longname = "mean intensity during MHW", prec = "double")
  def_max_int <- ncvar_def(name = "max_int", units = "deg_C", dim = list(endim, latdim, londim), 
                           missval = fillvalue, longname = "maximum intensity during MHW", prec = "double")
  def_cum_int <- ncvar_def(name = "cum_int", units = "deg_C days", dim = list(endim, latdim, londim), 
                           missval = fillvalue, longname = "cumulative intensity during MHW", prec = "double")
  def_onset <- ncvar_def(name = "rate_onset", units = "deg_C", dim = list(endim, latdim, londim), 
                           missval = fillvalue, longname = "daily increase in temperature from start to peak of MHW", prec = "double")
  def_decline <- ncvar_def(name = "rate_decline", units = "deg_C", dim = list(endim, latdim, londim),
                           missval = fillvalue, longname = "daily increase in temperature from start to peak of MHW", prec = "double")
  def_cat <- ncvar_def(name = "category", units = "Category", dim = list(endim, latdim, londim),
                       missval = 0, longname = "highest category reached during MHW", prec = "integer")
  
  # create netCDF file and put arrays
  ncout <- nc_create(ncfname, list(def_start, def_peak, def_end, def_dur, 
                                   def_mean_int, def_max_int, def_cum_int,
                                   def_onset, def_decline, def_cat), force_v4 = TRUE)
  
  # put variables
  print(paste0("Began filling NetCDF file at ", Sys.time()))
  ncvar_put(ncout, def_start, prep_start)
  ncvar_put(ncout, def_peak, prep_peak)
  ncvar_put(ncout, def_end, prep_end)
  ncvar_put(ncout, def_dur, prep_dur)
  ncvar_put(ncout, def_mean_int, prep_mean_int)
  ncvar_put(ncout, def_max_int, prep_max_int)
  ncvar_put(ncout, def_cum_int, prep_cum_int)
  ncvar_put(ncout, def_onset, prep_onset)
  ncvar_put(ncout, def_decline, prep_decline)
  ncvar_put(ncout, def_cat, prep_cat)
  
  # put additional attributes into dimension and data variables
  ncatt_put(ncout, "lon", "axis", "X") #,verbose=FALSE) #,definemode=FALSE)
  ncatt_put(ncout, "lat", "axis", "Y")
  ncatt_put(ncout, "event_no", "axis", "event_no")
  
  # add global attributes
  ncatt_put(ncout, 0, "title", paste0("MHW results from lon: ",
                                      min(lon)," to ",max(lon),
                                      " and lat: ",min(lat)," to ",max(lat)))
  ncatt_put(ncout, 0, "institution", "Institut de la Mer de Villefranche")
  ncatt_put(ncout, 0, "source", "NOAA Med SST ~4km")
  ncatt_put(ncout,0, "references", "Banzon et al. (2020) J. Atmos. Oce. Tech. 37:341-349")
  history <- paste0("Robert W Schlegel, ", date())
  ncatt_put(ncout, 0, "history", history)
  ncatt_put(ncout, 0, "Conventions", "Hobday et al. (2016; 2018)") # Assuming one has used the Hobday definition
  
  # Get a summary of the created file:
  ncout
}


# Convert files -----------------------------------------------------------

# Load the full brick of MHW event+cat results
registerDoParallel(cores = 14)
system.time(
MHW_event <- plyr::ldply(MHW_files, Rds_to_event_cat, .parallel = T)
) # 78 seconds

# Create the NetCDF file
system.time(
event_to_NetCDF(MHW_event, lon_lat = med_coords)
) # xxx seconds


# Test the output ---------------------------------------------------------

# Convenience function for comparing files
quick_grid <- function(df, var_choice){
  df %>% 
    filter(event_no == 13) %>% 
    ggplot(aes(x = lon, y = lat)) +
    geom_raster(aes_string(fill = var_choice)) +
    coord_cartesian(expand = F) +
    scale_fill_viridis_c() +
    theme(legend.position = "bottom")
}

# Thanks to the tidync package, loading the gridded data is very simple
MHW_res_nc <- tidync("~/pCloudDrive/MHWmed_data/MHW_event.nc") %>% 
  hyper_tibble() #%>% 
  # mutate(date_peak = as.Date(date_peak, origin = "1970-01-01"))

# Plot the duration results
quick_grid(MHW_event, "duration")
quick_grid(MHW_res_nc, "duration")

# Cumulative intensity
quick_grid(MHW_event, "intensity_cumulative")
quick_grid(MHW_res_nc, "cum_int")

# Rate of onset
quick_grid(MHW_event, "rate_onset")
quick_grid(MHW_res_nc, "rate_onset")

# Category
MHW_event %>% 
  mutate(category = as.integer(category)) %>% 
  quick_grid("category")
quick_grid(MHW_res_nc, "category")
