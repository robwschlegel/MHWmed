# code/MHW_NetCDF.R
# This script contains the code used to convert MHW results to NetCDF files


# Setup -------------------------------------------------------------------

# Needed packages
library(tidyverse)
library(tidync)
library(ncdf4)
library(abind)

# For test visuals
library(lattice)
library(RColorBrewer)


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

# Function for creating arrays from data.frames
# df <- filter(df_step, t == 5786)
df_acast <- function(df){
  
  # Ensure correct grid size
  med_coords_sub <- med_coords %>% 
    filter(lat == df$lat[1])
  
  # Round data for massive file size reduction
  df$temp <- round(df$temp, 2)
  
  # Force grid
  res <- df %>%
    right_join(med_coords_sub, by = c("lon", "lat")) %>% 
    arrange(lon, lat)
  
  # Create array
  res_array <- base::array(res$temp, dim = c(1305,1,1))
  dimnames(res_array) <- list(lon = unique(med_coords_sub$lon),
                              lat = unique(med_coords_sub$lat),
                              t = unique(na.omit(res$t)))
  return(res_array)
}

# Wrapper function for last step before data are entered into NetCDF files
# df <- all_clim
df_proc <- function(df){
  
  # Filter NA and convert dates to integer
  df_step <- df[,1:5] %>% 
    mutate(temp = ifelse(is.na(temp), NA, temp),
           t = as.integer(t)) %>% 
    na.omit()
  
  # Acast
  dfa <- df_step %>%
    mutate(t2 = t) %>% 
    group_by(t2) %>%
    nest() %>%
    mutate(data2 = purrr::map(data, df_acast)) %>%
    select(-data)
  
  # Final form
  dfa_res <- abind(dfa$data2, along = 3, hier.names = T)
  rm(df_step, dfa); gc()
  # dimnames(dfa_res)
  return(dfa_res)
}

# Load a Rds file and save it as NetCDF
# file_name <- MHW_files[1]
Rds_to_NetCDF <- function(file_name){

  # Load R file
  MHW_res <- readRDS(file_name)
  
  # Extract the four different data.frames
  event_clim <- MHW_res %>% 
    dplyr::select(-cat) %>% 
    unnest(event) %>% 
    filter(row_number() %% 2 == 1) %>% 
    unnest(event) %>% 
    ungroup()
  event_event <- MHW_res %>% 
    dplyr::select(-cat) %>% 
    unnest(event) %>% 
    filter(row_number() %% 2 == 0) %>% 
    unnest(event) %>% 
    ungroup()
  cat_clim <- MHW_res %>% 
    dplyr::select(-event) %>% 
    unnest(cat) %>% 
    filter(row_number() %% 2 == 1) %>% 
    unnest(cat) %>% 
    ungroup()
  cat_event <- MHW_res %>% 
    dplyr::select(-event) %>% 
    unnest(cat) %>% 
    filter(row_number() %% 2 == 0) %>% 
    unnest(cat) %>% 
    ungroup()
  
  # Combine same shaped dataframes
  all_clim <- left_join(event_clim, cat_clim, by = c("lon", "lat", "t", "event_no")) #%>% 
    # filter(lon < -18) # For testing
  all_event <- left_join(event_event, cat_event, 
                         by = c("lon", "lat", "event_no", "duration", 
                                "intensity_max" = "i_max", "date_peak" = "peak_date")) #%>% 
    # filter(lon < -18) # For testing
  
  # Process data for adding to NetCDF file
  # Convert single columns in data.frames into arrays
  proc_clim <- df_proc(all_clim)
  
  # Check array
  dim(proc_clim)
  dimnames(proc_clim)
  plot(proc_clim[,,1])
  proc_test <- data.frame(proc_clim[,1,1])
  
  # Get file attributes
  lon <- unique(med_coords$lon)
  lat <- unique(med_coords$lat)[1]
  time <- as.integer(unique(all_clim$t))
  tunits <- "days since 1970-01-01"
  doy <- unique(all_clim$doy)
  eventno <- unique(all_event$event_no)
  
  # Length of each attribute
  nlon <- length(lon)
  nlat <- length(lat)
  ndoy <- length(doy)
  nt <- length(time)
  nen <- length(eventno)
  
  # path and file name, set dname
  ncpath <- "~/Desktop/"
  ncname <- "test"  
  ncfname <- paste0(ncpath, ncname, ".nc")
  dname <- "tmp"
  
  # create and write the netCDF file -- ncdf4 version
  # define dimensions
  londim <- ncdim_def("lon","degrees_east",as.double(lon)) 
  latdim <- ncdim_def("lat","degrees_north",as.double(lat)) 
  timedim <- ncdim_def("time",tunits,as.double(time))
  
  # define variables
  fillvalue <- 1e32
  dlname <- "sea surface temperature"
  tmp_def <- ncvar_def("tmp","deg_C", list(londim, latdim, timedim), fillvalue, dlname, prec = "single")
  
  # create netCDF file and put arrays
  ncout <- nc_create(ncfname, list(tmp_def), force_v4 = TRUE)
  
  # put variables
  ncvar_put(ncout, tmp_def, proc_clim)
  
  # put additional attributes into dimension and data variables
  ncatt_put(ncout,"lon","axis","X") #,verbose=FALSE) #,definemode=FALSE)
  ncatt_put(ncout,"lat","axis","Y")
  ncatt_put(ncout,"time","axis","T")
  
  # add global attributes
  # ncatt_put(ncout,0,"title",title$value)
  # ncatt_put(ncout,0,"institution",institution$value)
  # ncatt_put(ncout,0,"source",datasource$value)
  # ncatt_put(ncout,0,"references",references$value)
  history <- paste("Party Panda", date(), sep=", ")
  # ncatt_put(ncout,0,"history",history)
  # ncatt_put(ncout,0,"Conventions",Conventions$value)
  
  # Get a summary of the created file:
  ncout
}


# Convert files -----------------------------------------------------------



# Test the output ---------------------------------------------------------

test_nc <- tidync("~/Desktop/test.nc") %>% 
  hyper_tibble()


# Full tutorial on dataframe to array -------------------------------------

# Set variable lengths etc.
lon <- as.array(unique(med_lon$lon))
lat <- as.array(unique(med_lat$lat))[1]
time <- as.array(as.integer(unique(all_clim$t)))
tunits <- "days since 1970-01-01"
nlon <- length(lon)
nlat <- length(lat)
nt <- length(time)

# Convert columns to arrays
temp_mat <- as.matrix(all_clim$temp)
temp_array <- array(temp_mat, dim = c(nlon, nlat, nt))
dim(temp_array)

# Check output
cutpts <- seq(10, 30, by = 2)
levelplot(temp_array[,,200] ~ lon * lat, data = med_coords[med_coords$lat == lat,], at = cutpts, cuts = 11, pretty = T, 
          col.regions = (rev(brewer.pal(10,"RdBu"))), main = "Mean July Temperature (C)")

# Convert natural dataframes missing pixels due to land etc.
j2 <- sapply(all_clim$lon, function(x) which.min(abs(lon-x)))
k2 <- sapply(all_clim$lat, function(x) which.min(abs(lat-x)))

# partial loop avoidance for tmp_array
fillvalue <- 1e32
temp_array <- array(fillvalue, dim = c(nlon, nlat))
tmp_array <- array(fillvalue, dim = c(nlon, nlat, nt))
nobs <- dim(all_clim)[1]
for (l in 1:nt) {
  temp_array[cbind(j2,k2)] <- as.matrix(all_clim[1:nobs,l+2]) 
  tmp_array[,,l] <- temp_array
}

# loop avoidance for all of the variables
nobs <- dim(all_clim)[1]
l <- rep(1:nt, each = nobs)
tmp_array[cbind(j2,k2,l)] <- as.matrix(all_clim[1:nobs,3])

