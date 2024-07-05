# RCMEMS
# The purpose of this script is to download data from the CMEMS portal
# This data portal is designed to work with python
# So gettting it to run with R is a bit of a hack, but it is supported


# Setup -------------------------------------------------------------------

# Get the development version of reticulate
remotes::install_github('rstudio/reticulate')

# The packages we will use
library(tidyverse) # A staple for most modern data management in R
library(tidync) # For easily dealing with NetCDF data
library(reticulate) # For using Python in R

# NB: Another way to do this, but not recommended. Easier to use CLI
# Set up virtual environment
# virtualenv_create(envname = "CopernicusMarine", python = install_python())
# virtualenv_install("CopernicusMarine", packages = c("copernicusmarine"))
# reticulate::use_virtualenv("CopernicusMarine", required = TRUE)
# Store the python package to use the functions in R
# cmt <- import("copernicusmarine")
# Add login crednetials if necessary
# I.e. https://www.copernicus.eu/en
# cmt$login("<username>", "<password>")

# Install copernicus marine toolbox directly via python via the terminal
system("python -m pip install copernicusmarine")


# Download ----------------------------------------------------------------

# Wrapper to pass to CLI
subset_CMEMS <- function(df, data_ID, var_ID, out_dir){
  start_date <- df$start_date; end_date <- df$end_date
  out_name <- paste0(data_ID,"_",start_date,"_", end_date)
  system(paste("copernicusmarine subset -i", data_ID, "--start-datetime", start_date, "--end-datetime", end_date, 
              "-v", var_ID, "-o", out_dir, "-f", out_name, "--force-download"))
}

# Dataframe of start and end dates to cycle through
dates <- data.frame(start_date = paste0(1982:2023, "-01-01"), end_date = paste0(1982:2023, "-12-31")) |> 
  mutate(row_idx = 1:n())

# Ply the downloads by year
## Takes about 20 seconds per year
plyr::d_ply(dates, c("row_idx"), subset_CMEMS,
            data_ID = "cmems_SST_MED_SST_L4_REP_OBSERVATIONS_010_021", 
            var_ID = "analysed_sst", out_dir = "~/data/CMEMS_Med/")


# Test --------------------------------------------------------------------

ncdf4::nc_open("~/data/CMEMS_Med/cmems_SST_MED_SST_L4_REP_OBSERVATIONS_010_021_1982-01-01_1982-12-31.nc")

sst_test <- tidync::tidync("~/data/CMEMS_Med/cmems_SST_MED_SST_L4_REP_OBSERVATIONS_010_021_1982-01-01_1982-12-31.nc") |> hyper_tibble()
