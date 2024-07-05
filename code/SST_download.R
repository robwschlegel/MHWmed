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

# Set up virtual environment
virtualenv_create(envname = "CopernicusMarine", python = install_python())
virtualenv_install("CopernicusMarine", packages = c("copernicusmarine"))
reticulate::use_virtualenv("CopernicusMarine", required = TRUE)

# Store the python package to use the functions in R
cmt <- import("copernicusmarine")

# Add login crednetials if necessary
# I.e. https://www.copernicus.eu/en
cmt$login("<username>", "<password>")

# Or install directly via python via the terminal
system("python -m pip install copernicusmarine")


# Download ----------------------------------------------------------------

# NB: For whatever reason, this hangs and does not download
cmt$subset(
  dataset_id="cmems_SST_MED_SST_L4_REP_OBSERVATIONS_010_021",
  variables=list("analysed_sst"),
  minimum_longitude=-18.125,
  maximum_longitude=36.32500076293945,
  minimum_latitude=30.125,
  maximum_latitude=46.025001525878906,
  start_datetime="2024-06-04T00:00:00",
  end_datetime="2024-06-04T00:00:00",
  output_directory = "~/pCloudDrive/data/CMEMS_Med/"
)

# Or use CLI directly
# NB: One would create a function or loop to create the CLI calls desired
system('copernicusmarine get -i cmems_SST_MED_SST_L4_REP_OBSERVATIONS_010_021 --filter "*2023012[1]*" -o "/home/robert/pCloudDrive/data/CMEMS_Med/" --force-download')

# Subset via CLI
# NB: Subset can also be used to compile files into a single download
system('copernicusmarine subset -i cmems_SST_MED_SST_L4_REP_OBSERVATIONS_010_021 --start-datetime 2022-01-01 --end-datetime 2022-01-31 -o "/home/robert/pCloudDrive/data/CMEMS_Med/" --force-download')
