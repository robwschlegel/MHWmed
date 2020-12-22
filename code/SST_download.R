# SST_download.R
# This script contains the code to FTP access the hi-res SST product for the Med


# Setup -------------------------------------------------------------------

# Packages not available via CRAN
remotes::install_github("skgrange/threadr")

# The packages we will use
library(tidyverse) # A staple for most modern data management in R
# library(RCurl) # For helping R to make sense of URLs for web hosted data
# library(XML) # For reading the HTML tables created by RCurl
library(tidync) # For easily dealing with NetCDF data
# library(doParallel) # For parallel processing
library(threadr) # For downloading from FTP sites that require user credentials
# library(RCMEMS) # For subsetting CMEMS data before download


# Download ----------------------------------------------------------------

# ftp://nrt.cmems-du.eu/Core/

# As the file was not currently mounted Nathaniel sent it to me