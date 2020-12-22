# code/MHW_detect.R
# This script detects MHWs in the temperature data


# Setup -------------------------------------------------------------------

# Necessary packages
library(tidyverse)
library(heatwaveR)
library(tidync)

# The data locations
med_SST_files <- dir("data/", pattern = ".nc", full.names = T)


# Load data ---------------------------------------------------------------

sst_test <- tidync(med_SST_files[1])
