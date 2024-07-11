# RCMEMS
# The purpose of this script is to download data from the CMEMS portal
# This data portal is designed to work with python
# So getting it to run with R is a bit of a hack, but it is supported


# Setup -------------------------------------------------------------------

# This script does not call any packages explicitly
# But does require that some be installed, as may be seen below

# Install Copernicus marine toolbox directly via python via the terminal
## NB: This requires that python and pip are installed locally
# system("python -m pip install copernicusmarine")

# Load credentials
## NB: This file is not uploaded to GitHub to avoid security breach
## Create account here: https://marine.copernicus.eu/
# And then uncomment and run this code:
# CMEMS_cred <- data.frame(username = "username", password = "password")
# readr::write_csv(CMEMS_cred, "metadata/CMEMS_cred.csv")
CMEMS_cred <- readr::read_csv("metadata/CMEMS_cred.csv")

# NB: This stopped working on RStudio, but works on Positron...
# I assume it's something to do with not being able to fins the python path correctly in RStudio


# Download ----------------------------------------------------------------

# Wrapper to pass to CLI
subset_CMEMS <- function(df, data_ID, var_ID, out_dir){
  start_date <- df$start_date; end_date <- df$end_date
  out_name <- paste0(data_ID,"_",start_date,"_", end_date)
  if(!file.exists(paste0(out_dir, out_name))){
    system(paste("copernicusmarine subset -i", data_ID,
    "--username", CMEMS_cred$username, "--password", CMEMS_cred$password,
    "--start-datetime", start_date, "--end-datetime", end_date, 
    "-v", var_ID, "-o", out_dir, "-f", out_name, "--force-download"))
  }

}

# Dataframe of start and end dates to cycle through
dates <- data.frame(start_date = paste0(1982:2023, "-01-01"), 
                    end_date = paste0(1982:2023, "-12-31"),
                    row_idx = seq_len(length(1982:2023)))

# Ply the downloads by year
## Takes about 20 seconds per year
plyr::d_ply(dates, c("row_idx"), subset_CMEMS, .parallel = FALSE,
            data_ID = "cmems_SST_MED_SST_L4_REP_OBSERVATIONS_010_021", 
            var_ID = "analysed_sst", out_dir = "~/data/MED_REP/")


# Test --------------------------------------------------------------------

ncdf4::nc_open("~/data/MED_REP/cmems_SST_MED_SST_L4_REP_OBSERVATIONS_010_021_2023-01-01_2023-12-31.nc")

sst_test <- tidync::tidync("~/data/MED_REP/cmems_SST_MED_SST_L4_REP_OBSERVATIONS_010_021_2023-01-01_2023-12-31.nc") |> 
  tidync::hyper_tibble()

