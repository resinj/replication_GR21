# Download and load data used in data examples
# Warning: Some files will be saved in the current working directory!
library(openxlsx)

# Bank of England CPI inflation forecasts
temp = tempfile()
download.file("https://www.bankofengland.co.uk/-/media/boe/files/monetary-policy-report/2021/february/mpr-february-2021-chart-slides-and-data.zip",temp)
unzip(zipfile=temp, files = "Parameters for MPC CPI inflation projections from February 2004.xlsx")
raw_fcasts = read.xlsx("Parameters for MPC CPI inflation projections from February 2004.xlsx",sheet = "CPI Forecast")
unlink(temp)

# ONS CPI inflation rates
temp = tempfile()
download.file("https://www.ons.gov.uk/generator?format=csv&uri=/economy/inflationandpriceindices/timeseries/d7g7/mm23",temp)
raw_obs = read.csv(temp,header = FALSE)
unlink(temp)

# Tredennick et al. (2020) butterfly population forecasts and observations
temp = tempfile()
download.file("https://zenodo.org/record/4311358/files/modsel-guide-v1.zip",temp)
unzip(zipfile=temp, files = "modsel-guide-v1/results/obs_vs_pred.csv",junkpaths = TRUE)
obs_pred <- read.csv("obs_vs_pred.csv",header=T)
unlink(temp)

