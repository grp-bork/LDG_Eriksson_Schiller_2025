## ====================================================================================================
## R script to calculate monthly and annual mean and standard deviation of species richness 
## projections from CEPHALOPOD SDMs. Processes multiple model output folders using parallel 
## computation and saves ensemble raster stacks for downstream analyses.
##
## Author:       Dominic Eriksson
## Date:         8th of October 2025
## Affiliation:  Environmental Physics Group, UP, ETH Zürich, Switzerland
## Contact:      deriksson@ethz.ch
##
## Input files:  Cephalopod model output folders containing CALL.RData, QUERY.RData, and MODEL.RData
##
## Output files: Raster stacks of ensemble monthly and annual mean and standard deviation, 
##               including predictive power and number of algorithms contributing to the ensemble.
##
## Strategy:
##   This script compiles model projections for each taxon, algorithm, and bootstrap replicate
##   into multi-layer raster stacks. Monthly predictions are combined into mean and standard
##   deviation rasters, and annual statistics are calculated across all months. Ensemble
##   predictions are augmented with predictive performance metrics and the number of algorithms
##   included. Parallel processing is employed to efficiently process multiple taxa and 
##   subfolders in the CEPHALOPOD output.
##
## Required R packages (tested versions):
##   - abind      1.4.5
##   - doParallel 1.0.17
##   - foreach    1.5.2
##   - raster     3.6.26
##   - dplyr      1.1.2
##   - sp         2.1.4
##   - tidyr      1.3.1
## ====================================================================================================


# Clear workspace
rm(list = ls())

# Load required libraries
library(abind)
library(doParallel)
library(foreach)
library(raster)
library(dplyr)
library(sp)

# Use here::here() for portable paths
if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
library(here)

# Source configuration file (relative path) (https://github.com/alexschickele/CEPHALOPOD)
setwd("/net/sea/work/deriksson/Projects/CEPHALOPOD") # adjust to your working directory
source("./code/00_config.R")

# Log and Non-log model outputs
trans <- c("Non_log", "Log")
t <- 1

# Depth levels
depth_level <- c("MLD", "0_10Meters") # we only use MLD in our manuscript, d = 1
d <- 1

# Input and output directories --> Save the cephalopod output data in Data_generated folder
project_wd <- "Data/Data_generated/Cephalopod_output/"
wd_out <- paste0("Code/2_Extract_model_output/3_Output/", trans[t])
if (!dir.exists(wd_out)) dir.create(wd_out, recursive = TRUE)

# List top-level folders
FOLDER_NAMES <- list.files(file.path(project_wd, trans[t]), full.names = TRUE, pattern = "chunk")

# Loop over folders
for (f in seq_along(FOLDER_NAMES)) {
  SUBFOLDER_NAMES <- list.files(FOLDER_NAMES[f], full.names = TRUE, pattern = "1000|5000")

  # Parallelization setup
  n.cores <- 10
  cl <- makeCluster(n.cores, outfile = "")
  registerDoParallel(cl)

  foreach::foreach(s = seq_along(SUBFOLDER_NAMES), .packages = c("abind", "sp", "raster", "dplyr")) %dopar% {
    SUBFOLDER_NAME <- SUBFOLDER_NAMES[s]
    print(paste("Processing folder", f, "subfolder", s))

    CALL <- get(load(file.path(FOLDER_NAMES[f], "CALL.RData")))
    QUERY <- get(load(file.path(SUBFOLDER_NAME, "QUERY.RData")))

    model_file_path <- file.path(SUBFOLDER_NAME, "MODEL.RData")
    if (!file.exists(model_file_path)) {
      message("MODEL.RData not found in ", SUBFOLDER_NAME, ". Skipping.")
      return(NULL)
    } else {
      MODEL <- get(load(model_file_path))
    }

    if ((length(MODEL$MODEL_LIST) == 0) & CALL$FAST == TRUE) return(NULL)

    # Raster template
    r0 <- CALL$ENV_DATA[[1]][[1]]
    r0_values <- getValues(r0)
    land <- setValues(r0, ifelse(is.na(r0_values), 9999, NA))

    # Filter recommendations
    rec <- qc_recommandations(QUERY = QUERY, MODEL = MODEL, DATA_TYPE = CALL$DATA_TYPE)
    MODEL$recommandations <- rec
    filtered_rows <- rec %>% filter(PRE_VIP == 1, FIT == 1, CUM_VIP == 1, DEV == 1)
    if (nrow(filtered_rows) == 0) return(NULL)
    loop_over <- rownames(filtered_rows)

    # Build ensemble arrays
    list_array <- lapply(loop_over, function(i) {
      if (CALL$DATA_TYPE == "proportions") {
        MODEL[["MBTR"]][["proj"]][["y_hat"]][,,i,]
      } else {
        MODEL[[i]][["proj"]][["y_hat"]]
      }
    })

    combined_array <- abind::abind(list_array, along = 4)
    combined_array[combined_array < 0] <- 0

    mean_array <- apply(combined_array, c(1, 3), mean, na.rm = TRUE)
    sd_array <- apply(combined_array, c(1, 3), sd, na.rm = TRUE)

    r_mean_stack <- raster::stack(lapply(1:12, function(i) setValues(r0, mean_array[, i])))
    r_sd_stack <- raster::stack(lapply(1:12, function(i) setValues(r0, sd_array[, i])))
    names(r_mean_stack) <- names(r_sd_stack) <- month.name

    annual_mean_array <- apply(combined_array, c(1, 2, 4), mean, na.rm = TRUE)
    r_annual_mean <- setValues(r0, apply(annual_mean_array, 1, mean, na.rm = TRUE))
    r_annual_sd <- setValues(r0, apply(annual_mean_array, 1, sd, na.rm = TRUE))
    names(r_annual_mean) <- "Annual_Mean"
    names(r_annual_sd) <- "Annual_Standard_Deviation"

    r_ensemble_means <- raster::stack(r_mean_stack, r_annual_mean)
    r_ensemble_sds <- raster::stack(r_sd_stack, r_annual_sd)

    wgs84 <- sp::CRS("+proj=longlat +datum=WGS84 +no_defs")
    raster::projection(r_ensemble_means) <- wgs84
    raster::projection(r_ensemble_sds) <- wgs84

    predictive_power <- ifelse(length(MODEL$MODEL_LIST) == 1, MODEL[[MODEL$MODEL_LIST]]$eval$R2, MODEL$ENSEMBLE$eval$R2)
    n_algorithms <- length(loop_over)

    r_ensemble_means <- raster::addLayer(r_ensemble_means, setValues(r_annual_mean, rep(predictive_power, ncell(r0))), setValues(r_annual_mean, rep(n_algorithms, ncell(r0))))
    names(r_ensemble_means)[(nlayers(r_ensemble_means)-2):nlayers(r_ensemble_means)] <- c("Annual_Mean", "Predictive_Power", "Number_of_Algorithms")

    # Save results
    l_final <- list(r_ensemble_means, r_ensemble_sds)
    save(l_final, file = file.path(wd_out, paste0(basename(SUBFOLDER_NAME), "_v2.RData")))
    print(paste("Finished folder", f, "subfolder", s))
  }

  parallel::stopCluster(cl)
}

## ====================================================================================================
## END OF SCRIPT
## ====================================================================================================

