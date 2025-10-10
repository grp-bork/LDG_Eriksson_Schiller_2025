## ====================================================================================================
## R script to build ensemble raster projections for CEPHALOPOD SDMs
## Processes model output folders, constructs ensemble raster stacks for each algorithm
## and bootstrap replicate, and saves the results for downstream analysis.
##
## Author:       Dominic Eriksson
## Date:         8th of October 2025
## Affiliation:  Environmental Physics Group, UP, ETH Zürich, Switzerland
## Contact:      deriksson@ethz.ch
##
## Input files:  Cephalopod model output folders containing CALL.RData, QUERY.RData, and MODEL.RData
##
## Output files: Ensemble raster stacks saved as RData files. File names indicate whether log 
##               or non-log predictors were used and contain information on ensemble members.
##
## Strategy:
##   This script compiles individual model projections from CEPHALOPOD outputs into ensemble raster
##   stacks. For each modeled taxon, it iterates over validated algorithms and bootstrap replicates,
##   converts predicted values into raster layers aligned with environmental grids, and stacks them
##   by month and bootstrap. Negative predictions are set to zero, and spatial reference is applied.
##   Ensemble stacks are saved for each taxon and subfolder to enable downstream analysis of 
##   spatial biodiversity patterns and prediction uncertainty.
##
## Required R packages (tested versions):
##   - abind      1.4.5
##   - doParallel 1.0.17
##   - foreach    1.5.2
##   - raster     3.6.26
##   - dplyr      1.1.2
##   - sp         2.1.4
##   - tidyr      1.3
##
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
library(tidyr)


# Load configuration file from the Cephalopod modeling pipeline (https://github.com/alexschickele/CEPHALOPOD)
setwd("/net/sea/work/deriksson/Projects/CEPHALOPOD")
source("./code/00_config.R")

# Choose log or non-log predictors
trans <- c("Non_log", "Log")
t <- 2

# Depth levels
depth_level <- c("MLD", "0_10Meters")
d <- 1

# Input and output directories --> Save the cephalopod output data in Data_generated folderproject_wd <- "Data/Data_generated/Cephalopod_output"
project_wd <- "Data/Data_generated/Cephalopod_output"
wd_out <- "Code/2_Extract_model_output/2_Output/"
# Create output directory if it doesn't exist
if (!dir.exists(wd_out)) {
  dir.create(wd_out, recursive = TRUE)
}

# List top-level folders to loop over - pattern needs to match the output file name of the cepaholpod pipeline. This can be individual depending on how 
# you set up the modeling procedure. Normaly the folder name contains the name of the target variable (e.g. "name_of_organism_modeled")
FOLDER_NAMES <- list.files(file.path(project_wd, trans[t]), full.names = TRUE, pattern = "chunk")

# Parallelization setup
n.cores <- 20
cl <- makeCluster(n.cores, outfile = "")
registerDoParallel(cl)

# Loop over folders
for (f in seq_along(FOLDER_NAMES)) {
  SUBFOLDER_NAMES <- list.files(FOLDER_NAMES[f], full.names = TRUE, pattern = "1000|5000")

  foreach::foreach(s = seq_along(SUBFOLDER_NAMES), .packages = c("abind", "sp", "raster", "dplyr")) %dopar% {
    SUBFOLDER_NAME <- SUBFOLDER_NAMES[s]
    print(paste("Processing folder", f, "of", length(FOLDER_NAMES), "subfolder", s, "of", length(SUBFOLDER_NAMES)))

    CALL <- get(load(file.path(FOLDER_NAMES[f], "CALL.RData")))
    QUERY <- get(load(file.path(SUBFOLDER_NAME, "QUERY.RData")))

    model_file_path <- file.path(SUBFOLDER_NAME, "MODEL.RData")
    if (!file.exists(model_file_path)) {
      message("MODEL.RData not found in ", SUBFOLDER_NAME, ". Skipping.")
      return(NULL)
    } else {
      MODEL <- get(load(model_file_path))
    }

    if ((length(MODEL$MODEL_LIST) == 0) & CALL$FAST == TRUE) {
      message("No validated algorithms for projections.")
      return(NULL)
    }

    # Initialize raster template
    r0 <- CALL$ENV_DATA[[1]][[1]]
    r0_values <- getValues(r0)
    land <- setValues(r0, ifelse(is.na(r0_values), 9999, NA))

    # Filter recommendations
    rec <- qc_recommandations(QUERY = QUERY, MODEL = MODEL, DATA_TYPE = CALL$DATA_TYPE)
    MODEL$recommandations <- rec
    filtered_rows <- rec %>% filter(PRE_VIP == 1, FIT == 1, CUM_VIP == 1, DEV == 1)
    if (nrow(filtered_rows) == 0) return(NULL)
    loop_over <- rownames(filtered_rows)

    # Build ensemble projections
    list_array <- lapply(loop_over, function(i) {
      if (CALL$DATA_TYPE == "proportions") {
        MODEL[["MBTR"]][["proj"]][["y_hat"]][,,i,]
      } else {
        MODEL[[i]][["proj"]][["y_hat"]]
      }
    })
    names(list_array) <- loop_over

    list_algorithms <- list()
    for (loop in loop_over) {
      d <- list_array[[loop]]
      raster_layers <- list()

      for (bootstrap in 1:dim(d)[2]) {
        for (month in 1:dim(d)[3]) {
          layer_data <- d[, bootstrap, month]
          raster_layer <- setValues(r0, layer_data)
          names(raster_layer) <- paste0("Bootstrap_", bootstrap, "_Month_", month)
          raster_layers[[names(raster_layer)]] <- raster_layer
        }
      }

      r <- raster::stack(raster_layers)
      r[r < 0] <- 0
      raster::projection(r) <- sp::CRS("+proj=longlat +datum=WGS84 +no_defs")
      list_algorithms[[loop]] <- r
    }

    # Save ensemble raster stack
    save(list_algorithms, file = file.path(wd_out, paste0(basename(SUBFOLDER_NAME), "_ensembleMembers_", trans[t], "_v2.RData")))
    print(paste("Finished folder", f, "subfolder", s))
  }
}

# Stop cluster
parallel::stopCluster(cl)

## ====================================================================================================
## END OF SCRIPT
## ====================================================================================================
