
## ====================================================================================================
## This R script extracts individual ensemble members from CEPHALOPOD habitat model
## outputs for prokaryotic taxa. It processes model predictions from parallel-computed chunks,
## extracting bootstrap replicates and temporal predictions to create comprehensive raster stacks
## for downstream uncertainty and variability analyses.
##
## Author:       Dominic Eriksson
## Date:         26th of February 2026
## Affiliation:  Environmental Physics Group UP, ETH Zürich, Switzerland
## Contact:      deriksson@ethz.ch
##
## Input files:  
##   - MODEL.RData, CALL.RData, and QUERY.RData files from chunked CEPHALOPOD outputs
##   - Environmental data layers used as spatial templates
##
## Output files: 
##   - RData files containing raster stacks with individual ensemble members
##   - Each file contains bootstrap replicates across temporal dimensions
##   - Spatial projections formatted for uncertainty quantification
##
## Strategy:
##   The script processes CEPHALOPOD model outputs in parallel, extracting individual ensemble
##   members (bootstrap replicates) for each algorithm and temporal step. Quality control filters
##   ensure only validated model outputs are processed. The resulting raster stacks preserve
##   the full uncertainty structure needed for robust biodiversity assessments and enable
##   detailed analysis of model variability and prediction confidence.
##
## Required R packages (tested versions):
##   - abind      1.4.5
##   - doParallel 1.0.17
##   - foreach    1.5.2
##   - raster     3.6.26
##   - dplyr      1.1.2
##   - sp         1.6.1
##   - tidyr      1.3.0
## ====================================================================================================

# Clear working space
rm(list = ls())

### Load Required Libraries
library(abind)      # For combining multidimensional arrays
library(doParallel) # For parallel processing
library(foreach)    # For parallel loops
library(raster)     # For raster data manipulation
library(dplyr)      # For data manipulation
library(sp)         # For spatial data handling
library(tidyr)      # For data tidying

# Load configuration file from the Cephalopod modeling pipeline (https://github.com/alexschickele/CEPHALOPOD)
source("../../CEPHALOPOD/code/00_config.R")
# Directories
project_wd <- "Code/6_Coastal_influence/Cephalopod_output"
wd_out <- "Code/6_Coastal_influence/2_Output/"
# Create output directory if it doesn't exist
if (!dir.exists(wd_out)) {
  dir.create(wd_out, recursive = TRUE)
}

# Get folder names to loop over
FOLDER_NAMES <- list.files(project_wd, full.names = TRUE, pattern = "chunk")

### Parallelization Setup
n.cores <- 10
cl <- makeCluster(n.cores, outfile = "")
registerDoParallel(cl)

### Process Each Folder and Subfolder
for (f in seq_along(FOLDER_NAMES)) {
# for (f in 5:length(FOLDER_NAMES)) {
  
  # Get subfolder names
  SUBFOLDER_NAMES <- list.files(FOLDER_NAMES[f], full.names = TRUE, pattern = "1000|5000")
  
  # Parallel loop over subfolders
  foreach::foreach(s = seq_along(SUBFOLDER_NAMES), .packages = c("abind", "sp", "raster", "dplyr")) %dopar% {
    
    # Print progress
    print(paste("Processing folder", f, "of", length(FOLDER_NAMES), "and subfolder", s, "of", length(SUBFOLDER_NAMES)))
    
    # Load parameter and query data
    SUBFOLDER_NAME <- SUBFOLDER_NAMES[s]
    CALL <- get(load(paste0(FOLDER_NAMES[f], "/CALL.RData")))
    QUERY <- get(load(paste0(SUBFOLDER_NAME, "/QUERY.RData")))
    
    # Check for MODEL.RData and load if exists
    model_file_path <- paste0(SUBFOLDER_NAME, "/MODEL.RData")
    if (!file.exists(model_file_path)) {
      message("MODEL.RData not found in ", SUBFOLDER_NAME, ". Skipping.")
      return(NULL)
    } else {
      MODEL <- get(load(model_file_path))
    }
    
    # Check if projections exist
    if ((length(MODEL$MODEL_LIST) == 0) & CALL$FAST == TRUE) {
      message("No validated algorithms for projections.")
      return(NULL)
    }
    
    ### Set Initial Parameters and Masks
    r0 <- CALL$ENV_DATA[[1]][[1]]
    # land <- raster::setValues(r0, ifelse(is.na(r0), 9999, NA)) # Not working anymore
    # Extract values from the raster
    r0_values <- getValues(r0)
    # Check for NA values and replace as desired
    new_values <- ifelse(is.na(r0_values), 9999, NA)
    # Set the modified values back to the raster
    land <- setValues(r0, new_values)
    
    ### Filter Recommendations
    rec <- qc_recommandations(QUERY = QUERY, MODEL = MODEL, DATA_TYPE = CALL$DATA_TYPE)
    MODEL$recommandations <- rec
    filtered_rows <- rec %>% filter(PRE_VIP == 1, FIT == 1, CUM_VIP == 1, DEV == 1)
    if (nrow(filtered_rows) == 0) {
      message("No rows meet criteria. Skipping.")
      return(NULL)
    }
    loop_over <- rownames(filtered_rows)
    
    ### Build Ensembles
    list_array <- lapply(loop_over, function(i) {
      if (CALL$DATA_TYPE == "proportions") {
        MODEL[["MBTR"]][["proj"]][["y_hat"]][,,i,]
      } else {
        MODEL[[i]][["proj"]][["y_hat"]]
      }
    })
    
   # Name list_array
   names(list_array) <- loop_over
   
   list_algorithms <- list()
    for( loop in loop_over){   

        # Extract model data
        d <- list_array[[loop]]

            # Initialize a list to store raster layers for each bootstrap-month combination
            raster_layers <- list()

            # Loop through each bootstrap
            for (bootstrap in 1:dim(d)[2]) {  # Loop over the 10 bootstraps
              # Loop through each month
              for (month in 1:dim(d)[3]) {  # Loop over the 12 months
                  # Extract the data for the current bootstrap and month
                  layer_data <- d[, bootstrap, month]
                  
                  # Create a raster layer using the template r0
                  raster_layer <- setValues(r0, layer_data)
                  
                  # Assign a meaningful name to the raster layer
                  layer_name <- paste0("Bootstrap_", bootstrap, "_Month_", month)
                  names(raster_layer) <- layer_name
                  
                  # Store the raster layer in the list
                  raster_layers[[layer_name]] <- raster_layer
              }
            }

        # Stack all in on
        r <- raster::stack(raster_layers)

        # Replace all negative values with zero in the raster stack
        r[r < 0] <- 0

        # Set the projection to WGS84
        wgs84 <- sp::CRS("+proj=longlat +datum=WGS84 +no_defs")
        raster::projection(r) <- wgs84
        
        # Save in the list
        list_algorithms[[loop]] <- r
    }

# Save list object 
save(
    list_algorithms,   
    file = paste0(wd_out, basename(SUBFOLDER_NAME), "_ensembleMembers_v2.RData")
)
    
    # Print progress
    print(paste("Finished folder", f, "subfolder", s))
  }
}

### Stop Cluster
parallel::stopCluster(cl)

## ====================================================================================================
## END OF SCRIPT
## ====================================================================================================