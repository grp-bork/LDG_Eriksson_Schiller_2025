##====================================================================================================
## Ensemble Means with Predictive Power Extraction for Genus-Level Prokaryotic Taxa
##
## This R script calculates monthly and annual mean values and standard deviations of species
## richness projections from CEPHALOPOD models for genus-level prokaryotic taxa analysis.
## The script processes ensemble predictions across multiple algorithms and bootstraps, computes
## latitudinal diversity gradient uncertainties, and merges results for comprehensive analysis.
## Features parallel processing for computational efficiency and includes predictive power metrics.
##
## Author:      Dominic Eriksson
##              Environmental Physics Group, UP
##              ETH Zürich, Switzerland
## Contact:     deriksson@ethz.ch
## Date:        February 27th, 2026
## Affiliation: ETH Zürich, Environmental Physics Group, UP
##
## Input files:
## - ./Code/Cephalopod_output/Non_log_genusLevel/chunks*/*/MODEL.RData (Genus-level model outputs)
## - ./Code/Cephalopod_output/Non_log/chunks*/*/MODEL.RData (Class-level model outputs for uncertainties)
## - ./Code/Cephalopod_output/*/CALL.RData (Model configuration files)
## - ./Code/Cephalopod_output/*/*/QUERY.RData (Query parameter files)
##
## Output files:
## - ./Code/2_Extract_model_outputs/6_Output/GenusLevel/*.RData (Monthly/annual ensemble means and SDs)
## - ./Code/2_Extract_model_outputs/6_Output/Non_log/LDG_uncertainties/*_monthly_v3.csv (Monthly LDG uncertainties)
## - ./Code/2_Extract_model_outputs/6_Output/Non_log/LDG_uncertainties/*_annual_v3.csv (Annual LDG uncertainties)
## - Merged_LDG_uncertainties_monthly_v3.csv (Consolidated monthly uncertainties)
## - Merged_LDG_uncertainties_annual_v2.csv (Consolidated annual uncertainties)
##
## Strategy:
## 1. Process genus-level ensemble members with parallel computing for efficiency
## 2. Calculate monthly richness means and standard deviations across algorithms
## 3. Compute annual richness from monthly means before ensemble averaging
## 4. Extract predictive power (R²) and algorithm counts as metadata layers
## 5. Generate latitudinal diversity gradient uncertainty quantification
## 6. Process 4D arrays (space × bootstrap × month × algorithm) into summary statistics
## 7. Merge class-level results for comprehensive uncertainty analysis
## 8. Export results in standardized formats for downstream visualization
##
## Required R packages:
## - abind (1.4.5): Multidimensional array operations and binding along specified dimensions
## - doParallel (1.0.17): Parallel processing backend for foreach loops with cluster management
## - foreach (1.5.2): Enhanced looping constructs supporting parallel execution
## - raster (3.5.15): Spatial raster data manipulation, stacking, and projection handling
## - dplyr (1.0.9): Data manipulation with grouping, summarizing, and filtering operations
## - sp (1.4.7): Spatial data classes and coordinate reference system definitions
## - parallel (4.2.2): Core parallel computing functionality for cluster operations
##====================================================================================================

# Clear working space
rm(list = ls())

### Load Required Libraries
library(abind)      # For combining multidimensional arrays
library(doParallel) # For parallel processing
library(foreach)    # For parallel loops
library(raster)     # For raster data manipulation
library(dplyr)      # For data manipulation
library(sp)         # For spatial data handling

# Load specific functions
setwd("../../CEPHALOPOD")
source(file = "../../CEPHALOPOD/code/00_config.R")

## Log and Non-log model outputs ---------------------------------------------------------------

# # Choose between log and non-log transformed predictors
# trans <- c("Non_log", "Log")
# t <- 1

# ### Set Parameters and Directories
# depth_level <- c("MLD", "0_10Meters")  # Define depth levels: Mixed Layer Depth (MLD) or 0-10 meters
# d <- 1                                 # Specify the depth level index

# Directories
project_wd <- "./Code/Cephalopod_output/Non_log_genusLevel"
wd_out <- paste0("./Code/2_Extract_model_outputs/6_Output/GenusLevel/")
# Create output directory if it does not exist
if (!dir.exists(wd_out)) {
  dir.create(wd_out, recursive = TRUE)
}

# Get folder names to loop over
FOLDER_NAMES <- list.files(project_wd, full.names = TRUE, pattern = "chunk")


### Process Each Folder and Subfolder
# for (f in seq_along(FOLDER_NAMES)) {
for (f in 18:length(FOLDER_NAMES)) {
  
  # Get subfolder names
  SUBFOLDER_NAMES <- list.files(FOLDER_NAMES[f], full.names = TRUE, pattern = "1000|5000")
  
  ### Parallelization Setup
  n.cores <- 10
  cl <- makeCluster(n.cores, outfile = "")
  registerDoParallel(cl)

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
    
    ### Monthly Ensembles
    combined_array <- abind::abind(list_array, along = 4)
    # Set all negative values to zero
    combined_array[combined_array < 0] <- 0

    mean_array <- apply(combined_array, c(1, 3), mean, na.rm = TRUE)
    sd_array <- apply(combined_array, c(1, 3), sd, na.rm = TRUE)
    
    r_mean_stack <- raster::stack(lapply(1:12, function(i) setValues(r0, mean_array[, i])))
    r_sd_stack <- raster::stack(lapply(1:12, function(i) setValues(r0, sd_array[, i])))
    names(r_mean_stack) <- names(r_sd_stack) <- month.name
    
    ### Annual Ensembles
    annual_mean_array <- apply(combined_array, c(1, 2, 4), mean, na.rm = TRUE)
    r_annual_mean <- setValues(r0, apply(annual_mean_array, c(1), mean, na.rm = TRUE))
    r_annual_sd <- setValues(r0, apply(annual_mean_array, c(1), sd, na.rm = TRUE))
    names(r_annual_mean) <- "Annual_Mean"
    names(r_annual_sd) <- "Annual_Standard_Deviation"
    
    ### Merge Results
    r_ensemble_means <- raster::stack(r_mean_stack, r_annual_mean)
    r_ensemble_sds <- raster::stack(r_sd_stack, r_annual_sd)
    
    ### Add Metadata
    wgs84 <- sp::CRS("+proj=longlat +datum=WGS84 +no_defs")
    raster::projection(r_ensemble_means) <- wgs84
    raster::projection(r_ensemble_sds) <- wgs84
    predictive_power <- ifelse(length(MODEL$MODEL_LIST) == 1, 
                                MODEL[[MODEL$MODEL_LIST]]$eval$R2, 
                                MODEL$ENSEMBLE$eval$R2)
    n_algorithms <- length(loop_over)
    
    r_ensemble_means <- raster::addLayer(r_ensemble_means, setValues(r_annual_mean, rep(predictive_power, ncell(r0))), setValues(r_annual_mean, rep(n_algorithms, ncell(r0))))
    names(r_ensemble_means)[nlayers(r_ensemble_means) - 2] <- "Annual_Mean"
    names(r_ensemble_means)[nlayers(r_ensemble_means) - 1] <- "Predictive_Power"
    names(r_ensemble_means)[nlayers(r_ensemble_means)] <- "Number_of_Algorithms"
    
    ### Save Results
    l_final <- list(r_ensemble_means, r_ensemble_sds)
    save(l_final, file = paste0(wd_out, basename(SUBFOLDER_NAME), ".RData"))
    print(paste("Finished folder", f, "subfolder", s))
  }

  ### Stop Cluster
  parallel::stopCluster(cl)

} # Close loop across folder names


## ====================================================================================================
## LDG uncertainties
## ====================================================================================================


# Load libraries
library(abind)
library(sp)
library(raster)
library(dplyr)
library(parallel)
library(doParallel)

# Choose between log and non-log transformed predictors
trans <- c("Non_log", "Log")
t <- 1

### Set Parameters and Directories
depth_level <- c("MLD", "0_10Meters")  # Define depth levels: Mixed Layer Depth (MLD) or 0-10 meters
d <- 1                                 # Specify the depth level index

# Output directory
wd_out <- "./Code/2_Extract_model_outputs/6_Output/"


# Directories
project_wd <- "./Code/Cephalopod_output"

# Get folder names to loop over
FOLDER_NAMES <- list.files(paste0(project_wd, "/", trans[t]), full.names = TRUE, pattern = "chunk")

### Compute uncertainty of the annual LDG

### Parallelization Setup
n.cores <- 20
cl <- makeCluster(n.cores, outfile = "")
registerDoParallel(cl)

### Process Each Folder and Subfolder
for (f in seq_along(FOLDER_NAMES)) {
# for (f in 1:11) {
  
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
    
    ### Check Dimensions Before Abind
    dims <- sapply(list_array, dim)
    print(dims)  # This will show the dimensions of each array in the list_array

    # If there is a dimension mismatch, skip the subfolder
    if (any(apply(dims, 1, function(x) length(unique(x))) > 1)) {
      message("Dimension mismatch detected in ", SUBFOLDER_NAME, ". Skipping.")
      return(NULL)
    }

    ### Monthly Ensembles
    combined_array <- tryCatch({
      abind::abind(list_array, along = 4)
    }, error = function(e) {
      message("Error in abind for ", SUBFOLDER_NAME, ": ", e$message)
      return(NULL)
    })

    # Replace negative values with zero
    combined_array[combined_array < 0] <- 0
    
    # If abind fails, skip the subfolder
    if (is.null(combined_array)) {
      return(NULL)
    }
    
    # Extract latitudes from raster dimensions - Be careful to check if the sequence is in the right order
    # latitudes <- seq(-89.5, 89.5, by = 1)  # Assuming 180 rows from -90 to 90 with 1° resolution
    latitudes <- seq(89.5, -89.5, by = -1) 

    # Convert the 4D array into a data frame
    # df <- data.frame(
    #   lat = rep(latitudes, each = 360),  # Repeat each latitude across 360 longitudes
    #   value = as.vector(combined_array), # Flatten the array to a vector
    #   bootstrap = rep(1:10, each = 64800 * 12 * 2),
    #   month = rep(rep(1:12, each = 64800 * 2), times = 10),
    #   algorithm = rep(rep(1:length(loop_over), each = 64800), times = 10 * 12)
    # )

      df <- data.frame(
        lat = rep(latitudes, each = 360),  # Repeat each latitude across 360 longitudes
        value = as.vector(combined_array), # Flatten the array to a vector
        bootstrap = rep(1:10, each = 64800 * 12 * length(loop_over)),
        month = rep(rep(1:12, each = 64800 * length(loop_over)), times = 10),
        algorithm = rep(rep(1:length(loop_over), each = 64800), times = 10 * 12)
    )

    # Compute latitudinal mean per bootstrap, month, and algorithm
    lat_summary <- df %>%
      group_by(lat, bootstrap, month, algorithm) %>%
      summarise(lat_mean = mean(value, na.rm = TRUE), .groups = "drop")

    # Now summarize across bootstrap and algorithm
    final_summary <- lat_summary %>%
      group_by(lat, month) %>%
      summarise(
        min = min(lat_mean, na.rm = TRUE),
        max = max(lat_mean, na.rm = TRUE),
        mean = mean(lat_mean, na.rm = TRUE),
        median = median(lat_mean, na.rm = TRUE),
        sd = sd(lat_mean, na.rm = TRUE),
        iqr = IQR(lat_mean, na.rm = TRUE),
        .groups = "drop"
      )

    # Step 1: Compute the annual mean for each latitude, bootstrap, and algorithm
    annual_lat_summary <- lat_summary %>%
      group_by(lat, bootstrap, algorithm) %>%
      summarise(annual_mean = mean(lat_mean, na.rm = TRUE), .groups = "drop")

    # Step 2: Compute summary statistics across bootstrap and algorithm
    annual_final_summary <- annual_lat_summary %>%
      group_by(lat) %>%
      summarise(
        min = min(annual_mean, na.rm = TRUE),
        max = max(annual_mean, na.rm = TRUE),
        mean = mean(annual_mean, na.rm = TRUE),
        median = median(annual_mean, na.rm = TRUE),
        sd = sd(annual_mean, na.rm = TRUE),
        iqr = IQR(annual_mean, na.rm = TRUE),
        .groups = "drop"
      )

    # Write results to CSV
    write.csv(final_summary, file = paste0(wd_out, trans[t], "/LDG_uncertainties/", basename(SUBFOLDER_NAME), "_statsAcrossEnsembleMembers_LDG_monthly_v3.csv"))
    write.csv(annual_final_summary, file = paste0(wd_out, trans[t], "/LDG_uncertainties/", basename(SUBFOLDER_NAME), "_statsAcrossEnsembleMembers_LDG_annual_v3.csv"))
    
    print(paste("Finished folder", f, "subfolder", s))
  }
}

### Stop Cluster
parallel::stopCluster(cl)

### ==============================================================================================
### Merge for Jonas
### ===========================================================================================


## Monthly ---------------------------------------------

# Get filenames
fnames <- list.files("/net/sea/work/deriksson/Projects/Prokaryotes_LDG_v3/Code/2_Extract_model_outputs/6_Output/Non_log/LDG_uncertainties", pattern = "_v3.csv", full.names = TRUE)

# Grep richness 5000 files
fnames <- fnames[grep("Richness_5000", fnames)]
fnames <- fnames[grep("monthly", fnames)]

## Add Desulfobacterota since it is not in the non-log files
# Check for directories containing "Desulfobacterota"
base_dir <- "./Code/2_Extract_model_outputs/6_Output/"
# List all directories and subdirectories
all_dirs <- list.dirs(base_dir, recursive = TRUE, full.names = TRUE)
# Filter for those containing "Desulfobacterota"
desulfo_dirs <- all_dirs[grep("Desulfobacterota", all_dirs)]
# Display results
if (length(desulfo_dirs) > 0) {
  cat("Found directories containing 'Desulfobacterota':\n")
  print(desulfo_dirs)
} else {
  cat("No directories containing 'Desulfobacterota' found.\n")
}

# Subset for class levels

## Define the class-level clades - Dependent on the input data 
# Option 1:
class_level_clades <- c(
  "d__Archaea_p__Thermoproteota_c__Nitrososphaeria",
  "d__Bacteria_p__Fibrobacterota_c__Fibrobacteria",
  "d__Bacteria_p__Bacteroidota_c__Bacteroidia",
  "d__Archaea_p__Thermoplasmatota_c__Poseidoniia",
  "d__Bacteria_p__Verrucomicrobiota_c__Verrucomicrobiae",
  "d__Bacteria_p__Desulfobacterota_D_c__UBA1144",
  "d__Bacteria_p__Marinisomatota_c__Marinisomatia",
  "d__Bacteria_p__Chloroflexota_c__Dehalococcoidia",
  "d__Bacteria_p__Pseudomonadota_c__Gammaproteobacteria",
  "d__Bacteria_p__SAR324_c__SAR324",
  "d__Bacteria_p__Chlamydiota_c__Chlamydiia",
  "d__Bacteria_p__Bacteroidota_c__Rhodothermia",
  "d__Bacteria_p__Verrucomicrobiota_c__Kiritimatiellia",
  "d__Bacteria_p__Actinomycetota_c__Acidimicrobiia",
  "d__Bacteria_p__Cyanobacteriota_c__Cyanobacteriia",
  "d__Bacteria_p__Pseudomonadota_c__Alphaproteobacteria",
  "d__Archaea_p__Asgardarchaeota_c__Heimdallarchaeia",
  "d__Bacteria_p__Bdellovibrionota_c__Bacteriovoracia",
  "d__Bacteria_p__Myxococcota_c__UBA796",
  "d__Bacteria_p__Myxococcota_c__UBA4151",
  "d__Bacteria_p__Myxococcota_c__XYA12-FULL-58-9"
)
# Option 2:


# Filter the filenames that contain any of the class-level clades
filtered_fnames <- fnames[grep(paste(class_level_clades, collapse = "|"), fnames)]

# View the filtered filenames
print(filtered_fnames)

# Loop
l_sd <- list()
for( f in seq_along(filtered_fnames) ){

  # Print progress
  print(paste("Processing file", f, "of", length(filtered_fnames)))

  # Open data
  df <- read.csv(filtered_fnames[f], header = TRUE, sep = ",", stringsAsFactors = FALSE)

  # Remove first column
  df <- df[, -1]

  # Format clade name
  clade <- basename(filtered_fnames[f])
  clade <- gsub("_statsAcrossEnsembleMembers_LDG_monthly.csv", "", clade)
  clade <- gsub("Richness_5000_", "", clade)

  df$clade <- clade

  # Save in list
  l_sd[[f]] <- df

}

# Merge all dataframes in the list into one dataframe
merged_df <- do.call(rbind, l_sd)

# Write the merged dataframe to a CSV file
write.csv(merged_df, file = "./Code/2_Extract_model_outputs/6_Output/Non_log/LDG_uncertainties/Merged_LDG_uncertainties_monthly_v3.csv", row.names = FALSE)




## Annual ---------------------------------------------

# Get filenames
fnames <- list.files("./Code/2_Extract_model_outputs/6_Output/Non_log/LDG_uncertainties", pattern = "_v3.csv", full.names = TRUE)

# Grep richness 5000 files
fnames <- fnames[grep("Richness_5000", fnames)]
fnames <- fnames[grep("annual", fnames)]

# Subset for class levels


## Define the class-level clades - Dependent on the input data
# Option 1:
class_level_clades <- c(
  "d__Archaea_p__Thermoproteota_c__Nitrososphaeria",
  "d__Bacteria_p__Fibrobacterota_c__Fibrobacteria",
  "d__Bacteria_p__Bacteroidota_c__Bacteroidia",
  "d__Archaea_p__Thermoplasmatota_c__Poseidoniia",
  "d__Bacteria_p__Verrucomicrobiota_c__Verrucomicrobiae",
  "d__Bacteria_p__Desulfobacterota_D_c__UBA1144",
  "d__Bacteria_p__Marinisomatota_c__Marinisomatia",
  "d__Bacteria_p__Chloroflexota_c__Dehalococcoidia",
  "d__Bacteria_p__Pseudomonadota_c__Gammaproteobacteria",
  "d__Bacteria_p__SAR324_c__SAR324",
  "d__Bacteria_p__Chlamydiota_c__Chlamydiia",
  "d__Bacteria_p__Bacteroidota_c__Rhodothermia",
  "d__Bacteria_p__Verrucomicrobiota_c__Kiritimatiellia",
  "d__Bacteria_p__Actinomycetota_c__Acidimicrobiia",
  "d__Bacteria_p__Cyanobacteriota_c__Cyanobacteriia",
  "d__Bacteria_p__Pseudomonadota_c__Alphaproteobacteria",
  "d__Archaea_p__Asgardarchaeota_c__Heimdallarchaeia",
  "d__Bacteria_p__Bdellovibrionota_c__Bacteriovoracia",
  "d__Bacteria_p__Myxococcota_c__UBA796",
  "d__Bacteria_p__Myxococcota_c__UBA4151",
  "d__Bacteria_p__Myxococcota_c__XYA12-FULL-58-9"
)
# Extract just the class names from the full taxonomic paths
class_names_only <- gsub(".*_c__(.+)", "\\1", class_level_clades)

# Filter the filenames that contain any of the class names
filtered_fnames <- fnames[grep(paste(class_names_only, collapse = "|"), fnames)]

# View the filtered filenames
print(filtered_fnames)

# Loop
l_sd <- list()
for( f in seq_along(filtered_fnames) ){

  # Print progress
  print(paste("Processing file", f, "of", length(filtered_fnames)))

  # Open data
  df <- read.csv(filtered_fnames[f], header = TRUE, sep = ",", stringsAsFactors = FALSE)

  # Remove first column
  df <- df[, -1]

  # Format clade name
  clade <- basename(filtered_fnames[f])
  clade <- gsub("_statsAcrossEnsembleMembers_LDG_annual.csv", "", clade)
  clade <- gsub("Richness_5000_", "", clade)

  df$clade <- clade

  # Save in list
  l_sd[[f]] <- df

}

# Merge all dataframes in the list into one dataframe
merged_df <- do.call(rbind, l_sd)

# Write the merged dataframe to a CSV file
write.csv(merged_df, file = "./Code/2_Extract_model_outputs/6_Output/Non_log/LDG_uncertainties/Merged_LDG_uncertainties_annual_v2.csv", row.names = FALSE)


## ====================================================================================================
## Percentage increase - NOT NEEDED
## ====================================================================================================


# Load libraries
library(abind)
library(sp)
library(raster)
library(dplyr)
library(parallel)
library(doParallel)

# Choose between log and non-log transformed predictors
trans <- c("Non_log", "Log")
t <- 1

### Set Parameters and Directories
depth_level <- c("MLD", "0_10Meters")  # Define depth levels: Mixed Layer Depth (MLD) or 0-10 meters
d <- 1                                 # Specify the depth level index

# Output directory
wd_out <- "./Code/2_Extract_model_outputs/6_Output/"

### Compute uncertainty of the annual LDG

### Parallelization Setup
n.cores <- 20
cl <- makeCluster(n.cores, outfile = "")
registerDoParallel(cl)

### Process Each Folder and Subfolder
for (f in seq_along(FOLDER_NAMES)) {
# for (f in 1:11) {
  
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
    
    ### Check Dimensions Before Abind
    dims <- sapply(list_array, dim)
    print(dims)  # This will show the dimensions of each array in the list_array

    # If there is a dimension mismatch, skip the subfolder
    if (any(apply(dims, 1, function(x) length(unique(x))) > 1)) {
      message("Dimension mismatch detected in ", SUBFOLDER_NAME, ". Skipping.")
      return(NULL)
    }

    ### Monthly Ensembles
    combined_array <- tryCatch({
      abind::abind(list_array, along = 4)
    }, error = function(e) {
      message("Error in abind for ", SUBFOLDER_NAME, ": ", e$message)
      return(NULL)
    })
    
    # If abind fails, skip the subfolder
    if (is.null(combined_array)) {
      return(NULL)
    }
    
    # Extract latitudes from raster dimensions - Be careful to check if the sequence is in the right order
    # latitudes <- seq(-89.5, 89.5, by = 1)  # Assuming 180 rows from -90 to 90 with 1° resolution
    latitudes <- seq(89.5, -89.5, by = -1) 

    # Convert the 4D array into a data frame
    # df <- data.frame(
    #   lat = rep(latitudes, each = 360),  # Repeat each latitude across 360 longitudes
    #   value = as.vector(combined_array), # Flatten the array to a vector
    #   bootstrap = rep(1:10, each = 64800 * 12 * 2),
    #   month = rep(rep(1:12, each = 64800 * 2), times = 10),
    #   algorithm = rep(rep(1:length(loop_over), each = 64800), times = 10 * 12)
    # )

        df <- data.frame(
      lat = rep(latitudes, each = 360),  # Repeat each latitude across 360 longitudes
      value = as.vector(combined_array), # Flatten the array to a vector
      bootstrap = rep(1:10, each = 64800 * 12 * length(loop_over)),
      month = rep(rep(1:12, each = 64800 * length(loop_over)), times = 10),
      algorithm = rep(rep(1:length(loop_over), each = 64800), times = 10 * 12)
    )

    # Compute latitudinal mean per bootstrap, month, and algorithm
    lat_summary <- df %>%
      group_by(lat, bootstrap, month, algorithm) %>%
      summarise(lat_mean = mean(value, na.rm = TRUE), .groups = "drop")

    # Now summarize across bootstrap and algorithm
    final_summary <- lat_summary %>%
      group_by(lat, month) %>%
      summarise(
        min = min(lat_mean, na.rm = TRUE),
        max = max(lat_mean, na.rm = TRUE),
        mean = mean(lat_mean, na.rm = TRUE),
        median = median(lat_mean, na.rm = TRUE),
        sd = sd(lat_mean, na.rm = TRUE),
        iqr = IQR(lat_mean, na.rm = TRUE),
        .groups = "drop"
      )

    # Step 1: Compute the annual mean for each latitude, bootstrap, and algorithm
    annual_lat_summary <- lat_summary %>%
      group_by(lat, bootstrap, algorithm) %>%
      summarise(annual_mean = mean(lat_mean, na.rm = TRUE), .groups = "drop")

    # Step 2: Compute summary statistics across bootstrap and algorithm
    annual_final_summary <- annual_lat_summary %>%
      group_by(lat) %>%
      summarise(
        min = min(annual_mean, na.rm = TRUE),
        max = max(annual_mean, na.rm = TRUE),
        mean = mean(annual_mean, na.rm = TRUE),
        median = median(annual_mean, na.rm = TRUE),
        sd = sd(annual_mean, na.rm = TRUE),
        iqr = IQR(annual_mean, na.rm = TRUE),
        .groups = "drop"
      )

    # Write results to CSV
    write.csv(final_summary, file = paste0(wd_out, trans[t], "/", depth_level[d], "/LDG_uncertainties/", basename(SUBFOLDER_NAME), "_statsAcrossEnsembleMembers_LDG_monthly.csv"))
    write.csv(annual_final_summary, file = paste0(wd_out, trans[t], "/", depth_level[d], "/LDG_uncertainties/", basename(SUBFOLDER_NAME), "_statsAcrossEnsembleMembers_LDG_annual.csv"))
    
    print(paste("Finished folder", f, "subfolder", s))
  }
}

### Stop Cluster
parallel::stopCluster(cl)

## ====================================================================================================
## END OF SCRIPT
## ====================================================================================================
