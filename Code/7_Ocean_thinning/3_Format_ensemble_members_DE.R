## ====================================================================================================
## This R script formats individual ensemble member files from basin-stratified subsampled CEPHALOPOD
## model outputs into structured data tables. Due to the large size of datasets from 1,100 sampling
## scenarios (11 cap values × 100 iterations), it processes files individually using parallel computing
## to avoid memory constraints while preserving the full ensemble structure for sampling bias analysis.
##
## Author:       Dominic Eriksson
## Date:         26th of February 2026
## Affiliation:  Environmental Physics Group UP, ETH Zürich, Switzerland
## Contact:      deriksson@ethz.ch
##
## Input files:  
##   - Individual ensemble member RData files from basin-stratified subsampling model runs
##   - Raster stacks containing bootstrap replicates across sampling scenarios and temporal dimensions
##
## Output files: 
##   - Individual CSV files for each sampling scenario with ensemble members in long format
##   - Files structured for memory-efficient processing of sampling bias sensitivity analysis
##   - Optimized format for calculating ensemble statistics across different sample size caps
##
## Strategy:
##   The script employs parallel processing to handle the computationally intensive task of converting
##   raster-based ensemble outputs to tabular format across 1,100 different sampling scenarios.
##   Each file is processed independently and saved immediately to prevent memory overflow while
##   preserving the complete uncertainty structure needed for robust sampling bias assessment.
##   This approach enables efficient downstream analysis of how basin-specific sampling intensity
##   affects biodiversity prediction reliability.
##
## Required R packages (tested versions):
##   - raster     3.6.26
##   - tidyr      1.3.0
##   - dplyr      1.1.2
##   - data.table 1.14.8
##   - parallel   4.2.2
## ====================================================================================================

# Clear workspace
rm(list = ls())

# Libraries
library(raster)     # For raster data manipulation
library(tidyr)      # For data tidying
library(dplyr)      # For data manipulation
library(data.table) # For efficient data handling
library(parallel)   # For parallel processing

# Directories
wd_in <- "Code/7_Ocean_thinning/2_Output" # for loading ensemble member files from basin-stratified subsampling
wd_out <- "Code/7_Ocean_thinning/3_Output/"
# Create output directory if it does not exist
if (!dir.exists(wd_out)) {
    dir.create(wd_out, recursive = TRUE)
}


# Get filenames
fnames <- list.files(wd_in, full.names = TRUE)

# Function to process each file
process_file <- function(f) {
    # Print progress
    cat(paste0("Loading file ", which(fnames == f), ", out of ", length(fnames), ".\n"))
    cat(paste0("File name: ", f, "\n"))
    
    # Subset
    data <- get(load(f))
    
    # Get index
    idx <- length(data)
    
    # Save in dataframe
    l2 <- list()
    for(i in 1:idx) {
        df <- raster::as.data.frame(data[[i]], xy = TRUE)
        df$alg <- names(data)[i]
        
        # Add clade indicator
        df$clade <- gsub("_ensembleMembers_v2.RData", "", basename(f))
        
        # Store in list object
        l2[[i]] <- df
    }
    
    # Merge into one dataframe
    df <- do.call("rbind", l2)
    
    # Extract clade name for filename
    clade_name <- gsub("_ensembleMembers_v2.RData", "", basename(f))
    
    # Save dataframe directly to avoid memory issues
    data.table::fwrite(
        df, 
        file = paste0(wd_out, "ensembleMembers_monthly_", clade_name, "_v2.csv")
    )
    
    cat(paste0("Saved: ensembleMembers_monthly_", clade_name, "_v2.csv\n"))
}

# Process all ensemble member files using parallel processing to avoid memory issues
# Use parallel processing with 20 cores
mclapply(fnames, process_file, mc.cores = 20)

cat("\nAll ensemble member files have been processed and saved individually.\n")
cat(paste0("Files saved in: ", wd_out, "\n"))

# Note: Individual files are saved separately to avoid memory issues.
# Each file contains ensemble members for one sampling scenario and is named:
# "ensembleMembers_monthly_{scenario_name}_v2.csv"
# 
# If you need to combine scenarios later, read them back individually
# and process them in chunks to avoid memory constraints.

## ====================================================================================================
## END OF SCRIPT
## ====================================================================================================