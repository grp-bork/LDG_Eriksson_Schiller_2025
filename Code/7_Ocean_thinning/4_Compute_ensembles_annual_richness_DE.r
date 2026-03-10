
## ====================================================================================================
## This R script computes annual ensemble mean richness from monthly predictions across basin-
## stratified subsampling scenarios. It processes individual ensemble member files from 1,100
## different sampling scenarios (11 cap values × 100 iterations), aggregating temporal predictions
## and bootstrap replicates to create annual richness estimates for sampling bias sensitivity analysis.
##
## Author:       Dominic Eriksson
## Date:         26th of February 2026
## Affiliation:  Environmental Physics Group UP, ETH Zürich, Switzerland
## Contact:      deriksson@ethz.ch
##
## Input files:  
##   - Individual CSV files with monthly ensemble members from basin-stratified subsampling
##   - Each file contains bootstrap replicates across temporal dimensions for one sampling scenario
##
## Output files: 
##   - Annual richness CSV files for each sampling scenario
##   - Aggregated predictions ready for sampling bias uncertainty quantification
##   - Files structured for comparative analysis across different sample size caps
##
## Strategy:
##   The script employs parallel processing to handle the computationally intensive aggregation
##   across 1,100 sampling scenarios. For each scenario, monthly predictions are averaged to
##   annual values, then bootstrap replicates and algorithms are aggregated to produce final
##   ensemble means. This hierarchical approach preserves the sampling structure needed for
##   robust assessment of how basin-specific sampling intensity affects biodiversity predictions
##   while maintaining computational efficiency through memory-optimized processing.
##
## Required R packages (tested versions):
##   - data.table 1.14.8
##   - tidyr      1.3.0
##   - dplyr      1.1.2
##   - parallel   4.2.2
## ====================================================================================================

# Clear workspace
rm(list = ls())

# Load libraries
library(data.table) # For efficient data handling
library(tidyr)      # For data tidying
library(dplyr)      # For data manipulation
library(parallel)   # For parallel processing

# Directories
wd_in <- "Code/7_Ocean_thinning/3_Output/"
wd_out <- "Code/7_Ocean_thinning/4_Output/"

# Create output directory if it does not exist
if (!dir.exists(wd_out)) {
    dir.create(wd_out, recursive = TRUE)
}

# Get filenames
fnames <- list.files(wd_in, full.names = TRUE)

# Function to process each sampling scenario file
process_annual_richness <- function(f) {
    # Print progress
    cat(paste0("Processing file ", which(fnames == f), " out of ", length(fnames), ": ", basename(f), "\n"))
    
    # Load dataframe
    df_final <- data.table::fread(f)
    
    # Extract scenario name from filename
    clade_name <- gsub("ensembleMembers_monthly_|_v2\\.csv", "", basename(f))
    
    # Reformat table
    df_long <- df_final %>%
      tidyr::pivot_longer(
        cols = starts_with("Bootstrap_"),
        names_to = "bootstrap_month",
        values_to = "richness"
      ) %>%
      tidyr::extract(
        col = bootstrap_month,
        into = c("bootstrap", "month"),
        regex = "Bootstrap_(\\d+)_Month_(\\d+)"
      ) %>%
      mutate(
        bootstrap = as.integer(bootstrap),
        month = as.integer(month)
      )
    
    # Convert back to data.table for efficient aggregation
    df_long <- data.table::as.data.table(df_long)
    
    # Compute annual richness for each ensemble member
    df_agg <- df_long[
      , .(annual_richness = mean(richness, na.rm = TRUE)),
      by = .(x, y, bootstrap, alg, clade)
    ]
    
    # Compute annual means for each distance group
    df_agg <- df_agg[
      , .(annual_richness = mean(annual_richness, na.rm = TRUE)),  
      by = .(x, y, clade)
    ]
    
    # Save result with clade name in filename
    output_file <- paste0(wd_out, "annual_richness_", clade_name, "_v2.csv")
    data.table::fwrite(df_agg, file = output_file)
    
    cat(paste0("Saved: annual_richness_", clade_name, "_v2.csv\n"))
    
    return(NULL)  # Don't return data to save memory
}

# Process all files in parallel using available cores for efficient computation
cat("Starting parallel processing of annual richness calculations...\n")
mclapply(fnames, process_annual_richness, mc.cores = 20)

cat("\nAll files processed successfully!\n")
cat(paste0("Output files saved in: ", wd_out, "\n"))

## ====================================================================================================
## END OF SCRIPT
## ====================================================================================================

