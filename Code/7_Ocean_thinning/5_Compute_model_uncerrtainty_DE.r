
## ====================================================================================================
## This R script computes model uncertainty metrics from basin-stratified subsampling sensitivity
## analysis by calculating coefficient of variation in prokaryotic richness predictions across
## different sampling scenarios. It quantifies how variations in sample size caps per ocean basin
## affect the reliability and consistency of biodiversity predictions at each spatial location.
##
## Author:       Dominic Eriksson
## Date:         26th of February 2026
## Affiliation:  Environmental Physics Group UP, ETH Zürich, Switzerland
## Contact:      deriksson@ethz.ch
##
## Input files:  
##   - Annual richness CSV files from 1,100 basin-stratified subsampling scenarios
##   - Each file contains ensemble mean predictions for one sampling configuration
##   - Files represent different sample size caps (10-350 samples per basin) across iterations
##
## Output files: 
##   - CSV file with coefficient of variation metrics for each spatial coordinate
##   - Includes mean richness, standard deviation, and CV across sampling scenarios
##   - Spatial uncertainty map for assessing sampling bias effects on predictions
##
## Strategy:
##   The script employs memory-efficient batch processing to handle 1,100 scenarios simultaneously.
##   For each spatial coordinate, it calculates the coefficient of variation in richness predictions
##   across all sampling scenarios, providing a robust measure of how sensitive biodiversity
##   predictions are to sampling bias. High CV values indicate locations where sampling
##   heterogeneity strongly affects model predictions, while low CV values suggest robust
##   predictions regardless of sampling intensity variations.
##
## Required R packages (tested versions):
##   - data.table 1.14.8
## ====================================================================================================

# Clear workspace
rm(list = ls())

# Libraries
library(data.table) # For efficient data handling

# Directory
wd_in <- "Code/7_Ocean_thinning/4_Output/"
wd_out <- "Code/7_Ocean_thinning/5_Output/"

# Create output directory if it does not exist
if (!dir.exists(wd_out)) {
    dir.create(wd_out, recursive = TRUE)
}

# Get filenames
fnames <- list.files(wd_in, full.names = TRUE)

# Function to calculate coefficient of variation across sampling scenarios for each coordinate
calculate_richness_cv <- function(file_list, wd_out, output_file = "richness_cv_by_coordinates.csv", 
                                  batch_size = 50, verbose = TRUE) {
  
  if(verbose) cat("Processing", length(file_list), "files in batches of", batch_size, "\n")
  
  # Initialize storage for all coordinates and their richness values
  all_data <- list()
  
  # Process files in batches to manage memory
  n_batches <- ceiling(length(file_list) / batch_size)
  
  for(batch in 1:n_batches) {
    if(verbose) cat("Processing batch", batch, "of", n_batches, "\n")
    
    # Determine file indices for this batch
    start_idx <- (batch - 1) * batch_size + 1
    end_idx <- min(batch * batch_size, length(file_list))
    batch_files <- file_list[start_idx:end_idx]
    
    # Process files in current batch
    batch_data <- rbindlist(lapply(batch_files, function(file) {
      if(verbose && (which(batch_files == file) %% 10 == 1)) {
        cat("  Reading file", which(batch_files == file), "of", length(batch_files), "\n")
      }
      
      tryCatch({
        # Read the file
        df <- fread(file)
        
        # Ensure we have the expected columns
        if(!all(c("x", "y", "annual_richness") %in% names(df))) {
          warning("Missing expected columns in file: ", file)
          return(data.table())
        }
        
        # Remove rows with NA richness values
        df <- df[!is.na(annual_richness)]
        
        # Keep only necessary columns
        df[, .(x, y, annual_richness)]
        
      }, error = function(e) {
        warning("Error reading file ", file, ": ", e$message)
        return(data.table())
      })
    }))
    
    # Combine with existing data
    if(nrow(batch_data) > 0) {
      all_data[[batch]] <- batch_data
    }
    
    # Clean up memory
    gc()
  }
  
  if(verbose) cat("Combining all batches...\n")
  
  # Combine all batches
  combined_data <- rbindlist(all_data)
  
  if(verbose) cat("Calculating coefficient of variation by coordinates...\n")
  
  # Calculate CV for each coordinate pair
  cv_results <- combined_data[, {
    richness_values <- annual_richness[!is.na(annual_richness)]
    if(length(richness_values) >= 2) {
      mean_val <- mean(richness_values)
      sd_val <- sd(richness_values)
      cv <- if(mean_val != 0) (sd_val / mean_val) * 100 else NA
      list(
        n_scenarios = length(richness_values),
        mean_richness = mean_val,
        sd_richness = sd_val,
        cv_richness = cv
      )
    } else {
      list(
        n_scenarios = length(richness_values),
        mean_richness = ifelse(length(richness_values) == 1, richness_values[1], NA),
        sd_richness = NA,
        cv_richness = NA
      )
    }
  }, by = .(x, y)]
  
  # Sort by coordinates for easier interpretation
  setorder(cv_results, x, y)
  
  if(verbose) {
    cat("Results summary:\n")
    cat("  Total coordinate pairs:", nrow(cv_results), "\n")
    cat("  Pairs with valid CV:", sum(!is.na(cv_results$cv_richness)), "\n")
    cat("  Mean number of scenarios per coordinate:", mean(cv_results$n_scenarios), "\n")
  }
  
  # Save results
  fwrite(cv_results, file.path(wd_out, output_file))
  if(verbose) cat("Results saved to:", output_file, "\n")
  
  return(cv_results)
}

# Execute coefficient of variation calculation across all sampling scenarios
results <- calculate_richness_cv(fnames, wd_out, output_file = "prokaryote_richness_cv.csv")

## ====================================================================================================
## END OF SCRIPT
## ====================================================================================================
