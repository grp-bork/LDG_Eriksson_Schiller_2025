## ====================================================================================================
## This R script computes global annual richness and Shannon diversity ensemble means and standard
## deviations for genus-level prokaryotic taxa. It processes ensemble member datasets to generate
## spatially explicit monthly and annual diversity statistics, creates visualization outputs, and
## handles both richness and Shannon diversity indices with parallel processing for efficiency.
##
## Author:       Dominic Eriksson
## Date:         6th of March 2026
## Affiliation:  Environmental Physics Group UP, ETH Zürich, Switzerland
## Contact:      deriksson@ethz.ch
##
## Input files:
##   - Ensemble member files: ./Code/3_Compute_ensembles_and_uncertainties/1_Output/*.csv
##   - Monthly and annual richness data from CEPHALOPOD model outputs
##   - Shannon diversity ensemble member files (optional)
##
## Output files:
##   - Ensemble means (monthly/annual): ./Code/3_Compute_ensembles_and_uncertainties/2_Output/ensemble_means_*.csv
##   - Visualization plots: ./Code/3_Compute_ensembles_and_uncertainties/2_Output/ensemble_mean_annual_*.png
##   - Shannon diversity outputs: ./Code/3_Compute_ensembles_and_uncertainties/2_Output/ensemble_means_*_shannon_*.csv
##
## Strategy:
##   This script performs parallel processing of ensemble datasets to compute diversity statistics.
##   For each taxonomic clade, it calculates ensemble means and standard deviations across algorithms
##   and bootstrap replicates for all 12 months plus annual means. The workflow includes: (1) parallel
##   processing of richness ensemble files with mean and SD calculations for monthly and annual data;
##   (2) generation of spatial visualization plots showing annual ensemble means; (3) similar processing
##   pipeline for Shannon diversity indices; and (4) comprehensive error handling and progress reporting
##   with detailed summaries of successful and failed processing attempts.
##
## Required R packages (tested versions):
##   - data.table    (version 1.15.4)
##   - dplyr         (version 1.1.2)
##   - ggplot2       (version 3.4.2)
##   - parallel      (version 4.3.2)
##   - foreach       (version 1.5.2)
##   - doParallel    (version 1.0.17)
##
## ====================================================================================================

# Clear working space
rm(list = ls())

### ===========================
### Libraries
### ===========================
library(data.table)  # For fast data manipulation and file I/O
library(dplyr)       # For data manipulation and grouping operations
library(ggplot2)     # For spatial visualization and plotting
library(parallel)    # For parallel computing infrastructure
library(foreach)     # For parallel foreach loops
library(doParallel)  # For parallel processing backend

### ===========================
### Configurable Project Paths - adjust as needed
### ===========================
# Input and output directories
wd_in <- "./Code/3_Compute_ensembles_and_uncertainties/1_Output/"
wd_out <- "./Code/3_Compute_ensembles_and_uncertainties/2_Output/"

# Create output directory if it does not exist
if (!dir.exists(wd_out)) {
    dir.create(wd_out, recursive = TRUE)
}

## ====================================================================================================
## SECTION 1: Species Richness Ensemble Processing
## ====================================================================================================

# Setup parallel processing
n_cores <- min(detectCores() - 1, 20)  # Use available cores minus 1, max 20
cat("Using", n_cores, "cores for parallel processing\n")
cat("Total files to process:", length(fnames), "\n")

cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Export necessary variables and functions to cluster
clusterExport(cl, c("fnames", "wd_out"))
clusterEvalQ(cl, {
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

# Define function to process one file
process_clade_file <- function(fname_path) {
  # Load data
  df <- data.table::fread(fname_path)
  
  # Compute ensemble means for all monthly columns and annual_mean
  # Average across alg and bootstrap for each x, y, clade combination
  df_ensemble_means <- df[
    , .(
      mean_month_1 = mean(month_1, na.rm = TRUE),
      mean_month_2 = mean(month_2, na.rm = TRUE),
      mean_month_3 = mean(month_3, na.rm = TRUE),
      mean_month_4 = mean(month_4, na.rm = TRUE),
      mean_month_5 = mean(month_5, na.rm = TRUE),
      mean_month_6 = mean(month_6, na.rm = TRUE),
      mean_month_7 = mean(month_7, na.rm = TRUE),
      mean_month_8 = mean(month_8, na.rm = TRUE),
      mean_month_9 = mean(month_9, na.rm = TRUE),
      mean_month_10 = mean(month_10, na.rm = TRUE),
      mean_month_11 = mean(month_11, na.rm = TRUE),
      mean_month_12 = mean(month_12, na.rm = TRUE),
      mean_annual = mean(annual_mean, na.rm = TRUE),
      
      # Also compute standard deviations for uncertainty quantification
      sd_month_1 = sd(month_1, na.rm = TRUE),
      sd_month_2 = sd(month_2, na.rm = TRUE),
      sd_month_3 = sd(month_3, na.rm = TRUE),
      sd_month_4 = sd(month_4, na.rm = TRUE),
      sd_month_5 = sd(month_5, na.rm = TRUE),
      sd_month_6 = sd(month_6, na.rm = TRUE),
      sd_month_7 = sd(month_7, na.rm = TRUE),
      sd_month_8 = sd(month_8, na.rm = TRUE),
      sd_month_9 = sd(month_9, na.rm = TRUE),
      sd_month_10 = sd(month_10, na.rm = TRUE),
      sd_month_11 = sd(month_11, na.rm = TRUE),
      sd_month_12 = sd(month_12, na.rm = TRUE),
      sd_annual = sd(annual_mean, na.rm = TRUE)
    ),
    by = .(x, y, clade)
  ]
  
  # Get clade name for file saving
  current_clade <- unique(df_ensemble_means$clade)
  safe_clade_name <- gsub("[^A-Za-z0-9_-]", "_", current_clade)
  
  # Save the ensemble means data
  output_file <- file.path(wd_out, paste0("ensemble_means_monthly_annual_", safe_clade_name, ".csv"))
  data.table::fwrite(df_ensemble_means, output_file)
  
  # Create and save visualization
  tryCatch({
    p <- ggplot(df_ensemble_means) +
      geom_raster(aes(x = x, y = y, fill = mean_annual)) +
      scale_fill_gradient2(
        low = "blue",
        mid = "green", 
        high = "red",
        midpoint = median(df_ensemble_means$mean_annual, na.rm = TRUE),
        na.value = "grey80"
      ) +
      labs(
        title = paste("Annual Ensemble Mean for", current_clade),
        x = "Longitude",
        y = "Latitude", 
        fill = "Mean Richness"
      ) +
      theme_minimal()
    
    plot_file <- file.path(wd_out, paste0("ensemble_mean_annual_", safe_clade_name, ".png"))
    ggsave(plot_file, plot = p, width = 10, height = 6)
  }, error = function(e) {
    # If plot fails, continue without stopping
    warning(paste("Plot failed for", current_clade, ":", e$message))
  })
  
  # Return summary info
  list(
    clade = current_clade,
    file = basename(fname_path),
    rows = nrow(df_ensemble_means),
    output_file = output_file,
    mean_annual_range = range(df_ensemble_means$mean_annual, na.rm = TRUE)
  )
}

# Export the function to cluster
clusterExport(cl, "process_clade_file")

# Process all files in parallel
cat("Starting parallel processing...\n")
start_time <- Sys.time()

results <- foreach(i = seq_along(fnames), 
                   .packages = c("data.table", "dplyr", "ggplot2"),
                   .errorhandling = "pass") %dopar% {
  
  fname_path <- fnames[i]
  
  # Process the file
  result <- process_clade_file(fname_path)
  
  # Add progress info
  result$file_index <- i
  result$total_files <- length(fnames)
  
  return(result)
}

# Stop cluster
stopCluster(cl)

end_time <- Sys.time()
cat("Parallel processing completed in", round(as.numeric(end_time - start_time, units = "mins"), 2), "minutes\n")

### Processing Results Summary and Output Validation

# Summarize results
successful_results <- results[sapply(results, function(x) !inherits(x, "try-error"))]
failed_results <- results[sapply(results, function(x) inherits(x, "try-error"))]

cat("\n=== PROCESSING SUMMARY ===\n")
cat("Total files processed:", length(fnames), "\n")
cat("Successful:", length(successful_results), "\n")
cat("Failed:", length(failed_results), "\n")

if (length(failed_results) > 0) {
  cat("\nFailed files:\n")
  for (i in which(sapply(results, function(x) inherits(x, "try-error")))) {
    cat(" -", basename(fnames[i]), ":", as.character(results[[i]]), "\n")
  }
}

# Check output files
actual_csv_files <- list.files(wd_out, pattern = "ensemble_means_monthly_annual_.*\\.csv", full.names = FALSE)
actual_png_files <- list.files(wd_out, pattern = "ensemble_mean_annual_.*\\.png", full.names = FALSE)

cat("\nOutput files created:\n")
cat("CSV files:", length(actual_csv_files), "\n")
cat("PNG files:", length(actual_png_files), "\n")

# Sample of successful results
if (length(successful_results) > 0) {
  cat("\nSample of processed clades:\n")
  sample_results <- head(successful_results, 5)
  for (result in sample_results) {
    cat(" -", result$clade, ": ", result$rows, "rows, annual range:", 
        round(result$mean_annual_range[1], 3), "-", round(result$mean_annual_range[2], 3), "\n")
  }
}

cat("\nAll ensemble means saved to:", wd_out, "\n")

## ====================================================================================================
## SECTION 2: Shannon Diversity Ensemble Processing
## ====================================================================================================

### Shannon Diversity File Processing

# Get Shannon diversity filenames to load
fnames_shannon <- list.files(
  path = wd_in,
  pattern = "monthly_and_annual_shannon_ensemble_members_.*\\.csv",
  full.names = TRUE
)

cat("\n=== STARTING SHANNON DIVERSITY PROCESSING ===\n")
cat("Shannon files found:", length(fnames_shannon), "\n")

if (length(fnames_shannon) > 0) {
  
  # Setup parallel processing for Shannon
  cl_shannon <- makeCluster(n_cores)
  registerDoParallel(cl_shannon)
  
  # Export necessary variables and functions to cluster
  clusterExport(cl_shannon, c("fnames_shannon", "wd_out"))
  clusterEvalQ(cl_shannon, {
    library(data.table)
    library(dplyr)
    library(ggplot2)
  })
  
  # Define function to process one Shannon file
  process_clade_file_shannon <- function(fname_path) {
    # Load data
    df <- data.table::fread(fname_path)
    
    # Compute ensemble means for all monthly columns and annual_mean
    # Average across alg and bootstrap for each x, y, clade combination
    df_ensemble_means <- df[
      , .(
        mean_month_1 = mean(month_1, na.rm = TRUE),
        mean_month_2 = mean(month_2, na.rm = TRUE),
        mean_month_3 = mean(month_3, na.rm = TRUE),
        mean_month_4 = mean(month_4, na.rm = TRUE),
        mean_month_5 = mean(month_5, na.rm = TRUE),
        mean_month_6 = mean(month_6, na.rm = TRUE),
        mean_month_7 = mean(month_7, na.rm = TRUE),
        mean_month_8 = mean(month_8, na.rm = TRUE),
        mean_month_9 = mean(month_9, na.rm = TRUE),
        mean_month_10 = mean(month_10, na.rm = TRUE),
        mean_month_11 = mean(month_11, na.rm = TRUE),
        mean_month_12 = mean(month_12, na.rm = TRUE),
        mean_annual = mean(annual_mean, na.rm = TRUE),
        
        # Also compute standard deviations for uncertainty quantification
        sd_month_1 = sd(month_1, na.rm = TRUE),
        sd_month_2 = sd(month_2, na.rm = TRUE),
        sd_month_3 = sd(month_3, na.rm = TRUE),
        sd_month_4 = sd(month_4, na.rm = TRUE),
        sd_month_5 = sd(month_5, na.rm = TRUE),
        sd_month_6 = sd(month_6, na.rm = TRUE),
        sd_month_7 = sd(month_7, na.rm = TRUE),
        sd_month_8 = sd(month_8, na.rm = TRUE),
        sd_month_9 = sd(month_9, na.rm = TRUE),
        sd_month_10 = sd(month_10, na.rm = TRUE),
        sd_month_11 = sd(month_11, na.rm = TRUE),
        sd_month_12 = sd(month_12, na.rm = TRUE),
        sd_annual = sd(annual_mean, na.rm = TRUE)
      ),
      by = .(x, y, clade)
    ]
    
    # Get clade name for file saving
    current_clade <- unique(df_ensemble_means$clade)
    safe_clade_name <- gsub("[^A-Za-z0-9_-]", "_", current_clade)
    
    # Save the ensemble means data
    output_file <- file.path(wd_out, paste0("ensemble_means_monthly_annual_shannon_", safe_clade_name, ".csv"))
    data.table::fwrite(df_ensemble_means, output_file)
    
    # Create and save visualization
    tryCatch({
      p <- ggplot(df_ensemble_means) +
        geom_raster(aes(x = x, y = y, fill = mean_annual)) +
        scale_fill_gradient2(
          low = "blue",
          mid = "green", 
          high = "red",
          midpoint = median(df_ensemble_means$mean_annual, na.rm = TRUE),
          na.value = "grey80"
        ) +
        labs(
          title = paste("Annual Shannon Diversity Ensemble Mean for", current_clade),
          x = "Longitude",
          y = "Latitude", 
          fill = "Mean Shannon\nDiversity"
        ) +
        theme_minimal()
      
      plot_file <- file.path(wd_out, paste0("ensemble_mean_annual_shannon_", safe_clade_name, ".png"))
      ggsave(plot_file, plot = p, width = 10, height = 6)
    }, error = function(e) {
      # If plot fails, continue without stopping
      warning(paste("Shannon plot failed for", current_clade, ":", e$message))
    })
    
    # Return summary info
    list(
      clade = current_clade,
      file = basename(fname_path),
      rows = nrow(df_ensemble_means),
      output_file = output_file,
      mean_annual_range = range(df_ensemble_means$mean_annual, na.rm = TRUE)
    )
  }
  
  # Export the function to cluster
  clusterExport(cl_shannon, "process_clade_file_shannon")
  
  # Process all Shannon files in parallel
  cat("Starting Shannon parallel processing...\n")
  start_time_shannon <- Sys.time()
  
  results_shannon <- foreach(i = seq_along(fnames_shannon), 
                     .packages = c("data.table", "dplyr", "ggplot2"),
                     .errorhandling = "pass") %dopar% {
    
    fname_path <- fnames_shannon[i]
    
    # Process the file
    result <- process_clade_file_shannon(fname_path)
    
    # Add progress info
    result$file_index <- i
    result$total_files <- length(fnames_shannon)
    
    return(result)
  }
  
  # Stop cluster
  stopCluster(cl_shannon)
  
  end_time_shannon <- Sys.time()
  cat("Shannon parallel processing completed in", round(as.numeric(end_time_shannon - start_time_shannon, units = "mins"), 2), "minutes\n")
  
  # Summarize Shannon results
  successful_results_shannon <- results_shannon[sapply(results_shannon, function(x) !inherits(x, "try-error"))]
  failed_results_shannon <- results_shannon[sapply(results_shannon, function(x) inherits(x, "try-error"))]
  
  cat("\n=== SHANNON PROCESSING SUMMARY ===\n")
  cat("Total Shannon files processed:", length(fnames_shannon), "\n")
  cat("Successful:", length(successful_results_shannon), "\n")
  cat("Failed:", length(failed_results_shannon), "\n")
  
  if (length(failed_results_shannon) > 0) {
    cat("\nFailed Shannon files:\n")
    for (i in which(sapply(results_shannon, function(x) inherits(x, "try-error")))) {
      cat(" -", basename(fnames_shannon[i]), ":", as.character(results_shannon[[i]]), "\n")
    }
  }
  
  # Check Shannon output files
  actual_csv_files_shannon <- list.files(wd_out, pattern = "ensemble_means_monthly_annual_shannon_.*\\.csv", full.names = FALSE)
  actual_png_files_shannon <- list.files(wd_out, pattern = "ensemble_mean_annual_shannon_.*\\.png", full.names = FALSE)
  
  cat("\nShannon output files created:\n")
  cat("CSV files:", length(actual_csv_files_shannon), "\n")
  cat("PNG files:", length(actual_png_files_shannon), "\n")
  
  # Sample of successful Shannon results
  if (length(successful_results_shannon) > 0) {
    cat("\nSample of processed Shannon clades:\n")
    sample_results_shannon <- head(successful_results_shannon, 5)
    for (result in sample_results_shannon) {
      cat(" -", result$clade, ": ", result$rows, "rows, annual range:", 
          round(result$mean_annual_range[1], 3), "-", round(result$mean_annual_range[2], 3), "\n")
    }
  }
  
  cat("\nAll Shannon ensemble means saved to:", wd_out, "\n")
  
} else {
  cat("No Shannon diversity files found in:", wd_in, "\n")
  cat("Looking for files matching pattern: monthly_and_annual_shannon_ensemble_members_*.csv\n")
}

## ====================================================================================================
## END OF SCRIPT
## ====================================================================================================