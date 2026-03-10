
## ====================================================================================================
## This R script removes marginal seas from prokaryotic diversity ensemble datasets using IHO ocean
## region shapefiles. It processes both richness and Shannon diversity data by identifying points in
## marginal seas and setting their diversity values to NA, while preserving spatial coordinates and
## metadata for downstream analysis. The script handles large datasets efficiently using parallel
## processing and generates comprehensive summaries and visualizations.
##
## Author:       Dominic Eriksson
## Date:         6th of March 2026
## Affiliation:  Environmental Physics Group UP, ETH Zürich, Switzerland
## Contact:      deriksson@ethz.ch
##
## Input files:
##   - Ensemble means: ./Code/3_Compute_ensembles_and_uncertainties/2_Output/ensemble_means_*.csv
##   - Shannon diversity: ./Code/3_Compute_ensembles_and_uncertainties/2_Output/ensemble_means_*_shannon_*.csv
##   - IHO ocean regions shapefile function: ./Functions/iho_ocean_regions.R
##
## Output files:
##   - Filtered richness: ./Code/3_Compute_ensembles_and_uncertainties/2c_Output/marginal_seas_filtered_*.csv
##   - Filtered Shannon: ./Code/3_Compute_ensembles_and_uncertainties/2c_Output/marginal_seas_filtered_shannon_*.csv
##   - Visualization plots: ./Code/3_Compute_ensembles_and_uncertainties/2c_Output/marginal_seas_filtering_summary*.png
##   - Processing summaries: ./Code/3_Compute_ensembles_and_uncertainties/2c_Output/processing_summary*.csv
##
## Strategy:
##   This script performs spatial filtering of marine diversity data to exclude marginal seas from
##   global ocean analyses. The workflow includes: (1) loading ensemble datasets and applying IHO
##   ocean region classification to identify marginal seas (Mediterranean, Red Sea, Baltic Sea, etc.);
##   (2) parallel processing of multiple taxonomic files by setting diversity values to NA for points
##   in marginal seas while preserving coordinate and metadata information; (3) generation of summary
##   visualizations showing spatial distribution of filtered vs. valid data points; and (4) comprehensive
##   reporting of processing statistics including point counts, filtering percentages, and file summaries
##   for both richness and Shannon diversity datasets.
##
## Required R packages (tested versions):
##   - raster        (version 3.6-26)
##   - sf            (version 1.0-12)
##   - data.table    (version 1.15.4)
##   - ggplot2       (version 3.4.2)
##   - parallel      (version 4.3.2)
##   - dplyr         (version 1.1.2)
##
## ====================================================================================================

# Clear working space
rm(list = ls())

### ===========================
### Libraries
### ===========================
library(raster)      # For spatial raster data handling
library(sf)          # For spatial features and shapefile operations
library(data.table)  # For fast data manipulation and file I/O
library(ggplot2)     # For spatial visualization and plotting
library(parallel)    # For parallel computing infrastructure
library(dplyr)       # For data manipulation and filtering operations

### ===========================
### Configurable Project Paths - adjust as needed
### ===========================
# Load custom functions
source("./Functions/iho_ocean_regions.R")

# Define the marginal seas to remove from global ocean analysis
marginal_seas <- c(
    "Mediterranean Sea - Eastern Basin",
    "Mediterranean Sea - Western Basin", 
    "Red Sea",
    "Baltic Sea",
    "Gulf of Finland",    # Part of Baltic system
    "Gulf of Bothnia",    # Part of Baltic system
    "Gulf of Riga",       # Part of Baltic system
    "Black Sea",
    "Sea of Azov",        # Connected to Black Sea
    "Persian Gulf",
    "Gulf of Aden",       # Extension of Red Sea
    "Gulf of Aqaba",      # Part of Red Sea
    "Caspian Sea",        # Note: might not be in ocean shapefile as it's landlocked
    "Adriatic Sea",       # Part of Mediterranean
    "Aegean Sea",         # Part of Mediterranean
    "Tyrrhenian Sea",     # Part of Mediterranean
    "Ionian Sea",         # Part of Mediterranean
    "Sea of Marmara",     # Part of Mediterranean system
    "Gulf of Suez",       # Part of Red Sea system
    "Ligurian Sea",       # Part of Mediterranean
    "Balearic (Iberian Sea)",  # Part of Mediterranean
    "Alboran Sea"         # Part of Mediterranean
)

# Input and output directories
wd_in <- "./Code/3_Compute_ensembles_and_uncertainties/2_Output/"
wd_out <- "./Code/3_Compute_ensembles_and_uncertainties/2c_Output/"

# Create output directory if it doesn't exist
if(!dir.exists(wd_out)) {
    dir.create(wd_out, recursive = TRUE)
}

## ====================================================================================================
## SECTION 1: Species Richness Marginal Seas Filtering
## ====================================================================================================

# Get all CSV files to process
fnames <- list.files(wd_in, pattern = "csv$", full.names = TRUE)

cat("Found", length(fnames), "CSV files to process:\n")
for(i in seq_along(fnames)) {
    cat(i, ":", basename(fnames[i]), "\n")
}

# Define richness columns to set to NA
richness_columns <- c("mean_month_1", "mean_month_2", "mean_month_3", "mean_month_4", 
                     "mean_month_5", "mean_month_6", "mean_month_7", "mean_month_8",
                     "mean_month_9", "mean_month_10", "mean_month_11", "mean_month_12",
                     "mean_annual", "sd_month_1", "sd_month_2", "sd_month_3", "sd_month_4",
                     "sd_month_5", "sd_month_6", "sd_month_7", "sd_month_8", "sd_month_9",
                     "sd_month_10", "sd_month_11", "sd_month_12", "sd_annual")

# Define function to process a single file
process_single_file <- function(current_file, marginal_seas, richness_columns, wd_out) {
    
    # Load required libraries in worker
    library(data.table)
    library(dplyr)
    library(sf)
    
    # Load functions in worker
    source("./Functions/iho_ocean_regions.R")
    
    # Load data
    df <- fread(current_file)
    original_points <- nrow(df)
    
    # Rename coordinates for function
    names(df)[1:2] <- c("decimalLongitude", "decimalLatitude")
    
    # Apply function to get ocean regions
    df_with_regions <- iho_ocean_regions(df)
    
    # Create filtered dataset by setting marginal seas to NA
    df_filtered <- df_with_regions
    
    # Identify marginal seas points
    marginal_mask <- df_filtered$marine_region %in% marginal_seas
    marginal_points <- sum(marginal_mask)
    valid_points <- sum(!marginal_mask)
    
    # Set richness values to NA for marginal seas
    df_filtered[marginal_mask, richness_columns] <- NA
    
    # Rename coordinates back to x, y
    df_filtered <- df_filtered %>%
        rename(x = decimalLongitude, y = decimalLatitude)
    
    # Generate output filename
    output_filename <- gsub("ensemble_means_monthly_annual_", "marginal_seas_filtered_", basename(current_file))
    output_path <- file.path(wd_out, output_filename)
    
    # Save filtered data
    fwrite(df_filtered, output_path)
    
    # Return summary for this file
    return(data.frame(
        file = basename(current_file),
        original_points = original_points,
        marginal_points = marginal_points,
        valid_points = valid_points,
        output_file = output_filename
    ))
}

# Determine number of cores to use
n_cores <- 20 # Leave one core free
n_cores <- min(n_cores, length(fnames))  # Don't use more cores than files

cat("\nUsing", n_cores, "cores for parallel processing...\n")
cat(paste(rep("=", 80), collapse=""), "\n")
cat("PROCESSING FILES IN PARALLEL\n")
cat(paste(rep("=", 80), collapse=""), "\n")

# Create cluster
cl <- makeCluster(n_cores)

# Export necessary objects to cluster
clusterExport(cl, c("marginal_seas", "richness_columns", "wd_out"))

# Process files in parallel
start_time <- Sys.time()

processing_results <- parLapply(cl, fnames, process_single_file, 
                               marginal_seas = marginal_seas,
                               richness_columns = richness_columns,
                               wd_out = wd_out)

# Stop cluster
stopCluster(cl)

end_time <- Sys.time()
processing_time <- end_time - start_time

# Combine results
processing_summary <- do.call(rbind, processing_results)

# Calculate totals
total_files_processed <- nrow(processing_summary)
total_marginal_points <- sum(processing_summary$marginal_points)
total_original_points <- sum(processing_summary$original_points)
total_valid_points <- sum(processing_summary$valid_points)

# Create summary visualization using the first file
cat("\nCreating summary visualization...\n")
df_example <- fread(file.path(wd_out, processing_summary$output_file[1]))

# Create status column for visualization
marginal_mask_viz <- df_example$marine_region %in% marginal_seas
df_example$data_status <- ifelse(marginal_mask_viz, "Marginal seas (NA)", "Valid ocean data")

p_filtered <- ggplot(df_example, aes(x = x, y = y, color = data_status)) +
    geom_point(size = 0.3, alpha = 0.7) +
    scale_color_manual(values = c("Marginal seas (NA)" = "grey60", "Valid ocean data" = "blue")) +
    labs(
        title = "Marginal Seas Filtering Results",
        subtitle = paste("Example from", processing_summary$file[1]),
        x = "Longitude",
        y = "Latitude", 
        color = "Data Status"
    ) +
    theme_minimal() +
    theme(
        legend.position = "bottom",
        plot.title = element_text(size = 14)
    )

ggsave(file.path(wd_out, "marginal_seas_filtering_summary.png"), 
       plot = p_filtered, width = 16, height = 10, dpi = 300)

# Print final summary
cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("PARALLEL PROCESSING COMPLETE\n") 
cat(paste(rep("=", 80), collapse=""), "\n")
cat("Processing time:", round(processing_time, 2), attr(processing_time, "units"), "\n")
cat("Cores used:", n_cores, "\n")
cat("Total files processed:", total_files_processed, "\n")
cat("Total original points:", total_original_points, "\n")
cat("Total marginal seas points set to NA:", total_marginal_points, "\n")
cat("Total valid ocean points:", total_valid_points, "\n")
cat("Percentage of points filtered:", round(100 * total_marginal_points / total_original_points, 2), "%\n")

cat("\nProcessing summary by file:\n")
print(processing_summary)

# Save processing summary
write.csv(processing_summary, file.path(wd_out, "processing_summary.csv"), row.names = FALSE)

cat("\nFiles saved in:", wd_out, "\n")
cat("- Filtered CSV files: marginal_seas_filtered_*.csv\n")
cat("- Summary visualization: marginal_seas_filtering_summary.png\n") 
cat("- Processing summary: processing_summary.csv\n")

cat("\nMarginal seas removed:\n")
for(sea in marginal_seas) {
    cat("- ", sea, "\n")
}

# Check length of csv files created
created_files <- list.files(wd_out, pattern = "marginal_seas_filtered_.*\\.csv$", full.names = TRUE)
cat("\nNumber of filtered CSV files created:", length(created_files), "\n")

## ====================================================================================================
## SECTION 2: Shannon Diversity Marginal Seas Filtering
## ====================================================================================================

# Get Shannon diversity CSV files to process
fnames_shannon <- list.files(wd_in, pattern = "ensemble_means_monthly_annual_shannon_.*\\.csv$", full.names = TRUE)

cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("SHANNON DIVERSITY MARGINAL SEAS FILTERING\n")
cat(paste(rep("=", 80), collapse=""), "\n")
cat("Found", length(fnames_shannon), "Shannon CSV files to process:\n")
for(i in seq_along(fnames_shannon)) {
    cat(i, ":", basename(fnames_shannon[i]), "\n")
}

if(length(fnames_shannon) > 0) {
    
    # Define Shannon diversity columns to set to NA (same structure as richness)
    shannon_columns <- c("mean_month_1", "mean_month_2", "mean_month_3", "mean_month_4", 
                         "mean_month_5", "mean_month_6", "mean_month_7", "mean_month_8",
                         "mean_month_9", "mean_month_10", "mean_month_11", "mean_month_12",
                         "mean_annual", "sd_month_1", "sd_month_2", "sd_month_3", "sd_month_4",
                         "sd_month_5", "sd_month_6", "sd_month_7", "sd_month_8", "sd_month_9",
                         "sd_month_10", "sd_month_11", "sd_month_12", "sd_annual")
    
    # Define function to process a single Shannon file
    process_single_file_shannon <- function(current_file, marginal_seas, shannon_columns, wd_out) {
        
        # Load required libraries in worker
        library(data.table)
        library(dplyr)
        library(sf)
        
        # Load functions in worker
        source("./Functions/iho_ocean_regions.R")
        
        # Load data
        df <- fread(current_file)
        original_points <- nrow(df)
        
        # Rename coordinates for function
        names(df)[1:2] <- c("decimalLongitude", "decimalLatitude")
        
        # Apply function to get ocean regions
        df_with_regions <- iho_ocean_regions(df)
        
        # Create filtered dataset by setting marginal seas to NA
        df_filtered <- df_with_regions
        
        # Identify marginal seas points
        marginal_mask <- df_filtered$marine_region %in% marginal_seas
        marginal_points <- sum(marginal_mask)
        valid_points <- sum(!marginal_mask)
        
        # Set Shannon diversity values to NA for marginal seas
        df_filtered[marginal_mask, shannon_columns] <- NA
        
        # Rename coordinates back to x, y
        df_filtered <- df_filtered %>%
            rename(x = decimalLongitude, y = decimalLatitude)
        
        # Generate output filename for Shannon
        output_filename <- gsub("ensemble_means_monthly_annual_shannon_", "marginal_seas_filtered_shannon_", basename(current_file))
        output_path <- file.path(wd_out, output_filename)
        
        # Save filtered data
        fwrite(df_filtered, output_path)
        
        # Return summary for this file
        return(data.frame(
            file = basename(current_file),
            original_points = original_points,
            marginal_points = marginal_points,
            valid_points = valid_points,
            output_file = output_filename
        ))
    }
    
    # Determine number of cores to use for Shannon
    n_cores_shannon <- 20 # Leave one core free
    n_cores_shannon <- min(n_cores_shannon, length(fnames_shannon))  # Don't use more cores than files
    
    cat("\nUsing", n_cores_shannon, "cores for Shannon parallel processing...\n")
    cat(paste(rep("=", 80), collapse=""), "\n")
    cat("PROCESSING SHANNON FILES IN PARALLEL\n")
    cat(paste(rep("=", 80), collapse=""), "\n")
    
    # Create cluster for Shannon
    cl_shannon <- makeCluster(n_cores_shannon)
    
    # Export necessary objects to cluster
    clusterExport(cl_shannon, c("marginal_seas", "shannon_columns", "wd_out"))
    
    # Process Shannon files in parallel
    start_time_shannon <- Sys.time()
    
    processing_results_shannon <- parLapply(cl_shannon, fnames_shannon, process_single_file_shannon, 
                                           marginal_seas = marginal_seas,
                                           shannon_columns = shannon_columns,
                                           wd_out = wd_out)
    
    # Stop cluster
    stopCluster(cl_shannon)
    
    end_time_shannon <- Sys.time()
    processing_time_shannon <- end_time_shannon - start_time_shannon
    
    # Combine Shannon results
    processing_summary_shannon <- do.call(rbind, processing_results_shannon)
    
    # Calculate totals for Shannon
    total_files_processed_shannon <- nrow(processing_summary_shannon)
    total_marginal_points_shannon <- sum(processing_summary_shannon$marginal_points)
    total_original_points_shannon <- sum(processing_summary_shannon$original_points)
    total_valid_points_shannon <- sum(processing_summary_shannon$valid_points)
    
    # Create summary visualization using the first Shannon file
    cat("\nCreating Shannon summary visualization...\n")
    df_example_shannon <- fread(file.path(wd_out, processing_summary_shannon$output_file[1]))
    
    # Create status column for visualization
    marginal_mask_viz_shannon <- df_example_shannon$marine_region %in% marginal_seas
    df_example_shannon$data_status <- ifelse(marginal_mask_viz_shannon, "Marginal seas (NA)", "Valid ocean data")
    
    p_filtered_shannon <- ggplot(df_example_shannon, aes(x = x, y = y, color = data_status)) +
        geom_point(size = 0.3, alpha = 0.7) +
        scale_color_manual(values = c("Marginal seas (NA)" = "grey60", "Valid ocean data" = "darkgreen")) +
        labs(
            title = "Shannon Diversity Marginal Seas Filtering Results",
            subtitle = paste("Example from", processing_summary_shannon$file[1]),
            x = "Longitude",
            y = "Latitude", 
            color = "Data Status"
        ) +
        theme_minimal() +
        theme(
            legend.position = "bottom",
            plot.title = element_text(size = 14)
        )
    
    ggsave(file.path(wd_out, "marginal_seas_filtering_summary_shannon.png"), 
           plot = p_filtered_shannon, width = 16, height = 10, dpi = 300)
    
    # Print final summary for Shannon
    cat("\n", paste(rep("=", 80), collapse=""), "\n")
    cat("SHANNON PARALLEL PROCESSING COMPLETE\n") 
    cat(paste(rep("=", 80), collapse=""), "\n")
    cat("Shannon processing time:", round(processing_time_shannon, 2), attr(processing_time_shannon, "units"), "\n")
    cat("Cores used:", n_cores_shannon, "\n")
    cat("Total Shannon files processed:", total_files_processed_shannon, "\n")
    cat("Total original points:", total_original_points_shannon, "\n")
    cat("Total marginal seas points set to NA:", total_marginal_points_shannon, "\n")
    cat("Total valid ocean points:", total_valid_points_shannon, "\n")
    cat("Percentage of Shannon points filtered:", round(100 * total_marginal_points_shannon / total_original_points_shannon, 2), "%\n")
    
    cat("\nShannon processing summary by file:\n")
    print(processing_summary_shannon)
    
    # Save Shannon processing summary
    write.csv(processing_summary_shannon, file.path(wd_out, "processing_summary_shannon.csv"), row.names = FALSE)
    
    cat("\nShannon files saved in:", wd_out, "\n")
    cat("- Filtered Shannon CSV files: marginal_seas_filtered_shannon_*.csv\n")
    cat("- Shannon summary visualization: marginal_seas_filtering_summary_shannon.png\n") 
    cat("- Shannon processing summary: processing_summary_shannon.csv\n")
    
    # Check length of Shannon csv files created
    created_files_shannon <- list.files(wd_out, pattern = "marginal_seas_filtered_shannon_.*\\.csv$", full.names = TRUE)
    cat("\nNumber of filtered Shannon CSV files created:", length(created_files_shannon), "\n")
    
} else {
    cat("No Shannon diversity files found in:", wd_in, "\n")
    cat("Looking for files matching pattern: ensemble_means_monthly_annual_shannon_*.csv\n")
}

## ====================================================================================================
## END OF SCRIPT
## ====================================================================================================