## ====================================================================================================
## This R script creates comprehensive Latitudinal Diversity Gradients (LDG) for prokaryotic genera
## by processing ensemble datasets with marginal seas filtering and parallel computing optimization.
## The script aggregates spatial richness and Shannon diversity data across longitudes to create
## latitudinal patterns, computes ensemble statistics with uncertainties, and generates publication-ready
## visualizations showing diversity-latitude relationships for all genera.
##
## Author:       Dominic Eriksson
## Date:         6th of March 2026
## Affiliation:  Environmental Physics Group UP, ETH Zürich, Switzerland
## Contact:      deriksson@ethz.ch
##
## Input files:
##   - Richness ensemble data: ./Code/3_Compute_ensembles_and_uncertainties/1b_Output/monthly_and_annual_richness_ensemble_members_*.csv
##   - Shannon ensemble data: ./Code/3_Compute_ensembles_and_uncertainties/1b_Output/monthly_and_annual_shannon_ensemble_members_*.csv
##   - IHO ocean regions function: ./Functions/iho_ocean_regions.R
##
## Output files:
##   - Individual LDG ensemble members: ./Code/3_Compute_ensembles_and_uncertainties/2e_Output/LDG_ensemble_members_*.csv
##   - Final aggregated LDG data: ./Code/3_Compute_ensembles_and_uncertainties/2e_Output/LDG_final_aggregated_*.csv
##   - Shannon LDG ensemble data: ./Code/3_Compute_ensembles_and_uncertainties/2e_Output/LDG_shannon_*
##   - LDG visualization plots: ./Code/3_Compute_ensembles_and_uncertainties/2e_Output/LDG_*_plot.png
##   - Processing summaries: ./Code/3_Compute_ensembles_and_uncertainties/2e_Output/LDG_*processing_summary.csv
##   - Coordinate mappings: ./Code/3_Compute_ensembles_and_uncertainties/2e_Output/unique_coords_with_regions.csv
##
## Strategy:
##   This script implements a comprehensive LDG analysis pipeline using parallel computing for efficient
##   processing of large ensemble datasets. The workflow includes: (1) loading ensemble member datasets
##   and applying IHO ocean region mapping with marginal seas filtering to exclude enclosed water bodies
##   that may bias global diversity patterns; (2) creating individual LDG profiles by aggregating richness
##   and Shannon diversity values across longitudes for each latitude band, maintaining ensemble member
##   identity for uncertainty quantification; (3) computing final ensemble statistics including means,
##   standard deviations, and sample sizes for robust diversity gradient characterization; (4) generating
##   publication-quality visualizations with confidence intervals showing latitudinal diversity patterns;
##   and (5) implementing parallel processing architecture to handle 272+ genera efficiently with
##   comprehensive quality control and performance monitoring throughout the analysis pipeline.
##
## Required R packages (tested versions):
##   - data.table    (version 1.15.4)
##   - dplyr         (version 1.1.2)
##   - sf            (version 1.0-14)
##   - ggplot2       (version 3.4.2)
##   - viridis       (version 0.6.3)
##   - parallel      (version 4.3.2)
##
## ====================================================================================================

# Clear working space
rm(list = ls())

## ====================================================================================================
## SECTION 1: Richness Latitudinal Diversity Gradients
## ====================================================================================================

### ===========================
### Libraries
### ===========================
library(data.table)  # For fast data manipulation and aggregation
library(dplyr)       # For data transformation and pipeline operations
library(sf)          # For spatial data handling and ocean region processing
library(ggplot2)     # For creating publication-quality LDG visualizations
library(viridis)     # For colorblind-friendly color palettes
library(parallel)    # For parallel processing of multiple genera

### ===========================
### Load Functions
### ===========================
source("./Functions/iho_ocean_regions.R")

### ===========================
### Configurable Project Paths - adjust as needed
### ===========================
wd_in <- "./Code/3_Compute_ensembles_and_uncertainties/1b_Output/"
wd_out <- "./Code/3_Compute_ensembles_and_uncertainties/2e_Output/"

# Create output directory if it doesn't exist
if(!dir.exists(wd_out)) {
    dir.create(wd_out, recursive = TRUE)
    cat("Created output directory:", wd_out, "\n")
}

# Get all long-format CSV files
fnames <- list.files(wd_in, pattern = "monthly_and_annual_richness_ensemble_members_.*\\.csv$", full.names = TRUE)

cat("Found", length(fnames), "ensemble files to process:\n")
for(i in seq_along(fnames)) {
    cat(i, ":", basename(fnames[i]), "\n")
}

if(length(fnames) == 0) {
    stop("No ensemble files found in input directory!")
}

# Define marginal seas to remove (same as before)
marginal_seas <- c(
    "Mediterranean Sea - Eastern Basin",
    "Mediterranean Sea - Western Basin", 
    "Red Sea",
    "Baltic Sea",
    "Gulf of Finland",
    "Gulf of Bothnia",
    "Gulf of Riga", 
    "Black Sea",
    "Sea of Azov",
    "Persian Gulf",
    "Gulf of Aden",
    "Gulf of Aqaba",
    "Caspian Sea",
    "Adriatic Sea",
    "Aegean Sea", 
    "Tyrrhenian Sea",
    "Ionian Sea",
    "Sea of Marmara",
    "Gulf of Suez",
    "Ligurian Sea",
    "Balearic (Iberian Sea)",
    "Alboran Sea"
)

# Define columns to set to NA for marginal seas
richness_columns <- c("month_1", "month_2", "month_3", "month_4", "month_5", "month_6",
                     "month_7", "month_8", "month_9", "month_10", "month_11", "month_12", 
                     "annual_mean")

# Define function to process a single file
process_single_genus <- function(current_file, marginal_seas, richness_columns, wd_out) {
    
    # Load required libraries in worker
    library(data.table)
    library(dplyr)
    library(sf)
    library(ggplot2)
    library(viridis)
    
    # Load functions in worker
    source("./Functions/iho_ocean_regions.R")
    
    # Extract genus name from filename
    genus_name <- gsub("monthly_and_annual_richness_ensemble_members_", "", basename(current_file))
    genus_name <- gsub("\\.csv$", "", genus_name)
    
    cat("Processing:", genus_name, "\n")
    
    # Load the long-format data
    df_long <- fread(current_file)
    original_dims <- dim(df_long)
    
    # Create a unique coordinate dataframe for the IHO function
    unique_coords <- df_long[, .(x, y)] %>% unique()
    
    # Only run IHO function once per coordinate set if not already done
    coords_file <- file.path(wd_out, "unique_coords_with_regions.csv")
    if(file.exists(coords_file)) {
        # Load existing coordinate mapping
        coords_with_regions <- fread(coords_file)
    } else {
        # Create coordinate mapping for the first time
        unique_coords_renamed <- copy(unique_coords)
        names(unique_coords_renamed) <- c("decimalLongitude", "decimalLatitude")
        
        # Apply IHO ocean regions function
        coords_with_regions <- iho_ocean_regions(unique_coords_renamed)
        
        # Rename back and select needed columns
        coords_with_regions <- coords_with_regions %>%
            rename(x = decimalLongitude, y = decimalLatitude) %>%
            dplyr::select(x, y, marine_region)
        
        # Save for reuse
        fwrite(coords_with_regions, coords_file)
    }
    
    # Merge marine regions with original data
    df_long_with_regions <- merge(df_long, coords_with_regions, by = c("x", "y"), all.x = TRUE)
    
    # Identify marginal seas
    marginal_mask <- df_long_with_regions$marine_region %in% marginal_seas
    
    # Set marginal seas richness values to NA
    df_filtered <- copy(df_long_with_regions)
    df_filtered[marginal_mask, (richness_columns) := NA]
    
    # Remove marine_region column since it's no longer needed
    df_filtered[, marine_region := NULL]
    
    # Step 1: Create LDG for each ensemble member (aggregate by latitude, algorithm, and bootstrap)
    month_cols <- paste0("month_", 1:12)
    all_richness_cols <- c(month_cols, "annual_mean")
    
    # Aggregate across longitude (x) for each y, alg, and bootstrap combination
    ldg_ensemble_members <- df_filtered[, lapply(.SD, function(x) mean(x, na.rm = TRUE)), 
                                       by = .(y, alg, bootstrap, clade), 
                                       .SDcols = all_richness_cols]
    
    # Replace NaN with NA (occurs when all values are NA)
    for(col in all_richness_cols) {
        ldg_ensemble_members[is.nan(get(col)), (col) := NA]
    }
    
    # Save individual ensemble member LDGs
    ensemble_filename <- paste0("LDG_ensemble_members_", genus_name, ".csv")
    ensemble_path <- file.path(wd_out, ensemble_filename)
    fwrite(ldg_ensemble_members, ensemble_path)
    
    # Step 2: Create final aggregated LDG
    ldg_final <- ldg_ensemble_members[, .(
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
        sd_annual = sd(annual_mean, na.rm = TRUE),
        n_ensemble_members = .N
    ), by = .(y, clade)]
    
    # Replace NaN with NA
    mean_cols <- paste0("mean_", c(month_cols, "annual"))
    sd_cols <- paste0("sd_", c(month_cols, "annual"))
    all_summary_cols <- c(mean_cols, sd_cols)
    
    for(col in all_summary_cols) {
        ldg_final[is.nan(get(col)), (col) := NA]
    }
    
    # Save final aggregated LDG
    final_filename <- paste0("LDG_final_aggregated_", genus_name, ".csv")
    final_path <- file.path(wd_out, final_filename)
    fwrite(ldg_final, final_path)
    
    # Create visualization (save one plot per genus)
    if(!all(is.na(ldg_final$mean_annual))) {
        # Clean genus name for display
        genus_display <- gsub("genus_", "", unique(ldg_final$clade))
        genus_display <- gsub("_", " ", genus_display)
        
        # Plot annual mean LDG
        p_ldg <- ggplot(ldg_final, aes(x = y, y = mean_annual)) +
            geom_line(size = 1, color = "darkblue") +
            geom_ribbon(aes(ymin = pmax(0, mean_annual - sd_annual), 
                           ymax = mean_annual + sd_annual), 
                        alpha = 0.3, fill = "lightblue") +
            geom_point(size = 0.8, color = "darkblue") +
            labs(
                title = "Latitudinal Diversity Gradient - Annual Mean",
                subtitle = paste("Genus:", genus_display),
                x = "Latitude (°)",
                y = "Mean Annual Richness",
                caption = "Ribbon shows ±1 standard deviation"
            ) +
            theme_minimal() +
            theme(
                plot.title = element_text(size = 14, hjust = 0.5),
                plot.subtitle = element_text(size = 12, hjust = 0.5),
                axis.title = element_text(size = 10),
                axis.text = element_text(size = 8)
            )
        
        plot_filename <- paste0("LDG_", genus_name, "_plot.png")
        ggsave(file.path(wd_out, plot_filename), p_ldg, width = 10, height = 6, dpi = 300)
    }
    
    # Return summary for this genus
    return(data.frame(
        genus = genus_name,
        original_rows = original_dims[1],
        filtered_rows = nrow(df_filtered),
        marginal_seas_rows = sum(marginal_mask),
        ldg_latitudes = nrow(ldg_final),
        peak_richness = ifelse(all(is.na(ldg_final$mean_annual)), NA, max(ldg_final$mean_annual, na.rm = TRUE)),
        peak_latitude = ifelse(all(is.na(ldg_final$mean_annual)), NA, ldg_final[which.max(mean_annual), y]),
        ensemble_file = ensemble_filename,
        final_file = final_filename,
        stringsAsFactors = FALSE
    ))
}

# Determine number of cores to use
n_cores <- detectCores() - 2  # Leave two cores free
n_cores <- min(n_cores, length(fnames))  # Don't use more cores than files

cat("\nUsing", n_cores, "cores for parallel processing...\n")
cat(paste(rep("=", 80), collapse=""), "\n")
cat("PROCESSING ALL GENERA IN PARALLEL\n")
cat(paste(rep("=", 80), collapse=""), "\n")

# Create cluster
cl <- makeCluster(n_cores)

# Export necessary objects to cluster
clusterExport(cl, c("marginal_seas", "richness_columns", "wd_out"))

# Process files in parallel
start_time <- Sys.time()

processing_results <- parLapply(cl, fnames, process_single_genus, 
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
total_genera <- nrow(processing_summary)
total_original_rows <- sum(processing_summary$original_rows)
total_marginal_seas_rows <- sum(processing_summary$marginal_seas_rows)

# Save processing summary
write.csv(processing_summary, file.path(wd_out, "LDG_processing_summary.csv"), row.names = FALSE)

# Print final summary
cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("PARALLEL LDG PROCESSING COMPLETE\n") 
cat(paste(rep("=", 80), collapse=""), "\n")
cat("Processing time:", round(processing_time, 2), attr(processing_time, "units"), "\n")
cat("Cores used:", n_cores, "\n")
cat("Total genera processed:", total_genera, "\n")
cat("Total original data rows:", total_original_rows, "\n")
cat("Total marginal seas rows filtered:", total_marginal_seas_rows, "\n")
cat("Percentage of data in marginal seas:", round(100 * total_marginal_seas_rows / total_original_rows, 2), "%\n")

# Show top 10 genera by peak richness
top_genera <- processing_summary[order(-peak_richness, na.last = TRUE)[1:min(10, nrow(processing_summary))], ]
cat("\nTop 10 genera by peak richness:\n")
print(top_genera[, c("genus", "peak_richness", "peak_latitude")])

cat("\nFiles created in", wd_out, ":\n")
cat("- LDG ensemble member files: LDG_ensemble_members_*.csv\n")
cat("- Final aggregated LDG files: LDG_final_aggregated_*.csv\n")
cat("- LDG visualization plots: LDG_*_plot.png\n")
cat("- Processing summary: LDG_processing_summary.csv\n")
cat("- Coordinate mapping (reused): unique_coords_with_regions.csv\n")

if(total_genera != 272) {
    cat("\nWARNING: Expected 272 genera but processed", total_genera, "\n")
} else {
    cat("\nSUCCESS: All 272 genera processed successfully!\n")
}

# Clean up the coordinate mapping file if desired (optional)
# file.remove(file.path(wd_out, "unique_coords_with_regions.csv"))

## ====================================================================================================
## SECTION 2: Shannon Diversity Latitudinal Diversity Gradients
## ====================================================================================================

# Get all Shannon ensemble files
fnames_shannon <- list.files(wd_in, pattern = "monthly_and_annual_shannon_ensemble_members_.*\\.csv$", full.names = TRUE)

cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("SHANNON DIVERSITY LDG PROCESSING\n")
cat(paste(rep("=", 80), collapse=""), "\n")
cat("Found", length(fnames_shannon), "Shannon ensemble files to process:\n")
for(i in seq_along(fnames_shannon)) {
    cat(i, ":", basename(fnames_shannon[i]), "\n")
}

if(length(fnames_shannon) > 0) {
    
    # Define Shannon columns to set to NA for marginal seas
    shannon_columns <- c("month_1", "month_2", "month_3", "month_4", "month_5", "month_6",
                         "month_7", "month_8", "month_9", "month_10", "month_11", "month_12", 
                         "annual_mean")
    
    # Define function to process a single Shannon genus
    process_single_genus_shannon <- function(current_file, marginal_seas, shannon_columns, wd_out) {
        
        # Load required libraries in worker
        library(data.table)
        library(dplyr)
        library(sf)
        library(ggplot2)
        library(viridis)
        
        # Load functions in worker
        source("./Functions/iho_ocean_regions.R")
        
        # Extract genus name from filename
        genus_name <- gsub("monthly_and_annual_shannon_ensemble_members_", "", basename(current_file))
        genus_name <- gsub("\\.csv$", "", genus_name)
        
        cat("Processing Shannon:", genus_name, "\n")
        
        # Load the long-format Shannon data
        df_long <- fread(current_file)
        original_dims <- dim(df_long)
        
        # Create a unique coordinate dataframe for the IHO function
        unique_coords <- df_long[, .(x, y)] %>% unique()
        
        # Use the existing coordinate mapping file if available
        coords_file <- file.path(wd_out, "unique_coords_with_regions.csv")
        if(file.exists(coords_file)) {
            # Load existing coordinate mapping
            coords_with_regions <- fread(coords_file)
        } else {
            # Create coordinate mapping for Shannon (same coordinates as richness)
            unique_coords_renamed <- copy(unique_coords)
            names(unique_coords_renamed) <- c("decimalLongitude", "decimalLatitude")
            
            # Apply IHO ocean regions function
            coords_with_regions <- iho_ocean_regions(unique_coords_renamed)
            
            # Rename back and select needed columns
            coords_with_regions <- coords_with_regions %>%
                rename(x = decimalLongitude, y = decimalLatitude) %>%
                dplyr::select(x, y, marine_region)
            
            # Save for reuse
            fwrite(coords_with_regions, coords_file)
        }
        
        # Merge marine regions with original Shannon data
        df_long_with_regions <- merge(df_long, coords_with_regions, by = c("x", "y"), all.x = TRUE)
        
        # Identify marginal seas
        marginal_mask <- df_long_with_regions$marine_region %in% marginal_seas
        
        # Set marginal seas Shannon values to NA
        df_filtered <- copy(df_long_with_regions)
        df_filtered[marginal_mask, (shannon_columns) := NA]
        
        # Remove marine_region column since it's no longer needed
        df_filtered[, marine_region := NULL]
        
        # Step 1: Create Shannon LDG for each ensemble member
        month_cols <- paste0("month_", 1:12)
        all_shannon_cols <- c(month_cols, "annual_mean")
        
        # Aggregate across longitude (x) for each y, alg, and bootstrap combination
        ldg_ensemble_members <- df_filtered[, lapply(.SD, function(x) mean(x, na.rm = TRUE)), 
                                           by = .(y, alg, bootstrap, clade), 
                                           .SDcols = all_shannon_cols]
        
        # Replace NaN with NA (occurs when all values are NA)
        for(col in all_shannon_cols) {
            ldg_ensemble_members[is.nan(get(col)), (col) := NA]
        }
        
        # Save individual Shannon ensemble member LDGs
        ensemble_filename <- paste0("LDG_shannon_ensemble_members_", genus_name, ".csv")
        ensemble_path <- file.path(wd_out, ensemble_filename)
        fwrite(ldg_ensemble_members, ensemble_path)
        
        # Step 2: Create final aggregated Shannon LDG
        ldg_final <- ldg_ensemble_members[, .(
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
            sd_annual = sd(annual_mean, na.rm = TRUE),
            n_ensemble_members = .N
        ), by = .(y, clade)]
        
        # Replace NaN with NA
        mean_cols <- paste0("mean_", c(month_cols, "annual"))
        sd_cols <- paste0("sd_", c(month_cols, "annual"))
        all_summary_cols <- c(mean_cols, sd_cols)
        
        for(col in all_summary_cols) {
            ldg_final[is.nan(get(col)), (col) := NA]
        }
        
        # Save final aggregated Shannon LDG
        final_filename <- paste0("LDG_shannon_final_aggregated_", genus_name, ".csv")
        final_path <- file.path(wd_out, final_filename)
        fwrite(ldg_final, final_path)
        
        # Create Shannon visualization (save one plot per genus)
        if(!all(is.na(ldg_final$mean_annual))) {
            # Clean genus name for display
            genus_display <- gsub("genus_", "", unique(ldg_final$clade))
            genus_display <- gsub("_", " ", genus_display)
            
            # Plot annual mean Shannon LDG
            p_ldg <- ggplot(ldg_final, aes(x = y, y = mean_annual)) +
                geom_line(size = 1, color = "darkorchid") +
                geom_ribbon(aes(ymin = pmax(0, mean_annual - sd_annual), 
                               ymax = mean_annual + sd_annual), 
                            alpha = 0.3, fill = "plum") +
                geom_point(size = 0.8, color = "darkorchid") +
                labs(
                    title = "Latitudinal Shannon Diversity Gradient - Annual Mean",
                    subtitle = paste("Genus:", genus_display),
                    x = "Latitude (°)",
                    y = "Mean Annual Shannon Diversity",
                    caption = "Ribbon shows ±1 standard deviation"
                ) +
                theme_minimal() +
                theme(
                    plot.title = element_text(size = 14, hjust = 0.5),
                    plot.subtitle = element_text(size = 12, hjust = 0.5),
                    axis.title = element_text(size = 10),
                    axis.text = element_text(size = 8)
                )
            
            plot_filename <- paste0("LDG_shannon_", genus_name, "_plot.png")
            ggsave(file.path(wd_out, plot_filename), p_ldg, width = 10, height = 6, dpi = 300)
        }
        
        # Return summary for this Shannon genus
        return(data.frame(
            genus = genus_name,
            original_rows = original_dims[1],
            filtered_rows = nrow(df_filtered),
            marginal_seas_rows = sum(marginal_mask),
            ldg_latitudes = nrow(ldg_final),
            peak_shannon = ifelse(all(is.na(ldg_final$mean_annual)), NA, max(ldg_final$mean_annual, na.rm = TRUE)),
            peak_latitude = ifelse(all(is.na(ldg_final$mean_annual)), NA, ldg_final[which.max(mean_annual), y]),
            ensemble_file = ensemble_filename,
            final_file = final_filename,
            stringsAsFactors = FALSE
        ))
    }
    
    # Determine number of cores to use for Shannon
    n_cores_shannon <- detectCores() - 2  # Leave two cores free
    n_cores_shannon <- min(n_cores_shannon, length(fnames_shannon))  # Don't use more cores than files
    
    cat("\nUsing", n_cores_shannon, "cores for Shannon parallel processing...\n")
    cat(paste(rep("=", 80), collapse=""), "\n")
    cat("PROCESSING ALL SHANNON GENERA IN PARALLEL\n")
    cat(paste(rep("=", 80), collapse=""), "\n")
    
    # Create cluster for Shannon
    cl_shannon <- makeCluster(n_cores_shannon)
    
    # Export necessary objects to cluster
    clusterExport(cl_shannon, c("marginal_seas", "shannon_columns", "wd_out"))
    
    # Process Shannon files in parallel
    start_time_shannon <- Sys.time()
    
    processing_results_shannon <- parLapply(cl_shannon, fnames_shannon, process_single_genus_shannon, 
                                           marginal_seas = marginal_seas,
                                           shannon_columns = shannon_columns,
                                           wd_out = wd_out)
    
    # Stop cluster
    stopCluster(cl_shannon)
    
    end_time_shannon <- Sys.time()
    processing_time_shannon <- end_time_shannon - start_time_shannon
    
    # Combine Shannon results
    processing_summary_shannon <- do.call(rbind, processing_results_shannon)
    
    # Calculate Shannon totals
    total_genera_shannon <- nrow(processing_summary_shannon)
    total_original_rows_shannon <- sum(processing_summary_shannon$original_rows)
    total_marginal_seas_rows_shannon <- sum(processing_summary_shannon$marginal_seas_rows)
    
    # Save Shannon processing summary
    write.csv(processing_summary_shannon, file.path(wd_out, "LDG_shannon_processing_summary.csv"), row.names = FALSE)
    
    # Print final Shannon summary
    cat("\n", paste(rep("=", 80), collapse=""), "\n")
    cat("PARALLEL SHANNON LDG PROCESSING COMPLETE\n") 
    cat(paste(rep("=", 80), collapse=""), "\n")
    cat("Shannon processing time:", round(processing_time_shannon, 2), attr(processing_time_shannon, "units"), "\n")
    cat("Cores used:", n_cores_shannon, "\n")
    cat("Total Shannon genera processed:", total_genera_shannon, "\n")
    cat("Total original Shannon data rows:", total_original_rows_shannon, "\n")
    cat("Total Shannon marginal seas rows filtered:", total_marginal_seas_rows_shannon, "\n")
    cat("Percentage of Shannon data in marginal seas:", round(100 * total_marginal_seas_rows_shannon / total_original_rows_shannon, 2), "%\n")
    
    # Show top 10 Shannon genera by peak diversity
    top_genera_shannon <- processing_summary_shannon[order(-peak_shannon, na.last = TRUE)[1:min(10, nrow(processing_summary_shannon))], ]
    cat("\nTop 10 Shannon genera by peak diversity:\n")
    print(top_genera_shannon[, c("genus", "peak_shannon", "peak_latitude")])
    
    cat("\nShannon files created in", wd_out, ":\n")
    cat("- Shannon LDG ensemble member files: LDG_shannon_ensemble_members_*.csv\n")
    cat("- Final aggregated Shannon LDG files: LDG_shannon_final_aggregated_*.csv\n")
    cat("- Shannon LDG visualization plots: LDG_shannon_*_plot.png\n")
    cat("- Shannon processing summary: LDG_shannon_processing_summary.csv\n")
    
    if(total_genera_shannon != 272) {
        cat("\nWARNING: Expected 272 Shannon genera but processed", total_genera_shannon, "\n")
    } else {
        cat("\nSUCCESS: All 272 Shannon genera processed successfully!\n")
    }
    
} else {
    cat("No Shannon diversity files found in:", wd_in, "\n")
    cat("Looking for files matching pattern: monthly_and_annual_shannon_ensemble_members_*.csv\n")
}

## ====================================================================================================
## END OF SCRIPT
## ====================================================================================================