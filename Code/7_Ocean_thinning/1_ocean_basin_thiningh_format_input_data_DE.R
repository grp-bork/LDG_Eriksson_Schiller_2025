## ====================================================================================================
## This R script performs basin-stratified subsampling of prokaryotic diversity data to evaluate
## the effects of sampling bias on biodiversity patterns. It creates multiple downsampled datasets
## with different sample size caps across ocean basins, enabling sensitivity analysis of how
## uneven sampling effort affects species distribution modeling outcomes.
##
## Author:       Dominic Eriksson
## Date:         26th of February 2026
## Affiliation:  Environmental Physics Group UP, ETH Zürich, Switzerland
## Contact:      deriksson@ethz.ch
##
## Input files:  
##   - Alpha diversity indices CSV with prokaryotic richness data
##   - Metadata CSV with sample locations and environmental information
##   - Ocean regions classification function
##
## Output files: 
##   - RDS file with all subsampled datasets for different cap values
##   - Combined CSV file formatted for CEPHALOPOD modeling pipeline
##   - Comparative visualization maps showing subsampling effects
##
## Strategy:
##   The script applies basin-stratified subsampling with multiple cap values (10-350 samples
##   per basin) across 100 iterations each, creating 1,100 subsampled datasets. This approach
##   enables systematic evaluation of sampling bias effects while maintaining representative
##   coverage across ocean basins. The resulting datasets facilitate robust assessment of
##   how sampling heterogeneity influences biodiversity pattern detection.
##
## Required R packages (tested versions):
##   - dplyr         1.1.2
##   - data.table    1.14.8
##   - ggplot2       3.4.2
##   - sf            1.0.12
##   - rnaturalearth 0.3.3
## ====================================================================================================

# Clear workspace
rm(list = ls())

# Load libraries
library(dplyr)      # For data manipulation
library(data.table) # For efficient data handling

# Set working directories
wd_in <- "Data/Jonas_Richter/Alpha diversity indices/alpha_diversity_prokaryotes_domain_phylum_class_order_family_genus_rarefied5000_V3.csv"

# Load data
df <- data.table::fread(wd_in)

# Subset clade of interest - biosample, strings with prokaryotes and strings with domain - We only do prokaryotes at the moment
df_diversity <- df %>%
    select(biosample, prokaryotes_Richness_5000)

## Load data - NOTE: Adjust metadata version
df_metadata <- data.table::fread("Data/Jonas_Richter/Metadaten/metadata_V23.csv")
# Subset data of interest in metadata
df_metadata <- df_metadata[, c("biosample", "depth", "month", "year", "lon", "lat", "analysis_type", "distance_to_coast")]

## Merge Jonas diversity with metadata
df <- base::merge(df_diversity, df_metadata, by = "biosample", all.x = TRUE)

# Columns for CEPHALOPOD
cols <- c(
    "scientificname",
    "worms_id",
    "decimallatitude",
    "decimallongitude",
    "depth",
    "year",
    "month",
    "measurementvalue",
    "measurementunit",
    "taxonrank"
)

## Surface to MLD --------------------------------------------------------------------------------------------

# Set output directory
wd_out <- "Code/7_Ocean_thinning/1_Output/"

# Filter for surface data
unique(df$analysis_type)
df_surface <- df[which(df$analysis_type != "not_included"), ]
df_surface <- df_surface[analysis_type %in% c("MLD")]

# Rename coordinate columns
setnames(df_surface, c("lon", "lat"), c("decimalLongitude", "decimalLatitude"))

# Load function to add ocean regions
source("../../Functions/subset_ocean_regions.R")

# Add ocean regions
df_surface_regions <- subset_ocean_regions(df_surface)

# Convert to data.table
df_surface_regions <- data.table::as.data.table(df_surface_regions)

# Count number of samples per marine region
sample_counts <- df_surface_regions[, .N, by = marine_region][order(-N)]
print("Number of samples per marine region:")
print(sample_counts)

sum(sample_counts$N)  # Total number of samples with MLD data

### ====================================================================================================
### Basin-stratified subsampling sensitivity analysis
### ====================================================================================================

# Set reproducible seed
set.seed(123)

# Define subsampling parameters
caps <- c(350, 300, 250, 200, 150, 100, 50, 25, 20, 15, 10)
n_iterations <- 100

# Initialize list to store results
subsampled_datasets <- list()

# Create progress counter
total_combinations <- length(caps) * n_iterations
current_combination <- 0

cat("Starting basin-stratified subsampling...\n")
cat("Caps:", caps, "\n")
cat("Iterations per cap:", n_iterations, "\n")
cat("Total combinations:", total_combinations, "\n\n")

# Create all subsampled datasets
for(cap in caps) {
  cat("Processing cap =", cap, "samples per basin...\n")
  
  for(iter in 1:n_iterations) {
    current_combination <- current_combination + 1
    
    # Create subsample_id
    subsample_id <- sprintf("cap%d_iter%03d", cap, iter)
    
    # Subsample each marine region according to rules
    subsampled_data <- df_surface_regions[
      , .SD[if(.N > cap) sample(.N, cap) else 1:.N], 
      by = marine_region
    ]
    
    # Add identifier columns
    subsampled_data[, subsample_cap := cap]
    subsampled_data[, subsample_iter := iter]
    subsampled_data[, subsample_id := subsample_id]
    
    # Store in list with descriptive name
    subsampled_datasets[[subsample_id]] <- subsampled_data
    
    # Progress update every 50 iterations
    if(current_combination %% 50 == 0) {
      cat("  Completed", current_combination, "of", total_combinations, "combinations\n")
    }
  }
}

# Final summary
cat("\nSubsampling completed successfully!\n")
cat("Created", length(subsampled_datasets), "subsampled datasets\n")

# Validation: Check structure of first few datasets
cat("\nValidation - sample sizes in first dataset:\n")
first_dataset <- subsampled_datasets[[1]]
validation_counts <- first_dataset[, .N, by = .(marine_region, subsample_cap)][order(-N)]
print(validation_counts)

cat("\nValidation - identifier columns in first dataset:\n")
cat("subsample_cap:", unique(first_dataset$subsample_cap), "\n")
cat("subsample_iter:", unique(first_dataset$subsample_iter), "\n") 
cat("subsample_id:", unique(first_dataset$subsample_id), "\n")

# Save the subsampled datasets
saveRDS(subsampled_datasets, file.path(wd_out, "subsampled_datasets_basin_stratified.rds"))
cat("\nSubsampled datasets saved to:", file.path(wd_out, "subsampled_datasets_basin_stratified.rds"), "\n")

### ====================================================================================================
### Plot maps comparing highest and lowest cap subsampling
### ====================================================================================================

# Load required libraries for mapping
library(ggplot2)       # For visualization
library(sf)            # For spatial data handling
library(rnaturalearth) # For world map data

# Extract first (highest cap) and last (lowest cap) datasets
first_dataset <- subsampled_datasets[[1]]  # cap350_iter001
last_dataset <- subsampled_datasets[[length(subsampled_datasets)]]  # cap200_iter100

# Get world map
world <- ne_countries(scale = "medium", returnclass = "sf")

# Add dataset identifier for plotting
first_dataset[, dataset := paste0("Cap ", unique(subsample_cap), " (", nrow(first_dataset), " samples)")]
last_dataset[, dataset := paste0("Cap ", unique(subsample_cap), " (", nrow(last_dataset), " samples)")]

# Combine datasets for plotting
combined_data <- rbind(first_dataset, last_dataset)

# Create side-by-side comparison map
p1 <- ggplot(combined_data) +
  geom_sf(data = world, fill = "grey90", color = "white", size = 0.2) +
  geom_point(aes(x = decimalLongitude, y = decimalLatitude, color = marine_region), 
             size = 0.8, alpha = 0.7) +
  facet_wrap(~dataset, ncol = 2) +
  labs(
    title = "Basin-stratified subsampling comparison",
    subtitle = "Sample distribution between highest and lowest cap values",
    color = "Marine Region",
    x = "Longitude", 
    y = "Latitude"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    strip.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12)
  ) +
  guides(color = guide_legend(override.aes = list(size = 2), ncol = 4))

# Print the plot
print(p1)

# Save the comparison plot
ggsave(
  filename = file.path(wd_out, "subsampling_comparison_map.png"),
  plot = p1,
  width = 16, height = 8, 
  dpi = 300
)

# Summary statistics for comparison
cat("\nSample size comparison:\n")
cat("Highest cap (", unique(first_dataset$subsample_cap), "):", nrow(first_dataset), "total samples\n")
cat("Lowest cap (", unique(last_dataset$subsample_cap), "):", nrow(last_dataset), "total samples\n")
cat("Reduction:", nrow(first_dataset) - nrow(last_dataset), "samples (",
    round((1 - nrow(last_dataset)/nrow(first_dataset)) * 100, 1), "% reduction)\n")

# Basin-wise comparison
cat("\nSamples per basin comparison:\n")
first_counts <- first_dataset[, .N, by = marine_region][order(-N)]
last_counts <- last_dataset[, .N, by = marine_region][order(-N)]
comparison <- merge(first_counts, last_counts, by = "marine_region", suffixes = c("_cap350", "_cap200"))
print(comparison)

# Format for Cephalopod
list_dataframes <- list()
for(i in seq_along(subsampled_datasets)) {

    # Print progress
    message("Formatting dataset: ", names(subsampled_datasets)[i])
  
    # Subset
    df_sub <- subsampled_datasets[[i]]

    # Create dataframe for CEPHALOPOD
    df_final <- data.frame(
        scientificname = paste0("prokaryotes_Richness_5000_", names(subsampled_datasets)[i]),
        worms_id = paste0("prokaryotes_Richness_5000_", names(subsampled_datasets)[i]),
        decimallatitude = df_sub$decimalLatitude,
        decimallongitude = df_sub$decimalLongitude,
        depth = df_sub$depth,
        year = df_sub$year,
        month = df_sub$month,
        measurementvalue = df_sub$prokaryotes_Richness_5000,
        measurementunit = paste0("prokaryotes_Richness_5000_", names(subsampled_datasets)[i]),
        taxonrank = paste0("prokaryotes_Richness_5000_", names(subsampled_datasets)[i])
    ) 
    list_dataframes[[i]] <- df_final
}

# Merge all dataframes into one
df_all_subsampled <- do.call(rbind, list_dataframes)
dim(df_all_subsampled)
head(df_all_subsampled)

# Save combined data
write.csv(
    df_all_subsampled, 
    paste0(wd_out, "prokaryotes_Richness_5000_all_basin_stratified_subsamples_MLD.csv"), 
    row.names = FALSE)

## ====================================================================================================
## END OF SCRIPT
## ====================================================================================================  