## ====================================================================================================
## This script computes the diversity hotspots for each clade and adds the environmental conditions at these locations.
## Diversity hotspots are defined above the 75th percentile of the global diversity values.
##
## Author:       Dominic Eriksson
## Date:         8th of October 2025
## Affiliation:  Environmental Physics Group UP, ETH Zürich, Switzerland
## Contact:      deriksson@ethz.ch
##
## Input files:  
##   - Files in folder: /Code/3_Compute_ensembles_and_uncertainties/2_Output/Global_annual_richness_ensemble_mean_richness.csv
##   - Files in folder: /Data/Environmental_parameters/Surface_predictors/
##
## Output files: 
##   - CSV files containing the richness hotspot information and the environmental conditions at these locations.
##
## Strategy:
##   This script identifies environmental conditions associated with high-diversity hotspots
##   across microbial clades using global alpha-diversity ensemble data (here, richness).
##   1. Load long-format ensemble data containing annual mean richness for each grid cell.
##   2. Load environmental predictor rasters and compute annual means per variable.
##   3. For each microbial clade, compute clade-specific richness percentiles (65th, 75th, 85th, 95th).
##   4. Identify grid cells exceeding each percentile and convert their coordinates to spatial points.
##   5. Extract environmental variable values at these high-diversity locations from raster data.
##   6. Combine extracted environmental data with richness and percentile information.
##   7. Save the compiled dataset for further analysis of diversity-environment relationships.
##   8. Generate visual checks of high-diversity hotspots and associated environmental variables for each clade and percentile.
##
## Required R packages (tested versions):
##   - raster       (version 3.6.26)
##   - dplyr        (version 1.1.2)
##   - data.table   (version 1.15.4)
##   - sp           (version 2.1.4)
##
## ====================================================================================================


# Clear workspace
rm(list = ls())

# Libraries
library(raster)
library(dplyr)
library(data.table)
library(sp)


# Directories
wd_div <- "Code/3_Compute_ensembles_and_uncertainties/2_Output/Global_annual_richness_ensemble_mean_richness.csv"
wd_meta <- "Environmental_parameters/Surface_predictors/" # navigate to where env. predictor files are saved
wd_out <- "Code/4_Hotspots_diversity/1_Output/"

# Create wd_out directory if not existing
if (!dir.exists(wd_out)) {
    dir.create(wd_out, recursive = TRUE)
}

# Load metadata
fnames_meta <- list.files(wd_meta, full.names = TRUE)

raster_list <- list()
for(i in seq_along( fnames_meta )){

    # Print progress
    print(paste("Processing file", i, "of", length(fnames_meta)))

    # Load daata    
    meta <- brick(fnames_meta[i])

    # Calculate annual mean
    meta_mean <- calc(meta, mean, na.rm = TRUE)

    # Name raster layer
    names(meta_mean) <- gsub(".nc", "", basename(fnames_meta[i]))

    # Save in list
    raster_list[[i]] <- meta_mean
    
}

# Stack all rasters
raster_stack_meta <- stack(raster_list)


# Load data
df_div <- data.table::fread(wd_div)
df_long <- df_div


# Extract environmental conditions for each percentile
results_list <- list()

# Looping setting/vectors
clades <- unique(df_long$clade)
percentiles <- c(0.65, 0.75, 0.85, 0.95)
results_list <- list()
idx <- 1

# Loop though
for (cl in clades) {
  df_clade <- df_long %>% filter(clade == cl)
  quantile_values <- quantile(df_clade$mean_richness, probs = percentiles, na.rm = TRUE)
  
  for (i in seq_along(percentiles)) {
    p <- percentiles[i]
    q <- quantile_values[i]
    
    # Get locations above clade-specific percentile
    high_div <- df_clade %>% filter(mean_richness > q)
    if (nrow(high_div) == 0) next
    
    # Convert to SpatialPoints
    sp_points <- sp::SpatialPoints(high_div[, .(x, y)], proj4string = CRS("+proj=longlat +datum=WGS84"))
    
    # Extract environmental data
    env_vals <- raster::extract(raster_stack_meta, sp_points)
    
    # Combine with coordinates and clade
    env_df <- cbind(high_div, as.data.frame(env_vals))
    env_df$percentile <- p
    
    results_list[[idx]] <- env_df
    idx <- idx + 1
  }
}

# Combine all results
results_df <- bind_rows(results_list)

# Check dataframe
head(results_df)

# Visual checks
library(ggplot2)
library(patchwork) # For multiplot

# Filter for d__Archaea and each percentile
plots <- list()
percentiles_to_plot <- c(0.65, 0.75, 0.85, 0.95)

for (p in percentiles_to_plot) {
  df_archaea <- results_df %>%
    filter(clade == "domain_Archaea", percentile == p)
  
  p_plot <- ggplot(df_archaea, aes(x = x, y = y)) +
    borders("world", colour = "grey70", fill = "grey95") +
    geom_point(color = "red", size = 1.2, alpha = 0.7) +
    coord_fixed() +
    labs(
      title = paste("domain_Archaea Richness >", p, "percentile"),
      x = "Longitude",
      y = "Latitude"
    ) +
    theme_minimal()
  
  plots[[as.character(p)]] <- p_plot
}

# Show all plots in a grid
wrap_plots(plots, ncol = 2)

# Plot temperature for d__Archaea and all percentiles
# List to save plots
plots <- list()
percentiles_to_plot <- c(0.65, 0.75, 0.85, 0.95)

for (p in percentiles_to_plot) {
  df_archaea <- results_df %>%
    filter(clade == "domain_Archaea", percentile == p)
  
  p_plot <- ggplot(df_archaea, aes(x = x, y = y, color = climatology_t_0_50)) +
    borders("world", colour = "grey70", fill = "grey95") +
    geom_point(size = 1.2, alpha = 0.8) +
    scale_color_viridis_c(option = "plasma", na.value = "grey80") +
    coord_fixed() +
    labs(
      title = paste("domain_Archaea >", p, "percentile\nclimatology_t_0_50"),
      x = "Longitude",
      y = "Latitude",
      color = "t_0_50"
    ) +
    theme_minimal()
  
  plots[[as.character(p)]] <- p_plot
}

wrap_plots(plots, ncol = 2)

# Save the results
data.table::fwrite(
  results_df, 
  file.path(wd_out, "EnvStatistics_inDivHotspots_Percentiles_v6.csv"), 
  row.names = FALSE
)

## ====================================================================================================
## END OF SCRIPT
## ====================================================================================================
