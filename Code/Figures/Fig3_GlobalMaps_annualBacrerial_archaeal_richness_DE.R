## ====================================================================================================
## This script generates global maps of modeled bacterial and archaeal richness with uncertainty representation.
##
## Author:       Dominic Eriksson
## Date:         8th of October 2025
## Affiliation:  Environmental Physics Group UP, ETH Zürich, Switzerland
## Contact:      deriksson@ethz.ch
##
## Input files:  
##   - Global ensemble mean richness CSV: Code/3_Compute_ensembles_and_uncertainties/2_Output/Global_annual_richness_ensemble_mean_richness.csv
##   - Ocean regions shapefile: GOaS_v1_20211214/goas_v01.shp
##   - Model statistics RData files: 2_Extract_model_output/3_Output/Non_log/domain_Bacteria_Richness_5000_v2.RData
##                                   2_Extract_model_output/3_Output/Non_log/domain_Archaea_Richness_5000_v2.RData
##
## Output files: 
##   - Global maps of bacterial and archaeal richness
##
## Strategy:
##   This script visualizes global maps in microbial richness and model uncertainty:
##   1. Load ensemble mean richness data for bacterial and archaeal clades.
##   2. Transform spatial data into Pacific-centered projection and read ocean basin shapefiles.
##   3. Compute the coefficient of variation (CV) for each spatial grid cell:
##       - CV = (Standard Deviation / Mean) * 100
##       - Apply minimum mean richness threshold to exclude low-richness regions.
##   4. Define uncertainty thresholds:
##       - Absolute CV threshold (default 20%)
##       - Relative threshold: 70th percentile of CV values per clade
##       - Final CV threshold = maximum of absolute and relative thresholds
##   5. Generate stippling masks for cells exceeding the CV threshold and above the mean richness threshold.
##   6. Convert stippling points into Pacific-centered CRS for plotting.
##   7. Create contour intervals based on quantiles of richness values to represent spatial gradients.
##   8. Combine raster, contour, and stippling layers into a ggplot figure.
##   9. Annotate plots with model predictive power (e.g., R^2).
##  10. Save high-resolution SVG maps for publication.
##
## Required R packages (tested versions):
##   - ggplot2      (version 3.5.2)
##   - dplyr        (version 1.1.2)
##   - sf           (version 1.0.14)
##   - tidyterra    (version 0.4.0)
##   - cowplot      (version 1.1.1)
##   - terra        (version 1.7.78)
##   - hexbin       (version 1.28.3)
##   - data.table   (version 1.15.4)
##
## ====================================================================================================



### ==============================================================================
### Global maps annual bacterial and archaeal richness
### ==============================================================================

# Clear workspace
rm(list = ls())

# Load libraries
library(ggplot2)
library(dplyr)
library(sf)
library(tidyterra)
library(cowplot)
library(terra)
library(hexbin)
library(data.table)

# Working directories - adjust as needed
wd_in <- "Code/3_Compute_ensembles_and_uncertainties/2_Output/Global_annual_richness_ensemble_mean_richness.csv"
wd_out <- "Code/Figures/Fig3/"

# Create output directory if it doesn't exist
if (!dir.exists(wd_out)) {
  dir.create(wd_out, recursive = TRUE)
}

# Filenames for model stats - From the script 1_Extract_model_output/3_Extract_ensemble_means_with_predictivePower.R
fnames_model_stats <- c(
    "Code/2_Extract_model_output/3_Output/Non_log/domain_Bacteria_Richness_5000_v2.RData",
    "Code/2_Extract_model_output/3_Output/Non_log/domain_Archaea_Richness_5000_v2.RData"
)


## Required data & data formatting
# Read shapefile with geometries of ocean basins - Note: Specific shapefile needed for Pacific-centered map
sf.ocean <- sf::st_read("Data/Ocean_regions_shapefile/GOaS_v1_20211214/goas_v01.shp") # direct to the folder where the ocean shapefile is saved, downloadable from marine.org (Global Oceans and seas)
# Pacific centered world map: Source: https://stackoverflow.com/questions/56146735/visual-bug-when-changing-robinson-projections-central-meridian-with-ggplot2
worldMap <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
  sf::st_make_valid()

# Set projection pacific centered
lon_pacific <- '200' 
target_crs <- sf::st_crs( paste0("+proj=eqc +x_0=0 +y_0=0 +lat_0=0 +lon_0=", lon_pacific) )

# define a long & slim polygon that overlaps the meridian line & set its CRS to match
# that of world
# Centered in lon 200
offset <- 180 - as.numeric(lon_pacific)


polygon <- st_polygon(x = list(rbind(
  c(-0.0001 - offset, 90),
  c(0 - offset, 90),
  c(0 - offset, -90),
  c(-0.0001 - offset, -90),
  c(-0.0001 - offset, 90)
))) %>%
  st_sfc() %>%
  st_set_crs(4326)


# modify world dataset to remove overlapping portions with world's polygons
world2 <- worldMap %>% sf::st_difference(polygon)
#> Warning: attribute variables are assumed to be spatially constant throughout all
#> geometries
# Transform
sf.ocean_pacific <- world2 %>% sf::st_transform(crs = target_crs)

# Define the file names and titles - Switch
df_div <- fread(wd_in)

# Choose index - SWITCH
i <- 1
index_div_model <- c("domain_Bacteria", "domain_Archaea")
# index_obs <- c("d_Bacteria_Richness_5000", "d__Archaea_Richness_5000")
titles <- c("Modeled bacterial richness", "Modeled archaeal richness")

# Extract title of the current clade
clade <- titles[i]


# Convert to spatraster object
df_div_sub <- filter(df_div, clade == index_div_model[i])
dat <- terra::rast(df_div_sub[, c("x", "y", "mean_richness")])
# Set reference system
terra::crs(dat) <- 'epsg:4326'

# Uncertainty stippling  
mean_min_threshold <- 5      # For richness
cv_abs_threshold <- 20       # Absolute CV threshold (%)

# Compute CV
df_div <- df_div %>%
  dplyr::mutate(
    mean_richness = ifelse(mean_richness < mean_min_threshold, NA, mean_richness),
    CV = 100 * sd_richness / mean_richness,
    CV = ifelse(is.na(mean_richness), NA, CV)
  )

# Calculate 70th percentile of CV for each layer
cv_percentile <- df_div %>%
  dplyr::group_by(clade) %>%
  dplyr::summarise(cv_70th = quantile(CV, probs = 0.7, na.rm = TRUE))

# Merge
df_div <- dplyr::left_join(df_div, cv_percentile, by = "clade")

# Compute final stippling threshold
df_div <- df_div %>%
  dplyr::mutate(
    final_cv_threshold = pmax(cv_abs_threshold, cv_70th, na.rm = TRUE),
    stipple = as.integer(!is.na(CV) & CV >= final_cv_threshold & !is.na(mean_richness))
  )

# Only keep stippling points
uncertainty <- df_div %>%
  dplyr::filter(clade == index_div_model[i], stipple == 1)

# Convert to spatraster object
uncertainty <- sf::st_as_sf(uncertainty, coords = c("x", "y"), crs = 4326)
  
# Extract longitude and latitude from the geometry column
uncertainty_coords <- cbind(
  st_drop_geometry(uncertainty),  # Drop the geometry column but keep other data
  st_coordinates(uncertainty)    # Extract coordinates as separate columns
)

# Rename the coordinate columns for clarity
colnames(uncertainty_coords)[ncol(uncertainty_coords) - 1:0] <- c("longitude", "latitude")

# View the updated data
head(uncertainty_coords)

# Create hexbin object for stippling
hb  <- erode(hexbin(uncertainty_coords$longitude, uncertainty_coords$latitude, xbins = 50))
df2 <- as.data.frame(hcell2xy(hb))

# Convert stippling points to the Pacific-centered CRS
stipple_points_sf <- st_as_sf(df2, coords = c("x", "y"), crs = 4326)
stipple_points_sf <- st_transform(stipple_points_sf, crs = target_crs)

  
# Get the minimum and maximum values of the raster
min_val <- terra::global(dat, fun = "min", na.rm = TRUE)[1, 1]
max_val <- terra::global(dat, fun = "max", na.rm = TRUE)[1, 1]
  
# Calculate quantiles to define contour intervals - To have contour lines where most of the data is centered!
quantiles <- quantile(values(dat), probs = seq(0, 1, by = 0.1), na.rm = TRUE)
  
# Define the contour step intervals based on quantiles
contour_intervals <- unique(round(quantiles, 0))
  
# Ensure the intervals cover the entire range
contour_intervals <- sort(unique(c(min_val, contour_intervals, max_val)))
  
# Create annotations for plot
model_stats <- get(load(fnames_model_stats[i]))
predictive_power <- round(unique(model_stats[[1]]$Predictive_Power), 2)
plot_title <- paste0(clade)
  
# Example coordinates
text_coords <- data.frame(lon = 50, lat = -80)

# Convert coordinates to sf object and transform CRS
text_sf <- sf::st_as_sf(text_coords, coords = c("lon", "lat"), crs = 4326)
text_sf <- sf::st_transform(text_sf, st_crs(sf.ocean_pacific))
text_coords_transformed <- sf::st_coordinates(text_sf)
  
# Generate world map
plot <- ggplot(sf.ocean_pacific) +
      
  # Base map with ocean regions
  geom_sf(fill = "grey", size = 0.1) +
      
    # Contour plot with filled contours
    tidyterra::geom_spatraster_contour_filled(
      data = dat, 
      show.legend = TRUE, 
      breaks = contour_intervals
    ) +
      
    # Add stippling upper 70th percentile bound
    geom_sf(data = stipple_points_sf, fill = "white", alpha = 1, size = 0.5, shape = 21, stroke = 0.1) +
      
    # Annotate predictive power (e.g., R^2)
    annotate(
      "text", 
      x = text_coords_transformed[1, "X"], 
      y = text_coords_transformed[1, "Y"],
      label = bquote("R"^2 * ": " * .(predictive_power)), 
      size = 8 / .pt, 
      color = "red"
    ) +
      
    # Add title and axis labels
    labs(
      title = plot_title,
      size = "Within sample richness",
      fill = "SDM richness"
    ) +
    # Customize theme for small figure
    theme_minimal(base_size = 8) +
    theme(
      text = element_text(size = 8),
      legend.position = "right",
      legend.key.size = unit(0.2, 'cm'),
      legend.key.height = unit(0.2, 'cm'),
      legend.key.width = unit(0.2, 'cm'),
      legend.title = element_text(size = 8, hjust = 0.5),
      legend.text = element_text(size = 8),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      plot.title = element_text(size = 8, face = "bold", hjust = 0.5)
    ) +
    xlab("") +
    ylab("") +
    guides(
      size = guide_legend(position = "bottom"),
      fill = guide_legend(position = "left")
    )

# Calculate 75th percentile
contour_75 <- quantile(values(dat), probs = 0.75, na.rm = TRUE)

# Add white contour line to the plot
plot_withContour <- plot +
  tidyterra::geom_spatraster_contour(
    data = dat,
    breaks = contour_75,
    color = "white",
    size = 0.35,
    linetype = "solid"
)

# Print the plot
print(plot_withContour)

# Save as high resolution image
ggsave(
  filename = paste0(wd_out, "GlobalMap_annual_richness_", index_div_model[i], "_v2.svg"),
  plot = plot_withContour,
  width = 3, height = 3, 
  dpi = 300
)

## ====================================================================================================
## END OF SCRIPT
## ====================================================================================================
