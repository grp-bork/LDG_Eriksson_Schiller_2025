
##====================================================================================================
## Supplementary Figure S13 Panel C: Global Maps and Correlation Analysis of Genus-Level Richness
##
## This R script generates global distribution maps of annual mean richness for Prochlorococcus and 
## Pelagibacter based on species distribution modeling ensemble predictions. The analysis includes 
## sophisticated uncertainty quantification using coefficient of variation (CV) thresholding, 
## Pacific-centered cartographic projections, and hexagonal stippling for model uncertainty visualization.
## Additionally, the script performs correlation analysis between total prokaryotic richness and 
## genus-specific richness patterns using Pearson correlation with hexagonal binning scatterplots.
##
## Author:      Dominic Eriksson
##              Environmental Physics Group, UP
##              ETH Zürich, Switzerland
## Contact:     deriksson@ethz.ch
## Date:        February 27th, 2026
## Affiliation: ETH Zürich, Environmental Physics Group, UP
##
## Input files:
## - /net/sea/work/deriksson/Projects/Prokaryotes_LDG_v3/Code/4_Compute_ensembles_uncertainties/2c_Output/marginal_seas_filtered_genus_Prochlorococcus_A.csv
## - /net/sea/work/deriksson/Projects/Prokaryotes_LDG_v3/Code/4_Compute_ensembles_uncertainties/2c_Output/marginal_seas_filtered_genus_Pelagibacter.csv
## - /net/sea/work/deriksson/Projects/Prokaryotes_LDG_v3/Code/4_Compute_ensembles_uncertainties/2_Output/Global_annual_richness_ensemble_mean_richness_v2.csv
## - /net/sea/work/deriksson/Projects/Prokaryotes_LDG_v3/Code/2_Extract_model_outputs/6_Output/GenusLevel/genus_Prochlorococcus_A_Richness_5000.RData
## - /net/sea/work/deriksson/Projects/Prokaryotes_LDG_v3/Code/2_Extract_model_outputs/6_Output/GenusLevel/genus_Pelagibacter_Richness_5000.RData
## - /net/sea/work/deriksson/Projects/Notion_DE_reviewNatureComms/Data/Ocean_regions_shapefile/GOaS_v1_20211214/goas_v01.shp
##
## Output files:
## - GlobalMap_annual_richness_Prochlorococcus_richness.svg (Global Prochlorococcus distribution map)
## - GlobalMap_annual_richness_Pelagibacter_richness.svg (Global Pelagibacter distribution map)
## - Correlation_Prokaryotes_Prochlorococcus.svg (Prokaryotes vs Prochlorococcus correlation plot)
## - Correlation_Prokaryotes_Pelagibacter.svg (Prokaryotes vs Pelagibacter correlation plot)
##
## Strategy:
## 1. Load required spatial analysis, visualization, and statistical libraries
## 2. Set up Pacific-centered map projection (lon_0 = 200°) for optimal marine visualization
## 3. Load genus-specific richness data and model performance statistics
## 4. Calculate coefficient of variation (CV) for uncertainty quantification
## 5. Apply dual-threshold uncertainty detection (absolute 20% + 70th percentile)
## 6. Exclude low-richness regions (minimum threshold = 5) to avoid statistical artifacts
## 7. Create hexagonal stippling masks for high-uncertainty regions
## 8. Generate contour-filled global maps with 75th percentile highlighting
## 9. Perform Pearson correlation analysis between prokaryotic and genus-specific richness
## 10. Create hexagonal binning scatterplots with statistical annotations and regression lines
##
## Required R packages:
## - ggplot2 (3.3.6): Advanced data visualization and mapping framework
## - dplyr (1.0.9): Data manipulation and transformation operations
## - sf (1.0.7): Simple features spatial data handling and coordinate transformations
## - tidyterra (0.3.1): Tidy integration of terra raster objects with ggplot2
## - cowplot (1.1.1): Advanced plot arrangement and publication layouts
## - terra (1.6.7): Modern raster data processing and spatial analysis
## - hexbin (1.28.2): Hexagonal binning for spatial point pattern analysis
## - data.table (1.14.2): High-performance data manipulation and CSV processing
## - raster (3.5.15): Traditional raster data processing and coordinate operations
## - sp (1.4.7): Legacy spatial data classes and coordinate system handling
## - rnaturalearth (0.1.0): Natural Earth map data access for cartographic visualization
##====================================================================================================

# Key Methodology:
# 1. **Coefficient of Variation (CV) Calculation**:
#    - CV is calculated for each spatial unit as:
#      CV = (Standard Deviation / Mean) * 100
#    - NaN values (e.g., from division by zero or missing data) are replaced with 0 to ensure
#      a clean and stable dataset.

# 2. **Threshold Determination**:
#    - Two thresholds are combined to determine regions of high uncertainty:
#      a) An absolute CV threshold (`cv_abs_threshold`, default = 20%), ensures that only
#         regions with CV above a minimum variability level are considered.
#      b) A relative threshold, calculated as the 70th percentile of CV values within each
#         layer, dynamically adapts to the distribution of variability in the data.
#    - The final threshold for each layer is the maximum of these two values, ensuring
#      both minimum and relative variability criteria are met.

# 3. **Exclusion of Low-Richness Regions**:
#    - To avoid inflating CV in regions where mean richness is low, a minimum mean richness
#      threshold (`mean_min_threshold`, default = 5 for richness and 0.5 for shannon) is applied. Only regions where the mean
#      richness exceeds this threshold are included in the analysis, ensuring the CV is not
#      disproportionately high due to small mean values, which would lead to misleading uncertainty interpretations.

# 4. **Stippling Mask Creation**:
#    - For each spatial layer, a binary stippling mask is created where:
#      a) CV exceeds the final threshold, and
#      b) Mean richness is greater than or equal to the minimum threshold.
#    - These stippling masks are stacked to provide a multi-layer representation of areas
#      with significant variability across spatial and temporal dimensions.

# Rationale:
# This method ensures that the uncertainty index highlights areas of genuinely high
# variability while avoiding overrepresentation of noise or artifacts caused by
# low richness values. The combination of absolute, relative, and minimum richness
# thresholds provides a balanced and adaptive framework for analyzing spatial
# uncertainty in ecological datasets.


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
library(raster)
library(sp)

# Working directories
wd_in_pro <- "/net/sea/work/deriksson/Projects/Prokaryotes_LDG_v3/Code/4_Compute_ensembles_uncertainties/2c_Output/marginal_seas_filtered_genus_Prochlorococcus_A.csv"
wd_in_pelagi <- "/net/sea/work/deriksson/Projects/Prokaryotes_LDG_v3/Code/4_Compute_ensembles_uncertainties/2c_Output/marginal_seas_filtered_genus_Pelagibacter.csv"
wd_out <- "/net/sea/work/deriksson/Projects/Prokaryotes_LDG_v3/Revisions_Cell_Host_and_Microbes/Code/1_Output/"

# Filenames for model stats
fnames_model_stats <- c(
    "/net/sea/work/deriksson/Projects/Prokaryotes_LDG_v3/Code/2_Extract_model_outputs/6_Output/GenusLevel/genus_Prochlorococcus_A_Richness_5000.RData",
    "/net/sea/work/deriksson/Projects/Prokaryotes_LDG_v3/Code/2_Extract_model_outputs/6_Output/GenusLevel/genus_Pelagibacter_Richness_5000.RData"
)

## Required data & data formatting
# Read shapefile with geometries of ocean basins - Note: Specific shapefile needed for Pacific-centered map
sf.ocean <- sf::st_read("/net/sea/work/deriksson/Projects/Notion_DE_reviewNatureComms/Data/Ocean_regions_shapefile/GOaS_v1_20211214/goas_v01.shp")
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
df_div <- fread(wd_in_pro) # index i = 1
df_div <- fread(wd_in_pelagi)# index i = 2

# Choose index - SWITCH
i <- 1 # change index accordingly: 1 = Prochlorococcus, 2 = Pelagibacter
# index_div_model <- c("domain_Bacteria", "domain_Archaea")
# index_obs <- c("d_Bacteria_Richness_5000", "d__Archaea_Richness_5000")
titles <- c("Prochlorococcus richness", "Pelagibacter richness")

# Subset input data from models
clade <- titles[i]
  
# Convert to spatraster object
dat <- terra::rast(df_div[, c("x", "y", "mean_annual")])
# Set reference system
terra::crs(dat) <- 'epsg:4326'

# Uncertainty stippling  
mean_min_threshold <- 5      # For richness
cv_abs_threshold <- 20       # Absolute CV threshold (%)

df_div <- df_div %>%
  dplyr::mutate(
    mean_richness = ifelse(mean_annual < mean_min_threshold, NA, mean_annual),
    CV = 100 * sd_annual / mean_annual,
    CV = ifelse(is.na(mean_annual), NA, CV)
  )

cv_percentile <- df_div %>%
  dplyr::group_by(clade) %>%
  dplyr::summarise(cv_70th = quantile(CV, probs = 0.7, na.rm = TRUE))

# Merge
df_div <- dplyr::left_join(df_div, cv_percentile, by = "clade")

df_div <- df_div %>%
  dplyr::mutate(
    final_cv_threshold = pmax(cv_abs_threshold, cv_70th, na.rm = TRUE),
    stipple = as.integer(!is.na(CV) & CV >= final_cv_threshold & !is.na(mean_annual))
)

# Only keep stippling points
uncertainty <- df_div %>%
  dplyr::filter(stipple == 1)

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

# # Replace 166 with 150 for better interpretation --> Only for prochlorococcus
# contour_intervals[contour_intervals == 166] <- 150
  


# Create annotations for plot
model_stats <- get(load(fnames_model_stats[i]))
predictive_power <- round(unique(model_stats[[1]]$Predictive_Power), 2)
# number_samples <- nrow(df_sub)
plot_title <- paste0(clade)
  
# Example coordinates
text_coords <- data.frame(
  lon = 50, lat = -80
)
  
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
    label = paste0("R² = ", predictive_power),
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
  filename = paste0(wd_out, "GlobalMap_annual_richness_", gsub(" ", "_", clade), ".svg"),
  plot = plot_withContour,
  width = 3, height = 3, 
  dpi = 300
)

### ==============================================================================
### PUBLICATION-READY FIGURE CAPTION
### ==============================================================================

# Figure Caption:
# Global distribution of annual mean richness for [Prochlorococcus/Pelagibacter] based on 
# species distribution modeling ensemble predictions. The color gradient represents predicted 
# richness values derived from an ensemble of species distribution models, with darker colors 
# indicating higher richness. White stippling marks regions of high model uncertainty, defined 
# as areas where the coefficient of variation (CV) exceeds both an absolute threshold (20%) 
# and the 70th percentile of CV values within the dataset, while maintaining a minimum mean 
# richness threshold (≥5) to exclude regions with artificially inflated variability due to 
# low mean values. The white contour line delineates the 75th percentile of richness values, 
# highlighting regions of exceptionally high diversity. The R² value indicates the predictive 
# performance of the species distribution models based on cross-validation. The map uses a 
# Pacific-centered equidistant cylindrical projection (lon_0 = 200°) to optimize visualization 
# of global oceanic patterns. Model predictions are based on ensemble means from 30 model 
# iterations, incorporating environmental predictors and accounting for spatial autocorrelation 
# in marine prokaryotic communities. Uncertainty analysis ensures that only regions with 
# genuine ecological variability (rather than statistical artifacts) are highlighted as 
# uncertain, providing a robust assessment of model reliability across different oceanic regions.



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
library(raster)
library(sp)

# Working directories
wd_in_pro <- "/net/sea/work/deriksson/Projects/Prokaryotes_LDG_v3/Code/4_Compute_ensembles_uncertainties/2c_Output/marginal_seas_filtered_genus_Prochlorococcus_A.csv"
wd_in_pelagi <- "/net/sea/work/deriksson/Projects/Prokaryotes_LDG_v3/Code/4_Compute_ensembles_uncertainties/2c_Output/marginal_seas_filtered_genus_Pelagibacter.csv"
wd_in_prokaryotes <- "/net/sea/work/deriksson/Projects/Prokaryotes_LDG_v3/Code/4_Compute_ensembles_uncertainties/2_Output/Global_annual_richness_ensemble_mean_richness_v2.csv"
wd_out <- "/net/sea/work/deriksson/Projects/Prokaryotes_LDG_v3/Revisions_Cell_Host_and_Microbes/Code/5_Output/"

# Load data
df_prokaryotes <- fread(wd_in_prokaryotes)
# Subset for clade == prokaryotes
df_prokaryotes <- df_prokaryotes %>% filter(clade == "prokaryotes")

df_pro <- fread(wd_in_pro)
df_pro <- df_pro %>% dplyr::select(x, y, mean_annual)

df_pelagi <- fread(wd_in_pelagi)
df_pelagi <- df_pelagi %>% dplyr::select(x, y, mean_annual)

## Correlation Analysis: Prokaryotes vs Prochlorococcus ----------------------

# Merge prokaryotes and prochlorococcus data
merged_pro <- df_prokaryotes %>%
  dplyr::left_join(df_pro, by = c("x", "y"), suffix = c("_prokaryotes", "_prochlorococcus")) %>%
  dplyr::filter(!is.na(mean_richness) & !is.na(mean_annual)) %>%
  dplyr::select(x, y, mean_richness, mean_annual)

# Calculate Pearson correlation
cor_pro <- cor.test(merged_pro$mean_richness, merged_pro$mean_annual, method = "pearson")

# Print results
cat("=== PROKARYOTES vs PROCHLOROCOCCUS ===\n")
cat("Number of valid data points:", nrow(merged_pro), "\n")
cat("Pearson's r:", round(cor_pro$estimate, 4), "\n")
cat("p-value:", format(cor_pro$p.value, scientific = TRUE, digits = 3), "\n")
cat("95% CI:", round(cor_pro$conf.int[1], 4), "to", round(cor_pro$conf.int[2], 4), "\n")
cat("Significance:", ifelse(cor_pro$p.value < 0.001, "***", 
                           ifelse(cor_pro$p.value < 0.01, "**",
                                  ifelse(cor_pro$p.value < 0.05, "*", "ns"))), "\n\n")

# Create scatterplot with hexbins
p_pro <- ggplot(merged_pro, aes(x = mean_richness, y = mean_annual)) +
  geom_hex(bins = 30) +
  scale_fill_viridis_c(
    option = "C",
    name = "Count"
  ) +
  geom_smooth(method = "lm", color = "blue", se = FALSE, size = 1) +
  annotate("text", 
           x = -Inf, y = Inf, 
           hjust = -0.1, vjust = 1.2,
           label = paste0("r = ", round(cor_pro$estimate, 3), "\n",
                         "p ", ifelse(cor_pro$p.value < 0.001, "< 0.001***", 
                                     ifelse(cor_pro$p.value < 0.01, "< 0.01**",
                                            ifelse(cor_pro$p.value < 0.05, "< 0.05*", 
                                                   paste("=", round(cor_pro$p.value, 3)))))),
           color = "red", size = 8 / .pt, fontface = "bold") +
  labs(
    title = "Prokaryote vs Prochlorococcus richness",
    x = "Prokaryote richness",
    y = "Prochlorococcus richness"
  ) +
  theme_minimal(base_size = 8) +
  theme(
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8, face = "bold"),
    plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.key.size = unit(0.3, 'cm'),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8)
  )

# Display and save plot
print(p_pro)
ggsave(paste0(wd_out, "Correlation_Prokaryotes_Prochlorococcus.svg"), 
       p_pro, width = 3, height = 3, dpi = 300)

## Correlation Analysis: Prokaryotes vs Pelagibacter -------------------------

# Merge prokaryotes and pelagibacter data
merged_pelagi <- df_prokaryotes %>%
  dplyr::left_join(df_pelagi, by = c("x", "y"), suffix = c("_prokaryotes", "_pelagibacter")) %>%
  dplyr::filter(!is.na(mean_richness) & !is.na(mean_annual)) %>%
  dplyr::select(x, y, mean_richness, mean_annual)

# Calculate Pearson correlation
cor_pelagi <- cor.test(merged_pelagi$mean_richness, merged_pelagi$mean_annual, method = "pearson")

# Print results
cat("=== PROKARYOTES vs PELAGIBACTER ===\n")
cat("Number of valid data points:", nrow(merged_pelagi), "\n")
cat("Pearson's r:", round(cor_pelagi$estimate, 4), "\n")
cat("p-value:", format(cor_pelagi$p.value, scientific = TRUE, digits = 3), "\n")
cat("95% CI:", round(cor_pelagi$conf.int[1], 4), "to", round(cor_pelagi$conf.int[2], 4), "\n")
cat("Significance:", ifelse(cor_pelagi$p.value < 0.001, "***", 
                           ifelse(cor_pelagi$p.value < 0.01, "**",
                                  ifelse(cor_pelagi$p.value < 0.05, "*", "ns"))), "\n\n")

# Create scatterplot with hexbins
p_pelagi <- ggplot(merged_pelagi, aes(x = mean_richness, y = mean_annual)) +
  geom_hex(bins = 30) +
  scale_fill_viridis_c(
    option = "C",
    name = "Count"
  ) +
  geom_smooth(method = "lm", color = "blue", se = FALSE, size = 1) +
  annotate("text", 
           x = -Inf, y = Inf, 
           hjust = -0.1, vjust = 1.2,
           label = paste0("r = ", round(cor_pelagi$estimate, 3), "\n",
                         "p ", ifelse(cor_pelagi$p.value < 0.001, "< 0.001***", 
                                     ifelse(cor_pelagi$p.value < 0.01, "< 0.01**",
                                            ifelse(cor_pelagi$p.value < 0.05, "< 0.05*", 
                                                   paste("=", round(cor_pelagi$p.value, 3)))))),
           color = "red", size = 8 / .pt, fontface = "bold") +
  labs(
    title = "Prokaryote vs Pelagibacter richness",
    x = "Prokaryote richness",
    y = "Pelagibacter richness"
  ) +
  theme_minimal(base_size = 8) +
  theme(
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8, face = "bold"),
    plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.key.size = unit(0.3, 'cm'),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8)
  )

# Display and save plot
print(p_pelagi)
ggsave(paste0(wd_out, "Correlation_Prokaryotes_Pelagibacter.svg"), 
       p_pelagi, width = 3, height = 3, dpi = 300)

### ====================================================================================================================
### END OF SCRIPT
### ====================================================================================================================