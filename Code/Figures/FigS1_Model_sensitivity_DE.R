## ====================================================================================================
## Supplementary Figure S1: Model Sensitivity Analysis - Coastal Distance and Basin-Stratified 
## Subsampling Uncertainty in Global Prokaryotic Richness Predictions
##
## This R script creates publication-quality visualizations comparing two types of model uncertainty
## for prokaryotic richness predictions: (1) uncertainty from coastal distance filtering and 
## (2) uncertainty from basin-stratified subsampling sensitivity analysis. Both analyses generate
## Pacific-centered world maps showing coefficient of variation patterns to assess prediction
## reliability under different methodological assumptions.
##
## Author:       Dominic Eriksson
## Date:         27th of February 2026
## Affiliation:  Environmental Physics Group UP, ETH Zürich, Switzerland
## Contact:      deriksson@ethz.ch
##
## Input files:  
##   - CSV with coastal distance filtering uncertainty metrics (no marginal seas)
##   - CSV with basin-stratified subsampling uncertainty metrics (no marginal seas)
##   - Ocean basin shapefiles for cartographic context
##
## Output files: 
##   - Two high-resolution SVG maps for supplementary publication
##   - Panel A: Coastal distance filtering coefficient of variation
##   - Panel B: Basin-stratified subsampling coefficient of variation
##
## Strategy:
##   The script employs advanced cartographic techniques including Pacific-centered projections
##   and spatraster visualization to create publication-ready comparative maps. Both analyses
##   use white-to-navy color scales for coefficient of variation, enabling direct comparison
##   of prediction sensitivity across different methodological approaches. This comparison
##   reveals spatial patterns in model robustness and identifies regions most affected by
##   sampling bias and coastal proximity assumptions.
##
## Required R packages (tested versions):
##   - raster        3.6.26
##   - ggplot2       3.4.2
##   - dplyr         1.1.2
##   - ggpubr        0.6.0
##   - betapart      1.5.6
##   - sf            1.0.12
##   - rnaturalearth 0.3.3
##   - terra         1.7.29
##   - tidyterra     0.4.0
## ====================================================================================================

# Clear workspace
rm(list = ls())

## Load libraries
library(raster)        # For raster data manipulation
library(ggplot2)       # For data visualization
library(dplyr)         # For data manipulation
library(ggpubr)        # For publication-ready plots
library(betapart)      # For biodiversity analysis
library(sf)            # For spatial data handling
library(rnaturalearth) # For world map data
library(terra)         # For modern raster handling
library(tidyterra)     # For terra integration with ggplot2

# Directories
wd_out <- "Code/Figures/Output/"

# Create output directory if it doesn't exist
if (!dir.exists(wd_out)) {
    dir.create(wd_out, recursive = TRUE)
}

## ====================================================================================================
## PANEL A: Coastal Distance Filtering Uncertainty Analysis
## ====================================================================================================

# Load coastal influence uncertainty data
wd_in_coastal <- "Code/6_Coastal_influence/4b_Output/Model_uncertainty_across_distance_groups_no_marginal_seas_prokaryotes.csv"

### ====================================================================================================================
### Panel A: Coastal distance filtering model uncertainty
### ====================================================================================================================


## Required data & data formatting
# Read shapefile with geometries of ocean basins
sf.ocean <- st_read("Data/Ocean_regions_shapefile/GOaS_v1_20211214/goas_v01.shp")
# Pacific centered world map: Source: https://stackoverflow.com/questions/56146735/visual-bug-when-changing-robinson-projections-central-meridian-with-ggplot2
worldMap <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_make_valid()
# Set projection pacific centered
lon_pacific <- '200' 
target_crs <- st_crs( paste0("+proj=eqc +x_0=0 +y_0=0 +lat_0=0 +lon_0=", lon_pacific) )

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
world2 <- worldMap %>% st_difference(polygon)
#> Warning: attribute variables are assumed to be spatially constant throughout all
#> geometries
# Transform
sf.ocean_pacific <- world2 %>% st_transform(crs = target_crs)


# Convert coastal data to spatraster object
dat_coastal <- read.csv(wd_in_coastal)
dat <- data.frame(
  x = dat$decimalLongitude,
  y = dat$decimalLatitude,
  cv_richness  = dat$cv_richness
)

# Convert to spatraster object
dat <- terra::rast(dat)
# Set reference system
crs(dat) <- 'epsg:4326'

# Plot world map - Publication ready version
gg_plot <-  ggplot(sf.ocean_pacific) + 
    
    geom_sf(fill = "grey90", color = "grey70", linewidth = 0.2) +
    
    # Plot continuous raster with white to red color scale
    tidyterra::geom_spatraster(data = dat) + 

    # Add contour lines at CV = 10 and 20 to highlight low uncertainty areas
    # tidyterra::geom_spatraster_contour(
    #   data = dat, 
    #   breaks = c(10, 20),
    #   color = "black", 
    #   linewidth = 0.6,
    #   alpha = 0.8
    # ) +

    # White to blue color scale for model uncertainty (alternative to red)
    scale_fill_gradient(
      low = "white", 
      high = "navy",
      na.value = "transparent",
      name = "CV (%)",
      guide = guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = 12,
        barheight = 0.8,
        frame.colour = "black",
        ticks.colour = "black"
      )
    ) +

    # Publication quality theme
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10),
      legend.margin = margin(t = 10, b = 5),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )

# Print Panel A plot
print(gg_plot_coastal)

# Save Panel A as SVG
ggsave(
    filename = paste0(wd_out, "FigS1_PanelA_CoastalDistance_ModelUncertainty.svg"), 
    plot = gg_plot_coastal, 
    width = 7, height = 4.5,
    units = "in",
    dpi = 300
)

## ====================================================================================================
## PANEL B: Basin-Stratified Subsampling Uncertainty Analysis
## ====================================================================================================

# Load basin-stratified subsampling uncertainty data
wd_in_thinning <- "Code/7_Ocean_thinning/6_Output/prokaryote_richness_cv_no_marginal_seas.csv"

### ====================================================================================================================
### Panel B: Basin-stratified subsampling model uncertainty
### ====================================================================================================================

# Convert basin thinning data to spatraster object (reusing same projection setup)
dat_thinning <- read.csv(wd_in_thinning)
dat_thinning <- data.frame(
  x = dat_thinning$decimalLongitude,
  y = dat_thinning$decimalLatitude,
  cv_richness  = dat_thinning$cv_richness
)

# Convert to spatraster object
dat_thinning <- terra::rast(dat_thinning)
# Set reference system
crs(dat_thinning) <- 'epsg:4326'

# Plot Panel B - Basin-stratified subsampling uncertainty
gg_plot_thinning <-  ggplot(sf.ocean_pacific) + 
    geom_sf(fill = "grey90", color = "grey70", linewidth = 0.2) +
    
    # Plot continuous raster
    tidyterra::geom_spatraster(data = dat_thinning) + 

    # White to navy color scale for model uncertainty
    scale_fill_gradient(
      low = "white", 
      high = "navy",
      na.value = "transparent",
      name = "CV (%)",
      guide = guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = 12,
        barheight = 0.8,
        frame.colour = "black",
        ticks.colour = "black"
      )
    ) +

    # Publication quality theme
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10),
      legend.margin = margin(t = 10, b = 5),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )

# Print Panel B plot
print(gg_plot_thinning)

# Save Panel B as SVG
ggsave(
    filename = paste0(wd_out, "FigS1_PanelB_BasinThinning_ModelUncertainty.svg"), 
    plot = gg_plot_thinning, 
    width = 7, height = 4.5,
    units = "in",
    dpi = 300
)

# **Supplementary Figure S1. Model sensitivity analysis comparing two sources of uncertainty in 
# global prokaryotic richness predictions.** 
# Panel A shows coefficient of variation (CV, %) from coastal distance filtering sensitivity, 
# representing uncertainty arising from different minimum distance-to-coast thresholds applied 
# to sample selection. Panel B shows CV from basin-stratified subsampling sensitivity across 
# 1,100 scenarios with sample size caps ranging from 10-350 samples per ocean basin across 
# 100 iterations each. Both panels use Pacific-centered projections with white-to-navy color 
# scales where higher CV values indicate greater prediction sensitivity to methodological choices. 
# Continental areas are shown in grey. Marginal seas were excluded to focus on open ocean 
# dynamics. This comparison reveals spatial patterns in model robustness and identifies regions 
# most affected by sampling bias and coastal proximity assumptions.

## ====================================================================================================
## END OF SCRIPT
## ====================================================================================================
