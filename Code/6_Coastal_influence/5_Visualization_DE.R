## ====================================================================================================
## This R script creates visualizations of model uncertainty for prokaryotic
## richness predictions across coastal to open ocean gradients. It generates Pacific-centered world
## maps showing coefficient of variation patterns, highlighting regions where coastal distance
## filtering creates the greatest prediction uncertainty.
##
## Author:       Dominic Eriksson
## Date:         26th of February 2026
## Affiliation:  Environmental Physics Group UP, ETH Zürich, Switzerland
## Contact:      deriksson@ethz.ch
##
## Input files:  
##   - CSV with model uncertainty metrics (no marginal seas)
##   - Ocean basin shapefiles for cartographic context
##
## Output files: 
##   - High-resolution SVG map for publication
##   - Coefficient of variation visualization across global oceans
##
## Strategy:
##   The script employs advanced cartographic techniques including Pacific-centered projections
##   and spatraster visualization to create publication-ready maps. The coefficient of variation
##   from coastal distance filtering is displayed using a white-to-navy color scale, enabling
##   clear identification of regions where model predictions are most sensitive to coastal
##   proximity assumptions.
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
wd_in <- "Code/6_Coastal_influence/4b_Output/Model_uncertainty_across_distance_groups_no_marginal_seas_prokaryotes.csv"
wd_out <- "Code/6_Coastal_influence/5_Output/"

### ====================================================================================================================
### MAP - Showing variance of model outputs regarding the inclusion of coastal samples
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


# Convert to spatraster object
dat <- read.csv(wd_in)
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

# Print plot
print(gg_plot)

# Save as SVG for publication (one column of two-row multiplot)
ggsave(
    filename = paste0(wd_out, "Map_coefficientOfVariation_richnessProkaryotes_CoastToOpenOcean_publication_v2.svg"), 
    plot = gg_plot, 
    width = 7, height = 4.5,
    units = "in",
    dpi = 300
)

## ====================================================================================================
## END OF SCRIPT
## ====================================================================================================