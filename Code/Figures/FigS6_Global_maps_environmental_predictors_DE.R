##====================================================================================================
## Supplementary Figure S6: Global Maps of Environmental Predictors
##
## This R script creates global maps of annual mean environmental predictors used in prokaryotic diversity modeling.
## The generated maps are Pacific-centered and visualize key oceanographic variables including nutrient concentrations 
## (nitrate, phosphate, silicate), physical parameters (temperature, salinity, mixed layer depth), biological 
## variables (chlorophyll-a, primary production), and dynamic features (eddy kinetic energy, Lyapunov exponents).
## Each environmental layer is displayed with appropriate color scaling and geographic projections optimized 
## for marine ecosystem visualization.
##
## Author:      Dominic Eriksson
##              Environmental Physics Group, UP
##              ETH Zürich, Switzerland
## Contact:     deriksson@ethz.ch
## Date:        February 27th, 2026
## Affiliation: ETH Zürich, Environmental Physics Group, UP
##
## Input files:
## - ../../Data/Ocean_regions_shapefile/GOaS_v1_20211214/goas_v01.shp (Ocean basins shapefile)
## - ../../Data/Environmental_parameters/Surface_predictors/*.nc (NetCDF environmental parameter files)
##
## Output files:
## - ./Output/Annual_surface_*_Pacific_centered_v4.png (Individual environmental predictor maps)
## - ./Output/Annual_surface_predictors_Pacific_centered_v4.png (Combined multiplot of all predictors)
##
## Strategy:
## 1. Load required spatial analysis and visualization libraries
## 2. Set up Pacific-centered map projection and coordinate system
## 3. Load ocean basins shapefile and natural earth country boundaries
## 4. Create meridian polygon for Pacific-centered projection transformation
## 5. Transform geographic data to Pacific-centered coordinate reference system
## 6. Load and process NetCDF environmental parameter files
## 7. Calculate annual means for each environmental variable
## 8. Apply log10 transformation to chlorophyll-a for better visualization
## 9. Generate individual maps with turbo color palette and appropriate scaling
## 10. Combine all maps into multiplot and save publication-ready figures
##
## Required R packages:
## - ggplot2 (3.3.6): Advanced data visualization and plotting framework
## - dplyr (1.0.9): Data manipulation and transformation operations
## - ggpubr (0.4.0): Arranging multiple ggplot figures into publication grids
## - sf (1.0.7): Simple features spatial data handling and manipulation
## - rnaturalearth (0.1.0): Access to natural earth map data for visualization
## - terra (1.6.7): Modern raster data handling and spatial analysis
## - tidyterra (0.3.1): Tidy integration of terra raster objects with ggplot2
## - raster (3.5.15): Traditional raster data processing and analysis
## - cowplot (1.1.1): Advanced plot arrangement and publication-ready layouts
##====================================================================================================


### ====================================================================================================================
### Preparation
### ====================================================================================================================

# Clear system
rm(list = ls())

## Load libraries
library(ggplot2)        # For data visualization, including the maps with ocean regions and contour plots
library(dplyr)          # For data manipulation, filtering, and selecting relevant data
library(ggpubr)         # For arranging multiple ggplot figures into a grid
library(sf)             # For handling and manipulating spatial data in the form of simple features
library(rnaturalearth)  # For accessing map data for countries, oceans, and geographical regions
library(terra)          # For handling and processing raster data, essential for spatial analysis
library(tidyterra)      # For working with raster data in a tidy format that integrates well with ggplot2
library(raster)
library(cowplot)

### ====================================================================================================================
### Visualize global map of annual environmental predictors - Pacific centered
### ====================================================================================================================

## Required data & data formatting
# Read shapefile with geometries of ocean basins - Note: Specific shapefile needed for Pacific-centered map
sf.ocean <- sf::st_read("../../Data/Ocean_regions_shapefile/GOaS_v1_20211214/goas_v01.shp")
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


# Directory
wd_in <- "../../Data/Environmental_parameters/Surface_predictors/"
wd_out <- "./Output/"


# Adjust predictor names
column_name_map <- c(
  "climatology_A_0_50" = "Apparent Oxygen Utilization (micromoles/kg)",
  "climatology_A_CHLA_regridded" = "Chlorophyll-a (milligrams/m^3)",
  "climatology_A_PAR_regridded" = "Photosynthetically Active Radiation (Ein/m^2/day)",
  "climatology_M_0_0" = "Mixed Layer Depth (meters)",
  "climatology_Nstar_0_50" = "N* (micromoles/kg)",
  "climatology_O_0_50" = "Oxygen Saturation (%)",
  "log10_climatology_S_PP_regridded" = "Primary Production (milligrams of carbon/m^2/day)",
  "log10_climatology_TOT_POC_CMEMS" = "Total POC (mg C/m^3)",
  "climatology_o_0_50" = "Dissolved O2 (micromoles/kg)",
  "log10_climatology_eke_aviso" = "Eddy Kinetic Energy (cm^2/s^2)",
  "climatology_fsle_aviso_2001_2020" = "Lyap. Exp. (days^-1)",
  "climatology_i_0_50" = "Silicate (micromoles/kg)",
  "climatology_n_0_50" = "Nitrate (micromoles/kg)",
  "climatology_p_0_50" = "Phosphate (micromoles/kg)",
  "climatology_s_0_50" = "Salinity",
  "climatology_t_0_50" = "Temperature (°Celsius)"
)

# Get filenames to loop over
fnames <- list.files(wd_in, full.names = TRUE)

# Initialize a list to store the plots
plot_list <- list()

for(f in seq_along(fnames)){

    # Subset input data from models
    env_name <- gsub(".nc", "", basename(fnames[f]))
    r <- brick(fnames[f])

    # Calculate the annual mean
    r <- calc(r, fun = mean, na.rm = TRUE)

    # Name the layer
    names(r) <- env_name

    # Convert to spatraster object
    dat <- raster::as.data.frame(r, xy = TRUE)
    dat <- terra::rast(dat[, c("x", "y", env_name)])
    # Set reference system
    terra::crs(dat) <- 'epsg:4326'

    # Update env_name using column_name_map if not found skip to next iteration
    if (env_name %in% names(column_name_map)) {
        env_name_title <- column_name_map[[env_name]]
    } else {
        next
    }

    # If chlorophyll, visualize log10
    if( basename(fnames[f]) == "climatology_A_CHLA_regridded.nc"){
        dat <- log10(dat)}

    # Define min and max
    min_val <- min(values(dat), na.rm = TRUE)
    max_val <- max(values(dat), na.rm = TRUE)

    # Calculate percentiles
    percentiles <- quantile(values(dat), probs = seq(0, 1, 0.1), na.rm = TRUE)

    # Create breaks with more intervals in the range where most observations fall
    breaks <- c(
        seq(min_val, percentiles[4], length.out = 3),  # More breaks in the lower range
        seq(percentiles[4], percentiles[7], length.out = 4),  # More breaks in the middle range
        seq(percentiles[7], max_val, length.out = 3)  # More breaks in the upper range
    )

    # Define legend breaks and labels
    legend_breaks <- c(min_val, percentiles[5], max_val)
    legend_labels <- round(legend_breaks, 2)

    # Plot world map
    gg_plot <- ggplot() + 
        geom_sf(data = sf.ocean_pacific, fill = "grey") +
        
        # Plot raster data with turbo color palette
        tidyterra::geom_spatraster(data = dat) +
        scale_fill_viridis_c(option = "turbo", name = "Value") +
        
        # Add boundaries for continents and countries
        geom_sf(data = sf.ocean_pacific, fill = NA, color = "black", size = 0.5) +
        
        # Legend settings
        theme(
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            legend.position = "right",
            # legend.key.size = unit(10, 'cm'), 
            # legend.key.height = unit(0.3, 'cm'), 
            # legend.key.width = unit(0.3, 'cm'),
            legend.title = element_text(hjust = 0.5, size = unit(10, "pt")), 
            legend.text = element_text(size = unit(8, "pt")), 

            # Plot background color
            panel.background = element_rect(fill = "white"),
            plot.background = element_rect(fill = "white"),
            
            # Axis text and title size
            axis.text = element_text(size = unit(10, "pt")),
            axis.title = element_text(size = unit(10, "pt")),
            plot.title = element_text(hjust = 0.5, size = unit(12, "pt"))
        ) +
        
        # Axis labels
        xlab("") +
        ylab("") +
        ggtitle(column_name_map[[env_name]])
        
    # Add the plot to the list
    plot_list[[f]] <- gg_plot

    # Save as high resolution plot
    ggsave(
        filename = paste0(wd_out, "Annual_surface_", gsub( ".nc", "", basename(fnames[f]) ) ,"_Pacific_centered_v4.png"), 
        plot = gg_plot, 
        width = 18,
        height = 15,
        dpi = 300)
}

# Remove NULL elements from the list
plot_list <- plot_list[!sapply(plot_list, is.null)]

# Create a multiplot using cowplot
multiplot <- plot_grid(plotlist = plot_list, ncol = 2)

# Print the multiplot
print(multiplot)

# Save as high resolution plot
ggsave(
    filename = paste0(wd_out, "Annual_surface_predictors_Pacific_centered_v4.png"), 
    plot = multiplot, 
    width = 12,
    height = 12,
    dpi = 300)

# Figure caption:
# Global maps of annual mean environmental predictors in the surface ocean, centered 
# on the Pacific basin. Each panel represents a key environmental variable, 
# including nutrient concentrations (e.g., nitrate, phosphate), 
# chlorophyll-a, oxygen saturation, temperature, and other physical and 
# biological parameters. 

### ====================================================================================================================
### END OF SCRIPT
### ====================================================================================================================
