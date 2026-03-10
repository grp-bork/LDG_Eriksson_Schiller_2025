
## ====================================================================================================
## This R script generates comprehensive global maps showing annual richness patterns across 15 major
## prokaryotic clades (2 Archaeal and 13 Bacterial classes). The analysis creates Pacific-centered
## global maps with filled contour plots, uncertainty visualization through stippling, and model
## performance annotations (R² values). Each clade is visualized individually with species richness
## predictions displayed using quantile-based contour intervals and coefficient of variation thresholds
## for uncertainty representation.
##
## Author:       Dominic Eriksson
## Date:         27th of February 2026
## Affiliation:  Environmental Physics Group UP, ETH Zürich, Switzerland
## Contact:      deriksson@ethz.ch
##
## Input files:
##   - Global annual richness ensemble mean CSV with clade-specific richness predictions
##   - Ocean regions shapefile for Pacific-centered map projection
##   - Model statistics RData files containing predictive power metrics for each clade
##   - Natural Earth country boundaries for continental outlines
##
## Output files:
##   - Multi-panel SVG figure with 15 global richness maps (3×5 grid layout)
##   - Multi-panel PNG figure for presentations and web display
##   - Individual map plots for each prokaryotic clade with uncertainty visualization
##
## Strategy:
##   The script processes 15 prokaryotic clades spanning major Archaeal (Poseidoniia, Nitrososphaeria)
##   and Bacterial lineages (Acidimicrobiia, Bacteroidia, Cyanobacteriia, Alphaproteobacteria, etc.).
##   For each clade, quantile-based contour intervals optimize visualization of richness ranges while
##   coefficient of variation thresholds (70th percentile or 20% absolute threshold) identify high
##   uncertainty regions for white stippling. Pacific-centered projection optimizes marine ecosystem
##   visualization with R² model performance metrics integrated into each map.
##
## Required R packages (tested versions):
##   - ggplot2        3.4.2
##   - dplyr          1.1.2
##   - sf             1.0.12
##   - tidyterra      0.4.0
##   - cowplot        1.1.1
##   - terra          1.7.29
##   - hexbin         1.28.3
##   - viridis        0.6.3
##   - rnaturalearth  0.3.3
##   - stringr        1.5.0
## ====================================================================================================


### ====================================================================================================================
### Preparation
### ====================================================================================================================

# Clear workspace
rm(list = ls())

# Load libraries
library(ggplot2)        # For creating publication-quality plots and visualizations
library(dplyr)          # For data manipulation and filtering operations
library(sf)             # For handling spatial data and geographic projections
library(tidyterra)      # For working with terra raster objects in ggplot2
library(cowplot)        # For combining multiple plots into publication-ready figures
library(terra)          # For handling and processing raster data
library(hexbin)         # For hexagonal binning used in uncertainty visualization
library(viridis)        # For perceptually uniform color palettes
library(rnaturalearth)  # For accessing natural earth map data (loaded via sf functions)
library(stringr)        # For string manipulation and text wrapping functions

### ====================================================================================================================
### Data Loading and Configuration
### ====================================================================================================================

# Working directories
wd_in <- "../3_Compute_ensembles_uncertainties/2_Output/Global_annual_richness_ensemble_mean_richness.csv"
wd_out <- "./Output/"

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


# Load data
df_rich <- read.csv(wd_in)

# Define taxa of interest
unique(df_rich$clade)
# Clades of interest
clade <- c(
  "class_Poseidoniia",
  "class_Nitrososphaeria",
  "class_Acidimicrobiia",
  "class_Bacteroidia",
  "class_Chlamydiia",
  "class_Cyanobacteriia",
  "class_UBA1144",
  "class_Fibrobacteria",
  "class_Marinisomatia",
  "class_UBA4151",
  "class_UBA796",
  "class_XYA12-FULL-58-9",
  "class_Alphaproteobacteria",
  "class_Gammaproteobacteria",
  "class_SAR324"
)

# Model stats
fnames_model_stats <- c(
  # Archaea
  "../2_Extract_model_outputs/3_Output/Non_log/class_Poseidoniia_Richness_5000_v2.RData",
  "../2_Extract_model_outputs/3_Output/Non_log/class_Nitrososphaeria_Richness_5000_v2.RData",
  # Bacteria
  "../2_Extract_model_outputs/3_Output/Non_log/class_Acidimicrobiia_Richness_5000_v2.RData",
  "../2_Extract_model_outputs/3_Output/Log/class_Bacteroidia_Richness_5000_v2.RData",
  "../2_Extract_model_outputs/3_Output/Non_log/class_Chlamydiia_Richness_5000_v2.RData",
  "../2_Extract_model_outputs/3_Output/Non_log/class_Cyanobacteriia_Richness_5000_v2.RData",
  "../2_Extract_model_outputs/3_Output/Non_log/class_UBA1144_Richness_5000_v2.RData",
  "../2_Extract_model_outputs/3_Output/Non_log/class_Fibrobacteria_Richness_5000_v2.RData",
  "../2_Extract_model_outputs/3_Output/Non_log/class_Marinisomatia_Richness_5000_v2.RData",
  "../2_Extract_model_outputs/3_Output/Non_log/class_UBA4151_Richness_5000_v2.RData",
  "../2_Extract_model_outputs/3_Output/Non_log/class_UBA796_Richness_5000_v2.RData",
  "../2_Extract_model_outputs/3_Output/Non_log/class_XYA12-FULL-58-9_Richness_5000_v2.RData",
  "../2_Extract_model_outputs/3_Output/Non_log/class_Alphaproteobacteria_Richness_5000_v2.RData",
  "../2_Extract_model_outputs/3_Output/Non_log/class_Gammaproteobacteria_Richness_5000_v2.RData",
  "../2_Extract_model_outputs/3_Output/Non_log/class_SAR324_Richness_5000_v2.RData"
)

# Reorder fnames_model_stats to match clade order
fnames_model_stats_ordered <- sapply(
  clade,
  function(cl) fnames_model_stats[grep(cl, fnames_model_stats)]
)

# Titles for each clade
titles <- c(
  "Archaea - Thermoplasmatota - Poseidoniia", 
  "Archaea - Thermoproteota - Nitrososphaeria", 
  "Bacteria - Actinomycetota - Acidimicrobiia", 
  "Bacteria - Bacteroidota - Bacteroidia", 
  "Bacteria - Chlamydiota - Chlamydiia", 
  "Bacteria - Cyanobacteriota - Cyanobacteriia", 
  "Bacteria - Desulfobacterota D - UBA1144", 
  "Bacteria - Fibrobacterota - Fibrobacteria", 
  "Bacteria - Marinisomatota - Marinisomatia", 
  "Bacteria - Myxococcota - UBA4151", 
  "Bacteria - Myxococcota - UBA796", 
  "Bacteria - Myxococcota - XYA12-FULL-58-9", 
  "Bacteria - Pseudomonadota - Alphaproteobacteria",
  "Bacteria - Pseudomonadota - Gammaproteobacteria",
  "Bacteria - SAR324 - SAR324"
)

### ====================================================================================================================
### Generate Individual Clade Maps
### ====================================================================================================================

# Loop through each file and create the plots
plot_list <- list()
for (i in seq_along(clade)){

  # Print progress
  print(paste0("Processing clade: ", clade[i], " (", i, "/", length(clade), ")"))
  
  # Convert to spatraster object
  df_rich$clade <- as.character(df_rich$clade)  # Ensure character type
  target_clade <- trimws(clade[i])
  df_div_sub <- dplyr::filter(df_rich, trimws(clade) == target_clade)
  dat <- terra::rast(df_div_sub[, c("x", "y", "mean_richness")])
  # Set reference system
  terra::crs(dat) <- 'epsg:4326'

  # Uncertainty stippling  
  mean_min_threshold <- 1     # For richness
  cv_abs_threshold <- 20       # Absolute CV threshold (%)

# Compute coefficient of variation
  df_div <- df_rich %>%
    dplyr::mutate(
      mean_richness = ifelse(mean_richness < mean_min_threshold, NA, mean_richness),
      CV = 100 * sd_richness / mean_richness,
      CV = ifelse(is.na(mean_richness), NA, CV)
    )

  # Compute the 70th percentile of CV for each clade
  cv_percentile <- df_div %>%
  dplyr::group_by(clade) %>%
  dplyr::summarise(cv_70th = quantile(CV, probs = 0.7, na.rm = TRUE))
  df_div <- dplyr::left_join(df_div, cv_percentile, by = "clade")

  df_div <- df_div %>%
    dplyr::mutate(
      final_cv_threshold = pmax(cv_abs_threshold, cv_70th, na.rm = TRUE),
      stipple = as.integer(!is.na(CV) & CV >= final_cv_threshold & !is.na(mean_richness))
    )

  # Only keep stippling points and subset for clade
  uncertainty <- df_div %>%
    dplyr::filter(clade == clade[i], stipple == 1)

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

  if( length(contour_intervals) < 7 ){
    # If richness low in general
    contour_intervals <- unique(round(quantiles, 2))
  }
  
  # Ensure the intervals cover the entire range
  contour_intervals <- sort(unique(c(min_val, contour_intervals, max_val)))
  
  # Create annotations for plot
  model_stats <- get(load(fnames_model_stats[i]))
  predictive_power <- round(unique(model_stats[[1]]$Predictive_Power), 2)
  # number_samples <- nrow(df_sub_obs)
  plot_title <- titles[i]
    wrap_title <- function(title, width = 30) {
    stringr::str_wrap(title, width = width)
    }
    plot_title_wrapped <- wrap_title(plot_title)

  
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
    geom_sf(fill = "grey", size = 0.1) +
    tidyterra::geom_spatraster_contour_filled(
      data = dat, 
      show.legend = TRUE, 
      breaks = contour_intervals
    ) +
    geom_sf(data = stipple_points_sf, fill = "white", alpha = 1, size = 0.5, shape = 21, stroke = 0.1) +
    # geom_sf(data = df_sub_sf, aes(size = measurementvalue), shape = 21, fill = "black", alpha = 0.6) +
    # scale_size_continuous(range = c(0.1, 1)) +
    annotate(
      "text", 
      x = text_coords_transformed[1, "X"], 
      y = text_coords_transformed[1, "Y"],
      label = bquote("R"^2 * ": " * .(predictive_power)), 
      size = 6 / .pt, 
      color = "red"
    ) +
    labs(
      title = plot_title_wrapped,
      size = "Within sample richness",
      fill = "Modeled richness"
    ) +
    theme_minimal(base_size = 8) +
    theme(
      text = element_text(size = 7),
      legend.position = "right",
      legend.direction = "vertical",
      legend.key.size = unit(0.12, 'cm'),
      legend.key.height = unit(0.12, 'cm'),
      legend.key.width = unit(0.18, 'cm'),
      legend.title = element_text(size = 7, hjust = 0.5),
      legend.text = element_text(size = 6),
      legend.box.margin = margin(0, 0, 0, 0, "mm"),
      legend.box.spacing = unit(0, "mm"),
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
      size = guide_legend(position = "bottom", keywidth = unit(0.18, "cm"), keyheight = unit(0.12, "cm")),
      fill = guide_legend(position = "right", keywidth = unit(0.18, "cm"), keyheight = unit(0.12, "cm"))
    )

  # Calculate the 75th percentile threshold
  contour_75 <- quantile(values(dat), probs = 0.75, na.rm = TRUE)

  # Add the contour line to your plot
  plot <- plot +
    tidyterra::geom_spatraster_contour(
      data = dat,
      breaks = contour_75,
      color = "white",
      size = 0.3,
      linetype = "solid"
    )

  # Remove the last value (upper bound) for labeling
  fill_labels <- as.character(contour_intervals[-length(contour_intervals)])


  # Generate a color vector with the same length as fill_labels
  colors <- viridis(length(fill_labels))

  # Then use in your scale_fill_manual
  plot <- plot +
    scale_fill_manual(
      values = colors,
      labels = fill_labels,
      name = "Modeled richness"
    )

  # Save in plot_list
  plot_list[[i]] <- plot

} # close loop across clades

### ====================================================================================================================
### Create Multiplot and Save Outputs
### ====================================================================================================================

# Arrange multiplot: 3 columns, adjust rows as needed
ncol <- 3
nrow <- ceiling(length(plot_list) / ncol)
multiplot <- cowplot::plot_grid(plotlist = plot_list, ncol = ncol, nrow = nrow, align = "hv")


# Print multiplot
print(multiplot)

# Save multiplot (A4 is 8.27 x 11.69 inches, so use less)
ggsave(
  filename = paste0(wd_out, "GlobalMap_annualRichness_acrossClades_Multiplot_v2.svg"),
  plot = multiplot,
  width = 8, height = 10,  # Less than A4
  dpi = 300
)

# Save multiplot (A4 is 8.27 x 11.69 inches, so use less)
ggsave(
  filename = paste0(wd_out, "GlobalMap_annualRichness_acrossClades_Multiplot_v2.png"),
  plot = multiplot,
  width = 8, height = 10,  # Less than A4
  dpi = 300
)

## ====================================================================================================
## END OF SCRIPT
## ====================================================================================================
