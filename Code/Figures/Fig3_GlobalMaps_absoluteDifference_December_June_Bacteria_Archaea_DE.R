## ====================================================================================================
## This script generates global maps of seasonal differences in bacterial and archaeal richness (December – June),
## highlighting areas of increase or decrease and visualizing regions of high model uncertainty.
##
## Author:       Dominic Eriksson
## Date:         8th of October, 2025
## Affiliation:  Environmental Physics Group UP, ETH Zürich, Switzerland
## Contact:      deriksson@ethz.ch
##
## Input files:  
##   - Seasonal richness difference CSV: Code/3_Compute_ensembles_and_uncertainties/3_Output/Global_map_difference_December_June_ensemble__mean_sd_richness.csv
##   - Ocean regions shapefile: GOaS_v1_20211214/goas_v01.shp
##
## Output files: 
##   - Global maps of Δ richness for bacterial and archaeal clades
##
## Strategy:
##   This script visualizes seasonal changes in microbial richness while accounting for model uncertainty:
##   1. Load seasonal Δ richness data for bacterial and archaeal clades.
##   2. Transform spatial data into Pacific-centered projection and read ocean basin shapefiles.
##   3. Convert Δ richness values into a raster object for plotting.
##   4. Define symmetric color bins around zero to highlight increases (red) and decreases (blue):
##       - Include a narrow white bin around zero for visual balance.
##   5. Identify grid cells with model agreement below 75% and generate stippling masks.
##   6. Convert stippling points into Pacific-centered CRS for plotting.
##   7. Create contour intervals based on Δ richness values to represent spatial gradients.
##   8. Combine raster, contour, and stippling layers into ggplot figures for each clade.
##   9. Apply tight, publication-ready styling with legends and color scales.
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
### Richness
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

# Working directories
wd_in <- "Code/3_Compute_ensembles_and_uncertainties/3_Output/Global_map_difference_December_June_ensemble__mean_sd_richness.csv"
wd_out <- "Code/Figures/Fig3/"

# Create output directory if it doesn't exist
if (!dir.exists(wd_out)) {
  dir.create(wd_out, recursive = TRUE)
}


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

# Load the data
df <- data.table::fread(wd_in)

# Convert mean_pct_increase to numeric
df$mean_diff <- as.numeric(df$mean_diff)
df$sd_diff <- as.numeric(df$sd_diff)
  

### Bacteria -------------------------------------------------

# Filter for d__Bacteria or d__Archaea in clade column
df_sub <- df %>%
    dplyr::filter(clade == "domain_Bacteria") %>%
    dplyr::select(clade, x, y, mean_diff, 
                  sd_diff)

# Convert to spatraster object
dat <- terra::rast(df_sub[, c("x", "y", "mean_diff")])

# Set reference system
terra::crs(dat) <- 'epsg:4326'
  
# Automatic binning with 10 bins and white centered at zero
data_values <- terra::values(dat, na.rm = TRUE)
data_min <- min(data_values, na.rm = TRUE)
data_max <- max(data_values, na.rm = TRUE)

# Calculate symmetric range around zero for balanced color scheme
max_abs <- max(abs(data_min), abs(data_max))

# Define white bin around zero (narrow range)
white_range <- max_abs * 0.05  # 5% of max absolute value for white bin
zero_lower <- -white_range
zero_upper <- white_range

# Create 4 bins below zero, 1 white bin at zero, 4 bins above zero (total 9 intervals = 10 bins)
neg_bins <- seq(data_min, zero_lower, length.out = 5)  # 4 negative bins
pos_bins <- seq(zero_upper, data_max, length.out = 5)  # 4 positive bins

# Combine all intervals
contour_intervals <- c(neg_bins, zero_upper, pos_bins[-1])

# Remove any duplicates
contour_intervals <- unique(sort(contour_intervals))

print(contour_intervals)

# Calculate bins and colors
n_bins <- length(contour_intervals) - 1
zero_bin <- findInterval(0, contour_intervals, all.inside = TRUE)

# Count bins below and above zero
n_below <- zero_bin - 1  # bins below zero
n_above <- n_bins - zero_bin  # bins above zero

# Create colors: blue → white → red
colors_neg <- colorRampPalette(c("blue", "white"))(n_below + 1)
colors_pos <- colorRampPalette(c("white", "red"))(n_above + 1)

# Combine colors, avoiding duplicate white
colors <- c(colors_neg, colors_pos[-1])

# Get model uncertainties
uncertainty <- df %>%
    dplyr::filter(clade == "domain_Bacteria") %>%
    dplyr::select(clade, x, y, agreement_diff)

# Remove NAs
uncertainty <- na.omit(uncertainty)

# Stipple areas with model agreement < 75%
uncertainty <- uncertainty[uncertainty$agreement_diff < 75, ]

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

# Generate world map
plot <- ggplot(sf.ocean_pacific) +

  # Add basemap
  geom_sf(fill = "grey", size = 0.1) +

  # Contour plot with filled contours
  tidyterra::geom_spatraster_contour_filled(
    data = dat,
    show.legend = TRUE,
    breaks = contour_intervals
  ) +
  scale_fill_manual(
    values = colors,
    name = "Richness Δ"
  ) +

  # Add stippling upper 70th percentile bound
  geom_sf(data = stipple_points_sf, fill = "black", alpha = 0.7, size = 0.5, shape = 21, stroke = 0.1) +

  coord_sf(expand = FALSE) +
  theme_minimal(base_size = 8) +
  theme(
    text = element_text(size = 8),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.key.size = unit(0.2, 'cm'),
    legend.key.height = unit(0.2, 'cm'),
    legend.key.width = unit(0.6, 'cm'),
    legend.title = element_text(size = 8, hjust = 0.5),
    legend.text = element_text(size = 8),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
    plot.margin = margin(2, 2, 2, 2, "mm")
  ) +
  xlab("") +
  ylab("") +
guides(
  fill = guide_legend(
    title.position = "top",
    title.hjust = 0.5,
    nrow = 1,
    byrow = TRUE,
    keywidth = unit(0.18, "cm"),     # very small keys
    keyheight = unit(0.15, "cm"),
    label.position = "bottom",
    label.hjust = 0.5,
    label.theme = element_text(size = 4, angle = 90, vjust = 0.5, hjust = 1)  # tiny font
  )
) +
theme(
  legend.key.size = unit(0.08, 'cm'),   # minimal key size
  legend.title = element_text(size = 5, hjust = 0.5), 
  plot.margin = margin(1, 2, 2, 1, "mm")  # very tight margins
)

# Print plot
print(plot)

# Save as high resolution image
    ggsave(
      filename = paste0(wd_out, "GlobalMap_percentageIncrease_Bacteria_richness.svg"),
      plot = plot,
      width = 3, height = 3, 
      dpi = 300
    )


### Archaea --------------------------------

# Filter Archaea in clade column
df_sub <- df %>%
    dplyr::filter(clade == "domain_Archaea") %>%
    dplyr::select(clade, x, y, mean_diff)
  
# Convert to spatraster object
dat <- terra::rast(df_sub[, c("x", "y", "mean_diff")])

# Set reference system
terra::crs(dat) <- 'epsg:4326'
  
# Automatic binning with 10 bins and white centered at zero
data_values <- terra::values(dat, na.rm = TRUE)
data_min <- min(data_values, na.rm = TRUE)
data_max <- max(data_values, na.rm = TRUE)

# Calculate symmetric range around zero for balanced color scheme
max_abs <- max(abs(data_min), abs(data_max))

# Define white bin around zero (narrow range)
white_range <- max_abs * 0.05  # 5% of max absolute value for white bin
zero_lower <- -white_range
zero_upper <- white_range

# Create 4 bins below zero, 1 white bin at zero, 4 bins above zero (total 9 intervals = 10 bins)
neg_bins <- seq(data_min, zero_lower, length.out = 5)  # 4 negative bins
pos_bins <- seq(zero_upper, data_max, length.out = 5)  # 4 positive bins

# Combine all intervals
contour_intervals <- c(neg_bins, zero_upper, pos_bins[-1])

# Remove any duplicates
contour_intervals <- unique(sort(contour_intervals))

print(contour_intervals)

# Calculate bins and colors
n_bins <- length(contour_intervals) - 1
zero_bin <- findInterval(0, contour_intervals, all.inside = TRUE)

# Count bins below and above zero
n_below <- zero_bin - 1  # bins below zero
n_above <- n_bins - zero_bin  # bins above zero

# Create colors: blue → white → red
colors_neg <- colorRampPalette(c("blue", "white"))(n_below + 1)
colors_pos <- colorRampPalette(c("white", "red"))(n_above + 1)

# Combine colors, avoiding duplicate white
colors <- c(colors_neg, colors_pos[-1])

# Get model uncertainties
uncertainty <- df %>%
    dplyr::filter(clade == "domain_Archaea") %>%
    dplyr::select(clade, x, y, agreement_diff)

# Remove NAs
uncertainty <- na.omit(uncertainty)

# Only keep rows with agreement_pct < 75
uncertainty <- uncertainty[uncertainty$agreement_diff < 75, ]
  
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

# Generate world map
plot <- ggplot(sf.ocean_pacific) +

  # Add basemap
  geom_sf(fill = "grey", size = 0.1) +

  # Contour plot with filled contours
  tidyterra::geom_spatraster_contour_filled(
    data = dat,
    show.legend = TRUE,
    breaks = contour_intervals
  ) +
  scale_fill_manual(
    values = colors,
    name = expression("Richness" ~ Delta)
  ) +

  # Add stippling upper 70th percentile bound
  geom_sf(data = stipple_points_sf, fill = "black", alpha = 0.7, size = 0.5, shape = 21, stroke = 0.1) +

  coord_sf(expand = FALSE) +
  theme_minimal(base_size = 8) +
  theme(
    text = element_text(size = 8),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.key.size = unit(0.2, 'cm'),
    legend.key.height = unit(0.2, 'cm'),
    legend.key.width = unit(0.6, 'cm'),
    legend.title = element_text(size = 8, hjust = 0.5),
    legend.text = element_text(size = 8),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
    plot.margin = margin(2, 2, 2, 2, "mm")
  ) +
  xlab("") +
  ylab("") +
guides(
  fill = guide_legend(
    title.position = "top",
    title.hjust = 0.5,
    nrow = 1,
    byrow = TRUE,
    keywidth = unit(0.18, "cm"),     # very small keys
    keyheight = unit(0.15, "cm"),
    label.position = "bottom",
    label.hjust = 0.5,
    label.theme = element_text(size = 4, angle = 90, vjust = 0.5, hjust = 1)  # tiny font
  )
) +
theme(
  legend.key.size = unit(0.08, 'cm'),   # minimal key size
  legend.title = element_text(size = 5, hjust = 0.5), 
  plot.margin = margin(1, 2, 2, 1, "mm")  # very tight margins
)

# Print plot
print(plot)

# Save as high resolution image
    ggsave(
      filename = paste0(wd_out, "GlobalMap_percentageIncrease_Archaea_v2.svg"),
      plot = plot,
      width = 3, height = 3, 
      dpi = 300
    )

## ====================================================================================================
## END OF SCRIPT
## ====================================================================================================
