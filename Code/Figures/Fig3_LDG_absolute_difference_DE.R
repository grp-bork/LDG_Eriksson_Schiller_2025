##====================================================================================================
## Figure 3: Global Maps of Absolute Richness Differences Between December and June
##
## This R script generates global maps showing absolute prokaryotic richness
## differences between December and June for both Bacteria and Archaea. The analysis uses
## absolute values to highlight winter impact as suggested by reviewers, employing a white-to-red
## color scheme instead of the traditional blue-white-red diverging scale. Maps are generated
## with Pacific-centered projection, uncertainty stippling, and publication-quality formatting.
##
## Author:      Dominic Eriksson
##              Environmental Physics Group, UP
##              ETH Zürich, Switzerland
## Contact:     deriksson@ethz.ch
## Date:        February 27th, 2026
## Affiliation: ETH Zürich, Environmental Physics Group, UP
##
## Input files:
## - ./Code/4_Compute_ensembles_uncertainties/5_Output/Global_map_difference_December_June_ensemble__mean_sd_richness_v2.csv
## - ./Data/Ocean_regions_shapefile/GOaS_v1_20211214/goas_v01.shp (Ocean basin geometries)
## - ./Code/4_Compute_ensembles_uncertainties/5_Output/Global_map_difference_December_June_ensemble_mean_sd_shannon_v2.csv
## - ./Code/4_Compute_ensembles_uncertainties/5_Output/Global_map_difference_December_June_ensemble_mean_sd_chao1_v2.csv
##
## Output files:
## - GlobalMap_absolute_richnessIncrease_Bacteria_richness.svg/.png (Bacteria absolute richness difference map)
## - GlobalMap_absolute_richnessIncrease_Archaea.svg/.png (Archaea absolute richness difference map)
##
## Strategy:
## 1. Load and prepare environmental data with Pacific-centered world map projection
## 2. Process richness difference data for Bacteria and Archaea domains
## 3. Calculate absolute values of richness differences to highlight winter impact
## 4. Create white-to-red color scheme for absolute value visualization
## 5. Generate uncertainty stippling for model agreement < 75%
## 6. Create publication-quality global maps with proper legends and formatting
## 7. Calculate hemisphere-specific statistics (mean, SD, median, IQR)
## 8. Perform correlation analysis between diversity indices (richness, Shannon, Chao1)
##
## Required R packages:
## - ggplot2 (3.3.6): Advanced data visualization and publication-quality plotting
## - dplyr (1.0.9): Data manipulation, filtering, and transformation operations
## - sf (1.0.7): Spatial data handling and coordinate reference system management
## - tidyterra (0.3.1): Integration between terra and ggplot2 for raster visualization
## - cowplot (1.1.1): Publication-ready plot themes and multi-panel arrangements
## - terra (1.6.17): Modern spatial data analysis and raster processing
## - hexbin (1.28.2): Hexagonal binning for uncertainty stippling visualization
## - data.table (1.14.2): Fast and efficient data manipulation and file reading
## - rnaturalearth (0.1.0): World map data for cartographic visualization
##====================================================================================================

### ==============================================================================================
### Richness Absolute Difference Analysis (Reviewer Suggested)
### ==============================================================================================

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

# Working directories
wd_in <- "./Code/4_Compute_ensembles_uncertainties/5_Output/Global_map_difference_December_June_ensemble__mean_sd_richness_v2.csv"
wd_out <- "./Output/"
# Create output directory if it doesn't exist
if (!dir.exists(wd_out)) {
  dir.create(wd_out, recursive = TRUE)
}

## Required data & data formatting
# Read shapefile with geometries of ocean basins - Note: Specific shapefile needed for Pacific-centered map
sf.ocean <- sf::st_read("./Data/Ocean_regions_shapefile/GOaS_v1_20211214/goas_v01.shp")
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
df$mean_december_june_richnessDiff <- as.numeric(df$mean_december_june_richnessDiff)
df$sd_december_june_richnessDiff <- as.numeric(df$sd_december_june_richnessDiff)
  

### Bacteria - Absolute Richness Difference (Reviewer Suggested Approach) -------

# Filter for d__Bacteria in clade column
df_sub <- df %>%
    dplyr::filter(clade == "domain_Bacteria") %>%
    dplyr::select(clade, x, y, mean_december_june_richnessDiff, 
                  sd_december_june_richnessDiff)

# Convert to spatraster object and take absolute values directly
dat <- terra::rast(df_sub[, c("x", "y", "mean_december_june_richnessDiff")])
# Set reference system
terra::crs(dat) <- 'epsg:4326'

# Take absolute values to highlight winter impact (reviewer suggestion)
dat_abs <- abs(dat)

# Get data values for binning (absolute values only)
data_values <- terra::values(dat_abs, na.rm = TRUE)
data_min <- min(data_values, na.rm = TRUE)  # This will be 0 or close to 0
data_max <- max(data_values, na.rm = TRUE)

# Create bins from 0 to max (only positive values)
contour_intervals <- seq(data_min, data_max, length.out = 10)

print(contour_intervals)

# Create colors: white → red (no blue)
n_bins <- length(contour_intervals) - 1
colors <- colorRampPalette(c("white", "red"))(n_bins)

# Get model uncertainties
uncertainty <- df %>%
    dplyr::filter(clade == "domain_Bacteria") %>%
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

  # Contour plot with filled contours using absolute values
  tidyterra::geom_spatraster_contour_filled(
    data = dat_abs,  # Use absolute values
    show.legend = TRUE,
    breaks = contour_intervals
  ) +
  scale_fill_manual(
    values = colors,
    name = "|Richness Δ|"  # Updated legend title to show absolute values
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
      filename = paste0(wd_out, "GlobalMap_absolute_richnessIncrease_Bacteria_richness.svg"),
      plot = plot,
      width = 3, height = 3, 
      dpi = 300
    )
#
    ggsave(
      filename = paste0(wd_out, "GlobalMap_absolute_richnessIncrease_Bacteria_richness.png"),
      plot = plot,
      width = 3, height = 3, 
      dpi = 300
    )



# Convert raster to data frame with coordinates
df_map <- as.data.frame(dat, xy = TRUE, na.rm = TRUE)
names(df_map)[3] <- "richness_diff"

# Northern Hemisphere (y > 0)
north <- df_map[df_map$y > 0 & !is.na(df_map$richness_diff), ]
mean_north <- mean(north$richness_diff)
sd_north <- sd(north$richness_diff)
median_north <- median(north$richness_diff)
iqr_north <- quantile(north$richness_diff, probs = c(0.25, 0.75))

# Southern Hemisphere (y < 0)
south <- df_map[df_map$y < 0 & !is.na(df_map$richness_diff), ]
mean_south <- mean(south$richness_diff)
sd_south <- sd(south$richness_diff)
median_south <- median(south$richness_diff)
iqr_south <- quantile(south$richness_diff, probs = c(0.25, 0.75))

cat("Northern Hemisphere:\n",
    "Mean =", mean_north, "\n",
    "SD =", sd_north, "\n",
    "Median =", median_north, "\n",
    "IQR =", iqr_north[1], "to", iqr_north[2], "\n\n")

cat("Southern Hemisphere:\n",
    "Mean =", mean_south, "\n",
    "SD =", sd_south, "\n",
    "Median =", median_south, "\n",
    "IQR =", iqr_south[1], "to", iqr_south[2], "\n")

### Archaea - Absolute Richness Difference (Reviewer Suggested Approach) ------

# Filter Archaea in clade column
df_sub <- df %>%
    dplyr::filter(clade == "domain_Archaea") %>%
    dplyr::select(clade, x, y, mean_december_june_richnessDiff)
  
# Convert to spatraster object and take absolute values directly
dat <- terra::rast(df_sub[, c("x", "y", "mean_december_june_richnessDiff")])
# Set reference system
terra::crs(dat) <- 'epsg:4326'

# Take absolute values to highlight winter impact (reviewer suggestion)
dat_abs <- abs(dat)

# Get data values for binning (absolute values only)
data_values <- terra::values(dat_abs, na.rm = TRUE)
data_min <- min(data_values, na.rm = TRUE)  # This will be 0 or close to 0
data_max <- max(data_values, na.rm = TRUE)

# Create bins from 0 to max (only positive values)
contour_intervals <- seq(data_min, data_max, length.out = 10)

print(contour_intervals)

# Create colors: white → red (no blue)
n_bins <- length(contour_intervals) - 1
colors <- colorRampPalette(c("white", "red"))(n_bins)

# Get model uncertainties for Archaea
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
      filename = paste0(wd_out, "GlobalMap_richnessIncrease_Archaea_v3.svg"),
      plot = plot,
      width = 3, height = 3, 
      dpi = 300
    )

# As suggested by reviewer, use absolute values to highlight winter impact better

# Take absolute values of the data for Archaea
dat_abs <- abs(dat)

# Get data values for binning
data_values <- terra::values(dat_abs, na.rm = TRUE)
data_min <- min(data_values, na.rm = TRUE)  # This will be 0 or close to 0
data_max <- max(data_values, na.rm = TRUE)

# Create bins from 0 to max (only positive values)
contour_intervals <- seq(data_min, data_max, length.out = 10)

print(contour_intervals)

# Create colors: white → red (no blue)
n_bins <- length(contour_intervals) - 1
colors <- colorRampPalette(c("white", "red"))(n_bins)

# Generate world map
plot <- ggplot(sf.ocean_pacific) +

  # Add basemap
  geom_sf(fill = "grey", size = 0.1) +

  # Contour plot with filled contours using absolute values
  tidyterra::geom_spatraster_contour_filled(
    data = dat_abs,  # Use absolute values
    show.legend = TRUE,
    breaks = contour_intervals
  ) +
  scale_fill_manual(
    values = colors,
    name = expression("|Richness" ~ Delta ~ "|")  # Updated legend title for absolute values
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
  filename = paste0(wd_out, "GlobalMap_absolute_richnessIncrease_Archaea.svg"),
  plot = plot,
  width = 3, height = 3, 
  dpi = 300
)

ggsave(
  filename = paste0(wd_out, "GlobalMap_absolute_richnessIncrease_Archaea.png"),
  plot = plot,
  width = 3, height = 3, 
  dpi = 300
)

# Convert raster to data frame with coordinates
df_map <- as.data.frame(dat, xy = TRUE, na.rm = TRUE)
names(df_map)[3] <- "richness_diff"

# Northern Hemisphere (y > 0)
north <- df_map[df_map$y > 0 & !is.na(df_map$richness_diff), ]
mean_north <- mean(north$richness_diff)
sd_north <- sd(north$richness_diff)
median_north <- median(north$richness_diff)
iqr_north <- quantile(north$richness_diff, probs = c(0.25, 0.75))

# Southern Hemisphere (y < 0)
south <- df_map[df_map$y < 0 & !is.na(df_map$richness_diff), ]
mean_south <- mean(south$richness_diff)
sd_south <- sd(south$richness_diff)
median_south <- median(south$richness_diff)
iqr_south <- quantile(south$richness_diff, probs = c(0.25, 0.75))

cat("Northern Hemisphere:\n",
    "Mean =", mean_north, "\n",
    "SD =", sd_north, "\n",
    "Median =", median_north, "\n",
    "IQR =", iqr_north[1], "to", iqr_north[2], "\n\n")

cat("Southern Hemisphere:\n",
    "Mean =", mean_south, "\n",
    "SD =", sd_south, "\n",
    "Median =", median_south, "\n",
    "IQR =", iqr_south[1], "to", iqr_south[2], "\n")
    


### ============================================================
### Correlation across diversity indices
### ============================================================

# Load files
richness <- fread("./Code/4_Compute_ensembles_uncertainties/5_Output/Global_map_difference_December_June_ensemble__mean_sd_richness_v2.csv")
shannon <- fread("./Code/4_Compute_ensembles_uncertainties/5_Output/Global_map_difference_December_June_ensemble_mean_sd_shannon_v2.csv")
chao1 <- fread("./Code/4_Compute_ensembles_uncertainties/5_Output/Global_map_difference_December_June_ensemble_mean_sd_chao1_v2.csv")

# Quick Shannon map for Bacteria
shannon_sub <- shannon %>%
    dplyr::filter(clade == "domain_Bacteria") %>%
    dplyr::select(x, y, mean_december_june_shannonDiff) %>%
    na.omit()

# Convert to spatraster object
shannon_rast <- terra::rast(shannon_sub)
terra::crs(shannon_rast) <- 'epsg:4326'

# Quick plot
ggplot(sf.ocean_pacific) +
  geom_sf(fill = "grey", size = 0.1) +
  tidyterra::geom_spatraster(
    data = shannon_rast,
    show.legend = TRUE
  ) +
  scale_fill_gradient2(
    low = "blue", 
    mid = "white", 
    high = "red",
    midpoint = 0,
    name = "Shannon Δ",
    na.value = "transparent"
  ) +
  coord_sf(expand = FALSE) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  ) +
  ggtitle("Shannon Diversity Difference - Bacteria")

# Quick Chao1 map for Bacteria
chao1_sub <- chao1 %>%
    dplyr::filter(clade == "domain_Bacteria") %>%
    dplyr::select(x, y, mean_december_june_chao1Diff) %>%
    na.omit()

# Convert to spatraster object
chao1_rast <- terra::rast(chao1_sub)
terra::crs(chao1_rast) <- 'epsg:4326'

# Quick plot
ggplot(sf.ocean_pacific) +
  geom_sf(fill = "grey", size = 0.1) +
  tidyterra::geom_spatraster(
    data = chao1_rast,
    show.legend = TRUE
  ) +
  scale_fill_gradient2(
    low = "blue", 
    mid = "white", 
    high = "red",
    midpoint = 0,
    name = "Chao1 Δ",
    na.value = "transparent"
  ) +
  coord_sf(expand = FALSE) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  ) +
  ggtitle("Chao1 Diversity Difference - Bacteria")


# Correlation analysis between diversity indices

# Merge datasets for Bacteria
bacteria_richness <- richness %>%
    dplyr::filter(clade == "domain_Bacteria") %>%
    dplyr::select(x, y, mean_december_june_richnessDiff) %>%
    na.omit()

bacteria_shannon <- shannon %>%
    dplyr::filter(clade == "domain_Bacteria") %>%
    dplyr::select(x, y, mean_december_june_shannonDiff) %>%
    na.omit()

bacteria_chao1 <- chao1 %>%
    dplyr::filter(clade == "domain_Bacteria") %>%
    dplyr::select(x, y, mean_december_june_chao1Diff) %>%
    na.omit()

# Merge all three datasets
combined_data <- bacteria_richness %>%
    dplyr::inner_join(bacteria_shannon, by = c("x", "y")) %>%
    dplyr::inner_join(bacteria_chao1, by = c("x", "y"))

# Calculate correlations and significance tests
# Richness vs Shannon
cor_richness_shannon <- cor.test(combined_data$mean_december_june_richnessDiff, 
                                combined_data$mean_december_june_shannonDiff, 
                                method = "pearson")

# Richness vs Chao1
cor_richness_chao1 <- cor.test(combined_data$mean_december_june_richnessDiff, 
                              combined_data$mean_december_june_chao1Diff, 
                              method = "pearson")

# Shannon vs Chao1
cor_shannon_chao1 <- cor.test(combined_data$mean_december_june_shannonDiff, 
                             combined_data$mean_december_june_chao1Diff, 
                             method = "pearson")

# Print results
cat("Correlation Results for Bacteria:\n\n")

cat("Richness vs Shannon:\n")
cat("Pearson r =", round(cor_richness_shannon$estimate, 3), "\n")
cat("p-value =", format(cor_richness_shannon$p.value, scientific = TRUE, digits = 3), "\n")
cat("95% CI: [", round(cor_richness_shannon$conf.int[1], 3), ", ", 
    round(cor_richness_shannon$conf.int[2], 3), "]\n")
cat("n =", nrow(combined_data), "\n\n")

cat("Richness vs Chao1:\n")
cat("Pearson r =", round(cor_richness_chao1$estimate, 3), "\n")
cat("p-value =", format(cor_richness_chao1$p.value, scientific = TRUE, digits = 3), "\n")
cat("95% CI: [", round(cor_richness_chao1$conf.int[1], 3), ", ", 
    round(cor_richness_chao1$conf.int[2], 3), "]\n")
cat("n =", nrow(combined_data), "\n\n")

cat("Shannon vs Chao1:\n")
cat("Pearson r =", round(cor_shannon_chao1$estimate, 3), "\n")
cat("p-value =", format(cor_shannon_chao1$p.value, scientific = TRUE, digits = 3), "\n")
cat("95% CI: [", round(cor_shannon_chao1$conf.int[1], 3), ", ", 
    round(cor_shannon_chao1$conf.int[2], 3), "]\n")
cat("n =", nrow(combined_data), "\n")



## Archaea correlation analysis -------------
# Merge datasets for Archaea
archaea_richness <- richness %>%
    dplyr::filter(clade == "domain_Archaea") %>%
    dplyr::select(x, y, mean_december_june_richnessDiff) %>%
    na.omit()

archaea_shannon <- shannon %>%
    dplyr::filter(clade == "domain_Archaea") %>%
    dplyr::select(x, y, mean_december_june_shannonDiff) %>%
    na.omit()

archaea_chao1 <- chao1 %>%
    dplyr::filter(clade == "domain_Archaea") %>%
    dplyr::select(x, y, mean_december_june_chao1Diff) %>%
    na.omit()

# Merge all three datasets
combined_data_archaea <- archaea_richness %>%
    dplyr::inner_join(archaea_shannon, by = c("x", "y")) %>%
    dplyr::inner_join(archaea_chao1, by = c("x", "y"))

# Calculate correlations and significance tests
# Richness vs Shannon
cor_richness_shannon_archaea <- cor.test(combined_data_archaea$mean_december_june_richnessDiff, 
                                        combined_data_archaea$mean_december_june_shannonDiff, 
                                        method = "pearson")

# Richness vs Chao1
cor_richness_chao1_archaea <- cor.test(combined_data_archaea$mean_december_june_richnessDiff, 
                                      combined_data_archaea$mean_december_june_chao1Diff, 
                                      method = "pearson")

# Shannon vs Chao1
cor_shannon_chao1_archaea <- cor.test(combined_data_archaea$mean_december_june_shannonDiff, 
                                     combined_data_archaea$mean_december_june_chao1Diff, 
                                     method = "pearson")

# Print results
cat("Correlation Results for Archaea:\n\n")

cat("Richness vs Shannon:\n")
cat("Pearson r =", round(cor_richness_shannon_archaea$estimate, 3), "\n")
cat("p-value =", format(cor_richness_shannon_archaea$p.value, scientific = TRUE, digits = 3), "\n")
cat("95% CI: [", round(cor_richness_shannon_archaea$conf.int[1], 3), ", ", 
    round(cor_richness_shannon_archaea$conf.int[2], 3), "]\n")
cat("n =", nrow(combined_data_archaea), "\n\n")

cat("Richness vs Chao1:\n")
cat("Pearson r =", round(cor_richness_chao1_archaea$estimate, 3), "\n")
cat("p-value =", format(cor_richness_chao1_archaea$p.value, scientific = TRUE, digits = 3), "\n")
cat("95% CI: [", round(cor_richness_chao1_archaea$conf.int[1], 3), ", ", 
    round(cor_richness_chao1_archaea$conf.int[2], 3), "]\n")
cat("n =", nrow(combined_data_archaea), "\n\n")

cat("Shannon vs Chao1:\n")
cat("Pearson r =", round(cor_shannon_chao1_archaea$estimate, 3), "\n")
cat("p-value =", format(cor_shannon_chao1_archaea$p.value, scientific = TRUE, digits = 3), "\n")
cat("95% CI: [", round(cor_shannon_chao1_archaea$conf.int[1], 3), ", ", 
    round(cor_shannon_chao1_archaea$conf.int[2], 3), "]\n")
cat("n =", nrow(combined_data_archaea), "\n")

### ====================================================================================================================
### END OF SCRIPT
### ====================================================================================================================