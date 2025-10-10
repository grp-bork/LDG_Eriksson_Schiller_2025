## ====================================================================================================
## This script visualizes the global distribution of oceanic metagenomic samples and their environmental context.
##
## Author:       Dominic Eriksson
## Date:         8th of October, 2024
## Affiliation:  Environmental Physics Group UP, ETH Zürich, Switzerland
## Contact:      deriksson@ethz.ch
##
## Input files:  
##   - Metadata CSV: /Data/External/samples_contextual.xlsx
##   - Alpha diversity data in folder: /Data/Internal/01_alpha_diversity_clade_level_JS/alpha_diversity_prokaryotes_domain_phylum_class_rarefied5000.csv
##   - Ocean regions shapefile: GOaS_v1_20211214/goas_v01.shp for Pacific centered map
##
## Output files: 
##   - Global maps of sampling effort
##   - Scatterplots of depth vs longitude and latitude
##   - Hovmoeller plots of samples by month and season
##
## Strategy:
##   This script provides a comprehensive visualization of the spatial, temporal and vertical distribution of oceanic
##   metagenomic samples highlighting sampling coverage.
##   1. Load sample metadata and alpha diversity indices
##   2. Prepare Pacific-centered global maps of sample locations using sf and rnaturalearth data.
##   3. Merge alpha diversity metrics with metadata and clean for missing values.
##   4. Compute summary statistics of samples per ocean basin and visualize distributions with treemaps.
##   5. Generate scatterplots of sample depth vs longitude and latitude, annotated by analysis type (MLD or mesopelagic).
##   6. Aggregate samples into latitude bins and calculate seasonal counts for Hovmoeller plots.
##
## Required R packages (tested versions):
##   - dplyr        (version 1.1.2)
##   - ggplot2      (version 3.5.2)
##   - data.table   (version 1.15.4)
##   - sf           (version 1.0.14)
##   - rnaturalearth (version 0.3.9)
##   - viridis      (version 0.6.2)
##   - raster       (version 3.6.26)
##   - terra        (version 1.7.78)
##   - tidyterra    (version 0.4.0)
##   - scales       (version 1.2.1)
##   - ncdf4        (version 1.6.0)
##   - treemap      (version 2.4-3)
##   - plyr         (version 1.8.8)
##
## ====================================================================================================


### ====================================================
### PREPARE ENVIRONMENT
### ====================================================
# Clear workspace
rm(list = ls())

# Libraries

# Libraries
library(dplyr)
library(ggplot2)
library(treemap)
library(data.table)
library(sf)
library(rnaturalearth)
library(viridis)
library(raster)
library(terra)
library(tidyterra)
library(scales)
library(ncdf4)
library(plyr)
library(readxl)


# Directories - Set your working dirctory accordingly
wd_out <- "Code/Figures/Fig1/"

# Check if folder exists, if not create it (including all parent directories)
if (!dir.exists(wd_out)) {
  dir.create(wd_out, recursive = TRUE, showWarnings = FALSE)
  message("Created directory: ", wd_out)
} else {
  message("Directory already exists: ", wd_out)
}

## Pacific centered global map
# Read shapefile with geometries of ocean basins
sf.ocean <- sf::st_read("Data/Ocean_regions_shapefile/GOaS_v1_20211214/goas_v01.shp") # direct to the folder where the ocean shapefile is saved, downloadable from marine.org (Global Oceans and seas)

# Pacific centered world map: Source: https://stackoverflow.com/questions/56146735/visual-bug-when-changing-robinson-projections-central-meridian-with-ggplot2
worldMap <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
  sf::st_make_valid()

# Set projection pacific centered
lon_pacific <- '200' 
target_crs <- sf::st_crs( paste0("+proj=eqc +x_0=0 +y_0=0 +lat_0=0 +lon_0=", lon_pacific) )

# define a long & slim polygon that overlaps the meridian line & set its CRS to match that of world Centered in lon 200
offset <- 180 - as.numeric(lon_pacific)

polygon <- sf::st_polygon(x = list(rbind(
  c(-0.0001 - offset, 90),
  c(0 - offset, 90),
  c(0 - offset, -90),
  c(-0.0001 - offset, -90),
  c(-0.0001 - offset, 90)
))) %>%
  sf::st_sfc() %>%
  sf::st_set_crs(4326)


# modify world dataset to remove overlapping portions with world's polygons
world2 <- worldMap %>% sf::st_difference(polygon)
#> Warning: attribute variables are assumed to be spatially constant throughout all
#> geometries
# Transform
sf.ocean_pacific <- world2 %>% sf::st_transform(crs = target_crs)


### ====================================================
### GLOBAL MAP OF SAMPLING EFFORT
### ====================================================

# Load data (Computed alpha diversity on samples and metadata)
df = as.data.table(read_excel("../Data/diversity_metrics.xlsx", sheet = 2))
# meta data
sheet2 = as.data.table(read_excel("../Data/External/samples_contextual.xlsx", sheet = 2))
sheet4 = as.data.table(read_excel("../Data/External/samples_contextual.xlsx", sheet = 4))
df_metadata = merge(sheet2, sheet4, by = "biosample", all = TRUE)


# Add lon and lat columns from the metadata to df
df <- df%>%
  dplyr::left_join( df_metadata %>% dplyr::select(biosample, lat, lon, analysis_type, depth), by = "biosample")

# Remove NAs in lon and lat columns
df <- df[!is.na(df$lat) & !is.na(df$lon), ]

# Subset for analysis_type
df <- df %>% dplyr::filter(analysis_type != "not_included")

# Remove analysis type
df <- df %>% dplyr::select(-analysis_type)

# Remove NAs in prokaryotes_Richness_5000
df <- df[!is.na(df$prokaryotes_Richness_5000), ]

## Get statistics of fraction of samples in ocean basins

# Load function to assign additional column for ocean regions
source("Data/Functions/subset_ocean_regions.R")

# We show sample distribution of the complete prokaryotic dataset
fractions <- df
fractions <- df[, c("lon", "lat", "prokaryotes_Richness_5000")]
names(fractions) <- c("decimalLongitude", "decimalLatitude", "prokaryotes_Richness_5000")

# Apply function
fractions <- subset_ocean_regions(fractions)

# Rename column
names(fractions)[1] <- "prokaryotic_richness"
summary(fractions)

# Count the number of samples in each marine region
region_counts <- fractions %>%
  dplyr::group_by(marine_region) %>%
  dplyr::summarise(n_samples = n()) %>%
  dplyr::mutate(percentage = 100 * n_samples / sum(n_samples)) %>%
  dplyr::arrange(desc(n_samples))
print(region_counts)

# Convert to spatial object for plotting
df_sf <- sf::st_as_sf(df, coords = c("lon", "lat"), crs = 4326)

# Subset the data if neccassary
df_sub <- df_sf["prokaryotes_Richness_5000"]

# Define custom colors
custom_colors <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d")

# Print world map
gg_plot <- ggplot(sf.ocean_pacific) + 
  geom_sf(fill = "grey") +
  coord_sf(expand = FALSE) +
  geom_sf(data = df_sub, fill = "black", shape = 21, size = 0.25, color = "transparent") +
  scale_fill_viridis(option = "viridis") +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 6, face = "bold"),   # Smaller legend title
    legend.text = element_text(size = 6),                   # Smaller legend text
    legend.key.size = unit(0.4, "cm"),                      # Smaller legend key
    legend.position = "bottom",
    plot.margin = margin(10, 10, 10, 10)
  ) +
  guides(fill = guide_legend(
    title = "Prok. richness",
    title.position = "top",
    title.hjust = 0.5,
    keywidth = unit(0.8, "cm"),   # Decrease key width
    keyheight = unit(0.4, "cm"),  # Decrease key height
    label.position = "bottom",
    label.hjust = 0.5,
    nrow = 1
  )) +
  ggtitle(paste0("Global distribution of samples (n = ", nrow(df_sub), ")" ))

# Print the plot
print(gg_plot)

# Save as before
ggsave(
  paste0(wd_out, "Global_distribution_of_samples_and_Prok_richness_5000_v3.svg"),
  gg_plot,
  width = 3, height = 3,
  dpi = 300
)

### ====================================================
### SCATTERPLOT - DEPTH VS. LONGITUDE
### ====================================================

# Load data
df = as.data.table(read_excel("../Data/diversity_metrics.xlsx", sheet = 2))
# meta data
sheet2 = as.data.table(read_excel("../Data/External/samples_contextual.xlsx", sheet = 2))
sheet4 = as.data.table(read_excel("../Data/External/samples_contextual.xlsx", sheet = 4))
df_metadata = merge(sheet2, sheet4, by = "biosample", all = TRUE)

# Add lon and lat columns from the metadata to df
df <- df %>%
  left_join( df_metadata %>% dplyr::select(biosample, lat, lon, analysis_type, depth), by = "biosample")

# Remove NAs in lon and lat columns
df <- df[!is.na(df$lat) & !is.na(df$lon), ]

# Subset for analysis_type
df <- df %>% filter(analysis_type != "not_included")

# Remove NAs in prokaryotes_Richness_5000
df <- df[!is.na(df$prokaryotes_Richness_5000), ]

# Subset dataframe
df_sub <- df[, c("lon", "lat", "depth", "analysis_type")]


# Calculate the number of points for each analysis type
mld <- nrow(filter(df, analysis_type == "MLD"))
meso <- nrow(filter(df, analysis_type == "mesopelagic"))

# Create annotation data manually
annotation_data <- data.frame(
  analysis_type = c("Surface", "Mesopelagic"),
  n = c(mld, meso),
  lon = c(-150, -156),  # Position at the right edge of the plot
  depth = c(900, 1000)  # Position at the bottom of the plot
)

# Define custom colors for analysis types
custom_colors <- c("MLD" = "#bc2128", "mesopelagic" = "#146db4")

# Shift longitude for Pacific-centered plotting (centered at 200°E)
shift_longitude <- function(lon) {
  shifted <- lon
  shifted[lon < 0] <- shifted[lon < 0] + 360
  shifted <- shifted - 200
  # Wrap to -180:180 for plotting
  shifted[shifted > 180] <- shifted[shifted > 180] - 360
  shifted
}

df_sub$lon_pacific <- shift_longitude(df_sub$lon)

# Plot
gg_plot <- ggplot(data = df_sub, aes(x = lon_pacific, y = depth, color = analysis_type)) + 
  geom_point(size = 0.5) +
  scale_y_reverse() +
  theme_minimal(base_size = 8) +
  theme(
    axis.title.x = element_text(size = 8, face = "bold"),
    axis.title.y = element_text(size = 8, face = "bold"),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.ticks = element_line(size = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.line = element_line(size = 0.5),
    axis.title.x.top = element_text(size = 8, face = "bold"),
    axis.text.x.top = element_text(size = 8),
    axis.title.x.bottom = element_blank(),
    axis.text.x.bottom = element_blank(),
    legend.title = element_text(size = 6, face = "bold"),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.3, "cm"),
    legend.position = "bottom",
    plot.margin = margin(5, 5, 5, 5),
    panel.grid = element_blank(),  # Remove grid lines
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()  # Remove minor grid lines
  ) +
  labs(
    x = "Longitude (Pacific-centered, °)",
    y = "Depth (meters)",
    color = "Analysis Type"
  ) +
  scale_x_continuous(
    position = "top",
    limits = c(-180, 180),
    breaks = seq(-180, 180, by = 60)
  ) +
  scale_color_manual(values = custom_colors) +
  guides(color = guide_legend(
    title.position = "top",
    title.hjust = 0.5,
    keywidth = unit(0.5, "cm"),
    keyheight = unit(0.3, "cm"),
    label.position = "bottom",
    label.hjust = 0.5,
    nrow = 1
  )) +
  geom_text(
    data = annotation_data,
    aes(x = shift_longitude(lon), y = depth, label = paste0("n = ", n), color = analysis_type),
    size = 8 / .pt, hjust = 1.1, vjust = -1.1, show.legend = FALSE
  )

# Print the plot
print(gg_plot)

# Save high resolution plot
ggsave(
  paste0(wd_out, "Depth_distribution_with_custom_y_axis_intervals_longitudinal_v3.svg"),
  gg_plot,
  width =5, height = 1.5,
  dpi = 300
)


### ====================================================
### SCATTERPLOT - DEPTH VS. LATITUDE
### ====================================================

# Load data
df = as.data.table(read_excel("../Data/diversity_metrics.xlsx", sheet = 2))
# meta data
sheet2 = as.data.table(read_excel("../Data/External/samples_contextual.xlsx", sheet = 2))
sheet4 = as.data.table(read_excel("../Data/External/samples_contextual.xlsx", sheet = 4))
df_metadata = merge(sheet2, sheet4, by = "biosample", all = TRUE)

# Add lon and lat columns from the metadata to df
df <- df %>%
  left_join( df_metadata %>% dplyr::select(biosample, lat, lon, analysis_type, depth), by = "biosample")

# Remove NAs in lon and lat columns
df <- df[!is.na(df$lat) & !is.na(df$lon), ]

# Subset for analysis_type
df <- df %>% filter(analysis_type != "not_included")

# Remove NAs in prokaryotes_Richness_5000
df <- df[!is.na(df$prokaryotes_Richness_5000), ]

# Subset dataframe
df_sub <- df[, c("lon", "lat", "depth", "analysis_type")]


# Calculate the number of points for each analysis type
mld <- nrow(filter(df, analysis_type == "MLD"))
meso <- nrow(filter(df, analysis_type == "mesopelagic"))

# Create annotation data manually
annotation_data <- data.frame(
  analysis_type = c("Surface", "Mesopelagic"),
  n = c(mld, meso),
  lat = c(-70, -70),  # Position at the right edge of the plot
  depth = c(900, 1000)  # Position at the bottom of the plot
)

# Define custom colors for analysis types
custom_colors <- c("MLD" = "#bc2128", "mesopelagic" = "#146db4")

# Plot the data with flipped y-axis and x-axis at the bottom
gg_plot <- ggplot(data = df_sub, aes(x = lat, y = depth, color = analysis_type)) + 
  geom_point(size = 0.5) +
  scale_y_reverse() +
  scale_x_continuous(position = "bottom", limits = c(-65, 82), breaks = c(-50, 0, 50)) +
  theme_minimal(base_size = 8) +
  theme(
    axis.title.x = element_text(size = 8, face = "bold"),
    axis.title.y = element_text(size = 8, face = "bold"),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.ticks = element_line(size = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.line = element_line(size = 0.5),
    axis.title.x.top = element_blank(),
    axis.text.x.top = element_blank(),
    axis.title.x.bottom = element_text(size = 8, face = "bold"),
    axis.text.x.bottom = element_text(size = 8),
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 8),
    legend.position = "bottom",
    plot.margin = margin(5, 5, 5, 5),
    panel.grid = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Latitude (°)",
    y = "Depth (meters)",
    color = "Analysis Type"
  ) +
  scale_color_manual(values = custom_colors) +
  geom_text(
    data = annotation_data,
    aes(x = lat, y = depth, label = paste0("n = ", n), color = analysis_type),
    size = 8 / .pt, hjust = 1.1, vjust = -1.1, show.legend = FALSE
  )

print(gg_plot)

ggsave(
  paste0(wd_out, "Depth_distribution_with_custom_y_axis_intervals_latitudinal_v4.svg"),
  gg_plot,
  width = 2.6, height = 1.5,
  dpi = 300
)


### ====================================================
### Hovmoeller plots by seasons 
### ====================================================

# Load data
df = as.data.table(read_excel("../Data/diversity_metrics.xlsx", sheet = 2))
# meta data
sheet2 = as.data.table(read_excel("../Data/External/samples_contextual.xlsx", sheet = 2))
sheet4 = as.data.table(read_excel("../Data/External/samples_contextual.xlsx", sheet = 4))
df_metadata = merge(sheet2, sheet4, by = "biosample", all = TRUE)

# Add lon and lat columns from the metadata to df
df <- df%>%
  dplyr::left_join( df_metadata %>% dplyr::select(biosample, lat, lon, analysis_type, depth, month), by = "biosample")

# Remove NAs in lon and lat columns
df <- df[!is.na(df$lat) & !is.na(df$lon), ]

# Subset for analysis_type
df <- df %>% dplyr::filter(analysis_type != "not_included")

# Remove analysis type
df <- df %>% dplyr::select(-analysis_type)

# Remove NAs in prokaryotes_Richness_5000
df <- df[!is.na(df$prokaryotes_Richness_5000), ]
summary(df$prokaryotes_Richness_5000) # No NAs left

# Grid latitude
df$y <- round(df$lat + 0.5) - 0.5

# Subset data
df_sub <- df[, c("y", "month")]

# Define the size of the new latitude bins
latitude_bin_size <- 5  # 5° bins

# Aggregate the latitude bins
df_sub <- df_sub %>%
  dplyr::mutate(lat_bin = cut(y, breaks = seq(-80, 90, by = latitude_bin_size), include.lowest = TRUE))

# Define the seasons, accounting for Northern and Southern Hemispheres
df_sub <- df_sub %>%
  dplyr::mutate(season = case_when(
    (y >= 0 & month %in% c(12, 1, 2)) | (y < 0 & month %in% c(6, 7, 8)) ~ "Winter",
    (y >= 0 & month %in% c(3, 4, 5)) | (y < 0 & month %in% c(9, 10, 11)) ~ "Spring",
    (y >= 0 & month %in% c(6, 7, 8)) | (y < 0 & month %in% c(12, 1, 2)) ~ "Summer",
    (y >= 0 & month %in% c(9, 10, 11)) | (y < 0 & month %in% c(3, 4, 5)) ~ "Fall"
  ))


# Count the number of samples for each season in each aggregated latitude bin
df_binned_season <- df_sub %>%
  dplyr::group_by(lat_bin, season) %>%
  dplyr::summarise(sample_count = dplyr::n()) %>%
  dplyr::ungroup()

# Define the order of seasons
season_names <- c("Winter", "Spring", "Summer", "Fall")

# Remove NAs in season column
head(df_binned_season)
df_binned_season <- na.omit(df_binned_season)

hovmoeller_plot_season <- ggplot(df_binned_season, aes(x = season, y = lat_bin, fill = sample_count)) +
  geom_tile() +
  scale_fill_gradient(low = "lightgrey", high = "black", name = "Sample Count") +
  scale_x_discrete(limits = season_names) +
  labs(
    title = "Number of Samples by Latitude and Season",
    x = "Season",
    y = "Latitude Bin (°)"
  ) +
  theme_minimal(base_size = 8) +
  theme(
    text = element_text(size = 8),
    axis.text.x = element_text(size = 8, angle = 90, hjust = 1),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 8, face = "bold"),
    plot.title = element_text(size = 8, face = "bold"),
    legend.title = element_text(size = 5, face = "bold"),
    legend.text = element_text(size = 5),
    legend.key.size = unit(0.18, "cm"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Print the plot
print(hovmoeller_plot_season)

# Save output file
output_file <- paste0(wd_out, "Hovmoller_Plot_Latitude_Season_v6.svg")
ggsave(
  filename = output_file,
  plot = hovmoeller_plot_season,
  width = 1.5,
  height = 3,
  dpi = 300
)

## ====================================================================================================
## END OF SCRIPT
## ====================================================================================================
