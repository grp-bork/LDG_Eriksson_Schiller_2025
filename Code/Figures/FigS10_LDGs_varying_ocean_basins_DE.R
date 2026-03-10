## ====================================================================================================
## This R script performs comprehensive analysis of latitudinal diversity gradients (LDGs) for
## prokaryotic richness across different ocean basins. The analysis examines prokaryotic diversity
## patterns along longitude transects in four major ocean regions: Eastern Pacific (130°W),
## Western Pacific (160°E), Eastern Atlantic (20°W), and Central Atlantic (30°W). The script generates
## publication-ready visualizations including a global map showing transect locations and
## latitudinal richness profiles for each ocean basin.
##
## Author:       Dominic Eriksson
## Date:         27th of February 2026
## Affiliation:  Environmental Physics Group UP, ETH Zürich, Switzerland
## Contact:      deriksson@ethz.ch
##
## Input files:
##   - Ensemble members richness CSV with monthly prokaryotic diversity predictions
##   - Natural Earth country boundaries for global map visualization
##   - Longitude transect configuration defining ocean basin locations
##
## Output files:
##   - Multi-panel SVG figure combining global map and latitudinal richness profiles
##   - Individual SVG files for map and richness plots (optimized for publication)
##   - Basin summary statistics CSV with richness metrics by ocean region
##   - Comparative analysis data for prokaryotic diversity patterns
##
## Strategy:
##   The script applies modular configuration for easy modification of longitude bins and basin
##   parameters across four major ocean transects. Ensemble statistics calculate mean richness
##   and uncertainty bands for each latitude-longitude combination, enabling robust assessment
##   of latitudinal diversity gradients. Global map visualization uses colored longitude transects
##   while latitudinal profiles display richness patterns with uncertainty ribbons. Publication-ready
##   output optimized for 5×6 inch dimensions facilitates comparative analysis of prokaryotic
##   diversity patterns across ocean basins.
##
## Required R packages (tested versions):
##   - data.table      1.14.8
##   - dplyr           1.1.2
##   - tidyr           1.3.0
##   - ggplot2         3.4.2
##   - sf              1.0.12
##   - cowplot         1.1.1
##   - rnaturalearth   0.3.3
##   - rnaturalearthdata 0.1.0
## ====================================================================================================


### ====================================================================================================================
### Preparation
### ====================================================================================================================

# Clear workspace
rm(list = ls())

# Libraries
library(data.table)     # For fast data reading and manipulation
library(dplyr)          # For data manipulation and filtering operations
library(tidyr)          # For data reshaping and tidying functions
library(ggplot2)        # For creating publication-quality plots and visualizations
library(sf)             # For handling spatial data and geographic projections
library(cowplot)        # For combining multiple plots into publication-ready figures
library(rnaturalearth)  # For accessing natural earth map data
library(rnaturalearthdata) # Additional natural earth datasets

### ====================================================================================================================
### Configuration Section - Easily Change These Parameters
### ====================================================================================================================

# Define longitude bins to analyze (CHANGE HERE for different bins)
LONGITUDE_BINS <- c(-130, 160, -20, -30)

# Define basin names and colors (ADD/MODIFY HERE for new bins)
BASIN_CONFIG <- data.frame(
    longitude = c(-130, 160, -20, -30),
    basin_name = c("Eastern_Pacific", "Western_Pacific", "Eastern_Atlantic", "Central_Atlantic"),
    basin_label = c("Eastern Pacific", "Western Pacific", "Eastern Atlantic", "Central Atlantic"),
    lon_label = c("130W", "160E", "20W", "30W"),  # Fixed to match actual longitudes
    color = c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"),  # Blue, Orange, Green, Red
    stringsAsFactors = FALSE
)

# Order basins from west to east for consistent plotting
BASIN_ORDER <- c("Eastern_Pacific", "Central_Atlantic", "Eastern_Atlantic", "Western_Pacific")

# Directories
wd_out <- "./Output/"
wd_in <- "../4_Compute_ensembles_uncertainties/1_Output/df_long_monthly_ensembleMembers_richness_5000_v2.csv"

# Create output directory
if (!dir.exists(wd_out)) {
    dir.create(wd_out, recursive = TRUE)
}

### ====================================================================================================================
### Helper Functions
### ====================================================================================================================

# Function to create basin labels and colors dynamically
create_basin_mapping <- function(longitude_bins, basin_config) {
    # Filter config for selected longitude bins
    config_subset <- basin_config[basin_config$longitude %in% longitude_bins, ]
    
    # Create named vectors for colors and labels
    colors <- setNames(config_subset$color, config_subset$basin_name)
    labels <- setNames(config_subset$basin_label, config_subset$basin_name)
    
    return(list(colors = colors, labels = labels, config = config_subset))
}

# Function to add basin information to data
add_basin_info <- function(data, basin_config) {
    data %>%
        dplyr::left_join(
            basin_config %>% dplyr::select(longitude, basin_name, basin_label, lon_label),
            by = c("x_rounded_1" = "longitude")
        ) %>%
        dplyr::mutate(
            # Fallback for any unmatched longitudes
            basin_name = ifelse(is.na(basin_name), paste0("Lon_", x_rounded_1), basin_name),
            basin_label = ifelse(is.na(basin_label), paste0("Longitude ", x_rounded_1, "°"), basin_label),
            lon_label = ifelse(is.na(lon_label), 
                             paste0(abs(x_rounded_1), ifelse(x_rounded_1 < 0, "W", "E")), 
                             lon_label),
            # Order basins consistently
            basin_name = factor(basin_name, levels = BASIN_ORDER)
        )
}

### ====================================================================================================================
### Data Loading and Processing
### ====================================================================================================================

cat("Loading and processing data...\n")

# Load world map
world <- ne_countries(scale = "medium", returnclass = "sf")

# Create longitude lines for map
longitude_lines <- data.frame(
    lon = seq(-180, 180, by = 10),
    lat_start = -90,
    lat_end = 90
)

# Load data
df_long <- data.table::fread(wd_in)

# Round longitude to nearest whole degree
df_long$x_rounded_1 <- round(df_long$x)

# Get basin mapping
basin_mapping <- create_basin_mapping(LONGITUDE_BINS, BASIN_CONFIG)

# Filter data for selected longitude bins
df_long_filtered <- df_long %>%
    filter(x_rounded_1 %in% LONGITUDE_BINS)

cat("Available longitude bins in filtered data:\n")
print(sort(unique(df_long_filtered$x_rounded_1)))

### ====================================================================================================================
### Calculate Ensemble Statistics
### ====================================================================================================================

cat("Calculating ensemble statistics...\n")

# Function to calculate statistics for a domain
calculate_ldg_stats <- function(data, domain_name) {
    data %>%
        dplyr::filter(clade == paste0("domain_", domain_name)) %>%
        dplyr::group_by(y, x_rounded_1) %>%
        dplyr::summarise(
            mean_richness = mean(richness, na.rm = TRUE),
            sd_richness = sd(richness, na.rm = TRUE),
            n_models = n(),
            .groups = "drop"
        ) %>%
        dplyr::mutate(
            lower_bound = mean_richness - sd_richness,
            upper_bound = mean_richness + sd_richness
        ) %>%
        add_basin_info(basin_mapping$config) %>%
        dplyr::mutate(domain = domain_name)
}


# Function to calculate statistics for prokaryotes
calculate_prokaryotes_ldg_stats <- function(data) {
    data %>%
        dplyr::filter(clade == "prokaryotes") %>%
        dplyr::group_by(y, x_rounded_1) %>%
        dplyr::summarise(
            mean_richness = mean(richness, na.rm = TRUE),
            sd_richness = sd(richness, na.rm = TRUE),
            n_models = n(),
            .groups = "drop"
        ) %>%
        dplyr::mutate(
            lower_bound = mean_richness - sd_richness,
            upper_bound = mean_richness + sd_richness
        ) %>%
        add_basin_info(basin_mapping$config) %>%
        dplyr::mutate(domain = "Prokaryotes")
}

# Calculate for combined prokaryotes
prokaryotes_ldg_by_longitude <- calculate_prokaryotes_ldg_stats(df_long_filtered)

# Summary statistics
basin_summary <- prokaryotes_ldg_by_longitude %>%
    dplyr::group_by(domain, basin_name, basin_label, x_rounded_1) %>%
    dplyr::summarise(
        n_latitudes = n(),
        mean_richness_overall = mean(mean_richness, na.rm = TRUE),
        max_richness = max(mean_richness, na.rm = TRUE),
        min_richness = min(mean_richness, na.rm = TRUE),
        latitude_max_richness = y[which.max(mean_richness)],
        .groups = "drop"
    )

cat("\nSummary by ocean basin for prokaryotes:\n")
print(basin_summary)

### ====================================================================================================================
### Create Optimized Plots for Publication
### ====================================================================================================================

cat("Creating optimized plots...\n")

# Create data for colored longitude lines
long_bins_colored <- basin_mapping$config %>%
    dplyr::select(longitude, basin_name, basin_label) %>%
    dplyr::mutate(
        lon = longitude,
        lat_start = -90,
        lat_end = 90,
        basin_name = factor(basin_name, levels = BASIN_ORDER)
    )

# Optimized map plot
p_map_optimized <- ggplot() +
    # Add world map
    geom_sf(data = world, fill = "lightgray", color = "white", size = 0.1) +
    # Add ALL vertical longitude lines (every 10°) in light gray
    geom_segment(
        data = longitude_lines,
        aes(x = lon, xend = lon, y = lat_start, yend = lat_end),
        color = "lightgray",
        size = 0.2,
        alpha = 0.5
    ) +
    # Add colored longitude lines for selected basins
    geom_segment(
        data = long_bins_colored,
        aes(x = lon, xend = lon, y = lat_start, yend = lat_end, color = basin_name),
        size = 1,
        alpha = 1
    ) +
    # Add longitude labels
    geom_text(
        data = long_bins_colored,
        aes(x = lon, y = -100, label = paste0(lon, "°"), color = basin_name),
        size = 2,  # Reduced text size
        hjust = 0.5,
        vjust = 1,
        fontface = "bold"
    ) +
    # Apply colors and labels
    scale_color_manual(
        values = basin_mapping$colors,
        labels = basin_mapping$labels,
        name = "Ocean Basin"
    ) +
    # Map settings
    coord_sf(
        xlim = c(-180, 180),
        ylim = c(-90, 90),
        crs = st_crs(4326),
        expand = FALSE
    ) +
    # Styling optimized for small dimensions
    theme_minimal(base_size = 8) +
    theme(
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
        panel.grid.major = element_line(color = "lightblue", size = 0.1),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "lightblue", color = NA),
        legend.position = "bottom",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.3, "cm"),
        plot.margin = margin(2, 2, 2, 2, "pt")
    ) +
    labs(
        title = "Ocean Basin Longitude Transects",
        x = "Longitude (°)",
        y = "Latitude (°)"
    )

# Optimized Prokaryotes plot
p_prokaryotes_optimized <- ggplot(prokaryotes_ldg_by_longitude, aes(x = mean_richness, y = y)) +
    # Add ribbons for standard deviation
    geom_ribbon(
        aes(xmin = lower_bound, xmax = upper_bound, fill = basin_name),
        alpha = 0.4,
        color = NA
    ) +
    # Add mean richness lines
    geom_path(
        aes(color = basin_name),
        size = 0.5
    ) +
    # Facet by basin (ordered west to east)
    facet_wrap(~ basin_name, scales = "free_x", ncol = 4) +
    # Apply colors
    scale_color_manual(values = basin_mapping$colors) +
    scale_fill_manual(values = basin_mapping$colors) +
    # Y-axis settings
    scale_y_continuous(
        breaks = c(-60, -30, 0, 30, 60),  # Fewer breaks for compact display
        labels = c("-60", "-30", "0", "30", "60")
    ) +
    # Labels and theme optimized for small size
    labs(
        title = "Prokaryotic Richness",
        x = "Richness",
        y = "Latitude (°)"
    ) +
    theme_minimal(base_size = 8) +
    theme(
        plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 8, face = "bold"),
        strip.text = element_text(size = 7, face = "bold"),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey90", size = 0.2),
        plot.margin = margin(2, 2, 2, 2, "pt"),
        panel.spacing = unit(0.1, "cm")
    )

### ====================================================================================================================
### Create Multiplot
### ====================================================================================================================

cat("Creating multiplot...\n")

# Create final optimized multiplot
multiplot_final <- plot_grid(
    p_map_optimized,
    p_prokaryotes_optimized,
    ncol = 1,
    nrow = 2,
    rel_heights = c(0.7, 1),
    align = "v",
    axis = "lr",
    labels = c("A", "B"),
    label_size = 8,
    label_fontface = "bold"
)

# Display plots
print(p_map_optimized)
print(p_prokaryotes_optimized)
print(multiplot_final)

### ====================================================================================================================
### Save Outputs
### ====================================================================================================================

cat("Saving outputs...\n")

# Save optimized multiplot with exact dimensions (adjusted for 2-panel layout)
ggsave(
    paste0(wd_out, "Multiplot_Prokaryotes_Optimized_5x6_v3.svg"), 
    multiplot_final, 
    width = 5, 
    height = 6, 
    units = "in", 
    dpi = 300
)

# Save individual optimized plots
ggsave(
    paste0(wd_out, "Map_Ocean_Basin_Transects_optimized_v3.svg"), 
    p_map_optimized, 
    width = 5, height = 2, units = "in", dpi = 300
)

ggsave(
    paste0(wd_out, "Prokaryotes_LDG_optimized_v3.svg"), 
    p_prokaryotes_optimized, 
    width = 5, height = 3, units = "in", dpi = 300
)

# Save summary data
write.csv(basin_summary, 
          paste0(wd_out, "Basin_Summary_Statistics_v3.csv"), 
          row.names = FALSE)

### ====================================================================================================================
### Print Summary
### ====================================================================================================================

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("ANALYSIS COMPLETE\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("Longitude bins analyzed:", paste(LONGITUDE_BINS, collapse = ", "), "\n")
cat("Number of basins:", nrow(basin_mapping$config), "\n")
cat("Output directory:", wd_out, "\n")
cat("\nFiles created:\n")
cat("- Multiplot_Prokaryotes_Optimized_5x6_v2.svg (MAIN OUTPUT - 5x6 inches)\n")
cat("- Map_Ocean_Basin_Transects_optimized_v2.svg\n")
cat("- Prokaryotes_LDG_optimized_v2.svg\n")
cat("- Basin_Summary_Statistics_v2.csv\n")
cat("\nOptimized for 8pt font and 5x6 inch dimensions\n")
cat("To change longitude bins, modify LONGITUDE_BINS and BASIN_CONFIG at the top.\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# Print prokaryotic richness statistics
cat("\nProkaryotic richness summary by basin:\n")
basin_stats <- prokaryotes_ldg_by_longitude %>%
    group_by(basin_name, basin_label, lon_label) %>%
    summarise(
        n_latitudes = n(),
        mean_richness_basin = mean(mean_richness, na.rm = TRUE),
        max_richness_basin = max(mean_richness, na.rm = TRUE),
        min_richness_basin = min(mean_richness, na.rm = TRUE),
        latitude_max_richness = y[which.max(mean_richness)],
        .groups = "drop"
    ) %>%
    arrange(desc(mean_richness_basin))

print(basin_stats)

cat("\nOptimized multiplot saved as 'Multiplot_Prokaryotes_Optimized_5x6_v2.svg'\n")
cat("Dimensions: 5 inches wide x 6 inches tall\n")
cat("Font size: 8pt base with smaller elements at 6-7pt\n")


## ====================================================================================================
## END OF SCRIPT
## ====================================================================================================