## ====================================================================================================
## This script generates global consensus hotspot maps for microbial clades. 
## Hotspots are defined as cells exceeding the 75th percentile of mean richness per clade, and 
## consensus maps show discrete counts of hotspots in regard to clades present per grid cell.
##
## Author:       Dominic Eriksson
## Date:         8th of October 2025
## Affiliation:  Environmental Physics Group, UP, ETH Zürich, Switzerland
## Contact:      deriksson@ethz.ch
##
## Input files:
##   - Global_annual_richness_ensemble_mean_richness.csv
##       Code/3_Compute_ensembles_and_uncertainties/2_Output/
##   - Ocean regions shapefile (GOaS_v1_20211214/goas_v01.shp)
##       Data/Ocean_regions_shapefile/
##
## Output files:
##       Pacific-centered consensus hotspot maps for four functional groups, 
##       with discrete count of taxa present per grid cell.
##
## Strategy:
##   1. Load global richness data per clade and ocean shapefiles.
##   2. Prepare Pacific-centered projection by shifting longitudes and trimming polygons.
##   3. For each clade:
##       - Compute the 75th percentile threshold of mean richness.
##       - Define hotspot cells above this threshold.
##       - Store hotspot status along with coordinates and clade.
##   4. Assign clades to four functional groups.
##   5. Count number of taxa present in hotspot cells per group (discrete counts).
##   6. Prepare discrete color scales for counts and generate maps for each group using ggplot2.
##   7. Combine all four maps with patchwork, maintaining individual legends.
##   8. Save final combined map as high-resolution SVG.
##
## Required R packages (tested versions):
##   - ggplot2          (version 3.5.2)
##   - dplyr            (version 1.1.2)
##   - sf               (version 1.1-9)
##   - terra             (version 1.7-43)
##   - maps             (version 3.4.1)
##   - rnaturalearth    (version 0.3.5)
##   - rnaturalearthdata (version 0.1.0)
##   - tidyterra        (version 0.5.3)
##   - patchwork        (version 1.1.2)
##
## ====================================================================================================


# Clear workspace
rm(list = ls())

# Libraries
library(ggplot2)
library(maps)
library(dplyr)
library(ggplot2)
library(sf)  # for spatial plotting, if you want coastlines etc.

# Directories
wd_out <- "Code/Figures/Fig5/"

# Create output directory if it doesn't exist
if (!dir.exists(wd_out)) {
  dir.create(wd_out, recursive = TRUE)
}

# Load data
wd_in <- "Code/3_Compute_ensembles_and_uncertainties/2_Output/Global_annual_richness_ensemble_mean_richness.csv"
df_rich <- read.csv(wd_in)

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

# Define taxa of interest
unique(df_rich$clade)
clade <- c(
    "domain_Bacteria",
    "domain_Archaea",
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

# Create a data frame to store hotspot info for all clades
hotspot_df <- data.frame()

for (i in seq_along(clade)) {
  # Print progress
  print(paste0("Processing clade: ", clade[i], " (", i, "/", length(clade), ")"))
  
  # Convert to spatraster object
  df_rich$clade <- as.character(df_rich$clade)
  target_clade <- trimws(clade[i])
  df_div_sub <- dplyr::filter(df_rich, trimws(clade) == target_clade)
  dat <- terra::rast(df_div_sub[, c("x", "y", "mean_richness")])
  terra::crs(dat) <- 'epsg:4326'
  
  # Convert raster to data frame
  df_coords <- as.data.frame(dat, xy = TRUE, na.rm = TRUE)
  names(df_coords)[3] <- "mean_richness"
  
  # Calculate 75th percentile threshold
  threshold_75 <- quantile(df_coords$mean_richness, probs = 0.75, na.rm = TRUE)
  
  # Create hotspot column (1 if above 75th percentile, else 0)
  df_coords$hotspot <- as.integer(df_coords$mean_richness > threshold_75)
  
  # Add clade information
  df_coords$clade <- target_clade
  
  # Append to main data frame
  hotspot_df <- rbind(hotspot_df, df_coords[, c("x", "y", "hotspot", "clade")])
}

# Check the final hotspot data frame
head(hotspot_df)
summary(hotspot_df)
unique(hotspot_df$clade)

# Create group column based on clade names
hotspot_df$group <- NA

# Group 1: Nitrososphaeria & Fibrobacteria
hotspot_df$group[grepl("Nitrososphaeria", hotspot_df$clade)] <- "group1"
hotspot_df$group[grepl("Fibrobacteria", hotspot_df$clade)] <- "group1"

# Group 2: Bacteria, Acidimicrobiia, Cyanobacteriia, Alphaproteobacteria
hotspot_df$group[grepl("Bacteria", hotspot_df$clade)] <- "group2"
hotspot_df$group[grepl("Acidimicrobiia", hotspot_df$clade)] <- "group2"
hotspot_df$group[grepl("Cyanobacteriia", hotspot_df$clade)] <- "group2"
hotspot_df$group[grepl("Alphaproteobacteria", hotspot_df$clade)] <- "group2"

# Group 3: Chlamydiia, UBA4151, UBA796, XYA12-FULL-58-9
hotspot_df$group[grepl("Chlamydiia", hotspot_df$clade)] <- "group3"
hotspot_df$group[grepl("UBA4151", hotspot_df$clade)] <- "group3"
hotspot_df$group[grepl("UBA796", hotspot_df$clade)] <- "group3"
hotspot_df$group[grepl("XYA12-FULL-58-9", hotspot_df$clade)] <- "group3"

# Group 4: Archaea, Poseidoniia, Bacteroidia, UBA1144, Marinisomatia, Gammaproteobacteria, SAR324
hotspot_df$group[grepl("Archaea", hotspot_df$clade)] <- "group4"
hotspot_df$group[grepl("Poseidoniia", hotspot_df$clade)] <- "group4"
hotspot_df$group[grepl("Bacteroidia", hotspot_df$clade)] <- "group4"
hotspot_df$group[grepl("UBA1144", hotspot_df$clade)] <- "group4"
hotspot_df$group[grepl("Marinisomatia", hotspot_df$clade)] <- "group4"
hotspot_df$group[grepl("Gammaproteobacteria", hotspot_df$clade)] <- "group4"
hotspot_df$group[grepl("SAR324", hotspot_df$clade)] <- "group4"

# Check the group assignments
table(hotspot_df$group, useNA = "always")

# Calculate the total number of taxa in each group
taxa_per_group <- hotspot_df %>%
  filter(!is.na(group)) %>%
  group_by(group) %>%
  summarise(total_taxa = n_distinct(clade), .groups = "drop")

print("Taxa per group:")
print(taxa_per_group)

# Create consensus hotspot map with DISCRETE COUNT calculation (not percentages)
consensus_map_groups <- hotspot_df %>%
  filter(hotspot == 1) %>%
  group_by(group, x, y) %>%
  summarise(hotspot_count = n(), .groups = "drop") %>%
  filter(!is.na(group)) %>%
  left_join(taxa_per_group, by = "group")

# Check the count distribution
print("Count distribution:")
print(summary(consensus_map_groups$hotspot_count))

# Create breaks based on the maximum possible count for each group
max_counts <- taxa_per_group$total_taxa
overall_max <- max(max_counts)

# Create discrete breaks based on counts (0 to max possible taxa)
count_breaks <- 0:overall_max
count_colors <- colorRampPalette(c("#f7fbff", "#08306b"))(overall_max)

# Function to create plots with discrete count data
create_count_plot <- function(group_name, group_title) {
  plot_data <- consensus_map_groups %>% filter(group == group_name)
  
  if (nrow(plot_data) == 0) {
    return(NULL)
  }
  
  # Get the maximum count for this specific group
  group_max <- taxa_per_group$total_taxa[taxa_per_group$group == group_name]
  
  # Create specific breaks and colors for this group
  group_breaks <- 0:group_max
  # Use a more visible color palette that starts with light blue instead of white
  group_colors <- colorRampPalette(c("#c6dbef", "#08306b"))(group_max + 1)
  
  plot_raster <- terra::rast(plot_data[, c("x", "y", "hotspot_count")])
  terra::crs(plot_raster) <- "EPSG:4326"
  
  ggplot(sf.ocean_pacific) +
    geom_sf(fill = "grey", size = 0.1) +
    tidyterra::geom_spatraster_contour_filled(
      data = plot_raster, 
      show.legend = TRUE,
      breaks = group_breaks
    ) +
    scale_fill_manual(
      values = group_colors,
      name = "# Taxa Present",
      labels = as.character(group_breaks[-1]),  # Remove the 0 break for labels
      drop = FALSE
    ) +
    theme_minimal(base_size = 8) +
    theme(
      text = element_text(size = 8),
      legend.position = "right",
      legend.key.size = unit(0.3, 'cm'),
      legend.key.height = unit(0.25, 'cm'),
      legend.key.width = unit(0.25, 'cm'),
      legend.title = element_text(size = 8, hjust = 0.5),
      legend.text = element_text(size = 7),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      plot.title = element_text(size = 8, face = "bold", hjust = 0.5)
    ) +
    labs(title = group_title) +
    xlab("") + ylab("")
}

# Create plots for each group (each with its own legend)
plot1 <- create_count_plot("group1", "Group 1: Nitrososphaeria & Fibrobacteria")
plot2 <- create_count_plot("group2", "Group 2: Bacteria, Acidimicrobiia, Cyanobacteriia, Alphaproteobacteria")
plot3 <- create_count_plot("group3", "Group 3: Chlamydiia, UBA4151, UBA796, XYA12-FULL-58-9")
plot4 <- create_count_plot("group4", "Group 4: Archaea, Poseidoniia, Bacteroidia, UBA1144, Marinisomatia, Gammaproteobacteria, SAR324")

# Print individual plots
print(plot1)
print(plot2)
print(plot3)
print(plot4)

# Combine all 4 plots - each keeps its own legend
library(patchwork)
combined_4groups_counts <- plot1 / plot2 / plot3 / plot4

print(combined_4groups_counts)

# Save as high resolution svg image
ggsave(filename = paste0(wd_out, "Global_consensus_hotspot_maps_4groups_counts_v1.svg"),
       plot = combined_4groups_counts,
       width = 10, height = 14, dpi = 300, device = "svg")


## ====================================================================================================
## END OF SCRIPT
## ====================================================================================================
