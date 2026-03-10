## ====================================================================================================
## This script generates spider plots of environmental characteristics for microbial clades. 
##
## Author:       Dominic Eriksson
## Date:         8th of October 2025
## Affiliation:  Environmental Physics Group, UP, ETH Zürich, Switzerland
## Contact:      deriksson@ethz.ch
##
## Input files:
##   - EnvStatistics_inDivHotspots_Percentiles_v6.csv:
##       Code/4_Hotspots_diversity/1_Output/EnvStatistics_inDivHotspots_Percentiles_v6.csv
##
## Output files:
##   - Radar plots for each functional group (SVG) with clade-specific colors and standardized scales:
##       Radar_plot_group1_ggradar.svg
##       Radar_plot_group2_ggradar.svg
##       Radar_plot_group3_ggradar.svg
##       Radar_plot_group4_ggradar.svg
##
## Strategy:
##   1. Load environmental percentile data for hotspots per clade.
##   2. Filter for clades of interest.
##   3. Assign clades to four functional groups and define short names.
##   4. Define taxon-specific colors for plotting and clockwise predictor order.
##   5. Compute median values of selected environmental variables per clade.
##   6. Scale environmental variables to 0–1 for radar plotting.
##   7. Prepare data for ggradar: rename variables with friendly names and reorder columns.
##   8. Generate radar plots for each functional group using ggradar:
##       - Scaled axis values (0–1)
##       - Line and point sizes adjusted
##       - Group-specific colors
##       - Legends and titles with publication-ready fonts
##   9. Save plots as SVG files.
##
## Required R packages (tested versions):
##   - dplyr       (version 1.1.2)
##   - scales      (version 1.2.1)
##   - fmsb        (version 0.7.3)
##   - gridExtra   (version 2.3)
##   - grid        (base)
##   - ggplot2     (version 3.5.2)
##   - rnaturalearth (version 0.3.5)
##   - rnaturalearthdata (version 0.1.0)
##   - ggradar     (version 0.2.2)
##   - tibble      (version 3.2.1)
##
## ====================================================================================================


### ==============================================================================================
### Radar plots using Jonas groupings and ordered env. predictors clockwise according to physiochemical effects
### ===============================================================================================
# Load libraries
library(dplyr)
library(scales)
library(fmsb)
library(gridExtra)
library(grid)
library(ggplot2)
library(dplyr)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggradar)
library(tibble)

# Directories
wd_out <- "Code/Figures/Fig5/"

# Load data
df <- read.csv("Code/4_Hotspots_diversity/1_Output/EnvStatistics_inDivHotspots_Percentiles_v6.csv")

# Clades of interest - Based on figure 4a
clades <- c(
  "domain_Bacteria", # Group 2
  "domain_Archaea", # Group 4
  "class_Poseidoniia", # group 4
  "class_Nitrososphaeria", # group 1
  "class_Acidimicrobiia", # group 2
  "class_Bacteroidia", # group 4
  "class_Chlamydiia", # group 3
  "class_Cyanobacteriia", # group 2
  "class_UBA1144", # group 4
  "class_Fibrobacteria", # group 1  
  "class_Marinisomatia", # group 4
  "class_UBA4151", # group 3
  "class_UBA796", # group 3
  "class_XYA12-FULL-58-9", # group 3
  "class_Alphaproteobacteria", # group 2
  "class_Gammaproteobacteria", # group 4
  "class_SAR324" # group 4
)

# Create group assignments based on figure 4a
group_assignments <- data.frame(
  clade = clades,
  group = c(
    "group2",  # domain_Bacteria
    "group4",  # domain_Archaea
    "group4",  # Poseidoniia
    "group1",  # Nitrososphaeria
    "group2",  # Acidimicrobiia
    "group4",  # Bacteroidia
    "group3",  # Chlamydiia
    "group2",  # Cyanobacteriia
    "group4",  # UBA1144
    "group1",  # Fibrobacteria
    "group4",  # Marinisomatia
    "group3",  # UBA4151
    "group3",  # UBA796
    "group3",  # XYA12-FULL-58-9
    "group2",  # Alphaproteobacteria
    "group4",  # Gammaproteobacteria
    "group4"   # SAR324
  ),
  short_name = c(
    "Bacteria", "Archaea", "Poseidoniia", "Nitrososphaeria", "Acidimicrobiia", "Bacteroidia",
    "Chlamydiia", "Cyanobacteriia", "UBA1144", "Fibrobacteria",
    "Marinisomatia", "UBA4151", "UBA796", "XYA12-FULL-58-9",
    "Alphaproteobacteria", "Gammaproteobacteria", "SAR324"
  )
)

# Group labels based on new groupings
group_labels <- c(
  "group1" = "Group 1: Nitrososphaeria & Fibrobacteria",
  "group2" = "Group 2: Bacteria, Acidimicrobiia, Cyanobacteriia, Alphaproteobacteria",
  "group3" = "Group 3: Chlamydiia, UBA4151, UBA796, XYA12-FULL-58-9",
  "group4" = "Group 4: Archaea, Poseidoniia, Bacteroidia, UBA1144, Marinisomatia, Gammaproteobacteria, SAR324"
)

# Define distinct colors for each taxon within groups
taxon_colors <- c(
  # Group 1 colors (2 taxa)
  "Nitrososphaeria" = "#FF6B6B", "Fibrobacteria" = "#4ECDC4",
  # Group 2 colors (4 taxa)  
  "Bacteria" = "#000000fe", "Acidimicrobiia" = "#F1948A", "Cyanobacteriia" = "#60a278f7", "Alphaproteobacteria" = "#cfb445f7",
  # Group 3 colors (4 taxa)
  "Chlamydiia" = "#82E0AA", "UBA4151" = "#BB8FCE", "UBA796" = "#85C1E9", "XYA12-FULL-58-9" = "#F8C471",
  # Group 4 colors (7 taxa)
  "Archaea" = "#d27373ff", "Poseidoniia" = "#876f6fff", "Bacteroidia" = "#DDA0DD", "UBA1144" = "#96CEB4", 
  "Marinisomatia" = "#FFEAA7", "Gammaproteobacteria" = "#a6923bff", "SAR324" = "#F7DC6F"
)

# Environmental parameters: Temperature at top (12 o'clock), then anticlockwise
env_vars <- c(
  "climatology_t_0_50", # Temperature (12 o'clock)
  rev(c(
    "climatology_A_PAR_regridded",           # PAR
    "climatology_p_0_50",                    # Phosphate
    "climatology_n_0_50",                    # Nitrate
    "climatology_i_0_50",                    # Silicate
    "climatology_Nstar_0_50",                # N*
    "log10_climatology_TOT_POC_CMEMS",       # Total POC
    "climatology_A_CHLA_regridded",          # Chlorophyll-a
    "log10_climatology_S_PP_regridded",      # Primary Production
    "climatology_o_0_50",                    # Dissolved O2
    "climatology_O_0_50",                    # Oxygen Saturation
    "climatology_A_0_50",                    # AOU
    "climatology_M_0_0",                     # Mixed Layer Depth
    "log10_climatology_eke_aviso",           # Eddy Kinetic Energy
    "climatology_fsle_aviso_2001_2020"       # FSLE
  ))
)

# Adjust predictor names (keeping same friendly names)
column_name_map <- c(
  "climatology_t_0_50" = "Temperature",
  "climatology_A_PAR_regridded" = "PAR",
  "climatology_p_0_50" = "Phosphate",
  "climatology_n_0_50" = "Nitrate",
  "climatology_i_0_50" = "Silicate",
  "climatology_Nstar_0_50" = "N*",
  "log10_climatology_TOT_POC_CMEMS" = "Total POC",
  "climatology_A_CHLA_regridded" = "Chlorophyll-a",
  "log10_climatology_S_PP_regridded" = "Primary Production",
  "climatology_o_0_50" = "Dissolved O2",
  "climatology_O_0_50" = "Oxygen Saturation",
  "climatology_A_0_50" = "AOU",
  "climatology_M_0_0" = "Mixed Layer Depth",
  "log10_climatology_eke_aviso" = "Eddy Kinetic Energy",
  "climatology_fsle_aviso_2001_2020" = "FSLE"
)

# Filter data for clades of interest and add group information
df_filtered <- df %>%
  dplyr::filter(clade %in% clades) %>%
  left_join(group_assignments, by = "clade")

# Compute median values per individual clade for selected env variables
env_means_clades <- df_filtered %>%
  group_by(clade, group, short_name) %>%
  summarise(across(all_of(env_vars), median, na.rm = TRUE), .groups = "drop")

# Filter for Group 1
group1_data <- df_filtered %>% filter(group == "group1")

# Get world map data
world <- ne_countries(scale = "medium", returnclass = "sf")

# Plot
ggplot() +
  geom_sf(data = world, fill = "grey95", color = "grey60") +
  geom_point(data = group1_data, 
             aes(x = x, y = y, color = short_name), 
             size = 2, alpha = 0.7) +
  coord_sf(expand = FALSE) +
  scale_color_manual(values = c("Nitrososphaeria" = "#FF6B6B",
                                "Fibrobacteria" = "#4ECDC4")) +
  theme_minimal() +
  labs(title = "Richness Hotspots for Group 1 Clades",
       color = "Clade") +
  theme(legend.position = "bottom")


# Scale env variables (0-1) column-wise for radar plot
scaled_env_clades <- env_means_clades
scaled_env_clades[env_vars] <- lapply(env_means_clades[env_vars], scales::rescale)


# Prepare data for ggradar: use simplified names and correct order
friendly_names <- unname(column_name_map[env_vars])
ggradar_data <- scaled_env_clades %>%
  dplyr::select(group, short_name, dplyr::all_of(env_vars))
colnames(ggradar_data)[3:ncol(ggradar_data)] <- friendly_names

# ggradar expects the group column first, then label (for legend), then variables
ggradar_data <- ggradar_data %>%
  dplyr::rename(label = short_name)

# Plot for each group and save as SVG with appropriate legend and font sizes
group_order <- paste0("group", 1:4)
for (group_name in group_order) {
  group_df <- ggradar_data %>% filter(group == group_name)
  if (nrow(group_df) == 0) next
  # Move label column to first position, then variables in correct order
  plot_df <- group_df %>% dplyr::select(label, dplyr::all_of(friendly_names))
  # Get colors for this group in the order of label
  group_colors <- taxon_colors[as.character(plot_df$label)]
  # ggradar expects label column first
  p <- ggradar(
    plot_df,
    values.radar = c("0", "0.5", "1"),
    grid.min = 0, grid.mid = 0.5, grid.max = 1,
    group.line.width = 0.8,
    group.point.size = 1,
    legend.position = "right",
    legend.text.size = 10,
    axis.label.size = 4,
    font.radar = "Arial",
    group.colours = group_colors
  ) +
    ggtitle(paste0("Environmental Characteristics: ", group_labels[group_name])) +
    theme(
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10),
      plot.title = element_text(size = 13, face = "bold"),
      text = element_text(family = "Arial")
    )
  ggsave(filename = paste0(wd_out, "/Radar_plot_", group_name, "_ggradar.svg"), plot = p, width = 5, height = 5)
}

message("Saved ggradar radar charts for each group with proper scale labels.")

## ====================================================================================================
## END OF SCRIPT
## ====================================================================================================
