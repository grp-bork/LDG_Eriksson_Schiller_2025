## ====================================================================================================
## This script generates latitudinal diversity gradients (LDG) for modeled bacterial and archaeal richness
## comparing two months (June vs December) and highlights latitude bins with statistically significant differences.
##
## Author:       Dominic Eriksson
## Date:         8th of October 2025
## Affiliation:  Environmental Physics Group, UP, ETH Zürich, Switzerland
## Contact:      deriksson@ethz.ch
##
## Input files:
##   - Monthly alpha-diversity CSV: Code/3_Compute_ensembles_and_uncertainties/5_Output/df_latGradients_alphaDiv_month_v5.csv
##   - Significance summary CSV: Code/5_Wilcox_test/1_Output/bin_significance_summary_all_bins.csv
##
## Output files:
##   - Lineplots of latitudinal richness for Bacteria and Archaea (June vs December)
##
## Strategy:
##   1. Load monthly richness data for bacterial and archaeal clades.
##   2. Reshape data to obtain mean and standard deviation per latitude and month.
##   3. Compute upper and lower bounds (mean ± SD) for ribbon plots.
##   4. Identify maximum mean richness for each clade and annotate on the plot.
##   5. Load significance summary and highlight latitude bins with >50% model significance.
##   6. Create LDG plots with:
##       - Colored ribbons for June and December
##       - Lines for mean richness per month
##       - Skyblue rectangles for significant bins
##   7. Save plots as high-resolution SVG files.
##
## Required R packages (tested versions):
##   - data.table   (version 1.15.4)
##   - dplyr        (version 1.1.2)
##   - tidyr        (version 1.3.1)
##   - ggplot2      (version 3.5.2)
##
## ====================================================================================================


# Libraries
# Required libraries
library(data.table)  # for fread()
library(dplyr)       # for %>%, filter(), select(), mutate()
library(tidyr)       # for pivot_longer(), pivot_wider()
library(ggplot2)     # for ggplot(), geom_ribbon(), geom_path(), etc.

# Directories - adjust as needed
wd_out <- "Code/Figures/Fig3/"

# Create output directory if it doesn't exist
if (!dir.exists(wd_out)) {
  dir.create(wd_out, recursive = TRUE)
}

# Load data
df <- fread("Code/3_Compute_ensembles_and_uncertainties/5_Output/df_latGradients_alphaDiv_month_v5.csv")
segment_data <- data.table::fread("Code/5_Wilcox_test/1_Output/bin_significance_summary_all_bins.csv")

## LDG plot for domain_Bacteria with significance segments ----------------

# Extract domain_Bacteria data for June and December
bacteria_data <- df %>%
    filter(clade == "domain_Bacteria") %>%
    select(y, stat_type, June, December) %>%
    pivot_longer(
        cols = c(June, December),
        names_to = "Month",
        values_to = "richness"
    ) %>%
    filter(!is.na(richness)) %>%
    pivot_wider(
        names_from = stat_type,
        values_from = richness
    ) %>%
    filter(!is.na(mean) & !is.na(sd))

# Get maximum richness for bacteria
max_mean_diversity_bacteria <- max(bacteria_data$mean, na.rm = TRUE)

# Prepare segments for Bacteria (>50% significance)
bacteria_segments <- segment_data %>%
    filter(clade == "domain_Bacteria" & fraction_models_significant > 0.5) %>%
    mutate(
        segment_x = max_mean_diversity_bacteria * 1.02,
        segment_xend = max_mean_diversity_bacteria * 1.05
    )

# Define y-axis labels
y_labels_to_show <- c(-75, -50, -25, 0, 25, 50, 75)

# Plot for Bacteria
p_bacteria <- ggplot(bacteria_data, aes(x = mean, y = y, group = Month)) +
  # Add skyblue rectangles for significant bins FIRST (behind other elements)
  geom_rect(
    data = bacteria_segments,
    aes(xmin = segment_x, xmax = segment_xend, ymin = y - 0.5, ymax = y + 0.5),
    fill = "skyblue",
    color = NA,
    inherit.aes = FALSE
  ) +
  # December ribbon
  geom_ribbon(
    data = bacteria_data %>% filter(Month == "December"),
    aes(
      y = y,
      xmin = mean - sd,
      xmax = mean + sd
    ),
    fill = "darkgrey",
    alpha = 0.5,
    inherit.aes = FALSE,
    color = NA
  ) +
  # June ribbon
  geom_ribbon(
    data = bacteria_data %>% filter(Month == "June"),
    aes(
      y = y,
      xmin = mean - sd,
      xmax = mean + sd
    ),
    fill = "deeppink",
    alpha = 0.3,
    inherit.aes = FALSE,
    color = NA
  ) +
  # December line
  geom_path(
    data = bacteria_data %>% filter(Month == "December"),
    aes(color = Month),
    size = 0.5
  ) +
  # June line
  geom_path(
    data = bacteria_data %>% filter(Month == "June"),
    aes(color = Month),
    size = 0.5
  ) +
  scale_color_manual(
    values = c("December" = "#000000", "June" = "#FF00FF"),
    name = "Month"
  ) +
  scale_y_continuous(
    breaks = y_labels_to_show,
    labels = as.character(y_labels_to_show)
  ) +
  labs(
    title = "Bacterial richness - December vs June",
    x = "Richness",
    y = "Latitude (°)",
    color = "Month"
  ) +
  theme_minimal(base_size = 8) +
  theme(
    plot.title = element_text(size = 8, face = "bold", hjust = 0.5, vjust = 1),
    axis.text.x = element_text(size = 8, angle = 90, hjust = 1),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.3, "cm"),
    legend.key.height = unit(0.3, "cm"),
    legend.key.width = unit(0.3, "cm"),
    panel.grid = element_blank(),
    legend.position = "none"
  )

# Display plot for Bacteria
print(p_bacteria)

# Save plot
ggsave(
  paste0(wd_out, "LatGradients_DecemberJune_Bacteria_withSegments_v2.svg"), 
  p_bacteria, 
  width = 1,
  height = 2.2, 
  units = "in",
  dpi = 300
)


## LDG plot for domain_Archaea with significance segments ----------------

# Extract domain_Archaea data for June and December
archaea_data <- df %>%
    filter(clade == "domain_Archaea") %>%
    select(y, stat_type, June, December) %>%
    pivot_longer(
        cols = c(June, December),
        names_to = "Month",
        values_to = "richness"
    ) %>%
    filter(!is.na(richness)) %>%
    pivot_wider(
        names_from = stat_type,
        values_from = richness
    ) %>%
    filter(!is.na(mean) & !is.na(sd))

# Get maximum richness for archaea
max_mean_diversity_archaea <- max(archaea_data$mean, na.rm = TRUE)

# Prepare segments for Archaea (>50% significance)
archaea_segments <- segment_data %>%
    filter(clade == "domain_Archaea" & fraction_models_significant > 0.5) %>%
    mutate(
        segment_x = max_mean_diversity_archaea * 1.02,
        segment_xend = max_mean_diversity_archaea * 1.05
    )

# Plot for Archaea
p_archaea <- ggplot(archaea_data, aes(x = mean, y = y, group = Month)) +
  # Add skyblue rectangles for significant bins FIRST (behind other elements)
  geom_rect(
    data = archaea_segments,
    aes(xmin = segment_x, xmax = segment_xend, ymin = y - 0.5, ymax = y + 0.5),
    fill = "skyblue",
    color = NA,
    inherit.aes = FALSE
  ) +
  # December ribbon
  geom_ribbon(
    data = archaea_data %>% filter(Month == "December"),
    aes(
      y = y,
      xmin = mean - sd,
      xmax = mean + sd
    ),
    fill = "darkgrey",
    alpha = 0.5,
    inherit.aes = FALSE,
    color = NA
  ) +
  # June ribbon
  geom_ribbon(
    data = archaea_data %>% filter(Month == "June"),
    aes(
      y = y,
      xmin = mean - sd,
      xmax = mean + sd
    ),
    fill = "deeppink",
    alpha = 0.3,
    inherit.aes = FALSE,
    color = NA
  ) +
  # December line
  geom_path(
    data = archaea_data %>% filter(Month == "December"),
    aes(color = Month),
    size = 0.5
  ) +
  # June line
  geom_path(
    data = archaea_data %>% filter(Month == "June"),
    aes(color = Month),
    size = 0.5
  ) +
  scale_color_manual(
    values = c("December" = "#000000", "June" = "#FF00FF"),
    name = "Month"
  ) +
  scale_y_continuous(
    breaks = y_labels_to_show,
    labels = as.character(y_labels_to_show)
  ) +
  labs(
    title = "Archaeal richness - December vs June",
    x = "Richness",
    y = "Latitude (°)",
    color = "Month"
  ) +
  theme_minimal(base_size = 8) +
  theme(
    plot.title = element_text(size = 8, face = "bold", hjust = 0.5, vjust = 1),
    axis.text.x = element_text(size = 8, angle = 90, hjust = 1),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.3, "cm"),
    legend.key.height = unit(0.3, "cm"),
    legend.key.width = unit(0.3, "cm"),
    panel.grid = element_blank(),
    legend.position = "none"
  )

# Display plot for Archaea
print(p_archaea)

# Save plot
ggsave(
  paste0(wd_out, "LatGradients_DecemberJune_Archaea_withSegments_v2.svg"), 
  p_archaea, 
  width = 1,
  height = 2.2, 
  units = "in",
  dpi = 300
)


## ====================================================================================================
## END OF SCRIPT
## ====================================================================================================
