## ====================================================================================================
## This script generates global latitudinal diversity gradients (LDG) for modeled bacterial and archaeal annual richness.
## It visualizes the mean richness and standard deviation along latitudes and across ensemble members.
##
## Author:       Dominic Eriksson
## Date:         29th of September, 2025
## Affiliation:  Environmental Physics Group, UP, ETH Zürich, Switzerland
## Contact:      deriksson@ethz.ch
##
## Input files:  
##   - Annual alpha-diversity CSV: Code/3_Compute_ensembles_and_uncertainties/4_Output/df_latGradients_alphaDiv_annual__jonas_v8.csv
##
## Output files: 
##   - Lineplots of latitudinal richness for bacterial and archaeal clades
##
## Strategy:
##   This script visualizes the latitudinal distribution of microbial richness while incorporating variability:
##   1. Load annual alpha-diversity data for bacterial and archaeal clades.
##   2. Pivot the data to have separate columns for mean and standard deviation.
##   3. Compute upper and lower bounds (mean ± SD) for each latitude.
##   4. Identify the maximum mean richness and annotate it on the plot.
##   5. Generate ribbon plots for SD and overlay mean richness lines.
##   6. Add global labels and annotations for clarity.
##   7. Save high-resolution SVG plots for each clade.
##
## Required R packages (tested versions):
##   - data.table   (version 1.15.4)
##   - dplyr        (version 1.1.2)
##   - tidyr        (version 1.3.1)
##   - ggplot2      (version 3.5.2)
##
## ====================================================================================================


# Load necessary libraries
# Required libraries
library(data.table)  # for fread()
library(dplyr)       # for %>%, select(), filter(), mutate()
library(tidyr)       # for pivot_wider()
library(ggplot2)     # for ggplot(), geom_ribbon(), geom_path(), etc.

# Directory
wd_out <- "Code/Figures/Fig3/"

# Create output directory if it doesn't exist
if (!dir.exists(wd_out)) {
  dir.create(wd_out, recursive = TRUE)
}

# Load data
df <- fread("Code/3_Compute_ensembles_and_uncertainties/4_Output/df_latGradients_alphaDiv_annual__jonas_v8.csv")

## LDG for Bacteria --------------------------------------------------------

# Extract domain_Bacteria data from the dataframe
bacteria_data <- df %>%
    dplyr::select(y, stat_type, domain_Bacteria) %>%
    dplyr::filter(!is.na(domain_Bacteria)) %>%
    tidyr::pivot_wider(
        names_from = stat_type,
        values_from = domain_Bacteria
    ) %>%
    filter(!is.na(mean) & !is.na(sd)) %>%
    mutate(
        lower_bound = mean - sd,
        upper_bound = mean + sd
    )

# Find the maximum mean diversity
max_mean_diversity <- round(max(bacteria_data$mean, na.rm = TRUE), 0)

# Define y-axis labels
y_labels_to_show <- c(-75, -50, -25, 0, 25, 50, 75)

# Create the plot
p <- ggplot() +
    # Global richness ribbon and line
    geom_ribbon(
        data = bacteria_data, 
        aes(y = y, xmin = lower_bound, xmax = upper_bound),
        fill = "darkgrey", alpha = 0.7
    ) +
    geom_path(
        data = bacteria_data, 
        aes(x = mean, y = y, group = 1), 
        color = "black", 
        size = 0.5
    ) +
    # Vertical line at maximum richness
    geom_vline(
        xintercept = max_mean_diversity, 
        linetype = "dashed", 
        color = "red",
        size = 0.2
    ) +
    # Annotation for maximum richness
    annotate(
        "text", 
        x = Inf, 
        y = Inf, 
        label = paste("Max. mean =", round(max_mean_diversity, 2)), 
        color = "red", 
        vjust = 1.5, 
        hjust = 1.1, 
        size = 8 / .pt
    ) +
    # Global label
    annotate(
        "text", 
        x = Inf, 
        y = -Inf, 
        label = "Global", 
        color = "black", 
        vjust = -1.5, 
        hjust = 1.1, 
        size = 8 / .pt, 
        fontface = "bold"
    ) +
    # Axis and theme settings
    scale_y_continuous(breaks = y_labels_to_show) +
    theme_minimal(base_size = 8) +
    theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
        legend.position = "none",
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        panel.grid = element_blank()
    ) +
    labs(
        title = "Modeled bacterial richness",
        x = "Richness",
        y = "Latitude (°)"
    )

# Display plot
print(p)

# Save plot
ggsave(
  paste0(wd_out, "LatGradients_annualRichness_Bacteria.svg"), 
  p, 
  width = 1,
  height = 2.2, 
  units = "in",
  dpi = 300
)

## LDG for Archaea --------------------------------------------------------

# Extract domain_Archaea data from the dataframe
archaea_data <- df %>%
    dplyr::select(y, stat_type, domain_Archaea) %>%
    filter(!is.na(domain_Archaea)) %>%
   tidyr::pivot_wider(
        names_from = stat_type,
        values_from = domain_Archaea
    ) %>%
    filter(!is.na(mean) & !is.na(sd)) %>%
    mutate(
        lower_bound = mean - sd,
        upper_bound = mean + sd
    )

# Find the maximum mean diversity
max_mean_diversity <- round(max(archaea_data$mean, na.rm = TRUE), 0)

# Define y-axis labels
y_labels_to_show <- c(-75, -50, -25, 0, 25, 50, 75)

# Create the plot
p_archaea <- ggplot() +
    # Global richness ribbon and line
    geom_ribbon(
        data = archaea_data, 
        aes(y = y, xmin = lower_bound, xmax = upper_bound),
        fill = "darkgrey", alpha = 0.7
    ) +
    geom_path(
        data = archaea_data, 
        aes(x = mean, y = y, group = 1), 
        color = "black", 
        size = 0.5
    ) +
    # Vertical line at maximum richness
    geom_vline(
        xintercept = max_mean_diversity, 
        linetype = "dashed", 
        color = "red",
        size = 0.2
    ) +
    # Annotation for maximum richness
    annotate(
        "text", 
        x = Inf, 
        y = Inf, 
        label = paste("Max. mean =", round(max_mean_diversity, 2)), 
        color = "red", 
        vjust = 1.5, 
        hjust = 1.1, 
        size = 8 / .pt
    ) +
    # Global label
    annotate(
        "text", 
        x = Inf, 
        y = -Inf, 
        label = "Global", 
        color = "black", 
        vjust = -1.5, 
        hjust = 1.1, 
        size = 8 / .pt, 
        fontface = "bold"
    ) +
    # Axis and theme settings
    scale_y_continuous(breaks = y_labels_to_show) +
    theme_minimal(base_size = 8) +
    theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
        legend.position = "none",
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        panel.grid = element_blank()
    ) +
    labs(
        title = "Modeled archaeal richness",
        x = "Richness",
        y = "Latitude (°)"
    )

# Display plot
print(p_archaea)

# Save plot
ggsave(
  paste0(wd_out, "LatGradients_annualRichness_Archaea_V8.svg"), 
  p_archaea, 
  width = 1,
  height = 2.2, 
  units = "in",
  dpi = 300
)

## ====================================================================================================
## END OF SCRIPT
## ====================================================================================================
