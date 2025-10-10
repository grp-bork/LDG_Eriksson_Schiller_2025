# ====================================================================================================
## Monthly latitudinal Richness Gradients for Bacteria and Archaea
##
## Author:       Dominic Eriksson
## Date:         8th of October 2025
## Affiliation:  Environmental Physics Group UP, ETH Zürich, Switzerland
## Contact:      deriksson@ethz.ch
##
## Input files:  
##   - Files in folder: Code/3_Compute_ensembles_and_uncertainties/1_Output/df_long_monthly_ensembleMembers_richness_5000.csv
##
## Output files: 
##   - CSV files containing the 1° binned monthly latitudinal diversity profiles for all clades
##
## Strategy:
##   This script computes latitudinal diversity gradients (LDGs) for microbial clades
##   using monthly ensemble data of alpha-diversity indices (here, richness).
##   1. Load the long-format ensemble data containing monthly richness values for each grid cell.
##   2. Aggregate across models and months to calculate the mean and standard deviation of richness per latitude.
##   3. Compute global ensemble statistics (mean, SD, lower and upper bounds) for each clade.
##   4. Focus on specific months (e.g., June and December) for seasonal comparisons.
##   5. Visualize LDGs
##
## Required R packages (tested versions):
##   - data.table   (version 1.15.4)
##   - dplyr        (version 1.1.2)
##   - ggplot2      (version 3.5.2)
##   - tidyr        (version 1.3.1)
##
## ====================================================================================================

# Clear workspace
rm(list = ls())

# Load libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)

# Directories - adjust as needed
wd_in <- "Code/3_Compute_ensembles_and_uncertainties/1_Output/df_long_monthly_ensembleMembers_richness_5000.csv"
wd_out <- "Code/3_Compute_ensembles_and_uncertainties/5_Output/"

# Create output directory if it doesn't exist
if (!dir.exists(wd_out)) {
  dir.create(wd_out, recursive = TRUE)
  message("Created directory: ", wd_out)
}

# Load data
df_long <- data.table::fread(wd_in)

# Compute ensemble statistics for December and June
ldg_ensemble <- df_long %>%
    dplyr::group_by(y, clade, month) %>%
    dplyr::summarise(
        mean_richness = mean(richness, na.rm = TRUE),
        sd_richness = sd(richness, na.rm = TRUE),
        n_models = n(),
        .groups = "drop"
    ) %>%
    dplyr::mutate(
        lower_bound = mean_richness - sd_richness,
        upper_bound = mean_richness + sd_richness
    )

# Check data
head(ldg_ensemble)

# Name months
# Add month names to ldg_ensemble
ldg_ensemble <- ldg_ensemble %>%
    mutate(
        month = month.name[month]  # Convert numeric month to month name
    )

# Check the result
head(ldg_ensemble)

# Visual check
# Filter for domain_Bacteria, June and December
bacteria_june_dec <- ldg_ensemble %>%
    filter(clade == "class_SAR324", month %in% c("June", "December"))

# Create the plot
p_bacteria_june_dec <- ggplot(bacteria_june_dec, aes(x = y, y = mean_richness, color = month)) +
    # Add ribbons for uncertainty (±sd)
    geom_ribbon(
        aes(ymin = lower_bound, ymax = upper_bound, fill = month),
        alpha = 0.3,
        color = NA
    ) +
    # Add lines for mean richness
    geom_line(size = 1) +
    # Color and fill scales
    scale_color_manual(
        values = c("June" = "#FF00FF", "December" = "#000000"),
        name = "Month"
    ) +
    scale_fill_manual(
        values = c("June" = "#FF00FF", "December" = "#000000"),
        name = "Month"
    ) +
    # Labels and theme
    labs(
        title = "Bacterial Richness LDG: June vs December",
        x = "Latitude (°)",
        y = "Mean Richness",
        color = "Month",
        fill = "Month"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom"
    )

# Display the plot
print(p_bacteria_june_dec)

#* Convert to alternative formatting with stat_type column
# Convert ldg_ensemble to Jonas format
jonas_format_new <- ldg_ensemble %>%
    # Create separate rows for mean and sd
    tidyr::pivot_longer(
        cols = c(mean_richness, sd_richness),
        names_to = "stat_type",
        values_to = "value"
    ) %>%
    # Clean up stat_type names
    dplyr::mutate(
        stat_type = case_when(
            stat_type == "mean_richness" ~ "mean",
            stat_type == "sd_richness" ~ "sd"
        )
    ) %>%
    # Select only needed columns
    dplyr::select(y, clade, month, stat_type, value) %>%
    # Pivot wider to get months as columns
    pivot_wider(
        names_from = month,
        values_from = value
    ) %>%
    # Reorder columns to match Jonas format
    dplyr::select(y, January, February, March, April, May, June, 
           July, August, September, October, November, December, 
           clade, stat_type)

# Check the result
head(jonas_format_new)

# Quick check - plot LDG bacteria June December from jonas_format_new
bacteria_check <- jonas_format_new %>%
    dplyr::filter(grepl("domain_Bacteria", clade)) %>%
    dplyr::select(y, June, December, clade, stat_type) %>%
    tidyr::pivot_longer(
        cols = c(June, December),
        names_to = "Month",
        values_to = "richness"
    ) %>%
    tidyr::pivot_wider(
        names_from = stat_type,
        values_from = richness
    ) %>%
    dplyr::mutate(
        lower_bound = mean - sd,
        upper_bound = mean + sd
    )

# Quick plot
p_check <- ggplot(bacteria_check, aes(x = y, y = mean, color = Month)) +
    geom_ribbon(
        aes(ymin = lower_bound, ymax = upper_bound, fill = Month),
        alpha = 0.3,
        color = NA
    ) +
    geom_line(size = 1) +
    scale_color_manual(
        values = c("June" = "#FF00FF", "December" = "#000000")
    ) +
    scale_fill_manual(
        values = c("June" = "#FF00FF", "December" = "#000000")
    ) +
    labs(
        title = "Quick Check: Bacteria June vs December from Jonas Format",
        x = "Latitude (°)",
        y = "Mean Richness"
    ) +
    theme_minimal()

print(p_check)

# Save the data in Jonas format
data.table::fwrite(
    jonas_format_new,
    paste0(wd_out, "df_latGradients_alphaDiv_month_v5.csv"),
    row.names = FALSE
)

## ====================================================================================================
## END OF SCRIPT
## ====================================================================================================
