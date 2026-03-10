## ====================================================================================================
## Annual latitudinal Richness Gradients for Bacteria and Archaea
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
##   - CSV files containing the 1° binned annual latitudinal diversity profiles in two different formats,
##     including one long format style. 
##   - Data frame contains ensemble means and standard deviation across ensemble members and longitude.
##
## Strategy:
##   This script computes latitudinal diversity gradients (LDGs) for microbial clades
##   using the monthly ensemble data of alpha-diversity indices (here, richness).
##   1. Load the long-format ensemble data containing monthly richness values for each grid cell.
##   2. Aggregate across models and months to calculate the mean and standard deviation of richness per latitude.
##   3. Compute global ensemble statistics (mean, SD, lower and upper bounds) for each clade.
##   4. Filter for specific clades if needed for targeted analyses.
##   5. Save results in two formats: 
##       - Standard format with one row per latitude and clade
##       - Alternative “Jonas” format with stat_type column for plotting or downstream analyses.
##
## Required R packages (tested versions):
##   - ggplot2      (version 3.5.2)
##   - dplyr        (version 1.1.2)
##   - sf           (version 1.0.14)
##   - tidyterra    (version 0.4.0)
##   - cowplot      (version 1.1.1)
##   - terra        (version 1.7.78)
##   - data.table   (version 1.15.4)
##   - tidyr        (version 1.3.1)
##
## ====================================================================================================

# Clear workspace
rm(list = ls())

# Load libraries
library(ggplot2)
library(dplyr)
library(sf)
library(tidyterra)
library(cowplot)
library(terra)
library(data.table)
library(tidyr)

# Directories - adjust as needed
wd_in <- "Code/3_Compute_ensembles_and_uncertainties/1_Output/df_long_monthly_ensembleMembers_richness_5000.csv"
wd_out <- "Code/3_Compute_ensembles_and_uncertainties/4_Output/"

# Create output directory if it doesn't exist
if (!dir.exists(wd_out)) {
  dir.create(wd_out, recursive = TRUE)
  message("Created directory: ", wd_out)
}

# Load data
df_long <- data.table::fread(wd_in)

## ====================================================================================================
## Bacteria LDG
## ====================================================================================================

# Compute ensemble statistics for Global (all regions)
bacteria_ldg_global <- df_long %>%
    dplyr::group_by(y, clade) %>%
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

# Quick check of data
head(bacteria_ldg_global)

# Filter for specific clade
bacteria_only <- bacteria_ldg_global %>%
    filter(clade == "class_SAR324")

# Save
write.csv(bacteria_ldg_global,
          file = paste0(wd_out, "df_latGradients_alphaDiv_annual_v6.csv"),
          row.names = FALSE)

## ====================================================================================================
## Convert to Alternative Format
## ====================================================================================================

# Alternative formatting with stat_type column
jonas_format_new <- bacteria_ldg_global %>%
    tidyr::pivot_longer(cols = c(mean_richness, sd_richness),
                 names_to = "stat_type", values_to = "value") %>%
    dplyr::mutate(stat_type = case_when(
        stat_type == "mean_richness" ~ "mean",
        stat_type == "sd_richness" ~ "sd"
    )) %>%
    dplyr::select(y, clade, stat_type, value) %>%
    tidyr::pivot_wider(names_from = clade, values_from = value) %>%
    dplyr::select(y, stat_type, everything()) %>%
    dplyr::arrange(y, stat_type)

# Check
head(jonas_format_new)

# Save Jonas-format CSV
write.csv(jonas_format_new,
          file = paste0(wd_out, "df_latGradients_alphaDiv_annual__jonas_v8.csv"),
          row.names = FALSE)

## ====================================================================================================
## END OF SCRIPT
## ====================================================================================================
