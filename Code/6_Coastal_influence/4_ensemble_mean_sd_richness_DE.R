## ====================================================================================================
## This R script computes ensemble statistics (means, standard deviations, and coefficients of 
## variation) from prokaryotic richness predictions across coastal distance groups. It quantifies
## model uncertainty by analyzing variability between different distance thresholds and creates
## spatial visualizations of prediction confidence for biodiversity assessments.
##
## Author:       Dominic Eriksson
## Date:         26th of February 2026
## Affiliation:  Environmental Physics Group UP, ETH Zürich, Switzerland
## Contact:      deriksson@ethz.ch
##
## Input files:  
##   - Long-format CSV with ensemble members from all coastal distance groups
##   - Contains bootstrap replicates across spatial and temporal dimensions
##
## Output files: 
##   - CSV file with ensemble statistics and uncertainty metrics
##   - PNG visualizations of standard deviation and coefficient of variation maps
##   - Spatial uncertainty assessments for model validation
##
## Strategy:
##   The script aggregates ensemble predictions across multiple levels: first computing annual
##   means from monthly predictions, then averaging across bootstrap replicates and algorithms,
##   and finally calculating uncertainty metrics across distance groups. This hierarchical
##   approach quantifies different sources of model uncertainty while preserving spatial
##   structure for robust biodiversity pattern analysis.
##
## Required R packages (tested versions):
##   - data.table 1.14.8
##   - dplyr      1.1.2
##   - ggplot2    3.4.2
## ====================================================================================================

# Clear workspace
rm(list = ls())

# Load libraries
library(data.table) # For efficient data handling
library(dplyr)      # For data manipulation
library(ggplot2)    # For visualization

# Directories
wd_in <- "Code/6_Coastal_influence/3_Output/df_long_monthly_ensembleMembers_richness_5000_v2.csv"
wd_out <- "Code/6_Coastal_influence/4_Output/"
# Create output directory if it does not exist
if (!dir.exists(wd_out)) {
    dir.create(wd_out, recursive = TRUE)
}

# Load data
df_long <- data.table::fread(wd_in)

# Compute annual richness for each ensemble member
df_agg <- df_long[
  , .(annual_richness = mean(richness, na.rm = TRUE)),
  by = .(x, y, bootstrap, alg, clade)
]

# Compute annual means for each distance group
df_agg <- df_agg[
  , .(annual_richness = mean(annual_richness, na.rm = TRUE)),  
  by = .(x, y, clade)
]

# Compute model uncertainty as standard deviation and coefficient of variation across distance groups/clades
df_agg <- df_agg[
  , .(
      mean_richness = mean(annual_richness, na.rm = TRUE),
      sd_richness = sd(annual_richness, na.rm = TRUE),
      cv_richness = (sd(annual_richness, na.rm = TRUE) / mean(annual_richness, na.rm = TRUE)*100)
    ),
  by = .(x, y)
]

# Visual check standard deviation - using appropriate sequential color scheme
ggplot(df_agg) +
  geom_raster(aes(x = x, y = y, fill = sd_richness)) +
  scale_fill_viridis_c(
    option = "plasma",  # Good for uncertainty: purple/pink to yellow
    na.value = "grey80",
    name = "Standard\nDeviation"
  ) +
  labs(
    title = "Model uncertainty (SD) across coastal distance groups for prokaryotes",
    x = "Longitude", 
    y = "Latitude",
    fill = "SD Richness"
  ) +
  theme_minimal()

# Save this plot
ggsave(
  file.path(wd_out, "Model_standard_deviation_across_distance_groups_annualRichness_prokaryotes.png"),
  width = 10, height = 6
) 

# Visual check on coefficient of variation
ggplot(df_agg) +
  geom_raster(aes(x = x, y = y, fill = cv_richness)) +
  scale_fill_viridis_c(
    option = "inferno",  # Good for CV: dark purple to yellow/white
    na.value = "grey80",
    name = "CV\n(%)"
  ) +
  labs(
    title = "Model uncertainty (CV) across coastal distance groups for prokaryotes",
    x = "Longitude",
    y = "Latitude",
    fill = "CV Richness (%)"
  ) +
  theme_minimal()

# Save this plot
ggsave(
  file.path(wd_out, "Model_coefficient_of_variation_across_distance_groups_annualRichness_prokaryotes.png"),
  width = 10, height = 6
)

# Save output data
fwrite(df_agg, file.path(wd_out, "Model_uncertainty_across_distance_groups_annualRichness_prokaryotes.csv"))

## ====================================================================================================
## END OF SCRIPT
## ====================================================================================================