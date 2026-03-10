## ====================================================================================================
## We compute annual ensemble means/sd for each diversity index (Richness, Shannon, Chao1).
##
## Author:       Dominic Eriksson
## Date:         8th of October 2025
## Affiliation:  Environmental Physics Group UP, ETH Zürich, Switzerland
## Contact:      deriksson@ethz.ch
##
## Input files:  
##   - Files in folder: /Code/3_Compute_ensembles_and_uncertainties/1_Output/df_long_monthly_ensembleMembers_richness_5000.csv
##
## Output files: 
##   - CSV files containing the long format data frames for each ensemble member.
##   - File names indicate which alpha diversity is saved in the file
##
## Strategy:
##   This script aggregates monthly ensemble predictions of prokaryotic diversity indices 
##   (Richness, Shannon, Chao1) into annual ensemble means and standard deviations. 
##   For each index, spatially resolved predictions are averaged across bootstrap replicates 
##   and algorithms for each clade. The aggregated data are saved as CSVs, and visual checks 
##   of the annual ensemble means are generated using raster plots for Bacteria and Archaea.
##
## Required R packages (tested versions):
##   - data.table  (version 1.15.4)
##   - dplyr       (version 1.1.2)
##   - ggplot2     (version 3.5.2)
##
## ====================================================================================================


# Clear workspace
rm(list = ls())

### =========================================================================
### Richness

# Load libraries
library(data.table)
library(dplyr)
library(ggplot2)

# Directories - adjust as needed
wd_in <- "Code/3_Compute_ensembles_and_uncertainties/1_Output/df_long_monthly_ensembleMembers_richness_5000.csv"
wd_out <- "Code/3_Compute_ensembles_and_uncertainties/2_Output/"
# Create output directory if it does not exist
if (!dir.exists(wd_out)) {
    dir.create(wd_out, recursive = TRUE)
}

# Load data
df_long <- data.table::fread(wd_in)

# Compute annual richness ensemble mean and sd
df_agg <- df_long[
  , .(mean_richness = mean(richness, na.rm = TRUE),
      sd_richness = sd(richness, na.rm = TRUE)),
  by = .(x, y, clade)
]

# Visual check on annual richness d__Bacteria
ggplot(df_agg %>% filter(clade == "domain_Bacteria")) +
  geom_raster(aes(x = x, y = y, fill = mean_richness)) +
  scale_fill_gradient2(
    low = "blue",
    mid = "green",
    high = "red",
    midpoint = median(df_agg$mean_richness, na.rm = TRUE),  # or set a specific midpoint
    na.value = "grey80"
  ) +
  labs(
    title = "Annual Richness Ensemble Mean for d__Bacteria",
    x = "Longitude",
    y = "Latitude",
    fill = "Richness"
  ) +
  theme_minimal()

# # Save this plot
# ggsave(
#   file.path(wd_out, "Global_annual_richness_ensemble_mean_Bacteria.png"),
#   width = 10, height = 6
# )

# Visual check on annual richness d__Archaea
ggplot(df_agg %>% filter(clade == "domain_Archaea")) +
  geom_raster(aes(x = x, y = y, fill = mean_richness)) +
  scale_fill_gradient2(
    low = "blue",
    mid = "green",
    high = "red",
    midpoint = median(df_agg$mean_richness, na.rm = TRUE),  # or set a specific midpoint
    na.value = "grey80"
  ) +
  labs(
    title = "Annual Richness Ensemble Mean for d__Archaea",
    x = "Longitude",
    y = "Latitude",
    fill = "Richness"
  ) +
  theme_minimal()

# # Save this plot
# ggsave(
#   file.path(wd_out, "Global_annual_richness_ensemble_mean_Archaea.png"),
#   width = 10, height = 6
# )

# Save the aggregated data
fwrite(
    df_agg, 
    file.path(wd_out, "Global_annual_richness_ensemble_mean_richness.csv")
)


### =========================================================================
### Shannon

# Load libraries
library(data.table)
library(dplyr)
library(ggplot2)

# Directories
wd_in_shannon <- "Code/3_Compute_ensembles_uncertainties/1_Output/df_long_monthly_ensembleMembers_shannon_5000.csv"

# Load data
df_long_shannon <- data.table::fread(wd_in_shannon)

# Compute annual shannon ensemble mean and sd
df_agg_shannon <- df_long_shannon[
  , .(mean_shannon = mean(shannon, na.rm = TRUE),
      sd_shannon = sd(shannon, na.rm = TRUE)),
  by = .(x, y, clade)
]

# Visual check on annual shannon d__Bacteria
ggplot(df_agg_shannon %>% filter(clade == "domain_Bacteria")) +
  geom_raster(aes(x = x, y = y, fill = mean_shannon)) +
  scale_fill_gradient2(
    low = "blue",
    mid = "green",
    high = "red",
    midpoint = median(df_agg_shannon$mean_shannon, na.rm = TRUE),
    na.value = "grey80"
  ) +
  labs(
    title = "Annual Shannon Ensemble Mean for d__Bacteria",
    x = "Longitude",
    y = "Latitude",
    fill = "Shannon"
  ) +
  theme_minimal()

# # Save this plot
# ggsave(
#   file.path(wd_out, "Global_annual_shannon_ensemble_mean_Bacteria.png"),
#   width = 10, height = 6
# )

# Visual check on annual shannon d__Archaea
ggplot(df_agg_shannon %>% filter(clade == "domain_Archaea")) +
  geom_raster(aes(x = x, y = y, fill = mean_shannon)) +
  scale_fill_gradient2(
    low = "blue",
    mid = "green",
    high = "red",
    midpoint = median(df_agg_shannon$mean_shannon, na.rm = TRUE),
    na.value = "grey80"
  ) +
  labs(
    title = "Annual Shannon Ensemble Mean for d__Archaea",
    x = "Longitude",
    y = "Latitude",
    fill = "Shannon"
  ) +
  theme_minimal()

# # Save this plot
# ggsave(
#   file.path(wd_out, "Global_annual_shannon_ensemble_mean_Archaea.png"),
#   width = 10, height = 6
# )

# Save the aggregated data
fwrite(
    df_agg_shannon, 
    file.path(wd_out, "Global_annual_shannon_ensemble_mean_shannon.csv")
)

### =========================================================================
### Chao1

# Directories
wd_in_chao1 <- "Code/3_Compute_ensembles_uncertainties/1_Output/df_long_monthly_ensembleMembers_chao1_5000.csv"

# Load data
df_long_chao1 <- data.table::fread(wd_in_chao1)

# Compute annual chao1 ensemble mean and sd
df_agg_chao1 <- df_long_chao1[
  , .(mean_chao1 = mean(chao1, na.rm = TRUE),
      sd_chao1 = sd(chao1, na.rm = TRUE)),
  by = .(x, y, clade)
]

# Visual check on annual chao1 d__Bacteria
ggplot(df_agg_chao1 %>% filter(clade == "domain_Bacteria")) +
  geom_raster(aes(x = x, y = y, fill = mean_chao1)) +
  scale_fill_gradient2(
    low = "blue",
    mid = "green",
    high = "red",
    midpoint = median(df_agg_chao1$mean_chao1, na.rm = TRUE),
    na.value = "grey80"
  ) +
  labs(
    title = "Annual Chao1 Ensemble Mean for d__Bacteria",
    x = "Longitude",
    y = "Latitude",
    fill = "Chao1"
  ) +
  theme_minimal()

# # Save this plot
# ggsave(
#   file.path(wd_out, "Global_annual_chao1_ensemble_mean_Bacteria.png"),
#   width = 10, height = 6
# )

# Visual check on annual chao1 d__Archaea
ggplot(df_agg_chao1 %>% filter(clade == "domain_Archaea")) +
  geom_raster(aes(x = x, y = y, fill = mean_chao1)) +
  scale_fill_gradient2(
    low = "blue",
    mid = "green",
    high = "red",
    midpoint = median(df_agg_chao1$mean_chao1, na.rm = TRUE),
    na.value = "grey80"
  ) +
  labs(
    title = "Annual Chao1 Ensemble Mean for d__Archaea",
    x = "Longitude",
    y = "Latitude",
    fill = "Chao1"
  ) +
  theme_minimal()

# # Save this plot
# ggsave(
#   file.path(wd_out, "Global_annual_chao1_ensemble_mean_Archaea.png"),
#   width = 10, height = 6
# )

# Save the aggregated data
fwrite(
    df_agg_chao1, 
    file.path(wd_out, "Global_annual_chao1_ensemble_mean_chao1.csv")
)

## ====================================================================================================
## END OF SCRIPT
## ====================================================================================================
