## ====================================================================================================
## R script to compute the difference between December and June
##
## Author:       Dominic Eriksson
## Date:         8th of October 2025
## Affiliation:  Environmental Physics Group UP, ETH Zürich, Switzerland
## Contact:      deriksson@ethz.ch
##
## Input files:  
##   - Files in folder: /Code/3_Compute_ensembles_and_uncertainties/3_Output/
##
## Output files: 
##   - CSV files containing the difference between December and June of the corresponding diversity index.
##
## Strategy:
##   This script computes the difference between December and June ensemble predictions 
##   for prokaryotic diversity indices (Richness, Shannon, Chao1) across spatial grids. 
##   Monthly ensemble predictions are filtered for June and December, aggregated by clade, 
##   algorithm, and bootstrap, and reshaped to calculate differences and directional changes. 
##   Summary statistics—including mean, standard deviation, coefficient of variation, and 
##   agreement among directional change computed for each grid cell. The results are saved as CSVs, 
##   and raster plots visualize the spatial patterns of December–June differences for Bacteria 
##   and Archaea.
##
## Required R packages (tested versions):
## Required R packages (tested versions):
##   - data.table  (version 1.15.4)
##   - dplyr       (version 1.1.2)
##   - ggplot2     (version 3.5.2)
##
## ====================================================================================================

# Clear work space
rm(list = ls())

# Load libraries
library(data.table)
library(dplyr)
library(ggplot2)

# Base project directory - adjust as needed
project_dir <- "LDG_github_script_and_data/" # Adjust to your project directory 

# Input/Output directories
wd_out <- paste0(project_dir, "/Code/3_Compute_ensembles_and_uncertainties/3_Output/")
if(!dir.exists(wd_out)) dir.create(wd_out, recursive = TRUE)

# Function to compute difference and summary
compute_dec_june_diff <- function(file_path, value_col, clade_label) {
  df_long <- fread(file_path)
  
  # Filter June (6) and December (12)
  df_sub <- df_long[month %in% c(6, 12)]
  
  # Step 1: Aggregate by grouping variables
  df_agg <- df_sub[, .(val = first(get(value_col))), 
                   by = .(x, y, clade, alg, bootstrap, month)]
  
  # Step 2: Pivot to wide format
  df_wide <- dcast(df_agg, x + y + clade + alg + bootstrap ~ month, value.var = "val")
  setDT(df_wide)
  
  # Step 3: Calculate difference and create final dataset
  diff_col_name <- paste0("dec_june_", value_col, "Diff")
  june_col_name <- paste0(value_col, "_june")
  dec_col_name <- paste0(value_col, "_dec")
  #
  df_diff <- df_wide[, c(diff_col_name) := `12` - `6`
                     ][, c("x", "y", "clade", "alg", "bootstrap", 
                           june_col_name, dec_col_name, diff_col_name) := 
                         .(x, y, clade, alg, bootstrap, `6`, `12`, get(diff_col_name))]
  
  # Keep only the columns we want
  keep_cols <- c("x", "y", "clade", "alg", "bootstrap", june_col_name, dec_col_name, diff_col_name)
  df_diff <- df_diff[, ..keep_cols]
  
  # Add direction
  df_diff <- df_diff %>%
    mutate(direction = ifelse(get(dec_col_name) > get(june_col_name), 1,
                              ifelse(get(dec_col_name) < get(june_col_name), 0, NA)))
  
  # Compute ensembles and other statistics
  df_summary <- df_diff %>%
    group_by(x, y, clade) %>%
    summarise(
      mean_diff = mean(get(diff_col_name), na.rm = TRUE),
      sd_diff = sd(get(diff_col_name), na.rm = TRUE),
      cv_diff = ifelse(mean_diff == 0, NA, sd_diff / mean_diff * 100),
      n_models = sum(!is.na(direction)),
      n_increase = sum(direction == 1, na.rm = TRUE),
      n_decrease = sum(direction == 0, na.rm = TRUE),
      agreement_diff = 100 * pmax(n_increase, n_decrease) / n_models,
      .groups = "drop"
    )
  
  list(diff = df_diff, summary = df_summary)
}

# Function to plot raster difference
plot_diff <- function(df_summary, value_col, clade_name, output_file, title_text) {
  ggplot(df_summary %>% filter(clade == clade_name), 
         aes(x = x, y = y, fill = mean_diff)) +
    geom_raster() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, na.value = "grey80") +
    coord_fixed() +
    labs(x = "Longitude", y = "Latitude", fill = "Mean Difference", title = title_text) +
    theme_minimal() -> p
  ggsave(output_file, plot = p, width = 10, height = 6)
}

# =========================
# Richness
# =========================
# Apply function
res <- compute_dec_june_diff(
  file_path = paste0(project_dir, "/Code/3_Compute_ensembles_and_uncertainties/1_Output/df_long_monthly_ensembleMembers_richness_5000.csv"),
  value_col = "richness",
  clade_label = "Richness"
)

# Check res$summary
summary(res$summary)

# Save
fwrite(res$summary, paste0(wd_out, "Global_map_difference_December_June_ensemble__mean_sd_richness.csv"))

# Check some quick plots
plot_diff(res$summary, "richness", "domain_Bacteria", 
          file.path(wd_out, "Global_richness_difference_December_June_Bacteria.png"), 
          "Richness Difference (Dec - Jun): d__Bacteria")
plot_diff(res$summary, "richness", "domain_Archaea", 
          file.path(wd_out, "Global_richness_difference_December_June_Archaea.png"), 
          "Richness Difference (Dec - Jun): d__Archaea")

# =========================
# Shannon
# =========================
# Apply function
res_sh <- compute_dec_june_diff(
  file_path = paste0(project_dir, "/Code/3_Compute_ensembles_and_uncertainties/1_Output/df_long_monthly_ensembleMembers_shannon_5000.csv"),
  value_col = "shannon",
  clade_label = "Shannon"
)

# Save file
fwrite(res_sh$summary, paste0(wd_out, "Global_map_difference_December_June_ensemble_mean_sd_shannon.csv"))

plot_diff(res_sh$summary, "shannon", "domain_Bacteria", 
          file.path(wd_out, "Global_shannon_difference_December_June_Bacteria.png"), 
          "Shannon Difference (Dec - Jun): d__Bacteria")
plot_diff(res_sh$summary, "shannon", "domain_Archaea", 
          file.path(wd_out, "Global_shannon_difference_December_June_Archaea.png"), 
          "Shannon Difference (Dec - Jun): d__Archaea")

# =========================
# Chao1
# =========================
# Apply function
res_chao1 <- compute_dec_june_diff(
  file_path = paste0(project_dir, "/Code/3_Compute_ensembles_and_uncertainties/1_Output/df_long_monthly_ensembleMembers_chao1_5000.csv"),
  value_col = "chao1",
  clade_label = "Chao1"
)

# Save file
fwrite(res_chao1$summary, file.path(wd_out, "Global_map_difference_December_June_ensemble_mean_sd_chao1.csv"))

plot_diff(res_chao1$summary, "chao1", "domain_Bacteria", 
          file.path(wd_out, "Global_chao1_difference_December_June_Bacteria.png"), 
          "Chao1 Difference (Dec - Jun): d__Bacteria")
plot_diff(res_chao1$summary, "chao1", "domain_Archaea", 
          file.path(wd_out, "Global_chao1_difference_December_June_Archaea.png"), 
          "Chao1 Difference (Dec - Jun): d__Archaea")

## ====================================================================================================
## END OF SCRIPT
## ====================================================================================================
