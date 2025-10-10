## ====================================================================================================
## This script uses Wilcoxon tests to compare richness between June and December for each 1° latitudinal bin.
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
##   - CSV files containing the richness hotspot information and the environmental conditions at these locations.
##
## Strategy:
##   This script evaluates latitudinal gradients of microbial richness between June and December.
##   1. Load the monthly ensemble richness data for all grid cells and clades.
##   2. Aggregate latitudes into 1° bins and filter for June (month 6) and December (month 12).
##   3. Generate a unique model ID per bootstrap and algorithm combination.
##   4. Perform a Wilcoxon rank-sum test for each latitude bin, clade, and model to compare June vs December richness.
##   5. Determine the fraction of significant bins per model and clade and summarize across models for each latitude.
##   6. Calculate mean richness per clade and month for visualization.
##   7. Visualize the latitudinal gradients with colored segments indicating bins where >50% of models show significant seasonal differences.
##   8. Save all results including:
##       - Wilcoxon test results by model and bin
##       - Fraction of significant bins per model
##       - Summary of bins with model significance
##       - Bins exceeding the 50% significance threshold for downstream analyses.

##
## Required R packages (tested versions):
##   - dplyr        (version 1.1.2)
##   - data.table   (version 1.15.4)
##   - ggplot2      (version 3.5.2)
##   - tidyr        (version 1.3.1)
##   - ggnewscale   (version 0.4.8)
##
## ====================================================================================================

# Load libraries
library(dplyr)
library(data.table)
library(ggplot2)
library(tidyr)
library(ggnewscale)

# Directories - adjust as needed
wd_out <- "Code/5_Wilcox_test/1_Output/"

# Create output directory if it doesn't exist
if (!dir.exists(wd_out)) {
  dir.create(wd_out, recursive = TRUE)
  message("Created directory: ", wd_out)
}

# Load data
df_ensemble <- data.table::fread("Code/3_Compute_ensembles_and_uncertainties/1_Output/df_long_monthly_ensembleMembers_richness_5000.csv")

# Save original loaded data
original <- df_ensemble

# Bin latitude into 1° bins if not already done
df_ensemble <- df_ensemble %>%
  mutate(y = floor(y))  # y is latitude bin

# Filter for June and December
df_test <- df_ensemble %>%
  filter(month %in% c(6, 12))

# Create model ID
df_test <- df_test %>%
  mutate(model_id = paste(bootstrap, alg, sep = "_"))

# Wilcoxon test for each clade, latitude bin, and model
wilcox_results_by_model <- df_test %>%
  group_by(y, clade, model_id, bootstrap, alg) %>%
  summarise(
    wilcox_pairwise_p_value = tryCatch({
      dec_vals <- richness[month == 12]
      jun_vals <- richness[month == 6]
      if(length(dec_vals) > 0 & length(jun_vals) > 0 & 
         !all(is.na(dec_vals)) & !all(is.na(jun_vals))) {
        wilcox.test(jun_vals, dec_vals, exact = FALSE)$p.value
      } else {
        NA_real_
      }
    }, error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  mutate(significant = !is.na(wilcox_pairwise_p_value) & wilcox_pairwise_p_value < 0.001)

# Calculate fraction of significant bins per model and clade
fraction_per_model <- wilcox_results_by_model %>%
  group_by(clade, model_id, bootstrap, alg) %>%
  summarise(
    total_bins = sum(!is.na(wilcox_pairwise_p_value)),
    significant_bins = sum(significant, na.rm = TRUE),
    fraction_significant = ifelse(total_bins > 0, significant_bins / total_bins, 0),
    .groups = "drop"
  )

# Calculate mean fraction across models for each latitude bin and clade
bin_significance_summary <- wilcox_results_by_model %>%
  group_by(y, clade) %>%
  summarise(
    models_with_data = sum(!is.na(wilcox_pairwise_p_value)),
    models_significant = sum(significant, na.rm = TRUE),
    fraction_models_significant = ifelse(models_with_data > 0, models_significant / models_with_data, 0),
    .groups = "drop"
  ) %>%
  # Only keep bins where at least some models show significance
  filter(fraction_models_significant > 0)

# Calculate mean richness for plotting
mean_richness <- df_test %>%
  group_by(clade, y, month) %>%
  summarise(mean_richness = mean(richness, na.rm = TRUE), .groups = "drop")

# Create plot data
plot_df <- mean_richness %>%
  filter(clade %in% c("domain_Bacteria", "domain_Archaea"))

# Get max richness for segment positioning
max_richness_per_clade <- plot_df %>%
  group_by(clade) %>%
  summarise(max_richness = max(mean_richness, na.rm = TRUE), .groups = "drop")


# Create segment data for ALL clades (not just bacteria and archaea)
segment_data_all <- wilcox_results_by_model %>%
  dplyr::group_by(y, clade) %>%
  dplyr::summarise(
    models_with_data = sum(!is.na(wilcox_pairwise_p_value)),
    models_significant = sum(significant, na.rm = TRUE),
    fraction_models_significant = ifelse(models_with_data > 0, models_significant / models_with_data, 0),
    .groups = "drop"
  ) %>%
  # Remove the filter to include all clades
  # filter(clade %in% c("domain_Bacteria", "domain_Archaea")) %>%
  left_join(max_richness_per_clade, by = "clade") %>%
  mutate(
    y_start = max_richness * 1.02,
    y_end = max_richness * 1.08,
    # Create color categories
    color_category = ifelse(fraction_models_significant > 0.50, "above_50", "below_50")
  )

# Separate data for different coloring
segment_data_above_50 <- segment_data_all %>% filter(color_category == "above_50")
segment_data_below_50 <- segment_data_all %>% filter(color_category == "below_50" & fraction_models_significant > 0)

# Create the plot with >50% threshold
p_50 <- ggplot(plot_df, aes(x = y, y = mean_richness)) +
  # Add the line plots for June and December
  geom_line(
    data = plot_df %>% filter(month == 6),
    aes(group = 1), 
    color = "#1f77b4", 
    size = 0.8
  ) +
  geom_line(
    data = plot_df %>% filter(month == 12),
    aes(group = 1), 
    color = "#ff7f0e", 
    size = 0.8
  ) +
  # Add WHITE segments for bins with <50% significance
  geom_segment(
    data = segment_data_below_50,
    aes(x = y, xend = y, y = y_start, yend = y_end),
    color = "white",
    size = 2
  ) +
  # Add SKYBLUE segments for bins with >50% significance
  geom_segment(
    data = segment_data_above_50,
    aes(x = y, xend = y, y = y_start, yend = y_end),
    color = "skyblue",
    size = 2
  ) +
  facet_wrap(~ clade, scales = "free_y") +
  labs(
    title = "Latitudinal Gradient of Richness (June vs December)\nwith Model Significance Fractions (>75% threshold)",
    x = "Latitude bin (y)",
    y = "Mean richness"
  ) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "right") +
  # Add manual legend for the months
  annotate("text", x = Inf, y = Inf, label = "June", color = "#1f77b4", 
           hjust = 1.1, vjust = 2, size = 3) +
  annotate("text", x = Inf, y = Inf, label = "December", color = "#ff7f0e", 
           hjust = 1.1, vjust = 3.5, size = 3)

# Display the plot
print(p_50)


# Save the dataframes with wilcox results and fractions
data.table::fwrite(
  wilcox_results_by_model, 
  paste0(wd_out, "wilcox_results_by_model_and_bin.csv"), 
  row.names = FALSE
)

data.table::fwrite(
  fraction_per_model, 
  paste0(wd_out, "fraction_significant_bins_per_model.csv"), 
  row.names = FALSE
)

data.table::fwrite(
  segment_data_all, 
  paste0(wd_out, "bin_significance_summary_all_bins.csv"), 
  row.names = FALSE
)

# Also save just the bins above 50% threshold
data.table::fwrite(
  segment_data_above_50, 
  paste0(wd_out, "bins_above_50pct_significant.csv"), 
  row.names = FALSE
)

## ====================================================================================================
## END OF SCRIPT
## ====================================================================================================
