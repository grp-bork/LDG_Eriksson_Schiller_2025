## ====================================================================================================
## This script generates bubble plots showing the top 3 environmental predictors per microbial clade.
## It visualizes the median importance along with asymmetric error bars (IQR) for each predictor, grouped by taxonomic clades.
##
## Author:       Dominic Eriksson
## Date:         8th of October, 2025
## Affiliation:  Environmental Physics Group, UP, ETH Zürich, Switzerland
## Contact:      deriksson@ethz.ch
##
## Input files:  
##   - Final selected model files: 
##       Code/3_Compute_ensembles_and_uncertainties/1_Output/final_selected_files_Richness_5000.txt
##   - Predictor ranking CSVs: 
##       Code/2_Extract_model_outputs/1_Output/Non_log/Predictor_rankings_MLD_v5.csv
##       Code/2_Extract_model_outputs/1_Output/Log/Predictor_rankings_MLD_v5.csv
##
## Output files: 
##   - Bubble plots with asymmetric error bars for top 3 predictors per taxon:
##
## Strategy:
##   This script visualizes the top predictors for each taxon while accounting for variability:
##   1. Load the final selected model files for log and non-log clades.
##   2. Load predictor ranking CSVs and subset to relevant clades.
##   3. Merge log and non-log clades, filtering only target clades and removing unwanted ones.
##   4. Assign groups to clades for faceted plotting.
##   5. Recode predictor names to human-readable labels.
##   6. Compute median and IQR (Q25-Q75) for each predictor per taxon.
##   7. Identify top 3 predictors per taxon and calculate asymmetric distances for error bars.
##   8. Define color mapping for environmental predictor categories.
##   9. Generate bubble plots with asymmetric error bars, overlay median values, and annotate IQR.
##  10. Save high-resolution SVG plots for publication-quality figures.
##
## ====================================================================================================================


### ==============================================================================================
### Predictor ranking boxplot for clades
### ==============================================================================================

# Clear workspace
rm(list = ls())

# Load libraries
library(stringr)
library(ggplot2)

# Directories
wd_out <- "Code/Figures/Fig5/"

# Create output directory if it doesn't exist
if (!dir.exists(wd_out)) {
  dir.create(wd_out, recursive = TRUE)
}

### ==============================================================================================
### Bubble plots for top 3 predictors per taxon
### ==============================================================================================

# Load the final selected files from the text file
fnames_model_stats <- readLines("Code/3_Compute_ensembles_and_uncertainties/1_Output/final_selected_files_Richness_5000.txt")

# Define which clades are log and non-log based on the file paths
log_clades <- fnames_model_stats[grepl("_Log_v2\\.RData$", fnames_model_stats)]
nonlog_clades <- fnames_model_stats[grepl("_Non_log_v2\\.RData$", fnames_model_stats)]

# Extract clade names from file paths (remove path and extension)
extract_clade <- function(x) {
  # Extract the clade name from the file path
  basename_file <- basename(x)
  # Remove the suffix after the clade name
  clade_name <- sub("_Richness_5000_ensembleMembers_(Log|Non_log)_v2\\.RData$", "_Richness_5000", basename_file)
  return(clade_name)
}

# Check if log files were used
log_clade_names <- extract_clade(log_clades)
nonlog_clade_names <- extract_clade(nonlog_clades)

# Print to check what we extracted
cat("Log clades:\n")
print(log_clade_names)
cat("\nNon-log clades:\n")
print(nonlog_clade_names)

# Load data
df_nonLog <- read.csv("Code/2_Extract_model_outputs/1_Output/Non_log/Predictor_rankings_MLD_v5.csv")
df_Log <- read.csv("Code/2_Extract_model_outputs/1_Output/Log/Predictor_rankings_MLD_v5.csv")

# Subset each dataframe for the relevant clades
df_nonLog_sub <- df_nonLog[df_nonLog$clade %in% nonlog_clade_names, ]
df_Log_sub <- df_Log[df_Log$clade %in% log_clade_names, ]

# Check how many rows we got from each
cat("Rows in df_nonLog_sub:", nrow(df_nonLog_sub), "\n")
cat("Rows in df_Log_sub:", nrow(df_Log_sub), "\n")

# Merge the two dataframes
df_all <- dplyr::bind_rows(df_nonLog_sub, df_Log_sub)

# Filter for clades that start with "class_" or "domain_"
df_all <- df_all %>%
  filter(grepl("^(class_|domain_)", clade))

# Remove clade class_Calditrichia_Richness_5000
df_all <- df_all %>%
  filter(clade != "class_Calditrichia_Richness_5000")

cat("Total rows in merged df_all:", nrow(df_all), "\n")
cat("Unique clades in df_all:", length(unique(df_all$clade)), "\n")

# Create group assignments based on clade names (updated with domains)
group_assignments <- data.frame(
  clade = c(
    "domain_Bacteria_Richness_5000",
    "domain_Archaea_Richness_5000",
    "class_Poseidoniia_Richness_5000",
    "class_Nitrososphaeria_Richness_5000",
    "class_Acidimicrobiia_Richness_5000",
    "class_Bacteroidia_Richness_5000",
    "class_Chlamydiia_Richness_5000",
    "class_Cyanobacteriia_Richness_5000",
    "class_UBA1144_Richness_5000",
    "class_Fibrobacteria_Richness_5000",
    "class_Marinisomatia_Richness_5000",
    "class_UBA4151_Richness_5000",
    "class_UBA796_Richness_5000",
    "class_XYA12-FULL-58-9_Richness_5000",
    "class_Alphaproteobacteria_Richness_5000",
    "class_Gammaproteobacteria_Richness_5000",
    "class_SAR324_Richness_5000"
  ),
  group = c( # Grouping found in main figure 4
    "group2",  # domain_Bacteria (added to Group 2)
    "group4",  # domain_Archaea (added to Group 4)
    "group4",  # Poseidoniia (moved to Group 4)
    "group1",  # Nitrososphaeria
    "group2",  # Acidimicrobiia (moved to Group 2)
    "group4",  # Bacteroidia (moved to Group 4)
    "group3",  # Chlamydiia (moved to Group 3)
    "group2",  # Cyanobacteriia (moved to Group 2)
    "group4",  # UBA1144 (moved to Group 4)
    "group1",  # Fibrobacteria
    "group4",  # Marinisomatia (moved to Group 4)
    "group3",  # UBA4151 (moved to Group 3)
    "group3",  # UBA796 (moved to Group 3)
    "group3",  # XYA12-FULL-58-9 (moved to Group 3)
    "group2",  # Alphaproteobacteria (moved to Group 2)
    "group4",  # Gammaproteobacteria (moved to Group 4)
    "group4"   # SAR324 (moved to Group 4)
  ),
  simple_name = c(
    "Bacteria", "Archaea", "Poseidoniia", "Nitrososphaeria", "Acidimicrobiia", "Bacteroidia",
    "Chlamydiia", "Cyanobacteriia", "UBA1144", "Fibrobacteria",
    "Marinisomatia", "UBA4151", "UBA796", "XYA12-FULL-58-9",
    "Alphaproteobacteria", "Gammaproteobacteria", "SAR324"
  )
)

# Updated Group labels (now 4 groups instead of 5)
group_labels <- c(
  "group1" = "Group 1: Nitrososphaeria & Fibrobacteria",
  "group2" = "Group 2: Bacteria, Acidimicrobiia, Cyanobacteriia, Alphaproteobacteria",
  "group3" = "Group 3: Chlamydiia, UBA4151, UBA796, XYA12-FULL-58-9",
  "group4" = "Group 4: Archaea, Poseidoniia, Bacteroidia, UBA1144, Marinisomatia, Gammaproteobacteria, SAR324"
)

# Add group information to the data
df_all <- df_all %>%
  left_join(group_assignments, by = "clade") %>%
  filter(!is.na(group))  # Remove any clades not in our groups

# Adjust predictor names
column_name_map <- c(
  "climatology_A_0_50" = "AOU",
  "climatology_A_CHLA_regridded" = "Chlorophyll-a",
  "climatology_A_PAR_regridded" = "PAR",
  "climatology_M_0_0" = "Mixed Layer Depth",
  "climatology_Nstar_0_50" = "N*",
  "climatology_O_0_50" = "Oxygen Saturation",
  "log10_climatology_S_PP_regridded" = "Primary Production",
  "log10_climatology_TOT_POC_CMEMS" = "Total POC",
  "climatology_o_0_50" = "Dissolved O2",
  "log10_climatology_eke_aviso" = "Eddy Kinetic Energy",
  "climatology_fsle_aviso_2001_2020" = "FSLE",
  "climatology_i_0_50" = "Silicate",
  "climatology_n_0_50" = "Nitrate",
  "climatology_p_0_50" = "Phosphate",
  "climatology_s_0_50" = "Salinity",
  "climatology_t_0_50" = "Temperature",
  "log10_climatology_i_0_50" = "Silicate",
  "log10_climatology_n_0_50" = "Nitrate",
  "log10_climatology_p_0_50" = "Phosphate"
)

# Add a new column with the adjusted predictor names
df_all <- df_all %>%
  mutate(adjusted_variable = dplyr::recode(variable, !!!column_name_map))

# Step 1: Get top 3 predictors for each individual taxon (simple_name) with IQR
top3_per_taxon_with_iqr <- df_all %>%
  group_by(simple_name, adjusted_variable) %>%
  summarize(
    median_value = median(value, na.rm = TRUE),
    q25 = quantile(value, 0.25, na.rm = TRUE),
    q75 = quantile(value, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(simple_name, desc(median_value)) %>%
  group_by(simple_name) %>%
  slice_head(n = 3) %>%
  mutate(rank = row_number()) %>%
  ungroup()

# Add group information to the data
top3_per_taxon_with_iqr <- top3_per_taxon_with_iqr %>%
  left_join(group_assignments, by = "simple_name") %>%
  filter(!is.na(group))  # Remove any clades not in our groups

# Calculate the distances for asymmetric error bars
top3_per_taxon_with_iqr <- top3_per_taxon_with_iqr %>%
  mutate(
    lower_distance = median_value - q25,  # Distance from median to Q25
    upper_distance = q75 - median_value,  # Distance from median to Q75
    max_distance = max(c(lower_distance, upper_distance), na.rm = TRUE)  # For scaling
  )

# Check the updated dataframe
head(top3_per_taxon_with_iqr)

# Define colors for adjusted variables based on categories
color_map <- c(
  "Mixed Layer Depth" = "#6c79dc",   # Blue (Oceanographic)
  "FSLE" = "#0925f6",                # Greenish-Blue (Oceanographic)
  "Eddy Kinetic Energy" = "#20128f", # Dark Blue (Oceanographic)
  
  "Temperature" = "#E41A1C",         # Red (Energy-related)
  "PAR" = "#FF4500",                 # Orange-Red (Energy-related)
  "Primary Production" = "#FF6347",   # Tomato (Energy-related)
  
  "Total POC" = "#4DAF4A",           # Green (Biological)
  "Chlorophyll-a" = "#00A087",       # Teal Green (Biological)
  
  "Silicate" = "#F781BF",            # Pink (Nutrient-related)
  "Nitrate" = "#E377C2",             # Soft Pink (Nutrient-related)
  "Phosphate" = "#FF69B4",           # Hot Pink (Nutrient-related)
  "N*" = "#D147A3",                  # Magenta (Nutrient-related)
  "AOU" = "#DDA0DD",                 # Plum (Nutrient-related)
  
  "Salinity" = "#A9A9A9",            # Grey (Salinity)
  "Oxygen Saturation" = "#FF8C00",   # Orange (Oxygen)
  "Dissolved O2" = "#FFA500"         # Gold (Oxygen)
)

# Bubble plot with asymmetric error bars
bubble_plot_errorbar_asymmetric <- ggplot(top3_per_taxon_with_iqr, aes(x = rank, y = simple_name)) +
  # Asymmetric error bars - left side shows distance to Q25, right side shows distance to Q75
  geom_errorbarh(aes(xmin = rank - (lower_distance / max_distance) * 0.3,
                     xmax = rank + (upper_distance / max_distance) * 0.3,
                     color = adjusted_variable),
                 height = 0.15, alpha = 0.6, linewidth = 1.5) +
  geom_point(aes(color = adjusted_variable, size = median_value), 
             alpha = 0.6) +  # Lower alpha for better text contrast
  # Add median value text inside the bubbles with black text
  geom_text(aes(label = round(median_value, 0)),
            size = 2.5, 
            hjust = 0.5, 
            vjust = 0.5,
            color = "black",
            fontface = "bold") +
  # Annotate IQR range at the end of error bars
  geom_text(aes(x = rank + (upper_distance / max_distance) * 0.3 + 0.08,
                y = simple_name,
                label = paste0(round(q25, 0), "-", round(q75, 0)),
                color = adjusted_variable),
            size = 1.8,
            hjust = 0,
            vjust = 0.5,
            fontface = "bold",
            show.legend = FALSE) +
  facet_wrap(~ group, scales = "free_y", ncol = 1) +
  scale_x_continuous(breaks = c(1, 2, 3), labels = c("1st", "2nd", "3rd")) +
  scale_size_continuous(
    name = "Median\nImportance",
    range = c(3, 10),
    breaks = c(5, 10, 15, 20),
    labels = c("5", "10", "15", "20")
  ) +
  scale_color_manual(
    name = "Environmental\nPredictor",
    values = color_map,
    na.value = "grey50"
  ) +
  theme_classic() +
  labs(
    x = "Predictor Rank",
    y = "Taxon",
    title = "Top 3 Environmental Predictors per Taxon (Asymmetric error bars show skew)"
  ) +
  theme(
    text = element_text(size = 7),
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(size = 7),
    axis.title.x = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    strip.text = element_text(size = 6, face = "bold"),
    plot.background = element_rect(fill = "white"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
    legend.position = "bottom",
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 5),
    legend.key.size = unit(0.3, "cm"),
    legend.key.height = unit(0.25, "cm"),
    legend.key.width = unit(0.25, "cm"),
    legend.margin = margin(t = 2, r = 2, b = 2, l = 2),
    legend.box = "horizontal",
    plot.title = element_text(size = 7, hjust = 0.5, face = "bold"),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
    panel.spacing = unit(0.3, "lines")
  ) +
  guides(
    color = guide_legend(
      override.aes = list(size = 4, alpha = 1),
      ncol = 4,
      title.position = "top",
      title.hjust = 0.5,
      keywidth = unit(0.3, "cm"),
      keyheight = unit(0.25, "cm")
    ),
    size = guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      ncol = 4,
      keywidth = unit(0.3, "cm"),
      keyheight = unit(0.25, "cm")
    )
  )

print(bubble_plot_errorbar_asymmetric)

# Save the asymmetric version
ggsave(
  filename = paste0(wd_out, "/Top3_predictors_bubble_plot_errorbar__enhanced_asymmetric_v2.svg"),
  plot = bubble_plot_errorbar_asymmetric,
  width = 4,
  height = 8,  # Reduced height since we now have 4 groups instead of 5
  dpi = 300,
  device = "svg"
)

## ====================================================================================================
## END OF SCRIPT
## ====================================================================================================
