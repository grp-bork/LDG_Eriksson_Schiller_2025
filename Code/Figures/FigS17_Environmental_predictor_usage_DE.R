##====================================================================================================
## Supplementary Figure S17: Environmental Predictor Usage Analysis
##
## This R script performs comprehensive analysis of environmental predictor usage patterns across 
## prokaryotic taxonomic groups in species distribution modeling. The analysis extracts predictor 
## selection from final random forest models, quantifies predictor importance rankings, and creates 
## publication-ready visualizations showing which environmental variables are most important for 
## different prokaryotic clades. The script generates multiple plot types including heatmaps and 
## dot plots with importance scaling to reveal predictor usage patterns and ecological insights.
##
## Author:      Dominic Eriksson
##              Environmental Physics Group, UP
##              ETH Zürich, Switzerland
## Contact:     deriksson@ethz.ch
## Date:        February 27th, 2026
## Affiliation: ETH Zürich, Environmental Physics Group, UP
##
## Input files:
## - /net/sea/work/deriksson/Projects/Prokaryotes_LDG_v3/Code/4_Compute_ensembles_uncertainties/1_Output/final_selected_files_richness_5000_v2.txt
## - /net/sea/work/deriksson/Projects/Prokaryotes_LDG_v3/Code/2_Extract_model_outputs/1_Output/Non_log/Predictor_rankings_MLD_v6.csv
## - /net/sea/work/deriksson/Projects/Prokaryotes_LDG_v3/Code/2_Extract_model_outputs/1_Output/Log/Predictor_rankings_MLD_v6.csv
## - /net/sea/work/deriksson/Projects/Prokaryotes_LDG_v3/Code/Cephalopod_output/*/MODEL.RData (Multiple model files)
##
## Output files:
## - predictor_usage_heatmap.png (Binary heatmap showing predictor usage by taxonomic group)
## - predictor_usage_dotplot.png (Colored dot plot with predictor categories and importance scaling)
## - predictor_usage_dotplot_black.png/.pdf/.svg (Publication-ready black dot plot with importance scaling)
##
## Strategy:
## 1. Load and process predictor ranking data from Log and Non_log model versions
## 2. Filter for domain and class level taxonomic groups of interest
## 3. Create taxonomic group assignments and apply readable predictor name mapping
## 4. Calculate summary statistics (median, IQR) for predictor importance by taxon
## 5. Extract final model formulas from CEPHALOPOD model output files
## 6. Map predictor variables to informative names using standardized nomenclature
## 7. Create presence/absence matrix showing predictor usage across groups
## 8. Generate heatmap visualization of predictor selection patterns
## 9. Create colored dot plots with categorical predictor coloring and importance scaling
## 10. Generate publication-ready black dot plots with statistical annotations
##
## Required R packages:
## - stringr (1.4.0): String manipulation and pattern matching operations
## - ggplot2 (3.3.6): Advanced data visualization and publication-quality plotting
## - dplyr (1.0.9): Data manipulation, filtering, and transformation operations
## - tidyr (1.1.3): Data tidying and reshaping for visualization workflows
##====================================================================================================

### ==============================================================================================
### Predictor ranking boxplot for clades
### ==============================================================================================

# Load libraries
library(stringr)
library(ggplot2)
library(dplyr)

# Clear workspace
rm(list = ls())

# Directories
wd_out <- "./Output/" #your output directory here

### ==============================================================================================
### Bubble plots for top 3 predictors per taxon
### ==============================================================================================

# Load the final selected files from the text file
fnames_model_stats <- readLines("./Code/4_Compute_ensembles_uncertainties/1_Output/final_selected_files_richness_5000_v2.txt")

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

log_clade_names <- extract_clade(log_clades)
nonlog_clade_names <- extract_clade(nonlog_clades)

# Print to check what we extracted
cat("Log clades:\n")
print(log_clade_names)
cat("\nNon-log clades:\n")
print(nonlog_clade_names)

# Load data
df_nonLog <- read.csv("./Code/2_Extract_model_outputs/1_Output/Non_log/Predictor_rankings_MLD_v6.csv")
df_Log <- read.csv("./Code/2_Extract_model_outputs/1_Output/Log/Predictor_rankings_MLD_v6.csv")

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
  group = c(
    "group3",  # domain_Bacteria (added to Group 2)
    "group2",  # domain_Archaea (added to Group 4)
    "group2",  # Poseidoniia (moved to Group 4)
    "group1",  # Nitrososphaeria
    "group3",  # Acidimicrobiia (moved to Group 2)
    "group2",  # Bacteroidia (moved to Group 4)
    "group4",  # Chlamydiia (moved to Group 3)
    "group3",  # Cyanobacteriia (moved to Group 2)
    "group2",  # UBA1144 (moved to Group 4)
    "group1",  # Fibrobacteria
    "group2",  # Marinisomatia (moved to Group 4)
    "group4",  # UBA4151 (moved to Group 3)
    "group4",  # UBA796 (moved to Group 3)
    "group4",  # XYA12-FULL-58-9 (moved to Group 3)
    "group3",  # Alphaproteobacteria (moved to Group 2)
    "group2",  # Gammaproteobacteria (moved to Group 4)
    "group2"   # SAR324 (moved to Group 4)
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


# Step 1: Get summary statistics for ALL predictors for each individual taxon (simple_name) with IQR
all_predictors_with_iqr <- df_all %>%
  group_by(simple_name, adjusted_variable) %>%
  summarize(
    median_value = median(value, na.rm = TRUE),
    q25 = quantile(value, 0.25, na.rm = TRUE),
    q75 = quantile(value, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(simple_name, desc(median_value)) %>%
  group_by(simple_name) %>%
  mutate(rank = row_number()) %>%
  ungroup()


# Open file
model <- get(load("./Code/Cephalopod_output/Non_log/chunks11_2025_07_14_07_17/domain_Archaea_Richness_5000/MODEL.RData"))

sort(names(model$ENSEMBLE))

# Directories
wd_in <- "./Code/Cephalopod_output/"
wd_out <- "../Prokaryotes_LDG_v3/Revisions_Cell_Host_and_Microbes/Code/2_Output/"

# Get filenames of groups of interest
groups <- c(
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
    )

# Find complete paths for directories matching the group names
group_paths <- list()

# Search for directories recursively
all_dirs <- list.dirs(wd_in, recursive = TRUE, full.names = TRUE)

# Find matches for each group
for(group in groups) {
    # Find directories that contain the group name
    matching_dirs <- all_dirs[grepl(group, basename(all_dirs))]
    
    if(length(matching_dirs) > 0) {
        group_paths[[group]] <- matching_dirs
        cat("Found", length(matching_dirs), "path(s) for", group, ":\n")
        for(path in matching_dirs) {
            cat("  ", path, "\n")
        }
        cat("\n")
    } else {
        cat("No paths found for", group, "\n\n")
    }
}

# Create a named vector with all found paths
all_group_paths <- unlist(group_paths)
names(all_group_paths) <- rep(names(group_paths), sapply(group_paths, length))

# Add 05_standard_maps.pdf to the paths
all_group_paths <- paste0(all_group_paths, "/05_standard_maps.pdf")

# Check file existence and prefer Non_log over Log versions
final_paths <- c()

# Debug: Check the names
# cat("Names of all_group_paths:\n")
# print(names(all_group_paths))

# Use the groups vector directly
for(group in groups) {
    # Get all paths for this specific group - search in the paths themselves, not the names
    group_paths_all <- all_group_paths[grepl(group, all_group_paths)]
    
    # Check which files exist for this group
    group_existing <- group_paths_all[file.exists(group_paths_all)]
    
    if(length(group_existing) > 0) {
        # Check if there's a Non_log version
        nonlog_file <- group_existing[grepl("/Non_log/", group_existing)]
        log_file <- group_existing[grepl("/Log/", group_existing)]
        
        if(length(nonlog_file) > 0) {
            # Use Non_log version
            selected_file <- nonlog_file[1]
            cat("Using Non_log file for", group, "\n")
        } else if(length(log_file) > 0) {
            # Use Log version
            selected_file <- log_file[1]
            cat("Using Log file for", group, "\n")
        } else {
            # Use first available
            selected_file <- group_existing[1]
            cat("Using file for", group, "\n")
        }
        
        final_paths <- c(final_paths, selected_file)
        names(final_paths)[length(final_paths)] <- group
    } else {
        cat("No existing files found for", group, "\n")
    }
}

cat("\nFinal selected paths:\n")
print(final_paths)
class(final_paths)
length(final_paths)

# Replace the 5_standard_maps.pdf with MODEL.RData
final_paths <- sub("05_standard_maps.pdf", "MODEL.RData", final_paths)

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

# Extract formulas from all model files
group_formulas <- list()

cat("\n=== EXTRACTING FORMULAS FOR ALL GROUPS ===\n\n")

for(i in 1:length(final_paths)) {
    group_name <- names(final_paths)[i]
    model_path <- final_paths[i]
    
    cat("Processing", group_name, "...\n")
    cat("Model path:", model_path, "\n")
    
    # Check if file exists
    if(file.exists(model_path)) {
        # Load the model
        tryCatch({
            model_data <- get(load(model_path))
            
            # Extract the formula
            if(!is.null(model_data$RF$final_wf$fit$actions$model$formula)) {
                rf_formula <- model_data$RF$final_wf$fit$actions$model$formula
                
                # Extract predictor names from the formula
                predictor_vars <- all.vars(rf_formula)[-1]  # Remove response variable (first one)
                
                # Map to informative names
                informative_names <- sapply(predictor_vars, function(x) {
                    if(x %in% names(column_name_map)) {
                        return(column_name_map[x])
                    } else {
                        return(x)  # Keep original name if not in map
                    }
                })
                
                # Store results
                group_formulas[[group_name]] <- list(
                    formula = rf_formula,
                    predictors = predictor_vars,
                    informative_names = as.character(informative_names)
                )
                
                cat("  SUCCESS: Found", length(predictor_vars), "predictors\n")
                cat("  Predictors:", paste(informative_names, collapse = ", "), "\n")
                
            } else {
                cat("  WARNING: No formula found in RF model\n")
            }
        }, error = function(e) {
            cat("  ERROR loading model:", e$message, "\n")
        })
    } else {
        cat("  ERROR: File does not exist\n")
    }
    cat("\n")
}

# Print summary of all formulas
cat("\n=== SUMMARY OF PREDICTORS BY GROUP ===\n\n")
for(group_name in names(group_formulas)) {
    cat(group_name, ":\n")
    predictors <- group_formulas[[group_name]]$informative_names
    cat("  ", paste(predictors, collapse = ", "), "\n\n")
}

# Create visualization of predictor usage across groups
library(ggplot2)
library(dplyr)
library(tidyr)

# Create a matrix showing which predictors are used by which groups
all_predictors <- unique(unlist(lapply(group_formulas, function(x) x$informative_names)))
groups_clean <- names(group_formulas)

# Create presence/absence matrix
predictor_matrix <- matrix(0, nrow = length(groups_clean), ncol = length(all_predictors))
rownames(predictor_matrix) <- groups_clean
colnames(predictor_matrix) <- all_predictors

# Fill the matrix
for(i in 1:length(groups_clean)) {
    group <- groups_clean[i]
    if(group %in% names(group_formulas)) {
        predictors <- group_formulas[[group]]$informative_names
        predictor_matrix[i, predictors] <- 1
    }
}

# Convert to data frame for ggplot
plot_data <- as.data.frame(predictor_matrix) %>%
    mutate(Group = rownames(.)) %>%
    pivot_longer(cols = -Group, names_to = "Predictor", values_to = "Used") %>%
    mutate(
        Group = factor(Group),
        Predictor = factor(Predictor, levels = all_predictors),
        Used = as.logical(Used)
    ) %>%
    # Add predictor colors
    mutate(Predictor_Color = case_when(
        Predictor %in% names(color_map) ~ color_map[Predictor],
        TRUE ~ "#808080"  # Default gray for unmapped predictors
    ))

# Clean group names for better display
plot_data$Group_Clean <- gsub("domain_|class_", "", plot_data$Group)
plot_data$Group_Clean <- gsub("_Richness_5000", "", plot_data$Group_Clean)
plot_data$Group_Clean <- gsub("_", " ", plot_data$Group_Clean)

# Join with median importance values
plot_data <- plot_data %>%
    left_join(all_predictors_with_iqr, 
              by = c("Group_Clean" = "simple_name", "Predictor" = "adjusted_variable")) %>%
    # For predictors not in all_predictors_with_iqr, set median_value to 0
    mutate(median_value = ifelse(is.na(median_value), 0, median_value))

# Create the heatmap
p1 <- ggplot(plot_data, aes(x = Predictor, y = Group_Clean, fill = Used)) +
    geom_tile(color = "white", size = 0.5) +
    scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "#2E86AB"), 
                     labels = c("FALSE" = "Not used", "TRUE" = "Used")) +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 11, face = "bold"),
        panel.grid = element_blank(),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    ) +
    labs(
        title = "Environmental Predictors Used in Species Distribution Models",
        subtitle = "Prokaryotic Groups (Bacteria & Archaea)",
        x = "Environmental Predictors",
        y = "Taxonomic Groups",
        fill = "Predictor Usage"
    )

# Save the plot
ggsave(paste0(wd_out, "predictor_usage_heatmap.png"), p1, 
       width = 14, height = 8, dpi = 300, bg = "white")

# Alternative: Dot plot with colored predictors by category (size scaled by importance)
p2 <- plot_data %>%
    filter(Used == TRUE) %>%
    ggplot(aes(x = Predictor, y = Group_Clean)) +
    geom_point(aes(color = Predictor, size = median_value), alpha = 0.8) +
    scale_color_manual(values = color_map, guide = "none") +
    scale_size_continuous(range = c(1, 8), name = "Median\nImportance") +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.position = "right"
    ) +
    labs(
        title = "Environmental Predictors by Taxonomic Group",
        subtitle = "Colors represent predictor categories, size represents median importance",
        x = "Environmental Predictors",
        y = "Taxonomic Groups"
    )

# Save the dot plot
ggsave(paste0(wd_out, "predictor_usage_dotplot.png"), p2, 
       width = 14, height = 8, dpi = 300, bg = "white")

# Alternative: Dot plot with black points only - Publication quality (size scaled by importance)
p3 <- plot_data %>%
    filter(Used == TRUE) %>%
    ggplot(aes(x = Predictor, y = Group_Clean, size = median_value)) +
    geom_point(alpha = 0.9, color = "black") +
    scale_size_continuous(range = c(1, 6), name = "Median\nImportance") +
    theme_minimal(base_size = 12) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 11, color = "black"),
        axis.text.y = element_text(size = 11, color = "black"),
        axis.title = element_text(size = 13, face = "bold", color = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey90", size = 0.3),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5, margin = margin(b = 10)),
        plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 15)),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10)),
        plot.margin = margin(15, 15, 15, 15),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        legend.position = "right"
    ) +
    labs(
        title = "Environmental predictors used in prokaryotic habitat models",
        subtitle = "Dot size represents median predictor importance",
        x = "Environmental predictors",
        y = "Taxonomic groups"
    )

# Save the black dot plot in multiple formats
ggsave(paste0(wd_out, "predictor_usage_dotplot_black.png"), p3, 
       width = 12, height = 8, dpi = 300, bg = "white")

ggsave(paste0(wd_out, "predictor_usage_dotplot_black.pdf"), p3, 
       width = 12, height = 8, device = "pdf", bg = "white")

ggsave(paste0(wd_out, "predictor_usage_dotplot_black.svg"), p3, 
       width = 12, height = 8, device = "svg", bg = "white")

# Create a summary table
predictor_summary <- plot_data %>%
    group_by(Predictor) %>%
    summarise(
        Groups_Using = sum(Used),
        Percentage = round(Groups_Using / length(unique(plot_data$Group_Clean)) * 100, 1),
        Color = first(Predictor_Color),
        .groups = "drop"
    ) %>%
    arrange(desc(Groups_Using))

cat("\n=== PREDICTOR USAGE SUMMARY ===\n")
print(predictor_summary)

cat("\nPlots saved to:", wd_out, "\n")
cat("- predictor_usage_heatmap.png\n")
cat("- predictor_usage_dotplot.png\n")

### ====================================================================================================================
### END OF SCRIPT
### ====================================================================================================================


