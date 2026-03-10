
## ====================================================================================================
## R script to extract successfully modeled taxa from CEPHALOPOD SDM outputs. CEPHALOPOD creates 
## output folders for all targeted taxa, but some modeling runs may fail; log files contain 
## information on why a model did not complete successfully. 
##
## Author:       Dominic Eriksson
## Date:         8th of October 2025
## Affiliation:  Environmental Physics Group, UP, ETH Zürich, Switzerland
## Contact:      deriksson@ethz.ch
##
## Input files:  directories_with_maps.txt - list of directories containing CEPHALOPOD model 
##               outputs with successful modeling runs (filtered to rarefaction depth 5000)
##
## Output files: modeling_status_heatmap.png - heatmap visualization of modeling status by
##               taxonomic clade and diversity metric
##               clade_diversity_modeling_status.csv - summary table of modeling status
##               presence_absence_table_clades_metrics.csv - detailed presence/absence table
##
## Strategy:
##   The script identifies which taxa have successfully generated diversity metric models by 
##   parsing directory names from the directories_with_maps.txt file. It filters for 5000 
##   rarefaction depth directories, extracts clade and diversity metric information, creates
##   presence/absence matrices, and generates both tabular summaries and heatmap visualizations
##   showing modeling success across all taxonomic groups and diversity indices.
##
## ====================================================================================================

# Directories
wd_out <- "./Code/2_Extract_model_outputs/5_Output/"

# Load txt file
directories <- readLines("./Code/2_Extract_model_outputs/5_Output/directories_with_maps.txt")

# Remove all directories containing "Archive" in the path
directories <- directories[!grepl("Archive", directories)]

# Only keep the 5000 rarefaction depth directories
directories <- directories[grepl("_5000$", directories)]

# Check the results
cat("Number of directories after filtering:", length(directories), "\n")
head(directories)

# Create a table of successfully modeled clades and diversity metrics
library(dplyr)

# Extract the basename (last part of the path)
basenames <- basename(directories)

# Parse the clade and diversity metric information
parse_info <- data.frame(
    full_path = directories,
    basename = basenames,
    stringsAsFactors = FALSE
)

# Extract clade and diversity metric from basename
# Pattern: clade_DiversityMetric_rarefactionDepth
parse_info$clade <- gsub("_[^_]*_[0-9]+$", "", parse_info$basename)
parse_info$diversity_metric <- gsub(".*_([^_]+)_[0-9]+$", "\\1", parse_info$basename)
parse_info$rarefaction <- gsub(".*_([0-9]+)$", "\\1", parse_info$basename)

# Check the parsing
head(parse_info)

# Get unique clades and metrics
unique_clades <- sort(unique(parse_info$clade))
unique_metrics <- sort(unique(parse_info$diversity_metric))

cat("Unique clades found:", length(unique_clades), "\n")
cat("Unique metrics found:", length(unique_metrics), "\n")
print(unique_metrics)

# Create presence/absence table
presence_table <- matrix(0, 
                        nrow = length(unique_clades), 
                        ncol = length(unique_metrics),
                        dimnames = list(unique_clades, unique_metrics))

# Fill in the table
for(i in 1:nrow(parse_info)) {
    clade <- parse_info$clade[i]
    metric <- parse_info$diversity_metric[i]
    presence_table[clade, metric] <- 1
}

# Convert to data frame for better display
presence_df <- as.data.frame(presence_table)
presence_df$clade <- rownames(presence_df)
presence_df <- presence_df[, c("clade", unique_metrics)]

# Display the table
cat("\nPresence/Absence Table for Diversity Metrics:\n")
print(presence_df)

# Summary statistics
cat("\nSummary:\n")
cat("Total clades:", nrow(presence_df), "\n")
for(metric in unique_metrics) {
    count <- sum(presence_df[, metric])
    cat(paste0(metric, ": ", count, " clades modeled\n"))
}

# Create visualization
library(ggplot2)
library(tidyr)

# Reshape data for plotting
plot_data <- presence_df %>%
    pivot_longer(cols = -clade, names_to = "diversity_metric", values_to = "modeled") %>%
    mutate(
        modeled_factor = factor(modeled, levels = c(0, 1), labels = c("Not Modeled", "Modeled")),
        clade_clean = gsub("domain_|class_|phylum_", "", clade)
    )

# Create heatmap visualization
p_heatmap <- ggplot(plot_data, aes(x = diversity_metric, y = clade_clean, fill = modeled_factor)) +
    geom_tile(color = "white", size = 0.5) +
    scale_fill_manual(
        values = c("Not Modeled" = "lightgray", "Modeled" = "darkblue"),
        name = "Status"
    ) +
    labs(
        title = "Modeling Status of Diversity Metrics by Taxonomic Clade",
        subtitle = "Overview of successfully modeled diversity indices",
        x = "Diversity Metric",
        y = "Taxonomic Clade",
        caption = "Blue = Successfully modeled, Gray = Not modeled"
    ) +
    theme_minimal(base_size = 12, base_family = "Arial") +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5, color = "gray30"),
        plot.caption = element_text(size = 10, color = "gray50"),
        legend.position = "bottom",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "white"),
        strip.background = element_rect(fill = "gray95", color = "black")
    )

# Display the plot
print(p_heatmap)

# Save the plot
ggsave(
    filename = paste0(wd_out, "modeling_status_heatmap.png"),
    plot = p_heatmap,
    width = 10,
    height = max(8, nrow(presence_df) * 0.3),
    dpi = 300
)

# Save the table as CSV
write.csv(presence_df, paste0(wd_out, "clade_diversity_modeling_status.csv"), row.names = FALSE)

cat("\nFiles saved:\n")
cat("- Visualization: modeling_status_heatmap.png\n")
cat("- Data table: clade_diversity_modeling_status.csv\n")
# Save the table to a CSV file
write.csv(presence_df, file = paste0(wd_out, "presence_absence_table_clades_metrics.csv"), row.names = FALSE)
cat("\nPresence/Absence table saved to:", paste0(wd_out, "presence_absence_table_clades_metrics.csv\n"))

# ===============================================================================
# END OF SCRIPT
# ===============================================================================

