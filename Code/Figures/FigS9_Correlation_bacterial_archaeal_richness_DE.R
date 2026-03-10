## ====================================================================================================
## This R script performs comprehensive correlation analysis between bacterial and archaeal diversity
## patterns across global ocean basins. The analysis examines three diversity indices (richness,
## Shannon diversity, and Chao1 estimates) to assess the relationship between bacterial and archaeal
## communities in both global oceans (0-90° latitude) and tropical/subtropical regions (0-40° latitude).
## The script generates hexagonal density plots, calculates Pearson correlations with statistical
## significance testing, and creates publication-ready visualizations to examine potential coupling
## or decoupling of bacterial and archaeal diversity patterns.
##
## Author:       Dominic Eriksson
## Date:         27th of February 2026
## Affiliation:  Environmental Physics Group UP, ETH Zürich, Switzerland
## Contact:      deriksson@ethz.ch
##
## Input files:
##   - Domain Archaea richness RData with ensemble mean predictions
##   - Domain Bacteria richness RData with ensemble mean predictions
##   - Domain Archaea Shannon diversity RData with ensemble statistics
##   - Domain Bacteria Shannon diversity RData with ensemble statistics
##   - Domain Archaea Chao1 estimates RData with diversity predictions
##   - Domain Bacteria Chao1 estimates RData with diversity predictions
##
## Output files:
##   - SVG hexagonal density plot showing bacterial-archaeal richness correlations
##   - Correlation table CSV with statistics for all diversity indices and regions
##   - Console output with comprehensive correlation analysis results
##
## Strategy:
##   The script loads ensemble mean diversity predictions for both Bacteria and Archaea domains
##   across three diversity indices (richness, Shannon, Chao1). For each index, spatial data
##   is joined by coordinates and filtered for finite values. Regional analysis compares
##   global patterns (0-90° latitude) with tropical/subtropical patterns (0-40° latitude).
##   Pearson correlations with significance testing examine potential coupling between domains.
##   Hexagonal density plots visualize richness relationships with linear regression fits
##   and statistical annotations. Comprehensive correlation tables enable comparison across
##   all diversity indices and regions to assess bacterial-archaeal diversity coupling.
##
## Required R packages (tested versions):
##   - dplyr           1.1.2
##   - ggplot2         3.4.2
##   - ggpmisc         0.5.2
##   - hexbin          1.28.3
##   - raster          3.6.26
## ====================================================================================================

### ====================================================================================================================
### Preparation
### ====================================================================================================================

# Clear workspace
rm(list = ls())

# Load required libraries
library(dplyr)          # For data manipulation and filtering operations
library(ggplot2)        # For creating publication-quality plots and visualizations
library(ggpmisc)        # For statistical annotations on ggplot2 (stat_poly_eq)
library(hexbin)         # For hexagonal binning in density plots
library(raster)         # For handling and processing raster data

# Directories
wd_out <- "./Output/"

## Load all diversity indices data --------------------------------------------

# Load Richness data
r_archaea <- get(load("../2_Extract_model_outputs/6_Output/Non_log/domain_Archaea_Richness_5000_v3.RData"))
r_bacteria <- get(load("../2_Extract_model_outputs/6_Output/Non_log/domain_Bacteria_Richness_5000_v3.RData"))

# Load Shannon data
s_archaea <- get(load("../2_Extract_model_outputs/6_Output/Non_log/domain_Archaea_Shannon_5000_v3.RData"))
s_bacteria <- get(load("../2_Extract_model_outputs/6_Output/Non_log/domain_Bacteria_Shannon_5000_v3.RData"))

# Load Chao1 data
c_archaea <- get(load("../2_Extract_model_outputs/6_Output/Non_log/domain_Archaea_Chao1_5000_v3.RData"))
c_bacteria <- get(load("../2_Extract_model_outputs/6_Output/Non_log/domain_Bacteria_Chao1_5000_v3.RData"))

# Extract ensemble means
r_archaea <- r_archaea[[1]]
r_bacteria <- r_bacteria[[1]]
s_archaea <- s_archaea[[1]]
s_bacteria <- s_bacteria[[1]]
c_archaea <- c_archaea[[1]]
c_bacteria <- c_bacteria[[1]]

## Prepare data for each diversity index ------------------------------------

# Richness
df_richness <- inner_join(
  as.data.frame(r_archaea$Annual_Mean, xy = TRUE) %>% rename(Archaea = Annual_Mean),
  as.data.frame(r_bacteria$Annual_Mean, xy = TRUE) %>% rename(Bacteria = Annual_Mean),
  by = c("x", "y")
) %>%
  filter(is.finite(Archaea), is.finite(Bacteria))

# Shannon
df_shannon <- inner_join(
  as.data.frame(s_archaea$Annual_Mean, xy = TRUE) %>% rename(Archaea = Annual_Mean),
  as.data.frame(s_bacteria$Annual_Mean, xy = TRUE) %>% rename(Bacteria = Annual_Mean),
  by = c("x", "y")
) %>%
  filter(is.finite(Archaea), is.finite(Bacteria))

# Chao1
df_chao1 <- inner_join(
  as.data.frame(c_archaea$Annual_Mean, xy = TRUE) %>% rename(Archaea = Annual_Mean),
  as.data.frame(c_bacteria$Annual_Mean, xy = TRUE) %>% rename(Bacteria = Annual_Mean),
  by = c("x", "y")
) %>%
  filter(is.finite(Archaea), is.finite(Bacteria))

### ====================================================================================================================
### Correlation Analysis and Visualization
### ====================================================================================================================

## Calculate correlations for all indices and regions ----------------------

# Function to calculate correlation and format results
calc_correlation <- function(df, diversity_index, region_name) {
  cor_result <- cor.test(df$Archaea, df$Bacteria)
  
  # Significance stars
  p_stars <- if (cor_result$p.value < 0.001) "***"
            else if (cor_result$p.value < 0.01) "**"
            else if (cor_result$p.value < 0.05) "*"
            else ""
  
  data.frame(
    Diversity_Index = diversity_index,
    Region = region_name,
    Correlation = round(cor_result$estimate, 3),
    P_value = round(cor_result$p.value, 6),
    Significance = p_stars,
    N_points = nrow(df)
  )
}

# Calculate correlations for all combinations
correlation_results <- bind_rows(
  # Richness
  calc_correlation(df_richness, "Richness", "Global (0-90°)"),
  calc_correlation(df_richness %>% filter(abs(y) <= 40), "Richness", "Tropics/Subtropics (0-40°)"),
  
  # Shannon
  calc_correlation(df_shannon, "Shannon", "Global (0-90°)"),
  calc_correlation(df_shannon %>% filter(abs(y) <= 40), "Shannon", "Tropics/Subtropics (0-40°)"),
  
  # Chao1
  calc_correlation(df_chao1, "Chao1", "Global (0-90°)"),
  calc_correlation(df_chao1 %>% filter(abs(y) <= 40), "Chao1", "Tropics/Subtropics (0-40°)")
)

## Create and display the correlation table ---------------------------------

cat("Correlation between Archaeal and Bacterial Diversity Indices\n")
cat("=============================================================\n\n")
print(correlation_results, row.names = FALSE)

# Create richness plot (as before)
df_plot <- bind_rows(
  df_richness %>% mutate(Region = "Global (0-90°)"),
  df_richness %>% filter(abs(y) <= 40) %>% mutate(Region = "Tropics/Subtropics (0-40°)")
)

facet_labels <- data.frame(
  Region = c("Global (0-90°)", "Tropics/Subtropics (0-40°)"),
  label = c(
    sprintf("r = %.3f %s", 
            correlation_results$Correlation[correlation_results$Diversity_Index == "Richness" & correlation_results$Region == "Global (0-90°)"],
            correlation_results$Significance[correlation_results$Diversity_Index == "Richness" & correlation_results$Region == "Global (0-90°)"]),
    sprintf("r = %.3f %s", 
            correlation_results$Correlation[correlation_results$Diversity_Index == "Richness" & correlation_results$Region == "Tropics/Subtropics (0-40°)"],
            correlation_results$Significance[correlation_results$Diversity_Index == "Richness" & correlation_results$Region == "Tropics/Subtropics (0-40°)"])
  )
)

# Plot 
p <- ggplot(df_plot, aes(x = Archaea, y = Bacteria)) +
  geom_hex(bins = 30, alpha = 0.8) +
  scale_fill_viridis_c(
    name = "Count",
    option = "inferno",
    trans = "log10"
  ) +
  geom_smooth(method = "lm", color = "blue", se = FALSE, size = 0.7) +
  facet_wrap(~Region, scales = "free") +
  geom_text(
    data = facet_labels,
    aes(x = -Inf, y = Inf, label = label),
    hjust = -0.1, vjust = 1.2, color = "red", size = 8 / .pt, inherit.aes = FALSE
  ) +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold", size = 8),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.position = "bottom"
  ) +
  labs(
    x = "Archaeal Richness",
    y = "Bacterial Richness",
    title = "Correlation between Archaeal and Bacterial Richness"
  )

# Print the plot
print(p)

# Save plot and correlation table
ggsave(
    paste0(wd_out, "correlation_richness_bacteria_archaea_v3.svg"), 
    plot = p, width = 3, height = 3, 
    units = "in", dpi = 300)

# Save correlation table as CSV
write.csv(correlation_results, 
          paste0(wd_out, "correlation_table_all_diversity_indices_v2.csv"), 
          row.names = FALSE)

cat("\n\nCorrelation table saved to:", paste0(wd_out, "correlation_table_all_diversity_indices.csv"))


## ====================================================================================================
## END OF SCRIPT
## ====================================================================================================
