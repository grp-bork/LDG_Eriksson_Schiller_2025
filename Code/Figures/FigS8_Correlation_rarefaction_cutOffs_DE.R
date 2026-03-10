## ====================================================================================================
## This R script generates Figure S4, illustrating the correlation between
## diversity indices (Richness, Shannon, and Chao1) across varying rarefaction cut-offs (1000 vs. 5000).
##
## Author:       Dominic Eriksson
## Date:         20th of January 2026
## Affiliation:  Environmental Physics Group UP, ETH Zürich, Switzerland
## Contact:      deriksson@ethz.ch
##
## Input files:
##   - RData files containing annual ensemble mean predictions for Archaea and Bacteria diversity indices:
##   - Files saved in folder: Code/Figures/
##
## Output files:
##   - SVG figure showing hexbin scatter plots of Archaeal vs. Bacterial richness correlations.
##   - CSV table summarizing correlation coefficients, p-values, significance, and sample sizes for 
##     all indices and regions.
##
## Required R packages (tested versions):
##   - dplyr    1.1.2
##   - ggplot2  3.5.2
##   - GGally   2.1.2
##   - reshape  1.4.4
##   - raster   3.6.26
## ====================================================================================================

# Clear the workspace
rm(list = ls())

# Load libraries
library(ggplot2)
library(raster)
library(dplyr)
library(GGally)
library(reshape2)

# We use non-normalized data for the output, the aim is to compare the same diversity metric across different rarefied outputs. We 
# compare the annual ensemble

# Depth levels
depth_levels <- c("MLD", "0_10Meters")
d <- 1 # We only check the MLD level

# Directories
wd_in <- "Code/2_Extract_model_outputs/6_Output/Non_log/"
wd_out <- "Code/Figures/FigS4/"

### ===========================================================================
### Data formatting
### ===========================================================================

# Get filenames
fnames <- list.files( wd_in, full.names = TRUE, recursive = TRUE, pattern = ".RData", ignore.case = TRUE)

# Print the filenames for verification
print(fnames)

# Initialize an empty list to store the raster bricks
raster_list <- list()

# Loop through each file and read it as a brick
for (f in seq_along(fnames)) {
  
  # Open raster file
  r <- get(load(fnames[f]))

  # Extract annual ensemble mean
  r <- r[[1]]
  r <- r$Annual_Mean

  # If the actual clade name contains "." replace it with "_" - SECURITY STEP
  names(r) <- gsub(".RData", "", basename(fnames[f]), perl = TRUE)

  # Save in list object
  raster_list[[f]] <- r
}

# Stack all raster into one stack
raster_all <- do.call("stack", raster_list)

### Scatter matrix correlation between 5000 and 1000 rarefactioin levels diversity across domains

# Load necessary libraries
library(raster)
library(ggplot2)
library(dplyr)
library(reshape2)

# Stack all raster into one stack
raster_all <- do.call("stack", raster_list)

# Remove raster layers that contain a "."
# layers_to_keep <- grep("\\.", names(raster_all), invert = TRUE, value = TRUE)
# raster_all <- raster_all[[layers_to_keep]]

# Extract layers for 1000 and 5000 rarefied data
layers_1000 <- grep("_1000_v3$", names(raster_all), value = TRUE)
layers_5000 <- grep("_5000_v3$", names(raster_all), value = TRUE)

# Remove the rarefaction indication
clean_layers_1000 <- gsub("_1000_v3$", "", layers_1000)
clean_layers_5000 <- gsub("_5000_v3$", "", layers_5000)

# Extract common layers
common_layers <- intersect(clean_layers_1000, clean_layers_5000)

# Add the rarefaction indication back to the common layers
common_layers_1000 <- paste0(common_layers, "_1000_v3")
common_layers_5000 <- paste0(common_layers, "_5000_v3")

# Ensure the layers are in the same order for comparison
layers_1000 <- sort(common_layers)
layers_5000 <- sort(common_layers)

# Extract the corresponding layers
raster_1000 <- raster_all[[common_layers_1000]]
raster_5000 <- raster_all[[common_layers_5000]]

# Convert raster layers to dataframes
df_1000 <- as.data.frame(raster_1000, xy = TRUE)
df_5000 <- as.data.frame(raster_5000, xy = TRUE)

# Remove the rarefaction indication from the column names
colnames(df_1000) <- gsub("_1000_v3$", "", colnames(df_1000))
colnames(df_5000) <- gsub("_5000_v3$", "", colnames(df_5000))
 
# Melt dataframes for plotting
df_1000_melt <- melt(df_1000, id.vars = c("x", "y"), variable.name = "Layer", value.name = "Diversity_1000")
df_5000_melt <- melt(df_5000, id.vars = c("x", "y"), variable.name = "Layer", value.name = "Diversity_5000")

# Check the head of the melted dataframes
head(df_1000_melt)
head(df_5000_melt)

# Merge dataframes based on Layer, x, and y
merged_data <- merge(df_1000_melt, df_5000_melt, by = c("Layer", "x", "y"))

# Check the head of the merged dataframe
head(merged_data)

# Remove NAs
merged_data <- na.omit(merged_data)

# Remove rows that include Pielou in the clade column
merged_data <- merged_data[!grepl("Pielou", merged_data$Layer), ]

# Calculate Spearman correlations and statistical significance
cor_results <- merged_data %>%
  group_by(Layer) %>%
  summarize(
    spearman_cor = cor(Diversity_1000, Diversity_5000, method = "spearman"),
    p_value = cor.test(Diversity_1000, Diversity_5000, method = "spearman")$p.value
  ) %>%
  mutate(
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ ""
    ),
    label = paste0("rho = ", round(spearman_cor, 2), significance)
  )


# Load necessary libraries
library(ggplot2)
library(dplyr)
library(reshape2)

# Subset the data for Archaea and Bacteria
archaea_data <- merged_data[grepl("Archaea", merged_data$Layer), ]
bacteria_data <- merged_data[grepl("Bacteria", merged_data$Layer), ]

# Combine the subsets
combined_data <- rbind(archaea_data, bacteria_data)

# Remove rows with diversity values equal to 0
combined_df <- combined_data[combined_data$Diversity_1000 != 0 & combined_data$Diversity_5000 != 0, ]


# Extract the domain and diversity index from the Layer column
combined_df <- combined_df %>%
  mutate(Domain = ifelse(grepl("Archaea", Layer), "Archaea", "Bacteria"),
         Diversity_Index = case_when(
           grepl("Richness", Layer) ~ "Richness",
           grepl("Shannon", Layer) ~ "Shannon",
           grepl("Chao1", Layer) ~ "Chao1"
         ))

# Ensure cor_results has Domain and Diversity_Index columns
cor_results <- cor_results %>%
  mutate(Domain = ifelse(grepl("Archaea", Layer), "Archaea", "Bacteria"),
         Diversity_Index = case_when(
           grepl("Richness", Layer) ~ "Richness",
           grepl("Shannon", Layer) ~ "Shannon",
           grepl("Chao1", Layer) ~ "Chao1"
         ))

# Calculate annotation positions for each facet
positions <- combined_df %>%
  group_by(Domain, Diversity_Index) %>%
  summarize(
    x_pos = (max(Diversity_1000, na.rm = TRUE)/2 ),
    y_pos = max(Diversity_5000, na.rm = TRUE)
  )


# Join positions to cor_results_summary
cor_results_summary <- cor_results %>%
  group_by(Domain, Diversity_Index) %>%
  summarize(
    spearman_cor = mean(spearman_cor),
    p_value = mean(p_value),
    significance = first(significance),
    label = first(label),
    .groups = 'drop'
  ) %>%
  left_join(positions, by = c("Domain", "Diversity_Index"))

# Find the maximum count in the hexbin (after plotting once, or estimate)
# Or use a high value like 100000 if you want to cap the legend

max_count <- max(ggplot_build(gg_plot)$data[[1]]$count, na.rm = TRUE)
# Or set a fixed upper break, e.g. 10000 or 50000

breaks <- c(1, 10, 100, 1000, 10000, max_count)
labels <- c("1", "10", "100", "1,000", "10,000", paste0(">", 10000))

# Calculate Pearson correlations and statistical significance for each facet
cor_results_pearson <- combined_df %>%
  group_by(Domain, Diversity_Index) %>%
  summarize(
    pearson_cor = cor(Diversity_1000, Diversity_5000, method = "pearson"),
    p_value = cor.test(Diversity_1000, Diversity_5000, method = "pearson")$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ ""
    ),
    label = paste0("r = ", round(pearson_cor, 2), significance)
  )

# Calculate annotation positions for each facet
positions <- combined_df %>%
  dplyr::group_by(Domain, Diversity_Index) %>%
  dplyr::summarize(
    x_pos = max(Diversity_1000, na.rm = TRUE),
    y_pos = max(Diversity_5000, na.rm = TRUE),
    .groups = "drop"
  )

# Join positions to correlation results
cor_results_pearson <- dplyr::left_join(cor_results_pearson, positions, by = c("Domain", "Diversity_Index"))

# Plot with annotation
gg_plot <- ggplot(combined_df, aes(x = Diversity_1000, y = Diversity_5000)) +
  geom_hex(bins = 30) +
  scale_fill_viridis_c(
    option = "C",
    name = "log₁₀(Count)",
    trans = "log10",
    breaks = breaks,
    labels = labels,
    limits = c(1, max_count)
  ) +
  geom_smooth(method = "lm", se = FALSE, color = "blue", linewidth = 0.5) +
  facet_wrap(Domain ~ Diversity_Index, scales = "free") +
  labs(
    title = "Correlation Between Rarefied 1000 and 5000 Reads",
    x = "Diversity (Rarefied 1000 Reads)",
    y = "Diversity (Rarefied 5000 Reads)"
  ) +
  geom_text(
    data = cor_results_pearson,
    aes(x = x_pos, y = y_pos, label = label),
    hjust = 1.1, vjust = 1.1, size = 8 / .pt, color = "red", fontface = "bold", inherit.aes = FALSE
  ) +
  theme_bw(base_size = 8, base_family = "Arial") +
  theme(
    axis.text.x = element_text(size = 8, color = "black", angle = 90, hjust = 1),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title = element_text(size = 8, face = "bold"),
    strip.text = element_text(size = 8, face = "bold"),
    plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 8, family = "Arial"),
    legend.text = element_text(size = 8, family = "Arial"),
    legend.position = "bottom",
    legend.key.height = unit(0.3, "cm"),
    legend.key.width = unit(0.3, "cm"),
    panel.grid = element_blank()
  )

print(gg_plot)

ggsave(
  filename = paste0(wd_out, "/Correlation_Matrix_Rarefied_1000_5000_domains_v3_hexbin.svg"),
  plot = gg_plot,
  width = 3,
  height = 4,
  dpi = 300
)

### =========================================================================
### End of script
### =========================================================================