### This script analyses the richness on genus-level and compares it to the species level richness, particularly with respect to the contribution of Alphaproteobacteria and Cyanobacteriia taxa to prokryotic richness
### Author: Jonas Schiller

# Load dependencies
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggtext)
library(readxl)

# Data
dir.create("../../Data") #should already contain diversity_metrics.xlsx and samples_contextual.xlsx
dir.create("../../Data/Internal")
dir.create("../../Data/Data_generated") #should already contain alpha_div_genus.csv
dir.create("../../Data/Internal/FigS1_alpha_diversity_clade_level_GENUS_analysis_JS")

# Output
dir.create("../Output")
dir.create("../Output/FigS1_alpha_diversity_clade_level_GENUS_analysis_JS")

# rarefaction level
rarefying_level = 5000
diversity_metric = "Richness"

## read data
## GTDB genus diversity & GTDB species diversity
div_gtdb_genus = fread("../../Data/Data_generated/alpha_div_genus.csv") # or if script already run from: "../../Data/Internal/03.2_richness_rarefied_genus/alpha_div_genus.csv"
# alpha diversity
div_gtdb_species = as.data.table(read_excel("../../Data/diversity_metrics.xlsx", sheet = ifelse(rarefying_level == 5000,
                                                                                                     2, # sheet 2 contains alpha diversity metrics from a rarefaction level of 5,000 mOTU counts
                                                                                                     1) # sheet 1 contains alpha diversity metrics from a rarefaction level of 1,000 mOTU counts
))

# meta data
sheet2 = as.data.table(read_excel("../../Data/samples_contextual.xlsx", sheet = 2))
sheet4 = as.data.table(read_excel("../../Data/samples_contextual.xlsx", sheet = 4))
meta = merge(sheet2, sheet4, by = "biosample", all = TRUE)

# Generate dynamic column names
genus_col         = paste0(diversity_metric, "_", rarefying_level)
species_col       = paste0("prokaryotes_", diversity_metric, "_", rarefying_level)
cyano_genus_col   = paste0(diversity_metric, "_Cyano_", rarefying_level)
alpha_genus_col   = paste0(diversity_metric, "_Alphaproteo_", rarefying_level)
cyano_species_col = paste0("d__Bacteria_p__Cyanobacteriota_c__Cyanobacteriia_", diversity_metric, "_", rarefying_level)
alpha_species_col = paste0("d__Bacteria_p__Pseudomonadota_c__Alphaproteobacteria_", diversity_metric, "_", rarefying_level)

# Join genus and species diversity by biosample, then join with metadata
div_gtdb_all = left_join(div_gtdb_genus, div_gtdb_species, by = "biosample")
meta_div_gtdb = left_join(meta, div_gtdb_all, by = "biosample")

# Create 5° latitude bins
meta_div_gtdb = meta_div_gtdb %>%
  mutate(lat_bin = cut(lat, breaks = seq(-90, 90, by = 5), include.lowest = TRUE))

# Column names
genus_col        = paste0("Richness_", rarefying_level)
cyano_genus_col  = paste0("Richness_Cyano_", rarefying_level)
alpha_genus_col  = paste0("Richness_Alphaproteo_", rarefying_level)

species_col        = paste0("prokaryotes_", diversity_metric, "_", rarefying_level)
cyano_species_col  = paste0("class_Cyanobacteriia_", diversity_metric, "_", rarefying_level)
alpha_species_col  = paste0("class_Alphaproteobacteria_", diversity_metric, "_", rarefying_level)

#  Subset 
mld = meta_div_gtdb %>% filter(analysis_type == "MLD")
mesopelagic = meta_div_gtdb %>% filter(analysis_type == "mesopelagic")

#  Helper function 
print_mean_sd = function(data, cyano_col, alpha_col, total_col, label) {
  cat(label, "\n")
  cat("  Cyanobacteriia:       ", 
      round(mean(data[[cyano_col]], na.rm = TRUE), 1), " ± ",
      round(sd(data[[cyano_col]], na.rm = TRUE), 1), "\n")
  cat("  Alphaproteobacteria:  ", 
      round(mean(data[[alpha_col]], na.rm = TRUE), 1), " ± ",
      round(sd(data[[alpha_col]], na.rm = TRUE), 1), "\n")
  cat("  Total:                ", 
      round(mean(data[[total_col]], na.rm = TRUE), 1), " ± ",
      round(sd(data[[total_col]], na.rm = TRUE), 1), "\n\n")
}

## boxplots with overlaid points, faceted by analysis_type
div_col = paste0(diversity_metric, "_", rarefying_level)

# Filter data
df = meta_div_gtdb[meta_div_gtdb$analysis_type != "not_included", ]

# Bin latitudes into 5° bins
lat_breaks = seq(-90, 90, by = 5)
df$lat_bin = cut(df$lat, breaks = lat_breaks, include.lowest = TRUE, right = FALSE)
all_bins = levels(cut(lat_breaks[-1], breaks = lat_breaks, include.lowest = TRUE, right = FALSE))
df$lat_bin = factor(df$lat_bin, levels = all_bins)

# Format bin labels
clean_labels = gsub("\\[|\\)|\\]", "", levels(df$lat_bin))
clean_labels = gsub(",", " to ", clean_labels)
clean_labels = gsub("(-?\\d+)", "\\1°", clean_labels)
levels(df$lat_bin) = clean_labels

# Style labels: black for 0–40°, gray otherwise
label_colors = sapply(levels(df$lat_bin), function(label) {
  nums = as.numeric(gsub("°", "", unlist(strsplit(label, " to "))))
  if (all(nums >= -40 & nums < 40)) {
    paste0("<span style='color:black;'>", label, "</span>")
  } else {
    paste0("<span style='color:#696969;'>", label, "</span>")
  }
})
levels(df$lat_bin) = label_colors

# Custom x-axis label skipping every other tick
lat_bin_levels = levels(df$lat_bin)
x_labels_custom = setNames(
  object = ifelse(seq_along(lat_bin_levels) %% 2 == 1, lat_bin_levels, ""),
  nm = lat_bin_levels
)

# Add lat region info
df = df %>%
  mutate(abs_lat = abs(lat),
         lat_region = case_when(
           abs_lat <= 40 ~ "0–40°",
           abs_lat > 40  ~ "41–90°"
         ))

# Wilcoxon + fold-change stats
lat_region_stats = df %>%
  mutate(abs_lat = abs(lat),
         lat_region = case_when(
           abs_lat <= 40 ~ "0–40°",
           abs_lat > 40  ~ "41–90°"
         )) %>%
  group_by(analysis_type) %>%
  group_map(~{
    median_vals = .x %>%
      group_by(lat_region) %>%
      summarise(median_val = median(.data[[div_col]], na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = lat_region, values_from = median_val, names_prefix = "med_")
    
    test = wilcox.test(.x[[div_col]] ~ .x$lat_region)
    
    tibble(
      analysis_type = .y$analysis_type,
      label_region = paste0(
        "Wilcoxon:\n0–40° vs 41–90°\n",
        "FC = ", round(median_vals[["med_0–40°"]] / median_vals[["med_41–90°"]], 1),
        ", p = ", signif(test$p.value, 4)
      )
    )
  }) %>%
  bind_rows()

# Plot
p = ggplot(df, aes(x = lat_bin, y = .data[[div_col]])) +
  geom_boxplot(outlier.shape = NA, fill = "lightgray") +
  geom_jitter(width = 0.2, alpha = 0.4, size = 0.3) +
  facet_wrap(~analysis_type, nrow = 2, scales = "free_y") +
  geom_text(
    data = lat_region_stats,
    aes(x = levels(df$lat_bin)[1], y = Inf, label = label_region),
    inherit.aes = FALSE,
    hjust = 0, vjust = 1, size = 3.2, color = "blue", fontface = "bold"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_markdown(angle = 90, vjust = 0.5, hjust = 1),
    axis.ticks.y = element_line(),
    axis.ticks.x = element_line(),
    panel.grid = element_blank()
  ) +
  labs(
    x = "Latitude bin (5° step)",
    y = paste("Prokaryotic", diversity_metric),
    title = ""
  ) +
  scale_x_discrete(labels = x_labels_custom, drop = FALSE)
p
# Save
ggsave("../Output/FigS1_alpha_diversity_clade_level_GENUS_analysis_JS/LDG_genus.png", plot = p, width = 4, height = 5)
ggsave("../Output/FigS1_alpha_diversity_clade_level_GENUS_analysis_JS/LDG_genus.svg", plot = p, width = 4, height = 5)

## contribution of alphaproteobacteria and cyanobacteriia taxa to overall prokaryotic richness (pie chart)
# Filter for MLD samples only
df_mld = meta_div_gtdb %>% filter(analysis_type == "MLD")

# Helper function 
to_percent_pie = function(cyano, alpha, total) {
  cyano[is.na(cyano)] = 0
  alpha[is.na(alpha)] = 0
  total[is.na(total)] = 0
  
  other = total - cyano - alpha
  other = max(0, other)
  
  values = c(cyano, alpha, other)
  props = round(100 * values / total, 1)
  
  data.frame(
    Group = c("Cyanobacteriia", "Alphaproteobacteria", "Others"),
    Richness = values,
    Percent = props
  )
}

# Genus level column names 
genus_total_col = paste0("Richness_", rarefying_level)
genus_cyano_col = paste0("Richness_Cyano_", rarefying_level)
genus_alpha_col = paste0("Richness_Alphaproteo_", rarefying_level)

# Species level column names 
species_total_col = paste0("prokaryotes_", diversity_metric, "_", rarefying_level)
species_cyano_col = paste0("class_Cyanobacteriia_", diversity_metric, "_", rarefying_level)
species_alpha_col = paste0("class_Alphaproteobacteria_", diversity_metric, "_", rarefying_level)

# Genus level means 
genus_total = mean(df_mld[[genus_total_col]], na.rm = TRUE)
genus_cyano = mean(df_mld[[genus_cyano_col]], na.rm = TRUE)
genus_alpha = mean(df_mld[[genus_alpha_col]], na.rm = TRUE)

pie_genus = to_percent_pie(genus_cyano, genus_alpha, genus_total)
pie_genus$Level = "Genus"

# Species level means 
species_total = mean(df_mld[[species_total_col]], na.rm = TRUE)
species_cyano = mean(df_mld[[species_cyano_col]], na.rm = TRUE)
species_alpha = mean(df_mld[[species_alpha_col]], na.rm = TRUE)

pie_species = to_percent_pie(species_cyano, species_alpha, species_total)
pie_species$Level = "Species"

# Combine and plot 
pie_combined = rbind(pie_genus, pie_species)

ggplot(pie_combined, aes(x = "", y = Richness, fill = Group)) +
  geom_col(width = 1, color = "white") +
  geom_text(aes(label = paste0(Percent, "%")), 
            position = position_stack(vjust = 0.5), size = 5) +
  coord_polar(theta = "y") +
  facet_wrap(~Level, scales = "free") + 
  theme_void(base_size = 14) +
  labs(title = paste0("Mean Prokaryotic ", diversity_metric, 
                      " Contribution (MLD Samples, ", rarefying_level, "-rarefied)")) +
  theme(legend.title = element_blank())

# Save 
ggsave("../Output/FigS1_alpha_diversity_clade_level_GENUS_analysis_JS/fraction_richness_MLD.png", width = 7, height = 4.5)
ggsave("../Output/FigS1_alpha_diversity_clade_level_GENUS_analysis_JS/fraction_richness_MLD.svg", width = 7, height = 4.5)