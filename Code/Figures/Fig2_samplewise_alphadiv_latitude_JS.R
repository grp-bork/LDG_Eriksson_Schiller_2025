### This script plots sample-wise alpha diversity against latitude, either for epipelagic samples (defined by the mixed layer depth = MLD) or mesopelagic samples (200 - 1000 m depth). A statistical analysis is applied to determine significant differences in alpha diversity between (sub)-tropical regions (0 - 40° absolute latitude) and temperate to polar regions (>40° absolute latitude)
### Author: Jonas Schiller

# load dependencies
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggtext)
library(readxl)

# Metrics
analysis = "MLD" # "mesopelagic" "MLD"
rarefying_sample_size = 5000
alpha_diversity_metric = "Richness" # Richness Shannon Chao1

## create dirs
dir = paste0(analysis, "_",
             rarefying_sample_size, "_",
             alpha_diversity_metric)
# Data
dir.create("../../Data") #should already contain samples_contextual.xlsx and diversity_metrics.xlsx
dir.create("../../Data/Internal")
dir.create("../../Data/Data_generated")
dir.create("../../Data/Internal/Fig2_samplewise_alphadiv_latitude_JS")
dir.create(paste0("../../Data/Internal/Fig2_samplewise_alphadiv_latitude_JS/",dir))

# Output
dir.create("../Output")
dir.create("../Output/Fig2_samplewise_alphadiv_latitude_JS")
dir.create(paste0("../Output/Fig2_samplewise_alphadiv_latitude_JS/",dir))

## read data
# meta data
sheet2 = as.data.table(read_excel("../../Data/samples_contextual.xlsx", sheet = 2))
sheet4 = as.data.table(read_excel("../../Data/samples_contextual.xlsx", sheet = 4))
meta = merge(sheet2, sheet4, by = "biosample", all = TRUE)
# alpha diversity
alpha_diversity = as.data.table(read_excel("../../Data/diversity_metrics.xlsx", sheet = ifelse(rarefying_sample_size == 5000,
                                                                                                                          2, # sheet 2 contains alpha diversity metrics from a rarefaction level of 5,000 mOTU counts
                                                                                                                          1) # sheet 1 contains alpha diversity metrics from a rarefaction level of 1,000 mOTU counts
                                           ))

# join contextual and alpha diversity data
meta_alpha_div = left_join(meta,
                           alpha_diversity)

# filter for water layer (MLD or mesopelagic)
meta_alpha_div = meta_alpha_div %>%
  filter(analysis_type == analysis)


# create latitude bins
meta_alpha_div$lat_group_5 = cut(meta_alpha_div$lat,
                                 breaks = seq(-90, 90, by = 5),
                                 include.lowest = TRUE,
                                 right = FALSE,
                                 labels = paste(seq(-90, 85, by = 5), seq(-85, 90, by = 5), sep = "-"))


### Correlation of Bacterial and Archaeal alpha diversity with overall prokaryotic alpha diversity
## Bacteria with Prokaryotes
# Compute Pearson correlation
cor_test = cor.test(
  meta_alpha_div[[paste0("prokaryotes_", alpha_diversity_metric, "_", rarefying_sample_size)]],
  meta_alpha_div[[paste0("domain_Bacteria_", alpha_diversity_metric, "_", rarefying_sample_size)]],
  method = "pearson"
)
# Prepare correlation label
cor_label = paste0(
  "Pearson r = ", round(cor_test$estimate, 2),
  "\np = ", format.pval(cor_test$p.value, digits = 3, eps = 0.001)
)
# Create plot
# Calculate y position
y_vals = meta_alpha_div[[paste0("domain_Bacteria_", alpha_diversity_metric, "_", rarefying_sample_size)]]
y_pos = 0.95 * max(meta_alpha_div[[paste0("domain_Bacteria_", alpha_diversity_metric, "_", rarefying_sample_size)]], na.rm = TRUE)
p = ggplot(
  meta_alpha_div,
  aes_string(
    x = paste0("prokaryotes_", alpha_diversity_metric, "_", rarefying_sample_size),
    y = paste0("domain_Bacteria_", alpha_diversity_metric, "_", rarefying_sample_size)
  )
) +
  geom_point(size = 1) +
  geom_smooth(method = "lm", color = "blue", se = FALSE) +
  annotate("text", x = 50, y = y_pos, label = cor_label, hjust = 0, vjust = 1, size = 6, color = "red") + 
  theme_minimal() +
  theme(
    text = element_text(size = 18),
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    panel.grid = element_blank(),
    axis.ticks = element_line()
  ) +
  labs(
    title = "",
    x = paste0("Prokaryotic ", alpha_diversity_metric),
    y = paste0("Bacterial ", alpha_diversity_metric)
  )
p
ggsave(plot = p, file = paste0("../Output/Fig2_samplewise_alphadiv_latitude_JS/", dir, "/Prokaryotes_Bacteria.png"), width = 8, height = 8, units = "cm")
ggsave(plot = p, file = paste0("../Output/Fig2_samplewise_alphadiv_latitude_JS/", dir, "/Prokaryotes_Bacteria.svg"), width = 8, height = 8, units = "cm", device = "svg")

## Archaea with Prokaryotes
# Compute Pearson correlation
cor_test = cor.test(
  meta_alpha_div[[paste0("prokaryotes_", alpha_diversity_metric, "_", rarefying_sample_size)]],
  meta_alpha_div[[paste0("domain_Archaea_", alpha_diversity_metric, "_", rarefying_sample_size)]],
  method = "pearson"
)
# Prepare correlation label
cor_label = paste0(
  "Pearson r = ", round(cor_test$estimate, 2),
  "\np = ", format.pval(cor_test$p.value, digits = 3, eps = 0.001)
)
# Create plot
# Calculate y position
y_vals = meta_alpha_div[[paste0("domain_Bacteria_", alpha_diversity_metric, "_", rarefying_sample_size)]]
y_pos = 0.95 * max(meta_alpha_div[[paste0("domain_Archaea_", alpha_diversity_metric, "_", rarefying_sample_size)]], na.rm = TRUE)
p = ggplot(
  meta_alpha_div,
  aes_string(
    x = paste0("prokaryotes_", alpha_diversity_metric, "_", rarefying_sample_size),
    y = paste0("domain_Archaea_", alpha_diversity_metric, "_", rarefying_sample_size)
  )
) +
  geom_point(size = 1) +
  geom_smooth(method = "lm", color = "blue", se = FALSE) +
  annotate("text", x = 50, y = y_pos, label = cor_label, hjust = 0, vjust = 1, size = 6, color = "red") + 
  theme_minimal() +
  theme(
    text = element_text(size = 18),
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    axis.ticks = element_line(),
    panel.grid = element_blank()
  ) +
  labs(
    title = "",
    x = paste0("Prokaryotic ", alpha_diversity_metric),
    y = paste0("Archaeal ", alpha_diversity_metric)
  )
p
ggsave(plot = p, file = paste0("../Output/Fig2_samplewise_alphadiv_latitude_JS/", dir, "/Prokaryotes_Archaea.png"), width = 8, height = 8, units = "cm")
ggsave(plot = p, file = paste0("../Output/Fig2_samplewise_alphadiv_latitude_JS/", dir, "/Prokaryotes_Archaea.svg"), width = 8, height = 8, units = "cm", device = "svg")


### Wilcoxon groups
# preparation: long format and subsetting to domain level
df = meta_alpha_div[!is.na(meta_alpha_div$lat_group_5), ]
df_long = df %>%
  dplyr::select(lat_group_5, lat, 
                bacteria = paste0("domain_Bacteria_", alpha_diversity_metric, "_", rarefying_sample_size),
                archaea  = paste0("domain_Archaea_", alpha_diversity_metric, "_", rarefying_sample_size)) %>%
  pivot_longer(cols = c(bacteria, archaea), names_to = "domain", values_to = "richness") %>%
  mutate(
    domain = factor(ifelse(domain == "bacteria", "Bacteria", "Archaea"), levels = c("Bacteria", "Archaea")),
    abs_lat = abs(lat),
    lat_region = case_when(
      abs_lat <= 40 ~ "0-40°",
      abs_lat > 40 ~ "41-90°"
    )
  )

# Wilcoxon + Fold Change for 0–40° vs 41–90° absolute latitude
lat_region_stats = df_long %>%
  group_by(domain, lat_region) %>%
  summarise(
    median_richness = median(richness, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = lat_region, values_from = median_richness, names_prefix = "med_") %>%
  rowwise() %>%
  mutate(
    fold_change = `med_0-40°` / `med_41-90°`,
    p_value_region = wilcox.test(
      richness ~ lat_region,
      data = df_long[df_long$domain == domain, ]
    )$p.value,
    label_region = paste0(
      "Wilcoxon:\n0–40° vs 41–90°\n",
      "FC = ", round(fold_change, 2),
      ", p = ", format.pval(p_value_region, digits = 3, eps = 0.001)
    )
  ) %>%
  ungroup() %>%
  dplyr::select(domain, label_region)

# Create lat bins
lat_breaks = seq(-90, 90, by = 5)
df_long$lat_bin = cut(df_long$lat, 
                       breaks = lat_breaks, 
                       include.lowest = TRUE, 
                       right = FALSE)

# Set all possible bins
all_bins = levels(cut(lat_breaks[-1], breaks = lat_breaks, include.lowest = TRUE, right = FALSE))
df_long$lat_bin = factor(df_long$lat_bin, levels = all_bins)

# Clean up labels: "[x,y)" -> "x to y"
clean_labels = gsub("\\[|\\)|\\]", "", levels(df_long$lat_bin))
clean_labels = gsub(",", " to ", clean_labels)
clean_labels = gsub("(-?\\d+)", "\\1°", clean_labels)
levels(df_long$lat_bin) = clean_labels

# Determine color per label
label_colors = sapply(levels(df_long$lat_bin), function(label) {
  # Extract numeric bounds
  nums = as.numeric(gsub("°", "", unlist(strsplit(label, " to "))))
  if (all(nums >= -40 & nums < 40)) {
    paste0("<span style='color:black;'>", label, "</span>")
  } else {
    paste0("<span style='color:#696969;'>", label, "</span>")
  }
})
levels(df_long$lat_bin) = label_colors


# Plot
lat_bin_levels = levels(df_long$lat_bin)
x_labels_custom = setNames(
  object = ifelse(seq_along(lat_bin_levels) %% 2 == 1, lat_bin_levels, ""),
  nm = lat_bin_levels
)
p = ggplot(df_long, aes(x = lat_bin, y = richness)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = .1, width = 0.2, height = 0) +
  facet_wrap(~domain, nrow = 1, scales = "free_y") +
  geom_text(
    data = lat_region_stats, 
    aes(x = levels(df_long$lat_bin)[1], y = Inf, label = label_region),
    inherit.aes = FALSE,
    hjust = 0, vjust = 1, 
    size = 4, color = "blue", fontface = "bold"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 13),
    axis.text.x = ggtext::element_markdown(angle = 90, hjust = 1, vjust = 0.5),
    axis.ticks.y = element_line(),
    axis.ticks.x = element_line(),  
    panel.grid = element_blank()
  ) +
  labs(
    title = "",
    x = "Latitude bin (5° step)",
    y = "Richness"
  ) + 
  scale_x_discrete(labels = x_labels_custom, drop = FALSE)
p
# Save plots
ggsave(plot = p, file = paste0("../Output/Fig2_samplewise_alphadiv_latitude_JS/", dir, "/boxplot_latitude_diversity_Bacteria_Archaea_stat_wilcoxon.png"), width = 30, height = 8, units = "cm")
ggsave(plot = p, file = paste0("../Output/Fig2_samplewise_alphadiv_latitude_JS/", dir, "/boxplot_latitude_diversity_Bacteria_Archaea_stat_wilcoxon.svg"), width = 30, height = 8, units = "cm", device = "svg")