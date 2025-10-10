### This script randomly subsamples the surface samples to the count of mesopelagic samples to test the robustness of the observed surface LDG
### Author: Jonas Schiller

# load dependencies
library(data.table)
library(ggplot2)
library(readxl)
library(dplyr)
library(tidyr)
library(purrr)

# Metrics
rarefying_sample_size = 5000
alpha_diversity_metric = "Richness" # Richness Shannon Chao1

## create dirs
dir = paste0(rarefying_sample_size, "_",
             alpha_diversity_metric)

# Data
dir.create("../../Data")#should already contain samples_contextual.xlsx and diversity_metrics.xlsx
dir.create("../../Data/Internal")
dir.create("../../Data/Data_generated") 
dir.create("../../Data/Internal/FigS2_LDG_surface_mesopelagic_stat_JS")
dir.create(paste0("../../Data/Internal/FigS2_LDG_surface_mesopelagic_stat_JS/",dir))

# Output
dir.create("../Output")
dir.create("../Output/FigS2_LDG_surface_mesopelagic_stat_JS")
dir.create(paste0("../Output/FigS2_LDG_surface_mesopelagic_stat_JS/",dir))

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

### Prepare data
# create latitude bins
meta_alpha_div$lat_group_5 = cut(meta_alpha_div$lat,
                                       breaks = seq(-90, 90, by = 5),
                                       include.lowest = TRUE,
                                       right = FALSE,
                                       labels = paste(seq(-90, 85, by = 5), seq(-85, 90, by = 5), sep = "-"))

# water layer subset
meta_alpha_div_mld = meta_alpha_div[meta_alpha_div$analysis_type == "MLD",]
meta_alpha_div_mesopelagic = meta_alpha_div[meta_alpha_div$analysis_type == "mesopelagic",]

# subset surface samples randomly to count of mesopelagic samples (multiple iterations) --> calculte Wilcoxon test (low vs high latitudes)
# Setup
set.seed(123)
n_mesopelagic = sum(!is.na(meta_alpha_div_mesopelagic$prokaryotes_Richness_5000))
n_iterations = 1000

# Fold change + Wilcoxon function
calculate_fc_wilcox = function(sampled_df) {
  df_long = sampled_df %>%
    select(lat, 
           bacteria = paste0("domain_Bacteria_", alpha_diversity_metric, "_", rarefying_sample_size),
           archaea  = paste0("domain_Archaea_", alpha_diversity_metric, "_", rarefying_sample_size)) %>%
    pivot_longer(cols = c(bacteria, archaea), names_to = "domain", values_to = "richness") %>%
    mutate(
      domain = factor(ifelse(domain == "bacteria", "Bacteria", "Archaea"), 
                      levels = c("Bacteria", "Archaea")),
      abs_lat = abs(lat),
      lat_region = case_when(
        abs_lat <= 40 ~ "0-40°",
        abs_lat > 40 ~ "41-90°"
      )
    )
  
  stats_df = df_long %>%
    group_by(domain, lat_region) %>%
    summarise(median_richness = median(richness, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = lat_region, values_from = median_richness, names_prefix = "med_") %>%
    rowwise() %>%
    mutate(
      fold_change = `med_0-40°` / `med_41-90°`,
      p_value = wilcox.test(
        richness ~ lat_region,
        data = df_long[df_long$domain == domain, ]
      )$p.value
    ) %>%
    ungroup() %>%
    select(domain, fold_change, p_value)
  
  return(stats_df)
}

# Run bootstrap iterations
results = purrr::map_dfr(1:n_iterations, ~{
  sampled = meta_alpha_div_mld %>%
    dplyr::slice_sample(n = n_mesopelagic)
  calculate_fc_wilcox(sampled) %>%
    mutate(iteration = .x)
})

# Plot fold change distribution as boxplot
p1_box = ggplot(results, aes(x = domain, y = fold_change)) +
  geom_boxplot(fill = "steelblue", color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.3, color = "darkblue", size = .1) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  theme_minimal() +
  theme(
    text = element_text(size = 15),
    panel.grid = element_blank(),
    axis.ticks = element_line()
  ) +
  labs(
    title = alpha_diversity_metric,
    x = "Domain",
    y = "Fold Change"
  )

ggsave(paste0("../Output/FigS2_LDG_surface_mesopelagic_stat_JS/", dir, "/foldchange_boxplot.png"), plot = p1_box,
       width = 8, height = 8, units = "cm")

ggsave(paste0("../Output/FigS2_LDG_surface_mesopelagic_stat_JS/", dir, "/foldchange_boxplot.svg"), plot = p1_box,
       width = 8, height = 8, units = "cm", device = "svg")


# Plot Wilcoxon p-value distribution as boxplot
p2_box = ggplot(results, aes(x = domain, y = p_value)) +
  geom_boxplot(fill = "salmon", color = "black", outlier.shape = NA) +
  geom_jitter(aes(color = p_value < 0.05), width = 0.2, alpha = 0.3, size = 0.1) +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "darkred")) +
  scale_y_log10() +
  theme_minimal() +
  theme(
    text = element_text(size = 15),
    panel.grid = element_blank(),
    axis.ticks = element_line(),
    legend.position = "none"
  ) +
  labs(
    title = alpha_diversity_metric,
    x = "Domain",
    y = "p-value"
  )

ggsave(paste0("../Output/FigS2_LDG_surface_mesopelagic_stat_JS/", dir, "/p_wilcoxon_boxplot.png"), plot = p2_box,
       width = 8, height = 8, units = "cm")

ggsave(paste0("../Output/FigS2_LDG_surface_mesopelagic_stat_JS/", dir, "/p_wilcoxon_boxplot.svg"), plot = p2_box,
       width = 8, height = 8, units = "cm", device = "svg")