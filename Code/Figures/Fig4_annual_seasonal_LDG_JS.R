### This script compares and plots the similarity in class-level (+ Domain) LDGs (pair-wise Pearson correlation of modeled 1° latitude binned richness)
### Author: Jonas Schiller

# load dependencies
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(Hmisc)
library(ggdendro)
library(RColorBrewer)
library(ggtext)

# Metrics
analysis = "MLD" # "mesopelagic" "MLD"
rarefying_sample_size = 5000
alpha_diversity_metric = "Richness" # Richness Shannon Chao1

## create dirs
dir = paste0(analysis, "_",
             rarefying_sample_size, "_",
             alpha_diversity_metric)
# Data
dir.create("../Data")
dir.create("../Data/Internal")
dir.create("../Data/External") #should already contain df_latGradients_alphaDiv_annual.csv and mOTUsv4.0.gtdb.taxonomy.80mv.tsv and df_latGradients_alphaDiv_month.csv and Wilcox_pairwise_test_june_December.csv
dir.create("../Data/Internal/Fig4_annual_seasonal_LDG_JS")
dir.create(paste0("../Data/Internal/Fig4_annual_seasonal_LDG_JS/",dir))
# Output
dir.create("./Output")
dir.create("./Output/Fig4_annual_seasonal_LDG_JS")
dir.create(paste0("./Output/Fig4_annual_seasonal_LDG_JS/",dir))

# Read data
## modeled, binned (1°) diversity
diversity = fread("../Data/External/df_latGradients_alphaDiv_annual.csv")
## GTDB
gtdb = fread("../Data/External/mOTUsv4.0.gtdb.taxonomy.80mv.tsv")
# if a cell contains "Unknown", replace entire cell content with an NA
gtdb_na = gtdb[, lapply(.SD, function(x) ifelse(grepl("Unknown", x), NA, x))]

# Diversity per month
diversity_month = fread("../Data/External/df_latGradients_alphaDiv_month.csv")
# seasonal significance
wilcoxon = fread("../Data/External/Wilcox_pairwise_test_june_December.csv")

## Format
# Define the columns to keep (excluding unwanted ones)
cols_to_keep = setdiff(names(diversity), c("stat_type"))
# Select only those columns (keep all stat_type values for now)
div_temp = diversity[, ..cols_to_keep]
# Melt to long format
div_temp = melt(div_temp, id.vars = "y", variable.name = "clade", value.name = "Annual")

# taxonomy lookup tables
class_lookup = gtdb_na %>%
  distinct(class, phylum, domain)
phylum_lookup = gtdb_na %>%
  distinct(phylum, domain)
domain_lookup = gtdb_na %>%
  distinct(domain)

# identification of taxonomic level and extraction of taxon name
div_temp = div_temp %>%
  mutate(
    level = case_when(
      grepl("^class_", clade)  ~ "class",
      grepl("^phylum_", clade) ~ "phylum",
      grepl("^domain_", clade) ~ "domain",
      TRUE ~ NA_character_
    ),
    taxon_name = sub("^(class|phylum|domain)_", "", clade)
  ) %>%
  # Join for class-level entries
  left_join(class_lookup, by = c("taxon_name" = "class")) %>%
  # Join for phylum-level entries
  left_join(phylum_lookup, by = c("taxon_name" = "phylum"), suffix = c("", "_phylum")) %>%
  # Join for domain-level entries
  left_join(domain_lookup, by = c("taxon_name" = "domain"), suffix = c("", "_domain")) %>%
  mutate(
    clade = case_when(
      level == "class"  ~ paste0("domain__", domain, "_phylum__", phylum, "_class__", taxon_name),
      level == "phylum" ~ paste0("domain__", domain_phylum, "_phylum__", taxon_name),
      level == "domain" ~ paste0("domain__", taxon_name),
      TRUE ~ NA_character_
    )
  ) %>%
  select(-level, -taxon_name, -phylum, -domain, -domain_phylum)
# Join back with stat_type column to allow filtering
div_temp = cbind(stat_type = diversity$stat_type, div_temp)
# Final filtering for mean and renaming
div_long_mean = div_temp %>%
  filter(stat_type == "mean") %>%
  select(-stat_type)


### Pairwise correlation of Clade-level alpha diversities
make_dist_from_cor = function(C, method = c("abs", "signed")) {
  method = match.arg(method)
  C = as.matrix(C)
  diag(C) = 1
  
  if (method == "abs") {
    D = 1 - abs(C)
  } else {
    D = 1 - C
  }
  
  as.dist(D)
}

extract_base_name = function(x) {
  sub("_class__.*", "", x)
}

tax_level = function(x) {
  ifelse(grepl("_class__", x), "class", "phylum")
}

# wide matrix
diversity_wide = div_long_mean %>%
  select(y, clade, Annual) %>%
  pivot_wider(names_from = clade, values_from = Annual)
X = diversity_wide %>%
  select(-y) %>%
  as.data.frame()

# remove clades with too few observations
min_non_na = 10
keep_cols = names(X)[colSums(!is.na(X)) >= min_non_na]
X = X[, keep_cols, drop = FALSE]
# remove near-constant clades
nzv = sapply(X, function(v) sd(v, na.rm = TRUE))
X = X[, nzv > 0, drop = FALSE]

# correlation
cor_type = "pearson"
cor_results = rcorr(as.matrix(X), type = cor_type)
correlation_matrix = cor_results$r
p_value_matrix = cor_results$P

# BH/FDR correction on off-diagonal p-values
pvec = p_value_matrix[upper.tri(p_value_matrix)]
padj = rep(NA_real_, length(pvec))
padj[!is.na(pvec)] = p.adjust(pvec[!is.na(pvec)], method = "BH")
padj_matrix = matrix(
  NA_real_,
  nrow = nrow(p_value_matrix),
  ncol = ncol(p_value_matrix),
  dimnames = dimnames(p_value_matrix)
)
padj_matrix[upper.tri(padj_matrix)] = padj
padj_matrix = t(padj_matrix)
padj_matrix[upper.tri(padj_matrix)] = padj
diag(padj_matrix) = 0
alpha = 0.05
masked_correlation_matrix = correlation_matrix
masked_correlation_matrix[padj_matrix > alpha] = NA

# removing phyla if the finer taxonomic level is already representing it (corr > 0.9)
lab = colnames(correlation_matrix)
base = extract_base_name(lab)
lvl = tax_level(lab)

bases_with_both = intersect(base[lvl == "phylum"], base[lvl == "class"])
phylum_idx = which(lvl == "phylum" & base %in% bases_with_both)

to_remove = logical(length(lab))
for (i in phylum_idx) {
  same_base_class = which(lvl == "class" & base == base[i])
  if (any(correlation_matrix[i, same_base_class] > 0.9, na.rm = TRUE)) {
    to_remove[i] = TRUE
  }
}

labels_to_keep = lab[!to_remove]

correlation_matrix = correlation_matrix[labels_to_keep, labels_to_keep, drop = FALSE]
masked_correlation_matrix = masked_correlation_matrix[labels_to_keep, labels_to_keep, drop = FALSE]
padj_matrix = padj_matrix[labels_to_keep, labels_to_keep, drop = FALSE]


# clustering
dist_method = "signed"   # use "abs" to cluster by shape regardless of sign
distance_matrix = make_dist_from_cor(correlation_matrix, method = dist_method)
hc = hclust(distance_matrix, method = "average")
clades_order = hc$labels[hc$order]

# dendrogram
dendro = as.dendrogram(hc)
dendro_data = dendro_data(dendro)

# prepare labels
dendro_data$labels$label = gsub(
  pattern = paste0("_", alpha_diversity_metric, "_", rarefying_sample_size),
  replacement = "",
  x = dendro_data$labels$label
)
dendro_data$labels$label = gsub("domain__|phylum__|class__", "", dendro_data$labels$label)
dendro_data$labels$label = gsub("_", " ", dendro_data$labels$label)

# word frequencies
first_word_freq = table(sapply(strsplit(dendro_data$labels$label, " "), `[`, 1))
second_word_freq = table(sapply(strsplit(dendro_data$labels$label, " "), `[`, 2))
pal1 = brewer.pal(8, "Dark2")
pal2 = brewer.pal(12, "Paired")
dendro_data$labels = dendro_data$labels %>%
  mutate(
    first_word = sapply(strsplit(label, " "), `[`, 1),
    second_word = sapply(strsplit(label, " "), `[`, 2),
    remaining_text = sapply(strsplit(label, " "), function(x) paste(x[-c(1, 2)], collapse = " ")),
    first_word_color = ifelse(
      first_word_freq[first_word] > 1,
      pal1[(as.integer(as.factor(first_word)) - 1) %% length(pal1) + 1],
      "black"
    ),
    second_word_color = ifelse(
      !is.na(second_word) & second_word_freq[second_word] > 1,
      pal2[(as.integer(as.factor(second_word)) - 1) %% length(pal2) + 1],
      "black"
    ),
    html_label = paste0(
      "<span style='color:", first_word_color, "'>", first_word, "</span> ",
      ifelse(
        !is.na(second_word),
        paste0("<span style='color:", second_word_color, "'>", second_word, "</span> "),
        ""
      ),
      remaining_text
    )
  )

ggplot(segment(dendro_data)) +
  geom_segment(aes(x = y, y = x, xend = yend, yend = xend)) +
  geom_richtext(
    data = dendro_data$labels,
    aes(x = 0, y = x, label = paste0("<b>", html_label, "</b>")),
    hjust = 0,
    vjust = 0.5,
    angle = 0,
    size = 5,
    fill = NA,
    label.color = NA
  ) +
  labs(title = "", x = "", y = "") +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    plot.title = element_blank()
  ) +
  scale_y_reverse() +
  xlim(c(6, -15))
ggsave(paste0("./Output/Fig4_annual_seasonal_LDG_JS/",dir, "/pairwise_pearson_dendrogram.svg"),
       width = 21, height = 40, units = "cm")

cladenames = data.frame(clade = hc$labels[hc$order])

## Plot LDGs
### ...Individual with SD
# Set latitude bin
div_temp[, latitude_bin := y]

# Pivot clade and value columns to long format (now 'clade' holds taxon name, 'Annual' the value)
div_temp_long = melt(
  div_temp,
  id.vars = c("latitude_bin", "stat_type", "clade"),
  measure.vars = "Annual",
  variable.name = "metric",
  value.name = "value"
)

# Filter to include only clades in cladenames$clade
div_temp_long = div_temp_long[clade %in% cladenames$clade]

# Reshape to wide so each clade has columns for mean and sd
div_temp_wide = dcast(
  div_temp_long,
  latitude_bin + clade ~ stat_type,
  value.var = "value"
)

# Compute lower and upper bounds
div_temp_wide[, `:=`(
  lower = pmax(0, mean - sd),
  upper = mean + sd
)]

# Clean clade names
div_temp_wide[, clade := gsub(
  paste0("_", alpha_diversity_metric, "_", rarefying_sample_size, "_median"),
  "",
  clade
)]

# Convert clade to factor with defined order
div_temp_wide[, clade := factor(clade, levels = cladenames$clade)]

# Plotting
x_breaks = seq(-80, 80, by = 10)

p = ggplot(div_temp_wide[!is.na(mean)], aes(x = latitude_bin, y = mean, group = clade)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "darkgrey", alpha = 0.5) +
  geom_line(size = 1) +
  theme_minimal() +
  theme(
    text = element_text(size = 25),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 17),
    plot.title = element_text(hjust = 0.5),
    strip.text = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_line(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2),
    panel.spacing = unit(0.3, "lines"),
    strip.background = element_blank()
  ) +
  labs(x = "Latitude", y = "", title = "") +
  scale_x_continuous(breaks = x_breaks, limits = c(-80, 80)) +
  scale_y_continuous(limits = c(0, NA), expand = c(0, 0))+
  facet_wrap(~clade, ncol = 1, scales = "free_y")
p
ggsave(paste0("./Output/Fig4_annual_seasonal_LDG_JS/", dir, "/LDG_per_clade.png"),
       width = 15, height = 50, units = "cm")
ggsave(paste0("./Output/Fig4_annual_seasonal_LDG_JS/", dir, "/LDG_per_clade.svg"),
       width = 15, height = 50, units = "cm")





### Seasonal
## Create taxonomic level and stripped name
# for wilcoxon
wilcoxon = wilcoxon %>%
  mutate(
    level = case_when(
      grepl("^class_", clade)  ~ "class",
      grepl("^phylum_", clade) ~ "phylum",
      grepl("^domain_", clade) ~ "domain",
      TRUE ~ NA_character_
    ),
    taxon_name = sub("^(class|phylum|domain)_", "", clade)
  ) %>%
  # join taxonomy info
  left_join(class_lookup, by = c("taxon_name" = "class")) %>%
  left_join(phylum_lookup, by = c("taxon_name" = "phylum"), suffix = c("", "_phylum")) %>%
  left_join(domain_lookup, by = c("taxon_name" = "domain"), suffix = c("", "_domain")) %>%
  # Reconstruct clade names
  mutate(
    clade = case_when(
      level == "class"  ~ paste0("domain__", domain, "_phylum__", phylum, "_class__", taxon_name),
      level == "phylum" ~ paste0("domain__", domain_phylum, "_phylum__", taxon_name),
      level == "domain" ~ paste0("domain__", taxon_name),
      TRUE ~ NA_character_
    )
  ) %>%
  select(-level, -taxon_name, -phylum, -domain, -domain_phylum)
# for diversity_month
diversity_month = diversity_month %>%
  mutate(
    level = case_when(
      grepl("^class_", clade)  ~ "class",
      grepl("^phylum_", clade) ~ "phylum",
      grepl("^domain_", clade) ~ "domain",
      TRUE ~ NA_character_
    ),
    taxon_name = sub("^(class|phylum|domain)_", "", clade)
  ) %>%
  # Join to taxonomic info
  left_join(class_lookup, by = c("taxon_name" = "class")) %>%
  left_join(phylum_lookup, by = c("taxon_name" = "phylum"), suffix = c("", "_phylum")) %>%
  left_join(domain_lookup, by = c("taxon_name" = "domain"), suffix = c("", "_domain")) %>%
  # Create unified clade names
  mutate(
    clade = case_when(
      level == "class"  ~ paste0("domain__", domain, "_phylum__", phylum, "_class__", taxon_name),
      level == "phylum" ~ paste0("domain__", domain_phylum, "_phylum__", taxon_name),
      level == "domain" ~ paste0("domain__", taxon_name),
      TRUE ~ NA_character_
    )
  ) %>%
  select(-level, -taxon_name, -phylum, -domain, -domain_phylum)

## December/June LDGs per clade
# Define color scheme
color_months = c("June" = "deeppink", "December" = "black")
fill_months  = c("June" = scales::alpha("deeppink", 0.3), "December" = scales::alpha("darkgray", 0.5))
#Add latitude bin
diversity_month[, latitude_bin := y]
# Melt June & December columns
div_month_long = melt(
  diversity_month,
  id.vars = c("latitude_bin", "clade", "stat_type"),
  measure.vars = c("June", "December"),
  variable.name = "Month",
  value.name = "value"
)
div_month_long = div_month_long[div_month_long$clade %in% cladenames$clade]
# Cast wide (mean/sd)
div_month_wide = dcast(
  div_month_long,
  latitude_bin + clade + Month ~ stat_type,
  value.var = "value"
)
# Compute bounds
div_month_wide[, `:=`(
  lower = pmax(0, mean - sd),
  upper = mean + sd
)]
#Clean clade names
div_month_wide[, taxon := gsub(
  paste0("_", alpha_diversity_metric, "_", rarefying_sample_size, "_median"),
  "",
  clade
)]
div_month_wide[, taxon := factor(taxon, levels = clades_order)]
div_month_wide[, Month := factor(Month, levels = c("June", "December"))]
# Base plot
x_breaks = seq(-80, 80, by = 10)
p = ggplot(div_month_wide[!is.na(mean)], aes(x = latitude_bin, y = mean, group = interaction(taxon, Month))) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Month), alpha = 0.5) +
  geom_line(aes(color = Month), size = 1) +
  theme_minimal() +
  theme(
    text = element_text(size = 25),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 17),
    plot.title = element_text(hjust = 0.5),
    strip.text = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_line(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2),
    panel.spacing = unit(0.3, "lines"),
    strip.background = element_blank(),
    legend.position = "bottom"
  ) +
  labs(
    x = "Latitude", y = "", title = "",
    color = "Month", fill = "Month"
  ) +
  scale_color_manual(values = color_months, drop = FALSE) +
  scale_fill_manual(values = fill_months, drop = FALSE) +
  scale_x_continuous(breaks = x_breaks, limits = c(-80, 80)) +
  scale_y_continuous(limits = c(0, NA), expand = c(0, 0))

# Significant bins from Wilcoxon
# Clean clade names and factor levels
wilcoxon[, taxon := clade]
wilcoxon[, taxon := factor(taxon, levels = clades_order)]
# Filter: only bins where fraction significant > 0.5
sig_bins = wilcoxon[color_category == "above_50"]
# Max y for each clade
clade_max_y = div_month_wide[, .(max_y = max(upper, na.rm = TRUE)), by = taxon]
# Prepare blue bar segments
sig_line_data = merge(
  sig_bins[, .(latitude_bin = y, taxon)],
  clade_max_y,
  by = "taxon",
  all.x = TRUE
)[, `:=`(
  x = latitude_bin - 0.5,
  xend = latitude_bin + 0.5,
  y = max_y
)]
sig_line_data = sig_line_data[!is.na(taxon)]

## % significant labels
# Count total and "above_50" bins by taxon
total_bins = wilcoxon[, .N, by = taxon]
sig_counts = sig_bins[, .N, by = taxon]
# Compute percentage
sig_percent = merge(total_bins, sig_counts, by = "taxon", all.x = TRUE)[
  , .(taxon, perc_sig = round(100 * coalesce(N.y, 0) / N.x, 1))
]
# Set label positions (x fixed, y = 90% of upper range)
label_data = div_month_wide[, .(
  x = 72,
  y = 0.9 * max(upper, na.rm = TRUE)
), by = taxon]
label_data = merge(label_data, sig_percent, by = "taxon", all.x = TRUE)
label_data[is.na(perc_sig), perc_sig := 0]
label_data[, label := paste0(perc_sig, "%")]
label_data[, taxon := factor(taxon, levels = clades_order)]
# plot with free y axis
p_freeY = p + facet_wrap(~taxon, ncol = 1, scales = "free_y") +
  geom_segment(
    data = sig_line_data,
    aes(x = x, xend = xend, y = y, yend = y),
    inherit.aes = FALSE,
    color = "skyblue",
    alpha = 1,
    linewidth = 4
  ) +
  geom_label(
    data = label_data,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    fill = "gray",
    alpha = 0.5,
    color = "blue",
    size = 6,
    label.size = 0.2,
    label.r = unit(0.25, "lines")
  )
ggsave(paste0("./Output/Fig4_annual_seasonal_LDG_JS/", dir, "/LDG_month6_12_per_taxon_freeY_wilcoxon_withPerc.png"),
       plot = p_freeY, width = 15, height = 50, units = "cm")
ggsave(paste0("./Output/Fig4_annual_seasonal_LDG_JS/", dir, "/LDG_month6_12_per_taxon_freeY_wilcoxon_withPerc.svg"),
       plot = p_freeY, width = 15, height = 50, units = "cm")



## Mean richness difference Winter/Summer
# Copy input data and keep only 'mean'
dm = copy(diversity_month)
dm = dm[stat_type == "mean"]
# Assign hemisphere
dm[, hemisphere := ifelse(latitude_bin < 0, "SH", "NH")]
# Assign summer and winter values per hemisphere
dm[, summer := fifelse(
  hemisphere == "NH", June,
  fifelse(hemisphere == "SH", December, NA_real_)
)]
dm[, winter := fifelse(
  hemisphere == "NH", December,
  fifelse(hemisphere == "SH", June, NA_real_)
)]
# Remove rows with missing values
dm = dm[!is.na(summer) & !is.na(winter)]
# First calculate means and SDs
seasonal_stats = dm[, .(
  mean_summer = mean(summer, na.rm = TRUE),
  mean_winter = mean(winter, na.rm = TRUE),
  sd_summer = sd(summer, na.rm = TRUE),
  sd_winter = sd(winter, na.rm = TRUE)
), by = clade]
# Then calculate % difference and SD of % difference across samples
dm[, pct_diff_sample := 100 * (winter - summer) / summer]
sd_pct_diff = dm[, .(sd_pct_diff = sd(pct_diff_sample, na.rm = TRUE)), by = clade]
# Merge SD of % difference into main table
seasonal_stats = merge(seasonal_stats, sd_pct_diff, by = "clade")
# Now compute mean percent difference
seasonal_stats[, pct_diff := 100 * (mean_winter - mean_summer) / mean_summer]
# Clean taxon names
seasonal_stats[, taxon := gsub(
  paste0("_", alpha_diversity_metric, "_", rarefying_sample_size, "_median"),
  "",
  clade
)]
seasonal_stats = seasonal_stats[!is.na(taxon)]
seasonal_stats[, taxon := factor(taxon, levels = rev(clades_order))]
# Label
seasonal_stats[, label := paste0(round(pct_diff, 0), "±", round(sd_pct_diff, 0), "%")]
# Plot
ggplot(seasonal_stats[!is.na(taxon)], aes(x = "", y = taxon, fill = pct_diff)) +
  geom_tile(width = 0.9, height = 0.9, alpha = 0.5, color = "white") +
  geom_text(aes(label = label), color = "black", size = 1.2, fontface = "bold") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(x = NULL, y = NULL, title = NULL) +
  theme_void(base_size = 9) +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(-1, -1, -1, -1),
    legend.position = "none"
  )
ggsave(
  sprintf("./Output/Fig4_annual_seasonal_LDG_JS/%s/percent_change_mean_winter_vs_summer.png", dir),
  width = 0.7, height = 10, units = "cm"
)
ggsave(
  sprintf("./Output/Fig4_annual_seasonal_LDG_JS/%s/percent_change_mean_winter_vs_summer.svg", dir),
  width =0.7, height = 10, units = "cm"
)