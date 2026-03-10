### This script quantifies the fraction of taxonomy and diversity captured per sample.
### Author: Jonas Schiller

library(data.table)
library(ggplot2)
library(dplyr)
library(vegan)
library(tidyr)
library(readxl)

# create dir

# create output dir
## Plots
dir.create("./Output")
dir.create("./Output/FigS5_taxonomic_representation")
dir.create("../Data/External") # store samples_contextual.xlsx and mOTUsv4.profiles_subset_to_OMDB.tsv.gz here


## Data
dir.create("../Data/Internal/FigS5_taxonomic_representation")

## read

# meta
sheet2 = as.data.table(read_excel("../Data/External/samples_contextual.xlsx", sheet = 2))
sheet4 = as.data.table(read_excel("../Data/External/samples_contextual.xlsx", sheet = 4))
meta = merge(sheet2, sheet4, by = "biosample", all = TRUE)
# subset
meta = meta[meta$analysis_type != "not_included", ]

# motu
motu = fread("../Data/External/01_motus_meta/mOTUsv4.profiles_subset_to_OMDB.tsv.gz")
# columns are sample IDs --> clean them
colnames(motu) = make.names(colnames(motu))
colnames(motu) = gsub("_METAG","",colnames(motu))
colnames(motu) = sub("^[^_]+_", "",colnames(motu))
colnames(motu) = gsub("\\.", "-", colnames(motu))
## rm the motu col (store mOTU IDs as a vector)
motu_ids = motu$V1

# Unassigned fraction
## calculate unassigned fraction

# biosamples present in BOTH motu and meta
biosamples <- intersect(
  setdiff(names(motu), "V1"),
  meta$biosample
)

# extract unassigned counts
unassigned_counts <- as.numeric(
  motu[V1 == "mOTUv4.0_unassigned", ..biosamples]
)

# total counts per biosample
all_counts <- colSums(motu[, ..biosamples])

# build data.frame
unassigned <- data.frame(
  biosample = biosamples,
  unassigned_counts = unassigned_counts,
  all_counts = all_counts,
  unassigned_fraction = unassigned_counts / all_counts,
  row.names = NULL
)

## join to meta

meta_unassigned = dplyr::left_join(meta, unassigned)

# filter to keep only >5000 counts
meta_unassigned = meta_unassigned[meta_unassigned$rarefaction_included == "rarefaction_5000 & rarefaction_1000", ]


## plot

# lat bins
# Convert latitude to absolute values
meta_unassigned$abs_lat <- abs(meta_unassigned$lat)

# Create bins in 20° steps
meta_unassigned$lat_bin <- cut(
  meta_unassigned$abs_lat,
  breaks = c(0, 20, 40, 60, 80),
  labels = c("0-20", "21-40", "41-60", "61-80"),
  right = TRUE, include.lowest = TRUE
)

dodge_w <- 0.6

meta_unassigned <- meta_unassigned %>% filter(!is.na(lat_bin))
meta_na    <- meta_unassigned %>% filter(is.na(lat_bin))

median_labels <- meta_unassigned %>%
  group_by(analysis_type, lat_bin) %>%
  summarize(
    median_val = median(unassigned_fraction, na.rm = TRUE),
    y_pos = boxplot.stats(unassigned_fraction)$stats[5],
    .groups = "drop"
  )

ggplot(meta_unassigned, aes(x = analysis_type, y = unassigned_fraction, fill = lat_bin)) +
  geom_boxplot(
    data = meta_unassigned,
    outlier.shape = NA,
    position = position_dodge(width = dodge_w),
    color = "black"
  ) +
  geom_jitter(
    data = meta_unassigned,
    size = 0.3,
    alpha = 0.3,
    position = position_jitterdodge(
      jitter.width = 0.1,
      dodge.width = dodge_w
    )
  ) +
  geom_jitter(
    data = meta_na,
    aes(x = analysis_type, y = unassigned_fraction),
    inherit.aes = FALSE,
    size = 0.3,
    alpha = 0.3,
    color = "grey70",
    position = position_jitter(width = 0.12, height = 0)
  ) +
  geom_text(
    data = median_labels,
    aes(
      x = analysis_type,
      y = y_pos * 1.05,
      label = paste0(round(median_val * 100, 1), "%")
    ),
    position = position_dodge(width = dodge_w),
    angle = 90,
    vjust = 0.5,
    hjust = -0.5,
    size = 3.5,
    color = "black",
    show.legend = FALSE
  ) +
  theme_minimal() +
  labs(x = "Ocean layer", y = "Unassigned fraction", fill = "Latitude bin") +
  scale_fill_brewer(palette = "Spectral", na.value = "grey70", na.translate = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.15)))

ggsave("./Output/FigS5_taxonomic_representation/unassigned_fraction_layer.png", width = 10, height = 10, units = "cm")
ggsave("./Output/FigS5_taxonomic_representation/unassigned_fraction_layer.svg", width = 10, height = 10, units = "cm")

# Rarefaction

## calculate slopes

# biosamples present in BOTH motu and meta
biosamples <- intersect(
  setdiff(names(motu), "V1"),
  meta$biosample
)

# total counts per biosample (column sums, since columns are samples)
motus_sum <- colSums(motu[, ..biosamples])

# store results
slope_results <- data.table(
  biosample  = biosamples,
  motus_sum  = as.numeric(motus_sum),
  slope_5000 = NA_real_
)

# compute slope at 5000 for samples with enough counts
for (j in seq_along(biosamples)) {
  bs <- biosamples[j]
  sample_counts <- motu[[bs]]
  total_counts  <- slope_results[j, motus_sum]
  
  if (!is.na(total_counts) && total_counts >= 5000) {
    slope_results[j, slope_5000 := rareslope(sample_counts, sample = 5000)]
  }
}

## join

# left join slopes onto meta
meta_slopes <- meta %>%
  left_join(slope_results)

# filter to keep only >5000 counts
meta_slopes = meta_slopes[meta_slopes$rarefaction_included == "rarefaction_5000 & rarefaction_1000", ]


## plot

# Create bins in 20° steps
meta_slopes$abs_lat <- abs(meta_slopes$lat)
meta_slopes$lat_bin <- cut(
  meta_slopes$abs_lat,
  breaks = c(0, 20, 40, 60, 80),
  labels = c("0-20", "21-40", "41-60", "61-80"),
  right = TRUE, include.lowest = TRUE
)

# reshape to long format for plotting
df_long <- meta_slopes %>%
  dplyr::select(analysis_type, lat_bin, slope_5000) %>%
  pivot_longer(
    cols = starts_with("slope_"),
    names_to = "motu_count",
    values_to = "slope"
  ) %>%
  mutate(
    motu_count = gsub("slope_", "", motu_count),
    motu_count = paste0("Slope at ", motu_count, " reads")
  ) %>%
  filter(!is.na(lat_bin))

median_labels <- df_long %>%
  group_by(analysis_type, lat_bin, motu_count) %>%
  summarize(
    median_val = median(slope, na.rm = TRUE),
    y_pos = boxplot.stats(slope)$stats[5],
    .groups = "drop"
  )

# plot
dodge_w <- 0.6

ggplot(df_long, aes(x = analysis_type, y = slope, fill = lat_bin)) +
  geom_boxplot(
    outlier.shape = NA,
    position = position_dodge(width = dodge_w),
    color = "black"
  ) +
  geom_jitter(
    size = 0.3,
    alpha = 0.3,
    position = position_jitterdodge(
      jitter.width = 0.1,
      dodge.width = dodge_w
    )
  ) +
  geom_text(
    data = median_labels,
    aes(
      x = analysis_type,
      y = y_pos * 1.05,
      label = round(median_val, 3)
    ),
    position = position_dodge(width = dodge_w),
    angle = 90,
    vjust = 0.5,
    hjust = -0.4,
    size = 3.5,
    color = "black",
    show.legend = FALSE
  ) +
  facet_wrap(~ motu_count) +
  theme_minimal() +
  theme(
    strip.text = element_blank()   # <-- removes "Slope at 5000 reads"
  ) +
  labs(
    x = "Ocean layer",
    y = "Rarefaction slope (first derivative)",
    fill = "Latitude bin"
  ) +
  scale_fill_brewer(palette = "Spectral") +
  scale_color_brewer(palette = "Spectral") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.15)))
ggsave("./Output/FigS5_taxonomic_representation/slope_5000.png", width = 10, height = 10, units = "cm")
ggsave("./Output/FigS5_taxonomic_representation/slope_5000.svg", width = 10, height = 10, units = "cm")

### curves

set.seed(1)

meta_curves <- meta_slopes %>%
  filter(
    biosample %in% biosamples,
    !is.na(lat_bin),
    rarefaction_included == "rarefaction_5000 & rarefaction_1000"
  )

# keep all mesopelagic and a random 10% of MLD
picked_meso <- meta_curves %>%
  filter(analysis_type == "mesopelagic") %>%
  pull(biosample)

picked_mld_pool <- meta_curves %>%
  filter(analysis_type == "MLD") %>%
  pull(biosample)

n_mld <- max(1, ceiling(0.10 * length(picked_mld_pool)))
picked_mld <- if (length(picked_mld_pool) > 0) sample(picked_mld_pool, size = n_mld) else character(0)

picked <- unique(c(picked_meso, picked_mld))

picked_meta <- meta_curves %>%
  filter(biosample %in% picked) %>%
  select(biosample, analysis_type, lat_bin)

step_size <- 250
max_depth_global <- 5000

curve_list <- lapply(picked, function(bs) {
  sample_counts <- motu[[bs]]
  tot <- sum(sample_counts, na.rm = TRUE)
  if (is.na(tot) || tot < 2) return(NULL)
  
  max_depth <- min(tot, max_depth_global)
  
  if (max_depth < step_size) {
    depths <- max_depth
  } else {
    depths <- seq(from = step_size, to = max_depth, by = step_size)
    if (tail(depths, 1) != max_depth) depths <- c(depths, max_depth)
    depths <- unique(depths)
  }
  
  rich <- vegan::rarefy(sample_counts, sample = depths)
  
  tibble(
    biosample = bs,
    depth = depths,
    richness = as.numeric(rich)
  )
})

df_curve <- bind_rows(curve_list) %>%
  left_join(picked_meta, by = "biosample") %>%
  filter(!is.na(lat_bin), !is.na(richness))

ggplot(df_curve, aes(x = depth, y = richness, group = biosample, color = lat_bin)) +
  geom_line(alpha = 0.8, linewidth = 0.6) +
  facet_wrap(~analysis_type) +
  theme_minimal() +
  theme(text = element_text(size = 14))+
  labs(x = "Rarefied reads", y = "Expected richness", color = "Latitude bin") +
  scale_color_brewer(palette = "Spectral")
ggsave("./Output/FigS5_taxonomic_representation/curve_5000.png", width = 20, height = 10, units = "cm")
ggsave("./Output/FigS5_taxonomic_representation/curve_5000.svg", width = 20, height = 10, units = "cm")