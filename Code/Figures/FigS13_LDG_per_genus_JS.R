### This script compares the shapes of the annual species-level LDGs per genus (alpha diversity metric = richness) 
### Author: Jonas Schiller

library(data.table)
library(vegan)
library(dplyr)
library(readxl)
library(data.table)
library(tidyr)
library(Hmisc)
library(ggplot2)
library(ggdendro)
library(dynamicTreeCut)
library(dendextend)
library(scales)
library(ggrepel)

# create dir
dir.create("../Data/Internal/FigS13_LDG_per_genus")
dir.create("./Output/FigS13_LDG_per_genus")

# modeled annual species richness per genus (mean per 1° latitude bin)
mean_annual_dt = fread("../Data/External/mean_annual_ldg_genus.csv")

# mOTU taxonomic profile
motu = fread("../Data/External/01_motus_meta/mOTUsv4.profiles_subset_to_OMDB.tsv.gz")
colnames(motu) = make.names(colnames(motu))
colnames(motu) = gsub("_METAG","",colnames(motu))
colnames(motu) = sub("^[^_]+_", "",colnames(motu))
colnames(motu) <- gsub("\\.", "-", colnames(motu)) # replace points with dashes

# GTDB
gtdb = fread("../Data/External/01_motus_meta/mOTUsv4.0.gtdb.taxonomy.80mv.tsv")
# if a cell contains "Unknown", replace entire cell content with an NA
gtdb_na <- gtdb[, lapply(.SD, function(x) ifelse(grepl("Unknown", x), NA, x))]
# consider only mOTUs that are in the taxonomic profile
index = gtdb_na$motu %in% motu[[1]]
gtdb_na_select <- gtdb_na[index, ]

# meta data
sheet2 = as.data.table(read_excel("../Data/External/samples_contextual.xlsx", sheet = 2))
sheet4 = as.data.table(read_excel("../Data/External/samples_contextual.xlsx", sheet = 4))
meta = merge(sheet2, sheet4, by = "biosample", all = TRUE)


## relative abundance genus - full taxonomic profile
### summary rel. abundance 
# select relevant samples
samples_5000 <- unlist(meta[rarefaction_included == "rarefaction_5000 & rarefaction_1000" & analysis_type != "MLD", biosample])
motu_counts  <- motu[, ..samples_5000]
# attach motu ids
motu[, motu_id := V1]

# join GTDB taxonomy (all ranks)
motu_tax <- merge(
  motu,
  gtdb_na_select[, .(motu, domain, phylum, class, order, family, genus)],
  by.x = "motu_id",
  by.y = "motu",
  all.x = TRUE
)
# identify sample columns
sample_cols <- setdiff(
  names(motu_tax),
  c("motu_id", "V1", "domain", "phylum", "class", "order", "family", "genus")
)
## Average relative abundance of taxa (for panel A in Fig. S13)
# function to summarize one taxonomic rank
summarize_rank <- function(dt, rank_col, sample_cols) {
  
  # split assigned and unassigned
  assigned   <- dt[!is.na(get(rank_col))]
  unassigned <- dt[is.na(get(rank_col))]
  
  # sum counts per taxon
  tax_counts <- assigned[
    , lapply(.SD, sum),
    by = rank_col,
    .SDcols = sample_cols
  ]
  
  # add unassigned as its own taxon
  if (nrow(unassigned) > 0) {
    unassigned_counts <- unassigned[, lapply(.SD, sum), .SDcols = sample_cols]
    unassigned_counts[, (rank_col) := "unassigned"]
    tax_counts <- rbind(tax_counts, unassigned_counts, fill = TRUE)
  }
  
  # convert counts to relative abundance per sample
  rel_abund <- copy(tax_counts)
  rel_abund[, (sample_cols) := lapply(.SD, function(x) x / sum(x)), .SDcols = sample_cols]
  
  # long format for stats
  rel_long <- melt(
    rel_abund,
    id.vars = rank_col,
    measure.vars = sample_cols,
    value.name = "rel_abundance"
  )
  
  # compute summary statistics
  summary_dt <- rel_long[
    ,
    .(
      median_rel_abundance = median(rel_abundance, na.rm = TRUE),
      mean_rel_abundance   = mean(rel_abundance, na.rm = TRUE),
      SD                   = sd(rel_abundance, na.rm = TRUE)
    ),
    by = rank_col
  ]
  
  # standardize column names
  setnames(summary_dt, rank_col, "taxon")
  summary_dt[, taxonomic_rank := rank_col]
  
  summary_dt[]
}
# taxonomic ranks to process
taxonomic_ranks <- c("domain", "phylum", "class", "order", "family", "genus")
# compute summaries for all ranks
rel_abundance_summary <- rbindlist(
  lapply(taxonomic_ranks, summarize_rank, dt = motu_tax, sample_cols = sample_cols),
  fill = TRUE
)
# order output for readability
setorder(rel_abundance_summary, taxonomic_rank, taxon)



## Dendrogram representing similarity in annual LDG shape (species richness per genus)
# function for dist from correlation
make_dist_from_cor <- function(C) {
  C <- as.matrix(C)
  diag(C) <- 1
  
  # signed distance
  D <- 1 - C
  
  as.dist(D)
}

# wide matrix: rows = y, cols = genus, values = mean_annual
diversity_wide <- mean_annual_dt %>%
  select(y, genus, mean_annual) %>%
  pivot_wider(names_from = genus, values_from = mean_annual)

# numeric matrix for correlation/clustering
X <- diversity_wide %>%
  select(-y) %>%
  as.data.frame()

# correlations across genera
cor_type <- "pearson"
cor_results <- rcorr(as.matrix(X), type = cor_type)

C <- cor_results$r
P <- cor_results$P

# clean correlation matrix for distance conversion
C_for_dist <- C
C_for_dist[is.na(C_for_dist)] <- 0
diag(C_for_dist) <- 1

# correlation distance and hierarchical clustering
D  <- make_dist_from_cor(C_for_dist)
linkage <- "average"
hc <- hclust(D, method = linkage)

# dendrogram objects for ggplot/dendextend workflows
dend <- as.dendrogram(hc)
dd <- ggdendro::dendro_data(dend)

# keep genus naming consistent end-to-end
genus_names <- colnames(X)
stopifnot(length(genus_names) == ncol(X))

# enforce dimnames so labels propagate through dist/hclust
C_for_dist <- as.matrix(C_for_dist)
rownames(C_for_dist) <- colnames(C_for_dist) <- genus_names
diag(C_for_dist) <- 1

# check dist carries labels from matrix dimnames
D <- make_dist_from_cor(C_for_dist)
stopifnot(identical(attr(D, "Labels"), genus_names))

# check hclust labels match genus names
hc <- hclust(D, method = linkage)
stopifnot(identical(hc$labels, genus_names))

# full distance matrix for methods that need matrix input
distM <- as.matrix(D)
stopifnot(identical(rownames(distM), genus_names))
stopifnot(identical(colnames(distM), genus_names))

# cut tree
clusters_raw = cutree(hc, k = 7)

# attach genus names in the same order as hclust labels
clusters <- setNames(as.integer(clusters_raw), hc$labels)
stopifnot(length(clusters) == length(genus_names))
stopifnot(all(names(clusters) %in% genus_names))

# leaf x positions + cluster membership
leaf_df <- dd$labels %>%
  as_tibble() %>%
  mutate(cluster = clusters[label]) %>%
  filter(!is.na(cluster))

# relabel clusters clockwise for consistent plotting
cluster_map <- leaf_df %>%
  group_by(cluster) %>%
  summarise(x_mean = mean(x), .groups = "drop") %>%
  arrange(x_mean) %>%
  mutate(cluster_new = row_number()) %>%
  { setNames(.$cluster_new, as.character(.$cluster)) }

leaf_df <- leaf_df %>%
  mutate(cluster_new = unname(cluster_map[as.character(cluster)]))

# color branches by cluster subtrees
K <- length(unique(leaf_df$cluster_new))
cols <- setNames(hue_pal()(K), as.character(1:K))

# map original cut clusters to the clockwise-renamed ids
clusters_new <- setNames(
  unname(cluster_map[as.character(clusters)]),
  names(clusters)
)

# apply colors to dendrogram branches using the renamed cluster ids
dend_col <- color_branches(
  dend,
  clusters = clusters_new[labels(dend)],
  col = cols
)

dd_col <- as.ggdend(dend_col)

# segment table used by ggplot
seg_df <- dd_col$segments %>%
  as_tibble() %>%
  mutate(cluster_new = as.integer(names(cols)[match(col, cols)]))

# place one label per cluster just outside the tips
outer_r <- max(c(seg_df$y, seg_df$yend), na.rm = TRUE)
tip_offset <- 0.05 * outer_r

cluster_labels <- leaf_df %>%
  group_by(cluster_new) %>%
  summarise(x = mean(x), .groups = "drop") %>%
  mutate(
    y = -tip_offset,
    label = letters[as.integer(cluster_new)],
    angle = 90 - 360 * (x - min(dd$labels$x)) / (max(dd$labels$x) - min(dd$labels$x)),
    hjust = ifelse(angle < -90, 1, 0),
    angle = ifelse(angle < -90, angle + 180, angle)
  )

# pick specific tips to annotate and strip the genus_ prefix
tip_labels_df <- dd$labels %>%
  as_tibble() %>%
  filter(label %in% c("genus_Prochlorococcus_A", "genus_Pelagibacter")) %>%
  mutate(label_clean = sub("^genus_", "", label))

# plot 1: colored circular dendrogram + cluster ids + two labeled tips
p_dend <- ggplot() +
  geom_segment(
    data = seg_df,
    aes(x = x, y = y, xend = xend, yend = yend, color = factor(cluster_new)),
    linewidth = 0.4
  ) +
  geom_text(
    data = cluster_labels,
    aes(x = x, y = y, label = label, color = factor(cluster_new), angle = angle, hjust = hjust),
    size = 6,
    fontface = "bold"
  ) +
  coord_polar(theta = "x", clip = "off") +
  scale_y_reverse(expand = expansion(mult = c(0.18, 0.02))) +
  scale_color_manual(values = cols) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = margin(20, 20, 20, 20),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

p_dend <- p_dend +
  geom_text(
    data = tip_labels_df,
    aes(x = x, y = 0, label = label_clean),
    color = "black",
    size = 5,
    hjust = 0,
    vjust = 0.5
  )

print(p_dend)

# write dendrogram as png and svg
ggsave("./Output/FigS13_LDG_per_genus/dendrogram_genus.png", plot = p_dend, width = 30, height = 20, units = "cm", dpi = 300)
ggsave("./Output/FigS13_LDG_per_genus/dendrogram_genus.svg", plot = p_dend, width = 30, height = 20, units = "cm", bg = "white")

# apply renamed cluster ids to genus table
cluster_df <- tibble(
  genus = names(clusters),
  cluster = unname(cluster_map[as.character(unname(clusters))])
)


## Per Genus: species richness LDG
# plot 2: scaled LDGs per renamed cluster
ldg_clustered <- mean_annual_dt %>%
  inner_join(cluster_df, by = "genus") %>%
  group_by(genus) %>%
  mutate(richness_scaled = as.numeric(scale(mean_annual))) %>%
  ungroup()

# attach the same letter labels used in the dendrogram (a, b, c, ...)
ldg_clustered <- ldg_clustered %>%
  mutate(cluster_label = letters[as.integer(cluster)])

# count genera per cluster for facet annotations
n_df <- ldg_clustered %>%
  distinct(cluster_label, genus) %>%
  count(cluster_label, name = "n")

# join relative abundance table
ldg_clustered$taxon = gsub("genus_", "", ldg_clustered$genus)
ldg_clustered = dplyr::left_join(ldg_clustered, rel_abundance_summary)

# one row per facet for placing axis titles as text
axis_lab_df <- ldg_clustered %>%
  distinct(cluster_label)

# color facet strip titles using the same palette as the dendrogram clusters
cluster_cols <- cols
facet_labeller <- labeller(
  cluster_label = function(x) {
    idx <- match(x, letters)
    sprintf("<span style='color:%s'><b>%s</b></span>", cluster_cols[as.character(idx)], x)
  }
)

p_ldg <- ggplot(ldg_clustered, aes(x = y, y = richness_scaled, group = genus, color = mean_rel_abundance)) +
  geom_line(alpha = 0.6) +
  facet_wrap(~ cluster_label, scales = "free", labeller = facet_labeller, nrow = 1) +
  geom_text(
    data = n_df,
    aes(label = paste0("n = ", n)),
    x = -Inf, y = Inf,
    hjust = -0.1, vjust = 1.2,
    inherit.aes = FALSE,
    size = 4
  ) +
  scale_x_continuous(
    breaks = seq(-75, 75, by = 25),
    limits = c(-80, 80)
  ) +
  scale_color_viridis_c(
    name = "Mean relative\nabundance",
    trans = "log10",
    na.value = "grey80",
    guide = guide_colorbar(direction = "horizontal", title.position = "top", label.position = "bottom")
  ) +
  labs(x = NULL, y = NULL) +
  geom_text(
    data = axis_lab_df,
    aes(x = -Inf, y = 0.5, label = "Richness (scaled)"),
    inherit.aes = FALSE,
    angle = 90,
    hjust = 0.5,
    vjust = -3
  ) +
  geom_text(
    data = axis_lab_df,
    aes(x = 0, y = -Inf, label = "Latitude (°)"),
    inherit.aes = FALSE,
    vjust = 4
  ) +
  theme_minimal() +
  ylim(c(-3, 4)) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background  = element_rect(fill = "transparent", color = NA),
    strip.background  = element_rect(fill = "transparent", color = NA),
    axis.ticks = element_line(),
    axis.text  = element_text(size = 12),
    axis.title = element_text(size = 11),
    strip.text = ggtext::element_markdown(face = "bold", size = 14),
    plot.margin = margin(14, 14, 32, 28),
    panel.spacing = unit(1.2, "lines"),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.text = element_text(angle = 45, hjust = 1)
  ) +
  coord_cartesian(clip = "off")

print(p_ldg)

# write LDG plot as png and svg
ggsave("./Output/FigS13_LDG_per_genus/ldg_genus.png", plot = p_ldg, width = 50, height = 10, units = "cm", dpi = 300)
ggsave("./Output/FigS13_LDG_per_genus/ldg_genus.svg", plot = p_ldg, width = 50, height = 10, units = "cm")

## Distribution of mean rel. abundance
# one dot per genus (mean rel. abundance), keep cluster letters and colors consistent
abund_genus <- ldg_clustered %>%
  distinct(cluster_label, genus, mean_rel_abundance) %>%
  mutate(x = 1)

# color map for a..z using your existing 'cols'
cluster_cols <- cols
col_map <- setNames(cluster_cols[as.character(seq_along(letters))], letters)

# fixed y position for the "x-axis letter" in log space (same height in every facet)
min_y <- min(abund_genus$mean_rel_abundance, na.rm = TRUE)
y_axis_label <- 0

axis_df <- abund_genus %>%
  distinct(cluster_label) %>%
  mutate(x = 1, y = y_axis_label)

# label only these two genera and strip the genus_ prefix
label_df <- abund_genus %>%
  filter(genus %in% c("genus_Prochlorococcus_A", "genus_Pelagibacter")) %>%
  mutate(label = sub("^genus_", "", genus))

# deterministic label positions per facet (keep left offset like before)
label_pos <- label_df %>%
  mutate(
    x_lab = 0.85,
    y_lab = mean_rel_abundance * 1.4
  )

p_abund <- ggplot(abund_genus, aes(x = x, y = mean_rel_abundance)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, linewidth = 0.4, color = "grey30") +
  geom_point(size = 0.5, alpha = 0.9) +
  facet_wrap(~ cluster_label, nrow = 1) +
  
  scale_y_continuous(
    trans = "log10",
    expand = expansion(mult = c(0.45, 0.25))  # room for the letters below
  ) +
  
  # ONE tick per facet (no manual tick drawing)
  scale_x_continuous(
    limits = c(0.7, 1.3),
    breaks = 1,
    labels = NULL,
    minor_breaks = NULL
  ) +
  
  labs(
    title = "Mean relative abundance of taxa within Latitudinal Diversity Gradient groups",
    subtitle = "Species-level richness per genus",
    x = NULL,
    y = "Mean relative abundance (log10)"
  ) +
  
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background  = element_rect(fill = "transparent", color = NA),
    
    strip.text = element_blank(),
    strip.background = element_blank(),
    
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    
    axis.text.y  = element_text(size = 12),
    axis.title   = element_text(size = 11),
    
    plot.title    = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 12),
    
    plot.margin = margin(14, 14, 32, 28),
    panel.spacing = unit(1.2, "lines"),
    legend.position = "none"
  ) +
  coord_cartesian(clip = "off") +
  
  # one colored letter per facet at a fixed "x-axis" height
  geom_text(
    data = axis_df,
    aes(x = x, y = y, label = cluster_label, color = cluster_label),
    fontface = "bold",
    size = 7,
    vjust = 0,          # anchor text downward
    inherit.aes = FALSE
  )+
  scale_color_manual(values = col_map, guide = "none") +
  
  # leader line: label → original point (unchanged offset)
  geom_segment(
    data = label_pos,
    aes(
      x = x_lab,
      y = y_lab,
      xend = 1,
      yend = mean_rel_abundance
    ),
    color = "blue",
    linewidth = 0.35,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = label_pos,
    aes(x = x_lab, y = y_lab, label = label),
    color = "blue",
    size = 3.5,
    hjust = 1,
    inherit.aes = FALSE
  )
print(p_abund)

ggsave("./Output/FigS13_LDG_per_genus/abundance_boxplot_by_cluster.png",
       plot = p_abund, width = 23, height = 10, units = "cm", dpi = 300)
ggsave("./Output/FigS13_LDG_per_genus/abundance_boxplot_by_cluster.svg",
       plot = p_abund, width = 23, height = 10, units = "cm")