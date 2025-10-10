### This script compares the similarity in class-level LDGs (pair-wise Pearson correlation of modeled 1° latitude binned richness)
### Author: Jonas Schiller

# load dependencies
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggdendro)
library(dendextend)
library(svglite)

# Metrics
rarefying_sample_size = 5000
alpha_diversity_metric = "Richness" # Richness Shannon Chao1
average_metric = "mean"

## create dirs
dir = paste0(rarefying_sample_size, "_",
             alpha_diversity_metric, "_",
             average_metric)

# Data
dir.create("../../Data")
dir.create("../../Data/Internal")
dir.create("../../Data/Data_generated") #should already contain df_latGradients_alphaDiv_annual_class.csv and dist_bac.csv and dist_arc.csv
dir.create("../../Data/Internal/FigS8_tanglegram_LDG_phylogeny_JS")
dir.create(paste0("../../Data/Internal/FigS8_tanglegram_LDG_phylogeny_JS/",dir))

# Output
dir.create("../Output")
dir.create("../Output/FigS8_tanglegram_LDG_phylogeny_JS")
dir.create(paste0("../Output/FigS8_tanglegram_LDG_phylogeny_JS/",dir))

# Read data
# alpha diversity
diversity = fread("../../Data/Data_generated/df_latGradients_alphaDiv_annual.csv")
# phylogenetic distances
dist_bac = fread("../../Data/Data_generated/dist_bac.csv")
dist_arc = fread("../../Data/Data_generated/dist_arc.csv")

# clades_order (from tanglegram)
clades_order <- c(
  "class_Fibrobacteria",
  "class_Nitrososphaeria",
  "class_Acidimicrobiia",
  "class_Cyanobacteriia",
  "class_Alphaproteobacteria",
  "class_Chlamydiia",
  "class_UBA796",
  "class_UBA4151",
  "class_XYA12-FULL-58-9",
  "class_SAR324",
  "class_Gammaproteobacteria",
  "class_Marinisomatia",
  "class_UBA1144",
  "class_Bacteroidia",
  "class_Poseidoniia"
) 

# ensure required columns are selected
diversity = diversity %>% dplyr::select(c("y", "stat_type",all_of(clades_order)))

# Adjust names: Remove everything before and including the last "__"
dist_bac$Var1 = sub(".*__", "", dist_bac$Var1)
dist_bac$Var2 = sub(".*__", "", dist_bac$Var2)
dist_arc$Var1 = sub(".*__", "", dist_arc$Var1)
dist_arc$Var2 = sub(".*__", "", dist_arc$Var2)

# Convert pairwise distance tables into separate matrices
dist_matrix_bac <- dcast(dist_bac, Var1 ~ Var2, value.var = "value", fill = 0)
rownames(dist_matrix_bac) <- dist_matrix_bac$Var1
dist_matrix_bac <- as.matrix(dist_matrix_bac[, -1])

dist_matrix_arc <- dcast(dist_arc, Var1 ~ Var2, value.var = "value", fill = 0)
rownames(dist_matrix_arc) <- dist_matrix_arc$Var1
dist_matrix_arc <- as.matrix(dist_matrix_arc[, -1])

# Ensure symmetry
dist_matrix_bac[lower.tri(dist_matrix_bac)] <- t(dist_matrix_bac)[lower.tri(dist_matrix_bac)]
dist_matrix_arc[lower.tri(dist_matrix_arc)] <- t(dist_matrix_arc)[lower.tri(dist_matrix_arc)]

# Convert to distance objects
dist_obj_bac <- as.dist(dist_matrix_bac)
dist_obj_arc <- as.dist(dist_matrix_arc)

# Perform hierarchical clustering separately
hc_bac <- hclust(dist_obj_bac, method = "ward.D2")
hc_arc <- hclust(dist_obj_arc, method = "ward.D2")

## connection of the 2 trees
# Extract original labels
labels_bac <- hc_bac$labels
labels_arc <- hc_arc$labels

# Find highest merge points in each tree
max_bac <- max(hc_bac$height)
max_arc <- max(hc_arc$height)
artificial_height <- max(max_bac, max_arc)

# Create an artificial linkage distance matrix
link_matrix <- matrix(artificial_height, 
                      nrow = length(labels_bac), 
                      ncol = length(labels_arc),
                      dimnames = list(labels_bac, labels_arc))

# Combine distance matrices with artificial links
full_matrix <- rbind(
  cbind(dist_matrix_bac, link_matrix), 
  cbind(t(link_matrix), dist_matrix_arc)
)

# Ensure symmetry
full_matrix[lower.tri(full_matrix)] <- t(full_matrix)[lower.tri(full_matrix)]

# Convert to dist object
dist_obj_combined <- as.dist(full_matrix)

# Perform final hierarchical clustering
hc_combined <- hclust(dist_obj_combined, method = "ward.D2")

# Assign correct labels from full matrix rownames
hc_combined$labels <- rownames(full_matrix)

# base R plot
plot(hc_combined, main = "Combined Bacteria & Archaea Dendrogram", cex = 0.6)

# ggplot plot
dendro_data_combined <- dendro_data(as.dendrogram(hc_combined))

# Prepare labels for ggplot
labels_data_combined <- data.frame(
  label = hc_combined$labels[hc_combined$order],
  x = seq_along(hc_combined$labels)
)

# Plot using ggplot
ggplot() +
  geom_segment(data = dendro_data_combined$segments, 
               aes(x = y, y = x, xend = yend, yend = xend)) +
  geom_text(data = labels_data_combined, 
            aes(x = min(dendro_data_combined$segments$y) - 0.2, y = x, label = label), 
            hjust = 1, size = 3) +  
  expand_limits(x = min(dendro_data_combined$segments$y) - 1.5) +  
  theme_minimal() +
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        panel.grid = element_blank()) +
  labs(title = "Combined Bacteria & Archaea Dendrogram")

## LDG similarity dendrogram
# Filter only the rows where stat_type is "mean"
diversity_mean <- diversity[stat_type == "mean"]

# Drop the 'y' and 'stat_type' columns to get only the clade columns
diversity_wide <- diversity_mean[, !c("y", "stat_type")]

# remove everything until the last "_", including the last "_"
colnames(diversity_wide) = sub(".*_(?!.*_)", "", colnames(diversity_wide), perl = TRUE)

# Compute correlation matrix using Pearson
cor_results <- Hmisc::rcorr(as.matrix(diversity_wide), type = "pearson")

# Extract correlation and p-value matrices
correlation_matrix <- cor_results$r
p_value_matrix <- cor_results$P

# Perform hierarchical clustering on correlation distances
distance_matrix <- as.dist(1 - correlation_matrix)
hc <- hclust(distance_matrix, method = "ward.D2")
clades_order <- hc$labels[hc$order]

# Plot the dendrogram
plot(hc)

## Tanglegram
# plot
tanglegram(untangle(as.dendrogram(hc), as.dendrogram(hc_combined)), highlight_distinct_edges = FALSE)
score <- entanglement(hc, hc_combined)
print(paste("Entanglement Score:", score))


# Create the tanglegram with custom labels
png(filename = paste0("../Output/FigS8_tanglegram_LDG_phylogeny_JS/",dir, "/tanglegram.png"), width = 20, height = 10, units = "cm", res = 1000)
tanglegram(untangle(rev(as.dendrogram(hc)), rev(as.dendrogram(hc_combined))),
           highlight_distinct_edges = FALSE,
           common_subtrees_color_lines = FALSE,
           margin_inner = 9,
           main_left = "Similarity in surface ocean annual LDG",   
           main_right = "Phylogenetic distance",
           cex_main = 1,
           cex_main_left = .9,
           cex_main_right = .9,
           edge.lwd = 1,
           main = paste("Entanglement Score:", round(score, 3)))
dev.off()

svg(filename = paste0("../Output/FigS8_tanglegram_LDG_phylogeny_JS/", dir, "/tanglegram.svg"),
    width = 20 / 2.54, height = 10 / 2.54)
tanglegram(
  untangle(rev(as.dendrogram(hc)), rev(as.dendrogram(hc_combined))),
  highlight_distinct_edges = FALSE,
  common_subtrees_color_lines = FALSE,
  margin_inner = 9,
  main_left = "Similarity in surface ocean annual LDG",   
  main_right = "Phylogenetic distance",
  cex_main = 1,
  cex_main_left = .9,
  cex_main_right = .9,
  edge.lwd = 1,
  main = paste("Entanglement Score:", round(score, 3))
)
dev.off()