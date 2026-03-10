### This script compares the similarity in class-level LDGs (pair-wise Pearson correlation of modeled 1° latitude binned richness) to the phylogenetic distance of the classes
### Author: Jonas Schiller, Anna Mankowski

# load dependencies
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggdendro)
library(dendextend)
library(svglite)
library(ape)
library(phytools)

### create dirs
# Data
dir.create("../Data/Internal/FigS14_tanglegram_LDG_phylogeny_JS")
dir.create("../Data/Internal/FigS14_tanglegram_LDG_phylogeny_JS/")

# Output
dir.create("./FigS14_tanglegram_LDG_phylogeny_JS")
dir.create("./FigS14_tanglegram_LDG_phylogeny_JS/")

# Clades
clades_order = c(
  "domain__Bacteria_phylum__Fibrobacterota_class__Fibrobacteria",
  "domain__Archaea_phylum__Thermoproteota_class__Nitrososphaeria",
  "domain__Bacteria_phylum__Pseudomonadota_class__Gammaproteobacteria",
  "domain__Bacteria_phylum__SAR324_class__SAR324",
  "domain__Bacteria_phylum__Bacteroidota_class__Bacteroidia",
  "domain__Archaea_phylum__Thermoplasmatota_class__Poseidoniia",
  "domain__Bacteria_phylum__Marinisomatota_class__Marinisomatia",
  "domain__Bacteria_phylum__Desulfobacterota_D_class__UBA1144",
  "domain__Bacteria_phylum__Actinomycetota_class__Acidimicrobiia",
  "domain__Bacteria_phylum__Cyanobacteriota_class__Cyanobacteriia",
  "domain__Bacteria_phylum__Pseudomonadota_class__Alphaproteobacteria",
  "domain__Bacteria_phylum__Chlamydiota_class__Chlamydiia",
  "domain__Bacteria_phylum__Myxococcota_class__UBA796",
  "domain__Bacteria_phylum__Myxococcota_class__UBA4151",
  "domain__Bacteria_phylum__Myxococcota_class__XYA12-FULL-58-9"
)
# extract only the class__x part
clades_order = paste0("class_", stringr::str_extract(clades_order, "(?<=class__)\\S+"))


## LDG similarity dendrogram
hc = readRDS("../Data/External/hc.rds")
# Convert to dendrogram
dend = as.dendrogram(hc)
# Labels to remove
drop_labels = c("domain__Bacteria", "domain__Archaea")

# Prune dendrogram
dend_pruned = prune(dend, drop_labels)
# Convert back to hclust
hc_pruned = as.hclust(dend_pruned)
hc_pruned$labels = sub(".*_", "", hc_pruned$labels)
plot(hc_pruned)


## Phylogenetic distances
# pruned bacterial backbone
bacttree = read.tree("../Data/External/bac120_r220_pruned.newick")
plot(bacttree)
# pruned archaeal backbone
archtree = read.tree("../Data/External/arc53_r220_pruned.newick")
plot(archtree)

# concatenate for tanglegram
super_tree = bind.tree(
  bacttree,
  archtree,
  where = "root",
  position = 0
)
plot(super_tree)

# label format
super_tree$tip.label = sub(".*_", "", super_tree$tip.label)
# Resolving polytomies deterministically
super_tree_bin = multi2di(super_tree, random = FALSE)
# Coercing to ultrametric
super_tree_u = force.ultrametric(super_tree_bin, method = "extend")

# preserving phylogenetic distance
hc_phylo_exact = as.hclust.phylo(super_tree_u)
dend_phylo = as.dendrogram(hc_phylo_exact)


### Tangelgram
# Correlation-based dendrogram from your LDG similarity clustering
dend_corr = as.dendrogram(hc_pruned)

# Make sure labels are in the same format on both sides
labels(dend_phylo) = sub(".*_", "", labels(dend_phylo))
labels(dend_corr)  = sub(".*_", "", labels(dend_corr))

# If one side has extra/missing tips, restrict both to the shared set
common_labs = intersect(labels(dend_corr), labels(dend_phylo))
dend_corr  = prune(dend_corr,  setdiff(labels(dend_corr),  common_labs))
dend_phylo = prune(dend_phylo, setdiff(labels(dend_phylo), common_labs))

stopifnot(setequal(labels(dend_corr), labels(dend_phylo)))

# Untangle using dendextend (not ape!) to reduce crossings by rotating subtrees
dends_unt = dendextend::untangle(
  dendextend::dendlist(rev(dend_corr), rev(dend_phylo)),
  method = "step2side"
)
# Entanglement is a “how tangled are the connecting lines?” score (0 = perfectly aligned order, 1 = very tangled)
score = dendextend::entanglement(dends_unt)
print(paste("Entanglement Score:", score))
# On-screen tanglegram
dendextend::tanglegram(
  dends_unt,
  highlight_distinct_edges = FALSE,
  common_subtrees_color_lines = FALSE,
  margin_inner = 9,
  main_left = "Similarity in surface mixed layer depth annual LDG",
  main_right = "Phylogenetic distance",
  cex_main = 1,
  cex_main_left = .9,
  cex_main_right = .9,
  edge.lwd = 1,
  main = paste("Entanglement Score:", round(score, 3))
)

# png
png(
  filename = paste0("./FigS14_tanglegram_LDG_phylogeny_JS/tanglegram.png"),
  width = 20, height = 10, units = "cm", res = 1000
)
dendextend::tanglegram(
  dends_unt,
  highlight_distinct_edges = FALSE,
  common_subtrees_color_lines = FALSE,
  margin_inner = 9,
  main_left = "Similarity in surface mixed layer depth annual LDG",
  main_right = "Phylogenetic distance",
  cex_main = 1,
  cex_main_left = .9,
  cex_main_right = .9,
  edge.lwd = 1,
  main = paste("Entanglement Score:", round(score, 3))
)
dev.off()

# svg
svg(
  filename = paste0("./FigS14_tanglegram_LDG_phylogeny_JS/tanglegram.svg"),
  width = 20 / 2.54, height = 10 / 2.54
)
dendextend::tanglegram(
  dends_unt,
  highlight_distinct_edges = FALSE,
  common_subtrees_color_lines = FALSE,
  margin_inner = 9,
  main_left = "Similarity in surface mixed layer depth annual LDG",
  main_right = "Phylogenetic distance",
  cex_main = 1,
  cex_main_left = .9,
  cex_main_right = .9,
  edge.lwd = 1,
  main = paste("Entanglement Score:", round(score, 3))
)
dev.off()