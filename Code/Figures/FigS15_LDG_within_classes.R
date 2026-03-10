### This script shows the LDGs of the up to five genera within a class.
### Author: Jonas Schiller

library(data.table)
library(tidyr)
library(dplyr)
library(Hmisc)
library(ggplot2)
library(ggdendro)
library(Hmisc)
library(dynamicTreeCut)
library(dendextend)
library(scales)
library(ggrepel)
library(data.table)
library(stringr)
library(ggrepel)
library(grid)
library(RColorBrewer)

# Preparation
### create dirs
# Data
dir.create("../Data/Internal")
dir.create(paste0("../Data/Internal/04.2_genus_LDG"))

# Output
dir.create("./FigS15_LDG_within_classes")
dir.create(paste0("./FigS15_LDG_within_classes"))


### Read data
mean_annual_dt = fread("../Data/External/mean_annual_ldg_genus.csv")
# clades
clades = c(
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

# gtdb
gtdb = fread("../Data/External/mOTUsv4.0.gtdb.taxonomy.80mv.tsv")
# if a cell contains "Unknown", replace entire cell content with an NA
gtdb_na = gtdb[, lapply(.SD, function(x) ifelse(grepl("Unknown", x), NA, x))]

## plot
norm_genus = function(x) {
  x = as.character(x)
  x = str_trim(x)
  x = str_replace(x, "^g__", "")
  x = str_replace(x, "^genus_", "")
  x = str_replace(x, ";.*$", "")
  x
}

# LDG table
setnames(mean_annual_dt, old = c("y", "mean_annual"), new = c("Latitude", "richness"), skip_absent = TRUE)
mean_annual_dt[, genus_key := norm_genus(genus)]

# Taxonomy table
tax_dt = unique(as.data.table(gtdb_na)[, .(domain, phylum, class, order, family, genus)])
tax_dt = tax_dt[!is.na(genus) & genus != ""]
tax_dt[, genus_key := norm_genus(genus)]

# Merge on genus_key (NOT genus)
ldg_tax = merge(
  mean_annual_dt,
  tax_dt,
  by = "genus_key",
  all = FALSE,
  suffixes = c(".ldg", ".tax")
)

message("Rows after merge: ", nrow(ldg_tax))
message("Unique genera (LDG): ", uniqueN(ldg_tax$genus.ldg))
message("Unique classes in merged data: ", uniqueN(ldg_tax$class))

# 4) Parse clades and keep ONLY explicit classes in the provided vector (ignore domain-only entries)
clade_dt = data.table(clade = clades)
clade_dt[, domain_clade := str_match(clade, "domain__([^_]+)")[,2]]
clade_dt[, class_clade  := str_match(clade, "class__([^_]+)$")[,2]]

# keep only those with an explicit class
clade_classes = clade_dt[!is.na(class_clade), .(domain = domain_clade, class = class_clade)]

# enforce facet order = order in your input vector (and de-duplicate while preserving order)
facet_levels = unique(paste0(clade_classes$domain, " | ", clade_classes$class))

# filter strictly to those class/domain pairs
ldg_keep = ldg_tax[clade_classes, on = .(domain, class), nomatch = 0L]

# label facets exactly as desired
ldg_keep[, class_label := paste0(domain, " | ", class)]
ldg_keep[, class_label := factor(class_label, levels = facet_levels)]

# plot extended
# settings
x_min = -80
x_max =  80
x_extend = 60
x_breaks = seq(-80, 80, by = 20)
N = 5

outdir = "./Output/FigS15_LDG_within_classes/"
outfilepng = file.path(outdir, "LDG_genera_within_classes.png")
outfilesvg = file.path(outdir, "LDG_genera_within_classes.svg")

ldg_keep[, Latitude := as.numeric(Latitude)]

# rightmost point per genus per facet (NO x-range restriction)
anchor_dt = ldg_keep[
  !is.na(richness) & !is.na(Latitude),
  .SD[which.max(Latitude)],
  by = .(class_label, genus_key)
]

# top N genera per facet
if (is.finite(N)) {
  top_genera = ldg_keep[
    !is.na(richness),
    .(m = mean(richness, na.rm = TRUE)),
    by = .(class_label, genus_key)
  ][order(-m)][, .SD[1:min(.N, N)], by = class_label]
  
  anchor_dt = anchor_dt[top_genera, on = .(class_label, genus_key), nomatch = 0L]
}

# label info (labels placed in extended region)
label_dt = copy(anchor_dt)
label_dt[, genus_lab := genus_key]
label_dt[, label_x := x_max + x_extend * 0.6]

# keep only genera that will be drawn/labelled (for correct palette sizing)
keep_keys = unique(label_dt$genus_key)
plot_dt = ldg_keep[genus_key %in% keep_keys]

# palette: Spectral without yellow
# Spectral has a yellow-ish center; remove it by filtering out yellow-like hex codes.
spectral_full = brewer.pal(11, "Spectral")
is_yellowish = function(hex) {
  rgb = grDevices::col2rgb(hex)
  r = rgb[1,] / 255; g = rgb[2,] / 255; b = rgb[3,] / 255
  (r > 0.85 & g > 0.80 & b < 0.55) | (r > 0.85 & g > 0.85 & b > 0.55) # catches the typical Spectral yellow
}
spectral_no_yellow = spectral_full[!is_yellowish(spectral_full)]

# map colors to genera (consistent across facets)
need = length(keep_keys)
cols = rep(spectral_no_yellow, length.out = need)
names(cols) = keep_keys

# plot
p = ggplot(plot_dt, aes(x = Latitude, y = richness, group = genus_key, color = genus_key)) +
  geom_line(alpha = 0.6, linewidth = 0.5, na.rm = TRUE) +
  ggrepel::geom_label_repel(
    data = label_dt,
    aes(
      x = label_x,
      y = richness,
      label = genus_lab,
      fill = genus_key
    ),
    color = "black",            
    alpha = 0.55,               
    direction = "y",
    hjust = 0,
    box.padding = 0.35,        
    point.padding = 0.25,
    force = 2,                  
    force_pull = 0.1,
    max.overlaps = Inf,
    segment.color = NA,         
    size = 4,
    label.size = 0.6,
    label.r = unit(0.15, "lines"),
    seed = 1,
    show.legend = FALSE
  ) +
  facet_wrap(~ class_label, scales = "free_y", drop = TRUE, nrow = 2) +
  coord_cartesian(clip = "off") +
  scale_x_continuous(
    limits = c(x_min, x_max + x_extend),
    breaks = x_breaks,
    minor_breaks = NULL
  ) +
  scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) +
  scale_color_manual(values = cols, guide = "none") +
  scale_fill_manual(values = cols, guide = "none") +
  theme_minimal() +
  theme(
    text = element_text(size = 25),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 17),
    axis.text.y = element_text(size = 17),
    strip.text = element_text(size = 16),
    panel.grid = element_blank(),
    axis.ticks = element_line(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2),
    panel.spacing = unit(0.3, "lines"),
    strip.background = element_blank(),
    legend.position = "none",
    plot.margin = margin(5.5, 160, 5.5, 5.5)
  ) +
  labs(x = "Latitude (°)", y = "", title = "")

p

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
ggsave(outfilepng, p, width = 55, height = 20, units = "cm", dpi = 300)
ggsave(outfilesvg, p, width = 55, height = 20, units = "cm", dpi = 300)