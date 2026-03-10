### This script shows the median richness across taxonomic ranks, either for epipelagic samples (defined by the mixed layer depth = MLD) or mesopelagic samples (200 - 1000 m depth).
### Author: Jonas Schiller

# load dependencies
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readxl)
library(htmlwidgets)
library(webshot)
library(networkD3)

# Metrics
analysis = "MLD" # "mesopelagic" "MLD"

## create dirs
dir = analysis

# Data
dir.create("../Data")
dir.create("../Data/Internal")
dir.create("../Data/External") #should already contain samples_contextual.xlsx and diversity_metrics.xlsx and mOTUsv4.0.gtdb.taxonomy.80mv.tsv
dir.create("../Data/Internal/FigS2_Sankey")
dir.create(paste0("../Data/Internal/FigS2_Sankey/",dir))

# Output
dir.create("./Output")
dir.create("./Output/FigS2_Sankey")
dir.create(paste0("./Output/FigS2_Sankey/",dir))

## read data
# meta data
sheet2 = as.data.table(read_excel("../Data/External/samples_contextual.xlsx", sheet = 2))
sheet4 = as.data.table(read_excel("../Data/External/samples_contextual.xlsx", sheet = 4))
meta = merge(sheet2, sheet4, by = "biosample", all = TRUE)
# alpha diversity
file = "../Data/External/species_diversity_taxonomic_rank.xlsx"
sheets = excel_sheets(file)
alpha_diversity = NULL
for (s in sheets) {
  dt = as.data.table(read_excel(file, sheet = s))
  if (is.null(alpha_diversity)) {
    alpha_diversity = dt
  } else {
    # keep only columns not already present
    new_cols = setdiff(names(dt), names(alpha_diversity))
    if (length(new_cols) > 0) {
      alpha_diversity = merge(
        alpha_diversity,
        dt[, c("biosample", new_cols), with = FALSE],
        by = "biosample",
        all = TRUE
      )
    }
  }
}
# mOTUs - GTDB link
gtdb = fread("../Data/External/mOTUsv4.0.gtdb.taxonomy.80mv.tsv") # download @ https://omdb.microbiomics.io/repository/ocean/motus-cols
gtdb_na <- gtdb[, lapply(.SD, function(x) ifelse(grepl("Unknown", x), NA, x))]


# filter for Richness
richness_cols = c("biosample",colnames(alpha_diversity)[grep("Richness", colnames(alpha_diversity))])
alpha_diversity = alpha_diversity %>% dplyr::select(all_of(richness_cols))

# join contextual and alpha diversity data
meta_alpha_div = left_join(meta,
                           alpha_diversity)

# filter for water layer (MLD or mesopelagic)
meta_alpha_div = meta_alpha_div %>%
  filter(analysis_type == analysis) %>%
  filter(rarefaction_included == "rarefaction_5000 & rarefaction_1000")

## Average diversity per Clade
## summary statistics
diversity_stats = data.frame(
  median_alpha_diversity = t(meta_alpha_div[,
                                                  lapply(.SD, median, na.rm = TRUE),
                                                  .SDcols = patterns("Richness")]),
  
  mean_alpha_diversity = t(meta_alpha_div[,
                                                lapply(.SD, mean, na.rm = TRUE),
                                                .SDcols = patterns("Richness")]),
  
  sd_alpha_diversity = t(meta_alpha_div[,
                                              lapply(.SD, sd, na.rm = TRUE),
                                              .SDcols = patterns("Richness")]),
  
  non_na_count = t(meta_alpha_div[,
                                        lapply(.SD, function(x) sum(!is.na(x))),
                                        .SDcols = patterns("Richness")]),
  
  non_zero_count = t(meta_alpha_div[,
                                          lapply(.SD, function(x) sum(x != 0, na.rm = TRUE)),
                                          .SDcols = patterns("Richness")])
)

# format output
diversity_stats$taxon = rownames(diversity_stats)
colnames(diversity_stats) = c(
  "median_alpha_diversity",
  "mean_alpha_diversity",
  "sd_alpha_diversity",
  "non_na_count",
  "non_zero_count",
  "taxon"
)
rownames(diversity_stats) = NULL


#### Proportion of mean richness with GTDB ref
prokaryotes_row = "prokaryotes_Richness_5000"
  
  # Taxonomic level
  diversity_stats$taxonomic_level = dplyr::case_when(
    grepl("genus_",  diversity_stats$taxon) ~ "Genus",
    grepl("family_", diversity_stats$taxon) ~ "Family",
    grepl("order_",  diversity_stats$taxon) ~ "Order",
    grepl("class_",  diversity_stats$taxon) ~ "Class",
    grepl("phylum_", diversity_stats$taxon) ~ "Phylum",
    grepl("^domain_", diversity_stats$taxon) ~ "Domain",
    grepl("^prokaryotes_", diversity_stats$taxon) ~ "Prokaryotes",
    TRUE ~ "Unknown"
  )
  
  # prokaryotes richness
  prokaryotes_value = diversity_stats %>%
    filter(taxon == prokaryotes_row) %>%
    pull(mean_alpha_diversity)
  
  # proportion unknown
  proportion_unknown = diversity_stats %>%
    group_by(taxonomic_level) %>%
    summarise(
      sum_mean_richness = sum(mean_alpha_diversity, na.rm = TRUE),
      .groups = "drop"
    )
  proportion_unknown = proportion_unknown %>%
    mutate(
      proportion_unknown = sum_mean_richness /
        sum_mean_richness[taxonomic_level == "Prokaryotes"]
    )
  
  plot_data = proportion_unknown %>%
    filter(taxonomic_level != "Prokaryotes") %>%
    mutate(
      taxonomic_level = factor(
        taxonomic_level,
        levels = c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
      ),
      prokaryotic_richness_coverage = proportion_unknown
    )
  
  # Plot
  x_left = 0.77 - 0.002
  ggplot(plot_data,
         aes(x = taxonomic_level, y = prokaryotic_richness_coverage)) +
    geom_col(fill = "white", color = "black", linewidth = 0.4, alpha = .5) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
    coord_flip(ylim = c(0.8, 1.2), clip = "off") +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 14),
      panel.grid = element_blank(),
      title = element_text(size = 12)
    ) +
    labs(
      x = NULL,
      y = "",
      title = "Fraction of\nprokaryotic richness\nwith GTDB reference"
    ) +
    geom_text(
      aes(x = taxonomic_level, y = x_left, label = taxonomic_level),
      inherit.aes = FALSE,
      hjust = 0, fontface = "bold", color = "black", size = 5
    ) +
    geom_text(
      aes(label = sprintf("%.3f%%", 100 * prokaryotic_richness_coverage)),
      hjust = -0.05,
      size = 5,
      color = "black"
    )
  ggsave(
    file = paste0("./Output/FigS2_Sankey/", dir, "/fraction_GTDB_ref.png"),
    width = 5, height = 22, units = "cm"
  )
  ggsave(
    file = paste0("./Output/FigS2_Sankey/", dir, "/fraction_GTDB_ref.svg"),
    width = 5, height = 22, units = "cm"
  )


#### Sankey
# keep the original column names around
original_colnames = colnames(meta_alpha_div)

# find the metric columns such as _Richness_5000
metric_cols = grep(paste0("_", "Richness"), original_colnames, value = TRUE)

# rename taxon columns using full GTDB paths where possible
class_cols  = metric_cols[grepl("^class_", metric_cols)]
class_names <- sub("^class_(.*)_Richness.*$", "\\1", class_cols)

class_lut = unique(gtdb_na[, .(class, phylum, domain)])
matched_class = class_lut[match(class_names, class_lut$class)]

valid_class = !is.na(matched_class$class) & !is.na(matched_class$phylum) & !is.na(matched_class$domain)
class_cols = class_cols[valid_class]
class_names = class_names[valid_class]
matched_class = matched_class[valid_class, ]

new_class_colnames = paste0(
  "domain__", matched_class$domain,
  "_phylum__", matched_class$phylum,
  "_class__", class_names,
  "_", "Richness"
)

order_cols  = metric_cols[grepl("^order_", metric_cols)]
order_names <- sub("^order_(.*)_Richness.*$", "\\1", order_cols)

order_lut = unique(gtdb_na[, .(order, class, phylum, domain)])
matched_order = order_lut[match(order_names, order_lut$order)]

valid_order = !is.na(matched_order$order) & !is.na(matched_order$class) &
  !is.na(matched_order$phylum) & !is.na(matched_order$domain)
order_cols = order_cols[valid_order]
order_names = order_names[valid_order]
matched_order = matched_order[valid_order, ]

new_order_colnames = paste0(
  "domain__", matched_order$domain,
  "_phylum__", matched_order$phylum,
  "_class__", matched_order$class,
  "_order__", order_names,
  "_", "Richness"
)

family_cols  = metric_cols[grepl("^family_", metric_cols)]
family_names <- sub("^family_(.*)_Richness.*$", "\\1", family_cols)

family_lut = unique(gtdb_na[, .(family, order, class, phylum, domain)])
matched_family = family_lut[match(family_names, family_lut$family)]

valid_family = !is.na(matched_family$family) & !is.na(matched_family$order) & !is.na(matched_family$class) &
  !is.na(matched_family$phylum) & !is.na(matched_family$domain)
family_cols = family_cols[valid_family]
family_names = family_names[valid_family]
matched_family = matched_family[valid_family, ]

new_family_colnames = paste0(
  "domain__", matched_family$domain,
  "_phylum__", matched_family$phylum,
  "_class__", matched_family$class,
  "_order__", matched_family$order,
  "_family__", family_names,
  "_", "Richness"
)

genus_cols  = metric_cols[grepl("^genus_", metric_cols)]
genus_names <- sub("^genus_(.*)_Richness.*$", "\\1", genus_cols)

genus_lut = unique(gtdb_na[, .(genus, family, order, class, phylum, domain)])
matched_genus = genus_lut[match(genus_names, genus_lut$genus)]

valid_genus = !is.na(matched_genus$genus) & !is.na(matched_genus$family) & !is.na(matched_genus$order) &
  !is.na(matched_genus$class) & !is.na(matched_genus$phylum) & !is.na(matched_genus$domain)
genus_cols = genus_cols[valid_genus]
genus_names = genus_names[valid_genus]
matched_genus = matched_genus[valid_genus, ]

new_genus_colnames = paste0(
  "domain__", matched_genus$domain,
  "_phylum__", matched_genus$phylum,
  "_class__", matched_genus$class,
  "_order__", matched_genus$order,
  "_family__", matched_genus$family,
  "_genus__", genus_names,
  "_", "Richness"
)

phylum_cols  = metric_cols[grepl("^phylum_", metric_cols)]
phylum_names <- sub("^phylum_(.*)_Richness.*$", "\\1", phylum_cols)

phylum_lut = unique(gtdb_na[, .(phylum, domain)])
matched_phylum = phylum_lut[match(phylum_names, phylum_lut$phylum)]

valid_phylum = !is.na(matched_phylum$phylum) & !is.na(matched_phylum$domain)
phylum_cols = phylum_cols[valid_phylum]
phylum_names = phylum_names[valid_phylum]
matched_phylum = matched_phylum[valid_phylum, ]

new_phylum_colnames = paste0(
  "domain__", matched_phylum$domain,
  "_phylum__", phylum_names,
  "_", "Richness"
)

domain_cols  = metric_cols[grepl("^domain_", metric_cols)]
domain_names <- sub("^domain_(.*)_Richness.*$", "\\1", domain_cols)

new_domain_colnames = paste0(
  "domain__", domain_names,
  "_", "Richness"
)

# apply the renaming map
new_colnames_map = setNames(
  c(new_class_colnames, new_order_colnames, new_family_colnames, new_genus_colnames, new_phylum_colnames, new_domain_colnames),
  c(class_cols,       order_cols,       family_cols,       genus_cols,       phylum_cols,       domain_cols)
)
setnames(meta_alpha_div, old = names(new_colnames_map), new = new_colnames_map)

  
cladenames_modeled = c()
  
  # compute mean richness per taxon column
  diversity_stats = data.frame(
    mean_alpha_diversity = t(meta_alpha_div[,
                                                     lapply(.SD, mean, na.rm = TRUE),
                                                     .SDcols = patterns("Richness")])
  )
  diversity_stats$taxon = rownames(diversity_stats)
  rownames(diversity_stats) = NULL
  
  diversity_stats$taxonomic_level = dplyr::case_when(
    grepl("_genus__",  diversity_stats$taxon) ~ "Genus",
    grepl("_family__", diversity_stats$taxon) ~ "Family",
    grepl("_order__",  diversity_stats$taxon) ~ "Order",
    grepl("_class__",  diversity_stats$taxon) ~ "Class",
    grepl("_phylum__", diversity_stats$taxon) ~ "Phylum",
    grepl("^domain__", diversity_stats$taxon) ~ "Domain",
    grepl("^prokaryotes_", diversity_stats$taxon) ~ "Prokaryotes",
    TRUE ~ "Unknown"
  )
  
  # strip suffix so parsing works with underscores
  diversity_stats$taxon = gsub(
    pattern = paste0("_", "Richness"),
    replacement = "",
    x = diversity_stats$taxon
  )
  
  mean_tbl = diversity_stats %>%
    dplyr::select(taxon, mean_alpha_diversity, taxonomic_level)
  
  # helper functions to walk the hierarchy safely
  domain_key = function(x) ifelse(grepl("^domain__", x), sub("(_phylum__.*)$", "", x), NA_character_)
  phylum_key = function(x) ifelse(grepl("_phylum__", x), sub("(_class__.*)$", "", x), NA_character_)
  class_key  = function(x) ifelse(grepl("_class__", x),  sub("(_order__.*)$", "", x), NA_character_)
  order_key  = function(x) ifelse(grepl("_order__", x),  sub("(_family__.*)$", "", x), NA_character_)
  family_key = function(x) ifelse(grepl("_family__", x), sub("(_genus__.*)$", "", x), NA_character_)
  
  # prokaryotes row name used as the root node
  prokaryotes_row = paste0("prokaryotes_", "Richness")
  prokaryotes_row = gsub(paste0("_", "Richness"), "", prokaryotes_row)
  
  # pick the main branches to keep
  keep_domains = mean_tbl %>% filter(taxonomic_level == "Domain") %>% pull(taxon)
  keep_phyla   = mean_tbl %>% filter(taxonomic_level == "Phylum", mean_alpha_diversity >= 2) %>% pull(taxon)
  
  classes_tbl = mean_tbl %>%
    filter(taxonomic_level == "Class") %>%
    mutate(pk = phylum_key(taxon))
  
  orders_tbl = mean_tbl %>%
    filter(taxonomic_level == "Order") %>%
    mutate(ck = class_key(taxon))
  
  families_tbl = mean_tbl %>%
    filter(taxonomic_level == "Family") %>%
    mutate(ok = order_key(taxon))
  
  genera_tbl = mean_tbl %>%
    filter(taxonomic_level == "Genus") %>%
    mutate(fk = family_key(taxon))
  
  top_class_per_phylum = classes_tbl %>%
    filter(pk %in% keep_phyla) %>%
    group_by(pk) %>%
    slice_max(order_by = mean_alpha_diversity, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  keep_classes = top_class_per_phylum$taxon
  
  top_order_per_class = orders_tbl %>%
    filter(ck %in% keep_classes) %>%
    group_by(ck) %>%
    slice_max(order_by = mean_alpha_diversity, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  keep_orders = top_order_per_class$taxon
  
  top_family_per_order = families_tbl %>%
    filter(ok %in% keep_orders) %>%
    group_by(ok) %>%
    slice_max(order_by = mean_alpha_diversity, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  keep_families = top_family_per_order$taxon
  
  top_genus_per_family = genera_tbl %>%
    filter(fk %in% keep_families) %>%
    group_by(fk) %>%
    slice_max(order_by = mean_alpha_diversity, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  keep_genera = top_genus_per_family$taxon
  
  # keep all taxa with high mean and keep their strongest downstream chain
  keep_high = mean_tbl %>%
    filter(!is.na(mean_alpha_diversity), mean_alpha_diversity > 10) %>%
    pull(taxon)
  
  keep_high_levels = mean_tbl$taxonomic_level[match(keep_high, mean_tbl$taxon)]
  keep_high_domains = keep_high[keep_high_levels == "Domain"]
  keep_high_phyla   = keep_high[keep_high_levels == "Phylum"]
  keep_high_classes = keep_high[keep_high_levels == "Class"]
  keep_high_orders  = keep_high[keep_high_levels == "Order"]
  keep_high_fams    = keep_high[keep_high_levels == "Family"]
  keep_high_genus   = keep_high[keep_high_levels == "Genus"]
  
  biggest_class_in_phylum = function(ph){
    classes_tbl %>%
      filter(pk == ph) %>%
      slice_max(order_by = mean_alpha_diversity, n = 1, with_ties = FALSE) %>%
      pull(taxon)
  }
  biggest_order_in_class = function(cl){
    orders_tbl %>%
      filter(ck == cl) %>%
      slice_max(order_by = mean_alpha_diversity, n = 1, with_ties = FALSE) %>%
      pull(taxon)
  }
  biggest_family_in_order = function(or){
    families_tbl %>%
      filter(ok == or) %>%
      slice_max(order_by = mean_alpha_diversity, n = 1, with_ties = FALSE) %>%
      pull(taxon)
  }
  biggest_genus_in_family = function(fa){
    genera_tbl %>%
      filter(fk == fa) %>%
      slice_max(order_by = mean_alpha_diversity, n = 1, with_ties = FALSE) %>%
      pull(taxon)
  }
  
  high_chain_from_phyla = unlist(lapply(keep_high_phyla, function(ph){
    cl = biggest_class_in_phylum(ph)
    if(length(cl) == 0) return(character(0))
    or = biggest_order_in_class(cl)
    fa = if(length(or) > 0) biggest_family_in_order(or) else character(0)
    ge = if(length(fa) > 0) biggest_genus_in_family(fa) else character(0)
    unique(c(cl, or, fa, ge))
  }))
  
  high_chain_from_classes = unlist(lapply(keep_high_classes, function(cl){
    or = biggest_order_in_class(cl)
    fa = if(length(or) > 0) biggest_family_in_order(or) else character(0)
    ge = if(length(fa) > 0) biggest_genus_in_family(fa) else character(0)
    unique(c(or, fa, ge))
  }))
  
  high_chain_from_orders = unlist(lapply(keep_high_orders, function(or){
    fa = biggest_family_in_order(or)
    ge = if(length(fa) > 0) biggest_genus_in_family(fa) else character(0)
    unique(c(fa, ge))
  }))
  
  high_chain_from_families = unlist(lapply(keep_high_fams, function(fa){
    ge = biggest_genus_in_family(fa)
    unique(c(ge))
  }))
  
  keep_high_downstream = unique(c(high_chain_from_phyla, high_chain_from_classes, high_chain_from_orders, high_chain_from_families))
  keep_high_all = unique(c(keep_high, keep_high_downstream))
  
  # include ancestors for high taxa so the tree stays connected
  keep_all_anc_domains = unique(domain_key(keep_high_all)); keep_all_anc_domains = keep_all_anc_domains[!is.na(keep_all_anc_domains)]
  keep_all_anc_phyla   = unique(phylum_key(keep_high_all)); keep_all_anc_phyla   = keep_all_anc_phyla[!is.na(keep_all_anc_phyla)]
  keep_all_anc_classes = unique(class_key(keep_high_all));  keep_all_anc_classes = keep_all_anc_classes[!is.na(keep_all_anc_classes)]
  keep_all_anc_orders  = unique(order_key(keep_high_all));  keep_all_anc_orders  = keep_all_anc_orders[!is.na(keep_all_anc_orders)]
  keep_all_anc_fams    = unique(family_key(keep_high_all)); keep_all_anc_fams    = keep_all_anc_fams[!is.na(keep_all_anc_fams)]
  
  # clean up modeled taxa names
  cladenames_modeled_clean = gsub(
    pattern = paste0("_", "Richness"),
    replacement = "",
    x = cladenames_modeled
  )
  
  modeled_domains = unique(domain_key(cladenames_modeled_clean)); modeled_domains = modeled_domains[!is.na(modeled_domains)]
  modeled_phyla   = unique(phylum_key(cladenames_modeled_clean)); modeled_phyla   = modeled_phyla[!is.na(modeled_phyla)]
  modeled_classes = unique(class_key(cladenames_modeled_clean));  modeled_classes = modeled_classes[!is.na(modeled_classes)]
  modeled_orders  = unique(order_key(cladenames_modeled_clean));  modeled_orders  = modeled_orders[!is.na(modeled_orders)]
  modeled_fams    = unique(family_key(cladenames_modeled_clean)); modeled_fams    = modeled_fams[!is.na(modeled_fams)]
  
  # final set of taxa to keep
  keep_taxa = unique(c(
    prokaryotes_row,
    keep_domains, keep_phyla, keep_classes, keep_orders, keep_families, keep_genera,
    keep_high_all, keep_all_anc_domains, keep_all_anc_phyla, keep_all_anc_classes, keep_all_anc_orders, keep_all_anc_fams,
    cladenames_modeled_clean, modeled_domains, modeled_phyla, modeled_classes, modeled_orders, modeled_fams
  ))
  
  taxonomic_data = mean_tbl %>%
    filter(taxon %in% keep_taxa) %>%
    mutate(
      group = case_when(
        grepl("^domain__Archaea", taxon) ~ "Archaea",
        grepl("^domain__Bacteria", taxon) ~ "Bacteria",
        TRUE ~ NA_character_
      )
    )
  
  # build the links for the sankey
  links_prokaryotes = taxonomic_data %>%
    filter(taxonomic_level == "Domain") %>%
    mutate(source = "Prokaryotes") %>%
    dplyr::select(source, taxon, mean_alpha_diversity) %>%
    rename(target = taxon, value = mean_alpha_diversity)
  
  links_domain = taxonomic_data %>%
    filter(taxonomic_level == "Phylum") %>%
    mutate(source = sub("_phylum__.*", "", taxon)) %>%
    dplyr::select(source, taxon, mean_alpha_diversity) %>%
    rename(target = taxon, value = mean_alpha_diversity)
  
  links_phylum = taxonomic_data %>%
    filter(taxonomic_level == "Class") %>%
    mutate(source = sub("_class__.*", "", taxon)) %>%
    dplyr::select(source, taxon, mean_alpha_diversity) %>%
    rename(target = taxon, value = mean_alpha_diversity)
  
  links_class = taxonomic_data %>%
    filter(taxonomic_level == "Order") %>%
    mutate(source = sub("_order__.*", "", taxon)) %>%
    dplyr::select(source, taxon, mean_alpha_diversity) %>%
    rename(target = taxon, value = mean_alpha_diversity)
  
  links_order = taxonomic_data %>%
    filter(taxonomic_level == "Family") %>%
    mutate(source = sub("_family__.*", "", taxon)) %>%
    dplyr::select(source, taxon, mean_alpha_diversity) %>%
    rename(target = taxon, value = mean_alpha_diversity)
  
  links_family = taxonomic_data %>%
    filter(taxonomic_level == "Genus") %>%
    mutate(source = sub("_genus__.*", "", taxon)) %>%
    dplyr::select(source, taxon, mean_alpha_diversity) %>%
    rename(target = taxon, value = mean_alpha_diversity)
  
  all_links = bind_rows(links_prokaryotes, links_domain, links_phylum, links_class, links_order, links_family)
  
  # build rank tables used for consistent node ordering
  get_domain = domain_key
  get_phylum = phylum_key
  get_class  = class_key
  get_order  = order_key
  get_family = family_key
  
  phylum_rank_tbl = taxonomic_data %>%
    filter(taxonomic_level == "Phylum") %>%
    mutate(domain_k = get_domain(taxon)) %>%
    group_by(domain_k) %>%
    arrange(desc(mean_alpha_diversity), taxon) %>%
    mutate(phylum_rank = row_number()) %>%
    ungroup() %>%
    dplyr::select(taxon, phylum_rank)
  
  class_rank_tbl = taxonomic_data %>%
    filter(taxonomic_level == "Class") %>%
    mutate(phylum_k = get_phylum(taxon)) %>%
    group_by(phylum_k) %>%
    arrange(desc(mean_alpha_diversity), taxon) %>%
    mutate(class_rank = row_number()) %>%
    ungroup() %>%
    dplyr::select(taxon, class_rank)
  
  order_rank_tbl = taxonomic_data %>%
    filter(taxonomic_level == "Order") %>%
    mutate(class_k = get_class(taxon)) %>%
    group_by(class_k) %>%
    arrange(desc(mean_alpha_diversity), taxon) %>%
    mutate(order_rank = row_number()) %>%
    ungroup() %>%
    dplyr::select(taxon, order_rank)
  
  family_rank_tbl = taxonomic_data %>%
    filter(taxonomic_level == "Family") %>%
    mutate(order_k = get_order(taxon)) %>%
    group_by(order_k) %>%
    arrange(desc(mean_alpha_diversity), taxon) %>%
    mutate(family_rank = row_number()) %>%
    ungroup() %>%
    dplyr::select(taxon, family_rank)
  
  genus_rank_tbl = taxonomic_data %>%
    filter(taxonomic_level == "Genus") %>%
    mutate(family_k = get_family(taxon)) %>%
    group_by(family_k) %>%
    arrange(desc(mean_alpha_diversity), taxon) %>%
    mutate(genus_rank = row_number()) %>%
    ungroup() %>%
    dplyr::select(taxon, genus_rank)
  
  nodes = data.frame(taxon = unique(c(all_links$source, all_links$target))) %>%
    mutate(
      level = case_when(
        taxon == "Prokaryotes" ~ 1,
        grepl("^domain__", taxon) & !grepl("_phylum__", taxon) ~ 2,
        grepl("_phylum__", taxon) & !grepl("_class__", taxon) ~ 3,
        grepl("_class__", taxon) & !grepl("_order__", taxon) ~ 4,
        grepl("_order__", taxon) & !grepl("_family__", taxon) ~ 5,
        grepl("_family__", taxon) & !grepl("_genus__", taxon) ~ 6,
        grepl("_genus__", taxon) ~ 7,
        TRUE ~ NA_real_
      ),
      domain_k = get_domain(taxon),
      phylum_k = get_phylum(taxon),
      class_k  = get_class(taxon),
      order_k  = get_order(taxon),
      family_k = get_family(taxon)
    ) %>%
    left_join(taxonomic_data %>% dplyr::select(taxon, mean_alpha_diversity, group), by = "taxon") %>%
    left_join(phylum_rank_tbl, by = "taxon") %>%
    left_join(class_rank_tbl,  by = "taxon") %>%
    left_join(order_rank_tbl,  by = "taxon") %>%
    left_join(family_rank_tbl, by = "taxon") %>%
    left_join(genus_rank_tbl,  by = "taxon") %>%
    mutate(
      pr  = ifelse(!is.na(phylum_rank), sprintf("%05d", phylum_rank), "99999"),
      cr  = ifelse(!is.na(class_rank),  sprintf("%05d", class_rank),  "99999"),
      orr = ifelse(!is.na(order_rank),  sprintf("%05d", order_rank),  "99999"),
      fr  = ifelse(!is.na(family_rank), sprintf("%05d", family_rank), "99999"),
      gr  = ifelse(!is.na(genus_rank),  sprintf("%05d", genus_rank),  "99999"),
      
      pr_branch = case_when(
        level == 3 ~ pr,
        level %in% c(4,5,6,7) ~ {
          ph_node = phylum_rank_tbl$phylum_rank[match(phylum_k, phylum_rank_tbl$taxon)]
          ifelse(is.na(ph_node), "99999", sprintf("%05d", ph_node))
        },
        TRUE ~ "99999"
      ),
      
      cr_branch = case_when(
        level == 4 ~ cr,
        level %in% c(5,6,7) ~ {
          cl_node = class_rank_tbl$class_rank[match(class_k, class_rank_tbl$taxon)]
          ifelse(is.na(cl_node), "99999", sprintf("%05d", cl_node))
        },
        TRUE ~ "99999"
      ),
      
      or_branch = case_when(
        level == 5 ~ orr,
        level %in% c(6,7) ~ {
          od_node = order_rank_tbl$order_rank[match(order_k, order_rank_tbl$taxon)]
          ifelse(is.na(od_node), "99999", sprintf("%05d", od_node))
        },
        TRUE ~ "99999"
      ),
      
      fr_branch = case_when(
        level == 6 ~ fr,
        level == 7 ~ {
          fa_node = family_rank_tbl$family_rank[match(family_k, family_rank_tbl$taxon)]
          ifelse(is.na(fa_node), "99999", sprintf("%05d", fa_node))
        },
        TRUE ~ "99999"
      ),
      
      parent_sort = case_when(
        level == 1 ~ "0",
        level == 2 ~ paste0("1|", taxon),
        level == 3 ~ paste0("2|", domain_k, "|", pr),
        level == 4 ~ paste0("3|", domain_k, "|", pr_branch, "|", cr_branch),
        level == 5 ~ paste0("4|", domain_k, "|", pr_branch, "|", cr_branch, "|", orr),
        level == 6 ~ paste0("5|", domain_k, "|", pr_branch, "|", cr_branch, "|", or_branch, "|", fr),
        level == 7 ~ paste0("6|", domain_k, "|", pr_branch, "|", cr_branch, "|", or_branch, "|", fr_branch, "|", gr),
        TRUE ~ paste0("9|", taxon)
      )
    ) %>%
    arrange(level, parent_sort, desc(mean_alpha_diversity), taxon) %>%
    mutate(node_id = row_number() - 1)
  
  links = all_links %>%
    left_join(nodes %>% dplyr::select(taxon, node_id), by = c("source" = "taxon")) %>%
    rename(source_id = node_id) %>%
    left_join(nodes %>% dplyr::select(taxon, node_id), by = c("target" = "taxon")) %>%
    rename(target_id = node_id) %>%
    dplyr::select(source_id, target_id, value)
  
  # build the display labels
  nodes$display_taxon = dplyr::case_when(
    grepl("_genus__",  nodes$taxon) ~ sub(".*_genus__",  "", nodes$taxon),
    grepl("_family__", nodes$taxon) ~ sub(".*_family__", "", nodes$taxon),
    grepl("_order__",  nodes$taxon) ~ sub(".*_order__",  "", nodes$taxon),
    grepl("_class__",  nodes$taxon) ~ sub(".*_class__",  "", nodes$taxon),
    grepl("_phylum__", nodes$taxon) ~ sub(".*_phylum__", "", nodes$taxon),
    grepl("^domain__", nodes$taxon) ~ sub("^domain__", "", nodes$taxon),
    TRUE ~ nodes$taxon
  )
  
  nodes$taxon = ifelse(
    !is.na(nodes$mean_alpha_diversity),
    paste0(nodes$display_taxon, " [", round(nodes$mean_alpha_diversity, 0), "]"),
    nodes$display_taxon
  )
  
  # render the sankey
  colour_scale = 'd3.scaleOrdinal()
    .domain(["Archaea", "Bacteria"])
    .range(["darkgray", "black"])'
  
  sankey = sankeyNetwork(
    Links = links,
    Nodes = nodes,
    Source = "source_id",
    Target = "target_id",
    Value  = "value",
    NodeID = "taxon",
    NodeGroup = "group",
    colourScale = colour_scale,
    units = "mean_alpha_diversity",
    fontSize = 20,
    nodeWidth = 20,
    nodePadding = 18,
    iterations = 0,
    fontFamily = "Sarif",
    width = 1400,
    height = 800
  )
  
  sankey = htmlwidgets::onRender(
    sankey,
    '
    function(el, x) {
      d3.select(el).selectAll(".node text")
        .attr("x", -10)
        .attr("text-anchor", "end");
    }
    '
  )
  
  sankey
  
  htmlwidgets::saveWidget(
    sankey,
    paste0("./Output/FigS2_Sankey/", dir, "/mean_alpha_diversity_richness_clade_of_prokaryotes_sankey.html")
  )
  
  try(webshot::install_phantomjs(), silent = TRUE)
  
  try(
    webshot::webshot(
      paste0("./Output/FigS2_Sankey/", dir, "/mean_alpha_diversity_richness_clade_of_prokaryotes_sankey.html"),
      file = paste0("./Output/FigS2_Sankey/", dir, "/mean_alpha_diversity_richness_clade_of_prokaryotes_sankey.png")
    ),
    silent = TRUE
  )