### This script calculates richness (after rarefaction) for overall prokaryotes and classes Alphaproteobacteria and Cyanobacteriia from a taxonomic profile that is being summarized on genus level
### Author: Jonas Schiller

# load dependencies
library(parallel)
library(data.table)
library(Matrix)
library(ggplot2)
library(vegan)
library(dplyr)
library(reshape2)

# Data
dir.create("../Data")
dir.create("../Data/Internal")
dir.create("../Data/External") #should already contain mOTUsv4.profiles_subset_to_OMDB.tsv.gz and mOTUsv4.0.gtdb.taxonomy.80mv.tsv and samples_contextual.xlsx
dir.create("../Data/Internal/FigS4_alpha_diversity_clade_level_GENUS_JS")

# Output
dir.create("./Output")
dir.create("./Output/FigS4_alpha_diversity_clade_level_GENUS_JS")

## Read
# mOTU
motu = fread("../Data/External/01_motus_meta/mOTUsv4.profiles_subset_to_OMDB.tsv.gz")
colnames(motu) = make.names(colnames(motu))
colnames(motu) = gsub("_METAG","",colnames(motu))
colnames(motu) = sub("^[^_]+_", "",colnames(motu))
colnames(motu) <- gsub("\\.", "-", colnames(motu))
# GTDB
gtdb = fread("../Data/External/mOTUsv4.0.gtdb.taxonomy.80mv.tsv")
# if a cell contains "Unknown", replace entire cell content with an NA
gtdb_na <- gtdb[, lapply(.SD, function(x) ifelse(grepl("Unknown", x), NA, x))]
# consider only mOTUs that are in the taxonomic profile
index = gtdb_na$motu %in% motu[[1]]
gtdb_na_select <- gtdb_na[index, ]
# meta data
sheet2 = as.data.table(read_excel("../Data/External/samples_contextual.xlsx", sheet = 2))
sheet4 = as.data.table(read_excel("../Data/External/samples_contextual.xlsx", sheet = 4))
meta = merge(sheet2, sheet4, by = "biosample", all = TRUE)

## Generate a genus-level taxonomic profile
# Extract the count columns (everything except V1)
motu_counts <- motu[, -1, with = FALSE]
motu_ids <- motu[[1]]
# Join taxonomy to motu table
motu[, motu_id := V1]
motu_tax <- merge(motu, gtdb_na_select[, .(motu, genus)], by.x = "motu_id", by.y = "motu", all.x = TRUE)
# Identify sample columns (all except taxonomic columns and motu_id)
sample_cols <- setdiff(names(motu_tax), c("motu_id", "V1", "genus"))
# Separate assigned and unassigned
motu_assigned <- motu_tax[!is.na(genus)]
motu_unassigned <- motu_tax[is.na(genus)]
# Sum counts by genus
motu_genus <- motu_assigned[, lapply(.SD, sum), by = genus, .SDcols = sample_cols]
# Sum unassigned counts and add to the 'unassigned' genus
unassigned_counts <- motu_unassigned[, lapply(.SD, sum), .SDcols = sample_cols]
unassigned_counts[, genus := "unassigned"]
# Combine assigned and unassigned genus counts
motu_genus <- rbind(motu_genus, unassigned_counts, fill = TRUE)
# order rows alphabetically or by total abundance
motu_genus <- motu_genus[order(genus)]

# Build proper taxonomic label with rank names
gtdb_na_select[, full_tax_string := paste0(
  "domain__", domain,
  "_phylum__", phylum,
  "_class__", class,
  "_order__", order,
  "_family__", family,
  "_genus__", genus
)]

## Replace genus name with full taxonomy label
# Map unique genus names to full taxonomy string
genus_to_tax <- unique(gtdb_na_select[!is.na(genus), .(genus, full_tax_string)])
# Merge with motu_genus
motu_genus <- merge(
  motu_genus,
  genus_to_tax,
  by = "genus",
  all.x = TRUE,
  sort = FALSE
)
# Replace genus name with full taxonomy label, fallback to original if no match
motu_genus[, genus := ifelse(!is.na(full_tax_string), full_tax_string, genus)]
# Drop helper column
motu_genus[, full_tax_string := NULL]
# Save updated table
fwrite(motu_genus, "../Data/Internal/FigS4_alpha_diversity_clade_level_GENUS_JS/motu_genus.csv")

## compare unassigned fraction between motu and motu_genus
# Define sample columns
sample_cols <- setdiff(names(motu), c("V1", "motu_id"))
# Compute unassigned fraction in original motu table
motu_total <- motu[, lapply(.SD, sum), .SDcols = sample_cols]
motu_unassigned <- motu[V1 == "mOTUv4.0_unassigned", ..sample_cols]
frac_motu <- as.data.table(t(motu_unassigned / motu_total))
setnames(frac_motu, "V1", "fraction")
frac_motu[, sample := rownames(motu_unassigned)]
frac_motu[, source := "motu"]
# Compute unassigned fraction in genus-level table
genus_total <- motu_genus[, lapply(.SD, sum), .SDcols = sample_cols]
genus_unassigned <- motu_genus[genus == "unassigned", ..sample_cols]
frac_genus <- as.data.table(t(genus_unassigned / genus_total))
setnames(frac_genus, "V1", "fraction")
frac_genus[, sample := rownames(genus_unassigned)]
frac_genus[, source := "motu_genus"]
# Combine for plotting
frac_combined <- rbind(frac_motu, frac_genus)
# plot
# Compute medians
medians <- frac_combined[, .(median_frac = round(median(fraction), 3)), by = source]
# Plot with labels
ggplot(frac_combined, aes(x = source, y = fraction)) +
  geom_boxplot() +
  geom_text(data = medians, aes(x = source, y = median_frac + 0.05, label = median_frac), 
            vjust = 0, fontface = "bold") +
  ylab("Fraction of Unassigned Counts") +
  xlab("") +
  ggtitle("Comparison of Unassigned Fractions per Sample") +
  theme_minimal()
ggsave("./Output/FigS4_alpha_diversity_clade_level_GENUS_JS/fraction_unassigned_motu_species_genus.png")

## Calculate alpha diversity of classes Alphaproteobacteria and Cyanobacteriia
# Transpose: samples as rows, taxa as columns
motu_genus_t <- as.data.table(t(motu_genus[, -1]))
colnames(motu_genus_t) <- motu_genus$genus
motu_genus_t[, biosample := colnames(motu_genus)[-1]]
setcolorder(motu_genus_t, "biosample")
# Drop unassigned genus if present
motu_genus_assigned <- motu_genus_t[, !"unassigned", with = FALSE]
# Identify genus columns containing target classes
alpha_cols <- grep("class__Alphaproteobacteria", names(motu_genus_assigned), value = TRUE)
cyano_cols <- grep("class__Cyanobacteriia", names(motu_genus_assigned), value = TRUE)
# Initialize result table with capitalized names
alpha_div <- data.frame(
  biosample = motu_genus_assigned$biosample,
  
  Shannon_5000 = NA_real_,
  Richness_5000 = NA_integer_,
  Chao1_5000 = NA_real_,
  
  Shannon_1000 = NA_real_,
  Richness_1000 = NA_integer_,
  Chao1_1000 = NA_real_,
  
  Shannon_Alphaproteo_5000 = NA_real_,
  Richness_Alphaproteo_5000 = NA_integer_,
  Chao1_Alphaproteo_5000 = NA_real_,
  
  Shannon_Cyano_5000 = NA_real_,
  Richness_Cyano_5000 = NA_integer_,
  Chao1_Cyano_5000 = NA_real_,
  
  Shannon_Alphaproteo_1000 = NA_real_,
  Richness_Alphaproteo_1000 = NA_integer_,
  Chao1_Alphaproteo_1000 = NA_real_,
  
  Shannon_Cyano_1000 = NA_real_,
  Richness_Cyano_1000 = NA_integer_,
  Chao1_Cyano_1000 = NA_real_
)
# Loop over samples
for (i in 1:nrow(motu_genus_assigned)) {
  counts <- as.numeric(motu_genus_assigned[i, -1])
  total_counts <- sum(counts)
  
  if (total_counts >= 1000) {
    rare_1000 <- rrarefy(matrix(counts, nrow = 1), sample = 1000)
    rare_1000_vec <- as.numeric(rare_1000[1, ])
    
    alpha_div$Shannon_1000[i] <- diversity(rare_1000_vec)
    alpha_div$Richness_1000[i] <- sum(rare_1000_vec > 0)
    alpha_div$Chao1_1000[i] <- as.numeric(estimateR(rare_1000_vec)[2])
    
    alpha_1000 <- na.omit(rare_1000_vec[match(alpha_cols, names(motu_genus_assigned)[-1])])
    cyano_1000 <- na.omit(rare_1000_vec[match(cyano_cols, names(motu_genus_assigned)[-1])])
    
    if (sum(alpha_1000) > 0) {
      alpha_div$Shannon_Alphaproteo_1000[i] <- diversity(alpha_1000)
      alpha_div$Richness_Alphaproteo_1000[i] <- sum(alpha_1000 > 0)
      alpha_div$Chao1_Alphaproteo_1000[i] <- as.numeric(estimateR(alpha_1000)[2])
    }
    
    if (sum(cyano_1000) > 0) {
      alpha_div$Shannon_Cyano_1000[i] <- diversity(cyano_1000)
      alpha_div$Richness_Cyano_1000[i] <- sum(cyano_1000 > 0)
      alpha_div$Chao1_Cyano_1000[i] <- as.numeric(estimateR(cyano_1000)[2])
    }
  }
  
  if (total_counts >= 5000) {
    rare_5000 <- rrarefy(matrix(counts, nrow = 1), sample = 5000)
    rare_5000_vec <- as.numeric(rare_5000[1, ])
    
    alpha_div$Shannon_5000[i] <- diversity(rare_5000_vec)
    alpha_div$Richness_5000[i] <- sum(rare_5000_vec > 0)
    alpha_div$Chao1_5000[i] <- as.numeric(estimateR(rare_5000_vec)[2])
    
    alpha_5000 <- na.omit(rare_5000_vec[match(alpha_cols, names(motu_genus_assigned)[-1])])
    cyano_5000 <- na.omit(rare_5000_vec[match(cyano_cols, names(motu_genus_assigned)[-1])])
    
    if (sum(alpha_5000) > 0) {
      alpha_div$Shannon_Alphaproteo_5000[i] <- diversity(alpha_5000)
      alpha_div$Richness_Alphaproteo_5000[i] <- sum(alpha_5000 > 0)
      alpha_div$Chao1_Alphaproteo_5000[i] <- as.numeric(estimateR(alpha_5000)[2])
    }
    
    if (sum(cyano_5000) > 0) {
      alpha_div$Shannon_Cyano_5000[i] <- diversity(cyano_5000)
      alpha_div$Richness_Cyano_5000[i] <- sum(cyano_5000 > 0)
      alpha_div$Chao1_Cyano_5000[i] <- as.numeric(estimateR(cyano_5000)[2])
    }
  }
}
# Save result
fwrite(alpha_div, "../Data/Internal/FigS4_alpha_diversity_clade_level_GENUS_JS/alpha_div_genus.csv")