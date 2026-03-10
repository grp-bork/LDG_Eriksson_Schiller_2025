### This script calculates alpha diversity metrics (richness, Shannon, Chao1) after rarefaction (e.g., to 5,000 mOTU counts) on the species-corresponding mOTU level for across all prokaryotes and on taxonomic ranks (e.g., domain, phylum, class)
### Code can be run in parallel (choose number of cores in line 18)
### Author: Jonas Schiller

# load dependencies
library(parallel)
library(data.table)
library(vegan)
library(dplyr)

# create dir
dir.create("../Data")
dir.create("../Data/Internal")
dir.create("../Data/Internal/01_alpha_diversity_clade_level_JS") # directory for output files of this script
dir.create("../Data/External") # store mOTUsv4.profiles_subset_to_OMDB.tsv.gz and mOTUsv4.0.gtdb.taxonomy.80mv.tsv here

## Choose on how many cores you want to run the code in parallel (default = 1)
n_cores = 1

## Read taxonomic mOTU profile
motu = fread("../Data/External/mOTUsv4.profiles_subset_to_OMDB.tsv.gz")
# columns are sample IDs --> clean them
colnames(motu) = make.names(colnames(motu))
colnames(motu) = gsub("_METAG","",colnames(motu))
colnames(motu) = sub("^[^_]+_", "",colnames(motu))
colnames(motu) = gsub("\\.", "-", colnames(motu))
## rm the motu col (store mOTU IDs as a vector)
motu_ids = motu$V1
motu = motu %>% dplyr::select(-V1)


## Read GTDB taxonomic assignment per mOTU
gtdb = fread("../Data/External/mOTUsv4.0.gtdb.taxonomy.80mv.tsv")
# if a cell contains "Unknown", replace entire cell content with an NA
gtdb_na = gtdb[, lapply(.SD, function(x) ifelse(grepl("Unknown", x), NA, x))]
# consider only mOTUs that are in the taxonomic motu profile
index = gtdb_na$motu %in% motu[[1]]
gtdb_na_select = gtdb_na[index, ]


## Rarefaction and alpha diversity calculation
# seed for reproducibility
set.seed(1234)

# creation of a data frame to store alpha diversity results
richness_index = data.frame(biosample = colnames(motu))

# definition of the taxonomic levels to calculate alpha diversity for
taxonomic_levels = c("prokaryotes", "domain", "phylum", "class", "order", "family", "genus")

# definition of number of rarefaction iterations and rarefaction sample size
num_iterations = 200
sample_size = 5000

# loop through each taxonomic level
for (taxonomic_level in taxonomic_levels) {
  
  richness_index = data.frame(biosample = colnames(motu))
  
  # retrieving unique taxa for current taxonomic level
  if (taxonomic_level == "prokaryotes") {
    unique_values = "prokaryotes"
  } else {
    unique_values = na.omit(unique(gtdb_na_select[[taxonomic_level]]))
  }
  
  # initialization of columns in richness_index
  for (value in unique_values) {
    richness_index[[paste(taxonomic_level, value, "Richness", sample_size, sep = "_")]] = numeric(ncol(motu))
    richness_index[[paste(taxonomic_level, value, "Shannon", sample_size, sep = "_")]] = numeric(ncol(motu))
    richness_index[[paste(taxonomic_level, value, "Chao1", sample_size, sep = "_")]] = numeric(ncol(motu))
  }
  
  # loop through each sample
  for (i in 1:ncol(motu)) {
    
    tryCatch({
      
      # generating a unique set of seeds for this sample
      seeds = sample(1:1e6, num_iterations)
      
      # function to calculate alpha diversity metrics for one iteration
      iteration_function = function(rep, seed) {
        set.seed(seed)
        suppressMessages(tryCatch({
          
          rarefied_row = vegan::rrarefy(motu[, ..i], sample = sample_size)
          
          sapply(unique_values, function(value) {
            if(value == "prokaryotes"){
              clade_counts = rarefied_row[-length(rarefied_row)]
            } else{
              taxa_indices = which(gtdb_na_select[[taxonomic_level]] == value)
              clade_counts = rarefied_row[taxa_indices]
            }
            
            richness = sum(clade_counts > 0, na.rm = TRUE)
            
            shannon = if (sum(clade_counts, na.rm = TRUE) > 0) {
              vegan::diversity(clade_counts, index = "shannon")
            } else { NA }
            
            chao1 = if (sum(clade_counts, na.rm = TRUE) > 0) {
              as.numeric(vegan::estimateR(clade_counts)[2])
            } else { NA }
            
            c(Richness = richness, Shannon = shannon, Chao1 = chao1)
          })
          
        }, error = function(e) {
          setNames(rep(NA, length(unique_values) * 3),
                   c(paste0(unique_values, "_Richness"),
                     paste0(unique_values, "_Shannon"),
                     paste0(unique_values, "_Chao1")))
        }))
      }
      
      # running parallel rarefaction iterations
      diversity_values = mclapply(1:num_iterations, function(rep) {
        iteration_function(rep, seeds[rep])
      }, mc.cores = n_cores)
      
      # aggregation of results
      for (value in unique_values) {
        richness_col_name = paste(taxonomic_level, value, "Richness", sample_size, sep = "_")
        shannon_col_name = paste(taxonomic_level, value, "Shannon", sample_size, sep = "_")
        chao1_col_name = paste(taxonomic_level, value, "Chao1", sample_size, sep = "_")
        
        richness_index[[richness_col_name]][i] = mean(sapply(diversity_values, function(x) x["Richness", value]), na.rm = TRUE)
        richness_index[[shannon_col_name]][i] = mean(sapply(diversity_values, function(x) x["Shannon", value]), na.rm = TRUE)
        richness_index[[chao1_col_name]][i] = mean(sapply(diversity_values, function(x) x["Chao1", value]), na.rm = TRUE)
      }
      
    }, error = function(e) {
      for (value in unique_values) {
        richness_index[[paste(taxonomic_level, value, "Richness", sample_size, sep = "_")]][i] = NA
        richness_index[[paste(taxonomic_level, value, "Shannon", sample_size, sep = "_")]][i] = NA
        richness_index[[paste(taxonomic_level, value, "Chao1", sample_size, sep = "_")]][i] = NA
      }
    })
    
    print(paste0("Completed ", taxonomic_level, " ", round(i / ncol(motu) * 100, 3), "%"))
  }
  
  # avoid to report alpha diversity metrics if the counts are too low for a given sample (NA instead)
  rows_to_na = colSums(motu) < sample_size
  richness_index[rows_to_na, 2:ncol(richness_index)] = NA
  
  # clean: remove if no non-NA alpha diversity metric >0 in any sample for a given alpha diversity metric of a taxonomic level
  richness_index = richness_index[, c(TRUE, colSums(richness_index[, -1], na.rm = TRUE) > 0)]
  
  # write as .csv
  output_file = paste0("../Data/Internal/01_alpha_diversity_clade_level_JS/alpha_diversity_", taxonomic_level, "_rarefied", sample_size, ".csv")
  fwrite(richness_index, file = output_file, row.names = FALSE)
  print(paste0("Saved results to ", output_file))
}