## ====================================================================================================
## This R script creates chunked data files from basin-stratified subsampled prokaryotic diversity
## datasets to enable efficient parallel processing in the CEPHALOPOD species distribution modeling
## pipeline. It splits the large subsampled datasets by taxonomic groups into smaller, manageable
## files optimized for distributed computing across multiple ocean basin thinning scenarios.
##
## Author:       Dominic Eriksson
## Date:         26th of February 2026
## Affiliation:  Environmental Physics Group UP, ETH Zürich, Switzerland
## Contact:      deriksson@ethz.ch
##
## Input files:  
##   - Combined CSV file with all basin-stratified subsampled diversity data
##   - Contains multiple cap values and iteration combinations for sensitivity analysis
##
## Output files: 
##   - Multiple chunked CSV files (chunk1.csv, chunk2.csv, etc.)
##   - Each chunk contains data for up to 100 taxonomic groups for parallel processing
##   - Optimized for evaluating sampling bias effects across ocean basins
##
## Strategy:
##   The script processes the comprehensive subsampled dataset containing 1,100 different
##   sampling scenarios (11 cap values × 100 iterations) and splits them into parallel-
##   processable chunks. This chunking approach enables efficient distributed modeling of
##   sampling bias sensitivity across multiple compute cores or nodes, facilitating robust
##   assessment of how basin-specific sampling intensity affects biodiversity predictions.
##
## Required R packages (tested versions):
##   - data.table 1.14.8
## ==================================================================================================== 

# Clear workspace
rm(list = ls())

# Load libraries
library(data.table) # For efficient data handling

# Set directory
wd_in <- "Code/7_Ocean_thinning/1_Output/prokaryotes_Richness_5000_all_basin_stratified_subsamples_MLD.csv"
wd_out <- "Code/7_Ocean_thinning/1b_Output/"
#
# Check if folder exists, if not create it (including all parent directories)
if (!dir.exists(wd_out)) {
  dir.create(wd_out, recursive = TRUE, showWarnings = FALSE)
  message("Created directory: ", wd_out)
} else {
  message("Directory already exists: ", wd_out)
}

# Load data
df_final <- data.table::fread(wd_in)

## Split the dataframe into chunks

# Get unique taxonrank values
unique_taxonrank <- unique(df_final$taxonrank)

# Split the unique taxonrank values into chunks of 100
taxonrank_chunks <- split(unique_taxonrank, ceiling(seq_along(unique_taxonrank) / 100))

# Loop through each chunk and save the corresponding data
for (i in seq_along(taxonrank_chunks)) {
  
  # Subset the data for the current chunk of taxonrank values
  df_chunk <- df_final[df_final$taxonrank %in% taxonrank_chunks[[i]], ]
  
  # Save the chunk to a CSV file
  data.table::fwrite(
    df_chunk,
    paste0(wd_out, "chunk", i, ".csv"), # Adjust filename if necessary to not overwrite existing ones
    row.names = FALSE
  )
}

## ====================================================================================================
## END OF SCRIPT
## ====================================================================================================