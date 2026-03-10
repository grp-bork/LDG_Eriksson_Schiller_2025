## ====================================================================================================
## This R script creates chunked data files from combined prokaryotic diversity datasets to enable
## efficient parallel processing in the CEPHALOPOD species distribution modeling pipeline. It splits
## large datasets by taxonomic groups into smaller, manageable files for distributed computing.
##
## Author:       Dominic Eriksson
## Date:         26th of February 2026
## Affiliation:  Environmental Physics Group UP, ETH Zürich, Switzerland
## Contact:      deriksson@ethz.ch
##
## Input files:  
##   - Combined CSV file with all distance groups and diversity data
##
## Output files: 
##   - Multiple chunked CSV files (chunk1.csv, chunk2.csv, etc.)
##   - Each chunk contains data for up to 100 coastal groups for parallel processing
##
## Strategy:
##   The script reads the combined diversity dataset and splits it into smaller chunks based on
##   unique coastal distance identifiers. This chunking approach enables parallel processing of species
##   distribution models across multiple compute cores or nodes, significantly reducing computational
##   time for large-scale biodiversity modeling workflows.
##
## Required R packages (tested versions):
##   - data.table 1.14.8
## ==================================================================================================== 

# Clear workspace
rm(list = ls())

# Load necessary libraries
library(data.table)

# Set directory
wd_in <- "Code/6_Coastal_influence/1_Output/prokaryotes_Richness_5000_all_distance_groups_MLD.csv"
wd_out <- "Code/6_Coastal_influence/1b_Output/"
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
taxonrank_chunks <- split(unique_taxonrank, ceiling(seq_along(unique_taxonrank) / 10))

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
