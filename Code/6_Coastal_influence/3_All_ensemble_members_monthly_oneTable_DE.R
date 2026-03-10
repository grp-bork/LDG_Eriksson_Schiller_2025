## ====================================================================================================
## This R script compiles individual ensemble member files from CEPHALOPOD model outputs into a
## long-format table for efficient downstream analysis. It processes raster-based
## predictions across bootstrap replicates and temporal dimensions, creating a unified dataset
## that facilitates ensemble statistics calculation and comparative analyses.
##
## Author:       Dominic Eriksson
## Date:         26th of February 2026
## Affiliation:  Environmental Physics Group UP, ETH Zürich, Switzerland
## Contact:      deriksson@ethz.ch
##
## Input files:  
##   - Individual ensemble member RData files from successful model runs
##   - Raster stacks containing bootstrap replicates across monthly predictions
##
## Output files: 
##   - Comprehensive CSV table with all ensemble members in long format
##   - Columns include spatial coordinates, algorithm, clade, bootstrap, month, and richness values
##   - Optimized structure for calculating ensemble means, standard deviations, and uncertainty metrics
##
## Strategy:
##   The script iterates through all ensemble member files, extracts raster data into data frames,
##   and combines them into a single long-format table. This approach enables efficient computation
##   of ensemble statistics across diverse taxonomic indices while preserving the full uncertainty
##   structure from bootstrap replicates and temporal variability. The resulting table serves as
##   the foundation for robust biodiversity pattern analysis and model validation.
##
## Required R packages (tested versions):
##   - raster     3.6.26
##   - tidyr      1.3.0
##   - dplyr      1.1.2
##   - data.table 1.14.8
## ====================================================================================================

# Clear workspace
rm(list = ls())

# Libraries
library(raster)     # For raster data manipulation
library(tidyr)      # For data tidying
library(dplyr)      # For data manipulation
library(data.table) # For efficient data handling

# Directories
wd_in <- "Code/6_Coastal_influence/2_Output" # for loading ensemble member files from successfully modeled taxa
wd_out <- "Code/6_Coastal_influence/3_Output/"
# Create output directory if it does not exist
if (!dir.exists(wd_out)) {
    dir.create(wd_out, recursive = TRUE)
}


# Get filenames
fnames <- list.files(wd_in, full.names = TRUE)

# We open each ensemble member and save in one dataframe, which is helpful efficient calculation on ensemble means and sd's across diverse indices.
l <- list()
for(f in seq_along(fnames)){

    # Print progress
    print( paste0("Loading file ", f, ", out of ", length(fnames), ".") )
    print( paste0("File name: ", fnames[f]) )

    # Subset
    data <- get(load(fnames[f]))

    # Get index
    idx <- length(data)

  # Save in dataframe
  l2 <- list()
    for(i in 1:idx ){

        df <- raster::as.data.frame(data[[i]], xy = TRUE)
        df$alg <- names(data)[i]
        
        # Add clade indicator
        df$clade <- gsub( "_Richness.*", "", basename(fnames[f]) )
        df$clade <- gsub("ensembleMembers_Log_v2.RData","", basename(fnames[f]) )

        # Store in list object
        l2[[i]] <- df

    }

    # Merge into one dataframe
    df <- do.call("rbind", l2)
    l[[f]] <- df

}

# Merge into one dataframe
df_final <- do.call("rbind", l)

# Check names
unique(df_final$clade)

# Reformat table
df_long <- df_final %>%
  tidyr::pivot_longer(
    cols = starts_with("Bootstrap_"),
    names_to = "bootstrap_month",
    values_to = "richness"
  ) %>%
  tidyr::extract(
    col = bootstrap_month,
    into = c("bootstrap", "month"),
    regex = "Bootstrap_(\\d+)_Month_(\\d+)"
  ) %>%
  mutate(
    bootstrap = as.integer(bootstrap),
    month = as.integer(month)
  )

# Check clades
unique(df_long$clade)
head(df_long)

# Save df_long - Use fwrite for faster writing
data.table::fwrite(
    df_long, 
    file = paste0(wd_out, "df_long_monthly_ensembleMembers_richness_5000_v2.csv")
)

## ====================================================================================================
## END OF SCRIPT
## ====================================================================================================