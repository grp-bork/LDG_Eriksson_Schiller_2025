## ====================================================================================================
## This R script removes marginal seas from basin-stratified subsampling uncertainty analysis to
## focus the sampling bias assessment on open ocean regions. It uses IHO (International Hydrographic
## Organization) ocean region classifications to filter out enclosed and semi-enclosed seas that
## may have different sampling characteristics and oceanographic properties than open ocean environments.
##
## Author:       Dominic Eriksson
## Date:         26th of February 2026
## Affiliation:  Environmental Physics Group UP, ETH Zürich, Switzerland
## Contact:      deriksson@ethz.ch
##
## Input files:  
##   - CSV with coefficient of variation metrics from basin-stratified subsampling analysis
##   - IHO ocean regions function for spatial classification
##
## Output files: 
##   - Filtered CSV with marginal seas removed
##   - Dataset focused on open ocean regions for sampling bias uncertainty assessment
##
## Strategy:
##   The script applies IHO ocean region classification to identify and remove marginal seas
##   including the Mediterranean, Red Sea, Baltic Sea, Black Sea, Persian Gulf, and their
##   associated basins and gulfs. This filtering ensures that sampling bias uncertainty
##   analyses focus on truly oceanic environments and avoid potential confounding effects
##   from enclosed sea dynamics that may respond differently to sampling heterogeneity.
##
## Required R packages (tested versions):
##   - raster     3.6.26
##   - sf         1.0.12
##   - data.table 1.14.8
##   - ggplot2    3.4.2
##   - parallel   4.2.2
##   - dplyr      1.1.2
## ====================================================================================================

# Clear workspace
rm(list = ls())

# Libraries
library(raster)     # For raster data manipulation
library(sf)         # For spatial data handling
library(data.table) # For efficient data handling
library(ggplot2)    # For visualization
library(parallel)   # For parallel processing
library(dplyr)      # For data manipulation

# Load functions
source("../../Functions/iho_ocean_regions.R")

# Directories
wd_in <- "Code/7_Ocean_thinning/5_Output/prokaryote_richness_cv.csv"
wd_out <- "Code/7_Ocean_thinning/6_Output/"

# Create output directory if it doesn't exist
if(!dir.exists(wd_out)) {
    dir.create(wd_out, recursive = TRUE)
}

# Load data
df <- fread(wd_in)

# Rename coordinates for function
names(df)[1:2] <- c("decimalLongitude", "decimalLatitude")

# Apply function to get ocean regions
df_with_regions <- iho_ocean_regions(df)

# Define the marginal seas to remove
marginal_seas <- c(
    "Mediterranean Sea - Eastern Basin",
    "Mediterranean Sea - Western Basin", 
    "Red Sea",
    "Baltic Sea",
    "Gulf of Finland",    # Part of Baltic system
    "Gulf of Bothnia",    # Part of Baltic system
    "Gulf of Riga",       # Part of Baltic system
    "Black Sea",
    "Sea of Azov",        # Connected to Black Sea
    "Persian Gulf",
    "Gulf of Aden",       # Extension of Red Sea
    "Gulf of Aqaba",      # Part of Red Sea
    "Caspian Sea",        # Note: might not be in ocean shapefile as it's landlocked
    "Adriatic Sea",       # Part of Mediterranean
    "Aegean Sea",         # Part of Mediterranean
    "Tyrrhenian Sea",     # Part of Mediterranean
    "Ionian Sea",         # Part of Mediterranean
    "Sea of Marmara",     # Part of Mediterranean system
    "Gulf of Suez",       # Part of Red Sea system
    "Ligurian Sea",       # Part of Mediterranean
    "Balearic (Iberian Sea)",  # Part of Mediterranean
    "Alboran Sea"         # Part of Mediterranean
)

# Filter out marginal seas
df_filtered <- df_with_regions %>%
    filter(!marine_region %in% marginal_seas)

# Check how many rows were removed
cat("Original rows:", nrow(df_with_regions), "\n")
cat("Filtered rows:", nrow(df_filtered), "\n") 
cat("Removed rows:", nrow(df_with_regions) - nrow(df_filtered), "\n")

# Save the filtered dataset
fwrite(df_filtered, file.path(wd_out, "prokaryote_richness_cv_no_marginal_seas.csv"))

## ====================================================================================================
## END OF SCRIPT
## ====================================================================================================