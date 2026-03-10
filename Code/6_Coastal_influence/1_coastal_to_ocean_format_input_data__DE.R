## ====================================================================================================
## This R script formats prokaryotic diversity data for coastal influence analysis by processing
## richness indices and merging them with environmental metadata. It creates distance-based groups
## for modeling the effects of coastal proximity on marine prokaryotic diversity patterns.
##
## Author:       Dominic Eriksson
## Date:         26th of February 2026
## Affiliation:  Environmental Physics Group UP, ETH Zürich, Switzerland
## Contact:      deriksson@ethz.ch
##
## Input files:  
##   - Alpha diversity indices CSV file with prokaryotic richness data
##   - Metadata CSV file with sample location and environmental information
##
## Output files: 
##   - Formatted CSV file with combined diversity
##   - Metadata for different distance groups to the coast, structured for CEPHALOPOD modeling.
##
## Strategy:
##   The script processes alpha diversity data by merging it with metadata,
##   then creates distance-based subsets representing different levels of coastal influence.
##   Each distance group filters samples based on minimum distance to coast, enabling
##   analysis of how coastal proximity affects prokaryotic diversity patterns across
##   oceanic gradients.
##
## Required R packages (tested versions):
##   - dplyr      1.1.2
##   - data.table 1.14.8
## ====================================================================================================

# Load libraries
library(dplyr)
library(data.table)

# Set working directories
wd_in <- "Data/Jonas_Richter/Alpha diversity indices/alpha_diversity_prokaryotes_domain_phylum_class_order_family_genus_rarefied5000_V3.csv"


# Load data
df <- data.table::fread(wd_in)

# Subset clade of interest - biosample, strings with prokaryotes and strings with domain - We only do prokaryotes at the moment
df_diversity <- df %>%
    select(biosample, prokaryotes_Richness_5000)

## Load data - NOTE: Adjust metadata version
df_metadata <- data.table::fread("Data/Jonas_Richter/Metadaten/metadata_V23.csv")
# Subset data of interest in metadata
df_metadata <- df_metadata[, c("biosample", "depth", "month", "year", "lon", "lat", "analysis_type", "distance_to_coast")]

## Merge Jonas diversity with metadata
df <- base::merge(df_diversity, df_metadata, by = "biosample", all.x = TRUE)

# Columns for CEPHALOPOD
cols <- c(
    "scientificname",
    "worms_id",
    "decimallatitude",
    "decimallongitude",
    "depth",
    "year",
    "month",
    "measurementvalue",
    "measurementunit",
    "taxonrank"
)

## Surface to MLD --------------------------------------------------------------------------------------------

# Set output directory
wd_out <- "Code/6_Coastal_influence/1_Output/"

# Filter for surface data
unique(df$analysis_type)
df_surface <- df[which(df$analysis_type != "not_included"), ]
df_surface <- df_surface[analysis_type %in% c("MLD")]

# Define all distance groups
distance_groups <- unique(c(seq(5, 100, by = 5), seq(100, 500, by = 50), seq(500, 2500, by = 100)))
groups <- paste0(  "d_", distance_groups)

l_groups <- list()
for(i in groups){

    # Print progress
    message("Processing group: ", i)

    # Subset
    df_sub <- df_surface %>%
    filter(distance_to_coast >= as.numeric(gsub("d_", "", i)))

    # Create dataframe for CEPHALOPOD
    df_final <- data.frame(
        scientificname = paste0("prokaryotes_Richness_5000_d_", i),
        worms_id = paste0("prokaryotes_Richness_5000_d_", i),
        decimallatitude = df_sub$lat,
        decimallongitude = df_sub$lon,
        depth = df_sub$depth,
        year = df_sub$year,
        month = df_sub$month,
        measurementvalue = df_sub$prokaryotes_Richness_5000,
        measurementunit = paste0("prokaryotes_Richness_5000_d_", i),
        taxonrank = paste0("prokaryotes_Richness_5000_d_", i)
    ) 

    # # Save data
    # write.csv(
    #     df_final, 
    #     paste0(wd_out, "prokaryotes_Richness_5000_", i, "_MLD.csv"), 
    #     row.names = FALSE
    # )   

    # Store in list
    l_groups[[i]] <- df_final

}

# Combine all groups into one dataframe
df_all_groups <- do.call(rbind, l_groups)

# Save combined data
write.csv(
    df_all_groups, 
    paste0(wd_out, "prokaryotes_Richness_5000_all_distance_groups_MLD.csv"), 
    row.names = FALSE)

## ====================================================================================================
## END OF SCRIPT
## ====================================================================================================


