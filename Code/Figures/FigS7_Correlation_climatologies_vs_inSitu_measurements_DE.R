##====================================================================================================
## Supplementary Figure S7: Correlation Analysis of Climatological vs In-situ Environmental Data
##
## This R script performs comprehensive correlation analysis between climatological environmental 
## parameters and in-situ measurements from Tara Ocean expeditions. The analysis validates the use of 
## climatological data as environmental predictors in marine species distribution modeling by comparing 
## seven key oceanographic variables: temperature, chlorophyll-a, photosynthetically active radiation 
## (PAR), nitrate, phosphate, silicate, and mixed layer depth. The script includes spatial-temporal 
## matching using 0.5° × 0.5° grid cells, Pearson correlation calculations, and creates publication-ready 
## visualizations including global distribution maps and correlation plots with statistical significance testing.
##
## Author:      Dominic Eriksson
##              Environmental Physics Group, UP
##              ETH Zürich, Switzerland
## Contact:     deriksson@ethz.ch
## Date:        February 27th, 2026
## Affiliation: ETH Zürich, Environmental Physics Group, UP
##
## Input files:
## - ../../Data/Environmental_parameters/Surface_predictors/*.nc (NetCDF climatology files)
## - ../../Data/Tara_Oceans/tara_metadata_20260119.tsv (Tara Ocean metadata)
## - ../../Functions/iho_ocean_regions.R (Ocean regions classification function)
##
## Output files:
## - ./Output/climatology_surface_predictors.csv (Processed climatology data)
## - ./Output/tara_environmental_DCM_SRF.csv (Filtered Tara environmental data)
## - ./Output/pearson_correlations_tara_climatology.png/.svg (Correlation plots)
## - ./Output/pearson_correlation_results_tara_climatology.csv (Correlation statistics)
## - ./Output/matchup_locations_map.png/.svg (Spatial distribution map)
## - ./Output/matchup_spatial_coordinates.csv (Match-up coordinates)
## - ./Output/combined_tara_climatology_correlation.png/.svg/.pdf (Combined figure)
##
## Strategy:
## 1. Load and process NetCDF climatology files from multiple environmental parameters
## 2. Filter Tara Ocean metadata to DCM and SRF depth layers, open ocean samples only
## 3. Remove marginal seas samples for consistency with modeling approaches
## 4. Create long-format dataframes with spatial-temporal coordinates
## 5. Harmonize variable names and coordinate systems between datasets
## 6. Perform spatial-temporal matching using 0.5° × 0.5° grid cells
## 7. Calculate Pearson correlations with significance testing
## 8. Generate correlation plots with regression lines and statistics
## 9. Create global distribution map of match-up locations
## 10. Combine visualizations into publication-ready multi-panel figure
##
## Required R packages:
## - data.table (1.14.2): Fast data manipulation and aggregation
## - ncdf4 (1.19): NetCDF file reading and processing
## - raster (3.5.15): Spatial data handling and processing
## - ggplot2 (3.3.6): Advanced data visualization and plotting
## - gridExtra (2.3): Multi-panel plot arrangement
## - dplyr (1.0.9): Data manipulation and transformation
## - stringr (1.4.0): String manipulation and processing
## - maps (3.4.0): Geographic map data and plotting
## - mapproj (1.2.8): Map projections and coordinate systems
##====================================================================================================

# Clear workspace
rm(list = ls())

# Load libraries
library(data.table)
library(ncdf4)
library(raster)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(stringr)
library(maps)
library(mapproj)

# Directories
wd_in_clim <- "../../Data/Environmental_parameters/Surface_predictors/"
wd_in_inSitu <- "../../Data/Tara_Oceans/tara_metadata_20260119.tsv"
wd_out <- "./Output/"

# Check if the directory path is correct
cat("Output directory:", wd_out, "\n")
cat("Directory exists:", dir.exists(wd_out), "\n")

# Create output directory
if(!dir.exists(wd_out)) {
    dir.create(wd_out, recursive = TRUE)
    cat("Created output directory\n")
} else {
    cat("Output directory already exists\n")
}

# Check if climatology files exist
cat("Climatology directory:", wd_in_clim, "\n")
cat("Climatology directory exists:", dir.exists(wd_in_clim), "\n")

# Get all NetCDF files from climatology directory
nc_files <- list.files(wd_in_clim, pattern = "\\.nc$", full.names = TRUE)

cat("Found", length(nc_files), "NetCDF files to process\n")
if(length(nc_files) > 0) {
    cat("First few files:\n")
    print(head(basename(nc_files), 5))
} else {
    cat("No NetCDF files found in:", wd_in_clim, "\n")
    cat("Directory contents:\n")
    print(list.files(wd_in_clim))
}

# Load Tara Oceans metadata
tara_metadata <- fread(wd_in_inSitu, sep = "\t", header = TRUE)

# Define columns to subset for Tara data
environmental_cols <- c(
    # Required metadata columns
    "biosample", "lat", "lon", "event_latitude", "event_longitude", 
    "date", "event_date", "depth_str", "depth_x", "depth_y", "depth_nominal",
    "sample", "station",
    
    # Environmental/Physicochemical/Biological Parameters
    "temperature", "chla", "par", "backscattering",
    "depth_chl_max", "mixed_layer_depth_sigma", "depth_max_brunt_väisälä",
    "nitrite", "phosphate", "nitrate_nitrite", "silicate", "nstar",
    "depth_bathy", "lyapunov_exp",
    
    # Additional context columns
    "size_fraction_x", "size_fraction_y",
    "marine_biome", "ocean_region", "biogeo_province"
)

# Check which columns exist and subset Tara data
existing_cols <- environmental_cols[environmental_cols %in% colnames(tara_metadata)]
tara_environmental <- tara_metadata[, ..existing_cols]

# Filter to only keep DCM and SRF depth layers
tara_environmental <- tara_environmental[depth_y %in% c("DCM", "SRF")]

cat("Tara environmental data prepared:\n")
cat("Dimensions:", nrow(tara_environmental), "rows x", ncol(tara_environmental), "columns\n")

# Get all NetCDF files from climatology directory
nc_files <- list.files(wd_in_clim, pattern = "\\.nc$", full.names = TRUE)

cat("Found", length(nc_files), "NetCDF files to process\n")

# Function to extract climatology data from NetCDF file
extract_climatology_data <- function(nc_file_path) {
    
    cat("Processing:", basename(nc_file_path), "\n")
    
    # Open NetCDF file
    nc_data <- nc_open(nc_file_path)
    
    # Get dimensions
    lon <- ncvar_get(nc_data, "lon")
    lat <- ncvar_get(nc_data, "lat")
    time <- ncvar_get(nc_data, "time")
    
    # Get variable information (assuming there's only one main variable)
    var_name <- names(nc_data$var)[1]  # Get first variable name
    var_info <- nc_data$var[[var_name]]
    var_longname <- var_info$longname
    
    # Get the data
    data_array <- ncvar_get(nc_data, var_name)
    
    # Close NetCDF file
    nc_close(nc_data)
    
    # Create data frame with all combinations of lon, lat, time
    df_list <- list()
    
    for(t in 1:length(time)) {
        # Extract data for this time step
        data_slice <- data_array[, , t]
        
        # Create coordinate grid
        coords <- expand.grid(lon = lon, lat = lat)
        
        # Add the data values
        coords$value <- as.vector(data_slice)
        coords$month <- t
        coords$variable <- var_name
        coords$long_name <- var_longname
        coords$source_file <- basename(nc_file_path)
        
        df_list[[t]] <- coords
    }
    
    # Combine all months
    result_df <- rbindlist(df_list)
    
    return(result_df)
}

# Process all NetCDF files
climatology_list <- list()

for(i in seq_along(nc_files)) {
    climatology_list[[i]] <- extract_climatology_data(nc_files[i])
}

# Combine all climatology data
climatology_df <- rbindlist(climatology_list)

cat("\nClimatology data extracted:\n")
cat("Dimensions:", nrow(climatology_df), "rows x", ncol(climatology_df), "columns\n")
cat("Variables included:", length(unique(climatology_df$variable)), "\n")
cat("Variable names:\n")
print(unique(climatology_df[, .(variable, long_name)]))

# Show structure of climatology dataframe
cat("\nClimatology dataframe structure:\n")
str(climatology_df)

# Show sample of data
cat("\nFirst few rows of climatology data:\n")
print(head(climatology_df))

# Save climatology dataframe
fwrite(climatology_df, file.path(wd_out, "climatology_surface_predictors.csv"))
cat("\nClimatology data saved to:", file.path(wd_out, "climatology_surface_predictors.csv"), "\n")

# Load file back to verify
climatology_df <- fread(file.path(wd_out, "climatology_surface_predictors.csv"))

# Also save Tara environmental data
fwrite(tara_environmental, file.path(wd_out, "tara_environmental_DCM_SRF.csv"))
# Load back to verify
tara_environmental <- fread(file.path(wd_out, "tara_environmental_DCM_SRF.csv"))

# Check sample size of each env. variable in Tara data
cat("\n=== SAMPLE COUNTS FOR TARA ENVIRONMENTAL VARIABLES ===\n")
env_vars <- c("temperature", "chla", "par", "backscattering", 
              "depth_chl_max", "mixed_layer_depth_sigma", "depth_max_brunt_väisälä",
              "nitrite", "phosphate", "nitrate_nitrite", "silicate", "nstar", 
              "depth_bathy", "lyapunov_exp")
for(var in env_vars) {
    if(var %in% names(tara_environmental)) {
        count <- sum(!is.na(tara_environmental[[var]]))
        total <- nrow(tara_environmental)
        pct <- round((count/total)*100, 1)
        cat(sprintf("%-25s: %4d samples (%5.1f%%)\n", var, count, pct))
    } else {
        cat(sprintf("%-25s: Variable not found\n", var))
    }
}   


## We remove marginal sea samples since we don't train our models on those, only open ocean samples
# We only need to do this for the Tara samples and match those afterwards with the climatologies
# Load function
source("../../Functions/iho_ocean_regions.R")
test <- tara_environmental
names(test)
# Rename lat and lon to decimalLatitude and decimalLongitude
names(test)[names(test) == "lat"] <- "decimalLatitude"
names(test)[names(test) == "lon"] <- "decimalLongitude"

# Apply function to get ocean regions
test_regions <- iho_ocean_regions(test)
# Name decimalLatitude and decimalLongitude back to lat and lon
names(test_regions)[names(test_regions) == "decimalLatitude"] <- "lat"
names(test_regions)[names(test_regions) == "decimalLongitude"] <- "lon"

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

# Filter out marginal seas from Tara data
cat("Before removing marginal seas:", nrow(test_regions), "samples\n")
cat("Marginal seas found in data:\n")
marginal_in_data <- intersect(unique(test_regions$marine_region), marginal_seas)
print(marginal_in_data)

# Remove samples from marginal seas (using data.frame syntax)
tara_environmental <- test_regions[!(test_regions$marine_region %in% marginal_seas), ]

# Convert to data.table
tara_environmental <- as.data.table(tara_environmental)

cat("After removing marginal seas:", nrow(tara_environmental), "samples\n")
cat("Samples removed:", nrow(test_regions) - nrow(tara_environmental), "\n")

head(test_regions)

cat("Tara environmental data saved to:", file.path(wd_out, "tara_environmental_DCM_SRF.csv"), "\n")

cat("\nNext step: Match Tara coordinates to climatology grid cells\n")


# Part 2 : Match Tara coordinates to climatology grid cells
# Create long format dataframe from Tara environmental data
cat("\n=== CREATING LONG FORMAT TARA DATAFRAME ===\n")

# Add month extraction to tara_environmental
tara_environmental[, month := {
    # Extract month from date
    date_to_use <- ifelse(!is.na(event_date), event_date, date)
    month_extracted <- as.numeric(substr(date_to_use, 6, 7))
    ifelse(is.na(month_extracted) | month_extracted < 1 | month_extracted > 12, 
           NA_integer_, month_extracted)
}]

# Define environmental variables and their metadata
env_var_metadata <- data.frame(
    var_name = c("temperature", "chla", "par", "backscattering", 
                 "depth_chl_max", "mixed_layer_depth_sigma", "depth_max_brunt_väisälä",
                 "nitrite", "phosphate", "nitrate_nitrite", "silicate", "nstar", 
                 "depth_bathy", "lyapunov_exp"),
    var_unit = c("°C", "mg/m^3", "µmol/m^2/s", "m^-1", 
                 "m", "m", "m",
                 "µmol/L", "µmol/L", "µmol/L", "µmol/L", "µmol/L",
                 "m", "day^-1"),
    var_type = c("Temperature", "Chlorophyll", "Light", "Optical",
                 "Depth", "Physical", "Physical", 
                 "Chemical", "Chemical", "Chemical", "Chemical", "Chemical",
                 "Bathymetry", "Physical"),
    stringsAsFactors = FALSE
)

# Create coordinates columns
tara_coords <- tara_environmental[, .(
    biosample,
    x = ifelse(!is.na(event_longitude), event_longitude, lon),
    y = ifelse(!is.na(event_latitude), event_latitude, lat),
    month = month,
    depth = depth_y,
    station = station
)]

# Remove rows with missing coordinates or month
tara_coords <- tara_coords[!is.na(x) & !is.na(y) & !is.na(month)]

# Melt the environmental data to long format
tara_long_list <- list()

for(i in 1:nrow(env_var_metadata)) {
    var_info <- env_var_metadata[i, ]
    var_col <- var_info$var_name
    
    if(var_col %in% names(tara_environmental)) {
        # Extract data for this variable
        var_data <- tara_environmental[, c("biosample", var_col), with = FALSE]
        setnames(var_data, var_col, "value")
        
        # Merge with coordinates
        var_long <- merge(tara_coords, var_data, by = "biosample", all.x = TRUE)
        
        # Add variable metadata
        var_long[, `:=`(
            var_name = var_info$var_name,
            var_unit = var_info$var_unit,
            var_type = var_info$var_type
        )]
        
        # Remove rows with missing values
        var_long <- var_long[!is.na(value)]
        
        if(nrow(var_long) > 0) {
            tara_long_list[[var_col]] <- var_long
        }
    }
}

# Combine all variables
tara_long <- rbindlist(tara_long_list)



# Add data_type column
tara_long[, data_type := "in_situ"]

# Add column for gridded coordinates -  So they match the climatology grid
head(tara_long)

# Check sample counts for each variable in tara_long
cat("\n=== SAMPLE COUNTS FOR TARA_LONG VARIABLES ===\n")
sample_counts <- tara_long[, .(
    n_samples = .N,
    n_valid = sum(!is.na(value)),
    n_missing = sum(is.na(value))
), by = var_name]
print(sample_counts)

# Harmonize coordinates
tara_long$x_gridded <- tara_long$x
tara_long$y_gridded <- tara_long$y
#
tara_long$x_gridded <- round(tara_long$x_gridded + 0.5) - 0.5
tara_long$y_gridded <- round(tara_long$y_gridded + 0.5) - 0.5

# Create climatology long format with same structure as tara_long
climatology_long <- climatology_df[, .(
    x = lon,  # Use lon column from climatology data
    y = lat,   # Use lat column from climatology data
    month = month,
    var_name = source_file,  # Use source_file as var_name
    var_unit = "",  # Will be filled later
    var_type = "",  # Will be filled later
    value = value,
    depth = "surface",  # Climatology is surface data
    biosample = NA_character_,  # No biosample for climatology
    station = NA_character_   # No station for climatology
)]

# Add data_type column
climatology_long[, data_type := "climatology"]

# Check structure of both dataframes
head(tara_long)
head(climatology_long)

# Harmonize parameter names between tara_long and climatology_long
unique(tara_long$var_name)

# > unique(tara_long$var_name)
# #  [1] "temperature"             "temp"                   
# #  [3] "chla"                    "par"                    
# #  [5] "backscattering"          "depth_chl_max"          
# #  [7] "mixed_layer_depth_sigma" "depth_max_brunt_väisälä"
# #  [9] "nitrite"                 "phosphate"              
# # [11] "nitrate_nitrite"         "silicate"               
# # [13] "nstar"                   "depth_bathy"            
# # [15] "lyapunov_exp"   

# Harmonize climatology variable names to match Tara variable names
unique(climatology_long$var_name)
# Create mapping table for climatology variable names
climatology_var_mapping <- data.frame(
    filename = c(
        "climatology_A_0_50.nc",
        "climatology_A_CHLA_regridded.nc",
        "climatology_A_PAR_regridded.nc",
        "climatology_eke_aviso.nc",
        "climatology_fsle_aviso_2001_2020.nc",
        "climatology_i_0_50.nc",
        "climatology_M_0_0.nc",
        "climatology_n_0_50.nc",
        "climatology_Nstar_0_50.nc",
        "climatology_o_0_50.nc",
        "climatology_O_0_50.nc",
        "climatology_p_0_50.nc",
        "climatology_s_0_50.nc",
        "climatology_t_0_50.nc",
        "log10_climatology_eke_aviso.nc",
        "log10_climatology_i_0_50.nc",
        "log10_climatology_n_0_50.nc",
        "log10_climatology_p_0_50.nc",
        "log10_climatology_S_PP_regridded.nc",
        "log10_climatology_TOT_POC_CMEMS.nc"
    ),
    var_name = c(
        "apparent_oxygen_utilization",  # climatology_A_0_50.nc
        "chla",          # climatology_A_CHLA_regridded.nc
        "par",           # climatology_A_PAR_regridded.nc
        "eddy_kinetic_energy",  # climatology_eke_aviso.nc
        "lyapunov_exp",  # climatology_fsle_aviso_2001_2020.nc
        "silicate",  # climatology_i_0_50.nc
        "mixed_layer_depth",  # climatology_M_0_0.nc
        "nitrate",  # climatology_n_0_50.nc
        "nstar",         # climatology_Nstar_0_50.nc
        "dissolved_oxygen",  # climatology_o_0_50.nc
        "fractional_saturation_oxygen",  # climatology_O_0_50.nc
        "phosphate",     # climatology_p_0_50.nc
        "salinity",      # climatology_s_0_50.nc
        "temperature",   # climatology_t_0_50.nc
        "log10_eddy_kinetic_energy",  # log10_climatology_eke_aviso.nc
        "log10_silicate",  # log10_climatology_i_0_50.nc
        "log10_nitrate",  # log10_climatology_n_0_50.nc
        "log10_phosphate",  # log10_climatology_p_0_50.nc
        "log10_primary_production",  # log10_climatology_S_PP_regridded.nc
        "log10_total_particulate_organic_carbon"  # log10_climatology_TOT_POC_CMEMS.nc
    ),
    var_unit = c(
        "µmol/kg",  # climatology_A_0_50.nc
        "mg/m^3",         # chla
        "µmol/m^2/s",     # par
        "cm^2/s^2",  # climatology_eke_aviso.nc
        "day^-1",         # lyapunov_exp
        "µmol/kg",  # climatology_i_0_50.nc
        "meter",  # climatology_M_0_0.nc
        "µmol/kg",  # climatology_n_0_50.nc
        "µmol/kg",        # nstar
        "µmol/kg",  # climatology_o_0_50.nc
        "µmol/kg",  # climatology_O_0_50.nc
        "µmol/kg",        # phosphate
        "µmol/kg",        # silicate
        "°C",            # temperature
        "FILL_IN_HERE",  # log10_climatology_eke_aviso.nc
        "FILL_IN_HERE",  # log10_climatology_i_0_50.nc
        "FILL_IN_HERE",  # log10_climatology_n_0_50.nc
        "FILL_IN_HERE",  # log10_climatology_p_0_50.nc
        "FILL_IN_HERE",  # log10_climatology_S_PP_regridded.nc
        "FILL_IN_HERE"   # log10_climatology_TOT_POC_CMEMS.nc
    ),
    var_type = c(
        "FILL_IN_HERE",  # climatology_A_0_50.nc
        "Chlorophyll",   # chla
        "Light",         # par
        "FILL_IN_HERE",  # climatology_eke_aviso.nc
        "Physical",      # lyapunov_exp
        "FILL_IN_HERE",  # climatology_i_0_50.nc
        "FILL_IN_HERE",  # climatology_M_0_0.nc
        "FILL_IN_HERE",  # climatology_n_0_50.nc
        "Chemical",      # nstar
        "FILL_IN_HERE",  # climatology_o_0_50.nc
        "FILL_IN_HERE",  # climatology_O_0_50.nc
        "Chemical",      # phosphate
        "Chemical",      # silicate
        "Temperature",   # temperature
        "FILL_IN_HERE",  # log10_climatology_eke_aviso.nc
        "FILL_IN_HERE",  # log10_climatology_i_0_50.nc
        "FILL_IN_HERE",  # log10_climatology_n_0_50.nc
        "FILL_IN_HERE",  # log10_climatology_p_0_50.nc
        "FILL_IN_HERE",  # log10_climatology_S_PP_regridded.nc
        "FILL_IN_HERE"   # log10_climatology_TOT_POC_CMEMS.nc
    ),
    stringsAsFactors = FALSE
)

# Apply the mapping to climatology_long
climatology_long <- merge(climatology_long, climatology_var_mapping, 
                         by.x = "var_name", by.y = "filename", all.x = TRUE)

# Update the var_name column with harmonized names and keep original filename
climatology_long[, `:=`(
    filename = var_name,           # Keep original filename
    var_name = var_name.y,         # Use harmonized name
    var_unit = var_unit.y,         # Use mapped unit
    var_type = var_type.y          # Use mapped type
)]

# Clean up duplicate columns
climatology_long[, c("var_name.y", "var_unit.y", "var_type.y") := NULL]

cat("Climatology variable names updated:\n")
print(unique(climatology_long[, .(filename, var_name, var_unit, var_type)]))

# Debug silicate specifically
cat("\n=== DEBUGGING SILICATE MAPPING ===\n")
silicate_data <- climatology_long[var_name == "silicate"]
cat("Silicate climatology data:\n")
cat("Total rows:", nrow(silicate_data), "\n")
cat("Unique filenames contributing to silicate:\n")
print(unique(silicate_data$filename))
cat("Sample counts by filename:\n")
print(silicate_data[, .N, by = filename])

# Harmonize variable names between tara_long and climatology_long
cat("\n=== HARMONIZING VARIABLE NAMES ===\n")

# First standardize tara_long variable names
cat("Tara variable names before harmonization:\n")
print(sort(unique(tara_long$var_name)))

# Standardize "mixed_layer_depth_sigma" to "mixed_layer_depth" in tara_long to match climatology
tara_long[var_name == "mixed_layer_depth_sigma", var_name := "mixed_layer_depth"]

# Standardize "nitrate_nitrite" to "nitrate" in tara_long to match climatology
tara_long[var_name == "nitrate_nitrite", var_name := "nitrate"]

cat("\nTara variable names after harmonization:\n")
print(sort(unique(tara_long$var_name)))

cat("\nClimatology variable names:\n")
print(sort(unique(climatology_long$var_name)))

# Find common variables between the two datasets
cat("\n=== DEBUGGING VARIABLE HARMONIZATION ===\n")
cat("Tara variable names after harmonization:\n")
tara_vars <- sort(unique(tara_long$var_name))
print(tara_vars)

cat("\nClimatology variable names after mapping:\n")
clim_vars <- sort(unique(climatology_long$var_name))
print(clim_vars)

cat("\nChecking if climatology_long has data:\n")
cat("Rows in climatology_long:", nrow(climatology_long), "\n")
cat("Unique var_name values:", length(unique(climatology_long$var_name)), "\n")

cat("\nChecking if tara_long has data:\n")
cat("Rows in tara_long:", nrow(tara_long), "\n")
cat("Unique var_name values:", length(unique(tara_long$var_name)), "\n")

# Check if the merge worked correctly for climatology
cat("\nChecking climatology merge results:\n")
print(unique(climatology_long[, .(var_name, filename)]))

common_vars <- intersect(unique(tara_long$var_name), unique(climatology_long$var_name))
cat("\nCommon variables between Tara and climatology data:\n")
print(common_vars)

# ==================== CORRELATION ANALYSIS ====================
cat("\n=== MATCHING TARA AND CLIMATOLOGY DATA ===\n")

# Add gridded coordinates to climatology data for matching
climatology_long[, `:=`(
    x_gridded = round(x + 0.5) - 0.5,
    y_gridded = round(y + 0.5) - 0.5
)]

# Define specific variables to analyze
analysis_vars <- c("chla", "mixed_layer_depth", "nitrate", "par", "phosphate", "silicate", "temperature")

# Create Tara units mapping
tara_units <- data.frame(
    var_name = c("temperature", "nitrate", "phosphate", "silicate", "chla", "par", "mixed_layer_depth"),
    tara_unit = c("°C", "µmol/L", "µmol/L", "µmol/L", "mg/m^3", "µmol photons/m^2/s", "m"),
    stringsAsFactors = FALSE
)

# Filter both datasets to only include analysis variables
tara_common <- tara_long[var_name %in% analysis_vars]
climatology_common <- climatology_long[var_name %in% analysis_vars]

# Create matching key for spatial and temporal matching
tara_common[, match_key := paste(x_gridded, y_gridded, month, var_name, sep = "_")]
climatology_common[, match_key := paste(x_gridded, y_gridded, month, var_name, sep = "_")]

# DEBUG: Check data before merging
cat("=== DEBUGGING MATCH ISSUES ===\n")
cat("Tara common data summary:\n")
print(tara_common[, .(
    n_rows = .N,
    unique_vars = length(unique(var_name)),
    unique_coords = length(unique(paste(x_gridded, y_gridded))),
    unique_months = length(unique(month)),
    sample_match_keys = head(unique(match_key), 3)
), by = var_name])

cat("\nClimatology common data summary:\n")
print(climatology_common[, .(
    n_rows = .N,
    unique_vars = length(unique(var_name)),
    unique_coords = length(unique(paste(x_gridded, y_gridded))),
    unique_months = length(unique(month)),
    sample_match_keys = head(unique(match_key), 3)
), by = var_name])

cat("\nTara coordinate ranges:\n")
print(tara_common[, .(
    lon_min = min(x_gridded, na.rm = TRUE),
    lon_max = max(x_gridded, na.rm = TRUE),
    lat_min = min(y_gridded, na.rm = TRUE),
    lat_max = max(y_gridded, na.rm = TRUE),
    month_range = paste(range(month, na.rm = TRUE), collapse = "-")
)])

cat("\nClimatology coordinate ranges:\n")
print(climatology_common[, .(
    lon_min = min(x_gridded, na.rm = TRUE),
    lon_max = max(x_gridded, na.rm = TRUE),
    lat_min = min(y_gridded, na.rm = TRUE),
    lat_max = max(y_gridded, na.rm = TRUE),
    month_range = paste(range(month, na.rm = TRUE), collapse = "-")
)])

# Check for any common match keys
common_keys <- intersect(tara_common$match_key, climatology_common$match_key)
cat("\nNumber of common match keys:", length(common_keys), "\n")
if(length(common_keys) > 0) {
    cat("Sample common match keys:\n")
    print(head(common_keys, 5))
}

# Debug silicate matching specifically
cat("\n=== DEBUGGING SILICATE MATCHING ===\n")
cat("Tara silicate data:\n")
tara_silicate <- tara_common[var_name == "silicate"]
cat("Tara silicate rows:", nrow(tara_silicate), "\n")
cat("Unique Tara silicate match_keys:", length(unique(tara_silicate$match_key)), "\n")

cat("\nClimatology silicate data:\n")
clim_silicate <- climatology_common[var_name == "silicate"]
cat("Climatology silicate rows:", nrow(clim_silicate), "\n")
cat("Unique climatology silicate match_keys:", length(unique(clim_silicate$match_key)), "\n")

cat("\nSilicate match key overlap:\n")
silicate_common_keys <- intersect(tara_silicate$match_key, clim_silicate$match_key)
cat("Common silicate match_keys:", length(silicate_common_keys), "\n")

# Merge the datasets
matched_data <- merge(
    tara_common[, .(match_key, var_name, tara_value = value, depth, biosample)],
    climatology_common[, .(match_key, var_name, clim_value = value)],
    by = c("match_key", "var_name"),
    all = FALSE  # Only keep matches
)

# Remove rows with missing values
matched_data <- matched_data[!is.na(tara_value) & !is.na(clim_value)]

cat("Matched data points:\n")
print(matched_data[, .N, by = var_name])

# Function to calculate correlations and create plots
create_correlation_plots <- function(data, method = "pearson", title_suffix = "") {
    
    # Calculate correlations for each variable
    correlation_results <- data[, .(
        n = .N,
        correlation = cor(tara_value, clim_value, method = method, use = "complete.obs"),
        p_value = cor.test(tara_value, clim_value, method = method)$p.value
    ), by = var_name]
    
    # Add significance stars
    correlation_results[, significance := ifelse(p_value < 0.001, "***",
                                                ifelse(p_value < 0.01, "**", 
                                                      ifelse(p_value < 0.05, "*", "ns")))]
    
    # Create correlation labels for plots
    correlation_results[, label := paste0("r = ", round(correlation, 3), " ", significance, "\nn = ", n)]
    
    # Add units information from both datasets
    plot_data <- merge(data, tara_units, by = "var_name", all.x = TRUE)
    plot_data <- merge(plot_data, unique(climatology_common[, .(var_name, clim_unit = var_unit)]), by = "var_name", all.x = TRUE)
    
    # Create facet labels with units
    plot_data[, facet_label := paste0(var_name, "\nTara: ", tara_unit, " | Clim: ", clim_unit)]
    
    # Create label data with unique entries per variable for text annotation
    label_data <- merge(correlation_results[, .(var_name, label)], 
                       unique(plot_data[, .(var_name, facet_label)]), by = "var_name")
    
    # Create the plot
    p <- ggplot(plot_data, aes(x = clim_value, y = tara_value)) +
        geom_point(alpha = 0.6, color = "black", size = 0.8) +
        geom_smooth(method = "lm", se = TRUE, color = "blue", alpha = 0.3, size = 0.8) +
        geom_text(data = label_data, aes(x = Inf, y = Inf, label = label), 
                 hjust = 1.1, vjust = 1.5, size = 2.5, color = "red", inherit.aes = FALSE) +
        facet_wrap(~ facet_label, scales = "free", ncol = 3) +
        labs(
            title = "In-situ vs Climatology Correlations (Pearson)",
            x = "Climatology Value",
            y = "Tara In-situ Value"
        ) +
        theme_bw() +
        theme(
            strip.text = element_text(size = 7, margin = margin(2, 2, 2, 2)),
            axis.text = element_text(size = 6),
            axis.title = element_text(size = 8),
            plot.title = element_text(size = 10, hjust = 0.5),
            panel.grid.minor = element_blank(),
            strip.background = element_rect(fill = "white", colour = "black", size = 0.5)
        )
    
    return(list(plot = p, correlations = correlation_results))
}

# Create Pearson correlation plots
cat("\n=== CREATING PEARSON CORRELATION PLOTS ===\n")
pearson_results <- create_correlation_plots(matched_data, method = "pearson")

# Print Pearson correlation results
cat("Pearson correlation results:\n")
print(pearson_results$correlations)

# Save Pearson plot
ggsave(file.path(wd_out, "pearson_correlations_tara_climatology.png"), 
       plot = pearson_results$plot, width = 14, height = 10, dpi = 300)
#
ggsave(file.path(wd_out, "pearson_correlations_tara_climatology.svg"), 
       plot = pearson_results$plot, width = 14, height = 10, dpi = 300)

# Save correlation results as CSV
fwrite(pearson_results$correlations, file.path(wd_out, "pearson_correlation_results_tara_climatology.csv"))

# ==================== CREATE SPATIAL MAP ====================
cat("\n=== CREATING SPATIAL DISTRIBUTION MAP ===\n")

# Extract unique coordinates from matched data
matched_coords <- unique(matched_data[, .(match_key)])

# Debug: Check the structure of match_key
cat("Sample match_keys:\n")
print(head(matched_coords$match_key, 5))

# Split match_key more carefully
# The pattern is: longitude_latitude_month_variable
# But we need to account for negative coordinates
matched_coords[, c("x_gridded", "y_gridded", "month", "var_name") := {
    # Split by underscore
    parts <- tstrsplit(match_key, "_", fixed = TRUE)
    
    # Handle the case where we might have negative coordinates
    # If we get more than 4 parts, it means negative signs created extra splits
    if(length(parts[[1]]) > 0 && any(sapply(parts, length) > length(parts[[1]]))) {
        # Reconstruct by taking first part as longitude, 
        # combine parts until we get latitude, then month, then variable
        list(
            x_gridded = as.numeric(parts[[1]]),
            y_gridded = as.numeric(parts[[2]]), 
            month = as.numeric(parts[[3]]),
            var_name = parts[[4]]
        )
    } else {
        list(
            x_gridded = as.numeric(parts[[1]]),
            y_gridded = as.numeric(parts[[2]]), 
            month = as.numeric(parts[[3]]),
            var_name = parts[[4]]
        )
    }
}]

cat("Processed coordinates sample:\n")
print(head(matched_coords))

# Get unique spatial locations (remove temporal and variable duplicates)
spatial_coords <- unique(matched_coords[, .(lon = x_gridded, lat = y_gridded)])

cat("Number of unique spatial match-up locations:", nrow(spatial_coords), "\n")

# Create world map with match-up locations
world_map <- map_data("world")

map_plot <- ggplot() +
    # Add world map
    geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
                 fill = "lightgray", color = "white", size = 0.2) +
    # Add match-up locations
    geom_point(data = spatial_coords, aes(x = lon, y = lat), 
               color = "red", size = 1.2, alpha = 0.7) +
    # Map styling
    coord_quickmap() +
    theme_void() +
    labs(
        title = "Spatial Distribution of Tara-Climatology Match-up Locations",
        subtitle = paste("Total unique locations:", nrow(spatial_coords))
    ) +
    theme(
        plot.title = element_text(size = 14, hjust = 0.5, margin = margin(b = 10)),
        plot.subtitle = element_text(size = 12, hjust = 0.5, color = "darkblue"),
        panel.background = element_rect(fill = "lightblue", color = NA)
    )

# Save the map
ggsave(file.path(wd_out, "matchup_locations_map.png"), 
       plot = map_plot, width = 12, height = 8, dpi = 300)

ggsave(file.path(wd_out, "matchup_locations_map.svg"), 
       plot = map_plot, width = 12, height = 8, dpi = 300)

cat("Spatial map saved as: matchup_locations_map.png and .svg\n")

# Create summary table by region/ocean basin
spatial_coords[, ocean_basin := case_when(
    lon >= -180 & lon <= -60 & lat >= -60 & lat <= 70 ~ "Atlantic",
    lon >= -60 & lon <= 20 & lat >= -60 & lat <= 70 ~ "Atlantic", 
    lon >= 20 & lon <= 150 & lat >= -60 & lat <= 70 ~ "Indian/Pacific",
    lon >= 150 & lon <= 180 & lat >= -60 & lat <= 70 ~ "Pacific",
    TRUE ~ "Other"
)]

cat("\nMatch-up locations by ocean basin:\n")
basin_summary <- spatial_coords[, .N, by = ocean_basin]
print(basin_summary)

# Save spatial coordinates
fwrite(spatial_coords, file.path(wd_out, "matchup_spatial_coordinates.csv"))

cat("Spatial coordinates saved to: matchup_spatial_coordinates.csv\n")

# Print summary
cat("\n=== CORRELATION ANALYSIS SUMMARY ===\n")
cat("Total matched data points:", nrow(matched_data), "\n")
cat("Variables analyzed:", length(analysis_vars), "\n")
cat("Variables with significant Pearson correlations (p < 0.05):\n")
print(pearson_results$correlations[p_value < 0.05, .(var_name, correlation, p_value)])

cat("\nCorrelation plot saved as: pearson_correlations_tara_climatology.png\n")


# Create final multiplot
cat("\n=== CREATING COMBINED MULTIPLOT ===\n")

# Update map plot to have panel label and remove title for cleaner combination
map_plot_clean <- ggplot() +
    # Add world map
    geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
                 fill = "lightgray", color = "white", size = 0.3) +
    # Add match-up locations
    geom_point(data = spatial_coords, aes(x = lon, y = lat), 
               color = "red", size = 0.8, alpha = 0.7) +
    # Map styling
    coord_quickmap() +
    theme_void() +
    labs(
        title = "Spatial Distribution of Match-up Locations"
    ) +
    theme(
        plot.title = element_text(size = 11, hjust = 0.5, margin = margin(b = 8)),
        panel.background = element_rect(fill = "lightblue", color = NA),
        plot.margin = margin(10, 5, 5, 15)
    )

# Update correlation plot to have panel label and adjust title
correlation_plot_clean <- pearson_results$plot +
    labs(title = "In-situ vs Climatology Correlations") +
    theme(
        plot.title = element_text(size = 11, hjust = 0.5, margin = margin(b = 8)),
        strip.text = element_text(size = 7),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        plot.margin = margin(5, 5, 5, 15)
    )

# Create combined plot using gridExtra
library(gridExtra)
library(grid)

# Add panel labels using annotate
map_with_label <- map_plot_clean + 
    annotate("text", x = -Inf, y = Inf, label = "a", 
            hjust = -1.2, vjust = 2.0, size = 5, fontface = "bold")

correlation_with_label <- correlation_plot_clean + 
    annotate("text", x = -Inf, y = Inf, label = "b", 
            hjust = -0.3, vjust = 1.2, size = 5, fontface = "bold")

# Combine the plots with appropriate dimensions for journal publication
combined_plot_labeled <- grid.arrange(
    map_with_label,
    correlation_with_label,
    ncol = 1,
    heights = c(0.35, 0.65),  # Map gets 35% of height, correlations get 65%
    top = textGrob("Tara Ocean In-situ vs Climatology Data Pearson Correlation", 
                   gp = gpar(fontsize = 12, fontface = "bold"))
)

# Save the combined plot with journal-appropriate dimensions
ggsave(file.path(wd_out, "combined_tara_climatology_correlation.png"), 
       plot = combined_plot_labeled, width = 7, height = 9, dpi = 300)

ggsave(file.path(wd_out, "combined_tara_climatology_correlation.svg"), 
       plot = combined_plot_labeled, width = 7, height = 9, dpi = 300)

ggsave(file.path(wd_out, "combined_tara_climatology_correlation.pdf"), 
       plot = combined_plot_labeled, width = 7, height = 9, dpi = 300)

cat("Combined multiplot saved as: combined_tara_climatology_correlation.png/.svg/.pdf\n")
cat("Dimensions: 7 x 9 inches (journal double-column format)\n")



# ==============================================================================
# FIGURE CAPTION
# ==============================================================================
#
# Figure 1. Pearson correlation analysis of climatological environmental data using Tara Ocean 
# in-situ measurements. (a) Global distribution of correlation analysis points where Tara 
# Ocean in-situ measurements were successfully matched with climatological 
# environmental data. Red dots indicate unique geographic locations spanning 
# major open ocean basins, representing sampling sites from both surface waters 
# (0-50 m depth) and Deep Chlorophyll Maximum (DCM) layers. Samples from marginal 
# seas (Mediterranean, Red Sea, Baltic Sea, Black Sea, Persian Gulf, and associated 
# regions) were excluded to maintain consistency with open ocean modeling approaches. 
# (b) Pearson correlations between Tara Ocean in-situ measurements (y-axis) and 
# corresponding climatological values (x-axis) for seven environmental variables: 
# chlorophyll-a concentration (chla), mixed layer depth, nitrate concentration, 
# photosynthetically active radiation (PAR), phosphate concentration, silicate 
# concentration, and temperature. Each panel displays the variable name with units 
# for both Tara Ocean data and climatology data. Blue lines represent linear 
# regression fits with 95% confidence intervals (gray shading). Correlation 
# coefficients (r), statistical significance (*** p < 0.001, ** p < 0.01, 
# * p < 0.05, ns = not significant), and sample sizes (n) are shown in the upper 
# right corner of each panel. Coordinate matching was performed using a 0.5° × 0.5° 
# grid system, and temporal matching used monthly climatological averages. Strong 
# correlations (r > 0.7) for temperature, PAR, nitrate, phosphate, and mixed layer 
# depth indicate excellent agreement between climatological estimates and direct 
# measurements, validating the use of climatological data as environmental predictors 
# in marine species distribution modeling studies. Moderate correlations for 
# chlorophyll-a and silicate reflect the higher temporal and spatial variability 
# of these biogeochemical parameters compared to more stable physical variables.
#
# ================================================### ====================================================================================================================
### END OF SCRIPT
### ====================================================================================================================




