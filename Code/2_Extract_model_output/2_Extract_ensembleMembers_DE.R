
##====================================================================================================
## Ensemble Members Extraction for Prokaryotic Species Distribution Models
##
## This R script extracts and processes ensemble members from CEPHALOPOD model outputs for both
## class-level and genus-level prokaryotic taxa. The analysis processes model predictions across
## multiple bootstraps and monthly time steps, creating comprehensive datasets for latitudinal
## diversity gradient analysis. The script handles parallel processing for computational efficiency
## and generates both individual ensemble member data and summary statistics.
##
## Author:      Dominic Eriksson
##              Environmental Physics Group, UP
##              ETH Zürich, Switzerland
## Contact:     deriksson@ethz.ch
## Date:        February 27th, 2026
## Affiliation: ETH Zürich, Environmental Physics Group, UP
##
## Input files:
## - ./Code/Cephalopod_output/Log/chunks*/*/MODEL.RData (Class-level model outputs)
## - ./Code/Cephalopod_output/Non_log_genusLevel/chunks*/*/MODEL.RData (Genus-level model outputs)
## - ./Code/Cephalopod_output/*/CALL.RData (Model configuration files)
## - ./Code/Cephalopod_output/*/*/QUERY.RData (Query parameter files)
##
## Output files:
## - *_ensembleMembers_Log_v2.RData/*_ensembleMembers_Non_log_v2.RData (Individual ensemble member rasters)
## - Statistics_acrossEnsembleMembers_selected_clades_v3.csv (Complete ensemble member statistics)
## - LDG_Statistics_acrossEnsembleMembers_selected_clades_v3.csv (Latitudinal diversity gradient summaries)
##
## Strategy:
## 1. Configure parallel processing environment for computational efficiency
## 2. Process class-level ensemble members from Log-transformed predictors
## 3. Extract model predictions across bootstraps (10) and months (12) for each algorithm
## 4. Apply quality control recommendations and filter validated algorithms
## 5. Convert model arrays to georeferenced raster stacks with WGS84 projection
## 6. Process genus-level ensemble members with similar workflow
## 7. Aggregate ensemble member data into comprehensive dataframes
## 8. Calculate summary statistics (mean, SD, median, IQR) by latitude and month
## 9. Generate latitudinal diversity gradient visualizations for validation
##
## Required R packages:
## - abind (1.4.5): Multidimensional array operations and binding
## - doParallel (1.0.17): Parallel processing backend for foreach loops
## - foreach (1.5.2): Enhanced looping constructs for parallel execution
## - raster (3.5.15): Spatial raster data manipulation and analysis
## - dplyr (1.0.9): Data manipulation, filtering, and transformation operations
## - sp (1.4.7): Spatial data classes and coordinate reference systems
## - tidyr (1.1.3): Data tidying and reshaping for analysis workflows
## - data.table (1.14.2): Fast and efficient data manipulation and file I/O
## - ggplot2 (3.3.6): Advanced data visualization for quality control plots
##====================================================================================================

# Clear working space
rm(list = ls())

### ==============================================================================================
### Class Level Ensemble Members Extraction
### ==============================================================================================

### Load Required Libraries
library(abind)      # For combining multidimensional arrays
library(doParallel) # For parallel processing
library(foreach)    # For parallel loops
library(raster)     # For raster data manipulation
library(dplyr)      # For data manipulation
library(sp)         # For spatial data handling
library(tidyr)      # For data tidying and reshaping

# Load specific functions
setwd("../../CEPHALOPOD")
source(file = "../../CEPHALOPOD/code/00_config.R")

# Choose between log and non-log transformed predictors
trans <- c("Non_log", "Log")
t <- 2  # Using Log-transformed predictors

### Set Parameters and Directories
depth_level <- c("MLD", "0_10Meters")  # Define depth levels: Mixed Layer Depth (MLD) or 0-10 meters
d <- 1                                 # Specify the depth level index

# Directories
project_wd <- "./Code/Cephalopod_output"
wd_out <- "./Code/2_Extract_model_outputs/3_Output/"
# Create output directory if it doesn't exist
if (!dir.exists(wd_out)) {
  dir.create(wd_out, recursive = TRUE)
}

# Get folder names to loop over
FOLDER_NAMES <- list.files(paste0(project_wd, "/", trans[t],"/"), full.names = TRUE, pattern = "chunk")

### Parallelization Setup
n.cores <- 20
cl <- makeCluster(n.cores, outfile = "")
registerDoParallel(cl)

### Process Each Folder and Subfolder
for (f in seq_along(FOLDER_NAMES)) {
# for (f in 5:length(FOLDER_NAMES)) {
  
  # Get subfolder names
  SUBFOLDER_NAMES <- list.files(FOLDER_NAMES[f], full.names = TRUE, pattern = "1000|5000")
  
  # Parallel loop over subfolders
  foreach::foreach(s = seq_along(SUBFOLDER_NAMES), .packages = c("abind", "sp", "raster", "dplyr")) %dopar% {
    
    # Print progress
    print(paste("Processing folder", f, "of", length(FOLDER_NAMES), "and subfolder", s, "of", length(SUBFOLDER_NAMES)))
    
    # Load parameter and query data
    SUBFOLDER_NAME <- SUBFOLDER_NAMES[s]
    CALL <- get(load(paste0(FOLDER_NAMES[f], "/CALL.RData")))
    QUERY <- get(load(paste0(SUBFOLDER_NAME, "/QUERY.RData")))
    
    # Check for MODEL.RData and load if exists
    model_file_path <- paste0(SUBFOLDER_NAME, "/MODEL.RData")
    if (!file.exists(model_file_path)) {
      message("MODEL.RData not found in ", SUBFOLDER_NAME, ". Skipping.")
      return(NULL)
    } else {
      MODEL <- get(load(model_file_path))
    }
    
    # Check if projections exist
    if ((length(MODEL$MODEL_LIST) == 0) & CALL$FAST == TRUE) {
      message("No validated algorithms for projections.")
      return(NULL)
    }
    
    ### Set Initial Parameters and Masks
    r0 <- CALL$ENV_DATA[[1]][[1]]
    # land <- raster::setValues(r0, ifelse(is.na(r0), 9999, NA)) # Not working anymore
    # Extract values from the raster
    r0_values <- getValues(r0)
    # Check for NA values and replace as desired
    new_values <- ifelse(is.na(r0_values), 9999, NA)
    # Set the modified values back to the raster
    land <- setValues(r0, new_values)
    
    ### Filter Recommendations
    rec <- qc_recommandations(QUERY = QUERY, MODEL = MODEL, DATA_TYPE = CALL$DATA_TYPE)
    MODEL$recommandations <- rec
    filtered_rows <- rec %>% filter(PRE_VIP == 1, FIT == 1, CUM_VIP == 1, DEV == 1)
    if (nrow(filtered_rows) == 0) {
      message("No rows meet criteria. Skipping.")
      return(NULL)
    }
    loop_over <- rownames(filtered_rows)
    
    ### Build Ensembles
    list_array <- lapply(loop_over, function(i) {
      if (CALL$DATA_TYPE == "proportions") {
        MODEL[["MBTR"]][["proj"]][["y_hat"]][,,i,]
      } else {
        MODEL[[i]][["proj"]][["y_hat"]]
      }
    })
    
   # Name list_array
   names(list_array) <- loop_over
   
   list_algorithms <- list()
    for( loop in loop_over){   

        # Extract model data
        d <- list_array[[loop]]

            # Initialize a list to store raster layers for each bootstrap-month combination
            raster_layers <- list()

            # Loop through each bootstrap
            for (bootstrap in 1:dim(d)[2]) {  # Loop over the 10 bootstraps
              # Loop through each month
              for (month in 1:dim(d)[3]) {  # Loop over the 12 months
                  # Extract the data for the current bootstrap and month
                  layer_data <- d[, bootstrap, month]
                  
                  # Create a raster layer using the template r0
                  raster_layer <- setValues(r0, layer_data)
                  
                  # Assign a meaningful name to the raster layer
                  layer_name <- paste0("Bootstrap_", bootstrap, "_Month_", month)
                  names(raster_layer) <- layer_name
                  
                  # Store the raster layer in the list
                  raster_layers[[layer_name]] <- raster_layer
              }
            }

        # Stack all in on
        r <- raster::stack(raster_layers)

        # Replace all negative values with zero in the raster stack
        r[r < 0] <- 0

        # Set the projection to WGS84
        wgs84 <- sp::CRS("+proj=longlat +datum=WGS84 +no_defs")
        raster::projection(r) <- wgs84
        
        # Save in the list
        list_algorithms[[loop]] <- r
    }

# Save list object 
save(
    list_algorithms,   
    file = paste0(wd_out, basename(SUBFOLDER_NAME), "_ensembleMembers_", trans[t], "_v2.RData")
)
    
    # Print progress
    print(paste("Finished folder", f, "subfolder", s))
  }
}

### Stop Cluster
parallel::stopCluster(cl)

### ===============================================================================
### Create one big dataframe with all ensemble members
### ===============================================================================

# Get filenames
# fnames <- list.files("/home/deriksson/Projects/Prokaryotes_LDG_v2/Code/2_Extract_model_outputs/4_Output", full.names = TRUE) 4_Output?
fnames <- list.files("./Code/2_Extract_model_outputs/3_Output", full.names = TRUE)

# Subset filenames for Richness_5000
fnames <- fnames[grepl("Richness_5000", fnames)]
fnames <- fnames[grepl("_v2", fnames)]

# Define selected clades
selected_clades <- c(
  "d__Archaea_Richness_5000",
  "d__Bacteria_Richness_5000",
  "d__Archaea_p__Thermoproteota_c__Nitrososphaeria_Richness_5000",
  "d__Bacteria_p__Fibrobacterota_c__Fibrobacteria_Richness_5000",
  "d__Bacteria_p__Bacteroidota_c__Bacteroidia_Richness_5000",
  "d__Archaea_p__Thermoplasmatota_c__Poseidoniia_Richness_5000",
  "d__Bacteria_p__Verrucomicrobiota_c__Verrucomicrobiae_Richness_5000",
  "d__Bacteria_p__Desulfobacterota_D_c__UBA1144_Richness_5000",
  "d__Bacteria_p__Marinisomatota_c__Marinisomatia_Richness_5000",
  "d__Bacteria_p__Chloroflexota_c__Dehalococcoidia_Richness_5000",
  "d__Bacteria_p__Pseudomonadota_c__Gammaproteobacteria_Richness_5000",
  "d__Bacteria_p__SAR324_c__SAR324_Richness_5000",
  "d__Bacteria_p__Chlamydiota_Richness_5000",
  "d__Bacteria_p__Bacteroidota_c__Rhodothermia_Richness_5000",
  "d__Bacteria_p__Verrucomicrobiota_c__Kiritimatiellia_Richness_5000",
  "d__Bacteria_p__Actinomycetota_c__Acidimicrobiia_Richness_5000",
  "d__Bacteria_p__Cyanobacteriota_c__Cyanobacteriia_Richness_5000",
  "d__Bacteria_p__Pseudomonadota_c__Alphaproteobacteria_Richness_5000",
  "d__Archaea_p__Asgardarchaeota_c__Heimdallarchaeia_Richness_5000",
  "d__Bacteria_p__Bdellovibrionota_c__Bacteriovoracia_Richness_5000",
  "d__Bacteria_p__Myxococcota_c__UBA796_Richness_5000",
  "d__Bacteria_p__Myxococcota_c__UBA4151_Richness_5000",
  "d__Bacteria_p__Myxococcota_c__XYA12-FULL-58-9_Richness_5000"
)

# Filter filenames based on selected clades
fnames <- fnames[grepl(paste(selected_clades, collapse = "|"), fnames)]

# Loop through files
final_results <- list()
for( f in seq_along(fnames) ){

    list_algorithms <- get(load(fnames[f]))


            list_df <- list()
            for(i in seq_along(list_algorithms)){

            # Extract data
            d <- as.data.frame( list_algorithms[[i]], xy = TRUE) 

            # Convert to long format
            df_long <- d %>%
                tidyr::pivot_longer(
                cols = starts_with("Bootstrap"),  # Select all columns starting with "Bootstrap"
                names_to = c("Bootstrap", "Month"),  # Split column names into "Bootstrap" and "Month"
                names_pattern = "Bootstrap_(\\d+)_Month_(\\d+)",  # Regex to extract bootstrap and month
                values_to = "richness"  # Name of the column for values
                )


            # Add algorithm info
            df_long$algorithms <- names(list_algorithms)[i]

            # Save in list
            list_df[[i]] <- df_long

            }

          # Combine all data frames into one
          df_long <- bind_rows(list_df)
          head(df_long)

        # Add clade name
        clade_name <- gsub(".RData", "", basename(fnames[f]) )
        df_long$clade <- clade_name

        # Print progress
        print(paste("Finished file", f, "of", length(fnames)))

        # Save in list
        final_results[[f]] <- df_long
}

# Merge into one dataframe
df_final <- do.call("rbind", final_results)

# Check dataframe
head(df_final)

# Save complete dataframe not LDG statistics
library(data.table)
fwrite(
  df_final,
  file = paste0("/home/deriksson/Projects/Prokaryotes_LDG_v3/Code/2_Extract_model_outputs/3_Output/", "Statistics_acrossEnsembleMembers_selected_clades_v3.csv"),
  row.names = FALSE
)

# Load library
library(dplyr)

# Compute mean, sd, median, and IQR grouped by x, Month, and clade
summary_stats <- df_final %>%
  group_by(y, Month, clade) %>%
  summarise(
    mean_richness = mean(richness, na.rm = TRUE),  # Compute mean
    sd_richness = sd(richness, na.rm = TRUE),      # Compute standard deviation
    median_richness = median(richness, na.rm = TRUE),  # Compute median
    iqr_richness = IQR(richness, na.rm = TRUE),    # Compute interquartile range
    .groups = "drop"  # Ungroup after summarising
  )

# Check summary statistics
head(summary_stats)

# Save summary stats
write.csv(
  summary_stats,
  file = paste0("/home/deriksson/Projects/Prokaryotes_LDG_v3/Code/2_Extract_model_outputs/3_Output/", "LDG_Statistics_acrossEnsembleMembers_selected_clades_v3.csv"),
  row.names = FALSE

)

# Quick plot
library(ggplot2)
library(dplyr)

# Load data
summary_stats <- read.csv(paste0("./Code/2_Extract_model_outputs/3_Output/", "LDG_Statistics_acrossEnsembleMembers_selected_clades_v3.csv"))

# Filter for June and December
df_june_december <- summary_stats %>%
  filter(Month %in% c("6", "12")) %>%
  mutate(Month = factor(Month, levels = c("6", "12"), labels = c("June", "December")))

# Plot the data
ggplot(df_june_december, aes(x = y, y = mean_richness, color = Month, group = Month)) +
  geom_line(size = 1.2) +  # Line plot for mean richness
  geom_ribbon(aes(ymin = mean_richness - sd_richness, ymax = mean_richness + sd_richness, fill = Month), alpha = 0.2) +  # Ribbon for ± sd
  facet_wrap(~ clade, scales = "free_y") +  # Facet by clade
  labs(
    title = "Latitudinal Gradients for June and December by Clade",
    x = "Latitude",
    y = "Mean Richness",
    color = "Month",
    fill = "Month"
  ) +
  scale_color_manual(values = c("June" = "blue", "December" = "red")) +  # Custom colors for lines
  scale_fill_manual(values = c("June" = "blue", "December" = "red")) +  # Custom colors for ribbons
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10),  # Adjust facet label size
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    legend.position = "top"  # Place legend at the top
  )

### ===============================================================================
### Check, plottting
### ===============================================================================

list_df <- list()
for(i in seq_along(list_algorithms)){

  # Extract data
  d <- as.data.frame( list_algorithms[[i]], xy = TRUE) 

  # Convert to long format
  df_long <- d %>%
    pivot_longer(
      cols = starts_with("Bootstrap"),  # Select all columns starting with "Bootstrap"
      names_to = c("Bootstrap", "Month"),  # Split column names into "Bootstrap" and "Month"
      names_pattern = "Bootstrap_(\\d+)_Month_(\\d+)",  # Regex to extract bootstrap and month
      values_to = "richness"  # Name of the column for values
    )


  # Add algorithm info
  df_long$algorithms <- names(list_algorithms)[i]

  # Save in list
  list_df[[i]] <- df_long

}

# Combine all data frames into one
df_long <- bind_rows(list_df)
head(df_long)

# Aggregate using the median of richness
df_aggregated <- df_long %>%
  group_by(y, Bootstrap, Month, algorithms) %>%
  summarise(
    median_richness = median(richness, na.rm = TRUE),  # Calculate the median, ignoring NA values
    .groups = "drop"  # Ungroup after summarising
  )

# Check the result
head(df_aggregated)

# Check lineplot
library(ggplot2)

# Plot line plots for each bootstrap and algorithm, faceted by month
ggplot(df_aggregated, aes(x = y, y = median_richness, color = algorithms, group = interaction(Bootstrap, algorithms))) +
  geom_line(alpha = 0.7) +  # Add lines with some transparency
  facet_wrap(~ Month, scales = "free_y") +  # Facet by month
  labs(
    title = "Latitudinal Gradients by Bootstrap and Algorithm",
    x = "Latitude",
    y = "Median Richness",
    color = "Algorithm"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10),  # Adjust facet label size
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    legend.position = "bottom"  # Place legend at the bottom
  )


library(dplyr)
library(ggplot2)

# Filter for June (Month == 6) and December (Month == 12)
df_june_december <- df_aggregated %>%
  filter(Month %in% c("6", "12")) %>%
  group_by(y, Month) %>%
  summarise(
    median_richness = median(median_richness, na.rm = TRUE),  # Compute the median across all bootstraps and algorithms
    .groups = "drop"
  )

# Plot the median richness for June and December
ggplot(df_june_december, aes(x = y, y = median_richness, color = Month, group = Month)) +
  geom_line(size = 1.2) +  # Add lines with increased thickness
  labs(
    title = "Median Latitudinal Gradients for June and December",
    x = "Latitude",
    y = "Median Richness",
    color = "Month"
  ) +
  scale_color_manual(values = c("6" = "blue", "12" = "red"), labels = c("June", "December")) +
  theme_minimal() +
  theme(
    legend.position = "top",  # Place legend at the top
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    strip.text = element_text(size = 10)  # Adjust facet label size
  )




### ==============================================================================================
### Genus Level Ensemble Members Extraction
### ==============================================================================================

### Load Required Libraries
library(abind)      # For combining multidimensional arrays
library(doParallel) # For parallel processing
library(foreach)    # For parallel loops
library(raster)     # For raster data manipulation
library(dplyr)      # For data manipulation
library(sp)         # For spatial data handling
library(tidyr)      # For data tidying and reshaping

# Load specific functions
setwd("../../CEPHALOPOD")
source(file = "../../CEPHALOPOD/code/00_config.R")

# Directories for genus-level analysis
project_wd <- "./Code/Cephalopod_output/Non_log_genusLevel/"
wd_out <- "./Code/2_Extract_model_outputs/3_Output/Genus_level/"
# Create output directory if it doesn't exist
if (!dir.exists(wd_out)) {
  dir.create(wd_out, recursive = TRUE)
}

# Get folder names to loop over
FOLDER_NAMES <- list.files(project_wd, full.names = TRUE, pattern = "chunk")

### Parallelization Setup
n.cores <- 10
cl <- makeCluster(n.cores, outfile = "")
registerDoParallel(cl)

### Process Each Folder and Subfolder
for (f in seq_along(FOLDER_NAMES)) {
# for (f in 5:length(FOLDER_NAMES)) {
# for (f in 1:4) {
  
  # Get subfolder names
  SUBFOLDER_NAMES <- list.files(FOLDER_NAMES[f], full.names = TRUE, pattern = "1000|5000")
  
  # Parallel loop over subfolders
  foreach::foreach(s = seq_along(SUBFOLDER_NAMES), .packages = c("abind", "sp", "raster", "dplyr")) %dopar% {
    
    # Print progress
    print(paste("Processing folder", f, "of", length(FOLDER_NAMES), "and subfolder", s, "of", length(SUBFOLDER_NAMES)))
    
    # Load parameter and query data
    SUBFOLDER_NAME <- SUBFOLDER_NAMES[s]
    CALL <- get(load(paste0(FOLDER_NAMES[f], "/CALL.RData")))
    QUERY <- get(load(paste0(SUBFOLDER_NAME, "/QUERY.RData")))
    
    # Check for MODEL.RData and load if exists
    model_file_path <- paste0(SUBFOLDER_NAME, "/MODEL.RData")
    if (!file.exists(model_file_path)) {
      message("MODEL.RData not found in ", SUBFOLDER_NAME, ". Skipping.")
      return(NULL)
    } else {
      MODEL <- get(load(model_file_path))
    }
    
    # Check if projections exist
    if ((length(MODEL$MODEL_LIST) == 0) & CALL$FAST == TRUE) {
      message("No validated algorithms for projections.")
      return(NULL)
    }
    
    ### Set Initial Parameters and Masks
    r0 <- CALL$ENV_DATA[[1]][[1]]
    # land <- raster::setValues(r0, ifelse(is.na(r0), 9999, NA)) # Not working anymore
    # Extract values from the raster
    r0_values <- getValues(r0)
    # Check for NA values and replace as desired
    new_values <- ifelse(is.na(r0_values), 9999, NA)
    # Set the modified values back to the raster
    land <- setValues(r0, new_values)
    
    ### Filter Recommendations
    rec <- qc_recommandations(QUERY = QUERY, MODEL = MODEL, DATA_TYPE = CALL$DATA_TYPE)
    MODEL$recommandations <- rec
    filtered_rows <- rec %>% filter(PRE_VIP == 1, FIT == 1, CUM_VIP == 1, DEV == 1)
    if (nrow(filtered_rows) == 0) {
      message("No rows meet criteria. Skipping.")
      return(NULL)
    }
    loop_over <- rownames(filtered_rows)
    
    ### Build Ensembles
    list_array <- lapply(loop_over, function(i) {
      if (CALL$DATA_TYPE == "proportions") {
        MODEL[["MBTR"]][["proj"]][["y_hat"]][,,i,]
      } else {
        MODEL[[i]][["proj"]][["y_hat"]]
      }
    })
    
   # Name list_array
   names(list_array) <- loop_over
   
   list_algorithms <- list()
    for( loop in loop_over){   

        # Extract model data
        d <- list_array[[loop]]

            # Initialize a list to store raster layers for each bootstrap-month combination
            raster_layers <- list()

            # Loop through each bootstrap
            for (bootstrap in 1:dim(d)[2]) {  # Loop over the 10 bootstraps
              # Loop through each month
              for (month in 1:dim(d)[3]) {  # Loop over the 12 months
                  # Extract the data for the current bootstrap and month
                  layer_data <- d[, bootstrap, month]
                  
                  # Create a raster layer using the template r0
                  raster_layer <- setValues(r0, layer_data)
                  
                  # Assign a meaningful name to the raster layer
                  layer_name <- paste0("Bootstrap_", bootstrap, "_Month_", month)
                  names(raster_layer) <- layer_name
                  
                  # Store the raster layer in the list
                  raster_layers[[layer_name]] <- raster_layer
              }
            }

        # Stack all in on
        r <- raster::stack(raster_layers)

        # Replace all negative values with zero in the raster stack
        r[r < 0] <- 0

        # Set the projection to WGS84
        wgs84 <- sp::CRS("+proj=longlat +datum=WGS84 +no_defs")
        raster::projection(r) <- wgs84
        
        # Save in the list
        list_algorithms[[loop]] <- r
    }

# Save list object 
save(
    list_algorithms,   
    file = paste0(wd_out, basename(SUBFOLDER_NAME), "_ensembleMembers.RData")
)
    
    # Print progress
    print(paste("Finished folder", f, "subfolder", s))
  }
}

### Stop Cluster
parallel::stopCluster(cl)

### ===============================================================================
### Create one big dataframe with all ensemble members
### ===============================================================================

# Get filenames
# fnames <- list.files("/home/deriksson/Projects/Prokaryotes_LDG_v2/Code/2_Extract_model_outputs/4_Output", full.names = TRUE) 4_Output?
fnames <- list.files("./Code/2_Extract_model_outputs/3_Output", full.names = TRUE)

# Subset filenames for Richness_5000
fnames <- fnames[grepl("Richness_5000", fnames)]
fnames <- fnames[grepl("_v2", fnames)]

# Define selected clades
selected_clades <- c(
  "d__Archaea_Richness_5000",
  "d__Bacteria_Richness_5000",
  "d__Archaea_p__Thermoproteota_c__Nitrososphaeria_Richness_5000",
  "d__Bacteria_p__Fibrobacterota_c__Fibrobacteria_Richness_5000",
  "d__Bacteria_p__Bacteroidota_c__Bacteroidia_Richness_5000",
  "d__Archaea_p__Thermoplasmatota_c__Poseidoniia_Richness_5000",
  "d__Bacteria_p__Verrucomicrobiota_c__Verrucomicrobiae_Richness_5000",
  "d__Bacteria_p__Desulfobacterota_D_c__UBA1144_Richness_5000",
  "d__Bacteria_p__Marinisomatota_c__Marinisomatia_Richness_5000",
  "d__Bacteria_p__Chloroflexota_c__Dehalococcoidia_Richness_5000",
  "d__Bacteria_p__Pseudomonadota_c__Gammaproteobacteria_Richness_5000",
  "d__Bacteria_p__SAR324_c__SAR324_Richness_5000",
  "d__Bacteria_p__Chlamydiota_Richness_5000",
  "d__Bacteria_p__Bacteroidota_c__Rhodothermia_Richness_5000",
  "d__Bacteria_p__Verrucomicrobiota_c__Kiritimatiellia_Richness_5000",
  "d__Bacteria_p__Actinomycetota_c__Acidimicrobiia_Richness_5000",
  "d__Bacteria_p__Cyanobacteriota_c__Cyanobacteriia_Richness_5000",
  "d__Bacteria_p__Pseudomonadota_c__Alphaproteobacteria_Richness_5000",
  "d__Archaea_p__Asgardarchaeota_c__Heimdallarchaeia_Richness_5000",
  "d__Bacteria_p__Bdellovibrionota_c__Bacteriovoracia_Richness_5000",
  "d__Bacteria_p__Myxococcota_c__UBA796_Richness_5000",
  "d__Bacteria_p__Myxococcota_c__UBA4151_Richness_5000",
  "d__Bacteria_p__Myxococcota_c__XYA12-FULL-58-9_Richness_5000"
)

# Filter filenames based on selected clades
fnames <- fnames[grepl(paste(selected_clades, collapse = "|"), fnames)]

# Loop through files
final_results <- list()
for( f in seq_along(fnames) ){

    list_algorithms <- get(load(fnames[f]))


            list_df <- list()
            for(i in seq_along(list_algorithms)){

            # Extract data
            d <- as.data.frame( list_algorithms[[i]], xy = TRUE) 

            # Convert to long format
            df_long <- d %>%
                tidyr::pivot_longer(
                cols = starts_with("Bootstrap"),  # Select all columns starting with "Bootstrap"
                names_to = c("Bootstrap", "Month"),  # Split column names into "Bootstrap" and "Month"
                names_pattern = "Bootstrap_(\\d+)_Month_(\\d+)",  # Regex to extract bootstrap and month
                values_to = "richness"  # Name of the column for values
                )


            # Add algorithm info
            df_long$algorithms <- names(list_algorithms)[i]

            # Save in list
            list_df[[i]] <- df_long

            }

          # Combine all data frames into one
          df_long <- bind_rows(list_df)
          head(df_long)

        # Add clade name
        clade_name <- gsub(".RData", "", basename(fnames[f]) )
        df_long$clade <- clade_name

        # Print progress
        print(paste("Finished file", f, "of", length(fnames)))

        # Save in list
        final_results[[f]] <- df_long
}

# Merge into one dataframe
df_final <- do.call("rbind", final_results)

# Check dataframe
head(df_final)

# Save complete dataframe not LDG statistics
library(data.table)
fwrite(
  df_final,
  file = paste0("/home/deriksson/Projects/Prokaryotes_LDG_v3/Code/2_Extract_model_outputs/3_Output/", "Statistics_acrossEnsembleMembers_selected_clades_v3.csv"),
  row.names = FALSE
)

# Load library
library(dplyr)

# Compute mean, sd, median, and IQR grouped by x, Month, and clade
summary_stats <- df_final %>%
  group_by(y, Month, clade) %>%
  summarise(
    mean_richness = mean(richness, na.rm = TRUE),  # Compute mean
    sd_richness = sd(richness, na.rm = TRUE),      # Compute standard deviation
    median_richness = median(richness, na.rm = TRUE),  # Compute median
    iqr_richness = IQR(richness, na.rm = TRUE),    # Compute interquartile range
    .groups = "drop"  # Ungroup after summarising
  )

# Check summary statistics
head(summary_stats)

# Save summary stats
write.csv(
  summary_stats,
  file = paste0("/home/deriksson/Projects/Prokaryotes_LDG_v3/Code/2_Extract_model_outputs/3_Output/", "LDG_Statistics_acrossEnsembleMembers_selected_clades_v3.csv"),
  row.names = FALSE

)

# Quick plot
library(ggplot2)
library(dplyr)

# Load data
summary_stats <- read.csv(paste0("./Code/2_Extract_model_outputs/3_Output/", "LDG_Statistics_acrossEnsembleMembers_selected_clades_v3.csv"))

# Filter for June and December
df_june_december <- summary_stats %>%
  filter(Month %in% c("6", "12")) %>%
  mutate(Month = factor(Month, levels = c("6", "12"), labels = c("June", "December")))

# Plot the data
ggplot(df_june_december, aes(x = y, y = mean_richness, color = Month, group = Month)) +
  geom_line(size = 1.2) +  # Line plot for mean richness
  geom_ribbon(aes(ymin = mean_richness - sd_richness, ymax = mean_richness + sd_richness, fill = Month), alpha = 0.2) +  # Ribbon for ± sd
  facet_wrap(~ clade, scales = "free_y") +  # Facet by clade
  labs(
    title = "Latitudinal Gradients for June and December by Clade",
    x = "Latitude",
    y = "Mean Richness",
    color = "Month",
    fill = "Month"
  ) +
  scale_color_manual(values = c("June" = "blue", "December" = "red")) +  # Custom colors for lines
  scale_fill_manual(values = c("June" = "blue", "December" = "red")) +  # Custom colors for ribbons
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10),  # Adjust facet label size
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    legend.position = "top"  # Place legend at the top
  )

### ===============================================================================
### Check, plottting
### ===============================================================================

list_df <- list()
for(i in seq_along(list_algorithms)){

  # Extract data
  d <- as.data.frame( list_algorithms[[i]], xy = TRUE) 

  # Convert to long format
  df_long <- d %>%
    pivot_longer(
      cols = starts_with("Bootstrap"),  # Select all columns starting with "Bootstrap"
      names_to = c("Bootstrap", "Month"),  # Split column names into "Bootstrap" and "Month"
      names_pattern = "Bootstrap_(\\d+)_Month_(\\d+)",  # Regex to extract bootstrap and month
      values_to = "richness"  # Name of the column for values
    )


  # Add algorithm info
  df_long$algorithms <- names(list_algorithms)[i]

  # Save in list
  list_df[[i]] <- df_long

}

# Combine all data frames into one
df_long <- bind_rows(list_df)
head(df_long)

# Aggregate using the median of richness
df_aggregated <- df_long %>%
  group_by(y, Bootstrap, Month, algorithms) %>%
  summarise(
    median_richness = median(richness, na.rm = TRUE),  # Calculate the median, ignoring NA values
    .groups = "drop"  # Ungroup after summarising
  )

# Check the result
head(df_aggregated)

# Check lineplot
library(ggplot2)

# Plot line plots for each bootstrap and algorithm, faceted by month
ggplot(df_aggregated, aes(x = y, y = median_richness, color = algorithms, group = interaction(Bootstrap, algorithms))) +
  geom_line(alpha = 0.7) +  # Add lines with some transparency
  facet_wrap(~ Month, scales = "free_y") +  # Facet by month
  labs(
    title = "Latitudinal Gradients by Bootstrap and Algorithm",
    x = "Latitude",
    y = "Median Richness",
    color = "Algorithm"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10),  # Adjust facet label size
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    legend.position = "bottom"  # Place legend at the bottom
  )


library(dplyr)
library(ggplot2)

# Filter for June (Month == 6) and December (Month == 12)
df_june_december <- df_aggregated %>%
  filter(Month %in% c("6", "12")) %>%
  group_by(y, Month) %>%
  summarise(
    median_richness = median(median_richness, na.rm = TRUE),  # Compute the median across all bootstraps and algorithms
    .groups = "drop"
  )

# Plot the median richness for June and December
ggplot(df_june_december, aes(x = y, y = median_richness, color = Month, group = Month)) +
  geom_line(size = 1.2) +  # Add lines with increased thickness
  labs(
    title = "Median Latitudinal Gradients for June and December",
    x = "Latitude",
    y = "Median Richness",
    color = "Month"
  ) +
  scale_color_manual(values = c("6" = "blue", "12" = "red"), labels = c("June", "December")) +
  theme_minimal() +
  theme(
    legend.position = "top",  # Place legend at the top
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    strip.text = element_text(size = 10)  # Adjust facet label size
  )


