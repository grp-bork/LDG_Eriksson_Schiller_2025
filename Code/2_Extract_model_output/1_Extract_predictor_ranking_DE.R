## ====================================================================================================
## This R script extracts and compiles predictor importance rankings from CEPHALOPOD habitat 
## model outputs for prokaryotic taxa. It processes individual ensemble model files, aggregates predictor 
## rankings by variable and clade, and saves both detailed and combined feature importance tables for 
## downstream analyses.
##
## Author:       Dominic Eriksson
## Date:         8th of October 2025
## Affiliation:  Environmental Physics Group UP, ETH Zürich, Switzerland
## Contact:      deriksson@ethz.ch
##
## Input files:  
##   - MODEL.RData, CALL.RData, and QUERY.RData files from CEPHALOPOD output folders for log and 
##     non-log environmental predictors.
##
## Output files: 
##   - CSV files containing predictor rankings for each ensemble member.
##   - Aggregated feature importance tables summarizing median, interquartile range, and bounds of 
##     predictor contributions across taxa and predictor transformations.
##
## Strategy:
##   The script automates the extraction and compilation of variable importance metrics from ensemble SDM 
##   outputs. For each clade, model outputs are filtered based on quality criteria, and predictor rankings 
##   from log-transformed and non-log-transformed environmental variables are merged, ensuring complete 
##   coverage across taxa. The resulting combined tables provide a robust overview of the environmental 
##   drivers shaping prokaryotic diversity and serve as input for downstream analyses of model interpretation 
##   and cross-taxa comparisons.
##
## Required R packages (tested versions):
##   - parallel 4.2.2
##   - raster   3.6.26
##   - dplyr    1.1.2
##   - tibble   3.3.0
## ====================================================================================================


# Clear workspace
rm(list = ls())

# Load necessary libraries
library(parallel)
library(raster)
library(dplyr)
library(tibble)

# Load configuration file from the Cephalopod modeling pipeline (https://github.com/alexschickele/CEPHALOPOD)
source("./CEPHALOPOD/code/00_config.R")

# SWITCH - Choose log or non-log predictors
trans <- c("Non_log", "Log") # Just in case you ran the model two times with and wiithout log transfomred predictors.
t <- 1

# Depth levels
depth_level <- c("MLD", "0_10Meters") # we only use MLD in our analysis in our manuscript --> d = 1
d <- 1

# Input and output directories --> Save the cephalopod output data in Data_generated folderdirectory <- paste0("Data/Data_generated/Cephalopod_output/", trans[t])
directory <- paste0("Data/Data_generated/Cephalopod_output/", trans[t])
wd_out <- paste0("Code/2_Extract_model_output/1_Output/", trans[t], "/")
if (!dir.exists(wd_out)) dir.create(wd_out, recursive = TRUE)

# List model files
model_files <- list.files(directory, full.names = TRUE, recursive = TRUE, pattern = "MODEL.RData")

# Function to process each model file
process_model_file <- function(model_file) {
  CALL <- get(load(file.path(dirname(dirname(model_file)), "CALL.RData")))
  model <- get(load(model_file))
  QUERY <- get(load(gsub("MODEL.RData", "QUERY.RData", model_file)))

  if ((length(model$MODEL_LIST) == 0) & CALL$FAST == TRUE) return(NULL)

  if (!"ENSEMBLE" %in% names(model)) {
    rec <- qc_recommandations(QUERY = QUERY, MODEL = MODEL, DATA_TYPE = CALL$DATA_TYPE)
    MODEL$recommandations <- rec
    filtered_rows <- rec %>% filter(PRE_VIP == 1, FIT == 1, CUM_VIP == 1, DEV == 1)
    if (nrow(filtered_rows) == 0) return(NULL)
    loop_over <- rownames(filtered_rows)
    df_predRanking <- model[[loop_over]]$vip
  } else {
    df_predRanking <- model$ENSEMBLE$vip
  }

  df_predRanking$clade <- basename(dirname(model_file))
  return(df_predRanking)
}

# Parallel processing
df_list <- parallel::mclapply(model_files, process_model_file, mc.cores = detectCores()/2)
df_list <- Filter(Negate(is.null), df_list)
df_list <- df_list[sapply(df_list, tibble::is_tibble)]
df <- do.call(rbind, df_list)

# Save merged predictor rankings
write.csv(df, file = file.path(wd_out, paste0("Predictor_rankings_", depth_level[d], ".csv")), row.names = FALSE)

# Aggregate predictor rankings by variable and clade
df_aggregated <- df %>%
  dplyr::group_by(variable, clade) %>%
  dplyr::summarize(var_explained = median(value, na.rm = TRUE),
            iqr = IQR(value, na.rm = TRUE),
            lower_bound = var_explained - iqr / 2,
            upper_bound = var_explained + iqr / 2) %>%
  dplyr::rename(species = clade)

# Save aggregated feature importance
write.csv(df_aggregated, file = file.path(wd_out, paste0("feature_importance_combined_", depth_level[d], ".csv")), row.names = FALSE)

## Combine non-log and log feature importance (portable paths) - Note that you need to have run t = 1 and t = 2 above first

# Input directories for log and non-log
wd_in_log <- "/net/sea/work/deriksson/Projects/Prokaryotes_LDG_v3_github/Code/1_Extract_model_output/1_Output/Log/"
wd_in_nonlog <- "/net/sea/work/deriksson/Projects/Prokaryotes_LDG_v3_github/Code/1_Extract_model_output/1_Output/Non_log/"

# Load files
df_nonlog <- read.csv(paste0(wd_in_nonlog, "feature_importance_combined_MLD.csv"))
df_log <- read.csv(paste0(wd_in_log, "feature_importance_combined_MLD.csv"))

# Let's add clades that have only been modeled via log transformed env. parameters
missing_entries <- df_log %>% anti_join(df_nonlog, by = "species")
df_combined <- bind_rows(df_nonlog, missing_entries)

# Save - We save outside of the log and non-log folders for clarity. 
write.csv(df_combined, file = paste0(wd_out, "feature_importance_combined_V4_MLD.csv"), row.names = FALSE)


## ====================================================================================================
## END OF SCRIPT
## ====================================================================================================
