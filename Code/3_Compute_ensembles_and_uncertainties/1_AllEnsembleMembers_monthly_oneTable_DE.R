## ====================================================================================================
## This R script compiles ensemble member outputs for species richness, Shannon, and Chao1 indices
## into long-format data frames for further analysis. It filters, reshapes, and saves monthly 
## ensemble data across all successfully modeled taxa.
##
## Author:       Dominic Eriksson
## Date:         8th of October 2025
## Affiliation:  Environmental Physics Group UP, ETH Zürich, Switzerland
## Contact:      deriksson@ethz.ch
##
## Input files:  
##   - Files in folder: Code/2_Extract_model_output/2_Output
##   - Files in folder: Code/2_Extract_model_output/4_Output
##
## Output files: 
##   - CSV files containing the long format data frames for each ensemble member.
##   - File names indicate which alpha diversity is saved in the file
##
## Strategy:
##   This script processes Cephalopod SDM ensemble outputs to compile long-format datasets 
##   for prokaryotic diversity indices (Richness, Shannon, Chao1). 
##   For each index, non-log-transformed files are preferred, with log-transformed files used as fallback. 
##   Ensemble raster stacks for each algorithm and bootstrap replicate are converted to tables 
##   with spatial coordinates, combined across algorithms and clades, reshaped to long format, 
##   annotated with bootstrap number, month, and clade, and filtered for undesired clades. 
##   The resulting tables are saved as CSVs for downstream analyses of monthly ensemble predictions 
##   and uncertainty quantification.
##
## Required R packages (tested versions):
##   - raster      (version 3.6-26)
##   - tidyr       (version 1.3.1)
##   - dplyr       (version 1.1.2)
##   - data.table  (version 1.15.4)
##
## ====================================================================================================

# Clear working space
rm(list = ls())

### ===========================
### Libraries
### ===========================
library(raster)
library(tidyr)
library(dplyr)
library(data.table)

### ===========================
### Configurable Project Paths - adjust as needed
### ===========================
wd_in  <- "Code/2_Extract_model_output/2_Output"
wd_out <- "Code/3_Compute_ensembles_and_uncertainties/1_Output"
clade_tables_dir <- file.path("Code/2_Extract_model_output/4_Output")

# Create output directory if it does not exist
if (!dir.exists(wd_out)) {
    dir.create(wd_out, recursive = TRUE)
}

### ===========================
### Function to select best file (Non-log preferred, fallback Log)
### ===========================
get_best_file <- function(non_log, log, fnames) {
  if (!is.na(non_log)) {
    pattern <- paste0(non_log, "_ensembleMembers_Non_log_v2.RData")
    file <- fnames[grepl(pattern, fnames, fixed = TRUE)]
    if (length(file) > 0) return(file[1])
  }
  if (!is.na(log)) {
    pattern <- paste0(log, "_ensembleMembers_Log_v2.RData")
    file <- fnames[grepl(pattern, fnames, fixed = TRUE)]
    if (length(file) > 0) return(file[1])
  }
  return(NA)
}

### ===========================
### Process Indices: Richness, Shannon, Chao1
### ===========================
indices <- c("Richness_5000", "Shannon_5000", "Chao1_5000")

for (index in indices) {

  # Load clade table
  clade_table <- read.csv(file.path(clade_tables_dir, paste0(index, "_taxa.csv")))
  
  # Get filenames
  fnames <- list.files(
    wd_in,
    full.names = TRUE,
    pattern = index
  )
  
  # Select best files based on non-log/log preference
  selected_files <- mapply(get_best_file, clade_table$non_log, clade_table$log, MoreArgs = list(fnames = fnames))
  selected_files <- na.omit(selected_files)
  
  fnames <- selected_files
  writeLines(fnames, con = file.path(wd_out, paste0("final_selected_files_", index, ".txt")))
  
  # Load and compile ensemble members
  l <- list()
  for (f in seq_along(fnames)) {
    message(paste0("Loading file ", f, " of ", length(fnames), ": ", basename(fnames[f])))
    
    data <- get(load(fnames[f]))
    idx <- length(data)
    
    l2 <- list()
    for (i in 1:idx) {
      df <- raster::as.data.frame(data[[i]], xy = TRUE)
      df$alg <- names(data)[i]
      df$clade <- gsub(paste0("_", index, ".*"), "", basename(fnames[f]))
      l2[[i]] <- df
    }
    
    l[[f]] <- do.call("rbind", l2)
  }
  
  df_final <- do.call("rbind", l)
  
  # Reformat table to long format
  df_long <- df_final %>%
    tidyr::pivot_longer(
      cols = starts_with("Bootstrap_"),
      names_to = "bootstrap_month",
      values_to = tolower(gsub("_5000", "", index))
    ) %>%
    tidyr::extract(
      col = bootstrap_month,
      into = c("bootstrap", "month"),
      regex = "Bootstrap_(\\d+)_Month_(\\d+)"
    ) %>%
    mutate(
      bootstrap = as.integer(bootstrap),
      month = as.integer(month)
    ) %>%
    filter(!clade %in% c("phylum_Calditrichia"))
  
  # Save long-format dataframe
  data.table::fwrite(
    df_long,
    file = file.path(wd_out, paste0("df_long_monthly_ensembleMembers_", tolower(gsub("_5000","",index)), "_5000.csv"))
  )
  
  message(paste0("Finished processing index: ", index))
}

## ====================================================================================================
## END OF SCRIPT
## ====================================================================================================
