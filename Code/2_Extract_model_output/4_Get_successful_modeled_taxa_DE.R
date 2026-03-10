## ====================================================================================================
## R script to extract successfully modeled taxa from CEPHALOPOD SDM outputs. CEPHALOPOD creates 
## output folders for all targeted taxa, but some modeling runs may fail; log files contain 
## information on why a model did not complete successfully. 
##
## Author:       Dominic Eriksson
## Date:         8th of October 2025
## Affiliation:  Environmental Physics Group, UP, ETH Zürich, Switzerland
## Contact:      deriksson@ethz.ch
##
## Input files:  Cephalopod output folders containing CALL.RData, QUERY.RData, MODEL.RData, 
##               and standard map log files (*.pdf)
##
## Output files: CSV tables listing all taxa for which successful modeling was achieved, 
##               separated by alpha diversity metric (Richness, Shannon, Chao1) and by 
##               predictor type (log or non-log)
##
## Strategy:
##   The script identifies which taxa have successfully generated standard maps by checking 
##   log and non-log output folders. Basenames of completed model folders are collected, 
##   matched across log types, and saved in structured CSV tables. This provides a concise 
##   record of taxa with completed projections for downstream analyses.
##
## ====================================================================================================


### ====================================================================
### Richness
### ====================================================================

## Check which taxa have a final map

# Directories
wd_out <- "Code/2_Extract_model_output/4_Output/"

# Input and output directories --> Save the cephalopod output data in Data_generated folder
log_dir <- "Data/Data_generated/Cephalopod_output/Log"
log_files <- list.files(log_dir, pattern = "05_standard_maps.pdf$", recursive = TRUE, full.names = TRUE)

# Non_log directory
nonlog_dir <- "Data/Data_generated/Cephalopod_output/Non_log"
nonlog_files <- list.files(nonlog_dir, pattern = "05_standard_maps.pdf$", recursive = TRUE, full.names = TRUE)

# Subset for Richness_5000
log_files_richness <- log_files[grep("Richness_5000", log_files)]
nonlog_files_richness <- nonlog_files[grep("Richness_5000", nonlog_files)]

# Create data frame with only the basename of the richness filenames
log_base <- basename(dirname(log_files_richness))
nonlog_base <- basename(dirname(nonlog_files_richness))

# Make sure both vectors are the same length for the table (fill with NA if needed)
max_len <- max(length(log_base), length(nonlog_base))
log_base <- c(log_base, rep(NA, max_len - length(log_base)))
nonlog_base <- c(nonlog_base, rep(NA, max_len - length(nonlog_base)))

# Combine into a data frame
richness_table <- data.frame(
  log = log_base,
  non_log = nonlog_base,
  stringsAsFactors = FALSE
)

print(richness_table)

# Get base names for matching
log_base <- basename(dirname(log_files_richness))
nonlog_base <- basename(dirname(nonlog_files_richness))

# Function to get full path by base name
get_path <- function(base, log_files, nonlog_files, log_base, nonlog_base) {
  idx_nonlog <- match(base, nonlog_base)
  idx_log <- match(base, log_base)
  if (!is.na(idx_nonlog)) {
    return(nonlog_files[idx_nonlog])
  } else if (!is.na(idx_log)) {
    return(log_files[idx_log])
  } else {
    return(NA)
  }
}

# All unique base names from both columns
all_bases <- unique(c(log_base, nonlog_base))

# Build the table
log_col <- sapply(all_bases, function(b) ifelse(b %in% log_base, b, NA))
nonlog_col <- sapply(all_bases, function(b) ifelse(b %in% nonlog_base, b, NA))
path_col <- sapply(all_bases, function(b) get_path(b, log_files_richness, nonlog_files_richness, log_base, nonlog_base))

richness_table <- data.frame(
  log = log_col,
  non_log = nonlog_col,
  path = path_col,
  stringsAsFactors = FALSE
)

# Print table 
head(richness_table)

# Save as CSV
write.csv(richness_table, file = paste0(wd_out, "Richness_5000_taxa.csv"), row.names = FALSE)

### ====================================================================
### Shannon
### ====================================================================

# Subset for Shannon_5000
log_files_Shannon <- log_files[grep("Shannon_5000", log_files)]
nonlog_files_Shannon <- nonlog_files[grep("Shannon_5000", nonlog_files)]

# Create data frame with only the basename of the richness filenames
log_base <- basename(dirname(log_files_Shannon))
nonlog_base <- basename(dirname(nonlog_files_Shannon))

# Make sure both vectors are the same length for the table (fill with NA if needed)
max_len <- max(length(log_base), length(nonlog_base))
log_base <- c(log_base, rep(NA, max_len - length(log_base)))
nonlog_base <- c(nonlog_base, rep(NA, max_len - length(nonlog_base)))

# Combine into a data frame
shannon_table <- data.frame(
  log = log_base,
  non_log = nonlog_base,
  stringsAsFactors = FALSE
)

# Print the table
print(shannon_table)

# Get base names for matching
log_base <- basename(dirname(log_files_Shannon))
nonlog_base <- basename(dirname(nonlog_files_Shannon))

# Function to get full path by base name
get_path <- function(base, log_files, nonlog_files, log_base, nonlog_base) {
  idx_nonlog <- match(base, nonlog_base)
  idx_log <- match(base, log_base)
  if (!is.na(idx_nonlog)) {
    return(nonlog_files[idx_nonlog])
  } else if (!is.na(idx_log)) {
    return(log_files[idx_log])
  } else {
    return(NA)
  }
}

# All unique base names from both columns
all_bases <- unique(c(log_base, nonlog_base))

# Build the table
log_col <- sapply(all_bases, function(b) ifelse(b %in% log_base, b, NA))
nonlog_col <- sapply(all_bases, function(b) ifelse(b %in% nonlog_base, b, NA))
path_col <- sapply(all_bases, function(b) get_path(b, log_files_Shannon, nonlog_files_Shannon, log_base, nonlog_base))

shannon_table <- data.frame(
  log = log_col,
  non_log = nonlog_col,
  path = path_col,
  stringsAsFactors = FALSE
)

# Print table 
head(shannon_table)

# Save as CSV
write.csv(shannon_table, file = paste0(wd_out, "Shannon_5000_taxa.csv"), row.names = FALSE)

### ====================================================================
### Chao1
### ====================================================================

# Subset for Chao1_5000
log_files_Chao1 <- log_files[grep("Chao1_5000", log_files)]
nonlog_files_Chao1 <- nonlog_files[grep("Chao1_5000", nonlog_files)]

# Create data frame with only the basename of the richness filenames
log_base <- basename(dirname(log_files_Chao1))
nonlog_base <- basename(dirname(nonlog_files_Chao1))

# Make sure both vectors are the same length for the table (fill with NA if needed)
max_len <- max(length(log_base), length(nonlog_base))
log_base <- c(log_base, rep(NA, max_len - length(log_base)))
nonlog_base <- c(nonlog_base, rep(NA, max_len - length(nonlog_base)))

# Combine into a data frame
chao1_table <- data.frame(
  log = log_base,
  non_log = nonlog_base,
  stringsAsFactors = FALSE
)

# Print the table
print(chao1_table)

# Get base names for matching
log_base <- basename(dirname(log_files_Chao1))
nonlog_base <- basename(dirname(nonlog_files_Chao1))

# Function to get full path by base name
get_path <- function(base, log_files, nonlog_files, log_base, nonlog_base) {
  idx_nonlog <- match(base, nonlog_base)
  idx_log <- match(base, log_base)
  if (!is.na(idx_nonlog)) {
    return(nonlog_files[idx_nonlog])
  } else if (!is.na(idx_log)) {
    return(log_files[idx_log])
  } else {
    return(NA)
  }
}

# All unique base names from both columns
all_bases <- unique(c(log_base, nonlog_base))

# Build the table
log_col <- sapply(all_bases, function(b) ifelse(b %in% log_base, b, NA))
nonlog_col <- sapply(all_bases, function(b) ifelse(b %in% nonlog_base, b, NA))
path_col <- sapply(all_bases, function(b) get_path(b, log_files_Chao1, nonlog_files_Chao1, log_base, nonlog_base))

chao1_table <- data.frame(
  log = log_col,
  non_log = nonlog_col,
  path = path_col,
  stringsAsFactors = FALSE
)

# Print table 
head(chao1_table)

# Save as CSV
write.csv(chao1_table, file = paste0(wd_out, "Chao1_5000_taxa.csv"), row.names = FALSE)

## ====================================================================================================
## END OF SCRIPT
## ====================================================================================================
