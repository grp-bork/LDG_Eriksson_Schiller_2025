## ====================================================================================================
## This R script creates comprehensive multi-plot visualizations of Latitudinal Diversity Gradients
## (LDG) for prokaryotic genera using aggregated ensemble data. The script processes final LDG datasets
## to generate publication-ready multi-panel figures showing richness and Shannon diversity patterns
## across latitude bands, with uncertainty visualization through confidence ribbons and systematic
## organization into paginated layouts for optimal presentation of large numbers of genera.
##
## Author:       Dominic Eriksson
## Date:         6th of March 2026
## Affiliation:  Environmental Physics Group UP, ETH Zürich, Switzerland
## Contact:      deriksson@ethz.ch
##
## Input files:
##   - Richness LDG data: ./Code/3_Compute_ensembles_and_uncertainties/2e_Output/LDG_final_aggregated_*.csv
##   - Shannon LDG data: ./Code/3_Compute_ensembles_and_uncertainties/2e_Output/LDG_shannon_final_aggregated_*.csv
##
## Output files:
##   - Richness LDG pages: ./Code/3_Compute_ensembles_and_uncertainties/2f_Output/LDG_page_*.png
##   - Shannon LDG pages: ./Code/3_Compute_ensembles_and_uncertainties/2f_Output/LDG_shannon_page_*.png
##   - Genus lookup tables: ./Code/3_Compute_ensembles_and_uncertainties/2f_Output/genus_lookup_table*.csv
##
## Strategy:
##   This script implements a comprehensive LDG visualization pipeline for large-scale prokaryotic
##   diversity analysis. The workflow includes: (1) loading aggregated LDG datasets with ensemble
##   statistics and filtering out genera without valid diversity data to focus on meaningful patterns;
##   (2) creating individual LDG plots for each genus with uncertainty ribbons showing standard
##   deviations, using distinct color schemes for richness (blue) and Shannon diversity (purple) to
##   enable clear differentiation between diversity metrics; (3) organizing plots into systematic
##   multi-panel layouts with 12 plots per page arranged in 3-column grids for optimal readability
##   and comparison; (4) generating comprehensive lookup tables mapping genera to specific pages and
##   files for efficient navigation; and (5) implementing quality control procedures to handle missing
##   data gracefully and ensure all valid genera are processed with consistent formatting standards.
##
## Required R packages (tested versions):
##   - ggplot2       (version 3.4.2)
##   - data.table    (version 1.15.4)
##   - patchwork     (version 1.1.3)
##
## ====================================================================================================

# Clear working space
rm(list = ls())

### ===========================
### Libraries
### ===========================
library(ggplot2)     # For creating publication-quality plots and visualizations
library(data.table)  # For fast data loading and manipulation
library(patchwork)   # For combining multiple plots into publication layouts

### ===========================
### Configurable Project Paths - adjust as needed
### ===========================
wd_in <- "./Code/3_Compute_ensembles_and_uncertainties/2e_Output/"
wd_out <- "./Code/3_Compute_ensembles_and_uncertainties/2f_Output/"

# Create output directory if it doesn't exist
if(!dir.exists(wd_out)) {
    dir.create(wd_out, recursive = TRUE)
    cat("Created output directory:", wd_out, "\n")
}

## ====================================================================================================
## SECTION 1: Richness Latitudinal Diversity Gradient Visualizations
## ====================================================================================================

# Get all aggregated LDG files
fnames <- list.files(wd_in, pattern = "LDG_final_aggregated.*\\.csv$", full.names = TRUE)

cat("Found", length(fnames), "aggregated LDG files to process\n")

if(length(fnames) == 0) {
    stop("No aggregated LDG files found in input directory!")
}

# Function to extract genus name from filename
extract_genus_name <- function(filename) {
    genus_name <- gsub("LDG_final_aggregated_", "", basename(filename))
    genus_name <- gsub("\\.csv$", "", genus_name)
    genus_name <- gsub("genus_", "", genus_name)
    return(genus_name)
}

# Function to create a single LDG plot
create_ldg_plot <- function(data, genus_name) {
    
    # Filter out NA values for plotting (both mean and sd)
    plot_data <- data[!is.na(data$mean_annual) & !is.na(data$sd_annual), ]
    
    # Skip if no valid data
    if(nrow(plot_data) == 0) {
        return(NULL)
    }
    
    # Sort by latitude to ensure proper line connection
    plot_data <- plot_data[order(plot_data$y), ]
    
    # Create the LDG plot with ribbon for standard deviation
    p <- ggplot(plot_data, aes(x = mean_annual, y = y)) +
        geom_ribbon(aes(xmin = pmax(0, mean_annual - sd_annual), 
                       xmax = mean_annual + sd_annual), 
                   alpha = 0.3, fill = "lightblue", na.rm = TRUE) +
        geom_point(color = "darkblue", size = 0.5) +
        geom_path(color = "darkblue", size = 0.4) +
        labs(
            title = genus_name,
            x = "Richness",
            y = "Latitude"
        ) +
        theme_minimal() +
        theme(
            plot.title = element_text(size = 9, hjust = 0.5, face = "bold"),
            axis.text = element_text(size = 7),
            axis.title = element_text(size = 8),
            panel.grid.minor = element_blank()
        )
    
    return(p)
}

# Settings
plots_per_page <- 12
n_files <- length(fnames)

# First, filter out files with no valid data
valid_files <- c()
valid_genera <- c()

cat("Checking files for valid data...\n")
for(i in seq_along(fnames)) {
    df <- fread(fnames[i])
    if(any(!is.na(df$mean_annual))) {
        valid_files <- c(valid_files, fnames[i])
        valid_genera <- c(valid_genera, extract_genus_name(fnames[i]))
        cat("Valid:", extract_genus_name(fnames[i]), "\n")
    } else {
        cat("Skipping (no data):", extract_genus_name(fnames[i]), "\n")
    }
}

cat("Found", length(valid_files), "genera with valid data\n")

n_valid_files <- length(valid_files)
n_pages <- ceiling(n_valid_files / plots_per_page)

cat("Creating", n_pages, "pages with up to", plots_per_page, "plots each\n")

# Create summary table for valid genera only
summary_table <- data.frame(
    genus = valid_genera,
    page_number = ceiling((1:n_valid_files) / plots_per_page),
    multiplot_file = paste0("LDG_page_", sprintf("%02d", ceiling((1:n_valid_files) / plots_per_page)), ".png"),
    stringsAsFactors = FALSE
)

# Process files in groups
for(page in 1:n_pages) {
    
    cat("Creating page", page, "of", n_pages, "\n")
    
    # Get file indices for this page
    start_idx <- (page - 1) * plots_per_page + 1
    end_idx <- min(page * plots_per_page, n_valid_files)
    current_files <- valid_files[start_idx:end_idx]
    
    # Create plots for this page
    plots <- list()
    
    for(i in seq_along(current_files)) {
        file_path <- current_files[i]
        genus_name <- extract_genus_name(file_path)
        
        cat("  Loading:", genus_name, "\n")
        
        # Load data and create plot
        df <- fread(file_path)
        plot_obj <- create_ldg_plot(df, genus_name)
        
        if(!is.null(plot_obj)) {
            plots[[length(plots) + 1]] <- plot_obj
        }
    }
    
    # Only create multiplot if we have plots
    if(length(plots) > 0) {
        
        # Add empty plots if needed to fill grid (optional)
        while(length(plots) < plots_per_page && length(plots) < 12) {
            empty_plot <- ggplot() + theme_void()
            plots[[length(plots) + 1]] <- empty_plot
        }
        
        # Combine plots using patchwork
        combined_plot <- wrap_plots(plots, ncol = 3)
        
        # Add title
        final_plot <- combined_plot + plot_annotation(
            title = paste("Prokaryotic Genus Latitudinal Diversity Gradients - Page", page),
            theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))
        )
        
        # Save the plot
        output_filename <- paste0("LDG_page_", sprintf("%02d", page), ".png")
        output_path <- file.path(wd_out, output_filename)
        
        ggsave(output_path, final_plot, width = 16, height = 20, dpi = 300)
        
        cat("  Saved:", output_filename, "\n\n")
    }
}

# Save summary table
write.csv(summary_table, file.path(wd_out, "genus_lookup_table.csv"), row.names = FALSE)

# Print summary
cat("COMPLETE!\n")
cat("Total files found:", n_files, "\n")
cat("Genera with valid data:", n_valid_files, "\n")
cat("Pages created:", n_pages, "\n")
cat("Summary saved as: genus_lookup_table.csv\n")

# Show first 10 entries
cat("\nFirst 10 entries:\n")
print(head(summary_table, 10))

cat("\nSUCCESS: All valid genera processed!\n")

## ====================================================================================================
## SECTION 2: Shannon Diversity Latitudinal Diversity Gradient Visualizations
## ====================================================================================================

# Get all aggregated Shannon LDG files
fnames_shannon <- list.files(wd_in, pattern = "LDG_shannon_final_aggregated.*\\.csv$", full.names = TRUE)

cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("SHANNON DIVERSITY LDG VISUALIZATION\n")
cat(paste(rep("=", 80), collapse=""), "\n")
cat("Found", length(fnames_shannon), "aggregated Shannon LDG files to process\n")

if(length(fnames_shannon) > 0) {
    
    # Function to extract genus name from Shannon filename
    extract_genus_name_shannon <- function(filename) {
        genus_name <- gsub("LDG_shannon_final_aggregated_", "", basename(filename))
        genus_name <- gsub("\\.csv$", "", genus_name)
        genus_name <- gsub("genus_", "", genus_name)
        return(genus_name)
    }
    
    # Function to create a single Shannon LDG plot
    create_ldg_plot_shannon <- function(data, genus_name) {
        
        # Filter out NA values for plotting (both mean and sd)
        plot_data <- data[!is.na(data$mean_annual) & !is.na(data$sd_annual), ]
        
        # Skip if no valid data
        if(nrow(plot_data) == 0) {
            return(NULL)
        }
        
        # Sort by latitude to ensure proper line connection
        plot_data <- plot_data[order(plot_data$y), ]
        
        # Create the Shannon LDG plot with ribbon for standard deviation
        p <- ggplot(plot_data, aes(x = mean_annual, y = y)) +
            geom_ribbon(aes(xmin = pmax(0, mean_annual - sd_annual), 
                           xmax = mean_annual + sd_annual), 
                       alpha = 0.3, fill = "plum", na.rm = TRUE) +
            geom_point(color = "darkorchid", size = 0.5) +
            geom_path(color = "darkorchid", size = 0.4) +
            labs(
                title = genus_name,
                x = "Shannon Diversity",
                y = "Latitude"
            ) +
            theme_minimal() +
            theme(
                plot.title = element_text(size = 9, hjust = 0.5, face = "bold"),
                axis.text = element_text(size = 7),
                axis.title = element_text(size = 8),
                panel.grid.minor = element_blank()
            )
        
        return(p)
    }
    
    # Settings for Shannon
    plots_per_page_shannon <- 12
    n_files_shannon <- length(fnames_shannon)
    
    # First, filter out Shannon files with no valid data
    valid_files_shannon <- c()
    valid_genera_shannon <- c()
    
    cat("Checking Shannon files for valid data...\n")
    for(i in seq_along(fnames_shannon)) {
        df <- fread(fnames_shannon[i])
        if(any(!is.na(df$mean_annual))) {
            valid_files_shannon <- c(valid_files_shannon, fnames_shannon[i])
            valid_genera_shannon <- c(valid_genera_shannon, extract_genus_name_shannon(fnames_shannon[i]))
            cat("Valid Shannon:", extract_genus_name_shannon(fnames_shannon[i]), "\n")
        } else {
            cat("Skipping Shannon (no data):", extract_genus_name_shannon(fnames_shannon[i]), "\n")
        }
    }
    
    cat("Found", length(valid_files_shannon), "Shannon genera with valid data\n")
    
    n_valid_files_shannon <- length(valid_files_shannon)
    n_pages_shannon <- ceiling(n_valid_files_shannon / plots_per_page_shannon)
    
    cat("Creating", n_pages_shannon, "Shannon pages with up to", plots_per_page_shannon, "plots each\n")
    
    # Create summary table for valid Shannon genera only
    summary_table_shannon <- data.frame(
        genus = valid_genera_shannon,
        page_number = ceiling((1:n_valid_files_shannon) / plots_per_page_shannon),
        multiplot_file = paste0("LDG_shannon_page_", sprintf("%02d", ceiling((1:n_valid_files_shannon) / plots_per_page_shannon)), ".png"),
        stringsAsFactors = FALSE
    )
    
    # Process Shannon files in groups
    for(page in 1:n_pages_shannon) {
        
        cat("Creating Shannon page", page, "of", n_pages_shannon, "\n")
        
        # Get file indices for this page
        start_idx <- (page - 1) * plots_per_page_shannon + 1
        end_idx <- min(page * plots_per_page_shannon, n_valid_files_shannon)
        current_files <- valid_files_shannon[start_idx:end_idx]
        
        # Create plots for this page
        plots <- list()
        
        for(i in seq_along(current_files)) {
            file_path <- current_files[i]
            genus_name <- extract_genus_name_shannon(file_path)
            
            cat("  Loading Shannon:", genus_name, "\n")
            
            # Load data and create Shannon plot
            df <- fread(file_path)
            plot_obj <- create_ldg_plot_shannon(df, genus_name)
            
            if(!is.null(plot_obj)) {
                plots[[length(plots) + 1]] <- plot_obj
            }
        }
        
        # Only create multiplot if we have plots
        if(length(plots) > 0) {
            
            # Add empty plots if needed to fill grid (optional)
            while(length(plots) < plots_per_page_shannon && length(plots) < 12) {
                empty_plot <- ggplot() + theme_void()
                plots[[length(plots) + 1]] <- empty_plot
            }
            
            # Combine plots using patchwork
            combined_plot <- wrap_plots(plots, ncol = 3)
            
            # Add title for Shannon
            final_plot <- combined_plot + plot_annotation(
                title = paste("Prokaryotic Genus Shannon Diversity Latitudinal Gradients - Page", page),
                theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))
            )
            
            # Save the Shannon plot
            output_filename <- paste0("LDG_shannon_page_", sprintf("%02d", page), ".png")
            output_path <- file.path(wd_out, output_filename)
            
            ggsave(output_path, final_plot, width = 16, height = 20, dpi = 300)
            
            cat("  Saved:", output_filename, "\n\n")
        }
    }
    
    # Save Shannon summary table
    write.csv(summary_table_shannon, file.path(wd_out, "genus_lookup_table_shannon.csv"), row.names = FALSE)
    
    # Print Shannon summary
    cat("SHANNON COMPLETE!\n")
    cat("Total Shannon files found:", n_files_shannon, "\n")
    cat("Shannon genera with valid data:", n_valid_files_shannon, "\n")
    cat("Shannon pages created:", n_pages_shannon, "\n")
    cat("Shannon summary saved as: genus_lookup_table_shannon.csv\n")
    
    # Show first 10 Shannon entries
    cat("\nFirst 10 Shannon entries:\n")
    print(head(summary_table_shannon, 10))
    
    cat("\nSUCCESS: All valid Shannon genera processed!\n")
    
} else {
    cat("No Shannon diversity files found in:", wd_in, "\n")
    cat("Looking for files matching pattern: LDG_shannon_final_aggregated_*.csv\n")
}

## ====================================================================================================
## END OF SCRIPT
## ====================================================================================================