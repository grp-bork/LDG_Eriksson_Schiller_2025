## ====================================================================================================
## This R script creates comprehensive multi-plot maps for prokaryotic genera showing annual mean
## richness and Shannon diversity patterns. It processes filtered ensemble datasets to generate
## publication-ready visualizations with consistent color schemes, spatial layouts, and comprehensive
## documentation. The script handles both richness and Shannon diversity data with different color
## palettes for clear distinction between diversity metrics.
##
## Author:       Dominic Eriksson
## Date:         6th of March 2026
## Affiliation:  Environmental Physics Group UP, ETH Zürich, Switzerland
## Contact:      deriksson@ethz.ch
##
## Input files:
##   - Filtered richness data: ./Code/3_Compute_ensembles_and_uncertainties/2c_Output/marginal_seas_filtered_*.csv
##   - Filtered Shannon data: ./Code/3_Compute_ensembles_and_uncertainties/2c_Output/marginal_seas_filtered_shannon_*.csv
##
## Output files:
##   - Richness multiplots: ./Code/3_Compute_ensembles_and_uncertainties/2d_Output/genus_maps_multiplot_*.png
##   - Shannon multiplots: ./Code/3_Compute_ensembles_and_uncertainties/2d_Output/genus_maps_shannon_multiplot_*.png
##   - Color scale references: ./Code/3_Compute_ensembles_and_uncertainties/2d_Output/*_colorscale_reference*.png
##   - Processing summaries: ./Code/3_Compute_ensembles_and_uncertainties/2d_Output/genus_mapping_summary*.csv
##
## Strategy:
##   This script creates publication-quality multi-plot visualizations of prokaryotic diversity patterns.
##   The workflow includes: (1) loading filtered ensemble datasets and organizing them into systematic
##   groups for multi-plot layouts (6 plots per figure for optimal readability); (2) generating richness
##   maps using viridis color scheme with square-root transformation for enhanced visualization of patterns;
##   (3) creating Shannon diversity maps using plasma color scheme to distinguish from richness data;
##   (4) producing comprehensive documentation including color scale references, processing summaries, and
##   genus mapping tables; and (5) implementing quality control checks to ensure all expected genera are
##   processed correctly with consistent formatting and naming conventions.
##
## Required R packages (tested versions):
##   - ggplot2       (version 3.4.2)
##   - data.table    (version 1.15.4)
##   - gridExtra     (version 2.3)
##   - grid          (version 4.3.2)
##   - scales        (version 1.2.1)
##   - viridis       (version 0.6.3)
##
## ====================================================================================================

# Clear working space
rm(list = ls())

### ===========================
### Libraries
### ===========================
library(ggplot2)     # For creating publication-quality plots and maps
library(data.table)  # For fast data manipulation and file I/O
library(gridExtra)   # For arranging multiple plots in grid layouts
library(grid)        # For low-level grid graphics and layout control
library(scales)      # For axis scaling and transformation functions
library(viridis)     # For colorblind-friendly color palettes

### ===========================
### Configurable Project Paths - adjust as needed
### ===========================
# Input and output directories
wd_in <- "./Code/3_Compute_ensembles_and_uncertainties/2c_Output/"
wd_out <- "./Code/3_Compute_ensembles_and_uncertainties/2d_Output/"

# Create output directory if it doesn't exist
if(!dir.exists(wd_out)) {
    dir.create(wd_out, recursive = TRUE)
    cat("Created output directory:", wd_out, "\n")
}

## ====================================================================================================
## SECTION 1: Species Richness Multi-Plot Generation
## ====================================================================================================

# Get all filtered CSV files
fnames <- list.files(wd_in, pattern = "marginal_seas_filtered_.*\\.csv$", full.names = TRUE)

cat("Found", length(fnames), "filtered CSV files to process\n")

if(length(fnames) == 0) {
    stop("No filtered CSV files found in input directory!")
}

# Function to create a single genus map
create_genus_map <- function(data, genus_name) {
    
    # Create the map with viridis color scheme
    p <- ggplot(data, aes(x = x, y = y, fill = mean_annual)) +
        geom_raster() +
        scale_fill_viridis_c(
            name = "Annual Mean\nRichness",
            na.value = "grey80",
            option = "viridis",  # Options: "viridis", "plasma", "inferno", "magma", "cividis"
            trans = "sqrt",      # Square root transformation for better visualization
            labels = function(x) sprintf("%.2f", x)
        ) +
        labs(
            title = "Annual Mean Richness",
            subtitle = genus_name,
            x = "Longitude",
            y = "Latitude"
        ) +
        theme_minimal() +
        theme(
            plot.title = element_text(size = 12, hjust = 0.5),
            plot.subtitle = element_text(size = 14, hjust = 0.5, face = "bold"),
            axis.text = element_text(size = 8),
            axis.title = element_text(size = 10),
            legend.title = element_text(size = 9),
            legend.text = element_text(size = 8),
            legend.key.size = unit(0.8, "cm"),
            plot.margin = margin(5, 5, 5, 5)
        )
    
    return(p)
}

# Function to create multi-plot with multiple genera
create_multiplot <- function(file_list, output_filename, plots_per_page = 6) {
    
    plot_list <- list()
    genus_names <- character()
    
    # Create individual plots
    for(i in seq_along(file_list)) {
        cat("  Loading", basename(file_list[i]), "\n")
        
        # Load data
        data <- fread(file_list[i])
        genus_name <- unique(data$clade)[1]
        genus_names[i] <- genus_name
        
        # Clean genus name for display
        genus_display <- gsub("genus_", "", genus_name)
        genus_display <- gsub("_", " ", genus_display)
        
        # Create map
        plot_list[[i]] <- create_genus_map(data, genus_display)
    }
    
    # Arrange plots in grid
    if(length(plot_list) <= 4) {
        # 2x2 grid for <= 4 plots
        ncol <- 2
        nrow <- 2
    } else if(length(plot_list) <= 6) {
        # 2x3 grid for 5-6 plots
        ncol <- 3
        nrow <- 2
    } else {
        # 3x3 grid for 7-9 plots (max)
        ncol <- 3
        nrow <- 3
    }
    
    # Create the multiplot
    multiplot <- do.call(grid.arrange, c(plot_list, ncol = ncol, nrow = nrow))
    
    # Save the multiplot
    output_path <- file.path(wd_out, paste0(output_filename, ".png"))
    ggsave(output_path, multiplot, width = 16, height = 12, dpi = 300)
    
    return(genus_names)
}

# Determine how to group files for multiplots
plots_per_multiplot <- 6  # Adjust this number based on desired layout
n_files <- length(fnames)
n_multiplots <- ceiling(n_files / plots_per_multiplot)

cat("\nCreating", n_multiplots, "multiplots with up to", plots_per_multiplot, "plots each\n")
cat("Using viridis color scheme for consistent, colorblind-friendly visualization\n")
cat(paste(rep("=", 80), collapse=""), "\n")

# Process files in groups
all_genera <- character()

for(i in 1:n_multiplots) {
    
    # Calculate file indices for this multiplot
    start_idx <- (i - 1) * plots_per_multiplot + 1
    end_idx <- min(i * plots_per_multiplot, n_files)
    
    current_files <- fnames[start_idx:end_idx]
    
    cat("Creating multiplot", i, "of", n_multiplots, "\n")
    cat("Processing files", start_idx, "to", end_idx, "\n")
    
    # Create output filename
    multiplot_name <- sprintf("genus_maps_multiplot_%02d", i)
    
    # Create multiplot
    genera_in_plot <- create_multiplot(current_files, multiplot_name, plots_per_multiplot)
    all_genera <- c(all_genera, genera_in_plot)
    
    cat("Saved:", paste0(multiplot_name, ".png"), "\n")
    cat("Genera included:", paste(genera_in_plot, collapse = ", "), "\n\n")
}

# Create a sample plot showing the viridis color scale
cat("Creating color scale reference...\n")

# Load first file as example
sample_data <- fread(fnames[1])
sample_genus <- gsub("genus_", "", unique(sample_data$clade)[1])
sample_genus <- gsub("_", " ", sample_genus)

p_colorscale <- ggplot(sample_data, aes(x = x, y = y, fill = mean_annual)) +
    geom_raster() +
    scale_fill_viridis_c(
        name = "Annual Mean Richness",
        na.value = "grey80",
        option = "viridis",
        trans = "sqrt"
    ) +
    labs(
        title = "Viridis Color Scale Reference",
        subtitle = paste("Example from", sample_genus),
        x = "Longitude",
        y = "Latitude"
    ) +
    theme_minimal() +
    theme(
        plot.title = element_text(size = 16, hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)
    )

ggsave(file.path(wd_out, "viridis_colorscale_reference.png"), 
       p_colorscale, width = 12, height = 8, dpi = 300)

# Create summary document
summary_df <- data.frame(
    File_Number = 1:length(fnames),
    Genus = all_genera,
    Input_File = basename(fnames),
    Multiplot_Number = ceiling((1:length(fnames)) / plots_per_multiplot),
    stringsAsFactors = FALSE
)

# Save summary
write.csv(summary_df, file.path(wd_out, "genus_mapping_summary.csv"), row.names = FALSE)

# Print final summary
cat(paste(rep("=", 80), collapse=""), "\n")
cat("MAPPING COMPLETE\n")
cat(paste(rep("=", 80), collapse=""), "\n")
cat("Total genera processed:", length(all_genera), "\n")
cat("Total multiplots created:", n_multiplots, "\n")
cat("Plots per multiplot:", plots_per_multiplot, "\n")
cat("Color scheme: Viridis (colorblind-friendly)\n")
cat("Output directory:", wd_out, "\n")

cat("\nFiles created:\n")
multiplot_files <- list.files(wd_out, pattern = "genus_maps_multiplot_.*\\.png$")
for(file in multiplot_files) {
    cat("- ", file, "\n")
}
cat("- viridis_colorscale_reference.png\n")
cat("- genus_mapping_summary.csv\n")

cat("\nFirst 10 genera processed:\n")
print(head(summary_df[, c("Genus", "Multiplot_Number")], 10))

if(length(all_genera) != 272) {
    cat("\nWARNING: Expected 272 genera but processed", length(all_genera), "\n")
} else {
    cat("\nSUCCESS: All 272 genera mapped with viridis color scheme!\n")
}

## ====================================================================================================
## SECTION 2: Shannon Diversity Multi-Plot Generation
## ====================================================================================================

# Get all filtered Shannon diversity CSV files
fnames_shannon <- list.files(wd_in, pattern = "marginal_seas_filtered_shannon_.*\\.csv$", full.names = TRUE)

cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("SHANNON DIVERSITY MAPPING\n")
cat(paste(rep("=", 80), collapse=""), "\n")
cat("Found", length(fnames_shannon), "filtered Shannon CSV files to process\n")

if(length(fnames_shannon) > 0) {
    
    # Function to create a single genus map for Shannon diversity
    create_genus_map_shannon <- function(data, genus_name) {
        
        # Create the map with plasma color scheme for Shannon diversity
        p <- ggplot(data, aes(x = x, y = y, fill = mean_annual)) +
            geom_raster() +
            scale_fill_viridis_c(
                name = "Annual Mean\nShannon Diversity",
                na.value = "grey80",
                option = "plasma",  # Different color scheme for Shannon
                trans = "sqrt",     # Square root transformation for better visualization
                labels = function(x) sprintf("%.2f", x)
            ) +
            labs(
                title = "Annual Mean Shannon Diversity",
                subtitle = genus_name,
                x = "Longitude",
                y = "Latitude"
            ) +
            theme_minimal() +
            theme(
                plot.title = element_text(size = 12, hjust = 0.5),
                plot.subtitle = element_text(size = 14, hjust = 0.5, face = "bold"),
                axis.text = element_text(size = 8),
                axis.title = element_text(size = 10),
                legend.title = element_text(size = 9),
                legend.text = element_text(size = 8),
                legend.key.size = unit(0.8, "cm"),
                plot.margin = margin(5, 5, 5, 5)
            )
        
        return(p)
    }
    
    # Function to create multi-plot with multiple genera for Shannon diversity
    create_multiplot_shannon <- function(file_list, output_filename, plots_per_page = 6) {
        
        plot_list <- list()
        genus_names <- character()
        
        # Create individual plots
        for(i in seq_along(file_list)) {
            cat("  Loading", basename(file_list[i]), "\n")
            
            # Load data
            data <- fread(file_list[i])
            genus_name <- unique(data$clade)[1]
            genus_names[i] <- genus_name
            
            # Clean genus name for display
            genus_display <- gsub("genus_", "", genus_name)
            genus_display <- gsub("_", " ", genus_display)
            
            # Create Shannon diversity map
            plot_list[[i]] <- create_genus_map_shannon(data, genus_display)
        }
        
        # Arrange plots in grid
        if(length(plot_list) <= 4) {
            # 2x2 grid for <= 4 plots
            ncol <- 2
            nrow <- 2
        } else if(length(plot_list) <= 6) {
            # 2x3 grid for 5-6 plots
            ncol <- 3
            nrow <- 2
        } else {
            # 3x3 grid for 7-9 plots (max)
            ncol <- 3
            nrow <- 3
        }
        
        # Create the multiplot
        multiplot <- do.call(grid.arrange, c(plot_list, ncol = ncol, nrow = nrow))
        
        # Save the multiplot
        output_path <- file.path(wd_out, paste0(output_filename, ".png"))
        ggsave(output_path, multiplot, width = 16, height = 12, dpi = 300)
        
        return(genus_names)
    }
    
    # Determine how to group Shannon files for multiplots
    plots_per_multiplot_shannon <- 6  # Adjust this number based on desired layout
    n_files_shannon <- length(fnames_shannon)
    n_multiplots_shannon <- ceiling(n_files_shannon / plots_per_multiplot_shannon)
    
    cat("\nCreating", n_multiplots_shannon, "Shannon multiplots with up to", plots_per_multiplot_shannon, "plots each\n")
    cat("Using plasma color scheme for Shannon diversity visualization\n")
    cat(paste(rep("=", 80), collapse=""), "\n")
    
    # Process Shannon files in groups
    all_genera_shannon <- character()
    
    for(i in 1:n_multiplots_shannon) {
        
        # Calculate file indices for this multiplot
        start_idx <- (i - 1) * plots_per_multiplot_shannon + 1
        end_idx <- min(i * plots_per_multiplot_shannon, n_files_shannon)
        
        current_files <- fnames_shannon[start_idx:end_idx]
        
        cat("Creating Shannon multiplot", i, "of", n_multiplots_shannon, "\n")
        cat("Processing files", start_idx, "to", end_idx, "\n")
        
        # Create output filename for Shannon
        multiplot_name <- sprintf("genus_maps_shannon_multiplot_%02d", i)
        
        # Create Shannon multiplot
        genera_in_plot <- create_multiplot_shannon(current_files, multiplot_name, plots_per_multiplot_shannon)
        all_genera_shannon <- c(all_genera_shannon, genera_in_plot)
        
        cat("Saved:", paste0(multiplot_name, ".png"), "\n")
        cat("Genera included:", paste(genera_in_plot, collapse = ", "), "\n\n")
    }
    
    # Create a sample plot showing the plasma color scale for Shannon
    cat("Creating Shannon color scale reference...\n")
    
    # Load first Shannon file as example
    sample_data_shannon <- fread(fnames_shannon[1])
    sample_genus_shannon <- gsub("genus_", "", unique(sample_data_shannon$clade)[1])
    sample_genus_shannon <- gsub("_", " ", sample_genus_shannon)
    
    p_colorscale_shannon <- ggplot(sample_data_shannon, aes(x = x, y = y, fill = mean_annual)) +
        geom_raster() +
        scale_fill_viridis_c(
            name = "Annual Mean Shannon Diversity",
            na.value = "grey80",
            option = "plasma",
            trans = "sqrt"
        ) +
        labs(
            title = "Plasma Color Scale Reference (Shannon Diversity)",
            subtitle = paste("Example from", sample_genus_shannon),
            x = "Longitude",
            y = "Latitude"
        ) +
        theme_minimal() +
        theme(
            plot.title = element_text(size = 16, hjust = 0.5),
            plot.subtitle = element_text(size = 14, hjust = 0.5),
            legend.title = element_text(size = 12),
            legend.text = element_text(size = 10)
        )
    
    ggsave(file.path(wd_out, "plasma_colorscale_reference_shannon.png"), 
           p_colorscale_shannon, width = 12, height = 8, dpi = 300)
    
    # Create Shannon summary document
    summary_df_shannon <- data.frame(
        File_Number = 1:length(fnames_shannon),
        Genus = all_genera_shannon,
        Input_File = basename(fnames_shannon),
        Multiplot_Number = ceiling((1:length(fnames_shannon)) / plots_per_multiplot_shannon),
        stringsAsFactors = FALSE
    )
    
    # Save Shannon summary
    write.csv(summary_df_shannon, file.path(wd_out, "genus_mapping_summary_shannon.csv"), row.names = FALSE)
    
    # Print final Shannon summary
    cat(paste(rep("=", 80), collapse=""), "\n")
    cat("SHANNON MAPPING COMPLETE\n")
    cat(paste(rep("=", 80), collapse=""), "\n")
    cat("Total Shannon genera processed:", length(all_genera_shannon), "\n")
    cat("Total Shannon multiplots created:", n_multiplots_shannon, "\n")
    cat("Plots per Shannon multiplot:", plots_per_multiplot_shannon, "\n")
    cat("Shannon color scheme: Plasma (distinguishable from richness)\n")
    cat("Output directory:", wd_out, "\n")
    
    cat("\nShannon files created:\n")
    multiplot_files_shannon <- list.files(wd_out, pattern = "genus_maps_shannon_multiplot_.*\\.png$")
    for(file in multiplot_files_shannon) {
        cat("- ", file, "\n")
    }
    cat("- plasma_colorscale_reference_shannon.png\n")
    cat("- genus_mapping_summary_shannon.csv\n")
    
    cat("\nFirst 10 Shannon genera processed:\n")
    print(head(summary_df_shannon[, c("Genus", "Multiplot_Number")], 10))
    
    if(length(all_genera_shannon) != 272) {
        cat("\nWARNING: Expected 272 Shannon genera but processed", length(all_genera_shannon), "\n")
    } else {
        cat("\nSUCCESS: All 272 Shannon genera mapped with plasma color scheme!\n")
    }
    
} else {
    cat("No Shannon diversity files found in:", wd_in, "\n")
    cat("Looking for files matching pattern: marginal_seas_filtered_shannon_*.csv\n")
}

## ====================================================================================================
## END OF SCRIPT
## ====================================================================================================