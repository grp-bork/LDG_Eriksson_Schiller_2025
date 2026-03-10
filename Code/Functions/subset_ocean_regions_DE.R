# This function is used to subset oceanic regions of interest.

# Available ocean regions to subset:
# [1] "Southern Ocean"
# [2] "South Atlantic Ocean"
# [3] "South Pacific Ocean"
# [4] "North Pacific Ocean"
# [5] "South China and Easter Archipelagic Seas"
# [6] "Indian Ocean"
# [7] "Mediterranean Region"
# [8] "Baltic Sea"
# [9] "North Atlantic Ocean"
# [10] "Arctic Ocean"

subset_ocean_regions <- function(
    data # path to the csv file or a dataframe with coordinates indicated by "decimalLongitude" and "decimalLatitude"
){
        # Load libraries
        library(sf)

        # Index and corresponding marine regions
        ind <- c(
        1, # "Southern Ocean"
        2, # "South Atlantic Ocean"
        3, # "South Pacific Ocean"
        4, # "North Pacific Ocean"
        5, # "South China and Easter Archipelagic Seas"
        6, # "Indian Ocean"
        7, # "Mediterranean Region"
        8, # "Baltic Sea"
        9, # "North Atlantic Ocean"
        10) # "Arctic Ocean"
        #
        ind_names <- c(
            "Southern Ocean",
            "South Atlantic Ocean",
            "South Pacific Ocean",
            "North Pacific Ocean",
            "South China and Easter Archipelagic Seas",
            "Indian Ocean",
            "Mediterranean Region",
            "Baltic Sea",
            "North Atlantic Ocean",
            "Arctic Ocean"
        )

        # Read shapefile with geometries of ocean basins
        sf.ocean <- st_read("Data/Ocean_regions_shapefile/GOaS_v1_20211214/goas_v01.shp")

        # Convert diazotroph dataset to sf object
        sf_data <- sf::st_as_sf(data, coords = c("decimalLongitude", "decimalLatitude")) # Convert foreign object to sf object

        # Get geometries from each ocean basin
        ocean_geoms <- st_geometry(sf.ocean)

        ## Loop through to subset each marine region
        l_regions <- list()
        for(i in seq_along(ind)){

            # Print progress
            message( paste0("Subsettin ", ind_names[i]) )

            # Select ocean region of interest
            region <- ocean_geoms[[ind[i]]]

            # Retrieve number of observation
            good_points <- st_filter(sf_data, region)

            # Convert back to dataframe
            coords <- st_coordinates(good_points)
            coords <- as.data.frame(coords)
            names(coords) <- c("decimalLongitude", "decimalLatitude")
            df <- as.data.frame(good_points)
            df <- data.frame(df[, -ncol(df)])
            if(nrow(df) > 0){
                df$marine_region <- ind_names[i]
            }
            df <- cbind(df, coords)

            # Save in list
            if(nrow(df) > 0){
                l_regions[[i]] <- df
            }

        } # close loop across ind

        # Merge list to dataframe
        df_final <- do.call("rbind", l_regions) 

        # Return object
        return(df_final)
}
