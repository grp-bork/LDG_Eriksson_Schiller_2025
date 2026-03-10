# This function is used to subset oceanic regions including marginal seas.

iho_ocean_regions <- function(
    data # path to the csv file or a dataframe with coordinates indicated by "decimalLongitude" and "decimalLatitude"
){
        # Load libraries
        library(sf)

        # Read shapefile with geometries of ocean basins
        sf.ocean <- st_read("/net/sea/work/deriksson/Projects/Prokaryotes_LDG_v3/Data/Shape_files/World_Seas_IHO_v3/World_Seas_IHO_v3.shp")

        # # Create Index and index names
        # ind <- unique(sf.ocean$NAME)
        # ind_names <- unique(sf.ocean$NAME)

        ind <- seq_len(nrow(sf.ocean))  # Use row numbers: 1, 2, 3, ..., 101
        ind_names <- sf.ocean$NAME      # Get the corresponding names

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
