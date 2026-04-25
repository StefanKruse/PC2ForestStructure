require(terra)
require(lidR)
require(sf)
require(concaveman)
require(dbscan)

input_path <- "..." # set your input path
output_path <- "..." # set your output path
if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)

las_files <- list.files(input_path, pattern = "\\.laz$", full.names = TRUE)
cat("Found", length(las_files), "LAS files\n")

for (i in 1:length(las_files)) {
  cat("\n=== File", i, "/", length(las_files), "===\n")
  all_rasters <- list()   
  las <- readLAS(las_files[i])
  
  base_name <- tools::file_path_sans_ext(basename(las_files[i]))
  base_name <- gsub("predicted", "species_structure", base_name)
  
  dtm <- rasterize_terrain(las, res = 0.5)
  las <- las - dtm
  las@data$Z[las@data$Z < 0] <- 0
  
  tree <- filter_poi(las, Species > 2)
  if (npoints(tree) == 0) {
    cat("No trees, skipping\n"); next
  }
  
  tree_crs <- st_crs(tree)
  tree_crs<-tree_crs$input
  tree_coords <- tree@data[, c("X", "Y")]
  full_extent <- ext(apply(tree_coords, 2, range))
  r_template <- rast(full_extent, resolution = 0.2, crs = tree_crs)
  
  z_min <- min(tree$Z); z_max <- max(tree$Z)
  slice_breaks <- seq(floor(z_min), ceiling(z_max), by = 1)
  slice_breaks <- slice_breaks[slice_breaks > 0]
  
  for (o in 1:(length(slice_breaks) - 1)) {
    lower_bound <- slice_breaks[o]
    upper_bound <- slice_breaks[o + 1]
    
    slice <- filter_poi(tree, Z >= lower_bound & Z < upper_bound)
    if (npoints(slice) == 0) next
    
    points <- slice@data
    points$Species <- as.factor(points$Species)
    
    coords <- points[, c("X", "Y")]
    coords_df <- data.frame(X = coords[, 1], Y = coords[, 2], 
                            Species = points$Species, Tree = points$Tree)
    pts_sf <- st_as_sf(coords_df, coords = c("X", "Y"), crs = tree_crs)
    
    # Use Tree field for individual trees
    unique_trees <- unique(points$Tree[!is.na(points$Tree)])
    
    all_tree_polys <- list()  # Collect all tree polygons for this slice
    
    for (tree_id in unique_trees) {
      tree_pts_mask <- pts_sf$Tree == tree_id & !is.na(pts_sf$Tree)
      tree_pts <- pts_sf[tree_pts_mask, ]
      
      if (nrow(tree_pts) < 3) next
      
      tree_coords <- st_coordinates(tree_pts)
      
      # FIXED: Higher concavity + length_threshold for smoother polygons
      hull <- concaveman(tree_coords, concavity = 1.5, length_threshold = 0.1)
      
      # Robust closure
      if (nrow(hull) > 3) {
        if (!all.equal(hull[1, ], hull[nrow(hull), ], tolerance = 1e-6)) {
          hull <- rbind(hull, hull[1, ])
        }
        
        poly <- tryCatch({
          st_polygon(list(hull))
        }, error = function(e) NULL)
        
        if (!is.null(poly)) {
          tree_poly_sf <- st_sfc(poly, crs = tree_crs) %>% 
            st_sf(species = tree_pts$Species[1], tree_id = tree_id) %>%
            st_make_valid()
          all_tree_polys[[length(all_tree_polys) + 1]] <- tree_poly_sf
        }
      }
    }
    
    if (length(all_tree_polys) > 0) {
      # Combine ALL tree polygons into one MultiPolygon per slice
      all_trees_sf <- do.call(rbind, all_tree_polys)
      
      # Rasterize with SPECIES values (one layer per slice)
      slice_raster <- rasterize(vect(all_trees_sf), r_template, field = "species")
      layer_name <- paste0("slice", o, "_", sprintf("%.1f-%.1fm", lower_bound, upper_bound))
      names(slice_raster) <- layer_name
      all_rasters[[layer_name]] <- slice_raster
      
      cat(sprintf("  Slice %d (%.1f-%.1fm): %d tree polygons\n", 
                  o, lower_bound, upper_bound, nrow(all_trees_sf)))
    }
  }
  
  if (length(all_rasters) > 0) {
    multi_layer_raster <- rast(all_rasters)
    outfile <- file.path(output_path, paste0(base_name, ".tif"))
    
    writeRaster(multi_layer_raster, outfile, 
                gdal = c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES"), 
                overwrite = TRUE)
    
    cat("Saved", nlyr(multi_layer_raster), "slice layers:", basename(outfile), "\n")
  }
}

