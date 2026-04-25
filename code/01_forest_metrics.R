library(lidR)
library(e1071)
library(terra)
library(sf)
library(dplyr)

# --- Configuration ---
output_path          <- ".../Forest_Structure"	# set your output path

if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)
target_dirs <- c( # set your input path(s)
  "...",
  "...",
  "..."
)

# 2. List all .las and .laz files across all directories
# pattern = "\\.la[sz]$": matches files ending in .las or .laz (case-insensitive)
# full.names = TRUE: returns the absolute path needed for processing
all_files <- list.files(
  path = target_dirs, 
  pattern = "\\.la[sz]$", 
  full.names = TRUE, 
  ignore.case = TRUE
)
las_files            <- all_files
base_res             <- 1
derived_res          <- c(2, 5, 10, 20)
height_threshold     <- 2
understory_threshold <- 5
dtm_resolution       <- 1
dtm_method           <- tin()

# --- Pixel-Metric-Functions ---------------------------------------------------
compute_forest_metrics <- function(Z, Intensity,
                                   ht = height_threshold,
                                   ut = understory_threshold) {
  if (length(Z) == 0) return(NULL)
  
  zmax  <- max(Z,  na.rm = TRUE)
  zmean <- mean(Z, na.rm = TRUE)
  zsd   <- sd(Z,   na.rm = TRUE)
  zmin  <- min(Z,  na.rm = TRUE)
  
  cv_h    <- ifelse(zmean == 0, NA_real_, zsd / zmean)
  c_cover <- sum(Z > ht) / length(Z)
  crr     <- ifelse(zmax == zmin, NA_real_, (zmean - zmin) / (zmax - zmin))  # CRR per pixel
  u_dens  <- ifelse(sum(Z > ht) == 0, 0,
                    sum(Z > ht & Z < ut) / sum(Z > ht))
  
  std <- stdmetrics_z(Z)
  
  c(
    list(
      cv_height          = cv_h,
      canopy_cover       = c_cover,
      gap_fraction       = 1 - c_cover,
      canopy_relief_ratio = crr,          # spatial, 1m resolution
      understory_density = u_dens
    ),
    std
  )
}

# --- Processing-Loop --------------------------------------------------
for (file_path in las_files) {
  
  las <- tryCatch(
    readLAS(file_path, filter = "-drop_z_below 0"),
    error = function(e) {
      cat(sprintf("Error on reading: %s — %s\n", basename(file_path), e$message))
      NULL
    }
  )
  if (is.null(las)) next
  
  base_name <- tools::file_path_sans_ext(basename(file_path))
  cat(sprintf("\n--- Starting Processing: %s ---\n", base_name))
  
  # 1. terrain correction
  las_ground <- filter_poi(las, Classification == 2)
  if (nrow(las_ground@data) == 0) {
    cat("    -> WARNING: No ground points. Continuing.\n")
    rm(las); gc(); next
  }
  
  dtm      <- rasterize_terrain(las_ground, res = dtm_resolution,
                                algorithm = dtm_method)
  las_norm <- normalize_height(las, dtm)
  las_norm <- filter_poi(las_norm, Z >= 0)
  
  if (nrow(las_norm@data) == 0) {
    cat("    -> WARNING: Point cloud empty after normalizing.\n")
    rm(las, las_ground, dtm); gc(); next
  }
  
  # 2. crop to exact 20×20m
  pc_xmin <- min(las_norm@data$X, na.rm = TRUE)
  pc_ymin <- min(las_norm@data$Y, na.rm = TRUE)
  
  perfect_bbox <- st_bbox(
    c(xmin = pc_xmin, ymin = pc_ymin,
      xmax = pc_xmin + 20, ymax = pc_ymin + 20),
    crs = st_crs(las_norm)
  )
  las_final <- clip_roi(las_norm, perfect_bbox)
  
  if (nrow(las_final@data) == 0) {
    cat("    -> WARNING: Point cloud empty after cropping.\n")
    rm(las, las_ground, dtm, las_norm); gc(); next
  }
  
  # 3. Pixel-metrics (1m raster)
  cat(sprintf("  -> Pixel-metrics at %sm resolution\n", base_res))
  
  base_raster <- pixel_metrics(
    las_final,
    func   = ~compute_forest_metrics(Z, Intensity),
    res    = base_res,
    origin = c(pc_xmin, pc_ymin)
  )
  
  e_base      <- terra::ext(pc_xmin, pc_xmin + 20, pc_ymin, pc_ymin + 20)
  base_raster <- terra::crop(base_raster, e_base)
  
  # 4. CHM-based metrics Metriken
  cat("  -> Calculate CHM-based metrics\n")
  
  chm <- tryCatch(
    grid_canopy(las_final, res = 1, algorithm = p2r()),
    error = function(e) {
      cat(sprintf("    -> WARNING: CHM could not be computed: %s\n", e$message))
      NULL
    }
  )
  
  # Rugosity: focal SD (3×3m Fenster) → spatial layer (1m resolution)
  rugosity_layer <- if (!is.null(chm)) {
    chm_terra <- terra::rast(chm)   # RasterLayer → SpatRaster
    rug <- terra::focal(chm_terra, w = matrix(1, 3, 3), fun = sd, na.rm = TRUE)
    terra::resample(rug, base_raster[[1]], method = "bilinear")
  } else {
    r <- terra::rast(base_raster[[1]])
    terra::values(r) <- NA_real_
    r
  }
  names(rugosity_layer) <- "rugosity_chm"
  
  # 5. Tile-level metrics as constant layer
  cat("  -> Calculate tile-level metrics\n")
  Z_all <- las_final@data$Z
  
  # Canopy relief ratio
  z_range <- max(Z_all) - min(Z_all)
  crr     <- ifelse(z_range == 0, NA_real_,
                    (mean(Z_all) - min(Z_all)) / z_range)
  
  # Proportion of returns above mean height
  prop_above_mean <- sum(Z_all > mean(Z_all)) / length(Z_all)
  
  tile_scalars <- list(
    prop_above_mean     = prop_above_mean
  )
  
  tile_layers <- terra::rast(lapply(tile_scalars, function(val) {
    r <- terra::rast(base_raster[[1]])
    terra::values(r) <- val
    r
  }))
  names(tile_layers) <- names(tile_scalars)
  
  # 6. Merge layer
  base_raster <- c(base_raster, rugosity_layer, tile_layers)
  
  # 7. Save base raster
  terra::writeRaster(
    base_raster,
    file.path(output_path, paste0(base_name, "_res", base_res, ".tif")),
    overwrite = TRUE
  )
  
  # 8. Aggregate raster
  pixel_layer_names <- names(base_raster)[!names(base_raster) %in%
                                            names(tile_scalars)]
  pixel_stack       <- base_raster[[pixel_layer_names]]
  
  for (res_agg in derived_res) {
    agg_raster <- terra::aggregate(pixel_stack,
                                   fact = res_agg / base_res,
                                   fun  = "mean", na.rm = TRUE)
    
    # Build tile-level layer on new resolution 
    tile_layers_agg <- terra::rast(lapply(tile_scalars, function(val) {
      r <- terra::rast(agg_raster[[1]])
      terra::values(r) <- val
      r
    }))
    names(tile_layers_agg) <- names(tile_scalars)
    
    terra::writeRaster(
      c(agg_raster, tile_layers_agg),
      file.path(output_path, paste0(base_name, "_res", res_agg, ".tif")),
      overwrite = TRUE
    )
  }
  
  rm(las, las_ground, las_norm, las_final, dtm, base_raster, chm, rugosity_layer)
  gc()
}

cat("\n Batch-Processing done.\n")
