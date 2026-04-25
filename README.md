# Forest Structure Metrics from Segmented LiDAR

This repository contains two R scripts to process airborne LiDAR tile data in which individual trees have been segmented and labeled by species. The first script computes a comprehensive set of forest structural metrics and exports them as multi‑layer GeoTIFF rasters at multiple pixel resolutions. The second script uses species‑labeled tree points to generate vertically stratified “species structure” layers, again as multi‑band GeoTIFFs.

## Input requirements

- **LiDAR format**: `.las` or `.laz` tiles (ASPRS LAS format).
- **Point labeling**:
  - Script 1: Generic ground / nonground classification (e.g., `Classification == 2` for ground).
  - Script 2: At least one numeric `Species` attribute and a `Tree` identifier per segmented tree crown.
- **Coordinate system**: All tiles should be in the same projected CRS (e.g., UTM).

## Script 1: Forest structural metrics

**File**: `01_forest_metrics.R`

### Purpose
Compute pixel‑wise and tile‑wise forest structural metrics from normalized LiDAR point clouds and export them as multi‑layer GeoTIFFs at base and aggregated resolutions (e.g., 1, 2, 5, 10, 20 m).

### Key metrics computed

- **Height distribution per pixel**
  - `cv_height`: coefficient of variation of Z within the pixel.
  - `canopy_cover`: proportion of points above a height threshold (e.g., 2 m).
  - `gap_fraction`: 1 − `canopy_cover`.
  - `canopy_relief_ratio`: `(zmean − zmin)/(zmax − zmin)` for the pixel.
  - `understory_density`: fraction of canopy returns that lie between the height threshold and a lower upper limit (e.g., 5 m).
- **CHM‑based layer**
  - `rugosity_chm`: local standard deviation of a 3×3 m window on the canopy height model (CHM).
- **Tile‑level scalar metrics**
  - `prop_above_mean`: proportion of all tree returns above the overall mean tree height for the 20×20 m tile.

### Output files

For each input tile, the script writes:

- `*_res1.tif` (1 m base resolution).
- `*_res2.tif`, `*_res5.tif`, `*_res10.tif`, `*_res20.tif` (aggregated resolutions).

Each GeoTIFF contains one band per metric, plus the CHM‑rugosity and tile‑scalar layers.

### Configuration

Edit the top block of `01_forest_metrics.R`:

```r
output_path          <- "path/to/Forest_Structure/"
target_dirs          <- c("path/to/tiles/dir1", "path/to/tiles/dir2")
base_res             <- 1
derived_res          <- c(2, 5, 10, 20)
height_threshold     <- 2
understory_threshold <- 5
dtm_resolution       <- 1
dtm_method           <- tin()
```

You can adjust:
- `target_dirs`: list of directories containing your `.las`/`.laz` files.
- `base_res` and `derived_res`: the raster resolutions (m) to export.
- `height_threshold` and `understory_threshold`: define canopy vs. understory in m.

## Script 2: Species‑explicit structure layers

**File**: `02_species_structure.R`

### Purpose
Leverage species‑labeled tree segments to build vertically stratified “species structure” rasters. Each layer corresponds to a vertical slice (e.g., 1 m height intervals) and encodes the dominant species (or presence) within raster cells.

### How it works

1. **Terrain normalization**
   - Fitting a DTM at 0.5 m resolution and normalizing tree heights.
2. **Tree polygonization**
   - For each tree segment (`Tree` ID), the script computes a concave hull around its 2D projection (X, Y) and represents each tree as a polygon.
3. **Vertical slicing**
   - Tree points are partitioned into height slices (e.g., 0–1 m, 1–2 m, etc.).
4. **Rasterization per slice**
   - For each slice, all tree polygons of that slice are rasterized into a single layer using the `Species` attribute as the pixel value.
5. **Multi‑layer output**
   - All slice layers are combined into one multi‑band GeoTIFF, one band per height interval.

### Output files

- `species_structure_*.tif` (one per input tile).
- Each band is named `sliceX_low-highm`, where `X` is the slice index and `low–high` is the height interval.

### Configuration

Edit the top block of `02_species_structure.R`:

```r
input_path  <- "path/to/segmented/tiles/"
output_path <- "path/to/Forest_Structure/Species/"
```

You can additionally adjust:
- `resolution = 0.5` for the DTM and `resolution = 0.2` for the raster template.
- The `slice_breaks` sequence if you want coarser or finer vertical bins.

## Dependencies

Both scripts require the following R packages:

- `lidR` (≥ 4.x) – LiDAR I/O and processing.
- `terra` – raster creation, aggregation, and export.
- `sf` – spatial features and polygons.
- `concaveman` – fast concave hulls from R.
- `dplyr` and `e1071` (used only in script 1, but left in for completeness).

Install via:

```r
install.packages(c("lidR", "terra", "sf", "dplyr", "concaveman"))
```

## Example usage

Assuming stored in a project folder with `01_forest_metrics.R` and `02_species_structure.R`:

```bash
cd /path/to/your/project
Rscript 01_forest_metrics.R
Rscript 02_species_structure.R
```

Adjust the `output_path` and `target_dirs`/`input_path` to match your tree‑segmented LiDAR structure on your filesystem.

## Notes

- Tiles are assumed to represent exactly 20×20 m patches aligned on a common grid.
- The species‑structure script assumes that the `Species` attribute is numeric and meaningful (e.g., species‑ or functional‑group IDs).
- If you change the slicing logic or resolutions, make sure to keep the overall height range and coordinate alignment consistent across tiles so that downstream mosaicking works.
- Please adapt file paths, project names, and citation details to your specific study or dataset.

## Project context
This code base is developed as part of **[POINTR](https://helmholtz-imaging.de/project/pointr/)**, a **[Helmholtz Imaging Project](https://helmholtz-imaging.de/projects/)** titled *Mapping Boreal Forest Change Using 3D Radar and Point Cloud Data*.

Within that context, the workflows in this repository support the extraction of forest-structure and species-explicit raster products from segmented LiDAR point clouds. These products contribute to the broader POINTR objective of detecting subtle structural changes in boreal forests by integrating point cloud information, radar remote sensing, and ecological analysis.
