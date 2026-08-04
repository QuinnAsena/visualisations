# Shared configuration and data loading for the Alaska fire visualisation.
#
# Everything more than one script needs lives here so the fire dataset has
# exactly ONE definition. It previously had three: ak_fire.r, ak_fire_anim.R and a
# since-deleted prep script each read, filtered and validated the perimeters
# independently, which meant the coastline clip below would have had to be
# repeated in three places to stay consistent.

library(terra)
library(here)

# ---- projection and resolution ----------------------------------------------

# All three AK vector layers are natively EPSG:3338, which is also the right
# projection for an Alaska map.
target_crs <- "EPSG:3338"

# Analysis resolution. AK spans ~1360 x 1800 km, so a 2k-wide frame resolves
# ~700 m/px; 250 m leaves headroom for zooms without an unusable raster.
dem_res <- 250

# Animation raster resolution: already finer than the video can resolve, and a
# quarter of the cells of dem_res.
anim_res <- 500

first_year <- 1980
last_year  <- 2020

# Vertical exaggeration applied to the DEM before slope/aspect. At 250 m real
# relief is heavily smoothed and shading reads flat at z = 1. Affects
# ak_hillshade.tif only -- ak_dem.tif stays true elevation.
hillshade_z <- 3

# ---- paths -------------------------------------------------------------------

dem_collection <- "arcticdem-mosaics-v4.1-32m"
dem_dir        <- here("ak_data", "dem", "arcticdem_32m")
dem_file       <- here("ak_data", "ak_dem.tif")
hillshade_file <- here("ak_data", "ak_hillshade.tif")
fire_dsn       <- here("ak_data", "historic_fire", "raw_data", "fire")
derived_dir    <- here("ak_data", "derived")
stac_url       <- "https://stac.pgc.umn.edu/api/v1/"

# PREDICTOR=3 is the floating-point predictor (elevation, hillshade);
# PREDICTOR=2 is the integer one (year, counts).
gtiff_flt <- c("COMPRESS=DEFLATE", "PREDICTOR=3", "TILED=YES")
gtiff_int <- c("COMPRESS=DEFLATE", "PREDICTOR=2", "TILED=YES")

# ---- loaders -----------------------------------------------------------------

read_vect <- function(path, crs_out = target_crs) {
  v <- terra::vect(path)
  if (!terra::same.crs(v, crs_out)) v <- terra::project(v, crs_out)
  v
}

load_ak       <- function() read_vect(here("ak_data", "AK_no_islands.shp"))
load_interior <- function() read_vect(here("ak_data", "AK_interior.shp"))

# Fire perimeters for the study period, clipped to the mapped landmass.
#
# 28 of the 3550 perimeters do not touch AK_no_islands. They are NOT bad
# coordinates: nearly all are real fires on Kodiak (Moser Bay, Larsen Bay,
# Karluk, Old Harbor), Afognak, St Lawrence Island and the Aleutians -- islands
# this basemap deliberately excludes -- plus one east of the 141st meridian, in
# Yukon. Median distance from the landmass is 130 km and none is an edge sliver,
# so nothing looks like a digitising error; the basemap simply has nowhere to put
# them, and they render as fire floating in the sea.
#
# crop() rather than a drop-the-polygon filter, because it also trims the coastal
# overhang of perimeters that straddle the shoreline, so reported areas match the
# pixels actually drawn. Cost: 30,885 ha, 0.17% of total burned area.
load_fires <- function(ak = load_ak(), clip = TRUE) {
  hf <- read_vect(file.path(fire_dsn, "AK_fire_location_polygons.shp"))
  hf$FIREYEAR <- as.numeric(hf$FIREYEAR)
  hf <- hf[!is.na(hf$FIREYEAR) &
             hf$FIREYEAR >= first_year & hf$FIREYEAR <= last_year, ]
  hf <- makeValid(hf)
  if (clip) hf <- crop(hf, ak)
  # Measured from geometry, not the stored Shape_Area, which is in the units of
  # whatever CRS the shapefile was authored in and is stale after any clip.
  hf$area_ha <- expanse(hf, unit = "ha")
  hf
}

# Annual burned area, dissolving within each year first so overlapping perimeters
# from the same year are not double-counted.
annual_burned_area <- function(hf) {
  yrs <- first_year:last_year
  ha <- vapply(yrs, function(y) {
    s <- hf[hf$FIREYEAR == y, ]
    if (nrow(s) == 0) return(0)
    sum(expanse(aggregate(s), unit = "ha"))
  }, numeric(1))
  data.frame(year = yrs, burned_ha = ha)
}

# Cached form: the 41 per-year dissolves take a couple of minutes, and the render
# needs the same numbers on every run.
get_annual_burned_area <- function(hf = NULL, cache = file.path(derived_dir,
                                                  "annual_burned_area.csv")) {
  if (file.exists(cache)) return(utils::read.csv(cache))
  if (is.null(hf)) hf <- load_fires()
  a <- annual_burned_area(hf)
  dir.create(dirname(cache), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(a, cache, row.names = FALSE)
  a
}
