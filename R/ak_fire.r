# Fetches ArcticDEM tiles and builds the 250 m DEM + hillshade for Alaska.
# Configuration, paths and data loading live in common.R.
#
# The ABoVE reference grid is deliberately NOT used as a projection template: it
# is 1 km and spans AK + western Canada, so passing it to project() would resample
# the DEM to 1 km over a mostly-empty continental extent. See common.R for the CRS
# actually used.

library(terra)
library(rstac)
library(here)

source(here("R", "common.R"))

# ---- Helpers -----------------------------------------------------------------

# Utility, not part of the pipeline: lists ArcticDEM collection ids so the
# resolution/version in dem_collection can be confirmed before a download.
list_arcticdem_collections <- function(url = stac_url) {
  cols <- stac(url) |> collections() |> get_request()
  vapply(cols$collections, function(x) x$id, character(1))
}

# assets_download() writes each asset to output_dir + the path part of its URL,
# so the local destination is predictable before anything is fetched. That is what
# makes resuming possible.
asset_local_path <- function(href, outdir) {
  file.path(outdir, sub("^[a-zA-Z][a-zA-Z0-9+.-]*://[^/]+/", "", href))
}

# A tile counts as present only if it opens AND its last block reads. Reading the
# bottom-right cell is what forces that block off disk, which is exactly what a
# connection dropped mid-transfer leaves behind. Catches truncation, not silent
# corruption.
tile_ok <- function(f) {
  if (!file.exists(f) || file.size(f) == 0) return(FALSE)
  r <- try(terra::rast(f), silent = TRUE)
  if (inherits(r, "try-error")) return(FALSE)
  probe <- try(r[terra::nrow(r), terra::ncol(r)], silent = TRUE)
  !inherits(probe, "try-error")
}

# Download ArcticDEM mosaic tiles covering landscape_poly.
#
# Only the DEM asset is fetched: the published hillshade/browse images are baked
# at a fixed illumination, whereas shade() lets the sun angle be tuned to suit
# the visualisation.
#
# Resumable. rstac cannot resume on its own -- asset_download() hands
# `overwrite` straight to httr::write_disk(), so overwrite = FALSE errors on the
# first file already on disk and overwrite = TRUE refetches everything. So the
# item list is filtered to the tiles actually missing, and they are fetched one
# at a time: a dropped connection then costs one tile rather than the whole run,
# and re-running picks up where it stopped.
download_dem <- function(landscape_poly,
                         collection      = dem_collection,
                         outdir          = dem_dir,
                         asset_names     = "dem",
                         verify          = TRUE,
                         retries         = 2,
                         clip_to_polygon = TRUE) {

  marker <- file.path(outdir, ".download_complete")
  if (file.exists(marker)) {
    message("ArcticDEM download already complete in ", outdir, ". Skipping.")
    return(invisible(outdir))
  }

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  # Buffer the real polygon (not its extent) before transforming: its many
  # vertices keep the lon/lat bounding box honest, where a 4-corner rectangle
  # would bow inwards and under-cover the north edge
  bbox_4326 <- terra::buffer(landscape_poly, width = 1000) |>
    terra::project("EPSG:4326") |>
    terra::ext() |>
    as.vector()
  # terra returns xmin, xmax, ymin, ymax; STAC wants xmin, ymin, xmax, ymax
  bbox_4326 <- unname(bbox_4326[c("xmin", "ymin", "xmax", "ymax")])

  # Query the ArcticDEM STAC catalogue for mosaic tiles intersecting the bbox
  stac_query <- stac_search(
    stac(stac_url),
    collections = collection,
    bbox = bbox_4326
  ) |>
    post_request() |>
    items_fetch()

  if (length(stac_query$features) == 0) {
    stop("No ", collection, " tiles returned for the requested bbox.")
  }

  # Asset keys vary between collections and versions, so check rather than assume
  available <- rstac::items_assets(stac_query)
  absent    <- setdiff(asset_names, available)
  if (length(absent) > 0) {
    warning("Asset(s) not in ", collection, ": ", paste(absent, collapse = ", "),
            "\nAvailable: ", paste(available, collapse = ", "))
    asset_names <- intersect(asset_names, available)
  }
  if (length(asset_names) == 0) {
    stop("None of the requested assets exist in ", collection)
  }

  feats <- stac_query$features

  # A bbox query over Alaska is a loose fit: the panhandle and Seward Peninsula
  # make the state a diagonal, so many returned tiles are entirely ocean or
  # Canada. Drop any tile whose footprint does not actually meet the polygon --
  # fewer tiles to fetch, and a smaller vrt to mosaic later.
  if (clip_to_polygon) {
    poly_ll <- terra::project(landscape_poly, "EPSG:4326")
    tile_polys <- do.call(rbind, lapply(feats, function(f) {
      b <- as.numeric(f$bbox)
      terra::as.polygons(terra::ext(b[1], b[3], b[2], b[4]), crs = "EPSG:4326")
    }))
    hits <- terra::is.related(tile_polys, poly_ll, "intersects")
    message(sprintf("bbox query returned %d tile(s); %d actually meet the polygon",
                    length(feats), sum(hits)))
    feats <- feats[hits]
    stac_query$features <- feats
    if (length(feats) == 0) stop("No tiles intersect the landscape polygon.")
  }

  # Which tiles are still needed
  dest <- lapply(feats, function(f) {
    vapply(asset_names, function(a) {
      href <- f$assets[[a]]$href
      if (is.null(href)) NA_character_ else asset_local_path(href, outdir)
    }, character(1))
  })

  needed <- vapply(dest, function(paths) {
    if (anyNA(paths)) return(TRUE)
    if (verify) any(!vapply(paths, tile_ok, logical(1)))
    else        any(!file.exists(paths))
  }, logical(1))

  message(sprintf("%d tile(s) in query | %d already present | %d to fetch",
                  length(feats), sum(!needed), sum(needed)))

  if (!any(needed)) {
    file.create(marker)
    message("Nothing to do.")
    return(invisible(outdir))
  }

  # Clear anything that exists but failed verification, or write_disk() will
  # refuse to write over it
  stale <- unlist(dest[needed])
  stale <- stale[!is.na(stale) & file.exists(stale)]
  if (length(stale) > 0) {
    message("Removing ", length(stale), " incomplete tile(s) before refetching.")
    unlink(stale)
  }

  idx <- which(needed)
  failed <- character(0)
  for (k in seq_along(idx)) {
    i  <- idx[k]
    id <- feats[[i]]$id
    one <- stac_query
    one$features <- feats[i]

    done <- FALSE
    for (attempt in seq_len(retries + 1)) {
      done <- tryCatch({
        assets_download(one, asset_names = asset_names, output_dir = outdir,
                        create_json = FALSE)
        TRUE
      }, error = function(e) {
        message(sprintf("  [%d/%d] %s attempt %d failed: %s",
                        k, length(idx), id, attempt, conditionMessage(e)))
        FALSE
      })
      if (done) break
    }

    if (done) {
      message(sprintf("  [%d/%d] %s", k, length(idx), id))
    } else {
      failed <- c(failed, id)
    }
  }

  if (length(failed) == 0) {
    file.create(marker)
    message("Download complete: ", length(idx), " tile(s) fetched.")
  } else {
    warning(length(failed), " tile(s) still missing after ", retries + 1,
            " attempts: ", paste(utils::head(failed, 10), collapse = ", "),
            if (length(failed) > 10) ", ..." else "",
            "\nRe-run download_dem() to retry only those.")
  }

  invisible(outdir)
}

# ---- Read vector data --------------------------------------------------------

ak <- load_ak()

# Restore this only if the outputs need to align with other ABoVE products, and
# pass crs(above_study_domain) -- never the raster itself -- to project():
# above_study_domain_file <- list.files(
#   "//10.60.2.10/FF_Lab/project_data/na_boreal/data_sets/ABoVE_reference_grid_v2_1527/data",
#   pattern = "\\.tif$", full.names = TRUE)[1]
# above_study_domain <- rast(above_study_domain_file)

# ---- DEM and hillshade -------------------------------------------------------

# The ArcticDEM tiles are ~3 GB of scaffolding; ak_dem.tif and ak_hillshade.tif
# are the durable products, so the tiles can be deleted once both exist. The
# download is therefore gated on BOTH stages: it runs only when neither the
# processed rasters nor a complete tile set is present. A re-run after the tiles
# have been deleted will not refetch them.
#
# "Complete tile set" is the marker AND at least one tile still on disk, so
# emptying the directory but leaving the marker behind is not mistaken for done.
dem_marker <- file.path(dem_dir, ".download_complete")
have_tiles <- file.exists(dem_marker) &&
  length(list.files(dem_dir, pattern = "_dem\\.tif$", recursive = TRUE)) > 0
have_processed <- file.exists(dem_file) && file.exists(hillshade_file)

if (!have_processed && !have_tiles) {
  download_dem(ak)
  have_tiles <- file.exists(dem_marker)
}

# Both products come out of one block so the hillshade can be derived from the
# unmasked DEM, which avoids an NA hairline along the coast
if (!have_processed) {

  dem_tiles <- list.files(dem_dir, pattern = "_dem\\.tif$",
                          full.names = TRUE, recursive = TRUE)
  if (length(dem_tiles) == 0) {
    # Reachable when the tiles have been deleted but only one of the two
    # processed rasters survives, so say what to do rather than just failing.
    stop("Need to build:\n",
         if (!file.exists(dem_file)) paste0("  - ", dem_file, "\n"),
         if (!file.exists(hillshade_file)) paste0("  - ", hillshade_file, "\n"),
         "but no ArcticDEM tiles are present under ", dem_dir, ".\n",
         "Either restore the tiles, or call download_dem(ak) to refetch them ",
         "(~3 GB). If a stale ", basename(dem_marker),
         " marker is present, delete it first.")
  }

  dem_vrt <- vrt(dem_tiles)
  names(dem_vrt) <- "elevation"

  # Crop in the DEM's native CRS *before* projecting, otherwise the whole of
  # Alaska is warped at source resolution first. The buffer keeps terrain()
  # from seeing an artificial edge inside the area of interest.
  ak_native <- terra::buffer(ak, width = 5000) |> terra::project(crs(dem_vrt))
  dem_crop  <- crop(dem_vrt, ak_native)

  # "average" resamples correctly for a large downsampling factor; bilinear
  # would point-sample and alias the terrain
  dem_proj <- project(dem_crop, target_crs, res = dem_res,
                      method = "average", threads = TRUE)

  if (!file.exists(dem_file)) {
    ak_dem <- mask(crop(dem_proj, ak), ak)
    writeRaster(ak_dem, dem_file, gdal = gtiff_flt, overwrite = TRUE)
  } else {
    ak_dem <- rast(dem_file)
  }

  if (!file.exists(hillshade_file)) {
    # Exaggerate before deriving slope, not after: terrain() reads slope from the
    # elevation-to-cell-size ratio, so scaling elevation is what steepens it.
    dem_z  <- if (hillshade_z == 1) dem_proj else dem_proj * hillshade_z
    slope  <- terrain(dem_z, v = "slope",  unit = "radians")
    aspect <- terrain(dem_z, v = "aspect", unit = "radians")

    # Multi-directional hillshade: averaging several sun azimuths reads better
    # than a single light source and avoids the flat-facing-away look.
    # A vector of directions gives one layer per azimuth, so mean() collapses
    # them. Do not pipe through rast() first -- that is the template constructor
    # and silently drops the values.
    # terrain() also yields "aspect"/"slope"/"TRI" directly if wanted as layers.
    ak_hillshade <- shade(slope, aspect, angle = 35,
                          direction = c(225, 270, 315, 360)) |>
      mean()
    names(ak_hillshade) <- "hillshade"

    ak_hillshade <- mask(crop(ak_hillshade, ak), ak)
    writeRaster(ak_hillshade, hillshade_file, gdal = gtiff_flt, overwrite = TRUE)
  } else {
    ak_hillshade <- rast(hillshade_file)
  }

} else {
  ak_dem       <- rast(dem_file)
  ak_hillshade <- rast(hillshade_file)
}

# ---- Fire perimeters ---------------------------------------------------------

# Read, year-filter, validate, clip to the mapped landmass and measure areas.
# One definition, in common.R, shared with ak_fire_anim.R.
histfire_yr <- load_fires(ak)

# Annual totals, dissolving within each year so overlapping perimeters from the
# same year are not double-counted. Cached under ak_data/derived/.
annual <- get_annual_burned_area(histfire_yr)
