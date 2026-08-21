# Shared loading, palette and ring geometry for the iLand tree-ring disc.
#
# Data: a LOCAL copy of the iLand project's processed summaries, under iland_data/.
# The source of truth is
#   //10.60.2.10/FF_Lab/personal_storage/quinn_storage/iLand-visualisation-ai/
#     data/processed/area_dom/
# copied locally so rendering never depends on the share. Re-copy after any
# preprocess.R run there. The 10 m map tiles (~1.9 GB) are deliberately NOT copied.

library(arrow)
library(here)

iland_dir  <- here("iland_data")
area_dom_f <- file.path(iland_dir, "area_dom", "dominant_species_area.parquet")

# Downscaled GCM climate, summarised to one value per (gcm, ssp, var, year) by
# esa-2026/R/process_gcms.R. Source:
#   //10.60.2.10/FF_Lab/personal_storage/quinn_storage/gcm_summaries
# 1089 files, 2.17 MB. Copied locally because reading them one-by-one over the
# share did not finish in 10 minutes (66 files alone took 44 s), then
# consolidated once into a single parquet.
climate_raw <- file.path(iland_dir, "climate", "raw")
climate_f   <- file.path(iland_dir, "climate", "gcm_annual.parquet")

# Alaska boundaries for the locator map. Deliberately re-declared here rather
# than sourcing R/common.R (the AK fire piece's module): that file sets
# first_year <- 1980 / last_year <- 2020, which would silently clash with this
# piece's 2014-2100. Two small loaders are cheaper than that trap.
ak_shp    <- here("ak_data", "AK_no_islands.shp")
intak_shp <- here("ak_data", "AK_interior.shp")

load_ak       <- function() terra::project(terra::vect(ak_shp), "EPSG:3338")
load_interior <- function() terra::project(terra::vect(intak_shp), "EPSG:3338")

# Landscape footprints for the locator map, derived once from the iLand env
# grids and cached. Same pattern as the climate consolidation: touch the share
# on first run only.
landscape_src <- paste0("//10.60.2.10/FF_Lab/personal_storage/quinn_storage/",
                        "landscape_init_ak_can")
footprint_f   <- file.path(iland_dir, "landscapes", "footprints.gpkg")

# ---- palette -----------------------------------------------------------------
#
# DIVERGES FROM esa-2026 ON PURPOSE. The inherited palette solved dark-mode
# contrast but had never been checked for colour-vision deficiency, and it fails
# badly: under Machado 2009 at full severity its worst all-pairs separation is
# dE 3.8 (white spruce vs mixed spruce) -- effectively one colour. Mixed spruce vs
# mixed deciduous is 4.6. A deuteranope saw roughly three groups, not six.
#
# Designed here against hard gates rather than by eye. Reasoning about it failed
# first: an even lightness ladder scored 0.9, WORSE than the legacy palette,
# because CVD simulation redistributes lightness according to chroma. A search
# over the achievable space showed six colours can reach dE 21.7 unconstrained,
# and 13.9 while keeping species-evocative hues -- so semantics cost almost
# nothing. This palette is the hill-climbed result inside those hue bands:
#
#   all-pairs CVD dE    9.8  (was 3.8)      | normal-vision dE 15.6
#   contrast on #0d0d0d 3.3 min             | max chroma 0.151, off-neon
#   nearest the temperature strip  dE 12.1
#
# Lightness is assigned by ABUNDANCE, not just semantics: the two smallest
# categories are the brightest, because a 2% wedge has to be findable.
#   Potr 12% -> gold at 13.1:1   Bene 2% -> coral at 8.6:1
#   Mixed.deciduous 21% -> dark brown at 3.3:1, deliberately receding
# Putting the big warm category dark and the tiny one bright is also what buys
# the two the lightness separation they need to survive CVD.
#
# Birch as coral rather than white is the one real semantic compromise: gold went
# to aspen, which owns that association more strongly, and coral is what remains
# in the warm band with enough separation from both gold and brown.

SPECIES_COLOURS_DARK <- c(
  Pima            = "#1C716D",  # Picea mariana        black spruce, dark teal
  Pigl            = "#4FC1BB",  # Picea glauca         white spruce, light teal
  Potr            = "#F6D14D",  # Populus tremuloides  aspen, autumn gold
  Bene            = "#F89177",  # Betula neoalaskana   Alaska birch, coral
  Mixed.spruce    = "#659423",  # olive green
  Mixed.deciduous = "#8C5920",  # dark bark brown
  `Not forested`  = "#78909c"   # unused in this piece; kept for completeness
)

SPECIES_LABELS <- c(
  Pima            = "Black spruce",
  Pigl            = "White spruce",
  Potr            = "Trembling aspen",
  Bene            = "Alaska birch",
  Mixed.spruce    = "Mixed spruce",
  Mixed.deciduous = "Mixed deciduous",
  `Not forested`  = "Not forested"
)

# esa-2026 ordering (conifer -> deciduous -> mixed), minus Not forested, which is
# dropped from this piece.
SPECIES_ORDER_6 <- c("Pima", "Pigl", "Potr", "Bene",
                     "Mixed.spruce", "Mixed.deciduous")

SPECIES_COLOURS_DARK_6 <- SPECIES_COLOURS_DARK[SPECIES_ORDER_6]

# GCM_COLOURS is deliberately NOT copied here: nothing in this piece colours by
# climate model. When the GCM x replicate grid arrives it is a one-line lift from
# esa-2026/R/theme.R, which keys both the NorESM2-MM and NorEsm2-MM spellings.

# Surface + chrome, matching ak_fire_anim.R so the two pieces sit together on the
# site. Darker than esa-2026's #191919, so the dark species palette can only gain
# contrast -- verified in the checks below rather than assumed.
# Text lifted one step from the ak_fire_anim values: light grey on near-black is
# harder to read than its contrast ratio suggests. ink_mut carries the caption
# body, and 5.41:1 cleared the 4.5:1 WCAG floor while still reading as strain.
# The three-step hierarchy is preserved -- OKLCH L 1.000 / 0.878 / 0.708.
bg      <- "#0d0d0d"
ink_pri <- "#ffffff"   # 21.0:1
ink_sec <- "#d9d7cc"   # 13.5:1  (was #c3c2b7, 10.9:1)
ink_mut <- "#a3a199"   #  7.5:1  (was #898781,  5.4:1)

# agg_png runs at res = 72 with pointsize = 12, so cex = points / 12.
pt <- function(points) points / 12

# ---- loading -----------------------------------------------------------------

# Scenario family only, one fire condition (default onlysimfalse -- the runs
# where fire kills trees), Not forested dropped, shares renormalised over the
# remaining six categories.
#
# The renormalisation is not cosmetic: total landscape area is a constant
# 58,396.3 ha across all years, but that constant includes Not forested. Dropping
# it leaves a subtotal that varies year to year, so shares MUST be recomputed
# against that subtotal or the rings will not close to a full circle.
load_area_dom <- function(landscape = "landscape_alaska_01_2015-2100scenario",
                          path = area_dom_f, fire = "onlysimfalse") {
  d <- as.data.frame(read_parquet(path))

  d <- d[d$landscape == landscape &
           grepl(fire, d$treatment, fixed = TRUE) &
           d$sp_dom %in% SPECIES_ORDER_6, ]
  if (nrow(d) == 0) stop("No rows for landscape ", landscape)

  key <- paste(d$treatment, d$replicate, d$year_actual, sep = "\r")
  tot <- tapply(d$area_ha, key, sum)
  d$share <- d$area_ha / as.numeric(tot[key])

  d$sp_dom <- factor(d$sp_dom, levels = SPECIES_ORDER_6)
  d$model  <- sub("ssp.*$", "", d$treatment)
  d[order(d$treatment, d$replicate, d$year_actual, d$sp_dom), ]
}

# Wide share matrix (years x species) for one replicate, in SPECIES_ORDER_6.
share_matrix <- function(d, treatment, replicate) {
  s <- d[d$treatment == treatment & d$replicate == replicate, ]
  if (nrow(s) == 0) stop("No rows for ", treatment, " rep ", replicate)
  yrs <- sort(unique(s$year_actual))
  m <- matrix(0, nrow = length(yrs), ncol = length(SPECIES_ORDER_6),
              dimnames = list(yrs, SPECIES_ORDER_6))
  m[cbind(match(s$year_actual, yrs), match(as.character(s$sp_dom), SPECIES_ORDER_6))] <-
    s$share
  m
}

# ---- disturbance -------------------------------------------------------------

# Total-variation distance between consecutive years' compositions: half the sum
# of absolute share changes, so it lands in [0, 1]. 0 = nothing moved, 1 = the
# whole landscape swapped category.
#
# Derived from the SAME numbers that draw the rings, deliberately. The fire tables
# come from separate `onlyfire` simulations, so their replicate 1 is not this
# replicate 1 -- using them here would look convincing and be wrong. This is a
# disturbance index, NOT area burned, and must be labelled as such.
disturbance_index <- function(m) {
  d <- c(0, rowSums(abs(diff(m))) / 2)
  names(d) <- rownames(m)
  d
}

# ---- climate -----------------------------------------------------------------

# The climate files spell it NorESM2-MM; the iLand treatments spell it
# NorEsm2-MM. Joining without normalising silently returns nothing.
canon_gcm <- function(x) sub("NorEsm2-MM", "NorESM2-MM", x, fixed = TRUE)

consolidate_climate <- function(raw = climate_raw, out = climate_f) {
  fs <- list.files(raw, pattern = "[.]parquet$", full.names = TRUE)
  if (length(fs) == 0) {
    stop("No climate parquets under ", raw,
         ".\nCopy them from //10.60.2.10/.../quinn_storage/gcm_summaries first.")
  }
  message("Consolidating ", length(fs), " climate files (one-off)...")
  d <- do.call(rbind, lapply(fs, function(f) as.data.frame(read_parquet(f))))
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  write_parquet(d, out)
  d
}

# Annual precipitation and mean temperature per GCM.
#
# CAVEAT: process_gcms.R takes a spatial mean over the whole downscaled raster,
# so this is a REGIONAL series -- the same climate applies to L01/L02/L03. It is
# not landscape-specific.
load_climate <- function(year_min = 2014, year_max = 2100) {
  if (!file.exists(climate_f)) consolidate_climate()
  d <- as.data.frame(read_parquet(climate_f))
  d <- d[d$year >= year_min & d$year <= year_max, ]
  w <- reshape(d[, c("gcm", "year", "var", "value")],
               idvar = c("gcm", "year"), timevar = "var", direction = "wide")
  names(w) <- sub("^value[.]", "", names(w))
  w$tmean <- (w$tasmax + w$tasmin) / 2
  w[order(w$gcm, w$year), c("gcm", "year", "pr", "tasmax", "tasmin", "tmean")]
}

# Normalisation constants computed ONCE across every GCM in play.
#
# This is the trap this function exists to close: if each disc normalised to its
# own range, NorESM2-MM's +2.0 degC rise would render with exactly the same
# colour spread as UKESM1-0-LL's +6.2, and a GCM grid would be silently
# incomparable while looking perfectly fine.
SCENARIO_LANDSCAPES <- sprintf("landscape_alaska_%02d_2015-2100scenario", 1:3)

# Every cross-run normalisation constant, computed ONCE and passed into each disc.
#
# Named shared_scales, not climate_scales, deliberately. The bug this closes was a
# channel deriving its own scale instead of reading the shared one, and disturbance
# slipped through precisely because a climate-only name gave it nowhere to live.
# One object, one answer to "where does my scale come from".
shared_scales <- function(cl = load_climate(),
                          landscapes = SCENARIO_LANDSCAPES) {
  sc <- climate_scales(cl)
  sc$di_max <- disturbance_max(landscapes)
  sc
}

# Largest disturbance-index value anywhere in the pool: every landscape x model x
# replicate, not just the run being drawn. Pooling one level too narrow is how the
# same bug would reappear when a grid adds landscapes 02 and 03.
disturbance_max <- function(landscapes = SCENARIO_LANDSCAPES) {
  mx <- 0
  for (ls in landscapes) {
    d <- load_area_dom(ls)
    for (tr in unique(d$treatment)) {
      for (rp in unique(d$replicate[d$treatment == tr])) {
        mx <- max(mx, disturbance_index(share_matrix(d, tr, rp)))
      }
    }
  }
  mx
}

climate_scales <- function(cl = load_climate()) {
  # Temperature bounds come from the SMOOTHED series, because the smoothed
  # series is what the strip actually draws. Taking them from raw annual values
  # left 35% of the ramp unreachable (raw span 10.30 degC vs 6.65 drawn) and
  # needlessly compressed every GCM's traverse -- NorESM2-MM covered 27% of the
  # ramp where it should cover 42%.
  sm <- unlist(lapply(split(cl, cl$gcm),
                      function(s) smooth_trend(s$tmean, s$year)))
  list(pr_lo = min(cl$pr), pr_hi = max(cl$pr),
       t_lo  = min(sm),    t_hi  = max(sm),
       # Percentile rank against the POOLED distribution -- see ring_widths().
       # Pooled, not per-GCM, so a systematically drier model keeps systematically
       # thinner rings (mean rank 0.39 / 0.50 / 0.61 for NorESM2 / TaiESM1 /
       # UKESM1) instead of every model being stretched to fill 0-1.
       pr_rank = stats::ecdf(cl$pr))
}

# Sequential cold -> warm across the whole shared range, NOT pivoted at 0 degC.
# Freezing is marked as a TICK on the key instead, so the absolute reference
# survives without the ramp being hostage to it.
#
# Blue -> magenta, and specifically NOT viridis, despite viridis being the
# clearer ramp in isolation. Measured against the species palette over the range
# all three GCMs actually use:
#
#   ramp      lightness traverse   worst dE vs a species colour
#   original       0.039                  10.5   collides
#   viridis        0.391                   1.7   collides badly - its green is
#                                                effectively a species green
#   cividis        0.435                   4.6   collides (yellow vs aspen)
#   plasma         0.415                   2.2   collides (orange vs deciduous)
#   this ramp      ~0.31                  20.4   clears
#
# Blues through magentas are the one region the species palette leaves free, so
# this buys viridis's actual mechanism -- MONOTONE LIGHTNESS, which is what makes
# a small traverse visible -- without the collision. The original ramp changed
# lightness by 0.002 across NorESM2-MM's span, which is why that disc read as a
# single flat colour: it was pure hue rotation.
#
# The cold end is deliberately not near-black. Once the scale was corrected to
# the smoothed range the darkest step became reachable, and against the bg-painted
# gutter it read as a GAP in the strip rather than a cold year (dE 3.7 from
# #0d0d0d). This cold end clears the background at dE 21.6.
TEMP_RAMP <- c("#1d3a6b", "#3560a8", "#8168c2", "#d07fc9", "#ffeef8")

temp_colour <- function(t, sc, n = 256) {
  ramp <- colorRampPalette(TEMP_RAMP)(n)
  f <- (t - sc$t_lo) / (sc$t_hi - sc$t_lo)
  ramp[pmax(1, pmin(n, round(f * (n - 1)) + 1))]
}

# Centred running mean with SHRINKING windows at the edges. Kept as the fallback
# when mgcv is unavailable; smooth_trend() is the one to call.
#
# The shrinking window matters more than it looks: padding the ends with raw
# values leaves 2014 and 2100 unsmoothed, and in this data both happen to be
# unrepresentative (-1.34 in a series reaching -3.99; +1.40 in one reaching
# +3.75). Those are the two rings a viewer compares first, so raw endpoints
# actively hid the warming.
running_mean <- function(x, k = 11) {
  n <- length(x)
  half <- (k - 1) %/% 2
  vapply(seq_len(n), function(i) mean(x[max(1, i - half):min(n, i + half)]),
         numeric(1))
}

# GAM smooth, matching the approach used in the esa-2026 deck.
#
# Year-to-year noise (sd of successive differences 1.65 degC) swamps the trend
# when read ring to ring. Both smoothers agree closely -- correlation 0.98-0.99,
# max difference 0.64 degC -- but the GAM is monotone rising in 100% of years
# for TaiESM1 and UKESM1 against 66% and 73% for the running mean. For a strip
# whose whole job is showing an increase, a series without local reversals is
# strictly better, so the GAM wins on more than consistency.
smooth_trend <- function(y, x = seq_along(y), k = 11) {
  if (!requireNamespace("mgcv", quietly = TRUE)) return(running_mean(y, k))
  as.numeric(mgcv::predict.gam(mgcv::gam(y ~ s(x)), data.frame(x = x)))
}

# Latewood hairline warms from near-black at the pith to a dark burnt brown at
# the rim. Atmospheric only -- it never touches the species fills, so identity
# is preserved.
latewood_colour <- function(t, sc, n = 256) {
  ramp <- colorRampPalette(c("#0d0d0d", "#301608"))(n)
  f <- (t - sc$t_lo) / (sc$t_hi - sc$t_lo)
  ramp[pmax(1, pmin(n, round(f * (n - 1)) + 1))]
}

# ---- landscape footprints (locator map) --------------------------------------

# Where landscapes 01-03 actually sit, taken from each run's env.grid.tif -- the
# same source esa-2026/R/map_plot.R uses. Cached to a local GeoPackage on first
# call so the render path never depends on the share.
#
# Projected to EPSG:3338 to match AK_no_islands.shp / AK_interior.shp, which are
# already local. The env grids are ~58 kha each, so at Alaska scale they are
# specks -- the locator draws them as marked points, not as shapes to be read.
load_footprints <- function(cache = footprint_f) {
  if (file.exists(cache)) return(terra::vect(cache))
  if (!requireNamespace("terra", quietly = TRUE)) stop("terra needed")
  message("Deriving landscape footprints from the share (one-off)...")
  ids <- sprintf("landscape_alaska_%02d", 1:3)
  polys <- lapply(ids, function(id) {
    f <- file.path(landscape_src, id, "gis", "env.grid.tif")
    if (!file.exists(f)) stop("Missing env grid: ", f)
    r <- terra::rast(f)
    p <- terra::as.polygons(!is.na(r), dissolve = TRUE)
    p <- p[as.data.frame(p)[[1]] == 1, ]           # keep the valid-data part
    p <- terra::project(p, "EPSG:3338")
    p$landscape <- id
    p[, "landscape"]
  })
  out <- do.call(rbind, polys)
  dir.create(dirname(cache), recursive = TRUE, showWarnings = FALSE)
  terra::writeVector(out, cache, overwrite = TRUE)
  out
}

# ---- ring geometry -----------------------------------------------------------

# Ring width is suppressed by BOTH a dry year and a disturbed one, as in a real
# tree. That makes a thin ring ambiguous on its own -- the scars are what
# disambiguate it: thin ring with a scar is disturbance, thin ring without is dry.
#
# Precipitation enters as a PERCENTILE RANK, not a linear min-max stretch.
#
# Annual precip is near-normal with a CV of 8.5%, so linear scaling spends most of
# the width range on rare extremes and leaves the bulk of years bunched: the
# middle 50% of years occupied just 0.236 of the available 0-1 range, which reads
# as no variation at all. Ranking against the pooled distribution doubles that
# span to 0.496 and spends the width channel on the years that actually exist.
# (Those two figures are fractions of the width range, so they hold regardless of
# how big the disc is rendered; the pixel figures they were first measured in do
# not, and have been removed rather than left to rot.)
#
# The scale is therefore relative AND non-linear -- the caption must say so.
# Dendrochronology does the same thing: ring width index is a normalised
# quantity, not millimetres.
#
# Defaults match the tuned values in iland_rings.R. They disagreed for a while,
# which is harmless only as long as every caller passes them explicitly.
ring_widths <- function(di, pr, sc, w_clim_min = 0.20, dist_pen = 0.45,
                        dist_gamma = 0.35) {
  pr_n <- sc$pr_rank(pr)
  w_clim <- w_clim_min + (1 - w_clim_min) * pr_n
  w_dist <- 1 - dist_pen * dist_amp(di, sc, dist_gamma)
  w_clim * w_dist
}

# Disturbance on the SHARED scale, compressed by dist_gamma.
#
# Two jobs in one place, because ring width and scar size must never disagree
# about how bad a year was. Previously each computed di / max(di) from the single
# run it held, so every disc rendered its own worst year at full strength and a
# grid could not be compared -- per-run maxima span 14.9x across the 108 runs (0.0137-0.2038).
#
# dist_gamma is what keeps a quiet run visible on a shared scale. At the old 0.7
# the quietest run's worst year narrowed its ring by only 8%, which reads as
# nothing; at 0.35 it narrows 18% while fidelity to true magnitude stays at 0.97
# (correlation between a run's real maximum and the narrowing it shows).
dist_amp <- function(di, sc, dist_gamma = 0.35) {
  if (is.null(sc$di_max)) {
    stop("sc has no di_max -- pass shared_scales(), not climate_scales(). ",
         "Without it di/NULL silently returns length zero.")
  }
  pmin(di / sc$di_max, 1)^dist_gamma
}

radii_from_widths <- function(w, r0 = 0.085, R = 1) {
  c(r0, r0 + cumsum(w) / sum(w) * (R - r0))   # n + 1 boundaries for n rings
}

# Which species lost the most share in each year -- where the wound landed in
# composition space. Used to place the scar wedge at a data-derived angle rather
# than an arbitrary one.
scar_species <- function(m) {
  dm <- rbind(0, diff(m))
  idx <- max.col(-dm, ties.method = "first")
  loss <- -dm[cbind(seq_len(nrow(dm)), idx)]
  data.frame(year = as.integer(rownames(m)),
             sp = SPECIES_ORDER_6[idx],
             loss = loss,
             stringsAsFactors = FALSE)
}
