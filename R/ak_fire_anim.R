# Animated 1980-2020 Alaska fire history: map + synchronised annual burned-area
# subplot, encoded to MP4 and WebM for web hosting.
#
# Runs with or without ak_hillshade.tif. Until the ArcticDEM tiles finish, land
# renders as a flat tone; once the hillshade exists it is picked up automatically
# and no other change is needed.

library(terra)
library(here)
library(ragg)

source(here("R", "common.R"))         # config, paths, fire loading
source(here("R", "fire_palette.R"))   # validated ramps (silent on source)

# ---- configuration -----------------------------------------------------------

# Frame is 1920 px wide across ~1360 km of Alaska, i.e. ~700 m/px, so raster work
# at anim_res (500 m, see common.R) is already finer than the video can resolve.

# Alaska in EPSG:3338 is PORTRAIT -- 1358 km wide by 1804 km tall (the panhandle
# drags the extent south-east). A 16:9 frame would letterbox it with ~40% dead
# width, so the frame is 4:3 and the leftover column carries the annotation and
# the timeline instead of being wasted. map_split is where that split falls.
width     <- 1920
height    <- 1440
map_split <- 0.56          # map occupies 0 -> map_split, panel the remainder
fps       <- 6

# Colour. `crimson` is the validated fire-age ramp, ordered OLDEST -> NEWEST.
fire_ramp <- fire_ramps$crimson
# Accent for the frame's own year. Not part of the sequential encoding -- it is a
# highlight, so it sits above the ramp's top step in lightness. It carries no
# legend key: the record ends in 2020, so any present-tense label would be wrong,
# and the huge year numeral plus the matching timeline row already identify it.
ignite <- "#FFF6E0"

bg        <- "#0d0d0d"   # background / ocean
land_flat <- "#232320"   # land tone used when no hillshade is available
# Terrain must stay SUBORDINATE to the fire. A lighter land_hi renders lovely
# relief and then competes with the data for attention -- at #57554f the
# hillshade was the loudest element in the frame. Keep the land in a narrow dark
# band; it still shows the ranges, it just stops shouting.
land_lo   <- "#161614"   # hillshade darkest (a touch above bg so coast reads)
land_hi   <- "#3b3934"   # hillshade lightest
intak_col <- "#4e4c45"   # interior outline: must clear land_hi to stay visible
ink_pri   <- "#ffffff"
ink_sec   <- "#c3c2b7"
ink_mut   <- "#898781"
hairline  <- "#2c2c2a"

label_years <- c(2004, 2005, 2015)   # direct-labelled in the subplot

frames_dir <- here("render", "frames")
out_dir    <- here("render")
dir.create(frames_dir, recursive = TRUE, showWarnings = FALSE)

# ---- data --------------------------------------------------------------------

ak    <- load_ak()
intak <- load_interior()
hf    <- load_fires(ak)      # year-filtered, valid, clipped to the landmass
cat("fire polygons", first_year, "-", last_year, ":", nrow(hf), "\n")

annual <- get_annual_burned_area(hf)
annual$cum_ha <- cumsum(annual$burned_ha)

grid <- rast(ext(ak), resolution = anim_res, crs = target_crs)
ak_mask <- rasterize(ak, grid)

have_hs <- file.exists(hillshade_file)
if (have_hs) {
  message("using hillshade: ", hillshade_file)
  hs <- resample(rast(hillshade_file), grid, method = "average")
  hs <- mask(hs, ak_mask)

  # shade() returns the cosine of the incidence angle, which goes negative on
  # slopes facing away from the light -- more so with hillshade_z exaggeration.
  # Below zero is all equally "in shadow", so clamp: left unclamped, a linear
  # colour ramp spends a third of its range separating shades of black.
  hs <- clamp(hs, 0, 1, values = TRUE)

  # Stretch on percentiles rather than min/max. Flat ground sits at sin(35 deg)
  # = 0.57, so the informative band is narrow and the tails would otherwise eat
  # most of the ramp.
  smp <- spatSample(hs, 2e5, na.rm = TRUE)[, 1]
  hs_range <- unname(quantile(smp, c(0.02, 0.98)))
  message(sprintf("hillshade stretch: %.3f -> %.3f", hs_range[1], hs_range[2]))
} else {
  message("no hillshade yet (", hillshade_file, ") -- rendering flat land")
  hs <- NULL
  hs_range <- NULL
}

years <- first_year:last_year
n <- length(years)

# ---- per-year burn masks -----------------------------------------------------
# Rasterised once up front so the frame loop is pure compositing.

burn_by_year <- vector("list", n)
for (i in seq_len(n)) {
  y <- years[i]
  s <- hf[hf$FIREYEAR == y, ]
  burn_by_year[[i]] <- if (nrow(s) == 0) NULL else rasterize(s, grid)
  cat(sprintf("\rrasterising %d (%d/%d)", y, i, n))
}
cat("\n")

# ---- frame renderer ----------------------------------------------------------

fmt_ha <- function(x) formatC(round(x), big.mark = ",", format = "d")

# Fixed colour lookup for age 0..40, so a given fire age is the same colour in
# every frame. Deriving it per frame would let the ramp rescale as the age range
# grows and the animation would drift.
#
# The ramp's palest step is held back from the age scale: it is nearly as light as
# the ignition accent, so including it made this year's fires indistinguishable
# from last year's and the "flash" did not read. Ages 1..40 therefore run gold ->
# crimson, leaving near-white to mean "burning now" and nothing else.
age_cols <- colorRampPalette(rev(fire_ramp[-length(fire_ramp)]))(n)

# Dissolve the interior polygon: it has 23 parts, and drawing every internal
# boundary put a grey web across the map that read as noise rather than context.
intak_outline <- aggregate(intak)

# Two regions via par(fig): terra's plot() manages its own aspect-preserving
# region, so annotation CANNOT be positioned in map user-coordinates. The panel
# gets its own plot.window(0..1) and everything in it is placed in plain 0-1
# fractions, which is both reliable and readable.
render_frame <- function(i, path) {
  y <- years[i]

  agg_png(path, width = width, height = height, background = bg)
  on.exit(invisible(dev.off()), add = TRUE)
  par(bg = bg)

  # ---------------- map region ----------------
  par(fig = c(0, map_split, 0, 1), mar = c(0, 0, 0, 0), new = FALSE)
  plot(ak, col = land_flat, border = NA, axes = FALSE, legend = FALSE,
       mar = c(0, 0, 0, 0))

  if (!is.null(hs)) {
    plot(hs, col = colorRampPalette(c(land_lo, land_hi))(128),
         add = TRUE, legend = FALSE, maxcell = Inf, range = hs_range)
  }

  # interior AK: an optional highlight, not a hard boundary
  plot(intak_outline, border = intak_col, lwd = 1.2, add = TRUE)

  # fire, coloured by age; `recent` carries the most recent burn year per pixel
  if (!is.null(recent_cache$r)) {
    age <- y - recent_cache$r
    plot(age, col = age_cols, add = TRUE, legend = FALSE, maxcell = Inf,
         range = c(0, n - 1))
  }
  # this year's fires on top, as the ignition accent
  if (!is.null(burn_by_year[[i]])) {
    plot(burn_by_year[[i]], col = ignite, add = TRUE, legend = FALSE,
         maxcell = Inf)
  }

  # ---------------- annotation panel ----------------
  par(fig = c(map_split, 1, 0, 1), mar = c(0, 0, 0, 0), new = TRUE)
  plot.new(); plot.window(c(0, 1), c(0, 1))

  L <- 0.10          # panel left margin
  R <- 0.94          # panel right margin

  text(L, 0.955, y, col = ignite, adj = c(0, 1), cex = 8.6)
  text(L + 0.005, 0.855, "Alaska wildfire, 1980-2020", col = ink_sec,
       adj = c(0, 1), cex = 2.0)

  text(L + 0.005, 0.800,
       paste0(fmt_ha(annual$cum_ha[annual$year == y]), " ha"),
       col = ink_pri, adj = c(0, 1), cex = 2.6)
  text(L + 0.005, 0.762, paste0("cumulative burned area since ", first_year),
       col = ink_mut, adj = c(0, 1), cex = 1.4)

  this_yr <- annual$burned_ha[annual$year == y]
  text(L + 0.005, 0.712, paste0(fmt_ha(this_yr), " ha this year"),
       col = ink_sec, adj = c(0, 1), cex = 1.6)

  # ---- fire-age scale legend (obliged: multi-hue semantic-heat ramp) ----
  # The gradient is built from age_cols itself, not from fire_ramp, so the legend
  # cannot drift out of step with the map the way it does if the two are written
  # independently. Reversed here because the bar reads oldest -> newest.
  ly1 <- 0.640; ly0 <- 0.615
  grad <- rev(age_cols)
  xs <- seq(L, L + 0.60, length.out = length(grad) + 1)
  for (j in seq_along(grad)) {
    rect(xs[j], ly0, xs[j + 1], ly1, col = grad[j], border = NA)
  }
  text(L, ly1 + 0.012, "years since fire", col = ink_mut, adj = c(0, 0), cex = 1.35)
  text(L, ly0 - 0.008, "40", col = ink_mut, adj = c(0, 1), cex = 1.2)
  text(L + 0.60, ly0 - 0.008, "1", col = ink_mut, adj = c(1, 1), cex = 1.2)

  # ---- annual burned area, as a vertical timeline ----
  # One series, so no legend box -- the heading names it. Horizontal bars down a
  # tall panel beat 41 vertical bars crammed into a strip, and reading downward
  # matches time moving forward. Linear scale on purpose: 2004 is 13.5x the
  # median year and that skew IS the story; a sqrt axis would flatten it.
  # Cumulative total is the number above, never a second axis.
  text(L, 0.545, "Annual area burned", col = ink_sec, adj = c(0, 0), cex = 1.6)

  ty1 <- 0.520; ty0 <- 0.055
  row <- (ty1 - ty0) / n
  bh  <- row * 0.66                 # gap between rows is surface, not a border
  bx0 <- L + 0.10
  # Hold back 22% of the track so the direct label on the longest bar has room.
  # At full width the "2.7M ha" label on 2004 ran off the edge of the frame.
  bmax <- (R - bx0) * 0.78
  ymax <- max(annual$burned_ha)

  # axis hairline, solid and one shade off the surface
  segments(bx0, ty0, bx0, ty1, col = hairline, lwd = 2)

  for (j in seq_len(n)) {
    yy <- ty1 - (j - 0.5) * row
    if (years[j] > y) {
      # A fixed-width stub, not the real value: it establishes that the timeline
      # continues without revealing data the frame has not reached yet.
      rect(bx0, yy - bh / 2, bx0 + 0.012, yy + bh / 2, col = "#211f1d",
           border = NA)
      next
    }
    len <- (annual$burned_ha[j] / ymax) * bmax
    age <- y - years[j]
    col <- if (age == 0) ignite else age_cols[min(age + 1, n)]
    rect(bx0, yy - bh / 2, bx0 + len, yy + bh / 2, col = col, border = NA)
  }

  # decade labels only; a label per row would be unreadable
  for (d in seq(1980, 2020, by = 10)) {
    j <- which(years == d)
    text(bx0 - 0.02, ty1 - (j - 0.5) * row, d, col = ink_mut,
         adj = c(1, 0.5), cex = 1.2)
  }
  # the current year always names itself
  text(bx0 - 0.02, ty1 - (which(years == y) - 0.5) * row, y, col = ignite,
       adj = c(1, 0.5), cex = 1.25)

  # selective direct labels: only the years that matter, once reached
  for (lyr in label_years) {
    if (lyr > y) next
    j <- which(years == lyr)
    len <- (annual$burned_ha[j] / ymax) * bmax
    text(bx0 + len + 0.015, ty1 - (j - 0.5) * row,
         paste0(format(round(annual$burned_ha[j] / 1e6, 1), nsmall = 1), "M ha"),
         col = ink_sec, adj = c(0, 0.5), cex = 1.2)
  }

  # Credit only what is actually in the frame: until the hillshade exists the
  # terrain source is not being used and must not be claimed.
  text(L, 0.022,
       paste0("Perimeters: AICC",
              if (!is.null(hs)) " | Terrain: ArcticDEM v4.1" else "",
              " | EPSG:3338"),
       col = ink_mut, adj = c(0, 0), cex = 1.15)
}

# ---- render ------------------------------------------------------------------

# `recent` accumulates forward in time: one raster carried across frames rather
# than 41 stacked layers, which keeps memory flat.
recent_cache <- new.env()
recent_cache$r <- NULL

render_years <- getOption("anim_years", years)   # subset for quick previews
paths <- character(0)

for (i in seq_len(n)) {
  y <- years[i]
  if (!is.null(burn_by_year[[i]])) {
    recent_cache$r <- if (is.null(recent_cache$r)) {
      ifel(is.na(burn_by_year[[i]]), NA, y)
    } else {
      ifel(!is.na(burn_by_year[[i]]), y, recent_cache$r)
    }
  }
  if (!(y %in% render_years)) next
  p <- file.path(frames_dir, sprintf("frame_%03d.png", i))
  render_frame(i, p)
  paths <- c(paths, p)
  cat(sprintf("\rrendered %d (%d/%d)", y, i, n))
}
cat("\n")

if (isTRUE(getOption("anim_encode", TRUE))) {
  library(av)
  mp4  <- file.path(out_dir, "ak_fire_1980_2020.mp4")
  webm <- file.path(out_dir, "ak_fire_1980_2020.webm")
  av_encode_video(paths, mp4, framerate = fps, codec = "libx264", verbose = FALSE)
  av_encode_video(paths, webm, framerate = fps, codec = "libvpx-vp9", verbose = FALSE)
  cat("mp4 :", mp4, round(file.size(mp4) / 1024^2, 2), "MB\n")
  cat("webm:", webm, round(file.size(webm) / 1024^2, 2), "MB\n")
}
