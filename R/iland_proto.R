# Four rough prototypes: ensemble-scale companions to the tree-ring disc.
#
# DELIBERATELY UNPOLISHED. Three of these four get thrown away, so none of them
# gets the treatment the disc got -- no CVD audit, no pixel-level legibility
# measurement, no reading guide, no web derivative. The question each one has to
# answer is only: does the ENSEMBLE read as the subject, is precipitation
# recoverable from the printed key, and does disturbance stay visibly separate
# from precipitation.
#
# What they all share, and why:
#
#   * Replicate spread is the subject. Measured on landscape 01, the 2100
#     deciduous share runs 0.149-0.534 with fire feedback and only 0.094-0.128
#     without it -- a 7x difference in standard deviation. Per-GCM MEANS differ
#     by just 0.02 (0.338 / 0.339 / 0.319) across +2.0 to +6.2 degC of warming.
#     So replicate variability dwarfs the climate model, and the disc -- which
#     varies by model -- was showing the least important axis.
#
#   * Every run starts at EXACTLY 0.3992 deciduous share, so a single trunk or
#     origin is literal, not a drafting convenience.
#
#   * Climate is REGIONAL and per-GCM, not per-replicate. Encoding it on all 12
#     replicates would repeat one series twelve times, so it is drawn once per
#     model as its own object. That is also what decouples it from disturbance.
#
#   * Precipitation is LINEAR with a printed calibration bar, not the disc's
#     percentile rank. See pr_width() in iland_ensemble.R.
#
#   * Disturbance is a SEPARATE channel from precipitation everywhere -- ticks,
#     notches or riffles, never width. That was the third of the three critiques.

library(ragg)
library(here)

source(here("R", "iland_ensemble.R"))

pw <- 1600
ph <- 1200

# Observed 2100 range is 0.149-0.534 and every run starts at 0.3992, so this
# spans the data with a little air at both ends.
DECID_LO <- 0.05
DECID_HI <- 0.60

# Top-N most disturbed years per run get a mark. Rank-based on purpose: it needs
# no shared magnitude scale, so these prototypes can use the fast
# climate_scales() instead of shared_scales(), which iterates all 108 runs.
DIST_TOP_N <- 8

# Two disturbance colours, because the mark lands on different grounds. The
# disc's charred brown reads as burnt tissue when it sits ON a species fill, but
# on the #0d0d0d background it is invisible -- which is exactly how the canopy's
# marks disappeared on the first render. A light ember is the reverse.
DIST_ON_FILL <- "#1c0f06"
DIST_ON_BG   <- "#ffb38a"

# End minus start, NOT diff(range()). The smoothed series dips below its 2014
# value before rising, so the range overstates warming: NorESM2-MM reads +2.8
# from the range and +2.0 from the endpoints, and +2.0 is the figure quoted
# everywhere else in this repo.
warming <- function(tmean) tmean[length(tmean)] - tmean[1]

# Linear mapping helper: value range -> position range.
mk_scale <- function(lo, hi, a, b) function(v) a + (b - a) * (v - lo) / (hi - lo)

# Vectorised in BOTH arguments. adjustcolor() cannot do this: it accepts a vector
# of colours or a vector of alphas but not one colour against many alphas, which
# is exactly what fading a thread along its length needs.
alpha <- function(col, a) {
  n <- max(length(col), length(a))
  m <- grDevices::col2rgb(rep_len(col, n)) / 255
  grDevices::rgb(m[1, ], m[2, ], m[3, ], alpha = rep_len(a, n))
}

# ---- 1. the canopy -----------------------------------------------------------

# Cross-section becomes silhouette: the disc's sibling. Time is vertical, and the
# one horizontal axis is spent on the spruce <-> deciduous tipping question,
# which is what these runs are actually asking. Every replicate is a branch from
# one trunk, so the crown spread IS the ensemble uncertainty.
#
# One panel per GCM, because the no-feedback control has to sit beside its own
# model. Each panel carries a single climate ribbon: thickness = precipitation
# (linear, calibrated), colour = temperature. Two climate channels in one object,
# drawn once rather than twelve times.
proto_canopy <- function(path = here("render", "proto_canopy.png"),
                         landscape = "landscape_alaska_01_2015-2100scenario",
                         cl = load_climate(), sc = climate_scales(cl)) {
  e  <- ensemble(landscape, "onlysimfalse")
  ct <- ensemble(landscape, "onlysimtrue")

  agg_png(path, width = pw, height = ph, background = bg)
  on.exit(invisible(dev.off()), add = TRUE)
  par(bg = bg, mar = c(0, 0, 0, 0))
  plot.new(); plot.window(c(0, 1), c(0, 1))

  ys <- mk_scale(2014, 2100, 0.085, 0.845)
  panel_w <- 0.205
  gap <- 0.017

  for (i in seq_along(GCM_ORDER)) {
    g  <- GCM_ORDER[i]
    x0 <- 0.045 + (i - 1) * (panel_w + gap)
    rib_x <- x0 + 0.020                       # climate ribbon centreline
    xs <- mk_scale(DECID_LO, DECID_HI, x0 + 0.052, x0 + panel_w)

    runs <- by_model(e, g)
    yrs  <- runs[[1]]$years
    clm  <- climate_for(g, yrs, cl)

    # --- climate ribbon: thickness = precip mm, colour = temperature ---
    hw <- pr_width(clm$pr, sc) * 0.0125
    tc <- temp_colour(clm$tmean, sc)
    for (k in seq_len(length(yrs) - 1)) {
      polygon(c(rib_x - hw[k], rib_x + hw[k], rib_x + hw[k + 1], rib_x - hw[k + 1]),
              c(ys(yrs[k]), ys(yrs[k]), ys(yrs[k + 1]), ys(yrs[k + 1])),
              col = tc[k], border = NA)
    }

    # --- the no-feedback control, behind: fire spreads but kills nothing ---
    for (r in by_model(ct, g)) {
      lines(xs(r$decid), ys(r$years), col = alpha(ink_mut, 0.45), lwd = 1)
    }

    # --- the real runs: 12 branches, coloured by what dominates ---
    for (r in runs) {
      x <- xs(r$decid); y <- ys(r$years)
      n <- length(x)
      segments(x[-n], y[-n], x[-1], y[-1],
               col = alpha(SPECIES_COLOURS_DARK[r$dom[-n]], 0.8), lwd = 2.6)
      # disturbance, as its own channel: a tick at the worst years
      k <- order(-r$di)[seq_len(DIST_TOP_N)]
      k <- k[r$di[k] > 0]
      points(x[k], y[k], pch = 18, col = alpha(DIST_ON_BG, 0.9), cex = 0.8)
    }

    text(x0 + panel_w / 2 + 0.010, 0.885, g, col = GCM_COLOURS[g],
         adj = c(0.5, 0), cex = pt(20), font = 2)
    text(x0 + panel_w / 2 + 0.010, 0.862,
         sprintf("%+.1f°C", warming(clm$tmean)),
         col = ink_mut, adj = c(0.5, 0), cex = pt(15))

    # deciduousness axis, once per panel
    for (v in c(0.1, 0.3, 0.5)) {
      text(xs(v), 0.062, v, col = ink_mut, adj = c(0.5, 1), cex = pt(13))
      segments(xs(v), 0.072, xs(v), 0.078, col = ink_mut, lwd = 1)
    }
  }

  # year axis on the left
  for (yy in seq(2020, 2100, 20)) {
    text(0.036, ys(yy), yy, col = ink_mut, adj = c(1, 0.5), cex = pt(14))
  }

  # ---- chrome ----
  text(0.045, 0.975, "Every future this forest has", col = ink_pri,
       adj = c(0, 1), cex = pt(34))
  text(0.045, 0.935,
       sprintf("iLand · Landscape %s · ssp245 · 12 replicates per model, from one starting forest",
               sub("^landscape_alaska_(\\d+)_.*$", "\\1", landscape)),
       col = ink_sec, adj = c(0, 1), cex = pt(18))

  L <- 0.735
  species_key(L, 0.790, 0.036, cex = pt(16))
  temp_key(L, 0.545, 0.20, 0.020, sc, cex = pt(15),
           label = "ribbon colour — temperature trend")

  # The calibration key, two anchors rather than one bar -- see width_key().
  width_key(L, 0.430, 0.0125, sc, cex = pt(15),
            label = "ribbon width — rainfall, linear not ranked")

  segments(L, 0.345, L + 0.030, 0.345, col = alpha(ink_mut, 0.7), lwd = 1)
  text(L + 0.042, 0.345, "no fire feedback (9 runs)", col = ink_mut,
       adj = c(0, 0.5), cex = pt(15))
  points(L + 0.015, 0.310, pch = 18, col = DIST_ON_BG, cex = 1.1)
  text(L + 0.042, 0.310, "most disrupted years", col = ink_mut,
       adj = c(0, 0.5), cex = pt(15))

  text(L, 0.258,
       paste0("Left to right within a panel is the share of\n",
              "landscape dominated by deciduous trees.\n\n",
              "All 36 runs begin at exactly 0.399. Where\n",
              "fire kills trees they fan to 0.15–0.53 by\n",
              "2100; where it does not they collapse\n",
              "together to 0.09–0.13. Fire is what makes\n",
              "the future uncertain — and it matters far\n",
              "more than which climate model you pick."),
       col = ink_mut, adj = c(0, 1), cex = pt(15))
  invisible(path)
}

# ---- 2. the braided river ----------------------------------------------------

# An aerial view of a braided glacial channel -- the Tanana, the Knik. Time flows
# left to right, each replicate is a channel, and the channels merge where
# replicates agree and braid apart where they diverge.
#
# Ribbon thickness is precipitation, which is regional and so identical across
# the 12 channels of a model. That is not a compromise: it makes the whole river
# swell and pinch together, which is exactly how discharge behaves.
proto_braid <- function(path = here("render", "proto_braid.png"),
                        landscape = "landscape_alaska_01_2015-2100scenario",
                        cl = load_climate(), sc = climate_scales(cl)) {
  e <- ensemble(landscape, "onlysimfalse")

  agg_png(path, width = pw, height = ph, background = bg)
  on.exit(invisible(dev.off()), add = TRUE)
  par(bg = bg, mar = c(0, 0, 0, 0))
  plot.new(); plot.window(c(0, 1), c(0, 1))

  xs <- mk_scale(2014, 2100, 0.055, 0.985)
  row_h <- 0.225
  row_y <- c(0.635, 0.395, 0.155)      # top to bottom, one per GCM

  for (i in seq_along(GCM_ORDER)) {
    g <- GCM_ORDER[i]
    runs <- by_model(e, g)
    yrs <- runs[[1]]$years
    clm <- climate_for(g, yrs, cl)
    ys <- mk_scale(DECID_LO, DECID_HI, row_y[i], row_y[i] + row_h)

    # thickness = precip, in the row's own units
    hw <- pr_width(clm$pr, sc) * row_h * 0.075

    for (r in runs) {
      y <- ys(r$decid)
      cols <- alpha(SPECIES_COLOURS_DARK[r$dom], 0.30)
      for (k in seq_len(length(yrs) - 1)) {
        polygon(c(xs(yrs[k]), xs(yrs[k + 1]), xs(yrs[k + 1]), xs(yrs[k])),
                c(y[k] - hw[k], y[k + 1] - hw[k + 1],
                  y[k + 1] + hw[k + 1], y[k] + hw[k]),
                col = cols[k], border = NA)
      }
      # riffles: disturbance as its own mark, not as width
      k <- order(-r$di)[seq_len(DIST_TOP_N)]
      k <- k[r$di[k] > 0]
      segments(xs(yrs[k]), y[k] - hw[k] * 1.6, xs(yrs[k]), y[k] + hw[k] * 1.6,
               col = alpha("#1c0f06", 0.85), lwd = 1.6)
    }

    # Temperature strip immediately BELOW the title and ABOVE its own channels,
    # so it is unambiguously attached. Below the band it sat next to the next
    # model's label and read as belonging to that one instead.
    tc <- temp_colour(clm$tmean, sc)
    sy <- row_y[i] + row_h + 0.006
    for (k in seq_len(length(yrs) - 1)) {
      rect(xs(yrs[k]), sy, xs(yrs[k + 1]), sy + 0.011, col = tc[k], border = NA)
    }
    text(0.055, sy + 0.020, g, col = GCM_COLOURS[g],
         adj = c(0, 0), cex = pt(19), font = 2)
    text(0.985, sy + 0.020, sprintf("%+.1f°C", warming(clm$tmean)),
         col = ink_mut, adj = c(1, 0), cex = pt(15))
  }

  for (yy in seq(2020, 2100, 20)) {
    text(xs(yy), 0.140, yy, col = ink_mut, adj = c(0.5, 1), cex = pt(14))
  }

  text(0.055, 0.975, "A braided century", col = ink_pri, adj = c(0, 1),
       cex = pt(34))
  text(0.055, 0.938,
       paste0("iLand · Landscape ",
              sub("^landscape_alaska_(\\d+)_.*$", "\\1", landscape),
              " · ssp245 · 12 replicate channels per climate model.  ",
              "Height in a band is deciduous share; channels braid where the runs disagree."),
       col = ink_sec, adj = c(0, 1), cex = pt(16))

  species_key(0.055, 0.088, 0.030, cex = pt(14), cols = 3, colgap = 0.175)
  width_key(0.500, 0.045, row_h * 0.075, sc, cex = pt(14),
            label = "channel width — rainfall")
  segments(0.500, 0.006, 0.500, 0.020, col = DIST_ON_FILL, lwd = 2)
  text(0.515, 0.013, "riffle = most disrupted years", col = ink_mut,
       adj = c(0, 0.5), cex = pt(14))
  temp_key(0.800, 0.040, 0.16, 0.014, sc, cex = pt(14),
           label = "strip above each band — temperature")
  invisible(path)
}

# ---- 3. the trajectory fan ---------------------------------------------------

# Composition space rather than time: where does this forest END UP. Neither the
# disc nor a faceted time series can show that, which makes this the most
# analytically novel of the four -- and the one that reads most like a chart
# rather than an object.
#
# The six categories split cleanly three-and-three, so a ternary is honest here:
# conifer (Pima, Pigl) / deciduous (Potr, Bene) / mixed (both Mixed classes).

# Endpoints of a constant-component line across the FULL simplex, for gridlines.
# Inverting tern_xy: decid = 2y/sqrt3, mixed = x - y/sqrt3, conifer = 1 - x - y/sqrt3.
tern_gridline <- function(kind, k) {
  h <- sqrt(3) / 2
  if (kind == "decid") {
    list(x = c(k / 2, 1 - k / 2), y = c(h * k, h * k))
  } else if (kind == "conifer") {
    d <- c(0, 1 - k)
    list(x = (1 - k) - 0.5 * d, y = h * d)
  } else {
    d <- c(0, 1 - k)
    list(x = k + 0.5 * d, y = h * d)
  }
}

proto_fan <- function(path = here("render", "proto_fan.png"),
                      landscape = "landscape_alaska_01_2015-2100scenario",
                      cl = load_climate(), sc = climate_scales(cl)) {
  e  <- ensemble(landscape, "onlysimfalse")
  ct <- ensemble(landscape, "onlysimtrue")
  P  <- lapply(e,  function(r) tern_xy(tern_coords(r$m)))
  P0 <- lapply(ct, function(r) tern_xy(tern_coords(r$m)))

  # ZOOMED. Drawn full-simplex first and it failed outright: the runs occupy a
  # small pocket, so roughly 90% of the triangle was empty and the fraying -- the
  # entire point of the figure -- collapsed into an illegible knot. That is a fact
  # about the data, not a drafting slip, so the honest fix is to zoom and keep a
  # locator inset showing how little of composition space is ever visited.
  allp <- do.call(rbind, c(P, P0))
  cx <- mean(range(allp[, "x"])); cy <- mean(range(allp[, "y"]))
  half <- max(diff(range(allp[, "x"])), diff(range(allp[, "y"]))) / 2 * 1.16
  xlim <- c(cx - half, cx + half); ylim <- c(cy - half, cy + half)

  agg_png(path, width = pw, height = ph, background = bg)
  on.exit(invisible(dev.off()), add = TRUE)
  par(bg = bg)

  # Its own region with the data in TRUE ternary coordinates, so base R clips the
  # gridlines to the window for free and asp = 1 keeps the geometry unsheared.
  par(fig = c(0.025, 0.66, 0.03, 0.97), mar = c(0, 0, 0, 0), new = FALSE)
  plot.new(); plot.window(xlim, ylim, asp = 1)

  # 5% gridlines in all three families
  for (k in seq(0, 1, 0.05)) {
    for (kind in c("decid", "conifer", "mixed")) {
      g <- tern_gridline(kind, k)
      lines(g$x, g$y, col = alpha(ink_mut, 0.20), lwd = 1)
    }
  }
  for (k in seq(0, 1, 0.05)) {
    g <- tern_gridline("decid", k)
    if (g$y[1] > ylim[1] && g$y[1] < ylim[2]) {
      text(xlim[1], g$y[1], sprintf("%.0f%% deciduous", k * 100), col = ink_mut,
           adj = c(0, -0.3), cex = pt(13))
    }
  }

  for (i in seq_along(P0)) {
    lines(P0[[i]][, "x"], P0[[i]][, "y"], col = alpha(ink_mut, 0.40), lwd = 1)
  }
  for (i in seq_along(P)) {
    p <- P[[i]]; n <- nrow(p)
    a <- seq(0.10, 0.90, length.out = n - 1)
    segments(p[-n, "x"], p[-n, "y"], p[-1, "x"], p[-1, "y"],
             col = alpha(GCM_COLOURS[e[[i]]$model], a), lwd = 1.7)
    points(p[n, "x"], p[n, "y"], pch = 21, bg = GCM_COLOURS[e[[i]]$model],
           col = bg, cex = 1.6, lwd = 1)
  }
  p0 <- P[[1]][1, , drop = FALSE]
  points(p0[, "x"], p0[, "y"], pch = 21, bg = ink_pri, col = bg, cex = 2.6,
         lwd = 1.5)
  text(p0[, "x"], p0[, "y"] + half * 0.045, "2014 — every run starts here",
       col = ink_pri, adj = c(0.5, 0), cex = pt(15), font = 2)

  # ---- locator inset: the whole simplex, with the zoom marked ----
  par(fig = c(0.055, 0.185, 0.055, 0.230), mar = c(0, 0, 0, 0), new = TRUE)
  plot.new(); plot.window(c(-0.06, 1.06), c(-0.06, 1.0), asp = 1)
  h <- sqrt(3) / 2
  polygon(c(0, 1, 0.5), c(0, 0, h), col = alpha(ink_mut, 0.10),
          border = alpha(ink_mut, 0.55), lwd = 1)
  rect(xlim[1], ylim[1], xlim[2], ylim[2], border = ink_pri, lwd = 1.6)
  text(0, -0.03, "conifer", col = ink_mut, adj = c(0.5, 1), cex = pt(12))
  text(1, -0.03, "mixed", col = ink_mut, adj = c(0.5, 1), cex = pt(12))
  text(0.5, h + 0.02, "deciduous", col = ink_mut, adj = c(0.5, 0), cex = pt(12))

  # ---- chrome ----
  par(fig = c(0.66, 1, 0, 1), mar = c(0, 0, 0, 0), new = TRUE)
  plot.new(); plot.window(c(0, 1), c(0, 1))
  L <- 0.06
  text(L, 0.955, "Where the forest lands", col = ink_pri, adj = c(0, 1),
       cex = pt(30))
  text(L, 0.915,
       paste0("iLand · Landscape ",
              sub("^landscape_alaska_(\\d+)_.*$", "\\1", landscape),
              " · ssp245 · 63 runs"),
       col = ink_sec, adj = c(0, 1), cex = pt(17))

  gcm_key(L, 0.855, 0.040, cex = pt(16))
  segments(L, 0.730, L + 0.055, 0.730, col = alpha(ink_mut, 0.6), lwd = 1)
  text(L + 0.075, 0.730, "no fire feedback", col = ink_mut, adj = c(0, 0.5),
       cex = pt(16))
  points(L + 0.027, 0.688, pch = 21, bg = ink_pri, col = bg, cex = 2, lwd = 1.5)
  text(L + 0.075, 0.688, "shared 2014 start", col = ink_mut, adj = c(0, 0.5),
       cex = pt(16))

  temp_key(L, 0.590, 0.72, 0.022, sc, cex = pt(15),
           label = "temperature trend, by model")
  text(L, 0.535,
       "climate is regional — one series per model,\nshared by all 12 replicates",
       col = ink_mut, adj = c(0, 1), cex = pt(14))

  text(L, 0.455,
       paste0("Composition space, not time. Each thread is\n",
              "one run travelling 2014 → 2100, fading in as\n",
              "time passes; all 63 begin at the white dot.\n\n",
              "The fraying is the point. Runs sharing a\n",
              "climate model still finish far apart, while the\n",
              "no-feedback runs stay bundled — so the\n",
              "spread is fire, not climate.\n\n",
              "This is a ZOOM. The inset shows how small a\n",
              "pocket of all possible compositions the forest\n",
              "ever visits — it never approaches any corner.\n",
              "Gridlines are 5% apart."),
       col = ink_mut, adj = c(0, 1), cex = pt(15))
  invisible(path)
}

# ---- 4. the core rack --------------------------------------------------------

# Keeps the disc's physical-specimen conceit but changes the specimen: a rack of
# vertical cores, one per replicate, read by scanning across. The most legible of
# the four and the least novel -- it is a stacked area chart in a good coat.
#
# Core width is precipitation, so the whole rack swells and pinches in step
# (climate is regional). Disturbance is a dark notch cut across a core.
proto_rack <- function(path = here("render", "proto_rack.png"),
                       landscape = "landscape_alaska_01_2015-2100scenario",
                       cl = load_climate(), sc = climate_scales(cl)) {
  e <- ensemble(landscape, "onlysimfalse")

  agg_png(path, width = pw, height = ph, background = bg)
  on.exit(invisible(dev.off()), add = TRUE)
  par(bg = bg, mar = c(0, 0, 0, 0))
  plot.new(); plot.window(c(0, 1), c(0, 1))

  ys <- mk_scale(2014, 2100, 0.100, 0.855)
  pitch <- 0.0225        # centre-to-centre spacing of cores
  half  <- 0.0088        # half-width of a core at maximum rainfall
  grpgap <- 0.030

  yrs <- e[[1]]$years
  ny <- length(yrs)
  ytop <- ys(yrs) + (ys(yrs[2]) - ys(yrs[1]))   # each year occupies one slab

  cx0 <- 0.055
  for (i in seq_along(GCM_ORDER)) {
    g <- GCM_ORDER[i]
    runs <- by_model(e, g)
    clm <- climate_for(g, yrs, cl)
    hw <- pr_width(clm$pr, sc) * half

    gx0 <- cx0 + (i - 1) * (12 * pitch + grpgap)

    for (j in seq_along(runs)) {
      r <- runs[[j]]
      cx <- gx0 + (j - 1) * pitch

      # species bands stacked across the core width
      cum <- t(apply(r$m, 1, function(v) c(0, cumsum(v))))
      for (k in seq_along(SPECIES_ORDER_6)) {
        rect(cx - hw + 2 * hw * cum[, k], ys(yrs),
             cx - hw + 2 * hw * cum[, k + 1], ytop,
             col = SPECIES_COLOURS_DARK_6[k], border = NA)
      }
      # Disturbance notch -- separate channel from width. Cut in from the SIDES
      # rather than across the whole core: a full-width block at 8 years per core
      # punched so much of the core out that the species bands stopped reading.
      kk <- order(-r$di)[seq_len(DIST_TOP_N)]
      kk <- kk[r$di[kk] > 0]
      rect(cx - hw * 1.30, ys(yrs[kk]), cx - hw * 0.35, ytop[kk],
           col = DIST_ON_FILL, border = NA)
      rect(cx + hw * 0.35, ys(yrs[kk]), cx + hw * 1.30, ytop[kk],
           col = DIST_ON_FILL, border = NA)

      text(cx, 0.088, r$replicate, col = ink_mut, adj = c(0.5, 1), cex = pt(12))
    }

    gc <- gx0 + (12 * pitch - pitch) / 2
    text(gc, 0.875, g, col = GCM_COLOURS[g], adj = c(0.5, 0), cex = pt(19),
         font = 2)
    # temperature strip beside each group
    tc <- temp_colour(clm$tmean, sc)
    sx <- gx0 + 12 * pitch - pitch * 0.35
    for (k in seq_len(ny - 1)) {
      rect(sx, ys(yrs[k]), sx + 0.005, ys(yrs[k + 1]), col = tc[k], border = NA)
    }
  }

  for (yy in seq(2020, 2100, 20)) {
    text(0.046, ys(yy), yy, col = ink_mut, adj = c(1, 0.5), cex = pt(14))
  }

  text(0.055, 0.975, "Thirty-six cores from the same forest", col = ink_pri,
       adj = c(0, 1), cex = pt(32))
  text(0.055, 0.938,
       paste0("iLand · Landscape ",
              sub("^landscape_alaska_(\\d+)_.*$", "\\1", landscape),
              " · ssp245 · one core per replicate. Bands across a core are species share; ",
              "read variability by scanning sideways."),
       col = ink_sec, adj = c(0, 1), cex = pt(16))

  species_key(0.055, 0.055, 0.026, cex = pt(14), cols = 3, colgap = 0.150)
  width_key(0.520, 0.010, half, sc, cex = pt(14), len = 0.030,
            label = "core width — rainfall")
  temp_key(0.700, 0.028, 0.15, 0.013, sc, cex = pt(14),
           label = "strip right of each group — temperature")
  rect(0.885, 0.024, 0.899, 0.042, col = DIST_ON_FILL, border = NA)
  text(0.907, 0.033, "notch = most disrupted years", col = ink_mut,
       adj = c(0, 0.5), cex = pt(14))
  invisible(path)
}

# ---- run ---------------------------------------------------------------------

PROTOS <- list(canopy = proto_canopy, braid = proto_braid,
               fan = proto_fan, rack = proto_rack)

if (!isTRUE(getOption("iland_proto_no_run", FALSE))) {
  dir.create(here("render"), recursive = TRUE, showWarnings = FALSE)
  cl <- load_climate()
  sc <- climate_scales(cl)
  for (nm in getOption("iland_proto_which", names(PROTOS))) {
    p <- PROTOS[[nm]](cl = cl, sc = sc)
    cat("wrote", p, "\n")
  }
}
