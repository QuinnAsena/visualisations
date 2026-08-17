# Tree-ring composition disc: iLand species composition 2014-2100 as a cross-section.
#
# Radius     = time (NOT linear -- ring width carries data, so year ticks are
#              mandatory rather than decorative)
# Angle      = share of landscape area by dominant species, in SPECIES_ORDER_6
# Ring width = precipitation percentile, narrowed again in disturbed years.
#              Both suppress growth in a real tree, so a thin ring is ambiguous
#              on its own -- the scars are what disambiguate it.
# Ring line  = a hairline at every year boundary, warmed by that year's
#              temperature. This is what makes the thing read as a core rather
#              than a wobbly pie chart: in a real core it is the latewood band
#              you count, and tight line spacing IS the signal.
# Scar       = radial wound at the sector of whichever species lost the most share
#              that year, with later rings deflected around it
# Gutter     = a wedge cut at 12 o'clock carrying the year axis and a per-year
#              temperature strip. Doubles as the slot an increment borer leaves.

library(ragg)
library(here)

source(here("R", "iland_common.R"))

# ---- configuration -----------------------------------------------------------

width      <- 1920
height     <- 1440
disc_split <- 0.64       # a bigger disc buys pixels per ring: at 87 rings the
                         # mean ring is only ~6 px, so every one counts, and
                         # those pixels are what let w_clim_min go lower

# Ring width = precipitation PERCENTILE (see ring_widths), then narrowed again in
# disturbed years. Both suppress growth in a real tree, which makes a thin ring
# ambiguous on its own -- the scars are what disambiguate it.
# Lowering this is the right lever for exaggerating precipitation. A multiplier
# on the rank looked equivalent but clipped 9-19 years to identical widths for
# the same pixel floor, which is information loss for no extra spread.
w_clim_min <- 0.20       # driest year's share of the wettest year's width.
                         # 4.9x ring ratio; the thinnest ring is 2.18 px, and
                         # below ~2 px rings dissolve under antialiasing, so
                         # this is close to the floor for a 602 px disc.
dist_pen   <- 0.45       # width lost at maximum disturbance
# Compresses the SHARED disturbance scale so a quiet run stays visible. At 0.7 the
# quietest of the 108 runs narrowed its worst ring by only 8%; at 0.35 it narrows
# 18% and fidelity to true magnitude is still 0.97. Feeds both ring width and
# scar size via dist_amp().
dist_gamma <- 0.35
r0         <- 0.085      # pith

# Pima and Pigl are adjacent in stack order and only 14.3 dE apart -- below the
# 15 floor -- so they must not abut. A surface gap is the prescribed remedy;
# re-stepping a colour esa-2026 already validated is not.
sector_gap <- 0.0025     # radians inset each side of a species sector

gutter <- 7 * pi / 180   # angular wedge at 12 o'clock: the climate strip + axis

ring_line_a   <- 0.32    # alpha of the year boundary hairline
ring_line_lwd <- 0.7

# Wounds displace boundaries by a multiple of the MEAN ring width, not by a
# fraction of absolute radius -- scaling to radius made the outer rings balloon
# into lobes because the same fraction is a far bigger distance out there.
scar_top_n      <- 8
scar_sigma      <- 0.13  # angular spread (radians)
scar_depth_rings <- 9.0  # displacement in units of mean ring width. Deep enough
                         # that fires visibly dent the outline -- the disc is not
                         # meant to stay circular. The cummax guard in
                         # boundary_field() is load-bearing at this depth.
                         #
                         # Raised from 6.0 when disturbance moved to the shared
                         # scale: at 6.0 every run reached amp 1.0 by normalising
                         # to its own worst year, so scars flattened everywhere
                         # once amp became the run's real standing in the pool.
                         #
                         # The OUTLINE IS TEXTURE, NOT A MEASUREMENT. Rendered rim
                         # dent correlates only 0.28 with a run's peak disturbance
                         # and 0.37 with its cumulative total, and narrowing
                         # scar_sigma does not improve it (0.19-0.28 at every value
                         # tried). The cause is timing: at scar_decay 0.99 a
                         # mid-century wound accumulates persistence over 40-odd
                         # rings while a 2090s wound of the same size gets a few.
                         # Authentic for a core; useless as a readout. The
                         # quantitative disturbance channel is ring WIDTH, which
                         # holds fidelity 0.976 -- and the caption claims nothing
                         # about the silhouette.
scar_decay      <- 0.99  # per-ring persistence of an older wound. Near 1 on
                         # purpose: a real catface never fully closes, and at
                         # 0.90 a mid-century scar retained 0.4% of its
                         # amplitude by the rim, so the outline stayed perfectly
                         # circular no matter how deep the wound was set.

# Deflection alone is not enough: a bent ring boundary is indistinguishable from
# a birch sector, because both read as radial features. A real fire scar shows as
# a dark lesion the tree grows over, so draw one -- charred brown rather than
# black, so it reads as burnt tissue and not as a hole punched in the disc.
#
# Width is a constant PHYSICAL half-width converted to an angle at each radius,
# not a constant angle. An angular width tapering outward produced little dark
# chevrons that read as floating debris; a constant-width radial crack reads as a
# wound. Length scales with severity, so the worst years cut furthest.
lesion_col    <- "#1c0f06"
lesion_alpha  <- 0.94
lesion_halfw  <- 0.0090  # half-width in disc-radius units
lesion_min    <- 4       # rings a scar runs outward, at minimum severity
lesion_max    <- 17      # ... and at maximum

# Locator inset. The land tone sits between the background and the interior
# outline so the silhouette reads without competing with the disc. Named rather
# than borrowing `ignite`, which belongs to ak_fire_anim.R -- a different piece
# in the same repo, and not in scope here.
# Both lifted for legibility: at #232320 the silhouette was only 1.2:1 on the
# background and the outline 2.3:1, which is present-but-unreadable. Land stays
# darker (L 0.360) than every species colour (L 0.43-0.89), so a brighter inset
# still cannot compete with the disc for attention.
locator_land <- "#3f3d37"   # 1.8:1 on bg (was #232320, 1.2:1)
intak_col    <- "#807d71"   # 4.7:1 on bg (was #4e4c45, 2.3:1); dE 22.9 from land
locator_mark <- "#FFF6E0"

n_theta <- 1441          # deflection grid resolution

sz_title <- pt(30)
sz_body  <- pt(21)
sz_label <- pt(18)
sz_axis  <- pt(16)
sz_year  <- pt(18)       # year markers: up 2 pt, bold, and on a plate

# ---- deflected boundary field ------------------------------------------------

# Radius of every ring boundary at every angle, after wounds.
#
# Computed as a full (n+1) x n_theta field and then forced monotone with cummax
# down the boundary axis. Without that, a deep wound on a thin ring can push one
# boundary past the next and the polygon self-intersects -- which shows up as a
# torn disc rather than an obvious error.
boundary_field <- function(r, scars, theta) {
  n1 <- length(r)
  mean_w <- mean(diff(r))
  fld <- matrix(r, nrow = n1, ncol = length(theta))

  if (nrow(scars) > 0) {
    for (i in seq_len(nrow(scars))) {
      dth <- (theta - scars$angle[i] + pi) %% (2 * pi) - pi
      shape <- exp(-(dth / scar_sigma)^2)
      k <- scars$ring[i]
      rows <- k:n1
      decay <- scar_decay^(rows - k)
      fld[rows, ] <- fld[rows, ] -
        outer(scars$amp[i] * decay * scar_depth_rings * mean_w, shape)
    }
  }

  fld <- pmax(fld, r0 * 0.5)
  fld <- apply(fld, 2, cummax)          # never let a boundary overtake the next
  fld
}

# Linear interpolation of one boundary's radius at arbitrary angles.
rad_at <- function(fld, k, theta, theta_grid) {
  th <- theta %% (2 * pi)
  approx(theta_grid, fld[k, ], xout = th, rule = 2)$y
}

# ---- drawing -----------------------------------------------------------------

draw_disc <- function(m, di, pr, tmean, sc, cx = 0.5, cy = 0.5, scale = 0.49) {
  n <- nrow(m)
  years <- as.integer(rownames(m))
  w <- ring_widths(di, pr, sc, w_clim_min = w_clim_min,
                   dist_pen = dist_pen, dist_gamma = dist_gamma)
  r <- radii_from_widths(w, r0 = r0)

  # Scars sit at the angular midpoint of the sector that lost the most share
  # that year -- the wound lands where the damage did.
  # NOT `sc` -- that is the climate-scales argument. Shadowing it made
  # sc$t_lo NULL, and arithmetic against NULL silently returns length zero
  # rather than erroring where the mistake was.
  scsp <- scar_species(m)
  keep <- order(-di)[seq_len(min(scar_top_n, n))]
  keep <- keep[di[keep] > 0]
  span <- 2 * pi - gutter
  a_start <- pi / 2 - gutter / 2

  scars <- data.frame(ring = integer(0), angle = numeric(0), amp = numeric(0))
  if (length(keep) > 0) {
    ang <- vapply(keep, function(k) {
      cum <- c(0, cumsum(m[k, ]))
      j <- match(scsp$sp[k], SPECIES_ORDER_6)
      a_start - span * (cum[j] + cum[j + 1]) / 2
    }, numeric(1))
    # Same shared, gamma-compressed amplitude the ring widths use, so scar size
    # and ring narrowing can never tell different stories about the same year.
    scars <- data.frame(ring = keep, angle = ang %% (2 * pi),
                        amp = dist_amp(di[keep], sc, dist_gamma))
  }

  theta_grid <- seq(0, 2 * pi, length.out = n_theta)
  fld <- boundary_field(r, scars, theta_grid)

  # Latewood warms with the year's temperature. Atmospheric only -- it colours
  # the boundary hairline, never the species fills, so identity is preserved.
  line_cols <- grDevices::adjustcolor(latewood_colour(tmean, sc),
                                      alpha.f = ring_line_a)

  for (k in seq_len(n)) {
    cum <- c(0, cumsum(m[k, ]))
    for (j in seq_along(SPECIES_ORDER_6)) {
      a0 <- a_start - span * cum[j]
      a1 <- a_start - span * cum[j + 1]
      if (a0 - a1 <= 2 * sector_gap) next
      th <- seq(a0 - sector_gap, a1 + sector_gap,
                length.out = max(3, ceiling((a0 - a1) / 0.012)))
      ri <- rad_at(fld, k,     th, theta_grid)
      # a hair of outward overlap covers the antialiasing seam between rings
      ro <- rad_at(fld, k + 1, th, theta_grid) * 1.0015
      polygon(c(cx + ri * cos(th) * scale, cx + rev(ro * cos(th)) * scale),
              c(cy + ri * sin(th) * scale, cy + rev(ro * sin(th)) * scale),
              col = SPECIES_COLOURS_DARK_6[j], border = NA)
    }
  }

  # Lesions: a narrow radial crack from the injury year outward, closing at the
  # tip as the callus grows over it.
  if (nrow(scars) > 0) {
    lcol <- grDevices::adjustcolor(lesion_col, alpha.f = lesion_alpha)
    for (i in seq_len(nrow(scars))) {
      k0 <- scars$ring[i]
      len <- round(lesion_min + (lesion_max - lesion_min) * scars$amp[i])
      k1 <- min(k0 + len, n + 1)
      if (k1 <= k0) next
      ks <- k0:k1
      # radius along the crack's centreline, and the angular half-width that
      # keeps its physical width constant
      rc <- vapply(ks, function(kk) rad_at(fld, kk, scars$angle[i], theta_grid),
                   numeric(1))
      close <- 1 - ((ks - k0) / (k1 - k0))^2.2      # tapers shut at the tip
      hw <- pmin(lesion_halfw * close / pmax(rc, 1e-3), 0.09)
      left  <- scars$angle[i] - hw
      right <- scars$angle[i] + hw
      polygon(c(cx + rc * cos(left) * scale, rev(cx + rc * cos(right) * scale)),
              c(cy + rc * sin(left) * scale, rev(cy + rc * sin(right) * scale)),
              col = lcol, border = NA)
    }
  }

  # Year boundaries last, over the fills: this is the ring structure itself.
  th <- seq(a_start - span, a_start, length.out = 900)
  for (k in 2:(n + 1)) {
    rr <- rad_at(fld, k, th, theta_grid)
    lines(cx + rr * cos(th) * scale, cy + rr * sin(th) * scale,
          col = line_cols[min(k, n)], lwd = ring_line_lwd)
  }

  # Pith
  tt <- seq(0, 2 * pi, length.out = 200)
  polygon(cx + r0 * 0.92 * cos(tt) * scale, cy + r0 * 0.92 * sin(tt) * scale,
          col = bg, border = NA)

  # ---- the gutter: a per-year temperature strip, doubling as the year axis ----
  gth <- seq(a_start, a_start + gutter, length.out = 24)
  outer_r <- max(fld[n + 1, ]) * 1.002
  polygon(c(cx + r0 * 0.9 * cos(gth) * scale, cx + rev(outer_r * cos(gth)) * scale),
          c(cy + r0 * 0.9 * sin(gth) * scale, cy + rev(outer_r * sin(gth)) * scale),
          col = bg, border = NA)

  # Strip and labels follow the DEFLECTED boundaries, not the plain radii.
  # scar_sigma is 7.4 deg and the gutter is 7 deg wide, so any scar within ~15
  # deg of 12 o'clock displaces the rings beside the slot while the strip stayed
  # put -- a visible seam between the strip and the ring it is supposed to mark.
  gmid <- a_start + gutter / 2
  rg <- vapply(seq_len(n + 1), function(k) rad_at(fld, k, gmid, theta_grid),
               numeric(1))

  tcols <- temp_colour(tmean, sc)
  for (k in seq_len(n)) {
    polygon(c(cx + rg[k] * cos(gth) * scale, cx + rev(rg[k + 1] * 1.02 * cos(gth)) * scale),
            c(cy + rg[k] * sin(gth) * scale, cy + rev(rg[k + 1] * 1.02 * sin(gth)) * scale),
            col = tcols[k], border = NA)
  }

  # Year labels sit on a background plate so they stay legible over whatever the
  # strip or the fills happen to be at that radius -- previously they were muted
  # grey over an orange wedge and effectively invisible.
  amid <- a_start + gutter * 1.35
  lab_years <- years[years %% 20 == 0]
  for (yy in c(years[1], lab_years)) {
    k <- match(yy, years)
    rr <- rg[k]
    x <- cx + rr * cos(amid) * scale
    y <- cy + rr * sin(amid) * scale
    lab <- as.character(yy)
    wpad <- strwidth(lab, cex = sz_year) / 2 + 0.006
    hpad <- strheight(lab, cex = sz_year) / 2 + 0.005
    rect(x - 2 * wpad, y - hpad, x + 0.004, y + hpad, col = bg, border = NA)
    text(x - 0.008, y, lab, col = ink_sec, adj = c(1, 0.5), cex = sz_year,
         font = 2)
    segments(x - 0.004, y, x + 0.010, y, col = ink_sec, lwd = 1.6)
  }

  invisible(list(r = r, w = w, scars = scars, fld = fld))
}

# ---- locator map -------------------------------------------------------------

# Draw a SpatVector's rings with base polygon(). terra's own plot() would need
# its own device window, which cannot be nested inside the disc's coordinates.
draw_rings <- function(v, tx, ty, col, border = NA, lwd = 1) {
  g <- terra::geom(v)
  key <- paste(g[, "geom"], g[, "part"], g[, "hole"], sep = "_")
  for (k in unique(key)) {
    i <- key == k
    polygon(tx(g[i, "x"]), ty(g[i, "y"]), col = col, border = border, lwd = lwd)
  }
}

# Alaska inset marking where the three iLand landscapes sit, with the one (or
# ones) on show highlighted. The single-disc frame puts it in the disc region's
# free bottom-left corner; the triptych puts it under the species key. Either
# way it is aspect-corrected for whichever region it lands in -- see below.
#
# The landscapes are ~60 kha each; Alaska is ~150 Mha. At this size that is under
# 4 px across, so they are drawn as MARKERS rather than shapes -- a footprint
# polygon here would be a sub-pixel smudge pretending to be geography.
draw_locator <- function(ak, intak, fp, current, x0, y0, w, h) {
  e <- as.vector(terra::ext(ak))
  ew <- e[["xmax"]] - e[["xmin"]]
  eh <- e[["ymax"]] - e[["ymin"]]

  # Separate x and y scales, derived from the CURRENT region's physical aspect.
  # A single shared scale is only correct where one user unit is the same size in
  # both directions -- true in the disc regions (asp = 1), false in the legend
  # region, which is 810 x 1140 px over a 0-1 window. Using one scale there
  # stretched Alaska vertically.
  usr <- par("usr"); pin <- par("pin")
  ppx <- pin[1] / (usr[2] - usr[1])     # inches per user unit, x
  ppy <- pin[2] / (usr[4] - usr[3])     # inches per user unit, y
  s_phys <- min(w * ppx / ew, h * ppy / eh)
  sx <- s_phys / ppx
  sy <- s_phys / ppy
  ox <- x0 + (w - ew * sx) / 2
  oy <- y0 + (h - eh * sy) / 2
  tx <- function(x) ox + (x - e[["xmin"]]) * sx
  ty <- function(y) oy + (y - e[["ymin"]]) * sy

  # Simplify to the resolution actually rendered: at ~190 px across 1358 km,
  # one pixel is ~7 km, so sub-3 km coastline detail is invisible work.
  draw_rings(terra::simplifyGeom(ak, tolerance = 3000), tx, ty,
             col = locator_land)
  draw_rings(terra::simplifyGeom(terra::aggregate(intak), tolerance = 3000),
             tx, ty, col = NA, border = intak_col, lwd = 1.1)

  # `current` may name one landscape or several. Naming several is how a
  # landscape comparison stays on ONE large map instead of three small ones:
  # each marker is numbered, and the discs are already titled to match.
  ct <- terra::crds(terra::centroids(fp))
  ids <- fp$landscape
  for (i in seq_len(nrow(ct))) {
    on <- ids[i] %in% current
    mark <- if (on) locator_mark else ink_mut
    points(tx(ct[i, 1]), ty(ct[i, 2]), pch = 21, bg = mark, col = mark,
           cex = if (on) 1.5 else 0.9, lwd = 1)
    if (on) {
      points(tx(ct[i, 1]), ty(ct[i, 2]), pch = 1, col = locator_mark, cex = 3.2,
             lwd = 1.2)
      text(tx(ct[i, 1]) + 0.028, ty(ct[i, 2]),
           sub("^landscape_alaska_", "", ids[i]),
           col = locator_mark, adj = c(0, 0.5), cex = sz_axis, font = 2)
    }
  }
  text(x0, y0 - 0.012,
       if (length(current) > 1) "landscape locations" else "landscape location",
       col = ink_mut, adj = c(0, 1), cex = sz_axis)
}

# ---- full frame --------------------------------------------------------------

# Everything one disc needs, in one place, so the single-disc and triptych
# renderers cannot drift apart on how the climate join is done.
disc_data <- function(landscape, model, replicate, cl) {
  d <- load_area_dom(landscape)
  tr <- unique(d$treatment[d$model == model])
  if (length(tr) != 1) stop("Expected one treatment for model ", model)
  m <- share_matrix(d, tr, replicate)
  years <- as.integer(rownames(m))

  # canon_gcm bridges NorEsm2-MM (iLand) and NorESM2-MM (climate). Without it
  # this join silently returns nothing.
  cg <- cl[cl$gcm == canon_gcm(model), ]
  idx <- match(years, cg$year)
  if (anyNA(idx)) stop("Climate missing for years: ",
                       paste(years[is.na(idx)], collapse = ", "))
  list(m = m, di = disturbance_index(m), years = years,
       pr = cg$pr[idx],
       # Smoothed for colour so the trend reads instead of the annual noise; the
       # caption says so. Precipitation stays raw -- its year-to-year jitter is
       # the texture we want in the ring widths.
       tmean = smooth_trend(cg$tmean[idx], years))
}

render_disc_frame <- function(path,
                              landscape = "landscape_alaska_01_2015-2100scenario",
                              model = "NorEsm2-MM", replicate = 1,
                              cl = load_climate(), sc = shared_scales(cl)) {
  dd <- disc_data(landscape, model, replicate, cl)
  m <- dd$m; di <- dd$di; years <- dd$years; pr <- dd$pr; tmean <- dd$tmean

  agg_png(path, width = width, height = height, background = bg)
  on.exit(invisible(dev.off()), add = TRUE)
  par(bg = bg)

  # asp = 1 makes one x unit equal one y unit in device space, so the disc is a
  # circle rather than an ellipse stretched by the region's aspect.
  par(fig = c(0, disc_split, 0, 1), mar = c(0, 0, 0, 0), new = FALSE)
  plot.new(); plot.window(c(0, 1), c(0, 1), asp = 1)
  draw_disc(m, di, pr, tmean, sc)

  # asp = 1 makes y extend past [0, 1] (the region is taller than wide), so the
  # corner below the disc is real drawable space rather than off-canvas.
  draw_locator(load_ak(), load_interior(), load_footprints(),
               current = sub("_\\d{4}-\\d{4}scenario$", "", landscape),
               x0 = 0.015, y0 = -0.060, w = 0.155, h = 0.205)

  par(fig = c(disc_split, 1, 0, 1), mar = c(0, 0, 0, 0), new = TRUE)
  plot.new(); plot.window(c(0, 1), c(0, 1))
  L <- 0.06

  # Title on one line: measured at 0.641 of the 0.91 usable panel width.
  text(L, 0.945, paste(length(years), "years of a boreal landscape"),
       col = ink_pri, adj = c(0, 1), cex = sz_title)
  text(L, 0.870,
       paste0("iLand · Landscape ",
              sub("^landscape_alaska_(\\d+)_.*$", "\\1", landscape),
              " · ", min(years), "-", max(years)),
       col = ink_sec, adj = c(0, 1), cex = sz_body)
  text(L, 0.827, paste0(model, " · ssp245 · replicate ", replicate),
       col = ink_mut, adj = c(0, 1), cex = sz_body)

  # Species key in two columns: the widest label needs 0.277 of a 0.455 column,
  # so it fits comfortably and buys three rows of vertical space.
  colw <- 0.455
  for (j in seq_along(SPECIES_ORDER_6)) {
    col_i <- (j - 1) %/% 3
    row_i <- (j - 1) %% 3
    xx <- L + col_i * colw
    yy <- 0.755 - row_i * 0.048
    rect(xx, yy - 0.014, xx + 0.048, yy + 0.014,
         col = SPECIES_COLOURS_DARK_6[j], border = NA)
    text(xx + 0.066, yy, SPECIES_LABELS[SPECIES_ORDER_6[j]], col = ink_sec,
         adj = c(0, 0.5), cex = sz_label)
  }

  # ---- temperature key (the gutter strip) ----
  ty <- 0.555
  # Short label: the number is right-aligned on the same line, and the fuller
  # "GAM trend" wording lives in the caption below.
  text(L, ty + 0.030, "Temperature trend", col = ink_sec,
       adj = c(0, 0), cex = sz_label)
  # Magnitude in numbers, on the key's own line so it cannot collide with the
  # severity block below. On a scale shared across GCMs a model that warms half
  # as much SHOULD look half as dramatic -- that is the scale working, not
  # failing -- but the reader still deserves the figure.
  text(L + 0.40, ty + 0.032, sprintf("%+.1f°C", tail(tmean, 1) - tmean[1]),
       col = ink_pri, adj = c(1, 0), cex = sz_axis)
  tseq <- seq(sc$t_lo, sc$t_hi, length.out = 160)
  xs <- seq(L, L + 0.40, length.out = 161)
  for (j in seq_along(tseq)) {
    rect(xs[j], ty - 0.012, xs[j + 1], ty + 0.012,
         col = temp_colour(tseq[j], sc), border = NA)
  }
  # Freezing is a REFERENCE tick, not a pivot -- the ramp is sequential across
  # the whole range. Guarded because a shared range that never crossed zero would
  # otherwise draw this tick outside the bar.
  if (sc$t_lo < 0 && sc$t_hi > 0) {
    xz <- L + 0.40 * (0 - sc$t_lo) / (sc$t_hi - sc$t_lo)
    segments(xz, ty - 0.019, xz, ty + 0.019, col = ink_pri, lwd = 2)
    text(xz, ty - 0.026, "0", col = ink_sec, adj = c(0.5, 1), cex = sz_axis)
  }
  text(L, ty - 0.026, paste0(round(sc$t_lo, 1), "°C"), col = ink_mut,
       adj = c(0, 1), cex = sz_axis)
  text(L + 0.40, ty - 0.026, paste0("+", round(sc$t_hi, 1), "°C"),
       col = ink_mut, adj = c(1, 1), cex = sz_axis)

  # ---- severity key ----
  # Scar LENGTH and depth scale with disturbance; width deliberately does not,
  # which is what keeps them reading as cracks rather than wedges.
  #
  # Severity is shown as a relative low-to-high scale only. An absolute figure
  # was tried and pulled: the honest one is "hectares that changed dominant
  # species", which is easily misread as hectares burned, and the fire tables
  # that would justify the latter belong to separate `onlyfire` simulations.
  # When matched fire output exists, this key is where the real number goes.
  sy <- 0.400
  text(L, sy + 0.030, "Scar severity", col = ink_sec, adj = c(0, 0), cex = sz_label)
  rect(L, sy - 0.036, L + 0.40, sy + 0.014,
       col = SPECIES_COLOURS_DARK_6["Mixed.spruce"], border = NA)
  for (j in 1:3) {
    xx <- L + 0.07 + (j - 1) * 0.13
    len <- 0.012 + 0.032 * c(0.25, 0.6, 1)[j]
    polygon(c(xx - 0.0035, xx + 0.0035, xx + 0.0012, xx - 0.0012),
            c(sy + 0.012, sy + 0.012, sy + 0.012 - len, sy + 0.012 - len),
            col = lesion_col, border = NA)
  }
  text(L + 0.07, sy - 0.046, "low", col = ink_mut, adj = c(0.5, 1), cex = sz_axis)
  text(L + 0.33, sy - 0.046, "high", col = ink_mut, adj = c(0.5, 1), cex = sz_axis)

  # ---- how to read ----
  # The panel bottom is a hard edge -- keep this to ~12 lines including blanks
  # or the last paragraph runs off the frame.
  # Plain language over precision of phrasing, but not at the cost of accuracy.
  # "A percentile rank across all three GCMs" was correct and unreadable; the
  # caveat it carried -- that width is ranked, not measured -- still has to
  # survive, because a wide ring genuinely does not mean "this many mm".
  # 18 lines including blanks. The panel bottom is a hard edge: at sz_axis a line
  # is ~0.0133 of panel height, so from 0.282 there is room for about 20.
  text(L, 0.320, "How to read it", col = ink_pri, adj = c(0, 1), cex = sz_body)
  text(L, 0.282,
       paste0("Each ring is one year, counted outward from the centre. The\n",
              "wedges show how much of the landscape each species\n",
              "dominates that year.\n\n",
              "Ring width is rainfall — wide rings wet, narrow rings dry — and\n",
              "narrows further when the forest changes species sharply. A\n",
              "narrow ring WITH a scar means disruption; one without means\n",
              "simply a dry year. Widths are ranked against all three climate\n",
              "models rather than measured, so wide means wetter than most\n",
              "years, not wetter by a set amount.\n\n",
              "The slot at the top is the year axis, coloured by temperature\n",
              "(smoothed, so warming shows through the noise). Because\n",
              "widths vary, distance from the centre is not proportional to\n",
              "time — read the years there.\n\n",
              "Scars mark the species that lost the most ground that year:\n",
              "position shows WHICH species, not where. Severity is species\n",
              "change, not area burned."),
       col = ink_mut, adj = c(0, 1), cex = sz_axis)
}

# ---- triptych ----------------------------------------------------------------

# Three discs in a triangle with one shared legend, for comparing along a single
# axis. ADDITIVE: render_disc_frame() is untouched and remains the primary output.
#
# A 3 x 3 grid was rejected on measurement, not taste. At 1920 px square a cell
# disc has a 301 px radius, which puts the mean ring at 3.2 px and the thinnest at
# 1.1 px -- under the ~2 px floor where antialiasing dissolves them, and the 87
# hairlines merge into grey. Matching today's single-disc quality in a 3 x 3 needs
# a 3840 px canvas. Three discs get there at 3000 px, and three panels is also
# about as much as one figure can be read at once.
#
# `vary` is limited to model / replicate / landscape because SSP is not an option:
# ssp370 exists only in the onlyfire runs, so the vegetation output is ssp245 only.
TRIPTYCH_VARY <- c("model", "replicate", "landscape")

render_triptych <- function(path, vary = "model",
                            landscape = "landscape_alaska_01_2015-2100scenario",
                            model = "TaiESM1", replicate = 1,
                            values = NULL,
                            cl = load_climate(), sc = shared_scales(cl),
                            width = 3000, height = 2850) {
  vary <- match.arg(vary, TRIPTYCH_VARY)
  if (is.null(values)) {
    values <- switch(vary,
      model     = c("NorEsm2-MM", "TaiESM1", "UKESM1-0-LL"),
      replicate = c(1, 2, 3),
      landscape = SCENARIO_LANDSCAPES)
  }
  if (length(values) != 3) stop("A triptych needs exactly 3 values, got ",
                                length(values))

  specs <- lapply(values, function(v) {
    s <- list(landscape = landscape, model = model, replicate = replicate)
    s[[vary]] <- if (vary == "replicate") as.integer(v) else v
    s
  })
  dd <- lapply(specs, function(s) disc_data(s$landscape, s$model, s$replicate, cl))

  lab_of <- function(s) switch(vary,
    model     = s$model,
    replicate = paste("replicate", s$replicate),
    landscape = paste("Landscape",
                      sub("^landscape_alaska_(\\d+)_.*$", "\\1", s$landscape)))

  agg_png(path, width = width, height = height, background = bg)
  on.exit(invisible(dev.off()), add = TRUE)
  par(bg = bg)

  # Equal square regions so all three discs render at identical scale -- a
  # comparison figure whose panels differ in size compares nothing.
  side <- 1270
  w <- side / width; h <- side / height
  regions <- list(c(0.05, 0.05 + w, 0.44, 0.44 + h),          # top left
                  c(0.95 - w, 0.95, 0.44, 0.44 + h),          # top right
                  c(0.5 - w / 2, 0.5 + w / 2, 0.005, 0.005 + h))  # bottom centre

  for (i in 1:3) {
    par(fig = regions[[i]], mar = c(0, 0, 0, 0), new = (i > 1))
    plot.new(); plot.window(c(0, 1), c(0, 1), asp = 1)
    draw_disc(dd[[i]]$m, dd[[i]]$di, dd[[i]]$pr, dd[[i]]$tmean, sc)
    # Top-left corner of the region. These regions are SQUARE, so asp = 1 gives
    # no vertical expansion and anything at negative y is clipped away -- unlike
    # the single-disc frame, whose taller region made that space real. The
    # corners of a square holding an inscribed circle are free, so use one.
    text(0.02, 0.98, lab_of(specs[[i]]), col = ink_pri, adj = c(0, 1),
         cex = sz_body * 1.2)
    text(0.02, 0.935,
         sprintf("%+.1f°C over %d-%d", tail(dd[[i]]$tmean, 1) - dd[[i]]$tmean[1],
                 min(dd[[i]]$years), max(dd[[i]]$years)),
         col = ink_mut, adj = c(0, 1), cex = sz_axis * 1.15)
  }

  # ---- title band ----
  par(fig = c(0.04, 0.96, 0.89, 1), mar = c(0, 0, 0, 0), new = TRUE)
  plot.new(); plot.window(c(0, 1), c(0, 1))
  text(0, 0.72, "87 years of a boreal landscape", col = ink_pri,
       adj = c(0, 0.5), cex = sz_title * 1.35)
  fixed <- switch(vary,
    model     = sprintf("Landscape %s · ssp245 · replicate %d",
                        sub("^landscape_alaska_(\\d+)_.*$", "\\1", landscape),
                        replicate),
    replicate = sprintf("Landscape %s · %s · ssp245",
                        sub("^landscape_alaska_(\\d+)_.*$", "\\1", landscape),
                        model),
    landscape = sprintf("%s · ssp245 · replicate %d", model, replicate))
  text(0, 0.26, paste0("iLand · ", fixed, "  —  compared by ", vary),
       col = ink_sec, adj = c(0, 0.5), cex = sz_body * 1.15)

  # ---- species key, bottom left ----
  par(fig = c(0.02, 0.29, 0.03, 0.43), mar = c(0, 0, 0, 0), new = TRUE)
  plot.new(); plot.window(c(0, 1), c(0, 1))
  for (j in seq_along(SPECIES_ORDER_6)) {
    yy <- 0.93 - (j - 1) * 0.085
    rect(0.02, yy - 0.028, 0.16, yy + 0.028,
         col = SPECIES_COLOURS_DARK_6[j], border = NA)
    text(0.21, yy, SPECIES_LABELS[SPECIES_ORDER_6[j]], col = ink_sec,
         adj = c(0, 0.5), cex = sz_label * 1.25)
  }
  # ONE large locator for every variant, always here. Three small per-disc maps
  # were tried for the landscape comparison and dropped: a single big map with
  # numbered markers is more legible, and it keeps the layout identical across
  # all three variants. The species key ends at y ~0.48, so the whole lower half
  # of this region is free.
  draw_locator(load_ak(), load_interior(), load_footprints(),
               current = sub("_\\d{4}-\\d{4}scenario$", "",
                             if (vary == "landscape") unlist(values) else landscape),
               x0 = 0.02, y0 = 0.035, w = 0.62, h = 0.40)

  # ---- keys and short caption, bottom right ----
  par(fig = c(0.71, 0.98, 0.03, 0.43), mar = c(0, 0, 0, 0), new = TRUE)
  plot.new(); plot.window(c(0, 1), c(0, 1))
  text(0, 0.96, "Temperature trend", col = ink_sec, adj = c(0, 1),
       cex = sz_label * 1.25)
  tseq <- seq(sc$t_lo, sc$t_hi, length.out = 120)
  xs <- seq(0, 0.86, length.out = 121)
  for (j in seq_along(tseq)) {
    rect(xs[j], 0.84, xs[j + 1], 0.90, col = temp_colour(tseq[j], sc),
         border = NA)
  }
  if (sc$t_lo < 0 && sc$t_hi > 0) {
    xz <- 0.86 * (0 - sc$t_lo) / (sc$t_hi - sc$t_lo)
    segments(xz, 0.825, xz, 0.915, col = ink_pri, lwd = 2)
  }
  text(0, 0.80, paste0(round(sc$t_lo, 1), "°C"), col = ink_mut, adj = c(0, 1),
       cex = sz_axis * 1.2)
  text(0.86, 0.80, paste0("+", round(sc$t_hi, 1), "°C"), col = ink_mut,
       adj = c(1, 1), cex = sz_axis * 1.2)

  text(0, 0.66, "Scar severity", col = ink_sec, adj = c(0, 1), cex = sz_label * 1.25)
  rect(0, 0.50, 0.86, 0.62, col = SPECIES_COLOURS_DARK_6["Mixed.spruce"],
       border = NA)
  for (j in 1:3) {
    xx <- 0.15 + (j - 1) * 0.28
    len <- 0.03 + 0.075 * c(0.25, 0.6, 1)[j]
    polygon(c(xx - 0.008, xx + 0.008, xx + 0.003, xx - 0.003),
            c(0.615, 0.615, 0.615 - len, 0.615 - len), col = lesion_col,
            border = NA)
  }
  text(0.15, 0.47, "low", col = ink_mut, adj = c(0.5, 1), cex = sz_axis * 1.2)
  text(0.71, 0.47, "high", col = ink_mut, adj = c(0.5, 1), cex = sz_axis * 1.2)

  # Deliberately shorter than the single-disc caption -- three discs already ask
  # a lot of the eye. The full reading guide lives in docs/tree-ring-disc.md.
  text(0, 0.36,
       paste0("One ring per year, outward from the centre; wedges are the\n",
              "share of landscape each species dominates.\n\n",
              "Width is rainfall, narrowed further where the forest changes\n",
              "species sharply. Cracks mark the most disrupted years.\n\n",
              "All three discs share one scale, so they can be compared\n",
              "directly. Radius is not proportional to time."),
       col = ink_mut, adj = c(0, 1), cex = sz_axis * 1.2)

  invisible(path)
}

if (!isTRUE(getOption("iland_rings_no_run", FALSE))) {
  dir.create(here("render"), recursive = TRUE, showWarnings = FALSE)
  cl <- load_climate(); sc <- shared_scales(cl)
  # Both GCMs on the SAME scales -- rendering them together is the check that
  # the shared normalisation works.
  for (g in getOption("iland_rings_models", c("NorEsm2-MM", "TaiESM1"))) {
    out <- here("render", paste0("iland_rings_L01_", g, "_rep1.png"))
    render_disc_frame(out, model = g, cl = cl, sc = sc)
    cat("wrote", out, "\n")
  }
  # Triptychs are opt-in so the default source() keeps producing exactly the
  # single-disc figures it always has:
  #   options(iland_rings_triptych = c("model", "landscape")); source(...)
  for (v in getOption("iland_rings_triptych", character(0))) {
    out <- here("render", paste0("iland_triptych_", v, ".png"))
    render_triptych(out, vary = v, cl = cl, sc = sc)
    cat("wrote", out, "\n")
  }
}
