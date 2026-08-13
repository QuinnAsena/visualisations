# Tree-ring composition disc: iLand species composition 2015-2101 as a cross-section.
#
# Radius     = time (NOT linear -- ring width carries data, so decade ticks are
#              mandatory rather than decorative)
# Angle      = share of landscape area by dominant species, in SPECIES_ORDER_6
# Ring width = inverse of a derived disturbance index, so a disrupted year lays
#              down a thin ring the way a stressed tree does
# Ring line  = a hairline at every year boundary. This is what makes the thing
#              read as a core rather than a wobbly pie chart: in a real core it is
#              the latewood band you count, and tight line spacing IS the signal.
# Scar       = radial wound at the sector of whichever species lost the most share
#              that year, with later rings deflected around it
# Gutter     = a wedge cut at 12 o'clock carrying the year axis. Doubles as the
#              slot a increment borer would leave.

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
gamma      <- 0.7
r0         <- 0.085      # pith

# Pima and Pigl are adjacent in stack order and only 14.3 dE apart -- below the
# 15 floor -- so they must not abut. A surface gap is the prescribed remedy;
# re-stepping a colour esa-2026 already validated is not.
sector_gap <- 0.0025     # radians inset each side of a species sector

gutter <- 7 * pi / 180   # angular wedge at 12 o'clock: the climate strip + axis

temp_smooth_k <- 11      # running-mean window for the temperature strip

ring_line_a   <- 0.32    # alpha of the year boundary hairline
ring_line_lwd <- 0.7

# Wounds displace boundaries by a multiple of the MEAN ring width, not by a
# fraction of absolute radius -- scaling to radius made the outer rings balloon
# into lobes because the same fraction is a far bigger distance out there.
scar_top_n      <- 8
scar_sigma      <- 0.13  # angular spread (radians)
scar_depth_rings <- 6.0  # displacement in units of mean ring width. Deep enough
                         # that fires visibly dent the outline -- the disc is not
                         # meant to stay circular. The cummax guard in
                         # boundary_field() is load-bearing at this depth.
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
                   dist_pen = dist_pen, gamma = gamma)
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
    scars <- data.frame(ring = keep, angle = ang %% (2 * pi),
                        amp = di[keep] / max(di))
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

  tcols <- temp_colour(tmean, sc)
  for (k in seq_len(n)) {
    polygon(c(cx + r[k] * cos(gth) * scale, cx + rev(r[k + 1] * 1.02 * cos(gth)) * scale),
            c(cy + r[k] * sin(gth) * scale, cy + rev(r[k + 1] * 1.02 * sin(gth)) * scale),
            col = tcols[k], border = NA)
  }

  # Year labels sit on a background plate so they stay legible over whatever the
  # strip or the fills happen to be at that radius -- previously they were muted
  # grey over an orange wedge and effectively invisible.
  amid <- a_start + gutter * 1.35
  lab_years <- years[years %% 20 == 0]
  for (yy in c(years[1], lab_years)) {
    k <- match(yy, years)
    rr <- r[k]
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

# ---- full frame --------------------------------------------------------------

render_disc_frame <- function(path,
                              landscape = "landscape_alaska_01_2015-2100scenario",
                              model = "NorEsm2-MM", replicate = 1,
                              cl = load_climate(), sc = climate_scales(cl)) {
  d <- load_area_dom(landscape)
  tr <- unique(d$treatment[d$model == model])
  if (length(tr) != 1) stop("Expected one treatment for model ", model)
  m <- share_matrix(d, tr, replicate)
  di <- disturbance_index(m)
  years <- as.integer(rownames(m))

  # canon_gcm bridges NorEsm2-MM (iLand) and NorESM2-MM (climate). Without it
  # this join silently returns nothing.
  cg <- cl[cl$gcm == canon_gcm(model), ]
  idx <- match(years, cg$year)
  if (anyNA(idx)) stop("Climate missing for years: ",
                       paste(years[is.na(idx)], collapse = ", "))
  pr <- cg$pr[idx]
  # Smoothed for colour so the trend reads instead of the annual noise; the
  # caption says so. Precipitation stays raw -- its year-to-year jitter is the
  # texture we want in the ring widths.
  tmean <- smooth_trend(cg$tmean[idx], years)

  agg_png(path, width = width, height = height, background = bg)
  on.exit(invisible(dev.off()), add = TRUE)
  par(bg = bg)

  # asp = 1 makes one x unit equal one y unit in device space, so the disc is a
  # circle rather than an ellipse stretched by the region's aspect.
  par(fig = c(0, disc_split, 0, 1), mar = c(0, 0, 0, 0), new = FALSE)
  plot.new(); plot.window(c(0, 1), c(0, 1), asp = 1)
  draw_disc(m, di, pr, tmean, sc)

  par(fig = c(disc_split, 1, 0, 1), mar = c(0, 0, 0, 0), new = TRUE)
  plot.new(); plot.window(c(0, 1), c(0, 1))
  L <- 0.06

  text(L, 0.945, paste(length(years), "years"), col = ink_pri,
       adj = c(0, 1), cex = sz_title)
  text(L, 0.882, "of a boreal landscape", col = ink_pri, adj = c(0, 1), cex = sz_title)
  text(L, 0.800,
       paste0("iLand · Landscape ",
              sub("^landscape_alaska_(\\d+)_.*$", "\\1", landscape),
              " · ", min(years), "-", max(years)),
       col = ink_sec, adj = c(0, 1), cex = sz_body)
  text(L, 0.757, paste0(model, " · ssp245 · replicate ", replicate),
       col = ink_mut, adj = c(0, 1), cex = sz_body)

  ky <- 0.700
  for (j in seq_along(SPECIES_ORDER_6)) {
    yy <- ky - (j - 1) * 0.044
    rect(L, yy - 0.014, L + 0.048, yy + 0.014,
         col = SPECIES_COLOURS_DARK_6[j], border = NA)
    text(L + 0.066, yy, SPECIES_LABELS[SPECIES_ORDER_6[j]], col = ink_sec,
         adj = c(0, 0.5), cex = sz_label)
  }

  # ---- temperature key (the gutter strip) ----
  ty <- 0.400
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
  # mark freezing, which is where the ramp pivots and the strip flips
  xz <- L + 0.40 * (0 - sc$t_lo) / (sc$t_hi - sc$t_lo)
  segments(xz, ty - 0.019, xz, ty + 0.019, col = ink_pri, lwd = 2)
  text(L, ty - 0.026, paste0(round(sc$t_lo, 1), "°C"), col = ink_mut,
       adj = c(0, 1), cex = sz_axis)
  text(xz, ty - 0.026, "0", col = ink_sec, adj = c(0.5, 1), cex = sz_axis)
  text(L + 0.40, ty - 0.026, paste0("+", round(sc$t_hi, 1), "°C"),
       col = ink_mut, adj = c(1, 1), cex = sz_axis)


  # ---- severity key ----
  # Scar LENGTH and depth scale with disturbance; width deliberately does not,
  # which is what keeps them reading as cracks rather than wedges.
  sy <- 0.290
  text(L, sy + 0.030, "Scar severity", col = ink_sec, adj = c(0, 0), cex = sz_label)
  rect(L, sy - 0.036, L + 0.40, sy + 0.014,
       col = SPECIES_COLOURS_DARK_6["Mixed.spruce"], border = NA)
  for (j in 1:3) {
    frac <- c(0.25, 0.6, 1)[j]
    xx <- L + 0.07 + (j - 1) * 0.13
    len <- 0.012 + 0.032 * frac
    polygon(c(xx - 0.0035, xx + 0.0035, xx + 0.0012, xx - 0.0012),
            c(sy + 0.012, sy + 0.012, sy + 0.012 - len, sy + 0.012 - len),
            col = lesion_col, border = NA)
  }
  text(L + 0.07, sy - 0.046, "low", col = ink_mut, adj = c(0.5, 1), cex = sz_axis)
  text(L + 0.33, sy - 0.046, "high", col = ink_mut, adj = c(0.5, 1), cex = sz_axis)

  # ---- how to read ----
  # The panel bottom is a hard edge -- keep this to ~11 lines including blanks
  # or the last paragraph runs off the frame.
  text(L, 0.205, "How to read it", col = ink_pri, adj = c(0, 1), cex = sz_body)
  text(L, 0.166,
       paste0("One ring per year, outward from the pith, split by the share\n",
              "of landscape area each species dominates.\n\n",
              "Width is precipitation — a percentile rank across all three\n",
              "GCMs (457–732 mm), so relative and non-linear — narrowed\n",
              "further in disrupted years. A thin ring WITH a scar is\n",
              "disturbance; a thin ring without one is a dry year.\n\n",
              "Slot colour is temperature as a GAM trend, so the warming\n",
              "shows through the year-to-year noise. Scars mark the\n",
              "species that lost most ground: the angle says WHICH species,\n",
              "not where. Radius is not proportional to time."),
       col = ink_mut, adj = c(0, 1), cex = sz_axis)
}

if (!isTRUE(getOption("iland_rings_no_run", FALSE))) {
  dir.create(here("render"), recursive = TRUE, showWarnings = FALSE)
  cl <- load_climate(); sc <- climate_scales(cl)
  # Both GCMs on the SAME scales -- rendering them together is the check that
  # the shared normalisation works.
  for (g in getOption("iland_rings_models", c("NorEsm2-MM", "TaiESM1"))) {
    out <- here("render", paste0("iland_rings_L01_", g, "_rep1.png"))
    render_disc_frame(out, model = g, cl = cl, sc = sc)
    cat("wrote", out, "\n")
  }
}
