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
disc_split <- 0.60

w_min <- 0.25            # thinnest ring as a fraction of the thickest
gamma <- 0.7
r0    <- 0.085           # pith

# Pima and Pigl are adjacent in stack order and only 14.3 dE apart -- below the
# 15 floor -- so they must not abut. A surface gap is the prescribed remedy;
# re-stepping a colour esa-2026 already validated is not.
sector_gap <- 0.0025     # radians inset each side of a species sector

gutter <- 6 * pi / 180   # angular wedge reserved at 12 o'clock for the axis

ring_line_col <- "#0d0d0d"
ring_line_a   <- 0.30    # alpha of the year boundary hairline
ring_line_lwd <- 0.7

# Wounds displace boundaries by a multiple of the MEAN ring width, not by a
# fraction of absolute radius -- scaling to radius made the outer rings balloon
# into lobes because the same fraction is a far bigger distance out there.
scar_top_n      <- 8
scar_sigma      <- 0.13  # angular spread (radians)
scar_depth_rings <- 2.4  # displacement in units of mean ring width
scar_decay      <- 0.88  # per-ring persistence of an older wound

# Deflection alone is not enough: a bent ring boundary is indistinguishable from
# a birch sector, because both read as radial features. A real fire scar shows as
# a dark lesion the tree grows over, so draw one -- charred brown rather than
# black, so it reads as burnt tissue and not as a hole punched in the disc.
#
# Width is a constant PHYSICAL half-width converted to an angle at each radius,
# not a constant angle. An angular width tapering outward produced little dark
# chevrons that read as floating debris; a constant-width radial crack reads as a
# wound. Length scales with severity, so the worst years cut furthest.
lesion_col    <- "#241408"
lesion_alpha  <- 0.90
lesion_halfw  <- 0.0055  # half-width in disc-radius units
lesion_min    <- 3       # rings a scar runs outward, at minimum severity
lesion_max    <- 11      # ... and at maximum

n_theta <- 1441          # deflection grid resolution

sz_title <- pt(30)
sz_body  <- pt(21)
sz_label <- pt(18)
sz_axis  <- pt(16)

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

draw_disc <- function(m, di, cx = 0.5, cy = 0.5, scale = 0.455) {
  n <- nrow(m)
  years <- as.integer(rownames(m))
  r <- ring_radii(di, w_min = w_min, gamma = gamma, r0 = r0)

  # Scars sit at the angular midpoint of the sector that lost the most share
  # that year -- the wound lands where the damage did.
  sc <- scar_species(m)
  keep <- order(-di)[seq_len(min(scar_top_n, n))]
  keep <- keep[di[keep] > 0]
  span <- 2 * pi - gutter
  a_start <- pi / 2 - gutter / 2

  scars <- data.frame(ring = integer(0), angle = numeric(0), amp = numeric(0))
  if (length(keep) > 0) {
    ang <- vapply(keep, function(k) {
      cum <- c(0, cumsum(m[k, ]))
      j <- match(sc$sp[k], SPECIES_ORDER_6)
      a_start - span * (cum[j] + cum[j + 1]) / 2
    }, numeric(1))
    scars <- data.frame(ring = keep, angle = ang %% (2 * pi),
                        amp = di[keep] / max(di))
  }

  theta_grid <- seq(0, 2 * pi, length.out = n_theta)
  fld <- boundary_field(r, scars, theta_grid)

  line_col <- grDevices::adjustcolor(ring_line_col, alpha.f = ring_line_a)

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
          col = line_col, lwd = ring_line_lwd)
  }

  # Pith
  tt <- seq(0, 2 * pi, length.out = 200)
  polygon(cx + r0 * 0.92 * cos(tt) * scale, cy + r0 * 0.92 * sin(tt) * scale,
          col = bg, border = NA)

  # ---- the gutter and its year axis ----
  gth <- seq(a_start, a_start + gutter, length.out = 40)
  outer_r <- max(fld[n + 1, ]) * 1.002
  polygon(c(cx + r0 * 0.9 * cos(gth) * scale, cx + rev(outer_r * cos(gth)) * scale),
          c(cy + r0 * 0.9 * sin(gth) * scale, cy + rev(outer_r * sin(gth)) * scale),
          col = bg, border = NA)

  amid <- a_start + gutter / 2
  dec <- which(years %% 20 == 0)
  for (k in dec) {
    rr <- r[k]
    lines(cx + c(rr, rr * 1.0) * cos(amid) * scale,
          cy + c(rr, rr) * sin(amid) * scale)
    points(cx + rr * cos(amid) * scale, cy + rr * sin(amid) * scale,
           pch = 16, cex = 0.55, col = ink_sec)
    text(cx + rr * cos(amid) * scale - 0.012, cy + rr * sin(amid) * scale,
         years[k], col = ink_sec, adj = c(1, 0.5), cex = sz_axis)
  }
  # First year anchors the pith. The last year is skipped when it is already a
  # decade tick -- 2100 was being drawn twice, on top of itself.
  text(cx + r[1] * cos(amid) * scale - 0.012, cy + r[1] * sin(amid) * scale,
       years[1], col = ink_mut, adj = c(1, 0.5), cex = sz_axis)
  if (!(years[n] %% 20 == 0)) {
    rr <- max(fld[n + 1, ])
    text(cx + rr * cos(amid) * scale - 0.012, cy + rr * sin(amid) * scale,
         years[n], col = ink_mut, adj = c(1, 0.5), cex = sz_axis)
  }

  invisible(list(r = r, scars = scars))
}

# ---- full frame --------------------------------------------------------------

render_disc_frame <- function(path,
                              landscape = "landscape_alaska_01_2015-2100scenario",
                              model = "NorEsm2-MM", replicate = 1) {
  d <- load_area_dom(landscape)
  tr <- unique(d$treatment[d$model == model])
  if (length(tr) != 1) stop("Expected one treatment for model ", model)
  m <- share_matrix(d, tr, replicate)
  di <- disturbance_index(m)
  years <- as.integer(rownames(m))

  agg_png(path, width = width, height = height, background = bg)
  on.exit(invisible(dev.off()), add = TRUE)
  par(bg = bg)

  # asp = 1 makes one x unit equal one y unit in device space, so the disc is a
  # circle rather than an ellipse stretched by the region's aspect.
  par(fig = c(0, disc_split, 0, 1), mar = c(0, 0, 0, 0), new = FALSE)
  plot.new(); plot.window(c(0, 1), c(0, 1), asp = 1)
  draw_disc(m, di)

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

  ky <- 0.665
  for (j in seq_along(SPECIES_ORDER_6)) {
    yy <- ky - (j - 1) * 0.050
    rect(L, yy - 0.016, L + 0.052, yy + 0.016,
         col = SPECIES_COLOURS_DARK_6[j], border = NA)
    text(L + 0.072, yy, SPECIES_LABELS[SPECIES_ORDER_6[j]], col = ink_sec,
         adj = c(0, 0.5), cex = sz_label)
  }

  text(L, 0.315, "How to read it", col = ink_pri, adj = c(0, 1), cex = sz_body)
  text(L, 0.272,
       paste0("One ring per year, counted outward from the\n",
              "pith. Each ring is divided by the share of\n",
              "landscape area each species dominates.\n\n",
              "Rings narrow when composition shifts hard, so\n",
              "tightly spaced lines mark disrupted decades —\n",
              "radius is not proportional to time.\n\n",
              "Wounds mark the ", scar_top_n, " most disrupted years, set at\n",
              "the species that lost the most ground. This is\n",
              "a derived index of composition change, not\n",
              "observed area burned."),
       col = ink_mut, adj = c(0, 1), cex = sz_axis)
}

if (!isTRUE(getOption("iland_rings_no_run", FALSE))) {
  out <- here("render", "iland_rings_L01_NorEsm2_rep1.png")
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  render_disc_frame(out)
  cat("wrote", out, "\n")
}
