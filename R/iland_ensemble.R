# Ensemble-scale loaders for iLand scenario output.
#
# The disc (R/iland_rings.R) is a SINGLE-REALISATION object: one landscape, one
# model, one replicate. That is deliberate, and it is also its one structural
# limit -- a disc spends angle on species and radius on time, so no axis is left
# for an ensemble. This module exists for the opposite kind of figure, whose
# subject IS the spread across replicates. Everything here returns the bundle.
#
# Reuses iland_common.R for the palette, the chrome and the per-run maths.
# Nothing here changes how a disc renders.

library(here)
source(here("R", "iland_common.R"))

# ---- fire condition ----------------------------------------------------------

# The naming is inverted from intuition, so state it once:
#
#   onlysimfalse -- the real runs. Fire spreads AND kills trees. 12 replicates.
#   onlysimtrue  -- iLand's <onlySimulation> test mode. Fire spreads but removes
#                   no biomass, so these are near-deterministic BY CONSTRUCTION.
#                   9 replicates.
#
# That makes onlysimtrue a free no-feedback CONTROL. Drawn beside the real runs
# it shows how much of the ensemble spread is fire's doing -- exactly the thing a
# single disc cannot show, and it costs nothing because the runs are already in
# the local cache. load_area_dom() filters them out by default.
FIRE_CONDITIONS <- c("onlysimfalse", "onlysimtrue")
FIRE_LABELS <- c(onlysimfalse = "fire kills trees",
                 onlysimtrue  = "no fire feedback")
FIRE_N_REP  <- c(onlysimfalse = 12L, onlysimtrue = 9L)

# ---- GCM identity ------------------------------------------------------------

# Keyed on BOTH spellings on purpose: the climate files say NorESM2-MM, the iLand
# treatments say NorEsm2-MM. This is the "one-line lift from esa-2026/R/theme.R"
# that iland_common.R's comment anticipated.
#
# NOT AUDITED for this piece. These were tuned against esa-2026's #191919
# surface, not the #0d0d0d used here, and never checked for colour-vision
# deficiency. Fine for a prototype; measure before anything ships.
GCM_COLOURS <- c(
  "NorESM2-MM"  = "#1baf7a",   # aqua
  "NorEsm2-MM"  = "#1baf7a",
  "TaiESM1"     = "#2a78d6",   # blue
  "UKESM1-0-LL" = "#eb6834"    # orange
)

# Coolest to hottest, by GAM-smoothed 2014-2100 warming: +2.0 / +4.8 / +6.2 degC.
GCM_ORDER <- c("NorEsm2-MM", "TaiESM1", "UKESM1-0-LL")

# ---- species groupings -------------------------------------------------------

# The six categories are exclusive and split three-and-three. Mixed.spruce is
# Pima/Pigl co-dominant and Mixed.deciduous is Bene/Potr co-dominant, so
# "deciduousness" is a clean sum rather than a weighted guess.
CONIFER_SP <- c("Pima", "Pigl")
DECID_SP   <- c("Potr", "Bene")
MIXED_SP   <- c("Mixed.spruce", "Mixed.deciduous")
DECID_ALL  <- c("Potr", "Bene", "Mixed.deciduous")

# The 1-D tipping axis: boreal spruce giving way to deciduous under fire is THE
# question these runs are asking, so it earns the one spatial axis available.
decid_share <- function(m) rowSums(m[, DECID_ALL, drop = FALSE])

# Three-part composition for a ternary projection.
tern_coords <- function(m) {
  cbind(conifer = rowSums(m[, CONIFER_SP, drop = FALSE]),
        decid   = rowSums(m[, DECID_SP,   drop = FALSE]),
        mixed   = rowSums(m[, MIXED_SP,   drop = FALSE]))
}

# Equilateral triangle: conifer bottom-left, mixed bottom-right, deciduous apex.
tern_xy <- function(tc) {
  cbind(x = tc[, "mixed"] + 0.5 * tc[, "decid"],
        y = (sqrt(3) / 2) * tc[, "decid"])
}

# Which category holds the most area each year -- used to colour a thread by what
# is actually dominant, so species survive a figure whose axes are spent on time
# and composition.
dom_sp <- function(m) SPECIES_ORDER_6[max.col(m, ties.method = "first")]

# ---- the ensemble ------------------------------------------------------------

# Every run for one landscape and one fire condition, as a flat list.
#
# The `d` argument lets a caller hand in an already-loaded frame:
# load_area_dom() re-reads a 129k-row parquet every call, and a figure wanting
# both fire conditions across three landscapes would otherwise read it six times.
ensemble <- function(landscape = "landscape_alaska_01_2015-2100scenario",
                     fire = "onlysimfalse", d = NULL) {
  fire <- match.arg(fire, FIRE_CONDITIONS)
  if (is.null(d)) d <- load_area_dom(landscape, fire = fire)
  out <- list()
  for (tr in unique(d$treatment)) {
    mdl <- sub("ssp.*$", "", tr)
    for (rp in sort(unique(d$replicate[d$treatment == tr]))) {
      m <- share_matrix(d, tr, rp)
      out[[paste(mdl, rp, sep = "_")]] <- list(
        m = m, model = mdl, replicate = rp, fire = fire,
        years = as.integer(rownames(m)),
        decid = decid_share(m), di = disturbance_index(m), dom = dom_sp(m))
    }
  }
  mdls <- vapply(out, function(r) r$model, character(1))
  reps <- vapply(out, function(r) r$replicate, integer(1))
  out[order(match(mdls, GCM_ORDER), reps)]
}

# Runs for one model, preserving replicate order.
by_model <- function(e, model) e[vapply(e, function(r) r$model == model, TRUE)]

# ---- climate -----------------------------------------------------------------

# Annual precipitation in mm/yr and the GAM-smoothed temperature trend, aligned
# to a run's years. Regional, not per-landscape -- process_gcms.R takes a spatial
# mean over the whole downscaled raster.
climate_for <- function(model, years, cl = load_climate()) {
  cg <- cl[cl$gcm == canon_gcm(model), ]
  idx <- match(years, cg$year)
  if (anyNA(idx)) stop("Climate missing for years: ",
                       paste(years[is.na(idx)], collapse = ", "))
  list(pr = cg$pr[idx], tmean = smooth_trend(cg$tmean[idx], years))
}

# ---- encoded channels, with a stated gain ------------------------------------

# LINEAR, unlike the disc's percentile rank -- and that reversal is the point.
# The disc ranked rainfall because a linear map put the middle 50% of years into
# 0.236 of the width range, invisible with no axis and no key to recover it. A
# printed calibration bar makes a linear map readable as a quantity, so the rank
# is no longer needed and "wide means wetter by this many mm" becomes true again.
PR_W_MIN <- 0.28   # driest year's share of the widest

pr_width <- function(pr, sc, w_min = PR_W_MIN) {
  f <- (pr - sc$pr_lo) / (sc$pr_hi - sc$pr_lo)
  w_min + (1 - w_min) * pmin(pmax(f, 0), 1)
}

# Width units per mm of rainfall -- the SLOPE of the mapping above.
pr_gain <- function(sc, w_min = PR_W_MIN) (1 - w_min) / (sc$pr_hi - sc$pr_lo)

# ---- keys --------------------------------------------------------------------

# The printed calibration key, chosen over a fully labelled axis panel.
#
# An encoded channel is only quantitative if the reader is told what one unit of
# it means. A bar spanning a round number of real units does that in the figure's
# own visual language, where a bare number would not.
calib_bar <- function(x, y, len, label, vertical = FALSE,
                      col = ink_sec, cex = pt(15), lwd = 2, tick = 0.006) {
  if (vertical) {
    segments(x, y, x, y + len, col = col, lwd = lwd)
    segments(x - tick, c(y, y + len), x + tick, c(y, y + len), col = col, lwd = lwd)
    text(x + tick * 2.5, y + len / 2, label, col = col, adj = c(0, 0.5), cex = cex)
  } else {
    segments(x, y, x + len, y, col = col, lwd = lwd)
    segments(c(x, x + len), y - tick, c(x, x + len), y + tick, col = col, lwd = lwd)
    text(x + len / 2, y - tick * 2.5, label, col = col, adj = c(0.5, 1), cex = cex)
  }
  invisible(NULL)
}

# Calibrating an OFFSET linear channel -- the honest version of the key.
#
# pr_width() has a non-zero floor (PR_W_MIN) so the driest year stays visible.
# That makes width linear in rainfall but NOT proportional to it: no width
# corresponds to 100 mm, so a single bar captioned "width = 100 mm/yr" is a lie,
# and it also comes out about 10 px long, which is how the error announced
# itself. Two reference cross-sections drawn at the real extremes calibrate it
# truthfully -- the reader has two labelled anchors and interpolates between them,
# which is what a scale bar does anyway.
#
# `wmax` is the HALF-width the figure uses at maximum rainfall, in that figure's
# own units, so the swatches are drawn at exactly the size they appear.
width_key <- function(x, y, wmax, sc, label = "ribbon width = rainfall",
                      cex = pt(15), len = 0.048, gap = 0.034, col = ink_sec) {
  for (i in 1:2) {
    v  <- if (i == 1) sc$pr_lo else sc$pr_hi
    yy <- y + (i - 1) * gap
    hw <- pr_width(v, sc) * wmax
    rect(x, yy - hw, x + len, yy + hw, col = col, border = NA)
    text(x + len + 0.010, yy, sprintf("%.0f mm/yr", v), col = ink_mut,
         adj = c(0, 0.5), cex = cex)
  }
  text(x, y + gap + wmax + 0.014, label, col = ink_sec, adj = c(0, 0), cex = cex)
  invisible(NULL)
}

temp_key <- function(x, y, w, h, sc, cex = pt(15), label = "temperature trend") {
  tseq <- seq(sc$t_lo, sc$t_hi, length.out = 120)
  xs <- seq(x, x + w, length.out = 121)
  for (j in seq_along(tseq)) {
    rect(xs[j], y, xs[j + 1], y + h, col = temp_colour(tseq[j], sc), border = NA)
  }
  # Freezing is a reference tick, not a pivot -- the ramp is sequential.
  if (sc$t_lo < 0 && sc$t_hi > 0) {
    xz <- x + w * (0 - sc$t_lo) / (sc$t_hi - sc$t_lo)
    segments(xz, y - h * 0.35, xz, y + h * 1.35, col = ink_pri, lwd = 2)
  }
  text(x, y - h * 0.6, sprintf("%.1f°C", sc$t_lo), col = ink_mut,
       adj = c(0, 1), cex = cex)
  text(x + w, y - h * 0.6, sprintf("+%.1f°C", sc$t_hi), col = ink_mut,
       adj = c(1, 1), cex = cex)
  text(x, y + h * 1.9, label, col = ink_sec, adj = c(0, 0), cex = cex)
  invisible(NULL)
}

species_key <- function(x, y, dy, sw = 0.022, cex = pt(15), cols = 1,
                        colgap = 0.26) {
  n <- length(SPECIES_ORDER_6)
  per <- ceiling(n / cols)
  for (j in seq_len(n)) {
    ci <- (j - 1) %/% per
    ri <- (j - 1) %% per
    xx <- x + ci * colgap
    yy <- y - ri * dy
    rect(xx, yy - dy * 0.3, xx + sw, yy + dy * 0.3,
         col = SPECIES_COLOURS_DARK_6[j], border = NA)
    text(xx + sw * 1.5, yy, SPECIES_LABELS[SPECIES_ORDER_6[j]], col = ink_sec,
         adj = c(0, 0.5), cex = cex)
  }
  invisible(NULL)
}

gcm_key <- function(x, y, dy, cex = pt(15)) {
  for (j in seq_along(GCM_ORDER)) {
    yy <- y - (j - 1) * dy
    segments(x, yy, x + 0.030, yy, col = GCM_COLOURS[GCM_ORDER[j]], lwd = 3)
    text(x + 0.042, yy, GCM_ORDER[j], col = ink_sec, adj = c(0, 0.5), cex = cex)
  }
  invisible(NULL)
}
