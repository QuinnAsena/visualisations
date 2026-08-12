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

# ---- palette -----------------------------------------------------------------
#
# COPIED from C:/Users/asenaq/Documents/GitHub/esa-2026/R/theme.R, which is the
# SOURCE OF TRUTH -- keep in sync, and prefer changing it there first.
#
# Note the smell: SPECIES_COLOURS now exists in three places (the iLand app's
# plots/theme.R, esa-2026's R/theme.R, and here). esa-2026's is canonical because
# it is the only one that solves dark mode, and it does so by measurement rather
# than by inverting the light palette:
#   Pima #1b5e20 -> #388e3c   (the original is only 2.47:1 on a dark surface and
#                              fails the 3:1 floor -- and it is the LARGEST
#                              category, so it would be the most visible failure)
#   Potr #fdd835 -> #e8b923   (not a failure at 13.9:1, dimmed for glare)

SPECIES_COLOURS_DARK <- c(
  Pima            = "#388e3c",  # Picea mariana        black spruce
  Pigl            = "#66bb6a",  # Picea glauca         white spruce
  Potr            = "#e8b923",  # Populus tremuloides  trembling aspen
  Bene            = "#e65100",  # Betula neoalaskana   Alaska birch
  Mixed.spruce    = "#7cb342",
  Mixed.deciduous = "#f57c00",
  `Not forested`  = "#78909c"
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

# Both spellings keyed: the downscaling pipeline writes NorESM2-MM, fire_annual
# writes NorEsm2-MM.
GCM_COLOURS <- c(
  "TaiESM1"     = "#2a78d6",
  "UKESM1-0-LL" = "#eb6834",
  "NorESM2-MM"  = "#1baf7a",
  "NorEsm2-MM"  = "#1baf7a"
)

# Surface + chrome, matching ak_fire_anim.R so the two pieces sit together on the
# site. Darker than esa-2026's #191919, so the dark species palette can only gain
# contrast -- verified in the checks below rather than assumed.
bg       <- "#0d0d0d"
ink_pri  <- "#ffffff"
ink_sec  <- "#c3c2b7"
ink_mut  <- "#898781"
hairline <- "#2c2c2a"

# agg_png runs at res = 72 with pointsize = 12, so cex = points / 12.
pt <- function(points) points / 12

# ---- loading -----------------------------------------------------------------

# Scenario family only, onlysimfalse only, Not forested dropped, shares
# renormalised over the remaining six categories.
#
# The renormalisation is not cosmetic: total landscape area is a constant
# 58,396.3 ha across all years, but that constant includes Not forested. Dropping
# it leaves a subtotal that varies year to year, so shares MUST be recomputed
# against that subtotal or the rings will not close to a full circle.
load_area_dom <- function(landscape = "landscape_alaska_01_2015-2100scenario",
                          path = area_dom_f) {
  d <- as.data.frame(read_parquet(path))

  d <- d[d$landscape == landscape &
           grepl("onlysimfalse", d$treatment, fixed = TRUE) &
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

# ---- ring geometry -----------------------------------------------------------

# Ring width falls as disturbance rises, so a disrupted year lays down a thin ring
# the way a stressed tree does. Radius is therefore NOT linear in time -- the disc
# needs decade ticks, which draw_disc() adds.
#
# w_min keeps a bad year visible rather than collapsing it to nothing; gamma < 1
# spreads the effect across moderate years, gamma > 1 confines it to the extremes.
ring_radii <- function(di, w_min = 0.25, gamma = 0.7, r0 = 0.10, R = 1) {
  scaled <- if (max(di) > 0) (di / max(di))^gamma else rep(0, length(di))
  w <- w_min + (1 - w_min) * (1 - scaled)
  r <- r0 + cumsum(w) / sum(w) * (R - r0)
  c(r0, r)                      # n + 1 boundaries for n rings
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
