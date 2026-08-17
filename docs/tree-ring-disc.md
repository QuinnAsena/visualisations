# The tree-ring disc — what it shows and every decision behind it

Plain-language documentation for `R/iland_rings.R` + `R/iland_common.R`.

This exists so the figure can be audited. Several encoding choices were mine
rather than specified, and **§4 lists every one of them** with the reasoning, the
alternative that was rejected, and what would have to change to undo it.

---

## 1. What you are looking at

A cross-section through a simulated boreal landscape, drawn like a tree core.
One disc = one iLand run: one landscape, one climate model, one replicate.

Each **ring is one year**, 2014 at the pith out to 2100 at the bark — 87 rings.
Within a ring, the landscape is split into **wedges by angle**, one per dominant
species category, sized by the share of landscape area that species dominates
that year. So reading outward shows composition changing; reading around a ring
shows the composition at one moment.

Five channels carry data:

| channel | what it encodes |
|---|---|
| angle within a ring | share of area dominated by each species |
| ring width | rainfall that year, narrowed further by species change |
| slot colour (12 o'clock) | temperature trend |
| dark radial cracks | the 8 most disrupted years |
| outline irregularity | accumulated scars, which never fully close |

---

## 2. Where the data comes from

Everything is copied locally first. Rendering never touches the network share,
because reading the climate files one at a time over it did not finish in ten
minutes.

| what | source | local cache |
|---|---|---|
| species composition | `iLand-visualisation-ai/data/processed/area_dom/dominant_species_area.parquet` | `iland_data/area_dom/` (1.1 MB) |
| climate | `quinn_storage/gcm_summaries/` — 1089 files | `iland_data/climate/gcm_annual.parquet` (2.2 MB) |
| landscape locations | `landscape_init_ak_can/{landscape}/gis/env.grid.tif` | `iland_data/landscapes/footprints.gpkg` (116 KB) |
| Alaska outlines | `ak_data/AK_no_islands.shp`, `AK_interior.shp` | already in the repo |

**Re-copy `area_dom` after any `preprocess.R` run** on the iLand project. The
climate and footprint caches rebuild themselves if deleted.

Filters applied: **scenario** family only (not spinup), **onlysimfalse** only
(the runs where fire actually kills trees), **ssp245** — which is the only SSP the
vegetation runs and the climate summaries have. `Not forested` is dropped, leaving
six categories.

---

## 3. How each channel works

### Angle — species share

Shares are **renormalised over the six kept categories**. This matters: total
landscape area is a constant 58,396.3 ha, but that constant includes
`Not forested`. Dropping it leaves a subtotal that moves year to year
(53,851–58,396 ha), so shares must be recomputed against that subtotal or the
rings do not close into a full circle.

Species are drawn in a fixed order from 12 o'clock, so sector boundaries can be
compared between neighbouring rings.

### Ring width — rainfall × disturbance

```
width  =  rainfall term  ×  disturbance term
```

Both suppress growth in a real tree, so a narrow ring is genuinely ambiguous —
**the scars are what disambiguate it.** Narrow with a scar means disruption;
narrow without one means a dry year. That is stated on the figure.

The consequence to keep in mind: **radius is no longer proportional to time.**
You cannot read a year off a distance; the year labels in the slot are the only
time axis.

### Slot colour — temperature

The wedge cut at 12 o'clock is both the year axis and a per-year temperature
strip. It doubles as the slot an increment borer would leave.

Temperature is **smoothed** before colouring, and the strip is labelled as a
trend. Raw annual values swamp the signal: year-to-year noise has a standard
deviation of 1.65 °C against a total rise of 2–6 °C.

### Scars — the disturbance index

A **derived index**, not observed fire. For each year:

```
index  =  ½ × Σ |share this year − share last year|   across the six categories
```

This is total variation distance. Halving is not arbitrary: shares sum to 1 in
both years, so gains and losses cancel exactly and the raw sum counts every
moved hectare twice. Halving gives the fraction that actually moved. An index of
0.044 means 4.4% of the landscape changed dominant species.

The eight highest-index years get a crack, positioned at the wedge of whichever
species lost the most share. Scars persist outward and dent the outline, because
a real catface never fully closes over.

**Why not real fire data?** The fire tables come from separate `onlyfire`
simulations. Their replicate 1 is a different simulation from the vegetation's
replicate 1, so their area burned belongs to another run — attaching it here
would look convincing and be wrong. This index is derived from the same numbers
that draw the rings, so it cannot disagree with them.

---

## 4. Decisions I made that you did not direct

### Rainfall enters as a rank, not a measurement

**What:** each year's rainfall is converted to its percentile position within the
pooled distribution across all three climate models, then mapped to ring width.

**Why:** annual rainfall is near-normal with a coefficient of variation of ~6%
(457–732 mm). Mapping that linearly to width made the middle 50% of years occupy
just 0.236 of the available width range — visually no variation at all. Ranking
doubles that to 0.496.

**Cost, and it is real:** width is now relative *and* non-linear. A wide ring
means "wetter than most years", not "wetter by this many mm". Precedent:
dendrochronology works in ring-width *index*, not millimetres. The figure says so.

**To undo:** use a linear min–max stretch in `ring_widths()`. Expect the rainfall
signal to become invisible again.

### The rank is pooled across models, not per-model

Pooling means a systematically drier model keeps systematically thinner rings
(mean rank 0.39 / 0.50 / 0.61 for NorESM2-MM / TaiESM1 / UKESM1-0-LL). Ranking
per-model would stretch every model to fill the full range and destroy exactly
the comparison a multi-model grid exists to make.

### Temperature is smoothed with a GAM

You suggested a GAM; the decision to smooth at all was mine. Both smoothers agree
closely (correlation 0.98–0.99, max difference 0.64 °C), but the GAM is monotone
rising in 100% of years for two of three models, against 66% and 73% for an
11-year running mean. For a strip whose only job is showing an increase, no local
reversals is strictly better.

### Every channel's scale is shared across runs

`shared_scales()` computes normalisation constants **once** across every model,
replicate and landscape, and that one object is passed into each disc. Without it,
a model warming +2.0 °C renders with the same colour spread as one warming
+6.2 °C, and a grid becomes silently incomparable while looking fine.

The temperature bounds come from the **smoothed** series, because that is what is
drawn. Taking them from raw annual values left 35% of the colour ramp unreachable.

It is called `shared_scales`, not `climate_scales`, for a reason: under the old
name **disturbance was left out for several rounds.** Both ring width and scar
size computed `di / max(di)` from the single run they held, so every disc rendered
its own worst year at full strength. Across the 108 runs, per-run maxima span
0.0154 to 0.2038 — a **13× spread** — so two discs showing identical narrowing
could differ thirteenfold in what actually happened. One object, one answer to
"where does my scale come from".

### Disturbance is compressed by `dist_gamma`, not left linear

A shared scale on its own is not enough. At the original exponent of 0.7 the
quietest of the 108 runs narrowed its worst ring by only **8%** — present in the
data, invisible on the page.

`dist_gamma = 0.35` lifts small relative values while keeping the order:

| | quiet run | loud run | fidelity |
|---|---|---|---|
| per-run (the bug) | 45% | 45% | none |
| shared, γ=0.7 | 8% | 45% | 1.00 |
| **shared, γ=0.35** | **17%** | **45%** | **0.976** |
| pooled rank | 41% | 45% | 0.59 |

"Fidelity" is the correlation between a run's true peak disturbance and the
narrowing it shows. Pooled rank — the treatment rainfall gets — was rejected here
because it flattens a 13× magnitude difference to almost nothing, which is the
opposite failure to the one being fixed.

`dist_amp()` is the single place this happens, feeding both ring width and scar
size, so the two can never tell different stories about the same year. It raises
an error if handed a scales object without `di_max`, because `di / NULL` silently
returns length zero — the exact class of bug that hid the original problem.

### The colour ramp is blue→magenta, not viridis

Viridis is the clearer ramp in isolation, and was tried. Its green sits ΔE 1.7
from the black spruce colour — effectively identical, in a strip drawn beside
spruce wedges. Cividis (4.6) and plasma (2.2) collide too, on their yellow and
orange ends. Blues through magentas are the one region the species palette leaves
free, and the ramp used gets ΔE 20.4 clear while keeping viridis's actual
mechanism, monotone lightness.

### Scars sit at the species that lost the most share

Nothing in the data says *where* a disturbance happened — it is a landscape
aggregate. Rather than an arbitrary angle, the crack is placed at the wedge of the
species that lost most that year, so **the angle says which species was hit, not
where.** The figure states this.

### A hairline at every year boundary

This is what makes it read as a core rather than a wobbly pie chart: in a real
core the latewood band is what you count, and tight line spacing *is* the signal.
Without it the disc looked like a decorated pie chart.

### The species palette was rebuilt for colour-vision deficiency

The inherited esa-2026 palette solved dark-mode *contrast* but had never been
checked for CVD. Under Machado 2009 at full severity its worst all-pairs
separation was **ΔE 3.8** (white spruce vs mixed spruce) — effectively one
colour; mixed spruce vs mixed deciduous was 4.6. A deuteranope saw roughly three
groups, not six.

Designing it by reasoning failed first: an even lightness ladder scored **0.9**,
*worse* than the legacy palette, because CVD simulation redistributes lightness
according to chroma. So the space was searched instead. Six colours can reach
ΔE 21.7 unconstrained and **13.9 while keeping species-evocative hues**, so
semantics cost almost nothing. The palette is the hill-climbed result inside
those hue bands:

| | | |
|---|---|---|
| all-pairs CVD ΔE | **9.8** | was 3.8 |
| all-pairs normal ΔE | 15.6 | ≥ 15 |
| min contrast on `#0d0d0d` | 3.3 | ≥ 3 |
| nearest the temperature strip | 12.1 | ≥ 12 |
| max chroma | 0.151 | off-neon |

Lightness is assigned by **abundance**, not just semantics — the two smallest
categories are the brightest, because a 2% wedge has to be findable:
aspen (12%) gold at 13.1:1, birch (2%) coral at 8.6:1, while mixed deciduous
(21%) recedes to a dark brown at 3.3:1. Putting the big warm category dark and
the tiny one bright is also what buys those two the lightness separation they
need to survive CVD.

Birch as coral rather than white is the one real semantic compromise: gold went
to aspen, which owns that association more strongly.

### Two small colour adjustments

- A gap between species wedges, because adjacent wedges physically touch and the
  closest adjacent pair is ΔE 9.8 under CVD.
- Text lifted one step from the fire animation's values (body text 5.4:1 → 7.5:1
  contrast). Light grey on near-black is harder to read than its ratio suggests.

---

## 5. What the figure cannot tell you

- **The index is not fire.** It is species change from any cause — fire,
  succession, mortality, establishment. A slow successional shift and a
  stand-replacing burn look identical if the same area crosses a category boundary.
- **It only counts changes that cross a category boundary.** If black spruce
  thins severely but stays dominant, the index registers nothing. Partial
  disturbance is under-counted.
- **It measures net change, not turnover.** A year where 3% moves spruce→
  deciduous and 3% moves back scores zero, though 6% of the landscape turned over.
  Treat it as a lower bound.
- **No spatial information.** Landscape aggregate only.
- **The lobed outline is texture, not a measurement.** Rendered rim dent
  correlates only **0.28** with a run's peak disturbance and 0.37 with its
  cumulative total, and narrowing the wounds does not improve it (0.19–0.28 at
  every `scar_sigma` tried). The cause is timing: at `scar_decay = 0.99` a
  mid-century wound accumulates persistence over 40-odd rings while a 2090s wound
  of the same size gets a few. Authentic for a core, useless as a readout. The two
  extremes do look visibly different, but you cannot rank two mid-range discs by
  their silhouette. **Ring width is the quantitative disturbance channel**, and the
  figure's caption claims nothing about the outline.
- **Climate is regional, not per-landscape.** `process_gcms.R` takes a spatial
  mean over the whole downscaled raster, so all three landscapes share one
  climate series.
- **One replicate per disc, deliberately.** Replicates are genuinely stochastic;
  averaging them would hide the variability that is part of the story.

---

## 6. Running it

```r
source("R/iland_rings.R")     # renders NorEsm2-MM and TaiESM1 for landscape 01
```

```r
# a single model
options(iland_rings_models = "UKESM1-0-LL"); source("R/iland_rings.R")

# define functions without rendering
options(iland_rings_no_run = TRUE); source("R/iland_rings.R")
```

### Triptychs

Three discs in a triangle sharing one legend, for comparing along a single axis.
Additive — `render_disc_frame()` is untouched and stays the primary output, so
triptychs are opt-in:

```r
options(iland_rings_triptych = c("model", "replicate", "landscape"))
source("R/iland_rings.R")

# or directly
render_triptych("render/x.png", vary = "landscape", model = "TaiESM1")
```

`vary` is one of `model`, `replicate`, `landscape`. **Not SSP** — `ssp370` exists
only in the `onlyfire` runs, so the vegetation output is ssp245 only.

A 3 × 3 grid was rejected on measurement. At 1920 px square a cell disc has a
301 px radius, putting the mean ring at 3.2 px and the thinnest at **1.1 px** —
below the ~2 px floor where antialiasing dissolves them, and the 87 hairlines
merge into grey. Matching today's single-disc quality in a 3 × 3 needs a 3840 px
canvas. Three discs reach it at 3000 px (~2.5 MB, ~10 s), and three panels is
also about as much as one figure can be read at once.

One caveat for `vary = "landscape"`: all three temperature strips will be
identical, because the climate series is regional rather than landscape-specific
(see §2). It is correct, but the strip carries no information in that variant.

Knobs worth touching, all at the top of `R/iland_rings.R`:

| knob | now | effect |
|---|---|---|
| `w_clim_min` | 0.20 | rainfall drama. **Near its floor** — the thinnest ring is ~2.2 px and below ~2 px rings dissolve under antialiasing |
| `dist_pen` | 0.45 | how much disturbance narrows a ring; competes with rainfall for the same channel |
| `dist_gamma` | 0.35 | how much a quiet run still shows on the shared scale. Raising it toward 1 is more faithful to magnitude and less visible; see §4 |
| `scar_top_n` | 8 | how many years get a crack — per run, so every disc always has 8 |
| `scar_depth_rings` | 9.0 | how far scars dent the outline. Calibrated against the current `di_max`; see §7 |
| `lesion_halfw` | 0.0090 | crack thickness |
| `disc_split` / `scale` | 0.64 / 0.49 | disc size, which is what buys pixels per ring |

---

## 7. Known open issues

None blocking. The per-run disturbance normalisation that used to live here is
fixed — see §4.

Two things to be aware of before building the grid:

- **`scar_depth_rings` is calibrated against the current pool.** It was raised
  6.0 → 9.0 when disturbance moved to the shared scale, tuned so the reference
  disc (L01 / TaiESM1 / rep 1) kept roughly the dent it had. If runs are added
  whose disturbance exceeds the current `di_max` of 0.2038, every existing disc
  quietly flattens, because amplitudes are relative to that maximum. Re-check the
  three reference discs after any new preprocessing run.
- **`w_clim_min` is at its floor.** The thinnest ring is ~2.2 px and below about
  2 px rings dissolve under antialiasing. More rainfall drama needs a bigger
  canvas, not a smaller `w_clim_min`.

---

## 8. Notes for future edits

- **Do not blanket-regex these files from PowerShell.** `Set-Content -Encoding
  UTF8` on Windows PowerShell 5.1 writes a byte-order mark that R cannot parse,
  and repairing it by byte round-trip double-encodes every `°`, `·` and em dash.
  Use targeted edits.
- **`draw_locator()` is aspect-aware; do not assume square units.** It derives
  separate x and y scales from `par("pin")` and `par("usr")` of the *current*
  region. Disc regions use `asp = 1` so both come out identical and the
  correction is a no-op, but the legend region has no `asp` and its units are
  6.667 vs 7.407 inches per unit — a 1.11× vertical stretch that rendered Alaska
  11% too tall until it was fixed. Any new region hosting the locator inherits
  the fix automatically; anything else drawn in a non-`asp` region does not.
- **Square regions clip negative coordinates.** The triptych's disc regions are
  square, so `asp = 1` gives no vertical expansion and labels placed below the
  disc vanish. The single-disc frame's region is taller than wide, which is why
  the same code worked there. Use the corners of the square — they are free
  around an inscribed circle.
- `R/common.R` is the **Alaska fire** piece's module, not this one. It sets
  `first_year <- 1980` / `last_year <- 2020`, which would clash with this piece's
  2014–2100. `iland_common.R` re-declares the few loaders it needs on purpose.
- Species colours exist in three places: the iLand app, `esa-2026/R/theme.R`, and
  here. **esa-2026 is canonical** — it is the only one that solves dark mode, and
  it does so by measurement.
