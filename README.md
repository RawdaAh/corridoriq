# AAM Corridor Lab

Code and simulation for a research program on **shared low-altitude Advanced Air Mobility (AAM) corridors**:

1. **Wind-aware routing** — *Highways in the Sky: UAV Routing under Dynamic Wind Profiles in Low-Altitude Urban Airspace* (TRR).
   Time-dependent Pareto label-correcting search on a 4D corridor network with wind projected onto every link.
2. **Shared-intent reservation (SIRE)** — *Shared Advanced Air Mobility Corridors: Configuration, Capacity, and Demand Management under a Reservation Protocol* (TRBAM-27-05211).
   Sequential limited-disclosure registry of link-time reservations for competing operators.
3. **Rolling-horizon scheduling** — priority-aware rolling-horizon corridor scheduling with full-mission commitments (in preparation).
   Two-stage priority ILP over space-time plans vs. FCFS vs. an offline benchmark.
4. **Bridge: value of wind information** — every mission is planned with four planning winds and flown in the true wind. It measures
   the energy lost to imperfect wind knowledge (Paper 1) and the **reservation drift** that a SIRE registry would have to absorb (Paper 2).

All four run on **one shared core**: the same modular corridor lattice, wind field, energy model and link-time reservation logic.

5. **Real city (downtown Dallas)** — the corridor is laid over OpenStreetMap buildings between Dallas Farmers Market and the
   American Airlines Center (2.46 km, 30–120 m). Every node and link that would hit a building (8 m clearance) is removed, a
   building-aware wind is computed around the real footprints, and missions are planned with the four planning winds and flown in
   the true wind, exactly as in GeoWind2Plan; many-operator SIRE traffic runs on the same city corridor.

**Interactive demos**
* `demo/city.html` — 3D Dallas: buildings, roads, your corridor network, wind streamlines, quadcopters flying each planned route,
  Start/Goal markers, wind dial, random tail/head/crosswind missions, energy/altitude/ground-speed/wind-on-track charts, and a
  corridor-traffic mode running the SIRE registry live in the browser.
* `demo/index.html` — replays the three paper experiments (abstract lattice, space-time ledgers, request boards).

## Real-city layer

```bash
python run_geo_mission.py --direction 165 --speed 4          # one headwind mission, printed
python run_geo_mission.py --missions 48 --speeds 2 4 8         # GeoWind2Plan-style table by wind class
python run_geo_mission.py --fetch --bbox <min_lon min_lat max_lon max_lat> --start <lat lon> --goal <lat lon> --osm data/osm/mycity.osm
python make_geo_demo.py                                         # rebuild demo/city.html
```
* `aamsim/geo.py` — OSM parsing (heights from `height`, `building:levels` × 3.5 m, or a type default; building parts), corridor
  placement along a start→goal axis, rasterised building heights, removal of blocked nodes/links, launch pads.
* `aamsim/urbanwind.py` — building-aware wind surrogate: per-altitude potential flow around the solid footprints (red-black SOR),
  a wake deficit behind buildings and a power-law shear profile. It reproduces the qualitative structure (canyon speed-up,
  wakes, turning) in about a second; it is **not CFD**. Replace it with CFD or LES output through `aamsim.wind.wind_from_arrays`
  for a paper result.
* `aamsim/geo_plan.py` — the four planning winds (ground truth; a coarse weather-model forecast with a 15° direction error and
  −15% speed bias and no buildings; the height-only profile; zero wind), energy-optimal planning, evaluation in the truth.
  Crosswinds cost energy through a crab angle (ground speed = √(v² − w⊥²) + w∥), an extension of the TRR projection.
* `demo/geo_engine.js` — the same wind solver, planner and SIRE registry in JavaScript, so the page recomputes live when the
  wind dial moves. `tests/test_geo.py` checks the browser engine against Python (node wind to 1e-12, identical routes up to
  exact ties, identical SIRE accept/reject decisions).

Map data © OpenStreetMap contributors, available under the Open Database License (ODbL).

---

## Installation

```bash
pip install -r requirements.txt
```

Python ≥ 3.10. `numba` is optional; without it SIRE falls back to a vectorized numpy search, which is slower but gives identical results (this is tested).

## Code

```
run_wind_routing.py       Paper 1 experiments (TRR Sec. 6.2-6.5)
run_sire_experiment.py    Paper 2 experiments (TRB Figs. 1-4) + equity and wind extensions
run_rolling_horizon.py    Paper 3 design grid (FCFS / rolling ILP / offline)
run_wind_belief.py        Bridge experiment: plan with a belief, fly in the truth
summarize_results.py      prints tables for any result CSV
make_demo.py              builds demo/index.html (paper replays)
make_geo_demo.py          builds demo/city.html (3D Dallas: OSM buildings, wind, drones, SIRE traffic)
run_geo_mission.py        real-city wind-information experiment (OSM buildings)
aamsim/
  corridor.py             modular prism lattice (edges, face and body diagonals), portals
  wind.py                 wind fields, TRR link projection (Eqs. 14-25), synthetic truth, planning beliefs
  energy.py               vehicle model, baseline and phase-dependent energy (Eqs. 10-13, 26-33)
  pareto.py               time-dependent bi-objective label-correcting search (Algorithms 1-2)
  sire.py                 SIRE registry, operator route search on (node, slot), KPIs (SR, Gini, cascade)
  rolling.py              port of the priority rolling-horizon ILP notebook, with independent audit
  belief.py               four planning winds, mission sampling, evaluation in the truth
  plotting.py             paper figure style (colour-blind-validated palette)
  geo.py                  OpenStreetMap buildings, corridor over a real city
  urbanwind.py            building-aware wind surrogate (potential flow + wakes + shear)
  geo_plan.py             four planning winds, energy-optimal routes, evaluation in the truth
demo/geo_engine.js        browser twin of urbanwind / geo_plan / sire
data/osm/dallas.osm.gz       OpenStreetMap extract of downtown Dallas (ODbL)
tests/test_core.py        geometry, reproduction of published values, reservation exclusivity, audit
demo/template.html        viewer; make_demo.py injects the data and writes demo/index.html
```

## Usage

### Paper 1 — wind-aware routing

```bash
python run_wind_routing.py --experiment all
python run_wind_routing.py --experiment orientation --phase-dependent   # climb + take-off energy terms
python run_wind_routing.py --experiment congestion --lanes 4 --reps 10
```
Output: `results/wind_routing/{orientation,spatial,temporal,congestion}.csv` and `figures/`. Each row reports the number of Pareto paths,
mean travel time and energy, ascents/descents, and the share of time and energy by altitude band and by longitudinal block (Tables 1–5).

### Paper 2 — SIRE

```bash
python run_sire_experiment.py --experiment demand    --reps 30 --workers 8
python run_sire_experiment.py --experiment geometry  --reps 30 --workers 8
python run_sire_experiment.py --experiment all       --reps 5            # quick look
```
| experiment | varies | paper |
|---|---|---|
| `demand` | 100 … 4,000 trips | Fig. 1 |
| `geometry` | lateral modules N_W 1–4 × altitude modules N_H 1–4 at 2,000 trips | Fig. 2 |
| `direction` | forward/reverse split × shared (undirected) vs direction-separated key | Fig. 3 |
| `dispatch` | dispatch window 5–60 min × uniform vs peak | Fig. 4 |
| `equity` | same-minute tie-break: operator id (published) vs lottery vs rotating | new |
| `wind` | link times from the TRR shear wind (tail vs head) | new |

Replicates use common random numbers across cells (`seed = 1000 + rep`). A 2,000-trip replicate takes about 1–2 s with numba.

### Paper 3 — rolling horizon

```bash
python run_rolling_horizon.py --design base  --seeds 3
python run_rolling_horizon.py --design quick --seeds 5  --workers 4
python run_rolling_horizon.py --design full  --seeds 30 --workers 8     # the notebook's full grid
```
Every schedule is rebuilt from geometry and audited (no shared resource token, request windows, freeze rule) before metrics are saved.
`summary.csv` holds paired ILP − FCFS and offline − ILP differences with bootstrap 95% CIs over seeds. Large loads (≥ 200 requests) solve many MILPs
with SciPy/HiGHS and can take a long time per seed, as in the original notebook.

### Bridge — value of wind information

```bash
python run_wind_belief.py --missions 256 --wind-speed 15 --workers 8
python run_wind_belief.py --mission 0                       # one mission, printed
python run_wind_belief.py --wind-speed 8 --select time       # pick the fastest Pareto route instead of the cheapest
python run_wind_belief.py --model continuous                 # crosswind and oblique components act on every link
```
| method | planning wind |
|---|---|
| Ground truth | the true field (perfect information reference) |
| Forecast | truth with a speed bias (±15%), a direction error (σ = 15°), smooth spatial error and a 6-min lag |
| Profile wind | one sounding at departure, averaged over each altitude layer and applied uniformly |
| Wind-agnostic | zero wind |

Output: `results/wind_belief/missions_<speed>mph_<rule>.csv` (one row per mission and method) and `summary_<...>.csv`, a table in the
style of GeoWind2Plan: Wh/km (mean, sd), extra energy vs ground truth, extra time, absolute arrival error, the reservation slot shift
and the share of links that leave their reserved slot, mean altitude and ground speed, reported for all, tailwind, headwind and crosswind missions.

### Demo

```bash
python make_demo.py              # runs small scenarios of every study (~2-3 min), writes demo/index.html
python make_demo.py --from-json  # rebuild the HTML only
```
Open `demo/index.html` in a browser. Tabs: wind routing, SIRE reservations, rolling horizon, wind information.

### Tests

```bash
python -m pytest -q tests
```

---

## Reproduction of published results

Values produced by this repository against the published or submitted numbers.

**Paper 1 (TRR), Table 1, altitude-varying wind**

| wind | time (min) repo / paper | energy (kWh) repo / paper |
|---|---|---|
| none | 9.74 / 9.74 | 0.1702 / 0.170 |
| 0° tail | 5.50 / 5.50 | 0.0539 / 0.054 |
| 45° | 6.32 / 6.30 | 0.0705 / 0.070 |
| 135° | 10.03 / 9.88 | 0.1772 / 0.174 |
| 180° head | 16.29 / 16.29 | 0.4775 / 0.477 |
| 225° | 9.95 / 10.04 | 0.1738 / 0.175 |
| phase-dependent: none / tail / head | 9.74 / 5.50 / 16.29 (same) | 0.187 / 0.071 / 0.495 vs 0.187 / 0.071 / 0.495 |

Spatial blocks (Table 3): tail 5.95 min / 0.0635 kWh with 2 Pareto paths (paper 5.95 / 0.063 / 2); 45° 6.73 / 0.080 (paper 6.71 / 0.079).
The 180° spatial case gives 15.86 min / 0.454 kWh against 17.51 / 0.559 in the paper; the block-to-node assignment at block boundaries
is the likely difference and is set in `aamsim.wind.trr_spatial_wind`. The congestion experiment (6.5) depends on the lane count, which
the paper does not state; use `--lanes`.

**Paper 2 (TRB), Figure 1** (10 replications; paper uses 30)

| demand | SR repo / paper | accepted repo / paper | Gini repo / paper |
|---|---|---|---|
| 200 | 99.8% / 99.57% | | |
| 600 | 89.5% / 89.68% | | |
| 1,000 | 69.9% / 69.97% | | 0.183 / 0.183 |
| 2,000 | 36.6% / 36.86% | 731 / 737 | 0.477 / 0.478 |
| 4,000 | 17.5% / 17.68% | 699 / 707 | 0.638 / 0.636 |

Accepted-trip travel time 585 s at low demand and 2.2% inflation at 2,000 trips (paper 585 s, 2.30%). Figure 3, 50/50 split: shared
44.3% vs separated 72.5% (paper 44.20% / 72.17%). The notebook's reproduction reference (D = 2,000, replicate 1, SR = 0.3645) is
reproduced exactly and is a unit test.

**Paper 3** — the 50-request pilot of the notebook is reproduced (ILP and offline mean delay 0.016 min, max 0.2 min; all requests served).

## Modelling notes and limits

* Link travel time and energy follow the TRR equations: the nodal wind is averaged on the link and projected on the link heading. Under the
  paper's rule, winds more than 45° off the track do not change ground speed; `model="continuous"` uses the full along-track component.
* SIRE routes depart exactly at the intended slot (no airborne waiting), keyed on the undirected physical link (or the directed arc for
  the separated rule); a link held for τ seconds blocks `max(1, ceil(τ/Δt))` slots.
* Routes in Papers 1 and 3 may not land at an intermediate ground node. SIRE keeps the TRB notebook's rule (only ground-to-ground links removed).
* The "ground truth" wind in the bridge experiment is **synthetic** (power-law shear, smooth building-wake bumps, veer, drift). It is a
  controlled stand-in so the planning beliefs have structure to get wrong. For a paper result, replace it with interpolated
  observations or model output (e.g. HRRR analysis vs. forecast, or CFD) through `aamsim.wind.wind_from_arrays`.
* Rolling horizon: `rolling.py` implements the 14 Oct 2026 notebook (priority ILP). The earlier draft's learned afterstate critic and Lagrangian
  decomposition are not implemented here.

## Citation

If you use this code, please cite the corresponding papers (Ahmed and Abdelghany, Southern Methodist University).

## License

MIT
