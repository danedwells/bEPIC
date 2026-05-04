# bEPIC — Bayesian EPIC Earthquake Early Warning Location Algorithm

bEPIC is a Bayesian extension of the ShakeAlert EPIC (Earthquake Point-source Integrated Code) system. It processes EPIC system log files to produce iteratively updated posterior estimates of earthquake location and magnitude, adding one station trigger at a time (each update is a new "version").

---

## Table of Contents

1. [Overview](#1-overview)
2. [Repository Structure](#2-repository-structure)
3. [Environment Setup](#3-environment-setup)
4. [The .run File Format](#4-the-run-file-format)
5. [Per-Event Directory Structure](#5-per-event-directory-structure)
6. [Core Pipeline — Step by Step](#6-core-pipeline--step-by-step)
7. [Module Reference: src/bEPIC/](#7-module-reference-srcbepic)
   - [bEPIC_main.py](#bepic_mainpy)
   - [locate.py](#locatepy)
   - [likelihood.py](#likelihoodpy)
   - [prior.py](#priorpy)
   - [posterior.py](#posteriorpy)
   - [magnitude.py](#magnitudepy)
   - [geospatial_util.py](#geospatial_utilpy)
   - [data_util.py](#data_utilpy)
   - [bEPIC_analysis.py](#bepic_analysispy)
   - [pull_logs.py](#pull_logspy)
   - [generate_synthetic_run_file.py](#generate_synthetic_run_filepy)
8. [Key Parameters and Constants](#8-key-parameters-and-constants)
9. [Output Files](#9-output-files)
10. [Testing](#10-testing)
11. [Zextra — EPIC C Locator Tools](#11-zextra--epic-c-locator-tools)
    - [EPIC_locate_prelim.py](#epic_locate_prelimpy)
    - [example_call_EPIC_locate.py](#example_call_epic_locatepy)
    - [download_usgs_event.py](#download_usgs_eventpy)
    - [plot_prior.py](#plot_priorpy)
12. [Catalog and Testing Data](#12-catalog-and-testing-data)
13. [Known Issues and Notes](#13-known-issues-and-notes)

---

## 1. Overview

### What it does

For a single seismic event, bEPIC reads EPIC's real-time output (station trigger times and peak displacements) and runs a Bayesian grid search. At each version (each new station trigger) it computes:

1. **Likelihood** — how well does each point on a spatial grid explain the observed P-wave arrival times?
2. **Prior** — how historically seismically active is each grid point, based on ANSS catalog?
3. **Posterior** — the normalized product; the MAP (maximum a posteriori) point is the estimated epicenter.
4. **Magnitude** — estimated from peak ground displacement (Pd) at each station using a fixed scaling relation.

### Relationship to EPIC

The real-time ShakeAlert EPIC system is written in C/C++. bEPIC is a Python Bayesian wrapper around the same inputs. The `Zextra/` directory also contains `EPIC_locate_prelim.py`, a direct Python transliteration of the EPIC C `searchGrid` routine, for comparison and testing.

---

## 2. Repository Structure

```
bEPIC/
├── src/bEPIC/                   # Core package (importable as `from bEPIC import ...`)
│   ├── __init__.py
│   ├── bEPIC_main.py            # Entry points: initialize_bEPIC_event(), run_bEPIC()
│   ├── locate.py                # 2-station initial grid center estimate
│   ├── likelihood.py            # Travel-time misfit likelihood surface
│   ├── prior.py                 # KDE seismicity prior; catalog download
│   ├── posterior.py             # Posterior = normalize(prior × likelihood)
│   ├── magnitude.py             # Pd-based magnitude scaling (Chung et al. 2020)
│   ├── geospatial_util.py       # Coordinate transforms, grid creation, distances
│   ├── data_util.py             # Log parsing, .run generation, USGS queries, TT table
│   ├── bEPIC_analysis.py        # Post-hoc trigger time misfit vs. USGS catalog
│   ├── pull_logs.py             # CLI script: download + parse EPIC logs from archive
│   └── generate_synthetic_run_file.py  # Generates a synthetic .run for testing
│
├── data/                        # Static data files (non-Python)
│   ├── h2p+ak135.080            # P-wave travel time table (distance km → time s)
│   └── prior_seismicity_catalog.txt  # ANSS catalog (M≥3, western US, 2000–present)
│
├── example/
│   └── run_example_event_01.py  # Minimal example: run event 126625
│
├── tests/
│   ├── test_locate.py           # pytest: get_two_station_location()
│   └── test_geospatial_util.py  # pytest: LL2cartd, ckm2LLd, get_dist_between_two_points_km
│
├── Zextra/                      # Supplementary scripts (EPIC C replication + USGS download)
│   ├── EPIC_locate_prelim.py    # Python transliteration of EPIC C searchGrid
│   ├── example_call_EPIC_locate.py   # Example usage of EPIC_locate_prelim
│   ├── download_usgs_event.py   # Download USGS event + phases → EPIC inputs
│   ├── plot_prior.py            # Plot prior_seis_grid_US_Canada.tt3
│   ├── prior_seis_grid_US_Canada.tt3  # Pre-built prior grid (US + Canada)
│   ├── phases.csv               # Example USGS phase download (Central America event)
│   ├── bEPIC_testing_catalog.txt     # EPIC vs. ANSS comparison catalog (~740 events)
│   └── run_files/               # ~600+ real .run files for batch testing
│
├── generate_synthetic_run_file.py    # Root-level copy: writes to bEPIC_test/999999/
├── test_1.py                    # Ad-hoc dev script (uses torch, postgres_id='9999')
├── pyproject.toml               # Package metadata (pip install -e .)
└── CLAUDE.md                    # Claude Code instructions for this repo
```

---

## 3. Environment Setup

### Required environment variable

```bash
export BEPIC=/path/to/bEPIC
```

Used by `data_util.travel_time_function()` and `pull_logs.py` to locate data files. Without it, those functions will fail with a `KeyError`.

> **Note:** `prior.py` and `bEPIC_main.py` currently have `bepic` hardcoded to `/home/daned/2024_NEHRP/bEPIC` rather than reading `$BEPIC`. This is a known inconsistency — the env var path is preferred.

### Installation

```bash
cd /path/to/bEPIC
pip install -e .
```

After this, `from bEPIC import bEPIC_main` works from anywhere.

### Dependencies

| Package | Usage |
|---------|-------|
| `numpy` | All numerical operations |
| `pandas` | CSV/TSV I/O, all tabular data |
| `scipy` | `interp1d` (travel times), `gaussian_kde` (prior) |
| `obspy` | `gps2dist_azimuth` (distances) |
| `requests` | USGS ComCat API queries |
| `libcomcat` | `prior.generate_prior_seismicity_catalog()` only |
| `torch` | Imported in `test_1.py` only (not in core pipeline) |
| `matplotlib` | `Zextra/plot_prior.py` only |
| `pytest` | Running unit tests |

---

## 4. The .run File Format

The `.run` file is the primary input to `run_bEPIC()`. It is a CSV with one row per (version, station) pair.

```
version,order,station,channel,network,location,longitude,latitude,trigger time,tterr,logPd
0,1,SAO,HNZ,BK,00,-121.4472,36.764,1538771380.09,-0.053,-1.89848
0,2,SAO,HHZ,BK,00,-121.4472,36.764,1538771380.11,-0.073,-1.875756
0,3,BSR,HNZ,NC,--,-121.5203,36.6674,1538771381.18,-0.035,-2.964038
1,1,SAO,HNZ,BK,00,-121.4472,36.764,1538771380.09,-0.037,-1.89848
1,6,1865,HNZ,NP,--,-121.3908,36.8831,1538771381.27,-0.080,-1.830575
...
```

| Column | Type | Description |
|--------|------|-------------|
| `version` | int | EPIC update index. Version 0 = first alert (2 stations). Each new station trigger is a new version. |
| `order` | int | Station's trigger order within this version (1 = first to trigger) |
| `station` | str | Station code |
| `channel` | str | Channel code (e.g. HHZ, HNZ) |
| `network` | str | FDSN network code |
| `location` | str | Location code (e.g. 00, --, blank) |
| `longitude` | float | Station longitude (degrees) |
| `latitude` | float | Station latitude (degrees) |
| `trigger time` | float | Unix timestamp of P-wave trigger at this station |
| `tterr` | float | Travel time residual from EPIC's location. Rows with `tterr <= -999` are flagged as bad picks and excluded from location. |
| `logPd` | float | log10 of peak displacement (Pd) in meters, used for magnitude estimation |

**Important:** The same station may appear in multiple versions with slightly different `tterr` values as EPIC refines its location. The station count at any version = the number of unique stations in that version's rows with `tterr > -999`.

---

## 5. Per-Event Directory Structure

```
{project_parent_directory}/{postgres_id}/
    {postgres_id}.run                          # Primary input (CSV)
    EPIC/
        {postgres_id}_event_summary_log.txt    # Per-version EPIC location + metadata (TSV)
        {postgres_id}_event_triggers_log.txt   # Per-version station trigger details (TSV)
        {postgres_id}_location_triggers_log.txt
        {postgres_id}_station_summary_log.txt
        {postgres_id}_station_counts_log.txt
        {postgres_id}_misc_log.txt
        {postgres_id}_epic_location_log.txt    # EPIC's own location output (TSV)
    bEPIC/
        {postgres_id}_bEPIC_log.txt            # bEPIC output: per-version location + magnitude (TSV)
    USGS/
        usgs_event_summary.txt                 # Matched USGS ComCat event (TSV)
        usgs_trigger_time_misfit.txt           # Station trigger residuals vs. USGS location (TSV)
    plots/                                     # (empty; reserved for output plots)
```

`postgres_id` is always a 6-character zero-padded string (e.g., integer `126625` → `'126625'`). This zero-padding is applied automatically inside all functions that accept an integer.

---

## 6. Core Pipeline — Step by Step

```
Raw EPIC archive log (.log.gz)
    └─► pull_logs.py                  # download + decompress
        └─► data_util.parse_log()     # extract 7 TSV files into EPIC/
            └─► data_util.generate_run_file()  # build {id}.run from event_triggers_log
                └─► run_bEPIC()
                    ├── version 0: locate.get_two_station_location() → CenterPoint
                    ├── for each version:
                    │   ├── likelihood.calculate_likelihood()  → like (p×p), misfit (p×p)
                    │   ├── prior.compute_prior()              → prior (p×p)
                    │   ├── posterior.compute_posterior()      → post (p×p), MAP lon/lat
                    │   └── magnitude.compute_magnitude()      → M (float)
                    └── save → bEPIC/{id}_bEPIC_log.txt
```

**Version skipping:** If the number of stations in version N equals the number in version N-1, the locate/likelihood/prior/posterior computations are skipped and the previous version's locations are reused. Only magnitude is recomputed. This happens because some EPIC versions update only travel time residuals or magnitude estimates without a new station triggering.

---

## 7. Module Reference: src/bEPIC/

### bEPIC_main.py

Top-level entry point. Two public functions:

---

#### `initialize_bEPIC_event(project_parent_directory, postgres_id)`

Sets up the directory structure and input files for an event. Run this before `run_bEPIC()`.

**Actions:**
1. Checks that `{project_parent_directory}/{postgres_id}/` exists (exits with a print if not).
2. Creates subdirectories `EPIC/`, `plots/`, `bEPIC/`, `USGS/` if missing.
3. Moves any EPIC log files found loose in the event directory into `EPIC/`.
4. If no `.run` file exists, calls `data_util.generate_run_file()`.
5. If no `USGS/usgs_event_summary.txt` exists, calls `data_util.search_for_USGS_event()` (wrapped in try/except; failure is non-fatal).

---

#### `run_bEPIC(project_parent_directory, postgres_id, velocity_model, GridSize, GridSpacing)`

Main Bayesian location loop.

**Args:**
- `project_parent_directory` (str): path ending in `/`
- `postgres_id` (str or int): auto zero-padded to 6 digits
- `velocity_model` (str): `'h2p+ak135'` or `'constant'`
- `GridSize` (float): grid half-width in km (default usage: 200 km)
- `GridSpacing` (float): grid node spacing in km (default usage: 2 km)

**Logic:**
1. Loads the `.run` file; adds a `sigma` column (all ones — uniform uncertainty weights).
2. Loops over versions 0 through max.
3. Filters rows by `version == v` and `tterr > -999`.
4. Version 0: calls `locate.get_two_station_location()` for the initial grid center.
5. Versions > 0: uses previous posterior lon/lat as the new grid center.
6. If station count unchanged from the previous version: skips locate/likelihood/prior/posterior, reuses previous locations.
7. Calls likelihood, prior, posterior, then magnitude (at both likelihood MAP and posterior MAP).
8. Appends one row per version to `bEPIC_df`.
9. Saves `bEPIC_df` to `{postgres_id}/bEPIC/{postgres_id}_bEPIC_log.txt` as TSV.

**Output columns of bEPIC_log.txt:**

| Column | Description |
|--------|-------------|
| `version` | Version index |
| `num stations` | Number of stations used |
| `likelihood lon/lat` | MAP of likelihood surface |
| `likelihood mag` | Magnitude at likelihood MAP |
| `prior lon/lat` | MAP of prior surface |
| `posterior lon/lat` | MAP of posterior (best estimate) |
| `posterior mag` | Magnitude at posterior MAP |

---

### locate.py

#### `get_two_station_location(sta_df) → [lon, lat]`

Estimates the initial grid center from the first two stations (order == 1 and order == 2).

**Method:** Weighted 1/3 – 2/3 interpolation between the two station locations. The earlier-triggering station gets 2/3 weight (the earthquake is closer to the station that triggered first).

```
If station 1 triggers first:
    eq_lat = (2*lat1 + lat2) / 3
    eq_lon = (2*lon1 + lon2) / 3

If station 2 triggers first (or tie goes to station 1 branch):
    eq_lat = (2*lat2 + lat1) / 3
    eq_lon = (2*lon2 + lon1) / 3
```

Returns `[eq_lon, eq_lat]`.

---

### likelihood.py

#### `calculate_likelihood(CenterPoint, sta_df, velocity_model, GridSize=200, GridSpacing=2)`

Computes the data-driven likelihood surface L(location | data) over a 2D spatial grid.

**Method:** Gaussian likelihood from P-wave travel time residuals.

**Step-by-step:**

1. Builds a square grid centered on `CenterPoint` via `geospatial_util.make_grid()`. The grid has `(2*GridSize/GridSpacing + 1)` nodes per side.

2. Converts station lon/lat to local Cartesian coordinates (km) relative to `CenterPoint` using `geospatial_util.LL2cartd()`.

3. Computes slant distance from every grid node to every station (including a fixed depth of **8.0 km**):
   ```
   distance[node, station] = sqrt(Δx² + Δy² + 8²)   [km]
   ```

4. Computes P-wave travel time from each grid node to each station:
   - `'h2p+ak135'`: table lookup via `data_util.travel_time_function()`, interpolated with `scipy.interp1d`
   - `'constant'`: `travel_time = distance / 6.0`

5. Estimates the average origin time at each grid node:
   ```
   OT[node] = mean( trigger_time[sta] - travel_time[node, sta] )   over all stations
   ```

6. Computes forward-predicted trigger times: `t_calc[node, sta] = OT[node] + travel_time[node, sta]`

7. Computes squared residuals: `tt_error[node, sta] = |t_calc - t_obs|²`

8. Sums residuals across stations → `misfit[node]` (p×p array)

9. Gaussian per-station likelihood: `rho[node, sta] = exp(-0.5 * tt_error / sigma²)`
   where `sigma` comes from the `.run` file's `sigma` column (initialized to 1.0 for all stations in `run_bEPIC`).

10. Combined likelihood: `like[node] = product over stations of rho[node, sta]` → p×p array

11. Finds the grid node with the maximum likelihood → `likelihood_lon`, `likelihood_lat`.

**Returns:** `(like, misfit, likelihood_lon, likelihood_lat)`

---

### prior.py

#### `generate_prior_seismicity_catalog()`

Downloads the ANSS historical seismicity catalog using `libcomcat.search` and saves it to `data/prior_seismicity_catalog.txt` (TSV).

**Region:** lon -135 to -112, lat 30 to 50 (western US)
**Filter:** M ≥ 3, from 2000-01-01 to present
**Columns:** `ANSS ID, date, timestamp, lon, lat, depth, mag`

Called automatically by `compute_prior()` if the catalog file doesn't exist.

---

#### `compute_prior(CenterPoint, GridSize, GridSpacing, ANSS_timestamp=None) → (prior_seis, prior_lon, prior_lat)`

Computes a 2D spatial prior representing historical seismicity density.

**Method:** Gaussian KDE over ANSS earthquake locations, evaluated on the search grid.

**Steps:**

1. Checks if `data/prior_seismicity_catalog.txt` exists; builds it if not.
2. Builds the search grid with `geospatial_util.make_grid()`.
3. If `ANSS_timestamp` is `None`, uses `datetime.now().timestamp()` (appropriate for real-time use; for replay events, pass the event's origin time to avoid including future seismicity in the prior).
4. Filters the catalog to the grid's bounding box, M ≥ 3, and timestamp < `ANSS_timestamp`.
5. Converts earthquake lon/lat to Cartesian km.
6. Fits `scipy.stats.kde.gaussian_kde` with Scott's bandwidth rule.
7. **Isotropizes the covariance:** replaces the KDE's covariance matrix with `I * max_eigenvalue`. This removes directional bias from anisotropic seismicity distributions.
8. Evaluates the KDE on the full grid → `prior_seis` (p×p array).
9. Finds the grid node with the highest prior density → `prior_lon`, `prior_lat`.

**Returns:** `(prior_seis, prior_lon, prior_lat)`

---

### posterior.py

#### `compute_posterior(CenterPoint, GridSize, GridSpacing, prior_function, likelihood_function) → (post, posterior_lon, posterior_lat)`

Combines prior and likelihood into a normalized posterior using Bayes' theorem.

**Method:**
```
k    = ∬ prior(x,y) × likelihood(x,y) dx dy     [normalization constant; computed via np.trapz]
post = (1/k) × likelihood × prior
```

The double integral is computed using `np.trapz` twice (once over x, once over y) on the grid coordinates.

Finds the MAP (maximum a posteriori) grid node → `posterior_lon`, `posterior_lat`.

**Returns:** `(post, posterior_lon, posterior_lat)`

---

### magnitude.py

#### `compute_magnitude(run_df, version, epicenter) → float`

Estimates earthquake magnitude from peak ground displacement (Pd) at each station.

**Scaling relation (Chung et al. 2020, BSSA):**
```
M = a·log10(Pd) + b·log10(R) + c
    where a=1.23, b=1.39, c=5.39
    Pd = peak displacement [m]
    R  = hypocentral distance [km] (epicentral, no depth term)
```

**Steps:**
1. Filters `run_df` to the given version.
2. Computes distance from the given epicenter to each station using `geospatial_util.get_dist_between_two_points_km()`.
3. Excludes stations at distance = 0 or > 200 km (set to NaN).
4. Converts `logPd` column: `Pd = 10^logPd`.
5. Handles any zero Pd values by setting them to NaN.
6. Computes per-station magnitude and returns the `nanmean`.

---

### geospatial_util.py

Coordinate transforms and grid construction. All functions assume a locally flat Earth approximation with an oblate spheroid correction.

---

#### `make_grid(CenterPoint, GridSize=200, GridSpacing=2) → (grid_lons, grid_lats, grid_x, grid_x, grid_x_ravel, grid_y_ravel)`

Builds a square 2D search grid centered on `CenterPoint`.

**Grid construction:**
```python
grid_inc = np.arange(GridSpacing, GridSize + GridSpacing, GridSpacing)
grid_x   = np.hstack((grid_inc[::-1]*-1, 0, grid_inc))  # e.g. [..., -4, -2, 0, 2, 4, ...] km
```
Number of nodes per axis: `2 * (GridSize / GridSpacing) + 1`  (e.g., 201 for GridSize=200, GridSpacing=2)

Converts the 1D km offsets to lon/lat via `ckm2LLd()`. Returns raveled meshgrid coordinates for vectorized operations.

**Returns:** `(grid_lons, grid_lats, grid_x, grid_x, grid_x_ravel, grid_y_ravel)`
> Note: the 3rd and 4th return values are the same array `grid_x` (both x and y are symmetric).

---

#### `LL2cartd(lon, lat, lon0, lat0, rot) → [x_rot, y_rot]`

Converts lon/lat arrays to local Cartesian coordinates (meters) centered on `(lon0, lat0)` with optional rotation by `rot` degrees.

- North = +y, East = +x before rotation.
- All inputs must be numpy arrays of equal length.
- Originally written by Amanda Thomas (2006), modified by Lujia Feng (2008), translated to Python by Amy Williamson (2018).

---

#### `ckm2LLd(xx, yy, lon0, lat0, rot) → [lon, lat]`

Inverse of `LL2cartd`. Converts Cartesian meters back to lon/lat. Used by `make_grid()`.

---

#### `get_dist_between_two_points_km(lon1, lat1, lon2, lat2) → float`

Returns the geodesic distance in km between two points using `obspy.geodetics.gps2dist_azimuth`. Note argument order: **lon first, then lat** (unlike obspy's native lat-first convention — the wrapping handles the swap internally).

---

### data_util.py

Handles all I/O: EPIC log parsing, `.run` file generation, USGS queries, and travel time interpolation.

---

#### `generate_run_file(project_parent_directory, postgres_id)`

Parses `EPIC/{postgres_id}_event_triggers_log.txt` and writes `{postgres_id}.run`.

**Source columns used from event_triggers_log.txt** (by column index, 0-based after skipping header):
- col 1: version
- col 3: order
- cols 4,5,6,7: sta, chan, net, loc
- cols 8,9: lat, lon
- col 10: trigger time (ISO string → converted to Unix timestamp)
- col 16: logPd
- col 34: tterr

---

#### `search_for_USGS_event(project_parent_directory, postgres_id)`

Searches USGS ComCat for a matching event by expanding a time+distance window around EPIC's reported location and time.

**Algorithm:**
- Starts with dt=0 seconds and dx=0 km.
- Iteratively expands: +5 seconds and +5 km each step if no event found.
- Shrinks by 1 if multiple events returned (to tighten the window).
- Stops when exactly 1 event is found or dx > 100 km (gives up, fills with NaN).

Saves result to `USGS/usgs_event_summary.txt` (TSV). Columns: `postgres id, USGS ID, USGS time, USGS lat, USGS lon, USGS depth, USGS mag`.

---

#### `search_comcat(starttime, endtime, maxradiuskm, latitude, longitude) → dict`

Low-level USGS ComCat API wrapper. Queries `earthquake.usgs.gov/fdsnws/event/1/query?format=geojson` with given spatiotemporal bounds. Returns the parsed GeoJSON dict. Timeout: 60 seconds.

---

#### `search_comcat_by_eventid(project_parent_directory, postgres_id, eventid) → DataFrame`

Fetches a specific USGS event by its ComCat event ID (e.g. `'nc73093981'`). Converts the millisecond timestamp in the GeoJSON to a Unix timestamp. Saves and returns a DataFrame with the same columns as `usgs_event_summary.txt`.

---

#### `travel_time_function(velocity_model) → interpolator`

Loads `$BEPIC/data/h2p+ak135.080` and returns a `scipy.interpolate.interp1d` function mapping distance (km) → travel time (seconds). The table has two columns: distance, time.

The `velocity_model` argument is accepted but currently unused — the h2p+ak135 table is always loaded regardless. The branching between table lookup and constant-velocity is done in `likelihood.py` itself.

---

#### `parse_log(project_parent_directory, log_file, event_id, epic_id)`

Parses a raw EPIC `.log` text file and extracts 7 DataFrames by scanning for tagged line prefixes:

| Tag | DataFrame | Description |
|-----|-----------|-------------|
| `E:I:` or `E:I:F:` | `event_summary_df` | Per-version event summary (44 fields) |
| `E:I:T:` | `event_triggers_df` | Per-version station trigger details (39 fields) |
| `L:T:` | `location_triggers_df` | Per-version location algorithm triggers (14 fields) |
| `E:S:` | `station_summary_df` | Per-version station count summary (13 fields) |
| `A:` | `a_df` (misc) | Station association details (32 fields) |
| `E:C:` | `station_counts_df` | Station count details (15 fields) |
| `L:E:` | `epic_location_df` | EPIC's grid search output (16 fields) |

Only lines where the event ID field matches `epic_id` are extracted. All 7 files are saved as TSV to `{event_id}/EPIC/`. Also creates the event and EPIC directories if they don't exist.

**Args:**
- `event_id`: the zero-padded postgres ID string (used for directory and filenames)
- `epic_id`: the EPIC internal event ID string (used to filter log lines)

---

### bEPIC_analysis.py

#### `compute_station_trigger_misfit(postgres_id, project_parent_directory)`

Post-processing: computes trigger time residuals relative to the USGS catalog location (not the EPIC/bEPIC estimate).

**Method:**
1. Loads `USGS/usgs_event_summary.txt` for the catalog origin.
2. Loads the `.run` file; uses only the **last version** (most complete station list).
3. Flags which stations have `tterr > -999` (those used in location).
4. Computes station distances from the USGS catalog epicenter (Cartesian, includes USGS depth).
5. Computes predicted trigger time: `USGS_time + distance/6.0` (constant 6 km/s — hardcoded).
6. Computes offset: `observed_trigger_time - predicted_trigger_time`.

Saves to `USGS/usgs_trigger_time_misfit.txt` (TSV). Columns: `station, network, channel, station lon, station lat, trigger offset, station distance, station in location`.

---

### pull_logs.py

CLI script to download and parse EPIC archive logs.

**Usage:**
```bash
python src/bEPIC/pull_logs.py <instance> <day> <postgres_id> <epic_id>
```

**Args:**
- `instance`: EPIC instance string, e.g. `epic@eew-nc-prod1` (the `epic@` prefix is stripped automatically)
- `day`: date string used to construct the log filename, e.g. `2018-10-05`
- `postgres_id`: integer event ID
- `epic_id`: EPIC internal event ID (used to filter log lines)

**What it does:**
1. Constructs the log URL: `http://131.215.66.120/{instance}/epic/epic_{day}.log.gz`
2. Downloads the `.log.gz` to a local cache directory (`/home/gcl/RA/williamson/EPIC_unprocessed_logs/`).
3. Decompresses with `gzip`.
4. Calls `data_util.parse_log()`.
5. Removes the `.gz` file (but keeps the decompressed `.log`).

> **Hardcoded paths:** `project_parent_directory` and `log_location` point to Caltech GCL servers. Must be updated for use on other systems.

---

### generate_synthetic_run_file.py

Generates a synthetic `.run` file for testing bEPIC without real data.

**Synthetic event:** M4.0 near Berkeley, CA
- EQ_LON = -122.27, EQ_LAT = 37.87
- EQ_DEPTH_KM = 8.0 (matches `calculate_likelihood` hardcoded depth)
- EQ_ORIGIN = 2024-01-15 10:30:00 UTC
- VELOCITY = 6.0 km/s

**6 synthetic stations:** BKS, OAKD, MCCM, SUTB, FARB, WENL (real NCEDC stations, BK and NC networks, Northern California)

**Output:** Trigger times computed as `OT + slant_distance/6.0 + N(0, 0.1)` noise. Each station appears in only its own version (version = station index), so version 0 has 1 station, version 5 has all 6.

There are two copies:
- `src/bEPIC/generate_synthetic_run_file.py` — writes to `/tmp/bEPIC_test/999999/`
- `generate_synthetic_run_file.py` (root) — writes to `./bEPIC_test/999999/`

---

## 8. Key Parameters and Constants

| Parameter | Location | Value | Description |
|-----------|----------|-------|-------------|
| `eq_depth` | `likelihood.py:37` | 8.0 km | Fixed assumed depth for all distance calculations |
| `velocity_model='constant'` | `likelihood.py:68` | 6.0 km/s | Constant P-wave velocity fallback |
| `mag_floor` | `prior.py:75` | M ≥ 3 | Minimum magnitude for ANSS prior catalog |
| `a, b, c` | `magnitude.py:27-29` | 1.23, 1.39, 5.39 | Pd magnitude scaling coefficients (Chung et al. 2020) |
| max station distance | `magnitude.py:42` | 200 km | Stations beyond this are excluded from magnitude |
| `GridSize` | `bEPIC_main` / caller | 200 km (typical) | Grid half-width |
| `GridSpacing` | `bEPIC_main` / caller | 2 km (typical) | Grid node spacing |
| KDE bandwidth | `prior.py:97` | Scott's rule | Used by `scipy.stats.kde.gaussian_kde` |
| Bad pick threshold | `bEPIC_main.py:120` | `tterr > -999` | Rows with tterr ≤ -999 are excluded |
| ComCat timeout | `data_util.py:285` | 60 seconds | HTTP timeout for USGS queries |
| ComCat limit | `data_util.py:284` | 20000 | Max events returned per ComCat query |

---

## 9. Output Files

### `{postgres_id}_bEPIC_log.txt` (main output)

Tab-separated, one row per version. Written by `run_bEPIC()`.

| Column | Description |
|--------|-------------|
| `version` | Version index (0 = first alert) |
| `num stations` | Number of stations used |
| `likelihood lon/lat` | MAP of likelihood surface alone |
| `likelihood mag` | Magnitude estimated at likelihood MAP |
| `prior lon/lat` | MAP of prior surface alone |
| `posterior lon/lat` | MAP of posterior (recommended best estimate) |
| `posterior mag` | Magnitude estimated at posterior MAP |

### `usgs_event_summary.txt`

Tab-separated, one row. Columns: `postgres id, USGS ID, USGS time, USGS lat, USGS lon, USGS depth, USGS mag`. Written by `search_for_USGS_event()` or `search_comcat_by_eventid()`.

### `usgs_trigger_time_misfit.txt`

Tab-separated, one row per station (last version only). Written by `bEPIC_analysis.compute_station_trigger_misfit()`. Columns: `station, network, channel, station lon, station lat, trigger offset, station distance, station in location`.

---

## 10. Testing

Tests use `pytest` and live in `tests/`.

### `tests/test_locate.py`

Tests `locate.get_two_station_location()`:

| Test | What it checks |
|------|----------------|
| `test_first_station_triggers_first_weighted_toward_it` | Output is 2/3 weighted toward the earlier-triggering station |
| `test_second_station_triggers_first_weighted_toward_it` | Same, with station 2 triggering first |
| `test_simultaneous_triggers_uses_station_1_weighting` | Ties go to the station-1 branch |
| `test_output_is_list_of_two_floats` | Return type and length |

### `tests/test_geospatial_util.py`

Tests coordinate transform functions:

| Test | What it checks |
|------|----------------|
| `test_LL2cartd_same_point_is_zero` | A point converted relative to itself is at origin |
| `test_LL2cartd_north_is_positive_y` | North = +y axis |
| `test_LL2cartd_east_is_positive_x` | East = +x axis |
| `test_round_trip_LL_to_cart_and_back` | `LL2cartd` → `ckm2LLd` recovers input coordinates to 1e-4° |
| `test_distance_same_point_is_zero` | Zero distance to self |
| `test_distance_LA_to_SF_approximate` | LA–SF ~560 km (530–600 km range check) |

### Running tests

```bash
cd /home/daned/2024_NEHRP/bEPIC
pytest tests/
```

---

## 11. Zextra — EPIC C Locator Tools

The `Zextra/` directory contains standalone scripts for (a) replicating the EPIC C grid search algorithm in Python and (b) downloading real USGS event data as inputs. These scripts are **not** part of the installed `bEPIC` package. They contain dependencies on other repositories for testing purposes only.

---

### EPIC_locate_prelim.py

A faithful Python transliteration of the EPIC C/C++ `searchGrid` function. Named "EPIC C" to indicate it mirrors the C source code structure rather than using idiomatic Python/NumPy vectorization.

#### Classes

**`Event`**
```python
Event(lat, lon, time, misfit_rms, misfit_ave, eventid, version)
```
Holds the current event estimate. Key fields: `.lat`, `.lon`, `.time`, `.depth` (fixed 8.0), `.trigs` (list of `TriggerManager`).

**`TriggerManager`**
```python
TriggerManager(lon, lat, trigger_time, sta, net, chan)
```
One per station trigger. Also stores `.dist`, `.tt`, `.tterror` (populated after location).

**`EPIC_PARAMS`**
```python
params = EPIC_PARAMS()
params.PriorGridFile = '/path/to/prior_seis_grid_US_Canada.tt3'
params.use_prior = True
params.GridSize  = 25      # nodes from center to edge
params.GridKm    = 50      # km from center to edge
params.method    = 'EPIC C'
```
Grid spacing = `GridKm / GridSize` km (e.g., 50/25 = 2 km).

**`PriorFile`**
Reads the `.tt3` binary prior grid format. Header (6 lines): `mx`, `my`, `xlower`, `ylower`, `dx`, `dy`. Data: `mx × my` floats, read with `np.loadtxt`, flipped vertically with `np.flipud`.

**`SearchOut`**
Output container. Key fields: `posterior_lon`, `posterior_lat`, `best_OT`, `best_location_post`, `best_like`, `best_prior`, `misfit_ave`, `best_misfit`.

#### Functions

**`latLonToXY(event) → event`**
Converts all station lon/lat to local Cartesian (km) relative to `event.lat/lon`. Populates `trig.stax` and `trig.stay`. Uses oblate spheroid radius correction.

**`E2Location_locate(params, event) → (t, output_df)`**
Main location function. Calls `latLonToXY()`, then `E2Location_searchGrid()`. After finding the best location, computes per-station distances and travel times from the MAP point. Prints formatted results to stdout.

**`E2Location_searchGrid(params, event) → (t, output_df)`**
The grid search. Iterates over every grid node (nested loop: y outer, x inner), then over every station, computing:
1. Distance from grid node to station (flat-earth, no depth term).
2. Travel time via the h2p+ak135 table (`ttf(dist)`).
3. Implied origin time per station: `trigger_time - travel_time`.
4. Average OT and RMS/mean misfit across stations.
5. Per-station Gaussian likelihood: `exp(-0.5 * tterror²)`.
6. Combined likelihood (product) × prior → posterior.

Tracks the MAP grid point. Returns `SearchOut` and a DataFrame with one row per grid node: `y, x, lat, lon, like, prior, post, misfitrms, misfitave`.

The `'python bypass'` method (vectorized NumPy) is present in the code but deprecated and removed from use.

**`get_dist_between_two_points_km(lon1, lat1, lon2, lat2)`**
Same as `geospatial_util.get_dist_between_two_points_km()` — a local copy using `obspy.geodetics.gps2dist_azimuth`.

#### Prior grid format (`.tt3`)

The `prior_seis_grid_US_Canada.tt3` file used by `PriorFile`:
- 6-line header: `mx` (1375), `my` (1613), `xlower` (-140.0), `ylower` (28.0), `dx` (~0.02037°), `dy` (~0.01798°)
- Coverage: lon -140 to ~-112, lat 28 to ~57 (US + Canada)
- Data: `mx × my` float values in row-major order, **y-flipped on load**

---

### example_call_EPIC_locate.py

Minimal working example for `EPIC_locate_prelim`. Reproduces event 126625 (M3.8, Gilroy CA, 2018-10-05) using 5 station triggers hand-copied from the `.run` file. Seeds the Event with the first station's location and trigger time, then calls `E2Location_locate()`.

---

### download_usgs_event.py

Downloads a real USGS event and its phase picks, then formats them as `Event` and `TriggerManager` objects for use with `EPIC_locate_prelim.E2Location_locate()`.

**Usage:**
```bash
python Zextra/download_usgs_event.py <usgs_eventid> [options]

# Examples:
python Zextra/download_usgs_event.py us7000n72h
python Zextra/download_usgs_event.py us7000n72h --max-dist 8.0 --phases P Pn Pg
python Zextra/download_usgs_event.py us7000n72h --run --prior-grid /path/to/prior_seis_grid_US_Canada.tt3
```

**Options:**
| Flag | Default | Description |
|------|---------|-------------|
| `--max-dist` | 5.0° | Max epicentral distance (degrees) for phase filtering |
| `--phases` | P Pn Pg Pb | Phase types to include |
| `--run` | False | Also call `E2Location_locate()` after building the event |
| `--prior-grid` | (hardcoded path) | Path to `.tt3` prior grid file for `--run` |

**API calls:**
1. `GET earthquake.usgs.gov/fdsnws/event/1/query?eventid={id}&format=geojson` — event origin
2. `GET {phases.csv url from GeoJSON product list}` — phase picks (same format as `Zextra/phases.csv`)
3. `GET service.iris.edu/fdsnws/station/1/query?net=...&sta=...` — station coordinates (cached per net+sta)

**Phase CSV format:** `Channel, Distance, Azimuth, Phase, Arrival Time, Status, Residual, Weight`
Channel string: `"NET STA CHA LOC"` (space-separated).

The USGS catalog origin is used to seed the Event's initial lat/lon/time, so the grid is pre-centered on the known answer. To simulate an EEW scenario, override with the first station's coordinates before calling locate.

---

### plot_prior.py

Visualizes the `prior_seis_grid_US_Canada.tt3` file. Reads the `.tt3` header and data using a local `PriorFile` class (same logic as in `EPIC_locate_prelim.py`), reconstructs the lon/lat grid with `np.linspace`, and plots as a scatter of squares colored by prior value on a log scale using `matplotlib`.

---

## 12. Catalog and Testing Data

### `Zextra/bEPIC_testing_catalog.txt`

Tab-separated. ~740 events from 2018–2022. Used for batch testing bEPIC performance.

Columns: `ANSS ID, ANSS date, ANSS mag, ANSS lat, ANSS lon, ANSS depth, postgres id, EPIC ID, EPIC instance, EPIC mag, EPIC lat, EPIC lon, EPIC OT, EPIC alert time, distance error`

- `postgres id`: links to `.run` files in `Zextra/run_files/`
- `distance error`: km between ANSS and EPIC locations
- Some rows have `nan` for EPIC OT, alert time, and distance error (events where EPIC data is incomplete)

### `Zextra/run_files/`

Contains ~600+ real `.run` files named by postgres_id (e.g. `126625.run`). Directly usable as input to `run_bEPIC()`. These cover the same events as `bEPIC_testing_catalog.txt`.

---

## 13. Known Issues and Notes

- **Hardcoded `bepic` path:** `prior.py` and `bEPIC_main.py` set `bepic = "/home/daned/2024_NEHRP/bEPIC"` directly instead of reading `$BEPIC`. If running on a different machine, update this or ensure the env var is set and change these lines to `bepic = os.environ['BEPIC']`.

- **Hardcoded archive paths in `pull_logs.py`:** `project_parent_directory` and `log_location` point to Caltech GCL servers (`/home/gcl/RA/williamson/...`). Must be updated for use elsewhere.

- **Velocity model argument to `travel_time_function()`:** The `velocity_model` parameter is accepted but unused — the function always loads `h2p+ak135.080`. The actual branching happens in `likelihood.py` where `velocity_model == 'constant'` uses 6 km/s inline.

- **`bEPIC_analysis.compute_station_trigger_misfit()` uses constant velocity (6 km/s)** regardless of the `velocity_model` used in `run_bEPIC()`.

- **`test_1.py` imports `torch`:** This is an ad-hoc development script and is not part of the production pipeline. The torch import is non-functional in this context.

- **`sigma` column:** Initialized to all-ones in `run_bEPIC()`. The `tterr` column from the `.run` file is available but not used to weight the likelihood. Future versions could use per-station uncertainties here.

- **2D grid, fixed depth:** All location is done in 2D (lon/lat only). Depth is fixed at 8.0 km everywhere in the likelihood computation and in `bEPIC_analysis`. This is a known simplification inherited from EPIC.
