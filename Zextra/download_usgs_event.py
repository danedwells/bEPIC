#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Download USGS event data and phase picks, then format them for use with
EPIC_locate_prelim.py (populating Event and TriggerManager objects).

Usage:
    python download_usgs_event.py <usgs_eventid>
    python download_usgs_event.py <usgs_eventid> --max-dist 5.0 --phases P Pn Pg

Arguments:
    eventid     USGS event ID (e.g., us7000n72h)
    --max-dist  Maximum epicentral distance in degrees to include (default: 5.0)
    --phases    Phase types to include (default: P Pn Pg Pb)
    --run       Also run EPIC locate after building event/triggers

The script:
    1. Queries USGS ComCat for event origin (lat, lon, depth, time, magnitude)
    2. Downloads the phases.csv product for that event
    3. Filters phases by type and distance
    4. Looks up station coordinates from IRIS FDSNWS
    5. Builds Event and TriggerManager objects ready for E2Location_locate()

Created: 2026-03-20
"""

import sys
import argparse
import requests
import pandas as pd
from io import StringIO
from datetime import datetime, timezone

import EPIC_locate_prelim


# ---------------------------------------------------------------------------
# USGS / IRIS data retrieval
# ---------------------------------------------------------------------------

def get_usgs_event(eventid):
    """Return the USGS ComCat GeoJSON for a given event ID."""
    url = "https://earthquake.usgs.gov/fdsnws/event/1/query"
    r = requests.get(url, params={"eventid": eventid, "format": "geojson"}, timeout=30)
    r.raise_for_status()
    return r.json()


def get_phases_df(geojson):
    """
    Download the phases.csv product from the event's USGS product list.

    Returns a DataFrame with columns:
        Channel, Distance, Azimuth, Phase, Arrival Time, Status, Residual, Weight
    or None if the product is unavailable.
    """
    products = geojson.get("properties", {}).get("products", {})
    phase_products = products.get("phase-data", [])
    if not phase_products:
        return None

    contents = phase_products[0].get("contents", {})
    if "phases.csv" not in contents:
        return None

    url = contents["phases.csv"]["url"]
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    return pd.read_csv(StringIO(r.text))


def get_station_coords(network, station, channel, origin_time_iso):
    """
    Query IRIS FDSNWS for the lat/lon of a station at the time of the event.

    Returns (lat, lon) or (None, None) if not found.
    """
    url = "https://service.iris.edu/fdsnws/station/1/query"
    params = {
        "net": network,
        "sta": station,
        "cha": channel,
        "level": "station",
        "format": "text",
        "endafter": origin_time_iso[:10],   # date only is enough
    }
    try:
        r = requests.get(url, params=params, timeout=15)
        r.raise_for_status()
        lines = [ln for ln in r.text.strip().split("\n") if not ln.startswith("#")]
        if lines:
            parts = lines[0].split("|")
            return float(parts[2]), float(parts[3])
    except Exception as e:
        pass
    return None, None


# ---------------------------------------------------------------------------
# Main builder
# ---------------------------------------------------------------------------

def build_event_and_triggers(eventid, max_dist_deg=5.0, phases_filter=None):
    """
    Download USGS data and return a populated EPIC_locate_prelim.Event object.

    Parameters
    ----------
    eventid : str
        USGS event ID, e.g. 'us7000n72h'
    max_dist_deg : float
        Keep only phases with epicentral distance <= this many degrees.
    phases_filter : list of str
        Phase types to include, e.g. ['P', 'Pn', 'Pg', 'Pb'].

    Returns
    -------
    event : EPIC_locate_prelim.Event  (or None on failure)
    """
    if phases_filter is None:
        phases_filter = ["P", "Pn", "Pg", "Pb"]

    # ------------------------------------------------------------------
    # 1. Event origin
    # ------------------------------------------------------------------
    print(f"\nFetching event {eventid} from USGS ComCat...")
    geojson = get_usgs_event(eventid)

    props  = geojson["properties"]
    coords = geojson["geometry"]["coordinates"]  # [lon, lat, depth_km]

    evlon   = coords[0]
    evlat   = coords[1]
    evdepth = coords[2]
    evmag   = props["mag"]
    evtime  = props["time"] / 1000.0   # ms → seconds since epoch

    origin_iso = datetime.fromtimestamp(evtime, tz=timezone.utc).isoformat()

    print(f"  Title : {props.get('title', eventid)}")
    print(f"  Origin: lat={evlat:.4f}  lon={evlon:.4f}  depth={evdepth} km  M={evmag}")
    print(f"  Time  : {origin_iso}")

    # ------------------------------------------------------------------
    # 2. Phase data
    # ------------------------------------------------------------------
    print("\nDownloading phases.csv from USGS...")
    phases_df = get_phases_df(geojson)
    if phases_df is None:
        print("  ERROR: No phases.csv product found for this event.")
        return None

    print(f"  Total phase picks : {len(phases_df)}")

    # Filter by phase type and distance
    phases_df = phases_df[phases_df["Phase"].isin(phases_filter)].copy()
    phases_df = phases_df[phases_df["Distance"] <= max_dist_deg].copy()
    print(f"  After filter (dist<={max_dist_deg}°, phases={phases_filter}): {len(phases_df)}")

    if phases_df.empty:
        print("  No phases remain after filtering. Try increasing --max-dist.")
        return None

    # Parse "NET STA CHA LOC" channel string
    def parse_channel(ch):
        parts = ch.strip().split()
        net = parts[0]
        sta = parts[1]
        cha = parts[2]
        loc = parts[3] if len(parts) > 3 else "--"
        return net, sta, cha, loc

    parsed = phases_df["Channel"].apply(parse_channel)
    phases_df["net"] = [p[0] for p in parsed]
    phases_df["sta"] = [p[1] for p in parsed]
    phases_df["cha"] = [p[2] for p in parsed]
    phases_df["loc"] = [p[3] for p in parsed]

    # ------------------------------------------------------------------
    # 3. Build Event (seeded with the USGS catalog origin)
    # ------------------------------------------------------------------
    event = EPIC_locate_prelim.Event(
        lat       = evlat,
        lon       = evlon,
        time      = evtime,
        misfit_rms= 0,
        misfit_ave= 0,
        eventid   = eventid,
        version   = 0,
    )

    # ------------------------------------------------------------------
    # 4. Fetch station coordinates and build TriggerManagers
    # ------------------------------------------------------------------
    print("\nFetching station coordinates from IRIS FDSNWS...")
    coord_cache = {}   # (net, sta) → (lat, lon)
    seen        = set()

    for _, row in phases_df.iterrows():
        key = (row["net"], row["sta"], row["cha"])
        if key in seen:
            continue
        seen.add(key)

        cache_key = (row["net"], row["sta"])
        if cache_key in coord_cache:
            sta_lat, sta_lon = coord_cache[cache_key]
        else:
            sta_lat, sta_lon = get_station_coords(row["net"], row["sta"], row["cha"], origin_iso)
            coord_cache[cache_key] = (sta_lat, sta_lon)

        if sta_lat is None:
            print(f"  SKIP  {row['net']}.{row['sta']}.{row['cha']} — coordinates not found")
            continue

        # Arrival time → Unix timestamp
        arrival_dt = datetime.fromisoformat(row["Arrival Time"].replace("Z", "+00:00"))
        arrival_ts = arrival_dt.timestamp()

        t = EPIC_locate_prelim.TriggerManager(
            lon          = sta_lon,
            lat          = sta_lat,
            sta          = row["sta"],
            net          = row["net"],
            chan         = row["cha"],
            trigger_time = arrival_ts,
        )
        event.trigs.append(t)
        print(f"  ADD   {row['net']}.{row['sta']}.{row['cha']:5s}  "
              f"lat={sta_lat:.4f}  lon={sta_lon:.4f}  "
              f"phase={row['Phase']}  dist={row['Distance']:.2f}°  "
              f"t={arrival_ts:.3f}")

    print(f"\nTotal triggers added: {len(event.trigs)}")
    return event


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Download a USGS event and build EPIC_locate_prelim inputs."
    )
    parser.add_argument("eventid",
                        help="USGS event ID (e.g. us7000n72h)")
    parser.add_argument("--max-dist", type=float, default=5.0,
                        help="Max epicentral distance in degrees (default: 5.0)")
    parser.add_argument("--phases", nargs="+", default=["P", "Pn", "Pg", "Pb"],
                        help="Phase types to include (default: P Pn Pg Pb)")
    parser.add_argument("--run", action="store_true",
                        help="Run EPIC locate after downloading (requires prior grid file)")
    parser.add_argument("--prior-grid",
                        default="/Users/amy/projects/container_bEPIC/data/prior_seis_grid_US_Canada.tt3",
                        help="Path to prior grid file (needed with --run)")
    args = parser.parse_args()

    event = build_event_and_triggers(
        args.eventid,
        max_dist_deg   = args.max_dist,
        phases_filter  = args.phases,
    )

    if event is None or len(event.trigs) < 2:
        print("\nNeed at least 2 triggers to run EPIC locate. Exiting.")
        sys.exit(1)

    if args.run:
        print("\n--- Running EPIC locate ---")
        params = EPIC_locate_prelim.EPIC_PARAMS()
        params.PriorGridFile = args.prior_grid
        params.use_prior     = True
        params.GridSize      = 25
        params.GridKm        = 50
        params.method        = "EPIC C"

        t, output_df = EPIC_locate_prelim.E2Location_locate(params, event)
    else:
        print("\nEvent and triggers ready. Pass --run to also execute EPIC locate.")
        print("Example code to run manually:\n")
        print("    params = EPIC_locate_prelim.EPIC_PARAMS()")
        print("    params.PriorGridFile = '/path/to/prior_seis_grid_US_Canada.tt3'")
        print("    params.use_prior = True")
        print("    params.GridSize  = 25")
        print("    params.GridKm    = 50")
        print("    params.method    = 'EPIC C'")
        print("    t, output_df = EPIC_locate_prelim.E2Location_locate(params, event)")
