#!/usr/bin/env python3
"""
catalog_fetch.py — ANSS catalog download for bEPIC's seismicity prior.

Replaces the libcomcat dependency in prior.py by querying the USGS FDSNWS
event API directly via requests.  Results are identical: a tab-separated
catalog file at $BEPIC/data/prior_seismicity_catalog.txt.

The FDSNWS API caps responses at 20,000 events per request, so the download
is chunked into one request per calendar year to stay well under that limit.
"""

import io
import os
from datetime import datetime

import pandas as pd
import requests

FDSNWS_URL = "https://earthquake.usgs.gov/fdsnws/event/1/query"
HEADERS    = {"User-Agent": "bEPIC catalog_fetch"}
TIMEOUT    = 120  # seconds per request


def generate_prior_seismicity_catalog():
    """
    Download the ANSS catalog and save it to $BEPIC/data/prior_seismicity_catalog.txt.

    Queries M≥3 earthquakes in the western US (lon -135 to -112, lat 30 to 50)
    from 2000 to the present, one calendar year per request to avoid the
    FDSNWS 20,000-event response limit.

    Output columns match the format expected by prior.compute_prior():
        ANSS ID | date | timestamp | lon | lat | depth | mag
    """
    bepic = os.path.dirname(os.path.abspath(__file__))
    region  = [-135, -112, 30, 50]   # lon_min, lon_max, lat_min, lat_max
    now     = datetime.now()
    chunks  = []

    for year in range(2000, now.year + 1):
        start = f"{year}-01-01"
        end   = f"{year + 1}-01-01" if year < now.year else now.strftime("%Y-%m-%d")

        params = {
            "format":       "csv",
            "starttime":    start,
            "endtime":      end,
            "minlongitude": region[0],
            "maxlongitude": region[1],
            "minlatitude":  region[2],
            "maxlatitude":  region[3],
            "minmagnitude": 3,
            "limit":        20000,
            "orderby":      "time-asc",
        }

        response = requests.get(FDSNWS_URL, params=params,
                                headers=HEADERS, timeout=TIMEOUT)
        response.raise_for_status()

        chunk = pd.read_csv(io.StringIO(response.text))
        if not chunk.empty:
            chunks.append(chunk)

    raw = pd.concat(chunks, ignore_index=True)

    df = pd.DataFrame({
        "ANSS ID":   raw["id"],
        "date":      raw["time"],
        "timestamp": pd.to_datetime(raw["time"], utc=True).apply(lambda t: t.timestamp()),
        "lon":       raw["longitude"],
        "lat":       raw["latitude"],
        "depth":     raw["depth"],
        "mag":       raw["mag"],
    })

    df.to_csv(bepic + "/data/prior_seismicity_catalog.txt", sep="\t", index=False)
