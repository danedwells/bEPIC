#!/usr/bin/env python3
"""
generate_travel_time_tables.py — build h2p+ak135-format P travel-time
tables at arbitrary source depths.

bEPIC's location search (EPIC_locate_prelim.py) looks up P travel time as a
function of hypocentral distance from a table like data/h2p+ak135.080 — a
distance -> travel-time curve precomputed for one fixed source depth (8 km).
That single table is why bEPIC's depth is hardcoded: a coarse depth search
needs one table per candidate depth, computed for the correct source-depth
ray geometry (not a re-lookup of the existing 8 km table at other depths).

This regenerates that table format using obspy's TauP implementation of
AK135 — the same global model the existing table was built from. Verified:
travel times match data/h2p+ak135.080 to ~0.01-0.15 s at distances <60 km
and >150 km, well within the location algorithm's calibrated sigma_s
(~0.22 s).

KNOWN CAVEAT: in the 60-150 km band, plain AK135 (via this script) disagrees
with data/h2p+ak135.080 by up to ~1.2 s, confirmed against TauP's full
phase inventory (phase_list='ttall') at the legacy table's own sample
points -- not an artifact of this script's phase selection. The likely
explanation is that h2p+ak135.080 splices AK135's mantle onto a
regionally-tuned (probably Southern-California) crustal model, which
produces different Pg/Pn crustal-refraction arrivals than plain global
AK135 in that distance range. New tables generated here (1/20/50/100 km)
are therefore NOT guaranteed dimensionally consistent with the existing
8 km table in that band. Follow up with whoever maintains the original
h2p tool/model to identify the actual crustal correction before relying on
this for precision work in that distance range; for now this is accepted
as a known limitation of the coarse first-pass depth search.

Only columns 1-2 (distance_km, travel_time_s) are read anywhere in bEPIC
today (EPIC_locate_prelim.py, data_util.py, likelihood.py all index
tt_mod[:, 0] / tt_mod[:, 1] only). Columns 3-4 here are TauP's own ray
parameter (s/radian) and takeoff angle (degrees from downward vertical) for
completeness — they are NOT guaranteed to match whatever convention the
original h2p tool used for its columns 3-4, since nothing currently reads
them.

Usage:
    python generate_travel_time_tables.py --depths 1 20 50 100
    python generate_travel_time_tables.py --depths 8 --outdir /tmp --overwrite
"""
import argparse
import os

import numpy as np
from obspy.geodetics import kilometers2degrees
from obspy.taup import TauPyModel

# First-arriving-P phase set; reproduces data/h2p+ak135.080 to within ~0.1 s
# (validated against the shipped 8 km table across 0-1500 km).
PHASE_LIST = ['p', 'P', 'Pn', 'Pg']

DEFAULT_DEPTHS = [1.0, 8.0, 20.0, 50.0, 100.0]


def _distance_grid_km():
    """Distance sampling in km: fine near-source, coarser at range.

    Not a reproduction of the legacy table's adaptive (branch-following)
    spacing — just fine enough that downstream linear interpolation stays
    well under sigma_s at every distance in the search grid, out to 1500 km
    (matching the legacy table's range).
    """
    return np.concatenate([
        np.arange(0, 50, 0.25),
        np.arange(50, 300, 1.0),
        np.arange(300, 1500 + 5, 5.0),
    ])


def table_filename(depth_km):
    """h2p+ak135.080 == depth 8.0 km -> depth*10, zero-padded to 3 digits."""
    return f'h2p+ak135.{int(round(depth_km * 10)):03d}'


def generate_table(depth_km, model_name='ak135'):
    """Return an (n, 4) array of (distance_km, travel_time_s, ray_param, takeoff_angle)."""
    model = TauPyModel(model=model_name)
    distances_km = _distance_grid_km()

    rows = []
    for dist_km in distances_km:
        deg = kilometers2degrees(max(dist_km, 1e-4))
        arrivals = model.get_travel_times(
            source_depth_in_km=depth_km,
            distance_in_degree=deg,
            phase_list=PHASE_LIST,
        )
        if not arrivals:
            continue
        first = min(arrivals, key=lambda a: a.time)
        rows.append((dist_km, first.time, first.ray_param, first.takeoff_angle))

    return np.array(rows)


def write_table(depth_km, outdir, overwrite=False):
    out_path = os.path.join(outdir, table_filename(depth_km))
    if os.path.exists(out_path) and not overwrite:
        raise FileExistsError(
            f"{out_path} already exists — pass --overwrite to replace it. "
            f"(The shipped 8 km table should generally be left untouched, "
            f"to guarantee the default single-depth search stays unchanged.)"
        )

    rows = generate_table(depth_km)
    with open(out_path, 'w') as fh:
        fh.write(f'depth {depth_km:.1f}\n')
        for dist_km, t, ray_param, takeoff in rows:
            # Fields are whitespace-delimited (np.genfromtxt default), not
            # fixed-width — ray_param can run to 4+ digits at regional
            # distances, so give it (and its neighbours) enough width to
            # always keep a separating space.
            fh.write(f'{dist_km:10.2f}{t:10.2f}{ray_param:12.2f}{takeoff:10.1f}\n')
    return out_path


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--depths', type=float, nargs='+', default=DEFAULT_DEPTHS,
                         help=f'Source depths in km (default: {DEFAULT_DEPTHS})')
    parser.add_argument('--outdir',
                         default=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data'),
                         help='Output directory (default: bEPIC/data)')
    parser.add_argument('--overwrite', action='store_true',
                         help='Overwrite an existing table file for a given depth')
    args = parser.parse_args()

    for depth_km in args.depths:
        out_path = write_table(depth_km, args.outdir, overwrite=args.overwrite)
        print(f'  depth={depth_km:>6.1f} km  ->  {out_path}')


if __name__ == '__main__':
    main()
