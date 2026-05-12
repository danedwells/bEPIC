
import numpy as np
import pytest
from scipy import interpolate
import os
import pandas as pd

from bEPIC.EPIC_locate_prelim import (
    Event, TriggerManager, EPIC_PARAMS, SearchOut, latLonToXY
)
import bEPIC.EPIC_locate_prelim as _ep

# Locate the bEPIC package root (needed for the travel-time table)
bepic = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(_ep.__file__))))


# ---------------------------------------------------------------------------
# Minimal mock prior
# ---------------------------------------------------------------------------

class _MockPrior:
    """
    Uniform prior over a coarse grid covering the test region.
    Gives use_prior=False tests a valid prior object without needing a .tt3 file.
    For use_prior=True tests, a non-uniform grid is used to exercise the lookup.
    """
    def __init__(self, uniform=True):
        self.lons = np.linspace(-125, -114, 60)
        self.lats = np.linspace(32,   43,   60)
        if uniform:
            self.grid = np.full((60, 60), 0.5)
        else:
            # Simple gradient so different grid points get different prior weights
            lon_idx, lat_idx = np.meshgrid(np.arange(60), np.arange(60))
            self.grid = (lat_idx + lon_idx + 1).astype(float)
            self.grid /= self.grid.max()


# ---------------------------------------------------------------------------
# Reference implementation: the original triple-nested loop
# (copied verbatim from the commented-out block in E2Location_searchGrid)
# ---------------------------------------------------------------------------

def _searchgrid_loop(event, trigs, params):
    """
    Reference: the original triple-nested loop from the commented-out section
    of E2Location_searchGrid (method='EPIC C').

    Copied verbatim for verification against the vectorized replacement.
    Row order matches the vectorized output: y-outer, x-inner.
    """
    evlat = event.lat
    evlon = event.lon

    num_trigs       = min(len(trigs), params.MAX_EVENT_TRIGS)
    trigs           = trigs[:num_trigs]
    trig_ot         = np.zeros(num_trigs)

    grid_size       = params.GridSize
    grid_km         = params.GridKm
    grid_spacing_km = grid_km / grid_size

    prior_info    = params.prior
    _prior_xlower = prior_info.lons[0]
    _prior_ylower = prior_info.lats[0]
    _prior_dx     = float(np.diff(prior_info.lons).mean())
    _prior_dy     = float(np.diff(prior_info.lats).mean())
    _prior_mx     = len(prior_info.lons)
    _prior_my     = len(prior_info.lats)

    vel_mod_filename = os.path.join(bepic, 'data', 'h2p+ak135.080')
    tt_mod = np.genfromtxt(vel_mod_filename, skip_header=1)
    ttf    = interpolate.interp1d(tt_mod[:, 0], tt_mod[:, 1])

    lat0  = evlat
    lon0  = evlon
    R     = 6378.137
    ff    = 1. / 298.257
    lat0r = lat0 * np.pi / 180.
    r     = R * (1 - ff * np.power(np.sin(lat0r), 2))
    mpd   = r * np.pi / 180.
    f     = mpd * np.cos(lat0r)

    ybeg   = -1 * grid_size
    yend   =  grid_size
    grid_y = grid_x = np.linspace(ybeg, yend, 2 * grid_size + 1)

    t    = SearchOut()
    rows = []
    j = i = 0

    for y in grid_y:
        ykm  = y * grid_spacing_km
        ylat = lat0 + ykm / mpd
        if params.use_prior == True:
            a = (ylat - _prior_ylower) / _prior_dy
            j = int(a + 0.5)
            if j < 0:        j = 0
            if j >= _prior_my: j = _prior_my - 1

        for x in grid_x:
            xkm   = x * grid_spacing_km
            sumOT = 0.0
            for it in range(num_trigs):
                tx  = trigs[it].stax - xkm
                ty  = trigs[it].stay - ykm
                dist = np.sqrt(tx * tx + ty * ty)
                stt  = float(ttf(dist))
                trig_ot[it] = trigs[it].time - stt
                sumOT += trig_ot[it]
            aveOT = sumOT / num_trigs

            rms = ttsum = 0.0
            for it in range(num_trigs):
                rms   += np.power(trig_ot[it] - aveOT, 2)
                ttsum += np.fabs(trig_ot[it] - aveOT)
            misfitsq   = rms   / num_trigs
            misfit_ave = ttsum / num_trigs

            like = 1.0
            for it in range(num_trigs):
                tterror = trig_ot[it] - aveOT
                like   *= np.exp(-0.5 * tterror * tterror)
            like = np.sqrt(like / num_trigs)

            Prior = 1.0
            xlon  = lon0 + xkm / f
            if params.use_prior:
                a = (xlon - _prior_xlower) / _prior_dx
                i = int(a + 0.5)
                if i < 0:        i = 0
                if i >= _prior_mx: i = _prior_mx - 1
                Prior = prior_info.grid[j, i]

            post = like * Prior

            rows.append({'y': y, 'x': x, 'lat': ylat, 'lon': xlon,
                         'like': like, 'prior': Prior, 'post': post,
                         'misfitrms': misfitsq, 'misfitave': misfit_ave})

            if post > t.best_location_post:
                t.best_location_post = post
                t.posterior_lon = xlon
                t.posterior_lat = ylat
                t.best_misfit   = misfitsq
                t.misfit_ave    = misfit_ave
                t.best_OT       = aveOT
                t.best_grid_x   = x
                t.best_grid_y   = y
                t.best_value    = post
                t.best_like     = like
                t.best_prior    = Prior

    return t, pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Shared fixture: event + triggers from example event 126625 (Gilroy M3.8)
# ---------------------------------------------------------------------------

def _make_event_and_trigs(use_prior):
    event = Event(
        lat=36.764, lon=-121.4472, time=1538771380.09,
        misfit_rms=0, misfit_ave=0, eventid=126625, version=0,
    )
    # Triggers copied directly from Amy's script (example_call_EPIC_locate.py)
    for lon, lat, tt, sta, net, chan in [
        (-121.4472, 36.764,  1538771380.09, 'SAO',  'BK', 'HNZ'),
        (-121.4472, 36.764,  1538771380.1,  'SAO', 'BK', 'HHZ'),
        (-121.5203, 36.6674, 1538771381.18, 'BSR',  'NC', 'HNZ'),
        (-121.287,  37.008,  1538771383.34, 'PACP', 'BK', 'HHZ'),
        (-121.287,  37.008,  1538771383.34, 'PACP', 'BK', 'HHZ')
    ]:
        event.trigs.append(
            TriggerManager(lon=lon, lat=lat, trigger_time=tt,
                           sta=sta, net=net, chan=chan)
        )
    event = latLonToXY(event)   # sets stax / stay on each trigger

    params = EPIC_PARAMS()
    params.prior           = _MockPrior(uniform=not use_prior)
    params.use_prior       = use_prior
    params.GridSize        = 5    # 11x11 grid — fast but non-trivial
    params.GridKm          = 10
    params.method          = 'EPIC C'
    params.MAX_EVENT_TRIGS = 10

    return event, event.trigs, params


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_grid_arrays_match_no_prior():
    """
    Per-grid-point values (like, prior, post, misfitrms, misfitave) from the
    loop and vectorized implementations must agree to floating-point precision
    when use_prior=False.
    """
    event, trigs, params = _make_event_and_trigs(use_prior=False)

    t_loop, df_loop = _searchgrid_loop(event, trigs, params)
    t_vec,  df_vec  = _ep.E2Location_searchGrid(event, trigs, params)

    # Sort both by (y, x) to guarantee row alignment
    df_loop = df_loop.sort_values(['y', 'x']).reset_index(drop=True)
    df_vec  = df_vec.sort_values(['y', 'x']).reset_index(drop=True)

    for col in ('like', 'prior', 'post', 'misfitrms', 'misfitave'):
        np.testing.assert_allclose(
            df_loop[col].values, df_vec[col].values,
            rtol=1e-10,
            err_msg=f"Column '{col}' differs between loop and vectorized (no prior)",
        )


def test_grid_arrays_match_with_prior():
    """
    Per-grid-point values must agree when use_prior=True with a non-uniform
    prior, exercising the prior lookup in both implementations.
    """
    event, trigs, params = _make_event_and_trigs(use_prior=True)

    t_loop, df_loop = _searchgrid_loop(event, trigs, params)
    t_vec,  df_vec  = _ep.E2Location_searchGrid(event, trigs, params)

    df_loop = df_loop.sort_values(['y', 'x']).reset_index(drop=True)
    df_vec  = df_vec.sort_values(['y', 'x']).reset_index(drop=True)

    for col in ('like', 'prior', 'post', 'misfitrms', 'misfitave'):
        np.testing.assert_allclose(
            df_loop[col].values, df_vec[col].values,
            rtol=1e-10,
            err_msg=f"Column '{col}' differs between loop and vectorized (with prior)",
        )


def test_best_location_scalars_match_no_prior():
    """
    The SearchOut best-location scalars (posterior_lon/lat, best_misfit,
    best_OT, best_like, best_prior, best_value) must match between the loop
    and vectorized implementations when use_prior=False.
    """
    event, trigs, params = _make_event_and_trigs(use_prior=False)

    t_loop, _ = _searchgrid_loop(event, trigs, params)
    t_vec,  _ = _ep.E2Location_searchGrid(event, trigs, params)

    for attr in ('posterior_lon', 'posterior_lat', 'best_misfit',
                 'misfit_ave', 'best_OT', 'best_like', 'best_prior', 'best_value'):
        assert abs(getattr(t_loop, attr) - getattr(t_vec, attr)) < 1e-10, (
            f"SearchOut.{attr} differs: loop={getattr(t_loop, attr)}, "
            f"vec={getattr(t_vec, attr)}"
        )


def test_best_location_scalars_match_with_prior():
    """
    The SearchOut best-location scalars must match when use_prior=True.
    """
    event, trigs, params = _make_event_and_trigs(use_prior=True)

    t_loop, _ = _searchgrid_loop(event, trigs, params)
    t_vec,  _ = _ep.E2Location_searchGrid(event, trigs, params)

    for attr in ('posterior_lon', 'posterior_lat', 'best_misfit',
                 'misfit_ave', 'best_OT', 'best_like', 'best_prior', 'best_value'):
        assert abs(getattr(t_loop, attr) - getattr(t_vec, attr)) < 1e-10, (
            f"SearchOut.{attr} differs: loop={getattr(t_loop, attr)}, "
            f"vec={getattr(t_vec, attr)}"
        )


# ---------------------------------------------------------------------------
# Activity mask tests
# ---------------------------------------------------------------------------

def _inventory_df(rows):
    """Build a station inventory DataFrame from (lon, lat, sta, net) tuples."""
    return pd.DataFrame(rows, columns=['longitude', 'latitude', 'station', 'network'])


def _params_with_inventory(inv_df, threshold=0.30):
    """Return a minimal EPIC_PARAMS with an activity mask inventory attached."""
    params = EPIC_PARAMS()
    params.prior                   = _MockPrior(uniform=True)
    params.use_prior               = False
    params.GridSize                = 5     # 11×11 grid, ±10 km
    params.GridKm                  = 10
    params.method                  = 'EPIC C'
    params.MAX_EVENT_TRIGS         = 10
    params.station_inventory       = inv_df
    params.activity_mask_threshold = threshold
    return params


def _single_trig_event():
    """One triggered station (SAO/BK) at the event epicentre."""
    event = Event(
        lat=36.764, lon=-121.4472, time=1538771380.09,
        misfit_rms=0, misfit_ave=0, eventid=999, version=0,
    )
    trig = TriggerManager(
        lon=-121.4472, lat=36.764, trigger_time=1538771380.09,
        sta='SAO', net='BK', chan='HNZ',
    )
    event.trigs = [trig]
    event = latLonToXY(event)
    return event, [trig]


def test_no_inventory_mask_column_is_ones():
    """With station_inventory=None the output_df must contain an all-ones activity_mask."""
    event, trigs, params = _make_event_and_trigs(use_prior=False)
    # station_inventory is None by default in EPIC_PARAMS
    _, df = _ep.E2Location_searchGrid(event, trigs, params)
    assert 'activity_mask' in df.columns
    assert (df['activity_mask'] == 1.0).all()


def test_inventory_triggered_only_mask_is_ones():
    """Inventory = triggered stations only → no untriggered → all-ones mask."""
    event, trigs, _ = _make_event_and_trigs(use_prior=False)
    rows = [(t.lon, t.lat, t.sta, t.net) for t in trigs]
    inv = _inventory_df(rows)
    params = _params_with_inventory(inv)
    _, df = _ep.E2Location_searchGrid(event, trigs, params)
    assert 'activity_mask' in df.columns
    assert (df['activity_mask'] == 1.0).all()


def test_dense_untriggered_masks_far_nodes():
    """
    Ten untriggered stations strung eastward from the epicentre should mask
    the eastern grid nodes but leave the centre and western nodes active.

    Setup: 1 triggered station at the event centre; 10 untriggered stations at
    Δlon = +0.01°, +0.02°, …, +0.10° (≈ 0.9, 1.8, …, 9.0 km to the east).

    Eastern edge node (XKM = +10 km):
      R_MAX = 10 km (distance to the single triggered station).
      All 10 untriggered stations lie within 10 km of this node.
      FRAC = 1 / 11 ≈ 0.09 < threshold → masked.

    Centre node (XKM = 0):
      R_MAX = 0 (triggered station is co-located with this node).
      No untriggered station lies within 0 km → FRAC = 1.0 → active.

    Western edge node (XKM = −10 km):
      R_MAX = 10 km, but all untriggered stations are ≥ 10.9 km away
      (east of centre) → none inside → FRAC = 1.0 → active.
    """
    event, trigs = _single_trig_event()

    utrig_rows = [(-121.4472 + 0.01 * k, 36.764, f'U{k:02d}', 'XX')
                  for k in range(1, 11)]
    inv = _inventory_df([(-121.4472, 36.764, 'SAO', 'BK')] + utrig_rows)
    params = _params_with_inventory(inv)

    _, df = _ep.E2Location_searchGrid(event, trigs, params)
    assert 'activity_mask' in df.columns

    centre = df[(df['x'] == 0.0) & (df['y'] == 0.0)]
    assert centre['activity_mask'].values[0] == 1.0, "Centre node should be active"

    east_edge = df[(df['x'] == 5.0) & (df['y'] == 0.0)]
    assert east_edge['activity_mask'].values[0] == 0.0, "Eastern edge node should be masked"

    west_edge = df[(df['x'] == -5.0) & (df['y'] == 0.0)]
    assert west_edge['activity_mask'].values[0] == 1.0, "Western edge node should be active"


def test_station_network_matching():
    """
    A station at the same coordinates as a triggered station but with a
    different (station, network) code must be treated as untriggered.
    Coordinate-based matching would incorrectly exclude it.
    """
    event, trigs = _single_trig_event()

    # SAO/BK is triggered. DXX/XX is at the exact same location but a
    # different code — it must count as untriggered.
    # Add 10 such duplicates so fraction drops below 0.30 at edge nodes.
    rows = [(-121.4472, 36.764, 'SAO', 'BK')]
    rows += [(-121.4472, 36.764, f'D{i:02d}', 'XX') for i in range(10)]
    inv = _inventory_df(rows)
    params = _params_with_inventory(inv)

    _, df = _ep.E2Location_searchGrid(event, trigs, params)
    assert 'activity_mask' in df.columns
    # Edge nodes (R_MAX > 0) should be masked because 10 untriggered stations
    # co-located with the triggered station lie within every non-zero radius.
    edge_nodes = df[df['x'].abs() == 5.0]
    assert (edge_nodes['activity_mask'] == 0.0).all(), (
        "Edge nodes should be masked when untriggered stations fill the circle"
    )


def test_mask_zeroes_posterior_and_preserves_unmasked():
    """
    Masked nodes (αm=0) must have post=0; active nodes must have
    post = like (since use_prior=False → prior=1 everywhere).
    """
    event, trigs = _single_trig_event()

    utrig_rows = [(-121.4472 + 0.01 * k, 36.764, f'U{k:02d}', 'XX')
                  for k in range(1, 11)]
    inv = _inventory_df([(-121.4472, 36.764, 'SAO', 'BK')] + utrig_rows)
    params = _params_with_inventory(inv)

    _, df = _ep.E2Location_searchGrid(event, trigs, params)

    masked = df[df['activity_mask'] == 0.0]
    assert len(masked) > 0, "Expected at least one masked node"
    assert (masked['post'] == 0.0).all(), "Masked nodes must have post = 0"

    active = df[df['activity_mask'] == 1.0]
    assert len(active) > 0, "Expected at least one active node"
    # use_prior=False → prior=1 → post = like * 1 * 1
    np.testing.assert_allclose(
        active['post'].values, active['like'].values,
        rtol=1e-10,
        err_msg="Active node post should equal like when use_prior=False",
    )
