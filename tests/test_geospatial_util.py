import numpy as np
import pytest
from bEPIC import geospatial_util


# ---------------------------------------------------------------------------
# LL2cartd
# ---------------------------------------------------------------------------

def test_LL2cartd_same_point_is_zero():
    """A point converted relative to itself should be at the origin."""
    lon, lat = np.array([-118.0]), np.array([34.0])
    x, y = geospatial_util.LL2cartd(lon, lat, lon0=-118.0, lat0=34.0, rot=0)
    assert abs(x[0]) < 1e-6
    assert abs(y[0]) < 1e-6


def test_LL2cartd_north_is_positive_y():
    """A point directly north of the origin should have positive y, zero x."""
    lon = np.array([-118.0])
    lat = np.array([35.0])   # 1 degree north of origin
    x, y = geospatial_util.LL2cartd(lon, lat, lon0=-118.0, lat0=34.0, rot=0)
    assert y[0] > 0
    assert abs(x[0]) < 1.0   # essentially zero (small floating point residual ok)


def test_LL2cartd_east_is_positive_x():
    """A point directly east of the origin should have positive x, zero y."""
    lon = np.array([-117.0])  # 1 degree east
    lat = np.array([34.0])
    x, y = geospatial_util.LL2cartd(lon, lat, lon0=-118.0, lat0=34.0, rot=0)
    assert x[0] > 0
    assert abs(y[0]) < 1.0


# ---------------------------------------------------------------------------
# ckm2LLd  (inverse of LL2cartd)
# ---------------------------------------------------------------------------

def test_round_trip_LL_to_cart_and_back():
    """Converting to Cartesian and back should recover the original coordinates."""
    lon0, lat0 = -118.0, 34.0
    lon_in = np.array([-119.5, -117.0, -118.0])
    lat_in = np.array([33.0,    35.0,   34.5])

    x, y = geospatial_util.LL2cartd(lon_in, lat_in, lon0=lon0, lat0=lat0, rot=0)
    lon_out, lat_out = geospatial_util.ckm2LLd(x, y, lon0=lon0, lat0=lat0, rot=0)

    np.testing.assert_allclose(lon_out, lon_in, atol=1e-4)
    np.testing.assert_allclose(lat_out, lat_in, atol=1e-4)


# ---------------------------------------------------------------------------
# get_dist_between_two_points_km
# ---------------------------------------------------------------------------

def test_distance_same_point_is_zero():
    """Distance from a point to itself should be zero."""
    d = geospatial_util.get_dist_between_two_points_km(-118.0, 34.0, -118.0, 34.0)
    assert d == 0.0


def test_distance_LA_to_SF_approximate():
    """LA to SF is roughly 560 km. Check we're in the right ballpark."""
    d = geospatial_util.get_dist_between_two_points_km(-118.24, 34.05, -122.45, 37.77)
    assert 530 < d < 600
