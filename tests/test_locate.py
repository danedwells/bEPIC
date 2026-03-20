import numpy as np
import pandas as pd
import pytest
from bEPIC import locate


def make_sta_df(lon1, lat1, t1, lon2, lat2, t2):
    """Helper: build a minimal sta_df with two stations."""
    return pd.DataFrame({
        'order':        [1, 2],
        'longitude':    [lon1, lon2],
        'latitude':     [lat1, lat2],
        'trigger time': [t1, t2],
    })


def test_first_station_triggers_first_weighted_toward_it():
    """When station 1 triggers first, the center point should be closer to it."""
    sta_df = make_sta_df(
        lon1=-120.0, lat1=36.0, t1=100.0,   # station 1 triggers first
        lon2=-118.0, lat2=34.0, t2=110.0,
    )
    center = locate.get_two_station_location(sta_df)
    # Expected: (2*sta1 + sta2) / 3
    expected_lon = (2 * -120.0 + -118.0) / 3
    expected_lat = (2 * 36.0   +  34.0)  / 3
    assert abs(center[0] - expected_lon) < 1e-9
    assert abs(center[1] - expected_lat) < 1e-9


def test_second_station_triggers_first_weighted_toward_it():
    """When station 2 triggers first, the center point should be closer to it."""
    sta_df = make_sta_df(
        lon1=-120.0, lat1=36.0, t1=110.0,
        lon2=-118.0, lat2=34.0, t2=100.0,   # station 2 triggers first
    )
    center = locate.get_two_station_location(sta_df)
    expected_lon = (2 * -118.0 + -120.0) / 3
    expected_lat = (2 *  34.0  +  36.0)  / 3
    assert abs(center[0] - expected_lon) < 1e-9
    assert abs(center[1] - expected_lat) < 1e-9


def test_simultaneous_triggers_uses_station_1_weighting():
    """Equal trigger times uses the station_01_OT <= station_02_OT branch."""
    sta_df = make_sta_df(
        lon1=-120.0, lat1=36.0, t1=100.0,
        lon2=-118.0, lat2=34.0, t2=100.0,
    )
    center = locate.get_two_station_location(sta_df)
    expected_lon = (2 * -120.0 + -118.0) / 3
    expected_lat = (2 *  36.0  +  34.0)  / 3
    assert abs(center[0] - expected_lon) < 1e-9
    assert abs(center[1] - expected_lat) < 1e-9


def test_output_is_list_of_two_floats():
    sta_df = make_sta_df(-120.0, 36.0, 100.0, -118.0, 34.0, 110.0)
    center = locate.get_two_station_location(sta_df)
    assert len(center) == 2
