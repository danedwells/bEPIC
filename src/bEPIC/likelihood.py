#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 13:34:37 2022

@author: amy
"""
import numpy as np
from bEPIC import geospatial_util,data_util
  

def calculate_likelihood(CenterPoint, sta_df, velocity_model, GridSize=200, GridSpacing=2):
    """
    Computes the likelihood surface L(m|d) for earthquake location.

    For each node on a spatial grid, computes the likelihood that the earthquake
    originated there, given the observed P-wave trigger times at each station.
    Uses a Gaussian likelihood: L = exp(-0.5 * misfit / sigma^2)

    Args:
        CenterPoint (list):     [lon, lat] center of the search grid.
        sta_df (DataFrame):     Active stations with columns: longitude, latitude, trigger time, sigma.
        velocity_model (str):   Either 'h2p+ak135' (tabulated) or 'constant' (6.0 km/s).
        GridSize (float):       Spatial extent of the search grid in km.
        GridSpacing (float):    Grid resolution in km.

    Returns:
        like (2D array):        Likelihood surface (p x p).
        misfit (2D array):      Sum of squared travel time residuals (p x p).
        likelihood_lon (float): Longitude of the peak likelihood location.
        likelihood_lat (float): Latitude of the peak likelihood location.
    """

    (grid_lons, grid_lats, grid_x, grid_y,
     grid_x_ravel, grid_y_ravel) = geospatial_util.make_grid(CenterPoint, GridSize, GridSpacing)

    eq_depth = 8.0  # assumed earthquake depth in km

    # Grid and station dimensions
    n = len(sta_df['longitude'])  # number of active stations
    p = len(grid_x)               # number of grid nodes along one axis
    m = len(grid_x_ravel)         # total grid nodes (p x p)

    # Convert station lon/lat to Cartesian coordinates (km) relative to CenterPoint
    stax, stay = geospatial_util.LL2cartd(
        np.array(sta_df['longitude']),
        np.array(sta_df['latitude']),
        CenterPoint[0], CenterPoint[1], 0
    )
    stax_km = stax / 1000
    stay_km = stay / 1000

    # Compute straight-line distance (km) from each grid node to each station
    # Includes a fixed depth term; result is an (m x n) matrix
    station_distance = np.sqrt(
        np.square(np.subtract(np.tile(grid_x_ravel, (n, 1)).T, np.tile(stax_km, (m, 1)))) +
        np.square(np.subtract(np.tile(grid_y_ravel, (n, 1)).T, np.tile(stay_km, (m, 1)))) +
        eq_depth**2
    )

    # Compute P-wave travel time from each grid node to each station — (m x n) matrix
    if velocity_model == 'h2p+ak135':
        # Use a pre-tabulated travel time function interpolated over distance
        ttf = data_util.travel_time_function(velocity_model)
        travel_time = ttf(np.ravel(station_distance)).reshape(m, n)
    elif velocity_model == 'constant':
        # Simple constant velocity approximation (6.0 km/s)
        travel_time = station_distance / 6.0

    # Estimate origin time at each grid node as the mean of (trigger time - travel time)
    # across all stations; result is a vector of length m
    average_OT = np.mean(
        np.subtract(np.tile(sta_df['trigger time'], (m, 1)), travel_time),
        axis=1
    )

    # Forward-model the expected trigger time at each station from each grid node
    # result is an (m x n) matrix
    trigger_time_calc = np.add(np.tile(average_OT, (n, 1)).T, travel_time)

    # Compute squared residuals between observed and modeled trigger times
    # result is an (m x n) matrix
    tt_error = np.abs(np.subtract(trigger_time_calc, np.tile(sta_df['trigger time'], (m, 1))))**2

    # Sum residuals across stations to get total misfit at each grid node — (p x p) matrix
    misfit = np.sum(tt_error, axis=1).reshape(p, p)

    # Compute per-station Gaussian likelihood: exp(-0.5 * residual / sigma^2) — (m x n) matrix
    rho = np.exp(-0.5 * (tt_error / (sta_df['sigma'].values**2)))

    # Combined likelihood at each grid node: product across all stations — (p x p) matrix
    like = np.prod(rho, axis=1).reshape(p, p)

    # Extract the grid coordinates of the peak likelihood
    best_location  = np.where(like == np.max(like))
    likelihood_lon = grid_lons[best_location[1][0]]
    likelihood_lat = grid_lats[best_location[0][0]]

    return (like, misfit, likelihood_lon, likelihood_lat)
