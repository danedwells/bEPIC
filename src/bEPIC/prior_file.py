#!/usr/bin/env python3
"""
prior_file.py — bEPIC adapter for SeismicPrior objects.

Provides compute_prior_from_model(), which has the same call signature and
return values as prior.compute_prior() so it can be used as a drop-in
replacement inside the run_bEPIC version loop.

The prior is bilinearly interpolated from the pre-built SeismicPrior grid
onto the dynamic bEPIC spatial grid centered on CenterPoint.  Grid points
that fall outside the prior's geographic extent are assigned zero probability.
If the interpolated region has no coverage at all (e.g. the event is far
outside California), the prior falls back to a uniform distribution.
"""

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from bEPIC import geospatial_util


def compute_prior_from_model(CenterPoint, GridSize, GridSpacing, prior_model):
    """
    Evaluate a SeismicPrior on the bEPIC spatial grid.

    Parameters
    ----------
    CenterPoint : list[float, float]
        [longitude, latitude] of the grid center.
    GridSize : float
        Grid radius in km.
    GridSpacing : float
        Grid resolution in km.
    prior_model : SeismicPrior
        Pre-built prior (from priors.prior_model.SeismicPrior).

    Returns
    -------
    prior_grid : np.ndarray
        2-D probability array on the bEPIC grid, shape (n, n), normalized
        to sum to 1.  prior_grid[row, col] corresponds to
        (grid_lats[row], grid_lons[col]).
    prior_lon : float
        Longitude of the prior's peak (MAP estimate).
    prior_lat : float
        Latitude of the prior's peak (MAP estimate).
    """
    grid_lons, grid_lats, _, _, _, _ = geospatial_util.make_grid(
        CenterPoint, GridSize, GridSpacing
    )

    # Build a bilinear interpolator over the pre-computed prior grid.
    # prior_model.grid has shape (len(lats), len(lons)), matching the
    # (points, points) convention expected by RegularGridInterpolator.
    interp = RegularGridInterpolator(
        (prior_model.lats, prior_model.lons),
        prior_model.grid,
        method='linear',
        bounds_error=False,
        fill_value=0.0,
    )

    # Construct the bEPIC lon/lat meshgrid and query the interpolator.
    # Result shape: (len(grid_lats), len(grid_lons)) = (n, n).
    LON, LAT = np.meshgrid(grid_lons, grid_lats)
    pts        = np.column_stack([LAT.ravel(), LON.ravel()]) # Technically hstack *could* work here too
    prior_grid = interp(pts).reshape(LON.shape)

    # Normalize; fall back to uniform if the prior has no coverage here.
    total = prior_grid.sum()
    if total > 0:
        prior_grid = prior_grid / total
    else:
        prior_grid = np.ones_like(prior_grid) / prior_grid.size

    best      = np.where(prior_grid == np.max(prior_grid))
    prior_lon = grid_lons[best[1][0]]
    prior_lat = grid_lats[best[0][0]]

    return prior_grid, prior_lon, prior_lat
