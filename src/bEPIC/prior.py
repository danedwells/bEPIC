#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 14:02:35 2022

@author: amy
"""
from libcomcat.search import search
from datetime import datetime
import pandas as pd
import os
from scipy.stats import kde
import pandas as pd
import numpy as np
from bEPIC import geospatial_util
import os
from datetime import datetime
#bepic=os.environ['BEPIC']
bepic = "/home/daned/2024_NEHRP/bEPIC"


def generate_prior_seismicity_catalog():
    """
    Args:

    Returns:
    
    """

    region = [-135, -112, 30, 50]
    earthquake = search(starttime=datetime(2000, 1, 1, 0, 0), endtime=datetime.now(),minlatitude=region[2], 
                        maxlatitude=region[3], minlongitude=region[0], maxlongitude=region[1],
                        minmagnitude=3)
    
    columns=['ANSS ID','date','timestamp','lon','lat','depth','mag']
    df = pd.DataFrame(columns=columns)
    for eq in earthquake:
        df.loc[len(df.index)] = [eq.id,str(eq.time),eq.time.timestamp(),eq.longitude,eq.latitude,eq.depth,eq.magnitude]
    

    df.to_csv(bepic+'/data/prior_seismicity_catalog.txt',sep='\t',index=False) 
    
    

def compute_prior(CenterPoint,GridSize,GridSpacing,ANSS_timestamp=None):
    """
    Compute the seismic prior, which is basically just a smoothed map of seismicity

    Args:
    CenterPoint: List of [float, float]
    GridSize: (int)
    GridSpacing: (int)
    ANSS_timestamp: datetime.timestamp()

    Returns:
    prior_seis: np.mgrid - 2D grid of seismicity (smoothed)
    prior_lon: float - longitude of max seismicity
    prior_lat: float - latitude of max seismicity

    """
    if os.path.exists(bepic+'/data/prior_seismicity_catalog.txt')==False:
        print('cannot find catalog... building a new one....')
        generate_prior_seismicity_catalog()
        

    # Generate a spatial grid of lon/lat points centered on CenterPoint
    (grid_lons, grid_lats, _, _,
     _, _) = geospatial_util.make_grid(CenterPoint, GridSize, GridSpacing)

    # Default to current time if no timestamp provided (i.e. not a replay event)
    if ANSS_timestamp is None:
        ANSS_timestamp = datetime.now().timestamp()

    # Load historical seismicity catalog and filter to grid bounds, mag floor, and before target time
    mag_floor = 3
    df = pd.read_csv(bepic+'/data/prior_seismicity_catalog.txt', sep='\t')
    df = df[
        (df['lon'].between(np.min(grid_lons), np.max(grid_lons))) &
        (df['lat'].between(np.min(grid_lats), np.max(grid_lats))) &
        (df['mag'] >= mag_floor) &
        (df['timestamp'] < ANSS_timestamp)
    ].reset_index(drop=True)

    # Convert lon/lat to Cartesian kilometers relative to grid origin
    df['x_km'], df['y_km'] = geospatial_util.LL2cartd(
        df['lon'].values, df['lat'].values, grid_lons[0], grid_lats[0], 0
    )
    df['x_km'] /= 1000
    df['y_km'] /= 1000

    # Convert grid nodes to Cartesian kilometers
    grid_x, grid_y = geospatial_util.LL2cartd(grid_lons, grid_lats, grid_lons[0], grid_lats[0], 0)
    grid_x = grid_x / 1000
    grid_y = grid_y / 1000

    # Fit a Gaussian KDE to the historical earthquake locations
    k = kde.gaussian_kde(df[['x_km', 'y_km']].values.T, bw_method='scott')

    # Isotropize the covariance: force equal spread in all directions
    spread = np.max(k._data_covariance)
    k._data_covariance = np.identity(2) * spread
    k._data_inv_cov    = np.linalg.inv(k._data_covariance)
    k.covariance       = k._data_covariance * k.factor**2
    k.inv_cov          = k._data_inv_cov / k.factor**2

    # Evaluate the KDE on the full grid to produce a 2D seismicity density map
    xi, yi = np.mgrid[np.min(grid_x):np.max(grid_x):len(grid_x)*1j,
                      np.min(grid_y):np.max(grid_y):len(grid_y)*1j]
    zi = k(np.vstack([xi.flatten(), yi.flatten()]))
    prior_seis = zi.reshape(xi.shape).T

    # Find the grid point with the highest prior seismicity density
    best_location_prior = np.where(prior_seis == np.max(prior_seis))
    prior_lon = grid_lons[best_location_prior[1][0]]
    prior_lat = grid_lats[best_location_prior[0][0]]\
    
    # Return the 2D prior seismicity grid and the coordinates of its peak
    return (prior_seis, prior_lon, prior_lat)





