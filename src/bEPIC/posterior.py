#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 14:57:47 2022

@author: amy
"""


def compute_posterior(CenterPoint,GridSize,GridSpacing,prior_function,likelihood_function):
    """
    Computes the posterior location function with appropriate normalization demoninator (k)
    Uses np.trapz to integrate over x and y coordinates

    Args:
    ---------------
    CenterPoint: 
    GridSize
    GridSpacing
    prior_function
    likelihood_function

    Returns:
    ---------------
    post
    posterior_lon
    posterior_lat
    
    """

    from bEPIC import geospatial_util
    import numpy as np
    
    (grid_lons,grid_lats,grid_x,grid_y,
     grid_x_ravel,grid_y_ravel) = geospatial_util.make_grid(CenterPoint,GridSize,GridSpacing)

    k = np.trapezoid(np.trapezoid(prior_function*likelihood_function,grid_x,axis=0),grid_y)
    post = (1/k)*likelihood_function*prior_function
     
     
    best_location_post  =  np.where(post == np.max(post))
    posterior_lon   =  grid_lons[best_location_post[1][0]]
    posterior_lat   =  grid_lats[best_location_post[0][0]]
    
    
    return(post,posterior_lon,posterior_lat)