#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 13:44:43 2022

@author: amy
"""


def compute_magnitude(run_df,version,epicenter):
    #-------------------------------------------------------------------------#
    #    Scaling relation 
    #        M = a*log10(Pd)+ b*log10(R)+ c
    #
    #  a = 1.23, b = 1.39, c = 5.39 (from Chung et al. 2020 BSSA)
    #
    #-------------------------------------------------------------------------#
    from bEPIC import geospatial_util
    import numpy as np
    
    
    idx    =   np.where(run_df['version']==version)[0]  
    sta_df =   run_df.iloc[idx].reset_index(drop=True)
    
    
    
    a = 1.23
    b=1.39
    c=5.39
    
    # get distance from each station

    
    num_stations = len(sta_df)
    station_distance = np.zeros(num_stations)
    for j in range(num_stations):
        station_distance[j]= geospatial_util.get_dist_between_two_points_km(epicenter[0],epicenter[1],
                                                                 sta_df['longitude'].iloc[j],sta_df['latitude'].iloc[j])
        if station_distance[j]==0:
            station_distance[j]=np.nan
            
        if station_distance[j]>200:
            station_distance[j]=np.nan
            
    Pd = 10**(sta_df['logPd'].values)
    
    if 0 in Pd:
        print(Pd)
        Pd[np.argwhere(Pd==0)[0][0]] = np.nan
        
        
    if 0 in   station_distance:
        print(station_distance)
        station_distance[np.argwhere(station_distance==0)[0][0]] = np.nan
        
    m = a*np.log10(Pd) + b*np.log10(station_distance) + c 
    mag = np.nanmean(m)
    return(mag)
