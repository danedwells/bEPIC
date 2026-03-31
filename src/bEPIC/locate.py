#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 11:55:10 2022

@author: amy
"""
import numpy as np
import pandas as pd

def get_two_station_location(sta_df):

    station_01_idx = np.where(sta_df['order']==1)[0][0]
    station_01_lon = sta_df['longitude'].iloc[station_01_idx]
    station_01_lat = sta_df['latitude'].iloc[station_01_idx]
    station_01_OT = sta_df['trigger time'].iloc[station_01_idx]

    station_02_idx = np.where(sta_df['order']==2)[0][0]
    station_02_lon = sta_df['longitude'].iloc[station_02_idx]
    station_02_lat = sta_df['latitude'].iloc[station_02_idx]
    station_02_OT = sta_df['trigger time'].iloc[station_02_idx]

    if station_01_OT <= station_02_OT:
        eq_lat = (station_01_lat*2 + station_02_lat)/3
        eq_lon = (station_01_lon*2 + station_02_lon)/3
    elif station_02_OT <= station_01_OT:
        eq_lat = (station_02_lat*2 + station_01_lat)/3
        eq_lon = (station_02_lon*2 + station_01_lon)/3

    grid_center = [ eq_lon, eq_lat]
    return(grid_center)