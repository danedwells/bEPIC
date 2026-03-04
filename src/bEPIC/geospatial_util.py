#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 11:59:27 2022

@author: amy
"""
import numpy as np
from obspy.geodetics import gps2dist_azimuth
from math import radians, sin, cos, pi


def make_grid(CenterPoint,GridSize=200,GridSpacing=2):

    grid_inc = np.arange(GridSpacing, GridSize + GridSpacing, GridSpacing)
    grid_x = np.hstack((grid_inc[::-1]*-1,0,grid_inc))

    grid_lons,grid_lats = ckm2LLd(grid_x*1000,grid_x*1000,CenterPoint[0],CenterPoint[1],0)

    xx,yy=np.meshgrid(grid_x ,grid_x )
    grid_x_ravel = np.ravel(xx)
    grid_y_ravel = np.ravel(yy)

    return(grid_lons,grid_lats,grid_x,grid_x,grid_x_ravel,grid_y_ravel)




def ckm2LLd(xx,yy,lon0,lat0,rot):
      # convert cartesian to lon lat
      # and rotate if wanted
      R = 6378137
      ff = 1/298.257
      r = R*(1-ff*np.sin(lat0* np.pi / 180.)**2) # r - radius at lat [m]

      mpd = r*np.pi/180
      cos_rot = np.cos(rot* np.pi / 180.);        # cos_rot - cos of rotation angle
      sin_rot = np.sin(rot* np.pi / 180.);        # sin_rot - sin of rotation angle

      if len(xx)==len(yy):
             x_rot = xx*cos_rot + yy*sin_rot
             y_rot =-xx*sin_rot + yy*cos_rot
             # transform from xx,yy to lon,lat using (lon0,lat0) as origin
             lat = lat0 + y_rot/mpd
             lon = lon0 + x_rot/mpd/np.cos(lat0* np.pi / 180.)
      else:
          print('xx and yy are not consistent!')
      return[lon,lat]


def LL2cartd(lon,lat,lon0,lat0,rot):
    """ This program
    (1) transforms lon,lat to local Cartesian coordinate using (lon0,lat0) as origin
    (2) rotates the Cartesian coordinates by rot [degree]
    Input: 
    (1) lon,lat [degree] MUST BE A NUMPY ARRAY
    (2) lon0,lat0 [degree] is used as the origin reference point
    (3) rot [degree] is the rotation angle form the old coord system (east is X+, north is y+)
        to the new one.
        Counterclockwise rotation is positive
        Clockwise rotation is negative
    Output:
    (1) x_rot, y_rot [m] 
    
    first written my Amanda Thomas, 2006
    modified by Lujia Feng, 2008
    translated to python by Amy Williamson, 2018
    last modified by ALW on July 27 2018 """


    R = 6378137          # R - Earth's radius at the equator [m]
    ff = 1/298.257       # ff - flattening factor
    lat0r=radians(lat0)
    r = R*(1-ff*sin(lat0r)**2)     # radius - radius at lat [m]

    mpd = r*pi/180;                # mpd - meters per degree
    cos_rot = cos(radians(rot))    # cos_rot - cos of rotation angle
    sin_rot = sin(radians(rot))    # sin_rot - sin of rotation angle

    if len(lon) == len(lat):
        #transform from lon,lat to xx,yy using lon0,lat0 as origin
       yy = (lat-lat0)*mpd
       xx = (lon-lon0)*mpd*cos(radians(lat0))
       # rotate the coordinate system by rot [degrees]
       x_rot =  xx*cos_rot + yy*sin_rot
       y_rot = -xx*sin_rot + yy*cos_rot
    else:
        print('lon and lat are not consistent')
    LL2cart_out=[x_rot,y_rot]
    return LL2cart_out


def get_dist_between_two_points_km(lon1,lat1,lon2,lat2):

    
    m,az1,az2 = gps2dist_azimuth(lat1, lon1, lat2, lon2)
    d = m/1000
    return(d)

