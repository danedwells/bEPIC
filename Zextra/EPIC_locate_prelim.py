#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

This script consolidates older versions of bEPIC into one script. It also best
mirrors what is being called on the ShakeAlert EPIC version of the code (E2Locate)

What this code needs as inputs:
    
    (1) an event object that contains the initial epicenter guess (lon, lat, origin time)
    
    event = Event(lat = 37,          <-----epicenter lat
                  lon = -122,        <-----epicenter lon
                  time = 100,        <-----dummy origin time (timestamp)
                  misfit_rms = 0,    <------ initial misfit (dummy var)
                  misfit_ave = 0,     <------ initial misfit (dummy var)
                  eventid = 100,     <------ event ID (dummy for keeping track of events)
                  version = 0)       <------ version ( for keeping track of events)

    (2) what stations trigger and when


    t = TriggerManager(lon = station_lon, lat = station_lat, sta=station_name, net=station_network, chan=station_channel,trigger_time = station_trigger_timestamp)
    event.trigs.append(t)    #(for each trigger, you add it to the event)

    

    (3) event replay parameters:
    
    params = EPIC_PARAMS()
    params.PriorGridFile =  '/Users/amy/projects/container_bEPIC/data/prior_seis_grid_US_Canada.tt3'
    params.use_prior = True   # can toggle True or False. False JUST uses the misfit grid, acting like old EPIC
    params.GridSize = 25  
    params.GridKm = 50
    params.method = 'EPIC C'     # keep this as EPIC C
    
    
    
    # GridSize and GridKm are grandfathered the nomenclature of EPIC, which is a bit wonky.
    # EPIC uses a square grid. GridKm is the distance from the center to the edge. So if
    # GridKm = 50, the dimensions of the grid search are 100 by 100 km. 
    # The grid spacing is 50/25 = 2, meaning that each grid node is 2 km apart.
    #
    # I know, it is unnecessarily confusing.
    
    
    
    # call the event using the params and event objects
    
    t,output_df = E2Location_locate(params,event)


    # ------------------------------------------------------
    # OUTPUTS
    
    t contains the best grid location, prior value, misfit location, etc etc
    
    output_df is a dataframe that is grid n x grid m in length containing the output
    information (misfit, likelihood, prior, post) for all grid nodes
    
    



Created on Fri Oct 31 13:24:40 2025
@author: amy
"""
import numpy as np
from scipy import interpolate
import pandas as pd

def latLonToXY(event):
    lat0 =event.lat
    lon0 = event.lon
    
    
    R = 6378137;
    ff = 1./298.257                            #// flattening factor
    lat0r = lat0*np.pi/180.
    r = R*(1 - ff*np.power(np.sin(lat0r), 2))  #// radius - radius at lat [m]
    mpd = r*np.pi/180.;                        #// mpd - meters per degree
    
    # // get the station coordinates on the local grid centered at lon0, lat0
    num_trigs = len(event.trigs)
    
    for i in range(num_trigs):
        event.trigs[i].stay = (event.trigs[i].lat-lat0)*mpd/1000
    
    f = mpd*np.cos(lat0r)/1000
    for i in range(num_trigs):
        event.trigs[i].stax = (event.trigs[i].lon-lon0)*f
    return(event)



class Event:  
    def __init__(self,lat,lon,time,misfit_rms,misfit_ave, eventid, version):
        self.lat = lat
        self.lon = lon
        self.time = time
        self.depth = 8.0
        self.misfit_rms = misfit_rms
        self.misfit_ave = misfit_ave
        self.eventid = eventid
        self.version = version
        
        self.trigs=[]


class LocInfo:  
    def __init__(self,initial_lat, initial_lon, initial_depth, initial_time):
        self.initial_lat = initial_lat
        self.initial_lon = initial_lon
        self.initial_time = initial_time
        self.initial_depth = 8.0

        
class TriggerManager:
    
    def __init__(self,lon,lat,trigger_time,sta,net,chan):
        self.lon          = lon
        self.lat          = lat
        self.time = trigger_time
        self.sta          = sta
        self.net          = net
        self.chan         = chan
        self.dist         = np.nan
        self.tt           = np.nan
        self.tterror      = np.nan

        


class PriorFile:
    def __init__(self,PriorGridFile):
        
        fid = open(PriorGridFile)
        self.mx     = int(fid.readline().split()[0]);       
        self.my     = int(fid.readline().split()[0])
        self.xlower = float(fid.readline().split()[0]);      
        self.ylower = float(fid.readline().split()[0])
        self.dx     = float(fid.readline().split()[0]);      
        self.dy     = float(fid.readline().split()[0])
        
        fid.close()
 
        
        prior_prior = np.loadtxt(PriorGridFile , skiprows=6, dtype = np.float64)
        self.prior = (np.flipud(prior_prior)).flatten()
        
        

class SearchOut:
    def __init__(self):
        self.best_location_post = 0
        self.posterior_lon = np.nan
        self.posterior_lat = np.nan
        self.best_misfit = 0
        self.misfit_ave = 0
        self.best_OT = np.nan
        self.best_grid_x = np.nan
        self.best_grid_y = np.nan
        self.best_value = 0
        self.best_like = 0
        self.best_prior = 0
        


class EPIC_PARAMS:
    def __init__(self):
        self.MAX_EVENT_TRIGS = 100
        self.LocationPVelocity = 6.0
        
        

def get_dist_between_two_points_km(lon1,lat1,lon2,lat2):
    from obspy.geodetics import gps2dist_azimuth
    
    m,az1,az2 = gps2dist_azimuth(lat1, lon1, lat2, lon2)
    d = m/1000
    return(d)
   
##########################################################################################


def E2Location_locate(params,event):
    


    #// Initialize search variables
    evlat  = event.lat                   # I think this is just to save the old location
    evlon  = event.lon
    evtime = event.time


   
    #// set new event location 
    event.lat = evlat
    event.lon = evlon
    event.time = evtime
    event.misfit_rms = 0
    event.misfit_ave = 0

    # ---------------------------------------------------
    # // MULTI STATION LOCATION

    #for all triggers in trigger object
    #// Find event and station coordinates on x-y grid (flat earth).
    event = latLonToXY(event)
    
    trigs = event.trigs
    
    
    t,output_df = E2Location_searchGrid(event,trigs, params)
    
    evlon = t.posterior_lon
    evlat = t.posterior_lat
    
    num_trigs = len(trigs)
    for i in range(num_trigs):
        trigs[i].dist = get_dist_between_two_points_km(evlon,evlat,trigs[i].lon,trigs[i].lat)
        trigs[i].tt = trigs[i].time - t.best_OT
    
    
    ddist = get_dist_between_two_points_km(event.lon,event.lat,evlon,evlat)
    
    print("L:E:H: \t eventid \t lat0 \t lon0 \t lat \t lon \t ddist \t avefit \t rmsfit \t nT")
    print("L:E:   \t"+str(event.eventid)
          +'\t'+str(np.round(event.lat,4))
          +'\t '+str(np.round(event.lon,4))
          +'\t'+str(np.round(evlat,4))
          +'\t'+str(np.round(evlon,4))
          +'\t'+str(np.round(ddist,4))
          +'\t'+str(np.round(t.misfit_ave,4))
          +'\t'+str(np.round(t.best_misfit,4))
          +'\t'+str(num_trigs))
    print('# --------------------------------------------------------')
    print("L:T:H: \t eventid \t nT \t sta \t chan \t net \t lat \t lon \t dist \t tt")
    for i in range(num_trigs):
        print('L:T:  \t'+str(event.eventid)
              +'\t'+str(num_trigs)
              +'\t'+trigs[i].sta +'\t'+trigs[i].chan  +'\t'+trigs[i].net
              +'\t'+str(np.round(trigs[i].lat,4))
              +'\t'+str(np.round(trigs[i].lon,4))
              +'\t'+str(np.round(trigs[i].dist,2))
              +'\t'+str(np.round(trigs[i].tt,2)))
    print('# --------------------------------------------------------')
    
    return(t,output_df)
        





def E2Location_searchGrid(event, trigs, params):
    
    # replicate searchGrid
    # running a full replicate is really slow in Python- python does really well at
    # vectorized functions. C++ is doing a lot of looping here over multiple threads
    
    
    #// Initialize search variables
    evlat  = event.lat                   # I think this is just to save the old location
    evlon  = event.lon
    evtime = event.time
    
   
    num_trigs = len(trigs)

    trig_ot = np.zeros(num_trigs)

    
   #// The location grid is a square with (2*GridSize + 1) grid-points on each side
   #// The grid point separation is (GridKm / GridSize)
    grid_size       = params.GridSize            # // 25
    grid_width      = 2*grid_size + 1;
    grid_km         = params.GridKm              #  // 50.
    grid_spacing_km = grid_km/grid_size;          # // 2.0
    p_velocity      = params.LocationPVelocity   # // 6.0
   
    t = SearchOut()
    prior_info = PriorFile(params.PriorGridFile)

   
    # velocity model
    vel_mod_filename = '/Users/amy/projects/EPICdb/testscripts/h2p+ak135.080'
    tt_mod = np.genfromtxt(vel_mod_filename, skip_header=1)
    ttf = interpolate.interp1d(tt_mod[:,0], tt_mod[:,1])


    lat0  = evlat;       
    lon0  = evlon
    R     = 6378.137;               
    ff    = 1./298.257                                # // flattening factor
    lat0r = lat0*np.pi/180.;        
    r     = R*(1 - ff*np.power(np.sin(lat0r), 2));    # // radius - radius at lat [m]
    mpd   = r*np.pi/180.;           
    f     = mpd*np.cos(lat0r)                         # // mpd - meters per degree

    j = 0   #// lat index in the prior grid
    i = 0   #// lon index in the prior grid


    ybeg = -1 *grid_size
    yend = grid_size
    grid_y =grid_x = np.linspace(ybeg,yend,2*grid_size + 1)


    
    

    if params.method == 'python bypass':
        
        ''' do not use the bypass it is being phased out'''
        
        
        n = len(trigs)              # number of active stations
        p = len(grid_x)                    # number of grid nodes in x (or y) direction
        m = p*p                            # total number of grid nodes  (p*p)
       
        xx,yy=np.meshgrid(grid_x ,grid_y )
    
        
        prior_val = []
        for y in grid_y:
            
            ykm = y*grid_spacing_km
            ylat = lat0 + ykm/mpd 
            
            if params.use_prior == True:
                a = (ylat - prior_info.ylower)/prior_info.dy
                j = (int)(a+0.5)
                if(j<0): j=0;
                if (j >= prior_info.my): j = prior_info.my - 1
                
                
            for x in grid_x:
                xkm = x*grid_spacing_km
                xlon = lon0 +xkm/f
                if params.use_prior == True:
                    a = (xlon - prior_info.xlower)/prior_info.dx
                    i = (int)(a+0.5)
                    if (i < 0): i = 0;
                    if (i >= prior_info.mx): i = prior_info.mx - 1
                
                prior_val = np.append(prior_val,prior_info.prior[j*prior_info.mx+i])
    
    
        
    
    
    
    
    elif params.method == 'EPIC C':
        
        # will want to save all values into dataframe
        output_df = pd.DataFrame(columns=['y','x','lat','lon','like','prior','post','misfitrms','misfitave'])
        
        
        print('grid spacing: '+str(grid_spacing_km))
        for y in grid_y:
            
            ykm = y*grid_spacing_km
            ylat = lat0 + ykm/mpd 
            
            if params.use_prior == True:
                a = (ylat - prior_info.ylower)/prior_info.dy
                j = (int)(a+0.5)
                if(j<0): j=0;
                if (j >= prior_info.my): j = prior_info.my - 1
                
                
            
            for x in grid_x:
                xkm = x*grid_spacing_km
                sumOT = 0
                
                for it in range(num_trigs):
                    tx   = trigs[it].stax - xkm
                    ty   = trigs[it].stay - ykm
                    dist = np.sqrt(tx*tx + ty*ty)
                    
                    #// d = distance from i'th station to grid point x,y
                    #// stt = travel time from this grid point to station
                    #// find distance in velocity model
                    
                    stt = ttf(dist)
                    
                    #// event origin time for this station is trigger_time - travel_time
                    trig_ot[it] = trigs[it].time - stt;
                    sumOT += trig_ot[it];
                #// Average origin time for this grid point
                aveOT = sumOT/num_trigs
                #// compute rms of trigger travel time errors
                rms   = 0
                ttsum = 0
                for it in range(num_trigs):
                    rms   += np.power(trig_ot[it]-aveOT,2)
                    ttsum += np.fabs(trig_ot[it]-aveOT)
    
                misfitsq = rms/num_trigs
                misfit_ave = ttsum/num_trigs
                
                like = 1
                for it in range(num_trigs):
                    tterror = trig_ot[it] - aveOT
                    d = tterror*tterror
                    
                    e = np.exp(-0.5*d)
                    like *= e
                    
                like = np.sqrt(like/num_trigs)
                post = like
                Prior = 1
                
                xlon = lon0 +xkm/f
                if params.use_prior == True:
                    a = (xlon - prior_info.xlower)/prior_info.dx
                    i = (int)(a+0.5)
                    if (i < 0): i = 0;
                    if (i >= prior_info.mx): i = prior_info.mx - 1
                    
                    Prior = prior_info.prior[j*prior_info.mx+i]
                    post = like*Prior
                
                
                output_df.loc[len(output_df)] = [y,x,ylat,xlon,like,Prior,post,misfitsq,misfit_ave]
                if post > t.best_location_post:
                    t.best_location_post = post
                    t.posterior_lon = xlon
                    t.posterior_lat = ylat
                    t.best_misfit = misfitsq
                    t.misfit_ave = misfit_ave
                    t.best_OT = aveOT
                    t.best_grid_x = x
                    t.best_grid_y = y
                    t.best_value = post
                    t.best_like = like
                    t.best_prior = Prior
                    
    
    print("best_location_post: "+str(np.round(t.best_location_post, 12)))
    print("posterior_lon: "+str(np.round(t.posterior_lon, 6)))
    print("posterior_lat: "+str(np.round(t.posterior_lat, 6)))
    print("bestOT: "+str(np.round(t.best_OT, 6)))        
    print("best_misfit: "+str(np.round(t.best_misfit, 6)))    
    
    print("best_grid_x: "+str(int(t.best_grid_x)))    
    print("best_grid_y: "+str(int(t.best_grid_y)))    
    
    print("best_value: "+str(np.round(t.best_value, 12)))   
    print("best_like: "+str(np.round(t.best_like, 12)))   
    print("best_prior: "+str(np.round(t.best_prior, 12)))   
    print('#################################################################')
      
    return(t,output_df)
                    
                






