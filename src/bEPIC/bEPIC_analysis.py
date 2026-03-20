#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 17:35:41 2022

@author: amy
"""
from bEPIC import geospatial_util
import pandas as pd
import numpy as np
import os



def compute_station_trigger_misfit(postgres_id,project_parent_directory):

    
    
   
    if os.path.exists(project_parent_directory+postgres_id+'/USGS/'+'usgs_event_summary.txt')==True:
        # open USGS file
        usgs_df = pd.read_csv(project_parent_directory+postgres_id+'/USGS/'+'usgs_event_summary.txt',sep='\t')
        
        # open run file
    
        run_df = pd.read_csv(project_parent_directory+postgres_id+'/'+postgres_id+'.run')
            
        last_version = np.unique(run_df['version'])[-1]
        idx = np.where(run_df['version'] ==last_version  )[0]
        run_df = run_df.iloc[idx].reset_index(drop=True)
    
        used_in_loc = np.zeros(len(run_df))
        idx = np.where(run_df['tterr']>-999)[0]
        used_in_loc[idx]=1
    
        CenterPoint=[usgs_df['USGS lon'].iloc[0],usgs_df['USGS lat'].iloc[0]]
        stax,stay =  geospatial_util.LL2cartd(np.array(run_df['longitude']),
                                   np.array(run_df['latitude']),
                                   CenterPoint[0],CenterPoint[1],0)
        station_distance = np.sqrt( (stax/1000)**2 + (stay/1000 )**2+ usgs_df['USGS depth'].iloc[0]**2)
        travel_time = station_distance / 6.0 # seconds
        usgs_predicted_time = usgs_df['USGS time'].iloc[0]  + travel_time 
        time_offset = run_df['trigger time']-usgs_predicted_time
    
    
    
        misift_df = pd.DataFrame({'station':run_df['station'],
                                  'network':run_df['network'],
                                  'channel':run_df['channel'],
                                  'station lon':run_df['longitude'],
                                  'station lat':run_df['latitude'],
                                  'trigger offset':time_offset,
                                  'station distance':station_distance,
                                  'station in location':used_in_loc})   
    
        misift_df.to_csv(project_parent_directory+postgres_id+'/USGS/'+'usgs_trigger_time_misfit.txt',sep='\t',index=False)
    
        
        
    else:
        print('no USGS file!!!')
        
        
