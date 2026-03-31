#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 11:18:41 2022

@author: amy
"""
import pandas as pd        
import numpy as np
from datetime import datetime,timedelta
from obspy.geodetics import gps2dist_azimuth

       
def generate_run_file(project_parent_directory,postgres_id):

    station_trigger_log_file = project_parent_directory+postgres_id+'/EPIC/'+postgres_id+'_event_triggers_log.txt'
    output_filename =  project_parent_directory+postgres_id+'/'+postgres_id+'.run'
    
    
    #-----------------------------------------------------------------------------#
    v =  np.genfromtxt(station_trigger_log_file,usecols=[1],skip_header=1)
    o =  np.genfromtxt(station_trigger_log_file,usecols=[3],skip_header=1)

    station_sta   =  np.genfromtxt(station_trigger_log_file,usecols=[4], dtype='str',skip_header=1)
    station_chan  =  np.genfromtxt(station_trigger_log_file,usecols=[5], dtype='str',skip_header=1)
    station_net   =  np.genfromtxt(station_trigger_log_file,usecols=[6], dtype='str',skip_header=1)
    station_loc   =  np.genfromtxt(station_trigger_log_file,usecols=[7], dtype='str',skip_header=1)
    station_t_str =  np.genfromtxt(station_trigger_log_file,usecols=[10],dtype='str',skip_header=1)
    station_lon   =  np.genfromtxt(station_trigger_log_file,usecols=[9],skip_header=1)
    station_lat   =  np.genfromtxt(station_trigger_log_file,usecols=[8],skip_header=1)
    station_t_str =  np.genfromtxt(station_trigger_log_file,usecols=[10],dtype='str',skip_header=1)
    station_t     = []
    version       = []
    order         = []


    for i in range(len(station_t_str)):
        station_t = np.append(station_t,datetime.strptime(station_t_str[i], '%Y-%m-%dT%H:%M:%S.%f').timestamp())
        version   = np.append(version,int(v[i]))
        order     = np.append(order,int(o[i]))
    station_logpd = np.genfromtxt(station_trigger_log_file,usecols=[16],skip_header=1)
    station_tterr = np.genfromtxt(station_trigger_log_file,usecols=[34],skip_header=1)

    #--------------------------------------------------#
    df = pd.DataFrame({'version':version,
                       'order':order,
                       'station':station_sta,
                       'channel':station_chan,
                       'network':station_net,
                       'location':station_loc,
                       'longitude':station_lon,
                       'latitude':station_lat,
                       'trigger time':station_t,
                       'tterr':station_tterr,
                       'logPd':station_logpd})
    df = df.astype({'version':'int'})
    df = df.astype({'order':'int'})
    df.to_csv(output_filename,index=False)
    #--------------------------------------------------#    
    
def search_for_USGS_event(project_parent_directory,postgres_id):
    
    
    event_summary_df = pd.read_csv(project_parent_directory+postgres_id+'/EPIC/'+postgres_id+'_event_summary_log.txt',sep='\t')
    EPIC_ot = datetime.strptime(event_summary_df['time'].iloc[-1], '%Y-%m-%dT%H:%M:%S.%f')
    dt=0
    dx=0
    found=False

    
    while found == False:
            starttime = (EPIC_ot-timedelta(seconds=dt))
            endtime   = (EPIC_ot+timedelta(seconds=dt))
            maxradiuskm=dx
            latitude=event_summary_df['event lat'].iloc[-1]
            longitude=event_summary_df['event lon'].iloc[-1]
            jdict = search_comcat(starttime,endtime,maxradiuskm,latitude,longitude)
            
            if len(jdict["features"]) ==1:
            
                print('found event with dx: ',dx,' and dt: ',dt)
                USGS_event_id =jdict["features"][0]['id']
                USGS_event_latitude=jdict["features"][0]['geometry']['coordinates'][1]
                USGS_event_longitude=jdict["features"][0]['geometry']['coordinates'][0]
                USGS_event_depth = jdict["features"][0]['geometry']['coordinates'][2]
                USGS_event_magnitude = jdict["features"][0]['properties']['mag']
                USGS_event_time=jdict["features"][0]['properties']['time']
                
                m,az1,az2 = gps2dist_azimuth(USGS_event_latitude, USGS_event_longitude, event_summary_df['event lat'].iloc[-1], event_summary_df['event lon'].iloc[-1])
                
                
                
                catalog_df = pd.DataFrame({'postgres id':postgres_id,'USGS ID':USGS_event_id,'USGS time':USGS_event_time,
                                                'USGS lat':USGS_event_latitude,'USGS lon':USGS_event_longitude,
                                                'USGS depth':USGS_event_depth,'USGS mag':USGS_event_magnitude},index=[0])   
                
                
                
                found = True
            elif len(jdict["features"])>1:
                # multiple events found
                dt=dt-1
                dx=dx-1
            else:
                dt=dt+5
                dx=dx+5
    
            if dx > 100:
                    print('this is getting out of hand, need a human to help')
    
                    found = True
                    catalog_df = pd.DataFrame({'postgres id':postgres_id,'USGS ID':np.nan,'USGS time':np.nan,
                                                    'USGS lat':np.nan,'USGS lon':np.nan,
                                                    'USGS depth':np.nan,'USGS mag':np.nan},index=[0])
    
    catalog_df.to_csv(project_parent_directory+postgres_id+'/USGS/'+'usgs_event_summary.txt',sep='\t',index=False)



def travel_time_function(velocity_model):
    from scipy import interpolate
    import numpy as np
    import os
    
    
    bepic = os.path.dirname(os.path.abspath(__file__))
    tt_file = bepic+'/data/h2p+ak135.080'
    
    tt_mod = np.genfromtxt(tt_file,skip_header=1)
    tt_mod_distance = tt_mod[:,0]
    tt_mod_time = tt_mod[:,1]
    ttf = interpolate.interp1d(tt_mod_distance, tt_mod_time)

    return(ttf)   


def parse_log(project_parent_directory,log_file,event_id,epic_id):
    import os
    import pandas as pd

    
    bepic = os.path.dirname(os.path.abspath(__file__))
    
    
    print('paring event long for epic event id '+epic_id)
    ee=0
    #-----------------------------------------------------------------------------#
    file1 = open(log_file, 'r')
    Lines = file1.readlines()


    event_summary_df = pd.DataFrame(columns=['eventid','version','event lat','event lon','depth','mag','time','latu','lonu','depu',
               'magu','timeu','lk','nTb','nSb','nT','nS','ave','rms','fitOK','splitOK','near',
               'statrig','active','inact','nsta','percent','prcntOK','mindist','maxdist','distOK',
               'azspan','MOK','nSOK','LOK','Tdif','tpave','pdave','TF','TOK','AZOK','AOK','Ast','alert_time'])


    event_triggers_df = pd.DataFrame(columns=['eventid','version','update','order','sta','chan','net','loc','lat','lon','trigger time','rsmp','tsmp',
               'log taup','taup snr','dsmp','log pd','pd snr','assoc','tpmag','utpm','pdmag','updm','uch','ukm','upd',
               'ups','utp','uts','tel','tsec','distkm','azimuth','TF','tterr','azerror','incid','plen','sps'])



    location_triggers_df = pd.DataFrame(columns=['eventid','version','nT','index',
                                                 'sta','chan','net','loc','lat','lon','U','dist','tt','tterr'])


    station_summary_df = pd.DataFrame(columns= ['eventid','version','event lat','event lon','time','mindist','maxdist','percent','near sta cnt','sta trig cnt',
               'active','inactive','nsta'])


    a_df = pd.DataFrame(columns=['eventid','version','event lat','event lon','depth','sta','chan','net','loc','lat','lon','date','tt',
               'dist','ttok','dok','nttok','ndok','ttmin','ttmax','relo','nlat','nlon','ndep','ntime','nrms','nave',
               'fitok','dmin','mdok','percent','pok'])
   
    station_counts_df = pd.DataFrame(columns=['eventid','version','sta','net','lat','lon','cluster','dist','tt','time','time check','active',
             'trig','clu','ctrig'])

    epic_location_df = pd.DataFrame(columns =['eventid','version','s','lat0','lon0','dep0','time0','lat','lon','depth','time','ddist','avefit','rmsfit','nT','nS'])


    for k in range(len(Lines)):
        l0 = list(filter(None, Lines[k].split('|')))[-1]
        l = list(filter(None, l0.split(' ')))
        #------------------------------------------------------------------------------------------------------------------------------------#
        if ('E:I:' in l) or ('E:I:F:' in l):
            # THIS IS THE EVENT INFORMATION
            if l[1]==epic_id:
                #print(l)
                ee=1
                
                
                event_summary_df.loc[len(event_summary_df)] =[int(l[1]),int(l[2]),float(l[3]),float(l[4]),
                                                            float(l[5]) ,float(l[6]),l[7],float(l[8]),float(l[9]),float(l[10]),
                                                            float(l[11]),float(l[12]),float(l[13]),int(l[14]),int(l[15]),int(l[16]),int(l[17]),
                                                            float(l[18]) ,float(l[19]),int(l[20]),int(l[21]),int(l[22]),
                                                            int(l[23]),int(l[24]),int(l[25]),int(l[26]),float(l[27]),
                                                            int(l[28]),float(l[29]),float(l[30]),int(l[31]),float(l[32]),int(l[33]),
                                                            int(l[34]),int(l[35]),float(l[36]),float(l[37]),float(l[38]),float(l[39]),
                                                            int(l[40]),int(l[41]),int(l[42]),l[43],l[44].split('\n')[0]]
        if ('E:I:T:' in l):
            #--------------------- EVENT INFO TRIGGER ------------------------#
            if l[1]==epic_id:
                
                event_triggers_df.loc[len(event_triggers_df)]=[l[1],int(l[2]),int(l[3]),int(l[4]),l[5],l[6],l[7],l[8],float(l[9]),float(l[10]),
                                                            l[11],int(l[12]),int(l[13]),float(l[14]),float(l[15]),int(l[16]),float(l[17]),float(l[18]),
                                                            l[19] ,float(l[20]),int(l[21]),float(l[22]),int(l[23]),int(l[24]),int(l[25]),
                                                            int(l[26]),int(l[27]),int(l[28]),int(l[29]),int(l[30]),float(l[31]),float(l[32]),float(l[33]),int(l[34]),
                                                            float(l[35]),float(l[36]),float(l[37]),int(l[38]),float(l[39])]
        if 'L:T:' in l:
            # THIS IS A LOCATION TRIGGER LINE
            if l[1]==epic_id:
                location_triggers_df.loc[len(location_triggers_df)]=[l[1],int(l[2]),int(l[3]),int(l[4]),l[5],l[6],l[7],l[8],float(l[9]),
                                                                   float(l[10]),int(l[11]),float(l[12]),float(l[13]),float(l[14])]
        if 'E:S:' in l:
            #STATION COUNT SUMMARY
            if l[1]==epic_id:
                station_summary_df.loc[len(station_summary_df)]=[l[1],int(l[2]),float(l[3]),float(l[4]),l[5],float(l[6]),
                                                                float(l[7]),float(l[8]),int(l[9]),int(l[10]),
                                                                int(l[11]),int(l[12]),int(l[13])]
        if "A:" in l:
            if l[1]==epic_id:
                a_df .loc[len(a_df)]=[l[1],int(l[2]),float(l[3]),float(l[4]),float(l[5]),l[6],l[7],l[8],l[9],
                                    float(l[10]),float(l[11]),l[12],float(l[13]),float(l[14]),int(l[15]),
                                    int(l[16]),int(l[17]),int(l[18]),float(l[19]),float(l[20]),
                                    int(l[21]),float(l[22]),float(l[23]),float(l[24]),l[25],float(l[26]),float(l[27]),
                                    int(l[28]),float(l[29]),int(l[30]),float(l[31]),int(l[32])]
        if 'E:C:' in l:
            #STATION COUNT DETAILS 
            if l[1]==epic_id:
                station_counts_df.loc[len(station_counts_df)]=[l[1],int(l[2]),l[3],l[4],float(l[5]),float(l[6]),
                                                             l[7],float(l[8]),float(l[9]),l[10],
                                                             float(l[11]),int(l[12]),int(l[13]),int(l[14]),int(l[15])]
        if 'L:E:' in l:
            # THIS IS A LOCATION ALGORITHM LINE
            if l[1]==epic_id:
                epic_location_df.loc[len(epic_location_df)]=[l[1],int(l[2]),int(l[3]),float(l[4]),float(l[5]),float(l[6]),l[7],float(l[8]),float(l[9]),
                                                            float(l[10]),l[11],float(l[12]),float(l[13]),float(l[14]),int(l[15]),int(l[16])]
                
                
        #------------------------------------------------------------------------------------------------------------------------------------#    



    # need to save the dataframes
    if ee ==1:
        
        
        if os.path.exists(project_parent_directory+event_id+'/')==False:
            os.mkdir(project_parent_directory+event_id+'/')
            os.mkdir(project_parent_directory+event_id+'/EPIC/')
            
            
        event_summary_df.to_csv(    project_parent_directory+event_id+'/EPIC/'+event_id+'_event_summary_log.txt',sep='\t',index=False)
        event_triggers_df.to_csv(   project_parent_directory+event_id+'/EPIC/'+event_id+'_event_triggers_log.txt',sep='\t',index=False)
        location_triggers_df.to_csv(project_parent_directory+event_id+'/EPIC/'+event_id+'_location_triggers_log.txt',sep='\t',index=False)
        station_summary_df.to_csv(  project_parent_directory+event_id+'/EPIC/'+event_id+'_station_summary_log.txt',sep='\t',index=False)
        a_df.to_csv(                project_parent_directory+event_id+'/EPIC/'+event_id+'_misc_log.txt',sep='\t',index=False)
        station_counts_df.to_csv(   project_parent_directory+event_id+'/EPIC/'+event_id+'_station_counts_log.txt',sep='\t',index=False)
        epic_location_df.to_csv(    project_parent_directory+event_id+'/EPIC/'+event_id+'_epic_location_log.txt',sep='\t',index=False)
        
        



def search_comcat(starttime,endtime,maxradiuskm,latitude,longitude):
    import requests
    from urllib.parse import urlencode
    ##############################################################################
    TIMEFMT = "%Y-%m-%dT%H:%M:%S"
    SEARCH_TEMPLATE = "https://earthquake.usgs.gov/fdsnws/event/1/query?format=geojson"
    libversion ='2.0.13'
    HEADERS = {"User-Agent": "libcomcat v%s" % libversion}
    TIMEOUT = 60
    
    
    template = SEARCH_TEMPLATE
    newargs = {}
    newargs['starttime']     = starttime.strftime(TIMEFMT)
    newargs['endtime']       = endtime.strftime(TIMEFMT)
    newargs['maxradiuskm']   = maxradiuskm
    newargs['longitude']     = longitude
    newargs['latitude']      = latitude
    newargs["limit"] = 20000
    
    paramstr = urlencode(newargs)
    url = template + "&" + paramstr
    
    response = requests.get(url, timeout=TIMEOUT, headers=HEADERS)
    jdict = response.json()
    ##############################################################################
    return(jdict)


def search_comcat_by_eventid(project_parent_directory,postgres_id,eventid):

    from datetime import datetime,timedelta
    import requests
    from urllib.parse import urlencode
    import pandas as pd
    ##############################################################################
    
    SEARCH_TEMPLATE = "https://earthquake.usgs.gov/fdsnws/event/1/query?format=geojson"
    libversion ='2.0.13'
    HEADERS = {"User-Agent": "libcomcat v%s" % libversion}
    TIMEOUT = 60
    
    
    template = SEARCH_TEMPLATE
    newargs = {}
    newargs['eventid']=eventid
    paramstr = urlencode(newargs)
    url = template + "&" + paramstr
    
    response = requests.get(url, timeout=TIMEOUT, headers=HEADERS)
    jdict = response.json()
    
    
    USGS_event_id=eventid
    USGS_event_latitude=jdict['geometry']['coordinates'][1]
    USGS_event_longitude=jdict['geometry']['coordinates'][0]
    USGS_event_depth = jdict['geometry']['coordinates'][2]
    USGS_event_magnitude =jdict['properties']['mag']
    time_in_msec = jdict['properties']['time']
    time_in_sec = time_in_msec // 1000
    msec = time_in_msec - (time_in_sec * 1000)
    dtime = datetime(1970, 1, 1) + timedelta(seconds=time_in_sec)
    dt = timedelta(milliseconds=msec)
    dtime = dtime + dt


    USGS_event_time=dtime.timestamp()
    
    
    catalog_df = pd.DataFrame({'postgres id':postgres_id,'USGS ID':USGS_event_id,'USGS time':USGS_event_time,
                                    'USGS lat':USGS_event_latitude,'USGS lon':USGS_event_longitude,
                                    'USGS depth':USGS_event_depth,'USGS mag':USGS_event_magnitude},index=[0])   
    
    catalog_df.to_csv(project_parent_directory+postgres_id+'/USGS/'+'usgs_event_summary.txt',sep='\t',index=False)
    
    
    return(catalog_df)

