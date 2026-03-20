#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 11:39:24 2026

@author: amy
"""

import EPIC_locate_prelim


#   set up params
params = EPIC_locate_prelim.EPIC_PARAMS()
params.PriorGridFile =  '/Users/amy/projects/container_bEPIC/data/prior_seis_grid_US_Canada.tt3'
params.use_prior = True
params.GridSize = 25
params.GridKm = 50
params.method = 'EPIC C'  




# set up starting event info


# this code assumes that I am inputting an event, so need to populate my event class
# with dummy values




# set up trigs

# initial location, the rest are just dummy variables
event = EPIC_locate_prelim.Event(lat = 36.764, 
              lon = -121.4472, 
              time = 1538771380.11, 
              misfit_rms = 0, 
              misfit_ave = 0, 
              eventid = 126625, 
              version = 0)


''' 126625.run
version,order,station,channel,network,location,longitude,latitude,trigger time,tterr,logPd
0,1,SAO,HNZ,BK,00,-121.4472,36.764,1538771380.09,-0.053,-1.89848
0,2,SAO,HHZ,BK,00,-121.4472,36.764,1538771380.11,-0.073,-1.875756
0,3,BSR,HNZ,NC,--,-121.5203,36.6674,1538771381.18,-0.035,-2.964038
0,5,PACP,HHZ,BK,00,-121.287,37.008,1538771383.34,0.081,-2.38656
0,4,PACP,HNZ,BK,00,-121.287,37.008,1538771383.34,0.081,-2.428661
'''
t = EPIC_locate_prelim.TriggerManager(lon = -121.4472, lat = 36.764, sta='SAO', net='BK', chan='HNZ',trigger_time = 1538771380.09)
event.trigs.append(t)

t = EPIC_locate_prelim.TriggerManager(lon = -121.4472, lat = 36.764, sta='SAO', net='BK', chan='HHZ',trigger_time = 1538771380.11)
event.trigs.append(t)

t = EPIC_locate_prelim.TriggerManager(lon = -121.5203, lat = 36.6674, sta='BSR', net='NC', chan='HNZ',trigger_time =1538771381.18)
event.trigs.append(t)

t = EPIC_locate_prelim.TriggerManager(lon =-121.287, lat = 37.008, sta='PACP', net='BK', chan='HHZ',trigger_time =1538771383.34)
event.trigs.append(t)

t = EPIC_locate_prelim.TriggerManager(lon =-121.287, lat = 37.008, sta='PACP', net='BK', chan='HNZ',trigger_time =1538771383.34)
event.trigs.append(t)

t,output_df = EPIC_locate_prelim.E2Location_locate(params,event)

