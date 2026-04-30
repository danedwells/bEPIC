#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 11:45:35 2022

@author: amy
"""
from bEPIC import bEPIC_main
postgres_id = 126625
#--------------------------------------#

init=0
run=1
#--------------------------------------#
# run variables
project_parent_directory = '/Users/amy/projects/bEPIC_event_files/'


velocity_model = 'h2p+ak135'   # 'constant
GridSize=200
GridSpacing=2


if init ==1:
    bEPIC_main.initialize_bEPIC_event(project_parent_directory,postgres_id)
    
if run ==1:
    bEPIC_main.run_bEPIC(project_parent_directory,postgres_id,velocity_model,GridSize,GridSpacing)

