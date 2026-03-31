#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 16:49:31 2022

@author: amy
"""
import sys
import os
import shutil # higher level version of os. Caution: can delete entire directories without issue.
import gzip # python library for interacting with gzip files (GNU compression software)
import urllib.request # Python library for opening and reading URLs
from bEPIC import data_util

#-------------------------------------------------------------#
# VARIABLE 1 needs to be instance
# VARIABLE 2 needs to be day
# VARIABLE 2 needs to be postges_d
bepic = os.path.dirname(os.path.abspath(__file__))



instance = sys.argv[1]
day = sys.argv[2]
postgres_id = sys.argv[3]
epic_id = sys.argv[4]

#-------------------------------------------------------------#
instance = instance.replace(" ", "").split('@')[-1]

project_parent_directory='/home/gcl/RA/williamson/bEPIC_events/'
log_location='/home/gcl/RA/williamson/EPIC_unprocessed_logs/'
log_file = log_location+instance+'_'+day+'.log'
event_id = str(postgres_id)

if os.path.exists(project_parent_directory+event_id+'/'+event_id+'_event_summary_log.txt')==False:
    try:
        url='http://131.215.66.120/'+instance+'/epic/epic_'+day+'.log.gz'
        if os.path.exists(log_location+instance+'_'+day+'.log')==False:
            print(' ... downloading log file from EPIC archive')
            urllib.request.urlretrieve(url,log_location+instance+'_'+day+'.log.gz')
            with gzip.open(log_location+instance+'_'+day+'.log.gz', 'rb') as f_in:
                with open(log_location+instance+'_'+day+'.log', 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
        else:
            print(' ... log file '+instance+'_'+day+'.log alread exsists')
        try:
            print(' ... parsing log file ...')
            data_util.parse_log(project_parent_directory,log_file,event_id,epic_id)
            #os.remove(log_file)
            os.remove(log_file+'.gz')
        except:
            a=1
    except:
        print(' ... erorr: file not found')

else:
    print(' ... event already has a local log file. exiting....')

