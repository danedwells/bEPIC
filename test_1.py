#%%

import numpy as np
import torch
torch.set_num_threads(1)
import matplotlib.pyplot as plt
from bEPIC import bEPIC_main

#%%

postgres_id = '9999'  # References some database (SQL) 


#postgres_id = 126625
#--------------------------------------#

init=1
run=1
#--------------------------------------#
# run variables
project_parent_directory = '.'


velocity_model = 'h2p+ak135'   # 'constant
GridSize=200
GridSpacing=2


if init ==1:
    bEPIC_main.initialize_bEPIC_event(project_parent_directory,postgres_id)

#%%
    
if run ==1:
    bEPIC_main.run_bEPIC(project_parent_directory,postgres_id,velocity_model,GridSize,GridSpacing)

