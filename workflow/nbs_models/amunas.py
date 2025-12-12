#!/usr/bin/env python
# -*- coding: utf-8 -*-

#import libraries
import xarray as xr
import rasterio
from rasterio.plot import show
from matplotlib import pyplot
import os
import rioxarray as rxr
import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import mapping
import matplotlib.colors as colors
import xesmf as xe

ephemeral = str(os.environ['EPHEMERAL'])
wd_hpc=str(os.environ['PBS_O_WORKDIR'])
home=str(os.environ['HOME'])
SCENARIO_SH = str(os.environ['SCENARIO_SH'])
#JULES_PATH = str(os.environ['JULES_PATH_SH'])
HIST_OUTPUT = str(os.environ['HIST_OUTPUT_SH'])
SCENARIO_SH_OUTPUT = str(os.environ['SCENARIO_OUTPUT_SH'])
RC = str(os.environ['RC'])
SC = str(os.environ['SCENARIO'])

## import jules runs (note this is 2000 to 2099)
## scenario rcp45_ACCESS1-0
#jules = xr.open_dataset(home+"/hydrodata/jules_runs/coupled_for_nbs/rcp85_ACCESS1-3_coupled_jules_oggm_00_99.nc",decode_coords="all")
jules = xr.open_dataset(os.path.join(SCENARIO_SH_OUTPUT, RC+"_"+SC+"_degradation_coupled_jules_oggm_00_99.nc"), decode_coords="all")

#mask >4000masl
mask=xr.open_dataset("/rds/general/user/cg2117/home/netcdf/mask_4000m.nc")['frac']

## load UH
uh = pd.read_csv(home+'/unit_hydro.csv').fillna(0)
uh["conv"]=uh["conv"]/uh["conv"].sum()

area=xr.zeros_like(jules.surf_roff)
ET0=xr.zeros_like(jules.surf_roff)

div_ratio=100 # 100% diversion scenario

## et area will vary linearly with the amount of water diverted
et_ratio = 1.2E6/75 #75l/s equated to 1.2 km2 of greened area 

## calculating diverted water
diversion=np.multiply(jules.surf_roff,mask)
diversion= div_ratio/100*diversion.where(diversion.time.dt.month.isin([1,2,3,4,12])) # kg/m2/s
area= np.multiply(diversion,cell_area[0,...])*et_ratio # this is the greened area
## Correct surface flow
jules["surf_roff"]= jules.surf_roff-diversion.fillna(0)

# Correct subsurface flow
# Calculate ET
ET0=np.divide(np.multiply(area,jules.fao_et0), #kg/s
              cell_area[0,...]) #kg/m2/s
# Calculate infiltration
infilt=diversion.fillna(0)-ET0.fillna(0) ##infiltration is now in kg/m2/s
## unit hydrograph i,j wise for shallow subsurface flow
shallow_asnp = infilt.values
shallow_asnp = np.append(shallow_asnp,np.zeros((365,75,102)),axis=0)
shallow_asnp1 = np.zeros_like(shallow_asnp)
for i in range(jules.sizes["lon"]):
    for j in range(jules.sizes["lat"]):
        for t in range(jules.sizes["time"]):
            if shallow_asnp[t,j,i] > 0:
                shallow_asnp1[t:t+uh["conv"].size,j,i]+= shallow_asnp[t,j,i]*uh["conv"].values
                


#replace flow
jules.sub_surf_roff.values += shallow_asnp1[:jules.sizes["time"],...]
jules.runoff.values = jules.surf_roff.values + jules.sub_surf_roff.values
jules.esoil_gb.values += ET0.values
# Save netcdf
jules.to_netcdf(os.path.join(SCENARIO_SH_OUTPUT, RC+"_"+SC+"_"+"amunas_coupled_jules_oggm_00_99.nc"))
