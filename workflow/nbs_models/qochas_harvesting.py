#!/usr/bin/env python
# -*- coding: utf-8 -*-

import xarray as xr
import rasterio
from rasterio.plot import show
from matplotlib import pyplot
import os
import rioxarray as rxr
import geopandas as gpd
import numpy as np
import pandas as pd
import scipy.stats as stats
from shapely.geometry import mapping
import xesmf as xe

ephemeral = str(os.environ['EPHEMERAL'])
wd_hpc=str(os.environ['PBS_O_WORKDIR'])
home=str(os.environ['HOME'])
SCENARIO_SH = str(os.environ['SCENARIO_SH'])
HIST_OUTPUT = str(os.environ['HIST_OUTPUT_SH'])
SCENARIO_SH_OUTPUT = str(os.environ['SCENARIO_OUTPUT_SH'])
RC= str(os.environ['RC'])
SC = str(os.environ['SCENARIO'])

## jules output. change accordingly
#jules = xr.open_dataset(home+"/hydrodata/jules_runs/coupled_for_nbs/rcp85_ACCESS1-3_coupled_jules_oggm_00_99.nc",decode_coords="all").sel(time=slice("2000-01-01", "2018-12-31"))
jules = xr.open_dataset(os.path.join(SCENARIO_SH_OUTPUT, RC+"_"+SC+"_degradation_coupled_jules_oggm_00_99.nc"), decode_coords="all")

#mask >4000masl
mask=xr.open_dataset("/rds/general/user/cg2117/home/netcdf/mask_4000m.nc")['frac']

### Qochas parameter relations derived from medians from a Sierra Azul small dataset.
# We will produce linear spaces for storage volume but the range of storage volume would be delimited by the
# corresponding contribution area as this cannot exceed the cell size.
# Storage to contribution area is 13729.72/202500 = 0.06558. The cell size is an average of 1.56E7 m2.
# then the storage volume will be 1023030 M3, so 1E6 m3.
# [1] 13729.72 ## storage volume m3
# [1] 202500 ## contribution area m2
# [1] 11652.64 ## qocha area m2
st_max=1E6
vol_area=11652.64/13729.72
vol_acc=202500/13729.72


# cell area
cell_area = xr.open_dataset("/rds/general/user/cg2117/home/netcdf/gridcell_area_rahu_v2.nc")['area']

## masking cells where there exists flows over the entire time
masking=jules.surf_roff.sum(dim="time")
masking=masking.where(masking>0)



qochas_cap=st_max
qochas_cap=np.multiply(xr.full_like(mask,st_max).where(masking>0),mask)
qochas_area = qochas_cap*vol_area
qochas_acc = qochas_cap*vol_acc

### Zero arrays for various intermediate variables
St=xr.zeros_like(jules.surf_roff) #Storage at time step t
Qav = xr.zeros_like(cell_area[0,...]) # storage volume + contributing runoff - losses (Et + drainage). or available storage
#R = xr.zeros_like(jules.surf_roff) # Potential recharge at time step t
Qin = xr.zeros_like(cell_area[0,...])
ET = xr.zeros_like(jules.surf_roff)


### Qocha model
# cell_area=15526711 #grid_area[0,j,i].values for grid area
for t in range(jules.sizes["time"]-1):
    #### Calculate losses

    ## calculate ET
    ET[t,...]=np.minimum(np.abs(jules.fao_et0[t,...]) * qochas_area * 86.4,St[t,...])
    St[t,...]=np.add(St[t,...],-ET[t,...])
    ## Calculate available storage
    Qav=np.add(qochas_cap,-St[t,...])
    # water that actually enters the qocha is limited by the storage capacity

    # Calculate inflow. This is considering overflow, i.e. cannot flow more than available storage
    Qin= np.minimum(qochas_acc*jules.surf_roff[t,...]*86.4,Qav)
    ## correcting flow
    jules.surf_roff[t,...]= (jules.surf_roff[t,...]-(np.divide(Qin.values/86.4,cell_area[0,...])))
    #update storage
    St[t+1,...] = np.add(Qin,St[t,...])

jules.runoff.values = jules.surf_roff.values + jules.sub_surf_roff.values
jules['ET_qocha']=ET
jules['St_qocha']=St
jules = jules.drop_vars(["melt","perc","rain","area","precip","rainfall","snowfall","fao_et0","esoil_gb","ecan_gb","ei_gb","elake"])

# Save netcdf
jules.to_netcdf(os.path.join(SCENARIO_SH_OUTPUT, RC+"_"+SC+"_"+"qochas_harvesting_coupled_jules_oggm_00_99.nc"))


#available_storage[counter,:]=St.sel(time=jules.surf_roff.time.dt.year.isin(np.arange(2080,2101)),drop=True).resample(time="ME").mean().groupby("time.month").mean().mean(dim=["lat","lon"])
#et_month[counter,:]=ET.sel(time=jules.surf_roff.time.dt.year.isin(np.arange(2080,2101)),drop=True).resample(time="ME").mean().groupby("time.month").mean().mean(dim=["lat","lon"])
