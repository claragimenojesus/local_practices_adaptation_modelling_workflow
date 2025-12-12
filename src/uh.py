#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import xarray as xr
import rioxarray as rio
import geopandas as gpd
from shapely.geometry import mapping
import pandas as pd
import os
from datetime import datetime, timedelta
import re

INPUT_DIRECTORY = str(os.environ['MOSART_INPUT'])
INPUT_FILENAME = str(os.environ['INPUT_FILENAME'])
OUTPUT_DIRECTORY = str(os.environ['MOSART_OUTPUT'])
NBS = str(os.environ['NBS'])
SCENARIO_SH = str(os.environ['SCENARIO_SH'])
SC= str(os.environ['SC'])

#SC = re.findall('5/(.*)',SCENARIO_SH)[0]
## unit hydrograph

def unit_hydro(V,k,n,t,er):
    from scipy.special import gamma
    import numpy as np
    from scipy.sparse import spdiags

    #inputs
    # V = volume of unit hydrograph
    # k,n = parameters from Nash (1957)
    # t = unit hydrograph time
    # er = effective rainfall time-series

    # outputs
    # DGW = deeper groundwater flow
    # uh = unit hydrograph

    uh = np.array(V/(k*gamma(n))*np.multiply(np.exp(-t/k),np.power((t/k),n-1)))
    uh = uh/np.sum(uh)

    # from Ana Mijic as from Anthony script
    m = er.size
    r = t.size
    n = r+m-1

    pmat = spdiags(np.multiply(np.ones(r),np.row_stack(er)),np.arange(0,-er.size,-1),n,r)

    DGW = pmat*uh
    DGW = DGW[0:er.size]

    return [DGW, uh]

# k0 = 18.59720523 for KM105, and k0= 22.32687392 for PISAC
# n0 = 22.57835457 for KM105, and n0= 18.8277447 for PISAC
k0_km105=18.59720523
n0_km105=22.57835457
k0_pisac=22.32687392
n0_pisac=18.8277447
V=1

## load jules results, change accordingly
#jules = xr.open_dataset(home+"/hydrodata/jules_runs/coupled/runoff_scaled_and_partitioned.nc",decode_coords="all")
jules = xr.open_dataset(os.path.join(INPUT_DIRECTORY,NBS, INPUT_FILENAME),decode_coords="all") # MODIFY THIS DEPENDING ON HOW WE LOOP FOR BARIABLES

## Change depending shapefile
catchment_km105= gpd.read_file("/rds/general/user/cg2117/home/catchments/km105.shp")
catchment_pisac = gpd.read_file("/rds/general/user/cg2117/home/catchments/PISAC.shp")

## JULES gridded results as timeseries
q_deep=jules.q_deep.rio.set_spatial_dims(x_dim="lon",y_dim="lat",inplace=True)
q_deep.rio.set_spatial_dims(x_dim="lon",y_dim="lat",inplace=True)
q_deep.rio.write_crs("epsg:4326",inplace=True)

grid_area = xr.open_dataset(os.path.join(INPUT_DIRECTORY,NBS,"mosart_new.nc"))['area']
grid_area = grid_area.expand_dims("time")
grid_area['lat']=q_deep['lat']
grid_area['lon']=q_deep['lon']
deep_total_ts = np.multiply(q_deep, grid_area[0,:,:]).sum(dim=("lat","lon"))/1000
#deep_total_ts = q_deep.sum(dim=["lat","lon"])
#deep_total_ts = deep_total_ts*83424939/1000 # think we should multiply by each gridcell area in the grid instead...

q_deep_km105=q_deep.rio.clip(catchment_km105.geometry.apply(mapping),catchment_km105.crs,drop=False,all_touched=True)
deep_km105_ts = np.multiply(q_deep_km105, grid_area[0,:,:]).sum(dim=("lat","lon"))/1000
#deep_km105_ts = q_deep_km105.sum(dim=["lat","lon"])
#deep_km105_ts=deep_km105_ts*83424939/1000 # multiplying by mean pixel area for 48x48 grid

q_deep_pisac=q_deep.rio.clip(catchment_pisac.geometry.apply(mapping),catchment_pisac.crs,drop=False,all_touched=True)
deep_pisac_ts = np.multiply(q_deep_pisac, grid_area[0,:,:]).sum(dim=("lat","lon"))/1000
#deep_pisac_ts = q_deep_pisac.sum(dim=["lat","lon"])
#deep_pisac_ts=deep_pisac_ts*83424939/1000 # multiplying by mean pixel area for 48x48 grid



## convolve timeseries
[deeper_mod_total,uh] = unit_hydro(V,k0_km105,n0_km105,np.arange(deep_total_ts.size),deep_total_ts) # Assumption here is we use the uh parameters from km105 derived by Jose to be applied for the entire catchment so that we can obtain a timeseries of deep contributions for the full catchment

[deeper_mod_km105,uh]= unit_hydro(V,k0_km105,n0_km105,np.arange(deep_km105_ts.size),deep_km105_ts)

[deeper_mod_pisac,uh]= unit_hydro(V,k0_pisac,n0_pisac,np.arange(deep_pisac_ts.size),deep_pisac_ts)

## Write deeper baseflows to csv files
t = np.arange(datetime(2000,1,1), datetime(2100,1,1), timedelta(days=1)).astype(datetime)

df_total = pd.DataFrame({"time" : t, "deeper_baseflow" : deeper_mod_total})
df_total.to_csv(os.path.join(OUTPUT_DIRECTORY, NBS, SC+"_"+NBS+"deeper_total.csv"), index=False)

df_km105 = pd.DataFrame({"time" : t, "deeper_baseflow" : deeper_mod_km105})
df_km105.to_csv(os.path.join(OUTPUT_DIRECTORY,NBS,SC+"_"+NBS,"deeper_km105.csv"), index=False) # MODIFY PATHS HERE AS NEEDED

df_pisac = pd.DataFrame({"time" : t, "deeper_baseflow" : deeper_mod_pisac})
df_pisac.to_csv(os.path.join(OUTPUT_DIRECTORY,NBS,SC+"_"+NBS,"deeper_pisac.csv"), index=False) # MODIFY PATHS HERE AS NEEDED
