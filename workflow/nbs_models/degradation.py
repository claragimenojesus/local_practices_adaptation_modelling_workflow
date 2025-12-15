#!/usr/bin/env python
# -*- coding: utf-8 -*-

import xarray as xr
import os
import numpy as np
import matplotlib.pyplot as plt

SCENARIO_OUTPUT = str(os.environ['SCENARIO_OUTPUT_SH'])
SCENARIO = str(os.environ['SCENARIO']) 
RC = str(os.environ['RC']) 
HIST_OUTPUT = str(os.environ['HIST_OUTPUT_SH'])
DEG_MASK = str(os.environ['DEG_MASK'])
LAND_FRAC_MASK = str(os.environ['LAND_FRAC_MASK'])

vars = ['runoff','surf_roff','sub_surf_roff','ei_gb','ecan_gb','esoil_gb','elake','fao_et0']
hydro_hist = xr.open_mfdataset(HIST_OUTPUT+"/*hydrology*20*2D.nc")[vars]
hydro_fut = xr.open_mfdataset(SCENARIO_OUTPUT+"/*hydrology*20*2D.nc")[vars]
meteo_hist = xr.open_mfdataset(HIST_OUTPUT+"/*meteo*20*2D.nc")[['precip','rainfall']]
meteo_fut = xr.open_mfdataset(SCENARIO_OUTPUT+"/*meteo*20*2D.nc")[['precip','rainfall']]
meteo_00_99 = xr.concat([meteo_hist, meteo_fut], dim="time")
meteo_00_99.to_netcdf(SCENARIO_OUTPUT+"/daily_meteo_00_99.nc")
hydro_00_99 = xr.concat([hydro_hist, hydro_fut], dim="time")
hydro_00_99.to_netcdf(SCENARIO_OUTPUT+"/daily_hydrology_00_99.nc")
hydro_hist.close()
hydro_fut.close()
hydro_00_99.close()

# Loading parameters derived using one-sided t-test
alpha_q = -35.35/100 
alpha_qss =  3.39/100

alpha_q_min = -57.92/100 
alpha_q_max = -12.77/100

alpha_qss_min = -15.62/100
alpha_qss_max = 22.40/100 

conservation_mask =  xr.open_dataset(DEG_MASK)['Band1']
land_frac = xr.open_dataset(LAND_FRAC_MASK)['land_frac']
conservation_mask['lat'] = land_frac['lat']
conservation_mask['lon'] = land_frac['lon']
conservation_mask = np.multiply(conservation_mask, land_frac)
conservation_mask = conservation_mask.fillna(0)

hydro_output = xr.open_dataset(os.path.join(SCENARIO_OUTPUT, "daily_hydrology_00_99.nc"))
meteo_output = xr.open_dataset(os.path.join(SCENARIO_OUTPUT, "daily_meteo_00_99.nc"))

conservation_mask = conservation_mask.expand_dims("time")
LAT = conservation_mask['lat']
LON = conservation_mask['lon']
NLAT = len(LAT)
NLON = len(LON)

vars = ['runoff','surf_roff','sub_surf_roff','ei_gb','ecan_gb','esoil_gb','elake']
conserved = np.multiply(hydro_output[vars],conservation_mask[0,:,:])#.where(conservation_mask[0,:,:]>0)
non_conserved = hydro_output[vars]-conserved[vars]

conserved_meteo = np.multiply(meteo_output, conservation_mask[0,:,:])#.where(conservation_mask[0,:,:]>0)
non_conserved_meteo = meteo_output-conserved_meteo

hydro_output.close()
meteo_output.close()

x_c = xr.open_dataset(os.path.join(SCENARIO_OUTPUT, "daily_hydrology_00_99.nc"))

# First, we modify the flows only taking into account the conserved area of each pixel (this will help for the operations)
x_c['runoff'] = np.maximum(conserved['runoff']*(1+alpha_q),0) # limit to positive flows

# Now, to modify the sub-surface runoff we need to ensure Qss<=Q in the conservation area
x_c['sub_surf_roff'] = np.minimum( np.maximum( conserved['sub_surf_roff']*(1+alpha_qss)*(1+alpha_q),0 ), x_c['runoff'])

# we modify only one of the ET fluxes, as anyways in the analysis we will add them all up so adding it to just one of the components is ok as it is a sum.
# We modify this flux before adding the non-conserved portions of the catchment. Here, ET-=increase in runoff
x_c['esoil_gb'] = x_c['esoil_gb'] + (conserved['runoff']-x_c['runoff'])#- np.minimum(conserved['runoff']*(1+alpha_q),0) # adding negative flows here.

x_c['runoff'] += non_conserved['runoff']
x_c['sub_surf_roff'] += non_conserved['sub_surf_roff']
x_c['surf_roff'] = x_c['runoff']-x_c['sub_surf_roff']

x_c.to_netcdf(os.path.join(SCENARIO_OUTPUT, "degradation_"+RC+"_"+SCENARIO+"_"+"00_99.nc"))
x_c.close()



# Now for the lower bound parameters

x_c = xr.open_dataset(os.path.join(SCENARIO_OUTPUT, "daily_hydrology_00_99.nc"))

# First, we modify the flows only taking into account the conserved area of each pixel (this will help for the operations)
x_c['runoff'] = np.maximum(conserved['runoff']*(1+alpha_q_min),0)

# Now, to modify the sub-surface runoff we need to ensure Qss<=Q in the conservation area
x_c['sub_surf_roff'] = np.minimum( np.maximum( conserved['sub_surf_roff']*(1+alpha_qss_min)*(1+alpha_q_min),0 ), x_c['runoff'])

# we modify only one of the ET fluxes, as anyways in the analysis we will add them all up so adding it to just one of the components is ok as it is a sum.
# We modify this flux before adding the non-conserved portions of the catchment. Here, ET-=increase in runoff
x_c['esoil_gb'] = x_c['esoil_gb'] + (conserved['runoff']-x_c['runoff']) #- np.minimum(conserved['runoff']*(1+alpha_q_min),0) # this needs redoing

x_c['runoff'] += non_conserved['runoff']
x_c['sub_surf_roff'] += non_conserved['sub_surf_roff']
x_c['surf_roff'] = x_c['runoff']-x_c['sub_surf_roff']

x_c.to_netcdf(os.path.join(SCENARIO_OUTPUT,"degradation_min_"+RC+"_"+SCENARIO+"_"+"00_99.nc"))
x_c.close()


# Now for the higher bound parameters

x_c = xr.open_dataset(os.path.join(SCENARIO_OUTPUT, "daily_hydrology_00_99.nc"))

# First, we modify the flows only taking into account the conserved area of each pixel (this will help for the operations)
x_c['runoff'] = np.maximum(conserved['runoff']*(1+alpha_q_max),0)

# Now, to modify the sub-surface runoff we need to ensure Qss<=Q in the conservation area
x_c['sub_surf_roff'] = np.minimum( np.maximum( conserved['sub_surf_roff']*(1+alpha_qss_max)*(1+alpha_q_max),0 ), x_c['runoff'])

# we modify only one of the ET fluxes, as anyways in the analysis we will add them all up so adding it to just one of the components is ok as it is a sum.
# We modify this flux before adding the non-conserved portions of the catchment. Here, ET-=increase in runoff
x_c['esoil_gb'] = x_c['esoil_gb'] + (conserved['runoff']-x_c['runoff']) #- np.minimum(conserved['runoff']*(1+alpha_q_max),0) # this needs redoing

x_c['runoff'] += non_conserved['runoff']
x_c['sub_surf_roff'] += non_conserved['sub_surf_roff']
x_c['surf_roff'] = x_c['runoff']-x_c['sub_surf_roff']

x_c.to_netcdf(os.path.join(SCENARIO_OUTPUT,"degradation_max_"+RC+"_"+SCENARIO+"_"+"00_99.nc"))
x_c.close()
