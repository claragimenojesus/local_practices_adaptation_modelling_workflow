#!/usr/bin/env python
# -*- coding: utf-8 -*-

import xarray as xr
import os
import numpy as np
import matplotlib.pyplot as plt

SCENARIO_OUTPUT = str(os.environ['SCENARIO_OUTPUT_SH'])
SCENARIO = str(os.environ['SCENARIO']) ## added by CGJ 16/07/2024
RC = str(os.environ['RC']) # added by CGJ 16/07/2024
HIST_OUTPUT = str(os.environ['HIST_OUTPUT_SH'])

#if SCENARIO != "GFDL-CM3":
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

alpha_q = -35.35/100 #10.77/100
alpha_qss =  3.39/100 #2.22/100

alpha_q_min = -57.92/100 #-0.60/100 #-43.85/100
alpha_q_max = -12.77/100 #8.16/100 #15/100

alpha_qss_min = -15.62/100 #-9.77/100 #52.73/100
alpha_qss_max = 22.40/100 #30.69 #0

conservation_mask =  xr.open_dataset("/rds/general/user/cg2117/home/cums_vub_2005.nc")['Band1']
land_frac = xr.open_dataset("/rds/general/user/cg2117/home/netcdf/jules_land_frac_ESA_rahu_clipped.nc")['land_frac']
conservation_mask['lat'] = land_frac['lat']
conservation_mask['lon'] = land_frac['lon']
conservation_mask = np.multiply(conservation_mask, land_frac)
conservation_mask = conservation_mask.fillna(0)

#hydro_output = xr.open_mfdataset("/rds/general/user/cg2117/home/project/rahu/jules-output/u-cz655/JULES_vn6.1.S2.daily_hydrology.20*.2D.nc")
#meteo_output = xr.open_mfdataset("/rds/general/user/cg2117/home/project/rahu/jules-output/u-cz655/JULES_vn6.1.S2.daily_meteo.20*.2D.nc")
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


#x_c = xr.open_mfdataset("/rds/general/user/cg2117/home/project/rahu/jules-output/u-cz655/JULES_vn6.1.S2.daily_hydrology.20*.2D.nc")
x_c = xr.open_dataset(os.path.join(SCENARIO_OUTPUT, "daily_hydrology_00_99.nc"))

# First, we modify the flows only taking into account the conserved area of each pixel (this will help for the operations)
x_c['runoff'] = np.maximum(conserved['runoff']*(1+alpha_q),0) # limit to positive flows

# NOT NEEDED FOR DEGRADATION RUN AS ALL ALPHA_Q<0, so shouldn't have to limit Q<=P
# To ensure balance over the entire timeseries, we adjust the total pixel runoff (which currently only takes into account conserved area) by a factor sum(Pconserved)/sum(Q).
#for i in range(len(LAT)):
#    for j in range(len(LON)):
#        if x_c['runoff'][...,i,j].sum(dim="time") > conserved_meteo['precip'][...,i,j].sum(dim="time"):
#            x_c['runoff'][...,i,j] *= np.divide(conserved_meteo['precip'][...,i,j].sum(dim="time"), x_c['runoff'][...,i,j].sum(dim="time"))

# Now, to modify the sub-surface runoff we need to ensure Qss<=Q in the conservation area
x_c['sub_surf_roff'] = np.minimum( np.maximum( conserved['sub_surf_roff']*(1+alpha_qss)*(1+alpha_q),0 ), x_c['runoff'])

# we modify only one of the ET fluxes, as anyways in the analysis we will add them all up so adding it to just one of the components is ok as it is a sum.
# We modify this flux before adding the non-conserved portions of the catchment. Here, ET-=increase in runoff
x_c['esoil_gb'] = x_c['esoil_gb'] + (conserved['runoff']-x_c['runoff'])- np.minimum(conserved['runoff']*(1+alpha_q),0) # adding negative flows here.

x_c['runoff'] += non_conserved['runoff']
x_c['sub_surf_roff'] += non_conserved['sub_surf_roff']
x_c['surf_roff'] = x_c['runoff']-x_c['sub_surf_roff']

#x_c.to_netcdf("/rds/general/user/cg2117/home/project/rahu/jules-output/degradation_maskv2_00_18.nc")
x_c.to_netcdf(os.path.join(SCENARIO_OUTPUT, "degradation_"+RC+"_"+SCENARIO+"_"+"00_99.nc"))
x_c.close()


#x_c = xr.open_mfdataset("/rds/general/user/cg2117/home/project/rahu/jules-output/u-cz655/JULES_vn6.1.S2.daily_hydrology.20*.2D.nc")
x_c = xr.open_dataset(os.path.join(SCENARIO_OUTPUT, "daily_hydrology_00_99.nc"))

# First, we modify the flows only taking into account the conserved area of each pixel (this will help for the operations)
x_c['runoff'] = np.maximum(conserved['runoff']*(1+alpha_q_min),0)
# To ensure balance over the entire timeseries, we adjust the total pixel runoff (which currently only takes into account conserved area) by a factor sum(Pconserved)/sum(Q).
#for i in range(len(LAT)):
#    for j in range(len(LON)):
#        if x_c['runoff'][...,i,j].sum(dim="time") > conserved_meteo['precip'][...,i,j].sum(dim="time"):
#            x_c['runoff'][...,i,j] *= np.divide(conserved_meteo['precip'][...,i,j].sum(dim="time"), x_c['runoff'][...,i,j].sum(dim="time"))

# Now, to modify the sub-surface runoff we need to ensure Qss<=Q in the conservation area
x_c['sub_surf_roff'] = np.minimum( np.maximum( conserved['sub_surf_roff']*(1+alpha_qss_min)*(1+alpha_q_min),0 ), x_c['runoff'])

# we modify only one of the ET fluxes, as anyways in the analysis we will add them all up so adding it to just one of the components is ok as it is a sum.
# We modify this flux before adding the non-conserved portions of the catchment. Here, ET-=increase in runoff
x_c['esoil_gb'] = x_c['esoil_gb'] + (conserved['runoff']-x_c['runoff']) - np.minimum(conserved['runoff']*(1+alpha_q_min),0) # this needs redoing

x_c['runoff'] += non_conserved['runoff']
x_c['sub_surf_roff'] += non_conserved['sub_surf_roff']
x_c['surf_roff'] = x_c['runoff']-x_c['sub_surf_roff']

#x_c.to_netcdf("/rds/general/user/cg2117/home/project/rahu/jules-output/degradation_min_maskv2_00_18.nc")
x_c.to_netcdf(os.path.join(SCENARIO_OUTPUT,"degradation_min_"+RC+"_"+SCENARIO+"_"+"00_99.nc"))
x_c.close()




#x_c = xr.open_mfdataset("/rds/general/user/cg2117/home/project/rahu/jules-output/u-cz655/JULES_vn6.1.S2.daily_hydrology.20*.2D.nc")
x_c = xr.open_dataset(os.path.join(SCENARIO_OUTPUT, "daily_hydrology_00_99.nc"))

# First, we modify the flows only taking into account the conserved area of each pixel (this will help for the operations)
x_c['runoff'] = np.maximum(conserved['runoff']*(1+alpha_q_max),0)

# To ensure balance over the entire timeseries, we adjust the total pixel runoff (which currently only takes into account conserved area) by a factor sum(Pconserved)/sum(Q).
#for i in range(len(LAT)):
#    for j in range(len(LON)):
#        if x_c['runoff'][...,i,j].sum(dim="time") > conserved_meteo['precip'][...,i,j].sum(dim="time"):
#            x_c['runoff'][...,i,j] *= np.divide(conserved_meteo['precip'][...,i,j].sum(dim="time"), x_c['runoff'][...,i,j].sum(dim="time"))

# Now, to modify the sub-surface runoff we need to ensure Qss<=Q in the conservation area
x_c['sub_surf_roff'] = np.minimum( np.maximum( conserved['sub_surf_roff']*(1+alpha_qss_max)*(1+alpha_q_max),0 ), x_c['runoff'])

# we modify only one of the ET fluxes, as anyways in the analysis we will add them all up so adding it to just one of the components is ok as it is a sum.
# We modify this flux before adding the non-conserved portions of the catchment. Here, ET-=increase in runoff
x_c['esoil_gb'] = x_c['esoil_gb'] + (conserved['runoff']-x_c['runoff']) - np.minimum(conserved['runoff']*(1+alpha_q_max),0) # this needs redoing

x_c['runoff'] += non_conserved['runoff']
x_c['sub_surf_roff'] += non_conserved['sub_surf_roff']
x_c['surf_roff'] = x_c['runoff']-x_c['sub_surf_roff']

#x_c.to_netcdf("/rds/general/user/cg2117/home/project/rahu/jules-output/degradation_max_maskv2_00_18.nc")
x_c.to_netcdf(os.path.join(SCENARIO_OUTPUT,"degradation_max_"+RC+"_"+SCENARIO+"_"+"00_99.nc"))
x_c.close()
