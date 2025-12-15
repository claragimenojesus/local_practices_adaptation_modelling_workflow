#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import netCDF4 as nc
import pandas as pd
import xarray as xr
import glob
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset, date2num, num2date
import datetime as dt
import pyreadr
import re

SCENARIO_SH = str(os.environ['SCENARIO_SH'])
#JULES_PATH = str(os.environ['JULES_PATH_SH'])
HIST_OUTPUT = str(os.environ['HIST_OUTPUT_SH'])
SCENARIO_SH_OUTPUT = str(os.environ['SCENARIO_OUTPUT_SH'])

melt_d = pyreadr.read_r(SCENARIO_SH+'/daily_melt_gridded.rda')
melt_d = melt_d['DFm']
melt_d = melt_d.rename(({'dim_0':'time','dim_1':'lat','dim_2':'lon'}))

rain_d = pyreadr.read_r(SCENARIO_SH+'/daily_rain_gridded.rda')
rain_d = rain_d['DFr']
rain_d = rain_d.rename(({'dim_0':'time','dim_1':'lat','dim_2':'lon'}))

MELT = melt_d
AREA = xr.open_dataset(".../gridcell_area.nc")['area']

nc = xr.open_dataset(os.path.join(SCENARIO_SH_OUTPUT, "JULES_vn6.1.S2.daily_hydrology.2019.2D.nc"))
LAT = nc['lat']
LAT = LAT.reindex(lat=list(reversed(LAT.lat)))
LON = nc['lon']
dates = pd.date_range(start=dt.datetime(2000,1,1,0),end=dt.datetime(2099,12,31,0),freq='D')
dates = dates.to_pydatetime()
ncfile = Dataset(SCENARIO_SH+"/melt_daily_kgm2s.nc", mode="w", format="NETCDF4")
lat_dim = ncfile.createDimension('lat', 75)
lon_dim = ncfile.createDimension('lon', 102)
time_dim = ncfile.createDimension('time', 36525)
lat = ncfile.createVariable('lat', np.float64, ('lat',))
lat.units = 'degrees_north'
lat.long_name = 'latitude'
lon = ncfile.createVariable('lon', np.float64, ('lon',))
lon.units = 'degrees_east'
lon.long_name = 'longitude'
time = ncfile.createVariable('time', np.float64, ('time',))
time.units = 'days since 2000-01-01'
time.long_name = 'time'
time.calendar = 'proleptic_gregorian'
TIME = date2num(dates, time.units)
melt = ncfile.createVariable('melt', np.float64, ('time','lat','lon'))
melt.units = 'kg m-2 s-1'
melt.standard_name = 'Daily gridbox glacier melt'
nlats = len(lat_dim)
nlons = len(lon_dim)
ntimes = len(time_dim)
lat[:] = LAT
lon[:] = LON
melt[:,:,:] = MELT*1000/AREA[0,:,:]
time[:] = TIME
ncfile.close()

nc_melt = xr.open_dataset(SCENARIO_SH+"/melt_daily_kgm2s.nc")
nc_melt.reindex(lat=list(reversed(nc_melt.lat))).to_netcdf(SCENARIO_SH+"/melt_daily_kgm2s_v2.nc")
nc_melt.close()

RAIN = rain_d
#nc = xr.open_dataset(JULES_PATH)
LAT = nc['lat']
LAT = LAT.reindex(lat=list(reversed(LAT.lat)))
LON = nc['lon']
dates = pd.date_range(start=dt.datetime(2000,1,1,0),end=dt.datetime(2099,12,31,0),freq='D')
dates = dates.to_pydatetime()
ncfile = Dataset(SCENARIO_SH+"/rain_daily_kgm2s.nc", mode="w", format="NETCDF4")
lat_dim = ncfile.createDimension('lat', 75)
lon_dim = ncfile.createDimension('lon', 102)
time_dim = ncfile.createDimension('time', 36525)
lat = ncfile.createVariable('lat', np.float64, ('lat',))
lat.units = 'degrees_north'
lat.long_name = 'latitude'
lon = ncfile.createVariable('lon', np.float64, ('lon',))
lon.units = 'degrees_east'
lon.long_name = 'longitude'
time = ncfile.createVariable('time', np.float64, ('time',))
time.units = 'days since 2000-01-01'
time.long_name = 'time'
time.calendar = 'proleptic_gregorian'
TIME = date2num(dates, time.units)
rain = ncfile.createVariable('rain', np.float64, ('time','lat','lon'))
rain.units = 'kg m-2 s-1'
rain.standard_name = 'Daily gridbox rain over glacier area directly contributing to runoff'
nlats = len(lat_dim)
nlons = len(lon_dim)
ntimes = len(time_dim)
lat[:] = LAT
lon[:] = LON
rain[:,:,:] = RAIN*1000/AREA[0,:,:]
time[:] = TIME
ncfile.close()

nc_rain = xr.open_dataset(SCENARIO_SH+"/rain_daily_kgm2s.nc")
nc_rain.reindex(lat=list(reversed(nc_rain.lat))).to_netcdf(SCENARIO_SH+"/rain_daily_kgm2s_v2.nc")
nc_rain.close()


# and create area and perc files for correct latitudes
nc_area = xr.open_dataset(SCENARIO_SH+"/glacier_area_daily_v2.nc")
nc_area.reindex(lat=list(reversed(nc_area.lat))).to_netcdf(SCENARIO_SH+"/glacier_area_daily_v3.nc")
nc_area.close()

nc_perc = xr.open_dataset(SCENARIO_SH+"/glacier_perc_daily_v2.nc")
nc_perc.reindex(lat=list(reversed(nc_perc.lat))).to_netcdf(SCENARIO_SH+"/glacier_perc_daily_v3.nc")
nc_perc.close()
AREA.close()

vars = ['runoff','surf_roff','sub_surf_roff','ei_gb','ecan_gb','esoil_gb','elake','fao_et0']
hydro_hist = xr.open_mfdataset(HIST_OUTPUT+"/*hydrology*20*2D.nc")[vars]
hydro_fut = xr.open_mfdataset(SCENARIO_SH_OUTPUT+"/*hydrology*20*2D.nc")[vars]
meteo_hist = xr.open_mfdataset(HIST_OUTPUT+"/*meteo*20*2D.nc")[['precip','rainfall']]
meteo_fut = xr.open_mfdataset(SCENARIO_SH_OUTPUT+"/*meteo*20*2D.nc")[['precip','rainfall']]
meteo_00_99 = xr.concat([meteo_hist, meteo_fut], dim="time")
meteo_00_99.to_netcdf(SCENARIO_SH_OUTPUT+"/daily_meteo_00_99.nc")
hydro_00_99 = xr.concat([hydro_hist, hydro_fut], dim="time")
hydro_00_99.to_netcdf(SCENARIO_SH_OUTPUT+"/daily_hydrology_00_99.nc")
hydro_hist.close()
hydro_fut.close()
hydro_00_99.close()


perc = xr.open_dataset(SCENARIO_SH+"/glacier_perc_daily_v3.nc")
daily_melt = xr.open_dataset(SCENARIO_SH+"/melt_daily_kgm2s_v2.nc")
daily_rain = xr.open_dataset(SCENARIO_SH+"/rain_daily_kgm2s_v2.nc")
daily_jules = xr.open_dataset(SCENARIO_SH_OUTPUT+"/daily_hydrology_00_99.nc")
daily_meteo = xr.open_dataset(SCENARIO_SH_OUTPUT+"/daily_meteo_00_99.nc")
daily_jules['time'] = daily_rain['time']
perc['time'] = daily_rain['time']
daily_meteo['time'] = daily_rain['time']
daily_area = xr.open_dataset(SCENARIO_SH+"/glacier_area_daily_v3.nc")


# Now let's create our new modified netcdf
LAT = daily_jules['lat']
LON = daily_jules['lon']
MELT = daily_melt['melt']
MELT = MELT.fillna(0)
PERC = perc['percentage']
PERC = PERC.fillna(0)
RUNOFF = daily_jules['runoff']
RUNOFF = RUNOFF.fillna(0)
SURF_ROFF = daily_jules['surf_roff']
SURF_ROFF = SURF_ROFF.fillna(0)
SUB_SURF_ROFF = daily_jules['sub_surf_roff']
SUB_SURF_ROFF = SUB_SURF_ROFF.fillna(0)
ET0 = daily_jules['fao_et0']
ET0 = ET0.fillna(0)
EI_GB = daily_jules['ei_gb']
EI_GB = EI_GB.fillna(0)
ELAKE = daily_jules['elake']
ELAKE = ELAKE.fillna(0)
ECAN_GB = daily_jules['ecan_gb']
ECAN_GB = ECAN_GB.fillna(0)
ESOIL_GB = daily_jules['esoil_gb']
ESOIL_GB = ESOIL_GB.fillna(0)
RAIN = daily_rain['rain']
RAIN = RAIN.fillna(0)
AREA = daily_area['area']
AREA = AREA.fillna(0)
PRECIP = daily_meteo['precip']
PRECIP = PRECIP.fillna(0)
RAINFALL = daily_meteo['rainfall']
RAINFALL = RAINFALL.fillna(0)

dates = pd.date_range(start=dt.datetime(2000,1,1,0),end=dt.datetime(2099,12,31,0),freq='D')
dates = dates.to_pydatetime()

RC = re.findall('(?<=\/JULES_OGGM\/)([^\/]+)',SCENARIO_SH)[0]
SC = re.findall('5/(.*)',SCENARIO_SH)[0]

ncfile = Dataset(os.path.join(SCENARIO_SH_OUTPUT+RC+"/"+SC+"/"+RC+"_"+SC+"_"+"coupled_jules_oggm_00_99.nc"), mode="w", format="NETCDF4")     

lat_dim = ncfile.createDimension('lat', 75)
lon_dim = ncfile.createDimension('lon', 102)
time_dim = ncfile.createDimension('time', 36525)
lat = ncfile.createVariable('lat', np.float64, ('lat',))
lat.units = 'degrees_north'
lat.long_name = 'latitude'
lon = ncfile.createVariable('lon', np.float64, ('lon',))
lon.units = 'degrees_east'
lon.long_name = 'longitude'
time = ncfile.createVariable('time', np.float64, ('time',))
time.units = 'days since 2019-01-01'
time.long_name = 'time'
time.calendar = 'proleptic_gregorian'

TIME = date2num(dates, time.units)

melt = ncfile.createVariable('melt', np.float64, ('time','lat','lon'))
melt.units = 'kg m-2 s-1'
melt.standard_name = 'Glacier melt per glacier area rate'

runoff = ncfile.createVariable('runoff', np.float64, ('time','lat','lon'))
runoff.units = 'kg m-2 s-1'
runoff.standard_name = 'Gridbox runoff rate'

surf_roff = ncfile.createVariable('surf_roff', np.float64, ('time','lat','lon'))
surf_roff.units = 'kg m-2 s-1'
surf_roff.standard_name = 'Gridbox surface runoff'

sub_surf_roff = ncfile.createVariable('sub_surf_roff', np.float64, ('time','lat','lon'))
sub_surf_roff.units = 'kg m-2 s-1'
sub_surf_roff.standard_name = 'Gridbox sub-surface runoff'

et0 = ncfile.createVariable('fao_et0', np.float64, ('time','lat','lon'))
et0.units = 'kg m-2 s-1'
et0.standard_name = 'FAO Penman-Monteith evapotranspiration for reference crop'

ei_gb = ncfile.createVariable('ei_gb', np.float64, ('time','lat','lon'))
ei_gb.units = 'kg m-2 s-1'
ei_gb.standard_name = 'Gridbox sublimation from lying snow or sea-ice'

elake = ncfile.createVariable('elake', np.float64, ('time','lat','lon'))
elake.units = 'kg m-2 s-1'
elake.standard_name = 'Gridbox mean evaporation from lakes'

esoil_gb = ncfile.createVariable('esoil_gb', np.float64, ('time','lat','lon'))
esoil_gb.units = 'kg m-2 s-1'
esoil_gb.standard_name = 'Gridbox surface evapotranspiration from soil moisture store'

ecan_gb = ncfile.createVariable('ecan_gb', np.float64, ('time','lat','lon'))
ecan_gb.units = 'kg m-2 s-1'
ecan_gb.standard_name = 'Gridbox mean evaporation from canopy/surface store'

perc = ncfile.createVariable('perc', np.float64, ('time','lat','lon'))
perc.units = "%"
perc.standard_name = "Percentage of glacier area per gridcell area"

rain = ncfile.createVariable('rain', np.float64, ('time','lat','lon'))
rain.units = 'kg m-2 s-1'
rain.standard_name = "Rainfall rate over glacier area contributing directly to runoff"

area = ncfile.createVariable('area', np.float64, ('time','lat','lon'))
area.units = 'm2'
area.standard_name = "Gridbox glacier area"

precip = ncfile.createVariable('precip', np.float64, ('time','lat','lon'))
precip.units = "kg m-2 s-1"
precip.standard_name = "Gridbox precipitation rate"

rainfall = ncfile.createVariable('rainfall', np.float64, ('time', 'lat', 'lon'))
rainfall.units = "kg m-2 s-1"
rainfall.standard_name = "Gridbox rainfall rate"

snowfall = ncfile.createVariable('snowfall', np.float64, ('time','lat','lon'))
snowfall.units = "kg m-2 s-1"
snowfall.standard_name = "Gridbox snowfall rate"

nlats = len(lat_dim)
nlons = len(lon_dim)
ntimes = len(time_dim)
lat[:] = LAT
lon[:] = LON
melt[:,:,:] = MELT[:,:,:]
perc[:,:,:] = PERC[:,:,:]
sub_surf_roff[:,:,:]=SUB_SURF_ROFF[:,:,:]*(1-PERC[:,:,:])
surf_roff[:,:,:]=SURF_ROFF[:,:,:]*(1-PERC[:,:,:]) + MELT[:,:,:] + RAIN[:,:,:]
runoff[:,:,:]=surf_roff[:,:,:]+sub_surf_roff[:,:,:]
et0[:,:,:]=ET0[:,:,:]
ei_gb[:,:,:]=EI_GB[:,:,:]
elake[:,:,:]=ELAKE[:,:,:]
ecan_gb[:,:,:]=ECAN_GB[:,:,:]
esoil_gb[:,:,:]=ESOIL_GB[:,:,:]
rain[:,:,:]=RAIN[:,:,:]
area[:,:,:]=AREA[:,:,:]
precip[:,:,:]=PRECIP[:,:,:]
rainfall[:,:,:]=RAINFALL[:,:,:]
snowfall[:,:,:]=PRECIP[:,:,:]-RAINFALL[:,:,:]
time[:] = TIME
ncfile.close()

daily_jules.close()
perc.close()
daily_area.close()
daily_melt.close()
daily_rain.close()
nc.close()
daily_meteo.close()
