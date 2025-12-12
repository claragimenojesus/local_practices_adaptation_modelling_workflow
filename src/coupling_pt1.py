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

perc_y = pyreadr.read_r(SCENARIO_SH+'/percentage_glacier_per_gridcell_area.rda')
perc_y = perc_y['perc']
perc_y = perc_y.rename(({'dim_0':'time','dim_1':'lat','dim_2':'lon'}))

area_y = pyreadr.read_r(SCENARIO_SH+'/yearly_areas.rda')
area_y = area_y['DF']
area_y = area_y.rename(({'dim_0':'time', 'dim_1':'lat', 'dim_2':'lon'}))

# Create nc file for % of glacier area per gridcell
AREA = area_y
nc = xr.open_dataset(os.path.join(SCENARIO_SH_OUTPUT,"JULES_vn6.1.S2.daily_hydrology.2019.2D.nc"))
LAT = nc['lat']
LAT = LAT.reindex(lat=list(reversed(LAT.lat)))
LON = nc['lon']
dates = pd.date_range(start=dt.datetime(2000,1,1,0),periods=101,freq='AS-JAN') # yearly datetime 1st of jan of 2000 to 2100
dates = dates.to_pydatetime()
ncfile = Dataset(SCENARIO_SH+"/glacier_area_yyyy.nc", mode="w", format="NETCDF4")
lat_dim = ncfile.createDimension('lat', 75)
lon_dim = ncfile.createDimension('lon', 102)
time_dim = ncfile.createDimension('time', 101)
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
area = ncfile.createVariable('area', np.float64, ('time','lat','lon'))
area.units = 'm2'
area.standard_name = 'Glacier area'
nlats = len(lat_dim)
nlons = len(lon_dim)
ntimes = len(time_dim)
lat[:] = LAT
lon[:] = LON
area[:,:,:] = AREA
time[:] = TIME
ncfile.close()
nc.close()

# Create nc file for % of glacier area per gridcell
PERC = perc_y
nc = xr.open_dataset(os.path.join(SCENARIO_SH_OUTPUT, "JULES_vn6.1.S2.daily_hydrology.2019.2D.nc"))
LAT = nc['lat']
LAT = LAT.reindex(lat=list(reversed(LAT.lat)))
LON = nc['lon']
dates = pd.date_range(start=dt.datetime(2000,1,1,0),periods=101,freq='AS-JAN')
dates = dates.to_pydatetime()
ncfile = Dataset(SCENARIO_SH+"/glacier_perc_per_gridcell.nc", mode="w", format="NETCDF4")
lat_dim = ncfile.createDimension('lat', 75)
lon_dim = ncfile.createDimension('lon', 102)
time_dim = ncfile.createDimension('time', 101)
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
percentage = ncfile.createVariable('percentage', np.float64, ('time','lat','lon'))
percentage.units = '%'
percentage.standard_name = 'glacier_percentage_per_gricell_area'
nlats = len(lat_dim)
nlons = len(lon_dim)
ntimes = len(time_dim)
lat[:] = LAT
lon[:] = LON
percentage[:,:,:] = PERC
time[:] = TIME
ncfile.close()
nc.close()
