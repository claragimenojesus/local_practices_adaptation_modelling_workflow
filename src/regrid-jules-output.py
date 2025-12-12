#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import numpy as np
import pandas as pd
import xarray
import click
import yaml
# for development:
# from importlib import reload
from constants import OUTPUT_VARS 

START = int(os.environ['JULES_START_YEAR'])
END = int(os.environ['JULES_END_YEAR'])
YEARS = np.arange(START, END + 1)
SUITE = str(os.environ['JULES_SUITE'])
ID_STEM = str(os.environ['JULES_ID_STEM'])
JOB_NAME = str(os.environ['JULES_JOB_NAME'])
PROFILE_NAME = str(os.environ['JULES_PROFILE_NAME'])
GRIDFILE = str(os.environ['JULES_GRIDFILE'])
DATADIR = str(os.environ['JULES_OUTPUT_DIRECTORY'])
OUTDIR = str(os.environ['REGRID_OUTPUT_DIRECTORY'])
FILE_SUFFIX = str(os.environ['REGRID_FILE_SUFFIX'])
Y_DIM_NAME = str(os.environ['JULES_Y_DIM_NAME'])
X_DIM_NAME = str(os.environ['JULES_X_DIM_NAME'])
MASK_VAR_NAME = str(os.environ['JULES_MASK_VAR_NAME'])
SOIL_DIM_NAME = str(os.environ['JULES_SOIL_DIM_NAME'])
TILE_DIM_NAME = str(os.environ['JULES_TILE_DIM_NAME'])
PFT_DIM_NAME = str(os.environ['JULES_PFT_DIM_NAME'])

# Open 2D land fraction file
y = xarray.open_dataset(GRIDFILE)
MASK = y[MASK_VAR_NAME].values
LAT = y[Y_DIM_NAME].values[:]
LON = y[X_DIM_NAME].values[:]
NLAT = len(LAT)
NLON = len(LON)
# Find out whether latitudes are N-S or S-N
NS_LAT = LAT[0] > LAT[1] 
y.close()

def aggregate_to_month(x):
    # Get the number of days in each month for the current dataset
    days_in_month = xarray.DataArray(
        np.ones(x.time.shape[0]),
        {'time' : x.time.values},
        dims=['time']
    ).groupby('time.month').sum(dim='time')
    x_month = (
        x.groupby('time.month').mean(dim='time')
        * 60 * 60 * 24 * days_in_month
    )
    return x_month

def convert_to_2d(x, variables):

    # Obtain length of dimensions
    if 'time' in x.dims:
        ntime = x.dims['time']
    if 'month' in x.dims:
        ntime = x.dims['month']
    if 'tile' in x.dims:
        ntile = x.dims['tile']
    if 'soil' in x.dims:                
        nsoil = x.dims['soil']
    if 'pft' in x.dims:
        npft = x.dims['pft']
    if 'soilt' in x.dims:
        nsoilt = x.dims['soilt']

    # Update mask to handle subgrid
    # This method hinges on np.isclose(...) being able to
    # reliably match coordinates 
    x_lat = x.latitude.values.squeeze()
    x_lon = x.longitude.values.squeeze()
    x_mask = np.zeros_like(MASK)
    for i in range(len(x_lat)):
        latv = x_lat[i]
        lonv = x_lon[i]
        try:
            latv_ix = np.where(np.isclose(LAT, latv))[0]
            lonv_ix = np.where(np.isclose(LON, lonv))[0]
        except:
            raise ValueError("1")
        if (len(latv_ix) != 1) | (len(lonv_ix) != 1):
            raise ValueError("2")
        x_mask[latv_ix, lonv_ix] = 1        
    
    # Get IDs of cells in the 1D JULES output
    arr = np.arange(NLAT * NLON).reshape((NLAT, NLON))
    if NS_LAT:
        arr = np.flipud(arr)
    y_index = arr[x_mask > 0]

    # Get coordinate pairs for grid
    lat_index = np.array([np.ones((NLON)) * i for i in range(NLAT)], dtype=int).flatten()
    lon_index = np.tile(np.arange(NLON), NLAT)    
    indices = np.array([[lat_index[i], lon_index[i]] for i in y_index], dtype=int)    
    def get_coords(soil=False, soilt=False, tile=False, pft=False):
        coords={}
        dims=[]
        if 'time' in x.dims:
            coords['time'] = x.time.values[:]
            dims.append('time')
        elif 'month' in x.dims:
            coords['month'] = x.month.values[:]
            dims.append('month')            
        if soil:
            coords[SOIL_DIM_NAME] = np.arange(1, nsoil + 1)
            dims.append(SOIL_DIM_NAME)
        if soilt:
            coords['soilt'] = x.soilt.values[:]
            dims.append('soilt')
        if tile:
            coords[TILE_DIM_NAME] = np.arange(1, ntile + 1)
            dims.append(TILE_DIM_NAME)
        if pft:
            coords[PFT_DIM_NAME] = np.arange(1, npft + 1)
            dims.append(PFT_DIM_NAME)
        coords['lat'] = LAT
        coords['lon'] = LON
        dims = dims + ['lat', 'lon']
        return coords, dims
    
    # Define functions to handle outputs with different dimensions
    def convert_gridbox_soilt_var(output_var):
        # (time, soilt, y, x)
        arr = np.zeros((ntime, nsoilt) + (NLAT, NLON)) * np.nan
        for i in range(ntime):
            for j in range(nsoilt):
                output_arr = output_var.values[:][i, j, ...].squeeze()
                arr_ij = arr[i, j, ...]
                arr_ij[tuple(indices.T)] = output_arr
                arr[i, j, ...] = arr_ij
        coords, dims = get_coords(soilt=True)
        xarr = xarray.DataArray(arr, coords, dims, attrs=output_var.attrs)
        return xarr

    def convert_gridbox_var(output_var):
        # (time, y, x)
        arr = np.zeros((ntime,) + (NLAT, NLON)) * np.nan
        for i in range(ntime):
            output_arr = output_var.values[:][i, ...].squeeze()
            arr_ij = arr[i, ...]
            arr_ij[tuple(indices.T)] = output_arr
            arr[i, ...] = arr_ij                
        # TODO: add mask?
        coords, dims = get_coords()
        xarr = xarray.DataArray(arr, coords, dims, attrs=output_var.attrs)
        return xarr

    def convert_tile_var(output_var):
        # (time, tile, y, x)
        arr = np.zeros((ntime, ntile) + (NLAT, NLON)) * np.nan
        for i in range(ntime):
            for j in range(ntile):
                output_arr = output_var.values[:][i, j, ...].squeeze()
                arr_ij = arr[i, j, ...]
                arr_ij[tuple(indices.T)] = output_arr
                arr[i, j, ...] = arr_ij
        coords, dims = get_coords(tile=True)
        xarr = xarray.DataArray(arr, coords, dims, attrs=output_var.attrs)
        return xarr

    def convert_soil_soilt_var(output_var):
        # (time, soil, soilt, y, x)
        arr = np.zeros((ntime, nsoil, nsoilt) + (NLAT, NLON)) * np.nan
        for i in range(ntime):
            for j in range(nsoil):
                for k in range(nsoilt):                            
                    output_arr = output_var.values[:][i, j, k, ...].squeeze()
                    arr_ijk = arr[i, j, k, ...]
                    arr_ijk[tuple(indices.T)] = output_arr
                    arr[i, j, k, ...] = arr_ijk
        coords, dims = get_coords(soil=True, soilt=True)
        xarr = xarray.DataArray(arr, coords, dims, attrs=output_var.attrs)
        return xarr

    def convert_soil_var(output_var):
        # (time, soil, y, x)
        arr = np.zeros((ntime, nsoil) + (NLAT, NLON)) * np.nan
        for i in range(ntime):
            for j in range(nsoil):
                output_arr = output_var.values[:][i, j, ...].squeeze()
                arr_ij = arr[i, j, ...]
                arr_ij[tuple(indices.T)] = output_arr
                arr[i, j, ...] = arr_ij
        coords, dims = get_coords(soil=True)
        xarr = xarray.DataArray(arr, coords, dims, attrs=output_var.attrs)
        return xarr

    def convert_pft_var(output_var):
        # (time, pft, y, x)
        arr = np.zeros((ntime, npft) + (NLAT, NLON)) * np.nan
        for i in range(ntime):
            for j in range(npft):
                output_arr = output_var.values[:][i, j, ...].squeeze()
                arr_ij = arr[i, j, ...]
                arr_ij[tuple(indices.T)] = output_arr
                arr[i, j, ...] = arr_ij
        coords, dims = get_coords(pft=True)
        xarr = xarray.DataArray(arr, coords, dims, attrs=output_var.attrs)
        return xarr

    # Loop through variables and create a list of DataArrays
    xarr_list = []
    for var in variables:
        if var in x.variables:
            if 'tile' in x[var].dims:
                xarr = convert_tile_var(x[var])
            elif 'soil' in x[var].dims:
                if 'soilt' in x[var].dims:
                    xarr = convert_soil_soilt_var(x[var])
                else:
                    xarr = convert_soil_var(x[var])
            elif 'pft' in x[var].dims:
                xarr = convert_pft_var(x[var])
            else:
                if 'soilt' in x[var].dims:
                    xarr = convert_gridbox_soilt_var(x[var])
                else:
                    xarr = convert_gridbox_var(x[var])
            xarr.name = var
            xarr_list.append(xarr)

    # Merge DataArray objects into single dataset
    ds = xarray.merge(xarr_list)
    return ds

def main():
    for yr in YEARS:
        print('Processing output profile ' + PROFILE_NAME + ' for year ' + str(yr) + '...')
        # Open JULES 1D output file
        job_name = JOB_NAME.format(year=yr)
        FN = ID_STEM + '.' + job_name + '.' + PROFILE_NAME + '.' + str(yr) + '.nc'
        x = xarray.open_dataset(os.path.join(DATADIR, FN))
        ds = convert_to_2d(x, OUTPUT_VARS[PROFILE_NAME])
        ds['lat'].attrs['standard_name'] = 'latitude'
        ds['lat'].attrs['units'] = 'degrees_north'
        ds['lon'].attrs['standard_name'] = 'longitude'
        ds['lon'].attrs['units'] = 'degrees_east'
        ds.to_netcdf(
            os.path.join(OUTDIR, os.path.splitext(FN)[0] + '.' + FILE_SUFFIX + '.nc'),
            format="NETCDF4"
        )
        x.close()
    
if __name__ == '__main__':
    main()
