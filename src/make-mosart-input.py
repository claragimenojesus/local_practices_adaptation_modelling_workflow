#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import geopandas as gpd
import pyflwdir
import rasterio
from rasterio import features
from rasterio.transform import Affine
import numpy as np
import netCDF4
import xarray as xr

# local libraries
from pyflwdir.gis_utils import affine_to_coords
from gislib.utils import xy_to_subidx, clip_bbox_global, build_vrt, vrt_props
from gislib import rio

# Default fill vals for netCDF 
F8_FILLVAL = netCDF4.default_fillvals['f8']
F4_FILLVAL = netCDF4.default_fillvals['f4']
I4_FILLVAL = netCDF4.default_fillvals['i4']
I8_FILLVAL = netCDF4.default_fillvals['i8']

AUX_DATADIR = str(os.environ['AUX_DATADIR'])
OUTPUT_DIRECTORY = str(os.environ['OUTPUT_DIRECTORY'])
MEAN_ANNUAL_RUNOFF = str(os.environ['MEAN_ANNUAL_RUNOFF'])
MEAN_ANNUAL_RUNOFF_VARNAME = str(os.environ['MEAN_ANNUAL_RUNOFF_VARNAME'])
MOSART_RES = str(os.environ['RES'])
                 
dir_file = os.path.join(AUX_DATADIR, MOSART_RES + '_flwdir_subrgn.tif')
lat_file = os.path.join(AUX_DATADIR, MOSART_RES + '_outlat_subrgn.tif')
lon_file = os.path.join(AUX_DATADIR, MOSART_RES + '_outlon_subrgn.tif')
# wth_file = os.path.join(AUX_DATADIR, MOSART_RES + '_rivwth_subrgn.tif')
len_file = os.path.join(AUX_DATADIR, MOSART_RES + '_rivlen_subrgn.tif')
elv_file = os.path.join(AUX_DATADIR, MOSART_RES + '_elevtn_subrgn.tif')
slp_file = os.path.join(AUX_DATADIR, MOSART_RES + '_rivslp_subrgn.tif')
hslp_file = os.path.join(AUX_DATADIR, 'merit_slope_' + MOSART_RES + '.tif')

with rasterio.open(dir_file, 'r') as src:
    flwdir = src.read(1)
    transform = src.transform
    crs = src.crs
    latlon = crs.to_epsg() == 4326

with rasterio.open(lat_file, 'r') as src:
    outlat = src.read(1)
    
with rasterio.open(lon_file, 'r') as src:
    outlon = src.read(1)
    
# with rasterio.open(wth_file, 'r') as src:
#     rivwth = src.read(1)

with rasterio.open(len_file, 'r') as src:
    rivlen = src.read(1)

with rasterio.open(elv_file, 'r') as src:
    elevtn = src.read(1)

with rasterio.open(slp_file, 'r') as src:
    rivslp = src.read(1)

with rasterio.open(hslp_file, 'r') as src:
    hslp = src.read(1)

flwdir = pyflwdir.from_array(
    flwdir, ftype='d8', transform=transform, latlon=latlon, cache=True
)
    
OUTLET_X = float(os.environ['OUTLET_X'])
OUTLET_Y = float(os.environ['OUTLET_Y'])
OUTLET_COORDS = (OUTLET_X, OUTLET_Y)

def main():

    # ================================== #
    # 1 - Flow direction
    # ================================== #
    
    # Retrieve validity flags - TODO
    uparea = flwdir.upstream_area(unit='m2')            
    # Extract river basin at coarse resolution, then use as a mask
    basins = flwdir.basins(xy=OUTLET_COORDS)
    # Get global cell id and downstream id
    ids = np.array([i for i in range(flwdir.size)]).reshape(flwdir.shape)
    ids = np.array(
        [i for i in range(flwdir.size)], dtype=np.int32
    ).reshape(flwdir.shape)
    # Start from 1, not zero
    ids += 1
    # Flip the array so that the index starts from the bottom left,
    # increasing from left to right fastest.
    ids = np.flipud(ids)
    # Get the flow direction map in terms of the NEXTXY global format
    nextxy = flwdir.to_array('nextxy')
    nextx = nextxy[0,...]
    nexty = nextxy[1,...]
    ny, nx = flwdir.shape    
    # Preallocate output
    dn_id = np.zeros((flwdir.shape))
    for i in range(ny):
        for j in range(nx):
            # account for zero indexing
            yi = nexty[i,j] - 1  
            xi = nextx[i,j] - 1
            if (yi >= 0) & (xi >= 0):
                idx = ids[yi,xi]
            else:
                idx = -9999
            dn_id[i,j] = idx

    # ================================== #
    # 2 - Area
    # ================================== #

    # Get longitude/latitude coordinates
    lon_vals, lat_vals = affine_to_coords(transform, flwdir.shape)

    # Compute grid area in radians
    area_rad = np.zeros((ny, nx), dtype=np.float64)
    area_m2 = np.zeros((ny, nx), dtype=np.float64)
    R = 6371007.2               # Radius of the earth
    for i in range(len(lon_vals)):
        lon0 = (lon_vals[i] - transform[0] / 2) * (np.pi / 180)
        lon1 = (lon_vals[i] + transform[0] / 2) * (np.pi / 180)
        for j in range(len(lat_vals)):
            lat0 = (lat_vals[j] + transform[4] / 2) * (np.pi / 180)
            lat1 = (lat_vals[j] - transform[4] / 2) * (np.pi / 180)
            area_rad[j,i] = (np.sin(lat1) - np.sin(lat0)) * (lon1 - lon0)
            area_m2[j,i] = (np.sin(lat1) - np.sin(lat0)) * (lon1 - lon0) * R ** 2

    # 1 = land; 0 = ocean
    land_frac = ~np.isnan(elevtn)
    
    # ================================== #
    # 3 - Accumulate mean annual Q 
    # ================================== #

    # We use this to estimate some channel parameters, based on
    # empirical relationships

    # TODO: check units
    runoff = xr.open_dataset(MEAN_ANNUAL_RUNOFF)[MEAN_ANNUAL_RUNOFF_VARNAME]
    # Convert m/y to m3/y
    runoff *= area_m2
    # Convert m3/y to m3/s
    runoff /= (365 * 24 * 60 * 60)
    Qmean = flwdir.accuflux(runoff, direction="up")
    Qmean = Qmean.astype(np.float64)
    
    # ================================== #
    # 3 - Compute channel bankfull depth/width
    # ================================== #
    
    # HyMAP (https://journals.ametsoc.org/view/journals/hydr/13/6/jhm-d-12-021_1.xml#bib8)
    Beta = 18.                  # TODO: check this value
    alpha = 3.73e-3
    width_main = np.clip(Beta * Qmean ** 0.5, 10., None).astype(np.float64)
    depth_main = np.clip(alpha * width_main, 2., None).astype(np.float64)
    width_main_floodplain = width_main * 3.

    # Alternatives
    # NB first two from https://ec-jrc.github.io/lisflood/pdfs/Dataset_hydro.pdf
    # width_main = uparea ** 0.0032    
    # width_main = Qmean_lr ** 0.539
    # width_main = flw.subgrid_rivavg(idxs_out, rivwth)
    
    # Tributary bankfull depth/width (use gridbox runoff)
    width_trib = np.clip(
        Beta * runoff.values ** 0.5, 10., None
    ).astype(np.float64)
    depth_trib = np.clip(
        alpha * width_trib, 2., None
    ).astype(np.float64)

    # Manning channel/overland (LISFLOOD)
    n_channel = (
        0.025 + 0.015
        * np.clip(50. / (uparea / 1000 / 1000), None, 1)
        + 0.030 * np.clip(elevtn / 2000., None, 1.)
    ).astype(np.float64)
    
    # Initial guess
    n_overland = (np.ones_like(n_channel) * 0.03).astype(np.float64)
    
    # ================================== #
    # 4 - Compute river length/slope
    # ================================== #

    length_main = rivlen
    slope_main = rivslp
    # MOSART makes tributary slope equal to main river slope
    slope_trib = slope_main
    # This is not implemented yet:    
    # rivslp2 = flw.subgrid_rivslp2(idxs_out, elevtn)
    
    # ================================== #
    # 5 - Write output
    # ================================== #

    mask = ~(basins.astype(bool))
    
    # First output file defines the model grid
    nco = netCDF4.Dataset(
        os.path.join(OUTPUT_DIRECTORY, 'land.nc'), 'w', format='NETCDF4'
    )
    nco.createDimension('lat', len(lat_vals))
    nco.createDimension('lon', len(lon_vals))

    var = nco.createVariable(
        'lon', 'f8', ('lon',)
    )
    var.units = 'degrees_east'
    var.standard_name = 'lon'
    var.long_name = 'longitude'
    var[:] = lon_vals.astype(np.float64)

    var = nco.createVariable(
        'lat', 'f8', ('lat',)
    )
    var.units = 'degrees_north'
    var.standard_name = 'lat'
    var.long_name = 'latitude'
    var[:] = lat_vals.astype(np.float64)

    # 0 = ocean; 1 = land.
    var = nco.createVariable('mask', 'i4', ('lat', 'lon'), fill_value=I4_FILLVAL)
    var.units = '1'
    var.long_name = 'land domain mask'
    var[:] = basins.astype(np.int32)

    var = nco.createVariable('frac', 'f8', ('lat', 'lon'), fill_value=F8_FILLVAL)
    var.units = '1'
    var.long_name = 'fraction of grid cell that is active'
    land_frac_masked = np.ma.masked_array(
        land_frac, mask=mask, dtype=np.float64
    )
    var[:] = land_frac_masked

    var = nco.createVariable('area', 'f8', ('lat', 'lon'), fill_value=F8_FILLVAL)
    var.units = 'radians^2'
    var.long_name = 'area of grid cell in radians squared'
    area_rad_masked = np.ma.masked_array(
        area_rad, mask=mask, dtype=np.float64
    )
    var[:] = area_rad_masked
    nco.close()

    # Second output file contains the model parameters
    nco = netCDF4.Dataset(
        os.path.join(OUTPUT_DIRECTORY, 'mosart.nc'), 'w', format='NETCDF4'
    )
    nco.createDimension('lat', len(lat_vals))
    nco.createDimension('lon', len(lon_vals))

    var = nco.createVariable(
        'lon', 'f8', ('lon',)
    )
    var.units = 'degrees_east'
    var.standard_name = 'lon'
    var.long_name = 'longitude'
    var[:] = lon_vals.astype(np.float64)

    var = nco.createVariable(
        'lat', 'f8', ('lat',)
    )
    var.units = 'degrees_north'
    var.standard_name = 'lat'
    var.long_name = 'latitude'
    var[:] = lat_vals.astype(np.float64)

    var = nco.createVariable('ID', 'f8', ('lat', 'lon'), fill_value=F8_FILLVAL)
    var.units = '1'
    var[:] = ids

    var = nco.createVariable('dnID', 'f8', ('lat', 'lon'), fill_value=F8_FILLVAL)
    var.units = '1'
    dn_id_masked = dn_id
    dn_id_masked[~(basins>0)] = -9999.  # why do this?
    var[:] = dn_id_masked

    var = nco.createVariable('fdir', 'f8', ('lat', 'lon'), fill_value=F8_FILLVAL)
    var.units = '1'
    var.long_name = 'flow direction based on D8 algorithm'
    fdir = flwdir.to_array('d8')
    fdir_masked = np.ma.masked_array(
        fdir, mask=mask, dtype=np.float64
    )
    var[:] = fdir_masked

    var = nco.createVariable('frac', 'f8', ('lat', 'lon'), fill_value=F8_FILLVAL)
    var.units = '1'
    var.long_name = 'fraction of the unit draining to the outlet'
    var[:] = basins
    basins_masked = np.ma.masked_array(
        basins, mask=mask, dtype=np.float64
    )
    var[:] = basins_masked

    # TODO: check slope units
    var = nco.createVariable('rslp', 'f8', ('lat', 'lon'), fill_value=F8_FILLVAL)
    var.units = '1'
    var.long_name = 'main channel slope'
    slope_main_masked = np.ma.masked_array(
        slope_main, mask=mask, dtype=np.float64
    )
    var[:] = slope_main_masked
    
    var = nco.createVariable('rlen', 'f8', ('lat', 'lon'), fill_value=F8_FILLVAL)
    var.units = 'm'
    var.long_name = 'main channel length'
    length_main_masked = np.ma.masked_array(
        length_main, mask=mask, dtype=np.float64
    )
    var[:] = length_main_masked

    var = nco.createVariable('rdep', 'f8', ('lat', 'lon'), fill_value=F8_FILLVAL)
    var.units = 'm'
    var.long_name = 'bankfull depth of main channel'
    depth_main_masked = np.ma.masked_array(
        depth_main, mask=mask, dtype=np.float64
    )
    var[:] = depth_main_masked

    var = nco.createVariable('rwid', 'f8', ('lat', 'lon'), fill_value=F8_FILLVAL)
    var.units = 'm'
    var.long_name = 'bankfull width of main channel'
    width_main_masked = np.ma.masked_array(
        width_main, mask=mask, dtype=np.float64
    )
    var[:] = width_main_masked

    var = nco.createVariable('rwid0', 'f8', ('lat', 'lon'), fill_value=F8_FILLVAL)
    var.units = 'm'
    var.long_name = 'floodplain width linked to main channel'
    width_main_floodplain_masked = np.ma.masked_array(
        width_main_floodplain, mask=mask, dtype=np.float64
    )
    var[:] = width_main_floodplain_masked

    # TODO: drainage density
    var = nco.createVariable('gxr', 'f8', ('lat', 'lon'), fill_value=F8_FILLVAL)
    var.units = 'm^-1'
    var.long_name = 'drainage density'
    gxr = np.full((ny, nx), 0.5, dtype=np.float64)
    gxr_masked = np.ma.masked_array(
        gxr, mask=mask, dtype=np.float64
    )
    var[:] = gxr_masked

    var = nco.createVariable('hslp', 'f8', ('lat', 'lon'), fill_value=F8_FILLVAL)
    var.units = '1'
    var.long_name = 'topographic slope'
    hslp_masked = np.ma.masked_array(
        hslp, mask=mask, dtype=np.float64
    )
    var[:] = hslp_masked

    var = nco.createVariable('twid', 'f8', ('lat', 'lon'), fill_value=F8_FILLVAL)
    var.units = 'm'
    var.long_name = 'bankfull width of local tributaries'
    width_trib_masked = np.ma.masked_array(
        width_trib, mask=mask, dtype=np.float64
    )
    var[:] = width_trib_masked

    var = nco.createVariable('tslp', 'f8', ('lat', 'lon'), fill_value=F8_FILLVAL)
    var.units = '1'
    var.long_name = 'mean tributary channel slope averaged through the unit'
    slope_trib_masked = np.ma.masked_array(
        slope_trib, mask=mask, dtype=np.float64
    )
    var[:] = slope_trib_masked

    var = nco.createVariable('area', 'f8', ('lat', 'lon'), fill_value=F8_FILLVAL)
    var.units = 'm^2'
    var.long_name = 'local drainage area'
    area_m2_masked = np.ma.masked_array(
        area_m2, mask=mask, dtype=np.float64
    )
    var[:] = area_m2_masked

    # NB This variable is not actually used in MOSART  (new feature?)
    var = nco.createVariable('areaTotal', 'f8', ('lat', 'lon'), fill_value=F8_FILLVAL)
    var.units = 'm^2'
    var.long_name = 'total upstream drainage area including local unit: multi-flow direction'
    uparea_masked = np.ma.masked_array(
        uparea, mask=mask, dtype=np.float64
    )
    var[:] = uparea_masked

    var = nco.createVariable('areaTotal2', 'f8', ('lat', 'lon'), fill_value=F8_FILLVAL)
    var.units = 'm^2'
    var.long_name = 'total upstream drainage area including local unit: single-flow direction'
    var[:] = uparea_masked
    
    var = nco.createVariable('nr', 'f8', ('lat', 'lon'), fill_value=F8_FILLVAL)
    var.units = '1'
    var.long_name = 'Manning\'s roughness coefficient for main channel flow'
    n_channel_masked = np.ma.masked_array(
        n_channel, mask=mask, dtype=np.float64
    )
    var[:] = n_channel_masked

    var = nco.createVariable('nt', 'f8', ('lat', 'lon'), fill_value=F8_FILLVAL)
    var.units = '1'
    var.long_name = 'Manning\'s roughness coefficient for tributary channel flow'    
    var[:] = n_channel_masked

    var = nco.createVariable('nh', 'f8', ('lat', 'lon'), fill_value=F8_FILLVAL)
    var.units = '1'
    var.long_name = 'Manning\'s roughness coefficient for overland flow'
    n_overland_masked = np.ma.masked_array(
        n_overland, mask=mask, dtype=np.float64
    )
    var[:] = n_overland_masked    
    nco.close()

##### If using reduced grid setup (suite u-cz655) where the latlon are inverted, then modify the grid for land.nc and mosart.nc files by running the following: (CGJ 11/01/2024)
land_file = xr.open_dataset(os.path.join(OUTPUT_DIRECTORY, 'land.nc'))
land_file_new = land_file.reindex(lat=list(reversed(land_file.lat)))
land_file_new.to_netcdf(os.path.join(OUTPUT_DIRECTORY, 'land_new.nc'), 'w', format="NETCDF4")

mosart_file = xr.open_dataset(os.path.join(OUTPUT_DIRECTORY, 'mosart.nc'))
mosart_file_new = mosart_file.reindex(lat=list(reversed(mosart_file.lat)))
mosart_file_new.to_netcdf(os.path.join(OUTPUT_DIRECTORY, 'mosart_new.nc'), 'w', format="NETCDF4")


#### end of added code by CGJ
        
if __name__ == "__main__":
    main()

