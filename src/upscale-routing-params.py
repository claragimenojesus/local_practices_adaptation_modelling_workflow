#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Flow direction data

# The `FlwdirRaster` object is at the core of the pyflwdir package.
# It contains gridded flow direction data, parsed to an actionable common format
# which describes the linear index of the next dowsntream cell.

# Currently `pyflwdir` supports two local flow direction (D8) data types
# according to the arcgis **D8** convention and pcraster **LDD** convention one
# global flow direction type according to the CaMa-Flood **NEXTXY** convention.
# Local flow direction data types describe the next downstream cell based on a
# local direction towards one of its neighboring cells, while global flow
# direction types describe the next downstream cell based on a global index.

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

# Default fill vals for netCDF 
F8_FILLVAL = netCDF4.default_fillvals['f8']
F4_FILLVAL = netCDF4.default_fillvals['f4']
I4_FILLVAL = netCDF4.default_fillvals['i4']
I8_FILLVAL = netCDF4.default_fillvals['i8']

AUX_DATADIR = str(os.environ['AUX_DATADIR'])
OUTPUT_DIRECTORY = str(os.environ['OUTPUT_DIRECTORY'])
MEAN_ANNUAL_RUNOFF = str(os.environ['MEAN_ANNUAL_RUNOFF'])

# We read the flow direction and elevation raster data, including meta-data,
# using, [rasterio](https://rasterio.readthedocs.io/en/latest/)
dir_file = os.path.join(AUX_DATADIR, 'merit_rgn_dir.tif')
# wth_file = os.path.join(AUX_DATADIR, 'merit_rgn_wth.tif')
# elv_file = os.path.join(AUX_DATADIR, 'merit_rgn_elv.tif')
# elv_lr_file = os.path.join(AUX_DATADIR, 'merit_rgn_elv_lr.tif')
# frc_lr_file = os.path.join(AUX_DATADIR, 'merit_rgn_land_fraction_lr.tif')
# # slp_lr_file = os.path.join(AUX_DATADIR, 'merit_basin_slope_lr.tif')

with rasterio.open(dir_file, 'r') as src:
    flwdir = src.read(1)
    transform = src.transform
    crs = src.crs
    latlon = crs.to_epsg() == 4326

# with rasterio.open(wth_file, 'r') as src:
#     rivwth = src.read(1)

# with rasterio.open(elv_file, 'r') as src:
#     elevtn = src.read(1)

# with rasterio.open(elv_lr_file, 'r') as src:
#     elevtn_lr = src.read(1)

# with rasterio.open(frc_lr_file, 'r') as src:
#     frac_lr = src.read(1)

# # with rasterio.open(slp_lr_file, 'r') as src:
# #     slp_lr = src.read(1)

# This is a requirement:
flwdir = flwdir.astype('uint8')

flw = pyflwdir.from_array(
    flwdir, ftype='d8', transform=transform, latlon=latlon, cache=True
)

SCALE_FACTOR = int(os.environ['SCALE_FACTOR'])
# # SCALE_FACTOR = 120 # 0.0008333 -> 0.1 [TODO: user-specified]
# OUTLET_X = float(os.environ['OUTLET_X'])
# OUTLET_Y = float(os.environ['OUTLET_Y'])
# # OUTLET_X=-72.6271
# # OUTLET_Y=-13.0045
# OUTLET_COORDS = (OUTLET_X, OUTLET_Y)

def main():

    # ================================== #
    # 1 - Upscale flow direction
    # ================================== #
    
    # Apply scale factor to upscale flow direction map
    # uparea = flw.upstream_area(unit='km2')
    flwdir_lr, idxs_out = flw.upscale(
        scale_factor=SCALE_FACTOR,
        # uparea=uparea,
        method='ihu'
    )
    
    # Retrieve validity flags
    valid = flw.upscale_error(flwdir_lr, idxs_out)

    # Write output
    ftype = flw.ftype
    shape_lr = flwdir_lr.shape
    dims=('y','x')
    transform_lr = Affine(
        transform[0] * SCALE_FACTOR, transform[1], transform[2],
        transform[3], transform[4] * SCALE_FACTOR, transform[5]
    )
    lon_vals, lat_vals = affine_to_coords(transform_lr, flwdir_lr.shape)
    coords={'y': lat_vals, 'x': lon_vals}

    # Flow direction map
    da_flwdir = xr.DataArray(
        name='flwdir',
        data=flwdir_lr.to_array(ftype),
        coords=coords,
        dims=dims,
        attrs=dict(
            long_name=f'{ftype} flow direction', 
            _FillValue=flw._core._mv)            
    )
    da_flwdir.rio.set_crs(4326)
    da_flwdir.rio.to_raster(
        os.path.join(
            AUX_DATADIR,
            "merit_flow_direction.tif"
        )
    )

    # Downstream IDs
    da_outidx = xr.DataArray(
        name='outidx',
        data=idxs_out,
        coords=coords,
        dims=dims,
        attrs=dict(
            long_name = 'subgrid outlet index', 
            # _FillValue = flw._core._mv
            _FillValue = I4_FILLVAL
        )
    )
    da_outidx.rio.set_crs(4326)
    da_outidx.rio.to_raster(
        os.path.join(
            AUX_DATADIR,
            "merit_subgrid_outlet_index.tif"
        )
    )

    # # NOT USED:
    
    # # flw = pyflwdir.from_array(
    # #     flwdir_lr.to_array(ftype), ftype='d8', transform=transform_lr, latlon=latlon, cache=True
    # # )
    
    # # # Translate outlet indices to global x,y coordinates
    # # # NB `dir_file` is arbitrary - the routine uses it to
    # # # extract the coordinates of each point
    # # with xr.open_rasterio(elv_file) as template:
    # #     x_out, y_out = template.rio.idx_to_xy(
    # #         idxs_out,
    # #         mask = idxs_out != flw._core._mv
    # #     )
        
    # # # Extract river basin at coarse resolution, then use as a mask
    # # basins_lr = flwdir_lr.basins(xy=OUTLET_COORDS)

    # # Recompute upstream area for the coarse resolution
    # uparea_lr = flwdir_lr.upstream_area(unit='m2')
    
    # # Get global cell id and downstream id    
    # ids = np.array([i for i in range(flwdir_lr.size)]).reshape(flwdir_lr.shape)
    # ids = np.array(
    #     [i for i in range(flwdir_lr.size)], dtype=np.int32
    # ).reshape(flwdir_lr.shape)

    # # Start from 1, not zero
    # ids += 1

    # # Flip the array so that the index starts from the bottom left,
    # # increasing from left to right fastest.
    # ids = np.flipud(ids)

    # # Get the flow direction map in terms of the NEXTXY global format
    # nextxy = flwdir_lr.to_array('nextxy')
    # nextx = nextxy[0,...]
    # nexty = nextxy[1,...]
    # ny, nx = flwdir_lr.shape
    
    # # Preallocate output
    # dn_id = np.zeros((flwdir_lr.shape))
    # for i in range(ny):
    #     for j in range(nx):
    #         # account for zero indexing
    #         yi = nexty[i,j] - 1  
    #         xi = nextx[i,j] - 1
    #         if (yi >= 0) & (xi >= 0):
    #             idx = ids[yi,xi]
    #         else:
    #             idx = -9999
    #         dn_id[i,j] = idx

    # # Write output
    

    # # # ================================== #
    # # # 2 - Area
    # # # ================================== #

    # # # Get longitude/latitude coordinates
    # # transform_lr = Affine(
    # #     transform[0] * SCALE_FACTOR, transform[1], transform[2],
    # #     transform[3], transform[4] * SCALE_FACTOR, transform[5]
    # # )
    # # lon_vals, lat_vals = affine_to_coords(transform_lr, flwdir_lr.shape)

    # # # Compute grid area in radians
    # # area_rad = np.zeros((ny, nx), dtype=np.float64)
    # # area_m2 = np.zeros((ny, nx), dtype=np.float64)
    # # R = 6371007.2               # Radius of the earth
    # # for i in range(len(lon_vals)):
    # #     lon0 = (lon_vals[i] - transform_lr[0] / 2) * (np.pi / 180)
    # #     lon1 = (lon_vals[i] + transform_lr[0] / 2) * (np.pi / 180)
    # #     for j in range(len(lat_vals)):
    # #         lat0 = (lat_vals[j] + transform_lr[4] / 2) * (np.pi / 180)
    # #         lat1 = (lat_vals[j] - transform_lr[4] / 2) * (np.pi / 180)
    # #         area_rad[j,i] = (np.sin(lat1) - np.sin(lat0)) * (lon1 - lon0)
    # #         area_m2[j,i] = (np.sin(lat1) - np.sin(lat0)) * (lon1 - lon0) * R ** 2
    
    # # # ================================== #
    # # # 3 - Accumulate mean annual Q 
    # # # ================================== #

    # # # We use this to estimate some channel parameters, based on
    # # # empirical relationships
    
    # # runoff_lr = xr.open_dataset(MEAN_ANNUAL_RUNOFF)['runoff']
    # # # Convert m/y to m3/y
    # # runoff_lr *= area_m2
    # # # Convert m3/y to m3/s
    # # runoff_lr /= (365 * 24 * 60 * 60)
    # # Qmean_lr = flwdir_lr.accuflux(runoff_lr, direction="up")
    # # Qmean_lr = Qmean_lr.astype(np.float64)
    
    # # # ================================== #
    # # # 3 - Compute channel bankfull depth/width
    # # # ================================== #
    
    # # # HyMAP (https://journals.ametsoc.org/view/journals/hydr/13/6/jhm-d-12-021_1.xml#bib8)
    # # Beta = 18.                  # TODO: check this value
    # # alpha = 3.73e-3
    # # width_main = np.clip(Beta * Qmean_lr ** 0.5, 10., None).astype(np.float64)
    # # depth_main = np.clip(alpha * width_main, 2., None).astype(np.float64)
    # # width_main_floodplain = width_main * 3.

    # # # Alternatives
    # # # NB first two from https://ec-jrc.github.io/lisflood/pdfs/Dataset_hydro.pdf
    # # # width_main = uparea ** 0.0032    
    # # # width_main = Qmean_lr ** 0.539
    # # # width_main = flw.subgrid_rivavg(idxs_out, rivwth)
    
    # # # Tributary bankfull depth/width (use gridbox runoff)
    # # width_trib = np.clip(Beta * runoff_lr.values ** 0.5, 10., None).astype(np.float64)
    # # depth_trib = np.clip(alpha * width_trib, 2., None).astype(np.float64)

    # # # Manning channel/overland (LISFLOOD)
    # # n_channel = (
    # #     0.025 + 0.015
    # #     * np.clip(50. / (uparea_lr / 1000 / 1000), None, 1)
    # #     + 0.030 * np.clip(elevtn_lr / 2000., None, 1.)
    # # ).astype(np.float64)
    
    # # # Initial guess
    # # n_overland = (np.ones_like(n_channel) * 0.03).astype(np.float64)
    
    # # # ================================== #
    # # # 4 - Compute river length/slope
    # # # ================================== #
    
    # # length_main = flw.subgrid_rivlen(idxs_out, direction="up").astype(np.float64)
    # # slope_main = flw.subgrid_rivslp(idxs_out, elevtn).astype(np.float64)
    # # # This is not implemented yet:    
    # # # rivslp2 = flw.subgrid_rivslp2(idxs_out, elevtn)

    # # # MOSART makes tributary slope equal to main river slope
    # # slope_trib = slope_main
    
    # # # ================================== #
    # # # 5 - Write output
    # # # ================================== #

    # # mask = ~(basins_lr.astype(bool))
    
    # # # First output file defines the model grid
    # # nco = netCDF4.Dataset(
    # #     os.path.join(OUTPUT_DIRECTORY, 'land.nc'), 'w', format='NETCDF4'
    # # )
    # # nco.createDimension('lat', len(lat_vals))
    # # nco.createDimension('lon', len(lon_vals))

    # # var = nco.createVariable(
    # #     'lon', 'f8', ('lon',)
    # # )
    # # var.units = 'degrees_east'
    # # var.standard_name = 'lon'
    # # var.long_name = 'longitude'
    # # var[:] = lon_vals.astype(np.float64)

    # # var = nco.createVariable(
    # #     'lat', 'f8', ('lat',)
    # # )
    # # var.units = 'degrees_north'
    # # var.standard_name = 'lat'
    # # var.long_name = 'latitude'
    # # var[:] = lat_vals.astype(np.float64)

    # # # 0 = ocean; 1 = land.
    # # var = nco.createVariable('mask', 'i4', ('lat', 'lon'), fill_value=I4_FILLVAL)
    # # var.units = '1'
    # # var.long_name = 'land domain mask'
    # # basins_lr_masked = basins_lr.astype(np.int32)
    # # var[:] = basins_lr_masked

    # # var = nco.createVariable('frac', 'f8', ('lat', 'lon'), fill_value=F8_FILLVAL)
    # # var.units = '1'
    # # var.long_name = 'fraction of grid cell that is active'
    # # frac_lr_masked = np.ma.masked_array(frac_lr, mask=mask, dtype=np.float64)
    # # var[:] = frac_lr_masked

    # # var = nco.createVariable('area', 'f8', ('lat', 'lon'), fill_value=F8_FILLVAL)
    # # var.units = 'radians^2'
    # # var.long_name = 'area of grid cell in radians squared'
    # # area_rad_masked = np.ma.masked_array(area_rad, mask=mask, dtype=np.float64)
    # # var[:] = area_rad_masked
    # # nco.close()

    # # # Second output file contains the model parameters
    # # nco = netCDF4.Dataset(
    # #     os.path.join(OUTPUT_DIRECTORY, 'mosart.nc'), 'w', format='NETCDF4'
    # # )
    # # nco.createDimension('lat', len(lat_vals))
    # # nco.createDimension('lon', len(lon_vals))

    # # var = nco.createVariable(
    # #     'lon', 'f8', ('lon',)
    # # )
    # # var.units = 'degrees_east'
    # # var.standard_name = 'lon'
    # # var.long_name = 'longitude'
    # # var[:] = lon_vals.astype(np.float64)

    # # var = nco.createVariable(
    # #     'lat', 'f8', ('lat',)
    # # )
    # # var.units = 'degrees_north'
    # # var.standard_name = 'lat'
    # # var.long_name = 'latitude'
    # # var[:] = lat_vals.astype(np.float64)

    # # var = nco.createVariable('ID', 'f8', ('lat', 'lon'), fill_value=F8_FILLVAL)
    # # var.units = '1'
    # # var[:] = ids

    # # var = nco.createVariable('dnID', 'f8', ('lat', 'lon'), fill_value=F8_FILLVAL)
    # # var.units = '1'
    # # dn_id_masked = dn_id
    # # dn_id_masked[~(basins_lr>0)] = -9999.
    # # var[:] = dn_id_masked

    # # var = nco.createVariable('fdir', 'f8', ('lat', 'lon'), fill_value=F8_FILLVAL)
    # # var.units = '1'
    # # var.long_name = 'flow direction based on D8 algorithm'
    # # fdir = flwdir_lr.to_array(flw.ftype)
    # # fdir_masked = np.ma.masked_array(fdir, mask=mask, dtype=np.float64)
    # # var[:] = fdir_masked

    # # # For now we set this to the basin map, but it would be worth
    # # # experimenting to see if having fractional values to account
    # # # for basin boundary pixels would worl.
    # # var = nco.createVariable('frac', 'f8', ('lat', 'lon'), fill_value=F8_FILLVAL)
    # # var.units = '1'
    # # var.long_name = 'fraction of the unit draining to the outlet'
    # # var[:] = basins_lr
    # # basins_lr_masked = np.ma.masked_array(basins_lr, mask=mask, dtype=np.float64)
    # # var[:] = basins_lr_masked

    # # # TODO: check slope units
    # # var = nco.createVariable('rslp', 'f8', ('lat', 'lon'), fill_value=F8_FILLVAL)
    # # var.units = '1'
    # # var.long_name = 'main channel slope'
    # # slope_main_masked = np.ma.masked_array(slope_main, mask=mask, dtype=np.float64)
    # # var[:] = slope_main_masked
    
    # # var = nco.createVariable('rlen', 'f8', ('lat', 'lon'), fill_value=F8_FILLVAL)
    # # var.units = 'm'
    # # var.long_name = 'main channel length'
    # # length_main_masked = np.ma.masked_array(length_main, mask=mask, dtype=np.float64)
    # # var[:] = length_main_masked

    # # var = nco.createVariable('rdep', 'f8', ('lat', 'lon'), fill_value=F8_FILLVAL)
    # # var.units = 'm'
    # # var.long_name = 'bankfull depth of main channel'
    # # depth_main_masked = np.ma.masked_array(depth_main, mask=mask, dtype=np.float64)
    # # var[:] = depth_main_masked

    # # var = nco.createVariable('rwid', 'f8', ('lat', 'lon'), fill_value=F8_FILLVAL)
    # # var.units = 'm'
    # # var.long_name = 'bankfull width of main channel'
    # # width_main_masked = np.ma.masked_array(width_main, mask=mask, dtype=np.float64)
    # # var[:] = width_main_masked

    # # var = nco.createVariable('rwid0', 'f8', ('lat', 'lon'), fill_value=F8_FILLVAL)
    # # var.units = 'm'
    # # var.long_name = 'floodplain width linked to main channel'
    # # width_main_floodplain_masked = np.ma.masked_array(width_main_floodplain, mask=mask, dtype=np.float64)
    # # var[:] = width_main_floodplain_masked

    # # # TODO: drainage density
    # # var = nco.createVariable('gxr', 'f8', ('lat', 'lon'), fill_value=F8_FILLVAL)
    # # var.units = 'm^-1'
    # # var.long_name = 'drainage density'
    # # gxr = np.full((ny, nx), 0.5, dtype=np.float64)
    # # gxr_masked = np.ma.masked_array(gxr, mask=mask, dtype=np.float64)
    # # var[:] = gxr_masked

    # # var = nco.createVariable('hslp', 'f8', ('lat', 'lon'), fill_value=F8_FILLVAL)
    # # var.units = '1'
    # # var.long_name = 'topographic slope'
    # # slp_lr_masked = np.ma.masked_array(slp_lr, mask=mask, dtype=np.float64)
    # # var[:] = slp_lr_masked

    # # var = nco.createVariable('twid', 'f8', ('lat', 'lon'), fill_value=F8_FILLVAL)
    # # var.units = 'm'
    # # var.long_name = 'bankfull width of local tributaries'
    # # width_trib_masked = np.ma.masked_array(width_trib, mask=mask, dtype=np.float64)
    # # var[:] = width_trib_masked

    # # var = nco.createVariable('tslp', 'f8', ('lat', 'lon'), fill_value=F8_FILLVAL)
    # # var.units = '1'
    # # var.long_name = 'mean tributary channel slope averaged through the unit'
    # # slope_trib_masked = np.ma.masked_array(slope_trib, mask=mask, dtype=np.float64)
    # # var[:] = slope_trib_masked

    # # var = nco.createVariable('area', 'f8', ('lat', 'lon'), fill_value=F8_FILLVAL)
    # # var.units = 'm^2'
    # # var.long_name = 'local drainage area'
    # # area_m2_masked = np.ma.masked_array(area_m2, mask=mask, dtype=np.float64)
    # # var[:] = area_m2_masked

    # # # NB This variable is not actually used in MOSART (I guess it may be a new feature?)
    # # var = nco.createVariable('areaTotal', 'f8', ('lat', 'lon'), fill_value=F8_FILLVAL)
    # # var.units = 'm^2'
    # # var.long_name = 'total upstream drainage area including local unit: multi-flow direction'
    # # uparea_lr_masked = np.ma.masked_array(uparea_lr, mask=mask, dtype=np.float64)
    # # var[:] = uparea_lr_masked

    # # var = nco.createVariable('areaTotal2', 'f8', ('lat', 'lon'), fill_value=F8_FILLVAL)
    # # var.units = 'm^2'
    # # var.long_name = 'total upstream drainage area including local unit: single-flow direction'
    # # var[:] = uparea_lr_masked
    
    # # var = nco.createVariable('nr', 'f8', ('lat', 'lon'), fill_value=F8_FILLVAL)
    # # var.units = '1'
    # # var.long_name = 'Manning\'s roughness coefficient for main channel flow'
    # # n_channel_masked = np.ma.masked_array(n_channel, mask=mask, dtype=np.float64)
    # # var[:] = n_channel_masked

    # # var = nco.createVariable('nt', 'f8', ('lat', 'lon'), fill_value=F8_FILLVAL)
    # # var.units = '1'
    # # var.long_name = 'Manning\'s roughness coefficient for tributary channel flow'    
    # # var[:] = n_channel_masked

    # # var = nco.createVariable('nh', 'f8', ('lat', 'lon'), fill_value=F8_FILLVAL)
    # # var.units = '1'
    # # var.long_name = 'Manning\'s roughness coefficient for overland flow'
    # # n_overland_masked = np.ma.masked_array(n_overland, mask=mask, dtype=np.float64)
    # # var[:] = n_overland_masked
    
    # # nco.close()
        
    # # # # setup output DataArray
    # # # ftype = flw.ftype
    # # # shape_lr = flwdir_lr.shape
    # # # dims=('y','x')
    # # # # apply scale factor to the transform
    # # # transform_lr = Affine(
    # # #     transform[0] * scale_factor, transform[1], transform[2],
    # # #     transform[3], transform[4] * scale_factor, transform[5]
    # # # )
    # # # xs, ys = affine_to_coords(transform_lr, shape_lr)    
    # # # coords={'y': ys, 'x': xs}

    # # # # Flow direction map
    # # # # ##################
    # # # da_flwdir = xr.DataArray(
    # # #     name='flwdir',
    # # #     data=flwdir_lr.to_array(ftype),
    # # #     coords=coords,
    # # #     dims=dims,
    # # #     attrs=dict(
    # # #         long_name=f'{ftype} flow direction', 
    # # #         _FillValue=flw._core._mv)            
    # # # )

    # # # # Add subgrid outlets to flow direction map
    # # # da_flwdir.coords['x_out'] = xr.Variable(
    # # #     dims=dims, 
    # # #     data=x_out, 
    # # #     attrs = dict(
    # # #         long_name = 'subgrid outlet x coordinate', 
    # # #         _FillValue = np.nan)
    # # # )    
    # # # da_flwdir.coords['y_out'] = xr.Variable(
    # # #     dims=dims, 
    # # #     data=y_out, 
    # # #     attrs = dict(
    # # #         long_name = 'subgrid outlet y coordinate', 
    # # #         _FillValue = np.nan)
    # # # )
    # # # # Add validity flag to dataset coordinates
    # # # da_flwdir.coords['flwerr'] = xr.Variable(
    # # #     dims=dims,
    # # #     data=valid,
    # # #     attrs=dict(
    # # #         long_name = 'valid flw dirs', 
    # # #         _FillValue = -1
    # # #     )
    # # # )

    # # # # Subgrid outlet index
    # # # # ####################
    # # # da_outidx = xr.DataArray(
    # # #     name='outidx',
    # # #     data=idxs_out,
    # # #     coords=coords,
    # # #     dims=dims,
    # # #     attrs=dict(
    # # #         long_name = 'subgrid outlet index', 
    # # #         _FillValue = flw._core._mv
    # # #     )
    # # # )

    # # # # Global low-resolution index
    # # # # ###########################
    # # # da_ids = xr.DataArray(
    # # #     name='ID',
    # # #     data=ids,
    # # #     coords=coords,
    # # #     dims=dims,
    # # #     attrs=dict(
    # # #         long_name = 'id',
    # # #         _FillValue=I4_FILLVAL
    # # #     )            
    # # # )

    # # # # Global low-resolution downstream IDs
    # # # # ####################################
    # # # da_dnids = xr.DataArray(
    # # #     name='dnID',
    # # #     data=dn_id,
    # # #     coords=coords,
    # # #     dims=dims,
    # # #     attrs=dict(
    # # #         long_name = 'dnID',
    # # #         _FillValue=I4_FILLVAL
    # # #     )            
    # # # )
    
    # # # # River basin map
    # # # # ###############
    # # # da_basins = xr.DataArray(
    # # #     name='basins',
    # # #     data=basins_lr,
    # # #     coords=coords,
    # # #     dims=dims,
    # # #     attrs=dict(
    # # #         long_name = 'river basins',
    # # #         _FillValue=flw._core._mv
    # # #     )            
    # # # )

    # # # # Upstream area map
    # # # # #################
    # # # da_uparea = xr.DataArray(
    # # #     name='uparea',
    # # #     data=uparea_lr,
    # # #     coords=coords,
    # # #     dims=dims,
    # # #     attrs=dict(
    # # #         long_name = 'upstream drainage area', 
    # # #         _FillValue=flw._core._mv
    # # #     )            
    # # # )

    # # # # Annual mean discharge
    # # # # #####################
    # # # da_runoff = xr.DataArray(
    # # #     name='Q',
    # # #     data=Qmean_lr,
    # # #     coords=coords,
    # # #     dims=dims,
    # # #     attrs=dict(
    # # #         long_name = 'annual total discharge', 
    # # #         _FillValue=F8_FILLVAL
    # # #     )
    # # # )

    # # # # River length
    # # # # ############
    # # # da_rivlen = xr.DataArray(
    # # #     name='rivlen',
    # # #     data=rivlen,
    # # #     coords=coords,
    # # #     dims=dims,
    # # #     attrs=dict(
    # # #         long_name = 'main river channel length', 
    # # #         _FillValue=F8_FILLVAL
    # # #     )            
    # # # )

    # # # # River slope
    # # # # ###########
    # # # da_rivslp = xr.DataArray(
    # # #     name='rivslp',
    # # #     data=rivslp,
    # # #     # data=rivslp2,
    # # #     coords=coords,
    # # #     dims=dims,
    # # #     attrs=dict(
    # # #         long_name = 'main river channel slope', 
    # # #         _FillValue=F8_FILLVAL
    # # #     )            
    # # # )

    # # # # Channel bankfull width
    # # # # ######################
    # # # da_wbmain = xr.DataArray(
    # # #     name='wbmain',
    # # #     data=wbmain,
    # # #     coords=coords,
    # # #     dims=dims,
    # # #     attrs=dict(
    # # #         long_name = 'main channel bankfull width',
    # # #         _FillValue=F8_FILLVAL
    # # #     )
    # # # )

    # # # # Channel bankfull depth
    # # # # ######################
    # # # da_hbmain = xr.DataArray(
    # # #     name='hbmain',
    # # #     data=hbmain,
    # # #     coords=coords,
    # # #     dims=dims,
    # # #     attrs=dict(
    # # #         long_name = 'main channel bankfull depth',
    # # #         _FillValue=F8_FILLVAL
    # # #     )
    # # # )

    # # # # Tributary bankfull width
    # # # # ########################
    # # # da_wbtrib = xr.DataArray(
    # # #     name='wbtrib',
    # # #     data=wbtrib,
    # # #     coords=coords,
    # # #     dims=dims,
    # # #     attrs=dict(
    # # #         long_name = 'tributary channel bankfull width',
    # # #         _FillValue=F8_FILLVAL
    # # #     )
    # # # )

    # # # # Tributary bankfull depth
    # # # # ########################
    # # # da_hbtrib = xr.DataArray(
    # # #     name='hbtrib',
    # # #     data=hbtrib,
    # # #     coords=coords,
    # # #     dims=dims,
    # # #     attrs=dict(
    # # #         long_name = 'tributary channel bankfull depth',
    # # #         _FillValue=F8_FILLVAL
    # # #     )
    # # # )

    # # # # Channel Manning's n
    # # # # ###################
    # # # da_n_channel = xr.DataArray(
    # # #     name='n_channel',
    # # #     data=n_channel,
    # # #     coords=coords,
    # # #     dims=dims,
    # # #     attrs=dict(
    # # #         long_name = 'Channel Manning n',
    # # #         _FillValue=F8_FILLVAL
    # # #     )
    # # # )

    # # # # Overland Manning's n
    # # # # ####################
    # # # da_n_overland = xr.DataArray(
    # # #     name='n_overland',
    # # #     data=n_overland,
    # # #     coords=coords,
    # # #     dims=dims,
    # # #     attrs=dict(
    # # #         long_name = 'Overland Manning n',
    # # #         _FillValue=F8_FILLVAL
    # # #     )
    # # # )
        
    # # # # Write output
    # # # da_flwdir.rio.to_raster(os.path.join(AUX_DATADIR, "merit_vilcanota_Fdir.tif"))
    # # # da_ids.rio.to_raster(os.path.join(AUX_DATADIR, "merit_vilcanota_ID.tif"))
    # # # da_dnids.rio.to_raster(os.path.join(AUX_DATADIR, "merit_vilcanota_dnID.tif"))
    # # # da_basins.rio.to_raster(os.path.join(AUX_DATADIR, "merit_vilcanota_basins.tif"))
    # # # da_uparea.rio.to_raster(os.path.join(AUX_DATADIR, "merit_vilcanota_Atotal.tif"))
    # # # da_runoff.rio.to_raster(os.path.join(AUX_DATADIR, "merit_vilcanota_Qanntot.tif"))
    # # # da_rivlen.rio.to_raster(os.path.join(AUX_DATADIR, "merit_vilcanota_Fdis.tif"))
    # # # da_rivslp.rio.to_raster(os.path.join(AUX_DATADIR, "merit_vilcanota_Schannel.tif"))
    # # # da_wbmain.rio.to_raster(os.path.join(AUX_DATADIR, "merit_vilcanota_Wmain.tif"))
    # # # da_hbmain.rio.to_raster(os.path.join(AUX_DATADIR, "merit_vilcanota_Hmain.tif"))
    # # # da_wbtrib.rio.to_raster(os.path.join(AUX_DATADIR, "merit_vilcanota_Wtrib.tif"))
    # # # da_hbtrib.rio.to_raster(os.path.join(AUX_DATADIR, "merit_vilcanota_Htrib.tif"))
    # # # da_n_channel.rio.to_raster(os.path.join(AUX_DATADIR, "merit_vilcanota_n_channel.tif"))
    # # # da_n_overland.rio.to_raster(os.path.join(AUX_DATADIR, "merit_vilcanota_n_overland.tif"))

if __name__ == "__main__":
    main()


# Not used:

# np.random.seed(seed=101)
# matplotlib.rcParams['savefig.bbox'] = 'tight'
# matplotlib.rcParams['savefig.dpi'] = 256
# plt.style.use('seaborn-whitegrid')

# def quickplot(gdfs=[], maps=[], hillshade=True, title='', filename='flw', save=False):
#     fig = plt.figure(figsize=(8, 15))
#     ax = fig.add_subplot(projection=ccrs.PlateCarree())
#     # plot hillshade background
#     if hillshade:
#         ls = matplotlib.colors.LightSource(azdeg=115, altdeg=45)
#         hillshade = ls.hillshade(
#             np.ma.masked_equal(elevtn, -9999), vert_exag=1e3)
#         ax.imshow(hillshade, origin='upper', extent=flw.extent,
#                   cmap='Greys', alpha=0.3, zorder=0)
#     # plot geopandas GeoDataFrame,
#     for gdf, kwargs in gdfs:
#         gdf.plot(ax=ax, **kwargs)
#     for data, nodata, kwargs in maps:
#         ax.imshow(np.ma.masked_equal(data, nodata),
#                   origin='upper', extent=flw.extent, **kwargs)
#     ax.set_aspect('equal')
#     ax.set_title(title, fontsize='large')
#     ax.text(0.01, 0.01, 'created with pyflwdir',
#             transform=ax.transAxes, fontsize='large')
#     if save:
#         plt.savefig(f'{filename}.png')
#     return ax
