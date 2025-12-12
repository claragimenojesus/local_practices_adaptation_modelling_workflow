#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import pyflwdir
import rasterio
from rasterio import features
from rasterio.transform import Affine
import numpy as np
import netCDF4
import xarray as xr
import rioxarray

# local libraries
from pyflwdir.gis_utils import affine_to_coords
# from gislib import rio

# Default fill vals for netCDF 
F8_FILLVAL = netCDF4.default_fillvals['f8']
F4_FILLVAL = netCDF4.default_fillvals['f4']
I4_FILLVAL = netCDF4.default_fillvals['i4']

AUX_DATADIR = str(os.environ['AUX_DATADIR'])
RES = str(os.environ['RES'])
OUTLET_X = float(os.environ['OUTLET_X'])
OUTLET_Y = float(os.environ['OUTLET_Y'])
OUTLET_COORDS = (OUTLET_X, OUTLET_Y)

def main():

    for res in ['30sec', '05min', '15min']:
        dir_file = os.path.join(AUX_DATADIR, res + '_flwdir_subrgn.tif')
        
        # Open flow direction map
        with rasterio.open(dir_file, 'r') as src:
            flwdir = src.read(1)
            transform = src.transform
            crs = src.crs
            latlon = crs.to_epsg() == 4326
        
        # Create a flow direction object
        flw = pyflwdir.from_array(
            flwdir,
            ftype='d8',
            transform=transform,
            latlon=latlon,
            cache=True
        )
        # Compute river basins
        basins = flw.basins(xy=OUTLET_COORDS).astype(np.int32)
        # Get coordinates
        xs, ys = affine_to_coords(transform, basins.shape)
        coords={'y': ys, 'x': xs}
        dims=coords.keys()
        # Create xarray.DataArray object, write raster
        da_basins = xr.DataArray(
            name='basins',
            data=basins,
            coords=coords,
            dims=dims,
            attrs=dict(
                long_name='basins', 
                _FillValue=I4_FILLVAL
            )            
        )
        da_basins = da_basins.rio.set_crs(crs)
        da_basins.rio.to_raster(
            os.path.join(AUX_DATADIR, res + "_ihu_basins.tif")
        )

if __name__ == "__main__":
    main()
