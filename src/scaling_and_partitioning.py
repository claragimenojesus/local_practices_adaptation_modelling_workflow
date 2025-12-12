#!/usr/bin/env python
# -*- coding: utf-8 -*-

import xarray as xr
import numpy as np
import os

INPUT_DIRECTORY = str(os.environ['INPUT_DIRECTORY'])
INPUT_FILENAME = str(os.environ['INPUT_FILENAME'])
# changing this to incorporate NBS
#MOSART_MASK = str(os.environ['MOSART_MASK'])
LAND_FRAC = str(os.environ['LAND_FRAC'])
GRID_AREA = str(os.environ['GRID_AREA'])
#MOSART_AREA = str(os.environ['MOSART_AREA']) # changes for NBS
HYD_COND = str(os.environ['HYD_COND'])
#OUTPUT_DIRECTORY = str(os.environ['OUTPUT_DIRECTORY']) # changes for NBS

# adapting script to add NBS scenarios - 22/10/2024 CGJ
NBS = str(os.environ['NBS'])
MOSART_MASK = os.path.join(INPUT_DIRECTORY, NBS, 'land_new.nc')
MOSART_AREA = os.path.join(INPUT_DIRECTORY, NBS, 'mosart_new.nc')
OUTPUT_DIRECTORY = os.path.join(INPUT_DIRECTORY, NBS, 'coupled_scpar_runoff.nc')
## end of changes - 22/10/2024

def main():
    IN_FN = os.path.join(INPUT_DIRECTORY, NBS, INPUT_FILENAME)
    print(IN_FN)
    runoff = xr.open_dataset(IN_FN)
    mosart_mask = xr.open_dataset(MOSART_MASK)['mask']
    mosart_mask['lon'] = runoff['lon']
    mosart_mask['lat'] = runoff['lat']
    runoff_input = np.multiply(runoff, mosart_mask)

    runoff.close()
    mosart_mask.close()

    land_frac = xr.open_dataset(LAND_FRAC)['land_frac']
    grid_area_102 = xr.open_dataset(GRID_AREA)['area']
    area_hydrosheds_total = np.multiply(grid_area_102[0,:,:], land_frac).sum(dim=("lat","lon"))
    mosart_area = xr.open_dataset(MOSART_AREA)['area']
    mosart_area['lat'] = runoff_input['lat']
    mosart_area['lon'] = runoff_input['lon']
    area_mosart_total = (np.multiply(runoff_input['runoff'].groupby("time.year").sum(dim="time").mean(dim="year"), mosart_area)/(runoff_input['runoff'].groupby("time.year").sum(dim="time").mean(dim="year"))).sum(dim=("lat","lon"))

    land_frac.close()
    grid_area_102.close()
    mosart_area.close()

    # Calculate the scaling parameter due to remapping method (it should be the same for all)
    scaling = area_hydrosheds_total/area_mosart_total

    # Apply scaling parameter
    new_scaled_regridded = runoff_input
    new_scaled_regridded['runoff'] = new_scaled_regridded['runoff']*scaling
    new_scaled_regridded['surf_roff'] = new_scaled_regridded['surf_roff']*scaling
    new_scaled_regridded['sub_surf_roff'] = new_scaled_regridded['sub_surf_roff']*scaling
    new_scaled_regridded['ei_gb'] = new_scaled_regridded['ei_gb']*scaling
    new_scaled_regridded['elake'] = new_scaled_regridded['elake']*scaling
    new_scaled_regridded['esoil_gb'] = new_scaled_regridded['esoil_gb']*scaling
    new_scaled_regridded['ecan_gb'] = new_scaled_regridded['ecan_gb']*scaling
    new_scaled_regridded['melt'] = new_scaled_regridded['melt']*scaling
    new_scaled_regridded['area'] = new_scaled_regridded['area']*scaling
    new_scaled_regridded['rain'] = new_scaled_regridded['rain']*scaling
    new_scaled_regridded['perc'] = new_scaled_regridded['perc']*scaling
    new_scaled_regridded['fao_et0'] = new_scaled_regridded['fao_et0']*scaling
    new_scaled_regridded['precip'] = new_scaled_regridded['precip']*scaling
    new_scaled_regridded['rainfall'] = new_scaled_regridded['rainfall']*scaling
    new_scaled_regridded['snowfall'] = new_scaled_regridded['snowfall']*scaling

    runoff_input.close()
    area_hydrosheds_total.close()
    area_mosart_total.close()

    # Load rock hydraulic conductivity
    hyd_cond = xr.open_dataset(HYD_COND).fillna(0) # this is in m/s. jules output is in kg/m2/s=mm/s
    # for some reasons rasterization produces NA in the lower limit

    ## I am copying the hydraulic conductivity values across time dimension. This allows us to do array operations in xarray without dimension problems
    hyd_cond0 = xr.zeros_like(new_scaled_regridded.sub_surf_roff)

    for t in range(new_scaled_regridded.sizes["time"]):
        hyd_cond0[t,...]=hyd_cond.Band1.values
    ### Create deep reservoir

    q_deep = xr.zeros_like(new_scaled_regridded.sub_surf_roff)
    # calculate deep percolation. Hydraulic conductivity is multiplied by 1000 (m/s to mm/s )
    q_deep[...]=np.minimum(new_scaled_regridded.sub_surf_roff,1000*hyd_cond0)

    ## correct shallow subsurface flow
    new_scaled_regridded.sub_surf_roff.values = np.add(new_scaled_regridded.sub_surf_roff,-q_deep)

    q_deep = q_deep.rename('q_deep')
    q_deep.attrs['standard_name'] = "Gridbox deep sub-surface runoff"
    new_scaled_regridded['q_deep'] = q_deep
    new_scaled_regridded.to_netcdf(OUTPUT_DIRECTORY)

    new_scaled_regridded.close()
    hyd_cond.close()


if __name__ == '__main__':
    main()
