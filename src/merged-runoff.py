#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import numpy as np
# import xarray as xr
import subprocess

START = int(os.environ['JULES_START_YEAR'])
END = int(os.environ['JULES_END_YEAR'])
YEARS = np.arange(START, END + 1)
SUITE = str(os.environ['JULES_SUITE'])
ID_STEM = str(os.environ['JULES_ID_STEM'])
JOB_NAME = str(os.environ['JULES_JOB_NAME'])
PROFILE_NAME = str(os.environ['JULES_PROFILE_NAME'])
INPUT_DIRECTORY = str(os.environ['INPUT_DIRECTORY'])
INPUT_FILE_SUFFIX = str(os.environ['INPUT_FILE_SUFFIX'])
OUTPUT_FILENAME = str(os.environ['OUTPUT_FILENAME'])
# OUTPUT_DIRECTORY = str(os.environ['OUTPUT_DIRECTORY'])

def main():    
    fs = []
    for yr in YEARS:
        job_name = JOB_NAME.format(year=yr)
        FN = ID_STEM + '.' + job_name + '.' + PROFILE_NAME + '.' + str(yr) + '.' + INPUT_FILE_SUFFIX + '.nc'
        NEWFN = os.path.join('/tmp', 'jules_runoff_' + str(yr) + '.nc')
        subprocess.run(['cdo', '-select,name=runoff,surf_roff,sub_surf_roff', os.path.join(INPUT_DIRECTORY, FN), NEWFN])
        # fs.append(os.path.join(INPUT_DIRECTORY, FN))
        fs.append(NEWFN)

    s = ''
    for f in fs: s = s + f + ' '
    print(s)
    print(OUTPUT_FILENAME)
    subprocess.run(['cdo', 'mergetime', s, OUTPUT_FILENAME])
    
        
    # x = xr.open_mfdataset(fs)
    # # # write this dataset as a single file - FAR TOO SLOW!!!
    # # x_runoff = x[['surf_roff', 'sub_surf_roff', 'runoff']]
    # # x_runoff.to_netcdf(os.path.join(OUTPUT_DIRECTORY, 'jules_runoff_' + str(START) + '_' + str(END) + '.nc'))
    # # compute annual averages for hydrological variables
    # jules_vars = [
    #     'surf_roff', 'sub_surf_roff', 'runoff',
    #     'precip', 'ecan_gb', 'elake', 'esoil_gb',
    #     'fao_et0'
    # ]        
    # for var in jules_vars:
    #     x[var] = x[var] * 60 * 60 * 24 / 1000  # kg m-2 s-1 -> m d-1
    # # m d-1 -> m y-1
    # x_year = x.groupby("time.year").sum(dim="time")
    # x_annual_mean = x_year.mean(dim="year")
    # x_annual_mean.to_netcdf(OUTPUT_FILENAME)
    # x.close()

if __name__ == '__main__':
    main()
