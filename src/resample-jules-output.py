#!/usr/bin/env python

import os
import sys
import tempfile
import subprocess
import numpy as np
import xarray

START = int(os.environ['JULES_START_YEAR'])
END = int(os.environ['JULES_END_YEAR'])
YEARS = np.arange(START, END + 1)
SUITE = str(os.environ['JULES_SUITE'])
ID_STEM = str(os.environ['JULES_ID_STEM'])
JOB_NAME = str(os.environ['JULES_JOB_NAME'])
PROFILE_NAME = str(os.environ['JULES_PROFILE_NAME'])
YMAX = float(os.environ['N'])
YMIN = float(os.environ['S'])
XMAX = float(os.environ['E'])
XMIN = float(os.environ['W'])
RES = str(os.environ['RES'])
INPUT_FILE_SUFFIX = str(os.environ['INPUT_FILE_SUFFIX'])
OUTPUT_FILE_SUFFIX = str(os.environ['OUTPUT_FILE_SUFFIX'])
INPUT_DIRECTORY = str(os.environ['INPUT_DIRECTORY'])
OUTPUT_DIRECTORY = str(os.environ['OUTPUT_DIRECTORY'])

if RES == '30sec':
    DRES = 0.008333333333
elif RES == '05min':
    DRES = 0.083333333333
elif RES == '15min':
    DRES = 0.25
else:
    stop()
    
def main():
    job_name = JOB_NAME.format(year=START)
    FN = os.path.join(
        INPUT_DIRECTORY, ID_STEM + '.' + job_name + '.' + PROFILE_NAME + '.' + str(START) + '.' + INPUT_FILE_SUFFIX + '.nc'
    )    
    with xarray.open_dataset(FN) as x:
        lat = x['lat'].values
        ns_lat = lat[0] > lat[1]
    # write CDO gridfile
    xsize = abs((XMAX - XMIN) / DRES)
    ysize = abs((YMAX - YMIN) / DRES)
    xfirst = XMIN + DRES / 2.
    xinc = abs(DRES)
    if ns_lat:
        yfirst = YMAX - DRES / 2.
        yinc = abs(DRES) * -1.
    else:
        yfirst = YMIN + DRES / 2.
        yinc = abs(DRES)
    with open('/tmp/gridfile.txt', 'w') as f:
        f.write('gridtype=lonlat\n')
        f.write('xsize=' + str(int(xsize)) + '\n')
        f.write('ysize=' + str(int(ysize)) + '\n')
        f.write('xfirst=' + str(xfirst) + '\n')
        f.write('xinc=' + str(xinc) + '\n')
        f.write('yfirst=' + str(yfirst) + '\n')
        f.write('yinc=' + str(yinc) + '\n')
        f.close()

    for yr in YEARS:
        job_name = JOB_NAME.format(year=yr)
        IN_FN = os.path.join(
            INPUT_DIRECTORY, ID_STEM + '.' + job_name + '.' + PROFILE_NAME + '.' + str(yr) + '.' + INPUT_FILE_SUFFIX + '.nc'            
        )
        TMP_FN = tempfile.NamedTemporaryFile(suffix='.nc')
        OUT_FN = os.path.join(
            OUTPUT_DIRECTORY, ID_STEM + '.' + job_name + '.' + PROFILE_NAME + '.' + str(yr) + '.' + OUTPUT_FILE_SUFFIX + '.nc'
        )
        # use bilinear interpolation for all continuous variables
        subprocess.run([
            'cdo',
            'remapbil,/tmp/gridfile.txt',
            IN_FN,
            TMP_FN.name
        ])
        # ensure that variables have datatype 'double', not 'short',
        # which seems to cause problems in JULES (not exactly sure why...)
        subprocess.run([
            'cdo',
            '-b', 'F64', 'copy',
            TMP_FN.name,
            OUT_FN
        ])
        
        
if __name__ == '__main__':
    main()

# echo "gridtype=$GRIDTYPE
# xsize=$XSIZE
# ysize=$YSIZE
# xfirst=$XFIRST
# xinc=$XINC
# yfirst=$YFIRST
# yinc=$YINC" > /tmp/gridfile.txt

# for YEAR in $(seq ${JULES_START_YEAR} ${JULES_END_YEAR})
# do
#     FN=$INPUT_DIRECTORY/$JULES_ID_STEM.S2.$JULES_PROFILE_NAME.$YEAR.2D.nc
#     # # 1-select variables
#     # OUTFN1=tmp1.nc
#     # cdo select,name=lat,lon,runoff,surf_roff,sub_surf_roff $FN $OUTFN1
#     # # 2-resample
#     OUTFN1=tmp1.nc
#     # cdo remapbil,/tmp/gridfile.txt $OUTFN1 $OUTFN2
#     cdo remapbil,/tmp/gridfile.txt $FN $OUTFN1
#     # 3-convert datatype
#     OUTFN2=$OUTPUT_DIRECTORY/$JULES_ID_STEM.S2.$JULES_PROFILE_NAME.$YEAR.2D.regrid.nc
#     cdo -b F64 copy $OUTFN1 $OUTFN2
#     # tidy up
#     rm -f tmp1.nc tmp2.nc
# done

    
