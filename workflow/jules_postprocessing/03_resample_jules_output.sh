#!/bin/bash

yaml() {
    python3 -c "import yaml;print(yaml.safe_load(open('$1'))$2)"
}

ceil() {
    python3 -c "import math;print(math.ceil($1*4)/4)" 
}

floor() {
    python3 -c "import math;print(math.floor($1*4)/4)"
}

CONFIG_FILE="$1"

export JULES_START_YEAR=$(yaml $CONFIG_FILE "['jules']['start_year']")
export JULES_END_YEAR=$(yaml $CONFIG_FILE "['jules']['end_year']")
export JULES_SUITE=$(yaml $CONFIG_FILE "['jules']['suite']")
export JULES_ID_STEM=$(yaml $CONFIG_FILE "['jules']['id_stem']")
export JULES_JOB_NAME=$(yaml $CONFIG_FILE "['jules']['job_name']")
export JULES_PROFILE_NAME=$(yaml $CONFIG_FILE "['jules']['profile_name']")
# export JULES_GRIDFILE=$(yaml $CONFIG_FILE "['jules']['gridfile']")
# export JULES_OUTPUT_DIRECTORY=$(yaml $CONFIG_FILE "['jules']['raw_output_directory']")
# export GRIDTYPE=$(yaml $CONFIG_FILE "['resample_jules_output']['gridtype']")
# export XSIZE=$(yaml $CONFIG_FILE "['resample_jules_output']['xsize']")
# export YSIZE=$(yaml $CONFIG_FILE "['resample_jules_output']['ysize']")
# export XFIRST=$(yaml $CONFIG_FILE "['resample_jules_output']['xfirst']")
# export YFIRST=$(yaml $CONFIG_FILE "['resample_jules_output']['yfirst']")
# export XINC=$(yaml $CONFIG_FILE "['resample_jules_output']['xinc']")
# export YINC=$(yaml $CONFIG_FILE "['resample_jules_output']['yinc']")
NORTH=$(yaml $CONFIG_FILE "['resample_jules_output']['ymax']")
SOUTH=$(yaml $CONFIG_FILE "['resample_jules_output']['ymin']")
EAST=$(yaml $CONFIG_FILE "['resample_jules_output']['xmax']")
WEST=$(yaml $CONFIG_FILE "['resample_jules_output']['xmin']")
export N=$(ceil $NORTH)
export S=$(floor $SOUTH)
export E=$(ceil $EAST)
export W=$(floor $WEST)
# Must be one of 05min 15min 30sec [TODO: check]
export RES=$(yaml $CONFIG_FILE "['resample_jules_output']['resolution']")
export INPUT_FILE_SUFFIX=$(yaml $CONFIG_FILE "['resample_jules_output']['input_file_suffix']")
export OUTPUT_FILE_SUFFIX=$(yaml $CONFIG_FILE "['resample_jules_output']['output_file_suffix']")
export INPUT_DIRECTORY=$(yaml $CONFIG_FILE "['resample_jules_output']['input_directory']")
export OUTPUT_DIRECTORY=$(yaml $CONFIG_FILE "['resample_jules_output']['output_directory']")
export SRC_DIR=$(pwd)/../src

CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh

# Deactivate
conda deactivate

conda activate mosart
python ${SRC_DIR}/resample-jules-output.py
conda deactivate 

# NOT USED:

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

# # We also need a grid area file, so make that here
# cdo gridarea $OUTPUT_DIRECTORY/JULES_vn6.1.S2.daily_hydrology.1980.2D.regrid.nc $OUTPUT_DIRECTORY/gridarea.nc

# CONDA_BASE=$(conda info --base)
# source $CONDA_BASE/etc/profile.d/conda.sh
# conda deactivate

# # # export DATADIR=$(pwd)/../data
# # export DATADIR=/home/sm510/projects/ganges_water_machine/jules-output/u-cg201/netcdf
# # export AUX_DATADIR=$(pwd)/../data/aux
# # export OUTDIR=$HOME/projects/rahu/mosart-input

# export SRC_DIR=$(pwd)/../src
# conda activate mosart
# # python ${SRC_DIR}/format-jules-output.py
# bash ${SRC_DIR}/regrid-runoff-data.sh
# # python ${SRC_DIR}/aggregate-jules.py
# conda deactivate 
