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

# export JULES_START_YEAR=$(yaml $CONFIG_FILE "['jules']['start_year']")
# export JULES_END_YEAR=$(yaml $CONFIG_FILE "['jules']['end_year']")
# export JULES_SUITE=$(yaml $CONFIG_FILE "['jules']['suite']")
# export JULES_ID_STEM=$(yaml $CONFIG_FILE "['jules']['id_stem']")
# export JULES_PROFILE_NAME=$(yaml $CONFIG_FILE "['jules']['profile_name']")
export MERIT_DATADIR=$(yaml $CONFIG_FILE "['mosart']['merit_hydro_directory']")
export MERIT_IHU_DATADIR=$(yaml $CONFIG_FILE "['mosart']['merit_ihu_directory']")
export GEOMORPHO_DATADIR=$(yaml $CONFIG_FILE "['mosart']['geomorpho90m_directory']")
export AUX_DATADIR=$(yaml $CONFIG_FILE "['mosart']['aux_directory']")
# export MEAN_ANNUAL_RUNOFF=$(yaml $CONFIG_FILE "['mosart']['mean_annual_runoff']")
NORTH=$(yaml $CONFIG_FILE "['mosart']['ymax']")
SOUTH=$(yaml $CONFIG_FILE "['mosart']['ymin']")
EAST=$(yaml $CONFIG_FILE "['mosart']['xmax']")
WEST=$(yaml $CONFIG_FILE "['mosart']['xmin']")
export N=$(ceil $NORTH)
export S=$(floor $SOUTH)
export E=$(ceil $EAST)
export W=$(floor $WEST)
# Must be one of 05min 15min 30sec [TODO: check]
export RES=$(yaml $CONFIG_FILE "['mosart']['resolution']")
# Identify river outlets using the 30sec map
export OUTLET_X=$(yaml $CONFIG_FILE "['mosart']['outlet_x']")
export OUTLET_Y=$(yaml $CONFIG_FILE "['mosart']['outlet_y']")
export LATLON_LOCATION=$(yaml $CONFIG_FILE "['mosart']['grass_location']")
export LATLON_MAPSET=${LATLON_LOCATION}/$(yaml $CONFIG_FILE "['mosart']['grass_mapset']")
# export SCALE_FACTOR=$(yaml $CONFIG_FILE "['mosart']['merit_scale_factor']")
# export OUTPUT_DIRECTORY=$(yaml $CONFIG_FILE "['mosart']['output_directory']")
# export INPUT_DIRECTORY=$(yaml $CONFIG_FILE "['mosart']['input_directory']")
export SRC_DIR=$(pwd)/../src

# =================================================================== #
# Enable Anaconda to be activated/deactivated within the script
# =================================================================== #

CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh

# Deactivate
conda deactivate

# =================================================================== #
# Crop MERIT-IHU data to the current region
# =================================================================== #

# Set GDAL_DATA here
# (On my PC at least, this has been causing issues - probably Anaconda-related)
export GDAL_DATA=$(/usr/bin/gdal-config --datadir)

if [[ ! -d ${AUX_DATADIR} ]]
then
    mkdir -p ${AUX_DATADIR}
fi

# for VAR in elevtn flwdir outlat outlon rivlen rivslp rivwth uparea
for VAR in elevtn flwdir outlat outlon rivlen rivslp uparea
do
    for RES in 30sec 05min 15min
    do	
	gdalwarp \
	    -overwrite \
	    -te $W $S $E $N \
	    ${MERIT_IHU_DATADIR}/${RES}_${VAR}.tif \
	    ${AUX_DATADIR}/${RES}_${VAR}_subrgn.tif 
    done
done

# =================================================================== #
# Mosaic Geomorpho90m [slope] data
# =================================================================== #

# It would be possible to calculate slope in GIS, but this is a bit
# tricky when dealing with lat/long format because the horizontal units
# (i.e. degrees) are different to the vertical units (m). As a result
# either a correction factor needs to be applied, or the map needs to
# be converted to a different projection. Both are difficult to
# generalise. Instead we use the precomputed maps from geomorpho90m.

VAR=slope
ptn="${VAR}_90M_*.tif"
find $GEOMORPHO_DATADIR -type f -iname "${ptn}" > /tmp/geomorpho90m_${VAR}_filenames.txt
gdalbuildvrt \
    -overwrite \
    -te $W $S $E $N \
    -tr 0.0008333333333 0.0008333333333 \
    -input_file_list /tmp/geomorpho90m_${VAR}_filenames.txt \
    ${AUX_DATADIR}/geomorpho90m_${VAR}.vrt

# =================================================================== #
# Run Python script to extract basins
# =================================================================== #

# NB We use `pyflwdir` to be consistent with the upscaling routines.

conda activate mosart
python ${SRC_DIR}/make-river-basins.py
conda deactivate 

# # =================================================================== #
# # Compute average slope in the coarse grid cell
# # =================================================================== #

# export GDAL_DATA=$(/usr/bin/gdal-config --datadir)

# # Create location/mapset if they do not already exist
# if [[ ! -d ${LATLON_LOCATION} ]]
# then
#     grass -c -e epsg:4326 ${LATLON_LOCATION}
# fi

# if [[ ! -d ${LATLON_MAPSET} ]]
# then
#     grass -c -e ${LATLON_MAPSET}
# fi

# chmod u+x ${SRC_DIR}/grass_process_merit.sh
# export GRASS_BATCH_JOB=${SRC_DIR}/grass_process_merit.sh
# grass76 ${LATLON_MAPSET}
# unset GRASS_BATCH_JOB

# NOT USED:

# # # =================================================================== #
# # # Mosaic MERIT Hydro maps for current study area
# # # =================================================================== #

# # if [[ ! -d ${AUX_DATADIR} ]]
# # then
# #     mkdir -p ${AUX_DATADIR}
# # fi

# # if [[ ! -d ${OUTPUT_DIRECTORY} ]]
# # then
# #     mkdir -p ${OUTPUT_DIRECTORY}
# # fi

# # for VAR in elv dir wth upa
# # do
# #     ptn="*_${VAR}.tif"
# #     find ${MERIT_DATADIR} -type f -iname "${ptn}" > /tmp/merit_${VAR}_filenames.txt
# #     # if [[ ! -f ${AUX_DATADIR}/merit_${VAR}.tif ]]
# #     # then	
# #     gdalbuildvrt \
# # 	-overwrite \
# # 	-te $W $S $E $N \
# # 	-tr 0.0008333333333 0.0008333333333 \
# # 	-input_file_list /tmp/merit_${VAR}_filenames.txt \
# # 	${AUX_DATADIR}/merit_${VAR}.vrt
# #     gdal_translate ${AUX_DATADIR}/merit_${VAR}.vrt ${AUX_DATADIR}/merit_${VAR}.tif
# #     # fi    
# # done

# # =================================================================== #
# # Run GRASS GIS scripts to process MERIT Hydro data
# # =================================================================== #

# # Create location/mapset if they do not already exist
# if [[ ! -d ${LATLON_LOCATION} ]]
# then
#     grass -c -e epsg:4326 ${LATLON_LOCATION}
# fi

# if [[ ! -d ${LATLON_MAPSET} ]]
# then
#     grass -c -e ${LATLON_MAPSET}
# fi


# chmod u+x ${SRC_DIR}/grass_process_merit.sh
# export GRASS_BATCH_JOB=${SRC_DIR}/grass_process_merit.sh
# grass76 ${LATLON_MAPSET}
# unset GRASS_BATCH_JOB

# # =================================================================== #
# # Run upscaling routines; write output data
# # =================================================================== #

# conda activate mosart
# python ${SRC_DIR}/upscale-routing-params.py
# conda deactivate 

# # =================================================================== #
# # Mosaic Geomorpho90m [slope] data
# # =================================================================== #

# # It would be possible to calculate slope in GIS, but this is a bit
# # tricky when dealing with lat/long format because the horizontal units
# # (i.e. degrees) are different to the vertical units (m). As a result
# # either a correction factor needs to be applied, or the map needs to
# # be converted to a different projection. Both are difficult to
# # generalise. Instead we use the precomputed maps from geomorpho90m.

# VAR=slope
# ptn="${VAR}_90M_*.tif"
# find $GEOMORPHO_DATADIR -type f -iname "${ptn}" > /tmp/geomorpho90m_${VAR}_filenames.txt
# # if [[ ! -f ${AUX_DATADIR}/geomorpho90m_${VAR}.tif ]]
# # then
# gdalbuildvrt \
#     -overwrite \
#     -te $W $S $E $N \
#     -tr 0.0008333333333 0.0008333333333 \
#     -input_file_list /tmp/geomorpho90m_${VAR}_filenames.txt \
#     ${AUX_DATADIR}/geomorpho90m_${VAR}.vrt
# gdal_translate ${AUX_DATADIR}/geomorpho90m_${VAR}.vrt ${AUX_DATADIR}/geomorpho90m_${VAR}.tif
# # fi

# # =================================================================== #
# # Run Python script to extract basins
# # =================================================================== #

# # NB We use `pyflwdir` to be consistent with the upscaling routines.

# conda activate mosart
# python ${SRC_DIR}/make-river-basins.py
# conda deactivate 










# NOT USED:

# Help()
# {
#     echo "Program to create input maps for mosartwmpy."
#     echo
#     echo "Syntax: create-mosartwm-input.sh [-h|o]"
#     echo "options:"
#     echo "-h | --help          Print this help message."
#     echo "-o | --overwrite     Overwrite existing database."
#     echo "--merit-datadir      Location of MERIT Hydro data."
#     echo "--geomorpho90m-datadir Location of Geomorpho90m data."
#     echo "--aux-datadir        Location to store intermediate outputs."
#     echo "--res                Resolution of model region, in DMS (e.g. 0:15, not 0.25)."
#     echo "--ext                Extent of model region (xmin, ymin, xmax, ymax), in DMS."
#     echo "-d | --destdir       Output directory."
#     echo
# }

# while [[ $# -gt 0 ]]
# do
#     key="$1"
#     case $key in
# 	-h|--help)
# 	    Help
# 	    # HELP="$2"
# 	    # shift
# 	    # shift
# 	    exit
# 	    ;;
# 	-o|--overwrite)
# 	    OVERWRITE='--overwrite'
# 	    shift
# 	    ;;
# 	--merit-datadir)
# 	    MERIT_DATADIR="$2"
# 	    shift
# 	    shift
# 	    ;;
# 	--geomorpho90m-datadir)
# 	    GEOMORPHO90M_DATADIR="$2"
# 	    shift
# 	    shift
# 	    ;;
#       --aux-datadir)
#           AUX_DATADIR="$2"
#           shift
#           shift
#           ;;
#       --res)
#           XRES="$2"
#           YRES="$3"
#           shift
#           shift
#           ;;
#       --ext)
#           XMIN="$2"
#           YMIN="$3"
#           XMAX="$4"
#           YMAX="$5"
#           shift
#           shift
#           ;;
# 	-d|--destdir)
# 	    OUTDIR="$2"
# 	    shift
# 	    shift
# 	    ;;
# 	*)  # unknown option
# 	    POSITIONAL+=("$1") # save it in an array for later
# 	    shift # past argument
# 	    ;;
#     esac
# done
