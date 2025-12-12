#!/bin/bash

yaml() {
    python3 -c "import yaml;print(yaml.safe_load(open('$1'))$2)"
}

CONFIG_FILE="$1"

export AUX_DATADIR=$(yaml $CONFIG_FILE "['mosart']['aux_directory']")
export MEAN_ANNUAL_RUNOFF=$(yaml $CONFIG_FILE "['mosart']['mean_annual_runoff']")
export MEAN_ANNUAL_RUNOFF_VARNAME=$(yaml $CONFIG_FILE "['mosart']['mean_annual_runoff_varname']")
# Must be one of 05min 15min 30sec [TODO: check]
export RES=$(yaml $CONFIG_FILE "['mosart']['resolution']")
# Identify river outlets using the 30sec map
export OUTLET_X=$(yaml $CONFIG_FILE "['mosart']['outlet_x']")
export OUTLET_Y=$(yaml $CONFIG_FILE "['mosart']['outlet_y']")
export OUTPUT_DIRECTORY=$(yaml $CONFIG_FILE "['mosart']['output_directory']")
export SRC_DIR=$(pwd)/../src

# Set GDAL_DATA here
# (On my PC at least, this has been causing issues - probably Anaconda-related)
export GDAL_DATA=$(gdal-config --datadir)

# =================================================================== #
# Enable Anaconda to be activated/deactivated within the script
# =================================================================== #

CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh

# Deactivate
conda deactivate

# =================================================================== #
# Run upscaling routines; write output data
# =================================================================== #

conda activate mosart
python ${SRC_DIR}/make-mosart-input.py
conda deactivate 

# # =================================================================== #
# # Mosaic MERIT Hydro maps for current study area
# # =================================================================== #

# if [[ ! -d ${AUX_DATADIR} ]]
# then
#     mkdir -p ${AUX_DATADIR}
# fi

# if [[ ! -d ${OUTPUT_DIRECTORY} ]]
# then
#     mkdir -p ${OUTPUT_DIRECTORY}
# fi

# for VAR in elv dir wth upa
# do
#     ptn="*_${VAR}.tif"
#     find ${MERIT_DATADIR} -type f -iname "${ptn}" > /tmp/merit_${VAR}_filenames.txt
#     # if [[ ! -f ${AUX_DATADIR}/merit_${VAR}.tif ]]
#     # then	
#     gdalbuildvrt \
# 	-overwrite \
# 	-te $W $S $E $N \
# 	-tr 0.0008333333333 0.0008333333333 \
# 	-input_file_list /tmp/merit_${VAR}_filenames.txt \
# 	${AUX_DATADIR}/merit_${VAR}.vrt
#     gdal_translate ${AUX_DATADIR}/merit_${VAR}.vrt ${AUX_DATADIR}/merit_${VAR}.tif
#     # fi    
# done

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
