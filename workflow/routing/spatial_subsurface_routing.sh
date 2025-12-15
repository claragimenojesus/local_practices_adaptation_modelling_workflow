#!/bin/bash
#PBS -l select=1:ncpus=23:mem=25gb:ompthreads=8
#PBS -l walltime=13:00:00
#PBS -J 1-30

eval "$(~/anaconda3/bin/conda shell.bash hook)"
source activate mosart

cd $PBS_O_WORKDIR

export CATCH_MASK_PATH = ".../cuenca_vilcanota.shp"
export DIST_TO_RIVER_PATH = ".../dist_to_river_20000.tif"
export JULES_OUTPUT = ".../jules-output/"
export RC = "rcp45" # or rcp85
export NBS = ""
export HYD_COND = ".../hydrogeo_k.nc"

# Extract a series of basins along river transect of choice in GRASS and store. In our case, 100 equidistant points along the longest glaciated river transect in the catchment were extracted.
export RESAMPLED_BASINS = ".../resampled_basins/" 
export OUTPUT_PATH = ".../"

export SRC_DIR=$(pwd)/../src
python ${SRC_DIR}/deep_convolution.py
