#!/bin/bash

CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda deactivate

export SRC_DIR=$(pwd)/../src
conda activate mosart
python ${SRC_DIR}/mosart-rahu.py # File mosart-rahu needs to be adapted for each scenario to point to the correct MOSART configuration file
conda deactivate
