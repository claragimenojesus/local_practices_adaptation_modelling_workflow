#!/bin/bash


yaml() {
    python3 -c "import yaml;print(yaml.safe_load(open('$1'))$2)"
}

CONFIG_FILE="$1"


export INPUT_DIRECTORY=$(yaml $CONFIG_FILE "['scaling_and_partitioning']['input_directory']")
export INPUT_FILENAME=$(yaml $CONFIG_FILE "['scaling_and_partitioning']['input_filename']")
export MOSART_MASK=$(yaml $CONFIG_FILE "['scaling_and_partitioning']['mosart_mask']")
export LAND_FRAC=$(yaml $CONFIG_FILE "['scaling_and_partitioning']['land_frac']")
export GRID_AREA=$(yaml $CONFIG_FILE "['scaling_and_partitioning']['grid_area']")
export MOSART_AREA=$(yaml $CONFIG_FILE "['scaling_and_partitioning']['mosart_area']")
export HYD_COND=$(yaml $CONFIG_FILE "['scaling_and_partitioning']['hyd_cond']")
export OUTPUT_DIRECTORY=$(yaml $CONFIG_FILE "['scaling_and_partitioning']['output_directory']")

CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda deactivate

export SRC_DIR=$(pwd)/../src
conda activate mosart
python ${SRC_DIR}/scaling_and_partitioning.py
conda deactivate
