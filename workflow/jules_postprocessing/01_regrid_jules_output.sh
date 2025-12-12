#!/bin/bash

yaml() {
    python3 -c "import yaml;print(yaml.safe_load(open('$1'))$2)"
}

CONFIG_FILE="$1"

export JULES_START_YEAR=$(yaml $CONFIG_FILE "['jules']['start_year']")
export JULES_END_YEAR=$(yaml $CONFIG_FILE "['jules']['end_year']")
export JULES_SUITE=$(yaml $CONFIG_FILE "['jules']['suite']")
export JULES_ID_STEM=$(yaml $CONFIG_FILE "['jules']['id_stem']")
export JULES_JOB_NAME=$(yaml $CONFIG_FILE "['jules']['job_name']")
export JULES_PROFILE_NAME=$(yaml $CONFIG_FILE "['jules']['profile_name']")
export JULES_GRIDFILE=$(yaml $CONFIG_FILE "['jules']['gridfile']")
export JULES_OUTPUT_DIRECTORY=$(yaml $CONFIG_FILE "['jules']['jules_output_directory']")
export JULES_Y_DIM_NAME=$(yaml $CONFIG_FILE "['jules']['y_dim_name']")
export JULES_X_DIM_NAME=$(yaml $CONFIG_FILE "['jules']['x_dim_name']")
export JULES_MASK_VAR_NAME=$(yaml $CONFIG_FILE "['jules']['mask_var_name']")
export JULES_SOIL_DIM_NAME=$(yaml $CONFIG_FILE "['jules']['soil_dim_name']")
export JULES_TILE_DIM_NAME=$(yaml $CONFIG_FILE "['jules']['tile_dim_name']")
export JULES_PFT_DIM_NAME=$(yaml $CONFIG_FILE "['jules']['pft_dim_name']")
export REGRID_OUTPUT_DIRECTORY=$(yaml $CONFIG_FILE "['regrid_jules_output']['output_directory']")
export REGRID_FILE_SUFFIX=$(yaml $CONFIG_FILE "['regrid_jules_output']['file_suffix']")

CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda deactivate

export SRC_DIR=$(pwd)/../src
conda activate mosart
python ${SRC_DIR}/regrid-jules-output.py
conda deactivate 
