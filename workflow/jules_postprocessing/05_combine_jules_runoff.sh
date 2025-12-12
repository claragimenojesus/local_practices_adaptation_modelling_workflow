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
export INPUT_DIRECTORY=$(yaml $CONFIG_FILE "['combined_runoff']['input_directory']")
export INPUT_FILE_SUFFIX=$(yaml $CONFIG_FILE "['combined_runoff']['input_file_suffix']")
export OUTPUT_FILENAME=$(yaml $CONFIG_FILE "['combined_runoff']['output_filename']")
# export OUTPUT_DIRECTORY=$(yaml $CONFIG_FILE "['mosart']['output_directory']")

# =================================================================== #
# Enable Anaconda to be activated/deactivated within the script
# =================================================================== #

CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh

# Deactivate
conda deactivate

# =================================================================== #
# Aggregate JULES output in time
# =================================================================== #

CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda deactivate

export SRC_DIR=$(pwd)/../src
conda activate mosart
python ${SRC_DIR}/merged-runoff.py
conda deactivate 
