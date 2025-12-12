#!/bin/bash

Help()
{
    echo "Script to..."
    echo
    echo "Syntax: 02_aggregate_jules_output.sh [-h] my-config.yaml"
    echo
    echo "[-h | --help]          Print this help message."
    echo "my-config.yaml         YAML configuration file."
    echo
}

POSITIONAL=()
while [[ $# -gt 0 ]]
do
    key="$1"
    case $key in
	-h|--help)
	    Help
	    exit
	    ;;
	*)  # unknown option
	    POSITIONAL+=("$1") # save it in an array for later
	    shift # past argument
	    ;;
    esac
done

CONFIG_FILE="${POSITIONAL[0]}"

yaml() {
    python3 -c "import yaml;print(yaml.safe_load(open('$1'))$2)"
}

export JULES_START_YEAR=$(yaml $CONFIG_FILE "['jules']['start_year']")
export JULES_END_YEAR=$(yaml $CONFIG_FILE "['jules']['end_year']")
export JULES_SUITE=$(yaml $CONFIG_FILE "['jules']['suite']")
export JULES_ID_STEM=$(yaml $CONFIG_FILE "['jules']['id_stem']")
export JULES_JOB_NAME=$(yaml $CONFIG_FILE "['jules']['job_name']")
export JULES_PROFILE_NAME=$(yaml $CONFIG_FILE "['jules']['profile_name']")
export INPUT_DIRECTORY=$(yaml $CONFIG_FILE "['aggregate_jules_output']['input_directory']")
export INPUT_FILE_SUFFIX=$(yaml $CONFIG_FILE "['aggregate_jules_output']['input_file_suffix']")
export OUTPUT_DIRECTORY=$(yaml $CONFIG_FILE "['aggregate_jules_output']['output_directory']")

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
python ${SRC_DIR}/aggregate-jules.py
conda deactivate 
