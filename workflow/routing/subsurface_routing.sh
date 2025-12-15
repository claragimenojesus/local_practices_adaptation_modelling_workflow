#!/bin/bash
#PBS -l select=1:ncpus=20:mem=20gb
#PBS -l walltime=03:00:00

#module load anaconda3/personal
#source activate mosart
eval "$(~/anaconda3/bin/conda shell.bash hook)"
source activate mosart

cd $PBS_O_WORKDIR

export INPUT_FILENAME="coupled_scpar_runoff.nc"
export MOSART_INPUT=$(ls -dv "$EPHEMERAL/"mosart-input/rcp45/* | head -n 2 | tail -n 1)
export MOSART_OUTPUT=$(ls -dv "$EPHEMERAL/"mosart-output/rcp45/* | head -n 2 | tail -n 1)
export RCP_DIR="$EPHEMERAL/mosart-output/rcp45"
export ANALYSIS_DIR="$HOME/analysis_files"
export RC="rcp45"
export SCENARIO_SH=$(ls -dv "$EPHEMERAL/"jules-output/rcp45/* | head -n 2 | tail -n 1)
export SC=$(basename $SCENARIO_SH)

export NBS="" # insert relevant NbS scenario
eval "$(~/anaconda3/bin/conda shell.bash hook)"
source activate mosart

export SRC_DIR=$(pwd)/src
python ${SRC_DIR}/uh.py

