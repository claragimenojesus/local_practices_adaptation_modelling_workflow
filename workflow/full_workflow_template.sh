#!/bin/bash
#PBS -l select=1:ncpus=23:mem=25gb:ompthreads=8
#PBS -l walltime=13:00:00
#PBS -J 1-30

eval "$(~/anaconda3/bin/conda shell.bash hook)"
source activate mosart

cd $PBS_O_WORKDIR

export NBS=""

CONFIG_FILE=$(ls -dv rahu-config_files_rcp45/* | head -n $PBS_ARRAY_INDEX | tail -n 1)
MOSART_FILE=$(ls -dv config_rahu_files_rcp45/mosart* | head -n $PBS_ARRAY_INDEX | tail -n 1)
CONFIG_FILE_METEO=$(ls -dv rahu-config_files_rcp45_meteo/* | head -n $PBS_ARRAY_INDEX | tail -n 1)

export SCENARIO_SH=$(ls -dv  "$EPHEMERAL/"JULES_OGGM/rcp45/* | head -n $PBS_ARRAY_INDEX | tail -n 1)
export SCENARIO_OUTPUT_SH=$(ls -dv "$EPHEMERAL/"jules-output/rcp45/* | head -n $PBS_ARRAY_INDEX | tail -n 1)
export HIST_OUTPUT_SH=$(ls -dv "$HOME/"project/rahu/jules-output/u-cz655)
export SCENARIO=$(basename $SCENARIO_SH)
export RC="rcp45"

echo "This script will pre-process rcp85 scenarios for MOSART routing"
echo "Regridding jules output"
./01_regrid_jules_output.sh $CONFIG_FILE
./01_regrid_jules_output.sh $CONFIG_FILE_METEO


export SRC_DIR=$(pwd)/../src
echo "Coupling with JULES-OGGM"
python ${SRC_DIR}/coupling_pt1.py


cd $SCENARIO_SH
cdo -yearadd -gtc,100000000 -inttime,2000-01-01,00:00:00,1day glacier_area_yyyy.nc glacier_area_yyyy.nc glacier_area_daily.nc
cdo -yearadd -gtc,100000000 -inttime,2000-01-01,00:00:00,1day glacier_perc_per_gridcell.nc glacier_perc_per_gridcell.nc glacier_perc_daily.nc

# then remove last day of data
cdo -delete,date="2100-01-01T00:00:00" glacier_area_daily.nc glacier_area_daily_v2.nc
cdo -delete,date="2100-01-01T00:00:00" glacier_perc_daily.nc glacier_perc_daily_v2.nc

# also remove files not needed anymore
rm glacier_area_yyyy.nc
rm glacier_area_yyyy_v2.nc
rm glacier_perc_per_gridcell.nc
rm glacier_perc_per_gridcell_v2.nc

cd $PBS_O_WORKDIR

python ${SRC_DIR}/coupling_pt2.py

cd $SCENARIO_SH
# Delete files not needed
rm glacier_area_daily.nc
rm glacier_area_daily_v2.nc
rm glacier_perc_daily.nc
rm glacier_perc_daily_v2.nc
rm rain_daily_kgm2s.nc

cd $PBS_O_WORKDIR
echo "Resampling jules output"
./03c_resample_jules_output.sh $CONFIG_FILE
echo "Resampling complete"
echo "Computing jules annual mean runoff"
./04c_compute_jules_annual_mean_runoff.sh $CONFIG_FILE
echo "Annual mean runoff complete"
echo "Combining jules runoff"
./05c_combine_jules_runoff.sh $CONFIG_FILE
echo "Combined jules runoff complete"
echo "Create mosart input"
./07_create-mosartwm-input.sh $CONFIG_FILE
echo "Mosart input complete"
./08_scaling_and_partitioning.sh $CONFIG_FILE

echo "Running Mosart"
NUMBA_NUM_THREADS=32 ./$MOSART_FILE
echo "Finished routing with Mosart"
