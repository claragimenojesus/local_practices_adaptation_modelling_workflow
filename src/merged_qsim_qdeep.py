#!/usr/bin/env python
# -*- coding: utf-8 -*-

import xarray as xr
import os
import pandas as pd
import re

INPUT_DIRECTORY = str(os.environ['MOSART_INPUT'])
INPUT_FILENAME = str(os.environ['INPUT_FILENAME'])
OUTPUT_DIRECTORY = str(os.environ['MOSART_OUTPUT'])
RCP_DIR = str(os.environ['RCP_DIR'])
ANALYSIS_DIR = str(os.environ['ANALYSIS_DIR'])
RC = str(os.environ['RC'])
NBS = str(os.environ['NBS'])
SCENARIO_SH = str(os.environ['SCENARIO_SH'])
#SC = re.findall('5/(.*)',SCENARIO_SH)[0]
SC = str(os.environ['SC'])

qsim_km105 = pd.read_csv(os.path.join(OUTPUT_DIRECTORY, NBS,SC+"_"+NBS, "output_pt_km105.csv"))
qdeep_km105 = pd.read_csv(os.path.join(OUTPUT_DIRECTORY, NBS,SC+"_"+NBS,"deeper_km105.csv"))

df_km105 = pd.DataFrame({"time" : qsim_km105['time'], "mosart_flow" : qsim_km105['RIVER_DISCHARGE_OVER_LAND_LIQ'], "deeper_baseflow" : qdeep_km105['deeper_baseflow'], "total_flow" : qsim_km105['RIVER_DISCHARGE_OVER_LAND_LIQ']+qdeep_km105['deeper_baseflow']})

df_km105.to_csv(os.path.join(OUTPUT_DIRECTORY, NBS,SC+"_"+NBS,"qtotal_km105.csv"), index=False)

qsim_pisac = pd.read_csv(os.path.join(OUTPUT_DIRECTORY,NBS,SC+"_"+NBS, "output_pt_Pisac.csv"))
qdeep_pisac = pd.read_csv(os.path.join(OUTPUT_DIRECTORY, NBS,SC+"_"+NBS,"deeper_pisac.csv"))

df_pisac = pd.DataFrame({"time" : qsim_pisac['time'], "mosart_flow" : qsim_pisac['RIVER_DISCHARGE_OVER_LAND_LIQ'], "deeper_baseflow" : qdeep_pisac['deeper_baseflow'], "total_flow" : qsim_pisac['RIVER_DISCHARGE_OVER_LAND_LIQ']+qdeep_pisac['deeper_baseflow']})

df_pisac.to_csv(os.path.join(OUTPUT_DIRECTORY, NBS,SC+"_"+NBS, "qtotal_pisac.csv"), index=False)

# Script to obtain all routed flows in a single csv file to then be able to plot it
# one file for rcp45 scenarios with dimensions timexscenario and one file for rcp85 scenarios
