#!/usr/bin/env python
# -*- coding: utf-8 -*-

import xarray as xr
import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd

NBS="degradation_"

GRID_FILE = xr.open_dataset("/mnt/test/clara/mosart-input/rcp45/ACCESS1-0/mosart_new.nc")['area']
GRID_FILE = GRID_FILE.expand_dims("time")

months = np.arange(12)+1
hist_qtm_45_df = pd.DataFrame({'month':months}).set_index("month")
fut_qtm_45_df = pd.DataFrame({'month':months}).set_index("month")
fut_qtm_85_df = pd.DataFrame({'month':months}).set_index("month")

hist_qgm_45_df = pd.DataFrame({'month':months}).set_index("month")
fut_qgm_45_df = pd.DataFrame({'month':months}).set_index("month")
fut_qgm_85_df = pd.DataFrame({'month':months}).set_index("month")

#hist_qngm_45_df = pd.DataFrame({'month':months}).set_index("month")
#fut_qngm_45_df = pd.DataFrame({'month':months}).set_index("month")
#fut_qngm_85_df = pd.DataFrame({'month':months}).set_index("month")

hist_qsm_45_df = pd.DataFrame({'month':months}).set_index("month")
fut_qsm_45_df = pd.DataFrame({'month':months}).set_index("month")
fut_qsm_85_df = pd.DataFrame({'month':months}).set_index("month")

hist_qssm_45_df = pd.DataFrame({'month':months}).set_index("month")
fut_qssm_45_df = pd.DataFrame({'month':months}).set_index("month")
fut_qssm_85_df = pd.DataFrame({'month':months}).set_index("month")

hist_qdm_45_df = pd.DataFrame({'month':months}).set_index("month")
fut_qdm_45_df = pd.DataFrame({'month':months}).set_index("month")
fut_qdm_85_df = pd.DataFrame({'month':months}).set_index("month")

INPUT_FILE = "/home/clara/project/NbS/mosart-input"
RC45 = "rcp45"
RC85 = "rcp85"
SCENARIOS = sorted(os.listdir(os.path.join(INPUT_FILE,RC45)))
FILENAME = "coupled_scpar_runoff.nc"
MOSART_OUTPUT = "/home/clara/project/NbS/mosart-output"

for scenario in SCENARIOS:
    file = xr.open_dataset(os.path.join(INPUT_FILE, RC45, scenario,NBS, FILENAME))
    file['lat'] = GRID_FILE['lat']
    file['lon'] = GRID_FILE['lon']

    #hist_qg_m = file['melt'].sel(time=slice("2000-01-01","2018-12-31")).resample(time="M").mean().groupby("time.month").mean()
    hist_qg_m = file['melt'].sel(time=slice("2002-01-01","2018-12-31")).resample(time="M").mean().groupby("time.month").mean() # Selecting 2002 to 2018 because uh deep has 2 year memory and needs to spinup
    hist_qg_m = np.multiply(hist_qg_m,GRID_FILE[0,:,:]).sum(dim=("lat","lon"))/1000 # kg/m2/s*m2/1000=m3/s

    #hist_qs_m = file['surf_roff'].sel(time=slice("2000-01-01","2018-12-31")).resample(time="M").mean().groupby("time.month").mean()
    hist_qs_m = file['surf_roff'].sel(time=slice("2002-01-01","2018-12-31")).resample(time="M").mean().groupby("time.month").mean()
    hist_qs_m = np.multiply(hist_qs_m,GRID_FILE[0,:,:]).sum(dim=("lat","lon"))/1000 # kg/m2/s*m2/1000=m3/s

    #hist_qss_m = file['sub_surf_roff'].sel(time=slice("2000-01-01","2018-12-31")).resample(time="M").mean().groupby("time.month").mean()
    hist_qss_m = file['sub_surf_roff'].sel(time=slice("2002-01-01","2018-12-31")).resample(time="M").mean().groupby("time.month").mean()
    hist_qss_m = np.multiply(hist_qss_m,GRID_FILE[0,:,:]).sum(dim=("lat","lon"))/1000 # kg/m2/s*m2/1000=m3/s

    #hist_qt_m = (file['surf_roff']+file['sub_surf_roff']).sel(time=slice("2000-01-01","2018-12-31")).resample(time="M").mean().groupby("time.month").mean()
    hist_qt_m = (file['surf_roff']+file['sub_surf_roff']).sel(time=slice("2002-01-01","2018-12-31")).resample(time="M").mean().groupby("time.month").mean()
    hist_qt_m = np.multiply(hist_qt_m,GRID_FILE[0,:,:]).sum(dim=("lat","lon"))/1000 # kg/m2/s*m2/1000=m3/s

#    km105 = pd.read_csv(os.path.join(MOSART_OUTPUT,RC45,scenario,NBS,scenario+"_"+NBS,"deeper_km105.csv"))
#    km105['time'] = pd.to_datetime(km105['time'])
#    km105 = km105.set_index("time")
    #hist_km105 = km105.loc['2000-01-01':'2018-12-31']
#    hist_km105 = km105.loc['2002-01-01':'2018-12-31']
#    hist_km105 = hist_km105.reset_index()
#    hist_km105 = hist_km105.groupby(hist_km105['time'].dt.month).mean()
#    hist_flow_m = pd.DataFrame({'month':months, 'flow':hist_qt_m.values}).set_index("month")

    # Instead since we take the contributing flows from entire catchment, let's take deep contribution derived for total catchment
    deep = pd.read_csv(os.path.join(MOSART_OUTPUT,RC45,scenario,NBS,scenario+"_"+NBS,"deeper_total.csv"))
    deep['time'] = pd.to_datetime(deep['time'])
    deep = deep.set_index("time")
    hist_deep = deep.loc['2002-01-01':'2018-12-31']
    hist_deep = hist_deep.reset_index()
    hist_deep = hist_deep.groupby(hist_deep['time'].dt.month).mean()
    hist_flow_m = pd.DataFrame){'month':months, 'flow':hist_qt_m.values}).set_index("month")

    hist_qgm_45_df[scenario] = hist_qg_m.values
    hist_qtm_45_df[scenario] = hist_flow_m['flow'] + hist_deep['deeper_baseflow'] # assumption of no monthly delay for surface and shallow subsurface
    #hist_qngm_45_df[scenario] = hist_qtm_45_df[scenario] - hist_qgm_45_df[scenario]
    hist_qsm_45_df[scenario] = hist_qs_m.values-hist_qg_m.values
    hist_qssm_45_df[scenario] = hist_qss_m.values
    hist_qdm_45_df[scenario] = hist_deep['deeper_baseflow']


    fut_qg_m = file['melt'].sel(time=slice("2080-01-01","2100-12-31")).resample(time="M").mean().groupby("time.month").mean()
    fut_qg_m = np.multiply(fut_qg_m,GRID_FILE[0,:,:]).sum(dim=("lat","lon"))/1000 # kg/m2/s*m2/1000=m3/s

    fut_qt_m = (file['surf_roff']+file['sub_surf_roff']).sel(time=slice("2080-01-01","2100-12-31")).resample(time="M").mean().groupby("time.month").mean()
    fut_qt_m = np.multiply(fut_qt_m,GRID_FILE[0,:,:]).sum(dim=("lat","lon"))/1000 # kg/m2/s*m2/1000=m3/s

    fut_qs_m = file['surf_roff'].sel(time=slice("2080-01-01","2100-12-31")).resample(time="M").mean().groupby("time.month").mean()
    fut_qs_m = np.multiply(fut_qs_m,GRID_FILE[0,:,:]).sum(dim=("lat","lon"))/1000 # kg/m2/s*m2/1000=m3/s

    fut_qss_m = file['sub_surf_roff'].sel(time=slice("2080-01-01","2100-12-31")).resample(time="M").mean().groupby("time.month").mean()
    fut_qss_m = np.multiply(fut_qss_m,GRID_FILE[0,:,:]).sum(dim=("lat","lon"))/1000 # kg/m2/s*m2/1000=m3/s

#    km105 = pd.read_csv(os.path.join(MOSART_OUTPUT,RC45,scenario,NBS,scenario+"_"+NBS,"deeper_km105.csv"))
#    km105['time'] = pd.to_datetime(km105['time'])
#    km105 = km105.set_index("time")
#    fut_km105 = km105.loc['2080-01-01':'2100-12-31']
#    fut_km105 = fut_km105.reset_index()
#    fut_km105 = fut_km105.groupby(fut_km105['time'].dt.month).mean()
#    fut_flow_m = pd.DataFrame({'month':months, 'flow':fut_qt_m.values}).set_index("month")

    deep = pd.read_csv(os.path.join(MOSART_OUTPUT,RC45,scenario,NBS,scenario+"_"+NBS,"deeper_total.csv"))
    deep['time'] = pd.to_datetime(deep['time'])
    deep = deep.set_index("time")
    fut_deep = deep.loc['2080-01-01':'2100-12-31']
    fut_deep = fut_deep.reset_index()
    fut_deep = fut_deep.groupby(fut_deep['time'].dt.month).mean()
    fut_flow_m = pd.DataFrame({'month':months, 'flow':fut_qt_m.values}).set_index("month")

    fut_qgm_45_df[scenario] = fut_qg_m.values
    fut_qtm_45_df[scenario] = fut_flow_m['flow'] + fut_deep['deeper_baseflow']
    #fut_qngm_45_df[scenario] = fut_qtm_45_df[scenario] - fut_qgm_45_df[scenario]
    fut_qsm_45_df[scenario] = fut_qs_m.values-fut_qg_m.values
    fut_qssm_45_df[scenario] = fut_qss_m.values
    fut_qdm_45_df[scenario] = fut_deep['deeper_baseflow']

    file = xr.open_dataset(os.path.join(INPUT_FILE, RC85, scenario, NBS,FILENAME))
    file['lat'] = GRID_FILE['lat']
    file['lon'] = GRID_FILE['lon']

    fut_qg_m = file['melt'].sel(time=slice("2080-01-01","2100-12-31")).resample(time="M").mean().groupby("time.month").mean()
    fut_qg_m = np.multiply(fut_qg_m,GRID_FILE[0,:,:]).sum(dim=("lat","lon"))/1000 # kg/m2/s*m2/1000=m3/s

    fut_qt_m = (file['surf_roff']+file['sub_surf_roff']).sel(time=slice("2080-01-01","2100-12-31")).resample(time="M").mean().groupby("time.month").mean()
    fut_qt_m = np.multiply(fut_qt_m,GRID_FILE[0,:,:]).sum(dim=("lat","lon"))/1000 # kg/m2/s*m2/1000=m3/s

    fut_qs_m = file['surf_roff'].sel(time=slice("2080-01-01","2100-12-31")).resample(time="M").mean().groupby("time.month").mean()
    fut_qs_m = np.multiply(fut_qs_m,GRID_FILE[0,:,:]).sum(dim=("lat","lon"))/1000 # kg/m2/s*m2/1000=m3/s

    fut_qss_m = file['sub_surf_roff'].sel(time=slice("2080-01-01","2100-12-31")).resample(time="M").mean().groupby("time.month").mean()
    fut_qss_m = np.multiply(fut_qss_m,GRID_FILE[0,:,:]).sum(dim=("lat","lon"))/1000 # kg/m2/s*m2/1000=m3/s

#    km105 = pd.read_csv(os.path.join(MOSART_OUTPUT,RC85,scenario,NBS,scenario+"_"+NBS,"deeper_km105.csv"))
#    km105['time'] = pd.to_datetime(km105['time'])
#    km105 = km105.set_index("time")
#    fut_km105 = km105.loc['2080-01-01':'2100-12-31']
#    fut_km105 = fut_km105.reset_index()
#    fut_km105 = fut_km105.groupby(fut_km105['time'].dt.month).mean()
#    fut_flow_m = pd.DataFrame({'month':months, 'flow':fut_qt_m.values}).set_index("month")

    deep = pd.read_csv(os.path.join(MOSART_OUTPUT,RC85,scenario,NBS,scenario+"_"+NBS,"deeper_total.csv"))
    deep['time'] = pd.to_datetime(deep['time'])
    deep = deep.set_index("time")
    fut_deep = deep.loc['2080-01-01':'2100-12-31']
    fut_deep = fut_deep.reset_index()
    fut_deep = fut_deep.groupby(fut_deep['time'].dt.month).mean()
    fut_flow_m = pd.DataFrame({'month':months, 'flow':fut_qt_m.values}).set_index("month")

    fut_qgm_85_df[scenario] = fut_qg_m.values
    fut_qtm_85_df[scenario] = fut_flow_m['flow'] + fut_deep['deeper_baseflow']
    #fut_qngm_85_df[scenario] = fut_qtm_85_df[scenario] - fut_qgm_85_df[scenario]
    fut_qsm_85_df[scenario] = fut_qs_m.values-fut_qg_m.values
    fut_qssm_85_df[scenario] = fut_qss_m.values
    fut_qdm_85_df[scenario] = fut_deep['deeper_baseflow']


fut_qgm_85_df.to_csv("/home/clara/NbS_paper/processed_flows/"+NBS+"/fut_qgm_85_df.csv")
fut_qtm_85_df.to_csv("/home/clara/NbS_paper/processed_flows/"+NBS+"/fut_qtm_85_df.csv")
fut_qsm_85_df.to_csv("/home/clara/NbS_paper/processed_flows/"+NBS+"fut_qsm_85_df.csv")
fut_qssm_85_df.to_csv("/home/clara/NbS_paper/processed_flows/"+NBS+"fut_qssm_85_df.csv")
fut_qdm_85_df.to_csv("/home/clara/NbS_paper/processed_flows/"+NBS+"fut_qdm_85_df.csv")
#fut_qngm_85_df.to_csv("/home/clara/NbS_paper/fut_qngm_85_df.csv")

fut_qgm_45_df.to_csv("/home/clara/NbS_paper/processed_flows/"+NBS+"/fut_qgm_45_df.csv")
fut_qtm_45_df.to_csv("/home/clara/NbS_paper/processed_flows/"+NBS+"/fut_qtm_45_df.csv")
fut_qsm_45_df.to_csv("/home/clara/NbS_paper/processed_flows/"+NBS+"/fut_qsm_45_df.csv")
fut_qssm_45_df.to_csv("/home/clara/NbS_paper/processed_flows/"+NBS+"/fut_qssm_45_df.csv")
fut_qdm_45_df.to_csv("/home/clara/NbS_paper/processed_flows/"+NBS+"/fut_qdm_45_df.csv")
#fut_qngm_45_df.to_csv("/home/clara/NbS_paper/fut_qngm_45_df.csv")

hist_qgm_45_df.to_csv("/home/clara/NbS_paper/processed_flows/"+NBS+"/hist_qgm_45_df.csv")
hist_qtm_45_df.to_csv("/home/clara/NbS_paper/processed_flows/"+NBS+"/hist_qtm_45_df.csv")
hist_qsm_45_df.to_csv("/home/clara/NbS_paper/processed_flows/"+NBS+"/hist_qsm_45_df.csv")
hist_qssm_45_df.to_csv("/home/clara/NbS_paper/processed_flows/"+NBS+"/hist_qssm_45_df.csv")
hist_qdm_45_df.to_csv("/home/clara/NbS_paper/processed_flows/"+NBS+"/hist_qdm_45_df.csv")
#hist_qngm_45_df.to_csv("/home/clara/NbS_paper/hist_qngm_45_df.csv")
