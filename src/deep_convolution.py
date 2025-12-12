#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import xarray as xr
import rioxarray
import glob
from tqdm import tqdm
import os
import matplotlib.pyplot as plt
import scipy.signal
import geopandas as gpd
from rasterio.mask import mask
import rasterio

# Import catchment mask
catchment = gpd.read_file("/home/clara/rahu_data/catchments/cuenca_vilcanota_final.shp")
catchment_utm = catchment.to_crs(epsg=32719)

# Import UH for temporal convolution
def unit_hydro(V,k,n,t,er):
    from scipy.special import gamma
    import numpy as np
    from scipy.sparse import spdiags

    #inputs
    # V = volume of unit hydrograph
    # k,n = parameters from Nash (1957)
    # t = unit hydrograph time
    # er = effective rainfall time-series

    # outputs
    # DGW = deeper groundwater flow
    # uh = unit hydrograph

    uh = np.array(V/(k*gamma(n))*np.multiply(np.exp(-t/k),np.power((t/k),n-1)))
    uh = uh/np.sum(uh)

    # from Ana Mijic as from Anthony script
    m = er.size
    r = t.size
    n = r+m-1

    pmat = spdiags(np.multiply(np.ones(r),np.row_stack(er)),np.arange(0,-er.size,-1),n,r)

    DGW = pmat*uh
    DGW = DGW[0:er.size]

    return [DGW, uh]

# k0 = 18.59720523 for KM105, and k0= 22.32687392 for PISAC
# n0 = 22.57835457 for KM105, and n0= 18.8277447 for PISAC
k0_km105=18.59720523
n0_km105=22.57835457
k0_pisac=22.32687392
n0_pisac=18.8277447
V=1

# Import PDF for spatial convolution
# Using 20,000 threshold contributing cells for stream network
# Load PDF for spatial convolution using 20k threshold

with rasterio.open("/home/clara/dist_to_river_20000.tif") as src:
    # 3. Mask (crop) raster to catchment
    out_image, out_transform = mask(src, catchment_utm.geometry, crop=True)
    out_meta = src.meta.copy()

# 4. Prepare data for plotting
dist_masked = out_image[0]
dist_masked = np.where(dist_masked < 0, np.nan, dist_masked)  # mask invalid values

dist_km = dist_masked[~np.isnan(dist_masked)] /1000
bin_edges = np.linspace(0, 346.7, 100)
hist_counts, bin_edges  = np.histogram(dist_km, bins=bin_edges)
pdf = hist_counts / hist_counts.sum()
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

# PErforming operations for future degradation rcp4.5 first

# FUTURE DEGRADATION RCP4.5 - TAKING MEDIAN ACROSS ALL FUTURE SCENARIOS

base_path = "/home/clara/project/NbS/jules-output/rcp45"
gcm_list = sorted(os.listdir(base_path))
filename_template = "rcp45_{gcm}_degradation_coupled_jules_oggm_00_99.nc"

hyd_cond = xr.open_dataset("/home/jec18/gis_layer/groundwater/hydrogeo_k.nc").fillna(0)

resampled_dir = "/home/clara/NbS_paper/resampled_basins/"
n_basins=100
area_m2 = 15675183.637074107 # m2

total_area_m2 = np.load(os.path.join(resampled_dir, "total_area_m2.npy"))
watershed_area = np.load(os.path.join(resampled_dir, "watershed_area_km2.npy"))

surf_matrix = np.zeros((n_basins, len(gcm_list)))
subsurface_matrix = np.zeros((n_basins, len(gcm_list)))
melt_matrix = np.zeros((n_basins, len(gcm_list)))
deep_matrix = np.zeros((n_basins, len(gcm_list)))

deep_dry_conv = np.zeros((n_basins,len(gcm_list)))

#time_index = pd.to_datetime(hist.time.values)
# Dry season mask
#dry_mask = (time_index >= "2080-01-01") & (time_index <= "2099-12-31") & (time_index.month >=5) & (time_index.month <=11)

for idx, gcm in enumerate(gcm_list):
    print(f"Processing scenario {idx+1}/30: {gcm}")

    nc_path = os.path.join(base_path, gcm, filename_template.format(gcm=gcm))

    file = xr.open_dataset(nc_path)[["surf_roff", "sub_surf_roff", "melt"]].sel(
        time=slice("2077-01-01", "2099-12-31"), drop=True
    )

    # Repeat hydraulic conductivity over time
    hyd_cond0 = xr.zeros_like(file.sub_surf_roff)
    hyd_vals = hyd_cond.Band1.values
    for t in range(file.sizes["time"]):
        hyd_cond0[t, ...] = hyd_vals

    q_deep = xr.zeros_like(file.sub_surf_roff)
    q_deep[...] = np.minimum(file.sub_surf_roff, 1000 * hyd_cond0)
    file["sub_surf_roff"].values -= q_deep

    # Subtract glacier melt from surface runoff
    file["surf_roff"] -= file["melt"]

    q_deep=q_deep.rio.set_spatial_dims(x_dim="lon",y_dim="lat",inplace=True)
    q_deep.rio.write_crs("epsg:4326",inplace=True)

    # convert fluxes to mm/d
    file = file*3600*24
    q_deep = q_deep*3600*24

    file = file.rio.write_crs("EPSG:4326")

    time_index = pd.to_datetime(file.time.values)
    # Dry season mask
    dry_mask = (time_index >= "2080-01-01") & (time_index <= "2099-12-31") & (time_index.month >=5) & (time_index.month <=11)

    for i in range(n_basins):
        basins_resampled = np.load(os.path.join(resampled_dir, f"basin_{i:03d}.npy"))

        # Extract daily timeseries per flux
        surf_ts = (file['surf_roff'].values * basins_resampled).sum(axis=(1, 2)) * area_m2 / total_area_m2[i]
        subsurf_ts = (file['sub_surf_roff'].values * basins_resampled).sum(axis=(1, 2)) * area_m2 / total_area_m2[i]
        melt_ts = (file['melt'].values * basins_resampled).sum(axis=(1, 2)) * area_m2 / total_area_m2[i]
        deep_ts = (q_deep.values * basins_resampled).sum(axis=(1,2)) * area_m2 / total_area_m2[i]
        # Temporal convolution for each sub-basin
        deep_daily_conv, uh = unit_hydro(V, k0_km105, n0_km105, np.arange(len(deep_ts)), deep_ts)

        surf_matrix[i, idx] = surf_ts[dry_mask].mean()
        subsurface_matrix[i, idx] = subsurf_ts[dry_mask].mean()
        melt_matrix[i, idx] = melt_ts[dry_mask].mean()
        deep_matrix[i, idx] = deep_daily_conv[dry_mask].mean()

    # Spatial convolution
    deep_dry_conv[:,idx] = scipy.signal.convolve(deep_matrix[:,idx], pdf, mode="full")[:n_basins]

# Get mean of all scenarios for historic
surface_vector = np.nanmedian(surf_matrix,axis=1)
subsurface_vector = np.nanmedian(subsurface_matrix,axis=1)
melt_vector = np.nanmedian(melt_matrix,axis=1)
deep_vector = np.nanmedian(deep_dry_conv, axis=1)

# Convolution of deep component
distance = (np.arange(1, 101) * 346.7) / 100  # same as (1:100) * 346.7 / 100 in R

# Build DataFrame
fut_45_df = pd.DataFrame({
    'watershed_area': watershed_area,
    'distance': distance,
    'surf_vector': surface_vector,
    'subsurface_vector': subsurface_vector,
    'melt_vector': melt_vector,
    'deep_vector': deep_vector
})

fut_45_df.to_csv("/home/clara/NbS_paper/transect_figures/fut45_dry_season_depth_df.csv", index=False)

distance = fut_45_df['distance']
subsurface = fut_45_df['subsurface_vector']
surface = fut_45_df['surf_vector']
melt = fut_45_df['melt_vector']
deep = fut_45_df['deep_vector']

# Compute cumulative sums for stacking
bottom = deep
middle = bottom + subsurface
top = middle + surface
top2 = top + melt

# Set up journal-compliant styling
plt.rcParams.update({
    "axes.labelsize": 10,
    "axes.titlesize": 9,
    "legend.fontsize": 9,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
})


# Create figure
fig, ax = plt.subplots(figsize=(12,7))  # Width = 180 mm, height adjustable

# Plot stacked areas
ax.fill_between(distance, 0, bottom, label='Deep sub-surface runoff', color='royalblue', alpha=0.8)
ax.fill_between(distance, bottom, middle, label='Shallow sub-surface runoff', color='dodgerblue', alpha=0.8)
ax.fill_between(distance, middle, top, label='Surface runoff', color='skyblue', alpha=0.8)
ax.fill_between(distance, top, top2, label='Glacier melt', color='lightcoral', alpha=0.8)

# Labels
ax.set_xlabel('Distance from headwater (km)', size=10, labelpad=10)
ax.set_ylabel('Runoff depth (mm d⁻¹)',size=10, labelpad=10)

# Axis and grid
ax.set_ylim(0, 1.1)
ax.grid(True, linestyle='--', linewidth=0.4, alpha=0.6)


# Legend
ax.legend(loc='upper right', frameon=False, fontsize=10)

# Save as high-resolution vector graphic
plt.tight_layout()
plt.title("BASELINE SCENARIO MEAN DRY SEASON FLUXES - RCP4.5 ENSEMBLE MEAN (2080-2100)", fontsize=11, weight="bold")
plt.savefig("/home/clara/NbS_paper/transect_figures/fut45_dry_season_depth.png", dpi=600)  # or .svg
