#!/bin/bash

export GRASS_MESSAGE_FORMAT=plain

# remove any existing mask
r.mask -r

# # Import data
# for VAR in elv dir wth upa
# do
#     r.in.gdal \
# 	input=${AUX_DATADIR}/merit_${VAR}.vrt \
# 	output=merit_${VAR} \
# 	--overwrite
# done

r.in.gdal \
    input=${AUX_DATADIR}/geomorpho90m_slope.vrt \
    output=merit_slope \
    --overwrite

# This map derived with pyflwdir
gdal_translate \
    -a_srs EPSG:4326 \
    ${AUX_DATADIR}/merit_basin.tif \
    ${AUX_DATADIR}/merit_basin_ll.tif

mv ${AUX_DATADIR}/merit_basin_ll.tif ${AUX_DATADIR}/merit_basin.tif
r.in.gdal \
    input=${AUX_DATADIR}/merit_basin.tif \
    output=merit_basin \
    --overwrite

# g.region raster=merit_elv
# g.region -p

# # Land fraction
# r.mapcalc \
#     "land_frac = if(isnull(merit_elv), 0, 1)" \
#     --overwrite

# for VAR in elv dir wth upa
# do    
#     r.out.gdal\
# 	input=merit_${VAR} \
# 	output=${AUX_DATADIR}/merit_rgn_${VAR}.tif \
# 	--overwrite
# done

# Set region to coarse model resolution
g.region raster=merit_slope res=0:06
g.region -p

# # Compute land fraction
# r.resamp.stats \
#     input=land_frac \
#     output=land_frac_resamp \
#     method=average \
#     --overwrite

# r.out.gdal \
#     input=land_frac_resamp \
#     output=${AUX_DATADIR}/merit_rgn_land_fraction_lr.tif \
#     --overwrite

# r.resamp.stats \
#     input=merit_elv \
#     output=merit_elv_resamp \
#     method=average \
#     --overwrite

# r.out.gdal \
#     input=merit_elv_resamp \
#     output=${AUX_DATADIR}/merit_rgn_elv_lr.tif \
#     --overwrite

# Set mask so that the following commands only take
# the average of pixels within the basin.
r.mask -r
r.mask raster=merit_basin
r.resamp.stats \
    input=merit_slope \
    output=merit_slope_resamp \
    method=average \
    --overwrite

r.out.gdal \
    input=merit_slope_resamp \
    output=${AUX_DATADIR}/merit_basin_slope_lr.tif \
    --overwrite

# Remove the mask, to calculate grid cell area

# TODO
# # Manning's n - overland 
# r.in.gdal \
#     input=/mnt/scratch/scratch/data/ESA_CCI_LC//ESACCI-LC-L4-LCCS-Map-300m-P1Y-2015-v2.0.7.tif \
#     output=esacci_lc \
#     --overwrite
# TODO: lookup table based on literature values (for now just use bulk value)
# r.mapcalc "n = 0.03" --overwrite
# r.out.gdal \
#     input=n \
#     output=${AUX_DATADIR}/manning_n_overland_0.100000Deg.tif \
#     --overwrite

# 11;Herbaceous cover
# 12;Tree or shrub cover
# 20;Cropland, irrigated or post-flooding
# 30;Mosaic cropland (>50%) / natural vegetation (tree, shrub, herbaceous cover) (<50%)
# 40;Mosaic natural vegetation (tree, shrub, herbaceous cover) (>50%) / cropland (<50%)
# 50;Tree cover, broadleaved, evergreen, closed to open (>15%)
# 60;Tree cover, broadleaved, deciduous, closed to open (>15%)
# 61;Tree cover, broadleaved, deciduous, closed (>40%)
# 62;Tree cover, broadleaved, deciduous, open (15-40%)
# 70;Tree cover, needleleaved, evergreen, closed to open (>15%)
# 71;Tree cover, needleleaved, evergreen, closed (>40%)
# 72;Tree cover, needleleaved, evergreen, open (15-40%)
# 80;Tree cover, needleleaved, deciduous, closed to open (>15%)
# 81;Tree cover, needleleaved, deciduous, closed (>40%)
# 82;Tree cover, needleleaved, deciduous, open (15-40%)
# 90;Tree cover, mixed leaf type (broadleaved and needleleaved)
# 100;Mosaic tree and shrub (>50%) / herbaceous cover (<50%)
# 110;Mosaic herbaceous cover (>50%) / tree and shrub (<50%)
# 120;Shrubland
# 121;Shrubland evergreen
# 122;Shrubland deciduous
# 130;Grassland
# 140;Lichens and mosses
# 150;Sparse vegetation (tree, shrub, herbaceous cover) (<15%)
# 151;Sparse tree (<15%)
# 152;Sparse shrub (<15%)
# 153;Sparse herbaceous cover (<15%)
# 160;Tree cover, flooded, fresh or brakish water
# 170;Tree cover, flooded, saline water
# 180;Shrub or herbaceous cover, flooded, fresh/saline/brakish water
# 190;Urban areas
# 200;Bare areas
# 201;Consolidated bare areas
# 202;Unconsolidated bare areas
# 210;Water bodies
# 220;Permanent snow and ice


# # Reproject topographic slope
# r.proj \
#     location=utm19s \
#     mapset=jules \
#     input=merit_slope \
#     output=merit_slope \
#     method=bilinear \
#     --overwrite

# ACCUMULATION:

# # This results in a segmentation fault:
# r.accumulate \
#     direction=merit_grass_dir \
#     format=45degree \
#     accumulation=merit_flow_accumulation \
#     --overwrite
# r.out.gdal \
#     input=merit_flow_accumulation \
#     output=merit_flow_accumulation.tif \
#     --overwrite

# # Alternative
# r.watershed \
#     -s \
#     elevation=merit_elv \
#     accumulation=merit_accum \
#     drainage=merit_draindir \
#     --overwrite

# r.out.gdal \
#     input=merit_accum \
#     output=merit_flow_accumulation.tif \
#     --overwrite

# # r.out.gdal \
# #     input=merit_flow_accumulation \
# #     output=merit_flow_accumulation.tif \
# #     --overwrite


# # Reclass MERIT to GRASS

# # MERIT:      GRASS:
# # 32 64 128   3 2 1 
# # 16    1     4   8
# # 8  4  2     5 6 7
# # NB: '-1' indicates inland depression in both maps
# #     '0' indicates a river mouth in MERIT (in GRASS?)
# echo "1 = 8
# 2       = 7      
# 4       = 6
# 8       = 5
# 16      = 4
# 32      = 3
# 64      = 2
# 128     = 1
# *       = -1
# " > /tmp/merit_to_grass_draindir_rcl.txt

# r.reclass \
#     input=merit_dir \
#     output=merit_grass_dir \
#     rules=/tmp/merit_to_grass_draindir_rcl.txt \
#     --overwrite

# # Delineate Vilcanota basin
# r.water.outlet \
#     input=merit_grass_dir \
#     output=merit_vilcanota_basin \
#     coordinates=-72.6271,-13.0045 \
#     --overwrite

# r.out.gdal \
#     input=merit_vilcanota_basin \
#     output=${AUX_DATADIR}/merit_vilcanota_basin.tif \
#     --overwrite
