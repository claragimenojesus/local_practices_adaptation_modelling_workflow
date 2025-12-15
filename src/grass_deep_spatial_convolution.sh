#!/bin/bash

###############################################################################
# GRASS GIS workflow to estimate flow distance to major river network
# and support downstream convolution of deep recharge contributions
#
# Study area: Vilcanota catchment
# Coordinate system: UTM Zone 19S (EPSG:32719)
###############################################################################

#-------------------------------#
# User-defined paths and inputs #
#-------------------------------#

GRASSDB="/path/to/grassdb"
LOCATION="vilcanota_utm"
MAPSET="PERMANENT"

INPUT_DEM_GRASS="dem_filled"                 # Original DEM in GRASS (geographic CRS)
WORKDIR="/path/to/workdir"
DEM_TIF="${WORKDIR}/dem_filled.tif"
DEM_UTM_TIF="${WORKDIR}/dem_filled_utm.tif"

# Stream accumulation thresholds to test
THRESHOLDS=(200 500 1000 10000 20000)

# Selected threshold used in final analysis
FINAL_THRESHOLD=20000

#------------------------------#
# 1. Export and reproject DEM  #
#------------------------------#

# Export depression-filled DEM from GRASS to GeoTIFF
r.out.gdal \
    input=${INPUT_DEM_GRASS} \
    output=${DEM_TIF} \
    format=GTiff \
    --overwrite

# Reproject DEM to UTM Zone 19S (meters required for r.stream.distance)
gdalwarp \
    -t_srs EPSG:32719 \
    ${DEM_TIF} \
    ${DEM_UTM_TIF}

#---------------------------------------------#
# 2. Create new GRASS location and import DEM #
#---------------------------------------------#

# Create a new GRASS location using the UTM DEM
grass -c ${DEM_UTM_TIF} ${GRASSDB}/${LOCATION}

# Launch GRASS session non-interactively
grass ${GRASSDB}/${LOCATION}/${MAPSET} --exec bash << 'EOF'

# Import UTM-projected DEM
r.in.gdal \
    input=/path/to/workdir/dem_filled_utm.tif \
    output=dem_filled_utm \
    --overwrite

# Set computational region to DEM
g.region raster=dem_filled_utm -p

#-------------------------------------------------------------#
# 3. Flow routing: direction, accumulation, stream extraction #
#-------------------------------------------------------------#

# Loop over candidate stream thresholds to assess sensitivity
for THR in 200 500 1000 10000 20000; do

    r.watershed -s \
        elevation=dem_filled_utm \
        threshold=${THR} \
        drainage=flow_dir_${THR} \
        accumulation=flow_accum_${THR} \
        stream=stream_${THR} \
        --overwrite

done

#------------------------------------------------#
# 4. Distance to downstream major river network #
#------------------------------------------------#

# Compute flow distance to selected major river network
r.stream.distance \
    stream_rast=stream_20000 \
    direction=flow_dir_20000 \
    distance=dist_to_river_20000 \
    --overwrite

EOF

###############################################################################
# 5. Post-processing (outside GRASS)
#
# - Export dist_to_river_20000 raster
# - Generate histogram of distances
# - Convolute deep recharge contributions using distance PDF
# - Perform downstream routing and transect aggregation
#
# These steps are implemented in Python (deep_convolution.py).
###############################################################################
