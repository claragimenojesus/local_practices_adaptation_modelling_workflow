#!/usr/bin/env Rscript

# Script to rasterize JULES-OGGM outputs for CMIP5 simulations to JULES VU catchment grid
# Created Feb 2024 - Clara Gimeno Jésus

library(raster)
library(dplyr)
library(magrittr)
library(MASS)

# Read JULES_OGGM glacier center coordinates
pts <- read.csv("/home/clara/rahu_data/JULES_OGGM/cmip_runs_combined_v3/latlon_532.csv")
pts2 <- pts[,2:3]
pts2 <- relocate(pts2, 'Longitude')
pts2 <- as.matrix(pts2)

dir_areas <- list.files('/home/clara/rahu_data/JULES_OGGM/cmip_runs_combined_v3', pattern="area_m2.csv", full.names=TRUE) # index 1 to 30 rcp45, and 31 to 60 rcp85
dir_rcp45 <- list.dirs('/mnt/scratch2/clara/JULES_OGGM/rcp45')[-1]
dir_rcp85 <- list.dirs('/mnt/scratch2/clara/JULES_OGGM/rcp85')[-1]

# Read JULES-hydro gridcell area size
r <- raster("/home/clara/rahu_data/netcdf/gridcell_area_rahu_v2.nc")

# Let's do the same for rcp85 scenarios
for (i in 1:2) {
    area_oggm <- read.csv(dir_areas[i+30])
    names(area_oggm) <- NULL

    # Preallocate 3D dataset for glacier areas
    DF <- array(data=NA, dim=c(101,75,102))
    for (year in 1:101) {
        # Read each year of data
        area <- as.numeric(area_oggm[year,])

        # Transfer values associated with 'object' type spatial data (points, lines, polygons) to raster cells. If x represents points, each point is assigned to a grid cell. Points that fall on a border between cells are placed in the cell to the right and/or in the cell below. The value of a grid cell is determined by the values associated with the points and function fun.
        # Sum values of glacier areas to contributing JULES-hydro gridcell.
        rnew <- rasterize(pts2, r, area, fun=sum, na.rm=TRUE)
        df <- as.matrix(rnew, xy=TRUE)

        # Store in 3D array
        DF[year,,] <- df
    }
    save(DF, file=paste(dir_rcp85[i],"yearly_areas.rda",sep="/"))

    # Now create % glacier area per gridcell area 3d array
    perc <- array(data=NA, dim=c(101,75,102))
    for (year in 1:101) {
        for (lat in 1:75) {
            for (lon in 1:102) {
                perc[year,lat,lon] = DF[year,lat,lon]/r[lat,lon,1]
            }
        }
    }

    save(perc, file=paste(dir_rcp85[i],"percentage_glacier_per_gridcell_area.rda",sep="/"))
}

# Doing glacier melt daily file
for (i in 1:30) {
    melt_oggm <- read.csv(paste(dir_rcp85[i],"daily_melt.csv",sep="/"))

    DFm <- array(data=NA, dim=c(36525,75,102))
    for (day in 1:36525) {
        melt <- melt_oggm[,day+1]
        rnew <- rasterize(pts2, r, melt, na.rm=TRUE, fun=sum)
        dfmelt <- as.matrix(rnew, xy=TRUE)
        DFm[day,,] <- dfmelt
    }

    save(DFm, file=paste(dir_rcp85[i],"daily_melt_gridded.rda",sep="/"))

    # Doing daily rain file
    rain_oggm <- read.csv(paste(dir_rcp85[i],"daily_rain.csv",sep="/"))

    DFr <- array(data=NA, dim=c(36525,75,102))
    for (day in 1:36525) {
        rain <- rain_oggm[,day+1]
        rnew_rain <- rasterize(pts2, r, rain, na.rm=TRUE, fun=sum)
        dfrain <- as.matrix(rnew_rain, xy=TRUE)
        DFr[day,,] <- dfrain
    }

    save(DFr, file=paste(dir_rcp85[i],"daily_rain_gridded.rda",sep="/"))
}
