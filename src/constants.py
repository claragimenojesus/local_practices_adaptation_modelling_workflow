#!/usr/bin/env python
# -*- coding: utf-8 -*-

OUTPUT_VARS = {
    'daily_hydrology' : [
        'ecan_gb',                  # Gridbox mean evaporation from canopy/surface store (kg m-2 s-1)
        'esoil_gb',                 # Gridbox surface evapotranspiration from soil moisture store (kg m-2 s-1)
        'et_stom_gb',               # Gridbox stomatal transpiration (kg m-2 s-1)
        'fqw_gb',                   # Gridbox moisture flux from surface (kg m-2 s-1)
        'precip',                   # Gridbox precipitation rate (kg m-2 s-1)
        'rainfall',                 # Gridbox rainfall rate (kg m-2 s-1)
        'snowfall',                 # Gridbox snowfall rate (kg m-2 s-1)
        'canopy_gb',                # Gridbox canopy water content (kg m-2 s-1)
        'drain',                    # Drainage from bottom (nshyd) soil layer (kg m-2 s-1)
        'drain_soilt',              # Drainage from bottom (nshyd) soil layer on soil tiles (kg m-2 s-1)
        'elake',                    # Gridbox mean evaporation from lakes (kg m-2 s-1)
        'fsat',                     # Surface saturated fraction (1)
        'fsat_soilt',               # Surface saturated fraction on soil tiles (1)
        'fwetl',                    # Wetland fraction at end of model timestep (1)
        'fwetl_soilt',              # Wetland fraction at end of model timestep on soil tiles (1)
        'qbase',                    # Baseflow (lateral subsurface runoff) (kg m-2 s-1)
        'qbase_soilt',              # Baseflow (lateral subsurface runoff) on soil tiles (kg m-2 s-1)
        'qbase_zw',                 # Baseflow from deep LSH/TOPMODEL layer (kg m-2 s-1)
        'qbase_zw_soilt',           # Baseflow from deep LSH/TOPMODEL layer on soil tiles (kg m-2 s-1)
        'runoff',                   # Gridbox runoff rate (kg m-2 s-1)
        'smc_avail_top',            # Gridbox available moisture in top 1.000000m of soil (kg m-2)
        'smc_avail_tot',            # Gridbox available moisture in soil column (kg m-2)
        'smc_tot',                  # Gridbox total soil moisture in column (kg m-2)
        'sthzw',                    # Soil wetness in deep LSH/TOPMODEL layer (1)
        'sthzw_soilt',              # Soil wetness in deep LSH/TOPMODEL layer, on soil tiles (1)
        'sub_surf_roff',            # Gridbox sub-surface runoff (kg m-2 s-1)
        'surf_roff',                # Gridbox surface runoff (kg m-2 s-1)
        'sat_excess_roff_soilt',
        'tfall',                    # Gridbox throughfall (kg m-2 s-1)
        'zw',                       # Gridbox mean depth to water table (m)
        'zw_soilt',                 # Mean depth to water table on soil tiles (m)
        'fao_et0',                  # FAO Penman-Monteith evapotranspiration for reference crop (kg m-2 s-1)
        'ecan',                     # Tile evaporation from canopy/surface store for snow-free land tiles (kg m-2 s-1)
        'esoil',                    # Tile surface evapotranspiration from soil moisture store for snow-free land tiles (kg m-2 s-1)
        'et_stom',                  # Tile stomatal transpiration (kg m-2 s-1)
        'smcl',                     # Gridbox moisture content of each soil layer (kg m-2)
        'smcl_soilt',               # Moisture content of each soil layer on soil tiles (kg m-2)
        'soil_wet',                 # Gridbox total moisture content of each soil layer, as fraction of saturation (1)
        'sthu',                     # Gridbox unfrozen moisture content of each soil layer as a fraction of saturation (1)
        'sthu_soilt',               # Unfrozen moisture content of each soil layer as a fraction of saturation on soil tiles (1)
        'sthu_irr',                 # Gridbox wetness of each soil layer over irrigation (1)
        'sthu_irr_soilt',           # Wetness of each soil layer over irrigation on soil tiles (1)
        'irrig_water',              # Irrigation water demand (kg m-2 s-1)
        'snow_melt',
    ],
    'daily_vegetation' : [
        'lai'
    ]
}
