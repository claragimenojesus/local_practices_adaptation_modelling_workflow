#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import xarray as xr
import sys

#OUTDIR = '/rds/general/user/cg2117/ephemeral/mosart-output/u-cz655_final'
OUTDIR = str(os.environ['MOSART_OUTPUT'])
NBS = str(os.environ['NBS'])
SCENARIO_SH = str(os.environ['SCENARIO_SH'])
SC = str(os.environ['SC'])
#SC = re.findall('5/(.*)',SCENARIO_SH)[0]

OUTLET_X_km105=-72.5336
OUTLET_Y_km105=-13.1831

#OUTLET_X_km105_sen=-72.564
#OUTLET_Y_km105_sen=-13.174

OUTLET_X_Chilca=-72.3406
OUTLET_Y_Chilca=-13.2236

OUTLET_X_Pisac = -71.8545
OUTLET_Y_Pisac = -13.4225

OUTLET_X_Salcca = -71.2319
OUTLET_Y_Salcca = -14.1699
# Waldo coordinates
# km105 x=-72.5336 y=-13.1831
# Chilca x=-72.3406 y=-13.2236
# Pisac x=-71.8545 y=-13.4225
# Salcca x=-71.2319 y=-14.1699

# Jose calculated outlet -72.62478 -13.00928
#OUTLET_X = -72.533592
#OUTLET_Y = -13.183103

# Use Senamhi coordinates, ANA coordinates wrong
# km105 x=-72.564 y=-13.174
# Paucartambo x=-71.597 y=-13.317 (out of catchment boundaries)
# Salcca x=-71.232 y=-14.170
#OUTLET_X = -71.841
#OUTLET_Y = -13.428
# Pisac x=-71.841 y=-13.428
# Chilca x=-72.341 y=-13.221
# Chacllabamba x=-71.720295 y=-13.106774 (out of catchment boundaries)


def main():
    #OUTDIR_S = sys.argv[1]
    fs = [f for f in os.listdir(os.path.join(OUTDIR,NBS,SC+"_"+NBS)) if re.match(r'.*_[0-9]{4}_[0-9]{2}\.nc', f)]
    fs.sort()
    fs = [os.path.join(OUTDIR,NBS,SC+"_"+NBS, f) for f in fs]
    x = xr.open_mfdataset(fs)

    xi_km105 = x.sel(lat=OUTLET_Y_km105, lon=OUTLET_X_km105, method='nearest')
    xi_km105 = xi_km105.to_dataframe()
    xi_km105.to_csv(os.path.join(OUTDIR,NBS,SC+"_"+NBS, 'output_pt_km105.csv'))

#    xi_km105_sen = x.sel(lat=OUTLET_Y_km105_sen, lon=OUTLET_X_km105_sen, method='nearest')
#    xi_km105_sen = xi_km105_sen.to_dataframe()
#    xi_km105_sen.to_csv(os.path.join(OUTDIR, 'output_pt_km105_sen.csv'))

    xi_chilca = x.sel(lat=OUTLET_Y_Chilca, lon=OUTLET_X_Chilca, method='nearest')
    xi_chilca = xi_chilca.to_dataframe()
    xi_chilca.to_csv(os.path.join(OUTDIR, NBS,SC+"_"+NBS, 'output_pt_Chilca.csv'))

    xi_pisac = x.sel(lat=OUTLET_Y_Pisac, lon=OUTLET_X_Pisac, method='nearest')
    xi_pisac = xi_pisac.to_dataframe()
    xi_pisac.to_csv(os.path.join(OUTDIR, NBS,SC+"_"+NBS, 'output_pt_Pisac.csv'))

    xi_salcca = x.sel(lat=OUTLET_Y_Salcca, lon=OUTLET_X_Salcca, method='nearest')
    xi_salcca = xi_salcca.to_dataframe()
    xi_salcca.to_csv(os.path.join(OUTDIR, NBS,SC+"_"+NBS, 'output_pt_Salcca.csv'))


if __name__ == '__main__':
    main()
