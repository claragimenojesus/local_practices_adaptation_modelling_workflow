#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import xarray as xr
import sys

OUTDIR = str(os.environ['MOSART_OUTPUT'])
NBS = str(os.environ['NBS'])
SCENARIO_SH = str(os.environ['SCENARIO_SH'])
SC = str(os.environ['SC'])

OUTLET_X_km105=-72.5336
OUTLET_Y_km105=-13.1831

OUTLET_X_Pisac = -71.8545
OUTLET_Y_Pisac = -13.4225


def main():
    #OUTDIR_S = sys.argv[1]i
    if NBS!="":
        fs = [f for f in os.listdir(os.path.join(OUTDIR,NBS,SC+"_"+NBS)) if re.match(r'.*_[0-9]{4}_[0-9]{2}\.nc', f)]
        fs.sort()
        fs = [os.path.join(OUTDIR,NBS,SC+"_"+NBS, f) for f in fs]
        x = xr.open_mfdataset(fs)

        xi_km105 = x.sel(lat=OUTLET_Y_km105, lon=OUTLET_X_km105, method='nearest')
        xi_km105 = xi_km105.to_dataframe()
        xi_km105.to_csv(os.path.join(OUTDIR,NBS,SC+"_"+NBS, 'output_pt_km105.csv'))

        xi_pisac = x.sel(lat=OUTLET_Y_Pisac, lon=OUTLET_X_Pisac, method='nearest')
        xi_pisac = xi_pisac.to_dataframe()
        xi_pisac.to_csv(os.path.join(OUTDIR, NBS,SC+"_"+NBS, 'output_pt_Pisac.csv'))

    else:
        fs = [f for f in os.listdir(os.path.join(OUTDIR,NBS)) if re.match(r'.*_[0-9]{4}_[0-9]{2}\.nc', f)]
        fs.sort()
        fs = [os.path.join(OUTDIR,NBS, f) for f in fs]
        x = xr.open_mfdataset(fs)

        xi_km105 = x.sel(lat=OUTLET_Y_km105, lon=OUTLET_X_km105, method='nearest')
        xi_km105 = xi_km105.to_dataframe()
        xi_km105.to_csv(os.path.join(OUTDIR,NBS, 'output_pt_km105.csv'))

        xi_pisac = x.sel(lat=OUTLET_Y_Pisac, lon=OUTLET_X_Pisac, method='nearest')
        xi_pisac = xi_pisac.to_dataframe()
        xi_pisac.to_csv(os.path.join(OUTDIR, NBS, 'output_pt_Pisac.csv'))


if __name__ == '__main__':
    main()
