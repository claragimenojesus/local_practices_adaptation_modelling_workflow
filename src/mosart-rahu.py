 #!/usr/bin/env python
# coding: utf-8

from mosartwmpy import Model
from mosartwmpy.utilities.download_data import download_data
from datetime import date, datetime
import xarray
import numpy as np

def main():
    mosart_wm = Model()
    mosart_wm.initialize('./config_rahu.yaml') # change path appropriately
    # mosart_wm.config["simulation.end_date"] = date(1981, 5, 14)
    mosart_wm.update_until(mosart_wm.get_end_time())

if __name__ == '__main__':
    main()


