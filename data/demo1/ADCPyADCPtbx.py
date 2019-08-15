# script to convert raw data to netcdf using python module
# you will want this installation of python:
# https://github.com/ioos/notebooks_demos/wiki/Installing-Conda-Python-with-the-IOOS-environment
# at the anaconda prompt in the data directory, with IOOS3 activated, run this script as
# ~\python>python ADCPyExample.py > ADCPyStatusOutput.txt
# tested with MMfavs environment on 5/24/2019

# import sys
# # if ADCPy is not on your python path, do this
# sys.path.append("c:\\projects\\python\\ADCPy")

import TRDIstuff.TRDIpd0tonetcdf as pd0
import datetime as dt
import numpy as np

pd0File = "9991wh000.000"
serialnum = "473"
delta_t_to_use = "900"  # seconds, as a string
# timetype = "CF"  # the time variable will have a CF time format
timetype = "EPIC"  # there will be time and time2, and cf_time variables

if timetype == "CF":
    cdfFile = "py9991whCFtime.cdf"
else:
    cdfFile = "py9991whEPICtime.cdf"

goodens = [265, np.inf]
# goodens = [0, np.inf]

maxens, ensLen, ensData, startOfData = pd0.analyzepd0file(pd0File)
# will print out two layers of nested dictionary
print("data found in file ", pd0File)
for key, value in sorted(ensData.items()):
    vtype = type(ensData[key])
    if vtype == dict:
        print("\n", key, ":")
        for key1, value1 in sorted(ensData[key].items()):
            print(key1, ":", value1)
        else:
            print(key, ":", value)

print("----")
print("Start file conversion at", dt.datetime.now())
# when using this module this way, the start and end ensembles are required.
# use Inf to indicate all ensembles.
ensCount, cdfIdx, ensError = pd0.convert_pd0_to_netcdf(pd0File, cdfFile, goodens, serialnum, timetype, delta_t_to_use)
print("Finished file conversion at", dt.datetime.now())
print("Ensemble error = {}".format(ensError))
