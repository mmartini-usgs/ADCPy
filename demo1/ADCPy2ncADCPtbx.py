# script to convert raw data to netcdf using python module
# you will want this installation of python:
# https://github.com/ioos/notebooks_demos/wiki/Installing-Conda-Python-with-the-IOOS-environment
# at the anaconda prompt in the data directory, with IOOS3 activated, run this script as
# ~\python>python ADCPyExample.py > ADCPyStatusOutput.txt

import sys
# if ADCPy is not on your python path, do this
sys.path.append('c:\projects\python\ADCPy')

import ADCPcdf2ncEPIC as cdf2nc
import datetime as dt
import numpy as np

ncFile = "py9991wh.nc"
cdfFile = "py9991wh.cdf"
attFile = "glob_att999.txt"
settings = {}
settings['good_ensembles'] = [0, np.inf]
#settings['good_ensembles'] = [300, np.inf] # one to be sure is in the water
settings['orientation'] = "UP" # uplooking ADCP, for downlookers, use DOWN
settings['transducer_offset_from_bottom'] = 1
settings['transformation'] = "EARTH" # | BEAM | INST

# here provide a dictionary of metadata needed for processing

# choice may be to trim by pressure where user gives the trim ensembles 
# or the minimum pressure

print('----')			
print('Start file conversion at', dt.datetime.now())
# when using this module this way, the start and end ensembles are required.
# use Inf to indicate all ensembles.
cdf2nc.doEPIC_ADCPfile(cdfFile, ncFile, attFile, settings)
print('Finished file conversion at', dt.datetime.now())
