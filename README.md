# ADCPy - code to work with ADCP data from the raw binary using python 3.x

[![Documentation Status](https://readthedocs.org/projects/adcpy/badge/?version=latest)](https://adcpy.readthedocs.io/en/latest/?badge=latest)
[![Build Status](https://travis-ci.org/mmartini-usgs/ADCPy.svg?branch=master)](https://travis-ci.org/mmartini-usgs/ADCPy)

The purpose of this code is to prepare large amounts of ADCP data from the raw binary for use with xarray by converting it to netCDF.  

Written by Marinna Martini, 1/12/2017

Use at  your own risk - this is a work in progress and a python learning project.

As the code stands now, a 3.5 GB, single ping Workhorse ADCP .pd0 file with 3 Million ensembles will take 4-5 hours to convert.  I live with this, because I can just let the conversion happen overnight on such large data sets, and once my data is in netCDF, everything else is convenient and fast.  I suspect that more speed might be acheived by making use of xarray and dask to write the netCDF output, and I may do this if time allows, and I invite an enterprising soul to beat me to it.

The code is written as a module of functions, rather than classes, in order to be more readable and to make the structure of the raw data (particularly the TRDI instruments) understandable.

At USGS Coastal and Marine Geology we use the PMEL EPIC convention for netCDF as we started doing this back in the early 1990's.  Downstream we do convert to more current CF conventions, however our diagnostic and other legacy code for processing instrument data from binary and other raw formats depends on the EPIC convention for time, so you will see a time (Time (UTC) in True Julian Days: 2440000 = 0000 h on May 23, 1968) and time2 (msec since 0:00 GMT) variable created as default.  This may confuse your code.  If you want the more python friendly CF time (seconds since 1970-01-01T00:00:00 UTC) set timetype to CF.
