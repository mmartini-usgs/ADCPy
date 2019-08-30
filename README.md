# ADCPy - code to work with ADCP data from the raw binary using python 3.x

[![Documentation Status](https://readthedocs.org/projects/adcpy/badge/?version=latest)](https://adcpy.readthedocs.io/en/latest/?badge=latest)
[![Build Status](https://travis-ci.org/mmartini-usgs/ADCPy.svg?branch=master)](https://travis-ci.org/mmartini-usgs/ADCPy)
### Purpose
This code prepares large amounts of single ping ADCP data from the raw binary for use with xarray by converting it to netCDF.  

### Motivation
The code was written for the TRDI ADCP when I discovered theat TRDI's Velocity software could not easily export single ping data.  While there are other packages out there, as the time of writing this code, I had yet to find one that saved the data in netCDF format (so it can be accessed with xarray and dask), could be run on linux, windows and mac, and did not load it into memory (the files I have are > 2GB)

The code is written as a module of functions, rather than classes, ensemble information is stored as nested dicts, in order to be more readable and to make the structure of the raw data (particularly the TRDI instruments) understandable.

### Status
As the code stands now, a 3.5 GB, single ping Workhorse ADCP .pd0 file with 3 Million ensembles will take 4-5 hours to convert.  I live with this, because I can just let the conversion happen overnight on such large data sets, and once my data is in netCDF, everything else is convenient and fast.  I suspect that more speed might be acheived by making use of xarray and dask to write the netCDF output, and I may do this if time allows, and I invite an enterprising soul to beat me to it.  I use this code myself on a routine basis in my work, and continue to make it better as I learn more about python.

At USGS Coastal and Marine Geology we use the PMEL EPIC convention for netCDF as we started doing this back in the early 1990's.  Downstream we do convert to more current CF conventions, however our diagnostic and other legacy code for processing instrument data from binary and other raw formats depends on the EPIC convention for time, so you will see a time (Time (UTC) in True Julian Days: 2440000 = 0000 h on May 23, 1968) and time2 (msec since 0:00 GMT) variable created as default.  This may confuse your code.  If you want the more python friendly CF time (seconds since 1970-01-01T00:00:00 UTC) set timetype to CF.

Use at  your own risk - this is a work in progress and a python learning project.

Enjoy,

Marinna
