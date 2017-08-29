# ADCPy

code to work with ADCP data from the raw binary using python 3.x

Written by Marinna Martini, 1/12/2017

I used this guide to set up this project:
http://docs.python-guide.org/en/latest/writing/structure/

Use at  your own risk - this is a work in progress and a python learning project.

Note that a 3.5 GB, single ping Workhorse ADCP .pd0 file with 3 Million ensembles will take 4-5 hours to convert.

you will want this installation of python:
https://github.com/ioos/notebooks_demos/wiki/Installing-Conda-Python-with-the-IOOS-environment

At USGS Coastal and Marine we use the PMEL EPIC convention for netCDF as we started doing this back in the early 1990's.  Downstream we do convert to more current CF conventions, however our diagnostic and other legacy code for processing instrument data from binary and other raw formats depends on the EPIC convention for time, so you will see a time (Time (UTC) in True Julian Days: 2440000 = 0000 h on May 23, 1968) and time2 (msec since 0:00 GMT) variable created as default.  This may confuse your code.  If you want the more python friendly CF time (seconds since 1970-01-01T00:00:00 UTC) set timetype to CF.
