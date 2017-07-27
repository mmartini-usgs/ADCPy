# ADCPy

code to work with ADCP data from the raw binary using python 3.x

In this directory is data to practice or test your installation using data from a TRDI ADCP.
FOr users of the old USGS ADCP toolbox for MATLAB, this is the same as "demo" in that toolbox.

The raw binary ADCP file is 9991wh000.000
The text file with standard USGS metadata is glob_att999.txt
	This file is read and the contents added as global attributes tot he netCDF file
	Some of these attributes (like water_depth) are necessary for processing and converting the data
	
Data are processed in several steps
	1. Conversion from raw binary to netCDF is done by ADCPyADCPtbx.py which calls TRDIpd0tonetcdf.py, and this file is exclusive to TRDI data.  Data are output in raw form, as an EPIC compliant netCDF file in py9991wh.cdf
	2. Rotation from beam or instrument coordinates to East, North and Up is done by ADCPy2ncADCPtbx.py which calls ADCPcdf2ncEPIC.py.  Data are output to a releaseable USGS netCDF file as py9991wh.nc, however these data are not cleaned or edited in any way to trim bins out of the water or noisy data.

Written by Marinna Martini, 1/12/2017
