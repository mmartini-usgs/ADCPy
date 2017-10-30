Example scripts to process multiple files from a Nortek Signature instrument

Used MIDAS.exe by Ocean Illumination to export from ad2cp to netCDF.  It takes two steps, first to the proprietary .ntk. then from ntk to netCDF.  Did not change the defaults.  It generates 24 files.  These are netCDF4 files with several groups, Config, Data and under data 4 data types, with different time bases and bin setups:  Burst (typical 4 beam Janues data), IBurstHR – the fifth beam, and two altimeter data sets.

Then, in python (saved in python2 subdirectory), working only with the Burst data to get profiles:
+ convert.py – script implements the following for the each of the MIDAS output netCDF files
    + doNortekRawFile – convert MIDAS netCDF output to USGS raw netCDF
    + doEPIC_ADCPfile – perform transformation to Earth and save as EPIC netCDF
+ combine.py – implements catEPIC to combine the individual .nc and .cdf files
+ resample.py – implement Xarray to reduce the data to hourly and 5 min averages
+ postproc.py – clean up from Xarray, reinstate EPIC convention items which are not compatible with or are lost in Xarray operations: time, time2, DELTA_T, etc.
+ unfill.py – clean up any NaNs

Known issues:
+ At 5 min intervals there are missing data and these are filled with very large values, code will be adjusted to prevent this.
+ Xarray introduces a time offset that needs to be adjusted for.
+ Incomplete bursts detected and confirmed in Contour – see email exchange below.
