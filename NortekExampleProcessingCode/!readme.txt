Example scripts to process multiple files from a Nortek Signature instrument

Convert.py 
- convert MIDAS.EXE netcdf output to USGS raw netcdf
- perform beam to earth rotations, add EPIC metadata

resample.py
- average data, hourly or 5 min, using xarray resample

postproc.py
- some EPIC compliance housekeeping

unfill.py
- remove _FillValue, especially NaNs from coordinate variables and dimensions, these being added by resample currently.  In work to figure out how to use resample to avoid this behavior