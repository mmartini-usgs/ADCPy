# script to process raw data to netcdf using python module
# at the anaconda prompt in the data directory, with IOOS3 activated, runt his script as
# E:\data\MVCO14\101003_ADCP767\python>python do767py.py > output.txt

# TODO - known issues
# NotImplementedError: skipna=True not yet implemented for mean with dtype uint32

import sys
# this is important in order to import my package which is not on the python path
sys.path.append('c:\projects\python\ADCPy')
#import pd0splitter
#import os
import datetime as dt
import xarray as xr
import netCDF4 as nc
from EPICstuff import resample_cleanup

datapath = "E:\\data\\Matanzas\\WellTest2017\\Signature\\python2\\"

datafiles = ['1108sigall.nc']
doHourly = True
do5min = True

# ------------------ beyond here the user should not need to edit anything
allstarttime = dt.datetime.now()

for filenum in range(len(datafiles)):

    ncFile = datapath + datafiles[filenum]
    
    allstarttime = dt.datetime.now()
    
    infile = ncFile
    # first lets see what kind of time we have for the main axis
    # if xarray - MUST set decode_time=False to see the corrdinate time variable units
    ds = xr.open_dataset(infile, decode_times=False)
    timeunits = ds['time'].attrs['units']
    ds.close()
    print(timeunits)
    if 'True Julian Day' in timeunits:
        timetype = 'EPIC'
    else:
        timetype = 'CF'
    
    print('in %s we have %s time, xarray operations will work' % (ncFile, timetype))
    
    # now that we know what time we have, we can allow timestamps to be interpreted here
    # time is expressed as decimal seconds since the first ping in the data set.  
    # Acoustic current profiler time stamps may have resolution of as little as 16 Hz, depending on the instrument
    # these netCDF files have both cf_time and time and time2.  To prevent xarray from trying to use time & time2 
    # as the time axis for this dataset, we ask xarray to omit them
    # do NOT use xarray.load_dataset for these files, they are very large
    vars2omit = {'EPIC_time','EPIC_time2'}
    ds = xr.open_dataset(infile,decode_times=True,drop_variables=vars2omit,mask_and_scale=False)
    
    if doHourly:
        print('Before resampling')
        print(ds.dims)
        #print(ds['lat'])
        print(ds['lat'].__dict__)
        starttime = dt.datetime.now()
        print("start resample at %s" % starttime)
        #ds_1h = ds.resample('H','time',how='mean',keep_attrs=True,skipna=True)
        ds_1h = ds.resample('H','time',how='mean',closed='left',keep_attrs=True)
        endtime = dt.datetime.now()
        print("finished resample at %s" % endtime)
        print("processing time was %s" % (endtime-starttime))
        print(ds_1h.dims)
        
        fileparts = infile.split('.')
        outfile = fileparts[0]+"_1h.nc"
        print(outfile)
        
        print('after resampling')
        print(ds_1h['lat'].__dict__)
        ds_1h.to_netcdf(outfile)
        
        resample_cleanup([outfile])
    
    if do5min:
        print(ds.dims)
        starttime = dt.datetime.now()
        print("start resample at %s" % starttime)
        ds_5min = ds.resample('5min','time',how='mean',closed='left',keep_attrs=True)
        endtime = dt.datetime.now()
        print("finished resample at %s" % endtime)
        print("processing time was %s" % (endtime-starttime))
        print(ds_5min.dims)
        
        fileparts = infile.split('.')
        outfile = fileparts[0]+"_5m.nc"
        print(outfile)
        ds_5min.to_netcdf(outfile)

        resample_cleanup([outfile])
    
    ds.close()
    
    # check _FillValue
    print('as netcdf Dataset')
    ncdata = nc.Dataset(outfile, mode="r", format='NETCDF4')
    print(ncdata['lat'].__dict__)
    ncdata.close()
    
allendtime = dt.datetime.now()
print("For all the operations, processing time was %s" % (allendtime-allstarttime))