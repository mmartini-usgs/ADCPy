# -*- coding: utf-8 -*-
"""
Conatiner for misc utilities
Created on Fri Oct 20 13:27:04 2017

catEPIC - concatenate like EPIC convention netCDF files, handle CF or EPIC time
        - files must be presented in the proper chronological order, this
          function will simply stick the data together in the order it is found

@author: mmartini USGS Woods Hole
"""

from netCDF4 import Dataset
from netCDF4 import num2date
import datetime as dt
import numpy as np
import math
import os
import sys
# this is important in order to import my package which is not on the python path
sys.path.append('c:\projects\python\ADCPy')
#from TRDIpd0tonetcdf import julian
from TRDIpd0tonetcdf import ajd
from ADCPcdf2ncEPIC import EPICtime2datetime


def catEPIC(datafiles, outfile):
    
    nfiles = len(datafiles)
    
    # use the first file to establish some information
    nc0 = Dataset(datafiles[0], mode = 'r', format = 'NETCDF4')
    varnames = nc0.variables.keys()
    if 'time2' not in varnames: 
        CFtime = True
        if 'calendar' not in nc0['time'].__dict__:
            print('No calendar specified, using gregorian')
            nccalendar = 'gregorian'
        else:
            nccalendar = nc0['time'].calendar
    else:
        CFtime = False
    
    nc0.close()
           
    # now glean time information from all the files
    alldt = np.array([])
    timelens = []
    for ifile in range(nfiles):
        print(datafiles[ifile])
        nc = Dataset(datafiles[ifile], mode = 'r', format = 'NETCDF4')
        timelens.append(nc.dimensions['time'].size)
        if CFtime: 
            tobj = num2date(nc['time'][:],nc['time'].units,calendar=nccalendar)
            alldt = np.concatenate((alldt,tobj))
        else:
            tobj = EPICtime2datetime(nc['time'][:],nc['time2'][:])
            alldt = np.concatenate((alldt,tobj))
        
        print('file %d is %s to %s' % (ifile, tobj[0], tobj[-1]))
        nc.close()
        
    # it was the case in the MATLAB version that the schema technique
    # would reorder the variables - not sure python does this
    # reordering the variables might not be a problem here since they are 
    # iterated by name
    
    # dimensions
    ncid_out = Dataset(outfile, "w", clobber=True, format="NETCDF4")
    ncid_in = Dataset(datafiles[0], mode = 'r', format = 'NETCDF4')
    for dimname in ncid_in.dimensions.keys():
        if 'time' in dimname:
            ncid_out.createDimension('time',len(alldt))
        else:
            ncid_out.createDimension(dimname,ncid_in.dimensions[dimname].size)
    
    # global attributes
    for attname in ncid_in.__dict__:
        ncid_out.setncattr(attname,ncid_in.getncattr(attname))
    
    # variables with their attributes
    for varname in ncid_in.variables.keys():
        ncid_out.createVariable(varname, ncid_in[varname].datatype, 
                            dimensions = ncid_in[varname].dimensions)
        for attname in ncid_in[varname].__dict__:
            ncid_out[varname].setncattr(attname, ncid_in[varname].getncattr(attname))
            
    ncid_out.close()
    ncid_in.close()
    
    # load the data
    ncid_out = Dataset(outfile, mode='r+', format="NETCDF4")
    print(timelens)
    for ifile in range(nfiles):
        if ifile == 0: 
            outidx_start = 0
            outidx_end = outidx_start+timelens[ifile]
        else:
            outidx_start = outidx_end
            outidx_end = outidx_start+timelens[ifile]
        print('getting data from file %s and writing %d to %d' % (datafiles[ifile],outidx_start,outidx_end))
        ncid_in = Dataset(datafiles[ifile], mode="r", format="NETCDF4")
        for varname in ncid_in.variables.keys():
            dimnames = ncid_in[varname].dimensions
            if 'time' in dimnames:
                s = outidx_start
                e = outidx_end
            else:
                s = 0
                e = len(ncid_in[varname])
            ndims = len(ncid_in[varname].dimensions)
            print('%s, %d' % (varname, ndims))
            print(len(ncid_in[varname]))
            if ndims == 1:
                ncid_out[varname][s:e] = ncid_in[varname][:]
            elif ndims == 2:
                ncid_out[varname][s:e,:] = ncid_in[varname][:,:]
            elif ndims == 3:
                ncid_out[varname][s:e,:,:] = ncid_in[varname][:,:,:]
            elif ndims == 4:
                ncid_out[varname][s:e,:,:,:] = ncid_in[varname][:,:,:,:]                        
        ncid_in.close()
    
    # finally, put the correct time span in th output file
    units = "seconds since %d-%d-%d %d:%d:%f 0:00" % (alldt[0].year, 
        alldt[0].month,alldt[0].day,alldt[0].hour,alldt[0].minute,
        alldt[0].second+alldt[0].microsecond/1000000)
    elapsed = alldt-alldt[0] # output is a numpy array of timedeltas
    # have to iterate to get at the timedelta objects in the numpy container
    # seems crazy, someone please show me a better trick!
    elapsed_sec = []
    for e in elapsed:  elapsed_sec.append(e.total_seconds())
    t = np.zeros((len(alldt),1))
    t2 = np.zeros((len(alldt),1))
    for i in range(len(alldt)): 
        jd = ajd(alldt[i])
        t[i] = int(math.floor(jd))
        t2[i] = int((jd - math.floor(jd))*(24*3600*1000))
    if CFtime:
        ncid_out['time'][:] = elapsed_sec[:]
        ncid_out['time'].units = units
        ncid_out['EPIC_time'][:] = t[:]
        ncid_out['EPIC_time2'][:] = t2[:]
    else:
        ncid_out['CFtime'][:] = elapsed_sec[:]
        ncid_out['CFtime'].units = units
        ncid_out['time'][:] = int(math.floor(jd))
        ncid_out['time2'][:] = int((jd - math.floor(jd))*(24*3600*1000))
    
    # recompute start_time and end_time
    ncid_out.start_time = '%s' % num2date(ncid_out['time'][0],ncid_out['time'].units)
    print(ncid_out.start_time)
    ncid_out.stop_time = '%s' % num2date(ncid_out['time'][-1],ncid_out['time'].units)
    print(ncid_out.stop_time)
    
    # TODO update history
    
    ncid_out.close()
            
    
def __main():
    print('%s running on python %s' % (sys.argv[0], sys.version))
	
    if len(sys.argv) < 2:
        print("%s useage:" % sys.argv[0])
        print("catEPIC file_list outfile\n" )
        sys.exit(1)
        
    try:
        datafiles = sys.argv[1]
    except:
        print('error - input file list missing')
        sys.exit(1)
        
    try:
        outfile = sys.argv[2]
    except:
        print('error - output file name missing')
        sys.exit(1)
        
    # some input testing
    if ~os.path.isfile(datafiles[0]):
        print('error - input file not found')
        sys.exit(1)
        
    if ~os.path.isfile(outfile):
        print('%s will be overwritten' % outfile)
    
    print('Start concatination at ',dt.datetime.now())
    
    catEPIC(datafiles, outfile)
    
    print('Finished file concatination at ',dt.datetime.now())
 
if __name__ == "__main__":
    __main()        
    