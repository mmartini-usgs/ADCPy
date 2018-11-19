# -*- coding: utf-8 -*-
"""
Processing for mooring 10811 TRDI workhorse V data
Created on Fri Sep 21 09:40:21 2018

@author: mmartini
"""
import netCDF4 as nc
from netCDF4 import num2date
import math
import sys
import datetime as dt
# this is important in order to import my package which is not on the python path
sys.path.append('c:\projects\python\ADCPy\EPICstuff')
import reshapeEPIC

# making the indeces
inpath = 'E:\\data\\Sandwich\\10811_V20784\\python\\'
outpath = 'E:\\data\\Sandwich\\10811_V20784\\python\\'

# continuous file
contFile = '10811whV.nc' # rotated file 
nfiles = 1
burstFileRoot = '10811whVcurrents.nc'
#burstFileRoot = '10811whVcurrentssmallb.nc'
#burstFileRoot = 'junkcurrents.nc'
indexFile = '10811whVcurrentindecesnc.txt'
#indexFile = 'junk.txt'

sample_rate = 2
burstlength = 150
burst_interval = 900 # 15 min interval
burst_start_offset = 0
# note we are dropping EPIC time as it is causing problems
vars2omit = {'EPIC_time','EPIC_time2'}
atts2omit = {'valid_range'} # this is in older converted files and needs to be removed
dryrun = False


# ------------------- the rest of this should be automatic, no user settings below
opstart = dt.datetime.now()
print('Start script run at ',opstart)
dim = 'time'

# ----------- execute
    
allslices = reshapeEPIC.generate_expected_start_times(inpath+contFile,dim,
            burst_start_offset, burst_interval, burstlength, sample_rate)

# here we limit the slices for testing
#print('** reducing the number of slices')
slices = allslices #[200:300]

contcdf= nc.Dataset(inpath+contFile,format="NETCDF4")

print('the last time is {} seconds from the start of the experiment'.format(contcdf['time'][-1]))
print('looking up the boundaries... this takes about 10 minutes on a 12 GB file')

edges = reshapeEPIC.find_boundaries(contcdf['time'][:], slices)
for x in edges[0:5]:  print('at indeces {} to {} we found times {} to {}'.format(x[0],x[1], 
    contcdf['time'][x[0]],contcdf['time'][x[1]]))
burstlengths = list(map(lambda t: t[1]-t[0], edges))
for x in burstlengths[0:5]:  print('bursts are {} long'.format(x))
    
contcdf.close()
print('elapsed time is {} min'.format((dt.datetime.now()-opstart).total_seconds()/60))

# TODO - this is not working when called from outside this file
reshapeEPIC.save_indexes_to_file(inpath+contFile, outpath+indexFile, edges)
    
nburstsperfile = int(math.floor(len(edges)/nfiles))
# now iterate throught he number of output files
#for ifile in range(1):
for ifile in range(nfiles):
    s = burstFileRoot.split('.')
    burstFile = s[0]+(f'%02d.' % ifile)+s[1]
    print('making burst file {}'.format(burstFile))
    
    startburst = ifile*nburstsperfile
    endburst = startburst+nburstsperfile
    
    if endburst > len(edges):
        enburst = len(edges)
        
    edges_this_file = edges[startburst:endburst]
    sampleseachburst = list(map(lambda t: t[1]-t[0],edges_this_file))
    
    # if there are no samples in a burst, we will skip the burst
    # skip them by removing them from this index list
    # this cleans up the tail end of the last file
    # TODO - use None to signal 
    idx_empty_bursts = list(map(lambda x: False if x==0 else True, sampleseachburst))
    print('Zeros samples in {} bursts, these will be ommitted'.format(idx_empty_bursts.count(0)))

    contcdf = nc.Dataset(inpath+contFile,format="NETCDF4")
    tunits = contcdf['time'].units
    x = list(range(5))+list(range(len(edges_this_file)-5,len(edges_this_file)))
    for i in x:
        print('burst {} will be {} samples from {} to {}'.format(
                i, sampleseachburst[i], 
                num2date(contcdf['time'][edges_this_file[i][0]],tunits),
                num2date(contcdf['time'][edges_this_file[i][1]],tunits)))

    contcdf.close()

    if not dryrun:
        reshapeEPIC.reshapeEPIC(inpath+contFile, outpath+burstFile, burstlength, 
                    dim='time', edges=edges_this_file, drop=vars2omit,
                    atts2drop=atts2omit)        
    

print('End script run at ',dt.datetime.now())
print('elapsed time is {} min'.format((dt.datetime.now()-opstart).total_seconds()/60))