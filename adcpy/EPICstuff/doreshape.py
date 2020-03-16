# -*- coding: utf-8 -*-
"""
Test file for reshapeEPIC
Created on Fri Sep 21 09:40:21 2018

@author: mmartini
"""
import netCDF4 as nc
from netCDF4 import num2date
import math
import sys
import datetime as dt
# this is important in order to import my package which is not on the python path
sys.path.append(r'c:\projects\python\ADCPy\EPICstuff')
from EPICstuff.reshapeEPIC import reshapeEPIC

# making the indeces
#contFile = 'E:\\data\\Sandwich\\10811_V20784\\python\\10811whVsubset00nbetterfill.cdf'
#burstFile = 'E:\\data\\Sandwich\\10811_V20784\\python\\junkbetterfill.cdf'
contFile = 'E:\\data\\Sandwich\\10811_V20784\\python\\10811whVsubset00n.cdf'
burstFile = 'E:\\data\\Sandwich\\10811_V20784\\python\\junk.cdf'
sample_rate = 2
burstlength = 2048
burst_interval = 3600
burst_start_offset = 0
vars2omit = {'HdgSTD','PtchSTD','RollSTD','S','xmitc','Ambient_Temp','Pressure+','Pressure-',
        'Attitude_Temp','EWD1','EWD2','EWD3','EWD4','PressVar'}
dryrun = False
nfiles = 1

# ------------------- the rest of this should be automatic, no user settings below
opstart = dt.datetime.now()
print('Start script run at ',opstart)
dim = 'time'

# make a function to identify the indeces
# and this is where tuples are nice
def find_boundaries(data,edges):
    # data is a list object of time stamps from the data
    # edges is a list of tuples of start and end times
    idx = []
    for edge in edges:
        s = None
        for t in enumerate(data):
            if t[1] >= edge[0]:
                s = t[0]
                break

        e = None
        for t in enumerate(data):
            if t[1] >= edge[1]:
                e = t[0]
                break

        # print((s,e))
        idx.append((s,e))
    
    return idx


contcdf = nc.Dataset(contFile, format="NETCDF4")

if dim in contcdf.dimensions:
    print('the dimension we are operating on is {}'.format(dim))
else:
    print('{} not found in {} file, aborting'.format(dim, contFile))
    contcdf.close()
    
# get the number of bursts based ont he elapsed time
tfirst = num2date(contcdf['time'][1], contcdf['time'].units)
tlast = num2date(contcdf['time'][-1], contcdf['time'].units)
nsec = (tlast-tfirst).total_seconds()
nbursts = int(nsec / burst_interval)
burst_start_times = []
for x in range(nbursts): 
    burst_start_times.append(burst_start_offset+x*burst_interval)
    
burst_duration = (1 / sample_rate) * burstlength # seconds

burst_end_times = list(map(lambda x: x+burst_duration,burst_start_times))

print('start times {} such as {}...'.format(len(burst_start_times),burst_start_times[0:5]))
print('end times {} such as {}...'.format(len(burst_end_times),burst_end_times[0:5]))

# it turns out later to be convenient to have this as a list of tuples of start and end
slices = list(map(lambda s,e: (s,e), burst_start_times, burst_end_times))
print('edge tuples {} such as {}...'.format(len(slices),slices[0:5]))

burstnum = 0
print('the last time is {} seconds from the start of the experiment'.format(contcdf['time'][-1]))
print('looking up the boundaries')
edges = find_boundaries(contcdf['time'][:], slices)
for x in edges[0:5]:  print('at indeces {} to {} we found times {} to {}'.format(x[0], x[1],
    contcdf['time'][x[0]],contcdf['time'][x[1]]))
# later we'll figure out how to do this as a list comprehension, for now it works
burstlengths = list(map(lambda t: t[1]-t[0], edges))
for x in burstlengths[0:5]:
    print('bursts are {} long'.format(x))
    
nburstsperfile = int(math.floor(nbursts/nfiles))

contcdf.close()

if not dryrun:
    print('Starting the reshape')
    reshapeEPIC(contFile, burstFile, burstlength, dim='time', edges=edges, drop=vars2omit)
    
print('Finished script run at ', dt.datetime.now())
print('elapsed time is {} min'.format((dt.datetime.now()-opstart).total_seconds()/60))
