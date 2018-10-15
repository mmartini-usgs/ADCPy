# -*- coding: utf-8 -*-
"""
function issueflags = reshapeEPIC(contFile, burstFile, burstlength, ...
    dim=None, edges=None, drop=None)

apportion a continuous time series file into bursts (e.g. reshape)

contFile = netcdf file with continuous data
burstFile = file to store the reshaped data, attributes will be copied
burstlength = maximum number of samples in each burst
dim = dimension along which we will split the data, usually 'time' or 'Rec'
edges = list of tuples (start, end) of edges defining the edges of each burst
drop = set of variable names to omit from the output file

the expected dimensions are [time, depth, lat, lon] for continuous EPIC
we are reshaping to [time, sample, depth, lat, lon]
for ADCP files, beams are expressed as variables vel1, vel2, ... veln

NOTE:  if there is a shape problem then the data file might not be
properly understood.  It might be that this code won't work
this problem will become evident if an error is produced and operation
returns to the keyboard in debug mode.  If this happens, check the shapes
of the variables.

WARNING:  time may not be monotomically increasing within bursts (e.g.
along the sample dimension)
this means that if the number of samples per burst is inconsistent, or if
there are gaps in time, the end of a burst may be fill_value, including
time values

Marinna Martini for the USGS in Woods Hole, 9/20/2018
originally coded for MATLAB as reshapeEPIC
https://cmgsoft.repositoryhosting.com/trac/cmgsoft_m-cmg/browser/trunk/MMstuff/reshapeEPIC.m

Created on Thu Sep 20 14:52:42 2018

@author: mmartini
https://github.com/mmartini-usgs
"""

import sys
import os
import datetime as dt
import netCDF4 as nc
from netCDF4 import num2date
import numpy as np
import math

def reshapeEPIC(*args, **kwargs):
    # the argument passing here works fine
    print('%s running on python %s' % (sys.argv[0], sys.version))
    print('Start file conversion at ',dt.datetime.now())
    
#    for s in args:
#        print(s)
#    for k in kwargs.keys():
#        print('{} = {}\n'.format(k,kwargs[k]))
        
    contFile = args[0]
    burstFile = args[1]
    burstlength = args[2]        
    if 'dim' in kwargs.keys():  
        dim = kwargs['dim']    
    else:
        dim = 'time'
    if 'edges' in kwargs.keys():  
        edges = kwargs['edges']   
    else:
        # make edges based on burstlength later, 
        # when we know the total number of samples in the file
        edges = None
    if 'drop' in kwargs.keys():  
        drop = kwargs['drop']  
    else:
        drop = None
    if 'verbose' in kwargs.keys():  
        verbose = kwargs['verbose']  
    else:
        verbose = False
    if 'atts2drop' in kwargs.keys():
        atts2drop = kwargs['atts2drop']
    else:
        atts2drop = {}
    
    print('Start file conversion at ',dt.datetime.now())

    # check for the output file's existence before we try to delete it.
    try:
        os.remove(burstFile)
        print('{} removed'.format(burstFile))
    except FileNotFoundError:
        pass

    contcdf = nc.Dataset(contFile,format="NETCDF4")
    if dim in contcdf.dimensions:
        print('the dimension we are operating on is {}'.format(dim))
    else:
        print('{} not found in input file, aborting'.format(dim))
        contcdf.close()
    
    # create the new file
    burstcdf = nc.Dataset(burstFile, mode="w", clobber=True, format='NETCDF4')
    # incoming data may be uneven, we need proper fill to be added
    burstcdf.set_fill_on()
    
    # copy the global attributes
    # first get a dict of them so that we can iterate
    gatts = {}
    for attr in contcdf.ncattrs():
        #print('{} = {}'.format(attr,getattr(contcdf,attr)))
        gatts[attr] = getattr(contcdf,attr)
        
    # add a few more important ones we will fill in later
    gatts['start_time'] = ""
    gatts['stop_time'] = ""
    gatts['DELTA_T'] = ""
    gatts['history'] = getattr(contcdf,'history')+'; converted to bursts by reshapeEPIC.py'
    
    burstcdf.setncatts(gatts)

    print('Finished copying global attributes\n')
    
    for item in contcdf.dimensions.items():
        print('Defining dimension {} which is {} long in continuous file'.format(item[0],len(item[1])))
        if item[0] == dim: 
            # this is the dimension along which we will reshape
            if len(edges) > 1:
                nbursts = len(edges)
            else:
                nbursts = math.floor(len(item[1])/burstlength)
            burstcdf.createDimension(dim,nbursts)
            print('Reshaped dimension {} created for {} bursts'.format(item[0],nbursts))
        else:
            burstcdf.createDimension(item[0],len(item[1]))
            
    burstcdf.createDimension('sample',burstlength)
     
    # ---------------- set up the variables
    # order of dimensions matters.
    # per https://cmgsoft.repositoryhosting.com/trac/cmgsoft_m-cmg/wiki/EPIC_compliant
    # for a burst file dimension order needs to be time, sample, depth, [lat], [lon]
    # TODO - check the above order and reorder dimensions if necessary
    for cvar in contcdf.variables.items():
        cvarobj = cvar[1]
        print('{} is data type {}'.format(cvarobj.name,cvarobj.dtype))
        try:
            fillValue = cvarobj.getncattr('_FillValue')
            if verbose:
                print('\tthe fill value is {}'.format(fillValue))
        except AttributeError:
            print('\tfailed to read the fill value')
            fillValue = False # do not use None here!!!
                
        if verbose:
            print('\tfillValue in burst file will be set to {} (if None, then False will be used)'.format(fillValue))
        
        if cvarobj.name not in drop:  # are we copying this variable?
            dtype = cvarobj.dtype
            
            if dim in cvarobj.dimensions:  # are we reshaping this variable?
                vdims_cont = cvarobj.dimensions
                vdims_burst = []
                for t in enumerate(vdims_cont):
                    vdims_burst.append(t[1])
    
                    if t[1] == dim:
                        vdims_burst.append('sample')
                        print('\tappending sample in {}'.format(cvarobj.name))
    
                varobj = burstcdf.createVariable(cvarobj.name, dtype, tuple(vdims_burst), 
                                        fill_value=fillValue)
            else:
                # for a normal copy, no reshape
                varobj = burstcdf.createVariable(cvarobj.name, dtype, cvarobj.dimensions, 
                                        fill_value=fillValue)
                
            # copy the variable attributes
            # first get a dict of them so that we can iterate
            vatts = {}
            for attr in cvarobj.ncattrs():
                #print('{} = {}'.format(attr,getattr(contcdf,attr)))
                if attr not in atts2drop:
                    vatts[attr] = getattr(cvarobj,attr)
    
            try:
                varobj.setncatts(vatts)
            except AttributeError:
                print('AttributeError for {}'.format(cvarobj.name))

    # not a coordinate but a fill value of None might cause problems
    burstcdf.createVariable('burst', 'uint16', ('time'), fill_value=False)
    # these are coordinates and thus cannot have fill as their values
    burstcdf.createVariable('sample', 'uint16', ('sample'), fill_value=False)
    try:
        burstcdf.createVariable('depth', 'float32', ('depth'), fill_value=False)
    except:
        pass # likely depth was already set up if this happens               
        
    # TODO - if edges is length 1, then we need to create the burst edges here
    
    # --------- populate the file
    # note that we don't have to change from a define to a read mode here
    
    # coordinate variables are small(er) and can be done at once, be sure to use
    # generative methods
    print('\nNow populating data for {} bursts'.format(nbursts))
    burstcdf['burst'][:] = list(range(nbursts))
    burstcdf['sample'][:] = list(range(burstlength))

    nbins = len(burstcdf['depth'])
    try: 
        binsize = contcdf['depth'].bin_size_m
    except AttributeError:
        try: 
            binsize = contcdf['depth'].bin_size
        except AttributeError: 
            # TODO - get user input here
            print('Warning:  no depth size information found, assuming 1 m')
            binsize = 1
    
    try: 
        bin1distance = contcdf['depth'].center_first_bin_m
    except AttributeError:
        try: 
            bin1distance = contcdf['depth'].center_first_bin
        except AttributeError: 
            # TODO - get user input here
            print('Warning:  no depth center of first bin information found, assuming 0.5 bins ')
            bin1distance = binsize/2            
    
    ranges_m = list(map(lambda ibin: bin1distance/100+ibin*binsize/100, range(nbins)))
    burstcdf['depth'][:] = ranges_m
       
    issueflags = {}
    diagnosticvars = {} # vars to generate diagnostic output
    for cvar in contcdf.variables.items():
        varname = cvar[1].name
        issueflags[varname] = []
        if varname not in drop:
            # the variable objects in Continuous and Burst files
            cvarobj = contcdf[varname]
            bvarobj = burstcdf[varname]
        
            vdims_cont = cvarobj.dimensions
            vshapes_cont = cvarobj.shape
            vndims_cont = len(cvarobj.dimensions)
            vdims_burst = burstcdf[varname].dimensions
            vshapes_burst = burstcdf[varname].shape
            vndims_burst = len(burstcdf[varname].dimensions)
            if verbose:
                print('{}\tin Continuous file is data type {} shape {}'.format(
                        varname,cvarobj.dtype,cvarobj.shape))
                print('\tin Burst file it is data type {} shape {}'.format(
                        bvarobj.dtype,bvarobj.shape))

            try:
                fillval_burst = burstcdf[varname].getncattr('_FillValue')
            except:
                if ('EPIC' in varname) and verbose:
                    # EPIC was ending up with odd fill values in the raw file
                    # this will avoid the typerror when EPIC_time is written 
                    # not sure it is the best solution, for now it works
                    fillval_burst = 0
                    print('\tfillval_burst {}'.format(fillval_burst))
                    # this will prevent EPIC_time from being written
                    #fillval_burst = None
                else:
                    fillval_burst = None
    
            
            if 'sample' not in vdims_burst:
                bvarobj[:] = contcdf[varname][:]
            else:            
                for iburst in range(nbursts):
                    contcorner = np.zeros(vndims_cont)
                    contedges = np.ones(vndims_cont)
                    # look up data in the continuous file according to the user's indeces
                    contcorner[vdims_cont.index('time')] = edges[iburst][0]
                    ndatasamples = edges[iburst][1]-edges[iburst][0]
                    contedges[vdims_cont.index('time')] = ndatasamples
                    
                    if 'depth' in vdims_cont:
                        contedges[vdims_cont.index('depth')] = vshapes_cont[vdims_cont.index('depth')]
                    
                    if (iburst == 0) and verbose:
                        print('\tconcorner = {}, contedges = {}'.format(
                                contcorner, contedges))
                          
                    # get the data, and this will be contingent on the number of dims
                    if vndims_cont == 1:
                        data = contcdf[varname][int(contcorner[0]):int(contcorner[0])+int(contedges[0])]
                    elif vndims_cont == 2:
                        if varname in diagnosticvars:
                            data = contcdf[varname]
                        data = contcdf[varname][int(contcorner[0]):int(contcorner[0])+int(contedges[0]), \
                            int(contcorner[1]):int(contcorner[1])+int(contedges[1])]
                    elif vndims_cont == 3:
                        data = contcdf[varname][int(contcorner[0]):int(contcorner[0])+int(contedges[0]), \
                            int(contcorner[1]):int(contcorner[1])+int(contedges[1]), \
                            int(contcorner[2]):int(contcorner[2])+int(contedges[2])]
                    elif vndims_cont == 4:
                        data = contcdf[varname][int(contcorner[0]):int(contcorner[0])+int(contedges[0]), \
                            int(contcorner[1]):int(contcorner[1])+int(contedges[1]), \
                            int(contcorner[2]):int(contcorner[2])+int(contedges[2]), \
                            int(contcorner[3]):int(contcorner[3])+int(contedges[3])]
                    else:
                        if iburst == 0:
                            print('did not read data')
                    
                    burstcorner = np.zeros(vndims_burst)
                    burstedges = np.ones(vndims_burst) 
                    burstcorner[vdims_burst.index('time')] = iburst
                    burstedges[vdims_burst.index('time')] = burstlength
                                            
                    # since we don't have regular and recurring indeces, we need to handle
                    # situations where the data read is not the maximum number of samples
                    # samples MUST be the second dimension!
                    if ndatasamples < burstlength:
                        issueflags[varname].append(ndatasamples)
                        if len(data.shape) == 1:
                            # start with a filled array
                            burstdata = np.full((1,vshapes_burst[1]),fillval_burst)
                            burstdata[:, 0:ndatasamples] = data[:]
                        elif len(data.shape) == 2:
                            # start with a filled array
                            burstdata = np.full((1,vshapes_burst[1],vshapes_burst[2]),fillval_burst)
                            burstdata[:, 0:ndatasamples] = data[:,:]
                        elif len(data.shape) == 3:
                            # start with a filled array
                            burstdata = np.full((1,vshapes_burst[1],vshapes_burst[2],vshapes_burst[3]),fillval_burst)
                            burstdata[:, 0:ndatasamples,: ] = data[:,:,:]
                        elif len(data.shape) == 4:
                            # start with a filled array
                            burstdata = np.full((1,vshapes_burst[1],vshapes_burst[2],
                                                 vshapes_burst[3], vshapes_burst[4]),fillval_burst)
                            burstdata[:, 0:ndatasamples,:,:] = data[:,:,:,:]
                        elif len(data.shape) == 5:
                            # start with a filled array
                            burstdata = np.full((1,vshapes_burst[1],vshapes_burst[2],
                                                 vshapes_burst[3], vshapes_burst[4],
                                                 vshapes_burst[5]),fillval_burst)
                            burstdata[:, 0:ndatasamples,:,:,:] = data[:,:,:,:,:]
                    else:
                        burstdata = data
                        
                    if ('EPIC' in varname) and (iburst == 0) and verbose :
                        print('\tdata {}'.format(data[1:10]))
                        print('\tburstdata {}'.format(burstdata[1:10]))
                        
                    if (iburst == 0) and verbose:
                        print('\tburstdata.shape = {} burst file dims {}'.format(
                                burstdata.shape, vdims_burst))
                        print('\tvndims_burst = {}'.format(vndims_burst))
                        print('\tdata.shape = {}'.format(data.shape))
                    
                    # TODO -- we can simplify this code my making a function object
                    # some day when I am better at python
                                                                                                               
                    if len(burstdata.shape) == 1:
                        try:
                            burstcdf[varname][iburst] = burstdata[:]
                        except TypeError:
                            # TypeError: int() argument must be a string, a bytes-like object or a number, not 'NoneType'
                            # EPIC_time was given a fill value in the raw file.
                            # this was solved by making sure coordinate variables had fill value set to False
                            if iburst == 0:
                                print('\t{} in Burst file is data type {}, burstdata is type {} and got a TypeError when writing'.format(
                                    varname, bvarobj.dtype, type(burstdata)))
                        except IndexError: # too many indices for array
                            if iburst == 0:
                                print('too many indices for array')
                                print('iburst = {}'.format(iburst))
                                print('burstdata = {}'.format(burstdata))
                        except ValueError: 
                            if iburst == 0:
                                print('ValueError ')
                    elif len(burstdata.shape) == 2:
                        try:
                            burstcdf[varname][iburst,:] = burstdata[:,:]
                        except TypeError:
                            # TypeError: int() argument must be a string, a bytes-like object or a number, not 'NoneType'
                            # EPIC_time was given a fill value in the raw file.
                            # this was solved by making sure coordinate variables had fill value set to False
                            if iburst == 0:
                                print('\t{} in Burst file is data type {}, burstdata is type {} and got a TypeError when writing'.format(
                                    varname, bvarobj.dtype, type(burstdata)))
                        except IndexError: # too many indices for array
                            if iburst == 0:
                                print('too many indices for array')
                                print('iburst = {}'.format(iburst))
                                print('burstdata = {}'.format(burstdata))
                        except ValueError: 
                            if iburst == 0:
                                print('ValueError ')
                    elif len(burstdata.shape) == 3:
                        try:
                            burstcdf[varname][iburst,:,:] = burstdata[:,:,:]
                        except TypeError:
                            if iburst == 0:
                                print('\t{} is data type {} and got a TypeError when writing'.format(
                                    varname, cvarobj.dtype))
                        except IndexError: # too many indices for array
                            if iburst == 0:
                                print('too many indices for array')
                                print('iburst = {}'.format(iburst))
                                print('burstdata = {}'.format(burstdata))
                        except ValueError: 
                            if iburst == 0:
                                print('ValueError cannot reshape array of size 1 into shape (1,150,1,1)')
                                # here we have shapes [time lat lon]
                    elif len(burstdata.shape) == 4:
                        try:
                            burstcdf[varname][iburst,:,:,:] = burstdata[:,:,:,:]
                        except TypeError:
                            if iburst == 0:
                                print('\t{} is data type {} and got a TypeError when writing'.format(
                                    varname, cvarobj.dtype))
                        except IndexError: # too many indices for array
                            if iburst == 0:
                                print('too many indices for array')
                                print('iburst = {}'.format(iburst))
                                print('burstdata = {}'.format(burstdata))
                        except ValueError: 
                            if iburst == 0:
                                print('ValueError cannot reshape array of size 1 into shape (1,150,1,1)')
                                # here we have shapes [time lat lon]
                    elif len(burstdata.shape) == 5:
                        try:
                            burstcdf[varname][iburst,:,:,:] = burstdata[:,:,:,:,:]
                        except TypeError:
                            if iburst == 0:
                                print('\t{} is data type {} and got a TypeError when writing'.format(
                                    varname, cvarobj.dtype))
                        except IndexError: # too many indices for array
                            if iburst == 0:
                                print('too many indices for array')
                                print('iburst = {}'.format(iburst))
                                print('burstdata = {}'.format(burstdata))
                        except ValueError:
                            if iburst == 0:
                                print('got a value error')
                            
                    else:
                        if iburst == 0:
                            print('\tnot set up to write {} dimensions to burst file'.format(
                                    len(burstdata.shape)))
                    
                # end of for iburst in range(nbursts):
            # end of if 'sample' not in vdims_burst:
        # end of if varname not in drop:
    # for cvar in contcdf.variables.items():
    
    burstcdf.start_time = str(num2date(burstcdf['time'][0,0],burstcdf['time'].units))
    burstcdf.stop_time = str(num2date(burstcdf['time'][-1,0],burstcdf['time'].units))
    # TODO compute datetime
    
    burstcdf.close()
    contcdf.close()
    
    print('Finished file conversion at ',dt.datetime.now())

    return issueflags


# utility functions for creating and managing indexes
# make a function to identify the indeces
# and this is where tuples are nice
def find_boundaries(data,edges):
    # using a list of start and end timestamps (edges) that
    # delineate the beginning times and ending times of burts of measurements
    # find the indeces into the data that correspond to these edges
    # the time base may be irregular, it does not matter
    # data is a list object of time stamps from the data
    # edges is a list of tuples of start and end times
    # make the data an numpy array for access to numpy's methods
    nparray = np.array(data)
    
    idx = []
    for edge in edges:
        s = np.where(nparray >= edge[0])
        e = np.where(nparray >= edge[1])
        
        idx.append((int(s[0][0]),int(e[0][0])))
        if not (len(idx) % 100):
            print('.', end='')
            
    print('\n')
    
    return idx


def generate_expected_start_times(cdffile, dim, burst_start_offset, 
                                  burst_interval, burstlength, sample_rate):
    # generate a regular and recurring set of start and end timestamps that
    # delineate the beginning times and ending times of burts of measurements
    # cdfname = name of a continuous time series data file
    # udim = the unlimited or time dimension which we will find the indeces to reshape
    
    # TODO - do this from a first burst first sample time stamp
    print('the file we are looking up times in is {}'.format(cdffile))

    cdf= nc.Dataset(cdffile,format="NETCDF4")
    
    if dim in cdf.dimensions:
        print('the dimension we are operating on is {}'.format(dim))
    else:
        print('{} not found in {} file, aborting'.format(dim, cdffile))
        cdf.close()

    # get the number of bursts based on the elapsed time
    tfirst = num2date(cdf['time'][1],cdf['time'].units)
    tlast = num2date(cdf['time'][-1],cdf['time'].units)
    nbursts = int((tlast-tfirst).total_seconds() / burst_interval)
    
    burst_start_times = []
    for x in range(nbursts): 
        burst_start_times.append(burst_start_offset+x*burst_interval)

    burst_duration = (1 / sample_rate) * burstlength # seconds
    burst_end_times = list(map(lambda x: x+burst_duration,burst_start_times))
    
    print('start times {} such as {}...'.format(len(burst_start_times),burst_start_times[0:5]))
    print('end times {} such as {}...'.format(len(burst_end_times),burst_end_times[0:5]))

    print('the last time is {} seconds from the start of the experiment'.format(cdf['time'][-1]))
    # it turns out later to be convenient to have this as a list of tuples of start and end
    slices = list(map(lambda s,e: (s,e), burst_start_times, burst_end_times))
    print('edge tuples {} such as {}...'.format(len(slices),slices[0:5]))
    
    cdf.close()

    return slices

# TODO -- this is not working. format string is failing
def save_indexes_to_file(cdffile, txtfile, edge_tuples):
    # write indexes to a file with the time stamps for QA/QC
    cdf = nc.Dataset(cdffile,format="NETCDF4")

    tunits = cdf['time'].units
    indexFile = cdffile.split('.')[0]+'indeces.txt'
    with open(indexFile,'w') as outfile:
        outfile.write('Burst Indexes for {}\n'.format(cdffile))
        outfile.write('Burst, start index, end index, number of samples, start time, end time\n')
        for x in enumerate(edge_tuples):
            t0 = num2date(cdf['time'][x[1][0]],tunits)
            t1 = num2date(cdf['time'][x[1][1]],tunits)
            try:
                s = '{}, {},  {}, {}, {}, {}\n'.format(x[0],x[1][0],x[1][1],
                     x[1][1]-x[1][0],t0,t1)
            except:
                s = '{}, {}, {}, , , \n'.format(x[0],x[1][0],x[1][1])
            outfile.write(s)
        
    cdf.close()
    print('Indexes written to {}'.format(indexFile))


if __name__ == "__main__":
    # then we have been run from the command line
    if len(sys.argv) < 3:
        print("%s \nuseage:" % sys.argv[0])
        print("reshapeEPIC(contFile, burstFile, burstlength, [dim2changename], [edges], [vars2omit])" )
        sys.exit(1)
        
    if len(sys.argv) != 3:
        # TODO the keywork pairs do not get into reshape EPIC correctly,
        reshapeEPIC(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:])
    else:
        reshapeEPIC(sys.argv[1], sys.argv[2], sys.argv[3])
        
else:
    # we have been imported
    # the argument passing here works fine
    pass

