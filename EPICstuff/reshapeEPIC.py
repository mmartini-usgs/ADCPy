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

TODO - write a second tool that distributes the samples according to their 
time stamps

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
    
    for s in args:
        print(s)
    for k in kwargs.keys():
        print('{} = {}\n'.format(k,kwargs[k]))
        
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
    
    print(dim)
    print(edges)
    print(burstlength)
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
    
    nbins = len(contcdf['depth'])
    ranges_m = list(map(lambda ibin: 
        contcdf.TRDI_Bin_1_distance_cm/100+ibin*contcdf.TRDI_VBeam_Vertical_Depth_Cell_Size_cm/100,
        range(nbins)))
        
    # ---------------- set up the variables
    
    for cvar in contcdf.variables.items():
        cvarobj = cvar[1]
        print('{} is data type {}'.format(cvarobj.name,cvarobj.dtype))
        try:
            fillValue = cvarobj.getncattr('_FillValue')
            print('\tthe fill value is {}'.format(fillValue))
        except AttributeError:
            print('\tfailed to read the fill value')
            coordset = {} #{'time','EPIC_time','EPIC_time2','depth','lat','lon'}
            if cvarobj.name in coordset:
                # this will avoid the typerror when EPIC_time is written and is
                # not a good solution because it causes zeros in a time variable
                fillValue = 0 # TODO this is not a great solution
            else:
                fillValue = None
        
        if fillValue == None:
            fillValue = False
                
        print('\tfillValue in burst file will be set to {} (if None, then False)'.format(fillValue))
        
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
        pass # likely depth exists if this happens               
        
    # TODO - if edges is length 1, then we need to create the burst edges here
    
    # --------- populate the file
    # note that we don't have to change from a define to a read mode here
    
    # coordinate variables are small(er) and can be done at once, be sure to use
    # generative methods
    print('\nNow populating data for {} bursts'.format(nbursts))
    burstcdf['burst'][:] = list(range(nbursts))
    burstcdf['sample'][:] = list(range(burstlength))
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
        
            print('{} in Continuous file is data type {}'.format(varname,cvarobj.dtype))
            print('\tin Burst file it is data type {}'.format(bvarobj.dtype))
            vdims_cont = cvarobj.dimensions
            vshapes_cont = cvarobj.shape
            vndims_cont = len(cvarobj.dimensions)
            vdims_burst = burstcdf[varname].dimensions
            vshapes_burst = burstcdf[varname].shape
            vndims_burst = len(burstcdf[varname].dimensions)
            try:
                fillval_burst = burstcdf[varname].getncattr('_FillValue')
            except:
                if 'EPIC' in varname:
                    # EPIC was ending up with odd fill values in the raw file
                    # TODO this will cause pain later but we just want to get
                    # data into the file
                    # this will avoid the typerror when EPIC_time is written and is
                    # not a good solution
                    fillval_burst = 0 
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
                    
                    #print('\tburst {} is {} samples from time index {} to {}'.format(
                    #    iburst,ndatasamples,contcorner,contedges))
                    if 'depth' in vdims_cont:
                        contedges[vdims_cont.index('depth')] = vshapes_cont[vdims_cont.index('depth')]

                    # this was necessary in the MATLAB version, not sure I need it here yet
                    '''
                    # are we going beyond the end of the data?  If so - adjust the number of samples
                    if contcorner[vdims_cont.index('time')]+contedges[vdims_cont.index('time')] > \
                            len(contcdf['time']):
                        contedges[vdims_cont.index('time')] = \
                            len(contcdf['time'])-contcorner[vdims_cont.index('time')]
                        print('not enough samples for burst {} - truncating'.format(iburst))
                    # are we past the end of the file?
                    if contcorner[vdims_cont.index('time')] > contedges[vdims_cont.index('time')]:
                        print('end of data reached at burst {}'.format(iburst))
                        break
                    '''
                    
                    # get the data, and this will be contingent on the number of dims
                    if vndims_cont == 1:
                        data = contcdf[varname][int(contcorner[0]):int(contcorner[0])+int(contedges[0])]
                    elif vndims_cont == 2:
                        if varname in diagnosticvars:
                            data = contcdf[varname]
                        data = contcdf[varname][int(contcorner[0]):int(contcorner[0])+int(contedges[0]), \
                            int(contcorner[1]):int(contcorner[1])+int(contedges[1])]
                    elif vndims_cont == 3:
                        data = contcdf[varname][int(contcorner[0]):int(contcorner[0])+int(contedges[0])] \
                            [int(contcorner[1]):int(contcorner[1])+int(contedges[1])] \
                            [int(contcorner[2]):int(contcorner[2])+int(contedges[2])]
                    elif vndims_cont == 4:
                        data = contcdf[varname][int(contcorner[0]):int(contcorner[0])+int(contedges[0])] \
                            [int(contcorner[1]):int(contcorner[1])+int(contedges[1])] \
                            [int(contcorner[2]):int(contcorner[2])+int(contedges[2])] \
                            [int(contcorner[3]):int(contcorner[3])+int(contedges[3])]
                    else:
                        pass
                    
                    burstcorner = np.zeros(vndims_burst)
                    burstedges = np.ones(vndims_burst) 
                    burstcorner[vdims_burst.index('time')] = iburst
                    burstedges[vdims_burst.index('time')] = burstlength
                                            
                    # since we don't have regular and recurring indeces, we need to handle
                    # situations where the data read is not the maximum number of samples   
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
                    else:
                        burstdata = data
                        
                    if 'EPIC' in varname and iburst==0 :
                        print('\tdata {}'.format(data[1:10]))
                        print('\tburstdata {}'.format(burstdata[1:5,1:10]))
                        print('\tvndims_cont {}'.format(vndims_cont))
                        print('\tvndims_burst {}'.format(vndims_burst))
                                                                                       
                    if len(burstcdf[varname].shape) == 1:
                        pass # no known scalars
                    elif vndims_burst == 2:
                        try:
                            burstcdf[varname][iburst,:] = burstdata[:,:]
                        except TypeError:
                            # TypeError: int() argument must be a string, a bytes-like object or a number, not 'NoneType'
                            # burstdata is object numpy ndarray
                            # varname is str
                            # iburst is int
                            # there were Nones in burstdata, why this is a problem only for EPIC_time I don't know
                            # EPIC_time was given a fill vale in the raw file.
                            if iburst == 0:
                                print('\t{} in Burst file is data type {}, burstdata is type {} and got a TypeError when writing'.format(
                                    varname, bvarobj.dtype, type(burstdata)))
                    elif vndims_burst == 3:
                        try:
                            burstcdf[varname][iburst,:,:] = burstdata[:,:,:]
                        except TypeError:
                            if iburst == 0:
                                print('\t{} is data type {} and got a TypeError when writing'.format(
                                    varname, cvarobj.dtype))
                    else:
                        pass
                    
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



if __name__ == "__main__":
    # then we have been run from the command line
    if len(sys.argv) < 3:
        print("%s \nuseage:" % sys.argv[0])
        print("reshapeEPIC(contFile, burstFile, burstlength, [dim2changename], [edges], [vars2omit])" )
        sys.exit(1)
        
    if len(sys.argv) != 3:
        # TODO the keywork pairs do not get into reshape EPIC correctly,
        # a command line call of 
        # python reshapeEPIC.py junk1 junk2 3 dim='time' edges=[1,2,3,4,5]
        # parses to
        # reshapeEPIC.py running on python 3.6.6 | packaged by conda-forge | (default, Jul 26 2018, 11:48:23) [MSC v.1900 64 bit (AMD64)]
        #junk1
        #junk2
        #3
        #["dim='time'", 'edges=[1,2,3,4,5]']
        reshapeEPIC(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:])
    else:
        reshapeEPIC(sys.argv[1], sys.argv[2], sys.argv[3])
        
else:
    # we have been imported
    # the argument passing here works fine
    pass

