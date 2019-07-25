"""
This code takes the raw currents portion of the adcp data from the splitter
program pd0.py and outputs raw current profile data to a netCDF4 file.
If you have a file with wave packets data, the splitter must be run first.

As a script:

python TRDIpd0tonetcdf.py [path] pd0File cdfFile

or dopd0file(pd0File, cdfFile, goodens, serialnum, timetype, delta_t_to_use)

where:
    path            is a path to prepend to the following
    pd0File         is path of raw PD0 format input file with current ensembles
    cdfFile         is path of a netcdf4 EPIC compliant output file
    goodens         [start, end]  ensembles to export.  end = -1 for all ensembles in file
    serialnum       serial number of the instrument, a string
    timetype        "CF" for CF conventions, "EPIC" for EPIC conventions
    delta_t_to_use  time between ensembles, in seconds, a string.  15 min profiles would be 900

note that file names and paths may not include spaces
    
As a module:
import TRDIpd0tonetcdf as pd0

Notes:
    time and time2, the EPIC convention for netCDF, is not used here so that
    the resulting very large files generated can be reduced using existing 
    python too ls such as xarrays

Programmed according to the TRDI Workhorse Commands and Output Data Format document, March 2005
"""

# 10/4/2018 remove valid_range as it causes too many downstream problems
# 1/25/2017 MM got this running on old Workhorse ADCP data

import sys
import struct
import math
import numpy as np 
from netCDF4 import Dataset
import datetime as dt
from EPICstuff.EPICmisc import cftime2EPICtime
from EPICstuff.EPICmisc import ajd


def dopd0file(pd0File, cdfFile, goodens, serialnum, timetype, delta_t_to_use):

    # TODO figure out a better way to handle this situation
    # need this check in case this function is used as a stand alone function

    # this is necessary so that this function does not change the value
    # in the calling function

    ens2process = goodens[:]
    verbose = True  # diagnostic, True = turn on output, False = silent
    
    maxens, ensLen, ensData, dataStartPosn = analyzepd0file(pd0File, verbose)
    
    infile = open(pd0File, 'rb')
    
    infile.seek(dataStartPosn)
    
    if (ens2process[1] < 0) or ens2process[1] == np.inf:
        ens2process[1] = maxens
           
    # we are good to go, get the output file ready
    print('Setting up netCDF file %s' % cdfFile)
    cdf, cf_units = setupCdf(cdfFile, ensData, ens2process, serialnum, timetype, delta_t_to_use)
    # we want to save the time stamp from this ensemble since it is the
    # time from which all other times in the file will be relative to
    t0 = ensData['VLeader']['dtobj']

    cdfIdx = 0
    ensCount = 0
    verbose = False  # diagnostic, True = turn on output, False = silent
    nslantbeams = 4
        
    # priming read - for the while loop
    # note that ensemble lengths can change in the middle of the file!
    # horribly inefficient, but here we go, one step backward, two forward...
    bookmark = infile.tell() # save beginning of next ensemble
    # need to read the header from the file to know the ensemble size
    Header = readTRDIHeader(infile)
    if Header['sourceID'] != b'\x7f':
        print('non-currents ensemble found at %d' % bookmark)
    
    if ensLen != Header['nbytesperens']+2:
        ensLen = Header['nbytesperens']+2 # update to what we have
    
    # go back to where this ensemble started before we checked the header
    infile.seek(bookmark)
    ens = infile.read(ensLen)   
        
    while len(ens) > 0:
        # print('-- ensemble %d length %g, file position %g' % (ensCount, len(ens), infile.tell()))
        # print(ensData['Header'])
        ensData, ensError = parseTRDIensemble(ens, verbose)
        
        if (ensError == 'None') and (ensCount >= ens2process[0]):
            # write to netCDF
            if cdfIdx == 0:
                print('--- first ensembles read at %s and TRDI #%d' % (              
                    ensData['VLeader']['timestr'], ensData['VLeader']['Ensemble_Number']))
                
            varobj = cdf.variables['Rec']
            try:
                varobj[cdfIdx] = ensData['VLeader']['Ensemble_Number']
            except:
                # here we have reached the end of the netCDF file
                cdf.close()
                infile.close()
                return

            # time calculations done when vleader is read
            if timetype == 'CF':
                varobj = cdf.variables['time']
                elapsed = ensData['VLeader']['dtobj']-t0 # timedelta
                elapsed_sec = elapsed.total_seconds()
                # TODO - suspect my EPIC_time woes may be caused here
                # is elapsed_sec rolling over?
                if elapsed_sec == 0:
                    print('elapsed seconds from ensemble {} is {}'.format(ensCount,elapsed_sec))
                    
                varobj[cdfIdx] = elapsed_sec   
                t1, t2 = cftime2EPICtime(elapsed_sec,cf_units)
                varobj = cdf.variables['EPIC_time']
                #varobj[cdfIdx] = ensData['VLeader']['EPIC_time']
                varobj[cdfIdx] = t1
                varobj = cdf.variables['EPIC_time2']
                #varobj[cdfIdx] = ensData['VLeader']['EPIC_time2']
                varobj[cdfIdx] = t2
                
            else:
                varobj = cdf.variables['time']
                varobj[cdfIdx] = ensData['VLeader']['EPIC_time']
                varobj = cdf.variables['time2']
                varobj[cdfIdx] = ensData['VLeader']['EPIC_time2']
                varobj = cdf.variables['cf_time']
                elapsed = ensData['VLeader']['dtobj']-t0 # timedelta
                elapsed_sec = elapsed.total_seconds()
                varobj[cdfIdx] = elapsed_sec              
            
            # diagnostic
            if (ens2process[1]-ens2process[0]-1)<100:
                print('%d %15.8f %s' % (ensData['VLeader']['Ensemble_Number'], 
                                        ensData['VLeader']['julian_day_from_julian'],
                                        ensData['VLeader']['timestr']))
    
            varobj = cdf.variables['sv']
            varobj[cdfIdx] = ensData['VLeader']['Speed_of_Sound']

            for i in range(nslantbeams):
                varname = "vel%d" % (i+1)
                varobj = cdf.variables[varname]
                varobj[cdfIdx,:] = ensData['VData'][i,:]

            for i in range(nslantbeams):
                varname = "cor%d" % (i+1)
                varobj = cdf.variables[varname]
                varobj[cdfIdx,:] = ensData['CData'][i,:]

            for i in range(nslantbeams):
                varname = "att%d" % (i+1)
                varobj = cdf.variables[varname]
                varobj[cdfIdx,:] = ensData['IData'][i,:]

            if ('GData' in ensData):
                for i in range(nslantbeams):
                    varname = "PGd%d" % (i+1)
                    varobj = cdf.variables[varname]
                    varobj[cdfIdx,:] = ensData['GData'][i,:]

            varobj = cdf.variables['Rec']
            varobj[cdfIdx] = ensData['VLeader']['Ensemble_Number']
            varobj = cdf.variables['Hdg']
            varobj[cdfIdx] = ensData['VLeader']['Heading']
            varobj = cdf.variables['Ptch']
            varobj[cdfIdx] = ensData['VLeader']['Pitch']
            varobj = cdf.variables['Roll']
            varobj[cdfIdx] = ensData['VLeader']['Roll']
            varobj = cdf.variables['HdgSTD']
            varobj[cdfIdx] = ensData['VLeader']['H/Hdg_Std_Dev']
            varobj = cdf.variables['PtchSTD']
            varobj[cdfIdx] = ensData['VLeader']['P/Pitch_Std_Dev']
            varobj = cdf.variables['RollSTD']
            varobj[cdfIdx] = ensData['VLeader']['R/Roll_Std_Dev']
            varobj = cdf.variables['Tx']
            varobj[cdfIdx] = ensData['VLeader']['Temperature']
            varobj = cdf.variables['S']
            varobj[cdfIdx] = ensData['VLeader']['Salinity']
            varobj = cdf.variables['xmitc']
            varobj[cdfIdx] = ensData['VLeader']['Xmit_Current']
            varobj = cdf.variables['xmitv']
            varobj[cdfIdx] = ensData['VLeader']['Xmit_Voltage']
            varobj = cdf.variables['Ambient_Temp']
            varobj[cdfIdx] = ensData['VLeader']['Ambient_Temp']
            varobj = cdf.variables['Pressure+']
            varobj[cdfIdx] = ensData['VLeader']['Pressure_(+)']
            varobj = cdf.variables['Pressure-']
            varobj[cdfIdx] = ensData['VLeader']['Pressure_(-)']
            varobj = cdf.variables['Attitude_Temp']
            varobj[cdfIdx] = ensData['VLeader']['Attitude_Temp']
            varobj = cdf.variables['EWD1']
            varobj[cdfIdx] = int(ensData['VLeader']['Error_Status_Word_Low_16_bits_LSB'])
            varobj = cdf.variables['EWD2']
            varobj[cdfIdx] = int(ensData['VLeader']['Error_Status_Word_Low_16_bits_MSB'])
            varobj = cdf.variables['EWD3']
            varobj[cdfIdx] = int(ensData['VLeader']['Error_Status_Word_High_16_bits_LSB'])
            varobj = cdf.variables['EWD4']
            varobj[cdfIdx] = int(ensData['VLeader']['Error_Status_Word_High_16_bits_MSB'])
            
            if ensData['FLeader']['Depth_sensor_available'] == 'Yes':
                varobj = cdf.variables['Pressure']
                varobj[cdfIdx] = ensData['VLeader']['Pressure_deca-pascals']
                varobj = cdf.variables['PressVar']
                varobj[cdfIdx] = ensData['VLeader']['Pressure_variance_deca-pascals']

            # add bottom track data write to cdf here
            if ('BTData' in ensData):
                if ensData['BTData']['Mode'] == 0:
                    varobj = cdf.variables['BTRmin']
                    varobj[cdfIdx] = ensData['BTData']['Ref_Layer_Min']
                    varobj = cdf.variables['BTRnear']
                    varobj[cdfIdx] = ensData['BTData']['Ref_Layer_Near']
                    varobj = cdf.variables['BTRfar']
                    varobj[cdfIdx] = ensData['BTData']['Ref_Layer_Far']

                varnames = ('BTWe','BTWu','BTWv','BTWd')
                for i in range(nslantbeams):
                    varname = "BTR%d" % (i+1)
                    varobj = cdf.variables[varname]
                    varobj[cdfIdx] = ensData['BTData']['BT_Range'][i]
                    if ensData['FLeader']['Coord_Transform'] == 'EARTH':
                        varobj = cdf.variables[varnames[i]]
                    else:
                        varname = "BTV%d" % (i+1)
                        varobj = cdf.variables[varname]
                    
                    varobj[cdfIdx] = ensData['BTData']['BT_Vel'][i]
                    varname = "BTc%d" % (i+1)
                    varobj = cdf.variables[varname]
                    varobj[cdfIdx] = ensData['BTData']['BT_Corr'][i]
                    varname = "BTe%d" % (i+1)
                    varobj = cdf.variables[varname]
                    varobj[cdfIdx] = ensData['BTData']['BT_Amp'][i]
                    varname = "BTp%d" % (i+1)
                    varobj = cdf.variables[varname]
                    varobj[cdfIdx] = ensData['BTData']['BT_PGd'][i]
                    varname = "BTRSSI%d" % (i+1)
                    varobj = cdf.variables[varname]
                    varobj[cdfIdx] = ensData['BTData']['RSSI_Amp'][i]
                    
                    if ensData['BTData']['Mode'] == 0:
                        varobj[cdfIdx] = ensData['BTData']['Ref_Layer_Vel'][i]
                        varname = "BTRc%d" % (i+1)
                        varobj = cdf.variables[varname]
                        varobj[cdfIdx] = ensData['BTData']['Ref_Layer_Corr'][i]
                        varname = "BTRi%d" % (i+1)
                        varobj = cdf.variables[varname]
                        varobj[cdfIdx] = ensData['BTData']['Ref_Layer_Amp'][i]
                        varname = "BTRp%d" % (i+1)
                        varobj = cdf.variables[varname]
                        varobj[cdfIdx] = ensData['BTData']['Ref_Layer_PGd'][i]
                        


            if ('VBeamVData' in ensData):
                if ensData['VBeamLeader']['Vertical_Depth_Cells'] == ensData['FLeader']['Number_of_Cells']:
                    varobj = cdf.variables['vel5']
                    varobj[cdfIdx,:] = ensData['VBeamVData']
                    varobj = cdf.variables['cor5']
                    varobj[cdfIdx,:] = ensData['VBeamCData']
                    varobj = cdf.variables['att5']
                    varobj[cdfIdx,:] = ensData['VBeamIData']
                    if ('VBeamGData' in ensData):
                        varobj = cdf.variables['PGd5']
                        varobj[cdfIdx,:] = ensData['VBeamGData']
                    
            if ('WaveParams' in ensData):
                # we can get away with this because the key names and var names are the same
                for key in ensData['WaveParams']:
                    varobj = cdf.variables[key]
                    varobj[cdfIdx] = ensData['WaveParams'][key]
    
            if ('WaveSeaSwell' in ensData):
                # we can get away with this because the key names and var names are the same
                for key in ensData['WaveSeaSwell']:
                    varobj = cdf.variables[key]
                    varobj[cdfIdx] = ensData['WaveSeaSwell'][key]

            cdfIdx += 1
            
        elif ensError == 'no ID':
            print('Stopping because ID tracking lost')
            infile.close()
            cdf.close()
            sys.exit(1)
            
        ensCount += 1
        
        if ensCount > maxens:
            print('stopping at estimated end of file ensemble %d' % ens2process[1])
            break        
        
        # if maxens < 100:  n=10
        # elif (maxens > 100) and (maxens < 1000): n=100
        # elif (maxens > 1000) and (maxens < 10000): n=1000
        # elif (maxens > 10000) and (maxens < 100000): n=10000
        # elif (maxens > 100000) and (maxens < 1000000): n=100000
        # else: n = 1000000
        n = 10000
        
        ensf, ensi = math.modf(ensCount/n)
        if ensf == 0:
            print('%d ensembles read at %s and TRDI #%d' % (ensCount,              
                ensData['VLeader']['dtobj'], ensData['VLeader']['Ensemble_Number']))
        
        if ensCount >= ens2process[1]-1:
            print('stopping at requested ensemble %d' % ens2process[1])
            break
        
        # note that ensemble lengths can change in the middle of the file!
        # TODO - is there a fastwe way to do this??
        bookmark = infile.tell() # save beginning of next ensemble
        # TODO - since we are jumping around, we should check here to see
        #   how close to the end of the file we are - if it is within one
        #   header length - we are done
        #   need to read the header from the file to know the ensemble size
        Header = readTRDIHeader(infile)

        if Header == None:
            # we presume this is the end of the file, since we don't have header info
            print('end of file reached with incomplete header')
            break
            
        if Header['sourceID'] != b'\x7f':
            print('non-currents ensemble found at %d' % bookmark)
        
        if ensLen != Header['nbytesperens']+2:
            ensLen = Header['nbytesperens']+2 # update to what we have
        
        # TODO - fix this so that we aren't going back and forth, it is probably really slow
        # go back to where this ensemble started before we checked the header
        infile.seek(bookmark)
        ens = infile.read(ensLen)
        
    else:  # while len(ens) > 0:
        print('end of file reached')
        
    if ensCount < maxens:
        print('end of file reached after %d ensembles, less than estimated in the file' % ensCount)
    elif ensCount > maxens:
        print('end of file reached after %d ensembles, more than estimated in the file' % ensCount)
    
    infile.close()
    cdf.close()
    
    print('%d ensembles read, %d records written' % (ensCount, cdfIdx))

    return ensCount, cdfIdx, ensError
    
    
def matrixTranspose( matrix ):
    if not matrix: return []
    return [ [ row[ i ] for row in matrix ] for i in range( len( matrix[ 0 ] ) ) ]

# write a dictionary to netCDF attributes
def writeDict2atts(cdfobj, d, tag):
    
    i = 0;
    # first, convert as many of the values in d to numbers as we can
    for key in iter(d):
        if type(d[key]) == str:
            try:
                d[key] = float(d[key])
            except ValueError:
                # we really don't need to print here, 
                # but python insists we do something
                #print('   can\'t convert %s to float' % key)
                i += 1

    for key in iter(d):
        newkey = tag + key
        try:
            cdfobj.setncattr(newkey,d[key])
        except:
            print('can\'t set %s attribute' % key)
    
    # return the dictionary it's numerized values 
    return d
    

def parseTRDIensemble(ensbytes, verbose):
    ensData = {}
    ensError = 'None'
    ensData['Header'] = parseTRDIHeader(ensbytes)
    
    for i in range(ensData['Header']['ndatatypes']):
        # go to each offset and parse depending on what we find
        offset = ensData['Header']['offsets'][i]
        #raw, val = __parseTRDIushort(ensbytes, offset)
        val = struct.unpack('<H',ensbytes[offset:offset+2])[0]
        if val == 0: # \x00\x00
            if verbose: print('Fixed Leader found at %g' % offset)
            ensData['FLeader'] = parseTRDIFixedLeader(ensbytes, offset)
            # we need this to decode the other data records
            ncells = int(ensData['FLeader']['Number_of_Cells'])
            nbeams = 4 # the 5th beam has it's own record
        elif val == 128: # \x80\x00
            if verbose: print('Variable Leader found at %g' % offset)
            ensData['VLeader'] = parseTRDIVariableLeader(ensbytes, offset)
            #print(VLeader)
        elif val == 256:  # raw == b'\x00\x01': 256
            if verbose: print('Velocity found at %g' % offset)
            ensData['VData'] = parseTRDIVelocity(ensbytes, offset, ncells, nbeams)
        elif val == 512: #raw == b'\x00\x02': 
            if verbose: print('Correlation found at %g' % offset)
            ensData['CData'] = parseTRDICorrelation(ensbytes, offset, ncells, nbeams)
        elif val == 768: #raw == b'\x00\x03':
            if verbose: print('Intensity found at %g' % offset)
            ensData['IData'] = parseTRDIIntensity(ensbytes, offset, ncells, nbeams)
        elif val == 1024: #raw == b'\x00\x04':
            if verbose: print('PGood found at %g' % offset)
            ensData['GData'] = parseTRDIPercentGood(ensbytes, offset, ncells, nbeams)
        elif val == 1280: #raw == b'\x00\x05':
            if verbose: print('Status profile found at %g' % offset)
        elif val == 1536: #raw == b'\x00\x06':
            if verbose: print('BT found at %g' % offset)
            ensData['BTData'] = parseTRDIBottomTrack(ensbytes, offset, nbeams)
        elif val == 1792: #raw == b'\x00\x07':
            # this not defined in TRDI docs
            pass
        elif val == 2048: #raw == b'\x00\x08':
            if verbose: print('MicroCAT data found at %g' % offset)
        elif val == 12800: #raw == b'\x00\x32': #12800
            if verbose: print('Instrument transformation found at %g' % offset)
            ensData['XformMatrix'] = parseTRDIxformMatrix(ensbytes, offset, nbeams)
        elif val == 28672: #raw == b'\x00\x70':
            if verbose: print('V Series sytem config found at %g' % offset)
            ensData['VSysConfig'] = parseTRDIVSysConfig(ensbytes, offset)
        elif val == 28673: #raw == b'\x01\x70':
            if verbose: print('V Series ping setup found at %g' % offset)
            ensData['VPingSetup'] = parseTRDIVPingSetup(ensbytes,offset)
        elif val == 28674: #raw == b'\x02\x70':
            if verbose: print('V Series ADC Data found at %g' % offset)
            # currently not defined well in TRDI docs
        elif val == 28675: #raw == b'\x03\x70':
            if verbose: print('V Series System Configuration Data found at %g' % offset)
            # currently not defined well in TRDI docs
        elif val == 3841: #raw == b'\x01\x0f':
            if verbose: print('Vertical Beam Leader Data found at %g' % offset)
            ensData['VBeamLeader'] = parseTRDIVBeamLeader(ensbytes, offset)
        elif val == 2560: #raw == b'\x00\x0a':
            if verbose: print('Vertical Beam Velocity Data found at %g' % offset)
            ensData['VBeamVData'] = parseTRDIVertVelocity(ensbytes, offset, 
                ensData['VBeamLeader']['Vertical_Depth_Cells'])
        elif val == 2816: #raw == b'\x00\x0b':
            if verbose: print('Vertical Beam Correlation Data found at %g' % offset)
            ensData['VBeamCData'] = parseTRDIVertCorrelation(ensbytes, offset, 
                ensData['VBeamLeader']['Vertical_Depth_Cells'])
        elif val == 3072: #raw == b'\x00\x0c':
            if verbose: print('Vertical Beam Amplitude Data found at %g' % offset)
            ensData['VBeamIData'] = parseTRDIVertIntensity(ensbytes, offset, 
                ensData['VBeamLeader']['Vertical_Depth_Cells'])
        elif val == 3328: #raw == b'\x00\x0d':
            if verbose: print('Vertical Beam Percent Good Data found at %g' % offset)
            ensData['VBeamGData'] = parseTRDIVertPercentGood(ensbytes, offset, 
                ensData['VBeamLeader']['Vertical_Depth_Cells'])
        elif val == 28676: #raw == b'\x40\x70':
            if verbose: print('V Series Event Log Data found at %g' % offset)
        elif val == 11: # raw == b'\x0b\x00':
            if verbose: print('Wavesmon 4 Wave Parameters found at %g' % offset)
            ensData['WaveParams'] = parseTRDIWaveParameters(ensbytes, offset)
        elif val == 12: #raw == b'\x0c\x00':
            if verbose: print('Wavesmon 4 Sea and Swell found at %g' % offset)
            ensData['WaveSeaSwell'] = parseTRDIWaveSeaSwell(ensbytes, offset)
        else:
            print('ID %d unrecognized at %g' % (val, offset))
            ensError = 'no ID'
        
    csum = __computeChecksum(ensbytes)
    if csum != (ensbytes[-2]+(ensbytes[-1]<<8)):
        ensError = 'checksum failure'
        
    return ensData, ensError


def setupCdf(fname, ensData, gens, serialnum, timetype, delta_t_to_use):
     
    # note that 
    # f4 = 4 byte, 32 bit float
    # maxfloat = 3.402823*10**38;
    intfill = -32768
    floatfill = 1E35

    nens = gens[1]-gens[0]-1
    print('creating netCDF file %s with %d records' % (fname, nens))
    
    cdf = Dataset(fname, "w", clobber=True, format="NETCDF4")
    
    # dimensions, in EPIC order
    cdf.createDimension('time', nens)
    cdf.createDimension('depth', ensData['FLeader']['Number_of_Cells'])
    cdf.createDimension('lat', 1)
    cdf.createDimension('lon', 1)
    
    # write global attributes
    cdf.history = "translated to netCDF by TRDIpd0tonetcdf.py"
    cdf.sensor_type = "TRDI"
    cdf.serial_number = serialnum
    cdf.DELTA_T = delta_t_to_use
    cdf.sample_rate = ensData['FLeader']['Time_Between_Ping Groups']
    
    writeDict2atts(cdf, ensData['FLeader'], "TRDI_")
    
    varobj = cdf.createVariable('Rec', 'u4', 'time', fill_value=intfill)
    varobj.units = "count"
    varobj.long_name = "Ensemble Number"
    # the ensemble number is a two byte LSB and a one byte MSB (for the rollover)
    # varobj.valid_range = [0, 2**23]

    # it's not yet clear which way to go with this.  python tools like xarray 
    # and panoply demand that time be a CF defined time.
    # USGS CMG MATLAB tools need time and time2
    if timetype == 'EPIC':
        # we include cf_time for cf compliance and use by python packages like xarray
        # if f8, 64 bit is not used, time is clipped
        # for ADCP fast sampled, single ping data, need millisecond resolution
        varobj = cdf.createVariable('cf_time', 'f8', 'time')
        # for cf convention, always assume UTC for now, and use the UNIX Epoch as the reference
        varobj.units = "seconds since %d-%d-%d %d:%d:%f 0:00" % (ensData['VLeader']['Year'],
            ensData['VLeader']['Month'], ensData['VLeader']['Day'],ensData['VLeader']['Hour'],
            ensData['VLeader']['Minute'], ensData['VLeader']['Second']+
            ensData['VLeader']['Hundredths']/100)
        varobj.standard_name = "time"
        varobj.axis = "T"    
        # we include time and time2 for EPIC compliance
        varobj = cdf.createVariable('time','u4',('time'))
        varobj.units = "True Julian Day"
        varobj.epic_code = 624
        varobj.datum = "Time (UTC) in True Julian Days: 2440000 = 0000 h on May 23, 1968"
        varobj.NOTE = "Decimal Julian day [days] = time [days] + ( time2 [msec] / 86400000 [msec/day] )"    
        varobj = cdf.createVariable('time2','u4',('time'))
        varobj.units = "msec since 0:00 GMT"
        varobj.epic_code = 624
        varobj.datum = "Time (UTC) in True Julian Days: 2440000 = 0000 h on May 23, 1968"
        varobj.NOTE = "Decimal Julian day [days] = time [days] + ( time2 [msec] / 86400000 [msec/day] )"
        cf_units = ""
    else:
        # cf_time for cf compliance and use by python packages like xarray
        # if f8, 64 bit is not used, time is clipped
        # for ADCP fast sampled, single ping data, need millisecond resolution
        varobj = cdf.createVariable('time','f8',('time'))
        # for cf convention, always assume UTC for now, and use the UNIX Epoch as the reference
        varobj.units = "seconds since %d-%d-%d %d:%d:%f 0:00" % (ensData['VLeader']['Year'],
            ensData['VLeader']['Month'],ensData['VLeader']['Day'],ensData['VLeader']['Hour'],
            ensData['VLeader']['Minute'],ensData['VLeader']['Second']+
            ensData['VLeader']['Hundredths']/100)
        cf_units = "seconds since %d-%d-%d %d:%d:%f 0:00" % (ensData['VLeader']['Year'],
            ensData['VLeader']['Month'],ensData['VLeader']['Day'],ensData['VLeader']['Hour'],
            ensData['VLeader']['Minute'],ensData['VLeader']['Second']+
            ensData['VLeader']['Hundredths']/100)
        varobj.standard_name = "time"
        varobj.axis = "T"
        varobj.type = "UNEVEN"
        # we include time and time2 for EPIC compliance
        # this statement resulted in a fill value of -1??
        #varobj = cdf.createVariable('EPIC_time','u4',('time'))
        varobj = cdf.createVariable('EPIC_time','u4',('time'),fill_value=False)
        varobj.units = "True Julian Day"
        varobj.epic_code = 624
        varobj.datum = "Time (UTC) in True Julian Days: 2440000 = 0000 h on May 23, 1968"
        varobj.NOTE = "Decimal Julian day [days] = time [days] + ( time2 [msec] / 86400000 [msec/day] )"    
        # this statement resulted in a fill value of -1??
        #varobj = cdf.createVariable('EPIC_time2','u4',('time'))
        varobj = cdf.createVariable('EPIC_time2','u4',('time'),fill_value=False)
        varobj.units = "msec since 0:00 GMT"
        varobj.epic_code = 624
        varobj.datum = "Time (UTC) in True Julian Days: 2440000 = 0000 h on May 23, 1968"
        varobj.NOTE = "Decimal Julian day [days] = time [days] + ( time2 [msec] / 86400000 [msec/day] )"  
        
    varobj = cdf.createVariable('bindist','f4',('depth'),fill_value=floatfill)
    # note name is one of the netcdf4 reserved attributes, use setncattr
    varobj.setncattr('name', "bindist")
    varobj.units = "m"
    varobj.long_name = "bin distance from instrument for slant beams"
    varobj.epic_code = 0
    #varobj.valid_range = [0 0]
    varobj.NOTE = "distance is calculated from center of bin 1 and bin size"
    bindist = []
    for idx in range(ensData['FLeader']['Number_of_Cells']):
        bindist.append(idx*(ensData['FLeader']['Depth_Cell_Length_cm']/100)+ensData['FLeader']['Bin_1_distance_cm']/100)
    varobj[:] = bindist[:]
        
    varobj = cdf.createVariable('depth','f4',('depth')) # no fill for ordinates
    varobj.units = "m"
    varobj.long_name = "distance from transducer, depth placeholder"
    varobj.center_first_bin_m = ensData['FLeader']['Bin_1_distance_cm']/100
    varobj.blanking_distance_m = ensData['FLeader']['Blank_after_Transmit_cm']/100
    varobj.bin_size_m = ensData['FLeader']['Depth_Cell_Length_cm']/100
    varobj.bin_count = ensData['FLeader']['Number_of_Cells']
    varobj[:] = bindist[:]

    varobj = cdf.createVariable('sv','f4',('time'),fill_value=floatfill)
    varobj.units = "m s-1"
    varobj.long_name = "sound velocity (m s-1)"
    #varobj.valid_range = [1400, 1600]
    
    for i in range(4):
        varname = "vel%d" % (i+1)
        varobj = cdf.createVariable(varname,'f4',('time','depth'),fill_value=floatfill)
        varobj.units = "mm s-1"
        varobj.long_name = "Beam %d velocity (mm s-1)" % (i+1)
        varobj.epic_code = 1277+i
        #varobj.valid_range = [-32767, 32767]
    
    for i in range(4):
        varname = "cor%d" % (i+1)
        varobj = cdf.createVariable(varname,'u2',('time','depth'),fill_value=intfill)
        varobj.units = "counts"
        varobj.long_name = "Beam %d correlation" % (i+1)
        varobj.epic_code = 1285+i
        #varobj.valid_range = [0, 255]

    for i in range(4):
        varname = "att%d" % (i+1)
        varobj = cdf.createVariable(varname,'u2',('time','depth'),fill_value=intfill)
        varobj.units = "counts"
        varobj.epic_code = 1281+i
        varobj.long_name = "ADCP attenuation of beam %d" % (i+1)
        #varobj.valid_range = [0, 255]

    if ('GData' in ensData):
        for i in range(4):
            varname = "PGd%d" % (i+1)
            varobj = cdf.createVariable(varname,'u2',('time','depth'),fill_value=intfill)
            varobj.units = "counts"
            varobj.long_name = "Percent Good Beam %d" % (i+1)
            varobj.epic_code = 1241+i
            #varobj.valid_range = [0, 100]

    varobj = cdf.createVariable('Hdg','f4',('time'),fill_value=floatfill)
    varobj.units = "hundredths of degrees"
    varobj.long_name = "INST Heading"
    varobj.epic_code = 1215
    varobj.heading_alignment = ensData['FLeader']['Heading_Alignment_Hundredths_of_Deg']
    varobj.heading_bias = ensData['FLeader']['Heading_Bias_Hundredths_of_Deg']
    #varobj.valid_range = [0, 36000]
    if ensData['FLeader']['Heading_Bias_Hundredths_of_Deg'] == 0:
        varobj.NOTE_9 = "no heading bias was applied by EB during deployment or by wavesmon"
    else:
        varobj.NOTE_9 = "a heading bias was applied by EB during deployment or by wavesmon"
    
    varobj = cdf.createVariable('Ptch','f4',('time'),fill_value=floatfill)
    varobj.units = "hundredths of degrees"
    varobj.long_name = "INST Pitch"
    varobj.epic_code = 1216
    #varobj.valid_range = [-18000, 18000] # physical limit, not sensor limit
    
    varobj = cdf.createVariable('Roll','f4',('time'),fill_value=floatfill)
    varobj.units = "hundredths of degrees"
    varobj.long_name = "INST Roll"
    varobj.epic_code = 1217
    #varobj.valid_range = [-18000, 18000] # physical limit, not sensor limit

    varobj = cdf.createVariable('HdgSTD','f4',('time'),fill_value=floatfill)
    varobj.units = "degrees"
    varobj.long_name = "Heading Standard Deviation"

    varobj = cdf.createVariable('PtchSTD','f4',('time'),fill_value=floatfill)
    varobj.units = "tenths of degrees"
    varobj.long_name = "Pitch Standard Deviation"

    varobj = cdf.createVariable('RollSTD','f4',('time'),fill_value=floatfill)
    varobj.units = "tenths of degrees"
    varobj.long_name = "Roll Standard Deviation"

    varobj = cdf.createVariable('Tx','f4',('time'),fill_value=floatfill)
    varobj.units = "hundredths of degrees"
    varobj.long_name = "ADCP Transducer Temperature"
    varobj.epic_code = 3017
    #varobj.valid_range = [-500, 4000]    

    varobj = cdf.createVariable('S','f4',('time'),fill_value=floatfill)
    varobj.units = "PPT"
    varobj.long_name = "SALINITY (PPT)"
    varobj.epic_code = 40
    #varobj.valid_range = [0, 40]    

    varobj = cdf.createVariable('xmitc','f4',('time'),fill_value=floatfill)
    varobj.units = "amps"
    varobj.long_name = "transmit current"

    varobj = cdf.createVariable('xmitv','f4',('time'),fill_value=floatfill)
    varobj.units = "volts"
    varobj.long_name = "transmit voltage"

    varobj = cdf.createVariable('Ambient_Temp','i2',('time'),fill_value=intfill)
    varobj.units = "C"
    varobj.long_name = "Ambient_Temp"

    varobj = cdf.createVariable('Pressure+','i2',('time'),fill_value=intfill)
    varobj.units = "unknown"
    varobj.long_name = "Pressure+"

    varobj = cdf.createVariable('Pressure-','i2',('time'),fill_value=intfill)
    varobj.units = "unknown"
    varobj.long_name = "Pressure-"

    varobj = cdf.createVariable('Attitude_Temp','i2',('time'),fill_value=intfill)
    varobj.units = "C"
    varobj.long_name = "Attitude_Temp"

    for i in range(4):
        varname = "EWD%d" % (i+1)
        varobj = cdf.createVariable(varname,'u2',('time'),fill_value=intfill)
        varobj.units = "binary flag"
        varobj.long_name = "Error Status Word %d" % (i+1)

    if  ensData['FLeader']['Depth_sensor_available'] == 'Yes':
        varobj = cdf.createVariable('Pressure','f4',('time'),fill_value=floatfill)
        varobj.units = "deca-pascals"
        varobj.long_name = "ADCP Transducer Pressure"
        varobj.epic_code = 4
        #varobj.valid_range = [0, maxfloat]
    
        varobj = cdf.createVariable('PressVar','f4',('time'),fill_value=floatfill)
        varobj.units = "deca-pascals"
        varobj.long_name = "ADCP Transducer Pressure Variance"
        #varobj.valid_range = [0, maxfloat]
    
    if 'BTData' in ensData: 
        # write globals attributable to BT setup
        cdf.setncattr('TRDI_BT_pings_per_ensemble',ensData['BTData']['Pings_per_ensemble'])
        cdf.setncattr('TRDI_BT_reacquire_delay',ensData['BTData']['delay_before_reacquire'])
        cdf.setncattr('TRDI_BT_min_corr_mag',ensData['BTData']['Corr_Mag_Min'])
        cdf.setncattr('TRDI_BT_min_eval_mag',ensData['BTData']['Eval_Amp_Min'])
        cdf.setncattr('TRDI_BT_min_percent_good',ensData['BTData']['PGd_Minimum'])
        cdf.setncattr('TRDI_BT_mode',ensData['BTData']['Mode'])
        cdf.setncattr('TRDI_BT_max_err_vel',ensData['BTData']['Err_Vel_Max'])
        #cdf.setncattr('TRDI_BT_max_tracking_depth',ensData['BTData'][''])
        #cdf.setncattr('TRDI_BT_shallow_water_gain',ensData['BTData'][''])

        for i in range(4):
            varname = "BTR%d" % (i+1)
            varobj = cdf.createVariable(varname,'u8',('time'),fill_value=intfill)
            varobj.units = "cm"
            varobj.long_name = "BT Range %d" % (i+1)
            #varobj.valid_range = [0, 65536*16777215]

        for i in range(4):
            varnames = ('BTWe','BTWu','BTWv','BTWd')
            longnames = ('BT Error Velocity','BT Eastward Velocity','BT Northward Velocity','BT Vertical Velocity')
            if ensData['FLeader']['Coord_Transform'] == 'EARTH':
                varobj = cdf.createVariable(varnames[i+1],'i2',('time'),fill_value=intfill)
                varobj.units = "mm s-1"
                varobj.long_name = "%s, mm s-1" % longnames[i+1]
                #varobj.valid_range = [-32768, 32767]
                
            else:
                varname = "BTV%d" % (i+1)
                varobj = cdf.createVariable(varname,'i2',('time'),fill_value=intfill)
                varobj.units = "mm s-1"
                varobj.long_name = "BT velocity, mm s-1 %d" % (i+1)
                #varobj.valid_range = [-32768, 32767]
                
        for i in range(4):
            varname = "BTc%d" % (i+1)
            varobj = cdf.createVariable(varname,'u2',('time'),fill_value=intfill)
            varobj.units = "counts"
            varobj.long_name = "BT correlation %d" % (i+1)
            #varobj.valid_range = [0, 255]
                
        for i in range(4):
            varname = "BTe%d" % (i+1)
            varobj = cdf.createVariable(varname,'u2',('time'),fill_value=intfill)
            varobj.units = "counts"
            varobj.long_name = "BT evaluation amplitude %d" % (i+1)
            #varobj.valid_range = [0, 255]
            
        for i in range(4):
            varname = "BTp%d" % (i+1)
            varobj = cdf.createVariable(varname,'u2',('time'),fill_value=intfill)
            varobj.units = "percent"
            varobj.long_name = "BT percent good %d" % (i+1)
            #varobj.valid_range = [0, 100]

        for i in range(4):
            varname = "BTRSSI%d" % (i+1)
            varobj = cdf.createVariable(varname,'u2',('time'),fill_value=intfill)
            varobj.units = "counts"
            varobj.long_name = "BT Receiver Signal Strength Indicator %d" % (i+1)
            #varobj.valid_range = [0, 255]

        if ensData['BTData']['Mode'] == 0: # water reference layer was used
            varobj = cdf.createVariable('BTRmin','f4',('time'),fill_value=floatfill)
            varobj.units = 'dm'
            varobj.long_name = "BT Ref. min"
            #varobj.valid_range = [0,999]
            varobj = cdf.createVariable('BTRnear','f4',('time'),fill_value=floatfill)
            varobj.units = 'dm'
            varobj.long_name = "BT Ref. near"
            #varobj.valid_range = [0,9999]
            varobj = cdf.createVariable('BTRfar','f4',('time'),fill_value=floatfill)
            varobj.units = 'dm'
            varobj.long_name = "BT Ref. far"
            #varobj.valid_range = [0,9999]
                
            for i in range(4):
                varname = "BTRv%d" % (i+1)
                varobj = cdf.createVariable(varname,'i2',('time'),fill_value=intfill)
                varobj.units = "mm s-1"
                varobj.long_name = "BT Ref. velocity, mm s-1 %d" % (i+1)
                #varobj.valid_range = [-32768, 32767]

            for i in range(4):
                varname = "BTRc%d" % (i+1)
                varobj = cdf.createVariable(varname,'u2',('time'),fill_value=intfill)
                varobj.units = "counts"
                varobj.long_name = "BT Ref. correlation %d" % (i+1)
                #varobj.valid_range = [0, 255]
                    
            for i in range(4):
                varname = "BTRi%d" % (i+1)
                varobj = cdf.createVariable(varname,'u2',('time'),fill_value=intfill)
                varobj.units = "counts"
                varobj.long_name = "BT Ref. intensity %d" % (i+1)
                #varobj.valid_range = [0, 255]
                
            for i in range(4):
                varname = "BTRp%d" % (i+1)
                varobj = cdf.createVariable(varname,'u2',('time'),fill_value=intfill)
                varobj.units = "percent"
                varobj.long_name = "BT Ref. percent good %d" % (i+1)
                varobj.epic_code = 1269+i
                #varobj.valid_range = [0, 100]
        
        
    
    if 'VPingSetup' in ensData:
        writeDict2atts(cdf, ensData['VPingSetup'], "TRDI_VBeam_")

    if 'VBeamLeader' in ensData:
        writeDict2atts(cdf, ensData['VBeamLeader'], "TRDI_VBeam_")

    if ('VBeamVData' in ensData):
        if ensData['VBeamLeader']['Vertical_Depth_Cells'] == ensData['FLeader']['Number_of_Cells']:
            varobj = cdf.createVariable("vel5",'f4',('time','depth'),fill_value=floatfill)
            varobj.units = "mm s-1"
            varobj.long_name = "Beam 5 velocity (mm s-1)"
            #varobj.valid_range = [-32767, 32767]
            varobj = cdf.createVariable("cor5",'u2',('time','depth'),fill_value=intfill)
            varobj.units = "counts"
            varobj.long_name = "Beam 5 correlation"
            #varobj.valid_range = [0, 255]
            varobj = cdf.createVariable("att5",'u2',('time','depth'),fill_value=intfill)
            varobj.units = "counts"
            varobj.long_name = "ADCP attenuation of beam 5"
            #varobj.valid_range = [0, 255]
            if ('VBeamGData' in ensData):
                varobj = cdf.createVariable("PGd5",'u2',('time','depth'),fill_value=intfill)
                varobj.units = "counts"
                varobj.long_name = "Percent Good Beam 5"
                #varobj.valid_range = [0, 100]
            else:
                cdf.TRDI_VBeam_note1 = 'Vertical beam data found without Percent Good'
        else:
            print('Vertical beam data found with different number of cells.')
            cdf.TRDI_VBeam_note = 'Vertical beam data found with different number of cells. Vertical beam data not exported to netCDF'
            print('Vertical beam data not exported to netCDF')

    if ('WaveParams' in ensData):
        # no units given for any of these in the TRDI docs
        varobj = cdf.createVariable("Hs",'f4',('time'),fill_value=floatfill)
        varobj.units = "m"
        varobj.long_name = "Significant Wave Height (m)"
        varobj = cdf.createVariable("Tp",'f4',('time'),fill_value=floatfill)
        varobj.units = "s"
        varobj.long_name = "Peak Wave Period (s)"
        varobj = cdf.createVariable("Dp",'f4',('time'),fill_value=floatfill)
        varobj.units = "Deg."
        varobj.long_name = "Peak Wave Direction (Deg.)"
        #varobj.valid_range = [0, 360]
        varobj = cdf.createVariable("Dm",'f4',('time'),fill_value=floatfill)
        varobj.units = "Deg."
        varobj.long_name = "Mea Peak Wave Direction (Deg.)"
        #varobj.valid_range = [0, 360]
        varobj = cdf.createVariable("SHmax",'f4',('time'),fill_value=floatfill)
        varobj.units = "m"
        varobj.long_name = "Maximum Wave Height (m)"
        varobj.note = "from zero crossing anaylsis of surface track time series"
        varobj = cdf.createVariable("SH13",'f4',('time'),fill_value=floatfill)
        varobj.units = "m"
        varobj.long_name = "Significant Wave Height of the largest 1/3 of the waves (m)"
        varobj.note = "in the field from zero crossing anaylsis of surface track time series"
        varobj = cdf.createVariable("SH10",'f4',('time'),fill_value=floatfill)
        varobj.units = "m"
        varobj.long_name = "Significant Wave Height of the largest 1/10 of the waves (m)"
        varobj.note = "in the field from zero crossing anaylsis of surface track time series"
        varobj = cdf.createVariable("STmax",'f4',('time'),fill_value=floatfill)
        varobj.units = "s"
        varobj.long_name = "Maximum Peak Wave Period (s)"
        varobj.note = "from zero crossing anaylsis of surface track time series"
        varobj = cdf.createVariable("ST13",'f4',('time'),fill_value=floatfill)
        varobj.units = "s"
        varobj.long_name = "Period associated with the peak wave height of the largest 1/3 of the waves (s)"
        varobj.note = "in the field from zero crossing anaylsis of surface track time series"
        varobj = cdf.createVariable("ST10",'f4',('time'),fill_value=floatfill)
        varobj.units = "s"
        varobj.long_name = "Period associated with the peak wave height of the largest 1/10 of the waves (s)"
        varobj.note = "in the field from zero crossing anaylsis of surface track time series"
        varobj = cdf.createVariable("T01",'f4',('time'),fill_value=floatfill)
        varobj.units = " " 
        varobj = cdf.createVariable("Tz",'f4',('time'),fill_value=floatfill)
        varobj.units = " "
        varobj = cdf.createVariable("Tinv1",'f4',('time'),fill_value=floatfill)
        varobj.units = " "
        varobj = cdf.createVariable("S0",'f4',('time'),fill_value=floatfill)
        varobj.units = " "
        varobj = cdf.createVariable("Source",'f4',('time'),fill_value=floatfill)
        varobj.units = " "

    if ('WaveSeaSwell' in ensData):
        # no units given for any of these in the TRDI docs
        varobj = cdf.createVariable("HsSea",'f4',('time'),fill_value=floatfill)
        varobj.units = "m"
        varobj.long_name = "Significant Wave Height (m)"
        varobj.note = "in the sea region of the power spectrum"
        varobj = cdf.createVariable("HsSwell",'f4',('time'),fill_value=floatfill)
        varobj.units = "m"
        varobj.long_name = "Significant Wave Height (m)"
        varobj.note = "in the swell region of the power spectrum"
        varobj = cdf.createVariable("TpSea",'f4',('time'),fill_value=floatfill)
        varobj.units = "s"
        varobj.long_name = "Peak Wave Period (s)"
        varobj.note = "in the sea region of the power spectrum"
        varobj = cdf.createVariable("TpSwell",'f4',('time'),fill_value=floatfill)
        varobj.units = "s"
        varobj.long_name = "Peak Wave Period (s)"
        varobj.note = "in the swell region of the power spectrum"
        varobj = cdf.createVariable("DpSea",'f4',('time'),fill_value=floatfill)
        varobj.units = "Deg."
        varobj.long_name = "Peak Wave Direction (Deg.)"
        varobj.note = "in the sea region of the power spectrum"
        varobj = cdf.createVariable("DpSwell",'f4',('time'),fill_value=floatfill)
        varobj.units = "Deg."
        varobj.long_name = "Peak Wave Direction (Deg.)"
        varobj.note = "in the swell region of the power spectrum"
        varobj = cdf.createVariable("SeaSwellPeriod",'f4',('time'),fill_value=floatfill)
        varobj.units = "s"
        varobj.long_name = "Transition Period between Sea and Swell (s)"     
        
    return cdf, cf_units

def bitstrLE(byte): # make a bit string from little endian byte
    # surely there's a better way to do this!!
    bits = ""
    for i in [7,6,5,4,3,2,1,0]:
        if (byte >> i) & 1:
            bits+="1"
        else:
            bits+="0"
    return bits

def bitstrBE(byte): # make a bit string from big endian byte
    # surely there's a better way to do this!!
    bits = ""
    for i in range(8): # Big Endian
        if (byte[0] >> i) & 1:
            bits+="1"
        else:
            bits+="0"
    return bits

# read header directly from a file pointer
# tests for end of file are here
def readTRDIHeader(infile):
    HeaderData = {}
    try:
        HeaderData['headerID'] = infile.read(1)
    except:
        return None
    try:
        HeaderData['sourceID'] = infile.read(1)
    except:
        return None
    try:
        HeaderData['nbytesperens'] = struct.unpack('<H',infile.read(2))[0]
    except:
        return None
    infile.read(1) # spare, skip it
    HeaderData['ndatatypes'] = infile.read(1)[0] # remember, bytes objects are arrays
    offsets = [0]*HeaderData['ndatatypes'] # predefine a list of ints to fill
    for i in range(HeaderData['ndatatypes']):
        offsets[i] = struct.unpack('<H',infile.read(2))[0]
        
    HeaderData['offsets'] = offsets

    return HeaderData

# read header data from an emsemble
def parseTRDIHeader(bstream):
    HeaderData = {}
    HeaderData['headerID'] = bstream[0] # byte 1
    HeaderData['sourceID'] = bstream[1] # byte 2
    HeaderData['nbytesperens'] = struct.unpack('<H',bstream[2:4])[0]
    # spare, skip it, byte 5
    HeaderData['ndatatypes'] = bstream[5] # byte 6
    offsets = [0]*HeaderData['ndatatypes'] # predefine a list of ints to fill
    for i in range(HeaderData['ndatatypes']):
        offsets[i] = struct.unpack('<H',bstream[6+i*2:6+i*2+2])[0]
        
    HeaderData['offsets'] = offsets

    return HeaderData
    
def parseTRDIFixedLeader(bstream, offset):
    # bstream is a bytes object that contains an entire ensemble
    # offset is the location in the bytes object of the first byte of this data format
    FLeaderData = {}
    leaderID = struct.unpack('<H',bstream[offset:offset+2])[0]
    if leaderID != 0:
        print("expected fixed leader ID, instead found %g",leaderID)
        return -1
    FLeaderData['CPU_Version'] = "%s.%s" % (bstream[offset+2],bstream[offset+4])
    
    FLeaderData['System_Configuration_LSB'] = bitstrLE(bstream[offset+4])
    # anyone who has a better way to convert these bits, pelase tell me!
    FLeaderData['System_Frequency'] = int(FLeaderData['System_Configuration_LSB'][5:8],2)
    SysFreqs = (75,150,300,600,1200,2400)  
    FLeaderData['System_Frequency'] = SysFreqs[FLeaderData['System_Frequency']]
    if FLeaderData['System_Configuration_LSB'][4] == "1":
        FLeaderData['Beam_Pattern'] = 'Convex'
    else:
        FLeaderData['Beam_Pattern'] = 'Concave'
    FLeaderData['Sensor_Configuration'] = int(FLeaderData['System_Configuration_LSB'][2:4],2) + 1
    if FLeaderData['System_Configuration_LSB'][1] == "1":
        FLeaderData['Transducer_Head_Is_Attached'] = 'Yes'
    else:
        FLeaderData['Transducer_Head_Is_Attached'] = 'No'
    if FLeaderData['System_Configuration_LSB'][0] == "1":
        FLeaderData['Orientation'] = 'Up-facing beams'
    else:
        FLeaderData['Orientation'] = 'Down-facing beams'
        
    FLeaderData['System_Configuration_MSB'] = bitstrLE(bstream[offset+5])
    FLeaderData['Beam_Angle'] = int(FLeaderData['System_Configuration_MSB'][5:8],2)
    # the angles 15, 20, and 30 are used by the Workhorse
    # the angle 25 is used by the Sentinel V, and so far, is always 25
    Angles = (15,20,30,0,0,0,0,25)  
    FLeaderData['Beam_Angle'] = Angles[FLeaderData['Beam_Angle']]
    FLeaderData['Beam_Configuration'] = int(FLeaderData['System_Configuration_MSB'][0:4],2)
    if FLeaderData['Beam_Configuration'] == 4:
        FLeaderData['Beam_Configuration'] = '4-bm janus'
    elif FLeaderData['Beam_Configuration'] == 5:
        FLeaderData['Beam_Configuration'] = '5-bm janus cfig demod'
    elif FLeaderData['Beam_Configuration'] == 15:
        FLeaderData['Beam_Configuration'] = '5-bm janus cfig (2 demd)'
    else: FLeaderData['Beam_Configuration'] = 'unknown'
    
    FLeaderData['Simulated_Data'] = bstream[offset+6]
    
    FLeaderData['Lag_Length'] = bstream[offset+7]
    FLeaderData['Number_of_Beams'] = bstream[offset+8]
    FLeaderData['Number_of_Cells'] = bstream[offset+9]
    FLeaderData['Pings_Per_Ensemble'] = struct.unpack('<h',bstream[offset+10:offset+12])[0]
    FLeaderData['Depth_Cell_Length_cm'] = struct.unpack('<h',bstream[offset+12:offset+14])[0]
    FLeaderData['Blank_after_Transmit_cm'] = struct.unpack('<h',bstream[offset+14:offset+16])[0]
    FLeaderData['Signal_Processing_Mode'] = bstream[offset+16]
    FLeaderData['Low_Corr_Threshold'] = bstream[offset+17]
    FLeaderData['No._Code_Reps'] = bstream[offset+18]
    FLeaderData['PGd_Minimum'] = bstream[offset+19]
    FLeaderData['Error_Velocity_Threshold'] = struct.unpack('<h',bstream[offset+20:offset+22])[0]
    # TODO ping group time needs to be formatted better
    FLeaderData['Time_Between_Ping Groups'] = "%03d:%02d:%02d" % (bstream[offset+22],bstream[offset+23],bstream[offset+24])

    FLeaderData['Coord_Transform_LSB'] = bitstrLE(bstream[offset+25])
    FLeaderData['Coord_Transform'] = int(FLeaderData['Coord_Transform_LSB'][3:5],2)
    Xforms = ('BEAM','INST','SHIP','EARTH')  
    FLeaderData['Coord_Transform'] = Xforms[FLeaderData['Coord_Transform']]
    if FLeaderData['Coord_Transform_LSB'][5] == '1':
        FLeaderData['Tilts_Used'] = 'Yes'
    else:
        FLeaderData['Tilts_Used'] = 'No'
    if FLeaderData['Coord_Transform_LSB'][6] == '1':
        FLeaderData['3-Beam_Solution_Used'] = 'Yes'
    else:
        FLeaderData['3-Beam_Solution_Used'] = 'No'
    if FLeaderData['Coord_Transform_LSB'][7] == '1':
        FLeaderData['Bin_Mapping_Used'] = 'Yes'
    else:
        FLeaderData['Bin_Mapping_Used'] = 'No'
        
    FLeaderData['Heading_Alignment_Hundredths_of_Deg'] = struct.unpack('<h',bstream[offset+26:offset+28])[0]
    FLeaderData['Heading_Bias_Hundredths_of_Deg'] = struct.unpack('<h',bstream[offset+28:offset+30])[0]
    
    FLeaderData['Sensor_Source_Byte'] = bitstrLE(bstream[offset+30])
    if FLeaderData['Sensor_Source_Byte'][1] == '1':
        FLeaderData['Calculate_EC_from_ED_ES_and_ET'] = 'Yes'
    else:
        FLeaderData['Calculate_EC_from_ED_ES_and_ET'] = 'No'
    if FLeaderData['Sensor_Source_Byte'][2] == '1':
        FLeaderData['Uses_ED_from_depth_sensor'] = 'Yes'
    else:
        FLeaderData['Uses_ED_from_depth_sensor'] = 'No'
    if FLeaderData['Sensor_Source_Byte'][3] == '1':
        FLeaderData['Uses_EH_from_transducer_heading_sensor'] = 'Yes'
    else:
        FLeaderData['Uses_EH_from_transducer_heading_sensor'] = 'No'
    if FLeaderData['Sensor_Source_Byte'][4] == '1':
        FLeaderData['Uses_EP_from_transducer_pitch_sensor'] = 'Yes'
    else:
        FLeaderData['Uses_EP_from_transducer_pitch sensor'] = 'No'
    if FLeaderData['Sensor_Source_Byte'][5] == '1':
        FLeaderData['Uses_ER_from_transducer_roll_sensor'] = 'Yes'
    else:
        FLeaderData['Uses_ER_from_transducer_roll_sensor'] = 'No'
    if FLeaderData['Sensor_Source_Byte'][6] == '1':
        FLeaderData['Uses_ES_from_conductivity_sensor'] = 'Yes'
    else:
        FLeaderData['Uses_ES_from_conductivity_sensor'] = 'No'
    if FLeaderData['Sensor_Source_Byte'][7] == '1':
        FLeaderData['Uses_ET_from_transducer_temperature_sensor'] = 'Yes'
    else:
        FLeaderData['Uses_ET_from_transducer_temperature_sensor'] = 'No'
        
    FLeaderData['Sensor_Avail_Byte'] = bitstrLE(bstream[offset+31])
    if FLeaderData['Sensor_Avail_Byte'][1] == '1':
        FLeaderData['Speed_of_sound_sensor_available'] = 'Yes'
    else:
        FLeaderData['Speed_of_sound_sensor_available'] = 'No'
    if FLeaderData['Sensor_Avail_Byte'][2] == '1':
        FLeaderData['Depth_sensor_available'] = 'Yes'
    else:
        FLeaderData['Depth_sensor_available'] = 'No'
    if FLeaderData['Sensor_Avail_Byte'][3] == '1':
        FLeaderData['Heading_sensor_available'] = 'Yes'
    else:
        FLeaderData['Heading_sensor_available'] = 'No'
    if FLeaderData['Sensor_Avail_Byte'][4] == '1':
        FLeaderData['Pitch_sensor_available'] = 'Yes'
    else:
        FLeaderData['Pitch_sensor_available'] = 'No'
    if FLeaderData['Sensor_Avail_Byte'][5] == '1':
        FLeaderData['Roll_sensor_available'] = 'Yes'
    else:
        FLeaderData['Roll_sensor_available'] = 'No'
    if FLeaderData['Sensor_Avail_Byte'][6] == '1':
        FLeaderData['Conductivity_sensor_available'] = 'Yes'
    else:
        FLeaderData['Conductivity_sensor_available'] = 'No'
    if FLeaderData['Sensor_Avail_Byte'][7] == '1':
        FLeaderData['Temperature_sensor_available'] = 'Yes'
    else:
        FLeaderData['Temperature_sensor_available'] = 'No'
        
    FLeaderData['Bin_1_distance_cm'] = struct.unpack('<h',bstream[offset+32:offset+34])[0]
    FLeaderData['Xmit_pulse_length_cm'] = struct.unpack('<h',bstream[offset+34:offset+36])[0]
    FLeaderData['Ref_Lyr_Avg_Starting_cell'] = bstream[offset+36]
    FLeaderData['Ref_Lyr_Avg_Ending_cell'] = bstream[offset+37]
    FLeaderData['False_Target_Threshold'] = bstream[offset+38]
    FLeaderData['Transmit_lag_distance_cm'] = struct.unpack('<h',bstream[offset+40:offset+42])[0]
    FLeaderData['CPU_Board_Serial_Number'] = ""
    for i in range(8):
        FLeaderData['CPU_Board_Serial_Number'] = FLeaderData['CPU_Board_Serial_Number'] + ("%x" % bstream[offset+42+i])
   
    FLeaderData['System_Bandwidth'] = struct.unpack('<h',bstream[offset+50:offset+52])[0]
    FLeaderData['System_Power'] = bstream[offset+52]
    FLeaderData['Base_Frequency_Index'] = bstream[offset+53]
    #TODO these two need to be interpreted as spare if WH ADCP
    #rawBytes, FLeaderData['Serial Number for Remus only'] = struct.unpack('<H',infile.read(2))[0]
    #FLeaderData['Beam Angle for H-ADCP only'] = "%g" % infile.read(1)[0]
    
    return FLeaderData

def parseTRDIVariableLeader(bstream, offset):
    # bstream is a bytes object that contains an entire ensemble
    # offset is the location in the bytes object of the first byte of this data format
    VLeaderData = {}
    leaderID = struct.unpack('<H',bstream[offset:offset+2])[0]
    if leaderID != 128:
        print("expected variable leader ID, instead found %g",leaderID)
        return -1
    VLeaderData['Ensemble_Number'] = struct.unpack('<H',bstream[offset+2:offset+4])[0]
    VLeaderData['Year'] = bstream[offset+4]
    if VLeaderData['Year'] < 50: # circa 2000
        VLeaderData['Year'] += 2000
    else:
        VLeaderData['Year'] += 1900
        
    VLeaderData['Month'] = bstream[offset+5]
    VLeaderData['Day'] = bstream[offset+6]
    VLeaderData['Hour'] = bstream[offset+7]
    VLeaderData['Minute'] = bstream[offset+8]
    VLeaderData['Second'] = bstream[offset+9]
    VLeaderData['Hundredths'] = bstream[offset+10]
    VLeaderData['Ensemble_#_MSB'] = bstream[offset+11]
    VLeaderData['Ensemble_Number'] = VLeaderData['Ensemble_Number']+(VLeaderData['Ensemble_#_MSB']<<16)

    VLeaderData['timestr'] = "%04d:%02d:%02d %02d:%02d:%02d.%03d" % (
        VLeaderData['Year'], VLeaderData['Month'],
        VLeaderData['Day'], VLeaderData['Hour'], VLeaderData['Minute'],
        VLeaderData['Second'], VLeaderData['Hundredths'])
    
    # compute time and time2
    jd = julian(VLeaderData['Year'],VLeaderData['Month'],VLeaderData['Day'],
                VLeaderData['Hour'],VLeaderData['Minute'],VLeaderData['Second'],
                VLeaderData['Hundredths'])        
    VLeaderData['dtobj'] = dt.datetime(VLeaderData['Year'], VLeaderData['Month'],
        VLeaderData['Day'], VLeaderData['Hour'], VLeaderData['Minute'],
        VLeaderData['Second'], VLeaderData['Hundredths']*10000)
    # centiseconds * 10000 = microseconds
    jddt = ajd(VLeaderData['dtobj'])
    VLeaderData['julian_day_from_as_datetime_object'] = jddt
    VLeaderData['julian_day_from_julian'] = jd
    #VLeaderData['time'] = jd
    VLeaderData['EPIC_time'] = int(math.floor(jd))
    VLeaderData['EPIC_time2'] = int((jd - math.floor(jd))*(24*3600*1000))
    
    VLeaderData['BIT_Result_Byte_13'] = bitstrLE(bstream[offset+12])
    VLeaderData['Demod_1_error_bit'] = int(VLeaderData['BIT_Result_Byte_13'][3])
    VLeaderData['Demod_0_error_bit'] = int(VLeaderData['BIT_Result_Byte_13'][4])
    VLeaderData['Timing_Card_error_bit'] = int(VLeaderData['BIT_Result_Byte_13'][6])

    VLeaderData['Speed_of_Sound'] = struct.unpack('<H',bstream[offset+14:offset+16])[0]
    VLeaderData['Depth_of_Transducer'] = struct.unpack('<H',bstream[offset+16:offset+18])[0]
    VLeaderData['Heading, Pitch, Roll units'] = "hundredths_of_a_degree"    
    VLeaderData['Heading'] = struct.unpack('<H',bstream[offset+18:offset+20])[0]
    VLeaderData['Pitch'] = struct.unpack('<h',bstream[offset+20:offset+22])[0]
    VLeaderData['Roll'] = struct.unpack('<h',bstream[offset+22:offset+24])[0]
    VLeaderData['Salinity'] = struct.unpack('<H',bstream[offset+24:offset+26])[0]
    VLeaderData['Temperature'] = struct.unpack('<H',bstream[offset+26:offset+28])[0]
    VLeaderData['MPT_minutes'] = bstream[offset+28]
    VLeaderData['MPT_seconds'] = bstream[offset+29]
    VLeaderData['MPT_hundredths'] = bstream[offset+30]
    VLeaderData['H/Hdg_Std_Dev'] = bstream[offset+31]
    VLeaderData['P/Pitch_Std_Dev'] = bstream[offset+32]
    VLeaderData['R/Roll_Std_Dev'] = bstream[offset+33]
    # the V Series PDO Output is different for the ADC channels        
    # V PD0 this is ADC Channel 0 not used    
    VLeaderData['Xmit_Current'] = bstream[offset+34] # ADC Channel 0
    # V PD0 this is ADC Channel 1 XMIT Voltage    
    VLeaderData['Xmit_Voltage'] = bstream[offset+35] # ADC Channel 1
    # V PD0 this is ADC Channel 2 not used    
    VLeaderData['Ambient_Temp'] = bstream[offset+36] #ADC Channel 2
    # V PD0 this is ADC Channel 3 not used    
    VLeaderData['Pressure_(+)'] = bstream[offset+37] #ADC Channel 3
    # V PD0 this is ADC Channel 4 not used    
    VLeaderData['Pressure_(-)'] = bstream[offset+38] #ADC Channel 4
    # V PD0 this is ADC Channel 5 not used    
    VLeaderData['Attitude_Temp'] = bstream[offset+39] #ADC Channel 5
    # V PD0 this is ADC Channel 6 not used    
    VLeaderData['Attitude'] = bstream[offset+40] #ADC Channel 6
    # V PD0 this is ADC Channel 7 not used    
    VLeaderData['Contamination_Sensor'] = bstream[offset+41] #ADC Channel 7

    VLeaderData['Error_Status_Word_Low_16_bits_LSB'] = bitstrLE(bstream[offset+42])
    VLeaderData['Bus_Error_exception'] = int(VLeaderData['Error_Status_Word_Low_16_bits_LSB'][7])
    VLeaderData['Address_Error_exception'] = int(VLeaderData['Error_Status_Word_Low_16_bits_LSB'][6])
    VLeaderData['Illegal_Instruction_exception'] = int(VLeaderData['Error_Status_Word_Low_16_bits_LSB'][5])
    VLeaderData['Zero_Divide_exception'] = int(VLeaderData['Error_Status_Word_Low_16_bits_LSB'][4])
    VLeaderData['Emulator_exception'] = int(VLeaderData['Error_Status_Word_Low_16_bits_LSB'][3])
    VLeaderData['Unassigned_exception'] = int(VLeaderData['Error_Status_Word_Low_16_bits_LSB'][2])
    VLeaderData['Watchdog_restart_occurred'] = int(VLeaderData['Error_Status_Word_Low_16_bits_LSB'][1])
    VLeaderData['Battery_Saver_power'] = int(VLeaderData['Error_Status_Word_Low_16_bits_LSB'][0])

    VLeaderData['Error_Status_Word_Low_16_bits_MSB'] = bitstrLE(bstream[offset+43])
    VLeaderData['Pinging'] = int(VLeaderData['Error_Status_Word_Low_16_bits_MSB'][7])
    VLeaderData['Cold_Wakeup_occurred'] = int(VLeaderData['Error_Status_Word_Low_16_bits_MSB'][1])
    VLeaderData['Unknown_Wakeup_occurred'] = int(VLeaderData['Error_Status_Word_Low_16_bits_MSB'][0])

    VLeaderData['Error_Status_Word_High_16_bits_LSB'] = bitstrLE(bstream[offset+44])
    VLeaderData['Clock_Read_error_occurred'] = int(VLeaderData['Error_Status_Word_High_16_bits_LSB'][7])
    VLeaderData['Unexpected_alarm'] = int(VLeaderData['Error_Status_Word_High_16_bits_LSB'][6])
    VLeaderData['Clock_jump_forward'] = int(VLeaderData['Error_Status_Word_High_16_bits_LSB'][5])
    VLeaderData['Clock_jump_backward'] = int(VLeaderData['Error_Status_Word_High_16_bits_LSB'][4])

    VLeaderData['Error_Status_Word_High_16_bits_MSB'] = bitstrLE(bstream[offset+42])
    VLeaderData['Power_Fail_(Unrecorded)'] = int(VLeaderData['Error_Status_Word_High_16_bits_MSB'][4])
    VLeaderData['Spurious_level_4_intr_(DSP)'] = int(VLeaderData['Error_Status_Word_High_16_bits_MSB'][3])
    VLeaderData['Spurious_level_5_intr_(UART)'] = int(VLeaderData['Error_Status_Word_High_16_bits_MSB'][2])
    VLeaderData['Spurious_level_6_intr_(CLOCK)'] = int(VLeaderData['Error_Status_Word_High_16_bits_MSB'][1])
    VLeaderData['Level_7_interrupt_occurred'] = int(VLeaderData['Error_Status_Word_High_16_bits_MSB'][0])

    # pressure of the water at the transducer head relative to one atmosphere (sea level)
    #VLeaderData['Pressure word byte 1'] = bitstrLE(bstream[offset+48])
    #VLeaderData['Pressure word byte 2'] = bitstrLE(bstream[offset+49])
    #VLeaderData['Pressure word byte 3'] = bitstrLE(bstream[offset+50])
    #VLeaderData['Pressure word byte 4'] = bitstrLE(bstream[offset+51])
    VLeaderData['Pressure_deca-pascals'] = bstream[offset+48]+(bstream[offset+49]<<8)+(bstream[offset+50]<<16)+(bstream[offset+51]<<24)
    VLeaderData['Pressure_variance_deca-pascals'] = bstream[offset+52]+(bstream[offset+53]<<8)+(bstream[offset+54]<<16)+(bstream[offset+55]<<24)

    VLeaderData['RTC_Century'] = bstream[offset+57]
    VLeaderData['RTC_Year'] = bstream[offset+58]
    VLeaderData['RTC_Month'] = bstream[offset+59]
    VLeaderData['RTC_Day'] = bstream[offset+60]
    VLeaderData['RTC_Hour'] = bstream[offset+61]
    VLeaderData['RTC_Minute'] = bstream[offset+62]
    VLeaderData['RTC_Second'] = bstream[offset+63]
    VLeaderData['RTC_Hundredths'] = bstream[offset+64]
    
    return VLeaderData

def parseTRDIVelocity(bstream, offset, ncells, nbeams):
    # each velocity value is stored as a two byte, twos complement integer 
    # [-32768 to 32767] with the LSB sent first.  Units are mm/s.
    # A value of -32768 = 0x8000 is a bad velocity value

    if bstream[offset+1] != 1:
        print("expected velocity ID, instead found %g",bstream[offset+1])
        return -1

    # start with a numpy array of bad values
    data = np.ones((nbeams,ncells),dtype=int) * -32768
    ibyte = 2
    for icell in range(ncells):
        for ibeam in range(nbeams):
            data[ibeam,icell] = struct.unpack('<h',bstream[offset+ibyte:offset+ibyte+2])[0]
            ibyte = ibyte+2   

    return data

def parseTRDICorrelation(bstream, offset, ncells, nbeams):
    if bstream[offset+1] != 2:
        print("expected correlation ID, instead found %g",bstream[offset+1])
        return -1
        
    # start with a numpy array of bad values
    data = np.ones((nbeams,ncells),dtype=int) * -32768
    ibyte = 2
    for icell in range(ncells):
        for ibeam in range(nbeams):
            data[ibeam,icell] = bstream[offset+ibyte]
            ibyte = ibyte+1

    return data    
    
def parseTRDIIntensity(bstream, offset, ncells, nbeams):
    if bstream[offset+1] != 3:
        print("expected intensity ID, instead found %g",bstream[offset+1])
        return -1
        
    # start with a numpy array of bad values
    data = np.ones((nbeams,ncells),dtype=int) * -32768
    ibyte = 2
    for icell in range(ncells):
        for ibeam in range(nbeams):
            data[ibeam,icell] = bstream[offset+ibyte]
            ibyte = ibyte+1

    return data    

def parseTRDIPercentGood(bstream, offset, ncells, nbeams):
    if bstream[offset+1] != 4:
        print("expected intensity ID, instead found %g",bstream[offset+1])
        return -1
        
    # start with a numpy array of bad values
    data = np.ones((nbeams,ncells),dtype=int) * -32768
    ibyte = 2
    for icell in range(ncells):
        for ibeam in range(nbeams):
            data[ibeam,icell] = bstream[offset+ibyte]
            ibyte = ibyte+1

    return data    

def parseTRDIxformMatrix(bstream, offset, nbeams):
    if bstream[offset+1] != 50: # \x00\x32
        print("expected transformation matrix ID, instead found %g",bstream[offset+1])
        return -1
        
    # start with a numpy array of bad values
    data = np.zeros((nbeams, 3), dtype=int)
    ibyte = 2
    
    for iaxis in range(3):
        for ibeam in range(nbeams):
            data[ibeam,iaxis] = struct.unpack('<h',bstream[offset+ibyte:offset+ibyte+2])[0]
            ibyte = ibyte+2
            
    return data

def parseTRDIVPingSetup(bstream, offset):
    # bstream is a bytes object that contains an entire ensemble
    # offset is the location in the bytes object of the first byte of this data format
    VPingSetupData = {}
    leaderID = struct.unpack('<H',bstream[offset:offset+2])[0]
    if leaderID != 28673: #\x70\x01 stored little endian
        print("expected V Series Ping Setup ID, instead found %g" % leaderID)
        return -1
    VPingSetupData['Ensemble_Interval_ms'] = bstream[offset+4]+(bstream[offset+5]<<8)+(bstream[offset+6]<<16)+(bstream[offset+7]<<24)
    VPingSetupData['Number_of_Pings'] = struct.unpack('<H',bstream[offset+10:offset+12])[0]
    VPingSetupData['Time_Between_Pings_ms'] = bstream[offset+10]+(bstream[offset+11]<<8)+(bstream[offset+12]<<16)+(bstream[offset+13]<<24)
    VPingSetupData['Offset_Between_Ping_Groups_ms'] = bstream[offset+14]+(bstream[offset+15]<<8)+(bstream[offset+16]<<16)+(bstream[offset+17]<<24)
    VPingSetupData['Ping_Sequence_Number'] = struct.unpack('<h',bstream[offset+22:offset+24])[0]
    VPingSetupData['Ambiguity_Velocity'] = struct.unpack('<h',bstream[offset+24:offset+26])[0]
    VPingSetupData['RX_Gain'] = bstream[offset+26]
    VPingSetupData['RX_Beam_Mask'] = bstream[offset+27]  
    VPingSetupData['TX_Beam_Mask'] = bstream[offset+28]  
    VPingSetupData['Ensemble_Offset'] = bstream[offset+30]+(bstream[offset+31]<<8)+(bstream[offset+32]<<16)+(bstream[offset+33]<<24)
    VPingSetupData['Ensemble_Count'] = bstream[offset+34]+(bstream[offset+35]<<8)

    VPingSetupData['Deployment_Start_Century'] = bstream[offset+36]
    VPingSetupData['Deployment_Start_Year'] = bstream[offset+37]
    VPingSetupData['Deployment_Start_Month'] = bstream[offset+38]
    VPingSetupData['Deployment_Start_Day'] = bstream[offset+39]
    VPingSetupData['Deployment_Start_Hour'] = bstream[offset+40]
    VPingSetupData['Deployment_Start_Minute'] = bstream[offset+41]
    VPingSetupData['Deployment_Start_Second'] = bstream[offset+42]
    VPingSetupData['Deployment_Start_Hundredths'] = bstream[offset+43]

    return VPingSetupData

def parseTRDIVSysConfig(bstream, offset):
    # bstream is a bytes object that contains an entire ensemble
    # offset is the location in the bytes object of the first byte of this data format
    VSysConfigData = {}
    leaderID = struct.unpack('<H',bstream[offset:offset+2])[0]
    if leaderID != 28672: #\x70\x00 stored little endian
        print("expected V Series System Config ID, instead found %g" % leaderID)
        return -1
    VSysConfigData['Firmware_Version'] = "%02d:%02d:%02d:%02d" % (bstream[offset+2],bstream[offset+3],bstream[offset+4],bstream[offset+5])
    VSysConfigData['System_Frequency'] = bstream[offset+6]+(bstream[offset+7]<<8)+(bstream[offset+8]<<16)+(bstream[offset+9]<<24)
    VSysConfigData['Pressure_Rating'] = struct.unpack('<H',bstream[offset+10:offset+12])[0]

    return VSysConfigData
    
def parseTRDIVBeamLeader(bstream, offset):
    # bstream is a bytes object that contains an entire ensemble
    # offset is the location in the bytes object of the first byte of this data format
    VBeamLeaderData = {}
    leaderID = struct.unpack('<H',bstream[offset:offset+2])[0]
    if leaderID != 3841: #\x0f\x01 stored little endian
        print("expected Vertical Beam Leader ID, instead found %g" % leaderID)
        return -1
    VBeamLeaderData['Vertical_Depth_Cells'] = struct.unpack('<H',bstream[offset+2:offset+4])[0]
    VBeamLeaderData['Vertical_Pings'] = struct.unpack('<H',bstream[offset+4:offset+6])[0]
    VBeamLeaderData['Vertical_Depth_Cell_Size_cm'] = struct.unpack('<H',bstream[offset+6:offset+8])[0]
    VBeamLeaderData['Vertical_First_Cell_Range_cm'] = struct.unpack('<H',bstream[offset+8:offset+10])[0]
    VBeamLeaderData['Vertical_Mode'] = struct.unpack('<H',bstream[offset+10:offset+12])[0]
    # 1 = low resolution slant beam cells = vertical beam cells
    # 2 = High resolution, dedicated surface tracking ping with 4:1 transmit/receive ratio or larger
    VBeamLeaderData['Vertical_Transmit_cm'] = struct.unpack('<H',bstream[offset+12:offset+14])[0]
    VBeamLeaderData['Vertical_Lag_Length_cm'] = struct.unpack('<H',bstream[offset+14:offset+16])[0]
    VBeamLeaderData['Transmit_Code_Elements'] = struct.unpack('<H',bstream[offset+16:offset+18])[0]
    VBeamLeaderData['Ping_Offset_Time'] = struct.unpack('<H',bstream[offset+30:offset+32])[0]

    return VBeamLeaderData
    
def parseTRDIVertVelocity(bstream, offset, ncells):
    leaderID = struct.unpack('<H',bstream[offset:offset+2])[0]
    if leaderID != 2560: #\x0a\x00 stored little endian
        print("expected Vertical Beam velocity ID, instead found %g" % leaderID)
        return -1

    # start with a numpy array of bad values
    data = np.ones((ncells),dtype=int) * -32768
    ibyte = 2
    for icell in range(ncells):
        data[icell] = struct.unpack('<h',bstream[offset+ibyte:offset+ibyte+2])[0]
        ibyte += 2   

    return data

def parseTRDIVertCorrelation(bstream, offset, ncells):
    leaderID = struct.unpack('<H',bstream[offset:offset+2])[0]
    if leaderID != 2816: #\x0b\x00 stored little endian
        print("expected Vertical Beam correlation ID, instead found %g" % leaderID)
        return -1

    # start with a numpy array of bad values
    data = np.ones((ncells),dtype=int) * -32768
    ibyte = 2
    for icell in range(ncells):
        data[icell] = bstream[offset+ibyte]
        ibyte += 1

    return data

def parseTRDIVertIntensity(bstream, offset, ncells):
    leaderID = struct.unpack('<H',bstream[offset:offset+2])[0]
    if leaderID != 3072: #\x0c\x00 stored little endian
        print("expected Vertical Beam intensity ID, instead found %g" % leaderID)
        return -1

    # start with a numpy array of bad values
    data = np.ones((ncells),dtype=int) * -32768
    ibyte = 2
    for icell in range(ncells):
        data[icell] = bstream[offset+ibyte]
        ibyte += 1

    return data

def parseTRDIVertPercentGood(bstream, offset, ncells):
    leaderID = struct.unpack('<H',bstream[offset:offset+2])[0]
    if leaderID != 3328: #\x0d\x00 stored little endian
        print("expected Vertical Beam percent good ID, instead found %g" % leaderID)
        return -1

    # start with a numpy array of bad values
    data = np.ones((ncells),dtype=int) * -32768
    ibyte = 2
    for icell in range(ncells):
        data[icell] = bstream[offset+ibyte]
        ibyte += 1

    return data
    
def parseTRDIVEventLog(bstream, offset):
    VEventLogData = {}
    leaderID = struct.unpack('<H',bstream[offset:offset+2])[0]
    if leaderID != 28676: #\x70\x04 stored little endian
        print("expected V Series Event Log ID, instead found %g" % leaderID)
        return -1
        
    VEventLogData['Fault_Count'] = struct.unpack('<H',bstream[offset+2:offset+4])[0]
    # TODO read the fault codes and output to a text file
        
    return VEventLogData

def parseTRDIWaveParameters(bstream, offset):
    # bstream is a bytes object that contains an entire ensemble
    # offset is the location in the bytes object of the first byte of this data format
    data = {}
    leaderID = struct.unpack('<H',bstream[offset:offset+2])[0]
    if leaderID != 11: #\x00\x0b stored little endian
        print("expected Wave Parameters ID, instead found %g" % leaderID)
        return -1
    data['Hs'] = struct.unpack('<H',bstream[offset+2:offset+4])[0]
    data['Tp'] = struct.unpack('<H',bstream[offset+4:offset+6])[0]
    data['Dp'] = struct.unpack('<H',bstream[offset+6:offset+8])[0]
    data['Dm'] = struct.unpack('<H',bstream[offset+16:offset+18])[0]
    data['SHmax'] = struct.unpack('<H',bstream[offset+30:offset+32])[0]
    data['SH13'] = struct.unpack('<H',bstream[offset+32:offset+34])[0]
    data['SH10'] = struct.unpack('<H',bstream[offset+34:offset+36])[0]
    data['STmax'] = struct.unpack('<H',bstream[offset+36:offset+38])[0]
    data['ST13'] = struct.unpack('<H',bstream[offset+38:offset+40])[0]
    data['ST10'] = struct.unpack('<H',bstream[offset+40:offset+42])[0]
    data['T01'] = struct.unpack('<H',bstream[offset+42:offset+44])[0]
    data['Tz'] = struct.unpack('<H',bstream[offset+44:offset+46])[0]
    data['Tinv1'] = struct.unpack('<H',bstream[offset+46:offset+48])[0]
    data['S0'] = struct.unpack('<H',bstream[offset+48:offset+50])[0]
    data['Source'] = bstream[offset+52]

    return data
     
def parseTRDIWaveSeaSwell(bstream, offset):
    # bstream is a bytes object that contains an entire ensemble
    # offset is the location in the bytes object of the first byte of this data format
    data = {}
    leaderID = struct.unpack('<H',bstream[offset:offset+2])[0]
    if leaderID != 12: #\x00\x0c stored little endian
        print("expected Wave Sea and Swell ID, instead found %g" % leaderID)
        return -1
    data['HsSea'] = struct.unpack('<H',bstream[offset+2:offset+4])[0]
    data['HsSwell'] = struct.unpack('<H',bstream[offset+4:offset+6])[0]
    data['TpSea'] = struct.unpack('<H',bstream[offset+6:offset+8])[0]
    data['TpSwell'] = struct.unpack('<H',bstream[offset+8:offset+10])[0]
    data['DpSea'] = struct.unpack('<H',bstream[offset+10:offset+12])[0]
    data['DpSwell'] = struct.unpack('<H',bstream[offset+12:offset+14])[0]
    data['SeaSwellPeriod'] = struct.unpack('<H',bstream[offset+44:offset+46])[0]

    return data

def parseTRDIBottomTrack(bstream, offset, nbeams):
    # bstream is a bytes object that contains an entire ensemble
    # offset is the location in the bytes object of the first byte of this data format
    data = {}
    leaderID = struct.unpack('<H',bstream[offset:offset+2])[0]
    if leaderID != 1536: #\x00\x06 stored little endian
        print("expected Bottom Track ID, instead found %g" % leaderID)
        return -1
    data['Pings_per_ensemble'] = struct.unpack('<H',bstream[offset+2:offset+4])[0]
    data['delay_before_reacquire'] = struct.unpack('<H',bstream[offset+4:offset+6])[0]
    data['Corr_Mag_Min'] = bstream[offset+6]
    data['Eval_Amp_Min'] = bstream[offset+7]
    data['PGd_Minimum'] = bstream[offset+8]
    data['Mode'] = bstream[offset+9]
    data['Err_Vel_Max'] = struct.unpack('<H',bstream[offset+10:offset+12])[0]
    data['BT_Range_LSB'] = np.ones((nbeams),dtype=int) * -32768
    ibyte = 16
    for ibeam in range(nbeams):
        data['BT_Range_LSB'][ibeam] = struct.unpack('<h',bstream[offset+ibyte:offset+ibyte+2])[0]
        ibyte = ibyte+2
    # the meaning and direction depends on the coordinate system used
    data['BT_Vel'] = np.ones((nbeams),dtype=float) * 1e35
    ibyte = 24
    for ibeam in range(nbeams):
        data['BT_Vel'][ibeam] = struct.unpack('<h',bstream[offset+ibyte:offset+ibyte+2])[0]
        ibyte = ibyte+2
    data['BT_Corr'] = np.ones((nbeams),dtype=int) * -32768
    ibyte = 32
    for ibeam in range(nbeams):
        data['BT_Corr'][ibeam] = bstream[offset+ibyte]
        ibyte = ibyte+1
    data['BT_Amp'] = np.ones((nbeams),dtype=int) * -32768
    ibyte = 36
    for ibeam in range(nbeams):
        data['BT_Amp'][ibeam] = bstream[offset+ibyte]
        ibyte = ibyte+1
    data['BT_PGd'] = np.ones((nbeams),dtype=int) * -32768
    ibyte = 40
    for ibeam in range(nbeams):
        data['BT_PGd'][ibeam] = bstream[offset+ibyte]
        ibyte = ibyte+1
    data['Ref_Layer_Min'] = struct.unpack('<H', bstream[offset+44:offset+46])[0]
    data['Ref_Layer_Near'] = struct.unpack('<H', bstream[offset+46:offset+48])[0]
    data['Ref_Layer_Far'] = struct.unpack('<H', bstream[offset+48:offset+50])[0]
    data['Ref_Layer_Vel'] = np.ones((nbeams), dtype=float) * 1e35
    ibyte = 50
    for ibeam in range(nbeams):
        data['Ref_Layer_Vel'][ibeam] = struct.unpack('<h', bstream[offset+ibyte:offset+ibyte+2])[0]
        ibyte = ibyte+2
    data['Ref_Layer_Corr'] = np.ones(nbeams, dtype=int) * -32768
    ibyte = 58
    for ibeam in range(nbeams):
        data['Ref_Layer_Corr'][ibeam] = bstream[offset+ibyte]
        ibyte = ibyte+1
    data['Ref_Layer_Amp'] = np.ones(nbeams, dtype=int) * -32768
    ibyte = 62
    for ibeam in range(nbeams):
        data['Ref_Layer_Amp'][ibeam] = bstream[offset+ibyte]
        ibyte = ibyte+1
    data['Ref_Layer_PGd'] = np.ones(nbeams, dtype=int) * -32768
    ibyte = 66
    for ibeam in range(nbeams):
        data['Ref_Layer_PGd'][ibeam] = bstream[offset+ibyte]
        ibyte = ibyte+1
    data['BT_Max_Depth'] = struct.unpack('<H', bstream[offset+70:offset+72])[0]
    data['RSSI_Amp'] = np.ones(nbeams, dtype=int) * -32768
    ibyte = 72
    for ibeam in range(nbeams):
        data['RSSI_Amp'][ibeam] = bstream[offset+ibyte]
        ibyte = ibyte+1
    data['GAIN'] = bstream[offset+76]
    data['BT_Range_MSB'] = np.ones(nbeams, dtype=int) * -32768
    ibyte = 77
    for ibeam in range(nbeams):
        data['BT_Range_MSB'][ibeam] = bstream[offset+ibyte]
        ibyte = ibyte+1
    data['BT_Range'] = np.ones(nbeams, dtype=int) * -32768
    for ibeam in range(nbeams):
        data['BT_Range'][ibeam] = data['BT_Range_LSB'][ibeam]+(data['BT_Range_MSB'][ibeam]<<16)

    return data


def __computeChecksum(ensemble):
    """Compute a checksum from header, length, and ensemble"""
    cs = 0    
    for byte in range(len(ensemble)-2):
        cs += ensemble[byte]
    return cs & 0xffff


# these date conversion functions came from
# http://stackoverflow.com/questions/31142181/calculating-julian-date-in-python/41769526#41769526
# AND OLD rps matlab code
def julian(year,month,day,hour,mn,sec,hund):
    # from julian.m and hms2h.m
    # convert hours, minutes and seconds to decimal hours
    decimalsec = sec+hund/100
    decimalhrs = hour+mn/60+decimalsec/3600
    mo=month+9
    yr=year-1
    
    if month > 2:
        mo -= 3
        yr = year
        
    c = math.floor(yr/100)
    yr = yr - c*100
    d = day
    j = math.floor((146097*c)/4)+math.floor((1461*yr)/4)+ \
        math.floor((153*mo +2)/5)+d+1721119

    # If you want julian days to start and end at noon, 
    # replace the following line with:
    # j=j+(decimalhrss-12)/24;
    j=j+decimalhrs/24
    
    return j

# moved to EPICstuff.py
"""
def jdn(dto):
    # Given datetime object returns Julian Day Number
    year = dto.year
    month = dto.month
    day = dto.day

    not_march = month < 3
    if not_march:
        year -= 1
        month += 12

    fr_y = math.floor(year / 100)
    reform = 2 - fr_y + math.floor(fr_y / 4)
    jjs = day + (
        math.floor(365.25 * (year + 4716)) + math.floor(30.6001 * (month + 1)) + reform - 1524)

    return jjs

def ajd(dto):
    #Given datetime object returns Astronomical Julian Day.
    #Day is from midnight 00:00:00+00:00 with day fractional
    #value added.
    jdd = jdn(dto)
    day_fraction = dto.hour / 24.0 + dto.minute / 1440.0 + dto.second / 86400.0
    return jdd + day_fraction - 0.5

def cftime2EPICtime(timecount, timeunits):
    # take a CF time variable and convert to EPIC time and time2
    # timecountis the integer count of minutes (for instance) since the time stamp
    # given in timeunits
    buf = timeunits.split()
    t0 = dt.datetime.strptime(buf[2]+' '+buf[3], '%Y-%m-%d %H:%M:%S.%f')
    t0j = ajd(t0)
    # julian day for EPIC is the beginning of the day e.g. midnight
    t0j = t0j+0.5 # add 0.5 because ajd() subtracts 0.5 
    
    if buf[0] == 'hours':
        tj = timecount/(24)
    elif buf[0] == 'minutes':
        tj = timecount/(24*60)
    elif buf[0] == 'seconds':
        tj = timecount/(24*60*60)
    elif buf[0] == 'milliseconds':
        tj = timecount/(24*60*60*1000)
    elif buf[0] == 'microseconds':
        tj = timecount/(24*60*60*1000*1000)
        
    tj = t0j+tj
    
    time = np.floor(tj)
    time2 = np.floor((tj-time)*(24*3600*1000))
    
    return time, time2
"""
"""
# this does not work
def cf2EPICtime(cftime, cfunits, cfcalendar):
    tobj = num2date(cftime,cfunits,calendar=cfcalendar)
    elapsed_sec = []
    for idx in range(len(tobj)):
        tdelta = tobj[idx]-tobj[0] # timedelta
        elapsed_sec.append(tdelta.total_seconds())
    # from the datetime object convert to time and time2
    jd = []
    time = []
    time2 = []
    for idx in range(len(tobj)):
        j = julian(tobj[idx].year,tobj[idx].month,tobj[idx].day, \
                   tobj[idx].hour,tobj[idx].minute,tobj[idx].second,\
                   floor(tobj[idx].microsecond/1000000))
        jd.append(j)
        time.append(int(floor(j)))
        time2.append(int((j - floor(j))*(24*3600*1000)))
        
    return time, time2
"""

    
def analyzepd0file(pd0File, verbose=False):
    # determine the input file size
    # read some ensembles, make an estimate of the number of ensembles within
    infile = open(pd0File, 'rb')
    
    while infile.tell() < 3000:
        b1 = infile.read(1)
        if b1 == b'\x7f':
            b2 = infile.read(1)
            if b2 == b'\x7f':
                break
    else:
        print('Desired TRDI 7f7f ID not found within 3 kB from beginning of the file')
        infile.close()
        sys.exit(1)
    
    startofdata = infile.tell()-2
    if startofdata != 0:
        print('data starts %d bytes into the file' % startofdata)
    infile.seek(startofdata)
    
    # need to read the header from the file to know the ensemble size
    Header = readTRDIHeader(infile)
    
    if Header['sourceID'] != b'\x7f':
        print('error - this is not a currents file')
        infile.close()
        
    # number of bytes per ensemble in the header does not include the checksum
    ensLen = Header['nbytesperens']+2
    print('ensemble length = %g' % ensLen)
    print(Header)
    # it is faster to define the netCDF file with a known length
    # for this we need to estimate how many ensembles we will be reading
    # for some reason, sys.getsizeof(infile) does not report the true length 
    # of the input file, so we will go to the end and see how far we have gone
    # there is a problem though.  While TRDI's documentation says the V Series
    # System Configuration data is always sent, this is not the case, so reading
    # only the first ensemble will not give the ensemble size typical over the
    # entire file
    # rewind and read this several ensembles because further in the ensemble
    # length can change on files output from Velocity
    infile.seek(startofdata)
    
    nens2check = 5
    nbytesperens = [0 for i in range(nens2check)]
    ndatatypes = [0 for i in range(nens2check)]
    
    for i in range(nens2check):
        fileposn = infile.tell()
        Header = readTRDIHeader(infile)
        ensLen = Header['nbytesperens']+2
        infile.seek(fileposn)
        ensData, ensError = parseTRDIensemble(infile.read(ensLen), verbose)
        if ensError != 'None':
            print('problem reading the first ensemble: ' + ensError)
            # infile.close()
            # sys.exit(1)
        
        if i == 0:
            firstEnsData = ensData
        print('ensemble %d has %d bytes and %d datatypes' % (ensData['VLeader']['Ensemble_Number'], 
            ensData['Header']['nbytesperens'], ensData['Header']['ndatatypes']))
        nbytesperens[i] = ensData['Header']['nbytesperens']+2
        ndatatypes[i] = ensData['Header']['ndatatypes']
        
    # the guess here is that if the first two ensembles are not the same,
    # it's the second ensemble that is representative of the data
    if nbytesperens[0] != nbytesperens[1]:
        ensLen = nbytesperens[1]
    else:
        ensLen = nbytesperens[0]      
    
    infile.seek(0,2)
    nbytesinfile = infile.tell()
    maxens = (nbytesinfile/ensLen)-1
    print('estimating %g ensembles in file using a %d ensemble size' % (maxens, ensLen))
        
    infile.close()
    
    print(ensData['Header'])
    print('ensemble length = %g' % ensLen)
    print('estimating %g ensembles in file' % maxens)
    
    # return maxens, ensLen, ensData, startofdata
    return maxens, ensLen, firstEnsData, startofdata


def __main():

    # TODO add - and -- types of command line arguments
    print('%s running on python %s' % (sys.argv[0], sys.version))

    if len(sys.argv) < 2:
        print("%s usage:" % sys.argv[0])
        print("TRDIpd0tonetcdf infilename outfilename [startingensemble endingensemble] " )
        print("[serialnum] [timetype] [delta_t_to_use] " )
        sys.exit(1)
    
    try:
        infileName = sys.argv[1]
    except:
        print('error - input file name missing')
        sys.exit(1)
        
    try:
        outfileName = sys.argv[2]
    except:
        print('error - output file name missing')
        sys.exit(1)
        
    print('Converting %s to %s' % (infileName, outfileName))

    try:
        goodens = [int(sys.argv[3]), int(sys.argv[4])]
    except:
        print('No starting and ending ensembles specified, processing entire file')
        goodens = [0, -1]
        
    try:
        serialnum = sys.argv[5]
    except:
        print('No serial number provided')
        serialnum = "unknown"      
    
    try:
        timetype = sys.argv[6]
    except:
        print('Time type will be CF')
        timetype = "CF"      
    
    try:
        delta_t_to_use = sys.argv[7]
    except:
        print('delta_t_to_use will be NONE')
        delta_t_to_use = "NONE"

    print('Start file conversion at ', dt.datetime.now())
    dopd0file(infileName, outfileName, goodens, serialnum, timetype, delta_t_to_use)
    
    print('Finished file conversion at ', dt.datetime.now())

    
if __name__ == "__main__":
    __main()
