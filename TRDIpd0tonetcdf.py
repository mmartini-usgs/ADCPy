"""
This code takes the raw currents portion of the adcp data from the splitter
program pd0.py and outputs raw current profile data to a netCDF4 file.
If you have a file with wave packets data, the splitter must be run first.

As a script:

python TRDIpd0tonetcdf.py [path] pd0File cdfFile

where:
    path         is a path to prepend to the following
    pd0File      is path of raw PD0 format input file with current ensembles
    cdfFile      is path of a netcdf4 EPIC compliant output file
    start        ensemble at which to start exporting
    end          ensemble at which to stop exporting
    
As a module:
import TRDIpd0tonetcdf as pd0

Notes:
    time and time2, the EPIC convention for netCDF, is not used here so that
    the resulting very large files generated can be reduced using existing 
    python too ls such as xarrays

Programmed according to the TRDI Workhorse Commands and Output Data Format document, March 2005
"""

# 1/25/2017 MM got this running on old Workhorse ADCP data

import sys, struct, math
import numpy as np 
from netCDF4 import Dataset
import datetime as dt

def dopd0file(pd0File, cdfFile, goodens):
    
    maxens, ensLen, ensData = analyzepd0file(pd0File)
    
    infile = open(pd0File, 'rb')
    
    if goodens[1] == np.inf:
        goodens[1] = maxens
           
    # we are good to go, get the output file ready
    print('Setting up netCDF file %s' % cdfFile)
    cdf = setupCdf(cdfFile, ensData, goodens)

    cdfIdx = 0
    ensCount = 0
    verbose = 0 # diagnostic, 1 = turn on output, 0 = silent
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
        #print('-- ensemble %d length %g, file position %g' % (ensCount, len(ens), infile.tell()))
        #print(ensData['Header'])        
        ensData, ensError = parseTRDIensemble(ens, verbose)
        
        if (ensError == 'None') and (ensCount >= goodens[0]):
            # write to netCDF
            if cdfIdx == 0:
                print('--- first ensembles read at %s and TRDI #%d' % (              
                    ensData['VLeader']['timestr'], ensData['VLeader']['Ensemble_Number']))
                
            varobj = cdf.variables['rec']
            try:
                varobj[cdfIdx] = ensData['VLeader']['Ensemble_Number']
            except:
                # here we have reached the end of the netCDF file
                cdf.close()
                infile.close()
                return

            # time calculations done when vleader is read
            varobj = cdf.variables['time']
            varobj[cdfIdx] = ensData['VLeader']['julian_day_from_julian']
            # diagnostic
            if (goodens[1]-goodens[0]-1)<100:
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
                varname = "AGC%d" % (i+1)
                varobj = cdf.variables[varname]
                varobj[cdfIdx,:] = ensData['IData'][i,:]

            if ('GData' in ensData):
                for i in range(nslantbeams):
                    varname = "PGd%d" % (i+1)
                    varobj = cdf.variables[varname]
                    varobj[cdfIdx,:] = ensData['GData'][i,:]

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
            varobj = cdf.variables['dac']
            varobj[cdfIdx] = ensData['VLeader']['Ambient_Temp']
            varobj = cdf.variables['VDD3']
            varobj[cdfIdx] = ensData['VLeader']['Pressure_(+)']
            varobj = cdf.variables['VDD1']
            varobj[cdfIdx] = ensData['VLeader']['Pressure_(-)']
            varobj = cdf.variables['VDC']
            varobj[cdfIdx] = ensData['VLeader']['Attitude_Temp']
            varobj = cdf.variables['EWD1']
            varobj[cdfIdx] = int(ensData['VLeader']['Error_Status_Word_Low_16_bits_LSB'])
            varobj = cdf.variables['EWD2']
            varobj[cdfIdx] = int(ensData['VLeader']['Error_Status_Word_Low_16_bits_MSB'])
            varobj = cdf.variables['EWD3']
            varobj[cdfIdx] = int(ensData['VLeader']['Error_Status_Word_High_16_bits_LSB'])
            varobj = cdf.variables['EWD4']
            varobj[cdfIdx] = int(ensData['VLeader']['Error_Status_Word_High_16_bits_MSB'])
            varobj = cdf.variables['Pressure']
            varobj[cdfIdx] = ensData['VLeader']['Pressure_deca-pascals']
            varobj = cdf.variables['PressVar']
            varobj[cdfIdx] = ensData['VLeader']['Pressure_variance_deca-pascals']

            if ('VBeamVData' in ensData):
                if ensData['VBeamLeader']['Vertical_Depth_Cells'] == ensData['FLeader']['Number_of_Cells']:
                    varobj = cdf.variables['vel5']
                    varobj[cdfIdx,:] = ensData['VBeamVData']
                    varobj = cdf.variables['cor5']
                    varobj[cdfIdx,:] = ensData['VBeamCData']
                    varobj = cdf.variables['AGC5']
                    varobj[cdfIdx,:] = ensData['VBeamIData']
                    if ('VBeamGData' in ensData):
                        varobj = cdf.variables['PGd5']
                        varobj[cdfIdx,:] = ensData['VBeamGData']
                    
            if ('WaveParams' in ensData):
                # we can get away with this because the key names and var names are the same
                for key, value in ensData['WaveParams']:
                    varobj = cdf.variables[key]
                    varobj[cdfIdx] = ensData['WaveParams'][key]
    
            if ('WaveSeaSwell' in ensData):
                # we can get away with this because the key names and var names are the same
                for key, value in ensData['WaveSeaSwell']:
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
            print('stopping at estimated end of file ensemble %d' % goodens[1])
            break        
        
        if maxens < 100:  n=10
        elif (maxens > 100) and (maxens < 1000): n=100
        elif (maxens > 1000) and (maxens < 10000): n=1000
        elif (maxens > 10000) and (maxens < 100000): n=10000
        elif (maxens > 100000) and (maxens < 1000000): n=100000
        else: n = 1000000
        
        ensf, ensi = math.modf(ensCount/n)
        if ensf == 0:
            #print('%d ensembles read at %s and TRDI #%d' % (ensCount,              
            #    ensData['VLeader']['timestr'], ensData['VLeader']['Ensemble_Number']))
            print('%d ensembles read at %s and TRDI #%d' % (ensCount,              
                ensData['VLeader']['dtobj'], ensData['VLeader']['Ensemble_Number']))
        
        if ensCount >= goodens[1]-1:
            print('stopping at requested ensemble %d' % goodens[1])
            break
        
        # note that ensemble lengths can change in the middle of the file!
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
        
    else:
        print('end of file reached')
        
    if ensCount < maxens:
        print('end of file reached after %d ensembles, less than estimated in the file' % ensCount)
    elif ensCount > maxens:
        print('end of file reached after %d ensembles, more than estimated in the file' % ensCount)
    
    infile.close()
    cdf.close()
    
    print('%d ensembles read, %d records written' % (ensCount, cdfIdx))
    
    
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
        elif val == 256: #raw == b'\x00\x01': 256
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
        elif val == 1792: #raw == b'\x00\x07':
            # this not defined in TRDI docs
            donothing()
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
    
def setupCdf(fname, ensData, gens):
    
    # TODO - detect BT data
    
    # note that 
    # f4 = 4 byte, 32 bit float
    maxfloat = 3.402823*10**38;
    intfill = -32768
    floatfill = 1E35
    
    nens = gens[1]-gens[0]-1
    print('creating netCDF file %s with %d records' % (fname, nens))
    
    cdf = Dataset(fname, "w", clobber=True, format="NETCDF4")
    
    # dimensions, in EPIC order
    cdf.createDimension('rec',nens)
    cdf.createDimension('depth',ensData['FLeader']['Number_of_Cells'])
    cdf.createDimension('lat',1)
    cdf.createDimension('lon',1)
    
    # write global attributes
    cdf.history = "translated to netCDF by adcpcurrents2cdf.py"
    
    writeDict2atts(cdf, ensData['FLeader'], "TRDI_")
    
    varobj = cdf.createVariable('rec','u4',('rec'),fill_value=intfill)
    varobj.units = "count"
    varobj.long_name = "Ensemble Number"
    # the ensemble number is a two byte LSB and a one byte MSB (for the rollover)
    varobj.valid_range = [0, 2**23]

    # this is a CF convention, and if f8, 64 bit is not used, time is clipped
    # for ADCP fast sampled, single ping data, need millisecond resolution
    varobj = cdf.createVariable('time','f8',('rec'))
    varobj.units = "milliseconds since 1986-5-23 00:00:00.0 0:00"
    varobj.NOTE = "Decimal Julian day [days] = time [days] + ( time2 [msec] / 86400000 [msec/day] )"    

    # we are not using these EPIC definitions yet.  They are here for reference
    #varobj = cdf.createVariable('rec','int',('time'))
    #varobj = cdf.createVariable('time','int',('rec'))
    #varobj.units = "True Julian Day"
    #varobj.epic_code = 624
    #varobj.datum = "Time (UTC) in True Julian Days: 2440000 = 0000 h on May 23, 1968"
    #varobj.NOTE = "Decimal Julian day [days] = time [days] + ( time2 [msec] / 86400000 [msec/day] )"    
    #varobj = cdf.createVariable('time2','int',('rec'))
    #varobj.units = "msec since 0:00 GMT"
    #varobj.epic_code = 624
    #varobj.datum = "Time (UTC) in True Julian Days: 2440000 = 0000 h on May 23, 1968"
    #varobj.NOTE = "Decimal Julian day [days] = time [days] + ( time2 [msec] / 86400000 [msec/day] )"    

    varobj = cdf.createVariable('sv','f4',('rec'),fill_value=floatfill)
    varobj.units = "m s-1"
    varobj.long_name = "sound velocity (m s-1)"
    varobj.valid_range = [1400, 1600]
    
    for i in range(4):
        varname = "vel%d" % (i+1)
        varobj = cdf.createVariable(varname,'f4',('rec','depth'),fill_value=floatfill)
        varobj.units = "mm s-1"
        varobj.long_name = "Beam %d velocity (mm s-1)" % (i+1)
        varobj.epic_code = 1280+i
        varobj.valid_range = [-32767, 32767]
    
    for i in range(4):
        varname = "cor%d" % (i+1)
        varobj = cdf.createVariable(varname,'u2',('rec','depth'),fill_value=intfill)
        varobj.units = "counts"
        varobj.long_name = "Beam %d correlation" % (i+1)
        varobj.epic_code = 1294+i
        varobj.valid_range = [0, 255]

    for i in range(4):
        varname = "AGC%d" % (i+1)
        varobj = cdf.createVariable(varname,'u2',('rec','depth'),fill_value=intfill)
        varobj.units = "counts"
        varobj.epic_code = 1221+i
        varobj.long_name = "Echo Intensity (AGC) Beam %d" % (i+1)
        varobj.valid_range = [0, 255]

    if ('GData' in ensData):
        for i in range(4):
            varname = "PGd%d" % (i+1)
            varobj = cdf.createVariable(varname,'u2',('rec','depth'),fill_value=intfill)
            varobj.units = "counts"
            varobj.long_name = "Percent Good Beam %d" % (i+1)
            varobj.epic_code = 1241+i
            varobj.valid_range = [0, 100]

    varobj = cdf.createVariable('Hdg','f4',('rec'),fill_value=floatfill)
    varobj.units = "hundredths of degrees"
    varobj.long_name = "INST Heading"
    varobj.epic_code = 1215
    varobj.heading_alignment = ensData['FLeader']['Heading_Alignment_Hundredths_of_Deg.']
    varobj.heading_bias = ensData['FLeader']['Heading_Bias_Hundredths_of_Deg.']
    varobj.valid_range = [0, 36000]
    if ensData['FLeader']['Heading_Bias_Hundredths_of_Deg.'] == 0:
        varobj.NOTE_9 = "no heading bias was applied by EB during deployment or by wavesmon"
    else:
        varobj.NOTE_9 = "a heading bias was applied by EB during deployment or by wavesmon"
    
    varobj = cdf.createVariable('Ptch','f4',('rec'),fill_value=floatfill)
    varobj.units = "hundredths of degrees"
    varobj.long_name = "INST Pitch"
    varobj.epic_code = 1216
    varobj.valid_range = [-18000, 18000] # physical limit, not sensor limit
    
    varobj = cdf.createVariable('Roll','f4',('rec'),fill_value=floatfill)
    varobj.units = "hundredths of degrees"
    varobj.long_name = "INST Roll"
    varobj.epic_code = 1217
    varobj.valid_range = [-18000, 18000] # physical limit, not sensor limit

    varobj = cdf.createVariable('HdgSTD','f4',('rec'),fill_value=floatfill)
    varobj.units = "degrees"
    varobj.long_name = "Heading Standard Deviation"

    varobj = cdf.createVariable('PtchSTD','f4',('rec'),fill_value=floatfill)
    varobj.units = "tenths of degrees"
    varobj.long_name = "Pitch Standard Deviation"

    varobj = cdf.createVariable('RollSTD','f4',('rec'),fill_value=floatfill)
    varobj.units = "tenths of degrees"
    varobj.long_name = "Roll Standard Deviation"

    varobj = cdf.createVariable('Tx','f4',('rec'),fill_value=floatfill)
    varobj.units = "hundredths of degrees"
    varobj.long_name = "ADCP Transducer Temperature"
    varobj.epic_code = 3017
    varobj.valid_range = [-500, 4000]    

    varobj = cdf.createVariable('S','f4',('rec'),fill_value=floatfill)
    varobj.units = "PPT"
    varobj.long_name = "SALINITY (PPT)"
    varobj.epic_code = 40
    varobj.valid_range = [0, 40]    

    varobj = cdf.createVariable('xmitc','f4',('rec'),fill_value=floatfill)
    varobj.units = "amps"
    varobj.long_name = "transmit current"

    varobj = cdf.createVariable('xmitv','f4',('rec'),fill_value=floatfill)
    varobj.units = "volts"
    varobj.long_name = "transmit voltage"

    varobj = cdf.createVariable('dac','i2',('rec'),fill_value=intfill)
    varobj.units = "counts"
    varobj.long_name = "DAC output"

    varobj = cdf.createVariable('VDD3','i2',('rec'),fill_value=intfill)
    varobj.units = "volts"
    varobj.long_name = "battery voltage 3"

    varobj = cdf.createVariable('VDD1','i2',('rec'),fill_value=intfill)
    varobj.units = "volts"
    varobj.long_name = "battery voltage 1"

    varobj = cdf.createVariable('VDC','i2',('rec'),fill_value=intfill)
    varobj.units = "volts"
    varobj.long_name = "VDC"

    for i in range(4):
        varname = "EWD%d" % (i+1)
        varobj = cdf.createVariable(varname,'u2',('rec'),fill_value=intfill)
        varobj.units = "binary flag"
        varobj.long_name = "Error Status Word %d" % (i+1)

    varobj = cdf.createVariable('Pressure','f4',('rec'),fill_value=floatfill)
    varobj.units = "deca-pascals"
    varobj.long_name = "ADCP Transducer Pressure"
    varobj.epic_code = 4
    varobj.valid_range = [0, maxfloat]

    varobj = cdf.createVariable('PressVar','f4',('rec'),fill_value=floatfill)
    varobj.units = "deca-pascals"
    varobj.long_name = "ADCP Transducer Pressure Variance"
    varobj.valid_range = [0, 2**31]
    
    if 'VPingSetup' in ensData:
        writeDict2atts(cdf, ensData['VPingSetup'], "TRDI_VBeam_")

    if 'VBeamLeader' in ensData:
        writeDict2atts(cdf, ensData['VBeamLeader'], "TRDI_VBeam_")

    if ('VBeamVData' in ensData):
        if ensData['VBeamLeader']['Vertical_Depth_Cells'] == ensData['FLeader']['Number_of_Cells']:
            varobj = cdf.createVariable("vel5",'f4',('rec','depth'),fill_value=floatfill)
            varobj.units = "mm s-1"
            varobj.long_name = "Beam 5 velocity (mm s-1)"
            varobj.valid_range = [-32767, 32767]
            varobj = cdf.createVariable("cor5",'u2',('rec','depth'),fill_value=intfill)
            varobj.units = "counts"
            varobj.long_name = "Beam 5 correlation"
            varobj.valid_range = [0, 255]
            varobj = cdf.createVariable("AGC5",'u2',('rec','depth'),fill_value=intfill)
            varobj.units = "counts"
            varobj.long_name = "Echo Intensity (AGC) Beam 5"
            varobj.valid_range = [0, 255]
            if ('VBeamGData' in ensData):
                varobj = cdf.createVariable("PGd5",'u2',('rec','depth'),fill_value=intfill)
                varobj.units = "counts"
                varobj.long_name = "Percent Good Beam 5"
                varobj.valid_range = [0, 100]
            else:
                cdf.TRDI_VBeam_note1 = 'Vertical beam data found without Percent Good'
        else:
            print('Vertical beam data found with different number of cells.')
            cdf.TRDI_VBeam_note = 'Vertical beam data found with different number of cells. Vertical beam data not exported to netCDF'
            print('Vertical beam data not exported to netCDF')

    if ('WaveParams' in ensData):
        # no units given for any of these in the TRDI docs
        varobj = cdf.createVariable("Hs",'f4',('rec'),fill_value=floatfill)
        varobj.units = "m"
        varobj.long_name = "Significant Wave Height (m)"
        varobj = cdf.createVariable("Tp",'f4',('rec'),fill_value=floatfill)
        varobj.units = "s"
        varobj.long_name = "Peak Wave Period (s)"
        varobj = cdf.createVariable("Dp",'f4',('rec'),fill_value=floatfill)
        varobj.units = "Deg."
        varobj.long_name = "Peak Wave Direction (Deg.)"
        varobj.valid_range = [0, 360]
        varobj = cdf.createVariable("Dm",'f4',('rec'),fill_value=floatfill)
        varobj.units = "Deg."
        varobj.long_name = "Mea Peak Wave Direction (Deg.)"
        varobj.valid_range = [0, 360]
        varobj = cdf.createVariable("SHmax",'f4',('rec'),fill_value=floatfill)
        varobj.units = "m"
        varobj.long_name = "Maximum Wave Height (m)"
        varobj.note = "from zero crossing anaylsis of surface track time series"
        varobj = cdf.createVariable("SH13",'f4',('rec'),fill_value=floatfill)
        varobj.units = "m"
        varobj.long_name = "Significant Wave Height of the largest 1/3 of the waves (m)"
        varobj.note = "in the field from zero crossing anaylsis of surface track time series"
        varobj = cdf.createVariable("SH10",'f4',('rec'),fill_value=floatfill)
        varobj.units = "m"
        varobj.long_name = "Significant Wave Height of the largest 1/10 of the waves (m)"
        varobj.note = "in the field from zero crossing anaylsis of surface track time series"
        varobj = cdf.createVariable("STmax",'f4',('rec'),fill_value=floatfill)
        varobj.units = "s"
        varobj.long_name = "Maximum Peak Wave Period (s)"
        varobj.note = "from zero crossing anaylsis of surface track time series"
        varobj = cdf.createVariable("ST13",'f4',('rec'),fill_value=floatfill)
        varobj.units = "s"
        varobj.long_name = "Period associated with the peak wave height of the largest 1/3 of the waves (s)"
        varobj.note = "in the field from zero crossing anaylsis of surface track time series"
        varobj = cdf.createVariable("ST10",'f4',('rec'),fill_value=floatfill)
        varobj.units = "s"
        varobj.long_name = "Period associated with the peak wave height of the largest 1/10 of the waves (s)"
        varobj.note = "in the field from zero crossing anaylsis of surface track time series"
        varobj = cdf.createVariable("T01",'f4',('rec'),fill_value=floatfill)
        varobj.units = " " 
        varobj = cdf.createVariable("Tz",'f4',('rec'),fill_value=floatfill)
        varobj.units = " "
        varobj = cdf.createVariable("Tinv1",'f4',('rec'),fill_value=floatfill)
        varobj.units = " "
        varobj = cdf.createVariable("S0",'f4',('rec'),fill_value=floatfill)
        varobj.units = " "
        varobj = cdf.createVariable("Source",'f4',('rec'),fill_value=floatfill)
        varobj.units = " "

    if ('WaveSeaSwell' in ensData):
        # no units given for any of these in the TRDI docs
        varobj = cdf.createVariable("HsSea",'f4',('rec'),fill_value=floatfill)
        varobj.units = "m"
        varobj.long_name = "Significant Wave Height (m)"
        varobj.note = "in the sea region of the power spectrum"
        varobj = cdf.createVariable("HsSwell",'f4',('rec'),fill_value=floatfill)
        varobj.units = "m"
        varobj.long_name = "Significant Wave Height (m)"
        varobj.note = "in the swell region of the power spectrum"
        varobj = cdf.createVariable("TpSea",'f4',('rec'),fill_value=floatfill)
        varobj.units = "s"
        varobj.long_name = "Peak Wave Period (s)"
        varobj.note = "in the sea region of the power spectrum"
        varobj = cdf.createVariable("TpSea",'f4',('rec'),fill_value=floatfill)
        varobj.units = "s"
        varobj.long_name = "Peak Wave Period (s)"
        varobj.note = "in the swell region of the power spectrum"
        varobj = cdf.createVariable("DpSea",'f4',('rec'),fill_value=floatfill)
        varobj.units = "Deg."
        varobj.long_name = "Peak Wave Direction (Deg.)"
        varobj.note = "in the sea region of the power spectrum"
        varobj = cdf.createVariable("DpSwell",'f4',('rec'),fill_value=floatfill)
        varobj.units = "Deg."
        varobj.long_name = "Peak Wave Direction (Deg.)"
        varobj.note = "in the swell region of the power spectrum"
        varobj = cdf.createVariable("SeaSwellPeriod",'f4',('rec'),fill_value=floatfill)
        varobj.units = "s"
        varobj.long_name = "Transition Period between Sea and Swell (s)"     
        
    # TODO add bottom track

    return cdf

def donothing():
    # this function is for placeholders where pylint is whining about
    # missing indented blocks
    i = 0
    return i

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
def readTRDIHeader(infile):
    HeaderData = {}
    HeaderData['headerID'] = infile.read(1)
    HeaderData['sourceID'] = infile.read(1)
    HeaderData['nbytesperens'] = struct.unpack('<H',infile.read(2))[0]
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
    FLeaderData['Beam_Angle'] = int(FLeaderData['System_Configuration_MSB'][6:8],2)
    Angles = (15,20,30,0)  
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
        
    FLeaderData['Heading_Alignment_Hundredths_of_Deg.'] = struct.unpack('<h',bstream[offset+26:offset+28])[0]
    FLeaderData['Heading_Bias_Hundredths_of_Deg.'] = struct.unpack('<h',bstream[offset+28:offset+30])[0]
    
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
        VLeaderData['Second'], VLeaderData['Hundredths']*1000)
    jddt = ajd(VLeaderData['dtobj'])
    VLeaderData['julian_day_from_as_datetime_object'] = jddt
    VLeaderData['julian_day_from_julian'] = jd
    #VLeaderData['time'] = int(math.floor(jd))
    VLeaderData['time'] = jd
    #VLeaderData['time2'] = int((jd - math.floor(jd))*(24*3600*1000))
        
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
    data['Sea_Swell_Period'] = struct.unpack('<H',bstream[offset+44:offset+46])[0]

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
    # conver hours, minutes and seconds to decimal hours
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
    
def analyzepd0file(pd0File):
    # determine the input file size
    # read some ensembles, make an estimate of the number of ensembles within
    infile = open(pd0File, 'rb')
    
    while (infile.read(1) != b'\x7f') & (infile.read(1) != b'\x7f'):
        idx = infile.tell()
        if idx < 3000:
            print('Desired TRDI 7f7f ID not found within 3 kB from beginning of the file')
            infile.close()
            sys.exit(1)
        
    infile.seek(0)
    
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
    # there is aproblem though.  WHile TRDI's documentation says the V Series 
    # Sytem Configuration data is always sent, this is not the case, so reading
    # only the first ensemble will not give the ensemble size typical over the
    # entire file
    # rewind and read this several ensembles because further in the ensemble
    # length can change on files output from Velocity
    infile.seek(0)
    
    nens2check = 5
    nbytesperens = [0 for i in range(nens2check)]
    ndatatypes = [0 for i in range(nens2check)]
    
    for i in range(nens2check):
        fileposn = infile.tell()
        Header = readTRDIHeader(infile)
        ensLen = Header['nbytesperens']+2
        infile.seek(fileposn)
        ensData, ensError = parseTRDIensemble(infile.read(ensLen), 0)
        if ensError != 'None':
            print('error - problem reading the first ensemble')
            infile.close()
            sys.exit(1)
            
        print('ensemble %d has %d bytes and %d datatypes' % (ensData['VLeader']['Ensemble_Number'], 
            ensData['Header']['nbytesperens'],ensData['Header']['ndatatypes']))
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
    maxens = nbytesinfile/ensLen
    print('estimating %g ensembles in file using a %d ensemble size' % (maxens, ensLen))
        
    infile.close()
    
    print(ensData['Header'])
    print('ensemble length = %g' % ensLen)
    print('estimating %g ensembles in file' % maxens)
    
    return maxens, ensLen, ensData

    
def __main():
# TODO add - and -- types of command line arguments
    print('%s running on python %s' % (sys.argv[0], sys.version))
	
    if len(sys.argv) < 2:
        print("%s useage:" % sys.argv[0])
        print("TRDIpd0tonetcdf infilename outfilename (startingensemble endingensemble)" )
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
        print('No starting or ending ensembles specfied, processing entire file')
        goodens = [0,np.inf]
    
    print('Start file conversion at ',dt.datetime.now())
    dopd0file(infileName, outfileName, goodens)
    
    print('Finished file conversion at ',dt.datetime.now())

    
if __name__ == "__main__":
    __main()