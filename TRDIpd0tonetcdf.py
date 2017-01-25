"""
This code takes the raw currents portion of the adcp data from the splitter
program pd0.py and outputs raw current profile data to a netCDF4 file.

As a script:

python adcpcurrents2cdf.py [path] pd0File cdfFile

where:
    path         is a path to prepend to the following
    pd0File      is path of raw PD0 format input file with current ensembles
    cdfFile      is path of a netcdf4 EPIC compliant output file

Programmed according to the TRDI Workhorse Commands and Output Data Format document, March 2005
"""

import sys, struct
import numpy as np 
from netCDF4 import Dataset
import datetime as dt
import math

def dopd0file(pd0File, cdfFile, goodens):
    
    ensCount = 0
    verbose = 0 # diagnostic, 1 = turn on output, 0 = silent
    
    infile = open(pd0File, 'rb')
    
    while (infile.read(1) != b'\x7f') & (infile.read(1) != b'\x7f'):
        print('ID not found')
        
    infile.seek(0)
    
    # need to read the header from the file to know the ensemble size
    Header = readTRDIHeader(infile)
    
    if Header['sourceID'] != b'\x7f':
        print('error - this is not a currents file')
        infile.close()
        
    # number of bytes per ensemble in the header does not include the checksum
    ensLen = Header['nbytesperens']+2
    print('ensemble length = %g' % ensLen)
    # it is faster to define the netCDF file with a known length
    # for this we need to estimate how many ensembles we will be reading
    # for some reason, sys.getsizeof(infile) does not report the true length 
    # of the input file, so we will go to the end and see how far we have gone
    infile.seek(0,2)
    nbytesinfile = infile.tell()
    maxens = nbytesinfile/ensLen
    print('estimating %g ensembles in file' % maxens)
    
    if goodens[1] == np.inf:
        goodens[1] = maxens
        
    # rewind and read this first ensemble because we need things 
    # to set up the output file    
    infile.seek(0)
    ensData, ensError = parseTRDIensemble(infile.read(ensLen), verbose)
    if ensError == 'None':
        print('Setting up netCDF file %s' % cdfFile)
        # we are good to go, get the output file ready
        cdf = setupCdf(cdfFile, ensData, goodens)
        cdfIdx = 0
        #ncells = ensData['FLeader']['Number_of_Cells']
        nslantbeams = 4
    else:
        print('error - problem reading the first ensemble')
        sys.exit(1)
        
    # rewind to start to do the full file
    infile.seek(0)
    # priming read - for the while loop
    ens = infile.read(ensLen)
    
    verbose = 0

    while len(ens) > 0:
        #print('ensemble %d length %g, file position %g' % (ensCount, len(ens), infile.tell()))
        ensData, ensError = parseTRDIensemble(ens, verbose)
        
        if (ensError == 'None') & (ensCount >= goodens[0]):
            # write to netCDF
            if cdfIdx == 0:
                print('first ensembles read at %s and TRDI #%d' % (              
                    ensData['VLeader']['timestr'], ensData['VLeader']['Ensemble_Number']))
                
            # TODO add VLeaderData['Ensemble_#_MSB'] 
            varobj = cdf.variables['Rec']
            try:
                varobj[cdfIdx] = ensData['VLeader']['Ensemble_Number']
            except:
                print(ensData['VLeader']['Ensemble_Number'])
                print(cdfIdx)
                print((goodens[1]-goodens[0]-1))
                cdf.close()
                sys.exit(1)  

            # time calculations done when vleader is read
            varobj = cdf.variables['time']
            varobj[cdfIdx] = ensData['VLeader']['time']
            varobj = cdf.variables['time2']
            varobj[cdfIdx] = ensData['VLeader']['time2']
            
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

            cdfIdx += 1
            
        ensCount += 1
        
        if ensCount > maxens:
            print('stopping at estimated end of file ensemble %d' % goodens[1])
            break        
            
        ensf, ensi = math.modf(ensCount/100)
        if ensf == 0:
            #print('%d ensembles read at %s and TRDI #%d' % (ensCount,              
            #    ensData['VLeader']['timestr'], ensData['VLeader']['Ensemble_Number']))
            print('%d ensembles read at %s and TRDI #%d' % (ensCount,              
                ensData['VLeader']['dtobj'], ensData['VLeader']['Ensemble_Number']))
        
        if ensCount >= goodens[1]-1:
            print('stopping at requested ensemble %d' % goodens[1])
            break
        
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
    
def setupCdf(fname, ensData, gens):
    
    nens = gens[1]-gens[0]-1
    print('creating netCDF file %s with %d records' % (fname, nens))
    
    cdf = Dataset(fname, "w", clobber=True, format="NETCDF4")
    
    # dimensions, in EPIC order
    cdf.createDimension('time',nens)
    cdf.createDimension('depth',int(ensData['FLeader']['Number_of_Cells']))
    cdf.createDimension('lat',1)
    cdf.createDimension('lon',1)
    
    # write global attributes
    cdf.history = "translated to netCDF by adcpcurrents2cdf.py"
    
    writeDict2atts(cdf, ensData['FLeader'])
    
    varobj = cdf.createVariable('time','int',('time'))
    varobj.units = "True Julian Day"
    varobj.epic_code = 624
    varobj.datum = "Time (UTC) in True Julian Days: 2440000 = 0000 h on May 23, 1968"
    varobj.NOTE = "Decimal Julian day [days] = time [days] + ( time2 [msec] / 86400000 [msec/day] )"    
    
    varobj = cdf.createVariable('time2','int',('time'))
    varobj.units = "msec since 0:00 GMT"
    varobj.epic_code = 624
    varobj.datum = "Time (UTC) in True Julian Days: 2440000 = 0000 h on May 23, 1968"
    varobj.NOTE = "Decimal Julian day [days] = time [days] + ( time2 [msec] / 86400000 [msec/day] )"    

    varobj = cdf.createVariable('Rec','int',('time'),fill_value=-32768)
    varobj.units = "count"
    varobj.long_name = "Ensemble Number"
    varobj.valid_range = [0, 2**15]

    varobj = cdf.createVariable('sv','float',('time'),fill_value=1E35)
    varobj.units = "m s-1"
    varobj.long_name = "sound velocity (m s-1)"
    varobj.valid_range = [1400, 1600]
    
    for i in range(4):
        varname = "vel%d" % (i+1)
        varobj = cdf.createVariable(varname,'float',('time','depth'),fill_value=-32768)
        varobj.units = "mm s-1"
        varobj.long_name = "Beam %d velocity (mm s-1)" % (i+1)
        varobj.epic_code = 1280+i
        varobj.valid_range = [-32767, 32767]
    
    for i in range(4):
        varname = "cor%d" % (i+1)
        varobj = cdf.createVariable(varname,'float',('time','depth'),fill_value=-32768)
        varobj.units = "counts"
        varobj.long_name = "Beam %d correlation" % (i+1)
        varobj.epic_code = 1294+i
        varobj.valid_range = [0, 255]

    for i in range(4):
        varname = "AGC%d" % (i+1)
        varobj = cdf.createVariable(varname,'float',('time','depth'),fill_value=-32768)
        varobj.units = "counts"
        varobj.epic_code = 1221+i
        varobj.long_name = "Echo Intensity (AGC) Beam %d" % (i+1)
        varobj.valid_range = [0, 255]

    for i in range(4):
        varname = "PGd%d" % (i+1)
        varobj = cdf.createVariable(varname,'float',('time','depth'),fill_value=-32768)
        varobj.units = "counts"
        varobj.long_name = "Percent Good Beam %d" % (i+1)
        varobj.epic_code = 1241+i
        varobj.valid_range = [0, 100]

    varobj = cdf.createVariable('Hdg','float',('time'),fill_value=1E35)
    varobj.units = "hundredths of degrees"
    varobj.long_name = "INST Heading"
    varobj.epic_code = 1215
    varobj.heading_alignment = ensData['FLeader']['Heading_Alignment_Hundredths_of_Deg.']
    varobj.heading_bias = ensData['FLeader']['Heading_Bias_Hundredths_of_Deg.']
    varobj.valid_range = [0, 360]
    if ensData['FLeader']['Heading_Bias_Hundredths_of_Deg.'] == 0:
        varobj.NOTE_9 = "no heading bias was applied by EB during deployment or by wavesmon"
    else:
        varobj.NOTE_9 = "a heading bias was applied by EB during deployment or by wavesmon"
    
    varobj = cdf.createVariable('Ptch','float',('time'),fill_value=1E35)
    varobj.units = "hundredths of degrees"
    varobj.long_name = "INST Pitch"
    varobj.epic_code = 1216
    varobj.valid_range = [-180, 180] # physical limit, not sensor limit
    
    varobj = cdf.createVariable('Roll','float',('time'),fill_value=1E35)
    varobj.units = "hundredths of degrees"
    varobj.long_name = "INST Roll"
    varobj.epic_code = 1217
    varobj.valid_range = [-180, 180] # physical limit, not sensor limit

    varobj = cdf.createVariable('HdgSTD','float',('time'),fill_value=1E35)
    varobj.units = "degrees"
    varobj.long_name = "Heading Standard Deviation"

    varobj = cdf.createVariable('PtchSTD','float',('time'),fill_value=1E35)
    varobj.units = "tenths of degrees"
    varobj.long_name = "Pitch Standard Deviation"

    varobj = cdf.createVariable('RollSTD','float',('time'),fill_value=1E35)
    varobj.units = "tenths of degrees"
    varobj.long_name = "Roll Standard Deviation"

    varobj = cdf.createVariable('Tx','float',('time'),fill_value=1E35)
    varobj.units = "hundredths of degrees"
    varobj.long_name = "ADCP Transducer Temperature"
    varobj.epic_code = 3017
    varobj.valid_range = [-5, 40]    

    varobj = cdf.createVariable('S','float',('time'),fill_value=1E35)
    varobj.units = "PPT"
    varobj.long_name = "SALINITY (PPT)"
    varobj.epic_code = 40
    varobj.valid_range = [0, 40]    

    varobj = cdf.createVariable('xmitc','float',('time'),fill_value=1E35)
    varobj.units = "amps"
    varobj.long_name = "transmit current"

    varobj = cdf.createVariable('xmitv','float',('time'),fill_value=1E35)
    varobj.units = "volts"
    varobj.long_name = "transmit voltage"

    varobj = cdf.createVariable('dac','float',('time'),fill_value=1E35)
    varobj.units = "counts"
    varobj.long_name = "DAC output"

    varobj = cdf.createVariable('VDD3','float',('time'),fill_value=1E35)
    varobj.units = "volts"
    varobj.long_name = "battery voltage 3"

    varobj = cdf.createVariable('VDD1','float',('time'),fill_value=1E35)
    varobj.units = "volts"
    varobj.long_name = "battery voltage 1"

    varobj = cdf.createVariable('VDC','float',('time'),fill_value=1E35)
    varobj.units = "volts"
    varobj.long_name = "VDC"

    for i in range(4):
        varname = "EWD%d" % (i+1)
        varobj = cdf.createVariable(varname,'float',('time'),fill_value=1E35)
        varobj.units = "binary flag"
        varobj.long_name = "Error Status Word %d" % (i+1)

    varobj = cdf.createVariable('Pressure','float',('time'),fill_value=1E35)
    varobj.units = "deca-pascals"
    varobj.long_name = "ADCP Transducer Pressure"
    varobj.epic_code = 4
    varobj.valid_range = [0, 2**31]

    varobj = cdf.createVariable('PressVar','float',('time'),fill_value=1E35)
    varobj.units = "deca-pascals"
    varobj.long_name = "ADCP Transducer Pressure Variance"
    varobj.valid_range = [0, 2**31]

    return cdf
    
def matrixTranspose( matrix ):
    if not matrix: return []
    return [ [ row[ i ] for row in matrix ] for i in range( len( matrix[ 0 ] ) ) ]

# write a dictionary to netCDF attributes
def writeDict2atts(cdfobj, d):
    
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
        newkey = "TRDI_" + key
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
        raw, val = __parsenext2TRDIbytes(ensbytes, offset)
        if raw == b'\x00\x00':
            if verbose: print('Fixed Leader found at %g' % offset)
            ensData['FLeader'] = parseTRDIFixedLeader(ensbytes, offset)
            # we need this to decode the other data records
            ncells = int(ensData['FLeader']['Number_of_Cells'])
            nbeams = 4 # the 5th beam has it's own record
            #print(FLeader)
        elif raw == b'\x80\x00':
            if verbose: print('Variable Leader found at %g' % offset)
            ensData['VLeader'] = parseTRDIVariableLeader(ensbytes, offset)
            #print(VLeader)
        elif raw == b'\x00\x01':
            if verbose: print('Velocity found at %g' % offset)
            ensData['VData'] = parseTRDIVelocity(ensbytes, offset, ncells, nbeams)
        elif raw == b'\x00\x02':
            if verbose: print('Correlation found at %g' % offset)
            ensData['CData'] = parseTRDICorrelation(ensbytes, offset, ncells, nbeams)
        elif raw == b'\x00\x03':
            if verbose: print('Intensity found at %g' % offset)
            ensData['IData'] = parseTRDIIntensity(ensbytes, offset, ncells, nbeams)
        elif raw == b'\x00\x04':
            if verbose: print('PGood found at %g' % offset)
            ensData['GData'] = parseTRDIPercentGood(ensbytes, offset, ncells, nbeams)
        elif raw == b'\x00\x32':
            print('Instrument transformation found at %g' % offset)
        elif raw == b'\x00\x70':
            print('V Series sytem config found at %g' % offset)
        elif raw == b'\x01\x70':
            print('V Series ping setup found at %g' % offset)
        elif raw == b'\x02\x70':
            print('V Series ADC Data found at %g' % offset)
        elif raw == b'\x03\x70':
            print('V Series Features Header Data found at %g' % offset)
        elif raw == b'\x01\x0f':
            print('Vertical Beam Leader Data found at %g' % offset)
        elif raw == b'\x00\x0a':
            print('Vertical Beam Velocity Data found at %g' % offset)
        elif raw == b'\x00\x0b':
            print('Vertical Beam Correlation Data found at %g' % offset)
        elif raw == b'\x00\x0c':
            print('Vertical Beam Amplitude Data found at %g' % offset)
        elif raw == b'\x00\x0d':
            print('Vertical Beam Percent Good Data found at %g' % offset)
        elif raw == b'\x40\x70':
            print('V Series Event Log Data found at %g' % offset)
        elif raw == b'\x0b\x00':
            print('Wavesmon 4 Wave Parameters found at %g' % offset)
        elif raw == b'\x0c\x00':
            print('Wavesmon 4 Sea and Swell found at %g' % offset)
        else:
            print('no known ID found')
            ensError = 'no ID'
        
    csum = __computeChecksum(ensbytes)
    if csum != (ensbytes[-2]+(ensbytes[-1]<<8)):
        ensError = 'checksum failure'
        
    return ensData, ensError
    

def bitstrLE(byte):
    bits = ""
    for i in [7,6,5,4,3,2,1,0]:
    #for i in range(7,0,-1):
        #print(i)
        if (byte >> i) & 1:
            bits+="1"
        else:
            bits+="0"
    return bits

def bitstrBE(byte): # make a bit string from big endian byte
    # surely there's a better way to do this!!
    bits = ""
    for i in range(8): # Big Endian
        #print(i)
        if (byte[0] >> i) & 1:
            bits+="1"
        else:
            bits+="0"
    return bits

#TODO - break this into a separate file for module-wide use
def __next2TRDIbytes(file):
    # read the next two bytes as little Endian, Unsigned Short
    raw = file.read(2)
    return (raw, struct.unpack('<H',raw)[0])
    
def __parsenext2TRDIbytes(bstream, i):
    # read the next two bytes as little Endian, Unsigned Short
    raw = bstream[i:i+2]
    try:
        data = struct.unpack('<H',raw)[0]
        
    except:
        print('could not unpack %g', raw)
        data = ''
        
    return (raw, data)

#TODO - break this into a separate file for module-wide use
def readTRDIHeader(infile):
    HeaderData = {}
    HeaderData['headerID'] = infile.read(1)
    HeaderData['sourceID'] = infile.read(1)
    rawBytes, HeaderData['nbytesperens'] = __next2TRDIbytes(infile)
    infile.read(1) # spare, skip it
    HeaderData['ndatatypes'] = infile.read(1)[0] # remember, bytes objects are arrays
    offsets = [0]*HeaderData['ndatatypes'] # predefine a list of ints to fill
    for i in range(HeaderData['ndatatypes']):
        rawBytes, data = __next2TRDIbytes(infile)
        offsets[i] = data
        
    HeaderData['offsets'] = offsets

    return HeaderData

#TODO - break this into a separate file for module-wide use
def parseTRDIHeader(bstream):
    HeaderData = {}
    HeaderData['headerID'] = bstream[0] # byte 1
    HeaderData['sourceID'] = bstream[1] # byte 2
    rawBytes, HeaderData['nbytesperens'] = __parsenext2TRDIbytes(bstream, 2) # byte 3-4
    # spare, skip it, byte 5
    HeaderData['ndatatypes'] = bstream[5] # byte 6
    offsets = [0]*HeaderData['ndatatypes'] # predefine a list of ints to fill
    for i in range(HeaderData['ndatatypes']):
        rawBytes, data = __parsenext2TRDIbytes(bstream, 6+i*2)
        offsets[i] = data
        
    HeaderData['offsets'] = offsets

    return HeaderData
    
#TODO - break this into a separate file for module-wide use
def parseTRDIFixedLeader(bstream, offset):
    # bstream is a bytes object that contains an entire ensemble
    # offset is the location in the bytes object of the first byte of this data format
    FLeaderData = {}
    rawBytes, leaderID = __parsenext2TRDIbytes(bstream,offset)
    if leaderID != 0:
        print("expected fixed leader ID, instead found %b",leaderID)
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
    
    FLeaderData['Simulated_Data'] = "%g" % bstream[offset+6]
    
    FLeaderData['Lag_Length'] = "%g" % bstream[offset+7]
    FLeaderData['Number_of_Beams'] = "%g" % bstream[offset+8]
    FLeaderData['Number_of_Cells'] = "%g" % bstream[offset+9]
    rawBytes, FLeaderData['Pings_Per_Ensemble'] = __parsenext2TRDIbytes(bstream,offset+10)
    rawBytes, FLeaderData['Depth_Cell_Length_cm'] = __parsenext2TRDIbytes(bstream,offset+12)
    rawBytes, FLeaderData['Blank_after_Transmit_cm'] = __parsenext2TRDIbytes(bstream,offset+14)
    FLeaderData['Signal_Processing_Mode'] = "%g" % bstream[offset+16]
    FLeaderData['Low_Corr_Threshold'] = "%g" % bstream[offset+17]
    FLeaderData['No._Code_Reps'] = "%g" % bstream[offset+18]
    FLeaderData['PGd_Minimum'] = "%g" % bstream[offset+19]
    FLeaderData['Error_Velocity_Threshold'] = __parsenext2TRDIbytes(bstream,offset+20)
    # TODO ping group time needs to be formatted better
    FLeaderData['Time_Between_Ping Groups'] = "%3d:%2d:%2d" % (bstream[offset+22],bstream[offset+23],bstream[offset+24])

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
        
    rawBytes, FLeaderData['Heading_Alignment_Hundredths_of_Deg.'] = __parsenext2TRDIbytes(bstream,offset+26)
    rawBytes, FLeaderData['Heading_Bias_Hundredths_of_Deg.'] = __parsenext2TRDIbytes(bstream,offset+28)
    
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
        
    rawBytes, FLeaderData['Bin_1_distance_cm'] = __parsenext2TRDIbytes(bstream,offset+32)
    rawBytes, FLeaderData['Xmit_pulse_length_cm'] = __parsenext2TRDIbytes(bstream,offset+34)
    FLeaderData['Ref_Lyr_Avg_Starting_cell'] = "%g" % bstream[offset+36]
    FLeaderData['Ref_Lyr_Avg_Ending_cell'] = "%g" % bstream[offset+37]
    FLeaderData['False_Target_Threshold'] = "%g" % bstream[offset+38]
    rawBytes, FLeaderData['Transmit_lag_distance_cm'] = __parsenext2TRDIbytes(bstream,offset+40)
    FLeaderData['CPU_Board_Serial_Number'] = ""
    for i in range(8):
        FLeaderData['CPU_Board_Serial_Number'] = FLeaderData['CPU_Board_Serial_Number'] + ("%x" % bstream[offset+42+i])
   
    rawBytes, FLeaderData['System_Bandwidth'] = __parsenext2TRDIbytes(bstream,offset+50)
    FLeaderData['System_Power'] = "%g" % bstream[offset+52]
    FLeaderData['Base_Frequency_Index'] = "%g" % bstream[offset+53]
    #TODO these two need to be interpreted as spare if WH ADCP
    #rawBytes, FLeaderData['Serial Number for Remus only'] = __next2TRDIbytes(infile)
    #FLeaderData['Beam Angle for H-ADCP only'] = "%g" % infile.read(1)[0]
    
    return FLeaderData

#TODO - break this into a separate file for module-wide use
def parseTRDIVariableLeader(bstream, offset):
    # bstream is a bytes object that contains an entire ensemble
    # offset is the location in the bytes object of the first byte of this data format
    VLeaderData = {}
    rawBytes, leaderID = __parsenext2TRDIbytes(bstream,offset)
    if leaderID != 128:
        print("expected variable leader ID, instead found %b",leaderID)
        return -1
    rawBytes, VLeaderData['Ensemble_Number'] = __parsenext2TRDIbytes(bstream,offset+2)
    VLeaderData['Year'] = int("%g" % bstream[offset+4])
    if VLeaderData['Year'] < 50: # circa 2000
        VLeaderData['Year'] += 2000
    else:
        VLeaderData['Year'] += 1900
        
    VLeaderData['Month'] = int("%g" % bstream[offset+5])
    VLeaderData['Day'] = int("%g" % bstream[offset+6])
    VLeaderData['Hour'] = int("%g" % bstream[offset+7])
    VLeaderData['Minute'] = int("%g" % bstream[offset+8])
    VLeaderData['Second'] = int("%g" % bstream[offset+9])
    VLeaderData['Hundredths'] = int("%g" % bstream[offset+10])
    VLeaderData['Ensemble_#_MSB'] = "%g" % bstream[offset+11]

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
    jd = ajd(VLeaderData['dtobj'])
    VLeaderData['julian day'] = jd
    VLeaderData['time'] = int(math.floor(jd))
    VLeaderData['time2'] = int((jd - math.floor(jd))*(24*3600*1000))
        
    VLeaderData['BIT_Result_Byte_13'] = bitstrLE(bstream[offset+12])
    VLeaderData['Demod_1_error_bit'] = int(VLeaderData['BIT_Result_Byte_13'][3])
    VLeaderData['Demod_0_error_bit'] = int(VLeaderData['BIT_Result_Byte_13'][4])
    VLeaderData['Timing_Card_error_bit'] = int(VLeaderData['BIT_Result_Byte_13'][6])

    rawBytes, VLeaderData['Speed_of_Sound'] = __parsenext2TRDIbytes(bstream,offset+14)
    rawBytes, VLeaderData['Depth_of_Transducer'] = __parsenext2TRDIbytes(bstream,offset+16)
    VLeaderData['Heading, Pitch, Roll units'] = "hundredths_of_a_degree"    
    rawBytes, VLeaderData['Heading'] = __parsenext2TRDIbytes(bstream,offset+18)
    rawBytes, VLeaderData['Pitch'] = __parsenext2TRDIbytes(bstream,offset+20)
    rawBytes, VLeaderData['Roll'] = __parsenext2TRDIbytes(bstream,offset+22)
    rawBytes, VLeaderData['Salinity'] = __parsenext2TRDIbytes(bstream,offset+24)
    rawBytes, VLeaderData['Temperature'] = __parsenext2TRDIbytes(bstream,offset+26)
    VLeaderData['MPT_minutes'] = "%g" % bstream[offset+28]
    VLeaderData['MPT_seconds'] = "%g" % bstream[offset+29]
    VLeaderData['MPT_hundredths'] = "%g" % bstream[offset+30]
    VLeaderData['H/Hdg_Std_Dev'] = int("%g" % bstream[offset+31])
    VLeaderData['P/Pitch_Std_Dev'] = int("%g" % bstream[offset+32])
    VLeaderData['R/Roll_Std_Dev'] = int("%g" % bstream[offset+33])
    # the V Series PDO Output is different for the ADC channels        
    # V PD0 this is ADC Channel 0 not used    
    VLeaderData['Xmit_Current'] = "%g" % bstream[offset+34] # ADC Channel 0
    # V PD0 this is ADC Channel 1 XMIT Voltage    
    VLeaderData['Xmit_Voltage'] = "%g" % bstream[offset+35] # ADC Channel 1
    # V PD0 this is ADC Channel 2 not used    
    VLeaderData['Ambient_Temp'] = "%g" % bstream[offset+36] #ADC Channel 2
    # V PD0 this is ADC Channel 3 not used    
    VLeaderData['Pressure_(+)'] = "%g" % bstream[offset+37] #ADC Channel 3
    # V PD0 this is ADC Channel 4 not used    
    VLeaderData['Pressure_(-)'] = "%g" % bstream[offset+38] #ADC Channel 4
    # V PD0 this is ADC Channel 5 not used    
    VLeaderData['Attitude_Temp'] = "%g" % bstream[offset+39] #ADC Channel 5
    # V PD0 this is ADC Channel 6 not used    
    VLeaderData['Attitude'] = "%g" % bstream[offset+40] #ADC Channel 6
    # V PD0 this is ADC Channel 7 not used    
    VLeaderData['Contamination_Sensor'] = "%g" % bstream[offset+41] #ADC Channel 7

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

    VLeaderData['RTC_Century'] = "%g" % bstream[offset+57]
    VLeaderData['RTC_Year'] = "%g" % bstream[offset+58]
    VLeaderData['RTC_Month'] = "%g" % bstream[offset+59]
    VLeaderData['RTC_Day'] = "%g" % bstream[offset+60]
    VLeaderData['RTC_Hour'] = "%g" % bstream[offset+61]
    VLeaderData['RTC_Minute'] = "%g" % bstream[offset+62]
    VLeaderData['RTC_Second'] = "%g" % bstream[offset+63]
    VLeaderData['RTC_Hundredths'] = "%g" % bstream[offset+64]
    
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
            value = struct.unpack('<h',bstream[offset+ibyte:offset+ibyte+2])
            data[ibeam,icell] = value[0]
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
    
    # factored for readability
def __computeChecksum(ensemble):
    """Compute a checksum from header, length, and ensemble"""
    cs = 0    
    for byte in range(len(ensemble)-2):
        # since the for loop returns an int to byte, use as-is in python 3x
        #value = struct.unpack('B', byte)[0]
        #cs += value
        cs += ensemble[byte]
    return cs & 0xffff
    
# these date conversion functions came from
# http://stackoverflow.com/questions/31142181/calculating-julian-date-in-python/41769526#41769526
# AND OLD rps matlab code
# I don't know what the Italy is.
# TODO - check this, or replace with better solution
def julian(year,month,day,hour,mn,sec,hund):
    # from RPS' julian.m and hms2h.m
    s = sec+hund/100
    h=hour+(mn+s/60)/60
    mo=month+9
    yr=year-1
    
    if month > 2:
        mo -= 3
        yr = year
        
    c = math.floor(yr/100)
    yr = yr - c*100
    d = day
    j = math.floor((146097*c)/4)+math.floor((1461*yr)/4)+math.floor((153*mo +2)/5)+d+1721119

    # If you want julian days to start and end at noon, 
    # replace the following line with:
    # j=j+(h-12)/24;
    j=j+h/24
    
    return j
    
def jdn(dto):
    """
    Given datetime object returns Julian Day Number
    """
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
    # Waht is this?
    #if jjs < ITALY:
    #    jjs -= reform

    return jjs
    # end jdn

def ajd(dto):
    """
    Given datetime object returns Astronomical Julian Day.
    Day is from midnight 00:00:00+00:00 with day fractional
    value added.
    """
    jdd = jdn(dto)
    day_fraction = dto.hour / 24.0 + dto.minute / 1440.0 + dto.second / 86400.0
    return jdd + day_fraction - 0.5
    # end ajd
    
def __main():
    
    print('%s running on python %s' % (sys.argv[0], sys.version))
    
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
        
    try:
        goodens = [int(sys.argv[3]), int(sys.argv[4])]
    except:
        goodens = [0,np.inf]
        
    print('Converting %s to %s' % (infileName, outfileName))

    dopd0file(infileName, outfileName, goodens)
    
if __name__ == "__main__":
    __main()