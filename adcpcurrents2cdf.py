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

# adapted from Greg Dusek's adcpwaves.py by Marinna Martini, USGS Woods Hole

import sys, struct
#import os
#import io
#from netCDF4 import Dataset

#print(sys.path)

def dopd0file(pd0File, cdfFile):
    
    infile = open(pd0File, 'rb')
    
    while (infile.read(1) != b'\x7f') & (infile.read(1) != b'\x7f'):
        print('ID not found')
        
    infile.seek(0)
    
    # need to read the header from the file to know the ensemble size
    Header = readTRDIHeader(infile)
    
    if Header['sourceID'] != b'\x7f':
        print('error - this is not a currents file')
        infile.close()
        
    infile.seek(0)
    
    # read one ensemble as a chunk to asses what kind of data we have
    # note the header's number of bytes per ensemble does not include 
    # the two checksum bytes
    ensData = readTRDIensemble(infile.read(Header['nbytesperens']+2));
    
    
    #print(ensData)
    
    infile.close()
    #cdf.close()
    print('processing complete')

def readTRDIensemble(ensbytes):
    ensData = {}
    ensError = 'None'
    ensData['Header'] = parseTRDIHeader(ensbytes)
    
    for i in range(ensData['Header']['ndatatypes']):
        # go to each offset and parse depending on what we find
        offset = ensData['Header']['offsets'][i]
        raw, val = __parsenext2TRDIbytes(ensbytes, offset)
        if raw == b'\x00\x00':
            #print('Fixed Leader found at %g' % offset)
            ensData['FLeader'] = parseTRDIFixedLeader(ensbytes, offset)
            # we need this to decode the other data records
            ncells = int(ensData['FLeader']['Number of Cells'])
            nbeams = 4 # the 5th beam has it's own record
            #print(FLeader)
        elif raw == b'\x80\x00':
            #print('Variable Leader found at %g' % offset)
            ensData['VLeader'] = parseTRDIVariableLeader(ensbytes, offset)
            #print(VLeader)
        elif raw == b'\x00\x01':
            #print('Velocity found at %g' % offset)
            ensData['VData'] = parseTRDIVelocity(ensbytes, offset, ncells, nbeams)
        elif raw == b'\x00\x02':
            #print('Correlation')
            ensData['CData'] = parseTRDICorrelation(ensbytes, offset, ncells, nbeams)
        elif raw == b'\x00\x03':
            #print('Intensity found at %g' % offset)
            ensData['IData'] = parseTRDIIntensity(ensbytes, offset, ncells, nbeams)
        elif raw == b'\x00\x04':
            #print('PGood found at %g' % offset)
            ensData['GData'] = parseTRDIPercentGood(ensbytes, offset, ncells, nbeams)
        
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
    return (raw, struct.unpack('<H',raw)[0])

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
    FLeaderData['CPU Version'] = "%s.%s" % (bstream[offset+2],bstream[offset+4])
    
    FLeaderData['System Configuration LSB'] = bitstrLE(bstream[offset+4])
    # anyone who has a better way to convert these bits, pelase tell me!
    FLeaderData['System Frequency'] = int(FLeaderData['System Configuration LSB'][5:8],2)
    SysFreqs = (75,150,300,600,1200,2400)  
    FLeaderData['System Frequency'] = SysFreqs[FLeaderData['System Frequency']]
    if FLeaderData['System Configuration LSB'][4] == "1":
        FLeaderData['Beam Pattern'] = 'Convex'
    else:
        FLeaderData['Beam Pattern'] = 'Concave'
    FLeaderData['Sensor Configuration'] = int(FLeaderData['System Configuration LSB'][2:4],2) + 1
    if FLeaderData['System Configuration LSB'][1] == "1":
        FLeaderData['Transducer Head Is Attached'] = 'Yes'
    else:
        FLeaderData['Transducer Head Is Attached'] = 'No'
    if FLeaderData['System Configuration LSB'][0] == "1":
        FLeaderData['Orientation'] = 'Up-facing beams'
    else:
        FLeaderData['Orientation'] = 'Down-facing beams'
        
    FLeaderData['System Configuration MSB'] = bitstrLE(bstream[offset+5])
    FLeaderData['Beam Angle'] = int(FLeaderData['System Configuration MSB'][6:8],2)
    Angles = (15,20,30,0)  
    FLeaderData['Beam Angle'] = Angles[FLeaderData['Beam Angle']]
    FLeaderData['Beam Configuration'] = int(FLeaderData['System Configuration MSB'][0:4],2)
    if FLeaderData['Beam Configuration'] == 4:
        FLeaderData['Beam Configuration'] = '4-bm janus'
    elif FLeaderData['Beam Configuration'] == 5:
        FLeaderData['Beam Configuration'] = '5-bm janus cfig demod'
    elif FLeaderData['Beam Configuration'] == 15:
        FLeaderData['Beam Configuration'] = '5-bm janus cfig (2 demd)'
    else: FLeaderData['Beam Configuration'] = 'unknown'
    
    FLeaderData['Simulated Data'] = "%g" % bstream[offset+6]
    
    FLeaderData['Lag Length'] = "%g" % bstream[offset+7]
    FLeaderData['Number of Beams'] = "%g" % bstream[offset+8]
    FLeaderData['Number of Cells'] = "%g" % bstream[offset+9]
    rawBytes, FLeaderData['Pings Per Ensemble'] = __parsenext2TRDIbytes(bstream,offset+10)
    rawBytes, FLeaderData['Depth Cell Length cm'] = __parsenext2TRDIbytes(bstream,offset+12)
    rawBytes, FLeaderData['Blank after Transmit cm'] = __parsenext2TRDIbytes(bstream,offset+14)
    FLeaderData['Signal Processing Mode'] = "%g" % bstream[offset+16]
    FLeaderData['Low Corr Threshold'] = "%g" % bstream[offset+17]
    FLeaderData['No. Code Reps'] = "%g" % bstream[offset+18]
    FLeaderData['%Gd Minimum'] = "%g" % bstream[offset+19]
    FLeaderData['Error Velocity Threshold'] = __parsenext2TRDIbytes(bstream,offset+20)
    # TODO ping group time needs to be formatted better
    FLeaderData['Time Between Ping Groups'] = "%3d:%2d:%2d" % (bstream[offset+22],bstream[offset+23],bstream[offset+24])

    FLeaderData['Coord Transform LSB'] = bitstrLE(bstream[offset+25])
    FLeaderData['Coord Transform'] = int(FLeaderData['Coord Transform LSB'][3:5],2)
    Xforms = ('BEAM','INST','SHIP','EARTH')  
    FLeaderData['Coord Transform'] = Xforms[FLeaderData['Coord Transform']]
    if FLeaderData['Coord Transform LSB'][5] == '1':
        FLeaderData['Tilts Used'] = 'Yes'
    else:
        FLeaderData['Tilts Used'] = 'No'
    if FLeaderData['Coord Transform LSB'][6] == '1':
        FLeaderData['3-Beam Solution Used'] = 'Yes'
    else:
        FLeaderData['3-Beam Solution Used'] = 'No'
    if FLeaderData['Coord Transform LSB'][7] == '1':
        FLeaderData['Bin Mapping Used'] = 'Yes'
    else:
        FLeaderData['Bin Mapping Used'] = 'No'
        
    rawBytes, FLeaderData['Heading Alignment Hundredths of Deg.'] = __parsenext2TRDIbytes(bstream,offset+26)
    rawBytes, FLeaderData['Heading Bias Hundredths of Deg.'] = __parsenext2TRDIbytes(bstream,offset+28)
    
    FLeaderData['Sensor Source Byte'] = bitstrLE(bstream[offset+30])
    if FLeaderData['Sensor Source Byte'][1] == '1':
        FLeaderData['Calculate EC from ED, ES, and ET'] = 'Yes'
    else:
        FLeaderData['Calculate EC from ED, ES, and ET'] = 'No'
    if FLeaderData['Sensor Source Byte'][2] == '1':
        FLeaderData['Uses ED from depth sensor'] = 'Yes'
    else:
        FLeaderData['Uses ED from depth sensor'] = 'No'
    if FLeaderData['Sensor Source Byte'][3] == '1':
        FLeaderData['Uses EH from transducer heading sensor'] = 'Yes'
    else:
        FLeaderData['Uses EH from transducer heading sensor'] = 'No'
    if FLeaderData['Sensor Source Byte'][4] == '1':
        FLeaderData['Uses EP from transducer pitch sensor'] = 'Yes'
    else:
        FLeaderData['Uses EP from transducer pitch sensor'] = 'No'
    if FLeaderData['Sensor Source Byte'][5] == '1':
        FLeaderData['Uses ER from transducer roll sensor'] = 'Yes'
    else:
        FLeaderData['Uses ER from transducer roll sensor'] = 'No'
    if FLeaderData['Sensor Source Byte'][6] == '1':
        FLeaderData['Uses ES from conductivity sensor'] = 'Yes'
    else:
        FLeaderData['Uses ES from conductivity sensor'] = 'No'
    if FLeaderData['Sensor Source Byte'][7] == '1':
        FLeaderData['Uses ET from transducer temperature sensor'] = 'Yes'
    else:
        FLeaderData['Uses ET from transducer temperature sensor'] = 'No'
        
    FLeaderData['Sensor Avail Byte'] = bitstrLE(bstream[offset+31])
    if FLeaderData['Sensor Avail Byte'][1] == '1':
        FLeaderData['Speed of sound sensor available'] = 'Yes'
    else:
        FLeaderData['Speed of sound sensor available'] = 'No'
    if FLeaderData['Sensor Avail Byte'][2] == '1':
        FLeaderData['Depth sensor available'] = 'Yes'
    else:
        FLeaderData['Depth sensor available'] = 'No'
    if FLeaderData['Sensor Avail Byte'][3] == '1':
        FLeaderData['Heading sensor available'] = 'Yes'
    else:
        FLeaderData['Heading sensor available'] = 'No'
    if FLeaderData['Sensor Avail Byte'][4] == '1':
        FLeaderData['Pitch sensor available'] = 'Yes'
    else:
        FLeaderData['Pitch sensor available'] = 'No'
    if FLeaderData['Sensor Avail Byte'][5] == '1':
        FLeaderData['Roll sensor available'] = 'Yes'
    else:
        FLeaderData['Roll sensor available'] = 'No'
    if FLeaderData['Sensor Avail Byte'][6] == '1':
        FLeaderData['Conductivity sensor available'] = 'Yes'
    else:
        FLeaderData['Conductivity sensor available'] = 'No'
    if FLeaderData['Sensor Avail Byte'][7] == '1':
        FLeaderData['Temperature sensor available'] = 'Yes'
    else:
        FLeaderData['Temperature sensor available'] = 'No'
        
    rawBytes, FLeaderData['Bin 1 distance cm'] = __parsenext2TRDIbytes(bstream,offset+32)
    rawBytes, FLeaderData['Xmit pulse length cm'] = __parsenext2TRDIbytes(bstream,offset+34)
    FLeaderData['Ref Lyr Avg Starting cell'] = "%g" % bstream[offset+36]
    FLeaderData['Ref Lyr Avg Ending cell'] = "%g" % bstream[offset+37]
    FLeaderData['False Target Threshold'] = "%g" % bstream[offset+38]
    rawBytes, FLeaderData['Transmit lag distance cm'] = __parsenext2TRDIbytes(bstream,offset+40)
    FLeaderData['CPU Board Serial Number'] = ""
    for i in range(8):
        FLeaderData['CPU Board Serial Number'] = FLeaderData['CPU Board Serial Number'] + ("%x" % bstream[offset+42+i])
   
    rawBytes, FLeaderData['System Bandwidth'] = __parsenext2TRDIbytes(bstream,offset+50)
    FLeaderData['System Power'] = "%g" % bstream[offset+52]
    FLeaderData['Base Frequency Index'] = "%g" % bstream[offset+53]
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
    rawBytes, VLeaderData['Ensemble Number'] = __parsenext2TRDIbytes(bstream,offset+2)
    VLeaderData['Year'] = "%g" % bstream[offset+4]
    VLeaderData['Month'] = "%g" % bstream[offset+5]
    VLeaderData['Day'] = "%g" % bstream[offset+6]
    VLeaderData['Hour'] = "%g" % bstream[offset+7]
    VLeaderData['Minute'] = "%g" % bstream[offset+8]
    VLeaderData['Second'] = "%g" % bstream[offset+9]
    VLeaderData['Hundredths'] = "%g" % bstream[offset+10]
    VLeaderData['Ensemble # MSB'] = "%g" % bstream[offset+11]

    VLeaderData['BIT Result Byte 13'] = bitstrLE(bstream[offset+12])
    VLeaderData['Demod 1 error bit'] = int(VLeaderData['BIT Result Byte 13'][3])
    VLeaderData['Demod 0 error bit'] = int(VLeaderData['BIT Result Byte 13'][4])
    VLeaderData['Timing Card error bit'] = int(VLeaderData['BIT Result Byte 13'][6])

    rawBytes, VLeaderData['Speed of Sound'] = __parsenext2TRDIbytes(bstream,offset+14)
    rawBytes, VLeaderData['Depth of Transducer'] = __parsenext2TRDIbytes(bstream,offset+16)
    VLeaderData['Heading, Pitch, Roll units'] = "hundredths of a degree"    
    rawBytes, VLeaderData['Heading'] = __parsenext2TRDIbytes(bstream,offset+18)
    rawBytes, VLeaderData['Pitch'] = __parsenext2TRDIbytes(bstream,offset+20)
    rawBytes, VLeaderData['Roll'] = __parsenext2TRDIbytes(bstream,offset+22)
    rawBytes, VLeaderData['Salinity'] = __parsenext2TRDIbytes(bstream,offset+24)
    rawBytes, VLeaderData['Temperature'] = __parsenext2TRDIbytes(bstream,offset+26)
    VLeaderData['MPT minutes'] = "%g" % bstream[offset+28]
    VLeaderData['MPT seconds'] = "%g" % bstream[offset+29]
    VLeaderData['MPT hundredths'] = "%g" % bstream[offset+30]
    VLeaderData['H/Hdg Std Dev'] = "%g" % bstream[offset+31]
    VLeaderData['P/Pitch Std Dev'] = "%g" % bstream[offset+32]
    VLeaderData['R/Roll Std Dev'] = "%g" % bstream[offset+33]
    VLeaderData['Xmit Current'] = "%g" % bstream[offset+34] # ADC Channel 0
    VLeaderData['Xmit Voltage'] = "%g" % bstream[offset+35] # ADC Channel 1
    VLeaderData['Ambient Temp'] = "%g" % bstream[offset+36] #ADC Channel 2
    VLeaderData['Pressure (+)'] = "%g" % bstream[offset+37] #ADC Channel 3
    VLeaderData['Pressure (-)'] = "%g" % bstream[offset+38] #ADC Channel 4
    VLeaderData['Attitude Temp'] = "%g" % bstream[offset+39] #ADC Channel 5
    VLeaderData['Attitude'] = "%g" % bstream[offset+40] #ADC Channel 6
    VLeaderData['Contaminatino Sensor'] = "%g" % bstream[offset+41] #ADC Channel 7

    VLeaderData['Error Status Word Low 16 bits LSB'] = bitstrLE(bstream[offset+42])
    VLeaderData['Bus Error exception'] = int(VLeaderData['Error Status Word Low 16 bits LSB'][7])
    VLeaderData['Address Error exception'] = int(VLeaderData['Error Status Word Low 16 bits LSB'][6])
    VLeaderData['Illegal Instruction exception'] = int(VLeaderData['Error Status Word Low 16 bits LSB'][5])
    VLeaderData['Zero Divide exception'] = int(VLeaderData['Error Status Word Low 16 bits LSB'][4])
    VLeaderData['Emulator exception'] = int(VLeaderData['Error Status Word Low 16 bits LSB'][3])
    VLeaderData['Unassigned exception'] = int(VLeaderData['Error Status Word Low 16 bits LSB'][2])
    VLeaderData['Watchdog restart occurred'] = int(VLeaderData['Error Status Word Low 16 bits LSB'][1])
    VLeaderData['Battery Saver power'] = int(VLeaderData['Error Status Word Low 16 bits LSB'][0])

    VLeaderData['Error Status Word Low 16 bits MSB'] = bitstrLE(bstream[offset+43])
    VLeaderData['Pinging'] = int(VLeaderData['Error Status Word Low 16 bits MSB'][7])
    VLeaderData['Cold Wakeup occurred'] = int(VLeaderData['Error Status Word Low 16 bits MSB'][1])
    VLeaderData['Unknown Wakeup occurred'] = int(VLeaderData['Error Status Word Low 16 bits MSB'][0])

    VLeaderData['Error Status Word High 16 bits LSB'] = bitstrLE(bstream[offset+44])
    VLeaderData['Clock Read error occurred'] = int(VLeaderData['Error Status Word High 16 bits LSB'][7])
    VLeaderData['Unexpected alarm'] = int(VLeaderData['Error Status Word High 16 bits LSB'][6])
    VLeaderData['Clock jump forward'] = int(VLeaderData['Error Status Word High 16 bits LSB'][5])
    VLeaderData['Clock jump backward'] = int(VLeaderData['Error Status Word High 16 bits LSB'][4])

    VLeaderData['Error Status Word High 16 bits MSB'] = bitstrLE(bstream[offset+42])
    VLeaderData['Power Fail (Unrecorded)'] = int(VLeaderData['Error Status Word High 16 bits MSB'][4])
    VLeaderData['Spurious level 4 intr (DSP)'] = int(VLeaderData['Error Status Word High 16 bits MSB'][3])
    VLeaderData['Spurious level 5 intr (UART)'] = int(VLeaderData['Error Status Word High 16 bits MSB'][2])
    VLeaderData['Spurious level 6 intr (CLOCK)'] = int(VLeaderData['Error Status Word High 16 bits MSB'][1])
    VLeaderData['Level 7 interrupt occurred'] = int(VLeaderData['Error Status Word High 16 bits MSB'][0])

    # pressure of the water at the transducer head relative to one atmosphere (sea level)
    #VLeaderData['Pressure word byte 1'] = bitstrLE(bstream[offset+48])
    #VLeaderData['Pressure word byte 2'] = bitstrLE(bstream[offset+49])
    #VLeaderData['Pressure word byte 3'] = bitstrLE(bstream[offset+50])
    #VLeaderData['Pressure word byte 4'] = bitstrLE(bstream[offset+51])
    VLeaderData['Pressure, deca-pascals'] = bstream[offset+48]+(bstream[offset+49]<<8)+(bstream[offset+50]<<16)+(bstream[offset+51]<<24)
    VLeaderData['Pressure variance, deca-pascals'] = bstream[offset+52]+(bstream[offset+53]<<8)+(bstream[offset+54]<<16)+(bstream[offset+55]<<24)

    VLeaderData['RTC Century'] = "%g" % bstream[offset+57]
    VLeaderData['RTC Year'] = "%g" % bstream[offset+58]
    VLeaderData['RTC Month'] = "%g" % bstream[offset+59]
    VLeaderData['RTC Day'] = "%g" % bstream[offset+60]
    VLeaderData['RTC Hour'] = "%g" % bstream[offset+61]
    VLeaderData['RTC Minute'] = "%g" % bstream[offset+62]
    VLeaderData['RTC Second'] = "%g" % bstream[offset+63]
    VLeaderData['RTC Hundredths'] = "%g" % bstream[offset+64]
    
    return VLeaderData

def parseTRDIVelocity(bstream, offset, ncells, nbeams):
    # each velocity value is stored as a two byte, twos complement integer 
    # [-32768 to 32767] with the LSB sent first.  Units are mm/s.
    # A value of -32768 = 0x8000 is a bad velocity value
    # this will create a cell x beam matrix
    # where VelocityData[:][0] will get all the cells for beam 1
    VelocityData = [[0 for x in range(ncells)] for y in range(nbeams)] 
    if bstream[offset+1] != 1:
        print("expected velocity ID, instead found %g",bstream[offset+1])
        return -1
        
    ibyte = 2
    for ibeam in range(nbeams):
        for icell in range(ncells):
            VelocityData[ibeam][icell] = struct.unpack('<h',bstream[offset+ibyte:offset+ibyte+2])
            ibyte = ibyte+2            
            
    return VelocityData    

def parseTRDICorrelation(bstream, offset, ncells, nbeams):
    CorrelationData = [[0 for x in range(ncells)] for y in range(nbeams)] 
    if bstream[offset+1] != 2:
        print("expected correlation ID, instead found %g",bstream[offset+1])
        return -1
        
    ibyte = 2
    for ibeam in range(nbeams):
        for icell in range(ncells):
            #CorrelationData[ibeam][icell] = struct.unpack('<h',bstream[offset+ibyte:offset+ibyte+2])
            CorrelationData[ibeam][icell] = bstream[offset+ibyte]
            ibyte = ibyte+1            
            
    return CorrelationData    
    
def parseTRDIIntensity(bstream, offset, ncells, nbeams):
    IntensityData = [[0 for x in range(ncells)] for y in range(nbeams)] 
    if bstream[offset+1] != 3:
        print("expected intensity ID, instead found %g",bstream[offset+1])
        return -1
        
    ibyte = 2
    for ibeam in range(nbeams):
        for icell in range(ncells):
            #CorrelationData[ibeam][icell] = struct.unpack('<h',bstream[offset+ibyte:offset+ibyte+2])
            IntensityData[ibeam][icell] = bstream[offset+ibyte]
            ibyte = ibyte+1            
            
    return IntensityData    

def parseTRDIPercentGood(bstream, offset, ncells, nbeams):
    PercentGoodData = [[0 for x in range(ncells)] for y in range(nbeams)] 
    if bstream[offset+1] != 4:
        print("expected intensity ID, instead found %g",bstream[offset+1])
        return -1
        
    ibyte = 2
    for ibeam in range(nbeams):
        for icell in range(ncells):
            #CorrelationData[ibeam][icell] = struct.unpack('<h',bstream[offset+ibyte:offset+ibyte+2])
            PercentGoodData[ibeam][icell] = bstream[offset+ibyte]
            ibyte = ibyte+1            
            
    return PercentGoodData    
    
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
    
def __main():
    
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
        
    dopd0file(infileName, outfileName)
    
if __name__ == "__main__":
    __main()