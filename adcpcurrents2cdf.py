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

# TODO change byte strings to something that appears correctly in the netCDF

import sys, struct
#import os
#import io
import netCDF4 as nc
from netCDF4 import Dataset

#print(sys.path)

def dopd0file(pd0File, cdfFile):
    
    wavesID=0x797f
    currentsID=127 #0x7f7f
    ensCount = 0
    
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
    nens = nbytesinfile/ensLen
    print('estimating %g ensembles in file' % nens)

    # rewind and read this first ensemble because we need things 
    # to set up the output file    
    infile.seek(0)
    ensData, ensError = parseTRDIensemble(infile.read(ensLen))
    if ensError == 'None':
        # we are good to go, get the output file ready
        cdf = setupCdf(cdfFile, ensData, nens)
    else:
        print('error - problem reading the first ensemble')
        sys.exit(1)
        
    # rewind to start to do the full file
    infile.seek(0)
    # priming read - for the while loop
    ens = infile.read(ensLen)

    while len(ens) > 0:
        print('ensemble %d length %g, file position %g' % (ensCount, len(ens), infile.tell()))
        ensData, ensError = parseTRDIensemble(ens)
        # write to netCDF
        p = cdf.variables['Pressure']
        p[ensCount] = ensData['VLeader']['Depth_of_Transducer']
        
        ensCount += 1
        if ensCount > nens:
            break;        
        ens = infile.read(ensLen)
    else:
        print('end of file reached')
        
    if ensCount < nens:
        print('end of file reached after %d ensembles, less than expected' % ensCount)
    elif ensCount > nens:
        print('end of file reached after %d ensembles, more than expected' % ensCount)
    
    infile.close()
    cdf.close()
    
    print('processing complete')
    
def setupCdf(fname, ensData, nens):
    cdf = Dataset(fname, "w", clobber=True, format="NETCDF4")
    
    # dimensions, in EPIC order
    cdf.createDimension('time',nens)
    cdf.createDimension('depth',int(ensData['FLeader']['Number_of_Cells']))
    cdf.createDimension('lat',1)
    cdf.createDimension('lon',1)
    
    # write global attributes
    cdf.history = "translated to netCDF by adcpcurrents2cdf.py"
    
    writeDict2atts(cdf, ensData['FLeader'])
    
    varobj = cdf.createVariable('Pressure','float',('time'),fill_value=1E35)
    varobj.units = "deca-pascals"
    varobj.long_name = "ADCP Transducer Pressure"
    varobj.valid_range = [0, 2**31]
    
    #print(len(ensData['FLeader']))
    #print(int(ensData['FLeader']['Number_of_Beams']))

    return cdf

# write a dictionary to netCDF attributes
def writeDict2atts(cdfobj, d):
    
    # first, convert as many of the values in d to numbers as we can
    for key in iter(d):
        if type(d[key]) == str:
            try:
                d[key] = float(d[key])
            except ValueError:
                # we really don't need to print here, 
                # but python insists we do something
                print('   can\'t convert %s to float' % key)

    for key in iter(d):
        newkey = "TRDI_" + key
        try:
            cdfobj.setncattr(newkey,d[key])
        except:
            print('can\'t set %s attribute' % key)
    
    # return the dictionary it's numerized values 
    return d
    

def parseTRDIensemble(ensbytes):
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
            ncells = int(ensData['FLeader']['Number_of_Cells'])
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
    VLeaderData['Year'] = "%g" % bstream[offset+4]
    VLeaderData['Month'] = "%g" % bstream[offset+5]
    VLeaderData['Day'] = "%g" % bstream[offset+6]
    VLeaderData['Hour'] = "%g" % bstream[offset+7]
    VLeaderData['Minute'] = "%g" % bstream[offset+8]
    VLeaderData['Second'] = "%g" % bstream[offset+9]
    VLeaderData['Hundredths'] = "%g" % bstream[offset+10]
    VLeaderData['Ensemble_#_MSB'] = "%g" % bstream[offset+11]

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
    VLeaderData['H/Hdg_Std_Dev'] = "%g" % bstream[offset+31]
    VLeaderData['P/Pitch_Std_Dev'] = "%g" % bstream[offset+32]
    VLeaderData['R/Roll_Std_Dev'] = "%g" % bstream[offset+33]
    VLeaderData['Xmit_Current'] = "%g" % bstream[offset+34] # ADC Channel 0
    VLeaderData['Xmit_Voltage'] = "%g" % bstream[offset+35] # ADC Channel 1
    VLeaderData['Ambient_Temp'] = "%g" % bstream[offset+36] #ADC Channel 2
    VLeaderData['Pressure_(+)'] = "%g" % bstream[offset+37] #ADC Channel 3
    VLeaderData['Pressure_(-)'] = "%g" % bstream[offset+38] #ADC Channel 4
    VLeaderData['Attitude_Temp'] = "%g" % bstream[offset+39] #ADC Channel 5
    VLeaderData['Attitude'] = "%g" % bstream[offset+40] #ADC Channel 6
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
    
    print('adcpcurrents2cdf.py running on python %s' % sys.version)
    
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
    dopd0file(infileName, outfileName)
    
if __name__ == "__main__":
    __main()