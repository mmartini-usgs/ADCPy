# -*- coding: utf-8 -*-
"""
This code takes a raw netcdf file containing data from any 4 beam Janus 
acoustic doppler profiler, with or without a center beam, and transforms the
data into Earth coordinates.  Data are output to netCDF using controlled
vocabulary for the variable names, following the EPIC convention wherever
possible.  This means given along beam velocity for each 
beam are used to compute East, North, and two redundant vertical velocities.
The difference between the two vertical velocities can be considered as the
error velocity, as originally introduced for the RD Instruments workhorse ADCP.

Output is stored in a netcdf file structured according to PMEL EPIC conventions.

Depth dependent attributes are compute from the mean Pressure found in the raw
data file.  So it is best to have the time series trimmed to the in water
time or to provide the good ensemble indeces for in water time

Created on Tue May 16 13:33:31 2017

@author: mmartini
"""
import sys, struct, math, os
import numpy as np 
from netCDF4 import Dataset
import datetime as dt
from datetime import datetime, timedelta
import xarray as xr

def doEPIC_ADCPfile(cdfFile, ncFile, attFile, goodens):
    
    # this function will operate on the files using the netCDF package
    nc = setupEPICnc(ncFile, cdfFile, attFile, goodens)
    
    # at the end of writing data, add these attributes
    # DELTA_T
    # VAR_DESC
    
    return
    

def setupEPICnc(fname, rawcdf, attfile, goodens):
     
    # note that 
    # f4 = 4 byte, 32 bit float
    maxfloat = 3.402823*10**38;
    intfill = -32768
    floatfill = 1E35
    
    rawcdf = Dataset(cdfFile, mode='r',format='NETCDF4')
    
    # check the ensemble limits asked for by the user
    nens = rawcdf.variables['rec'].size
    if goodens[0] < 0:
        goodens[0] = 0

    if goodens[1] > nens:
        goodens[1] = nens-1
                 
    print('creating netCDF file %s with %d records' % (fname, nens))
    
    nbins = len(rawcdf.dimensions['depth'])
    
    cdf = Dataset(fname, "w", clobber=True, format="NETCDF4")
    
    # dimensions, in EPIC order
    cdf.createDimension('time',nens)
    cdf.createDimension('depth',nbins)
    cdf.createDimension('lat',1)
    cdf.createDimension('lon',1)
    
    # write global attributes
    cdf.history = rawcdf.history+"rotations calculated and converted to EPIC format by ADCPcdf2ncEPIC.py"
    
    # these get returned as a dictionary
    gatts = read_globalatts(attfile)
    
    writeDict2atts(cdf, gatts, "")
    
    # some time things we will need and need to check
    t = gregorian(rawcdf.variables['time'][0])
    print("data starts at %s" % t['dtobj'])
    cdf.start_time = "%s" % t['dtobj']
    lastens = len(rawcdf.variables['time'])
    t = gregorian(rawcdf.variables['time'][lastens-1])
    print("data ends at %s" % t['dtobj'])
    cdf.stop_time = "%s" % t['dtobj']
    
    # more standard attributes
    cdf.latitude_units = "degree_north"
    cdf.longitude_units = "degree_east"
    cdf.initial_isntrument_height_note = "height in meters above bottom: accurate for tripod mounted instruments" 
    cdf.CREATION_DATE = "%s" % datetime.datetime.now()
    cdf.nominal_sensor_depth_note = "inst_depth = (water_depth - inst_height); nominal depth below surface, meters"
    cdf.DATA_TYPE = "ADCP"
    cdf.FILL_FLAG = 0
    cdf.COMPOSITE = 0
    
    # attributes that the names will vary depending on the ADCP vendor
    # TODO make a function template that feeds these to a dictionary depending
    # on whose ADCP this is
    # TRDI attributes
    cdf.bin_size = rawcdf.TRDI_Depth_Cell_Length_cm/100
    cdf.bin_count = rawcdf.TRDI_Number_of_Cells
    cdf.center_first_bin = rawcdf.TRDI_Bin_1_distance_cm/100
    cdf.blanking_distance = rawcdf.TRDI_Blank_after_Transmit_cm/100
    cdf.transform = rawcdf.TRDI_Coord_Transform
    cdf.beam_angle = rawcdf.TRDI_Beam_Angle
    # cdf.magnetic_variation_applied
    # cdf.magnetic_variation_applied_note = "as stated by user, not provided by Veloity processing"
    
    # attributes requiring user input
    # TODO create a way to get this information from the user
    # hardwired for now
    cdf.transducer_offset_from_bottom = 0
    cdf.orientation = "UP" #rawcdf.TRDI_Orientation # need translation to UP from "Up-facing beams"
    cdf.orientation_note = "hardwired by programmer"
    cdf.depth_note = "uplooking bin depths = WATER_DEPTH-transducer_offset_from_bottom-bindist"
    cdf.initial_instrument_height = 0
    cdf.serial_number = "TBD"
    cdf.magnetic_variation_at_site = 0
    #cdf.WATER_DEPTH = 
    #cdf.WATER_DEPTH_source = 
    #cdf.WATER_DEPTH_datum = "not yet assigned" 
    
    # computed
    # have to make things numpy to do what would be really simple in MATLAB
    #p  = np.array(cdf.variables['Pressure'])
    # we are going to assume here that the time series has been trimmed to a
    # valid in water period of time
    cdf.nominal_instrument_depth = np.mean(cdf.variables["Pressure"])
    cdf.nominal_sensor_depth = np.mean(cdf.variables["Pressure"])
    
    varobj = cdf.createVariable('Rec','u4',('time'),fill_value=intfill)
    varobj.units = "count"
    varobj.long_name = "Ensemble Number"
    # the ensemble number is a two byte LSB and a one byte MSB (for the rollover)
    varobj.valid_range = [0, 2**23]

    # if f8, 64 bit is not used, time is clipped
    # for ADCP fast sampled, single ping data, need millisecond resolution
    varobj = cdf.createVariable('time','f8',('time'))
    #varobj.units = "milliseconds since 1968-5-23 00:00:00.0 0:00" # UTC is understood
    #varobj.units = "days since 1968-5-23 00:00:00.0 0:00"
    #varobj.NOTE = "Decimal Julian day [days] = time [days] + ( time2 [msec] / 86400000 [msec/day] )"    
    # for cf convention, always assume UTC for now, and use the UNIX Epoch as the reference
    # note that ncBrowse can read netCDF files with this time convention
    varobj.units = "seconds since %d-%d-%d %d:%d:%f 0:00" % (ensData['VLeader']['Year'],
        ensData['VLeader']['Month'],ensData['VLeader']['Day'],ensData['VLeader']['Hour'],
        ensData['VLeader']['Minute'],ensData['VLeader']['Second']+
        ensData['VLeader']['Hundredths']/100)
    
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

    # TODO need to figure out when this depth information gets put in
    varobj = cdf.createVariable('depth','f4',('depth'),fill_value=floatfill)
    varobj.units = "m"
    varobj.long_name = "DEPTH (M)"
    varobj.epic_code = 3
    varobj.center_first_bin = cdf.center_first_bin
    varobj.blanking_distance = cdf.blanking_distance
    varobj.bin_size = cdf.bin_size
    varobj.bin_count = nbins
    #varobj.WATER_DEPTH = cdf.WATER_DEPTH
    #varobj.WATER_DEPTH_source = cdf.WATER_DEPTH_source
    #varobj.WATER_DEPTH_datum = "not yet assigned" 
    
    varobj = cdf.createVariable('lat','f4',('lat'),fill_value=floatfill)
    varobj.units = "degree_north"
    varobj.epic_code = 500
    varobj.name = "LAT"
    varobj.long_name = "LATITUDE"
    varobj.generic_name = "lat"
    varobj.datum = "NAD83"

    varobj = cdf.createVariable('lon','f4',('lon'),fill_value=floatfill)
    varobj.units = "degree_east"
    varobj.epic_code = 502
    varobj.name = "LON"
    varobj.long_name = "LONGITUDE"
    varobj.generic_name = "lon"
    varobj.datum = "NAD83"

    varobj = cdf.createVariable('sv','f4',('time'),fill_value=floatfill)
    varobj.units = "m s-1"
    varobj.long_name = "sound velocity (m s-1)"
    varobj.valid_range = [1400, 1600]
    
    rawcdf.close()
    cdf.close()


"""

    
    for i in range(4):
        varname = "vel%d" % (i+1)
        varobj = cdf.createVariable(varname,'f4',('time','depth'),fill_value=floatfill)
        varobj.units = "mm s-1"
        varobj.long_name = "Beam %d velocity (mm s-1)" % (i+1)
        varobj.epic_code = 1280+i
        varobj.valid_range = [-32767, 32767]
    
    for i in range(4):
        varname = "cor%d" % (i+1)
        varobj = cdf.createVariable(varname,'u2',('time','depth'),fill_value=intfill)
        varobj.units = "counts"
        varobj.long_name = "Beam %d correlation" % (i+1)
        varobj.epic_code = 1294+i
        varobj.valid_range = [0, 255]

    for i in range(4):
        varname = "AGC%d" % (i+1)
        varobj = cdf.createVariable(varname,'u2',('time','depth'),fill_value=intfill)
        varobj.units = "counts"
        varobj.epic_code = 1221+i
        varobj.long_name = "Echo Intensity (AGC) Beam %d" % (i+1)
        varobj.valid_range = [0, 255]

    if ('GData' in ensData):
        for i in range(4):
            varname = "PGd%d" % (i+1)
            varobj = cdf.createVariable(varname,'u2',('time','depth'),fill_value=intfill)
            varobj.units = "counts"
            varobj.long_name = "Percent Good Beam %d" % (i+1)
            varobj.epic_code = 1241+i
            varobj.valid_range = [0, 100]

    varobj = cdf.createVariable('Hdg','f4',('time'),fill_value=floatfill)
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
    
    varobj = cdf.createVariable('Ptch','f4',('time'),fill_value=floatfill)
    varobj.units = "hundredths of degrees"
    varobj.long_name = "INST Pitch"
    varobj.epic_code = 1216
    varobj.valid_range = [-18000, 18000] # physical limit, not sensor limit
    
    varobj = cdf.createVariable('Roll','f4',('time'),fill_value=floatfill)
    varobj.units = "hundredths of degrees"
    varobj.long_name = "INST Roll"
    varobj.epic_code = 1217
    varobj.valid_range = [-18000, 18000] # physical limit, not sensor limit

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
    varobj.valid_range = [-500, 4000]    

    varobj = cdf.createVariable('S','f4',('time'),fill_value=floatfill)
    varobj.units = "PPT"
    varobj.long_name = "SALINITY (PPT)"
    varobj.epic_code = 40
    varobj.valid_range = [0, 40]    

    varobj = cdf.createVariable('xmitc','f4',('time'),fill_value=floatfill)
    varobj.units = "amps"
    varobj.long_name = "transmit current"

    varobj = cdf.createVariable('xmitv','f4',('time'),fill_value=floatfill)
    varobj.units = "volts"
    varobj.long_name = "transmit voltage"

    varobj = cdf.createVariable('dac','i2',('time'),fill_value=intfill)
    varobj.units = "counts"
    varobj.long_name = "DAC output"

    varobj = cdf.createVariable('VDD3','i2',('time'),fill_value=intfill)
    varobj.units = "volts"
    varobj.long_name = "battery voltage 3"

    varobj = cdf.createVariable('VDD1','i2',('time'),fill_value=intfill)
    varobj.units = "volts"
    varobj.long_name = "battery voltage 1"

    varobj = cdf.createVariable('VDC','i2',('time'),fill_value=intfill)
    varobj.units = "volts"
    varobj.long_name = "VDC"

    for i in range(4):
        varname = "EWD%d" % (i+1)
        varobj = cdf.createVariable(varname,'u2',('time'),fill_value=intfill)
        varobj.units = "binary flag"
        varobj.long_name = "Error Status Word %d" % (i+1)

    varobj = cdf.createVariable('Pressure','f4',('time'),fill_value=floatfill)
    varobj.units = "deca-pascals"
    varobj.long_name = "ADCP Transducer Pressure"
    varobj.epic_code = 4
    varobj.valid_range = [0, maxfloat]

    varobj = cdf.createVariable('PressVar','f4',('time'),fill_value=floatfill)
    varobj.units = "deca-pascals"
    varobj.long_name = "ADCP Transducer Pressure Variance"
    varobj.valid_range = [0, 2**31]
    
    # TODO test this BT stuff with bottom track datas
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
            # TODO need to find a 64 bit unsigned integer - or double - designation here
            varobj = cdf.createVariable(varname,'f4',('time'),fill_value=floatfill)
            varobj.units = "cm"
            varobj.long_name = "BT Range %d" % (i+1)
            varobj.valid_range = [0, 65536*16777215]

        for i in range(4):
            varnames = ('BTWe','BTWu','BTWv','BTWd')
            longnames = ('BT Error Velocity','BT Eastward Velocity','BT Northward Velocity','BT Vertical Velocity')
            if ensData['FLeader']['Coord_Transform'] == 'EARTH':
                # TODO need to find u2 signed integer, 16 or 32 bit equivalent for this variable's declaration
                varobj = cdf.createVariable(varnames[i+1],'u2',('time'),fill_value=intfill)
                varobj.units = "mm s-1"
                varobj.long_name = "%s, mm s-1" % longnames[i+1]
                varobj.valid_range = [-32768, 32767]
                
            else:
                # TODO need to find u2 signed integer, 16 or 32 bit equivalent for this variable's declaration
                varname = "BTV%d" % (i+1)
                varobj = cdf.createVariable(varname,'u2',('time'),fill_value=intfill)
                varobj.units = "mm s-1"
                varobj.long_name = "BT velocity, mm s-1 %d" % (i+1)
                varobj.valid_range = [-32768, 32767]
                
        for i in range(4):
            varname = "BTc%d" % (i+1)
            varobj = cdf.createVariable(varname,'u2',('time'),fill_value=intfill)
            varobj.units = "counts"
            varobj.long_name = "BT correlation %d" % (i+1)
            varobj.valid_range = [0, 255]
                
        for i in range(4):
            varname = "BTe%d" % (i+1)
            varobj = cdf.createVariable(varname,'u2',('time'),fill_value=intfill)
            varobj.units = "counts"
            varobj.long_name = "BT evaluation amplitude %d" % (i+1)
            varobj.valid_range = [0, 255]
            
        for i in range(4):
            varname = "BTp%d" % (i+1)
            varobj = cdf.createVariable(varname,'u2',('time'),fill_value=intfill)
            varobj.units = "percent"
            varobj.long_name = "BT percent good %d" % (i+1)
            varobj.valid_range = [0, 100]

        for i in range(4):
            varname = "BTRSSI%d" % (i+1)
            varobj = cdf.createVariable(varname,'u2',('time'),fill_value=intfill)
            varobj.units = "counts"
            varobj.long_name = "BT Receiver Signal Strength Indicator %d" % (i+1)
            varobj.valid_range = [0, 255]

        if ensData['BTData']['Mode'] == 0: # water reference layer was used
            varobj = cdf.createVariable('BTRmin','f4',('time'),fill_value=floatfill)
            varobj.units = 'dm'
            varobj.long_name = "BT Ref. min"
            varobj.valid_range = [0,999]
            varobj = cdf.createVariable('BTRnear','f4',('time'),fill_value=floatfill)
            varobj.units = 'dm'
            varobj.long_name = "BT Ref. near"
            varobj.valid_range = [0,9999]
            varobj = cdf.createVariable('BTRfar','f4',('time'),fill_value=floatfill)
            varobj.units = 'dm'
            varobj.long_name = "BT Ref. far"
            varobj.valid_range = [0,9999]
                
            for i in range(4):
                # TODO need to find u2 signed integer, 16 or 32 bit equivalent for this variable's declaration
                varname = "BTRv%d" % (i+1)
                varobj = cdf.createVariable(varname,'f4',('time'),fill_value=floatfill)
                varobj.units = "mm s-1"
                varobj.long_name = "BT Ref. velocity, mm s-1 %d" % (i+1)
                varobj.valid_range = [-32768, 32767]

            for i in range(4):
                varname = "BTRc%d" % (i+1)
                varobj = cdf.createVariable(varname,'u2',('time'),fill_value=intfill)
                varobj.units = "counts"
                varobj.long_name = "BT Ref. correlation %d" % (i+1)
                varobj.valid_range = [0, 255]
                    
            for i in range(4):
                varname = "BTRi%d" % (i+1)
                varobj = cdf.createVariable(varname,'u2',('time'),fill_value=intfill)
                varobj.units = "counts"
                varobj.long_name = "BT Ref. intensity %d" % (i+1)
                varobj.valid_range = [0, 255]
                
            for i in range(4):
                varname = "BTRp%d" % (i+1)
                varobj = cdf.createVariable(varname,'u2',('time'),fill_value=intfill)
                varobj.units = "percent"
                varobj.long_name = "BT Ref. percent good %d" % (i+1)
                varobj.epic_code = 1269+i
                varobj.valid_range = [0, 100]
        
        
    
    if 'VPingSetup' in ensData:
        writeDict2atts(cdf, ensData['VPingSetup'], "TRDI_VBeam_")

    if 'VBeamLeader' in ensData:
        writeDict2atts(cdf, ensData['VBeamLeader'], "TRDI_VBeam_")

    if ('VBeamVData' in ensData):
        if ensData['VBeamLeader']['Vertical_Depth_Cells'] == ensData['FLeader']['Number_of_Cells']:
            varobj = cdf.createVariable("vel5",'f4',('time','depth'),fill_value=floatfill)
            varobj.units = "mm s-1"
            varobj.long_name = "Beam 5 velocity (mm s-1)"
            varobj.valid_range = [-32767, 32767]
            varobj = cdf.createVariable("cor5",'u2',('time','depth'),fill_value=intfill)
            varobj.units = "counts"
            varobj.long_name = "Beam 5 correlation"
            varobj.valid_range = [0, 255]
            varobj = cdf.createVariable("AGC5",'u2',('time','depth'),fill_value=intfill)
            varobj.units = "counts"
            varobj.long_name = "Echo Intensity (AGC) Beam 5"
            varobj.valid_range = [0, 255]
            if ('VBeamGData' in ensData):
                varobj = cdf.createVariable("PGd5",'u2',('time','depth'),fill_value=intfill)
                varobj.units = "counts"
                varobj.long_name = "Percent Good Beam 5"
                varobj.valid_range = [0, 100]
            else:
                cdf.TRDI_VBeam_note1 = 'Vertical beam data found without Percent Good'
        else:
            print('Vertical beam data found with different number of cells.')
            cdf.TRDI_VBeam_note = 'Vertical beam data found with different number of cells. Vertical beam data not exported to netCDF'
            print('Vertical beam data not exported to netCDF')
    return cdf
"""


# this code taken in part from RPS' gregorian.m, hundredths of sec added
def gregorian(jd):
    
    # jd is a julian decimal day
    print("jd = %f" % jd)
    secs = (jd % 1)*24*3600
    msecs = (jd % 1)*24*3600*1000
    j = math.floor(jd)-1721119
    i = (4*j)-1
    y = math.floor(i/146097)
    j = i - 146097*y
    i = math.floor(j/4)
    i = 4*i +3
    j = math.floor(i/1461)
    d = math.floor(((i - 1461*j) +4)/4)
    i = 5*d -3
    m = math.floor(i/153)
    d = math.floor(((i - 153*m) +5)/5)
    y = y*100 +j
    mo=m-9
    yr=y+1
    if m<10:
        mo = m+3
        yr = y
    
    hr = msecs/(3600*1000)
    hour = math.floor(msecs/(3600*1000))
    mn = (hr - hour)*60
    minu = math.floor(mn)
    sc = (mn - minu)*60
    sec = math.floor(sc)
    hund = (sc-sec)*100
    
    g = {}
    g['ymdhmsh'] = [yr, mo, d, hour, minu, sec, hund]
    # why multiplying hundredths by 10000? because datetime wants microseconds
    g['dtobj'] = dt.datetime(yr, mo, d, hour, minu, sec, math.floor(hund*10000))
    return g
    
def read_globalatts(fname):        
    # read_globalatts: read in file of metadata for a tripod or mooring
    # usage : gatt=read_globalatts(fname)
    #
    # reads global attributes for an experiment from a text file (fname)
    # called by all data processing programs to get uniform metadata input
    #  one argument is required- the name of the file to read- it should have
    #  this form:
    # SciPi;C. Sherwood
    # PROJECT; ONR
    # EXPERIMENT; RIPPLES DRI
    # DESCRIPTION; Stress, SSC, and Bedforms at MVCO 12-m fine/coarse transition site
    # DATA_SUBTYPE; MOORED
    # DATA_ORIGIN; USGS WHFS Sed Trans Group
    # COORD_SYSTEM; GEOGRAPHIC
    # Conventions; PMEL/EPIC
    # MOORING; 836
    # WATER_DEPTH; 10.99
    # latitude; 41.336063
    # longitude; -70.559615
    # magnetic_variation; -15
    # Deployment_date; 27-Aug-2007
    # Recovery_date;  ?
    gatts = {}
    f = open(fname,'r')
    for line in f:
        line = line.strip()
        cols = line.split(";")
        gatts[cols[0]] = cols[1].strip()
        
    f.close()
    return gatts
    
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

    
def __main():
    # TODO add - and -- types of command line arguments
    print('%s running on python %s' % (sys.argv[0], sys.version))
	
    if len(sys.argv) < 2:
        print("%s useage:" % sys.argv[0])
        print("ADCPcdf2ncEPIC rawcdfname ncEPICname USGSattfile [startingensemble endingensemble]\n" )
        print("starting and ending ensemble are netcdf file indeces, NOT TRDI ensemble numbers")
        print("USGSattfile is a file containing EPIC metadata")
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
        
    try:
        attfileName = sys.argv[3]
    except:
        print('error - global attribute file name missing')
        sys.exit(1)
        
    # TODO - fix these so they are looking int he correct path
    """    
    # some input testing
    if ~os.path.isfile(infileName):
        print('error - input file not found')
        sys.exit(1)
        
    if ~os.path.isfile(attfileName):
        print('error - attribute file not found')
        sys.exit(1)
    """

    print('Converting %s to %s' % (infileName, outfileName))

    try:
        goodens = [int(sys.argv[4]), int(sys.argv[5])]
    except:
        print('No starting and ending ensembles specfied, processing entire file')
        goodens = [0,np.inf]
    
    print('Start file conversion at ',dt.datetime.now())
    doEPIC_ADCPfile(infileName, outfileName, attfileName, goodens)
    
    print('Finished file conversion at ',dt.datetime.now())


if __name__ == "__main__":
    __main()