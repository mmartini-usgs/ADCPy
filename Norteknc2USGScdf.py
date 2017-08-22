"""
This code takes Nortek exported netcdf files of Signature data and
outputs raw current profile data to a netCDF4 file.
Data are taken from the "Burst" group of "Data"

As a script:

python Norteknc2USGScdf.py [path] infileName outfileName

where:
    path         is a path to prepend to the following
    infileName   is path of netCDF4 input file from a Nortek Signature
    outfileName      is path of a netcdf4 output file
    start        ensemble at which to start exporting
    end          ensemble at which to stop exporting
    
note that file names and paths may not include spaces
    
As a module:
import Norteknc2USGScdf as pd0

Notes:
    time and time2, the EPIC convention for netCDF, is not used here so that
    the resulting very large files generated can be reduced using existing 
    python too ls such as xarrays
    
    Files exported by the Contour program have differently named attributes
    and this program may be updated to handle them at a later date.

"""

# 1/25/2017 MM got this running on old Workhorse ADCP data

import sys, math
import numpy as np 
from netCDF4 import Dataset
from netCDF4 import num2date
import datetime as dt
from TRDIpd0tonetcdf import julian

def doNortekRawFile(infileName, outfileName, goodens, timetype):
    
    nc = Dataset(infileName, mode='r', format='NETCDF4')
    
    maxens = len(nc['Data']['Burst']['time'])
        
    if goodens[1] == np.inf:
        goodens[1] = maxens
           
    # we are good to go, get the output file ready
    print('Setting up netCDF output file %s' % outfileName)

    # set up some pointers to the netCDF groups
    config = nc['Config']
    data = nc['Data']['Burst']
    idata = nc['Data']['IBurst']
    
    # TODO - pay attention to the possible number of bursts.
    # we are assuming here that the first burst is the primary sample set of
    # the four slant beams
     
    # note that 
    # f4 = 4 byte, 32 bit float
    maxfloat = 3.402823*10**38;
    intfill = -32768
    floatfill = 1E35
    
    nens = goodens[1]-goodens[0]
    print('creating netCDF file %s with %d records' % (outfileName, nens))
    
    cdf = Dataset(outfileName, 'w', clobber=True, format='NETCDF4')
    
    # dimensions, in EPIC order
    cdf.createDimension('time',nens)
    cdf.createDimension('depth',config.burst_nCells)
    cdf.createDimension('lat',1)
    cdf.createDimension('lon',1)
    
    # write global attributes
    cdf.history = "translated to USGS netCDF by Norteknc2USGScdf.py"
    cdf.sensor_type = 'Nortek'
    cdf.serial_number = config.serialNumberDoppler
    
    # TODO - reduce the number if attributes we copy from the nc file
    # build a dictionary of global attributes is a faster way to load attributes
    # into a netcdf file http://unidata.github.io/netcdf4-python/#netCDF4.Dataset.setncatts

    # put the "sensor_type" in front of the attributes that come directly
    # from the instrument data        
    Nortek_config = dictifyatts(config, 'Nortek_')
    cdf.setncatts(Nortek_config)
        
    # it's not yet clear which way to go with this.  python tools like xarray 
    # and panoply demand that time be a CF defined time.
    # USGS CMG MATLAB tools need time and time2
    # create the datetime object from the CF time
    tobj = num2date(data['time'][:],data['time'].units,calendar=data['time'].calendar)
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
                   math.floor(tobj[idx].microsecond/10))
        jd.append(j)
        time.append(int(math.floor(j)))
        time2.append(int((j - math.floor(j))*(24*3600*1000)))
    if timetype=='EPIC':
        # we include cf_time for cf compliance and use by python packages like xarray
        # if f8, 64 bit is not used, time is clipped
        # for ADCP fast sampled, single ping data, need millisecond resolution
        varobj = cdf.createVariable('cf_time','f8',('time'))
        varobj.setncatts(dictifyatts(data['time'],''))
        varobj[:] = data['time'][:]
        # we include time and time2 for EPIC compliance
        varobj = cdf.createVariable('time','u4',('time'))
        varobj.units = "True Julian Day"
        varobj.epic_code = 624
        varobj.datum = "Time (UTC) in True Julian Days: 2440000 = 0000 h on May 23, 1968"
        varobj.NOTE = "Decimal Julian day [days] = time [days] + ( time2 [msec] / 86400000 [msec/day] )"    
        varobj[:] = time[:]
        varobj = cdf.createVariable('time2','u4',('time'))
        varobj.units = "msec since 0:00 GMT"
        varobj.epic_code = 624
        varobj.datum = "Time (UTC) in True Julian Days: 2440000 = 0000 h on May 23, 1968"
        varobj.NOTE = "Decimal Julian day [days] = time [days] + ( time2 [msec] / 86400000 [msec/day] )"    
        varobj[:] = time2[:]
    else:
        # cf_time for cf compliance and use by python packages like xarray
        # if f8, 64 bit is not used, time is clipped
        # for ADCP fast sampled, single ping data, need millisecond resolution
        varobj = cdf.createVariable('time','f8',('time'))
        varobj.setncatts(dictifyatts(data['time'],''))
        varobj[:] = data['time'][:]
        varobj = cdf.createVariable('EPIC_time','u4',('time'))
        varobj.units = "True Julian Day"
        varobj.epic_code = 624
        varobj.datum = "Time (UTC) in True Julian Days: 2440000 = 0000 h on May 23, 1968"
        varobj.NOTE = "Decimal Julian day [days] = time [days] + ( time2 [msec] / 86400000 [msec/day] )"   
        varobj[:] = time[:]
        varobj = cdf.createVariable('EPIC_time2','u4',('time'))
        varobj.units = "msec since 0:00 GMT"
        varobj.epic_code = 624
        varobj.datum = "Time (UTC) in True Julian Days: 2440000 = 0000 h on May 23, 1968"
        varobj.NOTE = "Decimal Julian day [days] = time [days] + ( time2 [msec] / 86400000 [msec/day] )"  
        varobj[:] = time2[:]
        
    cdf.start_time = '%s' % num2date(data['time'][0],data['time'].units)
    cdf.stop_time = '%s' % num2date(data['time'][-1],data['time'].units)
    print('times from the input file')
    print(cdf.start_time)
    print(cdf.stop_time)

    print('times from the output file')
    print('%s' % num2date(cdf['time'][0],cdf['time'].units))
    print('%s' % num2date(cdf['time'][-1],cdf['time'].units))

    varobj = cdf.createVariable('Rec','u4',('time'),fill_value=intfill)
    varobj.units = "count"
    varobj.long_name = "Ensemble Count for each burst"
    varobj.valid_range = [0, 2**23]
    varobj[:] = data['EnsembleCount'][:]

    varobj = cdf.createVariable('sv','f4',('time'),fill_value=floatfill)
    varobj.units = "m s-1"
    varobj.long_name = "sound velocity (m s-1)"
    varobj.valid_range = [1400, 1600]
    varobj[:] = data['SpeedOfSound'][:]
    
    # there are separate Amplitude_Range, Correlation_Range and Velocity_Range
    # we will pass on Velocity_Range as bindist
    varobj = cdf.createVariable('bindist','f4',('depth'),fill_value=floatfill)
    # note name is one of the netcdf4 reserved attributes, use setncattr
    varobj.setncattr('name', "bindist")
    varobj.units = "m"
    varobj.long_name = "bin distance from instrument for slant beams"
    varobj.epic_code = 0
    #varobj.valid_range = [0 0]
    varobj.NOTE = "distance is not specified by Nortek as along beam or vertical"
    varobj[:] = data['Velocity Range'][:]
    
    # map the Nortek beams onto TRDI order since later code expects TRDI order
    TRDInumber = [3,1,4,2]
    for i in range(4):
        varname = "vel%d" % TRDInumber[i]
        key = 'VelocityBeam%d' % (i+1)
        varobj = cdf.createVariable(varname,'f4',('time','depth'),fill_value=floatfill)
        varobj.units = "m s-1"
        varobj.long_name = "Beam %d velocity (m s-1)" % TRDInumber[i]
        varobj.epic_code = 1277+i
        varobj.NOTE = 'beams reordered from Nortek 1-2-3-4 to TRDI 3-1-4-2, as viewed clockwise from compass 0 degree reference, when instrument is up-looking'
        #varobj.valid_range = [-32767, 32767]
        varobj[:,:] = data[key][:,:]
    
    for i in range(4):
        varname = "cor%d" % (i+1)
        key = 'CorrelationBeam%d' % (i+1)
        varobj = cdf.createVariable(varname,'u2',('time','depth'),fill_value=intfill)
        varobj.units = "percent"
        varobj.long_name = "Beam %d correlation" % (i+1)
        #varobj.epic_code = 1285+i
        varobj.valid_range = [0, 100]
        varobj[:,:] = data[key][:,:]

    for i in range(4):
        varname = "att%d" % (i+1)
        key = 'AmplitudeBeam%d' % (i+1)
        varobj = cdf.createVariable(varname,'f4',('time','depth'),fill_value=intfill)
        varobj.units = "dB"
        #varobj.epic_code = 1281+i
        varobj.long_name = "ADCP amplitude of beam %d" % (i+1)
        #varobj.valid_range = [0, 255]
        varobj[:,:] = data[key][:,:]

    varname = 'Heading'
    varobj = cdf.createVariable('Hdg','f4',('time'),fill_value=floatfill)
    varobj.units = "degrees"
    varobj.long_name = "INST Heading"
    varobj.epic_code = 1215
    varobj.valid_range = [0, 360]
    # TODO can we tell on a Signature if a magvar was applied at deployment?
    # no metadata found in the .nc file global attributes
    #varobj.NOTE_9 = "no heading bias was applied during deployment"
    varobj[:] = data[varname][:]
    
    varname = 'Pitch'
    varobj = cdf.createVariable('Ptch','f4',('time'),fill_value=floatfill)
    varobj.units = "degrees"
    varobj.long_name = "INST Pitch"
    varobj.epic_code = 1216
    varobj.valid_range = [-18000, 18000] # physical limit, not sensor limit
    varobj[:] = data[varname][:]
    
    varname = 'Roll'
    varobj = cdf.createVariable('Roll','f4',('time'),fill_value=floatfill)
    varobj.units = "degrees"
    varobj.long_name = "INST Roll"
    varobj.epic_code = 1217
    varobj.valid_range = [-18000, 18000] # physical limit, not sensor limit
    varobj[:] = data[varname][:]

    # The Signature records magnetometer data we are not converting at this time

    varname = 'WaterTemperature'
    varobj = cdf.createVariable('Tx','f4',('time'),fill_value=floatfill)
    varobj.units = "degrees"
    varobj.long_name = "Water temperature at ADCP"
    # TODO - verify if Signature is IPTS-1990
    #  20:T  :TEMPERATURE (C)          :temp:C:f10.2:IPTS-1990 standard
    #varobj.epic_code = 28
    varobj.valid_range = [-500, 4000]    
    varobj[:] = data[varname][:]

    varname = 'Pressure'
    varobj = cdf.createVariable('Pressure','f4',('time'),fill_value=floatfill)
    varobj.units = "dbar"
    varobj.long_name = "ADCP Transducer Pressure"
    varobj.epic_code = 4
    varobj.valid_range = [0, maxfloat]
    varobj[:] = data[varname][:]

    # TODO - Signature can bottom track, and we don't have an example yet
    """
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
            varobj.valid_range = [0, 65536*16777215]

        for i in range(4):
            varnames = ('BTWe','BTWu','BTWv','BTWd')
            longnames = ('BT Error Velocity','BT Eastward Velocity','BT Northward Velocity','BT Vertical Velocity')
            if ensData['FLeader']['Coord_Transform'] == 'EARTH':
                varobj = cdf.createVariable(varnames[i+1],'i2',('time'),fill_value=intfill)
                varobj.units = "mm s-1"
                varobj.long_name = "%s, mm s-1" % longnames[i+1]
                varobj.valid_range = [-32768, 32767]
                
            else:
                varname = "BTV%d" % (i+1)
                varobj = cdf.createVariable(varname,'i2',('time'),fill_value=intfill)
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
                varname = "BTRv%d" % (i+1)
                varobj = cdf.createVariable(varname,'i2',('time'),fill_value=intfill)
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
    """   
        
    # it is possible in a Signature for the vertical beam data to be on a
    # different time base.  Test for this.  If it is the same time base we can 
    # include it now.  If it isn't we will have to add it later by some other
    # code.  5th beam Signature data is stored under the IBurst group
    # it is also possible for the number of bins to be different
    
    if (data['time'].size == idata['time'].size):
        if (nc['Data']['Burst']['Velocity Range'].size == nc['Data']['IBurst']['Velocity Range'].size):
            varobj = cdf.createVariable("vel5",'f4',('time','depth'),fill_value=floatfill)
            varobj.units = "m s-1"
            varobj.long_name = "Beam 5 velocity (m s-1)"
            #varobj.valid_range = [-32767, 32767]
            varobj[:,:] = idata['VelocityBeam5'][:,:]
            
            varobj = cdf.createVariable("cor5",'u2',('time','depth'),fill_value=intfill)
            varobj.units = "percent"
            varobj.long_name = "Beam 5 correlation"
            varobj.valid_range = [0, 100]
            varobj[:,:] = idata['CorrelationBeam5'][:,:]
            
            varobj = cdf.createVariable("att5",'u2',('time','depth'),fill_value=intfill)
            varobj.units = "dB"
            varobj.long_name = "ADCP amplitude of beam 5"
            #varobj.valid_range = [0, 255]
            varobj[:,:] = idata['AmplitudeBeam5'][:,:]

        else:
            print('Vertical beam data found with different number of cells.')
            cdf.Nortek_VBeam_note = 'Vertical beam data found with different number of cells. Vertical beam data not exported to netCDF'
            print('Vertical beam data not exported to netCDF')
    else:
        print('Vertical beam data found with different number of ensembles.')
        cdf.Nortek_VBeam_note = 'Vertical beam data found with different number of ensembles. Vertical beam data not exported to netCDF'
        print('Vertical beam data not exported to netCDF')

    nc.close()
    cdf.close()
    
    print('%d ensembles copied' % maxens)

def dictifyatts(varptr, tag):
    theDict = {}
    for key in varptr.ncattrs():
        newkey = tag+key
        theDict[newkey] = varptr.getncattr(key)

    return theDict
    
def __main():
# TODO add - and -- types of command line arguments
    print('%s running on python %s' % (sys.argv[0], sys.version))
	
    if len(sys.argv) < 2:
        print("%s useage:" % sys.argv[0])
        print("Norteknc2USGScdf infilename outfilename [startingensemble endingensemble]" )
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
        print('No starting and ending ensembles specfied, processing entire file')
        goodens = [0,np.inf]
        
    try:
        timetype = sys.argv[5]
    except:
        print('Time type will be CF')
        timetype = "CF"      
    
    print('Start file conversion at ',dt.datetime.now())
    doNortekRawFile(infileName, outfileName, goodens, timetype)
    
    print('Finished file conversion at ',dt.datetime.now())

    
if __name__ == "__main__":
    __main()