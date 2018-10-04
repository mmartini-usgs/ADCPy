"""
This code takes Nortek exported netcdf files of Signature data and
outputs raw current profile data to a netCDF4 file.
Data are taken from the "Burst" group of "Data"

As a script:

python Norteknc2USGScdf.py [path] infileBName outfileName

where:
    path         is a path to prepend to the following
    infileBName   is path of netCDF4 input Burst file from a Nortek Signature
    infileIName   is path of netCDF4 input IBurst file from a Nortek Signature
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

# TODO note that the ADCPcdf2ncEPIC keys off attributes for "burst" and
# there may be burst1, avg, avg1 etc. in a Nortek file.  This may need to
# be dealt with as separate .cdf files

# 10/4/2018 MM remove valid_range as it causes too many downstream problems

import sys, math
from netCDF4 import Dataset
from netCDF4 import num2date
import datetime as dt
from TRDIstuff.TRDIpd0tonetcdf import julian
from EPICstuff.EPICmisc import cftime2EPICtime

def doNortekRawFile(infileBName, infileIName, outfileName, goodens, timetype):
    
    if infileIName == '':
        midasdata = True
        
    # this is necessary so that this function does not change the value
    # in the calling function
    ens2process = goodens[:]
    
    nc = Dataset(infileBName, mode='r', format='NETCDF4')
    # so the MIDAS output and the Contour output are different, here we hope
    # to handle both, since MIDAS has been more tolerant of odd data
    if midasdata:
        ncI = nc
    else:
        ncI = Dataset(infileIName, mode='r', format='NETCDF4')
    
    maxens = len(nc['Data']['Burst']['time'])
    print('%s has %d ensembles' % (infileBName,maxens))
    
    # TODO - ens2process[1] has the file size from the previous file run when multiple files are processed!
    if ens2process[1] < 0:
        ens2process[1] = maxens
           
    # we are good to go, get the output file ready
    print('Setting up netCDF output file %s' % outfileName)

    # set up some pointers to the netCDF groups
    config = nc['Config']
    data = nc['Data']['Burst']
    if 'IBurstHR' in ncI['Data'].groups:
        idata = ncI['Data']['IBurstHR']
        HRdata = True
    else:
        idata = ncI['Data']['IBurst']
        HRdata = False
    
    # TODO - pay attention to the possible number of bursts.
    # we are assuming here that the first burst is the primary sample set of
    # the four slant beams
     
    # note that 
    # f4 = 4 byte, 32 bit float
    #maxfloat = 3.402823*10**38;
    intfill = -32768
    floatfill = 1E35
    
    nens = ens2process[1]-ens2process[0]
    print('creating netCDF file %s with %d records' % (outfileName, nens))
    
    cdf = Dataset(outfileName, 'w', clobber=True, format='NETCDF4')
    
    # dimensions, in EPIC order
    cdf.createDimension('time',nens)
    if midasdata: 
        cdf.createDimension('depth',config.burst_nCells)
    else:
        cdf.createDimension('depth',config.Instrument_burst_nCells)
    cdf.createDimension('lat',1)
    cdf.createDimension('lon',1)
    
    # write global attributes
    cdf.history = "translated to USGS netCDF by Norteknc2USGScdf.py"
    cdf.sensor_type = 'Nortek'
    if midasdata:
        cdf.serial_number = config.serialNumberDoppler
    else:
        cdf.serial_number = config.Instrument_serialNumberDoppler
    
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
    CFcount = data['time'][:]
    CFunits = data['time'].units
    EPICtime, EPICtime2 = cftime2EPICtime(CFcount,CFunits)
    print('CFcount[0] = %f, CFunits = %s' % (CFcount[0], CFunits))
    print('EPICtime[0] = %f, EPICtime2[0] = %f' % (EPICtime[0], EPICtime2[0]))
    elapsed_sec = []
    for idx in range(len(tobj)):
        tdelta = tobj[idx]-tobj[0] # timedelta
        elapsed_sec.append(tdelta.total_seconds())
    # from the datetime object convert to time and time2
    jd = []
    time = []
    time2 = []
    # using u2 rather than u4 here because when EPIC time is written from this
    # cdf to the nc file, it was getting messed up
    #file is 1108sig001.cdf
    #EPIC first time stamp = 25-Sep-2017 15:00:00
    #seconds since 1970-01-01T00:00:00 UTC
    #CF first time stamp = 25-Sep-2017 15:00:00
    #file is 1108sig001.nc
    #EPIC first time stamp = 08-Oct-5378 00:01:04
    #seconds since 1970-01-01T00:00:00 UTC
    #CF first time stamp = 25-Sep-2017 15:00:00
    timevartype = 'u2'
    for idx in range(len(tobj)):
        j = julian(tobj[idx].year,tobj[idx].month,tobj[idx].day, \
                   tobj[idx].hour,tobj[idx].minute,tobj[idx].second,\
                   math.floor(tobj[idx].microsecond/10))
        jd.append(j)
        time.append(int(math.floor(j)))
        time2.append(int((j - math.floor(j))*(24*3600*1000)))
    if timetype=='CF':
        # cf_time for cf compliance and use by python packages like xarray
        # if f8, 64 bit is not used, time is clipped
        # TODO test this theory, because downstream 64 bit time is a problem
        # for ADCP fast sampled, single ping data, need millisecond resolution
        #cf_time = data['time'][:]
        #cdf.createVariable('time','f8',('time'))
        #cdf['time'].setncatts(dictifyatts(data['time'],''))
        #cdf['time'][:] = cf_time[:]
        varobj = cdf.createVariable('time','f8',('time'))
        varobj.setncatts(dictifyatts(data['time'],''))
        varobj[:] = data['time'][:]
        varobj = cdf.createVariable('EPIC_time',timevartype,('time'))
        varobj.units = "True Julian Day"
        varobj.epic_code = 624
        varobj.datum = "Time (UTC) in True Julian Days: 2440000 = 0000 h on May 23, 1968"
        varobj.NOTE = "Decimal Julian day [days] = time [days] + ( time2 [msec] / 86400000 [msec/day] )"   
        #varobj[:] = time[:]
        varobj[:] = EPICtime[:]
        varobj = cdf.createVariable('EPIC_time2',timevartype,('time'))
        varobj.units = "msec since 0:00 GMT"
        varobj.epic_code = 624
        varobj.datum = "Time (UTC) in True Julian Days: 2440000 = 0000 h on May 23, 1968"
        varobj.NOTE = "Decimal Julian day [days] = time [days] + ( time2 [msec] / 86400000 [msec/day] )"  
        #varobj[:] = time2[:]
        varobj[:] = EPICtime2[:]
    else:
        # we include cf_time for cf compliance and use by python packages like xarray
        # if f8, 64 bit is not used, time is clipped
        # for ADCP fast sampled, single ping data, need millisecond resolution
        varobj = cdf.createVariable('cf_time','f8',('time'))
        varobj.setncatts(dictifyatts(data['time'],''))
        varobj[:] = data['time'][:]
        # we include time and time2 for EPIC compliance
        varobj = cdf.createVariable('time',timevartype,('time'))
        varobj.units = "True Julian Day"
        varobj.epic_code = 624
        varobj.datum = "Time (UTC) in True Julian Days: 2440000 = 0000 h on May 23, 1968"
        varobj.NOTE = "Decimal Julian day [days] = time [days] + ( time2 [msec] / 86400000 [msec/day] )"    
        #varobj[:] = time[:]
        varobj[:] = EPICtime[:]
        varobj = cdf.createVariable('time2',timevartype,('time'))
        varobj.units = "msec since 0:00 GMT"
        varobj.epic_code = 624
        varobj.datum = "Time (UTC) in True Julian Days: 2440000 = 0000 h on May 23, 1968"
        varobj.NOTE = "Decimal Julian day [days] = time [days] + ( time2 [msec] / 86400000 [msec/day] )"    
        #varobj[:] = time2[:]
        varobj[:] = EPICtime2[:]
        
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
    #varobj.valid_range = [0, 2**23]
    varobj[:] = data['EnsembleCount'][:]

    varobj = cdf.createVariable('sv','f4',('time'),fill_value=floatfill)
    varobj.units = "m s-1"
    varobj.long_name = "sound velocity (m s-1)"
    #varobj.valid_range = [1400, 1600]
    varobj[:] = data['SpeedOfSound'][:]
    
    # get the number
    
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
    if midasdata:
        vardata = data.variables['Velocity Range'][:] # in raw data
    else:        
        # because this is a coordinate variable, one can't just say data['Burst_Velocity_Beam_Range'][:]
        try: 
            vardata = data.variables['Burst Velocity_Range'][:] # in raw data
        except:
            vardata = data.variables['Burst Velocity Beam_Range'][:] # in processed data
        
    varobj[:] = vardata
    nbbins = vardata.size
    
    # map the Nortek beams onto TRDI order since later code expects TRDI order
    TRDInumber = [3,1,4,2]
    for i in range(4):
        varname = "vel%d" % TRDInumber[i]
        if midasdata:
            key = 'VelocityBeam%d' % (i+1)
        else:
            key = 'Vel_Beam%d' % (i+1)
        varobj = cdf.createVariable(varname,'f4',('time','depth'),fill_value=floatfill)
        varobj.units = "m s-1"
        varobj.long_name = "Beam %d velocity (m s-1)" % TRDInumber[i]
        varobj.epic_code = 1277+i
        varobj.NOTE = 'beams reordered from Nortek 1-2-3-4 to TRDI 3-1-4-2, as viewed clockwise from compass 0 degree reference, when instrument is up-looking'
        #varobj.valid_range = [-32767, 32767]
        varobj[:,:] = data[key][:,:]
    
    for i in range(4):
        varname = "cor%d" % (i+1)
        if midasdata:
            key = 'CorrelationBeam%d' % (i+1)
        else:
            key = 'Cor_Beam%d' % (i+1)
        varobj = cdf.createVariable(varname,'u2',('time','depth'),fill_value=intfill)
        varobj.units = "percent"
        varobj.long_name = "Beam %d correlation" % (i+1)
        #varobj.epic_code = 1285+i
        #varobj.valid_range = [0, 100]
        varobj[:,:] = data[key][:,:]

    for i in range(4):
        varname = "att%d" % (i+1)
        if midasdata:
            key = 'AmplitudeBeam%d' % (i+1)
        else:
            key = 'Amp_Beam%d' % (i+1)
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
    #varobj.valid_range = [0, 360]
    # TODO can we tell on a Signature if a magvar was applied at deployment?
    # no metadata found in the .nc file global attributes
    #varobj.NOTE_9 = "no heading bias was applied during deployment"
    varobj[:] = data[varname][:]
    
    varname = 'Pitch'
    varobj = cdf.createVariable('Ptch','f4',('time'),fill_value=floatfill)
    varobj.units = "degrees"
    varobj.long_name = "INST Pitch"
    varobj.epic_code = 1216
    #varobj.valid_range = [-18000, 18000] # physical limit, not sensor limit
    varobj[:] = data[varname][:]
    
    varname = 'Roll'
    varobj = cdf.createVariable('Roll','f4',('time'),fill_value=floatfill)
    varobj.units = "degrees"
    varobj.long_name = "INST Roll"
    varobj.epic_code = 1217
    #varobj.valid_range = [-18000, 18000] # physical limit, not sensor limit
    varobj[:] = data[varname][:]

    # The Signature records magnetometer data we are not converting at this time

    varname = 'WaterTemperature'
    varobj = cdf.createVariable('Tx','f4',('time'),fill_value=floatfill)
    varobj.units = "degrees"
    varobj.long_name = "Water temperature at ADCP"
    # TODO - verify if Signature is IPTS-1990
    #  20:T  :TEMPERATURE (C)          :temp:C:f10.2:IPTS-1990 standard
    #varobj.epic_code = 28
    #varobj.valid_range = [-500, 4000]    
    varobj[:] = data[varname][:]

    varname = 'Pressure'
    varobj = cdf.createVariable('Pressure','f4',('time'),fill_value=floatfill)
    varobj.units = "dbar"
    varobj.long_name = "ADCP Transducer Pressure"
    varobj.epic_code = 4
    #varobj.valid_range = [0, maxfloat]
    varobj[:] = data[varname][:]

    # TODO - Signature can bottom track, and we don't have an example yet
    # we will want to model it on the TRDI ADCP format that laredy exists
        
    # it is possible in a Signature for the vertical beam data to be on a
    # different time base.  Test for this.  If it is the same time base we can 
    # include it now.  If it isn't we will have to add it later by some other
    # code.  5th beam Signature data is stored under the IBurst group
    # it is also possible for the number of bins to be different
    if midasdata: 
        vrkey = 'Velocity Range'
    else:
        vrkey = 'IBurstHR Velocity_Range'
        
    if (data['time'].size == idata['time'].size):
        if (nbbins == idata.variables[vrkey].size):
            varobj = cdf.createVariable("vel5",'f4',('time','depth'),fill_value=floatfill)
            varobj.units = "m s-1"
            varobj.long_name = "Beam 5 velocity (m s-1)"
            #varobj.valid_range = [-32767, 32767]
            varobj[:,:] = idata['VelocityBeam5'][:,:]
            #varobj[:,:] = idata['Vel_Beam5'][:,:] # contour version
            
            varobj = cdf.createVariable("cor5",'u2',('time','depth'),fill_value=intfill)
            varobj.units = "percent"
            varobj.long_name = "Beam 5 correlation"
            #varobj.valid_range = [0, 100]
            varobj[:,:] = idata['CorrelationBeam5'][:,:]
            #varobj[:,:] = idata['Cor_Beam5'][:,:] # contour version
            
            varobj = cdf.createVariable("att5",'u2',('time','depth'),fill_value=intfill)
            varobj.units = "dB"
            varobj.long_name = "ADCP amplitude of beam 5"
            #varobj.valid_range = [0, 255]
            varobj[:,:] = idata['AmplitudeBeam5'][:,:]
            #varobj[:,:] = idata['Amp_Beam5'][:,:] # contour version

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
        if key.startswith('Instrument_'):
            # we need to strip the 'Instrument_' off the beginning
            n = key.find('_')
            newkey = tag+key[n+1:]
        else:
            newkey = tag+key
        theDict[newkey] = varptr.getncattr(key)

    return theDict
    
def __main():
# TODO add - and -- types of command line arguments
    print('%s running on python %s' % (sys.argv[0], sys.version))
	
    if len(sys.argv) < 3:
        print("%s useage:" % sys.argv[0])
        print("Norteknc2USGScdf infileBName infileIName outfilename [startingensemble endingensemble]" )
        sys.exit(1)
    
    try:
        infileBName = sys.argv[1]
    except:
        print('error - Burst input file name missing')
        sys.exit(1)
        
    try:
        infileIName = sys.argv[2]
    except:
        print('error - IBurst input file name missing')
        sys.exit(1)

    try:
        outfileName = sys.argv[3]
    except:
        print('error - output file name missing')
        sys.exit(1)
        
    print('Converting %s to %s' % (infileBName, outfileName))

    try:
        goodens = [int(sys.argv[4]), int(sys.argv[5])]
    except:
        print('No starting and ending ensembles specfied, processing entire file')
        goodens = [0,-1]
        
    try:
        timetype = sys.argv[6]
    except:
        print('Time type will be CF')
        timetype = "CF"      
    
    print('Start file conversion at ',dt.datetime.now())
    doNortekRawFile(infileBName, infileIName, outfileName, goodens, timetype)
    
    print('Finished file conversion at ',dt.datetime.now())

    
if __name__ == "__main__":
    __main()