# post rotation and resample adjustments
# we update thigns like DELTA_T, add minima and maxima to our new dataset
from netCDF4 import Dataset
from netCDF4 import num2date
import datetime as dt
import numpy as np
import sys
# this is important in order to import my package which is not on the python path
sys.path.append('c:\projects\python\ADCPy')
#from ADCPcdf2ncEPIC import cf2EPICtime
#from ADCPcdf2ncEPIC import EPICtime2datetime
from TRDIpd0tonetcdf import ajd

datapath = "E:\\data\\MVCO15\\10573_Signature\\Nortekncsmall\\"

datafiles = ['10573sig000_1h.nc',\
             '10573sig001_1h.nc',\
             '10573sig002_1h.nc',\
             '10573sig003_1h.nc']

# ------------------ beyond here the user should not need to edit anything
def cftime2EPICtime(time, timeunits):
    timecount = np.array(time)
    # take a CF time variable and convert to EPIC time and time2
    # timecountis the integer count of minutes (for instance) since the time stamp
    # given in timeunits
    buf = timeunits.split()
    t0 = dt.datetime.strptime(buf[2]+' '+buf[3], '%Y-%m-%d %H:%M:%S')
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

allstarttime = dt.datetime.now()

for filenum in range(len(datafiles)):
#for filenum in range(0,1):

    pydata = datapath+datafiles[filenum]
    
    
    # re-open the dataset for numerical operations such as min and max
    # we have to make attribute changes, etc. so need to open with the netCDF package
    pyd = Dataset(pydata, mode="r+", format='NETCDF4')
    

    vkeys = pyd.variables.keys()
    dkeys = pyd.dimensions.keys()
    
    # add minimum and maximum attribures
    for key in pyd.variables.keys():
        if (key not in dkeys) & (key not in {'time','EPIC_time','EPIC_time2'}):
            mn = np.min(pyd[key][:])
            # why doesn't this detect NaN?
            if np.isnan(mn): mn = 1E35
            mx = np.max(pyd[key][:])
            if np.isnan(mx): mx = 1E35
            print('%s min = %f, max = %f' % (key, mn, mx))
            pyd[key].minimum = mn
            pyd[key].maximum = mx
           
    # Add back EPIC time
    timecount = pyd['time']
    timeunits = pyd['time'].units
    timecal = pyd['time'].calendar

    #time, time2 = cf2EPICtime(timecount,timeunits,'proleptic_gregorian')
    time, time2 = cftime2EPICtime(timecount,timeunits)
    print('first time = %7d and %8d' % (time[0],time2[0]))
    
    try:
        varobj = pyd.createVariable('EPIC_time','u4',('time'))
    except:
        print('EPIC_time exists, updating.')
        varobj = pyd['EPIC_time'] 
        
    varobj.units = "True Julian Day"
    varobj.epic_code = 624
    varobj.datum = "Time (UTC) in True Julian Days: 2440000 = 0000 h on May 23, 1968"
    varobj.NOTE = "Decimal Julian day [days] = time [days] + ( time2 [msec] / 86400000 [msec/day] )"  
    varobj[:] = time[:]
    
    try:
        varobj = pyd.createVariable('EPIC_time2','u4',('time'))
    except:
        print('EPIC_time2 exists, updating.')
        varobj = pyd['EPIC_time2']
            
    varobj.units = "msec since 0:00 GMT"
    varobj.epic_code = 624
    varobj.datum = "Time (UTC) in True Julian Days: 2440000 = 0000 h on May 23, 1968"
    varobj.NOTE = "Decimal Julian day [days] = time [days] + ( time2 [msec] / 86400000 [msec/day] )"   
    varobj[:] = time2[:]
    
    # re-compute DELTA_T in seconds
    dtime = np.diff(pyd['time'][:])
    dtm = dtime.mean().astype('float').round()
    u = timeunits.split()
    if u[0] == 'minutes':
        dtm = dtm*60
    elif u[0] == 'hours':
        dtm = dtm*60*60
    elif u[0] == 'milliseconds':
        dtm = dtm/1000
    elif u[0] == 'microseconds':
        dtm = dtm/1000000
    
    DELTA_T = '%d' % int(dtm)
        
    pyd.DELTA_T = DELTA_T
    print(DELTA_T)
    
    # check start and stop time
    pyd.start_time = '%s' % num2date(pyd['time'][0],pyd['time'].units)
    pyd.stop_time = '%s' % num2date(pyd['time'][-1],pyd['time'].units)
    print('cf start time %s' % pyd.start_time)
    print('cf stop time %s' % pyd.stop_time)
    
    # add the variable descriptions
    vardesc = ''
    for key in pyd.variables.keys():
        if key not in dkeys:
            vardesc = vardesc+':'+key
    
    vardesc = vardesc[1:]
    print(vardesc)
    pyd.VAR_DESC = vardesc
        
    pyd.close()
    
    
    
