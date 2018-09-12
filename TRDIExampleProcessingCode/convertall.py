# script to process raw data to netcdf using python module
# at the anaconda prompt in the data directory, with IOOS3 activated, runt his script as
# E:\data\MVCO14\101003_ADCP767\python>python do767py.py > output.txt

import sys
# this is important in order to import my package which is not on the python path
sys.path.append('c:\projects\python\ADCPy')
import TRDIstuff.TRDIpd0tonetcdf as pd0
import EPICstuff.ADCPcdf2ncEPIC as cdf2nc
import datetime as dt

datapath = 'C:\\data\\10811_V20784\\'

pd0File = datapath + 'Waves_and_Currents_CCB_20161213T204403.pd0'
cdfFile = datapath + '10811whV.cdf'
ncFile = datapath + '10811whV.nc'
# test subset
#pd0File = datapath + '10811subset00n.000'
#cdfFile = datapath + '10811whVsubset00n.cdf'
#ncFile = datapath + '10811whVsubset00n.nc'
attFile = datapath + 'glob_att1081.txt'
serialnum = '20784'
goodens = [0,-1]
timetype = 'CF'

settings = {}
settings['good_ensembles'] = [0, -1]
settings['orientation'] = 'UP' # uplooking ADCP, for downlookers, use DOWN
settings['transducer_offset_from_bottom'] = 0.91
settings['transformation'] = 'EARTH' # | BEAM | INST
settings['use_pressure_for_WATER_DEPTH'] = True

do_part_one = True
do_part_two = True

print(pd0File)

allstarttime = dt.datetime.now()

if do_part_one:
    # this command will display basics about the raw binary data
    maxens, ensLen, ensData, startOfData = pd0.analyzepd0file(pd0File, 1)
    
    # will print out two layers of nested dictionary
    print('data found in file ',pd0File)
    for key, value in sorted(ensData.items()):
        vtype = type(ensData[key])
        if vtype == dict:
            print('\n',key,':')
            for key1, value1 in sorted(ensData[key].items()):
                print(key1,':',value1)
            else:
                print(key,':',value)
        else:
            print('\nData ---> %s' % key)
    
    starttime = dt.datetime.now()
    print("start binary to raw cdf conversion at %s" % starttime)
    # when using this module this way, the start and end ensembles are required.
    # use Inf to indicate all ensembles.
    print('Converting from %s\n to %s\n' % (pd0File,cdfFile))
    pd0.dopd0file(pd0File, cdfFile, goodens, serialnum, timetype)
    endtime = dt.datetime.now()
    print("finished binary to raw cdf conversion at %s" % starttime)
    print("processing time was %s" % (endtime-starttime))

if do_part_two:
    settings['timetype'] = timetype
    
    starttime = dt.datetime.now()
    print("start raw cdf to EPIC conversion at %s" % starttime)
    # when using this module this way, the start and end ensembles are required.
    # use Inf to indicate all ensembles.
    cdf2nc.doEPIC_ADCPfile(cdfFile, ncFile, attFile, settings)
    endtime = dt.datetime.now()
    print("finished raw cdf to EPIC conversion at %s" % endtime)
    print("processing time was %s" % (endtime-starttime))

allendtime = dt.datetime.now()
print("For all the operations, processing time was %s" % (allendtime-allstarttime))