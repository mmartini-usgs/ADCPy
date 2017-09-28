# script to process raw data to netcdf using python module
# at the anaconda prompt in the data directory, with IOOS3 activated, run this script as
# E:\data\MVCO14\101003_ADCP767\python>python do767py.py > output.txt

import sys
# this is important in order to import my package which is not on the python path
sys.path.append('c:\projects\python\ADCPy')
import Norteknc2USGScdf as ntk
import ADCPcdf2ncEPIC as cdf2nc
#import os
import datetime as dt

datapath = "E:\\data\\MVCO15\\10573_Signature\\Nortekncsmall\\"

datafiles = ['S100076A006_MVCO15_00.ad2cp.00000_1.nc',\
             'S100076A006_MVCO15_00.ad2cp.00000_2.nc',\
             'S100076A006_MVCO15_01.ad2cp.00000_1.nc',\
             'S100076A006_MVCO15.ad2cp.00000.nc']
outfileroot = '10573sig'
attFile = datapath + 'glob_att1057.txt'
good_ens = [0,-1]

settings = {}
# add things to the settings dictionary to add them to global attributes
settings['good_ensembles'] = [0,-1]
#settings['good_ensembles'] = [0, 2048*10]
settings['orientation'] = 'UP' # uplooking ADCP, for downlookers, use DOWN
settings['transducer_offset_from_bottom'] = 2.02
# turn this off if you do not have data that is trimmed to in water only
settings['use_pressure_for_WATER_DEPTH'] = False 
settings['adjust_to_UTC'] = 5 # for EST to UTC, if no adjustment, set to 0 or omit

do_part_one = False # Nortek nc to USGS cdf
do_part_two = True # apply rotations, output EPIC file
    
# --------------  beyond here the user should not need to change things
for filenum in range(len(datafiles)):
#for filenum in range(1,len(datafiles)):
#for filenum in range(0,1):

    print('\n--------------\n')
    NortekncFile = datapath + datafiles[filenum]
    print(NortekncFile)
    rawcdfFile = datapath + ('%s%03d.cdf' % (outfileroot,filenum))
    print(rawcdfFile)
    EPICFile = datapath + ('%s%03d.nc' % (outfileroot,filenum))
    print(EPICFile)
    print('\n')
    
    timetype = 'CF'
    
    allstarttime = dt.datetime.now()
    
    if do_part_one:
        starttime = dt.datetime.now()
        print("start binary to raw cdf conversion at %s" % starttime)
        # when using this module this way, the start and end ensembles are required.
        # use Inf to indicate all ensembles.
        print('Converting from %s\n to %s\n' % (NortekncFile,rawcdfFile))
        try:
            ntk.doNortekRawFile(NortekncFile, rawcdfFile, good_ens, timetype)
        except:
            sys.exit(1)
        endtime = dt.datetime.now()
        print("finished binary to raw cdf conversion at %s" % starttime)
        print("processing time was %s\n" % (endtime-starttime))
    
    if do_part_two:
        settings['timetype'] = timetype
        settings['transformation'] = 'EARTH' # | BEAM | INST
        
        starttime = dt.datetime.now()
        print("start raw cdf to EPIC conversion at %s" % starttime)
        print('Converting %s\n to %s\n' % (rawcdfFile,EPICFile))
        # when using this module this way, the start and end ensembles are required.
        # use Inf to indicate all ensembles.
        cdf2nc.doEPIC_ADCPfile(rawcdfFile, EPICFile, attFile, settings)
        endtime = dt.datetime.now()
        print("finished raw cdf to EPIC conversion at %s" % endtime)
        print("processing time was %s\n" % (endtime-starttime))
    
    allendtime = dt.datetime.now()
    print("For all the operations, processing time was %s" % (allendtime-allstarttime))


