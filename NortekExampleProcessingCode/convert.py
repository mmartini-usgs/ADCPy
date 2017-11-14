# script to process raw data to netcdf using python module
# at the anaconda prompt in the data directory, with IOOS3 activated, run this script as
# E:\data\MVCO14\101003_ADCP767\python>python do767py.py > output.txt

# ---- this script handles output from MIDAS

import sys
# this is important in order to import my package which is not on the python path
sys.path.append('c:\projects\python\ADCPy')
import Norteknc2USGScdf as ntk
import ADCPcdf2ncEPIC as cdf2nc
#import os
import datetime as dt

datapathin = "E:\\data\\Matanzas\\WellTest2017\\Signature\\S100593A010_WHOI_well\\"
datapathout = "E:\\data\\Matanzas\\WellTest2017\\Signature\\python2\\"

# this is output from the MIDAS program
datafileroot = 'S100593A010_WHOI_well.ad2cp.00000_'
nfiles = 24

outfileroot = '1108sig'
attFile = 'E:\\data\\Matanzas\\WellTest2017\\glob_att1108.txt'
good_ens = [0,-1]

settings = {}
# add things to the settings dictionary to add them to global attributes
settings['good_ensembles'] = [0,-1]
#settings['good_ensembles'] = [0, 2048*10]
settings['orientation'] = 'DOWN' # uplooking ADCP, for downlookers, use DOWN
settings['transducer_offset_from_bottom'] = 2.02
# turn this off if you do not have data that is trimmed to in water only
settings['use_pressure_for_WATER_DEPTH'] = False 
settings['adjust_to_UTC'] = 0 # for EST to UTC, if no adjustment, set to 0 or omit

do_part_one = True # Nortek nc to USGS cdf
do_part_two = True # apply rotations, output EPIC file
    
# --------------  beyond here the user should not need to change things
timetype = 'CF'
allstarttime = dt.datetime.now()
    
#for filenum in range(1,nfiles+1):
for filenum in range(1,2):

    print('\n--------------\n')
    NortekncFile = datapathin + ('%s%d.nc' % (datafileroot,filenum))
    print(NortekncFile)
    rawcdfFile = datapathout + ('%s%03d.cdf' % (outfileroot,filenum))
    print(rawcdfFile)
    EPICFile = datapathout + ('%s%03d.nc' % (outfileroot,filenum))
    print(EPICFile)
    print('\n')
    
    if do_part_one:
        starttime = dt.datetime.now()
        print("start binary to raw cdf conversion at %s" % starttime)
        # when using this module this way, the start and end ensembles are required.
        # use Inf to indicate all ensembles.
        print('Converting from %s\n and %s\n to %s\n' % (NortekncFile,'',rawcdfFile))
        ntk.doNortekRawFile(NortekncFile, '', rawcdfFile, good_ens, timetype)

        # try:
        #     ntk.doNortekRawFile(NortekncBFile, NortekncIFile, rawcdfFile, good_ens, timetype)
        # except:
        #     sys.exit(1)
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


