# -*- coding: utf-8 -*-
"""
Created on Wed Oct  3 15:53:39 2018

@author: mmartini
"""

import sys
import glob
import datetime as dt
# this is important in order to import my package which is not on the python path
sys.path.append('c:\projects\python\ADCPy\EPICstuff')
import repopulateEPIC

# making the indeces
datapath = 'E:\\data\\Sandwich\\10811_V20784\\python\\'

#burstFileRoot = '10811whVwave*.cdf'
burstFileRoot = '10811whVwave00.cdf'

sample_rate = 2

dryrun = False

# ------------------- the rest of this should be automatic, no user settings below
opstart = dt.datetime.now()
print('Start script run at ',opstart)

# ----------- execute

files = glob.glob(datapath+burstFileRoot)

# remove existing repo files
for file in files:
    if 'repo' in file:
        #print('got repo '+file)
        files.remove(file)
        
if len(files) == 1:
    files = [files[0]]
    
print(files)

for burstFile in files:
    s = burstFile.split('.')
    newFile = s[0]+'repo.'+s[1]
    print('making repopulated file {}'.format(newFile))

    if not dryrun:
        repopulateEPIC.repopulateEPIC(burstFile, newFile, sample_rate)

print('End script run at ',dt.datetime.now())
print('elapsed time is {} min'.format((dt.datetime.now()-opstart).total_seconds()/60))