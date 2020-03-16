# -*- coding: utf-8 -*-
"""
Created on Wed Oct  3 15:53:39 2018

@author: mmartini
"""

import sys
import glob
import datetime as dt
import adcpy.EPICstuff.repopulateEPIC as repop

# making the indeces
datapath = 'E:\\data\\Sandwich\\10811_V20784\\python\\'

burstFileRoot = '10811whVcurrents00.nc'
#burstFileRoot = '10811whVcurrentssmallb00.nc'

sample_rate = 2

dryrun = False

# ------------------- the rest of this should be automatic, no user settings below
opstart = dt.datetime.now()
print('Start script run at ',opstart)

# ----------- execute

files = glob.glob(datapath+burstFileRoot)
if len(files) == 1:
    files = [files[0]]

for burstFile in files:
    s = burstFile.split('.')
    newFile = s[0]+'repo.'+s[1]
    print('making repopulated file {}'.format(newFile))

    if not dryrun:
        repop(burstFile, newFile, sample_rate)

print('End script run at ',dt.datetime.now())
print('elapsed time is {} min'.format((dt.datetime.now()-opstart).total_seconds()/60))