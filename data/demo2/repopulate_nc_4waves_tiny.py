# -*- coding: utf-8 -*-
"""
Created on Wed Oct  3 15:53:39 2018

@author: mmartini
"""

import glob
import datetime as dt
import adcpy.EPICstuff.repopulateEPIC as repo

# making the indices
data_path = ''

burstFileRoot = '10631whVwaves*.nc'

sample_rate = 2

dry_run = False  # allows checking of file names

# ------------------- the rest of this should be automatic, no user settings below
run_start_time = dt.datetime.now()
print('Start script run at ', run_start_time)

# ----------- execute

files = glob.glob(data_path + burstFileRoot)

# remove existing repo files
for file in files:
    if 'repo' in file:
        print('got repo '+file)
        files.remove(file)
        
if len(files) == 1:
    files = [files[0]]
    
print(files)

for burstFile in files:
    s = burstFile.split('.')
    newFile = s[0]+'repo.'+s[1]
    print('making repopulated file {}'.format(newFile))

    if not dry_run:
        repo.repopulateEPIC(burstFile, newFile, sample_rate)

print('End script run at ', dt.datetime.now())
print('elapsed time is {} min'.format((dt.datetime.now() - run_start_time).total_seconds() / 60))
