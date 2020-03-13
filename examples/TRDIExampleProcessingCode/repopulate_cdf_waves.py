# -*- coding: utf-8 -*-
"""
Created on Wed Oct  3 15:53:39 2018

@author: mmartini
"""

import glob
import datetime as dt
import os
import adcpy.EPICstuff.repopulateEPIC as repo

# making the indices
input_path = r'E:\data\Matanzas\V23881\candidatefiles'
output_path = r'E:\data\Matanzas\V23881\pythonls'

burstFileRoot = r'\11101whVwaves*.cdf'

sample_rate = 2

dry_run = False  # allows checking of file names

# ------------------- the rest of this should be automatic, no user settings below
run_start_time = dt.datetime.now()
print('Start script run at ', run_start_time)

# ----------- execute

files = glob.glob(input_path + burstFileRoot)
print(output_path + burstFileRoot)
print(files)

# remove existing repo files
for path_file in files:
    file = r'\\' + os.path.split(path_file)[1]

    if 'repo' in file:
        print('got repo '+output_path +file)
        files.remove(output_path + file)
        
if len(files) == 1:
    files = [files[0]]
    
print(files)

for path_file in files:
    burstFile = r'\\' + os.path.split(path_file)[1]

    s = burstFile.split('.')
    newFile = s[0]+'repo.'+s[1]
    print('making repopulated file {}'.format(newFile))

    if not dry_run:
        repo.repopulateEPIC(input_path + '\\' + burstFile, output_path + newFile, sample_rate)

print('End script run at ', dt.datetime.now())
print('elapsed time is {} min'.format((dt.datetime.now() - run_start_time).total_seconds() / 60))