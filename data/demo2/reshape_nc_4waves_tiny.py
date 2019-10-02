"""
Processing for mooring 10631 TRDI workhorse V data

@author: mmartini
"""
import netCDF4 as nc
from netCDF4 import num2date
import math
import sys
import datetime as dt
import adcpy.EPICstuff.reshapeEPIC as reshape

# making the indices
input_path = r''
output_path = r''

# continuous file
continuous_file = r'10631whV.nc'  # rotated file with a 1D, continuous time series
number_of_output_files = 1
burst_file_name = r'10631whVwaves.nc'
index_file_name = r'10631whVwavesindecesnc.txt'

sample_rate = 2
burst_length = 2048
burst_interval = 3600  # 60 min interval
burst_start_offset = 0  # at which ensemble do we start the first wave burst?
# note we are dropping EPIC time (if present) as we are moving on to CF time
variables_to_omit = {'EPIC_time', 'EPIC_time2'}
attributes_to_omit = {'valid_range'}  # this is in older converted files and needs to be removed
dry_run = False

# ------------------- the rest of this should be automatic, no user settings below
operation_start = dt.datetime.now()
print('Start script run at ', operation_start)
dim = 'time'

# ----------- execute
    
all_slices = reshape.generate_expected_start_times(input_path + continuous_file, dim,
                                                   burst_start_offset, burst_interval, burst_length, sample_rate)

# here we limit the slices for testing
# print('** reducing the number of slices')
slices = all_slices  # [200:300]

continuous_netcdf_object = nc.Dataset(input_path + continuous_file, format="NETCDF4")

print('the last time is {} seconds from the start of the experiment'.format(continuous_netcdf_object['time'][-1]))
print('looking up the boundaries... this takes about 10 minutes on a 12 GB file')

edges = reshape.find_boundaries(continuous_netcdf_object['time'][:], slices)
for x in edges[0:5]:
    print('at indices {} to {} we found times {} to {}'.format(x[0], x[1], continuous_netcdf_object['time'][x[0]],
                                                               continuous_netcdf_object['time'][x[1]]))
burst_lengths = list(map(lambda t: t[1] - t[0], edges))
for x in burst_lengths[0:5]:
    print('bursts are {} long'.format(x))
    
continuous_netcdf_object.close()
print('elapsed time is {} min'.format((dt.datetime.now() - operation_start).total_seconds() / 60))

# TODO - this is not working when called from outside this file
reshape.save_indexes_to_file(input_path + continuous_file, edges, output_path + index_file_name)
    
number_of_bursts_per_file = int(math.floor(len(edges) / number_of_output_files))
# now iterate through the number of output files
# for ifile in range(1):
for file_index in range(number_of_output_files):
    s = burst_file_name.split('.')
    burstFile = s[0] + (f'%02d.' % file_index) + s[1]
    print('making burst file {}'.format(burstFile))
    
    burst_start_index = file_index * number_of_bursts_per_file
    burst_end_index = burst_start_index + number_of_bursts_per_file
    
    if burst_end_index > len(edges):
        enburst = len(edges)
        
    edges_this_file = edges[burst_start_index:burst_end_index]
    samples_in_each_burst = list(map(lambda t: t[1] - t[0] + 1, edges_this_file))
    
    # if there are no samples in a burst, we will skip the burst
    # skip them by removing them from this index list
    # this cleans up the tail end of the last file
    # TODO - use None to signal 
    idx_empty_bursts = list(map(lambda y: False if x == 0 else True, samples_in_each_burst))
    print('Zeros samples in {} bursts, these will be omitted'.format(idx_empty_bursts.count(0)))

    continuous_netcdf_object = nc.Dataset(input_path + continuous_file, format="NETCDF4")
    time_units = continuous_netcdf_object['time'].units
    number_to_display = 5
    if number_of_bursts_per_file < number_to_display or number_of_bursts_per_file < number_to_display*2:
        number_to_display = number_of_bursts_per_file
        x = list(range(number_to_display))
    else:
        x = list(range(number_to_display))+list(range(len(edges_this_file)-number_to_display-1, len(edges_this_file)-1))
    for i in x:
        print('burst {} will be {} samples from {} to {}'.format(
                i, samples_in_each_burst[i],
                num2date(continuous_netcdf_object['time'][edges_this_file[i][0]], time_units),
                num2date(continuous_netcdf_object['time'][edges_this_file[i][1]], time_units)))

    continuous_netcdf_object.close()

    if not dry_run:
        reshape.reshapeEPIC(input_path + continuous_file, output_path + burstFile, burst_length,
                            dim='time', edges=edges_this_file, drop=variables_to_omit,
                            variable_attributes_to_omit=attributes_to_omit)

print('End script run at ', dt.datetime.now())
print('elapsed time is {} min'.format((dt.datetime.now() - operation_start).total_seconds() / 60))
