Currently in work:
- docstrings, trying to write legible docstrings for sphinx and docstring parsing so that the cods are only in one place, at the function definition or file comments.  The result is OK, and not perfect for either.
- writing more tests
- bringing code up to PEP8 standards

Known issues:

It's possible for pressure standard deviation to be zero as computed by the instrument.  Detect this and omit it from the final .nc file if so.

This code will likely run very slowly when making the inital conversions.  May need to restructure to make use of xarray.  Try to do matrix transformations with numpy.einsum?  Either way, it is possible for ensemble types to vary in the middle of the file, and thus, each ensemble must be examined wthout assumption, this is a slow process.  This code also keeps only one ensemble at a time in memory.

It is possible for the 5th beam of a V or Signature to have a different number of cells.  This is trapped, but not dealt with, in the code.  In this case 5th beam data will not be output.

There may be a case where there is only wave data in a Velocity output file.  The program will probably not handle this well.

Learn how to implement and output a log file.

Time values that are not unique are still being output for very fast sampled time series or where there are bottom track pings or other combined sampling.  This needs to be sorted out, possibly by implementing groups.

With some data sets, when time is written to netCDF, there are invalid values.  This seems to happen if Inf is used as the last ensemble to read, with the intent to read the whole file.  In these cases it is best to use an explicit ensemble count. 

Need perform on raw data before performing rotations:
-- QA/QC thresholds 
-- masking
-- bin mapping

During rotations:
-- 3 beam solution

Post processing cleanup things to implement
-- post resample, shift time base to center of burst (resample default is probably beginning of hour or time interval)

Efficiency
-- read in buffers of data.  Since the unpacking module is the first thing I wrote in python, it is inefficient because it reads by the byte.  Change this to read with bigger buffers - though this may not be possible if there a multiple data types in each file with unknown sampling schemes.  Dolfyn (https://github.com/lkilcher/dolfyn) has reading modules that do most, but not all of the data reading this package needs, however does not save to netCDF and  think keeps all data in memory.  Might be a good model for an update.
