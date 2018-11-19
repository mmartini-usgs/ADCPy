Known issues:

Have posted a fix for the following, will be testing on further data sets:  Still having trouble converting between time formats in python.  EPIC time therefore is not correct.  The problem has something to do with if EPIC_time are stored as u2 or u4.  In one case, the EPIC time is good until about 1/3 of the way through a long time series, and then for no apparent reason, jumps ahead about 6 months.  Need to discuss this with someone who knows more about netCDF data types as they are treated in python than I do.

It's possible for pressure standard deviation to be zero as computed by the instrument.  Detect this and omit it from the final .nc file if so.

This code will likely run very slowly when making the inital conversions.  May need to restructure to make use of xarray.  Try to do matrix transformations with numpy.einsum?

It is possible for the 5th beam of a V or Signature to have a different number of cells.  This is trapped, but not dealt with, in the code.  In this case 5th beam data will not be output.

There may be a case where there is only wave data in a Velocity output file.  The program will probably not handle this well.

Learn how to implement and output a log file.

Time values that are not unique are still being output for very fast sampled time series or where there are bottom track pings or other combined sampling.  This needs to be sorted out, possibly by implementing groups.

With some data sets, when time is written to netCDF, there are invalid values.  This seems to happen if Inf is used as the last ensemble to read, with the intent to read the whole file.  In these cases it is best to use an explicit ensemble count. 

EPIC time
-- is computing bad dates, I have been unable to find the problem.
-- needs to be removed.  It is an impediment for xarray and probably incompatible with stglib
-- CF time is OK for most uses, including tools such as ncbrowse, panoply and MATLAB

Help docstrings need to be improved

Need perform on raw data before performing rotations:
-- QA/QC thresholds 
-- masking
-- bin mapping

During rotations:
-- 3 beam solution

Post processing cleanup things to implement
-- post resample, shift time base to center of burst (resample default is probably beginning of hour or time interval)

Efficiency
-- read in buffers of data.  Since the unpacking module is the first thing I wrote in python, it is inefficient because it reads by the byte.  Change this to read with bigger buffers - though this may not be possible if there a multiple data types in each file with unknown sampling schemes.  Dolfyn (https://github.com/lkilcher/dolfyn) has reading modules that do most, but not all of the data reading this package needs.  Might be a good model for an update.
