Known issues:

Still having trouble converting between time formats in python.  EPIC time therefore is not correct.  The problem has something to do with if EPIC_time are stored as u2 or u4.  In one case, the EPIC time is good until about 1/3 of the way through a long time series, and then for no apparent reason, jumps ahead about 6 months.  Need to discuss this with someone who knows more about netCDF data types as they are treated in python than I do.

It's possible for pressure standard deviation to be zero as computed by the instrument.  Detect this and omit it from the final .nc file if so.

This code will likely run very slowly when making the inital conversions.  May need to restructure to make use of xarray.  Try to do matrix transformations with numpy.einsum?

It is possible for the 5th beam of a V or Signature to have a different number of cells.  This is trapped, but not dealt with, in the code.  In this case 5th beam data will not be output.

There may be a case where there is only wave data in a Velocity output file.  The program will probably not handle this well.

Learn how to impliment and output a log file.

Time values that are not unique are still being output for very fast sampled time series or where there are bottom track pings or other combined sampling.  This needs to be sorted out, possibly by implmenting groups.

With some data sets, when time is written to netCDF, there are invalid values.  This seems to happen if Inf is used as the last ensemble to read, with the intent to read the whole file.  In these cases it is best to use an explicit ensemble count. 

Need perform on raw data before performing rotations:
-- QA/QC thresholds 
-- masking
-- bin mapping

During rotations:
-- 3 beam solution

Post processing cleanup things to implement
-- trim to surface, by following surface or by bin
-- trim in time
-- swap between EPIC and CF time base
-- post resample, shift time base to center of burst (resample default is probably beginning of hour or time interval)

For wave processing
-- reshape into bursts