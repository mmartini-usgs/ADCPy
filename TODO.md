Known issues:

This code will likely run very slowly.  Need to restructure to make use of xarray.  Need to do matrix transformations with numpy.einsum

It is possible for the 5th beam of a V to have a different number of cells.  This is trapped, but not dealt with, in the code.  In this case 5th beam data will nto be output.

There may be a case where there is only wave data in a Velocity output file.  The program will probably not handle this well.

Learn how to impliment and output a log file.

Time values that are not unique are still being output for very fast sampled time series or where there are bottom track pings or other combined sampling.  This needs to be sorted out, possibly by implmenting groups.

With some data sets, when time is written to netCDF, there are invalid values.  This seems to happen if Inf is used as the last ensemble to read, with the intent to read the whole file.  In these cases it is best to use an explicit ensemble count. 