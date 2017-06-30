6/30/2017 Caught that I was 10 sec off in my hundredths to microsecond conversion when manupulating time
5/9/2017 Put in code to handle bottom track ping data, have not yet tested on real data
2/28/2017 Needed time expressed as CF convention that xarray could understand.  Did this by expressing as seconds since the first timestamp in the file.
		Saving a version of the code before re-writing for another xarray requirement, that time be the index.
		Saving another version to include BT data and with time successfully written as an ordinate dimension.
		Now with time as an index, xarray can do its time magic and resample can be used.
2/21/2017 Got time figured out.  Using time as a variable, ensemble as the record variable.  Had to define the time variable in netCDF as f8 (64 bit floating point) to get teh necessary resolution for milliseconds in a day.
2/7/2017 Still having issues getting time right.  Went back to time and time2 output with rec as the longest dimension.  Worked with very large (>3mil enembles, 6 GB netcdf file size) files.
2/2/2017 Change longest dimension to 'Rec' rather than time.  Some ADCP records do not have unique time values.  Remove time2.  Make time a julian day.
2/1/2017 Added ability to parse Velocity output data, 5 beam
1/30/2017 Fixed tilts, was reading unsigned, rather than signed shorts
1/26/2017 Can now read through, though not yet output, data exported by Velocity with the new V series output, though somethign causes the program to exit ungracefully
1/25/2017 Got a completed version of a program that reads old workhorse data, currents only, pretty well.  Need to make some decisions about what units to use - convert or leave as found in the TRDI format