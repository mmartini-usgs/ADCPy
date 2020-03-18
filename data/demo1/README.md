Workhorse Sentinel data

  - The raw binary ADCP file is 9991wh000.000
  - The text file with standard USGS metadata is glob_att999.txt This file is read and the contents added as global 
    attributes to the netCDF file.  Some of these attributes (like water_depth) are necessary for processing and 
    converting the data.

You will need to decide what type of time convention to use, this code will support EPIC or CF.
	
Data are processed in several steps

  - Part 1: Conversion from raw binary to netCDF is done by convert_pd0_to_netcdf() in adcpy.TRDIstuff.TRDIpd0tonetcdf,
    and this file is exclusive to TRDI data.  Data are output in raw form, as a netCDF file in py9991wh.cdf
  - Part 2: Calculate rotation from beam or instrument coordinates to East, North and Up, with doEPIC_ADCPfile() in
    adcpy.EPICstuff.ADCPcdf2ncEPIC.  Data are output to a netCDF in EPIC conventions file as py9991wh.nc.  These data 
    are not cleaned or edited in any way to trim bins out of the water or noisy data.

The files py9991wh.cdf and py9991.nc are provided as output examples.

The file 9991adcp.doc is the summary of field notes and metadata for this instrument's deployment.

The output files can be examined and data plotted in puython using the netCDF4 package.  If CF time was specified, 
the xarray package can be used. 


