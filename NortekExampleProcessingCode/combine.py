# simply run catEPIC
import sys
# this is important in order to import my package which is not on the python path
sys.path.append('c:\projects\python\ADCPy')
from EPICstuff import catEPIC

path = '.\\'
rootname = '1108sig'
nfiles = 24
do_nc_files = True
do_cdf_files = True

if do_nc_files:
    # combine the .nc files
    datafiles = []
    for i in range(1,nfiles+1):
        datafiles.append(path + ('%s%03d.nc' % (rootname,i)))
    
    print(datafiles)
    
    outfile = path + rootname + 'all.nc'

    catEPIC(datafiles, outfile)
    
if do_cdf_files:
    # combine the .cdf files
    datafiles = []
    for i in range(1,nfiles+1):
        datafiles.append(path + ('%s%03d.cdf' % (rootname,i)))
    
    print(datafiles)
    
    outfile = path + rootname + 'all.cdf'
    
    catEPIC(datafiles, outfile)
