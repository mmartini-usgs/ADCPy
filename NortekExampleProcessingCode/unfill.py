# We clean up cruft from resample.  No _FillValues for dimensions or 
# coordinate variables, especially not NaN!
#import xarray as xr
from netCDF4 import Dataset
import numpy as np

datapath = "E:\\data\\MVCO15\\10573_Signature\\Nortekncsmall\\"

datafiles = ['10573sig000_1h.nc',\
             '10573sig001_1h.nc',\
             '10573sig002_1h.nc',\
             '10573sig003_1h.nc']

for filenum in range(len(datafiles)):
    
    pydata = datapath+datafiles[filenum]
        
    # check _FillValue
    print('Before, as netcdf Dataset')
    ncdata = Dataset(pydata, mode="r", format='NETCDF4')
    print(ncdata['lat'].__dict__)
    ncdata.close()

    # re-open the dataset for numerical operations such as min and max
    # we have to make attribute changes, etc. so need to open with the netCDF package
    pyd = Dataset(pydata, mode="r+", format='NETCDF4')

    # fill value cleanup
    # remove any fill values from dimensions and coordinate variables
    print('Dimensions:')
    for key in pyd.dimensions.keys():
        fillok = True
        try:
            fill = np.asarray(pyd[key].getncattr('_FillValue'))
            print('%s has fill %g, removing attribute' % (key,fill))
            pyd[key].delncattr('_FillValue')
        except:
            print('%s has no _FillValue defined' % key) 
            fillok = False
    
    # make sure that all Fill Values match variable type and are not NaN
    print('Coordinate Variables:')
    for key in pyd.dimensions.keys():
        fillok = True
        try:
            fill = np.asarray(pyd[key].getncattr('_FillValue'))
            print('%s has fill %g' % (key,fill))
            if np.isnan(fill):
                print('nan detected')
                np.asarray(pyd[key].getncattr('_FillValue'))
                newfill = 1e35
                #print(type(newfill.astype(np.float)))
                print(newfill)
                print(type(newfill))
                pyd[key].attrs.update({'_FillValue':1E35})
        except:
            print('%s has no _FillValue defined' % key) 
            fillok = False
    
    pyd.close()
    
    # check _FillValue
    print('After, as netcdf Dataset')
    ncdata = Dataset(pydata, mode="r", format='NETCDF4')
    print(ncdata['lat'].__dict__)
    ncdata.close()
