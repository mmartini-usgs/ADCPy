# Here we test the conversion of a pd0 file to netcdf format.
#
# these tests are designed to follow the example here
# http://doc.pytest.org/en/latest/fixture.html
import pytest
from pathlib import Path

here = Path(__file__).parent

@pytest.fixture
def trdipd0_convert_pd0file():
    import adcpy.TRDIstuff.TRDIpd0tonetcdf as TRDIpd0
    import numpy as np
    # this is a test data file, very simple ADCP file
    # of older Workhorse ADCP data with only profiler
    # you will use assert ensemble_count == 499 for this file
    # pd0_input_file = here.joinpath("9991wh000.000")
    # netcdf_output = here.joinpath("9991wh000.cdf")
    # serial_number = "473"
    # delta_t = "900"  # seconds, as a string
    # this is a test data file of a larger data set of Sentinel V data
    # you will use assert ensemble_count == 8430 for this file
    pd0_input_file = here.joinpath("11121whV23857profiles.pd0")
    serial_number = "23857"
    delta_t = "3600"  # seconds, as a string
    netcdf_output = here.joinpath("11121whV23857profiles.cdf")
    time_type = "CF"  # the time variable will have a CF time format
    ensembles_to_process = [265, np.inf]
    ensemble_count, netcdf_index, read_error = TRDIpd0.convert_pd0_to_netcdf(pd0_input_file, netcdf_output,
                                                                             ensembles_to_process,
                                                                             serial_number, time_type, delta_t)
    print(ensemble_count)
    return ensemble_count, netcdf_index, read_error


def test__output(trdipd0_convert_pd0file):
    ensemble_count, netcdf_index, read_error = trdipd0_convert_pd0file
    assert ensemble_count == 8430
    assert read_error is None

