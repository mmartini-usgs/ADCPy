# this setup modelled on 
# https://github.com/pypa/sampleproject/blob/master/setup.py
import setuptools

# Get the long description from the README file
with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='ADCPy',
    version='0.0.3',
    author='Marinna Martini',
    author_email='mmartini@usgs.gov',
    description='read ADCP data from TRDI and Nortek instruments',
    long_description=long_description, # read from README.md above
    long_description_content_type='text/markdown',
    url='https://github.com/mmartini-usgs/ADCPy',
    packages=setuptools.find_packages(exclude=('tests', 'docs')),
    classifiers=['Programming Language :: Python :: 3',
                 'License :: Public Domain',
                 'Operating System :: OS Independent',
                 'Development Status :: 3 - Alpha',
                 'Intended Audience :: Science/Research',
                 ],
    python_requires='>=3.5',
    keywords='acoustic doppler profiler ADCP',
)


# TODO - include data for demos
#    # If there are data files included in your packages that need to be
#    # installed, specify them here.
#    #
#    # If using Python 2.6 or earlier, then these have to be included in
#    # MANIFEST.in as well.
#    package_data={  # Optional
#        'sample': ['package_data.dat'],
#    },
#
#    # Although 'package_data' is the preferred approach, in some case you may
#    # need to place data files outside of your packages. See:
#    # http://docs.python.org/3.4/distutils/setupscript.html#installing-additional-files
#    #
#    # In this case, 'data_file' will be installed into '<sys.prefix>/my_data'
#    data_files=[('my_data', ['data/data_file'])],  # Optional


