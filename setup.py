# this setup modelled on 
# https://github.com/pypa/sampleproject/blob/master/setup.py
from setuptools import setup, find_packages
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    readme = f.read()

with open('LICENSE.md') as f:
    license = f.read()

setup(
    name='ADCPy',
    version='0.0.0',
	# TODO implement versioneer here
	# version=versioneer.get_version(),
    # cmdclass=versioneer.get_cmdclass(),
    description=('code to work with ADCP data from the raw binary' 
	    ' in python 3x, a learning project'),
    long_description=readme, # read from README.md above
    long_description_content_type='text/markdown',
    url='https://github.com/mmartini-usgs/ADCPy',
    author='Marinna Martini',
    author_email='mmartini@usgs.gov',
	license='Public domain',
    classifiers=['Development Status :: 3 - Alpha',
                 'Intended Audience :: Science/Research',
                 'Programming Language :: Python',
                 'Programming Language :: Python :: 3.7'
				 ],
    keywords='acoustic doppler profiler ADCP',  
    packages=find_packages(exclude=('tests', 'docs')),
	python_requires='>=3.5',
	# install_requires=['peppercorn'],  # Optional, here as a reminder
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


