# -*- coding: utf-8 -*-

from setuptools import setup, find_packages


with open('README.md') as f:
    readme = f.read()

with open('LICENSE.md') as f:
    license = f.read()

setup(
    name='ADCPy',
    version='0.0.0',
    description='code to work with ADCP data from the raw binary in python 3x',
    long_description=readme,
    author='Marinna Martini',
    author_email='mmartini@usgs.gov',
    url='https://github.com/mmartini-usgs/ADCPy',
    license=license,
    packages=find_packages(exclude=('tests', 'docs'))
)
