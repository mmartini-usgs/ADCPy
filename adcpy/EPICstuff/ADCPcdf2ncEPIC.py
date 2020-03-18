"""
This code takes a raw netcdf file containing data from any 4 beam Janus 
acoustic doppler profiler, with or without a center beam, and transforms the
data into Earth coordinates.  Data are output to netCDF using controlled
vocabulary for the variable names, following the EPIC convention wherever
possible.  

ADCPcdf2ncEPIC.doEPIC_ADCPfile(cdfFile, ncFile, attFile, settings)

cdfFile = path to a USGS raw netCDF ADCP data file

ncFile = a netcdf file structured according to PMEL EPIC conventions

attFile = a file containing global attributes (metadata) for the data.  See below

settings = a dictionary of preferences for the processing::

    'good_ensembles': [0, np.inf] # starting and ending indices of the input file. For all data use [0,np.inf]
    'orientation': 'UP' # uplooking ADCP, for downlooking, use DOWN
    'transducer_offset_from_bottom': 1.0 # a float in meters
    'transformation': 'EARTH' # | BEAM | INST
    'adjust_to_UTC': 5 # for EST to UTC, if no adjustment, set to 0 or omit
    
Depth dependent attributes are compute from the mean Pressure found in the raw
data file.  So it is best to have the time series trimmed to the in water
time or to provide the good ensemble indices for in water time

Note that file names and paths may not include spaces

Example contents of a Global Attribute file::

    SciPi; J.Q. Scientist
    PROJECT; USGS Coastal Marine Geology Program
    EXPERIMENT; MVCO 2015 Stress Comparison
    DESCRIPTION; Quadpod 13.9m
    DATA_SUBTYPE; MOORED
    COORD_SYSTEM; GEOGRAPHIC + SAMPLE
    Conventions; PMEL/EPIC
    MOORING; 1057
    WATER_DEPTH; 13.9
    WATER_DEPTH_NOTE; (meters), nominal
    WATER_DEPTH_source; ship fathometer
    latitude; 41.3336633
    longitude; -70.565877
    magnetic_variation; -14.7
    Deployment_date; 17-Nov-2015
    Recovery_date; 14-Dec-2015
    DATA_CMNT;  
    platform_type; USGS aluminum T14 quadpod
    DRIFTER; 0
    POS_CONST; 0
    DEPTH_CONST; 0
    Conventions; PMEL/EPIC
    institution; United States Geological Survey, Woods Hole Coastal and Marine Science Center
    institution_url; http://woodshole.er.usgs.gov

Created on Tue May 16 13:33:31 2017

@author: mmartini
"""

import os
import sys
import numpy as np 
from netCDF4 import Dataset
import netCDF4 as netcdf
import datetime as dt
from datetime import datetime


# noinspection PyPep8Naming
def doEPIC_ADCPfile(cdfFile, ncFile, attFile, settings):
    """
    Convert a raw netcdf file containing data from any 4 beam Janus acoustic doppler profiler,
    with or without a center beam, and transforms the data into Earth coordinates.  Data are output to netCDF
    using controlled vocabulary for the variable names, following the EPIC convention wherever possible.

    :param str cdfFile:  raw netCDF input data file name
    :param str ncFile:   output file name
    :param str attFile:  text file containing metadata
    :param dict settings: a dict of settings as follows::

        'good_ensembles': [] # starting and ending indices of the input file. For all data use [0,np.inf]
        'orientation': 'UP' # uplooking ADCP, for downlooking, use DOWN
        'transducer_offset_from_bottom': 2.02 # in meters
        'transformation': 'EARTH' # | BEAM | INST
        'adjust_to_UTC': 5 # for EST to UTC, if no adjustment, set to 0 or omit

    """
    # check some of the settings we can't live without
    # set flags, then remove from the settings list if we don't want them in metadata
    if 'good_ensembles' not in settings.keys():
        settings['good_ensembles'] = [0, np.inf]  # nothing from user, do them all
        print('No starting and ending ensembles specfied, processing entire file')
    if 'orientation' not in settings.keys():
        settings['orientation'] = "UP"
        settings['orientation_note'] = "assumed by program"
        print('No orientation specfied, assuming up-looking')
    else:
        settings['orientation_note'] = "user provided orientation"
    if 'transducer_offset_from_bottom' not in settings.keys():
        settings['transducer_offset_from_bottom'] = 0
        print('No transducer_offset_from_bottom, assuming 0')
    if 'transformation' not in settings.keys():
        settings['transformation'] = "EARTH"
    if 'adjust_to_UTC' not in settings.keys():
        settings['adjust_to_UTC'] = 0
    # TODO implement this time_type_out selection, right now does what is in the raw netCDF
    # if 'time_type_out' not in settings.keys():
    #    settings['time_type_out'] = "CF"
    # if 'time_type_out' in settings.keys():
    #    time_type_out = settings['time_type_out']
    # else:
    #    time_type_out = "CF"
    if 'use_pressure_for_WATER_DEPTH' in settings.keys():
        if settings['use_pressure_for_WATER_DEPTH']:
            usep4waterdepth = True
            settings.pop('use_pressure_for_WATER_DEPTH')
        else:
            usep4waterdepth = False
    else:
        usep4waterdepth = True
        
    rawcdf = Dataset(cdfFile, mode='r', format='NETCDF4')
    rawvars = []
    for key in rawcdf.variables.keys():
        rawvars.append(key)
      
    # this function will operate on the files using the netCDF package
    nc = setupEPICnc(ncFile, rawcdf, attFile, settings)
    
    nbeams = nc.number_of_slant_beams  # what if this isn't 4?
    nbins = len(rawcdf.dimensions['depth'])
    nens = len(rawcdf.dimensions['time'])
    ncvars = []
    for key in nc.variables.keys():
        ncvars.append(key)

    declination = nc.magnetic_variation_at_site

    # start and end indices
    s = settings['good_ensembles'][0]
    if settings['good_ensembles'][1] < 0:
        e = nens
    else:
        e = settings['good_ensembles'][1]
    print('Converting from index %d to %d of %s' % (s, e, cdfFile))

    # many variables do not need processing and can just be copied to the
    # new EPIC convention
    varlist = {'sv': 'SV_80', 'Rec': 'Rec'}
    
    for key in varlist:
        varobj = nc.variables[varlist[key]]
        varobj[:] = rawcdf.variables[key][s:e]  
        
    # check the time zone, Nortek data are usually set to UTC, no matter what
    # the actual time zone of deployment might have been
    if abs(settings['adjust_to_UTC']) > 0:
        nc.time_zone_change_applied = settings['adjust_to_UTC']
        nc.time_zone_change_applied_note = "adjust time to UTC requested by user"    
    toffset = settings['adjust_to_UTC']*3600

    # determine what kind of time setup we have in the raw file
    timevars = ['time', 'time2', 'EPIC_time', 'EPIC_time2', 'cf_time']
    timevars_in_file = [item for item in timevars if item in rawvars]

    if timevars_in_file == ['time', 'time2']:
        time_type = "EPIC"
    elif timevars_in_file == ['time', 'time2', 'cf_time']:
        time_type = "EPIC_with_CF"
    elif timevars_in_file == ['time', 'EPIC_time', 'EPIC_time2']:
        time_type = "CF_with_EPIC"
    elif timevars_in_file == ['time']:
        time_type = "CF"
    else:
        time_type = None
        print("Unrecognized time arrangement, known variables found: {}".format(timevars_in_file))

    print("The raw netCDF file has time_type {}".format(time_type))

    # raw variable name : EPIC variable name
    if settings['time_type_out'] == 'EPIC':
        varlist = {'time': 'time', 'time2': 'time2'}
    elif settings['time_type_out'] == 'CF_with_EPIC':
        varlist = {'time': 'time', 'EPIC_time': 'EPIC_time', 'EPIC_time2': 'EPIC_time2'}
    elif settings['time_type_out'] == 'EPIC_with_CF':
        varlist = {'time': 'time', 'time2': 'time2', 'cf_time': 'cf_time'}
    else:  # only CF time, the default
        varlist = {'time': 'time'}

    # TODO let user select type of time output, right now it uses what is in the netCDF file
    for key in varlist:
        print(key)
        varobj = nc.variables[varlist[key]]
        varobj[:] = rawcdf.variables[key][s:e]+toffset
        
    # TRDI instruments have heading, pitch, roll and temperature in hundredths of degrees
    if rawcdf.sensor_type == "TRDI":
        degree_factor = 100
    else:
        degree_factor = 1
    
    varlist = {'Ptch': 'Ptch_1216', 'Roll': 'Roll_1217', 'Tx': 'Tx_1211'}
    for key in varlist:
        varobj = nc.variables[varlist[key]]
        varobj[:] = rawcdf.variables[key][s:e]/degree_factor 
              
    #  TODO will need an instrument dependent methodology to check for any previous adjustments to heading
    #  prior to this correction. for instance, with TRDI instruments, Velocity or the EB command might have applied
    #  a correction.  If EB is set, then that value was applied to the raw data seen by TRDIpd0tonetcdf.py
    nc.magnetic_variation_applied = declination
    nc.magnetic_variation_applied_note = "as provided by user"
    heading = rawcdf.variables['Hdg'][s:e]/degree_factor + declination
    heading[heading >= 360] = heading[heading >= 360] - 360
    heading[heading < 0] = heading[heading < 0] + 360
    nc['Hdg_1215'][:] = heading
    
    # pressure needs to be in db or m
    if 'Pressure' in rawvars:
        pconvconst = 1  # when in doubt, do nothing
        punits = rawcdf['Pressure'].units
        if 'deca-pascals' in punits:
            pconvconst = 1000  # decapascals to dbar = /1000
            print('Pressure in deca-pascals will be converted to db')
        nc['P_1'][:] = rawcdf.variables['Pressure'][s:e]/pconvconst
    
    # check units of current velocity and convert to cm/s
    vunits = rawcdf['vel1'].units
    vconvconst = 1  # when in doubt, do nothing
    if (vunits == 'mm s-1') | (vunits == 'mm/s'):
        vconvconst = 0.1  # mm/s to cm/s
    elif (vunits == 'm s-1') | (vunits == 'm/s'):
        vconvconst = 100  # m/s to cm/s

    print('Velocity in {} will be converted using a multiplier of {}'.format(vunits, vconvconst))

    if 'vel5' in rawvars:
        nc['Wvert'][:] = rawcdf.variables['vel5'][s:e, :] * vconvconst

    if 'cor5' in rawvars:
        nc['corvert'][:] = rawcdf.variables['cor5'][s:e, :]
    
    if 'att5' in rawvars:
        nc['AGCvert'][:] = rawcdf.variables['att5'][s:e, :]

    if 'PGd4' in rawvars:
        nc['PGd_1203'][:, :, 0, 0] = rawcdf.variables['PGd4'][s:e, :]
        
    if 'PressVar' in rawvars:
        nc['SDP_850'][:, 0, 0] = rawcdf.variables['PressVar'][s:e]
    
    bindist = np.arange(len(nc['bindist']))
    bindist = bindist*nc.bin_size+nc.center_first_bin
    nc['bindist'][:] = bindist

    # figure out DELTA_T - we need to use the cf_time, more convenient
    if settings['time_type_out'] == 'CF':
        # we will set the main "time' variable to CF convention
        timekey = 'time'
    else:
        timekey = 'cf_time'
    
    # this calculation uses for CF time
    dtime = np.diff(nc[timekey][:])
    delta_t = '%s' % int((dtime.mean().astype('float')).round())  # needs to be a string
    nc.DELTA_T = delta_t
    
    # depths and heights
    nc.initial_instrument_height = settings['transducer_offset_from_bottom']
    nc.initial_instrument_height_note = "height in meters above bottom: accurate for tripod mounted instruments" 
    # compute depth, make a guess we want to average all depths recorded 
    # deeper than user supplied water depth
    # idx is returned as a tuple, the first of which is the actual index values

    # set the water depth here, this will be used throughout
    # the user may have put units next to the depth
    if type(nc.WATER_DEPTH) is str:
        water_depth = nc.WATER_DEPTH.split()
        water_depth = float(water_depth[0])
    else:
        water_depth = nc.WATER_DEPTH
    if ('Pressure' in rawvars) and usep4waterdepth:
        idx = np.where(nc['P_1'][:] > water_depth/2)
        # now for the mean of only on bottom pressure measurements
        if len(idx[0]) > 0:
            pmean = nc['P_1'][idx[0]].mean()
        else:
            pmean = 0  # this could be if the ADCP is in air the whole time
        print('Site WATER_DEPTH given is %f' % water_depth)
        print('Calculated mean water level from P_1 is %f m' % pmean)
        print('Updating site WATER_DEPTH to %f m' % pmean)
        nc.WATER_DEPTH = pmean+nc.transducer_offset_from_bottom
        nc.WATER_DEPTH_source = "water depth = MSL from pressure sensor, (meters), nominal"
        nc.WATER_DEPTH_NOTE = nc.WATER_DEPTH_source
        nc.nominal_sensor_depth = nc.WATER_DEPTH-settings['transducer_offset_from_bottom']
        nc.nominal_sensor_depth_note = "inst_depth = (water_depth - inst_height); nominal depth below surface, meters"
        varnames = ['bindist', 'depth']
        # WATER_DEPTH_datum is not used in this circumstance.
    else:
        print('Site WATER_DEPTH given is %f' % water_depth)
        print('No pressure data available, so no adjustment to water depth made')
        nc.WATER_DEPTH_source = "water depth as given by user, (meters), nominal"
        nc.WATER_DEPTH_NOTE = nc.WATER_DEPTH_source
        nc.nominal_sensor_depth = water_depth-settings['transducer_offset_from_bottom']
        nc.nominal_sensor_depth_note = "inst_depth = (water_depth - inst_height); nominal depth below surface, meters"
        varnames = ['bindist', 'depth']
        # WATER_DEPTH_datum is not used in this circumstance.

    for varname in varnames:
        nc[varname].WATER_DEPTH = water_depth
        nc[varname].WATER_DEPTH_source = nc.WATER_DEPTH_source
        nc[varname].transducer_offset_from_bottom = nc.transducer_offset_from_bottom
    
    # update depth variable for location of bins based on WATER_DEPTH information
    if "UP" in nc.orientation:
        depths = water_depth-nc.transducer_offset_from_bottom-nc['bindist']
    else:
        depths = -1 * (water_depth-nc.transducer_offset_from_bottom+nc['bindist'])
    
    nc['depth'][:] = depths        
    
    nc.start_time = '%s' % netcdf.num2date(nc[timekey][0], nc[timekey].units)
    nc.stop_time = '%s' % netcdf.num2date(nc[timekey][-1], nc[timekey].units)
    
    # some of these repeating attributes depended on depth calculations
    # these are the same for all variables because all sensors are in the same
    # package, as of now, no remote sensors being logged by this ADCP
    ncvarnames = []
    for key in nc.variables.keys():
        ncvarnames.append(key)
    omitnames = []
    for key in nc.dimensions.keys():
        omitnames.append(key)
    omitnames.append("Rec")
    omitnames.append("depth")
    for varname in ncvarnames:
        if varname not in omitnames:
            varobj = nc.variables[varname]
            varobj.sensor_type = nc.INST_TYPE
            varobj.sensor_depth = nc.nominal_sensor_depth
            varobj.initial_sensor_height = nc.initial_instrument_height 
            varobj.initial_sensor_height_note = "height in meters above bottom:  " +\
                                                "accurate for tripod mounted instruments"
            varobj.height_depth_units = "m"

    print('finished copying data, starting computations at %s' % (dt.datetime.now()))

    print('averaging cor at %s' % (dt.datetime.now()))
    # this will be a problem - it loads all into memory
    cor = (rawcdf.variables['cor1'][s:e, :] + rawcdf.variables['cor2'][s:e, :] +
           rawcdf.variables['cor3'][s:e, :] + rawcdf.variables['cor4'][s:e, :]) / 4
    nc['cor'][:, :, 0, 0] = cor[:, :]

    print('averaging AGC at %s' % (dt.datetime.now()))
    # this will be a problem - it loads all into memory
    agc = (rawcdf.variables['att1'][s:e, :] + rawcdf.variables['att2'][s:e, :] +
           rawcdf.variables['att3'][s:e, :]+rawcdf.variables['att4'][s:e, :]) / 4
    nc['AGC_1202'][:, :, 0, 0] = agc[:, :]
    
    print('converting %d ensembles from beam to earth %s' % (len(nc[timekey]), dt.datetime.now()))
    
    # check our indexing
    print('magnetic variation at site = %f' % nc.magnetic_variation_at_site)
    print('magnetic variation applied = %f' % nc.magnetic_variation_applied)
    print('magnetic variation applied note = %s' % nc.magnetic_variation_applied_note)
    n = int(len(heading)/2)
    print('From the middle of the time series at ensemble #%d, we have:' % n)
    print('heading variable in this python process = %f' % heading[n])
    print('rawcdf Hdg[n] = %f' % rawcdf['Hdg'][n])
    print('nc Hdg_1215[n] = %f' % nc['Hdg_1215'][n, 0, 0])
    
    # TODO add depth bin mapping
    
    # this beam arrangement is for TRDI Workhorse and V, other instruments 
    # should be re-ordered to match
    rawvarnames = ["vel1", "vel2", "vel3", "vel4"]
    ncidx = 0
    if settings['transformation'].upper() == "BEAM":
        ncvarnames = ["Beam1", "Beam2", "Beam3", "Beam4"]        
        for idx in range(s, e):
            for beam in range(nbeams):
                nc[ncvarnames[beam]][ncidx, :, 0, 0] = \
                    rawcdf.variables[rawvarnames[beam]][idx, :] * vconvconst
            ncidx = ncidx + 1
            
    elif (settings['transformation'].upper() == "INST") or (settings['transformation'].upper() == "EARTH"):
        ncvarnames = ["X", "Y", "Z", "Error"]
        # the dolfyn way (https://github.com/lkilcher/dolfyn)
        # load the ADCP data object - we have converted this from a class object to nested dictionaries for use here
        adcpo = {
            'props': {
                'coord_sys': "beam",
                'inst2earth:fixed': False,
            },
            'config': {
                'beam_angle': nc.beam_angle,
                'beam_pattern': nc.beam_pattern,
                'orientation': nc.orientation,
            },
            # note declination is applied immediately when heading is read from the raw data file
            'declination_in_heading': True,
            # dolfyn shape for ensemble data is [bins x beams x ens]
            'vel': np.ones([nbins, nbeams], dtype='float') * np.nan,
        }
        # vels has to be pre-defined to get the shapes to broadcast
        # noinspection PyUnusedLocal
        vels = np.ones([nbins, 1], dtype='float') * np.nan
        
        # Nortek and TRDI do their along beam velocity directions opposite for
        # slant beams.  Vertical beam directions are the same.
        if rawcdf.sensor_type == 'Nortek':
            beam_vel_multiplier = -1
        else:
            beam_vel_multiplier = 1
        
        for idx in range(s, e):
            for beam in range(nbeams):
                # load data of one ensemble to dolfyn shape, in cm/s
                # adcpo['vel'][:,beam,0] = rawcdf.variables[rawvarnames[beam]][idx,:] * 0.1
                vels = rawcdf.variables[rawvarnames[beam]][idx, :] * vconvconst * beam_vel_multiplier
                adcpo['vel'][:, beam] = vels
            
            # need to keep setting this with new beam data since we are iterating
            adcpo['props']['coord_sys'] = "beam" 
            beam2inst(adcpo)  # adcpo['vel'] is returned in inst coordinates
        
            if settings['transformation'].upper() == "EARTH":
                ncvarnames = ["u_1205", "v_1206", "w_1204", "Werr_1201"]
                adcpo['heading_deg'] = nc.variables['Hdg_1215'][ncidx]
                adcpo['pitch_deg'] = nc.variables['Ptch_1216'][ncidx]
                adcpo['roll_deg'] = nc.variables['Roll_1217'][ncidx]
                inst2earth(adcpo)
                
            for beam in range(nbeams):
                nc[ncvarnames[beam]][ncidx, :, 0, 0] = adcpo['vel'][:, beam]
                   
            ncidx = ncidx + 1
            
            # immediate - then less feedback
            ensf, ensi = np.modf(ncidx/1000)
            if (ensf == 0) and (ncidx < 10000):
                print('%d of %d ensembles read' % (ncidx, nens))
            else:
                ensf, ensi = np.modf(ncidx/10000)
                if ensf == 0:
                    print('%d of %d ensembles read' % (ncidx, nens))
            
    nc.transform = settings['transformation'].upper()
    
    print('closing files at %s' % (dt.datetime.now()))

    rawcdf.close()
    nc.close()            


def cal_earth_rotmatrix(heading=0, pitch=0, roll=0, declination=0):
    """
    this transformation matrix is from the R.D. Instruments Coordinate Transformation booklet.
    It presumes the beams are in the same position as RDI Workhorse ADCP beams, where,
    when looking down on the transducers::

        Beam 3 is in the direction of the compass' zero reference
        Beam 1 is to the right
        Beam 2 is to the left
        Beam 4 is opposite beam 3
        Pitch is about the beam 2-1 axis and is positive when beam 3 is raised
        Roll is about the beam 3-4 axis and is positive when beam 2 is raised
        Heading increases when beam 3 is rotated towards beam 1

    Nortek Signature differs in these ways::

        TRDI beam 3 = Nortek beam 1
        TRDI beam 1 = Nortek beam 2
        TRDI beam 4 = Nortek beam 3
        TRDI beam 2 = Nortek beam 4
        Heading, pitch and roll behave the same as TRDI

    :param float heading: ADCP heading in degrees
    :param float pitch: ADCP pitch in degrees
    :param float roll: ADCP roll in degrees
    :param float declination: heading offset from true, Westerly is negative
    :return:
    """
    heading = heading + declination
    ch = np.cos(heading)
    sh = np.sin(heading)
    cp = np.cos(pitch)
    sp = np.sin(pitch)
    cr = np.cos(roll)
    sr = np.sin(roll)
    
    return np.asmatrix(np.array([
                                [(ch*cr+sh*sp*sr),  (sh*cp), (ch*sr-sh*sp*cr)],
                                [(-sh*cr+ch*sp*sr), (ch*cp), (-sh*sr-ch*sp*cr)],
                                [(-cp*sr),          sp,      (cp*cr)]
                                ]))

# start code from dolfyn


def calc_beam_rotmatrix(theta=20, convex=True, degrees=True):
    """
    Calculate the rotation matrix from beam coordinates to
    instrument head coordinates.
    per dolfyn rotate.py code here: https://github.com/lkilcher/dolfyn

    :param float theta: is the angle of the heads (usually 20 or 30 degrees)
    :param int convex: is a flag for convex or concave head configuration.
    :param bool degrees: is a flag which specifies whether theta is in degrees or radians (default: degrees=True)

    """
    deg2rad = np.pi / 180.
    if degrees:
        theta = theta * deg2rad
    if convex == 0 or convex == -1:
        c = -1
    else:
        c = 1
    a = 1 / (2. * np.sin(theta))
    b = 1 / (4. * np.cos(theta))
    d = a / (2. ** 0.5)
    return np.array([[c * a, -c * a, 0, 0],
                     [0, 0, -c * a, c * a],
                     [b, b, b, b],
                     [d, d, -d, -d]])


def _cat4rot(tpl):
    # TODO do we need this function _cat4rot
    """
    helper function

    :param tpl:
    :return: numpy array
    """
    tmp = []
    for vl in tpl:
        tmp.append(vl[:, None, :])
    return np.concatenate(tuple(tmp), axis=1)


def beam2inst(adcpo, reverse=False, force=False):
    """
    Rotate velocities from beam to instrument coordinates.

    :param dict adcpo: containing the beam velocity data.
    :param bool reverse: If True, this function performs the inverse rotation (inst->beam).
    :param bool force: When true do not check which coordinate system the data is in prior to performing this rotation.
    """
    if not force:
        if not reverse and adcpo['props']['coord_sys'] != 'beam':
            raise ValueError('The input must be in beam coordinates.')
        if reverse and adcpo['props']['coord_sys'] != 'inst':
            raise ValueError('The input must be in inst coordinates.')
    if 'rotmat' in adcpo['config'].keys():  
        rotmat = adcpo['config']['rotmat']
    else:
        rotmat = calc_beam_rotmatrix(adcpo['config']['beam_angle'],
                                     adcpo['config']['beam_pattern'].lower()
                                     == 'convex')
    cs = 'inst'
    if reverse:
        # Can't use transpose because rotation is not between
        # orthogonal coordinate systems
        rotmat = np.linalg.inv(rotmat)
        cs = 'beam'
    # raw = adcpo['vel'].transpose()
    raw = np.asmatrix(adcpo['vel'])
    # here I end up with an extra dimension of 4
    # vels = np.einsum('ij,jkl->ikl', rotmat, raw)
    # vels = np.einsum('ij,jk->ik', rotmat, raw)
    vels = np.array(np.asmatrix(rotmat)*raw.transpose())
    # vels = np.einsum('ij,jkl->ikl', rotmat, adcpo['vel'])
    # ValueError: operands could not be broadcast together with remapped
    # shapes [original->remapped]: (4,4)->(4,newaxis,newaxis,4) (16,4,1)->(4,1,16) 
    # adcpo['vel'] = np.einsum('ij,jkl->ikl', rotmat, adcpo['vel'])
    adcpo['vel'] = vels.transpose()
    adcpo['props']['coord_sys'] = cs


def inst2earth(adcpo, reverse=False, fixed_orientation=False, force=False):
    """
    Rotate velocities from the instrument to earth coordinates.

    :param dict adcpo: containing the data in instrument coordinates
    :param bool reverse: If True, this function performs the inverse rotation (earth->inst).
    :param bool fixed_orientation: When true, take the average orientation and apply it over the whole record.
    :param bool force: When true do not check which coordinate system the data is in prior to performing this rotation.

    Notes
    -----
    The rotation matrix is taken from the Teledyne RDI ADCP Coordinate Transformation manual January 2008

    When performing the forward rotation, this function sets the 'inst2earth:fixed' flag to the value of
    `fixed_orientation. When performing the reverse rotation, that value is 'popped' from the props dict and the input
    value to this function`fixed_orientation` has no effect. If `'inst2earth:fixed'` is not in the props dict then
    the input value *is* used.
    """
    deg2rad = np.pi / 180.
    if not force:
        if not reverse and adcpo['props']['coord_sys'] != 'inst':
            raise ValueError('The input must be in inst coordinates.')
        if reverse and adcpo['props']['coord_sys'] != 'earth':
            raise ValueError('The input must be in earth coordinates.')
    if not reverse and 'declination' in adcpo['props'].keys() and not adcpo['props']['declination_in_heading']:
        # Only do this if making the forward rotation.
        adcpo['heading_deg'] += adcpo['props']['declination']
        adcpo['props']['declination_in_heading'] = True
    r = adcpo['roll_deg'] * deg2rad
    p = np.arctan(np.tan(adcpo['pitch_deg'] * deg2rad) * np.cos(r))
    h = adcpo['heading_deg'] * deg2rad
    if adcpo['config']['orientation'].lower() == 'up': 
        r += np.pi
    ch = np.cos(h)
    sh = np.sin(h)
    cr = np.cos(r)
    sr = np.sin(r)
    cp = np.cos(p)
    sp = np.sin(p)
    rotmat = np.empty((3, 3, len(r)))
    rotmat[0, 0, :] = ch * cr + sh * sp * sr
    rotmat[0, 1, :] = sh * cp
    rotmat[0, 2, :] = ch * sr - sh * sp * cr
    rotmat[1, 0, :] = -sh * cr + ch * sp * sr
    rotmat[1, 1, :] = ch * cp
    rotmat[1, 2, :] = -sh * sr - ch * sp * cr
    rotmat[2, 0, :] = -cp * sr
    rotmat[2, 1, :] = sp
    rotmat[2, 2, :] = cp * cr
    # Only operate on the first 3-components, b/c the 4th is err_vel
    # ess = 'ijk,jlk->ilk'
    cs = 'earth'
    if reverse:
        cs = 'inst'
        fixed_orientation = adcpo['props'].pop('inst2earth:fixed', fixed_orientation)
        # ess = ess.replace('ij', 'ji')
    else:
        adcpo['props']['inst2earth:fixed'] = fixed_orientation
    if fixed_orientation:
        # ess = ess.replace('k,', ',')
        rotmat = rotmat.mean(-1)
    # todo is the einsum method better?  If so, uncomment the ess statements above
    # vels = np.einsum(ess, rotmat, adcpo['vel'][:,:3])
    vels = np.asmatrix(rotmat) * np.asmatrix(adcpo['vel'][:, :3].transpose())
    adcpo['vel'][:, :3] = vels.transpose()
    adcpo['props']['coord_sys'] = cs

# end code from dolfyn


def setupEPICnc(fname, rawcdf, attfile, settings):
    """
    Construct an empty netCDF output file to EPIC conventions

    :param str fname: output netCDF file name
    :param Dataset rawcdf: input netCDF raw data file object
    :param str attfile: metadata text file
    :param dict settings: settings as follows::

        'good_ensembles': [] # starting and ending indices of the input file. For all data use [0,np.inf]
        'orientation': 'UP' # uplooking ADCP, for downlooking, use DOWN
        'transducer_offset_from_bottom': 2.02 # in meters
        'transformation': 'EARTH' # | BEAM | INST
        'adjust_to_UTC': 5 # for EST to UTC, if no adjustment, set to 0 or omit

    :return: netCDF file object
    """
    # note that 
    # f4 = 4 byte, 32 bit float
    # maximum value for 32 bit float = 3.402823*10**38;
    intfill = -32768
    floatfill = 1E35
    
    # check the ensemble limits asked for by the user
    nens = rawcdf.variables['Rec'].size
    if settings['good_ensembles'][1] < 0:
        settings['good_ensembles'][1] = nens
    if settings['good_ensembles'][0] < 0:
        settings['good_ensembles'][0] = 0
    if settings['good_ensembles'][1] > nens:
        settings['good_ensembles'][1] = nens-1
    nens2write = settings['good_ensembles'][1]-settings['good_ensembles'][0]
                 
    print('creating netCDF file %s with %d records' % (fname, nens2write))
    
    rawvars = []
    for key in rawcdf.variables.keys():
        rawvars.append(key)
    
    nbins = len(rawcdf.dimensions['depth'])
    
    cdf = Dataset(fname, "w", clobber=True, format="NETCDF4")
    
    # dimensions, in EPIC order
    cdf.createDimension('time', nens2write)
    cdf.createDimension('depth', nbins)
    cdf.createDimension('lat', 1)
    cdf.createDimension('lon', 1)
    
    # write global attributes
    cdf.history = rawcdf.history + "rotations calculated and converted to EPIC format by ADCPcdf2ncEPIC.py"
    
    # these get returned as a dictionary
    gatts = read_globalatts(attfile)
    
    if 'WATER_DEPTH' not in gatts.keys():
        # noinspection PyTypeChecker
        gatts['WATER_DEPTH'] = 0.0  # nothing from user
        print('No WATER_DEPTH found, check depths of bins and WATER_DEPTH!')
    gatts['orientation'] = settings['orientation'].upper()
    
    if 'serial_number' not in gatts.keys():
        gatts['serial_number'] = "unknown"

    if 'magnetic_variation' not in gatts.keys():
        # noinspection PyTypeChecker
        gatts['magnetic_variation_at_site'] = 0.0
        print('No magnetic_variation, assuming magnetic_variation_at_site = 0')
    else:
        gatts['magnetic_variation_at_site'] = gatts['magnetic_variation']
        gatts.pop('magnetic_variation')
    if type(gatts['MOORING']) != str:
        gatts['MOORING'] = str(int(np.floor(gatts['MOORING'])))

    writeDict2atts(cdf, gatts, "")
    
    # more standard attributes
    cdf.latitude_units = "degree_north"
    cdf.longitude_units = "degree_east"
    cdf.CREATION_DATE = "%s" % datetime.now()
    cdf.DATA_TYPE = "ADCP"
    cdf.FILL_FLAG = 0
    cdf.COMPOSITE = 0
    
    # attributes that the names will vary depending on the ADCP vendor
    if rawcdf.sensor_type == "TRDI":
        # TRDI attributes
        if any('VBeam' in item for item in rawcdf.ncattrs()):
            cdf.INST_TYPE = "TRDI Workhorse V"
        else:
            cdf.INST_TYPE = "TRDI Workhorse"

        cdf.bin_size = rawcdf.TRDI_Depth_Cell_Length_cm/100
        cdf.bin_count = rawcdf.TRDI_Number_of_Cells
        cdf.center_first_bin = rawcdf.TRDI_Bin_1_distance_cm/100
        cdf.blanking_distance = rawcdf.TRDI_Blank_after_Transmit_cm/100
        cdf.transform = rawcdf.TRDI_Coord_Transform
        cdf.beam_angle = rawcdf.TRDI_Beam_Angle
        cdf.number_of_slant_beams = rawcdf.TRDI_Number_of_Beams
        cdf.heading_bias_applied_EB = rawcdf.TRDI_Heading_Bias_Hundredths_of_Deg
        cdf.beam_angle = rawcdf.TRDI_Beam_Angle
        cdf.beam_pattern = rawcdf.TRDI_Beam_Pattern
    elif rawcdf.sensor_type == "Nortek":
        # Nortek attributes
        # TODO - what to do about multiple sampling schemes?
        cdf.INST_TYPE = "Nortek Signature"
        cdf.bin_size = rawcdf.Nortek_burst_cellSize
        cdf.bin_count = rawcdf.Nortek_burst_nCells
        # Nortek Signature does not seem to have an odd offset to center of
        # bin 1.  This value comes from the Velocity Range provided by Nortek
        cdf.center_first_bin = rawcdf['bindist'][0]
        cdf.blanking_distance = rawcdf.Nortek_burst_blankingDistance
        cdf.transform = rawcdf.Nortek_burst_coordSystem
        # Nortek provides two angles for each beam, theta being from the vertical
        cdf.beam_angle = rawcdf.Nortek_beamConfiguration1_theta
        cdf.number_of_slant_beams = 4
        # there's no indication from Nortek's metadata that magvar is applied or not
        # TODO have to include information from the user
        cdf.heading_bias_applied_EB = 0
        # hard coded based on known Signature design
        # Could be deduced from the theta and phi beam angles
        cdf.beam_pattern = "Convex"

    # attributes requiring user input
    cdf.transducer_offset_from_bottom = settings['transducer_offset_from_bottom']
    cdf.initial_instrument_height = settings['transducer_offset_from_bottom']
    # TODO check on orientation, using user input for now
    # rawcdf.TRDI_Orientation 
    # need translation to UP from "Up-facing beams"
    cdf.orientation = settings['orientation'].upper()
    cdf.orientation_note = settings['orientation_note']
    if settings['orientation'].upper() == 'UP':
        cdf.depth_note = "uplooking bin depths = WATER_DEPTH-transducer_offset_from_bottom-bindist"
    else:
        cdf.depth_note = "downlooking bin depths = WATER_DEPTH-transducer_offset_from_bottom+bindist"
    cdf.serial_number = rawcdf.serial_number
    
    # TODO consider using a float for time since less common integers are causing issues
    # the problem is, CF time is a count, so integer is appropriate
    timetype = 'u2'  # u4 may be causing downstream problems with NCO
    # u2 caused rollover problems when EPIC time was stored or read:
    # file is 1108sig001.nc
    # EPIC first time stamp = 08-Oct-5378 00:01:04
    # seconds since 1970-01-01T00:00:00 UTC
    # CF first time stamp = 25-Sep-2017 15:00:00
    # the bad EPIC time is because a u4 datatype in the cdf file
    # is being sent to a u2 datatype in hte nc file.  Changing u2 to u4, etc.
    # causes other problems
    # timetype = 'u4' # u4 causes problems downstream in catEPIC with fill values
    # for now choosing to live with the u2 problems
    varobj = cdf.createVariable('Rec', timetype, ('time',), fill_value=intfill)
    varobj.units = "count"
    varobj.long_name = "Ensemble Number"

    if settings['timetype'] == 'CF':
        # if f8, 64 bit is not used, time is clipped
        # TODO test this theory, because downstream 64 bit time is a problem
        # for ADCP fast sampled, single ping data, need millisecond resolution
        # for CF convention
        varobj = cdf.createVariable('time', 'f8', ('time',))
        varobj.units = rawcdf.variables['time'].units
        # for EPIC convention
        varobj = cdf.createVariable('EPIC_time', timetype, ('time',))
        varobj.units = "True Julian Day"
        varobj.epic_code = 624
        varobj.datum = "Time (UTC) in True Julian Days: 2440000 = 0000 h on May 23, 1968"
        varobj.NOTE = "Decimal Julian day [days] = time [days] + ( time2 [msec] / 86400000 [msec/day] )"    
        varobj = cdf.createVariable('EPIC_time2', timetype, ('time',))
        varobj.units = "msec since 0:00 GMT"
        varobj.epic_code = 624
        varobj.datum = "Time (UTC) in True Julian Days: 2440000 = 0000 h on May 23, 1968"
        varobj.NOTE = "Decimal Julian day [days] = time [days] + ( time2 [msec] / 86400000 [msec/day] )"    
    else:
        # if f8, 64 bit is not used, time is clipped
        # for ADCP fast sampled, single ping data, need millisecond resolution
        # for CF convention
        varobj = cdf.createVariable('cf_time', 'f8', ('time',))
        varobj.units = rawcdf.variables['cf_time'].units
        # for EPIC convention
        varobj = cdf.createVariable('time', timetype, ('time',))
        varobj.units = "True Julian Day"
        varobj.epic_code = 624
        varobj.datum = "Time (UTC) in True Julian Days: 2440000 = 0000 h on May 23, 1968"
        varobj.NOTE = "Decimal Julian day [days] = time [days] + ( time2 [msec] / 86400000 [msec/day] )"    
        varobj = cdf.createVariable('time2', timetype, ('time',))
        varobj.units = "msec since 0:00 GMT"
        varobj.epic_code = 624
        varobj.datum = "Time (UTC) in True Julian Days: 2440000 = 0000 h on May 23, 1968"
        varobj.NOTE = "Decimal Julian day [days] = time [days] + ( time2 [msec] / 86400000 [msec/day] )"    

    varobj = cdf.createVariable('depth', 'f4', ('depth',))
    varobj.units = "m"
    varobj.long_name = "DEPTH (M)"
    varobj.epic_code = 3
    varobj.center_first_bin = cdf.center_first_bin
    varobj.blanking_distance = cdf.blanking_distance
    varobj.bin_size = cdf.bin_size
    varobj.bin_count = nbins
    varobj.transducer_offset_from_bottom = cdf.transducer_offset_from_bottom
    
    varobj = cdf.createVariable('lat', 'f8', ('lat',))
    varobj.units = "degree_north"
    varobj.epic_code = 500
    # note name is one of the netcdf4 reserved attributes, use setncattr
    varobj.setncattr('name', "LAT")
    varobj.long_name = "LATITUDE"
    varobj.datum = "NAD83"
    varobj[:] = float(gatts['latitude'])    

    varobj = cdf.createVariable('lon', 'f8', ('lon',))
    varobj.units = "degree_east"
    varobj.epic_code = 502
    # note name is one of the netcdf4 reserved attributes, use setncattr
    varobj.setncattr('name', "LON")
    varobj.long_name = "LONGITUDE"
    varobj.datum = "NAD83"
    varobj[:] = float(gatts['longitude'])    
    
    varobj = cdf.createVariable('bindist', 'f4', ('depth',), fill_value=floatfill)
    # note name is one of the netcdf4 reserved attributes, use setncattr
    varobj.setncattr('name', "bindist")
    varobj.units = "m"
    varobj.long_name = "bin distance from instrument"
    varobj.epic_code = 0
    varobj.center_first_bin = cdf.center_first_bin
    varobj.blanking_distance = cdf.blanking_distance
    varobj.bin_size = cdf.bin_size
    varobj.bin_count = nbins
    varobj.transducer_offset_from_bottom = cdf.transducer_offset_from_bottom
    varobj.NOTE = "distance is along profile from instrument head to center of bin"
    
    varobj = cdf.createVariable('SV_80', 'f4', ('time', 'lat', 'lon'), fill_value=floatfill)
    varobj.units = "m s-1"
    varobj.epic_code = 80
    # note name is one of the netcdf4 reserved attributes, use setncattr
    varobj.setncattr('name', "SV")
    varobj.long_name = "SOUND VELOCITY (M/S)"
    
    varobj = cdf.createVariable('Hdg_1215', 'f4', ('time', 'lat', 'lon'), fill_value=floatfill)
    varobj.units = "degrees"
    # note name is one of the netcdf4 reserved attributes, use setncattr
    varobj.setncattr('name', "Hdg")
    varobj.long_name = "INST Heading"
    varobj.epic_code = 1215
    # varobj.heading_alignment = rawvarobj.heading_alignment
    # varobj.heading_bias = rawvarobj.heading_bias
    
    varobj = cdf.createVariable('Ptch_1216', 'f4', ('time', 'lat', 'lon'), fill_value=floatfill)
    varobj.units = "degrees"
    varobj.long_name = "INST Pitch"
    varobj.epic_code = 1216
    
    varobj = cdf.createVariable('Roll_1217', 'f4', ('time', 'lat', 'lon'), fill_value=floatfill)
    varobj.units = "degrees"
    varobj.long_name = "INST Roll"
    varobj.epic_code = 1217
    
    varobj = cdf.createVariable('Tx_1211', 'f4', ('time', 'lat', 'lon'), fill_value=floatfill)
    varobj.units = "C"
    # note name is one of the netcdf4 reserved attributes, use setncattr
    varobj.setncattr('name', "T")
    varobj.long_name = "instrument Transducer Temp."
    varobj.epic_code = 1211

    if 'Pressure' in rawvars:
        # rawvarobj = rawcdf.variables['Pressure']
        varobj = cdf.createVariable('P_1', 'f4', ('time', 'lat', 'lon'), fill_value=floatfill)
        varobj.units = "dbar"
        # note name is one of the netcdf4 reserved attributes, use setncattr
        varobj.setncattr('name', "P")
        varobj.long_name = "PRESSURE (DB)"
        varobj.epic_code = 1
        
    if 'PressVar' in rawvars:
        varobj = cdf.createVariable('SDP_850', 'f4', ('time', 'lat', 'lon'), fill_value=floatfill)
        varobj.setncattr('name', 'SDP')
        varobj.long_name = "STAND. DEV. (PRESS)"
        varobj.units = "mbar"
        varobj.epic_code = 850

    varobj = cdf.createVariable('cor', 'u2', ('time', 'depth', 'lat', 'lon'), fill_value=intfill)
    varobj.setncattr('name', 'cor')
    varobj.long_name = "Slant Beam Average Correlation (cor)"
    varobj.units = "counts"
    varobj.epic_code = 1202
    varobj.NOTE = "Calculated from the slant beams"

    if 'PGd4' in rawvars:
        varobj = cdf.createVariable('PGd_1203', 'u2', ('time', 'depth', 'lat', 'lon'), fill_value=intfill)
        varobj.setncattr('name', 'Pgd')
        varobj.long_name = "Percent Good Pings"
        varobj.units = "percent"
        varobj.epic_code = 1203
        varobj.NOTE = "Percentage of good 4-bem solutions (Field #4)"

    varobj = cdf.createVariable('AGC_1202', 'u2', ('time', 'depth', 'lat', 'lon'), fill_value=intfill)
    varobj.setncattr('name', 'AGC')
    varobj.long_name = "Average Echo Intensity (AGC)"
    varobj.units = "counts"
    varobj.epic_code = 1202
    varobj.NOTE = "Calculated from the slant beams"

    if 'cor5' in rawvars:
        varobj = cdf.createVariable('corvert', 'u2', ('time', 'depth', 'lat', 'lon'), fill_value=intfill)
        varobj.setncattr('name', 'cor')
        varobj.long_name = "Vertical Beam Correlation (cor)"
        varobj.units = "counts"
        varobj.epic_code = 1202
        varobj.NOTE = "From the center vertical beam"

    if 'att5' in rawvars:
        varobj = cdf.createVariable('AGCvert', 'u2', ('time', 'depth', 'lat', 'lon'), fill_value=intfill)
        varobj.setncattr('name', 'AGC')
        varobj.long_name = "Vertical Beam Echo Intensity (AGC)"
        varobj.units = "counts"
        varobj.epic_code = 1202
        varobj.NOTE = "From the center vertical beam"
        
    # repeating attributes that do not depend on height or depth calculations
    cdfvarnames = []
    for key in cdf.variables.keys():
        cdfvarnames.append(key)
    omitnames = []
    for key in cdf.dimensions.keys():
        omitnames.append(key)
    omitnames.append("Rec")
    for varname in cdfvarnames:
        if varname not in omitnames:
            varobj = cdf.variables[varname]
            varobj.serial_number = cdf.serial_number

    if settings['transformation'].upper() == "BEAM":
        varnames = ["Beam1", "Beam2", "Beam3", "Beam4"]
        codes = [0, 0, 0, 0]
    elif settings['transformation'].upper() == "INST":
        varnames = ["X", "Y", "Z", "Error"]
        codes = [0, 0, 0, 0]
    else:
        varnames = ["u_1205", "v_1206", "w_1204", "Werr_1201"]
        codes = [1205, 1206, 1204, 1201]
        
    for i in range(4):
        varname = varnames[i]
        varobj = cdf.createVariable(varname, 'f4', ('time', 'depth', 'lat', 'lon'), fill_value=floatfill)
        varobj.units = "cm s-1"
        varobj.long_name = "%s velocity (cm s-1)" % varnames[i]
        varobj.epic_code = codes[i]

    if 'vel5' in rawvars:
        varobj = cdf.createVariable('Wvert', 'f4', ('time', 'depth', 'lat', 'lon'), fill_value=floatfill)
        varobj.units = "cm s-1"
        varobj.long_name = "Vertical velocity (cm s-1)" 

    # TODO do we do bottom track data here?  Later?  Or as a separate thing?

    add_VAR_DESC(cdf)
        
    return cdf


# noinspection PyUnresolvedReferences
def add_VAR_DESC(cdf):
    """
    add the VAR_DESC global attribute constructed from variable names found in the file

    :param object cdf: netCDF file object
    """
    # cdf is an netcdf file object (e.g. pointer to open netcdf file)

    varkeys = cdf.variables.keys()  # get the names
    dimkeys = cdf.dimensions.keys()
    varnames = []
    for key in varkeys:
        varnames.append(key)
    dimnames = []
    for key in dimkeys:
        dimnames.append(key)
    buf = ""
    for varname in varnames:
        if varname not in dimnames:
            buf = "%s:%s" % (buf, varname)

    cdf.VAR_DESC = buf
    
    
def read_globalatts(fname):
    """
    read_globalatts: read in file of metadata for a tripod or mooring

    reads global attributes for an experiment from a text file (fname) called by all data processing programs
    to get uniform metadata input one argument is required- the name of the file to read- it should have this form::

        SciPi; J.Q. Scientist
        PROJECT; USGS Coastal Marine Geology Program
        EXPERIMENT; MVCO 2015 Stress Comparison
        DESCRIPTION; Quadpod 13.9m
        DATA_SUBTYPE; MOORED
        COORD_SYSTEM; GEOGRAPHIC + SAMPLE
        Conventions; PMEL/EPIC
        MOORING; 1057
        WATER_DEPTH; 13.9
        WATER_DEPTH_NOTE; (meters), nominal
        WATER_DEPTH_source; ship fathometer
        latitude; 41.3336633
        longitude; -70.565877
        magnetic_variation; -14.7
        Deployment_date; 17-Nov-2015
        Recovery_date; 14-Dec-2015
        DATA_CMNT;
        platform_type; USGS aluminum T14 quadpod
        DRIFTER; 0
        POS_CONST; 0
        DEPTH_CONST; 0
        Conventions; PMEL/EPIC
        institution; United States Geological Survey, Woods Hole Coastal and Marine Science Center
        institution_url; http://woodshole.er.usgs.gov

    :param str fname: input file name
    :return: dict of metadata
    """
    gatts = {}
    f = open(fname, 'r')
    for line in f:
        line = line.strip()
        cols = line.split(";")
        gatts[cols[0]] = cols[1].strip()

    f.close()
    return gatts


# noinspection PyUnresolvedReferences
def writeDict2atts(cdfobj, d, tag):
    """
    write a dictionary to netCDF attributes

    :param object cdfobj: netcdf file object
    :param dict d: metadata
    :param str tag: tag to add before each atrribute name
    :return: dict of metadata as written to file
    """
    i = 0
    # first, convert as many of the values in d to numbers as we can
    for key in iter(d):
        if type(d[key]) == str:
            try:
                d[key] = float(d[key])
            except ValueError:
                i += 1

    for key in iter(d):
        newkey = tag + key
        try:
            cdfobj.setncattr(newkey, d[key])
        except:
            print('can\'t set %s attribute' % key)
    
    return d


def floor(dec):
    """
    convenience function to round down
    provided to avoid loading the math package
    and because np.floor was causing unexpected behavior w.r.t ints

    :param float dec:
    :return: rounded number
    """
    return int(dec - (dec % 1))
    

def __main():
    print('%s running on python %s' % (sys.argv[0], sys.version))

    if len(sys.argv) < 2:
        print("%s useage:" % sys.argv[0])
        print("ADCPcdf2ncEPIC rawcdfname ncEPICname USGSattfile [startingensemble endingensemble]\n")
        print("starting and ending ensemble are netcdf file indeces, NOT TRDI ensemble numbers")
        print("USGSattfile is a file containing EPIC metadata")
        sys.exit(1)
        
    try:
        infile_name = sys.argv[1]
    except:
        print('error - input file name missing')
        sys.exit(1)
        
    try:
        outfile_name = sys.argv[2]
    except:
        print('error - output file name missing')
        sys.exit(1)
        
    try:
        attfile_name = sys.argv[3]
    except:
        print('error - global attribute file name missing')
        sys.exit(1)
        
    try:
        settings = sys.argv[4]
    except:
        print('error - settings missing - need dictionary of:')
        print('settings[\'good_ensembles\'] = [0, np.inf] # use np.inf for all ensembles or omit')
        print('settings[\'transducer_offset_from_bottom\'] = 1 # m')
        print('settings[\'transformation\'] = "EARTH" # | BEAM | INST')
        sys.exit(1)

    # some input testing
    if ~os.path.isfile(infile_name):
        print('error - input file not found')
        sys.exit(1)
        
    if ~os.path.isfile(attfile_name):
        print('error - attribute file not found')
        sys.exit(1)

    print('Converting %s to %s' % (infile_name, outfile_name))
    
    print('Start file conversion at ', dt.datetime.now())

    # noinspection PyTypeChecker
    doEPIC_ADCPfile(infile_name, outfile_name, attfile_name, settings)
    
    print(f'Finished file conversion at {dt.datetime.now()}')


if __name__ == "__main__":
    __main()
