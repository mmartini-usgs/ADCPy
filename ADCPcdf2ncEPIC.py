# -*- coding: utf-8 -*-
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
timetype = 'EPIC' for EPIC time and time2 convention, 
           'CF' for the time variable to be CF convention
           if CF is specified, EPIC time will be included as EPIC_time and EPIC_time2
settings = a dictionary of preferences for the processing
    settings['good_ensembles'] = starting and ending indeces of the input file.  
        For all data use [0,np.inf]
    settings['orientation'] = 'UP' # uplooking ADCP, for downlookers, use DOWN
    settings['transducer_offset_from_bottom'] = 2.02
    settings['transformation'] = 'EARTH' # | BEAM | INST
    # this code depnds on time being in UTC because of various time standards
    settings['adjust_to_UTC'] = 5 # for EST to UTC, if no adjustment, set to 0 or omit
    
Depth dependent attributes are compute from the mean Pressure found in the raw
data file.  So it is best to have the time series trimmed to the in water
time or to provide the good ensemble indeces for in water time

note that file names and paths may not include spaces

Example contents of a Global Attribute file:
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
#from math import floor
#import netCDF4 as nc4
from netCDF4 import Dataset
#from netCDF4 import num2date
import netCDF4 as netcdf
import datetime as dt
from datetime import datetime
#from TRDIpd0tonetcdf import julian
from TRDIpd0tonetcdf import ajd

def doEPIC_ADCPfile(cdfFile, ncFile, attFile, settings):
    
    # check some of the settings we can't live without
    # set flags, then remove fromt he settings list if we don't want them in 
    # metadata
    if 'good_ensembles' not in settings.keys():
        settings['good_ensembles'] = [0, np.inf] # nothing from user, do them all
        print('No starting and ending ensembles specfied, processing entire file')
    if 'orientation' not in settings.keys():
        settings['orientation'] = "UP"
        settings['orientation_note'] = "assumed by program"
        print('No orientation specfied, assuming up-looking')
    else: settings['orientation_note'] = "user provided orientation"        
    if 'transducer_offset_from_bottom' not in settings.keys():
        settings['transducer_offset_from_bottom'] = 0
        print('No transducer_offset_from_bottom, assuming 0')
    if 'transformation' not in settings.keys():
        settings['transformation'] = "EARTH"
    if 'timetype' not in settings.keys():
        settings['timetype'] = "CF"
    if 'adjust_to_UTC' not in settings.keys():
        settings['adjust_to_UTC'] = 0
    if 'beam_velocity_multiplier' in settings.keys():
        beamvelmultiplier = settings['beam_velocity_multiplier']
        settings.pop('beam_velocity_multiplier') 
    else: beamvelmultiplier = 1
    if 'use_pressure_for_WATER_DEPTH' in settings.keys():
        if settings['use_pressure_for_WATER_DEPTH']:
            usep4waterdepth = True
            settings.pop('use_pressure_for_WATER_DEPTH')
        else: usep4waterdepth = False
        
    rawcdf = Dataset(cdfFile, mode='r',format='NETCDF4')
    rawvars = []
    for key in rawcdf.variables.keys(): rawvars.append(key)
      
    # this function will operate on the files using the netCDF package
    nc = setupEPICnc(ncFile, rawcdf, attFile, settings)
    
    nbeams = nc.number_of_slant_beams #what if this isn't 4?
    nbins = len(rawcdf.dimensions['depth'])
    nens = len(rawcdf.dimensions['time'])
    ncvars = []
    for key in nc.variables.keys(): ncvars.append(key)

    declination = nc.magnetic_variation_at_site

    # start and end indices
    s = settings['good_ensembles'][0]
    if settings['good_ensembles'][1] < 0:
        e = nens
    else:
        e = settings['good_ensembles'][1]
    print('Converting from index %d to %d of %s' % (s,e,cdfFile))

    # many variables do not need processing and can just be copied to the
    # new EPIC convention
    varlist = {'sv':'SV_80','Rec':'Rec'}
    
    for key in varlist:
        varobj = nc.variables[varlist[key]]
        varobj[:] = rawcdf.variables[key][s:e]  
        
    # check the time zone, Nortek data are usually set to UTC, no matter what
    # the actual time zone of deployment might have been
    if abs(settings['adjust_to_UTC']) > 0:
        nc.time_zone_change_applied = settings['adjust_to_UTC']
        nc.time_zone_change_applied_note = "adjust time to UTC requested by user"    
    toffset = settings['adjust_to_UTC']*3600    
    # raw variable name : EPIC variable name
    # EPIC time (time, time2) is the default time convention
    if settings['timetype'] == 'CF':
        varlist = {'time':'time','EPIC_time':'EPIC_time','EPIC_time2':'EPIC_time2','sv':'SV_80','Rec':'Rec'}
    else:
        varlist = {'cf_time':'cf_time','time':'time','time2':'time2','sv':'SV_80','Rec':'Rec'}
    
    for key in varlist:
        varobj = nc.variables[varlist[key]]
        varobj[:] = rawcdf.variables[key][s:e]+toffset 
        
    # TRDI instruments have heading, pitch, roll and temperature in hundredths of degrees
    if rawcdf.sensor_type == "TRDI":
        degree_factor = 100
    else:
        degree_factor = 1
    
    varlist = {'Ptch':'Ptch_1216','Roll':'Roll_1217',
               'Tx':'Tx_1211'}
    
    for key in varlist:
        varobj = nc.variables[varlist[key]]
        varobj[:] = rawcdf.variables[key][s:e]/degree_factor 
              
    # TODO will need an isntrument dependent methodology to check for any
    # previous adjustments to heading prior to this correction.
    # for instance, with TRDI instruments, Velocity or the EB command
    # might have applied a correction.  If EB is set, then that value was 
    # applied to the raw data seen by TRDIpd0tonetcdf.py
    nc.magnetic_variation_applied = declination
    nc.magnetic_variation_applied_note = "as provided by user"
    heading = rawcdf.variables['Hdg'][s:e]/degree_factor + declination
    heading[heading >= 360] = heading[heading >= 360] - 360
    heading[heading < 0] = heading[heading < 0] + 360
    nc['Hdg_1215'][:] = heading
    
    # pressure needs to be in db or m
    if 'Pressure' in rawvars:
        pconvconst = 1 # when in doubt, do nothing
        punits = rawcdf['Pressure'].units
        if 'deca-pascals' in punits:
            pconvconst = 1000 # decapascals to dbar = /1000
            print('Pressure in deca-pascals will be converted to db' )
        nc['P_1'][:] = rawcdf.variables['Pressure'][s:e]/pconvconst
    
    # check units of current velocity and convert to cm/s
    vunits = rawcdf['vel1'].units
    vconvconst = 1 # when in doubt, do nothing
    if ('mm s-1' in vunits) | ('mm/s' in vunits):
        vconvconst = 0.1 # mm/s to cm/s
        print('Velocity in mm s-1 will be converted to cm s-1')
    if ('m s-1' in vunits) | ('m/s' in vunits):
        vconvconst = 100 # m/s to cm/s
        print('Velocity in m s-1 will be converted to cm s-1')

    if 'vel5' in rawvars:
        nc['Wvert'][:] = rawcdf.variables['vel5'][s:e,:] * vconvconst 

    if 'cor5' in rawvars:
        nc['corvert'][:] = rawcdf.variables['cor5'][s:e,:]
    
    if 'att5' in rawvars:
        nc['AGCvert'][:] = rawcdf.variables['att5'][s:e,:]

    if 'PGd4' in rawvars:
        nc['PGd_1203'][:,:,0,0] = rawcdf.variables['PGd4'][s:e,:]
        
    if 'PressVar' in rawvars:
        nc['SDP_850'][:,0,0] = rawcdf.variables['PressVar'][s:e]
    
    varobj = nc.variables['bindist']
    bindist = np.arange(len(nc['bindist']))
    bindist = bindist*nc.bin_size+nc.center_first_bin
    nc['bindist'][:] = bindist

    # figure out DELTA_T - we need to use the cf_time, more convenient
    if settings['timetype'] == 'CF':
        # we will set the main "time' variable to CF convention
        timekey = 'time'
    else:
        timekey = 'cf_time'
    
    # this calculation uses for CF time
    dtime = np.diff( nc[timekey][:])
    DELTA_T = '%s' % int((dtime.mean().astype('float')).round())
    nc.DELTA_T = DELTA_T
    
    # depths and heights
    nc.initial_instrument_height = settings['transducer_offset_from_bottom']
    nc.initial_instrument_height_note = "height in meters above bottom: accurate for tripod mounted instruments" 
    # compute depth, make a guess we want to average all depths recorded 
    # deeper than user supplied water depth
    # idx is returned as a tuple, the first of which is the actual index values

    # set the water depth here, this will be used throughout
    if ('Pressure' in rawvars) and usep4waterdepth:
        idx = np.where(nc['P_1'] > nc.WATER_DEPTH/2)
        # now for the mean of only on bottom pressure measurements
        pmean = nc['P_1'][idx[0]].mean()
        print('Site WATER_DEPTH given is %f' % nc.WATER_DEPTH)
        print('Calculated mean water level from P_1 is %f m' % pmean)
        print('Updating site WATER_DEPTH to %f m' % pmean)
        nc.WATER_DEPTH = pmean+nc.transducer_offset_from_bottom
        nc.WATER_DEPTH_source = "water depth = MSL from pressure sensor, (meters), nominal"
        nc.WATER_DEPTH_NOTE = nc.WATER_DEPTH_source
        nc.nominal_sensor_depth = nc.WATER_DEPTH-settings['transducer_offset_from_bottom']
        nc.nominal_sensor_depth_note = "inst_depth = (water_depth - inst_height); nominal depth below surface, meters"
        varnames = ['P_1','bindist','depth']
        # WATER_DEPTH_datum is not used in this circumstance.
    else:
        print('Site WATER_DEPTH given is %f' % nc.WATER_DEPTH)
        print('No pressure data available, so no adjustment to water depth made')
        nc.WATER_DEPTH_source = "water depth as given by user, (meters), nominal"
        nc.WATER_DEPTH_NOTE = nc.WATER_DEPTH_source
        nc.nominal_sensor_depth = nc.WATER_DEPTH-settings['transducer_offset_from_bottom']
        nc.nominal_sensor_depth_note = "inst_depth = (water_depth - inst_height); nominal depth below surface, meters"
        varnames = ['bindist','depth']
        # WATER_DEPTH_datum is not used in this circumstance.

    for varname in varnames:
        nc[varname].WATER_DEPTH = nc.WATER_DEPTH
        nc[varname].WATER_DEPTH_source = nc.WATER_DEPTH_source
        nc[varname].transducer_offset_from_bottom = nc.transducer_offset_from_bottom
    
    # update depth variable for location of bins based on WATER_DEPTH information
    if "UP" in nc.orientation:
        depths = nc.WATER_DEPTH-nc.transducer_offset_from_bottom-nc['bindist']
    else:
        depths = -1 * (nc.WATER_DEPTH-nc.transducer_offset_from_bottom+nc['bindist'])
    
    nc['depth'][:] = depths        
    
    nc.start_time = '%s' % netcdf.num2date(nc[timekey][0],nc[timekey].units)
    nc.stop_time = '%s' % netcdf.num2date(nc[timekey][-1],nc[timekey].units)
    
    # some of these repeating attributes depended on depth calculations
    # these are the same for all variables because all sensors are in the same
    # package, as of now, no remote sensors being logged by this ADCP
    ncvarnames = []
    for key in nc.variables.keys(): ncvarnames.append(key)
    omitnames = []
    for key in nc.dimensions.keys(): omitnames.append(key)
    omitnames.append("Rec")
    omitnames.append("depth")
    for varname in ncvarnames:
        if varname not in omitnames:
            varobj = nc.variables[varname]
            varobj.sensor_type = nc.INST_TYPE
            varobj.sensor_depth = nc.nominal_sensor_depth
            varobj.initial_sensor_height = nc.initial_instrument_height 
            varobj.initial_sensor_height_note = "height in meters above bottom:  accurate for tripod mounted instruments"
            varobj.height_depth_units = "m"

    print('finished copying data, starting computations at %s' % (dt.datetime.now()))

    print('averaging cor at %s' % (dt.datetime.now()))
    # this will be a problem - it loads all into memory
    cor = (rawcdf.variables['cor1'][s:e,:]+rawcdf.variables['cor2'][s:e,:]+ \
        rawcdf.variables['cor3'][s:e,:]+rawcdf.variables['cor4'][s:e,:]) / 4
    #varobj = nc.variables['cor']
    nc['cor'][:,:,0,0] = cor[:,:]

    print('averaging AGC at %s' % (dt.datetime.now()))
    # this will be a problem - it loads all into memory
    agc = (rawcdf.variables['att1'][s:e,:]+rawcdf.variables['att2'][s:e,:]+ \
        rawcdf.variables['att3'][s:e,:]+rawcdf.variables['att4'][s:e,:]) / 4
    #varobj = nc.variables['cor']
    nc['AGC_1202'][:,:,0,0] = agc[:,:]
    
    print('converting %d ensembles from beam to earth %s' % (len(nc[timekey]), dt.datetime.now()))
    
    # check our indexing
    print('magnetic variation at site = %f' % nc.magnetic_variation_at_site)
    print('magnetic variation applied = %f' % nc.magnetic_variation_applied)
    print('magnetic variation applied note = %s' % nc.magnetic_variation_applied_note)
    n = int(len(heading)/2)
    print('From the middle of the time series at ensemble #%d, we have:' % n)
    print('heading variable in this python process = %f' % heading[n])
    print('rawcdf Hdg[n] = %f' % rawcdf['Hdg'][n])
    print('nc Hdg_1215[n] = %f' % nc['Hdg_1215'][n,0,0])
    
    # TODO add depth bin mapping
    
    # this beam arrangement is for TRDI Workhorse and V, other instruments 
    # should be re-ordered to match
    rawvarnames = ["vel1", "vel2", "vel3", "vel4"]
    ncidx = 0
    if settings['transformation'].upper() == "BEAM":
        ncvarnames = ["Beam1", "Beam2", "Beam3", "Beam4"]        
        for idx in range(s,e):
            for beam in range(nbeams):
                nc[ncvarnames[beam]][ncidx,:,0,0] = \
                    rawcdf.variables[rawvarnames[beam]][idx,:] * vconvconst
            ncidx = ncidx + 1
            
    elif (settings['transformation'].upper() == "INST") or \
        (settings['transformation'].upper() == "EARTH"):
        ncvarnames = ["X","Y","Z","Error"]            
        # the dolfyn way
        # load the ADCP data object - we have converted this from a class object
        # to nested dictionaries for use here
        adcpo = {}
        adcpo['props'] = {}
        adcpo['props']['coord_sys'] = "beam"
        adcpo['props']['inst2earth:fixed'] = False
        adcpo['config'] = {}
        adcpo['config']['beam_angle'] = nc.beam_angle
        adcpo['config']['beam_pattern'] = nc.beam_pattern
        adcpo['config']['orientation'] = nc.orientation
        adcpo['declination_in_heading'] = True
        # dolfyn shape for ensemble data is [bins x beams x ens]
        adcpo['vel'] = np.ones([nbins,nbeams],dtype='float')*np.NAN
        # vels has to be pre-defined to get the shapes to broadcast
        vels = np.ones([nbins,1],dtype='float')*np.NAN
        
        # Nortek and TRDI do their along beam velocity directions opposite for
        # slant beams.  Vertical beam directions are the same.
        if rawcdf.sensor_type == 'Nortek':
            beamvelmultiplier = -1
        else:
            beamvelmultiplier = 1
        
        for idx in range(s,e):
            for beam in range(nbeams):
                # load data of one ensemble to dolfyn shape, in cm/s
                #adcpo['vel'][:,beam,0] = rawcdf.variables[rawvarnames[beam]][idx,:] * 0.1
                vels = rawcdf.variables[rawvarnames[beam]][idx,:] \
                    * vconvconst * beamvelmultiplier
                adcpo['vel'][:,beam] = vels
            
            # need to keep setting this with new beam data since we are iterating
            adcpo['props']['coord_sys'] = "beam" 
            beam2inst(adcpo) # adcpo['vel'] is returned in inst coordinates
        
            if settings['transformation'].upper() == "EARTH":
                ncvarnames = ["u_1205","v_1206","w_1204","Werr_1201"]
                adcpo['heading_deg'] = nc.variables['Hdg_1215'][ncidx]
                adcpo['pitch_deg'] = nc.variables['Ptch_1216'][ncidx]
                adcpo['roll_deg']  = nc.variables['Roll_1217'][ncidx]
                inst2earth(adcpo)
                
            for beam in range(nbeams):
                nc[ncvarnames[beam]][ncidx,:,0,0] = adcpo['vel'][:,beam] 
                   
            ncidx = ncidx + 1
            
            # immediate - then less feedback
            ensf, ensi = np.modf(ncidx/1000)
            if (ensf == 0) and (ncidx < 10000):
                print('%d of %d ensembles read' % (ncidx, nens))
            else:
                ensf, ensi = np.modf(ncidx/10000)
                if ensf == 0:
                    print('%d of %d ensembles read' % (ncidx, nens))
            
    # minima and maxima to be added as a separate operation after averaging

    print('closing files at %s' % (dt.datetime.now()))

    rawcdf.close()
    nc.close()            

   
""" this transformation matrix is from the R.D. Instruments Coordinate 
    Transformation booklet.  It presumes the beams are in the same position as
    RDI Workhorse ADCP beams, where, when looking down on the transducers:
        Beam 3 is in the direction of the compass' zero reference
        Beam 1 is to the right
        Beam 2 is to the left
        Beam 4 is opposite beam 3
        Pitch is about the beam 2-1 axis and is positive when beam 3 is raised
        Roll is about the beam 3-4 axis and is positive when beam 2 is raised
        Heading increases when beam 3 is rotated towards beam 1
    Nortek Signature differs in these ways:
        TRDI beam 3 = Nortek beam 1
        TRDI beam 1 = Nortek beam 2
        TRDI beam 4 = Nortek beam 3
        TRDI beam 2 = Nortek beam 4
        Heading, pitch and roll behave the same as TRDI
"""
# as per the R.D. Instruments Coordinate Transformation booklet                     
def cal_earth_rotmatrix(heading=0,pitch=0,roll=0,declination=0):
    # for declination, West is negative
    # heading, pitch and roll are in degrees
    ch = np.cos(heading)
    sh = np.sin(heading)
    cp = np.cos(pitch)
    sp = np.sin(pitch)
    cr = np.cos(roll)
    sr = np.sin(roll)
    
    return np.asmatrix(np.array([[(ch*cr+sh*sp*sr), (sh*cp), (ch*sr-sh*sp*cr)],
                                [(-sh*cr+ch*sp*sr), (ch*cp), (-sh*sr-ch*sp*cr)],
                                [(-cp*sr), sp, (cp*cr)]]))

# as per dolfyn rotate.py code here: https://github.com/lkilcher/dolfyn

def calc_beam_rotmatrix(theta=20, convex=True, degrees=True):
    """Calculate the rotation matrix from beam coordinates to
    instrument head coordinates.

    Parameters
    ----------
    theta : is the angle of the heads (usually 20 or 30 degrees)

    convex : is a flag for convex or concave head configuration.

    degrees : is a flag which specifies whether theta is in degrees
        or radians (default: degrees=True)
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
    tmp = []
    for vl in tpl:
        tmp.append(vl[:, None, :])
    return np.concatenate(tuple(tmp), axis=1)

def beam2inst(adcpo, reverse=False, force=False):
    """Rotate velocitiesfrom beam to instrument coordinates.

    Parameters
    ----------
    adpo : The ADP object containing the data.

    reverse : bool (default: False)
           If True, this function performs the inverse rotation
           (inst->beam).
    force : bool (default: False)
        When true do not check which coordinate system the data is in
        prior to performing this rotation.
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
    #raw = adcpo['vel'].transpose()
    raw = np.asmatrix(adcpo['vel'])
    # here I end up with an extra dimension of 4
    #vels = np.einsum('ij,jkl->ikl', rotmat, raw)
    #vels = np.einsum('ij,jk->ik', rotmat, raw)
    vels = np.array(np.asmatrix(rotmat)*raw.transpose())
    #vels = np.einsum('ij,jkl->ikl', rotmat, adcpo['vel'])
    #ValueError: operands could not be broadcast together with remapped 
    # shapes [original->remapped]: (4,4)->(4,newaxis,newaxis,4) (16,4,1)->(4,1,16) 
    #adcpo['vel'] = np.einsum('ij,jkl->ikl', rotmat, adcpo['vel'])
    adcpo['vel'] = vels.transpose()
    adcpo['props']['coord_sys'] = cs

def inst2earth(adcpo, reverse=False, fixed_orientation=False, force=False):
    """Rotate velocities from the instrument to earth coordinates.

    Parameters
    ----------
    adpo : The ADP object containing the data.

    reverse : bool (default: False)
           If True, this function performs the inverse rotation
           (earth->inst).
    fixed_orientation : bool (default: False)
        When true, take the average orientation and apply it over the
        whole record.
    force : bool (default: False)
        When true do not check which coordinate system the data is in
        prior to performing this rotation.

    Notes
    -----
    The rotation matrix is taken from the Teledyne RDI ADCP Coordinate
    Transformation manual January 2008

    When performing the forward rotation, this function sets the
    'inst2earth:fixed' flag to the value of `fixed_orientation`. When
    performing the reverse rotation, that value is 'popped' from the
    props dict and the input value to this function
    `fixed_orientation` has no effect. If `'inst2earth:fixed'` is not
    in the props dict than the input value *is* used.
    """
    deg2rad = np.pi / 180.
    if not force:
        if not reverse and adcpo['props']['coord_sys'] != 'inst':
            raise ValueError('The input must be in inst coordinates.')
        if reverse and adcpo['props']['coord_sys'] != 'earth':
            raise ValueError('The input must be in earth coordinates.')
    if (not reverse and 'declination' in adcpo['props'].keys() and not
            adcpo['props']['declination_in_heading']):
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
    ess = 'ijk,jlk->ilk'
    cs = 'earth'
    if reverse:
        cs = 'inst'
        fixed_orientation = adcpo['props'].pop('inst2earth:fixed', fixed_orientation)
        ess = ess.replace('ij', 'ji')
    else:
        adcpo['props']['inst2earth:fixed'] = fixed_orientation
    if fixed_orientation:
        ess = ess.replace('k,', ',')
        rotmat = rotmat.mean(-1)
    #vels = np.einsum(ess, rotmat, adcpo['vel'][:,:3])
    vels = np.asmatrix(rotmat) * np.asmatrix(adcpo['vel'][:,:3].transpose())
    adcpo['vel'][:,:3] = vels.transpose()
    adcpo['props']['coord_sys'] = cs

def inst2earth_heading(adpo):
    deg2rad = np.pi / 180.
    h = adpo.heading_deg[:] * deg2rad
    if 'heading_offset' in adpo.props.keys():
        h += adpo['props']['heading_offset'] * deg2rad
    if 'declination' in adpo['props'].keys():
        h += adpo['props']['declination'] * deg2rad
    return np.exp(-1j * h)
    
# end code from dolfyn

def setupEPICnc(fname, rawcdf, attfile, settings):
     
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
    for key in rawcdf.variables.keys(): rawvars.append(key)
    
    nbins = len(rawcdf.dimensions['depth'])
    
    cdf = Dataset(fname, "w", clobber=True, format="NETCDF4")
    
    # dimensions, in EPIC order
    cdf.createDimension('time',nens2write)
    cdf.createDimension('depth',nbins)
    cdf.createDimension('lat',1)
    cdf.createDimension('lon',1)
    
    # write global attributes
    cdf.history = rawcdf.history+"rotations calculated and converted to EPIC format by ADCPcdf2ncEPIC.py"
    
    # these get returned as a dictionary
    gatts = read_globalatts(attfile)
    
    if 'WATER_DEPTH' not in gatts.keys():
        gatts['WATER_DEPTH'] = 0 # nothing from user
        print('No WATER_DEPTH found, check depths of bins and WATER_DEPTH!')
    gatts['orientation'] = settings['orientation'].upper()
    
    if 'serial_number' not in gatts.keys():
        gatts['serial_number'] = "unknown"
    
    if 'magnetic_variation' not in gatts.keys():
        gatts['magnetic_variation_at_site'] = 0
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
        else: cdf.INST_TYPE = "TRDI Workhorse"
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
        cdf.beam_angle = rawcdf.Nortek_BeamCfg1_theta
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
    
      
    varobj = cdf.createVariable('Rec','u4',('time'),fill_value=intfill)
    varobj.units = "count"
    varobj.long_name = "Ensemble Number"

    if settings['timetype'] == 'CF':
        # if f8, 64 bit is not used, time is clipped
        # for ADCP fast sampled, single ping data, need millisecond resolution
        # for CF convention
        varobj = cdf.createVariable('time','f8',('time'))
        varobj.units = rawcdf.variables['time'].units
        # for EPIC convention
        varobj = cdf.createVariable('EPIC_time','u4',('time'))
        varobj.units = "True Julian Day"
        varobj.epic_code = 624
        varobj.datum = "Time (UTC) in True Julian Days: 2440000 = 0000 h on May 23, 1968"
        varobj.NOTE = "Decimal Julian day [days] = time [days] + ( time2 [msec] / 86400000 [msec/day] )"    
        varobj = cdf.createVariable('EPIC_time2','u4',('time'))
        varobj.units = "msec since 0:00 GMT"
        varobj.epic_code = 624
        varobj.datum = "Time (UTC) in True Julian Days: 2440000 = 0000 h on May 23, 1968"
        varobj.NOTE = "Decimal Julian day [days] = time [days] + ( time2 [msec] / 86400000 [msec/day] )"    
    else:
        # if f8, 64 bit is not used, time is clipped
        # for ADCP fast sampled, single ping data, need millisecond resolution
        # for CF convention
        varobj = cdf.createVariable('cf_time','f8',('time'))
        varobj.units = rawcdf.variables['cf_time'].units
        # for EPIC convention
        varobj = cdf.createVariable('time','u4',('time'))
        varobj.units = "True Julian Day"
        varobj.epic_code = 624
        varobj.datum = "Time (UTC) in True Julian Days: 2440000 = 0000 h on May 23, 1968"
        varobj.NOTE = "Decimal Julian day [days] = time [days] + ( time2 [msec] / 86400000 [msec/day] )"    
        varobj = cdf.createVariable('time2','u4',('time'))
        varobj.units = "msec since 0:00 GMT"
        varobj.epic_code = 624
        varobj.datum = "Time (UTC) in True Julian Days: 2440000 = 0000 h on May 23, 1968"
        varobj.NOTE = "Decimal Julian day [days] = time [days] + ( time2 [msec] / 86400000 [msec/day] )"    

    varobj = cdf.createVariable('depth','f4',('depth'))
    varobj.units = "m"
    varobj.long_name = "DEPTH (M)"
    varobj.epic_code = 3
    varobj.center_first_bin = cdf.center_first_bin
    varobj.blanking_distance = cdf.blanking_distance
    varobj.bin_size = cdf.bin_size
    varobj.bin_count = nbins
    varobj.transducer_offset_from_bottom = cdf.transducer_offset_from_bottom
    
    varobj = cdf.createVariable('lat','f8',('lat'))
    varobj.units = "degree_north"
    varobj.epic_code = 500
    # note name is one of the netcdf4 reserved attributes, use setncattr
    varobj.setncattr('name', "LAT")
    varobj.long_name = "LATITUDE"
    varobj.datum = "NAD83"
    varobj[:] = float(gatts['latitude'])    

    varobj = cdf.createVariable('lon','f8',('lon'))
    varobj.units = "degree_east"
    varobj.epic_code = 502
    # note name is one of the netcdf4 reserved attributes, use setncattr
    varobj.setncattr('name', "LON")
    varobj.long_name = "LONGITUDE"
    varobj.datum = "NAD83"
    varobj[:] = float(gatts['longitude'])    
    
    varobj = cdf.createVariable('bindist','f4',('depth'),fill_value=floatfill)
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
    
    varobj = cdf.createVariable('SV_80','f4',('time','lat','lon'),fill_value=floatfill)
    varobj.units = "m s-1"
    varobj.epic_code = 80
    # note name is one of the netcdf4 reserved attributes, use setncattr
    varobj.setncattr('name', "SV")
    varobj.long_name = "SOUND VELOCITY (M/S)"
    
    varobj = cdf.createVariable('Hdg_1215','f4',('time','lat','lon'),fill_value=floatfill)
    varobj.units = "degrees"
    # note name is one of the netcdf4 reserved attributes, use setncattr
    varobj.setncattr('name', "Hdg")
    varobj.long_name = "INST Heading"
    varobj.epic_code = 1215
    #varobj.heading_alignment = rawvarobj.heading_alignment
    #varobj.heading_bias = rawvarobj.heading_bias
    
    varobj = cdf.createVariable('Ptch_1216','f4',('time','lat','lon'),fill_value=floatfill)
    varobj.units = "degrees"
    varobj.long_name = "INST Pitch"
    varobj.epic_code = 1216
    
    varobj = cdf.createVariable('Roll_1217','f4',('time','lat','lon'),fill_value=floatfill)
    varobj.units = "degrees"
    varobj.long_name = "INST Roll"
    varobj.epic_code = 1217
    
    varobj = cdf.createVariable('Tx_1211','f4',('time','lat','lon'),fill_value=floatfill)
    varobj.units = "C"
    # note name is one of the netcdf4 reserved attributes, use setncattr
    varobj.setncattr('name', "T")
    varobj.long_name = "instrument Transducer Temp."
    varobj.epic_code = 1211

    if 'Pressure' in rawvars:
        #rawvarobj = rawcdf.variables['Pressure']
        varobj = cdf.createVariable('P_1','f4',('time','lat','lon'),fill_value=floatfill)
        varobj.units = "dbar"
        # note name is one of the netcdf4 reserved attributes, use setncattr
        varobj.setncattr('name', "P")
        varobj.long_name = "PRESSURE (DB)"
        varobj.epic_code = 1
        
    if 'PressVar' in rawvars:
        varobj = cdf.createVariable('SDP_850','f4',('time','lat','lon'),fill_value=floatfill)
        varobj.setncattr('name','SDP')
        varobj.long_name = "STAND. DEV. (PRESS)"
        varobj.units = "mbar"
        varobj.epic_code = 850

    varobj = cdf.createVariable('cor','u2',('time','depth','lat','lon'),fill_value=intfill)
    varobj.setncattr('name','cor')
    varobj.long_name = "Slant Beam Average Correlation (cor)"
    varobj.units = "counts"
    varobj.epic_code = 1202
    varobj.NOTE = "Calculated from the slant beams"

    if 'PGd4' in rawvars:
        varobj = cdf.createVariable('PGd_1203','u2',('time','depth','lat','lon'),fill_value=intfill)
        varobj.setncattr('name','Pgd')
        varobj.long_name = "Percent Good Pings"
        varobj.units = "percent"
        varobj.epic_code = 1203
        varobj.NOTE = "Percentage of good 4-bem solutions (Field #4)"

    varobj = cdf.createVariable('AGC_1202','u2',('time','depth','lat','lon'),fill_value=intfill)
    varobj.setncattr('name','AGC')
    varobj.long_name = "Average Echo Intensity (AGC)"
    varobj.units = "counts"
    varobj.epic_code = 1202
    varobj.NOTE = "Calculated from the slant beams"

    if 'cor5' in rawvars:
        varobj = cdf.createVariable('corvert','u2',('time','depth','lat','lon'),fill_value=intfill)
        varobj.setncattr('name','cor')
        varobj.long_name = "Vertical Beam Correlation (cor)"
        varobj.units = "counts"
        varobj.epic_code = 1202
        varobj.NOTE = "From the center vertical beam"

    if 'att5' in rawvars:
        varobj = cdf.createVariable('AGCvert','u2',('time','depth','lat','lon'),fill_value=intfill)
        varobj.setncattr('name','AGC')
        varobj.long_name = "Vertical Beam Echo Intensity (AGC)"
        varobj.units = "counts"
        varobj.epic_code = 1202
        varobj.NOTE = "From the center vertical beam"
        
    # repeating attributes that do not depend on height or depth calculations
    cdfvarnames = []
    for key in cdf.variables.keys(): cdfvarnames.append(key)
    omitnames = []
    for key in cdf.dimensions.keys(): omitnames.append(key)
    omitnames.append("Rec")
    for varname in cdfvarnames:
        if varname not in omitnames:
            varobj = cdf.variables[varname]
            varobj.serial_number = cdf.serial_number

    if settings['transformation'].upper() == "BEAM":
        varnames = ["Beam1", "Beam2", "Beam3", "Beam4"]
        codes = [0,0,0,0]
    elif settings['transformation'].upper() == "INST":
        varnames = ["X","Y","Z","Error"]            
        codes = [0,0,0,0]
    else:
        varnames = ["u_1205","v_1206","w_1204","Werr_1201"]
        codes = [1205, 1206, 1204, 1201]
        
    for i in range(4):
        varname = varnames[i]
        varobj = cdf.createVariable(varname,'f4',('time','depth','lat','lon'),fill_value=floatfill)
        varobj.units = "cm s-1"
        varobj.long_name = "%s velocity (cm s-1)" % varnames[i]
        varobj.epic_code = codes[i]

    if 'vel5' in rawvars:
        varobj = cdf.createVariable('Wvert','f4',('time','depth','lat','lon'),fill_value=floatfill)
        varobj.units = "cm s-1"
        varobj.long_name = "Vertical velocity (cm s-1)" 

    # TODO do we do bottom track data here?  Later?  Or as a separate thing?

    add_VAR_DESC(cdf)
        
    return cdf


def add_VAR_DESC(cdf):
    # cdf is an netcdf file object (e.g. pointer to open netcdf file)

    varkeys = cdf.variables.keys() #get the names
    dimkeys = cdf.dimensions.keys()
    varnames = []
    for key in varkeys: varnames.append(key)
    dimnames = []
    for key in dimkeys: dimnames.append(key)
    buf = ""
    for varname in varnames:
        if varname not in dimnames:
            buf = "%s:%s" % (buf,varname)

    cdf.VAR_DESC = buf  
    
    return
    
    
def read_globalatts(fname):        
    # read_globalatts: read in file of metadata for a tripod or mooring
    # usage : gatt=read_globalatts(fname)
    #
    # reads global attributes for an experiment from a text file (fname)
    # called by all data processing programs to get uniform metadata input
    #  one argument is required- the name of the file to read- it should have
    #  this form:
    # SciPi;C. Sherwood
    # PROJECT; ONR
    # EXPERIMENT; RIPPLES DRI
    # DESCRIPTION; Stress, SSC, and Bedforms at MVCO 12-m fine/coarse transition site
    # DATA_SUBTYPE; MOORED
    # DATA_ORIGIN; USGS WHFS Sed Trans Group
    # COORD_SYSTEM; GEOGRAPHIC
    # Conventions; PMEL/EPIC
    # MOORING; 836
    # WATER_DEPTH; 10.99
    # latitude; 41.336063
    # longitude; -70.559615
    # magnetic_variation; -15
    # Deployment_date; 27-Aug-2007
    # Recovery_date;  ?
    gatts = {}
    f = open(fname,'r')
    for line in f:
        line = line.strip()
        cols = line.split(";")
        gatts[cols[0]] = cols[1].strip()
        
        
    f.close()
    return gatts
    
# write a dictionary to netCDF attributes
def writeDict2atts(cdfobj, d, tag):
    
    i = 0;
    # first, convert as many of the values in d to numbers as we can
    for key in iter(d):
        if type(d[key]) == str:
            try:
                d[key] = float(d[key])
            except ValueError:
                # we really don't need to print here, 
                # but python insists we do something
                #print('   can\'t convert %s to float' % key)
                i += 1

    for key in iter(d):
        newkey = tag + key
        try:
            cdfobj.setncattr(newkey,d[key])
        except:
            print('can\'t set %s attribute' % key)
    
    # return the dictionary it's numerized values 
    return d

def floor(dec):
    # for convenience and to avoid loading math
    # anf then because np.floor was causing unexpected behavior w.r.t ints
    return int(dec - (dec % 1))
    
def s2hms(secs):
    hour = floor(secs/3600)
    mn = floor((secs % 3600)/60)
    sec = secs % 60
    return hour, mn, sec

def EPICtime2datetime(time,time2):
    
    dtos = []
    gtime = []
    for idx in range(len(time)):
        # time and time2 are the julian day and milliseconds 
        # in the day as per PMEL EPIC convention for netcdf
        jd = time[idx]+(time2[idx]/(24*3600*1000))
        secs = (jd % 1)*(24*3600)
        
        j = floor(jd) - 1721119
        in1 = 4*j-1
        y = floor(in1/146097)
        j = in1 - 146097*y
        in1 = floor(j/4)
        in1 = 4*in1 +3
        j = floor(in1/1461)
        d = floor(((in1 - 1461*j) +4)/4)
        in1 = 5*d -3
        m = floor(in1/153)
        d = floor(((in1 - 153*m) +5)/5)
        y = y*100 +j
        mo=m-9
        yr=y+1
        if m<10:
            mo = m+3
            yr = y
        hour, mn, sec = s2hms(secs)
        ss = floor(sec)
        
        hund = floor((sec-ss)*100)

        gtime.append([yr, mo, d, hour, mn, ss, hund])

        # centiseconds * 10000 = microseconds
        dto = dt.datetime(yr, mo, d, hour, mn, ss, int(hund*10000))
        
        dtos.append(dto)

    return gtime, dtos

def cftime2EPICtime(timecount, timeunits):
    # take a CF time variable and convert to EPIC time and time2
    # timecountis the integer count of minutes (for instance) since the time stamp
    # given in timeunits
    buf = timeunits.split()
    t0 = dt.datetime.strptime(buf[2]+' '+buf[3], '%Y-%m-%d %H:%M:%S')
    t0j = ajd(t0)
    # julian day for EPIC is the beginning of the day e.g. midnight
    t0j = t0j+0.5 # add 0.5 because ajd() subtracts 0.5 
    
    if buf[0] == 'hours':
        tj = timecount/(24)
    elif buf[0] == 'minutes':
        tj = timecount/(24*60)
    elif buf[0] == 'seconds':
        tj = timecount/(24*60*60)
    elif buf[0] == 'milliseconds':
        tj = timecount/(24*60*60*1000)
    elif buf[0] == 'microseconds':
        tj = timecount/(24*60*60*1000*1000)
        
    tj = t0j+tj
    
    time = np.floor(tj)
    time2 = np.floor((tj-time)*(24*3600*1000))
    
    return time, time2

"""
# this does not work
def cf2EPICtime(cftime, cfunits, cfcalendar):
    tobj = num2date(cftime,cfunits,calendar=cfcalendar)
    elapsed_sec = []
    for idx in range(len(tobj)):
        tdelta = tobj[idx]-tobj[0] # timedelta
        elapsed_sec.append(tdelta.total_seconds())
    # from the datetime object convert to time and time2
    jd = []
    time = []
    time2 = []
    for idx in range(len(tobj)):
        j = julian(tobj[idx].year,tobj[idx].month,tobj[idx].day, \
                   tobj[idx].hour,tobj[idx].minute,tobj[idx].second,\
                   floor(tobj[idx].microsecond/1000000))
        jd.append(j)
        time.append(int(floor(j)))
        time2.append(int((j - floor(j))*(24*3600*1000)))
        
    return time, time2
"""

def __main():
    print('%s running on python %s' % (sys.argv[0], sys.version))
	
    if len(sys.argv) < 2:
        print("%s useage:" % sys.argv[0])
        print("ADCPcdf2ncEPIC rawcdfname ncEPICname USGSattfile [startingensemble endingensemble]\n" )
        print("starting and ending ensemble are netcdf file indeces, NOT TRDI ensemble numbers")
        print("USGSattfile is a file containing EPIC metadata")
        sys.exit(1)
        
    try:
        infileName = sys.argv[1]
    except:
        print('error - input file name missing')
        sys.exit(1)
        
    try:
        outfileName = sys.argv[2]
    except:
        print('error - output file name missing')
        sys.exit(1)
        
    try:
        attfileName = sys.argv[3]
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
    if ~os.path.isfile(infileName):
        print('error - input file not found')
        sys.exit(1)
        
    if ~os.path.isfile(attfileName):
        print('error - attribute file not found')
        sys.exit(1)

    print('Converting %s to %s' % (infileName, outfileName))
    
    print('Start file conversion at ',dt.datetime.now())
    
    doEPIC_ADCPfile(infileName, outfileName, attfileName, settings)
    
    print('Finished file conversion at ',dt.datetime.now())
    


if __name__ == "__main__":
    __main()