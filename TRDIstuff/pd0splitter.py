"""
pd0splitter
===========

Functions to handle Acoustic Doppler Current Profiler data Teledyne RD Instruments pd0 format

As a script (to split a raw PD0 file into waves packets and currents binary files)

Usage:
    python pd0splitter.py pd0file packets_file currents_file [first_ensemble] [last_ensemble]

    :param str pd0file: is path of raw PD0 format input file
    :param str packets_file: is path of waves PD0 format output file
    :param str currents_file: is path of currents PD0 format output file
    :param int first_ensemble: is a integer ensemble number
    :param int last_ensemble: is a integer ensemble number

The pd0file is assumed to be in PD0 format.  PD0 format assumes the file is a succession of ensembles.
Each ensemble starts with a two byte header identifying the type of data contained in the ensemble.

Following the header is a two byte length field specifying the length of the header, length field, and data combined

Following the length field is raw data for the number of bytes indicated by the length field

Following the raw data is a checksum field which is the two least significant bytes of the sum of the byte values
of the header, length field, and raw data.

updated to run in python 3x, Marinna Martini 1/12/2017
adapted from pd0.py by Gregory P. Dusek http://trac.nccoos.org/dataproc/wiki/DPWP/docs
"""

import sys
import struct


def split(pd0file, packets_file, currents_file, first_ensemble, last_ensemble):
    """
    split ADCP data in pd0 format into current profiles and wave packets

    :param str pd0file: path and name of raw PD0 format input file
    :param str packets_file: path and name of waves PD0 format output file
    :param str currents_file: path and name of currents PD0 format output file
    :param int first_ensemble: ensemble number of the first ensemble to read
    :param int last_ensemble: ensemble number of the last ensemble to read
    """
    try:
        pd0file = open(pd0file, 'rb')
    except:
        print('Cannot open %s' % pd0file)
    try:
        packets_file = open(packets_file, 'wb')
    except:
        print('Cannot open %s' % packets_file)
    try:
        currents_file = open(currents_file, 'wb')
    except:
        print('Cannot open %s' % currents_file)
                
    # header IDs
    waves_id = 0x797f
    currents_id = 0x7f7f
    
    if last_ensemble < 0:
        last_ensemble = 1E35
        print('Reading from %d to the last ensemble found\n' % first_ensemble)
    else:
        print('Reading from ensemble %d to %d\n' % (first_ensemble, last_ensemble))

    # find the first instance of a waves or currents header
    raw_data = pd0file.read()
    first_waves = raw_data.find(struct.pack('<H', waves_id))
    first_currents = raw_data.find(struct.pack('<H', currents_id))

    # bail if neither waves nor currents found
    if (first_waves < 0) and (first_currents < 0):
        # raise IOError, "Neither waves nor currents header found"
        raise IOError('Neither waves nor currents header found')

    # get the starting point by throwing out unknown headers
    # and selecting the minimum
    first_file_ensemble = min([x for x in (first_waves, first_currents) if x >= 0])

    # seeks to the first occurrence of a waves or currents data
    pd0file.seek(first_file_ensemble)

    # loop through raw data
    raw_header, header = __nextLittleEndianUnsignedShort(pd0file)
    
    wave_count = 0
    current_count = 0

    while (header == waves_id) or (header == currents_id):
        # get ensemble length
        raw_length, length = __nextLittleEndianUnsignedShort(pd0file)
        # read up to the checksum
        raw_ensemble = pd0file.read(length-4)
        # get checksum
        raw_checksum, checksum = __nextLittleEndianUnsignedShort(pd0file)

        computed_checksum = __computeChecksum(raw_header, raw_length, raw_ensemble)

        if checksum != computed_checksum:
            raise IOError('Checksum error')

        # append to output stream
        if header == waves_id:
            wave_count = wave_count+1
            packets_file.write(raw_header)
            packets_file.write(raw_length)
            packets_file.write(raw_ensemble)
            packets_file.write(raw_checksum)
        elif header == currents_id:
            current_count = current_count+1
            if (current_count >= first_ensemble) & (current_count < last_ensemble):
                currents_file.write(raw_header)
                currents_file.write(raw_length)
                currents_file.write(raw_ensemble)
                currents_file.write(raw_checksum)
            elif current_count > last_ensemble:
                break

        try:
            raw_header, header = __nextLittleEndianUnsignedShort(pd0file)
        except struct.error:
            break
        
        if (current_count > 0) & ((current_count % 100) == 0):
            print('%d current ensembles read' % current_count)
        if (wave_count > 0) & ((wave_count % 1000) == 0):
            print('%d wave ensembles read' % wave_count)

    print('wave Ensemble count = %d\n' % wave_count)
    print('current Ensemble count = %d\n' % current_count)

    currents_file.close()
    packets_file.close()
    pd0file.close()


# convenience function reused for header, length, and checksum
def __nextLittleEndianUnsignedShort(file):
    """
    Get next little endian unsigned short from file

    :param file: file object open for reading as binary
    :return: a tuple of raw bytes and unpacked data
    """
    raw = file.read(2)
    # for python 3.5, struct.unpack('<H', raw)[0] needs to return a byte, not an int

    return raw, struct.unpack('<H', raw)[0]


# factored for readability
def __computeChecksum(file_header, ensemble_length, ensemble):
    """
    Compute a checksum

    :param file_header: file header
    :param ensemble_length: ensemble ensemble_length
    :param ensemble: ensemble raw data
    :return: checksum for ensemble
    """
    cs = 0
    for byte in file_header:
        # since the for loop returns an int to byte, use as-is
        # value = struct.unpack('B', byte)[0]
        # cs += value
        cs += byte
    for byte in ensemble_length:
        # value = struct.unpack('B', byte)[0]
        # cs += value
        cs += byte
    for byte in ensemble:
        # value = struct.unpack('B', byte)[0]
        # cs += value
        cs += byte
    return cs & 0xffff


def __main():

    print('%s running on python %s' % (sys.argv[0], sys.version))
    if len(sys.argv) < 3:
        print(__doc__)
        return

    try:
        pd0file = sys.argv[1]
    except:
        print('error - pd0 input file name missing')
        sys.exit(1)

    try:
        packets_file = sys.argv[2]
    except:
        print('error - packets output file name missing')
        sys.exit(1)

    try:
        currents_file = sys.argv[3]
    except:
        print('error - current profile output file name missing')
        sys.exit(1)

    print('Splitting %s to %s and %s' % (pd0file, packets_file, currents_file))

    try:
        first_ensemble = sys.argv[4]
    except:
        first_ensemble = 1

    try:
        last_ensemble = sys.argv[5]
    except:
        last_ensemble = -1

    split(pd0file, packets_file, currents_file, first_ensemble, last_ensemble)


if __name__ == "__main__":
    __main()
