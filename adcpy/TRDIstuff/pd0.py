#!/usr/bin/python
"""
Functions to handle Acoustic Doppler Current Profiler data from
a Teledyne RD Instruments instrument that is in pd0 format

As a script (to split a raw PD0 file into waves and currents):

python pd0.py [path] rawFile wavesFile currentsFile
python pd0.py [-p path] -r rawFile -w wavesFile -c currentsFile
python pd0.py [--path=path] \\
               --raw=rawFile \\
               --waves=wavesFile \\
               --currents=currentsFile

where:
    path         is a path to prepend to the following
    rawFile      is path of raw PD0 format input file
    wavesFile    is path of waves PD0 format output file
    currentsFile is path of currents PD0 format output file

or (to run the test suite):

python pd0.py -t
python pd0.py --test

or (to see this help message):

python pd0.py -h
python pd0.py --help

As a module:

import adcp.pd0
adcp.pd0.split(rawFile,wavesFile,currentsFile)

where:
    rawFile      is a file object representing the raw PD0 format input
    wavesFile    is a file object representing the waves PD0 format output
    currentsFile is a file object representing the currents PD0 format output
"""

import struct

#
# The rawfile is assumed to be in PD0 format.
#
# PD0 format assumes the file is a succession of ensembles.
#
# Each ensemble starts with a two byte header identifying the type of data
# contained in the ensemble.
#
# Following the header is a two byte length field specifying the length of
# the header, length field, and data combined
#
# Following the length field is raw data for the number of bytes indicated by
# the length field
#
# Following the raw data is a checksum field which is the two least
# significant bytes of the sum of the byte values of the header, length field,
# and raw data.
#

# updated to run in python 3x, Marinna Martini 1/12/2017
# adapted from pd0.py by Gregory P. Dusek
# http://trac.nccoos.org/dataproc/wiki/DPWP/docs


def split(raw_file, waves_file, currents_file):
    """
    Split PD0 format data into seperate waves and currents files

    :param binaryIO raw_file:
    :param binaryIO waves_file:
    :param binaryIO currents_file:
    :return:
    """

    # header IDs
    waves_id = 0x797F
    currents_id = 0x7F7F

    # convenience function reused for header, length, and checksum
    def __nextLittleEndianUnsignedShort(file):
        """Get next little endian unsigned short from file"""
        raw = file.read(2)
        """for python 3.5, struct.unpack('<H', raw)[0] needs to return a
           byte, not an int
        """
        return raw, struct.unpack("<H", raw)[0]

    # factored for readability
    def __computeChecksum(data, nbytes, ensemble):
        """Compute a checksum from header, length, and ensemble"""
        cs = 0
        for byte in data:
            # since the for loop returns an int to byte, use as-is
            # value = struct.unpack('B', byte)[0]
            # cs += value
            cs += byte
        for byte in nbytes:
            # value = struct.unpack('B', byte)[0]
            # cs += value
            cs += byte
        for byte in ensemble:
            # value = struct.unpack('B', byte)[0]
            # cs += value
            cs += byte
        return cs & 0xFFFF

    # find the first instance of a waves or currents header
    raw_data = raw_file.read()
    first_waves = raw_data.find(struct.pack("<H", waves_id))
    first_currents = raw_data.find(struct.pack("<H", currents_id))

    # bail if neither waves nor currents found
    if (first_waves < 0) and (first_currents < 0):
        #    raise IOError, "Neither waves nor currents header found"
        raise IOError("Neither waves nor currents header found")

    # get the starting point by throwing out unfound headers
    # and selecting the minimum
    first_ensemble = min([x for x in (first_waves, first_currents) if x >= 0])

    # seeks to the first occurrence of a waves or currents data
    raw_file.seek(first_ensemble)

    # loop through raw data
    raw_header, header = __nextLittleEndianUnsignedShort(raw_file)

    while (header == waves_id) or (header == currents_id):
        # get ensemble length
        raw_length, length = __nextLittleEndianUnsignedShort(raw_file)
        # read up to the checksum
        raw_ensemble = raw_file.read(length - 4)
        # get checksum
        raw_checksum, checksum = __nextLittleEndianUnsignedShort(raw_file)

        computed_checksum = __computeChecksum(raw_header, raw_length, raw_ensemble)

        if checksum != computed_checksum:
            raise IOError("Checksum error")

        # append to output stream
        if header == waves_id:
            waves_file.write(raw_header)
            waves_file.write(raw_length)
            waves_file.write(raw_ensemble)
            waves_file.write(raw_checksum)
        elif header == currents_id:
            currents_file.write(raw_header)
            currents_file.write(raw_length)
            currents_file.write(raw_ensemble)
            currents_file.write(raw_checksum)

        try:
            raw_header, header = __nextLittleEndianUnsignedShort(raw_file)
        except struct.error:
            break


def test():
    """Execute test suite"""
    try:
        import adcp.tests.runalltests as runalltests
    except:
        # possible if executed as script
        import sys
        import os

        sys.path.append(os.path.join(os.path.dirname(__file__), "tests"))
        import runalltests
    runalltests.runalltests(subset="pd0")


class __TestException(Exception):
    """Flow control for running as script"""

    pass


# wrapper function
def __test():
    """Execute test suite from command line"""
    test()
    # raise __TestException, 'Wrapper function for command line testing only'
    raise __TestException("Wrapper function for command line testing only")


def __main():
    """Process as script from command line"""
    import getopt
    import os
    import sys

    # get the command line options and arguments
    path = ""
    raw_name, waves_name, currents_name = 3 * [None]

    try:
        opts, args = getopt.gnu_getopt(
            sys.argv[1:],
            "htp:r:w:c:",
            ["help", "test", "path=", "raw=", "waves=", "currents="],
        )
        for opt, arg in opts:
            if opt in ["-h", "--help"]:
                raise getopt.GetoptError("")
            if opt in ["-t", "--test"]:
                __test()
            elif opt in ["-p", "--path"]:
                path = arg
            elif opt in ["-r", "--raw"]:
                raw_name = arg
            elif opt in ["-w", "--waves"]:
                waves_name = arg
            elif opt in ["-c", "--currents"]:
                currents_name = arg
            else:
                raise getopt.GetoptError("")
        if (raw_name is None) or (waves_name is None) or (currents_name is None):
            if len(args) not in [3, 4]:
                raise getopt.GetoptError("")
            else:
                if (
                    (raw_name is not None)
                    or (waves_name is not None)
                    or (currents_name is not None)
                ):
                    raise getopt.GetoptError("")
                else:
                    if len(args) == 4:
                        path = args[0]
                        del args[0]
                    raw_name = args[0]
                    waves_name = args[1]
                    currents_name = args[2]
        elif len(args) != 0:
            raise getopt.GetoptError("")
    except getopt.GetoptError:
        print(__doc__)
        return
    except __TestException:
        return

    # split a raw PD0 file
    raw_name = os.path.join(path, raw_name)
    print(("Raw file path:", raw_name))
    waves_name = os.path.join(path, waves_name)
    print(("Waves file path:", waves_name))
    currents_name = os.path.join(path, currents_name)
    print(("Currents file path:", currents_name))
    raw_file = open(raw_name, "rb")
    try:
        waves_file = open(waves_name, "wb")
        try:
            currents_file = open(currents_name, "wb")
            try:
                split(raw_file, waves_file, currents_file)
            finally:
                currents_file.close()
        finally:
            waves_file.close()
    finally:
        raw_file.close()


if __name__ == "__main__":
    __main()
