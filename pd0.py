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

def split(rawFile,wavesFile,currentsFile):
    """Split PD0 format data into seperate waves and currents

    split()rawFile,wavesFile,currentsFile -> None
    """

    # header IDs
    wavesId=0x797f
    currentsId=0x7f7f

    # convenience function reused for header, length, and checksum
    def __nextLittleEndianUnsignedShort(file):
        """Get next little endian unsigned short from file"""
        raw = file.read(2)
        """for python 3.5, struct.unpack('<H', raw)[0] needs to return a
           byte, not an int
        """
        return (raw, struct.unpack('<H', raw)[0])

    # factored for readability
    def __computeChecksum(header, length, ensemble):
        """Compute a checksum from header, length, and ensemble"""
        cs = 0    
        for byte in header:
            # since the for loop returns an int to byte, use as-is
            #value = struct.unpack('B', byte)[0]
            #cs += value
            cs += byte
        for byte in length:
            #value = struct.unpack('B', byte)[0]
            #cs += value
            cs += byte
        for byte in ensemble:
            #value = struct.unpack('B', byte)[0]
            #cs += value
            cs += byte
        return cs & 0xffff

    # find the first instance of a waves or currents header
    rawData=rawFile.read()
    firstWaves = rawData.find(struct.pack('<H',wavesId))
    firstCurrents= rawData.find(struct.pack('<H',currentsId))

    # bail if neither waves nor currents found
    if (firstWaves < 0) and (firstCurrents < 0):
    #    raise IOError, "Neither waves nor currents header found"
        raise IOError('Neither waves nor currents header found')

    # get the starting point by throwing out unfound headers
    # and selecting the minumum
    firstEnsemble = min([x for x in (firstWaves,firstCurrents) if x >= 0])

    #seeks to the first occurence of a waves or currents data
    rawFile.seek(firstEnsemble)

    # loop through raw data
    rawHeader, header = __nextLittleEndianUnsignedShort(rawFile)

    while (header == wavesId) or (header == currentsId):
        # get ensemble length
        rawLength, length = __nextLittleEndianUnsignedShort(rawFile)
        # read up to the checksum
        rawEnsemble = rawFile.read(length-4)
        # get checksum
        rawChecksum, checksum = __nextLittleEndianUnsignedShort(rawFile)

        computedChecksum = __computeChecksum(rawHeader, rawLength, rawEnsemble)

        if checksum != computedChecksum:
            raise IOError('Checksum error')

        # append to output stream
        if header == wavesId: 
            wavesFile.write(rawHeader)
            wavesFile.write(rawLength)
            wavesFile.write(rawEnsemble)
            wavesFile.write(rawChecksum)
        elif header == currentsId:
            currentsFile.write(rawHeader)
            currentsFile.write(rawLength)
            currentsFile.write(rawEnsemble)
            currentsFile.write(rawChecksum)

        try:
            rawHeader, header = __nextLittleEndianUnsignedShort(rawFile)
        except struct.error:
            break


def test():
    """Execute test suite"""
    try:
        import adcp.tests.runalltests as runalltests
    except:
        # possible if executed as script
        import sys,os
        sys.path.append(os.path.join(os.path.dirname(__file__),'tests'))
        import runalltests
    runalltests.runalltests(subset='pd0')

class __TestException(Exception):
    """Flow control for running as script"""
    pass

# wrapper function
def __test():
   """Execute test suite from command line"""
   test()
   # raise __TestException, 'Wrapper function for command line testing only'
   raise __TestException('Wrapper function for command line testing only')

def __main():
    """Process as script from command line"""
    import getopt,os,sys

    # get the command line options and arguments
    path = ''
    rawName,wavesName,currentsName = 3*[None]

    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:],
                                       'htp:r:w:c:',
                                       ['help',
                                        'test',
                                        'path=',
                                        'raw=',
                                        'waves=',
                                        'currents='])
        for opt,arg in opts:
            if opt in ['-h','--help']:
                raise getopt.GetoptError('')
            if opt in ['-t','--test']:
                __test()
            elif opt in ['-p','--path']:
                path = arg
            elif opt in ['-r','--raw']:
                rawName = arg
            elif opt in ['-w','--waves']:
                wavesName = arg
            elif opt in ['-c','--currents']:
                currentsName = arg
            else:
                raise getopt.GetoptError('')
        if (rawName is None) or \
           (wavesName is None) or \
           (currentsName is None):
            if len(args) not in [3, 4]:
                raise getopt.GetoptError('')
            else:
                if (rawName is not None) or \
                   (wavesName is not None) or \
                   (currentsName is not None):
                    raise getopt.GetoptError('')
                else:
                    if len(args) == 4:
                        path = args[0]
                        del args[0]
                    rawName = args[0]
                    wavesName = args[1]
                    currentsName = args[2]
        elif len(args) != 0:
            raise getopt.GetoptError('')
    except getopt.GetoptError:
        print (__doc__)
        return
    except __TestException:
        return

    # split a raw PD0 file
    rawName = os.path.join(path, rawName)
    print(('Raw file path:', rawName))
    wavesName = os.path.join(path, wavesName)
    print(('Waves file path:', wavesName))
    currentsName = os.path.join(path, currentsName)
    print(('Currents file path:', currentsName))
    rawFile = open(rawName, 'rb')
    try:
        wavesFile = open(wavesName, 'wb')
        try:
            currentsFile = open(currentsName, 'wb')
            try:
                split(rawFile, wavesFile, currentsFile)
            finally:
                currentsFile.close()
        finally:
            wavesFile.close()
    finally:
        rawFile.close()

if __name__ == "__main__":
    __main()

