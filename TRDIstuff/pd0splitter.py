"""
pd0splitter
===========

Functions to handle Acoustic Doppler Current Profiler data Teledyne RD Instruments pd0 format

As a script (to split a raw PD0 file into waves and currents)

Command line usage:
    python pd0splitter.py [path] rawFile wavesFile currentsFile

    python pd0splitter.py [-p path] -r rawFile -w wavesFile -c currentsFile -f firstEnsemble -l lastEnsemble

    python pd0splitter.py [--path=path] --raw=rawFile --waves=wavesFile --currents=currentsFile --first=firstEnsemble \\
        --last=lastEnsemble

    python pd0splitter.py -h

    python pd0splitter.py --help

Parameters:
    path:         is a path to prepend to the following
    rawFile:      is path of raw PD0 format input file
    wavesFile:    is path of waves PD0 format output file
    currentsFile: is path of currents PD0 format output file
    firstEnsemble:   is a integer ensemble number
    lastEnsemble:   is a integer ensemble number

The rawfile is assumed to be in PD0 format.  PD0 format assumes the file is a succession of ensembles.
Each ensemble starts with a two byte header identifying the type of data contained in the ensemble.

Following the header is a two byte length field specifying the length of the header, length field, and data combined

Following the length field is raw data for the number of bytes indicated by the length field

Following the raw data is a checksum field which is the two least significant bytes of the sum of the byte values
of the header, length field, and raw data.

updated to run in python 3x, Marinna Martini 1/12/2017
adapted from pd0.py by Gregory P. Dusek http://trac.nccoos.org/dataproc/wiki/DPWP/docs
"""

import getopt,os,sys
import struct


def split(rawFile, wavesFile, currentsFile, firstEnsemble, lastEnsemble):
    """
    split ADCP data in pd0 format into current profiles and wave packets

    :param str rawFile: path and name of raw PD0 format input file
    :param str wavesFile: path and name of waves PD0 format output file
    :param str currentsFile: path and name of currents PD0 format output file
    :param int firstEnsemble: ensemble number of the first ensemble to read
    :param int lastEnsemble: ensemble number of the last ensemble to read
    """
    try:
        rawFile = open(rawFile, 'rb')
    except:
        print('Cannot open %s' % rawFile)
    try:
        wavesFile = open(wavesFile, 'wb')
    except:
        print('Cannot open %s' % wavesFile)
    try:
        currentsFile = open(currentsFile, 'wb')
    except:
        print('Cannot open %s' % currentsFile)
                
    # header IDs
    wavesId = 0x797f
    currentsId = 0x7f7f
    
    if lastEnsemble < 0:
        lastEnsemble = 1E35
    
    print('Reading from %d to %d ensembles\n' % (firstEnsemble, lastEnsemble))

    # TODO move these function definitions out of this function
    # convenience function reused for header, length, and checksum
    def __nextLittleEndianUnsignedShort(file):
        """
        Get next little endian unsigned short from file

        :param file: file object open for reading as binary
        :return: a tuple of raw bytes and unpacked data
        """
        raw = file.read(2)
        # for python 3.5, struct.unpack('<H', raw)[0] needs to return a byte, not an int

        return (raw, struct.unpack('<H', raw)[0])

    # factored for readability
    def __computeChecksum(header, length, ensemble):
        """
        Compute a checksum

        :param header: file header
        :param length: ensemble length
        :param ensemble: ensemble raw data
        :return: checksum for ensemble
        """
        cs = 0    
        for byte in header:
            # since the for loop returns an int to byte, use as-is
            # value = struct.unpack('B', byte)[0]
            # cs += value
            cs += byte
        for byte in length:
            # value = struct.unpack('B', byte)[0]
            # cs += value
            cs += byte
        for byte in ensemble:
            # value = struct.unpack('B', byte)[0]
            # cs += value
            cs += byte
        return cs & 0xffff

    # find the first instance of a waves or currents header
    rawData=rawFile.read()
    firstWaves = rawData.find(struct.pack('<H', wavesId))
    firstCurrents= rawData.find(struct.pack('<H', currentsId))

    # bail if neither waves nor currents found
    if (firstWaves < 0) and (firstCurrents < 0):
        # raise IOError, "Neither waves nor currents header found"
        raise IOError('Neither waves nor currents header found')

    # get the starting point by throwing out unfound headers
    # and selecting the minimum
    firstFileEnsemble = min([x for x in (firstWaves, firstCurrents) if x >= 0])

    # seeks to the first occurence of a waves or currents data
    rawFile.seek(firstFileEnsemble)

    # loop through raw data
    rawHeader, header = __nextLittleEndianUnsignedShort(rawFile)
    
    waveCount = 0
    currentCount = 0

    while (header == wavesId) or (header == currentsId) or currentFlag:
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
            waveCount = waveCount+1
            wavesFile.write(rawHeader)
            wavesFile.write(rawLength)
            wavesFile.write(rawEnsemble)
            wavesFile.write(rawChecksum)
        elif header == currentsId:
            currentCount = currentCount+1
            if (currentCount >= firstEnsemble) & (currentCount < lastEnsemble):
                currentsFile.write(rawHeader)
                currentsFile.write(rawLength)
                currentsFile.write(rawEnsemble)
                currentsFile.write(rawChecksum)
            elif currentCount > lastEnsemble:
                break

        try:
            rawHeader, header = __nextLittleEndianUnsignedShort(rawFile)
        except struct.error:
            break
        
        if (currentCount > 0) & ((currentCount % 100) == 0):
            print('%d current ensembles read' % currentCount)
        if (waveCount > 0) & ((waveCount % 1000) == 0):
            print('%d wave ensembles read' % waveCount)

    print('wave Ensemble count = %d\n' % waveCount)
    print('current Ensemble count = %d\n' % currentCount)

    currentsFile.close()
    wavesFile.close()
    rawFile.close()


def __main():

    # get the command line options and arguments
    path = ''
    rawFile, wavesFile, currentsFile, firstEnsemble, lastEnsemble = 5*[None]

    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:],
                                       'htp:r:w:c:f:l',
                                       ['help',
                                        'path=',
                                        'raw=',
                                        'waves=',
                                        'currents=',
                                        'first=',
                                        'last-'])
        for opt,arg in opts:
            if opt in ['-h', '--help']:
                raise getopt.GetoptError('')
            elif opt in ['-p', '--path']:
                path = arg
            elif opt in ['-r', '--raw']:
                rawFile = arg
            elif opt in ['-w', '--waves']:
                wavesFile = arg
            elif opt in ['-c', '--currents']:
                currentsFile = arg
            elif opt in ['-f', '--first']:
                firstEnsemble = arg
            elif opt in ['-l', '--last']:
                lastEnsemble = arg
            else:
                raise getopt.GetoptError('')
        if (rawFile is None) or \
           (wavesFile is None) or \
           (currentsFile is None):
            if len(args) not in [3, 4, 5, 6]:
                raise getopt.GetoptError('')
            else:
                if (rawFile is not None) or \
                   (wavesFile is not None) or \
                   (currentsFile is not None):
                    raise getopt.GetoptError('')
                else:
                    if len(args) in [4, 5, 6]:
                        path = args[0]
                        del args[0]
                    rawFile = args[0]
                    wavesFile = args[1]
                    currentsFile = args[2]
                    if len(args) in [5, 6]:
                        firstEnsemble = args[3]
                        lastEnsemble = args[4]
        elif len(args) != 0:
            raise getopt.GetoptError('')
    except getopt.GetoptError:
        print (__doc__)
        return

    # split a raw PD0 file
    rawFile = os.path.join(path, rawFile)
    print(('Raw file path:', rawFile))
    wavesFile = os.path.join(path, wavesFile)
    print(('Waves file path:', wavesFile))
    currentsFile = os.path.join(path, currentsFile)
    print(('Currents file path:', currentsFile))

    split(rawFile, wavesFile, currentsFile, firstEnsemble, lastEnsemble)


if __name__ == "__main__":
    __main()

