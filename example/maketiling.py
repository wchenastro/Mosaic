#!/usr/bin/env python

import numpy as np
import os, sys, datetime
import argparse
import logging

import matplotlib.pyplot as plt
loggerFormat = '%(asctime)-15s %(filename)s  %(message)s'
logging.basicConfig(format = loggerFormat, level=logging.WARNING)
logger = logging.getLogger()

import katpoint
from mosaic import PsfSim, generate_nbeams_tiling, generate_radius_tiling
from mosaic.coordinate import convertBoresightToDegree, convert_equatorial_coordinate_to_pixel
from mosaic.plot import plot_overlap, plot_interferometry



def makeKatPointAntenna(antennaString):
    antennaKat = []

    for antenna in antennaString:
        antkat = katpoint.Antenna(antenna)
        antennaKat.append(antkat)

    return antennaKat

def creatBeamMatrix(antennaCoords, sourceCoord, observeTime, frequencies,
        duration, overlap, beamNum, subarray, update_interval, hourAngle,
        retile, retileRate, retile_threshold, overlay_source, plotting):

    if subarray != []:
        antennaKat = makeKatPointAntenna(
                [antennaCoords[int(ant)] for ant in subarray])
    else:
        antennaKat = makeKatPointAntenna(antennaCoords)

    # boresight = sourceCoord
    boresight = katpoint.Target('boresight, radec, {}, {}'.format(
                sourceCoord[0], sourceCoord[1]))
    reference = makeKatPointAntenna(
            ["ref, -30:42:39.8, 21:26:38.0, 1035.0",])[0]
    psf = PsfSim(antennaKat, frequencies[0])
    beamShape = psf.get_beam_shape(boresight, observeTime)
    beamShape.plot_psf("beamWithFit.png", shape_overlay=True)


    tiling = generate_nbeams_tiling(
            beamShape, beamNum, overlap = overlap)
    tiling_coordinats = tiling.get_equatorial_coordinates()
    np.savetxt("tilingCoord", tiling_coordinats)

    if plotting is not True:
        exit(0)
    else:
        outputPath = 'plots'
        if(not os.path.isdir(outputPath)):
            os.mkdir(outputPath)

    counter = []
    step = update_interval/60.
    hour = duration/3600.
    grids = []
    horizons = []
    retileTime = []
    startTime = observeTime
    localHourAngle = []
    retileStamp = []
    for minute in np.arange(0, hour*60, step):
        offsetTime = startTime + datetime.timedelta(minutes=minute)
        newBeamShape = psf.get_beam_shape(sourceCoord, offsetTime)
        tiling.beam_shape = newBeamShape
        extra_names = []
        extra_coords_pixel = []
        if overlay_source is not None:
            extra_coords = []
            for source in overlay_source:
                extra_names.append(source[0])
                extra_coords.append([source[1],source[2]])
            extra_coords_pixel = convert_equatorial_coordinate_to_pixel(
                    extra_coords, newBeamShape.bore_sight.equatorial)
        tiling.plot_tiling("plots/tiling%03d.svg" % minute,
            HD= True, index = True,
            extra_coordinates = extra_coords_pixel,
            extra_coordinates_text = extra_names)


def parseOptions(parser):
    # enable interpolation when ploting
    parser.add_argument('--inte', action='store_true', help=argparse.SUPPRESS)
    # plot the beam images
    parser.add_argument('--plot', action='store_true', help='plot the result')
    parser.add_argument('--ants', nargs=1, metavar="file", help='antenna coodinates files')
    parser.add_argument('--resolution', nargs=1, metavar="asec", help='resolution in arcsecond')
    parser.add_argument('--size', nargs=1, metavar="num", help='width in pixels')
    parser.add_argument('--freqrange', nargs=3, metavar=('s', 't', 'i'), help='freqencies range as start stop interval')
    parser.add_argument('--freq', nargs='+',  help='multiple freqencies')
    parser.add_argument('--frame', nargs=1, metavar="RADEC/AziAlt", help='source coordinate frame')
    parser.add_argument('--source', nargs=2, metavar=('RA/Azi', 'DEC/Alt'), help='source position in RADEC or AziAlt')
    parser.add_argument('--datetime', nargs=2, metavar=('date', 'time'), help='observation time in format: 03/10/2015 15:23:10.000001')
    parser.add_argument('--duration', nargs=1, metavar="num", help='offset to observeTime in second')
    parser.add_argument('--update_interval', nargs=1, metavar="num", help='interval to re-calculate overstep in second')
    parser.add_argument('--overlap', nargs=1, metavar="ratio", help='overlap point between beams')
    parser.add_argument('--retile', nargs=1, metavar="second", help='retile interval, default: 0(no retile)')
    parser.add_argument('--overlay_source', nargs=1, metavar="file", help='extra overlay sources')
    parser.add_argument('--retile_threshold', nargs=1, metavar="ratio", help='retile threshold, default: 0.5')
    parser.add_argument('--hour_angle',action='store_true', help='show hour angle in the x axis instead of time, default: false')
    parser.add_argument('--retile_rate',action='store_true', help='show rate of retile occurence instead of overlap rate, default: false')
    parser.add_argument('--beamnum', nargs=1, metavar="num", help='beam number of tiling')
    parser.add_argument('--subarray', nargs='+', metavar="num", help='list of antennas, saperated by comma')
    parser.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true")



    args = parser.parse_args()

    interpolation = True
    frequencies = [1.4e9,]
    paras = None
    plotting = False
    resolution = 10 #in arcsecond
    size = 20
    zoom = 1
    duration = 1
    frame = 'RADEC'
    overlap = 0.5
    beamnum = 400
    update_interval = 3600
    retile = 0
    hourAngle = False
    retile_rate = False
    overlay_source = None
    retile_threshold = 0.5
    if args.plot == True:
        plotting = True
        if args.inte == True:
            interpolation = True

    if args.hour_angle == True:
        hourAngle = True

    if args.retile_rate == True:
        retile_rate = True

    if args.ants is not None:
        with open(args.ants[0], 'r') as antFile:
            antennaCoords = antFile.readlines()
    else:
        parser.error("no antennas file, try --ants file")

    if args.overlay_source is not None:
        overlay_source = np.genfromtxt(args.overlay_source[0], dtype=None)
        """
        structured array with one row can not be itered.
        """
        if overlay_source.size == 1:
            overlay_source = [overlay_source.tolist(),]
    else:
        overlay_source = None


    if args.datetime is not None:
        observeTime=datetime.datetime.strptime(args.datetime[0] + " "
                + args.datetime[1], '%Y.%m.%d %H:%M:%S.%f')
    else:
        parser.error("no time specified, try --datetime date time")

    if args.subarray is not None:
        arrayString = "".join(args.subarray)
        subarray = arrayString.split(",")
    else:
        subarray = []


    if args.source is not None:
        sourceCoord = args.source
    else:
        parser.error("no source specifed, try --source RA DEC")

    if args.frame is not None:
        frame = args.frame[0].upper()
        if frame != 'RADEC' and frame != 'AZIALT':
            parser.error("frame not recongnized, should be RADEC or AziAlt")
    else:
        # logging.warning("frame not specified, default to RADEC")
        frame = 'RADEC'

    if args.duration is not None:
        duration = int(args.duration[0])
    if args.update_interval is not None:
        update_interval = int(args.update_interval[0])
    if args.retile is not None:
        retile = int(args.retile[0])
    if args.retile_threshold is not None:
        retile_threshold = float(args.retile_threshold[0])

    if args.beamnum is not None:
        beamnum = int(args.beamnum[0])

    if args.overlap is not None:
        overlap = float(args.overlap[0])
    else:
        overlap = 0.5

    if args.resolution is not None:
        resolution = float(args.resolution[0])
    if args.size is not None:
        size = int(args.size[0])
    if args.freqrange is None and args.freq is None:
        parser.error("no freqencies or frequency range specified")
    elif args.freq is not None:
        frequencies = [float(freq) for freq in args.freq]
    elif args.freqrange is not None:
        frequencies = np.arange(float(args.freqrange[0]),
                float(args.freqrange[1]), float(args.freqrange[2]))

    if args.verbose:
        logger.setLevel(logging.INFO)


    # paras = antennaCoords, sourceCoord, frame, observeTime, resolution, size, zoom
    paras  = {"antennaCoords": antennaCoords,
        "sourceCoord": sourceCoord,
        "observeTime": observeTime,
        "frequencies":frequencies,
        "duration":duration,
        "overlap":overlap,
        "beamNum":beamnum,
        "subarray":subarray,
        "update_interval":update_interval,
        "retile":retile,
        "retileRate":retile_rate,
        "hourAngle":hourAngle,
        "retile_threshold":retile_threshold,
        "overlay_source":overlay_source,
        "plotting":plotting}

    creatBeamMatrix(**paras)

def captureNegetiveNumber():
    for i, arg in enumerate(sys.argv):
          if (arg[0] == '-') and arg[1].isdigit(): sys.argv[i] = ' ' + arg

def main():
    captureNegetiveNumber()
    parser = argparse.ArgumentParser()
    parseOptions(parser)

if __name__ == "__main__":
    main()

'''
./maketiling.py --ants antenna.csv --freq 1.284e9 --source 00:24:05.67 -72:04:52.60 --overlay_source overlay_source_file --datetime 2020.05.02 06:02:13.663903 --beamnum 263 --verbose --overlap 0.7 --subarray 000, 001, 002, 003, 004, 006, 007, 008, 009, 010, 011, 012, 013, 014, 015, 016, 017, 018, 019, 020, 021, 022, 024, 026, 027, 028, 029, 030, 031, 033, 034, 035, 036, 038, 039, 040, 041, 042, 044, 045, 046, 047, 048, 049, 050, 051, 052, 053, 054, 055, 056, 057, 058, 059, 060, 061 --duration 7200 --update_interval 1800 --plot

the format of the overlay_source_file is:
source_name_1 ra dec
source_name_2 ra dec
.
.
.
'''
