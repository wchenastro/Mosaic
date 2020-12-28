#!/usr/bin/env python

import numpy as np
import sys, datetime
import argparse
import logging

loggerFormat = '%(asctime)-15s %(filename)s  %(message)s'
logging.basicConfig(format = loggerFormat, level=logging.WARNING)
logger = logging.getLogger()

import katpoint
from mosaic import PsfSim, generate_nbeams_tiling, generate_radius_tiling
from mosaic.coordinate import convertBoresightToDegree
from mosaic.plot import plot_overlap, plot_interferometry, plot_all

import matplotlib.pyplot as plt

def makeKatPointAntenna(antennaString):
    antennaKat = []

    for antenna in antennaString:
        antkat = katpoint.Antenna(antenna)
        antennaKat.append(antkat)

    return antennaKat

def creatBeamMatrix(antennaCoords, sourceCoord, observeTime, frequencies,
        duration, overlap, beamNum, subarray, update_interval, hourAngle,
        retile, retileRate, retile_threshold):

    if subarray != []:
        antennaKat = makeKatPointAntenna(
                [antennaCoords[int(ant)] for ant in subarray])
    else:
        antennaKat = makeKatPointAntenna(antennaCoords)

    boresight = katpoint.Target('boresight, radec, {}, {}'.format(
                sourceCoord[0], sourceCoord[1]))
    reference = makeKatPointAntenna(
            ["ref, -30:42:39.8, 21:26:38.0, 1035.0",])[0]
    psf = PsfSim(antennaKat, frequencies[0])
    beamShape = psf.get_beam_shape(boresight, observeTime, beam_number = 10000, beam_size = 1)



    tiling = generate_nbeams_tiling(
            beamShape, beamNum, overlap = overlap)
    tiling_coordinats = tiling.get_equatorial_coordinates()

    counter = []
    step = update_interval/60.
    hour = duration/3600.
    grids = []
    horizons = []
    retileTime = []
    startTime = observeTime
    localHourAngle = []


    antennasGeo = np.array([antenna.geo for antenna in beamShape.antennas])
    referenceGeo = beamShape.reference_antenna.geo

    for minute in np.arange(0, hour*60+1, step):
        offsetTime = startTime + datetime.timedelta(minutes=minute)
        newBeamShape = psf.get_beam_shape(sourceCoord, offsetTime, beam_number = 10000, beam_size = 1)
        overlap = tiling.calculate_overlap("counter", newBeamShape)

        normalizedCounts = overlap.calculate_fractions()

        if retile != 0 and normalizedCounts[1] < retile_threshold:
            tiling = generate_nbeams_tiling(newBeamShape, 400, overlap = 0.5)
            overlap = tiling.calculate_overlap("counter", newBeamShape)

        tiling.beam_shape = newBeamShape
        tiling.plot_tiling("plots/overlap/tiling%03d.png" % minute)

        counter.append(normalizedCounts)
        grids.append(overlap.metrics)
        horizons.append(newBeamShape.horizon)

        inteferometry_paras = [antennasGeo , referenceGeo, newBeamShape.horizon]
        beamshap_paras = [newBeamShape.psf.image,
                          newBeamShape.psf.bore_sight.equatorial,
                          newBeamShape.psf.image_range,
                          newBeamShape.axisH/newBeamShape.resolution,
                          newBeamShape.axisV/newBeamShape.resolution,
                          newBeamShape.angle,
                          None, True]
        overlap_paras = [overlap.metrics, "counter", None, "Trend"]
        retile_paras = [counter, np.arange(len(counter))*step, hour*60 +1, retile_threshold]

        plot_all(inteferometry_paras, beamshap_paras,
                overlap_paras, retile_paras,
                "plots/overlap/overview%03d.png" % minute)

    counter = np.array(counter)
    x = np.arange(len(counter))*step


def parseOptions(parser):
    # enable interpolation when ploting
    parser.add_argument('--inte', action='store_true', help=argparse.SUPPRESS)
    # plot the beam images
    parser.add_argument('--plot', action='store_true', help=argparse.SUPPRESS)
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
    duration = 0
    frame = 'RADEC'
    overlap = 0.5
    beamnum = 400
    update_interval = 10
    retile = 0
    hourAngle = False
    retile_rate = False
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
        "retile_threshold":retile_threshold}

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
./overview.py --ants antenna.csv --freq 1.284e9 --source 1:18:1.704509 -30:43:20.1839 --datetime 2017.03.19 7:03:21.000 --beamnum  400 --verbose --overlap 0.5 --subarray 000,001, 002, 003, 004, 005, 006, 007, 008, 009, 010, 011, 012, 013, 014, 015, 016, 017, 018, 019, 020, 021, 022, 023, 024, 025, 026, 027, 028, 029, 030, 031, 032, 034, 035, 036, 037, 038, 039, 040, 041, 042, 043, 047 --duration 43200 --update_interval 3600
'''
