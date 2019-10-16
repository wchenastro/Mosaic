#!/usr/bin/env python

import numpy as np
import sys, datetime
import argparse
import logging

loggerFormat = '%(asctime)-15s %(filename)s  %(message)s'
logging.basicConfig(format = loggerFormat, level=logging.WARNING)
logger = logging.getLogger()

import katpoint
from mosaic import PsfSim, DelayPolynomial, generate_nbeams_tiling, generate_radius_tiling
from mosaic.coordinate import convertBoresightToDegree


def makeKatPointAntenna(antennaString):
    antennaKat = []

    for antenna in antennaString:
        antkat = katpoint.Antenna(antenna)
        antennaKat.append(antkat)

    return antennaKat

def creatBeamMatrix(antennaCoords, sourceCoord, observeTime,
        frequencies, duration, overlap, beamNum, subarray):

    if subarray != []:
        antennaKat = makeKatPointAntenna([antennaCoords[int(ant)] for ant in subarray])
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
    # beamShape.plot_psf("beam.png", shape_overlay=False)
    # beamShape.psf.write_fits('psf.fits')
    # beamShape.plot_interferometry("interferometry.png")


    # margin = max(int(beamNum * 0.25), 16)
    tiling = generate_nbeams_tiling(
            beamShape, beamNum, overlap = overlap)
    tiling_coordinats = tiling.get_equatorial_coordinates()
    tiling.plot_tiling("tiling.svg", index = True)
    np.savetxt("tilingCoord", tiling_coordinats)
    np.savetxt("tilingCoord_pixel", tiling.coordinates)

    # targets = DelayPolynomial.make_katpoint_target(tiling_coordinats)
    # delay = DelayPolynomial(antennaKat, beamShape.bore_sight, targets, reference)
    # polynomials = delay.get_delay_polynomials(observeTime, duration)
    # np.save("polynomials", polynomials)


def parseOptions(parser):
    # enable interpolation when ploting
    parser.add_argument('--inte', action='store_true', help=argparse.SUPPRESS)
    # plot the beam images
    parser.add_argument('--plot', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--last', action='store_true', help='use last parameters')
    parser.add_argument('--ants', nargs=1, metavar="file", help='antenna coodinates files')
    parser.add_argument('--resolution', nargs=1, metavar="asec", help='resolution in arcsecond')
    parser.add_argument('--size', nargs=1, metavar="num", help='width in pixels')
    parser.add_argument('--zoom', nargs=1, metavar="num", help=argparse.SUPPRESS)
    parser.add_argument('--freqrange', nargs=3, metavar=('s', 't', 'i'), help='freqencies range as start stop interval')
    parser.add_argument('--freq', nargs='+',  help='multiple freqencies')
    parser.add_argument('--frame', nargs=1, metavar="RADEC/AziAlt", help='source coordinate frame')
    parser.add_argument('--source', nargs=2, metavar=('RA/Azi', 'DEC/Alt'), help='source position in RADEC or AziAlt')
    parser.add_argument('--datetime', nargs=2, metavar=('date', 'time'), help='observation time in format: 03/10/2015 15:23:10.000001')
    parser.add_argument('--duration', nargs=1, metavar="num", help='minutes offset to observeTime')
    parser.add_argument('--overlap', nargs=1, metavar="ratio", help='overlap point between beams')
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
    if args.plot == True:
        plotting = True
        if args.inte == True:
            interpolation = True

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
    if args.zoom is not None:
        zoom = int(args.zoom[0])

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
        "subarray":subarray}

    creatBeamMatrix(**paras)

def saveParas(paras):
    with open('bsParas', 'wb') as paraFile:
        np.save(paraFile, paras)

def loadParas():
    with open('bsParas', 'rb') as paraFile:
        paras = np.load(paraFile)
        return paras

def captureNegetiveNumber():
    for i, arg in enumerate(sys.argv):
          if (arg[0] == '-') and arg[1].isdigit(): sys.argv[i] = ' ' + arg


def main():
    captureNegetiveNumber()
    parser = argparse.ArgumentParser()
    parseOptions(parser)

if __name__ == "__main__":
    main()
