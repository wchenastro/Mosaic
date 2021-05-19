#!/usr/bin/env python

import numpy as np
import sys, datetime
import argparse
import logging

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
loggerFormat = '%(asctime)-15s %(filename)s  %(message)s'
logging.basicConfig(format = loggerFormat, level=logging.WARNING)
logger = logging.getLogger()

import katpoint
from mosaic.beamforming import PsfSim, generate_nbeams_tiling
from mosaic.coordinate import convertBoresightToDegree, createTilingRegion
from mosaic.plot import plot_overlap, plot_interferometry



def makeKatPointAntenna(antennaString):
    antennaKat = []

    for antenna in antennaString:
        antkat = katpoint.Antenna(antenna)
        antennaKat.append(antkat)

    return antennaKat

def creatBeamMatrix(antennaCoords, sourceCoord, observeTime, frequencies, beamNum,
        duration, overlap, subarray, update_interval, overlay_source, size, resolution,
        tilingMethod, tilingShape, tilingParameter, tilingParameterCoordinateType):

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

    counter = []
    step = update_interval/60.
    hour = duration/3600.
    grids = []
    horizons = []
    retileTime = []
    startTime = observeTime
    localHourAngle = []
    retileStamp = []
    for minute in np.arange(0, hour*60+1, step):
        offsetTime = startTime + datetime.timedelta(minutes=minute)
        newBeamShape = psf.get_beam_shape(sourceCoord, offsetTime, size, resolution)
        newBeamShape.plot_psf("plots/beamshape/beamshape{:03d}.png".format(int(minute)),
            overlap = overlap, shape_overlay=True)
        # newBeamShape.psf.write_fits(
                # "plots/beamshape/beamshapebeamshape{:03d}.fits".format(int(minute)))

        tiling = generate_nbeams_tiling(newBeamShape, beamNum, overlap, tilingMethod,
                tilingShape, tilingParameter, tilingParameterCoordinateType)
        tiling.plot_tiling("plots/tiling/tiling%03d.png" % minute,
            HD=True, edge=True)
        equatorial_coordinates = tiling.get_equatorial_coordinates()
        # np.savetxt("coordinate", equatorial_coordinates)
        actualShape = tiling.meta["axis"]
        createTilingRegion(equatorial_coordinates, actualShape,
                "plots/tiling/tiling%03d.reg" % minute)



def readPolygonRegion(filename):
    with open(filename, 'r') as regionFile:
        polygonLine = regionFile.readlines()[3]
        pointString = polygonLine.split(")")[0].split("(")[1]
        points = np.array(pointString.split(",")).astype(float)

        return points.reshape((-1, 2))

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
    parser.add_argument('--overlay_source', nargs=1, metavar="file", help='extra overlay sources')
    parser.add_argument('--beamnum', nargs=1, metavar="num", help='beam number of tiling')
    parser.add_argument('--subarray', nargs='+', metavar="num", help='list of antennas, saperated by comma')
    parser.add_argument('--tiling_method', nargs=1, metavar="method",
            help='tiling method, such as \"variable_overlap\" or \"variable_size\".')
    parser.add_argument('--tiling_shape', nargs=1, metavar="shape",
            help='shape of the tiling boundary, such as \"circle\", \"hexagon\", \"ellipse\", \"polygon\", \"annulus\".')
    parser.add_argument('--tiling_parameter', nargs='+', metavar='parameters', help='parameters for the tiling')
    parser.add_argument('--tiling_parameter_file', nargs=1, metavar='parameter_file', help='parameter_file for the tiling')
    parser.add_argument('--tiling_parameter_coordinate_type', nargs=1, metavar='coordinate_type',
            help='type of the coordinate of the tiling parameter, such as \"equatorial\", default is image(pixel) coordinate')
    parser.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true")



    args = parser.parse_args()

    interpolation = True
    frequencies = [1.4e9,]
    paras = None
    plotting = False
    resolution = None #in arcsecond
    size = 400
    zoom = 1
    duration = 0
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

    if args.tiling_shape is not None:
        tilingShape = args.tiling_shape[0]
    else:
        tilingShape = "circle"

    if args.tiling_parameter_coordinate_type is not None:
        tilingParameterCoordinateType = args.tiling_parameter_coordinate_type[0]
    else:
        tilingParameterCoordinateType  = 'equatorial'

    if args.tiling_method is not None:
        tilingMethod = args.tiling_method[0]
        if args.tiling_method[0] == "variable_overlap":
            if args.tiling_parameter is None and tilingShape != "polygon":
                parser.error("no parameter for \"variable_overlap\" method!")
                exit(-1)
            if tilingShape == "circle":
                tilingParameter = float(args.tiling_parameter[0])
            elif tilingShape == "hexagon":
                tilingParameter = [float(args.tiling_parameter[0]),
                                    float(args.tiling_parameter[1])]
            elif tilingShape == "ellipse":
                tilingParameter = [float(args.tiling_parameter[0]),
                                    float(args.tiling_parameter[1]),
                                    float(args.tiling_parameter[2])]
            elif tilingShape == "polygon":
                if args.tiling_parameter is not None:
                    parameterString = "".join(args.tiling_parameter).split(",")
                    tilingParameter = [float(coord) for coord in parameterString]
                    tilingParameter = np.array(tilingParameter).reshape(-1,2).tolist()
                elif args.tiling_parameter_file is not None:
                    tilingParameter = readPolygonRegion(
                            args.tiling_parameter_file[0]).tolist()
                else:
                    parser.error("no parameter for polygon tiling!")
                    exit(-1)
        else:
            tilingParameter = None
    else:
        tiling_method = "variable_size"
        tilingParameter = None

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
        "size":size,
        "resolution":resolution,
        "tilingMethod":tilingMethod,
        "tilingShape":tilingShape,
        "tilingParameter":tilingParameter,
        "tilingParameterCoordinateType":tilingParameterCoordinateType,
        "overlay_source":overlay_source}

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
