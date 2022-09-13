#!/usr/bin/env python

import numpy as np
import sys, datetime
import argparse
import logging

import matplotlib
# matplotlib.use('Agg')

import matplotlib.pyplot as plt
loggerFormat = '%(asctime)-15s %(filename)s  %(message)s'
logging.basicConfig(format = loggerFormat, level=logging.WARNING)
logger = logging.getLogger()

import katpoint
from mosaic.beamforming import PsfSim, generate_nbeams_tiling
from mosaic.coordinate import createTilingRegion, readPolygonRegion, convert_sexagesimal_to_degree, convert_equatorial_coordinate_to_pixel
from mosaic import __version__



def makeKatPointAntenna(antennaString):
    antennaKat = []

    for antenna in antennaString:
        antkat = katpoint.Antenna(antenna)
        antennaKat.append(antkat)

    return antennaKat

def createBeamMatrix(antennaCoords, sourceCoord, observeTime, frequencies, beamNum,
        duration, overlap, subarray, update_interval, overlay_source, overlay_source_name,
        size, resolution, tilingMethod, tilingShape, tilingParameter, tilingParameterCoordinateType,
        weights, interpolation, output):

    if subarray != []:
        antennaKat = makeKatPointAntenna(
                [antennaCoords[ant] for ant in subarray])
    else:
        antennaKat = makeKatPointAntenna(antennaCoords)


    # boresight = sourceCoord
    boresight = katpoint.Target('boresight, radec, {}, {}'.format(
                sourceCoord[0], sourceCoord[1]))
    reference = makeKatPointAntenna(
            ["ref, -30:42:39.8, 21:26:38.0, 1035.0",])[0]
    psf = PsfSim(antennaKat, frequencies[0])

    newBeamShape = psf.get_beam_shape(sourceCoord, observeTime, size, resolution, weights)
    if "psf_plot" in output:
        newBeamShape.plot_psf(output["psf_plot"][0], overlap = overlap,
                shape_overlay=True, interpolation=interpolation)
    if "psf_fits" in output:
        newBeamShape.psf.write_fits(output["psf_fits"][0])
    if ("tiling_plot" in output) or ("tiling_coordinate" in output) or ("tiling_region" in output):

        tiling = generate_nbeams_tiling(newBeamShape, beamNum, overlap, tilingMethod,
                tilingShape, tilingParameter, tilingParameterCoordinateType)
        if "tiling_plot" in output:
            tiling.plot_tiling(output["tiling_plot"][0], HD=True, edge=True,
                    extra_coordinates = overlay_source,
                    extra_coordinates_text = overlay_source_name)
        if "tiling_coordinate" in output:
            equatorial_coordinates = tiling.get_equatorial_coordinates()
            np.savetxt(output["tiling_coordinate"][0], equatorial_coordinates)
        if "tiling_region" in output:
            equatorial_coordinates = tiling.get_equatorial_coordinates()
            actualShape = tiling.meta["axis"]
            createTilingRegion(equatorial_coordinates, actualShape, output["tiling_region"])

def parseOptions(parser):
    # enable interpolation when ploting
    parser.add_argument('--no_interpolation', action='store_true', help=argparse.SUPPRESS)
    # plot the beam images
    parser.add_argument('--psf_plot', nargs='+', metavar="file [paras]", help='filename for the psf plot')
    parser.add_argument('--psf_fits', nargs='+', metavar="file [paras]", help='name for the psf fits file')
    parser.add_argument('--tiling_plot', nargs='+', metavar="file [paras]", help='filename for the tiling plot')
    parser.add_argument('--tiling_coordinate', nargs='+', metavar="file [paras]", help='filename for the tiling coordinate')
    parser.add_argument('--tiling_region', nargs=1, metavar="file", help='filename for the tiling region')
    parser.add_argument('--ants', nargs=1, metavar="file", help='antenna coodinates files')
    parser.add_argument('--resolution', nargs=1, metavar="asec", help='resolution in arcsecond')
    parser.add_argument('--size', nargs=1, metavar="num", help='width in pixels')
    parser.add_argument('--freqrange', nargs=3, metavar=('s', 't', 'i'), help='freqencies range as start stop interval')
    parser.add_argument('--freq', nargs='+',  help='multiple freqencies')
    parser.add_argument('--frame', nargs=1, metavar="RADEC/AziAlt", help='source coordinate frame')
    parser.add_argument('--source', nargs=2, metavar=('RA/Azi', 'DEC/Alt'), help='source position in RADEC or AziAlt')
    parser.add_argument('--datetime', nargs=2, metavar=('date', 'time'), help='observation time in format: 03/10/2015 15:23:10.000001')
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
    parser.add_argument("--version", help="show the version of this package", action="store_true")
    parser.add_argument("--weight", action="store_true",
            help='apply weights to individual antenna, attach weight after the item in --subarray, e.g., 0:0.5, 1:0.7, 2:0.5 ')




    args = parser.parse_args()

    if args.verbose:
        logger.setLevel(logging.INFO)

    if args.version is True:
        print(__version__)
        exit(0)
    else:
        logger.info("Mosaic " + __version__)

    if args.weight is True:
        weights = []
    else:
        weights = None

    frequencies = [1.4e9,]
    paras = None
    resolution = None #in arcsecond
    size = 400
    duration = 0
    overlap = 0.5
    beamnum = 400
    update_interval = 3600
    overlay_source = None
    output = {}
    if args.psf_plot is not None:
        output["psf_plot"] = args.psf_plot
    if args.psf_fits is not None:
        output["psf_fits"] = args.psf_fits
    if args.tiling_plot is not None:
        output["tiling_plot"] = args.tiling_plot
    if args.tiling_coordinate is not None:
        output["tiling_coordinate"] = args.tiling_coordinate
    if args.tiling_region is not None:
        output["tiling_region"] = args.tiling_region[0]


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
        subarray = []
        arrayString = "".join(args.subarray)
        ant_weights = arrayString.split(",")
        for ant_weight in ant_weights:
            ant_weight_pair = ant_weight.split(':')
            if weights is not None:
                subarray.append(int(ant_weight_pair[0]))
                if len(ant_weight_pair) > 1:
                    complex_weight = complex(ant_weight_pair[1])
                    if complex_weight.imag == 0:
                        weights.append(float(ant_weight_pair[1]))
                    else:
                        weights.append(complex_weight)
                else:
                    weights.append(1.0)
            else:
                subarray.append(int(ant_weight_pair[0]))
    else:
        subarray = []


    if args.source is not None:
        sourceCoord = args.source
    else:
        parser.error("no source specifed, try --source RA DEC")

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
        tilingMethod = "variable_size"
        tilingParameter = None

    if args.no_interpolation:
        interpolation = False
    else:
        interpolation = True

    if args.overlay_source is not None:
        overlay_coords = np.genfromtxt(args.overlay_source[0], dtype=None)
        if len(overlay_coords.shape) == 1:
            overlay_coords = overlay_coords.reshape(1, -1)
        bore_sight = convert_sexagesimal_to_degree([sourceCoord,])[0]
        overlay_source_degree = convert_sexagesimal_to_degree(overlay_coords[:, 1:])
        overlay_coordinates = convert_equatorial_coordinate_to_pixel(
                overlay_source_degree, bore_sight)

        overlay_source = overlay_coordinates
        overlay_source_name = overlay_coords[:, 0].astype('str')

    else:
        overlay_source = []
        overlay_source_name = []


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
        "overlay_source":overlay_source,
        "overlay_source_name":overlay_source_name,
        "weights":weights,
        "interpolation":interpolation,
        "output":output}

    createBeamMatrix(**paras)

def captureNegetiveNumber():
    for i, arg in enumerate(sys.argv):
          if (arg[0] == '-') and arg[1].isdigit(): sys.argv[i] = ' ' + arg

def main():
    captureNegetiveNumber()
    parser = argparse.ArgumentParser()
    parseOptions(parser)

if __name__ == "__main__":
    main()
