#!/usr/bin/env python

import numpy as np
import datetime
import argparse
import h5py

from astropy import wcs
from astropy.io import fits

import coordinate as coord
from interferometer import InterferometryObservation


class h5Writer:

    fileObj = None
    dataset = None
    dataIdx = 0

    def __init__(self, fileName):
        self.fileObj = h5py.File(fileName, 'w')
        self.header = self.fileObj.create_group('header')
        self.dataset = self.fileObj.create_group('dataset')
        self.dataIdx = 0


    def writeBeams(self, attr, data):
        newData = self.dataset.create_dataset(str(self.dataIdx), data=data)
        newData.attrs['waveLength'] = attr
        self.dataIdx += 1

    def close(self):
        self.fileObj.close()

def writeToFits(WCS, data):
    w = wcs.WCS(naxis=2)

    # "Airy's zenithal" projection
    w.wcs.crpix = WCS['crpix']
    w.wcs.cdelt = WCS['cdelt']
    w.wcs.crval = WCS['crval']
    w.wcs.ctype = WCS['ctype']
    # w.wcs.set_pv([(2, 1, 45.0)])

    header = w.to_header()

    hdu = fits.PrimaryHDU(header=header, data=data)
    # Save to FITS file
    hdu.writeto('psfsim.fits', overwrite=True)

def creatBeamMatrix(paras, freqencies, plotting = False, interpolation=True):

    antennaCoords, sourceCoord, frame, observeTime, resolution, size, zoom = paras

    # the values provided by http://public.ska.ac.za/meerkat
    arrayRefereceGEODET = (21.44389, -30.71317, 0)
    '''speed Of Light'''
    sol = 299792458
    waveLengths = sol/freqencies

    defaultBeamSizeFactor = 1
    defaultBeamNumber = 400
    defaultBoreSight = (21.411, -30.721)

    observation = InterferometryObservation(arrayRefereceGEODET,
            None, waveLengths)
    observation.setBoreSight(defaultBoreSight)
    observation.setBeamSizeFactor(zoom)
    observation.setGridNumber(size)
    observation.setBeamNumber(size*size)
    observation.setResolution(resolution)
    observation.setInterpolating(interpolation)
    observation.setAutoZoom(False)

    if frame == "RADEC":
        observation.setBoreSight(sourceCoord)
        observation.setInputType(InterferometryObservation.equatorialInput)
    elif frame == "AZIALT":
        observation.setHorizontal(sourceCoord)
        observation.setInputType(InterferometryObservation.horizontalInput)

    observation.setObserveTime(observeTime)

    saveParas(paras)

    # writer = h5Writer('beams.hdf5')
    images = observation.createPSF(antennaCoords, waveLengths, None, plotting)
    # writer.close()
    wcs = observation.getWCS()
    writeToFits(wcs, images)

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
    parser.add_argument('--freq', nargs=3, metavar=('s', 't', 'i'), required=True, help='freqencies range as start stop interval')
    parser.add_argument('--frame', nargs=1, metavar="RADEC/AziAlt", help='source coordinate frame')
    parser.add_argument('--source', nargs=2, metavar=('RA/Azi', 'DEC/Alt'), help='source position in RADEC or AziAlt')
    parser.add_argument('--datetime', nargs=2, metavar=('date', 'time'), help='observation time in format: 03/10/2015 15:23:10.000001')

    args = parser.parse_args()

    interpolation = False
    paras = None
    plotting = False
    resolution = 10 #in arcsecond
    size = 20
    zoom = 1
    frame = 'RADEC'
    if args.plot == True:
        plotting = True
        if args.inte == True:
            interpolation = True
    if args.last == True:
        paras = loadParas()
    else:
        if args.ants is not None:
            antennaCoords = np.loadtxt(args.ants[0], delimiter=',')
        else:
            parser.error("no antennas file, try --ants file")

        if args.datetime is not None:
            observeTime=datetime.datetime.strptime(args.datetime[0]+args.datetime[1], '%d/%m/%Y%H:%M:%S.%f')
        else:
            parser.error("no time specified, try --datetime date time")

        if args.source is not None:
            sourceCoord = (float(args.source[0]), float(args.source[1]))
        else:
            parser.error("no source specifed, try --source RA DEC")

        if args.frame is not None:
            frame = args.frame[0].upper()
            if frame != 'RADEC' and frame != 'AZIALT':
                parser.error("frame not recongnized, should be RADEC or AziAlt")
        else:
            print("frame not specified, default to RADEC")
            frame = 'RADEC'

        if args.resolution is not None:
            resolution = float(args.resolution[0])
        if args.size is not None:
            size = int(args.size[0])
        if args.zoom is not None:
            zoom = int(args.zoom[0])

        paras = antennaCoords, sourceCoord, frame, observeTime, resolution, size, zoom

    freqencies = np.arange(float(args.freq[0]), float(args.freq[1]), float(args.freq[2]))
    creatBeamMatrix(paras, freqencies, plotting, interpolation)

def saveParas(paras):
    with open('bsParas', 'wb') as paraFile:
        np.save(paraFile, paras)

def loadParas():
    with open('bsParas', 'rb') as paraFile:
        paras = np.load(paraFile)
        return paras


def main():
    parser = argparse.ArgumentParser()
    parseOptions(parser)

if __name__ == "__main__":
    main()
