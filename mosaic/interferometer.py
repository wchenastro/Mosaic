#! /usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

import mosaic.coordinate as coord
from mosaic.plot import plotBeamContour
from mosaic.utilities import normSigma, normInverse
from mosaic.beamshape import calculateBeamSize, createBeamshapeModel

import inspect, pickle, datetime, logging

logger = logging.getLogger(__name__)

class PointSpreadFunction(object):
    """
    class for point spread function
    """

    def __init__(self, image, bore_sight, width, wcs_header, image_range):
        self.image = image
        self.bore_sight = bore_sight
        self.width = width
        self.wcs_header = wcs_header
        self.image_range = image_range
        self.beamModel = None

    def write_fits(self, filename):
        coord.writeFits(self.wcs_header, self.image, filename)

class InterferometryObservation:

    def __init__(self, arrayReferece, waveLength):
        self.arrayReferece = arrayReferece
        self.waveLength = waveLength
        self.beamSizeFactor = 1
        self.beamNumber = 400
        self.interpolating = True
        self.autoZoom = True
        self.gridNumOfDFT = 100000.0
        self.imageDensity = 20
        self.resolution = 1/3600.0
        self.weights = None
        self.boresightFrame = coord.Boresight.EquatorialFrame
        self.updateBeamCoordinates(self.beamSizeFactor, self.imageDensity, self.gridNumOfDFT)


    def setInterpolating(self, state):
        if state == True:
            self.interpolating = True
        else:
            self.interpolating = False

    #def setInputType(self, inputType):
    #   self.inputType = inputType

    def setObserveTime(self, dateTime):
        if type(dateTime) != datetime.datetime:
            dateTime = coord.epochToDatetime(dateTime)
        self.observeTime = dateTime

    def setAntennaWeights(self, weights):
        self.antennaWeights = weights

    def getObserveTime(self):
        return self.observeTime

    def getPointSpreadFunction(self):
        return self.psf

    def setBeamSizeFactor(self, size):
        if size != self.beamSizeFactor:
            isUpdated = self.updateBeamCoordinates(size, self.imageDensity, self.gridNumOfDFT)
            if isUpdated:
                self.beamSizeFactor = size
                return True
            else:
                logger.warning("failed to set resolution")
                return False

    def setAutoZoom(self, switch):
        self.autoZoom = switch

    def getBeamSizeFactor(self):
        return self.beamSizeFactor

    def getBeamNumber(self):
        return self.beamNumber

    def setGridNumber(self, number):
        self.gridNumOfDFT = number


    def setBeamNumber(self, number):
        if number != self.beamNumber:
            oldBeamNumber = self.beamNumber
            self.beamNumber = number
            density = int(np.sqrt(number))
            if density % 2 != 0: density += 1
            oldDensity = self.imageDensity
            self.setImageDensity(density)
            isUpdated = self.updateBeamCoordinates(self.beamSizeFactor, density, self.gridNumOfDFT)
            if isUpdated:
                return True
            else:
                logger.warning("failed to set image density")
                self.beamNumber = oldBeamNumber
                self.setImageDensity(oldDensity)

    def setResolution(self, resolution):
        '''resolution default in arc second deg for now'''
        self.resolution = resolution/3600.0

    def getResolution(self):
        '''in degree'''
        return self.resolution * self.beamSizeFactor

    def getBaselines(self):
        return self.baselines

    def getAmplitude(self):
        return self.amplitude

    def getBeamCoordinates(self):
        return self.beamCoordinates

    def getBaselinesNumber(self):
        return len(self.baselines)

    def saveParas(self, fileName):
        coordinates = self.getAntCoordinates()
        observeTime = self.getObserveTime()
        source = self.getBoreSight()

        with open(fileName, 'wb') as paraFile:
            pickle.dump([coordinates, source, observeTime], paraFile)

    def getsynthesizedBeam(self):
        return self.beamSynthesized

    def setBoreSight(self, boresight=None, frame=None):
        if boresight is not None:
            self.boresightInput = boresight
            if type(boresight[0]) != type("str"):
                self.boresightCoord = boresight
            elif len(boresight[0].split(":")) > 1:
                self.boresightCoord  = coord.convertBoresightToDegree(boresight)
            else:
                self.boresightCoord = (float(boresight[0]), float(boresight[1]))
        if frame is not None:
            self.boresightFrame = frame

    def getBoreSight(self):
        return self.boresight

    #def getHorizontal(self):
    #   return self.boreSight.horizontal

    #def setHorizontal(self, horizontal):
    #   self.boresightInput = np.deg2rad(horizontal)

    def getProjectedBaselines(self):
        #return self.baselines
        return self.projectedBaselines
        #uvCoord = self.projectedBaselines/self.waveLength
        #return np.concatenate((uvCoord, -uvCoord))

    def getBeamAxis(self):
        beamshapeModel = self.fitContour()
        # get the 50% overlap beam shape
        bottomOverlap = beamshapeModel[0, 0]
        topOverlap = beamshapeModel[-1, 0]
        index = (0.5 - bottomOverlap) / (topOverlap-bottomOverlap) * (beamshapeModel.shape[0] - 1)
        axis1, axis2, angle = beamshapeModel[int(np.round(index))][1:]
        self.beamAxis[0:3] = [axis1, axis2, angle]

        width1, width2 = coord.convert_pixel_length_to_equatorial(axis1, axis2,
                angle, self.boresight.equatorial)

        # self.beamSize = [width1.degree, width2.degree]

        # logger.info("axis1: {:.3g}, axis2: {:.3g}, angle: {:.3f} in pixel plane"
                # .format(axis1, axis2, angle))

        logger.info("beamshape: width1: {:.5g} arcsec, width2: {:.5g} arcsec, angle: {:.5g} degree"
                .format(width1.arcsecond, width2.arcsecond, angle))
        """
        logger.info("beamshape: width1: {:.13g} arcsec, width2: {:.13g}, angle: {:.13g} arcsec in equatorial plane"
                .format(width1.degree, width2.degree, angle))
        """

        self.beamshapeModel = beamshapeModel
        return self.beamAxis

    def getBeamshapeModel(self):
        if self.beamshapeModel is None:
            self.beamshapeModel = self.fitContour()
        return self.beamshapeModel

    def getImageLength(self):
        return self.imageLength

    def getImageDensity(self):
        return self.imageDensity

    def setImageDensity(self, density):
        self.imageDensity = density

    def getImageData(self):
        return self.imageData

    def getAntCoordinates(self):
        return [antenna.geo for antenna in self.array.getAntennas]

    def createBaselines(self, antCoordinatesENU):
        baselines = []
        index = 1
        for antenna1 in antCoordinatesENU:
            for antenna2 in antCoordinatesENU[index:]:
                baselines.append([antenna1[0]-antenna2[0],
                        antenna1[1]-antenna2[1],antenna1[2]-antenna2[2]])
            index += 1

        return baselines

    def updateBeamCoordinates(self, interval, imageDensity, DFTSideLength):
        # interval = self.beamSizeFactor
        # halfLength = self.beamSizeFactor * 10
        halfLength = imageDensity/2 * interval
        # print halfLength, self.imageDensity/2, interval
        if halfLength > DFTSideLength/2:
            logger.warning('failed to update DFT grid')
            return False
        else:
            self.partialDFTGrid =self.createDFTGrid(
                    DFTSideLength, halfLength, interval)
            return True



    # def getAltAziFromRADEC(self, beamCoordinates, LSTDeg, arrayRefereceLatitude):
        # RA = np.deg2rad(beamCoordinates[:,0])
        # DEC = np.deg2rad(beamCoordinates[:,1])
        # LST = np.deg2rad(LSTDeg)


        # altitude, azimuth = coord.convertEquatorialToHorizontal(
                # RA, DEC, LST, arrayRefereceLatitude)

        # return altitude, azimuth


    def createDFTGrid(self, gridNum, halfLength, interval):
        ul = np.mgrid[0:halfLength:interval, 0:halfLength:interval]
        ur = np.mgrid[0:halfLength:interval, gridNum-halfLength:gridNum:interval]
        bl = np.mgrid[gridNum-halfLength:gridNum:interval, 0:halfLength:interval]
        br = np.mgrid[gridNum-halfLength:gridNum:interval, gridNum-halfLength:gridNum:interval]
        imagesCoord = np.array([
                        np.concatenate((
                        np.concatenate((ul[0].T, ur[0].T)).T,
                        np.concatenate((bl[0].T,br[0].T)).T)).flatten(),
                        np.concatenate((
                        np.concatenate((ul[1].T, ur[1].T)).T,
                        np.concatenate((bl[1].T, br[1].T)).T)).flatten()])

        # imagesCoord = np.array([
                        # np.vstack((
                        # np.hstack((ul[0], ur[0])),
                        # np.hstack((bl[0], br[0])))).flatten(),
                        # np.vstack((
                        # np.hstack((ul[1], ur[1])),
                        # np.hstack((bl[1], br[1])))).flatten()])

        return imagesCoord


    def calculateImageLength(self, rotatedProjectedBaselines, waveLength,
            zoomIn, density, gridNum, fixRange = None):

        if fixRange is None:
            uMax = np.amax(np.abs(rotatedProjectedBaselines[:,0]))/waveLength
            vMax = np.amax(np.abs(rotatedProjectedBaselines[:,1]))/waveLength
            uvMax = (uMax if uMax > vMax else vMax) * 2
        else:
            uvMax = fixRange

        # imageLength = 1/(1/(gridNum/uvMax))
        imageLength = gridNum/uvMax

        return imageLength

    def partialDFT(self, partialDFTGrid, rotatedProjectedBaselines, waveLength,
            imageLength, density, gridNum):
        "step: how many grids per uv unit"
        step = np.deg2rad(imageLength)
        halfGridNum = gridNum/2.
        uvSamples = []
        for baseline in rotatedProjectedBaselines:
            # print baseline
            uSlot = int(np.round(baseline[0]/waveLength*step + halfGridNum - 1))
            vSlot = int(np.round(halfGridNum - baseline[1]/waveLength*step - 1))

            uvSamples.append([uSlot, vSlot])
            uvSamples.append([gridNum - uSlot, gridNum - vSlot])
            # uvSamples.append([gridNum - uSlot - 1, gridNum - vSlot - 1])
            # uvGrids[vSlot][uSlot] = 1
            # uvGrids[-vSlot-1][-uSlot-1] = 1


        imagesCoord = partialDFTGrid
        # print imagesCoord
        fringeSum = np.zeros(density*density)
        for uv in uvSamples:
            fringeSum = fringeSum + np.exp(1j*np.pi*2*(imagesCoord[1]*uv[0] + imagesCoord[0]*uv[1])/gridNum)

        # fringeSum = np.sum(np.exp(1j*np.pi*2./gridNum*np.dot(np.fliplr(uvSamples), imagesCoord)), axis=0)

        fringeSum = 1.0 * fringeSum.reshape(density,density)/(len(uvSamples))

        # base = np.min(fringeSum.real)
        image = np.fft.fftshift(np.abs(fringeSum))
        # image = np.fft.fftshift(fringeSum.real - base)

        return image

    def calculateWeight(self, sidelength, antennaENU, resolution, boresight,
            waveLength, antennaWeights):

        speedOfLight = 299792458.

        # step = np.deg2rad(resolution)
        step = resolution
        margin = sidelength/2. * step
        # print(sidelength, margin*2, step)
        intervals = np.linspace(-margin-step , margin, sidelength)
        # lines = np.array(boresight).reshape(2,1) + intervals

        grids = np.meshgrid(intervals, intervals)
        coords = np.array(grids).reshape(2, -1).T

        equatorialCoordinates = coord.convert_pixel_coordinate_to_equatorial(
                coords, boresight.equatorial)
        altitude, azimuth = boresight.getAltAziFromRADEC(equatorialCoordinates,
                boresight.localSidereTime, np.deg2rad(boresight.arrayreference.geo[0]))

        theta = np.pi/2. - altitude
        phi = np.pi/2. - azimuth

        theta0 = np.pi/2. - boresight.horizontal[1]
        phi0 = np.pi/2.- boresight.horizontal[0]

        sourcePosition = np.array([np.sin(theta)*np.cos(phi),
            np.sin(theta)*np.sin(phi),
            np.cos(theta)])

        sourcePosition0 = np.array([np.sin(theta0)*np.cos(phi0),
        np.sin(theta0)*np.sin(phi0),
        np.cos(theta0)])

        receiverLocations = antennaENU
        kall = np.dot(receiverLocations, sourcePosition)/waveLength
        k0 = np.dot(receiverLocations, sourcePosition0)/waveLength
        k = (kall.T-k0).T


        weights = np.exp(1j * 2.0 * np.pi * k)
        # [antN, beams]
        if antennaWeights is not None and len(antennaWeights) == len(receiverLocations):
            weights = weights * np.array(antennaWeights).reshape(-1, 1)
        weights = weights.sum(axis=0)
        image = np.abs(weights).reshape(sidelength, sidelength)
        image = np.square(image/np.amax(image))
        # image = np.square(image/len(receiverLocations))
        # from matplotlib import pyplot as plt
        # plt.imshow(image, cmap="jet", origin="lower", vmin=0, vmax=1,
                # interpolation="nearest")
        # plt.tight_layout()
        # plt.savefig("weight.png")

        return image


    def performFFT(self, rotatedProjectedBaselines, waveLength, imageLength, gridNum):

        def Gaussian2DPDF(x, y, xMean, yMean, xSigma, ySigma):
            return np.exp(-((x-xMean)**2/(2*(xSigma**2)) + (y-yMean)**2/(2*(ySigma**2))))

        def normal(x, mu, sigma):
            return np.exp(-(x-mu)**2/(2*(sigma**2)))/(sigma*np.sqrt(2*np.pi))


        def FIRFilter(gridIndex, gridPosition, halfWidth):
            indexes = np.mgrid[gridIndex[1] - halfWidth:gridIndex[1] + halfWidth + 1:1,
                               gridIndex[0] - halfWidth:gridIndex[0] + halfWidth + 1:1]

            offsets = indexes - np.array([gridPosition[1], gridPosition[0]]).reshape(2,1,1)

            # dists = np.sqrt(np.sum(np.square(offsets.transpose(1,2,0)), axis = 2))

            values = np.exp(-1j*np.pi*offsets) * np.sinc(offsets)

            # gaussianValues = Gaussian2DPDF(offsets[1], offsets[0], 0, 0, 10.60, 10.60)
            # values = np.exp(-1j*np.pi*dists) * np.sinc(dists)

            # values = np.sum(values, axis=0) * gaussianValues

            # values = np.sum(values * normal(offsets, 0, 1), axis=0)
            values = np.sum(values, axis=0)

            return values

        step = imageLength
        # originalGridNum = gridNum
        gridNum = int(round(gridNum * 1.0))
        halfGridNum = gridNum/2.
        uvSamples = []
        uvGrids = np.zeros((gridNum, gridNum), dtype=np.complex)
        # uvGrids2 = np.zeros((gridNum, gridNum), dtype=np.complex)
        for baseline in rotatedProjectedBaselines:
            # uSlot = int((baseline[0]/waveLength*step + halfGridNum - 1))
            # vSlot = int((halfGridNum - baseline[1]/waveLength*step - 1))

            # uvGrids[vSlot][uSlot] += 1
            # uvGrids[-vSlot-1][-uSlot-1] += 1

            uSlot = baseline[0]/waveLength*step + halfGridNum - 1
            vSlot = halfGridNum - baseline[1]/waveLength*step - 1

            uGrid = int(round(uSlot))
            vGrid = int(round(vSlot))

            halfWidth = 5/2

            values = FIRFilter([uGrid, vGrid], [uSlot, vSlot], halfWidth)
            # valuesConj = FIRFilter([-uGrid, -vGrid], [-uSlot, -vSlot], halfWidth)
            valuesConj = np.rot90(values, 2)
            uvGrids[vGrid-halfWidth:vGrid+halfWidth+1, uGrid-halfWidth:uGrid+halfWidth+1] += values
            uvGrids[-(vGrid+halfWidth+1):-(vGrid-halfWidth), -(uGrid+halfWidth+1):-(uGrid-halfWidth)] += valuesConj



        np.savetxt('uvValues', np.abs(uvGrids))
        # psf = np.fft.ifft2(uvGrids, paddingToWidth)
        psf = np.fft.ifft2(uvGrids)
        # psf2 = np.fft.ifft2(uvGrids2, paddingToWidth)
        # image = np.fft.fftshift(np.abs(psf) + np.abs(psf2))
        image = np.fft.fftshift(np.abs(psf))
        image = image / image[int(halfGridNum)][int(halfGridNum)]

        # fullHalfWidth = paddingToWidth[0]/2
        # halfGridNum = int(originalGridNum/2)
        # trimmedImage = image[fullHalfWidth-1-halfGridNum:fullHalfWidth+halfGridNum-1, fullHalfWidth-1-halfGridNum:fullHalfWidth+halfGridNum-1]

        # return trimmedImage
        return image

    def getWCS(self):

        return self.WCS

    def calculateBeamScaleFromBaselines(self, rotatedProjectedBaselines, waveLength):
        baselineLengths = coord.distances(rotatedProjectedBaselines)
        baselineMax = np.amax(baselineLengths)
        baselineMin = np.amin(baselineLengths)
        # print baselineMax, baselineMin
        indexOfMaximum = np.argmax(baselineLengths)
        maxBaselineVector = rotatedProjectedBaselines[indexOfMaximum]
        # rotate vector on a surface https://math.stackexchange.com/questions/1830695/
        perpendicularOfMaxBaselineVector = coord.projectedRotate(
                np.pi/2., 0, maxBaselineVector, np.pi/2.)
        perpendicularUnitVector = perpendicularOfMaxBaselineVector/coord.distances(
                perpendicularOfMaxBaselineVector)
        perpendicularBaselines = np.dot(rotatedProjectedBaselines, perpendicularUnitVector)
        perpendicularBaselineMax = np.amax(np.abs(perpendicularBaselines))

        # print baselineMax
        minorAxis = np.rad2deg(1.22*waveLength/baselineMax/2.)
        majorAxis = np.rad2deg(1.22*waveLength/perpendicularBaselineMax/2.)

        return majorAxis, minorAxis

    def constructFitsHeader(self, density, step, boresight):
        """
        https://heasarc.gsfc.nasa.gov/docs/fcg/standard_dict.html
        CRVAL: coordinate system value at reference pixel
        CRPIX: coordinate system reference pixel
        CDELT: coordinate increment along axis
        CTYPE: name of the coordinate axis
        """
        self.WCS = {}
        # self.WCS['crpix'] = [density/2 -1, density/2 -1]
        # self.WCS['cdelt'] = [-step, step]
        # self.WCS['crval'] = [boresight[0] - abs(self.WCS['cdelt'][0]),
                             # boresight[1] - self.WCS['cdelt'][1]]
        # self.WCS['ctype'] = ["RA---TAN", "DEC--TAN"]

        self.WCS['crpix'] = [density/2. + 1, density/2. + 1]
        self.WCS['cdelt'] = [-step, step]
        self.WCS['crval'] = [boresight[0], boresight[1]]
        self.WCS['ctype'] = ["RA---TAN", "DEC--TAN"]

    def fitContour(self):
        # axis1, axis2, angle, overstep = calculateBeamSize(self.imageData,
                # self.imageDensity, self.imageLength, None, fit=True)
        # if overstep != 0:
            # elevation = np.rad2deg(self.boresight.horizontal[1])
            # if abs(elevation) < 20.:
                # logger.warning("Elevation is low %f" % elevation)
            # logger.warning("Beam shape probably is not correct. "
                           # "overstep at the power of {:.3}.".format(overstep))
        # angle = angle % 360. if abs(angle) > 360. else angle
        # self.beamAxis[0:3] = [axis1, axis2, angle]

        beamshapeModel = createBeamshapeModel(self.imageData, self.imageDensity,
                self.imageLength)

        return beamshapeModel

    def createContour(self, antennas, fileName=None, minAlt=0):


        self.array = coord.Array("main", antennas, self.arrayReferece)

        self.baselines = self.array.getBaselines()

        self.boresight = coord.Boresight("main", self.boresightCoord, self.observeTime,
                self.arrayReferece, self.boresightFrame)

        self.projectedBaselines = self.array.getRotatedProjectedBaselines(self.boresight)
        rotatedProjectedBaselines = self.projectedBaselines

        beamMajorAxisScale, beamMinorAxisScale = self.calculateBeamScaleFromBaselines(
                self.projectedBaselines, self.waveLength)
        # print self.waveLength, beamMinorAxisScale
        baselineMax = 1.22*self.waveLength/(beamMinorAxisScale*2)

        baselineNum = len(self.baselines)
        density = self.imageDensity
        gridNum = self.gridNumOfDFT
        imageLength = gridNum * self.resolution

        sidelength = density * self.beamSizeFactor
        windowLength = self.resolution * sidelength

        if baselineNum > 2 and self.autoZoom == True:
            axisRatio = beamMajorAxisScale/beamMinorAxisScale
            newBeamSizeFactor = axisRatio * beamMajorAxisScale*2*1.3 / (self.resolution * density)
            # print newBeamSizeFactor,
            if baselineMax > 2e3:
                logger.debug('larger')
                # print 'larger'
                newBeamSizeFactor = 6 if newBeamSizeFactor < 6 else int(round(newBeamSizeFactor))
            else:
                # print 'smaller', newBeamSizeFactor
                logger.debug('smaller')
                if newBeamSizeFactor > 6:
                    newBeamSizeFactor = 6
                elif newBeamSizeFactor > 1.:
                    newBeamSizeFactor = int(round(newBeamSizeFactor))
                else:
                    newBeamSizeFactor = 1

            # if newBeamSizeFactor < 3:
                # newBeamSizeFactor = 3 if baselineMax < 1e3 else 6
            # elif newBeamSizeFactor > 12:
                # newBeamSizeFactor = 10 if baselineMax < 1e3 else 20
            # else:
                # newBeamSizeFactor = int(round(newBeamSizeFactor))

            # print newBeamSizeFactor,
            # print newBeamSizeFactor
            self.setBeamSizeFactor(newBeamSizeFactor)

            sidelength = density * self.beamSizeFactor
            windowLength = self.resolution * sidelength
            imageLength = gridNum * self.resolution

            # image = self.partialDFT(self.partialDFTGrid, rotatedProjectedBaselines,
                    # self.waveLength, imageLength, density, gridNum)

            image  = self.calculateWeight(density, self.array.getENU(),
                self.resolution*self.beamSizeFactor,
                self.boresight, self.waveLength, self.antennaWeights)


            #if fileName != None:
            #    plotBeamContour(image, (0,0), windowLength,
            #        interpolation = self.interpolating, fileName='contourTest.png')

            sizeInfo = calculateBeamSize(image, density, windowLength, np.rad2deg(beamMajorAxisScale))
            # self.beamAxis = [sizeInfo[0], sizeInfo[1], sizeInfo[2]]
            # print sizeInfo[3]
            if sizeInfo[3] != 0:
                sigmaTest = normSigma(windowLength/2., 0, sizeInfo[3])
                majorAxis = normInverse(0.4, 0, sigmaTest)
                # majorAxis =  beamMajorAxisScale / (sizeInfo[3]/0.5)
            else:
                majorAxis, minorAxis, angle = sizeInfo[0], sizeInfo[1], sizeInfo[2]
            # print np.deg2rad(majorAxis), beamMajorAxisScale
            newBeamSizeFactor = 2*majorAxis*1.7 / (self.resolution *  density)
            # print newBeamSizeFactor
            overstep = sizeInfo[3]
            if overstep != 0:
                newBeamSizeFactor += 3.5
            # print newBeamSizeFactor
            # print majorAxis, minorAxis, angle
            if newBeamSizeFactor < 1.:
                sidelength = density * newBeamSizeFactor
                windowLength = self.resolution * sidelength
                if windowLength / (2.*majorAxis*1.4) < 1.1:
                    newBeamSizeFactor = 2
                else:
                    newBeamSizeFactor = 1
            else:
                newBeamSizeFactor = int(round(newBeamSizeFactor))
            self.setBeamSizeFactor(newBeamSizeFactor)
            # imageLength = self.calculateImageLength(rotatedProjectedBaselines,
                # self.waveLength, self.beamSizeFactor, density, gridNum, fixRange=width)


            # image = self.partialDFT(self.partialDFTGrid, rotatedProjectedBaselines,
                    # self.waveLength, imageLength, density, gridNum)

            image  = self.calculateWeight(density, self.array.getENU(),
                self.resolution*self.beamSizeFactor,
                self.boresight, self.waveLength, self.antennaWeights)


        else:
            # image = self.partialDFT(self.partialDFTGrid, rotatedProjectedBaselines,
                    # self.waveLength, imageLength, density, gridNum)
            image  = self.calculateWeight(density, self.array.getENU(),
                self.resolution*self.beamSizeFactor,
                self.boresight, self.waveLength, self.antennaWeights)


        sidelength = density * self.beamSizeFactor
        windowLength = self.resolution * sidelength

        self.imageLength = windowLength

        upper_left_pixel = [-windowLength/2., windowLength/2.] # x,y
        bottom_right_pixel = [windowLength/2., -windowLength/2.] # x,y

        coordinates_equatorial = coord.convert_pixel_coordinate_to_equatorial(
            [upper_left_pixel, bottom_right_pixel], self.boresight.equatorial)
        equatorial_range = [coordinates_equatorial[0][0], coordinates_equatorial[1][0], # left, right
                            coordinates_equatorial[0][1], coordinates_equatorial[1][1]] # up, bottom

        self.beamAxis = [None, None, None, equatorial_range]
        self.constructFitsHeader(density, self.resolution*self.beamSizeFactor, self.boresight.equatorial)

        self.psf = PointSpreadFunction(image, self.boresight,
                windowLength, self.WCS, equatorial_range)

        if fileName is not None:
            plotBeamContour(image, self.boresight.equatorial, equatorial_range,
                    interpolation = self.interpolating, fileName = fileName, colormap = True)

        # if baselineNum > 2:
            # resolution = windowLength/density
            # sizeInfo = calculateBeamSize(image, density, windowLength, beamMajorAxisScale, fit=True)
            # if sizeInfo[3] != 0:
                # elevation = np.rad2deg(self.boresight.horizontal[1])
                # if abs(elevation) < 20.:
                    # logger.warning("Elevation is low %f" % elevation)
                # logger.warning("Beam shape probably is not correct.")
            # self.beamAxis[0:3] = [sizeInfo[0], sizeInfo[1], sizeInfo[2]]

        self.imageData = image
        logger.info("beamshape simulation inputs, freq: {:.5g}, source: {}, "
                    "time: {}, subarray: {}".format(299792458./self.waveLength,
                        self.boresightInput, self.observeTime,
                        [ant.name for ant in antennas] if self.antennaWeights is None
                        else ["{}:{}".format(ant.name, antWeight) for ant, antWeight in zip(antennas, self.antennaWeights)]))
        return image
