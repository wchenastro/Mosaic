#! /usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

import coordinate as coord
from plot import plotBeamContour
from utilities import normSigma, normInverse
from beamshape import calculateBeamSize, trackBorder

import inspect, pickle, datetime, logging

loggerFormat = '%(asctime)-15s  %(filename)s  %(message)s'
logging.basicConfig(format = loggerFormat, level=logging.INFO)
logger = logging.getLogger(__name__)


class PointSpreadFunction(object):
    """
    class for point spread function
    """

    def __init__(self, image, bore_sight, width):
        self.image = image
        self.bore_sight = bore_sight
        self.width = width

class InterferometryObservation:

    equatorialInput = 0
    horizontalInput = 1

    def __init__(self, arrayRefereceGEODET, observeTime, waveLength):
        self.arrayRefereceGEODET = arrayRefereceGEODET
        self.observeTime = observeTime
        self.waveLength = waveLength
        self.baselines = []
        self.boreSightHorizontal = []
        self.projectedBaselines = []
        self.baselineIndex = []
        self.beamAxis = []
        self.beamSizeFactor = 1
        self.beamCoordinates = []
        self.beamNumber = 400
        self.localSiderealTime = 0
        self.beamSize = 1.22*self.waveLength/13.5
        self.interpolating = True
        self.amplitude = []
        self.autoZoom = True
        self.gridNumOfDFT = 100000.0
        self.imageDensity = 20
        self.projectedBaselines = []
        self.resolution = 1/3600.0
        self.inputType = self.equatorialInput
        self.WCS = {}



    def setInterpolating(self, state):
        if state == True:
            self.interpolating = True
        else:
            self.interpolating = False

    def setInputType(self, inputType):
        self.inputType = inputType

    def setObserveTime(self, dateTime):
        if type(dateTime) != datetime.datetime:
            dateTime = coord.epochToDatetime(dateTime)
        self.observeTime = dateTime

    def getObserveTime(self):
        return self.observeTime

    def getPointSpreadFunction(self):
        return self.psf

    def setBeamSizeFactor(self, size, autoZoom = True):
        if size != self.beamSizeFactor:
            oldSize = self.beamSizeFactor
            self.beamSizeFactor = size
            isUpdated = self.updateBeamCoordinates()
            if isUpdated:
                return True
            else:
                self.beamSizeFactor = oldSize
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
            isUpdated = self.updateBeamCoordinates()
            if isUpdated:
                return True
            else:
                self.beamNumber = oldBeamNumber
                self.setImageDensity(oldDensity)

    def setResolution(self, resolution):
        '''resolution default in arc second deg for now'''
        self.resolution = resolution/3600.0

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

    def setBoreSight(self, beamBoreSight):
        self.boreSight = beamBoreSight
        self.updateBeamCoordinates()

    def getBoreSight(self):
        return self.boreSight

    def getHorizontal(self):
        return self.boreSightHorizontal

    def setHorizontal(self, horizontal):
        self.boreSightHorizontal = np.deg2rad(horizontal)

    def getProjectedBaselines(self):
        uvCoord = self.projectedBaselines/self.waveLength
        return np.concatenate((uvCoord, -uvCoord))

    def getBeamAxis(self):
        return self.beamAxis

    def getImageLength(self):
        return self.imageLength

    def getImageDensity(self):
        return self.imageDensity

    def setImageDensity(self, density):
        self.imageDensity = density

    def getAntCoordinates(self):
        return self.antCoordinatesGEODET.tolist()

    def createBaselines(self, antCoordinatesENU):
        baselines = []
        index = 1
        for antenna1 in antCoordinatesENU:
            for antenna2 in antCoordinatesENU[index:]:
                baselines.append([antenna1[0]-antenna2[0],
                        antenna1[1]-antenna2[1],antenna1[2]-antenna2[2]])
            index += 1

        return baselines

    def updateBeamCoordinates(self):
        interval = self.beamSizeFactor
        # halfLength = self.beamSizeFactor * 10
        halfLength = self.imageDensity/2 * interval
        # print halfLength, self.imageDensity/2, interval
        if halfLength > self.gridNumOfDFT/2:
            return False
        else:
            self.partialDFTGrid = self.createDFTGrid(self.gridNumOfDFT, halfLength, interval)
            return True



    def getAltAziFromRADEC(self, beamCoordinates, LSTDeg, arrayRefereceLatitude):
        RA = np.deg2rad(beamCoordinates[:,0])
        DEC = np.deg2rad(beamCoordinates[:,1])
        LST = np.deg2rad(LSTDeg)


        altitude, azimuth = coord.convertEquatorialToHorizontal(
                RA, DEC, LST, arrayRefereceLatitude)

        return altitude, azimuth


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
            uSlot = int((baseline[0]/waveLength*step + halfGridNum - 1))
            vSlot = int((halfGridNum - baseline[1]/waveLength*step - 1))

            uvSamples.append([uSlot, vSlot])
            uvSamples.append([gridNum - uSlot - 1, gridNum - vSlot - 1])
            # uvGrids[vSlot][uSlot] = 1
            # uvGrids[-vSlot-1][-uSlot-1] = 1


        imagesCoord = partialDFTGrid
        fringeSum = np.zeros(density*density)
        for uv in uvSamples:
            fringeSum = fringeSum + np.exp(1j*np.pi*2*(imagesCoord[1]*uv[0] + imagesCoord[0]*uv[1])/gridNum)

        # fringeSum = np.sum(np.exp(1j*np.pi*2./gridNum*np.dot(np.fliplr(uvSamples), imagesCoord)), axis=0)

        fringeSum = fringeSum.reshape(density,density)/(len(uvSamples))

        image = np.fft.fftshift(np.abs(fringeSum))

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


    def createContour(self, antennacoor, fileName=None, minAlt=0):


        self.antCoordinatesGEODET = np.array(antennacoor)

        antCoordinatesECEF = coord.convertGodeticToECEF(self.antCoordinatesGEODET)

        arrayRefereceECEF = coord.convertGodeticToECEF([self.arrayRefereceGEODET])[0]
        antCoordinatesENU = coord.convertECEFToENU(antCoordinatesECEF, arrayRefereceECEF, self.arrayRefereceGEODET)


        self.baselines = self.createBaselines(antCoordinatesENU)

        arrayRefereceLatitude = np.deg2rad(self.arrayRefereceGEODET[0])
        LSTDeg =  coord.calculateLocalSiderealTime(self.observeTime, self.arrayRefereceGEODET[1])

        if self.inputType == self.equatorialInput:

            altitude, azimuth = self.getAltAziFromRADEC(np.array([self.boreSight]),
                        LSTDeg, arrayRefereceLatitude)
            self.boreSightHorizontal = (azimuth[0], altitude[0])

        elif self.inputType == self.horizontalInput:

            azimuth, altitude  = self.boreSightHorizontal
            RA, DEC = coord.convertHorizontalToEquatorial(azimuth, altitude, np.deg2rad(LSTDeg), arrayRefereceLatitude)
            self.boreSight = np.rad2deg([RA, DEC])
            print np.rad2deg([RA, DEC])


        self.saveParas('paras')

        # if abs(altitude[0]) < minAlt:
            # return
        # LHA = LST - RA
        LHA = np.deg2rad(LSTDeg) - np.deg2rad(self.boreSight[0])

        rotatedENU = coord.rotateENUToEquatorialPlane(self.baselines,
                arrayRefereceLatitude, self.boreSightHorizontal[0], self.boreSightHorizontal[1])



        # self.projectedBaselines = sv.projectedBaselinesDROP(
                # np.array([altitude[0]]), np.array([azimuth[0]]), self.baselines)

        # rotatedProjectedBaselines = sv.rotateCoordinateDROP(self.projectedBaselines,
                # np.pi/2.0 - altitude[0], azimuth[0])

        rotatedProjectedBaselines = coord.projectBaselines(rotatedENU, LHA, np.deg2rad(self.boreSight[1]))

        self.projectedBaselines = rotatedProjectedBaselines

        beamMajorAxisScale, beamMinorAxisScale = self.calculateBeamScaleFromBaselines(
                rotatedProjectedBaselines, self.waveLength)
        # print self.waveLength, beamMinorAxisScale
        baselineMax = 1.22*self.waveLength/(beamMinorAxisScale*2)

        baselineNum = len(self.baselines)
        density = self.imageDensity
        gridNum = self.gridNumOfDFT

        # width = 1/self.resolution
        # imageLength = self.calculateImageLength(rotatedProjectedBaselines,
                # self.waveLength, self.beamSizeFactor, density, gridNum, fixRange=width)
        imageLength = gridNum * self.resolution

        # newBeamSizeFactor = beamMajorAxisScale*1.3 / (imageLength/gridNum) / density
        # newBeamSizeFactor = beamMinorAxisScale*4*1.3 / (imageLength/gridNum) / density


        # if newBeamSizeFactor < 1:
            # newBeamSizeFactor = 1
        # else:
            # newBeamSizeFactor = int(round(newBeamSizeFactor))

        # if baselineNum > 2 and self.autoZoom == True and abs(newBeamSizeFactor - self.beamSizeFactor) > 0:

            # self.beamSizeFactor = newBeamSizeFactor
            # print "new beam factor:", newBeamSizeFactor
            # self.setBeamSizeFactor(newBeamSizeFactor)
        sidelength = density * self.beamSizeFactor
        windowLength = self.resolution * sidelength

        # image = self.partialDFT(self.partialDFTGrid, rotatedProjectedBaselines,
                # self.waveLength, imageLength, self.beamSizeFactor, density, gridNum)
        if baselineNum > 2 and self.autoZoom == True:
            # newBeamSizeFactor = beamMajorAxisScale*10*1.3 / (imageLength/gridNum) / density
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

            image = self.partialDFT(self.partialDFTGrid, rotatedProjectedBaselines,
                    self.waveLength, imageLength, density, gridNum)

            if fileName != None:
                plotBeamContour(image, (self.boreSight), windowLength,
                    interpolation = self.interpolating, fileName='contourTest.png')

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


            imageLength = gridNum * self.resolution
            image = self.partialDFT(self.partialDFTGrid, rotatedProjectedBaselines,
                    self.waveLength, imageLength, density, gridNum)

        else:
            image = self.partialDFT(self.partialDFTGrid, rotatedProjectedBaselines,
                    self.waveLength, imageLength, density, gridNum)


        sidelength = density * self.beamSizeFactor
        windowLength = self.resolution * sidelength

        self.imageLength = windowLength
        self.psf = PointSpreadFunction(image, self.boreSight, windowLength)
        if fileName != None:
            plotBeamContour(image, self.boreSight, windowLength,
                    interpolation = self.interpolating)

        if baselineNum > 2:
            sizeInfo = calculateBeamSize(image, density, windowLength, beamMajorAxisScale)
            if sizeInfo[3] != 0:
                elevation = np.rad2deg(self.boreSightHorizontal[1])
                if elevation < 20.:
                    logger.warning("Elevation is low %f" % elevation)
                logger.warning("Beam shape probably is not correct.")
            self.beamAxis = [sizeInfo[0], sizeInfo[1], sizeInfo[2]]

    def createPSF(self, antennacoor, waveLengths, writer, plotting):

        self.antCoordinatesGEODET = np.array(antennacoor)

        antCoordinatesECEF = coord.convertGodeticToECEF(self.antCoordinatesGEODET)

        arrayRefereceECEF = coord.convertGodeticToECEF([self.arrayRefereceGEODET])[0]
        antCoordinatesENU = coord.convertECEFToENU(antCoordinatesECEF, arrayRefereceECEF, self.arrayRefereceGEODET)

        self.baselines = self.createBaselines(antCoordinatesENU)

        arrayRefereceLatitude = np.deg2rad(self.arrayRefereceGEODET[0])
        LSTDeg =  coord.calculateLocalSiderealTime(self.observeTime, self.arrayRefereceGEODET[1])

        if self.inputType == self.equatorialInput:

            altitude, azimuth = self.getAltAziFromRADEC(np.array([self.boreSight]),
                        LSTDeg, arrayRefereceLatitude)
            self.boreSightHorizontal = (azimuth[0], altitude[0])
            # print("source in Azi, Alt is %f, %f" % (np.rad2deg(self.boreSightHorizontal[0]), np.rad2deg(self.boreSightHorizontal[1])))

        elif self.inputType == self.horizontalInput:

            azimuth, altitude  = self.boreSightHorizontal
            RA, DEC = coord.convertHorizontalToEquatorial(azimuth, altitude, np.deg2rad(LSTDeg), arrayRefereceLatitude)
            self.boreSight = np.rad2deg([RA, DEC])
            print("source in RA, DEC is %f, %f" % (self.boreSight[0], self.boreSight[1]))


        rotatedENU = coord.rotateENUToEquatorialPlane(self.baselines,
                arrayRefereceLatitude, self.boreSightHorizontal[0], self.boreSightHorizontal[1])




        # LHA = LST - RA
        LHA = np.deg2rad(LSTDeg) - np.deg2rad(self.boreSight[0])


        rotatedProjectedBaselines = coord.projectBaselines(rotatedENU, LHA, np.deg2rad(self.boreSight[1]))

        self.projectedBaselines = rotatedProjectedBaselines


        density = self.imageDensity
        gridNum = self.gridNumOfDFT

        width = np.rad2deg(1/self.resolution)
        imageLength = self.calculateImageLength(rotatedProjectedBaselines,
                None, self.beamSizeFactor, density, gridNum, fixRange=width)

        sidelength = density * self.beamSizeFactor
        windowLength = imageLength/gridNum*sidelength


        """
        https://heasarc.gsfc.nasa.gov/docs/fcg/standard_dict.html
        CRVAL: coordinate system value at reference pixel
        CRPIX: coordinate system reference pixel
        CDELT: coordinate increment along axis
        CTYPE: name of the coordinate axis
        """
        self.WCS['crpix'] = [density/2 -1, density/2 -1]
        self.WCS['cdelt'] = np.rad2deg([imageLength/gridNum, imageLength/gridNum])
        self.WCS['crval'] = [self.boreSight[0] - self.WCS['cdelt'][0],
                             self.boreSight[1] - self.WCS['cdelt'][1]]
        self.WCS['ctype'] = ["RA---TAN", "DEC--TAN"]

        chanIdx = 0
        images = []
        for waveLength in waveLengths:

            # image = self.partialDFT(self.partialDFTGrid, rotatedProjectedBaselines,
                    # waveLength, imageLength, self.beamSizeFactor, density, gridNum)

            # paddingToWidth = gridNum * width/(gridnum/imageLength)
            # fullWidth = int(round(width * imageLength))
            # paddingToWidth = [fullWidth, fullWidth]

            image = self.performFFT(rotatedProjectedBaselines, waveLength, imageLength, gridNum)

            # flippedImage = np.fliplr(image)
            if plotting == True:
                self.imageLength = windowLength
                psfName = 'psf' + str(chanIdx) + '.png'
                chanIdx += 1
                plotBeamContour(image, np.deg2rad(self.boreSight), windowLength,
                        fileName = psfName, interpolation = self.interpolating)
            # writer(waveLength, image)
            images.append(image)

        return image

