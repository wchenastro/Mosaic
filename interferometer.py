#! /usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

# import steervec as sv
import coordinate as coord
import beamShape as bs
# import createBeam as cb

import inspect
from scipy import interpolate


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
        self.resolution = np.deg2rad(1/3600.0)
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
        self.observeTime = dateTime

    def getObserveTime(self):
        return self.observeTime

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
        self.resolution = np.deg2rad(resolution/3600.0)

    def getBaselines(self):
        return self.baselines

    def getAmplitude(self):
        return self.amplitude

    def getBeamCoordinates(self):
        return self.beamCoordinates

    def getBaselinesNumber(self):
        return len(self.baselines)

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

    def calculateBeamOverlaps(self, ellipseCenters, radius, majorAxis, minorAxis, rotation):

        def isInsideEllips(center, majorAxis, minorAxis, rotation, testPointX, testPointY):

            xOffset = testPointX - center[0]
            yOffset = testPointY - center[1]
            cosa = np.cos(rotation)
            sina = np.sin(rotation)
            # majorSquare = np.power(majorAxis,2)
            # minorSquare = np.power(minorAxis,2)

            result = np.power((cosa*xOffset + sina*yOffset)/majorAxis,2) +\
                     np.power((sina*xOffset - cosa*yOffset)/minorAxis,2)


            # result = 1/majorAxis**2 * ((testPointX - center[0])*np.cos(rotation) +\
                    # (testPointY - center[1])*np.sin(rotation))**2 +\
                    # 1/minorAxis**2 * ((testPointX - center[0])*np.sin(rotation) -\
                    # (testPointY - center[1])*np.cos(rotation))**2
            return result

        longAxis = majorAxis if majorAxis > minorAxis else minorAxis
        gridNum = 1000
        # print 0.3*radius, longAxis
        halfSidelength = 0.3 * radius
        offsetCenter = [radius, radius]
        # np.savetxt('tmp/oldcenter', ellipseCenters)
        ellipseCenters = ellipseCenters + offsetCenter
        # np.savetxt('tmp/newcenter', ellipseCenters)
        innerEllipses = []
        for ellipseCenter in ellipseCenters:
            if (ellipseCenter[0] > (offsetCenter[0]-halfSidelength) and\
                ellipseCenter[0] < (offsetCenter[0]+halfSidelength)) and\
               (ellipseCenter[1] > (offsetCenter[1]-halfSidelength) and\
                ellipseCenter[1] < (offsetCenter[1]+halfSidelength)):
                innerEllipses.append(ellipseCenter)

        halfSidelength += longAxis
        # np.savetxt('tmp/innercenter', innerEllipses)
        step = 2*halfSidelength/gridNum
        width = longAxis
        squareEdgeX = [offsetCenter[0] - halfSidelength, offsetCenter[0] + halfSidelength]
        squareEdgeY = [offsetCenter[1] - halfSidelength, offsetCenter[1] + halfSidelength]
        # print squareEdgeX, squareEdgeY

        grids = np.mgrid[squareEdgeY[0]:squareEdgeY[1]:step, squareEdgeX[0]:squareEdgeX[1]:step]

        # nopoint = []
        gridLength = grids.shape[1]
        overlapCounter = np.zeros((gridLength, gridLength))
        for ellipseCenter in innerEllipses:
            horizontalBoarder = [ellipseCenter[0]-width, ellipseCenter[0]+width]
            verticalBoarder = [ellipseCenter[1]-width, ellipseCenter[1]+width]

            horizontalIndex = np.round([(horizontalBoarder[0]-squareEdgeX[0])/(2.0*halfSidelength)*gridNum,
                    (horizontalBoarder[1]-squareEdgeX[0])/(2.0*halfSidelength)*gridNum]).astype(int)
            verticalIndex = np.round([(verticalBoarder[0]-squareEdgeY[0])/(2.0*halfSidelength)*gridNum,
                    (verticalBoarder[1]-squareEdgeY[0])/(2.0*halfSidelength)*gridNum]).astype(int)
            # print verticalIndex, horizontalIndex
            gridX = grids[1][verticalIndex[0]: verticalIndex[1], horizontalIndex[0]: horizontalIndex[1]]
            gridY = grids[0][verticalIndex[0]: verticalIndex[1], horizontalIndex[0]: horizontalIndex[1]]
            # print gridX
            # print gridY
            result = isInsideEllips(ellipseCenter, majorAxis, minorAxis, rotation, gridX, gridY)
            mask = result<1
            result[mask] = 1
            result[~mask] = 0
            # print result.shape
            overlapCounter[verticalIndex[0]: verticalIndex[1], horizontalIndex[0]: horizontalIndex[1]] += result
            # print ellipseCenter, majorAxis, minorAxis, rotation
            # np.save('tmp/grid', [gridX, gridY])
            # np.savetxt('tmp/pointm', result)
            # exit(0)
            # if np.amin(result) > 1.:
                # np.savetxt('tmp/grid', [gridX,gridY])
                # print ellipseCenter
                # exit()
            # print len(gridX)
            # print result[result<1]
            # print len(gridY), np.amin(result), result[result<1]
        # np.savetxt('tmp/nopoint', nopoint)
        print np.count_nonzero(overlapCounter > 1), np.count_nonzero(overlapCounter == 1), np.count_nonzero(overlapCounter == 0)
        np.save('overlapCounter', overlapCounter)



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
            imageLength, zoomIn, density, gridNum):

        step = imageLength
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
        minorAxis = 1.22*waveLength/baselineMax/2.
        majorAxis = 1.22*waveLength/perpendicularBaselineMax/2.

        return majorAxis, minorAxis


    def trackBorder(self, image, threshold = 0.3, density = 20, interpolatedLength = 800):

        interpolater = interpolate.interp2d(range(density), range(density), image ,kind='cubic')
        image = interpolater(np.linspace(0, density - 1, interpolatedLength),
                np.linspace(0, density - 1, interpolatedLength))
        # np.savetxt('interpolate', image)
        interpolatedGrid = interpolatedLength*1.0/density


        imageCenter = [interpolatedLength/2, interpolatedLength/2]
        rowIdx = int(imageCenter[0])
        colIdx = int(imageCenter[1])

        trueCenterIndex = [int(rowIdx + interpolatedGrid/2 + 1), int(colIdx + interpolatedGrid/2 + 1)]

        class State:
            move, addRow, findBorder, findBorderReverse, upsideDown, end = range(6)

        border = []
        state = State.move
        rowStep = -1
        colStep = 1
        colBorder = interpolatedLength - 1
        rowBorder = 0
        bottom = False
        overstep = False
        offset = 0.1
        # closestToCenter = 1
        # closestToCenterIndex = []

        while state != State.end:
            if state == State.move:
                if colIdx != colBorder:
                    colIdx += colStep
                    if image[rowIdx, colIdx] < threshold:
                        border.append([rowIdx, colIdx])
                        state = State.addRow
                    # else:
                        # distToCenter = 1 - image[rowIdx, colIdx]
                        # if distToCenter  < closestToCenter:
                            # closestToCenter = distToCenter
                            # closestToCenterIndex = [rowIdx, colIdx]
                else:
                    overstep = True
                    state = State.addRow

            if state == State.addRow:
                if rowIdx != rowBorder:
                    rowIdx += rowStep
                    state = State.findBorder
                else:
                    overstep = True
                    state = State.upsideDown


            if state == State.findBorder:
                if image[rowIdx, colIdx] < threshold:
                    colStep *= -1
                    colBorder = interpolatedLength - 1 - colBorder
                    state = State.findBorderReverse
                else:
                    if colIdx != colBorder:
                        colIdx += colStep
                    else:
                        overstep = True
                        colStep *= -1
                        colBorder = interpolatedLength - 1 - colBorder
                        state = State.move

            if state == State.findBorderReverse:
                try:
                    if image[rowIdx, colIdx] < threshold:
                        if colIdx != colBorder:
                            colIdx += colStep
                        else:
                            state = State.upsideDown

                    else:
                        if len(border) > 2:
                            distLastSQ = (border[-1][0] - border[-2][0])**2 + (border[-1][1] - border[-2][1])**2
                            distNowSQ = (border[-1][0] - rowIdx)**2 + (border[-1][1] - colIdx)**2
                            if (distNowSQ - distLastSQ*1.5) > 0:
                                state = State.upsideDown
                                continue
                        border.append([rowIdx, colIdx])
                        state = State.move
                except IndexError:
                    print rowIdx, colIdx
                    raise

            if state == State.upsideDown:
                if bottom == True:
                    state = State.end
                else:
                    bottom = True
                    rowStep = 1
                    colStep = -1
                    colBorder = 0
                    rowBorder = interpolatedLength - 1
                    rowIdx = imageCenter[0]
                    colIdx = imageCenter[1]
                    state = State.move


        return border, trueCenterIndex, overstep

    def calculateBeamSize(self, image, density, windowLength, beamMajorAxisScale, interpolatedLength = 800, threshold = 0.4):
        interpolatedLength = interpolatedLength
        threshold = threshold
        border, closestToCenterIndex, overstep = self.trackBorder(
                image, threshold, density, interpolatedLength)

        if overstep == True: print 'overstep'
        np.savetxt('border', border)
        if len(border) < 10:
            print 'too little points'
            return None

        imageArray = np.array(border) - [0, closestToCenterIndex[1]]
        imageArray[:, 0] = closestToCenterIndex[0] - imageArray[:, 0]
        # np.savetxt('border', imageArray)
        distancesSQ = np.sum(np.square(imageArray), axis=1)
        minDistIndex = np.argmin(distancesSQ)
        minDist = np.sqrt(distancesSQ[minDistIndex])
        minDistVector = imageArray[minDistIndex]
        if overstep == False:
            maxDistIndex = np.argmax(distancesSQ)
            maxDist = np.sqrt(distancesSQ[maxDistIndex])
            maxDistVector = imageArray[maxDistIndex]
            angle = np.arctan2(maxDistVector[0], maxDistVector[1])

            minorAxis  = np.rad2deg(windowLength/interpolatedLength*minDist)
            majorAxis  = np.rad2deg(windowLength/interpolatedLength*maxDist)
        else:
            angle = np.pi/2. + np.arctan2(minDistVector[0], minDistVector[1])


            majorAxis = np.rad2deg(beamMajorAxisScale)
            minorAxis  = np.rad2deg(windowLength/interpolatedLength*minDist)

        return majorAxis, minorAxis, angle

    def createContour(self, antennacoor, fileName='contour.png', minAlt=0):


        self.antCoordinatesGEODET = np.array(antennacoor)

        antCoordinatesECEF = coord.convertGodeticToECEF(self.antCoordinatesGEODET)

        arrayRefereceECEF = coord.convertGodeticToECEF([self.arrayRefereceGEODET])[0]
        antCoordinatesENU = coord.convertECEFToENU(antCoordinatesECEF, arrayRefereceECEF, self.arrayRefereceGEODET)


        self.baselines = bs.createBaselines(antCoordinatesENU)

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
            newBeamSizeFactor = axisRatio * 2 * beamMajorAxisScale*2*1.3 / (self.resolution * density)
            if newBeamSizeFactor < 3:
                newBeamSizeFactor = 3
            elif newBeamSizeFactor > 12:
                newBeamSizeFactor = 12
            else:
                newBeamSizeFactor = int(round(newBeamSizeFactor))

            print newBeamSizeFactor
            self.setBeamSizeFactor(newBeamSizeFactor)

            sidelength = density * self.beamSizeFactor
            windowLength = self.resolution * sidelength

            image = self.partialDFT(self.partialDFTGrid, rotatedProjectedBaselines,
                    self.waveLength, imageLength, newBeamSizeFactor, density, gridNum)

            bs.plotBeamContour3(image, np.deg2rad(self.boreSight), windowLength,
                interpolation = self.interpolating, fileName='contourTest.png')

            sizeInfo = self.calculateBeamSize(image, density, windowLength, beamMajorAxisScale)
            # self.beamAxis = [sizeInfo[0], sizeInfo[1], sizeInfo[2]]
            majorAxis, minorAxis, angle = sizeInfo[0], sizeInfo[1], sizeInfo[2]
            # print np.deg2rad(majorAxis), beamMajorAxisScale
            newBeamSizeFactor = 2*np.deg2rad(majorAxis)*1.4 / (self.resolution *  density)
            if newBeamSizeFactor < 1:
                    newBeamSizeFactor = 1
            else:
                newBeamSizeFactor = int(round(newBeamSizeFactor))
            self.setBeamSizeFactor(newBeamSizeFactor)
            # imageLength = self.calculateImageLength(rotatedProjectedBaselines,
                # self.waveLength, self.beamSizeFactor, density, gridNum, fixRange=width)


            image = self.partialDFT(self.partialDFTGrid, rotatedProjectedBaselines,
                    self.waveLength, imageLength, self.beamSizeFactor, density, gridNum)

        else:
            image = self.partialDFT(self.partialDFTGrid, rotatedProjectedBaselines,
                    self.waveLength, imageLength, self.beamSizeFactor, density, gridNum)


        sidelength = density * self.beamSizeFactor
        windowLength = self.resolution * sidelength

        self.imageLength = windowLength
        bs.plotBeamContour3(image, np.deg2rad(self.boreSight), windowLength,
                interpolation = self.interpolating)

        if baselineNum > 2:
            sizeInfo = self.calculateBeamSize(image, density, windowLength, beamMajorAxisScale)
            self.beamAxis = [sizeInfo[0], sizeInfo[1], sizeInfo[2]]

    def createPSF(self, antennacoor, waveLengths, writer, plotting):

        self.antCoordinatesGEODET = np.array(antennacoor)

        antCoordinatesECEF = coord.convertGodeticToECEF(self.antCoordinatesGEODET)

        arrayRefereceECEF = coord.convertGodeticToECEF([self.arrayRefereceGEODET])[0]
        antCoordinatesENU = coord.convertECEFToENU(antCoordinatesECEF, arrayRefereceECEF, self.arrayRefereceGEODET)

        self.baselines = bs.createBaselines(antCoordinatesENU)

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

        width = 1/self.resolution
        imageLength = self.calculateImageLength(rotatedProjectedBaselines,
                None, self.beamSizeFactor, density, gridNum, fixRange=width)

        sidelength = density * self.beamSizeFactor
        windowLength = imageLength/gridNum*sidelength


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
                bs.plotBeamContour3(image, np.deg2rad(self.boreSight), windowLength,
                        fileName = psfName, interpolation = self.interpolating)
            # writer(waveLength, image)
            images.append(image)

        return image

