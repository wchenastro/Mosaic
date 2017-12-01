#! /usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

import steervec as sv
import coordinate as coord
import beamShape as bs
import createBeam as cb

import inspect
from scipy import interpolate


class InterferometryObservation:
    class Baseline:
        def __init__(longitude1, altitude1, longitude2, altitude2):
            self.longitude1 = longitude1
            self.altitude1 = altitude1
            self.longitude2 = longitude2
            self.altitude2 = altitude2


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
        self.beamSize = 1.22*self.waveLength/13.5
        self.interpolating = True
        self.amplitude = []
        self.autoZoom = True
        self.gridNumOfDFT = 1000.0

    def setInterpolating(self, state):
        if state == True:
            self.interpolating = True
        else:
            self.interpolating = False

    def setObserveTime(self, dateTime):
        self.observeTime = dateTime

    def getObserveTime(self):
        return self.observeTime

    def setBeamSizeFactor(self, size, autoZoom = True):
        if size != self.beamSizeFactor:
            self.beamSizeFactor = size
            self.updateBeamCoordinates()
            # self.autoZoom = autoZoom

    def setAutoZoom(self, switch):
        self.autoZoom = switch

    def getBeamSizeFactor(self):
        return self.beamSizeFactor

    def setBeamNumber(self, number):
        if number != self.beamNumber:
            self.beamNumber = number
            self.updateBeamCoordinates()

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

    def getProjectedBaselines(self):
        projectedBaselines = np.array(self.projectedBaselines)
        uvCoord = sv.rotateCoordinate(projectedBaselines, np.pi/2.0 -
                self.boreSightHorizontal[1], self.boreSightHorizontal[0])
        return np.concatenate((uvCoord, -uvCoord))

    def getBeamAxis(self):
        return self.beamAxis

    def getImageLength(self):
        return self.imageLength

    def getAntCoordinates(self):
        return self.antCoordinatesGEODET.tolist()

    def updateBeamCoordinates(self):
        halfLength = self.beamSizeFactor * 10
        interval = self.beamSizeFactor
        self.partialDFTGrid = self.createDFTGrid(self.gridNumOfDFT, halfLength, interval)


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
            zoomIn, density = 20, gridNum = 1000.0):

        halfGridNum = gridNum/2.
        # sidelength = density * zoomIn
        uMax = np.amax(np.abs(rotatedProjectedBaselines[:,0]))/waveLength
        vMax = np.amax(np.abs(rotatedProjectedBaselines[:,1]))/waveLength
        uvMax = uMax if uMax > vMax else vMax
        step = halfGridNum/(uvMax/waveLength)
        # imageLength = 1/(1/(halfGridNum/(uvMax/self.waveLength)))
        imageLength = step
        # windowLength = imageLength/gridNum*sidelength
        # print step

        return imageLength

    def partialDFT(self, partialDFTGrid, rotatedProjectedBaselines, waveLength,
            imageLength, zoomIn, density = 20, gridNum = 1000.0):

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
            fringeSum = fringeSum + np.exp(1j*np.pi*2*imagesCoord[1]*uv[0]/gridNum)*np.exp(1j*np.pi*2*imagesCoord[0]*uv[1]/gridNum)
        fringeSum = fringeSum.reshape(density,density)/(len(uvSamples))

        image = np.fft.fftshift(np.abs(fringeSum))

        return image

    def calculateBeamScaleFromBaselines(self, rotatedProjectedBaselines):
        baselineLengths = sv.distances(rotatedProjectedBaselines)
        baselineMax = np.amax(baselineLengths)
        baselineMin = np.amin(baselineLengths)
        indexOfMaximum = np.argmax(baselineLengths)
        maxBaselineVector = rotatedProjectedBaselines[indexOfMaximum]
        # rotate vector on a surface https://math.stackexchange.com/questions/1830695/
        perpendicularOfMaxBaselineVector = sv.projectedRotate(
                np.pi/2., 0, maxBaselineVector, np.pi/2.)
        perpendicularUnitVector = perpendicularOfMaxBaselineVector/sv.distances(
                perpendicularOfMaxBaselineVector)
        perpendicularBaselines = np.dot(rotatedProjectedBaselines, perpendicularUnitVector)
        perpendicularBaselineMax = np.amax(np.abs(perpendicularBaselines))

        minorAxis = 1.22*self.waveLength/baselineMax/2.
        majorAxis = 1.22*self.waveLength/perpendicularBaselineMax/2.

        return majorAxis, minorAxis


    def trackBorder(self, image, threshold = 0.3, density = 20, interpolatedLength = 800):

        interpolater = interpolate.interp2d(range(density), range(density), image ,kind='cubic')
        image = interpolater(np.linspace(0, density - 1, interpolatedLength),
                np.linspace(0, density - 1, interpolatedLength))
        # np.savetxt('interpolate', image)


        imageCenter = [interpolatedLength/2, interpolatedLength/2]
        rowIdx = int(imageCenter[0])
        colIdx = int(imageCenter[1])

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
        closestToCenter = 1
        closestToCenterIndex = []

        while state != State.end:
            # print rowIdx, colIdx, image[rowIdx, colIdx]
            if state == State.move:
                if colIdx != colBorder:
                    colIdx += colStep
                    if image[rowIdx, colIdx] < threshold:
                        border.append([rowIdx, colIdx])
                        state = State.addRow
                    else:
                        distToCenter = 1 - image[rowIdx, colIdx]
                        if distToCenter  < closestToCenter:
                            closestToCenter = distToCenter
                            closestToCenterIndex = [rowIdx, colIdx]
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

        return border, closestToCenterIndex, overstep


    def createContour(self, antennacoor, fileName='contour.png', minAlt=0):

        def getAltAziFromRADEC(beamCoordinates, LSTDeg, arrayRefereceLatitude):
            RA = np.deg2rad(beamCoordinates[:,0])
            DEC = np.deg2rad(beamCoordinates[:,1])
            LST = np.deg2rad(LSTDeg)

            altitude, azimuth = coord.convertEquatorialToHorizontal(
                    RA, DEC, LST, arrayRefereceLatitude)

            return altitude, azimuth

        self.antCoordinatesGEODET = np.array(antennacoor)

        antCoordinatesECEF = coord.convertGodeticToECEF(self.antCoordinatesGEODET)

        arrayRefereceECEF = coord.convertGodeticToECEF([self.arrayRefereceGEODET])[0]
        antCoordinatesENU = coord.convertECEFToENU(antCoordinatesECEF, arrayRefereceECEF, self.arrayRefereceGEODET)

        self.baselines = bs.createBaselines(antCoordinatesENU)

        LSTDeg =  coord.calculateLocalSiderealTime(self.observeTime, self.arrayRefereceGEODET[0])

        arrayRefereceLatitude = np.deg2rad(self.arrayRefereceGEODET[1])

        altitude, azimuth = getAltAziFromRADEC(np.array([self.boreSight]),
                LSTDeg, arrayRefereceLatitude)

        self.boreSightHorizontal = (azimuth[0], altitude[0])

        if abs(altitude[0]) < minAlt:
            return

        baselineNum = len(self.baselines)
        self.projectedBaselines = sv.projectedBaselines(
                np.array([altitude[0]]), np.array([azimuth[0]]), self.baselines)

        rotatedProjectedBaselines = sv.rotateCoordinate(self.projectedBaselines,
                np.pi/2.0 - altitude[0], azimuth[0])

        beamMajorAxisScale, beamMinorAxisScale = self.calculateBeamScaleFromBaselines(
                rotatedProjectedBaselines)

        density = 20
        gridNum = 1000.0

        imageLength = self.calculateImageLength(rotatedProjectedBaselines,
                self.waveLength, self.beamSizeFactor, density)

        newBeamSizeFactor = beamMajorAxisScale*1.3 / (imageLength/gridNum) / density
        # print "new factor", newBeamSizeFactor


        if newBeamSizeFactor < 1:
            newBeamSizeFactor = 1
        else:
            newBeamSizeFactor = int(newBeamSizeFactor)

        if baselineNum > 2 and self.autoZoom == True and abs(newBeamSizeFactor - self.beamSizeFactor) > 0:

            # self.beamSizeFactor = newBeamSizeFactor
            print "new beam factor:", newBeamSizeFactor
            self.setBeamSizeFactor(newBeamSizeFactor)

        sidelength = density * self.beamSizeFactor
        windowLength = imageLength/gridNum*sidelength

        image = self.partialDFT(self.partialDFTGrid, rotatedProjectedBaselines,
                self.waveLength, imageLength, self.beamSizeFactor, density)

        self.imageLength = windowLength
        bs.plotBeamContour3(image, self.boreSightHorizontal, windowLength,
                interpolation = self.interpolating)

        angle = 0
        if baselineNum > 2:
            interpolatedLength = 800
            threshold = 0.5
            border, closestToCenterIndex, overstep = self.trackBorder(
                    image, threshold, density, interpolatedLength)
            # np.savetxt('border', border)


            if len(border) < 10:
                majorAixs = np.rad2deg(beamMajorAxisScale)
                minorAixs = np.rad2deg(beamMajorAxisScale)
                angle = 0
                self.beamAxis = [majorAixs, minorAixs, angle]
                return

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
                # angleVert = np.arctan2(maxDistVector[0], maxDistVector[1])
                # angleHonr = np.pi/2. + np.arctan2(minDistVector[0], minDistVector[1])

                minorAixs  = np.rad2deg(windowLength/interpolatedLength*minDist)
                majorAixs  = np.rad2deg(windowLength/interpolatedLength*maxDist)
                # print np.rad2deg(majorAixs1), np.rad2deg(minorAixs1)
            else:
                angle = np.pi/2. + np.arctan2(minDistVector[0], minDistVector[1])


                majorAixs = np.rad2deg(beamMajorAxisScale)
                minorAixs  = np.rad2deg(windowLength/interpolatedLength*minDist)
                # self.beamAxis = [majorAixs, minorAixs, angle]
            self.beamAxis = [majorAixs, minorAixs, angle]

            print majorAixs, minorAixs, np.rad2deg(angle)



def sphericalToCartesian(theta, phi):
    x = np.sin(theta)*np.cos(phi)
    y = np.sin(theta)*np.sin(phi)
    z = np.cos(theta)

    return np.array([x, y, z]).T

def cartesianToSpherical(xyz):
    theta = np.arccos(xyz[:,2])
    phi = np.arctan2(xyz[:,1], xyz[:,0])

    return phi, theta

def rotateCoordinate(coordinates, theta, phi):

    rotationMatrix = np.array([
        [ np.cos(phi)*np.cos(theta),  -np.sin(phi),  np.cos(phi)*np.sin(theta)],
        [ np.sin(phi)*np.cos(theta),   np.cos(phi),  np.sin(phi)*np.sin(theta)],
        [-np.sin(theta),                    0,       np.cos(theta)]
    ])

    rotatedCoordinates = np.dot(coordinates, rotationMatrix.T.tolist())
    # print rotatedCoordinates
    return rotatedCoordinates

