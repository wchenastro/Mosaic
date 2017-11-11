#! /usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

import steervec as sv
import coordinate as coord
import beamShape as bs
import createBeam as cb

import inspect


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
        self.beamSizeFactor = 0.1
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
        # print inspect.stack()[1][3]
        if size != self.beamSizeFactor:
            self.beamSizeFactor = size
            self.updateBeamCoordinates()
            self.autoZoom = autoZoom

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

    # def getProjectedBaselines(self):
        # return np.array(self.projectedBaselines)

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
        # beamCoordinates, subBeamRadius = cb.optimizeGrid(
                # self.beamNumber, np.rad2deg(self.beamSize/self.beamSizeFactor)/2.,
                # cb.recGrid, 50, np.array(self.boreSight))
        # self.beamCoordinates = np.array(beamCoordinates)
        self.partialDFTGrid = self.createDFTGrid(self.gridNumOfDFT, 10, 1)


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

    def partialDFT(self):
        pass


    def createContour(self, antennacoor, fileName='contour.png', minAlt=0):

        def getAltAziFromRADEC(beamCoordinates, LSTDeg, arrayRefereceLatitude):
            RA = np.deg2rad(beamCoordinates[:,0])
            DEC = np.deg2rad(beamCoordinates[:,1])
            LST = np.deg2rad(LSTDeg)

            altitude, azimuth = coord.convertEquatorialToHorizontal(
                    RA, DEC, LST, arrayRefereceLatitude)

            return altitude, azimuth

        # antennasCoordinateFile = 'antennacoor'
        # beamCoordinateFile = 'inCoord'

        # beamCoordinates = coord.readCoordinates(beamCoordinateFile)
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

        self.projectedBaselines = sv.projectedBaselines(
                np.array([altitude[0]]), np.array([azimuth[0]]), self.baselines)

        rotatedProjectedBaselines = sv.rotateCoordinate(self.projectedBaselines,
                np.pi/2.0 - altitude[0], azimuth[0])

        halfGridNum = 500.0
        gridNum = halfGridNum*2.0
        # pixelNum = imageLength * imageLength
        # uvGrids = np.zeros((gridNum, gridNum), dtype=np.int)
        uMax = np.amax(np.abs(rotatedProjectedBaselines[:,0]))/self.waveLength
        vMax = np.amax(np.abs(rotatedProjectedBaselines[:,1]))/self.waveLength
        uvMax = uMax if uMax > vMax else vMax
        step = halfGridNum/(uvMax/self.waveLength)
        # imageLength = 1/(1/(halfGridNum/(uvMax/self.waveLength)))
        imageLength = step
        windowLength = imageLength/1000.0*20.0
        self.imageLength = windowLength
        # print step
        uvSamples = []
        for baseline in rotatedProjectedBaselines:
            # print baseline
            uSlot = int((baseline[0]/self.waveLength*step + halfGridNum - 1))
            vSlot = int((halfGridNum - baseline[1]/self.waveLength*step - 1))
            # print uSlot, vSlot

            uvSamples.append([uSlot, vSlot])
            uvSamples.append([gridNum - uSlot - 1, gridNum - vSlot - 1])
            # uvGrids[vSlot][uSlot] = 1
            # uvGrids[-vSlot-1][-uSlot-1] = 1


        imagesCoord = self.partialDFTGrid
        imagesValue = []
        # for coord in imagesCoord:
        fringeSum = np.zeros(20*20)
        for uv in uvSamples:
            fringeSum = fringeSum +  np.exp(1j*np.pi*2*imagesCoord[1]*uv[0]/gridNum)*np.exp(1j*np.pi*2*imagesCoord[0]*uv[1]/gridNum)
        fringeSum = fringeSum.reshape(20,20)/(len(uvSamples))


        # psf = np.fft.ifft2(uvGrids)

        # np.savetxt('psf', np.abs(psf))
        # np.savetxt('psf2', np.abs(fringeSum))
        # np.savetxt('vuGrid', uvGrids, fmt='%.0f')

        # print self.boreSightHorizontal
        bs.plotBeamContour3(np.fft.fftshift(np.abs(fringeSum)), self.boreSightHorizontal, windowLength)



        # projectedEastNorth = sv.projectedBaselines(
                # np.array([altitude[0]]), np.array([azimuth[0]]),
                # [[1,0,0], [0,1,0]])

        baselineLengths = sv.distances(rotatedProjectedBaselines)
        baselineMax = np.amax(baselineLengths)
        baselineMin = np.amin(baselineLengths)
        indexOfMaximum = np.argmax(baselineLengths)
        maxBaselineVector = rotatedProjectedBaselines[indexOfMaximum]
        # rotate vector on a surface https://math.stackexchange.com/questions/1830695/
        perpendicularOfMaxBaselineVector = sv.projectedRotate(
                np.pi/2., 0, maxBaselineVector, np.pi/2.)
                # altitude[0], azimuth[0], maxBaselineVector, np.pi/2.)
        perpendicularUnitVector = perpendicularOfMaxBaselineVector/sv.distances(perpendicularOfMaxBaselineVector)
        perpendicularBaselines = np.dot(rotatedProjectedBaselines, perpendicularUnitVector)
        perpendicularBaselineMax = np.amax(np.abs(perpendicularBaselines))
        baseNum = len(antCoordinatesENU)
        # factor = (0.140845*baseNum + 28.7887)/(baseNum + 11.6056)
        # Hyperbola function
        # factor = (50.5223 - 0.382166*baseNum)/(baseNum + 22.879)
        # factor = factor * (1 - abs(altitude[0])/(np.pi/2.)*0.5)

        # angle = np.arccos(np.dot(zenithUnitVector, maxBaselineUnitVector)/1)
        # phi, theta = cartesianToSpherical(np.array([maxBaselineVector/sv.distances(maxBaselineVector)]))
        # perpendicularRA, perpendicularDEC = coord.convertHorizontalToEquatorial(
                # np.pi/2. - phi, np.pi/2. - theta, np.deg2rad(LSTDeg), arrayRefereceLatitude)
        # print perpendicularRA, perpendicularDEC

        # angle = np.arctan2(np.rad2deg(perpendicularDEC[0]), np.rad2deg(perpendicularRA[0]))

        # rotatedUV = sv.rotateCoordinate(maxBaselineVector, -np.pi/2.0 + altitude[0], -np.pi/2.0 + azimuth[0])
        # rotatedUV = sv.rotateCoordinate(maxBaselineVector, np.pi/2.0 - altitude[0], azimuth[0])
        print baselineMin, baselineMax
        # angle = np.pi/2.0 + np.arccos(np.dot(maxBaselineVector, projectedEastNorth[0])/baselineMax*sv.distances(projectedEastNorth[0]))

        vectorSum = np.sum(rotatedProjectedBaselines, axis=0)

        # angle =  np.arctan2(maxBaselineVector[1], maxBaselineVector[0])
        # minorAixs = np.rad2deg(1.22*self.waveLength/baselineMax/2.)
        # majorAixs = np.rad2deg(1.22*self.waveLength/perpendicularBaselineMax/2.)
        # self.beamAxis = [majorAixs, minorAixs, np.rad2deg(angle)]

        angle = np.pi/2.0 +  np.arctan2(vectorSum[1], vectorSum[0])
        summedLength = sv.distances(vectorSum)
        print summedLength
        minorAixs = np.rad2deg(1.22*self.waveLength/summedLength/2.)
        majorAixs = np.rad2deg(1.22*self.waveLength/(summedLength*2)/2.)
        self.beamAxis = [majorAixs, minorAixs, np.rad2deg(angle)]

        print minorAixs, majorAixs, np.rad2deg(angle)



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

