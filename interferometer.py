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
        self.beamSizeFactor = 1
        self.beamCoordinates = []
        self.beamNumber = 400
        self.beamSize = 1.22*self.waveLength/13.5
        self.interpolating = True
        self.amplitude = []
        self.autoZoom = True

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

    def getAntCoordinates(self):
        return self.antCoordinatesGEODET.tolist()

    def updateBeamCoordinates(self):
        beamCoordinates, subBeamRadius = cb.optimizeGrid(
                self.beamNumber, np.rad2deg(self.beamSize/self.beamSizeFactor)/2.,
                cb.recGrid, 50, np.array(self.boreSight))
        self.beamCoordinates = np.array(beamCoordinates)


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

        altitude, azimuth = getAltAziFromRADEC(self.beamCoordinates,
                LSTDeg, arrayRefereceLatitude)

        if abs(altitude[0]) < minAlt:
            self.boreSightHorizontal = (azimuth[0], altitude[0])
            return

        # print 'autoZoom: ', self.autoZoom
        self.baselineIndex, self.projectedBaselines = sv.projectedBaselines(
                np.array([altitude[0]]), np.array([azimuth[0]]), self.baselines)

        drop, projectedEastNorth = sv.projectedBaselines(
                np.array([altitude[0]]), np.array([azimuth[0]]),
                [[1,0,0], [0,1,0]])

        # print projectedEastNorth

        # print sphericalToCartesian(np.pi/2.-altitude[0], np.pi/2-azimuth[0])

        if self.autoZoom == True:
            # np.savetxt('baselines0421', projectedBaselines)
            baselineLengths = sv.distances(self.projectedBaselines)
            baselineMax = np.amax(baselineLengths)
            indexOfMaximum = np.argmax(baselineLengths)
            maxBaselineVector = self.projectedBaselines[indexOfMaximum]
            maxBaselineVectorOriginal = self.baselineIndex[indexOfMaximum]
            # rotate vector on a surface https://math.stackexchange.com/questions/1830695/
            perpendicularOfMaxBaselineVector = sv.projectedRotate(
                    altitude[0], azimuth[0], maxBaselineVector, np.pi/2.)
            perpendicularUnitVector = perpendicularOfMaxBaselineVector/sv.distances(perpendicularOfMaxBaselineVector)
            perpendicularBaselines = np.dot(self.projectedBaselines, perpendicularUnitVector)
            perpendicularBaselineMax = np.amax(np.abs(perpendicularBaselines))
            # maxBaseline = np.amax(baselineLengths)
            # minBaseline = np.amin(baselineLengths)
            baseNum = len(antCoordinatesENU)
            # factor = (0.140845*baseNum + 28.7887)/(baseNum + 11.6056)
            # Hyperbola function
            factor = (50.5223 - 0.382166*baseNum)/(baseNum + 22.879)
            factor = factor * (1 - abs(altitude[0])/(np.pi/2.)*0.5)
            # factor = 1
            # print maxBaselineVector , perpendicularUnitVector
            # print self.baselines
            # print self.baselines
            # print baselineMax , perpendicularBaselineMax
            # projectedZenithVector = sv.projectedBaselines(
                    # np.array([altitude[0]]), np.array([azimuth[0]]), [[0,0,1],])[0]
            # zenithUnitVector = projectedZenithVector/sv.distances(projectedZenithVector)
            # maxBaselineUnitVector = maxBaselineVector/baselineMax

            # rotatedMaxBaselineVector = rotateCoordinate([maxBaselineVector,],
                    # (np.pi/2. - altitude[0]), (np.pi/2. - azimuth[0]))[0]

            # maxBaselineUnitVector = rotatedMaxBaselineVector/sv.distances(rotatedMaxBaselineVector)

            # angle = np.arccos(np.dot(zenithUnitVector, maxBaselineUnitVector)/1)
            phi, theta = cartesianToSpherical(np.array([maxBaselineVector/sv.distances(maxBaselineVector)]))
            perpendicularRA, perpendicularDEC = coord.convertHorizontalToEquatorial(
                    np.pi/2. - phi, np.pi/2. - theta, np.deg2rad(LSTDeg), arrayRefereceLatitude)
            # print perpendicularRA, perpendicularDEC

            # angle = np.arctan2(np.rad2deg(perpendicularDEC[0]), np.rad2deg(perpendicularRA[0]))

            # rotatedUV = sv.rotateCoordinate(maxBaselineVector, -np.pi/2.0 + altitude[0], -np.pi/2.0 + azimuth[0])
            # rotatedUV = sv.rotateCoordinate(maxBaselineVector, np.pi/2.0 - altitude[0], azimuth[0])
            print "rotate:"
            angle = np.pi/2.0 + np.arccos(np.dot(maxBaselineVector, projectedEastNorth[0])/baselineMax*sv.distances(projectedEastNorth[0]))
            # angle =  np.arctan2(rotatedUV [0], rotatedUV [1])
            print np.rad2deg(angle)
            self.beamAxis = [np.rad2deg(1.22*self.waveLength/perpendicularBaselineMax/2.),
                np.rad2deg(1.22*self.waveLength/baselineMax/2.), np.rad2deg(angle)]
            # print zenithUnitVector, perpendicularUnitVector
            # print 'aixsLengthCal: ', self.beamAxis, np.rad2deg(angle)
            requireBeamSizeFactor = factor *  self.beamSize/(1.22*self.waveLength/perpendicularBaselineMax)
            if abs(requireBeamSizeFactor - self.beamSizeFactor) > 1:
                rounedFactor = round(requireBeamSizeFactor)
                self.setBeamSizeFactor(rounedFactor if rounedFactor != 0 else 1)
                altitude, azimuth = getAltAziFromRADEC(self.beamCoordinates,
                    LSTDeg, arrayRefereceLatitude)
                self.baselineIndex, self.projectedBaselines =sv.projectedBaselines(
                        np.array([altitude[0]]), np.array([azimuth[0]]), self.baselines)

        else:
            self.autoZoom = True

        # print self.projectedBaselines
        waveNumbers = sv.waveNumber(altitude, azimuth, self.waveLength, False)

        weights = sv.weightVector(waveNumbers, self.baselines)
        # maxOriginalBaseline = [maxBaselineVectorOriginal,]
        # weightMaxOriginal = sv.weightVector(waveNumbers, maxOriginalBaseline )



        self.boreSightHorizontal = (azimuth[0], altitude[0])
        self.amplitude = bs.fringePlot(np.column_stack((np.rad2deg(altitude), np.rad2deg(azimuth))), weights, self.baselines,
        # self.amplitude = bs.fringePlot(self.beamCoordinates, weights, self.baselines,
                self.boreSight, np.rad2deg(self.beamSize),
                self.interpolating, fileName=fileName)

        # originalAmplitude = bs.fringePlot(self.beamCoordinates,
                # weightMaxOriginal,maxOriginalBaseline ,
                # self.boreSight, np.rad2deg(self.beamSize),
                # self.interpolating, fileName='maxBaseline.png')

        # print maxBaselineVector
        # rotatedUV = sv.rotateCoordinate(maxBaselineVector, -np.pi/2.0 + altitude[0], -np.pi/2.0 + azimuth[0])
        # print rotatedUV
        # print np.rad2deg(np.arctan2(maxBaselineVector[1], maxBaselineVector[0]))

        # np.savetxt('beamOrignal', originalAmplitude)
        # getCentralLine(self.beamCoordinates, originalAmplitude)


        # print "max baseline: ", maxBaseline

        # np.savetxt('beamCoor', np.c_[altitude, azimuth])
        # np.savetxt('antennaCoor', self.baselines)
        # np.savetxt('waveNumberPy', waveNumbers )
        # np.savetxt('baselines0421', projectedBaselines)
        # np.savetxt('weights0421', weights)


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

