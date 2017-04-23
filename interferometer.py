#! /usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

import steervec as sv
import coordinate as coord
import beamShape as bs
import createBeam as cb


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
        self.beamSizeFactor = 1
        self.beamCoordinates = []
        self.beamNumber = 400
        self.beamSize = 1.22*self.waveLength/13.5
        self.interpolating = True
        self.amplitude = []

    def setInterpolating(self, state):
        if state == True:
            self.interpolating = True
        else:
            self.interpolating = False

    def setObserveTime(self, dateTime):
        self.observeTime = dateTime

    def setBeamSizeFactor(self, size):
        if size != self.beamSizeFactor:
            self.beamSizeFactor = size
            self.updateBeamCoordinates()

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

    def getHorizontal(self):
        return self.boreSightHorizontal

    def updateBeamCoordinates(self):
        beamCoordinates, subBeamRadius = cb.optimizeGrid(
                self.beamNumber, np.rad2deg(self.beamSize/self.beamSizeFactor)/2.,
                cb.recGrid, 50, self.boreSight)
        self.beamCoordinates = np.array(beamCoordinates)


    def createContour(self, antennacoor):

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
        antCoordinatesGEODET = np.array(antennacoor)

        antCoordinatesECEF = coord.convertGodeticToECEF(antCoordinatesGEODET)

        arrayRefereceECEF = coord.convertGodeticToECEF([self.arrayRefereceGEODET])[0]
        antCoordinatesENU = coord.convertECEFToENU(antCoordinatesECEF, arrayRefereceECEF, self.arrayRefereceGEODET)

        self.baselines = bs.createBaselines(antCoordinatesENU)

        LSTDeg =  coord.calculateLocalSiderealTime(self.observeTime, self.arrayRefereceGEODET[0])

        arrayRefereceLatitude = np.deg2rad(self.arrayRefereceGEODET[1])

        # RA = np.deg2rad(self.beamCoordinates[:,0])
        # DEC = np.deg2rad(self.beamCoordinates[:,1])
        # LST = np.deg2rad(LSTDeg)
        # longitude = np.deg2rad(antCoordinatesGEODET[:,0])
        # latitude = np.deg2rad(antCoordinatesGEODET[:,1])

        altitude, azimuth = getAltAziFromRADEC(self.beamCoordinates,
                LSTDeg, arrayRefereceLatitude)




        projectedBaselines = sv.projectedBaselines(altitude, azimuth, self.baselines)


        # np.savetxt('baselines0421', projectedBaselines)
        maxBaseline = np.amax(np.abs(projectedBaselines))
        requireBeamSizeFactor = self.beamSize/(1.22*self.waveLength/maxBaseline)
        print "req beamfactor: ", requireBeamSizeFactor, "current beamFactor: ", self.beamSizeFactor
        if abs(requireBeamSizeFactor - self.beamSizeFactor) > 1:
            self.setBeamSizeFactor(round(requireBeamSizeFactor))
            altitude, azimuth = getAltAziFromRADEC(self.beamCoordinates,
                LSTDeg, arrayRefereceLatitude)
            projectedBaselines = sv.projectedBaselines(altitude, azimuth, self.baselines)

        waveNumbers = sv.waveNumber(altitude, azimuth, self.waveLength)

        weights = sv.weightVector(waveNumbers, self.baselines)



        self.boreSightHorizontal = (azimuth[0], altitude[0])
        self.amplitude = bs.fringePlot(self.beamCoordinates, weights, self.baselines,
                self.boreSight, np.rad2deg(self.beamSize), self.interpolating)

        print "max baseline: ", maxBaseline

        # np.savetxt('beamCoor', np.c_[altitude, azimuth])
        # np.savetxt('antennaCoor', self.baselines)
        # np.savetxt('waveNumberPy', waveNumbers )
        # np.savetxt('baselines0421', projectedBaselines)
        # np.savetxt('weights0421', weights)

