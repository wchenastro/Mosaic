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

    def setBeamNumber(self, number):
        if number != self.beamNumber:
            self.beamNumber = number
            self.updateBeamCoordinates()

    def getBaselines(self):
        return self.baselines

    def setBoreSight(self, beamBoreSight):
        self.boreSight = beamBoreSight
        self.updateBeamCoordinates()

    def getHorizontal(self):
        return self.boreSightHorizontal

    def updateBeamCoordinates(self):
        beamCoordinates, subBeamRadius = cb.optimizeGrid(
                self.beamNumber, self.beamSize/self.beamSizeFactor/2.,
                cb.recGrid, 50, self.boreSight)
        self.beamCoordinates = np.array(beamCoordinates)

    def createContour(self, antennacoor):

        antennasCoordinateFile = 'antennacoor'
        beamCoordinateFile = 'inCoord'

        # beamCoordinates = coord.readCoordinates(beamCoordinateFile)
        antCoordinatesGEODET = np.array(antennacoor)

        antCoordinatesECEF = coord.convertGodeticToECEF(antCoordinatesGEODET)

        arrayRefereceECEF = coord.convertGodeticToECEF([self.arrayRefereceGEODET])[0]
        antCoordinatesENU = coord.convertECEFToENU(antCoordinatesECEF, arrayRefereceECEF, self.arrayRefereceGEODET)

        self.baselines = bs.createBaselines(antCoordinatesENU)

        LSTDeg =  coord.calculateLocalSiderealTime(self.observeTime, self.arrayRefereceGEODET[0])

        RA = np.deg2rad(self.beamCoordinates[:,0])
        DEC = np.deg2rad(self.beamCoordinates[:,1])
        LST = np.deg2rad(LSTDeg)
        # longitude = np.deg2rad(antCoordinatesGEODET[:,0])
        # latitude = np.deg2rad(antCoordinatesGEODET[:,1])
        arrayRefereceLatitude = np.deg2rad(self.arrayRefereceGEODET[1])

        altitude, azimuth = coord.convertEquatorialToHorizontal(
                RA, DEC, LST, arrayRefereceLatitude)

        self.boreSightHorizontal = (azimuth[0], altitude[0])

        waveNumbers = sv.waveNumber(altitude, azimuth, self.waveLength)

        weights = sv.weightVector(waveNumbers, self.baselines)

        bs.fringePlot(self.beamCoordinates, weights, self.baselines, self.boreSight, self.beamSize, self.interpolating)

