import numpy as np
import datetime

import katpoint
import coordinate as coord
from interferometer import InterferometryObservation
from tile import ellipseCompact, ellipseGrid
from plot import plotPackedBeam
from beamshape import calculateBeamOverlaps

class psfsim:
    """
    Class for simulation of beam shape.

    Keyword arguments:
    antenna -- a list of antenna coordinates whose element is
        in the order of latitude(deg), longitude(deg), altitude(meter)
    source -- the boresight location or telescope pointing
        in the order of RA(deg), DEC(deg)
    time -- the observation time in datatime object or epoch seconds
    freqencies -- the central freqency of the observation

    """

    arrayRefereceGEODET = (-30.71106, 21.44389, 1035)
    '''speed Of Light'''
    sol = 299792458
    # self.waveLengths = sol/freqencies

    def __init__(self, antennas, source, time, freqencies):
        """
        constructor of the psfsim class.

        """

        self.pixels = 400
        self.boreSight = (21.44389, -30.71106)
        self.observation = None
        self.antennas = []

        waveLengths = float(self.sol)/np.array(freqencies)
        if waveLengths.shape == (1,): waveLengths = waveLengths[0]

        observation = InterferometryObservation(self.arrayRefereceGEODET,
                            None, waveLengths)
        observation.setBoreSight(self.boreSight)
        observation.setBeamSizeFactor(1)
        observation.setBeamNumber(self.pixels)
        observation.setInterpolating(True)
        observation.setAutoZoom(True)
        observation.setInputType(InterferometryObservation.equatorialInput)

        observation.setObserveTime(time)
        observation.setBoreSight(source)

        self.observation = observation
        self.antennas = antennas

    def setSource(self, source):
        """
        set the boreSight location or the pointing of the observation

        Keyword arguments:
        source -- -- the boresight location or telescope pointing
            in the order of RA(deg), DEC(deg)

        """

        self.observation.setBoreSight(source)

    def setTime(self, observeTime):
        """
        set the observation time

        Keyword arguments:
        time -- the observation time in datatime object or epoch seconds

        """

        self.observation.setObserveTime(observeTime)

    def setFreqencies(self, freqencies):
        self.waveLengths = sol/np.array(freqencies)

    def setBeamPixelNum(self, num):
        self.pixels = num

    def setAntennas(self, antennas):
        self.antenn = antennas

    def getBeamShape(self):

        self.observation.createContour(self.antennas, 'contour.png')
        shape = self.observation.getBeamAxis()
        print shape

        return shape

    def getHorizontal(self):

        return self.observation.getHorizontal()



class tilesim:


    def __init__(self, axisH, axisV, angle, beamNum, overlap):
        self.beamNum = beamNum
        self.axisH = axisH
        self.axisV = axisV
        self.angle = angle
        self.overlapPoint = overlap
        self.tilingRadius = 0
        self.tiling = []
        self.widthH = 0
        self.widthV = 0

    # def __init__(self, args):
        # pass

    def setBeamNumber(self, num):
        self.beamNum = num

    def setBeamShape(self, axisH, axisV, angle):
        self.axisH, self.axisV, self.angle = axisH, axisV, angle

    def getTilingRadius(self):
        return self.tilingRadius

    def inverseGaussian(self, p, mu, sigma):
        x = np.sqrt(-2.*sigma**2.*np.log(p))+mu
        return x

    def getTiling(self):
        self.sigmaH = self.axisH * (2./2.3556)
        self.sigmaV = self.axisV * (2./2.3556)

        widthH = self.inverseGaussian(self.overlapPoint, 0, self.sigmaH)
        widthV = self.inverseGaussian(self.overlapPoint, 0, self.sigmaV)

        self.coordinates, self.tilingRadius = ellipseCompact(
                self.beamNum, widthH, widthV, self.angle, 10)

        return self.coordinates

    def getTilingWithinRadius(self):
        self.sigmaH = self.axisH * (2./2.3556)
        self.sigmaV = self.axisV * (2./2.3556)

        widthH = self.inverseGaussian(self.overlapPoint, 0, self.sigmaH)
        widthV = self.inverseGaussian(self.overlapPoint, 0, self.sigmaV)
        self.widthH = widthH
        self.widthV = widthV


        tilingWithinRadius = ellipseGrid(
                self.tilingRadius, widthH, widthV, self.angle)
        self.tiling = tilingWithinRadius

        return tilingWithinRadius


    def plotSkyBrightness(self, fileName):
        overlapCounter = calculateBeamOverlaps(
                self.coordinates, self.tilingRadius,
                self.sigmaH,self.sigmaV, self.angle, fileName)

    def plotTiling(self, fileName):
        plotPackedBeam(self.tiling.T, self.angle, self.widthH,self.widthV,
                self.tilingRadius, fileName=fileName)


    def setTilingRadius(self, radius):
        self.tilingRadius = radius


    def setOverlapPoint(self, overlap):
        self.overlapPoint = overlap

class DelayPolynomial:


    # def __init__(self, args):
        # pass

    def __init__(self, antennas, targets,
            timestamp, duration, freqencies, reference):

        self.antennas = antennas
        self.timestamp = self.checkTime(timestamp)
        self.targets = targets
        self.freqencies = freqencies
        self.reference = reference
        self.duration = duration

    def checkTargets(self, targets):
        if isinstance(targets[0], katpoint.Target):
            return targets
        else:
            return makeKapointTarget(targets)

    def setFreqencies(self, freqencies):
        self.freqencies = freqencies

    def setAntennas(self, antennas):
        self.antennas = antennas

    def setTimestamp(self, timestamp):
        self.timestamp = self.checkTime(timestamp)

    def setTimeDuration(self, duration):
        self.duration = duration

    def setTargets(self, targets):
        self.targets = targets

    def setReference(self, reference):
        self.reference = reference

    def checkTime(self, time):
        if type(time) != int and type(time) != float:
            return coord.datetimeToEpoch(time)
        else:
            return time

    def dictToOrderList(self, dictObj):
        orderList = []
        for key in sorted(dictObj.iterkeys()):
            orderList.append(dictObj[key])

        return orderList

    @staticmethod
    def makeKapointTarget(sources):
        targets = []
        for source in sources:
            targetString = ",".join(['radec',
                            coord.angleToHour(source[0]),
                            coord.angleToDEC(source[1])])
            targets.append(katpoint.Target(targetString))

        return targets

    @staticmethod
    def makeAntenna(antennaString):
        antennaGeo = []
        antennaKat = []

        for antenna in antennaString:
            antkat = katpoint.Antenna(antenna)
            geo = antkat.position_wgs84
            antennaKat.append(antkat)
            antennaGeo.append(
                    [np.rad2deg(geo[0]), np.rad2deg(geo[1]), geo[2]])

        return antennaKat, antennaGeo



    def getPolynomial(self):
        antennaObjectList = self.antennas
        targets = self.targets
        ref = self.reference
        timestamp = (self.timestamp, self.timestamp + self.duration)
        sky_centre_freq = self.freqencies

        freqArray = []
        for freq in sky_centre_freq:
            targetArray = []
            for target in targets:
                dc = katpoint.DelayCorrection(antennaObjectList , ref, freq)
                delay, phase = dc.corrections(target, timestamp)
                phaseArray = np.array(self.dictToOrderList(phase))
                targetArray.append(phaseArray[::2][:,0,:])
            targetArray = np.array(targetArray)
            targetArray = targetArray - targetArray[0, :, :]
            freqArray.append(targetArray)

        freqArray = np.array(freqArray)

        # freq, beam, antenna [,pol], time, (phase, rate)

        return freqArray


