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
    freqencies -- the central frequency of the observation

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
        """
        set the central frequency of the observation

        Keyword arguments:
        freqencies -- the central frequency of the observation

        """

        self.waveLengths = sol/np.array(freqencies)

    def setBeamPixelNum(self, num):
        """
        set the number of pixels of the image, default is 400,
        that will generate a image with a dimension of 20x20

        Keyword arguments:
        num -- the number of pixels of the image

        """

        self.observation.setBeamNumber(num)

    def setAntennas(self, antennas):
        """
        set a list of antenna.

        Keyword arguments:
        antennas -- a list of coordinates whose element is
            in the order of latitude(deg), longitude(deg), altitude(meter)

        """

        self.antennas = antennas

    def getBeamShape(self):
        """
        return the beamshape of current oservation parameters
        assuming the beam is roughly a ellipse.


        return:
        axisH, semi-major axis of the ellipse
        axisV, semi-minor axis of the ellipse
        angle, orientation of the ellipse
        """


        self.observation.createContour(self.antennas, 'contour.png')
        shape = self.observation.getBeamAxis()
        # print shape

        return shape

    def getHorizontal(self):
        """
        return the horizontal coordinate of the pointing.
        this method  MUST be run after calling getBeamShape method


        return:
        azimuth, in degree
        elevation, in degree
        """


        return self.observation.getHorizontal()



class tilesim:
    """
    Class for generation of  tiling

    Keyword arguments:
    axisH -- semi-major axis of the beam of a ellipse shape in degree
    axisV -- semi-minor axis of the beam of a ellipse shape in degree
    angle -- orientation of the beam of a ellipse shape in radian
    overlap -- how much overlap between two beams, range in (0, 1)
    beamNum -- how many beams are there in the tiling

    """


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


    def setBeamNumber(self, num):
        """
        set the beam number of the tiling

        Keyword arguments:
        num -- beam number of the tiling

        """

        self.beamNum = num

    def setBeamShape(self, axisH, axisV, angle):
        """
        set the shape of the beam

        Keyword arguments:
        axisH -- semi-major axis of the beam of a ellipse shape in degree
        axisV -- semi-minor axis of the beam of a ellipse shape in degree
        angle -- orientation of the beam of a ellipse shape in radian

        """

        self.axisH, self.axisV, self.angle = axisH, axisV, angle

    def getTilingRadius(self):
        """
        return the radius of the generated tiling.
        This method SHOULD be called after calling getTiling()

        return:
        radius -- the radius of the generated tiling in Degree

        """

        return self.tilingRadius

    def inverseGaussian(self, p, mu, sigma):
        """
        return the offset to the Gaussian center given the Gaussian informaton

        Keyword arguments:
        p -- the Gaussian probabilty
        mu -- the center of the Gaussian
        sigma -- the deviation  of the Gaussian

        return:
        offset -- the offset to the Gaussian center given the parameters

        """
        x = np.sqrt(-2.*sigma**2.*np.log(p))+mu
        return x

    def getTiling(self):
        """
        return the tiling.

        return:
        tiling -- tiling coordinates in a list of [RA, DEC] pairs in degree

        """

        self.sigmaH = self.axisH * (2./2.3556)
        self.sigmaV = self.axisV * (2./2.3556)

        self.widthH = self.inverseGaussian(self.overlapPoint, 0, self.sigmaH)
        self.widthV = self.inverseGaussian(self.overlapPoint, 0, self.sigmaV)


        self.tiling, self.tilingRadius = ellipseCompact(
                self.beamNum, self.widthH, self.widthV, self.angle, 10)

        return self.tiling

    def getTilingWithinRadius(self):
        """
        return the tiling inside a specified region. this method will ignore
        the beamNum parameter.
        the tiling region MUST be specified first using setTilingRadius method

        return:
        tiling -- tiling coordinates in a list of [RA, DEC] pairs in degree

        """

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
        """
        plot the tiling pattern as sky temperature with specified file name
        This method SHOULD be called after calling getTiling()
        The format and directory can be specified in the file name such as
        "plot/pattern.png" or "pattern.pdf"

        """

        overlapCounter = calculateBeamOverlaps(
                self.coordinates, self.tilingRadius,
                self.axisH, self.axisV, self.angle, "heater", fileName)

    def plotTiling(self, fileName):
        """
        plot the tiling pattern with specified file name.
        This method SHOULD be called after calling getTiling()

        Keyword arguments:
        filename --  filename of the plot, the format and directory can be
        specified in the file name such as  "plot/pattern.png" or "pattern.pdf"

        """

        plotPackedBeam(self.tiling.T, self.angle, self.widthH,self.widthV,
                self.tilingRadius, fileName=fileName)

    def calculateOverlap(self, mode, fileName):
        """
        calculate overlap of the tiling pattern.

        Keyword arguments:
        mode -- mode of the calculation,
                "heater" will calculate the tiling pattern as sky temperature,
                "counter" will calculate the counts of the overlap regions,
                          non-overlap regions and empty regions.
                "both" will calculate both sky temperature and counter
        filename --  filename of the plot, the format and directory can be
                specified in the file name such as  "plot/pattern.png"

        return:
        overlapCounter -- counter when choose "counter" mode
        overlapHeater -- heater when choose "heater" mode
        overlapCounter, overlap heater -- both result will returenwhen choose
                       "both" mode


        """
        overlapCounter = calculateBeamOverlaps(
                self.tiling, self.tilingRadius,
                self.axisH, self.axisV, self.angle, mode, fileName)

        return overlapCounter


    def setTilingRadius(self, radius):
        """
        set the radius of the tiling before calling getTilingWithinRadius

        Keyword arguments:
        radius -- radius of the tiling

        """

        self.tilingRadius = radius


    def setOverlapPoint(self, overlap):
        """
        set how much overlap between two beams.

        Keyword arguments:
        radius -- float number range in (0, 1)

        """

        self.overlapPoint = overlap

class DelayPolynomial:

    """
    Class for generation of  delay polynomial

    Keyword arguments:
    antennas -- a list of antenna objects or coordinate in csv format
    timestamp -- the observation time in datatime object or epoch seconds
    duration -- the duration in which the polynomial is calcuated
    frequencies -- a list of frequencies on which the polynomail is calculated
    reference -- the reference antenna for delay calculation

    """


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