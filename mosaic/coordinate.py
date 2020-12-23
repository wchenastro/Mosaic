#!/usr/bin/env python

import datetime
import numpy as np
from astropy import wcs
from astropy.io import fits
import nvector as nv

from astropy import units as u
from astropy.coordinates import SkyCoord

class Antenna(object):
    def __init__(self, name, geo=None):
        self.name = name
        self.geo = geo
        self.enu = None

        grs80 = nv.FrameE(name='GRS80')
        lat, lon, height = geo
        geoPoint = grs80.GeoPoint(latitude=lat,
            longitude=lon, z=-height, degrees=True)
        self.ecef = geoPoint.to_ecef_vector().pvector.ravel()

class Baseline(object):
    def __init__(self, ant0, ant1):
        self.ant0 = ant0
        self.ant1 = ant1
        self.enu = np.array(ant0.enu) - np.array(ant1.enu)


class Array(object):
    def __init__(self, name, antennas, reference):
        self.name = name
        self.antennas = antennas
        self.reference = reference
        self.baselines = None

        ecefCoordinates = [antenna.ecef for antenna in antennas]
        enuCoordiantes = convertECEFToENU(ecefCoordinates, reference.ecef, reference.geo)
        for i in range(len(self.antennas)): self.antennas[i].enu = enuCoordiantes[i]

        self.baselines = self.createBaselines(self.antennas)

    def createBaselines(self, antennas):
         pairs = []
         index = 1
         for i in antennas:
             for j in antennas[index:]:
                 pairs.append(Baseline(i,j))
             index += 1
         return pairs

    def getBaselines(self):
        return self.baselines

    def getAntennas(self):
        return self.antennas

    def getRotatedProjectedBaselines(self, boresight):

        localHourAngle = (np.deg2rad(boresight.localSidereTime)
            - np.deg2rad(boresight.equatorial[0]))

        baselineCoordiantes = [baseline.enu for baseline in self.baselines]

        rotatedENU = rotateENUToEquatorialPlane(baselineCoordiantes,
            np.deg2rad(self.reference.geo[0]), boresight.horizontal[0],
            boresight.horizontal[1])

        rotatedProjectedBaselines = projectBaselines(rotatedENU,
            localHourAngle, np.deg2rad(boresight.equatorial[1]))

        return  rotatedProjectedBaselines

class Boresight(object):

    EquatorialFrame = 0
    HorizontalFrame = 1

    def __init__(self, name, boresight, time, arrayReference, frame):
        self.name = name
        self.time = time
        self.arrayreference = arrayReference

        arrayRefereceLatitude = np.deg2rad(arrayReference.geo[0])
        LSTDeg =  calculateLocalSiderealTime(time, arrayReference.geo[1])


        if frame == Boresight.EquatorialFrame:
            altitude, azimuth = self.getAltAziFromRADEC([boresight],
                LSTDeg, arrayRefereceLatitude)
            self.horizontal = (azimuth[0], altitude[0])
            self.equatorial = boresight
        elif frame == Boresight.HorizontalFrame:
            azimuth, altitude  = np.deg2rad(boresight)
            RA, DEC = convertHorizontalToEquatorial(azimuth, altitude,
                    np.deg2rad(LSTDeg), arrayRefereceLatitude)
            self.equatorial = np.rad2deg([RA, DEC])
            self.horizontal = np.deg2rad(boresight)
            #print("source in RA, DEC is %f, %f" % (self.boreSight[0], self.boreSight[1]))

        self.localSidereTime = LSTDeg
        localHourAngle = LSTDeg - self.equatorial[0]
        self.localHourAngle = (
                localHourAngle if localHourAngle > 0 else localHourAngle + 360.)

    def getAltAziFromRADEC(self, beamCoordinates, LSTDeg, arrayRefereceLatitude):
        beamCoordinates = np.array(beamCoordinates)
        RA = np.deg2rad(beamCoordinates[:,0])
        DEC = np.deg2rad(beamCoordinates[:,1])
        LST = np.deg2rad(LSTDeg)

        altitude, azimuth = convertEquatorialToHorizontal(
                RA, DEC, LST, arrayRefereceLatitude)

        return altitude, azimuth

def readCoordinates(coordinateFileName, delimiter=None):
    coordinateMatrix = None
    with open(coordinateFileName, 'r') as coordFile:
        coordinateArray = np.loadtxt(coordFile, delimiter=delimiter)

    return  coordinateArray

def convertEquatorialToHorizontal(RA, DEC, LST, latitude):
    '''local hour angle'''
    LHA = LST  - RA
    altitude = np.arcsin(np.sin(latitude)*np.sin(DEC) + np.cos(latitude)*np.cos(DEC)*np.cos(LHA))
    azimuth = np.arccos((np.sin(DEC) - np.sin(altitude)*np.sin(latitude))/(np.cos(altitude)*np.cos(latitude)))

    mask = np.sin(LHA) > 0.0
    azimuth[mask] = np.pi*2 - azimuth[mask]

    return altitude, azimuth

def convertHorizontalToEquatorial(azimuth, altitude, LST, latitude):

    DEC = np.arcsin(np.sin(altitude)*np.sin(latitude) + np.cos(altitude)*np.cos(latitude)*np.cos(azimuth))
    # print altitude, DEC, latitude
    cosLHA = (np.sin(altitude)-np.sin(DEC)*np.sin(latitude))/(np.cos(DEC)*np.cos(latitude))
    if cosLHA > 1 and (cosLHA - 1) < 0.0000000000000003:
        print("cos(LHA) corrected form %.20f to 1.0" % cosLHA)
        cosLHA = 1.
    # LHA = np.arccos((np.sin(altitude)-np.sin(DEC)*np.sin(latitude))/(np.cos(DEC)*np.cos(latitude)))
    LHA = np.arccos(cosLHA)

    if np.sin(azimuth) > 0:
        LHA = np.pi*2 - LHA

    RAorig = LST - LHA
    RA = RAorig if RAorig > 0 else RAorig + 2*np.pi

    return RA, DEC

def getHourAngle(RA, LST):

    LHA = LST - RA
    if LHA < 0:
        LHA + 2*np.pi

    return LHA

def projectBaselines(rotatedENU, HA, DEC):

    rotationMatrix = np.array([
            [ np.sin(HA),              np.cos(HA),                  0     ],
            [-np.sin(DEC)*np.cos(HA),  np.sin(DEC)*np.sin(HA), np.cos(DEC)],
            [ np.cos(DEC)*np.cos(HA), -np.cos(DEC)*np.sin(HA), np.sin(DEC)]])

    projectedBaselines = np.dot(rotatedENU, rotationMatrix.T)
    # print rotatedENU
    # print rotationMatrix.T
    # print projectedBaselines

    return projectedBaselines


def rotateENUToEquatorialPlane(ENU, latitude, azimuth, elevation):

    # rotationMatrix = np.array([
            # [0,     -np.sin(latitude), np.cos(latitude)],
            # [1,             0        ,        0        ],
            # [0,      np.cos(latitude), np.sin(latitude)]])

    # rotationMatrix = np.array([
            # [np.cos(latitude),     0, np.sin(latitude) ],
            # [0,             1        ,        0        ],
            # [-np.sin(latitude), 0,      np.cos(latitude)]])

    lengths = np.sqrt(np.sum(np.square(ENU), axis = 1))

    l = latitude

    ENU = np.array(ENU)
    epsilon = 0.000000000001
    azimuths = np.arctan2(ENU[:,0], (ENU[:,1] + epsilon))
    elevations = np.arcsin(ENU[:,2]/(lengths + epsilon))

    a = azimuths
    e = elevations


    rotationMatrix = np.array([
             np.cos(l)*np.sin(e) - np.sin(l)*np.cos(e)*np.cos(a),
                       np.cos(e) * np.sin(a),
             np.sin(l)*np.sin(e) + np.cos(l)*np.cos(e)*np.cos(a)])

    # rotatedENU = lengths * rotationMatrix
    rotatedENU = []

    # print ENU[:,2]
    # print (lengths + epsilon)

    for length, rotates in zip(lengths, rotationMatrix.T):
        rotatedENU.append([length*rotates[0],length*rotates[1],length*rotates[2]])



    # print lengths
    # print rotationMatrix
    # print rotatedENU

    return rotatedENU

def convertENUToUVW(ENU, waveLength, RA, DEC, LST):

    BaseLineInWavelength = ENU/waveLength
    HA = getHourAngle(RA, LST)

    rotationMatrix = np.array([
            [ np.sin(HA),              np.cos(HA),                  0     ],
            [-np.sin(DEC)*np.cos(HA),  np.sin(DEC)*np.sin(HA), np.cos(DEC)],
            [ np.cos(DEC)*np.cos(HA), -np.cos(DEC)*np.sin(HA), np.sin(DEC)]])

    UVW = np.dot(rotationMatrix, BaseLineInWavelength)

    return UVW

def calculateLocalSiderealTime(TimeInUTC, longitude, displayHour=False):

    def calculateDaysSinceJ2000(TimeInUTC):

        '''https://en.wikipedia.org/wiki/Epoch_(astronomy)'''
        J2000UTCTime=datetime.datetime(2000, 1, 1, 11,58,56)
        # compareTime=datetime.datetime(TimeInUTC[0], TimeInUTC[1], TimeInUTC[2],
                # TimeInUTC[3], TimeInUTC[4], TimeInUTC[5])
        timeDiff = TimeInUTC - J2000UTCTime
        timeDiffInDays = timeDiff.days + timeDiff.seconds/(60.0*60*24) + timeDiff.microseconds/(1000000.0*60.0*60*24)

        return timeDiffInDays


    if type(TimeInUTC) != datetime.datetime:
        TimeInUTC = epochToDatetime(TimeInUTC)
    daysSinceJ2000 = calculateDaysSinceJ2000(TimeInUTC)
    # print daysSinceJ2000
    '''http://www.stargazing.net/kepler/altaz.html'''
    LocalSiderealTime = 100.46 + 0.985647 * daysSinceJ2000 +longitude +\
                    15.0*(TimeInUTC.hour + TimeInUTC.minute/60.0 +\
                    TimeInUTC.second/3600.0 + TimeInUTC.microsecond/1000000.0/3600.0)

    LocalSiderealTime %= 360.0

    if displayHour == True:
        hour = LocalSiderealTime/360*24.
        minute = (hour - int(hour))*60
        second = (minute - int(minute))*60

        # print(int(hour), ":", int(minute), ":", second)
        print("{}:{}:{}".format(int(hour), int(minute), second))


    return LocalSiderealTime

def datetimeToEpoch(datetimeObjs):
    observeSeconds = []
    epoch = datetime.datetime.utcfromtimestamp(0)
    if type(datetimeObjs) == list:
        for datetimeObj in datetimeObjs:
            observeSeconds.append((datetimeObj - epoch).total_seconds())
    else:
        observeSeconds = (datetimeObjs - epoch).total_seconds()

    return observeSeconds

def epochToDatetime(epoches):
    epochDatetime = datetime.datetime.utcfromtimestamp(0)
    observeDatetime = []
    if type(epoches) == list:
        for epoch in epoches:
            timeDelter = datetime.timedelta(second = epoch)
            observeDatetime.append(epochDatetime + timeDelta)
    else:
        timeDelta = datetime.timedelta(seconds = epoches)
        observeDatetime = epochDatetime + timeDelta

    return observeDatetime


def convertGodeticToECEF(geodetics):

    '''http://itrf.ensg.ign.fr/faq.php?type=answer#question2'''
    grs80 = nv.FrameE(name='GRS80')
    ecefPoints = np.empty((0,3))
    for lat, lon, height in geodetics:
        geoPoint = grs80.GeoPoint(latitude=lat,
                longitude=lon, z=-height, degrees=True)
        ecefPoint = geoPoint.to_ecef_vector().pvector.ravel()
        ecefPoints = np.append(ecefPoints, [ecefPoint], axis = 0)

    return ecefPoints

def convertECEFToENU(ECEF, ECEFReference, GeodeticReference):
    offset = ECEF - ECEFReference
    lon = np.deg2rad(GeodeticReference[1])
    lat = np.deg2rad(GeodeticReference[0])
    rotationMatrix = np.array([
            [-np.sin(lon),              np.cos(lon),                  0     ],
            [-np.sin(lat)*np.cos(lon), -np.sin(lat)*np.sin(lon), np.cos(lat)],
            [ np.cos(lat)*np.cos(lon),  np.cos(lat)*np.sin(lon), np.sin(lat)]])

    return np.dot(offset, rotationMatrix.T)

def distances(vector):
    squares = np.square(vector)
    if len(squares.shape) > 1:
        elementWiseSum = np.sum(squares, axis=1)
    else:
        elementWiseSum = np.sum(squares)
    squareRoots = np.sqrt(elementWiseSum)
    return squareRoots

def angleToHour(angle, strfmt=True):
    hour = angle/360.*24.
    minute = (hour - int(hour))*60
    second = (minute - int(minute))*60

    hour = int(hour)
    minute = int(minute)

    if strfmt == True:
        return ":".join([str(hour), str(minute), str(second)])
    else:
        return hour, minute, second

def angleToDEC(angle, strfmt=True):
    if angle < 0:
        sign = -1
        angle = abs(angle)
    else:
        sign = 1
    # degree = int(angle)
    minute = (angle - int(angle))*60
    second = (minute - int(minute))*60

    degree = sign * int(angle)
    minute = int(minute)

    if strfmt == True:
        return ":".join([str(degree), str(minute), str(second)])
    else:
        return degree, minute, second

def projectedRotate(altitude, azimuth, baseline, angle):
    # rotate vector on a surface
    # https://math.stackexchange.com/questions/1830695/
    # https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
    '''different between spherical coordinate and horizontal coordinate system'''
    theta = np.pi/2. - altitude
    phi = np.pi/2.- azimuth
    # phi = azimuth
    sourcePosition = [np.sin(theta)*np.cos(phi),
                   np.sin(theta)*np.sin(phi),
                   np.cos(theta)]
    # sourcePosition = np.array([np.sin(theta)*np.cos(phi) + np.cos(theta)*np.cos(phi) - np.sin(phi),
               # np.sin(theta)*np.sin(phi) + np.cos(theta)*np.sin(phi),
               # np.cos(theta) - np.sin(phi)])


    # projectedRotated = np.cos(angle)*baseline + np.sin(angle)*np.cross(sourcePosition.T, baseline)
    # print sourcePosition, baseline
    projectedRotated = np.cos(angle)*baseline + (1-np.cos(angle))*np.cross(sourcePosition, baseline)

    return projectedRotated

def convert_pixel_coordinate_to_equatorial(pixel_coordinates, bore_sight):

    """
    https://heasarc.gsfc.nasa.gov/docs/fcg/standard_dict.html
    CRVAL: coordinate system value at reference pixel
    CRPIX: coordinate system reference pixel
    CDELT: coordinate increment along axis
    CTYPE: name of the coordinate axis
    """
    step = 1/10000000000.

    wcs_properties = wcs.WCS(naxis=2)
    wcs_properties.wcs.crpix = [0, 0]
    wcs_properties.wcs.cdelt = [step, step]
    wcs_properties.wcs.crval = bore_sight
    wcs_properties.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    scaled_pixel_coordinats = np.array(pixel_coordinates)/step

    equatorial_coodinates = wcs_properties.wcs_pix2world(scaled_pixel_coordinats, 0)

    return equatorial_coodinates


def convert_pixel_length_to_equatorial(axis1, axis2, angle, boresight):

    angle1Radian = np.deg2rad(angle)
    angle2Radian = np.deg2rad(angle + 90)
    pixel1 = [axis1 * np.cos(angle1Radian), axis1 * np.sin(angle1Radian)]
    pixel2 = [axis2 * np.cos(angle2Radian), axis2 * np.sin(angle2Radian)]
    equatorials = convert_pixel_coordinate_to_equatorial(
        [pixel1, pixel2], boresight)

    equatorial1 = SkyCoord(equatorials[0][0], equatorials[0][1], unit="deg")
    equatorial2 = SkyCoord(equatorials[1][0], equatorials[1][1], unit="deg")
    equatorialB = SkyCoord(boresight[0], boresight[1], unit="deg")

    width1 = equatorial1.separation(equatorialB)
    width2 = equatorial2.separation(equatorialB)

    # width1 = np.sqrt(np.sum(np.square(equatorials[0]-boresight)))
    # width2 = np.sqrt(np.sum(np.square(equatorials[1]-boresight)))


    return width1, width2

def convertBoresightToHour(boresight):
    boresight = SkyCoord(boresight[0], boresight[1], frame='icrs', unit='deg');
    RA = ":".join([str(int(boresight.ra.hms.h)), str(int(boresight.ra.hms.m)), str(boresight.ra.hms.s)])
    DEC = ":".join([str(int(boresight.dec.dms.d)), str(abs(int(boresight.dec.dms.m))), str(abs(boresight.dec.dms.s))])
    return (RA, DEC)

def convertBoresightToDegree(boresight):
    boresight_coord = SkyCoord(boresight[0]+" "+boresight[1], unit=(u.hourangle, u.deg))
    return(boresight_coord.ra.deg, boresight_coord.dec.deg)


def writeFits(wcsHeader, data, fileName):
    w = wcs.WCS(naxis=2)

    # "Airy's zenithal" projection
    w.wcs.crpix = wcsHeader['crpix']
    w.wcs.cdelt = wcsHeader['cdelt']
    w.wcs.crval = wcsHeader['crval']
    w.wcs.ctype = wcsHeader['ctype']
    # w.wcs.set_pv([(2, 1, 45.0)])

    header = w.to_header()

    hdu = fits.PrimaryHDU(header=header, data=data)
    # Save to FITS file
    hdu.writeto(str(fileName), overwrite=True)


