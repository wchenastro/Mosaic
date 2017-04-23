#!/usr/bin/env python

import datetime
import numpy as np

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

    # azimuthFix = []
    # for lha, azi in zip(LHA, azimuth):
        # if np.sin(lha) > 0.0:
            # azi = np.pi*2 - azi
        # azimuthFix.append(azi)
    # azimuth = np.array(azimuthFix)
    mask = np.sin(LHA) > 0.0
    azimuth[mask] = np.pi*2 - azimuth[mask]

    return altitude, azimuth

def getHourAngle(RA, LST):

    LHA = LST - RA
    if LHA < 0:
        LHA + 2*np.pi

    return LHA


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
        timeDiffInDays = timeDiff.days + timeDiff.seconds/(60.0*60*24)

        return timeDiffInDays


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

        print int(hour), ":", int(minute), ":", second


    return LocalSiderealTime


def convertGodeticToECEF(geodetics):

    import nvector as nv
    '''http://itrf.ensg.ign.fr/faq.php?type=answer#question2'''
    grs80 = nv.FrameE(name='GRS80')
    ecefPoints = np.empty((0,3))
    for lon, lat, height in geodetics:
        geoPoint = grs80.GeoPoint(latitude=lat,
                longitude=lon, z=-height, degrees=True)
        ecefPoint = geoPoint.to_ecef_vector().pvector.ravel()
        ecefPoints = np.append(ecefPoints, [ecefPoint], axis = 0)

    return ecefPoints


def convertECEFToENU(ECEF, ECEFReference, GeodeticReference):
    offset = ECEF - ECEFReference
    lon = np.deg2rad(GeodeticReference[0])
    lat = np.deg2rad(GeodeticReference[1])
    rotationMatrix = np.array([
            [-np.sin(lon),              np.cos(lon),                  0     ],
            [-np.sin(lat)*np.cos(lon), -np.sin(lat)*np.sin(lon), np.cos(lat)],
            [ np.cos(lat)*np.cos(lon),  np.cos(lat)*np.sin(lon), np.sin(lat)]])

    return np.dot(offset, rotationMatrix.T)

