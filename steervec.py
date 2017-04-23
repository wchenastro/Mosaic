#!/usr/bin/env python

import numpy as np

def waveNumber(altitude, azimuth, waveLength):

    '''different between spherical coordinate and horizontal coordinate system'''
    theta = np.pi/2. - altitude
    phi = np.pi/2.- azimuth
    # phi = azimuth

    u = np.array([np.sin(theta)*np.cos(phi),
               np.sin(theta)*np.sin(phi),
               np.cos(theta)])

    '''the nagetive sign indicates the direction(conjugate weight)'''
    waveNumbers = (-1) * u * 2 * np.pi / waveLength

    return waveNumbers

def projectedBaselines(altitude, azimuth, baselines):
    '''different between spherical coordinate and horizontal coordinate system'''
    theta = np.pi/2. - altitude
    phi = np.pi/2.- azimuth
    # phi = azimuth

    sourcePosition = np.array([np.sin(theta)*np.cos(phi),
               np.sin(theta)*np.sin(phi),
               np.cos(theta)])
    # sourcePosition = np.array([np.sin(theta)*np.cos(phi) + np.cos(theta)*np.cos(phi) - np.sin(phi),
               # np.sin(theta)*np.sin(phi) + np.cos(theta)*np.sin(phi),
               # np.cos(theta) - np.sin(phi)])


    projectedBaselines = []
    for source in sourcePosition.T:
        for baseline in baselines:
            projectedBaselines.append(np.absolute(np.cross(baseline, source)))
    # projectedBaselines  = np.absolute(np.cross(baselines, sourcePosition.T))

    return projectedBaselines


def weightVector(waveNumbers, receiverLocations):

    weights  = np.exp(-1j*np.dot(receiverLocations, waveNumbers))

    return  np.array(weights)

def waveNumberFreq(altitude, azimuth, frequencies):
    speedOfLight = 299792458

    '''different between spherical coordinate and horizontal coordinate system'''
    theta = np.pi/2. - altitude
    phi = np.pi/2.- azimuth
    # phi = azimuth

    u = np.array([np.sin(theta)*np.cos(phi),
               np.sin(theta)*np.sin(phi),
               np.cos(theta)])

    waveNumbers = []
    for frequency in frequencies:
        '''the nagetive sign indicates the direction(conjugate weight)'''
        waveNumbers.append((-1) * u * 2 * np.pi * frequency / speedOfLight)

    return np.array(waveNumbers)



def offsetWeight(waveNumber, receiverLocation):

    '''assume the first coordinate is the boreSight'''
    boreSightWaveNumber = waveNumber[:,[0]]
    boreSightWeight = np.exp(-1j*np.dot(receiverLocation, boreSightWaveNumber))
    '''delete the boreSight coordinate from coordinates array'''
    # waveNumberWithoutBoreSight = np.delete(waveNumber, 0, 0)
    waveNumberWithoutBoreSight = waveNumber[:,1:]

    '''delta of wave number between bore sight and other coordinates'''
    deltaWaveNumber = waveNumberWithoutBoreSight - boreSightWaveNumber

    delterWeights = np.exp(-1j*np.dot(receiverLocation, deltaWaveNumber))

    offsetWeights = boreSightWeight * delterWeights

    weightOfSite = np.concatenate(([boreSightWeight] + [offsetWeights]), axis = 1)

    return weightOfSite


def weightVector2(W):
    '''weight vector from uvw coordinates'''

    return np.exp(-1j*W)



if __name__ == "__main__":
    pass


