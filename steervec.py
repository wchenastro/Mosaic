#!/usr/bin/env python

import numpy as np

def waveNumber(altitude, azimuth, waveLength):

    '''different between spherical coordinate and horizontal coordinate system'''
    altitude = 90 - altitude

    u = np.array([np.sin(altitude)*np.cos(azimuth),
               np.sin(altitude)*np.sin(azimuth),
               np.cos(altitude)])

    '''the nagetive sign indicates the direction(conjugate weight)'''
    waveNumbers = (-1) * u * 2 * np.pi / waveLength

    return waveNumbers


def weightVector(waveNumbers, receiverLocations):

    weights  = np.exp(-1j*np.dot(receiverLocations, waveNumbers))

    return  np.array(weights)



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


