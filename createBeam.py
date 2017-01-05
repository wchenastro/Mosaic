#!/usr/bin/env python
from math import sqrt, pi, ceil
import sys
import random
import numpy as np


def hexagonGrid(beamNumber, beamRadius, subBeamRadius=None):
    sideLength = beamRadius*2
    if subBeamRadius == None:
        subBeamRadius = sqrt(sideLength**2 / beamNumber / pi )

    horizontalNumber = int(ceil(sideLength/(2*subBeamRadius)))


    if horizontalNumber % 2 == 0:
        horizontalNumber += 1

    coordinates = []
    evenLine = True

    oddLine = []
    evenLine = []

    # print 'horizon number', horizontalNumber

    oddLine.append([0,0])
    for i in range((horizontalNumber - 1)/2):
        oddLine.append([+(i+1)*2*subBeamRadius , 0])
        oddLine.append([-(i+1)*2*subBeamRadius , 0])
    for i in range((horizontalNumber - 1)/2):
        evenLine.append([+subBeamRadius + i*2*subBeamRadius , 0])
        evenLine.append([-subBeamRadius - i*2*subBeamRadius , 0])

    coordinates += oddLine
    twoLines = [evenLine, oddLine]
    lineType = 0
    verticalOffset = subBeamRadius * sqrt(3) * 0.5
    for i in range((horizontalNumber - 1)/2):
        coordinates += [[x, y+(i+1)*2*verticalOffset] for x, y in twoLines[lineType]]
        coordinates += [[x, y-(i+1)*2*verticalOffset] for x, y in twoLines[lineType]]
        lineType = abs(lineType - 1)


    # with open('coord', 'w') as coordFile:
        # for x, y in coordinates:
            # coordFile.write(' '.join([str(x), str(y)]) + '\n')


    inCircleCounter = 0;
    inCircleCoordinates = []
    distanceSqure = beamRadius*beamRadius
    for x, y in coordinates:
        if (x**2 + y**2) <= distanceSqure:
            inCircleCounter += 1
            inCircleCoordinates.append([x,y])

    # print inCircleCounter

    return inCircleCoordinates, subBeamRadius

def randomGrid(beamNumber, beamRadius, subBeamRadius=None):

    axisMin = 0.
    axisMax = beamRadius*2

    subBeamRadius = beamRadius**2 / beamNumber / 2.

    coordinates =  np.random.uniform(axisMin, axisMax, (beamNumber, 2))


    return coordinates, subBeamRadius

def recGrid(beamNumber, beamRadius, subBeamRadius=None):
    sideLength = 2*beamRadius
    if subBeamRadius == None:
        subBeamRadius = sqrt(sideLength*sideLength/beamNumber)/2.

    gridDivider = int(sideLength/2/subBeamRadius)
    if gridDivider % 2 == 0:
        gridDivider +=1

    coordinates = []

    singleLine = [[0,0]]
    for i in range((gridDivider - 1)/2):
        singleLine.append([+(i+1)*2*subBeamRadius , 0])
        singleLine.append([-(i+1)*2*subBeamRadius , 0])

    coordinates += singleLine
    for i in range((gridDivider - 1)/2):
        coordinates += [[x, y+(i+1)*2*subBeamRadius] for x, y in singleLine]
        coordinates += [[x, y-(i+1)*2*subBeamRadius] for x, y in singleLine]


    return coordinates, subBeamRadius



def squareGrid(beamNumber, beamRadius, subBeamRadius=None):
    sideLength = 2*beamRadius
    if subBeamRadius == None:
        subBeamRadius = sqrt(sideLength*sideLength/beamNumber)/2.

    gridDivider = int(sideLength/2/subBeamRadius)
    if gridDivider % 2 == 0:
        gridDivider +=1

    coordinates = []

    singleLine = [[0,0]]
    for i in range((gridDivider - 1)/2):
        singleLine.append([+(i+1)*2*subBeamRadius , 0])
        singleLine.append([-(i+1)*2*subBeamRadius , 0])

    coordinates += singleLine
    for i in range((gridDivider - 1)/2):
        coordinates += [[x, y+(i+1)*2*subBeamRadius] for x, y in singleLine]
        coordinates += [[x, y-(i+1)*2*subBeamRadius] for x, y in singleLine]

    inCircleCounter = 0;
    inCircleCoordinates = []
    for x, y in coordinates:
        if (x**2 + y**2) <= beamRadius*beamRadius:
            inCircleCounter += 1
            inCircleCoordinates.append([x,y])

    # print inCircleCounter

    # with open('coord', 'w') as coordFile:
        # for x, y in coordinates:
            # coordFile.write(' '.join([str(x), str(y)]) + '\n')

    return inCircleCoordinates, subBeamRadius

def optimizeGrid(beamNumber, beamRadius, beamPattern, error, boreSight = None):

    inCircleCoordinates, subBeamRadius = beamPattern(beamNumber, beamRadius)

    inCircleCount = len(inCircleCoordinates)
    trialCount = 0
    while(abs( inCircleCount - beamNumber) > error):

        if inCircleCount <= beamNumber:
            factor = random.uniform(0.8, 1.0)
        else:
            factor = random.uniform(1.0, 1.2)

        inCircleCoordinates, subBeamRadius = beamPattern(beamNumber, beamRadius, subBeamRadius*factor)
        inCircleCount = len(inCircleCoordinates)
        trialCount += 1
        if trialCount > 150:
            print('maximum trials')
            break
        # print(inCircleCount, subBeamRadius)

    if boreSight != None:
        inCircleCoordinates = [[x + boreSight[0], y + boreSight[1]] for x, y in inCircleCoordinates]
        # inCircleCoordinates.insert(0, boreSight)

    with open('inCoord', 'w') as inCoordFile:
        for x, y in inCircleCoordinates:
            inCoordFile.write(' '.join([str(x), str(y)]) + '\n')

    return inCircleCoordinates, subBeamRadius



def main():
    parameterLength = len(sys.argv)
    if parameterLength < 1+2:
        print 'not enough parameters'
        exit()

    beamNumber = int(sys.argv[1])
    beamRadius = float(sys.argv[2])

    beamBoreSight = None
    if len(sys.argv) == 1+4:
        beamBoreSight = [float(sys.argv[3]), float(sys.argv[4])]
    elif parameterLength == 1+3:
        print 'boreSight coordinates is not valid'
        exit()

    coordinates, subBeamRadius = optimizeGrid(beamNumber, beamRadius, hexagonGrid, beamBoreSight)

if __name__ == '__main__':
    main()
