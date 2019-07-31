#!/usr/bin/env python
from math import sqrt, pi, ceil
import sys
import random, logging
import numpy as np

logger = logging.getLogger(__name__)

def hexagonGrid(beamNumber, beamRadius, subBeamRadius=None):
    sideLength = beamRadius*2
    if subBeamRadius is None:
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

# def GaussianQuntileFunc(mean, sigma, probability, inverfFunc=None):

    # if inverfFunc == None:
        # try:
            # from scipy.special import erfinv
        # except ImportError:
            # print("no erfinv function available")
            # raise
        # inverfFunc = erfinv

    # return mean + sigma*np.sqrt(2)*inverfFunc(2*probability - 1)

class TilingNotOptimizedError(Exception):
    pass

def ellipseCompact(beamNumber, axisH, axisV, angle, error, seed=None, write=False):


    #angle = 180-angle
    area = beamNumber*np.pi*axisH*axisV
    beamRadius = np.sqrt(area/np.pi)*1.
    error = int(round(error/2.))
    if seed is None:
        random.seed(axisH)
    else:
        random.seed(seed)

    inCircleCoordinates = ellipseGrid(beamRadius, axisH, axisV, angle, write=False)

    inCircleCount = inCircleCoordinates.shape[1]
    trialCount = 0
    "the total number of beams should not less than required, so a padding error"
    "is added to the actual required beam number"
    while(abs(inCircleCount - (beamNumber+error)) > error):

        if inCircleCount <= beamNumber:
            factor = random.uniform(1.0, 1.1)
        else:
            factor = random.uniform(0.9, 1.0)
        beamRadius = beamRadius*factor
        inCircleCoordinates = ellipseGrid(beamRadius, axisH, axisV, angle)
        inCircleCount = inCircleCoordinates.shape[1]
        trialCount += 1
        if trialCount > 150:
            if abs(inCircleCount - beamNumber) < int(beamNumber*0.1):
                logger.warning("maximum trials reached in the tiling process, "
                    "the tiling is not well optimized, the difference in number "
                    "of beams between requried and generated is less than 10%")
                break
            else:
                logger.critical("maximum trials reached in the tiling process, "
                    "the tiling is not optimized. the number of requried beams and "
                    "generated beams is %d and %d, please consider increasing "
                    "the margin threshold if needed." % (beamNumber, inCircleCount))
                break
                # raise TilingNotOptimizedError({
                        # "message":"maximum trials reached in the tiling process",
                        # "required_beam_number":beamNumber,
                        # "generated_beam_number":inCircleCount,
                        # })
    logger.info("tiling result, required_beam_number: {}, generate_beam_number: {}, "
                "trial counter: {}".format(beamNumber, inCircleCount, trialCount))

    random.seed()
    # print inCircleCount
    if(write==True):
        with open('ellipsePack', 'w') as inCoordFile:
            for x, y in inCircleCoordinates.T:
                inCoordFile.write(' '.join([str(x), str(y)]) + '\n')

    return inCircleCoordinates.T, beamRadius



def ellipseGrid(beamRadius, axisH, axisV, angle, write=False):
    sideLength = beamRadius*2

    horizontalNumber = int(ceil(sideLength/(2*axisH)))
    verticalOffset = axisV * sqrt(3) * 0.5
    verticalNumber = int(ceil(sideLength/(2*verticalOffset)))



    if horizontalNumber % 2 == 0:
        horizontalNumber += 1
    if verticalNumber % 2 == 0:
        verticalNumber += 1

    coordinates = []
    evenLine = True

    oddLine = []
    evenLine = []

    # print 'horizon number', horizontalNumber

    oddLine.append([0,0])
    for i in range(int((horizontalNumber - 1)/2)):
        oddLine.append([+(i+1)*2*axisH , 0])
        oddLine.append([-(i+1)*2*axisH , 0])
    for i in range(int((horizontalNumber - 1)/2)):
        evenLine.append([+axisH + i*2*axisH , 0])
        evenLine.append([-axisH - i*2*axisH , 0])

    coordinates += oddLine
    twoLines = [evenLine, oddLine]
    lineType = 0
    # verticalOffset = subBeamRadius * sqrt(3) * 0.5
    maxAxis = axisH if axisH > axisV else axisV
    if maxAxis > beamRadius * 0.7:
        for i in range(int((verticalNumber - 1)/2)):
            coordinates += [[x, y+(i+1)*2*verticalOffset] for x, y in oddLine]
            coordinates += [[x, y-(i+1)*2*verticalOffset] for x, y in oddLine]
    else:
        for i in range(int((verticalNumber - 1)/2)):
            coordinates += [[x, y+(i+1)*2*verticalOffset] for x, y in twoLines[lineType]]
            coordinates += [[x, y-(i+1)*2*verticalOffset] for x, y in twoLines[lineType]]
            lineType = abs(lineType - 1)



    # with open('coord', 'w') as coordFile:
        # for x, y in coordinates:
            # coordFile.write(' '.join([str(x), str(y)]) + '\n')

    # from sympy import Ellipse, Circle, Point

    # from scipy.optimize import fsolve


    inCircleCoordinates = []
    counter = 0
    # primaryBeam = Circle(Point(0,0), beamRadius)
    # def makeFuns(x0, y0):
        # def funs(p):
            # x,y = p
            # return (x**2 + y**2 - beamRadius**2, (x-x0)**2/axisH**2 + (y-y0)**2/axisV**2 - 1)
        # return funs
    for x0, y0 in coordinates:
        if (x0**2 + y0**2) >= beamRadius*beamRadius: continue
        inCircleCoordinates.append([x0,y0])
        # result = fsolve(makeFuns(x0, y0), (1,1), full_output=True)
        # print(result[2], x0,y0)
        # if result[2] > 2:
            # print(result[3])
            # inCircleCoordinates.append([x0,y0])
        # ellipse = Ellipse(Point(x,y), axisH, axisV)
        # if len(Ellipse(Point(x,y), axisH, axisV).intersection(primaryBeam)) == 0:
        # print(x,y)
        # print(counter)

    # print len(coordinates), len(inCircleCoordinates)

    angle = np.deg2rad(angle)

    if angle == 0:
        inCircleCoordinatesRotated = np.array(inCircleCoordinates).T
    else:
        rotationMatrix = [[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]]
        inCircleCoordinatesRotated = np.dot(np.array(rotationMatrix), np.array(inCircleCoordinates).T)

    if write == True:
        with open('ellipsePack', 'w') as coordFile:
            for x, y in inCircleCoordinatesRotated.T:
                coordFile.write(' '.join([str(x), str(y)]) + '\n')



    # print inCircleCoordinatesRotated.shape
    return inCircleCoordinatesRotated


def randomGrid(beamNumber, beamRadius, subBeamRadius=None):

    axisMin = 0.
    axisMax = beamRadius*2

    subBeamRadius = beamRadius**2 / beamNumber / 2.

    coordinates =  np.random.uniform(axisMin, axisMax, (beamNumber, 2))


    return coordinates, subBeamRadius

def recGrid(beamNumber, subBeamRadius):

    area = beamNumber*(subBeamRadius**2)
    beamRadius = np.sqrt(area)

    sideLength = 2*beamRadius

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


    return np.array(coordinates), beamRadius



def squareGrid(beamNumber, beamRadius, subBeamRadius=None):
    sideLength = 2*beamRadius
    if subBeamRadius is None:
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

def optimizeGrid(beamNumber, beamRadius, beamPattern, error, boreSight = [0,0]):

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

    if (boreSight != [0,0]).all():
        inCircleCoordinates = [[x + boreSight[0], y + boreSight[1]] for x, y in inCircleCoordinates]
        # inCircleCoordinates.insert(0, boreSight)

    # with open('inCoord', 'w') as inCoordFile:
        # for x, y in inCircleCoordinates:
            # inCoordFile.write(' '.join([str(x), str(y)]) + '\n')

    return inCircleCoordinates, subBeamRadius



def main():
    parameterLength = len(sys.argv)
    if parameterLength < 1+2:
        print('not enough parameters')
        exit()

    beamNumber = int(sys.argv[1])
    beamRadius = float(sys.argv[2])

    beamBoreSight = None
    if len(sys.argv) == 1+4:
        beamBoreSight = [float(sys.argv[3]), float(sys.argv[4])]
    elif parameterLength == 1+3:
        print('boreSight coordinates is not valid')
        exit()

    coordinates, subBeamRadius = optimizeGrid(beamNumber, beamRadius, hexagonGrid, 5, beamBoreSight)
    print(np.rad2deg(subBeamRadius))

if __name__ == '__main__':
    main()
