#!/usr/bin/env python
from math import sqrt, pi, ceil
import random, logging
import numpy as np
from functools import partial

logger = logging.getLogger(__name__)

def createTiling(method, beamNumber, beamshape, overlap, tilingShape, parameter, error, seed=None):

    if method == "variable_size":
        axisH, axisV, angle = beamshape.width_at_overlap(overlap)
        area = beamNumber*np.pi*axisH*axisV
        if tilingShape == "circle":
            scale = np.sqrt(area/np.pi)*1.05
            gridding = partial(isInsideCircle, radius=scale)
        elif tilingShape == "hexagon":
            scale = np.sqrt(area/np.pi)*1.15
            gridding = partial(isInsideHexagon, circumradius=scale)

    elif method == "variable_overlap":
        if tilingShape == "circle":
            scale = parameter
            gridding = partial(isInsideCircle, radius=scale)

        elif tilingShape  == "ellipse":
            a, b , orientation = parameter
            # scale = a*1.2 if a > b else b*1.2
            scale = a if a > b else b
            gridding = partial(isInsideEllipse,
                    majorAxis=a, minorAxis=b, orientation=np.deg2rad(orientation))

        elif tilingShape  == "hexagon":
            scale, angle = parameter
            gridding = partial(isInsideHexagon,
                    circumradius=scale, angle = np.deg2rad(angle))

        elif tilingShape == "polygon":
            vertices = parameter
            scale = np.max(np.sqrt(np.square(vertices).sum(axis=1)))
            gridding = partial(isInsidePolygon, vertices=vertices)

        elif tilingShape == "annulus":
            annuluses = []
            scale = None
            for annulus in parameter:
                if annulus[0] == "polygon":
                    vertices = annulus[1]
                    annulusPart = partial(isInsidePolygon, vertices=vertices)
                    if scale == None:
                        scale = np.max(np.sqrt(np.square(vertices).sum(axis=1)))
                elif annulus[0] == "ellipse":
                    a, b , orientation = annulus[1]
                    annulusPart = partial(isInsideEllipse, majorAxis=a,
                        minorAxis=b, orientation=np.deg2rad(orientation))
                    if scale == None:
                        scale = a if a > b else b

                annuluses.append(annulusPart)

            gridding = partial(isInsideAnnulus, annuluses=annuluses)

        beamshapeModel = beamshape.beamshapeModel
        areaMatrix = np.pi * (beamshapeModel[:, 1] * beamshapeModel[:, 2])
        betterBeamArea = np.pi*scale*scale/beamNumber
        beamshapeIndex = np.argmin(np.square(areaMatrix - betterBeamArea))
        maxBeamshapeIndex = len(beamshapeModel) - 1
        if beamshapeIndex <= 0 or beamshapeIndex >= maxBeamshapeIndex:
            logger.warning("variable overlap outside (0, 1).")
            condition = {"optimized":False, "trial_count":-1, "cause": "off limit"}
            return [[0, 0]], scale, (0, 0, 0, 0), condition
        overlap, axisH, axisV, angle = beamshapeModel[beamshapeIndex]

    error = int(np.round(error/2.))
    if seed is None:
        np.random.seed(0)
    else:
        np.random.seed(int(seed))

    trialCount = 0
    betterTiling = [np.Inf]
    maxInsideCount  = -1
    maxTrial = 150
    indexList = []
    search, refine = list(range(2))
    state = search
    beamshapeIndexDelta = 0
    maxSearchCount = 3
    isReachValleyFast = False
    # smallError = max(int(error/3.), 1)
    smallError = max(int(error/2.), 1)
    diffRatio = 0
    condition = {"optimized":False, "trial_count":0}
    while(trialCount < maxTrial):

        insideCoordinates = createGrid(scale, axisH, axisV, angle, gridding)
        insideCount = len(insideCoordinates)
        beamNumDiff = abs(insideCount - (beamNumber))
        # if method == "variable_overlap":
            # print("{:4d}".format(insideCount), "{:4d}".format(beamNumDiff), beamshapeIndex,
                    # "{:.4f}".format(overlap), state, trialCount, beamshapeIndexDelta)
        # else:
            # print("{:3d}".format(insideCount), "{:2d}".format(beamNumDiff),
                    # "{:.7f}".format(beamshapeIndexDelta), "{:.4f}".format(scale), state, trialCount)

        # if beamNumber >= insideCount and (beamNumber - insideCount) <= smallError:
        if beamNumber == insideCount:
            condition["optimized"] = True;
            break
        # elif beamNumber < insideCount and (insideCount - beamNumber) <= int(smallError/3.):
            # break
        elif state == search and (beamNumDiff <= error or trialCount >= maxSearchCount):
            state = refine

        if insideCount < beamNumber and beamNumDiff < betterTiling[0]:
        # if beamNumDiff < smallError and beamNumDiff < betterTiling[0]:
            betterTiling = [beamNumDiff, insideCoordinates, scale, (axisH, axisV, angle, overlap)]


        if method == "variable_overlap":
            if state == search:
                beamArea = np.pi * (axisH * axisV)
                if insideCount > maxInsideCount:
                    maxInsideCount = insideCount
                    tilingArea = beamArea * max(insideCount, 1)
                betterBeamArea = np.random.uniform(0.95, 1.05) * tilingArea / beamNumber
                beamshapeIndex = np.argmin(np.square(areaMatrix - betterBeamArea))
                overlap, axisH, axisV, angle = beamshapeModel[beamshapeIndex]
            if state == refine:
                if trialCount < maxSearchCount:
                    isReachValleyFast = True
                    beamshapeIndexDelta = 40
                elif trialCount > maxSearchCount:
                    if insideCount < beamNumber and indexList[-1] >= indexList[-2] and indexList[-2] >= indexList[-3]:
                        beamshapeIndexDelta *= 2
                    elif insideCount > beamNumber and indexList[-1] <= indexList[-2] and indexList[-2] <= indexList[-3]:
                        beamshapeIndexDelta *= 2
                    else:
                        # if trialCount % 30 == 0 and beamshapeIndexDelta == 1:
                            # logger.warning("the optimization runs into a loop, "
                                # "but the number of generated beams is within the rough threshold.")
                        if beamshapeIndexDelta == 1:
                            condition["optimized"] = True;
                            break
                        beamshapeIndexDelta = max(int(beamshapeIndexDelta/2.), 1)
                elif trialCount == maxSearchCount and isReachValleyFast == False:
                    beamshapeIndexDelta = max(int(np.abs(np.round((indexList[-1] - indexList[-2])/2.))), 1)
                if beamshapeIndexDelta >= 10:
                  beamshapeIndexDelta = int(np.round(beamshapeIndexDelta * np.random.uniform(0.9, 1.1)))

                if insideCount < beamNumber:
                    beamshapeIndex += beamshapeIndexDelta
                else:
                    beamshapeIndex -= beamshapeIndexDelta
                if beamshapeIndex <= 0 or beamshapeIndex >= maxBeamshapeIndex:
                    logger.warning("variable overlap outside (0, 1).")
                    condition = {"optimized":False, "trial_count":trialCount, "cause": "off limit"}
                    break
                overlap, axisH, axisV, angle = beamshapeModel[beamshapeIndex]
            indexList.append(beamshapeIndex)

        elif method == "variable_size":
            if state == search:
                diffRatio = 5.0*beamNumDiff/beamNumber
                if insideCount <= beamNumber:
                    factor = random.uniform(1, 1+0.1*diffRatio)
                else:
                    factor = random.uniform(1-0.1*diffRatio, 1)
                scale *= factor
            if state == refine:
                if trialCount < maxSearchCount:
                    isReachValleyFast = True
                    diffRatio = 5.0*beamNumDiff/beamNumber
                    beamshapeIndexDelta = scale * 0.05*diffRatio
                elif trialCount > maxSearchCount:
                    if insideCount < beamNumber and indexList[-1] >= indexList[-2] and indexList[-2] >= indexList[-3]:
                        beamshapeIndexDelta *= 2
                    elif insideCount > beamNumber and indexList[-1] <= indexList[-2] and indexList[-2] <= indexList[-3]:
                        beamshapeIndexDelta *= 2
                    else:
                        if beamshapeIndexDelta == minBeamshapeIndexDelta:
                            condition["optimized"] = True;
                            break
                        beamshapeIndexDelta = max(beamshapeIndexDelta/2., minBeamshapeIndexDelta)
                elif trialCount == maxSearchCount:
                    if isReachValleyFast == False:
                        beamshapeIndexDelta = np.abs((indexList[-1] - indexList[-2])/2.)
                    minBeamshapeIndexDelta = beamshapeIndexDelta / 10.
                if insideCount < beamNumber:
                    scale += beamshapeIndexDelta
                else:
                    scale -= beamshapeIndexDelta
            indexList.append(scale)
            if tilingShape == "circle":
                gridding = partial(isInsideCircle, radius=scale)
            elif tilingShape == "hexagon":
                gridding = partial(isInsideHexagon, circumradius=scale)

        """
        if trialCount == 20:
           print("re-interpolate!")
           print(beamshapeModel.shape)
           lowerBound = min(indexList[5:])
           upperBound = max(indexList[5:])
           interval = 0.000000001
           level = beamshapeModel[lowerBound:upperBound+1][:,0]
           levelInterp = np.arange(level[0] , level[-1]+interval, interval)
           axis1Interp = np.interp(levelInterp, level, beamshapeModel[lowerBound:upperBound+1][:,1])
           axis2Interp = np.interp(levelInterp, level, beamshapeModel[lowerBound:upperBound+1][:,2])
           angleInterp = np.interp(levelInterp, level, beamshapeModel[lowerBound:upperBound+1][:,3])
           beamshapeModel = np.array([levelInterp, axis1Interp, axis2Interp, angleInterp]).T
           print(len(level), len(levelInterp), beamshapeModel.shape, lowerBound, upperBound)
           beamshapeIndex = int(len(levelInterp)/2)
        """


        trialCount += 1
        if trialCount > 149:
            condition["optimized"] = False;
            condition["cause"] = "maximum trials";
            actualDiff = np.Inf
            break
            # if abs(insideCount - beamNumber) < int(beamNumber*0.1):
                # logger.warning("maximum trials reached in the tiling process, "
                    # "the tiling is not well optimized, the difference in number "
                    # "of beams between requried and generated is less than 10%")
                # break
            # else:
                # logger.critical("maximum trials reached in the tiling process, "
                    # "the tiling is not optimized. the number of requried beams and "
                    # "generated beams is %d and %d, please consider increasing "
                    # "the margin threshold if needed." % (beamNumber, insideCount))
                # break

    if betterTiling != [np.Inf] and beamNumber < insideCount or (beamNumber - insideCount) > betterTiling[0]:
    # if  beamNumDiff > betterTiling[0]:
        # print("use better tiling")
        insideCoordinates, scale, (axisH, axisV, angle, overlap) = betterTiling[1:]
    # logger.info("tiling: required_beam_number: {}, generate_beam_number: {}, "
                # "trial counter: {}".format(beamNumber, len(insideCoordinates), trialCount))

    condition["trial_count"] = trialCount
    return insideCoordinates, scale, (axisH, axisV, angle, overlap), condition


def createGrid(scale, axisH, axisV, angle, isInsideBoundary):

    sideLength = scale*2

    horizontalNumber = int(np.ceil(sideLength/(2*axisH)))
    verticalOffset = axisV * sqrt(3) * 0.5
    verticalNumber = int(np.ceil(sideLength/(2*verticalOffset)))

    if horizontalNumber % 2 == 0:
        horizontalNumber += 1
    if verticalNumber % 2 == 0:
        verticalNumber += 1

    coordinates = []
    evenLine = True

    oddLine = []
    evenLine = []

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
    maxAxis = axisH if axisH > axisV else axisV
    if maxAxis > scale * 0.7:
        for i in range(int((verticalNumber - 1)/2)):
            coordinates += [[x, y+(i+1)*2*verticalOffset] for x, y in oddLine]
            coordinates += [[x, y-(i+1)*2*verticalOffset] for x, y in oddLine]
    else:
        for i in range(int((verticalNumber - 1)/2)):
            coordinates += [[x, y+(i+1)*2*verticalOffset] for x, y in twoLines[lineType]]
            coordinates += [[x, y-(i+1)*2*verticalOffset] for x, y in twoLines[lineType]]
            lineType = abs(lineType - 1)

    angle = np.deg2rad(angle)

    if angle != 0:
        rotationMatrix = [[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]]
        coordinatesRotated = np.dot(np.array(rotationMatrix), np.array(coordinates).T).T
    else:
        coordinatesRotated = coordinates

    insideCoordinatesRotated = isInsideBoundary(coordinates = coordinatesRotated)

    return insideCoordinatesRotated


def isInsideHexagon(coordinates, circumradius, angle=0):

    insideCoordinates = []

    if angle == 0:
        for coord in coordinates:
            x, y = abs(coord[0]), abs(coord[1])
            if x < 3**0.5 * min(circumradius - y, circumradius / 2.):
                insideCoordinates.append(coord)
        return insideCoordinates
    else:
        rotationMatrix = [[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]]
        coordinatesRotated = np.dot(np.array(rotationMatrix), np.array(coordinates).T).T

        for coord in coordinatesRotated:
            x, y = abs(coord[0]), abs(coord[1])
            if x < 3**0.5 * min(circumradius - y, circumradius / 2.):
                insideCoordinates.append(coord)
        rotationMatrix = [[np.cos(-angle), -np.sin(-angle)], [np.sin(-angle), np.cos(-angle)]]
        insideCoordinatesRotatedBack = np.dot(np.array(rotationMatrix), np.array(insideCoordinates).T).T
        return insideCoordinatesRotatedBack

def isInsidePolygon(coordinates, vertices):

    insideCoordinates = []

    verticesLength = len(vertices)
    for x,y in coordinates:
        c = False
        j = verticesLength - 1
        for i in range(verticesLength):
            if ((vertices[i][1] > y) is not (vertices[j][1] > y)) and (
                    x < (vertices[j][0] - vertices[i][0]) * (y - vertices[i][1]) / (vertices[j][1] - vertices[i][1]) + vertices[i][0]):
                c = not c;
            j = i
            i += 1

        if c is True:
            insideCoordinates.append([x,y])

    return insideCoordinates

def isInsideCircle(coordinates, radius):

    radiusSquare = radius**2
    distance = np.sum(np.square(coordinates), axis=1)

    return coordinates[distance < radiusSquare]


def isInsideEllipse(coordinates, majorAxis, minorAxis, orientation):

    cosa = np.cos(orientation)
    sina = np.sin(orientation)

    coordinates = np.array(coordinates)
    x = coordinates[:, 0]
    y = coordinates[:, 1]

    distance = np.power((cosa*x + sina*y)/majorAxis,2) +\
             np.power((sina*x - cosa*y)/minorAxis,2)

    return coordinates[distance < 1.0]

def isInsideAnnulus(coordinates, annuluses):

    outer = annuluses[0]
    inner = annuluses[1]

    outerCoord = outer(coordinates = coordinates)
    innerCoord = inner(coordinates = coordinates)

    inSideCoordinates = set([tuple(o) for o in outerCoord]) - set([tuple(i) for i in innerCoord])

    inSideCoordinates = [[t[0], t[1]] for t in list(inSideCoordinates)]

    return inSideCoordinates
