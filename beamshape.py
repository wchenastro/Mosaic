import numpy as np
from scipy import interpolate
from plot import plotOverlap


def calculateBeamOverlaps(ellipseCenters, radius, majorAxis, minorAxis, rotation, mode, fileName = None):

    def RotatedGaussian2DPDF(x, y, xMean, yMean, xSigma, ySigma, angle):
        angle = -(angle - np.pi)
        a = np.power(np.cos(angle), 2)/(2*xSigma**2) + np.power(np.sin(angle), 2)/(2*ySigma**2)
        b = - np.sin(2*angle)/(4*xSigma**2) + np.sin(2*angle)/(4*ySigma**2)
        c = np.power(np.sin(angle), 2)/(2*xSigma**2) + np.power(np.cos(angle), 2)/(2*ySigma**2)

        return np.exp(-(a*np.power(x-xMean, 2) + 2*b*(x-xMean)*(y-yMean) + c*np.power(y-yMean, 2)))

    def isInsideEllips(center, majorAxis, minorAxis, rotation, testPointX, testPointY):

        xOffset = testPointX - center[0]
        yOffset = testPointY - center[1]
        cosa = np.cos(rotation)
        sina = np.sin(rotation)
        # majorSquare = np.power(majorAxis,2)
        # minorSquare = np.power(minorAxis,2)

        result = np.power((cosa*xOffset + sina*yOffset)/majorAxis,2) +\
                 np.power((sina*xOffset - cosa*yOffset)/minorAxis,2)


        # result = 1/majorAxis**2 * ((testPointX - center[0])*np.cos(rotation) +\
                # (testPointY - center[1])*np.sin(rotation))**2 +\
                # 1/minorAxis**2 * ((testPointX - center[0])*np.sin(rotation) -\
                # (testPointY - center[1])*np.cos(rotation))**2
        return result

    if mode == "counter":
        mode = 1
    elif mode == "heater":
        mode = 2
    elif mode == "both":
        mode = 3

    # rotated = 0/(3600*24.) * 2 * np.pi
    # rotation += rotated
    longAxis = majorAxis if majorAxis > minorAxis else minorAxis
    gridNum = 1000
    # print 0.3*radius, longAxis
    halfSidelength = 0.15 * radius
    offsetCenter = [radius, radius]
    # np.savetxt('tmp/oldcenter', ellipseCenters)
    ellipseCenters = ellipseCenters + offsetCenter
    # np.savetxt('tmp/newcenter', ellipseCenters)
    innerEllipses = []
    for ellipseCenter in ellipseCenters:
        if (ellipseCenter[0] > (offsetCenter[0]-halfSidelength) and\
            ellipseCenter[0] < (offsetCenter[0]+halfSidelength)) and\
           (ellipseCenter[1] > (offsetCenter[1]-halfSidelength) and\
            ellipseCenter[1] < (offsetCenter[1]+halfSidelength)):
            innerEllipses.append(ellipseCenter)
    paddingRatio = 2*longAxis/halfSidelength
    halfSidelength *= 1 + paddingRatio
    # np.savetxt('tmp/innercenter', innerEllipses)
    # step = 2*halfSidelength/gridNum
    width = longAxis*2
    squareEdgeX = [offsetCenter[0] - halfSidelength, offsetCenter[0] + halfSidelength]
    squareEdgeY = [offsetCenter[1] - halfSidelength, offsetCenter[1] + halfSidelength]
    # print squareEdgeX, squareEdgeY

    # grids = np.mgrid[squareEdgeY[0]:squareEdgeY[1]:step, squareEdgeX[0]:squareEdgeX[1]:step]

    grids = np.meshgrid(np.linspace(squareEdgeX[0], squareEdgeX[1], gridNum),
            np.linspace(squareEdgeX[1], squareEdgeX[0], gridNum))


    # nopoint = []
    # gridLength = grids.shape[1]
    gridLength = gridNum
    overlapCounter = np.zeros((gridLength, gridLength))
    overlapHeater = np.zeros((gridLength, gridLength))
    # overlapCounter = np.full((gridLength, gridLength), np.inf)
    sigmaH, sigmaV = majorAxis * (2./2.3556),  minorAxis  * (2./2.3556)
    for ellipseCenter in innerEllipses:
        horizontalBoarder = [ellipseCenter[0]-width, ellipseCenter[0]+width]
        verticalBoarder = [ellipseCenter[1]-width, ellipseCenter[1]+width]

        horizontalIndex = np.round([(horizontalBoarder[0]-squareEdgeX[0])/(2.0*halfSidelength)*gridNum,
                (horizontalBoarder[1]-squareEdgeX[0])/(2.0*halfSidelength)*gridNum]).astype(int)
        verticalIndex = gridNum - np.round([(verticalBoarder[0]-squareEdgeY[0])/(2.0*halfSidelength)*gridNum,
                (verticalBoarder[1]-squareEdgeY[0])/(2.0*halfSidelength)*gridNum]).astype(int)
        # print verticalIndex, horizontalIndex
        insideThisBorder = np.s_[verticalIndex[1]: verticalIndex[0], horizontalIndex[0]: horizontalIndex[1]]
        gridX = grids[0][insideThisBorder]
        gridY = grids[1][insideThisBorder]

        #heat
        if mode == 2 or mode == 3:
            probability = RotatedGaussian2DPDF(gridX, gridY,ellipseCenter[0],
                    ellipseCenter[1], sigmaH, sigmaV, rotation)

            probabilityMask = (overlapHeater[insideThisBorder] < probability)
            overlapHeater[insideThisBorder][probabilityMask] = probability[probabilityMask]

        #counter
        if mode == 1 or mode == 3:
            counts = isInsideEllips(ellipseCenter, majorAxis, minorAxis, rotation, gridX, gridY)
            countMask = counts<1
            counts[countMask] = 1
            counts[~countMask] = 0
            overlapCounter[insideThisBorder] += counts

        # print ellipseCenter, majorAxis, minorAxis, rotation
        # np.save('tmp/grid', [gridX, gridY])
        # np.savetxt('tmp/pointm', result)
        # exit(0)
        # if np.amin(result) > 1.:
            # np.savetxt('tmp/grid', [gridX,gridY])
            # print ellipseCenter
            # exit()
        # print len(gridX)
        # print result[result<1]
        # print len(gridY), np.amin(result), result[result<1]
    # np.savetxt('tmp/nopoint', nopoint)
    trimmedGridLength = int(np.round(gridLength / (1 + 2*paddingRatio)))
    halfPaddingCount =  int(np.round((gridLength - trimmedGridLength) / 2.))
    overlapCounter = overlapCounter[halfPaddingCount:-halfPaddingCount, halfPaddingCount:-halfPaddingCount]
    overlapHeater = overlapHeater[halfPaddingCount:-halfPaddingCount, halfPaddingCount:-halfPaddingCount]
    # print np.count_nonzero(overlapCounter > 1), np.count_nonzero(overlapCounter == 1), np.count_nonzero(overlapCounter == 0)

    # unique, counts = np.unique(overlapCounter, return_counts=True)
    # print dict(zip(unique, counts))

    if fileName != None:
        prefix, suffix= fileName.split('.')
        if mode == 1 or mode == 3:
            plotOverlap(overlapCounter, fileName = prefix+"Counter."+suffix)
        if mode == 2 or mode == 3:
            plotOverlap(overlapHeater, fileName = prefix+"Heater."+suffix)

    if mode == 1:
        return overlapCounter
    elif mode == 2:
        return overlapHeater
    elif mode == 3:
        return overlapCounter, overlapHeater
    # np.save('overlapCounter', overlapCounter)


def trackBorder(image, threshold = 0.3, density = 20, interpolatedLength = 800):

    interpolater = interpolate.interp2d(range(density), range(density), image ,kind='cubic')
    image = interpolater(np.linspace(0, density - 1, interpolatedLength),
            np.linspace(0, density - 1, interpolatedLength))
    # np.savetxt('interpolate', image)
    interpolatedGrid = interpolatedLength*1.0/density


    imageCenter = [interpolatedLength/2, interpolatedLength/2]
    rowIdx = int(imageCenter[0])
    colIdx = int(imageCenter[1])

    trueCenterIndex = [int(rowIdx + interpolatedGrid/2 + 1), int(colIdx + interpolatedGrid/2 + 1)]

    class State:
        move, addRow, findBorder, findBorderReverse, upsideDown, end = range(6)

    border = []
    state = State.move
    rowStep = -1
    colStep = 1
    colBorder = interpolatedLength - 1
    rowBorder = 0
    bottom = False
    overstep = False
    maxOverstepValue = 0
    offset = 0.1
    # closestToCenter = 1
    # closestToCenterIndex = []

    while state != State.end:
        if state == State.move:
            if colIdx != colBorder:
                colIdx += colStep
                if image[rowIdx, colIdx] < threshold:
                    border.append([rowIdx, colIdx])
                    state = State.addRow
                # else:
                    # distToCenter = 1 - image[rowIdx, colIdx]
                    # if distToCenter  < closestToCenter:
                        # closestToCenter = distToCenter
                        # closestToCenterIndex = [rowIdx, colIdx]
            else:
                overstep = True
                if image[rowIdx, colIdx] > maxOverstepValue:
                    maxOverstepValue = image[rowIdx, colIdx]
                state = State.addRow

        if state == State.addRow:
            if rowIdx != rowBorder:
                rowIdx += rowStep
                state = State.findBorder
            else:
                overstep = True
                if image[rowIdx, colIdx] > maxOverstepValue:
                    maxOverstepValue = image[rowIdx, colIdx]
                state = State.upsideDown


        if state == State.findBorder:
            if image[rowIdx, colIdx] < threshold:
                colStep *= -1
                colBorder = interpolatedLength - 1 - colBorder
                state = State.findBorderReverse
            else:
                if colIdx != colBorder:
                    colIdx += colStep
                else:
                    overstep = True
                    if image[rowIdx, colIdx] > maxOverstepValue:
                        maxOverstepValue = image[rowIdx, colIdx]

                    colStep *= -1
                    colBorder = interpolatedLength - 1 - colBorder
                    state = State.move

        if state == State.findBorderReverse:
            try:
                if image[rowIdx, colIdx] < threshold:
                    if colIdx != colBorder:
                        colIdx += colStep
                    else:
                        state = State.upsideDown

                else:
                    if len(border) > 2:
                        distLastSQ = (border[-1][0] - border[-2][0])**2 + (border[-1][1] - border[-2][1])**2
                        distNowSQ = (border[-1][0] - rowIdx)**2 + (border[-1][1] - colIdx)**2
                        if (distNowSQ - distLastSQ*1.5) > 0:
                            state = State.upsideDown
                            continue
                    border.append([rowIdx, colIdx])
                    state = State.move
            except IndexError:
                print rowIdx, colIdx
                raise

        if state == State.upsideDown:
            if bottom == True:
                state = State.end
            else:
                bottom = True
                rowStep = 1
                colStep = -1
                colBorder = 0
                rowBorder = interpolatedLength - 1
                rowIdx = imageCenter[0]
                colIdx = imageCenter[1]
                state = State.move


    return border, trueCenterIndex, maxOverstepValue


def calculateBeamSize(image, density, windowLength,
        beamMajorAxisScale, interpolatedLength = 800, threshold = 0.5):
    border, closestToCenterIndex, overstep = trackBorder(
            image, threshold, density, interpolatedLength)

    if overstep != 0: print 'overstep'
    # np.savetxt('border', border)
    if len(border) < 10:
        print 'too little points:', len(border)
        return 0, 0, 0, overstep

    imageArray = np.array(border) - [0, closestToCenterIndex[1]]
    imageArray[:, 0] = closestToCenterIndex[0] - imageArray[:, 0]
    # np.savetxt('border', imageArray)
    distancesSQ = np.sum(np.square(imageArray), axis=1)
    minDistIndex = np.argmin(distancesSQ)
    minDist = np.sqrt(distancesSQ[minDistIndex])
    minDistVector = imageArray[minDistIndex]

    maxDistIndex = np.argmax(distancesSQ)
    maxDist = np.sqrt(distancesSQ[maxDistIndex])
    maxDistVector = imageArray[maxDistIndex]
    angle = np.arctan2(maxDistVector[0], maxDistVector[1])

    minorAxis  = (windowLength/interpolatedLength*minDist)
    majorAxis  = (windowLength/interpolatedLength*maxDist)


    # if overstep == False:
        # maxDistIndex = np.argmax(distancesSQ)
        # maxDist = np.sqrt(distancesSQ[maxDistIndex])
        # maxDistVector = imageArray[maxDistIndex]
        # angle = np.arctan2(maxDistVector[0], maxDistVector[1])

        # minorAxis  = np.rad2deg(windowLength/interpolatedLength*minDist)
        # majorAxis  = np.rad2deg(windowLength/interpolatedLength*maxDist)
    # else:
        # angle = np.pi/2. + np.arctan2(minDistVector[0], minDistVector[1])


        # majorAxis = np.rad2deg(beamMajorAxisScale)
        # minorAxis  = np.rad2deg(windowLength/interpolatedLength*minDist)

    # print majorAxis, minorAxis, threshold
    # majorSigma = normSigma(majorAxis, 0, threshold)
    # minorSigma = normSigma(minorAxis, 0, threshold)

    # majorAxis = normInverse(0.5, 0, majorSigma)
    # minorAxis = normInverse(0.5, 0, minorSigma)

    # print majorAxis, minorAxis
    return majorAxis, minorAxis, angle, overstep


