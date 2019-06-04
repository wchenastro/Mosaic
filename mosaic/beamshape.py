import numpy as np
from scipy import interpolate
from scipy.optimize import curve_fit
from utilities import normInverse
import logging

loggerFormat = '%(asctime)-15s  %(filename)s  %(message)s'
logging.basicConfig(format = loggerFormat, level=logging.WARNING)
logger = logging.getLogger(__name__)


def calculateBeamOverlaps(ellipseCenters, radius, majorAxis, minorAxis, rotation, overlap, mode):

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

    rotation = np.deg2rad(rotation)
    # rotated = 0/(3600*24.) * 2 * np.pi
    # rotation += rotated
    longAxis = majorAxis if majorAxis > minorAxis else minorAxis
    gridNum = 1000
    # print 0.3*radius, longAxis
    # halfSidelength = 0.15 * radius
    halfSidelength = longAxis*3.
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
    widthH = normInverse(overlap, 0, sigmaH)
    widthV = normInverse(overlap, 0, sigmaV)


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
            counts = isInsideEllips(ellipseCenter, widthH, widthV, rotation, gridX, gridY)
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

    if mode == 1:
        return overlapCounter
    elif mode == 2:
        return overlapHeater
    elif mode == 3:
        return overlapCounter, overlapHeater
    # np.save('overlapCounter', overlapCounter)


def trackBorder(image_orig, threshold = 0.3, density = 20, interpolatedLength = 800):

    interpolater = interpolate.interp2d(range(density), range(density), image_orig ,kind='cubic')
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
    right_border = int(interpolatedLength - 1)
    left_border = 0
    colBorder = right_border
    rowBorder = 0
    # horizonDirection = 1 # left
    bottom = False
    overstep = False
    maxOverstepValue = 0
    offset = 0.1
    filling = threshold * 0.66666666666666
    # closestToCenter = 1
    # closestToCenterIndex = []
    states = {
        State.move:"move",
        State.addRow:"addRow",
        State.findBorder:"findBorder",
        State.findBorderReverse:"findBorderReverse",
        State.upsideDown:"upsideDown",
        State.end:"end"}

    logs = []
    leftStep=True
    rowCenter = colIdx
    edges = [0,0]
    maxHalfWdith = interpolatedLength/2
    first = True
    stepCounter = 0
    while state != State.end:
        # log = "%s %d %d %d %d %s\n" % (states[state], rowIdx, colIdx, colStep, rowCenter, str(edges))
        # logs.append(log)
        if state == State.move:
            colIdx = rowCenter
            stepCounter = 0
            while colIdx != colBorder:
                if image[rowIdx, colIdx] < threshold or stepCounter >= maxHalfWdith:
                    border.append([rowIdx, colIdx])
                    if colStep == 1:
                        image[rowIdx, colIdx:] = filling
                        edges[1] = colIdx
                    else:
                        image[rowIdx, :colIdx] = filling
                        edges[0] = colIdx
                    break
                colIdx += colStep
                stepCounter += abs(colStep)
            if colIdx == colBorder:
                overstep = True
                if image[rowIdx, colIdx] > maxOverstepValue:
                    maxOverstepValue = image[rowIdx, colIdx]
                if colStep == 1:
                    edges[1] = colIdx
                else:
                    edges[0] = colIdx

            if leftStep == True:
                leftStep = False
            else:
                leftStep = True
                if first == True:
                    first = False
                    widthLeft = edges[1] - rowCenter
                    widthRight = rowCenter - edges[0]
                    smallerWidth = widthRight if widthLeft > widthRight else widthLeft
                    maxHalfWdith = smallerWidth
                else:
                    halfWidth = int(round((edges[1] - edges[0])/2.0))
                    if halfWidth < maxHalfWdith:
                        maxHalfWdith = halfWidth

                state = State.addRow
                if edges[0] != edges[1]:
                    rowCenter = edges[0] + np.argmax(image[rowIdx, edges[0]:edges[1]])
                    if image[rowIdx, edges[0]:edges[1]].max() < (threshold + 0.07):
                        state = State.upsideDown
            colStep *= -1
            if colStep == 1:
                colBorder = right_border
            else:
                colBorder = left_border

        elif state == State.addRow:
            if edges[0] == edges[1]:
                state = State.upsideDown
                border.pop()
            elif rowIdx != rowBorder:
                rowIdx += rowStep
                state = State.move
            else:
                overstep = True
                if image[rowIdx, colIdx] > maxOverstepValue:
                    maxOverstepValue = image[rowIdx, colIdx]
                state = State.upsideDown


        elif state == State.upsideDown:
            if bottom == True:
                state = State.end
            else:
                bottom = True
                rowStep = 1
                colStep = 1
                colBorder = right_border
                rowBorder = interpolatedLength - 1
                rowIdx = int(imageCenter[0])
                colIdx = int(imageCenter[1])
                rowCenter = colIdx
                first = True
                maxHalfWdith = interpolatedLength/2
                state = State.move

    # with open("/tmp/stateLog", 'w') as stateLogFile:
        # stateLogFile.writelines(logs)

    border = np.array(border)
    if border != []:
        topRow = border[:,0].max()
        bottomRow = border[:,0].min()
        image[topRow+1:, :] = filling
        image[:bottomRow, :] = filling

        # np.save("/tmp/trackimage", image)


    return border, trueCenterIndex, maxOverstepValue, image


def calculateBeamSize(image, density, windowLength,
        beamMajorAxisScale, interpolatedLength = 800, threshold = 0.2, fit=False):
    border, closestToCenterIndex, overstep, iterpolatedImage = trackBorder(
            image, threshold, density, interpolatedLength)

    # if overstep != 0: print 'overstep'
    # np.savetxt('border', border)
    if len(border) < 10:
        logger.info('less then 10 points in the border tracking:')
        return 0, 0, 0, overstep

    if fit == False:
        imageArray = np.array(border) - [0, closestToCenterIndex[1]]
        imageArray[:, 0] = closestToCenterIndex[0] - imageArray[:, 0]
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
    else:
        widthH, widthV, angle = fitEllipse(np.flipud(iterpolatedImage))
        minorAxis  = (windowLength/interpolatedLength*widthH)
        majorAxis  = (windowLength/interpolatedLength*widthV)
        angle = angle + np.pi/2.0
        if abs(angle) > 360.:
            angle = angle % 360.
        # angle = np.pi - (angle + np.pi/2.0)

    # print majorAxis, minorAxis
    return majorAxis , minorAxis, np.rad2deg(angle), overstep


def fitEllipse(image):
    def RotatedGaussian2DPDF((x, y), xMean, yMean, xSigma, ySigma, angle):
        angle = -(angle - np.pi)
        a = np.power(np.cos(angle), 2)/(2*xSigma**2) + np.power(np.sin(angle), 2)/(2*ySigma**2)
        b = - np.sin(2*angle)/(4*xSigma**2) + np.sin(2*angle)/(4*ySigma**2)
        c = np.power(np.sin(angle), 2)/(2*xSigma**2) + np.power(np.cos(angle), 2)/(2*ySigma**2)

        values = np.exp(-(a*np.power(x-xMean, 2) + 2*b*(x-xMean)*(y-yMean) + c*np.power(y-yMean, 2)))
        return values.ravel()

    yData = image
    dataShape = yData.shape
    yData = yData.ravel()

    grids = np.meshgrid(np.linspace(0, dataShape[0]-1, dataShape[0]),
            np.linspace(0, dataShape[1]-1, dataShape[1]))

    initial_guess = (dataShape[1]/2, dataShape[0]/2, 20,40,0)

    paras_bounds = ([300, 300, 10, 10, -2*np.inf], [500, 500, dataShape[0], dataShape[0], 2*np.inf])
    popt, pcov = curve_fit(RotatedGaussian2DPDF, grids, yData,
            p0=initial_guess, bounds=paras_bounds)

    centerX, centerY, sigmaH, sigmaV, angle = popt

    widthH = normInverse(0.5, 0, sigmaH)
    widthV = normInverse(0.5, 0, sigmaV)

    return widthH, widthV, angle


