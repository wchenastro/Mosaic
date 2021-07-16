import numpy as np
from scipy import interpolate
from scipy.optimize import curve_fit
from mosaic.utilities import normInverse
import logging
from matplotlib import pyplot as plt
from mosaic.fitEllipse import fit_ellipse
from matplotlib.patches import Ellipse

loggerFormat = '%(asctime)-15s  %(filename)s  %(message)s'
logging.basicConfig(format = loggerFormat, level=logging.WARNING)
logger = logging.getLogger(__name__)


def calculateBeamOverlaps(ellipseCenters, radius, majorAxis, minorAxis, rotation, overlap, mode, sideLength=None):

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
    if sideLength == None:
        halfSidelength = longAxis*3.
    else:
        halfSidelength = sideLength/2.
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
        move, addRow, findBorder, findBorderReverse, upsideDown, end = list(range(6))

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
    # filling = threshold * 0.66666666666666
    filling = 0
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

    # np.save("/tmp/trackimage", image)
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
        axis2  = (windowLength/interpolatedLength*minDist)
        axis1  = (windowLength/interpolatedLength*maxDist)
    else:
        widthH, widthV, angle = fitEllipse(iterpolatedImage)
        axis1  = (windowLength/interpolatedLength*widthH)
        axis2  = (windowLength/interpolatedLength*widthV)
        # print("fit angle: %.2f" % np.rad2deg(angle))
        #angle = angle + np.pi/2.0
        #angle = np.pi - angle
        # if abs(angle) > 360.:
          # angle = angle % 360.
        # angle = np.pi - (angle + np.pi/2.0)

    # print majorAxis, minorAxis
    return axis1, axis2, np.rad2deg(angle), overstep

def fitEllipse(image):
    def RotatedGaussian2DPDF(xy, xMean, yMean, xSigma, ySigma, angle):
        x,y = xy
        xSigma_square = xSigma**2
        ySigma_square = ySigma**2
        cos_angle_square = (np.cos(angle))**2
        sin_angle_square = (np.sin(angle))**2
        a = cos_angle_square/(2*xSigma_square) + sin_angle_square/(2*ySigma_square)
        b = - np.sin(2*angle)/(4*xSigma_square) + np.sin(2*angle)/(4*ySigma_square)
        c = sin_angle_square/(2*xSigma_square) + cos_angle_square/(2*ySigma_square)

        xMinusxMean = x-xMean
        yMinusyMean = y-yMean
        values = np.exp(-(a*xMinusxMean**2 + 2*b*(xMinusxMean)*(yMinusyMean) + c*yMinusyMean**2))
        return values

    yData = image
    dataShape = yData.shape

    X,Y = np.meshgrid(np.linspace(0, dataShape[0]-1, dataShape[0]),
            np.linspace(0, dataShape[1]-1, dataShape[1]))

    X = X[yData != 0]
    Y = Y[yData != 0]
    yData = yData[yData != 0]

    initial_guess = (dataShape[1]/2, dataShape[0]/2, 160, 160, 0)

    # paras_bounds = ([300, 300, 10, 10, -2*np.inf], [500, 500, dataShape[0], dataShape[0], 2*np.inf])
    paras_bounds = ([360, 360, 10, 10, -2*np.inf], [440, 440, dataShape[0], dataShape[0], 2*np.inf])
    popt, pcov = curve_fit(RotatedGaussian2DPDF, (X,Y), yData,
            p0=initial_guess, bounds=paras_bounds)

    centerX, centerY, sigmaH, sigmaV, angle = popt

    widthH = normInverse(0.5, 0, sigmaH)
    widthV = normInverse(0.5, 0, sigmaV)

    return widthH, widthV, angle

def createBeamshapeModel(originalImage, density, windowLength, interpolatedLength = 800):

    interpolater = interpolate.interp2d(range(density), range(density), originalImage, kind='cubic')
    image = interpolater(np.linspace(0, density, interpolatedLength, endpoint = False),
                np.linspace(0, density, interpolatedLength, endpoint = False))

    samples = []

    # levels = np.arange(0.1, 0.9 + 0.025, 0.025).tolist()
    levels = np.linspace(0.1, 0.9, 1 + int((0.9 - 0.1)/0.025)).tolist()
    # print(levels)
    # levels = [0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975]
    # levels = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
    plot = False
    if plot == True:
        thisDpi = 96
        fig = plt.figure(figsize=(1200./thisDpi,1200./thisDpi), dpi=thisDpi)
        ax0 = plt.subplot2grid((6, 6), (3, 0), colspan=3, rowspan=3)
        ax1 = plt.subplot2grid((6, 6), (3, 3), colspan=3, rowspan=3)
        ax2 = plt.subplot2grid((6, 6), (0, 0), colspan=6)
        ax3 = plt.subplot2grid((6, 6), (1, 0), colspan=6)
        ax4 = plt.subplot2grid((6, 6), (2, 0), colspan=6)
        ax0.imshow(image, cmap="jet")
        ax0.title.set_text('Contour')
        ax1.imshow(image, cmap="jet")
        ax1.title.set_text('Fit')
    else:
        fig, ax0 = plt.subplots()
    contour = ax0.contour(image, levels)
    dataShape = image.shape
    center = np.unravel_index(image.argmax(), dataShape)
    count = 0.1
    for segs, coll in zip(contour.allsegs, contour.collections):
        count += 0.025
        powerline = None
        fulllength = 0
        segLength = len(segs)
        if segLength > 1:
            for seg in segs:
                fulllength += len(seg)
            minimalLength = int(fulllength * 0.2)
            paths = coll.get_paths()
            pathIndexToDel = []
            for segIndex, seg in enumerate(segs):
                if len(seg) < minimalLength:
                    pathIndexToDel.append(segIndex)
                    continue
                if powerline is None:
                    powerline = seg
                else:
                    powerline = np.concatenate((powerline, seg))
            pathIndexToDel.reverse()
            for index in pathIndexToDel:
                del(paths[index])
                # print("small contour of length {} at level {} delelted".format(
                #           len(segs[index]), levels[len(samples)]))
        elif segLength == 1:
            powerline = np.array(segs[0])
        else:
            powerline = None

        if powerline is None:
            if count > 0.2:
                logger.warning('level {} countour is None!'.format(count))
            para = [np.nan, np.nan, 0, 0, np.nan]
        else:
            para = fitContour(powerline)
        samples.append(para)
        if plot == True:
            ellipse = Ellipse(center, width=2*para[0], height=2*para[1], angle=np.rad2deg(para[4]))
            ellipse.fill = False
            ax1.add_artist(ellipse)

    ### check NaN ##
    samples = np.array(samples).T
    for sampleIndex in np.arange(len(samples)):
        nans = np.isnan(samples[sampleIndex])
        if np.any(nans):
            levelNumber = levels[np.squeeze(np.argwhere(nans)[0][0])]
            if levelNumber > 0.2:
                logger.warning('level {} have NaN value!'.format(levelNumber))
            sample = samples[sampleIndex]
            sample[nans]= np.interp(np.array(levels)[nans], np.array(levels)[~nans], sample[~nans])
            samples[sampleIndex] = sample
    samples = samples.T

    ## extrapolate the overlap ratio of 1
    # levels.append(1.0)
    # samples.append([0, 0, samples[-1][2], samples[-1][3], samples[-1][4]])
    ## extrapolate the overlap ratio of 0
    # levels.insert(0, 0.0)
    # samples.insert(0, samples[0])
    # samples = np.array(samples)
    # samples[:, 0] = (1.*windowLength/interpolatedLength*samples[:, 0])
    # samples[:, 1] = (1.*windowLength/interpolatedLength*samples[:, 1])
    # samples[:, 4] = np.rad2deg(samples[:, 4])
    # samples[0, 0:2] = [0.89/2, 0.89/2]


    # samples = np.array(samples)
    samples[:, 0] = (1.*windowLength/interpolatedLength*samples[:, 0])
    samples[:, 1] = (1.*windowLength/interpolatedLength*samples[:, 1])

    angles = samples[:, 4]
    for i in range(1, len(angles)):
        if angles[i] - angles[i-1] > np.pi * 0.8:
            angles[i] -= np.pi
        elif angles[i] - angles[i-1] < -np.pi * 0.8:
            angles[i] += np.pi

    samples[:, 4] = np.rad2deg(angles)

    interpMethod ="cubic" # "quadratic, slinear, cubic"
    majorInterp = interpolate.interp1d(levels, samples[:, 0], kind=interpMethod)
    minorInterp = interpolate.interp1d(levels, samples[:, 1], kind=interpMethod)
    angleInterp = interpolate.interp1d(levels, samples[:, 4], kind=interpMethod)

    interval = 0.001
    # levelInterp = np.arange(levels[0], levels[-1]+interval, interval).tolist()
    levelInterp  = np.linspace(levels[0], levels[-1],
            1 + int((levels[-1] - levels[0])/interval)).tolist()
    # levelInterp[-1] = np.round(levelInterp[-1], 7)

    majorInterped = majorInterp(levelInterp).tolist()
    minorInterped = minorInterp(levelInterp).tolist()
    angleInterped = angleInterp(levelInterp).tolist()
    # angleInterped = interpolate.splev(levelInterp, tck, der=0)

    ### extrapolate the overlap ratio of 1
    levelInterp.append(1.0)
    majorInterped.append(0)
    minorInterped.append(0)
    angleInterped.append(angleInterped[-1])
    ### extrapolate the overlap ratio of 0
    levelInterp.insert(0, 0.0)
    majorInterped.insert(0, 0.89/2)
    minorInterped.insert(0, 0.89/2)
    angleInterped.insert(0, angleInterped[0])

    interval = 0.00001
    levelLinearInterp = np.arange(levelInterp[0], levelInterp[-1]+interval, interval)
    majorInterped = np.interp(levelLinearInterp, levelInterp, majorInterped)
    minorInterped = np.interp(levelLinearInterp, levelInterp, minorInterped)
    angleInterped = np.interp(levelLinearInterp, levelInterp, angleInterped)
    # print(len(levelInterp), len(levelLinearInterp))


    if plot == True:
        ax2.plot(levels, samples[:, 0])
        ax3.plot(levels, samples[:, 1])
        ax4.plot(levels, samples[:, 4])

        ax2.plot(levelLinearInterp, majorInterped)
        ax3.plot(levelLinearInterp, minorInterped )
        # ax4.plot(levelInterp, np.rad2deg(angleInterp(levelInterp)))
        ax4.plot(levelLinearInterp, angleInterped)

        ax2.set_ylabel('Major', size=15)
        ax2.set_xticklabels([])
        ax2.tick_params(axis='y', labelsize=15 )
        ax2.set_ylim(0, max(samples[1:, 0])*1.2)
        ax3.set_ylabel('Minor', size=15)
        ax3.tick_params(axis='y', labelsize=15 )
        ax3.set_xticklabels([])
        ax3.set_ylim(0, max(samples[1:, 1])*1.2)
        ax4.set_ylabel('Orientation', size=15)
        ax4.set_xlabel('Overlap ratio', size=15)
        ax4.tick_params(axis='x', labelsize=15 )

        fig.subplots_adjust(hspace=0)
        fig.tight_layout()
        fig.savefig("plots/levels.png")

    plt.close()

    beamshapeModel = np.array([levelLinearInterp, majorInterped, minorInterped, angleInterped]).T

    return beamshapeModel

def fitContour(contour):
    major, minor, centerX, centerY, angle = fit_ellipse(contour[:,0], contour[:,1])
    return major, minor, centerX, centerY, angle
