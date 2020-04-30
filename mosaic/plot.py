import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Circle
import matplotlib.animation as animation
from matplotlib.ticker import FormatStrFormatter, FixedLocator, FixedFormatter

def plotBeamContour(array, center, sideLength, fileName='contour.png', interpolation = True):
    thisDpi = 96.
    matplotlib.rcParams.update({'font.size': 8})
    fig = plt.figure(figsize=(400./thisDpi, 300./thisDpi), dpi=thisDpi)
    # plt.locator_params(axis='x', nbins=3)
    # plt.locator_params(axis='y', nbins=3)
    if type(sideLength) == list:
        plotRange = sideLength
    else:
        halfSideLength = sideLength/2.0
        xStart = (center[0] - halfSideLength)
        xEnd = (center[0] + halfSideLength)
        yStart = (center[1] - halfSideLength)
        yEnd = (center[1] + halfSideLength)
        plotRange = [xStart, xEnd, yStart, yEnd]
    interpolateOption = 'bicubic' if interpolation == True else 'nearest'
    extent = [plotRange[0],plotRange[1],plotRange[3],plotRange[2]]
    plt.imshow(np.fliplr(array),cmap=plt.cm.jet, vmin=0, vmax=1,
            interpolation=interpolateOption, extent=extent, origin='bottom')
    plt.colorbar()
    fig.gca().set_aspect('auto')
    fig.gca().xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    fig.gca().yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    plt.xticks((plotRange[0], center[0], plotRange[1]))
    plt.yticks((plotRange[3], center[1], plotRange[2]), rotation=90, va='center')


    plt.subplots_adjust(left=0.20, right=1.00)
    # axes = plt.gca()
    # axes.set_xlim([x.min(),x.max()])
    # axes.set_ylim([y.min(),y.max()])
    plt.xlabel('Right Ascension')
    plt.ylabel('Declination')
    plt.savefig(fileName, dpi=thisDpi)
    plt.close()

def plotPackedBeam(coordinates, angle, axis1, axis2, boresight, equaltorial_range,
        pixel_range, beamRadius, fileName='pack.png', scope=1.,
        extra_coordinates = [], extra_coordinates_text = [],
        transparency = 1., index = False, HD = True):
    #angle = 180 - angle
    # print("tiling angle: %.2f" % angle)
    # angle = angle - 90
    thisDpi = 96
    matplotlib.rcParams.update({'font.size': 8})
    plt.clf()
    if HD == True:
        fig = plt.figure(figsize=(1600./thisDpi, 1600./thisDpi), dpi=thisDpi)
    else:
        fig = plt.figure(figsize=(400./thisDpi, 400./thisDpi), dpi=thisDpi)
    # plt.axes().set_aspect('equal', 'datalim')
    axis = fig.add_subplot(111, aspect='equal')
    # center = coordinates[0]
    for idx in range(len(coordinates)):
        coord = coordinates[idx]
        ellipse = Ellipse(xy=coord, width=2*axis1, height=2*axis2, angle=angle,
                alpha = transparency)
        ellipse.fill = False
        axis.add_artist(ellipse)
        if index == True:
            axis.text(coord[0], coord[1], idx, size=6, ha='center', va='center')
    # print(extra_coordinates)
    # print(coordinates)
    for idx in range(len(extra_coordinates)):
        coord = extra_coordinates[idx]
        circle = Circle(xy=coord, radius = 0.0006)
        axis.add_artist(circle)
        if extra_coordinates_text != []:
            axis.text(coord[0], coord[1], extra_coordinates_text[idx],
                    size=5, ha='center', va='center')

    gridCenter = [0,0]
    circle = Ellipse(xy=gridCenter, width=2*beamRadius, height=2*beamRadius, angle=0)
    circle.fill = False
    axis.add_artist(circle)
    margin = beamRadius*1.3*scope
    axis.set_xlim(gridCenter[0]-margin, gridCenter[0]+margin)
    axis.set_ylim(gridCenter[1]-margin, gridCenter[1]+margin)
    beamNum = len(coordinates)

    xTicks = FixedLocator([pixel_range[0], 0, pixel_range[1]])
    xTicksLabel = FixedFormatter(["{:.2f}".format(equaltorial_range[0]),
                          "{:.2f}".format(boresight[0]),
                          "{:.2f}".format(equaltorial_range[1])])
    axis.xaxis.set_major_formatter(xTicksLabel)
    axis.xaxis.set_major_locator(xTicks)
    axis.xaxis.set_tick_params(tickdir="out", tick2On=False)
    axis.set_xlabel("Right Ascension")

    yTicks = FixedLocator([pixel_range[2], 0, pixel_range[3]])
    yTicksLabel = FixedFormatter(["{:.2f}".format(equaltorial_range[2]),
                          "{:.2f}".format(boresight[1]),
                          "{:.2f}".format(equaltorial_range[3])])
    axis.yaxis.set_major_formatter(yTicksLabel)
    axis.yaxis.set_major_locator(yTicks)
    axis.yaxis.set_tick_params(tickdir="out", tick2On=False)
    axis.set_ylabel("Declination")
    plt.yticks(axis.get_yticks(), visible=True, rotation="vertical")

    fig.tight_layout()

    plt.savefig(fileName, dpi=thisDpi)
    plt.close()

def plotBeamWithFit(array, center, sideLength, widthH, widthV, angle,
        fileName='contourfit.png', interpolation = True):
    thisDpi = 96.
    matplotlib.rcParams.update({'font.size': 8})
    fig = plt.figure(figsize=(400./thisDpi, 300./thisDpi), dpi=thisDpi)
    axis = fig.gca()
    if type(sideLength) == list:
        plotRange = sideLength
    else:
        halfSideLength = sideLength/2.0
        xStart = (center[0] - halfSideLength)
        xEnd = (center[0] + halfSideLength)
        yStart = (center[1] - halfSideLength)
        yEnd = (center[1] + halfSideLength)
        plotRange = [xStart, xEnd, yStart, yEnd]
    interpolateOption = 'bicubic' if interpolation == True else 'nearest'
    ims = axis.imshow(np.fliplr(array), cmap=plt.cm.jet, vmin=0, vmax=1,
            # interpolation=interpolateOption, extent=plotRange)
            interpolation=interpolateOption, aspect = 'equal', origin='bottom')

    imageShape = array.shape
    gridCenter = ((imageShape[1]/2.0 - 1), (imageShape[0]/2.0))
    # print("plot angle: %.2f" % angle)
    ellipse = Ellipse(gridCenter, width=2*widthH, height=2*widthV, angle= angle)
    ellipse.fill = False
    axis.add_artist(ellipse)

    fig.colorbar(ims)
    axis.set_aspect('auto')
    xTicks = FixedLocator([0, gridCenter[0], imageShape[1]-1])
    xTicksLabel = FixedFormatter(["{:.2f}".format(plotRange[0]),
                          "{:.2f}".format(center[0]),
                          "{:.2f}".format(plotRange[1])])
    axis.xaxis.set_major_formatter(xTicksLabel)
    axis.xaxis.set_major_locator(xTicks)
    axis.set_xlabel("RA")

    yTicks = FixedLocator([0, gridCenter[1], imageShape[0]-1])
    yTicksLabel = FixedFormatter(["{:.2f}".format(plotRange[3]),
                          "{:.2f}".format(center[1]),
                          "{:.2f}".format(plotRange[2])])
    axis.yaxis.set_major_formatter(yTicksLabel)
    axis.yaxis.set_major_locator(yTicks)
    axis.yaxis.set_tick_params(rotation=90)
    axis.set_ylabel("DEC")




    # plt.subplots_adjust(left=0.20, right=1.00)
    # axes = plt.gca()
    # axes.set_xlim([x.min(),x.max()])
    # axes.set_ylim([y.min(),y.max()])
    # plt.xlabel('Right Ascension')
    # plt.ylabel('Declination')
    plt.savefig(fileName, dpi=thisDpi)
    plt.close()


def plotBeamFit(sideLength, center, ellipseCenter, angle, axis1, axis2, fileName='fit.png'):
    halfSideLength = sideLength/2.0
    xMin = center[0] - halfSideLength
    xMax = center[0] + halfSideLength
    yMin = center[1] - halfSideLength
    yMax = center[1] + halfSideLength
    thisDpi = 96
    matplotlib.rcParams.update({'font.size': 8})
    plt.clf()
    fig = plt.figure(figsize=(400./thisDpi, 300./thisDpi), dpi=thisDpi)
    plt.subplots_adjust(right=0.92)
    ax = fig.add_subplot(111, aspect='equal')
    # step=sideLength/20.0
    # trueCenter = [center[0] + step/2., center[1] - step/2.]
    ellipse = Ellipse(xy=ellipseCenter, width=2*axis1, height=2*axis2, angle=angle)
    ellipse.fill = False
    ax.add_artist(ellipse)
    radius = max(axis1, axis2)
    ax.set_xlim(xMin, xMax)
    ax.set_ylim(yMin, yMax)
    # ax.get_xaxis().set_visible(False)
    # ax.get_yaxis().set_visible(False)
    # ax.patch.set_visible(False)
    # fig.patch.set_visible(False)
    ax.axis('off')
    plt.savefig(fileName, transparent=True, dpi=thisDpi)
    plt.close()

def plot_overlap(overlapTables, mode, fileName):
    overlapTables = np.array(overlapTables)
    if isinstance(overlapTables[0][0], np.ndarray):
        animate = True
    else:
        animate = False

    pointSize = 10
    thisDpi = 96
    matplotlib.rcParams.update({'font.size': 8})
    plt.clf()
    fig = plt.figure(figsize=(640./thisDpi, 480./thisDpi), dpi=thisDpi)
    plt.title("overlap changes with time")
    if animate == True:
        overlapTable0 = overlapTables[0]
    else:
        overlapTable0 = overlapTables
    if mode == "counter":
        image = plt.imshow(overlapTable0, cmap=plt.cm.jet, vmin=0, vmax=2)
    else:
        image = plt.imshow(overlapTable0, cmap=plt.cm.jet, vmin=0)
    axis = fig.gca()
    axis.set_aspect('equal', 'datalim')
    plt.colorbar()

    if animate == True:
        def animator(data):
            image.set_data(data)
            image.set_clim(0, 2)

        mov = animation.FuncAnimation(fig, animator, frames=overlapTables)

        mov.save(fileName, dpi=thisDpi, writer='imagemagick')

    else:
        plt.savefig(fileName, dpi=thisDpi)

    plt.close()

def plot_interferometry(antennas, center, horizons, fileName='horizon.gif'):

    def angleToCartesian(angleDeg, radius):
        angle = np.deg2rad(angleDeg)
        if angle < np.pi/2:
            x = np.sin(angle)*radius
            y = np.cos(angle)*radius
        elif angle > np.pi/2 and angle < np.pi:
            x = np.sin(angle)*radius
            y = np.cos(angle)*radius
        elif angle > np.pi and angle < np.pi*1.5:
            x = np.sin(angle)*radius
            y = np.cos(angle)*radius
        elif angle > np.pi*1.5 and angle < np.pi*2:
            x = np.sin(angle)*radius
            y = np.cos(angle)*radius

        return x, y

    antennas = np.array(antennas)
    antX = antennas[:,1] - center[1]
    antY = antennas[:,0] - center[0]

    antXMax = np.abs(antX).max()
    antYMax = np.abs(antY).max()

    radius = antXMax*1.2 if antXMax > antYMax else antYMax*1.2

    pointSize = 3
    thisDpi = 96
    # matplotlib.rcParams.update({'font.size': 8})
    plt.clf()

    fig = plt.figure(figsize=(640./thisDpi, 480./thisDpi), dpi=thisDpi)
    ax = fig.gca()
    "cartisian"
    ax.plot([-radius, radius, 0, 0, 0], [0, 0, 0, radius, -radius], c="k", alpha = 0.1)
    "antennas"
    ants = ax.scatter(antX,antY, s=pointSize)
    "horizon"
    horizon = plt.Circle((0,0), radius, fill=False)
    ax.add_artist(horizon )

    horizons = np.array(horizons)
    if isinstance(horizons[0], np.ndarray):
        azi0, alt0 = horizons[0]
    else:
        azi0, alt0 = horizons
    lineStyle = '--' if alt0 < 0 else '-'
    "altitude"
    altRadius = (90-np.abs(alt0))/90.*radius
    altCircle = plt.Circle((0,0), altRadius, fill=False, alpha=0.1, ls=lineStyle)
    ax.add_artist(altCircle)
    "azimuth"
    x, y = angleToCartesian(azi0, radius)
    aziLine, = ax.plot([0,x], [0, y], ls=lineStyle, alpha=0.1)
    "source"
    starColor = 'black' if alt0 < 0 else 'blue'
    starX, starY = angleToCartesian(azi0, altRadius)
    star, = ax.plot(starX, starY, marker='*')
    star.set_color(starColor)

    ax.set_aspect('equal', 'datalim')
    ax.set_xlim([-radius*1.3, +radius*1.3])
    ax.set_ylim([-radius*1.3, +radius*1.3])

    plt.legend([star, ants, horizon], ["source", "antennas", "horizon"])


    if isinstance(horizons[0], np.ndarray):
        def animator(horizon):
            azi, alt = horizon
            lineStyle = '--' if alt < 0 else '-'
            altRadius = (90-np.abs(alt))/90.*radius
            altCircle.set_radius(altRadius)
            x, y = angleToCartesian(azi, radius)
            aziLine.set_data([0,x], [0, y])
            starX, starY = angleToCartesian(azi, altRadius)
            starColor = 'black' if alt < 0 else 'blue'
            star.set_data(starX, starY)
            star.set_color(starColor)
            # return (star)
            return (aziLine, altCircle, star)

        mov = animation.FuncAnimation(fig, animator, frames=horizons)
        mov.save(fileName, dpi=thisDpi, writer='imagemagick')
    else:
        plt.savefig(fileName, dpi=thisDpi)

    plt.close()

def plot_beam_shape(array, center, sideLength, ellipseCenter, axis1, axis2, angle, fileName='contour.png', interpolation = True):
    thisDpi = 96.
    matplotlib.rcParams.update({'font.size': 8})
    fig = plt.figure(figsize=(400./thisDpi, 300./thisDpi), dpi=thisDpi)
    halfSideLength = sideLength/2.0
    xStart = (center[0] - halfSideLength)
    xEnd = (center[0] + halfSideLength)
    yStart = (center[1] - halfSideLength)
    yEnd = (center[1] + halfSideLength)
    plotRange = [xStart, xEnd, yStart, yEnd]
    interpolateOption = 'bicubic' if interpolation == True else 'nearest'
    plt.imshow(array,cmap=plt.cm.jet, vmin=0, vmax=1, interpolation=interpolateOption, extent=plotRange)
    plt.colorbar()
    ax = fig.gca()
    ax.set_aspect('equal', 'datalim')
    ellipse = Ellipse(xy=ellipseCenter, width=2*axis1, height=2*axis2, angle=angle)
    ellipse.fill = False
    ax.add_artist(ellipse)
    plt.savefig(fileName, dpi=thisDpi)
    plt.close()