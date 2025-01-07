import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Circle, RegularPolygon, Polygon
import matplotlib.animation as animation
from matplotlib.ticker import FormatStrFormatter, FixedLocator, FixedFormatter
from astropy import wcs

def plotBeamContour(array, center, sideLength, fileName='contour.png',
        colormap = False, interpolation = True):
    thisDpi = 96.
    matplotlib.rcParams.update({'font.size': 8})
    if(colormap == True):
        fig = plt.figure(figsize=(400./thisDpi, 300./thisDpi), dpi=thisDpi)
    else:
        fig = plt.figure(figsize=(400./thisDpi, 400./thisDpi), dpi=thisDpi)
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
    if(colormap == True):
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

def plotPackedBeam(coordinates, angle, axis1, axis2, boresight,
        tiling_meta, fileName='pack.png', scope=1., show_axis = True,
        extra_coordinates = [], extra_coordinates_text = [], subGroup = [],
        transparency = 1., edge = True, index = False, HD = True, raw = False,
        output_format='png'):


    step = 1/10000000000.
    wcs_properties = wcs.WCS(naxis=2)
    wcs_properties.wcs.crpix = [0, 0]
    wcs_properties.wcs.cdelt = [-step, step]
    wcs_properties.wcs.crval = boresight
    wcs_properties.wcs.ctype = ["RA---TAN", "DEC--TAN"]

    center = boresight
    resolution = step


    thisDpi = 96
    matplotlib.rcParams.update({'font.size': 8})
    plt.clf()
    if HD == True:
        width = 1600.
        extra_source_size = 20
        extra_source_text_size = 30
        ticklabel_size = 20
    else:
        width = 400.
        extra_source_size = 3
        extra_source_text_size = 12
        ticklabel_size = 7


    fig = plt.figure(figsize=(width/thisDpi, width/thisDpi), dpi=thisDpi)

    axis = fig.add_subplot(111,aspect='equal', projection=wcs_properties)

    beam_coordinate = np.array(coordinates)/resolution
    # beam_coordinate += [center[0] - 1, 0]
    # beam_coordinate = [0, center[1] - 1] + [1, -1]*beam_coordinate


    for idx in range(len(beam_coordinate)):
        coord = beam_coordinate[idx]
        ellipse = Ellipse(xy=coord,
                width=2.*axis1/resolution,height=2.*axis2/resolution, angle=angle,
                alpha = transparency)
        ellipse.fill = False
        axis.add_artist(ellipse)
        if index == True:
            axis.text(coord[0], coord[1], idx, size=6, ha='center', va='center')

    for group in subGroup:
        subCoordinates, subAngle, subAxis1, subAxis2 = group
        for idx in range(len(subCoordinates)):
            coord = subCoordinates[idx]/resolution
            ellipse = Ellipse(xy=coord,
                    width=2*subAxis1/resolution, height=2*subAxis2/resolution,
                    angle=subAngle, alpha = transparency)
            ellipse.fill = False
            axis.add_artist(ellipse)

    """
    for idx in range(len(extra_coordinates)):
        coord = extra_coordinates[idx]/resolution
        circle = Circle(xy=coord, radius = 0.0006)
        axis.add_artist(circle)
    """
    extra_coordinates = np.array(extra_coordinates)
    if extra_coordinates.size != 0:
        coord = extra_coordinates/resolution
        axis.scatter(coord[:,0], coord[:,1], s=extra_source_size, color='k')

        for idx in range(len(extra_coordinates_text)):
            axis.text(coord[idx][0], coord[idx][1], extra_coordinates_text[idx],
                    size=extra_source_text_size)

    gridCenter = center
    tilingScale = tiling_meta["scale"]
    if edge == True:
        shape = tiling_meta["shape"]
        if shape == "circle":
            edge = Circle(xy=gridCenter, radius = tilingScale/resolution,
                    alpha=0.1, linewidth = 3)
        elif shape == "ellipse":
            a, b, angle = tiling_meta["parameter"]
            edge = Ellipse(xy=gridCenter,
                    width=2*a/resolution, height=2*b/resolution,
                    angle=angle, alpha=0.1, linewidth = 3)
        elif shape == "hexagon":
            if tiling_meta["parameter"] is not None:
                angle = tiling_meta["parameter"][1]
            else:
                angle = 0
            edge = RegularPolygon(xy=gridCenter, numVertices=6,
                    radius = tilingScale/resolution, orientation=np.deg2rad(60.-angle),
                    alpha=0.1, linewidth = 3)
        elif shape == "polygon":
            vertices = tiling_meta["parameter"]
            edge = Polygon(xy=np.array(vertices)/resolution,
                    closed=True, alpha=0.1, linewidth = 3)

        if shape == "annulus":
            for annulus in tiling_meta["parameter"]:
                annulus_type = annulus[0]
                if annulus_type == "polygon":
                    vertices= annulus[1]
                    edge = Polygon(xy=np.array(vertices)/resolution,
                            closed=True, alpha=0.1, linewidth = 3)
                elif annulus_type == "ellipse":
                    a, b , angle = annulus[1]
                    edge = Ellipse(xy=gridCenter,
                            width=2*a/resolution, height=2*b/resolution,
                            angle=angle, alpha=0.1, linewidth = 3)
                edge.fill = False
                axis.add_artist(edge)
        else:
            edge.fill = False
            axis.add_artist(edge)

    margin = tilingScale/step*1.3*scope
    axis.set_xlim(gridCenter[0]-margin, gridCenter[0]+margin)
    axis.set_ylim(gridCenter[1]-margin, gridCenter[1]+margin)

    fig.tight_layout()

    if show_axis != True:
        plt.xticks([])
        plt.yticks([])
    else:
        ra = axis.coords[0]
        dec = axis.coords[1]
        ra.set_ticklabel(size=ticklabel_size)
        dec.set_ticklabel(size=ticklabel_size, rotation="vertical")
        dec.set_ticks_position('l')
        ra.set_major_formatter('hh:mm:ss')
        ra.set_ticks_position('b')
        ra.set_axislabel("RA", size=20)
        dec.set_major_formatter('dd:mm:ss')
        dec.set_axislabel("DEC", size=20)
        plt.subplots_adjust(left=0.04, bottom=0.05, right=0.98, top=0.96,
                wspace=0, hspace=0)

    if isinstance(fileName, str) is True:
        plt.savefig(fileName, dpi=thisDpi)
    else:
        plt.savefig(fileName, dpi=thisDpi, format=output_format)

    plt.close()



def rotatedEllipseParametric(center, major, minor, angle, parameter):

    x = center[0] + major*np.cos(angle)*np.cos(parameter) - minor*np.sin(angle)*np.sin(parameter)
    y = center[1] + major*np.sin(angle)*np.cos(parameter) + minor*np.cos(angle)*np.sin(parameter)

    return np.array([x,y]).T

def plotBeamWithFit2(array, center, sideLength, widthH, widthV, angle,
        fileName='contourfit.png', colormap = False, interpolation = True):
    thisDpi = 96.
    matplotlib.rcParams.update({'font.size': 8})
    if(colormap == True):
        fig = plt.figure(figsize=(400./thisDpi, 300./thisDpi), dpi=thisDpi)
    else:
        fig = plt.figure(figsize=(400./thisDpi, 400./thisDpi), dpi=thisDpi)

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
            interpolation=interpolateOption, aspect = 'equal', origin='lower')

    imageShape = array.shape
    gridCenter = ((imageShape[1]/2.0 - 1), (imageShape[0]/2.0))
    # print("plot angle: %.2f" % angle)
    ellipse = Ellipse(gridCenter, width=2*widthH, height=2*widthV, angle= angle)
    ellipse.fill = False
    axis.add_artist(ellipse)

    if(colormap == True):
        fig.colorbar(ims)
    axis.set_aspect('auto')
    xTicks = FixedLocator([0, gridCenter[0], imageShape[1]-1])
    xTicksLabel = FixedFormatter(["{:.2f}".format(plotRange[0]),
                          "{:.2f}".format(center[0]),
                          "{:.2f}".format(plotRange[1])])
    axis.xaxis.set_major_locator(xTicks)
    axis.xaxis.set_major_formatter(xTicksLabel)
    axis.set_xlabel("RA")

    yTicks = FixedLocator([0, gridCenter[1], imageShape[0]-1])
    yTicksLabel = FixedFormatter(["{:.2f}".format(plotRange[3]),
                          "{:.2f}".format(center[1]),
                          "{:.2f}".format(plotRange[2])])
    axis.yaxis.set_major_locator(yTicks)
    axis.yaxis.set_major_formatter(yTicksLabel)
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

def plotBeamWithFit(array, center, sideLength, widthH, widthV, angle, resolution,
        fileName='contourfit.png', colormap = False, interpolation = True,
        shapeOverlay = False, output_format = 'png'):
    thisDpi = 96.
    matplotlib.rcParams.update({'font.size': 8})

    shape = np.array(array).shape
    step = resolution
    wcs_properties = wcs.WCS(naxis=2)
    wcs_properties.wcs.crpix = [shape[1]/2.-0.5, shape[0]/2.-0.5]
    wcs_properties.wcs.cdelt = [-step, step]
    wcs_properties.wcs.crval = center
    wcs_properties.wcs.ctype = ["RA---TAN", "DEC--TAN"]

    if(colormap == True):
        fig = plt.figure(figsize=(400./thisDpi, 300./thisDpi), dpi=thisDpi)
    else:
        fig = plt.figure(figsize=(400./thisDpi, 400./thisDpi), dpi=thisDpi)

    axis = fig.add_subplot(111,aspect='equal', projection=wcs_properties)
    interpolateOption = 'bicubic' if interpolation == True else 'nearest'
    ims = axis.imshow(array, cmap=plt.cm.jet, vmin=0, vmax=1,
            interpolation=interpolateOption, aspect = 'equal', origin='lower')

    imageShape = array.shape
    if shapeOverlay == True:
        gridCenter = ((imageShape[1]/2.0), (imageShape[0]/2.0))
        ellipse = Ellipse(gridCenter,
                width=2*widthH/resolution, height=2*widthV/resolution, angle=angle)
        ellipse.fill = False
        axis.add_artist(ellipse)

    if(colormap == True):
        fig.colorbar(ims)

    fig.tight_layout()
    axis.set_aspect('auto')
    ra = axis.coords[0]
    dec = axis.coords[1]
    ra.set_major_formatter('hh:mm:ss')
    ra.set_ticklabel(size=8)
    dec.set_ticklabel(size=8, rotation="vertical")
    dec.set_ticks_position('l')
    ra.set_ticks_position('b')
    ra.set_axislabel("RA", size=8)
    dec.set_major_formatter('dd:mm:ss')
    dec.set_axislabel("DEC", size=8)
    plt.subplots_adjust(left=0.10, bottom=0.10, right=0.98, top=0.96,
            wspace=0, hspace=0)

    if isinstance(fileName, str) is True:
        plt.savefig(fileName, dpi=thisDpi)
    else:
        plt.savefig(fileName, dpi=thisDpi, format=output_format)

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

def plot_overlap(overlapTables, mode, fileName, scope = 1., axis = True, title = None):
    overlapTables = np.array(overlapTables)
    if isinstance(overlapTables[0][0], np.ndarray):
        animate = True
    else:
        animate = False

    pointSize = 10
    thisDpi = 96
    matplotlib.rcParams.update({'font.size': 8})
    plt.clf()
    fig = plt.figure(figsize=(640./thisDpi, 640./thisDpi), dpi=thisDpi)
    if title != None:
        plt.title(title)
    if animate == True:
        overlapTable0 = overlapTables[0]
    else:
        overlapTable0 = overlapTables
    if mode == "counter":
        image = plt.imshow(overlapTable0, cmap=plt.cm.jet, vmin=0, vmax=2)
    else:
        image = plt.imshow(overlapTable0, cmap=plt.cm.jet, vmin=0)
    axis = fig.gca()
    if axis != True:
        plt.axis('off')

    length = len(overlapTable0)
    gridCenter = [length/2,length/2]
    margin = int(length/2*scope)
    axis.set_xlim(gridCenter[0]-margin, gridCenter[0]+margin)
    axis.set_ylim(gridCenter[1]-margin, gridCenter[1]+margin)
    axis.set_aspect('auto')


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

def plot_all(interferometry, beamshape, overlap, retile, fileName):
    thisDpi = 96
    matplotlib.rcParams.update({'font.size': 8})
    fig = plt.figure(figsize=(1024./thisDpi, 1024./thisDpi), dpi=thisDpi)
    fig, axs = plt.subplots(2, 2)
    """
    array configuraton
    parameters: interferometry
    """
    def angleToCartesian(angleDeg, radius):
        angle = np.deg2rad(angleDeg)
        if angle <= np.pi/2:
            x = np.sin(angle)*radius
            y = np.cos(angle)*radius
        elif angle > np.pi/2 and angle <= np.pi:
            x = np.sin(angle)*radius
            y = np.cos(angle)*radius
        elif angle > np.pi and angle <= np.pi*1.5:
            x = np.sin(angle)*radius
            y = np.cos(angle)*radius
        elif angle > np.pi*1.5 and angle <= np.pi*2:
            x = np.sin(angle)*radius
            y = np.cos(angle)*radius

        return x, y

    antennas = interferometry[0]
    center = interferometry[1]
    horizons = interferometry[2]
    antennas = np.array(antennas)
    antX = antennas[:,1] - center[1]
    antY = antennas[:,0] - center[0]

    antXMax = np.abs(antX).max()
    antYMax = np.abs(antY).max()

    radius = antXMax*1.2 if antXMax > antYMax else antYMax*1.2

    pointSize = 3
    thisDpi = 96

    ax = axs[0, 0]
    ax.set_title("Array and Source")
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
    ax.legend([star, ants, horizon], ["source", "antennas", "horizon"], loc="upper right")

    """
    beam shape with fit
    parameters: beamshape
    """

    axis = axs[0, 1]
    axis.set_title("Beam shape simulation")
    array, center, sideLength, widthH, widthV, angle, drop, interpolation = beamshape

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
    ims = axis.imshow(array, cmap=plt.cm.jet, vmin=0, vmax=1,
            # interpolation=interpolateOption, extent=plotRange)
            interpolation=interpolateOption, aspect = 'equal', origin = 'lower')

    imageShape = array.shape
    # gridCenter = ((imageShape[1]/2.0 - 1), (imageShape[0]/2.0))
    gridCenter = ((imageShape[1]/2.0), (imageShape[0]/2.0))
    # print("plot angle: %.2f" % angle)
    ellipse = Ellipse(gridCenter, width=2*widthH, height=2*widthV, angle = angle)
    ellipse.fill = False
    axis.add_artist(ellipse)

    # fig.colorbar(ims)
    axis.set_aspect('equal')
    xTicks = FixedLocator([0, gridCenter[0], imageShape[1]-1])
    xTicksLabel = FixedFormatter(["{:.2f}".format(plotRange[0]),
                          "{:.2f}".format(center[0]),
                          "{:.2f}".format(plotRange[1])])
    axis.xaxis.set_major_locator(xTicks)
    axis.xaxis.set_major_formatter(xTicksLabel)
    axis.xaxis.set_tick_params(tickdir="out", tick2On=False)
    axis.set_xlabel("Right Ascension")

    yTicks = FixedLocator([0, gridCenter[1], imageShape[0]-1])
    yTicksLabel = FixedFormatter(["{:.2f}".format(plotRange[3]),
                          "{:.2f}".format(center[1]),
                          "{:.2f}".format(plotRange[2])])
    axis.yaxis.set_major_locator(yTicks)
    axis.yaxis.set_major_formatter(yTicksLabel)
    axis.yaxis.set_tick_params(tickdir="out", tick2On=False)
    axis.set_ylabel("Declination")
    # plt.yticks(axis.get_yticks(), visible=True, rotation="vertical")

    """
    overlap
    parameters: overlap
    """
    axis = axs[1, 0]
    axis.set_title("Tiling Overlap")
    overlapTables, mode, drop, title = overlap
    overlapTables = np.array(overlapTables)
    if isinstance(overlapTables[0][0], np.ndarray):
        animate = True
    else:
        animate = False

    pointSize = 10
    thisDpi = 96
    if title != None:
        pass
        # plt.title(title)
    else:
        pass
    overlapTable0 = overlapTables
    if mode == "counter":
        image = axis.imshow(overlapTable0, cmap=plt.cm.jet, vmin=0, vmax=2)
    else:
        image = axis.imshow(overlapTable0, cmap=plt.cm.jet, vmin=0)
    axis.set_aspect('equal', 'datalim')

    # fig.tight_layout()

    """
    reitle
    parameters: retile
    """
    axis = axs[1, 1]
    axis.set_title("Tiling Efficiency")
    counter, x, xMax, threshold = retile
    counter = np.array(counter)
    axis.set_xlim([x[0], xMax])
    axis.set_ylim([-0.05, 1.05])
    axis.plot(x, counter[:,0], c='r', ls='-',label='overlap')
    axis.plot(x, counter[:,1], c='g', ls='-',label='non-overlap')
    axis.plot(x, counter[:,2], c='b', ls='-',label='empty')
    axis.axhline(threshold, ls='--', lw = 0.5)
    axis.legend(loc="center right")
    # axis.ylabel('overlap fraction')
    fig.tight_layout()
    plt.savefig(fileName)
    fig.clf()
    plt.close(fig)
    plt.clf()
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
