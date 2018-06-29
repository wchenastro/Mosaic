import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

def plotBeamContour(array, center, sideLength, fileName='contour.png', interpolation = True):
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
    fig.gca().set_aspect('equal', 'datalim')
    # axes = plt.gca()
    # axes.set_xlim([x.min(),x.max()])
    # axes.set_ylim([y.min(),y.max()])
    plt.savefig(fileName, dpi=thisDpi)
    plt.close()

def plotOverlap(overlapTable, fileName='overlap.png'):
    pointSize = 10
    thisDpi = 96
    matplotlib.rcParams.update({'font.size': 8})
    plt.clf()
    fig = plt.figure(figsize=(400./thisDpi, 300./thisDpi), dpi=thisDpi)
    plt.imshow(overlapTable, cmap=plt.cm.jet, vmin=0)
    # plt.contour(overlapTable, cmap=plt.cm.jet, vmin=0, vmax=1)
    axis = fig.gca()
    axis.set_aspect('equal', 'datalim')
    plt.colorbar()
    # axis.axis('off')
    plt.savefig(fileName, dpi=thisDpi)
    plt.close()

def plotPackedBeam(coordinates, angle, axis1, axis2, beamRadius, fileName='pack.png', scope=1.):
    thisDpi = 96
    matplotlib.rcParams.update({'font.size': 8})
    plt.clf()
    fig = plt.figure(figsize=(400./thisDpi, 300./thisDpi), dpi=thisDpi)
    # plt.axes().set_aspect('equal', 'datalim')
    ax = fig.add_subplot(111, aspect='equal')
    for coord in coordinates:
        ellipse = Ellipse(xy=coord, width=2*axis1, height=2*axis2, angle=angle)
        ellipse.fill = False
        ax.add_artist(ellipse)
    circle = Ellipse(xy=(0,0), width=2*beamRadius, height=2*beamRadius, angle=0)
    circle.fill = False
    ax.add_artist(circle)
    ax.set_xlim(-beamRadius*1.3*scope, beamRadius*1.3*scope)
    ax.set_ylim(-beamRadius*1.3*scope, beamRadius*1.3*scope)
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
    plt.subplots_adjust(right=0.75)
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

