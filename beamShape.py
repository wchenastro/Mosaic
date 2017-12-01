#!/usr/bin/env python

import numpy as np
import scipy.interpolate
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib
import os

import fitellipse


def createBaselines(antCoordinatesENU):
    baselines = []
    index = 1
    for antenna1 in antCoordinatesENU:
        for antenna2 in antCoordinatesENU[index:]:
            baselines.append([antenna1[0]-antenna2[0],
                    antenna1[1]-antenna2[1],
                    antenna1[2]-antenna2[2]])
        index += 1

    return baselines

def fringePlot(beamCoordinates, weights, baselines, boresight, beamAperture, interpolating=True, outputFormat='png', allFringe=False, fileName='contour.png'):
    # beamPattern = primaryBeamPattern(beamAperture, boresight)
    # fringesOfAllBaselines = []
    phaseShifters = np.conj(weights[:,0][:,None])
    fringes = weights * phaseShifters
    sumOfFringes = np.sum(fringes.real, axis=0)
    # for weightOfEachBaseline in weights:
        # fringesOfEachBaseline = []
        # phaseShifter = np.conj(weightOfEachBaseline[0])
        # for beamCoord, weightOfBeam in zip(beamCoordinates, weightOfEachBaseline):
            # beamFactor = beamPattern.beamFactor(beamCoord)
            # beamFactor = 1
            # amplitude = beamFactor*(weightOfBeam*phaseShifter).real
            # fringesOfEachBaseline.append(amplitude)
        # fringesOfAllBaselines.append(fringesOfEachBaseline)
    # sumOfFringes = np.sum(fringesOfAllBaselines, axis=0)

    blocks = len(weights)
    if interpolating == True:
        plotBeamContour2(beamCoordinates[:,0], beamCoordinates[:,1], np.abs(sumOfFringes), blocks, fileName)
    else:
        plotBeamScatter(beamCoordinates[:,0], beamCoordinates[:,1], np.abs(sumOfFringes), blocks, fileName)

    return sumOfFringes

def readMarginalPoints(coordinates, amplitude, margin):
    offset = margin*0.1
    marginPoints = []

    for point, amp in zip(coordinates, amplitude):
        if abs(np.abs(amp)  - margin) < offset:
            marginPoints.append(point)

    return np.array(marginPoints)


def readMarginalPointsFromFile(margin):

    offset = margin*0.1
    marginPoints = []

    with open('beam', 'r') as coordFile:
        allPoints = np.loadtxt(coordFile)

    with open('marginal', 'w') as marginalFile:
        for point in allPoints:
            if abs(point[2]  - margin) < offset:
                marginPoints.append(point)
                marginalFile.write(' '.join([str(point[0]), str(point[1]), str(point[2])]) + '\n')

    return np.array(marginPoints)

def fitContour(points):
    xy = [points[:,0], points[:,1]]
    try:
        ellipse = fitellipse.fitellipse(xy)
    except RuntimeError as e:
        print str(e)
        return [], None, None, None
    center = ellipse[0]
    axis1 = ellipse[1]
    axis2 = ellipse[2]
    angle = ellipse[3]
    # print(np.rad2deg(angle), axis1, axis2)

    return center, angle, axis1, axis2


def fitEllipseBeam(coordinates, amplitude, marginal):
    points = readMarginalPoints(coordinates, amplitude, marginal)
    center, angle, axis1, axis2 = fitContour(points)

    return center, angle, axis1, axis2



class primaryBeamPattern:

    def __init__(self, aperture, center):
        self.aperture = aperture*0.32
        self.center = center
        self.normalization = self.normal(0, self.aperture, 0)

    def normal(self, mu, sigma, x):

        return np.exp(-(x-mu)**2/(2*(sigma**2)))/(sigma*np.sqrt(2*np.pi))

    def beamFactor(self, offset):

        distance = np.sqrt((offset[0] - self.center[0])**2
                + (offset[1] - self.center[1])**2)

        normalFactor = self.normal(0, self.aperture, distance)

        return normalFactor/self.normalization

def plotBeamContour2(x, y, z, maxValue, fileName='contour.png'):
    thisDpi = 96.
    matplotlib.rcParams.update({'font.size': 8})
    xi = np.linspace(x.min(), x.max(), len(x))
    yi = np.linspace(y.min(), y.max(), len(y))
    zi = scipy.interpolate.griddata((x, y), z, (xi[None,:], yi[:,None]), method='cubic')
    plt.clf()
    plt.figure(figsize=(400./thisDpi, 300./thisDpi), dpi=thisDpi)
    plt.contourf(xi,yi,zi,15,cmap=plt.cm.jet, vmin=0, vmax=maxValue)
    plt.colorbar()
    plt.axes().set_aspect('equal', 'datalim')
    axes = plt.gca()
    axes.set_xlim([x.min(),x.max()])
    axes.set_ylim([y.min(),y.max()])
    plt.savefig(fileName, dpi=thisDpi)
    plt.close()

def plotBeamContour3(array, center, sideLength, fileName='contour.png', interpolation = True):
    thisDpi = 96.
    matplotlib.rcParams.update({'font.size': 8})
    plt.figure(figsize=(400./thisDpi, 300./thisDpi), dpi=thisDpi)
    halfSideLength = sideLength/2.0
    xStart = np.rad2deg(center[0] - halfSideLength)
    xEnd = np.rad2deg(center[0] + halfSideLength)
    yStart = np.rad2deg(center[1] - halfSideLength)
    yEnd = np.rad2deg(center[1] + halfSideLength)
    plotRange = [xStart, xEnd, yStart, yEnd]
    interpolateOption = 'bicubic' if interpolation == True else 'nearest' 
    plt.imshow(array,cmap=plt.cm.jet,  interpolation=interpolateOption, extent=plotRange)
    plt.colorbar()
    plt.axes().set_aspect('equal', 'datalim')
    axes = plt.gca()
    # axes.set_xlim([x.min(),x.max()])
    # axes.set_ylim([y.min(),y.max()])
    plt.savefig(fileName, dpi=thisDpi)
    plt.close()


def plotBeamScatter(x, y, z, maxValue, fileName='contour.png'):
    pointSize = 10
    thisDpi = 96
    matplotlib.rcParams.update({'font.size': 8})
    plt.clf()
    plt.figure(figsize=(400./thisDpi, 300./thisDpi), dpi=thisDpi)
    plt.scatter(x,y,c=z,s=pointSize, cmap=plt.cm.jet, vmin=0, vmax=maxValue, edgecolor='none')
    plt.colorbar()
    plt.axes().set_aspect('equal', 'datalim')
    axes = plt.gca()
    axes.set_xlim([x.min(),x.max()])
    axes.set_ylim([y.min(),y.max()])
    plt.savefig(fileName, dpi=thisDpi)
    plt.close()

def plotAntennas(lons, lats, fileName='antennas.png'):
    pointSize = 10
    thisDpi = 96
    matplotlib.rcParams.update({'font.size': 8})
    plt.clf()
    plt.figure(figsize=(400./thisDpi, 300./thisDpi), dpi=thisDpi)
    plt.scatter(lons,lats)
    plt.axes().set_aspect('equal', 'datalim')
    axes = plt.gca()
    axes.set_xlim([lons.min(), lons.max()])
    axes.set_ylim([lats.min(), lats.max()])
    plt.savefig(fileName, dpi=thisDpi)
    plt.close()

def plotAntsNSource(lons, lats, center, alt, azi, fileName='AntsNSource.png'):

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

    # radius = 0.04

    antX = lons - center[0]
    antY = lats - center[1]

    antXMax = np.abs(antX).max()
    antYMax = np.abs(antY).max()

    radius = antXMax*1.2 if antXMax > antYMax else antYMax*1.2

    pointSize = 3
    thisDpi = 96
    matplotlib.rcParams.update({'font.size': 8})
    plt.clf()

    plt.figure(figsize=(400./thisDpi, 300./thisDpi), dpi=thisDpi)
    plt.scatter(antX,antY, s=pointSize)
    lineStyle = '--' if alt < 0 else '-'
    ax = plt.gcf().gca()
    ax.add_artist(plt.Circle((0,0), radius, fill=False))
    ax.add_artist(plt.Circle((0,0), (90-np.abs(np.rad2deg(alt)))/90.*radius, fill=False, ls=lineStyle))
    x, y = angleToCartesian(np.rad2deg(azi), radius)
    print 'azi, alt: ', np.rad2deg(azi), np.rad2deg(alt)
    plt.plot([0,x], [0, y], ls=lineStyle)

    plt.axes().set_aspect('equal', 'datalim')
    ax.set_xlim([-radius*1.3, +radius*1.3])
    ax.set_ylim([-radius*1.3, +radius*1.3])
    plt.savefig(fileName, dpi=thisDpi)
    plt.close()

def plotPackedBeam(coordinates, angle, axis1, axis2, beamRadius, fileName='pack.png'):
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
    ax.set_xlim(-beamRadius*1.3, beamRadius*1.3)
    ax.set_ylim(-beamRadius*1.3, beamRadius*1.3)
    plt.savefig(fileName, dpi=thisDpi)
    plt.close()

def plotBeamFit(sideLength, center, angle, axis1, axis2, fileName='fit.png'):
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
    step=sideLength/20.0
    trueCenter = [center[0] + step/2., center[1] - step/2.]
    ellipse = Ellipse(xy=trueCenter, width=2*axis1, height=2*axis2, angle=angle)
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


def plotBeamContourDot(dataFile, blockNumber, outputFormat, angle, axis1, axis2):

    plotScript="\"\
                set encoding utf8;\
                outputFormat='%s';\
                set terminal outputFormat font ',7' size 430, 330;\
                set output 'contour.'.outputFormat;\
                set size ratio -1;\
                set bmargin 5;\
                set lmargin 3;\
                set rmargin 0;\
                set tmargin 2;\
                set xlabel 'RA' font ',7' offset 0,1;\
                set ylabel 'DEC' font ',7' offset 2;\
                dataFile='%s';\
                stats dataFile every ::::0 using 1 name 'x' nooutput;\
                stats dataFile every ::::0 using 2 name 'y' nooutput;\
                boresightX = x_min;\
                boresightY = y_min;\
                set xtics rotate by -60;\
                set cbrange [0:%s];\
                axis1=%s;\
                axis2=%s;\
                if (axis1 != 0) {\
                set object 1 ellipse at 0,0 size 2*axis1,2*axis2 angle %s;}\
                plot dataFile using ((\$1)-boresightX):((\$2)-boresightY):(abs(\$3)) notitle with points palette pt 7 ps 0.7,\
                'marginal' using ((\$1)-boresightX):((\$2)-boresightY) notitle with points pt 7 ps 0.7;\
                set output;\
                \""

    os.system("gnuplot -e " + plotScript % (outputFormat, dataFile, blockNumber, axis1, axis2, angle));


def plotBeamContour(dataFile, blockNumber, outputFormat='png'):

    plotScript="\"\
                set encoding utf8;\
                outputFormat='%s';\
                set terminal outputFormat font ',7'  size 450, 350;\
                set output 'contour.'.outputFormat;\
                set print 'debuglog';\
                set size ratio -1;\
                set bmargin 3;\
                set lmargin 1;\
                set rmargin 0;\
                set tmargin 2;\
                set xlabel 'RA' font ',7';\
                set ylabel 'DEC' font ',7' offset 1;\
                dataFile='%s';\
                stats dataFile every ::::0 using 1 name 'x' nooutput;\
                stats dataFile every ::::0 using 2 name 'y' nooutput;\
                boresightX = x_min;\
                boresightY = y_min;\
                set dgrid3d 70,70 gauss 0.001;\
                unset surface;\
                set contour base;\
                set cbrange [0:%s];\
                set cntrlabel onecolor font ',4' start 15 interval -1;\
                set cntrparam levels auto;\
                set view map;\
                set pm3d;\
                set style textbox margins  0.3,  0.3 noborder;\
                set xtics rotate by -60;\
                set multiplot;\
                splot dataFile using ((\$1)-boresightX):((\$2)-boresightY):(abs(\$3)) notitle with line lc rgb '#702e3f54';\
                unset pm3d;\
                unset xtics;\
                unset ytics;\
                unset xlabel;\
                unset ylabel;\
                splot dataFile  using ((\$1)-boresightX):((\$2)-boresightY):(abs(\$3)) notitle with labels center boxed offset 0,-0.25;\
                unset multiplot;\
                set output;\
                \""

    os.system("gnuplot -e " + plotScript % (outputFormat, dataFile, blockNumber));

def plotPackedBeam2(dataFile, angle, axis1, axis2, beamRadius, outputFormat, ):

    plotScript="\"\
                set encoding utf8;\
                outputFormat='%s';\
                set terminal outputFormat.'cairo' font ',7' size 430, 330;\
                set output 'pack.'.outputFormat;\
                set size ratio -1;\
                set bmargin 5;\
                set lmargin 0;\
                set rmargin 0;\
                set tmargin 3;\
                set xlabel 'RA' font ',7' ;\
                set ylabel 'DEC' font ',7' ;\
                dataFile='%s';\
                beamRadius=%s;\
                axis1=%s;\
                axis2=%s;\
                set xtics rotate by -60;\
                set object 1 circle at 0,0 size beamRadius;\
                plot dataFile u 1:2:(2*axis1):(2*axis2):(%s) with ellipse notitle;\
                set output;\
                \""

    os.system("gnuplot -e " + plotScript % (outputFormat, dataFile, beamRadius, axis1, axis2, angle));



def plotFringesContour(dataFile, blockNumber, text, segmentLength, outputFormat='png'):

    plotScript="\"\
                set encoding utf8;\
                outputFormat='%s';\
                set terminal outputFormat font ',5';\
                set output 'fringeContour.'.outputFormat;\
                set size ratio -1;\
                set xlabel 'RA' font ',5';\
                set ylabel 'DEC' font ',5' offset 5;\
                dataFile='%s';\
                blockNumber=%s - 1;\
                titlesString ='%s';\
                titleLength = %s;\
                stats dataFile using 1 name 'x' nooutput;\
                stats dataFile using 2 name 'y' nooutput;\
                tics = 0.005;\
                set xtics x_min-tics, tics, x_max+tics;\
                set ytics y_min-tics, tics, y_max+tics rotate by -90;\
                set dgrid3d 70,70 gauss 0.00005;\
                unset surface;\
                set cntrparam levels auto;\
                set view map;\
                set pm3d;\
                set style textbox margins  0.3,  0.3 noborder;\
                set style data lines;\
                set print 'debuglog';\
                set multiplot layout 3,5;\
                do for [i=0:blockNumber] {\
                    fringeTitle = titlesString[i*titleLength+1:(i+1)*titleLength+1];\
                    set title fringeTitle;\
                    splot dataFile every :::i::i using 1:2:3 with line notitle;\
                }\
                unset multiplot;\
                set output;\
                \""

    os.system("gnuplot -e " + plotScript % (outputFormat, dataFile, blockNumber, text, segmentLength));
