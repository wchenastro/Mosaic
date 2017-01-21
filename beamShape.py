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

def fringePlot(beamCoordinates, weights, baselines, boresight, beamAperture, interpolating=True, outputFormat='png', allFringe=False):
    beamPattern = primaryBeamPattern(beamAperture, boresight)
    # fringeFile = open('fringes', 'w')
    # factorFile = open('factor', 'w')
    # beamFile = open('beam', 'w')
    fringeSum = [0]*len(beamCoordinates)
    for weightRow in weights:
        rowSum = []
        offset = np.angle(weightRow[0])
        for coord, weight, fringe in zip(beamCoordinates, weightRow, fringeSum):
            beamFactor = beamPattern.beamFactor(coord)
            # beamFactor = 1
            # factorFile.write(str(beamFactor) + '\n')
            amplitude = beamFactor*np.cos(np.angle(weight) - offset)
            rowSum.append(fringe + amplitude)
            # fringeFile.write(' '.join([str(coord[0]), str(coord[1]), str(amplitude)]) + '\n')
        fringeSum = rowSum[:]
        # fringeFile.write('\n')
    # fringeFile.close()
    # factorFile.close()
    baselinesString = ''
    for baseline in baselines:
        baselineINT = '{:7.0f}'.format(baseline[0]) + '{:7.0f}'.format(baseline[1]) + '{:7.0f}'.format(baseline[2])
        baselinesString += baselineINT + " " + '{:7.0f}'.format(np.linalg.norm(baseline)) + ' meters '
    stringSegmentLength = len(baselinesString)/len(weights)
    # plotFringesViaGnuplot('fringes', len(weights))
    if(allFringe == True):
        plotFringesContour('fringes', len(weights), baselinesString, stringSegmentLength, outputFormat)

    # beamSynthesized = []
    # for coord, weight in zip(beamCoordinates, fringeSum):
        # beamSynthesized.append([coord[0], coord[1], weight])
        # beamFile.write(' '.join([str(coord[0]), str(coord[1]), str(weight)]) + '\n')

    # beamFile.close()

    angle = 0;
    axis1 = 0;
    axis2 = 0;
    # if len(weights) >= 6:
        # try:
            # points = readMarginalPoints(len(weights)*0.4)
            # angle, axis1, axis2 = fitContour(points)
        # except Exception as e:
            # print(e)

    # plotBeamViaGnuplot('beam')
    blocks = len(weights)
    if interpolating == True:
        # plotBeamContour('beam', len(weights), outputFormat)
        plotBeamContour2(beamCoordinates[:,0], beamCoordinates[:,1], np.abs(fringeSum), blocks)
    else:
        # plotBeamContourDot('beam', len(weights), outputFormat, np.rad2deg(angle), axis1, axis2)
        plotBeamScatter(beamCoordinates[:,0], beamCoordinates[:,1], np.abs(fringeSum), blocks)

    # with open('beam', 'w') as fringeFile:
        # fringeSum = [0]*len(beamCoordinates)
        # for weightRow in weights:
            # rowSum = []
            # offset = np.angle(weightRow[0])
            # for weight, fring, coord in zip(weightRow, fringeSum, beamCoordinates):
                # beamFactor = beamPattern.beamFactor(coord)
                # rowSum.append(fring + np.cos(np.angle(weight) - offset)*beamFactor)
            # fringeSum = rowSum[:]

        # for coord, weight in zip(beamCoordinates, fringeSum):
            # fringeFile.write(' '.join([str(coord[0]), str(coord[1]), str(weight)]) + '\n')

    # plotBeamViaGnuplot('beam')

    return fringeSum

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
    ellipse = fitellipse.fitellipse(xy)
    center = ellipse[0]
    axis1 = ellipse[1]
    axis2 = ellipse[2]
    angle = ellipse[3]
    print(np.rad2deg(angle), axis1, axis2)

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

def plotBeamContour2(x, y, z, maxValue):
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
    plt.savefig('contour.png', dpi=thisDpi)
    plt.close()

def plotBeamScatter(x, y, z, maxValue):
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
    plt.savefig('contour.png', dpi=thisDpi)
    plt.close()


def plotPackedBeam(coordinates, angle, axis1, axis2, beamRadius):
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
    plt.savefig('pack.png', dpi=thisDpi)
    plt.close()

def plotBeamFit(coordinates, center, angle, axis1, axis2):
    xMin = coordinates[:,0].min()
    xMax = coordinates[:,0].max()
    yMin = coordinates[:,1].min()
    yMax = coordinates[:,1].max()
    thisDpi = 96
    matplotlib.rcParams.update({'font.size': 8})
    plt.clf()
    fig = plt.figure(figsize=(400./thisDpi, 300./thisDpi), dpi=thisDpi)
    plt.subplots_adjust(right=0.75)
    ax = fig.add_subplot(111, aspect='equal')
    ellipse = Ellipse(xy=center, width=2*axis1, height=2*axis2, angle=angle)
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
    plt.savefig('fit.png', transparent=True, dpi=thisDpi)
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
