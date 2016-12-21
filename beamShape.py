#!/usr/bin/env python

import numpy as np
import os


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

def fringePlot(beamCoordinates, weights, baselines, boresight, beamAperture, outputFormat='png', allFringe=False):
    beamPattern = primaryBeamPattern(beamAperture, boresight)
    fringeFile = open('fringes', 'w')
    beamFile = open('beam', 'w')
    fringeSum = [0]*len(beamCoordinates)
    for weightRow in weights:
        rowSum = []
        offset = np.angle(weightRow[0])
        for coord, weight, fringe in zip(beamCoordinates, weightRow, fringeSum):
            beamFactor = beamPattern.beamFactor(coord)
            # beamFactor = 1
            amplitude = beamFactor*np.cos(np.angle(weight) - offset)
            rowSum.append(fringe + amplitude)
            fringeFile.write(' '.join([str(coord[0]), str(coord[1]), str(amplitude)]) + '\n')
        fringeSum = rowSum[:]
        fringeFile.write('\n')
    fringeFile.close()
    baselinesString = ''
    for baseline in baselines:
        baselineINT = '{:7.0f}'.format(baseline[0]) + '{:7.0f}'.format(baseline[1]) + '{:7.0f}'.format(baseline[2])
        baselinesString += baselineINT + " " + '{:7.0f}'.format(np.linalg.norm(baseline)) + ' meters '
    stringSegmentLength = len(baselinesString)/len(weights)
    # plotFringesViaGnuplot('fringes', len(weights))
    if(allFringe == True):
        plotFringesContour('fringes', len(weights), baselinesString, stringSegmentLength, outputFormat)

    for coord, weight in zip(beamCoordinates, fringeSum):
        beamFile.write(' '.join([str(coord[0]), str(coord[1]), str(weight)]) + '\n')

    beamFile.close()

    # plotBeamViaGnuplot('beam')
    plotBeamContour('beam', len(weights), outputFormat)

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



class primaryBeamPattern:

    def __init__(self, aperture, center):
        self.aperture = aperture
        self.center = center
        self.normalization = self.normal(0, aperture, 0)

    def normal(self, mu, sigma, x):

        return np.exp(-(x-mu)**2/(2*sigma**2))/(sigma*np.sqrt(2*np.pi))

    def beamFactor(self, offset):

        distance = np.sqrt((offset[0] - self.center[0])**2
                + (offset[1] - self.center[1])**2)

        normalFactor = self.normal(0, self.aperture, distance)

        return normalFactor/self.normalization

def plotFringesViaGnuplot(dataFile, blockNumber):

    plotScript="\"\
                set encoding utf8;\
                set terminal pdf font ',15';\
                set output 'fringes.pdf';\
                set size ratio -1;\
                set xlabel 'RA';\
                set ylabel 'DEC';\
                set cblabel 'phase' font ',19';\
                set xtics rotate by -60;\
                dataFile='%s';\
                blknum=%s - 1;\
                do for [i=0:blknum] { plot dataFile every :::i::i u 1:2:3 with points  pt 7 ps 1 lt palette notitle};\
                set output;\
                \""

    os.system("gnuplot -e " + plotScript % (dataFile, blockNumber));

def plotBeamViaGnuplot(dataFile):

    plotScript="\"\
                set encoding utf8;\
                set terminal pdf font ',15';\
                set output 'beam.pdf';\
                set size ratio -1;\
                set xlabel 'RA';\
                set ylabel 'DEC';\
                set xtics rotate by -60;\
                dataFile='%s';\
                plot dataFile u 1:2:3 with points pt 7 ps 1 lt palette notitle;\
                set output;\
                \""

    os.system("gnuplot -e " + plotScript % (dataFile));

def plotBeamContour(dataFile, blockNumber, outputFormat='png'):

    plotScript="\"\
                set encoding utf8;\
                outputFormat='%s';\
                set terminal outputFormat font ',17';\
                set output 'contour.'.outputFormat;\
                set print 'debuglog';\
                set size ratio -1;\
                set xlabel 'RA' font ',11';\
                set ylabel 'DEC' font ',11';\
                dataFile='%s';\
                stats dataFile using 1 name 'x' nooutput;\
                stats dataFile using 2 name 'y' nooutput;\
                tics = 0.003;\
                set xtics x_min-tics, tics, x_max+tics rotate by -60;\
                set ytics y_min-tics, tics, y_max+tics;\
                set dgrid3d 70,70 gauss 0.001;\
                unset surface;\
                set contour base;\
                set cbrange [0:%s];\
                set cntrlabel onecolor font ',4' start 15 interval -1;\
                set cntrparam levels auto;\
                set view map;\
                set pm3d;\
                set style textbox margins  0.3,  0.3 noborder;\
                set style data lines;\
                set multiplot;\
                splot dataFile using 1:2:(abs(\$3)) notitle with line lc rgb '#702e3f54';\
                unset pm3d;\
                unset xtics;\
                splot dataFile  using 1:2:(abs(\$3)) notitle with labels center boxed offset 0,-0.25;\
                unset multiplot;\
                set output;\
                \""

    os.system("gnuplot -e " + plotScript % (outputFormat, dataFile, blockNumber));

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
