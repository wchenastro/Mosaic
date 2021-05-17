#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os, pickle
import numpy as np

from PyQt4.QtGui import *
from PyQt4.QtCore import *

from mosaic.interferometer import InterferometryObservation
from mosaic.plot import plotPackedBeam, plotBeamFit, plotBeamWithFit
from mosaic.tile import ellipseGrid, ellipseCompact
from mosaic.beamshape import calculateBeamOverlaps
from mosaic.coordinate import Antenna, Boresight, convertBoresightToHour
from mosaic.coordinate import convertBoresightToDegree, convert_pixel_coordinate_to_equatorial

import argparse
import logging, tempfile


class Cartesian(QWidget):


    def __init__(self, parent=None):
        QWidget.__init__(self, parent)

        self.painter = None
        self.dots = []
        self.highLightDots = []
        self.azimuth = 0
        self.elevation = 0
        self.width = 300.
        self.height = 300.
        self.halfWidth = 150.
        self.halfHeight = 150.
        self.xStart = 0
        self.xEnd = 0
        self.yStart = 0
        self.yEnd = 0


    def pixelCoordinateConv(self, value, direction):
        xRange = self.xEnd-self.xStart
        yRange = self.yEnd-self.yStart
        if direction == 'toCoord':
            x = (value[0]/self.width)*xRange+self.xStart
            y = ((self.height - value[1])/self.height)*yRange + self.yStart
        elif direction == 'toPixel':
            x = int(round((value[0] - self.xStart)/xRange*self.width))
            y = self.height - int(round((value[1] - self.yStart)/yRange*self.height))

        return [x,y]



    def paintEvent(self, event):
        self.setMouseTracking(True)

        painter = QPainter()
        self.painter = painter
        dashPen = QPen(Qt.DashLine)
        solidPen = QPen(Qt.SolidLine)
        painter.begin(self)
        painter.drawRect(0, 0, self.width, self.height)
        '''centered horizontal line'''
        painter.drawLine(0, self.halfHeight, self.width, self.halfHeight)
        '''centered vertical line'''
        painter.drawLine(self.halfWidth, self.height, self.halfWidth, 0)
        '''horizon'''
        painter.drawEllipse(QPoint(self.halfWidth, self.halfHeight),
                self.halfWidth, self.halfHeight)
        '''latitude'''
        if self.elevation < 0:
            painter.setPen(dashPen)
        painter.drawEllipse(QPoint(self.halfWidth, self.halfHeight),
                (90.-np.abs(self.elevation))/90.*self.halfWidth,
                (90.-np.abs(self.elevation))/90.*self.halfHeight)
        angleX, angleY = self.angleToCartesian(self.azimuth, self.halfWidth)
        painter.drawLine(self.halfWidth, self.halfHeight, angleX, angleY)
        painter.setPen(solidPen)

        for dot in self.dots:
            painter.drawEllipse(dot[0], dot[1],3,3)
            # painter.drawPoint(dot[0], dot[1])

        painter.setPen(Qt.blue)
        for dot in self.highLightDots:
            painter.drawEllipse(dot[0], dot[1],5,5)
            # painter.drawPoint(dot[0], dot[1])
        painter.setPen(Qt.black)
        self.highLightDots[:] = []

        painter.end()

    def sizeHint(self):
        return QSize(301, 301)

    def mouseMoveEvent(self, event):
        longitude, latitude = self.pixelCoordinateConv((event.x(), event.y()), 'toCoord')
        longitudeCoord.setPlaceholderText('{:6.4f}'.format(longitude))
        latitudeCoord.setPlaceholderText('{:6.4f}'.format(latitude))

    def addDots(self, dots):
        for latitude, longitude in dots:
            dot = self.pixelCoordinateConv((longitude, latitude), 'toPixel')
            self.dots.append(dot)
        self.update()

    def addHighLightDots(self, dots):
        for latitude, longitude in dots:
            dot = self.pixelCoordinateConv((longitude, latitude), 'toPixel')
            self.highLightDots.append(dot)
        self.update()

    def mouseReleaseEvent(self, event):
        height = 1035
        self.dots.append([event.x(), event.y()])
        longitude, latitude = self.pixelCoordinateConv((event.x(), event.y()), 'toCoord')
        position = coordinateList.rowCount()
        addRowToCoordinateList('*', position, latitude, longitude, height)
        self.update()
        updateContour()
        # updateBaselineList(observation.getBaselines())

    def angleToCartesian(self, angleDeg, radius):
        angle = np.deg2rad(angleDeg)
        if angle < np.pi/2:
            x = np.sin(angle)*radius + radius
            y = radius - np.cos(angle)*radius
        elif angle >= np.pi/2 and angle < np.pi:
            x = radius + np.sin(angle)*radius
            y = radius - np.cos(angle)*radius
        elif angle >= np.pi and angle < np.pi*1.5:
            x = radius + np.sin(angle)*radius
            y = -np.cos(angle)*radius + radius
        elif angle >= np.pi*1.5 and angle < np.pi*2:
            x = radius + np.sin(angle)*radius
            y = radius - np.cos(angle)*radius

        return x, y

    def clearDots(self):
        self.dots = []
        self.update()

    def removeDots(self, dotsToRemove):
        for latitude, longitude in dotsToRemove:
            dot = self.pixelCoordinateConv((longitude, latitude), 'toPixel')
            self.dots.remove(dot)
        self.update()

    def setAzAlt(self, horizontalCoord):
        self.azimuth = horizontalCoord[0]
        self.elevation = horizontalCoord[1]
        self.update()

    def setCenter(self, center, radius, swapXY=False):
        if swapXY == True:
            center = center[1], center[0]

        self.xStart = center[0] - radius
        self.xEnd = center[0] + radius
        self.yStart = center[1] - radius
        self.yEnd = center[1] + radius

        # print (center, self.xStart, self.xEnd, self.yStart, self.yEnd)

class miniCartesian(Cartesian):

    def paintEvent(self, event):
        self.setMouseTracking(True)

        painter = QPainter()
        self.painter = painter
        dashPen = QPen(Qt.DashLine)
        solidPen = QPen(Qt.SolidLine)
        painter.begin(self)
        painter.drawRect(0, 0, self.width, self.height)
        '''centered horizontal line'''
        painter.drawLine(0, self.halfHeight, self.width, self.halfHeight)
        '''centered vertical line'''
        painter.drawLine(self.halfWidth, self.height, self.halfWidth, 0)

        for dot in self.dots:
            painter.drawEllipse(dot[0], dot[1],1,1)

        # for dot in self.highLightDots:
            # painter.drawEllipse(dot[0], dot[1],5,5)
        # self.highLightDots[:] = []

        painter.end()

    def mouseReleaseEvent(self, event):
        pass

    def mouseMoveEvent(self, event):
        pass

def updateHorizontal(horizontalCoord):
    if horizontalCoord == []: return
    azimuth = np.rad2deg(horizontalCoord[0])
    elevation = np.rad2deg(horizontalCoord[1])
    axis.setAzAlt([azimuth, elevation])

    azimuthCoord.blockSignals(True)
    elevationCoord.blockSignals(True)
    azimuthCoord.setText('{:6.4f}'.format(azimuth))
    elevationCoord.setText('{:6.4f}'.format(elevation))
    azimuthCoord.blockSignals(False)
    elevationCoord.blockSignals(False)

def updateBoreSight(boreSightCoord):
    if boreSightCoord == []: return
    RA, DEC = convertBoresightToHour(boreSightCoord);
    RACoord.blockSignals(True)
    DECCoord.blockSignals(True)
    # RACoord.setText('{:6.5f}'.format(RA))
    # DECCoord.setText('{:6.5f}'.format(DEC))
    RACoord.setText(RA)
    DECCoord.setText(DEC)
    RACoord.blockSignals(False)
    DECCoord.blockSignals(False)


def onBoreSightUpdated():
    beamBoreSight = convertBoresightToDegree((str(RACoord.text()), str(DECCoord.text())))
    # beamBoreSight = (float(RACoord.text()), float(DECCoord.text()))
    preBoresight = observation.getBoreSight().equatorial
    if ((abs(beamBoreSight[0] - preBoresight[0]) < 1e-4) and (abs(beamBoreSight[1] - preBoresight[1]) < 1e-4)):
        return
    print('boresight edited')
    observation.setBoreSight(beamBoreSight, Boresight.EquatorialFrame)
    updateContour()
    updateHorizontal(observation.getBoreSight().horizontal)

def onHorizontalUpdated():
    horizontalCoord = (float(azimuthCoord.text()), float(elevationCoord.text()))
    preHorizontal = np.rad2deg(observation.getBoreSight().horizontal)
    if ((abs(horizontalCoord[0] - preHorizontal[0]) < 1e-4) and (abs(horizontalCoord[1] - preHorizontal[1]) < 1e-4)):
        return
    print('horizontal edited')
    observation.setBoreSight(horizontalCoord, Boresight.HorizontalFrame)
    updateContour()
    axis.setAzAlt(horizontalCoord)
    updateBoreSight(observation.getBoreSight().equatorial)


def onBeamSizeChanged(newValue):
    autoZoomCheckbox.setCheckState(Qt.Unchecked)
    isSet = observation.setBeamSizeFactor(newValue)
    if isSet:
        updateContour()
    else:
        beamSizeEdit.setValue(observation.getBeamSizeFactor())

def onRotationChanged():
    resetPackState()


def onBeamNumberChanged():
    isSet = observation.setBeamNumber(float(beamNumberEdit.text()))
    if isSet:
        updateContour()
    else:
        beamNumberEdit.setText(str(int(observation.getBeamNumber())))

def onInterpolateOptionChanged(state):
    if state == Qt.Checked:
        observation.setInterpolating(True)
    else:
        observation.setInterpolating(False)
    updateContour()

def onAutoZoomOptionChanged(state):
    if state == Qt.Checked:
        observation.setAutoZoom(True)
    else:
        observation.setAutoZoom(False)
    updateContour()

def onDateTimeChanged(dateTime):
    observation.setBoreSight(frame = Boresight.EquatorialFrame)
    observation.setObserveTime(dateTime.toPyDateTime())
    updateContour()
    updateHorizontal(observation.getBoreSight().horizontal)

def addRowToCoordinateList(name, position, latitude, longitude, height, hidden=''):
    coordinateList.insertRow(position)
    coordinateList.setVerticalHeaderItem(position, QTableWidgetItem(str(name)))
    '''
    coordinateList.setItem(position, 0, QTableWidgetItem('{:6.4f}'.format(latitude)))
    coordinateList.setItem(position, 1, QTableWidgetItem('{:6.4f}'.format(longitude)))
    coordinateList.setItem(position, 2, QTableWidgetItem('{:6.4f}'.format(height)))
    '''
    coordinateList.setItem(position, 0, QTableWidgetItem(str(latitude)))
    coordinateList.setItem(position, 1, QTableWidgetItem(str(longitude)))
    coordinateList.setItem(position, 2, QTableWidgetItem(str(height)))

    coordinateList.setItem(position, 3, QTableWidgetItem(hidden))
    coordinateList.setItem(position, 4, QTableWidgetItem('-'))

def resetPackState():
    onClickedPackButton2.newData = True
    onClickedPackButton2.state = 0
    onClickedPackButton2.fittedImage = None

def collectCoordinates():
    rowCount = coordinateList.rowCount()
    coordinates = []
    hidden = []
    for row in range(rowCount):
        #item = coordinateList.item(row, 3)
        #if item != None and item.text() == 'hidden': continue
        longitude = float(str(coordinateList.item(row, 1).text()))
        latitude = float(str(coordinateList.item(row, 0).text()))
        height = float(str(coordinateList.item(row, 2).text()))
        hiddenFlag = str(coordinateList.item(row, 3).text())
        if hiddenFlag == 'hidden':
            hidden.append(row)
        coordinates.append([latitude, longitude, height])

    return coordinates, hidden

def updateContour():
    # print('updateContour')

    # if not hasattr(updateContour, "lastObservetime"):
        # updateContour.lastObservetime = None

    # observeTime = dateTimeEdit.dateTime().toPyDateTime()
    # if(updateContour.lastObservetime != observeTime):
        # observation.setObserveTime(observeTime)
        # updateContour.lastObservetime = observeTime

    #FitsButton.disconnect()
    modifiers = QApplication.keyboardModifiers()
    if modifiers == Qt.ShiftModifier:
        observation.setAutoZoom(False)

    rowCount = coordinateList.rowCount()
    if rowCount < 2:
        label.setPixmap(blankImage)
        return
    coordinates = []
    for row in range(rowCount):
        item = coordinateList.item(row, 3)
        if item is not None and item.text() == 'hidden': continue
        longitude = float(str(coordinateList.item(row, 1).text()))
        latitude = float(str(coordinateList.item(row, 0).text()))
        height = float(str(coordinateList.item(row, 2).text()))
        coordinates.append([latitude, longitude, height])

    antennas = [Antenna("%03d" % idx, geo = ant)
        for idx, ant in enumerate(coordinates)]

    imageFileName = tempfile.gettempdir() + '/' + 'contour.png'
    image = observation.createContour(antennas, imageFileName)
    pixmap = QPixmap(imageFileName)
    label.setPixmap(pixmap.scaledToHeight(pixmap.height()))
    updateHorizontal(observation.getBoreSight().horizontal)
    # updateBaselineList(observation.getBaselines())
    projectedBaselines = observation.getProjectedBaselines()
    baselineCoordinates = np.array([pbl for pbl in projectedBaselines])
    uvCoord = baselineCoordinates / waveLength
    uvSamples = np.concatenate((uvCoord, -uvCoord))
    # np.save("/tmp/uvSamples", uvSamples)
    updateUVPlane(uvSamples)
    resetPackState()
    beamSizeFactor = observation.getBeamSizeFactor()
    beamSizeEdit.blockSignals(True)
    beamSizeEdit.setValue(beamSizeFactor)
    beamSizeEdit.blockSignals(False)


    #if os.path.exists(imageFileName):
    #    os.remove(imageFileName)
def onClickedFitsButton():

    from astropy import wcs
    from astropy.io import fits

    WCS = observation.getWCS()
    w = wcs.WCS(naxis=2)
    data = observation.getImageData()

    # "Airy's zenithal" projection
    w.wcs.crpix = WCS['crpix']
    w.wcs.cdelt = WCS['cdelt']
    w.wcs.crval = WCS['crval']
    w.wcs.ctype = WCS['ctype']
    # w.wcs.set_pv([(2, 1, 45.0)])

    header = w.to_header()

    hdu = fits.PrimaryHDU(header=header, data=data)
    # Save to FITS file

    fileName = QFileDialog.getSaveFileName(parent=None, caption='Save File')

    hdu.writeto(str(fileName), overwrite=True)



def onClickedAddGeoButton():
    longitude = longitudeCoord.text()
    latitude = latitudeCoord.text()
    height = "1035"
    if longitude == '' or latitude == '':
        return
    position = coordinateList.rowCount()
    addRowToCoordinateList('*', position, latitude, longitude, height)
    axis.addDots([[float(latitude), float(longitude)],])
    updateContour()

def onClickedOutputButton():
    coordinates, hidden = collectCoordinates()
    observeTime = observation.getObserveTime()
    source = observation.getBoreSight().equatorial

    fileName = QFileDialog.getSaveFileName(parent=None, caption='Save File')
    with open(fileName, 'wb') as paraFile:
        pickle.dump([coordinates, source, observeTime, hidden], paraFile)


def onClickedImportButton():

    fileName = QFileDialog.getOpenFileName()
    if str(fileName).endswith('.csv'):
        antennaCoords = np.loadtxt(str(fileName), delimiter=',')
        onlyAntenna = True
    else:
        with open(fileName, 'rb') as paraFile:
            paras = pickle.load(paraFile)
        antennaCoords = paras[0]
        source = paras[1]
        observeTime = paras[2]
        if len(paras) > 3:
            hidden = paras[3]
        else:
            hidden = []
        onlyAntenna = False

    onClickedDelAllButton()
    dateTimeEdit.blockSignals(True)
    azimuthCoord.blockSignals(True)
    elevationCoord.blockSignals(True)
    RACoord.blockSignals(True)
    DECCoord.blockSignals(True)

    position = 0
    dots = []
    for latitude, longitude, height in antennaCoords:
        if position in hidden:
            hiddenFlag = 'hidden'
        else:
            hiddenFlag = ''
            dots.append([float(latitude), float(longitude),])
        addRowToCoordinateList(position, position, latitude, longitude, height, hiddenFlag)
        position += 1

    axis.addDots(dots)

    if onlyAntenna == False:
        if source[0] != 'str':
            sourceDeg = source
            source = convertBoresightToHour(source)

        RACoord.setText(str(source[0]))
        RACoord.setCursorPosition(0)
        DECCoord.setText(str(source[1]))
        DECCoord.setCursorPosition(0)

        observation.setBoreSight(sourceDeg, Boresight.EquatorialFrame)

        newDateTime = QDateTime(observeTime.year,
                observeTime.month, observeTime.day,
                observeTime.hour, observeTime.minute,
                observeTime.second, observeTime.microsecond/1000)
        dateTimeEdit.setDateTime(newDateTime)
        observation.setObserveTime(observeTime)

    updateContour()
    updateHorizontal(observation.getBoreSight().horizontal)

    dateTimeEdit.blockSignals(False)
    dateTimeEdit.blockSignals(False)
    azimuthCoord.blockSignals(False)
    RACoord.blockSignals(False)
    DECCoord.blockSignals(False)


def onClickedDelAllButton():
    coordinateList.setRowCount(0)
    # baselineList.setRowCount(0)
    label.setPixmap(blankImage)
    axis.clearDots()
    UVPlane.clearDots()

def onClickedPackButton2():
    if not hasattr(onClickedPackButton2, "state"):
        onClickedPackButton2.state = 0
    if not hasattr(onClickedPackButton2, "fittedImage"):
        onClickedPackButton2.fittedImage = None
    if not hasattr(onClickedPackButton2, "newData"):
        onClickedPackButton2.newData = True

    tempdir = tempfile.gettempdir()

    if onClickedPackButton2.state == 1:
        onClickedPackButton2.state = 0
        fittedImage = onClickedPackButton2.fittedImage
        if(fittedImage is not None):
            label.setPixmap(fittedImage.scaledToHeight(fittedImage.height()))
        else:
            pixmap = QPixmap(tempdir + '/contour.png')
            label.setPixmap(pixmap.scaledToHeight(pixmap.height()))
        return

    if  onClickedPackButton2.newData == False:
        onClickedPackButton2.state = 1
        pixmap = QPixmap(tempdir + '/pack.png')
        label.setPixmap(pixmap.scaledToHeight(pixmap.height()))
        return

    number = observation.getBaselinesNumber()
    # amplitude = observation.getAmplitude()
    # coordinates = observation.getBeamCoordinates()
    # center, angle, axisH, axisV = fitEllipseBeam(coordinates, amplitude, number*0.4)
    # center = np.rad2deg(observation.getHorizontal())
    center = observation.getBoreSight().equatorial
    imageLength = observation.getImageLength()
    sizeInfo = observation.getBeamAxis()
    if sizeInfo is None:
        return
    else:
        axisH2, axisV2, angle2, plotRange = sizeInfo
    # print sizeInfo
    if center == []: return
    # print 'axisLengthFit: ', axisH, axisV, np.rad2deg(angle), center
    # plotBeamFit(coordinates, center, np.rad2deg(angle), axisH, axisV)
    # ==============================
    # imageDensity = observation.getImageDensity()
    # step=imageLength*1.0/imageDensity
    # ellipseCenter = [center[0] + step/2., center[1] - step/2.]
    # plotBeamFit(imageLength, center, ellipseCenter, angle2, axisH2, axisV2,
            # fileName = tempdir + '/fit.png')
    # bottomImage = QImage(tempdir + '/contour.png')
    # topImage = QImage(tempdir + '/fit.png')
    # fittedImage = QPixmap.fromImage(overlayImage(bottomImage, topImage))
    # ===============================
    imageFileName = tempfile.gettempdir() + '/' + 'contourfit.png'
    imageData = observation.getImageData()
    boresight = observation.getBoreSight().equatorial
    resolution = observation.getResolution()
    widthH = axisH2/resolution
    widthV = axisV2/resolution
    # boresight = (10, 9)
    # widthH = axisH2
    # widthV = axisV2
    plotBeamWithFit(imageData, center, plotRange, widthH, widthV, angle2,
                    fileName=imageFileName, colormap = True)
    pixmap = QPixmap(imageFileName)
    label.setPixmap(pixmap.scaledToHeight(pixmap.height()))
    fittedImage = pixmap
    # ===============================
    onClickedPackButton2.fittedImage = fittedImage
    label.setPixmap(fittedImage.scaledToHeight(fittedImage.height()))
    coordinates, beamRadius = ellipseCompact(400, axisH2, axisV2, angle2, 10)
    # ======================
    beamNumberTidal = coordinates.shape[0]
    # factor = float(packSizeEdit.text())
    factor = 1
    primaryBeamRadius = np.rad2deg(1.22*waveLength/13.5)/2. * factor
    coordinatesPrimary = ellipseGrid(primaryBeamRadius, beamRadius, beamRadius, 0).T
    beamNumberPrimary = len(coordinatesPrimary)


    # plotPackedBeam(coordinatesPrimary, 0, beamRadius, beamRadius, (0,0), primaryBeamRadius,
            # fileName = tempdir + '/primaryPack.png')
    # print("tidal: num:%d, axisH:%f, axisV:%f, radius:%f" % (beamNumberTidal, axisH2, axisV2, beamRadius))
    print("tidal: num:%d, axisH:%f, axisV:%f, radius:%f" % (beamNumberTidal, axisH2, axisV2, beamRadius)),
    print(", angle: %f" % angle2)
    # print("primary: num:%d" % beamNumberPrimary)
    # beamArea = np.pi*axisH*axisV
    # primaryBeamArea = np.pi*(beamRadius**2)
    # ratio = beamNumber*beamArea/primaryBeamArea
    # print(beamRadius)
    # print("%dx%f/%f=%f" % (beamNumber, beamArea, primaryBeamArea, ratio))
    # ======================
    maxDistance = np.max(
            np.sqrt(np.sum(np.square(coordinates), axis = 1)))
    upper_left_pixel = [-maxDistance, maxDistance] # x,y
    bottom_right_pixel = [maxDistance, -maxDistance] # x,y
    coordinates_equatorial = convert_pixel_coordinate_to_equatorial(
        [upper_left_pixel, bottom_right_pixel], center)

    equatorial_range = [
        coordinates_equatorial[0][0], coordinates_equatorial[1][0], # left, right
        coordinates_equatorial[0][1], coordinates_equatorial[1][1]] # up, bottom
    pixel_range = [upper_left_pixel[0], bottom_right_pixel[0], # left, right
                   upper_left_pixel[1], bottom_right_pixel[1]] # up, bottom

    plotPackedBeam(coordinates, angle2, axisH2, axisV2, center,
            equatorial_range, pixel_range, beamRadius,
            fileName = tempdir + '/pack.png', HD = False)
    pixmap = QPixmap(tempdir + '/pack.png')
    label.setPixmap(pixmap.scaledToHeight(pixmap.height()))
    onClickedPackButton2.state = 1
    onClickedPackButton2.newData = False

    #========parallactic
    # try:
        # parallacticData = np.load('parallactic.npy')
    # except:
        # parallacticData = np.array([])
    # parallacticData = np.append(parallacticData, [coordinates, np.rad2deg(angle2), axisH2, axisV2, beamRadius], axis=0)
    # np.save('parallactic.npy', parallacticData)

    #=========== overlaps
    # rotationOffset = float(packSizeEdit.text())
    # angleOffset = rotationOffset*60/(3600.*24) * 2 * np.pi
    # overlapCounter = calculateBeamOverlaps(
            # coordinates, beamRadius, axisH2, axisV2, angle2 + angleOffset)

    #=========== heapMap
    # tidalOffset = coordinates
    # mapGrid = []

    # np.savetxt('setcenter', coordinatesPrimary)
    # for beamSetCenter in coordinatesPrimary:
        # mapGrid.append(tidalOffset + beamSetCenter)

    # mapGrid = np.array(mapGrid)
    # shapes = mapGrid.shape
    # mapSum = np.zeros((shapes[0], shapes[1]))
    # for beamSetCenter in coordinatesPrimary:
        # mapSum += Gaussian2DPDF(mapGrid[:,:,0], mapGrid[:,:,1], beamSetCenter[0], beamSetCenter[1], primaryBeamRadius, primaryBeamRadius)
    # plotSkyHeatMap(mapGrid[:,:,0], mapGrid[:,:,1], mapSum)
    #=========== heapMap
    # heats = generateSkyHeatMap(primaryBeamRadius*1.5, coordinatesPrimary, beamRadius, primaryBeamRadius)
    # gridX = np.arange(1000)
    # gridY = np.arange(1000)
    # plotSkyHeatMap(heats)

def generateSkyHeatMap(mapRadius, beamSetCenters, beamSetRadius, primaryBeamRadius):
    gridNum = 1000.
    step = mapRadius/(gridNum/2)
    # grids = np.mgrid[0:2*mapRadius:step, 0:2*mapRadius:step]
    grids = np.mgrid[0:gridNum:1, 0:gridNum:1]
    beamSetGridLen = int(round(beamSetRadius/step))
    heapMap = np.zeros((1000, 1000))
    diameter = beamSetGridLen*2 + 1
    for beamSetCenter in beamSetCenters:
        xCenterIdx = int(round((beamSetCenter[0] + mapRadius)/(2*mapRadius)*gridNum))
        yCenterIdx = int(round((mapRadius - beamSetCenter[1])/(2*mapRadius)*gridNum))
        xStart = xCenterIdx - beamSetGridLen
        yStart = yCenterIdx - beamSetGridLen
        xEnd = xCenterIdx + beamSetGridLen
        yEnd = yCenterIdx + beamSetGridLen
        xGrid = grids[1][yStart:yEnd+1, xStart:xEnd+1].reshape(diameter*diameter)
        yGrid = grids[0][yStart:yEnd+1, xStart:xEnd+1].reshape(diameter*diameter)
        heat = Gaussian2DPDF(xGrid, yGrid, xCenterIdx, yCenterIdx, beamSetGridLen, beamSetGridLen)
        block = heapMap[yStart:yEnd+1, xStart:xEnd+1].reshape(diameter*diameter)
        blockLength = len(block)
        for i in range(blockLength):
            block[i] = heat[i] if heat[i] > block[i] else block[i]
        heapMap[yStart:yEnd+1, xStart:xEnd+1] = block.reshape(diameter,diameter)
        # np.savetxt("heapMap", heapMap)

    return heapMap


def Gaussian2DPDF(x, y, xMean, yMean, xSigma, ySigma):
    return np.exp(-((x-xMean)**2/(2*(xSigma**2)) + (y-yMean)**2/(2*(ySigma**2))))


def overlayImage(bottom, top):
    painter = QPainter()
    painter.begin(bottom)
    painter.drawImage(0, 0, top)
    painter.end()

    return bottom


def onClickedAtCoordinateList(row, column):
    x = float(str(coordinateList.item(row, 0).text()))
    y = float(str(coordinateList.item(row, 1).text()))

    if column == 4:
        coordinateList.removeRow(row)
        axis.removeDots([[x,y],])
        updateContour()
        return
    elif column == 3:
        item = coordinateList.item(row, column)
        if item is not None and str(item.text()) == 'hidden':
            coordinateList.setItem(row, column, QTableWidgetItem(''))
            axis.addDots([[x,y],])
        else:
            coordinateList.setItem(row, column, QTableWidgetItem('hidden'))
            axis.removeDots([[x,y],])
        updateContour()
        return

    items = coordinateList.selectedItems()
    dots = []
    for item in items:
        x = float(str(coordinateList.item(item.row(), 0).text()))
        y = float(str(coordinateList.item(item.row(), 1).text()))
        dots.append([x,y])

    axis.addHighLightDots(dots)
    return

def onPackSizeChanged():
    resetPackState()

def onCoordinateListSelectionChanged():
    items = coordinateList.selectedItems()
    dots = []
    for item in items:
        if item.column() == 3: continue
        x = float(str(coordinateList.item(item.row(), 0).text()))
        y = float(str(coordinateList.item(item.row(), 1).text()))
        dots.append([x,y])

    axis.addHighLightDots(dots)

# class onItemChangedAtCoordinateList(QObject):
    # def eventFilter(self, receiver, event):
        # print('fired')
        # if(event.type() == QEvent.Enter):
            # print('enter pressed')
            # return True
        # else:
            # return super(MyEventFilter,self).eventFilter(receiver, event)


# def updateBaselineList(baselines):
    # baselineList.setRowCount(0)
    # if baselines == None:return
    # index = 0
    # for baseline in baselines:
        # length = '{:6.2f}'.format(np.linalg.norm(baseline))
        # baselineList.insertRow(index)
        # vector = ' '.join(['{: 9.1f}'.format(i) for i in baseline])
        # baselineList.setItem(index, 0, QTableWidgetItem(vector))
        # baselineList.setItem(index, 1, QTableWidgetItem(str(length)))
        # index += 1

def updateUVPlane(uvSamples):
    # print baselines
    # if baselines == None:return

    UVPlane.clearDots()
    #swap x,y because the Cartesian class use (y,x) index.
    UVPlane.addDots(uvSamples[:,[1,0]])


np.set_printoptions(precision=3)
loggerFormat = '%(asctime)-15s %(filename)s  %(message)s'
logging.basicConfig(format = loggerFormat, level=logging.WARNING)
logger = logging.getLogger()

parser = argparse.ArgumentParser()
parser.add_argument("-v", "--verbose", help="increase output verbosity",
        action="store_true")

args = parser.parse_args()
if args.verbose:
    logger.setLevel(logging.INFO)

'''MeerKAT coordinates'''
# the values provided by http://public.ska.ac.za/meerkat
# Longitude, Latitude, Height
arrayReferece = Antenna('ref', [-30.71106, 21.44389, 1035])
'''observation time in UTC'''
observationTime = QDateTime.currentDateTime().toPyDateTime()
'''observation waveLength in meter'''
waveLength = 0.21

defaultBeamSizeFactor = 1
defaultBeamNumber = 400
defaultBoreSight = ("1:25:46.5336", "-30:42:39.815999")

observation = InterferometryObservation(arrayReferece, waveLength)
observation.setBoreSight(convertBoresightToDegree(defaultBoreSight), Boresight.EquatorialFrame)
# observation.setBeamSizeFactor(defaultBeamSizeFactor)
# observation.setBeamNumber(defaultBeamNumber)
observation.setObserveTime(observationTime)
observation.setInterpolating(True)
observation.setAutoZoom(True)

imageData = None



a = QApplication(sys.argv)

w = QWidget()
w.setWindowTitle("WaveRider")

axis =Cartesian(w)
axis.setCenter(arrayReferece.geo, 0.04, swapXY=True)
axis.move(500, 10)

label = QLabel(w)
blankImage = QPixmap(400, 300)
blankImage.fill(Qt.white)
label.setPixmap(blankImage)
label.setScaledContents(True)
label.move(10, 10)

longitudeCoordLabel = QLabel(w)
longitudeCoordLabel.setText('Longitude')
longitudeCoordLabel.move(100, 380)
latitudeCoordLabel = QLabel(w)
latitudeCoordLabel.setText('Latitude')
latitudeCoordLabel.move(10, 380)

longitudeCoord = QLineEdit(w)
latitudeCoord = QLineEdit(w)
longitudeCoord.resize(80,30)
latitudeCoord.resize(80,30)
longitudeCoord.move(100, 400)
latitudeCoord.move(10, 400)

addGeoButton = QPushButton('Add', w)
addGeoButton.clicked.connect(onClickedAddGeoButton)
addGeoButton.resize(40, 30)
addGeoButton.move(180, 400)

importButton = QPushButton('Import', w)
importButton.resize(60, 30)
importButton.clicked.connect(onClickedImportButton)
importButton.move(230, 400)

outputButton = QPushButton('Output', w)
outputButton.resize(60, 30)
outputButton.move(230, 360)
outputButton.clicked.connect(onClickedOutputButton)


DeleteAllButton = QPushButton('Delete All', w)
DeleteAllButton.clicked.connect(onClickedDelAllButton)
DeleteAllButton.move(300, 400)

FitsButton = QPushButton('Fits', w)
FitsButton.resize(50, 30)
FitsButton.move(400, 360)
FitsButton.clicked.connect(onClickedFitsButton)

PackButton = QPushButton('Pack', w)
PackButton.clicked.connect(onClickedPackButton2)
PackButton.resize(50, 30)
PackButton.move(400, 400)

# packSizeLabel = QLabel(w)
# packSizeLabel.setText('Div')
# packSizeLabel.move(450, 380)
# packSizeEdit = QSpinBox(w)
# packSizeEdit.move(450, 400)
# packSizeEdit.setValue(0)
# packSizeEdit.setMinimum(0)
# packSizeEdit.setMaximum(999)
# packSizeEdit.valueChanged.connect(onRotationChanged)



dateTimeLabel = QLabel(w)
dateTimeLabel.setText('UTC Time')
dateTimeLabel.move(500, 320)
dateTimeEdit = QDateTimeEdit(w)
dateTimeEdit.move(500, 340)
dateTimeEdit.setDisplayFormat("yyyy.MM.dd hh:mm:ss.zzz")
dateTimeEdit.setDateTime(QDateTime.currentDateTime())
dateTimeEdit.dateTimeChanged.connect(onDateTimeChanged)
dateTimeEdit.setWrapping(True)

beamSizeLabel = QLabel(w)
beamSizeLabel.setText('Zoom')
beamSizeLabel.move(710, 320)
beamSizeEdit = QSpinBox(w)
beamSizeEdit.move(710, 340)
# beamSizeEdit.setValue(defaultBeamSizeFactor)
beamSizeEdit.setValue(observation.getBeamSizeFactor())
beamSizeEdit.setMinimum(1)
beamSizeEdit.setMaximum(99)
beamSizeEdit.valueChanged.connect(onBeamSizeChanged)

autoZoomLabel = QLabel(w)
autoZoomLabel.setText(u"\u25F1")
autoZoomLabel.setToolTip('Auto Zoom')
autoZoomLabel.move(763, 320)
autoZoomCheckbox = QCheckBox(w)
autoZoomCheckbox.move(760, 345)
autoZoomCheckbox.setToolTip('Auto Zoom')
autoZoomCheckbox.setCheckState(Qt.Checked)
autoZoomCheckbox.stateChanged.connect(onAutoZoomOptionChanged)



beamNumberLabel = QLabel(w)
beamNumberLabel.setText('Pixels')
beamNumberLabel.move(780, 320)
beamNumberEdit = QLineEdit(w)
beamNumberEdit.move(785, 340)
beamNumberEdit.resize(40,30)
# beamNumberEdit.setText(str(defaultBeamNumber))
beamNumberEdit.setText(str(int(observation.getBeamNumber())))
beamNumberEdit.editingFinished.connect(onBeamNumberChanged)

interpolateLabel = QLabel(w)
interpolateLabel.setText(u"\u25A6")
interpolateLabel.setToolTip('Interpolation')
interpolateLabel.move(833, 320)
interpolateCheckbox = QCheckBox(w)
interpolateCheckbox.move(830, 345)
interpolateCheckbox.setToolTip('Interpolation')
interpolateCheckbox.setCheckState(Qt.Checked)
interpolateCheckbox.stateChanged.connect(onInterpolateOptionChanged)


RACoordLabel = QLabel(w)
RACoordLabel.setText('RA')
RACoordLabel.move(500, 380)
DECCoordLabel = QLabel(w)
DECCoordLabel.setText('DEC')
DECCoordLabel.move(590, 380)


RACoord = QLineEdit(w)
RACoord.move(500, 400)
RACoord.resize(80,30)
DECCoord = QLineEdit(w)
DECCoord.resize(80,30)
DECCoord.move(590, 400)
RACoord.setText(str(defaultBoreSight[0]))
DECCoord.setText(str(defaultBoreSight[1]))
RACoord.home(True)
DECCoord.home(True)
DECCoord.editingFinished.connect(onBoreSightUpdated)
RACoord.editingFinished.connect(onBoreSightUpdated)


azimuthCoordLabel = QLabel(w)
azimuthCoordLabel.setText('Azimuth')
azimuthCoordLabel.move(680, 380)
elevationCoordLabel = QLabel(w)
elevationCoordLabel.setText('Elevation')
elevationCoordLabel.move(770, 380)



azimuthCoord = QLineEdit(w)
azimuthCoord.move(680, 400)
azimuthCoord.resize(80,30)
azimuthCoord.editingFinished.connect(onHorizontalUpdated)
elevationCoord = QLineEdit(w)
elevationCoord.resize(80,30)
elevationCoord.move(770, 400)
elevationCoord.editingFinished.connect(onHorizontalUpdated)

coordinateListLabel = QLabel(w)
coordinateListLabel.setText('Antennas')
coordinateListLabel.move(10, 440)

coordinateList = QTableWidget(w)
coordinateList.setColumnCount(5)
coordinateList.resize(480, 300)
coordinateList.move(10, 460)
coordinateList.setHorizontalHeaderItem(0, QTableWidgetItem('latitude'))
coordinateList.setHorizontalHeaderItem(1, QTableWidgetItem('longitude'))
coordinateList.setHorizontalHeaderItem(2, QTableWidgetItem('height'))
coordinateList.setHorizontalHeaderItem(3, QTableWidgetItem('hide'))
coordinateList.setHorizontalHeaderItem(4, QTableWidgetItem('delete'))
coordinateList.cellClicked.connect(onClickedAtCoordinateList)
coordinateList.itemSelectionChanged.connect(onCoordinateListSelectionChanged)
# coordinateList.itemChanged.connect(onItemChangedAtCoordinateList)
# coordinateList.installEventFilter(onItemChangedAtCoordinateList())
coordinateList.setFocus()

UVPlaneLabel = QLabel(w)
UVPlaneLabel.setText('UV plane')
UVPlaneLabel.move(500, 440)

UVPlane = miniCartesian(w)
UVPlane.setCenter((0,0), 40000)
UVPlane.move(500, 460)


w.resize(860, 860)

w.show()


sys.exit(a.exec_())
