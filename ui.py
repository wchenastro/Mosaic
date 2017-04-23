#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os
import numpy as np

from PyQt4.QtGui import *
from PyQt4.QtCore import *

from interferometer import InterferometryObservation
from beamShape import fitEllipseBeam, plotPackedBeam, plotBeamFit
from createBeam import ellipseGrid, ellipseCompact

class Cartesian(QWidget):

    painter = None
    dots = []
    azimuth = 0
    elevation = 0
    width = 300.
    height = 300.
    halfWidth = 150.
    halfHeight = 150.
    xStart = 21.391
    xEnd = 21.431
    yStart = -30.741
    yEnd = -30.701

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
        painter.begin(self)
        painter.drawRect(0, 0, self.width, self.height)
        '''centered horzontal line'''
        painter.drawLine(0, self.halfHeight, self.width, self.halfHeight)
        '''centered vertical line'''
        painter.drawLine(self.halfWidth, self.height, self.halfWidth, 0)
        painter.drawEllipse(QPoint(self.halfWidth, self.halfHeight),
                self.halfWidth, self.halfHeight)
        painter.drawEllipse(QPoint(self.halfWidth, self.halfHeight),
                (90.-self.elevation)/90.*self.halfWidth,
                (90.-self.elevation)/90.*self.halfHeight)
        angleX, angleY = self.angleToCartesian(self.azimuth, self.halfWidth)
        painter.drawLine(self.halfWidth, self.halfHeight, angleX, angleY)

        for dot in self.dots:
            painter.drawEllipse(dot[0], dot[1],3,3)

        painter.end()

    def sizeHint(self):
        return QSize(301, 301)

    def mouseMoveEvent(self, event):
        longitude, altitude = self.pixelCoordinateConv((event.x(), event.y()), 'toCoord')
        longitudeCoord.setPlaceholderText('{:6.4f}'.format(longitude))
        altitudeCoord.setPlaceholderText('{:6.4f}'.format(altitude))

    def addDots(self, dots):
        for longitude, altitude in dots:
            dot = self.pixelCoordinateConv((longitude, altitude), 'toPixel')
            self.dots.append(dot)
        self.update()

    def mouseReleaseEvent(self, event):
        self.dots.append([event.x(), event.y()])
        longitude, altitude = self.pixelCoordinateConv((event.x(), event.y()), 'toCoord')
        addRowToCoordinateList(longitude, altitude)
        self.update()
        updateCountour()
        # updateBaselineList(observation.getBaselines())

    def angleToCartesian(self, angleDeg, radius):
        angle = np.deg2rad(angleDeg)
        if angle < np.pi/2:
            x = np.sin(angle)*radius + radius
            y = radius - np.cos(angle)*radius
        elif angle > np.pi/2 and angle < np.pi:
            x = radius + np.sin(angle)*radius
            y = radius - np.cos(angle)*radius
        elif angle > np.pi and angle < np.pi*1.5:
            x = radius + np.sin(angle)*radius
            y = -np.cos(angle)*radius + radius
        elif angle > np.pi*1.5 and angle < np.pi*2:
            x = radius + np.sin(angle)*radius
            y = radius - np.cos(angle)*radius

        return x, y

    def clearDots(self):
        self.dots = []
        self.update()

    def removeDots(self, dotsToRemove):
        for longitude, altitude in dotsToRemove:
            dot = self.pixelCoordinateConv((longitude, altitude), 'toPixel')
            self.dots.remove(dot)
        self.update()

    def setAzAlt(self, horizontalCoord):
        self.azimuth = horizontalCoord[0]
        self.elevation = np.abs(horizontalCoord[1])
        self.update()

def updateHorizontal(horizontalCoord):
    if horizontalCoord == []: return
    azimuth = np.rad2deg(horizontalCoord[0])
    elevation = np.rad2deg(horizontalCoord[1])
    axis.setAzAlt([azimuth, elevation])
    azimuthCoord.setText('{:6.4f}'.format(azimuth))
    elevationCoord.setText('{:6.4f}'.format(elevation))

def onBoreSightUpdated():
    beamBoreSight = (float(RACoord.text()), float(DECCoord.text()))
    observation.setBoreSight(beamBoreSight)
    updateCountour()
    updateHorizontal(observation.getHorizontal())

def onBeamSizeChanged():
    observation.setBeamSizeFactor(beamSizeEdit.value())
    updateCountour()

def onBeamNumberChanged():
    observation.setBeamNumber(float(beamNumberEdit.text()))
    updateCountour()

def onInterpolateOptionChanged(state):
    if state == Qt.Checked:
        observation.setInterpolating(True)
    else:
        observation.setInterpolating(False)
    updateCountour()

def onDateTimeChanged(dateTime):
    observation.setObserveTime(dateTime.toPyDateTime())
    updateCountour()
    updateHorizontal(observation.getHorizontal())

def addRowToCoordinateList(longitude, altitude):
    rowCount = coordinateList.rowCount()
    coordinateList.insertRow(rowCount)
    coordinateList.setItem(rowCount, 0, QTableWidgetItem('{:6.4f}'.format(longitude)))
    coordinateList.setItem(rowCount, 1, QTableWidgetItem('{:6.4f}'.format(altitude)))
    coordinateList.setItem(rowCount, 3, QTableWidgetItem('-'))

def resetPackState():
    onClickedPackButton2.newData = True
    onClickedPackButton2.state = 0
    onClickedPackButton2.fittedImage = None

def updateCountour():
    # print('updateCountour')

    # if not hasattr(updateCountour, "lastObservetime"):
        # updateCountour.lastObservetime = None

    # observeTime = dateTimeEdit.dateTime().toPyDateTime()
    # if(updateCountour.lastObservetime != observeTime):
        # observation.setObserveTime(observeTime)
        # updateCountour.lastObservetime = observeTime

    rowCount = coordinateList.rowCount()
    if rowCount < 2:
        label.setPixmap(blankImage)
        return
    coordinates = []
    for row in range(rowCount):
        longitude = float(str(coordinateList.item(row, 0).text()))
        altitude = float(str(coordinateList.item(row, 1).text()))
        coordinates.append([longitude, altitude, 0])

    observation.createContour(coordinates)
    pixmap = QPixmap(os.getcwd() + '/contour.png')
    label.setPixmap(pixmap.scaledToHeight(pixmap.height()))
    updateBaselineList(observation.getBaselines())
    updateHorizontal(observation.getHorizontal())
    resetPackState()
    beamSizeFactor = observation.getBeamSizeFactor()
    beamSizeEdit.blockSignals(True)
    beamSizeEdit.setValue(beamSizeFactor)
    beamSizeEdit.blockSignals(False)

def onClickedAddGeoButton():
    longitude = longitudeCoord.text()
    altitude = altitudeCoord.text()
    if longitude == '' or altitude == '':
        return
    rowCount = coordinateList.rowCount()
    coordinateList.insertRow(rowCount)
    coordinateList.setItem(rowCount, 0, QTableWidgetItem(longitude))
    coordinateList.setItem(rowCount, 1, QTableWidgetItem(altitude))
    coordinateList.setItem(rowCount, 3, QTableWidgetItem('-'))
    axis.addDots([[float(longitude), float(altitude)],])
    updateCountour()

def onClickedDelAllButton():
    coordinateList.setRowCount(0)
    baselineList.setRowCount(0)
    label.setPixmap(blankImage)
    axis.clearDots()

def onClickedPackButton():
    if not hasattr(onClickedPackButton, "state"):
        onClickedPackButton.state = 0


    if onClickedPackButton.state == 1:
        onClickedPackButton.state = 0
        pixmap = QPixmap(os.getcwd() + '/contour.png')
        label.setPixmap(pixmap.scaledToHeight(pixmap.height()))
        return

    onClickedPackButton.state = 1
    number = observation.getBaselinesNumber()
    angle, axisH, axisV = fitEllipseBeam(number*0.4)
    divider = float(packSizeEdit.value())
    beamRadius = np.rad2deg(1.22*waveLength/13.5)/2/divider
    coordinates = ellipseGrid(beamRadius, axisH, axisV, angle, write=False)
    # plotPackedBeam('ellipsePack', np.rad2deg(angle), axisH, axisV, beamRadius, 'png')
    plotPackedBeam(coordinates, np.rad2deg(angle), axisH, axisV, beamRadius, 'png')
    pixmap = QPixmap(os.getcwd() + '/pack.png')
    label.setPixmap(pixmap.scaledToHeight(pixmap.height()))

def onClickedPackButton2():
    if not hasattr(onClickedPackButton2, "state"):
        onClickedPackButton2.state = 0
    if not hasattr(onClickedPackButton2, "fittedImage"):
        onClickedPackButton2.fittedImage = None
    if not hasattr(onClickedPackButton2, "newData"):
        onClickedPackButton2.newData = True


    if onClickedPackButton2.state == 1:
        onClickedPackButton2.state = 0
        fittedImage = onClickedPackButton2.fittedImage
        if(fittedImage != None):
            label.setPixmap(fittedImage.scaledToHeight(fittedImage.height()))
        else:
            pixmap = QPixmap(os.getcwd() + '/contour.png')
            label.setPixmap(pixmap.scaledToHeight(pixmap.height()))
        return

    if  onClickedPackButton2.newData == False:
        onClickedPackButton2.state = 1
        pixmap = QPixmap(os.getcwd() + '/pack.png')
        label.setPixmap(pixmap.scaledToHeight(pixmap.height()))
        return

    number = observation.getBaselinesNumber()
    amplitude = observation.getAmplitude()
    coordinates = observation.getBeamCoordinates()
    center, angle, axisH, axisV = fitEllipseBeam(coordinates, amplitude, number*0.4)
    plotBeamFit(coordinates, center, np.rad2deg(angle), axisH, axisV)
    bottomImage = QImage(os.getcwd() + '/contour.png')
    topImage = QImage(os.getcwd() + '/fit.png')
    fittedImage = QPixmap.fromImage(overlayImage(bottomImage, topImage))
    onClickedPackButton2.fittedImage = fittedImage
    label.setPixmap(fittedImage.scaledToHeight(fittedImage.height()))
    coordinates, beamRadius = ellipseCompact(400, axisH, axisV, angle, 10)
    beamNumber = coordinates.shape[0]
    beamArea = np.pi*axisH*axisV
    primaryBeamArea = np.pi*(beamRadius**2)
    ratio = beamNumber*beamArea/primaryBeamArea
    # print(beamRadius)
    # print("%dx%f/%f=%f" % (beamNumber, beamArea, primaryBeamArea, ratio))
    plotPackedBeam(coordinates, np.rad2deg(angle), axisH, axisV, beamRadius)
    pixmap = QPixmap(os.getcwd() + '/pack.png')
    label.setPixmap(pixmap.scaledToHeight(pixmap.height()))
    onClickedPackButton2.state = 1
    onClickedPackButton2.newData = False

def overlayImage(bottom, top):
    painter = QPainter()
    painter.begin(bottom)
    painter.drawImage(0, 0, top)
    painter.end()

    return bottom


def onClickedAtCoordinateList(row, column):
    if column == 3:
        x = float(str(coordinateList.item(row, 0).text()))
        y = float(str(coordinateList.item(row, 1).text()))
        coordinateList.removeRow(row)
        axis.removeDots([[x,y],])
        updateCountour()


# class onItemChangedAtCoordinateList(QObject):
    # def eventFilter(self, receiver, event):
        # print('fired')
        # if(event.type() == QEvent.Enter):
            # print('enter pressed')
            # return True
        # else:
            # return super(MyEventFilter,self).eventFilter(receiver, event)


def updateBaselineList(baselines):
    baselineList.setRowCount(0)
    if baselines == None:return
    index = 0
    for baseline in baselines:
        length = '{:6.2f}'.format(np.linalg.norm(baseline))
        baselineList.insertRow(index)
        vector = ' '.join(['{: 9.1f}'.format(i) for i in baseline])
        baselineList.setItem(index, 0, QTableWidgetItem(vector))
        baselineList.setItem(index, 1, QTableWidgetItem(str(length)))
        index += 1

np.set_printoptions(precision=3)

'''MeerKAT coordinates'''
# the values provided by http://public.ska.ac.za/meerkat
arrayRefereceGEODET = (21.44389,-30.71317, 0)
'''observation time in UTC'''
observationTime = QDateTime.currentDateTime().toPyDateTime()
'''observation waveLength in meter'''
waveLength = 0.21

defaultBeamSizeFactor = 90
defaultBeamNumber = 400
defaultBoreSight = (21.411, -30.721)

observation = InterferometryObservation(arrayRefereceGEODET,
        observationTime, waveLength)
observation.setBoreSight(defaultBoreSight)
observation.setBeamSizeFactor(defaultBeamSizeFactor)
observation.setBeamNumber(defaultBeamNumber)
observation.setInterpolating(True)

a = QApplication(sys.argv)

w = QWidget()
w.setWindowTitle("WaveRider")

axis =Cartesian(w)
axis.move(500, 10)

label = QLabel(w)
blankImage = QPixmap(400, 300)
blankImage.fill(Qt.white)
label.setPixmap(blankImage)
label.move(10, 10)

longitudeCoordLabel = QLabel(w)
longitudeCoordLabel.setText('Longitude')
longitudeCoordLabel.move(10, 380)
altitudeCoordLabel = QLabel(w)
altitudeCoordLabel.setText('Altitude')
altitudeCoordLabel.move(100, 380)

longitudeCoord = QLineEdit(w)
altitudeCoord = QLineEdit(w)
longitudeCoord.resize(80,30)
altitudeCoord.resize(80,30)
longitudeCoord.move(10, 400)
altitudeCoord.move(100, 400)

addGeoButton = QPushButton('Add', w)
addGeoButton.clicked.connect(onClickedAddGeoButton)
addGeoButton.move(200, 400)

DeleteAllButton = QPushButton('Delete All', w)
DeleteAllButton.clicked.connect(onClickedDelAllButton)
DeleteAllButton.move(300, 400)

PackButton = QPushButton('Pack', w)
PackButton.clicked.connect(onClickedPackButton2)
PackButton.resize(50, 30)
PackButton.move(400, 400)

packSizeLabel = QLabel(w)
packSizeLabel.setText('Div')
packSizeLabel.move(450, 380)
packSizeEdit = QSpinBox(w)
packSizeEdit.move(450, 400)
packSizeEdit.setValue(1)
packSizeEdit.setMinimum(1)
# packSizeEdit.valueChanged.connect(onPackSizeChanged)



dateTimeLabel = QLabel(w)
dateTimeLabel.setText('UTC Time')
dateTimeLabel.move(500, 320)
dateTimeEdit = QDateTimeEdit(w)
dateTimeEdit.move(500, 340)
dateTimeEdit.setDisplayFormat("dd.MM.yyyy hh:mm:ss.zzz")
dateTimeEdit.setDateTime(QDateTime.currentDateTime())
dateTimeEdit.dateTimeChanged.connect(onDateTimeChanged)
dateTimeEdit.setWrapping(True)

beamSizeLabel = QLabel(w)
beamSizeLabel.setText('Zoom')
beamSizeLabel.move(710, 320)
beamSizeEdit = QSpinBox(w)
beamSizeEdit.move(710, 340)
beamSizeEdit.setValue(defaultBeamSizeFactor)
beamSizeEdit.setMinimum(1)
beamSizeEdit.setMaximum(1000)
beamSizeEdit.valueChanged.connect(onBeamSizeChanged)

beamNumberLabel = QLabel(w)
beamNumberLabel.setText('Beams')
beamNumberLabel.move(760, 320)
beamNumberEdit = QLineEdit(w)
beamNumberEdit.move(760, 340)
beamNumberEdit.resize(50,30)
beamNumberEdit.setText(str(defaultBeamNumber))
beamNumberEdit.editingFinished.connect(onBeamNumberChanged)

interpolateLabel = QLabel(w)
interpolateLabel.setText('INTL.')
interpolateLabel.setToolTip('Interpolation')
interpolateLabel.move(820, 320)
interpolateCheckbox = QCheckBox(w)
interpolateCheckbox.move(820, 345)
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
RACoord.setText(str(defaultBoreSight[0]))
RACoord.editingFinished.connect(onBoreSightUpdated)
DECCoord = QLineEdit(w)
DECCoord.resize(80,30)
DECCoord.move(590, 400)
DECCoord.setText(str(defaultBoreSight[1]))
DECCoord.editingFinished.connect(onBoreSightUpdated)


azimuthCoordLabel = QLabel(w)
azimuthCoordLabel.setText('Azimuth')
azimuthCoordLabel.move(680, 380)
elevationCoordLabel = QLabel(w)
elevationCoordLabel.setText('Elevation')
elevationCoordLabel.move(770, 380)



azimuthCoord = QLineEdit(w)
azimuthCoord.move(680, 400)
azimuthCoord.resize(80,30)
elevationCoord = QLineEdit(w)
elevationCoord.resize(80,30)
elevationCoord.move(770, 400)

coordinateListLabel = QLabel(w)
coordinateListLabel.setText('Antennas')
coordinateListLabel.move(10, 440)

coordinateList = QTableWidget(w)
coordinateList.setColumnCount(4)
coordinateList.resize(480, 300)
coordinateList.move(10, 460)
coordinateList.setHorizontalHeaderItem(0, QTableWidgetItem('longitude'))
coordinateList.setHorizontalHeaderItem(1, QTableWidgetItem('latitude'))
coordinateList.setHorizontalHeaderItem(2, QTableWidgetItem(''))
coordinateList.setHorizontalHeaderItem(3, QTableWidgetItem('delete'))
coordinateList.cellClicked.connect(onClickedAtCoordinateList)
# coordinateList.itemChanged.connect(onItemChangedAtCoordinateList)
# coordinateList.installEventFilter(onItemChangedAtCoordinateList())
coordinateList.setFocus()

baselineListLabel = QLabel(w)
baselineListLabel.setText('Baselines')
baselineListLabel.move(500, 440)


baselineList = QTableWidget(w)
baselineList.setColumnCount(2)
baselineList.resize(350, 300)
baselineList.move(500, 460)
baselineList.setHorizontalHeaderItem(0, QTableWidgetItem('vector'))
baselineList.setHorizontalHeaderItem(1, QTableWidgetItem('length'))
baselineHeader = baselineList.horizontalHeader()
baselineHeader.setResizeMode(0, QHeaderView.ResizeToContents)

w.resize(860, 800)

w.show()


sys.exit(a.exec_())
