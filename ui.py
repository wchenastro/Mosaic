#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os
import numpy as np


import steervec as sv
import coordinate as coord
import beamShape as bs
import createBeam as cb


from PyQt4.QtGui import *
from PyQt4.QtCore import *

class InterferometryObservation:
    class Baseline:
        def __init__(longitude1, altitude1, longitude2, altitude2):
            self.longitude1 = longitude1
            self.altitude1 = altitude1
            self.longitude2 = longitude2
            self.altitude2 = altitude2


    def __init__(self, arrayRefereceGEODET, observeTime, waveLength):
        self.arrayRefereceGEODET = arrayRefereceGEODET
        self.observeTime = observeTime
        self.waveLength = waveLength
        self.baselines = []
        self.boreSightHorizontal = []
        self.lastBeamBoreSight = ()
        self.beamCoordinates = []

    def setObserveTime(self, dateTime):
        self.observeTime = dateTime

    def getBaselines(self):
        return self.baselines

    def getHorizontal(self):
        return self.boreSightHorizontal

    def createContour(self, antennacoor, beamBoreSight):

        antennasCoordinateFile = 'antennacoor'
        beamCoordinateFile = 'inCoord'


        beamNumber = 400
        beamSize = 1.22*self.waveLength/13.5

        if self.lastBeamBoreSight != beamBoreSight:
            # print('BoreSight Changed')
            self.lastBeamBoreSight = beamBoreSight
            beamCoordinates, subBeamRadius = cb.optimizeGrid(beamNumber, beamSize/2, cb.hexagonGrid, beamBoreSight)
            self.beamCoordinates = np.array(beamCoordinates)

        # beamCoordinates = coord.readCoordinates(beamCoordinateFile)
        antCoordinatesGEODET = np.array(antennacoor)

        antCoordinatesECEF = coord.convertGodeticToECEF(antCoordinatesGEODET)

        arrayRefereceECEF = coord.convertGodeticToECEF([arrayRefereceGEODET])[0]
        antCoordinatesENU = coord.convertECEFToENU(antCoordinatesECEF, arrayRefereceECEF, arrayRefereceGEODET)

        self.baselines = bs.createBaselines(antCoordinatesENU)

        LSTDeg =  coord.calculateLocalSiderealTime(self.observeTime, self.arrayRefereceGEODET[1])

        RA = np.deg2rad(self.beamCoordinates[:,0])
        DEC = np.deg2rad(self.beamCoordinates[:,1])
        LST = np.deg2rad(LSTDeg)
        latitude = np.deg2rad(antCoordinatesGEODET[:,0])
        longitude = np.deg2rad(antCoordinatesGEODET[:,1])
        arrayRefereceLatitude = np.deg2rad(arrayRefereceGEODET[0])

        altitude, azimuth = coord.convertEquatorialToHorizontal(
                RA, DEC, LST, arrayRefereceLatitude)

        self.boreSightHorizontal = (azimuth[0], altitude[0])

        waveNumbers = sv.waveNumber(altitude, azimuth, self.waveLength)

        weights = sv.weightVector(waveNumbers, self.baselines)

        bs.fringePlot(self.beamCoordinates, weights, self.baselines, beamBoreSight, beamSize)

def pixelCorrdinateConv(value, direction):
    if direction == 'toCoord':
        x = (value[0]/300.)*(21.431 - 21.391) + 21.391
        y = ((300.0 - value[1])/300.)*(-30.701 - -30.741) + -30.741
    elif direction == 'toPixel':
        x = int(round((value[0] - 21.391)/(21.431 - 21.391)*300.))
        y = 300.0 - int(round((value[1] - -30.741)/(-30.701 - -30.741)*300.))

    return [x,y]

class Cartesian(QWidget):

    painter = None
    dots = []


    def paintEvent(self, event):
        self.setMouseTracking(True)

        painter = QPainter()
        self.painter = painter
        painter.begin(self)
        painter.drawRect(0, 0, 300, 300)
        painter.drawLine(0, 150, 300, 150)
        painter.drawLine(150, 300, 150, 0)

        for dot in self.dots:
            painter.drawEllipse(dot[0], dot[1],3,3)

        painter.end()

    def sizeHint(self):
        return QSize(301, 301)

    def mouseMoveEvent(self, event):
        longitude, altitude = pixelCorrdinateConv((event.x(), event.y()), 'toCoord')
        longitudeCoord.setPlaceholderText('{:6.4f}'.format(longitude))
        altitudeCoord.setPlaceholderText('{:6.4f}'.format(altitude))

    def addDots(self, dots):
        self.dots += dots
        self.update()

    def mouseReleaseEvent(self, event):
        self.dots.append([event.x(), event.y()])
        longitude, altitude = pixelCorrdinateConv((event.x(), event.y()), 'toCoord')
        addRowToCoordinateList(longitude, altitude)
        self.update()
        updateCountour()
        # updateBaselineList(observation.getBaselines())

    def clearDots(self):
        self.dots = []
        self.update()

    def removeDots(self, dotsToRemove):
        for dot in dotsToRemove:
            self.dots.remove(dot)
        self.update()

def updateHorizontal(horizontalCoord):
    if horizontalCoord == []: return
    azimuthCoord.setText('{:6.4f}'.format(np.rad2deg(horizontalCoord[0])))
    elevationCoord.setText('{:6.4f}'.format(np.rad2deg(horizontalCoord[1])))

def onBoreSightUpdated():
    updateCountour()
    updateHorizontal(observation.getHorizontal())

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
        coordinates.append([longitude,altitude, 0])

    beamBoreSight = (float(RACoord.text()), float(DECCoord.text()))
    observation.createContour(coordinates, beamBoreSight)
    pixmap = QPixmap(os.getcwd() + '/contour.png')
    label.setPixmap(pixmap.scaledToHeight(pixmap.height()*0.7))
    updateBaselineList(observation.getBaselines())
    updateHorizontal(observation.getHorizontal())

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
    dot = pixelCorrdinateConv((float(longitude), float(altitude)), 'toPixel')
    axis.addDots([dot,])
    updateCountour()

def onClickedDelAllButton():
    coordinateList.setRowCount(0)
    baselineList.setRowCount(0)
    label.setPixmap(blankImage)
    axis.clearDots()

def onClickedAtCoordinateList(row, column):
    if column == 3:
        x = float(str(coordinateList.item(row, 0).text()))
        y = float(str(coordinateList.item(row, 1).text()))
        pixels = pixelCorrdinateConv((x,y), 'toPixel')
        coordinateList.removeRow(row)
        axis.removeDots([pixels,])
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
arrayRefereceGEODET = (-30.721, 21.411, 0)
'''observation time in UTC'''
observationTime = QDateTime.currentDateTime().toPyDateTime()
'''observation waveLength in meter'''
waveLength = 0.21

observation = InterferometryObservation(arrayRefereceGEODET,
        observationTime, waveLength)

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

addGeoButton = QPushButton('Delete All', w)
addGeoButton.clicked.connect(onClickedDelAllButton)
addGeoButton.move(320, 400)



dateTimeLabel = QLabel(w)
dateTimeLabel.setText('UTC Time')
dateTimeLabel.move(500, 320)
dateTimeEdit = QDateTimeEdit(w)
dateTimeEdit.move(500, 340)
dateTimeEdit.setDisplayFormat("dd.MM.yyyy hh:mm:ss.zzz")
dateTimeEdit.setDateTime(QDateTime.currentDateTime())
dateTimeEdit.dateTimeChanged.connect(onDateTimeChanged)

RACoordLabel = QLabel(w)
RACoordLabel.setText('RA')
RACoordLabel.move(500, 380)
DECCoordLabel = QLabel(w)
DECCoordLabel.setText('DEC')
DECCoordLabel.move(590, 380)


RACoord = QLineEdit(w)
RACoord.move(500, 400)
RACoord.resize(80,30)
RACoord.setText('21.411')
RACoord.editingFinished.connect(onBoreSightUpdated)
DECCoord = QLineEdit(w)
DECCoord.resize(80,30)
DECCoord.move(590, 400)
DECCoord.setText('-30.721')
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
coordinateList.setHorizontalHeaderItem(1, QTableWidgetItem('altitude'))
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
