from PySide6.QtWidgets import  QStyledItemDelegate, QSpinBox, QLabel, QMenu, QApplication
from PySide6.QtCore import QAbstractTableModel, QMargins, Qt, QByteArray, QPointF, QRectF
from PySide6.QtGui import QColor, QImage, QPainter, QPainterPath, QPen, QPixmap, QAction
from PySide6.QtCharts import QChart, QChartView, QScatterSeries, QValueAxis, QLineSeries

import numpy as np
from scipy.optimize import fsolve
from rdkit.Chem.Draw import rdMolDraw2D
import base64
import math

import Wizard.Draw as dM

class SpectraScatter(QChartView):
    '''
    Scatterplot with custom marker,
    made to appear as mass spectra.
    '''
    def __init__(self, parent=None, drawStructure=True):
        super().__init__(QChart(), parent=parent)

        self.drawStructure = drawStructure
        self.curveDisplayed = False

        self.setStyleSheet("""
            QChartView {
            background-color: white;
            border: 1px solid #DDDDDD;
            border-radius: 1px;}         
            """)

        self.xaxis = QValueAxis()
        self.xaxis.setGridLineVisible(False)
        self.chart().addAxis(self.xaxis, Qt.AlignBottom)
        self.xaxis.setTickCount(2)
        self.xaxis.setLabelsVisible(False)

        self.yaxis = QValueAxis()
        self.yaxis.setGridLineVisible(False)
        self.yaxis.setRange(0, 100) 
        self.yaxis.setMin(0)
        self.chart().addAxis(self.yaxis, Qt.AlignLeft)
        self.yaxis.setLabelsAngle(-90)
        self.yaxis.setLabelsVisible(True)    
        self.yaxis.setLabelFormat("%.0f")
        self.yaxis.setTickCount(3)

        self.spectra = QScatterSeries()
        self.spectra.setMarkerSize(10000) # Big marker to look like peak
        self.spectra.setMarkerShape(self.spectra.MarkerShapeRectangle)
        self.chart().setMargins(QMargins(-5, 0, -5, -15))
        self.chart().addSeries(self.spectra)
        self.spectra.attachAxis(self.xaxis)
        self.spectra.attachAxis(self.yaxis)
        self.chart().legend().hide()

        self.curve = QLineSeries()
        self.chart().addSeries(self.curve)
        self.curve.attachAxis(self.xaxis)
        self.curve.attachAxis(self.yaxis)

        self.tooltip = QLabel(self)
        self.tooltip.setAutoFillBackground(True)
        self.tooltip.setFrameShape(QLabel.Shape.StyledPanel)
        self.tooltip.setFrameShadow(QLabel.Shadow.Raised)
        self.tooltip.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.tooltip.hide()

        if self.drawStructure:
            self.label = QLabel(self)
            self.label.setVisible(False)
            self.label.mouseDoubleClickEvent = self.mouseDoubleClickEvent
            self.label.setContextMenuPolicy(Qt.CustomContextMenu)
            self.label.customContextMenuRequested.connect(self.showContextMenu)

        self.setRubberBand(QChartView.HorizontalRubberBand)
        self.chart().axes(Qt.Horizontal)[0].rangeChanged.connect(self.updateYaxisRange)
        self.setRenderHint(QPainter.Antialiasing)
        self.setMouseTracking(True)

    def updateYaxisRange(self, minX, maxX):

        maxY = 0
        for peak in self.spectra.points():
            if minX <= peak.x() <= maxX:
                maxY = min(math.ceil(max(maxY, 1.5*peak.y())), 100)
        if maxY < 50:
            self.yaxis.setLabelFormat("%.1f")
        else: self.yaxis.setLabelFormat("%.0f")
        self.chart().axes(Qt.Vertical)[0].setRange(0, maxY)


    def mouseMoveEvent(self, event):
        super().mouseMoveEvent(event)
        mouseLocation = self.chart().mapToValue(event.pos())
        minDist = float('inf')
        chosenPeak = QPointF()

        for peak in self.spectra.points():
            distX = abs(peak.x() - mouseLocation.x())
            distY = abs(peak.y() - mouseLocation.y())
            if distX < minDist and distY < 0.1*(self.yaxis.max() - self.yaxis.min()):
                minDist = distX
                chosenPeak = peak

        if minDist < 5 and chosenPeak.y() > 0:
            self.tooltip.setText(f"{chosenPeak.x()}")
            tooltipPos = self.chart().mapToPosition(chosenPeak)
            tooltipPos.setX(tooltipPos.x() - self.tooltip.width() / 2)
            tooltipPos.setY(tooltipPos.y() - self.tooltip.height())
            self.tooltip.move(tooltipPos.toPoint())
            self.tooltip.show()
        else:
            self.tooltip.hide()

    def leaveEvent(self, event):
        self.tooltip.hide()
        super().leaveEvent(event)

    def mousePressEvent(self, event):
        self.tooltip.hide()
        super().mousePressEvent(event)
        if event.button() == Qt.RightButton:
            self.chart().zoomReset()

    def mouseReleaseEvent(self, event):
        super().mouseReleaseEvent(event)
        if event.button() == Qt.RightButton:
            self.chart().zoomReset()

    def mouseDoubleClickEvent(self, event):
        if event.button() == Qt.LeftButton and self.drawStructure:
            self.label.setVisible(not self.label.isVisible())

    def showContextMenu(self, pos):
        contextMenu = QMenu()
        copyAction = QAction("Copy Image", self)
        copyAction.triggered.connect(self.copyImage)
        contextMenu.addAction(copyAction)
        contextMenu.exec_(self.mapToGlobal(pos))

    def copyImage(self):
        if self.pixmap:
            clipboard = QApplication.clipboard()
            image = self.pixmap.toImage()
            clipboard.setImage(image)

    def resizeEvent(self, event):
        if self.drawStructure:
            self.label.resize(self.size().width(), self.size().height())
            pix = dM.drawMolecule(self.smiles, self.size().width()-2, self.size().height()-2)
            self.pixmap = dM.framePixmap(pix)
            self.label.setPixmap(self.pixmap)
        super().resizeEvent(event)

    def paintEvent(self, event):
        super().paintEvent(event)


    def setSpectra(self, precursorName, precursorMass, spectra, smiles = 'OCC(O)CO', isotopeSpectra=False):
        '''
        Updates displayed spectra with new spectra.
        '''
        self.precursorName = precursorName 
        self.chart().setTitle(self.precursorName)
        
        if isotopeSpectra:
            masses = [peak[0] for peak in spectra]
            self.binMin = min(masses)-1
            self.binMax = max(masses)+1
        else:
            try:
                self.binMax = int(spectra[0].mass+100)
            except: self.binMax = int(precursorMass+100)
            self.binMin = 0

        self.spectra.clear()

        linePath = QPainterPath()
        linePath.moveTo(5000, 5001) #  Makes a big square
        linePath.lineTo(5000, 10000) #  draws a line from
        linePath.closeSubpath() # middle to bottom. 'peak'
        image = QImage(10000, 10000, QImage.Format_ARGB32)
        painter = QPainter(image)
        pen = QPen(QColor(0, 0, 0), 1, Qt.SolidLine)
        painter.setPen(pen)
        painter.setBrush(painter.pen().color())
        painter.drawPath(linePath)
        painter.end()

        self.spectra.setBrush(image)
        self.spectra.setPen(QColor(Qt.transparent))

        for peak in spectra:
            try:
                x = peak.mass
                y = peak.intensity
            except:
                x = peak[0]
                y = peak[1] 
            self.spectra.append(x, y)

        self.xaxis.setRange(self.binMin, self.binMax)
        self.xaxis.setMin(self.binMin)
        self.xaxis.setMax(self.binMax)

        self.smiles = smiles
        self.resizeEvent(None) # <- redraws molecule
        #self.pixmap = dM.drawMolecule(self.smiles, self.size().width(), self.size().height())

    def drawBellCurves(self, spectra, res=5000):

        self.curve.clear()
        self.curveDisplayed = False

        pen = QPen(QColor(0, 0, 0), 1, Qt.SolidLine)
        self.curve.setPen(pen)

        sigma = spectra[0][0]/(2.355*res)

        xVals = np.linspace(math.floor(self.binMin), math.ceil(self.binMax), int((math.ceil(self.binMax) - math.floor(self.binMin)) * int(res/50)))
        yVals = np.zeros_like(xVals)

        for peak in spectra:
            yVals += np.round(peak[1] * np.exp(-(xVals - peak[0])**2 / (2 * sigma**2)), 3)

        self.curve.append(self.binMin, 0.001)
        for x, y in zip(xVals, yVals):
            if y > 0: self.curve.append(x, y)
        self.curve.append(self.binMax, 0.001)
        self.curveDisplayed = True





class SpectraTableModel(QAbstractTableModel):
    '''
    Custom table type to organise and display
    lipid spectra data for inspection/modification
    using SpectraEditWindow dialogue.
    '''
    def __init__(self, data=[], parent=None):
        QAbstractTableModel.__init__(self, parent)
        self.tdata = data
        self.headers = ['m/z (Da)', 'Abn. (%)']
    
    def flags(self, index):
        '''
        Forbids editing of mz value displayed in column 0.
        Allows editing of abundance value displayed in column 1.
        '''
        if index.column() == 1:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable
        else:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable 

    def rowCount(self, parent):
        return len(self.tdata)

    def columnCount(self, parent):
        return len(self.headers)
    
    def headerData(self, column, orientation, role=Qt.DisplayRole):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return self.headers[column]
    
    def data(self, index, role=Qt.DisplayRole):
        row = index.row()
        column = index.column()
        if role == Qt.DisplayRole:
            if column == 0:
                try:
                    val = str(self.tdata[row].mass)
                except: val = self.tdata[row][0]
                return val
            elif column == 1:
                try:
                    val = str(self.tdata[row].intensity)
                except: val = self.tdata[row][1]                
                return val
        elif role == Qt.TextAlignmentRole:
            return Qt.AlignCenter
        elif role == Qt.ToolTipRole:
            return str(self.tdata[row].Comment())

    def setData(self, index, value, role=Qt.EditRole):
        if role == Qt.EditRole:
            row = index.row()
            column = index.column()
            self.tdata[row].intensity = value
            self.dataChanged.emit(index, index)
            return True
        return QAbstractTableModel.setData(self, index, value, role)

    def resizeEvent(self, event):
        pass

    def insertRow(self):
        pass

    def removeRow(self, position):
        pass


class SpinBoxDelegate(QStyledItemDelegate):
    def __init__(self, parent=None):
        QStyledItemDelegate.__init__(self, parent)

    def createEditor(self, parent, option, index):
        editor = QSpinBox(parent)
        editor.setMinimum(0)
        editor.setMaximum(100)
        editor.setSingleStep(1)
        return editor
            
    def setEditorData(self, spinBox, index):
        value = int(index.model().data(index, Qt.DisplayRole))
        spinBox.setValue(value)

    def setModelData(self, spinBox, model, index):
        value = spinBox.value()
        model.setData(index, value)

    def updateEditorGeometry(self, editor, option, index):
        editor.setGeometry(option.rect)