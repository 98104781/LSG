from PySide6.QtWidgets import  QStyledItemDelegate, QSpinBox, QLabel, QMenu, QApplication
from PySide6.QtCore import QAbstractTableModel, QMargins, Qt, QByteArray, QPointF
from PySide6.QtGui import QColor, QImage, QPainter, QPainterPath, QPen, QPixmap, QAction
from PySide6.QtCharts import QChart, QChartView, QScatterSeries, QValueAxis

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
import base64

import Wizard.Draw as dM

class SpectraScatter(QChartView):
    '''
    Scatterplot with custom marker,
    made to appear as mass spectra.
    '''
    def __init__(self, parent=None):
        super().__init__(QChart(), parent=parent)

        self.xaxis = QValueAxis()
        self.xaxis.setGridLineVisible(False)
        self.chart().addAxis(self.xaxis, Qt.AlignBottom)
        self.xaxis.setTickCount(2)
        self.xaxis.setLabelsVisible(False)

        self.yaxis = QValueAxis()
        self.yaxis.setGridLineVisible(False)
        self.chart().addAxis(self.yaxis, Qt.AlignLeft)
        self.yaxis.setLabelsVisible(False)    

        self.spectra = QScatterSeries()
        self.spectra.setMarkerSize(10000) # Big marker to look like peak
        self.spectra.setMarkerShape(self.spectra.MarkerShapeRectangle)
        self.chart().setMargins(QMargins(-20, 0, -10, -10))
        self.chart().addSeries(self.spectra)
        self.spectra.attachAxis(self.xaxis)
        self.spectra.attachAxis(self.yaxis)
        self.chart().legend().hide()

        self.tooltip = QLabel(self)
        self.tooltip.setAutoFillBackground(True)
        self.tooltip.setFrameShape(QLabel.Shape.StyledPanel)
        self.tooltip.setFrameShadow(QLabel.Shadow.Raised)
        self.tooltip.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.tooltip.hide()

        self.label = QLabel(self)
        self.label.setVisible(False)
        self.label.mousePressEvent = self.mousePressEvent
        self.label.setContextMenuPolicy(Qt.CustomContextMenu)
        self.label.customContextMenuRequested.connect(self.showContextMenu)

    def mouseMoveEvent(self, event):
        mouseLocation = self.chart().mapToValue(event.pos())
        minDist = float('inf')
        chosenPeak = QPointF()

        for peak in self.spectra.points():
            dist = QPointF(peak.x() - mouseLocation.x(), 0).manhattanLength()
            if dist < minDist:
                minDist = dist
                chosenPeak = peak

        if minDist < 5 and chosenPeak.y() > 0:
            self.tooltip.setText(f"{chosenPeak.x()}")
            tooltipPos = self.chart().mapToPosition(chosenPeak)
            #tooltipPos.setY(0.8 * tooltipPos.y())
            tooltipPos.setX(tooltipPos.x() - self.tooltip.width() / 2)
            tooltipPos.setY(tooltipPos.y() - self.tooltip.height())
            self.tooltip.move(tooltipPos.toPoint())
            self.tooltip.show()
        else:
            self.tooltip.hide()

    def mousePressEvent(self, event):
        if event.button() == Qt.LeftButton:
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
        self.label.resize(self.size().width(), self.size().height())
        pix = dM.drawMolecule(self.smiles, self.size().width()-2, self.size().height()-2)
        self.pixmap = dM.framePixmap(pix)
        self.label.setPixmap(self.pixmap)
        super().resizeEvent(event)

    def setSpectra(self, precursorName, precursorMass, spectra, smiles = 'OCC(O)CO'):
        '''
        Updates displayed spectra with new spectra.
        '''
        self.precursorName = precursorName 
        self.chart().setTitle(self.precursorName)
        try:
            binMax = int(spectra[0].mass+100)
        except: binMax = int(precursorMass+100)

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
            x = peak.mass
            y = peak.intensity
            self.spectra.append(x, y)

        self.yaxis.setRange(0, 100)
        self.xaxis.setRange(0, binMax)

        self.smiles = smiles
        self.resizeEvent(None) # <- redraws molecule
        #self.pixmap = dM.drawMolecule(self.smiles, self.size().width(), self.size().height())

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
            return Qt.ItemIsSelectable

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
                return str(self.tdata[row].mass)
            elif column == 1:
                return str(self.tdata[row].intensity)
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