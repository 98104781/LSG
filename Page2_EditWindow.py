
import GenerateLipids as GL

from itertools import combinations_with_replacement as cwr
from PySide6.QtCharts import QChart, QChartView, QScatterSeries, QValueAxis
from PySide6.QtCore import QAbstractTableModel, QMargins, QModelIndex, Qt
from PySide6.QtGui import QColor, QImage, QPainter, QPainterPath, QPen
from PySide6.QtWidgets import QComboBox, QDialog, QHeaderView, QStyledItemDelegate, QTableView, QPushButton, QVBoxLayout, QHBoxLayout, QLabel, QSpinBox

class NewWindow(QDialog): # Opened from SpectraSetupPage
    '''
    Popup window to modify selected class/adduct spectra.
    '''
    def __init__(self, treeData):
        super().__init__()

        self.setWindowTitle('LSG3')
        self.setFixedSize(600, 510)
        self.vLayout = QVBoxLayout(self)
        self.vLayout2 = QVBoxLayout(self)
        self.hLayout = QHBoxLayout(self)

        # Example spectra will be GPL 16:0_18:1
        self.tails = [GL.sn(c=16, d=0, type='Acyl'),
                      GL.sn(c=18, d=1, type='Acyl')]

        self.tableView = QTableView()
        self.comboBox = QComboBox()
        self.spectra = SpectraScatter()
        self.spectra.setFixedHeight(200)
        self.label = QLabel("If generated spectra don't respect sn isomerism, all isomer dependent\n"
                            "fragments will be equal in intensity.\n\nIf the observed fragment intensities"
                            " differ from the default provided, they may\nbe manually updated on the left.")
        self.comboBox.currentTextChanged.connect(self.buildList) # Lipid data will be stored in comboboxs
        for lipidClass, adductList in treeData.items():  # Unpack lipid/adduct data from Treeview for use
            for adduct in adductList: # Unpack all adducts from adduct list for use in separate comboboxs
                self.comboBox.addItem(lipidClass.text(0)+' '+adduct.text(0), [lipidClass, adduct]) # here

        self.applyBTN = QPushButton('Apply Changes')
        self.applyBTN.clicked.connect(self.applyChanges)
        #self.resetBTN = QPushButton('Reset Changes')
        #self.applyBTN.clicked.connect(self.resetChanges)

        self.vLayout.addWidget(self.comboBox)
        self.hLayout.addWidget(self.tableView)

        self.vLayout2.addWidget(self.spectra)
        self.vLayout2.addWidget(self.label)
        self.vLayout2.addStretch()
        self.vLayout2.addWidget(self.applyBTN)
        self.hLayout.addLayout(self.vLayout2)

        self.vLayout.addLayout(self.hLayout)

    def buildLipid(self, data):
        '''
        Returns an example lipid for use.
        Lipid will be GPL (16:0_) 16:0_18:1.
        '''
        GL.Glycerolipid.instances = []
        cls = data[0].lipidClass # Example lipid
        _, comb, *_ = cwr(self.tails, cls.No_Tails)
        example = cls(*comb) # comb will be (16:0_) 16:0_18:1
        example.resolve_spectra(data[1].text(0), data[1].fragmentList)
        return example 
    
    def buildList(self):
        '''
        Updates selected data,
        Refreshes fragment list and spectra.
        '''
        data = self.comboBox.currentData()
        if data is not None:
            example = self.buildLipid(data)
            adduct = data[1].text(0)
            fragments = example.spectra[adduct]
            mz = GL.MA(example, adduct, 0).mass

            self.spectra.setSpectra(example.name+' '+adduct, mz, fragments)
            self.tableData = SpectraTableModel(fragments)
            self.tableView.setModel(self.tableData)
            self.tableView.setItemDelegateForColumn(1, SpinBoxDelegate(self.tableView))
            self.tableView.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
            self.tableView.verticalHeader().hide()
            self.tableData.dataChanged.connect(self.updateData)
        else: pass

    def updateData(self, index: QModelIndex, *other):
        '''
        Updates adduct spectra in list
        with new intensity.
        '''
        row = int(index.row())
        fragmentList = self.comboBox.currentData()[1].fragmentList
        data = self.tableData.table_data[row]
        fragmentList[data.fragmentType] = data.intensity
        self.comboBox.currentData()[1].fragmentList = fragmentList
        self.buildList()

    def applyChanges(self):
        '''
        Overwrites adduct spectra in memory
        with new adduct spectra.
        '''
        for x in range(self.comboBox.count()):
            data = self.comboBox.itemData(x)
            lipidClass = data[0].lipidClass
            adduct = data[1].text(0)
            fragmentList = data[1].fragmentList
            for fragment in list(fragmentList):
                if fragmentList[fragment] == 0:
                    fragmentList.pop(fragment)
            lipidClass.adducts[adduct] = fragmentList
        self.buildList()

    def resetChanges(self):
        pass

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
        self.spectra.setMarkerSize(10000)
        self.spectra.setMarkerShape(self.spectra.MarkerShapeRectangle)
        self.chart().setMargins(QMargins(-20, 0, -10, -10))
        self.chart().addSeries(self.spectra)
        self.spectra.attachAxis(self.xaxis)
        self.spectra.attachAxis(self.yaxis)
        self.chart().legend().hide()

    def setSpectra(self, precursorName, precursorMass, spectra):
        '''
        Updates displayed spectra with new spectra.
        '''
        self.chart().setTitle(precursorName)
        binMax = int(precursorMass+100)

        self.spectra.clear()

        linePath = QPainterPath()
        linePath.moveTo(5000, 5001)
        linePath.lineTo(5001, 10000)
        linePath.closeSubpath()
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
        self.xaxis.setRange(100, binMax)

class SpectraTableModel(QAbstractTableModel):
    '''
    Custom table type to organise and display
    lipid spectra data for inspection/modification
    using SpectraEditWindow dialogue.
    '''
    def __init__(self, data=[], parent=None):
        QAbstractTableModel.__init__(self, parent)
        self.table_data = data
        self.headers = ['m/z (Da)', 'Abn. (%)']
    
    def flags(self, index):
        '''
        Forbids editing of mz value displayed in column 0.
        Allows editing of abundance value displayed in column 1.
        '''
        if index.column() > 0:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable
        else:
            return Qt.ItemIsSelectable

    def rowCount(self, parent):
        return len(self.table_data)

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
                return str(self.table_data[row].mass)
            elif column == 1:
                return str(self.table_data[row].intensity)

    def setData(self, index, value, role=Qt.EditRole):
        if role == Qt.EditRole:
            row = index.row()
            column = index.column()
            self.table_data[row].intensity = value
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