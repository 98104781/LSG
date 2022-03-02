import Page2B_EditWindow as P2BEW

from PySide6.QtGui import QPixmap
from PySide6.QtCore import Qt, QAbstractTableModel, Property, Signal
from PySide6.QtWidgets import QPushButton, QTableView, QVBoxLayout, QHBoxLayout, QWizard, QWizardPage, QHeaderView

class Page(QWizardPage):
    '''
    Alternate second page. Used for defining specific
    lipids to generate.
    '''

    def __init__(self, parent=None):
        super(Page, self).__init__(parent)

        self.setTitle("Create specific lipids")
        self.setSubTitle(" ")
        self.setPixmap(QWizard.WatermarkPixmap, QPixmap('Images\GPLs.png'))
        self.vLayout = QVBoxLayout(self)
        self.hLayout = QHBoxLayout(self)
        self.setCommitPage(True)

        self.lipidList = []
        self.tableView = QTableView()
        self.tableView.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.buildList()
        self.vLayout.addWidget(self.tableView)

        self.addLipid = QPushButton('New Lipid')
        self.addLipid.clicked.connect(self.addNewLipid)
        self.hLayout.addWidget(self.addLipid)
        self.editLipid = QPushButton('Edit Selected')
        self.editLipid.clicked.connect(self.editSelectedLipid)
        self.hLayout.addWidget(self.editLipid)
        self.removeLipid = QPushButton('Remove Selected')
        self.removeLipid.clicked.connect(self.removeSelectedLipid)
        self.hLayout.addWidget(self.removeLipid)
        self.vLayout.addLayout(self.hLayout)

        self.registerField("lipidList", self, "tableProperty")

    def buildList(self):
        self.tableModel = LipidTableModel(self.lipidList)
        self.tableView.setModel(self.tableModel)
        self.completeChanged.emit()

    def openLipidEditor(self, selection=None):
        editspectrawindow = P2BEW.NewWindow(self, selection)
        if editspectrawindow.exec() > 0:
            lipid = editspectrawindow.lipid
            adduct = editspectrawindow.lipidAdduct.currentData()
            return lipid, adduct

    def addNewLipid(self):
        try: 
            lipid, adduct = self.openLipidEditor()
            self.lipidList.append([lipid, adduct])
        except: pass
        self.buildList()

    def editSelectedLipid(self):
        try: 
            row = self.tableView.selectedIndexes()[0].row()
            lipid, adduct = self.openLipidEditor(self.lipidList[row])
            if lipid: self.lipidList[row] = [lipid, adduct]
        except: pass
        self.buildList()

    def removeSelectedLipid(self):
        try:
            row = self.tableView.selectedIndexes()[0].row()
            self.lipidList.pop(row)
        except: pass
        self.buildList()

    def isComplete(self):
        if self.tableModel.rowCount(self) < 1:
            return False
        return super().isComplete()

    def setLipidList(self, data):
        self.lipidList = data

    def getLipidList(self):
        return self.lipidList

    def nextId(self):
        return 3 # Page 3

    tableProperty = Property(list, getLipidList, setLipidList)
    lipidListChanged = Signal()


class LipidTableModel(QAbstractTableModel):
    '''
    Custom table type to organise and display
    lipid data for inspection/modification.
    '''
    def __init__(self, data=[], parent=None):
        QAbstractTableModel.__init__(self, parent)
        self.tdata = data
        self.headers = ['Lipid']
    
    def flags(self, index):
        '''
        Forbids editing of value displayed in column 0.
        '''
        if index.column() == 0:
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
        if role == Qt.DisplayRole:
                return str(self.tdata[row][0].name+' '+
                self.tdata[row][1])
 
    def setData(self, index, value, role=Qt.EditRole):
        if role == Qt.EditRole:
            row = index.row()
            column = index.column()
            return True

        return QAbstractTableModel.setData(self, index, value, role)

    def resizeEvent(self, event):
        pass

    def insertRow(self):
        pass

    def removeRow(self, position):
        pass