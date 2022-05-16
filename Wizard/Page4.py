import Wizard.ResourcePath as RP
import Wizard.EditLipidAdduct as LAEW

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

        self.parent = parent

        self.setTitle("Create specific lipids")
        self.setSubTitle("Define a list of specific lipids to generate")
        image_Path = RP.resource_path('Images\GPLs.png')
        self.setPixmap(QWizard.WatermarkPixmap, QPixmap(image_Path))
        self.vLayout = QVBoxLayout(self)
        self.hLayout = QHBoxLayout(self)
        self.setCommitPage(True)
        self.tableView = QTableView()
        self.tableView.doubleClicked.connect(self.editSelectedLipid)
        self.tableView.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
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

        self.lipidList = []
        self.buildList()
        self.registerField("lipidList", self, "tableProperty")

    def buildList(self):
        self.tableModel = LipidTableModel(self.lipidList)
        self.tableView.setModel(self.tableModel)
        self.completeChanged.emit()

    def openLipidEditor(self, selection=None):

        if selection is None: # ie not editing a lipid
            lipidClass = self.parent.classes_to_generate
        else: lipidClass = [type(selection[0])] # else edit selection

        editspectrawindow = LAEW.LipidWindow(self, lipidClass, classEdit=False, selection=selection)
        if editspectrawindow.exec() > 0:
            lipid = editspectrawindow.lipid
            adduct = editspectrawindow.adductBox.currentText()
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

    def removeSelectedLipid(self): # Reverse list of indexes so ones towards the end are deleted first
        indexesToDelete = sorted([row.row() for row in self.tableView.selectedIndexes()], reverse=True)
        for index in indexesToDelete:
            self.lipidList.pop(index)
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
        return 5

    tableProperty = Property(list, getLipidList, setLipidList)
    lipidListChanged = Signal()

    def initializePage(self):
        self.lipidList = self.tableModel.tdata
        return super().initializePage()

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