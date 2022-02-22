import inspect

import Classes
import Classes_isomers
import SpectraEditWindow as SEW

from PySide6.QtCore import Property, Qt, Signal, QAbstractTableModel
from PySide6.QtGui import QPixmap
from PySide6.QtWidgets import QPushButton, QTableView, QVBoxLayout, QHBoxLayout, QWizard, QWizardPage, QHeaderView

# This gives me the classes, but not in the order that I set them... Really irks me !
gplClassList = [cls for _, cls in inspect.getmembers(Classes) if inspect.isclass(cls)]
gplClassList_Isomers = [cls for _, cls in inspect.getmembers(Classes_isomers) if inspect.isclass(cls)]

class Page(QWizardPage):
    '''
    Alternate second page. Used for defining specific
    lipids to generate.
    '''

    def __init__(self, parent=None):
        super(Page, self).__init__(parent)

        self.setTitle("Create specific lipids")
        self.setSubTitle("Incomplete ")
        self.setPixmap(QWizard.WatermarkPixmap, QPixmap('Images\GPLs.png'))
        self.vLayout = QVBoxLayout(self)
        self.hLayout = QHBoxLayout(self)

        self.tableView = QTableView()
        self.tableView.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.tableData = LipidTableModel()
        self.tableView.setModel(self.tableData)
        self.vLayout.addWidget(self.tableView)

        self.addLipid = QPushButton('New Lipid')
        self.addLipid.clicked.connect(self.addRow)
        self.hLayout.addWidget(self.addLipid)
        self.editLipid = QPushButton('Edit Selected')
        self.hLayout.addWidget(self.editLipid)
        self.removeLipid = QPushButton('Remove Selected')
        self.hLayout.addWidget(self.removeLipid)
        self.vLayout.addLayout(self.hLayout)

    def addRow(self):
        rows = self.tableData.rowCount(self)
        self.tableData.setRowCount(rows+1)

    def editRow(self, table, row):
        pass

    def removeRow(self, table, row):
        pass
    
    def isComplete(self):

        if self.tableData.rowCount(self) < 1:
            return False

        return super().isComplete()




    def initializePage(self) -> None:

        if self.field('isomerism') == False:
            classes_to_generate = gplClassList
        else:
            classes_to_generate = gplClassList_Isomers

        return super().initializePage()

    def nextId(self):
        return 3 # Page 3



class LipidTableModel(QAbstractTableModel):
    '''
    Custom table type to organise and display
    lipid data for inspection/modification.
    '''
    def __init__(self, data=[], parent=None):
        QAbstractTableModel.__init__(self, parent)
        self.table_data = data
        self.headers = ['Lipid']
    
    def flags(self, index):
        '''
        Forbids editing of value displayed in column 0.
        '''
        if index.column() == 0:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable

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
                return str('Test')
 
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