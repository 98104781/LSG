import os
import re
import random

import Wizard.Draw as dM
import Wizard.EditTail as ET
import ResourcePath as RP
import Lipids.GenerateLipids as GL

from PySide6.QtGui import QPixmap
from PySide6.QtCore import Qt, QAbstractTableModel, Property, Signal
from PySide6.QtWidgets import QPushButton, QTableView, QVBoxLayout, QHBoxLayout, QWizard, QWizardPage, QHeaderView, QFileDialog

class Page(QWizardPage):
    '''
    Define Specific Tails
    '''

    def __init__(self, parent=None):
        super(Page, self).__init__(parent)

        self.setTitle("Generate a range of lipids using specific tails")
        self.setSubTitle("Please define a list of tails to use.\nSome common tails are preset.")
        image_Path = RP.resource_path('Images\FAs.png')
        
        #image = dM.drawMolecule(smiles='OC(=O)'+random.randint(5,18)*'C', width=395, height=130)
        #rotatedImage = dM.rotatePixmap(image, 90)
        #framedImage = dM.framePixmap(rotatedImage)

        self.setPixmap(QWizard.WatermarkPixmap, QPixmap(image_Path))

        self.vLayout = QVBoxLayout(self)
        self.hLayout = QHBoxLayout(self)
        self.hLayout2 = QHBoxLayout(self)

        self.importButton = QPushButton('Import Tail List')
        self.hLayout.addWidget(self.importButton)
        self.importButton.clicked.connect(self.importTailList)
        self.exportButton = QPushButton('Export Tail List')
        self.exportButton.clicked.connect(self.exportTailList)
        self.hLayout.addWidget(self.exportButton)
        self.vLayout.addLayout(self.hLayout)

        self.tableView = QTableView()
        self.tableView.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.vLayout.addWidget(self.tableView)

        self.addLipid = QPushButton('New Tail')
        self.addLipid.clicked.connect(self.addNewTail)
        self.hLayout2.addWidget(self.addLipid)
        self.removeLipid = QPushButton('Remove Selected')
        self.removeLipid.clicked.connect(self.removeSelectedTail)
        self.hLayout2.addWidget(self.removeLipid)
        self.vLayout.addLayout(self.hLayout2)

        self.registerField("tailList", self, "tableProperty")

        # Common fatty acids from various sources
        # https://doi.org/10.1016/B978-0-12-809521-8.00002-7
        # https://doi.org/10.1016/j.pharmthera.2021.107972
        self.tailList = [GL.sn( 2, 0, type='Acyl'), GL.sn( 3, 0, type='Acyl'), GL.sn( 4, 0, type='Acyl'), GL.sn( 5, 0, type='Acyl'),
                         GL.sn( 6, 0, type='Acyl'), GL.sn( 7, 0, type='Acyl'), GL.sn( 8, 0, type='Acyl'), GL.sn( 9, 0, type='Acyl'),
                         GL.sn(10, 0, type='Acyl'), GL.sn(11, 0, type='Acyl'), GL.sn(12, 0, type='Acyl'), GL.sn(12, 1, type='Acyl'), 
                         GL.sn(13, 0, type='Acyl'), GL.sn(14, 0, type='Acyl'), GL.sn(14, 1, type='Acyl'), GL.sn(15, 0, type='Acyl'), 
                         GL.sn(15, 1, type='Acyl'), GL.sn(16, 0, type='Acyl'), GL.sn(16, 1, type='Acyl'), GL.sn(16, 2, type='Acyl'),
                         GL.sn(17, 0, type='Acyl'), GL.sn(17, 1, type='Acyl'), GL.sn(17, 2, type='Acyl'), GL.sn(18, 0, type='Acyl'), 
                         GL.sn(18, 1, type='Acyl'), GL.sn(18, 2, type='Acyl'), GL.sn(18, 3, type='Acyl'), GL.sn(18, 4, type='Acyl'), 
                         GL.sn(19, 0, type='Acyl'), GL.sn(19, 1, type='Acyl'), GL.sn(20, 0, type='Acyl'), GL.sn(20, 1, type='Acyl'), 
                         GL.sn(20, 2, type='Acyl'), GL.sn(20, 3, type='Acyl'), GL.sn(20, 4, type='Acyl'), GL.sn(20, 5, type='Acyl'), 
                         GL.sn(21, 0, type='Acyl'), GL.sn(21, 1, type='Acyl'), GL.sn(21, 2, type='Acyl'), GL.sn(22, 0, type='Acyl'), 
                         GL.sn(22, 1, type='Acyl'), GL.sn(22, 2, type='Acyl'), GL.sn(22, 3, type='Acyl'), GL.sn(22, 4, type='Acyl'),
                         GL.sn(22, 5, type='Acyl'), GL.sn(22, 6, type='Acyl'), GL.sn(23, 0, type='Acyl'), GL.sn(24, 0, type='Acyl'), 
                         GL.sn(24, 1, type='Acyl'), GL.sn(24, 4, type='Acyl'), GL.sn(25, 0, type='Acyl'), GL.sn(26, 0, type='Acyl'),
                         GL.sn(26, 1, type='Acyl'), GL.sn(26, 2, type='Acyl')]
        self.tailList.extend([GL.sn( 16, 0, type='Ether'),      GL.sn( 18, 0, type='Ether'),      GL.sn( 20, 0, type='Ether')])   
        self.tailList.extend([GL.sn( 16, 0, type='Vinyl'),      GL.sn( 18, 0, type='Vinyl'),      GL.sn( 20, 0, type='Vinyl')])      
        self.tailList.extend([GL.sn( 12, 0, type='Acyl', oh=1), GL.sn( 14, 0, type='Acyl', oh=1), GL.sn( 16, 0, type='Acyl', oh=1),
                              GL.sn( 17, 0, type='Acyl', oh=1), GL.sn( 18, 0, type='Acyl', oh=1), GL.sn( 18, 1, type='Acyl', oh=1),
                              GL.sn( 18, 2, type='Acyl', oh=1), GL.sn( 20, 1, type='Acyl', oh=1), GL.sn( 20, 2, type='Acyl', oh=1)]) 

        self.tailList.extend([GL.base( 16, type='Sphingosine'), GL.base( 17, type='Sphingosine'), GL.base( 18, type='Sphingosine'), GL.base( 19, type='Sphingosine'), GL.base( 20, type='Sphingosine')])
        self.tailList.extend([GL.base( 17, type='Dihydrodeoxysphinganine'), GL.base( 18, type='Dihydrodeoxysphinganine'), GL.base( 19, type='Dihydrodeoxysphinganine')])
        self.tailList.extend([GL.base( 18, type='Sphinganine')])
        self.tailList.extend([GL.base( 18, type='Phytosphingosine'), GL.base( 20, type='Phytosphingosine')])
        self.tailList.extend([GL.base( 18, type='Deoxysphinganine')])
        self.tailList.extend([GL.base( 18, type='Sphingadiene')])

        self.buildList()

    def importTailList(self):
        file_name, filter = QFileDialog.getOpenFileName(filter="Tail List (*.txt)", selectedFilter='')

        if file_name and os.path.exists(file_name):
            open_file = open(file_name, 'r')
            lines = [line.strip() for line in open_file.readlines()]
            lines = [re.split(r'\s{1,}', line) for line in lines]
            for line in lines:
                type=line[2]
                if type in ['Acyl', 'Ether', 'Vinyl']:
                    try:
                        tail = GL.sn(c=int(line[0]), d=int(line[1]),
                                     type=line[2],   me=int(line[3]),
                                     oh=int(line[4]),dt=int(line[5]))
                        self.getTail(tail)
                    except: continue
                elif type in GL.baseTypes:
                    try:
                        base = GL.base(c=int(line[0]),
                                       type=line[2],
                                       dt=int(line[5]))
                        self.getTail(base)
                    except: continue        
            open_file.close()

    def exportTailList(self):
        file_name, filter = QFileDialog.getSaveFileName(filter="Tail List (*.txt)", selectedFilter='')

        if file_name:
            if os.path.exists(file_name):
                try:os.remove(file_name) # Removes if exists
                except PermissionError:
                    pass
            save_file = open(file_name, 'x', newline='')
            outputStr = '\n'.join(f"{sn.c} {sn.d} {sn.type} {sn.me} {sn.oh} {sn.dt}" for sn in self.tailList)
            save_file.write(outputStr)
            save_file.close()

    def buildList(self):
        self.tailList = sorted(self.tailList)
        self.tableModel = TailTableModel(self.tailList)
        self.tableView.setModel(self.tableModel)
        self.completeChanged.emit()

    def getTail(self, tail):
        if tail not in self.tailList:
            self.tailList.append(tail)
            self.buildList()

    def openTailEditor(self, selection=None):
        tailWindow = ET.TailWindow()
        if tailWindow.exec() > 0:
            tail = ''
            return tail

    def addNewTail(self):
        tailWindow = ET.TailWindow('X')
        tailWindow.output.connect(self.getTail)
        tailWindow.exec()

    def removeSelectedTail(self): # Reverse list of indexes so ones towards the end are deleted first
        indexesToDelete = sorted([row.row() for row in self.tableView.selectedIndexes()], reverse=True)
        for index in indexesToDelete:
            self.tailList.pop(index)
        self.buildList()

    def isComplete(self):
        if self.tableModel.rowCount(self) < 1:
            return False
        return super().isComplete()

    def setTailList(self, data):
        self.tailList = data

    def getTailList(self):
        return self.tailList

    def nextId(self):
        return 3

    tableProperty = Property(list, getTailList, setTailList)
    lipidListChanged = Signal()


class TailTableModel(QAbstractTableModel):
    '''
    Custom table type to organise and display
    lipid data for inspection/modification.
    '''
    def __init__(self, data=[], parent=None):
        QAbstractTableModel.__init__(self, parent)
        self.tdata = data
        self.headers = ['Tail']
    
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
            if self.tdata[row].type in GL.baseTypes:
                return f"{self.tdata[row].name} {self.tdata[row].type}"
            return str(self.tdata[row].name)
 
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