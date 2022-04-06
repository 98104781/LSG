import inspect
from collections import Counter
import Wizard.Spectra as Spectra
import Lipids.GenerateLipids as GL
import Wizard.EditTail as ET
from itertools import combinations_with_replacement as cwr, product

from PySide6.QtCore import Signal, QModelIndex, QAbstractTableModel, QSize
from PySide6.QtGui import Qt, QCursor, QDoubleValidator, QIntValidator
from PySide6.QtWidgets import  QDialog, QPushButton, QLineEdit, QComboBox, QTableView, QVBoxLayout, QHBoxLayout, QHeaderView, QMenu

class NewWindow(QDialog):
    '''
    Window showing lipid classes, permits editing
    '''

    def __init__(self, lipid, adduct, fragtype):
        super().__init__()

        self.lipid = lipid
        self.adduct = adduct
        self.fragtype = fragtype 

        self.setWindowTitle('LSG3')
        self.setFixedSize(400, 510)
        self.vLayout = QVBoxLayout(self)
        self.hLayout = QHBoxLayout(self)
        self.hLayout2 = QHBoxLayout(self)

        self.tableView = QTableView()
        self.tableView.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.vLayout.addWidget(self.tableView)
        self.tableView.doubleClicked.connect(self.selectedFragment)

        self.buildList()

    def find_subclasses(self, module, cls):
        for name in dir(module):
            o = getattr(module, name)
            try:
                if (o != cls) and issubclass(o, cls):
                    yield name, o
            except TypeError: pass

    def buildList(self):

        totalfragList = list(self.find_subclasses(GL, GL.Fragment))

        fragList = []
        for frag in totalfragList:
            if frag[0][-1] == 'x':
                f = getattr(GL, frag[0][0:-1])
            else: f = frag[1]
            if f not in list(self.lipid.adducts[self.adduct].keys()):

                try: fgmt = f(self.lipid, self.adduct, 0)
                except: continue
                
                try: fragList.extend(fgmt)
                except AttributeError: pass # May be trying to make a headgroup fragment on lipid without headgroup
                except: fragList.append(fgmt)
                
        self.finalFragList = sorted(set(fragList), reverse=True)

        self.tableModel = ClassTableModel(self.finalFragList)
        self.tableView.setModel(self.tableModel)

    def selectedFragment(self):
        index = self.tableView.selectedIndexes()[0].row()
        self.selection = self.finalFragList[index].fragmentType
        self.done(1)
        self.close()






class ClassTableModel(QAbstractTableModel):
    def __init__(self, data=[], parent=None):
        QAbstractTableModel.__init__(self, parent)
        self.tdata = data
        self.headers = ['Fragment Types']

    def flags(self, index):
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
        frag = self.tdata[row]
        if role == Qt.DisplayRole:
            try: return f"{frag.mass}, {frag.Comment()}"
            except: print(frag)

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




