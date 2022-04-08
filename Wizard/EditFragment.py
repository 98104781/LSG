import Lipids.GenerateLipids as GL

from PySide6.QtGui import Qt
from PySide6.QtCore import QAbstractTableModel
from PySide6.QtWidgets import  QDialog, QTableView, QVBoxLayout, QHBoxLayout, QHeaderView

class PredefinedFragment(QDialog):
    '''
    Window showing lipid classes, permits editing
    '''

    def __init__(self, lipid, adduct):
        super().__init__()

        self.lipid = lipid
        self.adduct = adduct

        self.setWindowTitle('LSG3 - Predefined Fragments')
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
            if f not in [key for key in self.lipid.adducts[self.adduct]]:

                try: 
                    # First try if fragment can be made
                    fgmt = f(self.lipid, self.adduct, 0)

                    try:
                        fgmt = list(fgmt) # If it's a generator, turn to a list and assert for
                        for x in fgmt:    # every fragment in the list
                            x.Validate()
                            assert ((x.Charge() > 0) - (x.Charge() < 0) # Check correct polarity
                            == (GL.adducts[self.adduct][2] > 0) - (GL.adducts[self.adduct][2] < 0))   
                            fragList.append(x) # If it's all fine, add to the list.

                    except: # If it can't be turned to a list, probably a single fragment
                        fgmt.Validate()
                        assert ((fgmt.Charge() > 0) - (fgmt.Charge() < 0) # Check correct polarity
                        == (GL.adducts[self.adduct][2] > 0) - (GL.adducts[self.adduct][2] < 0)) 
                        fragList.append(fgmt) # If it's all fine, add to the list. 

                except AttributeError: pass # Failed to make fragment due to missing part. E.G. no headgroup.
                except AssertionError: pass # Assertion failed. E.G. Lipid formula lacks atom from fragment.
                
        self.finalFragList = sorted(set(fragList), reverse=True)

        self.tableModel = ClassTableModel(self.finalFragList)
        self.tableView.setModel(self.tableModel)

    def selectedFragment(self):
        index = self.tableView.selectedIndexes()[0].row()
        self.selection = self.finalFragList[index].fragmentType
        self.done(1)
        self.close()




class CustomisedFragment(QDialog):
    '''
    Window showing lipid classes, permits editing
    '''

    def __init__(self, lipid, adduct):
        super().__init__()

        self.lipid = lipid
        self.adduct = adduct

        self.setWindowTitle('LSG3 - Custom Fragments')
        self.setFixedSize(400, 510)
        self.vLayout = QVBoxLayout(self)
        self.hLayout = QHBoxLayout(self)
        self.hLayout2 = QHBoxLayout(self)

        self.tableView = QTableView()
        self.tableView.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.vLayout.addWidget(self.tableView)




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




