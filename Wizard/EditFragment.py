import re
import copy
from types import MethodType
from collections import Counter
import Lipids.GenerateLipids as GL

from PySide6.QtGui import Qt, QRegularExpressionValidator
from PySide6.QtCore import QAbstractTableModel, QRegularExpression
from PySide6.QtWidgets import  QDialog, QTableView, QVBoxLayout, QHBoxLayout, QHeaderView, QLineEdit, QPushButton

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
        self.tableView.doubleClicked.connect(self.fragmentSelected)

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
            
            # methods share name of class with suffix x. If method, remove x, call class.
            f = getattr(GL, frag[0][0:-1]) if frag[0][-1] == 'x' else frag[1]
            
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

                except IndexError: pass # Might be missing a tail
                except AttributeError: pass # Failed to make fragment due to missing part. E.G. no headgroup.
                except AssertionError: pass # Assertion failed. E.G. Lipid formula lacks atom from fragment.
                
        self.finalFragList = sorted(set(fragList), reverse=True)
        self.tableModel = ClassTableModel(self.finalFragList)
        self.tableView.setModel(self.tableModel)

    def fragmentSelected(self):
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
        frag = self.tdata[index.row()]
        if role == Qt.DisplayRole:
            try: return f"{frag.mass}, {frag.Comment()}"
            except: pass

    def setData(self, index, value, role=Qt.EditRole):
        if role == Qt.EditRole: return True
        return QAbstractTableModel.setData(self, index, value, role)

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

class CustomisedFragment(QDialog):
    '''
    Window showing lipid classes, permits editing
    '''

    def __init__(self, lipid, adduct):
        super().__init__()

        self.lipid = lipid
        self.adduct = adduct

        self.fragment = GL.Fragment(self.lipid, self.adduct, 0)

        self.setWindowTitle('LSG3 - Custom Fragments')
        self.setFixedSize(400, 125)
        self.vLayout = QVBoxLayout(self)
        self.hLayout = QHBoxLayout(self)
        self.hLayout2 = QHBoxLayout(self)

        self.massBox = QLineEdit()
        self.massBox.setDisabled(True)
        self.massBox.setToolTip('Fragment M/Z')
        self.formulaBox = QLineEdit()
        self.formulaBox.setDisabled(True)
        self.formulaBox.setToolTip('Fragment Formula')
        self.hLayout.addWidget(self.massBox)
        self.hLayout.addWidget(self.formulaBox)
        self.vLayout.addLayout(self.hLayout)

        self.fragmentString = QLineEdit()
        self.fragmentString.setToolTip('Special Terms:\n'
                                       'M: Lipid\n'
                                       'A: Adduct\n'
                                       'sn1 - sn3')
        self.chargeString = QLineEdit('1')
        self.chargeString.setToolTip('Fragment Charge')
        self.chargeString.setFixedWidth(25)
        self.rEx = QRegularExpression('[1-9]')
        self.chargeString.setAlignment(Qt.AlignCenter)
        self.chargeString.setValidator(QRegularExpressionValidator(self.rEx))
        self.fragmentString.textChanged.connect(self.updateFragmentString)
        self.chargeString.textChanged.connect(self.updateFragmentString)
        self.hLayout2.addWidget(self.fragmentString)
        self.hLayout2.addWidget(self.chargeString)
        self.vLayout.addLayout(self.hLayout2)

        self.acceptButton = QPushButton('Accept Fragment')
        self.acceptButton.clicked.connect(self.acceptFragment)
        #self.acceptButton.setDisabled(True)
        self.vLayout.addWidget(self.acceptButton)

    def acceptFragment(self):
        if not self.massBox.text() or not self.formulaBox.text():
            pass
        else:
            self.lipid.adducts[self.adduct][self.fragment] = 0
            self.done(1)
            self.close()

    def updateBoxes(self):
        self.massBox.setText(str(self.fragment.MZ()))
        self.formulaBox.setText(''.join(''.join((key, str(val))) for (key, val ) in self.fragment.Formula().items() if val !=0))
        #if self.fragment.Formula() <= self.lipid.formula:
        #    self.acceptButton.setDisabled(False)
        #else:self.acceptButton.setDisabled(True)

    def updateFragmentString(self):
        
        self.charge = int(self.chargeString.text() or 1)*int(GL.adducts[self.adduct][2]/abs(GL.adducts[self.adduct][2]))
        self.comment = self.fragmentString.text()

        # 'M+H-H2SO4' -> ['M', '-H', '+H2SO4']
        fragUnits = re.findall('(.[^+-]*)', self.fragmentString.text())
        self.fragment.fragUnits = fragUnits
        self.fragment.charge = self.charge
        self.fragment.comment = f"[{self.comment}] {self.charge}"

        fragTerms = GL.TermList([GL.defineFragmentTerm(self, s) for s in fragUnits])
        fragTerms.charge = self.charge
        fragTerms.comment = f"[{self.comment}] ({self.charge})"

        self.fragment.MZ = fragTerms.returnMass
        self.fragment.Formula = fragTerms.returnFormula
        self.fragment.Charge = fragTerms.returnCharge
        self.fragment.Comment = fragTerms.returnComment

        self.updateBoxes()


        
