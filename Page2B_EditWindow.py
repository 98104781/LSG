

import inspect
import Classes
import Classes_isomers
import GenerateLipids as GL

from itertools import combinations_with_replacement as cwr
from PySide6.QtGui import QIntValidator
from PySide6.QtWidgets import QDialog, QPushButton, QVBoxLayout, QHBoxLayout, QLineEdit, QComboBox

class NewWindow(QDialog):
    '''
    Popup window to add/modify selected lipid.
    '''
    def __init__(self, parent, lipidList, row):
        super().__init__()

        self.lipidList = lipidList
        self.row = row

        self.setWindowTitle('LSG3')
        self.setFixedSize(600, 510)
        self.vLayout = QVBoxLayout(self)
        self.vLayout2 = QVBoxLayout(self)
        self.hLayout = QHBoxLayout(self)
        self.hLayout2 = QHBoxLayout(self)

        self.lipidClass = QComboBox()
        self.hLayout.addWidget(self.lipidClass)
        self.lipidAdduct = QComboBox()
        self.hLayout.addWidget(self.lipidAdduct)
        self.vLayout.addLayout(self.hLayout)
        self.vLayout.addLayout(self.hLayout2)

        if parent.field('isomerism') == False:
            self.classes_to_generate = [cls for _, cls in inspect.getmembers(Classes) if inspect.isclass(cls)]
        else:
            self.classes_to_generate = [cls for _, cls in inspect.getmembers(Classes_isomers) if inspect.isclass(cls)]

        self.tails = []
        self.tailButton = {}
        self.lipidClass.currentTextChanged.connect(self.updateAdducts)
        self.lipidClass.currentTextChanged.connect(self.updateTailButtons)      
        for cls in self.classes_to_generate:
            self.lipidClass.addItem(cls.__name__, cls)

        self.acceptButton = QPushButton('Accept')
        self.acceptButton.clicked.connect(self.acceptLipid)
        self.vLayout.addWidget(self.acceptButton)

    def updateAdducts(self):
        self.lipidAdduct.clear()
        for adduct in self.lipidClass.currentData().adducts.keys():
            self.lipidAdduct.addItem(adduct, adduct)

    def updateTailButtons(self):
        for button in self.tailButton.keys():
            self.tailButton[button].setParent(None)
        self.tails = []
        self.tailButton = {}
        for button in range(self.lipidClass.currentData().No_Tails):
            self.tailButton[button] = QPushButton()
            self.tailButton[button].clicked.connect(self.open_specifyTailWindow)
            self.hLayout2.addWidget(self.tailButton[button])

    def open_specifyTailWindow(self):
        '''
        Opens external window to add/modify lipid.
        '''
        button = self.sender()
        editspectrawindow = TailWindow(button)
        editspectrawindow.exec()
        
    def acceptLipid(self):
        
        self.tails = []
        try:
            for button in self.tailButton.values():
                self.tails.append(button.tail)
            lipid = self.lipidClass.currentData()
            lipid = lipid(*self.tails)
            string = lipid.name+' '+self.lipidAdduct.currentData()
            if self.row is None:
                self.lipidList.insert(len(self.lipidList)+1, string)
            else: 
                self.lipidList[self.row] = string
            self.close()
        except: pass
    
class TailWindow(QDialog):
    '''
    Popup window to specify tail stats.
    '''
    def __init__(self, button):
        super().__init__()

        self.setWindowTitle('LSG3')
        self.setFixedSize(600, 125)
        self.vLayout = QVBoxLayout(self)
        self.hLayout = QHBoxLayout(self)
        self.button = button
    
        self.c = QLineEdit()
        self.c.setPlaceholderText('Chain Length')
        self.c.setValidator(QIntValidator(1, 100))
        self.hLayout.addWidget(self.c)
        self.d = QLineEdit()
        self.d.setPlaceholderText('Desaturation')
        self.d.setValidator(QIntValidator(1, 100))
        self.hLayout.addWidget(self.d)
        self.oh = QLineEdit()
        self.oh.setPlaceholderText('# of -OH groups')
        self.oh.setValidator(QIntValidator(1, 100))
        self.hLayout.addWidget(self.oh)
        self.dt = QLineEdit()
        self.dt.setPlaceholderText('# of deuterium labels')
        self.dt.setValidator(QIntValidator(1, 100))
        self.hLayout.addWidget(self.dt)
        self.vLayout.addLayout(self.hLayout)

        self.tail = None
        self.acceptButton = QPushButton('Accept')
        self.acceptButton.clicked.connect(self.accept)
        self.vLayout.addWidget(self.acceptButton)

    def accept(self):

        c = self.c.text()
        d = self.d.text()
        oh= self.oh.text()
        dt= self.dt.text()

        try: # Tail must be possible
            if (int(d or 0)+int(oh or 0)) > int(c or 1)/2: return # D + -OH < C
            if int(dt or 0) > (2*int(c or 1)-2*int(d or 0)-1): return # Deuterium <= Hydrogens
            self.tail = GL.sn(int(c), int(d or 0), type='Acyl', oh=int(oh or 0), dt=int(dt or 0))
            self.button.tail = self.tail
            self.button.setText(self.tail.name)
            self.close()
        except: pass
