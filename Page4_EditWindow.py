import copy
import inspect
import Spectra
import Classes
import Classes_isomers
import GenerateLipids as GL

import TailEditWindow as TEW
from PySide6.QtCore import QModelIndex
from PySide6.QtWidgets import QDialog, QPushButton, QVBoxLayout, QHBoxLayout, QLineEdit, QComboBox, QTableView, QHeaderView

class NewWindow(QDialog):
    '''
    Popup window to add/modify selected lipid.
    '''
    def __init__(self, parent, selection):
        super().__init__()

        if parent.field('isomerism') == False:
            self.classes_to_generate = [cls for cls in GL.Glycerolipid.__subclasses__() if inspect.getmodule(cls) == Classes]
            self.classes_to_generate.extend([cls for cls in GL.Sphingolipid.__subclasses__() if inspect.getmodule(cls) == Classes])
            self.classes_to_generate.extend([cls for cls in GL.OtherLipid.__subclasses__() if inspect.getmodule(cls) == Classes])
        else:
            self.classes_to_generate = [cls for cls in GL.Glycerolipid.__subclasses__() if inspect.getmodule(cls) == Classes_isomers]

        self.setWindowTitle('LSG3')
        self.setFixedSize(600, 510)
        self.hLayout = QHBoxLayout(self)
        self.vLayout = QVBoxLayout(self)
        self.hLayout2 = QHBoxLayout(self)
        self.hLayout3 = QHBoxLayout(self)

        self.tableView = QTableView()
        self.hLayout.addWidget(self.tableView)

        self.spectra = Spectra.SpectraScatter()
        self.spectra.setFixedHeight(200)
        self.spectra.setFixedWidth(400)
        self.vLayout.addWidget(self.spectra)
        self.hLayout.addLayout(self.vLayout)

        self.lipidClass = QComboBox()
        self.lipidAdduct = QComboBox()
        self.hLayout2.addWidget(self.lipidClass)
        self.hLayout2.addWidget(self.lipidAdduct)
        self.vLayout.addLayout(self.hLayout2)
        self.vLayout.addLayout(self.hLayout3)
        self.vLayout.addStretch()

        self.tailButton = {} # Buttons for tails stored in this dict
        self.lipidClass.currentTextChanged.connect(self.updateAdducts)
        self.lipidClass.currentTextChanged.connect(self.updateTailButtons)
        for cls in self.classes_to_generate:
            self.lipidClass.addItem(cls.__name__, cls)  

        if selection: # Selection exists if list item edited
            lipid = selection[0]  # Taken as list of [lipid,
            adduct = selection[1] # adduct]
            i = self.lipidClass.findText(lipid.lipid_class)
            # if i!= -1: it exists in the combo box.
            # if it exists, set combo box to that option.
            if i != -1: self.lipidClass.setCurrentIndex(i)
            i = self.lipidAdduct.findText(adduct)
            if i != -1: self.lipidAdduct.setCurrentIndex(i)

            if issubclass(self.lipidClass.currentData(), GL.Sphingolipid): j = 1
            else: j = 0 # Offset to let the base be first button 
            for i in range(lipid.No_Tails+j): # assign tails to buttons
                self.tailButton[i].tail = lipid.tails[i]
                self.tailButton[i].setText(lipid.tails[i].name)
            self.updateSpectra(lipid)

        self.lipidAdduct.currentTextChanged.connect(self.buildLipid)      

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
        self.tailButton = {}
        
        if issubclass(self.lipidClass.currentData(), GL.Sphingolipid):
            i = 1 # Create a button for the base first!
            self.tailButton[0] = QPushButton()
            self.tailButton[0].clicked.connect(self.specifyBase)
            self.hLayout3.addWidget(self.tailButton[0])
        else: i = 0 # Offset to let the base be first button 

        for button in range(self.lipidClass.currentData().No_Tails):
            self.tailButton[button+i] = QPushButton()
            self.tailButton[button+i].clicked.connect(self.specifyTail)
            self.hLayout3.addWidget(self.tailButton[button+i])

    def getTail(self, tail):
        # Signal captured from window
        self.tail = tail

    def specifyTail(self):
        # keep sender (button clicked)
        button = self.sender()
        tailWindow = TEW.TailWindow()
        tailWindow.output.connect(self.getTail)
        if tailWindow.exec() > 0:
        # tail stored in button for lipid
            button.tail = self.tail
        # button text updated to show tail
            button.setText(self.tail.name)
            self.buildLipid()

    def specifyBase(self):
        # keep sender (button clicked)
        button = self.sender()
        # sphingolipids can have several basetypes, use needed.
        types = self.lipidClass.currentData().base_types
        baseWindow = TEW.BaseWindow(types)
        baseWindow.output.connect(self.getTail)
        if baseWindow.exec() > 0:
        # tail stored in button for lipid
            button.tail = self.tail
        # button text updated to show tail
            button.setText(self.tail.name)
            self.buildLipid()

    def buildLipid(self):
        self.tails = []
        try: # for all the buttons on screen, grab the
             # tails that are assigned to them.
            for button in self.tailButton.values():
                self.tails.append(button.tail)
            lipidClass = self.lipidClass.currentData()
            self.lipid = lipidClass(*self.tails)
            self.updateSpectra(self.lipid)
        except: self.updateSpectra()

    def updateSpectra(self, lipid=None):
        try:
            adduct = self.lipidAdduct.currentData()
            lipid.resolve_spectra(adduct, lipid.adducts[adduct])
            formula = ''.join(''.join((key, str(val))) for (key, val) in lipid.formula.items())
            name = lipid.name+' '+adduct+' '+formula
            mz = GL.MA(lipid, adduct, 0).mass
            frags = lipid.spectra[adduct]
        except: name, mz, frags = '', 0, []
        self.spectra.setSpectra(name, mz, frags)
        self.table = Spectra.SpectraTableModel(frags)
        self.tableView.setModel(self.table)
        self.tableView.setItemDelegateForColumn(1, Spectra.SpinBoxDelegate(self.tableView))
        self.tableView.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.tableView.verticalHeader().hide()
        self.table.dataChanged.connect(self.updateLipid)

    def updateLipid(self, index: QModelIndex, *other):
        adduct = self.lipidAdduct.currentData()
        row = int(index.row())
        
        adducts = copy.deepcopy(self.lipid.adducts)
        data = self.table.tdata[row]
        adducts[adduct][data.fragmentType] = data.intensity
        self.lipid.adducts = adducts
        self.updateSpectra(self.lipid)

    def acceptLipid(self):
        self.done(1)
        self.close()
