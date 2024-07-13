import copy
from collections import Counter
import Wizard.Spectra as Spectra
import Lipids.GenerateLipids as GL
import Wizard.EditTail as ET
import Wizard.EditFragment as EF
from itertools import combinations_with_replacement as cwr, product

from PySide6.QtCore import Signal, QModelIndex, QAbstractTableModel, QMimeData
from PySide6.QtGui import Qt, QCursor, QDoubleValidator, QIntValidator, QAction
from PySide6.QtWidgets import  QDialog, QPushButton, QLineEdit, QComboBox, QTableView, QVBoxLayout, QHBoxLayout, QHeaderView, QMenu, QWidget, QSplitter, QApplication

class NewWindow(QDialog):
    '''
    Window showing lipid classes, permits editing
    '''

    def __init__(self, parent, classes):
        super().__init__()

        self.classList = classes 

        self.setWindowTitle('LSG3')
        self.setFixedSize(300, 510)
        self.vLayout = QVBoxLayout(self)
        self.hLayout = QHBoxLayout(self)
        self.hLayout2 = QHBoxLayout(self)

        self.tableView = QTableView()
        self.tableView.doubleClicked.connect(self.editClass)
        self.tableView.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.vLayout.addWidget(self.tableView)

        self.editLipid = QPushButton('Edit Class')
        self.editLipid.clicked.connect(self.editClass)
        self.hLayout2.addWidget(self.editLipid)
        self.vLayout.addLayout(self.hLayout2)

        self.buildList()

    def buildList(self):
        self.tableModel = ClassTableModel(self.classList)
        self.tableView.setModel(self.tableModel)

    def getClass(self, cls):
        if cls not in self.classList:
            self.classList.append(cls)
            self.buildList()

    def editClass(self):
        try:
            selection = self.tableView.selectedIndexes()[0].row()
            clsWindow = LipidWindow(self, [self.classList[selection]])
            clsWindow.output.connect(self.getClass)
            self.close()
            clsWindow.exec()
        except: pass # Invalid selection

    def removeSelectedClass(self):
        indexesToDelete = [row.row() for row in self.tableView.selectedIndexes()]
        indexesToDelete.reverse() # Reverse list of indexes so ones towards the end are deleted first
        for index in indexesToDelete:
            self.classList.pop(index)
        self.buildList()






class ClassTableModel(QAbstractTableModel):
    def __init__(self, data=[], parent=None):
        QAbstractTableModel.__init__(self, parent)
        self.tdata = data
        self.headers = ['Lipid Classes']

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
        cls = self.tdata[row]
        if role == Qt.DisplayRole:
            try: return str(cls.givenName)
            except: return str(cls.__name__)
        if role == Qt.ToolTipRole:
            try: toolTip = cls.tooltip
            except: 
                try: return str(cls.givenName)
                except: return str(cls.__name__)
            return toolTip

    def setData(self, index, value, role=Qt.EditRole):
        if role == Qt.EditRole:
            row = index.row()
            column = index.column()
            return True

        return QAbstractTableModel.setData(self, index, value, role)






class LipidWindow(QDialog):

    output = Signal(GL.Lipid)

    def __init__(self, parent, lipidClasses=[], classEdit=True, selection=None):
        super().__init__()
    
        self.setWindowFlag(Qt.WindowMinimizeButtonHint, True)
        self.setWindowFlag(Qt.WindowMaximizeButtonHint, True)

        self.parent = parent
        self.lipidClasses = lipidClasses
        self.classEdit = classEdit
        self.selection = selection
        self.defaultAdducts = [a for a in GL.adducts]

        self.setWindowTitle('LSG3')
        self.setMinimumWidth(600)
        self.setMinimumHeight(510)
        self.resize(600, 510)

        self.windowLayout = QVBoxLayout(self)
        self.windowLayout.setContentsMargins(0,0,0,0)

        self.splitter = QSplitter(Qt.Vertical)
        self.windowLayout.addWidget(self.splitter)
        self.splitter.setHandleWidth(3)
        self.splitter.setStyleSheet("""QSplitter::handle {
                                    background: qlineargradient(x1:0, y1:0, x2:1, y2:1,
                                        stop:0 #ccc, stop:1 #ccc);
                                    border: 1px solid #777;
                                    width: 15px;
                                    margin-top: 1px;
                                    margin-bottom: 1px;
                                    margin-left: 14px;
                                    margin-right: 14px;
                                    border-radius: 4px;
                                    }""")

        self.topWidget = QWidget(self.splitter)
        self.splitter.setCollapsible(0, False)
        self.vLayout = QVBoxLayout(self.topWidget)
        self.hLayout = QHBoxLayout(self.topWidget)
        self.vLayout.addLayout(self.hLayout)

        self.lipidBox = QComboBox()
        self.lipidBox.setToolTip('Lipid Class')
        self.lipidBox.currentIndexChanged.connect(self.updateSelectedClass)
        self.hLayout.addWidget(self.lipidBox)

        self.toolTipBox = QLineEdit()
        self.toolTipBox.setToolTip('Rename Lipid Class')
        self.toolTipBox.textChanged.connect(self.updateToolTip)
        self.hLayout.addWidget(self.toolTipBox)

        self.smilesBox = QLineEdit()
        self.smilesBox.setToolTip('Lipid SMILES')
        self.hLayout.addWidget(self.smilesBox)

        self.adductBox = QComboBox()
        self.adductBox.setToolTip('Current Adduct')
        self.adductBox.currentIndexChanged.connect(self.displayAdductDetails)
        self.adductBox.currentTextChanged.connect(self.buildLipid)
        self.hLayout2 = QHBoxLayout(self.topWidget)
        self.hLayout2.addWidget(self.adductBox)

        self.addAdduct = QPushButton('Add Adduct')
        self.addAdduct.clicked.connect(self.newAdduct)
        self.addAdduct.setAutoDefault(False)

        self.removeAdduct = QPushButton('Remove Adduct')
        self.removeAdduct.clicked.connect(self.delAdduct)
        self.removeAdduct.setAutoDefault(False)

        self.hLayout3 = QHBoxLayout(self.topWidget)
        self.hLayout3.addWidget(self.addAdduct)
        self.hLayout3.addWidget(self.removeAdduct)
        self.hLayout2.addLayout(self.hLayout3)
        self.vLayout.addLayout(self.hLayout2)

        self.adductName = QLineEdit()
        self.adductName.setToolTip('Adduct Name')
        self.adductName.textEdited.connect(self.updateAdductDetails)

        self.adductMass = QLineEdit() # Shouldn't be editable directly
        self.adductMass.setToolTip('Adduct Mass')

        self.adductCharge = QLineEdit()
        self.adductCharge.setToolTip('Adduct Charge')
        self.adductCharge.setValidator(QIntValidator())
        self.adductCharge.textEdited.connect(self.updateAdductDetails)

        self.adductFormula = QLineEdit()
        self.adductFormula.setToolTip('Adduct Formula\n'
                                      'eg: N:1,H:4')
        self.adductFormula.textEdited.connect(self.updateAdductDetails)

        self.hLayout4 = QHBoxLayout(self.topWidget)
        self.hLayout4.addWidget(self.adductName)
        self.hLayout4.addWidget(self.adductMass)
        self.hLayout4.addWidget(self.adductCharge)
        self.hLayout4.addWidget(self.adductFormula)
        self.vLayout.addLayout(self.hLayout4)

        self.spectra = Spectra.SpectraScatter(self)
        #self.spectra.setFixedSize(580, 200)
        self.spectra.setMinimumWidth(580)
        self.spectra.setMinimumHeight(200)
        #self.spectra.resize(580, 200)
        self.vLayout.addWidget(self.spectra)

        self.tailButton = {} # Buttons for tails stored in this dict
        self.hLayout5 = QHBoxLayout(self.topWidget)
        self.vLayout.addLayout(self.hLayout5)

        self.bottomWidget = QWidget(self.splitter)
        self.splitter.setCollapsible(1, False)
        self.splitter.setStretchFactor(1, 1)
        self.vLayout2 = QVBoxLayout(self.bottomWidget)

        self.tableView = QTableView()
        self.vLayout2.addWidget(self.tableView)
        self.tableView.setContextMenuPolicy(Qt.CustomContextMenu)
        self.tableView.customContextMenuRequested.connect(self.editAdducts)

        self.addFragmentButton = QPushButton('Add Fragment')
        self.addFragmentButton.clicked.connect(self.editAdducts)
        self.hLayout6 = QHBoxLayout(self.bottomWidget)
        self.hLayout6.addWidget(self.addFragmentButton)

        for cls in self.lipidClasses:
            try: self.lipidBox.addItem(cls.givenName, cls)
            except: self.lipidBox.addItem(cls.__name__, cls)
        if self.classEdit: # If the page is set up to modify the lipid class
            self.buildDefault()
        else: # If the page is set up to create / modify an individual lipid
            self.acceptButton = QPushButton('Accept')
            self.acceptButton.clicked.connect(self.acceptCreatedLipid)
            self.hLayout6.addWidget(self.acceptButton)

        self.vLayout2.addLayout(self.hLayout6)

    def updateSelectedClass(self):

        try: toolTip = self.lipidBox.currentData().tooltip
        except: toolTip = self.lipidBox.currentData().__name__
        self.toolTipBox.setText(toolTip)
        self.refreshAdductList()
        self.displayAdductDetails()
        self.createTailBoxes()

        if self.classEdit: # If the page is set up to modify the lipid class
            self.lipid = self.buildExampleLipid()
            self.updateSpectra(self.lipid)
        elif self.selection:
            self.lipid = self.selection[0]
            adduct = self.selection[1]
            i = self.adductBox.findText(adduct)
            if i != -1: self.adductBox.setCurrentIndex(i)
            for i in range(len(''.join(self.lipid.tailOrganisation))): # assign tails to buttons
                self.tailButton[i].tail = self.lipid.tails[i]
                self.tailButton[i].setText(self.lipid.tails[i].name)
            self.updateSpectra(self.lipid)

    def buildDefault(self):
        self.lipid = self.buildExampleLipid()
        self.updateSpectra(self.lipid)

    def acceptCreatedLipid(self):
        self.done(1)
        self.close()



    def createTailBoxes(self):
        for button in self.tailButton.keys():
            self.tailButton[button].setParent(None)
        self.tailButton = {}

        tailTypes = ''.join(self.lipidBox.currentData().tailOrganisation)
        i = 0
        for tailType in tailTypes:
            if tailType not in [',']:
                self.tailButton[i] = QPushButton()
                self.tailButton[i].tailType = tailType
                self.tailButton[i].setContextMenuPolicy(Qt.CustomContextMenu)
                self.tailButton[i].customContextMenuRequested.connect(self.tailButtonContextMenu)
                if tailType == 'B': 
                    self.tailButton[i].clicked.connect(self.specifyBase)
                else: 
                    self.tailButton[i].clicked.connect(self.specifyTail)
                self.hLayout5.addWidget(self.tailButton[i])
                i += 1

    def tailButtonContextMenu(self, pos):
        contextMenu = QMenu()
        defaultAction = QAction("Default", self)
        defaultAction.triggered.connect(self.buildDefault)
        contextMenu.addAction(defaultAction)
        contextMenu.exec_(self.sender().mapToGlobal(pos))


    def getTail(self, tail):
        # Signal captured from window
        self.tail = tail

    def specifyTail(self):
        # keep sender (button clicked)
        button = self.sender()
        tailWindow = ET.TailWindow(button.tailType)
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
        types = self.lipidBox.currentData().base_types
        baseWindow = ET.BaseWindow(types)
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
            lipidClass = self.lipidBox.currentData()
            self.lipid = lipidClass(*self.tails)
            try: self.lipid.name = self.lipid.specificname
            except: pass # For lipids with acyl in headgroup
            self.updateSpectra(self.lipid)
        except: self.updateSpectra()

    def updateLipid(self, index: QModelIndex, *other):

        adductType = self.adductBox.currentText() # Adduct
        lipidClass = self.lipidBox.currentData() # Lipid
        row = int(index.row()) # Location of edit in table
        adducts = copy.deepcopy(self.lipid.adducts) # Spectra
        data = self.table.tdata[row] # Edited data from table
        adducts[adductType][data.fragmentType] = data.intensity
        if self.classEdit:
            lipidClass.adducts = adducts # Overwrite Class
            self.buildLipid() # Rebuild Lipid
        else:
            self.lipid.adducts = adducts # Overwrite Lipid
            self.updateSpectra(self.lipid)

    def updateSpectra(self, lipid=None):

        try:
            adduct = self.adductBox.currentText()
            lipid.resolve_spectra(adduct, lipid.adducts[adduct])
            formula = ''.join(''.join((key, str(val))) for (key, val) in lipid.formula.items())
            name = lipid.name+' '+adduct+' '+formula
            mz = GL.MA(lipid, adduct, 0).mass
            frags = lipid.spectra[adduct]
            smiles = lipid.smiles
            self.smilesBox.setText(self.lipid.smiles)
        except: name, mz, frags, smiles = '', 0, [], ''
        self.spectra.setSpectra(name, mz, frags, smiles)
        for frag in frags:
            print(frag, frag.mass)
        self.table = Spectra.SpectraTableModel(frags)
        self.tableView.setModel(self.table)
        self.tableView.setItemDelegateForColumn(1, Spectra.SpinBoxDelegate(self.tableView))
        self.tableView.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.tableView.verticalHeader().hide()
        self.table.dataChanged.connect(self.updateLipid)




    def buildExampleLipid(self):
        '''
        Returns an example lipid for use.
        Lipid will be GPL (16:0_) 16:0_18:1\n
        or Sphingolipid with 18:0.
        '''

        constituents = self.lipidBox.currentData().tailOrganisation
        constituentList = [] # a combination of two tails, and another tail independent of the previous combination is needed.

        self.acyls =  [GL.sn(c=16, d=0, type='Acyl') ,
                       GL.sn(c=18, d=1, type='Acyl') ]
        self.ethers = [GL.sn(c=16, d=0, type='Ether'),
                       GL.sn(c=18, d=1, type='Ether') ]
        self.vinyls = [GL.sn(c=16, d=0, type='Vinyl'),
                       GL.sn(c=18, d=1, type='Vinyl') ]

        for x in constituents: # ie, for 'B', 'AA', 'A' in ['B', 'AA', 'A']
            group = dict(Counter(x)) # ie, {'B':1}, {'A':2}, {'A':1}
            for key in group:
                if key == 'B': # In the case of 'B', it indicates a base is needed. Get base list!
                    constituentList.append([GL.base(18, self.lipidBox.currentData().base_types[0])])
                elif key == 'A': # In the case of 'A', it indicates an acyl tail is needed. generate tail combination!
                    constituentList.append(cwr(self.acyls, r=group[key]))
                elif key == 'O': # In the case of 'O', it indicates an ether tail is needed. generate tail combination!
                    constituentList.append(cwr(self.ethers, r=group[key]))
                elif key == 'P': # In the case of 'P', it indicates a vinyl tail is needed. generate tail combination!
                    constituentList.append(cwr(self.vinyls, r=group[key]))
                else: pass

        _, comb, *_ = product(*constituentList)
        comb = self.flatten(comb)
        example =  self.lipidBox.currentData()(*comb)

        for i in range(len(''.join(constituents))): # assign tails to buttons
            self.tailButton[i].tail = example.tails[i]
            self.tailButton[i].setText(example.tails[i].name)

        example.resolve_spectra(self.adductBox.currentText(), self.adductBox.currentData())
        return example

    def flatten(self, data):
        if isinstance(data, tuple):
            for x in data: yield from self.flatten(x)
        else: yield data



    def editAdducts(self):
        menu = QMenu()
        menu.addAction('Add Predefined Fragment')
        menu.addSeparator()
        menu.addAction('Add Customised Fragment')
        menu.addSeparator()
        menu.addAction('Copy Value')
        menu.addSeparator()
        menu.addAction('Copy Entire Table')
        menu.addAction('... as .msp')
        selection = menu.exec(QCursor.pos())
        #try:

        if selection.text() == 'Add Predefined Fragment':
            fragment = EF.PredefinedFragment(self.lipid, self.adductBox.currentText())
            if fragment.exec() > 0:
                self.lipid.adducts[self.adductBox.currentText()][fragment.selection] = 0
                self.updateSpectra(self.lipid)

        elif selection.text() == 'Add Customised Fragment':
            fragment = EF.CustomisedFragment(self.lipid, self.adductBox.currentText())
            if fragment.exec() > 0:
                self.updateSpectra(self.lipid)

        elif selection.text() == 'Copy Value':
            selection = self.tableView.selectedIndexes()[0]
            if selection.isValid():
                cell_value = self.table.data(selection, Qt.DisplayRole)
                clipboard = QApplication.clipboard()
                clipboard.setText(str(cell_value))

        elif selection.text() == 'Copy Entire Table':
            rows = self.table.rowCount(self.tableView)
            columns = self.table.columnCount(self.tableView)
            data = self.lipid.name + '\t' + self.adductName.text() + '\n' + '\t'.join(self.table.headers) + '\n'

            for row in range(rows):
                intensity = self.table.data(self.table.index(row, 1))
                if intensity == '0':
                    continue
                for col in range(columns):
                    index = self.table.index(row, col)
                    data += str(self.table.data(index)) + '\t'
                data = data.strip() + '\n'

            data = data.strip()

            clipboard = QApplication.clipboard()
            mime_data = QMimeData()
            mime_data.setText(data)
            clipboard.setMimeData(mime_data)

        elif selection.text() == '... as .msp':

            adduct = self.adductBox.currentText()

            spectrum = [peak for peak in self.lipid.spectra[adduct] if peak.intensity !=0]

            if self.checklipidAmbiguity(spectrum):
                self.lipid.ambiguousName = f"{self.lipid.lipid_class} {self.addTailNames(self.lipid)}"
                self.lipid.ambiguoussmiles = ' '
            else: 
                self.lipid.ambiguousName = False
                self.lipid.ambiguoussmiles = False

            data = (f"NAME: {self.lipid.ambiguousName if self.lipid.ambiguousName else self.lipid.name} {adduct}\n"
                    f"IONMODE: {GL.adducts[adduct][1]}\n"
                    f"MW: {self.lipid.mass}\n"
                    f"PRECURSORMZ: {GL.MA(self.lipid, adduct, 0).mass}\n"
                    f"COMPOUNDCLASS: {self.lipid.lipid_class}\n"
                    f"FORMULA: {''.join(''.join((key, str(val))) for (key, val) in self.lipid.formula.items())}\n"
                    f"SMILES: {self.lipid.ambiguoussmiles if self.lipid.ambiguoussmiles else self.lipid.smiles}\n"
                    f"COMMENT: LSG in-silico\n" 
                    f"RETENTIONTIME: 0.00\n" # Pointless
                    f"PRECURSORTYPE: {adduct}\n"
                    f"Num Peaks: {len(spectrum)}\n")
            for peak in spectrum:
                data += f'{peak.mass} {peak.intensity} "{peak.Comment()}" \n'

            clipboard = QApplication.clipboard()
            mime_data = QMimeData()
            mime_data.setText(data)
            clipboard.setMimeData(mime_data)

        #except: pass # Invalid selection


    def addTailNames(self, lipid):
        c =  sum(snx.c  for snx in lipid.tails if snx.type != 'Headgroup')
        d =  sum(snx.d  for snx in lipid.tails if snx.type != 'Headgroup')
        me = sum(snx.me for snx in lipid.tails if snx.type != 'Headgroup')
        oh = sum(snx.oh for snx in lipid.tails if snx.type != 'Headgroup')
        dt = sum(snx.dt for snx in lipid.tails if snx.type != 'Headgroup')

        if 'Ether' in [snx.type for snx in lipid.tails if snx.type != 'Headgroup']:
            name = f"O-{c}:{d}"
        elif 'Vinyl' in [snx.type for snx in lipid.tails if snx.type != 'Headgroup']:
            name = f"P-{c}:{d}"
        else: name = f"{c}:{d}"
        
        if me > 0: # Methyl branching of fatty acid
            name += f";{me}-M" 
        if oh > 0: # Hydroxy functionalisation of fatty acid
            name += f";O{oh}"
        if dt > 0: # deuterium labelled fatty acids
            name += f"(D{dt})" # Deuterium doesn't update smiles currently.

        return name

    def checklipidAmbiguity(self, spectra):
        if not next((k.__class__.__name__ for k in spectra if 'FA' in k.__class__.__name__ or 'Cer' in k.__class__.__name__), False):
            return True
        else: return False

    def updateToolTip(self, text):
        try: self.lipidBox.currentData().tooltip = text
        except: pass

    def refreshAdductList(self):
        try: self.adductBox.clear()
        except: pass
        for adduct in self.lipidBox.currentData().adducts:
            self.adductBox.addItem(adduct, self.lipidBox.currentData().adducts[adduct])

    def newAdduct(self):
        menu = QMenu()
        adducts = GL.adducts
        for adduct in adducts:
            if adduct not in self.lipidBox.currentData().adducts:
                menu.addAction(adduct)
        menu.addSeparator()
        menu.addAction("Custom ...")
        selection = menu.exec(QCursor.pos())
        if selection and selection.text() != "Custom ...":
            self.newSpectra(selection.text())
        if selection and selection.text() == "Custom ...":
            self.customSpectra()

    def newSpectra(self, adduct):
        self.lipidBox.currentData().adducts[adduct] = {}
        self.refreshAdductList()
        self.adductBox.setCurrentText(adduct)

    def customSpectra(self):
        self.lipidBox.currentData().adducts['[M]'] = {}
        GL.adducts['[M]'] = [0,'Positive',1,{},'']

        self.refreshAdductList()
        self.adductBox.setCurrentText('[M]')

    def delAdduct(self):
        menu = QMenu()
        for adduct in self.lipidBox.currentData().adducts:
            menu.addAction(adduct)
        selection = menu.exec(QCursor.pos())
        try: 
            self.lipidBox.currentData().adducts.pop(selection.text())
            self.refreshAdductList()
        except: pass

    def displayAdductDetails(self):
        selectedAdduct = self.adductBox.currentText()

        try:
            self.adductName.setText(selectedAdduct)
            self.adductMass.setText(str(GL.adducts[selectedAdduct][0]))
            self.adductCharge.setText(str(GL.adducts[selectedAdduct][2]))
            formula = ','.join(':'.join((key, str(val))) for (key, val) in GL.adducts[selectedAdduct][3].items())
            self.adductFormula.setText(formula)
        except: pass

        if selectedAdduct in self.defaultAdducts:
            self.adductName.setDisabled(True)
            self.adductMass.setDisabled(True)
            self.adductCharge.setDisabled(True)
            self.adductFormula.setDisabled(True)
        else:
            self.adductName.setDisabled(False)
            self.adductMass.setDisabled(True)
            self.adductCharge.setDisabled(False)
            self.adductFormula.setDisabled(False)
    
    def updateAdductDetails(self):
        selectedAdduct = self.adductBox.currentText()
        try:
            
            if self.adductName.text() not in GL.adducts.keys():
                newAdduct = self.adductName.text()
                GL.adducts[newAdduct] = GL.adducts.pop(selectedAdduct)
                self.lipidBox.currentData().adducts[newAdduct] = self.lipidBox.currentData().adducts.pop(selectedAdduct)
                self.refreshAdductList()
                self.adductBox.setCurrentText(newAdduct)
                self.adductName.setFocus()

            temp = {}
            formula = self.adductFormula.text()
            formula = formula.split(',')
            formula = [element.strip() for element in formula]
            mass = 0
            for f in formula:
                key, val = f.split(':')
                mass += GL.elements[key][0][0]*int(val)
                temp[str(key)] = int(val)
            GL.adducts[selectedAdduct][3] = temp
            mass -= int(self.adductCharge.text())*GL.masses['e']

            self.adductMass.setText(str(round(mass, 6)))
            GL.adducts[selectedAdduct][0] = mass
            GL.adducts[selectedAdduct][2] = int(self.adductCharge.text())
            if GL.adducts[selectedAdduct][2] > 0:
                GL.adducts[selectedAdduct][1] = 'Positive'
            elif GL.adducts[selectedAdduct][2] < 0:
                GL.adducts[selectedAdduct][1] = 'Negative'
            elif GL.adducts[selectedAdduct][2] == 0:
                GL.adducts[selectedAdduct][1] = 'Positive'
                GL.adducts[selectedAdduct][2] = 1

        except: pass

