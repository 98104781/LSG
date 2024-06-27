
from scipy.stats import multinomial

from collections import Counter
import Wizard.ResourcePath as RP
import Lipids.GenerateLipids as GL
from Lipids.GenerateLipids import elements
import Wizard.EditLipidAdduct as LAEW
import Wizard.Spectra as Spectra

from PySide6.QtCore import QObject, Signal, QThread, Qt, Signal, Property, QMimeData
from PySide6.QtGui import QDoubleValidator, QCursor, QColor, QStandardItem, QStandardItemModel
from PySide6.QtWidgets import QProgressBar,  QPlainTextEdit, QPushButton, QVBoxLayout, QHBoxLayout, QLineEdit, QComboBox, QWizardPage, QTreeWidget, QTreeWidgetItem, QMenu, QSplitter, QWidget, QApplication, QDialog, QTableView, QHeaderView, QSlider, QLabel

class Page(QWizardPage):
    '''
    Final page of GUI. Allows for generation.
    '''
    treeDataChanged = Signal()

    def __init__(self, parent=None):
        super(Page, self).__init__(parent)

        self.wiz = parent
        self.classQbox = {} # GPL QCheckBoxes        
        self.classAdductQbox = {} # Adduct QCheckBoxes
        self.candidatesString = []
        self.candidatesData = []

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
        self.hLayout1 = QHBoxLayout(self.topWidget)
        self.vLayout.addLayout(self.hLayout1)

        self.ionType = QComboBox()
        self.ionType.addItems(['+/-', '+', '-'])
        self.hLayout1.addWidget(self.ionType)

        self.massInput = QLineEdit()
        self.massInput.setPlaceholderText('Target Mass')
        self.massInput.setValidator(QDoubleValidator(44, 5000, 6))
        self.hLayout1.addWidget(self.massInput)

        self.hLayout1.addStretch()

        self.errorInput = QLineEdit()
        self.errorInput.setPlaceholderText('Permitted Error')
        self.errorInput.setValidator(QDoubleValidator(0, 100, 6))
        self.errorInput.setText('15')
        self.hLayout1.addWidget(self.errorInput)

        self.errorType = QComboBox()
        self.errorType.addItems(['ppm', 'Da'])
        self.hLayout1.addWidget(self.errorType)

        self.guesstimateButton = QPushButton("Guesstimate")
        self.guesstimateButton.setStyleSheet("""QPushButton {
                                                    background-color: #c7ecee;
                                                    border: 1px solid #000000;
                                                    border-radius: 2px;
                                                    padding: 2px 4px;
                                                }
                                                
                                                QPushButton:hover {
                                                    background-color: #dcdcdc;
                                                }
                                                
                                                QPushButton:pressed {
                                                    background-color: #b6b6b6;
                                                } """)
        self.guesstimateButton.clicked.connect(self.guesstimate)
        self.vLayout.addWidget(self.guesstimateButton)

        self.output_console = QPlainTextEdit()
        self.output_console.setStyleSheet("""QPlainTextEdit {
                                                background-color: #000000;
                                                color: #f0f0f0;
                                                font-family: Courier, monospace;
                                                font-size: 10px;
                                                font-weight: bold;
                                                border: 1px solid #000000;}
                                                                    
                                            QPlainTextEdit:focus {
                                                border: 1px solid #00afff;} """)
        self.output_console.setReadOnly(True)
        self.vLayout.addWidget(self.output_console)

        self.bottomWidget = QWidget(self.splitter)
        self.splitter.setCollapsible(1, False)
        self.vLayout2 = QVBoxLayout(self.bottomWidget)

        self.treeView = QTreeWidget()
        self.treeView.setContextMenuPolicy(Qt.CustomContextMenu)
        self.treeView.customContextMenuRequested.connect(self.editLipidContextMenu)
        self.treeView.doubleClicked.connect(self.editLipid)
        self.treeView.setHeaderHidden(True)
        self.registerField("tree2", self, "tree_property")
        self.treeView.itemChanged.connect(self.itemChanged)
        self.treeView.itemSelectionChanged.connect(self.updateCheckedState)
        self.vLayout2.addWidget(self.treeView)

        self.hLayout2 = QHBoxLayout()
        self.vLayout2.addLayout(self.hLayout2)

        self.reviewButton = QPushButton("Review Candidates")
        self.reviewButton.clicked.connect(self.reviewCandidates)
        self.reviewButton.setEnabled(False)
        self.hLayout2.addWidget(self.reviewButton)

        self.modifybutton = QPushButton("Modify Lipid Templates")
        self.modifybutton.clicked.connect(self.editLipid)
        self.hLayout2.addWidget(self.modifybutton)



        self.progress_bar = QProgressBar()
        self.progress_bar.setStyleSheet(""" QProgressBar {
                                                border: 1px solid grey;
                                                border-radius: 1px;
                                                background-color: #ffffff;
                                                text-align: center;}
                                                                
                                            QProgressBar::chunk {
                                                background-color: #c7ecee;}""")      
        self.vLayout2.addWidget(self.progress_bar)
        self.progress_bar.reset()

        self.setFinalPage(True)

        self.white =    "<span style=' font-size:8pt; font-weight:600; color:white;' >"
        self.red =      "<span style=' font-size:8pt; font-weight:600; color:red;' >"
        self.blue =     "<span style=' font-size:8pt; font-weight:600; color:blue;' >"
        self.yellow =   "<span style=' font-size:8pt; font-weight:800; color:yellow;' >"
        self.green =    "<span style=' font-size:8pt; font-weight:600; color:limegreen;' >"
        self.orange =   "<span style=' font-size:8pt; font-weight:600; color:orange;' >"

        self.generatorThread = QThread() # Generator on second thread so GUI doesn't lag!

        return super().initializePage()

    def itemChanged(self, item):
        if item.checkState(0) == Qt.Unchecked:
            item.setBackground(0, Qt.white)
        else: 
            item.setBackground(0, QColor('#c7ecee'))
        if self.pageInitialised:
            self.updateConsole()

    def updateCheckedState(self):
        if QApplication.mouseButtons() == Qt.LeftButton:
            items = self.treeView.selectedItems()
            for item in items:
                if item.checkState(0) == Qt.Unchecked:
                    statusUpdate = Qt.Checked
                else: statusUpdate = Qt.Unchecked
                item.setCheckState(0, statusUpdate)
                self.itemChanged(item)
                #if statusUpdate == Qt.Checked:
                #    item.setBackground(0, QColor(180, 230, 180))
                #else: item.setBackground(0, Qt.white)

    def generatorComplete(self):

        self.guesstimateButton.setEnabled(True)
        self.treeView.setEnabled(True)
        self.modifybutton.setEnabled(True)
        self.errorType.setEnabled(True)
        self.errorInput.setEnabled(True)
        self.massInput.setEnabled(True)
        self.ionType.setEnabled(True)

        self.progress_bar.setMaximum(1)
        self.progress_bar.setValue(1)
        if self.output_console.toPlainText() == '':
            self.output_console.appendHtml(f"{self.yellow}No candidates found.")
        else: self.reviewButton.setEnabled(True)

        self.hasGenerated = True
        self.generatorThread.quit()

    def appendToConsole(self, string):
        self.candidatesString.append(string)
        self.output_console.appendHtml(string)

    def storeCandidate(self, data):
        self.candidatesData.append(data)

    def initializePage(self) -> None:

        self.hasGenerated = False
        self.pageInitialised = False

        self.guesstimateButton.setEnabled(True)
        try:
            if self.generatorThread.isRunning():
                self.progress_bar.setValue(0)
                self.progress_bar.setMaximum(0)
                self.guesstimateButton.setEnabled(False)
        except: pass

        self.treeView.clear()

        self.classList = self.wiz.classes_to_generate.copy()

        for cls in self.classList: #  Make boxes for Treeview
            self.classQbox[cls]  =  QTreeWidgetItem(self.treeView)
            root = self.classQbox[cls] # Creates tickbox for class
            root.lipidClass = cls # Custom variable to store class
            try: root.setText(0, cls.givenName) # Checks for names 
            except: root.setText(0, cls.__name__) # Else make name 
            root.setCheckState(0, Qt.Unchecked) #   Untick tickbox
            root.setFlags(root.flags() | Qt.ItemIsAutoTristate | Qt.ItemIsUserCheckable)
            try: root.setToolTip(0, cls.tooltip) #  Set tooltip to 
            except: root.setToolTip(0, cls.__name__) #  lipid name
            self.classAdductQbox[cls] = {} # Open dict for adducts

            for adduct in cls.adducts: # Make sub-box for treeview
                self.classAdductQbox[cls][adduct] = QTreeWidgetItem(self.classQbox[cls])
                child = self.classAdductQbox[cls][adduct] # Assign
                child.fragmentList = cls.adducts[adduct] # Adducts
                child.setText(0, adduct) #  Gives name for tickbox
                try: child.setCheckState(0, self.aBoxStatus[cls][adduct])
                except: child.setCheckState(0, Qt.Unchecked) # Tickstate
                child.setFlags(child.flags() | Qt.ItemIsUserCheckable) 

        self.classes_to_generate = [] # Empty lists initialised here
        self.tails_to_generate = [] # otherwise error thrown when
        self.bases_to_generate = [] # coming from page 2B

        # Prepare information to generate lipids:
        # List with limits for generated tails
        self.tails_to_generate = [int(self.field('cmin') or 2), int(self.field('cmax') or 30),
                                  int(self.field('dmin') or 0), int(self.field('dmax') or 12),
                                  int(self.field('omax') or 0), int(self.field('Umax') or 0)]

        base_types = [] # List with types of sphingoid bases to generate
        for cls in self.classList: # Sphingolipids can have one of many base type, which the lipid
            if issubclass(cls, GL.Sphingolipid): # is centered around. The possible types are defined in the
                base_types.extend(cls.base_types) # Sphingolipid.base_types list. Collect all unique base types 
            base_types = list(set(base_types)) # used so that they can be generated with the lipids.

        if self.field('ceramideVariability') is False:
            self.bases_to_generate = [18, 18, base_types]
        else: self.bases_to_generate = [max(int(self.field('cmin') or 0), 7), max(int(self.field('cmax') or 0), 7), base_types]

        self.acyls = []
        self.ethers = []
        self.vinyls = []
        self.bases = []

        self.output_console.clear() # Clear console, update with tails and lipids chosen
        self.progress_bar.reset()

        self.updateConsole()

        if (self.candidatesString != []):
            self.output_console.setPlainText('')
            for string in self.candidatesString:
                self.output_console.appendHtml(string)

        self.pageInitialised = True

        return super().initializePage()

    def editLipidContextMenu(self):
        menu = QMenu()
        menu.addAction("Edit Lipid Template")
        menu.addSeparator()
        menu.addAction("Clear Selection")
        selection = menu.exec(QCursor.pos())
        if selection != None:
            if selection.text() == "Edit Lipid Template":
                self.editLipid()
            elif selection.text() == "Clear Selection":
                self.aBoxStatus={}
                self.initializePage()

    def editLipid(self):

        self.aBoxStatus={} # If edit lipid is called, recall all tickstates for adduct boxes
        for cls in self.classList:
            self.aBoxStatus[cls] = {}
            for adduct in cls.adducts:
                self.aBoxStatus[cls][adduct] = self.classAdductQbox[cls][adduct].checkState(0)

        try:
            cls = self.treeView.selectedItems()[0]
            if cls.parent():
                cls = cls.parent().lipidClass
            else: cls = cls.lipidClass
            clsWindow = LAEW.LipidWindow(self, [cls])
            clsWindow.exec()
            self.initializePage()
        except: 
            try:
                classes_to_generate = self.wiz.classes_to_generate
                editLipidWindow = LAEW.NewWindow(self, classes_to_generate)
                editLipidWindow.exec()
                self.initializePage()
            except: pass

        pass

    def reviewCandidates(self):
        
        if self.candidatesData != []:
            #try:
            newWindow = NewWindow(self, self.candidatesData)
            newWindow.exec()
            #except: pass


    def updateConsole(self):
        
        self.classes_to_generate = {}
        self.selected_class_adducts = self.field('tree2')

        if (self.selected_class_adducts != {}):

            self.output_console.setPlainText('') # Clear console, begin listing classes / adducts selected.
            self.output_console.appendHtml(f"<ul>{self.yellow}> {self.white}The following will be assessed:")

            caString = "" # Generates string, for class and adduct list
            for item, item2 in self.selected_class_adducts.items():
                caString += f"<li>{self.yellow}> {self.white}{item.text(0)}{self.orange}"
                self.classes_to_generate.update({item.lipidClass:[]})
                adducts_to_generate = []

                for adduct in item2: # Update the selected adducts
                    adducts_to_generate.append(adduct.text(0))
                    caString += f" {adduct.text(0)}"
                    item.lipidClass.adducts[adduct.text(0)] = adduct.fragmentList
                caString+="</li>"

                self.classes_to_generate[item.lipidClass] = adducts_to_generate
            self.output_console.appendHtml(f"{caString}</ul>")
            self.reviewButton.setEnabled(False)

        else: self.output_console.setPlainText('')

    def guesstimate(self):

        massValid, mass, _ = self.massInput.validator().validate(self.massInput.text(), 0)
        errorValid, error, _ = self.errorInput.validator().validate(self.errorInput.text(), 0)
        self.output_console.setPlainText('')
        self.candidatesString = []
        self.candidatesData = []

        if massValid != QDoubleValidator.Acceptable:
            self.output_console.appendHtml(f"{self.yellow}Please enter valid mass ( 44 - 5000 ).\n")
        elif errorValid != QDoubleValidator.Acceptable or (float(error) > 2 and self.errorType.currentText() == "Da") or (float(error) > 100 and self.errorType.currentText() == "ppm"):
            self.output_console.appendHtml(f"{self.yellow}Please specify valid error ( 0 - 100 ppm, 0 - 2 Da).\n")
        elif len(self.classes_to_generate) < 1:
            self.output_console.appendHtml(f"{self.yellow}Please select at least one lipid or adduct.\n")
        else:

            self.guesstimateButton.setEnabled(False)
            self.treeView.setEnabled(False)
            self.modifybutton.setEnabled(False)
            self.errorType.setEnabled(False)
            self.errorInput.setEnabled(False)
            self.massInput.setEnabled(False)
            self.ionType.setEnabled(False)
            self.reviewButton.setEnabled(False)
            self.hasGenerated = False

            self.output_console.setPlainText('')
            self.targetMass = float(mass)
            self.candidateMinMass, self.candidateMaxMass = self.determineMassBounds(self.targetMass, float(error), self.errorType.currentText())

            self.progress_bar.setValue(0) # Sets loading animation
            self.progress_bar.setMaximum(0) # for progress bar

            self.generatorObject = Generator(self.candidateMinMass, self.candidateMaxMass, self.targetMass, self.classes_to_generate,
                                              self.ionType.currentText(), self.tails_to_generate, self.bases_to_generate)
            self.generatorObject.moveToThread(self.generatorThread)
            self.generatorThread.started.connect(self.generatorObject.run)
            self.generatorObject.consoleOutput.connect(self.appendToConsole)
            self.generatorObject.candidateData.connect(self.storeCandidate)
            self.generatorObject.finished.connect(self.generatorComplete)
            self.generatorObject.finished.connect(self.generatorObject.deleteLater)
            #self.generatorThread.finished.connect(self.generatorThread.deleteLater)
            self.generatorThread.start()
            self.completeChanged.emit()

    def determineMassBounds(self, mass, error, errorType):
        if errorType == 'ppm':
            massShift = mass*(error / 1000000)
            return  mass-massShift, mass+massShift
        elif errorType == 'Da':
            return mass-error, mass+error

    def nextId(self):
        return -1

    def treeData(self):
        '''
        Custom storage for tree data.
        Allows later determination of selected class/adducts.
        '''
        checkedBoxes = {}
        root = self.treeView.invisibleRootItem()

        childs_1 = root.childCount()
        for i in range(childs_1):
            item = root.child(i)
            if bool(item.checkState(0)):
                checkedBoxes[item] = []
                #item.lipidClass
                
                adducts = []
                childs_2 = item.childCount()
                for j in range(childs_2):
                    item2 = item.child(j)
                    if bool(item2.checkState(0)):
                        adducts.append(item2)
                        #item2.fragmentList
                checkedBoxes[item] = adducts

        return checkedBoxes

    def setTreeData(self, data):
        pass

    tree_property = Property("QVariant", treeData, setTreeData, treeDataChanged)

    def isComplete(self):
        if self.hasGenerated:
            try: # Disable finish button if running.
                if self.generatorThread.isRunning():
                    return False
            except: pass # Thread not yet created.
        else: return False
        return super().isComplete()



class Generator(QObject):
    
    consoleOutput = Signal(str)
    candidateData = Signal(list)
    finished = Signal()

    def __init__(self, candidateMinMass, candidateMaxMass, targetMass, classes_to_generate, ionMode, tails_to_generate, bases_to_generate):
        super().__init__()

        # self.tails_to_generate = [int(self.field('cmin') or 2), int(self.field('cmax') or 30),
        #                           int(self.field('dmin') or 0), int(self.field('dmax') or 12),
        #                           int(self.field('omax') or 0), int(self.field('Umax') or 0)]

        self.classes_to_generate = classes_to_generate
        self.cmin, self.cmax, self.dmin, self.dmax, self.omax, self.Umax = tails_to_generate
        self.b2G = bases_to_generate
        self.ionMode = ionMode

        self.candidateMinMass = candidateMinMass
        self.candidateMaxMass = candidateMaxMass
        self.targetMass = targetMass

        self.white =    "<span style=' font-size:8pt; color:white;' >"
        self.red =      "<span style=' font-size:8pt; color:red;' >"
        self.blue =     "<span style=' font-size:8pt; color:blue;' >"
        self.yellow =   "<span style=' font-size:8pt; color:yellow;' >"
        self.green =    "<span style=' font-size:8pt; color:limegreen;' >"
        self.orange =   "<span style=' font-size:8pt; color:orange;' >"

    def run(self):
        
        #import pydevd;pydevd.settrace(suspend=False)

        prevNumKeys = 0

        for lipidClass in self.classes_to_generate:

            constituents = lipidClass.tailOrganisation # List of tails and organisation: ['B', 'AA', 'A'], indicates a sphingoid Base, -

            numKeys = 0
            for x in constituents: # ie, for 'B', 'AA', 'A' in ['B', 'AA', 'A']
                group = dict(Counter(x)) # {'B':1}, {'A':2}, {'A':1}

                for key in group:
                    numKeys += group[key]

            constituentList = [] # a combination of two tails, and another tail independent of the previous combination is needed.
            minimumUnits = []

            for x in constituents:
                group = dict(Counter(x))

                for key in group:
                    if key == 'B':
                        self.bases = GL.generate_base_tails(self.b2G)
                        basesToUse = list(self.flatten([self.bases[basetype] if basetype in self.bases.keys() else self.bases[basetype] for basetype in lipidClass.base_types]))
                        basesToUse.sort()
                        [constituentList.append(basesToUse) for n in range(group[key])]
                        [minimumUnits.append(basesToUse) for n in range(group[key])]
                    elif key == 'A':
                        if numKeys != prevNumKeys:
                            self.acyls = GL.generate_tails([self.cmin, (self.cmax*numKeys), self.dmin, self.dmax*numKeys, self.omax*numKeys, self.Umax*numKeys], 'Acyl')
                            self.acyls.sort()
                        [constituentList.append(self.acyls) for n in range(group[key])]
                        [minimumUnits.append([GL.sn(1, type='Acyl')]) for n in range(group[key])]
                    elif key == 'O':                       
                        self.ethers = GL.generate_tails([self.cmin, (self.cmax*numKeys), self.dmin, self.dmax*numKeys, self.omax*numKeys, self.Umax*numKeys], 'Ether')
                        self.ethers.sort()
                        [constituentList.append(self.ethers) for n in range(group[key])]
                        [minimumUnits.append([GL.sn(1, type='Ether')]) for n in range(group[key])]
                    elif key == 'P':                      
                        self.vinyls = GL.generate_tails([self.cmin, (self.cmax*numKeys), self.dmin, self.dmax*numKeys, self.omax*numKeys, self.Umax*numKeys], 'Vinyl')
                        self.vinyls.sort()
                        [constituentList.append(self.vinyls) for n in range(group[key])]
                        [minimumUnits.append([GL.sn(1, type='Vinyl')]) for n in range(group[key])]
                    else: pass

            tailComb = [tailType[0] for tailType in constituentList]
            minLipid = lipidClass(*tailComb)
            tailComb[-1] = constituentList[-1][-1]
            maxLipid = lipidClass(*tailComb)

            for adduct in self.classes_to_generate[lipidClass]:

                if self.ionMode == '-' and GL.adducts[adduct][2] > 0:
                    continue
                if self.ionMode == '+' and GL.adducts[adduct][2] < 0:
                    continue

                minLipidMass = GL.MA(minLipid, adduct, 0).mass
                maxLipidMass = GL.MA(maxLipid, adduct, 0).mass
                lD = self.candidateMaxMass - minLipidMass # distance to lower threshold from candidate
                uD = maxLipidMass - self.candidateMinMass # distance to upper threshold from candidate

                if lD < 0 or uD < 0: # negative distance means candidate is outside of thresholds.
                    continue
                else:
                    tailComb = [tailType[0] for tailType in minimumUnits] # Gotta 
                    candidates = self.binSearchBetweenVals(self.candidateMinMass, self.candidateMaxMass, lipidClass, constituentList, minimumUnits, adduct)
                    for candidate in candidates:
                        
                        tailComb[-1] = constituentList[-1][candidate]
                        candidate = lipidClass(*tailComb)
                        candidateMass = GL.MA(candidate, adduct, 0).mass
                        candidateFormula = candidate.formula.copy()
                        candidateFormula += GL.adducts[adduct][3]
                        cleanedFormula = ''.join(''.join((key, str(val if val > 1 else ''))) for (key, val) in candidateFormula.items()) 
                        if GL.adducts[adduct][2] > 0: candidateCharge = f"{abs(GL.adducts[adduct][2])}+" if abs(GL.adducts[adduct][2]) > 1 else "+"
                        else: candidateCharge = f"{abs(GL.adducts[adduct][2])}-" if abs(GL.adducts[adduct][2]) > 1 else "-"

                        className = candidate.lipid_class
                        tailName = self.addTailNames(candidate)
                        massDiff = self.diffPPM(candidateMass, self.targetMass)

                        self.consoleOutput.emit(f"{self.yellow}> {self.white}{className} {tailName}, "
                                                        f"{self.orange}{adduct}, {self.yellow}[{cleanedFormula}]{candidateCharge}, " 
                                                        f"{self.green}{candidateMass}, {self.white}({massDiff} ppm)")
                        
                        self.candidateData.emit([f"{className} {tailName}", f"{adduct}", f"[{cleanedFormula}]", candidateFormula, candidateCharge, abs(GL.adducts[adduct][2]), candidateMass, massDiff])

        self.finished.emit()
        return

    def flatten(self, data):
        if isinstance(data, tuple) or isinstance(data, list):
            for x in data: yield from self.flatten(x)
        else: yield data

    def binSearchBetweenVals(self, leftLimit, rightLimit, lipidClass, constituents, minimum, adduct):

        def binSearch(arr, target, low):
            
            high = len(arr)
            while low < high:
                mid = (low + high) // 2
                tailComb[-1] = constituents[-1][mid]
                candidate = lipidClass(*tailComb)
                candidateMass = GL.MA(candidate, adduct, 0).mass
                if candidateMass < target:
                    low = mid + 1
                else:
                    high = mid
            return low

        tailComb = [tailType[0] for tailType in minimum]
        maxKey = len(constituents[-1])
        arr = list(range(maxKey))

        # Find first element >= lhs
        leftIndex = binSearch(arr, leftLimit, 0)

        # Find first element > rhs
        rightIndex = binSearch(arr, rightLimit, leftIndex)

        # Return elements between lhs and rhs
        return arr[leftIndex:rightIndex]

    def addTailNames(self, lipid):
        c =  sum(snx.c  for snx in lipid.tails if snx.type != 'Headgroup')
        d =  sum(snx.d  for snx in lipid.tails if snx.type != 'Headgroup')
        me = sum(snx.me for snx in lipid.tails if snx.type != 'Headgroup')
        oh = sum(snx.oh for snx in lipid.tails if snx.type != 'Headgroup')
        dt = sum(snx.dt for snx in lipid.tails if snx.type != 'Headgroup')

        name = f"{c}:{d}"
        
        if me > 0: # Methyl branching of fatty acid
            name += f";{me}-M" 
        if oh > 0: # Hydroxy functionalisation of fatty acid
            name += f";O{oh}"
        if dt > 0: # deuterium labelled fatty acids
            name += f"(D{dt})" # Deuterium doesn't update smiles currently.
        return name
    
    def diffPPM(self, candidate, target):
        return round(((candidate - target)/(target)) * 1000000, 2)

class NewWindow(QDialog):
    '''
    Window showing lipid classes, permits editing
    '''

    def __init__(self, parent, candidates):
        super().__init__()

        self.setWindowTitle('LSG3')
        self.setMinimumSize(600, 510)

        self.windowLayout = QVBoxLayout(self)
        self.windowLayout.setContentsMargins(0,0,0,0)

        self.splitter = QSplitter(Qt.Horizontal)
        self.windowLayout.addWidget(self.splitter)
        self.splitter.setHandleWidth(3)
        self.splitter.setStyleSheet("""QSplitter::handle {
                                    background: qlineargradient(x1:0, y1:0, x2:1, y2:1,
                                        stop:0 #ccc, stop:1 #ccc);
                                    border: 1px solid #777;
                                    width: 15px;
                                    margin-top: 14px;
                                    margin-bottom: 14px;
                                    margin-left: 1px;
                                    margin-right: 1px;
                                    border-radius: 4px;
                                    }""")

        self.leftWidget = QWidget(self.splitter)
        self.splitter.setCollapsible(0, False)
        self.vLayout = QVBoxLayout(self.leftWidget)

        self.spectra = Spectra.SpectraScatter(self, False)
        self.spectra.xaxis.setGridLineVisible(True)
        self.spectra.yaxis.setGridLineVisible(True)
        self.vLayout.addWidget(self.spectra)

        self.hLayout = QHBoxLayout(self)
        self.vLayout.addLayout(self.hLayout)

        self.valueLabel = QLabel()
        self.valueLabel.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.valueLabel.setFixedWidth(40)
        self.hLayout.addWidget(self.valueLabel)

        self.resSlider = QSlider(Qt.Horizontal)
        self.resSlider.setMinimum(5000)
        self.resSlider.setMaximum(300000)
        self.resSlider.setTickPosition(QSlider.TicksBelow)
        self.resSlider.setTickInterval(14750)     
        self.resSlider.setValue(5000)
        self.resSlider.sliderReleased.connect(self.changeResolution)
        self.hLayout.addWidget(self.resSlider)

        self.valueLabel.setText(str(self.resSlider.sliderPosition()))

        self.candidateTable = QTableView()
        tableData = QStandardItemModel(len(candidates), 2, parent=None)
        headers = ['Candidate', 'Diff (ppm)']
        tableData.setHorizontalHeaderLabels(headers)
        for row in range(len(candidates)):
            for col, value in enumerate([0, -1]):
                item = QStandardItem(str(candidates[row][value]))
                item.setFlags(Qt.ItemIsEnabled | Qt.ItemIsSelectable)
                item.setData(candidates[row])
                tableData.setItem(row, col, item)
        self.candidateTable.setModel(tableData)
        self.candidateTable.setSortingEnabled(True)
        self.candidateTable.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.candidateTable.setSelectionBehavior(QTableView.SelectRows)
        self.candidateTable.setSelectionMode(QTableView.SingleSelection)
        self.candidateTable.selectionModel().selectionChanged.connect(self.updateSpectra)

        self.vLayout.addWidget(self.candidateTable)

        self.rightWidget = QWidget(self.splitter)
        self.splitter.setCollapsible(1, False)
        self.vLayout2 = QVBoxLayout(self.rightWidget)

        self.spectraTable = QTableView()
        self.spectraTable.setContextMenuPolicy(Qt.CustomContextMenu)
        self.spectraTable.customContextMenuRequested.connect(self.copyData)
        self.vLayout2.addWidget(self.spectraTable) 


    def updateSpectra(self, item):
        index = item.indexes()[0]
        selection = self.candidateTable.model().itemFromIndex(index)
        rowData = selection.data()

        self.isotopeSpectra = self.determineIsotopeSpectra(rowData[-5], rowData[-3], .0001)

        self.spectra.setSpectra(f"{rowData[0]} {rowData[1]} {rowData[2]}{rowData[4]}", rowData[-2], self.isotopeSpectra, isotopeSpectra=True)
        self.spectra.drawBellCurves(self.isotopeSpectra, self.resSlider.sliderPosition())

        tdata = sorted(self.isotopeSpectra, key=lambda x: x[0])

        self.tableData = QStandardItemModel(len(tdata), 2, parent=None)
        self.spectraTable.setModel(self.tableData)
        self.spectraTable.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.spectraTable.verticalHeader().hide()
        headers = ['m/z (Da)', 'Abn. (%)']
        self.tableData.setHorizontalHeaderLabels(headers)
        for row in range(len(tdata)):
            for col in [0, 1]:
                if col == 1: item = QStandardItem(str(round(tdata[row][col], 2)))
                else: item = QStandardItem(str(tdata[row][col]))
                item.setFlags(Qt.ItemIsEnabled | Qt.ItemIsSelectable)
                item.setData(tdata[row])
                self.tableData.setItem(row, col, item)
                  
    def changeResolution(self):
            self.valueLabel.setText(str(self.resSlider.sliderPosition()))
            if self.spectra.curveDisplayed:
                self.spectra.drawBellCurves(self.isotopeSpectra, self.resSlider.sliderPosition())

    def sumContributions(self, elementContributions):
        
        def recursionLoop(depth, sum, intensity):
            if depth == len(elementContributions):
                contributions.append([sum, intensity])
                return
            for isotope in elementContributions[depth]:
                recursionLoop(depth + 1,
                                    sum + isotope[0],
                                    intensity * isotope[1])
        contributions = []
        recursionLoop(0, 0, 1)
        return contributions

    def findArrangements(self, count, numIsotopes):
        
        arrangements = []
        stack = [(0, [], count)]

        while stack:
            index, current, remaining = stack.pop()

            if len(current) == numIsotopes:
                if remaining == 0:
                    arrangements.append(current)
                continue

            for i in range(0, remaining + 1):
                stack.append((index + 1, current + [i], remaining - i))

        return arrangements

    def determineIsotopeSpectra(self, formula, charge, threshold):

        isotopeSpectra = []
        
        for element, count in formula.items():

            elementContributions = [] # Probability dist of isotope combs for element.
            mn = multinomial(count, [i[1] for i in elements[element]])
            uniqueArrangements = self.findArrangements(count, len(elements[element]))

            for comb in uniqueArrangements:
                intensity = mn.pmf(comb)
            
                if intensity > threshold:
                    mass = sum(val * isotope[0] for val, isotope in zip(comb, elements[element]))
                    elementContributions.append([mass, intensity])
            
            isotopeSpectra.append(elementContributions)

        isotopeSpectra = self.sumContributions(isotopeSpectra)
        eMassDiff = -elements['e'][0][0]*charge # Difference in mass caused by charge
        normalisationFactor = max([i[1] for i in isotopeSpectra]) # Most abundant mass should be 100% abn.
        isotopeSpectra = [[round((i[0] + eMassDiff)/abs(charge), 6), 100*(i[1]/normalisationFactor)] for i in isotopeSpectra if i[1] > threshold]

        return isotopeSpectra
    
    def copyData(self):
            
            menu = QMenu()
            menu.addAction('Copy Value')
            menu.addSeparator()
            menu.addAction('Copy Entire Table')
            selection = menu.exec(QCursor.pos())

            try:

                if selection.text() == 'Copy Value':
                    selection = self.spectraTable.selectedIndexes()[0]
                    if selection.isValid():
                        cell_value = self.tableData.data(selection, Qt.DisplayRole)
                        clipboard = QApplication.clipboard()
                        clipboard.setText(str(cell_value))

                elif selection.text() == 'Copy Entire Table':
                    rows = self.tableData.rowCount()
                    columns = self.tableData.columnCount()
                    data = '\t'.join(['m/z (Da)', 'Abn. (%)']) + '\n'

                    for row in range(rows):
                        intensity = self.tableData.data(self.tableData.index(row, 1))
                        if intensity == '0':
                            continue
                        for col in range(columns):
                            index = self.tableData.index(row, col)
                            data += str(self.tableData.data(index)) + '\t'
                        data = data.strip() + '\n'

                    data = data.strip()

                    clipboard = QApplication.clipboard()
                    mime_data = QMimeData()
                    mime_data.setText(data)
                    clipboard.setMimeData(mime_data)

            except: pass # Invalid selection