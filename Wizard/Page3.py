import Wizard.ResourcePath as RP
import Wizard.EditLipidAdduct as LAEW

from PySide6.QtGui import QPixmap, QCursor
from PySide6.QtCore import Property, Qt, Signal
from PySide6.QtWidgets import QPushButton, QCheckBox, QTreeWidget, QTreeWidgetItem, QVBoxLayout, QWizard, QWizardPage, QMenu

class Page(QWizardPage):
    '''
    Displays currently supported lipid classes
    with their currently supported adducts for generation.
    '''
    treeDataChanged = Signal()

    def __init__(self, parent=None):
        super(Page, self).__init__(parent)

        self.wiz = parent

        self.setTitle("Select lipid classes to generate")
        self.setSubTitle("Select from the list of available lipid classes below.\n"
                         "Spectra will be generated for the selected classes using the tails previously defined.")
        image_Path = RP.resource_path('Images\GPLs.png')
        self.setPixmap(QWizard.WatermarkPixmap, QPixmap(image_Path))
        self.setCommitPage(True)
        self.vLayout = QVBoxLayout(self)

        self.classQbox = {} # GPL QCheckBoxes        
        self.classAdductQbox = {} # Adduct QCheckBoxes

        self.specificOrganisation = QCheckBox('Respect headgroup-acyl isomerism')
        self.specificOrganisation.setToolTip('Some lipids contain fatty-acyls on headgroup, e.g. AC3PIM2. Duplicate fatty-acyl'
                                             '\ncombinations are required to generate all species e.g. 16:0/18:0_18:1 and 18:1/16:0_18:1.'
                                             '\nHowever this can lead to excessively large libraries. If unticked, only unique'
                                             '\nand location-unspecific combinations generated e.g. 16:0_18:0_18:1. However, ensure'
                                             '\nfatty-acyl location-specific fragments are excluded in the spectra selected below!')
        self.specificOrganisation.setChecked(True)
        self.registerField('specificOrganisation', self.specificOrganisation)

        self.treeView = QTreeWidget()
        self.treeView.setContextMenuPolicy(Qt.CustomContextMenu)
        self.treeView.customContextMenuRequested.connect(self.editLipidContextMenu)
        self.treeView.doubleClicked.connect(self.editLipid)
        self.treeView.setHeaderHidden(True)
        self.registerField("tree", self, "tree_property")
        
        self.modifybutton = QPushButton("Modify Lipid")
        self.modifybutton.clicked.connect(self.editLipid)

        self.vLayout.addWidget(self.specificOrganisation)
        self.vLayout.addWidget(self.treeView)
        self.vLayout.addWidget(self.modifybutton)

    def editLipidContextMenu(self):
        menu = QMenu()
        menu.addAction("Edit Lipid")
        selection = menu.exec(QCursor.pos())
        if selection:
            self.editLipid()

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

    def initializePage(self) -> None:

        self.treeView.clear()

        self.classList = self.wiz.classes_to_generate

        for cls in self.classList: #  Make boxes for Treeview
            self.classQbox[cls]  =  QTreeWidgetItem(self.treeView)
            root = self.classQbox[cls] # Creates tickbox for class
            root.lipidClass = cls # Custom variable to store class
            root.setText(0, cls.__name__) # Gives name for tickbox
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
                try:child.setCheckState(0, self.aBoxStatus[cls][adduct])
                except:child.setCheckState(0, Qt.Unchecked) # Tickstate
                child.setFlags(child.flags() | Qt.ItemIsUserCheckable) 

        return super().initializePage()
            
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

    def nextId(self):
        return 5