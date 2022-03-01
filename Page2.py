import inspect

import Classes
import Classes_isomers
from GenerateLipids import Glycerolipid, Sphingolipid
import Page2_EditWindow as P2EW

from PySide6.QtCore import Property, Qt, Signal
from PySide6.QtGui import QPixmap
from PySide6.QtWidgets import QPushButton, QTreeWidget, QTreeWidgetItem, QVBoxLayout, QWizard, QWizardPage

class Page(QWizardPage):
    '''
    Second page. Displays currently supported lipid classes
    with their currently supported adducts for generation.
    '''
    treeDataChanged = Signal()

    def __init__(self, parent=None):
        super(Page, self).__init__(parent)

        self.setTitle("Select lipid classes")
        self.setSubTitle("Select from the list of available glycerophospholipid classes below.\n"
                         "Spectra will be generated for the selected classes using the fatty acid ranges.")
        self.setPixmap(QWizard.WatermarkPixmap, QPixmap('Images\GPLs.png'))
        self.vLayout = QVBoxLayout(self)

        self.classQbox = {} # GPL QCheckBoxes        
        self.classAdductQbox = {} # Adduct QCheckBoxes

        self.treeView = QTreeWidget()
        self.treeView.setHeaderHidden(True)
        self.registerField("tree", self, "tree_property")

        self.modifybutton = QPushButton("Modify fragmentation spectra for selected adduct(s)")
        self.modifybutton.clicked.connect(self.open_editspectrawindow)

        self.vLayout.addWidget(self.treeView)
        self.vLayout.addWidget(self.modifybutton)
     

    def initializePage(self) -> None:

        self.treeView.clear()

        if self.field('isomerism') == False:
            classes_to_generate = [cls for cls in Glycerolipid.__subclasses__() if inspect.getmodule(cls) == Classes]
            classes_to_generate.extend([cls for cls in Sphingolipid.__subclasses__() if inspect.getmodule(cls) == Classes])
        else:
            classes_to_generate = [cls for cls in Glycerolipid.__subclasses__() if inspect.getmodule(cls) == Classes_isomers]

        for cls in classes_to_generate: #  Make boxes for Treeview
            self.classQbox[cls]  =  QTreeWidgetItem(self.treeView)
            root = self.classQbox[cls] # Creates tickbox for class
            root.lipidClass = cls # Custom variable to store class
            root.setText(0, cls.__name__) # Gives name for tickbox
            root.setCheckState(0, Qt.Unchecked) #   Untick tickbox
            root.setFlags(root.flags() | Qt.ItemIsAutoTristate | Qt.ItemIsUserCheckable)
            self.classAdductQbox[cls] = {} # Open dict for adducts

            for adduct in cls.adducts: # Make sub-box for treeview
                self.classAdductQbox[cls][adduct] = QTreeWidgetItem(self.classQbox[cls])
                child = self.classAdductQbox[cls][adduct] # Assign
                child.fragmentList = cls.adducts[adduct] # Adducts
                child.setText(0, adduct) #  Gives name for tickbox
                child.setCheckState(0, Qt.Unchecked) #  Untick box
                child.setFlags(child.flags() | Qt.ItemIsUserCheckable) 

        return super().initializePage()

    def open_editspectrawindow(self):
        '''
        Opens external window to modify selected spectra.
        Should display an example lipid spectra of GPL
        16:0_18:1 with appropriate masses and intensity.
        '''
        editspectrawindow = P2EW.NewWindow(self.field('tree'))
        editspectrawindow.exec()
            
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
        return 3 # Page 3