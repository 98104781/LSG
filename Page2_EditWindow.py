
import GenerateLipids as GL
import Spectra

from itertools import combinations_with_replacement as cwr
from PySide6.QtCore import QModelIndex
from PySide6.QtWidgets import QComboBox, QDialog, QHeaderView, QTableView, QPushButton, QVBoxLayout, QHBoxLayout, QLabel

class NewWindow(QDialog): # Opened from SpectraSetupPage
    '''
    Popup window to modify selected class/adduct spectra.
    '''
    def __init__(self, treeData):
        super().__init__()

        self.setWindowTitle('LSG3')
        self.setFixedSize(600, 510)
        self.vLayout = QVBoxLayout(self)
        self.vLayout2 = QVBoxLayout(self)
        self.hLayout = QHBoxLayout(self)

        # Example spectra will be GPL 16:0_18:1
        self.tails = [GL.sn(c=16, d=0, type='Acyl'),
                      GL.sn(c=18, d=1, type='Acyl')]

        self.tableView = QTableView()
        self.comboBox = QComboBox()
        self.spectra = Spectra.SpectraScatter()
        self.spectra.setFixedHeight(200)
        self.spectra.setFixedWidth(400)
        self.label = QLabel("If generated spectra don't respect sn isomerism, all isomer dependent\n"
                            "fragments will be equal in intensity.\n\nIf the observed fragment intensities"
                            " differ from the default provided, they may\nbe manually updated on the left.")
        self.comboBox.currentTextChanged.connect(self.buildList) # Lipid data will be stored in comboboxs
        for lipidClass, adductList in treeData.items():  # Unpack lipid/adduct data from Treeview for use
            for adduct in adductList: # Unpack all adducts from adduct list for use in separate comboboxs
                self.comboBox.addItem(lipidClass.text(0)+' '+adduct.text(0), [lipidClass, adduct]) # here

        self.applyBTN = QPushButton('Apply Changes')
        self.applyBTN.clicked.connect(self.applyChanges)
        #self.resetBTN = QPushButton('Reset Changes')
        #self.applyBTN.clicked.connect(self.resetChanges)

        self.vLayout.addWidget(self.comboBox)
        self.hLayout.addWidget(self.tableView)

        self.vLayout2.addWidget(self.spectra)
        self.vLayout2.addWidget(self.label)
        self.vLayout2.addStretch()
        self.vLayout2.addWidget(self.applyBTN)
        self.hLayout.addLayout(self.vLayout2)

        self.vLayout.addLayout(self.hLayout)

    def buildLipid(self, data):
        '''
        Returns an example lipid for use.
        Lipid will be GPL (16:0_) 16:0_18:1.
        '''
        GL.Glycerolipid.instances = []
        cls = data[0].lipidClass # Example lipid
        _, comb, *_ = cwr(self.tails, cls.No_Tails)
        example = cls(*comb) # comb will be (16:0_) 16:0_18:1
        example.resolve_spectra(data[1].text(0), data[1].fragmentList)
        return example 
    
    def buildList(self):
        '''
        Updates selected data,
        Refreshes fragment list and spectra.
        '''
        data = self.comboBox.currentData()
        if data is not None:
            example = self.buildLipid(data)
            adduct = data[1].text(0)
            fragments = example.spectra[adduct]
            mz = GL.MA(example, adduct, 0).mass

            self.spectra.setSpectra(example.name+' '+adduct, mz, fragments)
            self.table = Spectra.SpectraTableModel(fragments)
            self.tableView.setModel(self.table)
            self.tableView.setItemDelegateForColumn(1, Spectra.SpinBoxDelegate(self.tableView))
            self.tableView.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
            self.tableView.verticalHeader().hide()
            self.table.dataChanged.connect(self.updateData)
        else: pass

    def updateData(self, index: QModelIndex, *other):
        '''
        Updates adduct spectra in list
        with new intensity.
        '''
        row = int(index.row())
        fragmentList = self.comboBox.currentData()[1].fragmentList
        data = self.table.tdata[row]
        fragmentList[data.fragmentType] = data.intensity
        self.comboBox.currentData()[1].fragmentList = fragmentList
        self.buildList()

    def applyChanges(self):
        '''
        Overwrites adduct spectra in memory
        with new adduct spectra.
        '''
        for x in range(self.comboBox.count()):
            data = self.comboBox.itemData(x)
            lipidClass = data[0].lipidClass
            adduct = data[1].text(0)
            fragmentList = data[1].fragmentList
            for fragment in list(fragmentList):
                if fragmentList[fragment] == 0:
                    fragmentList.pop(fragment)
            lipidClass.adducts[adduct] = fragmentList
        self.buildList()

    def resetChanges(self):
        pass
