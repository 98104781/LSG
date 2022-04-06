import Wizard.EditLipidAdduct as LAEW

from PySide6.QtGui import QPixmap
from PySide6.QtWidgets import QRadioButton, QButtonGroup, QVBoxLayout, QHBoxLayout, QLabel, QWizard, QWizardPage, QPushButton

class Page(QWizardPage):
    '''
    First page. Cmin, cmax, dmin and dmax used to
    generate fatty acids. Later, to determine lipid spectra.
    '''
    def __init__(self, parent=None):
        super(Page, self).__init__(parent)

        self.parent = parent

        self.setTitle(" ")
        self.setSubTitle(" "
                         " "
                         " ")
        self.setPixmap(QWizard.WatermarkPixmap, QPixmap('Images\LSG.png'))
        self.vLayout =  QVBoxLayout(self)
        self.hLayout =  QHBoxLayout(self)
        self.hLayout2 = QHBoxLayout(self)

        # Determines whether a range of lipids will be generated, ie
        # from AA to ZZ, or whether specific lipids are to be generated.
        self.textBox1 = QLabel('Lipid data can be generated in one of two ways:\n')
        self.buttonGroup1 = QButtonGroup()
        self.lipidRange = QRadioButton('Generate range of lipids')
        self.buttonGroup1.addButton(self.lipidRange)
        self.lipidRange.setChecked(True)
        self.lipidRange.toggled.connect(self.enableOptions)
        self.lipidSpecific = QRadioButton('Generate specific lipids')
        self.buttonGroup1.addButton(self.lipidSpecific)
        self.lipidSpecific.toggled.connect(self.disableOptions)
        self.registerField('lipidSpecific', self.lipidSpecific)
 
        self.textBox2 = QLabel('\nGenerating a range of lipids will produce all combinations (without reptition)\n'
                               'using either a generated range of lipid tails, or tails that are specifically defined.\n\n')
        self.textBox3 = QLabel('\nGenerating specific lipids will produce only those lipids specifically defined.')

        self.buttonGroup2 = QButtonGroup()
        self.tailRange = QRadioButton('Using a range of tails')
        self.tailRange.setChecked(True)
        self.tailRange.clicked.connect(self.completeChanged)
        self.buttonGroup2.addButton(self.tailRange)
        self.tailSpecific = QRadioButton('Using specific tails')
        self.registerField('tailSpecific', self.tailSpecific)
        #self.tailSpecific.setDisabled(True) # Disabled for now, incomplete!
        self.tailSpecific.clicked.connect(self.completeChanged)
        self.buttonGroup2.addButton(self.tailSpecific)

        self.editLipids = QPushButton('Edit Lipid Classes')
        self.editLipids.clicked.connect(self.editLipidClasses)

        self.vLayout.addWidget(self.textBox1)
        self.hLayout.addWidget(self.lipidRange)
        self.hLayout.addWidget(self.lipidSpecific)
        self.vLayout.addLayout(self.hLayout)
        self.vLayout.addWidget(self.textBox2)
        self.hLayout2.addWidget(self.tailRange)
        self.hLayout2.addWidget(self.tailSpecific)
        self.vLayout.addLayout(self.hLayout2)
        self.vLayout.addWidget(self.textBox3)
        self.vLayout.addStretch()
        self.vLayout.addWidget(self.editLipids)

    def enableOptions(self):
        self.tailSpecific.setDisabled(False)
        self.tailRange.setDisabled(False)
        self.completeChanged.emit()

    def disableOptions(self):
        self.tailSpecific.setDisabled(True)
        self.tailSpecific.setChecked(False)
        self.tailRange.setDisabled(True)
        self.tailRange.setChecked(False)
        self.completeChanged.emit()

    def isComplete(self):

        if self.lipidRange.isChecked():
            if self.tailRange.isChecked() is False and self.tailSpecific.isChecked() is False:
                return False

        return super().isComplete()

    def nextId(self):
        if self.lipidSpecific.isChecked():
            return 4
        elif self.tailRange.isChecked():
            return 1
        elif self.tailSpecific.isChecked():
            return 2

    def editLipidClasses(self):
        classes_to_generate = self.parent.classes_to_generate
        editLipidWindow = LAEW.NewWindow(self, classes_to_generate)
        editLipidWindow.exec()
