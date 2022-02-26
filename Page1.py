from PySide6.QtGui import QIntValidator, QPixmap
from PySide6.QtWidgets import QRadioButton, QCheckBox, QVBoxLayout, QHBoxLayout, QLineEdit, QWizard, QWizardPage

class Page(QWizardPage):
    '''
    First page. Cmin, cmax, dmin and dmax used to
    generate fatty acids. Later, to determine lipid spectra.
    '''
    def __init__(self, parent=None):
        super(Page, self).__init__(parent)

        self.parent = parent

        self.setTitle("Fatty acid range")
        self.setSubTitle("Lipids will be generated according to the fatty acid range:\n"
                          "C min and C max determine chain lengths. "
                          "D min and D max determine the range of desaturation.")
        self.setPixmap(QWizard.WatermarkPixmap, QPixmap('Images\FAs.png'))
        self.vLayout = QVBoxLayout(self)
        self.hLayout = QHBoxLayout(self)

        # Determines whether a range of lipids will be generated, ie
        # from AA to ZZ, or whether specific lipids are to be generated.
        self.radButton1 = QRadioButton('Generate range of lipids')
        self.radButton2 = QRadioButton('Generate specific lipids')
        self.registerField('radButton2', self.radButton2)
        self.radButton1.setChecked(True)
        self.hLayout.addWidget(self.radButton1)
        self.hLayout.addWidget(self.radButton2)
        self.vLayout.addLayout(self.hLayout)

        # Optional tickbox to consider every fatty acid combination, ie AA, AB, BA, BB, ...
        # instead of only the unique combinations, ie AA, AB, BB, ...
        # Generation takes a very long time when activated !!
        self.isomerism = QCheckBox('Respect sn isomerism', self)
        self.registerField('isomerism', self.isomerism)
        self.vLayout.addWidget(self.isomerism)

        # Text and inputs on Page 1
        # Input box for C min value, which determines the minimum fatty acid chain length
        # Fatty acids will be generated from C min to C max
        self.cmin = QLineEdit()
        self.cmin.setPlaceholderText(' C min:                                                                                                                               (1 -> 40)')
        self.cmin.setValidator(QIntValidator(1, 40)) # Only allow ints between 1 - 40
        self.registerField('cmin*', self.cmin) # Register field to make it mandatory
        self.vLayout.addWidget(self.cmin)

        # Input box for C max value, which determines the maximum fatty acid chain length
        # Fatty acids will be generated from C min to C max
        self.cmax = QLineEdit()
        self.cmax.setPlaceholderText(' C max:                                                                                                                     (C min -> 40)')
        self.cmax.setValidator(QIntValidator(1, 40))
        self.registerField('cmax*', self.cmax)
        self.vLayout.addWidget(self.cmax)

        # Input box for D min value, which determines the minimum desaturation of a fatty acid
        # Fatty acids of a given length will be generated from D min to D max where D < C/2
        self.dmin = QLineEdit()
        self.dmin.setPlaceholderText(' D min:                                                                                                               (0 -> 40 < C min)')
        self.dmin.setValidator(QIntValidator(0, 40))
        self.registerField('dmin*', self.dmin)
        self.vLayout.addWidget(self.dmin)

        # Input box for D min value, which determines the maximum desaturation of a fatty acid
        # Fatty acids of a given length will be generated from D min to D max where D < C/2
        self.dmax = QLineEdit()     
        self.dmax.setPlaceholderText(' D max:                                                                                                     (D min -> 40 < C max)')
        self.dmax.setValidator(QIntValidator(0, 40))
        self.registerField('dmax*', self.dmax)
        self.vLayout.addWidget(self.dmax)

        # Optional input box for O max value, which determines maximum number of hydroxy groups
        # to add to the fatty acid. By default 0.
        self.hydroxytickbox = QCheckBox('Include hydroxy-functionalised tails', self)
        self.registerField('hydroxytickbox', self.hydroxytickbox)
        self.vLayout.addWidget(self.hydroxytickbox)
        self.omax = QLineEdit()
        self.omax.setPlaceholderText(' O max:                                                                                                                               (1 -> 10)')
        self.omax.setDisabled(True)
        self.omax.setValidator(QIntValidator(1, 10))
        self.registerField('omax', self.omax)
        self.hydroxytickbox.toggled.connect(self.omax.setEnabled)
        self.hydroxytickbox.toggled.connect(self.omax.clear)
        self.hydroxytickbox.toggled.connect(self.completeChanged)
        self.omax.textEdited.connect(self.completeChanged)
        self.vLayout.addWidget(self.omax)

        # Toggling radiobuttons 1 or 2 will toggle all entry fields
        self.radButton1.toggled.connect(self.radButton1Enabled)
        self.radButton2.toggled.connect(self.radButton2Enabled)
        self.radButton2.toggled.connect(self.completeChanged)
        # Enable everything for button 1
    def radButton1Enabled(self):
        self.cmin.setDisabled(False)
        self.cmax.setDisabled(False)
        self.dmin.setDisabled(False)
        self.dmax.setDisabled(False)
        self.hydroxytickbox.setDisabled(False)
        # Disable everything for button 2
    def radButton2Enabled(self):
        self.cmin.setDisabled(True)
        self.cmin.clear()
        self.cmax.setDisabled(True)
        self.cmax.clear()
        self.dmin.setDisabled(True)
        self.dmin.clear()
        self.dmax.setDisabled(True)
        self.dmax.clear()
        self.hydroxytickbox.setChecked(False)
        self.hydroxytickbox.setDisabled(True)

    def isComplete(self):
        '''
        Restricts acceptable values for cmin, cmax, dmin and dmax as int.
        Values for c should be at minimum 1, arbitrary limit at 100
        Values for d should be from 0 to cmax.
        Restrictions also defined above (self.registerField()).
        '''
        if self.radButton2.isChecked():
            return True
        elif int(self.cmin.text() or 1) > int(self.cmax.text() or 0):
            return False
        elif int(self.dmin.text() or 1) > int(self.dmax.text() or 0):         
            return False
        elif int(self.dmin.text() or 0) > int(self.cmin.text() or 0)-1:
            return False
        elif int(self.dmax.text() or 0) > int(self.cmax.text() or 0)-1:
            return False
        elif self.hydroxytickbox.isChecked() and int(self.omax.text() or 0) not in range(1, 11):
            return False  
        return super().isComplete()

    def nextId(self):
        if self.radButton2.isChecked():
            return 2 # Page 2B
        else:
            return 1 # Page 2

    
