from PySide6.QtGui import QIntValidator, QPixmap
from PySide6.QtWidgets import QCheckBox, QVBoxLayout, QLabel, QLineEdit, QWizard, QWizardPage

class Page(QWizardPage):
    '''
    First page. Cmin, cmax, dmin and dmax used to
    generate fatty acids. Later, to determine lipid spectra.
    '''
    def __init__(self, parent=None):
        super(Page, self).__init__(parent)

        self.setTitle("Fatty acid range")
        self.setSubTitle("Lipids will be generated according to the fatty acid range:\n"
                          "C min and C max determine chain lengths. "
                          "D min and D max determine the range of desaturation.")
        self.setPixmap(QWizard.WatermarkPixmap, QPixmap('Images\FAs.png'))
        self.vLayout = QVBoxLayout(self)

        # Text and inputs on Page 1
        # Input box for C min value, which determines the minimum fatty acid chain length
        # Fatty acids will be generated from C min to C max
        self.cminLabel = QLabel('C min:                              (1 -> 100)')
        self.cmin = QLineEdit()
        self.cmin.setValidator(QIntValidator(1, 100)) # Only allow ints between 1 - 100
        self.registerField('cmin*', self.cmin) # Register field to make it mandatory


        # Input box for C max value, which determines the maximum fatty acid chain length
        # Fatty acids will be generated from C min to C max
        self.cmaxLabel = QLabel('C max:                     (C min -> 100)')
        self.cmax = QLineEdit()
        self.cmax.setValidator(QIntValidator(1, 100))
        self.registerField('cmax*', self.cmax)

        # Input box for D min value, which determines the minimum desaturation of a fatty acid
        # Fatty acids of a given length will be generated from D min to D max where D < C/2
        self.dminLabel = QLabel('D min:               (0 -> 100 < C min)')
        self.dmin = QLineEdit()
        self.dmin.setValidator(QIntValidator(0, 100))
        self.registerField('dmin*', self.dmin)

        # Input box for D min value, which determines the maximum desaturation of a fatty acid
        # Fatty acids of a given length will be generated from D min to D max where D < C/2
        self.dmaxLabel = QLabel('D max:     (D min -> 100 < C max)')
        self.dmax = QLineEdit()
        self.dmax.setValidator(QIntValidator(0, 100))
        self.registerField('dmax*', self.dmax)



        self.omaxLabel = QLabel('O max:     (1 -> 10)')
        self.omax = QLineEdit()
        self.omax.setDisabled(True)
        self.omax.setValidator(QIntValidator(1, 10))
        self.registerField('omax', self.omax)
        self.omax.textEdited.connect(self.completeChanged)

        self.isomerism = QCheckBox('Respect sn isomerism', self)
        self.registerField('isomerism', self.isomerism)

        self.hydroxytickbox = QCheckBox('Include hydroxy-functionalised tails', self)
        self.hydroxytickbox.toggled.connect(self.omax.setEnabled)
        self.hydroxytickbox.toggled.connect(self.completeChanged)
        self.registerField('hydroxytickbox', self.hydroxytickbox)

        # Alternate text and inputs vertically
        self.vLayout.addWidget(self.cminLabel) # Place label
        self.vLayout.addWidget(self.cmin) # Place input box

        self.vLayout.addWidget(self.cmaxLabel)
        self.vLayout.addWidget(self.cmax)

        self.vLayout.addWidget(self.dminLabel)
        self.vLayout.addWidget(self.dmin)

        self.vLayout.addWidget(self.dmaxLabel)
        self.vLayout.addWidget(self.dmax)

        self.vLayout.addWidget(self.isomerism)
        self.vLayout.addWidget(self.hydroxytickbox)

        self.vLayout.addWidget(self.omaxLabel)
        self.vLayout.addWidget(self.omax)

    def isComplete(self):
        '''
        Restricts acceptable values for cmin, cmax, dmin and dmax as int.
        Values for c should be at minimum 1, arbitrary limit at 100
        Values for d should be from 0 to cmax.
        Restrictions also defined above (self.registerField()).
        '''
        if int(self.cmin.text() or 1) > int(self.cmax.text() or 0):
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
