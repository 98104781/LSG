import Wizard.ResourcePath as RP
from PySide6.QtGui import QIntValidator, QPixmap
from PySide6.QtWidgets import QCheckBox, QVBoxLayout, QHBoxLayout, QLineEdit, QWizard, QWizardPage

class Page(QWizardPage):
    '''
    First page. Cmin, cmax, dmin and dmax used to
    generate fatty acids. Later, to determine lipid spectra.
    '''
    def __init__(self, parent=None):
        super(Page, self).__init__(parent)

        self.parent = parent

        self.setTitle("Generate a range of lipids using a range of tails")
        self.setSubTitle("Please define the limits of the range to use.   Be aware, large ranges can produce large libraries!\n"
                         "C min and C max determine chain lengths.   "
                         "D min and D max determine the range of desaturation.")
        image_Path = RP.resource_path('Images\FAs.png')
        self.setPixmap(QWizard.WatermarkPixmap, QPixmap(image_Path))
        self.vLayout = QVBoxLayout(self)
        self.hLayout = QHBoxLayout(self)

        # Optional tickbox to consider every fatty acid combination, ie AA, AB, BA, BB, ...
        # instead of only the unique combinations, ie AA, AB, BB, ...
        # Generation takes a very long time when activated !!
        self.isomerism = QCheckBox('Respect sn isomerism', self)
        self.registerField('isomerism', self.isomerism)
        self.isomerism.hide()
        #self.vLayout.addWidget(self.isomerism) # Lipid isomers incomplete

        # Text and inputs on Page 1
        # Input box for C min value, which determines the minimum fatty acid chain length
        # Fatty acids will be generated from C min to C max
        self.cmin = QLineEdit()
        self.cmin.setPlaceholderText(' C min:'+86*' '+'(2 -> C max )')
        self.cmin.setValidator(QIntValidator(2, 30)) # Only allow ints between 1 - 40
        self.registerField('cmin*', self.cmin) # Register field to make it mandatory
        self.vLayout.addWidget(self.cmin)

        # Input box for C max value, which determines the maximum fatty acid chain length
        # Fatty acids will be generated from C min to C max
        self.cmax = QLineEdit()
        self.cmax.setPlaceholderText(' C max:'+85*' '+'(C min -> 30)')
        self.cmax.setValidator(QIntValidator(1, 30))
        self.registerField('cmax*', self.cmax)
        self.vLayout.addWidget(self.cmax)

        # Input box for D min value, which determines the minimum desaturation of a fatty acid
        # Fatty acids of a given length will be generated from D min to D max where D < C/2
        self.dmin = QLineEdit()
        self.dmin.setPlaceholderText(' D min:'+86*' '+'(0 -> D max )')# < C min)')
        self.dmin.setValidator(QIntValidator(0, 12))
        self.registerField('dmin*', self.dmin)
        self.vLayout.addWidget(self.dmin)

        # Input box for D min value, which determines the maximum desaturation of a fatty acid
        # Fatty acids of a given length will be generated from D min to D max where D < C/2
        self.dmax = QLineEdit()     
        self.dmax.setPlaceholderText(' D max:'+85*' '+'(D min -> 12)')# < C max)')
        self.dmax.setValidator(QIntValidator(0, 12))
        self.registerField('dmax*', self.dmax)
        self.vLayout.addWidget(self.dmax)

        # Optional input box for O max value, which determines maximum number of hydroxy groups
        # to add to the fatty acid. By default 0.
        self.hydroxytickbox = QCheckBox('Include hydroxy-functionalised tails', self)
        self.registerField('hydroxytickbox', self.hydroxytickbox)
        self.vLayout.addWidget(self.hydroxytickbox)
        self.omax = QLineEdit()
        self.omax.setPlaceholderText(' O max:'+85*' '+'(1 -> 5)')
        self.omax.setDisabled(True)
        self.omax.setValidator(QIntValidator(1, 8))
        self.registerField('omax', self.omax)
        self.hydroxytickbox.toggled.connect(self.omax.setEnabled)
        self.hydroxytickbox.toggled.connect(self.omax.clear)
        self.hydroxytickbox.toggled.connect(self.completeChanged)
        self.omax.textEdited.connect(self.completeChanged)
        self.vLayout.addWidget(self.omax)

        self.ceramideVariability = QCheckBox('Vary ceramide base length with fatty acids', self)
        self.registerField('ceramideVariability', self.ceramideVariability)
        self.vLayout.addWidget(self.ceramideVariability)


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
        elif self.hydroxytickbox.isChecked() and int(self.omax.text() or 0) not in range(1, 6):
            return False  
        return super().isComplete()

    def nextId(self):
        return 3

    
