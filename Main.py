import sys
import inspect
import Lipids.Classes as Classes
from Lipids.GenerateLipids import Glycerolipid, OtherLipid, Sphingolipid

from PySide6.QtWidgets import QApplication, QWizard

import Wizard.Page0 as Page0
import Wizard.Page1 as Page1
import Wizard.Page2 as Page2
import Wizard.Page3 as Page3
import Wizard.Page4 as Page4
import Wizard.Page5 as Page5

class CreateWindow(QWizard):
    '''
    Wizard to guide through spectra generation
    '''
    def __init__(self):
        super().__init__()

        self.setWizardStyle(QWizard.ModernStyle)
        self.setWindowTitle('Lipid Spectrum Generator')
        self.setFixedSize(600, 510)

        # Glycerolipids
        self.classes_to_generate =      [cls for cls in Glycerolipid.__subclasses__() if inspect.getmodule(cls) == Classes]
        # Sphingolipids
        self.classes_to_generate.extend([cls for cls in Sphingolipid.__subclasses__() if inspect.getmodule(cls) == Classes])
        # ETC Lipids, Cholesterol ester
        self.classes_to_generate.extend([cls for cls in OtherLipid.__subclasses__()   if inspect.getmodule(cls) == Classes])

        #print(len(self.classes_to_generate)) # Used to check number of classes available in console

        # Add Wizard Pages
        self.setPage(0, Page0.Page(self)) # Define range for lipid tails   or   choose to generate specific lipids.
        self.setPage(1, Page1.Page(self))
        self.setPage(2, Page2.Page(self))
        self.setPage(3, Page3.Page(self)) # Define lipid classes and spectra (if range is chosen).
        self.setPage(4, Page4.Page(self)) # Define which specific lipids to generate (if specific is chosen).
        self.setPage(5, Page5.Page(self)) # Generate the specified lipids.

if __name__ == '__main__':
    app = QApplication(sys.argv)
    mainWindow = CreateWindow()
    mainWindow.show()
    sys.exit(app.exec())