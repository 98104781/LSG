import sys
from PySide6.QtWidgets import QApplication, QWizard

import Page0
import Page1
import Page2
import Page3
import Page4
import Page5

class CreateWindow(QWizard):
    '''
    Wizard to guide through spectra generation
    '''
    def __init__(self):
        super().__init__()

        self.setWizardStyle(QWizard.ModernStyle)
        self.setWindowTitle('LSG3')
        self.setFixedSize(600, 510)

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