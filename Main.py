import sys
from PySide6.QtWidgets import QApplication, QWizard

import Page1
import Page2
import Page2B
import Page3

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
        self.setPage(0, Page1.Page(self)) # Define range for lipid tails   or   choose to generate specific lipids.
        self.setPage(1, Page2.Page(self)) # Define lipid classes and spectra (if range is chosen).
        self.setPage(2, Page2B.Page(self))# Define which specific lipids to generate (if specific is chosen).
        self.setPage(3, Page3.Page(self)) # Generate the specified lipids.

if __name__ == '__main__':
    app = QApplication(sys.argv)
    mainWindow = CreateWindow()
    mainWindow.show()
    sys.exit(app.exec())