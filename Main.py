import sys
import psutil
import inspect
import Lipids.Classes as Classes
from Lipids.GenerateLipids import Glycerolipid, OtherLipid, Sphingolipid

from PySide6.QtCore import QTimer
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
        self.setOption(QWizard.HaveCustomButton1, True)
        self.setButtonText(QWizard.CustomButton1, "RAM Usage: 0 MB")
        self.custom_button1 = self.button(QWizard.CustomButton1)
        self.custom_button1.setEnabled(False)
        self.custom_button1.setStyleSheet("""
            QPushButton[enabled="false"] {
                background-color: #ffffff;
                color: #000000;
                border: 1px solid #d0d0d0;            
                min-width: 120px;
            }
        """)

        self.setOption(QWizard.HaveCustomButton2, True)
        self.setButtonText(QWizard.CustomButton2, "CPU Usage: 0 %")
        self.custom_button2 = self.button(QWizard.CustomButton2)
        self.custom_button2.setEnabled(False)
        self.custom_button2.setStyleSheet("""
            QPushButton[enabled="false"] {
                background-color: #ffffff;
                color: #000000;
                border: 1px solid #d0d0d0;            
                min-width: 120px;
            }
        """)

        self.timer = QTimer(self)
        self.timer.timeout.connect(self.update_ram_usage)
        self.timer.timeout.connect(self.update_cpu_usage)
        self.timer.start(1000)  # Update every second

        # Glycerolipids
        self.classes_to_generate =      [cls for cls in Glycerolipid.__subclasses__() if inspect.getmodule(cls) == Classes]
        # Sphingolipids
        self.classes_to_generate.extend([cls for cls in Sphingolipid.__subclasses__() if inspect.getmodule(cls) == Classes])
        # ETC Lipids, Cholesterol ester
        self.classes_to_generate.extend([cls for cls in OtherLipid.__subclasses__()   if inspect.getmodule(cls) == Classes])

        #print('Classes: ', len(self.classes_to_generate)) # Used to check number of classes available in console

        #x = 0
        #for cls in self.classes_to_generate: # Used to check number of adducts available in console
        #    x += len(cls.adducts)
        #print('Adducts: ', x)

        # Add Wizard Pages
        self.setPage(0, Page0.Page(self)) # Define range for lipid tails   or   choose to generate specific lipids.
        self.setPage(1, Page1.Page(self))
        self.setPage(2, Page2.Page(self))
        self.setPage(3, Page3.Page(self)) # Define lipid classes and spectra (if range is chosen).
        self.setPage(4, Page4.Page(self)) # Define which specific lipids to generate (if specific is chosen).
        self.setPage(5, Page5.Page(self)) # Generate the specified lipids.

    def update_ram_usage(self):
        ram_usage = psutil.Process().memory_info().rss / (1024 ** 2)  # RAM usage in MB
        self.setButtonText(QWizard.CustomButton1, f"RAM Usage: {ram_usage:.0f} MB")

    def update_cpu_usage(self):
        cpu_usage = psutil.cpu_percent()
        self.setButtonText(QWizard.CustomButton2, f"CPU Usage: {cpu_usage:.2f}%")

if __name__ == '__main__':
    app = QApplication(sys.argv)
    mainWindow = CreateWindow()
    mainWindow.show()
    sys.exit(app.exec())
