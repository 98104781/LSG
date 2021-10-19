import os
import csv
import sys
import time
import Classes
import GenerateLipids as GL
from PySide6.QtWidgets import QFileDialog, QProgressBar
from itertools import combinations_with_replacement as cwr

from PySide6.QtCore import Property, Qt, Signal
from PySide6.QtGui import QIntValidator, QPixmap
from PySide6.QtWidgets import QApplication, QPlainTextEdit, QPushButton, QTreeWidget, QTreeWidgetItem, QVBoxLayout, QLabel, QLineEdit, QWizard, QWizardPage

gplClassList = [cls for cls in GL.Glycerolipid.__subclasses__()]




class CreateWindow(QWizard):

    def __init__(self):
        super().__init__()

        self.setWizardStyle(QWizard.ModernStyle)
        self.setWindowTitle('LSG3')
        self.setFixedSize(600, 510)

        # Add Wizard Pages
        self.addPage(Page1(self))
        self.addPage(Page2(self))
        self.addPage(Page3(self))




class Page1(QWizardPage):
    def __init__(self, parent=None):
        super(Page1, self).__init__(parent)

        self.setTitle("Fatty acid range")
        self.setSubTitle("Lipids will be generated according to the fatty acid range:\n"
                          "C min and C max determine chain lengths. "
                          "D min and D max determine the range of desaturation.")
        self.setPixmap(QWizard.WatermarkPixmap, QPixmap('Images\FAs.png'))
        self.vLayout = QVBoxLayout(self)

        # Text and inputs on Page 1
        self.cminLabel = QLabel('C min:                              (1 -> 100)') # Create label
        self.cmin = QLineEdit() # Create input box
        self.cmin.setValidator(QIntValidator(1, 100)) # Only allow ints between 1 - 100
        self.registerField('cmin*', self.cmin) # Register field to make it mandatory

        self.cmaxLabel = QLabel('C max:                     (C min -> 100)')
        self.cmax = QLineEdit()
        self.cmax.setValidator(QIntValidator(1, 100))
        self.registerField('cmax*', self.cmax)

        self.dminLabel = QLabel('D min:               (0 -> 100 < C min)')
        self.dmin = QLineEdit()
        self.dmin.setValidator(QIntValidator(0, 100))
        self.registerField('dmin*', self.dmin)

        self.dmaxLabel = QLabel('D max:     (D min -> 100 < C max)')
        self.dmax = QLineEdit()
        self.dmax.setValidator(QIntValidator(0, 100))
        self.registerField('dmax*', self.dmax)

        # Alternate text and inputs vertically
        self.vLayout.addWidget(self.cminLabel) # Place label
        self.vLayout.addWidget(self.cmin) # Place input box

        self.vLayout.addWidget(self.cmaxLabel)
        self.vLayout.addWidget(self.cmax)

        self.vLayout.addWidget(self.dminLabel)
        self.vLayout.addWidget(self.dmin)

        self.vLayout.addWidget(self.dmaxLabel)
        self.vLayout.addWidget(self.dmax)

    def isComplete(self):
        if int(self.cmin.text() or 0) > int(self.cmax.text() or 0):
            return False
        elif int(self.dmin.text() or 0) > int(self.dmax.text() or 0):
            return False
        elif int(self.dmin.text() or 0) > int(self.cmin.text() or 0)-1:
            return False
        elif int(self.dmax.text() or 0) > int(self.cmax.text() or 0)-1:
            return False
        return super().isComplete()




class Page2(QWizardPage):
    treeDataChanged = Signal()

    def __init__(self, parent=None):
        super(Page2, self).__init__(parent)

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
        
        for cls in gplClassList:
            self.classQbox[cls] = QTreeWidgetItem(self.treeView)
            root = self.classQbox[cls]
            # Custom variable to store lipid class
            root.lipidClass = cls
            root.setText(0, cls.__name__)
            root.setCheckState(0, Qt.Unchecked)
            root.setFlags(root.flags() | Qt.ItemIsAutoTristate | Qt.ItemIsUserCheckable)
            self.classAdductQbox[cls] = {} # Open class for adducts to add to

            for adduct in cls.adducts:
                self.classAdductQbox[cls][adduct] = QTreeWidgetItem(self.classQbox[cls])
                child = self.classAdductQbox[cls][adduct]
                # Custom variable to store fragment list
                child.fragmentList = cls.adducts[adduct]
                child.setText(0, adduct)
                child.setCheckState(0, Qt.Unchecked)
                child.setFlags(child.flags() | Qt.ItemIsUserCheckable)
                
        self.vLayout.addWidget(self.treeView)

    def treeData(self):
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




class Page3(QWizardPage):
    def __init__(self, parent=None):
        super(Page3, self).__init__(parent)

        self.setTitle("Select filetype to generate")
        self.setSubTitle("Press 'Generate' to create file \n"
                         "Currently .msp spectral library and .csv QE+ inclusion list supported")
        self.setPixmap(QWizard.WatermarkPixmap, QPixmap('Images\ADs.png'))
        self.vLayout = QVBoxLayout(self)

        self.generatebutton = QPushButton("Generate")
        self.generatebutton.clicked.connect(self.save_as)

        self.vLayout.addWidget(self.generatebutton)

        self.output_console = QPlainTextEdit()
        self.output_console.setReadOnly(True)
        self.vLayout.addWidget(self.output_console)

        self.progress_bar = QProgressBar()
        self.vLayout.addWidget(self.progress_bar)
    
    def initializePage(self) -> None:

        self.selected_class_adducts = self.field('tree')

        self.output_console.clear()
        self.output_console.appendPlainText('Tails will be generated from '
                                +self.field('cmin')+':'+self.field('dmin')+
                          ' -> '+self.field('cmax')+':'+self.field('dmax')+
                             '\nThe following classes will be generated:')
        
        caString = '' # Generates string, and class list
        self.classes_to_generate = []
        for item, item2 in self.selected_class_adducts.items():

            caString += '- '+item.text(0)+'\n'
            self.classes_to_generate.append(item.lipidClass)

            self.adducts_to_generate = {}
            for adduct in item2: # Update the selected adducts
                self.adducts_to_generate.update({adduct.text(0):adduct.fragmentList})
            item.lipidClass.adducts = self.adducts_to_generate

        self.output_console.appendPlainText(caString)
        self.tails_to_generate = [int(self.field('cmin')), int(self.field('cmax')),
                                  int(self.field('dmin')), int(self.field('dmax'))]
        
        self.progress_bar.reset()

        return super().initializePage()

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

    def as_msp(self, save_file, lipid_data):

        count = 0
        for lipid in lipid_data:
            lipid.generate_spectra()
        
            for spectrum in lipid.spectra:

                save_file.write(f"NAME: {lipid.name} {spectrum}\n")
                save_file.write(f"IONMODE: {GL.Masses[spectrum][1]}\n")
                save_file.write(f"PRECURSORMZ: {round(GL.MA(lipid, GL.Masses[spectrum]), 6)}\n")
                save_file.write(f"PRECURSORTYPE: {spectrum}\n")
                save_file.write(f"Num Peaks: {len(lipid.spectra[spectrum])}\n")
                save_file.writelines(f"{peak[0]} {peak[1]}\n" for peak in lipid.spectra[spectrum])
                save_file.write("\n")
                count += 1
                self.progress_bar.setValue(count)

        return count

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

    def as_csv(self, save_file, lipid_data): # For orbitrap inclusion list!

        count = 0
        writer = csv.writer(save_file)

        writer.writerow(['Mass [m/z]','Formula [M]','Formula type','Species','CS [z]','Polarity','Start [min]','End [min]','(N)CE','(N)CE type','MSX ID','Comment'])
        unique_mass = []
        for lipid in lipid_data:
            for adduct in lipid.adducts:
                MA = round(GL.Adduct_Spectra(lipid, adduct,{GL.MA:1}).spectrum[0][0], 6)
                if MA not in unique_mass: # This takes a lot of time!
                    unique_mass.append(MA) # Removes all the duplicates
                    writer.writerow([MA,'','',type(lipid).__name__ ,GL.Masses[adduct][2],GL.Masses[adduct][1],'','','','','',adduct])
                    count += 1
                    self.progress_bar.setValue(count)

        return count

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

    def save_as(self):

        # Generate lipids
        def Generate_Lipids():
            GL.Glycerolipid.instances = []
            tails = GL.generate_tails(self.tails_to_generate)
            for cls in self.classes_to_generate:
                for sn1, sn2 in cwr(tails, 2):
                    cls(sn2, sn1)
            return GL.Glycerolipid.instances
        lipid_list = Generate_Lipids()
        self.progress_bar.setMinimum(0)
        self.progress_bar.setMaximum(len(lipid_list))
        
        # Create save location
        file_name, _ = QFileDialog.getSaveFileName(filter="MSP (*.msp);;CSV (*.csv)")

        if file_name:

            t0 = time.time()

            # If save location exists, override
            if os.path.exists(file_name):
                os.remove(file_name) # Removes if exists
                self.output_console.appendPlainText('Overwriting {}'.format(file_name))
            else: self.output_console.appendPlainText('Creating {}'.format(file_name))

            save_file = open(file_name, 'x')

            # Based on file extension, export differently
            if '.msp' in file_name:
                count = self.as_msp(save_file, lipid_list)
            elif '.csv' in file_name:
                count = self.as_csv(save_file, lipid_list)
            else: self.output_console.appendPlainText('Unsupported file type')

            t1 = time.time()        
            self.output_console.appendPlainText(f"Generated {count} spectra in {t1-t0:.4f} seconds!")

            save_file.close()

        else: pass

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #



if __name__ == '__main__':
    app = QApplication(sys.argv)
    mainWindow = CreateWindow()
    mainWindow.show()
    status = app.exec()