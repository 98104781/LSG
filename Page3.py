import os
import time

import GenerateLipids as GL
import SaveAs

from itertools import combinations_with_replacement as cwr, product
from PySide6.QtWidgets import QFileDialog

from PySide6.QtGui import QPixmap
from PySide6.QtWidgets import QProgressBar
from PySide6.QtWidgets import QPlainTextEdit, QPushButton, QVBoxLayout, QWizard, QWizardPage

class Page(QWizardPage):
    '''
    Final page of GUI. Allows for generation.
    '''
    def __init__(self, parent=None):
        super(Page, self).__init__(parent)

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
        '''
        When the page is opened, update console information.
        '''
        self.selected_class_adducts = self.field('tree')

        self.output_console.clear()
        self.output_console.appendPlainText('Tails will be generated from '
                                +self.field('cmin')+':'+self.field('dmin')+
                          ' -> '+self.field('cmax')+':'+self.field('dmax'))
        if self.field('hydroxytickbox'):
            self.output_console.appendPlainText('Hydroxy-functionalised tails included.\n')
        self.output_console.appendPlainText('The following classes will be generated:')

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

    def save_as(self):
        '''
        Popup 'Save as' dialogue box.
        '''
        def Generate_Lipids():

            GL.Glycerolipid.instances = []

            self.progress_bar.setMinimum(0)
            self.progress_bar.setMaximum(len(self.classes_to_generate)+2)

            tails = GL.generate_acyl_tails(self.tails_to_generate)
            if self.field('hydroxytickbox'):
                hydroxytails = GL.generate_hydroxy_tails(self.tails_to_generate, int(self.field('omax')))
                tails.extend(hydroxytails)

            self.progress_bar.setValue(self.progress_bar.value()+1)

            if self.field('isomerism') == True:

                for cls in self.classes_to_generate:
                    for comb in product(tails, repeat = cls.No_Tails):
                        cls(*comb)
                    self.progress_bar.setValue(self.progress_bar.value()+1)
                return GL.Glycerolipid.instances
            else:
                for cls in self.classes_to_generate:
                    for comb in cwr(tails, cls.No_Tails):
                        cls(*comb)
                    self.progress_bar.setValue(self.progress_bar.value()+1)
                return GL.Glycerolipid.instances

        # Create save location
        file_name, filter = QFileDialog.getSaveFileName(filter="MSP (*.msp);;Orbitrap Inclusion (*.csv);;Skyline Transition (*.csv)", selectedFilter='')
        if file_name:

            t0 = time.time()

            # If save location exists, override
            if os.path.exists(file_name):
                self.output_console.appendPlainText('Overwriting {}...'.format(file_name))
                try:os.remove(file_name) # Removes if exists
                except PermissionError:
                    self.output_console.appendPlainText("Could not overwrite file. File may be in use.")
                    pass
            else: self.output_console.appendPlainText('Creating {}'.format(file_name))

            save_file = open(file_name, 'x', newline='')
            # Based on file extension, export differently
            if filter == "MSP (*.msp)":
                count = SaveAs.as_msp(self, save_file, Generate_Lipids())
            elif filter == "Orbitrap Inclusion (*.csv)":
                count = SaveAs.as_orb(self, save_file, Generate_Lipids())
            elif filter =="Skyline Transition (*.csv)":
                count = SaveAs.as_sky(self, save_file, Generate_Lipids())
                pass
            else: self.output_console.appendPlainText('Unsupported file type')
            save_file.close()

            t1 = time.time()
            self.progress_bar.setValue(self.progress_bar.value()+1)        
            self.output_console.appendPlainText(f"Generated {count} spectra in {t1-t0:.4f} seconds!")

        else: pass

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #