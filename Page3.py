from asyncio import base_tasks
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

        self.output_console.clear() # Clear console, update with tails and lipids chosen
        if self.field('radButton2'): # If generate specific lipids is selected:
            self.output_console.appendPlainText('The following specific lipids will be generated:')
            for lipid in self.field('lipidList'):
                self.output_console.appendPlainText(lipid[0].name+' '+lipid[1])
        else: # If generate lipid range is selected:
            self.selected_class_adducts = self.field('tree')
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

            # Prepare information to generate lipids:

            # List with limits for generated tails
            self.tails_to_generate = [int(self.field('cmin') or 0), int(self.field('cmax') or 0),
                                      int(self.field('dmin') or 0), int(self.field('dmax') or 0),
                                      int(self.field('omax') or 0)]
            

            base_types = [] # List with types of sphingoid bases to generate
            for cls in self.classes_to_generate: # Sphingolipids can have one of many base type, which the lipid
                if issubclass(cls, GL.Sphingolipid): # is centered around. The possible types are defined in the
                    base_types.extend(cls.base_types) # Sphingolipid.base_types list. Collect all unique base types 
                base_types = list(set(base_types)) # used so that they can be generated with the lipids.

            if self.field('ceramideVariability') is False:
                self.bases_to_generate = [18, 18, base_types]
            else: self.bases_to_generate = [int(self.field('cmin') or 0), int(self.field('cmax') or 0), base_types]
            
        self.progress_bar.reset()


        return super().initializePage()
    
    # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

    def save_as(self):
        '''
        Popup 'Save as' dialogue box.
        '''
        def generate_range():
            
            self.progress_bar.setMinimum(0)
            self.progress_bar.setMaximum(len(self.classes_to_generate)+2)
            self.tails = GL.generate_acyl_tails(self.tails_to_generate)
            self.bases = GL.generate_base_tails(self.bases_to_generate)
            self.progress_bar.setValue(self.progress_bar.value()+1)

            if self.field('isomerism') == True:
                for cls in self.classes_to_generate:
                    for adduct in cls.adducts: cls.adducts[adduct] = {k: v for k, v in cls.adducts[adduct].items() if v != 0}
                    if issubclass(cls, GL.Glycerolipid):
                        for comb in product(self.tails, repeat = cls.No_Tails):
                            yield cls(*comb)
                    elif issubclass(cls, GL.Sphingolipid):
                        for type in cls.base_types:
                            for comb in product(self.bases[type], self.tails):
                                yield cls(*comb)
                    else: pass
                    self.progress_bar.setValue(self.progress_bar.value()+1)

            else:
                for cls in self.classes_to_generate:
                    for adduct in cls.adducts: cls.adducts[adduct] = {k: v for k, v in cls.adducts[adduct].items() if v != 0}
                    if issubclass(cls, GL.Glycerolipid):
                        for comb in cwr(self.tails, cls.No_Tails):
                            yield cls(*comb)
                    elif issubclass(cls, GL.Sphingolipid):
                        for type in cls.base_types:
                            for comb in product(self.bases[type], self.tails):
                                yield cls(*comb)
                    else: pass
                    self.progress_bar.setValue(self.progress_bar.value()+1)

        def generate_specific():

            lipidList = self.field('lipidList')
            self.progress_bar.setMinimum(0)
            self.progress_bar.setMaximum(len(lipidList)+1)
            
            for lipid, selected_adduct in lipidList:
                adducts = lipid.adducts
                lipid.adducts = {k: v for k, v in adducts.items() if k is selected_adduct}
                lipid.adducts[selected_adduct] = {k: v for k, v in adducts[selected_adduct].items() if v != 0}
                self.progress_bar.setValue(self.progress_bar.value()+1)
            return [x[0] for x in lipidList]


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

            if self.field('radButton2'):
                generate = generate_specific
            else: generate = generate_range

            try: 
                save_file = open(file_name, 'x', newline='')
                # Based on file extension, export differently
                if filter == "MSP (*.msp)":
                    count = SaveAs.as_msp(self, save_file, generate())
                elif filter == "Orbitrap Inclusion (*.csv)":
                    count = SaveAs.as_orb(self, save_file, generate())
                elif filter =="Skyline Transition (*.csv)":
                    count = SaveAs.as_sky(self, save_file, generate())
                    pass
                else: self.output_console.appendPlainText('Unsupported file type')
                save_file.close()
            except: self.output_console.appendPlainText('File unavailable')

            t1 = time.time()
            self.progress_bar.setValue(self.progress_bar.value()+1)        
            self.output_console.appendPlainText(f"Generated {count} spectra in {t1-t0:.4f} seconds!")

        else: pass

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #