import os
import time

import GenerateLipids as GL
import SaveAs

from PySide6.QtWidgets import QFileDialog

from PySide6.QtCore import QThread
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
        self.setSubTitle("Press 'Generate' to create and export file \n"
                         ".msp spectral libraries, .csv QE+ inclusion list and .csv skyline transition lists supported")
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

    def unsupported_fileType(self):
        self.generatorThread.exit()
        self.progress_bar.setMaximum(1)
        self.generatebutton.setEnabled(True)
        self.output_console.appendPlainText('Unsupported file type')

    def completionText(self):
        self.generatorThread.exit()
        self.t1 = time.time()
        self.progress_bar.setMaximum(1)
        self.progress_bar.setValue(1)
        self.generatebutton.setEnabled(True)
        self.output_console.appendPlainText(f"Generated {self.generatorObject.count} {self.generatorObject.noun} in {self.t1-self.t0:.4f} seconds!")

    def classCompleted(self, cls):
        consoleText = self.output_console.toPlainText()
        consoleText = consoleText.replace('- '+cls.__name__+' ', '- '+cls.__name__+' - Completed')
        self.output_console.setPlainText(consoleText)

    def save_as(self):
        '''Popup 'Save as' dialogue box'''
        # Create save location
        file_name, filter = QFileDialog.getSaveFileName(filter="MSP (*.msp);;Orbitrap Inclusion (*.csv);;Skyline Transition (*.csv)", selectedFilter='')

        if file_name:
            if os.path.exists(file_name): # If save location exists, override
                self.output_console.appendPlainText('Overwriting {}...'.format(file_name))
                try:os.remove(file_name) # Removes if exists
                except PermissionError:
                    self.output_console.appendPlainText("Could not overwrite file. File may be in use.")
                    pass
            else: self.output_console.appendPlainText('Creating {}'.format(file_name))

            try:
                self.progress_bar.setValue(0) # Sets loading animation
                self.progress_bar.setMaximum(0) # for progress bar
                consoleText = self.output_console.toPlainText()
                consoleText = consoleText.replace(' - Completed', ' ')
                self.output_console.setPlainText(consoleText)

                self.generatorThread = QThread() # Generator on second thread so GUI doesn't lag!
                self.generatorObject = SaveAs.Generator(file_name, filter,
                    self.classes_to_generate, self.tails_to_generate, self.bases_to_generate, 
                    self.field('isomerism'), self.field('lipidSpecific'), self.field('lipidList'),
                    self.field('tailSpecific'), self.field('tailList'), self.field('specificOrganisation'))
                self.generatorObject.moveToThread(self.generatorThread)
                self.generatorObject.fileError.connect(self.unsupported_fileType)
                self.generatorThread.started.connect(self.generatorObject.run)
                self.generatorObject.finished.connect(self.completionText)
                self.generatorObject.finished.connect(self.generatorObject.deleteLater)
                self.generatorObject.finished.connect(self.completeChanged)
                self.generatorThread.finished.connect(self.generatorThread.deleteLater)
                self.generatorObject.progress.connect(self.classCompleted)
                self.generatorThread.start()
                self.generatebutton.setEnabled(False)
                self.completeChanged.emit()
                self.t0 = time.time()
            except: self.output_console.appendPlainText('Save location unavailable')
        else: pass

    def initializePage(self) -> None:
        '''
        When the page is opened, update console information.
        '''
        self.output_console.clear() # Clear console, update with tails and lipids chosen

        try:
            if self.generatorThread.isRunning():
                self.progress_bar.setValue(0)
                self.progress_bar.setMaximum(0)
                self.generatebutton.setEnabled(False)
        except: pass

        self.classes_to_generate = [] # Empty lists initialised here
        self.tails_to_generate = [] # otherwise error thrown when
        self.bases_to_generate = [] # coming from page 2B

        if self.field('lipidSpecific'): # If generate specific lipids is selected:
            self.output_console.appendPlainText('The following specific lipids will be generated:')
            for lipid in self.field('lipidList'):
                self.output_console.appendPlainText(lipid[0].name+' '+lipid[1])
       
        else: # If generate lipid range is selected:
            # Print to console the range of tails to be generated
            if self.field('tailSpecific'):
                self.output_console.appendPlainText('The following tails will be used:')
                self.output_console.appendPlainText(', '.join(tail.name for tail in self.field('tailList'))+'\n')
            else:
                self.output_console.appendPlainText('Tails will be generated from '
                                        +self.field('cmin')+':'+self.field('dmin')+
                                ' -> '+self.field('cmax')+':'+self.field('dmax'))
                if self.field('hydroxytickbox'): # Include any hydroxy tails in console too!
                    self.output_console.appendPlainText('Hydroxy-functionalised tails included.\n')

            self.selected_class_adducts = self.field('tree')
            self.output_console.appendPlainText('The following classes will be generated:')
            caString = '' # Generates string, and class list
            for item, item2 in self.selected_class_adducts.items():
                caString += '- '+item.text(0)+' \n'
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
    
    def isComplete(self):
        try: # Disable finish button if running.
            if self.generatorThread.isRunning():
                return False
        except: pass # Thread not yet created.
        return super().isComplete()