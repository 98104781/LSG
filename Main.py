import os
import csv
import sys
import time

import Classes
#import Classes_isomers
import GenerateLipids as GL

from PySide6.QtCharts import QBarSeries, QBarSet, QChart, QChartView, QScatterSeries, QValueAxis
from PySide6 import QtWidgets
from PySide6.QtGui import QColor, QImage, QIntValidator, QPainter, QPainterPath, QPen, QPixmap
from itertools import combinations_with_replacement as cwr
from PySide6.QtCore import Property, QAbstractTableModel, QMargins, QModelIndex, QRect, Qt, Signal
from PySide6.QtWidgets import QComboBox, QDialog, QFileDialog, QHeaderView, QProgressBar, QStyledItemDelegate, QTableView
from PySide6.QtWidgets import QApplication, QPlainTextEdit, QPushButton, QTreeWidget, QTreeWidgetItem, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QWizard, QWizardPage

gplClassList = [cls for cls in GL.Glycerolipid.__subclasses__()]

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
        self.addPage(FASetupPage(self))
        self.addPage(SpectraSetupPage(self))
        self.addPage(GenerateSpectraPage(self))

class FASetupPage(QWizardPage):
    '''
    First page. Cmin, cmax, dmin and dmax used to
    generate fatty acids. Later, to determine lipid spectra.
    '''
    def __init__(self, parent=None):
        super(FASetupPage, self).__init__(parent)

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
        return super().isComplete()


class SpectraSetupPage(QWizardPage):
    '''
    Second page. Displays currently supported lipid classes
    with their currently supported adducts for generation.
    '''
    treeDataChanged = Signal()

    def __init__(self, parent=None):
        super(SpectraSetupPage, self).__init__(parent)

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
            self.classAdductQbox[cls] = {} # Open for adducts to add to

            for adduct in cls.adducts:
                self.classAdductQbox[cls][adduct] = QTreeWidgetItem(self.classQbox[cls])
                child = self.classAdductQbox[cls][adduct]
                # Custom variable to store fragment list
                child.fragmentList = cls.adducts[adduct]
                child.setText(0, adduct)
                child.setCheckState(0, Qt.Unchecked)
                child.setFlags(child.flags() | Qt.ItemIsUserCheckable)

        self.modifybutton = QPushButton("Modify selected adduct spectra")
        self.modifybutton.clicked.connect(self.open_editspectrawindow)

        self.vLayout.addWidget(self.treeView)
        self.vLayout.addWidget(self.modifybutton)

    def open_editspectrawindow(self):
        '''
        Opens external window to modify selected spectra.
        Should display an example lipid spectra of GPL
        16:0_18:1 with appropriate masses and intensity.
        '''
        editspectrawindow = SpectraEditWindow(self.field('tree'))
        editspectrawindow.exec()
            
    def treeData(self):
        '''
        Custom storage for tree data.
        Allows later determination of selected class/adducts.
        '''
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



class SpectraEditWindow(QDialog): # Opened from SpectraSetupPage
    '''
    Popup window to modify selected class/adduct spectra.
    '''
    def __init__(self, treeData):
        super().__init__()

        self.setWindowTitle('LSG3')
        self.setFixedSize(600, 510)
        self.vLayout = QVBoxLayout(self)
        self.vLayout2 = QVBoxLayout(self)
        self.hLayout = QHBoxLayout(self)

        # Example spectra will be GPL 16:0_18:1
        self.tails = [GL.sn(c=16, d=0, type='Acyl'),
                      GL.sn(c=18, d=1, type='Acyl')]

        self.tableView = QTableView()
        self.tableView.setFixedWidth(137)
        self.comboBox = QComboBox()
        self.spectra = SpectraHistogram()
        self.spectra.setFixedHeight(200)
        self.comboBox.currentTextChanged.connect(self.buildList)
        for item, item2 in treeData.items():
            for adduct in item2:
                self.comboBox.addItem(item.text(0)+' '+adduct.text(0), [item, adduct])

        self.applyBTN = QPushButton('Apply Changes')
        self.applyBTN.clicked.connect(self.applyChanges)
        self.resetBTN = QPushButton('Reset Changes')
        self.applyBTN.clicked.connect(self.resetChanges)

        self.vLayout.addWidget(self.comboBox)
        self.hLayout.addWidget(self.tableView)

        self.vLayout2.addWidget(self.spectra)
        self.vLayout2.addStretch()
        self.vLayout2.addWidget(self.applyBTN)
        self.hLayout.addLayout(self.vLayout2)

        self.vLayout.addLayout(self.hLayout)

    def buildLipid(self, data):
        GL.Glycerolipid.instances = []
        cls = data[0].lipidClass # Example lipid
        _, comb, *_ = cwr(self.tails, cls.No_Tails)
        example = cls(*comb) # comb will be (16:0_) 16:0_18:1
        example.adducts = {data[1].text(0):data[1].fragmentList}
        example.resolve_spectra(example.adducts)
        return example
    
    def buildList(self):

        data = self.comboBox.currentData()

        example = self.buildLipid(data)
        adduct = data[1].text(0)
        fragments = example.spectra[adduct]
        mz = GL.MA(example, GL.Masses[adduct])

        self.spectra.setSpectra(example.name+' '+adduct, mz, fragments)
        self.tableData = SpectraTableModel(fragments)
        self.tableView.setModel(self.tableData)
        self.tableView.setItemDelegateForColumn(1, SpinBoxDelegate(self.tableView))
        self.tableView.horizontalHeader().setSectionResizeMode(QHeaderView.Fixed)
        self.tableView.resizeColumnsToContents()
        self.tableView.verticalHeader().hide()
        self.tableData.dataChanged.connect(self.updateData)

    def updateData(self, index: QModelIndex, *other):
        fragment = int(index.row())
        fragmentList = self.comboBox.currentData()[1].fragmentList
        data = self.tableData.table_data[fragment]
        fragmentList[data[2]] = data[1]
        self.comboBox.currentData()[1].fragmentList = fragmentList
        self.buildList()

    def applyChanges(self):
        for x in range(self.comboBox.count()):
            data = self.comboBox.itemData(x)
            lipidClass = data[0].lipidClass
            adduct = data[1].text(0)
            fragmentList = data[1].fragmentList
            for fragment in list(fragmentList):
                if fragmentList[fragment] == 0:
                    fragmentList.pop(fragment)
            lipidClass.adducts[adduct] = fragmentList
        self.buildList()

    def resetChanges(self):
        pass

class SpectraHistogram(QChartView):
    def __init__(self, parent=None):
        super().__init__(QChart(), parent=parent)

        self.xaxis = QValueAxis()
        self.xaxis.setGridLineVisible(False)
        self.chart().addAxis(self.xaxis, Qt.AlignBottom)
        self.xaxis.setTickCount(2)
        self.xaxis.setLabelsVisible(False)

        self.yaxis = QValueAxis()
        self.yaxis.setGridLineVisible(False)
        self.chart().addAxis(self.yaxis, Qt.AlignLeft)
        self.yaxis.setLabelsVisible(False)    

        self.spectra = QScatterSeries()
        self.spectra.setMarkerSize(10000)
        self.spectra.setMarkerShape(self.spectra.MarkerShapeRectangle)
        self.chart().setMargins(QMargins(-20, 0, -10, -10))
        self.chart().addSeries(self.spectra)
        self.spectra.attachAxis(self.xaxis)
        self.spectra.attachAxis(self.yaxis)
        self.chart().legend().hide()

    def setSpectra(self, precursorName, precursorMass, spectra):

        self.chart().setTitle(precursorName)
        binMax = int(precursorMass+100)

        self.spectra.clear()

        linePath = QPainterPath()
        linePath.moveTo(5000, 5001)
        linePath.lineTo(5001, 10000)
        linePath.closeSubpath()
        image = QImage(10000, 10000, QImage.Format_ARGB32)
        painter = QPainter(image)
        pen = QPen(QColor(0, 0, 0), 1, Qt.SolidLine)
        painter.setPen(pen)
        painter.setBrush(painter.pen().color())
        painter.drawPath(linePath)
        painter.end()

        self.spectra.setBrush(image)
        self.spectra.setPen(QColor(Qt.transparent))

        for peak in spectra:
            x = peak[0]
            y = peak[1]
            self.spectra.append(x, y)

        self.yaxis.setRange(0, 100)
        self.xaxis.setRange(100, binMax)




class SpectraTableModel(QAbstractTableModel):
    '''
    Custom table type to organise and display
    lipid spectra data for inspection/modification
    using SpectraEditWindow dialogue.
    Data should be list of lists in form:
    [[mz, abundance], [mz, abundance]...]
    '''
    def __init__(self, data=[], parent=None):
        QAbstractTableModel.__init__(self, parent)
        self.table_data = data
        self.headers = ['m/z (Da)', 'Abn. (%)']
    
    def flags(self, index):
        '''
        Forbids editing of mz value displayed in column 0.
        Allows editing of abundance value displayed in column 1.
        '''
        if index.column() > 0:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable
        else:
            return Qt.ItemIsSelectable

    def rowCount(self, parent):
        return len(self.table_data)

    def columnCount(self, parent):
        return len(self.headers)
    
    def headerData(self, column, orientation, role=Qt.DisplayRole):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return self.headers[column]
    
    def data(self, index, role=Qt.DisplayRole):
        row = index.row()
        column = index.column()
        if role == Qt.DisplayRole:
            return str(self.table_data[row][column])

    def setData(self, index, value, role=Qt.EditRole):
        if role == Qt.EditRole:
            row = index.row()
            column = index.column()
            self.table_data[row][column] = value
            self.dataChanged.emit(index, index)
            return True
        return QAbstractTableModel.setData(self, index, value, role)

    def resizeEvent(self, event):
        pass

    def insertRow(self):
        pass

    def removeRow(self, position):
        pass

class SpinBoxDelegate(QStyledItemDelegate):
    def __init__(self, parent=None):
        QStyledItemDelegate.__init__(self, parent)

    def createEditor(self, parent, option, index):
        editor = QtWidgets.QSpinBox(parent)
        editor.setMinimum(0)
        editor.setMaximum(100)
        editor.setSingleStep(1)
        return editor
            
    def setEditorData(self, spinBox, index):
        value = int(index.model().data(index, Qt.DisplayRole))
        spinBox.setValue(value)

    def setModelData(self, spinBox, model, index):
        value = spinBox.value()
        model.setData(index, value)

    def updateEditorGeometry(self, editor, option, index):
        editor.setGeometry(option.rect)

class GenerateSpectraPage(QWizardPage):
    '''
    Final page of GUI. Allows for generation.
    '''
    def __init__(self, parent=None):
        super(GenerateSpectraPage, self).__init__(parent)

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
        '''
        Defines how to export data when saved as .MSP.
        Contains lipid fragmentation informaiton.
        '''
        count = 0
        for lipid in lipid_data:
            lipid.resolve_spectra(lipid.adducts)
        
            for adduct in lipid.spectra:

                save_file.write(f"NAME: {lipid.name} {adduct}\n")
                save_file.write(f"IONMODE: {GL.Masses[adduct][1]}\n")
                save_file.write(f"PRECURSORMZ: {round(GL.MA(lipid, GL.Masses[adduct]), 6)}\n")
                save_file.write(f"COMPOUNDCLASS: {lipid.lipid_class}\n")
                save_file.write(f"FORMULA: {''.join(''.join((key, str(val))) for (key, val) in lipid.formula.items())} \n")
                save_file.write(f"RETENTIONTIME: 0.00\n") # Pointless
                save_file.write(f"PRECURSORTYPE: {adduct}\n")
                save_file.write(f"Num Peaks: {len(lipid.spectra[adduct])}\n")
                save_file.writelines(f"{frag[0]} {frag[1]}\n" for frag in lipid.spectra[adduct] if frag[1] > 0)
                save_file.write("\n")
                count += 1

        self.progress_bar.setValue(self.progress_bar.value()+1)
        return count

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

    def as_csv(self, save_file, lipid_data):
        '''
        Defines how to export data when saved as .CSV.
        Specifically for use in QE 'Orbitrap' inclusion list.
        '''
        count = 0
        writer = csv.writer(save_file)

        writer.writerow(['Mass [m/z]','Formula [M]','Formula type','Species','CS [z]','Polarity','Start [min]','End [min]','(N)CE','(N)CE type','MSX ID','Comment'])
        unique_mass = []

        for lipid in lipid_data:
            for adduct in lipid.adducts:
                ma = round(GL.MA(lipid, GL.Masses[adduct]), 6) # lipid.resolve_frag(adduct, {GL.MA:1})[0][0]
                if ma not in unique_mass: # This takes a lot of time!
                    unique_mass.append(ma) # Removes all the duplicates
                    writer.writerow([ma,'','',type(lipid).__name__ ,GL.Masses[adduct][2],GL.Masses[adduct][1],'','','','','',adduct])
                    count +=1
                else: continue

        self.progress_bar.setValue(self.progress_bar.value()+1)
        return count

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

    def save_as(self):
        '''
        Popup 'Save as' dialogue box.
        '''
        def Generate_Lipids():

            GL.Glycerolipid.instances = []

            self.progress_bar.setMinimum(0)
            self.progress_bar.setMaximum(len(self.classes_to_generate)+2)

            tails = GL.generate_tails(self.tails_to_generate)
            self.progress_bar.setValue(self.progress_bar.value()+1)

            for cls in self.classes_to_generate:
                for comb in cwr(tails, cls.No_Tails):
                    cls(*comb)
                self.progress_bar.setValue(self.progress_bar.value()+1)

            return GL.Glycerolipid.instances


        # Create save location
        file_name, _ = QFileDialog.getSaveFileName(filter="MSP (*.msp);;CSV (*.csv)")
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
            if '.msp' in file_name:
                count = self.as_msp(save_file, Generate_Lipids())
            elif '.csv' in file_name:
                count = self.as_csv(save_file, Generate_Lipids())
            else: self.output_console.appendPlainText('Unsupported file type')
            save_file.close()

            t1 = time.time()
            self.progress_bar.setValue(self.progress_bar.value()+1)        
            self.output_console.appendPlainText(f"Generated {count} spectra in {t1-t0:.4f} seconds!")

        else: pass

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #



if __name__ == '__main__':
    app = QApplication(sys.argv)
    mainWindow = CreateWindow()
    mainWindow.show()
    status = app.exec()