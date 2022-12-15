import Lipids.GenerateLipids as GL

from PySide6.QtCore import Signal
from PySide6.QtGui import QIntValidator
from PySide6.QtWidgets import QDialog, QPushButton, QVBoxLayout, QHBoxLayout, QLineEdit, QComboBox

class TailWindow(QDialog):
    '''
    Popup window to specify tail stats.
    '''
    output = Signal(GL.sn)

    def __init__(self, type='A'):
        super().__init__()

        self.setWindowTitle('LSG3')
        self.setFixedSize(600, 125)
        self.vLayout = QVBoxLayout(self)
        self.hLayout = QHBoxLayout(self)

        self.ty = QComboBox()
        self.hLayout.addWidget(self.ty)
        self.c = QLineEdit()
        self.c.setPlaceholderText('Chain Length')
        self.hLayout.addWidget(self.c)        
        self.d = QLineEdit()
        self.d.setPlaceholderText('Desaturation')
        self.hLayout.addWidget(self.d)
        self.oh = QLineEdit()  
        self.oh.setPlaceholderText('# Oxidation')
        self.hLayout.addWidget(self.oh)
        self.dt = QLineEdit()
        self.hLayout.addWidget(self.dt)
        self.dt.setPlaceholderText('# of deuterium labels')

        self.ty.currentIndexChanged.connect(self.updateComboBox)

        if type == 'A': self.ty.addItem('Acyl')
        elif type == 'O': self.ty.addItem('Ether')
        elif type == 'P': self.ty.addItem('Vinyl')
        elif type == 'X':
            self.ty.addItem('Acyl')
            self.ty.addItem('Ether')
            self.ty.addItem('Vinyl')
            for item in GL.baseTypes:
                self.ty.addItem(item)

        self.tail = None
        self.acceptButton = QPushButton('Accept')
        self.acceptButton.clicked.connect(self.acceptTail)
        self.vLayout.addLayout(self.hLayout)
        self.vLayout.addWidget(self.acceptButton)

    def updateComboBox(self):

        if self.ty.currentText() in ['Acyl', 'Ether', 'Vinyl']:
            self.c.setValidator(QIntValidator(1, 100))
            self.d.setEnabled(True)
            self.d.setValidator(QIntValidator(1, 100))
            self.oh.setEnabled(True)
            self.oh.setValidator(QIntValidator(1, 100))
            self.dt.setValidator(QIntValidator(1, 100))

        elif self.ty.currentText() in GL.baseTypes:
            self.c.setValidator(QIntValidator(6, 100))
            self.d.setEnabled(False)
            self.d.setText('')
            self.oh.setEnabled(False)
            self.oh.setText('')
            self.dt.setValidator(QIntValidator(1, 100))

    def acceptTail(self):

        ty = self.ty.currentText()
        c = self.c.text()
        d = self.d.text()
        oh= self.oh.text()
        dt= self.dt.text()

        if ty in ['Acyl', 'Ether', 'Vinyl']:
            try: # Tail must be possible
                if (int(d or 0)+int(oh or 0)) > int(c or 1)/2: return # D + -OH < C
                if int(dt or 0) > (2*int(c or 1)-2*int(d or 0)-1): return # Deuterium <= Hydrogens
                if int(c or 0) < 2: return
                self.tail = GL.sn(int(c), int(d or 0), type=ty, oh=int(oh or 0), dt=int(dt or 0))
                self.output.emit(self.tail)
                self.done(1)
                self.close()
            except: pass

        elif ty in GL.baseTypes:
            try: # Base must be possible
                if int(dt or 0) > (2*int(c or 1)): return # Deuterium <= Hydrogens
                if int(c or 0) < 7: return
                self.tail = GL.base(int(c), type=ty, dt=int(dt or 0))
                self.output.emit(self.tail)
                self.done(1)
                self.close()
            except: pass

            
class BaseWindow(QDialog):
    '''
    Popup window to specify base stats.
    '''

    output = Signal(GL.base)

    def __init__(self, types):
        super().__init__()

        self.setWindowTitle('LSG3')
        self.setFixedSize(600, 125)
        self.vLayout = QVBoxLayout(self)
        self.hLayout = QHBoxLayout(self)
        self.types = types
    
        self.ty = QComboBox()
        for type in self.types:
            self.ty.addItem(type)
        self.hLayout.addWidget(self.ty)

        self.c = QLineEdit()
        self.c.setPlaceholderText('Chain Length')
        self.c.setValidator(QIntValidator(6, 100))
        self.hLayout.addWidget(self.c)

        self.dt = QLineEdit()
        self.dt.setPlaceholderText('# of deuterium labels')
        self.dt.setValidator(QIntValidator(1, 100))
        self.hLayout.addWidget(self.dt)
        self.vLayout.addLayout(self.hLayout)

        self.tail = None
        self.acceptButton = QPushButton('Accept')
        self.acceptButton.clicked.connect(self.acceptTail)
        self.vLayout.addWidget(self.acceptButton)

    def acceptTail(self):

        ty= self.ty.currentText()
        c = self.c.text()
        dt= self.dt.text()

        try: # Base must be possible
            if int(dt or 0) > (2*int(c or 1)): return # Deuterium <= Hydrogens
            self.tail = GL.base(int(c), type=ty, dt=int(dt or 0))
            self.output.emit(self.tail)
            self.done(1)
            self.close()
        except: pass