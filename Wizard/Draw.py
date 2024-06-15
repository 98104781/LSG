from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
import base64

from PySide6.QtCore import QByteArray, Qt
from PySide6.QtGui import QPixmap, QTransform, QPainter, QColor, QFont, QPen 

def drawMolecule(smiles = 'OCC(O)CO', width = 0, height = 0):

    drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
    drawer.drawOptions().addStereoAnnotation = True
    drawer.drawOptions().unspecifiedStereoIsUnknown = True
    molStructure = Chem.MolFromSmiles(smiles)
    drawer.DrawMolecule(molStructure)
    drawer.FinishDrawing()

    molText = base64.b64encode(drawer.GetDrawingText()).decode('utf8')
    pixmap = QPixmap()
    pixmap.loadFromData(QByteArray.fromBase64(molText.encode()))
    
    return pixmap

def rotatePixmap(pixmap, rotation):

    pix = pixmap.transformed(QTransform().rotate(rotation))

    return pix

def framePixmap(pixmap, borderColour=QColor(0, 0, 0)):

    borderThickness = 1
    frameWidth = pixmap.width() + 2 * borderThickness
    frameHeight = pixmap.height() + 2 * borderThickness
    
    borderPixmap = QPixmap(frameWidth, frameHeight)
    borderPixmap.fill(borderColour)
    painter = QPainter(borderPixmap)
    
    painter.drawPixmap(borderThickness, borderThickness, pixmap)
    painter.end()
    
    return borderPixmap

def drawText(text = '', width=0, height=0, textColour=QColor(0, 0, 0)):

    font = QFont()
    font.setPointSize(24)
    font.setFamily("Times New Roman")
    
    pixmap = QPixmap(width, height)
    pixmap.fill(QColor(255, 255, 255))
    painter = QPainter(pixmap)
    painter.setRenderHint(QPainter.Antialiasing)
    
    painter.setFont(font)
    painter.setPen(textColour)
    painter.drawText(pixmap.rect(), Qt.AlignCenter, text)
    
    painter.end()

    return pixmap