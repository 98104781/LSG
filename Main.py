import sys
import Export
from PySide6.QtWidgets import QApplication, QMainWindow

class CreateWindow(QMainWindow):

    def __init__(self):
        super().__init__()

        self.setWindowTitle('LSG3')
        self.setGeometry(300, 100, 600, 300)

        self.menu_Bar()

    def menu_Bar(self):
        menu_bar = self.menuBar()
        file = menu_bar.addMenu('&File')
        
        file.addAction('Export', Export.save_as)

app = QApplication(sys.argv)
mainWindow = CreateWindow()
mainWindow.show()
status = app.exec()