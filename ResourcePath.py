import os
import sys

def resource_path(relative_path):
    base_path = getattr(sys, '_MEIPASS', os.path.dirname(os.path.abspath(__file__)))
    return os.path.join(base_path, relative_path)

def exe_path(relative_path):
    if getattr(sys, 'frozen', False):  # check if the application is bundled
        base_path = os.path.dirname(sys.executable)  # as Pyinstaller .exe
    else:
        base_path = getattr(sys, '_MEIPASS', os.path.dirname(os.path.abspath(__file__)))
    return os.path.join(base_path, relative_path)