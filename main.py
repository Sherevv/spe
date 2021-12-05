#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from program_gui import Ui_MainWindow
from planet_satellite_system import *
from PyQt4 import QtGui

if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    MainWindow = QtGui.QMainWindow()
    w3d = Widget3D()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    w3d.connect(MainWindow, ui)
    MainWindow.show()
    MainWindow.iren.Initialize()
    MainWindow.iren.Start()
    sys.exit(app.exec_())
