# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'program_gui.ui'
#
# Created: Fri Mar  3 21:46:23 2017
#      by: PyQt4 UI code generator 4.10.4
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName(_fromUtf8("MainWindow"))
        MainWindow.resize(1191, 835)
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.gridLayout_2 = QtGui.QGridLayout(self.centralwidget)
        self.gridLayout_2.setContentsMargins(0, 0, 0, -1)
        self.gridLayout_2.setObjectName(_fromUtf8("gridLayout_2"))
        self.scrollArea = QtGui.QScrollArea(self.centralwidget)
        self.scrollArea.setWidgetResizable(True)
        self.scrollArea.setObjectName(_fromUtf8("scrollArea"))
        self.scrollAreaWidgetContents = QtGui.QWidget()
        self.scrollAreaWidgetContents.setGeometry(QtCore.QRect(0, 0, 1189, 775))
        self.scrollAreaWidgetContents.setMinimumSize(QtCore.QSize(1180, 760))
        self.scrollAreaWidgetContents.setObjectName(_fromUtf8("scrollAreaWidgetContents"))
        self.controlFrame = QtGui.QFrame(self.scrollAreaWidgetContents)
        self.controlFrame.setGeometry(QtCore.QRect(10, 500, 1171, 271))
        self.controlFrame.setFrameShape(QtGui.QFrame.StyledPanel)
        self.controlFrame.setFrameShadow(QtGui.QFrame.Raised)
        self.controlFrame.setObjectName(_fromUtf8("controlFrame"))
        self.infoBox = QtGui.QGroupBox(self.controlFrame)
        self.infoBox.setGeometry(QtCore.QRect(990, 20, 181, 241))
        self.infoBox.setObjectName(_fromUtf8("infoBox"))
        self.gridLayoutWidget_2 = QtGui.QWidget(self.infoBox)
        self.gridLayoutWidget_2.setGeometry(QtCore.QRect(10, 30, 160, 98))
        self.gridLayoutWidget_2.setObjectName(_fromUtf8("gridLayoutWidget_2"))
        self.infoGridLayout_1 = QtGui.QGridLayout(self.gridLayoutWidget_2)
        self.infoGridLayout_1.setMargin(0)
        self.infoGridLayout_1.setObjectName(_fromUtf8("infoGridLayout_1"))
        self.evol_i = QtGui.QLineEdit(self.gridLayoutWidget_2)
        self.evol_i.setObjectName(_fromUtf8("evol_i"))
        self.infoGridLayout_1.addWidget(self.evol_i, 2, 1, 1, 1)
        self.label_info_e = QtGui.QLabel(self.gridLayoutWidget_2)
        self.label_info_e.setObjectName(_fromUtf8("label_info_e"))
        self.infoGridLayout_1.addWidget(self.label_info_e, 0, 0, 1, 1)
        self.label_info_i = QtGui.QLabel(self.gridLayoutWidget_2)
        self.label_info_i.setObjectName(_fromUtf8("label_info_i"))
        self.infoGridLayout_1.addWidget(self.label_info_i, 2, 0, 1, 1)
        self.label_info_a = QtGui.QLabel(self.gridLayoutWidget_2)
        self.label_info_a.setObjectName(_fromUtf8("label_info_a"))
        self.infoGridLayout_1.addWidget(self.label_info_a, 1, 0, 1, 1)
        self.evol_a = QtGui.QLineEdit(self.gridLayoutWidget_2)
        self.evol_a.setObjectName(_fromUtf8("evol_a"))
        self.infoGridLayout_1.addWidget(self.evol_a, 1, 1, 1, 1)
        self.evol_e = QtGui.QLineEdit(self.gridLayoutWidget_2)
        self.evol_e.setReadOnly(True)
        self.evol_e.setObjectName(_fromUtf8("evol_e"))
        self.infoGridLayout_1.addWidget(self.evol_e, 0, 1, 1, 1)
        self.manageBox = QtGui.QGroupBox(self.controlFrame)
        self.manageBox.setGeometry(QtCore.QRect(10, 20, 161, 241))
        self.manageBox.setObjectName(_fromUtf8("manageBox"))
        self.verticalLayoutWidget = QtGui.QWidget(self.manageBox)
        self.verticalLayoutWidget.setGeometry(QtCore.QRect(10, 30, 141, 95))
        self.verticalLayoutWidget.setObjectName(_fromUtf8("verticalLayoutWidget"))
        self.buttonVerticalLayout = QtGui.QVBoxLayout(self.verticalLayoutWidget)
        self.buttonVerticalLayout.setMargin(0)
        self.buttonVerticalLayout.setObjectName(_fromUtf8("buttonVerticalLayout"))
        self.startButton = QtGui.QPushButton(self.verticalLayoutWidget)
        self.startButton.setObjectName(_fromUtf8("startButton"))
        self.buttonVerticalLayout.addWidget(self.startButton)
        self.stopButton = QtGui.QPushButton(self.verticalLayoutWidget)
        self.stopButton.setObjectName(_fromUtf8("stopButton"))
        self.buttonVerticalLayout.addWidget(self.stopButton)
        self.solveButton = QtGui.QPushButton(self.verticalLayoutWidget)
        self.solveButton.setObjectName(_fromUtf8("solveButton"))
        self.buttonVerticalLayout.addWidget(self.solveButton)
        self.paramsBox = QtGui.QGroupBox(self.controlFrame)
        self.paramsBox.setGeometry(QtCore.QRect(180, 20, 791, 241))
        self.paramsBox.setObjectName(_fromUtf8("paramsBox"))
        self.gridLayoutWidget = QtGui.QWidget(self.paramsBox)
        self.gridLayoutWidget.setGeometry(QtCore.QRect(420, 30, 371, 95))
        self.gridLayoutWidget.setObjectName(_fromUtf8("gridLayoutWidget"))
        self.paramsGridLayout_2 = QtGui.QGridLayout(self.gridLayoutWidget)
        self.paramsGridLayout_2.setMargin(0)
        self.paramsGridLayout_2.setObjectName(_fromUtf8("paramsGridLayout_2"))
        self.label_a = QtGui.QLabel(self.gridLayoutWidget)
        self.label_a.setObjectName(_fromUtf8("label_a"))
        self.paramsGridLayout_2.addWidget(self.label_a, 1, 0, 1, 1)
        self.label_ecc = QtGui.QLabel(self.gridLayoutWidget)
        self.label_ecc.setObjectName(_fromUtf8("label_ecc"))
        self.paramsGridLayout_2.addWidget(self.label_ecc, 0, 0, 1, 1)
        self.label_i = QtGui.QLabel(self.gridLayoutWidget)
        self.label_i.setObjectName(_fromUtf8("label_i"))
        self.paramsGridLayout_2.addWidget(self.label_i, 2, 0, 1, 1)
        self.label_i_2 = QtGui.QLabel(self.gridLayoutWidget)
        self.label_i_2.setObjectName(_fromUtf8("label_i_2"))
        self.paramsGridLayout_2.addWidget(self.label_i_2, 2, 2, 1, 1)
        self.label_a_raz = QtGui.QLabel(self.gridLayoutWidget)
        self.label_a_raz.setTextFormat(QtCore.Qt.AutoText)
        self.label_a_raz.setObjectName(_fromUtf8("label_a_raz"))
        self.paramsGridLayout_2.addWidget(self.label_a_raz, 1, 2, 1, 1)
        self.edit_e = QtGui.QDoubleSpinBox(self.gridLayoutWidget)
        self.edit_e.setDecimals(6)
        self.edit_e.setMaximum(1.0)
        self.edit_e.setObjectName(_fromUtf8("edit_e"))
        self.paramsGridLayout_2.addWidget(self.edit_e, 0, 1, 1, 1)
        self.edit_a = QtGui.QDoubleSpinBox(self.gridLayoutWidget)
        self.edit_a.setDecimals(4)
        self.edit_a.setMaximum(999999999.99)
        self.edit_a.setObjectName(_fromUtf8("edit_a"))
        self.paramsGridLayout_2.addWidget(self.edit_a, 1, 1, 1, 1)
        self.edit_i = QtGui.QDoubleSpinBox(self.gridLayoutWidget)
        self.edit_i.setDecimals(4)
        self.edit_i.setMinimum(-180.0)
        self.edit_i.setMaximum(180.0)
        self.edit_i.setObjectName(_fromUtf8("edit_i"))
        self.paramsGridLayout_2.addWidget(self.edit_i, 2, 1, 1, 1)
        self.gridLayoutWidget_3 = QtGui.QWidget(self.paramsBox)
        self.gridLayoutWidget_3.setGeometry(QtCore.QRect(10, 30, 401, 146))
        self.gridLayoutWidget_3.setObjectName(_fromUtf8("gridLayoutWidget_3"))
        self.paramsGridLayout_1 = QtGui.QGridLayout(self.gridLayoutWidget_3)
        self.paramsGridLayout_1.setMargin(0)
        self.paramsGridLayout_1.setObjectName(_fromUtf8("paramsGridLayout_1"))
        self.label_Tp_raz = QtGui.QLabel(self.gridLayoutWidget_3)
        self.label_Tp_raz.setObjectName(_fromUtf8("label_Tp_raz"))
        self.paramsGridLayout_1.addWidget(self.label_Tp_raz, 2, 2, 1, 1)
        self.label_mp = QtGui.QLabel(self.gridLayoutWidget_3)
        self.label_mp.setObjectName(_fromUtf8("label_mp"))
        self.paramsGridLayout_1.addWidget(self.label_mp, 0, 0, 1, 1)
        self.label_Tp = QtGui.QLabel(self.gridLayoutWidget_3)
        self.label_Tp.setObjectName(_fromUtf8("label_Tp"))
        self.paramsGridLayout_1.addWidget(self.label_Tp, 2, 0, 1, 1)
        self.edit_ms = QtGui.QDoubleSpinBox(self.gridLayoutWidget_3)
        self.edit_ms.setDecimals(4)
        self.edit_ms.setMaximum(999999999.99)
        self.edit_ms.setObjectName(_fromUtf8("edit_ms"))
        self.paramsGridLayout_1.addWidget(self.edit_ms, 3, 1, 1, 1)
        self.label_mp_raz = QtGui.QLabel(self.gridLayoutWidget_3)
        self.label_mp_raz.setObjectName(_fromUtf8("label_mp_raz"))
        self.paramsGridLayout_1.addWidget(self.label_mp_raz, 0, 2, 1, 1)
        self.label_r0 = QtGui.QLabel(self.gridLayoutWidget_3)
        self.label_r0.setObjectName(_fromUtf8("label_r0"))
        self.paramsGridLayout_1.addWidget(self.label_r0, 1, 0, 1, 1)
        self.edit_mp = QtGui.QDoubleSpinBox(self.gridLayoutWidget_3)
        self.edit_mp.setDecimals(4)
        self.edit_mp.setMaximum(999999999.99)
        self.edit_mp.setObjectName(_fromUtf8("edit_mp"))
        self.paramsGridLayout_1.addWidget(self.edit_mp, 0, 1, 1, 1)
        self.edit_Tp = QtGui.QDoubleSpinBox(self.gridLayoutWidget_3)
        self.edit_Tp.setDecimals(4)
        self.edit_Tp.setMaximum(999999999.99)
        self.edit_Tp.setObjectName(_fromUtf8("edit_Tp"))
        self.paramsGridLayout_1.addWidget(self.edit_Tp, 2, 1, 1, 1)
        self.edit_r0 = QtGui.QDoubleSpinBox(self.gridLayoutWidget_3)
        self.edit_r0.setDecimals(4)
        self.edit_r0.setMaximum(999999999.99)
        self.edit_r0.setObjectName(_fromUtf8("edit_r0"))
        self.paramsGridLayout_1.addWidget(self.edit_r0, 1, 1, 1, 1)
        self.label_ms = QtGui.QLabel(self.gridLayoutWidget_3)
        self.label_ms.setObjectName(_fromUtf8("label_ms"))
        self.paramsGridLayout_1.addWidget(self.label_ms, 3, 0, 1, 1)
        self.label_r0_raz = QtGui.QLabel(self.gridLayoutWidget_3)
        self.label_r0_raz.setObjectName(_fromUtf8("label_r0_raz"))
        self.paramsGridLayout_1.addWidget(self.label_r0_raz, 1, 2, 1, 1)
        self.label_ms_raz = QtGui.QLabel(self.gridLayoutWidget_3)
        self.label_ms_raz.setObjectName(_fromUtf8("label_ms_raz"))
        self.paramsGridLayout_1.addWidget(self.label_ms_raz, 3, 2, 1, 1)
        self.horizontalLayoutWidget = QtGui.QWidget(self.paramsBox)
        self.horizontalLayoutWidget.setGeometry(QtCore.QRect(10, 190, 415, 41))
        self.horizontalLayoutWidget.setObjectName(_fromUtf8("horizontalLayoutWidget"))
        self.speedHorizontalLayout = QtGui.QHBoxLayout(self.horizontalLayoutWidget)
        self.speedHorizontalLayout.setMargin(0)
        self.speedHorizontalLayout.setObjectName(_fromUtf8("speedHorizontalLayout"))
        self.label_speed = QtGui.QLabel(self.horizontalLayoutWidget)
        self.label_speed.setObjectName(_fromUtf8("label_speed"))
        self.speedHorizontalLayout.addWidget(self.label_speed)
        self.edit_speed = QtGui.QDoubleSpinBox(self.horizontalLayoutWidget)
        self.edit_speed.setMinimum(-9999990.0)
        self.edit_speed.setMaximum(999999999.99)
        self.edit_speed.setProperty("value", 0.0)
        self.edit_speed.setObjectName(_fromUtf8("edit_speed"))
        self.speedHorizontalLayout.addWidget(self.edit_speed)
        self.horizontalLayoutWidget_2 = QtGui.QWidget(self.paramsBox)
        self.horizontalLayoutWidget_2.setGeometry(QtCore.QRect(530, 190, 241, 41))
        self.horizontalLayoutWidget_2.setObjectName(_fromUtf8("horizontalLayoutWidget_2"))
        self.timeHorizontalLayout = QtGui.QHBoxLayout(self.horizontalLayoutWidget_2)
        self.timeHorizontalLayout.setMargin(0)
        self.timeHorizontalLayout.setObjectName(_fromUtf8("timeHorizontalLayout"))
        self.label_time_integr = QtGui.QLabel(self.horizontalLayoutWidget_2)
        self.label_time_integr.setObjectName(_fromUtf8("label_time_integr"))
        self.timeHorizontalLayout.addWidget(self.label_time_integr)
        self.edit_time_integration = QtGui.QSpinBox(self.horizontalLayoutWidget_2)
        self.edit_time_integration.setMinimum(-99)
        self.edit_time_integration.setMaximum(99)
        self.edit_time_integration.setObjectName(_fromUtf8("edit_time_integration"))
        self.timeHorizontalLayout.addWidget(self.edit_time_integration)
        self.gridLayoutWidget_4 = QtGui.QWidget(self.scrollAreaWidgetContents)
        self.gridLayoutWidget_4.setGeometry(QtCore.QRect(10, 10, 1161, 451))
        self.gridLayoutWidget_4.setObjectName(_fromUtf8("gridLayoutWidget_4"))
        self.gridLayout = QtGui.QGridLayout(self.gridLayoutWidget_4)
        self.gridLayout.setContentsMargins(1, -1, -1, -1)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.widget_3 = QtGui.QWidget(self.gridLayoutWidget_4)
        self.widget_3.setObjectName(_fromUtf8("widget_3"))
        self.frame3D = QtGui.QFrame(self.widget_3)
        self.frame3D.setGeometry(QtCore.QRect(100, 0, 381, 221))
        self.frame3D.setObjectName(_fromUtf8("frame3D"))
        self.ySlider = QtGui.QSlider(self.widget_3)
        self.ySlider.setGeometry(QtCore.QRect(530, 10, 21, 191))
        self.ySlider.setMaximum(180)
        self.ySlider.setProperty("value", 90)
        self.ySlider.setOrientation(QtCore.Qt.Vertical)
        self.ySlider.setObjectName(_fromUtf8("ySlider"))
        self.gridLayout.addWidget(self.widget_3, 1, 2, 1, 1)
        self.tab_i = QtGui.QWidget(self.gridLayoutWidget_4)
        self.tab_i.setObjectName(_fromUtf8("tab_i"))
        self.gridLayout.addWidget(self.tab_i, 1, 1, 1, 1)
        self.tab_a = QtGui.QWidget(self.gridLayoutWidget_4)
        self.tab_a.setObjectName(_fromUtf8("tab_a"))
        self.gridLayout.addWidget(self.tab_a, 0, 2, 1, 1)
        self.tab_e = QtGui.QWidget(self.gridLayoutWidget_4)
        self.tab_e.setObjectName(_fromUtf8("tab_e"))
        self.gridLayout.addWidget(self.tab_e, 0, 1, 1, 1)
        self.timeSlider = QtGui.QSlider(self.scrollAreaWidgetContents)
        self.timeSlider.setGeometry(QtCore.QRect(20, 470, 533, 21))
        self.timeSlider.setMaximum(1000)
        self.timeSlider.setSingleStep(1)
        self.timeSlider.setPageStep(10)
        self.timeSlider.setOrientation(QtCore.Qt.Horizontal)
        self.timeSlider.setObjectName(_fromUtf8("timeSlider"))
        self.xSlider = QtGui.QSlider(self.scrollAreaWidgetContents)
        self.xSlider.setGeometry(QtCore.QRect(710, 470, 351, 21))
        self.xSlider.setMaximum(359)
        self.xSlider.setProperty("value", 180)
        self.xSlider.setOrientation(QtCore.Qt.Horizontal)
        self.xSlider.setObjectName(_fromUtf8("xSlider"))
        self.label = QtGui.QLabel(self.scrollAreaWidgetContents)
        self.label.setGeometry(QtCore.QRect(20, 0, 161, 41))
        self.label.setObjectName(_fromUtf8("label"))
        self.label_2 = QtGui.QLabel(self.scrollAreaWidgetContents)
        self.label_2.setGeometry(QtCore.QRect(20, 220, 161, 41))
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.label_3 = QtGui.QLabel(self.scrollAreaWidgetContents)
        self.label_3.setGeometry(QtCore.QRect(600, 0, 161, 41))
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.scrollArea.setWidget(self.scrollAreaWidgetContents)
        self.gridLayout_2.addWidget(self.scrollArea, 0, 0, 1, 1)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtGui.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1191, 27))
        self.menubar.setObjectName(_fromUtf8("menubar"))
        self.menu = QtGui.QMenu(self.menubar)
        self.menu.setObjectName(_fromUtf8("menu"))
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtGui.QStatusBar(MainWindow)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))
        MainWindow.setStatusBar(self.statusbar)
        self.action = QtGui.QAction(MainWindow)
        self.action.setObjectName(_fromUtf8("action"))
        self.action_save = QtGui.QAction(MainWindow)
        self.action_save.setObjectName(_fromUtf8("action_save"))
        self.action_load = QtGui.QAction(MainWindow)
        self.action_load.setObjectName(_fromUtf8("action_load"))
        self.menu.addAction(self.action_save)
        self.menu.addAction(self.action_load)
        self.menubar.addAction(self.menu.menuAction())

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow", None))
        self.infoBox.setTitle(_translate("MainWindow", "????????", None))
        self.label_info_e.setText(_translate("MainWindow", "e", None))
        self.label_info_i.setText(_translate("MainWindow", "i", None))
        self.label_info_a.setText(_translate("MainWindow", "a", None))
        self.manageBox.setTitle(_translate("MainWindow", "????????????????????", None))
        self.startButton.setText(_translate("MainWindow", "??????????????????", None))
        self.stopButton.setText(_translate("MainWindow", "????????????????????", None))
        self.solveButton.setText(_translate("MainWindow", "????????????????????", None))
        self.paramsBox.setTitle(_translate("MainWindow", "??????????????????", None))
        self.label_a.setText(_translate("MainWindow", "?????????????? ?????????????? (a)", None))
        self.label_ecc.setText(_translate("MainWindow", "???????????????????????????? (e)", None))
        self.label_i.setText(_translate("MainWindow", "???????????????????? (i)", None))
        self.label_i_2.setText(_translate("MainWindow", "????????.", None))
        self.label_a_raz.setText(_translate("MainWindow", "?? 1e+6 ??", None))
        self.label_Tp_raz.setText(_translate("MainWindow", "??????", None))
        self.label_mp.setText(_translate("MainWindow", "?????????? ??????????????", None))
        self.label_Tp.setText(_translate("MainWindow", "???????????? ???????????????? (Tp)", None))
        self.label_mp_raz.setText(_translate("MainWindow", "?? 1e+24 ????", None))
        self.label_r0.setText(_translate("MainWindow", "???????????? ??????????????", None))
        self.label_ms.setText(_translate("MainWindow", "?????????? ????????????????", None))
        self.label_r0_raz.setText(_translate("MainWindow", "?? 1e+6 ??", None))
        self.label_ms_raz.setText(_translate("MainWindow", "?? 1e+20 ????", None))
        self.label_speed.setText(_translate("MainWindow", "????????. ???????????????? ???????????????? 86400 ?? 10^", None))
        self.label_time_integr.setText(_translate("MainWindow", "??????????. ?????????????? 10^", None))
        self.label.setText(_translate("MainWindow", "????????????????????????????", None))
        self.label_2.setText(_translate("MainWindow", "????????????????????", None))
        self.label_3.setText(_translate("MainWindow", "?????????????? ??????????????", None))
        self.menu.setTitle(_translate("MainWindow", "????????", None))
        self.action.setText(_translate("MainWindow", "?????????????????? ??????????????????", None))
        self.action_save.setText(_translate("MainWindow", "?????????????????? ??????????????????", None))
        self.action_load.setText(_translate("MainWindow", "?????????????????? ??????????????????", None))

