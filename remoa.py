# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'remoa.ui'
#
# Created by: PyQt5 UI code generator 5.15.2
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtWidgets
from pyqtgraph import PlotWidget, ImageView


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1000, 800)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setMouseTracking(False)
        self.centralwidget.setAutoFillBackground(True)
        self.centralwidget.setObjectName("centralwidget")
        self.graphicsView = ImageView(self.centralwidget)
        self.graphicsView.setGeometry(QtCore.QRect(10, 10, 480, 340))
        self.graphicsView.setObjectName("graphicsView")
        self.graphicsView.ui.roiBtn.hide()
        self.graphicsView.ui.menuBtn.hide()
        self.graphicsView_2 = ImageView(self.centralwidget)
        self.graphicsView_2.setGeometry(QtCore.QRect(510, 10, 480, 340))
        self.graphicsView_2.setObjectName("graphicsView_2")
        self.graphicsView_2.ui.roiBtn.hide()
        self.graphicsView_2.ui.menuBtn.hide()
        self.graphicsView_3 = PlotWidget(self.centralwidget)
        self.graphicsView_3.setGeometry(QtCore.QRect(118, 380, 371, 250))
        self.graphicsView_3.setObjectName("graphicsView_3")
        self.graphicsView_4 = PlotWidget(self.centralwidget)
        self.graphicsView_4.setGeometry(QtCore.QRect(510, 380, 371, 250))
        self.graphicsView_4.setObjectName("graphicsView_4")
        self.horizontalSlider = QtWidgets.QSlider(self.centralwidget)
        self.horizontalSlider.setGeometry(QtCore.QRect(200, 350, 601, 30))
        self.horizontalSlider.setOrientation(QtCore.Qt.Horizontal)
        self.horizontalSlider.setObjectName("horizontalSlider")
        self.tableView = QtWidgets.QTableView(self.centralwidget)
        self.tableView.setGeometry(QtCore.QRect(120, 640, 760, 120))
        self.tableView.setObjectName("tableView")
        self.verticalLayoutWidget = QtWidgets.QWidget(self.centralwidget)
        self.verticalLayoutWidget.setGeometry(QtCore.QRect(10, 390, 101, 361))
        self.verticalLayoutWidget.setObjectName("verticalLayoutWidget")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.verticalLayoutWidget)
        self.verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout.setObjectName("verticalLayout")
        self.pushButton = QtWidgets.QPushButton(self.verticalLayoutWidget)
        self.pushButton.setObjectName("pushButton")
        self.verticalLayout.addWidget(self.pushButton)
        self.pushButton_3 = QtWidgets.QPushButton(self.verticalLayoutWidget)
        self.pushButton_3.setObjectName("pushButton_3")
        self.verticalLayout.addWidget(self.pushButton_3)
        self.pushButton_4 = QtWidgets.QPushButton(self.verticalLayoutWidget)
        self.pushButton_4.setObjectName("pushButton_4")
        self.verticalLayout.addWidget(self.pushButton_4)
        self.pushButton_5 = QtWidgets.QPushButton(self.verticalLayoutWidget)
        self.pushButton_5.setObjectName("pushButton_5")
        self.verticalLayout.addWidget(self.pushButton_5)
        self.pushButton_6 = QtWidgets.QPushButton(self.verticalLayoutWidget)
        self.pushButton_6.setObjectName("pushButton_6")
        self.verticalLayout.addWidget(self.pushButton_6)

        self.pushButton_7 = QtWidgets.QPushButton(self.verticalLayoutWidget)
        self.pushButton_7.setObjectName("pushButton_7")
        self.verticalLayout.addWidget(self.pushButton_7)

        # added right buttons
        self.verticalLayoutWidget2 = QtWidgets.QWidget(self.centralwidget)
        self.verticalLayoutWidget2.setGeometry(
            QtCore.QRect(890, 390, 101, 361))
        self.verticalLayoutWidget2.setObjectName("verticalLayoutWidget")
        self.verticalLayout2 = QtWidgets.QVBoxLayout(
            self.verticalLayoutWidget2)
        self.verticalLayout2.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout2.setObjectName("verticalLayout")
        self.pushButton_r1 = QtWidgets.QPushButton(self.verticalLayoutWidget2)
        self.pushButton_r1.setCheckable(True)
        self.pushButton_r1.setObjectName("pushButton_r1")
        self.pushButton_r1.setStyleSheet("background-color: green")
        self.verticalLayout2.addWidget(self.pushButton_r1)
        self.pushButton_r2 = QtWidgets.QPushButton(self.verticalLayoutWidget2)
        self.pushButton_r2.setCheckable(True)
        self.pushButton_r2.setStyleSheet("background-color: cyan")
        self.pushButton_r2.setObjectName("pushButton_r2")
        self.verticalLayout2.addWidget(self.pushButton_r2)
        self.pushButton_r3 = QtWidgets.QPushButton(self.verticalLayoutWidget2)
        self.pushButton_r3.setCheckable(True)
        self.pushButton_r3.setStyleSheet("background-color: magenta")
        self.pushButton_r3.setObjectName("pushButton_r3")
        self.verticalLayout2.addWidget(self.pushButton_r3)
        self.pushButton_r4 = QtWidgets.QPushButton(self.verticalLayoutWidget2)
        self.pushButton_r4.setCheckable(True)
        self.pushButton_r4.setStyleSheet("background-color: yellow")
        self.pushButton_r4.setObjectName("pushButton_r4")
        self.verticalLayout2.addWidget(self.pushButton_r4)

        self.pushButton_2 = QtWidgets.QPushButton(
            self.verticalLayoutWidget)
        self.pushButton_2.setObjectName("pushButton_2")
        self.verticalLayout.addWidget(self.pushButton_2)
        MainWindow.setCentralWidget(self.centralwidget)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)
        self.actionResp_Motion_Analyzer = QtWidgets.QAction(MainWindow)
        self.actionResp_Motion_Analyzer.setObjectName(
            "actionResp_Motion_Analyzer")

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.pushButton.setText(_translate("MainWindow", "Open"))
        self.pushButton_3.setText(_translate("MainWindow", "Calc."))
        self.pushButton_4.setText(_translate("MainWindow", "Save"))
        self.pushButton_5.setText(_translate("MainWindow", "Report"))
        self.pushButton_6.setText(_translate("MainWindow", "Clear"))
        self.pushButton_2.setText(_translate("MainWindow", "Close"))
        self.pushButton_7.setText(_translate("MainWindow", "Patient Dir."))
        self.pushButton_r1.setText(_translate("MainWindow", "Edit #1"))
        self.pushButton_r2.setText(_translate("MainWindow", "Edit #2"))
        self.pushButton_r3.setText(_translate("MainWindow", "Edit #3"))
        self.pushButton_r4.setText(_translate("MainWindow", "Edit #4"))
        self.actionResp_Motion_Analyzer.setText(
            _translate("MainWindow", "Resp. Motion Analyzer"))
