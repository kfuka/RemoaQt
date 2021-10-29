# - * - coding: utf8 - * -
from os import remove
from PyQt5.QtCore import center
import pydicom
import sys
import pyqtgraph as pg
import numpy as np
from PyQt5.QtCore import Qt
import matplotlib.pyplot as plt

from PyQt5.QtWidgets import QAction, QFileDialog, QHBoxLayout, QMainWindow, QApplication, QAction, QPushButton, QSlider, QVBoxLayout, QWidget, qApp, QGridLayout
from pathlib import Path

from pyqtgraph.graphicsItems.ROI import ROI


pg.setConfigOptions(imageAxisOrder='row-major')

class ImageView(pg.ImageView):
    def __init__(self, *args, **kwargs):
        pg.ImageView.__init__(self, *args, **kwargs)
        

class MyWindow(QMainWindow): # QWidgetクラスを使用します。

    def __init__(self):
        super().__init__()
        self.title = 'Qmain window'
        self.width = 1000
        self.height = 800
        self.initUI()
        self.roi_corrections = {} # dict in dictにしたい[slice][roinum]にできたらいいな
        self.current_slice = 0
        self.rois = []


    def initUI(self):
        centerWidget = QWidget()
        self.setCentralWidget(centerWidget)
        self.layout = QGridLayout()
        centerWidget.setLayout(self.layout)

        exitAction = QAction(' &Exit', self)
        exitAction.setShortcut('Ctrl+Q')
        exitAction.setStatusTip('Close GUI')
        exitAction.triggered.connect(qApp.quit)

        fileOpen = QAction(' &Open', self)
        fileOpen.triggered.connect(self.selectfile)

        # sliderをつけるところ
        self.sld = QSlider(Qt.Horizontal, self)
        self.sld.setRange(0,29)
        self.sld.setFocusPolicy(Qt.NoFocus)
        self.sld.valueChanged.connect(self.updateslice)


        menubar = self.menuBar()
        fileMenu = menubar.addMenu(' &File')
        fileMenu.addAction(fileOpen)
        fileMenu.addAction(exitAction)

        self.setWindowTitle(self.title)
        self.setGeometry(0, 0, self.width, self.height)
        self.statusBar()

        # ImageViewをつけたところ
        self.imv = ImageView()
        self.imv.imageItem.mouseClickEvent = self.mouse_click
        self.layout.addWidget(self.imv, 0,0)
        self.imv.ui.roiBtn.hide()
        self.imv.ui.menuBtn.hide() 

        self.imv2 = ImageView()
        self.imv2.imageItem.mouseClickEvent = self.mouse_click
        self.layout.addWidget(self.imv2, 0,1)
        self.imv2.ui.roiBtn.hide()
        self.imv2.ui.menuBtn.hide()



        self.layout.addWidget(self.sld, 1,0)

        # 計算ボタンをつけた．これいる？
        btn = QPushButton(text="calc dicom")
        btn.clicked.connect(self.load_dicom)
        self.layout.addWidget(btn, 1,1)

        self.show()

    def update_roi_correction(self):
        print("handle")
    
    def updateslice(self):
        # ここが2回呼ばれている？
        current_slice = self.sld.value()

        self.imv.setImage(self.array[int(current_slice), :, :],
                          autoRange=False,
                          autoLevels=False,
                          autoHistogramRange=False)

        # ROIがあれば
        if current_slice in self.roi_corrections:
            for n, a_roi in enumerate(self.rois):
                a_roi.setPos(np.array(self.roi_corrections[current_slice][n]) - 15, update=False)              

        


    def selectfile(self):
        directory = Path("")
        fname = QFileDialog.getOpenFileName(self, "Open file")
        self.file = fname[0]
        # print(self.file)
        self.load_dicom()
    
    def load_dicom(self):
        dicom = pydicom.dcmread(self.file)
        self.array = np.array(dicom.pixel_array)
        self.imv.setImage(self.array[0, :, :],
                          autoRange=False,
                          autoLevels=True,
                          autoHistogramRange=False,
                          levelMode="mono")
        
        # self.imv.autoLevels()
        hist = self.imv.getImageItem().getHistogram()
        mid_x = hist[0][np.argmax(hist[1][:len(hist[1])//2])]
        min_x = hist[0][0]
        max_x = (mid_x - min_x) * 3 + mid_x
        self.imv.setLevels(min_x,max_x)

    def mouse_click(self, event):
        current_slice = self.sld.value()

        # self.roi_correctionsに「今のスライス」のキーがあればROIを追加
        # なければキーを追加して値を入れる． [x, y]のリスト．
        # self.imvにROIを追加する．イメージの更新はしなくてよさそう，

        if current_slice in self.roi_corrections:
            self.roi_corrections[current_slice].append([event.pos().x(), event.pos().y()])
        else:
            self.roi_corrections[current_slice] = [[event.pos().x(), event.pos().y()]]
        
        roi = pg.RectROI(np.array([event.pos().x(), event.pos().y()])- 15,
                        [30,30], pen=(len(self.roi_corrections[current_slice]),9)) 
        roi.resizable = False
        roi.rotatable = False   
        roi.sigRegionChangeFinished.connect(self.roiMove)    
        self.imv.addItem(roi)
        self.rois.append(roi)
        # for the_roi in self.rois:
        #     the_roi.sigRegionChanged.connect(self.roiMove)

    def roiMove(self):
        # もしROIの動きを感知したら
        current_slice = self.sld.value()
        if current_slice in self.roi_corrections:
            self.roi_corrections.pop(current_slice)
        roi_update = []
        # self.roi_corrections[current_slice] = []
        for roi in self.rois:
            # そのスライスにある
            # self.roi_corrections[current_slice].append([roi.pos().x() + 15, roi.pos().y() + 15])
            roi_update.append([roi.pos().x() + 15, roi.pos().y() + 15])
        self.roi_corrections[current_slice] = roi_update



def main():
    app = QApplication(sys.argv)
    gui = MyWindow()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
