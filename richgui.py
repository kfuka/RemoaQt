import sys
import time
import os
import tkinter
from PyQt5 import QtWidgets, QtCore
from PyQt5.QtGui import QBrush, QKeySequence
from PyQt5.QtCore import QVariant, Qt
import pydicom
from remoa import Ui_MainWindow
import numpy as np
from open_file import openDirectory, get_wave_dicoms, openDicoms
import pyqtgraph as pg
import pyqtgraph.exporters
import wave_analysis
from csrt_multi_tracker import track_MIL
from scipy import signal
import configparser
import subprocess
from unite_class import build_report
# import tkinter as tk
from tkinter import filedialog
from PIL import Image

pg.setConfigOptions(imageAxisOrder='row-major')
spacing_correction_factor = 1.5 * 1550 / 2111.63

# edited

"""
.uiファイルをつくったら
pyuic5 remoa.ui -o remoa.py
でpythonファイルにしてしまう
"""

config = configparser.ConfigParser()
config.read("./config.ini")
data_base_folder = config.get("settings", "database_path")
color_waves = ("g", "c", "m", "y", "b", "k", "w", "g", "c", "m", "y", "b", "k")
roi_size = int(30)


class MainWindow(QtWidgets.QMainWindow, Ui_MainWindow):

    def __init__(self, parent=None, *args, obj=None, **kwargs):
        super(MainWindow, self).__init__(parent, *args, **kwargs)
        self.sub1 = None
        self.iscalculated = False
        self.setupUi(self)
        self.pushButton_2.clicked.connect(self.close)
        self.pushButton.clicked.connect(self.openfile)
        self.pushButton_6.clicked.connect(self.clear)
        self.pushButton_3.clicked.connect(self.calc)
        self.pushButton_4.clicked.connect(self.save_items)
        self.pushButton_7.clicked.connect(self.open_folder)
        self.pushButton_5.clicked.connect(self.create_report)
        self.horizontalSlider.valueChanged.connect(self.update)
        self.pushButton_r1.clicked.connect(self.edit1)
        self.pushButton_r2.clicked.connect(self.edit2)
        self.pushButton_r3.clicked.connect(self.edit3)
        self.pushButton_r4.clicked.connect(self.edit4)
        self.graphicsView.imageItem.mouseClickEvent = self.mouse_click1
        self.graphicsView_2.imageItem.mouseClickEvent = self.mouse_click2
        self.tableheaders = ["SI : V", "LR : V", "SI : H", "AP : H", "3D (mm)"]

        self.pushButton_3.setEnabled(False)
        self.pushButton_4.setEnabled(False)
        self.pushButton_7.setEnabled(False)
        self.pushButton_r1.setEnabled(False)
        self.pushButton_r2.setEnabled(False)
        self.pushButton_r3.setEnabled(False)
        self.pushButton_r4.setEnabled(False)

        self.roi_corrections1 = {}
        self.roi_corrections2 = {}
        self.rois1 = []
        self.rois2 = []
        self.marker_chase1 = []
        self.marker_chase2 = []
        self.mwave1 = []
        self.mwave2 = []

        self.inedit = False
        self.edit_num = 0

        self.key_left = QtWidgets.QShortcut(QKeySequence("Left"), self)
        self.key_left.activated.connect(self.left_press)
        self.key_rt = QtWidgets.QShortcut(QKeySequence("Right"), self)
        self.key_rt.activated.connect(self.right_press)

    def edit1(self, state):
        if state:
            self.inedit = True
            self.edit_num = 1
        else:
            self.inedit = False

    def edit2(self, state):
        if state:
            self.inedit = True
            self.edit_num = 2
        else:
            self.inedit = False

    def edit3(self, state):
        if state:
            self.inedit = True
            self.edit_num = 3
        else:
            self.inedit = False

    def edit4(self, state):
        if state:
            self.inedit = True
            self.edit_num = 4
        else:
            self.inedit = False

    def open_folder(self):
        if self.id:
            subprocess.run('explorer {}'.format(data_base_folder + str(self.id)), cwd=os.getcwd())
        else:
            subprocess.run('explorer {}'.format(data_base_folder), cwd=os.getcwd())

    def left_press(self):
        current_slice = int(self.horizontalSlider.value())
        if current_slice == 0:
            new_slice = len(self.array1) - 2
        else:
            new_slice = current_slice - 1
        self.horizontalSlider.setValue(new_slice)
        self.update()

    def right_press(self):
        current_slice = int(self.horizontalSlider.value())
        if current_slice == len(self.array1) - 2:
            new_slice = 0
        else:
            new_slice = current_slice + 1
        self.horizontalSlider.setValue(new_slice)
        self.update()

    def calc(self):
        # ret = multi_csrt(self.roi_corrections1[0], self.array1)
        self.iscalculated = True
        for i in range(len(self.rois1)):
            trackroi = [int(self.roi_corrections1[0][i][0] - roi_size / 2), int(self.roi_corrections1[0][i][1] - roi_size / 2), roi_size, roi_size]
            ret = track_MIL(trackroi, self.array1, "aaa", "bbb", 0)
            for a_slice, data in enumerate(ret):
                if a_slice in self.roi_corrections1:
                    self.roi_corrections1[a_slice][i] = [data[0], data[1]]
                else:
                    self.roi_corrections1[a_slice] = {}
                    self.roi_corrections1[a_slice][i] = [data[0], data[1]]
            maxx, maxy, wave, chased = self.wave_calc(ret, self.pixel_spacing1)
            self.marker_chase1.append(chased)
            self.mwave1.append(wave)

        for i in range(len(self.rois2)):
            trackroi = [int(self.roi_corrections2[0][i][0] - roi_size / 2), int(self.roi_corrections2[0][i][1] - roi_size / 2), roi_size, roi_size]
            ret = track_MIL(trackroi, self.array2, "aaa", "bbb", 0)
            for a_slice, data in enumerate(ret):
                if a_slice in self.roi_corrections2:
                    self.roi_corrections2[a_slice][i] = [data[0], data[1]]
                else:
                    self.roi_corrections2[a_slice] = {}
                    self.roi_corrections2[a_slice][i] = [data[0], data[1]]
            maxx, maxy, wave, chased = self.wave_calc(ret, self.pixel_spacing2)
            self.marker_chase2.append(chased)
            self.mwave2.append(wave)

        if self.mwave1 != []:
            for i in range(len(self.mwave1)):
                self.graphicsView_3.plot(self.wave_time1, self.mwave1[i], pen=pg.mkPen(color=color_waves[i]), antialias=True, name=str(i + 1))

        if self.mwave2 != []:
            for i in range(len(self.mwave2)):
                self.graphicsView_4.plot(self.wave_time2, self.mwave2[i], pen=pg.mkPen(color=color_waves[i]), antialias=True, name=str(i + 1))

        self.generate_table()

        self.pushButton_4.setEnabled(True)
        self.pushButton_6.setEnabled(True)
        self.pushButton_7.setEnabled(True)
        self.pushButton_r1.setEnabled(True)
        self.pushButton_r2.setEnabled(True)
        self.pushButton_r3.setEnabled(True)
        self.pushButton_r4.setEnabled(True)

    def generate_table(self):
        self.max_x1, self.max_x2, self.max_y1, self.max_y2, self.max_3d = self.get_maxshift_from_marker_chase()
        
        data = []
        for i in range(len(self.max_x1)):
            data.append([str(round(self.max_y1[i], 2)), str(round(self.max_x1[i], 2)), str(round(self.max_y2[i], 2)), str(round(self.max_x2[i], 2)), str(round(self.max_3d[i], 2))])
        self.model = TableModel(data, self.tableheaders)
        self.tableView.setModel(self.model)
        self.tableView.resizeRowsToContents()

    def get_maxshift_from_marker_chase(self):
        max_x1, max_y1, max_x2, max_y2 = [], [], [], []
        for j,i in enumerate(self.marker_chase1):
            x_shifts, y_shifts = i
            x2_shifts, y2_shifts = self.marker_chase2[j]
            max_x1.append((np.max(x_shifts) - np.min(x_shifts)) * self.pixel_spacing1[0])
            max_y1.append((np.max(y_shifts) - np.min(y_shifts)) * self.pixel_spacing1[1])
            max_x2.append((np.max(x2_shifts) - np.min(x2_shifts)) * self.pixel_spacing1[0])
            max_y2.append((np.max(y2_shifts) - np.min(y2_shifts)) * self.pixel_spacing1[1])

        LR_from_center = []
        AP_from_center = []
        for x1 in self.marker_chase1:
            LR_from_center.append((756.0 - np.average(x1[0])) * self.pixel_spacing1[0])
        for x2 in self.marker_chase2:
            AP_from_center.append((np.average(x2[0]) - 756.0) * self.pixel_spacing2[0])

        c_max_x1, c_max_y1, c_max_x2, c_max_y2, c_max_3d = [], [], [], [], []
        for n, the_x1 in enumerate(max_x1):
            c_max_x1.append(self.LRAP_correction(the_x1, AP_from_center[n]))
            c_max_y1.append(self.LRAP_correction(max_y1[n], AP_from_center[n]))
            c_max_x2.append(self.LRAP_correction(max_x2[n], LR_from_center[n]))
            c_max_y2.append(self.LRAP_correction(max_y2[n], LR_from_center[n]))
            c_max_3d.append(np.sqrt((self.LRAP_correction(the_x1, AP_from_center[n])) ** 2 + (self.LRAP_correction(max_x2[n], LR_from_center[n])) ** 2 + (max([self.LRAP_correction(max_y1[n], AP_from_center[n]), self.LRAP_correction(max_y2[n], LR_from_center[n])])) **2))

        return c_max_x1, c_max_x2, c_max_y1, c_max_y2, c_max_3d

    def LRAP_correction(self, val_iso, distance):
        return (1550.0 + distance) / 1550.0 * val_iso

    def wave_calc(self, tracker_returned, pixel_spacing):
        x_shifts = np.array([it[0] for it in tracker_returned])
        y_shifts = np.array([it[1] for it in tracker_returned])
        x_shifts = signal.savgol_filter(x_shifts, 7, 5, mode="nearest")
        y_shifts = signal.savgol_filter(y_shifts, 7, 5, mode="nearest")
        max_x_move = (np.max(x_shifts) - np.min(x_shifts)) * pixel_spacing[0]
        max_y_move = (np.max(y_shifts) - np.min(y_shifts)) * pixel_spacing[1]
        wave = np.array(y_shifts - np.min(y_shifts)) / np.max(y_shifts - np.min(y_shifts)) * 100
        return max_x_move, max_y_move, np.delete(wave, 0), [x_shifts, y_shifts]

    def clear(self):
        for roi in self.rois1:
            self.graphicsView.removeItem(roi)

        for roi in self.rois2:
            self.graphicsView_2.removeItem(roi)

        self.roi_corrections1 = {}
        self.roi_corrections2 = {}
        self.rois1 = []
        self.rois2 = []
        self.marker_chase1 = []
        self.marker_chase2 = []
        self.mwave1, self.mwave2 = [], []
        self.iscalculated = False
        self.pushButton_4.setEnabled(False)
        self.pushButton_7.setEnabled(False)
        self.pushButton_r1.setEnabled(False)
        self.pushButton_r2.setChecked(False)
        self.pushButton_r3.setChecked(False)
        self.pushButton_r4.setChecked(False)
        self.pushButton_r1.setChecked(False)
        self.pushButton_r2.setEnabled(False)
        self.pushButton_r3.setEnabled(False)
        self.pushButton_r4.setEnabled(False)
        self.init_graphics()
        self.generate_table()

    def mouse_click1(self, event):
        current_slice = self.horizontalSlider.value()

        # self.roi_correctionsに「今のスライス」のキーがあればROIを追加
        # なければキーを追加して値を入れる． [x, y]のリスト．
        # self.imvにROIを追加する．イメージの更新はしなくてよさそう，
        if current_slice == 0 and not self.iscalculated:
            if self.roi_corrections1 == {}:
                roi_num = 0
                self.roi_corrections1[current_slice] = {}
            else:
                roi_num = len(self.roi_corrections1[current_slice])

            self.roi_corrections1[current_slice][roi_num] = [event.pos().x(),
                                                             event.pos().y()]

            roi = pg.RectROI(np.array([event.pos().x(), event.pos().y()]) - 15,
                             [30, 30],
                             pen=pg.mkPen(color_waves[len(self.roi_corrections1[current_slice]) - 1]),
                             rotatable=False, resizable=False)
            roi.resizable = False
            roi.rotatable = False
            roi.sigRegionChangeFinished.connect(self.roiMove1)
            self.graphicsView.addItem(roi)
            self.rois1.append(roi)
        else:
            if self.inedit and self.iscalculated:
                for n, roi in enumerate(self.rois1):
                    if n == self.edit_num - 1:
                        self.roi_corrections1[current_slice + 1][n] = \
                             [event.pos().x() - 15, event.pos().y() - 15]
            self.graphicsView_3.getPlotItem().clearPlots()
            # self.update()
            self.mwave1, maxx, maxy, self.marker_chase1 = \
                self.update_mwave(self.rois1, self.mwave1,
                                  self.roi_corrections1, self.pixel_spacing1)
            self.plot_wave1()
            if self.mwave1 != []:
                for i in range(len(self.mwave1)):
                    self.graphicsView_3.plot(self.wave_time1, self.mwave1[i],
                                             pen=pg.mkPen(color=color_waves[i]),
                                             antialias=True, name=str(i + 1))
            self.generate_table()
            if current_slice == 28:
                self.update()
                self.inedit = False

            else:
                self.update(_=0, plusnum=1)
                self.horizontalSlider.setValue(current_slice + 1)

    def roiMove1(self):
        # もしROIの動きを感知したら
        current_slice = self.horizontalSlider.value() + 1
        for n, roi in enumerate(self.rois1):
            self.roi_corrections1[current_slice][n] = [roi.pos().x(),
                                                       roi.pos().y()]
        self.graphicsView_3.getPlotItem().clearPlots()
        self.mwave1, maxx, maxy, self.marker_chase1 = \
            self.update_mwave(self.rois1,
                              self.mwave1,
                              self.roi_corrections1,
                              self.pixel_spacing1)
        self.plot_wave1()
        if self.mwave1 != []:
            for i in range(len(self.mwave1)):
                self.graphicsView_3.plot(self.wave_time1, self.mwave1[i],
                                         pen=pg.mkPen(color=color_waves[i]),
                                         antialias=True,
                                         name=str(i + 1))
        self.generate_table()

    def update_mwave(self, rois, mwave, roi_corrections, pixel_spacing):
        mwave = []
        maxx, maxy = [], [], 
        chase = []
        for n, roi in enumerate(rois):
            x_shifts, y_shifts = [], []
            for cslice in roi_corrections:
                x_shifts.append(roi_corrections[cslice][n][0])
                y_shifts.append(roi_corrections[cslice][n][1])
            x_shifts = signal.savgol_filter(x_shifts, 7, 5, mode="nearest")
            y_shifts = signal.savgol_filter(y_shifts, 7, 5, mode="nearest")
            chase.append([x_shifts, y_shifts])
            max_x_move = (np.max(x_shifts) - np.min(x_shifts))\
                * pixel_spacing[0]
            max_y_move = (np.max(y_shifts) - np.min(y_shifts))\
                * pixel_spacing[1]
            maxx.append(max_x_move)
            maxy.append(max_y_move)
            wave = np.array(y_shifts - np.min(y_shifts))\
                / np.max(y_shifts - np.min(y_shifts)) * 100
            mwave.append(np.delete(wave, 0))
        return mwave, maxx, maxy, chase

    def mouse_click2(self, event):
        current_slice = self.horizontalSlider.value()

        # self.roi_correctionsに「今のスライス」のキーがあればROIを追加
        # なければキーを追加して値を入れる． [x, y]のリスト．
        # self.imvにROIを追加する．イメージの更新はしなくてよさそう，
        if current_slice == 0 and not self.iscalculated:
            if self.roi_corrections2 == {}:
                roi_num = 0
                self.roi_corrections2[current_slice] = {}
            else:
                roi_num = len(self.roi_corrections2[current_slice])

            self.roi_corrections2[current_slice][roi_num] = [event.pos().x(), event.pos().y()]
            
            roi = pg.RectROI(np.array([event.pos().x(), event.pos().y()]) - 15,
                            [30,30], pen=pg.mkPen(color_waves[len(self.roi_corrections2[current_slice]) - 1])) 
            roi.resizable = False
            roi.rotatable = False   
            roi.sigRegionChangeFinished.connect(self.roiMove2)    
            self.graphicsView_2.addItem(roi)
            self.rois2.append(roi)
        else:
            if self.inedit and self.iscalculated:
                for n, roi in enumerate(self.rois2):
                    if n == self.edit_num - 1:
                        self.roi_corrections2[current_slice + 1][n] = \
                             [event.pos().x() - 15, event.pos().y() - 15]
            self.graphicsView_4.getPlotItem().clearPlots()
            # self.update()
            self.mwave2, maxx, maxy, self.marker_chase2 = \
                self.update_mwave(self.rois2, self.mwave2,
                                  self.roi_corrections2, self.pixel_spacing2)
            self.plot_wave2()
            if self.mwave2 != []:
                for i in range(len(self.mwave2)):
                    self.graphicsView_4.plot(self.wave_time2, self.mwave2[i],
                                             pen=pg.mkPen(color=color_waves[i]),
                                             antialias=True, name=str(i + 1))
            self.generate_table()
            if current_slice == 28:
                self.update()
                self.inedit = False

            else:
                self.update(_=0, plusnum=1)
                self.horizontalSlider.setValue(current_slice + 1)

    def roiMove2(self):
        # # もしROIの動きを感知したら

        current_slice = self.horizontalSlider.value() + 1
        for n, roi in enumerate(self.rois2):
            self.roi_corrections2[current_slice][n] = [roi.pos().x(), roi.pos().y()]
        # if self.mwave1 != []:
        #     for i in range(len(self.mwave1)):
        #         self.graphicsView_3.setData(self.wave_time1, self.mwave1[i], pen=pg.mkPen(color=color_waves[i]), antialias=True, name=str(i))
        self.graphicsView_4.getPlotItem().clearPlots()
        self.mwave2, maxx, maxy, self.marker_chase2 = self.update_mwave(self.rois2, self.mwave2, self.roi_corrections2, self.pixel_spacing2)
        self.plot_wave2()
        if self.mwave2 != []:
            for i in range(len(self.mwave2)):
                self.graphicsView_4.plot(self.wave_time2, self.mwave2[i], pen=pg.mkPen(color=color_waves[i]), antialias=True, name=str(i + 1))
        self.generate_table()

    def openfile(self):
        try:
            self.clear()
            self.__init__()
        except:
            pass
        getdicomclass = openDirectory()
        dicom_dir = getdicomclass.get()
        folder_dict = get_wave_dicoms(dicom_dir)
        self.sub1 = openDicoms(folder_dict, parent=self)
        self.sub1.show()
        self.pushButton_3.setEnabled(True)

    def get_open_dicoms(self, dicoms):
        # ここで開くべきDicom Filesを返してもらう．
        # dictionalyになっている．
        # [file, [time, UID, beam direction]]
        # 表示するところまでこの関数で行う
        # 子クラスから呼び出されている

        self.open_dicom_list = dicoms
        for a in dicoms:
            if a[1][2] == "V":
                self.dicom_file_1 = a[0]
            elif a[1][2] == "H":
                self.dicom_file_2 = a[0]
            else:
                print(a)

        self.init_graphics()
        # self.key_left = QtWidgets.QShortcut(QKeySequence("Left"), self)
        # self.key_left.activated.connect(self.left_press)
        # self.key_rt = QtWidgets.QShortcut(QKeySequence("Right"), self)
        # self.key_rt.activated.connect(self.right_press)

    def init_graphics(self):
        self.dicom1 = pydicom.dcmread(self.dicom_file_1)
        self.dicom2 = pydicom.dcmread(self.dicom_file_2)
        self.id = self.dicom1.PatientID
        self.SOPUID = [self.dicom1[0x0008, 0x0018].value,
                       self.dicom2[0x0008, 0x0018].value]
        self.actime = [self.dicom1.AcquisitionTime,
                       self.dicom2.AcquisitionTime]
        self.study_date = [self.dicom1.StudyDate, self.dicom2.StudyDate]
        self.array1 = self.dicom1.pixel_array
        self.graphicsView.setImage(self.array1[1],
                                   autoRange=False,
                                   autoLevels=True,
                                   autoHistogramRange=False,
                                   levelMode="mono")
        self.array2 = self.dicom2.pixel_array
        self.graphicsView_2.setImage(self.array2[1],
                                     autoRange=False,
                                     autoLevels=True,
                                     autoHistogramRange=False,
                                     levelMode="mono")

        self.pixel_spacing1 = np.array([float(self.dicom1.PixelSpacing[0]),
                                        float(self.dicom1.PixelSpacing[1])]) \
                                        * spacing_correction_factor
        self.pixel_spacing2 = np.array([float(self.dicom2.PixelSpacing[0]),
                                        float(self.dicom2.PixelSpacing[1])]) * spacing_correction_factor

        self.horizontalSlider.setRange(0, len(self.array1) - 2)
        hist = self.graphicsView.getImageItem().getHistogram()
        mid_x = hist[0][np.argmax(hist[1][:len(hist[1])//2])]
        min_x = hist[0][0]
        max_x = (mid_x - min_x) * 3 + mid_x
        self.graphicsView.setLevels(min_x,max_x)
        hist = self.graphicsView_2.getImageItem().getHistogram()
        mid_x = hist[0][np.argmax(hist[1][:len(hist[1])//2])]
        min_x = hist[0][0]
        max_x = (mid_x - min_x) * 3 + mid_x
        self.graphicsView_2.setLevels(min_x,max_x)

        # wave をplotする
        self.wave1, self.wave_time1 = wave_analysis.wave_analysis(self.dicom1)
        self.wave2, self.wave_time2 = wave_analysis.wave_analysis(self.dicom2)
        self.plot_wave()
        self.line1.setValue(self.wave_time1[0])
        self.line2.setValue(self.wave_time2[0])

    def plot_wave(self):
        self.graphicsView_3.plotItem.clear()
        self.graphicsView_4.plotItem.clear()
        self.graphicsView_3.addLegend()
        self.graphicsView_4.addLegend()
        self.graphicsView_3.plot(self.wave_time1, self.wave1, pen=pg.mkPen(color="r"), antialias=True, name="Resp.", symbolBrush = "r", symbol="o", symbolPen="r", symbolSize=5)
        self.graphicsView_3.setLabel("bottom", text="time (sec.)")
        self.graphicsView_3.setLabel("left", text="Resp. Phase (%)")
        self.graphicsView_3.setXRange(np.min(self.wave_time1), np.max(self.wave_time1), padding=0.05)
        self.graphicsView_3.setYRange(0, 100, padding=0.05)
        self.graphicsView_4.plot(self.wave_time2, self.wave2, pen=pg.mkPen(color="r"), antialias=True, name="Resp.", symbolBrush = "r", symbol="o", symbolPen="r", symbolSize=5)
        self.graphicsView_4.setLabel("bottom", text="time (sec.)")
        self.graphicsView_4.setLabel("left", text="Resp. Phase (%)")
        self.graphicsView_4.setXRange(np.min(self.wave_time2), np.max(self.wave_time2), padding=0.05)
        self.graphicsView_4.setYRange(0, 100, padding=0.05)
        self.line1 = pg.InfiniteLine(pos=0, angle=90, pen=pg.mkPen(color="w", width=0.5))
        self.line2 = pg.InfiniteLine(pos=0, angle=90, pen=pg.mkPen(color="w", width=0.5))
        self.graphicsView_3.addItem(self.line1)
        self.graphicsView_4.addItem(self.line2)

    def plot_wave1(self):
        self.graphicsView_3.plot(self.wave_time1, self.wave1,
                                 pen=pg.mkPen(color="r"),
                                 antialias=True, name="Resp.",
                                 symbolBrush="r", symbol="o",
                                 symbolPen="r", symbolSize=5)

    def plot_wave2(self):
        self.graphicsView_4.plot(self.wave_time2,
                                 self.wave2,
                                 pen=pg.mkPen(color="r"),
                                 antialias=True, name="Resp.",
                                 symbolBrush="r", symbol="o",
                                 symbolPen="r", symbolSize=5)

    def update(self, _=0, plusnum=0):
        current_slice = int(self.horizontalSlider.value() + plusnum)

        self.graphicsView.setImage(self.array1[int(current_slice) + 1, :, :],
                                   autoRange=False,
                                   autoLevels=False,
                                   autoHistogramRange=False)

        self.graphicsView_2.setImage(self.array2[int(current_slice) + 1, :, :],
                                     autoRange=False,
                                     autoLevels=False,
                                     autoHistogramRange=False)
        
        self.line1.setValue(self.wave_time1[current_slice])
        self.line2.setValue(self.wave_time2[current_slice])

        for n, a_roi in enumerate(self.rois1):
            if current_slice in self.roi_corrections1:
                a_roi.setPos(np.array(self.roi_corrections1[current_slice + 1][n]),
                             y=None, update=False, finish=False)

        for n, a_roi in enumerate(self.rois2):
            if current_slice in self.roi_corrections1:
                a_roi.setPos(np.array(self.roi_corrections2[current_slice + 1][n]),
                             y=None, update=False, finish=False)

    def save_items(self):
        ut = time.time()
        directory_path = data_base_folder + str(self.id)
        if not os.path.isdir(directory_path):
            os.mkdir(directory_path)
        np.savez(directory_path + "/" + str(self.SOPUID[0]),
                time1=self.wave_time1,
                x1=self.wave1,
                y1=self.mwave1,
                time2=self.wave_time2,
                x2=self.wave2,
                y2=self.mwave2,
                UID=self.SOPUID,
                max_x1=self.max_x1,
                max_x2=self.max_x2,
                max_y1=self.max_y1,
                max_y2=self.max_y2,
                actime=self.actime,
                study_date=self.study_date,
                roi_center_1=self.marker_chase1,
                roi_center_2=self.marker_chase2,
                img1=self.array1,
                img2=self.array2
                )

        self.graphicsView.grab().save(directory_path + "/" + str(self.SOPUID[0]) + "AP.png")
        self.graphicsView_2.grab().save(directory_path + "/" + str(self.SOPUID[0]) + "LR.png")
        im1 = Image.open(directory_path + "/" + str(self.SOPUID[0]) + "AP.png")
        im2 = Image.open(directory_path + "/" + str(self.SOPUID[0]) + "LR.png")
        im1 = im1.crop((0,0,im1.height, im1.height))
        im2 = im2.crop((0,0,im2.height, im2.height))
        dst = Image.new('RGB', (im1.width + im2.width, im1.height))
        dst.paste(im1, (0, 0))
        dst.paste(im2, (im1.width, 0))
        dst.save(directory_path + "/" + str(self.SOPUID[0]) + "fig1.png")

    def create_report(self):
        tkroot = tkinter.Tk()
        tkroot.withdraw()
        if self.id:
            patient_save_folder = filedialog.askdirectory(initialdir=data_base_folder + str(self.id))
        else:
            patient_save_folder = filedialog.askdirectory(initialdir="./data_base/")
        tkroot.destroy()
        patient = build_report(patient_save_folder)
        print("Report created.")


class TableModel(QtCore.QAbstractTableModel):
    def __init__(self, data, headers=[]):
        super(TableModel, self).__init__()
        self._data = data
        self.headers = headers

    def data(self, index, role):
        if role == Qt.DisplayRole:
            return self._data[index.row()][index.column()]

        if role == Qt.BackgroundColorRole:
            if float(self._data[index.row()][index.column()]) > 5.0:
                if float(self._data[index.row()][index.column()]) > 10.0:
                    return QBrush(Qt.red)
                else:
                    return QBrush(Qt.yellow)

    def rowCount(self, index):
        return len(self._data)

    def columnCount(self, index):
        try:
            x = len(self._data[0])
        except ValueError:
            x = 0
        return x

    def headerData(self, col, orientation, role):
        if role == Qt.DisplayRole:
            if orientation == Qt.Horizontal:
                return QVariant(self.headers[col])
            else:
                return "marker %d" % (col + 1)
        return QVariant()


def main():
    app = QtWidgets.QApplication(sys.argv)

    window = MainWindow()
    window.show()
    app.exec()


if __name__ == "__main__":
    main()
