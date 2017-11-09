##!/usr/bin/python
# -*- coding: utf-8 -*-
################ DEVA software: D2AM Edf images Visualisation and Analysis ##############
###### Dependencies: numpy, scipy, matplotlib, lmfit (sudo easy_install -U lmfit), pyFAI, fabio, Basemap
#from gi.repository import GObject

from PyQt5.QtCore import *
from PyQt5.QtGui  import *
from PyQt5.QtWidgets import *
from PyQt5 import *
import PyQt5
from functools import partial
import xrayutilities
import threading, time
import sys
from sys import stdout
import os
import re, operator
import gc
from os import listdir
from os.path import isfile,join
import tempfile
import math
import numpy as N
from numpy import unravel_index
from scipy import ndimage, stats
from scipy.fftpack import fft, fftfreq, fftshift
from lmfit import Parameters, minimize
from PACKAGES.utilities import Combination_edf_by_translationXZ as EDF_XZ_combination
from PACKAGES.PopUpWindows import *
from PACKAGES.ReadSpecD1 import *
from PACKAGES.CommonFunctions import *
from PACKAGES.GettingData_MultiThreading import *
from PACKAGES.Geometry_Correction import *
from PACKAGES.PlotPolarGrid import *
from multiprocessing import Pool

from PACKAGES.save_h5 import App# as layout4_app
from PACKAGES.data_fits import App5# as layout4_app

from PACKAGES.new_PF import *

##############
## Graphic library ##
##############
import matplotlib as mpl
#mpl.use('GtkAgg')
from matplotlib.figure import Figure
#from matplotlib.axes import Subplot
#from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg as FigureCanvas
#from matplotlib.backends.backend_gtkagg import NavigationToolbar2GTKAgg as NavigationToolbar
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.cm import jet # colormap
from matplotlib.widgets import Cursor
from matplotlib.patches import Rectangle
from matplotlib.ticker import MaxNLocator
import colormaps as cmaps

import matplotlib.pyplot as plt
plt.register_cmap(name='magma', cmap=cmaps.magma)
plt.register_cmap(name='viridis', cmap=cmaps.viridis)
plt.register_cmap(name='inferno', cmap=cmaps.inferno)
plt.register_cmap(name='plasma', cmap=cmaps.plasma)
plt.set_cmap(cmaps.viridis)
cmap = cmaps.viridis


import fabio
import pyFAI

__author__="Tra NGUYEN THANH"
__email__ = "thanhtra0104@gmail.com"
__version__ = "3.1"
__date__="04/04/2017"

g_zoom_tb = 0
g_is_image_loaded=False

#mpl.rcParams['font.size'] = 18.0
mpl.rcParams['axes.labelsize'] = 'large'
mpl.rcParams['legend.fancybox'] = True
mpl.rcParams['legend.handletextpad'] = 0.5
mpl.rcParams['legend.fontsize'] = 'medium'
mpl.rcParams['figure.subplot.bottom'] = 0.13
mpl.rcParams['figure.subplot.top'] = 0.93
mpl.rcParams['figure.subplot.left'] = 0.14
mpl.rcParams['figure.subplot.right'] = 0.915
mpl.rcParams['image.cmap'] = 'jet'
mpl.rcParams['savefig.dpi'] = 300


#print( os.getcwd()
_SPEC_IMG_COL         = ["img", "Image_num"] #column containing image number in spec file
_IMG_NUMBER_SEPARATOR = "_"  #Please prefer to use the underscore to separate image number from prefix
_IMG_FORMAT           = ".edf"
# ********************************************************************************************************************************************
_DETECTOR_NAME    = "Maxipix"
_DETECTOR_PIXELS  = [516, 516] #Number of pixels in the first and second direction. First means vertical, second means horizontal direction
_PIXEL_SIZE        = 55e-6 #m
# ********************************************************************************************************************************************
_DETECTOR_ROTATIONS           = 1 #Number of rotations
_DETECTOR_MOTORS              = ["tth"] #Detector motors, from outermost to innermost if having more than 1 motor
_DETECTOR_ROTATION_DIRECTIONS = ['z+'] #positive rotation direction of the detector's arm, outermost --> inner most
_DETECTOR_PIXEL_DIRECTIONS    = ['z+', 'y+'] #Positive directions of the pixels, first (vertical) and second (horizontal) direction
# ********************************************************************************************************************************************
_SAMPLE_ROTATIONS             = 3 #Number of sample rotations
_SAMPLE_ROTATION_DIRECTIONS   = ['z+', 'x+', 'z+'] #Positive rotation directions of the sample motors, from outermost circle to innermost circle
_SAMPLE_ROTATION_MOTORS       = ["th", "chi", "phi"] #Sample rotation motors, from outermost to innermost circle
_XRAY_INCIDENT_BEAM_DIRECTION = 'x+' #Direction of incident X-ray beam. In general we choose this direction as x+

from os.path import expanduser
DEVA_dir = join(expanduser("~"),"DEVA")
s=open(join(DEVA_dir,"Instrument_configuration.DEVA"), "r").read()
exec(s)

    
class MyMainWindow(QMainWindow):
    def menuD5(self):
         self.set_detector(None, "D5")

    def menuD1(self):
         self.set_detector(None, "D1")

    def menuS70(self):
         self.set_detector(None, "S70")

    def menugonio(self):
        self.set_manip(None, "gonio")
    
    def menugisaxs(self):
        self.set_manip(None, "gisaxs")
    
    def menu_inst_config(self):
        self.instrument_config(None)

    def updateDetector(self):
        det_name = self.detector_lst_cb.currentText()
        self.set_detector(None, det_name)
    

    def showDetectorDialog(self):
        d = QDialog(self)
        b1 = QPushButton("Select",d)
        lbl = QLabel("Select Detector", d)
        self.detector_lst_cb = QComboBox(d)
        substrate_list = ["Maxipix", "XPAD D1", "XPAD D5", "S70", "PIXIRAD", "WOS"]
        self.detector_lst_cb.addItems(substrate_list)
        lbl.move(10, 30)
        self.detector_lst_cb.move(100, 30)
        b1.move(100, 70)
        d.setWindowTitle("Dialog")
        b1.clicked.connect( self.updateDetector)
        d.setWindowModality(Qt.ApplicationModal)
        d.exec_()

    def aboutDialog(self):
        d = QDialog(self)
        text = "Diffraction toolkit is principally used to visualize the images from ESRF or any image from synchrotron \n \
                You can view the images from different detectors provided we know the dimensions of the image. \n \
                This software is fully written in Python language. The graphical user interface is PyQt5 \n \
                and it works on Linux and Windows. \n \n \n \
                "

        lbl = QLabel(text, d)
        lbl2 = QLabel("", d)
        lbl.move(10,20)
        lbl2.move(150,600)
        d.setWindowModality(Qt.ApplicationModal)
        d.exec_()

    def __init__(self):
        super(MyMainWindow, self).__init__()
        self.setWindowTitle("Diffraction toolkit v1.0")
        #DEVA - D2AM EDF Visualisation and Analysis - version %s - Last update: %s"%(__version__, __date__))
        #self.set_size_request(1000, 850)
        #self.set_position(gtk.WIN_POS_CENTER)
        #self.set_border_width(5)
        self.setStyleSheet("QMessageBox { \
                           background-color: #222222; \
                            }")
        self.ctrlw = QWidget()
        root = os.path.dirname(os.path.realpath('__file__'))
        self.detector_type = _DETECTOR_NAME # My detector by default is XPAD D5
        self.menubar =  self.menuBar()#.addMenu(_DETECTOR_NAME)#gtk.MenuBar()

        #self.detector_menu = self.menuBar().addMenu(_DETECTOR_NAME) #gtk.Menu()
        self.detm = self.menuBar().addMenu(_DETECTOR_NAME)
        #gtk.MenuItem(_DETECTOR_NAME)
        #self.detm.set_submenu(self.detector_menu)

        self.detector_D5 = QAction(QIcon(root + '/new.png'), "&D5", self, statusTip="XPAD D5", triggered=self.menuD5)
        self.detector_D1 = QAction(QIcon(root + '/new.png'), "&D1", self, statusTip="XPAD D1", triggered=self.menuD1)
        self.detector_S70 = QAction(QIcon(root + '/new.png'), "&S70", self, statusTip="S70", triggered=self.menuS70)

        #self.inWard = self.menuBar().addMenu("&Outward Entry")
        self.detm.addAction(self.detector_D5) 
        self.detm.addAction(self.detector_D1) 
        self.detm.addAction(self.detector_S70) 

        #Menu manip
        self.experiment_type = "GONIO"
        #self.manip_menu = gtk.Menu()
        self.manipm = self.menuBar().addMenu("Experiment") #gtk.MenuItem("Experiment")
        #self.manipm.set_submenu(self.manip_menu)

        self.manip_gonio = QAction(QIcon(root + '/new.png'), "&Gonio", self, statusTip="Gonio", triggered=self.menugonio)
        self.manip_gisaxs = QAction(QIcon(root + '/new.png'), "&GISAXS", self, statusTip="GISAXS", triggered=self.menugisaxs)

        self.manipm.addAction(self.manip_gonio) 
        self.manipm.addAction(self.manip_gisaxs)

        self.instrument_menu = self.menuBar().addMenu("Instrument")
        self.instrument_item = QAction(QIcon(root + '/new.png'), "&Intrument Config", self, statusTip="Intrument Config", triggered=self.menu_inst_config)
        self.instrument_menu.addAction(self.instrument_item)

        #self.menubar.append(self.instrument_item)
        #self.menubar.append(self.detm)
        #self.menubar.append(self.manipm)

        self.toolbar = self.addToolBar("Home") #gtk.Toolbar()
        self.toolbar.setMovable(False)
        #self.toolbar.set_style(gtk.TOOLBAR_ICONS)

        self.opentb =  QAction(QIcon(root + '/stock_open_24.bmp'), "Open a folder containing EDF images", self, statusTip="Open a folder containing EDF images", triggered=self.choose_folder)
        self.refreshtb = QAction(QIcon(root + '/stock_refresh_24.bmp'), "&Reload data files", self, statusTip="Reload data files", triggered=self.folder_update)
        self.savetb =  QAction(QIcon(root + '/stock_save_24.bmp'), "Save image", self, statusTip="Save image", triggered=self.save_image)
        self.zoomtb = QAction(QIcon(root + '/stock_zoom_in_24.bmp'), "Zoom on", self, statusTip="Zoom on", triggered=self.zoom_on)
        self.hometb = QAction(QIcon(root + '/stock_home_24.bmp'), "Reset the Image", self, statusTip="Reset the Image", triggered=self.reset_image) 
        self.aspecttb = QAction(QIcon(root + '/stock_page_setup_24.png'), "Change the graph's aspect ratio", self, statusTip="Change the graph's aspect ratio", triggered=self.change_aspect_ratio) 
        self.loadcalibtb = QAction(QIcon(root + '/stock_convert_24.png'), "Load a calibration file (PONI file)", self, statusTip="Load a calibration file (PONI file)", triggered=self.load_calibration) 
        self.use_dark_tb = QAction(QIcon(root + '/stock_dialog_info_48.png'), "Use this image as dark data. Click again to cancel dark substraction.", self, statusTip="Use this image as dark data. Click again to cancel dark substraction.", triggered=self.get_dark) 
        self.quittb = QAction(QIcon(root + '/stock_exit_24.bmp'), "Quit the program", self, statusTip="Quit the program", triggered=sys.exit) 
        self.abouttb = QAction(QIcon(root + '/stock_about_24.png'), "About the program", self, statusTip="About the program", triggered=self.aboutDialog) 
        
        self.toolbar.addAction(self.opentb)#, 0)
        self.toolbar.addAction(self.refreshtb)#, 1)

        #self.toolbar.insert(self.sep)#, 2)
        self.toolbar.addAction(self.savetb)#, 3)
        self.toolbar.addAction(self.zoomtb)#, 4)
        self.toolbar.addAction(self.hometb)#, 5)
        #self.toolbar.addAction(self.aspecttb)#, 6)
        self.toolbar.addAction(self.loadcalibtb)#, 7)
        self.toolbar.addAction(self.use_dark_tb)#, 8)

        #self.toolbar.insert(self.sep2, 9)
        self.toolbar.addAction(self.quittb)#, 10)
        self.toolbar.addAction(self.abouttb)#, 10)

        #self.tooltips = gtk.Tooltips()
        #self.tooltips.set_tip(self.refreshtb,"Reload data files")
        #self.tooltips.set_tip(self.opentb,"Open a folder containing EDF images")
        self.graph_aspect = False
        self.DARK_CORRECTION = False

        ############################# BOXES ###############################################
        vbox = QVBoxLayout(self) #gtk.VBox()
        #vbox.pack_start(self.menubar,False,False,0)
        #vbox.pack_start(self.toolbar,False,False,0)
        hbox=QHBoxLayout(self) #gtk.HBox()
        ######################### TREE VIEW #############################################
        self.sw = QScrollArea()   #gtk.ScrolledWindow()
        #self.sw.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        #self.sw.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.sw.setWidgetResizable(True)
        #self.sw.setFixedHeight(50)
        self.sw.resize(450,450)
        #self.sw.set_shadow_type(gtk.SHADOW_ETCHED_IN)
        #self.sw.set_policy(gtk.POLICY_NEVER, gtk.POLICY_AUTOMATIC)

        self.store={}
        self.MODEL = QFileSystemModel()#gtk.TreeStore(str)
        self.MODEL.setFilter(QDir.AllDirs | QDir.NoDotAndDotDot | QDir.AllEntries)
        #self.MODEL.setNameFilters("*.edf")
        view_root = self.MODEL.setRootPath(os.getcwd())
        print("Model is ",self.MODEL)
        #filter1 = QStringList("")
        filter1 = ["EDF (*.edf)","EDF GZ(*.edf.gz)"] #["PONI (*.poni)", "All Files(*)"]
        self.MODEL.setNameFilters(filter1)

        #self.MODEL.setFillters()
        self.treeView = QTreeView()#gtk.TreeView(self.MODEL)
        self.treeView.resize(450, 400)
        self.treeView.setModel(self.MODEL)
        self.treeView.setSortingEnabled(True)
        self.treeView.sortByColumn(0,0)
        self.treeView.setRootIndex(view_root)
        self.treeView.hideColumn(1)
        self.treeView.hideColumn(2) 
        self.treeView.hideColumn(3) 
        self.treeView.hideColumn(4) 
        ## --------------------------------------------------------------------------------------
        ##self.treeView.connect("row-activated",self.on_changed_edf)
        #self.treeView.selectionModel().selectionChanged.connect(self.treevViewSelection)
        self.treeView.doubleClicked.connect(self.treevViewSelection)

        #rendererText = gtk.CellRendererText()
        self.TVcolumn = "EDF images" 
        #gtk.TreeViewColumn("EDF images", rendererText, text=0)
        #self.TVcolumn.set_sort_column_id(0)    
        #self.treeView.append_column(self.TVcolumn)
        self.treeView.setWindowTitle(self.TVcolumn)

        self.sw.setWidget(self.treeView)
        hbox.addWidget(self.sw, 1)#, False, False, 0)
        self.threads = []  # This is to manage multithreading jobs

        self.current_folder = os.getcwd()
        self.edf_folder     = self.current_folder
        #if self.current_folder is not os.getcwd():
        #    glib.timeout_add_seconds(5, self.folder_update)
        ############################# NOTEBOOK ################################################
        self.notebook = QTabWidget() #gtk.Notebook()
        self.page_single_figure = QHBoxLayout() #gtk.HBox()
        self.midle_panel = QVBoxLayout() #gtk.VBox()
        self.page_pole_figure = QHBoxLayout() #gtk.HBox()
        self.page_batch_correction = QHBoxLayout() #gtk.HBox() #For batch correction of EDF
        
        #**************************** Geometry setup ******************************************
        xdirection = {'x+':[1,0,0], 'x-': [-1,0,0], 'y+': [0,1,0], 'y-': [0,-1,0], 'z+':[0,0,1], 'z-': [0,0,-1]}
        self.qconv = xrayutilities.experiment.QConversion(_SAMPLE_ROTATION_DIRECTIONS,_DETECTOR_ROTATION_DIRECTIONS,xdirection[_XRAY_INCIDENT_BEAM_DIRECTION])
        #print( "qconv"
        #print( _SAMPLE_ROTATION_DIRECTIONS,_DETECTOR_ROTATION_DIRECTIONS,xdirection[_XRAY_INCIDENT_BEAM_DIRECTION]
        self.geometry_setup_tbl = QGridLayout() #gtk.Table(3,4,True)
        self.geometry_manual    = QPushButton("VALIDATE the above parameters setup")
        self.geometry_manual.clicked.connect( self.manual_calibration)
        #************** Check if the Geo_config.DEVA file exist **************
        geo_distance = ""
        geo_energy   = ""
        geo_direct   = ""
        geo_has_substrate= False
        geo_substrate= ""
        geo_substrate_other = ""
        geo_inplane  = "1 1 0"
        geo_outplane = "0 0 1"
        tmp_dir = DEVA_dir
        geo_file_name = join(tmp_dir,"Geo_config.DEVA")
        if isfile(geo_file_name):
            geoconfig_file = open(geo_file_name,'r+')
            for line in geoconfig_file:
                if line.startswith("DISTANCE"):
                    geo_distance = line.split("=")[-1].split('\n')[0]
                elif line.startswith("DIRECT_BEAM_XY"):
                    geo_direct   = line.split("=")[-1].split('\n')[0]
                elif line.startswith("ENERGY"):
                    geo_energy   = line.split("=")[-1].split('\n')[0]
                elif line.startswith("HAS_SUBSTRATE"):
                    geo_has_substrate = line.split("=")[-1].split('\n')[0]
                    if geo_has_substrate=="NO":
                        geo_has_substrate = False
                    elif geo_has_substrate=="YES":
                        geo_has_substrate = True
                elif line.startswith("SUBSTRATE"):
                    geo_substrate = line.split("=")[-1].split('\n')[0]
                elif line.startswith("INPLANE"):
                    geo_inplane = line.split("=")[-1].split('\n')[0]
                elif line.startswith("OUTPLANE"):
                    geo_outplane = line.split("=")[-1].split('\n')[0]
                
            geoconfig_file.close()
            
        self.geometry_distance_txt = QLabel("Distance Samp-Det (m):")
        self.geometry_distance     = QLineEdit()
        #self.geometry_distance.set_usize(30,0)
        self.geometry_distance.setText(geo_distance)
        self.geometry_direct_beam_txt = QLabel("Direct beam X,Y:")
        self.geometry_direct_beam     = QLineEdit()
        #self.geometry_direct_beam.set_usize(30,0)
        self.geometry_direct_beam.setText(geo_direct)
        self.geometry_energy_txt      = QLabel("Energy (eV):")
        self.geometry_energy          = QLineEdit()
        #self.geometry_energy.set_usize(30,0)
        self.geometry_energy.setText(geo_energy)
        #self.geometry_distance_txt.set_alignment(0,0.5)
        #self.geometry_direct_beam_txt.set_alignment(0,0.5)
        #self.geometry_energy_txt.set_alignment(0,0.5)
        
        self.geometry_UB_txt  = QLabel("Import a PONI file:")
        self.geometry_browse_UB = QPushButton("Browse PONI file")
        self.geometry_browse_UB.setCheckable(True)
        #self.geometry_browse_UB.set_usize(30,0)
        self.geometry_browse_UB.clicked[bool].connect(self.load_UBfile)
        self.UB_MATRIX_LOAD = False #By default, no UB matrix file is loaded
        
        self.geometry_substrate_txt  = QLabel("Substrate material:")
        self.geometry_substrate_other_txt  = QLabel("If other:")
        self.geometry_substrate_inplane_txt= QLabel("In-plane direction:")
        self.geometry_substrate_outplane_txt= QLabel("Normal direction:")
        
        self.geometry_direct_beam_txt.setToolTip("Position (in pixel) of the direct beam when all motors are at zero. X and Y position are separated by a comma, e.g. 300.5,650.7")
        self.geometry_substrate_txt.setToolTip("Substrate material")
        self.geometry_substrate_other_txt.setToolTip("The substrate material, i.e. Al, SiO2, CdTe, GaN,...")
        self.geometry_substrate_inplane_txt.setToolTip( "The substrate in-plane direction - separated by space - for calculation of the orientation matrix.")
        self.geometry_substrate_outplane_txt.setToolTip( "The substrate out-of-plane direction - separated by space - for calculation of the orientation matrix.")
        
        #self.geometry_substrate_txt.set_alignment(0,0.5)
        #self.geometry_substrate_other_txt.set_alignment(0,0.5)
        #self.geometry_substrate_inplane_txt.set_alignment(0,0.5)
        #self.geometry_substrate_outplane_txt.set_alignment(0,0.5)
        self.geometry_substrate = QComboBox()
        substrate_list = ["-- other", "Si", "Ge", "GaAs", "GaN", "GaP", "GaSb", "InAs", "InP", "InSb", "Al2O3"]
        self.geometry_substrate.addItems(substrate_list) #set_usize(30,0)
        #for s in range(len(substrate_list)):
        #    self.geometry_substrate.append_text(substrate_list[s])
        #self.geometry_substrate.setChecked(0)
        self.has_substrate = False
        
        if geo_has_substrate:
            if geo_substrate in substrate_list:
                geo_substrate_other = ""
                for i in range(len(substrate_list)):
                    if geo_substrate == substrate_list[i]:
                        self.geometry_substrate.setCurrentIndex(i) #setChecked(i)
                        break
            else:
                geo_substrate_other = geo_substrate
        
        self.geometry_substrate_other = QLineEdit()
        #self.geometry_substrate_other.set_usize(30,0)
        self.geometry_substrate_other.setText(geo_substrate_other)
        self.geometry_substrate_inplane = QLineEdit()
        #self.geometry_substrate_inplane.set_usize(30,0)
        self.geometry_substrate_inplane.setText(geo_inplane)
        self.geometry_substrate_outplane = QLineEdit()
        #self.geometry_substrate_outplane.set_usize(30,0)
        self.geometry_substrate_outplane.setText(geo_outplane)
        
        self.geometry_UB_txt.setToolTip("Import a UB matrix which is a text file with a 3x3 matrix (3 lines, 3 colunms)")
        
        #self.geometry_UB_txt.set_alignment(0,0.5)
                
        self.geometry_setup_tbl.addWidget(self.geometry_UB_txt, 1,1)
        self.geometry_setup_tbl.addWidget(self.geometry_browse_UB, 1,2)
        self.geometry_setup_tbl.addWidget(self.geometry_energy_txt, 2,1)
        self.geometry_setup_tbl.addWidget(self.geometry_energy, 2,2)
        self.geometry_setup_tbl.addWidget(self.geometry_distance_txt, 3,1)
        self.geometry_setup_tbl.addWidget(self.geometry_distance, 3,2)
        self.geometry_setup_tbl.addWidget(self.geometry_direct_beam_txt, 4,1)
        self.geometry_setup_tbl.addWidget(self.geometry_direct_beam, 4,2)
        
        self.geometry_setup_tbl.addWidget(self.geometry_substrate_txt, 1,3)
        self.geometry_setup_tbl.addWidget(self.geometry_substrate,1,4)
        self.geometry_setup_tbl.addWidget(self.geometry_substrate_other_txt,2,3)
        self.geometry_setup_tbl.addWidget(self.geometry_substrate_other,2,4)
        self.geometry_setup_tbl.addWidget(self.geometry_substrate_inplane_txt,3,3)
        self.geometry_setup_tbl.addWidget(self.geometry_substrate_inplane,3,4)
        self.geometry_setup_tbl.addWidget(self.geometry_substrate_outplane_txt,4,3)
        self.geometry_setup_tbl.addWidget(self.geometry_substrate_outplane,4,4)
        
        ############################# PAGE1: FIGURES ##########################################
        self.edf = ""
        self.edf_choosen = ""
        self.my_notes = []
        self.lines = []
        self.points=[]
        self.arb_lines_X=[]
        self.arb_lines_Y=[]
        self.arb_line_points=0
        self.SELECTED_IMG_NUM = None
        self.IMG_INIT = False

        self.header = {}
        self.calibrated = False #The detector need to be calibrated using pyFAI program
        self.calibrated_quantitative = False
        self.fig=Figure(dpi=100)
        #self.ax  = self.fig.add_subplot(111)
        self.ax  = self.fig.add_axes([0.12,0.2,0.7,0.7])
        self.MAIN_XLABEL = self.ax.set_xlabel("X (pixel)")
        self.MAIN_YLABEL = self.ax.set_ylabel("Y (pixel)")
        #self.fig.text(0.5, 0.92, self.edf, ha='center', fontsize=26)
        self.xDim0 = 0
        self.xDim1 = _DETECTOR_PIXELS[1]
        self.yDim0 = 0
        self.yDim1 = _DETECTOR_PIXELS[0]
        self.MAIN_EXTENT = (self.xDim0,self.xDim1,self.yDim0,self.yDim1)
        self.fig.subplots_adjust(left=0.1,bottom=0.20, top=0.90)
        self.data = N.zeros(shape=(self.yDim1,self.xDim1))
        self.vmin = 0
        self.vmax = 1000
        self.vmax_range = self.vmax

        #self.init_image()
        self.img = self.ax.imshow(self.data,vmin=self.vmin, vmax=self.vmax, cmap=cmap, interpolation='nearest',aspect='auto')

        self.canvas  = FigureCanvas(self.fig)
        self.main_figure_navBar = NavigationToolbar(self.canvas, self)
        self.cursor = Cursor(self.ax, color='k', linewidth=1, useblit=True)
        #Global color bar
        self.cax = self.fig.add_axes([0.85, 0.20, 0.03, 0.70])#left,bottom,width,height
        self.cb  = self.fig.colorbar(self.img, cax = self.cax, format="%.f")#format=fm
        self.cb.set_label(r'$Intensity\ (Counts\ per\ second)$', fontsize=12)
        self.cb.locator = MaxNLocator(nbins=6)
        #Horizontal colorbar
        self.cax2 = self.fig.add_axes([0.12, 0.1, 0.7, 0.03])#left,bottom,width,height
        self.cb2  = self.fig.colorbar(self.img, cax = self.cax2, orientation='horizontal', format="%.f")#format=fm
        self.cb2.set_label(r'$Intensity\ (Counts\ per\ second)$', fontsize=12)
        self.cb2.locator = MaxNLocator(nbins=6)
        self.cb2.ax.set_visible(False)
        ########### ZOOM ACTION #################################
        self.zoom = False #Press the Zoom button on the tool bar
        self.zoom_press = False #Button pressed in the axe to draw a rectangle
        self.ROI_ON = False
        self.ROI_press = False
        self.rect = Rectangle((0,0),0,0)
        self.roi_rect = Rectangle((0,0),0,0)
        self.x0 =0
        self.y0 =0
        self.x1 =0
        self.y1 =0
        self.ROI_x0 =0
        self.ROI_y0 =0
        self.ROI_x1 =0
        self.ROI_y1 =0
        self.ax.add_patch(self.rect)
        self.ax.add_patch(self.roi_rect)
        self.canvas.mpl_connect("motion_notify_event",self.on_motion)
        self.canvas.mpl_connect("button_press_event",self.on_press)
        self.canvas.mpl_connect("button_release_event",self.on_release)
        self.mouse_moved = False #If click without move: donot zoom the image

        #self.midle_panel.addWidget(self.main_figure_navBar)#, False,False, 2)
        #*********************************** SCAN SLIDER *****************************************
        
        #Variables for scan slider:
        self.SCAN_IMG = []
        self.SPEC_IMG = []#List of all spec images
        self.SPEC_FILE = ""
        self.SELECTED_IMG_NUM = None
        self.SPEC_IS_LOADED = False
        self.SPEC_DATA = None #SPEC object from xrayutilities.io.SPECFile(specfile)
        self.SPEC_SCAN_LIST = None #list of all spec scan, each spec scan is an object from spec_data.scan_list 
        self.SPEC_ACTUAL_SCAN = None #actual scan number
        self.SPEC_ACTUAL_SCAN_IMG = 0
        self.DATA_IS_LOADED = False
        self.SPEC_ACTUAL_SCAN_DATA = []
        self.SPEC_ALL_MOTORS_LIST  = []#List of all motor_scan, in uppercase (PHI,TSZ,ETA,...)
        self.SPEC_SKIPPED_MOTORS   = []#List of the scanning motors that we don't want to search
        
        self.scan_slider_frame = QFrame()#gtk.Frame()
        self.scan_slider_table_align = QGridLayout() #gtk.Alignment(0,0.5,1,1)
        #self.scan_slider_table_align.set_padding(10,5,5,5)
        skip_box = QHBoxLayout()
        
        #self.scan_slider_frame.setTile("Scan Slider")
        #self.scan_slider_frame.set_label_align(0.5,0.5)
        self.scan_slider_table = QGridLayout() #gtk.Table(3,3,False)
        self.scan_slider_specfile_txt = QLabel("Spec file")
        #self.scan_slider_specfile_txt.set_alignment(0,0.5)
        self.scan_slider_scanNumber_txt = QLabel("Scan #")
        #self.scan_slider_scanNumber_txt.set_alignment(0,0.5)
        self.scan_slider_path = QLineEdit()
        self.scan_slider_browseSpec = QPushButton("Browse spec file")
        #self.scan_slider_browseSpec.set_size_request(80,-1)
        self.scan_slider_browseSpec.clicked.connect(self.load_specFile)
        #scan_slider_spin_adj = gtk.Adjustment(1,0,1,1,10,0)#actual, min, max, step increment, page increment, page size
        
        self.scan_slider_spinButton = QSpinBox()
        #gtk.SpinButton(scan_slider_spin_adj, 1,0)
        self.scan_slider_imgSlider  = QSlider(Qt.Horizontal) #gtk.HScale()
        self.scan_slider_imgSlider.setRange(0,1)
        self.scan_slider_imgSlider.setValue(1)
        #self.scan_slider_imgSlider.set_digits(0)
        #self.scan_slider_imgSlider.set_increments(1,1)
        #self.scan_slider_spinButton.set_numeric(True)
        self.scan_slider_spinButton.valueChanged.connect(self.update_scan_slider)
        self.scan_slider_imgSlider.valueChanged[int].connect(self.slider_plot_scan)
        
        #Skip scans:
        self.scan_slider_skip_scans = QLabel("Skipped scans:")
        #self.scan_slider_skip_scans.set_alignment(0,0.5)
        self.scan_slider_skip_tsz   = QCheckBox("Tsz")
        self.scan_slider_skip_eta   = QCheckBox("Eta")
        self.scan_slider_skip_del   = QCheckBox("Del")
        self.scan_slider_skip_chi   = QCheckBox("Chi")
        self.scan_slider_skip_phi   = QCheckBox("Phi")
        self.scan_slider_skip_rox   = QCheckBox("Rox")
        self.scan_slider_skip_roy   = QCheckBox("Roy")
        self.scan_slider_skip_tox   = QCheckBox("Tox")
        self.scan_slider_skip_toy   = QCheckBox("Toy")
        self.scan_slider_skip_rien  = QCheckBox("Rien")
        #skip_box.pack_start(self.scan_slider_skip_scans, False, False, 0)
        #Plot roi
        plot_roi_txt = QLabel("Plot ROI:")
        plot_roi_txt.setToolTip( "Check this box to draw a ROI. Click on the image and drag the mouse to draw. Click again to ecrase the ROI")
        #plot_roi_txt.set_alignment(0,0.5)
        self.draw_roi_btn = QCheckBox("Draw a ROI")
        self.draw_roi_btn.stateChanged.connect(self.Enable_draw_roi)
        self.plot_3D_scan_btn = QPushButton("Plot 3D HKL of this scan")
        self.plot_3D_scan_btn.clicked.connect(self.plot_3D_scan)
        roi_box = QHBoxLayout() #gtk.HBox()
        roi_box.addWidget(self.draw_roi_btn)#, False, False, 0)
        self.use_spec_UB = QCheckBox("Use UB from Spec (if applicable)")
        roi_box.addWidget(self.use_spec_UB)#, False, False, 0)
        
        skip_box.addWidget(self.scan_slider_skip_tsz)#, False, False, 0)
        skip_box.addWidget(self.scan_slider_skip_eta)#, False, False, 0)
        skip_box.addWidget(self.scan_slider_skip_del)#, False, False, 0)
        skip_box.addWidget(self.scan_slider_skip_chi)#, False, False, 0)
        skip_box.addWidget(self.scan_slider_skip_phi)#, False, False, 0)
        skip_box.addWidget(self.scan_slider_skip_rox)#, False, False, 0)
        skip_box.addWidget(self.scan_slider_skip_roy)#, False, False, 0)
        skip_box.addWidget(self.scan_slider_skip_tox)#, False, False, 0)
        skip_box.addWidget(self.scan_slider_skip_toy)#, False, False, 0)
        skip_box.addWidget(self.scan_slider_skip_rien)#, False, False, 0)
        
        self.scan_slider_table.addWidget(self.scan_slider_specfile_txt, 1,1)
        self.scan_slider_table.addWidget(self.scan_slider_path, 1,2)
        self.scan_slider_table.addWidget(self.scan_slider_browseSpec, 1,3)
        self.scan_slider_table.addWidget(self.scan_slider_scanNumber_txt, 2,1)
        self.scan_slider_table.addWidget(self.scan_slider_spinButton, 2,2)
        self.scan_slider_table.addWidget(self.scan_slider_imgSlider, 2,3)
        #self.scan_slider_table.attach(self.scan_slider_skip_scans,0,1,2,3)
        #self.scan_slider_table.attach(skip_box, 1,3,2,3)
        self.scan_slider_table.addWidget(plot_roi_txt, 3,1)
        self.scan_slider_table.addLayout(roi_box, 3,2)
        self.scan_slider_table.addWidget(self.plot_3D_scan_btn, 3,3)

        self.detector_disposition_horizontal = QPushButton("Rotate detector")
        self.linear_scale_btn = QPushButton("Log scale")
        self.save_adj_btn = QPushButton("Save Corrected EDF")

        self.scan_slider_table.addWidget(self.linear_scale_btn, 4,1)
        self.scan_slider_table.addWidget(self.detector_disposition_horizontal, 4,2)
        self.scan_slider_table.addWidget(self.save_adj_btn, 4,3)
        
        self.scan_slider_skip_tsz.setChecked(True)
        self.scan_slider_skip_rien.setChecked(True)
        
        self.scan_slider_table_align.addLayout(self.scan_slider_table, 1, 1)
        self.scan_slider_frame.setLayout(self.scan_slider_table_align)

        #self.midle_panel_grid = QGridLayout() 

        self.midle_panel.addLayout(self.geometry_setup_tbl)
        self.midle_panel.addWidget(self.geometry_manual)
        self.midle_panel.addWidget(self.main_figure_navBar)
        self.midle_panel.addWidget(self.canvas)
        self.midle_panel.addWidget(self.scan_slider_frame)
        self.page_single_figure.addLayout(self.midle_panel)

        ########################################## Check Buttons RIGHT PANEL ###################

        self.right_panel = QGridLayout()#VBoxLayout() #gtk.VBox(False,0)

        # self.detector_disposition_horizontal = gtk.ToggleButton("Rotate detector")
        #self.detector_disposition_horizontal = QPushButton("Rotate detector")
        self.rotate_detector_n = 0
        # self.detector_disposition_horizontal.set_sensitive(False)
        self.detector_disposition_horizontal.clicked.connect(self.detector_disposition)
        self.horizontal_detector = False #By default, the detector is in the vertical position, i.e. 960 rows x 560 cols

        #self.linear_scale_btn = QPushButton("Log scale")
        self.linear_scale_btn.setCheckable(True)
        self.linear_scale_btn.clicked[bool].connect(self.log_update)

        self.log_scale=0

        self.adj_btn = QCheckBox("Geometry correction")
        self.adj_btn.stateChanged.connect(self.plot_update)
        
        self.detector_space_btn = QRadioButton("Detector map")
        self.detector_space_btn.setChecked(True)
        self.detector_space_btn.toggled.connect(self.plot_update)
        
        self.tth_chi_space_btn = QRadioButton("2Theta-Chi map")
        #gtk.RadioButton(self.detector_space_btn, "2Theta-Chi map")
        self.tth_chi_space_btn.toggled.connect(self.plot_update)
        
        self.hk_space_btn = QRadioButton("HK map")
        #gtk.RadioButton(self.detector_space_btn, "HK map")
        self.hk_space_btn.toggled.connect( self.plot_update)
        
        self.hl_space_btn = QRadioButton("HL map")
        #gtk.RadioButton(self.detector_space_btn, "HL map")
        self.hl_space_btn.toggled.connect(self.plot_update)
        
        self.kl_space_btn = QRadioButton("KL map")
        #gtk.RadioButton(self.detector_space_btn, "KL map")
        self.kl_space_btn.toggled.connect(self.plot_update)
        
        self.q_space_btn = QRadioButton("Q map")
        #gtk.RadioButton(self.detector_space_btn, "Q map")
        self.q_space_btn.toggled.connect(self.plot_update)

        #self.save_adj_btn = QPushButton("Save Corrected EDF")
        self.save_adj_btn.clicked.connect(self.save_adjust)

        self.separator = QSplitter(Qt.Horizontal)#gtk.HSeparator()
        #self.plotXYprofiles_btn = gtk.CheckButton("Plot X,Y profiles") #Plot a cross profile of X and Y data
        self.plotXYprofiles_btn = QRadioButton("X,Y profiles")
        self.plotXYprofiles_btn.setChecked(True)
        self.arbitrary_profiles_btn = QRadioButton("Arbitrary profiles")

        self.show_chi_delta_btn = QCheckBox("Show 2Theta,Chi")
        self.show_chi_delta_btn.stateChanged.connect(self.show_chi_delta)
        #self.show_chi_delta_btn.toggle()
        self.show_chi_delta_flag = False
        self.show_chi_txt   = QLabel()
        #self.show_chi_txt.set_alignment(0,0)
        self.show_delta_txt = QLabel()
        #self.show_delta_txt.set_alignment(0,0)
        self.show_d_txt = QLabel()
        #self.show_d_txt.set_alignment(0,0)
        #### Pack these options in a table
        self.option_table = QGridLayout()
        #gtk.Table(5,3,False) #5 rows, 3 cols, homogeneous
        
        #self.option_table.addWidget(self.linear_scale_btn, 1,1)
        #self.option_table.addWidget(self.detector_disposition_horizontal, 1,2)
        #self.option_table.addWidget(self.save_adj_btn, 1,3)
        
        #self.option_table.addWidget(self.detector_space_btn, 2,1)
        self.option_table.addWidget(self.adj_btn, 3,1)
        self.option_table.addWidget(self.plotXYprofiles_btn,4,1)
        self.option_table.addWidget(self.arbitrary_profiles_btn,5,1)
        
        self.option_table.addWidget(self.tth_chi_space_btn,3,2)
        # self.option_table.addWidget(self.hk_space_btn, 1,2,2,3)
        self.option_table.addWidget(self.q_space_btn, 4,2)
        self.option_table.addWidget(self.kl_space_btn, 5,2)
        # self.option_table.addWidget(self.kl_space_btn, 1,2,4,5)
        
        self.option_table.addWidget(self.show_chi_delta_btn,3,3)
        self.option_table.addWidget(self.show_delta_txt,2,3,2,3)
        self.option_table.addWidget(self.show_chi_txt, 2,3,3,4)
        self.option_table.addWidget(self.show_d_txt, 2,3,4,5)


        ### Options for profile plots
        self.profiles_log_btn = QPushButton("Y-Log")
        self.profiles_log_btn.setCheckable(True)
        self.profiles_log_btn.clicked[bool].connect(self.profiles_update)
        self.profiles_export_data_btn = QPushButton("Export data")
        self.profiles_export_data_btn.clicked.connect(self.profiles_export)
        integration_width = QLabel(" Profile integration width: ")
        
        self.integration_width = QLineEdit()
        self.integration_width.setText("10")
        #self.integration_width.set_usize(40,0)
        
        self.integration_width.setToolTip("Integration width in pixel. This is not applied to arbitrary profile")

        self.profiles_option_box = QHBoxLayout() #gtk.HBox(False,0)
        self.profiles_option_box.addWidget(self.profiles_log_btn)
        self.profiles_option_box.addWidget(self.profiles_export_data_btn)
        #self.profiles_option_box.addWidget(integration_width)
        #self.profiles_option_box.addWidget(self.integration_width)
        ### Figure of profiles plot
        self.fig_profiles = Figure()
        self.profiles_fringes = []
        self.profiles_ax2 = self.fig_profiles.add_subplot(211)
        self.profiles_ax2.set_xlabel("X profile", size=12)
        self.profiles_ax1 = self.fig_profiles.add_subplot(212)
        self.profiles_ax1.set_xlabel("Y profile", size=12)
        
        self.fig_profiles.subplots_adjust(bottom=0.1, top=0.95, hspace=0.30)
        self.profiles_canvas = FigureCanvas(self.fig_profiles)
        #self.profiles_canvas.set_size_request(450,50)
        self.profiles_canvas.resize(450,50)
        self.profiles_canvas.mpl_connect("button_press_event",self.profile_press)
        self.profiles_navBar = NavigationToolbar(self.profiles_canvas, self)
        self.cursor_pro1 = Cursor(self.profiles_ax1, color='k', linewidth=1, useblit=True)
        self.cursor_pro2 = Cursor(self.profiles_ax2, color='k', linewidth=1, useblit=True)


        #### Results of fitted curves
        self.fit_results_table = QGridLayout() # gtk.Table(7,3, False)
        title = QLabel()#"Fitted results:")
        self.chi_title = QLabel("Y")
        self.tth_title = QLabel("X")
        y0 = QLabel("y0:")
        xc = QLabel("xc:")
        A = QLabel("A:")
        w = QLabel("FWHM:")
        mu = QLabel("mu:")
        #y0.set_alignment(0,0.5)
        #xc.set_alignment(0,0.5)
        #A.set_alignment(0,0.5)
        #w.set_alignment(0,0.5)
        #mu.set_alignment(0,0.5)

        self.chi_fitted_y0 = QLabel()
        self.chi_fitted_xc = QLabel()
        self.chi_fitted_A = QLabel()
        self.chi_fitted_w = QLabel()
        self.chi_fitted_mu = QLabel()

        self.tth_fitted_y0 = QLabel()
        self.tth_fitted_xc = QLabel()
        self.tth_fitted_A = QLabel()
        self.tth_fitted_w = QLabel()
        self.tth_fitted_mu = QLabel()

        self.fit_results_table.addWidget(title,1,1,3,1)
        #self.fit_results_table.addWidget(self.chi_title,2,3,1,1)
        #self.fit_results_table.addWidget(self.tth_title,2,2,1,2)
        #self.fit_results_table.addWidget(y0,3,1)
        #self.fit_results_table.addWidget(xc,4,1)
        #self.fit_results_table.addWidget(A,5,1)
        #self.fit_results_table.addWidget(w,6,1)
        #self.fit_results_table.addWidget(mu,7,1)

        #self.fit_results_table.addWidget(self.chi_fitted_y0,3,3)
        #self.fit_results_table.addWidget(self.chi_fitted_xc,4,3)
        #self.fit_results_table.addWidget(self.chi_fitted_A,5,3)
        #self.fit_results_table.addWidget(self.chi_fitted_w,6,3)
        #self.fit_results_table.addWidget(self.chi_fitted_mu,7,3)

        #self.fit_results_table.addWidget(self.tth_fitted_y0,3,2)
        #self.fit_results_table.addWidget(self.tth_fitted_xc,4,2)
        #self.fit_results_table.addWidget(self.tth_fitted_A,5,2)
        #self.fit_results_table.addWidget(self.tth_fitted_w,6,2)
        #self.fit_results_table.addWidget(self.tth_fitted_mu,7,2)
        #self.fit_results_table.hide()
        #### PACK the right panel
        self.right_panel.addLayout(self.option_table,1,1,6,4)
        self.right_panel.addLayout(self.profiles_option_box,7,1,1,4)
        self.right_panel.addWidget(self.profiles_navBar,8,1,2,4)
        self.right_panel.addWidget(self.profiles_canvas,10,1,6,4)
        self.right_panel.addLayout(self.fit_results_table,16,1,7,4)

        #hbox.pack_end(fixed,False,False,5)
        self.page_single_figure.addLayout(self.right_panel)
        self.page_single_figure_w = QWidget()
        self.page_single_figure_w.setLayout( self.page_single_figure )
        self.page_pole_figure_w = QWidget()
        self.page_pole_figure_w.setLayout( self.page_pole_figure )
        self.page_batch_correction_w = QWidget()        
        self.page_batch_correction_w.setLayout( self.page_batch_correction )
        self.l4_app = App() #layout4_app()
        self.l5_app = App5() #layout4_app()
        self.notebook.addTab(self.page_single_figure_w, "EDF viewer")
        self.notebook.addTab(self.l4_app, "Convert edf to h5")
        self.notebook.addTab(self.page_pole_figure_w, "Pole figure")
        #self.notebook.addTab(self.page_batch_correction_w, "EDF batch correction")
        self.notebook.addTab(self.l5_app, "Temperature analysis")

        ###################### PAGE 2: POLE FIGURES ###################################
        self.left_panel_polefigure = QVBoxLayout() 
        #gtk.VBox() #This box is to contain options for pole figure construction
        self.right_panel_polefigure = QVBoxLayout()
        #gtk.VBox()

        self.images_used = QLabel("") #"Use images")
        #self.images_used.set_alignment(0,0.5)
        self.images_from = QLabel("From #: ")
        self.images_from_nb = QLineEdit()

        self.images_to   = QLabel("To #: ")
        self.images_to_nb = QLineEdit()

        self.h5SrcLbl = QLabel("Sample")
        self.h5SrcTb = QLineEdit()
        self.h5SrcFldBtn = QPushButton("Browse")
        self.h5SrcFldBtn.clicked.connect( self.srcFileForPlot )
        self.logLbl = QLabel("")
        self.logTb = QCheckBox("Log")
        self.wavelengthLbl = QLabel("Wavelength")
        self.wavelengthTb = QLineEdit()
        self.startLbl = QLabel("Start")
        self.startTb = QLineEdit()
        self.startTb.setValidator(QIntValidator())
        self.endLbl = QLabel("End")
        self.endTb = QLineEdit()
        self.endTb.setValidator(QIntValidator())
        self.omegaLbl = QLabel("Omega_25")
        self.omegaTb = QLineEdit()
        self.stepLbl = QLabel("Step")
        self.setpTb = QLineEdit()
        #self.setpTb.setValidator(QIntValidator())
        self.destFldLbl = QLabel("Destination")
        self.destFldTb = QLineEdit()
        self.destFldBtn = QPushButton("Browse")
        self.destFldBtn.clicked.connect( self.destFldForPlot )

        self.pole_2theta_txt  = QLabel("2 Theta: ")
        self.pole_2theta_field = QLineEdit()
        
        self.pole_chi_min_txt  = QLabel("Chi minimum: ")
        self.PF_chi_min = QLineEdit()

        self.select_phi = QRadioButton("Plot PHI")
        self.select_kphi= QRadioButton("Plot KPHI")#(self.select_phi,"PlotKPHI")

        self.check_chi_positive = QCheckBox("Chi positive only")

        self.plot_pole_figure_btn = QPushButton("Plot")
        self.plot_pole_figure_btn.clicked[bool].connect(self.generatePlotFigure)
        self.plot_pole_all_figure_btn = QPushButton("Plot All")
        self.plot_pole_all_figure_btn.clicked[bool].connect(self.generatePlotAlllFigure)
        self.log_pole_figure_btn = QPushButton("Log scale")
        self.log_pole_figure_btn.setCheckable(True)
        self.log_pole_figure_btn.clicked[bool].connect(self.replot_PF)
        #self.log_pole_figure_btn.clicked[bool].connect(self.generatePlotFigure)
        self.data_loading = QProgressBar() #gtk.ProgressBar()

        self.input_table = QGridLayout() #gtk.Table(5,2,True) #rows, cols, homogeneous
        self.input_table.addWidget(self.h5SrcLbl, 1,1)
        self.input_table.addWidget(self.h5SrcTb, 1,2)
        #self.input_table.addWidget(self.logLbl, 2,1)
        self.input_table.addWidget(self.h5SrcFldBtn, 2,1)
        self.input_table.addWidget(self.logTb, 2,2)
        self.input_table.addWidget(self.wavelengthLbl, 3,1)
        self.input_table.addWidget(self.wavelengthTb, 3,2)
        self.input_table.addWidget(self.startLbl, 4,1)
        self.input_table.addWidget(self.startTb, 4,2)
        self.input_table.addWidget(self.endLbl, 5,1)
        self.input_table.addWidget(self.endTb, 5,2)
        self.input_table.addWidget(self.omegaLbl, 6,1)
        self.input_table.addWidget(self.omegaTb, 6,2)
        self.input_table.addWidget(self.stepLbl, 7,1)
        self.input_table.addWidget(self.setpTb, 7,2)
        self.input_table.addWidget(self.destFldLbl, 8,1)
        self.input_table.addWidget(self.destFldTb, 8,2)
        self.input_table.addWidget(self.destFldBtn, 9,2)
        self.input_table.addWidget(self.pole_2theta_txt, 10,1)
        self.input_table.addWidget(self.pole_2theta_field, 10,2)
        #self.input_table.addWidget(self.images_from, 1,1,1,2)
        #self.input_table.addWidget(self.images_from_nb, 1,3,1,2)
        #self.input_table.addWidget(self.images_to, 2,1,1,2)
        #self.input_table.addWidget(self.images_to_nb,2,3,1,2)
        #self.input_table.addWidget(self.select_phi, 4,1,1,2)
        #self.input_table.addWidget(self.select_kphi, 4,3,1,2)

        self.left_panel_polefigure.addWidget(self.images_used)#, False, False, 5)
        self.left_panel_polefigure.addLayout(self.input_table)#, False, False, 2)
        self.left_panel_polefigure.addWidget(self.plot_pole_figure_btn)#, False, False, 2)
        self.left_panel_polefigure.addWidget(self.plot_pole_all_figure_btn)
        self.left_panel_polefigure.addWidget(self.data_loading)#, False, False, 2)
        self.page_pole_figure.addLayout(self.left_panel_polefigure)#,False,False,10)

        self.polefig=Figure(figsize=(10,8),dpi=100)
        self.polar_ax  = self.polefig.add_axes([0.1,0.1,0.75,0.8],polar=True)
        #self.polar_ax  = self.polefig.add_axes([0.1,0.1,0.75,0.8])
        #self.polar_ax.set_frame_on(False)
        self.polar_cax = self.polefig.add_axes([0.85,0.1,0.03,0.8])
        self.pole_canvas  = FigureCanvas(self.polefig)
        #self.PF_cursor = Cursor(self.polar_ax, color='k', linewidth=1, useblit=True)
        self.PF_navBar = NavigationToolbar(self.pole_canvas, self)

        self.right_panel_polefigure.addWidget(self.PF_navBar)#,False,False,0)
        self.right_panel_polefigure.addWidget(self.pole_canvas)#,True,True,0)
        self.page_pole_figure.addLayout(self.right_panel_polefigure)#,True, True,10)
        
        #*****************************************************************************
        #                  PAGE 3: CORRECTION OF A BATCH EDF
        #*****************************************************************************
        self.table_1 = QGridLayout() #gtk.Table(3,3, False)
        self.table_2 = QGridLayout() #gtk.Table(7,5, False)
        self.t1_src_folder_txt= QLabel("Source folder:")
        #self.t1_src_folder_txt.set_alignment(0,0.5)
        self.t1_src_path  = QLineEdit()
        #self.t1_src_path.set_usize(100,0)
        
        self.t1_src_button= QPushButton("Browse")
        self.t1_src_button.clicked.connect(self.select_source_folder)
        
        self.t1_des_folder_txt = QLabel("Destination folder:")
        #self.t1_des_folder_txt.set_alignment(0,0.5)
        self.t1_des_path = QLineEdit()
        #self.t1_des_path.set_usize(100,0)
        self.t1_des_button = QPushButton("Browse")
        self.t1_des_button.clicked.connect(self.select_destination_folder)
        
        self.t1_dark_img_txt = QLabel("Dark image:")
        #self.t1_dark_img_txt.set_alignment(0,0.5)
        self.t1_dark_img_path = QLineEdit()
        #self.t1_dark_img_path.set_usize(100,0)
        self.t1_dark_img_button = QPushButton("Browse")
        self.t1_dark_img_button.clicked.connect(self.select_dark_image)
        self.batch_DARK_CORRECTION = False
        
        self.t2_detector_txt = QLabel("Detector: ")
        #self.t2_detector_txt.set_alignment(0,0.5)
        self.t2_det_combobox = QComboBox()   #gtk.combo_box_new_text()
        self.t2_det_combobox.addItems(["D1","D5", "S70"])
        #self.t2_det_combobox.append_text("D1")
        #self.t2_det_combobox.append_text("D5")
        #self.t2_det_combobox.append_text("S70")
        #self.t2_det_combobox.setChecked(1)
        self.t2_det_combobox.currentIndexChanged.connect(self.batch_change_detector)
        
        self.d1_specfile_txt = QLabel("Spec file:")
        #self.d1_specfile_txt.set_alignment(1,0.5)
        self.d1_specfile_path = QLineEdit()
        self.d1_specfile_browse = QPushButton("Browse")
        #-----------------------------------------------------------------------------
        self.d1_specfile_browse.clicked.connect(self.select_normalisation_file)#,self.d1_specfile_path)
        
        self.d1_center_txt = QLabel("Center, separated by comma (X,Y):")
        #self.d1_center_txt.set_alignment(1,0.5)
        self.d1_center = QLineEdit()
        
        self.d1_distance_txt = QLabel("Distance (mm):")
        #self.d1_distance_txt.set_alignment(1,0.5)
        self.d1_distance = QLineEdit()
        
        self.t2_monitor_txt = QLabel("Monitor counter: ")
        #self.t2_monitor_txt.set_alignment(0,0.5)
        self.t2_mon_combobox = QComboBox() #gtk.combo_box_new_text()
        self.t2_mon_combobox.addItems(["vct1", "vtc2","vct3","d0_cps","pm01","pm02","pmo","pm1","pm2","roi1","roi2","roi3","roi4"])
        #self.t2_mon_combobox.append_text("vct1")
        #self.t2_mon_combobox.append_text("vct2")
        #self.t2_mon_combobox.append_text("vct3")
        #self.t2_mon_combobox.append_text("d0_cps")
        #self.t2_mon_combobox.append_text("pm01")
        #self.t2_mon_combobox.append_text("pm02")
        #self.t2_mon_combobox.append_text("pm0")
        #self.t2_mon_combobox.append_text("pm1")
        #self.t2_mon_combobox.append_text("pm2")
        #self.t2_mon_combobox.append_text("roi1")
        #self.t2_mon_combobox.append_text("roi2")
        #self.t2_mon_combobox.append_text("roi3")
        #self.t2_mon_combobox.append_text("roi4")
        #self.t2_mon_combobox.setChecked(0)
        
        self.t2_ref_mon_txt = QLabel("Reference monitor (?):")
        #self.t2_ref_mon_txt.set_alignment(0,0.5)
        self.t2_ref_mon_entry = QLineEdit()
        self.t2_img_start_txt = QLabel("Image begin (?):")
        #self.t2_img_start_txt.set_alignment(0,0.5)
        self.t2_img_start_entry = QLineEdit()
        self.t2_img_end_txt = QLabel("Image end (?):")
        #self.t2_img_end_txt.set_alignment(0,0.5)
        self.t2_img_end_entry = QLineEdit()
        self.t2_combine_XY = QLabel("Combination of X,Y translated images (n & n+1)")
        #self.t2_combine_XY.set_alignment(0,0.5)
        self.combine_XY = QCheckBox() #gtk.CheckButton()
        self.t2_ascii_out = QLabel("Save data as ASCII file")
        #self.t2_ascii_out.set_alignment(0,0.5)
        self.ascii_out = QCheckBox() #gtk.CheckButton()
        
        #self.t2_tooltip = gtk.Tooltips()
        self.t2_ref_mon_txt.setToolTip("Reference value of the monitor which corresponds to the highest value of the current machine. This is for the data normalization. If empty, the data will not be normalized")
        self.t2_img_start_txt.setToolTip("Starting image to be corrected. If empty the whole folder will be proceeded.")
        self.t2_img_end_txt.setToolTip("Ending image to be corrected. If empty the whole folder will be proceeded.")
        
        self.run_button = QPushButton("RUN")
        self.run_button.clicked.connect(self.batch_transform)
        #self.run_button.set_usize(30,0)
        #self.run_button.set_alignment(0.5,0.5)
        separator = QSplitter(Qt.Horizontal)#gtk.HSeparator()
        self.show_process = QLabel("Processing:")
        #self.show_process.set_alignment(0,0.5)
        self.show_process_info = QLabel()
        self.progressbar = QProgressBar(self) #gtk.ProgressBar()
        self.progressbar.setRange(0,100)
        
        #**** Packaging *************************************************************
        self.table_1.addWidget(self.t1_src_folder_txt, 1,1)
        self.table_1.addWidget(self.t1_src_path, 1,2)
        self.table_1.addWidget(self.t1_src_button, 1,3)
        self.table_1.addWidget(self.t1_des_folder_txt, 2,1)
        self.table_1.addWidget(self.t1_des_path, 2,2)
        self.table_1.addWidget(self.t1_des_button, 2,3)
        self.table_1.addWidget(self.t1_dark_img_txt, 3,1)
        self.table_1.addWidget(self.t1_dark_img_path, 3,2)
        self.table_1.addWidget(self.t1_dark_img_button, 3,3)
        
        self.table_2.addWidget(self.t2_detector_txt, 1,1)
        self.table_2.addWidget(self.t2_det_combobox, 1,2)
        self.table_2.addWidget(self.t2_monitor_txt, 2,1)
        self.table_2.addWidget(self.t2_mon_combobox, 2,2)
        self.table_2.addWidget(self.t2_ref_mon_txt, 3,1)
        self.table_2.addWidget(self.t2_ref_mon_entry, 3,2)
        self.table_2.addWidget(self.t2_img_start_txt, 4,1)
        self.table_2.addWidget(self.t2_img_start_entry, 4,2)
        self.table_2.addWidget(self.t2_img_end_txt, 5,1)
        self.table_2.addWidget(self.t2_img_end_entry, 5,2)
        self.table_2.addWidget(self.t2_combine_XY, 6,1)
        self.table_2.addWidget(self.combine_XY, 6,2)
        self.table_2.addWidget(self.t2_ascii_out, 7,1)
        self.table_2.addWidget(self.ascii_out, 7,2)
        
        self.table_2.addWidget(self.d1_specfile_txt,1,3)
        self.table_2.addWidget(self.d1_specfile_path,1,4)
        self.table_2.addWidget(self.d1_specfile_browse,1,5)
        self.table_2.addWidget(self.d1_center_txt,2,3)
        self.table_2.addWidget(self.d1_center,2,4)
        self.table_2.addWidget(self.d1_distance_txt,3,3)
        self.table_2.addWidget(self.d1_distance,3,4)
        
        
        self.page_batch_correction_container = QVBoxLayout() #gtk.VBox()
        self.page_batch_correction_container.addLayout(self.table_1)#, False, False, 10)
        self.page_batch_correction_container.addLayout(self.table_2)#, False, False, 10)
        self.page_batch_correction_container.addWidget(self.run_button)#, False, False, 10)
        #self.page_batch_correction_container.pack_start(separator, False, False, 10)
        self.table_3 = QGridLayout() #gtk.Table(1,2,False)
        self.table_3.addWidget(self.show_process, 1,1)
        self.table_3.addWidget(self.show_process_info,1,2)
        self.page_batch_correction_container.addLayout(self.table_3)#, False, False, 10)
        self.page_batch_correction_container.addWidget(self.progressbar)#, False, False, 10)
        self.page_batch_correction.addLayout(self.page_batch_correction_container)#, False, False, 20)
        
        ###################### PACK THE NOTEBOOK #####################################
        hbox.addWidget(self.notebook, 3)
        vbox.addLayout(hbox)#,True,True,0)
        ############################### Sliders ######################################
        #sld_box = gtk.Fixed()
        sld_box = QHBoxLayout()  #gtk.HBox(False,2)

        self.vmin_txt = QLabel("Vmin")
        #self.vmin_txt.set_alignment(0,0.5)
        #self.vmin_txt.set_justify(gtk.JUSTIFY_CENTER)
        self.vmax_txt = QLabel("Vmax")
        #self.vmax_txt.set_alignment(0,0.5)
        #self.vmax_txt.set_justify(gtk.JUSTIFY_CENTER)
        self.sld_vmin = QSlider(Qt.Horizontal)#gtk.HScale()
        self.sld_vmax = QSlider(Qt.Horizontal)#gtk.HScale()

        #self.sld_vmin.set_size_request(200,25)
        #self.sld_vmax.set_size_request(200,25)
        self.sld_vmin.setRange(0,self.vmax)
        self.sld_vmax.setRange(0,self.vmax)
        self.sld_vmax.setValue(0)#self.vmax)
        self.sld_vmin.setValue(0)
        self.sld_vmin.valueChanged.connect(self.scale_update)
        self.sld_vmax.valueChanged.connect(self.scale_update)

        #vmax_spin_adj         = gtk.Adjustment(self.vmax, 0, self.vmax_range, 0.5, 10.0, 0.0)
        self.vmax_spin_btn    = QSpinBox() #QSlider(Qt.Horizontal)   #gtk.SpinButton(vmax_spin_adj,1,1)
        #self.vmax_spin_btn.setRange(0,self.vmax_range)
        #self.vmax_spin_btn.setValue(self.vmax)
        #self.vmax_spin_btn.setTickInterval(0.5)
        #self.vmax_spin_btn.setSingleStep(10)
        #self.vmax_spin_btn.hide()
        #self.vmax_spin_btn.set_numeric(True)
        #self.vmax_spin_btn.set_wrap(True)
        #self.vmax_spin_btn.set_update_policy(gtk.UPDATE_IF_VALID)
        #self.vmax_spin_btn.set_size_request(80,-1)
        #self.vmax_spin_btn.set_alignment(0,0.5)
        self.vmax_spin_btn.valueChanged.connect(self.scale_update_spin)

        vmin_spin_adj = QSlider(Qt.Horizontal) #gtk.Adjustment(self.vmin, 0, self.vmax_range, 0.5, 10.0, 0.0)
        vmin_spin_adj.setRange(0,self.vmax_range)
        vmin_spin_adj.setValue(self.vmin)
        vmin_spin_adj.setTickInterval(0.5)
        vmin_spin_adj.setSingleStep(10)
        self.vmin_spin_btn    = QSpinBox() #gtk.SpinButton(vmin_spin_adj,1,1)
        #self.vmin_spin_btn.set_numeric(True)
        #self.vmin_spin_btn.set_wrap(True)
        #self.vmin_spin_btn.set_update_policy(gtk.UPDATE_IF_VALID)
        #self.vmin_spin_btn.set_size_request(80,-1)
        #self.vmax_spin_btn.set_alignment(0,0.5)
        self.vmin_spin_btn.valueChanged.connect(self.scale_update_spin)

        self.slider_reset_btn = QPushButton("Reset")
        #self.slider_reset_btn.set_size_request(50,30)
        #self.slider_reset_btn.set_alignment(1,0.5)
        self.slider_reset_btn.clicked.connect(self.reset_scale)

        sld_box.addWidget(self.vmin_txt)#,False,False,0)
        sld_box.addWidget(self.sld_vmin)#,False,False,0)
        sld_box.addWidget(self.vmin_spin_btn)#,False,False,0)
        sld_box.addWidget(self.vmax_txt)#,False,False,0)
        sld_box.addWidget(self.sld_vmax)#,False,False,0)
        sld_box.addWidget(self.vmax_spin_btn)#,False,False,0)
        sld_box.addWidget(self.slider_reset_btn)#,False,False,0)

        vbox.addLayout(sld_box)#,False,False,5)

        ################# Status bar #################################################
        self.status_bar = QHBoxLayout() #QStatusBar() #gtk.Table(1,22,True)

        self.x_pos_txt = QLabel("X: ")
        self.y_pos_txt = QLabel("Y: ")
        self.z_pos_txt = QLabel("Z: ")
        self.x_pos = QLabel()
        self.y_pos = QLabel()
        self.z_pos = QLabel()
        
        self.del_pos_txt = QLabel()
        self.eta_pos_txt = QLabel()
        self.phi_pos_txt = QLabel()
        self.kphi_pos_txt= QLabel()
        self.chi_pos_txt = QLabel()
        self.nu_pos_txt  = QLabel()
        self.mu_pos_txt  = QLabel()
        self.time_pos_txt= QLabel()
        
        self.del_pos = QLabel()
        self.eta_pos = QLabel()
        self.phi_pos = QLabel()
        self.kphi_pos = QLabel()
        self.chi_pos = QLabel()
        self.nu_pos  = QLabel()
        self.mu_pos  = QLabel()
        self.time_pos  = QLabel()
        
        self.status_bar_grid = QGridLayout()
        self.status_bar_grid.addWidget(self.x_pos_txt,1,1)
        self.status_bar_grid.addWidget(self.x_pos,1,2)
        self.status_bar_grid.addWidget(self.y_pos_txt,1,3)
        self.status_bar_grid.addWidget(self.y_pos,1,4)
        self.status_bar_grid.addWidget(self.z_pos_txt,1,5)
        self.status_bar_grid.addWidget(self.z_pos,1,6)
        
        
        self.status_bar_grid.addWidget(self.eta_pos_txt,1,7)
        self.status_bar_grid.addWidget(self.eta_pos,1,8)
        self.status_bar_grid.addWidget(self.chi_pos_txt,1,9)
        self.status_bar_grid.addWidget(self.chi_pos,1,10)
        self.status_bar_grid.addWidget(self.phi_pos_txt,1,11)
        self.status_bar_grid.addWidget(self.phi_pos,1,12)
        self.status_bar_grid.addWidget(self.nu_pos_txt,1,13)
        self.status_bar_grid.addWidget(self.nu_pos,1,14)
        self.status_bar_grid.addWidget(self.del_pos_txt,1,15)
        self.status_bar_grid.addWidget(self.del_pos,1,16)
        self.status_bar_grid.addWidget(self.time_pos_txt,1,17)
        self.status_bar_grid.addWidget(self.time_pos,1,18)
        self.status_bar.addLayout( self.status_bar_grid )
        vbox.addLayout(self.status_bar)#,False,False,0)
        self.ctrlw.setLayout( vbox )
        self.setCentralWidget( self.ctrlw )
        #print("Came in displaying everything ")
        self.show()#_all()
        #self.popup_info("Info","Kindly select a detector first")
        #self.progressbar.hide()
        #self.data_loading.hide()
        #self.d1_center.hide()
        #self.d1_center_txt.hide()
        #self.d1_distance.hide()
        #self.d1_distance_txt.hide()
        #self.d1_specfile_browse.hide()
        #self.d1_specfile_path.hide()
        #self.d1_specfile_txt.hide()

######################################## Definitions ########################################################################
    def instrument_config(self, QWidget):
        global _SPEC_IMG_COL ,_IMG_NUMBER_SEPARATOR, _IMG_FORMAT ,_DETECTOR_NAME , _DETECTOR_PIXELS ,_PIXEL_SIZE ,_DETECTOR_ROTATIONS, _DETECTOR_MOTORS, _DETECTOR_ROTATION_DIRECTIONS,_DETECTOR_PIXEL_DIRECTIONS, _SAMPLE_ROTATIONS , _SAMPLE_ROTATION_DIRECTIONS, _SAMPLE_ROTATION_MOTORS , _XRAY_INCIDENT_BEAM_DIRECTION
        
        #self.popup_info('info', "When the instrument is configured, please quit the application and restart it to take into account the changes. This feature will be changed in future release.")
        
        filename = join(DEVA_dir, "Instrument_configuration.DEVA")
        os.system("gedit  %s >/dev/null < /dev/null " % (filename),)

    def pro_format_coord(self,x,y):
        if abs(x)>10000:
            fmtx = '%.2e'%x
        else:
            fmtx = '%.2f'%x
        if abs(y)>10000:
            fmty = '%.2e'%y
        else:
            fmty = '%.2f'%y
        form = 'x=%s, y=%s'%(str(fmtx),str(fmty))
        return form
        
    def pole_format_coord(self, x,y):
        r = N.sqrt(x**2 + y**2)
        chi = N.degrees(r)
        phi = N.arctan2(y,x)
        phi = N.degrees(phi)
        if phi <0:
            phi += 360
        out = "Chi = %.2f, Phi = %.2f"%(chi, phi)
        return out

    def init_image(self):
        self.ax.clear()
        self.cax.clear()
        self.cax2.clear()
        self.ax.add_patch(self.rect)
        self.ax.add_patch(self.roi_rect)
        self.IMG_ZOOMED = False
        self.img = self.ax.imshow(self.data,origin='lower',vmin=self.vmin, vmax=self.vmax, cmap=cmap, interpolation='nearest',aspect='auto', extent=self.MAIN_EXTENT)
        
        if self.log_scale == 0:
            clabel = r'$Intensity\ (Counts\ per\ second)$'
        else:
            clabel = r'$Log_{10}\ (Counts\ per\ second)\ [arb.\ units]$'
        self.cb  = self.fig.colorbar(self.img, cax = self.cax, format="%.f")#format=fm
        self.cb.set_label(clabel, fontsize=12)
        self.cb.locator = MaxNLocator(nbins=6)
        #Horizontal colorbar
        self.cb2  = self.fig.colorbar(self.img, cax = self.cax2, orientation='horizontal', format="%.f")#format=fm
        self.cb2.set_label(clabel, fontsize=12)
        self.cb2.locator = MaxNLocator(nbins=6)
        self.cb2.ax.set_visible(False)
        
        if self.tth_chi_space_btn.isChecked():
            self.MAIN_XLABEL.set_text(r"$2\theta\ (deg.)$")
            self.MAIN_YLABEL.set_text(r"$\chi\ (deg.)$")
        
        elif self.hk_space_btn.isChecked():
            self.MAIN_XLABEL.set_text(r"$H\ (R.L.U)$")
            self.MAIN_YLABEL.set_text(r"$K\ (R.L.U)$")
        elif self.hl_space_btn.isChecked():
            self.MAIN_XLABEL.set_text(r"$H\ (R.L.U)$")
            self.MAIN_YLABEL.set_text(r"$L\ (R.L.U)$")
        elif self.kl_space_btn.isChecked():
            self.MAIN_XLABEL.set_text(r"$K\ (R.L.U)$")
            self.MAIN_YLABEL.set_text(r"$L\ (R.L.U)$")
        elif self.q_space_btn.isChecked():
            self.MAIN_XLABEL.set_text(r"$Qy\ (nm^{-1})$")
            self.MAIN_YLABEL.set_text(r"$Qz\ (nm^{-1})$")
        else:
            self.MAIN_XLABEL.set_text("X (pixel)")
            self.MAIN_YLABEL.set_text("Y (pixel)")
        
        if self.graph_aspect == True:
            self.ax.set_aspect('equal')
        else:
            self.ax.set_aspect('auto')
            
        self.IMG_INIT = True
        #self.cursor = Cursor(self.ax, color='k', linewidth=1, useblit=True)
        
    def change_aspect_ratio(self,w):
        self.graph_aspect = not (self.graph_aspect)
        if self.graph_aspect == True:
            self.ax.set_aspect('equal')
        else:
            self.ax.set_aspect('auto')
        self.canvas.draw()

    def popup_info(self,info_type,text):
        """ info_type = WARNING, INFO, QUESTION, ERROR """
        self.warning=QMessageBox.about(self, info_type.upper(), text)

    def load_calibration(self, widget):
        """ load a pre-calib file , PONI file """
        if 1==1 or self.loadcalibtb.get_active():
            dialog = QFileDialog.getOpenFileName(self, "Select an EDF folder", self.current_folder, "Poni files (*.poni)")
            if dialog!='' and type(dialog)==type((1,2)): #dialog.exec_()==QDialog.Accepted: #response == gtk.RESPONSE_OK:
                dialog = dialog[0]
                if dialog=='' or len(dialog)==0:
                    return 
                self.ponifile = dialog#.get_filename().decode('utf8')
                print( "Calibration file is: ",self.ponifile)
                self.azimuthalIntegration = pyFAI.load(self.ponifile)

                self.calibrated = True
                self.calibrated_quantitative = True
                print( "is calibrated? ",self.calibrated)
                # self.geometry_manual.setChecked(False)
                self.canvas.draw()
                s = os.path.basename(self.ponifile)
                self.popup_info("info","This detector is calibrated with the PONI file %s!!!"%s)
            else:
                print("File is not selected. Kindly select a file")
                pass
            #dialog.destroy()
            self.wavelength = self.azimuthalIntegration.wavelength
            self.distance   = self.azimuthalIntegration.dist
            self.energy     = xrayutilities.lam2en(self.wavelength)/1e10
            self.direct_beam= [self.azimuthalIntegration.poni2/_PIXEL_SIZE, self.azimuthalIntegration.poni1/_PIXEL_SIZE]
            self.geometry_energy.setText(str(self.energy))
            self.geometry_distance.setText(str(self.distance))
            self.geometry_direct_beam.setText(str(self.direct_beam[0])+","+str(self.direct_beam[1]))
            
            tmp_dir  = DEVA_dir
            geo_name = join(tmp_dir,"Geo_config.DEVA")
            geo_file = open(geo_name,"w")
            in_plane = self.geometry_substrate_inplane.text()
            out_of_plane = self.geometry_substrate_outplane.text()
            if in_plane != "" and out_of_plane != "":
                in_plane = in_plane.split()
                self.in_plane = N.asarray([int(i) for i in in_plane])
                out_of_plane = out_of_plane.split()
                self.out_of_plane = N.asarray([int(i) for i in out_of_plane])
            content  =""
            content += "ENERGY="+str(self.energy)+"\n"
            content += "DISTANCE="+str(self.distance)+"\n"
            content += "DIRECT_BEAM_XY="+str(self.direct_beam[0])+","+str(self.direct_beam[1])+"\n"
            if self.has_substrate:
                content += "HAS_SUBSTRATE=YES\n"
                content += "SUBSTRATE="+substrate+"\n"
            else:
                content += "HAS_SUBSTRATE=NO\n"
            
            content += "INPLANE="+str(self.in_plane[0])+" "+str(self.in_plane[1])+" "+str(self.in_plane[2])+"\n"
            content += "OUTPLANE="+str(self.out_of_plane[0])+" "+str(self.out_of_plane[1])+" "+str(self.out_of_plane[2])+"\n"
            
            geo_file.write(content)
            geo_file.close()
            print( content )
        else:
            self.calibrated = False
            self.calibrated_quantitative = False

    def manual_calibration(self,widget):
        """Checking input data for geometry setup
        """
        distance = self.geometry_distance.text()
        energy   = self.geometry_energy.text()
        direct_beam = self.geometry_direct_beam.text()
        if len(distance)==0:
            distance = 0
        if len(energy)==0:
            energy = 0

        self.distance = float(distance)
        self.energy   = float(energy)
        direct_beam = direct_beam.split(",")
        self.direct_beam = [float(direct_beam[0]),float(direct_beam[1])]
        from scipy.constants import h,c,e
        poni1 = self.direct_beam[1]*_PIXEL_SIZE
        poni2 = self.direct_beam[0]*_PIXEL_SIZE
        self.wavelength = h*c/e/self.energy

        self.qconv.energy = self.energy
                
        if self.UB_MATRIX_LOAD:
            self.qconv.UB = self.UB_MATRIX
            self.experiment = xrayutilities.HXRD([1,0,0],[0,0,1], en=self.energy, qconv=self.qconv)
        else:
            substrate = self.geometry_substrate.currentText()
            if substrate == "-- other":
                if self.geometry_substrate_other.text() == "":
                    self.experiment = xrayutilities.HXRD([1,0,0],[0,0,1], en=self.energy, qconv=self.qconv)
                    self.has_substrate = False
                else:
                    self.has_substrate = True
                    substrate = self.geometry_substrate_other.text()
            else:
                self.has_substrate = True
            if self.has_substrate:
                command = "self.substrate = xrayutilities.materials."+substrate
                exec(command)
                in_plane = self.geometry_substrate_inplane.text()
                out_of_plane = self.geometry_substrate_outplane.text()
                if in_plane != "" and out_of_plane != "":
                    in_plane = in_plane.split()
                    self.in_plane = N.asarray([int(i) for i in in_plane])
                    out_of_plane = out_of_plane.split()
                    self.out_of_plane = N.asarray([int(i) for i in out_of_plane])
                    #self.has_orientation_matrix = True
                    self.experiment = xrayutilities.HXRD(self.substrate.Q(self.in_plane),self.substrate.Q(self.out_of_plane), en=self.energy, qconv=self.qconv)
                else:
                    #self.has_orientation_matrix = False
                    self.experiment = xrayutilities.HXRD(self.substrate.Q(1,0,0),self.substrate.Q(0,0,1), en=self.energy, qconv=self.qconv)
                    
        
        if not self.calibrated_quantitative:
            self.azimuthalIntegration = pyFAI.azimuthalIntegrator.AzimuthalIntegrator(dist=self.distance,
                                                                                    poni1=poni1,
                                                                                    poni2=poni2,
                                                                                    rot1=None,
                                                                                    rot2=None,
                                                                                    rot3=None,
                                                                                    pixel1=_PIXEL_SIZE,
                                                                                    pixel2=_PIXEL_SIZE,
                                                                                    wavelength=self.wavelength)
            #self.experiment.Ang2Q.init_area('z+','y-', cch1=direct_beam[1], cch2=direct_beam[0], Nch1=Npx1,Nch2=Npx2, pwidth1=px1,pwidth2=px2, distance=distance, detrot=detrot, tiltazimuth=0, tilt=0)
        
            self.calibrated=True
        MSSG = "Your parameters have been taken into account.\nEnergy = %s eV\nDistance = %s m\nDirect beam position: %s,%s\n"%(str(self.energy),str(self.distance),str(self.direct_beam[0]),str(self.direct_beam[1]))
        if self.UB_MATRIX_LOAD:
            MSSG+= "\nYou have imported a UB matrix. If you donot want to use this UB matrix anymore, click the browse button again.\n\nYour actual UB matrix is:\n%s"%str(self.UB_MATRIX)
        else:
            if self.has_substrate:
                MSSG+= "\nYou do not have a UB matrix, you have to define your substrate material and it's orientation. This information will be considered to calculate the orientation matrix\n"
                MSSG+= "\nYour actual choise:\nSubstrate material: %s\nIn-plane direction: %s\nOut-of-plane direction: %s"%(str(substrate), str(self.in_plane), str(self.out_of_plane))
            else:
                MSSG+= "\nYou donot have a UB matrix, nor a substrate. A UB equal to unity will be applied\n"
        self.popup_info("info",MSSG)
        tmp_dir  = DEVA_dir
        geo_name = join(tmp_dir,"Geo_config.DEVA")
        geo_file = open(geo_name,"w")
        in_plane = self.geometry_substrate_inplane.text()
        out_of_plane = self.geometry_substrate_outplane.text()
        if in_plane != "" and out_of_plane != "":
            in_plane = in_plane.split()
            self.in_plane = N.asarray([int(i) for i in in_plane])
            out_of_plane = out_of_plane.split()
            self.out_of_plane = N.asarray([int(i) for i in out_of_plane])
        else:
            self.in_plane = ""
            self.out_of_plane = ""
        content  =""
        content += "ENERGY="+str(self.energy)+"\n"
        content += "DISTANCE="+str(self.distance)+"\n"
        content += "DIRECT_BEAM_XY="+str(self.direct_beam[0])+","+str(self.direct_beam[1])+"\n"
        if self.has_substrate:
            content += "HAS_SUBSTRATE=YES\n"
            content += "SUBSTRATE="+substrate+"\n"
        else:
            content += "HAS_SUBSTRATE=NO\n"
        
        content += "INPLANE="+str(self.in_plane[0])+" "+str(self.in_plane[1])+" "+str(self.in_plane[2])+"\n"
        content += "OUTPLANE="+str(self.out_of_plane[0])+" "+str(self.out_of_plane[1])+" "+str(self.out_of_plane[2])+"\n"
        
        geo_file.write(content)
        geo_file.close()
            
    def calculation_angular_coordinates(self):
        # if self.geometry_corrected:
        self.tableTwoTheta = self.azimuthalIntegration.twoThetaArray((self.data.shape[0],self.data.shape[1]))
        self.tableChi      = self.azimuthalIntegration.chiArray((self.data.shape[0],self.data.shape[1]))
        self.tableChi      = N.degrees(self.tableChi)-90
        self.table_dSpace = self.azimuthalIntegration.wavelength / (2*N.sin(self.tableTwoTheta/2.0)) * 1e10 #d in Angstrom
        self.tableTwoTheta = N.degrees(self.tableTwoTheta)
        
    def set_detector(self,widget, det_type):
        if det_type.upper()=="D5":
            print( "D5 selected" )
            #self.detm.set_label("XPAD D5")
            self.detector_type = "D5"

        elif det_type.upper()=="S70":
            print( "S70 selected" )
            #self.detm.set_label("XPAD S70")
            self.detector_type = "S70"
            self.xDim1 = 560
            self.yDim1 = 0
            self.yDim0 = 120
        elif det_type.upper()=="Pixirad":
            print( "D1 selected" )
            #self.detm.set_label("XPAD D1")
            self.detector_type = "Pixirad"
            self.xDim1 = 476
            self.yDim1 = 4096
        elif det_type.upper()=="WOS":
            print( "D1 selected" )
            #self.detm.set_label("XPAD D1")
            self.detector_type = "WOS"
            self.xDim1 = 598
            self.yDim1 = 1154

        elif det_type.upper()=="D1":
            print( "D1 selected" )
            #self.detm.set_label("XPAD D1")
            self.detector_type = "D1"
            self.xDim1 = 577
            self.yDim1 = 913

        self.init_image()
        self.canvas.draw()
        return 0
    
    def set_manip(self,widget, manip_type):
        if manip_type.upper()=="GONIO":
            ## -----------------------------------------
            #self.manipm.setText("Gonio")
            self.experiment_type = "GONIO"

        elif manip_type.upper()=="GISAXS":
            ## -----------------------------------------
            #self.manipm.setText("GISAXS")
            self.experiment_type = "GISAXS"
        self.canvas.draw()
        return 0

    def check_azimuthal_integrator(self):
        if not self.calibrated_quantitative:
            if self.experiment_type=="GISAXS":
                if hasattr(self,'det_outer'): self.nu = 0.0
                self.delta = 0.0
            #if _DETECTOR_ROTATIONS==2:            
            rot1 = N.radians(self.pyFAI_rot1)*(-1.)
            rot2 = N.radians(self.pyFAI_rot2)*(-1.)
            rot3 = N.radians(90-self.chi)
            self.azimuthalIntegration.rot1 = rot1
            self.azimuthalIntegration.rot2 = rot2
            self.azimuthalIntegration.rot3 = rot3
        self.calculation_angular_coordinates()
        
    def get_dark(self,w):
        if self.use_dark_tb.isChecked():
            #self.DARK_CORRECTION = not self.DARK_CORRECTION
            if self.SELECTED_IMG_NUM != None:
                self.DARK_DATA = self.fabioIMG.data
                self.DARK_CORRECTION = True
            else:
                self.use_dark_tb.setChecked(False)
                self.DARK_CORRECTION = False
        else:
            self.DARK_CORRECTION = False
        print( "Use Dark image: ",self.DARK_CORRECTION )
        
    def what_is_this_manip(self):
        """ Check manip type by the spec file"""
        kappapsic = re.compile(r'kappapsic')
        fourc     = re.compile(r'kappa')
        gisaxs    = re.compile(r'gisaxs')
        saxsext   = re.compile(r'saxsext')
        manip     = [kappapsic, fourc, gisaxs, saxsext]
        manip_list= ['kappapsic', 'fourc', 'gisaxs', 'saxsext']
        spec_file = os.path.basename(self.SPEC_FILE)
        print( "Spec file: ",spec_file )
        for m in range(len(manip)):
            if manip[m].findall(spec_file) !=[]:
                self.manip = manip_list[m]
                break
            else:
                self.manip = "Unknown"
        
    def read_header(self,header):
        self.counter = get_counters(header)
        self.motor = get_motors(header)
        motor_mne = self.motor.keys()
        self.count_time = self.counter['sec']
        #import pdb; pdb.set_trace()
        self.det_inner = self.motor[_DETECTOR_MOTORS[-1]]
        self.delta = self.det_inner
        self.th    = self.motor[_SAMPLE_ROTATION_MOTORS[0]]
        self.chi   = self.motor[_SAMPLE_ROTATION_MOTORS[1]]
        self.phi   = self.motor[_SAMPLE_ROTATION_MOTORS[2]]
        if _DETECTOR_ROTATIONS==2:
            self.det_outer = self.motor[_DETECTOR_MOTORS[0]]
            self.nu = self.det_outer
        if len(_DETECTOR_ROTATION_DIRECTIONS) == 1:
            rotd1 = _DETECTOR_ROTATION_DIRECTIONS[0]
            if 'z' in rotd1:
                self.pyFAI_rot1 = self.det_inner
                self.pyFAI_rot2 = 0
            elif 'y' in rotd1:
                self.pyFAI_rot2 = self.det_inner
                self.pyFAI_rot1 = 0
        if len(_DETECTOR_ROTATION_DIRECTIONS) == 2:
            rotd1 = _DETECTOR_ROTATION_DIRECTIONS[0]
            rotd2 = _DETECTOR_ROTATION_DIRECTIONS[1]
            if 'z' in rotd1:
                self.pyFAI_rot1 = self.det_outer
            elif 'y' in rotd1:
                self.pyFAI_rot2 = self.det_outer
            if 'z' in rotd2:
                self.pyFAI_rot1 = self.det_inner
            elif 'y' in rotd2:
                self.pyFAI_rot2 = self.det_inner
       
    def treevViewSelection(self, signal):#, old, new):
        #import pdb; pdb.set_trace()
        global g_is_image_loaded
        file_path=self.MODEL.filePath(signal)
        print(file_path)
        #print("Came in treeviewselction")
        g_is_image_loaded = True
        self.on_changed_edf(file_path, '', '')

    def on_changed_edf(self,widget,row,col):
        """ Change EDF by double clicking on the file name """
        self.clear_notes()
        self.init_image()
        #model = widget.get_model()
        #print( "TRA CHOOSE: ",model
        #print( "ROW: ",row,len(row)
        #print( "MODEL[ROW]: ",model[row]
        if os.path.isdir(widget):
            self.popup_info("ERROR","This is a directory. No action will be taken")
            print("This is a directory. No action will be taken")
            return
        if widget[ len(widget)-4: ]!=".edf":
            self.popup_info("ERROR","This is not an image file. Kindly check it.")
            print("This is not an image file. Kindly check it.")
            return
        path = widget.split(os.sep)
        self.edf_choosen = path[len(path)-1] #model[row][0]
        #Getting prefix:
        tmp=self.edf_choosen.split("_")
        if len(tmp)>2:
            prefix=""
            for p in range(len(tmp)-1):
                prefix=prefix+tmp[p]+"_"
        else:
            prefix=tmp[0]+"_"
        self.current_prefix = prefix        
        self.edf_folder = self.current_folder
        if len(row)==2:
            #print( "Selected Directory: ",model[row[0]][0]
            self.edf_folder = join(self.current_folder,model[row[0]][0])
            #print( "edf_folder: ",self.edf_folder
        self.edf = join(self.edf_folder,self.edf_choosen)
        #self.edf_pos.set_text("EDF choosen:  %s"%self.edf_choosen)
        try:
            self.MAIN_TITLE.setText(self.edf_choosen)#,fontsize=18)
        except:
            self.MAIN_TITLE = self.ax.set_title(self.edf_choosen,fontsize=12)
        #self.MAIN_TITLE = self.fig.suptitle(self.edf_choosen)
        ### Data Loading #########
        print("  >>>>>>   ",self.edf)
        self.fabioIMG = fabio.open(self.edf)
        self.header = self.fabioIMG.header
        self.Dim1 = self.fabioIMG.data.shape[0]
        self.Dim2 = self.fabioIMG.data.shape[1]
        if self.detector_type == "S70":
            self.fabioIMG.data = N.flipud(self.fabioIMG.data)
        #print( self.header
        self.read_header(self.header)
        self.time_pos_txt.setText("Seconds: ")
        self.time_pos.setText("%d"%self.count_time)
        if hasattr(self,'det_outer'):
            self.nu_pos_txt.setText(_DETECTOR_MOTORS[0])
            self.nu_pos.setText("%.2f"%self.det_outer)
        self.del_pos_txt.setText(_DETECTOR_MOTORS[-1])
        self.eta_pos_txt.setText(_SAMPLE_ROTATION_MOTORS[0])
        self.chi_pos_txt.setText(_SAMPLE_ROTATION_MOTORS[1])
        self.phi_pos_txt.setText(_SAMPLE_ROTATION_MOTORS[2])

        self.del_pos.setText("%.2f"%self.det_inner)
        self.eta_pos.setText("%.2f"%self.th)
        self.chi_pos.setText("%.2f"%self.chi)
        self.phi_pos.setText("%.2f"%self.phi)
        gc.collect() # Clear unused variables
        #self.data = self.fabioIMG.data
        self.plot_data()
        spliter = _IMG_NUMBER_SEPARATOR
        num = self.edf_choosen.split(spliter)[-1]
        num = num.split(".")[0]
        try:
            self.SELECTED_IMG_NUM = int(num)
        except:
            self.SELECTED_IMG_NUM = 1
        if self.SPEC_IS_LOADED and len(self.SPEC_IMG) != 0:
            self.check_and_update_scan_slider()#Check the scan number and image number --> set the spin button for scan num and the slider for img num
        return

    def load_specFile(self,widget):
        #dialog = gtk.FileChooserDialog("Select spec file",None,gtk.FILE_CHOOSER_ACTION_OPEN,(gtk.STOCK_CANCEL,gtk.RESPONSE_CANCEL, gtk.STOCK_OPEN, gtk.RESPONSE_OK))
        dialog = QFileDialog.getOpenFileNames(self, "Select spec file", self.current_folder, "")# "Poni files (*.poni)")
        #dialog.set_current_folder(self.current_folder)
        #response = dialog.run()
        print("dialog values >>>>>>>> ",dialog)
        if dialog!='' and type(dialog)==type((1,2)): #response == gtk.RESPONSE_OK:
            dialog = dialog[0]
            if dialog=='' or len(dialog)==0:
                return 
            file_choosen = dialog#.get_filename().decode('utf8')
            self.scan_slider_path.setText(file_choosen)
            self.SPEC_FILE = file_choosen
        else:
            print("Kindly select a valid file")
            return
        #dialog.destroy()
        while gtk.events_pending():
            gtk.main_iteration()
        if self.detector_type=="D1":
            self.SPEC_DATA = Read_Spec_D1(self.SPEC_FILE)
        else:
            self.SPEC_DATA = xrayutilities.io.SPECFile(self.SPEC_FILE)
        self.update_spec_data()
        
        return
        
    def load_UBfile(self,widget):
        if self.geometry_browse_UB.isChecked():
            dialog = QFileDialog.getOpenFileName(self,"Select a PONI file", self.current_folder,"Python Files (*.poni);;")#, options=options)
            print(" >>>>> ",dialog)
            if dialog!='' and type(dialog)==type((1,2)): #response == gtk.RESPONSE_OK:
                dialog = dialog[0]
                if dialog=='' or len(dialog)==0:
                    return 
                file_choosen = dialog#.get_filename().decode('utf8')
                self.UB_FILE = file_choosen
            else:
                print("Kindly select a Poni file.")
                self.popup_info("Info", "Kindly select a Poni file.")
                return
            ## Need to replace the loadtxt file with poni file reading
            self.UB_MATRIX = "" #N.loadtxt(self.UB_FILE)
            self.UB_MATRIX_LOAD=True
            #self.qconv.UB=self.UB_MATRIX
            #self.experiment = xrayutilities.HXRD([1,0,0],[0,0,1], en=self.energy, qconv=self.qconv)
            self.geometry_browse_UB.setText("PONI imported")
            print( "PONI file: ",self.UB_FILE )
            print( "PONI matrix: \n",self.UB_MATRIX )
        else:
            self.UB_MATRIX_LOAD = False
            self.geometry_browse_UB.setText("Browse PONI file")
        return
        
    def update_spec_data(self):
        self.SPEC_IMG = []
        self.SPEC_SCAN_LIST = self.SPEC_DATA.scan_list
        for i in range(len(_SPEC_IMG_COL)):
            if _SPEC_IMG_COL[i] not in self.SPEC_SCAN_LIST[0].colnames:
                img_col_found = False
                continue
            else:
                img_col_found = True
                self.IMG_COL  = _SPEC_IMG_COL[i]
                break
        if not img_col_found:
            img_col = ''
            for imc in _SPEC_IMG_COL:
                img_col+='%s, '%imc
            msg = "Spec file does not containt the image field. The current image number counter is %s. Please check the instrument configuration again."%img_col
            self.popup_info("error",msg)    
            return
        else:
            first_scan_num = self.SPEC_SCAN_LIST[0].nr
            last_scan_num = self.SPEC_SCAN_LIST[-1].nr
            
            self.SPEC_SCAN_RANGE = (first_scan_num, last_scan_num)
            #print( first_scan_num, last_scan_num
            self.SPEC_SCAN_NUM_LIST = []
            self.SPEC_ALL_MOTORS_LIST = []
            for i in range(len(self.SPEC_SCAN_LIST)):
                #item = self.SPEC_SCAN_LIST[i]
                self.SPEC_SCAN_NUM_LIST.append(self.SPEC_SCAN_LIST[i].nr)
                self.SPEC_ALL_MOTORS_LIST.append(self.SPEC_SCAN_LIST[i].colnames[0].upper())
                if self.SPEC_SCAN_LIST[i].scan_status == 'OK':
                    self.SPEC_SCAN_LIST[i].ReadData()
                    this_img_list = self.SPEC_SCAN_LIST[i].data[self.IMG_COL]
                    this_img_list = this_img_list.astype('int')
                else:
                    this_img_list = []
                self.SPEC_IMG.append(this_img_list)
                
            self.SPEC_IS_LOADED = True
            self.scan_slider_spinButton.set_range(self.SPEC_SCAN_RANGE[0], self.SPEC_SCAN_RANGE[1])
            
            #print( "SPEC_IMG len: ",len(self.SPEC_IMG)
            #print( "SPEC SCAN NUM LIST: ",self.SPEC_SCAN_NUM_LIST
            self.check_and_update_scan_slider()
    
    def check_skipped_motors(self):
        self.SPEC_SKIPPED_MOTORS = []
        if self.scan_slider_skip_chi.isChecked():
            self.SPEC_SKIPPED_MOTORS.append("Chi".upper())
        if self.scan_slider_skip_del.isChecked():
            self.SPEC_SKIPPED_MOTORS.append("Del".upper())
        if self.scan_slider_skip_eta.isChecked():
            self.SPEC_SKIPPED_MOTORS.append("Eta".upper())
        if self.scan_slider_skip_phi.isChecked():
            self.SPEC_SKIPPED_MOTORS.append("Phi".upper())
        if self.scan_slider_skip_tsz.isChecked():
            self.SPEC_SKIPPED_MOTORS.append("Tsz".upper())
        if self.scan_slider_skip_tox.isChecked():
            self.SPEC_SKIPPED_MOTORS.append("Tox".upper())
        if self.scan_slider_skip_toy.isChecked():
            self.SPEC_SKIPPED_MOTORS.append("Toy".upper())
        if self.scan_slider_skip_rox.isChecked():
            self.SPEC_SKIPPED_MOTORS.append("Rox".upper())
        if self.scan_slider_skip_roy.isChecked():
            self.SPEC_SKIPPED_MOTORS.append("Roy".upper())
        if self.scan_slider_skip_rien.isChecked():
            self.SPEC_SKIPPED_MOTORS.append("Rien".upper())
    
    def update_scan_slider(self,widget):
        self.update_scan_slider_now()
    
    def get_scan_data(self,i,edf_base):
        edf    = join(self.edf_folder, edf_base)
        data   = fabio.open(edf).data
        header = fabio.open(edf).header
        stdout.write("\r Loading %s"%edf_base)
        stdout.flush()
        return (i,data,header)
        
    def update_scan_slider_now(self):
        actual_scan_num = self.scan_slider_spinButton.value()
        actual_scan_num = int(actual_scan_num)
        self.check_skipped_motors()#To get the list of skipped motors
        print( "Actual scan number: ",actual_scan_num )
        #print( "Skipped scans: ",self.SPEC_SKIPPED_MOTORS
        #print( "First scan: %d, Last scan: %d"%(self.SPEC_SCAN_LIST[0].nr, self.SPEC_SCAN_LIST[-1].nr)
        for i in range(len(self.SPEC_SCAN_LIST)):
            #scan_motor = self.SPEC_ACTUAL_SCAN.colnames[0].upper()
            if (actual_scan_num == self.SPEC_SCAN_LIST[i].nr) and (self.SPEC_SCAN_LIST[i].colnames[0].upper() not in self.SPEC_SKIPPED_MOTORS):
                self.SPEC_ACTUAL_SCAN = self.SPEC_SCAN_LIST[i]
                self.SPEC_SCAN_MOTOR_NAME = self.SPEC_ACTUAL_SCAN.colnames[0]
                #print( "This scan motor: ",self.SPEC_SCAN_MOTOR_NAME
                break
            else:
                continue
        #print( "Actual scan number: ", self.SPEC_ACTUAL_SCAN.nr
        
        #self.SPEC_ACTUAL_SCAN.ReadData()#All scan are Data Ready
        #self.SPEC_ACTUAL_SCAN_IMG = []
        self.SPEC_ACTUAL_SCAN_IMG = self.SPEC_ACTUAL_SCAN.data[self.IMG_COL]
        self.SPEC_ACTUAL_SCAN_IMG = self.SPEC_ACTUAL_SCAN_IMG.astype('int')
        #print( "Actual scan images: ",self.SPEC_ACTUAL_SCAN_IMG
        actual_img_num  = self.SELECTED_IMG_NUM
        #actual_img_num = self.scan_slider_imgSlider.get_value()
        #actual_img_num = int(actual_img_num)
        print( "Current edf folder: ",self.edf_folder)
        print( "Current edf prefix: ",self.current_prefix)
        self.current_img_list=self.store_img[self.edf_folder][self.current_prefix]
        
        if (actual_img_num == None) or (actual_img_num not in self.SPEC_ACTUAL_SCAN_IMG):
            actual_img_num = self.SPEC_ACTUAL_SCAN_IMG[0]
        
        for k in self.store.keys():
            if actual_img_num in self.current_img_list:
                self.edf_folder = k
                break
            else:
                continue
        #print( "EDF folder: ",self.edf_folder
        #while gtk.events_pending():
            #gtk.main_iteration()
        
        if self.SPEC_ACTUAL_SCAN_IMG[0] == self.SPEC_ACTUAL_SCAN_IMG[-1]:
            self.scan_slider_imgSlider.setRange(self.SPEC_ACTUAL_SCAN_IMG[0], self.SPEC_ACTUAL_SCAN_IMG[-1]+1)
        else:
            self.scan_slider_imgSlider.setRange(self.SPEC_ACTUAL_SCAN_IMG[0], self.SPEC_ACTUAL_SCAN_IMG[-1])
        
        self.scan_slider_imgSlider.setValue(actual_img_num)#This will call the plot_scan too
        #print( "Actual image number: ",actual_img_num
        #print( self.SPEC_ACTUAL_SCAN_IMG
        self.SPEC_ACTUAL_SCAN_DATA = []
        self.SPEC_ACTUAL_SCAN_HEADER=[]
        #try:
        self.SPEC_ACTUAL_SCAN_IMG_NAMES = select_files_from_table(self.TABLE_STORE[self.edf_folder], self.current_prefix, self.SPEC_ACTUAL_SCAN_IMG.min(), self.SPEC_ACTUAL_SCAN_IMG.max())
        
        img_list = self.SPEC_ACTUAL_SCAN_IMG_NAMES
        #print( "Actual images: ",img_list
        #"""
        print( "Loading data for this scan %d ..."%self.SPEC_ACTUAL_SCAN.nr)
        
        #Getting data by multithreading
        threads = []
        for i in range(len(img_list)):
            t = get_scan_data_thread(i,self.edf_folder, img_list[i])
            threads.append(t)
            t.start()
        output = []
        for t in threads:
            t.join()
            output.append(t.Data)
        while len(threads)>0:
            threads.pop()
       # gobject.timeout_add(100, self._callback)
        output.sort()
        self.SPEC_ACTUAL_SCAN_DATA   = [o[1] for o in output]
        self.SPEC_ACTUAL_SCAN_HEADER = [o[2] for o in output]
        # print( "Don't worry about the order print(ed above."
        print( "End.")
        #except:
            #pass
            #self.popup_info("warning","Attention: Data not found for this scan.")
        
        #self.SPEC_SCAN_MOTOR_NAME = self.SPEC_ACTUAL_SCAN.colnames[0]
        self.SPEC_SCAN_MOTOR_DATA = self.SPEC_ACTUAL_SCAN.data[self.SPEC_SCAN_MOTOR_NAME]
        #print( "Motor values - size: ",self.SPEC_SCAN_MOTOR_DATA.shape
        #print( "plot_scan"
        self.plot_scan()
        #"""
        if not self.IMG_INIT:
            self.init_image()
        gc.collect()
        del threads
        return
    
    def check_and_update_scan_slider(self):
        """Get the actual scan object corresponding to the scan number and image number selected"""
        scan_found = False
        if self.DATA_IS_LOADED and self.SELECTED_IMG_NUM != None:
            self.check_skipped_motors()#To get the list of skipped motors
            #This_Prefix = self.edf_choosen.split("_")[0]
            #prefix = re.compile(This_Prefix)
            for i in range(len(self.SPEC_IMG)):
                if (self.SELECTED_IMG_NUM in self.SPEC_IMG[i]) and (self.SPEC_SCAN_LIST[i].colnames[0].upper() not in self.SPEC_SKIPPED_MOTORS):
                    self.SPEC_ACTUAL_SCAN = self.SPEC_SCAN_LIST[i]
                    scan_found = True
                    break
                else:
                    continue
            
        else:
            self.SPEC_ACTUAL_SCAN = self.SPEC_SCAN_LIST[-1]#if no image is selected, the last scan will be displayed
            
        #print( "actual scan: ",self.SPEC_ACTUAL_SCAN.nr
        if scan_found:    
            self.scan_slider_spinButton.set_value(self.SPEC_ACTUAL_SCAN.nr)#This will call the update scan slider too
            #self.update_scan_slider_now()#NOT NECESSARY
            self.SCAN_ROI_INTEGRATION_X =[]
            self.SCAN_ROI_INTEGRATION_Y =[]
            return
        else:
            return
        
    def slider_plot_scan(self, widget):
        self.plot_scan()
    
    def get_roi_data(self):
        if self.detector_space_btn.isChecked():
            r = sorted([int(self.ROI_y0), int(self.ROI_y1)])
            c = sorted([int(self.ROI_x0), int(self.ROI_x1)])
            
        elif self.tth_chi_space_btn.isChecked():
            x1 = get_index(self.tth_pyFAI, self.ROI_x0)
            x2 = get_index(self.tth_pyFAI, self.ROI_x1)
            y1 = get_index(self.chi_pyFAI, self.ROI_y0)
            y2 = get_index(self.chi_pyFAI, self.ROI_y1)
            r  = sorted([y1,y2])
            c  = sorted([x1,x2])
        elif self.hk_space_btn.isChecked() or self.hl_space_btn.isChecked() or self.kl_space_btn.isChecked() or self.q_space_btn.isChecked():
            x1 = get_index(self.QGridder.xaxis, self.ROI_x0)
            x2 = get_index(self.QGridder.xaxis, self.ROI_x1)
            y1 = get_index(self.QGridder.yaxis, self.ROI_y0)
            y2 = get_index(self.QGridder.yaxis, self.ROI_y1)
            r  = sorted([y1,y2])
            c  = sorted([x1,x2])
        return self.data[r[0]:r[1],c[0]:c[1]].sum()

    def plot_roi(self):
        self.profiles_ax1.cla()
        #self.profiles_ax1.format_coord = self.pro_format_coord
        #self.profiles_ax2.format_coord = self.pro_format_coord
        self.profiles_ax1.plot(self.SCAN_ROI_INTEGRATION_X, self.SCAN_ROI_INTEGRATION_Y, "r-o", lw=2)
        self.profiles_ax1.set_xlabel(self.SPEC_SCAN_MOTOR_NAME, size=14)
        self.profiles_ax1.set_ylabel("ROI integration", size=14)
        self.profiles_refresh()
        self.canvas.draw()
        
    def read_scan_header_D1(self, img_index, scan_motor):
        #self.what_is_this_manip()
        motors = self.SPEC_ACTUAL_SCAN.motors
        self.count_time = self.SPEC_ACTUAL_SCAN.count_time
        
        this_scan_motor_value = self.SPEC_ACTUAL_SCAN.data[scan_motor][img_index]
        mot_del  = re.compile(r'DEL')
        mot_eta  = re.compile(r'ETA')
        mot_chi  = re.compile(r'CHI')
        mot_phi  = re.compile(r'PHI')
        mot_nu  = re.compile(r'NU')
        mot_mu  = re.compile(r'MU')
        mot_list = [mot_del,mot_eta,mot_chi,mot_phi,mot_nu,mot_mu]
        real_mot_list = ['del','eta','chi','phi','nu','mu']
        for i in range(len(mot_list)):
            if mot_list[i].findall(scan_motor.upper())!=[]:
                motors[real_mot_list[i]] = this_scan_motor_value
        self.delta = motors['del']
        self.eta   = motors['eta']
        self.chi   = motors['chi']
        self.phi   = motors['phi']
        self.nu    = motors['nu']
        self.mu    = motors['mu']
    
    def plot_scan(self):
        """ Plot when the scan slider changes """
        if len(self.SPEC_ACTUAL_SCAN_DATA)>0:
            img_num = self.scan_slider_imgSlider.value()
            img_num = int(img_num)
            #print( "Image number: ",img_num
                
            img_index = N.where(self.SPEC_ACTUAL_SCAN_IMG == img_num)
            
            #img_index = N.where(self.current_img_list == img_num)
            img_index = img_index[0][0]
            self.data = self.SPEC_ACTUAL_SCAN_DATA[img_index]
            
            self.header = self.SPEC_ACTUAL_SCAN_HEADER[img_index]
            if self.adj_btn.isChecked():
                if self.data.size % 9600 == 0:
                    self.data = correct_geometry(self.data, detector_type=self.detector_type)
            if self.horizontal_detector:
                self.data = N.rot90(self.data)
            this_title = self.SPEC_ACTUAL_SCAN_IMG_NAMES[img_index]
            scan_motor = self.SPEC_SCAN_MOTOR_NAME
            this_motor_value = self.SPEC_SCAN_MOTOR_DATA[img_index]
            
            this_title = this_title +" - %s = %s"%(scan_motor, this_motor_value)
            self.MAIN_TITLE.set_text(this_title)
            
            if self.detector_type =="D1":
                self.read_scan_header_D1(img_index, scan_motor)
            else:
                self.read_header(self.header)
            if self.detector_space_btn.isChecked():
                self.MAIN_EXTENT = (0, self.data.shape[1], 0, self.data.shape[0])
            elif self.tth_chi_space_btn.isChecked():
                self.Angular_space_plot()
            elif self.hk_space_btn.isChecked():
                self.Reciprocal_space_plot(space="HK")
                #print( "X shape: ",self.QGridder.xaxis.shape
            elif self.hl_space_btn.isChecked():
                self.Reciprocal_space_plot(space="HL")
            elif self.kl_space_btn.isChecked():
                self.Reciprocal_space_plot(space="KL")
            elif self.q_space_btn.isChecked():
                self.Reciprocal_space_plot(space="Q")
            #print( "Data shape: ",self.data.shape
            if self.ROI_ON and self.ROI_DRAWN:
                self.SCAN_ROI_INTEGRATION_X.append(this_motor_value)
                roi_data = self.get_roi_data()
                self.SCAN_ROI_INTEGRATION_Y.append(roi_data)
                self.plot_roi()
            self.scale_plot()
            if img_num == self.SPEC_ACTUAL_SCAN_IMG.min() or img_num == self.SPEC_ACTUAL_SCAN_IMG.max():
                self.SCAN_ROI_INTEGRATION_X = []
                self.SCAN_ROI_INTEGRATION_Y = []
            
            self.img.set_extent(self.MAIN_EXTENT)
            self.slider_update()
                        
            self.SELECTED_IMG_NUM = img_num
            if isfile(join(self.edf_folder, self.SPEC_ACTUAL_SCAN_IMG_NAMES[img_index])):
                self.fabioIMG = fabio.open(join(self.edf_folder, self.SPEC_ACTUAL_SCAN_IMG_NAMES[img_index]))
            else:
                pass
        gc.collect()
        return
        
    def plot_3D_scan(self,w):
        import matplotlib.pyplot as plt
        """ popup a mayavi window to visualize the 3D data """
        gc.collect()
        if self.use_spec_UB.isChecked():
            self.UB_MATRIX = self.get_UB_from_spec(self.SPEC_ACTUAL_SCAN)
            self.qconv.UB = self.UB_MATRIX
            N.savetxt("UB_scan_%d.dat"%self.SPEC_ACTUAL_SCAN.nr, self.UB_MATRIX)
        #self.what_is_this_manip()
        #print( "Manip ",self.manip
        B = UB2B(self.UB_MATRIX)
        U = self.UB_MATRIX.dot(np.linalg.inv(B))
        DATA = []
        scan_motors = {}
        cch1 = self.direct_beam[1]
        cch2 = self.direct_beam[0]
        dim = self.SPEC_ACTUAL_SCAN_DATA[0].shape
        Nch1 = dim[0]
        Nch2 = dim[1]
        # reduce data: number of pixels to average in each detector direction
        default_nav = [4,4]
        default_roi = [0,Nch1, 0,Nch2]  # region of interest on the detector
        contour_level = 50
        if self.detector_type == "S70" or "MAXIPIX" in _DETECTOR_NAME.upper():
            contour_level = 50
            default_nav = [1,1]
        if self.detector_space_btn.isChecked() and self.ROI_ON and self.ROI_DRAWN:
            r = sorted([int(self.ROI_y0), int(self.ROI_y1)])
            c = sorted([int(self.ROI_x0), int(self.ROI_x1)])
            default_roi = [r[0],r[1],c[0],c[1]]
            
        #default_roi = np.asarray(default_roi)
        #default_nav = np.asarray(default_nav)    
        has_det_rot_out = False
        det_rot_in = _DETECTOR_MOTORS[-1]
        if _DETECTOR_ROTATIONS ==2:
            det_rot_out = _DETECTOR_MOTORS[0]
            has_det_rot_out = True
        sam_rot1 = _SAMPLE_ROTATION_MOTORS[0]
        sam_rot2 = _SAMPLE_ROTATION_MOTORS[1]
        sam_rot3 = _SAMPLE_ROTATION_MOTORS[2]
        scan_motors[det_rot_in] = []
        if has_det_rot_out:
            scan_motors[det_rot_out] = []
        scan_motors[sam_rot1]     = []
        scan_motors[sam_rot2]     = []
        scan_motors[sam_rot3]     = []

        print( "samp_rot_1, 2, 3, tth:")
        print( sam_rot1, sam_rot2, sam_rot3, det_rot_in)
        print( "nbr of images: ",len(self.SPEC_ACTUAL_SCAN_DATA))
        print( "Data processing ...")
    
        for i in range(len(self.SPEC_ACTUAL_SCAN_DATA)):
            # print( "Loading image: ",self.SPEC_ACTUAL_SCAN_IMG_NAMES[i]
            motor=get_motors(self.SPEC_ACTUAL_SCAN_HEADER[i])            
            data = self.SPEC_ACTUAL_SCAN_DATA[i]
            scan_motors[det_rot_in].append(motor[det_rot_in])
            if has_det_rot_out:
                scan_motors[det_rot_out].append(motor[det_rot_out])
            scan_motors[sam_rot1].append(motor[sam_rot1])
            scan_motors[sam_rot2].append(motor[sam_rot2])
            scan_motors[sam_rot3].append(motor[sam_rot3])
            data = xrayutilities.blockAverage2D(data, default_nav[0], default_nav[1], roi=default_roi)
            DATA.append(data)        
        
        
        try:
            tth   = N.asarray(scan_motors[det_rot_in])
            th   = N.asarray(scan_motors[sam_rot1])
            chi  = N.asarray(scan_motors[sam_rot2])
            phi  = N.asarray(scan_motors[sam_rot3])
            DATA = N.asarray(DATA)
            
            this_experiment = xrayutilities.HXRD([1,0,0],[0,0,1], en=self.energy, qconv=self.qconv)
            this_experiment.Ang2Q.init_area(_DETECTOR_PIXEL_DIRECTIONS[0], _DETECTOR_PIXEL_DIRECTIONS[1], cch1=cch1, cch2=cch2, Nch1=Nch1,Nch2=Nch2, 
                                      pwidth1=_PIXEL_SIZE,pwidth2=_PIXEL_SIZE, distance=self.distance, 
                                     # tiltazimuth=self.xu_tiltazimuth, tilt=self.xu_tilt, detrot=self.xu_detrot, outerangle_offset=self.xu_outerangle_offset,
                                      Nav=default_nav, roi=default_roi)
            
            """
            if has_det_rot_out:
                h,k,l=this_experiment.Ang2HKL(th,chi,phi,nu,tth, dettype='area')
            else:
                h,k,l=this_experiment.Ang2HKL(th,chi,phi,tth, dettype='area')
            """
            if self.has_substrate:
                h,k,l=this_experiment.Ang2HKL(th,chi,phi,tth, dettype='area', mat=self.substrate)
            else:
                #print( "Should be here"            
                h,k,l=this_experiment.Ang2HKL(th,chi,phi, tth, dettype='area', U=U, B=B)
            
            nx = 101
            ny = 101
            nz = 101
            print( "Initializing data grid...")
            gridder = xrayutilities.Gridder3D(nx,ny,nz)
            print( "Data gridding on a regular map...")
            gridder(h,k,l,DATA)
            
            MAP_data  = gridder.data.T
            maxint    = N.log10(MAP_data.max())
            MAP_data  = xrayutilities.maplog(MAP_data,maxint*0.5,0)
            print( "Plotting 3D image ...")
            x,y,z = gridder.xaxis, gridder.yaxis, gridder.zaxis
            
            PopUpRSM(MAP_data, gridder.xaxis, gridder.yaxis, gridder.zaxis)
        except:
            exc_type, exc_value, exc_traceback = sys.exc_info()
            self.popup_info("warning", "ERROR: %s"%str(exc_value))
        
        gc.collect()
        return

    
    def get_UB_from_spec(self,actual_scan):
        header = actual_scan.header
        ub     = []
        for h in header:
            if h.startswith("#G3"):
                h = h.split()
                h = h[1:]
                for i in range(len(h)):
                    ub.append(float(h[i]))
                ub = N.asarray(ub)
                ub = ub.reshape((3,3))
                return ub
            else:
                continue
            
    def change_scale(self,button, data):
        if button.isChecked():
            button.setText("Linear scale")
            data = N.log10(data+1e-6)
            
            actual_vmin = self.sld_vmin.value()
            actual_vmax = self.sld_vmax.value()
            #print( data
            #print( "Log scale called. Actual VMAX: ",actual_vmax
            self.vmax = N.log10(actual_vmax) if self.log_scale == 0 else actual_vmax
            if actual_vmin == 0:
                self.vmin=0
            elif actual_vmin >0:
                self.vmin = N.log10(actual_vmin) if self.log_scale == 0 else actual_vmin

            #vmax_range = N.max(data[N.isfinite(data)])
            self.vmax_range = data.max()
            #self.vmax_range = self.vmax_range if self.vmax_range < 8*med_log else 8*med_log
            self.log_scale = 1

        else:
            button.setText("Log scale")
            actual_vmin = self.sld_vmin.value()
            actual_vmax = self.sld_vmax.value()
            #print( "Linear scale called, Actual vmax: ",actual_vmax
            if self.log_scale == 0:
                self.vmax = actual_vmax
            elif self.log_scale == 1:
                self.vmax = N.power(10.,actual_vmax)
            #print( data
            self.vmax_range = data.max()
            #vmax_range = N.max(data[N.isfinite(data)])
            #if self.vmax_range > 100*med_lin:
            #       self.vmax_range = 100*med_lin

            if actual_vmin ==0:
                self.vmin = 0
            elif actual_vmin>0:
                if self.log_scale == 0:
                    self.vmin = actual_vmin
                elif self.log_scale == 1:
                    self.vmin = N.power(10.,actual_vmin)
            self.log_scale = 0

        self.sld_vmax.setRange(0,self.vmax_range)
        self.sld_vmin.setRange(0,self.vmax_range)
        return data
    
    def scale_plot(self):
        data = self.data.copy()
        data = self.change_scale(self.linear_scale_btn, data)
        self.img.set_data(data)
        #self.img.set_extent(self.MAIN_EXTENT)
        #self.canvas.draw()
    
    def log_update(self,widget):
        self.scale_plot()
        if self.log_scale==1:
            self.cb.set_label(r'$Log_{10}\ (Counts\ per\ second)\ [arb.\ units]$',fontsize=12)
            self.cb2.set_label(r'$Log_{10}\ (Counts\ per\ second)\ [arb.\ units]$',fontsize=12)
        else:
            self.cb.set_label(r'$Intensity\ (Counts\ per\ second)$', fontsize=12)
            self.cb2.set_label(r'$Intensity\ (Counts\ per\ second)$', fontsize=12)
        self.slider_update()

    def scale_update(self,widget):
        global g_is_image_loaded
        if g_is_image_loaded==False:
            #self.popup_info("Info","No image has been loaded. Kindly load an image")
            return
        self.vmin = self.sld_vmin.value()
        self.vmax = self.sld_vmax.value()
        self.vmin_spin_btn.setValue(self.vmin)
        self.vmax_spin_btn.setValue(self.vmax)
        self.slider_update()

    def scale_update_spin(self,widget):
        global g_is_image_loaded
        if g_is_image_loaded==False:
            #self.popup_info("Info","No image has been loaded. Kindly load an image")
            return
        self.vmin = self.vmin_spin_btn.value()
        self.vmax = self.vmax_spin_btn.value()
        self.slider_update()

    def slider_update(self):
        global g_is_image_loaded
        if g_is_image_loaded==False:
            self.popup_info("Info","No image has been loaded. Kindly load an image")
            return
        if self.notebook.currentIndex() == 0:
            self.img.set_clim(self.vmin, self.vmax)
            self.sld_vmax.setValue(self.vmax)
            self.sld_vmin.setValue(self.vmin)
            if self.linear_scale_btn.isChecked():
                #self.vmin_spin_btn.set_adjustment(gtk.Adjustment(self.vmin, 0, self.vmax_range, 0.1, 1.0, 0))
                self.vmin_spin_btn.setValue(self.vmin)
                self.vmin_spin_btn.setMinimum(0)
                self.vmin_spin_btn.setMaximum(self.vmax_range) 
                self.vmin_spin_btn.setSingleStep(0.1)
                self.vmax_spin_btn.setValue(self.vmax)
                self.vmax_spin_btn.setMinimum(0)
                self.vmax_spin_btn.setMaximum(self.vmax_range)
                self.vmax_spin_btn.setSingleStep(0.1)
                #self.vmax_spin_btn.set_adjustment(gtk.Adjustment(self.vmax, 0, self.vmax_range, 0.1, 1.0, 0))
            else:
                #self.vmin_spin_btn.set_adjustment(gtk.Adjustment(self.vmin, 0, self.vmax_range, 10, 100, 0))
                #self.vmax_spin_btn.set_adjustment(gtk.Adjustment(self.vmax, 0, self.vmax_range, 10, 100, 0))
                self.vmin_spin_btn.setValue(self.vmin)
                self.vmin_spin_btn.setMinimum(0)
                self.vmin_spin_btn.setMaximum(self.vmax_range) 
                self.vmin_spin_btn.setSingleStep(10)
                self.vmax_spin_btn.setValue(self.vmax)
                self.vmax_spin_btn.setMinimum(0)
                self.vmax_spin_btn.setMaximum(self.vmax_range)
                self.vmax_spin_btn.setSingleStep(10)
            #self.vmax_spin_btn.update()
            
            #self.ax.relim()
            #self.img.set_extent(self.MAIN_EXTENT)#******* RECHECK THIS FOR 2THETA-CHI
            if self.IMG_ZOOMED == True:
                self.ax.set_xlim(self.ZOOM_EXTENT[0], self.ZOOM_EXTENT[1])
                self.ax.set_ylim(self.ZOOM_EXTENT[2], self.ZOOM_EXTENT[3])
            else:
                self.ax.set_xlim(self.MAIN_EXTENT[0], self.MAIN_EXTENT[1])
                self.ax.set_ylim(self.MAIN_EXTENT[2], self.MAIN_EXTENT[3])
            self.ax.relim()
            self.canvas.draw()
            #self.ax.figure.canvas.draw_idle()
        elif self.notebook.currentIndex() == 1:
            self.polar_img.set_clim(self.vmin, self.vmax)
            self.sld_vmax.setValue(self.vmax)
            self.sld_vmin.setValue(self.vmin)
            #self.vmin_spin_btn.set_adjustment(gtk.Adjustment(self.vmin, 0, self.vmax_range, 0.5, 10.0, 0))
            #self.vmax_spin_btn.set_adjustment(gtk.Adjustment(self.vmax, 0, self.vmax_range, 0.5, 10.0, 0))
            self.vmin_spin_btn.setValue(self.vmin)
            self.vmin_spin_btn.setMinimum(0)
            self.vmin_spin_btn.setMaximum(self.vmax_range) 
            self.vmin_spin_btn.setSingleStep(0.5)
            self.vmax_spin_btn.setValue(self.vmax)
            self.vmax_spin_btn.setMinimum(0)
            self.vmax_spin_btn.setMaximum(self.vmax_range)
            self.vmax_spin_btn.setSingleStep(0.5)
            #self.vmax_spin_btn.update()
            #self.polar_ax.relim()
            self.plot_PF()

    def detector_disposition(self,widget):
        """ Set detector vertically or horizontally """
        #data = self.data.copy()
        # if self.detector_disposition_horizontal.get_active():
            # self.horizontal_detector = True
            # self.detector_disposition_horizontal.set_label("Rotate back")
            # self.data = N.rot90(self.data)
        # else:
            # self.horizontal_detector = False
            # self.detector_disposition_horizontal.set_label("Rotate detector")
            # self.data = N.rot90(self.data,3)
        self.rotate_detector_n +=1
        if self.rotate_detector_n == 3:
            self.rotate_detector_n = 0
        
        self.data = N.rot90(self.data,self.rotate_detector_n)
        self.img.set_array(self.data)

        imshape = self.data.shape
        self.xDim1 = imshape[1]
        self.yDim1 = imshape[0]
        self.MAIN_EXTENT =(0, self.xDim1, 0,self.yDim1)
        if self.xDim1 > self.yDim1:
            self.cb.ax.set_visible(False)
            self.cb2.ax.set_visible(True)
        else:
            self.cb.ax.set_visible(True)
            self.cb2.ax.set_visible(False)
        self.img.set_extent(self.MAIN_EXTENT)
        # self.ax.set_xlim(0,self.xDim1)
        # self.ax.set_ylim(0,self.yDim1)
        self.slider_update()

    def calcul_chi_2theta_d(self,event):
        """ EVENT Calculate chi, 2theta of the corresponding points on the image """
        x = event.xdata
        y = event.ydata
        self.check_azimuthal_integrator()
        chi = self.tableChi[y,x]
        tth = self.tableTwoTheta[y,x]
        d   = self.table_dSpace[y,x]
        return chi, tth, d

    def show_chi_delta(self,widget):
        """ If the show_chi_delta_btn is checked and the calibration file is loaded """
        if self.show_chi_delta_btn.isChecked():
            if self.calibrated and self.geometry_corrected:
                self.show_chi_delta_flag = True
                self.show_chi_txt.show()#set_visible(True)
                self.show_delta_txt.show()#set_visible(True)
                self.show_d_txt.show()#set_visible(True)
            else:
                self.show_chi_delta_btn.setChecked(False)
                self.show_chi_delta_flag = False
                self.popup_info("warning","Please calibrate the detector and correct the detector's geometry before checking this box!")

        else:
            self.show_chi_delta_flag = False
            self.show_chi_txt.hide()#set_visible(False)
            self.show_delta_txt.hide()#set_visible(False)
            self.show_d_txt.hide()#set_visible(False)

    def status_update(self,event):
        if event.inaxes==self.ax:
            xdata = event.xdata
            ydata = event.ydata
            if self.detector_space_btn.isChecked():
                zdata = self.data[int(ydata),int(xdata)]
            elif self.tth_chi_space_btn.isChecked():
                x = get_index(self.tth_pyFAI, xdata)
                y = get_index(self.chi_pyFAI, ydata)
                zdata = self.data[y,x]
            elif self.hk_space_btn.isChecked() or self.hl_space_btn.isChecked() or self.kl_space_btn.isChecked() or self.q_space_btn.isChecked():
                x = get_index(self.QGridder.xaxis, xdata)
                y = get_index(self.QGridder.yaxis, ydata)
                #print( "x=%d y=%d"%(x,y)
                zdata = self.data[y,x]
            if self.detector_space_btn.isChecked():
                self.x_pos.setText("%d"%xdata)
                self.y_pos.setText("%d"%ydata)
            else:
                self.x_pos.setText("%.2f"%xdata)
                self.y_pos.setText("%.2f"%ydata)
            self.z_pos.setText("%d"%zdata)
            
            if self.show_chi_delta_flag == True:
                if self.detector_space_btn.isChecked():
                    chi,tth,d = self.calcul_chi_2theta_d(event)
                    self.show_chi_txt.setText("Chi = %.2f"%chi)
                    self.show_delta_txt.setText("2Theta = %.2f"%tth)
                    self.show_d_txt.setText("d = %.4f A"%d)
                elif self.tth_chi_space_btn.isChecked():
                    tth = event.xdata
                    chi = event.ydata
                    d   = self.azimuthalIntegration.wavelength / (2*N.sin(N.radians(tth/2.0))) * 1e10 #d in Angstrom
                    self.show_chi_txt.setText("Chi = %.2f"%chi)
                    self.show_delta_txt.setText("2Theta = %.2f"%tth)
                    self.show_d_txt.setText("d = %.4f A"%d)
                
    def zoom_on(self,widget):
        """For the Zoom button"""
        global g_zoom_tb
        if g_zoom_tb==0 or self.zoomtb.isChecked():
            self.zoom = True
            self.popup_info("info","The activity is set for zooming. Kindly zoom the image")
            #self.canvas.window.setShape(Qt.CrossCursor)  #gtk.gdk.Cursor(gtk.gdk.CROSS))
            #self.canvas.window.set_cursor(QCursor(Qt.CrossCursor))  #gtk.gdk.Cursor(gtk.gdk.CROSS))
            self.cursor.visible = False
            g_zoom_tb = 1
        else:
            self.zoom = False
            #self.canvas.window.set_cursor(None)
            #self.canvas.window.setShape(Qt.ArrowCursor)#set_cursor(None)
            self.cursor.visible = True
            g_zoom_tb = 0

    def zoom_image(self):
        tmp_x = [self.x0, self.x1]
        tmp_y = [self.y0, self.y1]
        tmp_x = N.sort(tmp_x)
        tmp_y = N.sort(tmp_y)
        #self.xDim0, self.xDim1 = tmp_x[0], tmp_x[1]
        #self.yDim0, self.yDim1 = tmp_y[0], tmp_y[1]
        extent_x0, extent_x1 = tmp_x[0],tmp_x[1]
        #self.ax.set_xlim(self.xDim0, self.xDim1)
        self.ax.set_xlim(tmp_x[0], tmp_x[1])
        #if self.detector_type not in ["D5", "D1"]:
            ##self.ax.set_ylim(self.yDim1, self.yDim0)
            #self.ax.set_ylim(tmp_y[1], tmp_y[0])
            #extent_y0, extent_y1 = tmp_y[1],tmp_y[0]
        #else:
            ##self.ax.set_ylim(self.yDim0, self.yDim1)
            #self.ax.set_ylim(tmp_y[0], tmp_y[1])
            #extent_y0, extent_y1 = tmp_y[0],tmp_y[1]
        self.ax.set_ylim(tmp_y[0], tmp_y[1])
        extent_y0, extent_y1 = tmp_y[0],tmp_y[1]
        self.cb.ax.set_visible(False)
        self.MAIN_EXTENT = (extent_x0, extent_x1, extent_y0, extent_y1)
        #self.img.set_extent(self.MAIN_EXTENT)
        #self.ax.set_xlim(self.MAIN_EXTENT[0], self.MAIN_EXTENT[1])
        #self.ax.set_ylim(self.MAIN_EXTENT[2], self.MAIN_EXTENT[3])
        self.IMG_ZOOMED = True
        self.ZOOM_EXTENT = (extent_x0, extent_x1, extent_y0, extent_y1)
        self.canvas.draw()
        
    def reset_scale(self,widget):
        try:
            if self.linear_scale_btn.isChecked():
                self.vmin = 0
                self.vmax = N.log10(10*self.med)
            else:
                self.vmin = 0
                self.vmax = 10*self.med
            self.slider_update()
        except:
            self.popup_info("Error", "Error has occured may be due to no image being plotted or .edf file is loaded yet")
            return

    def reset_image(self,widget):
        self.Do_reset_image()
        
    def Do_reset_image(self):
        """For the Home button"""
        self.IMG_ZOOMED = False
        self.xDim0 = 0
        self.xDim1 = self.data.shape[1]
        self.yDim0 = 0
        self.yDim1 = self.data.shape[0]
        if self.tth_chi_space_btn.isChecked():
            self.xDim0 = self.tth_pyFAI.min()
            self.xDim1 = self.tth_pyFAI.max()
            self.yDim0 = self.chi_pyFAI.min()
            self.yDim1 = self.chi_pyFAI.max()
        if self.hk_space_btn.isChecked() or self.hl_space_btn.isChecked() or self.kl_space_btn.isChecked() or self.q_space_btn.isChecked():
            self.xDim0 = self.QGridder.xaxis.min()
            self.xDim1 = self.QGridder.xaxis.max()
            self.yDim0 = self.QGridder.yaxis.min()
            self.yDim1 = self.QGridder.yaxis.max()
            
        self.MAIN_EXTENT = (self.xDim0,self.xDim1,self.yDim0, self.yDim1)
        #if self.detector_type=="S70":
            #self.MAIN_EXTENT = (self.xDim0,self.xDim1,self.yDim1, self.yDim0)
        
        self.ax.set_xlim(self.MAIN_EXTENT[0],self.MAIN_EXTENT[1])
        self.ax.set_ylim(self.MAIN_EXTENT[2],self.MAIN_EXTENT[3])
        self.img.set_extent(self.MAIN_EXTENT)
        if self.data.shape[0] < self.data.shape[1]:
            self.cb.ax.set_visible(False)
            self.cb2.ax.set_visible(True)
        else:
            self.cb.ax.set_visible(True)
            self.cb2.ax.set_visible(False)

        self.canvas.draw()

    
    def get_list_dir(self,this_dir):
        ''' Get the list of directories inside this_dir'''
        #self.MODEL.clear()
        main_dir = this_dir
        list_dir = listdir(main_dir)
        no_of_threads = len(list_dir)
        pool = Pool()
        for i in range(no_of_threads):
            #t = threading.Thread(target=self._thread_scanning,args=(main_dir,list_dir[i].decode('utf8'),))
            #self.threads.append(t)
            #t.start()
            # self._thread_scanning(main_dir,list_dir[i].decode('utf8'))
            pool.apply_async(self._thread_scanning,args=(main_dir,list_dir[i],))
        pool.close()
        pool.join()
    
    def add_parent(self,parent, name):
        '''parent: a model parent, name: text '''
        return self.MODEL.append(parent, name)

    def add_child(self, parent_item, child_name):
        ''' parent_item: parent returned by self.add_parent'''
        self.MODEL.append(parent_item, child_name)
    """
    def get_list_dir(self, this_dir):
        self.MODEL.clear()
        for main_dir, subdir, files in os.walk(this_dir):
            if main_dir == self.current_folder:
                parent = self.add_parent(None, None)
            else:
                md = os.path.basename(main_dir)
                parent = self.add_parent(parent,md)
            main_store = [i for i in files if i.endswith(_IMG_FORMAT) or i.endswith(_IMG_FORMAT +".gz")]
            main_store = list_to_table(main_store,sort_col=0)
            if len(main_store)>0:
                self.TABLE_STORE[str(main_dir)] = main_store
                self.store[str(main_dir)] = get_column_from_table(main_store,0)
                for f in main_store:
                    self.add_child(parent,[f[0]])            
    """        
    def _thread_scanning(self,main_d,list_d):
        path = os.sep.join((main_d, list_d))  # Made use of os's sep instead...
        path = path.decode('utf8')
        if os.path.isdir(path):
            main_store= [i for i in listdir(path) if isfile(join(path,i)) and i.endswith(_IMG_FORMAT) or i.endswith(_IMG_FORMAT +".gz")]
            main_store = list_to_table(main_store,sort_col=0)
            if len(main_store)>0:
                parent = self.MODEL.append(None,[list_d])
                self.TABLE_STORE[str(path)] = main_store
                self.store[str(path)] = get_column_from_table(main_store,0)
                for f in main_store:
                    self.MODEL.append(parent,[f[0]])    
    
    def choose_folder(self,w):
        #dialog = gtk.FileChooserDialog(title="Select an EDF folder",action=gtk.FILE_CHOOSER_ACTION_SELECT_FOLDER, buttons = (gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL, gtk.STOCK_OPEN, gtk.RESPONSE_OK))
        dialog = QFileDialog(self, "Select an EDF folder", self.current_folder)
        dialog.setFileMode(QFileDialog.Directory)
        #dialog.setOption(QtGui.QFileDialog.ShowDirsOnly, True)
        #dialog.set_current_folder(self.current_folder)
        #response=dialog.run()

        if dialog.exec_()==QDialog.Accepted: #response==gtk.RESPONSE_OK:
            folder=dialog.selectedFiles()[0] #getOpenFileName()  #get_filename()
            #folder=folder.decode('utf8')
            folder_basename = os.path.basename(os.path.dirname(folder))
            #print( folder
            main_store= [i for i in listdir(folder) if isfile(join(folder,i)) and i.endswith(_IMG_FORMAT) or i.endswith(_IMG_FORMAT+".gz")]
            self.store = {}       #{'folder_1':[list name], 'folder_2': [list name], ...}
            self.TABLE_STORE = {} #{'folder_1":[table: name-prefix-number], 'folder_2':[table:name-prefix-number],...}
            
            self.current_folder = folder
            #print( self.store
            self.threads = []
            self.get_list_dir(self.current_folder)
            for t in self.threads:
                t.join()
            while len(self.threads)>0:
                self.threads.pop()
            #gobject.timeout_add(100, self._callback)
            if len(main_store)>0:
                #main_store = sorted(main_store)
                main_store = list_to_table(main_store,sort_col=0)
                self.TABLE_STORE[str(folder)] = main_store
                self.store[str(folder)] = get_column_from_table(main_store,0)
                filter1 = ["EDF (*.edf)","EDF GZ(*.edf.gz)"]
                self.MODEL.setNameFilters( filter1 )
                self.MODEL.setRootPath( folder_basename )

                #for i in main_store:
                #    self.MODEL.append(None,[i[0]])
            else:
                pass
            #self.TVcolumn.set_title(folder_basename)
            self.treeView.setWindowTitle(folder_basename)
            self.DATA_IS_LOADED = True
            
        else:
            pass
        #dialog.destroy()
        self.threads=[]
        if self.DATA_IS_LOADED:
            self.store_img = {}#store_img = {'folder':{'prefix': [list img_num], 'prefix2':[list inum],...} }
            for k in self.store.keys():
                self.store_img[k] = group_image_from_table(self.TABLE_STORE[k])

        self.folder_update('')

    def folder_update(self,widget):
        folder = self.current_folder
        if folder is not os.getcwd():
            main_store= [i for i in listdir(folder) if isfile(join(folder,i)) and i.endswith(_IMG_FORMAT) or i.endswith(_IMG_FORMAT+".gz")]
            main_store = list_to_table(main_store,sort_col=0)
            self.store={}
            self.TABLE_STORE = {}
            #self.list_store.clear()
            self.TABLE_STORE[self.current_folder] = main_store
            self.store[self.current_folder] = get_column_from_table(main_store,0)
            self.threads = []
            self.get_list_dir(self.current_folder)
            for t in self.threads:
                t.join()
            while len(self.threads)>0:
                self.threads.pop()
            #gobject.timeout_add(100, self._callback)
            self.store_img = {}
            for k in self.store.keys():
                self.store_img[k] = group_image_from_table(self.TABLE_STORE[k])
            
            #for i in main_store:
            #    self.MODEL.append(None,[i[0]])

            self.treeView.setRootIndex( self.MODEL.setRootPath( folder ) )
            
            self.DATA_IS_LOADED = True
            self.threads =[]
            #if self.SPEC_DATA:
                #self.SPEC_DATA.Update()
                #self.update_spec_data()
        #return 1

    def save_image(self,widget):
        #dialog = gtk.FileChooserDialog(title="Save image", action=gtk.FILE_CHOOSER_ACTION_SAVE, buttons = (gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL, gtk.STOCK_SAVE, gtk.RESPONSE_OK))
        dialog = QFileDialog.getSaveFileName(self, "Save Image", self.current_folder, "Images (*.png *.jpg *.pdf *.ps *.eps);; Other (*)")#, self.current_folder)
        #dialog.setFileMode(QFileDialog.Directory)
        filename = self.edf_choosen.split(".")[0] if self.edf_choosen != "" else ""
        #dialog.set_current_name(filename+".pdf")
        #dialog.set_filename(filename)
        #dialog.set_current_folder(self.current_folder)
        #filtre = gtk.FileFilter()
        #filtre.set_name("images")
        #filtre.add_pattern("*.png")
        #filtre.add_pattern("*.jpg")
        #filtre.add_pattern("*.pdf")
        #filtre.add_pattern("*.ps")
        #filtre.add_pattern("*.eps")
        #dialog.setNameFilters(["Images (*.png *.jpg *.pdf *.ps *.eps)", "Other (*)"])
        #filtre = gtk.FileFilter()
        #filtre.set_name("Other")
        #filtre.add_pattern("*")
        #dialog.add_filter(filtre)
        #response = dialog.run()
        print("File name for saving is ",dialog)
        if dialog!="" and type(dialog)==type((1,2)):
            dialog = dialog[0]
            if len(dialog)<1:
                print("Cannot save the data. Provide a valid file name")
            print("File is saved as ",dialog)
            self.fig.savefig( dialog )
        #if response==gtk.RESPONSE_OK:
        #    self.fig.savefig(dialog.get_filename().decode('utf8'))
        #dialog.destroy()

    def save_adjust(self,widget):
        """ Save the current EDF image, the file name will be name+adjusted """
        global g_is_image_loaded
        if g_is_image_loaded==False:
            self.popup_info("Info","No image has been loaded. Kindly load an image the EDF Image")
            return
        self.fabioIMG.data = self.data
        basename = os.path.basename(self.edf)
        name = basename.split(".")[0]+"_corrected.edf"
        filename = join(os.path.dirname(self.edf), name)
        #filename = filename.decode('utf8')
        self.fabioIMG.write(filename)
        self.popup_info("info", "Image %s is successfully saved !"%filename)

    def plot_update(self,widget):
        if self.tth_chi_space_btn.isChecked():
            self.MAIN_XLABEL.set_text(r"$2\theta\ (deg.)$")
            self.MAIN_YLABEL.set_text(r"$\chi\ (deg.)$")
        
        elif self.hk_space_btn.isChecked():
            self.MAIN_XLABEL.set_text(r"$H\ (R.L.U)$")
            self.MAIN_YLABEL.set_text(r"$K\ (R.L.U)$")
        elif self.hl_space_btn.isChecked():
            self.MAIN_XLABEL.set_text(r"$H\ (R.L.U)$")
            self.MAIN_YLABEL.set_text(r"$L\ (R.L.U)$")
        elif self.kl_space_btn.isChecked():
            self.MAIN_XLABEL.set_text(r"$K\ (R.L.U)$")
            self.MAIN_YLABEL.set_text(r"$L\ (R.L.U)$")
        elif self.q_space_btn.isChecked():
            self.MAIN_XLABEL.set_text(r"$Qy\ (nm^{-1})$")
            self.MAIN_YLABEL.set_text(r"$Qz\ (nm^{-1})$")
            
        else:
            self.MAIN_XLABEL.set_text("X (pixel)")
            self.MAIN_YLABEL.set_text("Y (pixel)")
        #self.Do_reset_image()
        self.plot_data()
    
    def Angular_space_plot(self):
        if self.calibrated == False:
            self.popup_info('warning','Please calibrate the detector before checking this!')
            self.tth_chi_space_btn.setChecked(False)
        elif self.calibrated==True and self.geometry_corrected==True:
            
            self.show_chi_delta_btn.setEnabled(False)
            self.show_chi_delta_flag=self.show_chi_delta_btn.isChecked()
            self.check_azimuthal_integrator()
            
            self.azimuthalIntegration.setChiDiscAtZero()
            
            self.data,self.tth_pyFAI,self.chi_pyFAI = self.azimuthalIntegration.integrate2d(self.data,self.data.shape[0],self.data.shape[1],unit="2th_deg")
            
            self.chi_pyFAI = self.chi_pyFAI - 90
            self.MAIN_EXTENT = (self.tth_pyFAI.min(), self.tth_pyFAI.max(), self.chi_pyFAI.min(), self.chi_pyFAI.max())
            
        else:
            self.popup_info('warning',"Please correct the detector's geometry prior to proceed this operation!")
            self.tth_chi_space_btn.setChecked(False)

    def Reciprocal_space_plot(self,space="HK"):
        """ space should be HK, HL or KL """
        if self.calibrated == False:
            self.popup_info('warning','Please calibrate the detector before checking this!')
            self.hk_space_btn.setChecked(False)
            self.hl_space_btn.setChecked(False)
            self.kl_space_btn.setChecked(False)
            self.q_space_btn.setChecked(False)
            self.tth_chi_space_btn.setChecked(False)
        elif self.calibrated==True and self.geometry_corrected==True:
            self.show_chi_delta_btn.set_sensitive(False)
            self.show_chi_delta_flag=False
            self.check_azimuthal_integrator()
            distance = self.azimuthalIntegration.dist
            #print( "Distance: ",distance
            cch1     = self.direct_beam[1]
            cch2     = self.direct_beam[0]
            detrot   = self.azimuthalIntegration.rot3
            tiltazimuth=self.azimuthalIntegration.rot1
            tilt     = self.azimuthalIntegration.rot2
            if detrot != None:
                detrot = N.degrees(detrot)*(-1)
            else:
                detrot = 0
            if tiltazimuth != None:
                tiltazimuth = N.degrees(tiltazimuth)*(-1)
            else:
                tiltazimuth = 0
            if tilt !=None:
                tilt   = N.degrees(tilt)*(-1)
            else:
                tilt   = 0
            ETA = self.th
            CHI = self.chi
            PHI = self.phi 
            if hasattr(self,'det_outer'): NU  = self.det_outer
            DEL = self.det_inner
            
            Nch1 = self.data.shape[0]
            Nch2 = self.data.shape[1]
            
            #pixel binning to reduce memory consumption
            # if self.detector_type != "S70":
            dim1 = self.data.shape[0]/3
            dim2 = self.data.shape[1]/3

                
            self.experiment.Ang2Q.init_area(_DETECTOR_PIXEL_DIRECTIONS[0], _DETECTOR_PIXEL_DIRECTIONS[1], cch1=cch1, cch2=cch2, Nch1=Nch1,Nch2=Nch2, pwidth1=_PIXEL_SIZE,pwidth2=_PIXEL_SIZE, distance=distance, detrot=detrot)
            if hasattr(self,'det_outer'):
                self.H,self.K,self.L = self.experiment.Ang2HKL(ETA,CHI,PHI,NU,DEL,dettype='area')
                self.QX,self.QY,self.QZ = self.experiment.Ang2Q.area(ETA,CHI,PHI,NU,DEL)
            else:
                self.H,self.K,self.L = self.experiment.Ang2HKL(ETA,CHI,PHI,DEL,dettype='area')
                self.QX,self.QY,self.QZ = self.experiment.Ang2Q.area(ETA,CHI,PHI,DEL)
            if self.has_substrate:
                if hasattr(self,'det_outer'):
                    self.H,self.K,self.L = self.experiment.Ang2HKL(ETA,CHI,PHI,NU,DEL, mat=self.substrate, dettype='area')
                    self.QX,self.QY,self.QZ = self.experiment.Ang2Q.area(ETA,CHI,PHI,NU,DEL)
                else:
                    self.H,self.K,self.L = self.experiment.Ang2HKL(ETA,CHI,PHI,DEL, mat=self.substrate, dettype='area')
                    self.QX,self.QY,self.QZ = self.experiment.Ang2Q.area(ETA,CHI,PHI,DEL)
            else:
                if hasattr(self,'det_outer'):
                    self.H,self.K,self.L = self.experiment.Ang2HKL(ETA,CHI,PHI,NU,DEL, dettype='area')
                    self.QX,self.QY,self.QZ = self.experiment.Ang2Q.area(ETA,CHI,PHI,NU,DEL)
                else:
                    self.H,self.K,self.L = self.experiment.Ang2HKL(ETA,CHI,PHI,DEL, dettype='area')
                    self.QX,self.QY,self.QZ = self.experiment.Ang2Q.area(ETA,CHI,PHI,DEL)
            
            self.QGridder = xrayutilities.Gridder2D(dim1, dim2)
            if space=="HK":
                self.QGridder(self.H, self.K, self.data)
            elif space=="HL":
                self.QGridder(self.H, self.L, self.data)
            elif space=="KL":
                self.QGridder(self.K, self.L, self.data)
            elif space=="Q":
                self.QGridder(self.QY, self.QZ, self.data)
            self.data = self.QGridder.data.T
            self.MAIN_EXTENT = (self.QGridder.xaxis.min(), self.QGridder.xaxis.max(), self.QGridder.yaxis.min(), self.QGridder.yaxis.max())
            
                
        else:
            self.popup_info('warning','Please correct the geometry of the detector prior to proceed this operation!')
            self.hk_space_btn.setChecked(False)
            self.hl_space_btn.setChecked(False)
            self.kl_space_btn.setChecked(False)
            self.q_space_btn.setChecked(False)
                
            #print( self.MAIN_EXTENT
    
    def plot_data(self):
        """plot the selected edf image"""
        global g_is_image_loaded
        if g_is_image_loaded==False:
            self.popup_info("Info", "No Image has been loaded. Kindly load an image to plot data")
            return
        self.data = self.fabioIMG.data
        self.geometry_corrected = False
        if self.DARK_CORRECTION:
            try:
                self.data = self.data - self.DARK_DATA
                null = (self.data <0)
                self.data[null] = 0
            except:
                pass
        if self.adj_btn.isChecked():
            adjust = True
        else:
            adjust = False
        
        ### Calculate the median and the deviation of this data ###
        self.med, self.nMAD = median_stats(self.data)
        #If the loaded EDF image is already adjusted, the size of the image is not divided by 9600, i.e size%9600 !=0
        # 9600 = 120*80 = the size of one module
        if self.data.size%9600 == 0:
            self.adj_btn.setEnabled(True)
            adjust = self.adj_btn.isChecked()
        else:
            self.adj_btn.setEnabled(False)
            adjust = False
            self.geometry_corrected = True
        
        if adjust:
            self.data = correct_geometry(self.data,detector_type=self.detector_type)
            self.geometry_corrected = True

        # if self.horizontal_detector == True:
            # self.data = N.rot90(self.data) #rotation in the clock-wise direction - right rotation
        self.data = N.rot90(self.data,self.rotate_detector_n)
        self.MAIN_EXTENT = (0, self.data.shape[1], 0, self.data.shape[0])
        
        if self.tth_chi_space_btn.isChecked():
            self.Angular_space_plot()
        elif self.hk_space_btn.isChecked():
            self.Reciprocal_space_plot(space="HK")
        elif self.hl_space_btn.isChecked():
            self.Reciprocal_space_plot(space="HL")
        elif self.kl_space_btn.isChecked():
            self.Reciprocal_space_plot(space="KL")
        elif self.q_space_btn.isChecked():
            self.Reciprocal_space_plot(space="Q")
        else:
            self.show_chi_delta_btn.setEnabled(True)
            self.show_chi_delta_flag = self.show_chi_delta_btn.isChecked()        
        
        self.scale_plot()

        if self.data.shape[0] < self.data.shape[1]:
            self.cb.ax.set_visible(False)
            self.cb2.ax.set_visible(True)
        else:
            self.cb.ax.set_visible(True)
            self.cb2.ax.set_visible(False)
        self.img.set_extent(self.MAIN_EXTENT)
        self.slider_update()

    def plot_profiles(self, x, y, cross_line=True):
        """Line x = 2theta profile, Column y = Chi profile, if not transform in 2theta,chi space, the X,Y profiles are plotted in detector coordinates"""
        integration_width = self.integration_width.text()
        if integration_width=="":
            integration_width = 10
        else:
            integration_width = int(integration_width)
            
        if cross_line:
            x=x[0]
            y=y[0]
            #Draw the cross-lines in the main figure:
            hline = self.ax.axhline(y, color='k', ls='--', lw=1)
            self.lines.append(hline)
            vline = self.ax.axvline(x, color='k', ls='--', lw=1)
            self.lines.append(vline)
            
            if self.tth_chi_space_btn.isChecked():
                x = get_index(self.tth_pyFAI,x)
                y = get_index(self.chi_pyFAI,y)
            elif self.hk_space_btn.isChecked() or self.kl_space_btn.isChecked() or self.hl_space_btn.isChecked() or self.q_space_btn.isChecked():
                x = get_index(self.QGridder.xaxis,x)
                y = get_index(self.QGridder.yaxis,y)
            x= int(x)
            y= int(y)        

            self.profiles_data_X = self.data[y-integration_width/2:y+integration_width/2,:].sum(axis=0)
            self.profiles_data_X = self.profiles_data_X / integration_width

            self.profiles_data_Y = self.data[:, x-integration_width/2:x+integration_width/2].sum(axis=-1)
            self.profiles_data_Y = self.profiles_data_Y / integration_width

            if self.tth_chi_space_btn.isChecked():
                X_label = "2 Theta (deg)"
                Y_label = "Chi (deg)"
                yc = self.chi_pyFAI[y]
                xc = self.tth_pyFAI[x]
                coor_X = self.tth_pyFAI
                coor_Y = self.chi_pyFAI
                self.chi_title.setText("Chi")
                self.tth_title.setText("2 Theta")
            
            elif self.hk_space_btn.isChecked():
                X_label = "H (r.l.u)"
                Y_label = "K (r.l.u)"
                yc = self.QGridder.yaxis[y]
                xc = self.QGridder.xaxis[x]
                coor_X = self.QGridder.xaxis
                coor_Y = self.QGridder.yaxis
                self.chi_title.setText("K")
                self.tth_title.setText("H")
            
            elif self.hl_space_btn.isChecked():
                X_label = "H (r.l.u)"
                Y_label = "L (r.l.u)"
                yc = self.QGridder.yaxis[y]
                xc = self.QGridder.xaxis[x]
                coor_X = self.QGridder.xaxis
                coor_Y = self.QGridder.yaxis
                self.chi_title.setText("L")
                self.tth_title.setText("H")
            elif self.kl_space_btn.isChecked():
                X_label = "K (r.l.u)"
                Y_label = "L (r.l.u)"
                yc = self.QGridder.yaxis[y]
                xc = self.QGridder.xaxis[x]
                coor_X = self.QGridder.xaxis
                coor_Y = self.QGridder.yaxis
                self.chi_title.setText("L")
                self.tth_title.setText("K")
            elif self.q_space_btn.isChecked():
                X_label = r"$Qy\ (nm^{-1})$"
                Y_label = r"$Qz\ (nm^{-1})$"
                yc = self.QGridder.yaxis[y]
                xc = self.QGridder.xaxis[x]
                coor_X = self.QGridder.xaxis
                coor_Y = self.QGridder.yaxis
                self.chi_title.setText("Qz")
                self.tth_title.setText("Qy")
            else:
                X_label = "X (pixels)"
                Y_label = "Y (pixels)"
                xc = x
                yc = y
                coor_X = N.arange(self.profiles_data_X.shape[0])
                coor_Y = N.arange(self.profiles_data_Y.shape[0])
                self.chi_title.setText("Y")
                self.tth_title.setText("X")
        else:
            #This is for arbitrary profile
            if self.tth_chi_space_btn.isChecked():
                x[0] = get_index(self.tth_pyFAI,x[0])
                y[0] = get_index(self.chi_pyFAI,y[0])
                x[1] = get_index(self.tth_pyFAI,x[1])
                y[1] = get_index(self.chi_pyFAI,y[1])
            
            elif self.hk_space_btn.isChecked() or self.hl_space_btn.isChecked() or self.kl_space_btn.isChecked() or self.q_space_btn.isChecked():
                x[0] = get_index(self.QGridder.xaxis,x[0])
                y[0] = get_index(self.QGridder.yaxis,y[0])
                x[1] = get_index(self.QGridder.xaxis,x[1])
                y[1] = get_index(self.QGridder.yaxis,y[1])
                
                
            num = int(N.hypot(x[1]-x[0], y[1]-y[0])) #Number of points to be taken along the line
            #print( "Number of points selected: ",num
            xi, yi = N.linspace(x[0], x[1], num), N.linspace(y[0], y[1], num)
            #self.profiles_data_X = self.profiles_data_Y = self.data[yi.astype(N.int), xi.astype(N.int)]
            self.profiles_data_X = self.profiles_data_Y = ndimage.map_coordinates(self.data, N.vstack((yi,xi)))
            if self.tth_chi_space_btn.isChecked():
                X_label = "2 Theta (deg)"
                Y_label = "Chi (deg)"
                coor_X = N.linspace(self.tth_pyFAI[int(x[0])], self.tth_pyFAI[int(x[1])], num)
                coor_Y = N.linspace(self.chi_pyFAI[int(y[0])], self.chi_pyFAI[int(y[1])], num)

                self.chi_title.setText("Chi")
                self.tth_title.setText("2 Theta")
            
            elif self.hk_space_btn.isChecked() or self.hl_space_btn.isChecked() or self.kl_space_btn.isChecked() or self.q_space_btn.isChecked():
                coor_X = N.linspace(self.QGridder.xaxis[int(x[0])], self.QGridder.xaxis[int(x[1])], num)
                coor_Y = N.linspace(self.QGridder.yaxis[int(y[0])], self.QGridder.yaxis[int(y[1])], num)
                if self.hk_space_btn.isChecked():
                    X_label = "H (r.l.u)"
                    Y_label = "K (r.l.u)"
                    self.chi_title.setText("K")
                    self.tth_title.setText("H")
                elif self.hl_space_btn.isChecked():
                    X_label = "H (r.l.u)"
                    Y_label = "L (r.l.u)"
                    self.chi_title.setText("L")
                    self.tth_title.setText("H")
                elif self.kl_space_btn.isChecked():
                    X_label = "K (r.l.u)"
                    Y_label = "L (r.l.u)"
                    self.chi_title.setText("L")
                    self.tth_title.setText("K")
                elif self.q_space_btn.isChecked():
                    X_label = r"$Qy\ (nm^{-1})$"
                    Y_label = r"$Qz\ (nm^{-1})$"
                    self.chi_title.setText("Qz")
                    self.tth_title.setText("Qy")
                    
            
            
            else:
                X_label = "X (pixels)"
                Y_label = "Y (pixels)"
                coor_X = xi
                coor_Y = yi
                self.chi_title.setText("Y")
                self.tth_title.setText("X")
            tmpX = N.sort(coor_X)
            tmpY = N.sort(coor_Y)
            xc = tmpX[self.profiles_data_X.argmax()]
            yc = tmpY[self.profiles_data_Y.argmax()]
        self.coor_X_export = coor_X
        self.coor_Y_export = coor_Y
        self.profiles_ax1.cla()
        self.profiles_ax2.cla()
        #self.profiles_ax1.format_coord = self.pro_format_coord
        #self.profiles_ax2.format_coord = self.pro_format_coord
        # The CHI or Vertical (Y) profile (ax1):

        self.profiles_ax1.plot(coor_Y, self.profiles_data_Y, color='blue', lw=1.5)
        Y_fitted_params, Y_fitted_data = fit(coor_Y, self.profiles_data_Y, yc, arbitrary= not cross_line)
        self.profiles_ax1.plot(coor_Y, Y_fitted_data, color='red', lw=1, alpha=0.8)
        self.profiles_ax1.set_xlabel(Y_label, size=14)

        # The TTH or Horizontal (X) profile (ax2):
        self.profiles_ax2.plot(coor_X, self.profiles_data_X, color='blue', lw=1.5)
        X_fitted_params, X_fitted_data = fit(coor_X, self.profiles_data_X, xc, arbitrary= not cross_line)
        self.profiles_ax2.plot(coor_X, X_fitted_data, color='red', lw=1, alpha=0.8)
        self.profiles_ax2.set_xlabel(X_label, size=14)
        self.profiles_canvas.draw()
        # Show the fitted results
        self.chi_fitted_y0.setText("%.4f"%Y_fitted_params['y0'].value)
        self.chi_fitted_xc.setText("%.4f"%Y_fitted_params['xc'].value)
        self.chi_fitted_A.setText("%.4f"%Y_fitted_params['A'].value)
        self.chi_fitted_w.setText("%.4f"%Y_fitted_params['w'].value)
        self.chi_fitted_mu.setText("%.4f"%Y_fitted_params['mu'].value)

        self.tth_fitted_y0.setText("%.4f"%X_fitted_params['y0'].value)
        self.tth_fitted_xc.setText("%.4f"%X_fitted_params['xc'].value)
        self.tth_fitted_A.setText("%.4f"%X_fitted_params['A'].value)
        self.tth_fitted_w.setText("%.4f"%X_fitted_params['w'].value)
        self.tth_fitted_mu.setText("%.4f"%X_fitted_params['mu'].value)

        self.profiles_refresh()
        self.canvas.draw()

    def draw_pointed(self, x, y, finished=False):
        p=self.ax.plot(x,y,'ro')
        self.points.append(p[0])
        if finished:
            l=self.ax.plot(self.arb_lines_X, self.arb_lines_Y, '--',linewidth=1.5, color='white')
            self.lines.append(l[0])
        else:
            pass

        #if self.detector_type=="S70":
            #self.ax.axis([self.xDim0,self.xDim1, self.yDim1,self.yDim0])
        #else:
            #self.ax.axis([self.xDim0, self.xDim1, self.yDim0, self.yDim1])
        self.ax.axis([self.MAIN_EXTENT[0],self.MAIN_EXTENT[1], self.MAIN_EXTENT[2], self.MAIN_EXTENT[3]])
        self.canvas.draw()

    def srcFileForPlot(self):
        #options = QFileDialog.Options()
        #options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self, "QFileDialog.getOpenFileName()", "",
                                                  "h5 Files (*.h5)")#, options=options)
        if fileName:
            self.h5SrcTb.setText(fileName)


    def destFldForPlot(self):
        dirName = str(QFileDialog.getExistingDirectory())
        if dirName:
            self.destFldTb.setText(dirName)

    def profiles_refresh(self):
        """ """
        if self.profiles_log_btn.isChecked():
            if len(self.profiles_ax1.get_lines())>0:
                try:
                    self.profiles_ax1.set_yscale('log')
                except ValueError:
                    pass
                
            if len(self.profiles_ax2.get_lines())>0:
                try:
                    self.profiles_ax2.set_yscale('log')
                except ValueError:
                    pass
            else:
                self.profiles_log_btn.setChecked(False)
                self.popup_info("error","The graphs have no data!")

        else:
            self.profiles_ax1.set_yscale('linear')
            self.profiles_ax2.set_yscale('linear')
        self.profiles_ax1.format_coord = self.pro_format_coord
        self.profiles_ax2.format_coord = self.pro_format_coord
        self.profiles_canvas.draw()
        #return

    def profiles_update(self, widget):
        self.profiles_refresh()

    def profiles_export(self,widget):
        """ Export X,Y profiles data in the same folder as the EDF image """
        proX_fname = self.edf.split(".")[0]+"_X_profile.dat"
        proY_fname = self.edf.split(".")[0]+"_Y_profile.dat"
        data_x_export = False
        data_y_export = False
        if len(self.profiles_ax2.get_lines())>0:
            proX_export = self.profiles_ax2.get_lines()[0].get_xydata()
            N.savetxt(proX_fname, proX_export)
            data_x_export = True
        if len(self.profiles_ax1.get_lines())>0:
            proY_export = self.profiles_ax1.get_lines()[0].get_xydata()
            N.savetxt(proY_fname, proY_export)
            data_y_export = True
        MSSG = "Data exported: \n\n"
        if data_x_export:
            MSSG+=proX_fname+"\n\n"
        if data_y_export:
            MSSG+=proY_fname+"\n\n"
        if data_x_export or data_y_export:
            self.popup_info('info',MSSG)
        else:
            self.popup_info('error','ERROR! No data exported!')

    def profile_press(self, event):
        """ Calculate thickness fringes """
        if event.inaxes == self.profiles_ax1:
            draw_fringes = True
            ax = self.profiles_ax1
            ax_data = ax.get_lines()[0].get_xydata()
            X_data = ax_data[:,0]
            Y_data = ax_data[:,1]
            xlabel = 'Y'
            title = "Linear regression of Y fringes"
            title_FFT = "Fast Fourier Transform of Y profile"
            xlabel_FFT= "Frequency"
        elif event.inaxes == self.profiles_ax2:
            draw_fringes = True
            ax = self.profiles_ax2
            ax_data = ax.get_lines()[0].get_xydata()
            X_data = ax_data[:,0]
            Y_data = ax_data[:,1]
            xlabel = 'X'
            title = "Linear regression of X fringes"
            title_FFT = "Fast Fourier Transform of X profile"
            xlabel_FFT= "Frequency"
        else:
            draw_fringes = False
            
        if draw_fringes and (event.button==1):
            if len(self.profiles_fringes)>0:
                self.profiles_fringes = N.asarray(self.profiles_fringes)
                self.profiles_fringes = N.sort(self.profiles_fringes)
                fringes_popup = PopUpFringes(self.profiles_fringes, xlabel, "Fringes order", title)
                self.profiles_fringes=[]
                self.clear_notes()
        elif draw_fringes and (event.button == 3):
            vline=ax.axvline(event.xdata, linewidth=2, color="green")
            self.lines.append(vline)
            self.profiles_fringes.append(event.xdata)
        
        elif draw_fringes and event.button == 2:
            XF,YF = Fourier(X_data, Y_data)
            popup_window=PopUpImage(XF, YF, xlabel_FFT, "Normalized intensity", title_FFT)
            
        self.profiles_canvas.draw()
        
    def draw_rect(self):
        self.rect.set_width(self.x1 - self.x0)
        self.rect.set_height(self.y1 - self.y0)
        self.rect.set_xy((self.x0, self.y0))
        self.rect.set_linestyle('solid')
        self.rect.set_facecolor("white")
        self.rect.set_alpha(0.3)
        self.rect.set_edgecolor("black")
        self.rect.set_visible(self.zoom_press)
        self.canvas.draw()
        
    def Enable_draw_roi(self,w):
        if self.draw_roi_btn.isChecked():
            self.ROI_ON = True
        else:
            self.ROI_ON = False
            self.ROI_DRAWN = False
            
    def draw_roi(self):
        self.roi_rect.set_width(self.ROI_x1 - self.ROI_x0)
        self.roi_rect.set_height(self.ROI_y1 - self.ROI_y0)
        self.roi_rect.set_xy((self.ROI_x0, self.ROI_y0))
        self.roi_rect.set_linestyle('solid')
        self.roi_rect.set_linewidth(2)
        self.roi_rect.set_facecolor("none")
        #self.rect.set_alpha(0.3)
        self.roi_rect.set_edgecolor("red")
        self.roi_rect.set_visible(True)
        
        self.canvas.draw()

    def clear_notes(self):
        self.roi_rect.set_visible(False)
        if len(self.my_notes)>=1:
            for txt in self.my_notes:
                try:
                    txt.remove()
                except ValueError:
                    break
        if len(self.lines)>=1:
            for line in self.lines:
                try:
                    line.remove()
                except ValueError:
                    break
        if len(self.points)>=1:
            for p in self.points:
                try:
                    p.remove()
                except ValueError:
                    break

        self.canvas.draw()
        self.my_notes = []
        self.lines=[]
        self.points=[]
        self.arb_lines_X=[]
        self.arb_lines_Y=[]
        self.arb_line_points = 0

    def peak_max(self, event):
        x = int(event.xdata)
        y = int(event.ydata)
        chi,tth,d=self.calcul_chi_2theta_d(event)
        chi_text = r'$ \chi \ = \ %.2f$'%chi
        tth_text = r'$ 2\theta \ = \ %.2f$'%tth
        d_text = r'$ d \ = \ %.4f \ \AA$'%d
        txt1 = self.ax.text(x,y+50,tth_text, fontsize=14, color="black")
        txt2 = self.ax.text(x,y,chi_text, fontsize=14, color="black")
        txt3 = self.ax.text(x,y-50,d_text, fontsize=14, color="black")
        self.my_notes.append(txt1)
        self.my_notes.append(txt2)
        self.my_notes.append(txt3)
        self.canvas.draw()

    def on_press(self, event):
        #********** Zoom action ***************************************
        if (self.zoom) and (event.inaxes == self.ax) and (event.button==1):
            #print( 'press'
            self.zoom_press = True
            self.x0 = event.xdata
            self.y0 = event.ydata
        #******** ROI plot *********************************
        elif (self.ROI_ON) and (event.inaxes == self.ax) and (event.button==1):
            self.ROI_press = True
            self.roi_rect.set_visible(False)
            self.ROI_x0 = event.xdata
            self.ROI_y0 = event.ydata
            self.SCAN_ROI_INTEGRATION_X = []
            self.SCAN_ROI_INTEGRATION_Y = []
            
        #******** Plot cross profiles *********************************
        elif (event.inaxes == self.ax) and (event.button==3) and self.plotXYprofiles_btn.isChecked():
            x = event.xdata
            y = event.ydata
            xx=[]
            yy=[]
            xx.append(x)
            yy.append(y)
            try:
                self.clear_notes()
            except:
                pass
            self.plot_profiles(xx,yy,cross_line=True)
        #******** Plot arbitrary profiles  ****************************
        elif (event.inaxes == self.ax) and (event.button==1) and self.arbitrary_profiles_btn.isChecked():
            self.arb_line_points +=1
            #print( "Number of points clicked: ",self.arb_line_points
            if self.arb_line_points>2:
                self.clear_notes()
                self.arb_line_points=1

            x = event.xdata
            y = event.ydata
            self.arb_lines_X.append(x)
            self.arb_lines_Y.append(y)
            if len(self.arb_lines_X)<2:
                finished=False
            elif len(self.arb_lines_X)==2:
                finished = True

            self.draw_pointed(x,y,finished)#If finished clicking, connect the two points by a line
            if finished:
                self.plot_profiles(self.arb_lines_X, self.arb_lines_Y, cross_line=False)
                self.arb_lines_X=[]
                self.arb_lines_Y=[]
        #******** Mark the peak information ****************************
        elif (event.inaxes == self.ax) and (event.button==3) and (self.show_chi_delta_flag == True):
            self.peak_max(event)

        elif event.button==2:
            self.clear_notes()
        else:
            return

    def on_motion(self,event):
        if event.inaxes == self.ax:
            self.status_update(event)
            if self.zoom_press:
                self.mouse_moved = True
                self.x1 = event.xdata
                self.y1 = event.ydata
                self.draw_rect()
            elif self.ROI_press:
                self.mouse_moved = True
                self.ROI_x1 = event.xdata
                self.ROI_y1 = event.ydata
                self.draw_roi()
        else:
            return

    def on_release(self, event):
        if (self.zoom) and (event.inaxes == self.ax):
            #print( 'release'

            self.zoom_press = False
            self.x1 = event.xdata
            self.y1 = event.ydata
            self.draw_rect()
            if self.mouse_moved==True:
                self.zoom_image()
                self.mouse_moved = False
                
        if (self.ROI_ON) and (event.inaxes == self.ax):
            #print( 'release'

            self.ROI_press = False
            self.ROI_x1 = event.xdata
            self.ROI_y1 = event.ydata
            self.draw_roi()
            self.ROI_DRAWN = True

    def check_input_data(self):
        """ Check if the user has entered correctly the data needed to construct the pole figure
        We need the start and end number of images, and the 2 theta for which we want to plot """
        try:
            self.img_deb = int(self.images_from_nb.text())
            self.img_fin = int(self.images_to_nb.text())
            self.pole_2theta = float(self.pole_2theta_field.text())
            #self.plot_chi_min = float(self.PF_chi_min.text())
            return True
        except:
            self.popup_info("error","Values are not entered correctly. Please check again!")
            return False

    
           
    def load_data(self, img_deb, img_fin, pole_2theta):
        """
        folder: directory where edf images are found
        img_deb, img_fin: number of beginning and ending images for the maps
        pole_2theta: 2theta (degrees) to construct the pole figure
        logscale: (1 or 0) to present the pole figure in log scale or linear scale
        """
        gc.collect()
        self.data_loading.show()
        #Check where the images are
        for k in self.store_img.keys():
            if (img_deb in self.store_img[k]) and (img_fin in self.store_img[k]):
                self.edf_folder = k
                break
        
        img_list  = select_files_from_list(self.store[self.edf_folder], img_deb, img_fin)
        
        img = fabio.open(join(self.edf_folder, img_list[0]))
        this_motor = get_motors(img.header)
        if not self.calibrated_quantitative:
            rot1 = N.radians(this_motor['nu'])*(-1.)
            rot2 = N.radians(this_motor['del'])*(-1.)
            self.azimuthalIntegration.rot1 = rot1
            self.azimuthalIntegration.rot2 = rot2
            self.azimuthalIntegration.rot3 = 0
        total = len(img_list)
        processed = 0.
        
        threads = []
        if self.select_phi.isChecked():            
            plot_kphi = False
        elif self.select_kphi.isChecked():
            plot_kphi = True
        for img in range(total):
            edf_basename = img_list[img]
            t = Pole_Figure_load_data(img, self.edf_folder, edf_basename, pole_2theta, self.azimuthalIntegration, plot_kphi, detector_type=self.detector_type)
            threads.append(t)
            t.start()
            processed +=1.
            fr = processed/total
            self.data_loading.set_fraction(fr)
            self.data_loading.setText("Data loading: "+str(int(fr*100))+"%")
            self.data_loading.set_show_text(True)
            while gtk.events_pending():
                gtk.main_iteration()
        output = []
        for t in threads:
            t.join()
            output.append(t.Data)
            chi = t.chi
            self.sample_tilt = t.chi_gonio
            self.omega = t.omega
            
        while len(threads)>0:
            threads.pop()
        #gobject.timeout_add(200, self._callback)
        output.sort()
        intensity = [o[1] for o in output]
        phi_table = [o[2] for o in output]
        intensity = N.asarray(intensity)
        phi_table = N.asarray(phi_table)
        
        if self.detector_type=="S70":
            chi = chi + 90
        else:
            chi = chi - 90
        del threads
        return chi, phi_table, intensity

    def coordinates_transform(self, psi, phi, tth):
        #Transformer les coordonnees du detecteur en coordonnees polaires, avec phi et psi sont des matrices 2D, omega=angle incidente, tth=2theta bragg
        #Ref: Marie-Ingrid Richard, J. Appl. Cryst. 46,(2013),1842
        #tilt = self.sample_tilt
        tilt = 90 - self.sample_tilt
        # tilt = 0
        omega = self.omega
        phi   = N.radians(phi)
        psi   = N.radians(psi)
        
        cosPhi = N.cos(phi)
        sinPhi = N.sin(phi)
        cosPsi = N.cos(psi)
        sinPsi = N.sin(psi)
        cosTil = N.cos(N.radians(tilt))
        sinTil = N.sin(N.radians(tilt))
        cosTta = N.cos(N.radians(tth/2.))
        sinTta = N.sin(N.radians(tth/2.))
        cosOme = N.cos(N.radians(omega))
        sinOme = N.sin(N.radians(omega))
        
        #Phi: rotate around z
        #Omega: Rotate around y
        #Chi: rotate around -x
        Rx = [[1,0,0], [0,cosTil,sinTil], [0,-sinTil,cosTil]]
        Ry = [[cosOme,0,sinOme], [0,1,0], [-sinOme,0,cosOme]]
        Rz = [[cosPhi,-sinPhi,0], [sinPhi, cosPhi,0], [0,0,1]]
        qLab = [-sinTta, -cosTta*sinPsi, cosTta*cosPsi]
        
        Rx = N.asarray(Rx)
        Ry = N.asarray(Ry)
        Rz = N.asarray(Rz)
        #qLab=N.asarray(qLab)
        
        qSample = Rz.dot(Rx).dot(Ry).dot(qLab)
        psi_pole = N.arctan2(N.sqrt(qSample[0]**2 + qSample[1]**2), qSample[2])
        phi_pole = N.arctan2(qSample[1], qSample[0])
        return psi_pole, phi_pole
    
    def plot_pole_figure(self,widget):
        """ Click the PLOT button to plot the pole figure """
        check = self.check_input_data()
        if check:
            # if self.detector_type is not "D1":
            self.PF_psi, self.PF_phi, self.PF_intensity = self.load_data(self.img_deb, self.img_fin, self.pole_2theta)
            self.plot_PF()
        else:
            self.popup_info("error","Please make sure that you have correctly entered all of the fields!")
        self.data_loading.hide()

    def validatePoleWidgets(self):
        if self.h5SrcTb.text()=='':
            self.popup_info("error","Kindly provide h5 sample value!")
            self.h5SrcTb.setFocus()
            return True
        if self.wavelengthTb .text()=='':
            self.popup_info("error","Kindly provide wavelength values!")
            self.wavelengthTb.setFocus()
            return True
        if self.startTb.text()=='':
            self.popup_info("error","Kindly provide start values!")
            self.startTb.setFocus()
            return True
        if self.endTb.text()=='':
            self.popup_info("error","Kindly provide end values!")
            self.endTb.setFocus()
            return True
        if self.omegaTb.text()=='':
            self.popup_info("error","Kindly provide omega values!")
            self.omegaTb.setFocus()
            return True
        #h5SrcTb
        if self.h5SrcTb.text()=='':
            self.popup_info("error","Kindly provide source h5 files!")
            self.h5SrcTb.setFocus()
            return True
        if self.destFldTb.text()=='':
            self.popup_info("error","Kindly provide destination folder!")
            self.destFldBtn.setFocus()
            return True
        if self.setpTb.text()=='':
            self.popup_info("error","Kindly provide step values!")
            self.setpTb.setFocus()
            return True
        return False

    def validate1PoleWidgets(self):
        if self.h5SrcTb.text()=='':
            self.popup_info("error","Kindly provide h5 sample value!")
            self.h5SrcTb.setFocus()
            return True
        if self.wavelengthTb .text()=='':
            self.popup_info("error","Kindly provide wavelength values!")
            self.wavelengthTb.setFocus()
            return True

        if self.omegaTb.text()=='':
            self.popup_info("error","Kindly provide omega values!")
            self.omegaTb.setFocus()
            return True
        #h5SrcTb
        if self.h5SrcTb.text()=='':
            self.popup_info("error","Kindly provide source h5 files!")
            self.h5SrcTb.setFocus()
            return True
        if self.destFldTb.text()=='':
            self.popup_info("error","Kindly provide destination folder!")
            self.destFldBtn.setFocus()
            return True
        return False

    def generatePlotAlllFigure(self):
        print("Came in generating plot figures")
        if not self.validatePoleWidgets():
            #print("Came in generating plot figures")
            edf_folders = [""]
            samples = [self.h5SrcTb.text()]
            _vmin=self.sld_vmin.value()
            _vmax=self.sld_vmin.value()
            log=self.logTb.isChecked() #True
            wavelength =self.wavelengthTb.text() #0.7044e-10
            omega_25 = float(self.omegaTb.text()) #3.5#12.5
            start = int(self.startTb.text())
            end = int(self.endTb.text() )
            step = float(self.setpTb.text() ) 
            tth_list_del25 = N.arange(start, end, step)
            source = self.h5SrcTb.text() #"S159-160.h5"
            destination = self.destFldTb.text() #"polefigure/S159-160/"
            theta = int(self.pole_2theta_field.text())
            params = {'sample': samples, 'log': log, 'wl': wavelength, 'start': start, 'end': end, 'step': step, "tth": tth_list_del25, 'omega': omega_25, 'src': source, 'dest': destination,
                    'vmax': _vmax, 'vmin': _vmin, 'theta': theta}
            pole_2theta = theta
            delta = 3.5
            #print( params )
            msg = draw_all(source, destination, makemovie=True, params=params)
            self.popup_info("Info", msg)
            #draw_one(source, delta, pole_2theta, params=params)


    def generatePlotFigure(self):
        print("Came in generating plot figures")
        if not self.validate1PoleWidgets():
            #print("Came in generateing plot figures")
            edf_folders = [""]
            samples = [self.h5SrcTb.text()]
            _vmin=self.sld_vmin.value()
            _vmax=self.sld_vmin.value()
            log=self.logTb.isChecked() #True
            wavelength =self.wavelengthTb.text() #0.7044e-10
            omega_25 = float(self.omegaTb.text()) #3.5#12.5
  #          start = int(self.startTb.text())
   #         end = int(self.endTb.text() )
    #        step = float(self.setpTb.text() ) 
     #       tth_list_del25 = N.arange(start, end, step)
            source = self.h5SrcTb.text() #"S159-160.h5"
            destination = self.destFldTb.text() #"polefigure/S159-160/"
            theta = int(self.pole_2theta_field.text())
            params = {'sample': samples, 'log': log, 'wl': wavelength, 'omega': omega_25, 'src': source, 'dest': destination,
                    'vmax': _vmax, 'vmin': _vmin, 'theta': theta}
            pole_2theta = theta
            delta = 3.5
            res = draw_one(source, delta, pole_2theta, params=params)
            self.popup_info("Info", res)

    def replot_PF(self,widget):
        """ Replot the pole figure with logscale option"""
        self.PF_intensity = self.change_scale(widget,self.PF_intensity)
        
        
    def plot_PF(self):
        #********** Plot the pole figure ***************
        gc.collect()
        self.polar_ax.clear()
        self.polar_cax.clear()
        PF_data = self.PF_intensity.copy()
        if self.log_scale == 1:
            PF_data = flat_data(PF_data, self.vmin, self.vmax, True)
            PF_data = N.log10(PF_data+1e-6)
            clabel = r'$Log_{10}\ (Counts\ per\ second)\ [arb.\ units]$'
            fmt = "%.1f"
        else:
            clabel = r'$Intensity\ (Counts\ per\ second)$'
            fmt = "%d"
            dynhigh = N.where(PF_data > self.vmax)
            dynlow  = N.where(PF_data < self.vmin)
            PF_data[dynhigh] = self.vmax
            PF_data[dynlow]  = self.vmin
        
        psi,phi   = N.meshgrid(self.PF_psi,self.PF_phi)
        
        psi,phi = self.coordinates_transform(psi, phi, self.pole_2theta)
        X = psi*N.cos(phi)
        Y = psi*N.sin(phi)
        clevels = N.linspace(self.vmin, self.vmax,50)
        
        self.polar_img = self.polar_ax.contourf(X,Y, PF_data, clevels, vmin=self.vmin, vmax=self.vmax)
        self.polar_ax.set_aspect('equal')
        self.polar_ax.set_frame_on(False)
        self.polar_ax.get_xaxis().set_ticks([])
        self.polar_ax.get_yaxis().set_ticks([])
        self.polar_ax.format_coord = self.pole_format_coord
        self.polar_cb  = self.polefig.colorbar(self.polar_img,cax=self.polar_cax, format=fmt)
        self.polar_cb.set_label(clabel, fontsize=12)
        self.polar_cb.locator = MaxNLocator(nbins=6)
        title   = "2 Theta = %.2f Deg."%self.pole_2theta
        self.polar_ax.text(0.5, 1.08, title, horizontalalignment='center', transform = self.polar_ax.transAxes, verticalalignment="center", fontsize=12)
        plotPolarGrid(self.polar_ax, abs(X).min(), X.max(), 6, 8)
        self.polar_ax.relim()
        self.pole_canvas.draw()
    #***************************************************************************
    #     For batch images correction 
    #***************************************************************************
    def select_source_folder(self, widget):
        #dialog = gtk.FileChooserDialog(title="Select a sources folder",action=gtk.FILE_CHOOSER_ACTION_SELECT_FOLDER, buttons = (gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL, gtk.STOCK_OPEN, gtk.RESPONSE_OK))
        #dialog.set_current_folder(self.current_folder)
        dialog = QFileDialog.getOpenFileNames(self, "Select a sources folder", self.current_folder, "")#Poni files (*.poni)")
        #response = dialog.run()
        if dialog!='' and type(dialog)==type((1,2)): #response==gtk.RESPONSE_OK:
            dialog = dialog[0]
            if dialog=='':
                return 
            folder=dialog#.get_filename()
            self.t1_src_path.setText(folder)
            self.current_folder = folder
            self.src_folder = folder
        else:
            pass
        #dialog.destroy()
    
    def select_destination_folder(self, widget):
        #dialog = gtk.FileChooserDialog(title="Select the destination folder",action=gtk.FILE_CHOOSER_ACTION_SELECT_FOLDER, buttons = (gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL, gtk.STOCK_OPEN, gtk.RESPONSE_OK))
        #dialog.set_current_folder(self.current_folder)
        #response = dialog.run()
        dialog = QFileDialog.getOpenFileNames(self, "Select a destination folder", self.current_folder, "")#Poni files (*.poni)")
        if dialog!='' and type(dialog)==type((1,2)): #response==gtk.RESPONSE_OK:
            dialog = dialog[0]
            if dialog=='':
                return 
            folder=dialog#.get_filename()
            self.t1_des_path.setText(folder)
            self.current_folder = folder
            self.des_folder = folder
        else:
            pass
        #dialog.destroy()
        
    def select_dark_image(self,w):
        #dialog = gtk.FileChooserDialog("Select a dark image",None,gtk.FILE_CHOOSER_ACTION_OPEN,(gtk.STOCK_CANCEL,gtk.RESPONSE_CANCEL, gtk.STOCK_OPEN, gtk.RESPONSE_OK))
        #dialog.set_current_folder(self.current_folder)
        #response = dialog.run()
        dialog = QFileDialog.getOpenFileNames(self, "Select a destination folder", self.current_folder, "")#Poni files (*.poni)")
        if dialog!='' and type(dialog)==type((1,2)): # response == gtk.RESPONSE_OK:
            dialog = dialog[0]
            if dialog=='':
                return 
            file_choosen = dialog#.get_filename()
            self.t1_dark_img_path.setText(file_choosen)
            self.current_folder = os.path.dirname(file_choosen)
            self.dark_img_file = file_choosen#.decode('utf8')
            self.dark_img = fabio.open(self.dark_img_file)
            self.batch_DARK_CORRECTION = True
            
        else:
            self.batch_DARK_CORRECTION = False
            pass
        #dialog.destroy()
    
    def show_proc(self,out_info):
        self.show_process_info.set_label(out_info)
    
    def batch_change_detector(self, widget):
        detector = self.t2_det_combobox.currentText()
        if detector == "D1":
            self.d1_center.show()
            self.d1_center_txt.show()
            self.d1_distance.show()
            self.d1_distance_txt.show()
            self.d1_specfile_browse.show()
            self.d1_specfile_path.show()
            self.d1_specfile_txt.show()
        else:
            self.d1_center.hide()
            self.d1_center_txt.hide()
            self.d1_distance.hide()
            self.d1_distance_txt.hide()
            self.d1_specfile_browse.hide()
            self.d1_specfile_path.hide()
            self.d1_specfile_txt.hide()
        
    def select_normalisation_file(self,widget,path):
        #dialog = gtk.FileChooserDialog("Select spec file",None,gtk.FILE_CHOOSER_ACTION_OPEN,(gtk.STOCK_CANCEL,gtk.RESPONSE_CANCEL, gtk.STOCK_OPEN, gtk.RESPONSE_OK))
        #dialog.set_current_folder(self.current_folder)
        #response = dialog.run()
        dialog = QFileDialog.getOpenFileNames(self, "Select a spec folder", self.current_folder, "")#Poni files (*.poni)")
        if dialog!='' and type(dialog)==type((1,2)): #response == gtk.RESPONSE_OK:
            dialog = dialog[0]
            if dialog=='':
                return 
            file_choosen = dialog#.get_filename()
            path.set_text(file_choosen)
            self.current_folder = os.path.dirname(file_choosen)
            #if label == "A":
                #self.attenuation_file = file_choosen.decode('utf8')
            #elif label == "S":
            self.spec_file = file_choosen.decode('utf8')
            #elif label == "M":
                #self.mca_file = file_choosen.decode('utf8')
        else:
            pass
        #dialog.destroy()
        
    def geometric_correction_D1(self,edf,adjusted_folder):
        distance = self.d1_distance.text()
        center   = self.d1_center.text()
        normalize= self.spec_file
        program  = 'xpad3_geometry.py'
        cmd = program+" --geometry=D1 --mask=search --distance="+distance+" --center="+center+" --normalize="+normalize+" --monitor=d0_cps --multiply="+str(self.monitor_ref)+" "+edf
        print( cmd)
        os.system(cmd)
        menage = "mv *g.edf %s"%self.des_folder
        os.system(menage)
        
    def combine_edf(self, edf_1, edf_2, normalisation=False, dark_data=None, detector_type="D5"):
        #normalisation: True, False
        #dark_img: fabio object for dark image
        if detector_type == "D5":
            DET_SIZE_X = 1153
            DET_SIZE_Y = 578
        elif detector_type == "S70":
            DET_SIZE_X = 120
            DET_SIZE_Y = 578
            
        img_1     = fabio.open(edf_1)
        img_2     = fabio.open(edf_2)
        data_1    = img_1.data
        data_2    = img_2.data
        if dark_data != None:
            data_1 = data_1 - dark_data
            data_2 = data_2 - dark_data
            data_1[data_1<0] = 0
            data_2[data_2<0] = 0
        motor_1   = get_motors(img_1.header)
        motor_2   = get_motors(img_2.header)
        counter_1 = get_counters(img_1.header)
        counter_2 = get_counters(img_2.header)
        
        if data_1.size%9600==0:
            data_1 = EDF_XZ_combination.correct_geometry(data_1, detector_type = detector_type)
            #img_1.data = d
        if data_2.size%9600==0:
            data_2 = EDF_XZ_combination.correct_geometry(data_2, detector_type = detector_type)
            #img_2.data = d
        
        X1 = motor_1['Xdet']
        X2 = motor_2['Xdet']
        Z1 = motor_1['Zdet']
        Z2 = motor_2['Zdet']
        mon1= counter_1[self.t2_mon_combobox.isChecked_text()]
        mon2= counter_2[self.t2_mon_combobox.isChecked_text()]
        
        deltaX = X2-X1
        deltaZ = Z2-Z1
        pixX   = int(abs(deltaX)*1e-3/_PIXEL_SIZE)
        pixZ   = int(abs(deltaZ)*1e-3/_PIXEL_SIZE)
        mat_decal_Z = N.zeros(shape=(pixZ, DET_SIZE_Y))
        mat_decal_X = N.zeros(shape=(DET_SIZE_X+pixZ,pixX))
        
        if detector_type == "D5":
            data_1 = EDF_XZ_combination.image_correction(data_1)
        data_1 = N.rot90(data_1)
        if deltaZ < 0:
            data_1 = N.vstack((mat_decal_Z,data_1))
        else:
            data_1 = N.vstack((data_1, mat_decal_Z))
        if deltaX >0:
            data_1 = N.hstack((mat_decal_X,data_1))
        else:
            data_1 = N.hstack((data_1, mat_decal_X))
        norm_1 = data_1 == 0.

        if detector_type == "D5":
            data_2 = EDF_XZ_combination.image_correction(data_2)
        data_2 = N.rot90(data_2)
        if deltaZ < 0:
            data_2 = N.vstack((data_2,mat_decal_Z))
        else:
            data_2 = N.vstack((mat_decal_Z, data_2))
        if deltaX > 0:
            data_2 = N.hstack((data_2,mat_decal_X))
        else:
            data_2 = N.hstack((mat_decal_X, data_2))
        norm_2 = data_2 == 0.
        
        if normalisation:
            data_1 = data_1 * self.monitor_ref/mon1
            data_2 = data_2 * self.monitor_ref/mon2
        
        data = (data_1+data_2)/2.0
        norm = norm_1 + norm_2
        data[norm] = data[norm] * 2.0
        #Saving the final image
        img_1.data = data
        return img_1
        
    def transform(self,edf,adjusted_folder,normalisation):
        """ transformer edf en edf adjusted (en rajoutant les gaps)
        Corriger les lignes mortes, nettoyer les pixels chauds
        Normaliser avec le Imachine = 200 mA: Icorr = Imesure*monitor_ref(@200 mA)/monitor_mesure, cela donne les intensites en CPS
        """
        detector_type = self.t2_det_combobox.currentText()
        
        edf_img = fabio.open(edf)
        data = edf_img.data#edfFile.data
        header = edf_img.header#edfFile.header
        counter_name = header['counter_mne'].split()
        counter_value= header['counter_pos'].split()
        counters= {}
        for i in range(len(counter_name)):
            counters[counter_name[i]] = float(counter_value[i])
        monitor_col = self.t2_mon_combobox.isChecked_text()
        monitor = counters[monitor_col]
        #Ajouter les gaps #
        data = correct_geometry(data, detector_type=detector_type)
        #Image cleaning for XPAD D5
        if self.t2_det_combobox.currentText() == "D5":
            data = EDF_XZ_combination.image_correction(data)
        #Normalisation
        if normalisation:
            monitor_ref = self.monitor_ref
            norm_factor = monitor_ref/monitor
            data = data * norm_factor
        #sauver EDF apres correction
        name = os.path.basename(edf)
        name_adjusted = name.split(".")[0]+"_corrected"
        if self.ascii_out.isChecked():
            ext = "dat"
        else:
            ext  = name.split(".")[1]
        #filename = adjusted_folder+"/"+name_adjusted+"."+ext
        filename = join(adjusted_folder, name_adjusted+"."+ext)
        if self.ascii_out.isChecked():
            N.savetxt(filename, data, header=str(header))
        else:
            edf_img.data = data
            edf_img.write(filename)
        
    def processing(self):
        self.progressbar.show()
        #try:
        src_folder = self.src_folder
        des_folder = self.des_folder
        try:
            img_beg = int(self.t2_img_start_entry.text())
            img_end = int(self.t2_img_end_entry.text())
            beg_end = True
        except:
            beg_end = False
        try:
            self.monitor_ref = float(self.t2_ref_mon_entry.text())
            normalisation = True
        except:
            self.monitor_ref = 1
            normalisation = False
        
        filelist=[i for i in listdir(src_folder) if isfile(join(src_folder,i)) and i.endswith(".edf") or i.endswith(".edf.gz") or i.endswith(".dat") or i.endswith(".dat.gz")]
        if beg_end:
            edf_list = select_files_from_list(filelist, img_beg, img_end)
        else:
            edf_list = filelist
        total = len(edf_list)
        processed = 0
        img_number_list = N.arange(img_beg, img_end+1)
        if self.batch_DARK_CORRECTION:
            dark_data=self.dark_img.data
            if normalisation:
                dark_data = dark_data * self.monitor_ref/get_counters(self.dark_img.header)[self.t2_mon_combobox.isChecked_text()]
        else:
            dark_data = None
            
        for e in range(len(edf_list)):
            #try:
            if not self.combine_XY.isChecked():
                edf_base = edf_list[e]
                edf = join(src_folder, edf_base)
                if self.t2_det_combobox.currentText() == "D1":
                    self.geometric_correction_D1(edf, des_folder)
                else:
                    self.transform(edf, des_folder, normalisation)
                out_info = "Image %s saved successfully!"%edf_base.split(".")[0]
            else:
                i = 2*e+1
                if i<len(edf_list):
                    edf_1_base = edf_list[i-1]
                    edf_2_base = edf_list[i]
                    edf_1      = join(src_folder, edf_1_base)
                    edf_2      = join(src_folder, edf_2_base)
                    if self.t2_det_combobox.currentText() != "D1":
                        combined_img = self.combine_edf(edf_1, edf_2, normalisation=normalisation, dark_data=dark_data, detector_type = self.t2_det_combobox.currentText())
                        #sauver EDF apres correction
                        name = "Combined_%04d_%04d"%(img_number_list[i-1],img_number_list[i])
                        
                        if self.ascii_out.isChecked():
                            ext = "dat"
                        else:
                            ext  = "edf"
                        fname = name+"."+ext
                        filename = join(des_folder,fname)
                        if self.ascii_out.isChecked():
                            N.savetxt(filename, combined_img.data, header=str(combined_img.header))
                        else:
                            combined_img.write(filename)
                        out_info = "Image %s saved successfully!"%name
                    else:
                        out_info = "This is not applied to D1 detector"
                    
            #except:
                #out_info = "Image %s does not existe or cannot be corrected!"%edf
            processed +=1.
            fraction = processed/total *100.
            fraction = int(fraction)
            self.progressbar.setValue(fraction/100.)
            self.progressbar.setText(str(fraction)+"%")
            self.progressbar.setTextVisible(True)
            #while gtk.events_pending():
                #gtk.main_iteration()
            self.show_process_info.setText(out_info)
            yield True    
        #except:
            #self.popup_info("warning","Please check that you have carefully entered all of the fields.")
        #self.progressbar.hide()
        yield False
        
    def batch_transform(self,widget):
        task = self.processing()
       # gobject.idle_add(task.next)
    def _callback(self):
        if threading.active_count() == 1:  # If only one left, scanning is done
            return False  # False make callback stop
        #print( threading.active_count()
        return True
#***********************************************************************************************
#                                       FINISHED
#***********************************************************************************************
if __name__ == "__main__":
    #print("Came in main application")
    app = QApplication(sys.argv)
    #gtk.threads_init()
    m=MyMainWindow()
    #m.showFullScreen()
    m.setWindowState(Qt.WindowMaximized)
    #m.popup_info("Info","Kindly select a detector first")
    m.showDetectorDialog()
    #self.set_size_request(1000, 850)
    #m.setGeometry(50, 50, 500, 450)
    m.show()
    sys.exit(app.exec_())
    #gtk.main()
