import sys,os
from PyQt5.QtWidgets import *
#QApplication, QWidget, QPushButton,QFileDialog,QAction,QMainWindow, QLineEdit,QLabel,QGridLayout,QVBoxLayout,QGroupBox,QMessageBox,QInputDialog
from PyQt5.QtGui import QIcon,QIntValidator
from PyQt5.QtCore import pyqtSlot

from conv_h5_2 import H5Img
from PyQt5.QtCore import *
from PyQt5.QtGui import *
import PyQt5
from PyQt5.uic import loadUi
from PyQt5 import uic

from save_h5_mg import *

progress_status = 0

class App(QWidget):
    def __init__(self):
        super().__init__()
        #loadUi('layout4_css.ui', self)
        self.files = []
        self.title = 'Layout 4'
        self.left = 200
        self.top = 200
        self.width = 600
        self.height = 250
        self.initUI()



    def initUI(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)
        label_css = "QLabel { \
                        font-weight: bold; \
                        font-size: 12px; \
                    }"
        text_css = "QLineEdit { \
                        padding: 1px; \
                        border-style: solid; \
                        border: 2px solid gray; \
                        border-radius: 8px; \
                    }"
        btn_css = "QPushButton { \
                    color: black \
                    border-width: 1px; \
                    border-color: #339; \
                    border-style: solid; \
                    border-radius: 7; \
                    padding: 3px; \
                    font-size: 10px; \
                    padding-left: 5px; \
                    padding-right: 5px; \
                    min-width: 50px; \
                    max-width: 50px; \
                    min-height: 18px; \
                    max-height: 18px; \
                   }"
        # Create textbox and lable
        self.lable1 = QLabel(self)
        self.lable1.setText('First Image No: ')
        self.textbox1 = QLineEdit(self)
        self.textbox1.setValidator(QIntValidator())
        self.textbox1.setStyleSheet( text_css )
        self.lable1.setStyleSheet(label_css)
        #self.textbox1.setText('709')

        self.lable2 = QLabel(self)
        self.lable2.setText('Source Foler: ')
        self.lable2.setStyleSheet(label_css)
        self.textbox2 = QLineEdit(self)
        self.textbox2.setStyleSheet( text_css )
        #self.textbox2.setText("")

        self.lable3 = QLabel(self)
        self.lable3.setText('Number of Images: ')
        self.lable3.setStyleSheet(label_css)
        self.textbox3 = QLineEdit(self)
        self.textbox3.setValidator(QIntValidator())
        self.textbox3.setStyleSheet( text_css )
        #self.textbox3.setText('3')

        self.lable4 = QLabel(self)
        self.lable4.setText('Select Poni: ')
        self.lable4.setStyleSheet(label_css)
        self.textbox4 = QLineEdit(self)
        self.textbox4.setStyleSheet( text_css )
        #self.textbox4.setText("17Mar10D5_0075-rsz.poni")

        self.flatLbl = QLabel("Flat File")
        self.darkLbl = QLabel("Dark File")
        self.maskLbl = QLabel("Mask File")
        self.flatCB = QLineEdit()
        self.darkCB = QLineEdit()
        self.maskCB = QLineEdit()
        self.flatBtn = QPushButton('Select EDF for Flat', self)
        self.flatBtn.clicked.connect(self.flatFuncSel)
        self.darkBtn = QPushButton('Select EDF for Dark', self)
        self.darkBtn.clicked.connect(self.darkFuncSel)
        self.maskBtn = QPushButton('Select EDF for Mask', self)
        self.maskBtn.clicked.connect(self.maskFuncSel)
        self.flatLbl.setStyleSheet(label_css)
        self.darkLbl.setStyleSheet(label_css)
        self.maskLbl.setStyleSheet(label_css)
        self.flatCB.setStyleSheet(text_css)
        self.darkCB.setStyleSheet(text_css)
        self.maskCB.setStyleSheet(text_css)
        self.flatBtn.setStyleSheet(btn_css)
        self.darkBtn.setStyleSheet(btn_css)
        self.maskBtn.setStyleSheet(btn_css)


        #set buttons properties.
        self.button1 = QPushButton('Integrate and View', self)
        self.button1.clicked.connect(self.on_click1)

        self.button2 = QPushButton('Integrate and Save', self)
        self.button2.clicked.connect(self.on_click2)


        self.int2dViewBtn = QPushButton('Integrate2d and View', self)
        self.int2dViewBtn.clicked.connect(self.int2dViewFunc)

        self.int2dSaveBtn = QPushButton('Integrate2d and Save', self)
        self.int2dSaveBtn.clicked.connect(self.int2dSaveFun)

        self.int2dLogBtn = QPushButton('Integrate2d and ViewLog', self)
        self.int2dLogBtn.clicked.connect(self.int2dLogFun)

        self.button3 = QPushButton('Select Poni Name', self)
        self.button3.clicked.connect(self.on_click3)

        self.button4 = QPushButton('Select Source Path', self)
        self.button4.clicked.connect(self.on_click4)

        self.imposeBtn = QPushButton('Super Impost and View', self)
        self.imposeBtn.clicked.connect(self.superImposeViewFun)
        self.imposeSaveBtn = QPushButton('Super Impost and View', self)
        self.imposeSaveBtn.clicked.connect(self.superImposeSaveFun)

        self.button1.setStyleSheet( btn_css )
        self.button2.setStyleSheet( btn_css )
        self.button3.setStyleSheet( btn_css )
        self.button4.setStyleSheet( btn_css )
        self.imposeSaveBtn.setStyleSheet( btn_css )
        self.imposeBtn.setStyleSheet( btn_css )
        self.int2dLogBtn.setStyleSheet( btn_css )
        self.int2dSaveBtn.setStyleSheet( btn_css )
        self.int2dViewBtn.setStyleSheet( btn_css )

        self.prefixLbl = QLabel(self)
        self.prefixLbl.setText('Prefix : ')
        self.prefixLbl.setStyleSheet(label_css)
        self.prefixTB = QLineEdit(self)
        #self.prefixTB.setValidator(QIntValidator())
        self.prefixTB.setStyleSheet( text_css )


        #setup the layout: We will use Grid layout to display it. 4x4
        self.grid = QGridLayout()
        self.grid.setColumnStretch(1, 4)
        self.grid.setColumnStretch(2, 4)
        self.grid.addWidget(self.lable1,0,0)
        self.grid.addWidget(self.textbox1, 0, 1)
        self.grid.addWidget(self.lable2, 0, 2)
        self.grid.addWidget(self.textbox2, 0, 3)
        self.grid.addWidget(self.lable3, 1, 0)
        self.grid.addWidget(self.textbox3, 1, 1)
        self.grid.addWidget(self.button4, 1, 3)
        self.grid.addWidget(self.lable4, 2, 0)
        self.grid.addWidget(self.textbox4, 2, 1)
        self.grid.addWidget(self.button3, 2, 2)
        self.grid.addWidget(self.flatLbl, 3, 0)
        self.grid.addWidget(self.darkLbl, 4, 0)
        self.grid.addWidget(self.maskLbl, 5, 0)
        self.grid.addWidget(self.flatCB, 3, 1)
        self.grid.addWidget(self.darkCB, 4, 1)
        self.grid.addWidget(self.maskCB, 5, 1)
        self.grid.addWidget(self.flatBtn, 3, 2)
        self.grid.addWidget(self.darkBtn, 4, 2)
        self.grid.addWidget(self.maskBtn, 5, 2)
        self.grid.addWidget(self.prefixTB, 6, 2)
        self.grid.addWidget(self.prefixLbl,6, 1)
        
        self.grid.addWidget(self.button1, 7, 1)
        self.grid.addWidget(self.button2, 7, 2)
        self.grid.addWidget(self.int2dViewBtn, 8, 1)
        self.grid.addWidget(self.int2dSaveBtn, 8, 2)      
        self.grid.addWidget(self.int2dLogBtn, 8, 3)      


        self.grid.addWidget(self.imposeBtn, 9, 2)  
        self.grid.addWidget(self.imposeSaveBtn, 9, 3)  
        
       

        self.lablePoni4 = QLabel(self)
        self.lablePoni4.setText('Select Poni2: ')
        self.textboxPoni4 = QLineEdit(self)
        self.lableSrc2 = QLabel(self)
        self.lableSrc2.setText('Source Folder2: ')
        self.textboxSrc2 = QLineEdit(self)
        self.lableImgNo1 = QLabel(self)
        self.lableImgNo1.setText('First Image No: ')
        self.textboxImgNo1 = QLineEdit(self)
        self.textboxImgNo1.setValidator(QIntValidator())
        self.buttonPoni3 = QPushButton('Select Poni File 2', self)
        self.buttonPoni3.clicked.connect(self.selectPoni2)

        self.srcBtn2 = QPushButton('Select Source Folder 2', self)
        self.srcBtn2.clicked.connect(self.selectSrcFolder2)

        self.lablePoni4.setStyleSheet(label_css)
        self.lableSrc2.setStyleSheet(label_css)
        self.lableImgNo1.setStyleSheet(label_css)
        self.textboxPoni4.setStyleSheet(text_css)
        self.textboxSrc2.setStyleSheet(text_css)
        self.textboxImgNo1.setStyleSheet(text_css)
        self.buttonPoni3.setStyleSheet( btn_css )

        self.grid.addWidget(self.lablePoni4, 10,1)
        self.grid.addWidget(self.textboxPoni4, 10,2)
        self.grid.addWidget(self.buttonPoni3, 10,3)
        self.grid.addWidget(self.lableSrc2, 11, 1)
        self.grid.addWidget(self.textboxSrc2, 11, 2)
        self.grid.addWidget(self.srcBtn2, 11, 3)
        self.grid.addWidget(self.lableImgNo1, 12, 1)
        self.grid.addWidget(self.textboxImgNo1, 12, 2)


        #self.progress = QProgressBar(self)
        #self.progress.setGeometry(200, 80, 250, 20)
        #self.grid.addWidget(self.progress, 13,2) 

        self.grid.setAlignment(Qt.AlignTop | Qt.AlignLeft)
        self.setLayout(self.grid)
        #self.setFixedSize(500, 400)
        
        #.setsetStyleSheet(label_css)
        #uic.loadUi('layout4_css.ui', self)

        #self.show()

    def openFileNameDialog(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self, "QFileDialog.getOpenFileName()", "",
                                                  "All Files (*.poni);;Python Files (*.py)", options=options)
        if fileName:
            self.textbox4.setText(fileName)

    def openFileNameDialog2(self):
        dirName = str(QFileDialog.getExistingDirectory())
        if dirName:
            self.textbox2.setText(dirName)

    def layout_popup_info(self,info_type,text):
        self.warning=QMessageBox.about(self, info_type.upper(), text)

    def getCheckBoxData(self):
        f = self.flatCB.text()
        d = self.darkCB.text()
        m = self.maskCB.text()
        if not os.path.isfile(f):
            f = ""
        if not os.path.isfile(d):
            d = ""
        if not os.path.isfile(m):
            m = ""
        data = {"flat": f, "dark": d, "mask":m}
        return data

    @pyqtSlot()
    def on_click1(self):
        cbdata = self.getCheckBoxData()
        if(self.textbox1.text()=='' or self.textbox2.text()=='' or self.textbox3.text()=='' or self.textbox4.text()=='' or self.prefixTB.text()==''):
            QMessageBox.about(self, "Error", "All fields are required")
        else:
            h = H5Img(int(self.textbox1.text()),int(self.textbox3.text()),self.textbox2.text(),self.textbox4.text(), prefix=self.prefixTB.text())
            i, okPressed = QInputDialog.getInt(self, "Get Image Plot Number", "Your Image Plot Number", 1, 0, 100, 1)
            if okPressed:
                #plotImgNo = i  # number of images to be plotted
                status = h.Plot(i, h.ai_list, cbdata)
                if status:
                    self.layout_popup_info("Error", "Files does not exist. Kindly check for missing files")
                else:
                    #self.grid.addWidget( status, 4, 0, 3)
                    #self.setLayout( self.grid )
                    print("To add the plot layout to the window")
                   
            else:
                QMessageBox.about(self, "Error", "Please set the number of images to be plotted")

    @pyqtSlot()
    def on_click2(self):
        cbdata = self.getCheckBoxData()
        if(self.textbox1.text()=='' or self.textbox2.text()=='' or self.textbox3.text()=='' or self.textbox4.text()=='' or self.prefixTB.text()==''):
            QMessageBox.about(self, "Error", "All fields are required")
        else:
            h = H5Img(int(self.textbox1.text()),int(self.textbox3.text()),self.textbox2.text(),self.textbox4.text(), prefix=self.prefixTB.text())
            text, okPressed = QInputDialog.getText(self, "Get Sample Name", "Your Sample Name:", QLineEdit.Normal, "")
            if okPressed and text != '':
                h.sample = text #output h5 name
                h.save2HDF(h.source, h.destination, h.sample, cbdata)
                QMessageBox.about(self, "Message", "Sample File Generated Successfully")
            else:
                QMessageBox.about(self, "Error", "Please set the sample file name")


    def int2dViewFunc(self):
        cbdata = self.getCheckBoxData()
        if(self.textbox1.text()=='' or self.textbox2.text()=='' or self.textbox3.text()=='' or self.textbox4.text()=='' or self.prefixTB.text()==''):
            QMessageBox.about(self, "Error", "All fields are required")
        else:
            h = H5Img(int(self.textbox1.text()),int(self.textbox3.text()),self.textbox2.text(),self.textbox4.text(), prefix=self.prefixTB.text())
            i, okPressed = QInputDialog.getInt(self, "Get Image Plot Number", "Your Image Plot Number", 1, 0, 100, 1)
            if okPressed:
                status = h.Plot2d(i, h.ai_list, cbdata)
                if status:
                    self.layout_popup_info("Error", "Files does not exist. Kindly check for missing files")
                else:
                    print("To add the plot layout to the window")
            else:
                QMessageBox.about(self, "Error", "Please set the number of images to be plotted")

    def superImposeSaveFun(self):
        samples = self.superImposeFun()

    def superImposeViewFun(self):
        samples = self.superImposeFun()
        return False ## remove the return False
        h = H5Img(poni_name=self.textboxPoni4.text(), prefix=self.prefixTB.text())
        #h.Plot2dSample(samples)
        mg = MultiGeometry(ais, unit="2th_deg", radial_range=(0,50), azimuth_range=(-180,0))
        I, tth, chi = mg.integrate2d(imgs,500,500)#, radial_range=[10,_max],azimuth_range=[chi_min,chi_max])
        plt.imshow(N.log10(I))
        plt.show()

    def superImposeFun(self):
        samples = ""
        if self.validateFields():
            edf_prefix = "17Oct02WOS_"
            ponifile_1 = self.textbox4.text() 
            ponifile_2 = self.textboxPoni4.text()
            if os.path.isfile('17Mar10D5_0001.edf.gz'):
                mask = fabio.open('17Mar10D5_0001.edf.gz') #'maskWos.edf.gz')
                flat = fabio.open('17Mar10D5_0001.edf.gz') #'flatWOS_18keV_calib9p5_mai.edf.gz')
            #~ dark=fabio.open('17Mar10D5_0079-rsz.edf.gz')
            ai_1 = pyFAI.load(ponifile_1) #ponifile_del25_nu0)
            ai_2 = pyFAI.load(ponifile_1) #ponifile_del25_nu2)
            ais = []
            ais = [ai_1, ai_2]
            edf_folders = [self.textbox2.text() , self.textboxSrc2.text()] #source folder 1 and source folder 2
            samples   = 'S229-30'#['E38_ESRF02_1'] #save as name
            img_deb   = [int(self.textbox1.text()), int(self.textboxImgNo1.text())]#, 1071] #first number of edf1 and first number of edf2
            img_fin   = [int(self.textbox1.text())+int(self.textbox3.text()), int(self.textboxImgNo1.text())+int(self.textbox3.text())]#, 1251] #last number of edf1 and last number of edf2

            destination = ""
            prefix = self.prefixTB.text()
            extract_all(destination, edf_folders, samples, img_deb, img_fin, ais, prefix)
        return samples

    def validateFields(self):
        poni2 = self.textboxPoni4.text()
        src2 = self.textboxSrc2.text() 
        imgno2 = self.textboxImgNo1.text()
        prefix = self.prefixTB.text()
        if len(prefix)==0:
            self.layout_popup_info("Error", "Kindly provide prefix text")
            return False
        if len(poni2)==0:
            self.layout_popup_info("Error", "Kindly provide Poni 2 file")
            return False
        if len(src2)==0:
            self.layout_popup_info("Error", "Kindly provide source 2 file")
            return False
        if len(imgno2)==0:
            self.layout_popup_info("Error", "Kindly provide image number 2 value")
            return False
        
        return True

    def int2dSaveFun(self):
        cbdata = self.getCheckBoxData()
        if(self.textbox1.text()=='' or self.textbox2.text()=='' or self.textbox3.text()=='' or self.textbox4.text()=='' or self.prefixTB.text()==''):
            QMessageBox.about(self, "Error", "All fields are required")
        else:
            h = H5Img(int(self.textbox1.text()),int(self.textbox3.text()),self.textbox2.text(),self.textbox4.text(), prefix=self.prefixTB.text())
            text, okPressed = QInputDialog.getText(self, "Get Sample Name", "Your Sample Name:", QLineEdit.Normal, "")
            if okPressed and text != '':
                h.sample = text #output h5 name
                h.save2HDF2d(h.source, h.destination, h.sample, cbdata)
                QMessageBox.about(self, "Message", "Sample File Generated Successfully")
            else:
                QMessageBox.about(self, "Error", "Please set the sample file name")

    def int2dLogFun(self):
        cbdata = self.getCheckBoxData()
        if(self.textbox1.text()=='' or self.textbox2.text()=='' or self.textbox3.text()=='' or self.textbox4.text()=='' or self.prefixTB.text()==''):
            QMessageBox.about(self, "Error", "All fields are required")
        else:
            h = H5Img(int(self.textbox1.text()),int(self.textbox3.text()),self.textbox2.text(),self.textbox4.text(), prefix=self.prefixTB.text())
            i, okPressed = QInputDialog.getInt(self, "Get Image Plot Number", "Your Image Plot Number", 1, 0, 100, 1)
            if okPressed:
                status = h.Plot2dLog(i, h.ai_list, cbdata)
                if status:
                    self.layout_popup_info("Error", "Files does not exist. Kindly check for missing files")
                else:
                    print("To add the plot layout to the window")
            else:
                QMessageBox.about(self, "Error", "Please set the number of images to be plotted")

    def selectPoni2(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self, "QFileDialog.getOpenFileName()", "",
                                                  "All Files (*.poni);;Python Files (*.py)", options=options)
        if fileName:
            self.textboxPoni4 .setText(fileName)

    @pyqtSlot()
    def on_click3(self):
        self.openFileNameDialog()

    def flatFuncSel(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self, "QFileDialog.getOpenFileName()", "",
                                                  "All EDF GZ File (*.edf.gz);;All EDF Files (*.edf)", options=options)
        if fileName:
            self.flatCB.setText(fileName)

    def darkFuncSel(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self, "QFileDialog.getOpenFileName()", "",
                                                  "All EDF GZ File (*.edf.gz);;All EDF Files (*.edf)", options=options)
        if fileName:
            self.darkCB.setText(fileName)

    def maskFuncSel(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self, "QFileDialog.getOpenFileName()", "",
                                                  "All EDF GZ File (*.edf.gz);;All EDF Files (*.edf)", options=options)
        if fileName:
            self.maskCB.setText(fileName)

    def selectSrcFolder2(self):
        dirName = str(QFileDialog.getExistingDirectory())
        if dirName:
            self.textboxSrc2.setText(dirName)

    @pyqtSlot()
    def on_click4(self):
        self.openFileNameDialog2()

if __name__ == '__main1__':
    app = QApplication(sys.argv)
    ex = App()
    sys.exit(app.exec_())
