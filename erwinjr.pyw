#===============================================================================
# ErwinJr is a simulation program for quantum semiconductor lasers.
# Copyright (C) 2012 Kale J. Franz, PhD
#
# A portion of this code is Copyright (c) 2011, California Institute of 
# Technology ("Caltech"). U.S. Government sponsorship acknowledged.
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#===============================================================================


from __future__ import division
import os, sys
from PyQt4.QtCore import *
from PyQt4.QtGui import *
import PyQt4.Qwt5 as Qwt
from numpy import *
from functools import partial
import time

import matplotlib
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D

import settings
import SupportClasses
import ThePhysics

#===============================================================================
# Version
#===============================================================================
ejVersion = 120217
majorVersion = '2.1.0'

#===============================================================================
# Global Variables
#===============================================================================
pi = 3.14159265
e0 = 1.60217653e-19;  #electron charge
eps0 = 8.85e-12;
m0 = 9.10938188e-31;   #free electron mass (kg)
hbar = 6.6260693e-34/(2*pi); #Planck's constant (J s)
kb = 1.386505e-23 / e0; #eV/K


class MainWindow(QMainWindow):
    def __init__(self, fileName=None, parent=None):
        super(MainWindow, self).__init__(parent)
        
        self.colors = [(149,115,179), (110,124,190), (147,177,132), (174,199,82), 
                       (128,128,130), (218,189,63), (223,155,74), (210,87,71), 
                       (185,82,159), (105,189,69), (20,20,20), (110,205,222), 
                       (57,82,164)]
                       
        if os.name == "nt":
            self.newLine = '\n'
        elif os.name == "posix":
            self.newLine = '\r\n'
        
        self.qclayers = ThePhysics.QCLayers()
        self.strata = ThePhysics.Strata()
        
        self.stratumMaterialsList = ['Active Core', 
                                     'InP',
                                     'GaAs',
                                     'InGaAs', 
                                     'InAlAs', 
                                     'Au', 
                                     'SiNx', 
                                     'SiO2', 
                                     'Air']
        self.waveguideFacetsList = ['as-cleaved + as-cleaved',
                                    'as-cleaved + perfect HR',
                                    'as-cleaved + perfect AR',
                                    'perfect AR + perfect HR',
                                    'custom coating + as-cleaved',
                                    'custom coating + perfect HR',
                                    'custom coating + perfect AR']
        self.substratesList = ['InP', 'GaAs', 'GaSb', 'GaN']

        self.filename = fileName
        self.plotDirty = False
        self.solveType = None
        
        self.plotVX = False
        self.plotVL = False
        self.plotLH = False
        self.plotSO = False
        
        self.stateHolder = []

        self.create_Quantum_menu()
        self.create_main_frame()
        
        self.create_zoomer()

        self.update_inputBoxes()
        
        self.input_substrate('InP')
        self.layerTable_refresh()
        self.layerTable.selectRow(1)
        self.layerTable.setFocus()
        
        self.stratumTable_refresh()
        
        self.update_inputBoxes()
        self.update_stratum_inputBoxes()
        
        qsettings = QSettings()
        self.recentFiles = qsettings.value("RecentFiles").toStringList()
        self.restoreGeometry(
                qsettings.value("MainWindow/Geometry").toByteArray())
        self.restoreState(qsettings.value("MainWindow/State").toByteArray())
        self.updateFileMenu()

        self.dirty = False

        if self.filename:
            self.fileOpen(self.filename)
        else:
            QTimer.singleShot(0, self.loadInitialFile)
            
        self.dirty = False
        self.update_windowTitle()

        

        
        
#===============================================================================
# Create Main Frame        
#===============================================================================

    def create_main_frame(self):
        
        self.mainTabWidget = QTabWidget()
        
        # ##########################
        #
        # Quantum Tab
        #
        # ##########################
        
        #set up quantumCanvas for band structure plot
        self.quantumCanvas = Qwt.QwtPlot(self)
        self.quantumCanvas.setCanvasBackground(Qt.white)
        self.quantumCanvas.canvas().setCursor(Qt.ArrowCursor)
        
        #set up layerTable
        self.layerTable = QTableWidget()
        self.layerTable.setSelectionBehavior(QTableWidget.SelectRows)
        self.layerTable.setSelectionMode(QTableWidget.SingleSelection)
        self.layerTable.setMaximumWidth(340)
        self.layerTable.setMinimumWidth(340)
        self.connect(self.layerTable,SIGNAL("itemChanged(QTableWidgetItem*)"),self.layerTable_itemChanged)
        self.connect(self.layerTable,SIGNAL("itemSelectionChanged()"),self.layerTable_itemSelectionChanged)
        
        #set up buttons
        deleteLayerButton = QPushButton("Delete Layer")
        self.connect(deleteLayerButton, SIGNAL("clicked()"), self.delete_layer)
        insertLayerAboveButton = QPushButton("Insert Layer Above")
        self.connect(insertLayerAboveButton, SIGNAL("clicked()"), self.insert_layerAbove)
        self.goButton = QPushButton("Solve Whole")
        self.connect(self.goButton, SIGNAL("clicked()"), self.solve_whole)
        self.solveButton = QPushButton("Solve Basis")
        self.connect(self.solveButton, SIGNAL("clicked()"), self.solve)
        
        #set up left inputs
        self.substrateBox = QComboBox()
        self.substrateBox.addItems(self.substratesList)
        self.connect(self.substrateBox,SIGNAL("currentIndexChanged(const QString)"), self.input_substrate)
        
        inputEFieldLabel = QLabel('<center><b><i>E<sub>field</sub></i></b></center>')
        self.inputEFieldBox = QDoubleSpinBox()
        self.inputEFieldBox.setDecimals(1)
        self.inputEFieldBox.setSuffix(' kV/cm')
        self.inputEFieldBox.setRange(0.0,250.0)
        self.connect(self.inputEFieldBox, SIGNAL("valueChanged(double)"), self.input_EField)
        
        inputHorzResLabel = QLabel('<center><b>Horizontal<br>Resolution</b></center>')
        self.inputHorzResBox = QComboBox();
        self.inputHorzResBox.addItems(['1.0','0.5','0.25','0.2','0.1'])
        self.connect(self.inputHorzResBox, SIGNAL("currentIndexChanged(int)"), self.input_horzRes)
        
        inputVertResLabel = QLabel('<center><b>Vertical<br>Resolution</b></center>')
        self.inputVertResBox = QDoubleSpinBox()
        self.inputVertResBox.setDecimals(2)
        self.inputVertResBox.setValue(0.5)
        self.inputVertResBox.setRange(0.0,10.0)
        self.inputVertResBox.setSingleStep(0.1)
        self.inputVertResBox.setSuffix(' meV')
        self.connect(self.inputVertResBox, SIGNAL("valueChanged(double)"), self.input_vertRes)
        
        inputRepeatsLabel = QLabel('<center><b>Structure Repeats</b></center>')
        self.inputRepeatsBox = QSpinBox()
        self.inputRepeatsBox.setValue(1)
        self.inputRepeatsBox.setRange(1,5)
        self.connect(self.inputRepeatsBox, SIGNAL("valueChanged(int)"), self.input_repeats)
        
        self.inputARInjectorCheck = QCheckBox("AR->Injector")
        self.inputInjectorARCheck = QCheckBox("Injector->AR")
        self.inputARInjectorCheck.setChecked(True)
        self.inputInjectorARCheck.setChecked(True)
        basis_groupBox = QGroupBox("Basis Divisions")
        vbox = QVBoxLayout()
        vbox.addWidget(self.inputARInjectorCheck)
        vbox.addWidget(self.inputInjectorARCheck)
        basis_groupBox.setLayout(vbox)
        self.connect(self.inputARInjectorCheck, SIGNAL("stateChanged(int)"), self.input_basis)
        self.connect(self.inputInjectorARCheck, SIGNAL("stateChanged(int)"), self.input_basis)
        
        # Lp groupbox
        self.LpFirstSpinbox = QSpinBox()
        self.LpFirstSpinbox.setValue(1)
        self.LpFirstSpinbox.setRange(1,1)
        self.connect(self.LpFirstSpinbox, SIGNAL("valueChanged(int)"), self.update_inputBoxes)
        self.LpLastSpinbox  = QSpinBox()
        self.LpLastSpinbox.setValue(1)
        self.LpLastSpinbox.setRange(1,1)
        self.connect(self.LpLastSpinbox, SIGNAL("valueChanged(int)"), self.update_inputBoxes)
        self.LpStringBox = QTextEdit('')
        self.LpStringBox.setReadOnly(True)
        self.LpStringBox.setSizePolicy(QSizePolicy(QSizePolicy.Fixed,QSizePolicy.Fixed))
        self.LpStringBox.setMaximumHeight(62)
        self.LpStringBox.setMaximumWidth(97)
        LpLayout = QGridLayout()
        LpLayout.addWidget(QLabel('<b>first</b>'), 0,0)
        LpLayout.addWidget(QLabel('<b>last</b>'), 0,1)
        LpLayout.addWidget(self.LpFirstSpinbox, 1,0)
        LpLayout.addWidget(self.LpLastSpinbox, 1,1)
        LpLayout.addWidget(self.LpStringBox, 2,0, 1,2)
        LpLayout_groupBox = QGroupBox("Period Info")
        LpLayout_groupBox.setLayout(LpLayout)
        
        #set up material composition inputs
        self.mtrl_title   = QLabel('<center><b>Mole Fractions</b></center')
        self.mtrl_header1    = QLabel('<center><b>In<sub>x</sub>Ga<sub>1-x</sub>As</b></center')
        self.mtrl_header2    = QLabel('<center><b>Al<sub>1-x</sub>In<sub>x</sub>As</b></center')
        self.mtrl_header3    = QLabel('<center><b>Offset</b></center')
        self.mtrl_row1 = QLabel('<center><b>#1</b></center')
        self.mtrl_row2 = QLabel('<center><b>#2</b></center')
        self.mtrl_row3 = QLabel('<center><b>#3</b></center')
        self.mtrl_row4 = QLabel('<center><b>#4</b></center')
        self.inputMoleFrac1Box = QDoubleSpinBox()
        self.inputMoleFrac1Box.setDecimals(3)
        self.inputMoleFrac1Box.setValue(0.53)
        self.inputMoleFrac1Box.setRange(0.0,1.0)
        self.inputMoleFrac1Box.setSingleStep(0.001)
        self.connect(self.inputMoleFrac1Box, SIGNAL("editingFinished()"), partial(self.input_moleFrac,1))
        self.inputMoleFrac2Box = QDoubleSpinBox()
        self.inputMoleFrac2Box.setDecimals(3)
        self.inputMoleFrac2Box.setValue(0.52)
        self.inputMoleFrac2Box.setRange(0.0,1.0)
        self.inputMoleFrac2Box.setSingleStep(0.001)
        self.connect(self.inputMoleFrac2Box, SIGNAL("editingFinished()"), partial(self.input_moleFrac,2))
        self.inputMoleFrac3Box = QDoubleSpinBox()
        self.inputMoleFrac3Box.setDecimals(3)
        self.inputMoleFrac3Box.setValue(0.53)
        self.inputMoleFrac3Box.setRange(0.0,1.0)
        self.inputMoleFrac3Box.setSingleStep(0.001)
        self.connect(self.inputMoleFrac3Box, SIGNAL("editingFinished()"), partial(self.input_moleFrac,3))
        self.inputMoleFrac4Box = QDoubleSpinBox()
        self.inputMoleFrac4Box.setDecimals(3)
        self.inputMoleFrac4Box.setValue(0.52)
        self.inputMoleFrac4Box.setRange(0.0,1.0)
        self.inputMoleFrac4Box.setSingleStep(0.001)
        self.connect(self.inputMoleFrac4Box, SIGNAL("editingFinished()"), partial(self.input_moleFrac,4))
        self.inputMoleFrac5Box = QDoubleSpinBox()
        self.inputMoleFrac5Box.setDecimals(3)
        self.inputMoleFrac5Box.setValue(0.53)
        self.inputMoleFrac5Box.setRange(0.0,1.0)
        self.inputMoleFrac5Box.setSingleStep(0.001)
        self.connect(self.inputMoleFrac5Box, SIGNAL("editingFinished()"), partial(self.input_moleFrac,5))
        self.inputMoleFrac6Box = QDoubleSpinBox()
        self.inputMoleFrac6Box.setDecimals(3)
        self.inputMoleFrac6Box.setValue(0.52)
        self.inputMoleFrac6Box.setRange(0.0,1.0)
        self.inputMoleFrac6Box.setSingleStep(0.001)
        self.connect(self.inputMoleFrac6Box, SIGNAL("editingFinished()"), partial(self.input_moleFrac,6))
        self.inputMoleFrac7Box = QDoubleSpinBox()
        self.inputMoleFrac7Box.setDecimals(3)
        self.inputMoleFrac7Box.setValue(0.53)
        self.inputMoleFrac7Box.setRange(0.0,1.0)
        self.inputMoleFrac7Box.setSingleStep(0.001)
        self.connect(self.inputMoleFrac7Box, SIGNAL("valueChanged(double)"), partial(self.input_moleFrac,7))
        self.inputMoleFrac8Box = QDoubleSpinBox()
        self.inputMoleFrac8Box.setDecimals(3)
        self.inputMoleFrac8Box.setValue(0.52)
        self.inputMoleFrac8Box.setRange(0.0,1.0)
        self.inputMoleFrac8Box.setSingleStep(0.001)
        self.connect(self.inputMoleFrac8Box, SIGNAL("valueChanged(double)"), partial(self.input_moleFrac,8))
        self.offset1Label = QLabel('')
        self.offset2Label = QLabel('')
        self.offset3Label = QLabel('')
        self.offset4Label = QLabel('')
        self.strainDescription = QLabel('')
        #self.strainDescription.setTextAlignment(Qt.AlignHCenter)
        mtrl_grid = QGridLayout()
        mtrl_grid.addWidget(self.mtrl_title, 0,0, 1,4)
        mtrl_grid.addWidget(self.mtrl_header1, 1,1)
        mtrl_grid.addWidget(self.mtrl_header2, 1,2)
        mtrl_grid.addWidget(self.mtrl_header3, 1,3)
        mtrl_grid.addWidget(self.mtrl_row1, 2,0)
        mtrl_grid.addWidget(self.inputMoleFrac1Box, 2,1)
        mtrl_grid.addWidget(self.inputMoleFrac2Box, 2,2)
        mtrl_grid.addWidget(self.offset1Label, 2,3)
        mtrl_grid.addWidget(self.mtrl_row2, 3,0)
        mtrl_grid.addWidget(self.inputMoleFrac3Box, 3,1)
        mtrl_grid.addWidget(self.inputMoleFrac4Box, 3,2)
        mtrl_grid.addWidget(self.offset2Label, 3,3)
        mtrl_grid.addWidget(self.mtrl_row3, 4,0)
        mtrl_grid.addWidget(self.inputMoleFrac5Box, 4,1)
        mtrl_grid.addWidget(self.inputMoleFrac6Box, 4,2)
        mtrl_grid.addWidget(self.offset3Label, 4,3)
        mtrl_grid.addWidget(self.mtrl_row4, 5,0)
        mtrl_grid.addWidget(self.inputMoleFrac7Box, 5,1)
        mtrl_grid.addWidget(self.inputMoleFrac8Box, 5,2)
        mtrl_grid.addWidget(self.offset4Label, 5,3)
        mtrl_grid.addWidget(self.strainDescription, 6,0, 1,4)
        self.mtrl_groupBox = QGroupBox()
        self.mtrl_groupBox.setLayout(mtrl_grid)
        
        #set up description box
        self.DescriptionBox = QTextEdit('')
        self.DescriptionBox.setReadOnly(False)
        self.DescriptionBox.setSizePolicy(QSizePolicy(QSizePolicy.Fixed,QSizePolicy.Fixed))
        self.DescriptionBox.setMaximumHeight(40)
        self.DescriptionBox.setMaximumWidth(190)
        self.connect(self.DescriptionBox, SIGNAL("textChanged()"), self.input_description)
        DescLayout = QVBoxLayout()
        DescLayout.addWidget(self.DescriptionBox)
        DescLayout_groupBox = QGroupBox("Description")
        DescLayout_groupBox.setLayout(DescLayout)

        #set up plot control inputs
        self.zoomButton = QPushButton("Zoom")
        self.zoomButton.setCheckable(True)
        self.connect(self.zoomButton, SIGNAL("toggled(bool)"), self.zoom)
        zoomOutButton = QPushButton("Zoom Out")
        self.connect(zoomOutButton, SIGNAL("clicked()"), self.zoom_out)
        self.panButton = QPushButton("Pan")
        self.panButton.setCheckable(True)
        self.connect(self.panButton, SIGNAL("toggled(bool)"), self.pan)
        self.wellSelectButton = QPushButton("Layer Select")
        self.wellSelectButton.setCheckable(True)
        self.connect(self.wellSelectButton, SIGNAL("toggled(bool)"), self.well_select)
        clearWFsButton = QPushButton("Clear Wavefunctions")
        self.connect(clearWFsButton, SIGNAL("clicked()"), self.clear_WFs)
        plotControl_grid = QGridLayout()
        plotControl_grid.addWidget(self.wellSelectButton, 0,0, 1,2)
        plotControl_grid.addWidget(self.zoomButton, 1,0, 1,1)
        plotControl_grid.addWidget(zoomOutButton, 1,1, 1,1)
        plotControl_grid.addWidget(self.panButton, 2,0, 1,1)
        plotControl_grid.addWidget(clearWFsButton, 2,1, 1,1)
        plotControl_groupBox = QGroupBox("Plot Controls")
        plotControl_groupBox.setLayout(plotControl_grid)
        
        #set up Calculate controls
        self.pairSelectButton = QPushButton("Pair Select")
        self.pairSelectButton.setCheckable(True)
        self.connect(self.pairSelectButton, SIGNAL("toggled(bool)"), self.pair_select)
        self.FoMButton = QPushButton("Figure of Merit")
        self.FoMButton.setEnabled(False)
        self.connect(self.FoMButton, SIGNAL("clicked()"), self.figure_of_merit)
        self.transferOpticalParametersButton = QPushButton("-> Optical Params")
        self.transferOpticalParametersButton.setEnabled(False)
        self.connect(self.transferOpticalParametersButton, SIGNAL("clicked()"), self.transfer_optical_parameters)
        self.pairSelectString = QTextEdit('')
        self.pairSelectString.setReadOnly(True)
        self.pairSelectString.setSizePolicy(QSizePolicy(QSizePolicy.Fixed,QSizePolicy.Fixed))
        self.pairSelectString.setMaximumHeight(130)
        self.pairSelectString.setMaximumWidth(190)
        calculateControl_grid = QGridLayout()
        calculateControl_grid.addWidget(self.pairSelectButton, 0,0, 1,2)
        calculateControl_grid.addWidget(self.FoMButton, 1,0, 1,1)
        calculateControl_grid.addWidget(self.transferOpticalParametersButton, 1,1, 1,1)
        calculateControl_grid.addWidget(self.pairSelectString, 2,0, 1,2)
        calculateControl_groupBox = QGroupBox("Calculate")
        calculateControl_groupBox.setLayout(calculateControl_grid)
       
        #lay out GUI
        #vBox1
        vBox1 = QVBoxLayout()
        vBox1.addWidget(QLabel("<center><b>Substrate</b></center>"))
        vBox1.addWidget(self.substrateBox)
        vBox1.addWidget(inputEFieldLabel)
        vBox1.addWidget(self.inputEFieldBox)
        vBox1.addWidget(inputHorzResLabel)
        vBox1.addWidget(self.inputHorzResBox)
        vBox1.addWidget(inputVertResLabel)
        vBox1.addWidget(self.inputVertResBox)
        vBox1.addWidget(inputRepeatsLabel)
        vBox1.addWidget(self.inputRepeatsBox)
        vBox1.addWidget(basis_groupBox)
        vBox1.addWidget(LpLayout_groupBox)
        #vBox1.addWidget(designBy_groupBox)
        vBox1.addStretch()
        
        #vBox2
        vBox2 = QGridLayout()
        vBox2.addWidget(insertLayerAboveButton, 0,0)
        vBox2.addWidget(deleteLayerButton, 0,1)        
        vBox2.addWidget(self.layerTable, 1,0, 1,2)
        #vBox2.addWidget(updateButton)      
        
        #vBox3
        vBox3 = QVBoxLayout()
        vBox3.addWidget(self.solveButton)
        vBox3.addWidget(self.goButton)
        vBox3.addWidget(DescLayout_groupBox)
        vBox3.addWidget(self.mtrl_groupBox)
        vBox3.addWidget(plotControl_groupBox)
        vBox3.addWidget(calculateControl_groupBox)
        vBox3.addStretch()
                
        #vBox4
        vBox4 = QVBoxLayout()
        vBox4.addWidget(self.quantumCanvas)
        
        quantumLayout = QHBoxLayout()
        quantumLayout.addLayout(vBox1)
        quantumLayout.addLayout(vBox2)
        quantumLayout.addLayout(vBox3)
        quantumLayout.addLayout(vBox4)
        
        self.quantumWidget = QWidget()
        self.quantumWidget.setLayout(quantumLayout)
        self.quantumWidget.setAutoFillBackground(True)
        self.quantumWidget.setBackgroundRole(QPalette.Window)
        

        # ##########################
        #
        # Optical Tab
        #
        # ##########################
        
        #vBox1
        vBox1Grid = QGridLayout()
        
        self.editOpticalParametersBox = QCheckBox('Edit Parameters')
        self.editOpticalParametersBox.setChecked(False)
        self.connect(self.editOpticalParametersBox,SIGNAL('toggled(bool)'),self.edit_optical_parameters)
        vBox1Grid.addWidget(self.editOpticalParametersBox, 0,0, 1,2, Qt.AlignCenter)
        
        wlLabel = QLabel('<b><center>Wavelength</center</b>')
        vBox1Grid.addWidget(wlLabel, 1,0, 1,2)
        self.wavelengthBox = QDoubleSpinBox()
        self.wavelengthBox.setValue(self.strata.wavelength)
        self.wavelengthBox.setSuffix(u' \u03BCm')
        self.wavelengthBox.setDecimals(3)
        self.wavelengthBox.setSingleStep(0.1)
        self.wavelengthBox.setRange(0.0,30.0)
        self.wavelengthBox.setReadOnly(True)
        self.wavelengthBox.setStyleSheet('color:gray')
        self.connect(self.wavelengthBox, SIGNAL("valueChanged(double)"), self.input_wavelength)
        vBox1Grid.addWidget(self.wavelengthBox, 2,0, 1,2)
        
        vBox1Grid.addItem(QSpacerItem(20,20), 3,0, 1,2)
        
        vBox1Grid.addWidget(QLabel('<b><center>Operating<br>Field</center></b>'), 4,0, 1,1)
        self.operatingFieldBox = QDoubleSpinBox()
        self.operatingFieldBox.setDecimals(1)
        self.operatingFieldBox.setRange(0.,300.)
        self.operatingFieldBox.setSingleStep(1)
        self.operatingFieldBox.setSuffix(u' kV/cm')
        self.operatingFieldBox.setReadOnly(True)
        self.operatingFieldBox.setStyleSheet('color:gray')
        self.connect(self.operatingFieldBox, SIGNAL("valueChanged(double)"), self.input_operatingField)
        vBox1Grid.addWidget(self.operatingFieldBox, 5,0, 1,1)
        
        vBox1Grid.addWidget(QLabel('<b><center>Active Core<br>Period Length</center></b>'), 4,1, 1,1)
        self.ACPeriodLengthBox = QDoubleSpinBox()
        self.ACPeriodLengthBox.setDecimals(1)
        self.ACPeriodLengthBox.setRange(0.,10000.)
        self.ACPeriodLengthBox.setSingleStep(1)
        self.ACPeriodLengthBox.setSuffix(u' \u212B')
        self.connect(self.ACPeriodLengthBox, SIGNAL("valueChanged(double)"), self.input_ACPeriodLength)
        self.ACPeriodLengthBox.setReadOnly(True)
        self.ACPeriodLengthBox.setStyleSheet('color:gray')
        vBox1Grid.addWidget(self.ACPeriodLengthBox, 5,1, 1,1)
        
        vBox1Grid.addWidget(QLabel('<b><center>Active Core<br>Periods</center></b>'), 6,0, 1,1)
        self.ACPeriodsBox = QSpinBox()
        self.ACPeriodsBox.setRange(1,99)
        self.connect(self.ACPeriodsBox, SIGNAL("valueChanged(int)"), self.input_ACPeriods)
        vBox1Grid.addWidget(self.ACPeriodsBox, 7,0, 1,1)

        vBox1Grid.addWidget(QLabel('<b><center>Operating<br>Voltage</center></b>'), 6,1, 1,1)
        self.OperatingVoltageBox = QLineEdit()
        self.OperatingVoltageBox.setMaximumWidth(150)
        self.OperatingVoltageBox.setReadOnly(True)
        self.OperatingVoltageBox.setStyleSheet('color:gray')
        vBox1Grid.addWidget(self.OperatingVoltageBox, 7,1, 1,1)
        
        vBox1Grid.addItem(QSpacerItem(20,20), 8,0, 1,2)
        
        vBox1Grid.addWidget(QLabel(u'<b><center><i>\u03B1<sub>core</sub></i></center></b>'), 9,0, 1,1)
        self.aCoreBox = QLineEdit()
        self.aCoreBox.setMaximumWidth(150)
        self.aCoreBox.setReadOnly(True)
        self.aCoreBox.setStyleSheet('color:gray')
        self.connect(self.aCoreBox, SIGNAL("editingFinished()"), self.input_aCore)
        vBox1Grid.addWidget(self.aCoreBox, 10,0, 1,1)
        
        vBox1Grid.addWidget(QLabel(u'<center><b>\u00F1<i><sub>core</sub></i></b></center>'), 9,1, 1,1)
        self.nCoreBox = QLineEdit()
        self.nCoreBox.setMaximumWidth(150)
        self.nCoreBox.setReadOnly(True)
        self.nCoreBox.setStyleSheet('color:gray')
        vBox1Grid.addWidget(self.nCoreBox, 10,1, 1,1)
        
        vBox1Grid.addItem(QSpacerItem(20,20), 11,0, 1,2)
        
        vBox1Grid.addWidget(QLabel(u'<center><b>Transition Broadening</b></center>'), 12,0, 1,2)
        self.transitionBroadeningBox = QDoubleSpinBox()
        self.transitionBroadeningBox.setMaximumWidth(85)
        self.transitionBroadeningBox.setDecimals(1)
        self.transitionBroadeningBox.setRange(0.,1000.)
        self.transitionBroadeningBox.setSingleStep(1)
        self.transitionBroadeningBox.setSuffix(u' meV')
        self.transitionBroadeningBox.setReadOnly(True)
        self.transitionBroadeningBox.setStyleSheet('color:gray')
        self.connect(self.transitionBroadeningBox, SIGNAL("valueChanged(double)"), self.input_transitionBroadening)
        vBox1Grid.addWidget(self.transitionBroadeningBox, 13,0, 1,2, Qt.AlignCenter)
        
        vBox1subGrid1 = QGridLayout()
        
        vBox1subGrid1.addWidget(QLabel(u'<b><center><i>\u03C4<sub>upper</sub></i></center></b>'), 0,0, 1,1)
        self.tauUpperBox = QDoubleSpinBox()
        self.tauUpperBox.setDecimals(3)
        self.tauUpperBox.setRange(0.,99.)
        self.tauUpperBox.setSingleStep(1)
        self.tauUpperBox.setSuffix(u' ps')
        self.tauUpperBox.setReadOnly(True)
        self.tauUpperBox.setStyleSheet('color:gray')
        #self.tauUpperBox.setSizePolicy(QSizePolicy(QSizePolicy.Fixed,QSizePolicy.Fixed))
        self.connect(self.tauUpperBox, SIGNAL("valueChanged(double)"), self.input_tauUpper)
        vBox1subGrid1.addWidget(self.tauUpperBox, 1,0, 1,1)
        
        vBox1subGrid1.addWidget(QLabel(u'<b><center><i>\u03C4<sub>lower</sub></i></b></center></b>'), 0,1, 1,1)
        self.tauLowerBox = QDoubleSpinBox()
        self.tauLowerBox.setDecimals(3)
        self.tauLowerBox.setRange(0.,99.)
        self.tauLowerBox.setSingleStep(1)
        self.tauLowerBox.setSuffix(u' ps')
        self.tauLowerBox.setReadOnly(True)
        self.tauLowerBox.setStyleSheet('color:gray')
        #self.tauLowerBox.setSizePolicy(QSizePolicy(QSizePolicy.Fixed,QSizePolicy.Fixed))
        self.connect(self.tauLowerBox, SIGNAL("valueChanged(double)"), self.input_tauLower)
        vBox1subGrid1.addWidget(self.tauLowerBox, 1,1, 1,1)
        
        vBox1subGrid1.addWidget(QLabel(u'<b><center><i>\u03C4<sub>upper,lower</sub></i></b></center></b>'), 0,2, 1,1)
        self.tauUpperLowerBox = QDoubleSpinBox()
        self.tauUpperLowerBox.setDecimals(3)
        self.tauUpperLowerBox.setRange(0.,99.)
        self.tauUpperLowerBox.setSingleStep(1)
        self.tauUpperLowerBox.setSuffix(u' ps')
        self.tauUpperLowerBox.setReadOnly(True)
        self.tauUpperLowerBox.setStyleSheet('color:gray')
        #self.tauUpperLowerBox.setSizePolicy(QSizePolicy(QSizePolicy.Fixed,QSizePolicy.Fixed))
        self.connect(self.tauUpperLowerBox, SIGNAL("valueChanged(double)"), self.input_tauUpperLower)
        vBox1subGrid1.addWidget(self.tauUpperLowerBox, 1,2, 1,1)
        
        vBox1Grid.addLayout(vBox1subGrid1, 14,0, 1,2)
        
        vBox1Grid.addWidget(QLabel(u'<b><center>optical dipole</center></b>'), 15,0, 1,1)
        self.opticalDipoleBox = QDoubleSpinBox()
        self.opticalDipoleBox.setDecimals(1)
        self.opticalDipoleBox.setRange(0.,10000.)
        self.opticalDipoleBox.setSingleStep(1)
        self.opticalDipoleBox.setSuffix(u' \u212B')
        self.opticalDipoleBox.setReadOnly(True)
        self.opticalDipoleBox.setStyleSheet('color:gray')
        self.connect(self.opticalDipoleBox, SIGNAL("valueChanged(double)"), self.input_opticalDipole)
        vBox1Grid.addWidget(self.opticalDipoleBox, 16,0, 1,1)
        
        vBox1Grid.addWidget(QLabel(u'<center><b>Figure of Merit</b></center>'), 15,1, 1,1)
        self.FoMBox = QLineEdit()
        self.FoMBox.setMaximumWidth(150)
        self.FoMBox.setReadOnly(True)
        self.FoMBox.setStyleSheet('color:gray')
        vBox1Grid.addWidget(self.FoMBox, 16,1, 1,1)
        
        vBox1Grid.addItem(QSpacerItem(20,20), 17,0, 1,2)
        
        vBox1Grid.addWidget(QLabel(u'<center><b>Waveguide Facets</b></center>'), 18,0, 1,2)
        self.waveguideFacetsBox = QComboBox();
        self.waveguideFacetsBox.addItems(self.waveguideFacetsList)
        self.connect(self.waveguideFacetsBox, SIGNAL("currentIndexChanged(const QString &)"), self.input_waveguideFacets)
        vBox1Grid.addWidget(self.waveguideFacetsBox, 19,0, 1,2)
        
        vBox1Grid.addWidget(QLabel(u'<b><center>Waveguide<br>Length</center></b>'), 20,0, 1,1)
        self.waveguideLengthBox = QDoubleSpinBox()
        self.waveguideLengthBox.setDecimals(1)
        self.waveguideLengthBox.setRange(0.,20.)
        self.waveguideLengthBox.setSingleStep(1)
        self.waveguideLengthBox.setSuffix(u' mm')
        self.connect(self.waveguideLengthBox, SIGNAL("valueChanged(double)"), self.input_waveguideLength)
        vBox1Grid.addWidget(self.waveguideLengthBox, 21,0, 1,1)
        
        
        self.customFacetBoxLabel = QLabel(u'<b><center>Custom<br>Reflectivity</center></b>')
        self.customFacetBoxLabel.setStyleSheet('color:gray')
        vBox1Grid.addWidget(self.customFacetBoxLabel, 20,1, 1,1)
        self.customFacetBox = QDoubleSpinBox()
        self.customFacetBox.setDecimals(1)
        self.customFacetBox.setRange(0.,100.)
        self.customFacetBox.setSingleStep(1)
        self.customFacetBox.setSuffix(u'%')
        self.customFacetBox.setEnabled(False)
        self.connect(self.customFacetBox, SIGNAL("valueChanged(double)"), self.input_customFacet)
        vBox1Grid.addWidget(self.customFacetBox, 21,1, 1,1)
        
        
        vBox1GridWidget = QWidget()
        vBox1GridWidget.setLayout(vBox1Grid)
        vBox1GridWidget.setContentsMargins(0,0,0,0)
        vBox1GridWidget.setMaximumWidth(235)
        vBox1 = QVBoxLayout()
        vBox1.addWidget(vBox1GridWidget)
        vBox1.setSpacing(0)
        vBox1.addStretch()
        
        
        #set up stratumTable
        self.stratumTable = QTableWidget()
        self.stratumTable.setSelectionBehavior(QTableWidget.SelectRows)
        self.stratumTable.setSelectionMode(QTableWidget.SingleSelection)
        self.stratumTable.setMinimumWidth(450)
        self.stratumTable.setMaximumWidth(450)
        self.stratumTable.setMinimumHeight(450)
        self.stratumTable.setMaximumHeight(650)
        self.connect(self.stratumTable,SIGNAL("itemChanged(QTableWidgetItem*)"),self.stratumTable_itemChanged)
        self.connect(self.stratumTable,SIGNAL("itemSelectionChanged()"),self.stratumTable_itemSelectionChanged)
        

        insertStratumAboveButton = QPushButton("Insert Stratum Above")
        self.connect(insertStratumAboveButton, SIGNAL("clicked()"), self.insert_stratumAbove)
        insertStratumBelowButton = QPushButton("Insert Stratum Below")
        self.connect(insertStratumBelowButton, SIGNAL("clicked()"), self.insert_stratumBelow)
        deleteStratumButton = QPushButton("Delete Stratum")
        self.connect(deleteStratumButton, SIGNAL("clicked()"), self.delete_stratum)
        
        #vBox2
        vBox2Grid = QGridLayout()
        vBox2Grid.addWidget(insertStratumAboveButton, 1,0, 1,1)
        vBox2Grid.addWidget(insertStratumBelowButton, 1,1, 1,1)
        vBox2Grid.addWidget(deleteStratumButton, 1,2, 1,1)
        vBox2Grid.addWidget(self.stratumTable, 2,0, 1,3)
        
        #Optimization
        optiFrameLayout = QGridLayout()
        
        self.opti1DChoiceBox = QComboBox()
        self.opti1DChoiceBox.addItems(['Thickness','Doping'])
        self.connect(self.opti1DChoiceBox, SIGNAL("currentIndexChanged(const QString &)"), self.input_opti1DChoice)
        optiFrameLayout.addWidget(self.opti1DChoiceBox, 1,0, 1,1)
        
        opti1DLayerBoxLabel = QLabel('<b><center>1<sup>st</sup> Dimension<br>Strata Number(s)</center></b>')
        opti1DLayerBoxLabel.setToolTip('Ex: 6,8')
        optiFrameLayout.addWidget(opti1DLayerBoxLabel, 0,1, 1,1)
        self.opti1DLayerBox = QLineEdit()
        self.connect(self.opti1DLayerBox, SIGNAL('editingFinished()'), self.input_opti1DLayer)
        self.opti1DLayerBox.setToolTip('Ex: 6,8')
        optiFrameLayout.addWidget(self.opti1DLayerBox, 1,1, 1,1)
        
        opti1DRangeBoxLabel = QLabel('<b><center>Optimization<br>Range</center></b>')
        opti1DRangeBoxLabel.setToolTip('Ex: 1:0.1:3')
        optiFrameLayout.addWidget(opti1DRangeBoxLabel, 0,2, 1,1)
        self.opti1DRangeBox = QLineEdit()
        self.opti1DRangeBox.setToolTip('Ex: 1:0.1:3')
        self.connect(self.opti1DRangeBox, SIGNAL('editingFinished()'), self.input_opti1DRange)
        optiFrameLayout.addWidget(self.opti1DRangeBox, 1,2, 1,1)
        
        self.opti1DRunButton = QPushButton('Optimize 1D')
        self.connect(self.opti1DRunButton, SIGNAL('clicked(bool)'), self.run_opti1D)
        optiFrameLayout.addWidget(self.opti1DRunButton, 1,3, 1,1)

        self.opti2DChoiceBox = QComboBox()
        self.opti2DChoiceBox.addItems(['Thickness','Doping'])
        self.connect(self.opti2DChoiceBox, SIGNAL("currentIndexChanged(const QString &)"), self.input_opti2DChoice)
        optiFrameLayout.addWidget(self.opti2DChoiceBox, 3,0, 1,1)
        
        opti2DLayerBoxLabel = QLabel('<b><center>2<sup>nd</sup> Dimension<br>Strata Number(s)</center></b>')
        opti2DLayerBoxLabel.setToolTip('Ex: 2 5 7')
        optiFrameLayout.addWidget(opti2DLayerBoxLabel, 2,1, 1,1)
        self.opti2DLayerBox = QLineEdit()
        self.opti2DLayerBox.setToolTip('Ex: 2 5 7')
        self.connect(self.opti2DLayerBox, SIGNAL('editingFinished()'), self.input_opti2DLayer)
        optiFrameLayout.addWidget(self.opti2DLayerBox, 3,1, 1,1)
        
        opti2DRangeBoxLabel = QLabel('<b><center>Optimization<br>Range</center></b>')
        opti2DRangeBoxLabel.setToolTip('Ex: 1:5')
        optiFrameLayout.addWidget(opti2DRangeBoxLabel, 2,2, 1,1)
        self.opti2DRangeBox = QLineEdit()
        self.opti2DRangeBox.setToolTip('Ex: 1:5')
        self.connect(self.opti2DRangeBox, SIGNAL('editingFinished()'), self.input_opti2DRange)
        optiFrameLayout.addWidget(self.opti2DRangeBox, 3,2, 1,1)
        
        self.opti2DRunButton = QPushButton('Optimize 2D')
        self.connect(self.opti2DRunButton, SIGNAL('clicked(bool)'), self.run_opti2D)
        optiFrameLayout.addWidget(self.opti2DRunButton, 3,3, 1,1)
        
        self.optiFrame = QGroupBox('Optimization')
        self.optiFrame.setMaximumWidth(450)
        self.optiFrame.setMinimumWidth(450)
        self.optiFrame.setLayout(optiFrameLayout)
        vBox2Grid.addWidget(self.optiFrame, 3,0, 1,3)
        
        vBox2 = QVBoxLayout()
        vBox2.addLayout(vBox2Grid)
        vBox2.addStretch()
        
        
        
        #vBox3
        self.plotModeButton = QPushButton("Plot Mode")
        self.connect(self.plotModeButton, SIGNAL("clicked()"), self.solve_mode)
        vBox3 = QVBoxLayout()
        vBox3.addWidget(self.plotModeButton)
        
        vBox3.addWidget(QLabel(u'<center><b><i>\u03B2<sub>eff</sub></i></b></center>'))
        self.betaEffBox = QLineEdit()
        self.betaEffBox.setMaximumWidth(150)
        self.betaEffBox.setEnabled(False)
        vBox3.addWidget(self.betaEffBox)
        
        self.modeCalculationsBox = QTextEdit('')
        self.modeCalculationsBox.setReadOnly(True)
        self.modeCalculationsBox.setSizePolicy(QSizePolicy(QSizePolicy.Fixed,QSizePolicy.Fixed))
        self.modeCalculationsBox.setMaximumHeight(175)
        self.modeCalculationsBox.setMaximumWidth(150)
        vBox3.addWidget(self.modeCalculationsBox)
        
        vBox3.addStretch()
        
        
        #vBox4

        #set up opticalCanvas for stratum / mode plot
        self.opticalCanvas = Qwt.QwtPlot(self)
        self.opticalCanvas.setCanvasBackground(Qt.white)
        self.opticalCanvas.canvas().setCursor(Qt.ArrowCursor)
        
        #optical optimization canvas
        self.optimization1DCanvas = Qwt.QwtPlot(self)
        self.optimization1DCanvas.setCanvasBackground(Qt.white)
        self.optimization1DCanvas.canvas().setCursor(Qt.ArrowCursor)
        self.optimization1DCanvas.setVisible(False)
        
        #2D optical optimization canvas
        #optimization2DFig = Figure((5.0, 4.0), dpi=dpi)
        self.optimization2DFig = Figure()
        self.optimization2DCanvas = FigureCanvas(self.optimization2DFig)
        #self.optimization2DAxes = self.optimization2DFig.add_subplot(111, projection='3d')
        margins = [0.05,0.05,0.95,0.95]
        self.optimization2DAxes = self.optimization2DFig.add_axes(margins, projection='3d')
        self.optimization2DAxes.autoscale(enable=True, axis='both', tight=True)
        #get the background color of the central widget
        #bgColor = self.mainTabWidget.palette().brush(QPalette.Window).color().name()
        bgColorRed = self.mainTabWidget.palette().brush(QPalette.Window).color().red()
        bgColorBlue = self.mainTabWidget.palette().brush(QPalette.Window).color().blue()
        bgColorGreen = self.mainTabWidget.palette().brush(QPalette.Window).color().green()
        self.bgColor = (bgColorRed/255.0, bgColorGreen/255.0, bgColorBlue/255.0, 1)
        self.optimization2DAxes.patch.set_color(self.bgColor)
        self.optimization2DFig.patch.set_color(self.bgColor)
        self.optimization2DCanvas.setVisible(False)


        
        vBox4 = QVBoxLayout()
        vBox4.addWidget(self.opticalCanvas)
        vBox4.addWidget(self.optimization1DCanvas)
        vBox4.addWidget(self.optimization2DCanvas)
        
        opticalLayout = QHBoxLayout()
        opticalLayout.addLayout(vBox1)
        opticalLayout.addLayout(vBox2)
        opticalLayout.addLayout(vBox3)
        opticalLayout.addLayout(vBox4)  
        
        
        
        opticalWidget = QWidget()
        opticalWidget.setLayout(opticalLayout)
        opticalWidget.setAutoFillBackground(True)
        opticalWidget.setBackgroundRole(QPalette.Window)        
        
        # ###############################
        #
        # Thermal Tab
        #
        # ###############################
        
#        vBox1 = QVBoxLayout()
        
#        thermalTable = QTableWidget()
#        thermalTable.setSelectionBehavior(QTableWidget.SelectRows)
#        thermalTable.setSelectionMode(QTableWidget.SingleSelection)
#        thermalTable.setMaximumWidth(380)
#        thermalTable.setMinimumWidth(380)
#        self.connect(thermalTable,SIGNAL("itemChanged(QTableWidgetItem*)"),self.stratumTable_itemChanged)
#        self.connect(thermalTable,SIGNAL("itemSelectionChanged()"),self.stratumTable_itemSelectionChanged)
#        vBox1.addWidget(thermalTable)
        
#        thermalWidget = QWidget()
#        thermalWidget.setLayout(vBox1)
#        thermalWidget.setAutoFillBackground(True)
#        thermalWidget.setBackgroundRole(QPalette.Window)
        
        
        self.mainTabWidget.addTab(self.quantumWidget, 'Quantum')
        self.mainTabWidget.addTab(opticalWidget, 'Optical')
#        self.mainTabWidget.addTab(thermalWidget, 'Thermal')
        self.connect(self.mainTabWidget, SIGNAL('currentChanged(int)'), self.change_main_tab)
        
        self.setCentralWidget(self.mainTabWidget)
        
        self.layerTable.selectRow(0)
        self.layerTable.setFocus()




#===============================================================================
# Optical Tab Input Controls
#===============================================================================

    def update_stratum_inputBoxes(self):
        self.wavelengthBox.setValue(self.strata.wavelength)
        self.operatingFieldBox.setValue(self.strata.operatingField)
        self.ACPeriodLengthBox.setValue(self.strata.Lp)
        self.ACPeriodsBox.setValue(self.strata.Np)
        self.OperatingVoltageBox.setText('{0:.1f} V'.format(self.strata.Np*self.strata.operatingField/self.strata.Lp))
        self.aCoreBox.setText('{0:.3f} cm^-1'.format(self.strata.aCore))
        self.nCoreBox.setText('{0.real:2.3f}+{0.imag:1.3e}j'.format(self.strata.nCore))
        self.transitionBroadeningBox.setValue(self.strata.transitionBroadening * 1000) #display in meV
        self.tauUpperBox.setValue(self.strata.tauUpper)
        self.tauLowerBox.setValue(self.strata.tauLower)
        self.tauUpperLowerBox.setValue(self.strata.tauUpperLower)
        self.opticalDipoleBox.setValue(self.strata.opticalDipole)
        self.FoMBox.setText(u'{0:4.0f} ps \u212B^2'.format(self.strata.FoM))
        self.waveguideFacetsBox.setCurrentIndex(self.waveguideFacetsList.index(self.strata.waveguideFacets))
        self.waveguideLengthBox.setValue(self.strata.waveguideLength)
        self.customFacetBox.setValue(self.strata.customFacet * 100) #display in percent
        
        if self.strata.waveguideFacets == 'as-cleaved + as-cleaved':
            self.strata.frontFacet = ThePhysics.reflectivity(self.strata.beta)
            self.strata.backFacet = ThePhysics.reflectivity(self.strata.beta)
        elif self.strata.waveguideFacets == 'as-cleaved + perfect HR':
            self.strata.frontFacet = ThePhysics.reflectivity(self.strata.beta)
            self.strata.backFacet = 1
        elif self.strata.waveguideFacets == 'as-cleaved + perfect AR':
            self.strata.frontFacet = 1e-9
            self.strata.backFacet = ThePhysics.reflectivity(self.strata.beta)
        elif self.strata.waveguideFacets == 'perfect AR + perfect HR':
            self.strata.frontFacet = 1e-9
            self.strata.backFacet = 1
        elif self.strata.waveguideFacets == 'custom coating + as-cleaved':
            self.strata.frontFacet = self.strata.customFacet
            self.strata.backFacet = ThePhysics.reflectivity(self.strata.beta)
        elif self.strata.waveguideFacets == 'custom coating + perfect HR':
            self.strata.frontFacet = self.strata.customFacet
            self.strata.backFacet = 1
        elif self.strata.waveguideFacets == 'custom coating + perfect AR':
            self.strata.frontFacet = 1e-9
            self.strata.backFacet = self.strata.customFacet
        
    def update_modeCalculations_box(self):
        
        self.strata.calculate_performance_parameters()
        
        reportString = u""
        
        #confinement factor
        reportString += u"\u0393: <b>{0:3.1f}%</b><br>".format(self.strata.confinementFactor*100)
        
        #waveguide loss
        reportString += u"<i>\u03B1<sub>wg</sub></i> : {0:3.1f} cm<sup>-1</sup><br>".format(self.strata.waveguideLoss)
        
        #mirror loss
        reportString += u"<i>\u03B1<sub>mirror</sub></i> : {0:3.1f} cm<sup>-1</sup><br>".format(self.strata.mirrorLoss)
        
        #gain
        reportString += u"gain: {0:3.3f} cm/A<br>".format(self.strata.gain)
        
        #Jth0
        reportString += u"<i>J<sub>th0</sub></i> : <b>{0:3.3f} kA/cm<sup>2</sup></b><br>".format(self.strata.Jth0)
        
        #Ith0
        reportString += u"<i>I<sub>th0</sub></i> : {0:5.1f} mA<br>".format(self.strata.Ith0*1000)
        
        #Voltage
        reportString += u"<i>V<sub>op</sub></i> : {0:3.1f} V<br>".format(self.strata.operatingVoltage)
        
        #Voltage Efficiency
        reportString += u"<i>\u03B7<sub>voltage</sub></i> : {0:3.1f}%<br>".format(self.strata.voltageEfficiency*100)
        
        #Extraction Efficiency
        reportString += u"<i>\u03B7<sub>extraction</sub></i> : {0:3.1f}%<br>".format(self.strata.extractionEfficiency*100)
        
        #Inversion Efficiency
        reportString += u"<i>\u03B7<sub>inversion</sub></i> : {0:3.1f}%<br>".format(self.strata.inversionEfficiency*100)

        #Modal Efficiency
        reportString += u"<i>\u03B7<sub>modal</sub></i> : {0:3.1f}%<br>".format(self.strata.modalEfficiency*100)
      
        self.modeCalculationsBox.setText(reportString)
        
    def edit_optical_parameters(self, toggleState):
        if toggleState == True:
            self.wavelengthBox.setReadOnly(False)
            self.wavelengthBox.setStyleSheet('color:black')
            self.operatingFieldBox.setReadOnly(False)
            self.operatingFieldBox.setStyleSheet('color:black')
            self.ACPeriodLengthBox.setReadOnly(False)
            self.ACPeriodLengthBox.setStyleSheet('color:black')
            self.aCoreBox.setReadOnly(False)
            self.aCoreBox.setStyleSheet('color:black')
            self.tauUpperBox.setReadOnly(False)
            self.tauUpperBox.setStyleSheet('color:black')
            self.tauUpperLowerBox.setReadOnly(False)
            self.tauUpperLowerBox.setStyleSheet('color:black')
            self.tauLowerBox.setReadOnly(False)
            self.tauLowerBox.setStyleSheet('color:black')
            self.opticalDipoleBox.setReadOnly(False)
            self.opticalDipoleBox.setStyleSheet('color:black')
            self.transitionBroadeningBox.setReadOnly(False)
            self.transitionBroadeningBox.setStyleSheet('color:black')
        else:
            self.wavelengthBox.setReadOnly(True)
            self.wavelengthBox.setStyleSheet('color:gray')
            self.operatingFieldBox.setReadOnly(True)
            self.operatingFieldBox.setStyleSheet('color:gray')
            self.ACPeriodLengthBox.setReadOnly(True)
            self.ACPeriodLengthBox.setStyleSheet('color:gray')
            self.aCoreBox.setReadOnly(True)
            self.aCoreBox.setStyleSheet('color:gray')
            self.tauUpperBox.setReadOnly(True)
            self.tauUpperBox.setStyleSheet('color:gray')
            self.tauLowerBox.setReadOnly(True)
            self.tauLowerBox.setStyleSheet('color:gray')
            self.tauUpperLowerBox.setReadOnly(True)
            self.tauUpperLowerBox.setStyleSheet('color:gray')
            self.opticalDipoleBox.setReadOnly(True)
            self.opticalDipoleBox.setStyleSheet('color:gray')
            self.transitionBroadeningBox.setReadOnly(True)
            self.transitionBroadeningBox.setStyleSheet('color:gray')
            
    def input_wavelength(self, value):
        self.strata.wavelength = value
        
        self.dirty = True
        self.update_windowTitle()       
        self.update_stratum_inputBoxes()
        self.stratumTable_refresh()
            
    def input_operatingField(self, value):
        self.strata.operatingField = value
        
        self.dirty = True
        self.update_windowTitle()       
        self.update_stratum_inputBoxes()
        self.stratumTable_refresh()
    
    def input_ACPeriodLength(self, value):
        self.strata.Lp = value
        
        self.dirty = True
        self.update_windowTitle()       
        self.update_stratum_inputBoxes()
        self.stratumTable_refresh()
    
    def input_ACPeriods(self, value):
        self.strata.Np = value
        
        self.dirty = True
        self.update_windowTitle()       
        self.update_stratum_inputBoxes()
        self.stratumTable_refresh()
            
    def input_aCore(self):
        initialText = unicode(self.aCoreBox.text())
        txt = initialText.split()[0]
        try:
            value = float(txt)
            self.strata.aCore = value
            self.strata.nCore = self.strata.get_nCore(self.qclayers)
            
            self.dirty = True
            self.update_windowTitle()       
            self.update_stratum_inputBoxes()
            self.stratumTable_refresh()
        except ValueError:
            self.aCore.setText(initialText)            
            
    def input_transitionBroadening(self, value):
        self.strata.transitionBroadening = value / 1000
        
        self.dirty = True
        self.update_windowTitle()       
        self.update_stratum_inputBoxes()
        self.stratumTable_refresh()
            
    def input_tauUpper(self, value):
        self.strata.tauUpper = value
        self.strata.FoM = self.strata.opticalDipole**2 * self.strata.tauUpper * (1- self.strata.tauLower/self.strata.tauUpperLower)
        
        self.dirty = True
        self.update_windowTitle()       
        self.update_stratum_inputBoxes()
        self.stratumTable_refresh()
        
    def input_tauLower(self, value):
        self.strata.tauLower = value
        self.strata.FoM = self.strata.opticalDipole**2 * self.strata.tauUpper * (1- self.strata.tauLower/self.strata.tauUpperLower)
        
        self.dirty = True
        self.update_windowTitle()       
        self.update_stratum_inputBoxes()
        self.stratumTable_refresh()
        
    def input_tauUpperLower(self, value):
        self.strata.tauUpperLower = value
        self.strata.FoM = self.strata.opticalDipole**2 * self.strata.tauUpper * (1- self.strata.tauLower/self.strata.tauUpperLower)
        
        self.dirty = True
        self.update_windowTitle()       
        self.update_stratum_inputBoxes()
        self.stratumTable_refresh()
        
    def input_opticalDipole(self, value):
        self.strata.opticalDipole = value
        self.strata.FoM = self.strata.opticalDipole**2 * self.strata.tauUpper * (1- self.strata.tauLower/self.strata.tauUpperLower)
        
        self.dirty = True
        self.update_windowTitle()       
        self.update_stratum_inputBoxes()
        self.stratumTable_refresh()
        
    def input_waveguideFacets(self, selection):
        self.strata.waveguideFacets = selection
        if selection == 'as-cleaved + as-cleaved':
            self.customFacetBoxLabel.setStyleSheet('color:gray')
            self.customFacetBox.setEnabled(False)
        elif selection == 'as-cleaved + perfect HR':
            self.customFacetBoxLabel.setStyleSheet('color:gray')
            self.customFacetBox.setEnabled(False)
        elif selection == 'as-cleaved + perfect AR':
            self.customFacetBoxLabel.setStyleSheet('color:gray')
            self.customFacetBox.setEnabled(False)
        elif selection == 'perfect AR + perfect HR':
            self.customFacetBoxLabel.setStyleSheet('color:gray')
            self.customFacetBox.setEnabled(False)
        elif selection == 'custom coating + as-cleaved':
            self.customFacetBoxLabel.setStyleSheet('color:black')
            self.customFacetBox.setEnabled(True)
        elif selection == 'custom coating + perfect HR':
            self.customFacetBoxLabel.setStyleSheet('color:black')
            self.customFacetBox.setEnabled(True)
        elif selection == 'custom coating + perfect AR':
            self.customFacetBoxLabel.setStyleSheet('color:black')
            self.customFacetBox.setEnabled(True)
        
        self.dirty = True
        self.update_windowTitle()       
        self.update_stratum_inputBoxes()
        self.stratumTable_refresh()
      
    def input_waveguideLength(self, value):
        self.strata.waveguideLength = value
        
        self.dirty = True
        self.update_windowTitle()       
        self.update_stratum_inputBoxes()
        self.stratumTable_refresh()
        
    def input_customFacet(self, value):
        self.strata.customFacet = value / 100.0
        
        self.dirty = True
        self.update_windowTitle()       
        self.update_stratum_inputBoxes()
        self.stratumTable_refresh()
        
    def input_opti1DChoice(self, selectionString):
        pass
    
    def input_opti1DLayer(self):
        pass
    
    def input_opti1DRange(self):
        pass
    
    def input_opti2DChoice(self, selectionString):
        pass
    
    def input_opti2DLayer(self):
        pass
    
    def input_opti2DRange(self):
        pass




#===============================================================================
# Optical Tab Strata Table Control
#===============================================================================

    def stratumTable_refresh(self):
        #calculate index for each layer
        self.strata.populate_rIndexes()
        
        #set properties for Active Core Layer
        for q in xrange(self.strata.stratumDopings.size):
            if self.strata.stratumMaterials[q] == 'Active Core':
                self.strata.stratumThicknesses[q] = self.strata.Np * self.strata.Lp * 1e-4
                self.strata.stratumDopings[q] = self.strata.nD
                self.strata.stratumRIndexes[q] = self.strata.nCore
        
        #update table
        self.stratumTable.clear()
        self.stratumTable.setColumnCount(6)
        self.stratumTable.setRowCount(self.strata.stratumDopings.size) #need to change
        self.stratumTable.setHorizontalHeaderLabels(['Material', 'Mole Frac', 'Thickness', 'Doping', 'Refractive Index', 'Loss'])
        vertLabels = []
        for n in xrange(self.strata.stratumDopings.size):
            vertLabels.append(str(n+1))
        self.stratumTable.setVerticalHeaderLabels(vertLabels)
        
        for q in xrange(self.strata.stratumDopings.size):
            #Stratum Material Setup
            materialWidget = QComboBox()
            materialWidget.addItems(self.stratumMaterialsList)
            materialWidget.setCurrentIndex(self.stratumMaterialsList.index(self.strata.stratumMaterials[q]))
            self.connect(materialWidget, SIGNAL("currentIndexChanged(int)"), partial(self.stratumTable_materialChanged, q))
            self.stratumTable.setCellWidget(q, 0, materialWidget)
            
            #Stratum Composition Setup
            composition = QTableWidgetItem()
            if self.strata.stratumMaterials[q] not in self.strata.needsCompositionList:
                composition.setFlags(Qt.NoItemFlags)
            else:
                composition.setData(0,'{0:3.3f}'.format(self.strata.stratumCompositions[q]))
                composition.setTextAlignment(Qt.AlignCenter)
            self.stratumTable.setItem(q, 1, composition)
            
            #Stratum Thickness Setup
            thickness = QTableWidgetItem(unicode(self.strata.stratumThicknesses[q]))
            thickness.setTextAlignment(Qt.AlignCenter)
            self.stratumTable.setItem(q, 2, thickness)
            if self.strata.stratumMaterials[q] == 'Active Core':
                thickness.setFlags(Qt.NoItemFlags)
            
            #Stratum Doping Setup
            doping = QTableWidgetItem()
            if self.strata.stratumMaterials[q] in self.strata.notDopableList:
                doping.setFlags(Qt.NoItemFlags)
            else:
                doping.setData(0,'{0:3.2f}'.format(self.strata.stratumDopings[q]))
                doping.setTextAlignment(Qt.AlignCenter)
                if self.strata.stratumMaterials[q] == 'Active Core':
                    doping.setFlags(Qt.NoItemFlags)
            self.stratumTable.setItem(q, 3, doping)
            
            #Stratum RIndex Setup
            rIndex = QTableWidgetItem('{0.real:2.3f}+{0.imag:1.3e}j'.format(self.strata.stratumRIndexes[q]))
            rIndex.setTextAlignment(Qt.AlignCenter)
            rIndex.setFlags(Qt.NoItemFlags)
            self.stratumTable.setItem(q, 4, rIndex)
            
            #Stratum Loss Setup
            loss = self.strata.stratumRIndexes[q].imag*4*pi/self.strata.wavelength/1e-4
            alpha = QTableWidgetItem('{0:3.2f}'.format(loss))
            alpha.setTextAlignment(Qt.AlignCenter)
            alpha.setFlags(Qt.NoItemFlags)
            self.stratumTable.setItem(q, 5, alpha)
            
        self.stratumTable.resizeColumnsToContents()
        
        self.update_opticalCanvas()

    def stratumTable_itemChanged(self, item):
        column = self.stratumTable.currentColumn()
        row = self.stratumTable.currentRow()
        if column == -1: #column == -1 on GUI initialization
            return
        elif column == 0:
            return
        elif column == 1:
            xFrac = float(item.text())
            if xFrac < 0 or xFrac > 1:
                QMessageBox.warning(self,'ErwinJr Error','Mole Fraction must be between 0 and 1')
            else:
                self.strata.stratumCompositions[row] = xFrac
        elif column == 2:
            self.strata.stratumThicknesses[row] = float(item.text())
        elif column == 3:
            self.strata.stratumDopings[row] = float(item.text())
        
        self.stratumTable_refresh()
        self.stratumTable.selectRow(row)
        
        self.dirty = True
        self.update_windowTitle()   
            
    def stratumTable_itemSelectionChanged(self):
        self.strata.stratumSelected = self.stratumTable.currentRow()
        if self.strata.stratumSelected >= 0 and self.strata.stratumSelected < self.qclayers.layerWidths.size:
            self.strata.populate_x()
            self.update_opticalCanvas()

    def stratumTable_materialChanged(self, row, selection):
        self.strata.stratumMaterials[row] = self.stratumMaterialsList[selection]
        
        self.stratumTable_refresh()
        self.stratumTable.selectRow(row)

        self.dirty = True
        self.update_windowTitle()  

    def insert_stratumAbove(self):
        row = self.stratumTable.currentRow()
        if row == -1:
            return
        
        if row == 0:
            self.strata.stratumMaterials.insert(row, self.strata.stratumMaterials[row])
            self.strata.stratumCompositions = hstack([self.strata.stratumCompositions[row], self.strata.stratumCompositions[row:,]])
            self.strata.stratumThicknesses = hstack([self.strata.stratumThicknesses[row], self.strata.stratumThicknesses[row:,]])
            self.strata.stratumDopings = hstack([self.strata.stratumDopings[row], self.strata.stratumDopings[row:,]])
            self.strata.stratumRIndexes = hstack([self.strata.stratumRIndexes[row], self.strata.stratumRIndexes[row:,]])

        else:
            self.strata.stratumMaterials.insert(row, self.strata.stratumMaterials[row])
            self.strata.stratumCompositions = hstack([self.strata.stratumCompositions[0:row], self.strata.stratumCompositions[row], self.strata.stratumCompositions[row:,]])
            self.strata.stratumThicknesses = hstack([self.strata.stratumThicknesses[0:row], self.strata.stratumThicknesses[row], self.strata.stratumThicknesses[row:,]])
            self.strata.stratumDopings = hstack([self.strata.stratumDopings[0:row], self.strata.stratumDopings[row], self.strata.stratumDopings[row:,]])
            self.strata.stratumRIndexes = hstack([self.strata.stratumRIndexes[0:row], self.strata.stratumRIndexes[row], self.strata.stratumRIndexes[row:,]])
            

        self.stratumTable_refresh()
        self.stratumTable.selectRow(row)
        self.stratumTable.setFocus()
        
        self.dirty = True
        self.update_windowTitle()
    
    def insert_stratumBelow(self):
        row = self.stratumTable.currentRow()
        if row == -1:
            return
        
        if row == self.strata.stratumDopings.size-1:
            self.strata.stratumMaterials.insert(row, self.strata.stratumMaterials[row])
            self.strata.stratumCompositions = hstack([self.strata.stratumCompositions[:], self.strata.stratumCompositions[row]])
            self.strata.stratumThicknesses = hstack([self.strata.stratumThicknesses[:], self.strata.stratumThicknesses[row]])
            self.strata.stratumDopings = hstack([self.strata.stratumDopings[:], self.strata.stratumDopings[row]])
            self.strata.stratumRIndexes = hstack([self.strata.stratumRIndexes[:], self.strata.stratumRIndexes[row]])

        else:
            self.strata.stratumMaterials.insert(row, self.strata.stratumMaterials[row])
            self.strata.stratumCompositions = hstack([self.strata.stratumCompositions[0:row+1], self.strata.stratumCompositions[row], self.strata.stratumCompositions[row+1:,]])
            self.strata.stratumThicknesses = hstack([self.strata.stratumThicknesses[0:row+1], self.strata.stratumThicknesses[row], self.strata.stratumThicknesses[row+1:,]])
            self.strata.stratumDopings = hstack([self.strata.stratumDopings[0:row+1], self.strata.stratumDopings[row], self.strata.stratumDopings[row+1:,]])
            self.strata.stratumRIndexes = hstack([self.strata.stratumRIndexes[0:row+1], self.strata.stratumRIndexes[row], self.strata.stratumRIndexes[row+1:,]])

        self.stratumTable_refresh()
        self.stratumTable.selectRow(row+1)
        self.stratumTable.setFocus()
        
        self.dirty = True
        self.update_windowTitle()
    
    def delete_stratum(self):
        #don't delete last stratum
        if self.strata.stratumDopings.size == 1:
            return
        
        row = self.stratumTable.currentRow()
        if row == -1:
            return
        
        self.strata.stratumMaterials.pop(row)
        self.strata.stratumThicknesses = hstack([self.strata.stratumCompositions[0:row], self.strata.stratumCompositions[row+1:,]])
        self.strata.stratumThicknesses = hstack([self.strata.stratumThicknesses[0:row], self.strata.stratumThicknesses[row+1:,]])
        self.strata.stratumDopings = hstack([self.strata.stratumDopings[0:row], self.strata.stratumDopings[row+1:,]])
        self.strata.stratumRIndexes = hstack([self.strata.stratumRIndexes[0:row], self.strata.stratumRIndexes[row+1:,]])

        #if current row was last row (now deleted)
        if row+1 > self.strata.stratumThicknesses.size:
            self.strata.stratumSelected -= 1
            row -= 1

        self.stratumTable.selectRow(row)
        self.stratumTable_refresh()
        self.stratumTable.selectRow(row)
        self.stratumTable.setFocus()
        
        self.dirty = True
        self.update_windowTitle()        




#===============================================================================
# Optical Tab Plotting and Plot Control
#===============================================================================
        
    def update_opticalCanvas(self):
        self.strata.populate_x()
        
        self.opticalCanvas.clear()
        
        self.curvenR = Qwt.QwtPlotCurve()
        self.curvenR.setData(self.strata.xPoints,self.strata.xn.real)
        self.curvenR.setPen(QPen(Qt.black, 1.5))
        if settings.antialiased:
            self.curvenR.setRenderHint(Qwt.QwtPlotItem.RenderAntialiased)
        self.curvenR.attach(self.opticalCanvas)
        self.opticalCanvas.setAxisTitle(Qwt.QwtPlot.yLeft, 'Refractive Index')
        
        if self.strata.stratumSelected >= 0 and self.strata.stratumSelected < self.strata.stratumThicknesses.size:
            mask = ~isnan(self.strata.xStratumSelected)
            self.stratumSelection = SupportClasses.MaskedCurve(self.strata.xPoints,self.strata.xStratumSelected,mask)
            self.stratumSelection.setPen(QPen(Qt.blue, 2))
            if settings.antialiased:
                self.stratumSelection.setRenderHint(Qwt.QwtPlotItem.RenderAntialiased)
            self.stratumSelection.attach(self.opticalCanvas)
        
        #plot Intensity
        if hasattr(self.strata,'xI'):
            self.curvexI = Qwt.QwtPlotCurve()
            self.curvexI.setData(self.strata.xPoints, self.strata.xI*self.strata.stratumRIndexes.real.max())
            self.curvexI.setPen(QPen(Qt.red, 1.5))
            if settings.antialiased:
                self.curvexI.setRenderHint(Qwt.QwtPlotItem.RenderAntialiased)
            self.curvexI.attach(self.opticalCanvas)
            self.opticalCanvas.setAxisTitle(Qwt.QwtPlot.yLeft, 'Refractive Index, Mode Intensity')
            
        self.opticalCanvas.setAxisTitle(Qwt.QwtPlot.xBottom, u'Position (\u03BCm)')
        self.opticalCanvas.replot()

    def solve_mode(self):
        betaInit = self.betaEffBox.text()
        if betaInit == '':
            betaInit = None
        else:
            betaInit = complex(str(betaInit))
        self.strata.beta = self.strata.beta_find(betaInit)
        self.betaEffBox.setText('{0.real:2.3f}+{0.imag:1.3e}j'.format(self.strata.beta))
        self.strata.mode_plot()
        self.update_modeCalculations_box()
        self.update_opticalCanvas()
        
    def run_opti1D(self):
        #get initial parameters
        try:
            optiType1D  = self.opti1DChoiceBox.currentText()
            strata1D    = array(SupportClasses.matlab_range(self.opti1DLayerBox.text()), dtype=int)
            strata1D   -= 1 #indexing starts at 0
            optiRange1D = array(SupportClasses.matlab_range(self.opti1DRangeBox.text()))
        except ValueError:
            QMessageBox.warning(self,"ErwinJr Error", "Invalid entry.")
            return
        
        #set up GUI
        self.optiFrame.setEnabled(False)
        self.plotModeButton.setEnabled(False)
        
        Jth0Array = zeros(optiRange1D.size)*NaN
        ylabel = '<i>J<sub>th0</sub></i>'
        
        stratumThicknessesInitial = self.strata.stratumThicknesses.copy()
        stratumDopingsInitial = self.strata.stratumDopings.copy()
        
        for q, rangeValue in enumerate(optiRange1D):
            if optiType1D == 'Thickness':
                self.strata.stratumThicknesses[strata1D] = rangeValue
                xlabel = u'Thickness (\u03BCm)'
            elif optiType1D == 'Doping':
                self.strata.stratumDopings[strata1D] = rangeValue
                xlabel = 'Doping (x10<sup>17</sup> cm<sup>-3</sup>)'
            elif optiType1D == 'Active Core Periods':
                pass
            elif optiTyp1D == 'deltaE':
                pass
            elif optiTyp1D == 'E3c':
                pass
            elif optiTyp1D == 'Custom Facet':
                pass
            elif optiTyp1D == 'Waveguide Length':
                pass
            elif optiTyp1D == 'Ridge Width':
                pass
            elif optiType1D == 'Tsink':
                pass
            self.stratumTable_refresh()
            self.update_stratum_inputBoxes()
            self.solve_mode()
            Jth0Array[q] = self.strata.Jth0
            self.plot_on_optimization1DCanvas(optiRange1D, xlabel, Jth0Array, ylabel)
        
        #reset initial values
        self.strata.stratumThicknesses = stratumThicknessesInitial
        self.strata.stratumDopings = stratumDopingsInitial
        self.update_stratum_inputBoxes()
        self.stratumTable_refresh()
        self.solve_mode()
            
        #reset GUI
        self.optiFrame.setEnabled(True)
        self.plotModeButton.setEnabled(True)
            
    def plot_on_optimization1DCanvas(self, xVals, xlabel, yVals, ylabel):
        self.optimization2DCanvas.setVisible(False)
        self.optimization1DCanvas.setVisible(True)
        
        self.optimization1DCanvas.clear()
        
        mask = ~isnan(yVals)
        optiCurve = SupportClasses.MaskedCurve(xVals, yVals, mask)
        optiCurve.setPen(QPen(Qt.blue, 1.5))
        if settings.antialiased:
            optiCurve.setRenderHint(Qwt.QwtPlotItem.RenderAntialiased)
        optiCurve.attach(self.optimization1DCanvas)

        self.optimization1DCanvas.setAxisTitle(Qwt.QwtPlot.xBottom, xlabel)
        self.optimization1DCanvas.setAxisTitle(Qwt.QwtPlot.yLeft, ylabel)
        self.optimization1DCanvas.replot()
            
    def run_opti2D(self):
        #get initial parameters
        try:
            optiType1D  = self.opti1DChoiceBox.currentText()
            strata1D    = array(SupportClasses.matlab_range(self.opti1DLayerBox.text()), dtype=int)
            strata1D   -= 1 #indexing starts at 0
            optiRange1D = array(SupportClasses.matlab_range(self.opti1DRangeBox.text()))
            optiType2D  = self.opti2DChoiceBox.currentText()
            strata2D    = array(SupportClasses.matlab_range(self.opti2DLayerBox.text()), dtype=int)
            strata2D   -= 1 #indexing starts at 0
            optiRange2D = array(SupportClasses.matlab_range(self.opti2DRangeBox.text()))
        except ValueError:
            QMessageBox.warning(self,"ErwinJr Error", "Invalid entry.")
            return
        
        Jth0Array = NaN * zeros((optiRange1D.size, optiRange2D.size))
        zlabel = '$J_{th0}$'
        
        stratumThicknessesInitial = self.strata.stratumThicknesses.copy()
        stratumDopingsInitial = self.strata.stratumDopings.copy()
        
        for qq, rangeValue2D in enumerate(optiRange2D):
            if optiType2D == 'Thickness':
                self.strata.stratumThicknesses[strata2D] = rangeValue2D
                xlabel = u'Thickness ($\mu m$)'
            elif optiType2D == 'Doping':
                self.strata.stratumDopings[strata2D] = rangeValue2D
                xlabel = 'Doping ($x10^{17} cm^{-3}$)'
            for q, rangeValue1D in enumerate(optiRange1D):
                if optiType1D == 'Thickness':
                    self.strata.stratumThicknesses[strata1D] = rangeValue1D
                    ylabel = u'Thickness ($\mu m$)'
                elif optiType1D == 'Doping':
                    self.strata.stratumDopings[strata1D] = rangeValue1D
                    ylabel = 'Doping ($x10^{17} cm^{-3}$)'
    
                self.stratumTable_refresh()
                self.update_stratum_inputBoxes()
                self.solve_mode()
                Jth0Array[q,qq] = self.strata.Jth0
                self.plot_on_optimization2DCanvas(optiRange1D, xlabel, optiRange2D, ylabel, Jth0Array, zlabel)
                QCoreApplication.processEvents()
                
        #reset initial values
        self.strata.stratumThicknesses = stratumThicknessesInitial
        self.strata.stratumDopings = stratumDopingsInitial
        self.update_stratum_inputBoxes()
        self.stratumTable_refresh()
        self.solve_mode()
        
        #reset GUI
        self.optiFrame.setEnabled(True)
        self.plotModeButton.setEnabled(True)
    
    def plot_on_optimization2DCanvas(self, xVals, xlabel, yVals, ylabel, zVals, zlabel):
        self.optimization1DCanvas.setVisible(False)
        self.optimization2DCanvas.setVisible(True)
        
        X,Y = meshgrid(yVals, xVals)
        Z = zVals

        self.optimization2DAxes.cla()
        self.optimization2DAxes.patch.set_color(self.bgColor)
        self.optimization2DAxes.mouse_init()
        
        normd = matplotlib.colors.Normalize(nanmin(nanmin(Z)), nanmax(nanmax(Z)))
        self.optimization2DAxes.plot_surface(X, Y, Z, cstride=1, rstride=1, norm=normd, 
                                             cmap=matplotlib.cm.Blues_r, linewidth=0, 
                                             antialiased=False, shade=False)
        self.optimization2DAxes.set_zlim(0.95*nanmin(nanmin(Z)), 1.05*nanmax(nanmax(Z)))
        self.optimization2DAxes.set_xlabel(xlabel)
        self.optimization2DAxes.set_ylabel(ylabel)
        self.optimization2DAxes.set_zlabel(zlabel)
        self.optimization2DCanvas.draw()




#===============================================================================
# Quantum Tab Input Controls
#===============================================================================

    def update_inputBoxes(self):
        self.inputMoleFrac1Box.setValue(self.qclayers.moleFrac1)
        self.inputMoleFrac2Box.setValue(self.qclayers.moleFrac2)
        self.inputMoleFrac3Box.setValue(self.qclayers.moleFrac3)
        self.inputMoleFrac4Box.setValue(self.qclayers.moleFrac4)
        self.inputMoleFrac5Box.setValue(self.qclayers.moleFrac5)
        self.inputMoleFrac6Box.setValue(self.qclayers.moleFrac6)
        self.inputMoleFrac7Box.setValue(self.qclayers.moleFrac7)
        self.inputMoleFrac8Box.setValue(self.qclayers.moleFrac8)
        
        self.qclayers.update_alloys()
        self.qclayers.update_strain()
        self.qclayers.populate_x()
        self.update_quantumCanvas()
        
        self.offset1Label.setText("%6.0f meV" % ((self.qclayers.EcG[1]-self.qclayers.EcG[0])*1000))
        self.offset2Label.setText("%6.0f meV" % ((self.qclayers.EcG[3]-self.qclayers.EcG[2])*1000))
        self.offset3Label.setText("%6.0f meV" % ((self.qclayers.EcG[5]-self.qclayers.EcG[4])*1000))
        self.offset4Label.setText("%6.0f meV" % ((self.qclayers.EcG[7]-self.qclayers.EcG[6])*1000))
        
        strainString = "<center>Net Strain: <b>%6.3f%%</b></center>" % self.qclayers.mismatch
        self.strainDescription.setText(strainString)
        
        self.inputVertResBox.setValue(self.qclayers.vertRes)
        self.inputEFieldBox.setValue(self.qclayers.EField)
        self.inputRepeatsBox.setValue(self.qclayers.repeats)
        self.inputHorzResBox.setCurrentIndex(self.inputHorzResBox.findText(
                                              QString(unicode(self.qclayers.xres))))
                                              
        self.DescriptionBox.setText(self.qclayers.description)
        
        self.update_Lp_box()
        
        self.dirty = True
        self.update_windowTitle()
        
    def input_substrate(self, substrateType):
        if substrateType == 'InP':
            self.qclayers.substrate = 'InP'
            self.materialList = ['InGaAs/AlInAs #1', 'InGaAs/AlInAs #2', 'InGaAs/AlInAs #3', 'InGaAs/AlInAs #4']
            self.mtrl_header1.setText('<center><b>In<sub>x</sub>Ga<sub>1-x</sub>As</b></center')
            self.mtrl_header2.setText('<center><b>Al<sub>1-x</sub>In<sub>x</sub>As</b></center')
            
            self.quantumCanvas.clear()
            #self.layerTable_refresh()
            self.update_Lp_limits()
            self.update_inputBoxes()
            self.layerTable_refresh()
            self.qclayers.populate_x()
            
        elif substrateType == 'GaAs':
            self.qclayers.substrate = 'GaAs'
            self.materialList = ['AlGaAs/AlGaAs #1', 'AlGaAs/AlGaAs #2', 'AlGaAs/AlGaAs #3', 'AlGaAs/AlGaAs #4']
            self.mtrl_header1.setText('<center><b>Al<sub>x</sub>Ga<sub>1-x</sub>As</b></center')
            self.mtrl_header2.setText('<center><b>Al<sub>1</sub>Ga<sub>1-x</sub>As</b></center')
            
            self.quantumCanvas.clear()
            #self.layerTable_refresh()
            self.update_Lp_limits()
            self.update_inputBoxes()
            self.layerTable_refresh()
            self.qclayers.populate_x()
            
        elif substrateType == 'GaSb':
            self.qclayers.substrate = 'GaSb'
            self.materialList = ['InAsSb/AlGaSb #1', 'InAsSb/AlGaSb #2', 'InAsSb/AlGaSb #3', 'InAsSb/AlGaSb #4']
            self.mtrl_header1.setText('<center><b>InAs<sub>y</sub>Sb<sub>1-y</sub></b></center')
            self.mtrl_header2.setText('<center><b>Al<sub>x</sub>Ga<sub>1-x</sub>Sb</b></center')
            
            self.quantumCanvas.clear()
            #self.layerTable_refresh()
            self.update_Lp_limits()
            self.update_inputBoxes()
            self.layerTable_refresh()
            self.qclayers.populate_x()
            
        elif substrateType == 'GaN':
            QMessageBox.information(self, 'ErwinJr Error', 'III-Nitride substrates have not yet been implemented.')
            
        else:
            raise TypeError('substrate selection not allowed')
        
    def input_EField(self):
        self.qclayers.EField = float(self.inputEFieldBox.value())
        
        self.qclayers.populate_x()
        self.update_quantumCanvas()
        
        self.dirty = True
        self.update_windowTitle()
        
    def input_horzRes(self):
        horzRes = unicode(self.inputHorzResBox.currentText())
        horzRes = float(horzRes)
        self.qclayers.xres = horzRes

        self.qclayers.populate_x()
        self.update_quantumCanvas()
        
        self.dirty = True
        self.update_windowTitle()

    def input_vertRes(self):
        self.qclayers.vertRes = float(self.inputVertResBox.value())

        self.dirty = True
        self.update_windowTitle()

    def input_repeats(self):
        self.qclayers.repeats = int(self.inputRepeatsBox.value())
            
        self.qclayers.populate_x()
        self.update_quantumCanvas()
        
        self.dirty = True
        self.update_windowTitle()
        
    def input_basis(self):
        self.qclayers.basisARInjector = self.inputARInjectorCheck.isChecked()
        self.qclayers.basisInjectorAR = self.inputInjectorARCheck.isChecked()
        self.dirty = True
        self.update_windowTitle()

    def update_Lp_limits(self):
        self.LpFirstSpinbox.setRange(1,self.qclayers.layerWidths.size-1)
        self.LpFirstSpinbox.setValue(1)
        self.LpLastSpinbox.setRange(1,self.qclayers.layerWidths.size-1)
        self.LpLastSpinbox.setValue(self.qclayers.layerWidths.size-1)
        
    def update_Lp_box(self):
        #self.LpFirstSpinbox.setRange(1,self.qclayers.layerWidths.size-1)
        #self.LpFirstSpinbox.setValue(1)
        #self.LpLastSpinbox.setRange(1,self.qclayers.layerWidths.size-1)
        #self.LpLastSpinbox.setValue(self.qclayers.layerWidths.size-1)
        LpFirst = self.LpFirstSpinbox.value()
        LpLast = self.LpLastSpinbox.value()+1 #+1 because range is not inclusive of last value
        Lp_string  = u"Lp: %g \u212B<br>" % (sum(self.qclayers.layerWidths[LpFirst:LpLast]))
        Lp_string += u"wells: %6.1f%%<br>" % (sum((1-self.qclayers.layerBarriers[LpFirst:LpLast])*self.qclayers.layerWidths[LpFirst:LpLast]) / sum(self.qclayers.layerWidths[LpFirst:LpLast]) * 100)
        Lp_string += u"n<sub>D</sub>: %6.3f\u00D710<sup>17</sup><br>" % (sum(self.qclayers.layerDopings[LpFirst:LpLast]*self.qclayers.layerWidths[LpFirst:LpLast]) / sum(self.qclayers.layerWidths[LpFirst:LpLast]))
        Lp_string += u"n<sub>s</sub>: %6.3f\u00D710<sup>11</sup>" % (sum(self.qclayers.layerDopings[LpFirst:LpLast]*self.qclayers.layerWidths[LpFirst:LpLast])*1e-2)
        self.LpStringBox.setText(Lp_string)

    def input_description(self):
        self.qclayers.description = self.DescriptionBox.toPlainText()
        self.dirty = True
        self.update_windowTitle()

    def input_moleFrac(self, boxID):
        if boxID == 1:
            self.qclayers.moleFrac1 = float(self.inputMoleFrac1Box.value())
        elif boxID == 2:
            self.qclayers.moleFrac2 = float(self.inputMoleFrac2Box.value())
        elif boxID == 3:
            self.qclayers.moleFrac3 = float(self.inputMoleFrac3Box.value())
        elif boxID == 4:
            self.qclayers.moleFrac4 = float(self.inputMoleFrac4Box.value())
        elif boxID == 5:
            self.qclayers.moleFrac5 = float(self.inputMoleFrac5Box.value())
        elif boxID == 6:
            self.qclayers.moleFrac6 = float(self.inputMoleFrac6Box.value())
        elif boxID == 7:
            self.qclayers.moleFrac7 = float(self.inputMoleFrac7Box.value())
        elif boxID == 8:
            self.qclayers.moleFrac8 = float(self.inputMoleFrac8Box.value())
            
        self.dirty = True
        self.update_windowTitle()
        
        self.update_inputBoxes()
        self.qclayers.update_alloys()
        self.qclayers.update_strain()
        self.qclayers.populate_x()
        self.update_quantumCanvas()



        
#===============================================================================
# Quantum Tab Layer Table Control
#===============================================================================
        
    def layerTable_refresh(self):
        self.layerTable.clear()
        self.layerTable.setColumnCount(6)
        self.layerTable.setRowCount(self.qclayers.layerWidths.size+1)
        self.layerTable.setHorizontalHeaderLabels(['Width', 'ML', 'Brr', 'AR', 'Doping', 'Material'])
        vertLabels = []
        for n in xrange(self.qclayers.layerWidths.size+1):
            vertLabels.append(str(n))
        self.layerTable.setVerticalHeaderLabels(vertLabels)        
        
        #color for barrier layers
        gray = QColor(230,230,240)
        gray2 = QColor(230,230,230)
        
        for q, layerWidth in enumerate(self.qclayers.layerWidths):
            #Width Setup
            width = QTableWidgetItem("%5.1f" % layerWidth)
            width.setTextAlignment(Qt.AlignCenter)
            if bool(self.qclayers.layerBarriers[q]):
                width.setBackgroundColor(gray)
            self.layerTable.setItem(q, 0, width)
            if q == 0:
                width.setFlags(Qt.NoItemFlags)
                width.setBackgroundColor(gray2)
            
            #ML Setup
            numML = layerWidth/self.qclayers.MLThickness[q]
            item = QTableWidgetItem("%5.1f" % numML)
            item.setTextAlignment(Qt.AlignCenter)
            if bool(self.qclayers.layerBarriers[q]):
                item.setBackgroundColor(gray)
            self.layerTable.setItem(q, 1, item)
            if q == 0:
                item.setFlags(Qt.NoItemFlags)
                item.setBackgroundColor(gray2)
                
            #Barrier Layer Setup
            item = QTableWidgetItem()
            item.setCheckState(int(self.qclayers.layerBarriers[q])*2)
            if bool(self.qclayers.layerBarriers[q]):
                item.setBackgroundColor(gray)
            self.layerTable.setItem(q, 2, item)
            if q == 0:
                item.setFlags(Qt.NoItemFlags)
                item.setBackgroundColor(gray2)

            #Active Region Layer Setup
            item = QTableWidgetItem()
            item.setCheckState(int(self.qclayers.layerARs[q])*2)
            if bool(self.qclayers.layerBarriers[q]):
                item.setBackgroundColor(gray)
            self.layerTable.setItem(q, 3, item)
            if q == 0:
                item.setFlags(Qt.NoItemFlags)
                item.setBackgroundColor(gray2)
            
            #Layer Doping Setup
            doping = QTableWidgetItem(unicode(self.qclayers.layerDopings[q]))
            doping.setTextAlignment(Qt.AlignCenter)
            if bool(self.qclayers.layerBarriers[q]):
                doping.setBackgroundColor(gray)
            self.layerTable.setItem(q, 4, doping)
            if q == 0:
                doping.setFlags(Qt.NoItemFlags)
                doping.setBackgroundColor(gray2)
                
            #Material Setup
            if q == 0:
                item = QTableWidgetItem(unicode(self.materialList[int(self.qclayers.layerMaterials[q]-1)]))
                item.setBackgroundColor(gray2)
                item.setFlags(Qt.NoItemFlags)
                self.layerTable.setItem(q, 5, item)
            else:
                materialWidget = QComboBox()
                materialWidget.addItems(self.materialList)
                materialWidget.setCurrentIndex(self.qclayers.layerMaterials[q]-1)
                self.connect(materialWidget, SIGNAL("currentIndexChanged(int)"), partial(self.layerTable_materialChanged, q))
                self.layerTable.setCellWidget(q, 5, materialWidget)
        
        self.layerTable.resizeColumnsToContents()

    def layerTable_itemChanged(self, item):
        column = self.layerTable.currentColumn()
        row = self.layerTable.currentRow()
        if column == -1: #column == -1 on GUI initialization
            return
        elif column == 0: #column == 0 for item change in Widths column
            if mod(float(item.text()), self.qclayers.xres) != 0 and self.qclayers.xres != 0.1:
                QMessageBox.warning(self,"ErwinJr - Warning",
                             "You entered a width that is not compatible with the minimum horizontal resolution.")
                return
            if row == self.qclayers.layerWidths.size: #add row at end of list
                self.qclayers.layerWidths = hstack([self.qclayers.layerWidths, float(item.text())])
                if self.qclayers.layerBarriers[-1] == 1:
                    self.qclayers.layerBarriers = hstack([self.qclayers.layerBarriers, 0])
                else:
                    self.qclayers.layerBarriers = hstack([self.qclayers.layerBarriers, 1])
                self.qclayers.layerARs = hstack([self.qclayers.layerARs, self.qclayers.layerARs[-1]])
                self.qclayers.layerMaterials = hstack([self.qclayers.layerMaterials, self.qclayers.layerMaterials[-1]])
                self.qclayers.layerDopings = hstack([self.qclayers.layerDopings, self.qclayers.layerDopings[-1]])
                self.qclayers.layerDividers = hstack([self.qclayers.layerDividers, self.qclayers.layerDividers[-1]])
                row += 1 #used so that last (blank) row is again selected
                
                #make first item the same as last item
                self.qclayers.layerWidths[0] = self.qclayers.layerWidths[-1]
                self.qclayers.layerBarriers[0] = self.qclayers.layerBarriers[-1]
                self.qclayers.layerARs[0] = self.qclayers.layerARs[-1]
                self.qclayers.layerMaterials[0] = self.qclayers.layerMaterials[-1]
                self.qclayers.layerDopings[0] = self.qclayers.layerDopings[-1]
                self.qclayers.layerDividers[0] = self.qclayers.layerDividers[-1]
                
                self.update_Lp_limits()
                
            elif row == self.qclayers.layerWidths.size-1:
                self.qclayers.layerWidths[row] = float(item.text())
                
                #make first item the same as last item
                self.qclayers.layerWidths[0] = self.qclayers.layerWidths[-1]
                self.qclayers.layerBarriers[0] = self.qclayers.layerBarriers[-1]
                self.qclayers.layerARs[0] = self.qclayers.layerARs[-1]
                self.qclayers.layerMaterials[0] = self.qclayers.layerMaterials[-1]
                self.qclayers.layerDopings[0] = self.qclayers.layerDopings[-1]
                self.qclayers.layerDividers[0] = self.qclayers.layerDividers[-1]  
            else: #change Width of selected row in-place
                self.qclayers.layerWidths[row] = float(item.text())
        elif column == 1: #column == 1 for ML
            if self.qclayers.xres != 0.1:
                QMessageBox.warning(self,"ErwinJr - Warning",
                             u"Horizontal Resolution of 0.1 \u212B required when setting monolayer thicknesses.")
                return
            if row == self.qclayers.layerWidths.size: #add row at end of list
                pass
            elif row == self.qclayers.layerWidths.size-1:
                self.qclayers.layerWidths[row] = self.qclayers.MLThickness[row] * float(item.text())
                
                #make first item the same as last item
                self.qclayers.layerWidths[0] = self.qclayers.layerWidths[-1]
                self.qclayers.layerBarriers[0] = self.qclayers.layerBarriers[-1]
                self.qclayers.layerARs[0] = self.qclayers.layerARs[-1]
                self.qclayers.layerMaterials[0] = self.qclayers.layerMaterials[-1]
                self.qclayers.layerDopings[0] = self.qclayers.layerDopings[-1]
                self.qclayers.layerDividers[0] = self.qclayers.layerDividers[-1]
                
                self.update_Lp_limits()
                
            else: #change Width of selected row in-place
                self.qclayers.layerWidths[row] = self.qclayers.MLThickness[row] * float(item.text())
        elif column == 2: #column == 2 for item change in Barrier column
            if row == self.qclayers.layerWidths.size: #don't do anything if row is last row
                return
            self.qclayers.layerBarriers[row] = int(item.checkState())//2
        elif column == 3: #column == 3 for item change in AR column
            if row == self.qclayers.layerWidths.size: #don't do anything if row is last row
                return
            self.qclayers.layerARs[row] = int(item.checkState())//2
        elif column == 4: #column == 4 for item change in Doping column
            if row == self.qclayers.layerWidths.size: #don't do anything if row is last row
                return
            self.qclayers.layerDopings[row] = float(item.text())
        elif column == 5: #column == 5 for item change in Materials column
            #self.qclayers.layerWidths[row] = int(item.text[row])
           pass
        else:
            pass
        self.layerTable_refresh()
        #self.qclayers.populate_x()
        #self.layerTable.selectRow(row)
        self.layerTable.setCurrentCell(row,column)
        self.layerTable.setFocus()
        
        self.update_Lp_box()

        self.dirty = True
        self.update_windowTitle()        
        
    def layerTable_materialChanged(self, row, selection):
        self.qclayers.layerMaterials[row] = selection+1
        #self.layerTable_refresh()
        self.qclayers.populate_x()
        self.layerTable.selectRow(row)

        self.dirty = True
        self.update_windowTitle()        

    def layerTable_itemSelectionChanged(self):
        #This is the primary call to update_quantumCanvas
        self.qclayers.layerSelected = self.layerTable.currentRow()
        if self.qclayers.layerSelected >= 0 and self.qclayers.layerSelected < self.qclayers.layerWidths.size:
            self.qclayers.populate_x()
            self.update_quantumCanvas()

    def delete_layer(self):
        #don't delete last layer
        if self.qclayers.layerWidths.size == 1:
            return
        
        row = self.layerTable.currentRow()
        if row == -1:
            return
        
        self.qclayers.layerWidths = hstack([self.qclayers.layerWidths[0:row], self.qclayers.layerWidths[row+1:,]])
        self.qclayers.layerBarriers = hstack([self.qclayers.layerBarriers[0:row], self.qclayers.layerBarriers[row+1:,]])
        self.qclayers.layerARs = hstack([self.qclayers.layerARs[0:row], self.qclayers.layerARs[row+1:,]])
        self.qclayers.layerMaterials = hstack([self.qclayers.layerMaterials[0:row], self.qclayers.layerMaterials[row+1:,]])
        self.qclayers.layerDopings = hstack([self.qclayers.layerDopings[0:row], self.qclayers.layerDopings[row+1:,]])
        self.qclayers.layerDividers = hstack([self.qclayers.layerDividers[0:row], self.qclayers.layerDividers[row+1:,]])
        
        if row == self.qclayers.layerWidths.size: #if row == last_row
            #make first item the same as last item
            self.qclayers.layerWidths[0] = self.qclayers.layerWidths[-1]
            self.qclayers.layerBarriers[0] = self.qclayers.layerBarriers[-1]
            self.qclayers.layerARs[0] = self.qclayers.layerARs[-1]
            self.qclayers.layerMaterials[0] = self.qclayers.layerMaterials[-1]
            self.qclayers.layerDopings[0] = self.qclayers.layerDopings[-1]
            self.qclayers.layerDividers[0] = self.qclayers.layerDiviers[-1]

        self.update_Lp_limits()
        self.update_Lp_box()

        self.qclayers.update_strain()
        self.layerTable_refresh()
        self.qclayers.populate_x()
        self.layerTable.selectRow(row)
        self.layerTable.setFocus()
        
        self.dirty = True
        self.update_windowTitle()
    
    def insert_layerAbove(self):
        row = self.layerTable.currentRow()
        if row == -1:
            return
        
        if row == 0:
            self.qclayers.layerWidths = hstack([0, self.qclayers.layerWidths[row:,]])
            if self.qclayers.layerBarriers[row] == 1:
                self.qclayers.layerBarriers = hstack([0, self.qclayers.layerBarriers[row:,]])
            else:
                self.qclayers.layerBarriers = hstack([1, self.qclayers.layerBarriers[row:,]])
            self.qclayers.layerARs = hstack([self.qclayers.layerARs[row], self.qclayers.layerARs[row:,]])
            self.qclayers.layerMaterials = hstack([self.qclayers.layerMaterials[row], self.qclayers.layerMaterials[row:,]])
            self.qclayers.layerDopings = hstack([self.qclayers.layerDopings[row], self.qclayers.layerDopings[row:,]])
            self.qclayers.layerDividers = hstack([self.qclayers.layerDividers[row], self.qclayers.layerDividers[row:,]])

        else:
            self.qclayers.layerWidths = hstack([self.qclayers.layerWidths[0:row], 0, self.qclayers.layerWidths[row:,]])
            if self.qclayers.layerBarriers[row] == 1:
                self.qclayers.layerBarriers = hstack([self.qclayers.layerBarriers[0:row], 0, self.qclayers.layerBarriers[row:,]])
            else:
                self.qclayers.layerBarriers = hstack([self.qclayers.layerBarriers[0:row], 1, self.qclayers.layerBarriers[row:,]])
            self.qclayers.layerARs = hstack([self.qclayers.layerARs[0:row], self.qclayers.layerARs[row], self.qclayers.layerARs[row:,]])
            self.qclayers.layerMaterials = hstack([self.qclayers.layerMaterials[0:row], self.qclayers.layerMaterials[row], self.qclayers.layerMaterials[row:,]])
            self.qclayers.layerDopings = hstack([self.qclayers.layerDopings[0:row], self.qclayers.layerDopings[row], self.qclayers.layerDopings[row:,]])
            self.qclayers.layerDividers = hstack([self.qclayers.layerDividers[0:row], self.qclayers.layerDividers[row], self.qclayers.layerDividers[row:,]])

        self.update_Lp_limits()
        self.update_Lp_box()
        
        self.layerTable_refresh()
        self.qclayers.populate_x()
        self.layerTable.selectRow(row)
        self.layerTable.setFocus()
        
        self.dirty = True
        self.update_windowTitle()            
            



#===============================================================================
# Quantum Tab Buttons
#===============================================================================

    def solve_whole(self):  #solves whole structure
        self.goButton.setEnabled(False)
        self.goButton.repaint()
        
        self.qclayers.populate_x_full()
        try:
            self.qclayers.solve_psi()
            self.plotDirty = True
            self.solveType = 'whole'
            self.update_quantumCanvas()
        except (IndexError,TypeError) as err:
            QMessageBox.warning(self, 'ErwinJr - Error', str(err))

        self.goButton.setEnabled(True)
        
    def solve(self):  #solves structure with basis
        self.solveButton.setEnabled(False)
        self.solveButton.repaint()
        
        try:
            self.dCL = ThePhysics.basisSolve(self.qclayers)
            ThePhysics.convert_dCL_to_data(self.qclayers, self.dCL)
            self.solveType = 'basis'        
            self.plotDirty = True
            self.update_quantumCanvas()
        except (ValueError,IndexError) as err:
            QMessageBox.warning(self,"ErwinJr - Error", str(err))
        
        self.solveButton.setEnabled(True)

    def pair_select(self, on):
        if on:
            self.wellSelectButton.setChecked(False)
            self.panButton.setChecked(False)
            self.zoomButton.setChecked(False)
            self.selectedWF = []
            self.stateHolder = []
            self.picker.setEnabled(True)
            self.quantumCanvas.canvas().setCursor(Qt.PointingHandCursor)
        else:
            self.picker.setEnabled(False)
            self.stateHolder = []
            self.pairSelectString.setText('')
            for curve in self.selectedWF:
                try:
                    curve.detach()
                except RuntimeError:
                    #underlying C/C++ object has been deleted
                    pass
            self.quantumCanvas.replot()
            self.quantumCanvas.canvas().setCursor(Qt.ArrowCursor)
            
    def hinput(self, aQPointF):
        #x data is self.qclayers.xPointsPost
        x = aQPointF.x()
        
        xLayerNum = argmin((self.qclayers.xPoints-x)**2)
        layerNum = self.qclayers.xLayerNums[xLayerNum]
        
        self.layerTable.selectRow(layerNum)
        self.layerTable.setFocus()

    def figure_of_merit(self):
        if len(self.stateHolder) < 2:
            return
        
        self.FoMButton.setEnabled(False)
        self.FoMButton.repaint()
        
        upper = self.stateHolder[-1]
        lower = self.stateHolder[-2]        
        if upper < lower:
            temp = upper
            upper = lower
            lower = temp
        
        self.tauUpper = 0; self.tauLower = 0
        for q in xrange(upper):
            self.tauUpper += 1/ThePhysics.lo_phonon_time(self.qclayers, upper, q)
        self.tauUpper = 1/self.tauUpper
        for q in xrange(lower):
            self.tauLower += 1/ThePhysics.lo_phonon_time(self.qclayers, lower, q)
        self.tauLower = 1/self.tauLower
        
        self.FoM = self.opticalDipole**2 * self.tauUpper * (1- self.tauLower/self.tauUpperLower)
        
        self.alphaISB = ThePhysics.alphaISB(self.qclayers, upper, lower)
        
        energyString  = u"<i>\u03C4<sub>upper</sub></i> : %6.2f ps<br><i>\u03C4<sub>lower</sub></i> : %6.2f ps" % (self.tauUpper, self.tauLower)
        energyString += u"<br>FoM: <b>%6.0f ps \u212B<sup>2</sup></b>" % (self.FoM)
        energyString += u"<br><i>\u03B1<sub>ISB</sub></i> : {0:.2f} cm<sup>2</sup>".format(self.alphaISB)

        self.pairSelectString.append(energyString)

        self.FoMButton.setEnabled(True)
        self.transferOpticalParametersButton.setEnabled(True)

    def ginput(self, aQPointF):
        #x data is self.qclayers.xPointsPost
        #y data is self.qclayers.xyPsiPsi[:,q] + self.qclayers.EigenE[q]
        
        x = aQPointF.x()
        y = aQPointF.y()
        
        xData = tile(self.qclayers.xPointsPost,(self.qclayers.xyPsiPsi.shape[1],1)).T
        yData = self.qclayers.xyPsiPsi + self.qclayers.EigenE
        
        xScale = self.quantumCanvas.axisScaleDiv(Qwt.QwtPlot.xBottom).upperBound() - self.quantumCanvas.axisScaleDiv(Qwt.QwtPlot.xBottom).lowerBound()
        yScale = self.quantumCanvas.axisScaleDiv(Qwt.QwtPlot.yLeft).upperBound() - self.quantumCanvas.axisScaleDiv(Qwt.QwtPlot.yLeft).lowerBound()

        r = nanmin(sqrt( ((xData-x)/xScale)**2 + ((yData-y)/yScale)**2 ), axis=0)
        selectedState = nanargmin(r)
        self.stateHolder.append(selectedState)
        
        q = selectedState
        mask = ~isnan(self.qclayers.xyPsiPsi[:,q])
        curve = SupportClasses.MaskedCurve(self.qclayers.xPointsPost,self.qclayers.xyPsiPsi[:,q] + self.qclayers.EigenE[q],mask)
        curve.setPen(QPen(Qt.black, 3))
        curve.setRenderHint(Qwt.QwtPlotItem.RenderAntialiased)
        curve.attach(self.quantumCanvas)
        self.selectedWF.append(curve)
        self.quantumCanvas.replot()
        
        if mod(len(self.stateHolder),2) == 0:
            self.FoMButton.setEnabled(True)
            E_i = self.qclayers.EigenE[self.stateHolder[-2]]
            E_j = self.qclayers.EigenE[self.stateHolder[-1]]
            if E_i > E_j:
                upper = self.stateHolder[-2]
                lower = self.stateHolder[-1]
            else:
                upper = self.stateHolder[-1]
                lower = self.stateHolder[-2]

            self.eDiff = 1000*(E_i-E_j)
            
            if self.solveType is 'basis':
                couplingEnergy = ThePhysics.coupling_energy(self.qclayers, self.dCL, upper, lower)
                self.transitionBroadening = ThePhysics.broadening_energy(self.qclayers, upper, lower)
                self.opticalDipole = ThePhysics.dipole(self.qclayers, upper, lower)            
                self.tauUpperLower = ThePhysics.lo_phonon_time(self.qclayers, upper, lower)
                energyString  = u"selected: %d, %d<br>energy diff: <b>%6.1f meV</b>\
                                 <br>coupling: %6.1f meV<br>broadening: %6.1f meV\
                                 <br>dipole: <b>%6.1f \u212B</b><br>LO scattering: <b>%6.2g ps</b>" % \
                               (self.stateHolder[-2], self.stateHolder[-1], self.eDiff, couplingEnergy, self.transitionBroadening, self.opticalDipole, self.tauUpperLower)

            elif self.solveType is 'whole':
                self.opticalDipole = ThePhysics.dipole(self.qclayers, upper, lower)
                self.tauUpperLower = ThePhysics.lo_phonon_time(self.qclayers, upper, lower)
                self.transitionBroadening = 0.1 * self.eDiff
                energyString = u"selected: %d, %d<br>energy diff: <b>%6.1f meV</b>\
                                 <br>dipole: %6.1f \u212B<br>LO scattering: %6.2g ps" % \
                               (self.stateHolder[-2], self.stateHolder[-1], self.eDiff, self.opticalDipole, self.tauUpperLower)
            else:
                self.FoMButton.setEnabled(False)
            
            self.pairSelectString.clear()
            self.pairSelectString.setText(energyString)        

    def transfer_optical_parameters(self):
        #set wavelength
        self.strata.wavelength = 1.24/abs(self.eDiff)*1000
        
        #set operating field
        self.strata.operatingField = self.qclayers.EField

        #set Lp
        LpFirst = self.LpFirstSpinbox.value()
        LpLast = self.LpLastSpinbox.value()+1 #+1 because range is not inclusive of last value
        self.strata.Lp = sum(self.qclayers.layerWidths[LpFirst:LpLast])
        
        #set nD doping sheet density
        self.strata.nD = sum(self.qclayers.layerDopings[LpFirst:LpLast]*self.qclayers.layerWidths[LpFirst:LpLast]) / sum(self.qclayers.layerWidths[LpFirst:LpLast])
        
        #set aCore
        self.strata.aCore = self.alphaISB
        
        #set nCore
        self.strata.nCore = self.strata.get_nCore(self.qclayers)
        
        #set tauUpper
        self.strata.tauUpper = self.tauUpper
        
        #set tauLower
        self.strata.tauLower = self.tauLower
        
        #set tauUpperLower
        self.strata.tauUpperLower = self.tauUpperLower
        
        #set optical dipole
        self.strata.opticalDipole = self.opticalDipole
        
        #set figure of merit
        self.strata.FoM = self.FoM
        
        #2gamma transition broadening
        self.strata.transitionBroadening = self.transitionBroadening / 1000 #store in eV
        
        #GUI settings
        self.transferOpticalParametersButton.setEnabled(False)
        self.editOpticalParametersBox.setChecked(False)
        
        #update all the input boxes
        self.update_stratum_inputBoxes()
        
        #update the stratumTable
        self.stratumTable_refresh()




#===============================================================================
# Quantum Tab Plotting and Plot Control
#===============================================================================

    def update_quantumCanvas(self):
        #setup for layer outline

        #PyQwt code
        if self.plotDirty: #self.plotdirty is True when self.go is executed
            self.quantumCanvas.clear()
        
        #plot xVc
        try:
            self.curveVc.detach()
            self.curveAR.detach()
        except AttributeError:
            #self.curveVc has not yet been formed
            pass
        except RuntimeError:
            #self.curveVc deleted with self.quantumCanvas.clear()
            pass
        self.curveVc = Qwt.QwtPlotCurve()
        self.curveVc.setData(self.qclayers.xPoints,self.qclayers.xVc)
        self.curveVc.setPen(QPen(Qt.black, 1))
        if settings.antialiased:
            self.curveVc.setRenderHint(Qwt.QwtPlotItem.RenderAntialiased)
        self.curveVc.attach(self.quantumCanvas)
        
        #plot Conduction Band L-Valley
        if self.plotVL:
            self.curveVL = Qwt.QwtPlotCurve()
            self.curveVL.setData(self.qclayers.xPoints,self.qclayers.xVL)
            self.curveVL.setPen(QPen(Qt.green, 1, Qt.DashLine))
            if settings.antialiased:
                self.curveVL.setRenderHint(Qwt.QwtPlotItem.RenderAntialiased)
            self.curveVL.attach(self.quantumCanvas)
        else:
            try:
                self.curveVL.detach()
            except (AttributeError, RuntimeError):
                pass
        
        #plot Conduction Band X-Valley
        if self.plotVX:
            self.curveVX = Qwt.QwtPlotCurve()
            self.curveVX.setData(self.qclayers.xPoints,self.qclayers.xVX)
            self.curveVX.setPen(QPen(Qt.magenta, 1, Qt.DashDotLine))
            if settings.antialiased:
                self.curveVX.setRenderHint(Qwt.QwtPlotItem.RenderAntialiased)
            self.curveVX.attach(self.quantumCanvas)
        else:
            try:
                self.curveVX.detach()
            except (AttributeError, RuntimeError):
                pass
            
        #plot Light Hole Valence Band
        if self.plotLH:
            self.curveLH = Qwt.QwtPlotCurve()
            self.curveLH.setData(self.qclayers.xPoints,self.qclayers.xVLH)
            self.curveLH.setPen(QPen(Qt.black, 1))
            if settings.antialiased:
                self.curveLH.setRenderHint(Qwt.QwtPlotItem.RenderAntialiased)
            self.curveLH.attach(self.quantumCanvas)
        else:
            try:
                self.curveLH.detach()
            except (AttributeError, RuntimeError):
                pass
            
        #plot Split Off Valence Band
        if self.plotSO:
            self.curveSO = Qwt.QwtPlotCurve()
            self.curveSO.setData(self.qclayers.xPoints,self.qclayers.xVSO)
            self.curveSO.setPen(QPen(Qt.red, 1, Qt.DashLine))
            if settings.antialiased:
                self.curveSO.setRenderHint(Qwt.QwtPlotItem.RenderAntialiased)
            self.curveSO.attach(self.quantumCanvas)
        else:
            try:
                self.curveSO.detach()
            except (AttributeError, RuntimeError):
                pass

        #highlight selected layer & make AR layers bold
        try:
            self.curveSelection.detach()
        except AttributeError:
            #self.curveSelection has not yet been formed
            pass       
        except RuntimeError:
            #self.curveSelection deleted with self.quantumCanvas.clear()        
            pass
        if self.qclayers.layerSelected >= 0 and self.qclayers.layerSelected < self.qclayers.layerWidths.size:
            mask = ~isnan(self.qclayers.xARs)
            self.curveAR =  SupportClasses.MaskedCurve(self.qclayers.xPoints,self.qclayers.xARs,mask)
            self.curveAR.setPen(QPen(Qt.black, 2))
            if settings.antialiased:
                self.curveAR.setRenderHint(Qwt.QwtPlotItem.RenderAntialiased)
            self.curveAR.attach(self.quantumCanvas)
            
            mask = ~isnan(self.qclayers.xLayerSelected)
            self.curveSelection = SupportClasses.MaskedCurve(self.qclayers.xPoints,self.qclayers.xLayerSelected,mask)
            self.curveSelection.setPen(QPen(Qt.blue, 1.5))
            if settings.antialiased:
                self.curveSelection.setRenderHint(Qwt.QwtPlotItem.RenderAntialiased)
            self.curveSelection.attach(self.quantumCanvas)
        
        #plot wavefunctions
        if self.plotDirty and hasattr(self.qclayers, 'EigenE'):
            self.plotDirty=False
            self.curveWF = []
            for q in xrange(self.qclayers.EigenE.size):
                mask = ~isnan(self.qclayers.xyPsiPsi[:,q])
                curve = SupportClasses.MaskedCurve(self.qclayers.xPointsPost,self.qclayers.xyPsiPsi[:,q] + self.qclayers.EigenE[q],mask)
                r,g,b = self.colors[mod(q,13)]
                curve.setPen(QPen(QColor(r,g,b), 2))
                curve.setRenderHint(Qwt.QwtPlotItem.RenderAntialiased)
                self.curveWF.append(curve)
                curve.attach(self.quantumCanvas)
            
        self.quantumCanvas.setAxisTitle(Qwt.QwtPlot.xBottom, u'Position (\u212B)')
        self.quantumCanvas.setAxisTitle(Qwt.QwtPlot.yLeft, 'Energy (eV)')
        self.quantumCanvas.replot()
        self.zoomer.setZoomBase()
        
    def create_zoomer(self):
        self.setContextMenuPolicy(Qt.NoContextMenu)
        self.zoomer = Qwt.QwtPlotZoomer(
            Qwt.QwtPlot.xBottom,
            Qwt.QwtPlot.yLeft,
            Qwt.QwtPicker.DragSelection,
            Qwt.QwtPicker.AlwaysOff,
            self.quantumCanvas.canvas())
        self.zoomer.setEnabled(False)
        #self.zoomer.setRubberBandPen(QPen(Qt.green))
        pattern = [
            Qwt.QwtEventPattern.MousePattern(Qt.LeftButton,
                                             Qt.NoModifier),
            Qwt.QwtEventPattern.MousePattern(Qt.MidButton,
                                             Qt.NoModifier),
            Qwt.QwtEventPattern.MousePattern(Qt.RightButton,
                                             Qt.NoModifier),
            Qwt.QwtEventPattern.MousePattern(Qt.LeftButton,
                                             Qt.ShiftModifier),
            Qwt.QwtEventPattern.MousePattern(Qt.MidButton,
                                             Qt.ShiftModifier),
            Qwt.QwtEventPattern.MousePattern(Qt.RightButton,
                                             Qt.ShiftModifier),
            ]
        self.zoomer.setMousePattern(pattern)
        
        self.picker = Qwt.QwtPlotPicker(
            Qwt.QwtPlot.xBottom,
            Qwt.QwtPlot.yLeft,
            Qwt.QwtPicker.PointSelection | Qwt.QwtPicker.DragSelection,
            Qwt.QwtPlotPicker.CrossRubberBand,
            Qwt.QwtPicker.AlwaysOn,
            self.quantumCanvas.canvas())
        self.picker.setRubberBandPen(QPen(Qt.green))
        self.picker.setTrackerPen(QPen(Qt.black))
        self.picker.connect(self.picker, SIGNAL('selected(const QwtDoublePoint&)'), self.ginput)
        self.picker.setEnabled(False)
        
        self.picker2 = Qwt.QwtPlotPicker(
            Qwt.QwtPlot.xBottom,
            Qwt.QwtPlot.yLeft,
            Qwt.QwtPicker.PointSelection | Qwt.QwtPicker.DragSelection,
            Qwt.QwtPlotPicker.NoRubberBand,
            Qwt.QwtPicker.AlwaysOff,
            self.quantumCanvas.canvas())
        #self.picker2.setRubberBandPen(QPen(Qt.green))
        self.picker2.setTrackerPen(QPen(Qt.black))
        self.picker2.connect(self.picker2, SIGNAL('selected(const QwtDoublePoint&)'), self.hinput)
        self.picker2.setEnabled(False)
        
        self.panner = Qwt.QwtPlotPanner(self.quantumCanvas.canvas())
        self.panner.setEnabled(False)
        
        self.zoom(False)

    def well_select(self, on):
        if on:
            self.pairSelectButton.setChecked(False)
            self.panButton.setChecked(False)
            self.zoomButton.setChecked(False)
            self.picker2.setEnabled(True)
            self.quantumCanvas.canvas().setCursor(Qt.PointingHandCursor)
        else:
            self.picker2.setEnabled(False)
            self.quantumCanvas.canvas().setCursor(Qt.ArrowCursor)

    def zoom(self, on):
        if on:
            self.pairSelectButton.setChecked(False)
            self.wellSelectButton.setChecked(False)
            self.panButton.setChecked(False)
            self.zoomer.setEnabled(True)
            self.quantumCanvas.canvas().setCursor(Qt.CrossCursor)
        else:
            self.zoomer.setEnabled(False)
            self.quantumCanvas.canvas().setCursor(Qt.ArrowCursor)

    def zoom_out(self):
        """Auto scale and clear the zoom stack
        """
        self.pairSelectButton.setChecked(False)
        self.wellSelectButton.setChecked(False)
        self.zoomButton.setChecked(False)
        self.panButton.setChecked(False)
        self.quantumCanvas.setAxisAutoScale(Qwt.QwtPlot.xBottom)
        self.quantumCanvas.setAxisAutoScale(Qwt.QwtPlot.yLeft)
        self.quantumCanvas.replot()
        self.zoomer.setZoomBase()

    def pan(self, on):
        if on:
            self.pairSelectButton.setChecked(False)
            self.wellSelectButton.setChecked(False)
            self.zoomButton.setChecked(False)
            self.quantumCanvas.setCursor(Qt.OpenHandCursor)
            self.panner.setEnabled(True)
            self.quantumCanvas.canvas().setCursor(Qt.OpenHandCursor)
        else:
            self.panner.setEnabled(False)
            self.quantumCanvas.canvas().setCursor(Qt.ArrowCursor)

    def clear_WFs(self):
        self.plotDirty = False
        delattr(self.qclayers,'EigenE')
        self.quantumCanvas.clear()
        self.update_quantumCanvas()




#===============================================================================
# General Menu Functions
#===============================================================================

    def change_main_tab(self, tabIdx):
        self.menuBar().clear()
        if tabIdx == 0:
            self.create_Quantum_menu()
        elif tabIdx == 1:
            self.create_Optical_menu()
        elif tabIdx == 2:
            pass
        else:
            assert 1==2

    def add_actions(self, target, actions):
        for action in actions:
            if action is None:
                target.addSeparator()
            else:
                target.addAction(action)

    def create_action(self, text, slot=None, shortcut=None, 
                        icon=None, tip=None, checkable=False, 
                        signal="triggered()", ischecked=False):
        action = QAction(text, self)
        if icon is not None:
            action.setIcon(QIcon("images/%s.png" % icon))
        if shortcut is not None:
            action.setShortcut(shortcut)
        if tip is not None:
            action.setToolTip(tip)
            action.setStatusTip(tip)
        if slot is not None:
            self.connect(action, SIGNAL(signal), slot)
        if checkable:
            action.setCheckable(True)
        if ischecked:
            action.setChecked(True)
        return action

    def create_Quantum_menu(self):
        #file menu
        self.file_menu = self.menuBar().addMenu("&File")
        newFileAction      = self.create_action("&New...", self.fileNew, QKeySequence.New, "filenew", "New ErwinJr file")
        openFileAction     = self.create_action("&Open", shortcut="Ctrl+O", slot=self.fileOpen, tip="Open ErwinJr file", icon="fileopen")
        saveFileAction     = self.create_action("&Save", shortcut="Ctrl+S", slot=self.fileSave, tip="Save ErwinJr file", icon="filesave")
        saveAsFileAction   = self.create_action("S&ave As", shortcut="Ctrl+W", slot=self.fileSaveAs, tip="Save ErwinJr file as", icon="filesaveas")
        exportQuantumCanvasAction = self.create_action("Export Band Diagram Image", slot=self.export_quantumCanvas, tip="Export Band Diagram Image")
        exportBandCSVAction = self.create_action("Export Band Diagram Data", slot=self.export_band_diagram_data, tip="Export Band Diagram Data")
        quit_action = self.create_action("&Quit", slot=self.close, shortcut="Ctrl+Q", tip="Close the application", icon="filequit")
        self.fileMenuActions = (newFileAction, openFileAction, saveFileAction, saveAsFileAction, None, exportBandCSVAction, exportQuantumCanvasAction, None, quit_action)
        self.connect(self.file_menu, SIGNAL("aboutToShow()"), self.updateFileMenu)
        #self.add_actions(self.file_menu, (newFileAction, openFileAction, saveFileAction, saveAsFileAction, None, quit_action))

        #edit menu
        self.edit_menu = self.menuBar().addMenu("&Edit")
        temperatureAction = self.create_action("&Temperature", slot=self.set_temperature, tip="Set temperature")
        bumpLayerAction = self.create_action("&Bump First Layer", slot=self.bump_first_layer, tip="Move zeroth layer to first layer")
        copyStructureAction = self.create_action("&Copy Structure", slot=self.copy_structure, tip="Copy Layer Structure to Clipboard")
        self.add_actions(self.edit_menu, (temperatureAction, bumpLayerAction, None, copyStructureAction))

        #view menu
        self.view_menu = self.menuBar().addMenu("&View")
        VXBandAction = self.create_action("X Valley Conduction Band", checkable=True, ischecked=self.plotVX, slot=self.view_VXBand)
        VLBandAction = self.create_action("L Valley Conduction Band", checkable=True, ischecked=self.plotVL, slot=self.view_VLBand)
        LHBandAction = self.create_action("Light Hole Valence Band", checkable=True, ischecked=self.plotLH, slot=self.view_LHBand)
        SOBandAction = self.create_action("Split Off Valence Band", checkable=True, ischecked=self.plotSO, slot=self.view_SOBand)
        self.add_actions(self.view_menu, (VXBandAction,VLBandAction,LHBandAction,SOBandAction))        

        #help menu
        self.help_menu = self.menuBar().addMenu("&Help")
        about_action = self.create_action("&About",shortcut='F1', slot=self.on_about)
        licenses_action = self.create_action("&License", slot=self.on_licenses)
        tutorialAction = self.create_action("&Tutorial", slot=self.on_tutorial)
        self.add_actions(self.help_menu, (tutorialAction,about_action,licenses_action))
        
    def create_Optical_menu(self):
        #file menu
        self.file_menu = self.menuBar().addMenu("&File")
        newFileAction      = self.create_action("&New...", self.fileNew, QKeySequence.New, "filenew", "New ErwinJr file")
        openFileAction     = self.create_action("&Open", shortcut="Ctrl+O", slot=self.fileOpen, tip="Open ErwinJr file", icon="fileopen")
        saveFileAction     = self.create_action("&Save", shortcut="Ctrl+S", slot=self.fileSave, tip="Save ErwinJr file", icon="filesave")
        saveAsFileAction   = self.create_action("S&ave As", shortcut="Ctrl+W", slot=self.fileSaveAs, tip="Save ErwinJr file as", icon="filesaveas")
        quit_action = self.create_action("&Quit", slot=self.close, shortcut="Ctrl+Q", tip="Close the application", icon="filequit")
        self.fileMenuActions = (newFileAction, openFileAction, saveFileAction, saveAsFileAction, None, quit_action)
        self.connect(self.file_menu, SIGNAL("aboutToShow()"), self.updateFileMenu)
        #self.add_actions(self.file_menu, (newFileAction, openFileAction, saveFileAction, saveAsFileAction, None, quit_action))

        #help menu
        self.help_menu = self.menuBar().addMenu("&Help")
        about_action = self.create_action("&About",shortcut='F1', slot=self.on_about)
        licenses_action = self.create_action("&License", slot=self.on_licenses)
        tutorialAction = self.create_action("&Tutorial", slot=self.on_tutorial)
        self.add_actions(self.help_menu, (tutorialAction,about_action,licenses_action))




#===============================================================================
# File Menu Items
#===============================================================================

    def updateFileMenu(self):
        self.file_menu.clear()
        self.add_actions(self.file_menu, self.fileMenuActions[:-1])
        current = (QString(self.filename)
                   if self.filename is not None else None)
        recentFiles = []
        for fname in self.recentFiles:
            if fname != current and QFile.exists(fname):
                recentFiles.append(fname)
        if recentFiles:
            self.file_menu.addSeparator()
            for i, fname in enumerate(recentFiles):
                action = QAction(
                        "&{0}  {1}".format(i + 1, QFileInfo(
                        fname).fileName()), self)
                action.setData(QVariant(fname))
                self.connect(action, SIGNAL("triggered()"),
                             partial(self.fileOpen,fname))
                self.file_menu.addAction(action)
        self.file_menu.addSeparator()
        self.file_menu.addAction(self.fileMenuActions[-1])

    def addRecentFile(self, fname):
        if fname is None:
            return
        if not self.recentFiles.contains(fname):
            self.recentFiles.prepend(QString(fname))
            while self.recentFiles.count() > 9:
                self.recentFiles.takeLast()
        
    def loadInitialFile(self):
        qsettings = QSettings()
        fname = unicode(qsettings.value("LastFile").toString())
        if fname and QFile.exists(fname):
            if fname.split('.')[-1] == 'qcl':
                self.qclLoad(fname)

            self.zoomer.zoom(0)
            self.quantumCanvas.clear()
            self.update_Lp_limits()
            self.update_inputBoxes()
            self.layerTable_refresh()
            self.qclayers.populate_x()

            self.strata.populate_rIndexes()
            self.update_stratum_inputBoxes()
            self.stratumTable_refresh()
            self.opticalCanvas.clear()
            self.update_opticalCanvas()
        
        self.filename = fname
        self.addRecentFile(fname)
        self.dirty = False
        self.update_windowTitle()       

    def fileNew(self):
        if not self.okToContinue():
            return False
            
        self.filename = None
        self.plotDirty = False
        self.quantumCanvas.clear()
        self.opticalCanvas.clear()
        self.optimization1DCanvas.clear()
        
        self.qclayers = ThePhysics.QCLayers()
        self.strata = ThePhysics.Strata()
        
        self.zoomer.zoom(0)

        self.update_Lp_limits()
        self.update_inputBoxes()
        self.layerTable_refresh()
        self.qclayers.populate_x()

        self.layerTable_refresh()
        self.layerTable.selectRow(1)
        self.layerTable.setFocus()
        
        self.update_stratum_inputBoxes()
        
        self.update_opticalCanvas()
        self.stratumTable_refresh()
        
        self.dirty = False
        self.update_windowTitle()

        return True

    def okToContinue(self):
        if self.dirty:
            reply = QMessageBox.question(self,
                                         "ErwinJr " + str(majorVersion) + " - Unsaved Changes",
                                         "Save unsaved changes?",
                                         QMessageBox.Yes|QMessageBox.No|QMessageBox.Cancel)
            if reply == QMessageBox.Cancel:
                return False
            elif reply == QMessageBox.Yes:
                self.fileSave()
        return True

    def update_windowTitle(self):
        if self.filename is not None:
            self.setWindowTitle("ErwinJr " + str(majorVersion) + " - %s[*]" % os.path.basename(str(self.filename)))
        else:
            self.setWindowTitle("ErwinJr " + str(majorVersion) + "[*]")
        self.setWindowModified(self.dirty)

    def fileOpen(self, fname = None):
        #if not self.okToContinue():
        #    return False
        
        #clear all old data, also calls self.okToContinue()
        if not self.fileNew(): 
            return False
        if fname is None:
            dir = os.path.dirname(str(self.filename)) if self.filename is not None else "."
            fname =unicode(QFileDialog.getOpenFileName(self,"ErwinJr - Choose file", 
                                                    dir, "ErwinJr files (*.qcl)\nAll files (*.*)"))
        #open file and determine if it is from the Matlab version of ErwinJr
        filehandle = open(fname, 'r')
        firstLine = filehandle.readline()
        filehandle.close()
        if fname:
            if firstLine.split(':')[0] == 'Description':
                self.qclPtonLoad(fname)
            elif firstLine == 'ErwinJr Data File'+self.newLine:
                self.qclLoad(fname)
            else:
                QMessageBox.warning(self,'ErwinJr Error','Could not recognize input file.')
                return
            self.zoomer.zoom(0)
            self.quantumCanvas.clear()
            self.update_Lp_limits()
            self.update_inputBoxes()
            self.layerTable_refresh()
            self.qclayers.populate_x()
            
            #if firstLine == 'ErwinJr Data File\n':
            self.strata.populate_rIndexes()
            self.update_stratum_inputBoxes()
            self.stratumTable_refresh()
            self.opticalCanvas.clear()
            self.update_opticalCanvas()
        
        self.filename = fname
        self.addRecentFile(fname)
        self.dirty = False
        self.update_windowTitle()
        
        return True

    def qclPtonLoad(self, fname):
        try:
            filehandle = open(fname, 'r')
            self.qclayers.description = filehandle.readline().split(':')[1].strip()
            self.qclayers.EField      = float(filehandle.readline().split(':')[1].strip())
            self.qclayers.xres        = float(filehandle.readline().split(':')[1].strip())
            self.qclayers.vertRes     = float(filehandle.readline().split(':')[1].strip())
            self.qclayers.moleFrac1   = float(filehandle.readline().split(':')[1].strip())
            self.qclayers.moleFrac2   = float(filehandle.readline().split(':')[1].strip())
            self.qclayers.moleFrac3   = float(filehandle.readline().split(':')[1].strip())
            self.qclayers.moleFrac4   = float(filehandle.readline().split(':')[1].strip())
            self.qclayers.solver      = filehandle.readline().split(':')[1].strip()
            self.qclayers.Temperature = float(filehandle.readline().split(':')[1].strip())
            self.qclayers.TempFoM     = float(filehandle.readline().split(':')[1].strip())
            self.qclayers.repeats     = int(filehandle.readline().split(':')[1].strip())
            self.qclayers.diffLength  = float(filehandle.readline().split(':')[1].strip())
            
            filehandle.readline() #throw the column description line away
            lines = filehandle.readlines();
            rows = len(lines)
            variables = ['layerWidths', 'layerBarriers', 'layerARs', 'layerDopings', 
                         'layerMaterials', 'layerDividers']
            for item in variables:
                setattr(self.qclayers, item, zeros(rows))
            for q, line in enumerate(lines):
                line = line.split('\t')
                self.qclayers.layerWidths[q]    = float(line[1])
                self.qclayers.layerBarriers[q]  = float(line[2])
                self.qclayers.layerARs[q]       = float(line[3])
                self.qclayers.layerMaterials[q] = float(line[4])
                self.qclayers.layerDopings[q]   = float(line[5])
                self.qclayers.layerDividers[q]  = float(line[6])
            self.qclayers.layerMaterials[nonzero(self.qclayers.layerMaterials == 4)[0]] = 1
            self.qclayers.layerMaterials[nonzero(self.qclayers.layerMaterials == 5)[0]] = 2
            
            filehandle.close()
            
            self.bump_first_layer()
        except Exception as err:
            QMessageBox.warning(self,"ErwinJr - Warning",
                             "Could not load *.qcl file.\n"+str(err))

    def qclLoad(self, fname):
        try:
            valDict = {}
            filehandle = open(fname, 'r')
            filehandle.readline() #throw away 'ErwinJr Data File'
            while True:
                line = filehandle.readline()
                if line == '# QC layers #'+self.newLine:
                    break
                line = line.split(':')
                valDict[line[0]] = line[1].strip()
            
            self.qclayers.description = valDict['Description']
            self.qclayers.substrate   = valDict['Substrate']
            self.qclayers.EField      = float(valDict['Efield'])
            self.qclayers.xres        = float(valDict['xres'])
            self.qclayers.vertRes     = float(valDict['Eres'])
            self.qclayers.moleFrac1   = float(valDict['moleFrac1'])
            self.qclayers.moleFrac2   = float(valDict['moleFrac2'])
            self.qclayers.moleFrac3   = float(valDict['moleFrac3'])
            self.qclayers.moleFrac4   = float(valDict['moleFrac4'])
            self.qclayers.moleFrac5   = float(valDict['moleFrac5'])
            self.qclayers.moleFrac6   = float(valDict['moleFrac6'])
            self.qclayers.moleFrac7   = float(valDict['moleFrac7'])
            self.qclayers.moleFrac8   = float(valDict['moleFrac8'])
            self.qclayers.solver      = valDict['Solver']
            self.qclayers.Temperature = float(valDict['Temp'])
            self.qclayers.TempFoM     = float(valDict['TempFoM'])
            self.qclayers.repeats     = int(valDict['PlotPeriods'])
            self.qclayers.diffLength  = float(valDict['DiffLeng'])
            
            self.strata.wavelength           = float(valDict['Wavelength'])
            self.strata.operatingField       = float(valDict['StratumField'])
            self.strata.Lp                   = float(valDict['Lp'])
            self.strata.Np                   = float(valDict['Np'])
            self.strata.aCore                = float(valDict['alphaCore'])
            self.strata.nCore                = complex(valDict['nCore'])
            self.strata.nD                   = float(valDict['nD'])
            self.strata.transitionBroadening = float(valDict['transitionBroadening'])        
            self.strata.tauUpper             = float(valDict['tauUpper'])
            self.strata.tauLower             = float(valDict['tauLower'])
            self.strata.tauUpperLower        = float(valDict['tauUpperLower'])
            self.strata.opticalDipole        = float(valDict['opticalDipole'])
            self.strata.FoM                  = float(valDict['FoM'])
            self.strata.waveguideFacets      = valDict['waveguideFacets']
            self.strata.waveguideLength      = float(valDict['waveguideLength'])
            self.strata.customFacet          = float(valDict['customFacet'])
            
            lines = []
            while True:
                line = filehandle.readline()
                if line == '# Optical strata #'+self.newLine:
                    break
                lines.append(line)
            rows = len(lines)
            variables = ['layerWidths', 'layerBarriers', 'layerARs', 'layerDopings', 
                         'layerMaterials', 'layerDividers']
            for item in variables:
                setattr(self.qclayers, item, zeros(rows))
            for q, line in enumerate(lines):
                line = line.split('\t')
                self.qclayers.layerWidths[q]    = float(line[1])
                self.qclayers.layerBarriers[q]  = float(line[2])
                self.qclayers.layerARs[q]       = float(line[3])
                self.qclayers.layerMaterials[q] = float(line[4])
                self.qclayers.layerDopings[q]   = float(line[5])
                self.qclayers.layerDividers[q]  = float(line[6])
            
            lines = filehandle.readlines()
            rows = len(lines)
            variables = ['stratumCompositions', 'stratumThicknesses', 'stratumDopings']
            for item in variables:
                setattr(self.strata, item, zeros(rows))
            self.strata.stratumMaterials = []
            for q, line in enumerate(lines):
                line = line.split('\t')
                self.strata.stratumMaterials.append(str(line[1]))
                self.strata.stratumCompositions[q] = float(line[2])
                self.strata.stratumThicknesses[q] = float(line[3])
                self.strata.stratumDopings[q]     = float(line[4])
            
            filehandle.close()
        except Exception as err:
            QMessageBox.warning(self,"ErwinJr - Warning",
                             "Could not load *.qcl file.\n"+str(err))

    def fileSave(self):
        if self.filename is None:
            return self.fileSaveAs()
        else:
            #os.path.extsep
            if self.filename.split('.')[-1] == 'qcl':
                if self.qclSave(self.filename):
                    self.dirty = False
                    self.update_windowTitle()
                    return True
                else:
                    return False
            else:
                raise IOError('The *.' + self.filename.split('.')[-1] + ' extension is not supported.')
                return False

    def fileSaveAs(self):
        fname = self.filename if self.filename is not None else "."
        typeString = "ErwinJr 2.x file (*.qcl)\nAll files (*.*)"
        fname = unicode(QFileDialog.getSaveFileName(self,"ErwinJr - Save File", QString(fname), typeString))
        if fname:
            if "." not in fname:
                fname += ".qcl"
            self.addRecentFile(fname)
            self.filename = fname
            return self.fileSave()
        return False

    def qclSave(self, fname):
        filehandle = open(fname, 'w')
        filehandle.write("ErwinJr Data File\n")
        filehandle.write("Version:" + str(ejVersion) + '\n')
        filehandle.write("Description:" + self.qclayers.description + '\n')
        filehandle.write("Substrate:" + self.qclayers.substrate + '\n')
        filehandle.write("Efield:" + str(self.qclayers.EField) + '\n')
        filehandle.write("xres:" + str(self.qclayers.xres) + '\n')
        filehandle.write("Eres:" + str(self.qclayers.vertRes) + '\n')
        filehandle.write("moleFrac1:" + str(self.qclayers.moleFrac1) + '\n')
        filehandle.write("moleFrac2:" + str(self.qclayers.moleFrac2) + '\n')
        filehandle.write("moleFrac3:" + str(self.qclayers.moleFrac3) + '\n')
        filehandle.write("moleFrac4:" + str(self.qclayers.moleFrac4) + '\n')
        filehandle.write("moleFrac5:" + str(self.qclayers.moleFrac5) + '\n')
        filehandle.write("moleFrac6:" + str(self.qclayers.moleFrac6) + '\n')
        filehandle.write("moleFrac7:" + str(self.qclayers.moleFrac7) + '\n')
        filehandle.write("moleFrac8:" + str(self.qclayers.moleFrac8) + '\n')
        filehandle.write("Solver:" + self.qclayers.solver + '\n')
        filehandle.write("Temp:" + str(self.qclayers.Temperature) + '\n')
        filehandle.write("TempFoM:" + str(self.qclayers.TempFoM) + '\n')
        filehandle.write("PlotPeriods:" + str(self.qclayers.repeats) + '\n')
        filehandle.write("DiffLeng:" + str(self.qclayers.diffLength) + '\n')
        
        filehandle.write("Wavelength:" + str(self.strata.wavelength) + '\n')
        filehandle.write("StratumField:" + str(self.strata.operatingField) + '\n')
        filehandle.write("Lp:" + str(self.strata.Lp) + '\n')
        filehandle.write("Np:" + str(self.strata.Np) + '\n')
        filehandle.write("alphaCore:" + str(self.strata.aCore) + '\n')
        filehandle.write("nCore:" + str(self.strata.nCore) + '\n')
        filehandle.write("nD:" + str(self.strata.nD) + '\n')
        filehandle.write("transitionBroadening:" + str(self.strata.transitionBroadening) + '\n')
        filehandle.write("tauUpper:" + str(self.strata.tauUpper) + '\n')
        filehandle.write("tauLower:" + str(self.strata.tauLower) + '\n')
        filehandle.write("tauUpperLower:" + str(self.strata.tauUpperLower) + '\n')
        filehandle.write("opticalDipole:" + str(self.strata.opticalDipole) + '\n')
        filehandle.write("FoM:" + str(self.strata.FoM) + '\n')
        filehandle.write("waveguideFacets:" + self.strata.waveguideFacets + '\n')
        filehandle.write("waveguideLength:" + str(self.strata.waveguideLength) + '\n')
        filehandle.write("customFacet:" + str(self.strata.customFacet) + '\n')
        
        filehandle.write("# QC layers #\n")
        for row in xrange(self.qclayers.layerWidths.size):
            string = "%d\t%f\t%d\t%d\t%d\t%f\t%d\n" % (row+1, self.qclayers.layerWidths[row], 
                      self.qclayers.layerBarriers[row], self.qclayers.layerARs[row], 
                      self.qclayers.layerMaterials[row], self.qclayers.layerDopings[row], 
                      self.qclayers.layerDividers[row])
            filehandle.write(string)
            
        filehandle.write("# Optical strata #\n")
        for row in xrange(self.strata.stratumDopings.size):
            string = "%d\t%s\t%f\t%f\t%f\n" % (row+1, self.strata.stratumMaterials[row], 
                      self.strata.stratumCompositions[row], self.strata.stratumThicknesses[row], 
                      self.strata.stratumDopings[row])
            filehandle.write(string)
            
        filehandle.close()
        return True

    def qclPtonSave(self, fname):
        filehandle = open(fname, 'w')
        filehandle.write("Description:" + self.qclayers.description + '\n')
        filehandle.write("Efield:" + str(self.qclayers.EField) + '\n')
        filehandle.write("xres:" + str(self.qclayers.xres) + '\n')
        filehandle.write("Eres:" + str(self.qclayers.vertRes) + '\n')
        filehandle.write("InGaAsx:" + str(self.qclayers.moleFrac1) + '\n')
        filehandle.write("AlInAsx:" + str(self.qclayers.moleFrac2) + '\n')
        filehandle.write("InGaAsx2:" + str(self.qclayers.moleFrac3) + '\n')
        filehandle.write("AlInAsx2:" + str(self.qclayers.moleFrac4) + '\n')
        filehandle.write("Solver:" + self.qclayers.solver + '\n')
        filehandle.write("Temp:" + str(self.qclayers.Temp) + '\n')
        filehandle.write("TempFoM:" + str(self.qclayers.TempFoM) + '\n')
        filehandle.write("PlotPeriods:" + str(self.qclayers.repeats) + '\n')
        filehandle.write("DiffLeng:" + str(self.qclayers.diffLength) + '\n')
        
        filehandle.write("regionNum\twellWdiths\tbarrierSwitch\tarSwitch\tmaterial\tdoping\tdivider\n")
        for row in xrange(self.qclayers.layerWidths.size):
            string = "%d\t%f\t%d\t%d\t%d\t%f\t%d\n" % (row+1, self.qclayers.layerWidths[row], 
                      self.qclayers.layerBarriers[row], self.qclayers.layerARs[row], 
                      self.qclayers.layerMaterials[row], self.qclayers.layerDopings[row], 
                      self.qclayers.layerDividers[row])
            filehandle.write(string)
            
        filehandle.close()
        return True

    def closeEvent(self, event):
        if self.okToContinue():
            qsettings = QSettings()
            filename = QVariant(QString(self.filename)) if self.filename is not None else QVariant()
            qsettings.setValue("LastFile", filename)
            recentFiles = QVariant(self.recentFiles) if self.recentFiles else QVariant()
            qsettings.setValue("RecentFiles", recentFiles)
            qsettings.setValue("MainWindow/Geometry", QVariant(self.saveGeometry()))
            qsettings.setValue("MainWindow/State", QVariant(self.saveState()))
        else:
            event.ignore()




#===============================================================================
# Export Functions
#===============================================================================

    def export_quantumCanvas(self):
        fname = unicode(QFileDialog.getSaveFileName(self,"ErwinJr - Export Band Structure Image",self.filename.split('.')[0],
                "Portable Network Graphics file (*.png)"))
        if not fname:
            return
        
        try:
            self.curveAR.detach()
            self.curveSelection.detach()
            self.quantumCanvas.replot()
        except:
            pass

        #set background color to white and save presets
        bgColor = self.quantumCanvas.canvasBackground()
        bgRole = self.mainTabWidget.backgroundRole()
        self.mainTabWidget.setBackgroundRole(QPalette.Base)
        self.quantumWidget.setBackgroundRole(QPalette.Base)
        self.quantumCanvas.setCanvasBackground(Qt.white)
        self.quantumCanvas.setAutoFillBackground(True)

        #save image
        QPixmap.grabWidget(self.quantumCanvas).save(fname+'.png', 'PNG')

        self.quantumWidget.setBackgroundRole(QPalette.Window)
        self.quantumCanvas.setCanvasBackground(bgColor)
        self.mainTabWidget.setBackgroundRole(bgRole)
        
        try:
            self.curveAR.attach(self.quantumCanvas)
            self.curveSelection.attach(self.quantumCanvas)
            self.quantumCanvas.replot()
        except:
            pass
        
    def export_band_diagram_data(self):
        fname = unicode(QFileDialog.getSaveFileName(self,"ErwinJr - Export Band Structure Data",self.filename.split('.')[0],
                "Comma-Separated Value file (*.csv)"))
        if fname != '': #if user doesn't click cancel
            savetxt(fname.split('.')[0] + '_CB' + '.csv', column_stack([self.qclayers.xPoints,self.qclayers.xVc]), delimiter=',')
        
            try: self.qclayers.xyPsiPsi
            except AttributeError: pass #band structure hasn't been solved yet
            else:
                xyPsiPsiEig = zeros(self.qclayers.xyPsiPsi.shape)
                for q in xrange(self.qclayers.EigenE.size):
                    xyPsiPsiEig[:,q] = self.qclayers.xyPsiPsi[:,q] + self.qclayers.EigenE[q]
                savetxt(fname.split('.')[0] + '_States' + '.csv', column_stack([self.qclayers.xPointsPost, xyPsiPsiEig]), delimiter=',')




#===============================================================================
# Edit Menu Items
#===============================================================================

    def bump_first_layer(self):
        self.qclayers.layerWidths = hstack([self.qclayers.layerWidths[-1], self.qclayers.layerWidths])
        self.qclayers.layerBarriers = hstack([self.qclayers.layerBarriers[-1], self.qclayers.layerBarriers])
        self.qclayers.layerARs = hstack([self.qclayers.layerARs[-1], self.qclayers.layerARs])
        self.qclayers.layerMaterials = hstack([self.qclayers.layerMaterials[-1], self.qclayers.layerMaterials])
        self.qclayers.layerDopings = hstack([self.qclayers.layerDopings[-1], self.qclayers.layerDopings])
        self.qclayers.layerDividers = hstack([self.qclayers.layerDividers[-1], self.qclayers.layerDividers])
    
        self.update_inputBoxes()
        self.layerTable_refresh()
        self.layerTable.setCurrentCell(1,0)
        self.layerTable.setFocus()
        
        self.dirty = True
        self.update_windowTitle()

    def set_temperature(self):
        nowTemp = ThePhysics.c.Temperature
        newTemp, buttonResponse = QInputDialog.getDouble(self, 'ErwinJr Input Dialog', 'Set Temperature', value=nowTemp, min=0)
        if buttonResponse:
            ThePhysics.c.set_temperature(newTemp)
            self.qclayers.Temperature = newTemp
            self.qclayers.populate_x()
            self.qclayers.populate_x_full()
            self.update_quantumCanvas()
            
    def copy_structure(self):
        clipboard = QApplication.clipboard()
        string = ''
        for layer in self.qclayers.layerWidths[1:]:
            string += '%g\n' % layer
        clipboard.setText(string)




#===============================================================================
# View Menu Items
#===============================================================================

    def view_VXBand(self):
        if self.plotVX:
            self.plotVX = False
        else:
            self.plotVX = True
        self.plotDirty = True
        self.update_quantumCanvas()
        
    def view_VLBand(self):
        if self.plotVL:
            self.plotVL = False
        else:
            self.plotVL = True
        self.plotDirty = True
        self.update_quantumCanvas()
        
    def view_LHBand(self):
        if self.plotLH:
            self.plotLH = False
        else:
            self.plotLH = True
        self.plotDirty = True
        self.update_quantumCanvas()
        
    def view_SOBand(self):
        if self.plotSO:
            self.plotSO = False
        else:
            self.plotSO = True
        self.plotDirty = True
        self.update_quantumCanvas()




#===============================================================================
# Help Menu Items
#===============================================================================

    def on_about(self):
        msg = """ ErwinJr 2.x Authors and Contributors
        
         * Kale J. Franz, PhD (Jet Propulsion Laboratory)
            kfranz@alumni.princeton.edu
            www.kalefranz.com
            
With contributions from:
         * Yamac Dikmelik (Johns Hopkins University)
         * Yu Song (Princeton University)
        """
        QMessageBox.about(self, "ErwinJr " + str(ejVersion), msg.strip())
        
    def on_licenses(self):
        copyright1 = """
#=======================================
# ErwinJr is a simulation program for quantum semiconductor lasers.
# Copyright (C) 2012 Kale J. Franz, PhD
#
# A portion of this code is Copyright (c) 2011, California Institute of 
# Technology ("Caltech"). U.S. Government sponsorship acknowledged.
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#=======================================
"""
        QMessageBox.about(self, "ErwinJr " + str(ejVersion), copyright1.strip())
        
    def on_tutorial(self):
        if os.name == "nt":
            os.startfile("tutorial.pdf")
        elif os.name == "posix":
            os.system("/usr/bin/xdg-open tutorial.pdf")  





def main():
    app = QApplication(sys.argv)
    app.setOrganizationName("JPL")
    app.setOrganizationDomain("erwinjr.org")
    app.setApplicationName("ErwinJr")
    qsettingsSystem = QSettings(QSettings.SystemScope,"JPL","ErwinJr")
    installDirectory = str(qsettingsSystem.value('installDirectory').toString())
    if installDirectory:
        os.chdir(installDirectory)
    
        
    
    # Create and display the splash screen
    splash_pix = QPixmap('images/erwinjr_splash.png')
    splash = QSplashScreen(splash_pix, Qt.WindowStaysOnTopHint)
    splash.setMask(splash_pix.mask())
    splash.show()
    app.processEvents()
    
    time.sleep(1)
    
    app.setWindowIcon(QIcon('images/EJpng48x48.png'))

    #set to cleanlooks on macs and linux
    if os.name == "posix":
        app.setStyle('cleanlooks')
    
    #this block handles a filename passed in by command line
    try:
        fileName = sys.argv[1]
        name, ext = os.path.splitext(fileName)
        assert ext == ".qcl"
        assert os.path.exists(fileName)
        fileName = os.path.abspath(fileName)
    except (IndexError, AssertionError):
        fileName = None

    form = MainWindow(fileName)
    form.show()
    splash.finish(form)
    
    # Import Psyco if available
    try:
        import psyco
        psyco.full()
    except ImportError:
        noPsycoBox = QMessageBox(QMessageBox.Question, 'EwrinJr '+str(majorVersion), "Psyco could not be loaded.\nExecution will be slowed.")
        noPsycoBox.exec_()
    
    qsettings = QSettings()
    if not qsettings.value('firstRun').toInt()[1]:
        if not installDirectory:
            qsettingsSystem.setValue("installDirectory", QVariant(os.getcwd()))
        firstRunBox = QMessageBox(QMessageBox.Question, 'EwrinJr '+str(majorVersion), "Welcome to ErwinJr!\n\
            Since this is your first time running the program, would you like to open an example file or a blank file?",
            parent=form)
        firstRunBox.addButton("Blank File", QMessageBox.NoRole)
        firstRunBox.addButton("Example File", QMessageBox.YesRole)
        ansr =firstRunBox.exec_()
        if ansr:
            form.fileOpen('examples/NPhoton PQLiu.qcl')
        else:
            form.fileNew()
        qsettings.setValue("firstRun", 1)
    
    
    app.exec_()

main()

