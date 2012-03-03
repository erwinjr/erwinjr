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
from numpy import *
from PyQt4.QtCore import *
from PyQt4.QtGui import *
import PyQt4.Qwt5 as Qwt


#===============================================================================
# Global Variables
#===============================================================================
e0 = 1.60217653e-19;  #electron charge
eps0 = 8.85e-12;
m0 = 9.10938188e-31;   #free electron mass (kg)
hbar = 6.6260693e-34/(2*pi); #Planck's constant (J s)
kb = 1.386505e-23 / e0; #eV/K


class Spy(QObject):
    
    def __init__(self, parent):
        QObject.__init__(self, parent)
        parent.setMouseTracking(True)
        parent.installEventFilter(self)

    # __init__()

    def eventFilter(self, _, event):
        if event.type() == QEvent.MouseMove:
            self.emit(SIGNAL("MouseMove"), event.pos())
        return False

    # eventFilter()

# class Spy

class MaskedData(Qwt.QwtArrayData):

    def __init__(self, x, y, mask):
        Qwt.QwtArrayData.__init__(self, x, y)
        self.__mask = asarray(mask, bool)
        # keep a copy of x and y for boundingRect()
        self.__x = asarray(x)
        self.__y = asarray(y)

    # __init__()

    def copy(self):
        return self

    # copy()

    def mask(self):
        return self.__mask

    # mask()

    def boundingRect(self):
        """Return the bounding rectangle of the data, accounting for the mask.
        """
        try:
            xmax = self.__x[self.__mask].max()
            xmin = self.__x[self.__mask].min()
            ymax = self.__y[self.__mask].max()
            ymin = self.__y[self.__mask].min()
    
            return QRectF(xmin, ymin, xmax-xmin, ymax-ymin)
        except ValueError:
            return QRectF()           

    # boundingRect()

# class MaskedData


class MaskedCurve(Qwt.QwtPlotCurve):

    def __init__(self, x, y, mask):
        Qwt.QwtPlotCurve.__init__(self)
        self.setData(MaskedData(x, y, mask))

    # __init__()

    def draw(self, painter, xMap, yMap, rect):
        # When the array indices contains the indices of all valid data points,
        # a chunks of valid data is indexed by
        # indices[first], indices[first+1], .., indices[last].
        # The first index of a chunk of valid data is calculated by:
        # 1. indices[i] - indices[i-1] > 1
        # 2. indices[0] is always OK
        # The last index of a chunk of valid data is calculated by:
        # 1. index[i] - index[i+1] < -1
        # 2. index[-1] is always OK
        indices = arange(self.data().size())[self.data().mask()]
        fs = array(indices)
        try:        
            fs[1:] -= indices[:-1]
            fs[0] = 2
            fs = indices[fs > 1]
            ls = array(indices)
            ls[:-1] -= indices[1:]
            ls[-1] = -2
            ls = indices[ls < -1]
            for first, last in zip(fs, ls):
                Qwt.QwtPlotCurve.drawFromTo(self, painter, xMap, yMap, first, last)
        except IndexError:
            pass
# class MaskedCurve

def matlab_range(rangeString):
    #rangeString = '1:1rangeString = self.currentStringBox.text():10 2:2:20'
    rangeString = rangeString.replace(',', ' ') #make comma and space the same
    matlabValues = []
    rangeList = str(rangeString).split(' ')
    for series in rangeList:
        item = series.split(':')
        if len(item) == 1:
            matlabValues.append(float(item[0]))
        elif len(item) == 2:
            matlabValues.extend(arange(float(item[0]),float(item[1])+1,1))
        elif len(item) == 3:
            matlabValues.extend(arange(float(item[0]),float(item[1])+float(item[2]),float(item[1])))
    return matlabValues
#matlab_range