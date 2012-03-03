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
from scipy import interpolate
import copy
from multiprocessing import Process, Queue

#import pylab as plt

import settings

import MaterialConstants
#create global variable c for material constants
c = MaterialConstants.MaterialConstants()

from ctypes import *
try:
    cFunctions=CDLL('cFunctions.so')
except WindowsError:
    cFunctions=CDLL('cFunctions.dll')



#===============================================================================
# Global Variables
#===============================================================================
e0 = 1.60217653e-19;  #electron charge
eps0 = 8.854187e-12;
m0 = 9.10938188e-31;   #free electron mass (kg)
h = 6.6260693e-34;
hbar = 6.6260693e-34/(2*pi); #Planck's constant (J s)
kb = 1.386505e-23 / e0; #eV/K
c0 = 299792458;

class Strata(object):
    def __init__(self):
        self.stratumMaterials = ['InP']
        self.stratumCompositions = array([0.])
        self.stratumThicknesses = array([0.])
        self.stratumDopings = array([0.])
        
        self.wavelength = 4.7
        self.operatingField = 0
        self.Lp = 1
        self.Np = 1
        self.aCore = 0
        self.nCore = 3
        self.nD = 0
        
        self.tauUpper = 0.0
        self.tauLower = 0.0
        self.tauUpperLower = 1.0e-3
        self.opticalDipole = 0.0
        self.FoM = 0.0
        self.transitionBroadening = 1.0e-5
        self.waveguideFacets = 'as-cleaved + as-cleaved'
        self.customFacet = 0.0
        self.waveguideLength = 3.0
        
        self.frontFacet = 0
        self.backFacet = 0
        
        self.beta = 3+0j
        
        self.xres = 0.01 #um
        self.stratumSelected = 0
        
        self.notDopableList = ['Air', 'Au', 'SiO2', 'SiNx']
        self.needsCompositionList = ['InGaAs','InAlAs']

        self.populate_rIndexes()
        
    def populate_rIndexes(self):
        wl = self.wavelength
        n_GaAs = sqrt(c.C1_GaAs+c.C2_GaAs*wl**2/(wl**2-c.C3_GaAs**2)+c.C4_GaAs*wl**2/(wl**2-c.C5_GaAs**2))
        n_InAs = sqrt(c.C1_InAs+c.C2_InAs*wl**2/(wl**2-c.C3_InAs**2)+c.C4_InAs*wl**2/(wl**2-c.C5_InAs**2))
        n_AlAs = sqrt(c.C1_AlAs+c.C2_AlAs*wl**2/(wl**2-c.C3_AlAs**2)+c.C4_AlAs*wl**2/(wl**2-c.C5_AlAs**2))
        n_InP = sqrt(c.C1_InP+c.C2_InP*wl**2/(wl**2-c.C3_InP**2)+c.C4_InP*wl**2/(wl**2-c.C5_InP**2))
        
        self.stratumRIndexes = zeros(self.stratumDopings.size, dtype=complex)
        for q, material in enumerate(self.stratumMaterials):
            if material == 'Active Core':
                self.stratumRIndexes[q] = self.nCore
            elif material == 'InP':
                nue = 1
                me0 = c.me0_InP
                a = 8.97E-5*wl**2/me0*self.stratumDopings[q]
                eps = n_InP**2 - a / (1+1j*5.305e-3*wl**2*nue)
                n_InPd = sqrt( 0.5 *(abs(eps) + eps.real))
                k_InPd = sqrt( 0.5 *(abs(eps) - eps.real))
                self.stratumRIndexes[q] = n_InPd + 1j*k_InPd
            elif material == 'GaAs':
                nue = 1
                me0 = c.me0_GaAs
                a = 8.97E-5*wl**2/me0*self.stratumDopings[q]
                eps = n_GaAs**2 - a / (1+1j*5.305e-3*wl**2*nue)
                n_GaAsd = sqrt( 0.5 *(abs(eps) + eps.real))
                k_GaAsd = sqrt( 0.5 *(abs(eps) - eps.real))
                self.stratumRIndexes[q] = n_GaAsd + 1j*k_GaAsd
            elif material == 'InGaAs':
                xFrac = self.stratumCompositions[q]
                n_InGaAs = xFrac*n_InAs + (1-xFrac)*n_GaAs
                nue = 1
                me0 = xFrac*c.me0_InAs + (1-xFrac)*c.me0_GaAs
                a = 8.97E-5*wl**2/me0*self.stratumDopings[q]
                eps = n_InGaAs**2 - a / (1+1j*5.305e-3*wl**2*nue)
                n_InGaAs = sqrt( 0.5 *(abs(eps) + eps.real))
                k_InGaAs = sqrt( 0.5 *(abs(eps) - eps.real))
                self.stratumRIndexes[q] = n_InGaAs + 1j*k_InGaAs
            elif material == 'InAlAs':
                xFrac = self.stratumCompositions[q]
                n_AlInAs = (1-xFrac)*n_AlAs + xFrac*n_InAs
                nue = 1
                me0 = (1-xFrac)*c.me0_AlAs + xFrac*c.me0_InAs
                a = 8.97E-5*wl**2/me0*self.stratumDopings[q]
                eps = n_AlInAs**2 - a / (1+1j*5.305e-3*wl**2*nue)
                n_AlInAs = sqrt( 0.5 *(abs(eps) + eps.real))
                k_AlInAs = sqrt( 0.5 *(abs(eps) - eps.real))
                self.stratumRIndexes[q] = n_AlInAs + 1j*k_AlInAs
            elif material == 'Au':
                C1=-0.1933; C2=0.3321; C3=0.0938;
                D1=-0.382; D2=6.8522; D3=-0.1289;
                n_Au = C1+wl*C2+wl*C3**2
                k_Au = D1+wl*D2+wl*D3**2
                self.stratumRIndexes[q] = n_Au+k_Au*1j
            elif material == 'SiNx':
                #from Jean Nguyen's Thesis
                C1 = 2.0019336; C2 = 0.15265213; C3 = 4.0495557;
                D0=-0.00282; D1=0.003029; D2=-0.0006982; D3=-0.0002839; D4=0.0001816; D5=-3.948e-005; D6=4.276e-006; D7=-2.314e-007; D8=4.982e-009;
                n_SiNx = C1 + C2/wl**2 + C3/wl**4
                k_SiNx = D0+D1*wl+D2*wl**2+D3*wl**3+D4*wl**4+D5*wl**5+D6*wl**6+D7*wl**7+D8*wl**8
                k_SiNx *= 100
                self.stratumRIndexes[q] = n_SiNx+k_SiNx*1j
            elif material == 'SiO2':
                #from Jean Nguyen's Thesis
                C1 = 1.41870; C2 = 0.12886725; C3 = 2.7573641e-5
                n_SiO2 = C1 + C2/wl**2 + C3/wl**4
                
                #this is a 4 peak Lorentzian fit to her data
                y0=-797.4627
                xc1=2.83043; w1=6.083822; A1=10881.9438
                xc2=8.95338; w2=1.38389113; A2=9167.662815
                xc3=12.3845492; w3=3.9792077; A3=12642.72911
                xc4=15.6387213; w4=0.6057751177; A4=3292.325272
                alpha = y0 + 2*A1/pi*w1/(4*(wl-xc1)**2+w1**2) + 2*A2/pi*w2/(4*(wl-xc2)**2+w2**2) \
                        + 2*A3/pi*w3/(4*(wl-xc3)**2+w3**2) + 2*A4/pi*w4/(4*(wl-xc4)**2+w4**2)
                k_SiO2 = alpha * wl*1e-4 / (4*pi)
                self.stratumRIndexes[q] = n_SiO2 + k_SiO2*1j
            elif material == 'Air':
                self.stratumRIndexes[q] = 1
                
    def get_nCore(self, data):
        wl = self.wavelength
        n_GaAs = sqrt(c.C1_GaAs+c.C2_GaAs*wl**2/(wl**2-c.C3_GaAs**2)+c.C4_GaAs*wl**2/(wl**2-c.C5_GaAs**2))
        n_InAs = sqrt(c.C1_InAs+c.C2_InAs*wl**2/(wl**2-c.C3_InAs**2)+c.C4_InAs*wl**2/(wl**2-c.C5_InAs**2))
        n_AlAs = sqrt(c.C1_AlAs+c.C2_AlAs*wl**2/(wl**2-c.C3_AlAs**2)+c.C4_AlAs*wl**2/(wl**2-c.C5_AlAs**2))
        
        n=zeros(8)
        n[0] = data.moleFrac1*n_InAs + (1-data.moleFrac1)*n_GaAs;
        n[1] = data.moleFrac2*n_InAs + (1-data.moleFrac2)*n_AlAs;
        n[2] = data.moleFrac3*n_InAs + (1-data.moleFrac3)*n_GaAs;
        n[3] = data.moleFrac4*n_InAs + (1-data.moleFrac4)*n_AlAs;
        n[4] = data.moleFrac5*n_InAs + (1-data.moleFrac5)*n_GaAs;
        n[5] = data.moleFrac6*n_InAs + (1-data.moleFrac6)*n_AlAs;
        n[6] = data.moleFrac7*n_InAs + (1-data.moleFrac7)*n_GaAs;
        n[7] = data.moleFrac8*n_InAs + (1-data.moleFrac8)*n_AlAs;
        nCore = sum(data.h*n)/sum(data.h)
        
        kCore = 1/(4*pi) * self.aCore * wl*1e-4
        
        return nCore+kCore*1j
        
    def populate_x(self):
        #use rounding to work with selected resolution
        self.stratumThicknesses /= self.xres
        self.stratumThicknesses  = self.stratumThicknesses.round()
        self.stratumThicknesses *= self.xres
        
        #convert to int to prevent machine rounding errors
        self.xPoints = arange(0.,int(self.stratumThicknesses.sum()/self.xres),1)
        self.xPoints *= self.xres
        
        stratumThicknessesCumSum = concatenate([[0.],self.stratumThicknesses.cumsum()])
        self.xn = zeros(self.xPoints.size, dtype=complex)
        self.xAC = zeros(self.xPoints.size) #binary designation for Active Core

        #extend layer data for all xpoints
        for q in xrange(0,self.stratumThicknesses.size):
            self.xn[stratumThicknessesCumSum[q]/self.xres:stratumThicknessesCumSum[q+1]/self.xres] = self.stratumRIndexes[q]
            if self.stratumMaterials[q] == 'Active Core':
                self.xAC[stratumThicknessesCumSum[q]/self.xres:stratumThicknessesCumSum[q+1]/self.xres] = 1
                
        # make array to show selected stratum in mainCanvas
        self.xStratumSelected = zeros(self.xPoints.shape)*NaN
        sS = self.stratumSelected
        cs = stratumThicknessesCumSum
        if sS != -1: #if not no row selected
            if sS == 0: #row for first layer is selected
                self.xStratumSelected[cs[sS]/self.xres:cs[sS+1]/self.xres] = self.xn.real[cs[sS]/self.xres:cs[sS+1]/self.xres]
            else: 
                self.xStratumSelected[cs[sS]/self.xres:cs[sS+1]/self.xres] = self.xn.real[cs[sS]/self.xres:cs[sS+1]/self.xres]

    def chi_find(self, beta):
        z0 = 0.003768
        k = 2*pi/self.wavelength
        
        alpha = sqrt(self.stratumRIndexes**2-beta**2)
        if alpha[0].imag < 0:
            alpha[0] = conj(alpha[0])
        if alpha[-1].imag < 0:
            alpha[-1] = conj(alpha[-1])
        gamma = z0*alpha/self.stratumRIndexes**2
        phi   = k*self.stratumThicknesses*alpha
        #zeta  = k*self.stratumThicknesses/z0
        
        Mj = []
        M = array([[1+0j,0],[0,1+0j]])
        for q in xrange(alpha.size):
            Mj.append(array([[cos(phi[q]), -1j/gamma[q]*sin(phi[q])],[-1j*gamma[q]*sin(phi[q]), cos(phi[q])]]))
        
        Mj.reverse()
        for mj in Mj:
            M = dot(mj,M)
            
        gammas = gamma[0]
        gammac = gamma[-1]
        
        chi = gammac*M[0,0] + gammac*gammas*M[0,1] + M[1,0] + gammas*M[1,1]
        return chi
        
    def zero_find(self, xVals, yVals):
        tck = interpolate.splrep(xVals.real,yVals.real,s=0)
        xnew = linspace(min(xVals.real),max(xVals.real),5e5)
        ynew = interpolate.splev(xnew,tck,der=0)

        #This routine looks for all of the zero crossings, and then picks each one out
        gtz = ynew > 0
        ltz = ynew < 0
        overlap1 = bitwise_and(gtz[0:-1],ltz[1:])
        overlap2 = bitwise_and(gtz[1:],ltz[0:-1])
        overlap  = bitwise_or(overlap1, overlap2)
        idxs = nonzero(overlap == True)[0]
        zeroCrossings = zeros(idxs.size)
        
        cFunctions.inv_quadratic_interp(xnew.ctypes.data_as(c_void_p), ynew.ctypes.data_as(c_void_p), 
                                                idxs.ctypes.data_as(c_void_p), int(zeroCrossings.size), 
                                                zeroCrossings.ctypes.data_as(c_void_p))
        return zeroCrossings
    
    def beta_find(self, betaInit = None):
        if True: #betaInit == None:
            betaMax  = max(self.stratumRIndexes.real)
            betaMin  = min(self.stratumRIndexes.real)

            betas = arange(betaMin.real+0.01,betaMax.real,0.01)

            if True: #do chi_find in c
                chiImag = zeros(len(betas),dtype=float)
                betasReal = betas.real
                betasImag = betas.imag
                stratumRIndexesReal = self.stratumRIndexes.real.copy()
                stratumRIndexesImag = self.stratumRIndexes.imag.copy()
                cFunctions.chiImag_array(c_double(self.wavelength), self.stratumThicknesses.ctypes.data_as(c_void_p), 
                                         stratumRIndexesReal.ctypes.data_as(c_void_p), stratumRIndexesImag.ctypes.data_as(c_void_p), 
                                         int(self.stratumRIndexes.size), betasReal.ctypes.data_as(c_void_p), betasImag.ctypes.data_as(c_void_p),
                                         int(betasReal.size), chiImag.ctypes.data_as(c_void_p))
                beta0s = self.zero_find(betas.real, chiImag)
            else:
                chi=zeros(betas.size, dtype=complex)
                for p, beta in enumerate(betas):
                    chi[p] = self.chi_find(beta)
                beta0s = self.zero_find(betas.real, chi.imag)
            beta = max(beta0s)+1j*min(self.stratumRIndexes.imag)
        else:
            beta = betaInit
        
        beta_find_precision = 1e-5
        if True:
            betaIn = beta
            stratumRIndexesReal = self.stratumRIndexes.real.copy()
            stratumRIndexesImag = self.stratumRIndexes.imag.copy()
            betaOut = array([0.0, 0.0])
            beta = cFunctions.beta_find(c_double(self.wavelength), self.stratumThicknesses.ctypes.data_as(c_void_p), 
                             stratumRIndexesReal.ctypes.data_as(c_void_p), stratumRIndexesImag.ctypes.data_as(c_void_p), 
                             int(self.stratumRIndexes.size), c_double(betaIn.real), c_double(betaIn.imag),
                             c_double(beta_find_precision), betaOut.ctypes.data_as(c_void_p))
            beta = betaOut[0] + 1j*betaOut[1]
        else:
            rInc = 0.0001; iInc = 1j*1e-6
            abschiNew=1
            while True:
                betas = [beta, beta+rInc, beta-rInc, beta+iInc, beta-iInc, 
                         beta+rInc+iInc, beta-rInc-iInc, beta+rInc-iInc,
                         beta-rInc+iInc]
                if True:
                    chi = zeros(len(betas),dtype=complex)
                    for p, betaIn in enumerate(betas):
                        chi[p] = self.chi_find(betaIn)
                else: #do chi_find in c
                    chi = zeros(len(betas),dtype=float)
                    abschi_find = cFunctions.abschi_find
                    abschi_find.restype = c_double
                    for p, betaIn in enumerate(betas):
                        stratumRIndexesReal = self.stratumRIndexes.real.copy()
                        stratumRIndexesImag = self.stratumRIndexes.imag.copy()
                        chi[p] = abschi_find(c_double(self.wavelength), self.stratumThicknesses.ctypes.data_as(c_void_p), 
                                             stratumRIndexesReal.ctypes.data_as(c_void_p), stratumRIndexesImag.ctypes.data_as(c_void_p), 
                                             int(self.stratumRIndexes.size), c_double(betaIn.real), c_double(betaIn.imag))
                abschiOld = abschiNew
                abschiNew = min(abs(chi))
                idx=argmin(abs(chi))
                beta=betas[idx]
                if abs(abschiOld -abschiNew)/abschiOld < beta_find_precision:
                    break
        return beta
        
    def mode_plot(self):
        n=copy.copy(self.stratumRIndexes)[::-1]
        thicknesses = copy.copy(self.stratumThicknesses)[::-1]
        xres = self.xres        
        
        z0 = 0.003768
        #z0 = 376.8
        k = 2*pi/self.wavelength
        M = array([[1+0j,0],[0,1+0j]])
        
        alpha = sqrt(n**2-self.beta**2)
        if alpha[0].imag < 0:
            alpha[0] = conj(alpha[0])
        if alpha[-1].imag < 0:
            alpha[-1] = conj(alpha[-1])
        gamma = z0*alpha/n**2
        phi   = k*thicknesses*alpha
        #zeta  = k*thicknesses/z0
        
        cs = stratumThicknessesCumSum = concatenate([[0.],thicknesses.cumsum()])
        xI = zeros(self.xPoints.size, dtype=complex)

        for q in xrange(thicknesses.size-1,-1,-1):
            xvec  = copy.copy(self.xPoints[cs[q]/xres:cs[q+1]/xres])[::-1]
            if len(xvec) == 0: #make sure xvec isn't empty
                continue
            xvec -= min(xvec)
            if q == 0:
                field = dot(M,array([1,gamma[-1]]))
                U = field[0]
                xI[cs[q]/xres:cs[q+1]/xres] = real(U)*exp(1j*k*alpha[q]*xvec) / n[q]**2
            elif q == self.stratumThicknesses.size-1:
                field = dot(M,array([1,gamma[-1]]))
                U = field[0]
                xI[cs[q]/xres:cs[q+1]/xres] = real(U)*exp(1j*k*alpha[q]*xvec) / n[q]**2
            else:
                field = dot(M,array([1,gamma[-1]]))
                U = field[0]; V = field[1]
                xI[cs[q]/xres:cs[q+1]/xres]  = U*cos(-k*alpha[q]*xvec) + 1j/gamma[q] * V*sin(-k*alpha[q]*xvec)
                xI[cs[q]/xres:cs[q+1]/xres] /= n[q]**2
                Mj = array([[cos(phi[q]), -1j/gamma[q]*sin(phi[q])],[-1j*gamma[q]*sin(phi[q]), cos(phi[q])]])
                M = dot(Mj,M)
        
        xI = abs(xI)**2 / max(abs(xI)**2)
        self.xI = xI[::-1]
        
        #calculate confinement factor
        self.confinementFactor = sum(self.xI * self.xAC) / sum(self.xI)
        
    def calculate_performance_parameters(self):
        #waveguide loss
        self.waveguideLoss = 4 * pi * self.beta.imag / (self.wavelength * 1e-6) * 1e-2
        
        #mirror loss
        self.mirrorLoss = -1 / (2 * self.waveguideLength * 0.1) * log(self.frontFacet * self.backFacet)
        
        #transition cross-section
        Eph = h * c0 / (self.wavelength * 1e-6)
        neff = self.beta.real
        z = self.opticalDipole * 1e-10
        deltaE = 0.1*Eph
        sigma0 = 4*pi*e0**2 / (h*c0*eps0*neff) * Eph/deltaE * z**2
        
        #gain
        tauEff = self.tauUpper * (1 - self.tauLower / self.tauUpperLower) * 1e-12
        Lp = self.Lp * 1e-10
        self.gain = sigma0 * tauEff / (e0 * Lp) #m/A
        self.gain *= 100 #cm/A
        
        #threshold current density
        self.Jth0 = (self.waveguideLoss + self.mirrorLoss) / (self.gain * self.confinementFactor) #A/cm^2
        self.Jth0 *= 1e-3 #kA/cm^2
        
        #threshold current
        self.Ith0 = self.Jth0*1e3 * (self.Np * self.Lp*1e-8) * self.waveguideLength*1e-1
        
        #operating voltage
        self.operatingVoltage = self.operatingField*1e3 * self.Lp*1e-8 * self.Np
        
        #voltage efficiency
        self.voltageEfficiency = 1.24/self.wavelength * self.Np / self.operatingVoltage
        
        #extraction efficiency
        self.extractionEfficiency = self.mirrorLoss / (self.mirrorLoss + self.waveguideLoss)
        
        #population inverstion efficiency
        tauEff = self.tauUpper * (1 - self.tauLower / self.tauUpperLower)
        self.inversionEfficiency = tauEff / (tauEff + self.tauLower)
        
        #modal efficiency
        xI = self.xI/max(self.xI)
        U = xI[nonzero(self.xAC)[0]]
        
        #interoplate over U_AC for each Np at the point xbar for Ubar
        #this is Faist's version
        numACs = self.stratumMaterials.count('Active Core')
        try:
            xVals = arange(self.xres,numACs*self.Np*self.Lp*1e-4,self.xres)
            assert xVals.size == U.size
        except AssertionError:
            try:
                xVals = arange(0,numACs*self.Np*self.Lp*1e-4,self.xres)
                assert xVals.size == U.size
            except AssertionError:
                xVals = arange(self.xres,numACs*self.Np*self.Lp*1e-4-self.xres,self.xres)
                assert xVals.size == U.size
        tck = interpolate.splrep(xVals,U,s=0)
        minx = 0.5*self.Lp*1e-4
        maxx = numACs*self.Np*self.Lp*1e-4-0.5*self.Lp*1e-4
        xbar = linspace(minx, maxx, numACs*self.Np)
        Ubar = interpolate.splev(xbar,tck,der=0)
        self.modalEfficiency = sum(Ubar)**2 / (numACs * self.Np * sum(Ubar**2))
        
        #Kale's version
        modalEfficiency = sum(Ubar) / self.Np #since Ubar taken from normalized xI
        #I guess we'll go with Faist's version since he probably knows better than I do.
       

class QCLayers(object):
    def __init__(self):
        self.layerWidths = array([1.,1.])
        self.layerBarriers = array([0,0])
        self.layerARs = array([0,0])
        self.layerMaterials = array([1,1])
        self.layerDopings = array([0.,0.])
        self.layerDividers = array([0,0])
        
        self.xres = 0.5
        self.EField = 0
        self.layerSelected = -1
        self.vertRes = 0.5
        self.repeats = 1
        
        self.description = ""
        self.solver = "SolverH"
        self.Temperature = 300
        self.TempFoM = 300
        self.diffLength = 0
        self.basisARInjector = True
        self.basisInjectorAR = True
        self.designByAngs = True
        self.designByML = False
        self.substrate = 'InP'
        
        self.moleFrac1 = 0.53
        self.moleFrac2 = 0.52
        self.moleFrac3 = 0.53
        self.moleFrac4 = 0.52
        self.moleFrac5 = 0.53
        self.moleFrac6 = 0.52
        self.moleFrac7 = 0.53
        self.moleFrac8 = 0.52
        
        self.update_alloys()
        self.update_strain()
        self.populate_x()

    #create xValues for mainCanvas plot ONLY; rest of xValues populated in populate_x_full    
#    def __getattr__(self, attr):
#        return self[attr]

    def populate_x(self):
        #use rounding to work with selected resolution
        self.layerWidths /= self.xres
        self.layerWidths = self.layerWidths.round()
        self.layerWidths*= self.xres
        
        #convert to int to prevent machine rounding errors
        self.xPoints = arange(0.,int(self.layerWidths.sum()/self.xres),1)
        self.xPoints *= self.xres
        
        layerWidthsCumSum = concatenate([[0.],self.layerWidths.cumsum()])
        self.xBarriers = zeros(self.xPoints.shape)
        self.xARs = zeros(self.xPoints.shape)
        self.xMaterials = zeros(self.xPoints.shape)
        self.xDopings = zeros(self.xPoints.shape)
        self.xLayerNums = zeros(self.xPoints.shape)

        #extend layer data for all xpoints
        for q in xrange(0,self.layerWidths.size):
            self.xBarriers[ layerWidthsCumSum[q]/self.xres:round(layerWidthsCumSum[q+1]/self.xres)] = self.layerBarriers[q]
            if self.layerARs[q] == 1:
                self.xARs[layerWidthsCumSum[q]/self.xres-1:layerWidthsCumSum[q+1]/self.xres+1] = 1
            self.xMaterials[layerWidthsCumSum[q]/self.xres:layerWidthsCumSum[q+1]/self.xres] = self.layerMaterials[q]
            self.xDopings[layerWidthsCumSum[q]/self.xres:layerWidthsCumSum[q+1]/self.xres] = self.layerDopings[q]
            self.xLayerNums[layerWidthsCumSum[q]/self.xres:layerWidthsCumSum[q+1]/self.xres] = q

            
        #plt.plot(self.xPoints, self.xBarriers,'o')

        #duplicate layer based on user input repeats
        repeats = self.repeats
        if self.repeats >= 2:
            self.xPoints = arange(0,self.layerWidths.sum()+self.layerWidths[1:].sum()*(self.repeats-1),self.xres)
            self.xBarriers = hstack([self.xBarriers, tile(self.xBarriers[layerWidthsCumSum[1]/self.xres:], self.repeats-1)])
            self.xARs = hstack([self.xARs, tile(self.xARs[layerWidthsCumSum[1]/self.xres:], self.repeats-1)])
            self.xMaterials = hstack([self.xMaterials, tile(self.xMaterials[layerWidthsCumSum[1]/self.xres:], self.repeats-1)])
            self.xDopings = hstack([self.xDopings, tile(self.xDopings[layerWidthsCumSum[1]/self.xres:], self.repeats-1)])
            self.xLayerNums = hstack([self.xLayerNums, tile(self.xLayerNums[layerWidthsCumSum[1]/self.xres:], self.repeats-1)])
        
        
        #this hack is needed because sometimes self.xPoints is one element too big
        self.xPoints = self.xPoints[0:self.xBarriers.size]
        
        self.update_strain()
        indx1 = nonzero(self.xMaterials == 1)[0]
        indx2 = nonzero(self.xMaterials == 2)[0]
        indx3 = nonzero(self.xMaterials == 3)[0]
        indx4 = nonzero(self.xMaterials == 4)[0]
        self.xVc  = zeros(self.xPoints.size)
        self.xVX  = zeros(self.xPoints.size)
        self.xVL  = zeros(self.xPoints.size)
        self.xVLH = zeros(self.xPoints.size)
        self.xVSO = zeros(self.xPoints.size)
        if indx1.size != 0:
            self.xVc[indx1]  = self.EcG[1]*self.xBarriers[indx1] + self.EcG[0]*(1-self.xBarriers[indx1]) - self.xPoints[indx1] * self.EField * 1e-5
            self.xVX[indx1]  = self.EcX[1]*self.xBarriers[indx1] + self.EcX[0]*(1-self.xBarriers[indx1]) - self.xPoints[indx1] * self.EField * 1e-5
            self.xVL[indx1]  = self.EcL[1]*self.xBarriers[indx1] + self.EcL[0]*(1-self.xBarriers[indx1]) - self.xPoints[indx1] * self.EField * 1e-5
            self.xVLH[indx1] = self.EvLH[1]*self.xBarriers[indx1] + self.EvLH[0]*(1-self.xBarriers[indx1]) - self.xPoints[indx1] * self.EField * 1e-5
            self.xVSO[indx1] = self.EvSO[1]*self.xBarriers[indx1] + self.EvSO[0]*(1-self.xBarriers[indx1]) - self.xPoints[indx1] * self.EField * 1e-5
        if indx2.size != 0:
            self.xVc[indx2]  = self.EcG[3]*self.xBarriers[indx2] + self.EcG[2]*(1-self.xBarriers[indx2]) - self.xPoints[indx2] * self.EField * 1e-5
            self.xVX[indx2]  = self.EcX[3]*self.xBarriers[indx2] + self.EcX[2]*(1-self.xBarriers[indx2]) - self.xPoints[indx2] * self.EField * 1e-5
            self.xVL[indx2]  = self.EcL[3]*self.xBarriers[indx2] + self.EcL[2]*(1-self.xBarriers[indx2]) - self.xPoints[indx2] * self.EField * 1e-5
            self.xVLH[indx2] = self.EvLH[3]*self.xBarriers[indx2] + self.EvLH[2]*(1-self.xBarriers[indx2]) - self.xPoints[indx2] * self.EField * 1e-5
            self.xVSO[indx2] = self.EvSO[3]*self.xBarriers[indx2] + self.EvSO[2]*(1-self.xBarriers[indx2]) - self.xPoints[indx2] * self.EField * 1e-5
        if indx3.size != 0:
            self.xVc[indx3]  = self.EcG[5]*self.xBarriers[indx3] + self.EcG[4]*(1-self.xBarriers[indx3]) - self.xPoints[indx3] * self.EField * 1e-5
            self.xVX[indx3]  = self.EcX[5]*self.xBarriers[indx3] + self.EcX[4]*(1-self.xBarriers[indx3]) - self.xPoints[indx3] * self.EField * 1e-5
            self.xVL[indx3]  = self.EcL[5]*self.xBarriers[indx3] + self.EcL[4]*(1-self.xBarriers[indx3]) - self.xPoints[indx3] * self.EField * 1e-5
            self.xVLH[indx3] = self.EvLH[5]*self.xBarriers[indx3] + self.EvLH[4]*(1-self.xBarriers[indx3]) - self.xPoints[indx3] * self.EField * 1e-5
            self.xVSO[indx3] = self.EvSO[5]*self.xBarriers[indx3] + self.EvSO[4]*(1-self.xBarriers[indx3]) - self.xPoints[indx3] * self.EField * 1e-5
        if indx4.size != 0:
            self.xVc[indx4]  = self.EcG[7]*self.xBarriers[indx4] + self.EcG[6]*(1-self.xBarriers[indx4]) - self.xPoints[indx4] * self.EField * 1e-5
            self.xVX[indx4]  = self.EcX[7]*self.xBarriers[indx4] + self.EcX[6]*(1-self.xBarriers[indx4]) - self.xPoints[indx4] * self.EField * 1e-5
            self.xVL[indx4]  = self.EcL[7]*self.xBarriers[indx4] + self.EcL[6]*(1-self.xBarriers[indx4]) - self.xPoints[indx4] * self.EField * 1e-5
            self.xVLH[indx4] = self.EvLH[7]*self.xBarriers[indx4] + self.EvLH[6]*(1-self.xBarriers[indx4]) - self.xPoints[indx4] * self.EField * 1e-5
            self.xVSO[indx4] = self.EvSO[7]*self.xBarriers[indx4] + self.EvSO[6]*(1-self.xBarriers[indx4]) - self.xPoints[indx4] * self.EField * 1e-5
            
        #self.xVc = self.xBarriers * 0.520 - self.xPoints * self.EField * 1e-5
        
        # make array to show selected layer in mainCanvas
        try:
            self.xLayerSelected = zeros(self.xPoints.shape)*NaN
            layerSelected = self.layerSelected
            if layerSelected != -1:
                if layerSelected == 0: #row for first layer is selected
                    for repeat in xrange(1,repeats+1):
                        base = layerWidthsCumSum[-1]/self.xres * (repeat-1)
                        self.xLayerSelected[base+layerWidthsCumSum[layerSelected]/self.xres:base+layerWidthsCumSum[layerSelected+1]/self.xres+1] = self.xVc[base+layerWidthsCumSum[layerSelected]/self.xres:base+layerWidthsCumSum[layerSelected+1]/self.xres+1]
                elif self.layerSelected == self.layerWidths.size: #last (blank) layer row is selected
                    pass
                else: 
                    for repeat in xrange(1,repeats+1):
                        base = sum(self.layerWidths[1:])/self.xres*(repeat-1)
                        self.xLayerSelected[base+layerWidthsCumSum[layerSelected]/self.xres-1:base+layerWidthsCumSum[layerSelected+1]/self.xres+1] = self.xVc[base+layerWidthsCumSum[layerSelected]/self.xres-1:base+layerWidthsCumSum[layerSelected+1]/self.xres+1]
        except IndexError:
            #index error happens in SolveBasis when the selected layer is greater than the number of layers in the solve swath
            # however, xLayerSelected is not used for the SolveBasis function
            pass

        self.xARs[nonzero(self.xARs==0)[0]] = NaN
        self.xARs *= self.xVc

    def populate_x_full(self):
        indx1 = nonzero(self.xMaterials == 1)[0]
        indx2 = nonzero(self.xMaterials == 2)[0]
        indx3 = nonzero(self.xMaterials == 3)[0]
        indx4 = nonzero(self.xMaterials == 4)[0]
        self.xEg = zeros(self.xPoints.size)
        self.xMc = zeros(self.xPoints.size)
        self.xESO = zeros(self.xPoints.size)
        self.xEp = zeros(self.xPoints.size)
        self.xF = zeros(self.xPoints.size)
        if indx1.size != 0:
            self.xEg[indx1]  = self.EgLH[1]*self.xBarriers[indx1] + self.EgLH[0]*(1-self.xBarriers[indx1])
            self.xMc[indx1]  = self.me[1]*self.xBarriers[indx1] + self.me[0]*(1-self.xBarriers[indx1])
            self.xESO[indx1] = self.ESO[1]*self.xBarriers[indx1] + self.ESO[0]*(1-self.xBarriers[indx1])
            self.xEp[indx1]  = self.Ep[1]*self.xBarriers[indx1] + self.Ep[0]*(1-self.xBarriers[indx1])
            self.xF[indx1]   = self.F[1]*self.xBarriers[indx1] + self.F[0]*(1-self.xBarriers[indx1])
        if indx2.size != 0:
            self.xEg[indx2]  = self.EgLH[3]*self.xBarriers[indx2] + self.EgLH[2]*(1-self.xBarriers[indx2])
            self.xMc[indx2]  = self.me[3]*self.xBarriers[indx2] + self.me[2]*(1-self.xBarriers[indx2])
            self.xESO[indx2] = self.ESO[3]*self.xBarriers[indx2] + self.ESO[2]*(1-self.xBarriers[indx2])
            self.xEp[indx2]  = self.Ep[3]*self.xBarriers[indx2] + self.Ep[2]*(1-self.xBarriers[indx2])
            self.xF[indx2]   = self.F[3]*self.xBarriers[indx2] + self.F[2]*(1-self.xBarriers[indx2])                    
        if indx3.size != 0:
            self.xEg[indx3]  = self.EgLH[5]*self.xBarriers[indx3] + self.EgLH[4]*(1-self.xBarriers[indx3])
            self.xMc[indx3]  = self.me[5]*self.xBarriers[indx3] + self.me[4]*(1-self.xBarriers[indx3])
            self.xESO[indx3] = self.ESO[5]*self.xBarriers[indx3] + self.ESO[4]*(1-self.xBarriers[indx3])
            self.xEp[indx3]  = self.Ep[5]*self.xBarriers[indx3] + self.Ep[4]*(1-self.xBarriers[indx3])
            self.xF[indx3]   = self.F[5]*self.xBarriers[indx3] + self.F[4]*(1-self.xBarriers[indx3])
        if indx4.size != 0:
            self.xEg[indx4]  = self.EgLH[7]*self.xBarriers[indx4] + self.EgLH[6]*(1-self.xBarriers[indx4])
            self.xMc[indx4]  = self.me[7]*self.xBarriers[indx4] + self.me[6]*(1-self.xBarriers[indx4])
            self.xESO[indx4] = self.ESO[7]*self.xBarriers[indx4] + self.ESO[6]*(1-self.xBarriers[indx4])
            self.xEp[indx4]  = self.Ep[7]*self.xBarriers[indx4] + self.Ep[6]*(1-self.xBarriers[indx4])
            self.xF[indx4]   = self.F[7]*self.xBarriers[indx4] + self.F[6]*(1-self.xBarriers[indx4])              
        self.xMc = self.xMc

    def update_alloys(self):  # c is a Material_Constant class instance
        if self.substrate == 'InP':
            self.numMaterials = 8
            variables = ['EgG', 'EgL', 'EgX', 'VBO', 'DSO', 'me0', 'acG', 'acL', 'acX', 
                     'Ep', 'F', 'XiX', 'b', 'av', 'alG', 'beG', 'alL', 'beL', 'alX', 
                     'beX', 'epss', 'epsInf', 'hwLO', 'a0', 'c11', 'c12', 'h']
            for item in variables:
                setattr(self, item, zeros(self.numMaterials))
                
            self.EgG[0] = self.moleFrac1*c.EgG_InAs + (1-self.moleFrac1)*c.EgG_GaAs - self.moleFrac1*(1-self.moleFrac1)*c.EgG_InGaAs;
            self.EgL[0] = self.moleFrac1*c.EgL_InAs + (1-self.moleFrac1)*c.EgL_GaAs - self.moleFrac1*(1-self.moleFrac1)*c.EgL_InGaAs;
            self.EgX[0] = self.moleFrac1*c.EgX_InAs + (1-self.moleFrac1)*c.EgX_GaAs - self.moleFrac1*(1-self.moleFrac1)*c.EgX_InGaAs;
            self.VBO[0] = self.moleFrac1*c.VBO_InAs + (1-self.moleFrac1)*c.VBO_GaAs - self.moleFrac1*(1-self.moleFrac1)*c.VBO_InGaAs;
            self.DSO[0] = self.moleFrac1*c.DSO_InAs + (1-self.moleFrac1)*c.DSO_GaAs - self.moleFrac1*(1-self.moleFrac1)*c.DSO_InGaAs;
            self.me0[0] = self.moleFrac1*c.me0_InAs + (1-self.moleFrac1)*c.me0_GaAs - self.moleFrac1*(1-self.moleFrac1)*c.me0_InGaAs;
            self.acG[0] = self.moleFrac1*c.acG_InAs + (1-self.moleFrac1)*c.acG_GaAs - self.moleFrac1*(1-self.moleFrac1)*c.acG_InGaAs;
            self.acL[0] = self.moleFrac1*c.acL_InAs + (1-self.moleFrac1)*c.acL_GaAs - self.moleFrac1*(1-self.moleFrac1)*c.acL_InGaAs;
            self.acX[0] = self.moleFrac1*c.acX_InAs + (1-self.moleFrac1)*c.acX_GaAs - self.moleFrac1*(1-self.moleFrac1)*c.acX_InGaAs;
            self.Ep[0]  = self.moleFrac1*c.Ep_InAs  + (1-self.moleFrac1)*c.Ep_GaAs  - self.moleFrac1*(1-self.moleFrac1)*c.Ep_InGaAs;
            self.F[0]   = self.moleFrac1*c.F_InAs   + (1-self.moleFrac1)*c.F_GaAs   - self.moleFrac1*(1-self.moleFrac1)*c.F_InGaAs;
            self.XiX[0] = self.moleFrac1*c.XiX_InAs + (1-self.moleFrac1)*c.XiX_GaAs;
            self.b[0]   = self.moleFrac1*c.b_InAs   + (1-self.moleFrac1)*c.b_GaAs;
            self.av[0]  = self.moleFrac1*c.av_InAs  + (1-self.moleFrac1)*c.av_GaAs;
            self.alG[0] = self.moleFrac1*c.alG_InAs + (1-self.moleFrac1)*c.alG_GaAs;
            self.beG[0] = self.moleFrac1*c.beG_InAs + (1-self.moleFrac1)*c.beG_GaAs;
            self.alL[0] = self.moleFrac1*c.alL_InAs + (1-self.moleFrac1)*c.alL_GaAs;
            self.beL[0] = self.moleFrac1*c.beL_InAs + (1-self.moleFrac1)*c.beL_GaAs;
            self.alX[0] = self.moleFrac1*c.alX_InAs + (1-self.moleFrac1)*c.alX_GaAs;
            self.beX[0] = self.moleFrac1*c.beX_InAs + (1-self.moleFrac1)*c.beX_GaAs;
            self.epss[0]   = self.moleFrac1*c.epss_InAs   + (1-self.moleFrac1)*c.epss_GaAs;
            self.epsInf[0] = self.moleFrac1*c.epsInf_InAs + (1-self.moleFrac1)*c.epsInf_GaAs;
            self.hwLO[0]   = self.moleFrac1*c.hwLO_InAs   + (1-self.moleFrac1)*c.hwLO_GaAs;
            #self.Ped[0]   = 3*self.moleFrac1*(1-self.moleFrac1)*(-c.acG_InAs+c.acG_GaAs)*(c.alc_InAs-c.alc_GaAs)/c.alc_InP; #Ped is the third term from Walle eqn (21)
            self.a0[0]     = self.moleFrac1*c.alc_InAs    + (1-self.moleFrac1)*c.alc_GaAs;
            self.c11[0]    = self.moleFrac1*c.c11_InAs    + (1-self.moleFrac1)*c.c11_GaAs;
            self.c12[0]    = self.moleFrac1*c.c12_InAs    + (1-self.moleFrac1)*c.c12_GaAs;
            
            self.EgG[1] = self.moleFrac2*c.EgG_InAs + (1-self.moleFrac2)*c.EgG_AlAs - self.moleFrac2*(1-self.moleFrac2)*c.EgG_AlInAs;
            self.EgL[1] = self.moleFrac2*c.EgL_InAs + (1-self.moleFrac2)*c.EgL_AlAs - self.moleFrac2*(1-self.moleFrac2)*c.EgL_AlInAs;
            self.EgX[1] = self.moleFrac2*c.EgX_InAs + (1-self.moleFrac2)*c.EgX_AlAs - self.moleFrac2*(1-self.moleFrac2)*c.EgX_AlInAs;
            self.VBO[1] = self.moleFrac2*c.VBO_InAs + (1-self.moleFrac2)*c.VBO_AlAs - self.moleFrac2*(1-self.moleFrac2)*c.VBO_AlInAs;
            self.DSO[1] = self.moleFrac2*c.DSO_InAs + (1-self.moleFrac2)*c.DSO_AlAs - self.moleFrac2*(1-self.moleFrac2)*c.DSO_AlInAs;
            self.me0[1] = self.moleFrac2*c.me0_InAs + (1-self.moleFrac2)*c.me0_AlAs - self.moleFrac2*(1-self.moleFrac2)*c.me0_AlInAs;
            self.acG[1] = self.moleFrac2*c.acG_InAs + (1-self.moleFrac2)*c.acG_AlAs - self.moleFrac2*(1-self.moleFrac2)*c.acG_AlInAs;
            self.acL[1] = self.moleFrac2*c.acL_InAs + (1-self.moleFrac2)*c.acL_AlAs - self.moleFrac2*(1-self.moleFrac2)*c.acL_AlInAs;
            self.acX[1] = self.moleFrac2*c.acX_InAs + (1-self.moleFrac2)*c.acX_AlAs - self.moleFrac2*(1-self.moleFrac2)*c.acX_AlInAs;
            self.Ep[1]  = self.moleFrac2*c.Ep_InAs  + (1-self.moleFrac2)*c.Ep_AlAs  - self.moleFrac2*(1-self.moleFrac2)*c.Ep_AlInAs;
            self.F[1]   = self.moleFrac2*c.F_InAs   + (1-self.moleFrac2)*c.F_AlAs   - self.moleFrac2*(1-self.moleFrac2)*c.F_AlInAs;
            self.XiX[1] = self.moleFrac2*c.XiX_InAs + (1-self.moleFrac2)*c.XiX_AlAs;
            self.b[1]   = self.moleFrac2*c.b_InAs   + (1-self.moleFrac2)*c.b_AlAs;
            self.av[1]  = self.moleFrac2*c.av_InAs  + (1-self.moleFrac2)*c.av_AlAs;
            self.alG[1] = self.moleFrac2*c.alG_InAs + (1-self.moleFrac2)*c.alG_AlAs;
            self.beG[1] = self.moleFrac2*c.beG_InAs + (1-self.moleFrac2)*c.beG_AlAs;
            self.alL[1] = self.moleFrac2*c.alL_InAs + (1-self.moleFrac2)*c.alL_AlAs;
            self.beL[1] = self.moleFrac2*c.beL_InAs + (1-self.moleFrac2)*c.beL_AlAs;
            self.alX[1] = self.moleFrac2*c.alX_InAs + (1-self.moleFrac2)*c.alX_AlAs;
            self.beX[1] = self.moleFrac2*c.beX_InAs + (1-self.moleFrac2)*c.beX_AlAs;
            self.epss[1]   = self.moleFrac2*c.epss_InAs   + (1-self.moleFrac2)*c.epss_AlAs;
            self.epsInf[1] = self.moleFrac2*c.epsInf_InAs + (1-self.moleFrac2)*c.epsInf_AlAs;
            self.hwLO[1]   = self.moleFrac2*c.hwLO_InAs   + (1-self.moleFrac2)*c.hwLO_AlAs;
            self.a0[1]     = self.moleFrac2*c.alc_InAs    + (1-self.moleFrac2)*c.alc_AlAs;
            self.c11[1]    = self.moleFrac2*c.c11_InAs    + (1-self.moleFrac2)*c.c11_AlAs;
            self.c12[1]    = self.moleFrac2*c.c12_InAs    + (1-self.moleFrac2)*c.c12_AlAs;

            self.EgG[2] = self.moleFrac3*c.EgG_InAs + (1-self.moleFrac3)*c.EgG_GaAs - self.moleFrac3*(1-self.moleFrac3)*c.EgG_InGaAs;
            self.EgL[2] = self.moleFrac3*c.EgL_InAs + (1-self.moleFrac3)*c.EgL_GaAs - self.moleFrac3*(1-self.moleFrac3)*c.EgL_InGaAs;
            self.EgX[2] = self.moleFrac3*c.EgX_InAs + (1-self.moleFrac3)*c.EgX_GaAs - self.moleFrac3*(1-self.moleFrac3)*c.EgX_InGaAs;
            self.VBO[2] = self.moleFrac3*c.VBO_InAs + (1-self.moleFrac3)*c.VBO_GaAs - self.moleFrac3*(1-self.moleFrac3)*c.VBO_InGaAs;
            self.DSO[2] = self.moleFrac3*c.DSO_InAs + (1-self.moleFrac3)*c.DSO_GaAs - self.moleFrac3*(1-self.moleFrac3)*c.DSO_InGaAs;
            self.me0[2] = self.moleFrac3*c.me0_InAs + (1-self.moleFrac3)*c.me0_GaAs - self.moleFrac3*(1-self.moleFrac3)*c.me0_InGaAs;
            self.acG[2] = self.moleFrac3*c.acG_InAs + (1-self.moleFrac3)*c.acG_GaAs - self.moleFrac3*(1-self.moleFrac3)*c.acG_InGaAs;
            self.acL[2] = self.moleFrac3*c.acL_InAs + (1-self.moleFrac3)*c.acL_GaAs - self.moleFrac3*(1-self.moleFrac3)*c.acL_InGaAs;
            self.acX[2] = self.moleFrac3*c.acX_InAs + (1-self.moleFrac3)*c.acX_GaAs - self.moleFrac3*(1-self.moleFrac3)*c.acX_InGaAs;
            self.Ep[2]  = self.moleFrac3*c.Ep_InAs  + (1-self.moleFrac3)*c.Ep_GaAs  - self.moleFrac3*(1-self.moleFrac3)*c.Ep_InGaAs;
            self.F[2]   = self.moleFrac3*c.F_InAs   + (1-self.moleFrac3)*c.F_GaAs   - self.moleFrac3*(1-self.moleFrac3)*c.F_InGaAs;
            self.XiX[2] = self.moleFrac3*c.XiX_InAs + (1-self.moleFrac3)*c.XiX_GaAs;
            self.b[2]   = self.moleFrac3*c.b_InAs   + (1-self.moleFrac3)*c.b_GaAs;
            self.av[2]  = self.moleFrac3*c.av_InAs  + (1-self.moleFrac3)*c.av_GaAs;
            self.alG[2] = self.moleFrac3*c.alG_InAs + (1-self.moleFrac3)*c.alG_GaAs;
            self.beG[2] = self.moleFrac3*c.beG_InAs + (1-self.moleFrac3)*c.beG_GaAs;
            self.alL[2] = self.moleFrac3*c.alL_InAs + (1-self.moleFrac3)*c.alL_GaAs;
            self.beL[2] = self.moleFrac3*c.beL_InAs + (1-self.moleFrac3)*c.beL_GaAs;
            self.alX[2] = self.moleFrac3*c.alX_InAs + (1-self.moleFrac3)*c.alX_GaAs;
            self.beX[2] = self.moleFrac3*c.beX_InAs + (1-self.moleFrac3)*c.beX_GaAs;
            self.epss[2]   = self.moleFrac3*c.epss_InAs   + (1-self.moleFrac3)*c.epss_GaAs;
            self.epsInf[2] = self.moleFrac3*c.epsInf_InAs + (1-self.moleFrac3)*c.epsInf_GaAs;
            self.hwLO[2]   = self.moleFrac3*c.hwLO_InAs   + (1-self.moleFrac3)*c.hwLO_GaAs;
            self.a0[2]     = self.moleFrac3*c.alc_InAs    + (1-self.moleFrac3)*c.alc_GaAs;
            self.c11[2]    = self.moleFrac3*c.c11_InAs    + (1-self.moleFrac3)*c.c11_GaAs;
            self.c12[2]    = self.moleFrac3*c.c12_InAs    + (1-self.moleFrac3)*c.c12_GaAs;
           
            self.EgG[3] = self.moleFrac4*c.EgG_InAs + (1-self.moleFrac4)*c.EgG_AlAs - self.moleFrac4*(1-self.moleFrac4)*c.EgG_AlInAs;
            self.EgL[3] = self.moleFrac4*c.EgL_InAs + (1-self.moleFrac4)*c.EgL_AlAs - self.moleFrac4*(1-self.moleFrac4)*c.EgL_AlInAs;
            self.EgX[3] = self.moleFrac4*c.EgX_InAs + (1-self.moleFrac4)*c.EgX_AlAs - self.moleFrac4*(1-self.moleFrac4)*c.EgX_AlInAs;
            self.VBO[3] = self.moleFrac4*c.VBO_InAs + (1-self.moleFrac4)*c.VBO_AlAs - self.moleFrac4*(1-self.moleFrac4)*c.VBO_AlInAs;
            self.DSO[3] = self.moleFrac4*c.DSO_InAs + (1-self.moleFrac4)*c.DSO_AlAs - self.moleFrac4*(1-self.moleFrac4)*c.DSO_AlInAs;
            self.me0[3] = self.moleFrac4*c.me0_InAs + (1-self.moleFrac4)*c.me0_AlAs - self.moleFrac4*(1-self.moleFrac4)*c.me0_AlInAs;
            self.acG[3] = self.moleFrac4*c.acG_InAs + (1-self.moleFrac4)*c.acG_AlAs - self.moleFrac4*(1-self.moleFrac4)*c.acG_AlInAs;
            self.acL[3] = self.moleFrac4*c.acL_InAs + (1-self.moleFrac4)*c.acL_AlAs - self.moleFrac4*(1-self.moleFrac4)*c.acL_AlInAs;
            self.acX[3] = self.moleFrac4*c.acX_InAs + (1-self.moleFrac4)*c.acX_AlAs - self.moleFrac4*(1-self.moleFrac4)*c.acX_AlInAs;
            self.Ep[3]  = self.moleFrac4*c.Ep_InAs  + (1-self.moleFrac4)*c.Ep_AlAs  - self.moleFrac4*(1-self.moleFrac4)*c.Ep_AlInAs;
            self.F[3]   = self.moleFrac4*c.F_InAs   + (1-self.moleFrac4)*c.F_AlAs   - self.moleFrac4*(1-self.moleFrac4)*c.F_AlInAs;
            self.XiX[3] = self.moleFrac4*c.XiX_InAs + (1-self.moleFrac4)*c.XiX_AlAs;
            self.b[3]   = self.moleFrac4*c.b_InAs   + (1-self.moleFrac4)*c.b_AlAs;
            self.av[3]  = self.moleFrac4*c.av_InAs  + (1-self.moleFrac4)*c.av_AlAs;
            self.alG[3] = self.moleFrac4*c.alG_InAs + (1-self.moleFrac4)*c.alG_AlAs;
            self.beG[3] = self.moleFrac4*c.beG_InAs + (1-self.moleFrac4)*c.beG_AlAs;
            self.alL[3] = self.moleFrac4*c.alL_InAs + (1-self.moleFrac4)*c.alL_AlAs;
            self.beL[3] = self.moleFrac4*c.beL_InAs + (1-self.moleFrac4)*c.beL_AlAs;
            self.alX[3] = self.moleFrac4*c.alX_InAs + (1-self.moleFrac4)*c.alX_AlAs;
            self.beX[3] = self.moleFrac4*c.beX_InAs + (1-self.moleFrac4)*c.beX_AlAs;
            self.epss[3]   = self.moleFrac4*c.epss_InAs   + (1-self.moleFrac4)*c.epss_AlAs;
            self.epsInf[3] = self.moleFrac4*c.epsInf_InAs + (1-self.moleFrac4)*c.epsInf_AlAs;
            self.hwLO[3]   = self.moleFrac4*c.hwLO_InAs   + (1-self.moleFrac4)*c.hwLO_AlAs;
            self.a0[3]     = self.moleFrac4*c.alc_InAs    + (1-self.moleFrac4)*c.alc_AlAs;
            self.c11[3]    = self.moleFrac4*c.c11_InAs    + (1-self.moleFrac4)*c.c11_AlAs;
            self.c12[3]    = self.moleFrac4*c.c12_InAs    + (1-self.moleFrac4)*c.c12_AlAs;
           
            self.EgG[4] = self.moleFrac5*c.EgG_InAs + (1-self.moleFrac5)*c.EgG_GaAs - self.moleFrac5*(1-self.moleFrac5)*c.EgG_InGaAs;
            self.EgL[4] = self.moleFrac5*c.EgL_InAs + (1-self.moleFrac5)*c.EgL_GaAs - self.moleFrac5*(1-self.moleFrac5)*c.EgL_InGaAs;
            self.EgX[4] = self.moleFrac5*c.EgX_InAs + (1-self.moleFrac5)*c.EgX_GaAs - self.moleFrac5*(1-self.moleFrac5)*c.EgX_InGaAs;
            self.VBO[4] = self.moleFrac5*c.VBO_InAs + (1-self.moleFrac5)*c.VBO_GaAs - self.moleFrac5*(1-self.moleFrac5)*c.VBO_InGaAs;
            self.DSO[4] = self.moleFrac5*c.DSO_InAs + (1-self.moleFrac5)*c.DSO_GaAs - self.moleFrac5*(1-self.moleFrac5)*c.DSO_InGaAs;
            self.me0[4] = self.moleFrac5*c.me0_InAs + (1-self.moleFrac5)*c.me0_GaAs - self.moleFrac5*(1-self.moleFrac5)*c.me0_InGaAs;
            self.acG[4] = self.moleFrac5*c.acG_InAs + (1-self.moleFrac5)*c.acG_GaAs - self.moleFrac5*(1-self.moleFrac5)*c.acG_InGaAs;
            self.acL[4] = self.moleFrac5*c.acL_InAs + (1-self.moleFrac5)*c.acL_GaAs - self.moleFrac5*(1-self.moleFrac5)*c.acL_InGaAs;
            self.acX[4] = self.moleFrac5*c.acX_InAs + (1-self.moleFrac5)*c.acX_GaAs - self.moleFrac5*(1-self.moleFrac5)*c.acX_InGaAs;
            self.Ep[4]  = self.moleFrac5*c.Ep_InAs  + (1-self.moleFrac5)*c.Ep_GaAs  - self.moleFrac5*(1-self.moleFrac5)*c.Ep_InGaAs;
            self.F[4]   = self.moleFrac5*c.F_InAs   + (1-self.moleFrac5)*c.F_GaAs   - self.moleFrac5*(1-self.moleFrac5)*c.F_InGaAs;
            self.XiX[4] = self.moleFrac5*c.XiX_InAs + (1-self.moleFrac5)*c.XiX_GaAs;
            self.b[4]   = self.moleFrac5*c.b_InAs   + (1-self.moleFrac5)*c.b_GaAs;
            self.av[4]  = self.moleFrac5*c.av_InAs  + (1-self.moleFrac5)*c.av_GaAs;
            self.alG[4] = self.moleFrac5*c.alG_InAs + (1-self.moleFrac5)*c.alG_GaAs;
            self.beG[4] = self.moleFrac5*c.beG_InAs + (1-self.moleFrac5)*c.beG_GaAs;
            self.alL[4] = self.moleFrac5*c.alL_InAs + (1-self.moleFrac5)*c.alL_GaAs;
            self.beL[4] = self.moleFrac5*c.beL_InAs + (1-self.moleFrac5)*c.beL_GaAs;
            self.alX[4] = self.moleFrac5*c.alX_InAs + (1-self.moleFrac5)*c.alX_GaAs;
            self.beX[4] = self.moleFrac5*c.beX_InAs + (1-self.moleFrac5)*c.beX_GaAs;
            self.epss[4]   = self.moleFrac5*c.epss_InAs   + (1-self.moleFrac5)*c.epss_GaAs;
            self.epsInf[4] = self.moleFrac5*c.epsInf_InAs + (1-self.moleFrac5)*c.epsInf_GaAs;
            self.hwLO[4]   = self.moleFrac5*c.hwLO_InAs   + (1-self.moleFrac5)*c.hwLO_GaAs;
            #self.Ped[0]   = 3*self.moleFrac1*(1-self.moleFrac1)*(-c.acG_InAs+c.acG_GaAs)*(c.alc_InAs-c.alc_GaAs)/c.alc_InP; #Ped is the third term from Walle eqn (21)
            self.a0[4]     = self.moleFrac5*c.alc_InAs    + (1-self.moleFrac5)*c.alc_GaAs;
            self.c11[4]    = self.moleFrac5*c.c11_InAs    + (1-self.moleFrac5)*c.c11_GaAs;
            self.c12[4]    = self.moleFrac5*c.c12_InAs    + (1-self.moleFrac5)*c.c12_GaAs;
            
            self.EgG[5] = self.moleFrac6*c.EgG_InAs + (1-self.moleFrac6)*c.EgG_AlAs - self.moleFrac6*(1-self.moleFrac6)*c.EgG_AlInAs;
            self.EgL[5] = self.moleFrac6*c.EgL_InAs + (1-self.moleFrac6)*c.EgL_AlAs - self.moleFrac6*(1-self.moleFrac6)*c.EgL_AlInAs;
            self.EgX[5] = self.moleFrac6*c.EgX_InAs + (1-self.moleFrac6)*c.EgX_AlAs - self.moleFrac6*(1-self.moleFrac6)*c.EgX_AlInAs;
            self.VBO[5] = self.moleFrac6*c.VBO_InAs + (1-self.moleFrac6)*c.VBO_AlAs - self.moleFrac6*(1-self.moleFrac6)*c.VBO_AlInAs;
            self.DSO[5] = self.moleFrac6*c.DSO_InAs + (1-self.moleFrac6)*c.DSO_AlAs - self.moleFrac6*(1-self.moleFrac6)*c.DSO_AlInAs;
            self.me0[5] = self.moleFrac6*c.me0_InAs + (1-self.moleFrac6)*c.me0_AlAs - self.moleFrac6*(1-self.moleFrac6)*c.me0_AlInAs;
            self.acG[5] = self.moleFrac6*c.acG_InAs + (1-self.moleFrac6)*c.acG_AlAs - self.moleFrac6*(1-self.moleFrac6)*c.acG_AlInAs;
            self.acL[5] = self.moleFrac6*c.acL_InAs + (1-self.moleFrac6)*c.acL_AlAs - self.moleFrac6*(1-self.moleFrac6)*c.acL_AlInAs;
            self.acX[5] = self.moleFrac6*c.acX_InAs + (1-self.moleFrac6)*c.acX_AlAs - self.moleFrac6*(1-self.moleFrac6)*c.acX_AlInAs;
            self.Ep[5]  = self.moleFrac6*c.Ep_InAs  + (1-self.moleFrac6)*c.Ep_AlAs  - self.moleFrac6*(1-self.moleFrac6)*c.Ep_AlInAs;
            self.F[5]   = self.moleFrac6*c.F_InAs   + (1-self.moleFrac6)*c.F_AlAs   - self.moleFrac6*(1-self.moleFrac6)*c.F_AlInAs;
            self.XiX[5] = self.moleFrac6*c.XiX_InAs + (1-self.moleFrac6)*c.XiX_AlAs;
            self.b[5]   = self.moleFrac6*c.b_InAs   + (1-self.moleFrac6)*c.b_AlAs;
            self.av[5]  = self.moleFrac6*c.av_InAs  + (1-self.moleFrac6)*c.av_AlAs;
            self.alG[5] = self.moleFrac6*c.alG_InAs + (1-self.moleFrac6)*c.alG_AlAs;
            self.beG[5] = self.moleFrac6*c.beG_InAs + (1-self.moleFrac6)*c.beG_AlAs;
            self.alL[5] = self.moleFrac6*c.alL_InAs + (1-self.moleFrac6)*c.alL_AlAs;
            self.beL[5] = self.moleFrac6*c.beL_InAs + (1-self.moleFrac6)*c.beL_AlAs;
            self.alX[5] = self.moleFrac6*c.alX_InAs + (1-self.moleFrac6)*c.alX_AlAs;
            self.beX[5] = self.moleFrac6*c.beX_InAs + (1-self.moleFrac6)*c.beX_AlAs;
            self.epss[5]   = self.moleFrac6*c.epss_InAs   + (1-self.moleFrac6)*c.epss_AlAs;
            self.epsInf[5] = self.moleFrac6*c.epsInf_InAs + (1-self.moleFrac6)*c.epsInf_AlAs;
            self.hwLO[5]   = self.moleFrac6*c.hwLO_InAs   + (1-self.moleFrac6)*c.hwLO_AlAs;
            self.a0[5]     = self.moleFrac6*c.alc_InAs    + (1-self.moleFrac6)*c.alc_AlAs;
            self.c11[5]    = self.moleFrac6*c.c11_InAs    + (1-self.moleFrac6)*c.c11_AlAs;
            self.c12[5]    = self.moleFrac6*c.c12_InAs    + (1-self.moleFrac6)*c.c12_AlAs;
        
            self.EgG[6] = self.moleFrac7*c.EgG_InAs + (1-self.moleFrac7)*c.EgG_GaAs - self.moleFrac7*(1-self.moleFrac7)*c.EgG_InGaAs;
            self.EgL[6] = self.moleFrac7*c.EgL_InAs + (1-self.moleFrac7)*c.EgL_GaAs - self.moleFrac7*(1-self.moleFrac7)*c.EgL_InGaAs;
            self.EgX[6] = self.moleFrac7*c.EgX_InAs + (1-self.moleFrac7)*c.EgX_GaAs - self.moleFrac7*(1-self.moleFrac7)*c.EgX_InGaAs;
            self.VBO[6] = self.moleFrac7*c.VBO_InAs + (1-self.moleFrac7)*c.VBO_GaAs - self.moleFrac7*(1-self.moleFrac7)*c.VBO_InGaAs;
            self.DSO[6] = self.moleFrac7*c.DSO_InAs + (1-self.moleFrac7)*c.DSO_GaAs - self.moleFrac7*(1-self.moleFrac7)*c.DSO_InGaAs;
            self.me0[6] = self.moleFrac7*c.me0_InAs + (1-self.moleFrac7)*c.me0_GaAs - self.moleFrac7*(1-self.moleFrac7)*c.me0_InGaAs;
            self.acG[6] = self.moleFrac7*c.acG_InAs + (1-self.moleFrac7)*c.acG_GaAs - self.moleFrac7*(1-self.moleFrac7)*c.acG_InGaAs;
            self.acL[6] = self.moleFrac7*c.acL_InAs + (1-self.moleFrac7)*c.acL_GaAs - self.moleFrac7*(1-self.moleFrac7)*c.acL_InGaAs;
            self.acX[6] = self.moleFrac7*c.acX_InAs + (1-self.moleFrac7)*c.acX_GaAs - self.moleFrac7*(1-self.moleFrac7)*c.acX_InGaAs;
            self.Ep[6]  = self.moleFrac7*c.Ep_InAs  + (1-self.moleFrac7)*c.Ep_GaAs  - self.moleFrac7*(1-self.moleFrac7)*c.Ep_InGaAs;
            self.F[6]   = self.moleFrac7*c.F_InAs   + (1-self.moleFrac7)*c.F_GaAs   - self.moleFrac7*(1-self.moleFrac7)*c.F_InGaAs;
            self.XiX[6] = self.moleFrac7*c.XiX_InAs + (1-self.moleFrac7)*c.XiX_GaAs;
            self.b[6]   = self.moleFrac7*c.b_InAs   + (1-self.moleFrac7)*c.b_GaAs;
            self.av[6]  = self.moleFrac7*c.av_InAs  + (1-self.moleFrac7)*c.av_GaAs;
            self.alG[6] = self.moleFrac7*c.alG_InAs + (1-self.moleFrac7)*c.alG_GaAs;
            self.beG[6] = self.moleFrac7*c.beG_InAs + (1-self.moleFrac7)*c.beG_GaAs;
            self.alL[6] = self.moleFrac7*c.alL_InAs + (1-self.moleFrac7)*c.alL_GaAs;
            self.beL[6] = self.moleFrac7*c.beL_InAs + (1-self.moleFrac7)*c.beL_GaAs;
            self.alX[6] = self.moleFrac7*c.alX_InAs + (1-self.moleFrac7)*c.alX_GaAs;
            self.beX[6] = self.moleFrac7*c.beX_InAs + (1-self.moleFrac7)*c.beX_GaAs;
            self.epss[6]   = self.moleFrac7*c.epss_InAs   + (1-self.moleFrac7)*c.epss_GaAs;
            self.epsInf[6] = self.moleFrac7*c.epsInf_InAs + (1-self.moleFrac7)*c.epsInf_GaAs;
            self.hwLO[6]   = self.moleFrac7*c.hwLO_InAs   + (1-self.moleFrac7)*c.hwLO_GaAs;
            self.a0[6]     = self.moleFrac7*c.alc_InAs    + (1-self.moleFrac7)*c.alc_GaAs;
            self.c11[6]    = self.moleFrac7*c.c11_InAs    + (1-self.moleFrac7)*c.c11_GaAs;
            self.c12[6]    = self.moleFrac7*c.c12_InAs    + (1-self.moleFrac7)*c.c12_GaAs;
            
            self.EgG[7] = self.moleFrac8*c.EgG_InAs + (1-self.moleFrac8)*c.EgG_AlAs - self.moleFrac8*(1-self.moleFrac8)*c.EgG_AlInAs;
            self.EgL[7] = self.moleFrac8*c.EgL_InAs + (1-self.moleFrac8)*c.EgL_AlAs - self.moleFrac8*(1-self.moleFrac8)*c.EgL_AlInAs;
            self.EgX[7] = self.moleFrac8*c.EgX_InAs + (1-self.moleFrac8)*c.EgX_AlAs - self.moleFrac8*(1-self.moleFrac8)*c.EgX_AlInAs;
            self.VBO[7] = self.moleFrac8*c.VBO_InAs + (1-self.moleFrac8)*c.VBO_AlAs - self.moleFrac8*(1-self.moleFrac8)*c.VBO_AlInAs;
            self.DSO[7] = self.moleFrac8*c.DSO_InAs + (1-self.moleFrac8)*c.DSO_AlAs - self.moleFrac8*(1-self.moleFrac8)*c.DSO_AlInAs;
            self.me0[7] = self.moleFrac8*c.me0_InAs + (1-self.moleFrac8)*c.me0_AlAs - self.moleFrac8*(1-self.moleFrac8)*c.me0_AlInAs;
            self.acG[7] = self.moleFrac8*c.acG_InAs + (1-self.moleFrac8)*c.acG_AlAs - self.moleFrac8*(1-self.moleFrac8)*c.acG_AlInAs;
            self.acL[7] = self.moleFrac8*c.acL_InAs + (1-self.moleFrac8)*c.acL_AlAs - self.moleFrac8*(1-self.moleFrac8)*c.acL_AlInAs;
            self.acX[7] = self.moleFrac8*c.acX_InAs + (1-self.moleFrac8)*c.acX_AlAs - self.moleFrac8*(1-self.moleFrac8)*c.acX_AlInAs;
            self.Ep[7]  = self.moleFrac8*c.Ep_InAs  + (1-self.moleFrac8)*c.Ep_AlAs  - self.moleFrac8*(1-self.moleFrac8)*c.Ep_AlInAs;
            self.F[7]   = self.moleFrac8*c.F_InAs   + (1-self.moleFrac8)*c.F_AlAs   - self.moleFrac8*(1-self.moleFrac8)*c.F_AlInAs;
            self.XiX[7] = self.moleFrac8*c.XiX_InAs + (1-self.moleFrac8)*c.XiX_AlAs;
            self.b[7]   = self.moleFrac8*c.b_InAs   + (1-self.moleFrac8)*c.b_AlAs;
            self.av[7]  = self.moleFrac8*c.av_InAs  + (1-self.moleFrac8)*c.av_AlAs;
            self.alG[7] = self.moleFrac8*c.alG_InAs + (1-self.moleFrac8)*c.alG_AlAs;
            self.beG[7] = self.moleFrac8*c.beG_InAs + (1-self.moleFrac8)*c.beG_AlAs;
            self.alL[7] = self.moleFrac8*c.alL_InAs + (1-self.moleFrac8)*c.alL_AlAs;
            self.beL[7] = self.moleFrac8*c.beL_InAs + (1-self.moleFrac8)*c.beL_AlAs;
            self.alX[7] = self.moleFrac8*c.alX_InAs + (1-self.moleFrac8)*c.alX_AlAs;
            self.beX[7] = self.moleFrac8*c.beX_InAs + (1-self.moleFrac8)*c.beX_AlAs;
            self.epss[7]   = self.moleFrac8*c.epss_InAs   + (1-self.moleFrac8)*c.epss_AlAs;
            self.epsInf[7] = self.moleFrac8*c.epsInf_InAs + (1-self.moleFrac8)*c.epsInf_AlAs;
            self.hwLO[7]   = self.moleFrac8*c.hwLO_InAs   + (1-self.moleFrac8)*c.hwLO_AlAs;
            self.a0[7] = self.moleFrac8*c.alc_InAs + (1-self.moleFrac8)*c.alc_AlAs; 
            self.c11[7]    = self.moleFrac8*c.c11_InAs    + (1-self.moleFrac8)*c.c11_AlAs;
            self.c12[7]    = self.moleFrac8*c.c12_InAs    + (1-self.moleFrac8)*c.c12_AlAs;
            
        if self.substrate == 'GaAs':
            self.numMaterials = 8
            variables = ['EgG', 'EgL', 'EgX', 'VBO', 'DSO', 'me0', 'acG', 'acL', 'acX', 
                     'Ep', 'F', 'XiX', 'b', 'av', 'alG', 'beG', 'alL', 'beL', 'alX', 
                     'beX', 'epss', 'epsInf', 'hwLO', 'a0', 'c11', 'c12', 'h']
            for item in variables:
                setattr(self, item, zeros(self.numMaterials))
            
            EgG_AlGaAs = c.EgG_AlGaAs + 1.310*self.moleFrac1
            self.EgG[0] = self.moleFrac1*c.EgG_AlAs + (1-self.moleFrac1)*c.EgG_GaAs - self.moleFrac1*(1-self.moleFrac1)*EgG_AlGaAs;
            self.EgL[0] = self.moleFrac1*c.EgL_AlAs + (1-self.moleFrac1)*c.EgL_GaAs - self.moleFrac1*(1-self.moleFrac1)*c.EgL_AlGaAs;
            self.EgX[0] = self.moleFrac1*c.EgX_AlAs + (1-self.moleFrac1)*c.EgX_GaAs - self.moleFrac1*(1-self.moleFrac1)*c.EgX_AlGaAs;
            self.VBO[0] = self.moleFrac1*c.VBO_AlAs + (1-self.moleFrac1)*c.VBO_GaAs - self.moleFrac1*(1-self.moleFrac1)*c.VBO_AlGaAs;
            self.DSO[0] = self.moleFrac1*c.DSO_AlAs + (1-self.moleFrac1)*c.DSO_GaAs - self.moleFrac1*(1-self.moleFrac1)*c.DSO_AlGaAs;
            self.me0[0] = self.moleFrac1*c.me0_AlAs + (1-self.moleFrac1)*c.me0_GaAs - self.moleFrac1*(1-self.moleFrac1)*c.me0_AlGaAs;
            self.acG[0] = self.moleFrac1*c.acG_AlAs + (1-self.moleFrac1)*c.acG_GaAs - self.moleFrac1*(1-self.moleFrac1)*c.acG_AlGaAs;
            self.acL[0] = self.moleFrac1*c.acL_AlAs + (1-self.moleFrac1)*c.acL_GaAs - self.moleFrac1*(1-self.moleFrac1)*c.acL_AlGaAs;
            self.acX[0] = self.moleFrac1*c.acX_AlAs + (1-self.moleFrac1)*c.acX_GaAs - self.moleFrac1*(1-self.moleFrac1)*c.acX_AlGaAs;
            self.Ep[0]  = self.moleFrac1*c.Ep_AlAs  + (1-self.moleFrac1)*c.Ep_GaAs  - self.moleFrac1*(1-self.moleFrac1)*c.Ep_AlGaAs;
            self.F[0]   = self.moleFrac1*c.F_AlAs   + (1-self.moleFrac1)*c.F_GaAs   - self.moleFrac1*(1-self.moleFrac1)*c.F_AlGaAs;
            self.XiX[0] = self.moleFrac1*c.XiX_AlAs + (1-self.moleFrac1)*c.XiX_GaAs;
            self.b[0]   = self.moleFrac1*c.b_AlAs   + (1-self.moleFrac1)*c.b_GaAs;
            self.av[0]  = self.moleFrac1*c.av_AlAs  + (1-self.moleFrac1)*c.av_GaAs;
            self.alG[0] = self.moleFrac1*c.alG_AlAs + (1-self.moleFrac1)*c.alG_GaAs;
            self.beG[0] = self.moleFrac1*c.beG_AlAs + (1-self.moleFrac1)*c.beG_GaAs;
            self.alL[0] = self.moleFrac1*c.alL_AlAs + (1-self.moleFrac1)*c.alL_GaAs;
            self.beL[0] = self.moleFrac1*c.beL_AlAs + (1-self.moleFrac1)*c.beL_GaAs;
            self.alX[0] = self.moleFrac1*c.alX_AlAs + (1-self.moleFrac1)*c.alX_GaAs;
            self.beX[0] = self.moleFrac1*c.beX_AlAs + (1-self.moleFrac1)*c.beX_GaAs;
            self.epss[0]   = self.moleFrac1*c.epss_AlAs   + (1-self.moleFrac1)*c.epss_GaAs;
            self.epsInf[0] = self.moleFrac1*c.epsInf_AlAs + (1-self.moleFrac1)*c.epsInf_GaAs;
            self.hwLO[0]   = self.moleFrac1*c.hwLO_AlAs   + (1-self.moleFrac1)*c.hwLO_GaAs;
            #self.Ped[0]   = 3*self.moleFrac1*(1-self.moleFrac1)*(-c.acG_InAs+c.acG_GaAs)*(c.alc_InAs-c.alc_GaAs)/c.alc_InP; #Ped is the third term from Walle eqn (21)
            self.a0[0]     = self.moleFrac1*c.alc_AlAs    + (1-self.moleFrac1)*c.alc_GaAs;
            self.c11[0]    = self.moleFrac1*c.c11_AlAs    + (1-self.moleFrac1)*c.c11_GaAs;
            self.c12[0]    = self.moleFrac1*c.c12_AlAs    + (1-self.moleFrac1)*c.c12_GaAs;
            
            EgG_AlGaAs = c.EgG_AlGaAs + 1.310*self.moleFrac2
            self.EgG[1] = self.moleFrac2*c.EgG_AlAs + (1-self.moleFrac2)*c.EgG_GaAs - self.moleFrac2*(1-self.moleFrac2)*EgG_AlGaAs;
            self.EgL[1] = self.moleFrac2*c.EgL_AlAs + (1-self.moleFrac2)*c.EgL_GaAs - self.moleFrac2*(1-self.moleFrac2)*c.EgL_AlGaAs;
            self.EgX[1] = self.moleFrac2*c.EgX_AlAs + (1-self.moleFrac2)*c.EgX_GaAs - self.moleFrac2*(1-self.moleFrac2)*c.EgX_AlGaAs;
            self.VBO[1] = self.moleFrac2*c.VBO_AlAs + (1-self.moleFrac2)*c.VBO_GaAs - self.moleFrac2*(1-self.moleFrac2)*c.VBO_AlGaAs;
            self.DSO[1] = self.moleFrac2*c.DSO_AlAs + (1-self.moleFrac2)*c.DSO_GaAs - self.moleFrac2*(1-self.moleFrac2)*c.DSO_AlGaAs;
            self.me0[1] = self.moleFrac2*c.me0_AlAs + (1-self.moleFrac2)*c.me0_GaAs - self.moleFrac2*(1-self.moleFrac2)*c.me0_AlGaAs;
            self.acG[1] = self.moleFrac2*c.acG_AlAs + (1-self.moleFrac2)*c.acG_GaAs - self.moleFrac2*(1-self.moleFrac2)*c.acG_AlGaAs;
            self.acL[1] = self.moleFrac2*c.acL_AlAs + (1-self.moleFrac2)*c.acL_GaAs - self.moleFrac2*(1-self.moleFrac2)*c.acL_AlGaAs;
            self.acX[1] = self.moleFrac2*c.acX_AlAs + (1-self.moleFrac2)*c.acX_GaAs - self.moleFrac2*(1-self.moleFrac2)*c.acX_AlGaAs;
            self.Ep[1]  = self.moleFrac2*c.Ep_AlAs  + (1-self.moleFrac2)*c.Ep_GaAs  - self.moleFrac2*(1-self.moleFrac2)*c.Ep_AlGaAs;
            self.F[1]   = self.moleFrac2*c.F_AlAs   + (1-self.moleFrac2)*c.F_GaAs   - self.moleFrac2*(1-self.moleFrac2)*c.F_AlGaAs;
            self.XiX[1] = self.moleFrac2*c.XiX_AlAs + (1-self.moleFrac2)*c.XiX_GaAs;
            self.b[1]   = self.moleFrac2*c.b_AlAs   + (1-self.moleFrac2)*c.b_GaAs;
            self.av[1]  = self.moleFrac2*c.av_AlAs  + (1-self.moleFrac2)*c.av_GaAs;
            self.alG[1] = self.moleFrac2*c.alG_AlAs + (1-self.moleFrac2)*c.alG_GaAs;
            self.beG[1] = self.moleFrac2*c.beG_AlAs + (1-self.moleFrac2)*c.beG_GaAs;
            self.alL[1] = self.moleFrac2*c.alL_AlAs + (1-self.moleFrac2)*c.alL_GaAs;
            self.beL[1] = self.moleFrac2*c.beL_AlAs + (1-self.moleFrac2)*c.beL_GaAs;
            self.alX[1] = self.moleFrac2*c.alX_AlAs + (1-self.moleFrac2)*c.alX_GaAs;
            self.beX[1] = self.moleFrac2*c.beX_AlAs + (1-self.moleFrac2)*c.beX_GaAs;
            self.epss[1]   = self.moleFrac2*c.epss_AlAs   + (1-self.moleFrac2)*c.epss_GaAs;
            self.epsInf[1] = self.moleFrac2*c.epsInf_AlAs + (1-self.moleFrac2)*c.epsInf_GaAs;
            self.hwLO[1]   = self.moleFrac2*c.hwLO_AlAs   + (1-self.moleFrac2)*c.hwLO_GaAs;
            self.a0[1]     = self.moleFrac2*c.alc_AlAs    + (1-self.moleFrac2)*c.alc_GaAs;
            self.c11[1]    = self.moleFrac2*c.c11_AlAs    + (1-self.moleFrac2)*c.c11_GaAs;
            self.c12[1]    = self.moleFrac2*c.c12_AlAs    + (1-self.moleFrac2)*c.c12_GaAs;
            
            EgG_AlGaAs = c.EgG_AlGaAs + 1.310*self.moleFrac3
            self.EgG[2] = self.moleFrac3*c.EgG_AlAs + (1-self.moleFrac3)*c.EgG_GaAs - self.moleFrac3*(1-self.moleFrac3)*EgG_AlGaAs;
            self.EgL[2] = self.moleFrac3*c.EgL_AlAs + (1-self.moleFrac3)*c.EgL_GaAs - self.moleFrac3*(1-self.moleFrac3)*c.EgL_AlGaAs;
            self.EgX[2] = self.moleFrac3*c.EgX_AlAs + (1-self.moleFrac3)*c.EgX_GaAs - self.moleFrac3*(1-self.moleFrac3)*c.EgX_AlGaAs;
            self.VBO[2] = self.moleFrac3*c.VBO_AlAs + (1-self.moleFrac3)*c.VBO_GaAs - self.moleFrac3*(1-self.moleFrac3)*c.VBO_AlGaAs;
            self.DSO[2] = self.moleFrac3*c.DSO_AlAs + (1-self.moleFrac3)*c.DSO_GaAs - self.moleFrac3*(1-self.moleFrac3)*c.DSO_AlGaAs;
            self.me0[2] = self.moleFrac3*c.me0_AlAs + (1-self.moleFrac3)*c.me0_GaAs - self.moleFrac3*(1-self.moleFrac3)*c.me0_AlGaAs;
            self.acG[2] = self.moleFrac3*c.acG_AlAs + (1-self.moleFrac3)*c.acG_GaAs - self.moleFrac3*(1-self.moleFrac3)*c.acG_AlGaAs;
            self.acL[2] = self.moleFrac3*c.acL_AlAs + (1-self.moleFrac3)*c.acL_GaAs - self.moleFrac3*(1-self.moleFrac3)*c.acL_AlGaAs;
            self.acX[2] = self.moleFrac3*c.acX_AlAs + (1-self.moleFrac3)*c.acX_GaAs - self.moleFrac3*(1-self.moleFrac3)*c.acX_AlGaAs;
            self.Ep[2]  = self.moleFrac3*c.Ep_AlAs  + (1-self.moleFrac3)*c.Ep_GaAs  - self.moleFrac3*(1-self.moleFrac3)*c.Ep_AlGaAs;
            self.F[2]   = self.moleFrac3*c.F_AlAs   + (1-self.moleFrac3)*c.F_GaAs   - self.moleFrac3*(1-self.moleFrac3)*c.F_AlGaAs;
            self.XiX[2] = self.moleFrac3*c.XiX_AlAs + (1-self.moleFrac3)*c.XiX_GaAs;
            self.b[2]   = self.moleFrac3*c.b_AlAs   + (1-self.moleFrac3)*c.b_GaAs;
            self.av[2]  = self.moleFrac3*c.av_AlAs  + (1-self.moleFrac3)*c.av_GaAs;
            self.alG[2] = self.moleFrac3*c.alG_AlAs + (1-self.moleFrac3)*c.alG_GaAs;
            self.beG[2] = self.moleFrac3*c.beG_AlAs + (1-self.moleFrac3)*c.beG_GaAs;
            self.alL[2] = self.moleFrac3*c.alL_AlAs + (1-self.moleFrac3)*c.alL_GaAs;
            self.beL[2] = self.moleFrac3*c.beL_AlAs + (1-self.moleFrac3)*c.beL_GaAs;
            self.alX[2] = self.moleFrac3*c.alX_AlAs + (1-self.moleFrac3)*c.alX_GaAs;
            self.beX[2] = self.moleFrac3*c.beX_AlAs + (1-self.moleFrac3)*c.beX_GaAs;
            self.epss[2]   = self.moleFrac3*c.epss_AlAs   + (1-self.moleFrac3)*c.epss_GaAs;
            self.epsInf[2] = self.moleFrac3*c.epsInf_AlAs + (1-self.moleFrac3)*c.epsInf_GaAs;
            self.hwLO[2]   = self.moleFrac3*c.hwLO_AlAs   + (1-self.moleFrac3)*c.hwLO_GaAs;
            self.a0[2]     = self.moleFrac3*c.alc_AlAs    + (1-self.moleFrac3)*c.alc_GaAs;
            self.c11[2]    = self.moleFrac3*c.c11_AlAs    + (1-self.moleFrac3)*c.c11_GaAs;
            self.c12[2]    = self.moleFrac3*c.c12_AlAs    + (1-self.moleFrac3)*c.c12_GaAs;
            
            EgG_AlGaAs = c.EgG_AlGaAs + 1.310*self.moleFrac4
            self.EgG[3] = self.moleFrac4*c.EgG_AlAs + (1-self.moleFrac4)*c.EgG_GaAs - self.moleFrac4*(1-self.moleFrac4)*EgG_AlGaAs;
            self.EgL[3] = self.moleFrac4*c.EgL_AlAs + (1-self.moleFrac4)*c.EgL_GaAs - self.moleFrac4*(1-self.moleFrac4)*c.EgL_AlGaAs;
            self.EgX[3] = self.moleFrac4*c.EgX_AlAs + (1-self.moleFrac4)*c.EgX_GaAs - self.moleFrac4*(1-self.moleFrac4)*c.EgX_AlGaAs;
            self.VBO[3] = self.moleFrac4*c.VBO_AlAs + (1-self.moleFrac4)*c.VBO_GaAs - self.moleFrac4*(1-self.moleFrac4)*c.VBO_AlGaAs;
            self.DSO[3] = self.moleFrac4*c.DSO_AlAs + (1-self.moleFrac4)*c.DSO_GaAs - self.moleFrac4*(1-self.moleFrac4)*c.DSO_AlGaAs;
            self.me0[3] = self.moleFrac4*c.me0_AlAs + (1-self.moleFrac4)*c.me0_GaAs - self.moleFrac4*(1-self.moleFrac4)*c.me0_AlGaAs;
            self.acG[3] = self.moleFrac4*c.acG_AlAs + (1-self.moleFrac4)*c.acG_GaAs - self.moleFrac4*(1-self.moleFrac4)*c.acG_AlGaAs;
            self.acL[3] = self.moleFrac4*c.acL_AlAs + (1-self.moleFrac4)*c.acL_GaAs - self.moleFrac4*(1-self.moleFrac4)*c.acL_AlGaAs;
            self.acX[3] = self.moleFrac4*c.acX_AlAs + (1-self.moleFrac4)*c.acX_GaAs - self.moleFrac4*(1-self.moleFrac4)*c.acX_AlGaAs;
            self.Ep[3]  = self.moleFrac4*c.Ep_AlAs  + (1-self.moleFrac4)*c.Ep_GaAs  - self.moleFrac4*(1-self.moleFrac4)*c.Ep_AlGaAs;
            self.F[3]   = self.moleFrac4*c.F_AlAs   + (1-self.moleFrac4)*c.F_GaAs   - self.moleFrac4*(1-self.moleFrac4)*c.F_AlGaAs;
            self.XiX[3] = self.moleFrac4*c.XiX_AlAs + (1-self.moleFrac4)*c.XiX_GaAs;
            self.b[3]   = self.moleFrac4*c.b_AlAs   + (1-self.moleFrac4)*c.b_GaAs;
            self.av[3]  = self.moleFrac4*c.av_AlAs  + (1-self.moleFrac4)*c.av_GaAs;
            self.alG[3] = self.moleFrac4*c.alG_AlAs + (1-self.moleFrac4)*c.alG_GaAs;
            self.beG[3] = self.moleFrac4*c.beG_AlAs + (1-self.moleFrac4)*c.beG_GaAs;
            self.alL[3] = self.moleFrac4*c.alL_AlAs + (1-self.moleFrac4)*c.alL_GaAs;
            self.beL[3] = self.moleFrac4*c.beL_AlAs + (1-self.moleFrac4)*c.beL_GaAs;
            self.alX[3] = self.moleFrac4*c.alX_AlAs + (1-self.moleFrac4)*c.alX_GaAs;
            self.beX[3] = self.moleFrac4*c.beX_AlAs + (1-self.moleFrac4)*c.beX_GaAs;
            self.epss[3]   = self.moleFrac4*c.epss_AlAs   + (1-self.moleFrac4)*c.epss_GaAs;
            self.epsInf[3] = self.moleFrac4*c.epsInf_AlAs + (1-self.moleFrac4)*c.epsInf_GaAs;
            self.hwLO[3]   = self.moleFrac4*c.hwLO_AlAs   + (1-self.moleFrac4)*c.hwLO_GaAs;
            self.a0[3]     = self.moleFrac4*c.alc_AlAs    + (1-self.moleFrac4)*c.alc_GaAs;
            self.c11[3]    = self.moleFrac4*c.c11_AlAs    + (1-self.moleFrac4)*c.c11_GaAs;
            self.c12[3]    = self.moleFrac4*c.c12_AlAs    + (1-self.moleFrac4)*c.c12_GaAs;
            
            EgG_AlGaAs = c.EgG_AlGaAs + 1.310*self.moleFrac5
            self.EgG[4] = self.moleFrac5*c.EgG_AlAs + (1-self.moleFrac5)*c.EgG_GaAs - self.moleFrac5*(1-self.moleFrac5)*EgG_AlGaAs;
            self.EgL[4] = self.moleFrac5*c.EgL_AlAs + (1-self.moleFrac5)*c.EgL_GaAs - self.moleFrac5*(1-self.moleFrac5)*c.EgL_AlGaAs;
            self.EgX[4] = self.moleFrac5*c.EgX_AlAs + (1-self.moleFrac5)*c.EgX_GaAs - self.moleFrac5*(1-self.moleFrac5)*c.EgX_AlGaAs;
            self.VBO[4] = self.moleFrac5*c.VBO_AlAs + (1-self.moleFrac5)*c.VBO_GaAs - self.moleFrac5*(1-self.moleFrac5)*c.VBO_AlGaAs;
            self.DSO[4] = self.moleFrac5*c.DSO_AlAs + (1-self.moleFrac5)*c.DSO_GaAs - self.moleFrac5*(1-self.moleFrac5)*c.DSO_AlGaAs;
            self.me0[4] = self.moleFrac5*c.me0_AlAs + (1-self.moleFrac5)*c.me0_GaAs - self.moleFrac5*(1-self.moleFrac5)*c.me0_AlGaAs;
            self.acG[4] = self.moleFrac5*c.acG_AlAs + (1-self.moleFrac5)*c.acG_GaAs - self.moleFrac5*(1-self.moleFrac5)*c.acG_AlGaAs;
            self.acL[4] = self.moleFrac5*c.acL_AlAs + (1-self.moleFrac5)*c.acL_GaAs - self.moleFrac5*(1-self.moleFrac5)*c.acL_AlGaAs;
            self.acX[4] = self.moleFrac5*c.acX_AlAs + (1-self.moleFrac5)*c.acX_GaAs - self.moleFrac5*(1-self.moleFrac5)*c.acX_AlGaAs;
            self.Ep[4]  = self.moleFrac5*c.Ep_AlAs  + (1-self.moleFrac5)*c.Ep_GaAs  - self.moleFrac5*(1-self.moleFrac5)*c.Ep_AlGaAs;
            self.F[4]   = self.moleFrac5*c.F_AlAs   + (1-self.moleFrac5)*c.F_GaAs   - self.moleFrac5*(1-self.moleFrac5)*c.F_AlGaAs;
            self.XiX[4] = self.moleFrac5*c.XiX_AlAs + (1-self.moleFrac5)*c.XiX_GaAs;
            self.b[4]   = self.moleFrac5*c.b_AlAs   + (1-self.moleFrac5)*c.b_GaAs;
            self.av[4]  = self.moleFrac5*c.av_AlAs  + (1-self.moleFrac5)*c.av_GaAs;
            self.alG[4] = self.moleFrac5*c.alG_AlAs + (1-self.moleFrac5)*c.alG_GaAs;
            self.beG[4] = self.moleFrac5*c.beG_AlAs + (1-self.moleFrac5)*c.beG_GaAs;
            self.alL[4] = self.moleFrac5*c.alL_AlAs + (1-self.moleFrac5)*c.alL_GaAs;
            self.beL[4] = self.moleFrac5*c.beL_AlAs + (1-self.moleFrac5)*c.beL_GaAs;
            self.alX[4] = self.moleFrac5*c.alX_AlAs + (1-self.moleFrac5)*c.alX_GaAs;
            self.beX[4] = self.moleFrac5*c.beX_AlAs + (1-self.moleFrac5)*c.beX_GaAs;
            self.epss[4]   = self.moleFrac5*c.epss_AlAs   + (1-self.moleFrac5)*c.epss_GaAs;
            self.epsInf[4] = self.moleFrac5*c.epsInf_AlAs + (1-self.moleFrac5)*c.epsInf_GaAs;
            self.hwLO[4]   = self.moleFrac5*c.hwLO_AlAs   + (1-self.moleFrac5)*c.hwLO_GaAs;
            self.a0[4]     = self.moleFrac5*c.alc_AlAs    + (1-self.moleFrac5)*c.alc_GaAs;
            self.c11[4]    = self.moleFrac5*c.c11_AlAs    + (1-self.moleFrac5)*c.c11_GaAs;
            self.c12[4]    = self.moleFrac5*c.c12_AlAs    + (1-self.moleFrac5)*c.c12_GaAs;
            
            EgG_AlGaAs = c.EgG_AlGaAs + 1.310*self.moleFrac6
            self.EgG[5] = self.moleFrac6*c.EgG_AlAs + (1-self.moleFrac6)*c.EgG_GaAs - self.moleFrac6*(1-self.moleFrac6)*EgG_AlGaAs;
            self.EgL[5] = self.moleFrac6*c.EgL_AlAs + (1-self.moleFrac6)*c.EgL_GaAs - self.moleFrac6*(1-self.moleFrac6)*c.EgL_AlGaAs;
            self.EgX[5] = self.moleFrac6*c.EgX_AlAs + (1-self.moleFrac6)*c.EgX_GaAs - self.moleFrac6*(1-self.moleFrac6)*c.EgX_AlGaAs;
            self.VBO[5] = self.moleFrac6*c.VBO_AlAs + (1-self.moleFrac6)*c.VBO_GaAs - self.moleFrac6*(1-self.moleFrac6)*c.VBO_AlGaAs;
            self.DSO[5] = self.moleFrac6*c.DSO_AlAs + (1-self.moleFrac6)*c.DSO_GaAs - self.moleFrac6*(1-self.moleFrac6)*c.DSO_AlGaAs;
            self.me0[5] = self.moleFrac6*c.me0_AlAs + (1-self.moleFrac6)*c.me0_GaAs - self.moleFrac6*(1-self.moleFrac6)*c.me0_AlGaAs;
            self.acG[5] = self.moleFrac6*c.acG_AlAs + (1-self.moleFrac6)*c.acG_GaAs - self.moleFrac6*(1-self.moleFrac6)*c.acG_AlGaAs;
            self.acL[5] = self.moleFrac6*c.acL_AlAs + (1-self.moleFrac6)*c.acL_GaAs - self.moleFrac6*(1-self.moleFrac6)*c.acL_AlGaAs;
            self.acX[5] = self.moleFrac6*c.acX_AlAs + (1-self.moleFrac6)*c.acX_GaAs - self.moleFrac6*(1-self.moleFrac6)*c.acX_AlGaAs;
            self.Ep[5]  = self.moleFrac6*c.Ep_AlAs  + (1-self.moleFrac6)*c.Ep_GaAs  - self.moleFrac6*(1-self.moleFrac6)*c.Ep_AlGaAs;
            self.F[5]   = self.moleFrac6*c.F_AlAs   + (1-self.moleFrac6)*c.F_GaAs   - self.moleFrac6*(1-self.moleFrac6)*c.F_AlGaAs;
            self.XiX[5] = self.moleFrac6*c.XiX_AlAs + (1-self.moleFrac6)*c.XiX_GaAs;
            self.b[5]   = self.moleFrac6*c.b_AlAs   + (1-self.moleFrac6)*c.b_GaAs;
            self.av[5]  = self.moleFrac6*c.av_AlAs  + (1-self.moleFrac6)*c.av_GaAs;
            self.alG[5] = self.moleFrac6*c.alG_AlAs + (1-self.moleFrac6)*c.alG_GaAs;
            self.beG[5] = self.moleFrac6*c.beG_AlAs + (1-self.moleFrac6)*c.beG_GaAs;
            self.alL[5] = self.moleFrac6*c.alL_AlAs + (1-self.moleFrac6)*c.alL_GaAs;
            self.beL[5] = self.moleFrac6*c.beL_AlAs + (1-self.moleFrac6)*c.beL_GaAs;
            self.alX[5] = self.moleFrac6*c.alX_AlAs + (1-self.moleFrac6)*c.alX_GaAs;
            self.beX[5] = self.moleFrac6*c.beX_AlAs + (1-self.moleFrac6)*c.beX_GaAs;
            self.epss[5]   = self.moleFrac6*c.epss_AlAs   + (1-self.moleFrac6)*c.epss_GaAs;
            self.epsInf[5] = self.moleFrac6*c.epsInf_AlAs + (1-self.moleFrac6)*c.epsInf_GaAs;
            self.hwLO[5]   = self.moleFrac6*c.hwLO_AlAs   + (1-self.moleFrac6)*c.hwLO_GaAs;
            self.a0[5]     = self.moleFrac6*c.alc_AlAs    + (1-self.moleFrac6)*c.alc_GaAs;
            self.c11[5]    = self.moleFrac6*c.c11_AlAs    + (1-self.moleFrac6)*c.c11_GaAs;
            self.c12[5]    = self.moleFrac6*c.c12_AlAs    + (1-self.moleFrac6)*c.c12_GaAs;
            
            EgG_AlGaAs = c.EgG_AlGaAs + 1.310*self.moleFrac7
            self.EgG[6] = self.moleFrac7*c.EgG_AlAs + (1-self.moleFrac7)*c.EgG_GaAs - self.moleFrac7*(1-self.moleFrac7)*EgG_AlGaAs;
            self.EgL[6] = self.moleFrac7*c.EgL_AlAs + (1-self.moleFrac7)*c.EgL_GaAs - self.moleFrac7*(1-self.moleFrac7)*c.EgL_AlGaAs;
            self.EgX[6] = self.moleFrac7*c.EgX_AlAs + (1-self.moleFrac7)*c.EgX_GaAs - self.moleFrac7*(1-self.moleFrac7)*c.EgX_AlGaAs;
            self.VBO[6] = self.moleFrac7*c.VBO_AlAs + (1-self.moleFrac7)*c.VBO_GaAs - self.moleFrac7*(1-self.moleFrac7)*c.VBO_AlGaAs;
            self.DSO[6] = self.moleFrac7*c.DSO_AlAs + (1-self.moleFrac7)*c.DSO_GaAs - self.moleFrac7*(1-self.moleFrac7)*c.DSO_AlGaAs;
            self.me0[6] = self.moleFrac7*c.me0_AlAs + (1-self.moleFrac7)*c.me0_GaAs - self.moleFrac7*(1-self.moleFrac7)*c.me0_AlGaAs;
            self.acG[6] = self.moleFrac7*c.acG_AlAs + (1-self.moleFrac7)*c.acG_GaAs - self.moleFrac7*(1-self.moleFrac7)*c.acG_AlGaAs;
            self.acL[6] = self.moleFrac7*c.acL_AlAs + (1-self.moleFrac7)*c.acL_GaAs - self.moleFrac7*(1-self.moleFrac7)*c.acL_AlGaAs;
            self.acX[6] = self.moleFrac7*c.acX_AlAs + (1-self.moleFrac7)*c.acX_GaAs - self.moleFrac7*(1-self.moleFrac7)*c.acX_AlGaAs;
            self.Ep[6]  = self.moleFrac7*c.Ep_AlAs  + (1-self.moleFrac7)*c.Ep_GaAs  - self.moleFrac7*(1-self.moleFrac7)*c.Ep_AlGaAs;
            self.F[6]   = self.moleFrac7*c.F_AlAs   + (1-self.moleFrac7)*c.F_GaAs   - self.moleFrac7*(1-self.moleFrac7)*c.F_AlGaAs;
            self.XiX[6] = self.moleFrac7*c.XiX_AlAs + (1-self.moleFrac7)*c.XiX_GaAs;
            self.b[6]   = self.moleFrac7*c.b_AlAs   + (1-self.moleFrac7)*c.b_GaAs;
            self.av[6]  = self.moleFrac7*c.av_AlAs  + (1-self.moleFrac7)*c.av_GaAs;
            self.alG[6] = self.moleFrac7*c.alG_AlAs + (1-self.moleFrac7)*c.alG_GaAs;
            self.beG[6] = self.moleFrac7*c.beG_AlAs + (1-self.moleFrac7)*c.beG_GaAs;
            self.alL[6] = self.moleFrac7*c.alL_AlAs + (1-self.moleFrac7)*c.alL_GaAs;
            self.beL[6] = self.moleFrac7*c.beL_AlAs + (1-self.moleFrac7)*c.beL_GaAs;
            self.alX[6] = self.moleFrac7*c.alX_AlAs + (1-self.moleFrac7)*c.alX_GaAs;
            self.beX[6] = self.moleFrac7*c.beX_AlAs + (1-self.moleFrac7)*c.beX_GaAs;
            self.epss[6]   = self.moleFrac7*c.epss_AlAs   + (1-self.moleFrac7)*c.epss_GaAs;
            self.epsInf[6] = self.moleFrac7*c.epsInf_AlAs + (1-self.moleFrac7)*c.epsInf_GaAs;
            self.hwLO[6]   = self.moleFrac7*c.hwLO_AlAs   + (1-self.moleFrac7)*c.hwLO_GaAs;
            self.a0[6]     = self.moleFrac7*c.alc_AlAs    + (1-self.moleFrac7)*c.alc_GaAs;
            self.c11[6]    = self.moleFrac7*c.c11_AlAs    + (1-self.moleFrac7)*c.c11_GaAs;
            self.c12[6]    = self.moleFrac7*c.c12_AlAs    + (1-self.moleFrac7)*c.c12_GaAs;
            
            EgG_AlGaAs = c.EgG_AlGaAs + 1.310*self.moleFrac8
            self.EgG[7] = self.moleFrac8*c.EgG_AlAs + (1-self.moleFrac8)*c.EgG_GaAs - self.moleFrac8*(1-self.moleFrac8)*EgG_AlGaAs;
            self.EgL[7] = self.moleFrac8*c.EgL_AlAs + (1-self.moleFrac8)*c.EgL_GaAs - self.moleFrac8*(1-self.moleFrac8)*c.EgL_AlGaAs;
            self.EgX[7] = self.moleFrac8*c.EgX_AlAs + (1-self.moleFrac8)*c.EgX_GaAs - self.moleFrac8*(1-self.moleFrac8)*c.EgX_AlGaAs;
            self.VBO[7] = self.moleFrac8*c.VBO_AlAs + (1-self.moleFrac8)*c.VBO_GaAs - self.moleFrac8*(1-self.moleFrac8)*c.VBO_AlGaAs;
            self.DSO[7] = self.moleFrac8*c.DSO_AlAs + (1-self.moleFrac8)*c.DSO_GaAs - self.moleFrac8*(1-self.moleFrac8)*c.DSO_AlGaAs;
            self.me0[7] = self.moleFrac8*c.me0_AlAs + (1-self.moleFrac8)*c.me0_GaAs - self.moleFrac8*(1-self.moleFrac8)*c.me0_AlGaAs;
            self.acG[7] = self.moleFrac8*c.acG_AlAs + (1-self.moleFrac8)*c.acG_GaAs - self.moleFrac8*(1-self.moleFrac8)*c.acG_AlGaAs;
            self.acL[7] = self.moleFrac8*c.acL_AlAs + (1-self.moleFrac8)*c.acL_GaAs - self.moleFrac8*(1-self.moleFrac8)*c.acL_AlGaAs;
            self.acX[7] = self.moleFrac8*c.acX_AlAs + (1-self.moleFrac8)*c.acX_GaAs - self.moleFrac8*(1-self.moleFrac8)*c.acX_AlGaAs;
            self.Ep[7]  = self.moleFrac8*c.Ep_AlAs  + (1-self.moleFrac8)*c.Ep_GaAs  - self.moleFrac8*(1-self.moleFrac8)*c.Ep_AlGaAs;
            self.F[7]   = self.moleFrac8*c.F_AlAs   + (1-self.moleFrac8)*c.F_GaAs   - self.moleFrac8*(1-self.moleFrac8)*c.F_AlGaAs;
            self.XiX[7] = self.moleFrac8*c.XiX_AlAs + (1-self.moleFrac8)*c.XiX_GaAs;
            self.b[7]   = self.moleFrac8*c.b_AlAs   + (1-self.moleFrac8)*c.b_GaAs;
            self.av[7]  = self.moleFrac8*c.av_AlAs  + (1-self.moleFrac8)*c.av_GaAs;
            self.alG[7] = self.moleFrac8*c.alG_AlAs + (1-self.moleFrac8)*c.alG_GaAs;
            self.beG[7] = self.moleFrac8*c.beG_AlAs + (1-self.moleFrac8)*c.beG_GaAs;
            self.alL[7] = self.moleFrac8*c.alL_AlAs + (1-self.moleFrac8)*c.alL_GaAs;
            self.beL[7] = self.moleFrac8*c.beL_AlAs + (1-self.moleFrac8)*c.beL_GaAs;
            self.alX[7] = self.moleFrac8*c.alX_AlAs + (1-self.moleFrac8)*c.alX_GaAs;
            self.beX[7] = self.moleFrac8*c.beX_AlAs + (1-self.moleFrac8)*c.beX_GaAs;
            self.epss[7]   = self.moleFrac8*c.epss_AlAs   + (1-self.moleFrac8)*c.epss_GaAs;
            self.epsInf[7] = self.moleFrac8*c.epsInf_AlAs + (1-self.moleFrac8)*c.epsInf_GaAs;
            self.hwLO[7]   = self.moleFrac8*c.hwLO_AlAs   + (1-self.moleFrac8)*c.hwLO_GaAs;
            self.a0[7]     = self.moleFrac8*c.alc_AlAs    + (1-self.moleFrac8)*c.alc_GaAs;
            self.c11[7]    = self.moleFrac8*c.c11_AlAs    + (1-self.moleFrac8)*c.c11_GaAs;
            self.c12[7]    = self.moleFrac8*c.c12_AlAs    + (1-self.moleFrac8)*c.c12_GaAs;
            
        elif self.substrate == 'GaSb':
            self.numMaterials = 8
            variables = ['EgG', 'EgL', 'EgX', 'VBO', 'DSO', 'me0', 'acG', 'acL', 'acX', 
                     'Ep', 'F', 'XiX', 'b', 'av', 'alG', 'beG', 'alL', 'beL', 'alX', 
                     'beX', 'epss', 'epsInf', 'hwLO', 'a0', 'c11', 'c12', 'h']
            for item in variables:
                setattr(self, item, zeros(self.numMaterials))
                
            self.EgG[0] = self.moleFrac1*c.EgG_InAs + (1-self.moleFrac1)*c.EgG_InSb - self.moleFrac1*(1-self.moleFrac1)*c.EgG_InAsSb;
            self.EgL[0] = self.moleFrac1*c.EgL_InAs + (1-self.moleFrac1)*c.EgL_InSb - self.moleFrac1*(1-self.moleFrac1)*c.EgL_InAsSb;
            self.EgX[0] = self.moleFrac1*c.EgX_InAs + (1-self.moleFrac1)*c.EgX_InSb - self.moleFrac1*(1-self.moleFrac1)*c.EgX_InAsSb;
            self.VBO[0] = self.moleFrac1*c.VBO_InAs + (1-self.moleFrac1)*c.VBO_InSb - self.moleFrac1*(1-self.moleFrac1)*c.VBO_InAsSb;
            self.DSO[0] = self.moleFrac1*c.DSO_InAs + (1-self.moleFrac1)*c.DSO_InSb - self.moleFrac1*(1-self.moleFrac1)*c.DSO_InAsSb;
            self.me0[0] = self.moleFrac1*c.me0_InAs + (1-self.moleFrac1)*c.me0_InSb - self.moleFrac1*(1-self.moleFrac1)*c.me0_InAsSb;
            self.acG[0] = self.moleFrac1*c.acG_InAs + (1-self.moleFrac1)*c.acG_InSb - self.moleFrac1*(1-self.moleFrac1)*c.acG_InAsSb;
            self.acL[0] = self.moleFrac1*c.acL_InAs + (1-self.moleFrac1)*c.acL_InSb - self.moleFrac1*(1-self.moleFrac1)*c.acL_InAsSb;
            self.acX[0] = self.moleFrac1*c.acX_InAs + (1-self.moleFrac1)*c.acX_InSb - self.moleFrac1*(1-self.moleFrac1)*c.acX_InAsSb;
            self.Ep[0]  = self.moleFrac1*c.Ep_InAs  + (1-self.moleFrac1)*c.Ep_InSb  - self.moleFrac1*(1-self.moleFrac1)*c.Ep_InAsSb;
            self.F[0]   = self.moleFrac1*c.F_InAs   + (1-self.moleFrac1)*c.F_InSb   - self.moleFrac1*(1-self.moleFrac1)*c.F_InAsSb;
            self.XiX[0] = self.moleFrac1*c.XiX_InAs + (1-self.moleFrac1)*c.XiX_InSb;
            self.b[0]   = self.moleFrac1*c.b_InAs   + (1-self.moleFrac1)*c.b_InSb;
            self.av[0]  = self.moleFrac1*c.av_InAs  + (1-self.moleFrac1)*c.av_InSb;
            self.alG[0] = self.moleFrac1*c.alG_InAs + (1-self.moleFrac1)*c.alG_InSb;
            self.beG[0] = self.moleFrac1*c.beG_InAs + (1-self.moleFrac1)*c.beG_InSb;
            self.alL[0] = self.moleFrac1*c.alL_InAs + (1-self.moleFrac1)*c.alL_InSb;
            self.beL[0] = self.moleFrac1*c.beL_InAs + (1-self.moleFrac1)*c.beL_InSb;
            self.alX[0] = self.moleFrac1*c.alX_InAs + (1-self.moleFrac1)*c.alX_InSb;
            self.beX[0] = self.moleFrac1*c.beX_InAs + (1-self.moleFrac1)*c.beX_InSb;
            self.epss[0]   = self.moleFrac1*c.epss_InAs   + (1-self.moleFrac1)*c.epss_InSb;
            self.epsInf[0] = self.moleFrac1*c.epsInf_InAs + (1-self.moleFrac1)*c.epsInf_InSb;
            self.hwLO[0]   = self.moleFrac1*c.hwLO_InAs   + (1-self.moleFrac1)*c.hwLO_InSb;
            self.a0[0]     = self.moleFrac1*c.alc_InAs    + (1-self.moleFrac1)*c.alc_InSb;
            self.c11[0]    = self.moleFrac1*c.c11_InAs    + (1-self.moleFrac1)*c.c11_InSb;
            self.c12[0]    = self.moleFrac1*c.c12_InAs    + (1-self.moleFrac1)*c.c12_InSb;
            
            EgG_AlGaSb = c.EgG_AlGaSb + 1.22*self.moleFrac2
            self.EgG[1] = self.moleFrac2*c.EgG_AlSb + (1-self.moleFrac2)*c.EgG_GaSb - self.moleFrac2*(1-self.moleFrac2)*EgG_AlGaSb;
            self.EgL[1] = self.moleFrac2*c.EgL_AlSb + (1-self.moleFrac2)*c.EgL_GaSb - self.moleFrac2*(1-self.moleFrac2)*c.EgL_AlGaSb;
            self.EgX[1] = self.moleFrac2*c.EgX_AlSb + (1-self.moleFrac2)*c.EgX_GaSb - self.moleFrac2*(1-self.moleFrac2)*c.EgX_AlGaSb;
            self.VBO[1] = self.moleFrac2*c.VBO_AlSb + (1-self.moleFrac2)*c.VBO_GaSb - self.moleFrac2*(1-self.moleFrac2)*c.VBO_AlGaSb;
            self.DSO[1] = self.moleFrac2*c.DSO_AlSb + (1-self.moleFrac2)*c.DSO_GaSb - self.moleFrac2*(1-self.moleFrac2)*c.DSO_AlGaSb;
            self.me0[1] = self.moleFrac2*c.me0_AlSb + (1-self.moleFrac2)*c.me0_GaSb - self.moleFrac2*(1-self.moleFrac2)*c.me0_AlGaSb;
            self.acG[1] = self.moleFrac2*c.acG_AlSb + (1-self.moleFrac2)*c.acG_GaSb - self.moleFrac2*(1-self.moleFrac2)*c.acG_AlGaSb;
            self.acL[1] = self.moleFrac2*c.acL_AlSb + (1-self.moleFrac2)*c.acL_GaSb - self.moleFrac2*(1-self.moleFrac2)*c.acL_AlGaSb;
            self.acX[1] = self.moleFrac2*c.acX_AlSb + (1-self.moleFrac2)*c.acX_GaSb - self.moleFrac2*(1-self.moleFrac2)*c.acX_AlGaSb;
            self.Ep[1]  = self.moleFrac2*c.Ep_AlSb  + (1-self.moleFrac2)*c.Ep_GaSb  - self.moleFrac2*(1-self.moleFrac2)*c.Ep_AlGaSb;
            self.F[1]   = self.moleFrac2*c.F_AlSb   + (1-self.moleFrac2)*c.F_GaSb   - self.moleFrac2*(1-self.moleFrac2)*c.F_AlGaSb;
            self.XiX[1] = self.moleFrac2*c.XiX_AlSb + (1-self.moleFrac2)*c.XiX_GaSb;
            self.b[1]   = self.moleFrac2*c.b_AlSb   + (1-self.moleFrac2)*c.b_GaSb;
            self.av[1]  = self.moleFrac2*c.av_AlSb  + (1-self.moleFrac2)*c.av_GaSb;
            self.alG[1] = self.moleFrac2*c.alG_AlSb + (1-self.moleFrac2)*c.alG_GaSb;
            self.beG[1] = self.moleFrac2*c.beG_AlSb + (1-self.moleFrac2)*c.beG_GaSb;
            self.alL[1] = self.moleFrac2*c.alL_AlSb + (1-self.moleFrac2)*c.alL_GaSb;
            self.beL[1] = self.moleFrac2*c.beL_AlSb + (1-self.moleFrac2)*c.beL_GaSb;
            self.alX[1] = self.moleFrac2*c.alX_AlSb + (1-self.moleFrac2)*c.alX_GaSb;
            self.beX[1] = self.moleFrac2*c.beX_AlSb + (1-self.moleFrac2)*c.beX_GaSb;
            self.epss[1]   = self.moleFrac2*c.epss_AlSb   + (1-self.moleFrac2)*c.epss_GaSb;
            self.epsInf[1] = self.moleFrac2*c.epsInf_AlSb + (1-self.moleFrac2)*c.epsInf_GaSb;
            self.hwLO[1]   = self.moleFrac2*c.hwLO_AlSb   + (1-self.moleFrac2)*c.hwLO_GaSb;
            self.a0[1]     = self.moleFrac2*c.alc_AlSb    + (1-self.moleFrac2)*c.alc_GaSb;
            self.c11[1]    = self.moleFrac2*c.c11_AlSb    + (1-self.moleFrac2)*c.c11_GaSb;
            self.c12[1]    = self.moleFrac2*c.c12_AlSb    + (1-self.moleFrac2)*c.c12_GaSb;
                
            self.EgG[2] = self.moleFrac3*c.EgG_InAs + (1-self.moleFrac3)*c.EgG_InSb - self.moleFrac3*(1-self.moleFrac3)*c.EgG_InAsSb;
            self.EgL[2] = self.moleFrac3*c.EgL_InAs + (1-self.moleFrac3)*c.EgL_InSb - self.moleFrac3*(1-self.moleFrac3)*c.EgL_InAsSb;
            self.EgX[2] = self.moleFrac3*c.EgX_InAs + (1-self.moleFrac3)*c.EgX_InSb - self.moleFrac3*(1-self.moleFrac3)*c.EgX_InAsSb;
            self.VBO[2] = self.moleFrac3*c.VBO_InAs + (1-self.moleFrac3)*c.VBO_InSb - self.moleFrac3*(1-self.moleFrac3)*c.VBO_InAsSb;
            self.DSO[2] = self.moleFrac3*c.DSO_InAs + (1-self.moleFrac3)*c.DSO_InSb - self.moleFrac3*(1-self.moleFrac3)*c.DSO_InAsSb;
            self.me0[2] = self.moleFrac3*c.me0_InAs + (1-self.moleFrac3)*c.me0_InSb - self.moleFrac3*(1-self.moleFrac3)*c.me0_InAsSb;
            self.acG[2] = self.moleFrac3*c.acG_InAs + (1-self.moleFrac3)*c.acG_InSb - self.moleFrac3*(1-self.moleFrac3)*c.acG_InAsSb;
            self.acL[2] = self.moleFrac3*c.acL_InAs + (1-self.moleFrac3)*c.acL_InSb - self.moleFrac3*(1-self.moleFrac3)*c.acL_InAsSb;
            self.acX[2] = self.moleFrac3*c.acX_InAs + (1-self.moleFrac3)*c.acX_InSb - self.moleFrac3*(1-self.moleFrac3)*c.acX_InAsSb;
            self.Ep[2]  = self.moleFrac3*c.Ep_InAs  + (1-self.moleFrac3)*c.Ep_InSb  - self.moleFrac3*(1-self.moleFrac3)*c.Ep_InAsSb;
            self.F[2]   = self.moleFrac3*c.F_InAs   + (1-self.moleFrac3)*c.F_InSb   - self.moleFrac3*(1-self.moleFrac3)*c.F_InAsSb;
            self.XiX[2] = self.moleFrac3*c.XiX_InAs + (1-self.moleFrac3)*c.XiX_InSb;
            self.b[2]   = self.moleFrac3*c.b_InAs   + (1-self.moleFrac3)*c.b_InSb;
            self.av[2]  = self.moleFrac3*c.av_InAs  + (1-self.moleFrac3)*c.av_InSb;
            self.alG[2] = self.moleFrac3*c.alG_InAs + (1-self.moleFrac3)*c.alG_InSb;
            self.beG[2] = self.moleFrac3*c.beG_InAs + (1-self.moleFrac3)*c.beG_InSb;
            self.alL[2] = self.moleFrac3*c.alL_InAs + (1-self.moleFrac3)*c.alL_InSb;
            self.beL[2] = self.moleFrac3*c.beL_InAs + (1-self.moleFrac3)*c.beL_InSb;
            self.alX[2] = self.moleFrac3*c.alX_InAs + (1-self.moleFrac3)*c.alX_InSb;
            self.beX[2] = self.moleFrac3*c.beX_InAs + (1-self.moleFrac3)*c.beX_InSb;
            self.epss[2]   = self.moleFrac3*c.epss_InAs   + (1-self.moleFrac3)*c.epss_InSb;
            self.epsInf[2] = self.moleFrac3*c.epsInf_InAs + (1-self.moleFrac3)*c.epsInf_InSb;
            self.hwLO[2]   = self.moleFrac3*c.hwLO_InAs   + (1-self.moleFrac3)*c.hwLO_InSb;
            self.a0[2]     = self.moleFrac3*c.alc_InAs    + (1-self.moleFrac3)*c.alc_InSb;
            self.c11[2]    = self.moleFrac3*c.c11_InAs    + (1-self.moleFrac3)*c.c11_InSb;
            self.c12[2]    = self.moleFrac3*c.c12_InAs    + (1-self.moleFrac3)*c.c12_InSb;
            
            EgG_AlGaSb = c.EgG_AlGaSb + 1.22*self.moleFrac4
            self.EgG[3] = self.moleFrac4*c.EgG_AlSb + (1-self.moleFrac4)*c.EgG_GaSb - self.moleFrac4*(1-self.moleFrac4)*EgG_AlGaSb;
            self.EgL[3] = self.moleFrac4*c.EgL_AlSb + (1-self.moleFrac4)*c.EgL_GaSb - self.moleFrac4*(1-self.moleFrac4)*c.EgL_AlGaSb;
            self.EgX[3] = self.moleFrac4*c.EgX_AlSb + (1-self.moleFrac4)*c.EgX_GaSb - self.moleFrac4*(1-self.moleFrac4)*c.EgX_AlGaSb;
            self.VBO[3] = self.moleFrac4*c.VBO_AlSb + (1-self.moleFrac4)*c.VBO_GaSb - self.moleFrac4*(1-self.moleFrac4)*c.VBO_AlGaSb;
            self.DSO[3] = self.moleFrac4*c.DSO_AlSb + (1-self.moleFrac4)*c.DSO_GaSb - self.moleFrac4*(1-self.moleFrac4)*c.DSO_AlGaSb;
            self.me0[3] = self.moleFrac4*c.me0_AlSb + (1-self.moleFrac4)*c.me0_GaSb - self.moleFrac4*(1-self.moleFrac4)*c.me0_AlGaSb;
            self.acG[3] = self.moleFrac4*c.acG_AlSb + (1-self.moleFrac4)*c.acG_GaSb - self.moleFrac4*(1-self.moleFrac4)*c.acG_AlGaSb;
            self.acL[3] = self.moleFrac4*c.acL_AlSb + (1-self.moleFrac4)*c.acL_GaSb - self.moleFrac4*(1-self.moleFrac4)*c.acL_AlGaSb;
            self.acX[3] = self.moleFrac4*c.acX_AlSb + (1-self.moleFrac4)*c.acX_GaSb - self.moleFrac4*(1-self.moleFrac4)*c.acX_AlGaSb;
            self.Ep[3]  = self.moleFrac4*c.Ep_AlSb  + (1-self.moleFrac4)*c.Ep_GaSb  - self.moleFrac4*(1-self.moleFrac4)*c.Ep_AlGaSb;
            self.F[3]   = self.moleFrac4*c.F_AlSb   + (1-self.moleFrac4)*c.F_GaSb   - self.moleFrac4*(1-self.moleFrac4)*c.F_AlGaSb;
            self.XiX[3] = self.moleFrac4*c.XiX_AlSb + (1-self.moleFrac4)*c.XiX_GaSb;
            self.b[3]   = self.moleFrac4*c.b_AlSb   + (1-self.moleFrac4)*c.b_GaSb;
            self.av[3]  = self.moleFrac4*c.av_AlSb  + (1-self.moleFrac4)*c.av_GaSb;
            self.alG[3] = self.moleFrac4*c.alG_AlSb + (1-self.moleFrac4)*c.alG_GaSb;
            self.beG[3] = self.moleFrac4*c.beG_AlSb + (1-self.moleFrac4)*c.beG_GaSb;
            self.alL[3] = self.moleFrac4*c.alL_AlSb + (1-self.moleFrac4)*c.alL_GaSb;
            self.beL[3] = self.moleFrac4*c.beL_AlSb + (1-self.moleFrac4)*c.beL_GaSb;
            self.alX[3] = self.moleFrac4*c.alX_AlSb + (1-self.moleFrac4)*c.alX_GaSb;
            self.beX[3] = self.moleFrac4*c.beX_AlSb + (1-self.moleFrac4)*c.beX_GaSb;
            self.epss[3]   = self.moleFrac4*c.epss_AlSb   + (1-self.moleFrac4)*c.epss_GaSb;
            self.epsInf[3] = self.moleFrac4*c.epsInf_AlSb + (1-self.moleFrac4)*c.epsInf_GaSb;
            self.hwLO[3]   = self.moleFrac4*c.hwLO_AlSb   + (1-self.moleFrac4)*c.hwLO_GaSb;
            self.a0[3]     = self.moleFrac4*c.alc_AlSb    + (1-self.moleFrac4)*c.alc_GaSb;
            self.c11[3]    = self.moleFrac4*c.c11_AlSb    + (1-self.moleFrac4)*c.c11_GaSb;
            self.c12[3]    = self.moleFrac4*c.c12_AlSb    + (1-self.moleFrac4)*c.c12_GaSb;
            
            self.EgG[4] = self.moleFrac5*c.EgG_InAs + (1-self.moleFrac5)*c.EgG_InSb - self.moleFrac5*(1-self.moleFrac5)*c.EgG_InAsSb;
            self.EgL[4] = self.moleFrac5*c.EgL_InAs + (1-self.moleFrac5)*c.EgL_InSb - self.moleFrac5*(1-self.moleFrac5)*c.EgL_InAsSb;
            self.EgX[4] = self.moleFrac5*c.EgX_InAs + (1-self.moleFrac5)*c.EgX_InSb - self.moleFrac5*(1-self.moleFrac5)*c.EgX_InAsSb;
            self.VBO[4] = self.moleFrac5*c.VBO_InAs + (1-self.moleFrac5)*c.VBO_InSb - self.moleFrac5*(1-self.moleFrac5)*c.VBO_InAsSb;
            self.DSO[4] = self.moleFrac5*c.DSO_InAs + (1-self.moleFrac5)*c.DSO_InSb - self.moleFrac5*(1-self.moleFrac5)*c.DSO_InAsSb;
            self.me0[4] = self.moleFrac5*c.me0_InAs + (1-self.moleFrac5)*c.me0_InSb - self.moleFrac5*(1-self.moleFrac5)*c.me0_InAsSb;
            self.acG[4] = self.moleFrac5*c.acG_InAs + (1-self.moleFrac5)*c.acG_InSb - self.moleFrac5*(1-self.moleFrac5)*c.acG_InAsSb;
            self.acL[4] = self.moleFrac5*c.acL_InAs + (1-self.moleFrac5)*c.acL_InSb - self.moleFrac5*(1-self.moleFrac5)*c.acL_InAsSb;
            self.acX[4] = self.moleFrac5*c.acX_InAs + (1-self.moleFrac5)*c.acX_InSb - self.moleFrac5*(1-self.moleFrac5)*c.acX_InAsSb;
            self.Ep[4]  = self.moleFrac5*c.Ep_InAs  + (1-self.moleFrac5)*c.Ep_InSb  - self.moleFrac5*(1-self.moleFrac5)*c.Ep_InAsSb;
            self.F[4]   = self.moleFrac5*c.F_InAs   + (1-self.moleFrac5)*c.F_InSb   - self.moleFrac5*(1-self.moleFrac5)*c.F_InAsSb;
            self.XiX[4] = self.moleFrac5*c.XiX_InAs + (1-self.moleFrac5)*c.XiX_InSb;
            self.b[4]   = self.moleFrac5*c.b_InAs   + (1-self.moleFrac5)*c.b_InSb;
            self.av[4]  = self.moleFrac5*c.av_InAs  + (1-self.moleFrac5)*c.av_InSb;
            self.alG[4] = self.moleFrac5*c.alG_InAs + (1-self.moleFrac5)*c.alG_InSb;
            self.beG[4] = self.moleFrac5*c.beG_InAs + (1-self.moleFrac5)*c.beG_InSb;
            self.alL[4] = self.moleFrac5*c.alL_InAs + (1-self.moleFrac5)*c.alL_InSb;
            self.beL[4] = self.moleFrac5*c.beL_InAs + (1-self.moleFrac5)*c.beL_InSb;
            self.alX[4] = self.moleFrac5*c.alX_InAs + (1-self.moleFrac5)*c.alX_InSb;
            self.beX[4] = self.moleFrac5*c.beX_InAs + (1-self.moleFrac5)*c.beX_InSb;
            self.epss[4]   = self.moleFrac5*c.epss_InAs   + (1-self.moleFrac5)*c.epss_InSb;
            self.epsInf[4] = self.moleFrac5*c.epsInf_InAs + (1-self.moleFrac5)*c.epsInf_InSb;
            self.hwLO[4]   = self.moleFrac5*c.hwLO_InAs   + (1-self.moleFrac5)*c.hwLO_InSb;
            self.a0[4]     = self.moleFrac5*c.alc_InAs    + (1-self.moleFrac5)*c.alc_InSb;
            self.c11[4]    = self.moleFrac5*c.c11_InAs    + (1-self.moleFrac5)*c.c11_InSb;
            self.c12[4]    = self.moleFrac5*c.c12_InAs    + (1-self.moleFrac5)*c.c12_InSb;
            
            EgG_AlGaSb = c.EgG_AlGaSb + 1.22*self.moleFrac6
            self.EgG[5] = self.moleFrac6*c.EgG_AlSb + (1-self.moleFrac6)*c.EgG_GaSb - self.moleFrac6*(1-self.moleFrac6)*EgG_AlGaSb;
            self.EgL[5] = self.moleFrac6*c.EgL_AlSb + (1-self.moleFrac6)*c.EgL_GaSb - self.moleFrac6*(1-self.moleFrac6)*c.EgL_AlGaSb;
            self.EgX[5] = self.moleFrac6*c.EgX_AlSb + (1-self.moleFrac6)*c.EgX_GaSb - self.moleFrac6*(1-self.moleFrac6)*c.EgX_AlGaSb;
            self.VBO[5] = self.moleFrac6*c.VBO_AlSb + (1-self.moleFrac6)*c.VBO_GaSb - self.moleFrac6*(1-self.moleFrac6)*c.VBO_AlGaSb;
            self.DSO[5] = self.moleFrac6*c.DSO_AlSb + (1-self.moleFrac6)*c.DSO_GaSb - self.moleFrac6*(1-self.moleFrac6)*c.DSO_AlGaSb;
            self.me0[5] = self.moleFrac6*c.me0_AlSb + (1-self.moleFrac6)*c.me0_GaSb - self.moleFrac6*(1-self.moleFrac6)*c.me0_AlGaSb;
            self.acG[5] = self.moleFrac6*c.acG_AlSb + (1-self.moleFrac6)*c.acG_GaSb - self.moleFrac6*(1-self.moleFrac6)*c.acG_AlGaSb;
            self.acL[5] = self.moleFrac6*c.acL_AlSb + (1-self.moleFrac6)*c.acL_GaSb - self.moleFrac6*(1-self.moleFrac6)*c.acL_AlGaSb;
            self.acX[5] = self.moleFrac6*c.acX_AlSb + (1-self.moleFrac6)*c.acX_GaSb - self.moleFrac6*(1-self.moleFrac6)*c.acX_AlGaSb;
            self.Ep[5]  = self.moleFrac6*c.Ep_AlSb  + (1-self.moleFrac6)*c.Ep_GaSb  - self.moleFrac6*(1-self.moleFrac6)*c.Ep_AlGaSb;
            self.F[5]   = self.moleFrac6*c.F_AlSb   + (1-self.moleFrac6)*c.F_GaSb   - self.moleFrac6*(1-self.moleFrac6)*c.F_AlGaSb;
            self.XiX[5] = self.moleFrac6*c.XiX_AlSb + (1-self.moleFrac6)*c.XiX_GaSb;
            self.b[5]   = self.moleFrac6*c.b_AlSb   + (1-self.moleFrac6)*c.b_GaSb;
            self.av[5]  = self.moleFrac6*c.av_AlSb  + (1-self.moleFrac6)*c.av_GaSb;
            self.alG[5] = self.moleFrac6*c.alG_AlSb + (1-self.moleFrac6)*c.alG_GaSb;
            self.beG[5] = self.moleFrac6*c.beG_AlSb + (1-self.moleFrac6)*c.beG_GaSb;
            self.alL[5] = self.moleFrac6*c.alL_AlSb + (1-self.moleFrac6)*c.alL_GaSb;
            self.beL[5] = self.moleFrac6*c.beL_AlSb + (1-self.moleFrac6)*c.beL_GaSb;
            self.alX[5] = self.moleFrac6*c.alX_AlSb + (1-self.moleFrac6)*c.alX_GaSb;
            self.beX[5] = self.moleFrac6*c.beX_AlSb + (1-self.moleFrac6)*c.beX_GaSb;
            self.epss[5]   = self.moleFrac6*c.epss_AlSb   + (1-self.moleFrac6)*c.epss_GaSb;
            self.epsInf[5] = self.moleFrac6*c.epsInf_AlSb + (1-self.moleFrac6)*c.epsInf_GaSb;
            self.hwLO[5]   = self.moleFrac6*c.hwLO_AlSb   + (1-self.moleFrac6)*c.hwLO_GaSb;
            self.a0[5]     = self.moleFrac6*c.alc_AlSb    + (1-self.moleFrac6)*c.alc_GaSb;
            self.c11[5]    = self.moleFrac6*c.c11_AlSb    + (1-self.moleFrac6)*c.c11_GaSb;
            self.c12[5]    = self.moleFrac6*c.c12_AlSb    + (1-self.moleFrac6)*c.c12_GaSb;
                
            self.EgG[6] = self.moleFrac7*c.EgG_InAs + (1-self.moleFrac7)*c.EgG_InSb - self.moleFrac7*(1-self.moleFrac7)*c.EgG_InAsSb;
            self.EgL[6] = self.moleFrac7*c.EgL_InAs + (1-self.moleFrac7)*c.EgL_InSb - self.moleFrac7*(1-self.moleFrac7)*c.EgL_InAsSb;
            self.EgX[6] = self.moleFrac7*c.EgX_InAs + (1-self.moleFrac7)*c.EgX_InSb - self.moleFrac7*(1-self.moleFrac7)*c.EgX_InAsSb;
            self.VBO[6] = self.moleFrac7*c.VBO_InAs + (1-self.moleFrac7)*c.VBO_InSb - self.moleFrac7*(1-self.moleFrac7)*c.VBO_InAsSb;
            self.DSO[6] = self.moleFrac7*c.DSO_InAs + (1-self.moleFrac7)*c.DSO_InSb - self.moleFrac7*(1-self.moleFrac7)*c.DSO_InAsSb;
            self.me0[6] = self.moleFrac7*c.me0_InAs + (1-self.moleFrac7)*c.me0_InSb - self.moleFrac7*(1-self.moleFrac7)*c.me0_InAsSb;
            self.acG[6] = self.moleFrac7*c.acG_InAs + (1-self.moleFrac7)*c.acG_InSb - self.moleFrac7*(1-self.moleFrac7)*c.acG_InAsSb;
            self.acL[6] = self.moleFrac7*c.acL_InAs + (1-self.moleFrac7)*c.acL_InSb - self.moleFrac7*(1-self.moleFrac7)*c.acL_InAsSb;
            self.acX[6] = self.moleFrac7*c.acX_InAs + (1-self.moleFrac7)*c.acX_InSb - self.moleFrac7*(1-self.moleFrac7)*c.acX_InAsSb;
            self.Ep[6]  = self.moleFrac7*c.Ep_InAs  + (1-self.moleFrac7)*c.Ep_InSb  - self.moleFrac7*(1-self.moleFrac7)*c.Ep_InAsSb;
            self.F[6]   = self.moleFrac7*c.F_InAs   + (1-self.moleFrac7)*c.F_InSb   - self.moleFrac7*(1-self.moleFrac7)*c.F_InAsSb;
            self.XiX[6] = self.moleFrac7*c.XiX_InAs + (1-self.moleFrac7)*c.XiX_InSb;
            self.b[6]   = self.moleFrac7*c.b_InAs   + (1-self.moleFrac7)*c.b_InSb;
            self.av[6]  = self.moleFrac7*c.av_InAs  + (1-self.moleFrac7)*c.av_InSb;
            self.alG[6] = self.moleFrac7*c.alG_InAs + (1-self.moleFrac7)*c.alG_InSb;
            self.beG[6] = self.moleFrac7*c.beG_InAs + (1-self.moleFrac7)*c.beG_InSb;
            self.alL[6] = self.moleFrac7*c.alL_InAs + (1-self.moleFrac7)*c.alL_InSb;
            self.beL[6] = self.moleFrac7*c.beL_InAs + (1-self.moleFrac7)*c.beL_InSb;
            self.alX[6] = self.moleFrac7*c.alX_InAs + (1-self.moleFrac7)*c.alX_InSb;
            self.beX[6] = self.moleFrac7*c.beX_InAs + (1-self.moleFrac7)*c.beX_InSb;
            self.epss[6]   = self.moleFrac7*c.epss_InAs   + (1-self.moleFrac7)*c.epss_InSb;
            self.epsInf[6] = self.moleFrac7*c.epsInf_InAs + (1-self.moleFrac7)*c.epsInf_InSb;
            self.hwLO[6]   = self.moleFrac7*c.hwLO_InAs   + (1-self.moleFrac7)*c.hwLO_InSb;
            self.a0[6]     = self.moleFrac7*c.alc_InAs    + (1-self.moleFrac7)*c.alc_InSb;
            self.c11[6]    = self.moleFrac7*c.c11_InAs    + (1-self.moleFrac7)*c.c11_InSb;
            self.c12[6]    = self.moleFrac7*c.c12_InAs    + (1-self.moleFrac7)*c.c12_InSb;
            
            EgG_AlGaSb = c.EgG_AlGaSb + 1.22*self.moleFrac8
            self.EgG[7] = self.moleFrac8*c.EgG_AlSb + (1-self.moleFrac8)*c.EgG_GaSb - self.moleFrac8*(1-self.moleFrac8)*EgG_AlGaSb;
            self.EgL[7] = self.moleFrac8*c.EgL_AlSb + (1-self.moleFrac8)*c.EgL_GaSb - self.moleFrac8*(1-self.moleFrac8)*c.EgL_AlGaSb;
            self.EgX[7] = self.moleFrac8*c.EgX_AlSb + (1-self.moleFrac8)*c.EgX_GaSb - self.moleFrac8*(1-self.moleFrac8)*c.EgX_AlGaSb;
            self.VBO[7] = self.moleFrac8*c.VBO_AlSb + (1-self.moleFrac8)*c.VBO_GaSb - self.moleFrac8*(1-self.moleFrac8)*c.VBO_AlGaSb;
            self.DSO[7] = self.moleFrac8*c.DSO_AlSb + (1-self.moleFrac8)*c.DSO_GaSb - self.moleFrac8*(1-self.moleFrac8)*c.DSO_AlGaSb;
            self.me0[7] = self.moleFrac8*c.me0_AlSb + (1-self.moleFrac8)*c.me0_GaSb - self.moleFrac8*(1-self.moleFrac8)*c.me0_AlGaSb;
            self.acG[7] = self.moleFrac8*c.acG_AlSb + (1-self.moleFrac8)*c.acG_GaSb - self.moleFrac8*(1-self.moleFrac8)*c.acG_AlGaSb;
            self.acL[7] = self.moleFrac8*c.acL_AlSb + (1-self.moleFrac8)*c.acL_GaSb - self.moleFrac8*(1-self.moleFrac8)*c.acL_AlGaSb;
            self.acX[7] = self.moleFrac8*c.acX_AlSb + (1-self.moleFrac8)*c.acX_GaSb - self.moleFrac8*(1-self.moleFrac8)*c.acX_AlGaSb;
            self.Ep[7]  = self.moleFrac8*c.Ep_AlSb  + (1-self.moleFrac8)*c.Ep_GaSb  - self.moleFrac8*(1-self.moleFrac8)*c.Ep_AlGaSb;
            self.F[7]   = self.moleFrac8*c.F_AlSb   + (1-self.moleFrac8)*c.F_GaSb   - self.moleFrac8*(1-self.moleFrac8)*c.F_AlGaSb;
            self.XiX[7] = self.moleFrac8*c.XiX_AlSb + (1-self.moleFrac8)*c.XiX_GaSb;
            self.b[7]   = self.moleFrac8*c.b_AlSb   + (1-self.moleFrac8)*c.b_GaSb;
            self.av[7]  = self.moleFrac8*c.av_AlSb  + (1-self.moleFrac8)*c.av_GaSb;
            self.alG[7] = self.moleFrac8*c.alG_AlSb + (1-self.moleFrac8)*c.alG_GaSb;
            self.beG[7] = self.moleFrac8*c.beG_AlSb + (1-self.moleFrac8)*c.beG_GaSb;
            self.alL[7] = self.moleFrac8*c.alL_AlSb + (1-self.moleFrac8)*c.alL_GaSb;
            self.beL[7] = self.moleFrac8*c.beL_AlSb + (1-self.moleFrac8)*c.beL_GaSb;
            self.alX[7] = self.moleFrac8*c.alX_AlSb + (1-self.moleFrac8)*c.alX_GaSb;
            self.beX[7] = self.moleFrac8*c.beX_AlSb + (1-self.moleFrac8)*c.beX_GaSb;
            self.epss[7]   = self.moleFrac8*c.epss_AlSb   + (1-self.moleFrac8)*c.epss_GaSb;
            self.epsInf[7] = self.moleFrac8*c.epsInf_AlSb + (1-self.moleFrac8)*c.epsInf_GaSb;
            self.hwLO[7]   = self.moleFrac8*c.hwLO_AlSb   + (1-self.moleFrac8)*c.hwLO_GaSb;
            self.a0[7]     = self.moleFrac8*c.alc_AlSb    + (1-self.moleFrac8)*c.alc_GaSb;
            self.c11[7]    = self.moleFrac8*c.c11_AlSb    + (1-self.moleFrac8)*c.c11_GaSb;
            self.c12[7]    = self.moleFrac8*c.c12_AlSb    + (1-self.moleFrac8)*c.c12_GaSb;
            
#        else:
#            raise TypeError('substrate selection not allowed')

        #set this once the others are set        
        self.epsrho = 1 / (1/self.epsInf - 1/self.epss)     

    def update_strain(self):  # c is a Material_Constant class instance

        if self.substrate == 'InP':
            self.a_parallel = c.alc_InP
        elif self.substrate == 'GaSb':
            self.a_parallel = c.alc_GaSb
        elif self.substrate == 'GaAs':
            self.a_parallel = c.alc_GaAs
        else:
            assert 0 == 1 #obviously not ;)
            
        indx = nonzero(self.layerMaterials[1:] == 1)[0]  #[1:] because first layer doesn't count
        self.h[1] = sum(self.layerWidths[indx]*self.layerBarriers[indx])
        self.h[0] = sum(self.layerWidths[indx]) - self.h[1]
        indx = nonzero(self.layerMaterials[1:] == 2)[0]
        self.h[3] = sum(self.layerWidths[indx]*self.layerBarriers[indx])
        self.h[2] = sum(self.layerWidths[indx]) - self.h[3]
        indx = nonzero(self.layerMaterials[1:] == 3)[0]  #[1:] because first layer doesn't count
        self.h[5] = sum(self.layerWidths[indx]*self.layerBarriers[indx])
        self.h[4] = sum(self.layerWidths[indx]) - self.h[5]
        indx = nonzero(self.layerMaterials[1:] == 4)[0]  #[1:] because first layer doesn't count
        self.h[7] = sum(self.layerWidths[indx]*self.layerBarriers[indx])
        self.h[6] = sum(self.layerWidths[indx]) - self.h[7]
        
        self.eps_parallel = self.a_parallel / self.a0 - 1; #Walle eqn 1b
        
        self.a_perp   = self.a0 * (1 - 2* self.c12 / self.c11 * self.eps_parallel); #Walle eqn 2a
        self.eps_perp = self.a_perp/self.a0 - 1; #Walle eqn 2b
        
        self.mismatch = 100 * sum(self.h*self.eps_perp) / sum(self.h);
        
        self.MLThickness = zeros(self.layerMaterials.size)
        self.MLThickness[nonzero((self.layerMaterials == 1) & (self.layerBarriers == 0))[0]] = self.a_perp[0] / 2.0
        self.MLThickness[nonzero((self.layerMaterials == 1) & (self.layerBarriers == 1))[0]] = self.a_perp[1] / 2.0
        self.MLThickness[nonzero((self.layerMaterials == 2) & (self.layerBarriers == 0))[0]] = self.a_perp[2] / 2.0
        self.MLThickness[nonzero((self.layerMaterials == 2) & (self.layerBarriers == 1))[0]] = self.a_perp[3] / 2.0
        self.MLThickness[nonzero((self.layerMaterials == 3) & (self.layerBarriers == 0))[0]] = self.a_perp[4] / 2.0
        self.MLThickness[nonzero((self.layerMaterials == 3) & (self.layerBarriers == 1))[0]] = self.a_perp[5] / 2.0
        self.MLThickness[nonzero((self.layerMaterials == 4) & (self.layerBarriers == 0))[0]] = self.a_perp[6] / 2.0
        self.MLThickness[nonzero((self.layerMaterials == 4) & (self.layerBarriers == 1))[0]] = self.a_perp[7] / 2.0
    
        # for In0.53Ga0.47As, EcG = 0.22004154
        #    use this as a zero point baseline
        baseln = 0.22004154
        
        self.Pec = (2*self.eps_parallel+self.eps_perp) * (self.acG)
        self.Pe = 2*self.av * (self.c11-self.c12) / self.c11 * self.eps_parallel
        self.Qe = - self.b * (self.c11+2*self.c12) / self.c11 * self.eps_parallel
        self.ESO  = sqrt(9*self.Qe**2+2*self.Qe*self.DSO+self.DSO**2);
        self.Varsh = - self.alG*c.Temperature**2/(c.Temperature+self.beG)
        
        #These are wrong
        #self.EcG = self.VBO + self.EgG + self.Pec + self.Pe + self.Varsh - baseln; #This equation is wrong; no self.Pe
        #self.EcL = self.VBO + self.EgL + (2*self.eps_parallel+self.eps_perp) * (self.acL+self.av) - baseln;
        #self.EcX = self.VBO + self.EgX + (2*self.eps_parallel+self.eps_perp)*(self.acX+self.av) + 2/3 * self.XiX * (self.eps_perp-self.eps_parallel) - baseln;
        
        # calculations for effective mass
        #  in part following Sugawara, PRB 48, 8102 (1993)
        
        self.EgLH = self.EgG + self.Pec + self.Pe + self.Varsh - 1/2*(self.Qe - self.DSO + sqrt(9*self.Qe**2+2*self.Qe*self.DSO+self.DSO**2));
        self.EgSO = self.EgG + self.Pec + self.Pe + self.Varsh - 1/2*(self.Qe - self.DSO - sqrt(9*self.Qe**2+2*self.Qe*self.DSO+self.DSO**2));
        self.ESO  = sqrt(9*self.Qe**2+2*self.Qe*self.DSO+self.DSO**2);
        
        # self.alpha = (self.Qe + self.DSO + sqrt(9*self.Qe.^2+2*self.Qe*self.DSO+self.DSO.^2)) / (2*sqrt(9*self.Qe.^2+2*self.Qe*self.DSO+self.DSO.^2));
        # self.beta = (-self.Qe - self.DSO + sqrt(9*self.Qe.^2+2*self.Qe*self.DSO+self.DSO.^2)) / (2*sqrt(9*self.Qe.^2+2*self.Qe*self.DSO+self.DSO.^2));
        # self.alpha = self.alpha / sqrt(self.alpha.^2 + self.beta.^2);
        # self.beta = self.beta / sqrt(self.alpha.^2 + self.beta.^2);
        self.me = 1 / ((1+2*self.F) + self.Ep/self.EgLH*(self.EgLH+2/3*self.ESO)/(self.EgLH + self.ESO));
        # self.me2 = 1 / ((1+2*self.F) + self.Ep/3 *( (sqrt(2)*self.alpha-self.beta).^2 / self.EgLH + (sqrt(2)*self.beta+self.alpha).^2 / (self.EgLH + self.ESO))); #using Igor's Ep instead of the calculated P2
        # self.me3 = 1 / ((1+2*self.F) + self.Ep/3 *( 2 / self.EgLH + 1 / (self.EgLH + self.ESO)));
        # self.me4 = 1 / ((1+2*self.F) + self.Ep/3 *( 2 / self.EgLH + 1 / (self.EgSO)));
        
        #corrections to the method used to calculate band edges, thanks to Yu Song
        self.EcG = self.VBO + self.EgG + self.Pec - baseln;
        self.EcL = self.VBO + self.EgL + (2*self.eps_parallel+self.eps_perp) * (self.acL+self.av) - baseln;
        self.EcX = self.VBO + self.EgX + (2*self.eps_parallel+self.eps_perp)*(self.acX+self.av) + 2/3 * self.XiX * (self.eps_perp-self.eps_parallel) - baseln;
        
        self.EgLH = self.EgG + self.Pec + self.Pe - 1/2*(self.Qe - self.DSO + self.ESO);
        self.EgSO = self.EgG + self.Pec + self.Pe - 1/2*(self.Qe - self.DSO - self.ESO);
        
        #1st MAJOR assumption: Varshney contribution to band edge is in proportion to percent of band offset
        #2nd major assumption: Temperature affects sattelite valleys in the same way it does the Gamma valley
        barrs = array([1,3,5,7])
        wells = array([0,2,4,6])
        CBOffset = self.EcG[barrs] - self.EcG[wells]
        VBOffset = (self.EcG[barrs] - self.EgLH[barrs]) - (self.EcG[wells] - self.EgLH[wells])
        percentCB = CBOffset / (CBOffset + VBOffset)
        percentCB = column_stack([percentCB,percentCB]).flatten() #applies percent CV to both well and barrier slots
        
        self.EcG += percentCB * self.Varsh
        self.EcL += percentCB * self.Varsh
        self.EcX += percentCB * self.Varsh
        self.EvLH = self.EcG - self.EgLH - ((1-percentCB) * self.Varsh)
        self.EvSO = self.EcG - self.EgSO - ((1-percentCB) * self.Varsh)
        
        
        #other code that I've junked for some reason or other
        # self.EgAV = self.EgG - self.alG*c.Temperature.^2/(c.Temperature+self.beG) + (2*self.eps_parallel+self.eps_perp) * (self.acG+self.av);
        # self.EgEG = self.EgAV + 1/3*self.ESO;
        # self.P2 = 1/2 * (1/self.me0 - (1+2*self.F)/1) * self.EgG * (self.EgG+self.DSO) / (self.EgG+2/3*self.DSO); #strain independent
        # self.P2 = 1/2 * (1/self.me0 - (1+2*self.F)/1) * self.EgLH * (self.EgLH+self.ESO) / (self.EgLH+2/3*self.ESO);
        # Pe = 2*self.acG * (self.c11-self.c12) / self.c11 * self.eps_parallel;
        # Elh = -Pe + 1/2 * (self.Qe - self.DSO + sqrt(9*self.Qe.^2+2*self.Qe*self.DSO+self.DSO.^2))
        # Elh = -Pe + 1/2 * (self.Qe - self.DSO + sqrt(9*self.Qe.^2+2*self.Qe*self.DSO+self.DSO.^2));
        # Eso = -Pe + 1/2 * (self.Qe - self.DSO - sqrt(9*self.Qe.^2+2*self.Qe*self.DSO+self.DSO.^2));
        # self.me = 1 / ((1+self.F) + 2*self.P2/3 *( (sqrt(2)*self.alpha-self.beta).^2 / self.EgLH + (sqrt(2)*self.beta+self.alpha).^2 / self.EgSO));
        # self.me3 = 1 / ((2+2*self.F) + self.Ep*(self.EgLH + 2/3*self.ESO) / self.EgLH / (self.EgLH + self.ESO) ); #using Igor's effective mass formula
        # self.me4 = self.EgAV / self.Ep;

    def solve_psi(self):
        Epoints = arange(min(self.xVc),max(self.xVc-115*self.EField*1e-5),self.vertRes/1000)
        xMcE = zeros(self.xPoints.shape)
        xPsi = zeros(self.xPoints.shape)
        psiEnd = zeros(Epoints.size)
        
        #TODO: add adaptive spacing for Eq
        #TODO: convert nested for loop to C
        if True:
            cFunctions.psiFnEnd(Epoints.ctypes.data_as(c_void_p), int(Epoints.size), 
                                int(xPsi.size), c_double(self.xres), c_double(self.EField), self.xVc.ctypes.data_as(c_void_p),
                                self.xEg.ctypes.data_as(c_void_p), self.xF.ctypes.data_as(c_void_p),
                                self.xEp.ctypes.data_as(c_void_p), self.xESO.ctypes.data_as(c_void_p), 
                                self.xMc.ctypes.data_as(c_void_p), xMcE.ctypes.data_as(c_void_p),
                                xPsi.ctypes.data_as(c_void_p), psiEnd.ctypes.data_as(c_void_p))
            #psiFnEnd(double *eEq, int eEqSize, int xPsiSize, double xres, double *xVc, double *xEg, 
            #  double *xF, double *xEp, double *xESO, double *xMcE, double *xPsi, double *xPsiEnd)
            # my_sum.sum(a.ctypes.data_as(c_void_p), int(10))
        else:
            for p, Eq in enumerate(Epoints):
                if True:
                    xMcE = m0 / (1+2*self.xF + self.xEp/3 * (2 / ((Eq-self.xVc)+self.xEg) + 1 / ((Eq-self.xVc)+self.xEg+self.xESO) ))
                else:
                    xMcE = self.xMc * (1 - (self.xVc - Eq) / self.xEg)
                xMcE[0:-1] = 0.5 * (xMcE[0:-1]+xMcE[1:])
                xPsi[0] = 0
                xPsi[1] = 1
                for q in xrange(1,xPsi.size-1):
                    xPsi[q+1] = ((2 * self.xres*1e-10 * self.xres*1e-10 / hbar / hbar * (self.xVc[q] - Eq)*e0 + 1 / xMcE[q] + 1 / xMcE[q-1]) * xPsi[q] - xPsi[q-1] / xMcE[q-1]) * xMcE[q]
                psiEnd[p] = xPsi[-1]
                
        #interpolate between solved-for E points        
        tck = interpolate.splrep(Epoints,psiEnd,s=0)
        xnew = linspace(Epoints[0],Epoints[-1],Epoints.size*1e2) #adds 100 points per solved-for E point
        ynew = interpolate.splev(xnew,tck,der=0)
        
#        #plot interpolated points over solved-for E points
#        plt.figure()
#        plt.plot(Epoints,psiEnd,'x',xnew,ynew,'o-')
#        plt.show()
       
        #find Eigen Energies
        #This routine looks for all of the zero crossings, and then picks each one out
        gtz = ynew > 0
        ltz = ynew < 0
        overlap1 = bitwise_and(gtz[0:-1],ltz[1:])
        overlap2 = bitwise_and(gtz[1:],ltz[0:-1])
        overlap  = bitwise_or(overlap1, overlap2)
        idxs = nonzero(overlap == True)[0]
        self.EigenE = zeros(idxs.size)
        
        if True:
            cFunctions.inv_quadratic_interp(xnew.ctypes.data_as(c_void_p), ynew.ctypes.data_as(c_void_p), 
                                            idxs.ctypes.data_as(c_void_p), int(self.EigenE.size), 
                                            self.EigenE.ctypes.data_as(c_void_p))
        else:
            for q, idx in enumerate(idxs): # do quadratic interpolation
                x0=xnew[idx-1]; fx0=ynew[idx-1]
                x1=xnew[idx];   fx1=ynew[idx]
                x2=xnew[idx+1]; fx2=ynew[idx+1]
                d1=(fx1-fx0)/(x1-x0); d2=(fx2-fx1)/(x2-x1);
                #inverse quadratic interpolation
                x3 = x0*fx1*fx2/(fx0-fx1)/(fx0-fx2) + x1*fx0*fx2/(fx1-fx0)/(fx1-fx2) + x2*fx0*fx1/(fx2-fx0)/(fx2-fx1)
                self.EigenE[q] = x3
#                if abs(d1) > 1e15 and abs(d2) > 1e15:
#                    self.EigenE[q] = 0

        if True:
            for q in xrange(self.EigenE.size):
                x0=self.EigenE[q]-self.vertRes/100000
                x1=self.EigenE[q]
                x2=self.EigenE[q]+self.vertRes/100000
                
                cFunctions.psiFn(c_double(x0), int(1), 
                                int(xPsi.size), c_double(self.xres), self.xVc.ctypes.data_as(c_void_p),
                                self.xEg.ctypes.data_as(c_void_p), self.xF.ctypes.data_as(c_void_p),
                                self.xEp.ctypes.data_as(c_void_p), self.xESO.ctypes.data_as(c_void_p), 
                                self.xMc.ctypes.data_as(c_void_p), xMcE.ctypes.data_as(c_void_p),
                                xPsi.ctypes.data_as(c_void_p))
                fx0 = xPsi[-1]
                
                cFunctions.psiFn(c_double(x1), int(1), 
                                int(xPsi.size), c_double(self.xres), self.xVc.ctypes.data_as(c_void_p),
                                self.xEg.ctypes.data_as(c_void_p), self.xF.ctypes.data_as(c_void_p),
                                self.xEp.ctypes.data_as(c_void_p), self.xESO.ctypes.data_as(c_void_p), 
                                self.xMc.ctypes.data_as(c_void_p), xMcE.ctypes.data_as(c_void_p),
                                xPsi.ctypes.data_as(c_void_p))
                fx1 = xPsi[-1]
                
                cFunctions.psiFn(c_double(x2), int(1), 
                                int(xPsi.size), c_double(self.xres), self.xVc.ctypes.data_as(c_void_p),
                                self.xEg.ctypes.data_as(c_void_p), self.xF.ctypes.data_as(c_void_p),
                                self.xEp.ctypes.data_as(c_void_p), self.xESO.ctypes.data_as(c_void_p), 
                                self.xMc.ctypes.data_as(c_void_p), xMcE.ctypes.data_as(c_void_p),
                                xPsi.ctypes.data_as(c_void_p))
                fx2 = xPsi[-1]
                
                #psiFn(double Eq, int startpoint, int xPsiSize, double xres, double *xVc, double *xEg, 
                #      double *xF, double *xEp, double *xESO, double *xMc, double *xMcE, double *xPsi)
                
               
                d1=(fx1-fx0)/(x1-x0); d2=(fx2-fx1)/(x2-x1);
                #inverse quadratic interpolation
                x3 = x0*fx1*fx2/(fx0-fx1)/(fx0-fx2) + x1*fx0*fx2/(fx1-fx0)/(fx1-fx2) + x2*fx0*fx1/(fx2-fx0)/(fx2-fx1)
                self.EigenE[q] = x3
                
                
        if False:
            for q in xrange(self.EigenE.size):
                x0=self.EigenE[q]-self.vertRes/1000000
                x1=self.EigenE[q]
                x2=self.EigenE[q]+self.vertRes/1000000
                
                cFunctions.psiFn(c_double(x0), int(1), 
                                int(xPsi.size), c_double(self.xres), self.xVc.ctypes.data_as(c_void_p),
                                self.xEg.ctypes.data_as(c_void_p), self.xF.ctypes.data_as(c_void_p),
                                self.xEp.ctypes.data_as(c_void_p), self.xESO.ctypes.data_as(c_void_p), 
                                self.xMc.ctypes.data_as(c_void_p), xMcE.ctypes.data_as(c_void_p),
                                xPsi.ctypes.data_as(c_void_p))
                fx0 = xPsi[-1]
                
                cFunctions.psiFn(c_double(x1), int(1), 
                                int(xPsi.size), c_double(self.xres), self.xVc.ctypes.data_as(c_void_p),
                                self.xEg.ctypes.data_as(c_void_p), self.xF.ctypes.data_as(c_void_p),
                                self.xEp.ctypes.data_as(c_void_p), self.xESO.ctypes.data_as(c_void_p), 
                                self.xMc.ctypes.data_as(c_void_p), xMcE.ctypes.data_as(c_void_p),
                                xPsi.ctypes.data_as(c_void_p))
                fx1 = xPsi[-1]
                
                cFunctions.psiFn(c_double(x2), int(1), 
                                int(xPsi.size), c_double(self.xres), self.xVc.ctypes.data_as(c_void_p),
                                self.xEg.ctypes.data_as(c_void_p), self.xF.ctypes.data_as(c_void_p),
                                self.xEp.ctypes.data_as(c_void_p), self.xESO.ctypes.data_as(c_void_p), 
                                self.xMc.ctypes.data_as(c_void_p), xMcE.ctypes.data_as(c_void_p),
                                xPsi.ctypes.data_as(c_void_p))
                fx2 = xPsi[-1]
                
                #psiFn(double Eq, int startpoint, int xPsiSize, double xres, double *xVc, double *xEg, 
                #      double *xF, double *xEp, double *xESO, double *xMc, double *xMcE, double *xPsi)
                
               
                d1=(fx1-fx0)/(x1-x0); d2=(fx2-fx1)/(x2-x1);
                #inverse quadratic interpolation
                x3 = x0*fx1*fx2/(fx0-fx1)/(fx0-fx2) + x1*fx0*fx2/(fx1-fx0)/(fx1-fx2) + x2*fx0*fx1/(fx2-fx0)/(fx2-fx1)
                self.EigenE[q] = x3

        #make array for Psi and fill it in
        if True:
            self.xyPsi = zeros(self.xPoints.size*self.EigenE.size)
            cFunctions.psiFill(int(xPsi.size), c_double(self.xres),
                               int(self.EigenE.size), self.EigenE.ctypes.data_as(c_void_p),
                               self.xVc.ctypes.data_as(c_void_p),
                               self.xEg.ctypes.data_as(c_void_p), self.xF.ctypes.data_as(c_void_p),
                               self.xEp.ctypes.data_as(c_void_p), self.xESO.ctypes.data_as(c_void_p), 
                               self.xMc.ctypes.data_as(c_void_p), xMcE.ctypes.data_as(c_void_p),
                               self.xyPsi.ctypes.data_as(c_void_p))
            #self.xyPsi.resize(a.xPoints.size, a.EigenE.size)
            self.xyPsi = self.xyPsi.reshape(self.xPoints.size, self.EigenE.size, order='F')
        else:
            self.xyPsi = zeros((self.xPoints.size,self.EigenE.size))
            for p, Eq in enumerate(self.EigenE):
                if True:
                    xMcE = m0 / (1+2*self.xF + self.xEp/3 * (2 / ((Eq-self.xVc)+self.xEg) + 1 / ((Eq-self.xVc)+self.xEg+self.xESO) ))
                else:
                    xMcE = self.xMc * (1 - (self.xVc - Eq) / self.xEg)
                xMcE[0:-1] = 0.5 * (xMcE[0:-1]+xMcE[1:])
                xPsi[1] = 1
                for q in xrange(2,xPsi.size-1):
                    xPsi[q+1] = ((2 * self.xres*1e-10 * self.xres*1e-10 / hbar / hbar * (self.xVc[q] - Eq)*e0 + 1 / xMcE[q] + 1 / xMcE[q-1]) * xPsi[q] - xPsi[q-1] / xMcE[q-1]) * xMcE[q]
                psiInt = sum(xPsi**2 * (1+(Eq-self.xVc)/(Eq-self.xVc+self.xEg)))
                A = 1 / sqrt( self.xres * 1e-10 * psiInt)
                self.xyPsi[:,p] = A * xPsi
        
        #remove states that come from oscillating end points
        psiEnd = self.xyPsi[-1,:]
        idxs = nonzero(abs(psiEnd)<10)[0]
        self.EigenE = self.EigenE[idxs]
        self.xyPsi = self.xyPsi[:,idxs]

        self.xyPsiPsi = self.xyPsi*self.xyPsi*settings.wf_scale #4.5e-10 scales size of wavefunctions, arbitrary for nice plots

        #remove states that are smaller than minimum height
        # addresses states high above band edge
        idxs = nonzero(self.xyPsiPsi.max(0) > settings.wf_scale * settings.wf_min_height)[0] #0.014 is arbitrary; changes if 4.5e-10 changes
        self.EigenE = self.EigenE[idxs]
        self.xyPsi = self.xyPsi[:,idxs]
        self.xyPsiPsi = self.xyPsiPsi[:,idxs]
        
        self.xyPsiPsi2 = copy.deepcopy(self.xyPsiPsi)

        #implement pretty plot
        for q in xrange(self.EigenE.size):
            prettyIdxs = nonzero(self.xyPsiPsi[:,q] > settings.wf_scale * settings.pretty_plot_factor)[0] #0.0005 is arbitrary
            self.xyPsiPsi[0:prettyIdxs[0],q] = NaN
            self.xyPsiPsi[prettyIdxs[-1]:,q] = NaN

        #decimate plot points
        idxs = arange(0,self.xPoints.size, int(settings.plot_decimate_factor/self.xres), dtype=int)
        self.xyPsiPsiDec = zeros([idxs.size, self.EigenE.size])
        for q in xrange(self.EigenE.size):
            self.xyPsiPsiDec[:,q] = self.xyPsiPsi[idxs,q]
        self.xyPsiPsi = self.xyPsiPsiDec
        self.xPointsPost = self.xPoints[idxs]


def basisSolve(data):
    #self.data.basisInjectorAR is 0-to-1
    #self.data.basisARInjector is 1-to-0
        
    #find all 0-to-1 & 1-to-0 transitions
    zeroTOone = []
    oneTOzero = []
    layerAR = hstack([data.layerARs[-1], data.layerARs])
    for q in xrange(0,layerAR.size-1):
        if layerAR[q] == 0 and layerAR[q+1] == 1:
            zeroTOone.append(q-1)
        if layerAR[q] == 1 and layerAR[q+1] == 0:
            oneTOzero.append(q)

    dividers = [0, data.layerARs.size-1]
    if data.basisInjectorAR:
        dividers += zeroTOone
    if data.basisARInjector:
        dividers += oneTOzero
    dividers = list(set(dividers)) #This converts the list into a set, thereby removing duplicates, and then back into a list.
    dividers.sort()
    
    #the first region is always the "wrap around" region
    layers = []
    for q in xrange(len(dividers)-1):
        layers.append(range(dividers[q],dividers[q+1]+1))
    
    dCL = [] #this is dataClassesList. it holds all of the Data classes for each individual solve section

    #for first period only
    #this handles all of the solving        
    for n in xrange(len(layers)):
        dCL.append(copy.deepcopy(data))
        dCL[n].repeats = 1
        
        #substitute proper layer characteristics into dCL[n]
        dCL[n].layerWidths = hstack([100, data.layerWidths[layers[n]], 30])
        dCL[n].layerBarriers = hstack([1, data.layerBarriers[layers[n]], 1])
        dCL[n].layerARs = hstack([0, data.layerARs[layers[n]], 0])
        dCL[n].layerMaterials = hstack([data.layerMaterials[layers[n]][0], data.layerMaterials[layers[n]], data.layerMaterials[layers[n]][-1]])
        dCL[n].layerDopings = hstack([0, data.layerDopings[layers[n]], 0])
        dCL[n].layerDividers = hstack([0, data.layerDividers[layers[n]], 0])
       
        #update and solve
        dCL[n].update_alloys()
        dCL[n].update_strain()
        dCL[n].populate_x()
        dCL[n].populate_x_full()
        dCL[n].solve_psi()
        
        #caculate offsets
        dCL[n].widthOffset = sum(data.layerWidths[range(0,dividers[n])]) #- 100/data.xres
        dCL[n].fieldOffset = -(dCL[n].widthOffset-100) * dCL[n].EField * 1e-5            
    
    solvePeriods = len(dCL)
    
    #create dCL's and offsets for repeat periods
    counter = solvePeriods
    if data.repeats > 1:
        for q in xrange(1,data.repeats):
            for p in xrange(0,solvePeriods):
                dCL.append(copy.deepcopy(dCL[p]))
                dCL[counter].widthOffset = sum(data.layerWidths[1:])*q + dCL[p].widthOffset #- 100/data.xres
                dCL[counter].fieldOffset = -(dCL[counter].widthOffset-100) * dCL[counter].EField * 1e-5
                counter += 1
    return dCL

def convert_dCL_to_data(data, dCL):
    #count number of wavefunctions
    numWFs = 0
    for q in xrange(len(dCL)):
        numWFs += dCL[q].EigenE.size

    #convert to int to prevent machine rounding errors
    data.xPointsPost = arange(-100,data.xPoints[-1]+30+data.xres,data.xres)

    data.xyPsi = zeros([data.xPointsPost.size, numWFs])
    data.xyPsiPsi = NaN*zeros([data.xPointsPost.size, numWFs])
    data.EigenE = zeros(numWFs)
    data.moduleID = zeros(numWFs)
    counter = 0
    for n in xrange(len(dCL)):
        for q in xrange(dCL[n].EigenE.size):
            wf = NaN*zeros(data.xPointsPost.size)
            begin = int(dCL[n].widthOffset/data.xres)
            end = begin + dCL[n].xyPsiPsi2[:,q].size
            wf[begin:end] = dCL[n].xyPsiPsi2[:,q]
            data.xyPsiPsi[:,counter] = wf
            data.EigenE[counter] = dCL[n].EigenE[q] + dCL[n].fieldOffset
            
            data.moduleID[counter] = n
            wf = zeros(data.xPointsPost.size)
            begin = int(dCL[n].widthOffset/data.xres)
            end = begin + dCL[n].xyPsi[:,q].size
            wf[begin:end] = dCL[n].xyPsi[:,q]
            data.xyPsi[:,counter] = wf
            data.xyPsiPsi[:,counter] = wf**2 * settings.wf_scale
            
            counter += 1
    data.xPointsPost = data.xPointsPost[int(100/data.xres):int(-30/data.xres)]
    data.xyPsi = data.xyPsi[int(100/data.xres):int(-30/data.xres)]
    data.xyPsiPsi = data.xyPsiPsi[int(100/data.xres):int(-30/data.xres)]
    
    #implement pretty plot
    for q in xrange(data.EigenE.size):
        prettyIdxs = nonzero(data.xyPsiPsi[:,q] > settings.wf_scale * settings.pretty_plot_factor)[0] #0.0005 is arbitrary
        data.xyPsiPsi[0:prettyIdxs[0],q] = NaN
        data.xyPsiPsi[prettyIdxs[-1]:,q] = NaN
        
    #sort by ascending energy
    sortID = argsort(data.EigenE)
    data.EigenE = data.EigenE[sortID]
    data.xyPsi = data.xyPsi[:,sortID]
    data.xyPsiPsi = data.xyPsiPsi[:,sortID]
    data.moduleID = data.moduleID[sortID]

#    #decimate plot points
#    idxs = arange(0,data.xPoints.size, int(settings.plot_decimate_factor/data.xres), dtype=int)
#    data.xyPsiPsiDec = zeros([idxs.size, data.EigenE.size])
#    for q in xrange(data.EigenE.size):
#        data.xyPsiPsiDec[:,q] = data.xyPsiPsi[idxs,q]
#    data.xyPsiPsi = data.xyPsiPsiDec
#    data.xPointsPost = data.xPoints[idxs]

def lo_phonon_time(data, upper, lower):
    if upper < lower:
        temp = upper
        upper = lower
        lower = temp
        
    psi_i = data.xyPsi[:,upper]
    psi_j = data.xyPsi[:,lower]
    E_i = data.EigenE[upper]
    E_j = data.EigenE[lower]

    if E_i-E_j-data.hwLO[0] < 0:
        return 1e20
        
    idxs_i = nonzero(psi_i > settings.wf_scale * settings.phonon_integral_factor)[0]
    idxs_j = nonzero(psi_j > settings.wf_scale * settings.phonon_integral_factor)[0]
    idx_first = [idxs_i[0], idxs_j[0]]
    idx_last  = [idxs_i[-1], idxs_j[-1]]
    if max(idx_first) > min(idx_last):
        return 1e20
    else:
        idx_first = min([idxs_i[0], idxs_j[0]])
        idx_last  = max([idxs_i[-1], idxs_j[-1]])
        psi_i = psi_i[idx_first:idx_last]
        psi_j = psi_j[idx_first:idx_last]
        xPoints = data.xPoints[idx_first:idx_last]

        xMcE_j = data.xMc * (1 - (data.xVc - E_j) / data.xEg)        
        McE_j = m0*sum(xMcE_j[idx_first:idx_last] * psi_j**2) / sum(psi_j**2) #weight non-parabolic effective mass by probability density
        
        kl = sqrt(2*McE_j/hbar**2 * (E_i-E_j-data.hwLO[0])*e0)
        dIij = zeros(xPoints.size)
        for n in xrange(xPoints.size):
            x1 = xPoints[n]*1e-10
            x2 = xPoints*1e-10
            dIij[n] = sum(psi_i*psi_j * exp(-kl*abs(x1-x2)) * psi_i[n]*psi_j[n] * (data.xres*1e-10)**2)
        Iij = sum(dIij)
        inverse_tau = McE_j * e0**2 * data.hwLO[0]*e0/hbar * Iij / (4 * hbar**2 * data.epsrho[0]*eps0 * kl)
        tau = 1e12/inverse_tau
        return tau

def dipole(data, upper, lower):
    if upper < lower:
        temp = upper
        upper = lower
        lower = temp
    psi_i = data.xyPsi[:,upper]
    psi_j = data.xyPsi[:,lower]
    E_i = data.EigenE[upper]
    E_j = data.EigenE[lower]
    
    data.populate_x_full()
    xMcE_i = data.xMc * (1 - (data.xVc - E_i) / data.xEg)
    xMcE_j = data.xMc * (1 - (data.xVc - E_j) / data.xEg)
    xMcE_j_avg = 0.5 * (xMcE_j[0:-1]+xMcE_j[1:])
    psi_i_avg = 0.5 * (psi_i[0:-1]+psi_i[1:])
    z = sum(psi_i_avg * diff(psi_j/xMcE_i) + 1/xMcE_j_avg * (psi_i_avg * diff(psi_j)))
    z *= hbar**2/(2*(E_i-E_j)*e0) * 1e10/m0
    return z

def coupling_energy(data, dCL, upper, lower):
    #here, psi_i is the left-most wavefunction, not the wf with the highest energy
    module_i = data.moduleID[upper]
    module_j = data.moduleID[lower]
    psi_i = data.xyPsi[:,upper]
    psi_j = data.xyPsi[:,lower]
    if module_i > module_j:
        module_i = data.moduleID[lower]
        module_j = data.moduleID[upper]
        psi_i = data.xyPsi[:,lower]
        psi_j = data.xyPsi[:,upper]

    DeltaV = ones(data.xPointsPost.size)
    first = dCL[int(module_i)].widthOffset/data.xres
    last = first + dCL[int(module_i)].xBarriers[100/data.xres:].size
    DeltaV[first:last] = dCL[int(module_i)].xBarriers[100/data.xres:]
    DeltaV = 1 - DeltaV
    couplingEnergy = sum(psi_i * DeltaV * psi_j) * data.xres*1e-10 * abs(data.EcG[1] - data.EcG[0]) * 1000 #* (1-data.xBarriers)
    return couplingEnergy

def broadening_energy(data, upper, lower):
    if upper < lower:
        temp = upper
        upper = lower
        lower = temp
    psi_i = data.xyPsi[:,upper]
    psi_j = data.xyPsi[:,lower]
    
    transitions = bitwise_xor(data.xBarriers[0:-1].astype(bool),data.xBarriers[1:].astype(bool))
    transitionIdxs = nonzero(transitions == True)[0]# + int(100/data.xres)
    psi2int2 = sum((psi_i[transitionIdxs]**2-psi_j[transitionIdxs]**2)**2)
    DeltaLambda = 0.76 * 1e-9
    twogamma = pi*data.me[0]*m0*e0**2/hbar**2 * DeltaLambda**2 * (data.EcG[1] - data.EcG[0])**2 * psi2int2 *1e3
    return twogamma

def alphaISB(data, stateR, lower):
    statesQ = []
    dipoles = []
    gammas = []
    energies = []
    for q in xrange(stateR+1, data.EigenE.size):
        dp = dipole(data, q, stateR)
        if abs(dp) > 1e-6:
            statesQ.append(q)
            dipoles.append(dp)
    for q in statesQ:
        gamma = broadening_energy(data, q, stateR)/2
        gammas.append(gamma)
        energies.append(data.EigenE[q]-data.EigenE[stateR])
        
    dipoles = array(dipoles)*1e-10 #in m
    gammas = array(gammas)/1000 #from meV to eV
    energies = abs(array(energies)) #in eV
    
    neff = 3
    Lp = sum(data.layerWidths[1:]) * 1e-10 #in m
    Nq = sum(data.layerDopings[1:]*data.layerWidths[1:]) / sum(data.layerWidths[1:])
    Nq *= 100**3 #convert from cm^-3 to m^-3
    Ns = sum(data.layerDopings[1:]*1e17 * data.layerWidths[1:]*1e-8) #in cm^-2
    Ns *= 100**2 #from cm^-2 to m^-2
    hw = data.EigenE[stateR] - data.EigenE[lower]
    
#    hw = arange(0.15,0.5,0.01)
#    for enerG in hw:
#        alphaISB = sum(energies * dipoles**2 * gammas / ((energies - enerG)**2 + gammas**2))

    h = 4.1356675e-15 #eV s
    #eps0 = 5.527e7 #1/eV-m
    eps0 = 8.854e-12 #C/V-m
    #print energies/h/c0 * dipoles**2
    #print gammas / ((energies - hw)**2 + gammas**2)
    
    alphaISB = sum(energies/h/c0 * dipoles**2 * gammas / ((energies - hw)**2 + gammas**2))
    alphaISB *= 4*pi*e0**2 / (eps0*neff) * pi/(2*Lp) * Ns
    alphaISB /= e0*100
    
    if False: #plot loss diagram
        hw = arange(0.15,0.5,0.001)
        alphaISBw = zeros(hw.size)
        for q, enerG in enumerate(hw):
            alphaISBw[q] = sum(energies/h/c0 * dipoles**2 * gammas / ((energies - enerG)**2 + gammas**2))
        alphaISBw *= 4*pi*e0**2 / (eps0*neff) * pi/(2*Lp) * Ns
        alphaISBw /= e0*100
        plt.figure()
        plt.plot(hw,alphaISBw)
        plt.show()
        
    return alphaISB

def reflectivity(beta):
    beta = beta.real
    return ((beta - 1) / (beta + 1))**2



if __name__  == "__main__":
    print 'hi', cFunctions.returnme()

    
        