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


class MaterialConstants(object):
    def __init__(self):
        Temperature = 300
        self.set_constants(Temperature)
        
    def set_temperature(self, Temperature):
        self.set_constants(Temperature)
        
    def set_constants(self, Temperature):
        self.Temperature = Temperature
        
        # GaAs constants
         #from Vurgaftman
        self.alc_GaAs = 5.65325 + 3.88e-5*(Temperature-300) #Angs at 80 K
        self.c11_GaAs = 1221 #GPa
        self.c12_GaAs = 566  #GPa
        self.EgG_GaAs = 1.519 #eV at 0 K 
        self.EgL_GaAs = 1.815 #eV
        self.EgX_GaAs = 1.981 #eV
        self.VBO_GaAs = -0.80 #eV
        self.DSO_GaAs = 0.341 #eV Delta_SO from top of valence band
        self.acG_GaAs = -7.17 #eV
        self.acL_GaAs = -4.91 #eV, NextNano DB
        self.acX_GaAs = -0.16 #eV, NaxtNano DB
        self.av_GaAs  = -1.16 #eV, valence band deformation potential
        self.b_GaAs   = -2.0  #eV, Pikus Burr Uniaxial Strain deformation potential
        self.XiG_GaAs = 0 #eV, uniaxial deformation potential, NextNano DB
        self.XiL_GaAs = 6.5 #eV
        self.XiX_GaAs = 14.26 #eV
        self.me0_GaAs = 0.067 #1/m0
        self.Ep_GaAs  = 28.8 #eV
        self.F_GaAs   = -1.94
        self.alG_GaAs = 0.5405e-3 #eV/K, Varshni alpha(Gamma)
        self.beG_GaAs = 204 #K, Varshni beta(Gamma)
        self.alX_GaAs = 0.460e-3
        self.beX_GaAs = 204
        self.alL_GaAs = 0.605e-3
        self.beL_GaAs = 204
        self.epss_GaAs = 12.9
        self.epsInf_GaAs = 10.86
        self.hwLO_GaAs = 35.3*1e-3
        self.C1_GaAs = 3.5 #Handbook of Optics, 2nd edition, Vol. 2. McGraw-Hill 1994
        self.C2_GaAs = 7.4969
        self.C3_GaAs = 0.4082
        self.C4_GaAs = 1.9347
        self.C5_GaAs = 37.17
        
        # InAs constants
         #from Vurgaftman
        self.alc_InAs = 6.0583 + 2.74e-5*(Temperature-300)
        self.c11_InAs = 832.9
        self.c12_InAs = 452.6
        self.EgG_InAs = 0.417
        self.EgL_InAs = 1.133
        self.EgX_InAs = 1.433
        self.VBO_InAs = -0.59
        self.DSO_InAs = 0.39
        self.acG_InAs = -5.08
        self.acL_InAs = -3.89 #eV, NextNano DB
        self.acX_InAs = -0.08 #eV, NextNano DB
        self.av_InAs  = -1.00
        self.b_InAs   = -1.8
        self.XiG_InAs = 0
        self.XiL_InAs = 11.35
        self.XiX_InAs = 3.7
        self.me0_InAs = 0.026
        self.Ep_InAs  = 21.5
        self.F_InAs   = -2.9
        self.alG_InAs = 0.276e-3 
        self.beG_InAs = 93 
        self.alX_InAs = 0.276e-3
        self.beX_InAs = 93
        self.alL_InAs = 0.276e-3
        self.beL_InAs = 93
        self.epss_InAs = 14.3
        self.epsInf_InAs = 11.6
        self.hwLO_InAs = 29.93*1e-3
        self.C1_InAs = 11.1 #Handbook of Optics, 2nd edition, Vol. 2. McGraw-Hill 1994
        self.C2_InAs = 0.71
        self.C3_InAs = 2.551
        self.C4_InAs = 2.75
        self.C5_InAs = 45.66
          
         
        # AlAs constants
         #from Vurgaftman
        self.alc_AlAs = 5.6611 + 2.90e-5*(Temperature-300)
        self.c11_AlAs = 1250
        self.c12_AlAs = 534
        self.EgG_AlAs = 3.099
        self.EgL_AlAs = 2.46
        self.EgX_AlAs = 2.24
        self.VBO_AlAs = -1.33
        self.DSO_AlAs = 0.28
        self.acG_AlAs = -5.64
        self.acL_AlAs = -3.07 #NextNano DB
        self.acX_AlAs = 2.54  #NextNano DB
        self.av_AlAs  = -2.47
        self.b_AlAs   = -2.3
        self.XiG_AlAs = 0
        self.XiL_AlAs = 11.35
        self.XiX_AlAs = 6.11
        self.me0_AlAs = 0.15
        self.Ep_AlAs  = 21.1
        self.F_AlAs   = -0.48
        self.alG_AlAs = 0.855e-3
        self.beG_AlAs = 530
        self.alX_AlAs = 0.70e-3
        self.beX_AlAs = 530
        self.alL_AlAs = 0.605e-3
        self.beL_AlAs = 204
        self.epss_AlAs = 10.06
        self.epsInf_AlAs = 8.16
        self.hwLO_AlAs = 49.8*1e-3
        self.C1_AlAs = 2.0792 #Handbook of Optics, 2nd edition, Vol. 2. McGraw-Hill 1994
        self.C2_AlAs = 6.0840
        self.C3_AlAs = 0.2822
        self.C4_AlAs = 1.900
        self.C5_AlAs = 27.62
        
        
        # AlSb constants
        # from Vurgaftman
        self.alc_AlSb = 6.1355 + 2.60e-5*(Temperature-300)
        self.c11_AlSb = 876.9
        self.c12_AlSb = 434.1
        self.EgG_AlSb = 2.386
        self.EgL_AlSb = 2.329
        self.EgX_AlSb = 1.696
        self.VBO_AlSb = -0.41
        self.DSO_AlSb = 0.676
        self.acG_AlSb = -4.5
        self.acL_AlSb = 0 #NextNano DB
        self.acX_AlSb = 2.54  #NextNano DB
        self.av_AlSb  = -1.4
        self.b_AlSb   = -1.35
        self.XiG_AlSb = 0
        self.XiL_AlSb = 11.35
        self.XiX_AlSb = 6.11
        self.me0_AlSb = 0.14
        self.Ep_AlSb  = 18.7
        self.F_AlSb   = -0.56
        self.alG_AlSb = 0.42e-3
        self.beG_AlSb = 140
        self.alX_AlSb = 0.39e-3
        self.beX_AlSb = 140
        self.alL_AlSb = 0.58e-3
        self.beL_AlSb = 140
        self.epss_AlSb = 12.04 #ISBN 0849389127
        self.epsInf_AlSb = 10.24 #ISBN 0849389127
        self.hwLO_AlSb = 42.7 #http://prb.aps.org/pdf/PRB/v43/i9/p7231_1
        
        # GaSb constants
        # from Vurgaftman
        self.alc_GaSb = 6.0959 + 4.72e-5*(Temperature-300)
        self.c11_GaSb = 884.2
        self.c12_GaSb = 402.6
        self.EgG_GaSb = 0.812
        self.EgL_GaSb = 0.875
        self.EgX_GaSb = 1.141
        self.VBO_GaSb = -0.03
        self.DSO_GaSb = 0.76
        self.acG_GaSb = -7.5
        self.acL_GaSb = 0 #unknown
        self.acX_GaSb = 0 #unknown
        self.av_GaSb  = -0.8
        self.b_GaSb   = -2.0
        self.XiG_GaSb = 0 #unknown
        self.XiL_GaSb = 0 #unknown
        self.XiX_GaSb = 0 #unknown
        self.me0_GaSb = 0.039
        self.Ep_GaSb  = 27.0
        self.F_GaSb   = -1.63
        self.alG_GaSb = 0.417e-3
        self.beG_GaSb = 140
        self.alX_GaSb = 0.475e-3
        self.beX_GaSb = 94
        self.alL_GaSb = 0.597e-3
        self.beL_GaSb = 140
        self.epss_GaSb = 0 #unknown
        self.epsInf_GaSb = 0 #unknown
        self.hwLO_GaSb = 0 #unknown
        
        # InSb constants
        # from Vurgaftman
        self.alc_InSb = 6.4794 + 3.48e-5*(Temperature-300)
        self.c11_InSb = 684.7
        self.c12_InSb = 373.5
        self.EgG_InSb = 0.235
        self.EgL_InSb = 0.93
        self.EgX_InSb = 0.63
        self.VBO_InSb = 0
        self.DSO_InSb = 0.81
        self.acG_InSb = -6.94
        self.acL_InSb = 0 #unknown
        self.acX_InSb = 0 #unknown
        self.av_InSb  = -0.36
        self.b_InSb   = -2.0
        self.XiG_InSb = 0 #unknown
        self.XiL_InSb = 0 #unknown
        self.XiX_InSb = 0 #unknown
        self.me0_InSb = 0.0135
        self.Ep_InSb  = 23.3
        self.F_InSb   = -0.23
        self.alG_InSb = 0.32e-3
        self.beG_InSb = 170
        self.alX_InSb = 0e-3 #unknown
        self.beX_InSb = 0 #unknown
        self.alL_InSb = 0e-3 #unknown
        self.beL_InSb = 0 #unknown
        self.epss_InSb = 0 #unknown
        self.epsInf_InSb = 0 #unknown
        self.hwLO_InSb = 0 #unknown
        
        # InP constants
        self.alc_InP = 5.8697+2.79e-5*(Temperature-300)
        self.me0_InP = 0.0795
        self.C1_InP = 7.255 #Handbook of Optics, 2nd edition, Vol. 2. McGraw-Hill 1994
        self.C2_InP = 2.316
        self.C3_InP = 0.6263
        self.C4_InP = 2.765
        self.C5_InP = 32.935
        
        # InGaAs constants
        self.EgG_InGaAs = 0.477
        self.EgL_InGaAs = 0.33
        self.EgX_InGaAs = 1.4
        self.VBO_InGaAs = -0.38
        self.DSO_InGaAs = 0.15
        self.acG_InGaAs = 2.61
        self.acL_InGaAs = 2.61 #NextNano DB
        self.acX_InGaAs = 2.61 #NextNano DB
        self.me0_InGaAs = 0.0091
        self.Ep_InGaAs  = -1.48
        self.F_InGaAs   = 1.77
         
        # AlInAs constants
        self.EgG_AlInAs = 0.70
        self.EgL_AlInAs = 0
        self.EgX_AlInAs = 0
        self.VBO_AlInAs = -0.64
        self.DSO_AlInAs = 0.15
        self.acG_AlInAs = -1.4
        self.acL_AlInAs = -1.4 #NextNano DB
        self.acX_AlInAs = -1.4 #NextNano DB
        self.me0_AlInAs = 0.049
        self.Ep_AlInAs  = -4.81
        self.F_AlInAs   = -4.44
         
        # AlGaAs constants
        self.EgG_AlGaAs = -0.127 #+ 1.310*x
        self.EgL_AlGaAs = 0.055
        self.EgX_AlGaAs = 0
        self.VBO_AlGaAs = 0
        self.DSO_AlGaAs = 0
        self.acG_AlGaAs = 0
        self.acL_AlGaAs = 0
        self.acX_AlGaAs = 0
        self.me0_AlGaAs = 0
        self.Ep_AlGaAs  = 0
        self.F_AlGaAs   = 0
        
        # AlAsSb constants
        self.EgG_AlAsSb = 0.8
        self.EgL_AlAsSb = 0.28
        self.EgX_AlAsSb = 0.28
        self.DSO_AlAsSb = 0.15
        self.VBO_AlAsSb = -1.71
        
        # AlGaSb constants
        self.EgG_AlGaSb = -0.044 #+ 1.22*x
        self.EgL_AlGaSb = 0
        self.EgX_AlGaSb = 0
        self.VBO_AlGaSb = 0
        self.DSO_AlGaSb = 0.2
        self.acG_AlGaSb = 0
        self.acL_AlGaSb = 0 #NextNano DB
        self.acX_AlGaSb = 0 #NextNano DB
        self.me0_AlGaSb = 0
        self.Ep_AlGaSb  = 0
        self.F_AlGaSb   = 0
        
        # InAsSb constants
        self.EgG_InAsSb = 0.67
        self.EgL_InAsSb = 0.6
        self.EgX_InAsSb = 0.6
        self.VBO_InAsSb = 0
        self.DSO_InAsSb = 1.2
        self.acG_InAsSb = 0
        self.acL_InAsSb = 0 #NextNano DB
        self.acX_InAsSb = 0 #NextNano DB
        self.me0_InAsSb = 0.035
        self.Ep_InAsSb  = 0
        self.F_InAsSb   = 0
