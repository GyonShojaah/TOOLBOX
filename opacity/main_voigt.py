import numpy as np
import sys
import datetime
import ConfigParser

from scipy import constants
from scipy import special

PRES     = 1.e3
TEMP     = 300.0
airWidth = 0.05
Tdep     = 0.5
MASS     = 18

PRES0 = 1.e3
TEMP0 = 300.0
WN   = 10000.0
dWN  = 2.0
WN_NUM = 1000

keywords = ['pressure:', 'temperature:', 'molecule:']

CC = constants.c*1e2
HH = constants.h*1e7
KK = constants.k*1e7
GG = constants.g*1e2
AMU = 1.661e-24




##############################################################################
def voigt(TT, T_0, PP, P_0, wn, airWidth, Tdep, mass) :
# def voigt(TT, T_0, PP, P_0, wn, airWidth, selfWidth, Tdep, mass, mixratio) :
    """
    =======================================================================

    voigt (TT, T_0, PP, P_0, wn, airWidth, Tdep, mass)

    -------- USAGE --------------------------------------------------------
    compute broadening function (voigt profile)
    given temperature and pressure

    -------- INPUT --------------------------------------------------------
    TT        (float)  temperature                             [K]
    T_0       (float)  standard temperature                    [K]
    PP        (float)  total pressure                          [barye]
    P_0       (float)  standard pressure                       [barye]
    wn        (float)  wave number                             [cm-1]
    airWidth  (float)  air-broadening coefficients             [cm-1]
    Tdep      (float)  power of TT depedence of Lorentz width  [-]
    mass      (float)  molecular mass                          [AMU]

    -------- OUTPUT --------------------------------------------------------
    voigtfunc (func)   Voigt Function

    =======================================================================
    """

    # QUESTION!!!
    # p_total = PP
    # p_partial = PP*mixratio
    # gamma_L = (airWidth*((p_total-p_partial)/P_0)+selfWidth*(p_partial/P_0))*(T_0/TT)**Tdep
    #
#    print mass, wn
    gamma_L = airWidth*(PP/P_0)*(T_0/TT)**Tdep
    gamma_D = wn*np.sqrt(2*np.log(2)*KK*TT/(mass*AMU*CC**2))

    print "################################"
    print "gamma_L", gamma_L
    print "gamma_D", gamma_D
    print "y=", np.sqrt(np.log(2.))*gamma_L/gamma_D
    print "################################"
    #
    # Voigt function
    #
    # integral_-infty^+infty dx' g_L(x-x')*g_D(x'-x0)
    #   pressure broadening : g_L(x) = gamma_L/pi/((x-x0)^2+gamma_L^2)
    #   Doppler broadening  : g_D(x) = (1/gamma_D)*(ln(2)/pi)^(1/2)*exp(-ln(2)*((x-x0)/gamma_D))
    #
    # Voigt function is the real part of the Faddeeva function
    # w(x+iy) = V(x,y) + iL(x,y)
    #
    def return_voigtfunc(x) : 
        realvalue = np.sqrt(np.log(2))/gamma_D * (x-wn)
        imagvalue = np.sqrt(np.log(2))*gamma_L/gamma_D
        complexvalue = realvalue + imagvalue*1j
        faddeeva = np.sqrt(np.log(2))/gamma_D/np.sqrt(np.pi) * special.wofz(complexvalue)
        return faddeeva.real
    return return_voigtfunc


##############################################################################
def doppler(TT, T_0, PP, P_0, wn, airWidth, Tdep, mass) :

    """
    =======================================================================

    voigt (TT, T_0, PP, P_0, wn, airWidth, Tdep, mass)

    -------- USAGE --------------------------------------------------------
    compute broadening function (voigt profile)
    given temperature and pressure

    -------- INPUT --------------------------------------------------------
    TT        (float)  temperature                             [K]
    T_0       (float)  standard temperature                    [K]
    PP        (float)  total pressure                          [barye]
    P_0       (float)  standard pressure                       [barye]
    wn        (float)  wave number                             [cm-1]
    airWidth  (float)  air-broadening coefficients             [cm-1]
    Tdep      (float)  power of TT depedence of Lorentz width  [-]
    mass      (float)  molecular mass                          [AMU]

    -------- OUTPUT --------------------------------------------------------
    voigtfunc (func)   Voigt Function

    =======================================================================
    """

    gamma_D = wn*np.sqrt(2*np.log(2)*KK*TT/(mass*AMU*CC**2))

    def return_dopplerfunc(x) : 
        fac = 1./gamma_D*(np.log(2.)/np.pi)**0.5
        return fac * np.exp( - np.log(2.) * ( ( x - wn )/gamma_D )**2 )

    return return_dopplerfunc



##############################################################################
def lorentz(TT, T_0, PP, P_0, wn, airWidth, Tdep, mass) :
# def voigt(TT, T_0, PP, P_0, wn, airWidth, selfWidth, Tdep, mass, mixratio) :
    """
    =======================================================================

    voigt (TT, T_0, PP, P_0, wn, airWidth, Tdep, mass)

    -------- USAGE --------------------------------------------------------
    compute broadening function (voigt profile)
    given temperature and pressure

    -------- INPUT --------------------------------------------------------
    TT        (float)  temperature                             [K]
    T_0       (float)  standard temperature                    [K]
    PP        (float)  total pressure                          [barye]
    P_0       (float)  standard pressure                       [barye]
    wn        (float)  wave number                             [cm-1]
    airWidth  (float)  air-broadening coefficients             [cm-1]
    Tdep      (float)  power of TT depedence of Lorentz width  [-]
    mass      (float)  molecular mass                          [AMU]

    -------- OUTPUT --------------------------------------------------------
    voigtfunc (func)   Voigt Function

    =======================================================================
    """

    gamma_L = airWidth*(PP/P_0)*(T_0/TT)**Tdep
    def return_lorentzfunc(x) : 
        fac = gamma_L / np.pi
        return fac / ( ( x - wn )**2 + gamma_L**2 )
    return return_lorentzfunc




##############################################################################
if __name__ == "__main__":

    voigtfunc   = voigt( TEMP, TEMP0, PRES, PRES0, WN, airWidth, Tdep, MASS ) # unit: [cm]
    dopplerfunc = doppler( TEMP, TEMP0, PRES, PRES0, WN, airWidth, Tdep, MASS ) # unit: [cm]
    lorentzfunc = lorentz( TEMP, TEMP0, PRES, PRES0, WN, airWidth, Tdep, MASS ) # unit: [cm]

    wn_array = np.linspace(WN-dWN, WN+dWN, WN_NUM)

    for ii in xrange(WN_NUM):
        print wn_array[ii], voigtfunc(wn_array[ii]), dopplerfunc(wn_array[ii]), lorentzfunc(wn_array[ii])

