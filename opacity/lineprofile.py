import numpy as np
from cgs import *
from molecules import molecules
from scipy import special

##############################################################################
def voigtfunction(TT, T_0, PP, P_0, wn, airWidth, Tdep, mass) :
# def voigtfunction(TT, T_0, PP, P_0, wn, airWidth, selfWidth, Tdep, mass, mixratio) :
    """
    =======================================================================

    voigtfunction (TT, T_0, PP, P_0, wn)

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
    #
    # QUESTION!!!
    # p_total = PP
    # p_partial = PP*mixratio
    # gamma_L = (airWidth*((p_total-p_partial)/P_0)+selfWidth*(p_partial/P_0))*(T_0/TT)**Tdep
    #
#    print mass, wn
    gamma_L = airWidth*(PP/P_0)*(T_0/TT)**Tdep
    gamma_D = wn*np.sqrt(2*np.log(2)*KK*TT/(mass*AMU*CC**2))
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
def linestrength(molecule, TT, T_0, wn, strength, energy) :
    """
    =======================================================================

    linestrength(TT, T_0, wn, energy) :

    -------- USAGE --------------------------------------------------------
    compute line strength
    given temperature and pressure

    -------- INPUT --------------------------------------------------------
    TT        (float)  temperature                             [K]
    T_0       (float)  standard temperature                    [K]
    wn        (float)  wavenumber of the transition            [cm-1/cm-2]
    energy    (float)  energy of the transition                [cm-1/cm-2]
    strength  (float)  line strength at standard condition     [cm-1/cm-2]

    -------- OUTPUT --------------------------------------------------------
    new_strength  (float)  new line strength                   [-]

    -------- REFERENCE -----------------------------------------------------
    Norton and Rinsland (1991)

    =======================================================================
    """

    beta = molecules[molecule]['TempExpQR']
    freq = molecules[molecule]['VibFreq']
    dege = molecules[molecule]['NumDeg']

    fac_exponent = np.exp(-HH*CC*energy/(KK*TT))/np.exp(-HH*CC*energy/(KK*T_0))
#    print 'fac_exponent', fac_exponent
    fac_exponent *= fac_exponent * (1 - np.exp(-HH*CC*wn/(KK*TT)))/(1 - np.exp(-HH*CC*wn/(KK*T_0)))
    fac_Q_rot = (T_0/TT)**beta
    fac_Q_vib = 1.0
    for i in xrange(len(freq)) :
        fac_Q_vib *= (( 1 - np.exp(-(HH*CC*freq[i])/(KK*TT)) )/( 1 - np.exp(-(HH*CC*freq[i])/(KK*T_0)) ))**dege[i]

    new_strength = (fac_exponent*fac_Q_rot*fac_Q_vib)*strength
    return new_strength
